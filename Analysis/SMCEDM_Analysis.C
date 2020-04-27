#include "SMCEDM_Analysis.h"
#include "progressbar.hpp"

void SMCEDM_Analysis(int event_begin, int event_end) {
	gROOT->SetBatch();
	int nthreads = 24;
	ROOT::EnableImplicitMT(nthreads);

	auto file = new TFile("delphes.root");
	//auto file = new TFile("/mnt/harddisk2/wjets/WJets_200to400/Events/delphes.root");
	//auto file = new TFile("/mnt/harddisk1/ttbar/dtG_2/Events/delphes.root");
	TTreeReader myReader("Delphes", file);

	char outfile_name[50];
	sprintf(outfile_name, "output_cms_wjets.root", event_begin, event_end);
	auto outfile = new TFile(outfile_name, "RECREATE");

	// Create the TTree and branches
	TTree *outtree = new TTree("outtree", "outtree");

	int   n_jets;
	int n_leptons;
	int n_bjets;

	float scalar_ht;
	float met;
	float met_phi;
	float aplanarity;
	float sphericity;

	std::vector<float> jet_pt;
	std::vector<float> fox_wolfram;

	outtree->Branch("br_njets"       , &n_jets,      "n_jets/I");
	outtree->Branch("br_nleptons"    , &n_leptons,   "n_leptons/I");
	outtree->Branch("br_nbjets"      , &n_bjets,     "n_bjets/I");
	outtree->Branch("br_scalar_HT"   , &scalar_ht,   "scalar_ht/F");
	outtree->Branch("br_jet_pt"      , &jet_pt);
	outtree->Branch("br_MET"         , &met,         "MET/F");
	outtree->Branch("br_MET_Phi"     , &met_phi,     "MET_Phi/F");
	outtree->Branch("br_sphericity"  , &sphericity,  "sphericity/F");
	outtree->Branch("br_Fox_Wolfram" , &fox_wolfram);

	//Jet Definitions
	TTreeReaderValue<int>            rv_Jet_size(myReader, "Jet_size");
	TTreeReaderArray<float>          ra_Jet_pT(myReader,   "Jet.PT");
	TTreeReaderArray<float>          ra_Jet_Eta(myReader,  "Jet.Eta");
	TTreeReaderArray<float>          ra_Jet_Phi(myReader,  "Jet.Phi");
	TTreeReaderArray<float>          ra_Jet_Mass(myReader, "Jet.Mass");
	TTreeReaderArray<unsigned int>   ra_Jet_BTag(myReader, "Jet.BTag");

	//Electron Definitions
	TTreeReaderValue<int>            rv_Electron_size(myReader,  "Electron_size");
	TTreeReaderArray<float>          ra_Electron_pT(myReader,    "Electron.PT");
	TTreeReaderArray<float>          ra_Electron_Eta(myReader,    "Electron.Eta");
	TTreeReaderArray<float>          ra_Electron_Phi(myReader,    "Electron.Phi");
	TTreeReaderArray<float>          ra_Electron_Energy(myReader,   "Electron.T");
	TTreeReaderArray<int>            ra_Electron_Charge(myReader, "Electron.Charge");

	//Muon Definitions
	TTreeReaderValue<int>            rv_Muon_size(myReader,  "Muon_size");
	TTreeReaderArray<float>          ra_Muon_pT(myReader,    "Muon.PT");
	TTreeReaderArray<float>          ra_Muon_Eta(myReader,    "Muon.Eta");
	TTreeReaderArray<float>          ra_Muon_Phi(myReader,    "Muon.Phi");
	TTreeReaderArray<float>          ra_Muon_Energy(myReader,   "Muon.T");
	TTreeReaderArray<int>            ra_Muon_Charge(myReader, "Muon.Charge");

	//MET Definitions
	TTreeReaderArray<float>          ra_MissingET_MET(myReader,   "MissingET.MET");
	TTreeReaderArray<float>          ra_MissingET_Phi(myReader,   "MissingET.Phi");

	//HT Definitions
	TTreeReaderArray<float>          ra_scalar_HT(myReader,   "ScalarHT.HT");

	// Histogram Definitions
	auto histo_bjet_pT      = new TH1F("histo_bjet_pT", "histo_bjet_pT", 200, 0, 2000);
	auto histo_light_jet_pT = new TH1F("histo_light_jet_pT", "histo_light_jet_pT", 200, 0, 2000);
	auto histo_had_top_mass = new TH1F("histo_had_top_mass", "histo_had_top_mass", 1000, 0, 1000);
	auto histo_had_top_eta  = new TH1F("histo_had_top_eta", "histo_had_top_eta", 100, -5, 5);
	auto histo_had_top_phi  = new TH1F("histo_had_top_phi", "histo_had_top_phi", 30, -3.14, 3.14);

	auto n_events = myReader.GetEntries(1);

	int n_btag;
	int n_btag_cut = 0;
	int n_light_jet;
	int n_light_jet_cut = 0;
	int n_lepton_cut = 0;
	int n_MET_cut  = 0;
	int n_events_passed = 0;
	float MET_threshold = 0.0;

	// initialize the progress bar
	const int limit = myReader.GetEntries(1);
	ProgressBar progressBar(limit, 70);

	if (event_end > n_events) { event_end = n_events; }

	// Loop over Events
	for (int i_event = event_begin; i_event < event_end; ++i_event) {
		++progressBar;
		// display the bar
		progressBar.display();

		n_btag = 0;

		bool pass_lepton_criteria(1), pass_b_jet_criteria(1), pass_MET_criteria(1), pass_light_jet_criteria(1);

		myReader.SetEntry(i_event);

		// put all electrons to a container (vector of leptons object)
		std::vector<leptons> v_electrons;
		leptons dummy_lepton;

		//******* count the number of negative and positive electrons *******//
		for (int i_electron = 0; i_electron < *rv_Electron_size; ++i_electron) {
			dummy_lepton.setPt(ra_Electron_pT.At(i_electron));
			dummy_lepton.setEta(ra_Electron_Eta.At(i_electron));
			dummy_lepton.setPhi(ra_Electron_Phi.At(i_electron));
			dummy_lepton.setM(0.511 / 1000);
			dummy_lepton.setQ(ra_Electron_Charge.At(i_electron));
			dummy_lepton.setPDGID(11);

			v_electrons.push_back(dummy_lepton);
		}

		// put all muons to a container (vector of leptons object)
		std::vector<leptons> v_muons;
		leptons dummy_muon;
		//******* count the number of negative and positive muons *******//
		for (int i_muon = 0; i_muon < *rv_Muon_size; ++i_muon) {
			dummy_lepton.setPt(ra_Muon_pT.At(i_muon));
			dummy_lepton.setEta(ra_Muon_Eta.At(i_muon));
			dummy_lepton.setPhi(ra_Muon_Phi.At(i_muon));
			dummy_lepton.setM(105 / 1000);
			dummy_lepton.setQ(ra_Muon_Charge.At(i_muon));
			dummy_lepton.setPDGID(13);

			v_muons.push_back(dummy_lepton);
		}

		// put all jets to a container (vector of jets object)
		std::vector<jets> v_jets;
		jets dummy_jet;
		//******* count the number of negative and positive muons *******//
		for (int i_jet = 0; i_jet < *rv_Jet_size; ++i_jet) {
			dummy_jet.setPt(ra_Jet_pT.At(i_jet));
			dummy_jet.setEta(ra_Jet_Eta.At(i_jet));
			dummy_jet.setPhi(ra_Jet_Phi.At(i_jet));
			dummy_jet.setM(ra_Jet_Mass.At(i_jet));
			dummy_jet.setBTag(ra_Jet_BTag.At(i_jet));

			v_jets.push_back(dummy_jet);
		}

		if (v_electrons.size() + v_muons.size() < 1) {
			pass_lepton_criteria = 0;
			n_lepton_cut++;
		}

		auto v_light_jets = create_light_jet_collection(v_jets);
		auto v_b_jets     = create_b_jet_collection(v_jets);

		n_btag = v_b_jets.size() - v_light_jets.size();

		if (v_light_jets.size() < 2) {
			pass_light_jet_criteria = 0;
			n_light_jet_cut++;
		}

		if (n_btag < 2) {
			pass_b_jet_criteria = 0;
			n_btag_cut++; }

		if (ra_MissingET_MET[0] < MET_threshold) {
			pass_MET_criteria = 0;
			n_MET_cut++;
		}

		// EVENT SELECTION (skips event if not matching all the criteria below)
		if (!pass_lepton_criteria || !pass_MET_criteria || !pass_light_jet_criteria) {
			continue;
		}

		auto v_W_candidates = create_W_candidates(v_light_jets);

		if (n_btag > 0) {
			auto v_had_top_candidates = create_had_top_candidates(v_b_jets, v_W_candidates);
		}

		for (auto& jet : v_jets) {
			jet_pt.push_back(jet.getPt());
		}

		n_jets     = v_jets.size();
		n_bjets    = v_b_jets.size();
		n_leptons  = v_electrons.size() + v_muons.size();
		scalar_ht  = ra_scalar_HT[0];

		aplanarity  = get_event_shape_variable(v_jets, "aplanarity");
		sphericity  = get_event_shape_variable(v_jets, "sphericity");
		fox_wolfram = get_fox_wolfram_parameters(v_jets, 4);

		met     = ra_MissingET_MET[0];
		met_phi = ra_MissingET_Phi[0];

		outtree->Fill();
		n_events_passed++;
	}

	outfile->cd();
	outtree->Write();

	std::cout << "n_light_jet_cut: " << n_light_jet_cut << std::endl;
	std::cout << "n_btag_cut: " << n_btag_cut << std::endl;
	std::cout << "n_lepton_cut: " << n_lepton_cut << std::endl;
	std::cout << "n_MET_cut: " << n_MET_cut << std::endl;
	std::cout << "n_events_passed/n_events: " << n_events_passed << "/" << n_events << std::endl;
}
