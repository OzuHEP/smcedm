void get_number_of_events(TString file_name, TString tree_name)
{
	auto file     = new TFile(file_name);
	TTreeReader myReader(tree_name, file);
	auto n_events = myReader.GetEntries(1);
	ofstream cout("output.txt");
    cout << n_events << endl;
    gROOT->ProcessLine(".q");
}