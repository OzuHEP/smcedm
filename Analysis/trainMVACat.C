#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

void trainMVACat();

void trainMVACat()
{
    char name[1000];
    float xsec;
    float weight;
    int n_events;

    TFile *outf = new TFile("mva_SMCEDM.root", "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("factory_mva_SMCEDM", outf, "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");

    TMVA::DataLoader loader("dataset");

    TFile *sigSrc[1];
    TFile *bkgSrc[2];

    sigSrc[0] = TFile::Open("/mnt/harddisk4/scratch/signal_dtG1_delphes.root");

    bkgSrc[0] = TFile::Open("/mnt/harddisk4/scratch/dyjets_delphes.root");
    bkgSrc[1] = TFile::Open("/mnt/harddisk4/scratch/w_jets_delphes.root ");

    TTree *sigTree[1];
    TTree *bkgTree[2];

    for (int k = 0; k < 1; k++)
    {
        xsec = 0.;
        sigTree[k] = (TTree *)sigSrc[0]->FindObjectAny("outtree");
        sigTree[k]->SetBranchAddress("br_weight",&weight);
        n_events = sigTree[k]->GetEntries();

        for (int i_event = 0; i_event < n_events; ++i_event)
        {
            sigTree[k]->GetEntry(i_event);
            xsec += weight;
        }

        loader.AddSignalTree(sigTree[k], 1.);
    }


    for (int k = 0; k < 2; k++)
    {
        xsec = 0.;
        bkgTree[k] = (TTree *)bkgSrc[k]->FindObjectAny("outtree");
        bkgTree[k]->SetBranchAddress("br_weight",&weight);
        n_events = bkgTree[k]->GetEntries();

        for (int i_event = 0; i_event < n_events; ++i_event)
        {
            bkgTree[k]->GetEntry(i_event);
            xsec += weight;
        }

        loader.AddBackgroundTree(bkgTree[k], 0.1);
    }

    //const int NVAR = 21;
    //TString VAR[NVAR] = {"nJets", "ht", "jetPt[0]", "jetPt[1]", "jetPt[2]", "jetPt[3]", "jetPt[4]", "sphericity", "aplanarity", "foxWolfram[0]", "foxWolfram[1]", "foxWolfram[2]", "foxWolfram[3]", "mTop", "yTop", "ptTop", "met", "metPhi", "nLeptons", "mW", "nBJets"};
    //char TYPE[NVAR] = {"I", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "I", "F", "I"};

    const int NVAR = 15;
    TString VAR[NVAR] = {"br_njets", "br_nbjets", "br_scalar_HT", "br_jet_pt[0]", "br_jet_pt[1]", "br_jet_pt[2]", "br_jet_pt[3]", "br_MET", "br_MET_Phi", "br_sphericity", "br_aplanarity","br_Fox_Wolfram[0]", "br_Fox_Wolfram[1]", "br_Fox_Wolfram[2]", "br_Fox_Wolfram[3]"};

    //char TYPE[NVAR] = {"I", "I", "I", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F"};

    for (int i = 0; i < NVAR; i++)
    {
        loader.AddVariable(VAR[i]);
    }

    TCut mycuts;
    int nTrain_Signal     = 3000;
    int nTrain_Background = 3000;

    loader.PrepareTrainingAndTestTree(mycuts, Form("nTrain_Signal=%d:nTrain_Background=%d:SplitMode=Random:NormMode=NumEvents:!V", nTrain_Signal, nTrain_Background));

    //factory->BookMethod( &loader, TMVA::Types::kKNN, "KNN",
    //                     "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:CreateMVAPdfs:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );

    factory->BookMethod(&loader, TMVA::Types::kBDT, "BDT", "!V:NTrees=200:MinNodeSize=20.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

    //factory->BookMethod( &loader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

    // Use Deep Neural-Network

    //////General layout.

    TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

    //////Training strategies.
    TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=100,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
    TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=100,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=100,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString trainingStrategyString ("TrainingStrategy=");
    trainingStrategyString += training0 + "|" + training1 + "|" + training2;

    //////General Options.
    TString dnnOptions ("!H:V:ErrorStrategy=CROSSENTROPY:VarTransform=N:"
                        "WeightInitialization=XAVIERUNIFORM");
    dnnOptions.Append (":");
    dnnOptions.Append (layoutString);
    dnnOptions.Append (":");
    dnnOptions.Append (trainingStrategyString);

    TString cpuOptions = dnnOptions + ":Architecture=GPU:";

    factory->BookMethod( &loader, TMVA::Types::kDNN, "DNN GPU", cpuOptions);

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    outf->Close();
}
