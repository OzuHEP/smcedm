#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

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
    bkgSrc[1] = TFile::Open("/mnt/harddisk4/scratch/wjets_delphes.root");

    TTree *sigTree[1];
    TTree *bkgTree[2];

    int n_signal = 0;
    int n_bkg    = 0;

    for (int k = 0; k < 1; k++)
    {
        xsec = 0.;
        sigTree[k] = (TTree *)sigSrc[0]->FindObjectAny("outtree");
        sigTree[k]->SetBranchAddress("br_weight", &weight);
        n_events = sigTree[k]->GetEntries();
        n_signal += n_events;

        for (int i_event = 0; i_event < n_events; ++i_event)
        {
            sigTree[k]->GetEntry(i_event);
            xsec += weight;
        }

        std::cout << "xsec_signal[" << k << "]=" << xsec << std::endl;
        loader.AddSignalTree(sigTree[k], xsec/n_events);
    }

    for (int k = 0; k < 2; k++)
    {
        xsec     = 0.;
        bkgTree[k] = (TTree *)bkgSrc[k]->FindObjectAny("outtree");
        //bkgTree[k]->SetBranchAddress("br_weight", &weight);
        n_events = bkgTree[k]->GetEntries();
        n_bkg += n_events;

        for (int i_event = 0; i_event < n_events; ++i_event)
        {
            bkgTree[k]->GetEntry(i_event);
            xsec += weight;
        }

        std::cout << "xsec_background[" << k << "]=" << xsec << std::endl;
        loader.AddBackgroundTree(bkgTree[k], xsec/n_events);
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

    int nTrain_Background = (int)(n_bkg*0.9);
    int nTrain_Signal     = nTrain_Background;
    int nTest_Signal      = (int)(nTrain_Signal*0.1);

    loader.PrepareTrainingAndTestTree(mycuts, Form("nTrain_Signal=%d:nTrain_Background=%d:nTest_Signal=%d::SplitMode=Random:NormMode=NumEvents:!V", nTrain_Signal, nTrain_Background, nTest_Signal));


    //factory->BookMethod( &loader, TMVA::Types::kRuleFit, "RuleFit","H:!V:RuleFitModule=RFTMVA:Model=ModRuleLinear:MinImp=0.001:RuleMinDist=0.001:NTrees=20:fEventsMin=0.01:fEventsMax=0.5:GDTau=-1.0:GDTauPrec=0.01:GDStep=0.01:GDNSteps=10000:GDErrScale=1.02" );

    factory->BookMethod(&loader, TMVA::Types::kBDT, "BDT", "!V:NTrees=200:MinNodeSize=20.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

    factory->BookMethod( &loader, TMVA::Types::kBDT, "BDTG","!H:!V:NTrees=1000:MinNodeSize=2.5%:BoostType=Grad:Shrinkage=0.10:UseBaggedBoost:BaggedSampleFraction=0.5:nCuts=20:MaxDepth=2" );

    //factory->BookMethod( &loader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

    //factory->BookMethod( &loader, TMVA::Types::kMLP, "MLPBFGS", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:!UseRegulator" );

    factory->BookMethod( &loader, TMVA::Types::kMLP, "MLPBNN", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=60:HiddenLayers=N+5:TestRate=5:TrainingMethod=BFGS:UseRegulator" ); // BFGS training with bayesian regulators

    // Use Deep Neural-Network

    //////General layout.

    TString layoutString ("Layout=TANH|256,TANH|256,TANH|256,LINEAR");

    //////Training strategies.
    TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                      "ConvergenceSteps=200,BatchSize=32,TestRepetitions=1,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.2+0.2+0.2, Multithreading=True");
    TString training1("LearningRate=1e-2,Momentuwrapm=0.9,Repetitions=1,"
                      "ConvergenceSteps=200,BatchSize=32,TestRepetitions=1,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                      "ConvergenceSteps=200,BatchSize=32,TestRepetitions=1,"
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
