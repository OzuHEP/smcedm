#include "TMVA/Factory.h"
#include "TMVA/MethodCategory.h"
#include "TFile.h"
#include "TTree.h"
#include "TCut.h"
#include "TH1F.h"
using namespace TMVA;
using namespace TMath;

void trainMVACat();

void trainMVACat() {
    char name[1000];
    // backgrounds:[ttbar, wjets]
    float sig_XSEC[1] = {7.51504679e+01};
    float bkg_XSEC[1] = {5.44570384e+01};

    float sig_NORM[1] = {1000000};
    float bkg_NORM[2] = {1207103};

    TFile *bkgSrc[1];

    bkgSrc[0] = TFile::Open("../wjets_600_800.root");

    TFile *sigSrc = TFile::Open("../signal_dtG_10.root");

    TTree *sigTree = (TTree *)sigSrc->Get("outtree");
    TTree *bkgTree[1];


    TFile *outf = new TFile("mva_smcedm.root", "RECREATE");
    TMVA::Factory *factory = new TMVA::Factory("factory_mva_smcedm", outf, "!V:!Silent:Color:DrawProgressBar:Transformations=I:AnalysisType=Classification");

    TMVA::DataLoader loader("dataset");

    loader.AddSignalTree(sigTree, sig_XSEC[0] / sig_NORM[0]);
    //loader.AddSignalTree(sigTree, sig_XSEC[0] / sig_NORM[0]);

    for (int k = 0; k < 1; k++) {
        bkgTree[k] = (TTree *)bkgSrc[k]->FindObjectAny("outtree");
        loader.AddBackgroundTree(bkgTree[k], bkg_XSEC[k] / bkg_NORM[k]);
    }

    const int NVAR = 15;
    TString VAR[NVAR] = {"br_njets", "br_scalar_HT",  "br_sphericity", "br_jet_pt[0]", "br_jet_pt[1]", "br_jet_pt[2]", "br_jet_pt[3]", "br_Fox_Wolfram[0]", "br_Fox_Wolfram[1]", "br_Fox_Wolfram[2]", "br_Fox_Wolfram[3]", "br_MET", "br_MET_Phi", "br_nleptons", "br_nbjets"};

    /*const int NVAR = 3;
    TString VAR[NVAR] = {"n_jets", "n_leptons", "n_bjets"};
    //char TYPE[NVAR] = {"I", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "F", "I", "F", "I"};*/

    for (int i = 0; i < NVAR; i++) {
        loader.AddVariable(VAR[i]);
    }

    TCut mycuts;
    loader.PrepareTrainingAndTestTree(mycuts, "nTrain_Signal=7000:nTrain_Background=8:nTest_Signal=1000:nTest_Background=3:SplitMode=Random:NormMode=NumEvents:!V");

    /*factory->BookMethod( &loader, TMVA::Types::kKNN, "KNN",
                         "H:nkNN=20:ScaleFrac=0.8:SigmaFact=1.0:CreateMVAPdfs:Kernel=Gaus:UseKernel=F:UseWeight=T:!Trim" );*/

    factory->BookMethod(&loader, TMVA::Types::kBDT, "BDT", "!V:NTrees=200:MinNodeSize=20.5%:MaxDepth=2:BoostType=AdaBoost:AdaBoostBeta=0.5:UseBaggedBoost:BaggedSampleFraction=0.5:SeparationType=GiniIndex:nCuts=20");

    factory->BookMethod( &loader, TMVA::Types::kMLP, "MLP", "H:!V:NeuronType=tanh:VarTransform=N:NCycles=600:HiddenLayers=N+5:TestRate=5:!UseRegulator" );

    // Use Deep Neural-Network

    //////General layout.

    TString layoutString ("Layout=TANH|128,TANH|128,TANH|128,LINEAR");

    //////Training strategies.
    TString training0("LearningRate=1e-1,Momentum=0.9,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.5+0.5+0.5, Multithreading=True");
    TString training1("LearningRate=1e-2,Momentum=0.9,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
                      "WeightDecay=1e-4,Regularization=L2,"
                      "DropConfig=0.0+0.0+0.0+0.0, Multithreading=True");
    TString training2("LearningRate=1e-3,Momentum=0.0,Repetitions=1,"
                      "ConvergenceSteps=20,BatchSize=256,TestRepetitions=10,"
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

    TString cpuOptions = dnnOptions + ":Architecture=GPU";

    factory->BookMethod( &loader, TMVA::Types::kDNN, "DNN CPU", cpuOptions);

    factory->TrainAllMethods();
    factory->TestAllMethods();
    factory->EvaluateAllMethods();
    outf->Close();
}
