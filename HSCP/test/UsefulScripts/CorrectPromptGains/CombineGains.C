#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TPaveText.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "../../AnalysisCode/tdrstyle.C"
#include "../../AnalysisCode/Analysis_CommonFunction.h"

#include <vector>
#include <string>
#include <iostream>
#include <fstream>

void Macro_Core(string input, string input2, string moduleName, string output);

string PathPrefix = "file:/afs/cern.ch/cms/tracker/sistrvalidation/WWW/CalibrationValidation/ParticleGain";
string PathToRoot = "sqlite/Gains.root";

void CombineGains() {
  ifstream fin("gains.txt", ifstream::in);
  string GainsName, GainsPath, PromptPath;
  system("rm -rf Gains && mkdir Gains");
  while (!fin.eof()) {
    fin >> GainsName >> GainsPath >> PromptPath;

    cout << "Processing: " << GainsName << " " << GainsPath << " " << PromptPath << endl;
    Macro_Core(PathPrefix + "/" + GainsPath + "/" + PathToRoot,   //first one is the true gains
               PathPrefix + "/" + PromptPath + "/" + PathToRoot,  //2nd one is the prompt gains
               "SiStripCalib",
               GainsName);
  }
  cerr << "About to commit a segfault! Noooo!" << endl;
}

void Macro_Core(string input, string input2, string moduleName, string GainsName) {
  unsigned int tree1_Index;
  unsigned int tree1_DetId;
  unsigned char tree1_APVId;
  double tree1_Gain;

  unsigned int tree2_Index;
  unsigned int tree2_DetId;
  unsigned char tree2_APVId;
  double tree2_Gain;

  TFile *f1 = new TFile(input.c_str(), "READ");
  TTree *t1 = (TTree *)GetObjectFromPath(f1, moduleName + "/APVGain");
  if (!t1) {
    cerr << "Tree " << input << " missing!" << endl;
    cerr << "Skipping!" << endl;
    return;
  }
  t1->SetBranchAddress("Index", &tree1_Index);
  t1->SetBranchAddress("DetId", &tree1_DetId);
  t1->SetBranchAddress("APVId", &tree1_APVId);
  t1->SetBranchAddress("Gain", &tree1_Gain);

  TFile *f2 = new TFile(input2.c_str(), "READ");
  TTree *t2 = (TTree *)GetObjectFromPath(f2, moduleName + "/APVGain");
  if (!t2) {
    cerr << "Tree " << input2 << " missing!" << endl;
    cerr << "Skipping!" << endl;
    return;
  }
  t2->SetBranchAddress("Index", &tree2_Index);
  t2->SetBranchAddress("DetId", &tree2_DetId);
  t2->SetBranchAddress("APVId", &tree2_APVId);
  t2->SetBranchAddress("Gain", &tree2_Gain);

  unsigned int tree3_Index;
  unsigned int tree3_DetId;
  unsigned char tree3_APVId;
  double tree3_Gain;

  TFile *f3 = new TFile(("Gains/" + GainsName + ".root").c_str(), "NEW");
  TTree *t3 = new TTree(GainsName.c_str(), GainsName.c_str());
  TBranch *branch_Index = t3->Branch("Index", &tree3_Index, "Index/i");
  TBranch *branch_DetId = t3->Branch("DetId", &tree3_DetId, "DetId/i");
  TBranch *branch_APVId = t3->Branch("APVId", &tree3_APVId, "APVId/b");
  TBranch *branch_Gain = t3->Branch("Gain", &tree3_Gain, "Gain/D");

  printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Looping on the Tree          :");
  int TreeStep = t1->GetEntries() / 50;
  if (TreeStep == 0)
    TreeStep = 1;
  for (size_t ientry = 0; ientry < t1->GetEntries(); ientry++) {
    if (ientry % TreeStep == 0) {
      printf(".");
      fflush(stdout);
    }
    t1->GetEntry(ientry);
    t2->GetEntry(ientry);

    tree3_Index = tree1_Index;
    tree3_DetId = tree1_DetId;
    tree3_APVId = tree1_APVId;
    tree3_Gain = tree1_Gain / tree2_Gain;
    t3->Fill();
  }
  printf("\n");
  f3 = t3->GetCurrentFile();
  f3->Write();
  f3->Close();
}
