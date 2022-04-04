#include <iostream>
#include <string>
#include "TTree.h"
#include "TFile.h"
#include "TObject.h"
#include "TList.h"
#include "../../AnalysisCode/Analysis_CommonFunction.h"

void ReorganizeGains() {
  TFile* fin = new TFile("Data13TeVGains_v2.root", "READ");
  TList* ObjList = fin->GetListOfKeys();
  for (int i = 0; i < ObjList->GetSize(); i++) {
    TTree* tmp = (TTree*)GetObjectFromPath(fin, ObjList->At(i)->GetName(), false);
    if (tmp->InheritsFrom("TTree")) {
      std::string name = std::string(tmp->GetName());
      std::string title = std::string(tmp->GetTitle());

      if (name.find("276870_to_999999") != std::string::npos && title.find("276870_to_999999") != std::string::npos) {
        name = "Gains_276870_to_278405";
      }

      TFile* fout = new TFile(("Gains/" + name + ".root").c_str(), "RECREATE");
      TTree* renamed = new TTree(name.c_str(), name.c_str());

      unsigned int tree1_Index;
      unsigned int tree1_DetId;
      unsigned char tree1_APVId;
      double tree1_Gain;

      unsigned int tree2_Index;
      unsigned int tree2_DetId;
      unsigned char tree2_APVId;
      double tree2_Gain;

      tmp->SetBranchAddress("Index", &tree1_Index);
      tmp->SetBranchAddress("DetId", &tree1_DetId);
      tmp->SetBranchAddress("APVId", &tree1_APVId);
      tmp->SetBranchAddress("Gain", &tree1_Gain);

      TBranch* branch_Index = renamed->Branch("Index", &tree2_Index, "Index/i");
      TBranch* branch_DetId = renamed->Branch("DetId", &tree2_DetId, "DetId/i");
      TBranch* branch_APVId = renamed->Branch("APVId", &tree2_APVId, "APVId/b");
      TBranch* branch_Gain = renamed->Branch("Gain", &tree2_Gain, "Gain/D");

      unsigned long int i = 0ul;
      for (; i < tmp->GetEntries(); i++) {
        tmp->GetEntry(i);

        tree2_Index = tree1_Index;
        tree2_DetId = tree1_DetId;
        tree2_APVId = tree1_APVId;
        tree2_Gain = tree1_Gain;
        renamed->Fill();
      }

      renamed->Write();
      delete renamed;
      delete tmp;
      fout->Close();
    }
  }
}
