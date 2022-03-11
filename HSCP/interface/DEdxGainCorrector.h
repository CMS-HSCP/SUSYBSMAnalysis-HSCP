#ifndef DEdxGainCorrector_h
#define DEdxGainCorrector_h

#include <string>
#include <map>
#include <unordered_map>
#include "TFile.h"
#include "TKey.h"
#include "TList.h"

class DEdxGainCorrector {
private:
  std::map<unsigned int, std::unordered_map<unsigned int, float> > TrackerGainsPerRuns;

public:
  std::unordered_map<unsigned int, float>* TrackerGains;
  DEdxGainCorrector() { TrackerGains = NULL; }
  ~DEdxGainCorrector() {}

  void setRun(unsigned int currentRun) {
    if (TrackerGainsPerRuns.size() <= 0) {
      TrackerGains = NULL;
      return;
    }
    std::map<unsigned int, std::unordered_map<unsigned int, float> >::iterator it,
        itPrev = TrackerGainsPerRuns.begin();
    for (it = TrackerGainsPerRuns.begin(); it != TrackerGainsPerRuns.end(); it++) {
      if (it->first > currentRun) {
        TrackerGains = &(itPrev->second);
        return;
      }  //runs are ordered, so the previous iterator correspond to our run
      itPrev = it;
    }
    TrackerGains = &(itPrev->second);  //just in case we go beyond the list of run for which we have a correciton
  }

  void LoadDeDxCalibration(std::string path) {
    TrackerGainsPerRuns.clear();
    TrackerGains = NULL;

    TFile* InputFile = new TFile(path.c_str(), "r");
    TList* ObjList = InputFile->GetListOfKeys();
    for (int i = 0; i < ObjList->GetSize(); i++) {
      TObject* tmp = ((TKey*)ObjList->At(i))->ReadObj();
      if (tmp->InheritsFrom("TTree")) {
        std::string dirName = ObjList->At(i)->GetName();
        unsigned int FirstRun, LastRun;
        sscanf(dirName.c_str(), "Gains_%d_to_%d", &FirstRun, &LastRun);
        printf("Add a new gain starting at run %d\n", FirstRun);

        TTree* t1 = (TTree*)tmp;
        unsigned int tree_DetId;
        t1->SetBranchAddress("DetId", &tree_DetId);
        unsigned char tree_APVId;
        t1->SetBranchAddress("APVId", &tree_APVId);
        float tree_Gain;
        t1->SetBranchAddress("Gain", &tree_Gain);
        //               float        tree_PrevGain;t1->SetBranchAddress("PrevGain"          ,&tree_PrevGain   );

        TrackerGains = &TrackerGainsPerRuns[FirstRun];
        for (unsigned int ientry = 0; ientry < t1->GetEntries(); ientry++) {
          t1->GetEntry(ientry);
          (*TrackerGains)[tree_DetId << 3 | tree_APVId] = tree_Gain;
        }
      }
    }
    InputFile->Close();
    delete InputFile;
  }
};

#endif
