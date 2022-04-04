
#include <exception>
#include <vector>

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
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"

#include "DataFormats/FWLite/interface/LuminosityBlock.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/Luminosity/interface/LumiSummary.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

using namespace fwlite;

#endif
void analyzeEDM() {
  //does nothing... this function is used to compile the code (without providing any argument)
}

void analyzeEDM(std::string fileName, std::string cwdPath = "./") {
  double Total = 0;
  double Skimmed = 0;
  printf("process %s\n", fileName.c_str());
  fflush(stdout);
  TFile* file = TFile::Open((cwdPath + "/FARM_SkimEff/outputs/" + fileName + ".root").c_str());
  if (!file or file->IsZombie())
    return;

  fwlite::LuminosityBlock ls(file);
  for (ls.toBegin(); !ls.atEnd(); ++ls) {
    fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
    nEventsTotalCounter.getByLabel(ls, "nEventsBefSkim");
    if (!nEventsTotalCounter.isValid()) {
      printf("Invalid nEventsTotalCounterH\n");
      continue;
    }
    Total += nEventsTotalCounter->value;

    fwlite::Handle<edm::MergeableCounter> nEventsSkimmedCounter;
    nEventsSkimmedCounter.getByLabel(ls, "nEventsBefEDM");
    if (!nEventsSkimmedCounter.isValid()) {
      printf("Invalid nEventsSkimmedCounterH\n");
      continue;
    }
    Skimmed += nEventsSkimmedCounter->value;
  }
  double Eff = Skimmed / Total;
  double Size = file->GetSize() / (Skimmed * 1024);

  system((std::string("mkdir -p ") + cwdPath + "/Results/").c_str());
  FILE* pFile = fopen((cwdPath + "/Results/" + fileName + ".txt").c_str(), "w");
  fprintf(pFile,
          "%60s : Skimmed/Total = %6.0f/%6.0f = %6.2f%%  EventSize=%6.2fKB\n",
          fileName.c_str(),
          Skimmed,
          Total,
          100.0 * Eff,
          Size);
  fclose(pFile);
}
