
#include <exception>
#include <vector>
#include <string>

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
#include "TCutG.h"

// class stSignal;
// namespace edm {class TriggerResults; class TriggerResultsByName; class InputTag;}
// namespace reco { class Vertex; class Track; class GenParticle;}
// namespace susybsm {class HSCParticle;}
namespace fwlite {
  class ChainEvent;
}
// namespace trigger {class TriggerEvent;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

/*
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripMatchedRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"
#include "DataFormats/TrackerRecHit2D/interface/ProjectedSiStripRecHit2D.h"
#include "DataFormats/TrackerRecHit2D/interface/SiStripRecHit1D.h"
*/
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

using namespace fwlite;
// using namespace reco;
using namespace edm;
using namespace std;
// using namespace trigger;

#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"

#endif

void ReadXSecFromSample() {
  InitBaseDirectory();
  GetSampleDefinition(samples, "../../AnalysisCode/Analysis_Samples.txt");
  keepOnlyValidSamples(samples);

  for (unsigned int i = 0; i < samples.size(); i++) {
    if (samples[i].Type != 2)
      continue;  //only process signal samples

    vector<string> FileName;
    GetInputFiles(samples[i], BaseDirectory, FileName, 1);

    // loop on all the run blocks for an EDM file in order to count the number of events that are in a sample
    // this is useful to determine how to normalize the events (compute weight)
    double Total = 0;
    double Error = 0;
    int NRuns = 0;
    for (unsigned int f = 0; f < FileName.size(); f++) {
      TFile *file = TFile::Open(FileName[f].c_str());
      fwlite::Run run(file);
      for (run.toBegin(); !run.atEnd(); ++run) {
        fwlite::Handle<GenRunInfoProduct> genRunInfo;
        genRunInfo.getByLabel(run, "generator");
        if (!genRunInfo.isValid()) {
          printf("Invalid genRunInfo Handle\n");
          continue;
        }
        Total += genRunInfo->internalXSec().value();
        Error += genRunInfo->internalXSec().error();
        NRuns++;
      }
    }
    printf("%40s --> xsec= %8E +- %8E pb\n", samples[i].Name.c_str(), Total / NRuns, Error / NRuns);
  }
}
