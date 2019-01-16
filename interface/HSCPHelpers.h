/*
  HSCP Helpers, a set of useful functions
  "C'e` modo e modo di dare una mano, questo lo sai anche tu"
*/

#ifndef HSCPHelpers_H
#define HSCPHelpers_H

#include "TH1F.h"
#include "TH3F.h"
#include "TMath.h"

#include "DataFormats/TrackReco/interface/HitPattern.h"

namespace hscphelpers {
  TH3F* preprocess_dEdx_template(TH3F* histo, bool splitByModuleType=true);
  int  muonStations(const reco::HitPattern& hitPattern);
}

#endif
