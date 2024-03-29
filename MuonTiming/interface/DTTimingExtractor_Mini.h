#ifndef MuonIdentification_DTTimingExtractor_Mini_H
#define MuonIdentification_DTTimingExtractor_Mini_H

/**\class DTTimingExtractor_Mini
 *
 * Extracts timing information associated to a muon track
 *
*/
// Adapted from original code  from
// Original Author:  Traczyk Piotr
//         Created:  Thu Oct 11 15:01:28 CEST 2007
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h"

#include "DataFormats/Common/interface/Ref.h"

#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "Geometry/CommonDetUnit/interface/GlobalTrackingGeometry.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "RecoMuon/TrackingTools/interface/MuonSegmentMatcher.h"
#include "RecoMuon/MuonIdentification/interface/TimeMeasurementSequence.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"
#include "DataFormats/TrackReco/interface/TrackExtraFwd.h"

#include <vector>

namespace edm {
  class ParameterSet;
  class EventSetup;
  class InputTag;
}

class MuonServiceProxy;

class DTTimingExtractor_Mini {

public:
  
  /// Constructor
  DTTimingExtractor_Mini(const edm::ParameterSet&, MuonSegmentMatcher *segMatcher);
  
  /// Destructor
  ~DTTimingExtractor_Mini();

 class TimeMeasurement
  {
   public:
     bool isLeft;
     bool isPhi;
     float posInLayer;
     float distIP;
     float timeCorr;
     int station;
     DetId driftCell;
  };

 void fillTiming(TimeMeasurementSequence &tmSequence, 
		 const std::vector<const DTRecSegment4D*> &segments,
		 reco::TrackRef muonTrack,
		 const edm::Event& iEvent, const edm::EventSetup& iSetup);

 void fillTiming(TimeMeasurementSequence &tmSequence, reco::TrackRef muonTrack,
		 const edm::Event& iEvent, const edm::EventSetup& iSetup);

private:
  double fitT0(double &a, double &b, const std::vector<double>& xl, const std::vector<double>& yl, const std::vector<double>& xr, const std::vector<double>& yr );

  edm::InputTag DTSegmentTags_; 
  unsigned int theHitsMin_;
  double thePruneCut_;
  double theTimeOffset_;
  double theError_;
  bool useSegmentT0_;
  bool doWireCorr_;
  bool dropTheta_;
  bool requireBothProjections_;
  bool debug;
  
  std::unique_ptr<MuonServiceProxy> theService;
  MuonSegmentMatcher *theMatcher;

};

#endif
