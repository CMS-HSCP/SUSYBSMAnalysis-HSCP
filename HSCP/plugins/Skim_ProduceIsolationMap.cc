// -*- C++ -*-
//
// Package:    ProduceIsolationMap
// Class:      ProduceIsolationMap
//
/*\class ProduceIsolationMap ProduceIsolationMap.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Loic Quertenmont
//         Created:  Wed Nov 10 16:41:46 CDT 2010
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/TrackExtra.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHit.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/TrackAssociator/interface/TrackAssociatorParameters.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetMatchInfo.h"

#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/DetId/interface/DetIdCollection.h"

//#include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h"
#include "Geometry/Records/interface/CaloTopologyRecord.h"
#include "Geometry/CaloTopology/interface/CaloTopology.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/Math/interface/deltaR.h"
#include <iostream>

//
// class declaration
//

using namespace susybsm;
using namespace edm;

class ProduceIsolationMap : public edm::EDProducer {
public:
  explicit ProduceIsolationMap(const edm::ParameterSet&);
  ~ProduceIsolationMap();
  virtual void produce(edm::Event&, const edm::EventSetup&) override;

private:
  edm::EDGetTokenT<reco::TrackCollection> TKToken_;
  edm::EDGetTokenT<reco::TrackCollection> inputCollectionToken_;
  std::vector<double> TKIsolationPtcut_;
  std::vector<double> IsolationConeDR_;
  std::vector<std::string> Label_;
  double candMinPt_;
  TrackDetectorAssociator trackAssociator_;
  TrackAssociatorParameters parameters_;
  // ----------member data ---------------------------
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
ProduceIsolationMap::ProduceIsolationMap(const edm::ParameterSet& iConfig) {
  TKToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("TKLabel"));
  inputCollectionToken_ = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("inputCollection"));
  TKIsolationPtcut_ = iConfig.getParameter<std::vector<double>>("TkIsolationPtCut");
  IsolationConeDR_ = iConfig.getParameter<std::vector<double>>("IsolationConeDR");
  Label_ = iConfig.getParameter<std::vector<std::string>>("Label");
  candMinPt_ = iConfig.getParameter<double>("CandidateMinPt");
  if (TKIsolationPtcut_.size() != IsolationConeDR_.size() || TKIsolationPtcut_.size() != Label_.size()) {
    printf(
        "The size of the following vector must be equal in the HSCP Isolation producer: TkIsolationPtCut, "
        "IsolationConeDR, Label\nFix your configuration file\n");
    exit(0);
  }

  // TrackAssociator parameters
  edm::ParameterSet parameters = iConfig.getParameter<edm::ParameterSet>("TrackAssociatorParameters");
  edm::ConsumesCollector iC = consumesCollector();
  parameters_.loadParameters(parameters, iC);
  trackAssociator_.useDefaultPropagator();

  //register your products
  for (unsigned int i = 0; i < Label_.size(); i++) {
    produces<ValueMap<HSCPIsolation>>(Label_[i]);
  }
}

ProduceIsolationMap::~ProduceIsolationMap() {}

void ProduceIsolationMap::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;

  using reco::TrackCollection;

  Handle<TrackCollection> TKHandle;
  iEvent.getByToken(TKToken_, TKHandle);
  if (!TKHandle.isValid()) {
    edm::LogError("ProduceIsolationMap") << "TK Tracks collection not found";
    return;
  }

  //loop through tracks.
  Handle<TrackCollection> tkTracks;
  iEvent.getByToken(inputCollectionToken_, tkTracks);
  std::vector<std::vector<HSCPIsolation>> IsolationInfoColl(Label_.size());
  for (unsigned int i = 0; i < Label_.size(); i++) {
    IsolationInfoColl[i].resize(tkTracks->size());
  }

  int TkIndex = 0;
  for (TrackCollection::const_iterator itTrack = tkTracks->begin(); itTrack != tkTracks->end(); ++itTrack, TkIndex++) {
    std::vector<double> SumPt;
    for (unsigned int i = 0; i < Label_.size(); i++) {
      SumPt.push_back(0);
    }
    std::vector<double> Count;
    for (unsigned int i = 0; i < Label_.size(); i++) {
      Count.push_back(0);
    }
    std::vector<double> CountHighPt;
    for (unsigned int i = 0; i < Label_.size(); i++) {
      CountHighPt.push_back(0);
    }

    if (itTrack->pt() >= candMinPt_) {
      // version in 10_6_X
      //         TrackDetMatchInfo info = trackAssociator_.associate(iEvent, iSetup, trackAssociator_.getFreeTrajectoryState(iSetup, *itTrack), parameters_);
      // version in 12_1_0
      TrackDetMatchInfo info = trackAssociator_.associate(
          iEvent,
          iSetup,
          trackAssociator_.getFreeTrajectoryState(&iSetup.getData(parameters_.bFieldToken), *itTrack),
          parameters_);
      for (unsigned int i = 0; i < Label_.size(); i++) {
        if (info.ecalRecHits.size() > 0) {
          IsolationInfoColl[i][TkIndex].Set_ECAL_Energy(
              info.coneEnergy(IsolationConeDR_[i], TrackDetMatchInfo::EcalRecHits));
        }
        if (info.hcalRecHits.size() > 0) {
          IsolationInfoColl[i][TkIndex].Set_HCAL_Energy(
              info.coneEnergy(IsolationConeDR_[i], TrackDetMatchInfo::HcalRecHits));
        }
      }

      for (TrackCollection::const_iterator itTrack2 = TKHandle->begin(); itTrack2 != TKHandle->end(); ++itTrack2) {
        if (fabs(itTrack->pt() - itTrack2->pt()) < 0.1 && fabs(itTrack->eta() - itTrack2->eta()) < 0.05)
          continue;
        float dR = deltaR(itTrack->momentum(), itTrack2->momentum());
        for (unsigned int i = 0; i < Label_.size(); i++) {
          if (dR > IsolationConeDR_[i])
            continue;
          SumPt[i] += itTrack2->pt();
          Count[i]++;
          if (itTrack2->pt() >= TKIsolationPtcut_[i])
            CountHighPt[i]++;
        }
      }

      for (unsigned int i = 0; i < Label_.size(); i++) {
        IsolationInfoColl[i][TkIndex].Set_TK_CountHighPt(CountHighPt[i]);
        IsolationInfoColl[i][TkIndex].Set_TK_Count(Count[i]);
        IsolationInfoColl[i][TkIndex].Set_TK_SumEt(SumPt[i]);
      }
    } else {
      for (unsigned int i = 0; i < Label_.size(); i++) {
        IsolationInfoColl[i][TkIndex].Set_ECAL_Energy(-1);
        IsolationInfoColl[i][TkIndex].Set_HCAL_Energy(-1);
        IsolationInfoColl[i][TkIndex].Set_TK_CountHighPt(-1);
        IsolationInfoColl[i][TkIndex].Set_TK_Count(-1);
        IsolationInfoColl[i][TkIndex].Set_TK_SumEt(-1);
      }
    }
  }

  for (unsigned int i = 0; i < Label_.size(); i++) {
    //Create empty output collections
    unique_ptr<ValueMap<HSCPIsolation>> trackHSCPIsolMap(new ValueMap<HSCPIsolation>);
    ValueMap<HSCPIsolation>::Filler filler(*trackHSCPIsolMap);

    filler.insert(tkTracks, IsolationInfoColl[i].begin(), IsolationInfoColl[i].end());
    filler.fill();
    iEvent.put(std::move(trackHSCPIsolMap), Label_[i]);
  }
}
//define this as a plug-in
DEFINE_FWK_MODULE(ProduceIsolationMap);
