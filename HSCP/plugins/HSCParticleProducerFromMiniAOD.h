// -*- C++ -*-
//
// Package:    HSCParticleProducerFromMiniAOD
// Class:      HSCParticleProducerFromMiniAOD
//
/**\class HSCParticleProducerFromMiniAOD HSCParticleProducer.h SUSYBSMAnalysis/HSCParticleProducer/interface/HSCParticleProducer.h

 Description: Producer for HSCP candidates, merging tracker dt information and rpc information

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Loic Quertenmont
//         Created:  Wed Oct 10 12:01:28 CEST 2007

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDFilter.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "CommonTools/UtilAlgos/interface/DeltaR.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "SUSYBSMAnalysis/HSCP/interface/BetaCalculatorTK.h"
#include "SUSYBSMAnalysis/HSCP/interface/BetaCalculatorMUON.h"
#include "SUSYBSMAnalysis/HSCP/interface/BetaCalculatorRPC.h"
#include "SUSYBSMAnalysis/HSCP/interface/BetaCalculatorECAL.h"
#include "SUSYBSMAnalysis/HSCP/interface/CandidateSelector.h"

#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "Math/GenVector/VectorUtil.h"

#include "DataFormats/PatCandidates/interface/Muon.h"

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include <TNtuple.h>
#include <TF1.h>

#include <vector>
#include <iostream>

//
// class decleration
//
class HSCParticleProducerFromMiniAOD : public edm::EDFilter {
public:
  explicit HSCParticleProducerFromMiniAOD(const edm::ParameterSet&);
  ~HSCParticleProducerFromMiniAOD();

private:
  virtual void beginJob();
  virtual bool filter(edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  std::vector<susybsm::HSCParticle> getHSCPSeedCollection(edm::Handle<reco::TrackCollection>& trackCollectionHandle,
                                                          edm::Handle<pat::MuonCollection>& muonCollectionHandle,
                                                          edm::Handle<pat::MuonCollection>& MTmuonCollectionHandle);

  // ----------member data ---------------------------
  bool Filter_;

  edm::EDGetTokenT<reco::TrackCollection> m_trackToken;
  edm::EDGetTokenT<reco::TrackCollection> m_trackIsoToken;
  edm::EDGetTokenT<pat::MuonCollection> m_muonsToken;
  edm::EDGetTokenT<pat::MuonCollection> m_MTmuonsToken;

  bool useBetaFromTk;
  bool useBetaFromMuon;
  bool useBetaFromRpc;
  bool useBetaFromEcal;

  float minTkP;
  float maxTkChi2;
  unsigned int minTkHits;
  float minMuP;
  float minSAMuPt;
  float minMTMuPt;
  float minDR;
  float minMTDR;
  float maxInvPtDiff;

  BetaCalculatorTK* beta_calculator_TK;
  BetaCalculatorMUON* beta_calculator_MUON;
  BetaCalculatorRPC* beta_calculator_RPC;
  BetaCalculatorECAL* beta_calculator_ECAL;

  std::vector<CandidateSelector*> Selectors;
};
