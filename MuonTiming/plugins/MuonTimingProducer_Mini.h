#ifndef MuonTiming_MuonTimingProducer_Mini_h
#define MuonTiming_MuonTimingProducer_Mini_h 1

// -*- C++ -*-
//
// Package:    MuonTimingProducer_Mini
// Class:      MuonTimingProducer_Mini
// 
/**\class MuonTimingProducer_Mini MuonTimingProducer_Mini.h RecoMuon/MuonIdentification/interface/MuonTimingProducer_Mini.h

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
// Adapted from original code  from
// Original Author:  Piotr Traczyk, CERN
//         Created:  Mon Mar 16 12:27:22 CET 2009
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDProducer.h"

#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "SUSYBSMAnalysis/MuonTiming/interface/MuonTimingFiller_Mini.h"


//
// class decleration
//

class MuonTimingProducer_Mini : public edm::stream::EDProducer<> {
   public:
      explicit MuonTimingProducer_Mini(const edm::ParameterSet&);
      ~MuonTimingProducer_Mini() override;

   private:
      void produce(edm::Event&, const edm::EventSetup&) override;
      
      // ----------member data ---------------------------
      edm::InputTag m_muonCollection;
      edm::EDGetTokenT<pat::MuonCollection> muonToken_;

      MuonTimingFiller_Mini* theTimingFiller_;

};

#endif
