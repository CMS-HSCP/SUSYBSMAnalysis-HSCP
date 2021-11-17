#ifndef MuonIdentification_MuonTimingFiller_Mini_h
#define MuonIdentification_MuonTimingFiller_Mini_h 1

// -*- C++ -*-
//
// Package:    MuonTimingFiller_Mini
// Class:      MuonTimingFiller_Mini
// 
/**\class MuonTimingFiller_Mini MuonTimingFiller_Mini.h RecoMuon/MuonIdentification/interface/MuonTimingFiller_Mini.h

 Description: Class filling the DT, CSC and Combined MuonTimeExtra objects

 Implementation:
     <Notes on implementation>
*/
// Adapted from original  from
// Original Author:  Piotr Traczyk, CERN
//         Created:  Mon Mar 16 12:27:22 CET 2009
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

#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/PatCandidates/interface/Muon.h"
#include "SUSYBSMAnalysis/MuonTiming/interface/DTTimingExtractor_Mini.h"
#include "SUSYBSMAnalysis/MuonTiming/interface/CSCTimingExtractor_Mini.h"

//
// class decleration
//

class MuonTimingFiller_Mini {
   public:
      MuonTimingFiller_Mini(const edm::ParameterSet&, edm::ConsumesCollector&& iC);
      ~MuonTimingFiller_Mini();
      void fillTiming( const pat::Muon& muon, reco::MuonTimeExtra& dtTime, 
                    reco::MuonTimeExtra& cscTime, reco::MuonTime& rpcTime, 
                    reco::MuonTimeExtra& combinedTime, 
                    edm::Event& iEvent, const edm::EventSetup& iSetup );

   private:
      void fillTimeFromMeasurements( const TimeMeasurementSequence& tmSeq, reco::MuonTimeExtra &muTime );
      void fillRPCTime( const pat::Muon& muon, reco::MuonTime &muTime, edm::Event& iEvent );
      void rawFit(double &a, double &da, double &b, double &db, 
                  const std::vector<double>& hitsx, const std::vector<double>& hitsy);
      void addEcalTime( const pat::Muon& muon, TimeMeasurementSequence &cmbSeq );
      void combineTMSequences( const pat::Muon& muon, const TimeMeasurementSequence& dtSeq, 
                               const TimeMeasurementSequence& cscSeq, TimeMeasurementSequence &cmbSeq );
      
      std::unique_ptr<MuonSegmentMatcher> theMatcher_;
      std::unique_ptr<DTTimingExtractor_Mini> theDTTimingExtractor_Mini_;
      std::unique_ptr<CSCTimingExtractor_Mini> theCSCTimingExtractor_Mini_;
      double errorEB_,errorEE_,ecalEcut_;
      bool useDT_, useCSC_, useECAL_;

};

#endif
