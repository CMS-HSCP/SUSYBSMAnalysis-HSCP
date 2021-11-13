// -*- C++ -*-
//
// Package:    HitAnalyzer/MuonFilter
// Class:      MuonFilter
//
/**\class MuonFilter MuonFilter.cc HitAnalyzer/MuonFilter/plugins/MuonFilter.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Elisabetta Manca
//         Created:  Mon, 21 May 2018 07:31:04 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/stream/EDFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/StreamID.h"

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "functions.h"

//
// class declaration
//

class MuonFilter : public edm::stream::EDFilter<>
{
public:
  explicit MuonFilter(const edm::ParameterSet &);
  ~MuonFilter();

  static void fillDescriptions(edm::ConfigurationDescriptions &descriptions);

private:
  virtual void beginStream(edm::StreamID) override;
  virtual bool filter(edm::Event &, const edm::EventSetup &) override;
  virtual void endStream() override;

  //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
  //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
  //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

  // ----------member data ---------------------------

  edm::EDGetTokenT<std::vector<reco::GenParticle>> GenParticlesToken_;
  edm::EDGetTokenT<std::vector<reco::Muon>> MuonToken_;
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
MuonFilter::MuonFilter(const edm::ParameterSet &iConfig)
{
  //now do what ever initialization is needed

  GenParticlesToken_ = consumes<std::vector<reco::GenParticle>>(edm::InputTag("genParticles"));
  MuonToken_ = consumes<std::vector<reco::Muon>>(edm::InputTag("muons"));
}

MuonFilter::~MuonFilter()
{

  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called on each new Event  ------------
bool MuonFilter::filter(edm::Event &iEvent, const edm::EventSetup &iSetup)
{
  using namespace edm;

  Handle<std::vector<reco::GenParticle>> genPartCollection;
  iEvent.getByToken(GenParticlesToken_, genPartCollection);

  auto genParticles = *genPartCollection.product();

  Handle<std::vector<reco::Muon>> muonCollection;
  iEvent.getByToken(MuonToken_, muonCollection);

  auto muons = *muonCollection.product();

  // loop over gen particles
  for (std::vector<reco::GenParticle>::const_iterator g = genParticles.begin(); g != genParticles.end(); ++g)
  {

    if (g->pt() < 3.)
      return false;

    // loop over reco muons

    for (std::vector<reco::Muon>::const_iterator m = muons.begin(); m != muons.end(); ++m)
    {

      //if(m->eta() < -2.3 || m->eta() > -2.1) return false;
      //if(m->phi() > -2.74) return false;

      if (!(m->isPFMuon()))
        return false;
      if (!(m->isGlobalMuon() || m->isTrackerMuon()))
        return false;

      float dR = deltaR(g->phi(), m->phi(), g->eta(), m->eta());
      if (dR < 0.15)
      {

        if (m->innerTrack()->normalizedChi2() > 1.5)
          return false;
        if (m->innerTrack()->numberOfValidHits() < 14)
          return false;
        if (m->innerTrack()->ptError() / m->innerTrack()->pt() > 0.1)
          return false;
      }
    }
  }

  return true;
}

// ------------ method called once each stream before processing any runs, lumis or events  ------------
void
    MuonFilter::beginStream(edm::StreamID)
{
}

// ------------ method called once each stream after processing all runs, lumis and events  ------------
void MuonFilter::endStream()
{
}

// ------------ method called when starting to processes a run  ------------
/*
void
MuonFilter::beginRun(edm::Run const&, edm::EventSetup const&)
{ 
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void
MuonFilter::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void
MuonFilter::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void
MuonFilter::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void MuonFilter::fillDescriptions(edm::ConfigurationDescriptions &descriptions)
{
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}
//define this as a plug-in
DEFINE_FWK_MODULE(MuonFilter);
