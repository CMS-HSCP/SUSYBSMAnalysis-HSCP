// -*- C++ -*-
//
// Package:    SimHitShifterRun2
// Class:      SimHitShifterRun2
//
/**\class SimHitShifterRun2 SimHitShifterRun2.cc simhitshifter/SimHitShifterRun2/src/SimHitShifterRun2.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Camilo Andres Carrillo Montoya,40 2-B15,+41227671625,
//         Created:  Mon Aug 30 18:35:05 CEST 2010
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

#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/RPCDigi/interface/RPCDigi.h"
#include "DataFormats/RPCDigi/interface/RPCDigiCollection.h"
#include <DataFormats/RPCRecHit/interface/RPCRecHit.h>
#include "DataFormats/RPCRecHit/interface/RPCRecHitCollection.h"
#include <DataFormats/MuonDetId/interface/RPCDetId.h>
#include <Geometry/RPCGeometry/interface/RPCGeometry.h>
#include "Geometry/RPCGeometry/interface/RPCGeomServ.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "DataFormats/DetId/interface/DetId.h"
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandMatchMap.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "MagneticField/Engine/interface/MagneticField.h"
#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include <DataFormats/GeometrySurface/interface/LocalError.h>
#include <DataFormats/GeometryVector/interface/LocalPoint.h>
#include "DataFormats/GeometrySurface/interface/Surface.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"

#include "TrackingTools/GeomPropagators/interface/Propagator.h"
#include "TrackingTools/GeomPropagators/interface/AnalyticalPropagator.h"
#include "TrackPropagation/SteppingHelixPropagator/interface/SteppingHelixPropagator.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "TrackingTools/Records/interface/TrackingComponentsRecord.h"
#include "RecoMuon/TrackingTools/interface/MuonPatternRecoDumper.h"

#include "DataFormats/MuonDetId/interface/MuonSubdetId.h"

#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
#include "DataFormats/MuonDetId/interface/CSCIndexer.h"
#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"

#include "Geometry/RPCGeometry/interface/RPCGeometry.h"
#include <Geometry/RPCGeometry/interface/RPCRoll.h>
#include <Geometry/Records/interface/MuonGeometryRecord.h>
#include <Geometry/RPCGeometry/interface/RPCGeomServ.h>
#include <Geometry/CommonDetUnit/interface/GeomDet.h>

#include <cmath>

//Root
#include "TFile.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH1.h"
#include "TH2F.h"
#include "TROOT.h"
#include "TMath.h"
#include "TCanvas.h"

//Track
#include "DataFormats/TrajectorySeed/interface/TrajectorySeedCollection.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"

#include <fstream>

//
// class declaration
//

class SimHitShifterRun2 : public edm::EDProducer {
public:
  explicit SimHitShifterRun2(const edm::ParameterSet&);
  ~SimHitShifterRun2();
  //edm::ESHandle <RPCGeometry> rpcGeo;
  virtual void beginRun(const edm::Run&, const edm::EventSetup&) override;
  std::map<unsigned int, float> shiftinfo;

private:
  std::string ShiftFileName;
  virtual void beginJob(const edm::Run&, const edm::EventSetup&);
  virtual void produce(edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
};

SimHitShifterRun2::SimHitShifterRun2(const edm::ParameterSet& iConfig) {
  ShiftFileName =
      iConfig.getUntrackedParameter<std::string>("ShiftFileName",
                                                 "/afs/cern.ch/user/c/carrillo/simhits/CMSSW_3_5_8_patch2/src/"
                                                 "simhitshifter/SimHitShifterRun2/Merged_Muon_RawId_Shift.txt");

  //iSetup.get<MuonGeometryRecord>().get(rpcGeo);

  std::ifstream ifin(ShiftFileName.c_str());

  int rawId;
  float offset;

  if (!ifin)
    std::cout << "Problem reading the map rawId shift " << ShiftFileName.c_str() << std::endl;
  assert(ifin);

  while (ifin.good()) {
    ifin >> rawId >> offset;
    shiftinfo[rawId] = offset;
    //std::cout<<"Adding " << rawId << " to the map, with an offset = " << offset << "\n";
  }

  produces<edm::PSimHitContainer>("MuonCSCHits");
  produces<edm::PSimHitContainer>("MuonDTHits");
  produces<edm::PSimHitContainer>("MuonRPCHits");
}

SimHitShifterRun2::~SimHitShifterRun2() {}

void SimHitShifterRun2::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;

  //std::cout << " Getting the SimHits " <<std::endl;
  std::vector<edm::Handle<edm::PSimHitContainer> > theSimHitContainers;
  iEvent.getManyByType(theSimHitContainers);
  //std::cout << " The Number of sim Hits is  " << theSimHitContainers.size() <<std::endl;

  std::unique_ptr<edm::PSimHitContainer> pcsc(new edm::PSimHitContainer);
  std::unique_ptr<edm::PSimHitContainer> pdt(new edm::PSimHitContainer);
  std::unique_ptr<edm::PSimHitContainer> prpc(new edm::PSimHitContainer);

  std::vector<PSimHit> theSimHits;

  using std::dec;
  using std::oct;

  for (int i = 0; i < int(theSimHitContainers.size()); i++) {
    theSimHits.insert(theSimHits.end(), theSimHitContainers.at(i)->begin(), theSimHitContainers.at(i)->end());
  }

  for (std::vector<PSimHit>::const_iterator iHit = theSimHits.begin(); iHit != theSimHits.end(); iHit++) {
    DetId theDetUnitId((*iHit).detUnitId());
    DetId simdetid = DetId((*iHit).detUnitId());

    if (simdetid.det() != DetId::Muon)
      continue;

    float newtof = 0;

    if (simdetid.det() == DetId::Muon && simdetid.subdetId() == MuonSubdetId::RPC) {  //Only RPCs
      if (shiftinfo.find(simdetid.rawId()) == shiftinfo.end()) {
        std::cout << "RPC Warning the RawId = " << simdetid.det() << " | " << simdetid.rawId() << "is not in the map"
                  << std::endl;
        newtof = (*iHit).timeOfFlight();
      } else {
        newtof = (*iHit).timeOfFlight() + shiftinfo[simdetid.rawId()];
      }

      PSimHit hit((*iHit).entryPoint(),
                  (*iHit).exitPoint(),
                  (*iHit).pabs(),
                  newtof,
                  (*iHit).energyLoss(),
                  (*iHit).particleType(),
                  simdetid,
                  (*iHit).trackId(),
                  (*iHit).thetaAtEntry(),
                  (*iHit).phiAtEntry(),
                  (*iHit).processType());
      prpc->push_back(hit);
    } else if (simdetid.det() == DetId::Muon && simdetid.subdetId() == MuonSubdetId::DT) {  //Only DTs
      DTChamberId TheChamberDetId = DTChamberId(simdetid.rawId() & 0xFFC3FFFF);
      unsigned int RawId = TheChamberDetId.rawId();

      if (shiftinfo.find(RawId & 0xFFC3FFFF) == shiftinfo.end()) {
        std::cout << "DT Warning the RawId = " << RawId << "leading to condensed Id: " << (RawId & 0xFFC3FFFF)
                  << " is not in the map" << std::endl;
        newtof = (*iHit).timeOfFlight();
      } else {
        //std::cout<<"DT RawId = "<<RawId<<"leading to condensed Id: " << (RawId&0xFFC3FFFF)<< " is in the map"<<std::endl;
        newtof = (*iHit).timeOfFlight() + shiftinfo[RawId & 0xFFC3FFFF];
      }

      PSimHit hit((*iHit).entryPoint(),
                  (*iHit).exitPoint(),
                  (*iHit).pabs(),
                  newtof,
                  (*iHit).energyLoss(),
                  (*iHit).particleType(),
                  simdetid,
                  (*iHit).trackId(),
                  (*iHit).thetaAtEntry(),
                  (*iHit).phiAtEntry(),
                  (*iHit).processType());
      pdt->push_back(hit);
    } else if (simdetid.det() == DetId::Muon && simdetid.subdetId() == MuonSubdetId::CSC) {  //Only CSCs
      CSCDetId TheCSCDetId = CSCDetId(simdetid.rawId());
      CSCDetId TheChamberDetId = TheCSCDetId.chamberId();
      unsigned int RawId = TheChamberDetId.rawId();

      if (shiftinfo.find(RawId & 0xFFFFFE07) == shiftinfo.end()) {
        std::cout << "CSC Warning the RawId = " << RawId << "leading to condensed Id: " << (RawId & 0xFFFFFE07)
                  << " is not in the map" << std::endl;
        newtof = (*iHit).timeOfFlight();
      } else {
        //std::cout<<"CSC RawId = "<<RawId<<"leading to condensed Id: " << (RawId&0xFFFFFE07) << " is in the map"<<std::endl;
        newtof = (*iHit).timeOfFlight() + shiftinfo[RawId & 0xFFFFFE07];
      }
      PSimHit hit((*iHit).entryPoint(),
                  (*iHit).exitPoint(),
                  (*iHit).pabs(),
                  newtof,
                  (*iHit).energyLoss(),
                  (*iHit).particleType(),
                  simdetid,
                  (*iHit).trackId(),
                  (*iHit).thetaAtEntry(),
                  (*iHit).phiAtEntry(),
                  (*iHit).processType());
      pcsc->push_back(hit);
    }
  }

  iEvent.put(std::move(pcsc), "MuonCSCHits");
  iEvent.put(std::move(pdt), "MuonDTHits");
  iEvent.put(std::move(prpc), "MuonRPCHits");
}

void SimHitShifterRun2::beginRun(const edm::Run& run, const edm::EventSetup& iSetup) {}

// ------------ method called once each job just before starting event loop  ------------
void SimHitShifterRun2::beginJob(const edm::Run& run, const edm::EventSetup& iSetup) {}

// ------------ method called once each job just after ending the event loop  ------------
void SimHitShifterRun2::endJob() {}

//define this as a plug-in
DEFINE_FWK_MODULE(SimHitShifterRun2);
