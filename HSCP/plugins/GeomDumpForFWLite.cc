// -*- C++ -*-
//
// Package:    GeomDumpForFWLite
// Class:      GeomDumpForFWLite
//
/**\class GeomDumpForFWLite GeomDumpForFWLite.cc 

 Description: <one line class summary>

 Implementation:
     <Notes on implementation>
*/
//
// Original Author:  Loic QUERTENMONT
//         Created:  Fri Oct 26 07:22:12 CEST 2007 in the context of FROG
//        Modified:  Fri Sep 04 11:00:00 CEST 2015 to dump geometry in a root file for FWLite

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/CommonDetUnit/interface/GeomDet.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerNumberingBuilder/interface/GeometricDet.h"
#include "DataFormats/GeometrySurface/interface/BoundSurface.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "Geometry/CommonDetUnit/interface/TrackingGeometry.h"
#include "Geometry/Records/interface/MuonGeometryRecord.h"
#include "Geometry/DTGeometry/interface/DTGeometry.h"
#include "Geometry/CSCGeometry/interface/CSCGeometry.h"
#include "Geometry/RPCGeometry/interface/RPCGeometry.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "DataFormats/GeometrySurface/interface/TrapezoidalPlaneBounds.h"
#include "DataFormats/GeometrySurface/interface/RectangularPlaneBounds.h"

#include "DataFormats/MuonDetId/interface/DTLayerId.h"
#include "DataFormats/MuonDetId/interface/DTChamberId.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/RPCDetId.h"

#include "TTree.h"
#include "TVector3.h"

using namespace edm;
using namespace std;

//
// class decleration
//

class GeomDumpForFWLite : public edm::EDAnalyzer {
public:
  explicit GeomDumpForFWLite(const edm::ParameterSet&);
  ~GeomDumpForFWLite();

private:
  virtual void beginRun(edm::Run& run, const edm::EventSetup&);
  virtual void analyze(const edm::Event&, const edm::EventSetup&);
  virtual void endJob();

  bool isInitialized;

  // ----------member data ---------------------------

  //edm::ESHandle<TrackerGeometry> tkGeom;
  const edm::ESGetToken<TrackerGeometry,TrackerDigiGeometryRecord> tkGeomToken_;
  edm::ESHandle<DTGeometry> DtGeom;
  edm::ESHandle<CSCGeometry> CscGeom;
  edm::ESHandle<RPCGeometry> RpcGeom;
  //      edm::ESHandle<CaloGeometry>    CaloGeom;
};

//
// constructors and destructor
//
GeomDumpForFWLite::GeomDumpForFWLite(const edm::ParameterSet& iConfig) :
  tkGeomToken_(esConsumes<TrackerGeometry,TrackerDigiGeometryRecord>())
{ 
	isInitialized = false;
}

GeomDumpForFWLite::~GeomDumpForFWLite() {}

// ------------ method called once each job just before starting event loop  ------------
void GeomDumpForFWLite::beginRun(edm::Run& run, const edm::EventSetup& iSetup) {}

// ------------ method called once each job just after ending the event loop  ------------
void GeomDumpForFWLite::endJob() {}

// ------------ method called to for each event  ------------
void GeomDumpForFWLite::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  if (isInitialized)
    return;
  isInitialized = true;

  unsigned int rawId;
  float trapezeParam;
  TVector3* posV = new TVector3();
  TVector3* widthV = new TVector3();
  TVector3* lengthV = new TVector3();
  TVector3* thickV = new TVector3();

  edm::Service<TFileService> tfs;
  TTree* outtree = tfs->make<TTree>("geom", "geom");
  outtree->Branch("rawId", &rawId, "rawId/l");
  outtree->Branch("trapezeParam", &trapezeParam, "trapezeParam/f");
  outtree->Branch("pos", &posV, 32000, 0);
  outtree->Branch("width", &widthV, 32000, 0);
  outtree->Branch("length", &lengthV, 32000, 0);
  outtree->Branch("thick", &thickV, 32000, 0);

  //iSetup.get<TrackerDigiGeometryRecord>().get(tkGeom);
  const auto tkGeom = &iSetup.getData(tkGeomToken_);
  iSetup.get<MuonGeometryRecord>().get(DtGeom);
  iSetup.get<MuonGeometryRecord>().get(CscGeom);
  iSetup.get<MuonGeometryRecord>().get(RpcGeom);

  vector<const GeomDet*> TkDets = tkGeom->dets();
  vector<const GeomDet*> DtDets = DtGeom->dets();
  vector<const GeomDet*> CscDets = CscGeom->dets();
  vector<const GeomDet*> RpcDets = RpcGeom->dets();
  vector<const GeomDet*> MuonDets;
  for (unsigned int i = 0; i < TkDets.size(); i++) {
    MuonDets.push_back(TkDets[i]);
  }
  for (unsigned int i = 0; i < DtDets.size(); i++) {
    MuonDets.push_back(DtDets[i]);
  }
  for (unsigned int i = 0; i < CscDets.size(); i++) {
    MuonDets.push_back(CscDets[i]);
  }
  for (unsigned int i = 0; i < RpcDets.size(); i++) {
    MuonDets.push_back(RpcDets[i]);
  }

  for (unsigned int i = 0; i < MuonDets.size(); i++) {
    DetId Detid = MuonDets[i]->geographicalId();
    unsigned int SubDet = Detid.subdetId();
    if (Detid.det() == 1 && (SubDet < 1 || SubDet > 6))
      continue;

    const GeomDet* DetUnit = MuonDets[i];
    if (!DetUnit)
      continue;
    const BoundPlane plane = DetUnit->surface();
    const TrapezoidalPlaneBounds* trapezoidalBounds(dynamic_cast<const TrapezoidalPlaneBounds*>(&(plane.bounds())));
    const RectangularPlaneBounds* rectangularBounds(dynamic_cast<const RectangularPlaneBounds*>(&(plane.bounds())));

    rawId = Detid.rawId();
    trapezeParam = 0;

    float width = 0;
    float length = 0;
    float thickness = 0;
    if (trapezoidalBounds) {
      std::array<const float, 4> const& parameters = (*trapezoidalBounds).parameters();
      width = parameters[0] * 2;
      length = parameters[3] * 2;
      thickness = (*trapezoidalBounds).thickness();
      trapezeParam = parameters[1] / parameters[0];
    } else if (rectangularBounds) {
      width = DetUnit->surface().bounds().width();
      length = DetUnit->surface().bounds().length();
      thickness = DetUnit->surface().bounds().thickness();
      trapezeParam = 1;
    }

    Surface::GlobalPoint WidthVector = plane.toGlobal(LocalPoint(width / 2, 0, 0));
    Surface::GlobalPoint LengthVector = plane.toGlobal(LocalPoint(0, length / 2, 0));
    Surface::GlobalPoint ThickVector = plane.toGlobal(LocalPoint(0, 0, thickness / 2));
    GlobalVector Pos = GlobalVector(DetUnit->position().basicVector());

    posV->SetX(Pos.x());
    posV->SetY(Pos.y());
    posV->SetZ(Pos.z());
    widthV->SetX(WidthVector.x() - Pos.x());
    widthV->SetY(WidthVector.y() - Pos.y());
    widthV->SetZ(WidthVector.z() - Pos.z());
    lengthV->SetX(LengthVector.x() - Pos.x());
    lengthV->SetY(LengthVector.y() - Pos.y());
    lengthV->SetZ(LengthVector.z() - Pos.z());
    thickV->SetX(ThickVector.x() - Pos.x());
    thickV->SetY(ThickVector.y() - Pos.y());
    thickV->SetZ(ThickVector.z() - Pos.z());

    outtree->Fill();
  }
}

//define this as a plug-in
DEFINE_FWK_MODULE(GeomDumpForFWLite);
