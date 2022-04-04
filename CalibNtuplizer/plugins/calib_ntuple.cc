// -*- C++ -*-
//
// Package:    stage/ntuple
// Class:      calib_ntuple
//
/**\class calib_ntuple calib_ntuple.cc stage/ntuple/plugins/calib_ntuple.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Collard Caroline
//         Created:  Mon, 14 Jan 2019 15:48:08 GMT
//
//

// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
//#include "DataFormats/SiPixelDetId/interface/PixelBarrelName.h"
#include "DataFormats/TrackerCommon/interface/PixelBarrelName.h"
#include "DataFormats/SiPixelDetId/interface/PixelBarrelNameUpgrade.h"
//#include "DataFormats/SiPixelDetId/interface/PixelEndcapName.h"
#include "DataFormats/TrackerCommon/interface/PixelEndcapName.h"
#include "DataFormats/SiPixelDetId/interface/PixelEndcapNameUpgrade.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"
//#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"

#include "DataFormats/Common/interface/ValueMap.h"
#include "DataFormats/TrackReco/interface/DeDxData.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "RecoTracker/DeDx/interface/DeDxTools.h"
#include "DataFormats/TrackReco/interface/TrackToTrackMap.h"

#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonReco/interface/Muon.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"

#include "SimDataFormats/CrossingFrame/interface/CrossingFrame.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"
#include "SimDataFormats/TrackerDigiSimLink/interface/StripDigiSimLink.h"
#include "SimDataFormats/TrackingHit/interface/PSimHitContainer.h"
//#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/Scalers/interface/LumiScalers.h"
#include "DataFormats/Luminosity/interface/LumiDetails.h"

#include "HSCP_codeFromAnalysis.h"

#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"

#include "TH1.h"
#include <TTree.h>
#include <string.h>
//
// class declaration
//

// If the analyzer does not use TFileService, please remove
// the template argument to the base class so the class inherits
// from  edm::one::EDAnalyzer<> and also remove the line from
// constructor "usesResource("TFileService");"
// This will improve performance in multithreaded jobs.
//

const int nMaxTrack = 10000;
const int nMaxDeDxH = 100000;
const int nMaxStrip = 100000;
const int nMaxStripprim = 100000;
const int nMaxSimHit = 100000;
const int nMaxGen = 1000;
const int nMaxMuon = 10000;
const int nMaxHSCP = 10000;

class calib_ntuple : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit calib_ntuple(const edm::ParameterSet&);
  ~calib_ntuple();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;

  // ----------member data ---------------------------
  //       edm::EDGetTokenT<reco::VertexCollection> m_primaryVertexTag;
  edm::EDGetTokenT<reco::TrackCollection> m_tracksTag;
  edm::EDGetTokenT<reco::DeDxHitInfoAss> m_dedxTag;
  edm::EDGetTokenT<edm::ValueMap<int>> m_dedxPrescaleTag;

  int printOut_;
  //

  TH3F* dEdxTemplatesUncorr = NULL;
  TH3F* dEdxTemplatesCorr = NULL;

  //       TH1D * histo;
  TTree* smalltree;
  int tree_runNumber;
  uint32_t tree_event;
  //       int      tree_npv;
  //       int      tree_ngoodpv;

  int tree_ntracks;
  float tree_track_pt[nMaxTrack];
  float tree_track_pterr[nMaxTrack];
  float tree_track_p[nMaxTrack];
  float tree_track_eta[nMaxTrack];
  float tree_track_phi[nMaxTrack];
  float tree_track_charge[nMaxTrack];
  float tree_track_chi2[nMaxTrack];
  int tree_track_nvalidhits[nMaxTrack];
  int tree_track_npixhits[nMaxTrack];
  int tree_track_missing[nMaxTrack];
  float tree_track_validfrac[nMaxTrack];
  float tree_track_validlast[nMaxTrack];
  int tree_track_qual[nMaxTrack];
  float tree_track_dz[nMaxTrack];
  float tree_track_dxy[nMaxTrack];
  int tree_track_index_hit[nMaxTrack];
  int tree_track_nhits[nMaxTrack];
  int tree_track_prescale[nMaxTrack];
  float tree_track_ih_ampl[nMaxTrack];
  float tree_track_ih_ampl_corr[nMaxTrack];

  int tree_dedxhits;
  uint32_t tree_dedx_detid[nMaxDeDxH];
  int tree_dedx_subdetid[nMaxDeDxH];
  int tree_dedx_modulgeom[nMaxDeDxH];
  float tree_dedx_charge[nMaxDeDxH];
  float tree_dedx_pathlength[nMaxDeDxH];
  float tree_dedx_posx[nMaxDeDxH];
  float tree_dedx_posy[nMaxDeDxH];
  float tree_dedx_posz[nMaxDeDxH];
  bool tree_dedx_strip[nMaxDeDxH];
  bool tree_dedx_pixel[nMaxDeDxH];
  bool tree_dedx_insideTkMod[nMaxDeDxH];

  int tree_sclus_index_strip[nMaxDeDxH];
  int tree_sclus_nstrip[nMaxDeDxH];
  int tree_sclus_firstsclus[nMaxDeDxH];
  float tree_sclus_barycenter[nMaxDeDxH];
  float tree_sclus_charge[nMaxDeDxH];
  float tree_sclus_errorclus[nMaxDeDxH];
  bool tree_sclus_ismerged[nMaxDeDxH];
  bool tree_sclus_sat254[nMaxDeDxH];
  bool tree_sclus_sat255[nMaxDeDxH];
  bool tree_sclus_shape[nMaxDeDxH];
  int tree_sclus_index_strip_corr[nMaxDeDxH];
  int tree_sclus_nstrip_corr[nMaxDeDxH];
  float tree_sclus_charge_corr[nMaxDeDxH];
  bool tree_sclus_clusclean[nMaxDeDxH];
  bool tree_sclus_clusclean2[nMaxDeDxH];

  int tree_nstrips;
  int tree_strip_ampl[nMaxStrip];
  int tree_nstrips_corr;
  int tree_strip_ampl_corr[nMaxStripprim];
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
calib_ntuple::calib_ntuple(const edm::ParameterSet& iConfig)
/*
 :
  trackTags_(iConfig.getUntrackedParameter<edm::InputTag>("tracks"))

{
*/
{
  //   m_primaryVertexTag   = consumes<reco::VertexCollection>(iConfig.getParameter<edm::InputTag>("primaryVertexColl"));
  m_tracksTag = consumes<reco::TrackCollection>(iConfig.getParameter<edm::InputTag>("tracks"));
  m_dedxTag = consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dedx"));
  m_dedxPrescaleTag = consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("dEdxHitInfoPrescale"));
  printOut_ = iConfig.getUntrackedParameter<int>("printOut");

  //now do what ever initialization is needed
  usesResource("TFileService");
  edm::Service<TFileService> fs;

  smalltree = fs->make<TTree>("ttree", "ttree");
  smalltree->Branch("runNumber", &tree_runNumber);
  smalltree->Branch("event", &tree_event);
  //   smalltree -> Branch ( "npv",     &tree_npv ) ;
  //   smalltree -> Branch ( "ngoodpv",     &tree_ngoodpv ) ;

  smalltree->Branch("ntracks", &tree_ntracks);
  smalltree->Branch("track_pt", tree_track_pt, "track_pt[ntracks]/F");
  smalltree->Branch("track_pterr", tree_track_pterr, "track_pterr[ntracks]/F");
  smalltree->Branch("track_p", tree_track_p, "track_p[ntracks]/F");
  smalltree->Branch("track_eta", tree_track_eta, "track_eta[ntracks]/F");
  smalltree->Branch("track_phi", tree_track_phi, "track_phi[ntracks]/F");
  smalltree->Branch("track_charge", tree_track_charge, "track_charge[ntracks]/F");
  smalltree->Branch("track_chi2", tree_track_chi2, "track_chi2[ntracks]/F");
  smalltree->Branch("track_nvalidhits", tree_track_nvalidhits, "track_nvalidhits[ntracks]/I");
  smalltree->Branch("track_npixhits", tree_track_npixhits, "track_npixhits[ntracks]/I");
  smalltree->Branch("track_missing", tree_track_missing, "track_missing[ntracks]/I");
  smalltree->Branch("track_validfraction", tree_track_validfrac, "track_validfraction[ntracks]/F");
  smalltree->Branch("track_validlast", tree_track_validlast, "track_validlast[ntracks]/F");
  smalltree->Branch("track_qual", tree_track_qual, "track_qual[ntracks]/I");
  smalltree->Branch("track_dz", tree_track_dz, "track_dz[ntracks]/F");
  smalltree->Branch("track_dxy", tree_track_dxy, "track_dxy[ntracks]/F");
  smalltree->Branch("track_index_hit", tree_track_index_hit, "track_index_hit[ntracks]/I");
  smalltree->Branch("track_nhits", tree_track_nhits, "track_nhits[ntracks]/I");
  smalltree->Branch("track_prescale", tree_track_prescale, "track_prescale[ntracks]/I");
  smalltree->Branch("track_ih_ampl", tree_track_ih_ampl, "track_ih_ampl[ntracks]/F");
  smalltree->Branch("track_ih_ampl_corr", tree_track_ih_ampl_corr, "track_ih_ampl_corr[ntracks]/F");

  smalltree->Branch("ndedxhits", &tree_dedxhits);
  smalltree->Branch("dedx_detid", tree_dedx_detid, "dedx_detid[ndedxhits]/i");
  smalltree->Branch("dedx_subdetid", tree_dedx_subdetid, "dedx_subdetid[ndedxhits]/I");
  smalltree->Branch("dedx_modulgeom", tree_dedx_modulgeom, "dedx_modulgeom[ndedxhits]/I");
  smalltree->Branch("dedx_charge", tree_dedx_charge, "dedx_charge[ndedxhits]/F");
  smalltree->Branch("dedx_pathlength", tree_dedx_pathlength, "dedx_pathlength[ndedxhits]/F");
  smalltree->Branch("dedx_posx", tree_dedx_posx, "dedx_posx[ndedxhits]/F");
  smalltree->Branch("dedx_posy", tree_dedx_posy, "dedx_posy[ndedxhits]/F");
  smalltree->Branch("dedx_posz", tree_dedx_posz, "dedx_posz[ndedxhits]/F");
  smalltree->Branch("dedx_isstrip", tree_dedx_strip, "dedx_isstrip[ndedxhits]/O");
  smalltree->Branch("dedx_ispixel", tree_dedx_pixel, "dedx_ispixel[ndedxhits]/O");
  smalltree->Branch("dedx_insideTkMod", tree_dedx_insideTkMod, "dedx_insideTkMod[ndedxhits]/O");

  smalltree->Branch("sclus_firstsclus", tree_sclus_firstsclus, "sclus_firstsclus[ndedxhits]/I");
  smalltree->Branch("sclus_barycenter", tree_sclus_barycenter, "sclus_barycenter[ndedxhits]/F");
  smalltree->Branch("sclus_charge", tree_sclus_charge, "sclus_charge[ndedxhits]/F");
  smalltree->Branch("sclus_errorclus", tree_sclus_errorclus, "sclus_errorclus[ndedxhits]/F");
  smalltree->Branch("sclus_ismerged", tree_sclus_ismerged, "sclus_ismerged[ndedxhits]/O");
  smalltree->Branch("sclus_index_strip", tree_sclus_index_strip, "sclus_index_strip[ndedxhits]/I");
  smalltree->Branch("sclus_nstrip", tree_sclus_nstrip, "sclus_nstrip[ndedxhits]/I");
  smalltree->Branch("sclus_sat254", tree_sclus_sat254, "sclus_sat254[ndedxhits]/O");
  smalltree->Branch("sclus_sat255", tree_sclus_sat255, "sclus_sat255[ndedxhits]/O");
  smalltree->Branch("sclus_shape", tree_sclus_shape, "sclus_shape[ndedxhits]/O");
  smalltree->Branch("sclus_index_strip_corr", tree_sclus_index_strip_corr, "sclus_index_strip_corr[ndedxhits]/I");
  smalltree->Branch("sclus_nstrip_corr", tree_sclus_nstrip_corr, "sclus_nstrip_corr[ndedxhits]/I");
  smalltree->Branch("sclus_charge_corr", tree_sclus_charge_corr, "sclus_charge_corr[ndedxhits]/F");
  smalltree->Branch("sclus_clusclean", tree_sclus_clusclean, "sclus_clusclean[ndedxhits]/O");
  smalltree->Branch("sclus_clusclean2", tree_sclus_clusclean2, "sclus_clusclean2[ndedxhits]/O");

  smalltree->Branch("nstrips", &tree_nstrips);
  smalltree->Branch("strip_ampl", tree_strip_ampl, "strip_ampl[nstrips]/I");
  smalltree->Branch("nstrips_corr", &tree_nstrips_corr);
  smalltree->Branch("strip_ampl_corr", tree_strip_ampl_corr, "strip_ampl_corr[nstrips_corr]/I");
}

calib_ntuple::~calib_ntuple() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called for each event  ------------
void calib_ntuple::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace reco;

  using reco::TrackCollection;

  EventID myEvId = iEvent.id();
  tree_runNumber = myEvId.run();
  tree_event = myEvId.event();

  /*
    edm::Handle<reco::VertexCollection> primaryVertex ;
    iEvent.getByToken(m_primaryVertexTag,primaryVertex);
    tree_npv = primaryVertex->size();
    const std::vector<reco::Vertex>& vertexColl = *primaryVertex;
    int index_pv=-1;
    int goodVerts=0;
    bool firstpvfound=false;
    for(unsigned int i=0;i<vertexColl.size();i++){
      if(vertexColl[i].isFake() || fabs(vertexColl[i].z())>24 || vertexColl[i].position().rho()>2 || vertexColl[i].ndof()<=4)continue; //only consider good vertex
      if(!firstpvfound) {
          firstpvfound=true;
          index_pv=i;
      }
      goodVerts++;
    }
    tree_ngoodpv = goodVerts;
*/

  edm::Handle<reco::TrackCollection> trackCollectionHandle;
  iEvent.getByToken(m_tracksTag, trackCollectionHandle);

  // inspired by /opt/sbg/data/safe1/cms/ccollard/HSCP/CMSSW_9_4_3/src/SUSYBSMAnalysis-HSCP/test/UsefulScripts/DeDxStudy/DeDxStudy
  edm::Handle<reco::DeDxHitInfoAss> dedxCollH;
  iEvent.getByToken(m_dedxTag, dedxCollH);
  if (!dedxCollH.isValid())
    std::cout << " access problem to dedxHitInfo collection " << std::endl;

  // inspired by https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc
  edm::Handle<edm::ValueMap<int>> dedxHitInfoPrescale;
  iEvent.getByToken(m_dedxPrescaleTag, dedxHitInfoPrescale);

  edm::ESHandle<TrackerTopology> tTopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(tTopoHandle);
  const TrackerTopology* tTopo = tTopoHandle.product();

  tree_ntracks = 0;
  tree_dedxhits = 0;
  tree_nstrips = 0;
  tree_nstrips_corr = 0;

  unsigned int ntrack_in_coll = 0;
  ntrack_in_coll = trackCollectionHandle->size();

  for (unsigned int c = 0; c < ntrack_in_coll; c++) {
    //      reco::TrackRef track;

    float pt_tr = 0;
    float pterr_tr = 0;
    float p_tr = 0;
    float eta_tr = 0;
    float phi_tr = 0;
    float charge_tr = 0;
    float chi2_tr = 0;
    int nh_tr = 0;
    int nhpix_tr = 0;
    int nmisstilL_tr = 0;
    float frac_tr = -1;
    float frac2_tr = -1;
    int qual_tr = 0;
    float dz_tr = -1000;
    float dxy_tr = -1000;
    const reco::DeDxHitInfo* dedxHits = nullptr;
    reco::DeDxHitInfoRef dedxHitsRef;

    reco::TrackRef track = reco::TrackRef(trackCollectionHandle.product(), c);
    if (printOut_ > 0)
      std::cout << " track with pT =  " << track->pt() << std::endl;
    pt_tr = track->pt();
    pterr_tr = track->ptError();
    p_tr = track->p();
    eta_tr = track->eta();
    phi_tr = track->phi();
    charge_tr = track->charge();
    chi2_tr = track->chi2() / track->ndof();
    nh_tr = track->numberOfValidHits();
    nhpix_tr = track->hitPattern().numberOfValidPixelHits();
    nmisstilL_tr = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
                   track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
    frac_tr = track->validFraction();
    frac2_tr = track->found() <= 0 ? -1 : track->found() / float(track->found() + nmisstilL_tr);
    qual_tr = track->qualityMask();  // >=2 = (2 meaning HighPurity tracks)
                                     /*
          if (index_pv>-1) {
          dz_tr=track->dz(vertexColl[index_pv].position());
          dxy_tr=track->dxy(vertexColl[index_pv].position());
          }
*/
    dedxHitsRef = dedxCollH->get(track.key());
    if (!dedxHitsRef.isNull())
      dedxHits = &(*dedxHitsRef);

    if (tree_ntracks < nMaxTrack) {
      tree_track_pt[tree_ntracks] = pt_tr;
      tree_track_pterr[tree_ntracks] = pterr_tr;
      tree_track_p[tree_ntracks] = p_tr;
      tree_track_eta[tree_ntracks] = eta_tr;
      tree_track_phi[tree_ntracks] = phi_tr;
      tree_track_charge[tree_ntracks] = charge_tr;
      tree_track_chi2[tree_ntracks] = chi2_tr;
      tree_track_nvalidhits[tree_ntracks] = nh_tr;
      tree_track_npixhits[tree_ntracks] = nhpix_tr;
      tree_track_missing[tree_ntracks] = nmisstilL_tr;
      tree_track_validfrac[tree_ntracks] = frac_tr;
      tree_track_validlast[tree_ntracks] = frac2_tr;
      tree_track_qual[tree_ntracks] = qual_tr;
      tree_track_dz[tree_ntracks] = dz_tr;
      tree_track_dxy[tree_ntracks] = dxy_tr;

      //load dEdx informations
      if (!dedxHits)
        continue;

      tree_track_prescale[tree_ntracks] = (*dedxHitInfoPrescale)[dedxHitsRef];

      //hit level dEdx information (only done for MIPs)

      if (printOut_ > 0)
        std::cout << " with " << dedxHits->size() << " dedxHits info " << std::endl;

      tree_track_nhits[tree_ntracks] = dedxHits->size();
      tree_track_index_hit[tree_ntracks] = tree_dedxhits;

      double dEdxSF[2] = {1., 1.};       // 0 : Strip SF, 1 : Pixel to Strip SF
      double dEdxSF_corr[2] = {1., 1.};  // 0 : Strip SF, 1 : Pixel to Strip SF

      // Hits, SF, Templates, usePixel, useClusterCleaning, reverseProb, useTrunc,  TrackerGains, useStrips, mustBeInside, MaxStripNOM,
      // correctFEDSat, CrossTalkInv, dropLowerValue, ErrorDeDx
      tree_track_ih_ampl[tree_ntracks] = 0;
      // correction inverseXtalk = 0 --> take the raw amplitudes of the cluster;
      reco::DeDxData dedxMObjTmp =
          computedEdx(dedxHits, dEdxSF, NULL, true, true, false, false, NULL, true, true, 99, false, 0, 0.15, NULL);
      reco::DeDxData* dedxMObj = dedxMObjTmp.numberOfMeasurements() > 0 ? &dedxMObjTmp : NULL;
      if (dedxMObj)
        tree_track_ih_ampl[tree_ntracks] = dedxMObj->dEdx();

      tree_track_ih_ampl_corr[tree_ntracks] = 0;
      // correction inverseXtalk = 1 --> modify the amplitudes based on xtalk for non-sat cluster + correct for saturation
      reco::DeDxData dedxMObjTmp2 = computedEdx(
          dedxHits, dEdxSF_corr, NULL, true, true, false, false, NULL, true, true, 99, false, 1, 0.15, NULL);
      reco::DeDxData* dedxMObj2 = dedxMObjTmp2.numberOfMeasurements() > 0 ? &dedxMObjTmp2 : NULL;
      if (dedxMObj2)
        tree_track_ih_ampl_corr[tree_ntracks] = dedxMObj2->dEdx();

      for (unsigned int h = 0; h < dedxHits->size(); h++) {
        if (tree_dedxhits < nMaxDeDxH) {
          DetId detid(dedxHits->detId(h));
          if (printOut_ > 0)
            std::cout << " DetId " << (uint32_t)detid << std::endl;
          tree_dedx_detid[tree_dedxhits] = (uint32_t)detid;
          tree_dedx_subdetid[tree_dedxhits] = detid.subdetId();
          if (detid.subdetId() < 3) {
            tree_dedx_modulgeom[tree_dedxhits] = 15;
            if (detid.subdetId() == 1) {
              if (printOut_ > -1) {
                std::cout << " GeomDetEnumerators::P1PXB " << static_cast<int>(GeomDetEnumerators::P1PXB) << std::endl;
                std::cout << " PixelSubdetector::PixelBarrel " << static_cast<int>(PixelSubdetector::PixelBarrel)
                          << std::endl;
                std::cout << "layer in pix barrel " << detid.subdetId();
                std::cout << " for detid " << (uint32_t)detid;
                //                     std::cout << " --> detid "<<  int((detid >>16)&0xF);  // 2016
                std::cout << " --> detid " << int((detid >> 20) & 0xF);  // 2017
                std::cout << " with tTopo " << tTopo->pxbLayer(detid.rawId()) << std::endl;

                //                     std::cout << "layer number " << PixelBarrelName(detid,1).layerName() << " name "  << PixelBarrelName(detid,1).name() << std::endl;
                //                     std::cout << "layer number " << PixelBarrelNameUpgrade(detid).layerName() << " name "  << PixelBarrelNameUpgrade(detid).name() << std::endl;
              }
            } else {
              if (printOut_ > -1) {
                std::cout << " GeomDetEnumerators::P1PXEC " << static_cast<int>(GeomDetEnumerators::P1PXEC)
                          << std::endl;
                std::cout << " PixelSubdetector::PixelEndcap " << static_cast<int>(PixelSubdetector::PixelEndcap)
                          << std::endl;
                std::cout << "disk in pix endcap " << detid.subdetId();
                std::cout << " for detid " << (uint32_t)detid;
                //                     std::cout << " --> detid "<<  int((detid >>16)&0xF);  // 2016
                std::cout << " --> detid " << int((detid >> 18) & 0xF);  // 2017
                std::cout << " side " << int((detid >> 23) & 0x3);
                std::cout << " with tTopo " << tTopo->pxfDisk(detid.rawId()) << std::endl;

                //                     std::cout << "disk number "<<  PixelEndcapName(detid,1).diskName() <<  " name "  << PixelEndcapName(detid,1).name() << std::endl;
                //                     std::cout << "disk number "<<  PixelEndcapNameUpgrade(detid).diskName() <<  " name "  << PixelEndcapNameUpgrade(detid).name() << std::endl;
              }
            }
          } else {
            SiStripDetId SSdetId(detid);
            //               tree_dedx_modulgeom[tree_dedxhits]=SSdetId.moduleGeometry();
            int moduleGeometry = 0;
            if (SSdetId.moduleGeometry() == SiStripModuleGeometry::IB1)
              moduleGeometry = 1;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::IB2)
              moduleGeometry = 2;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::OB1)
              moduleGeometry = 3;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::OB2)
              moduleGeometry = 4;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W1A)
              moduleGeometry = 5;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W2A)
              moduleGeometry = 6;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W3A)
              moduleGeometry = 7;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W1B)
              moduleGeometry = 8;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W2B)
              moduleGeometry = 9;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W3B)
              moduleGeometry = 10;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W4)
              moduleGeometry = 11;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W5)
              moduleGeometry = 12;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W6)
              moduleGeometry = 13;
            else if (SSdetId.moduleGeometry() == SiStripModuleGeometry::W7)
              moduleGeometry = 14;
            tree_dedx_modulgeom[tree_dedxhits] = moduleGeometry;
          }
          tree_dedx_insideTkMod[tree_dedxhits] =
              isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId() >= 3 ? dedxHits->stripCluster(h) : NULL);
          tree_dedx_charge[tree_dedxhits] = dedxHits->charge(h);
          tree_dedx_pathlength[tree_dedxhits] = dedxHits->pathlength(h);
          const LocalPoint position_hit = dedxHits->pos(h);
          tree_dedx_posx[tree_dedxhits] = position_hit.x();
          tree_dedx_posy[tree_dedxhits] = position_hit.y();
          tree_dedx_posz[tree_dedxhits] = position_hit.z();

          if (printOut_ > 0)
            std::cout << " charge of " << h << " = " << dedxHits->charge(h) << std::endl;
          if (printOut_ > 0)
            std::cout << " path length = " << dedxHits->pathlength(h) << std::endl;
          if (printOut_ > 0)
            std::cout << " ChargeOverPathlength " << dedxHits->charge(h) / dedxHits->pathlength(h) << std::endl;
          if (printOut_ > 0)
            std::cout << " DetId " << detid.subdetId() << std::endl;

          if (detid.subdetId() >= 3) {
            tree_dedx_strip[tree_dedxhits] = true;
            tree_dedx_pixel[tree_dedxhits] = false;

            tree_sclus_firstsclus[tree_dedxhits] = (int)dedxHits->stripCluster(h)->firstStrip();
            tree_sclus_barycenter[tree_dedxhits] = dedxHits->stripCluster(h)->barycenter();
            tree_sclus_charge[tree_dedxhits] = dedxHits->stripCluster(h)->charge();
            tree_sclus_errorclus[tree_dedxhits] = dedxHits->stripCluster(h)->getSplitClusterError();
            tree_sclus_ismerged[tree_dedxhits] = dedxHits->stripCluster(h)->isMerged();
            tree_sclus_index_strip[tree_dedxhits] = tree_nstrips;
            tree_sclus_index_strip_corr[tree_dedxhits] = tree_nstrips_corr;
            tree_sclus_charge_corr[tree_dedxhits] = 0;
            tree_sclus_sat254[tree_dedxhits] = false;
            tree_sclus_sat255[tree_dedxhits] = false;
            tree_sclus_shape[tree_dedxhits] = DeDxTools::shapeSelection(*(dedxHits->stripCluster(h)));

            //               std::vector <uint8_t> amplis = dedxHits->stripCluster(h)->amplitudes();
            //               std::vector <int> amps = convert(amplis);
            auto const& amplls = dedxHits->stripCluster(h)->amplitudes();
            std::vector<int> amps;
            for (unsigned int i = 0; i < amplls.size(); i++) {
              amps.push_back((int)amplls[i]);
            }

            if (printOut_ > 0)
              std::cout << " amps.size() " << amps.size() << std::endl;

            tree_sclus_nstrip[tree_dedxhits] = amps.size();
            for (unsigned int iclu = 0; iclu < amps.size(); iclu++) {
              if (printOut_ > 0)
                std::cout << "   amplitude (cluster" << iclu << ") = " << amps[iclu] << std::endl;
              if (amps[iclu] >= 254)
                tree_sclus_sat254[tree_dedxhits] = true;
              if (amps[iclu] == 255)
                tree_sclus_sat255[tree_dedxhits] = true;
              if (tree_nstrips < nMaxStrip) {
                tree_strip_ampl[tree_nstrips] = amps[iclu];
                tree_nstrips++;
              }  // end if MaxStrips
              else {
                std::cout << "Limit reached for tree_nstrips " << tree_nstrips << std::endl;
              }  // end if MaxStrips
            }    // end loop strips

            std::vector<int> ampsprim = CrossTalkInv(amps, 0.10, 0.04, true);
            tree_sclus_nstrip_corr[tree_dedxhits] = ampsprim.size();
            for (unsigned int iclu = 0; iclu < ampsprim.size(); iclu++) {
              if (printOut_ > 0)
                std::cout << " corrected amplitude (cluster" << iclu << ") = " << (int)ampsprim[iclu] << std::endl;
              tree_sclus_charge_corr[tree_dedxhits] += ampsprim[iclu];
              if (tree_nstrips_corr < nMaxStripprim) {
                tree_strip_ampl_corr[tree_nstrips_corr] = ampsprim[iclu];
                tree_nstrips_corr++;
              }  // end if nMaxStripprim
              else {
                std::cout << "Limit reached for tree_nstrips_corr " << tree_nstrips_corr << std::endl;
              }  // end if nMaxStripprim
            }    // end loop strips

            tree_sclus_clusclean[tree_dedxhits] = clusterCleaning(amps, 0);
            tree_sclus_clusclean2[tree_dedxhits] = clusterCleaning(ampsprim, 1);

          } else {
            tree_dedx_strip[tree_dedxhits] = false;
            tree_dedx_pixel[tree_dedxhits] = true;
            tree_sclus_firstsclus[tree_dedxhits] = -1;
            tree_sclus_barycenter[tree_dedxhits] = -1;
            tree_sclus_charge[tree_dedxhits] = -1;
            tree_sclus_errorclus[tree_dedxhits] = -1;
            tree_sclus_ismerged[tree_dedxhits] = 0;
            tree_sclus_index_strip[tree_dedxhits] = -1;
            tree_sclus_nstrip[tree_dedxhits] = 0;
            tree_sclus_sat254[tree_dedxhits] = 0;
            tree_sclus_sat255[tree_dedxhits] = 0;
            tree_sclus_shape[tree_dedxhits] = 0;
            tree_sclus_clusclean[tree_dedxhits] = 0;
            tree_sclus_charge_corr[tree_dedxhits] = -1;
            tree_sclus_nstrip_corr[tree_dedxhits] = 0;
            tree_sclus_clusclean2[tree_dedxhits] = 0;

          }  // end if detID
          tree_dedxhits++;

        }  //end if MaxDeDx
        else {
          std::cout << "Limit reached for tree_dedxhits " << tree_dedxhits << std::endl;
        }  // end if MaxDeDx
      }    // end loop dEdx
      tree_ntracks++;
    }  // end if MaxTracks
    else {
      std::cout << "Limit reached for tree_ntracks " << tree_ntracks << std::endl;
    }  // end if MaxTracks

  }  // end loop TrackCollection

  smalltree->Fill();

#ifdef THIS_IS_AN_EVENT_EXAMPLE
  Handle<ExampleData> pIn;
  iEvent.getByLabel("example", pIn);
#endif

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
  ESHandle<SetupData> pSetup;
  iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just before starting event loop  ------------
void calib_ntuple::beginJob() {}

// ------------ method called once each job just after ending the event loop  ------------
void calib_ntuple::endJob() {}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void calib_ntuple::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(calib_ntuple);
