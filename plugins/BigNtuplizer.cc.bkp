// -*- C++ -*-
//
// Class:      BigNtuplizer
// 
/**\class BigNtuplizer BigNtuplizer.cc 

   Description: Make a big ntuple from slimmed AOD for further use in the analysis

   Implementation:
   [Notes on implementation]
*/
//
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/Run.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "TTree.h"

#include <map>
#include <boost/algorithm/string.hpp>

#include "SUSYBSMAnalysis/HSCP/interface/BigNtuple.h"
#include "SUSYBSMAnalysis/HSCP/interface/CutFlowTracker.h"

#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPDeDxInfo.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "FWCore/Common/interface/TriggerNames.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "TH1F.h"
#include "TH3F.h"

//shady parts copy-pasted from the previous code
#include "SUSYBSMAnalysis/HSCP/interface/HSCPHelpers.h"
#include "SUSYBSMAnalysis/HSCP/interface/ModuleGeom.h"
#include "SUSYBSMAnalysis/HSCP/interface/MuonTimingCalculator.h"
#include "SUSYBSMAnalysis/HSCP/interface/DEdxGainCorrector.h"
#include "SUSYBSMAnalysis/HSCP/interface/L1BugEmulator.h"
#include "SUSYBSMAnalysis/HSCP/interface/DEdxHIPEmulator.h"
#include "SUSYBSMAnalysis/HSCP/interface/HIPTrackLossEmulator.h"

using namespace std;
using namespace edm;
using namespace reco;

//
// class declaration
//

class BigNtuplizer : public edm::one::EDAnalyzer<edm::one::SharedResources> {
public:
  explicit BigNtuplizer(const edm::ParameterSet&);
  ~BigNtuplizer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);
  float get_1d_value(TH1F* h, double x) { //returns the content of histogram at position x maybe better in a class
    int xbin = h->GetXaxis()->FindFixBin(x);
    return h->GetBinContent(xbin);
  }

private:
  virtual void beginRun(edm::Run& r) {
    unsigned int run = r.runAuxiliary().run();
    tofCalculator_.setRun(run);
    dEdxCorrector_.setRun(run);

    dEdxPixelToStrip_ = get_1d_value(h_dEdxPixelToStrip_, run);
    dEdxK_ = get_1d_value(h_dEdxK_, run);
    dEdxC_ = get_1d_value(h_dEdxC_, run);
    new_run_ = true;
  } 
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void initialize(const edm::Event&);
  virtual void endJob() override;

  typedef std::vector<std::string> vstring;


  // ----------member data ---------------------------
  //const edm::EDGetTokenT<std::vector<pat::Jet> > jets_;
  const edm::EDGetTokenT< std::vector<reco::Muon>  > muonsToken_;
  const edm::EDGetTokenT< std::vector<reco::Vertex>  > vertecesToken_;
  const edm::EDGetTokenT< std::vector<reco::Track>  > stAlMuonTracksToken_;
  const edm::EDGetTokenT< reco::BeamSpot > beamSpotToken_;
  const edm::EDGetTokenT< susybsm::HSCPIsolationValueMap  > isoHSCPToken_;
  const edm::EDGetTokenT< std::vector<reco::GenParticle>  > genParticlesToken_;
  const edm::EDGetTokenT< susybsm::HSCParticleCollection  > collectionHSCPToken_;
  const edm::EDGetTokenT< reco::DeDxHitInfoAss  > dedxHitToken_;
  const edm::EDGetTokenT< reco::MuonTimeExtraMap  > muonTOFToken_;
  const edm::EDGetTokenT< reco::MuonTimeExtraMap  > muonTOFDTToken_;
  const edm::EDGetTokenT< reco::MuonTimeExtraMap  > muonTOFCSCToken_;
  const edm::EDGetTokenT< susybsm::MuonSegmentCollection  > segCollectionCSCToken_;
  const edm::EDGetTokenT< susybsm::MuonSegmentCollection  > segCollectionDTToken_;
  const edm::EDGetTokenT< edm::TriggerResults > triggerResultsToken_;

  int typeMode_;
  edm::Service<TFileService> fs;
  TTree *tree_;
  BigNtuple ntuple_;
  CutFlowTracker tracker_;

protected:

  edm::Handle< std::vector<reco::Muon>  > muonsHandle_;
  edm::Handle< std::vector<reco::Vertex>  > vertecesHandle_;
  edm::Handle< std::vector<reco::Track>  > stAlMuonTracksHandle_;
  edm::Handle< reco::BeamSpot > beamSpotHandle_;
  edm::Handle< susybsm::HSCPIsolationValueMap  > isoHSCPHandle_;
  edm::Handle< std::vector<reco::GenParticle>  > genParticlesHandle_;
  edm::Handle< susybsm::HSCParticleCollection  > collectionHSCPHandle_;
  edm::Handle< reco::DeDxHitInfoAss  > dedxHitHandle_;
  edm::Handle< reco::MuonTimeExtraMap  > muonTOFHandle_;
  edm::Handle< reco::MuonTimeExtraMap  > muonTOFCSCHandle_;
  edm::Handle< reco::MuonTimeExtraMap  > muonTOFDTHandle_;
  edm::Handle< susybsm::MuonSegmentCollection  > segCollectionCSCHandle_;
  edm::Handle< susybsm::MuonSegmentCollection  > segCollectionDTHandle_;
  edm::Handle< edm::TriggerResults > triggerResultsHandle_;

  //pre-selection cuts
  const double        cutMaxEta_       =   2.1;    // cut on inner tracker track eta
  const double        cutMaxV3D_       =   99999;  //0.50;   // cut on 3D distance (cm) to closest vertex
  const double        cutMaxDZ_        =   0.50;   // cut on 1D distance (cm) to closest vertex in "Z" direction
  const double        cutMaxDXY_       =   0.50;   // cut on 2D distance (cm) to closest vertex in "R" direction
  const double        cutMaxChi2_      =   5.0;    // cut on Track maximal Chi2/NDF
  const int           cutMinQual_      =   2;      // cut on track quality (2 meaning HighPurity tracks)
  const unsigned int  cutMinNOH_       =   8;//7AMSB;      // cut on number of (valid) track pixel+strip hits 
  const int           cutMinNOPH_      =   2;      // cut on number of (valid) track pixel hits 
  const double        cutMinFOVH_      =   0.8;//0.0AMSB;    // cut on fraction of valid track hits
  const int           cutMaxNOMHTillLast_ = 99999;//1AMSB;     // cut on the number of missing hits from IP till last hit (excluding hits behind the last hit)
  const double        cutMinFOVHTillLast_ =-99999;//0.85AMSB;   // cut on the fraction of valid hits divided by total expected hits until the last one
  const unsigned int  cutMinNOM_       =   6;//7AMSB;      // cut on number of dEdx hits (generally equal to #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used)
  const double        cutMinNDOF_      =   8;      // cut on number of     DegreeOfFreedom used for muon TOF measurement
  const double        cutMinNDOFDT_    =   6;      // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
  const double        cutMinNDOFCSC_   =   6;      // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
  const double        cutMaxTOFErr_    =   0.15;//0.07;   // cut on error on muon TOF measurement
  const double        cutMaxPterr_     =   0.25;//0.50;//0.25;   // cut on error on track pT measurement 
  const double        cutMaxTIsol_     =  50;      // cut on tracker isolation (SumPt)
  const double        cutMaxRelTIsol_  =  9999999; // cut on relative tracker isolation (SumPt/Pt)
  const double        cutMaxEIsol_     =  0.30;    // cut on calorimeter isolation (E/P)
  const double        cutMinPt_        =  55.00;   // cut on pT    at PRE-SELECTION
  const double        cutMinIs_        =   0.0;    // cut on dEdxS at PRE-SELECTION (dEdxS is generally a  discriminator)
  const double        cutMinIm_        =   0.0;    // cut on dEdxM at PRE-SELECTION (dEdxM is generally an estimator    )
  const int           cutMinMuStations_=  2; 
  //double        cutMinTOF_       =   1.0;    // cut on TOF   at PRE-SELECTION

  //calibrations
  float dEdxSF_; //https://github.com/jozzez1/cmssw/blob/Run2HSCP16_v4/SUSYBSMAnalysis/HSCP/test/AnalysisCode_Eta12/Analysis_Step1_EventLoop.C#L110-L113
  //run dependent
  TH1F* h_dEdxPixelToStrip_;
  TH1F* h_dEdxK_;
  TH1F* h_dEdxC_;
  //caches values for the current run to save time
  float dEdxPixelToStrip_;
  float dEdxK_;
  float dEdxC_;

  TH3F* dEdxTemplates_;

  //other shit
  MuonTimingCalculator tofCalculator_;
  DEdxGainCorrector dEdxCorrector_;

  //
  std::vector<std::string> hlt_paths_;
  std::vector<long int> hlt_positions_;
  bool new_run_ = true;
};

BigNtuplizer::BigNtuplizer(const edm::ParameterSet& iConfig):
  muonsToken_(consumes<std::vector<reco::Muon>  >(iConfig.getParameter<edm::InputTag>("muons"))),  //muons
  vertecesToken_(consumes< std::vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVertices"))), //offlinePrimaryVertices
  stAlMuonTracksToken_(consumes< std::vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("stAlMuonTracks"))), //refittedStandAloneMuons
  beamSpotToken_(consumes<reco::BeamSpot >(iConfig.getParameter<edm::InputTag>("beamSpot"))), //offlineBeamSpot
  isoHSCPToken_(consumes< susybsm::HSCPIsolationValueMap  >(iConfig.getParameter<edm::InputTag>("isoHSCP"))), //"HSCPIsolation", "R03"
  genParticlesToken_(consumes< std::vector<reco::GenParticle> >(iConfig.getParameter<edm::InputTag>("genParticles"))),  // genParticlePlusGeant or genParticlesSkimmed or genParticles
  collectionHSCPToken_(consumes<susybsm::HSCParticleCollection >(iConfig.getParameter<edm::InputTag>("collectionHSCP"))),//HSCParticleProducer
  dedxHitToken_(consumes< reco::DeDxHitInfoAss >(iConfig.getParameter<edm::InputTag>("dedxHit"))),  //dedxHitInfo
  muonTOFToken_(consumes< reco::MuonTimeExtraMap >(iConfig.getParameter<edm::InputTag>("muonTOF"))),  // "muons", ????
  muonTOFDTToken_(consumes< reco::MuonTimeExtraMap >(iConfig.getParameter<edm::InputTag>("muonTOFDT"))), // "muons", ???
  muonTOFCSCToken_(consumes< reco::MuonTimeExtraMap >(iConfig.getParameter<edm::InputTag>("muonTOFCSC"))), // "muons", ???
  segCollectionCSCToken_(consumes<susybsm::MuonSegmentCollection >(iConfig.getParameter<edm::InputTag>("segCollectionCSC"))), // cscSegments
  segCollectionDTToken_(consumes< susybsm::MuonSegmentCollection >(iConfig.getParameter<edm::InputTag>("segCollectionDT"))), //  dtSegments
  triggerResultsToken_{consumes< edm::TriggerResults >(iConfig.getParameter<edm::InputTag>("HLT"))},
  typeMode_{iConfig.getParameter<int>("analysisType")},
  tree_(0),
  ntuple_(),
  tracker_(),
  dEdxSF_{iConfig.getParameter<float>("dEdxSF")},
  h_dEdxPixelToStrip_(0),
  h_dEdxK_(0),
  h_dEdxC_(0),
  dEdxPixelToStrip_(0),
  dEdxK_(0),
  dEdxC_(0),
  dEdxTemplates_(0),
  hlt_paths_{iConfig.getParameter<vector<string> >("HLTPaths")}
{
  usesResource("TFileService");
  tree_ = fs->make<TTree>("tree", "tree");
  ntuple_.set_evtinfo(tree_);

  //TODO Take from ROOT file
  edm::FileInPath calibrations(iConfig.getParameter<edm::FileInPath>("calibrations"));
  TFile tf(calibrations.fullPath().c_str());
  h_dEdxPixelToStrip_ = (TH1F*) ((TH1F*) tf.Get("dEdxPixelToStrip"))->Clone("dEdxPixelToStrip_"); 
  h_dEdxK_ = (TH1F*) ((TH1F*) tf.Get("dEdxK"))->Clone("dEdxK_");
  h_dEdxC_ = (TH1F*) ((TH1F*) tf.Get("dEdxC"))->Clone("dEdxC_");
  
  h_dEdxPixelToStrip_->SetDirectory(0);
  h_dEdxK_->SetDirectory(0);
  h_dEdxC_->SetDirectory(0);

  dEdxTemplates_ = hscphelpers::preprocess_dEdx_template((TH3F*) tf.Get("Charge_Vs_Path"));

  ModuleGeom::loadGeometry(iConfig.getParameter<edm::FileInPath>("geometry").fullPath()); //TODO should be replaced
  tofCalculator_.loadTimeOffset(iConfig.getParameter<edm::FileInPath>("TOFCalibration").fullPath()); //TODO check nothing is shady Is it a real calibration??
  dEdxCorrector_.TrackerGains = NULL; //This for some reason is always null o.O?
}



BigNtuplizer::~BigNtuplizer()
{

  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)

}


void BigNtuplizer::initialize(const edm::Event& iEvent){

  iEvent.getByToken(muonsToken_, muonsHandle_);
  iEvent.getByToken(vertecesToken_, vertecesHandle_);
  iEvent.getByToken(stAlMuonTracksToken_, stAlMuonTracksHandle_);
  iEvent.getByToken(beamSpotToken_, beamSpotHandle_);
  iEvent.getByToken(isoHSCPToken_, isoHSCPHandle_);
  iEvent.getByToken(genParticlesToken_, genParticlesHandle_);
  iEvent.getByToken(collectionHSCPToken_, collectionHSCPHandle_);
  iEvent.getByToken(dedxHitToken_, dedxHitHandle_);
  iEvent.getByToken(muonTOFToken_, muonTOFHandle_);
  iEvent.getByToken(muonTOFCSCToken_, muonTOFCSCHandle_);
  iEvent.getByToken(muonTOFDTToken_, muonTOFDTHandle_);
  iEvent.getByToken(segCollectionCSCToken_, segCollectionCSCHandle_);
  iEvent.getByToken(segCollectionDTToken_, segCollectionDTHandle_);
  iEvent.getByToken(triggerResultsToken_, triggerResultsHandle_);
}
//
// member functions
//

// ------------ method called for each event  ------------
void BigNtuplizer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  using namespace edm;
  initialize(iEvent);

  ntuple_.reset();
  ntuple_.fill_evtinfo(iEvent.id());
  tracker_.track("START");

  //Trigger
  //if new run refresh the trigger names info, somehow not stored in run but in the event
  if(new_run_) {
    const edm::TriggerNames& trigNames  = iEvent.triggerNames(*triggerResultsHandle_);
    hlt_positions_.clear();
    for(size_t i=0; i<hlt_paths_.size(); i++) hlt_positions_.push_back(-1); //reset everything to -1
    
    for (size_t i = 0; i < trigNames.size(); ++i) {
      const std::string &name = trigNames.triggerName(i);
      for(size_t j=0; j<hlt_paths_.size(); j++) {
	if(name.find(hlt_paths_.at(j)) != std::string::npos) {
	  hlt_positions_[j] = i;
	}
      }
    }
  }
  new_run_ = false;
  
  bool trig_pass = false;
  for(const long int idx : hlt_positions_) {
    trig_pass |= triggerResultsHandle_->accept(idx);    
  }
  if(!trig_pass) return;

  //loop over hscp candidates  ---line 1312 old code ---- 

  for (auto& hscpCand: *collectionHSCPHandle_) {

    susybsm::HSCParticle hscp  = hscpCand;
    reco::MuonRef  muon  = hscp.muonRef();   //this is needed for the 'muon only analysis'
    
    reco::TrackRef track;
    if(typeMode_!=3) track = hscp.trackRef();
    else {
      if(muon.isNull()) continue;
      track = muon->standAloneMuon();
    }
    if(track.isNull())continue;

    
    //TODO ---


    //pass preselection

    ////if(typeMode_==1 && !(hscp.type() == HSCParticleType::trackerMuon || hscp.type() == HSCParticleType::globalMuon))return false;
    //// if( (typeMode_==2 || typeMode_==4) && hscp.type() != HSCParticleType::globalMuon)return false;

    float GenBeta = -1; //FIXME
    ////Cut Flow: total number of events  st->Total->Fill(0.0,Event_Weight);   (utilizzato in Step3)
    // if(GenBeta>=0)st->Beta_Matched->Fill(GenBeta, Event_Weight);  //Control Plot
    // BS_Eta->Fill(track->eta(),Event_Weight);  
   
    if(fabs(track->eta())>cutMaxEta_) continue;

    int count = hscphelpers::muonStations(track->hitPattern());     ////function in Analysis_Step1_EventLoop.C
    // BS_MatchedStations->Fill(count, Event_Weight);      //Control Plot

    if(typeMode_==3 && count<cutMinMuStations_) continue;
    //// Cut Flow: st->Stations->Fill(0.0, Event_Weight);

    if(vertecesHandle_->size()<1) continue;  //si puo` fare? //const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;

    reco::Vertex bestVertex;
    double dzMin=10000;

    for (auto& vertex: *vertecesHandle_){

      if(vertex.isFake() || fabs(vertex.z())>24 || vertex.position().rho()>2 || vertex.ndof()<=4)continue;
      // BS_dzAll->Fill( track->dz (vertexColl[i].position()),Event_Weight); //Control Plot
      // BS_dxyAll->Fill(track->dxy(vertexColl[i].position()),Event_Weight); // Control Plot
      
      if(fabs(track->dz (vertex.position())) < fabs(dzMin) ){
	dzMin = fabs(track->dz (vertex.position()));
	bestVertex = vertex;
      }
    }
    
    // BS_NVertex->Fill(verteces.size(), Event_Weight);  //control plot _NoEventWeight for no weight one

    float dz = track->dz(bestVertex.position());
    float dxy = track->dxy(bestVertex.position());

    bool PUA = (vertecesHandle_->size()<15);
    bool PUB = (vertecesHandle_->size()>=15);

    // BS_TNOH_PUA->Fill(track->found(),Event_Weight);      //control plot
    // BS_TNOH_PUB->Fill(track->found(),Event_Weight);    //control plot
    // BS_TNOHFraction->Fill(track->validFraction(),Event_Weight);  //control plot
    // BS_TNOPH->Fill(track->hitPattern().numberOfValidPixelHits(),Event_Weight); //control plot
  
    if (typeMode_!=3 && (
			 track->found()<cutMinNOH_ || 
			 track->hitPattern().numberOfValidPixelHits()<cutMinNOPH_ || 
			 track->validFraction()<cutMinFOVH_ ))  continue;
  
    //this has to be fixed
    int missingHitsTillLast = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);;
    double validFractionTillLast = track->found()<=0?-1:track->found() / float(track->found() + missingHitsTillLast);
    ///////

    // BS_TNOHFractionTillLast->Fill(validFractionTillLast,Event_Weight);  //control plot
    // BS_TNOMHTillLast->Fill(missingHitsTillLast,Event_Weight); //control plot
    
    if (typeMode_!=3 && (
			 missingHitsTillLast>cutMaxNOMHTillLast_ ||
			 validFractionTillLast<cutMinFOVHTillLast_)) continue;

    ////Cut Flow:   st->TNOH  ->Fill(0.0,Event_Weight);  

    if (dedxSObj){
      // st->BS_TNOM->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);  //control plot 
      // if(PUA)st->BS_TNOM_PUA->Fill(dedxSObj->numberOfMeasurements(),Event_Weight); //control plot 
      // if(PUB)st->BS_TNOM_PUB->Fill(dedxSObj->numberOfMeasurements(),Event_Weight); //control plot 

      if (dedxSObj->numberOfMeasurements()<cutMinNOM_) continue;
    }
    //// Cut Flow:  st->TNOM  ->Fill(0.0,Event_Weight);   usato in Analysis_Step3_MakePlots.C

    if (tof){
      // st->BS_nDof->Fill(tof->nDof(),Event_Weight);  //control plot
      if((typeMode_>1  && typeMode_!=5) && tof->nDof()<cutMinNDOF_ && (dttof->nDof()<cutMinNDOFDT_ || csctof->nDof()<cutMinNDOFCSC_) ) continue;
    }
    //// Cut Flow: st->nDof  ->Fill(0.0,Event_Weight); usati in Analysis_Step3_MakePlots.C
    // st->BS_Qual->Fill(track->qualityMask(),Event_Weight); //control plot

    if(typeMode_!=3 && track->qualityMask()<cutMinQual_ ) continue;
    //// CutFlow: st->Qual  ->Fill(0.0,Event_Weight);   used in Analysis_Step3_MakePlots.C
    // st->BS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight); //control plot
    
    if(typeMode_!=3 && track->chi2()/track->ndof()>cutMaxChi2_ )continue;
    //// CutFlow:  st->Chi2  ->Fill(0.0,Event_Weight);  used in Analysis_Step3_MakePlots.C 

    // if(GenBeta>=0)st->Beta_PreselectedA->Fill(GenBeta, Event_Weight);  //control plot
    // st->BS_MPt ->Fill(track->pt(),Event_Weight); //control plot

    bool pass_trk_pt_central = track->pt() > cutMinPt_;   //base cut
    bool pass_trk_pt_shift = hscphelpers::rescaledPt(NVTrack.hitPattern()) > cutMinPt_;  // systematic selection cut

    if(!(pass_trk_pt_central || pass_trk_pt_shift))  continue;

    //if (pass_trk_pt_central)      st->MPt   ->Fill(0.0,Event_Weight);  //cut flow
    // if(dedxSObj && pass_trk_pt_central) st->BS_MIs->Fill(dedxSObj->dEdx(),Event_Weight); //control plot
    // if(dedxMObj && pass_trk_pt_central) st->BS_MIm->Fill(dedxMObj->dEdx(),Event_Weight); //control plot
    
    ///per l'analisi: dovro` applicare:  track->pt() > cutMinPt_;

    bool pass_dEdX = 1;
    bool pass_dEdX_rescale = 0;
    if (dedxSObj){
      pass_dEdX = dedxSObj->dEdx()>cutMinIs_;
      pass_dEdX_rescale = dedxSObj->dEdx()+RescaleI>cutMinIs_;
      if (!(pass_dEdX ||pass_dEdX_rescale )) continue;
    }      

    ///per l'analisi: dovro` applicare: pass_dEdX= dedxSObj->dEdx()>cutMinIs_

    if(dedxMObj && ((typeMode_!=5 && dedxMObj->dEdx()<cutMinIm_) || (typeMode_==5 && dedxMObj->dEdx()>cutMinIm_)) ) continue; 
    //// CutFlow:  if(pass_trk_pt_central && pass_dEdX) {st->MI   ->Fill(0.0,Event_Weight);}  ///used after I think....
    
    if(tof){
      // if(pass_trk_pt_central && pass_dEdX){st->BS_MTOF ->Fill(tof->inverseBeta(),Event_Weight);} //control plot 
      // if(pass_trk_pt_central && pass_dEdX)st->BS_TOFError->Fill(tof->inverseBetaErr(),Event_Weight);//control plot
      if(typeMode_>1 && typeMode_!=5  && tof->inverseBetaErr()>cutMaxTOFErr_) continue;   

      // if(pass_trk_pt_central && pass_dEdX) st->BS_TimeAtIP->Fill(tof->timeAtIpInOut(),Event_Weight); //control plot 
      if(typeMode_==3 && min(min(fabs(tof->timeAtIpInOut()-100), fabs(tof->timeAtIpInOut()-50)), min(fabs(tof->timeAtIpInOut()+100), fabs(tof->timeAtIpInOut()+50)))<5) continue;
    }
    
    // if(pass_trk_pt_central && pass_dEdX) st->BS_dzMinv3d->Fill(dz,Event_Weight);  // control plot
    // if(pass_trk_pt_central && pass_dEdX) st->BS_dxyMinv3d->Fill(dxy,Event_Weight);  // control plot
    // if(pass_trk_pt_central && pass_dEdX) st->BS_PV->Fill(goodVerts,Event_Weight); // control plot
    // if(pass_trk_pt_central && pass_dEdX) st->BS_PV_NoEventWeight->Fill(goodVerts); //control plot
    // if(pass_trk_pt_central && pass_dEdX && dedxSObj) st->BS_NOMoNOHvsPV->Fill(goodVerts,dedxSObj->numberOfMeasurements()/(double)track->found(),Event_Weight); //control plot


    //Require at least one good vertex except if cosmic event
    //  if(typeMode_==3 && goodVerts<1 && (!st || st->Name.find("Cosmic")==string::npos)) return false;
    if(typeMode_==3 && goodVerts<1) continue;   ////put exception for cosmics

    //For TOF only analysis match to a SA track without vertex constraint for IP cuts                 

    if(typeMode_==3) {
      //Find closest NV track                                                                                                                    
      reco::Track NVTrack;
      double minDr=15;
      for(auto& stAlMuon: *stAlMuonTracksHandle_){
	double dR = deltaR(track->eta(), track->phi(), stAlMuon.eta(), stAlMuon.phi());
	if(dR<minDr) {minDr=dR;
	  NVTrack=stAlMuon;}
      }
      // if(pass_trk_pt_central && pass_dEdX) st->BS_dR_NVTrack->Fill(minDr,Event_Weight); //control plot
      if(minDr>0.4) continue;
      //// CutFlow: if(pass_trk_pt_central && pass_dEdX)st->NVTrack->Fill(0.0,Event_Weight);
      
      //Find displacement of tracks with respect to beam spot                                                                                                                
      dz  = NVTrack.dz (beamSpot.position());   ///overwrite the value of dz and dxy
      dxy = NVTrack.dxy(beamSpot.position());
      if(hscphelpers::muonStations(NVTrack.hitPattern())<cutMinMuStations_) continue;
    }

    //// CutFlow: if(pass_trk_pt_central && pass_dEdX){st->MTOF ->Fill(0.0,Event_Weight);
    // if(GenBeta>=0 && pass_trk_pt_central && pass_dEdX)st->Beta_PreselectedB->Fill(GenBeta, Event_Weight);  //control plot
  
    double v3d = sqrt(dz*dz+dxy*dxy);

    // if(pass_trk_pt_central && pass_dEdX){st->BS_V3D->Fill(v3d,Event_Weight);}
    if(v3d>cutMaxV3D_ )continue;
    ////Cut Flow: if(pass_trk_pt_central && pass_dEdX){st->V3D  ->Fill(0.0,Event_Weight);}
    // if(pass_trk_pt_central && pass_dEdX)st->BS_Dxy->Fill(dxy, Event_Weight);

    ////TreeDXY = dxy;  magari dxy mi serve nelle ntuple??
    bool DXYSB = false;
    if(typeMode_!=5 && fabs(dxy)>cutMaxDXY_) continue;
    if(typeMode_==5 && fabs(dxy)>4) continue;
    if(typeMode_==5 && fabs(dxy)>cutMaxDXY_) DXYSB = true;

    //// Cut Flow: if(pass_trk_pt_central && pass_dEdX){st->Dxy  ->Fill(0.0,Event_Weight);}

    if(typeMode_!=3) {
      const ValueMap<HSCPIsolation>& IsolationMap = *isoHSCP.product();
      HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
      // if(pass_trk_pt_central && pass_dEdX){st->BS_TIsol ->Fill(hscpIso.Get_TK_SumEt(),Event_Weight);} //control plot

      if(hscpIso.Get_TK_SumEt()>cutMaxTIsol_) continue;
      ////Cut Flow: if(pass_trk_pt_central && pass_dEdX){st->TIsol   ->Fill(0.0,Event_Weight);}

      double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
      // if(pass_trk_pt_central && pass_dEdX){st->BS_EIsol ->Fill(EoP,Event_Weight);} // control plot

      if(EoP>cutMaxEIsol_)continue;
      /////CutFlow: if(pass_trk_pt_central && pass_dEdX){st->EIsol   ->Fill(0.0,Event_Weight);}

      // relative tracker isolation                                                                                          
      // if (pass_trk_pt_central && pass_dEdX) {  st->BS_SumpTOverpT->Fill(hscpIso.Get_TK_SumEt()/track->pt(), Event_Weight); }
      if(hscpIso.Get_TK_SumEt()/track->pt()>cutMaxRelTIsol_)continue;
      ////CutFlow : if (pass_trk_pt_central && pass_dEdX) {  st->SumpTOverpT   ->Fill(0.0,Event_Weight);}
    }

    ////Plot: if(pass_trk_pt_central && pass_dEdX){st->BS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);}
    if(typeMode_!=3 && (track->ptError()/track->pt())>cutMaxPterr_ && std::max(0.0,track->pt())<cutMinPt_)continue;
    ////CutFlow:  if(pass_trk_pt_central && pass_dEdX){st->Pterr   ->Fill(0.0,Event_Weight);}

    //Find distance to nearest segment on opposite side of detector                                                                                           
    double minPhi, minEta;
    double segSep=SegSep(hscp, ev, minPhi, minEta);  //check function
   
    
    ////Plots
    // if(pass_trk_pt_central && pass_dEdX){
    //   st->BS_SegSep->Fill(segSep, Event_Weight);  //control plot
    //   st->BS_SegMinPhiSep->Fill(minPhi, Event_Weight); //control plot 
    //   st->BS_SegMinEtaSep->Fill(minEta, Event_Weight); //control plot 
    //   //Plotting segment separation depending on whether track passed dz cut                                                                                   
    //   if(fabs(dz)>cutMaxDZ_) {
    // 	st->BS_SegMinEtaSep_FailDz->Fill(minEta, Event_Weight); //control plot 
    //   }
    //   else {
    // 	st->BS_SegMinEtaSep_PassDz->Fill(minEta, Event_Weight); //control plot 
    //   }
    //   //Plots for tracking failing Eta Sep cut                                                                                                                         
    //   if(fabs(minEta)<minSegEtaSep) {
    // 	//Needed to compare dz distribution of cosmics in pure cosmic and main sample                                                                                 
    // 	st->BS_Dz_FailSep->Fill(dz); //control plot 
    //   }
    // }
    
    //Now cut Eta separation                                                                                                                                            
    ////Cut Flow:  if(pass_trk_pt_central && pass_dEdX){st->SegSep->Fill(0.0,Event_Weight);}
    

    ////Plots
    // if(pass_trk_pt_central && pass_dEdX) {
    //   //Plots for tracks in dz control region                                                                                                              
    //   if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz && !muon->isGlobalMuon()) {
    // 	st->BS_Pt_FailDz->Fill(track->pt(), Event_Weight);  //control plot
    // 	st->BS_TOF_FailDz->Fill(tof->inverseBeta(), Event_Weight); //control plot 
    // 	if(fabs(track->eta())>CSCRegion) {
    // 	  st->BS_TOF_FailDz_CSC->Fill(tof->inverseBeta(), Event_Weight); //control plot 
    // 	  st->BS_Pt_FailDz_CSC->Fill(track->pt(), Event_Weight);   //control plot 
    // 	}
    // 	else if(fabs(track->eta())<DTRegion) {
    // 	  st->BS_TOF_FailDz_DT->Fill(tof->inverseBeta(), Event_Weight);  //control plot 
    // 	  st->BS_Pt_FailDz_DT->Fill(track->pt(), Event_Weight);      //control plot 
    // 	}
    //   }
    //   //Plots of dz                                                                                                                                
    //   st->BS_Dz->Fill(dz, Event_Weight);//control plot 
    //   if(fabs(track->eta())>CSCRegion) st->BS_Dz_CSC->Fill(dz,Event_Weight);  //control plot 
    //   else if(fabs(track->eta())<DTRegion) st->BS_Dz_DT->Fill(dz,Event_Weight);//control plot 
    //   st->BS_EtaDz->Fill(track->eta(),dz,Event_Weight); //control plot 
    // }

    //At this poitn is when it does: st->H_D_DzSidebands->Fill(CutIndex, DzType);    

    // se faccio il fill qua??

    bool DZSB = false;
    if(typeMode_!=5 && fabs(dz)>cutMaxDZ_) continue;
    if(typeMode_==5 && fabs(dz)>4) continue ;
    if(typeMode_==5 && fabs(dz)>cutMaxDZ_) DZSB = true;
    ////CutFlow :  if(pass_trk_pt_central && pass_dEdX){st->Dz  ->Fill(0.0,Event_Weight);}

    if(typeMode_==3 && fabs(minEta)<minSegEtaSep) continue;
    ////Plot: if(pass_trk_pt_central && pass_dEdX)st->BS_Phi->Fill(track->phi(),Event_Weight);
    if(typeMode_==3 && fabs(track->phi())>1.2 && fabs(track->phi())<1.9) continue;
    if(typeMode_==5 && fabs(dz)>4) continue;
    if(typeMode_==5 && fabs(dz)>cutMaxDZ_) DZSB = true;

    //skip HSCP that are compatible with cosmics.                                                                                                                            
    ////Plot: if(pass_trk_pt_central && pass_dEdX)st->BS_OpenAngle->Fill(OpenAngle,Event_Weight);

    bool OASB = false;
    if(typeMode_==5 && OpenAngle>=2.8)OASB = true;
    isCosmicSB = DXYSB && DZSB && OASB;
    isSemiCosmicSB = (!isCosmicSB && (DXYSB || DZSB || OASB));


    TVector3 outerHit = getOuterHitPos(dedxHits);  ////to be checked                                                                                                    
    TVector3 vertex(bestVertex.position().x(), bestVertex.position().y(), bestVertex.position().z());

    //dxy, dz,  phi, eta, p, pt, inverseBeta, dEdx s, dttof->inverseBeta(), csctof->inverseBeta(), tof->inverseBeta(),

    
    //Muon only prediction binned depending on where in the detector the track is and how many muon stations it has                                          
    //Binning not used for other analyses                                                                                                                      
    //Rimosso, non so se serve
   
    //// Cut Flow :if(pass_trk_pt_central && pass_dEdX){st->Basic  ->Fill(0.0,Event_Weight);}
  

  }
}
// ------------ method called once each job just before starting event loop  ------------
void 
BigNtuplizer::beginJob()
{
}

// ------------ method called once each job just after ending the event loop  ------------
void 
BigNtuplizer::endJob() 
{
}

// ------------ method called when starting to processes a run  ------------
/*
  void 
  BigNtuplizer::beginRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a run  ------------
/*
  void 
  BigNtuplizer::endRun(edm::Run const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
  void 
  BigNtuplizer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
  void 
  BigNtuplizer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
  {
  }
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
BigNtuplizer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(BigNtuplizer);
