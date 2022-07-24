#ifndef SUSYBSMAnalysis_Analyzer_Analyzer_h
#define SUSYBSMAnalysis_Analyzer_Analyzer_h
// -*- C++ -*-
//
// Package:    SUSYBSMAnalysis/Analyzer
// Class:      Analyzer
//
/**\class Analyzer Analyzer.cc SUSYBSMAnalysis/Analyzer/plugins/Analyzer.cc
*/
//
// Original Author:  Emery Nibigira
//         Created:  Thu, 01 Apr 2021 07:04:53 GMT
// Modifications by Tamas Almos Vami
//

// ~~~~~~~~~c++ include files ~~~~~~~~~
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <exception>
#include <unordered_map>

// ~~~~~~~~~ ROOT include files ~~~~~~~~~
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TObject.h"
#include "TVector3.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TTree.h"
#include "TProfile.h"
#include "TLorentzVector.h"
//#include "TCanvas.h"

// ~~~~~~~~~ CMSSW include files ~~~~~~~~~
//#include "FWCore/Framework/interface/one/EDAnalyzer.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"
#include "DataFormats/Common/interface/View.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "DataFormats/MuonReco/interface/MuonSelectors.h"
// Muons CSC segments
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"
#include "DataFormats/PatCandidates/interface/PackedCandidate.h"
#include "DataFormats/ParticleFlowCandidate/interface/PFCandidate.h"
#include "DataFormats/PatCandidates/interface/IsolatedTrack.h"
#include "DataFormats/PatCandidates/interface/PFIsolation.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHit.h"
#include "TrackingTools/TrackAssociator/interface/TrackDetectorAssociator.h"
#include "TrackingTools/PatternTools/interface/Trajectory.h"
#include "TrackingTools/TrackFitters/interface/TrajectoryStateCombiner.h"
#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/PatternTools/interface/TrajTrackAssociation.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHit.h"

#include "DataFormats/METReco/interface/CaloMET.h"
#include "DataFormats/METReco/interface/PFMET.h"
#include "DataFormats/METReco/interface/PFMETFwd.h"
#include "DataFormats/JetReco/interface/CaloJet.h"
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"

#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"

#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

// ~~~~~~~~~ user include files ~~~~~~~~~
#define FWCORE
#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "SUSYBSMAnalysis/Analyzer/interface/DeDxUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TOFUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TupleMaker.h"
#include "SUSYBSMAnalysis/Analyzer/interface/SaturationCorrection.h"
#include "SUSYBSMAnalysis/Analyzer/interface/MCWeight.h"
#include "SUSYBSMAnalysis/Analyzer/interface/Regions.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"

#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"

using namespace std;

class TupleMaker;
class MCWeight;

class Analyzer : public edm::EDAnalyzer {
public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  float scaleFactor(float eta);

  void initializeCuts(edm::Service<TFileService>& fs,
                      vector<float>& CutPt,
                      vector<float>& CutI,
                      vector<float>& CutTOF,
                      vector<float>& CutPt_Flip,
                      vector<float>& CutI_Flip,
                      vector<float>& CutTOF_Flip);

  bool passPreselection(const reco::TrackRef track,
                        const reco::DeDxHitInfo* dedxHits,
                        const reco::DeDxData* dedxSObj,
                        const reco::DeDxData* dedxMObj,
                        const reco::MuonTimeExtra* tof,
                        const edm::Event& iEvent,
                        const edm::EventSetup& iSetup,
                        const float pixelProbs[],
                        const float Event_Weight,
                        Tuple* tuple,
                        const float GenBeta,
                        const bool RescaleP,
                        const float RescaleI,
                        const float RescaleT,
                        float MassErr,
                        const bool Ih_Iso_cut,
                        const float closestBackgroundPDGsIDs[]);

  bool passSelection(const reco::TrackRef track,
                     const reco::DeDxData* dedxSObj,
                     const reco::DeDxData* dedxMObj,
                     const reco::MuonTimeExtra* tof,
                     const edm::Event& iEvent,
                     const float Event_Weight,
                     const int& CutIndex,
                     Tuple* tuple,
                     const bool isFlip,
                     const float GenBeta,
                     const bool RescaleP,
                     const float RescaleI,
                     const float RescaleT);

  float RescaledPt(const float& pt, const float& eta, const float& phi, const int& charge);
  GlobalPoint getOuterHitPos(const edm::EventSetup& iSetup, const reco::DeDxHitInfo* dedxHits);
  float SegSep(const reco::TrackRef track, const edm::Event& iEvent, float& minPhi, float& minEta);
  float combineProbs(float probOnTrackWMulti, int numRecHits) const;
  bool isHSCPgenID(const reco::GenParticle& gen);
  void calculateSyst(const reco::TrackRef track,
                     const reco::DeDxHitInfo* dedxHits,
                     const reco::DeDxData* dedxSObj,
                     const reco::DeDxData* dedxMObj,
                     const reco::MuonTimeExtra* tof,
                     const edm::Event& iEvent,
                     const edm::EventSetup& iSetup,
                     const float pixelProbs[],
                     const float Event_Weight,
                     Tuple* tuple,
                     const float GenBeta,
                     float MassErr,
                     const bool Ih_Iso_cut,
                     const float closestBackgroundPDGsIDs[]);

private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
//  virtual void isPixelTrack(const edm::Ref<std::vector<Trajectory>>&, bool&, bool&);


  // ----------member data ---------------------------
  edm::EDGetTokenT<vector<susybsm::HSCParticle>> hscpToken_;
  edm::EDGetTokenT<edm::ValueMap<susybsm::HSCPIsolation>> hscpIsoToken_;
  edm::EDGetTokenT<susybsm::MuonSegmentCollection> muonSegmentToken_;
  edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> muonTimeToken_;  // for reading inverse beta
  edm::EDGetTokenT<reco::MuonTimeExtraMap> muonDtTimeToken_;
  edm::EDGetTokenT<reco::MuonTimeExtraMap> muonCscTimeToken_;
  edm::EDGetTokenT<DTRecSegment4DCollection> muonDtSegmentToken_;
  edm::EDGetTokenT<CSCSegmentCollection> muonCscSegmentToken_;
  edm::EDGetTokenT<vector<reco::Vertex>> offlinePrimaryVerticesToken_;
  edm::EDGetTokenT<vector<reco::Vertex>> inclusiveSecondaryVerticesToken_;
  edm::EDGetTokenT<LumiScalersCollection> lumiScalersToken_;
  edm::EDGetTokenT<vector<reco::Track>> refittedStandAloneMuonsToken_;
  edm::EDGetTokenT<reco::BeamSpot> offlineBeamSpotToken_;
  edm::EDGetTokenT<vector<reco::Muon>> muonToken_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> pfMETToken_;
  edm::EDGetTokenT<reco::PFJetCollection> pfJetToken_;
  edm::EDGetTokenT<std::vector<reco::CaloMET>> caloMETToken_;
  edm::EDGetTokenT<std::vector<reco::CaloJet>> caloJetToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfoToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> trackToGenToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventToken_; // for reading generator weight

  vector<string> trigger_met_, trigger_mu_;

  vector<float> CutPt_, CutI_, CutTOF_;
  vector<float> CutPt_Flip_, CutI_Flip_, CutTOF_Flip_;
  //map<string, vector<float>> VCuts;

  map<string, TProfile*> HCuts;

  bool* HSCPTk;
  bool* HSCPTk_SystP;
  bool* HSCPTk_SystI;
  bool* HSCPTk_SystT;
  bool* HSCPTk_SystM;
  bool* HSCPTk_SystPU;
  bool* HSCPTk_SystHUp;
  bool* HSCPTk_SystHDown;
  float* MaxMass;
  float* MaxMass_SystP;
  float* MaxMass_SystI;
  float* MaxMass_SystT;
  float* MaxMass_SystM;
  float* MaxMass_SystPU;
  float* MaxMass_SystHUp;
  float* MaxMass_SystHDown;

  const reco::MuonTimeExtra* tof;
  const reco::MuonTimeExtra* dttof;
  const reco::MuonTimeExtra* csctof;

  float OpenAngle = -1;  //global variable needed by PassPreselection... Ugly isn't it?!
  float TreeDXY = -1;
  float TreeDZ = -1;
  float TreeprobQonTrack = -1;
  float TreeprobQonTracknoL1 = -1;
  float TreeprobXYonTrack = -1;
  float TreeprobXYonTracknoL1 = -1;
  bool isCosmicSB = false;
  bool isSemiCosmicSB = false;

  unsigned int typeMode_;
  unsigned int sampleType_;
  string sampleName_;
  string period_;

  bool skipSelectionPlot_;

  // binning for the pT, mass, IP distributions
  float ptHistoUpperBound_ = 4000;
  float massHistoUpperBound_ = 4000;
  unsigned int massNBins_ = 400;
  float cutOnIPbound_ = 1.0;
  unsigned int predBins_ = 0;
  unsigned int etaBins_ = 60;
  
  // Ias quantiles and pT_cut used to validate the background estimate method in data
  float Ias_quantiles[5]={ 0.039, 0.045, 0.053, 0.064, 0.082 }; //data or signal
  float pT_cut = 60;
  
  // binning for eta, ih, p, mass distributions used to validate the background estimate method in data
  int reg_etabins_ = 120;
  int reg_ihbins_ = 200;
  int reg_pbins_ = 200;
  int reg_massbins_ = 50;
  
  
  float dEdxS_UpLim_ = 1.0;
  float dEdxM_UpLim_ = 30.0;
  unsigned int numDzRegions_ = 6;



  //Variables used in the TOF only HSCP search
  float DTRegion = 0.9;       //Define the dividing line between DT and
  float CSCRegion = 0.9;      //CSC regions of CMS
  float CosmicMinDz = 70.;    //Min dz displacement to be tagged as cosmic muon
  float CosmicMaxDz = 120.;   //Max dz displacement for cosmic tagged tracks
  float minSegEtaSep = 0.1;  //Minimum eta separation between SA track and muon segment on opposite side of detector

  // Thresholds for candidate preselection
  float globalMaxEta_, globalMinPt_;
  unsigned int globalMinNOPH_;
  float globalMinFOVH_;
  unsigned int globalMinNOM_;
  float globalMaxChi2_, globalMaxEoP_, globalMaxDZ_, globalMaxDXY_, globalMaxTIsol_, globalMiniRelIsoAll_, globalMinIh_, globalMinTrackProbQCut_, globalMaxTrackProbQCut_, globalMinTrackProbXYCut_;
  unsigned int minMuStations_;
  float globalMinIs_, globalMinTOF_;
  float GlobalMinNDOF = 8;            // cut on number of     DegreeOfFreedom used for muon TOF measurement
  float GlobalMinNDOFDT = 6;          // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
  float GlobalMinNDOFCSC = 6;         // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
  float GlobalMaxTOFErr = 0.15;       //0.07;   // cut on error on muon TOF measurement

  bool skipPixel_ = true;
  bool useTemplateLayer_ = false;

  // The maximum number of different bins prediction is done in for any of the analyses (defines array size)
  const int MaxPredBins = 6; 

  //=============================================================
  Tuple* tuple;
  TupleMaker* tuple_maker;
  //=============================================================

  TH3F* dEdxTemplates = nullptr;

  float dEdxSF_0_, dEdxSF_1_;
  float dEdxSF[2] = {dEdxSF_0_, dEdxSF_1_};
  float dEdxK_;
  float dEdxC_;

  dedxGainCorrector trackerCorrector;
  string dEdxTemplate_;  // "MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", "Data13TeV16_dEdxTemplate.root"
  bool enableDeDxCalibration_;
  string dEdxCalibration_, geometry_, timeOffset_;
  muonTimingCalculator tofCalculator;

  float theFMIPX_ = 4;

  unsigned int saveTree_, saveGenTree_;

  bool useClusterCleaning, isData, isBckg, isSignal;

  unsigned int CurrentRun_ = 0;

  MCWeight* mcWeight;

  float EventWeight_ = 1.;
  float GeneratorWeight_ = 1.;
  float GeneratorBinningValues_ = 1.;
  //double SampleWeight_ = 1.;
  float CrossSection_ = 1.;

  vector<float> PUSystFactor_;

  unsigned int TrigInfo_ = 0;  //1 -mu only, 2- met only, 3 mu and met

  TRandom3* RNG = nullptr;
  bool is2016;
  bool is2016G;

  bool isMCglobal = false;

  float IntegratedLuminosity_ = 33676.4;          //13TeV16

  const std::string pixelCPE_;
  const int debug_;
  const bool hasMCMatch_, doTriggering_,calcSyst_;

  static constexpr const char* const MOD = "Analyzer";

};
#endif
