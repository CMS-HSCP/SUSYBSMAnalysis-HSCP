#ifndef SUSYBSMAnalysis_Analyzer_Analyzer_h
#define SUSYBSMAnalysis_Analyzer_Analyzer_h

// -*- C++ -*-
// 
// 
//  Class : Analyzer
// 
// 
//  Inspired by class Analyzer SUSYBSMAnalysis/Analyzer/plugins/Analyzer.cc
// 
//  Author : Raphael Haeberle
//        Created : Wed, 23 February 2022 16:12:48 GMT+1
// 
// 

// ~~~~~~~~~c++ include files ~~~~~~~~~
#include <memory>
#include <vector>
#include <string>
#include <map>
#include <exception>
#include <unordered_map>
#include <iostream>
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

// ~~~~~~~~~ CMSSW include files ~~~~~~~~~
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
// Muons CSC segments
#include "DataFormats/CSCRecHit/interface/CSCSegment.h"

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
#include "DataFormats/JetReco/interface/PFJet.h"
#include "DataFormats/JetReco/interface/PFJetCollection.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "HLTrigger/HLTcore/interface/HLTPrescaleProvider.h"
#include "HLTrigger/HLTcore/interface/FractionalPrescale.h"

// ~~~~~~~~~ user include files ~~~~~~~~~

#define FWCORE
#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "SUSYBSMAnalysis/Analyzer/interface/DeDxUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TOFUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TupleMaker.h"
#include "SUSYBSMAnalysis/Analyzer/interface/SaturationCorrection.h"
#include "SUSYBSMAnalysis/Analyzer/interface/MCWeight.h"

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "DataFormats/L1Trigger/interface/BXVector.h"
#include "DataFormats/L1Trigger/interface/L1TObjComparison.h"
#include "DataFormats/L1TGlobal/interface/GlobalAlgBlk.h"
#include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"



#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"
#include "DataFormats/Provenance/interface/ParameterSetID.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Utilities/interface/InputTag.h"


using namespace std;

class TupleMaker;
class MCWeight;

class Analyzer : public edm::EDAnalyzer {

public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  void AssoGenID(const vector<bool>& l1decision);

  void initializeCuts(edm::Service<TFileService>& fs,
                      vector<double>& CutPt,
                      vector<double>& CutI,
                      vector<double>& CutTOF,
                      vector<double>& CutPt_Flip,
                      vector<double>& CutI_Flip,
                      vector<double>& CutTOF_Flip);


  bool passPreselection(const susybsm::HSCParticle& hscp,
                        const reco::DeDxHitInfo* dedxHits,
                        const reco::DeDxData* dedxSObj,
                        const reco::DeDxData* dedxMObj,
                        const reco::MuonTimeExtra* tof,
                        const edm::Event& iEvent,
                        float Event_Weight,
                        Tuple* tuple,
                        const double& GenBeta,
                        bool RescaleP,
                        const double& RescaleI,
                        const double& RescaleT,
                        double MassErr,
                        bool Ih_Iso_cut = true);

  bool passSelection(const susybsm::HSCParticle& hscp,
                     const reco::DeDxData* dedxSObj,
                     const reco::DeDxData* dedxMObj,
                     const reco::MuonTimeExtra* tof,
                     const edm::Event& iEvent,
                     float Event_Weight,
                     const int& CutIndex,
                     Tuple*& tuple,
                     const bool isFlip,
                     const double& GenBeta,
                     bool RescaleP,
                     const double& RescaleI,
                     const double& RescaleT);

  bool passTrigger(const edm::Event& iEvent, bool isData, bool isCosmic = false, L1BugEmulator* emul = nullptr);
  //double SegSep(const susybsm::HSCParticle& hscp, const edm::Event& iEvent, double& minPhi, double& minEta);
  static constexpr unsigned int maxPhysicsTriggers = 512;

private:
  virtual void beginJob() override;
  virtual void beginRun(const edm::Run& run,const edm::EventSetup& iSetup) override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
  //virtual void isPixelTrack(const edm::Ref<std::vector<Trajectory>>&, bool&, bool&);

  // ----------member data ---------------------------
  edm::EDGetTokenT<vector<susybsm::HSCParticle>> hscpToken_;
  edm::EDGetTokenT<edm::ValueMap<susybsm::HSCPIsolation>> hscpIsoToken_;
  edm::EDGetTokenT<BXVector<GlobalAlgBlk>> l1resToken_;
  edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxToken_;
  edm::EDGetTokenT<vector<reco::Vertex>> offlinePrimaryVerticesToken_;
  
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> pfMETToken_;
  edm::EDGetTokenT<reco::PFJetCollection> pfJetToken_;
  edm::EDGetTokenT<std::vector<reco::CaloMET>> CaloMETToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  
  vector<string> trigger_met_, trigger_mu_;
  vector<double> CutPt_, CutI_, CutTOF_;
  vector<double> CutPt_Flip_, CutI_Flip_, CutTOF_Flip_;

  map<string, TProfile*> HCuts;
  const int nbl1names = 12;
  //following is hard coded for my actual purpose, a given list of L1 seeds of interest
  const string ListL1Names[12] = {"L1_SingleMu22","L1_SingleMu25","L1_ETMHF90_HTT60er","L1_ETMHF100","L1_ETMHF100_HTT60er","L1_ETMHF110","L1_ETMHF110_HTT60er","L1_ETMHF120","L1_ETMHF120_HTT60er","L1_ETT2000","L1_ETM150","L1_ETMHF150"};
  vector<double> EffL1Seeds;
  vector<int> L1Num;
  vector<int> L1Denom;
  vector<bool> L1Dec; 

  bool* HSCPTk;

  const reco::MuonTimeExtra* tof;
  const reco::MuonTimeExtra* dttof;
  const reco::MuonTimeExtra* csctof;


  double OpenAngle = -1;  //global variable needed by PassPreselection... Ugly isn't it?!
  double TreeDXY = -1;
  double TreeDZ = -1;
  bool isCosmicSB = false;
  bool isSemiCosmicSB = false;

  int TypeMode_;
  int SampleType_;
  string SampleName_;
  string Period_;
  int nbpasspresel = 0;
  int nbtot = 0;


  bool SkipSelectionPlot_;

  // binning for the pT, mass, IP distributions

  //Variables used in the TOF only HSCP search
  float DTRegion = 0.9;       //Define the dividing line between DT and
  float CSCRegion = 0.9;      //CSC regions of CMS
  float CosmicMinDz = 70.;    //Min dz displacement to be tagged as cosmic muon
  float CosmicMaxDz = 120.;   //Max dz displacement for cosmic tagged tracks
  double minSegEtaSep = 0.1;  //Minimum eta separation between SA track and muon segment on opposite side of detector
  int minMuStations = 2;
  
  // Thresholds for candidate preselection
  double GlobalMaxEta = 2.1;      // cut on inner tracker track eta
  double GlobalMaxV3D = 99999;    //0.50 cuts away signal;   // cut on 3D distance (cm) to closest vertex
  double GlobalMaxDZ = 0.50;      // cut on 1D distance (cm) to closest vertex in "Z" direction
  double GlobalMaxDXY = 0.50;     // cut on 2D distance (cm) to closest vertex in "R" direction
  double GlobalMaxChi2 = 5.0;     // cut on Track maximal Chi2/NDF
  int GlobalMinQual = 2;          // cut on track quality (2 meaning HighPurity tracks)
  unsigned int GlobalMinNOH = 8;  //7AMSB;      // cut on number of (valid) track pixel+strip hits
  int GlobalMinNOPH = 2;          // cut on number of (valid) track pixel hits
  double GlobalMinFOVH = 0.8;     //0.0AMSB;    // cut on fraction of valid track hits
  unsigned int GlobalMaxNOMHTillLast =
      99999;  //1AMSB;     // cut on the number of missing hits from IP till last hit (excluding hits behind the last hit)
  double GlobalMinFOVHTillLast =
      -99999;  //0.85AMSB;   // cut on the fraction of valid hits divided by total expected hits until the last one
  unsigned int GlobalMinNOM =
      6;  //7AMSB;      // cut on number of dEdx hits (generally equal to #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used)
  double GlobalMinNDOF = 8;            // cut on number of     DegreeOfFreedom used for muon TOF measurement
  double GlobalMinNDOFDT = 6;          // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
  double GlobalMinNDOFCSC = 6;         // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
  double GlobalMaxTOFErr = 0.15;       //0.07;   // cut on error on muon TOF measurement
  double GlobalMaxPterr = 0.25;        //0.50;//0.25;   // cut on error on track pT measurement
  double GlobalMaxTIsol = 50;          // cut on tracker isolation (SumPt)
  double GlobalMaxRelTIsol = 9999999;  // cut on relative tracker isolation (SumPt/Pt)
  double GlobalMaxEIsol = 0.30;        // cut on calorimeter isolation (E/P)
  double GlobalMinPt = 55.00;          // cut on pT    at PRE-SELECTION
  double GlobalMinIs = 0.0;            // cut on dEdxS at PRE-SELECTION (dEdxS is generally a  discriminator)
  double GlobalMinIm = 0.0;            // cut on dEdxM at PRE-SELECTION (dEdxM is generally an estimator    )
  double GlobalMinTOF = 1.0;           // cut on TOF   at PRE-SELECTION

  bool skipPixel = true;
  bool useTemplateLayer = false;

  const int MaxPredBins =
      6;  // The maximum number of different bins prediction is done in for any of the analyses (defines array size)

  double dEdxK_Data = 2.580;
  double dEdxC_Data = 3.922;
  double dEdxK_MC = 2.935;
  double dEdxC_MC = 3.197;

  //=============================================================
  Tuple* tuple;
  TupleMaker* tuple_maker;
  //=============================================================

  TH3F* dEdxTemplates = nullptr;
  double DeDxSF_0 = 1.00000;  // [0]  unchanged
  double DeDxSF_1 = 1.41822;  // [1]  Pixel data to SiStrip data
  double dEdxSF[2] = {DeDxSF_0, DeDxSF_1};
  double DeDxK = 0.0;
  double DeDxC = 0.0;

  dedxGainCorrector trackerCorrector;
  string DeDxTemplate;  // "MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", "Data13TeV16_dEdxTemplate.root"
  bool enableDeDxCalibration;
  string DeDxCalibration;  //"Data13TeVGains_v2.root" if Data
  string Geometry;         //CMS_GeomTree.root
  string TimeOffset;       //MuonTimeOffset.txt
  muonTimingCalculator tofCalculator;

  double FMIPX = 4;

  unsigned int STree = 0;
  unsigned int SGTree = 0;
  bool useClusterCleaning;
  bool isData;
  bool isBckg;
  bool isSignal;

  unsigned int CurrentRun_ = 0;

  MCWeight* mcWeight;

  float EventWeight_ = 1.;
  //double SampleWeight_ = 1.;
  double CrossSection_ = 1.;
  vector<float> PUSystFactor_;
  
  unsigned int TrigInfo_ = 0;  //1 -mu only, 2- met only, 3 mu and met
  //
  TRandom3* RNG = nullptr;
  bool is2016;
  bool is2016G;
  bool isMCglobal = false;
  //
  double preTrackingChangeL1IntLumi_ = 29679.982;  // pb
  double IntegratedLuminosity_ = 33676.4;          //13TeV16
  //
  const std::string pixelCPE_;
  const double trackProbQCut_;
  const int debugLevel_;
  //
  const edm::ESGetToken<TrackerTopology, TrackerTopologyRcd> tTopoToken_;
  const edm::ESGetToken<TrackerGeometry,TrackerDigiGeometryRecord> tkGeomToken_;
  const edm::ESGetToken<PixelClusterParameterEstimator, TkPixelCPERecord> pixelCPEToken_;
  HLTPrescaleProvider hltPSProv_;
  std::string hltProcess_; 
  string scenarios_;
};
#endif
