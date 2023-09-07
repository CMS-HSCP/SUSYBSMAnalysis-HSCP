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
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "FWCore/Common/interface/TriggerNames.h"

#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/PatCandidates/interface/TriggerObjectStandAlone.h"
#include "DataFormats/PatCandidates/interface/PackedTriggerPrescales.h"
#include "DataFormats/PatCandidates/interface/Electron.h"

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
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h"
#include "L1Trigger/L1TGlobal/interface/L1TGlobalUtil.h"
//#include "DataFormats/L1TGlobalTrigger/interface/GlobalAlgBlk.h"
//#include "DataFormats/L1TGlobalTrigger/interface/L1TGlobalOutput.h"


//
// ~~~~~~~~~ user include files ~~~~~~~~~
#define FWCORE

#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "RecoLocalTracker/Records/interface/TkPixelCPERecord.h"
#include "RecoLocalTracker/ClusterParameterEstimator/interface/PixelClusterParameterEstimator.h"
#include "Geometry/Records/interface/TrackerTopologyRcd.h"
#include "DataFormats/TrackerCommon/interface/TrackerTopology.h"
#include "FWCore/ParameterSet/interface/ParameterSetDescription.h"
#include "FWCore/ParameterSet/interface/ConfigurationDescriptions.h"
#include "RecoEgamma/EgammaTools/interface/ConversionTools.h"

#include "DataFormats/ParticleFlowReco/interface/PFDisplacedVertex.h"

#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "SUSYBSMAnalysis/Analyzer/interface/DeDxUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TOFUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TupleMaker.h"
#include "SUSYBSMAnalysis/Analyzer/interface/SaturationCorrection.h"
#include "SUSYBSMAnalysis/Analyzer/interface/MCWeight.h"
#include "SUSYBSMAnalysis/Analyzer/interface/Regions.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TrigToolsFuncs.h"


using namespace std;
class TupleMaker;
class MCWeight;

class Analyzer : public edm::EDAnalyzer {
public:
  explicit Analyzer(const edm::ParameterSet&);
  ~Analyzer();

  static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

  float muonRecoSFsForTrackEta(float eta, int syst);
  float muonIdSFsForTrackEta(float eta, int syst);
  float muonTriggerSFsForTrackEta(float eta, int syst);
  float triggerSystFactor(float eta, float beta, int syst);

  void initializeCuts(edm::Service<TFileService>& fs,
                      vector<float>& CutPt,
                      vector<float>& CutI,
                      vector<float>& CutTOF,
                      vector<float>& CutPt_Flip,
                      vector<float>& CutI_Flip,
                      vector<float>& CutTOF_Flip);
  template <typename T, size_t n>
  bool passPreselection(T (&passedCutsArray)[n], bool verbose);

  bool passSelection(const reco::TrackRef track,
                     const reco::DeDxData* dedxSObj,
                     const reco::DeDxData* dedxMObj,
                     const reco::MuonTimeExtra* tof,
                     const int& CutIndex,
                     Tuple* tuple,
                     const bool isFlip,
                     const float GenBeta,
                     const bool RescaleP,
                     const float RescaleI,
                     const float RescaleT);

  float shiftForPt(const float& pt, const float& eta, const float& phi, const int& charge);
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
                     Tuple* tuple,
                     const float GenBeta,
                     float MassErr,
                     const float closestBackgroundPDGsIDs[]);
   const reco::Candidate* findFirstMotherWithDifferentID(const reco::Candidate *particle);
   const reco::Candidate* findOriginalMotherWithSameID(const reco::Candidate *particle);


private:
  virtual void beginJob() override;
  virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
  virtual void endJob() override;
//  virtual void isPixelTrack(const edm::Ref<std::vector<Trajectory>>&, bool&, bool&);


  // ----------member data ---------------------------
  edm::EDGetTokenT<vector<susybsm::HSCParticle>> hscpToken_;
  edm::EDGetTokenT<reco::TrackCollection> genTrackToken_;
  edm::EDGetTokenT<edm::ValueMap<susybsm::HSCPIsolation>> hscpIsoToken_;
  edm::EDGetTokenT<susybsm::MuonSegmentCollection> muonSegmentToken_;
  edm::EDGetTokenT<reco::DeDxHitInfoAss> dedxToken_;
  edm::EDGetTokenT<edm::ValueMap<int>> dedxPrescaleToken_;
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
  edm::EDGetTokenT<vector<reco::Conversion> > conversionsToken_;
  edm::EDGetTokenT<reco::GsfElectronCollection> electronToken_;

  edm::EDGetTokenT<edm::ValueMap<bool> > electron_cutbasedID_decisions_veto_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_cutbasedID_decisions_loose_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_cutbasedID_decisions_medium_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_cutbasedID_decisions_tight_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wp80_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wp90_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wpHZZ_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wpLoose_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_mvaNoIsoID_decisions_wp80_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_mvaNoIsoID_decisions_wp90_Token_;
  edm::EDGetTokenT<edm::ValueMap<bool> > electron_mvaNoIsoID_decisions_wpLoose_Token_;
  edm::EDGetTokenT<edm::TriggerResults> triggerResultsToken_;
  edm::EDGetTokenT<pat::PackedTriggerPrescales> triggerPrescalesToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> trigEventToken_ ;
  string filterName_;
  string pathName_;
  string triggerPathNamesFile_;
  string muonHLTFilterNamesFile_;
  static const int NTriggersMAX = 1201;
  string triggerPathNames[NTriggersMAX];
  static const int MAX_MuonHLTFilters = 69;
  string muonHLTFilterNames[MAX_MuonHLTFilters];


  bool matchToHLTTrigger_;
  edm::EDGetTokenT<std::vector<reco::PFMET>> pfMETToken_;
  edm::EDGetTokenT<reco::PFJetCollection> pfJetToken_;
  edm::EDGetTokenT<std::vector<reco::CaloMET>> caloMETToken_;
  edm::EDGetTokenT<std::vector<reco::CaloJet>> caloJetToken_;
  edm::EDGetTokenT<trigger::TriggerEvent> triggerSummaryToken_;
  edm::EDGetTokenT<std::vector<PileupSummaryInfo>> pileupInfoToken_;
  edm::EDGetTokenT<std::vector<reco::GenParticle>> genParticleToken_;
  edm::EDGetTokenT<edm::Association<reco::GenParticleCollection>> trackToGenToken_;
  edm::EDGetTokenT<reco::PFCandidateCollection> pfCandToken_;
  edm::EDGetTokenT<GenEventInfoProduct> genEventToken_; // for reading generator weight
  edm::EDGetTokenT<L1GlobalTriggerReadoutRecord> l1GtReadoutRecordToken_;


  edm::Handle<edm::ValueMap<bool> > electron_cutbasedID_decisions_veto;
  edm::Handle<edm::ValueMap<bool> > electron_cutbasedID_decisions_loose;
  edm::Handle<edm::ValueMap<bool> > electron_cutbasedID_decisions_medium;
  edm::Handle<edm::ValueMap<bool> > electron_cutbasedID_decisions_tight;
  edm::Handle<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wp80;
  edm::Handle<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wp90;
  edm::Handle<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wpHZZ;
  edm::Handle<edm::ValueMap<bool> > electron_mvaIsoID_decisions_wpLoose;
  edm::Handle<edm::ValueMap<bool> > electron_mvaNoIsoID_decisions_wp80;
  edm::Handle<edm::ValueMap<bool> > electron_mvaNoIsoID_decisions_wp90;
  edm::Handle<edm::ValueMap<bool> > electron_mvaNoIsoID_decisions_wpLoose;

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

  bool isCosmicSB = false;
  bool isSemiCosmicSB = false;

  unsigned int typeMode_;
  unsigned int sampleType_;
  string sampleName_;
  string period_;

  bool tapeRecallOnly_, doBefTrigPlots_, doBefPreSplots_, doPostPreSplots_, doSystsPlots_;

  // binning for the pT, mass, IP distributions
  float ptHistoUpperBound_ = 4000;
  float pHistoUpperBound_ = 10000;
  float massHistoUpperBound_ = 4000;
  unsigned int massNBins_ = 400;
  float cutOnIPbound_ = 1.0;
  unsigned int predBins_ = 0;
  unsigned int etaBins_ = 60;

  // Ias quantiles and pT_cut used to validate the background estimate method in data
  //float Ias_quantiles[5]={ 0.039, 0.045, 0.053, 0.064, 0.082 }; //data or signal
  //Ias-quantiles update
  //Ias-quantiles { 40%, 50%, 60%, 70%, 80%, 90%, 99%, 99.9% }
  float Ias_quantiles[8]={ 0.014565036, 0.017987774, 0.022399569, 0.028518069, 0.038047370, 0.056746799, 0.13331622, 0.22018057 }; //data or signal -- IAS STRIP ONLY NO FSTRIP CUT
  float Fpix_quantiles[12]={ 0., 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7,0.8,0.9,0.99,1.0 }; //data or signal -- F PIXEL ONLY
 
  //float Ias_quantiles[5]={ 0.037, 0.042, 0.048, 0.056, 0.066 }; //data or signal //WIP new quantiles determined with new preselection cuts
  //pT cut update 60-->70 GeV
  float pT_cut = 70;

  // binning for eta, ih, p, mass distributions used to validate the background estimate method in data
  //eta 120-->200
  int reg_etabins_ = 200;
  //ih 200-->2000
  int reg_ihbins_ = 2000;
  //p 200-->2000
  int reg_pbins_ = 2000;
  //mass 50-->200
  int reg_massbins_ = 200;


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
  // TAV: I think these could / should be const 
  float globalMaxEta_, globalMinPt_, globalMaxPt_;
  unsigned int globalMinNOPH_;
  float globalMinFOVH_;
  unsigned int globalMinNOM_;
  float globalMaxChi2_, globalMaxEoP_, globalMaxDZ_, globalMaxDXY_, globalMaxTIsol_, globalMinDeltaRminJet_, globalMaxMiniRelIsoAll_, globalMinIh_, globalMinTrackProbQCut_, globalMaxTrackProbQCut_, globalMinTrackProbXYCut_;
  float globalMaxTrackProbXYCut_;
  unsigned int minMuStations_;
  float globalMinIs_, globalMinTOF_;
  bool puTreatment_, createGiTemplates_, createAndExitGitemplates_;
  int NbPuBins_;
  vector<int> PuBins_;
 
  double GiSysParamOne_; 
  double GiSysParamTwo_; 
  vector<int> NominalEntries_;

  float GlobalMinNDOF = 8;            // cut on number of     DegreeOfFreedom used for muon TOF measurement
  float GlobalMinNDOFDT = 6;          // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
  float GlobalMinNDOFCSC = 6;         // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
  float GlobalMaxTOFErr = 0.15;       //0.07;   // cut on error on muon TOF measurement
  bool exitWhenGenMatchNotFound_;
  bool useTemplateLayer_ = false;

  // The maximum number of different bins prediction is done in for any of the analyses (defines array size)
  const int MaxPredBins = 6;

  //=============================================================
  Tuple* tuple;
  TupleMaker* tuple_maker;
  //=============================================================

  TH3F* dEdxTemplates = nullptr;
  vector<TH3F*> dEdxTemplatesPU;

  float dEdxSF_0_, dEdxSF_1_;
  float dEdxSF[2] = {dEdxSF_0_, dEdxSF_1_};
  float dEdxK_;
  float dEdxC_;
  float globalIas_;
  float globalFiStrips_;
  float globalIh_;

  dedxGainCorrector trackerCorrector;
  string dEdxTemplate_;
  bool enableDeDxCalibration_;
  string timeOffset_;
  muonTimingCalculator tofCalculator;

  unsigned int saveTree_;

  bool useClusterCleaning, isData, isBckg, isSignal;

  unsigned int currentRun_ = 0;

  MCWeight* mcWeight;

  float eventWeight_ = 1.;
  float GeneratorWeight_ = 1.;
  float GeneratorBinningValues_ = 1.;
  //double SampleWeight_ = 1.;
  float CrossSection_ = 1.;

  vector<float> PUSystFactor_;

  TRandom3* RNG = nullptr;
  TRandom3* RNG2 = nullptr;
  //TRandom3* RNG3 = nullptr;
  bool is2016;
  bool is2016G;

  bool isMCglobal = false;

  float IntegratedLuminosity_ = 33676.4;          //13TeV16

  const std::string pixelCPE_;
  const int debug_;
  const bool hasMCMatch_,calcSyst_;

  unsigned int trigInfo_;

  static constexpr const char* const MOD = "Analyzer";
  int totMu22;
  int totMu22or25;
  int totLastMu;

  int passMu22;
  int passMu25;
  int passLastMu;

  int passMu22PostS;
  int passMu25PostS;
  int passLastMuPostS;

  TEfficiency* effl1Mu22;
  TEfficiency* effl1Mu22or25;
  TEfficiency* effl1LastMu;
  TEfficiency* effHltMu50;

  TEfficiency* effl1Mu22PostS;
  TEfficiency* effl1Mu22or25PostS;
  TEfficiency* effl1LastMuPostS;
  TEfficiency* effHltMu50PostS;
};
#endif
