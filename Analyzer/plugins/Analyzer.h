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
//
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
//#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPDeDxInfo.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
// Muons CSC segments
//#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

#include "DataFormats/VertexReco/interface/Vertex.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/Common/interface/ValueMap.h"

#include "DataFormats/Common/interface/TriggerResults.h"
#include "FWCore/Common/interface/TriggerNames.h"

// ~~~~~~~~~ user include files ~~~~~~~~~
#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "SUSYBSMAnalysis/Analyzer/interface/Calibration.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TOFUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/TupleSaver.h"

using namespace std;

class TupleSaver;

class Analyzer : public edm::EDAnalyzer {
   public:
      explicit Analyzer(const edm::ParameterSet&);
      ~Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      double scaleFactor(double eta);

      void initializeCuts(edm::Service<TFileService> &fs);

      bool PassTrigger(const edm::Event& iEvent, bool isData, bool isCosmic=false, L1BugEmulator* emul=nullptr);

      //bool PassPreselection(const susybsm::HSCParticle& hscp, const reco::DeDxHitInfo* dedxHits,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const edm::Event&, stPlots* st, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT, double MassErr);

      bool passPreselection(
         const susybsm::HSCParticle& hscp, 
         const reco::DeDxHitInfo* dedxHits,  
         const reco::DeDxData* dedxSObj, 
         const reco::DeDxData* dedxMObj,
         const reco::MuonTimeExtra* tof,
         const edm::Event& iEvent, 
         Tuple* &tuple, 
         const double& GenBeta, 
         bool RescaleP, 
         const double& RescaleI, 
         const double& RescaleT, 
         double MassErr);

      bool passSelection(
         const susybsm::HSCParticle& hscp,  
         const reco::DeDxData* dedxSObj, 
         const reco::DeDxData* dedxMObj, 
         const reco::MuonTimeExtra* tof, 
         const edm::Event& iEvent,
         const int& CutIndex, 
         Tuple* &tuple, 
         const bool isFlip, 
         const double& GenBeta, 
         bool RescaleP, 
         const double& RescaleI, 
         const double& RescaleT);

      bool passTrigger(const edm::Event& iEvent, bool isData, bool isCosmic=false, L1BugEmulator* emul=nullptr);

      //int  muonStations(const reco::HitPattern& hitPattern);
      double RescaledPt(const double& pt, const double& eta, const double& phi, const int& charge);
      TVector3 getOuterHitPos(const reco::DeDxHitInfo* dedxHits);
      double SegSep(const susybsm::HSCParticle& hscp, const edm::Event& iEvent, double& minPhi, double& minEta);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      // ----------member data ---------------------------
      // HSCP, dEdx and TOF collections
      edm::EDGetTokenT<vector<susybsm::HSCParticle>>           hscpToken_;
      edm::EDGetTokenT<edm::ValueMap<susybsm::HSCPIsolation>>  hscpIsoToken_;
      edm::EDGetTokenT<susybsm::MuonSegmentCollection>         muonSegmentToken_;
      //edm::EDGetTokenT<vector<reco::DeDxHitInfo>>      _dedxToken;
      edm::EDGetTokenT<reco::DeDxHitInfoAss>     dedxToken_;
      edm::EDGetTokenT<reco::MuonTimeExtraMap>   muonTimeToken_; // for reading inverse beta
      edm::EDGetTokenT<reco::MuonTimeExtraMap>   muonDtTimeToken_;
      edm::EDGetTokenT<reco::MuonTimeExtraMap>   muonCscTimeToken_;
      edm::EDGetTokenT<DTRecSegment4DCollection> muonDtSegmentToken_;
      edm::EDGetTokenT<CSCSegmentCollection>     muonCscSegmentToken_;
      edm::EDGetTokenT<vector<reco::Vertex>>     offlinePrimaryVerticesToken_;
      edm::EDGetTokenT<vector<reco::Track>>      refittedStandAloneMuonsToken_;
      edm::EDGetTokenT<reco::BeamSpot>           offlineBeamSpotToken_;
      edm::EDGetTokenT<vector<reco::Muon>>       muonToken_;
      edm::EDGetTokenT<edm::TriggerResults>      triggerResultsToken_;

      //edm::EDGetTokenT<reco::Track>  _tracksToken;//edm::EDGetTokenT<vector<reco::Track>>  _tracksToken;
      //edm::EDGetTokenT<vector<reco::DeDxHitInfo>>  _dedxHitInfosToken; //DataFormats/TrackReco/interface/DeDxHitInfo.h

      vector<double>  CutPt,      CutI,       CutTOF;
      vector<double>  CutPt_Flip, CutI_Flip,  CutTOF_Flip;
      //map<string, vector<double>> VCuts;

      map<string, TProfile*>  HCuts;

      bool* HSCPTk              ;
      bool* HSCPTk_SystP        ;
      bool* HSCPTk_SystI        ;
      bool* HSCPTk_SystT        ;
      bool* HSCPTk_SystM        ;
      bool* HSCPTk_SystPU       ;
      bool* HSCPTk_SystHUp      ;
      bool* HSCPTk_SystHDown    ;
      double* MaxMass           ;
      double* MaxMass_SystP     ;
      double* MaxMass_SystI     ;
      double* MaxMass_SystT     ;
      double* MaxMass_SystM     ;
      double* MaxMass_SystPU    ;
      double* MaxMass_SystHUp   ;
      double* MaxMass_SystHDown ;

      const reco::MuonTimeExtra* tof;
      const reco::MuonTimeExtra* dttof;
      const reco::MuonTimeExtra* csctof;

      double OpenAngle = -1; //global variable needed by PassPreselection... Ugly isn't it?!
      double TreeDXY = -1;
      double TreeDZ = -1;
      bool isCosmicSB = false;
      bool isSemiCosmicSB = false;

      int                TypeMode           =   0;      // 0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1
      int                SampleType         =   0;      // 0:Data, 1:MC, >=2:Signal

      bool SkipSelectionPlot                = false;

      // binning for the pT, mass, IP distributions
      double             PtHistoUpperBound   = 4000;
      double             MassHistoUpperBound = 4000;
      int                MassNBins           = 400;
      double             IPbound             = 1.0;
      int                PredBins            =   0;      //  How many different bins the prediction is split in for analysis being run, sets how many histograms are actually initialized.
      int                EtaBins             =  60;      //  How many bins we use for the background prediction method in Eta -- impacts background prediction method -- histograms with the name of the form "Pred_Eta*" in Analysis_PlotStructure.h
      double             dEdxS_UpLim         = 1.0;
      double             dEdxM_UpLim         = 30.0;
      int                DzRegions           =  6;     //Number of different Dz side regions used to make cosmic background prediction

      //Variables used in the TOF only HSCP search
      float              DTRegion        =   0.9;  //Define the dividing line between DT and 
      float              CSCRegion       =   0.9;  //CSC regions of CMS
      float              CosmicMinDz     =   70.;  //Min dz displacement to be tagged as cosmic muon
      float              CosmicMaxDz     =   120.; //Max dz displacement for cosmic tagged tracks
      double             minSegEtaSep    =  0.1;   //Minimum eta separation between SA track and muon segment on opposite side of detector

      int                minMuStations   =  2;

      // Thresholds for candidate preselection 
      double             GlobalMaxEta       =   2.1;    // cut on inner tracker track eta
      double             GlobalMaxV3D       =   99999;  //0.50;   // cut on 3D distance (cm) to closest vertex
      double             GlobalMaxDZ        =   0.50;   // cut on 1D distance (cm) to closest vertex in "Z" direction
      double             GlobalMaxDXY       =   0.50;   // cut on 2D distance (cm) to closest vertex in "R" direction
      double             GlobalMaxChi2      =   5.0;    // cut on Track maximal Chi2/NDF
      int                GlobalMinQual      =   2;      // cut on track quality (2 meaning HighPurity tracks)
      unsigned int       GlobalMinNOH       =   8;//7AMSB;      // cut on number of (valid) track pixel+strip hits 
      int                GlobalMinNOPH      =   2;      // cut on number of (valid) track pixel hits 
      double             GlobalMinFOVH      =   0.8;//0.0AMSB;    // cut on fraction of valid track hits
      unsigned int       GlobalMaxNOMHTillLast = 99999;//1AMSB;     // cut on the number of missing hits from IP till last hit (excluding hits behind the last hit)
      double             GlobalMinFOVHTillLast =-99999;//0.85AMSB;   // cut on the fraction of valid hits divided by total expected hits until the last one
      unsigned int       GlobalMinNOM       =   6;//7AMSB;      // cut on number of dEdx hits (generally equal to #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used)
      double             GlobalMinNDOF      =   8;      // cut on number of     DegreeOfFreedom used for muon TOF measurement
      double             GlobalMinNDOFDT    =   6;      // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
      double             GlobalMinNDOFCSC   =   6;      // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
      double             GlobalMaxTOFErr    =   0.15;//0.07;   // cut on error on muon TOF measurement
      double             GlobalMaxPterr     =   0.25;//0.50;//0.25;   // cut on error on track pT measurement 
      double             GlobalMaxTIsol     =  50;      // cut on tracker isolation (SumPt)
      double             GlobalMaxRelTIsol  =  9999999; // cut on relative tracker isolation (SumPt/Pt)
      double             GlobalMaxEIsol     =  0.30;    // cut on calorimeter isolation (E/P)
      double             GlobalMinPt        =  55.00;   // cut on pT    at PRE-SELECTION
      double             GlobalMinIs        =   0.0;    // cut on dEdxS at PRE-SELECTION (dEdxS is generally a  discriminator)
      double             GlobalMinIm        =   0.0;    // cut on dEdxM at PRE-SELECTION (dEdxM is generally an estimator    )
      double             GlobalMinTOF       =   1.0;    // cut on TOF   at PRE-SELECTION

      const int          MaxPredBins        =   6;      // The maximum number of different bins prediction is done in for any of the analyses (defines array size)
      

      double             dEdxK_Data    = 2.580;
      double             dEdxC_Data    = 3.922;
      double             dEdxK_MC      = 2.935;
      double             dEdxC_MC      = 3.197;

      //=============================================================
      Tuple* tuple;
      TupleSaver* tuple_saver;
      //=============================================================

      TH3F* dEdxTemplates = nullptr;
      double DeDxSF_0 = 1.00000; // [0]  unchanged
      double DeDxSF_1 = 1.41822; // [1]  Pixel data to SiStrip data
      double dEdxSF [2] = {
         DeDxSF_0,   
         DeDxSF_1    
      };
      double DeDxK = 0.0;
      double DeDxC = 0.0;

      dedxGainCorrector trackerCorrector;
      string DeDxTemplate;    // "MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", "Data13TeV16_dEdxTemplate.root"
      bool enableDeDxCalibration;
      string DeDxCalibration; //"Data13TeVGains_v2.root" if Data
      string Geometry;        //CMS_GeomTree.root
      string TimeOffset;      //MuonTimeOffset.txt
      muonTimingCalculator tofCalculator;

      // Emulators
      /*dedxHIPEmulator      HIPemulator;
      dedxHIPEmulator      HIPemulatorUp;
      dedxHIPEmulator      HIPemulatorDown;
      L1BugEmulator        L1Emul;
      HIPTrackLossEmulator HIPTrackLossEmul;*/

      bool useClusterCleaning;
      bool isData;
      bool isMC;
      bool isSignal;

      unsigned int CurrentRun = 0;

      float Event_Weight = 1;

      unsigned int TrigInfo =0; //1 -mu only, 2- met only, 3 mu and met 

      TRandom3* RNG = nullptr;
      bool is2016;
      bool is2016G;

      bool isMCglobal = false;

      double preTrackingChangeL1IntLumi       = 29679.982; // pb
      double IntegratedLuminosity             = 33676.4; //13TeV16

};
#endif               