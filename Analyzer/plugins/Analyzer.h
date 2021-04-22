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
//#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
//#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPDeDxInfo.h"
// Muons
#include "DataFormats/MuonReco/interface/Muon.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtra.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
// Muons CSC segments
//#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include <DataFormats/CSCRecHit/interface/CSCSegment.h>

#include "DataFormats/HepMCCandidate/interface/GenParticle.h"

// ~~~~~~~~~ user include files ~~~~~~~~~
#include "SUSYBSMAnalysis/Analyzer/interface/Analysis_Calibration.h"
#include "SUSYBSMAnalysis/Analyzer/interface/Analysis_Samples.h"
#include "SUSYBSMAnalysis/Analyzer/interface/Analysis_TOFUtility.h"
#include "SUSYBSMAnalysis/Analyzer/interface/Analysis_Tool.h"
//#include "SUSYBSMAnalysis/Analyzer/interface/Analysis_Tree.h"

using namespace std;

class Analyzer : public edm::EDAnalyzer {
   public:
      explicit Analyzer(const edm::ParameterSet&);
      ~Analyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);

      double scaleFactor(double eta);

      //bool PassPreselection(const susybsm::HSCParticle& hscp, const reco::DeDxHitInfo* dedxHits,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const edm::Event&, stPlots* st, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT, double MassErr);

      bool PassPreselection(
         const susybsm::HSCParticle& hscp, 
         const reco::DeDxHitInfo* dedxHits,  
         const reco::DeDxData* dedxSObj, 
         const reco::DeDxData* dedxMObj,
         const edm::Event& iEvent, 
         TTree* &tree,//stPlots* st, 
         const double& GenBeta, 
         bool RescaleP, 
         const double& RescaleI, 
         const double& RescaleT, 
         double MassErr);

      void InitTree(TTree* &tree, edm::Service<TFileService> &fs);
      void FillTree(TTree* &tree, const edm::Event& iEvent);
      void InitHist(edm::Service<TFileService> &fs);

   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
      
      // ----------member data ---------------------------
      edm::InputTag theInput;
      // HSCP, dEdx and TOF collections
      edm::EDGetTokenT<vector<susybsm::HSCParticle>>   _hscpToken;
      //edm::EDGetTokenT<vector<reco::DeDxHitInfo>>      _dedxToken;
      edm::EDGetTokenT<reco::DeDxHitInfoAss>   _dedxToken;
      /*edm::EDGetTokenT<edm::ValueMap<reco::MuonTimeExtra>> _tofToken;       // for reading inverse beta
      edm::EDGetTokenT<edm::ValueMap<reco::MuonTimeExtra>> _tofDtToken;
      edm::EDGetTokenT<edm::ValueMap<reco::MuonTimeExtra>> _tofCscToken;*/
      edm::EDGetTokenT<reco::MuonTimeExtraMap>   muonTimeToken_;
      edm::EDGetTokenT<reco::MuonTimeExtraMap>   muonDtTimeToken_;
      edm::EDGetTokenT<reco::MuonTimeExtraMap>   muonCscTimeToken_;
      edm::EDGetTokenT<DTRecSegment4DCollection> muonDtSegmentToken_;
      edm::EDGetTokenT<CSCSegmentCollection>     muonCscSegmentToken_;

      //edm::EDGetTokenT<reco::Track>  _tracksToken;//edm::EDGetTokenT<vector<reco::Track>>  _tracksToken;
      //edm::EDGetTokenT<vector<reco::DeDxHitInfo>>  _dedxHitInfosToken; //DataFormats/TrackReco/interface/DeDxHitInfo.h

      vector<double>  CutPt,      CutI,       CutTOF;
      vector<double>  CutPt_Flip, CutI_Flip,  CutTOF_Flip;

      map<string, TProfile*>  HCuts;
      map<string, TH2F*>  HMass;

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

      /*
      reco::MuonTimeExtra tof_;
      reco::MuonTimeExtra dttof_;
      reco::MuonTimeExtra csctof_;
      */
      reco::MuonTimeExtra tof;
      reco::MuonTimeExtra dttof;
      reco::MuonTimeExtra csctof;

      double OpenAngle = -1; //global variable needed by PassPreselection... Ugly isn't it?!
      double TreeDXY = -1;
      double TreeDZ = -1;
      bool isCosmicSB = false;
      bool isSemiCosmicSB = false;
      
      //unsigned int _minTracks;
      //TH1D *demohisto;
      //map<string, TProfile*> Histos;

      bool Debug;
      bool AddTree;

      // binning for the pT, mass, and IP distributions
      double             PtHistoUpperBound   = 4000;
      double             MassHistoUpperBound = 4000;
      int                MassNBins           = 400;
      double             IPbound             = 1.0;


      // Thresholds for candidate preselection --> note that some of the followings can be replaced at the beginning of Analysis_Step1_EventLoop function
      int                SampleType         =   0;      // 0:Data, 1:MC, >=2:Signal
      int                TypeMode           =   0;      // Tracker Only, Tracker+TOF, ...
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
      int                PredBins           =   0;      //  How many different bins the prediction is split in for analysis being run, sets how many histograms are actually initialized.

      double             dEdxK_Data    = 2.580;
      double             dEdxC_Data    = 3.922;
      double             dEdxK_MC      = 2.935;
      double             dEdxC_MC      = 3.197;

      string SampleTxtFile; //"Analysis_Samples.txt"

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
      string DeDxCalibration; //"Data13TeVGains_v2.root" if Data
      string Geometry;        //CMS_GeomTree.root
      string TimeOffset;      //MuonTimeOffset.txt
      muonTimingCalculator tofCalculator;

      // classes from Analysis_CommonFunction.h
      /*dedxHIPEmulator HIPemulator;
      dedxHIPEmulator HIPemulatorUp;
      dedxHIPEmulator HIPemulatorDown;*/
      /*L1BugEmulator        L1Emul;
      HIPTrackLossEmulator HIPTrackLossEmul;*/

      bool useClusterCleaning;
      bool isData;
      bool isMC;
      bool isSignal;

      unsigned int CurrentRun = 0;

      float Event_Weight = 1;

      TRandom3* RNG = nullptr;
      bool is2016;
      bool is2016G;

      //=============================================================
      //      
      //      Create Trees \& Branches
      //                  
      //=============================================================
      TTree*       Tree;
      /*
      unsigned int Run;
      unsigned int Event;
      unsigned int lumi;
      unsigned int bunch;
      */
      unsigned int NCuts;
      unsigned int Tree_Trig;
      unsigned int Tree_Run;
      unsigned int Tree_Event;
      unsigned int Tree_Lumi;
      unsigned int Tree_Hscp;
      float        Tree_Charge;
      float        Tree_Pt;
      float        Tree_PtErr;
      float        Tree_I;
      float        Tree_Ih;
      float        Tree_Ick;
      float        Tree_TOF;
      float        Tree_Mass;
      float        Tree_dZ;
      float        Tree_dXY;
      float        Tree_dR;
      float        Tree_eta;
      float        Tree_phi;
      unsigned int Tree_NOH; //number of (valid) track pixel+strip hits 
      unsigned int Tree_NOPH;//number of (valid) track pixel hits
      float        Tree_FOVH;//fraction of valid track hits
      unsigned int Tree_NOMH;//number of missing hits from IP till last hit (excluding hits behind the last hit)
      float        Tree_FOVHD;//fraction of valid hits divided by total expected hits until the last one
      unsigned int Tree_NOM;//number of dEdx hits (= #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used)
      float        Tree_Weight;
      float        Tree_GenId;
      float        Tree_GenCharge;
      float        Tree_GenMass;
      float        Tree_GenPt;
      float        Tree_GenEta;
      float        Tree_GenPhi;
   
      TTree*       GenTree;
      unsigned int GenTree_Run;
      unsigned int GenTree_Event;
      unsigned int GenTree_Lumi;
      unsigned int GenTree_Hscp;
      float        GenTree_Weight;
      float        GenTree_GenId;
      float        GenTree_GenCharge;
      float        GenTree_GenMass;
      float        GenTree_GenPt;
      float        GenTree_GenEta;
      float        GenTree_GenPhi;

      //=============================================================
      //      
      //      Create Histograms
      //                  
      //=============================================================
      TH2F*  Mass;
      TH2F*  MassTOF;
      TH2F*  MassComb;
      TH2F*  MaxEventMass;

      TH2F*  Mass_SystP;
      TH2F*  MassTOF_SystP;
      TH2F*  MassComb_SystP;
      TH2F*  MaxEventMass_SystP;

      TH2F*  Mass_SystI;
      TH2F*  MassTOF_SystI;
      TH2F*  MassComb_SystI;
      TH2F*  MaxEventMass_SystI;

      TH2F*  Mass_SystM;
      TH2F*  MassTOF_SystM;
      TH2F*  MassComb_SystM;
      TH2F*  MaxEventMass_SystM;

      TH2F*  Mass_SystT;
      TH2F*  MassTOF_SystT;
      TH2F*  MassComb_SystT;
      TH2F*  MaxEventMass_SystT;

      TH2F*  Mass_SystPU;
      TH2F*  MassTOF_SystPU;
      TH2F*  MassComb_SystPU;
      TH2F*  MaxEventMass_SystPU;

      TH2F*  Mass_SystHUp;
      TH2F*  MassTOF_SystH;
      TH2F*  MassComb_SystHUp;
      TH2F*  MaxEventMass_SystHUp;

      TH2F*  Mass_SystHDown;
      TH2F*  MassComb_SystHDown;
      TH2F*  MaxEventMass_SystHDown;



      TH2F*  Mass_Flip;
      TH2F*  MassTOF_Flip;
      TH2F*  MassComb_Flip;

      TProfile* IntLumi;
      TH1F* TotalE;
      TH1F* TotalEPU; 
      TH1F* TotalTE;
      TH1F* Total;
      TH1F* V3D; 
      TH1F* Chi2;  
      TH1F* Qual; 
      TH1F* TNOH;
      TH1F* TNOM;
      TH1F* nDof;
      TH1F* tofError;
      TH1F* Pterr;
      TH1F* MPt; 
      TH1F* MI; 
      TH1F* MTOF; 
      TH1F* TIsol;
      TH1F* EIsol;
      TH1F* SumpTOverpT;
      TH1F* Pt;	
      TH1F* I;	
      TH1F* TOF;
      TH1F* HSCPE;
      TH1F* NVTrack;
      TH1F* Stations;
      TH1F* Dxy;
      TH1F* Dz;
      TH1F* SegSep;
      TH1F* FailDz;
      TH1F* Basic;
      
      TH1F* HSCPE_SystP;
      TH1F* HSCPE_SystI;
      TH1F* HSCPE_SystM;
      TH1F* HSCPE_SystT;
      TH1F* HSCPE_SystPU;
      TH1F* HSCPE_SystHUp;
      TH1F* HSCPE_SystHDown;

      TH1F* Gen_DecayLength;
      TH1F* Beta_Gen;
      TH1F* Beta_GenCharged;
      TH1F* Beta_Triggered;
      TH1F* Beta_Matched;
      TH1F* Beta_PreselectedA;
      TH1F* Beta_PreselectedB;
      TH1F* Beta_PreselectedC;
      TH2F* Beta_SelectedP;
      TH2F* Beta_SelectedI;
      TH2F* Beta_SelectedT;

      TH1F*  BS_V3D;
      TH1F*  BS_Chi2;
      TH1F*  BS_Qual;
      TH1F*  BS_TNOH;
      TH1F*  BS_TNOH_PUA;
      TH1F*  BS_TNOH_PUB;
      TH1F*  BS_TNOHFraction;
      TH1F*  BS_TNOPH;
      TH1F*  BS_TNOHFractionTillLast;
      TH1F*  BS_TNOMHTillLast;
      TH1F*  BS_Eta;
      TH1F*  BS_TNOM;
      TH1F*  BS_TNOM_PUA;
      TH1F*  BS_TNOM_PUB;
      TProfile*  BS_NOMoNOHvsPV;
      TH1F*  BS_nDof;
      TH1F*  BS_TOFError;
      TH1F*  BS_Pterr;
      TH1F*  BS_MPt; 
      TH1F*  BS_MIs; 
      TH1F*  BS_MIm; 
      TH1F*  BS_MTOF;
      TH1F*  BS_TIsol;
      TH1F*  BS_EIsol;
      TH1F*  BS_SumpTOverpT;
      TH1F*  BS_dR_NVTrack;
      TH1F*  BS_MatchedStations;
      TH1F*  BS_NVertex;
      TH1F*  BS_NVertex_NoEventWeight;
      TH1F*  BS_PV;
      TH1F*  BS_PV_NoEventWeight;
      TH1F*  BS_dzAll;
      TH1F*  BS_dxyAll;
      TH1F*  BS_dzMinv3d;
      TH1F*  BS_dxyMinv3d;
      TH1F*  BS_SegSep;
      TH1F*  BS_SegMinPhiSep;
      TH1F*  BS_SegMinEtaSep;
      TH1F*  BS_SegMinEtaSep_FailDz;
      TH1F*  BS_SegMinEtaSep_PassDz;
      TH1F*  BS_Dz_FailSep;
      TH1F*  BS_InnerInvPtDiff;
      TH1F*  BS_Phi;
      TH1F*  BS_TimeAtIP;
      TH1F*  BS_OpenAngle;
      TH1F*  BS_OpenAngle_Cosmic;

      TH1F*  BS_Pt_FailDz;
      TH1F*  BS_Pt_FailDz_DT;
      TH1F*  BS_Pt_FailDz_CSC;
      TH1F*  BS_TOF_FailDz;
      TH1F*  BS_TOF_FailDz_DT;
      TH1F*  BS_TOF_FailDz_CSC;
      TH1F*  BS_Dxy;
      TH1F*  BS_Dxy_Cosmic;
      TH1F*  BS_Dz;
      TH1F*  BS_Dz_Cosmic;
      TH1F*  BS_Dz_CSC;
      TH1F*  BS_Dz_DT;
      //TH1F*  BS_Pt_Binned[MaxPredBins];// error: invalid use of non-static data member 'Analyzer::MaxPredBins'
      //TH1F*  BS_TOF_Binned[MaxPredBins]; 

      TH1F*  BS_LastHitDXY;
      TH1F*  BS_LastHitD3D;

      TH2F* AS_Eta_RegionA;
      TH2F* AS_Eta_RegionB;
      TH2F* AS_Eta_RegionC;
      TH2F* AS_Eta_RegionD;
      TH2F* AS_Eta_RegionE;
      TH2F* AS_Eta_RegionF;
      TH2F* AS_Eta_RegionG;
      TH2F* AS_Eta_RegionH;

      TH1F*  BS_P; 	   TH2F*  AS_P;
      TH1F*  BS_Pt;	   TH2F*  AS_Pt;
      TH1F*  BS_Pt_PUA;
      TH1F*  BS_Pt_PUB;
      TH1F*  BS_Pt_DT;
      TH1F*  BS_Pt_CSC;
      TH1F*  BS_Is;	   TH2F*  AS_Is;
      TH1F*  BS_Is_PUA;
      TH1F*  BS_Is_PUB;
      TH1F*  BS_Im;           TH2F*  AS_Im;
      TH1F*  BS_Im_PUA;
      TH1F*  BS_Im_PUB;
      TH1F*  BS_TOF;          TH2F*  AS_TOF;
      TH1F*  BS_TOF_PUA;
      TH1F*  BS_TOF_PUB;
      TH1F*  BS_TOF_DT;
      TH1F*  BS_TOF_CSC;
      TH1F*  BS_Is_Cosmic;
      TH1F*  BS_Pt_Cosmic;



      TH2F*  BS_EtaIs;        //TH3F*  AS_EtaIs;
      TH2F*  BS_EtaIm;        //TH3F*  AS_EtaIm;
      TH2F*  BS_EtaP;	   //TH3F*  AS_EtaP;
      TH2F*  BS_EtaPt;	   //TH3F*  AS_EtaPt;
      TH2F*  BS_EtaTOF;       //TH3F*  AS_EtaTOF;
      TH2F*  BS_EtaDz;
      TH2F*  BS_EtaNBH;       // number of bad hits vs Eta


      TH2F*  BS_PIs;	   TH3F*  AS_PIs;
      TH2F*  BS_PImHD; 
      TH2F*  BS_PIm;          TH3F*  AS_PIm;
      TH2F*  BS_PtIs;         TH3F*  AS_PtIs;
      TH2F*  BS_PtIm;         TH3F*  AS_PtIm;
      TH2F*  BS_PtTOF;
      TH2F*  BS_TOFIs;        TH3F*  AS_TOFIs;  
      TH2F*  BS_TOFIm;        TH3F*  AS_TOFIm;   

      //Prediction histograms
      TH1D* H_A;
      TH1D* H_B;
      TH1D* H_C;
      TH1D* H_D;
      TH1D* H_E;
      TH1D* H_F;
      TH1D* H_G;
      TH1D* H_H;
};
#endif