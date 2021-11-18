
#include <exception>
#include <vector>
#include <algorithm>
#include <fstream>
#include <iostream>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TCutG.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TLorentzVector.h"

namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "DataFormats/Scalers/interface/LumiScalers.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;


//#include "../../AnalysisCode_NewSyst_Hybr_WOverflow/Analysis_Step1_EventLoop.C"
#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"
#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"


#endif


double DistToHSCP (const reco::TrackRef& track, const std::vector<reco::GenParticle>& genColl);
bool isCompatibleWithCosmic (const reco::TrackRef& track, const std::vector<reco::Vertex>& vertexColl);
double GetMass (double P, double I, double K, double C);

const double L_Min               = 1   ;
const double L_Max               = 22  ; 
const int    L_NBins             = 21  ; 
const double P_Min               = 1   ;
const double P_Max               = 31  ; // 1 + 14 + 1; final one is for pixel!
const int    P_NBins             = 30  ; // 15th bin = pixel; 0 is underflow
const double Path_Min            = 0.2 ;
const double Path_Max            = 1.6 ;
const int    Path_NBins          = 42  ;
const double Charge_Min          = 0   ;
const double Charge_Max          = 5000;
const int    Charge_NBins        = 500 ;
const double Eta_Min             = -5  ;
const double Eta_Max             = 5   ;
const int    Eta_NBins           = 1000;
const double Phi_Min             = -3.5 ;
const double Phi_Max             = 3.5  ;
const int    Phi_NBins           = 700;
const double Momentum_Min        = 0   ; 
const double Momentum_Max        = 10  ;
const int    Momentum_NBins      = 1000;



struct dEdxStudyObj
{
   string Name;
   bool isDiscrim;
   bool isEstim;
   bool useTrunc;
   bool isHit;

   bool usePixel;
   bool useStrip;

   bool mustBeInside;
   bool removeCosmics;
   bool correctFEDSat;
   bool useClusterCleaning;
   bool fakeHIP;

   double dropLowerDeDxValue;
   bool trimPixel;

   bool skipPixelL1;

   double Kconst;
   double Cconst;

   int crossTalkInvAlgo; // 0  -- do not use crossTalkInversion
                         // 1  -- use existing algorithm developed by Claude

   TH3D* MG_Charge_Vs_Path;
   TH3D* L_Charge_Vs_Path;
   TH1D* HdedxMIP;
   TH1D* HdedxMIP4;
   TH1D* HdedxMIP8;
   TH1D* HdedxMIP12;
   TH2D* HdedxVsP;
   TH2D* HdedxVsPSyst;
   TH1D* HStripADC;
   TH2D* HStripADCvsP;
//   TH2D* HdedxVsQP;
//   TProfile2D* HdedxVsP_NS;
   TProfile* HdedxVsPProfile;
   TProfile* HdedxVsEtaProfile;
   TProfile* HdedxVsNOH;
   TProfile* HNOMVsdEdxProfile;
   TH2D* HdedxVsEta;
   TH2D* HNOMVsdEdx;
   TProfile* HNOSVsEtaProfile;
   TProfile* HNOMVsEtaProfile;
   TProfile* HNOMSVsEtaProfile;
   TH1D* HMass;
   TH1D* HMassHSCP;
   TH1D* HP;
   TH1D* HHit;
   TH1D* HHit_U;
   TProfile* HHitProfile; 
   TProfile* HHitProfile_U; 
   TH1D* HProtonHitSO; 
   TH1D* HProtonHitPO; 
   TProfile* Charge_Vs_FS[16];
   TH2D* Charge_Vs_XYH[16];
   TH2D* Charge_Vs_XYL[16];
   TH2D* Charge_Vs_XYHN[16];
   TH2D* Charge_Vs_XYLN[16];

//   TH3D* HdedxVsPVsEta;

   TH2D* HdedxVsP_wPixels;
   TH2D* HdedxVsP_wPixels_noL1;

   TH2D* HdedxVsP_Eta0p91;
   TH2D* HdedxVsP_0p91Eta1p74;
   TH2D* HdedxVsP_1p74Eta;

   TH2D* HdedxVsP_PU_0_10;
   TH2D* HdedxVsP_PU_10_15;
   TH2D* HdedxVsP_PU_15_20;
   TH2D* HdedxVsP_PU_20_25;
   TH2D* HdedxVsP_PU_25_30;
   TH2D* HdedxVsP_PU_30_35;
   TH2D* HdedxVsP_PU_35_40;
   TH2D* HdedxVsP_PU_40_50;
   TH2D* HdedxVsP_PU_50_inf;

   TH2D* HdedxVsP_instLumi_0_1000;
   TH2D* HdedxVsP_instLumi_1000_2000;
   TH2D* HdedxVsP_instLumi_2000_3000;
   TH2D* HdedxVsP_instLumi_3000_4000;
   TH2D* HdedxVsP_instLumi_4000_5000;
   TH2D* HdedxVsP_instLumi_5000_6000;
   TH2D* HdedxVsP_instLumi_6000_7000;
   TH2D* HdedxVsP_instLumi_7000_8000;
   TH2D* HdedxVsP_instLumi_8000_9000;
   TH2D* HdedxVsP_instLumi_9000_10000;
   TH2D* HdedxVsP_instLumi_10000_11000;
   TH2D* HdedxVsP_instLumi_11000_12000;
   TH2D* HdedxVsP_instLumi_12000_13000;
   TH2D* HdedxVsP_instLumi_13000_14000;
   TH2D* HdedxVsP_instLumi_14000_15000;
   TH2D* HdedxVsP_instLumi_15000_16000;
   TH2D* HdedxVsP_instLumi_16000_inf;

   TH3D* MG_Charge_Vs_Path_instLumi_0_5000;
   TH3D* MG_Charge_Vs_Path_instLumi_5000_10000;
   TH3D* MG_Charge_Vs_Path_instLumi_10000_15000;
   TH3D* MG_Charge_Vs_Path_instLumi_15000_inf;

   TH3D* MG_Charge_Vs_Path_PU_0_25;
   TH3D* MG_Charge_Vs_Path_PU_25_30;
   TH3D* MG_Charge_Vs_Path_PU_30_35;
   TH3D* MG_Charge_Vs_Path_PU_35_40;
   TH3D* MG_Charge_Vs_Path_PU_40_inf;

   TH3D* MG_Charge_Vs_Path_wPixels;

   
//adds by Dylan - dec2020 janv2021
   TH3D* MG_PathlengthVsEta; //only Ias
   TH3D* L_PathlengthVsEta;
   
   TH2D* HdedxVsP_PUlower20; //Ih
   TH2D* HdedxVsP_PUhigher35;

/*   TH2D* HdedxVsP_Eta0p1;
   TH2D* HdedxVsP_0p1Eta0p2;
   TH2D* HdedxVsP_0p2Eta0p3;
   TH2D* HdedxVsP_0p3Eta0p4;
   TH2D* HdedxVsP_0p4Eta0p5;
   TH2D* HdedxVsP_0p5Eta0p6;
   TH2D* HdedxVsP_0p6Eta0p7;
   TH2D* HdedxVsP_0p7Eta0p8;
   TH2D* HdedxVsP_0p8Eta0p9;
   TH2D* HdedxVsP_0p9Eta1p0;
   TH2D* HdedxVsP_1p0Eta1p1;
   TH2D* HdedxVsP_1p1Eta1p2;
   TH2D* HdedxVsP_1p2Eta1p3;
   TH2D* HdedxVsP_1p3Eta1p4;
   TH2D* HdedxVsP_1p4Eta1p5;
   TH2D* HdedxVsP_1p5Eta1p6;
   TH2D* HdedxVsP_1p6Eta1p7;
   TH2D* HdedxVsP_1p7Eta1p8;
   TH2D* HdedxVsP_1p8Eta1p9;
   TH2D* HdedxVsP_1p9Eta2p0;
   TH2D* HdedxVsP_2p0Eta2p1;*/

   TH2F* HdedxVsEta_highdEdx;
   TH2F* HdedxVsEta_0p9P1p0;
   TH2F* HdedxVsEta_5p0P;
   TH2F* HdedxVsAbsEta_5p0P;


   TH3F* dEdxTemplates = NULL;
   std::unordered_map<unsigned int,double>* TrackerGains = NULL;

   dEdxStudyObj(string Name_, int type_, int subdet_, TH3F* dEdxTemplates_=NULL, double K_=2.7, double C_=3.2, std::unordered_map<unsigned int,double>* TrackerGains_=NULL, bool mustBeInside_=false, bool removeCosmics_=false, bool correctFEDSat_=false, int useClusterCleaning_=0, int crossTalkInvAlgo_=0, double dropLowerDeDxValue_ = 0, bool trimPixel_ = false, bool fakeHIP_=true, bool skipPixelL1_=false){
      Name = Name_;

      if     (type_==0){ isHit=true;  isEstim= false; isDiscrim = false; useTrunc = false;} // hit level only
      else if(type_==1){ isHit=false; isEstim= true;  isDiscrim = false; useTrunc = false;} // harm2
      else if(type_==2){ isHit=false; isEstim= false; isDiscrim = true;  useTrunc = false;} // Ias via harm2
      else if(type_==3){ isHit=false; isEstim= true;  isDiscrim = false; useTrunc = true; } // trunc40
      else             { isHit=false; isEstim= false; isDiscrim = false;}

           if(subdet_==1){ usePixel = true;  useStrip = false;}
      else if(subdet_==2){ usePixel = false; useStrip = true; }
      else               { usePixel = true;  useStrip = true; }

      dEdxTemplates      = dEdxTemplates_;
      TrackerGains       = TrackerGains_;
      mustBeInside       = mustBeInside_;
      removeCosmics      = removeCosmics_; 
      correctFEDSat      = correctFEDSat_;
      useClusterCleaning = static_cast<bool> (useClusterCleaning_);
      crossTalkInvAlgo   = crossTalkInvAlgo_;
      dropLowerDeDxValue = dropLowerDeDxValue_;
      trimPixel          = trimPixel_;
      fakeHIP            = fakeHIP_;
      skipPixelL1        = skipPixelL1_;

      Kconst = K_;
      Cconst = C_;

      string HistoName;
      //HitLevel plot      
      if(isHit){ 
         HistoName = Name + "_Hit";               HHit                  = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20); 
         HistoName = Name + "_Hit_U";             HHit_U                = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20); 
         HistoName = Name + "_HitProfile";        HHitProfile           = new TProfile(  HistoName.c_str(), HistoName.c_str(),  50, 0, 100); 
         HistoName = Name + "_HitProfile_U";      HHitProfile_U         = new TProfile(  HistoName.c_str(), HistoName.c_str(),  50, 0, 100);
         HistoName = Name + "_StripADC";          HStripADC             = new TH1D(      HistoName.c_str(), HistoName.c_str(), 256, 0, 256);
         HistoName = Name + "_StripADCvsP";       HStripADCvsP          = new TH2D(      HistoName.c_str(), HistoName.c_str(), 500, 0, 10, 256, 0, 256);
         if(usePixel && useStrip){
            for(unsigned int g=0;g<16;g++){
               char Id[255]; sprintf(Id, "%02i", g);
               HistoName = Name + "_ChargeVsFS"+Id;       Charge_Vs_FS[g]       = new TProfile  ( HistoName.c_str(), HistoName.c_str(),  769, 0, 769);
               HistoName = Name + "_ChargeVsXYH"+Id;      Charge_Vs_XYH[g]      = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -15, 15, 250, -15, 15);
               HistoName = Name + "_ChargeVsXYL"+Id;      Charge_Vs_XYL[g]      = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -15, 15, 250, -15, 15);
               HistoName = Name + "_ChargeVsXYHN"+Id;     Charge_Vs_XYHN[g]     = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -1.5, 1.5, 250, -1.5, 1.5);
               HistoName = Name + "_ChargeVsXYLN"+Id;     Charge_Vs_XYLN[g]     = new TH2D      ( HistoName.c_str(), HistoName.c_str(),  250, -1.5, 1.5, 250, -1.5, 1.5);
            }
         }
      }

      //Track Level plots
      if(isEstim || isDiscrim){
         HistoName = Name + "_MIP";               HdedxMIP              = new TH1D(      HistoName.c_str(), HistoName.c_str(), 1000, 0, isDiscrim?1.0:25);
         HistoName = Name + "_MIP4";              HdedxMIP4             = new TH1D(      HistoName.c_str(), HistoName.c_str(), 1000, 0, isDiscrim?1.0:25);
         HistoName = Name + "_MIP8";              HdedxMIP8             = new TH1D(      HistoName.c_str(), HistoName.c_str(), 1000, 0, isDiscrim?1.0:25);
         HistoName = Name + "_MIP12";             HdedxMIP12            = new TH1D(      HistoName.c_str(), HistoName.c_str(), 1000, 0, isDiscrim?1.0:25);
         HistoName = Name + "_dedxVsP";           HdedxVsP              = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsPSyst";       HdedxVsPSyst          = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
//       HistoName = Name + "_dedxVsQP";          HdedxVsQP             = new TH2D(      HistoName.c_str(), HistoName.c_str(), 6000, -30, 30,1500,0, isDiscrim?1.0:15);
//       HistoName = Name + "_dedxVsP_NS";        HdedxVsP_NS           = new TProfile2D(HistoName.c_str(), HistoName.c_str(), 3000, 0, 30,1500,0, isDiscrim?1.0:15);
         HistoName = Name + "_Profile";           HdedxVsPProfile       = new TProfile(  HistoName.c_str(), HistoName.c_str(),   50, 0,100);
         HistoName = Name + "_Eta";               HdedxVsEtaProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_dedxVsNOH";         HdedxVsNOH            = new TProfile(  HistoName.c_str(), HistoName.c_str(),   80, 0, 80);
         HistoName = Name + "_NOMVsdEdxProfile";  HNOMVsdEdxProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),   200, 0, isDiscrim?1.0:25);
         HistoName = Name + "_NOMVsdEdx";         HNOMVsdEdx            = new TH2D(      HistoName.c_str(), HistoName.c_str(), 200, 0, isDiscrim?1.0:25, 30, 0, 30);
         HistoName = Name + "_Eta2D";             HdedxVsEta            = new TH2D(      HistoName.c_str(), HistoName.c_str(),   60,-3,  3, 100,0, isDiscrim?1.0:5);
         HistoName = Name + "_NOS";               HNOSVsEtaProfile      = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_NOM";               HNOMVsEtaProfile      = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_NOMS";              HNOMSVsEtaProfile     = new TProfile(  HistoName.c_str(), HistoName.c_str(),   60,-3,  3);
         HistoName = Name + "_P";                 HP                    = new TH1D(      HistoName.c_str(), HistoName.c_str(),   50, 0, 100);  
         //HistoName = Name + "_dedxVsPVsEta"; HdedxVsPVsEta              = new TH3D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15,Eta_NBins, Eta_Min, Eta_Max);

         HistoName = Name + "_dedxVsP_Eta0p91"; HdedxVsP_Eta0p91        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_0p91Eta1p74"; HdedxVsP_0p91Eta1p74= new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_1p74Eta"; HdedxVsP_1p74Eta        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_MG_ChargeVsPath";      MG_Charge_Vs_Path        = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_L_ChargeVsPath";      L_Charge_Vs_Path        = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         

         HistoName = Name + "_MG_ChargeVsPath_instLumi_0_5000";      MG_Charge_Vs_Path_instLumi_0_5000 = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_MG_ChargeVsPath_instLumi_5000_10000";      MG_Charge_Vs_Path_instLumi_5000_10000 = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_MG_ChargeVsPath_instLumi_10000_15000";      MG_Charge_Vs_Path_instLumi_10000_15000 = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_MG_ChargeVsPath_instLumi_15000_inf";      MG_Charge_Vs_Path_instLumi_15000_inf = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         
         
         HistoName = Name + "_MG_ChargeVsPath_PU_0_25";      MG_Charge_Vs_Path_PU_0_25 = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_MG_ChargeVsPath_PU_25_30";      MG_Charge_Vs_Path_PU_25_30 = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_MG_ChargeVsPath_PU_30_35";      MG_Charge_Vs_Path_PU_30_35 = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_MG_ChargeVsPath_PU_35_40";      MG_Charge_Vs_Path_PU_35_40 = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
         HistoName = Name + "_MG_ChargeVsPath_PU_40_inf";      MG_Charge_Vs_Path_PU_40_inf = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
        

         HistoName = Name + "_MG_PathlengthVsEta"; MG_PathlengthVsEta = new TH3D(HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Eta_NBins, Eta_Min, Eta_Max);
         HistoName = Name + "_L_PathlengthVsEta"; L_PathlengthVsEta = new TH3D(HistoName.c_str(), HistoName.c_str(), L_NBins, L_Min, L_Max, Path_NBins, Path_Min, Path_Max, Eta_NBins, Eta_Min, Eta_Max);
           
         HistoName = Name + "_dedxVsP_PUlower20"; HdedxVsP_PUlower20    = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PUhigher35";HdedxVsP_PUhigher35   = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         

         HistoName = Name + "_dedxVsP_PU_0_10"; HdedxVsP_PU_0_10        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_10_15"; HdedxVsP_PU_10_15        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_15_20"; HdedxVsP_PU_15_20        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_20_25"; HdedxVsP_PU_20_25        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_25_30"; HdedxVsP_PU_25_30        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_30_35"; HdedxVsP_PU_30_35        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_35_40"; HdedxVsP_PU_35_40        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_40_50"; HdedxVsP_PU_40_50        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_PU_50_inf"; HdedxVsP_PU_50_inf        = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsEta_highdEdx";             
            HdedxVsEta_highdEdx            = new TH2F(      HistoName.c_str(), HistoName.c_str(),   Eta_NBins, Eta_Min, Eta_Max, 100,5, 35);
         HistoName = Name + "_dedxVsEta_0p9P1p0";             
            HdedxVsEta_0p9P1p0            = new TH2F(      HistoName.c_str(), HistoName.c_str(),   Eta_NBins,Eta_Min,  Eta_Max, 100,0,20);
         HistoName = Name + "_dedxVsEta_5p0P";             
            HdedxVsEta_5p0P            = new TH2F(      HistoName.c_str(), HistoName.c_str(),   Eta_NBins,Eta_Min,  Eta_Max, 100,0,20);
         HistoName = Name + "_dedxVsAbsEta_5p0P";             
            HdedxVsAbsEta_5p0P            = new TH2F(      HistoName.c_str(), HistoName.c_str(),   Eta_NBins,0,  Eta_Max, 100,0,20);


         HistoName = Name + "_dedxVsP_instLumi_0_1000"; HdedxVsP_instLumi_0_1000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_1000_2000"; HdedxVsP_instLumi_1000_2000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_2000_3000"; HdedxVsP_instLumi_2000_3000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_3000_4000"; HdedxVsP_instLumi_3000_4000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_4000_5000"; HdedxVsP_instLumi_4000_5000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_5000_6000"; HdedxVsP_instLumi_5000_6000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_6000_7000"; HdedxVsP_instLumi_6000_7000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_7000_8000"; HdedxVsP_instLumi_7000_8000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_8000_9000"; HdedxVsP_instLumi_8000_9000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_9000_10000"; HdedxVsP_instLumi_9000_10000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_10000_11000"; HdedxVsP_instLumi_10000_11000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_11000_12000"; HdedxVsP_instLumi_11000_12000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_12000_13000"; HdedxVsP_instLumi_12000_13000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_13000_14000"; HdedxVsP_instLumi_13000_14000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_14000_15000"; HdedxVsP_instLumi_14000_15000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_15000_16000"; HdedxVsP_instLumi_15000_16000 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
         HistoName = Name + "_dedxVsP_instLumi_16000_inf"; HdedxVsP_instLumi_16000_inf = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);



/*            HistoName = Name + "_dedxVsP_Eta0p1";           
                HdedxVsP_Eta0p1                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15); 
            HistoName = Name + "_dedxVsP_0p1Eta0p2";           
                HdedxVsP_0p1Eta0p2                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15); 
            HistoName = Name + "_dedxVsP_0p2Eta0p3";           
                HdedxVsP_0p2Eta0p3                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_0p3Eta0p4";           
                HdedxVsP_0p3Eta0p4                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_0p4Eta0p5";           
                HdedxVsP_0p4Eta0p5                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_0p5Eta0p6";           
                HdedxVsP_0p5Eta0p6                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_0p6Eta0p7";           
                HdedxVsP_0p6Eta0p7                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_0p7Eta0p8";           
                HdedxVsP_0p7Eta0p8                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_0p8Eta0p9";           
                HdedxVsP_0p8Eta0p9                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_0p9Eta1p0";           
                HdedxVsP_0p9Eta1p0                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p0Eta1p1";           
                HdedxVsP_1p0Eta1p1                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p1Eta1p2";           
                HdedxVsP_1p1Eta1p2                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p2Eta1p3";           
                HdedxVsP_1p2Eta1p3                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p3Eta1p4";           
                HdedxVsP_1p3Eta1p4                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p4Eta1p5";           
                HdedxVsP_1p4Eta1p5                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p5Eta1p6";           
                HdedxVsP_1p5Eta1p6                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p6Eta1p7";           
                HdedxVsP_1p6Eta1p7                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p7Eta1p8";           
                HdedxVsP_1p7Eta1p8                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p8Eta1p9";           
                HdedxVsP_1p8Eta1p9                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_1p9Eta2p0";           
                HdedxVsP_1p9Eta2p0                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
            HistoName = Name + "_dedxVsP_2p0Eta2p1";           
                HdedxVsP_2p0Eta2p1                = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);*/



            if(usePixel)
            {
                HistoName = Name + "_dedxVsP_wPixels"; HdedxVsP_wPixels = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
                HistoName = Name + "_dedxVsP_wPixels_noL1"; HdedxVsP_wPixels_noL1 = new TH2D(      HistoName.c_str(), HistoName.c_str(),  500, 0, 10,1000,0, isDiscrim?1.0:15);
                HistoName = Name + "_MG_ChargeVsPath_wPixels";      MG_Charge_Vs_Path_wPixels        = new TH3D(      HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
           
            }
      }

      //estimator plot only
      if(isEstim){
         HistoName = Name + "_Mass";              HMass                 = new TH1D(      HistoName.c_str(), HistoName.c_str(),  250, 0, 10);
         HistoName = Name + "_MassHSCP";          HMassHSCP             = new TH1D(      HistoName.c_str(), HistoName.c_str(),  300, 0, 3000);
//         HistoName = Name + "_ProtonHit";         HProtonHit            = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20);
         HistoName = Name + "_ProtonHitPO";       HProtonHitPO          = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20);
         HistoName = Name + "_ProtonHitSO";       HProtonHitSO          = new TH1D(      HistoName.c_str(), HistoName.c_str(),  200, 0, 20);
      }
      
         

         


   }
};

void DeDxStudy(string DIRNAME="COMPILE", string INPUT="dEdx.root", string OUTPUT="out.root")
{
  if(DIRNAME=="COMPILE") return;

   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.07);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505,"X");
   TH1::AddDirectory(kTRUE);

   bool isData   = !(INPUT.find("MC")!=string::npos);
   bool isSignal = false;
   bool removeCosmics = true;
   std::vector<string> FileName;
   if(INPUT.find(".root")<std::string::npos){
      char* pch=strtok(&INPUT[0],",");
      while (pch!=NULL){
         FileName.push_back(pch);    
         pch=strtok(NULL,",");
      }
   }else{
      string SampleId = INPUT;
      InitBaseDirectory();
      GetSampleDefinition(samples , DIRNAME+"/../../AnalysisCode/Analysis_Samples.txt");
      stSample& sample = samples[JobIdToIndex(SampleId, samples)];
      isData   = (sample.Type==0);
      isSignal = (sample.Type==2);
      GetInputFiles(sample, BaseDirectory, FileName, 0);
   }

   unsigned int CurrentRun = 1;
   dedxGainCorrector trackerCorrector;
//   dedxHIPEmulator   HIPEmulator (true, "pixelDeDxStudy", "stripDeDxStudy");
//   HIPEmulator.setPeriodHIPRate(true);
   HIPemulator.setPeriodHIPRate(true);
//   printf("I am here!\n");
   TH3F* dEdxTemplates      = NULL;
   TH3F* dEdxTemplatesIn    = NULL;
   TH3F* dEdxTemplatesInc   = NULL;
   TH3F* dEdxTemplatesCCC   = NULL;
   TH3F* dEdxTemplatesCCC16 = NULL;
   TH3F* dEdxTemplatesCC    = NULL;
   TH3F* dEdxTemplatesCI    = NULL;
//   TH3F* dEdxTemplatesnewCCC= NULL;
   double dEdx_U [2] = {1.0, 1.0};

//   double dEdxSF_OldCC[2] = {1.0,1.34475};
//   double dEdxSF_NewCC[2] = {1.0,1.21845};


   bool SuppressFakeHIP = false;
   if(isData){
         SuppressFakeHIP = true; // never fake HIPs for Data -- they are already present as it is
         dEdxSF [0] = 1.00000;
         dEdxSF [1] = 1.00000;
         //dEdxSF [1] = 1.41822*1.0371500*0.99565; //preG
	     //dEdxSF [1] = 1.41822*1.1265500*1.00815;  //postG
//         dEdxTemplates    = loadDeDxTemplate(DIRNAME + "/../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", true);
//         dEdxTemplatesInc = loadDeDxTemplate(DIRNAME + "/../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", false);
//         dEdxTemplates      = loadDeDxTemplate (DIRNAME+"/Templates/dEdxTemplate_hit_SP_Data.root"           , true);
//         dEdxTemplatesIn    = loadDeDxTemplate (DIRNAME+"/Templates/dEdxTemplate_hit_SP_in_noC_Data.root"    , true);
//         dEdxTemplatesInc   = loadDeDxTemplate (DIRNAME+"/Templates/dEdxTemplate_hit_SP_Data.root"           , false);
//         dEdxTemplatesCC    = loadDeDxTemplate (DIRNAME+"/../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CC.root" , true);
//         dEdxTemplatesCI    = loadDeDxTemplate (DIRNAME+"/../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CI.root" , true);
         dEdxTemplatesCCC   = loadDeDxTemplate (DIRNAME+"/../../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true);
         //dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/dEdxTemplate_hit_SP_in_noC_CCC_Run278018.root", true);// commented from Joze code
	 /////	 dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/Data13TeV16_dEdxTemplate.root", true);  //tentativo a caso
	 //////dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_Run278018.root", true);
	 dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_RunPostG.root", true);
	 //dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_RunPreG.root", true);

//         trackerCorrector.LoadDeDxCalibration  (DIRNAME+"/../../../data/Data13TeVGains_v2.root");
	         trackerCorrector.TrackerGains = NULL;
   }else{
         dEdxSF [0]      = 1.09711;
         dEdxSF [1]      = 1.09256;
//	 dEdxSF_OldCC[0] = 1.01705;
//	 dEdxSF_OldCC[1] = 1.09522;
//	 dEdxSF_NewCC[0] = 1.04345;
//	 dEdxSF_NewCC[1] = 1.16599;
//         dEdxTemplates    = loadDeDxTemplate(DIRNAME + "/../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", true);
//         dEdxTemplatesInc = loadDeDxTemplate(DIRNAME + "/../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd.root", false); 
//         dEdxTemplates      = loadDeDxTemplate (DIRNAME+"/Templates/dEdxTemplate_hit_SP_MCMinBias.root"           , true);
//         dEdxTemplatesIn    = loadDeDxTemplate (DIRNAME+"/Templates/dEdxTemplate_hit_SP_in_noC_MCMinBias.root"    , true);
//         dEdxTemplatesInc   = loadDeDxTemplate (DIRNAME+"/Templates/dEdxTemplate_hit_SP_MCMinBias.root"           , false);
//         dEdxTemplatesCC    = loadDeDxTemplate (DIRNAME+"/../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CC.root" , true);
//         dEdxTemplatesCI    = loadDeDxTemplate (DIRNAME+"/../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CI.root" , true);
         dEdxTemplatesCCC   = loadDeDxTemplate (DIRNAME+"/../../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true);
         dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/MC13TeV16_dEdxTemplate.root", true);
//         dEdxTemplatesnewCCC= loadDeDxTemplate (DIRNAME+"/Templates/dEdxTemplate_hit_SP_in_noC_newCCC_MCMinBias.root", true);
         trackerCorrector.TrackerGains = NULL;
   }

   if (!trackerCorrector.TrackerGains) std::cerr << "Could not load TrackerGains!!!" << std::endl;
   if (trackerCorrector.TrackerGains) std::cerr << "TrackerGains loaded!!!" << std::endl;

   ifstream CorrFact2017L1,CorrFact2017L2,CorrFact2017L3,CorrFact2017L4,CorrFact2017R1,CorrFact2017R2;
   CorrFact2017L1.open("./CorrFactHistory2017L1.txt");
   CorrFact2017L2.open("./CorrFactHistory2017L2.txt");
   CorrFact2017L3.open("./CorrFactHistory2017L3.txt");
   CorrFact2017L4.open("./CorrFactHistory2017L4.txt");
   CorrFact2017R1.open("./CorrFactHistory2017R1.txt");
   CorrFact2017R2.open("./CorrFactHistory2017R2.txt");
   

   TFile* OutputHisto = new TFile((OUTPUT).c_str(),"RECREATE");  //File must be opened before the histogram are created

 

   std::vector<dEdxStudyObj*> results;
//   results.push_back(new dEdxStudyObj("hit_PO"      , 0, 1, NULL, 2.7, 3.2, NULL, false, false, false, 1, 1) );
//   results.push_back(new dEdxStudyObj("hit_PO_noHIP", 0, 1, NULL, 2.7, 3.2, NULL, false, false, false, 1, 1, 0, false, false) );
//   results.push_back(new dEdxStudyObj("hit_SO_raw"  , 0, 2, NULL, 2.7, 3.2, NULL) );
//   results.push_back(new dEdxStudyObj("hit_SO_in_noC_CCC", 0, 2, NULL, 2.7, 3.2  , trackerCorrector.TrackerGains, true, true, false, 1, 1) );
//   results.push_back(new dEdxStudyObj("hit_SO_in_noC_CCC_noHIP", 0, 2, NULL, 2.7, 3.2  , trackerCorrector.TrackerGains, true, true, false, 1, 1, 0, false, false) );
//   results.push_back(new dEdxStudyObj("hit_SO_in_noC_newCCC", 0, 2, NULL, 2.7, 3.2  , trackerCorrector.TrackerGains, true, true, false, 2, 2) );
//   results.push_back(new dEdxStudyObj("hit_SO_in_noC_newCCC_noHIP", 0, 2, NULL, 2.7, 3.2  , trackerCorrector.TrackerGains, true, true, false, 2, 2, 0, false, false) );
//   results.push_back(new dEdxStudyObj("hit_SP"      , 0, 3, NULL, 2.7, 3.2  , trackerCorrector.TrackerGains) );
//   results.push_back(new dEdxStudyObj("hit_SO_in"   , 0, 2, NULL, 2.7, 3.2  , trackerCorrector.TrackerGains, true) );
//   results.push_back(new dEdxStudyObj("hit_SP_in_noC", 0, 3, NULL, 2.7, 3.2 , trackerCorrector.TrackerGains, true) );
//   results.push_back(new dEdxStudyObj("hit_SP_in_noC_CI" , 0, 3, NULL,2.7,3.2, trackerCorrector.TrackerGains, true, true, false, false, 1) );
//   results.push_back(new dEdxStudyObj("hit_SP_in_noC_CC" , 0, 3, NULL,2.7,3.2, trackerCorrector.TrackerGains, true, true, false, true,  0) );
//   results.push_back(new dEdxStudyObj("hit_SP_in_noC_CCC", 0, 3, NULL, 2.7, 3.2, trackerCorrector.TrackerGains, true, true, false, 1,  1) );
//   results.push_back(new dEdxStudyObj("hit_SP_in_noC_newCCC", 0, 3, NULL, 2.7, 3.2, trackerCorrector.TrackerGains, true, true, false, 2,  2) );
//////   results.push_back(new dEdxStudyObj("harm2_PO_raw", 1, 1, NULL, 2.7, 3.2 , NULL) );
//////   results.push_back(new dEdxStudyObj("Hybr2005_PO_raw", 1, 1, NULL, 2.7, 3.2 , NULL, false, false, false, true, 1, 0.05, true) );
//////   results.push_back(new dEdxStudyObj("Hybr2010_PO_raw", 1, 1, NULL, 2.7, 3.2 , NULL, false, false, false, true, 1, 0.10, true) );
//////   results.push_back(new dEdxStudyObj("Hybr2015_PO_raw", 1, 1, NULL, 2.7, 3.2 , NULL, false, false, false, true, 1, 0.15, true) );
//   results.push_back(new dEdxStudyObj("harm2_SO"    , 1, 2, NULL            , trackerCorrector.TrackerGains) );
//   results.push_back(new dEdxStudyObj("harm2_SO_FS" , 1, 2, NULL            , trackerCorrector.TrackerGains, false, false, true) );
//   results.push_back(new dEdxStudyObj("harm2_SO_in" , 1, 2, NULL            , trackerCorrector.TrackerGains, true) );
//   results.push_back(new dEdxStudyObj("harm2_SO_in_noC"       , 1, 2, NULL  , trackerCorrector.TrackerGains, true, true) );
//   results.push_back(new dEdxStudyObj("harm2_SO_in_noC_CI"    , 1, 2, NULL  , trackerCorrector.TrackerGains, true, true, false, false, 1) );
//   results.push_back(new dEdxStudyObj("harm2_SO_in_noC_CC"    , 1, 2, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  0) );
   results.push_back(new dEdxStudyObj("harm2_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.600:2.655,isData?3.249:3.371, trackerCorrector.TrackerGains, true, true, false, true,  1) );
//////   results.push_back(new dEdxStudyObj("Hybr2005_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.687:2.938, isData?2.600:2.655, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.05, true) );
//////   results.push_back(new dEdxStudyObj("Hybr2010_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.687:2.938, isData?2.600:2.655, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.10, true) );
/////   results.push_back(new dEdxStudyObj("Hybr2015_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.687:2.938, isData?2.600:2.655, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.15, true) );
//////////   results.push_back(new dEdxStudyObj("Hybr2010_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.399:2.812, isData?3.496:3.178, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.1, true) ); //harm2 with 0.10 low value drop
   results.push_back(new dEdxStudyObj("Hybr2015_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.399:2.812, isData?3.496:3.178, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.15, true) ); //harm2 with 0.15 low value drop
////////   results.push_back(new dEdxStudyObj("Hybr2020_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.399:2.812, isData?3.496:3.178, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.2, true) ); //harm2 with 0.20 low value drop
////////   results.push_back(new dEdxStudyObj("Hybr2025_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.399:2.812, isData?3.496:3.178, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.25, true) ); //harm2 with 0.25 low value drop
////////   results.push_back(new dEdxStudyObj("Hybr2030_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.399:2.812, isData?3.496:3.178, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.30, true) ); //harm2 with 0.30 low value drop
////////   results.push_back(new dEdxStudyObj("Hybr2035_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.399:2.812, isData?3.496:3.178, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.35, true) ); //harm2 with 0.35 low value drop
/////////   results.push_back(new dEdxStudyObj("Hybr2040_SO_in_noC_CCC"   , 1, 2, NULL, isData?2.399:2.812, isData?3.496:3.178, trackerCorrector.TrackerGains, true, true, false, true,  1, 0.40, true) ); //harm2 with 0.40 low value drop
//   results.push_back(new dEdxStudyObj("harm2_SP"    , 1, 3, NULL            , trackerCorrector.TrackerGains) );
//   results.push_back(new dEdxStudyObj("harm2_SP_in" , 1, 3, NULL            , trackerCorrector.TrackerGains, true) );
//   results.push_back(new dEdxStudyObj("harm2_SP_in_noC"       , 1, 3, NULL  , trackerCorrector.TrackerGains, true, true) );
//   results.push_back(new dEdxStudyObj("harm2_SP_in_noC_CI"    , 1, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, false, 1) );
//   results.push_back(new dEdxStudyObj("harm2_SP_in_noC_CC"    , 1, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  0) );
////////   results.push_back(new dEdxStudyObj("harm2_SP_in_noC_CCC"     , 1, 3, NULL, isData?2.447:2.839, isData?3.218:3.029  , trackerCorrector.TrackerGains, true, true, false, 1,  1) );
//   results.push_back(new dEdxStudyObj("harm2_SP_in_noC_newCCC"     , 1, 3, NULL, isData?2.355:2.759, isData?3.282:3.074  , trackerCorrector.TrackerGains, true, true, false, 2,  2) );
//   results.push_back(new dEdxStudyObj("hybr201_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.389:2.804, isData?3.413:3.134  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.1) );
  results.push_back(new dEdxStudyObj("hybr2015_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.399:2.812, isData?3.496:3.178  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.15) );
  results.push_back(new dEdxStudyObj("hybr2015_SP_in_noC_CCC_noL1Pxl"   , 1, 3, NULL, isData?2.399:2.812, isData?3.496:3.178  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.15, false, true, true) );
//   results.push_back(new dEdxStudyObj("hybr202_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.408:2.794, isData?3.534:3.288  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.2) );
//   results.push_back(new dEdxStudyObj("hybr2025_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.361:2.760, isData?3.641:3.391  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.25) );
//   results.push_back(new dEdxStudyObj("hybr203_SP_in_noC_CCC"   , 1, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.3) );
//   results.push_back(new dEdxStudyObj("hybr2035_SP_in_noC_CCC"   , 1, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.35) );
//   results.push_back(new dEdxStudyObj("hybr204_SP_in_noC_CCC"   , 1, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.4) );
//   results.push_back(new dEdxStudyObj("Hybr201_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.525:2.883, isData?3.389:3.114  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.1, true) );
//   results.push_back(new dEdxStudyObj("Hybr2015_SP_in_noC_CCC_"   , 1, 3, NULL, isData?2.580:2.935, isData?3.922:3.197  , trackerCorrector.TrackerGains, true, true, false, 1,  1, 0.15, true) );     //THIS IS THE TO KEEP
//   results.push_back(new dEdxStudyObj("Hybr2020_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.580:2.935, isData?3.922:3.197  , trackerCorrector.TrackerGains, true, true, false, 1,  1, 0.20, true) );
//   results.push_back(new dEdxStudyObj("Hybr2025_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.580:2.935, isData?3.922:3.197  , trackerCorrector.TrackerGains, true, true, false, 1,  1, 0.25, true) );
//   results.push_back(new dEdxStudyObj("Hybr2030_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.580:2.935, isData?3.922:3.197  , trackerCorrector.TrackerGains, true, true, false, 1,  1, 0.30, true) );
//   results.push_back(new dEdxStudyObj("Hybr2035_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.580:2.935, isData?3.922:3.197  , trackerCorrector.TrackerGains, true, true, false, 1,  1, 0.35, true) );
//   results.push_back(new dEdxStudyObj("Hybr2015_SP_in_noC_newCCC"   , 1, 3, NULL, isData?2.567:2.863, isData?3.426:3.268  , trackerCorrector.TrackerGains, true, true, false, 2,  2, 0.15, true) );
//   results.push_back(new dEdxStudyObj("Hybr202_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.739:2.969, isData?3.472:3.305  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.2, true) );
//   results.push_back(new dEdxStudyObj("Hybr2025_SP_in_noC_CCC"   , 1, 3, NULL, isData?2.869:3.075, isData?3.486:3.333  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.25, true) );
//   results.push_back(new dEdxStudyObj("harm2_SP_in_noC_CCC_noF"     , 1, 3, NULL, isData?2.448:2.866, isData?3.218:2.963  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.0, true, false) );
//   results.push_back(new dEdxStudyObj("Hybr201_SP_in_noC_CCC_noF"   , 1, 3, NULL, isData?2.525:2.883, isData?3.389:3.114  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.1, true, false) );
//   results.push_back(new dEdxStudyObj("Hybr2015_SP_in_noC_CCC_noF"   , 1, 3, NULL, isData?2.687:2.938, isData?3.383:3.199  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.15, true, false) );
//   results.push_back(new dEdxStudyObj("Hybr202_SP_in_noC_CCC_noF"   , 1, 3, NULL, isData?2.739:2.969, isData?3.472:3.305  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.2, true, false) );
//   results.push_back(new dEdxStudyObj("Hybr2025_SP_in_noC_CCC_noF"   , 1, 3, NULL, isData?2.869:3.075, isData?3.486:3.333  , trackerCorrector.TrackerGains, true, true, false, true,  1, 0.25, true, false) );
   /*results.push_back(new dEdxStudyObj("trunc40_PO_raw", 3, 1, NULL            , NULL) );
   results.push_back(new dEdxStudyObj("trunc40_SO"    , 3, 2, NULL            , trackerCorrector.TrackerGains) );
   results.push_back(new dEdxStudyObj("trunc40_SO_FS" , 3, 2, NULL            , trackerCorrector.TrackerGains, false, false, true) );
   results.push_back(new dEdxStudyObj("trunc40_SO_in" , 3, 2, NULL            , trackerCorrector.TrackerGains, true) );
   results.push_back(new dEdxStudyObj("trunc40_SO_in_noC"       , 3, 2, NULL  , trackerCorrector.TrackerGains, true, true) );
   results.push_back(new dEdxStudyObj("trunc40_SO_in_noC_CI"    , 3, 2, NULL  , trackerCorrector.TrackerGains, true, true, false, false, 1) );
   results.push_back(new dEdxStudyObj("trunc40_SO_in_noC_CC"    , 3, 2, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  0) );
   results.push_back(new dEdxStudyObj("trunc40_SO_in_noC_CCC"   , 3, 2, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  1) );
   results.push_back(new dEdxStudyObj("trunc40_SP"    , 3, 3, NULL            , trackerCorrector.TrackerGains) );
   results.push_back(new dEdxStudyObj("trunc40_SP_in" , 3, 3, NULL            , trackerCorrector.TrackerGains, true) );
   results.push_back(new dEdxStudyObj("trunc40_SP_in_noC"       , 3, 3, NULL  , trackerCorrector.TrackerGains, true, true) );
   results.push_back(new dEdxStudyObj("trunc40_SP_in_noC_CI"    , 3, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, false, 1) );
   results.push_back(new dEdxStudyObj("trunc40_SP_in_noC_CC"    , 3, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  0) );*/
//   results.push_back(new dEdxStudyObj("trunc40_SP_in_noC_CCC"   , 3, 3, NULL  , trackerCorrector.TrackerGains, true, true, false, true,  1) );

//   results.push_back(new dEdxStudyObj("Ias_PO"      , 2, 1, dEdxTemplates, 2.7, 3.2 , NULL) );
//   results.push_back(new dEdxStudyObj("Ias_SO_inc"  , 2, 2, dEdxTemplatesInc, NULL) );
//   results.push_back(new dEdxStudyObj("Ias_SO"      , 2, 2, dEdxTemplates   , NULL) );
//   results.push_back(new dEdxStudyObj("Ias_SO_in"   , 2, 2, dEdxTemplatesIn , NULL, true) );
//   results.push_back(new dEdxStudyObj("Ias_SO_in_noC"         , 2, 2, dEdxTemplatesIn , NULL, true, true) );
//   results.push_back(new dEdxStudyObj("Ias_SO_in_noC_CI"      , 2, 2, dEdxTemplatesCI , NULL, true, true, false, false, 1) );
//   results.push_back(new dEdxStudyObj("Ias_SO_in_noC_CC"      , 2, 2, dEdxTemplatesCC , NULL, true, true, false, true,  0) );
//   results.push_back(new dEdxStudyObj("Ias_SO_in_noC_CCC"     , 2, 2, dEdxTemplatesCCC, 2.7, 3.2, NULL, true, true, false, true,  1) );
//   results.push_back(new dEdxStudyObj("Ias_SP_inc"  , 2, 3, dEdxTemplatesInc, NULL) );
//   results.push_back(new dEdxStudyObj("Ias_SP"      , 2, 3, dEdxTemplates   , NULL) );
//   results.push_back(new dEdxStudyObj("Ias_SP_in"   , 2, 3, dEdxTemplatesIn , NULL, true) );
//   results.push_back(new dEdxStudyObj("Ias_SP_in_noC"         , 2, 3, dEdxTemplatesIn , NULL, true, true) );
//   results.push_back(new dEdxStudyObj("Ias_SP_in_noC_CI"      , 2, 3, dEdxTemplatesCI , NULL, true, true, false, false, 1) );
//   results.push_back(new dEdxStudyObj("Ias_SP_in_noC_CC"      , 2, 3, dEdxTemplatesCC , NULL, true, true, false, true,  0) );
//   results.push_back(new dEdxStudyObj("Ias_SP_in_noC_CCC"     , 2, 3, dEdxTemplatesCCC,   2.7, 3.2, NULL, true, true, false, 1,  1) );
/////////results.push_back(new dEdxStudyObj("Ias_SP_in_noC_CCC16"   , 2, 3, dEdxTemplatesCCC16, 2.7, 3.2, NULL, true, true, false, 1,  1) );
results.push_back(new dEdxStudyObj("Ias_SO_in_noC_CCC16"   , 2, 2, dEdxTemplatesCCC16, 2.7, 3.2, NULL, true, true, false, 1,  1) );
//  results.push_back(new dEdxStudyObj("Ias_PO_in_noC_CCC16"   , 2, 1, dEdxTemplatesCCC16, 2.7, 3.2, NULL, true, true, false, 1,  1) );
//   results.push_back(new dEdxStudyObj("Ias_SP_in_noC_newCCC"     , 2, 3, dEdxTemplatesnewCCC, 2.7, 3.2, NULL, true, true, false, 2,  2) );


   TH1F* h_event = new TH1F("h_event","nof events",1,0,1);
   TH1F* h_goodVertices = new TH1F("h_goodVertices","nof good vertices",100,0,100);
   TH1F* h_instLumi = new TH1F("h_instLumi","instantaneous luminosity",10000,0,20000);
   TH1F* h_bunchLumi = new TH1F("h_bunchLumi","bunchLumi",10000,0,10000);
   TH1F* h_pileup = new TH1F("h_pileup","pile up",100,0,100);
   TH2D* InstantLumiVsGoodVertices  = new TH2D("instLumiVsGoodVtx", ";instantaneous luminosity;good vertices",10000,0,20000,100,0,100);
   TH2D* InstantLumiVsBunchLumi = new TH2D("instLumiVsBunchLumi",";instantaneous luminosity;bunchLumi",10000,0,20000,10000,0,10000);
   TH2D* PileUpVsInstantLumi = new TH2D("puVsInstLumi",";pile up;instantaneous luminosity",100,0,100,10000,0,20000);
   TH2D* PileUpVsBunchLumi = new TH2D("puVsBunchLumi",";pile up;bunchLumi",100,0,100,10000,0,10000);
   TH1F* h_trackPhi = new TH1F("h_trackPhi",";phi;",Phi_NBins,Phi_Min,Phi_Max);
   
   TH1F* h_invariantMass_prescaleTrack = new TH1F("invMass_prescaleTrack",";m [GeV];",1000,0,10);
   TH1F* h_invariantMass_prescaleMulti = new TH1F("invMass_prescaleMulti",";m [GeV];",1000,0,10);

   TH2F* h2_pions = new TH2F("pion_line",";p [GeV];dE/dx [MeV/cm]",300,0,60,100,0,20);
   TH2F* h2_pions_ih = new TH2F("pion_line_ih",";p [GeV];Ih [MeV/cm]",300,0,60,100,0,20);
   TH2F* h2_kaons = new TH2F("kaon_line",";p [GeV];dE/dx [MeV/cm]",300,0,60,100,0,20);
   TH2F* h2_kaons_ih = new TH2F("kaon_line_ih",";p [GeV];Ih [MeV/cm]",300,0,60,100,0,20);
   TH2F* h2_muons = new TH2F("muon_line",";p [GeV];dE/dx [MeV/cm]",5000,0,1000,100,0,20);
   TH2F* h2_muons_ih = new TH2F("muon_line_ih",";p [GeV];Ih [MeV/cm]",5000,0,1000,100,0,20);

   TH2F* h2_muons_Eta0p83 = new TH2F("muon_line_Eta0p83","|#eta|<0.83;p [GeV];dE/dx [MeV/cm]",5000,0,1000,100,0,20);
   TH2F* h2_muons_0p83Eta1p44 = new TH2F("muon_line_0p83Eta1p44","0.83<=|#eta|<1.44;p [GeV];dE/dx [MeV/cm]",5000,0,1000,100,0,20);
   TH2F* h2_muons_1p44Eta2p1 = new TH2F("muon_line_1p44Eta2p1","1.44<=|#eta|<2.1;p [GeV];dE/dx [MeV/cm]",5000,0,1000,100,0,20);
   TH2F* h2_muons_2p1Eta = new TH2F("muon_line_2p1Eta","2.1<=|#eta|;p [GeV];dE/dx [MeV/cm]",5000,0,1000,100,0,20);

   TH2F* HdedxVsEta_muonPhigher70  = new TH2F("dedxVsEta_muonPhigher70 ","muon momentum > 70 GeV;#eta;dE/dx [MeV/cm]", 1000,-5,  5, 100,0, 20);

   TH3F* EtaPhiVsP_dEdx0p5          = new TH3F("EtaPhiVsP_dEdx0p5", ";eta;phi;momentum",   Eta_NBins,Eta_Min,  Eta_Max, Phi_NBins, Phi_Min, Phi_Max,500,0,10);
   TH2F* EtaPhi_dEdx0p5_Phigher5    = new TH2F("EtaPhi_dEdx0p5_Phigher5","p>5GeV;eta;phi",   Eta_NBins,Eta_Min,  Eta_Max, Phi_NBins, Phi_Min, Phi_Max);
   TH2F* EtaPhi_dEdx0p5             = new TH2F("EtaPhi_dEdx0p5",";eta;phi",   Eta_NBins,Eta_Min,  Eta_Max, Phi_NBins ,Phi_Min, Phi_Max);        
         


   h_trackPhi->Sumw2();

   h2_pions->Sumw2();
   h2_pions_ih->Sumw2();
   h2_kaons->Sumw2();
   h2_kaons_ih->Sumw2();
   h2_muons->Sumw2();
   h2_muons_ih->Sumw2();
  
   h2_muons_Eta0p83->Sumw2();
   h2_muons_0p83Eta1p44->Sumw2();
   h2_muons_1p44Eta2p1->Sumw2();
   h2_muons_2p1Eta->Sumw2();

  cout<<"results size: "<<results.size()<<endl;
  for (size_t r=0; r<results.size(); r++) cout<<"results["<<r<<"] name: "<<results[r]->Name<<endl;
   printf("Progressing Bar           :0%%       20%%       40%%       60%%       80%%       100%%\n");
   long Ncluster_lowP=0;
   long Nsatcluster_lowP=0;
   for(unsigned int f=0;f<FileName.size();f++){
     cout<<"opening file: "<<FileName[f].c_str()<<endl;
     TFile* file = TFile::Open(FileName[f].c_str() );
     if(file==NULL) continue;
     cout<<"file loaded"<<endl;
     fwlite::Event ev(file);
     cout<<"events loaded"<<endl;
     cout<<"ev size: "<<ev.size()<<endl;
     printf("Scanning the ntuple %2i/%2i :", (int)f+1, (int)FileName.size());
     long treeStep(ev.size()/50);
     long iev=0;
     for(ev.toBegin(); !ev.atEnd(); ++ev){ iev++;
         if(iev%treeStep==0){printf(".");fflush(stdout);}
         h_event->Fill(iev);
         //cout<<"in the event loop!"<<endl;
	 HIPemulator.setEventRate();
     //if(iev>100) continue;

         if(CurrentRun != ev.eventAuxiliary().run()){
            CurrentRun = ev.eventAuxiliary().run();
            trackerCorrector.setRun (CurrentRun);

	    /* if (isData){
	          switch (CurrentRun)
	      case 278018: dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_Run278018.root", true); dEdxSF [0] = 1.00000; dEdxSF [1] = 1.41822*1.04098; break;
	      case 278308: dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_Run278308.root", true); dEdxSF [0] = 1.00000; dEdxSF [1] = 1.41822*1.06009;break;
	      case 279931: dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_Run279931.root", true); dEdxSF [0] = 1.00000; dEdxSF [1] = 1.41822*1.12230;break;
	      default:  dEdxTemplatesCCC16 = loadDeDxTemplate (DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_Run278018.root", true); dEdxSF [0] = 1.00000; dEdxSF [1] = 1.41822;break;
              
	     for (size_t r=0; r<results.size(); r++)
	       {
		 if (results[r]->Name.find("Ias") != std::string::npos)
		   results[r]->dEdxTemplates = dEdxTemplatesCCC16;
	       }
	     }*/
         }

         float lumi = -1;
         float bunchLumi = -1;
         float pileup = 0;

         fwlite::Handle<LumiScalersCollection> lumiScalers;
         lumiScalers.getByLabel(ev, "scalersRawToDigi");
         if(lumiScalers.isValid() && !lumiScalers->empty())
         {
             LumiScalersCollection::const_iterator scalit = lumiScalers->begin();
             lumi = scalit->instantLumi();
             bunchLumi = scalit->bunchLumi();
             pileup = scalit->pileup();
         }

         h_bunchLumi->Fill(bunchLumi);
         h_instLumi->Fill(lumi);
         h_pileup->Fill(pileup);
         PileUpVsInstantLumi->Fill(pileup,lumi);
         PileUpVsBunchLumi->Fill(pileup,bunchLumi);
         InstantLumiVsBunchLumi->Fill(lumi,bunchLumi);

         fwlite::Handle<DeDxHitInfoAss> dedxCollH;
         dedxCollH.getByLabel(ev, "dedxHitInfo");
         if(!dedxCollH.isValid()){printf("Invalid dedxCollH\n");continue;}

         fwlite::Handle<edm::ValueMap<int>> dedxHitInfoPrescale;
         dedxHitInfoPrescale.getByLabel(ev, "dedxHitInfo", "prescale");
         if(!dedxHitInfoPrescale.isValid()){printf("Invalid dedxHitInfoPrescale\n");continue;}

         fwlite::Handle< std::vector<reco::Track> > trackCollHandle;
         trackCollHandle.getByLabel(ev,"RefitterForDeDx");
         if(!trackCollHandle.isValid()){
            trackCollHandle.getByLabel(ev,"generalTracks");
               if (!trackCollHandle.isValid()){
                  printf("Invalid trackCollHandle\n");
                  continue;
               }
         }

         fwlite::Handle< std::vector<reco::Muon> > muonCollHandle;
         muonCollHandle.getByLabel(ev, "muons");
         if(!muonCollHandle.isValid()){printf("Invalid muonCollHandle\n");continue;}

         fwlite::Handle < std::vector<reco::Vertex> > vertexCollHandle;
         vertexCollHandle.getByLabel(ev, "offlinePrimaryVertices");
         if(!vertexCollHandle.isValid()){printf("Vertex Collection not found!\n"); continue;}
         const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
         if(vertexColl.size()<1){printf("NO VERTICES\n"); continue;}

         int goodVertices=0;
         for(unsigned int i=0;i<vertexColl.size();i++){
             if(vertexColl[i].isFake() || fabs(vertexColl[i].z())>24 || vertexColl[i].position().rho()>2 || vertexColl[i].ndof()<=4) continue;
             goodVertices++;
         }
         h_goodVertices->Fill(goodVertices);

         InstantLumiVsGoodVertices->Fill(lumi,goodVertices);

         fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
         if(isSignal){
            //get the collection of generated Particles
            genCollHandle.getByLabel(ev, "genParticlesSkimmed");
            if(!genCollHandle.isValid()){
               genCollHandle.getByLabel(ev, "genParticles");
               if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound\n");continue;}
            }
         }

         //cout<<"all collections loaded"<<endl;

         for(unsigned int m=0;m<muonCollHandle->size();m++){
                reco::MuonRef muon = reco::MuonRef(muonCollHandle.product(),m);
                if(muon.isNull()) continue;
                if(!muon->isGlobalMuon()) continue;
                if(muon->track().isNull()) continue;
                const DeDxHitInfo* dedxHits = NULL;
                DeDxHitInfoRef dedxHitsRef = dedxCollH->get(muon->track().key());
                if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
                if(!dedxHits)continue;
                int prescale = 1;
                if(dedxHitsRef.isNonnull()){
                    prescale = (*dedxHitInfoPrescale)[dedxHitsRef];
                }
                for(unsigned int h=0;h<dedxHits->size();h++){
                    DetId detid(dedxHits->detId(h));
                    double SF = dEdxSF[0];
                    if(detid.subdetId()<3) continue;
                    double Norm = 3.61e-06*265;
                    h2_muons->Fill(muon->track()->p(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);
                    if(fabs(muon->track()->eta()) < 0.83)  h2_muons_Eta0p83->Fill(muon->track()->p(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);
                    if(0.83 <= fabs(muon->track()->eta()) && fabs(muon->track()->eta()) < 1.44) h2_muons_0p83Eta1p44->Fill(muon->track()->p(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);
                    if(1.44 <= fabs(muon->track()->eta()) && fabs(muon->track()->eta()) < 2.1) h2_muons_1p44Eta2p1->Fill(muon->track()->p(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);
                    if(2.1 <= fabs(muon->track()->eta())) h2_muons_2p1Eta->Fill(muon->track()->p(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);    
                    if(muon->track()->p()>=70) HdedxVsEta_muonPhigher70->Fill(muon->track()->eta(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);
                }
                DeDxData dedxObjMuon   = computedEdx(dedxHits, dEdxSF, 0, false, true, false, false, 0, true);
                h2_muons_ih->Fill(muon->track()->p(),dedxObjMuon.dEdx(),prescale);
        }



         for(unsigned int c=0;c<trackCollHandle->size();c++){
            //basic track quality cuts
            reco::TrackRef track = reco::TrackRef( trackCollHandle.product(), c );
            if(track.isNull())continue;
            if(track->chi2()/track->ndof()>5 )continue;
            if(track->found()<8)continue;

            if(isSignal){
               if(track->pt()<45)continue;
               const std::vector<reco::GenParticle>& genColl = *genCollHandle;
               if (DistToHSCP (track, genColl)>0.03) continue;
            }else{
               if(track->pt()>=45)continue;
            }

            //load dEdx informations
            const DeDxHitInfo* dedxHits = NULL;
            DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
            if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
            if(!dedxHits)continue;

            //load prescale for 2017 - inspired from /opt/sbg/cms/safe1/cms/dapparu/HSCP/CMSSW_10_6_2/src/stage/ntuple/plugins/ntuple.cc and from https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc
            int prescale = 1;
            if(dedxHitsRef.isNonnull()){
                prescale = (*dedxHitInfoPrescale)[dedxHitsRef];
            }
            
            bool isPion=false;
            bool isKaon=false;
            for(unsigned int d=c+1;d<trackCollHandle->size();d++){
                reco::TrackRef track2 = reco::TrackRef(trackCollHandle.product(),d);
                if(track2.isNull()) continue;
                if(track2->charge()==track->charge()) continue;
                if(track2->chi2()/track2->ndof()>5 )continue;
                if(track2->found()<8)continue;
                if(abs(track->dxy())<0.1) continue;
                if(abs(track->dz())<0.1) continue;
                if(abs(track2->dxy()-track->dxy())>0.002) continue;
                if(abs(track2->dz()-track->dz())>0.005) continue;
                //if(track->p()>10) continue;
                //if(track2->p()>10) continue;
                const DeDxHitInfo* dedxHits2 = NULL;
                DeDxHitInfoRef dedxHitsRef2 = dedxCollH->get(track2.key());
                if(!dedxHitsRef2.isNull())dedxHits2 = &(*dedxHitsRef2);
                if(!dedxHits2)continue;
                int prescale2 = 1;
                if(dedxHitsRef2.isNonnull()){
                    prescale2 = (*dedxHitInfoPrescale)[dedxHitsRef2];
                }
                //prescale/=100; //because track->p < 10
                //prescale2/=100; //beacause track2->p < 10
                int multiPrescale = prescale*prescale2;
                TLorentzVector TL1;
                TLorentzVector TL2;
                TL1.SetPtEtaPhiM(track->pt(),track->eta(),track->phi(),0.);
                TL2.SetPtEtaPhiM(track2->pt(),track2->eta(),track2->phi(),0.);
                double InvMass = (TL1+TL2).M();
                h_invariantMass_prescaleTrack->Fill(InvMass,prescale);
                h_invariantMass_prescaleMulti->Fill(InvMass,multiPrescale);
                if((InvMass>0.495 && InvMass<0.499) || (InvMass>0.774 && InvMass<0.776 )) isPion=true; //kaon mass = 497.6 MeV ; rho mass = 775.5 MeV
                if((InvMass>1018 && InvMass<1020)) isKaon=true; //kaon mass = 497.6 MeV ; rho mass = 775.5 MeV
                //std::cout << "invMass: " << InvMass << " isPion: " << isPion << std::endl;
            }
            for(unsigned int h=0;h<dedxHits->size();h++){
                DetId detid(dedxHits->detId(h));
                double SF = dEdxSF[0];
                if(detid.subdetId()<3) continue;
                double Norm = 3.61e-06*265;
                if(isPion) h2_pions->Fill(track->p(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);
                if(isKaon) h2_kaons->Fill(track->p(),SF*Norm*dedxHits->charge(h)/dedxHits->pathlength(h),prescale);
            }
            DeDxData dedxObjPion   = computedEdx(dedxHits, dEdxSF, 0, false, true, false, false, 0, true);
//computedEdx( dedxHits, scaleFactors, templateHisto, usePixel, useClusterCleaning, reverseProb, useTruncated, TrackerGains, useStrip, mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo, dropLowerDeDxValue, dedxHIPEmulator* hipEmulator, double* dEdxErr){
            if(isPion) h2_pions_ih->Fill(track->p(),dedxObjPion.dEdx(),prescale);
            if(isKaon) h2_kaons_ih->Fill(track->p(),dedxObjPion.dEdx(),prescale);

            if(track->p() <= 1){
                for(unsigned int h=0;h<dedxHits->size();h++){
                    bool sat=false;
                    DetId detid(dedxHits->detId(h));
                    if(detid.subdetId()<3) continue;
                    vector <int> amps = convert(dedxHits->stripCluster(h)->amplitudes());
                    for(int j=0;j<amps.size();j++){
                        if(amps[j]>253) sat=true;
                    }
                    Ncluster_lowP++;
                    if(sat) Nsatcluster_lowP++;
                }
            }

            //hit level dEdx information (only done for MIPs)
            //if(track->pt() > 5) //FIXME why this pt cut ?
            if(track->pt() > 0){ //FIXME why this pt cut ?
               for(unsigned int h=0;h<dedxHits->size();h++){
                   DetId detid(dedxHits->detId(h));
                   int moduleGeometry = 0; // underflow bin -- debug purposes
                   int layer = 0;
                   if(detid.subdetId()>=3){ SiStripDetId SSdetId(detid); moduleGeometry = SSdetId.moduleGeometry(); if (moduleGeometry==15) {cerr << "ERROR! There is no SiStrip geometry 15!" << endl; exit (EXIT_FAILURE);}}
                   else if(detid.subdetId()<3){moduleGeometry = 15;} // 15 is for pixel
                   if(detid.subdetId()==3){layer=((detid>>14)&0x7);} //TIB
                   if(detid.subdetId()==4){layer=((detid>>9)&0x3)+10;} //TID
                   if(detid.subdetId()==5){layer=((detid>>14)&0x7)+4;} //TOB
                   if(detid.subdetId()==6){layer=((detid>>5)&0x7)+13;} //TEC
                   for(unsigned int R=0;R<results.size();R++){
//                      if (results[R]->Name.find("newCCC")!=string::npos){dEdxSF[0] = dEdxSF_NewCC[0]; dEdxSF[1] = dEdxSF_NewCC[1];}
//                      else {dEdxSF[0] = dEdxSF_OldCC[0]; dEdxSF[1] = dEdxSF_OldCC[1];}
                      double scaleFactor = dEdxSF[0];
                      if (detid.subdetId()<3) scaleFactor *= dEdxSF[1];
                      double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;

                    
                      results[R]->MG_PathlengthVsEta->Fill(moduleGeometry,dedxHits->pathlength(h),track->eta(),prescale);
                      results[R]->L_PathlengthVsEta->Fill(layer,dedxHits->pathlength(h),track->eta(),prescale);

                      if(results[R]->isEstim || results[R]->isDiscrim){
                      int charge = dedxHits->charge(h);
                      if (detid.subdetId()>=3 && results[R]->crossTalkInvAlgo==1){ //in case of crossTalkInv, give the corrected cluster charge
                         vector <int> amps = convert(dedxHits->stripCluster(h)->amplitudes());
                         amps = CrossTalkInv(amps, 0.10, 0.04, true);
                         charge = std::accumulate(amps.begin(), amps.end(), 0);}

                      results[R]->MG_Charge_Vs_Path->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      results[R]->L_Charge_Vs_Path->Fill (layer, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      
                      if(lumi<5000) results[R]->MG_Charge_Vs_Path_instLumi_0_5000->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      if(lumi>=5000 && lumi<10000) results[R]->MG_Charge_Vs_Path_instLumi_5000_10000->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      if(lumi>=10000 && lumi<15000) results[R]->MG_Charge_Vs_Path_instLumi_10000_15000->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      if(lumi>=15000) results[R]->MG_Charge_Vs_Path_instLumi_15000_inf->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      
                      if(pileup<25) results[R]->MG_Charge_Vs_Path_PU_0_25->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      if(pileup>= 25 && pileup<30) results[R]->MG_Charge_Vs_Path_PU_25_30->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      if(pileup>= 30 && pileup<35) results[R]->MG_Charge_Vs_Path_PU_30_35->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      if(pileup>= 35 && pileup<40) results[R]->MG_Charge_Vs_Path_PU_35_40->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      if(pileup>= 40) results[R]->MG_Charge_Vs_Path_PU_40_inf->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 

                      if(results[R]->usePixel)
                      {
                          int pixelLayer = 0;
                          if(detid.subdetId()==1) //pixel barrel
                          {
                              pixelLayer = ((detid>>16)&0xF); 
                          }
                          if(detid.subdetId()==2) //pixel endcaps
                          {
                              pixelLayer = ((detid>>16)&0xF);
                              moduleGeometry = 16;
                          }
                          if(pixelLayer==1) moduleGeometry = 17;

                      
                      results[R]->MG_Charge_Vs_Path_wPixels->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1)),prescale); 
                      }
                      
                      }

                      
                      if(!results[R]->isHit) continue; //only consider results related to hit info here
                      if(!results[R]->usePixel && detid.subdetId() <3)continue; // skip pixels
                      if(!results[R]->useStrip && detid.subdetId()>=3)continue; // skip strips
                      if(results[R]->mustBeInside && !isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId()>=3?dedxHits->stripCluster(h):NULL) )continue;
//                      if(results[R]->removeCosmics){ if (isCompatibleWithCosmic(track, vertexColl))continue;} //don't consider hits, which belong to cosmic tracks

                      if(results[R]->useClusterCleaning && detid.subdetId()>=3 && !clusterCleaning(dedxHits->stripCluster(h), results[R]->crossTalkInvAlgo)) continue; //if it fails clusterCleaning, skip it!
//                      if(results[R]->CCFunction && detid.subdetId()>=3 && !results[R]->CCFunction(dedxHits->stripCluster(h), results[R]->crossTalkInvAlgo, NULL)) continue; //if it fails clusterCleaning, skip it!

                      int charge = dedxHits->charge(h);
                      if (detid.subdetId()>=3 && results[R]->crossTalkInvAlgo==1){ //in case of crossTalkInv, give the corrected cluster charge
                         vector <int> amps = convert(dedxHits->stripCluster(h)->amplitudes());
                         amps = CrossTalkInv(amps, 0.10, 0.04, true);
                         charge = std::accumulate(amps.begin(), amps.end(), 0);

                         if (results[R]->useStrip && detid.subdetId()>=3){
                            for (unsigned int i=0; i<amps.size(); i++){
                               if (amps[i] > 0){
                                  results[R]->HStripADC   ->Fill (amps[i],prescale);
                                  results[R]->HStripADCvsP->Fill (track->p(), amps[i],prescale);
                               }
                            }
                         }
                      }
                      

                      double ChargeOverPathlength   = scaleFactor*Norm*charge/dedxHits->pathlength(h);
                      double ChargeOverPathlength_U = 1.0*Norm*charge/dedxHits->pathlength(h);
                      if (results[R]->fakeHIP && !(SuppressFakeHIP))
                         ChargeOverPathlength = HIPemulator.fakeHIP (detid.subdetId(), ChargeOverPathlength);

                      results[R]->HHit->Fill(ChargeOverPathlength,prescale);
                      results[R]->HHit_U->Fill(ChargeOverPathlength_U,prescale);
                      if (fabs(track->eta())<0.4){
                         results[R]->HHitProfile->Fill(track->p(), ChargeOverPathlength,prescale);
                         results[R]->HHitProfile_U->Fill(track->p(), ChargeOverPathlength_U,prescale);
                      }
                      if(results[R]->usePixel && results[R]->useStrip){
                          
                          
                         if(detid.subdetId()>=3)results[R]->Charge_Vs_FS[moduleGeometry]->Fill(dedxHits->stripCluster(h)->firstStrip(), charge,prescale); 
                         results[R]->Charge_Vs_XYH[moduleGeometry]->Fill(dedxHits->pos(h).x(), dedxHits->pos(h).y(),prescale); 
                         if(ChargeOverPathlength<1.6)results[R]->Charge_Vs_XYL[moduleGeometry]->Fill(dedxHits->pos(h).x(), dedxHits->pos(h).y(),prescale); 
    
                         if(moduleGeometry>=1 && moduleGeometry<=14){ // FIXME we don't have the geometry information for Pixels yet (TkModGeom* arrays) !!!
                            double nx, ny;
                            if(moduleGeometry<=4){
                               ny = dedxHits->pos(h).y() /  TkModGeomLength[moduleGeometry];
                               nx = dedxHits->pos(h).x() /  TkModGeomWidthT[moduleGeometry];
                            }else{
                               double  offset = TkModGeomLength[moduleGeometry] * (TkModGeomWidthT[moduleGeometry]+TkModGeomWidthB[moduleGeometry]) / (TkModGeomWidthT[moduleGeometry]-TkModGeomWidthB[moduleGeometry]);  // check sign if GeomWidthT[moduleGeometry] < TkModGeomWidthB[moduleGeometry] !!! 
                               double  tan_a = TkModGeomWidthT[moduleGeometry] / std::abs(offset + TkModGeomLength[moduleGeometry]);
                               ny = dedxHits->pos(h).y() /  TkModGeomLength[moduleGeometry];
                               nx = dedxHits->pos(h).x() / (tan_a*std::abs(dedxHits->pos(h).y()+offset));
                            }
                            //printf("%i - %f - %f --> %f - %f\n", moduleGeometry, dedxHits->pos(h).x(), dedxHits->pos(h).y(), nx, ny);
                            results[R]->Charge_Vs_XYHN[moduleGeometry]->Fill(nx, ny,prescale); 
                            if(ChargeOverPathlength<1.6)results[R]->Charge_Vs_XYLN[moduleGeometry]->Fill(nx, ny,prescale);
                         }
                      }
                   }
                }
             }

//             bool isCosmic = isCompatibleWithCosmic(track, vertexColl);
             bool lockOnTrack=false;
             double dEdxDebug = 0;
             for(unsigned int R=0;R<results.size();R++){
                if(!results[R]->isEstim and !results[R]->isDiscrim) continue; //only consider results related to estimator/discriminator variables here
//                if(results[R]->removeCosmics && isCosmic)continue; //don't consider cosmic tracks

//                if (results[R]->Name.find("newCCC")!=string::npos){dEdxSF[0] = dEdxSF_NewCC[0]; dEdxSF[1] = dEdxSF_NewCC[1];}
//                else {dEdxSF[0] = dEdxSF_OldCC[0]; dEdxSF[1] = dEdxSF_OldCC[1];}
                DeDxData dedxObj   = computedEdx(dedxHits, dEdxSF, results[R]->dEdxTemplates, results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 99, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, (results[R]->fakeHIP && !(SuppressFakeHIP))?&HIPemulator:NULL,NULL,results[R]->skipPixelL1);
//                DeDxData dedxObj_U = computedEdx(dedxHits, dEdx_U, results[R]->dEdxTemplates, results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 99, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, results[R]->trimPixel, (results[R]->fakeHIP && !(SuppressFakeHIP))?&HIPEmulator:NULL);

//                DeDxData dedxObj   = computeDeDx (dedxHits, dEdxSF, results[R]->dEdxTemplates, true, results[R]->CCFunction, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 99, false, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, results[R]->trimPixel, (results[R]->fakeHIP && !(SuppressFakeHIP) && results[R]->dEdxTemplates)?&HIPEmulator:NULL);
                 //results[R]->HdedxVsPVsEta->Fill(track->p(),dedxObj.dEdx(),track->eta(),prescale);
                 //

                results[R]->HdedxVsP    ->Fill(track->p(), dedxObj.dEdx() ,prescale);

                if(fabs(track->eta()) < 0.91) results[R]->HdedxVsP_Eta0p91->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.91 <= fabs(track->eta()) && fabs(track->eta()) < 1.74) results[R]->HdedxVsP_0p91Eta1p74->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.74 <= fabs(track->eta())) results[R]->HdedxVsP_1p74Eta->Fill(track->p(),dedxObj.dEdx(),prescale);

		if (isSignal) results[R]->HdedxVsP->SetBins(1000, 0, 2400, results[R]->isDiscrim?1000:2000, 0, results[R]->isDiscrim?1.0:30); // if it's signal sample increase axis range

                if(pileup<20) results[R]->HdedxVsP_PUlower20->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(pileup>=35) results[R]->HdedxVsP_PUhigher35->Fill(track->p(), dedxObj.dEdx() ,prescale);

                if(pileup<10) results[R]->HdedxVsP_PU_0_10->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=10 && pileup<15) results[R]->HdedxVsP_PU_10_15->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=15 && pileup<20) results[R]->HdedxVsP_PU_15_20->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=20 && pileup<25) results[R]->HdedxVsP_PU_20_25->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=25 && pileup<30) results[R]->HdedxVsP_PU_25_30->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=30 && pileup<35) results[R]->HdedxVsP_PU_30_35->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=35 && pileup<40) results[R]->HdedxVsP_PU_35_40->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=40 && pileup<50) results[R]->HdedxVsP_PU_40_50->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(pileup>=50) results[R]->HdedxVsP_PU_50_inf->Fill(track->p(), dedxObj.dEdx(), prescale);


                if(lumi<1000) results[R]->HdedxVsP_instLumi_0_1000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=1000 && lumi<2000) results[R]->HdedxVsP_instLumi_1000_2000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=2000 && lumi<3000) results[R]->HdedxVsP_instLumi_2000_3000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=3000 && lumi<4000) results[R]->HdedxVsP_instLumi_3000_4000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=4000 && lumi<5000) results[R]->HdedxVsP_instLumi_4000_5000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=5000 && lumi<6000) results[R]->HdedxVsP_instLumi_5000_6000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=6000 && lumi<7000) results[R]->HdedxVsP_instLumi_6000_7000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=7000 && lumi<8000) results[R]->HdedxVsP_instLumi_7000_8000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=8000 && lumi<9000) results[R]->HdedxVsP_instLumi_8000_9000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=9000 && lumi<10000) results[R]->HdedxVsP_instLumi_9000_10000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=10000 && lumi<11000) results[R]->HdedxVsP_instLumi_10000_11000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=11000 && lumi<12000) results[R]->HdedxVsP_instLumi_11000_12000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=12000 && lumi<13000) results[R]->HdedxVsP_instLumi_12000_13000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=13000 && lumi<14000) results[R]->HdedxVsP_instLumi_13000_14000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=14000 && lumi<15000) results[R]->HdedxVsP_instLumi_14000_15000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=15000 && lumi<16000) results[R]->HdedxVsP_instLumi_15000_16000->Fill(track->p(), dedxObj.dEdx(), prescale);
                if(lumi>=16000) results[R]->HdedxVsP_instLumi_16000_inf->Fill(track->p(), dedxObj.dEdx(), prescale);

                h_trackPhi->Fill(track->phi(),prescale);
                //EtaPhi_dEdx0p5->Fill(track->eta(),track->phi(),prescale);
                //if(track->p()>5) EtaPhi_dEdx0p5_Phigher5->Fill(track->eta(),track->phi(),prescale);
                //EtaPhiVsP_dEdx0p5->Fill(track->eta(),track->phi(),track->p(),prescale);


                /*if(fabs(track->eta()) < 0.1) results[R]->HdedxVsP_Eta0p1->Fill(track->p(), dedxObj.dEdx() ,prescale);                
                if(0.1 <= fabs(track->eta()) && fabs(track->eta()) < 0.2) results[R]->HdedxVsP_0p1Eta0p2->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.2 <= fabs(track->eta()) && fabs(track->eta()) < 0.3) results[R]->HdedxVsP_0p2Eta0p3->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.3 <= fabs(track->eta()) && fabs(track->eta()) < 0.4) results[R]->HdedxVsP_0p3Eta0p4->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.4 <= fabs(track->eta()) && fabs(track->eta()) < 0.5) results[R]->HdedxVsP_0p4Eta0p5->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.5 <= fabs(track->eta()) && fabs(track->eta()) < 0.6) results[R]->HdedxVsP_0p5Eta0p6->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.6 <= fabs(track->eta()) && fabs(track->eta()) < 0.7) results[R]->HdedxVsP_0p6Eta0p7->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.7 <= fabs(track->eta()) && fabs(track->eta()) < 0.8) results[R]->HdedxVsP_0p7Eta0p8->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.8 <= fabs(track->eta()) && fabs(track->eta()) < 0.9) results[R]->HdedxVsP_0p8Eta0p9->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(0.9 <= fabs(track->eta()) && fabs(track->eta()) < 1.0) results[R]->HdedxVsP_0p9Eta1p0->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.0 <= fabs(track->eta()) && fabs(track->eta()) < 1.1) results[R]->HdedxVsP_1p0Eta1p1->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.1 <= fabs(track->eta()) && fabs(track->eta()) < 1.2) results[R]->HdedxVsP_1p1Eta1p2->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.2 <= fabs(track->eta()) && fabs(track->eta()) < 1.3) results[R]->HdedxVsP_1p2Eta1p3->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.3 <= fabs(track->eta()) && fabs(track->eta()) < 1.4) results[R]->HdedxVsP_1p3Eta1p4->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.4 <= fabs(track->eta()) && fabs(track->eta()) < 1.5) results[R]->HdedxVsP_1p4Eta1p5->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.5 <= fabs(track->eta()) && fabs(track->eta()) < 1.6) results[R]->HdedxVsP_1p5Eta1p6->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.6 <= fabs(track->eta()) && fabs(track->eta()) < 1.7) results[R]->HdedxVsP_1p6Eta1p7->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.7 <= fabs(track->eta()) && fabs(track->eta()) < 1.8) results[R]->HdedxVsP_1p7Eta1p8->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.8 <= fabs(track->eta()) && fabs(track->eta()) < 1.9) results[R]->HdedxVsP_1p8Eta1p9->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(1.9 <= fabs(track->eta()) && fabs(track->eta()) < 2.0) results[R]->HdedxVsP_1p9Eta2p0->Fill(track->p(), dedxObj.dEdx() ,prescale);
                if(2.0 <= fabs(track->eta()) && fabs(track->eta()) < 2.1) results[R]->HdedxVsP_2p0Eta2p1->Fill(track->p(), dedxObj.dEdx() ,prescale);*/
                
                
                
                results[R]->HdedxVsEta_highdEdx->Fill(track->eta(),dedxObj.dEdx(),prescale);
                if(GetMass(track->p(),dedxObj.dEdx(), results[R]->Kconst, results[R]->Cconst) > 0.938-0.20 
                        && GetMass(track->p(),dedxObj.dEdx(), results[R]->Kconst, results[R]->Cconst) < 0.938+0.20 
                        && track->p()>=0.9 && track->p()<1.0) 
                {
                    results[R]->HdedxVsEta_0p9P1p0->Fill(track->eta(),dedxObj.dEdx(),prescale);
                }


                if(track->p()>5) results[R]->HdedxVsEta_5p0P->Fill(track->eta(),dedxObj.dEdx(),prescale);
                if(track->p()>5) results[R]->HdedxVsAbsEta_5p0P->Fill(fabs(track->eta()),dedxObj.dEdx(),prescale);
                    
   //             results[R]->HdedxVsQP   ->Fill(track->p()*track->charge(), dedxObj.dEdx() );
   //             results[R]->HdedxVsP_NS ->Fill(track->p(), dedxObj.dEdx(), dedxObj.numberOfSaturatedMeasurements() );

                if(track->pt()>10 && track->pt()<45 && dedxObj.numberOfMeasurements()>=(results[R]->useStrip?7:3) ){
                  results[R]->HdedxVsEtaProfile->Fill(track->eta(), dedxObj.dEdx() ,prescale);
                  results[R]->HdedxVsEta->Fill(track->eta(), dedxObj.dEdx() ,prescale);
                  results[R]->HNOMVsEtaProfile->Fill(track->eta(),dedxObj.numberOfMeasurements() ,prescale);
                  results[R]->HNOSVsEtaProfile->Fill(track->eta(),dedxObj.numberOfSaturatedMeasurements() ,prescale);
                  results[R]->HNOMSVsEtaProfile->Fill(track->eta(),dedxObj.numberOfMeasurements() - dedxObj.numberOfSaturatedMeasurements() ,prescale);
                }

                if(fabs(track->eta())>2.1) continue;
                if((int)dedxObj.numberOfMeasurements()<(results[R]->useStrip?10:3))continue;
//                if(track->found()<10) continue; // we cut on total number of hits instead of valid measurements

                if(track->pt()>5){
                   results[R]->HdedxVsNOH->Fill(track->found(), dedxObj.dEdx(),prescale);
                   results[R]->HNOMVsdEdxProfile->Fill(dedxObj.dEdx(), (int)dedxObj.numberOfMeasurements(),prescale);
                   results[R]->HNOMVsdEdx->Fill(dedxObj.dEdx(), (int)dedxObj.numberOfMeasurements(),prescale);
                   results[R]->HdedxMIP  ->Fill(dedxObj.dEdx(),prescale);
                   results[R]->HP->Fill(track->p(),prescale);

                     DeDxData dedxObj4  = computedEdx(dedxHits, dEdxSF, results[R]->dEdxTemplates, results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside,  4, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, (results[R]->fakeHIP && !(SuppressFakeHIP))?&HIPemulator:NULL);
                     DeDxData dedxObj8  = computedEdx(dedxHits, dEdxSF, results[R]->dEdxTemplates, results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside,  8, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, (results[R]->fakeHIP && !(SuppressFakeHIP))?&HIPemulator:NULL);
                     DeDxData dedxObj12 = computedEdx(dedxHits, dEdxSF, results[R]->dEdxTemplates, results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 12, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, (results[R]->fakeHIP && !(SuppressFakeHIP))?&HIPemulator:NULL);
//                   DeDxData dedxObj4  = computeDeDx (dedxHits, dEdxSF, results[R]->dEdxTemplates, true, results[R]->CCFunction, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside,  4, false, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, results[R]->trimPixel, (results[R]->fakeHIP && !(SuppressFakeHIP) && results[R]->dEdxTemplates)?&HIPEmulator:NULL);
//                   DeDxData dedxObj8  = computeDeDx (dedxHits, dEdxSF, results[R]->dEdxTemplates, true, results[R]->CCFunction, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside,  8, false, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, results[R]->trimPixel, (results[R]->fakeHIP && !(SuppressFakeHIP) && results[R]->dEdxTemplates)?&HIPEmulator:NULL);
//                   DeDxData dedxObj12 = computeDeDx (dedxHits, dEdxSF, results[R]->dEdxTemplates, true, results[R]->CCFunction, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 12, false, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, results[R]->trimPixel, (results[R]->fakeHIP && !(SuppressFakeHIP) && results[R]->dEdxTemplates)?&HIPEmulator:NULL);

                   results[R]->HdedxMIP4 ->Fill(dedxObj4 .dEdx(),prescale);
                   results[R]->HdedxMIP8 ->Fill(dedxObj8 .dEdx(),prescale);
                   results[R]->HdedxMIP12->Fill(dedxObj12.dEdx(),prescale);
                }
                if(fabs(track->eta())<0.4){
                   results[R]->HdedxVsPProfile  ->Fill(track->p(), dedxObj  .dEdx() ,prescale);
                }


                DeDxData dedxObjEstim   = results[R]->isEstim?dedxObj:computedEdx(dedxHits, dEdxSF, NULL                     , results[R]->usePixel, results[R]->useClusterCleaning, false, results[R]->useTrunc, results[R]->TrackerGains, results[R]->useStrip, results[R]->mustBeInside, 99, results[R]->correctFEDSat, results[R]->crossTalkInvAlgo, results[R]->dropLowerDeDxValue, (results[R]->fakeHIP && !(SuppressFakeHIP))?&HIPemulator:NULL);
                double Mass = GetMass(track->p(),dedxObjEstim.dEdx(), results[R]->Kconst, results[R]->Cconst);
                //if (Mass > 1.8756-0.20 && Mass < 1.8756+0.20 && dedxObjEstim.dEdx() > 5){// proton candidates
                if (Mass > 0.93827-0.20 && Mass < 0.93827+0.20 && dedxObjEstim.dEdx() > 4){// proton candidates
                   results[R]->HdedxVsPSyst->Fill(track->p(), dedxObj.dEdx() ,prescale);
		   
                   if (results[R]->isEstim && dedxObj.dEdx() > 7){
                      for(unsigned int h=0;h<dedxHits->size();h++){
                         DetId detid(dedxHits->detId(h));
                         double scaleFactor = dEdxSF[0];
                         if (detid.subdetId()<3) scaleFactor *= dEdxSF[1];
                         double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;
                       
                         int moduleGeometry = 0; // underflow bin -- debug purposes
                         if(detid.subdetId()>=3){ SiStripDetId SSdetId(detid); moduleGeometry = SSdetId.moduleGeometry(); if (moduleGeometry==15) {cerr << "ERROR! There is no SiStrip geometry 15!" << endl; exit (EXIT_FAILURE);}}
                         else if(detid.subdetId()<3){moduleGeometry = 15;} // 15 is for pixel
                         if(!results[R]->usePixel && detid.subdetId() <3)continue; // skip pixels
                         if(!results[R]->useStrip && detid.subdetId()>=3)continue; // skip strips
                         if(results[R]->mustBeInside && !isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId()>=3?dedxHits->stripCluster(h):NULL) )continue;
//                         if(results[R]->removeCosmics){ if (isCompatibleWithCosmic(track, vertexColl))continue;} //don't consider hits, which belong to cosmic tracks
                         if(results[R]->useClusterCleaning && detid.subdetId()>=3 && !clusterCleaning(dedxHits->stripCluster(h), results[R]->crossTalkInvAlgo)) continue; //if it fails clusterCleaning, skip it!
                       
                         int charge = dedxHits->charge(h);
                         if (detid.subdetId()>=3 && results[R]->crossTalkInvAlgo!=0){
                            vector <int> amps = CrossTalkInv(convert(dedxHits->stripCluster(h)->amplitudes()), 0.10, 0.04, true);
                            charge = std::accumulate(amps.begin(), amps.end(), 0);
                         }
                         double ChargeOverPathlength = scaleFactor*Norm*charge/dedxHits->pathlength(h);

//                         results[R]->HProtonHit->Fill(ChargeOverPathlength);
			 if (detid.subdetId()<3) results[R]->HProtonHitPO->Fill(ChargeOverPathlength,prescale);
			 else                    results[R]->HProtonHitSO->Fill(ChargeOverPathlength,prescale);
                      }
                   }
                }

                if(results[R]->isEstim && dedxObj.dEdx()>results[R]->Cconst + 3.0){  //mass can only be computed for dEdx estimators
                   if(track->p()<3.0){      results[R]->HMass->Fill(Mass,prescale);
                   }else{                   results[R]->HMassHSCP->Fill(Mass,prescale);
                   }
                }
                // FIXME DEBUG -- only for Ias -- to check what's going on during cluster cleaning
                /*if(results[R]->isDiscrim){
                   if (!results[R]->useClusterCleaning && !results[R]->crossTalkInvAlgo && results[R]->Name=="Ias_SO_in_noC" && dedxObj.dEdx()>0.9){ lockOnTrack=true; dEdxDebug = dedxObj.dEdx(); }
                   if (lockOnTrack && results[R]->useClusterCleaning && !results[R]->crossTalkInvAlgo && results[R]->Name=="Ias_SO_in_noC_CC" && dedxObj.dEdx()<0.9){
                      std::cout << "Track Ias dropped below 0.9 after ClusterCleaning!" << std::endl;
                      std::cout << "Track index = " << c << "\tdEdx difference = " << dEdxDebug - dedxObj.dEdx() << std::endl;
                      std::cout << "Failed clusters:" << std::endl;
                      lockOnTrack=false; dEdxDebug=0; //reset to the original values
                      // have to know which are the clusters that are removed, and how they look like!
                      // print out rejected clusters
                      for (unsigned int h=0;h<dedxHits->size();h++){
                         uint8_t exitCode;
                         if (!clusterCleaning(dedxHits->stripCluster(h), 0, &exitCode)){
                            std::cout << "\t[ ";
                            for (unsigned int digi_i=0;digi_i<dedxHits->stripCluster(h)->amplitudes().size();digi_i++)
                               std::cout << static_cast<unsigned int>(dedxHits->stripCluster(h)->amplitudes()[digi_i]) << " ";
                            std::cout << "]\t"; printClusterCleaningMessage (exitCode);
                         }
https://github.com/HgithubNLto2L2Q/HNgithubL/blob/master/HeavyNeutralLeptonAnalysis/plugins/HeavyNeutralLeptonAnalysis.cc                      }
                   }
                }*/
             }
         }
      }printf("\n");
      delete file;
   }
   std::cout << "Ncluster_lowP(<1GeV): " << Ncluster_lowP 
       << " Nsatcluster_lowP(<1GeV): " << Nsatcluster_lowP 
       << " ratio: " << (float)Nsatcluster_lowP/(float)Ncluster_lowP 
       << std::endl;
   OutputHisto->Write();
   OutputHisto->Close();  
}

double DistToHSCP (const reco::TrackRef& track, const std::vector<reco::GenParticle>& genColl){
   if(track.isNull())return false; // FIXME does this make sense? returning false to a double function?

   double RMin = 9999;
   for(unsigned int g=0;g<genColl.size();g++){
      if(genColl[g].pt()<5)continue;
      if(genColl[g].status()!=1)continue;
      int AbsPdg=abs(genColl[g].pdgId());
      if(AbsPdg<1000000 && AbsPdg!=17)continue;
      double dR = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
      if(dR<RMin)RMin=dR;
   }
   return RMin;
}

bool isCompatibleWithCosmic (const reco::TrackRef& track, const std::vector<reco::Vertex>& vertexColl){
   for (unsigned int vertex_i=0;vertex_i<vertexColl.size();vertex_i++){
      if(fabs(track->dz (vertexColl[vertex_i].position())) < 0.5 && fabs(track->dxy(vertexColl[vertex_i].position())) < 0.2)return false;
   }
   return true;
}

double GetMass (double P, double I, double K, double C){
   if (I-C<0) return -1;
   return sqrt((I-C)/K)*P;
}


