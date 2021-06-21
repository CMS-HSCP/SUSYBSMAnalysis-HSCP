// Original Author:  Loic Quertenmont

#ifndef SUSYBSMAnalysis_Analyzer_Tuple_h
#define SUSYBSMAnalysis_Analyzer_Tuple_h

struct Tuple {
  //=============================================================
  //      Declare Trees & Branches
  //=============================================================
   TTree*       Tree;
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
  //      Declare Histograms
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
   std::map<std::string,TH1F*> BS_Pt_Binned; //TH1F*  BS_Pt_Binned[MaxPredBins];
   std::map<std::string,TH1F*> BS_TOF_Binned; //TH1F*  BS_TOF_Binned[MaxPredBins]; 

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

  //Prediction histograms for muon only analysis which is binned depending on eta nd number of muon stations
  std::map<std::string,TH1D*> H_B_Binned; //TH1D* H_B_Binned[MaxPredBins];
  std::map<std::string,TH1D*> H_D_Binned; //TH1D* H_D_Binned[MaxPredBins];
  std::map<std::string,TH1D*> H_F_Binned; //TH1D* H_F_Binned[MaxPredBins];
  std::map<std::string,TH1D*> H_H_Binned; //TH1D* H_H_Binned[MaxPredBins];

  TH1D*  HCuts_Pt;
  TH1D*  HCuts_Is;
  TH1D*  HCuts_TOF;

  TH1D*  Hist_Pt ;
  TH1D*  Hist_Is  ;
  TH1D*  Hist_TOF;

 //FIXME ------ To be modified for Number Of Hits (NOH)
  TH3D*  Pred_EtaP ;
  TH2D*  Pred_I    ;
  TH2D*  Pred_TOF  ;
  TH2D*  Pred_EtaB;
  TH2D*  Pred_EtaS;
  TH2D*  Pred_EtaS2;

  //pz

  TH2D* PDF_E_Eta;
  TH2D* PDF_A_Eta;
  TH3D* PDF_H_EtaMass;
  TH3D* PDF_G_EtaP;
  TH3D* PDF_C_EtaP;
  TH3D* PDF_F_EtaICK;
  TH3D* PDF_B_EtaICK;

  TH2D* PDF_E_Eta_Flip;
  TH2D* PDF_A_Eta_Flip;
  TH3D* PDF_H_EtaMass_Flip;
  TH3D* PDF_G_EtaP_Flip;
  TH3D* PDF_C_EtaP_Flip;
  TH3D* PDF_F_EtaICK_Flip;
  TH3D* PDF_B_EtaICK_Flip;


// end FIXME

  TH2D*  RegionD_P;
  TH2D*  RegionD_I;
  TH2D*  RegionD_Ias;
  TH2D*  RegionD_TOF;

  TH2D*  RegionH_Ias;

  TH1D* H_A_Flip;
  TH1D* H_B_Flip;
  TH1D* H_C_Flip;
  TH1D* H_D_Flip;
  TH1D* H_E_Flip;
  TH1D* H_F_Flip;
  TH1D* H_G_Flip;
  TH1D* H_H_Flip;

  std::map<std::string,TH1D*> H_B_Binned_Flip; //TH1D* H_B_Binned_Flip[MaxPredBins];
  std::map<std::string,TH1D*> H_D_Binned_Flip; //TH1D* H_D_Binned_Flip[MaxPredBins];
  std::map<std::string,TH1D*> H_F_Binned_Flip; //TH1D* H_F_Binned_Flip[MaxPredBins];
  std::map<std::string,TH1D*> H_H_Binned_Flip; //TH1D* H_H_Binned_Flip[MaxPredBins];


 //FIXME ------ To be modified for Number Of Hits (NOH)
  TH3D*  Pred_EtaP_Flip ;
  TH2D*  Pred_I_Flip    ;
  TH2D*  Pred_TOF_Flip  ;
  TH2D*  Pred_EtaB_Flip;
  TH2D*  Pred_EtaS_Flip;
  TH2D*  Pred_EtaS2_Flip;
// end FIXME

  TH2D*  RegionD_P_Flip;
  TH2D*  RegionD_I_Flip;
  TH2D*  RegionD_Ias_Flip;
  TH2D*  RegionD_TOF_Flip;

  TH2D*  RegionH_Ias_Flip;

  TH2D* H_D_DzSidebands;

  TH2F*  genrecopT;
  TH1F*  genlevelpT;
  TH1F*  genleveleta;
  TH1F*  genlevelbeta;

  TH1D*  CtrlPt_S1_Is;
  TH1D*  CtrlPt_S2_Is;
  TH1D*  CtrlPt_S3_Is;
  TH1D*  CtrlPt_S4_Is;

  TH1D*  CtrlIs_S1_TOF;
  TH1D*  CtrlIs_S2_TOF;
  TH1D*  CtrlIs_S3_TOF;
  TH1D*  CtrlIs_S4_TOF;

  TH1D*  CtrlIm_S1_TOF;
  TH1D*  CtrlIm_S2_TOF;
  TH1D*  CtrlIm_S3_TOF;
  TH1D*  CtrlIm_S4_TOF;

  TH1D*  CtrlPt_S1_Im;
  TH1D*  CtrlPt_S2_Im;
  TH1D*  CtrlPt_S3_Im;
  TH1D*  CtrlPt_S4_Im;

  TH1D*  CtrlPt_S1_TOF;
  TH1D*  CtrlPt_S2_TOF;
  TH1D*  CtrlPt_S3_TOF;
  TH1D*  CtrlPt_S4_TOF;

  std::map<std::string,TH1D*> CtrlPt_S1_TOF_Binned; //TH1D* CtrlPt_S1_TOF_Binned[MaxPredBins];
  std::map<std::string,TH1D*> CtrlPt_S2_TOF_Binned; //TH1D* CtrlPt_S2_TOF_Binned[MaxPredBins];
  std::map<std::string,TH1D*> CtrlPt_S3_TOF_Binned; //TH1D* CtrlPt_S3_TOF_Binned[MaxPredBins];
  std::map<std::string,TH1D*> CtrlPt_S4_TOF_Binned; //TH1D* CtrlPt_S4_TOF_Binned[MaxPredBins];
};

#endif