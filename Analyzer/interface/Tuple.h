// Original Author:  Loic Quertenmont
// Modification by Tamas Almos Vami

#ifndef SUSYBSMAnalysis_Analyzer_Tuple_h
#define SUSYBSMAnalysis_Analyzer_Tuple_h

struct Tuple {
  //=============================================================
  //      Declare Trees & Branches
  //=============================================================

  TTree* Tree;
  unsigned int NCuts;
  unsigned int Tree_Trig;
  unsigned int Tree_Run;
  unsigned long Tree_Event;
  unsigned int Tree_Lumi;
  unsigned int Tree_PileUp;
  unsigned int Tree_nofVertices;
  unsigned int Tree_Hscp;
  unsigned int Tree_nmuons;
  unsigned int Tree_njets;
  float Tree_Weight;
  float Tree_GeneratorWeight;
  bool Tree_HLT_Mu50;
  bool Tree_HLT_PFMET120_PFMHT120_IDTight;
  bool Tree_HLT_PFHT500_PFMET100_PFMHT100_IDTight;
  bool Tree_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  bool Tree_HLT_MET105_IsoTrk50;
  float Tree_CaloMET;
  float Tree_RecoPFMET;
  float Tree_RecoPFMHT;
  float Tree_HLTPFMET;
  float Tree_HLTPFMHT;
  float Tree_RecoPFMET_eta;
  float Tree_RecoPFMET_phi;
  float Tree_RecoPFMET_significance;
  float Tree_Muon1_Pt;
  float Tree_Muon1_eta;
  float Tree_Muon1_phi;
  float Tree_Muon2_Pt;
  float Tree_Muon2_eta;
  float Tree_Muon2_phi;

  std::vector<float> Tree_vect_mT;
  std::vector<bool> Tree_passCutPt55;
  std::vector<bool> Tree_passPreselection_noIsolation_noIh;
  std::vector<bool> Tree_passPreselection;
  std::vector<bool> Tree_passSelection;
  std::vector<float> Tree_Charge;
  std::vector<float> Tree_Pt;
  std::vector<float> Tree_PtErr;
  std::vector<float> Tree_Ias;
  std::vector<float> Tree_Ias_noTIBnoTIDno3TEC;
  std::vector<float> Tree_Ias_PixelOnly;
  std::vector<float> Tree_Ias_StripOnly;
  std::vector<float> Tree_Ih;
  std::vector<float> Tree_Ick;  //return (Ih-C)/K
  std::vector<float> Tree_Fmip;
  std::vector<float> Tree_ProbXY;
  std::vector<float> Tree_ProbXY_noL1;
  std::vector<float> Tree_ProbQ;
  std::vector<float> Tree_ProbQ_noL1;
  std::vector<float> Tree_ProbQ_dEdx;
  std::vector<float> Tree_Ndof;
  std::vector<float> Tree_Chi2;
  std::vector<bool>  Tree_isHighPurity;
  std::vector<float> Tree_muon_eta;
  std::vector<bool> Tree_isMuon;
  std::vector<int>  Tree_Muon_selector;
  std::vector<bool> Tree_isElectron;
  std::vector<bool> Tree_isChHadron;
  std::vector<bool> Tree_isNeutHadron;
  std::vector<float> Tree_ECAL_energy;
  std::vector<float> Tree_HCAL_energy;
  std::vector<float> Tree_TOF;
  std::vector<float> Tree_TOFErr;
  std::vector<unsigned int> Tree_TOF_ndof;
  std::vector<float> Tree_DTTOF;
  std::vector<float> Tree_DTTOFErr;
  std::vector<unsigned int> Tree_DTTOF_ndof;
  std::vector<float> Tree_CSCTOF;
  std::vector<float> Tree_CSCTOFErr;
  std::vector<unsigned int> Tree_CSCTOF_ndof;
  std::vector<float> Tree_Mass;
  std::vector<float> Tree_MassErr;
  std::vector<float> Tree_dZ;
  std::vector<float> Tree_dXY;
  std::vector<float> Tree_dR;
  std::vector<float> Tree_eta;
  std::vector<float> Tree_phi;
  std::vector<unsigned int> Tree_NOH;   //number of (valid) track pixel+strip hits
  std::vector<unsigned int> Tree_NOPH;  //number of (valid) track pixel hits
  std::vector<float> Tree_FOVH;         //fraction of valid track hits
  std::vector<unsigned int>
      Tree_NOMH;                  //number of missing hits from IP till last hit (excluding hits behind the last hit)
  std::vector<float> Tree_FOVHD;  //fraction of valid hits divided by total expected hits until the last one
  std::vector<unsigned int>
      Tree_NOM;  //number of dEdx hits (= #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used)
  std::vector<float> Tree_iso_TK;
  std::vector<float> Tree_iso_ECAL;
  std::vector<float> Tree_iso_HCAL;
  std::vector<float> Tree_track_PFIsolationR005_sumChargedHadronPt;
  std::vector<float> Tree_track_PFIsolationR005_sumNeutralHadronPt;
  std::vector<float> Tree_track_PFIsolationR005_sumPhotonPt;
  std::vector<float> Tree_track_PFIsolationR005_sumPUPt;

  std::vector<float> Tree_track_PFIsolationR01_sumChargedHadronPt;
  std::vector<float> Tree_track_PFIsolationR01_sumNeutralHadronPt;
  std::vector<float> Tree_track_PFIsolationR01_sumPhotonPt;
  std::vector<float> Tree_track_PFIsolationR01_sumPUPt;
 
  std::vector<float> Tree_track_PFIsolationR03_sumChargedHadronPt;
  std::vector<float> Tree_track_PFIsolationR03_sumNeutralHadronPt;
  std::vector<float> Tree_track_PFIsolationR03_sumPhotonPt;
  std::vector<float> Tree_track_PFIsolationR03_sumPUPt;

  std::vector<float> Tree_track_PFIsolationR05_sumChargedHadronPt;
  std::vector<float> Tree_track_PFIsolationR05_sumNeutralHadronPt;
  std::vector<float> Tree_track_PFIsolationR05_sumPhotonPt;
  std::vector<float> Tree_track_PFIsolationR05_sumPUPt;

  std::vector<float> Tree_muon_PFIsolationR03_sumChargedHadronPt;
  std::vector<float> Tree_muon_PFIsolationR03_sumNeutralHadronPt;
  std::vector<float> Tree_muon_PFIsolationR03_sumPhotonPt;
  std::vector<float> Tree_muon_PFIsolationR03_sumPUPt;

  std::vector<float> Tree_Ih_noL1;
  std::vector<float> Tree_Ih_15drop;
  std::vector<float> Tree_Ih_StripOnly;
  std::vector<float> Tree_Ih_StripOnly_15drop;
  std::vector<float> Tree_Ih_SaturationCorrectionFromFits;
  std::vector<std::vector<float>> Tree_clust_charge;  //dedx charge -> either strip or pixel
  std::vector<std::vector<float>> Tree_clust_pathlength;
  std::vector<std::vector<unsigned int>> Tree_clust_nstrip;
  std::vector<std::vector<bool>> Tree_clust_sat254;
  std::vector<std::vector<bool>> Tree_clust_sat255;
  std::vector<std::vector<uint32_t>> Tree_clust_detid;
  std::vector<std::vector<bool>> Tree_clust_isStrip;  //is it a SiStrip cluster?
  std::vector<std::vector<bool>> Tree_clust_isPixel;  //is it a Pixel hit?
  std::vector<float> Tree_GenId;
  std::vector<float> Tree_GenCharge;
  std::vector<float> Tree_GenMass;
  std::vector<float> Tree_GenPt;
  std::vector<float> Tree_GenEta;
  std::vector<float> Tree_GenPhi;

  TTree* GenTree;
  unsigned int GenTree_Run;
  unsigned int GenTree_Event;
  unsigned int GenTree_Lumi;
  unsigned int GenTree_Hscp;
  float GenTree_Weight;
  float GenTree_GeneratorWeight;
  std::vector<float> GenTree_GenId;
  std::vector<float> GenTree_GenCharge;
  std::vector<float> GenTree_GenMass;
  std::vector<float> GenTree_GenPt;
  std::vector<float> GenTree_GenEta;
  std::vector<float> GenTree_GenPhi;


  //=============================================================
  //      Declare Histograms
  //=============================================================

  TH1F* CutFlow_nHSCP;


  TH2F* Mass;
  TH2F* MassTOF;
  TH2F* MassComb;
  TH2F* MaxEventMass;

  TH2F* Mass_SystP;
  TH2F* MassTOF_SystP;
  TH2F* MassComb_SystP;
  TH2F* MaxEventMass_SystP;

  TH2F* Mass_SystI;
  TH2F* MassTOF_SystI;
  TH2F* MassComb_SystI;
  TH2F* MaxEventMass_SystI;

  TH2F* Mass_SystM;
  TH2F* MassTOF_SystM;
  TH2F* MassComb_SystM;
  TH2F* MaxEventMass_SystM;

  TH2F* Mass_SystT;
  TH2F* MassTOF_SystT;
  TH2F* MassComb_SystT;
  TH2F* MaxEventMass_SystT;

  TH2F* Mass_SystPU;
  TH2F* MassTOF_SystPU;
  TH2F* MassComb_SystPU;
  TH2F* MaxEventMass_SystPU;

  TH2F* Mass_SystHUp;
  TH2F* MassTOF_SystH;
  TH2F* MassComb_SystHUp;
  TH2F* MaxEventMass_SystHUp;

  TH2F* Mass_SystHDown;
  TH2F* MassComb_SystHDown;
  TH2F* MaxEventMass_SystHDown;

  TH2F* Mass_Flip;
  TH2F* MassTOF_Flip;
  TH2F* MassComb_Flip;

  TProfile* IntLumi;
  TProfile* XSection;
  TH1F* EventsTotal;
  TH1F* TotalE;
  TH1F* TotalEPU;
  TH1F* TotalTE;
  TH1F* Total;
  TH1F* V3D;
  TH1F* N1Eta;
  TH1F* N1Chi2PerNdof;
  TH1F* N1Qual;
  TH1F* N1TNOH;
  TH1F* N1TNOPH;
  TH1F* N1TNOHFraction;
  TH1F* N1TNOM;
  TH1F* nDof;
  TH1F* tofError;
  TH1F* N1MPt;
  TH1F* N1MIh;
  TH1F* MTOF;
  TH1F* TIsol;
  TH1F* N1EIsol;
  TH1F* N1SumpTOverpT;
  TH1F* Pt;
  TH1F* N1PtErrOverPt;
  TH1F* I;
  TH1F* TOF;
  TH1F* HSCPE;
  TH1F* NVTrack;
  TH1F* N1Stations;
  TH1F* N1Dxy;
  TH1F* N1Dz;
  TH1F* N1SegSep;
  TH1F* FailDz;
  TH1F* Basic;
  TH1F* N1ProbQ;
  TH1F* ProbQNoL1;
  TH1F* N1ProbXY;
  TH1F* ProbXYNoL1;

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

  TH1F* BS_ProbQ;
  TH1F* BS_ProbXY;
  TH1F* BS_ProbQNoL1;
  TH1F* BS_ProbXYNoL1;
  
  TH1F* Beta_Matched;
  TH1F* Beta_PreselectedA;
  TH1F* Beta_PreselectedB;
  TH1F* Beta_PreselectedC;
  TH2F* Beta_SelectedP;
  TH2F* Beta_SelectedI;
  TH2F* Beta_SelectedT;

  TH1F* BS_V3D;
  TH1F* BS_Chi2PerNdof;
  TH1F* BS_Qual;
  TH1F* BS_TNOH;
  TH1F* BS_TNOH_PUA;
  TH1F* BS_TNOH_PUB;
  TH1F* BS_TNOHFraction;
  TH1F* BS_TNOPH;
  TH1F* BS_TNOHFractionTillLast;
  TH1F* BS_TNOMHTillLast;
  TH1F* BS_Eta;
  TH1F* BS_TNOM;
  TH1F* BS_TNOM_PUA;
  TH1F* BS_TNOM_PUB;
  TProfile* BS_NOMoNOHvsPV;
  TH1F* BS_nDof;
  TH1F* BS_TOFError;
  TH1F* BS_PtErrOverPt;
  TH1F* BS_PtErrOverPt2;
  TH1F* BS_MPt;
  TH1F* BS_MIs;
  TH1F* BS_MIh;
  TH1F* BS_MTOF;
  TH1F* BS_TIsol;
  TH1F* BS_EIsol;
  TH1F* BS_SumpTOverpT;
  TH1F* BS_dR_NVTrack;
  TH1F* BS_MatchedStations;
  TH1F* BS_NVertex;
  TH1F* BS_NVertex_NoEventWeight;
  TH1F* BS_PV;
  TH1F* BS_PV_NoEventWeight;
  TH1F* BS_dzAll;
  TH1F* BS_dxyAll;
  TH1F* BS_dzMinv3d;
  TH1F* BS_dxyMinv3d;
  TH1F* BS_SegSep;
  TH1F* BS_SegMinPhiSep;
  TH1F* BS_SegMinEtaSep;
  TH1F* BS_SegMinEtaSep_FailDz;
  TH1F* BS_SegMinEtaSep_PassDz;
  TH1F* BS_Dz_FailSep;
  TH1F* BS_InnerInvPtDiff;
  TH1F* BS_Phi;
  TH1F* BS_TimeAtIP;
  TH1F* BS_OpenAngle;
  TH1F* BS_OpenAngle_Cosmic;

  TH1F* BS_Pt_FailDz;
  TH1F* BS_Pt_FailDz_DT;
  TH1F* BS_Pt_FailDz_CSC;
  TH1F* BS_TOF_FailDz;
  TH1F* BS_TOF_FailDz_DT;
  TH1F* BS_TOF_FailDz_CSC;
  TH1F* BS_Dxy;
  TH1F* BS_Dxy_Cosmic;
  TH1F* BS_Dz;
  TH1F* BS_Dz_Cosmic;
  TH1F* BS_Dz_CSC;
  TH1F* BS_Dz_DT;
  std::map<std::string, TH1F*> BS_Pt_Binned;   //TH1F*  BS_Pt_Binned[MaxPredBins];
  std::map<std::string, TH1F*> BS_TOF_Binned;  //TH1F*  BS_TOF_Binned[MaxPredBins];

  TH1F* BS_LastHitDXY;
  TH1F* BS_LastHitD3D;
  TH2F* BS_PtErrOverPtVsPtErrOverPt2;

  TH1F* CutFlow;
  TH1F* CutFlowProbQLast;
  TH1F* CutFlowProbQFirst;

  TH2F* AS_Eta_RegionA;
  TH2F* AS_Eta_RegionB;
  TH2F* AS_Eta_RegionC;
  TH2F* AS_Eta_RegionD;
  TH2F* AS_Eta_RegionE;
  TH2F* AS_Eta_RegionF;
  TH2F* AS_Eta_RegionG;
  TH2F* AS_Eta_RegionH;

  TH1F* BS_P;
  TH2F* AS_P;
  TH1F* BS_Pt;
  TH2F* AS_Pt;
  TH1F* BS_Pt_PUA;
  TH1F* BS_Pt_PUB;
  TH1F* BS_Pt_DT;
  TH1F* BS_Pt_CSC;
  TH1F* BS_Is;
  TH2F* AS_Is;
  TH1F* BS_Is_PUA;
  TH1F* BS_Is_PUB;
  TH1F* BS_Ih;
  TH2F* AS_Ih;
  TH1F* BS_Ih_PUA;
  TH1F* BS_Ih_PUB;
  TH1F* BS_TOF;
  TH2F* AS_TOF;
  TH1F* BS_TOF_PUA;
  TH1F* BS_TOF_PUB;
  TH1F* BS_TOF_DT;
  TH1F* BS_TOF_CSC;
  TH1F* BS_Is_Cosmic;
  TH1F* BS_Pt_Cosmic;

  TH2F* BS_EtaIs;   //TH3F*  AS_EtaIs;
  TH2F* BS_EtaIh;   //TH3F*  AS_EtaIh;
  TH2F* BS_EtaP;    //TH3F*  AS_EtaP;
  TH2F* BS_EtaPt;   //TH3F*  AS_EtaPt;
  TH2F* BS_EtaTOF;  //TH3F*  AS_EtaTOF;
  TH2F* BS_EtaDz;
  TH2F* BS_EtaNBH;  // number of bad hits vs Eta

  TH2F* BS_PIs;
  TH3F* AS_PIs;
  TH2F* BS_PIhHD;
  TH2F* BS_PIh;
  TH3F* AS_PIh;
  TH2F* BS_PtIs;
  TH3F* AS_PtIs;
  TH2F* BS_PtIh;
  TH3F* AS_PtIh;
  TH2F* BS_PtTOF;
  TH2F* BS_TOFIs;
  TH3F* AS_TOFIs;
  TH2F* BS_TOFIh;
  TH3F* AS_TOFIh;

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
  std::map<std::string, TH1D*> H_B_Binned;  //TH1D* H_B_Binned[MaxPredBins];
  std::map<std::string, TH1D*> H_D_Binned;  //TH1D* H_D_Binned[MaxPredBins];
  std::map<std::string, TH1D*> H_F_Binned;  //TH1D* H_F_Binned[MaxPredBins];
  std::map<std::string, TH1D*> H_H_Binned;  //TH1D* H_H_Binned[MaxPredBins];

  TH1D* HCuts_Pt;
  TH1D* HCuts_Is;
  TH1D* HCuts_TOF;

  TH1D* Hist_Pt;
  TH1D* Hist_Is;
  TH1D* Hist_TOF;

  //FIXME ------ To be modified for Number Of Hits (NOH)
  TH3D* Pred_EtaP;
  TH2D* Pred_I;
  TH2D* Pred_TOF;
  TH2D* Pred_EtaB;
  TH2D* Pred_EtaS;
  TH2D* Pred_EtaS2;

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

  TH2D* RegionD_P;
  TH2D* RegionD_I;
  TH2D* RegionD_Ias;
  TH2D* RegionD_TOF;

  TH2D* RegionH_Ias;

  TH1D* H_A_Flip;
  TH1D* H_B_Flip;
  TH1D* H_C_Flip;
  TH1D* H_D_Flip;
  TH1D* H_E_Flip;
  TH1D* H_F_Flip;
  TH1D* H_G_Flip;
  TH1D* H_H_Flip;

  std::map<std::string, TH1D*> H_B_Binned_Flip;  //TH1D* H_B_Binned_Flip[MaxPredBins];
  std::map<std::string, TH1D*> H_D_Binned_Flip;  //TH1D* H_D_Binned_Flip[MaxPredBins];
  std::map<std::string, TH1D*> H_F_Binned_Flip;  //TH1D* H_F_Binned_Flip[MaxPredBins];
  std::map<std::string, TH1D*> H_H_Binned_Flip;  //TH1D* H_H_Binned_Flip[MaxPredBins];

  //FIXME ------ To be modified for Number Of Hits (NOH)
  TH3D* Pred_EtaP_Flip;
  TH2D* Pred_I_Flip;
  TH2D* Pred_TOF_Flip;
  TH2D* Pred_EtaB_Flip;
  TH2D* Pred_EtaS_Flip;
  TH2D* Pred_EtaS2_Flip;
  // end FIXME

  TH2D* RegionD_P_Flip;
  TH2D* RegionD_I_Flip;
  TH2D* RegionD_Ias_Flip;
  TH2D* RegionD_TOF_Flip;

  TH2D* RegionH_Ias_Flip;

  TH2D* H_D_DzSidebands;

  TH2F* genrecopT;
  TH1F* genlevelpT;
  TH1F* genleveleta;
  TH1F* genlevelbeta;

  TH1D* CtrlPt_S1_Is;
  TH1D* CtrlPt_S2_Is;
  TH1D* CtrlPt_S3_Is;
  TH1D* CtrlPt_S4_Is;

  TH1D* CtrlIs_S1_TOF;
  TH1D* CtrlIs_S2_TOF;
  TH1D* CtrlIs_S3_TOF;
  TH1D* CtrlIs_S4_TOF;

  TH1D* CtrlIh_S1_TOF;
  TH1D* CtrlIh_S2_TOF;
  TH1D* CtrlIh_S3_TOF;
  TH1D* CtrlIh_S4_TOF;

  TH1D* CtrlPt_S1_Ih;
  TH1D* CtrlPt_S2_Ih;
  TH1D* CtrlPt_S3_Ih;
  TH1D* CtrlPt_S4_Ih;

  TH1D* CtrlPt_S1_TOF;
  TH1D* CtrlPt_S2_TOF;
  TH1D* CtrlPt_S3_TOF;
  TH1D* CtrlPt_S4_TOF;

  std::map<std::string, TH1D*> CtrlPt_S1_TOF_Binned;  //TH1D* CtrlPt_S1_TOF_Binned[MaxPredBins];
  std::map<std::string, TH1D*> CtrlPt_S2_TOF_Binned;  //TH1D* CtrlPt_S2_TOF_Binned[MaxPredBins];
  std::map<std::string, TH1D*> CtrlPt_S3_TOF_Binned;  //TH1D* CtrlPt_S3_TOF_Binned[MaxPredBins];
  std::map<std::string, TH1D*> CtrlPt_S4_TOF_Binned;  //TH1D* CtrlPt_S4_TOF_Binned[MaxPredBins];
};

#endif
