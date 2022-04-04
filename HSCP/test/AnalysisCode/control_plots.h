// Original Author:  Jessica Prisciandaro
#include "CommonTools/UtilAlgos/interface/TFileService.h"

//#include "Analysis_Global.h"
//#include "Analysis_Samples.h"

//define a container for all plots that should be produced per sample
struct controlPlots {
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

  TH1F* BS_V3D;
  TH1F* BS_Chi2;
  TH1F* BS_Qual;
  TH1F* BS_TNOH;
  TH1F* BS_TNOH_PUA;
  TH1F* BS_TNOH_PUB;
  TH1F* BS_TNOHFraction;
  TH1F* BS_TNOPH;
  TH1F* BS_TNOHFractionTillLast;
  TH1F* BS_TNOMHTillLast;
  TH1F* BS_Eta;  //J.
  TH1F* BS_TNOM;
  TH1F* BS_TNOM_PUA;
  TH1F* BS_TNOM_PUB;
  TProfile* BS_NOMoNOHvsPV;
  TH1F* BS_nDof;
  TH1F* BS_TOFError;
  TH1F* BS_Pterr;
  TH1F* BS_MPt;
  TH1F* BS_MIs;
  TH1F* BS_MIm;
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
  TH1F* BS_Pt_Binned[MaxPredBins];
  TH1F* BS_TOF_Binned[MaxPredBins];

  TH1F* BS_LastHitDXY;
  TH1F* BS_LastHitD3D;

  TH1F* BS_P;
  TH1F* BS_Pt;
  TH1F* BS_Pt_PUA;
  TH1F* BS_Pt_PUB;
  TH1F* BS_Pt_DT;
  TH1F* BS_Pt_CSC;
  TH1F* BS_Is;
  TH1F* BS_Is_PUA;
  TH1F* BS_Is_PUB;
  TH1F* BS_Im;
  TH1F* BS_Im_PUA;
  TH1F* BS_Im_PUB;
  TH1F* BS_TOF;
  TH1F* BS_TOF_PUA;
  TH1F* BS_TOF_PUB;
  TH1F* BS_TOF_DT;
  TH1F* BS_TOF_CSC;
  TH1F* BS_Is_Cosmic;
  TH1F* BS_Pt_Cosmic;

  TH2F* BS_EtaIs;   //TH3F*  AS_EtaIs;
  TH2F* BS_EtaIm;   //TH3F*  AS_EtaIm;
  TH2F* BS_EtaP;    //TH3F*  AS_EtaP;
  TH2F* BS_EtaPt;   //TH3F*  AS_EtaPt;
  TH2F* BS_EtaTOF;  //TH3F*  AS_EtaTOF;
  TH2F* BS_EtaDz;

  TH2F* BS_PIs;
  TH3F* AS_PIs;
  TH2F* BS_PImHD;
  TH2F* BS_PIm;
  TH3F* AS_PIm;
  TH2F* BS_PtIs;
  TH3F* AS_PtIs;
  TH2F* BS_PtIm;
  TH3F* AS_PtIm;
  TH2F* BS_PtTOF;
  TH2F* BS_TOFIs;
  TH3F* AS_TOFIs;
  TH2F* BS_TOFIm;
  TH3F* AS_TOFIm;

  //Prediction histograms for muon only analysis which is binned depending on eta nd number of muon stations

  TH2F* genrecopT;
  TH1F* genlevelpT;
  TH1F* genleveleta;
  TH1F* genlevelbeta;
};

// initialize all the plots but also the directory structure to save them in the file
// WARNING: if you decide to add some histograms to the container, mind the binning of the histograms and keep in mind that we have a very large
// number of samples in our analysis... The size of the file can easilly explode
void stPlots_Init(TFile* HistoFile,
                  stPlots& st,
                  std::string BaseName,
                  unsigned int NCuts,
                  bool SkipSelectionPlot = false,
                  bool isSignal = true,
                  unsigned int NCuts_Flip = 0) {
  st.SelPlot = !SkipSelectionPlot;
  st.Name = BaseName;
  st.NCuts = NCuts;

  std::string Name;
  Name = BaseName;
  st.Directory = HistoFile->mkdir(Name.c_str(), Name.c_str());
  //return 0 if the directory already exist, in that case just take the directory
  if (!st.Directory)
    HistoFile->GetDirectory(Name.c_str());
  st.Directory->cd();

  if (SkipSelectionPlot)
    return;
  Name = "Gen_DecayLength";
  st.Gen_DecayLength = new TH1F(Name.c_str(), Name.c_str(), 1000, 0, 1000);
  st.Gen_DecayLength->Sumw2();
  Name = "Beta_Gen";
  st.Beta_Gen = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.Beta_Gen->Sumw2();
  Name = "Beta_GenChaged";
  st.Beta_GenCharged = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.Beta_GenCharged->Sumw2();
  Name = "Beta_Triggered";
  st.Beta_Triggered = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.Beta_Triggered->Sumw2();
  Name = "Beta_Matched";
  st.Beta_Matched = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.Beta_Matched->Sumw2();
  Name = "Beta_PreselectedA";
  st.Beta_PreselectedA = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.Beta_PreselectedA->Sumw2();
  Name = "Beta_PreselectedB";
  st.Beta_PreselectedB = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.Beta_PreselectedB->Sumw2();
  Name = "Beta_PreselectedC";
  st.Beta_PreselectedC = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.Beta_PreselectedC->Sumw2();
  Name = "Beta_SelectedP";
  st.Beta_SelectedP = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 20, 0, 1);
  st.Beta_SelectedP->Sumw2();
  Name = "Beta_SelectedI";
  st.Beta_SelectedI = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 20, 0, 1);
  st.Beta_SelectedI->Sumw2();
  Name = "Beta_SelectedT";
  st.Beta_SelectedT = new TH2F(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 20, 0, 1);
  st.Beta_SelectedT->Sumw2();

  Name = "BS_V3D";
  st.BS_V3D = new TH1F(Name.c_str(), Name.c_str(), 150, 0, IPbound);
  st.BS_V3D->Sumw2();
  Name = "BS_Chi2";
  st.BS_Chi2 = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 20);
  st.BS_Chi2->Sumw2();
  Name = "BS_Qual";
  st.BS_Qual = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 20);
  st.BS_Qual->Sumw2();
  Name = "BS_TNOH";
  st.BS_TNOH = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 40);
  st.BS_TNOH->Sumw2();
  Name = "BS_TNOH_PUA";
  st.BS_TNOH_PUA = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 40);
  st.BS_TNOH_PUA->Sumw2();
  Name = "BS_TNOH_PUB";
  st.BS_TNOH_PUB = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 40);
  st.BS_TNOH_PUB->Sumw2();
  Name = "BS_TNOHFraction";
  st.BS_TNOHFraction = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 1);
  st.BS_TNOHFraction->Sumw2();
  Name = "BS_TNOPH";
  st.BS_TNOPH = new TH1F(Name.c_str(), Name.c_str(), 16, 0, 8);
  st.BS_TNOPH->Sumw2();
  Name = "BS_TNOHFractionTillLast";
  st.BS_TNOHFractionTillLast = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 1);
  st.BS_TNOHFractionTillLast->Sumw2();
  Name = "BS_TNOMHTillLast";
  st.BS_TNOMHTillLast = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 20);
  st.BS_TNOMHTillLast->Sumw2();
  Name = "BS_Eta";
  st.BS_Eta = new TH1F(Name.c_str(), Name.c_str(), 50, -2.6, 2.6);
  st.BS_Eta->Sumw2();
  Name = "BS_TNOM";
  st.BS_TNOM = new TH1F(Name.c_str(), Name.c_str(), 40, 0, 40);
  st.BS_TNOM->Sumw2();
  Name = "BS_TNOM_PUA";
  st.BS_TNOM_PUA = new TH1F(Name.c_str(), Name.c_str(), 40, 0, 40);
  st.BS_TNOM_PUA->Sumw2();
  Name = "BS_TNOM_PUB";
  st.BS_TNOM_PUB = new TH1F(Name.c_str(), Name.c_str(), 40, 0, 40);
  st.BS_TNOM_PUB->Sumw2();
  Name = "BS_nDof";
  st.BS_nDof = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 40);
  st.BS_nDof->Sumw2();
  Name = "BS_TOFError";
  st.BS_TOFError = new TH1F(Name.c_str(), Name.c_str(), 25, 0, 0.25);
  st.BS_TOFError->Sumw2();
  Name = "BS_PtErr";
  st.BS_Pterr = new TH1F(Name.c_str(), Name.c_str(), 40, 0, 1);
  st.BS_Pterr->Sumw2();
  Name = "BS_MPt";
  st.BS_MPt = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_MPt->Sumw2();
  Name = "BS_MIs";
  st.BS_MIs = new TH1F(Name.c_str(), Name.c_str(), 50, 0, dEdxS_UpLim);
  st.BS_MIs->Sumw2();
  Name = "BS_MIm";
  st.BS_MIm = new TH1F(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  st.BS_MIm->Sumw2();
  Name = "BS_MTOF";
  st.BS_MTOF = new TH1F(Name.c_str(), Name.c_str(), 50, -2, 5);
  st.BS_MTOF->Sumw2();
  Name = "BS_TIsol";
  st.BS_TIsol = new TH1F(Name.c_str(), Name.c_str(), 25, 0, 100);
  st.BS_TIsol->Sumw2();
  Name = "BS_EIsol";
  st.BS_EIsol = new TH1F(Name.c_str(), Name.c_str(), 25, 0, 1.5);
  st.BS_EIsol->Sumw2();
  Name = "BS_SumpTOverpT";
  st.BS_SumpTOverpT = new TH1F(Name.c_str(), Name.c_str(), 80, 0.0, 2.0);
  st.BS_SumpTOverpT->Sumw2();
  Name = "BS_LastHitDXY";
  st.BS_LastHitDXY = new TH1F(Name.c_str(), Name.c_str(), 75, 0, 150);
  st.BS_LastHitDXY->Sumw2();
  Name = "BS_LastHitD3D";
  st.BS_LastHitD3D = new TH1F(Name.c_str(), Name.c_str(), 175, 0, 350);
  st.BS_LastHitD3D->Sumw2();
  Name = "BS_P";
  st.BS_P = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_P->Sumw2();
  Name = "BS_Pt";
  st.BS_Pt = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt->Sumw2();
  Name = "BS_Pt_PUA";
  st.BS_Pt_PUA = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_PUA->Sumw2();
  Name = "BS_Pt_PUB";
  st.BS_Pt_PUB = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_PUB->Sumw2();
  Name = "BS_Pt_Cosmic";
  st.BS_Pt_Cosmic = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_Cosmic->Sumw2();
  Name = "BS_Pt_DT";
  st.BS_Pt_DT = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_DT->Sumw2();
  Name = "BS_Pt_CSC";
  st.BS_Pt_CSC = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_CSC->Sumw2();
  Name = "BS_Is";
  st.BS_Is = new TH1F(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  st.BS_Is->Sumw2();
  Name = "BS_Is_PUA";
  st.BS_Is_PUA = new TH1F(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  st.BS_Is_PUA->Sumw2();
  Name = "BS_Is_PUB";
  st.BS_Is_PUB = new TH1F(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  st.BS_Is_PUB->Sumw2();
  Name = "BS_Is_Cosmic";
  st.BS_Is_Cosmic = new TH1F(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  st.BS_Is_Cosmic->Sumw2();
  Name = "BS_Im";
  st.BS_Im = new TH1F(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  st.BS_Im->Sumw2();
  Name = "BS_Im_PUA";
  st.BS_Im_PUA = new TH1F(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  st.BS_Im_PUA->Sumw2();
  Name = "BS_Im_PUB";
  st.BS_Im_PUB = new TH1F(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  st.BS_Im_PUB->Sumw2();
  Name = "BS_TOF";
  st.BS_TOF = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF->Sumw2();
  Name = "BS_TOF_PUA";
  st.BS_TOF_PUA = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF_PUA->Sumw2();
  Name = "BS_TOF_PUB";
  st.BS_TOF_PUB = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF_PUB->Sumw2();
  Name = "BS_TOF_DT";
  st.BS_TOF_DT = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF_DT->Sumw2();
  Name = "BS_TOF_CSC";
  st.BS_TOF_CSC = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF_CSC->Sumw2();
  Name = "BS_dR_NVTrack";
  st.BS_dR_NVTrack = new TH1F(Name.c_str(), Name.c_str(), 40, 0, 1);
  st.BS_dR_NVTrack->Sumw2();
  Name = "BS_MatchedStations";
  st.BS_MatchedStations = new TH1F(Name.c_str(), Name.c_str(), 8, -0.5, 7.5);
  st.BS_MatchedStations->Sumw2();
  Name = "BS_InnerInvPtDiff";
  st.BS_InnerInvPtDiff = new TH1F(Name.c_str(), Name.c_str(), 120, -4, 4);
  st.BS_InnerInvPtDiff->Sumw2();
  Name = "BS_Phi";
  st.BS_Phi = new TH1F(Name.c_str(), Name.c_str(), 50, -3.14, 3.14);
  st.BS_Phi->Sumw2();
  Name = "BS_TimeAtIP";
  st.BS_TimeAtIP = new TH1F(Name.c_str(), Name.c_str(), 50, -100, 100);
  st.BS_TimeAtIP->Sumw2();
  Name = "BS_OpenAngle";
  st.BS_OpenAngle = new TH1F(Name.c_str(), Name.c_str(), 50, -0.3, 3.15);
  st.BS_OpenAngle->Sumw2();
  Name = "BS_OpenAngle_Cosmic";
  st.BS_OpenAngle_Cosmic = new TH1F(Name.c_str(), Name.c_str(), 50, -0.3, 3.15);
  st.BS_OpenAngle_Cosmic->Sumw2();

  Name = "BS_NVertex";
  st.BS_NVertex = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 50);
  st.BS_NVertex->Sumw2();
  Name = "BS_NVertex_NoEventWeight";
  st.BS_NVertex_NoEventWeight = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 50);
  st.BS_NVertex_NoEventWeight->Sumw2();
  Name = "BS_PV";
  st.BS_PV = new TH1F(Name.c_str(), Name.c_str(), 60, 0, 60);
  st.BS_PV->Sumw2();
  Name = "BS_PV_NoEventWeight";
  st.BS_PV_NoEventWeight = new TH1F(Name.c_str(), Name.c_str(), 60, 0, 60);
  st.BS_PV_NoEventWeight->Sumw2();
  Name = "BS_NOMoNOHvsPV";
  st.BS_NOMoNOHvsPV = new TProfile(Name.c_str(), Name.c_str(), 60, 0, 60);
  st.BS_NOMoNOHvsPV->Sumw2();
  Name = "BS_dzAll";
  st.BS_dzAll = new TH1F(Name.c_str(), Name.c_str(), 200, -10, 10);
  st.BS_dzAll->Sumw2();
  Name = "BS_dxyAll";
  st.BS_dxyAll = new TH1F(Name.c_str(), Name.c_str(), 200, -10, 10);
  st.BS_dxyAll->Sumw2();
  Name = "BS_dzMinv3d";
  st.BS_dzMinv3d = new TH1F(Name.c_str(), Name.c_str(), 200, -10, 10);
  st.BS_dzMinv3d->Sumw2();
  Name = "BS_dxyMinv3d";
  st.BS_dxyMinv3d = new TH1F(Name.c_str(), Name.c_str(), 200, -10, 10);
  st.BS_dxyMinv3d->Sumw2();

  Name = "BS_SegSep";
  st.BS_SegSep = new TH1F(Name.c_str(), Name.c_str(), 50, 0, 2.5);
  st.BS_SegSep->Sumw2();
  Name = "BS_SegMinEtaSep";
  st.BS_SegMinEtaSep = new TH1F(Name.c_str(), Name.c_str(), 50, -1., 1.);
  st.BS_SegMinEtaSep->Sumw2();
  Name = "BS_SegMinPhiSep";
  st.BS_SegMinPhiSep = new TH1F(Name.c_str(), Name.c_str(), 50, -3.3, 3.3);
  st.BS_SegMinPhiSep->Sumw2();
  Name = "BS_SegMinEtaSep_FailDz";
  st.BS_SegMinEtaSep_FailDz = new TH1F(Name.c_str(), Name.c_str(), 50, -1., 1.);
  st.BS_SegMinEtaSep_FailDz->Sumw2();
  Name = "BS_SegMinEtaSep_PassDz";
  st.BS_SegMinEtaSep_PassDz = new TH1F(Name.c_str(), Name.c_str(), 50, -1., 1.);
  st.BS_SegMinEtaSep_PassDz->Sumw2();
  Name = "BS_Dz_FailSep";
  st.BS_Dz_FailSep = new TH1F(Name.c_str(), Name.c_str(), 50, -150, 150);
  st.BS_Dz_FailSep->Sumw2();

  Name = "BS_Dxy";
  st.BS_Dxy = new TH1F(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  st.BS_Dxy->Sumw2();
  Name = "BS_Dxy_Cosmic";
  st.BS_Dxy_Cosmic = new TH1F(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  st.BS_Dxy_Cosmic->Sumw2();
  Name = "BS_Dz";
  st.BS_Dz = new TH1F(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  st.BS_Dz->Sumw2();
  Name = "BS_Dz_Cosmic";
  st.BS_Dz_Cosmic = new TH1F(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  st.BS_Dz_Cosmic->Sumw2();
  Name = "BS_Dz_CSC";
  st.BS_Dz_CSC = new TH1F(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  st.BS_Dz_CSC->Sumw2();
  Name = "BS_Dz_DT";
  st.BS_Dz_DT = new TH1F(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  st.BS_Dz_DT->Sumw2();
  Name = "BS_Pt_FailDz";
  st.BS_Pt_FailDz = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_FailDz->Sumw2();
  Name = "BS_Pt_FailDz_DT";
  st.BS_Pt_FailDz_DT = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_FailDz_DT->Sumw2();
  Name = "BS_Pt_FailDz_CSC";
  st.BS_Pt_FailDz_CSC = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.BS_Pt_FailDz_CSC->Sumw2();
  Name = "BS_TOF_FailDz";
  st.BS_TOF_FailDz = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF_FailDz->Sumw2();
  Name = "BS_TOF_FailDz_DT";
  st.BS_TOF_FailDz_DT = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF_FailDz_DT->Sumw2();
  Name = "BS_TOF_FailDz_CSC";
  st.BS_TOF_FailDz_CSC = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
  st.BS_TOF_FailDz_CSC->Sumw2();
  Name = "genrecopT";
  st.genrecopT = new TH2F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 50, 0, PtHistoUpperBound);
  st.genrecopT->Sumw2();

  Name = "genlevelpT";
  st.genlevelpT = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  st.genlevelpT->Sumw2();
  Name = "genleveleta";
  st.genleveleta = new TH1F(Name.c_str(), Name.c_str(), 60, -3, 3);
  st.genleveleta->Sumw2();
  Name = "genlevelbeta";
  st.genlevelbeta = new TH1F(Name.c_str(), Name.c_str(), 20, 0, 1);
  st.genlevelbeta->Sumw2();

  //Initialize histograms for number of bins.  For everything but muon only PredBins=0 so no histograms created
  for (int i = 0; i < PredBins; i++) {
    char Suffix[1024];
    sprintf(Suffix, "_%i", i);
    Name = "BS_Pt_Binned";
    Name.append(Suffix);
    st.BS_Pt_Binned[i] = new TH1F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
    st.BS_Pt_Binned[i]->Sumw2();
    Name = "BS_TOF_Binned";
    Name.append(Suffix);
    st.BS_TOF_Binned[i] = new TH1F(Name.c_str(), Name.c_str(), 150, -1, 5);
    st.BS_TOF_Binned[i]->Sumw2();
  }

  Name = "BS_EtaIs";
  st.BS_EtaIs = new TH2F(Name.c_str(), Name.c_str(), 50, -3, 3, 50, 0, dEdxS_UpLim);
  Name = "BS_EtaIm";
  st.BS_EtaIm = new TH2F(Name.c_str(), Name.c_str(), 50, -3, 3, 100, 0, dEdxM_UpLim);
  Name = "BS_EtaP";
  st.BS_EtaP = new TH2F(Name.c_str(), Name.c_str(), 50, -3, 3, 50, 0, PtHistoUpperBound);
  Name = "BS_EtaPt";
  st.BS_EtaPt = new TH2F(Name.c_str(), Name.c_str(), 50, -3, 3, 50, 0, PtHistoUpperBound);
  Name = "BS_EtaTOF";
  st.BS_EtaTOF = new TH2F(Name.c_str(), Name.c_str(), 50, -3, 3, 50, 0, 3);
  Name = "BS_EtaDz";
  st.BS_EtaDz = new TH2F(Name.c_str(), Name.c_str(), 50, -3, 3, 50, -IPbound, IPbound);
  Name = "BS_PIs";
  st.BS_PIs = new TH2F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0, dEdxS_UpLim);
  Name = "BS_PImHD";
  st.BS_PImHD = new TH2F(Name.c_str(), Name.c_str(), 500, 0, PtHistoUpperBound, 1000, 0, dEdxM_UpLim);
  Name = "BS_PIm";
  st.BS_PIm = new TH2F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  Name = "BS_PtIs";
  st.BS_PtIs = new TH2F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 50, 0, dEdxS_UpLim);
  Name = "BS_PtIm";
  st.BS_PtIm = new TH2F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  Name = "BS_PtTOF";
  st.BS_PtTOF = new TH2F(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 50, 0, 5);
  //   Name = "BS_TOFIs"; st.BS_TOFIs = new TH2F(Name.c_str(), Name.c_str(),                   100, 1, 5, 100, 0, dEdxS_UpLim);
  Name = "BS_TOFIs";
  st.BS_TOFIs = new TH2F(Name.c_str(), Name.c_str(), 50, 0, 5, 50, 0, dEdxS_UpLim);
  //   Name = "BS_TOFIm"; st.BS_TOFIm = new TH2F(Name.c_str(), Name.c_str(),                   100, 1, 5, 200, 0, dEdxM_UpLim);
  Name = "BS_TOFIm";
  st.BS_TOFIm = new TH2F(Name.c_str(), Name.c_str(), 50, 0, 5, 100, 0, dEdxM_UpLim);

  Name = "H_D_DzSidebands";
  st.H_D_DzSidebands = new TH2D(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, DzRegions, 0, DzRegions);
  st.H_D_DzSidebands->Sumw2();

  //Background prediction histograms don't need to be made for signal or individual MC samples
}

HistoFile->cd();
}

// load all the plots from an already existing file
bool stPlots_InitFromFile(TFile* HistoFile, stPlots& st, std::string BaseName) {
  st.Name = BaseName;
  std::string Name;
  Name = BaseName;

  st.Directory = new TDirectory((Name + "ReadFromFile").c_str(), (Name + "ReadFromFile").c_str());
  st.Directory->cd();
  TDirectory::AddDirectory(kTRUE);
  TH1::AddDirectory(kTRUE);

  if (HistoFile->GetDirectory(BaseName.c_str()) == 0) {
    printf("Can't find subdirectory %s in opened file\n", BaseName.c_str());
    return false;
  }

  st.Gen_DecayLength = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Gen_DecayLength");
  st.Beta_Gen = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_Gen");
  st.Beta_GenCharged = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_GenCharged");
  st.Beta_Triggered = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_Triggered");
  st.Beta_Matched = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_Matched");
  st.Beta_PreselectedA = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_PreselectedA");
  st.Beta_PreselectedB = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_PreselectedB");
  st.Beta_PreselectedC = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_PreselectedC");
  st.Beta_SelectedP = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_SelectedP");
  st.Beta_SelectedI = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_SelectedI");
  st.Beta_SelectedT = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/Beta_SelectedT");

  st.BS_V3D = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_V3D");
  st.BS_Chi2 = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Chi2");
  st.BS_Qual = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Qual");
  st.BS_TNOH = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOH");
  st.BS_TNOH_PUA = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOH_PUA");
  st.BS_TNOH_PUB = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOH_PUB");
  st.BS_TNOHFraction = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOHFraction");
  st.BS_TNOPH = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOPH");
  st.BS_TNOHFractionTillLast = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOHFractionTillLast");
  st.BS_TNOMHTillLast = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOMHTillLast");
  st.BS_Eta = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Eta");
  st.BS_TNOM = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOM");
  st.BS_TNOM_PUA = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOM_PUA");
  st.BS_TNOM_PUB = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TNOM_PUB");
  st.BS_nDof = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_nDof");
  st.BS_TOFError = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOFError");
  st.BS_Pterr = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PtErr");
  st.BS_MPt = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_MPt");
  st.BS_MIm = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_MIm");
  st.BS_MIs = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_MIs");
  st.BS_MTOF = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_MTOF");
  st.BS_TIsol = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TIsol");
  st.BS_EIsol = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_EIsol");
  st.BS_SumpTOverpT = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_SumpTOverpT");
  st.BS_dR_NVTrack = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_dR_NVTrack");
  st.BS_MatchedStations = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_MatchedStations");
  st.BS_NVertex = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_NVertex");
  st.BS_NVertex_NoEventWeight =
      (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_NVertex_NoEventWeight");
  st.BS_PV = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PV");
  st.BS_NOMoNOHvsPV = (TProfile*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_NOMoNOHvsPV");
  st.BS_SegSep = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_SegSep");
  st.BS_SegMinPhiSep = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_SegMinPhiSep");
  st.BS_SegMinEtaSep = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_SegMinEtaSep");
  st.BS_SegMinEtaSep_FailDz = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_SegMinEtaSep_FailDz");
  st.BS_SegMinEtaSep_PassDz = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_SegMinEtaSep_PassDz");
  st.BS_Dz_FailSep = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Dz_FailSep");
  st.BS_InnerInvPtDiff = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_InnerInvPtDiff");
  st.BS_Phi = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Phi");
  st.BS_TimeAtIP = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TimeAtIP");
  st.BS_OpenAngle = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_OpenAngle");
  st.BS_OpenAngle_Cosmic = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_OpenAngle_Cosmic");

  st.BS_Pt_FailDz = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_FailDz");
  st.BS_Pt_FailDz_DT = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_FailDz_DT");
  st.BS_Pt_FailDz_CSC = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_FailDz_CSC");
  st.BS_TOF_FailDz = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_FailDz");
  st.BS_TOF_FailDz_DT = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_FailDz_DT");
  st.BS_TOF_FailDz_CSC = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_FailDz_CSC");
  st.BS_Dxy = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Dxy");
  st.BS_Dxy_Cosmic = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Dxy_Cosmic");
  st.BS_Dz = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Dz");
  st.BS_Dz_Cosmic = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Dz_Cosmic");
  st.BS_Dz_CSC = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Dz_CSC");
  st.BS_Dz_DT = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Dz_DT");
  st.genrecopT = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/genrecopT");
  st.genlevelpT = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/genlevelpT");
  st.genleveleta = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/genleveleta");
  st.genlevelbeta = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/genlevelbeta");
  st.BS_LastHitDXY = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_LastHitDXY");
  st.BS_LastHitD3D = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_LastHitD3D");
  st.BS_P = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_P");
  st.AS_P = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/AS_P");
  st.BS_Pt = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt");
  st.BS_Pt_PUA = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_PUA");
  st.BS_Pt_PUB = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_PUB");
  st.BS_Pt_Cosmic = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_Cosmic");
  st.BS_Pt_DT = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_DT");
  st.BS_Pt_CSC = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_CSC");
  st.AS_Pt = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/AS_Pt");
  st.BS_Im = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Im");
  st.BS_Im_PUA = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Im_PUA");
  st.BS_Im_PUB = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Im_PUB");
  st.AS_Im = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/AS_Im");
  st.BS_Is = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Is");
  st.BS_Is_PUA = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Is_PUA");
  st.BS_Is_PUB = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Is_PUB");
  st.BS_Is_Cosmic = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Is_Cosmic");
  st.AS_Is = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/AS_Is");
  st.BS_TOF = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF");
  st.BS_TOF_PUA = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_PUA");
  st.BS_TOF_PUB = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_PUB");
  st.BS_TOF_DT = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_DT");
  st.BS_TOF_CSC = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_CSC");
  st.AS_TOF = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/AS_TOF");
  st.BS_EtaIs = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_EtaIs");
  //st.AS_EtaIs  = (TH3F*)GetObjectFromPath(st.Directory, HistoFile,  BaseName + "/AS_EtaIs");
  st.BS_EtaIm = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_EtaIm");
  //st.AS_EtaIm  = (TH3F*)GetObjectFromPath(st.Directory, HistoFile,  BaseName + "/AS_EtaIm");
  st.BS_EtaP = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_EtaP");
  //st.AS_EtaP   = (TH3F*)GetObjectFromPath(st.Directory, HistoFile,  BaseName + "/AS_EtaP");
  st.BS_EtaPt = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_EtaPt");
  //st.AS_EtaPt  = (TH3F*)GetObjectFromPath(st.Directory, HistoFile,  BaseName + "/AS_EtaPt");
  st.BS_EtaTOF = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_EtaTOF");
  //st.AS_EtaTOF  = (TH3F*)GetObjectFromPath(st.Directory, HistoFile,  BaseName + "/AS_EtaTOF");
  st.BS_PIs = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PIs");
  st.BS_PImHD = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PImHD");
  st.BS_PIm = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PIm");
  st.BS_PtIs = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PtIs");
  st.BS_PtIm = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PtIm");
  st.BS_PtTOF = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_PtTOF");
  st.BS_TOFIs = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOFIs");
  st.BS_TOFIm = (TH2F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOFIm");

  st.BS_Pt_Binned[0] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_Binned_0");
  st.BS_Pt_Binned[1] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_Binned_1");
  st.BS_Pt_Binned[2] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_Binned_2");
  st.BS_Pt_Binned[3] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_Binned_3");
  st.BS_Pt_Binned[4] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_Binned_4");
  st.BS_Pt_Binned[5] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_Pt_Binned_5");

  st.BS_TOF_Binned[0] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_Binned_0");
  st.BS_TOF_Binned[1] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_Binned_1");
  st.BS_TOF_Binned[2] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_Binned_2");
  st.BS_TOF_Binned[3] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_Binned_3");
  st.BS_TOF_Binned[4] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_Binned_4");
  st.BS_TOF_Binned[5] = (TH1F*)GetObjectFromPath(st.Directory, HistoFile, BaseName + "/BS_TOF_Binned_5");

  HistoFile->cd();
  return true;
}

// Write the histograms to the file on disk and properly clean the memory from all the histograms
void stPlots_Clear(stPlots* st, bool WriteFirst = false) {
  if (WriteFirst) {
    st->Tree->SetDirectory(st->Directory);
    st->Directory->Write();
  }
  delete st->Directory;
}

// draw all plots that are not meant for comparison with other samples (mostly 2D plots that can't be superimposed)
void stPlots_Draw(stPlots& st, std::string SavePath, std::string LegendTitle, unsigned int CutIndex) {
  TypeMode = TypeFromPattern(SavePath);
  char YAxisTitle[2048];

  TObject** Histos = new TObject*[10];
  std::vector<std::string> legend;
  TCanvas* c1;

  char CutIndexStr[255];
  sprintf(CutIndexStr, "_%03i", CutIndex);

  char InputHist[1024];
  sprintf(InputHist, "Results/Type%i/Histos.root", TypeMode);
  TFile* InputFile = new TFile(InputHist);

  char PtCutStr[1024];
  sprintf(PtCutStr, "%.0f GeV", HCuts_Pt->GetBinContent(CutIndex + 1));
  char ICutStr[1024];
  sprintf(ICutStr, "%.2f", HCuts_Is->GetBinContent(CutIndex + 1));
  char TOFCutStr[1024];
  sprintf(TOFCutStr, "%.3f", HCuts_TOF->GetBinContent(CutIndex + 1));
  InputFile->Close();

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_EtaIs;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "#eta", dEdxS_Legend.c_str(), 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "EtaIs_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_EtaIm;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "#eta", dEdxM_Legend.c_str(), 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "EtaIm_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_EtaP;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "#eta", "p (GeV)", 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "EtaP_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_EtaPt;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "#eta", "p_{T} (GeV)", 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "EtaPt_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_EtaTOF;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "#eta", "1/#beta", 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "EtaTOF_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_PIs;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p (GeV)", dEdxS_Legend.c_str(), 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "PIs_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_PIm;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p (GeV)", dEdxM_Legend.c_str(), 0, 0, 0, 15, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "PIm_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_PtIs;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV)", dEdxS_Legend.c_str(), 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "PtIs_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_PtIm;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV)", dEdxM_Legend.c_str(), 0, 0, 0, 15, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "PtIm_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_PtTOF;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "p_{T} (GeV)", "1/#beta", 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "PtTOF_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_TOFIs;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "1/#beta", dEdxS_Legend.c_str(), 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TOFIs_BS", true);
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)st.BS_TOFIm;
  legend.push_back("Before Cut");
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "1/#beta", dEdxM_Legend.c_str(), 0, 0, 0, 15, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TOFIm_BS", true);
  delete c1;

  if (st.BS_Pt_PUA != NULL) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Pt_PUA->Clone();
    Histos1D[0]->Rebin(1);
    legend.push_back("NVtx<15");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_Pt_PUB->Clone();
    Histos1D[1]->Rebin(1);
    legend.push_back("NVtx>15");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    sprintf(YAxisTitle, "Fraction of tracks/%0.2f", ((TH1D*)Histos1D[0])->GetBinWidth(1));
    DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", "p_{T} (GeV)", YAxisTitle, 0, 1250, 1E-6, 1);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.20, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Pt_PU", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_TOF_PUA->Clone();
    Histos1D[0]->Rebin(1);
    legend.push_back("NVtx<15");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_TOF_PUB->Clone();
    Histos1D[1]->Rebin(1);
    legend.push_back("NVtx>15");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    sprintf(YAxisTitle, "Fraction of tracks/%0.2f", ((TH1D*)Histos1D[0])->GetBinWidth(1));
    DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", "1/#beta", YAxisTitle, 0, 4, 1E-6, 2);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.20, 0.045);  //,0.35);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "TOF_PU");
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Is_PUA->Clone();
    Histos1D[0]->Rebin(1);
    legend.push_back("NVtx<15");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_Is_PUB->Clone();
    Histos1D[1]->Rebin(1);
    legend.push_back("NVtx>15");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    sprintf(YAxisTitle, "Fraction of tracks/%0.2f", ((TH1D*)Histos1D[0])->GetBinWidth(1));
    DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", dEdxS_Legend.c_str(), YAxisTitle, 0, 0, 1E-6, 2);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.20, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Is_BS");
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Im_PUA->Clone();
    Histos1D[0]->Rebin(1);
    legend.push_back("NVtx<15");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_Im_PUB->Clone();
    Histos1D[1]->Rebin(1);
    legend.push_back("NVtx>15");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", dEdxM_Legend.c_str(), "Fraction of tracks", 0, 20, 1E-6, 2);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.20, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Im_BS");
    delete c1;
  }

  if (st.Name.find("Cosmic") != string::npos) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Pt_FailDz->Clone();
    Histos1D[0]->Rebin(1);
    legend.push_back("abs(dz)>35");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_Pt->Clone();
    Histos1D[1]->Rebin(1);
    legend.push_back("abs(dz)<35");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    DrawSuperposedHistos(
        (TH1**)Histos1D, legend, "E1", "p_{T} (GeV)", "arbitrary units", 0, 600, 0.0005, 1, false, false, true, false);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.2, 0.1);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "_Pt_Dz_Comp", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Pt_FailDz_CSC->Clone();
    Histos1D[0]->Rebin(1);
    legend.push_back("abs(dz)>35");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_Pt_CSC->Clone();
    Histos1D[1]->Rebin(1);
    legend.push_back("abs(dz)<35");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    DrawSuperposedHistos(
        (TH1**)Histos1D, legend, "E1", "p_{T} (GeV)", "arbitrary units", 0, 600, 0.0005, 1, false, false, true, false);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.2, 0.1);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "_Pt_Dz_Comp_CSC", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Pt_FailDz_DT->Clone();
    Histos1D[0]->Rebin(1);
    legend.push_back("abs(dz)>35");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_Pt_DT->Clone();
    Histos1D[1]->Rebin(1);
    legend.push_back("abs(dz)<35");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    DrawSuperposedHistos(
        (TH1**)Histos1D, legend, "E1", "p_{T} (GeV)", "arbitrary units", 0, 600, 0.0005, 1, false, false, true, false);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.2, 0.1);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "_Pt_Dz_Comp_DT", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_TOF_FailDz->Clone();
    legend.push_back("abs(z)>35");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_TOF->Clone();
    legend.push_back("abs(dz)<35");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    DrawSuperposedHistos(
        (TH1**)Histos1D, legend, "E1", "1/#beta", "arbitrary units", -2, 4, 0.0005, 1, false, false, true, false);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.2, 0.1);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "_TOF_Dz_Comp", true);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_TOF_FailDz_CSC->Clone();
    legend.push_back("abs(z)>35");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_TOF_CSC->Clone();
    legend.push_back("abs(dz)<35");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    DrawSuperposedHistos(
        (TH1**)Histos1D, legend, "E1", "1/#beta", "arbitrary units", -2, 4, 0.0005, 1, false, false, true, false);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.2, 0.1);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "_TOF_Dz_CSC_Comp", true);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_TOF_FailDz_DT->Clone();
    legend.push_back("abs(z)>35");
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[1] = (TH1*)st.BS_TOF_DT->Clone();
    legend.push_back("abs(dz)<35");
    if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
      Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
    DrawSuperposedHistos(
        (TH1**)Histos1D, legend, "E1", "1/#beta", "arbitrary units", -2, 4, 0.0005, 1, false, false, true, false);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.2, 0.1);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "_TOF_Dz_DT_Comp", true);
    delete c1;
  }

  if (TypeMode == 3 && st.Name.find("Data") != string::npos) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Pt_Binned[0]->Clone();
    legend.push_back("Bar - 2 Sta");
    Histos1D[0]->Rebin(2);
    if (Histos1D[0]->Integral() > 0)
      Histos1D[0]->Scale(1. / Histos1D[0]->Integral());
    Histos1D[1] = (TH1*)st.BS_Pt_Binned[1]->Clone();
    legend.push_back("Bar - 3 Sta");
    Histos1D[1]->Rebin(2);
    if (Histos1D[1]->Integral() > 0)
      Histos1D[1]->Scale(1. / Histos1D[1]->Integral());
    Histos1D[2] = (TH1*)st.BS_Pt_Binned[2]->Clone();
    legend.push_back("Bar - 4 Sta");
    Histos1D[2]->Rebin(2);
    if (Histos1D[2]->Integral() > 0)
      Histos1D[2]->Scale(1. / Histos1D[2]->Integral());
    Histos1D[3] = (TH1*)st.BS_Pt_Binned[3]->Clone();
    legend.push_back("For - 2 Sta");
    Histos1D[3]->Rebin(2);
    if (Histos1D[3]->Integral() > 0)
      Histos1D[3]->Scale(1. / Histos1D[3]->Integral());
    Histos1D[4] = (TH1*)st.BS_Pt_Binned[4]->Clone();
    legend.push_back("For - 3 Sta");
    Histos1D[4]->Rebin(2);
    if (Histos1D[4]->Integral() > 0)
      Histos1D[4]->Scale(1. / Histos1D[4]->Integral());
    Histos1D[5] = (TH1*)st.BS_Pt_Binned[5]->Clone();
    legend.push_back("For - 4 Sta");
    Histos1D[5]->Rebin(2);
    if (Histos1D[5]->Integral() > 0)
      Histos1D[5]->Scale(1. / Histos1D[5]->Integral());
    DrawSuperposedHistos((TH1**)Histos1D, legend, "COLZ", "p_T", "arbitrary units", 0, 0, 0, 0, false);
    //DrawLegend((TObject**)Histos1D,legend,"","P", 0.79, 0.92, 0.2, 0.1);
    DrawLegend((TObject**)Histos1D, legend, "", "P");
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Pt_Binned_BS", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_TOF_Binned[0]->Clone();
    legend.push_back("Bar - 2 Sta");
    Histos1D[0]->Rebin(2);
    if (Histos1D[0]->Integral() > 0)
      Histos1D[0]->Scale(1. / Histos1D[0]->Integral());
    Histos1D[1] = (TH1*)st.BS_TOF_Binned[1]->Clone();
    legend.push_back("Bar - 3 Sta");
    Histos1D[1]->Rebin(2);
    if (Histos1D[1]->Integral() > 0)
      Histos1D[1]->Scale(1. / Histos1D[1]->Integral());
    Histos1D[2] = (TH1*)st.BS_TOF_Binned[2]->Clone();
    legend.push_back("Bar - 4 Sta");
    Histos1D[2]->Rebin(2);
    if (Histos1D[2]->Integral() > 0)
      Histos1D[2]->Scale(1. / Histos1D[2]->Integral());
    Histos1D[3] = (TH1*)st.BS_TOF_Binned[3]->Clone();
    legend.push_back("For - 2 Sta");
    Histos1D[3]->Rebin(2);
    if (Histos1D[3]->Integral() > 0)
      Histos1D[3]->Scale(1. / Histos1D[3]->Integral());
    Histos1D[4] = (TH1*)st.BS_TOF_Binned[4]->Clone();
    legend.push_back("For - 3 Sta");
    Histos1D[4]->Rebin(2);
    if (Histos1D[4]->Integral() > 0)
      Histos1D[4]->Scale(1. / Histos1D[4]->Integral());
    Histos1D[5] = (TH1*)st.BS_TOF_Binned[5]->Clone();
    legend.push_back("For - 4 Sta");
    Histos1D[5]->Rebin(2);
    if (Histos1D[5]->Integral() > 0)
      Histos1D[5]->Scale(1. / Histos1D[5]->Integral());
    DrawSuperposedHistos((TH1**)Histos1D, legend, "COLZ", "1/#beta", "arbitrary units", 0, 2, 0, 0, false);
    DrawLegend((TObject**)Histos1D, legend, "", "P");
    //DrawLegend((TObject**)Histos1D,legend,"","P", 0.79, 0.92, 0.2, 0.1);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "TOF_Binned_BS", false);
    delete c1;
  }

  if (TypeMode == 5) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos[0] = (TH1*)st.BS_OpenAngle_Cosmic->Clone();
    legend.push_back("|Dz|>0.5cm & |Dxy|>0.5cm");
    ((TH1D*)Histos[0])->Rebin(2);
    DrawSuperposedHistos(
        (TH1**)Histos, legend, "HIST", "#theta max", "Number of tracks", -0.5, 3.2, 0, 0, false, false, true, true);
    DrawLegend((TObject**)Histos, legend, "", "L", 0.93, 0.88, 0.38, 0.045);
    //   c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "OpenAngle_Cosmic_BS", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos[0] = (TH1*)st.BS_Dxy_Cosmic->Clone();
    legend.push_back("|Dz|>0.5cm & #theta>2.8rad");
    ((TH1D*)Histos[0])->Rebin(4);
    //   sprintf(YAxisTitle,"Fraction of tracks/%2.0f (cm)",Histos[0]->GetBinWidth(1));
    DrawSuperposedHistos(
        (TH1**)Histos, legend, "HIST", "Dxy (cm)", "Number of tracks", -2, 2, 0, 0, false, false, true, true);
    DrawLegend((TObject**)Histos, legend, "", "L", 0.93, 0.88, 0.38, 0.045);
    //   c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Dxy_Cosmic_BS", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos[0] = (TH1*)st.BS_Dz_Cosmic->Clone();
    legend.push_back("|Dxy|>0.5cm & #theta>2.8rad");
    ((TH1D*)Histos[0])->Rebin(4);
    //   sprintf(YAxisTitle,"Fraction of tracks/%2.0f (cm)",Histos[0]->GetBinWidth(1));
    DrawSuperposedHistos(
        (TH1**)Histos, legend, "HIST", "Dz (cm)", "Number of tracks", 0, 0, 0, 0, false, false, true, true);
    DrawLegend((TObject**)Histos, legend, "", "L", 0.93, 0.88, 0.38, 0.045);
    //   c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Dz_Cosmic_BS", false);
    delete c1;
  }

  if (TypeMode == 3) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Dz->Clone();
    legend.push_back("Dz");
    Histos1D[0]->Rebin(1);
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[0]->Fit("gaus");
    sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos1D[0]->GetBinWidth(1));
    DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", "Dz (cm)", YAxisTitle, 0, 0, 5E-4, 2, false, false, true, true);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.38, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "DzFit_BS", false);
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos1D[0] = (TH1*)st.BS_Dxy->Clone();
    legend.push_back("Dxy");
    Histos1D[0]->Rebin(1);
    if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
      Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
    Histos1D[0]->Fit("gaus");
    sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos1D[0]->GetBinWidth(1));
    DrawSuperposedHistos(
        (TH1**)Histos1D, legend, "E1", "Dxy (cm)", YAxisTitle, 0, 0, 5E-4, 2, false, false, true, true);
    DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.38, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "DxyFit_BS", false);
    delete c1;
  }
}

// draw all plots that meant for comparison with other samples (mostly 1D plots that can be superimposed)
void stPlots_DrawComparison(std::string SavePath,
                            std::string LegendTitle,
                            unsigned int CutIndex,
                            unsigned int CutIndexTight,
                            stPlots* st1,
                            stPlots* st2 = NULL,
                            stPlots* st3 = NULL,
                            stPlots* st4 = NULL,
                            stPlots* st5 = NULL,
                            stPlots* st6 = NULL,
                            stPlots* st7 = NULL) {
  char CutIndexStr[255];
  sprintf(CutIndexStr, "_%03i", CutIndex);

  //bool IsTkOnly = (SavePath.find("Type0",0)<std::string::npos);
  TypeMode = TypeFromPattern(SavePath);
  char YAxisTitle[2048];

  std::vector<std::string> lg;
  std::vector<stPlots*> st;
  if (st1)
    st.push_back(st1);
  if (st2)
    st.push_back(st2);
  if (st3)
    st.push_back(st3);
  if (st4)
    st.push_back(st4);
  if (st5)
    st.push_back(st5);
  if (st6)
    st.push_back(st6);
  if (st7)
    st.push_back(st7);

  Color[2] = 2;

  std::vector<stSample> samples;
  GetSampleDefinition(samples);
  for (unsigned int i = 0; i < st.size(); i++) {
    int Index = -1;
    for (unsigned int s = 0; s < samples.size(); s++) {
      if (samples[s].Name == st[i]->Name) {
        Index = s;
        break;
      }
    }
    if (st[i]->Name.find("MCTr") != string::npos) {
      lg.push_back("MC (SM)");
    } else if (st[i]->Name.find("Data7TeV") != string::npos) {
      lg.push_back("Data #sqrt{s} = 7.0 TeV");
    } else if (st[i]->Name.find("Data8TeV") != string::npos) {
      lg.push_back("Data #sqrt{s} = 8.0 TeV");
    } else if (st[i]->Name.find("Cosmic") != string::npos) {
      lg.push_back("Cosmic");
      Color[i] = 7;
    } else if (Index == -1) {
      lg.push_back(st[i]->Name);
    } else {
      lg.push_back(samples[Index].Legend);
    }
  }

  TH1** Histos = new TH1*[10];
  TH1** Histos1D = new TH1*[10];
  std::vector<std::string> legend;
  TCanvas* c1;

  for (unsigned int i = 0; i < st.size(); i++) {
    //      if(st[i]->Name=="Data")continue;
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    Histos[0] = (TH1*)st[i]->Beta_Gen;
    legend.push_back("Gen");
    //      Histos[1] = (TH1*)st[i]->Beta_GenCharged;                                           legend.push_back("Charged Gen");
    Histos[1] = (TH1*)st[i]->Beta_Triggered;
    legend.push_back("Triggered");
    DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1", "#beta", "# HSCP", 0, 0, -10, -10);
    DrawLegend((TObject**)Histos, legend, "", "P", 0.65, 0.88, 0.20, 0.04);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, st[i]->Name + "_GenBeta", true);
    delete c1;
  }

  for (unsigned int C = 0; C < 2; C++) {
    unsigned int CutIndex_ = C == 0 ? CutIndex : CutIndexTight;

    for (unsigned int i = 0; i < st.size(); i++) {
      //      if(st[i]->Name=="Data")continue;
      c1 = new TCanvas("c1", "c1,", 600, 600);
      legend.clear();
      Histos[0] = (TH1*)st[i]->Beta_Gen;
      legend.push_back("Gen");
      //    Histos[1] = (TH1*)st[i]->Beta_GenCharged;                                           legend.push_back("Charged Gen");
      Histos[1] = (TH1*)st[i]->Beta_Triggered;
      legend.push_back("Triggered");
      Histos[2] = (TH1*)st[i]->Beta_Matched;
      legend.push_back("Reconstructed");
      //    Histos[0] = (TH1*)st[i]->Beta_PreselectedA;                                         legend.push_back("PreselectedA");
      //    Histos[0] = (TH1*)st[i]->Beta_PreselectedB;                                         legend.push_back("PreselectedB");
      Histos[3] = (TH1*)st[i]->Beta_PreselectedC;
      legend.push_back("Preselected");
      Histos[4] = (TH1*)st[i]->Beta_SelectedP->ProjectionY("A", CutIndex_ + 1, CutIndex_ + 1);
      legend.push_back("p_{T}>Cut");
      Histos[5] = (TH1*)st[i]->Beta_SelectedI->ProjectionY("B", CutIndex_ + 1, CutIndex_ + 1);
      legend.push_back("I  >Cut");
      if (!(TypeMode == 0 || TypeMode == 5)) {
        Histos[6] = (TH1*)st[i]->Beta_SelectedT->ProjectionY("C", CutIndex_ + 1, CutIndex_ + 1);
        legend.push_back("ToF>Cut");
      }
      DrawSuperposedHistos((TH1**)Histos, legend, "HIST E1", "#beta", "# HSCP", 0, 0, -10, -10);
      DrawLegend((TObject**)Histos, legend, "", "P", 0.65, 0.88, 0.20, 0.025);
      c1->SetLogy(true);
      DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
      if (C == 0)
        SaveCanvas(c1, SavePath, st[i]->Name + "_Beta");
      else
        SaveCanvas(c1, SavePath, st[i]->Name + "_BetaTight");
      //for(int l=0;l<legend.size();l++){delete Histos[l];}
      delete c1;
    }
  }

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->Gen_DecayLength->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Decay Length (cm)", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Gen_DecayLength", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  for (unsigned int i = 0; i < st.size(); i++) {
    if (st[i]->BS_Pt_PUA != NULL) {
      c1 = new TCanvas("c1", "c1,", 600, 600);
      legend.clear();
      Histos1D[0] = (TH1*)st[i]->BS_Pt_PUA->Clone();
      Histos1D[0]->Rebin(1);
      legend.push_back("NVtx<15");
      if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
        Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
      Histos1D[1] = (TH1*)st[i]->BS_Pt_PUB->Clone();
      Histos1D[1]->Rebin(1);
      legend.push_back("NVtx>15");
      if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
        Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
      sprintf(YAxisTitle, "Fraction of tracks/%0.2f", ((TH1D*)Histos1D[0])->GetBinWidth(1));
      DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", "p_{T} (GeV)", YAxisTitle, 0, 1250, 1E-6, 3);
      DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.25, 0.045);
      c1->SetLogy(true);
      DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1, SavePath, st[i]->Name + "_Pt_PU", false);
      delete c1;

      c1 = new TCanvas("c1", "c1,", 600, 600);
      legend.clear();
      Histos1D[0] = (TH1*)st[i]->BS_TOF_PUA->Clone();
      Histos1D[0]->Rebin(1);
      legend.push_back("NVtx<15");
      if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
        Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
      Histos1D[1] = (TH1*)st[i]->BS_TOF_PUB->Clone();
      Histos1D[1]->Rebin(1);
      legend.push_back("NVtx>15");
      if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
        Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
      sprintf(YAxisTitle, "Fraction of tracks/%0.2f", ((TH1D*)Histos1D[0])->GetBinWidth(1));
      DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", "1/#beta", YAxisTitle, 0, 4, 1E-6, 3);
      DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.25, 0.045);  //,0.35);
      c1->SetLogy(true);
      DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1, SavePath, st[i]->Name + "_TOF_PU");
      delete c1;

      c1 = new TCanvas("c1", "c1,", 600, 600);
      legend.clear();
      Histos1D[0] = (TH1*)st[i]->BS_Is_PUA->Clone();
      Histos1D[0]->Rebin(1);
      legend.push_back("NVtx<15");
      if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
        Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
      Histos1D[1] = (TH1*)st[i]->BS_Is_PUB->Clone();
      Histos1D[1]->Rebin(1);
      legend.push_back("NVtx>15");
      if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
        Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
      sprintf(YAxisTitle, "Fraction of tracks/%0.2f", ((TH1D*)Histos1D[0])->GetBinWidth(1));
      DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", dEdxS_Legend.c_str(), YAxisTitle, 0, 0, 1E-6, 3);
      DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.25, 0.045);
      c1->SetLogy(true);
      DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1, SavePath, st[i]->Name + "_Is_PU");
      delete c1;

      c1 = new TCanvas("c1", "c1,", 600, 600);
      legend.clear();
      Histos1D[0] = (TH1*)st[i]->BS_Im_PUA->Clone();
      Histos1D[0]->Rebin(1);
      legend.push_back("NVtx<15");
      if (Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1) > 0)
        Histos1D[0]->Scale(1.0 / Histos1D[0]->Integral(0, Histos1D[0]->GetNbinsX() + 1));
      Histos1D[1] = (TH1*)st[i]->BS_Im_PUB->Clone();
      Histos1D[1]->Rebin(1);
      legend.push_back("NVtx>15");
      if (Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1) > 0)
        Histos1D[1]->Scale(1.0 / Histos1D[1]->Integral(0, Histos1D[1]->GetNbinsX() + 1));
      DrawSuperposedHistos((TH1**)Histos1D, legend, "E1", dEdxM_Legend.c_str(), "Fraction of tracks", 0, 20, 1E-6, 3);
      DrawLegend((TObject**)Histos1D, legend, "", "P", 0.93, 0.88, 0.25, 0.045);
      c1->SetLogy(true);
      DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1, SavePath, st[i]->Name + "_Im_PU");
      delete c1;
    }
  }

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_V3D->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "V3D (cm)", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "V3D_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Chi2->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "#chi^{2}/ndof", "Fraction of tracks", 0, 0, 1E-3, 3, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Chi2_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Qual->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "quality", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Quality_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TNOH->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#NOH", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NOH_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TNOHFraction->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Fraction of hits", "Fraction of tracks", 0, 0, 1E-3, 3);
  //   DrawLegend((TObject**)Histos,legend,"","P",0.62, 0.90, 0.38, 0.05);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NOHFraction_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TNOPH->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#NOPH", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NOPH_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TNOHFractionTillLast->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "Fraction of hits (till last)", "Fraction of tracks", 0, 0, 1E-3, 3);
  //   DrawLegend((TObject**)Histos,legend,"","P",0.62, 0.90, 0.38, 0.05);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NOHFractionTillLast_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TNOMHTillLast->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#NOH missing (till last)", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NOMHTillLast_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Eta->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#eta", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  //c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Eta_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TNOM->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#NOM", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NOM_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_NOMoNOHvsPV->Clone();
    legend.push_back(lg[i]);
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#Vertices", "#NOM/#NOH", 0, 30, 0, 1);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  //   c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NOMoNOHvsPV_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_nDof->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "TOF_{nDof}", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "nDof_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TOFError->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "1/#beta error", "Fraction of tracks", 0, 0, 1E-3, 3, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TOFError_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Pterr->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "p_{T} Err / p_{T}", "Fraction of tracks", 0, 0, 1E-3, 3, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Pterr_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_MPt->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "p_{T} (GeV)", "Fraction of tracks", 0, 1250, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "MPt_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_MIs->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", dEdxS_Legend.c_str(), "Fraction of tracks", 0, 0, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "MIs_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_MIm->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", dEdxM_Legend.c_str(), "Fraction of tracks", 0, 20, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "MIm_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_MTOF->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "1/#beta", "Fraction of tracks", -2, 5, 1E-6, 1);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "MTOF_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TIsol->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "Isolation: Track SumPt (GeV)", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "IsolT_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_EIsol->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "Isolation: (Ecal + Hcal) Energy / p", "Fraction of tracks", 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "IsolE_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_dR_NVTrack->Clone();
    Histos[i]->Rebin(1);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "dR", "Fraction of tracks", 0, 0.4, 1E-3, 3, false, false, true, false);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "dR_NVTrack_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_SegSep->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos,
                       legend,
                       "E1",
                       "dR to opp side segment",
                       "Fraction of tracks",
                       0,
                       2.5,
                       1E-3,
                       3,
                       false,
                       false,
                       true,
                       false);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "SegSep_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_SegMinPhiSep->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos,
                       legend,
                       "E1",
                       "dPhi to opp side segment",
                       "Fraction of tracks",
                       0,
                       0,
                       1E-3,
                       3,
                       false,
                       false,
                       true,
                       false);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "SegMinPhiSep_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos,
                       legend,
                       "E1",
                       "dEta to opp side segment",
                       "Fraction of tracks",
                       0,
                       0,
                       1E-3,
                       3,
                       false,
                       false,
                       true,
                       true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.6, 0.92, 0.38, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "SegMinEtaSep_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep_FailDz->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos,
                       legend,
                       "E1",
                       "dR to opp side segment",
                       "Fraction of tracks",
                       -0.5,
                       0.5,
                       1E-3,
                       3,
                       false,
                       false,
                       true,
                       true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "SegMinEtaSep_FailDz_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_SegMinEtaSep_PassDz->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos,
                       legend,
                       "E1",
                       "dEta to opp side segment",
                       "Fraction of tracks",
                       -1.5,
                       1.5,
                       1E-3,
                       3,
                       false,
                       false,
                       true,
                       true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "SegMinEtaSep_PassDz_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_MatchedStations->Clone();
    legend.push_back(lg[i]);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Muon stations", "Fraction of tracks", 0, 5, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "MatchedStations_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_PV->Clone();
    Histos[i]->Rebin(1);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Primary Vertices", "Fraction of events", 0, 0, 0, 0.3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(false);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "PV_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_InnerInvPtDiff->Clone();
    Histos[i]->Rebin(4);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Inner Inverse Pt Diff", "Fraction of tracks", -3, 3, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "InnerInvPtDiff_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Phi->Clone();
    Histos[i]->Rebin(2);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "#phi", "Fraction of tracks", -3.14, 3.14, -1.2, -1.2, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(false);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Phi_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TimeAtIP->Clone();
    Histos[i]->Rebin(1);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Time At Vertex (ns)", "Fraction of tracks", 0, 0, 1E-4, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TimeAtIP_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  if (TypeMode == 5) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    for (unsigned int i = 0; i < st.size(); i++) {
      Histos[i] = (TH1*)st[i]->BS_OpenAngle->Clone();
      Histos[i]->Rebin(1);
      legend.push_back(lg[i]);
      if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
        Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
    }
    DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#theta max", "Fraction of tracks", 0, 0, 1E-4, 3);
    DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "OpenAngle_BS", false);
    for (unsigned int i = 0; i < st.size(); i++) {
      delete Histos[i];
    }
    delete c1;

    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    for (unsigned int i = 0; i < st.size(); i++) {
      Histos[i] = (TH1*)st[i]->BS_OpenAngle_Cosmic->Clone();
      Histos[i]->Rebin(1);
      legend.push_back(lg[i]);
      if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
        Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
    }
    DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#theta max", "Fraction of tracks", 0, 0, 1E-4, 3);
    DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "OpenAngle_Cosmic_BS", false);
    for (unsigned int i = 0; i < st.size(); i++) {
      delete Histos[i];
    }
    delete c1;
  }

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Dz_FailSep->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Dz (cm)", "Fraction of tracks", 0, 0, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Dz_FailSep_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Dxy->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "Dxy (cm)", YAxisTitle, -0.5, 0.5, 5E-4, 3, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Dxy_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  if (TypeMode == 5) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    for (unsigned int i = 0; i < st.size(); i++) {
      Histos[i] = (TH1*)st[i]->BS_Dxy_Cosmic->Clone();
      legend.push_back(lg[i]);
      Histos[i]->Rebin(1);
      if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
        Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
    }
    sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos[0]->GetBinWidth(1));
    DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Dxy (cm)", YAxisTitle, 0, 0, 1E-3, 3, false, false, true, true);
    DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Dxy_Cosmic_BS", false);
    for (unsigned int i = 0; i < st.size(); i++) {
      delete Histos[i];
    }
    delete c1;
  }

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Dz->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos(
      (TH1**)Histos, legend, "E1", "Dz (cm)", YAxisTitle, -0.5, 0.5, 5E-4, 3, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Dz_BS", false);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  if (TypeMode == 5) {
    c1 = new TCanvas("c1", "c1,", 600, 600);
    legend.clear();
    for (unsigned int i = 0; i < st.size(); i++) {
      Histos[i] = (TH1*)st[i]->BS_Dz_Cosmic->Clone();
      legend.push_back(lg[i]);
      Histos[i]->Rebin(1);
      if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
        Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
    }
    sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos[0]->GetBinWidth(1));
    DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Dz (cm)", YAxisTitle, 0, 0, 1E-3, 3, false, false, true, true);
    DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
    c1->SetLogy(true);
    DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
    SaveCanvas(c1, SavePath, "Dz_Cosmic_BS", false);
    for (unsigned int i = 0; i < st.size(); i++) {
      delete Histos[i];
    }
    delete c1;
  }

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Dz_CSC->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Dz (cm)", YAxisTitle, 0, 0, 1E-3, 3, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Dz_CSC_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Dz_DT->Clone();
    legend.push_back(lg[i]);
    Histos[i]->Rebin(1);
    if (Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1) > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral(0, Histos[i]->GetNbinsX() + 1));
  }
  sprintf(YAxisTitle, "Fraction of tracks/%2.0f (cm)", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Dz (cm)", YAxisTitle, 0, 0, 1E-3, 3, false, false, true, true);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Dz_DT_BS", true);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_LastHitDXY;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "last hit DXY", YAxisTitle, 0, 0, 1E-6, 15);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.92, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "LastHitDXY_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_LastHitDXY;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "last hit D3D", YAxisTitle, 0, 0, 1E-6, 15);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.92, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "LastHitD3D_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Is;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", dEdxS_Legend.c_str(), YAxisTitle, 0, 0, 1E-6, 15);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.92, 0.45, 0.045);
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", dEdxS_Legend.c_str(), YAxisTitle, 0, 1, 1E-6, 15, false, true);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Is_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Im;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", dEdxM_Legend.c_str(), "Fraction of tracks", 0, 20, 1E-6, 15);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.92, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Im_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)(st[i]->AS_Is->ProjectionY((st[i]->Name + "AA").c_str(), CutIndex + 1, CutIndex + 1));
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", dEdxS_Legend.c_str(), "Fraction of tracks", 0, 0, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, std::string("Is_AS") + CutIndexStr);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->AS_Im->ProjectionY((st[i]->Name + "BB").c_str(), CutIndex + 1, CutIndex + 1);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", dEdxM_Legend.c_str(), "Fraction of tracks", 0, 20, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, std::string("Im_AS") + CutIndexStr);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Pt;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%2.0f GeV", Histos[0]->GetBinWidth(1));
  //DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} (GeV)", YAxisTitle, 0,1250, 0.000000001, 1.2);
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "p_{T} (GeV)", YAxisTitle, 0, 1250, 1E-6, 15);
  //if(IsTkOnly) DrawLegend((TObject**)Histos,legend,"","P", 0.45, 0.42, 0.26, 0.05);
  //else DrawLegend((TObject**)Histos,legend,"","P", 0.51, 0.39, 0.33, 0.05);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.92, 0.45, 0.045);
  //DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "p_{T} (GeV)", YAxisTitle, 0,1250, 0.000000001, 1.2, false, true);

  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Pt_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_Pt_FailDz;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%2.0f GeV", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "p_{T} (GeV)", YAxisTitle, 0, 1250, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "Pt_FailDz_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->AS_Pt->ProjectionY((st[i]->Name + "CC").c_str(), CutIndex + 1, CutIndex + 1);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "p_{T} (GeV)", "Fraction of tracks", 0, 1250, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, std::string("Pt_AS") + CutIndexStr);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TOF;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "1/#beta", YAxisTitle, 0, 4, 1E-6, 15);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.92, 0.45, 0.045);  //,0.35);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TOF_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TOF_DT;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "1/#beta", YAxisTitle, -1, 4, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TOF_DT_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TOF_CSC;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "1/#beta", YAxisTitle, -1, 4, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);  //,0.35);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TOF_CSC_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->AS_TOF->ProjectionY((st[i]->Name + "DD").c_str(), CutIndex + 1, CutIndex + 1);
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "1/#beta", "Fraction of tracks", 1, 4, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);  //, 0.35);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, std::string("TOF_AS") + CutIndexStr);
  for (unsigned int i = 0; i < st.size(); i++) {
    delete Histos[i];
  }
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_TOF_FailDz;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "1/#beta", YAxisTitle, 0, 4, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);  //,0.35);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "TOF_FailDz_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_NVertex;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.0f vertex", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "Number of reconstructed vertices", YAxisTitle, 0, 0, 1E-3, 0.30);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "NVertex_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->BS_SumpTOverpT;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.3f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "#sump_{T}/p_{T}", YAxisTitle, 0, 0, 1E-3, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "SumptOverpt_BS");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->genlevelpT;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "gen-p_{T}", YAxisTitle, 0, 1200, 1E-6, 3);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);  //,0.35);
  c1->SetLogy(true);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "genpT");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->genleveleta;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "gen-#eta", YAxisTitle, -3, 3, 0.01, 0.040);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);  //,0.35);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "geneta");
  delete c1;

  c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  for (unsigned int i = 0; i < st.size(); i++) {
    Histos[i] = (TH1*)st[i]->genlevelbeta;
    legend.push_back(lg[i]);
    if (Histos[i]->Integral() > 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
  }
  sprintf(YAxisTitle, "Fraction of tracks/%0.2f", Histos[0]->GetBinWidth(1));
  DrawSuperposedHistos((TH1**)Histos, legend, "E1", "gen-#beta", YAxisTitle, 0, 1, 0.0001, 0.5);
  DrawLegend((TObject**)Histos, legend, "", "P", 0.93, 0.88, 0.45, 0.045);  //,0.35);
  DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  SaveCanvas(c1, SavePath, "genbeta");
  delete c1;
}
