#include "../interface/Tuple.h"

struct Tuple;

class TuplePlotter {
public:
  TuplePlotter(std::string DirName, int TypeMode, std::string SQRTS) {
    DirName_ = DirName;
    TypeMode_ = TypeMode;
    SQRTS_ = SQRTS;
  };
  ~TuplePlotter();

  void getObjects(Tuple*& tuple, TDirectory* dir);
  void draw(Tuple* tuple, TDirectory* dir, std::string LegendTitle, uint CutIndex);
  void CutFlow(Tuple* tuple, TDirectory* dir, uint CutIndex);
  //void CutFlowPlot(Tuple* tuple, TDirectory* dir, uint CutIndex, double ylow, double yhigh, bool setLog);
  void MassPrediction(Tuple* tuple, TDirectory* dir, std::string path, uint CutIndex);

private:
  std::string DirName_;
  int TypeMode_;
  std::string SQRTS_;
  float weight_;
};

//=============================================================
//
//     get objects: histos, trees
//
//=============================================================

void TuplePlotter::getObjects(Tuple*& tuple, TDirectory* dir) {
  std::string path = DirName_;
  std::string Name;

  Name = "IntLumi";
  tuple->IntLumi = (TProfile*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "XSection";
  tuple->XSection = (TProfile*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "NumEvents";
  tuple->NumEvents = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "N1_Eta";
  tuple->N1_Eta = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_Chi2PerNdof";
  tuple->N1_Chi2PerNdof = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_Qual";
  tuple->N1_Qual = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_TNOH";
  tuple->N1_TNOH = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_TNOPH";
  tuple->N1_TNOPH = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_TNOM";
  tuple->N1_TNOM = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_TNOHFraction";
  tuple->N1_TNOHFraction = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "nDof";
  tuple->nDof = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "tofError";
  tuple->tofError = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_PtErrOverPt";
  tuple->N1_PtErrOverPt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "TIsol";
  tuple->TIsol = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_EIsol";
  tuple->N1_EIsol = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_SumpTOverpT";
  tuple->N1_SumpTOverpT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_MPt";
  tuple->N1_MPt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_MIh";
  tuple->N1_MIh = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_MiniRelIsoChg";
  tuple->N1_MiniRelIsoChg  = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MTOF";
  tuple->MTOF = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Pt";
  tuple->Pt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "I";
  tuple->I = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "TOF";
  tuple->TOF = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "NVTrack";
  tuple->NVTrack = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_Stations";
  tuple->N1_Stations = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_Dxy";
  tuple->N1_Dxy = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_Dz";
  tuple->N1_Dz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_SegSep";
  tuple->N1_SegSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "FailDz";
  tuple->FailDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "CutFlow";
  tuple->CutFlow = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "CutFlowProbQFirst";
  tuple->CutFlowProbQFirst = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "CutFlowEta";
  tuple->CutFlowEta = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);


  Name = "N1_ProbQ";
  tuple->N1_ProbQ = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "ProbQNoL1";
  tuple->ProbQNoL1 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_ProbXY";
  tuple->N1_ProbXY = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "N1_pfType";
  tuple->N1_pfType = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "ProbXYNoL1";
  tuple->ProbXYNoL1 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "PrePreS_pfType";
  tuple->PrePreS_pfType = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "HSCPE";
  tuple->HSCPE = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "HSCPE_SystP";
  tuple->HSCPE_SystP = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "HSCPE_SystI";
  tuple->HSCPE_SystI = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "HSCPE_SystM";
  tuple->HSCPE_SystM = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "HSCPE_SystT";
  tuple->HSCPE_SystT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "HSCPE_SystPU";
  tuple->HSCPE_SystPU = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "HSCPE_SystHUp";
  tuple->HSCPE_SystHUp = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "HSCPE_SystHDown";
  tuple->HSCPE_SystHDown = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass";
  tuple->Mass = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF";
  tuple->MassTOF = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb";
  tuple->MassComb = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass";
  tuple->MaxEventMass = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_SystP";
  tuple->Mass_SystP = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF_SystP";
  tuple->MassTOF_SystP = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_SystP";
  tuple->MassComb_SystP = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass_SystP";
  tuple->MaxEventMass_SystP = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_SystI";
  tuple->Mass_SystI = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF_SystI";
  tuple->MassTOF_SystI = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_SystI";
  tuple->MassComb_SystI = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass_SystI";
  tuple->MaxEventMass_SystI = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_SystM";
  tuple->Mass_SystM = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF_SystM";
  tuple->MassTOF_SystM = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_SystM";
  tuple->MassComb_SystM = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass_SystM";
  tuple->MaxEventMass_SystM = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_SystT";
  tuple->Mass_SystT = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF_SystT";
  tuple->MassTOF_SystT = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_SystT";
  tuple->MassComb_SystT = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass_SystT";
  (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_SystPU";
  tuple->Mass_SystPU = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF_SystPU";
  tuple->MassTOF_SystPU = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_SystPU";
  tuple->MassComb_SystPU = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass_SystPU";
  tuple->MaxEventMass_SystPU = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_SystHUp";
  tuple->Mass_SystHUp = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF_SystH";
  tuple->MassTOF_SystH = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_SystHUp";
  tuple->MassComb_SystHUp = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass_SystHUp";
  tuple->MaxEventMass_SystHUp = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_SystHDown";
  tuple->Mass_SystHDown = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_SystHDown";
  tuple->MassComb_SystHDown = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MaxEventMass_SystHDown";
  tuple->MaxEventMass_SystHDown = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "Mass_Flip";
  tuple->Mass_Flip = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassTOF_Flip";
  tuple->MassTOF_Flip = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "MassComb_Flip";
  tuple->MassComb_Flip = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  //UNUSED//if(SkipSelectionPlot)return;

  Name = "Gen_DecayLength";
  tuple->Gen_DecayLength = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_Gen";
  tuple->Beta_Gen = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_GenChaged";
  tuple->Beta_GenCharged = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_Triggered";
  tuple->Beta_Triggered = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);


  Name = "genrecopT";
  tuple->genrecopT = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "Beta_Matched";
  tuple->Beta_Matched = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_PreselectedA";
  tuple->Beta_PreselectedA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_PreselectedB";
  tuple->Beta_PreselectedB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_PreselectedC";
  tuple->Beta_PreselectedC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_SelectedP";
  tuple->Beta_SelectedP = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_SelectedI";
  tuple->Beta_SelectedI = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "Beta_SelectedT";
  tuple->Beta_SelectedT = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "PrePreS_massT";
  tuple->PrePreS_massT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MiniRelIsoAll";
  tuple->PrePreS_MiniRelIsoAll = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MiniRelIsoChg";
  tuple->PrePreS_MiniRelIsoChg = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_RecoPFMET";
  tuple->PrePreS_RecoPFMET = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);


  Name = "PrePreS_Chi2PerNdof";
  tuple->PrePreS_Chi2PerNdof = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Qual";
  tuple->PrePreS_Qual = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOH";
  tuple->PrePreS_TNOH = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOH_PUA";
  tuple->PrePreS_TNOH_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOH_PUB";
  tuple->PrePreS_TNOH_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOHFraction";
  tuple->PrePreS_TNOHFraction = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOPH";
  tuple->PrePreS_TNOPH = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOHFractionTillLast";
  tuple->PrePreS_TNOHFractionTillLast = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOMHTillLast";
  tuple->PrePreS_TNOMHTillLast = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Eta";
  tuple->PrePreS_Eta = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOM";
  tuple->PrePreS_TNOM = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOM_PUA";
  tuple->PrePreS_TNOM_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TNOM_PUB";
  tuple->PrePreS_TNOM_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_nDof";
  tuple->PrePreS_nDof = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOFError";
  tuple->PrePreS_TOFError = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PtErrOverPt";
  tuple->PrePreS_PtErrOverPt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PtErrOverPt2";
  tuple->PrePreS_PtErrOverPt2 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MPt";
  tuple->PrePreS_MPt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MIs";
  tuple->PrePreS_MIs = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MIh";
  tuple->PrePreS_MIh = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MTOF";
  tuple->PrePreS_MTOF = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TIsol";
  tuple->PrePreS_TIsol = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_EIsol";
  tuple->PrePreS_EIsol = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_SumpTOverpT";
  tuple->PrePreS_SumpTOverpT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_LastHitDXY";
  tuple->PrePreS_LastHitDXY = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_LastHitD3D";
  tuple->PrePreS_LastHitD3D = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_P";
  tuple->PrePreS_P = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt";
  tuple->PrePreS_Pt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_PUA";
  tuple->PrePreS_Pt_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_PUB";
  tuple->PrePreS_Pt_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_Cosmic";
  tuple->PrePreS_Pt_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_DT";
  tuple->PrePreS_Pt_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_CSC";
  tuple->PrePreS_Pt_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Is";
  tuple->PrePreS_Is = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Is_PUA";
  tuple->PrePreS_Is_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Is_PUB";
  tuple->PrePreS_Is_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Is_Cosmic";
  tuple->PrePreS_Is_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Ih";
  tuple->PrePreS_Ih = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Ih_PUA";
  tuple->PrePreS_Ih_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Ih_PUB";
  tuple->PrePreS_Ih_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF";
  tuple->PrePreS_TOF = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF_PUA";
  tuple->PrePreS_TOF_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF_PUB";
  tuple->PrePreS_TOF_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF_DT";
  tuple->PrePreS_TOF_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF_CSC";
  tuple->PrePreS_TOF_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_dR_NVTrack";
  tuple->PrePreS_dR_NVTrack = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MatchedStations";
  tuple->PrePreS_MatchedStations = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_InnerInvPtDiff";
  tuple->PrePreS_InnerInvPtDiff = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Phi";
  tuple->PrePreS_Phi = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TimeAtIP";
  tuple->PrePreS_TimeAtIP = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_OpenAngle";
  tuple->PrePreS_OpenAngle = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_OpenAngle_Cosmic";
  tuple->PrePreS_OpenAngle_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "PrePreS_NVertex";
  tuple->PrePreS_NVertex = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_NVertex_NoEventWeight";
  (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PV";
  tuple->PrePreS_PV = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PV_NoEventWeight";
  tuple->PrePreS_PV_NoEventWeight = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_NOMoNOHvsPV";
  tuple->PrePreS_NOMoNOHvsPV = (TProfile*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_dzAll";
  tuple->PrePreS_dzAll = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_dxyAll";
  tuple->PrePreS_dxyAll = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_dzMinv3d";
  tuple->PrePreS_dzMinv3d = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_dxyMinv3d";
  tuple->PrePreS_dxyMinv3d = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "PrePreS_SegSep";
  tuple->PrePreS_SegSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_SegMinEtaSep";
  tuple->PrePreS_SegMinEtaSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_SegMinPhiSep";
  tuple->PrePreS_SegMinPhiSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_SegMinEtaSep_FailDz";
  tuple->PrePreS_SegMinEtaSep_FailDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_SegMinEtaSep_PassDz";
  tuple->PrePreS_SegMinEtaSep_PassDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Dz_FailSep";
  tuple->PrePreS_Dz_FailSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "PrePreS_Dxy";
  tuple->PrePreS_Dxy = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Dxy_Cosmic";
  tuple->PrePreS_Dxy_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Dz";
  tuple->PrePreS_Dz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Dz_Cosmic";
  tuple->PrePreS_Dz_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Dz_CSC";
  tuple->PrePreS_Dz_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Dz_DT";
  tuple->PrePreS_Dz_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_FailDz";
  tuple->PrePreS_Pt_FailDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_FailDz_DT";
  tuple->PrePreS_Pt_FailDz_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_Pt_FailDz_CSC";
  tuple->PrePreS_Pt_FailDz_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF_FailDz";
  tuple->PrePreS_TOF_FailDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF_FailDz_DT";
  tuple->PrePreS_TOF_FailDz_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_TOF_FailDz_CSC";
  tuple->PrePreS_TOF_FailDz_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PtErrOverPtVsPtErrOverPt2";
  tuple->PrePreS_PtErrOverPtVsPtErrOverPt2 = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PtErrOverPtVsPt";
  tuple->PrePreS_PtErrOverPtVsPt = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "PrePreS_ProbQ";
  tuple->PrePreS_ProbQ = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_ProbXY";
  tuple->PrePreS_ProbXY = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_ProbQNoL1";
  tuple->PrePreS_ProbQNoL1 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_ProbXYNoL1";
  tuple->PrePreS_ProbXYNoL1 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_MassErr";
  tuple->PrePreS_MassErr = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "PostPreS_massT";
  tuple->PostPreS_massT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MiniRelIsoAll";
  tuple->PostPreS_MiniRelIsoAll = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MiniRelIsoChg";
  tuple->PostPreS_MiniRelIsoChg = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_RecoPFMET";
  tuple->PostPreS_RecoPFMET = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_pfType";
  tuple->PostPreS_pfType = (TH1F*)GetObjectFromPath(dir, path + "/" + Name); 
  Name = "PostPreS_Chi2PerNdof";
  tuple->PostPreS_Chi2PerNdof = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Qual";
  tuple->PostPreS_Qual = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOH";
  tuple->PostPreS_TNOH = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOH_PUA";
  tuple->PostPreS_TNOH_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOH_PUB";
  tuple->PostPreS_TNOH_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOHFraction";
  tuple->PostPreS_TNOHFraction = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOPH";
  tuple->PostPreS_TNOPH = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOHFractionTillLast";
  tuple->PostPreS_TNOHFractionTillLast = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOMHTillLast";
  tuple->PostPreS_TNOMHTillLast = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Eta";
  tuple->PostPreS_Eta = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOM";
  tuple->PostPreS_TNOM = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOM_PUA";
  tuple->PostPreS_TNOM_PUA = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TNOM_PUB";
  tuple->PostPreS_TNOM_PUB = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_nDof";
  tuple->PostPreS_nDof = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TOFError";
  tuple->PostPreS_TOFError = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PtErrOverPt";
  tuple->PostPreS_PtErrOverPt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PtErrOverPt2";
  tuple->PostPreS_PtErrOverPt2 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Pt";
  tuple->PostPreS_Pt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIs";
  tuple->PostPreS_MIs = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIs_NoEventWeight";
  tuple->PostPreS_MIs_NoEventWeight = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIh";
  tuple->PostPreS_MIh = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIh_NoEventWeight";
  tuple->PostPreS_MIh_NoEventWeight = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MTOF";
  tuple->PostPreS_MTOF = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TIsol";
  tuple->PostPreS_TIsol = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_EIsol";
  tuple->PostPreS_EIsol = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_SumpTOverpT";
  tuple->PostPreS_SumpTOverpT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_LastHitDXY";
  tuple->PostPreS_LastHitDXY = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_LastHitD3D";
  tuple->PostPreS_LastHitD3D = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_P";
  tuple->PostPreS_P = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Pt";
  tuple->PostPreS_Pt = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Ih";
  tuple->PostPreS_Ih = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_dR_NVTrack";
  tuple->PostPreS_dR_NVTrack = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MatchedStations";
  tuple->PostPreS_MatchedStations = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_InnerInvPtDiff";
  tuple->PostPreS_InnerInvPtDiff = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Phi";
  tuple->PostPreS_Phi = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TimeAtIP";
  tuple->PostPreS_TimeAtIP = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_OpenAngle";
  tuple->PostPreS_OpenAngle = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_OpenAngle_Cosmic";
  tuple->PostPreS_OpenAngle_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "PostPreS_NVertex";
  tuple->PostPreS_NVertex = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_NVertex_NoEventWeight";
  (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PV";
  tuple->PostPreS_PV = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PV_NoEventWeight";
  tuple->PostPreS_PV_NoEventWeight = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_NOMoNOHvsPV";
  tuple->PostPreS_NOMoNOHvsPV = (TProfile*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_dzAll";
  tuple->PostPreS_dzAll = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_dxyAll";
  tuple->PostPreS_dxyAll = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_dzMinv3d";
  tuple->PostPreS_dzMinv3d = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_dxyMinv3d";
  tuple->PostPreS_dxyMinv3d = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "PostPreS_SegSep";
  tuple->PostPreS_SegSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_SegMinEtaSep";
  tuple->PostPreS_SegMinEtaSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_SegMinPhiSep";
  tuple->PostPreS_SegMinPhiSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_SegMinEtaSep_FailDz";
  tuple->PostPreS_SegMinEtaSep_FailDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_SegMinEtaSep_PassDz";
  tuple->PostPreS_SegMinEtaSep_PassDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Dz_FailSep";
  tuple->PostPreS_Dz_FailSep = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "PostPreS_Dxy";
  tuple->PostPreS_Dxy = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Dxy_Cosmic";
  tuple->PostPreS_Dxy_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Dz";
  tuple->PostPreS_Dz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Dz_Cosmic";
  tuple->PostPreS_Dz_Cosmic = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Dz_CSC";
  tuple->PostPreS_Dz_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Dz_DT";
  tuple->PostPreS_Dz_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Pt_FailDz";
  tuple->PostPreS_Pt_FailDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Pt_FailDz_DT";
  tuple->PostPreS_Pt_FailDz_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_Pt_FailDz_CSC";
  tuple->PostPreS_Pt_FailDz_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TOF_FailDz";
  tuple->PostPreS_TOF_FailDz = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TOF_FailDz_DT";
  tuple->PostPreS_TOF_FailDz_DT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_TOF_FailDz_CSC";
  tuple->PostPreS_TOF_FailDz_CSC = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PtErrOverPtVsPtErrOverPt2";
  tuple->PostPreS_PtErrOverPtVsPtErrOverPt2 = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PtErrOverPtVsPt";
  tuple->PostPreS_PtErrOverPtVsPt = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "PostPreS_ProbQ";
  tuple->PostPreS_ProbQ = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_ProbXY";
  tuple->PostPreS_ProbXY = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_ProbQNoL1";
  tuple->PostPreS_ProbQNoL1 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_ProbXYNoL1";
  tuple->PostPreS_ProbXYNoL1 = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MassErr";
  tuple->PostPreS_MassErr = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "PostPreS_EtaPerGenID";
  tuple->PostPreS_EtaPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_ProbQPerGenID";
  tuple->PostPreS_ProbQPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_ProbXYPerGenID";
  tuple->PostPreS_ProbXYPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PtPerGenID";
  tuple->PostPreS_PtPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_EIsolPerGenID";
  tuple->PostPreS_EIsolPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIhPerGenID";
  tuple->PostPreS_MIhPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIsPerGenID";
  tuple->PostPreS_MIsPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_massTPerGenID";
  tuple->PostPreS_massTPerGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_EtaPerMomGenID";
  tuple->PostPreS_EtaPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_ProbQPerMomGenID";
  tuple->PostPreS_ProbQPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_ProbXYPerMomGenID";
  tuple->PostPreS_ProbXYPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_PtPerMomGenID";
  tuple->PostPreS_PtPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_EIsolPerMomGenID";
  tuple->PostPreS_EIsolPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIhPerMomGenID";
  tuple->PostPreS_MIhPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_MIsPerMomGenID";
  tuple->PostPreS_MIsPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PostPreS_massTPerMomGenID";
  tuple->PostPreS_massTPerMomGenID  = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  
  Name = "genlevelpT";
  tuple->genlevelpT = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "genleveleta";
  tuple->genleveleta = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "genlevelbeta";
  tuple->genlevelbeta = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);

  //Initialize histograms for number of bins.
  for (int i = 0; i < 6; i++) {
    char Suffix[1024];
    sprintf(Suffix, "_%i", i);
    Name = "PrePreS_Pt_Binned";
    Name.append(Suffix);
    tuple->PrePreS_Pt_Binned[std::to_string(i)] = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
    Name = "PrePreS_TOF_Binned";
    Name.append(Suffix);
    tuple->PrePreS_TOF_Binned[std::to_string(i)] = (TH1F*)GetObjectFromPath(dir, path + "/" + Name);
  }

  Name = "AS_Eta_RegionA";
  tuple->AS_Eta_RegionA = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Eta_RegionB";
  tuple->AS_Eta_RegionB = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Eta_RegionC";
  tuple->AS_Eta_RegionC = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Eta_RegionD";
  tuple->AS_Eta_RegionD = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Eta_RegionE";
  tuple->AS_Eta_RegionE = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Eta_RegionF";
  tuple->AS_Eta_RegionF = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Eta_RegionG";
  tuple->AS_Eta_RegionG = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Eta_RegionH";
  tuple->AS_Eta_RegionH = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "AS_P";
  tuple->AS_P = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Pt";
  tuple->AS_Pt = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Is";
  tuple->AS_Is = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_Ih";
  tuple->AS_Ih = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_TOF";
  tuple->AS_TOF = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  Name = "PrePreS_EtaIs";
  tuple->PrePreS_EtaIs = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_EtaIh";
  tuple->PrePreS_EtaIh = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_EtaP";
  tuple->PrePreS_EtaP = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_EtaPt";
  tuple->PrePreS_EtaPt = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_EtaTOF";
  tuple->PrePreS_EtaTOF = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_EtaNBH";
  tuple->PrePreS_EtaNBH = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_EtaDz";
  tuple->PrePreS_EtaDz = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PIs";
  tuple->PrePreS_PIs = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_IhIs";
  tuple->PrePreS_IhIs = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PIh";
  tuple->PrePreS_PIh = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PtIs";
  tuple->PrePreS_PtIs = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PtIh";
  tuple->PrePreS_PtIh = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "PrePreS_PtTOF";
  tuple->PrePreS_PtTOF = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  //   Name = "PrePreS_TOFIs"; tuple->PrePreS_TOFIs = (TH2F*)GetObjectFromPath(dir, path+"/"+Name);
  Name = "PrePreS_TOFIs";
  tuple->PrePreS_TOFIs = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);
  //   Name = "PrePreS_TOFIh"; PrePreS_TOFIm = (TH2F*)GetObjectFromPath(dir, path+"/"+Name);
  Name = "PrePreS_TOFIh";
  tuple->PrePreS_TOFIh = (TH2F*)GetObjectFromPath(dir, path + "/" + Name);

  //   Name = "AS_EtaIs"; AS_EtaIs = (TH3F*)GetObjectFromPath(dir, path+"/"+Name);
  //   Name = "AS_EtaIh"; AS_EtaIm = (TH3F*)GetObjectFromPath(dir, path+"/"+Name);
  //   Name = "AS_EtaP" ; AS_EtaP  = (TH3F*)GetObjectFromPath(dir, path+"/"+Name);
  //   Name = "AS_EtaPt"; AS_EtaPt = (TH3F*)GetObjectFromPath(dir, path+"/"+Name);
  //   Name = "AS_EtaTOF"; AS_EtaTOF = (TH3F*)GetObjectFromPath(dir, path+"/"+Name);
  Name = "AS_PIs";
  tuple->AS_PIs = (TH3F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_PIh";
  tuple->AS_PIh = (TH3F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_PtIs";
  tuple->AS_PtIs = (TH3F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_PtIh";
  tuple->AS_PtIh = (TH3F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_TOFIs";
  tuple->AS_TOFIs = (TH3F*)GetObjectFromPath(dir, path + "/" + Name);
  Name = "AS_TOFIh";
  tuple->AS_TOFIh = (TH3F*)GetObjectFromPath(dir, path + "/" + Name);
}

//=============================================================
//
//  Draw all plots that are not meant for comparison with other samples (mostly 2D plots that can't be superimposed)
//
//=============================================================
void TuplePlotter::draw(Tuple* tuple, TDirectory* dir, std::string LegendTitle, uint CutIndex) {
  TH1D* HCuts_Pt = (TH1D*)GetObjectFromPath(dir, "HCuts_Pt");
  TH1D* HCuts_Is = (TH1D*)GetObjectFromPath(dir, "HCuts_I");
  TH1D* HCuts_TOF = (TH1D*)GetObjectFromPath(dir, "HCuts_TOF");
  char PtCutStr[1024];
  sprintf(PtCutStr, "%.0f GeV", HCuts_Pt->GetBinContent(CutIndex + 1));
  char ICutStr[1024];
  sprintf(ICutStr, "%.2f", HCuts_Is->GetBinContent(CutIndex + 1));
  char TOFCutStr[1024];
  sprintf(TOFCutStr, "%.3f", HCuts_TOF->GetBinContent(CutIndex + 1));

  TObject** Histos = new TObject*[10];
  std::vector<std::string> legend;

  TCanvas* c1 = new TCanvas("c1", "c1,", 600, 600);
  legend.clear();
  Histos[0] = (TH1*)tuple->PrePreS_EtaIs;
  legend.push_back("Before Cut");
  std::string dEdxS_Legend = "I_{as}";
  DrawSuperposedHistos((TH1**)Histos, legend, "COLZ", "#eta", dEdxS_Legend.c_str(), 0, 0, 0, 0, false);
  c1->SetLogz(true);
  c1->SetRightMargin(0.15);
  //DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
  c1->SaveAs(TString(DirName_) + "/EtaIs_BS.png");  //SaveCanvas(c1,SavePath,"EtaIs_BS", true);
  delete c1;
}

/*SelectionPlot*/

void TuplePlotter::MassPrediction(Tuple* tuple, TDirectory* dir, std::string path, uint CutIndex) {
  /*bool IsTkOnly = (TypeMode_==0);
   TH1D  *Pred = nullptr, *Data = nullptr;*/
  TH1D* Pred = GetProjectionFromPath(dir, path, CutIndex, "TmpPredMass");
  //"Data/Pred_" + HistoSuffix
  /*Data  = GetProjectionFromPath(dir, path, CutIndex, "TmpDataMass");*/
  TCanvas* c1 = new TCanvas("c1", "c1,", 600, 600);
  if (!Pred)
    std::cout << "Error: Histogram " << path << " NOT found." << std::endl;
  else {
    Pred->Draw();
    c1->SaveAs(TString(DirName_) + "/Pred_Mass.png");
  }
  TCanvas* c2 = new TCanvas("c2", "c2,", 600, 600);
  TH2F* Pred2D = (TH2F*)GetObjectFromPath(dir, path);
  Pred2D->Draw("colz");
  c2->SaveAs(TString(DirName_) + "/Pred_Mass2D.png");
  /*TCanvas* c2 = new TCanvas("c2","c2,",600,600);
   Data->Draw();*/
}

void TuplePlotter::CutFlow(Tuple* tuple, TDirectory* dir, uint CutIndex) {
  TH1D* HCuts_Pt = (TH1D*)GetObjectFromPath(dir, "HCuts_Pt");
  TH1D* HCuts_I = (TH1D*)GetObjectFromPath(dir, "HCuts_I");
  TH1D* HCuts_TOF = (TH1D*)GetObjectFromPath(dir, "HCuts_TOF");

  char Buffer[1024];
  sprintf(Buffer,
          "%s/CutFlow_%03i_Pt%03.0f_I%05.3f_TOF%04.3f.txt",
          DirName_.c_str(),
          CutIndex,
          HCuts_Pt->GetBinContent(CutIndex + 1),
          HCuts_I->GetBinContent(CutIndex + 1),
          HCuts_TOF->GetBinContent(CutIndex + 1));
/*
  FILE* pFile = fopen(Buffer, "w");

  fprintf(pFile, "#################### %20s ####################\n", DirName_.c_str());
  fprintf(pFile, "#Events                       = %4.2E\n", tuple->TotalE->Integral(0,tuple->TotalE->GetNbinsX()+1));
  fprintf(pFile,
          "#Triggered Events             = %4.2E Eff=%4.3E\n",
          tuple->TotalTE->Integral(0,tuple->TotalTE->GetNbinsX()+1),
          tuple->TotalTE->Integral(0,tuple->TotalTE->GetNbinsX()+1) / tuple->TotalE->Integral(0,tuple->TotalE->GetNbinsX()+1));
  fprintf(pFile, "#Tracks                       = %4.2E\n", tuple->Total->Integral(0,tuple->Total->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing TNOH   cuts   = %4.2E Eff=%4.3E\n",
          tuple->TNOH->Integral(0,tuple->TNOH->GetNbinsX()+1),
          tuple->TNOH->Integral(0,tuple->TNOH->GetNbinsX()+1) / tuple->Total->Integral(0,tuple->Total->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing TNOM   cuts   = %4.2E Eff=%4.3E\n",
          tuple->TNOM->Integral(0,tuple->TNOM->GetNbinsX()+1),
          tuple->TNOM->Integral(0,tuple->TNOM->GetNbinsX()+1) / tuple->TNOH->Integral(0,tuple->TNOH->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing nDof   cuts   = %4.2E Eff=%4.3E\n",
          tuple->nDof->Integral(0,tuple->nDof->GetNbinsX()+1),
          tuple->nDof->Integral(0,tuple->nDof->GetNbinsX()+1) / tuple->TNOM->Integral(0,tuple->TNOM->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Qual   cuts   = %4.2E Eff=%4.3E\n",
          tuple->Qual->Integral(0,tuple->Qual->GetNbinsX()+1),
          tuple->Qual->Integral(0,tuple->Qual->GetNbinsX()+1) / tuple->nDof->Integral(0,tuple->nDof->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Chi2PerNdof   cuts   = %4.2E Eff=%4.3E\n",
          tuple->Chi2PerNdof->Integral(0,tuple->Chi2PerNdof->GetNbinsX()+1),
          tuple->Chi2PerNdof->Integral(0,tuple->Chi2PerNdof->GetNbinsX()+1) / tuple->Qual->Integral(0,tuple->Qual->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Min Pt cuts   = %4.2E Eff=%4.3E\n",
          tuple->MPt->Integral(0,tuple->MPt->GetNbinsX()+1),
          tuple->MPt->Integral(0,tuple->MPt->GetNbinsX()+1) / tuple->Chi2PerNdof->Integral(0,tuple->Chi2PerNdof->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Min I  cuts   = %4.2E Eff=%4.3E\n",
          tuple->MI->Integral(0,tuple->MI->GetNbinsX()+1),
          tuple->MI->Integral(0,tuple->MI->GetNbinsX()+1) / tuple->MPt->Integral(0,tuple->MPt->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Min TOFcuts   = %4.2E Eff=%4.3E\n",
          tuple->MTOF->Integral(0,tuple->MTOF->GetNbinsX()+1),
          tuple->MTOF->Integral(0,tuple->MTOF->GetNbinsX()+1) / tuple->MI->Integral(0,tuple->MI->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Dxy    cuts   = %4.2E Eff=%4.3E\n",
          tuple->Dxy->Integral(0,tuple->Dxy->GetNbinsX()+1),
          tuple->Dxy->Integral(0,tuple->Dxy->GetNbinsX()+1) / tuple->MTOF->Integral(0,tuple->MTOF->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing TIsol  cuts   = %4.2E Eff=%4.3E\n",
          tuple->TIsol->Integral(0,tuple->TIsol->GetNbinsX()+1),
          tuple->TIsol->Integral(0,tuple->TIsol->GetNbinsX()+1) / tuple->Dxy->Integral(0,tuple->Dxy->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing EIsol  cuts   = %4.2E Eff=%4.3E\n",
          tuple->EIsol->Integral(0,tuple->EIsol->GetNbinsX()+1),
          tuple->EIsol->Integral(0,tuple->EIsol->GetNbinsX()+1) / tuple->TIsol->Integral(0,tuple->TIsol->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing PtErr  cuts   = %4.2E Eff=%4.3E\n",
          tuple->Pterr->Integral(0,tuple->Pterr->GetNbinsX()+1),
          tuple->Pterr->Integral(0,tuple->Pterr->GetNbinsX()+1) / tuple->EIsol->Integral(0,tuple->EIsol->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Dz  cuts      = %4.2E Eff=%4.3E\n",
          tuple->Dz->Integral(0,tuple->Dz->GetNbinsX()+1),
          tuple->Dz->Integral(0,tuple->Dz->GetNbinsX()+1) / tuple->Pterr->Integral(0,tuple->Pterr->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Basic  cuts   = %4.2E Eff=%4.3E\n",
          tuple->Basic->Integral(0,tuple->Basic->GetNbinsX()+1),
          tuple->Basic->Integral(0,tuple->Basic->GetNbinsX()+1) / tuple->Total->Integral(0,tuple->Total->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing Pt     cuts   = %4.2E Eff=%4.3E\n",
          tuple->Pt->GetBinContent(CutIndex + 1),
          tuple->Pt->GetBinContent(CutIndex + 1) / tuple->Pterr->Integral(0,tuple->Pterr->GetNbinsX()+1));
  fprintf(pFile,
          "#Tracks passing I      cuts   = %4.2E Eff=%4.3E\n",
          tuple->I->GetBinContent(CutIndex + 1),
          tuple->I->GetBinContent(CutIndex + 1) / tuple->Pt->GetBinContent(CutIndex + 1));
  fprintf(pFile,
          "#Tracks passing TOF    cuts   = %4.2E Eff=%4.3E\n",
          tuple->TOF->GetBinContent(CutIndex + 1),
          tuple->TOF->GetBinContent(CutIndex + 1) / tuple->I->GetBinContent(CutIndex + 1));
  fprintf(pFile,
          "#Tracks passing selection     = %4.2E Eff=%4.3E\n",
          tuple->TOF->GetBinContent(CutIndex + 1),
          tuple->TOF->GetBinContent(CutIndex + 1) / tuple->Total->Integral(0,tuple->Total->GetNbinsX()+1));
  fprintf(pFile, "--------------------\n");
  fprintf(pFile,
          "HSCP Detection Efficiency Before Trigger                           Eff=%4.3E\n",
          tuple->TOF->GetBinContent(CutIndex + 1) / (2 * tuple->TotalE->Integral(0,tuple->TotalE->GetNbinsX()+1)));
  fprintf(pFile,
          "HSCP Detection Efficiency After  Trigger                           Eff=%4.3E\n",
          tuple->TOF->GetBinContent(CutIndex + 1) / (2 * tuple->TotalTE->Integral(0,tuple->TotalTE->GetNbinsX()+1)));
  fprintf(pFile,
          "#HSCPTrack per HSCPEvent (with at least one HSCPTrack)             Eff=%4.3E\n",
          tuple->TOF->GetBinContent(CutIndex + 1) / (tuple->HSCPE->GetBinContent(CutIndex + 1)));
  fprintf(pFile,
          "HSCP Event Efficiency                                              Eff=%4.3E\n",
          tuple->HSCPE->GetBinContent(CutIndex + 1) / (tuple->TotalE->Integral(0,tuple->TotalE->GetNbinsX()+1)));

  fprintf(pFile, "\n\n");

  fclose(pFile);
*/
}

/*void TuplePlotter::CutFlowPlot(Tuple* tuple, TDirectory* dir, uint CutIndex, double ylow, double yhigh, bool setLog){
   TH1D* HCuts_Pt  = (TH1D*) GetObjectFromPath (dir, "HCuts_Pt");
   TH1D* HCuts_Is  = (TH1D*) GetObjectFromPath (dir, "HCuts_I");
   TH1D* HCuts_TOF = (TH1D*) GetObjectFromPath (dir, "HCuts_TOF");
   char PtCutStr [1024]; sprintf (PtCutStr, "p_{T}>%.0f GeV", HCuts_Pt ->GetBinContent(CutIndex+1));
   char ICutStr  [1024]; sprintf (ICutStr,  "I_{as}>%.2f",    HCuts_Is ->GetBinContent(CutIndex+1));
   char TOFCutStr[1024]; sprintf (TOFCutStr,"1/#beta>%.3f",   HCuts_TOF->GetBinContent(CutIndex+1));
}*/
