#include "SUSYBSMAnalysis/Analyzer/interface/Tuple.h"

struct Tuple;

class TupleMaker {
public:
  TupleMaker();
  ~TupleMaker();

  void initializeTuple(Tuple *&tuple,
                       TFileDirectory &dir,
                       unsigned int saveTree,
                       unsigned int saveGenTree,
                       bool SkipSelectionPlot,
                       int TypeMode,
                       bool isSignal,
                       unsigned int NCuts,
                       unsigned int NCuts_Flip,
                       float PtHistoUpperBound,
                       float MassHistoUpperBound,
                       int MassNBins,
                       float IPbound,
                       int PredBins,
                       int EtaBins,
                       float dEdxS_UpLim,
                       float dEdxM_UpLim,
                       int DzRegions,
                       float GlobalMinPt,
                       float GlobalMinTOF);

  void initializeRegions(Tuple *&tuple,
                        TFileDirectory &dir,
                        int etabins,
                        int ihbins,
                        int pbins,
                        int massbins);

  void fillTreeBranches(Tuple *&tuple,
                        const unsigned int &Trig,
                        const unsigned int &Run,
                        const unsigned long &Event,
                        const unsigned int &Lumi,
                        const unsigned int &PileUp,
                        const unsigned int &nofVertices,
                        const unsigned int &Hscp,
                        const unsigned int &nmuons,
                        const unsigned int &njets,
                        const float &weight,
                        const float &generator_weight,
                        const float &generator_binning_values,
                        const bool &HLT_Mu50,
                        const bool &HLT_PFMET120_PFMHT120_IDTight,
                        const bool &HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                        const bool &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                        const bool &HLT_MET105_IsoTrk50,
                        const float &CaloMET,
                        const float &RecoPFMET,
                        const float &RecoPFMHT,
                        const float &HLTPFMET,
                        const float &HLTPFMHT,
                        const float &RecoPFMET_eta,
                        const float &RecoPFMET_phi,
                        const float &RecoPFMET_significance,
                        const float &Muon1_Pt,
                        const float &Muon1_eta,
                        const float &Muon1_phi,
                        const float &Muon2_Pt,
                        const float &Muon2_eta,
                        const float &Muon2_phi,
                        const std::vector<float> &vect_jet_pt,
                        const std::vector<float> &vect_jet_eta,
                        const std::vector<float> &vect_jet_phi,
                        const std::vector<float> &vect_jet_mass,
                        const std::vector<float> &vect_jet_energy,
                        const std::vector<float> &vect_jet_pdgId,
                        const std::vector<float> &vect_jet_et,
                        const std::vector<float> &vect_jet_chargedEmEnergyFraction,
                        const std::vector<float> &vect_jet_neutralEmEnergyFraction,
                        const std::vector<float> &vect_mT,
                        const std::vector<bool> &passCutPt55,
                        const std::vector<bool> &passPreselection,
                        const std::vector<bool> &passSelection,
                        const std::vector<float> &Charge,
                        const std::vector<float> &Pt,
                        const std::vector<float> &PtErr,
                        const std::vector<float> &Ias,
                        const std::vector<float> &Ias_noTIBnoTIDno3TEC,
                        const std::vector<float> &Ias_PixelOnly,
                        const std::vector<float> &Ias_StripOnly,
                        const std::vector<float> &Ias_PixelOnly_noL1,
                        const std::vector<float> &Ih,
                        const std::vector<float> &Ick,
                        const std::vector<float> &Fmip,
                        const std::vector<float> &ProbXY,
                        const std::vector<float> &ProbXY_noL1,
                        const std::vector<float> &ProbQ,
                        const std::vector<float> &ProbQ_noL1,
                        const std::vector<float> &Ndof,
                        const std::vector<float> &Chi2,
                        const std::vector<int>   &QualityMask,
                        const std::vector<bool>  &isHighPurity,
                        const std::vector<float>  &EoverP,
                        const std::vector<bool>  &isMuon,
                        const std::vector<int>   &MuonSelector,
                        const std::vector<bool>  &isElectron,
                        const std::vector<bool>  &isChHadron,
                        const std::vector<bool>  &isNeutHadron,
                        const std::vector<float> &ECAL_energy,
                        const std::vector<float> &HCAL_energy,
                        const std::vector<float> &TOF,
                        const std::vector<float> &TOFErr,
                        const std::vector<unsigned int> &TOF_ndof,
                        const std::vector<float> &DTTOF,
                        const std::vector<float> &DTTOFErr,
                        const std::vector<unsigned int> &DTTOF_ndof,
                        const std::vector<float> &CSCTOF,
                        const std::vector<float> &CSCTOFErr,
                        const std::vector<unsigned int> &CSCTOF_ndof,
                        const std::vector<float> &Mass,
                        const std::vector<float> &MassErr,
                        const std::vector<float> &dZ,
                        const std::vector<float> &dXY,
                        const std::vector<float> &dR,
                        const std::vector<float> &eta,
                        const std::vector<float> &phi,
                        const std::vector<unsigned int> &noh,
                        const std::vector<unsigned int> &noph,
                        const std::vector<float> &fovh,
                        const std::vector<unsigned int> &nomh,
                        const std::vector<float> &fovhd,
                        const std::vector<unsigned int> &nom,
                        const std::vector<float> &iso_TK,
                        const std::vector<float> &iso_ECAL,
                        const std::vector<float> &iso_HCAL,
                        const std::vector<float> &PFMiniIso_relative,
                        const std::vector<float> &PFMiniIso_wMuon_relative,
                        const std::vector<float> &track_PFIsolationR005_sumChargedHadronPt,
                        const std::vector<float> &track_PFIsolationR005_sumNeutralHadronPt,
                        const std::vector<float> &track_PFIsolationR005_sumPhotonPt,
                        const std::vector<float> &track_PFIsolationR005_sumPUPt,
                        const std::vector<float> &track_PFIsolationR01_sumChargedHadronPt,
                        const std::vector<float> &track_PFIsolationR01_sumNeutralHadronPt,
                        const std::vector<float> &track_PFIsolationR01_sumPhotonPt,
                        const std::vector<float> &track_PFIsolationR01_sumPUPt,
                        const std::vector<float> &track_PFIsolationR03_sumChargedHadronPt,
                        const std::vector<float> &track_PFIsolationR03_sumNeutralHadronPt,
                        const std::vector<float> &track_PFIsolationR03_sumPhotonPt,
                        const std::vector<float> &track_PFIsolationR03_sumPUPt,
                        const std::vector<float> &track_PFIsolationR05_sumChargedHadronPt,
                        const std::vector<float> &track_PFIsolationR05_sumNeutralHadronPt,
                        const std::vector<float> &track_PFIsolationR05_sumPhotonPt,
                        const std::vector<float> &track_PFIsolationR05_sumPUPt,
                        const std::vector<float> &muon_PFIsolationR03_sumChargedHadronPt,
                        const std::vector<float> &muon_PFIsolationR03_sumNeutralHadronPt,
                        const std::vector<float> &muon_PFIsolationR03_sumPhotonPt,
                        const std::vector<float> &muon_PFIsolationR03_sumPUPt,
                        const std::vector<float> &Ih_noL1,
                        const std::vector<float> &Ih_15drop,
                        const std::vector<float> &Ih_StripOnly,
                        const std::vector<float> &Ih_StripOnly_15drop,
                        const std::vector<float> &Ih_PixelOnly_noL1,
                        const std::vector<float> &Ih_SaturationCorrectionFromFits,
                        const std::vector<std::vector<float>> &clust_charge,
                        const std::vector<std::vector<float>> &clust_pathlength,
                        const std::vector<std::vector<unsigned int>> &clust_nstrip,
                        const std::vector<std::vector<bool>> &clust_sat254,
                        const std::vector<std::vector<bool>> &clust_sat255,
                        const std::vector<std::vector<uint32_t>> &clust_detid,
                        const std::vector<std::vector<bool>> &clust_isStrip,
                        const std::vector<std::vector<bool>> &clust_isPixel,
                        const std::vector<float> &genid,
                        const std::vector<float> &gencharge,
                        const std::vector<float> &genmass,
                        const std::vector<float> &genpt,
                        const std::vector<float> &geneta,
                        const std::vector<float> &genphi);

  void fillGenTreeBranches(Tuple *&tuple,
                           const unsigned int &Run,
                           const unsigned long &Event,
                           const unsigned int &Lumi,
                           const float &weight,
                           const float &generator_weight,
                           const float &generator_binning_values,
                           const std::vector<float> &genid,
                           const std::vector<float> &gencharge,
                           const std::vector<float> &genmass,
                           const std::vector<float> &genpt,
                           const std::vector<float> &geneta,
                           const std::vector<float> &genphi);

  void fillControlAndPredictionHist(const susybsm::HSCParticle &hscp,
                                    const reco::DeDxData *dedxSObj,
                                    const reco::DeDxData *dedxMObj,
                                    const reco::MuonTimeExtra *tof,
                                    Tuple *&tuple,
                                    int TypeMode,
                                    float GlobalMinTOF,
                                    float Event_Weight,
                                    bool isCosmicSB,
                                    float DTRegion,
                                    const int MaxPredBins,
                                    float DeDxK,
                                    float DeDxC,
                                    std::vector<float> CutPt,
                                    std::vector<float> CutI,
                                    std::vector<float> CutTOF,
                                    std::vector<float> CutPt_Flip,
                                    std::vector<float> CutI_Flip,
                                    std::vector<float> CutTOF_Flip);

  void fillRegions(Tuple *&tuple,
                   float pt_cut,
                   float Ias_quantiles[6],
                   float eta,
                   float p,
                   float pt,
                   float pterr,
                   float ih,
                   float ias,
                   float m,
                   float tof,
                   float w);

  void writeRegions(Tuple *&tuple,
                    TFileDirectory &dir);
};

//=============================================================
//
//     Create histograms, tree variables and branches
//
//=============================================================

void TupleMaker::initializeTuple(Tuple *&tuple,
                                 TFileDirectory &dir,
                                 unsigned int saveTree,
                                 unsigned int saveGenTree,
                                 bool SkipSelectionPlot,
                                 int TypeMode,
                                 bool isSignal,
                                 unsigned int NCuts,
                                 unsigned int NCuts_Flip,
                                 float PtHistoUpperBound,
                                 float MassHistoUpperBound,
                                 int MassNBins,
                                 float IPbound,
                                 int PredBins,
                                 int EtaBins,
                                 float dEdxS_UpLim,
                                 float dEdxM_UpLim,
                                 int DzRegions,
                                 float GlobalMinPt,
                                 float GlobalMinTOF) {
  std::string Name;

  Name = "IntLumi";
  tuple->IntLumi = dir.make<TProfile>(Name.c_str(), Name.c_str(), 1, 0, 1);
  Name = "XSection";
  tuple->XSection = dir.make<TProfile>(Name.c_str(), Name.c_str(), 1, 0, 1);
  Name = "NumEvents";
  tuple->NumEvents = dir.make<TH1F>(Name.c_str(), Name.c_str(), 3, 0., 3.);
  Name = "ErrorHisto";
  tuple->ErrorHisto = dir.make<TH1F>(Name.c_str(), Name.c_str(), 10, 0., 10.);
  tuple->ErrorHisto->Sumw2();
  Name = "BefPreS_TriggerType";
  tuple->BefPreS_TriggerType = dir.make<TH1F>(Name.c_str(), Name.c_str(), 3, 0., 3.);
  tuple->BefPreS_TriggerType->Sumw2();
  Name = "HSCPCandidateType";
  tuple->HSCPCandidateType = dir.make<TH1F>(Name.c_str(), Name.c_str(), 6, 0., 6.);

  tuple->CutFlow = dir.make<TH1F>("CutFlow", ";CutFlowIndex", 17, 0., 17.);
  tuple->CutFlowReverse = dir.make<TH1F>("CutFlowReverse", ";CutFlowIndex", 17, 0., 17.);
  tuple->CutFlowProbQ =  dir.make<TH2F>("CutFlowProbQ",";ProbQ;CutFlowIndex",10, 0., 1., 17, 0., 17.);
  tuple->CutFlowEta = dir.make<TH2F>("CutFlowEta", ";Eta;CutFlowIndex", 50, -2.6, 2.6, 17, 0., 17.);
  tuple->CutFlowPfType = dir.make<TH2F>("CutFlowPfType", ";PfType;CutFlowIndex", 9, 0., 9., 17, 0., 17.);

  Name = "N1_Eta";
  tuple->N1_Eta = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -2.6, 2.6);
  Name = "N1_Chi2oNdof";
  tuple->N1_Chi2oNdof = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  Name = "N1_Qual";
  tuple->N1_Qual = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  tuple->N1_TNOM = dir.make<TH1F>("N1_TNOM", "N1_TNOM;Number of measurments", 40, 0, 40);
  Name = "N1_TNOPH";
  tuple->N1_TNOPH = dir.make<TH1F>(Name.c_str(), Name.c_str(), 8, 0, 8);
  Name = "N1_TNOHFraction";
  tuple->N1_TNOHFraction = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 1);
  Name = "nDof";
  tuple->nDof = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  Name = "tofError";
  tuple->tofError = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 0.25);
  Name = "TIsol";
  tuple->TIsol = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 100);
  Name = "N1_EoP";
  tuple->N1_EoP = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 1.5);
  tuple->N1_dRMinPfJet= dir.make<TH1F>("N1_dRMinPfJet",";dRMinPfJet",100,0.,1.5);
  tuple->N1_dRMinPfJet->Sumw2();
  Name = "N1_SumpTOverpT";
  tuple->N1_SumpTOverpT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 80, 0, 2);
  Name = "N1_Pt";
  tuple->N1_Pt = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  Name = "N1_Ih";
  tuple->N1_Ih = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  Name = "MTOF";
  tuple->MTOF = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -2, 5);
  Name = "Pt";
  tuple->Pt = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  Name = "I";
  tuple->I = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  Name = "TOF";
  tuple->TOF = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->TOF->Sumw2();
  Name = "NVTrack";
  tuple->NVTrack = dir.make<TH1F>(Name.c_str(), Name.c_str(), 1, 0, 1);
  Name = "N1_Stations";
  tuple->N1_Stations = dir.make<TH1F>(Name.c_str(), Name.c_str(), 1, 0, 1);
  Name = "N1_Dxy";
  tuple->N1_Dxy = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.10, 0.10);
  Name = "N1_Dz";
  tuple->N1_Dz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.10, 0.10);
  Name = "N1_PtErrOverPt";
  tuple->N1_PtErrOverPt = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 1);
  tuple->N1_PtErrOverPt->Sumw2();
  Name = "N1_SegSep";
  tuple->N1_SegSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 1, 0, 1);
  Name = "FailDz";
  tuple->FailDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 1, 0, 1);

  Name = "N1_ProbQ";
  tuple->N1_ProbQ = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->N1_ProbQ->Sumw2();
  Name = "N1_ProbQVsIas";
  tuple->N1_ProbQVsIas = dir.make<TH2F>(Name.c_str(), Name.c_str(), 100, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->N1_ProbQVsIas->Sumw2();
  Name = "ProbQNoL1";
  tuple->ProbQNoL1 = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->ProbQNoL1->Sumw2();
  Name = "N1_ProbXY";
  tuple->N1_ProbXY = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->N1_ProbXY->Sumw2();
  Name = "N1_MiniRelIsoAll";
  tuple->N1_MiniRelIsoAll = dir.make<TH1F>(Name.c_str(), Name.c_str(),  200, 0.0, 10.0);
  tuple->N1_MiniRelIsoAll->Sumw2();
  Name = "ProbXYNoL1";
  tuple->ProbXYNoL1 = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->ProbXYNoL1->Sumw2();
  Name = "N1_pfType";
  tuple->N1_pfType = dir.make<TH1F>(Name.c_str(), Name.c_str(), 9, 0, 9);
  tuple->N1_pfType->Sumw2();

  Name = "HSCPE";
  tuple->HSCPE = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE->Sumw2();
  Name = "HSCPE_SystP";
  tuple->HSCPE_SystP = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE_SystP->Sumw2();
  Name = "HSCPE_SystI";
  tuple->HSCPE_SystI = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE_SystI->Sumw2();
  Name = "HSCPE_SystM";
  tuple->HSCPE_SystM = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE_SystM->Sumw2();
  Name = "HSCPE_SystT";
  tuple->HSCPE_SystT = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE_SystT->Sumw2();
  Name = "HSCPE_SystPU";
  tuple->HSCPE_SystPU = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE_SystPU->Sumw2();
  Name = "HSCPE_SystHUp";
  tuple->HSCPE_SystHUp = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE_SystHUp->Sumw2();
  Name = "HSCPE_SystHDown";
  tuple->HSCPE_SystHDown = dir.make<TH1F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
  tuple->HSCPE_SystHDown->Sumw2();

  Name = "Mass";
  tuple->Mass = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass->Sumw2();
  Name = "MassTOF";
  tuple->MassTOF = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF->Sumw2();
  Name = "MassComb";
  tuple->MassComb = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb->Sumw2();
  Name = "MaxEventMass";
  tuple->MaxEventMass = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass->Sumw2();

  Name = "Mass_SystP";
  tuple->Mass_SystP = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_SystP->Sumw2();
  Name = "MassTOF_SystP";
  tuple->MassTOF_SystP = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystP->Sumw2();
  Name = "MassComb_SystP";
  tuple->MassComb_SystP =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystP->Sumw2();
  Name = "MaxEventMass_SystP";
  tuple->MaxEventMass_SystP =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystP->Sumw2();

  Name = "Mass_SystI";
  tuple->Mass_SystI = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_SystI->Sumw2();
  Name = "MassTOF_SystI";
  tuple->MassTOF_SystI = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystI->Sumw2();
  Name = "MassComb_SystI";
  tuple->MassComb_SystI =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystI->Sumw2();
  Name = "MaxEventMass_SystI";
  tuple->MaxEventMass_SystI =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystI->Sumw2();

  Name = "Mass_SystM";
  tuple->Mass_SystM = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_SystM->Sumw2();
  Name = "MassTOF_SystM";
  tuple->MassTOF_SystM = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystM->Sumw2();
  Name = "MassComb_SystM";
  tuple->MassComb_SystM =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystM->Sumw2();
  Name = "MaxEventMass_SystM";
  tuple->MaxEventMass_SystM =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystM->Sumw2();

  Name = "Mass_SystT";
  tuple->Mass_SystT = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_SystT->Sumw2();
  Name = "MassTOF_SystT";
  tuple->MassTOF_SystT = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystT->Sumw2();
  Name = "MassComb_SystT";
  tuple->MassComb_SystT =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystT->Sumw2();
  Name = "MaxEventMass_SystT";
  tuple->MaxEventMass_SystT =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystT->Sumw2();

  Name = "Mass_SystPU";
  tuple->Mass_SystPU = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_SystPU->Sumw2();
  Name = "MassTOF_SystPU";
  tuple->MassTOF_SystPU =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystPU->Sumw2();
  Name = "MassComb_SystPU";
  tuple->MassComb_SystPU =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystPU->Sumw2();
  Name = "MaxEventMass_SystPU";
  tuple->MaxEventMass_SystPU =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystPU->Sumw2();

  Name = "Mass_SystHUp";
  tuple->Mass_SystHUp = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_SystHUp->Sumw2();
  Name = "MassTOF_SystH";
  tuple->MassTOF_SystH = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystH->Sumw2();
  Name = "MassComb_SystHUp";
  tuple->MassComb_SystHUp =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystHUp->Sumw2();
  Name = "MaxEventMass_SystHUp";
  tuple->MaxEventMass_SystHUp =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystHUp->Sumw2();

  Name = "Mass_SystHDown";
  tuple->Mass_SystHDown =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_SystHDown->Sumw2();
  Name = "MassComb_SystHDown";
  tuple->MassComb_SystHDown =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystHDown->Sumw2();
  Name = "MaxEventMass_SystHDown";
  tuple->MaxEventMass_SystHDown =
      dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystHDown->Sumw2();

  Name = "Mass_Flip";
  tuple->Mass_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->Mass_Flip->Sumw2();
  Name = "MassTOF_Flip";
  tuple->MassTOF_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_Flip->Sumw2();
  Name = "MassComb_Flip";
  tuple->MassComb_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_Flip->Sumw2();

  if (SkipSelectionPlot)
    return;

  Name = "Gen_DecayLength";
  tuple->Gen_DecayLength = dir.make<TH1F>(Name.c_str(), Name.c_str(), 1000, 0, 1000);
  tuple->Gen_DecayLength->Sumw2();
  Name = "Beta_Gen";
  tuple->Beta_Gen = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->Beta_Gen->Sumw2();
  Name = "Beta_GenChaged";
  tuple->Beta_GenCharged = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 1);
  tuple->Beta_GenCharged->Sumw2();
  Name = "Beta_Triggered";
  tuple->Beta_Triggered = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 1);
  tuple->Beta_Triggered->Sumw2();

  Name = "Beta_Matched";
  tuple->Beta_Matched = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 1);
  tuple->Beta_Matched->Sumw2();
  Name = "Beta_PreselectedA";
  tuple->Beta_PreselectedA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 1);
  tuple->Beta_PreselectedA->Sumw2();
  Name = "Beta_PreselectedB";
  tuple->Beta_PreselectedB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 1);
  tuple->Beta_PreselectedB->Sumw2();
  Name = "Beta_PreselectedC";
  tuple->Beta_PreselectedC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 1);
  tuple->Beta_PreselectedC->Sumw2();
  Name = "Beta_SelectedP";
  tuple->Beta_SelectedP = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 20, 0, 1);
  tuple->Beta_SelectedP->Sumw2();
  Name = "Beta_SelectedI";
  tuple->Beta_SelectedI = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 20, 0, 1);
  tuple->Beta_SelectedI->Sumw2();
  Name = "Beta_SelectedT";
  tuple->Beta_SelectedT = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 20, 0, 1);
  tuple->Beta_SelectedT->Sumw2();


  Name = "BefPreS_GenPtVsdRMinGen";
  tuple->BefPreS_GenPtVsdRMinGen = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0., 1.);
  tuple->BefPreS_GenPtVsdRMinGen->Sumw2();
  tuple->BefPreS_GendRMin = dir.make<TH1F>("BefPreS_GendRMin",";dR_min;Gen candidate",100,0.,5.);
  tuple->BefPreS_GendRMin->Sumw2();
  tuple->BefPreS_GenPtVsdRMinGenPostCut = dir.make<TH2F>("BefPreS_GenPtVsdRMinGenPostCut", ";GenPt (GeV);dRMinGen (after cut)", 50, 0, PtHistoUpperBound, 50, 0., 0.05);
  tuple->BefPreS_GenPtVsdRMinGenPostCut->Sumw2();
  Name = "BefPreS_GenPtVsGenMinPt";
  tuple->BefPreS_GenPtVsGenMinPt = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0, 1.);
  tuple->BefPreS_GenPtVsGenMinPt->Sumw2();
  Name = "BefPreS_GenPtVsRecoPt"; 
  tuple->BefPreS_GenPtVsRecoPt = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 50, 0, PtHistoUpperBound);
  tuple->BefPreS_GenPtVsRecoPt->Sumw2();

  Name = "BefPreS_pfType";
  tuple->BefPreS_pfType = dir.make<TH1F>(Name.c_str(), Name.c_str(), 9, 0, 9);
  tuple->BefPreS_pfType->Sumw2();

  Name = "BefPreS_massT";
  tuple->BefPreS_massT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0);
  tuple->BefPreS_massT->Sumw2();

  tuple->BefPreS_MiniRelIsoAll = dir.make<TH1F>("BefPreS_MiniRelIsoAll",";MiniRelIsoAll;Tracks/bin", 150, 0.0, 1.5);
  tuple->BefPreS_MiniRelIsoAll->Sumw2();
  tuple->BefPreS_MiniRelIsoChg = dir.make<TH1F>("BefPreS_MiniRelIsoChg",";MiniRelIsoChg;Tracks/bin",  150, 0.0, 1.5);
  tuple->BefPreS_MiniRelIsoChg->Sumw2();

  Name = "BefPreS_RecoPFMET";
  tuple->BefPreS_RecoPFMET = dir.make<TH1F>(Name.c_str(), Name.c_str(),  200, 0.0, 2000.0);
  tuple->BefPreS_RecoPFMET->Sumw2();
  Name = "BefPreS_RecoPFHT";
  tuple->BefPreS_RecoPFHT = dir.make<TH1F>(Name.c_str(), Name.c_str(),  200, 0.0, 2000.0);
  tuple->BefPreS_RecoPFHT->Sumw2();
  tuple->BefPreS_CaloNumJets = dir.make<TH1F>("BefPreS_CaloNumJets", ";Number of calo jets;Jets/bin",  70, 0.0, 70.0);
  tuple->BefPreS_CaloNumJets->Sumw2();

  Name = "BefPreS_Chi2oNdof";
  tuple->BefPreS_Chi2oNdof = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  tuple->BefPreS_Chi2oNdof->Sumw2();
  Name = "BefPreS_Qual";
  tuple->BefPreS_Qual = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  tuple->BefPreS_Qual->Sumw2();
  Name = "BefPreS_TNOH_PUA";
  tuple->BefPreS_TNOH_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->BefPreS_TNOH_PUA->Sumw2();
  Name = "BefPreS_TNOH_PUB";
  tuple->BefPreS_TNOH_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->BefPreS_TNOH_PUB->Sumw2();
  Name = "BefPreS_TNOHFraction";
  tuple->BefPreS_TNOHFraction = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 1);
  tuple->BefPreS_TNOHFraction->Sumw2();
  Name = "BefPreS_TNOPH";
  tuple->BefPreS_TNOPH = dir.make<TH1F>(Name.c_str(), Name.c_str(), 8, 0, 8.);
  tuple->BefPreS_TNOPH->Sumw2();
  Name = "BefPreS_TNOHFractionTillLast";
  tuple->BefPreS_TNOHFractionTillLast = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 1);
  tuple->BefPreS_TNOHFractionTillLast->Sumw2();
  Name = "BefPreS_TNOMHTillLast";
  tuple->BefPreS_TNOMHTillLast = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  tuple->BefPreS_TNOMHTillLast->Sumw2();
  Name = "BefPreS_Eta";
  tuple->BefPreS_Eta = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -2.6, 2.6);
  tuple->BefPreS_Eta->Sumw2();
  Name = "BefPreS_TNOM";
  tuple->BefPreS_TNOM = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->BefPreS_TNOM->Sumw2();
  Name = "BefPreS_TNOM_PUA";
  tuple->BefPreS_TNOM_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->BefPreS_TNOM_PUA->Sumw2();
  Name = "BefPreS_TNOM_PUB";
  tuple->BefPreS_TNOM_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->BefPreS_TNOM_PUB->Sumw2();
  Name = "BefPreS_nDof";
  tuple->BefPreS_nDof = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->BefPreS_nDof->Sumw2();
  Name = "BefPreS_TOFError";
  tuple->BefPreS_TOFError = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 0.25);
  tuple->BefPreS_TOFError->Sumw2();
  Name = "BefPreS_PtErrOverPt";
  tuple->BefPreS_PtErrOverPt = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 1);
  tuple->BefPreS_PtErrOverPt->Sumw2();
  Name = "BefPreS_PtErrOverPt2";
  tuple->BefPreS_PtErrOverPt2 = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 0.003);
  tuple->BefPreS_PtErrOverPt2->Sumw2();
  Name = "BefPreS_Pt";
  tuple->BefPreS_Pt = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt->Sumw2();
  Name = "BefPreS_Ias";
  tuple->BefPreS_Ias = dir.make<TH1F>(Name.c_str(), Name.c_str(), 10, 0., 1.);
  tuple->BefPreS_Ias->Sumw2();
  Name = "BefPreS_Ih";
  tuple->BefPreS_Ih = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih->Sumw2();
  Name = "BefPreS_MTOF";
  tuple->BefPreS_MTOF = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -2, 5);
  tuple->BefPreS_MTOF->Sumw2();
  Name = "BefPreS_TIsol";
  tuple->BefPreS_TIsol = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 100);
  tuple->BefPreS_TIsol->Sumw2();
  Name = "BefPreS_EoP";
  tuple->BefPreS_EoP = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 1.5);
  tuple->BefPreS_EoP->Sumw2();
  Name = "BefPreS_SumpTOverpT";
  tuple->BefPreS_SumpTOverpT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 80, 0.0, 2.0);
  tuple->BefPreS_SumpTOverpT->Sumw2();
  Name = "BefPreS_LastHitDXY";
  tuple->BefPreS_LastHitDXY = dir.make<TH1F>(Name.c_str(), Name.c_str(), 75, 0, 150);
  tuple->BefPreS_LastHitDXY->Sumw2();
  Name = "BefPreS_LastHitD3D";
  tuple->BefPreS_LastHitD3D = dir.make<TH1F>(Name.c_str(), Name.c_str(), 175, 0, 350);
  tuple->BefPreS_LastHitD3D->Sumw2();
  Name = "BefPreS_P";
  tuple->BefPreS_P = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_P->Sumw2();
  Name = "BefPreS_Pt";
  tuple->BefPreS_Pt = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt->Sumw2();
  Name = "BefPreS_Pt_PUA";
  tuple->BefPreS_Pt_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_PUA->Sumw2();
  Name = "BefPreS_Pt_PUB";
  tuple->BefPreS_Pt_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_PUB->Sumw2();
  Name = "BefPreS_Pt_Cosmic";
  tuple->BefPreS_Pt_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_Cosmic->Sumw2();
  Name = "BefPreS_Pt_DT";
  tuple->BefPreS_Pt_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_DT->Sumw2();
  Name = "BefPreS_Pt_CSC";
  tuple->BefPreS_Pt_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_CSC->Sumw2();
  Name = "BefPreS_Ias";
  tuple->BefPreS_Ias = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ias->Sumw2();
  Name = "BefPreS_Ias_PUA";
  tuple->BefPreS_Ias_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ias_PUA->Sumw2();
  Name = "BefPreS_Ias_PUB";
  tuple->BefPreS_Ias_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ias_PUB->Sumw2();
  Name = "BefPreS_Ias_Cosmic";
  tuple->BefPreS_Ias_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ias_Cosmic->Sumw2();
  Name = "BefPreS_Ih_Cosmic";
  tuple->BefPreS_Ih_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih_Cosmic->Sumw2();
  Name = "BefPreS_Ih";
  tuple->BefPreS_Ih = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih->Sumw2();
  Name = "BefPreS_Ih_PUA";
  tuple->BefPreS_Ih_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih_PUA->Sumw2();
  Name = "BefPreS_Ih_PUB";
  tuple->BefPreS_Ih_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih_PUB->Sumw2();
  Name = "BefPreS_TOF";
  tuple->BefPreS_TOF = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF->Sumw2();
  Name = "BefPreS_TOF_PUA";
  tuple->BefPreS_TOF_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF_PUA->Sumw2();
  Name = "BefPreS_TOF_PUB";
  tuple->BefPreS_TOF_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF_PUB->Sumw2();
  Name = "BefPreS_TOF_DT";
  tuple->BefPreS_TOF_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF_DT->Sumw2();
  Name = "BefPreS_TOF_CSC";
  tuple->BefPreS_TOF_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF_CSC->Sumw2();
  Name = "BefPreS_dR_NVTrack";
  tuple->BefPreS_dR_NVTrack = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 1);
  tuple->BefPreS_dR_NVTrack->Sumw2();
  Name = "BefPreS_MatchedStations";
  tuple->BefPreS_MatchedStations = dir.make<TH1F>(Name.c_str(), Name.c_str(), 8, -0.5, 7.5);
  tuple->BefPreS_MatchedStations->Sumw2();
  Name = "BefPreS_InnerInvPtDiff";
  tuple->BefPreS_InnerInvPtDiff = dir.make<TH1F>(Name.c_str(), Name.c_str(), 120, -4, 4);
  tuple->BefPreS_InnerInvPtDiff->Sumw2();
  Name = "BefPreS_Phi";
  tuple->BefPreS_Phi = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -3.14, 3.14);
  tuple->BefPreS_Phi->Sumw2();
  Name = "BefPreS_TimeAtIP";
  tuple->BefPreS_TimeAtIP = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -100, 100);
  tuple->BefPreS_TimeAtIP->Sumw2();
  Name = "BefPreS_OpenAngle";
  tuple->BefPreS_OpenAngle = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -0.3, 3.15);
  tuple->BefPreS_OpenAngle->Sumw2();
  Name = "BefPreS_OpenAngle_Cosmic";
  tuple->BefPreS_OpenAngle_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -0.3, 3.15);
  tuple->BefPreS_OpenAngle_Cosmic->Sumw2();

  Name = "BefPreS_NVertex";
  tuple->BefPreS_NVertex = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 50);
  tuple->BefPreS_NVertex->Sumw2();
  Name = "BefPreS_NVertex_NoEventWeight";
  tuple->BefPreS_NVertex_NoEventWeight = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 50);
  tuple->BefPreS_NVertex_NoEventWeight->Sumw2();
  Name = "BefPreS_PV";
  tuple->BefPreS_PV = dir.make<TH1F>(Name.c_str(), Name.c_str(), 60, 0, 60);
  tuple->BefPreS_PV->Sumw2();
  Name = "BefPreS_PV_NoEventWeight";
  tuple->BefPreS_PV_NoEventWeight = dir.make<TH1F>(Name.c_str(), Name.c_str(), 60, 0, 60);
  tuple->BefPreS_PV_NoEventWeight->Sumw2();
  tuple->BefPreS_NOMoNOH = dir.make<TH1F>("BefPreS_NOMoNOH",";Num of measurment / num of hits;Tracks/bin",10,0.,1.0);
  tuple->BefPreS_NOMoNOH->Sumw2();
  Name = "BefPreS_NOMoNOHvsPV";
  tuple->BefPreS_NOMoNOHvsPV = dir.make<TProfile>(Name.c_str(), Name.c_str(), 60, 0, 60);
  tuple->BefPreS_NOMoNOHvsPV->Sumw2();
  Name = "BefPreS_dzAll";
  tuple->BefPreS_dzAll = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -10, 10);
  tuple->BefPreS_dzAll->Sumw2();
  Name = "BefPreS_dxyAll";
  tuple->BefPreS_dxyAll = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.2, 0.2);
  tuple->BefPreS_dxyAll->Sumw2();
  Name = "BefPreS_Dz";
  tuple->BefPreS_Dz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.2, 0.2);
  tuple->BefPreS_Dz->Sumw2();
  Name = "BefPreS_Dxy";
  tuple->BefPreS_Dxy = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.2, 0.2);
  tuple->BefPreS_Dxy->Sumw2();

  Name = "BefPreS_SegSep";
  tuple->BefPreS_SegSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 2.5);
  tuple->BefPreS_SegSep->Sumw2();
  Name = "BefPreS_SegMinEtaSep";
  tuple->BefPreS_SegMinEtaSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -1., 1.);
  tuple->BefPreS_SegMinEtaSep->Sumw2();
  Name = "BefPreS_SegMinPhiSep";
  tuple->BefPreS_SegMinPhiSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -3.3, 3.3);
  tuple->BefPreS_SegMinPhiSep->Sumw2();
  Name = "BefPreS_SegMinEtaSep_FailDz";
  tuple->BefPreS_SegMinEtaSep_FailDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -1., 1.);
  tuple->BefPreS_SegMinEtaSep_FailDz->Sumw2();
  Name = "BefPreS_SegMinEtaSep_PassDz";
  tuple->BefPreS_SegMinEtaSep_PassDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -1., 1.);
  tuple->BefPreS_SegMinEtaSep_PassDz->Sumw2();
  Name = "BefPreS_Dz_FailSep";
  tuple->BefPreS_Dz_FailSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -150, 150);
  tuple->BefPreS_Dz_FailSep->Sumw2();

  Name = "BefPreS_Dxy";
  tuple->BefPreS_Dxy = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.1, 0.1);
  tuple->BefPreS_Dxy->Sumw2();
  Name = "BefPreS_Dxy_Cosmic";
  tuple->BefPreS_Dxy_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->BefPreS_Dxy_Cosmic->Sumw2();
  Name = "BefPreS_Dz";
  tuple->BefPreS_Dz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.1, 0.1);
  tuple->BefPreS_Dz->Sumw2();
  Name = "BefPreS_Dz_Cosmic";
  tuple->BefPreS_Dz_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->BefPreS_Dz_Cosmic->Sumw2();
  Name = "BefPreS_Dz_CSC";
  tuple->BefPreS_Dz_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->BefPreS_Dz_CSC->Sumw2();
  Name = "BefPreS_Dz_DT";
  tuple->BefPreS_Dz_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->BefPreS_Dz_DT->Sumw2();
  Name = "BefPreS_Pt_FailDz";
  tuple->BefPreS_Pt_FailDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_FailDz->Sumw2();
  Name = "BefPreS_Pt_FailDz_DT";
  tuple->BefPreS_Pt_FailDz_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_FailDz_DT->Sumw2();
  Name = "BefPreS_Pt_FailDz_CSC";
  tuple->BefPreS_Pt_FailDz_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_FailDz_CSC->Sumw2();
  Name = "BefPreS_TOF_FailDz";
  tuple->BefPreS_TOF_FailDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF_FailDz->Sumw2();
  Name = "BefPreS_TOF_FailDz_DT";
  tuple->BefPreS_TOF_FailDz_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF_FailDz_DT->Sumw2();
  Name = "BefPreS_TOF_FailDz_CSC";
  tuple->BefPreS_TOF_FailDz_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->BefPreS_TOF_FailDz_CSC->Sumw2();
  Name = "BefPreS_GenPtVsRecoPt";
  tuple->BefPreS_GenPtVsRecoPt = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 50, 0, PtHistoUpperBound);
  tuple->BefPreS_GenPtVsRecoPt->Sumw2();
  Name = "BefPreS_PtErrOverPtVsPtErrOverPt2";
  tuple->BefPreS_PtErrOverPtVsPtErrOverPt2 = dir.make<TH2F>(Name.c_str(), Name.c_str(),  40, 0., 1., 40, 0., 0.003);
  tuple->BefPreS_PtErrOverPtVsPtErrOverPt2->Sumw2();
  Name = "BefPreS_PtErrOverPtVsPt";
  tuple->BefPreS_PtErrOverPtVsPt = dir.make<TH2F>(Name.c_str(), Name.c_str(),  40, 0., 1., 40, 0., 4000);
  tuple->BefPreS_PtErrOverPtVsPt->Sumw2();
  
  Name = "BefPreS_ProbQ";
  tuple->BefPreS_ProbQ = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->BefPreS_ProbQ->Sumw2();
  Name = "BefPreS_ProbXY";
  tuple->BefPreS_ProbXY = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->BefPreS_ProbXY->Sumw2();
  Name = "BefPreS_ProbQNoL1";
  tuple->BefPreS_ProbQNoL1 = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->BefPreS_ProbQNoL1->Sumw2();
  Name = "BefPreS_ProbXYNoL1";
  tuple->BefPreS_ProbXYNoL1 = dir.make<TH1F>(Name.c_str(), Name.c_str(), 100, 0, 1);
  tuple->BefPreS_ProbXYNoL1->Sumw2();
  Name = "BefPreS_MassErr";
  tuple->BefPreS_MassErr = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0., 10.);
  tuple->BefPreS_MassErr->Sumw2();
  Name = "BefPreS_ProbQVsIas";
  tuple->BefPreS_ProbQVsIas = dir.make<TH2F>(Name.c_str(), Name.c_str(), 100, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->BefPreS_ProbQVsIas->Sumw2();

  tuple->BefPreS_CluProbQVsPixelLayer = dir.make<TH2F>("BefPreS_CluProbQVsPixelLayer",";CluProbQ;Layer",20,0.,1.,4,0.,4.);
  tuple->BefPreS_CluProbQVsPixelLayer->Sumw2();
  tuple->BefPreS_CluProbXYVsPixelLayer = dir.make<TH2F>("BefPreS_CluProbXYVsPixelLayer",";CluProbXY;Layer",100,0.,1.,4,0.,4.);
  tuple->BefPreS_CluProbXYVsPixelLayer->Sumw2();
  tuple->BefPreS_CluNormChargeVsPixelLayer = dir.make<TH2F>("BefPreS_CluNormChargeVsPixelLayer",";CluNormCharge;Layer",5,0.,600.,4,0.,4.);
  tuple->BefPreS_CluNormChargeVsPixelLayer->Sumw2();
  tuple->BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma",";CluNormCharge;Layer",5,0.,600.,4,0.,4.);
  tuple->BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma->Sumw2();
  tuple->BefPreS_CluSizeVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeVsPixelLayer",";CluSize;Layer",10,0.,10.,4,0.,4.);
  tuple->BefPreS_CluSizeVsPixelLayer->Sumw2();
  tuple->BefPreS_CluSizeXVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeXVsPixelLayer",";CluSizeX;Layer",10,0.,10.,4,0.,4.);
  tuple->BefPreS_CluSizeXVsPixelLayer->Sumw2();
  tuple->BefPreS_CluSizeYVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeYVsPixelLayer",";CluSizeY;Layer",10,0.,10.,4,0.,4.);
  tuple->BefPreS_CluSizeYVsPixelLayer->Sumw2();
  tuple->BefPreS_CluSpecInCPEVsPixelLayer = dir.make<TH2F>("BefPreS_CluSpecInCPEVsPixelLayer",";CluSpecInCPE;Layer",4,0.,4.,4,0.,4.);
  tuple->BefPreS_CluSpecInCPEVsPixelLayer->Sumw2();

  tuple->BefPreS_CluCotBetaVsPixelLayer_lowProbXY = dir.make<TH2F>("BefPreS_CluCotBetaVsPixelLayer_lowProbXY",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->BefPreS_CluCotBetaVsPixelLayer_lowProbXY->Sumw2();
  tuple->BefPreS_CluCotAlphaVsPixelLayer_lowProbXY = dir.make<TH2F>("BefPreS_CluCotAlphaVsPixelLayer_lowProbXY",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
  tuple->BefPreS_CluCotAlphaVsPixelLayer_lowProbXY->Sumw2();
  tuple->BefPreS_CluCotBetaVsPixelLayer = dir.make<TH2F>("BefPreS_CluCotBetaVsPixelLayer",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->BefPreS_CluCotBetaVsPixelLayer->Sumw2();
  tuple->BefPreS_CluCotAlphaVsPixelLayer = dir.make<TH2F>("BefPreS_CluCotAlphaVsPixelLayer",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
  tuple->BefPreS_CluCotAlphaVsPixelLayer->Sumw2();

  tuple->BefPreS_CluNormChargeVsStripLayer_lowBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_lowBetaGamma",";CluNormCharge;Layer",5,0.,600.,20,0.,20.);
  tuple->BefPreS_CluNormChargeVsStripLayer_lowBetaGamma->Sumw2();

  tuple->BefPreS_dRMinPfJet= dir.make<TH1F>("BefPreS_dRMinPfJet",";dRMinPfJet",100,0.,1.5);
  tuple->BefPreS_dRMinPfJet->Sumw2();
  tuple->BefPreS_dRMinPfJetVsIas =  dir.make<TH2F>("BefPreS_dRMinPfJetVsIas",";dRMinPfJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->BefPreS_dRMinPfJetVsIas->Sumw2();
  tuple->BefPreS_dRMinCaloJet= dir.make<TH1F>("BefPreS_dRMinCaloJet",";dRMinCaloJet",100,0.,1.5);
  tuple->BefPreS_dRMinCaloJet->Sumw2();
  tuple->BefPreS_dRMinCaloJetVsIas =  dir.make<TH2F>("BefPreS_dRMinCaloJetVsIas",";dRMinCaloJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->BefPreS_dRMinCaloJetVsIas->Sumw2();
  tuple->BefPreS_genGammaBetaVsProbXYNoL1 =  dir.make<TH2F>("BefPreS_genGammaBetaVsProbXYNoL1", ";#gamma #beta;ProbXYNoL1",10,0.,1.3,20,0.,1.);
  tuple->BefPreS_genGammaBetaVsProbXYNoL1->Sumw2();
  tuple->BefPreS_dRVsPtPfJet = dir.make<TH2F>("BefPreS_dRVsPtPfJet",";dR(cand,jet);p_{T}",100,0.,1.5,100,0.,1000.);
  tuple->BefPreS_dRVsPtPfJet->Sumw2();
  tuple->BefPreS_dRVsdPtPfCaloJet = dir.make<TH2F>("BefPreS_dRVsdPtPfCaloJet",";dRmin;dPtPfCaloJet",100,0.,1.5,20,0.,100.);
  tuple->BefPreS_dRVsdPtPfCaloJet->Sumw2();
  

  Name = "PostPreS_TriggerType";
  tuple->PostPreS_TriggerType = dir.make<TH1F>(Name.c_str(), Name.c_str(), 3, 0., 3.);
  tuple->PostPreS_TriggerType->Sumw2();

  Name = "PostPreS_pfType";
  tuple->PostPreS_pfType = dir.make<TH1F>(Name.c_str(), Name.c_str(), 9, 0, 9);
  tuple->PostPreS_pfType->Sumw2();
  tuple->PostPreS_pfTypeVsIas = dir.make<TH2F>("PostPreS_pfTypeVsIas","PostPreS_pfTypeVsIas", 9, 0, 9.,20,0.,1.);
  tuple->PostPreS_pfTypeVsIas->Sumw2();

  Name = "PostPreS_massT";
  tuple->PostPreS_massT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0);
  tuple->PostPreS_massT->Sumw2();
  tuple->PostPreS_massTVsIas = dir.make<TH2F>("PostPreS_massTVsIas","PostPreS_massTVsIas",50, 0.0, 250.0, 20, 0., 1.);
  tuple->PostPreS_massTVsIas->Sumw2();
  
  tuple->PostPreS_MiniRelIsoAll = dir.make<TH1F>("PostPreS_MiniRelIsoAll",";MiniRelIsoAll;Tracks/bin", 150, 0.0, 1.5);
  tuple->PostPreS_MiniRelIsoAll->Sumw2();
  tuple->PostPreS_MiniRelIsoAllVsIas =  dir.make<TH2F>("PostPreS_MiniRelIsoAllVsIas","PostPreS_MiniRelIsoAllVsIas", 150, 0.0, 1.5, 10, 0.,1.);
  tuple->PostPreS_MiniRelIsoAllVsIas->Sumw2();
  tuple->PostPreS_MiniRelIsoChg = dir.make<TH1F>("PostPreS_MiniRelIsoChg",";MiniRelIsoChg;Tracks/bin",  150, 0.0, 1.5);
  tuple->PostPreS_MiniRelIsoChg->Sumw2();
  
  Name = "PostPreS_RecoPFMET";
  tuple->PostPreS_RecoPFMET = dir.make<TH1F>(Name.c_str(), Name.c_str(),  200, 0.0, 2000.0);
  tuple->PostPreS_RecoPFMET->Sumw2();
  Name = "PostPreS_RecoPFHT";
  tuple->PostPreS_RecoPFHT = dir.make<TH1F>(Name.c_str(), Name.c_str(),  200, 0.0, 2000.0);
  tuple->PostPreS_RecoPFHT->Sumw2();
  tuple->PostPreS_CaloNumJets = dir.make<TH1F>("PostPreS_CaloNumJets",";Number of calo jets;Jets/bin",  30, 0.0, 30.0);
  tuple->PostPreS_CaloNumJets->Sumw2();
  
  Name = "PostPreS_Chi2oNdof";
  tuple->PostPreS_Chi2oNdof = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  tuple->PostPreS_Chi2oNdof->Sumw2();
  tuple->PostPreS_Chi2oNdofVsIas = dir.make<TH2F>("PostPreS_Chi2oNdofVsIas","PostPreS_Chi2oNdofVsIas",20, 0, 20,10,0.,1.);
  tuple->PostPreS_Chi2oNdofVsIas->Sumw2();
  Name = "PostPreS_Qual";
  tuple->PostPreS_Qual = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  tuple->PostPreS_Qual->Sumw2();
  Name = "PostPreS_TNOH_PUA";
  tuple->PostPreS_TNOH_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->PostPreS_TNOH_PUA->Sumw2();
  Name = "PostPreS_TNOH_PUB";
  tuple->PostPreS_TNOH_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->PostPreS_TNOH_PUB->Sumw2();
  Name = "PostPreS_TNOHFraction";
  tuple->PostPreS_TNOHFraction = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 1);
  tuple->PostPreS_TNOHFraction->Sumw2();
  tuple->PostPreS_TNOHFractionVsIas = dir.make<TH2F>("PostPreS_TNOHFractionVsIas","PostPreS_TNOHFractionVsIas",50, 0, 1,10,0.,1.);
  tuple->PostPreS_TNOHFractionVsIas->Sumw2();
  Name = "PostPreS_TNOPH";
  tuple->PostPreS_TNOPH = dir.make<TH1F>(Name.c_str(), Name.c_str(), 8, 0, 8);
  tuple->PostPreS_TNOPH->Sumw2();
  tuple->PostPreS_TNOPHVsIas = dir.make<TH2F>("PostPreS_TNOPHVsIas","PostPreS_TNOPHVsIas", 8, 0, 8,20,0.,1.);
  Name = "PostPreS_TNOHFractionTillLast";
  tuple->PostPreS_TNOHFractionTillLast = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 1);
  tuple->PostPreS_TNOHFractionTillLast->Sumw2();
  Name = "PostPreS_TNOMHTillLast";
  tuple->PostPreS_TNOMHTillLast = dir.make<TH1F>(Name.c_str(), Name.c_str(), 20, 0, 20);
  tuple->PostPreS_TNOMHTillLast->Sumw2();
  Name = "PostPreS_Eta";
  tuple->PostPreS_Eta = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -2.6, 2.6);
  tuple->PostPreS_Eta->Sumw2();
  tuple->PostPreS_EtaVsIas =  dir.make<TH2F>("PostPreS_EtaVsIas","PostPreS_EtaVsIas", 50, -2.6, 2.6, 20,0.,1.);
  tuple->PostPreS_EtaVsIas->Sumw2();
  Name = "PostPreS_TNOM";
  tuple->PostPreS_TNOM = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->PostPreS_TNOM->Sumw2();
  tuple->PostPreS_TNOMVsIas = dir.make<TH2F>("PostPreS_TNOMVsIas","PostPreS_TNOMVsIas", 40, 0, 40,20,0.,1.);
  tuple->PostPreS_TNOMVsIas->Sumw2();
  Name = "PostPreS_TNOM_PUA";
  tuple->PostPreS_TNOM_PUA = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->PostPreS_TNOM_PUA->Sumw2();
  Name = "PostPreS_TNOM_PUB";
  tuple->PostPreS_TNOM_PUB = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->PostPreS_TNOM_PUB->Sumw2();
  Name = "PostPreS_nDof";
  tuple->PostPreS_nDof = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 40);
  tuple->PostPreS_nDof->Sumw2();
  Name = "PostPreS_TOFError";
  tuple->PostPreS_TOFError = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 0.25);
  tuple->PostPreS_TOFError->Sumw2();
  Name = "PostPreS_PtErrOverPt";
  tuple->PostPreS_PtErrOverPt = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 1);
  tuple->PostPreS_PtErrOverPt->Sumw2();
  tuple->PostPreS_PtErrOverPtVsIas =  dir.make<TH2F>("PostPreS_PtErrOverPtVsIas","PostPreS_PtErrOverPtVsIas", 40, 0, 1, 20, 0.,1.);
  Name = "PostPreS_PtErrOverPt2";
  tuple->PostPreS_PtErrOverPt2 = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 0.003);
  tuple->PostPreS_PtErrOverPt2->Sumw2();
  Name = "PostPreS_Pt";
  tuple->PostPreS_Pt = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->PostPreS_Pt->Sumw2();
  tuple->PostPreS_PtVsIas = dir.make<TH2F>("PostPreS_PtVsIas","PostPreS_PtVsIas", 50, 0, PtHistoUpperBound, 20, 0., 1.);
  tuple->PostPreS_PtVsIas->Sumw2();
  Name = "PostPreS_Ias";
  tuple->PostPreS_Ias = dir.make<TH1F>(Name.c_str(), Name.c_str(), 10, 0, dEdxS_UpLim);
  tuple->PostPreS_Ias->Sumw2();
  Name = "PostPreS_Ias_NoEventWeight";
  tuple->PostPreS_Ias_NoEventWeight = dir.make<TH1F>(Name.c_str(), Name.c_str(), 10, 0, dEdxS_UpLim);
  tuple->PostPreS_Ias_NoEventWeight->Sumw2();
  Name = "PostPreS_Ih";
  tuple->PostPreS_Ih = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  tuple->PostPreS_Ih->Sumw2();
  tuple->PostPreS_IhVsIas = dir.make<TH2F>("PostPreS_IhVsIas","PostPreS_IhVsIas",200, 0, dEdxM_UpLim, 20, 0.,1.);
  tuple->PostPreS_IhVsIas->Sumw2();
  Name = "PostPreS_Ih_NoEventWeight";
  tuple->PostPreS_Ih_NoEventWeight = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim);
  tuple->PostPreS_Ih_NoEventWeight->Sumw2();
  Name = "PostPreS_MTOF";
  tuple->PostPreS_MTOF = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -2, 5);
  tuple->PostPreS_MTOF->Sumw2();
  Name = "PostPreS_TIsol";
  tuple->PostPreS_TIsol = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 100);
  tuple->PostPreS_TIsol->Sumw2();
  tuple->PostPreS_TIsolVsIas = dir.make<TH2F>("PostPreS_TIsolVsIas","PostPreS_TIsolVsIas",25, 0, 100, 10, 0., 1.);
  tuple->PostPreS_TIsolVsIas->Sumw2();
  Name = "PostPreS_EoP";
  tuple->PostPreS_EoP = dir.make<TH1F>(Name.c_str(), Name.c_str(), 25, 0, 1.5);
  tuple->PostPreS_EoP->Sumw2();
  tuple->PostPreS_EoPVsIas = dir.make<TH2F>("PostPreS_EoPVsIas","PostPreS_EoPVsIas;E/ p;Ias",25, 0, 1.5, 20, 0.,1.);
  tuple->PostPreS_EoPVsIas->Sumw2();
  Name = "PostPreS_SumpTOverpT";
  tuple->PostPreS_SumpTOverpT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 80, 0.0, 2.0);
  tuple->PostPreS_SumpTOverpT->Sumw2();
  tuple->PostPreS_SumpTOverpTVsIas = dir.make<TH2F>("PostPreS_SumpTOverpTVsIas","PostPreS_SumpTOverpTVsIas;SumpTOverpT;Ias", 80, 0.0, 2.0, 20, 0.,1.);
  tuple->PostPreS_SumpTOverpTVsIas->Sumw2();
  Name = "PostPreS_LastHitDXY";
  tuple->PostPreS_LastHitDXY = dir.make<TH1F>(Name.c_str(), Name.c_str(), 75, 0, 150);
  tuple->PostPreS_LastHitDXY->Sumw2();
  tuple->PostPreS_LastHitDXYVsEta  = dir.make<TH2F>("PostPreS_LastHitDXYVsEta","PostPreS_LastHitDXYVsEta;LastHitDXY;Eta", 75, 0, 150, 20, 0.,1.);
  tuple->PostPreS_LastHitDXYVsEta->Sumw2();
  Name = "PostPreS_LastHitD3D";
  tuple->PostPreS_LastHitD3D = dir.make<TH1F>(Name.c_str(), Name.c_str(), 175, 0, 350);
  tuple->PostPreS_LastHitD3D->Sumw2();
  tuple->PostPreS_LastHitD3DVsEta = dir.make<TH2F>("PostPreS_LastHitD3DVsEta","PostPreS_LastHitD3DVsEta;LastHitD3D;Eta", 175, 0, 350, 20, 0.,1.);
  tuple->PostPreS_LastHitD3DVsEta->Sumw2();
  Name = "PostPreS_P";
  tuple->PostPreS_P = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->PostPreS_P->Sumw2();
  Name = "PostPreS_dR_NVTrack";
  tuple->PostPreS_dR_NVTrack = dir.make<TH1F>(Name.c_str(), Name.c_str(), 40, 0, 1);
  tuple->PostPreS_dR_NVTrack->Sumw2();
  Name = "PostPreS_MatchedStations";
  tuple->PostPreS_MatchedStations = dir.make<TH1F>(Name.c_str(), Name.c_str(), 8, -0.5, 7.5);
  tuple->PostPreS_MatchedStations->Sumw2();
  Name = "PostPreS_InnerInvPtDiff";
  tuple->PostPreS_InnerInvPtDiff = dir.make<TH1F>(Name.c_str(), Name.c_str(), 120, -4, 4);
  tuple->PostPreS_InnerInvPtDiff->Sumw2();
  Name = "PostPreS_Phi";
  tuple->PostPreS_Phi = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -3.14, 3.14);
  tuple->PostPreS_Phi->Sumw2();
  Name = "PostPreS_TimeAtIP";
  tuple->PostPreS_TimeAtIP = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -100, 100);
  tuple->PostPreS_TimeAtIP->Sumw2();
  Name = "PostPreS_OpenAngle";
  tuple->PostPreS_OpenAngle = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -0.3, 3.15);
  tuple->PostPreS_OpenAngle->Sumw2();
  Name = "PostPreS_OpenAngle_Cosmic";
  tuple->PostPreS_OpenAngle_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -0.3, 3.15);
  tuple->PostPreS_OpenAngle_Cosmic->Sumw2();
  
  Name = "PostPreS_NVertex";
  tuple->PostPreS_NVertex = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 50);
  tuple->PostPreS_NVertex->Sumw2();
  Name = "PostPreS_NVertex_NoEventWeight";
  tuple->PostPreS_NVertex_NoEventWeight = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 50);
  tuple->PostPreS_NVertex_NoEventWeight->Sumw2();
  Name = "PostPreS_PV";
  tuple->PostPreS_PV = dir.make<TH1F>(Name.c_str(), Name.c_str(), 60, 0, 60);
  tuple->PostPreS_PV->Sumw2();
  Name = "PostPreS_PV_NoEventWeight";
  tuple->PostPreS_PV_NoEventWeight = dir.make<TH1F>(Name.c_str(), Name.c_str(), 60, 0, 60);
  tuple->PostPreS_PV_NoEventWeight->Sumw2();
  tuple->PostPreS_NOMoNOH = dir.make<TH1F>("PostPreS_NOMoNOH",";Num of measurment / num of hits;Tracks/bin",10,0.,1.0);
  tuple->PostPreS_NOMoNOH->Sumw2();
  Name = "PostPreS_NOMoNOHvsPV";
  tuple->PostPreS_NOMoNOHvsPV = dir.make<TProfile>(Name.c_str(), Name.c_str(), 60, 0, 60);
  tuple->PostPreS_NOMoNOHvsPV->Sumw2();
  Name = "PostPreS_dzAll";
  tuple->PostPreS_dzAll = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -10, 10);
  tuple->PostPreS_dzAll->Sumw2();
  Name = "PostPreS_dxyAll";
  tuple->PostPreS_dxyAll = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.2, 0.2);
  tuple->PostPreS_dxyAll->Sumw2();
  Name = "PostPreS_Dz";
  tuple->PostPreS_Dz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.1, 0.1);
  tuple->PostPreS_Dz->Sumw2();
  Name = "PostPreS_Dxy";
  tuple->PostPreS_Dxy = dir.make<TH1F>(Name.c_str(), Name.c_str(), 200, -0.1, 0.1);
  tuple->PostPreS_Dxy->Sumw2();
  
  Name = "PostPreS_SegSep";
  tuple->PostPreS_SegSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, 2.5);
  tuple->PostPreS_SegSep->Sumw2();
  Name = "PostPreS_SegMinEtaSep";
  tuple->PostPreS_SegMinEtaSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -1., 1.);
  tuple->PostPreS_SegMinEtaSep->Sumw2();
  Name = "PostPreS_SegMinPhiSep";
  tuple->PostPreS_SegMinPhiSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -3.3, 3.3);
  tuple->PostPreS_SegMinPhiSep->Sumw2();
  Name = "PostPreS_SegMinEtaSep_FailDz";
  tuple->PostPreS_SegMinEtaSep_FailDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -1., 1.);
  tuple->PostPreS_SegMinEtaSep_FailDz->Sumw2();
  Name = "PostPreS_SegMinEtaSep_PassDz";
  tuple->PostPreS_SegMinEtaSep_PassDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -1., 1.);
  tuple->PostPreS_SegMinEtaSep_PassDz->Sumw2();
  Name = "PostPreS_Dz_FailSep";
  tuple->PostPreS_Dz_FailSep = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, -150, 150);
  tuple->PostPreS_Dz_FailSep->Sumw2();
  
  Name = "PostPreS_Dxy_Cosmic";
  tuple->PostPreS_Dxy_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->PostPreS_Dxy_Cosmic->Sumw2();
  Name = "PostPreS_Dz_Cosmic";
  tuple->PostPreS_Dz_Cosmic = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->PostPreS_Dz_Cosmic->Sumw2();
  Name = "PostPreS_Dz_CSC";
  tuple->PostPreS_Dz_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->PostPreS_Dz_CSC->Sumw2();
  Name = "PostPreS_Dz_DT";
  tuple->PostPreS_Dz_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -IPbound, IPbound);
  tuple->PostPreS_Dz_DT->Sumw2();
  Name = "PostPreS_Pt_FailDz";
  tuple->PostPreS_Pt_FailDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->PostPreS_Pt_FailDz->Sumw2();
  Name = "PostPreS_Pt_FailDz_DT";
  tuple->PostPreS_Pt_FailDz_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->PostPreS_Pt_FailDz_DT->Sumw2();
  Name = "PostPreS_Pt_FailDz_CSC";
  tuple->PostPreS_Pt_FailDz_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
  tuple->PostPreS_Pt_FailDz_CSC->Sumw2();
  Name = "PostPreS_TOF_FailDz";
  tuple->PostPreS_TOF_FailDz = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->PostPreS_TOF_FailDz->Sumw2();
  Name = "PostPreS_TOF_FailDz_DT";
  tuple->PostPreS_TOF_FailDz_DT = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->PostPreS_TOF_FailDz_DT->Sumw2();
  Name = "PostPreS_TOF_FailDz_CSC";
  tuple->PostPreS_TOF_FailDz_CSC = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
  tuple->PostPreS_TOF_FailDz_CSC->Sumw2();
  Name = "PostPreS_PtErrOverPtVsPtErrOverPt2";
  tuple->PostPreS_PtErrOverPtVsPtErrOverPt2 = dir.make<TH2F>(Name.c_str(), Name.c_str(),  40, 0., 1., 40, 0., 0.003);
  tuple->PostPreS_PtErrOverPtVsPtErrOverPt2->Sumw2();
  Name = "PostPreS_PtErrOverPtVsPt";
  tuple->PostPreS_PtErrOverPtVsPt = dir.make<TH2F>(Name.c_str(), Name.c_str(),  40, 0., 1., 40, 0., 4000);
  tuple->PostPreS_PtErrOverPtVsPt->Sumw2();
  Name = "PostPreS_GenPtVsRecoPt"; 
  tuple->PostPreS_GenPtVsRecoPt = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 50, 0, PtHistoUpperBound);
  tuple->PostPreS_GenPtVsRecoPt->Sumw2();
  
  tuple->PostPreS_ProbQ = dir.make<TH1F>("PostPreS_ProbQ", ";ProbQ", 20, 0., 1.);
  tuple->PostPreS_ProbQ->Sumw2();
  tuple->PostPreS_ProbQVsIas = dir.make<TH2F>("PostPreS_ProbQVsIas",";ProbQ;I_{as}", 20, 0., 1., 20, 0., 1.);
  tuple->PostPreS_ProbQVsIas->Sumw2();
  tuple->PostPreS_ProbXY = dir.make<TH1F>("PostPreS_ProbXY", ";ProbXY", 100, 0, 1);
  tuple->PostPreS_ProbXY->Sumw2();
  tuple->PostPreS_ProbXY_highIas = dir.make<TH1F>("PostPreS_ProbXY_highIas",";ProbXY (I_{as} > 0.6)", 100, 0, 1);
  tuple->PostPreS_ProbXY_highIas->Sumw2();
  tuple->PostPreS_ProbXYVsIas = dir.make<TH2F>("PostPreS_ProbXYVsIas",";ProbXY;Ias",  100, 0, 1, 10, 0., 1.);
  tuple->PostPreS_ProbXYVsIas->Sumw2();
  tuple->PostPreS_ProbXYVsIas_highIas = dir.make<TH2F>("PostPreS_ProbXYVsIas_highIas",";ProbXY (I_{as} > 0.6);Ias (I_{as} > 0.6)",  100, 0, 1, 10, 0., 1.);
  tuple->PostPreS_ProbXYVsIas_highIas->Sumw2();
  tuple->PostPreS_ProbXYVsProbQ = dir.make<TH2F>("PostPreS_ProbXYVsProbQ",";ProbXY;ProbQ",  100, 0., 1., 10, 0., 1.);
  tuple->PostPreS_ProbXYVsProbQ->Sumw2();
  tuple->PostPreS_ProbXYVsProbQ_highIas = dir.make<TH2F>("PostPreS_ProbXYVsProbQ_highIas",";ProbXY (I_{as} > 0.6);ProbQ (I_{as} > 0.6)",  100, 0., 1., 10, 0., 1.);
  tuple->PostPreS_ProbXYVsProbQ_highIas->Sumw2();
  tuple->PostPreS_ProbQNoL1 = dir.make<TH1F>("PostPreS_ProbQNoL1",";ProbQNoL1", 20, 0., 1.);
  tuple->PostPreS_ProbQNoL1->Sumw2();
  tuple->PostPreS_ProbQNoL1VsIas = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas","PostPreS_ProbQNoL1VsIas;PostPreS_ProbQNoL1VsIas;PostPreS_ProbQNoL1VsIas",20, 0., 1., 20, 0., 1.);
  tuple->PostPreS_ProbQNoL1VsIas->Sumw2();
  tuple->PostPreS_ProbXYNoL1 = dir.make<TH1F>("PostPreS_ProbXYNoL1", ";ProbXYNoL1", 100, 0, 1);
  tuple->PostPreS_ProbXYNoL1->Sumw2();
  tuple->PostPreS_ProbXYNoL1_highIas = dir.make<TH1F>("PostPreS_ProbXYNoL1_highIas", ";ProbXYNoL1 (I_{as} > 0.6)", 100, 0, 1);
  tuple->PostPreS_ProbXYNoL1_highIas->Sumw2();
  tuple->PostPreS_ProbXYNoL1VsIas = dir.make<TH2F>("PostPreS_ProbXYNoL1VsIas",";ProbXYNoL1;Ias", 100, 0., 1., 20, 0.,1.);
  tuple->PostPreS_ProbXYNoL1VsIas->Sumw2();
  tuple->PostPreS_ProbXYNoL1VsIas_highIas = dir.make<TH2F>("PostPreS_ProbXYNoL1VsIas_highIas",";ProbXYNoL1 (I_{as} > 0.6);Ias (I_{as} > 0.6)", 100, 0., 1., 20, 0.,1.);
  tuple->PostPreS_ProbXYNoL1VsIas_highIas->Sumw2();
  tuple->PostPreS_ProbXYNoL1VsProbQNoL1  = dir.make<TH2F>("PostPreS_ProbXYNoL1VsProbQNoL1",";ProbXYNoL1;ProbQNoL1", 100, 0., 1., 20,0.,1.);
  tuple->PostPreS_ProbXYNoL1VsProbQNoL1->Sumw2();
  tuple->PostPreS_ProbXYNoL1VsProbQNoL1_highIas  = dir.make<TH2F>("PostPreS_ProbXYNoL1VsProbQNoL1_highIas",";ProbXYNoL1 (I_{as} > 0.6);ProbQNoL1 (I_{as} > 0.6)", 100, 0., 1., 20,0.,1.);
  tuple->PostPreS_ProbXYNoL1VsProbQNoL1_highIas->Sumw2();

  Name = "PostPreS_MassErr";
  tuple->PostPreS_MassErr = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0., 10.);
  tuple->PostPreS_MassErr->Sumw2();
  tuple->PostPreS_MassErrVsIas = dir.make<TH2F>("PostPreS_MassErrVsIas","PostPreS_MassErrVsIas;PostPreS_MassErrVsIas;PostPreS_MassErrVsIas",50, 0., 10.,10,0.,1.);
  tuple->PostPreS_MassErrVsIas->Sumw2();

  Name = "PostPreS_EtaVsGenID";
  tuple->PostPreS_EtaVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(),  50, -2.6, 2.6, 4000, 0.0, 4000.0);
  tuple->PostPreS_EtaVsGenID->Sumw2();
  Name = "PostPreS_ProbQVsGenID";
  tuple->PostPreS_ProbQVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 100, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbQVsGenID->Sumw2();
  tuple->PostPreS_ProbQVsGenEnviromentID = dir.make<TH2F>("PostPreS_ProbQVsGenEnviromentID",";ProbQ;GenEnviromentID",20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbQVsGenEnviromentID->Sumw2();
  Name = "PostPreS_ProbXYVsGenID";
  tuple->PostPreS_ProbXYVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 100, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbXYVsGenID->Sumw2();
  Name = "PostPreS_PtVsGenID";
  tuple->PostPreS_PtVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 4000, 0.0, 4000.0);
  tuple->PostPreS_PtVsGenID->Sumw2();
  Name = "PostPreS_EoPVsGenID";
  tuple->PostPreS_EoPVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 25, 0, 1.5, 4000, 0.0, 4000.0);
  tuple->PostPreS_EoPVsGenID->Sumw2();
  Name = "PostPreS_IhVsGenID";  
  tuple->PostPreS_IhVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
  tuple->PostPreS_IhVsGenID->Sumw2();
  Name = "PostPreS_IasVsGenID";
  tuple->PostPreS_IasVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsGenID->Sumw2();
  tuple->PostPreS_IasVsGenEnviromentID = dir.make<TH2F>("PostPreS_IasVsGenEnviromentID",";Ias;GenEnviromentID", 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsGenEnviromentID->Sumw2();
  Name = "PostPreS_massTVsGenID";
  tuple->PostPreS_massTVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_massTVsGenID->Sumw2();
  Name = "PostPreS_miniIsoChgVsGenID";
  tuple->PostPreS_miniIsoChgVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoChgVsGenID->Sumw2();
  Name = "PostPreS_miniIsoAllVsGenID";
  tuple->PostPreS_miniIsoAllVsGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoAllVsGenID->Sumw2();
  tuple->PostPreS_MassVsGenID = dir.make<TH2F>("PostPreS_MassVsGenID",";Mass (GeV);GenID",80,0.,4000.,4000, 0.0, 4000.0);
  tuple->PostPreS_MassVsGenID->Sumw2();

  Name = "PostPreS_EtaVsMomGenID";
  tuple->PostPreS_EtaVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(),  50, -2.6, 2.6, 4000, 0.0, 4000.0);
  tuple->PostPreS_EtaVsMomGenID->Sumw2();
  Name = "PostPreS_ProbQVsMomGenID";
  tuple->PostPreS_ProbQVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbQVsMomGenID->Sumw2();
  Name = "PostPreS_ProbXYVsMomGenID";
  tuple->PostPreS_ProbXYVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbXYVsMomGenID->Sumw2();
  Name = "PostPreS_PtVsMomGenID";
  tuple->PostPreS_PtVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 4000, 0.0, 4000.0);
  tuple->PostPreS_PtVsMomGenID->Sumw2();
  Name = "PostPreS_EoPVsMomGenID";
  tuple->PostPreS_EoPVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 25, 0, 1.5, 4000, 0.0, 4000.0);
  tuple->PostPreS_EoPVsMomGenID->Sumw2();
  Name = "PostPreS_IhVsMomGenID";
  tuple->PostPreS_IhVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
  tuple->PostPreS_IhVsMomGenID->Sumw2();
  Name = "PostPreS_IasVsMomGenID";
  tuple->PostPreS_IasVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsMomGenID ->Sumw2();
  Name = "PostPreS_massTVsMomGenID";
  tuple->PostPreS_massTVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_massTVsMomGenID->Sumw2();
  Name = "PostPreS_miniIsoChgVsMomGenID";
  tuple->PostPreS_miniIsoChgVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoChgVsMomGenID->Sumw2();
  Name = "PostPreS_miniIsoAllVsMomGenID";
  tuple->PostPreS_miniIsoAllVsMomGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoAllVsMomGenID->Sumw2();
  tuple->PostPreS_MassVsMomGenID = dir.make<TH2F>("PostPreS_MassVsMomGenID","PostPreS_MassVsMomGenID;Mass;MomGenID",80,0.,4000.,4000, 0.0, 4000.0);
  tuple->PostPreS_MassVsMomGenID->Sumw2();
  
  Name = "PostPreS_EtaVsSiblingGenID";
  tuple->PostPreS_EtaVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(),  50, -2.6, 2.6, 4000, 0.0, 4000.0);
  tuple->PostPreS_EtaVsSiblingGenID->Sumw2();
  Name = "PostPreS_ProbQVsSiblingGenID";
  tuple->PostPreS_ProbQVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbQVsSiblingGenID->Sumw2();
  Name = "PostPreS_ProbXYVsSiblingGenID";
  tuple->PostPreS_ProbXYVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbXYVsSiblingGenID->Sumw2();
  Name = "PostPreS_PtVsSiblingGenID";
  tuple->PostPreS_PtVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 4000, 0.0, 4000.0);
  tuple->PostPreS_PtVsSiblingGenID->Sumw2();
  Name = "PostPreS_EoPVsSiblingGenID";
  tuple->PostPreS_EoPVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 25, 0, 1.5, 4000, 0.0, 4000.0);
  tuple->PostPreS_EoPVsSiblingGenID->Sumw2();
  Name = "PostPreS_IhVsSiblingGenID";
  tuple->PostPreS_IhVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
  tuple->PostPreS_IhVsSiblingGenID->Sumw2();
  Name = "PostPreS_IasVsSiblingGenID";
  tuple->PostPreS_IasVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsSiblingGenID ->Sumw2();
  Name = "PostPreS_massTVsSiblingGenID";
  tuple->PostPreS_massTVsSiblingGenID = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_massTVsSiblingGenID->Sumw2();
  tuple->PostPreS_MassVsSiblingGenID = dir.make<TH2F>("PostPreS_MassVsSiblingGenID","PostPreS_MassVsSiblingGenID;Mass;SiblingGenID",80, 0.0, 4000.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_MassVsSiblingGenID->Sumw2();
  
  Name = "PostPreS_EtaVsGenAngle";
  tuple->PostPreS_EtaVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(),  50, -2.6, 2.6, 100, 0.0, 1.0);
  tuple->PostPreS_EtaVsGenAngle->Sumw2();
  Name = "PostPreS_ProbQVsGenAngle";
  tuple->PostPreS_ProbQVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0,1.0);
  tuple->PostPreS_ProbQVsGenAngle->Sumw2();
  Name = "PostPreS_ProbXYVsGenAngle";
  tuple->PostPreS_ProbXYVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_ProbXYVsGenAngle->Sumw2();
  Name = "PostPreS_PtVsGenAngle";
  tuple->PostPreS_PtVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0.0, 1.0);
  tuple->PostPreS_PtVsGenAngle->Sumw2();
  Name = "PostPreS_EoPVsGenAngle";
  tuple->PostPreS_EoPVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 25, 0, 1.5, 100, 0.0, 1.0);
  tuple->PostPreS_EoPVsGenAngle->Sumw2();
  Name = "PostPreS_IhVsGenAngle";
  tuple->PostPreS_IhVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim, 100, 0.0, 1.0);
  tuple->PostPreS_IhVsGenAngle->Sumw2();
  Name = "PostPreS_IasVsGenAngle";
  tuple->PostPreS_IasVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 10, 0., 1., 100, 0.0, 1.0);
  tuple->PostPreS_IasVsGenAngle ->Sumw2();
  Name = "PostPreS_massTVsGenAngle";
  tuple->PostPreS_massTVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0, 100, 0.0, 1.0);
  tuple->PostPreS_massTVsGenAngle->Sumw2();
  Name = "PostPreS_miniIsoChgVsGenAngle";
  tuple->PostPreS_miniIsoChgVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoChgVsGenAngle->Sumw2();
  Name = "PostPreS_miniIsoAllVsGenAngle";
  tuple->PostPreS_miniIsoAllVsGenAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoAllVsGenAngle->Sumw2();
  tuple->PostPreS_MassVsGenAngle = dir.make<TH2F>("PostPreS_MassVsGenAngle",";Mass;GenAngle",80, 0.0,4000.0,  100, 0.0, 1.0);
  tuple->PostPreS_MassVsGenAngle->Sumw2();
  
  Name = "PostPreS_EtaVsGenMomAngle";
  tuple->PostPreS_EtaVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(),  50, -2.6, 2.6, 100, 0.0, 1.0);
  tuple->PostPreS_EtaVsGenMomAngle->Sumw2();
  Name = "PostPreS_ProbQVsGenMomAngle";
  tuple->PostPreS_ProbQVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0,1.0);
  tuple->PostPreS_ProbQVsGenMomAngle->Sumw2();
  Name = "PostPreS_ProbXYVsGenMomAngle";
  tuple->PostPreS_ProbXYVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_ProbXYVsGenMomAngle->Sumw2();
  Name = "PostPreS_PtVsGenMomAngle";
  tuple->PostPreS_PtVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0.0, 1.0);
  tuple->PostPreS_PtVsGenMomAngle->Sumw2();
  Name = "PostPreS_EoPVsGenMomAngle";
  tuple->PostPreS_EoPVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 25, 0, 1.5, 100, 0.0, 1.0);
  tuple->PostPreS_EoPVsGenMomAngle->Sumw2();
  Name = "PostPreS_IhVsGenMomAngle";
  tuple->PostPreS_IhVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 100, 0, dEdxM_UpLim, 100, 0.0, 1.0);
  tuple->PostPreS_IhVsGenMomAngle->Sumw2();
  Name = "PostPreS_IasVsGenMomAngle";
  tuple->PostPreS_IasVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 10, 0., 1., 100, 0.0, 1.0);
  tuple->PostPreS_IasVsGenMomAngle ->Sumw2();
  Name = "PostPreS_massTVsGenMomAngle";
  tuple->PostPreS_massTVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0, 100, 0.0, 1.0);
  tuple->PostPreS_massTVsGenMomAngle->Sumw2();
  Name = "PostPreS_miniIsoChgVsGenMomAngle";
  tuple->PostPreS_miniIsoChgVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoChgVsGenMomAngle->Sumw2();
  Name = "PostPreS_miniIsoAllVsGenMomAngle";
  tuple->PostPreS_miniIsoAllVsGenMomAngle = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoAllVsGenMomAngle->Sumw2();
  tuple->PostPreS_MassVsGenMomAngle = dir.make<TH2F>("PostPreS_MassVsGenMomAngle",";Mass;GenMomAngle",80,0.,4000., 100, 0.0, 1.0);
  tuple->PostPreS_MassVsGenMomAngle->Sumw2();

  Name = "PostPreS_ProbQVsIas";
  tuple->PostPreS_ProbQVsIas = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0., 1., 20, 0., 1.);
  tuple->PostPreS_ProbQVsIas->Sumw2();

  Name = "PostPreS_EtaVsGenNumSibling";
  tuple->PostPreS_EtaVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(),  50, -2.6, 2.6, 100, 0.0, 10.0);
  tuple->PostPreS_EtaVsGenNumSibling->Sumw2();
  Name = "PostPreS_ProbQVsGenNumSibling";
  tuple->PostPreS_ProbQVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0., 1., 10, 0.0,10.);
  tuple->PostPreS_ProbQVsGenNumSibling->Sumw2();
  Name = "PostPreS_ProbXYVsGenNumSibling";
  tuple->PostPreS_ProbXYVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 10, 0.0, 10.0);
  tuple->PostPreS_ProbXYVsGenNumSibling->Sumw2();
  Name = "PostPreS_PtVsGenNumSibling";
  tuple->PostPreS_PtVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0.0, 10.0);
  tuple->PostPreS_PtVsGenNumSibling->Sumw2();
  Name = "PostPreS_EoPVsGenNumSibling";
  tuple->PostPreS_EoPVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 25, 0, 1.5, 100, 0.0, 10.0);
  tuple->PostPreS_EoPVsGenNumSibling->Sumw2();
  Name = "PostPreS_IhVsGenNumSibling";
  tuple->PostPreS_IhVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 200, 0, dEdxM_UpLim, 100, 0.0, 10.0);
  tuple->PostPreS_IhVsGenNumSibling->Sumw2();
  Name = "PostPreS_IasVsGenNumSibling";
  tuple->PostPreS_IasVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 10, 0., 1., 100, 0.0, 10.0);
  tuple->PostPreS_IasVsGenNumSibling ->Sumw2();
  Name = "PostPreS_massTVsGenNumSibling";
  tuple->PostPreS_massTVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0.0, 250.0, 100, 0.0, 10.0);
  tuple->PostPreS_massTVsGenNumSibling->Sumw2();
  Name = "PostPreS_miniIsoChgVsGenNumSibling";
  tuple->PostPreS_miniIsoChgVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 10.0);
  tuple->PostPreS_miniIsoChgVsGenNumSibling->Sumw2();
  Name = "PostPreS_miniIsoAllVsGenNumSibling";
  tuple->PostPreS_miniIsoAllVsGenNumSibling = dir.make<TH2F>(Name.c_str(), Name.c_str(), 20, 0.0, 1.0, 100, 0.0, 10.0);
  tuple->PostPreS_miniIsoAllVsGenNumSibling->Sumw2();

  Name = "PostPreS_EoPVsPfType";
  tuple->PostPreS_EoPVsPfType = dir.make<TH2F>(Name.c_str(), Name.c_str(), 25, 0.0, 1.5, 9, 0.0, 9.0);
  tuple->PostPreS_EoPVsPfType->Sumw2();
  tuple->PostPreS_Mass = dir.make<TH1F>("PostPreS_Mass","PostPreS_Mass;Mass (GeV);Tracks / 50 GeV", 80,0.,4000.);
  tuple->PostPreS_Mass->Sumw2();
  tuple->PostPreS_MassVsPfType = dir.make<TH2F>("PostPreS_MassVsPfType","PostPreS_MassVsPfType;Mass (GeV);PF ID",80,0.,4000.,9, 0.0, 9.0);
  tuple->PostPreS_MassVsPfType->Sumw2();
  tuple->PostPreS_MassVsPt = dir.make<TH2F>("PostPreS_MassVsPt", ";Mass (GeV);", 80,0.,4000.,80,0.,4000.);
  tuple->PostPreS_MassVsPt->Sumw2();
  tuple->PostPreS_MassVsP = dir.make<TH2F>("PostPreS_MassVsP", ";Mass (GeV);", 80,0.,4000.,80,0.,4000.);
  tuple->PostPreS_MassVsP->Sumw2();
  tuple->PostPreS_MassVsTNOHFraction = dir.make<TH2F>("PostPreS_MassVsTNOHFraction", ";Mass (GeV);", 80,0.,4000.,50, 0, 1);
  tuple->PostPreS_MassVsTNOHFraction->Sumw2();
  tuple->PostPreS_MassVsTNOPH = dir.make<TH2F>("PostPreS_MassVsTNOPH", ";Mass (GeV);", 80,0.,4000.,8, 0, 8);
  tuple->PostPreS_MassVsTNOPH->Sumw2();
  tuple->PostPreS_MassVsTNOM = dir.make<TH2F>("PostPreS_MassVsTNOM", ";Mass (GeV);", 80,0.,4000.,40, 0, 40);
  tuple->PostPreS_MassVsTNOM->Sumw2();
  tuple->PostPreS_MassVsProbQNoL1 = dir.make<TH2F>("PostPreS_MassVsProbQNoL1", ";Mass (GeV);", 80,0.,4000.,20,0.,1.);
  tuple->PostPreS_MassVsProbQNoL1->Sumw2();
  tuple->PostPreS_MassVsProbXYNoL1 = dir.make<TH2F>("PostPreS_MassVsProbXYNoL1", ";Mass (GeV);", 80,0.,4000.,20,0.,1.);
  tuple->PostPreS_MassVsProbXYNoL1->Sumw2();
  tuple->PostPreS_MassVsEoP = dir.make<TH2F>("PostPreS_MassVsEoP", ";Mass (GeV);", 80,0.,4000.,30, 0, 0.3);
  tuple->PostPreS_MassVsEoP->Sumw2();
  tuple->PostPreS_MassVsSumpTOverpT = dir.make<TH2F>("PostPreS_MassVsSumpTOverpT", ";Mass (GeV);", 80,0.,4000.,80, 0, 2);
  tuple->PostPreS_MassVsSumpTOverpT->Sumw2();
  tuple->PostPreS_MassVsPtErrOverPt = dir.make<TH2F>("PostPreS_MassVsPtErrOverPt", ";Mass (GeV);", 80,0.,4000.,40, 0, 1);
  tuple->PostPreS_MassVsPtErrOverPt->Sumw2();
  tuple->PostPreS_MassVsTIsol = dir.make<TH2F>("PostPreS_MassVsTIsol", ";Mass (GeV);", 80,0.,4000., 25, 0, 100);
  tuple->PostPreS_MassVsTIsol->Sumw2();
  tuple->PostPreS_MassVsIh = dir.make<TH2F>("PostPreS_MassVsIh", ";Mass (GeV);", 80,0.,4000.,200, 0, dEdxM_UpLim);
  tuple->PostPreS_MassVsIh->Sumw2();
  tuple->PostPreS_MassVsMassT = dir.make<TH2F>("PostPreS_MassVsMassT", ";Mass (GeV);", 80,0.,4000.,50, 0.0, 250.0);
  tuple->PostPreS_MassVsMassT->Sumw2();
  tuple->PostPreS_MassVsMiniRelIsoAll = dir.make<TH2F>("PostPreS_MassVsMiniRelIsoAll", ";Mass (GeV);", 80,0.,4000.,20, 0., 0.2);
  tuple->PostPreS_MassVsMiniRelIsoAll->Sumw2();
  tuple->PostPreS_MassVsMassErr = dir.make<TH2F>("PostPreS_MassVsMassErr", ";Mass (GeV);", 80,0.,4000.,50, 0., 10.);
  tuple->PostPreS_MassVsMassErr->Sumw2();
  
  Name = "PostPreS_IasAllIhVsLayer";
  tuple->PostPreS_IasAllIhVsLayer = dir.make<TH3F>(Name.c_str(), Name.c_str(), 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 35, 0.,35.);
  tuple->PostPreS_IasAllIhVsLayer->Sumw2();
  Name = "PostPreS_IasPixelIhVsLayer";
  tuple->PostPreS_IasPixelIhVsLayer = dir.make<TH3F>(Name.c_str(), Name.c_str(), 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 10, 0.,10.);
  tuple->PostPreS_IasPixelIhVsLayer->Sumw2();
  Name = "PostPreS_IasStripIhVsLayer";
  tuple->PostPreS_IasStripIhVsLayer = dir.make<TH3F>(Name.c_str(), Name.c_str(), 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 25, 0.,25.);
  tuple->PostPreS_IasStripIhVsLayer->Sumw2();

  tuple->PostPreS_CluProbQVsPixelLayer = dir.make<TH2F>("PostPreS_CluProbQVsPixelLayer",";CluProbQ;Layer",20,0.,1.,4,0.,4.);
  tuple->PostPreS_CluProbQVsPixelLayer->Sumw2();
  tuple->PostPreS_CluProbXYVsPixelLayer = dir.make<TH2F>("PostPreS_CluProbXYVsPixelLayer",";CluProbXY;Layer",100,0.,1.,4,0.,4.);
  tuple->PostPreS_CluProbXYVsPixelLayer->Sumw2();
  tuple->PostPreS_CluSizeVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeVsPixelLayer",";CluSize;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeVsPixelLayer->Sumw2();
  tuple->PostPreS_CluSizeXVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeXVsPixelLayer",";CluSizeX;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeXVsPixelLayer->Sumw2();
  tuple->PostPreS_CluSizeYVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeYVsPixelLayer",";CluSizeY;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeYVsPixelLayer->Sumw2();
  tuple->PostPreS_CluSpecInCPEVsPixelLayer = dir.make<TH2F>("PostPreS_CluSpecInCPEVsPixelLayer",";CluSpecInCPE;Layer",4,0.,4.,4,0.,4.);
  tuple->PostPreS_CluSpecInCPEVsPixelLayer->Sumw2();

  tuple->PostPreS_CluProbQVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluProbQVsPixelLayer_highIas",";CluProbQ;Layer",20,0.,1.,4,0.,4.);
  tuple->PostPreS_CluProbQVsPixelLayer_highIas->Sumw2();
  tuple->PostPreS_CluProbXYVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluProbXYVsPixelLayer_highIas",";CluProbXY;Layer",100,0.,1.,4,0.,4.);
  tuple->PostPreS_CluProbXYVsPixelLayer_highIas->Sumw2();
  tuple->PostPreS_CluSizeVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeVsPixelLayer_highIas",";CluSize;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeVsPixelLayer_highIas->Sumw2();
  tuple->PostPreS_CluSizeXVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeXVsPixelLayer_highIas",";CluSizeX;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeXVsPixelLayer_highIas->Sumw2();
  tuple->PostPreS_CluSizeYVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeYVsPixelLayer_highIas",";CluSizeY;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeYVsPixelLayer_highIas->Sumw2();
  tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSpecInCPEVsPixelLayer_highIas",";CluSpecInCPE;Layer",4,0.,4.,4,0.,4.);
  tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->Sumw2();

  tuple->PostPreS_CluCotBetaVsPixelLayer_lowProbXY = dir.make<TH2F>("PostPreS_CluCotBetaVsPixelLayer_lowProbXY",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->PostPreS_CluCotBetaVsPixelLayer_lowProbXY->Sumw2();
  tuple->PostPreS_CluCotAlphaVsPixelLayer_lowProbXY = dir.make<TH2F>("PostPreS_CluCotAlphaVsPixelLayer_lowProbXY",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
  tuple->PostPreS_CluCotAlphaVsPixelLayer_lowProbXY->Sumw2();
  tuple->PostPreS_CluCotBetaVsPixelLayer = dir.make<TH2F>("PostPreS_CluCotBetaVsPixelLayer",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->PostPreS_CluCotBetaVsPixelLayer->Sumw2();
  tuple->PostPreS_CluCotAlphaVsPixelLayer = dir.make<TH2F>("PostPreS_CluCotAlphaVsPixelLayer",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
  tuple->PostPreS_CluCotAlphaVsPixelLayer->Sumw2();

  tuple->PostPreS_dRMinPfJet = dir.make<TH1F>("PostPreS_dRMinPfJet",";dRMinPfJet",100,0.,1.5);
  tuple->PostPreS_dRMinPfJet->Sumw2();
  tuple->PostPreS_closestPfJetMuonFraction = dir.make<TH1F>("PostPreS_closestPfJetMuonFraction",":closestPfJetMuonFraction",20,0.,1.);
  tuple->PostPreS_closestPfJetMuonFraction->Sumw2();
  tuple->PostPreS_closestPfJetElectronFraction = dir.make<TH1F>("PostPreS_closestPfJetElectronFraction",";closestPfJetElectronFraction",20,0.,1.);
  tuple->PostPreS_closestPfJetElectronFraction->Sumw2();
  tuple->PostPreS_closestPfJetPhotonFraction = dir.make<TH1F>("PostPreS_closestPfJetPhotonFraction",";closestPfJetPhotonFraction",20,0.,1.);
  tuple->PostPreS_closestPfJetPhotonFraction->Sumw2();

  tuple->PostPreS_closestPfJetMuonFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetMuonFractionVsIas",":closestPfJetMuonFraction;Ias",20,0.,1.,20,0.,1.);
  tuple->PostPreS_closestPfJetMuonFractionVsIas->Sumw2();
  tuple->PostPreS_closestPfJetElectronFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetElectronFractionVsIas",";closestPfJetElectronFraction;Ias",20,0.,1.,20,0.,1.);
  tuple->PostPreS_closestPfJetElectronFractionVsIas->Sumw2();
  tuple->PostPreS_closestPfJetPhotonFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetPhotonFractionVsIas",";closestPfJetPhotonFraction;Ias",20,0.,1.,20,0.,1.);
  tuple->PostPreS_closestPfJetPhotonFractionVsIas->Sumw2();

  tuple->PostPreS_dRMinPfJetVsIas = dir.make<TH2F>("PostPreS_dRMinPfJetVsIas", ";dRMinPfJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->PostPreS_dRMinPfJetVsIas->Sumw2();
  tuple->PostPreS_dRMinCaloJet = dir.make<TH1F>("PostPreS_dRMinCaloJet",";dRMinCaloJet",100,0.,1.5);
  tuple->PostPreS_dRMinCaloJet->Sumw2();
  tuple->PostPreS_dRMinCaloJetVsIas =  dir.make<TH2F>("PostPreS_dRMinCaloJetVsIas",";dRMinCaloJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->PostPreS_dRMinCaloJetVsIas->Sumw2();
  
  tuple->GenLevelBinning = dir.make<TH1F>("GenLevelBinning","GenLevelBinning",1200,0.,1200.);
  tuple->GenLevelBinning->Sumw2();
  tuple->GenLevelpT = dir.make<TH1F>("GenLevelpT", "GenLevelpT;Generator p_{T} (GeV)", 50, 0, PtHistoUpperBound);
  tuple->GenLevelpT->Sumw2();
  tuple->GenLevelEta = dir.make<TH1F>("GenLevelEta",";Generator #eta", 60, -3, 3);
  tuple->GenLevelEta->Sumw2();
  tuple->GenLevelBeta = dir.make<TH1F>("GenLevelBeta",";Generator #beta", 20, 0, 1);
  tuple->GenLevelBeta->Sumw2();
  tuple->GenLevelBetaGamma = dir.make<TH1F>("GenLevelBetaGamma",";Generator #beta #gamma;Gen canidate",4500,0.,450.);
  tuple->GenLevelBetaGamma->Sumw2();

  //Initialize histograms for number of bins.  For everything but muon only PredBins=0 so no histograms created
  for (int i = 0; i < PredBins; i++) {
    char Suffix[1024];
    sprintf(Suffix, "_%i", i);
    Name = "BefPreS_Pt_Binned";
    Name.append(Suffix);
    tuple->BefPreS_Pt_Binned[std::to_string(i)] = dir.make<TH1F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_Binned[std::to_string(i)]->Sumw2();
    Name = "BefPreS_TOF_Binned";
    Name.append(Suffix);
    tuple->BefPreS_TOF_Binned[std::to_string(i)] = dir.make<TH1F>(Name.c_str(), Name.c_str(), 150, -1, 5);
    tuple->BefPreS_TOF_Binned[std::to_string(i)]->Sumw2();
  }

  Name = "AS_Eta_RegionA";
  tuple->AS_Eta_RegionA = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionA->Sumw2();
  Name = "AS_Eta_RegionB";
  tuple->AS_Eta_RegionB = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionB->Sumw2();
  Name = "AS_Eta_RegionC";
  tuple->AS_Eta_RegionC = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionC->Sumw2();
  Name = "AS_Eta_RegionD";
  tuple->AS_Eta_RegionD = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionD->Sumw2();
  Name = "AS_Eta_RegionE";
  tuple->AS_Eta_RegionE = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionE->Sumw2();
  Name = "AS_Eta_RegionF";
  tuple->AS_Eta_RegionF = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionF->Sumw2();
  Name = "AS_Eta_RegionG";
  tuple->AS_Eta_RegionG = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionG->Sumw2();
  Name = "AS_Eta_RegionH";
  tuple->AS_Eta_RegionH = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->AS_Eta_RegionH->Sumw2();

  Name = "AS_P";
  tuple->AS_P = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, PtHistoUpperBound);
  tuple->AS_P->Sumw2();
  Name = "AS_Pt";
  tuple->AS_Pt = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, PtHistoUpperBound);
  tuple->AS_Pt->Sumw2();
  Name = "AS_Ias";
  tuple->AS_Ias = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 10, 0., 1.);
  tuple->AS_Ias->Sumw2();
  Name = "AS_Ih";
  tuple->AS_Ih = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 100, 0, dEdxM_UpLim);
  tuple->AS_Ih->Sumw2();
  Name = "AS_TOF";
  tuple->AS_TOF = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 1, 5);
  tuple->AS_TOF->Sumw2();

  Name = "BefPreS_EtaVsIas";
  tuple->BefPreS_EtaVsIas = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, -3, 3, 10, 0., 1.);
  Name = "BefPreS_EtaVsIh";
  tuple->BefPreS_EtaVsIh = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, -3, 3, 100, 0, dEdxM_UpLim);
  Name = "BefPreS_EtaVsP";
  tuple->BefPreS_EtaVsP = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, -3, 3, 50, 0, PtHistoUpperBound);
  Name = "BefPreS_EtaVsPt";
  tuple->BefPreS_EtaVsPt = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, -3, 3, 50, 0, PtHistoUpperBound);
  Name = "BefPreS_EtaVsTOF";
  tuple->BefPreS_EtaVsTOF = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, -3, 3, 50, 0, 3);
  Name = "BefPreS_EtaVsNBH";
  tuple->BefPreS_EtaVsNBH = dir.make<TH2F>(Name.c_str(), Name.c_str(), 60, -3, 3, 24, 0, 24);
  Name = "BefPreS_EtaVsDz";
  tuple->BefPreS_EtaVsDz = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, -3, 3, 50, -IPbound, IPbound);
  Name = "BefPreS_PVsIas";
  tuple->BefPreS_PVsIas = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0, dEdxS_UpLim);
  Name = "BefPreS_IhVsIas";
  tuple->BefPreS_IhVsIas = dir.make<TH2F>(Name.c_str(), Name.c_str(), 100, 0, dEdxM_UpLim, 100, 0, dEdxS_UpLim);
  Name = "BefPreS_PVsIh";
  tuple->BefPreS_PVsIh = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  Name = "BefPreS_PtVsIas";
  tuple->BefPreS_PtVsIas = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 10, 0., 1.);
  Name = "BefPreS_PtVsIh";
  tuple->BefPreS_PtVsIh = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  Name = "BefPreS_PtTOF";
  tuple->BefPreS_PtTOF = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, PtHistoUpperBound, 50, 0, 5);
  //   Name = "BefPreS_TOFIs"; tuple->BefPreS_TOFIs = dir.make<TH2F>(Name.c_str(), Name.c_str(),                   100, 1, 5, 100, 0, dEdxS_UpLim);
  Name = "BefPreS_TOFIs";
  tuple->BefPreS_TOFIs = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, 5, 10, 0., 1.);
  //   Name = "BefPreS_TOFIh"; BefPreS_TOFIm = dir.make<TH2F>(Name.c_str(), Name.c_str(),                   100, 1, 5, 200, 0, dEdxM_UpLim);
  Name = "BefPreS_TOFIh";
  tuple->BefPreS_TOFIh = dir.make<TH2F>(Name.c_str(), Name.c_str(), 50, 0, 5, 100, 0, dEdxM_UpLim);

  //   Name = "AS_EtaIs"; AS_EtaIs = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 10, 0., 1.);
  //   Name = "AS_EtaIh"; AS_EtaIm = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3,100, 0, dEdxM_UpLim);
  //   Name = "AS_EtaP" ; AS_EtaP  = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, PtHistoUpperBound);
  //   Name = "AS_EtaPt"; AS_EtaPt = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, PtHistoUpperBound);
  //   Name = "AS_EtaTOF"; AS_EtaTOF = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0,  NCuts, 50,-3, 3, 50, 0, 3);
  Name = "AS_PIs";
  tuple->AS_PIs =
      dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 10, 0., 1.);
  Name = "AS_PIh";
  tuple->AS_PIh =
      dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  Name = "AS_PtIs";
  tuple->AS_PtIs =
      dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 10, 0., 1.);
  Name = "AS_PtIh";
  tuple->AS_PtIh =
      dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  Name = "AS_TOFIs";
  tuple->AS_TOFIs = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, 5, 10, 0., 1.);
  Name = "AS_TOFIh";
  tuple->AS_TOFIh = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 50, 0, 5, 100, 0, dEdxM_UpLim);

  Name = "H_D_DzSidebands";
  tuple->H_D_DzSidebands = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, DzRegions, 0, DzRegions);
  tuple->H_D_DzSidebands->Sumw2();

  // Background prediction histograms don't need to be made for signal or individual MC samples
  // if (!isSignal) {
  // Although not needed let's still do it for every input
  if (true) {
    Name = "H_A";
    tuple->H_A = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_A->Sumw2();
    Name = "H_B";
    tuple->H_B = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_B->Sumw2();
    Name = "H_C";
    tuple->H_C = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_C->Sumw2();
    Name = "H_D";
    tuple->H_D = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_D->Sumw2();
    Name = "H_E";
    tuple->H_E = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_E->Sumw2();
    Name = "H_F";
    tuple->H_F = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_F->Sumw2();
    Name = "H_G";
    tuple->H_G = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_G->Sumw2();
    Name = "H_H";
    tuple->H_H = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
    tuple->H_H->Sumw2();

    //Initialize histograms for number of bins.  For everything but muon only PredBins=0 so no histograms created
    for (int i = 0; i < PredBins; i++) {
      char Suffix[1024];
      sprintf(Suffix, "_%i", i);
      Name = "H_B_Binned";
      Name.append(Suffix);
      tuple->H_B_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_B_Binned[std::to_string(i)]->Sumw2();
      Name = "H_D_Binned";
      Name.append(Suffix);
      tuple->H_D_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_D_Binned[std::to_string(i)]->Sumw2();
      Name = "H_F_Binned";
      Name.append(Suffix);
      tuple->H_F_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_F_Binned[std::to_string(i)]->Sumw2();
      Name = "H_H_Binned";
      Name.append(Suffix);
      tuple->H_H_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_H_Binned[std::to_string(i)]->Sumw2();
    }

    Name = "Hist_Is";
    tuple->Hist_Is = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, dEdxS_UpLim);
    tuple->Hist_Is->Sumw2();
    Name = "Hist_Pt";
    tuple->Hist_Pt = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, PtHistoUpperBound);
    tuple->Hist_Pt->Sumw2();
    Name = "Hist_TOF";
    tuple->Hist_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -10, 20);
    tuple->Hist_TOF->Sumw2();
    //The following are only used to create the predicted mass spectrum.  Memory intensive so don't initialize for analyses not doing mass fits
    if (TypeMode < 3) {
      Name = "Pred_I";
      tuple->Pred_I = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 400, 0, dEdxM_UpLim);
      tuple->Pred_I->Sumw2();
      Name = "Pred_EtaI";
      tuple->Pred_EtaI = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 400, 0, dEdxM_UpLim);
      tuple->Pred_EtaI->Sumw2();
      Name = "Pred_EtaB";
      tuple->Pred_EtaB = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3);
      tuple->Pred_EtaB->Sumw2();
      Name = "Pred_EtaS";
      tuple->Pred_EtaS = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3);
      tuple->Pred_EtaS->Sumw2();
      Name = "Pred_EtaS2";
      tuple->Pred_EtaS2 = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3);
      tuple->Pred_EtaS2->Sumw2();
      Name = "Pred_EtaP";
      tuple->Pred_EtaP = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 200, GlobalMinPt, PtHistoUpperBound);
      tuple->Pred_EtaP->Sumw2();
      Name = "Pred_TOF";
      tuple->Pred_TOF = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 200, GlobalMinTOF, 5);
      tuple->Pred_TOF->Sumw2();
      //pz

      Name = "PDF_G_EtaP";
      tuple->PDF_G_EtaP = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 200, GlobalMinPt, PtHistoUpperBound);
      tuple->PDF_G_EtaP->Sumw2();
      Name = "PDF_C_EtaP";
      tuple->PDF_C_EtaP = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 200, GlobalMinPt, PtHistoUpperBound);
      tuple->PDF_C_EtaP->Sumw2();

      Name = "PDF_A_Eta";
      tuple->PDF_A_Eta = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3);
      tuple->PDF_A_Eta->Sumw2();
      Name = "PDF_E_Eta";
      tuple->PDF_E_Eta = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3);
      tuple->PDF_E_Eta->Sumw2();

      Name = "PDF_B_EtaICK";
      tuple->PDF_B_EtaICK = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 60, -2., 3.);
      tuple->PDF_B_EtaICK->Sumw2();
      Name = "PDF_F_EtaICK";
      tuple->PDF_F_EtaICK = dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 60, -2., 3.);
      tuple->PDF_F_EtaICK->Sumw2();

      Name = "PDF_H_EtaMass";
      tuple->PDF_H_EtaMass = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, MassNBins, 0, MassHistoUpperBound);
      tuple->PDF_H_EtaMass->Sumw2();

      //pz FLIP
      Name = "PDF_G_EtaP_Flip";
      tuple->PDF_G_EtaP_Flip = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 200, GlobalMinPt, PtHistoUpperBound);
      tuple->PDF_G_EtaP_Flip->Sumw2();
      Name = "PDF_C_EtaP_Flip";
      tuple->PDF_C_EtaP_Flip = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 200, GlobalMinPt, PtHistoUpperBound);
      tuple->PDF_C_EtaP_Flip->Sumw2();

      Name = "PDF_A_Eta_Flip";
      tuple->PDF_A_Eta_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3);
      tuple->PDF_A_Eta_Flip->Sumw2();
      Name = "PDF_E_Eta_Flip";
      tuple->PDF_E_Eta_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3);
      tuple->PDF_E_Eta_Flip->Sumw2();

      Name = "PDF_B_EtaICK_Flip";
      tuple->PDF_B_EtaICK_Flip =
          dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 60, -2., 3.);
      tuple->PDF_B_EtaICK_Flip->Sumw2();
      Name = "PDF_F_EtaICK_Flip";
      tuple->PDF_F_EtaICK_Flip =
          dir.make<TH3F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, 60, -2., 3.);
      tuple->PDF_F_EtaICK_Flip->Sumw2();

      Name = "PDF_H_EtaMass_Flip";
      tuple->PDF_H_EtaMass_Flip = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts, 0, NCuts, EtaBins, -3, 3, MassNBins, 0, MassHistoUpperBound);
      tuple->PDF_H_EtaMass_Flip->Sumw2();
    }

    Name = "RegionD_I";
    tuple->RegionD_I = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 400, 0, dEdxM_UpLim);
    tuple->RegionD_I->Sumw2();
    Name = "RegionD_Ias";
    tuple->RegionD_Ias = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 100, 0, dEdxS_UpLim);
    tuple->RegionD_Ias->Sumw2();
    Name = "RegionD_P";
    tuple->RegionD_P = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 200, GlobalMinPt, PtHistoUpperBound);
    tuple->RegionD_P->Sumw2();
    Name = "RegionD_TOF";
    tuple->RegionD_TOF = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 200, GlobalMinTOF, 5);
    tuple->RegionD_TOF->Sumw2();

    Name = "RegionH_Ias";
    tuple->RegionH_Ias = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts, 100, 0, dEdxS_UpLim);
    tuple->RegionH_Ias->Sumw2();

    Name = "H_A_Flip";
    tuple->H_A_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_A_Flip->Sumw2();
    Name = "H_B_Flip";
    tuple->H_B_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_B_Flip->Sumw2();
    Name = "H_C_Flip";
    tuple->H_C_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_C_Flip->Sumw2();
    Name = "H_D_Flip";
    tuple->H_D_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_D_Flip->Sumw2();
    Name = "H_E_Flip";
    tuple->H_E_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_E_Flip->Sumw2();
    Name = "H_F_Flip";
    tuple->H_F_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_F_Flip->Sumw2();
    Name = "H_G_Flip";
    tuple->H_G_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_G_Flip->Sumw2();
    Name = "H_H_Flip";
    tuple->H_H_Flip = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip);
    tuple->H_H_Flip->Sumw2();

    for (int i = 0; i < PredBins; i++) {
      char Suffix[1024];
      sprintf(Suffix, "_%i", i);
      Name = "H_B_Binned_Flip";
      Name.append(Suffix);
      tuple->H_B_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_B_Binned_Flip[std::to_string(i)]->Sumw2();
      Name = "H_D_Binned_Flip";
      Name.append(Suffix);
      tuple->H_D_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_D_Binned_Flip[std::to_string(i)]->Sumw2();
      Name = "H_F_Binned_Flip";
      Name.append(Suffix);
      tuple->H_F_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_F_Binned_Flip[std::to_string(i)]->Sumw2();
      Name = "H_H_Binned_Flip";
      Name.append(Suffix);
      tuple->H_H_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), NCuts, 0, NCuts);
      tuple->H_H_Binned_Flip[std::to_string(i)]->Sumw2();
    }

    //The following are only used to create the predicted mass spectrum.  Memory intensive so don't initialize for analyses not doing mass fits
    if (TypeMode < 3) {
      Name = "Pred_I_Flip";
      tuple->Pred_I_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, 400, 0, dEdxM_UpLim);
      tuple->Pred_I_Flip->Sumw2();
      Name = "Pred_EtaB_Flip";
      tuple->Pred_EtaB_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, EtaBins, -3, 3);
      tuple->Pred_EtaB_Flip->Sumw2();
      Name = "Pred_EtaS_Flip";
      tuple->Pred_EtaS_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, EtaBins, -3, 3);
      tuple->Pred_EtaS_Flip->Sumw2();
      Name = "Pred_EtaS2_Flip";
      tuple->Pred_EtaS2_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, EtaBins, -3, 3);
      tuple->Pred_EtaS2_Flip->Sumw2();
      Name = "Pred_EtaP_Flip";
      tuple->Pred_EtaP_Flip = dir.make<TH3F>(
          Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, EtaBins, -3, 3, 200, GlobalMinPt, PtHistoUpperBound);
      tuple->Pred_EtaP_Flip->Sumw2();
      Name = "Pred_TOF_Flip";
      tuple->Pred_TOF_Flip =
          dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, 200, GlobalMinTOF, 5);
      tuple->Pred_TOF_Flip->Sumw2();
    }

    Name = "RegionD_I_Flip";
    tuple->RegionD_I_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, 400, 0, dEdxM_UpLim);
    tuple->RegionD_I_Flip->Sumw2();
    Name = "RegionD_Ias_Flip";
    tuple->RegionD_Ias_Flip =
        dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, 100, 0, dEdxS_UpLim);
    tuple->RegionD_Ias_Flip->Sumw2();
    Name = "RegionD_P_Flip";
    tuple->RegionD_P_Flip =
        dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, 200, GlobalMinPt, PtHistoUpperBound);
    tuple->RegionD_P_Flip->Sumw2();
    Name = "RegionD_TOF_Flip";
    tuple->RegionD_TOF_Flip = dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, 200, -3, 1);
    tuple->RegionD_TOF_Flip->Sumw2();

    Name = "RegionH_Ias_Flip";
    tuple->RegionH_Ias_Flip =
        dir.make<TH2F>(Name.c_str(), Name.c_str(), NCuts_Flip, 0, NCuts_Flip, 100, 0, dEdxS_UpLim);
    tuple->RegionH_Ias_Flip->Sumw2();

    Name = "CtrlPt_S1_Is";
    tuple->CtrlPt_S1_Is = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, dEdxS_UpLim);
    tuple->CtrlPt_S1_Is->Sumw2();
    Name = "CtrlPt_S2_Is";
    tuple->CtrlPt_S2_Is = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, dEdxS_UpLim);
    tuple->CtrlPt_S2_Is->Sumw2();
    Name = "CtrlPt_S3_Is";
    tuple->CtrlPt_S3_Is = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, dEdxS_UpLim);
    tuple->CtrlPt_S3_Is->Sumw2();
    Name = "CtrlPt_S4_Is";
    tuple->CtrlPt_S4_Is = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, dEdxS_UpLim);
    tuple->CtrlPt_S4_Is->Sumw2();

    Name = "CtrlPt_S1_Ih";
    tuple->CtrlPt_S1_Ih = dir.make<TH1D>(Name.c_str(), Name.c_str(), 400, 0, dEdxM_UpLim);
    tuple->CtrlPt_S1_Ih->Sumw2();
    Name = "CtrlPt_S2_Ih";
    tuple->CtrlPt_S2_Ih = dir.make<TH1D>(Name.c_str(), Name.c_str(), 400, 0, dEdxM_UpLim);
    tuple->CtrlPt_S2_Ih->Sumw2();
    Name = "CtrlPt_S3_Ih";
    tuple->CtrlPt_S3_Ih = dir.make<TH1D>(Name.c_str(), Name.c_str(), 400, 0, dEdxM_UpLim);
    tuple->CtrlPt_S3_Ih->Sumw2();
    Name = "CtrlPt_S4_Ih";
    tuple->CtrlPt_S4_Ih = dir.make<TH1D>(Name.c_str(), Name.c_str(), 400, 0, dEdxM_UpLim);
    tuple->CtrlPt_S4_Ih->Sumw2();

    Name = "CtrlIs_S1_TOF";
    tuple->CtrlIs_S1_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIs_S1_TOF->Sumw2();
    Name = "CtrlIs_S2_TOF";
    tuple->CtrlIs_S2_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIs_S2_TOF->Sumw2();
    Name = "CtrlIs_S3_TOF";
    tuple->CtrlIs_S3_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIs_S3_TOF->Sumw2();
    Name = "CtrlIs_S4_TOF";
    tuple->CtrlIs_S4_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIs_S4_TOF->Sumw2();

    Name = "CtrlIh_S1_TOF";
    tuple->CtrlIh_S1_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIh_S1_TOF->Sumw2();
    Name = "CtrlIh_S2_TOF";
    tuple->CtrlIh_S2_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIh_S2_TOF->Sumw2();
    Name = "CtrlIh_S3_TOF";
    tuple->CtrlIh_S3_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIh_S3_TOF->Sumw2();
    Name = "CtrlIh_S4_TOF";
    tuple->CtrlIh_S4_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, 0, 5);
    tuple->CtrlIh_S4_TOF->Sumw2();

    Name = "CtrlPt_S1_TOF";
    tuple->CtrlPt_S1_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
    tuple->CtrlPt_S1_TOF->Sumw2();
    Name = "CtrlPt_S2_TOF";
    tuple->CtrlPt_S2_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
    tuple->CtrlPt_S2_TOF->Sumw2();
    Name = "CtrlPt_S3_TOF";
    tuple->CtrlPt_S3_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
    tuple->CtrlPt_S3_TOF->Sumw2();
    Name = "CtrlPt_S4_TOF";
    tuple->CtrlPt_S4_TOF = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
    tuple->CtrlPt_S4_TOF->Sumw2();

    for (int i = 0; i < PredBins; i++) {
      char Suffix[1024];
      sprintf(Suffix, "_%i", i);
      Name = "CtrlPt_S1_TOF_Binned";
      Name.append(Suffix);
      tuple->CtrlPt_S1_TOF_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
      tuple->CtrlPt_S1_TOF_Binned[std::to_string(i)]->Sumw2();
      Name = "CtrlPt_S2_TOF_Binned";
      Name.append(Suffix);
      tuple->CtrlPt_S2_TOF_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
      tuple->CtrlPt_S2_TOF_Binned[std::to_string(i)]->Sumw2();
      Name = "CtrlPt_S3_TOF_Binned";
      Name.append(Suffix);
      tuple->CtrlPt_S3_TOF_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
      tuple->CtrlPt_S3_TOF_Binned[std::to_string(i)]->Sumw2();
      Name = "CtrlPt_S4_TOF_Binned";
      Name.append(Suffix);
      tuple->CtrlPt_S4_TOF_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), Name.c_str(), 200, -2, 7);
      tuple->CtrlPt_S4_TOF_Binned[std::to_string(i)]->Sumw2();
    }
  }

  //===================================================
  //
  //  HSCPCandidates ttree: different saving-levels
  //
  //===================================================
  tuple->Tree = dir.make<TTree>("HscpCandidates", "HscpCandidates");
  ///tuple->Tree->SetDirectory(0);
  if (saveTree > 0) {
    tuple->Tree->Branch("Trig", &tuple->Tree_Trig, "Trig/i");
    tuple->Tree->Branch("Run", &tuple->Tree_Run, "Run/i");
    tuple->Tree->Branch("Event", &tuple->Tree_Event, "Event/l");
    tuple->Tree->Branch("Lumi", &tuple->Tree_Lumi, "Lumi/i");
    tuple->Tree->Branch("PileUp", &tuple->Tree_PileUp, "PileUp/i");
    tuple->Tree->Branch("nofVtx", &tuple->Tree_nofVertices, "nofVtx/i");
    tuple->Tree->Branch("Hscp", &tuple->Tree_Hscp, "Hscp/i");
    tuple->Tree->Branch("nmuons", &tuple->Tree_nmuons, "nmuons/i");
    tuple->Tree->Branch("njets", &tuple->Tree_njets, "njets/i");
    tuple->Tree->Branch("Weight", &tuple->Tree_Weight, "Weight/F");
    tuple->Tree->Branch("GeneratorWeight", &tuple->Tree_GeneratorWeight, "GeneratorWeight/F");
    tuple->Tree->Branch("GeneratorBinningValues", &tuple->Tree_GeneratorBinningValues, "GeneratorBinningValues/F");
    tuple->Tree->Branch("HLT_Mu50", &tuple->Tree_HLT_Mu50, "HLT_Mu50/O");
    tuple->Tree->Branch(
        "HLT_PFMET120_PFMHT120_IDTight", &tuple->Tree_HLT_PFMET120_PFMHT120_IDTight, "HLT_PFMET120_PFMHT120_IDTight/O");
    tuple->Tree->Branch("HLT_PFHT500_PFMET100_PFMHT100_IDTight",
                        &tuple->Tree_HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                        "HLT_PFHT500_PFMET100_PFMHT100_IDTight/O");
    tuple->Tree->Branch("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60",
                        &tuple->Tree_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                        "HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60/O");
    tuple->Tree->Branch("HLT_MET105_IsoTrk50", &tuple->Tree_HLT_MET105_IsoTrk50, "HLT_MET105_IsoTrk50/O");
    tuple->Tree->Branch("CaloMET", &tuple->Tree_CaloMET, "CaloMET/F");
    tuple->Tree->Branch("RecoPFMET", &tuple->Tree_RecoPFMET, "RecoPFMET/F");
    tuple->Tree->Branch("RecoPFMHT", &tuple->Tree_RecoPFMHT, "RecoPFMHT/F");
    tuple->Tree->Branch("HLTPFMET", &tuple->Tree_HLTPFMET, "HLTPFMET/F");
    tuple->Tree->Branch("HLTPFMHT", &tuple->Tree_HLTPFMHT, "HLTPFMHT/F");
    tuple->Tree->Branch("RecoPFMET_eta", &tuple->Tree_RecoPFMET_eta, "RecoPFMET_eta/F");
    tuple->Tree->Branch("RecoPFMET_phi", &tuple->Tree_RecoPFMET_phi, "RecoPFMET_phi/F");
    tuple->Tree->Branch("RecoPFMET_significance", &tuple->Tree_RecoPFMET_significance, "RecoPFMET_significance/F");
    tuple->Tree->Branch("Muon1_Pt", &tuple->Tree_Muon1_Pt, "Muon1_Pt/F");
    tuple->Tree->Branch("Muon1_eta", &tuple->Tree_Muon1_eta, "Muon1_eta/F");
    tuple->Tree->Branch("Muon1_phi", &tuple->Tree_Muon1_phi, "Muon1_phi/F");
    tuple->Tree->Branch("Muon2_Pt", &tuple->Tree_Muon2_Pt, "Muon2_Pt/F");
    tuple->Tree->Branch("Muon2_eta", &tuple->Tree_Muon2_eta, "Muon2_eta/F");
    tuple->Tree->Branch("Muon2_phi", &tuple->Tree_Muon2_phi, "Muon2_phi/F");
    if (saveTree > 3) {
      tuple->Tree->Branch("Jet_pt", &tuple->Tree_jet_pt);
      tuple->Tree->Branch("Jet_eta", &tuple->Tree_jet_eta);
      tuple->Tree->Branch("Jet_phi", &tuple->Tree_jet_phi);
      tuple->Tree->Branch("Jet_mass", &tuple->Tree_jet_mass);
      tuple->Tree->Branch("Jet_energy", &tuple->Tree_jet_energy);
      tuple->Tree->Branch("Jet_pdgId", &tuple->Tree_jet_pdgId);
      tuple->Tree->Branch("Jet_et", &tuple->Tree_jet_et);
      tuple->Tree->Branch("Jet_chargedEmEnergyFraction", &tuple->Tree_jet_chargedEmEnergyFraction);
      tuple->Tree->Branch("Jet_neutralEmEnergyFraction", &tuple->Tree_jet_neutralEmEnergyFraction);
    }
    tuple->Tree->Branch("mT", &tuple->Tree_vect_mT);
    if (saveTree > 1) {
      tuple->Tree->Branch("passCutPt55", &tuple->Tree_passCutPt55);
      tuple->Tree->Branch("passPreselection", &tuple->Tree_passPreselection);
      tuple->Tree->Branch("passSelection", &tuple->Tree_passSelection);
    }
    tuple->Tree->Branch("Charge", &tuple->Tree_Charge);
    tuple->Tree->Branch("Pt", &tuple->Tree_Pt);
    tuple->Tree->Branch("PtErr", &tuple->Tree_PtErr);
    tuple->Tree->Branch("Ias", &tuple->Tree_Ias);
    tuple->Tree->Branch("Ias_noTIBnoTIDno3TEC", &tuple->Tree_Ias_noTIBnoTIDno3TEC);
    tuple->Tree->Branch("Ias_PixelOnly", &tuple->Tree_Ias_PixelOnly);
    tuple->Tree->Branch("Ias_StripOnly", &tuple->Tree_Ias_StripOnly);
    tuple->Tree->Branch("Ias_PixelOnly_noL1", &tuple->Tree_Ias_PixelOnly_noL1);
    tuple->Tree->Branch("Ih", &tuple->Tree_Ih);
    tuple->Tree->Branch("Ick", &tuple->Tree_Ick);
    tuple->Tree->Branch("Fmip", &tuple->Tree_Fmip);
    tuple->Tree->Branch("ProbXY", &tuple->Tree_ProbXY);
    tuple->Tree->Branch("ProbXY_noL1", &tuple->Tree_ProbXY_noL1);
    tuple->Tree->Branch("ProbQ", &tuple->Tree_ProbQ);
    tuple->Tree->Branch("ProbQ_noL1", &tuple->Tree_ProbQ_noL1);
    tuple->Tree->Branch("Ndof", &tuple->Tree_Ndof);
    tuple->Tree->Branch("Chi2", &tuple->Tree_Chi2);
    tuple->Tree->Branch("QualityMask", &tuple->Tree_QualityMask);
    tuple->Tree->Branch("isHighPurity", &tuple->Tree_isHighPurity);
    tuple->Tree->Branch("EoverP", &tuple->Tree_EoverP);
    tuple->Tree->Branch("isMuon", &tuple->Tree_isMuon);
    tuple->Tree->Branch("MuonSelector", &tuple->Tree_Muon_selector);
    tuple->Tree->Branch("isElectron", &tuple->Tree_isElectron);
    tuple->Tree->Branch("isChHadron", &tuple->Tree_isChHadron);
    tuple->Tree->Branch("isNeutHadron", &tuple->Tree_isNeutHadron);
    tuple->Tree->Branch("ECAL_energy", &tuple->Tree_ECAL_energy);
    tuple->Tree->Branch("HCAL_energy", &tuple->Tree_HCAL_energy);
    tuple->Tree->Branch("TOF", &tuple->Tree_TOF);
    tuple->Tree->Branch("TOFErr", &tuple->Tree_TOFErr);
    tuple->Tree->Branch("TOF_ndof", &tuple->Tree_TOF_ndof);
    tuple->Tree->Branch("DTTOF", &tuple->Tree_DTTOF);
    tuple->Tree->Branch("DTTOFErr", &tuple->Tree_DTTOFErr);
    tuple->Tree->Branch("DTTOF_ndof", &tuple->Tree_DTTOF_ndof);
    tuple->Tree->Branch("CSCTOF", &tuple->Tree_CSCTOF);
    tuple->Tree->Branch("CSCTOFErr", &tuple->Tree_CSCTOFErr);
    tuple->Tree->Branch("CSCTOF_ndof", &tuple->Tree_CSCTOF_ndof);
    tuple->Tree->Branch("Mass", &tuple->Tree_Mass);
    tuple->Tree->Branch("MassErr", &tuple->Tree_MassErr);
    tuple->Tree->Branch("dZ", &tuple->Tree_dZ);
    tuple->Tree->Branch("dXY", &tuple->Tree_dXY);
    tuple->Tree->Branch("dR", &tuple->Tree_dR);
    tuple->Tree->Branch("eta", &tuple->Tree_eta);
    tuple->Tree->Branch("phi", &tuple->Tree_phi);
    tuple->Tree->Branch("NOH", &tuple->Tree_NOH);
    tuple->Tree->Branch("NOPH", &tuple->Tree_NOPH);
    tuple->Tree->Branch("FOVH", &tuple->Tree_FOVH);
    tuple->Tree->Branch("NOMH", &tuple->Tree_NOMH);
    tuple->Tree->Branch("FOVHD", &tuple->Tree_FOVHD);
    tuple->Tree->Branch("NOM", &tuple->Tree_NOM);
    tuple->Tree->Branch("iso_TK", &tuple->Tree_iso_TK);
    tuple->Tree->Branch("iso_ECAL", &tuple->Tree_iso_ECAL);
    tuple->Tree->Branch("iso_HCAL", &tuple->Tree_iso_HCAL);
    if (saveTree > 1) {
      tuple->Tree->Branch("PFMiniIso_relative", &tuple->Tree_PFMiniIso_relative);
      tuple->Tree->Branch("PFMiniIso_wMuon_relative", &tuple->Tree_PFMiniIso_wMuon_relative);
      tuple->Tree->Branch("TrackPFIsolationR005_sumChargedHadronPt", &tuple->Tree_track_PFIsolationR005_sumChargedHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR005_sumNeutralHadronPt", &tuple->Tree_track_PFIsolationR005_sumNeutralHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR005_sumPhotonPt", &tuple->Tree_track_PFIsolationR005_sumPhotonPt);
      tuple->Tree->Branch("TrackPFIsolationR005_sumPUPt", &tuple->Tree_track_PFIsolationR005_sumPUPt);
      tuple->Tree->Branch("TrackPFIsolationR01_sumChargedHadronPt", &tuple->Tree_track_PFIsolationR01_sumChargedHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR01_sumNeutralHadronPt", &tuple->Tree_track_PFIsolationR01_sumNeutralHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR01_sumPhotonPt", &tuple->Tree_track_PFIsolationR01_sumPhotonPt);
      tuple->Tree->Branch("TrackPFIsolationR01_sumPUPt", &tuple->Tree_track_PFIsolationR01_sumPUPt);
      tuple->Tree->Branch("TrackPFIsolationR03_sumChargedHadronPt", &tuple->Tree_track_PFIsolationR03_sumChargedHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR03_sumNeutralHadronPt", &tuple->Tree_track_PFIsolationR03_sumNeutralHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR03_sumPhotonPt", &tuple->Tree_track_PFIsolationR03_sumPhotonPt);
      tuple->Tree->Branch("TrackPFIsolationR03_sumPUPt", &tuple->Tree_track_PFIsolationR03_sumPUPt);
      tuple->Tree->Branch("TrackPFIsolationR05_sumChargedHadronPt", &tuple->Tree_track_PFIsolationR05_sumChargedHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR05_sumNeutralHadronPt", &tuple->Tree_track_PFIsolationR05_sumNeutralHadronPt);
      tuple->Tree->Branch("TrackPFIsolationR05_sumPhotonPt", &tuple->Tree_track_PFIsolationR05_sumPhotonPt);
      tuple->Tree->Branch("TrackPFIsolationR05_sumPUPt", &tuple->Tree_track_PFIsolationR05_sumPUPt);
      tuple->Tree->Branch("MuonPFIsolationR03_sumChargedHadronPt", &tuple->Tree_muon_PFIsolationR03_sumChargedHadronPt);
      tuple->Tree->Branch("MuonPFIsolationR03_sumNeutralHadronPt", &tuple->Tree_muon_PFIsolationR03_sumNeutralHadronPt);
      tuple->Tree->Branch("MuonPFIsolationR03_sumPhotonPt", &tuple->Tree_muon_PFIsolationR03_sumPhotonPt);
      tuple->Tree->Branch("MuonPFIsolationR03_sumPUPt", &tuple->Tree_muon_PFIsolationR03_sumPUPt);
      tuple->Tree->Branch("Ih_noL1", &tuple->Tree_Ih_noL1);
      tuple->Tree->Branch("Ih_15drop", &tuple->Tree_Ih_15drop);
      tuple->Tree->Branch("Ih_StripOnly", &tuple->Tree_Ih_StripOnly);
      tuple->Tree->Branch("Ih_StripOnly_15drop", &tuple->Tree_Ih_StripOnly_15drop);
      tuple->Tree->Branch("Ih_PixelOnly_noL1", &tuple->Tree_Ih_PixelOnly_noL1);
      tuple->Tree->Branch("Ih_SaturationCorrectionFromFits", &tuple->Tree_Ih_SaturationCorrectionFromFits);
    }
    if (saveTree > 2) {
      tuple->Tree->Branch("clust_charge", &tuple->Tree_clust_charge);
      tuple->Tree->Branch("clust_pathlength", &tuple->Tree_clust_pathlength);
      tuple->Tree->Branch("clust_nstrip", &tuple->Tree_clust_nstrip);
      tuple->Tree->Branch("clust_sat254", &tuple->Tree_clust_sat254);
      tuple->Tree->Branch("clust_sat255", &tuple->Tree_clust_sat255);
      tuple->Tree->Branch("clust_detid", &tuple->Tree_clust_detid);
      tuple->Tree->Branch("clust_isStrip", &tuple->Tree_clust_isStrip);
      tuple->Tree->Branch("clust_isPixel", &tuple->Tree_clust_isPixel);
    }
    if (saveTree > 4) {
      tuple->Tree->Branch("GenId", &tuple->Tree_GenId);
      tuple->Tree->Branch("GenCharge", &tuple->Tree_GenCharge);
      tuple->Tree->Branch("GenMass", &tuple->Tree_GenMass);
      tuple->Tree->Branch("GenPt", &tuple->Tree_GenPt);
      tuple->Tree->Branch("GenEta", &tuple->Tree_GenEta);
      tuple->Tree->Branch("GenPhi", &tuple->Tree_GenPhi);
    }
  }

  tuple->GenTree = dir.make<TTree>("GenHscpCandidates", "GenHscpCandidates");
  ////tuple->GenTree->SetDirectory(0);
  if (saveGenTree > 0) {
    tuple->GenTree->Branch("Run", &tuple->GenTree_Run, "Run/i");
    tuple->GenTree->Branch("Event", &tuple->GenTree_Event, "Event/l");
    tuple->GenTree->Branch("Lumi", &tuple->GenTree_Lumi, "Lumi/i");
    /*tuple->GenTree->Branch("Hscp"    ,&tuple->GenTree_Hscp      ,"Hscp/i");*/
    tuple->GenTree->Branch("Weight", &tuple->GenTree_Weight, "Weight/F");
    tuple->GenTree->Branch("GeneratorWeight", &tuple->GenTree_GeneratorWeight, "GeneratorWeight/F");
    tuple->GenTree->Branch("GeneratorBinningValues", &tuple->GenTree_GeneratorBinningValues, "GeneratorBinningValues/F");
    tuple->GenTree->Branch("GenId", &tuple->GenTree_GenId);
    tuple->GenTree->Branch("GenCharge", &tuple->GenTree_GenCharge);
    tuple->GenTree->Branch("GenMass", &tuple->GenTree_GenMass);
    tuple->GenTree->Branch("GenPt", &tuple->GenTree_GenPt);
    tuple->GenTree->Branch("GenEta", &tuple->GenTree_GenEta);
    tuple->GenTree->Branch("GenPhi", &tuple->GenTree_GenPhi);
  }
}

//=============================================================
//
//      Initialize regions
//
//=============================================================

void TupleMaker::initializeRegions(Tuple *&tuple,
                                TFileDirectory &dir,
                                int etabins,
                                int ihbins,
                                int pbins,
                                int massbins) {
    tuple->rA_ias50.setSuffix("_regionA_ias50"); tuple->rA_ias50.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rC_ias50.setSuffix("_regionC_ias50"); tuple->rC_ias50.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rB_50ias60.setSuffix("_regionB_50ias60"); tuple->rB_50ias60.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rB_60ias70.setSuffix("_regionB_60ias70"); tuple->rB_60ias70.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rB_70ias80.setSuffix("_regionB_70ias80"); tuple->rB_70ias80.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rB_80ias90.setSuffix("_regionB_80ias90"); tuple->rB_80ias90.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rB_50ias90.setSuffix("_regionB_50ias90"); tuple->rB_50ias90.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rB_90ias100.setSuffix("_regionB_90ias100"); tuple->rB_90ias100.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rD_50ias60.setSuffix("_regionD_50ias60"); tuple->rD_50ias60.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rD_60ias70.setSuffix("_regionD_60ias70"); tuple->rD_60ias70.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rD_70ias80.setSuffix("_regionD_70ias80"); tuple->rD_70ias80.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rD_80ias90.setSuffix("_regionD_80ias90"); tuple->rD_80ias90.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rD_50ias90.setSuffix("_regionD_50ias90"); tuple->rD_50ias90.initHisto(dir, etabins, ihbins, pbins, massbins);
    tuple->rD_90ias100.setSuffix("_regionD_90ias100"); tuple->rD_90ias100.initHisto(dir, etabins, ihbins, pbins, massbins);
}

//=============================================================
//
//     Fill branches of Tree
//
//=============================================================

void TupleMaker::fillTreeBranches(Tuple *&tuple,
                                  const unsigned int &Trig,
                                  const unsigned int &Run,
                                  const unsigned long &Event,
                                  const unsigned int &Lumi,
                                  const unsigned int &PileUp,
                                  const unsigned int &nofVertices,
                                  const unsigned int &Hscp,
                                  const unsigned int &nmuons,
                                  const unsigned int &njets,
                                  const float &weight,
                                  const float &generator_weight,
                                  const float &generator_binning_values,
                                  const bool &HLT_Mu50,
                                  const bool &HLT_PFMET120_PFMHT120_IDTight,
                                  const bool &HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                                  const bool &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                                  const bool &HLT_MET105_IsoTrk50,
                                  const float &CaloMET,
                                  const float &RecoPFMET,
                                  const float &RecoPFMHT,
                                  const float &HLTPFMET,
                                  const float &HLTPFMHT,
                                  const float &RecoPFMET_eta,
                                  const float &RecoPFMET_phi,
                                  const float &RecoPFMET_significance,
                                  const float &Muon1_Pt,
                                  const float &Muon1_eta,
                                  const float &Muon1_phi,
                                  const float &Muon2_Pt,
                                  const float &Muon2_eta,
                                  const float &Muon2_phi,
                                  const std::vector<float> &Jet_pt,
                                  const std::vector<float> &Jet_eta,
                                  const std::vector<float> &Jet_phi,
                                  const std::vector<float> &Jet_mass,
                                  const std::vector<float> &Jet_energy,
                                  const std::vector<float> &Jet_pdgId,
                                  const std::vector<float> &Jet_et,
                                  const std::vector<float> &Jet_chargedEmEnergyFraction,
                                  const std::vector<float> &Jet_neutralEmEnergyFraction,
                                  const std::vector<float> &vect_mT,
                                  const std::vector<bool> &passCutPt55,
                                  const std::vector<bool> &passPreselection,
                                  const std::vector<bool> &passSelection,
                                  const std::vector<float> &Charge,
                                  const std::vector<float> &Pt,
                                  const std::vector<float> &PtErr,
                                  const std::vector<float> &Ias,
                                  const std::vector<float> &Ias_noTIBnoTIDno3TEC,
                                  const std::vector<float> &Ias_PixelOnly,
                                  const std::vector<float> &Ias_StripOnly,
                                  const std::vector<float> &Ias_PixelOnly_noL1,
                                  const std::vector<float> &Ih,
                                  const std::vector<float> &Ick,
                                  const std::vector<float> &Fmip,
                                  const std::vector<float> &ProbXY,
                                  const std::vector<float> &ProbXY_noL1,
                                  const std::vector<float> &ProbQ,
                                  const std::vector<float> &ProbQ_noL1,
                                  const std::vector<float> &Ndof,
                                  const std::vector<float> &Chi2,
                                  const std::vector<int>   &QualityMask,
                                  const std::vector<bool>  &isHighPurity,
                                  const std::vector<float>  &EoverP,
                                  const std::vector<bool>  &isMuon,
                                  const std::vector<int>   &MuonSelector,
                                  const std::vector<bool>  &isElectron,
                                  const std::vector<bool>  &isChHadron,
                                  const std::vector<bool>  &isNeutHadron,
                                  const std::vector<float> &ECAL_energy,
                                  const std::vector<float> &HCAL_energy,
                                  const std::vector<float> &TOF,  //equal to invBeta
                                  const std::vector<float> &TOFErr,
                                  const std::vector<unsigned int> &TOF_ndof,
                                  const std::vector<float> &DTTOF,
                                  const std::vector<float> &DTTOFErr,
                                  const std::vector<unsigned int> &DTTOF_ndof,
                                  const std::vector<float> &CSCTOF,
                                  const std::vector<float> &CSCTOFErr,
                                  const std::vector<unsigned int> &CSCTOF_ndof,
                                  const std::vector<float> &Mass,
                                  const std::vector<float> &MassErr,
                                  const std::vector<float> &dZ,
                                  const std::vector<float> &dXY,
                                  const std::vector<float> &dR,
                                  const std::vector<float> &eta,
                                  const std::vector<float> &phi,
                                  const std::vector<unsigned int> &noh,
                                  const std::vector<unsigned int> &noph,
                                  const std::vector<float> &fovh,
                                  const std::vector<unsigned int> &nomh,
                                  const std::vector<float> &fovhd,
                                  const std::vector<unsigned int> &nom,
                                  const std::vector<float> &iso_TK,
                                  const std::vector<float> &iso_ECAL,
                                  const std::vector<float> &iso_HCAL,
                                  const std::vector<float> &PFMiniIso_relative,
                                  const std::vector<float> &PFMiniIso_wMuon_relative,
                                  const std::vector<float> &track_PFIsolationR005_sumChargedHadronPt,
                                  const std::vector<float> &track_PFIsolationR005_sumNeutralHadronPt,
                                  const std::vector<float> &track_PFIsolationR005_sumPhotonPt,
                                  const std::vector<float> &track_PFIsolationR005_sumPUPt,
                                  const std::vector<float> &track_PFIsolationR01_sumChargedHadronPt,
                                  const std::vector<float> &track_PFIsolationR01_sumNeutralHadronPt,
                                  const std::vector<float> &track_PFIsolationR01_sumPhotonPt,
                                  const std::vector<float> &track_PFIsolationR01_sumPUPt,
                                  const std::vector<float> &track_PFIsolationR03_sumChargedHadronPt,
                                  const std::vector<float> &track_PFIsolationR03_sumNeutralHadronPt,
                                  const std::vector<float> &track_PFIsolationR03_sumPhotonPt,
                                  const std::vector<float> &track_PFIsolationR03_sumPUPt,
                                  const std::vector<float> &track_PFIsolationR05_sumChargedHadronPt,
                                  const std::vector<float> &track_PFIsolationR05_sumNeutralHadronPt,
                                  const std::vector<float> &track_PFIsolationR05_sumPhotonPt,
                                  const std::vector<float> &track_PFIsolationR05_sumPUPt,
                                  const std::vector<float> &muon_PFIsolationR03_sumChargedHadronPt,
                                  const std::vector<float> &muon_PFIsolationR03_sumNeutralHadronPt,
                                  const std::vector<float> &muon_PFIsolationR03_sumPhotonPt,
                                  const std::vector<float> &muon_PFIsolationR03_sumPUPt,
                                  const std::vector<float> &Ih_noL1,
                                  const std::vector<float> &Ih_15drop,
                                  const std::vector<float> &Ih_StripOnly,
                                  const std::vector<float> &Ih_StripOnly_15drop,
                                  const std::vector<float> &Ih_PixelOnly_noL1,
                                  const std::vector<float> &Ih_SaturationCorrectionFromFits,
                                  const std::vector<std::vector<float>> &clust_charge,
                                  const std::vector<std::vector<float>> &clust_pathlength,
                                  const std::vector<std::vector<unsigned int>> &clust_nstrip,
                                  const std::vector<std::vector<bool>> &clust_sat254,
                                  const std::vector<std::vector<bool>> &clust_sat255,
                                  const std::vector<std::vector<uint32_t>> &clust_detid,
                                  const std::vector<std::vector<bool>> &clust_isStrip,
                                  const std::vector<std::vector<bool>> &clust_isPixel,
                                  const std::vector<float> &genid,
                                  const std::vector<float> &gencharge,
                                  const std::vector<float> &genmass,
                                  const std::vector<float> &genpt,
                                  const std::vector<float> &geneta,
                                  const std::vector<float> &genphi) {
  tuple->Tree_Trig = Trig;
  tuple->Tree_Run = Run;
  tuple->Tree_Event = Event;
  tuple->Tree_Lumi = Lumi;
  tuple->Tree_PileUp = PileUp;
  tuple->Tree_nofVertices = nofVertices;
  tuple->Tree_Hscp = Hscp;
  tuple->Tree_nmuons = nmuons;
  tuple->Tree_njets = njets;
  tuple->Tree_Weight = weight;
  tuple->Tree_GeneratorWeight = generator_weight;
  tuple->Tree_GeneratorBinningValues = generator_binning_values;
  tuple->Tree_HLT_Mu50 = HLT_Mu50;
  tuple->Tree_HLT_PFMET120_PFMHT120_IDTight = HLT_PFMET120_PFMHT120_IDTight;
  tuple->Tree_HLT_PFHT500_PFMET100_PFMHT100_IDTight = HLT_PFHT500_PFMET100_PFMHT100_IDTight;
  tuple->Tree_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  tuple->Tree_HLT_MET105_IsoTrk50 = HLT_MET105_IsoTrk50;
  tuple->Tree_CaloMET = CaloMET;
  tuple->Tree_RecoPFMET = RecoPFMET;
  tuple->Tree_RecoPFMHT = RecoPFMHT;
  tuple->Tree_HLTPFMET = HLTPFMET;
  tuple->Tree_HLTPFMHT = HLTPFMHT;
  tuple->Tree_RecoPFMET_eta = RecoPFMET_eta;
  tuple->Tree_RecoPFMET_phi = RecoPFMET_phi;
  tuple->Tree_RecoPFMET_significance = RecoPFMET_significance;
  tuple->Tree_Muon1_Pt = Muon1_Pt;
  tuple->Tree_Muon1_eta = Muon1_eta;
  tuple->Tree_Muon1_phi = Muon1_phi;
  tuple->Tree_Muon2_Pt = Muon2_Pt;
  tuple->Tree_Muon2_eta = Muon2_eta;
  tuple->Tree_Muon2_phi = Muon2_phi;
  tuple->Tree_jet_pt = Jet_pt;
  tuple->Tree_jet_eta = Jet_eta;
  tuple->Tree_jet_phi = Jet_phi;
  tuple->Tree_jet_mass = Jet_mass;
  tuple->Tree_jet_energy = Jet_energy;
  tuple->Tree_jet_pdgId = Jet_pdgId;
  tuple->Tree_jet_et = Jet_et;
  tuple->Tree_jet_chargedEmEnergyFraction = Jet_chargedEmEnergyFraction;
  tuple->Tree_jet_neutralEmEnergyFraction = Jet_neutralEmEnergyFraction;
  tuple->Tree_vect_mT = vect_mT;
  tuple->Tree_passCutPt55 = passCutPt55;
  tuple->Tree_passPreselection = passPreselection;
  tuple->Tree_passSelection = passSelection;
  tuple->Tree_Charge = Charge;
  tuple->Tree_Pt = Pt;
  tuple->Tree_PtErr = PtErr;
  tuple->Tree_Ias = Ias;
  tuple->Tree_Ias_noTIBnoTIDno3TEC = Ias_noTIBnoTIDno3TEC;
  tuple->Tree_Ias_PixelOnly = Ias_PixelOnly;
  tuple->Tree_Ias_StripOnly = Ias_StripOnly;
  tuple->Tree_Ias_PixelOnly_noL1 = Ias_PixelOnly_noL1;
  tuple->Tree_Ih = Ih;
  tuple->Tree_Ick = Ick;
  tuple->Tree_Fmip = Fmip;
  tuple->Tree_ProbXY = ProbXY;
  tuple->Tree_ProbXY_noL1 = ProbXY_noL1;
  tuple->Tree_ProbQ = ProbQ;
  tuple->Tree_ProbQ_noL1 = ProbQ_noL1;
  tuple->Tree_Ndof = Ndof;
  tuple->Tree_Chi2 = Chi2;
  tuple->Tree_QualityMask = QualityMask;
  tuple->Tree_isHighPurity = isHighPurity;
  tuple->Tree_EoverP = EoverP;
  tuple->Tree_isMuon = isMuon;
  tuple->Tree_Muon_selector = MuonSelector;
  tuple->Tree_isElectron = isElectron;
  tuple->Tree_isChHadron = isChHadron;
  tuple->Tree_isNeutHadron = isNeutHadron;
  tuple->Tree_ECAL_energy = ECAL_energy;
  tuple->Tree_HCAL_energy = HCAL_energy;
  tuple->Tree_TOF = TOF;
  tuple->Tree_TOFErr = TOFErr;
  tuple->Tree_TOF_ndof = TOF_ndof;
  tuple->Tree_DTTOF = DTTOF;
  tuple->Tree_DTTOFErr = DTTOFErr;
  tuple->Tree_DTTOF_ndof = DTTOF_ndof;
  tuple->Tree_CSCTOF = CSCTOF;
  tuple->Tree_CSCTOFErr = CSCTOFErr;
  tuple->Tree_CSCTOF_ndof = CSCTOF_ndof;
  tuple->Tree_Mass = Mass;
  tuple->Tree_MassErr = MassErr;
  tuple->Tree_dZ = dZ;
  tuple->Tree_dXY = dXY;
  tuple->Tree_dR = dR;
  tuple->Tree_eta = eta;
  tuple->Tree_phi = phi;
  tuple->Tree_NOH = noh;
  tuple->Tree_NOPH = noph;
  tuple->Tree_FOVH = fovh;
  tuple->Tree_NOMH = nomh;
  tuple->Tree_FOVHD = fovhd;
  tuple->Tree_NOM = nom;
  tuple->Tree_iso_TK = iso_TK;
  tuple->Tree_iso_ECAL = iso_ECAL;
  tuple->Tree_iso_HCAL = iso_HCAL;
  tuple->Tree_PFMiniIso_relative = PFMiniIso_relative;
  tuple->Tree_PFMiniIso_wMuon_relative = PFMiniIso_wMuon_relative;
  tuple->Tree_track_PFIsolationR005_sumChargedHadronPt = track_PFIsolationR005_sumChargedHadronPt;
  tuple->Tree_track_PFIsolationR005_sumNeutralHadronPt = track_PFIsolationR005_sumNeutralHadronPt;
  tuple->Tree_track_PFIsolationR005_sumPhotonPt = track_PFIsolationR005_sumPhotonPt;
  tuple->Tree_track_PFIsolationR005_sumPUPt = track_PFIsolationR005_sumPUPt;
  tuple->Tree_track_PFIsolationR01_sumChargedHadronPt = track_PFIsolationR01_sumChargedHadronPt;
  tuple->Tree_track_PFIsolationR01_sumNeutralHadronPt = track_PFIsolationR01_sumNeutralHadronPt;
  tuple->Tree_track_PFIsolationR01_sumPhotonPt = track_PFIsolationR01_sumPhotonPt;
  tuple->Tree_track_PFIsolationR01_sumPUPt = track_PFIsolationR01_sumPUPt;
  tuple->Tree_track_PFIsolationR03_sumChargedHadronPt = track_PFIsolationR03_sumChargedHadronPt;
  tuple->Tree_track_PFIsolationR03_sumNeutralHadronPt = track_PFIsolationR03_sumNeutralHadronPt;
  tuple->Tree_track_PFIsolationR03_sumPhotonPt = track_PFIsolationR03_sumPhotonPt;
  tuple->Tree_track_PFIsolationR03_sumPUPt = track_PFIsolationR03_sumPUPt;
  tuple->Tree_track_PFIsolationR05_sumChargedHadronPt = track_PFIsolationR05_sumChargedHadronPt;
  tuple->Tree_track_PFIsolationR05_sumNeutralHadronPt = track_PFIsolationR05_sumNeutralHadronPt;
  tuple->Tree_track_PFIsolationR05_sumPhotonPt = track_PFIsolationR05_sumPhotonPt;
  tuple->Tree_track_PFIsolationR05_sumPUPt = track_PFIsolationR05_sumPUPt;
  tuple->Tree_muon_PFIsolationR03_sumChargedHadronPt = muon_PFIsolationR03_sumChargedHadronPt;
  tuple->Tree_muon_PFIsolationR03_sumNeutralHadronPt = muon_PFIsolationR03_sumNeutralHadronPt;
  tuple->Tree_muon_PFIsolationR03_sumPhotonPt = muon_PFIsolationR03_sumPhotonPt;
  tuple->Tree_muon_PFIsolationR03_sumPUPt = muon_PFIsolationR03_sumPUPt;
  tuple->Tree_Ih_noL1 = Ih_noL1;
  tuple->Tree_Ih_15drop = Ih_15drop;
  tuple->Tree_Ih_StripOnly = Ih_StripOnly;
  tuple->Tree_Ih_StripOnly_15drop = Ih_StripOnly_15drop;
  tuple->Tree_Ih_PixelOnly_noL1 = Ih_PixelOnly_noL1;
  tuple->Tree_Ih_SaturationCorrectionFromFits = Ih_SaturationCorrectionFromFits;
  tuple->Tree_clust_charge = clust_charge;
  tuple->Tree_clust_pathlength = clust_pathlength;
  tuple->Tree_clust_nstrip = clust_nstrip;
  tuple->Tree_clust_sat254 = clust_sat254;
  tuple->Tree_clust_sat255 = clust_sat255;
  tuple->Tree_clust_detid = clust_detid;
  tuple->Tree_clust_isStrip = clust_isStrip;
  tuple->Tree_clust_isPixel = clust_isPixel;
  tuple->Tree_GenId = genid;
  tuple->Tree_GenCharge = gencharge;
  tuple->Tree_GenMass = genmass;
  tuple->Tree_GenPt = genpt;
  tuple->Tree_GenEta = geneta;
  tuple->Tree_GenPhi = genphi;

  // Save in the tree
  tuple->Tree->Fill();
  //if (!SkipSelectionPlot_) tuple->Tree->Fill();
}

void TupleMaker::fillGenTreeBranches(Tuple *&tuple,
                                     const unsigned int &Run,
                                     const unsigned long &Event,
                                     const unsigned int &Lumi,
                                     /*const unsigned int &Hscp,*/
                                     const float &weight,
                                     const float &generator_weight,
                                     const float &generator_binning_values,
                                     const std::vector<float> &genid,
                                     const std::vector<float> &gencharge,
                                     const std::vector<float> &genmass,
                                     const std::vector<float> &genpt,
                                     const std::vector<float> &geneta,
                                     const std::vector<float> &genphi) {
  tuple->GenTree_Run = Run;
  tuple->GenTree_Event = Event;
  tuple->GenTree_Lumi = Lumi;
  /*tuple->GenTree_Hscp      = Hscp;*/
  tuple->GenTree_Weight = weight;
  tuple->GenTree_GeneratorWeight = generator_weight;
  tuple->GenTree_GeneratorBinningValues = generator_binning_values;
  tuple->GenTree_GenId = genid;
  tuple->GenTree_GenCharge = gencharge;
  tuple->GenTree_GenMass = genmass;
  tuple->GenTree_GenPt = genpt;
  tuple->GenTree_GenEta = geneta;
  tuple->GenTree_GenPhi = genphi;

  // Save in the tree
  tuple->GenTree->Fill();
}

//=============================================================
//
//     Fill of the ABCD related histograms -
//     -> this information will be used later in Step4 for the actual datadriven prediction
//
//=============================================================
void TupleMaker::fillControlAndPredictionHist(const susybsm::HSCParticle &hscp,
                                              const reco::DeDxData *dedxSObj,
                                              const reco::DeDxData *dedxMObj,
                                              const reco::MuonTimeExtra *tof,
                                              Tuple *&tuple,
                                              int TypeMode,
                                              float GlobalMinTOF,
                                              float Event_Weight,
                                              bool isCosmicSB,
                                              float DTRegion,
                                              const int MaxPredBins,
                                              float DeDxK,
                                              float DeDxC,
                                              std::vector<float> CutPt,
                                              std::vector<float> CutI,
                                              std::vector<float> CutTOF,
                                              std::vector<float> CutPt_Flip,
                                              std::vector<float> CutI_Flip,
                                              std::vector<float> CutTOF_Flip) {
  using namespace std;
  using namespace edm;


  reco::TrackRef track;
  if (TypeMode != 3)
    track = hscp.trackRef();
  else {
    reco::MuonRef muon = hscp.muonRef();
    if (muon.isNull())
      return;
    track = muon->standAloneMuon();
  }

  float MuonTOF = tof ? tof->inverseBeta() : GlobalMinTOF;

  float Is = 0;
  if (dedxSObj) {
    Is = dedxSObj->dEdx();
  }
  float Ih = 0;
  if (dedxMObj) {
    Ih = dedxMObj->dEdx();
  }

  if (!isCosmicSB) {
    tuple->Hist_Pt->Fill(track->pt(), Event_Weight);
    tuple->Hist_Is->Fill(Is, Event_Weight);
    tuple->Hist_TOF->Fill(MuonTOF, Event_Weight);
  }

  // std::cout << "Init After PT, Is and TOF plots are filled" << std::endl;
  //          /\ Ias
  //       /\  |----------------------------
  //        |  |   |           |             |
  //        |  |   |           |             |
  //        |  |   |    B      |     D       |
  //        |  |   |           |             |
  //        |  ------------------------------
  //        |  |   |           |             |
  //        |  |   |    A      |     C       |
  //        |  |   |           |             |
  //        |  |---|-----------|-------------|
  //        |  |   |           |             |
  //        |  /--------------------------------> PT
  //        | /       E       /    G
  //         /------------------------------->
  //        /
  //      TOF

  //Use different pt regions if using momentum from Stand Alone Muons
  std::vector<float> PtLimits;
  if (TypeMode != 3) {
    PtLimits.push_back(100);
    PtLimits.push_back(80);
    PtLimits.push_back(60);
  } else {
    PtLimits.push_back(240);
    PtLimits.push_back(170);
    PtLimits.push_back(120);
  }

  //Muon only prediction binned depending on where in the detector the track is and how many muon stations it has
  //Binning not used for other analyses
  int bin = -1;
  if (TypeMode == 3) {
    if (fabs(track->eta()) < DTRegion)
      bin = muonStations(track->hitPattern()) - 2;
    else
      bin = muonStations(track->hitPattern()) + 1;
  }

  if (!isCosmicSB) {
    if (track->pt() > PtLimits[0]) {
      tuple->CtrlPt_S4_Is->Fill(Is, Event_Weight);
      tuple->CtrlPt_S4_Ih->Fill(Ih, Event_Weight);
      if (tof)
        tuple->CtrlPt_S4_TOF->Fill(MuonTOF, Event_Weight);
      if (tof && bin >= 0 && bin < MaxPredBins)
        tuple->CtrlPt_S4_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
    } else if (track->pt() > PtLimits[1]) {
      tuple->CtrlPt_S3_Is->Fill(Is, Event_Weight);
      tuple->CtrlPt_S3_Ih->Fill(Ih, Event_Weight);
      if (tof)
        tuple->CtrlPt_S3_TOF->Fill(MuonTOF, Event_Weight);
      if (tof && bin >= 0 && bin < MaxPredBins)
        tuple->CtrlPt_S3_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
    } else if (track->pt() > PtLimits[2]) {
      tuple->CtrlPt_S2_Is->Fill(Is, Event_Weight);
      tuple->CtrlPt_S2_Ih->Fill(Ih, Event_Weight);
      if (tof)
        tuple->CtrlPt_S2_TOF->Fill(MuonTOF, Event_Weight);
      if (tof && bin >= 0 && bin < MaxPredBins)
        tuple->CtrlPt_S2_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
    } else {
      tuple->CtrlPt_S1_Is->Fill(Is, Event_Weight);
      tuple->CtrlPt_S1_Ih->Fill(Ih, Event_Weight);
      if (tof)
        tuple->CtrlPt_S1_TOF->Fill(MuonTOF, Event_Weight);
      if (tof && bin >= 0 && bin < MaxPredBins)
        tuple->CtrlPt_S1_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
    }

    if (Is > 0.2) {
      if (tof)
        tuple->CtrlIs_S4_TOF->Fill(MuonTOF, Event_Weight);
    } else if (Is > 0.1) {
      if (tof)
        tuple->CtrlIs_S3_TOF->Fill(MuonTOF, Event_Weight);
    } else if (Is > 0.05) {
      if (tof)
        tuple->CtrlIs_S2_TOF->Fill(MuonTOF, Event_Weight);
    } else {
      if (tof)
        tuple->CtrlIs_S1_TOF->Fill(MuonTOF, Event_Weight);
    }

    if (Ih > 4.4) {
      if (tof)
        tuple->CtrlIh_S4_TOF->Fill(MuonTOF, Event_Weight);
    } else if (Ih > 4.1) {
      if (tof)
        tuple->CtrlIh_S3_TOF->Fill(MuonTOF, Event_Weight);
    } else if (Ih > 3.8) {
      if (tof)
        tuple->CtrlIh_S2_TOF->Fill(MuonTOF, Event_Weight);
    } else {
      if (tof)
        tuple->CtrlIh_S1_TOF->Fill(MuonTOF, Event_Weight);
    }
  }

  //	 if(dedxMObj) Ih=dedxMObj->dEdx();
  float Ick = 0;
  if (dedxMObj)
    Ick = GetIck(Ih, DeDxK, DeDxC);  //GetIck(float I, bool MC, float dEdxK, float dEdxC)

  for (unsigned int CutIndex = 0; CutIndex < CutPt.size(); CutIndex++) {
    if (MuonTOF < GlobalMinTOF)
      continue;
    if (TypeMode == 5 && isCosmicSB)
      continue;
    bool PassPtCut = track->pt() >= CutPt[CutIndex];
    bool PassICut = (Is >= CutI[CutIndex]);
    bool PassTOFCut = MuonTOF >= CutTOF[CutIndex];

    if (PassTOFCut && PassPtCut && PassICut) {  //Region D
      tuple->H_D->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_D_Binned[to_string(bin)]->Fill(CutIndex, Event_Weight);
      tuple->RegionD_P->Fill(CutIndex, track->p(), Event_Weight);
      tuple->RegionD_I->Fill(CutIndex, Ih, Event_Weight);
      tuple->RegionD_Ias->Fill(CutIndex, Is, Event_Weight);
      tuple->RegionD_TOF->Fill(CutIndex, MuonTOF, Event_Weight);
      tuple->AS_Eta_RegionD->Fill(CutIndex, track->eta());
    } else if (PassTOFCut && PassPtCut && !PassICut) {  //Region C
      tuple->H_C->Fill(CutIndex, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
      tuple->PDF_C_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz
      //Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
      tuple->AS_Eta_RegionC->Fill(CutIndex, track->eta());
    } else if (PassTOFCut && !PassPtCut && PassICut) {  //Region B
      tuple->H_B->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_B_Binned[to_string(bin)]->Fill(CutIndex, Event_Weight);
      if (TypeMode < 2){
        tuple->Pred_I->Fill(CutIndex, Ih, Event_Weight);
        tuple->Pred_EtaI->Fill(CutIndex, track->eta(), Ih, Event_Weight);
      }
      if (TypeMode < 2)
        tuple->Pred_EtaS->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->PDF_B_EtaICK->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz
      //Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
      tuple->AS_Eta_RegionB->Fill(CutIndex, track->eta());
    } else if (PassTOFCut && !PassPtCut && !PassICut) {  //Region A
      tuple->H_A->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_TOF->Fill(CutIndex, MuonTOF, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaB->Fill(CutIndex, track->eta(), Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS2->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->AS_Eta_RegionA->Fill(CutIndex, track->eta());
      tuple->PDF_A_Eta->Fill(CutIndex, track->eta(), Event_Weight);  //pz

    } else if (!PassTOFCut && PassPtCut && PassICut) {  //Region H
      tuple->H_H->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_H_Binned[to_string(bin)]->Fill(CutIndex, Event_Weight);
      tuple->RegionH_Ias->Fill(CutIndex, Is, Event_Weight);
      if (TypeMode == 2 && Ick > 0)
        tuple->PDF_H_EtaMass->Fill(CutIndex, track->eta(), track->p() * sqrt(Ick), Event_Weight);  //pz
      //Pred_P->Fill(CutIndex,track->p(),        Event_Weight);
      //Pred_I->Fill(CutIndex,Ih,   Event_Weight);
      if (TypeMode == 2)
        tuple->AS_Eta_RegionH->Fill(CutIndex, track->eta());
    } else if (!PassTOFCut && PassPtCut && !PassICut) {  //Region G
      tuple->H_G->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
      tuple->AS_Eta_RegionG->Fill(CutIndex, track->eta());
      if (TypeMode == 2)
        tuple->PDF_G_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz
    } else if (!PassTOFCut && !PassPtCut && PassICut) {                             //Region F
      tuple->H_F->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_F_Binned[to_string(bin)]->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_I->Fill(CutIndex, Ih, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->AS_Eta_RegionF->Fill(CutIndex, track->eta());
      if (TypeMode == 2)
        tuple->PDF_F_EtaICK->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz

    } else if (!PassTOFCut && !PassPtCut && !PassICut) {  //Region E
      tuple->H_E->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaB->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->AS_Eta_RegionE->Fill(CutIndex, track->eta());
      if (TypeMode == 2)
        tuple->PDF_E_Eta->Fill(CutIndex, track->eta(), Event_Weight);  //pz
    }
  }

  //Use events with low TOF to check accuracy of background prediction
  for (unsigned int CutIndex = 0; CutIndex < CutPt_Flip.size(); CutIndex++) {
    if (TypeMode != 5 && MuonTOF >= GlobalMinTOF)
      continue;
    if (TypeMode == 5 && !isCosmicSB)
      continue;

    bool PassPtCut = track->pt() >= CutPt_Flip[CutIndex];
    bool PassICut = (Is >= CutI_Flip[CutIndex]);
    bool PassTOFCut = MuonTOF <= CutTOF_Flip[CutIndex];

    if (TypeMode == 5)
      PassTOFCut = true;

    if (PassTOFCut && PassPtCut && PassICut) {  //Region D
      tuple->RegionD_P_Flip->Fill(CutIndex, track->p(), Event_Weight);
      tuple->RegionD_I_Flip->Fill(CutIndex, Ih, Event_Weight);
      tuple->RegionD_Ias_Flip->Fill(CutIndex, Is, Event_Weight);
      tuple->RegionD_TOF_Flip->Fill(CutIndex, MuonTOF, Event_Weight);
      tuple->H_D_Flip->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_D_Binned_Flip[to_string(bin)]->Fill(CutIndex, Event_Weight);
    } else if (PassTOFCut && PassPtCut && !PassICut) {  //Region C
      tuple->H_C_Flip->Fill(CutIndex, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaP_Flip->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
      tuple->PDF_C_EtaP_Flip->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz
      //Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
    } else if (PassTOFCut && !PassPtCut && PassICut) {  //Region B
      tuple->H_B_Flip->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_B_Binned_Flip[to_string(bin)]->Fill(CutIndex, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_I_Flip->Fill(CutIndex, Ih, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaS_Flip->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->PDF_B_EtaICK_Flip->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz
      //Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
    } else if (PassTOFCut && !PassPtCut && !PassICut) {  //Region A
      tuple->H_A_Flip->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_TOF_Flip->Fill(CutIndex, MuonTOF, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaB_Flip->Fill(CutIndex, track->eta(), Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS2_Flip->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->PDF_A_Eta_Flip->Fill(CutIndex, track->eta(), Event_Weight);  //pz
    } else if (!PassTOFCut && PassPtCut && PassICut) {                    //Region H
      tuple->H_H_Flip->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_H_Binned_Flip[to_string(bin)]->Fill(CutIndex, Event_Weight);
      tuple->RegionH_Ias_Flip->Fill(CutIndex, Is, Event_Weight);
      if (TypeMode == 2 && Ick > 0)
        tuple->PDF_H_EtaMass_Flip->Fill(CutIndex, track->eta(), track->p() * sqrt(Ick), Event_Weight);  //pz

      //Pred_P_Flip->Fill(CutIndex,track->p(),        Event_Weight);
      //Pred_I_Flip->Fill(CutIndex,Ih,   Event_Weight);
    } else if (!PassTOFCut && PassPtCut && !PassICut) {  //Region G
      tuple->H_G_Flip->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaP_Flip->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
      if (TypeMode == 2)
        tuple->PDF_G_EtaP_Flip->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz

    } else if (!PassTOFCut && !PassPtCut && PassICut) {  //Region F
      tuple->H_F_Flip->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_F_Binned_Flip[to_string(bin)]->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_I_Flip->Fill(CutIndex, Ih, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS_Flip->Fill(CutIndex, track->eta(), Event_Weight);
      if (TypeMode == 2)
        tuple->PDF_F_EtaICK_Flip->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz
    } else if (!PassTOFCut && !PassPtCut && !PassICut) {                            //Region E
      tuple->H_E_Flip->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaB_Flip->Fill(CutIndex, track->eta(), Event_Weight);
      if (TypeMode == 2)
        tuple->PDF_E_Eta_Flip->Fill(CutIndex, track->eta(), Event_Weight);  //pz
    }
  }
}

//=============================================================
//
//  Fill regions used to validate the background estimate method
//
//=============================================================

void TupleMaker::fillRegions(Tuple *&tuple,
                             float pt_cut,
                             float Ias_quantiles[6],
                             float eta,
                             float p,
                             float pt,
                             float pterr,
                             float ih,
                             float ias,
                             float m,
                             float tof,
                             float w){
    if(pt<=pt_cut){
        if(ias<Ias_quantiles[0]) tuple->rA_ias50.fill(eta,p,pt,pterr,ih,ias,m,tof,w); 
        if(ias>=Ias_quantiles[0] && ias<Ias_quantiles[1]) tuple->rB_50ias60.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[1] && ias<Ias_quantiles[2]) tuple->rB_60ias70.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[2] && ias<Ias_quantiles[3]) tuple->rB_70ias80.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[3] && ias<Ias_quantiles[4]) tuple->rB_80ias90.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[0] && ias<Ias_quantiles[4]) tuple->rB_50ias90.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[4] && ias<Ias_quantiles[5]) tuple->rB_90ias100.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
    }else{
        if(ias<Ias_quantiles[0]) tuple->rC_ias50.fill(eta,p,pt,pterr,ih,ias,m,tof,w); 
        if(ias>=Ias_quantiles[0] && ias<Ias_quantiles[1]) tuple->rD_50ias60.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[1] && ias<Ias_quantiles[2]) tuple->rD_60ias70.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[2] && ias<Ias_quantiles[3]) tuple->rD_70ias80.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[3] && ias<Ias_quantiles[4]) tuple->rD_80ias90.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[0] && ias<Ias_quantiles[4]) tuple->rD_50ias90.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
        if(ias>=Ias_quantiles[4] && ias<Ias_quantiles[5]) {m<500?m=m:m=-1; tuple->rD_90ias100.fill(eta,p,pt,pterr,ih,ias,m,tof,w);} //blind in the last quantile
    }
}

void TupleMaker::writeRegions(Tuple *&tuple,
                              TFileDirectory &dir){
   dir.cd();
   tuple->rA_ias50.write();
   tuple->rC_ias50.write();
   tuple->rB_50ias60.write();
   tuple->rB_60ias70.write();
   tuple->rB_70ias80.write();
   tuple->rB_80ias90.write();
   tuple->rB_50ias90.write();
   tuple->rB_90ias100.write();
   tuple->rD_50ias60.write();
   tuple->rD_60ias70.write();
   tuple->rD_70ias80.write();
   tuple->rD_80ias90.write();
   tuple->rD_50ias90.write();
   tuple->rD_90ias100.write();
}
