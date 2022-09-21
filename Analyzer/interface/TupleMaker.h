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
  
  TH1::SetDefaultSumw2(kTRUE);

  tuple->IntLumi = dir.make<TProfile>("IntLumi", ";IntLumi", 1, 0, 1);
  tuple->XSection = dir.make<TProfile>("XSection", ";XSection", 1, 0, 1);
  
  tuple->NumEvents = dir.make<TH1F>("NumEvents",";;Number of events / category", 4, -0.5, 3.5);
  tuple->NumEvents->GetXaxis()->SetBinLabel(1,"(Re-weighted) events");
  tuple->NumEvents->GetXaxis()->SetBinLabel(2,"(Re-weighted) events w/ PU syst");
  tuple->NumEvents->GetXaxis()->SetBinLabel(3,"(Re-weighted) events w/ trig matching");
  
  tuple->dRMinHLTMuon = dir.make<TH1F>("dRMinHLTMuon",";#Delta R_{min,mu,HLT};Number of events/bin",100,0.,1.0);
  tuple->ErrorHisto = dir.make<TH1F>("ErrorHisto", ";;", 11, -0.5, 10.5);
  tuple->BefPreS_TriggerType = dir.make<TH1F>("BefPreS_TriggerType", ";;Events/category", 5, -0.5, 4.5);
  tuple->HSCPCandidateType = dir.make<TH1F>("HSCPCandidateType",";;Number of generator candidate / category", 6, -0.5, 5.5);
  tuple->HSCPCandidateType->GetXaxis()->SetBinLabel(1,"Neutral HSCP");
  tuple->HSCPCandidateType->GetXaxis()->SetBinLabel(2,"Single-charged");
  tuple->HSCPCandidateType->GetXaxis()->SetBinLabel(3,"Double-charged R-hadrons");
  tuple->HSCPCandidateType->GetXaxis()->SetBinLabel(4,"Tau-prime (1e or 2e)");
  tuple->HSCPCandidateType->GetXaxis()->SetBinLabel(5,"Else");
  
  tuple->BefPreS_RecoHSCParticleType = dir.make<TH1F>("BefPreS_RecoHSCParticleType",";;Track/category", 6, -0.5, 5.5);
  tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(1,"trackerMuon");
  tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(2,"matchedStandAloneMuon");
  tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(3,"standAloneMuon");
  tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(4,"innerTrack");
  tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(5,"unknown");
  // Can I do setBinLabel at this point?

  tuple->CutFlow = dir.make<TH1F>("CutFlow", ";CutFlowIndex", 17, 0., 17.);
  tuple->CutFlow->GetXaxis()->SetBinLabel(2,"Trigger");
  tuple->CutFlow->GetXaxis()->SetBinLabel(3,"p_{T}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(4,"#eta");
  tuple->CutFlow->GetXaxis()->SetBinLabel(5,"N_{no-L1 pixel hits}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(6,"N_{valid hit}/N_{all hits}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(7,"N_{dEdx hits}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(8,"HighPurity");
  tuple->CutFlow->GetXaxis()->SetBinLabel(9,"#chi^2 / N_{dof}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(10,"d_{z}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(11,"d_{xy}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(12,"MiniRelIsoAll");
  tuple->CutFlow->GetXaxis()->SetBinLabel(13,"MiniRelTkIso");
  tuple->CutFlow->GetXaxis()->SetBinLabel(14,"E/p");
  
  tuple->CutFlowReverse = dir.make<TH1F>("CutFlowReverse", ";CutFlowIndex", 17, -0.5, 16.5);
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(5,"N_{valid hit}/N_{all hits}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(8,"#chi^2 / N_{dof}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(13,"E/p");
  
  tuple->CutFlowProbQ =  dir.make<TH2F>("CutFlowProbQ",";ProbQ;",10, 0., 1.,17, -0.5, 16.5);
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(5,"N_{valid hit}/N_{all hits}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(8,"#chi^2 / N_{dof}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(13,"E/p");
  
  tuple->CutFlowEta = dir.make<TH2F>("CutFlowEta", ";#eta;", 50, -2.6, 2.6, 17, -0.5, 16.5);
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(5,"N_{valid hit}/N_{all hits}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(8,"#chi^2 / N_{dof}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(13,"E/p");
  
  tuple->CutFlowPfType = dir.make<TH2F>("CutFlowPfType", ";;", 9, -0.5, 8.5, 17, -0.5, 16.5);
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(5,"N_{valid hit}/N_{all hits}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(8,"#chi^2 / N_{dof}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(13,"E/p");
  
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(1,"AllTracks");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(2,"PFtracks");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(3,"isElectron");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(4,"isMuon");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(5,"isPhoton");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(6,"isChHadron");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(7,"#isNeutHadron");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(8,"isUndefined");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(9,"notPFtrack");

  tuple->N1_Eta = dir.make<TH1F>("N1_Eta",";#eta", 50, -2.6, 2.6);
  tuple->N1_Pt = dir.make<TH1F>("N1_Pt", ";p_{T} (GeV)", 50, 0, PtHistoUpperBound);
  tuple->N1_Pt_lowPt = dir.make<TH1F>("N1_Pt_lowPt", ";p_{T} (GeV)", 50, 0, 500);
  tuple->N1_Chi2oNdof = dir.make<TH1F>("N1_Chi2oNdof", ";#chi^2 / N_{dof}", 20, 0, 20);
  
  tuple->N1_Qual = dir.make<TH1F>("N1_Qual",";Tracks / category", 2, -0.5, 1.5);
  tuple->N1_Qual->GetXaxis()->SetBinLabel(1,"Not-HighPurity");
  tuple->N1_Qual->GetXaxis()->SetBinLabel(2,"HighPurity");
  
  tuple->N1_TNOM = dir.make<TH1F>("N1_TNOM", ";Number of measurments;Tracks / bin", 40, -0.5, 39.5);
  tuple->N1_TNOPH = dir.make<TH1F>("N1_TNOPH",";Number of pixel hits;Tracks / bin", 8, -0.5, 7.5);
  tuple->N1_TNOHFraction = dir.make<TH1F>("N1_TNOHFraction",";Number of valid hit fraction;Tracks / bin", 50, 0, 1);
//  tuple->N1_nDof = dir.make<TH1F>("nDof",";N_{dof}", 40, -0.5, 39.5);
//  tuple->N1_tofError = dir.make<TH1F>("tofError", ";tofError", 25, 0, 0.25);
  tuple->N1_TIsol = dir.make<TH1F>("TIsol", "TIsol", 25, 0, 100);
  tuple->N1_EoP = dir.make<TH1F>("N1_EoP",";Energy / Momentum", 25, 0, 1.5);
  tuple->N1_dRMinPfJet= dir.make<TH1F>("N1_dRMinPfJet",";dRMinPfJet",100,0.,1.5);
  tuple->N1_SumpTOverpT = dir.make<TH1F>("N1_SumpTOverpT",";SumpTOverpT", 80, 0, 2);
  tuple->N1_Ih = dir.make<TH1F>("N1_Ih",";I_{h}", 200, 0, dEdxM_UpLim);
  tuple->N1_MTOF = dir.make<TH1F>("N1_MTOF", ";TOF", 50, -2, 5);
  tuple->N1_I = dir.make<TH1F>("N1_I", ";I", NCuts, 0, NCuts);
  tuple->N1_TOF = dir.make<TH1F>("TOF", ";TOF", NCuts, 0, NCuts);

  tuple->NVTrack = dir.make<TH1F>("NVTrack", ";NVTrack", 1, 0, 1);
  tuple->N1_Stations = dir.make<TH1F>("N1_Stations", ";Stations", 1, 0, 1);
  tuple->N1_Dxy = dir.make<TH1F>("N1_Dxy",";d_{xy} (cm)", 200, -0.1, 0.1);
  tuple->N1_Dz = dir.make<TH1F>("N1_Dz",";d_{z} (cm)", 200, -0.3, 0.3);

  tuple->N1_PtErrOverPt = dir.make<TH1F>("N1_PtErrOverPt", ";#sigma_{p_{T}}/p_{T}", 40, 0, 1);
  tuple->N1_PtErrOverPt2 = dir.make<TH1F>("N1_PtErrOverPt2", ";#sigma_{p_{T}}/p_{T}^2", 40, 0, 1);
  tuple->N1_PtErrOverPtVsPt = dir.make<TH2F>("N1_PtErrOverPtVsPt", ";#sigma_{p_{T}}/p_{T};p_{T}",  40, 0., 1., 40, 0., 4000);

  tuple->N1_SegSep = dir.make<TH1F>("N1_SegSep", ";SegSep", 1, 0, 1);

  tuple->FailDz = dir.make<TH1F>("FailDz", ";FailDz", 1, 0, 1);

  tuple->N1_ProbQ = dir.make<TH1F>("N1_ProbQ", ";ProbQ", 100, 0, 1);
  tuple->N1_ProbQVsIas = dir.make<TH2F>("N1_ProbQVsIas",";ProbQ;I_{as}", 100, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->N1_ProbXY = dir.make<TH1F>("N1_ProbXY",";ProbXY", 100, 0, 1);
  tuple->N1_MiniRelIsoAll = dir.make<TH1F>("N1_MiniRelIsoAll", ";MiniRelIsoAll",  150, 0.0, 1.5);
  tuple->N1_MiniRelIsoAll_lowMiniRelIso = dir.make<TH1F>("N1_MiniRelIsoAll_lowMiniRelIso", ";MiniRelIsoAll",  100, 0.0, 0.1);
  tuple->N1_MiniRelTkIso = dir.make<TH1F>("N1_MiniRelTkIso", ";MiniRelTkIso",  150, 0.0, 1.5);
  tuple->N1_MiniRelTkIso_lowMiniRelIso = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso", ";MiniRelTkIso",  100, 0.0, 0.1);
  tuple->N1_MiniTkIso = dir.make<TH1F>("N1_MiniTkIso", ";MiniTkIso",  150, 0.0, 50.);
  tuple->N1_MiniRelTkIso_lowMiniRelIso_PUA = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso_PUA", ";MiniRelTkIso (PU < 15)",  100, 0.0, 0.1);
  tuple->N1_MiniRelTkIso_lowMiniRelIso_PUB = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso_PUB", ";MiniRelTkIso (15 =< PU < 30)",  100, 0.0, 0.1);
  tuple->N1_MiniRelTkIso_lowMiniRelIso_PUC = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso_PUC", ";MiniRelTkIso (PU >= 30)",  100, 0.0, 0.1);
  tuple->N1_MiniTkIso_PUA = dir.make<TH1F>("N1_MiniTkIso_PUA", ";MiniTkIso (PU < 15)",  150, 0.0, 50.);
  tuple->N1_MiniTkIso_PUB = dir.make<TH1F>("N1_MiniTkIso_PUB", ";MiniTkIso (15 =< PU < 30)",  150, 0.0, 50.);
  tuple->N1_MiniTkIso_PUC = dir.make<TH1F>("N1_MiniTkIso_PUC", ";MiniTkIso (PU >= 30)",  150, 0.0, 50.);

  tuple->N1_pfType = dir.make<TH1F>("N1_pfType", ";pfType", 9, 0, 9);

  tuple->HSCPE = dir.make<TH1F>("HSCPE", ";NCuts;HSCPE", NCuts, 0, NCuts);
  tuple->HSCPE_SystP = dir.make<TH1F>("HSCPE_SystP", ";NCuts;HSCPE_SystP", NCuts, 0, NCuts);
  tuple->HSCPE_SystI = dir.make<TH1F>("HSCPE_SystI", ";NCuts;HSCPE_SystI", NCuts, 0, NCuts);
  tuple->HSCPE_SystM = dir.make<TH1F>("HSCPE_SystM", ";NCuts;HSCPE_SystM", NCuts, 0, NCuts);
  tuple->HSCPE_SystT = dir.make<TH1F>("HSCPE_SystT", ";NCuts;HSCPE_SystT", NCuts, 0, NCuts);
  tuple->HSCPE_SystPU = dir.make<TH1F>("HSCPE_SystPU", ";NCuts;HSCPE_SystPU", NCuts, 0, NCuts);
  tuple->HSCPE_SystHUp = dir.make<TH1F>("HSCPE_SystHUp", ";NCuts;HSCPE_SystHUp", NCuts, 0, NCuts);
  tuple->HSCPE_SystHDown = dir.make<TH1F>("HSCPE_SystHDown", ";NCuts;HSCPE_SystHDown", NCuts, 0, NCuts);

  tuple->Mass = dir.make<TH2F>("Mass", ";NCuts;Mass", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF = dir.make<TH2F>("MassTOF", ";NCuts;MassTOF", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb = dir.make<TH2F>("MassComb", ";NCuts;MassComb", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass = dir.make<TH2F>("MaxEventMass", ";NCuts;MaxEventMass", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_SystP = dir.make<TH2F>("Mass_SystP", ";NCuts;Mass_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystP = dir.make<TH2F>("MassTOF_SystP", ";NCuts;MassTOF_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystP = dir.make<TH2F>("MassComb_SystP", ";NCuts;MassComb_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystP = dir.make<TH2F>("MaxEventMass_SystP", ";NCuts;MaxEventMass_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_SystI = dir.make<TH2F>("Mass_SystI", ";NCuts;Mass_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystI = dir.make<TH2F>("MassTOF_SystI", ";NCuts;MassTOF_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystI = dir.make<TH2F>("MassComb_SystI", ";NCuts;MassComb_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystI = dir.make<TH2F>("MaxEventMass_SystI", ";NCuts;MaxEventMass_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_SystM = dir.make<TH2F>("Mass_SystM", ";NCuts;Mass_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystM = dir.make<TH2F>("MassTOF_SystM", ";NCuts;MassTOF_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystM = dir.make<TH2F>("MassComb_SystM", ";NCuts;MassComb_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystM = dir.make<TH2F>("MaxEventMass_SystM", ";NCuts;MaxEventMass_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_SystT = dir.make<TH2F>("Mass_SystT", ";NCuts;Mass_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystT = dir.make<TH2F>("MassTOF_SystT", ";NCuts;MassTOF_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystT = dir.make<TH2F>("MassComb_SystT", ";NCuts;MassComb_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystT = dir.make<TH2F>("MaxEventMass_SystT", ";NCuts;MaxEventMass_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_SystPU = dir.make<TH2F>("Mass_SystPU", ";NCuts;Mass_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystPU = dir.make<TH2F>("MassTOF_SystPU", ";NCuts;MassTOF_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystPU = dir.make<TH2F>("MassComb_SystPU", ";NCuts;MassComb_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystPU = dir.make<TH2F>("MaxEventMass_SystPU", ";NCuts;MaxEventMass_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_SystHUp = dir.make<TH2F>("Mass_SystHUp", ";NCuts;Mass_SystHUp", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_SystH = dir.make<TH2F>("MassTOF_SystH", ";NCuts;MassTOF_SystH", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystHUp = dir.make<TH2F>("MassComb_SystHUp", ";NCuts;MassComb_SystHUp", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystHUp = dir.make<TH2F>("MaxEventMass_SystHUp", ";NCuts;MaxEventMass_SystHUp", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_SystHDown = dir.make<TH2F>("Mass_SystHDown", ";NCuts;Mass_SystHDown", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_SystHDown = dir.make<TH2F>("MassComb_SystHDown", ";NCuts;MassComb_SystHDown", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MaxEventMass_SystHDown = dir.make<TH2F>("MaxEventMass_SystHDown", ";NCuts;MaxEventMass_SystHDown", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  tuple->Mass_Flip = dir.make<TH2F>("Mass_Flip", ";NCuts;Mass_Flip", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassTOF_Flip = dir.make<TH2F>("MassTOF_Flip", ";NCuts;MassTOF_Flip", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
  tuple->MassComb_Flip = dir.make<TH2F>("MassComb_Flip", ";NCuts;MassComb_Flip", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);

  if (SkipSelectionPlot)
    return;

  tuple->Gen_DecayLength = dir.make<TH1F>("Gen_DecayLength", "Gen_DecayLength", 1000, 0, 1000);
  tuple->Beta_Gen = dir.make<TH1F>("Beta_Gen", "Beta_Gen", 100, 0, 1);
  tuple->Beta_GenCharged = dir.make<TH1F>("Beta_GenChaged", "Beta_GenChaged", 20, 0, 1);
  tuple->Beta_Triggered = dir.make<TH1F>("Beta_Triggered", "Beta_Triggered", 20, 0, 1);

  tuple->Beta_Matched = dir.make<TH1F>("Beta_Matched", "Beta_Matched", 20, 0, 1);
  tuple->Beta_PreselectedA = dir.make<TH1F>("Beta_PreselectedA", "Beta_PreselectedA", 20, 0, 1);
  tuple->Beta_PreselectedB = dir.make<TH1F>("Beta_PreselectedB", "Beta_PreselectedB", 20, 0, 1);
  tuple->Beta_PreselectedC = dir.make<TH1F>("Beta_PreselectedC", "Beta_PreselectedC", 20, 0, 1);


  tuple->BefPreS_GenPtVsdRMinGen = dir.make<TH2F>("BefPreS_GenPtVsdRMinGen", "BefPreS_GenPtVsdRMinGen", 50, 0, PtHistoUpperBound, 100, 0., 1.);
  tuple->BefPreS_GendRMin = dir.make<TH1F>("BefPreS_GendRMin",";dR_min;Gen candidate",100,0.,3.2);
  tuple->BefPreS_GenPtVsdRMinGenPostCut = dir.make<TH2F>("BefPreS_GenPtVsdRMinGenPostCut", ";GenPt (GeV);dRMinGen (after cut)", 50, 0, PtHistoUpperBound, 50, 0., 0.05);
  tuple->BefPreS_GenPtVsGenMinPt = dir.make<TH2F>("BefPreS_GenPtVsGenMinPt", "BefPreS_GenPtVsGenMinPt", 50, 0, PtHistoUpperBound, 100, 0, 1.);
  tuple->BefPreS_GenPtVsRecoPt = dir.make<TH2F>("BefPreS_GenPtVsRecoPt", "BefPreS_GenPtVsRecoPt", 50, 0, PtHistoUpperBound, 50, 0, PtHistoUpperBound);

  tuple->BefPreS_pfType = dir.make<TH1F>("BefPreS_pfType", "BefPreS_pfType", 9, 0, 9);

  tuple->BefPreS_massT = dir.make<TH1F>("BefPreS_massT", "BefPreS_massT", 50, 0.0, 250.0);

  tuple->BefPreS_MiniRelIsoAll = dir.make<TH1F>("BefPreS_MiniRelIsoAll",";MiniRelIsoAll;Tracks/bin", 150, 0.0, 1.5);
  tuple->BefPreS_MiniRelIsoChg = dir.make<TH1F>("BefPreS_MiniRelIsoChg",";MiniRelIsoChg;Tracks/bin",  150, 0.0, 1.5);
  tuple->BefPreS_MiniRelTkIso = dir.make<TH1F>("BefPreS_MiniRelTkIso",";MiniRelTkIso;Tracks/bin",  150, 0.0, 1.5);
  tuple->BefPreS_MiniTkIso = dir.make<TH1F>("BefPreS_MiniTkIso",";MiniTkIso;Tracks/bin",  150, 0.0, 50);

  tuple->BefPreS_RecoPFMET = dir.make<TH1F>("BefPreS_RecoPFMET", "BefPreS_RecoPFMET",  200, 0.0, 2000.0);
  tuple->BefPreS_RecoPFHT = dir.make<TH1F>("BefPreS_RecoPFHT", "BefPreS_RecoPFHT",  200, 0.0, 2000.0);
  tuple->BefPreS_CaloNumJets = dir.make<TH1F>("BefPreS_CaloNumJets", ";Number of calo jets;Jets/bin",  70, 0.0, 70.0);

  tuple->BefPreS_Chi2oNdof = dir.make<TH1F>("BefPreS_Chi2oNdof", "BefPreS_Chi2oNdof", 20, 0, 20);
  // This should just be a 2-bin plot where high-purity is not present = 0 or present = 1
  tuple->BefPreS_Qual = dir.make<TH1F>("BefPreS_Qual", "BefPreS_Qual", 2, -0.5, 1.5);
  tuple->BefPreS_TNOH_PUA = dir.make<TH1F>("BefPreS_TNOH_PUA", "BefPreS_TNOH_PUA",  40, -0.5, 39.5);
  tuple->BefPreS_TNOH_PUB = dir.make<TH1F>("BefPreS_TNOH_PUB", "BefPreS_TNOH_PUB", 40, -0.5, 39.5);
  tuple->BefPreS_TNOHFraction = dir.make<TH1F>("BefPreS_TNOHFraction", "BefPreS_TNOHFraction", 50, 0., 1.);
  tuple->BefPreS_TNOPH = dir.make<TH1F>("BefPreS_TNOPH",":Number of pixel hits", 8, -0.5, 7.5);
  tuple->BefPreS_TNOHFractionTillLast = dir.make<TH1F>("BefPreS_TNOHFractionTillLast", "BefPreS_TNOHFractionTillLast", 50, 0, 1);
  tuple->BefPreS_TNOMHTillLast = dir.make<TH1F>("BefPreS_TNOMHTillLast", "BefPreS_TNOMHTillLast", 20, 0, 20);
  tuple->BefPreS_Eta = dir.make<TH1F>("BefPreS_Eta", ";#eta", 50, -2.6, 2.6);
  tuple->BefPreS_TNOM = dir.make<TH1F>("BefPreS_TNOM",";Number of measurements", 40, -0.5, 39.5);
  tuple->BefPreS_TNOM_PUA = dir.make<TH1F>("BefPreS_TNOM_PUA",";Number of measurements (low PU)", 40, -0.5, 39.5);
  tuple->BefPreS_TNOM_PUB = dir.make<TH1F>("BefPreS_TNOM_PUB",";Number of measurements (high PU)", 40, -0.5, 39.5);
  tuple->BefPreS_nDof = dir.make<TH1F>("BefPreS_nDof",";Number of DF", 40, -0.5, 39.5);
  tuple->BefPreS_TOFError = dir.make<TH1F>("BefPreS_TOFError", "BefPreS_TOFError", 25, 0, 0.25);
  tuple->BefPreS_PtErrOverPt = dir.make<TH1F>("BefPreS_PtErrOverPt", ";#sigma_{p_{T}}/p_{T}", 40, 0, 1);
  tuple->BefPreS_PtErrOverPt2 = dir.make<TH1F>("BefPreS_PtErrOverPt2", ";#sigma_{p_{T}}/p_{T}^2", 40, 0, 0.003);
  tuple->BefPreS_Pt = dir.make<TH1F>("BefPreS_Pt", ";p_{T} (GeV)", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_lowPt = dir.make<TH1F>("BefPreS_Pt_lowPt", ";p_{T} (GeV)", 50, 0, 500);
  tuple->BefPreS_Ias = dir.make<TH1F>("BefPreS_Ias", "BefPreS_Ias", 10, 0., 1.);
  tuple->BefPreS_Ih = dir.make<TH1F>("BefPreS_Ih", "BefPreS_Ih", 200, 0, dEdxM_UpLim);
  tuple->BefPreS_MTOF = dir.make<TH1F>("BefPreS_MTOF", "BefPreS_MTOF", 50, -2, 5);
  tuple->BefPreS_TIsol = dir.make<TH1F>("BefPreS_TIsol", "BefPreS_TIsol", 25, 0, 100);
  tuple->BefPreS_EoP = dir.make<TH1F>("BefPreS_EoP", "BefPreS_EoP", 25, 0, 1.5);
  tuple->BefPreS_SumpTOverpT = dir.make<TH1F>("BefPreS_SumpTOverpT", "BefPreS_SumpTOverpT", 80, 0.0, 2.0);
  tuple->BefPreS_LastHitDXY = dir.make<TH1F>("BefPreS_LastHitDXY", "BefPreS_LastHitDXY", 75, 0, 150);
  tuple->BefPreS_LastHitD3D = dir.make<TH1F>("BefPreS_LastHitD3D", "BefPreS_LastHitD3D", 175, 0, 350);
  tuple->BefPreS_P = dir.make<TH1F>("BefPreS_P", "BefPreS_P", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt = dir.make<TH1F>("BefPreS_Pt", "BefPreS_Pt", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_PUA = dir.make<TH1F>("BefPreS_Pt_PUA", "BefPreS_Pt_PUA", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_PUB = dir.make<TH1F>("BefPreS_Pt_PUB", "BefPreS_Pt_PUB", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_Cosmic = dir.make<TH1F>("BefPreS_Pt_Cosmic", "BefPreS_Pt_Cosmic", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_DT = dir.make<TH1F>("BefPreS_Pt_DT", "BefPreS_Pt_DT", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_CSC = dir.make<TH1F>("BefPreS_Pt_CSC", "BefPreS_Pt_CSC", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Ias = dir.make<TH1F>("BefPreS_Ias", "BefPreS_Ias", 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ias_PUA = dir.make<TH1F>("BefPreS_Ias_PUA", "BefPreS_Ias_PUA", 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ias_PUB = dir.make<TH1F>("BefPreS_Ias_PUB", "BefPreS_Ias_PUB", 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ias_Cosmic = dir.make<TH1F>("BefPreS_Ias_Cosmic", "BefPreS_Ias_Cosmic", 100, 0, dEdxS_UpLim);
  tuple->BefPreS_Ih_Cosmic = dir.make<TH1F>("BefPreS_Ih_Cosmic", "BefPreS_Ih_Cosmic", 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih = dir.make<TH1F>("BefPreS_Ih", "BefPreS_Ih", 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih_PUA = dir.make<TH1F>("BefPreS_Ih_PUA", "BefPreS_Ih_PUA", 200, 0, dEdxM_UpLim);
  tuple->BefPreS_Ih_PUB = dir.make<TH1F>("BefPreS_Ih_PUB", "BefPreS_Ih_PUB", 200, 0, dEdxM_UpLim);
  tuple->BefPreS_TOF = dir.make<TH1F>("BefPreS_TOF", "BefPreS_TOF", 150, -1, 5);
  tuple->BefPreS_TOF_PUA = dir.make<TH1F>("BefPreS_TOF_PUA", "BefPreS_TOF_PUA", 150, -1, 5);
  tuple->BefPreS_TOF_PUB = dir.make<TH1F>("BefPreS_TOF_PUB", "BefPreS_TOF_PUB", 150, -1, 5);
  tuple->BefPreS_TOF_DT = dir.make<TH1F>("BefPreS_TOF_DT", "BefPreS_TOF_DT", 150, -1, 5);
  tuple->BefPreS_TOF_CSC = dir.make<TH1F>("BefPreS_TOF_CSC", "BefPreS_TOF_CSC", 150, -1, 5);
  tuple->BefPreS_dR_NVTrack = dir.make<TH1F>("BefPreS_dR_NVTrack", "BefPreS_dR_NVTrack", 40, 0, 1);
  tuple->BefPreS_MatchedStations = dir.make<TH1F>("BefPreS_MatchedStations", "BefPreS_MatchedStations", 8, -0.5, 7.5);
  tuple->BefPreS_InnerInvPtDiff = dir.make<TH1F>("BefPreS_InnerInvPtDiff", "BefPreS_InnerInvPtDiff", 120, -4, 4);
  tuple->BefPreS_Phi = dir.make<TH1F>("BefPreS_Phi", "BefPreS_Phi", 50, -3.14, 3.14);
  tuple->BefPreS_TimeAtIP = dir.make<TH1F>("BefPreS_TimeAtIP", "BefPreS_TimeAtIP", 50, -100, 100);
  tuple->BefPreS_OpenAngle = dir.make<TH1F>("BefPreS_OpenAngle", "BefPreS_OpenAngle", 50, -0.3, 3.15);
  tuple->BefPreS_OpenAngle_Cosmic = dir.make<TH1F>("BefPreS_OpenAngle_Cosmic", "BefPreS_OpenAngle_Cosmic", 50, -0.3, 3.15);

  tuple->BefPreS_NVertex = dir.make<TH1F>("BefPreS_NVertex", "BefPreS_NVertex", 50, -0.5, 49.5);
  tuple->BefPreS_NVertex_NoEventWeight = dir.make<TH1F>("BefPreS_NVertex_NoEventWeight", "BefPreS_NVertex_NoEventWeight", 50, -0.5, 49.5);
  tuple->BefPreS_PV = dir.make<TH1F>("BefPreS_PV", "BefPreS_PV", 60, 0, 60);
  tuple->BefPreS_PV_NoEventWeight = dir.make<TH1F>("BefPreS_PV_NoEventWeight", "BefPreS_PV_NoEventWeight", 60, 0, 60);
  tuple->BefPreS_NOMoNOH = dir.make<TH1F>("BefPreS_NOMoNOH",";Num of measurment / num of hits;Tracks/bin",10,0.,1.0);
  tuple->BefPreS_NOMoNOHvsPV = dir.make<TProfile>("BefPreS_NOMoNOHvsPV", "BefPreS_NOMoNOHvsPV", 60, 0, 60);
  tuple->BefPreS_dzAll = dir.make<TH1F>("BefPreS_dzAll", "BefPreS_dzAll", 200, -10, 10);
  tuple->BefPreS_dxyAll = dir.make<TH1F>("BefPreS_dxyAll", "BefPreS_dxyAll", 200, -0.2, 0.2);
  tuple->BefPreS_Dz = dir.make<TH1F>("BefPreS_Dz","d_{z} (cm)", 200, -0.3, 0.3);
  tuple->BefPreS_Dxy = dir.make<TH1F>("BefPreS_Dxy","d_{xy} (cm)", 200, -0.1, 0.1);

  tuple->BefPreS_SegSep = dir.make<TH1F>("BefPreS_SegSep", "BefPreS_SegSep", 50, 0, 2.5);
  tuple->BefPreS_SegMinEtaSep = dir.make<TH1F>("BefPreS_SegMinEtaSep", "BefPreS_SegMinEtaSep", 50, -1., 1.);
  tuple->BefPreS_SegMinPhiSep = dir.make<TH1F>("BefPreS_SegMinPhiSep", "BefPreS_SegMinPhiSep", 50, -3.3, 3.3);
  tuple->BefPreS_SegMinEtaSep_FailDz = dir.make<TH1F>("BefPreS_SegMinEtaSep_FailDz", "BefPreS_SegMinEtaSep_FailDz", 50, -1., 1.);
  tuple->BefPreS_SegMinEtaSep_PassDz = dir.make<TH1F>("BefPreS_SegMinEtaSep_PassDz", "BefPreS_SegMinEtaSep_PassDz", 50, -1., 1.);
  tuple->BefPreS_Dz_FailSep = dir.make<TH1F>("BefPreS_Dz_FailSep", "BefPreS_Dz_FailSep", 50, -150, 150);

  tuple->BefPreS_Dxy_Cosmic = dir.make<TH1F>("BefPreS_Dxy_Cosmic", "BefPreS_Dxy_Cosmic", 150, -IPbound, IPbound);
  tuple->BefPreS_Dz_Cosmic = dir.make<TH1F>("BefPreS_Dz_Cosmic", "BefPreS_Dz_Cosmic", 150, -IPbound, IPbound);
  tuple->BefPreS_Dz_CSC = dir.make<TH1F>("BefPreS_Dz_CSC", "BefPreS_Dz_CSC", 150, -IPbound, IPbound);
  tuple->BefPreS_Dz_DT = dir.make<TH1F>("BefPreS_Dz_DT", "BefPreS_Dz_DT", 150, -IPbound, IPbound);
  tuple->BefPreS_Pt_FailDz = dir.make<TH1F>("BefPreS_Pt_FailDz", "BefPreS_Pt_FailDz", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_FailDz_DT = dir.make<TH1F>("BefPreS_Pt_FailDz_DT", "BefPreS_Pt_FailDz_DT", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_Pt_FailDz_CSC = dir.make<TH1F>("BefPreS_Pt_FailDz_CSC", "BefPreS_Pt_FailDz_CSC", 50, 0, PtHistoUpperBound);
  tuple->BefPreS_TOF_FailDz = dir.make<TH1F>("BefPreS_TOF_FailDz", "BefPreS_TOF_FailDz", 150, -1, 5);
  tuple->BefPreS_TOF_FailDz_DT = dir.make<TH1F>("BefPreS_TOF_FailDz_DT", "BefPreS_TOF_FailDz_DT", 150, -1, 5);
  tuple->BefPreS_TOF_FailDz_CSC = dir.make<TH1F>("BefPreS_TOF_FailDz_CSC", "BefPreS_TOF_FailDz_CSC", 150, -1, 5);
  tuple->BefPreS_GenPtVsRecoPt = dir.make<TH2F>("BefPreS_GenPtVsRecoPt", "BefPreS_GenPtVsRecoPt", 50, 0, PtHistoUpperBound, 50, 0, PtHistoUpperBound);
  tuple->BefPreS_PtErrOverPtVsPt = dir.make<TH2F>("BefPreS_PtErrOverPtVsPt",  ";#sigma_{p_{T}}/p_{T};p_{T}",  40, 0., 1., 40, 0., 4000);
  tuple->BefPreS_PtErrOverPtVsPtErrOverPt2 = dir.make<TH2F>("BefPreS_PtErrOverPtVsPtErrOverPt2",  ";#sigma_{p_{T}}/p_{T};p_{T}^2",  40, 0., 1., 40, 0., 0.003);
  
  tuple->BefPreS_ProbQ = dir.make<TH1F>("BefPreS_ProbQ", "BefPreS_ProbQ", 100, 0, 1);
  tuple->BefPreS_ProbXY = dir.make<TH1F>("BefPreS_ProbXY", "BefPreS_ProbXY", 100, 0, 1);
  tuple->BefPreS_ProbQNoL1 = dir.make<TH1F>("BefPreS_ProbQNoL1", "BefPreS_ProbQNoL1", 100, 0, 1);
  tuple->BefPreS_ProbXYNoL1 = dir.make<TH1F>("BefPreS_ProbXYNoL1", "BefPreS_ProbXYNoL1", 100, 0, 1);
  tuple->BefPreS_MassErr = dir.make<TH1F>("BefPreS_MassErr", "BefPreS_MassErr", 50, 0., 10.);
  tuple->BefPreS_ProbQVsIas = dir.make<TH2F>("BefPreS_ProbQVsIas", "BefPreS_ProbQVsIas", 100, 0.0, 1.0, 100, 0.0, 1.0);
  
  tuple->BefPreS_EtaVsIas = dir.make<TH2F>("BefPreS_EtaVsIas", "BefPreS_EtaVsIas", 50, -3, 3, 10, 0., 1.);
  tuple->BefPreS_EtaVsIh = dir.make<TH2F>("BefPreS_EtaVsIh", "BefPreS_EtaVsIh", 50, -3, 3, 100, 0, dEdxM_UpLim);
  tuple->BefPreS_EtaVsP = dir.make<TH2F>("BefPreS_EtaVsP", "BefPreS_EtaVsP", 50, -3, 3, 50, 0, PtHistoUpperBound);
  tuple->BefPreS_EtaVsPt = dir.make<TH2F>("BefPreS_EtaVsPt", "BefPreS_EtaVsPt", 50, -3, 3, 50, 0, PtHistoUpperBound);
  tuple->BefPreS_EtaVsTOF = dir.make<TH2F>("BefPreS_EtaVsTOF", "BefPreS_EtaVsTOF", 50, -3, 3, 50, 0, 3);
  tuple->BefPreS_EtaVsNBH = dir.make<TH2F>("BefPreS_EtaVsNBH", "BefPreS_EtaVsNBH", 60, -3, 3, 24, 0, 24);
  tuple->BefPreS_EtaVsDz = dir.make<TH2F>("BefPreS_EtaVsDz", "BefPreS_EtaVsDz", 50, -3, 3, 50, -IPbound, IPbound);
  tuple->BefPreS_PVsIas = dir.make<TH2F>("BefPreS_PVsIas", "BefPreS_PVsIas", 50, 0, PtHistoUpperBound, 100, 0, dEdxS_UpLim);
  tuple->BefPreS_IhVsIas = dir.make<TH2F>("BefPreS_IhVsIas", "BefPreS_IhVsIas", 100, 0, dEdxM_UpLim, 100, 0, dEdxS_UpLim);
  tuple->BefPreS_PVsIh = dir.make<TH2F>("BefPreS_PVsIh", "BefPreS_PVsIh", 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  tuple->BefPreS_PtVsIas = dir.make<TH2F>("BefPreS_PtVsIas", "BefPreS_PtVsIas", 50, 0, PtHistoUpperBound, 10, 0., 1.);
  tuple->BefPreS_PtVsIh = dir.make<TH2F>("BefPreS_PtVsIh", "BefPreS_PtVsIh", 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  tuple->BefPreS_PtTOF = dir.make<TH2F>("BefPreS_PtTOF", "BefPreS_PtTOF", 50, 0, PtHistoUpperBound, 50, 0, 5);
    //tuple->BefPreS_TOFIs = dir.make<TH2F>("BefPreS_TOFIs", "BefPreS_TOFIs", 100, 1, 5, 100, 0, dEdxS_UpLim);
  tuple->BefPreS_TOFIs = dir.make<TH2F>("BefPreS_TOFIs", "BefPreS_TOFIs", 50, 0, 5, 10, 0., 1.);
    //tuple->BefPreS_TOFIm = dir.make<TH2F>("BefPreS_TOFIh", "BefPreS_TOFIh", 100, 1, 5, 200, 0, dEdxM_UpLim);
  tuple->BefPreS_TOFIh = dir.make<TH2F>("BefPreS_TOFIh", "BefPreS_TOFIh", 50, 0, 5, 100, 0, dEdxM_UpLim);

  tuple->BefPreS_CluProbQVsPixelLayer = dir.make<TH2F>("BefPreS_CluProbQVsPixelLayer",";CluProbQ;Layer",20,0.,1.,4,0.,4.);
  tuple->BefPreS_CluProbXYVsPixelLayer = dir.make<TH2F>("BefPreS_CluProbXYVsPixelLayer",";CluProbXY;Layer",100,0.,1.,4,0.,4.);
  tuple->BefPreS_CluNormChargeVsPixelLayer = dir.make<TH2F>("BefPreS_CluNormChargeVsPixelLayer",";CluNormCharge;Layer",100,0.,600.,4,0.,4.);
  tuple->BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma",";CluNormCharge;Layer",100,0.,600.,4,0.,4.);
  tuple->BefPreS_CluSizeVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeVsPixelLayer",";CluSize;Layer",10,0.,10.,4,0.,4.);
  tuple->BefPreS_CluSizeXVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeXVsPixelLayer",";CluSizeX;Layer",10,0.,10.,4,0.,4.);
  tuple->BefPreS_CluSizeYVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeYVsPixelLayer",";CluSizeY;Layer",10,0.,10.,4,0.,4.);
  tuple->BefPreS_CluSpecInCPEVsPixelLayer = dir.make<TH2F>("BefPreS_CluSpecInCPEVsPixelLayer",";CluSpecInCPE;Layer",4,0.,4.,4,0.,4.);

  tuple->BefPreS_CluCotBetaVsPixelLayer_lowProbXY = dir.make<TH2F>("BefPreS_CluCotBetaVsPixelLayer_lowProbXY",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->BefPreS_CluCotAlphaVsPixelLayer_lowProbXY = dir.make<TH2F>("BefPreS_CluCotAlphaVsPixelLayer_lowProbXY",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
  tuple->BefPreS_CluCotBetaVsPixelLayer = dir.make<TH2F>("BefPreS_CluCotBetaVsPixelLayer",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->BefPreS_CluCotAlphaVsPixelLayer = dir.make<TH2F>("BefPreS_CluCotAlphaVsPixelLayer",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);

  tuple->BefPreS_CluNormChargeVsStripLayer_lowBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_lowBetaGamma",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_Stat91 = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_Stat91",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatNot91 = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatNot91",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatHigherThan2 = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatHigherThan2",";CluNormCharge;Layer",600,0.,600.,20,0.,20.);

  tuple->BefPreS_dRMinPfJet= dir.make<TH1F>("BefPreS_dRMinPfJet",";dRMinPfJet",100,0.,1.5);
  tuple->BefPreS_dRMinPfJetVsIas =  dir.make<TH2F>("BefPreS_dRMinPfJetVsIas",";dRMinPfJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->BefPreS_dRMinCaloJet= dir.make<TH1F>("BefPreS_dRMinCaloJet",";dRMinCaloJet",100,0.,1.5);
  tuple->BefPreS_dRMinCaloJetVsIas =  dir.make<TH2F>("BefPreS_dRMinCaloJetVsIas",";dRMinCaloJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->BefPreS_genGammaBetaVsProbXYNoL1 =  dir.make<TH2F>("BefPreS_genGammaBetaVsProbXYNoL1", ";#gamma #beta;ProbXYNoL1",10,0.,1.3,20,0.,1.);
  tuple->BefPreS_dRVsPtPfJet = dir.make<TH2F>("BefPreS_dRVsPtPfJet",";dR(cand,jet);p_{T}",100,0.,1.5,100,0.,1000.);
  tuple->BefPreS_dRVsdPtPfCaloJet = dir.make<TH2F>("BefPreS_dRVsdPtPfCaloJet",";dRmin;dPtPfCaloJet",100,0.,1.5,20,0.,100.);
  

  tuple->PostPreS_TriggerType = dir.make<TH1F>("PostPreS_TriggerType", ";;Events/category", 4, -0.5, 3.5);
  tuple->PostPreS_RecoHSCParticleType = dir.make<TH1F>("PostPreS_RecoHSCParticleType",";;Track/category", 6, -0.5, 5.5);
  tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(1,"trackerMuon");
  tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(2,"matchedStandAloneMuon");
  tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(3,"standAloneMuon");
  tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(4,"innerTrack");
  tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(5,"unknown");

  tuple->PostPreS_pfType = dir.make<TH1F>("PostPreS_pfType", "PostPreS_pfType", 9, -0.5, 8.5);
  tuple->PostPreS_pfTypeVsIas = dir.make<TH2F>("PostPreS_pfTypeVsIas","PostPreS_pfTypeVsIas", 9, -0.5, 8.5,20,0.,1.);

  tuple->PostPreS_massT = dir.make<TH1F>("PostPreS_massT", "PostPreS_massT", 50, 0.0, 250.0);
  tuple->PostPreS_massTVsIas = dir.make<TH2F>("PostPreS_massTVsIas","PostPreS_massTVsIas",50, 0.0, 250.0, 20, 0., 1.);
  
  tuple->PostPreS_MiniRelIsoAll = dir.make<TH1F>("PostPreS_MiniRelIsoAll",";MiniRelIsoAll;Tracks/bin", 150, 0.0, 1.5);
  tuple->PostPreS_MiniRelIsoAllVsIas =  dir.make<TH2F>("PostPreS_MiniRelIsoAllVsIas","PostPreS_MiniRelIsoAllVsIas", 150, 0.0, 1.5, 10, 0.,1.);
  tuple->PostPreS_MiniRelIsoChg = dir.make<TH1F>("PostPreS_MiniRelIsoChg",";MiniRelIsoChg;Tracks/bin",  150, 0.0, 1.5);
  tuple->PostPreS_MiniTkIso = dir.make<TH1F>("PostPreS_MiniTkIso",";MiniTkIso;Tracks/bin", 150, 0.0, 50.);
  tuple->PostPreS_MiniRelTkIso = dir.make<TH1F>("PostPreS_MiniRelTkIso",";MiniRelTkIso;Tracks/bin", 150, 0.0, 1.5);
  
  tuple->PostPreS_RecoPFMET = dir.make<TH1F>("PostPreS_RecoPFMET", "PostPreS_RecoPFMET",  200, 0.0, 2000.0);
  tuple->PostPreS_RecoPFHT = dir.make<TH1F>("PostPreS_RecoPFHT", "PostPreS_RecoPFHT",  200, 0.0, 2000.0);
  tuple->PostPreS_CaloNumJets = dir.make<TH1F>("PostPreS_CaloNumJets",";Number of calo jets;Jets/bin",  30, 0.0, 30.0);
  
  tuple->PostPreS_Chi2oNdof = dir.make<TH1F>("PostPreS_Chi2oNdof", ";#chi^2/Ndof", 20, 0, 20);
  tuple->PostPreS_Chi2oNdofVsIas = dir.make<TH2F>("PostPreS_Chi2oNdofVsIas",";#chi^2/Ndof;I_{as}",20, 0, 20,10,0.,1.);
  tuple->PostPreS_Qual = dir.make<TH1F>("PostPreS_Qual","PostPreS_Qual", 2, -0.5, 1.5);
  tuple->PostPreS_TNOH_PUA = dir.make<TH1F>("PostPreS_TNOH_PUA", "Number of hits (low PU)", 40, -0.5, 39.5);
  tuple->PostPreS_TNOH_PUB = dir.make<TH1F>("PostPreS_TNOH_PUB", "Number of hits (high PU)", 40, -0.5, 39.5);
  tuple->PostPreS_TNOHFraction = dir.make<TH1F>("PostPreS_TNOHFraction", ";TNOHFraction", 20, 0, 1);
  tuple->PostPreS_TNOHFractionVsIas = dir.make<TH2F>("PostPreS_TNOHFractionVsIas","PostPreS_TNOHFractionVsIas",50, 0, 1,10,0.,1.);
  tuple->PostPreS_TNOPH = dir.make<TH1F>( "PostPreS_TNOPH",";Number of pixel hits", 8, -0.5, 7.5);
  tuple->PostPreS_TNOPHVsIas = dir.make<TH2F>("PostPreS_TNOPHVsIas",";_TNOPH;I_{as} ", 8,-0.5, 7.5, 20, 0., 1.);
  tuple->PostPreS_TNOHFractionTillLast = dir.make<TH1F>("PostPreS_TNOHFractionTillLast", ";TNOHFractionTillLast", 50, 0, 1);
  tuple->PostPreS_TNOMHTillLast = dir.make<TH1F>("PostPreS_TNOMHTillLast",";TNOMHTillLast", 20, -0.5, 19.5);
  tuple->PostPreS_Eta = dir.make<TH1F>("PostPreS_Eta", ";#eta", 50, -2.6, 2.6);
  tuple->PostPreS_EtaVsIas =  dir.make<TH2F>("PostPreS_EtaVsIas",";#eta;I_{as}", 50, -2.6, 2.6, 20,0.,1.);
  tuple->PostPreS_TNOM = dir.make<TH1F>("PostPreS_TNOM",";Number of measurement", 40, -0.5, 39.5);
  tuple->PostPreS_TNOMVsIas = dir.make<TH2F>("PostPreS_TNOMVsIas",";Number of measurement;I_{as}",  40, -0.5, 39.5, 20, 0., 1.);
  tuple->PostPreS_TNOM_PUA = dir.make<TH1F>("PostPreS_TNOM_PUA",  ";Number of measurement (low PU)", 40, -0.5, 39.5);
  tuple->PostPreS_TNOM_PUB = dir.make<TH1F>("PostPreS_TNOM_PUB",  ";Number of measurement (high PU)",  40, -0.5, 39.5);
  tuple->PostPreS_nDof = dir.make<TH1F>("PostPreS_nDof", ";nDof",  40, -0.5, 39.5);
  tuple->PostPreS_TOFError = dir.make<TH1F>("PostPreS_TOFError", "PostPreS_TOFError", 25, 0, 0.25);
  tuple->PostPreS_PtErrOverPt = dir.make<TH1F>("PostPreS_PtErrOverPt", ";#sigma_{p_{T}}/p_{T}", 40, 0, 1);
  tuple->PostPreS_PtErrOverPtVsIas =  dir.make<TH2F>("PostPreS_PtErrOverPtVsIas",";#sigma_{p_{T}}/p_{T};I_{as}", 40, 0, 1, 20, 0.,1.);
  tuple->PostPreS_PtErrOverPt2 = dir.make<TH1F>("PostPreS_PtErrOverPt2", ";#sigma_{p_{T}}/p_{T}^2", 40, 0, 0.003);
  tuple->PostPreS_Pt = dir.make<TH1F>("PostPreS_Pt", ";p_{T} (GeV)", 50, 0, PtHistoUpperBound);
  tuple->PostPreS_Pt_lowPt = dir.make<TH1F>("PostPreS_Pt_lowPt",";p_{T} (GeV)", 50, 0, 500);
  tuple->PostPreS_PtVsIas = dir.make<TH2F>("PostPreS_PtVsIas","PostPreS_PtVsIas", 50, 0, PtHistoUpperBound, 20, 0., 1.);
  tuple->PostPreS_Ias = dir.make<TH1F>("PostPreS_Ias", "PostPreS_Ias", 10, 0, dEdxS_UpLim);
  tuple->PostPreS_Ias_NoEventWeight = dir.make<TH1F>("PostPreS_Ias_NoEventWeight", "PostPreS_Ias_NoEventWeight", 10, 0, dEdxS_UpLim);
  tuple->PostPreS_Ih = dir.make<TH1F>("PostPreS_Ih", "PostPreS_Ih", 200, 0, dEdxM_UpLim);
  tuple->PostPreS_IhVsIas = dir.make<TH2F>("PostPreS_IhVsIas","PostPreS_IhVsIas",200, 0, dEdxM_UpLim, 20, 0.,1.);
  tuple->PostPreS_Ih_NoEventWeight = dir.make<TH1F>("PostPreS_Ih_NoEventWeight", "PostPreS_Ih_NoEventWeight", 200, 0, dEdxM_UpLim);
  tuple->PostPreS_MTOF = dir.make<TH1F>("PostPreS_MTOF", "PostPreS_MTOF", 50, -2, 5);
  tuple->PostPreS_TIsol = dir.make<TH1F>("PostPreS_TIsol", "PostPreS_TIsol", 25, 0, 100);
  tuple->PostPreS_TIsolVsIas = dir.make<TH2F>("PostPreS_TIsolVsIas","PostPreS_TIsolVsIas",25, 0, 100, 10, 0., 1.);
  tuple->PostPreS_EoP = dir.make<TH1F>("PostPreS_EoP", "PostPreS_EoP", 25, 0, 1.5);
  tuple->PostPreS_EoPVsIas = dir.make<TH2F>("PostPreS_EoPVsIas","PostPreS_EoPVsIas;E/ p;Ias",25, 0, 1.5, 20, 0.,1.);
  tuple->PostPreS_SumpTOverpT = dir.make<TH1F>("PostPreS_SumpTOverpT", "PostPreS_SumpTOverpT", 80, 0.0, 2.0);
  tuple->PostPreS_SumpTOverpTVsIas = dir.make<TH2F>("PostPreS_SumpTOverpTVsIas","PostPreS_SumpTOverpTVsIas;SumpTOverpT;Ias", 80, 0.0, 2.0, 20, 0.,1.);
  tuple->PostPreS_LastHitDXY = dir.make<TH1F>("PostPreS_LastHitDXY", "PostPreS_LastHitDXY", 75, 0, 150);
  tuple->PostPreS_LastHitDXYVsEta  = dir.make<TH2F>("PostPreS_LastHitDXYVsEta","PostPreS_LastHitDXYVsEta;LastHitDXY;Eta", 75, 0, 150, 20, 0.,1.);
  tuple->PostPreS_LastHitD3D = dir.make<TH1F>("PostPreS_LastHitD3D", "PostPreS_LastHitD3D", 175, 0, 350);
  tuple->PostPreS_LastHitD3DVsEta = dir.make<TH2F>("PostPreS_LastHitD3DVsEta","PostPreS_LastHitD3DVsEta;LastHitD3D;Eta", 175, 0, 350, 20, 0.,1.);
  tuple->PostPreS_P = dir.make<TH1F>("PostPreS_P", ";Momentum (GeV)", 50, 0, PtHistoUpperBound);
  tuple->PostPreS_dR_NVTrack = dir.make<TH1F>("PostPreS_dR_NVTrack", "PostPreS_dR_NVTrack", 40, 0, 1);
  tuple->PostPreS_MatchedStations = dir.make<TH1F>("PostPreS_MatchedStations",";MatchedStations",8, -0.5, 7.5);
  tuple->PostPreS_InnerInvPtDiff = dir.make<TH1F>("PostPreS_InnerInvPtDiff", "PostPreS_InnerInvPtDiff", 120, -4, 4);
  tuple->PostPreS_Phi = dir.make<TH1F>("PostPreS_Phi", "PostPreS_Phi", 50, -3.14, 3.14);
  tuple->PostPreS_TimeAtIP = dir.make<TH1F>("PostPreS_TimeAtIP", "PostPreS_TimeAtIP", 50, -100, 100);
  tuple->PostPreS_OpenAngle = dir.make<TH1F>("PostPreS_OpenAngle", "PostPreS_OpenAngle", 50, -0.3, 3.15);
  tuple->PostPreS_OpenAngle_Cosmic = dir.make<TH1F>("PostPreS_OpenAngle_Cosmic", "PostPreS_OpenAngle_Cosmic", 50, -0.3, 3.15);
  
  tuple->PostPreS_NVertex = dir.make<TH1F>("PostPreS_NVertex",";N_{vertex}", 50, -0.5, 49.5);
  tuple->PostPreS_NVertex_NoEventWeight = dir.make<TH1F>("PostPreS_NVertex_NoEventWeight",";N_{vertex} (NoEventWeight)", 50, -0.5, 49.5);
  tuple->PostPreS_PV = dir.make<TH1F>("PostPreS_PV",";PV", 60, -0.5, 59.5);
  tuple->PostPreS_PV_NoEventWeight = dir.make<TH1F>("PostPreS_PV_NoEventWeight", "PostPreS_PV_NoEventWeight", 60, 0, 60);
  tuple->PostPreS_NOMoNOH = dir.make<TH1F>("PostPreS_NOMoNOH",";Num of measurment / num of hits;Tracks/bin",10,0.,1.0);
  tuple->PostPreS_NOMoNOHvsPV = dir.make<TProfile>("PostPreS_NOMoNOHvsPV", "PostPreS_NOMoNOHvsPV", 60, 0, 60);
  tuple->PostPreS_dzAll = dir.make<TH1F>("PostPreS_dzAll", "PostPreS_dzAll", 200, -10, 10);
  tuple->PostPreS_dxyAll = dir.make<TH1F>("PostPreS_dxyAll", "PostPreS_dxyAll", 200, -0.2, 0.2);
  tuple->PostPreS_Dz = dir.make<TH1F>( "PostPreS_Dz", ";dz (cm)", 200, -0.3, 0.3);
  tuple->PostPreS_Dxy = dir.make<TH1F>("PostPreS_Dxy", "PostPreS_Dxy", 200, -0.1, 0.1);
  
  tuple->PostPreS_SegSep = dir.make<TH1F>("PostPreS_SegSep", "PostPreS_SegSep", 50, 0, 2.5);
  tuple->PostPreS_SegMinEtaSep = dir.make<TH1F>("PostPreS_SegMinEtaSep", "PostPreS_SegMinEtaSep", 50, -1., 1.);
  tuple->PostPreS_SegMinPhiSep = dir.make<TH1F>("PostPreS_SegMinPhiSep", "PostPreS_SegMinPhiSep", 50, -3.3, 3.3);
  tuple->PostPreS_SegMinEtaSep_FailDz = dir.make<TH1F>("PostPreS_SegMinEtaSep_FailDz", "PostPreS_SegMinEtaSep_FailDz", 50, -1., 1.);
  tuple->PostPreS_SegMinEtaSep_PassDz = dir.make<TH1F>("PostPreS_SegMinEtaSep_PassDz", "PostPreS_SegMinEtaSep_PassDz", 50, -1., 1.);
  tuple->PostPreS_Dz_FailSep = dir.make<TH1F>("PostPreS_Dz_FailSep", "PostPreS_Dz_FailSep", 50, -150, 150);
  
  tuple->PostPreS_Dxy_Cosmic = dir.make<TH1F>("PostPreS_Dxy_Cosmic", "PostPreS_Dxy_Cosmic", 150, -IPbound, IPbound);
  tuple->PostPreS_Dz_Cosmic = dir.make<TH1F>("PostPreS_Dz_Cosmic", "PostPreS_Dz_Cosmic", 150, -IPbound, IPbound);
  tuple->PostPreS_Dz_CSC = dir.make<TH1F>("PostPreS_Dz_CSC", "PostPreS_Dz_CSC", 150, -IPbound, IPbound);
  tuple->PostPreS_Dz_DT = dir.make<TH1F>("PostPreS_Dz_DT", "PostPreS_Dz_DT", 150, -IPbound, IPbound);
  tuple->PostPreS_Pt_FailDz = dir.make<TH1F>("PostPreS_Pt_FailDz", "PostPreS_Pt_FailDz", 50, 0, PtHistoUpperBound);
  tuple->PostPreS_Pt_FailDz_DT = dir.make<TH1F>("PostPreS_Pt_FailDz_DT", "PostPreS_Pt_FailDz_DT", 50, 0, PtHistoUpperBound);
  tuple->PostPreS_Pt_FailDz_CSC = dir.make<TH1F>("PostPreS_Pt_FailDz_CSC", "PostPreS_Pt_FailDz_CSC", 50, 0, PtHistoUpperBound);
  tuple->PostPreS_TOF_FailDz = dir.make<TH1F>("PostPreS_TOF_FailDz", "PostPreS_TOF_FailDz", 150, -1, 5);
  tuple->PostPreS_TOF_FailDz_DT = dir.make<TH1F>("PostPreS_TOF_FailDz_DT", "PostPreS_TOF_FailDz_DT", 150, -1, 5);
  tuple->PostPreS_TOF_FailDz_CSC = dir.make<TH1F>("PostPreS_TOF_FailDz_CSC", "PostPreS_TOF_FailDz_CSC", 150, -1, 5);
  tuple->PostPreS_PtErrOverPtVsPtErrOverPt2 = dir.make<TH2F>("PostPreS_PtErrOverPtVsPtErrOverPt2",  ";#sigma_{p_{T}}/p_{T};p_{T}^2",  40, 0., 1., 40, 0., 0.003);
  tuple->PostPreS_PtErrOverPtVsPt = dir.make<TH2F>("PostPreS_PtErrOverPtVsPt",  ";#sigma_{p_{T}}/p_{T};p_{T}",  40, 0., 1., 40, 0., 4000);
  tuple->PostPreS_GenPtVsRecoPt = dir.make<TH2F>("PostPreS_GenPtVsRecoPt", "PostPreS_GenPtVsRecoPt", 50, 0, PtHistoUpperBound, 50, 0, PtHistoUpperBound);
  
  tuple->PostPreS_ProbQ = dir.make<TH1F>("PostPreS_ProbQ", ";ProbQ", 20, 0., 1.);
  tuple->PostPreS_ProbQVsIas = dir.make<TH2F>("PostPreS_ProbQVsIas",";ProbQ;I_{as}", 20, 0., 1., 20, 0., 1.);
  tuple->PostPreS_ProbXY = dir.make<TH1F>("PostPreS_ProbXY", ";ProbXY", 100, 0, 1);
  tuple->PostPreS_ProbXY_highIas = dir.make<TH1F>("PostPreS_ProbXY_highIas",";ProbXY (I_{as} > 0.6)", 100, 0, 1);
  tuple->PostPreS_ProbXYVsIas = dir.make<TH2F>("PostPreS_ProbXYVsIas",";ProbXY;Ias",  100, 0, 1, 10, 0., 1.);
  tuple->PostPreS_ProbXYVsIas_highIas = dir.make<TH2F>("PostPreS_ProbXYVsIas_highIas",";ProbXY (I_{as} > 0.6);Ias (I_{as} > 0.6)",  100, 0, 1, 10, 0., 1.);
  tuple->PostPreS_ProbXYVsProbQ = dir.make<TH2F>("PostPreS_ProbXYVsProbQ",";ProbXY;ProbQ",  100, 0., 1., 10, 0., 1.);
  tuple->PostPreS_ProbXYVsProbQ_highIas = dir.make<TH2F>("PostPreS_ProbXYVsProbQ_highIas",";ProbXY (I_{as} > 0.6);ProbQ (I_{as} > 0.6)",  100, 0., 1., 10, 0., 1.);
  tuple->PostPreS_ProbQNoL1 = dir.make<TH1F>("PostPreS_ProbQNoL1",";ProbQNoL1", 20, 0., 1.);
  tuple->PostPreS_ProbQNoL1VsIas = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas",";ProbQNoL1;I_{as}",20, 0., 1., 20, 0., 1.);
  tuple->PostPreS_ProbXYNoL1 = dir.make<TH1F>("PostPreS_ProbXYNoL1", ";ProbXYNoL1", 100, 0, 1);
  tuple->PostPreS_ProbXYNoL1_highIas = dir.make<TH1F>("PostPreS_ProbXYNoL1_highIas", ";ProbXYNoL1 (I_{as} > 0.6)", 100, 0, 1);
  tuple->PostPreS_ProbXYNoL1VsIas = dir.make<TH2F>("PostPreS_ProbXYNoL1VsIas",";ProbXYNoL1;Ias", 100, 0., 1., 20, 0.,1.);
  tuple->PostPreS_ProbXYNoL1VsIas_highIas = dir.make<TH2F>("PostPreS_ProbXYNoL1VsIas_highIas",";ProbXYNoL1;I_{as} (I_{as} > 0.6)", 100, 0., 1., 20, 0.,1.);
  tuple->PostPreS_ProbXYNoL1VsProbQNoL1  = dir.make<TH2F>("PostPreS_ProbXYNoL1VsProbQNoL1",";ProbXYNoL1;ProbQNoL1", 100, 0., 1., 20,0.,1.);
  tuple->PostPreS_ProbXYNoL1VsProbQNoL1_highIas  = dir.make<TH2F>("PostPreS_ProbXYNoL1VsProbQNoL1_highIas",";ProbXYNoL1 (I_{as} > 0.6);ProbQNoL1 (I_{as} > 0.6)", 100, 0., 1., 20,0.,1.);

  tuple->PostPreS_MassErr = dir.make<TH1F>("PostPreS_MassErr", ";MassErr/Mass", 50, 0., 10.);
  tuple->PostPreS_MassErrVsIas = dir.make<TH2F>("PostPreS_MassErrVsIas",";MassErr/Mass;I_{as}",50, 0., 10.,10,0.,1.);

  tuple->PostPreS_EtaVsGenID = dir.make<TH2F>("PostPreS_EtaVsGenID", "PostPreS_EtaVsGenID",  50, -2.6, 2.6, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbQVsGenID = dir.make<TH2F>("PostPreS_ProbQVsGenID", "PostPreS_ProbQVsGenID", 100, 0.0, 1.0, 4000, 0.0, 4000.0);

  tuple->PostPreS_IasForStatus91 = dir.make<TH1F>("PostPreS_IasForStatus91",";I_{as} when status=91", 10, 0., 1.);
  tuple->PostPreS_IasForStatusNot91 = dir.make<TH1F>("PostPreS_IasForStatusNot91",";I_{as} when status!=91", 10, 0., 1.);

  tuple->PostPreS_ProbQVsGenEnviromentID = dir.make<TH2F>("PostPreS_ProbQVsGenEnviromentID",";ProbQ;GenEnviromentID",20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbXYVsGenID = dir.make<TH2F>("PostPreS_ProbXYVsGenID", "PostPreS_ProbXYVsGenID", 100, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_PtVsGenID = dir.make<TH2F>("PostPreS_PtVsGenID", "PostPreS_PtVsGenID", 50, 0, PtHistoUpperBound, 4000, 0.0, 4000.0);
  tuple->PostPreS_EoPVsGenID = dir.make<TH2F>("PostPreS_EoPVsGenID", "PostPreS_EoPVsGenID", 25, 0, 1.5, 4000, 0.0, 4000.0);
  tuple->PostPreS_IhVsGenID = dir.make<TH2F>("PostPreS_IhVsGenID", "PostPreS_IhVsGenID", 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsGenID = dir.make<TH2F>("PostPreS_IasVsGenID", "PostPreS_IasVsGenID", 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsGenEnviromentID = dir.make<TH2F>("PostPreS_IasVsGenEnviromentID",";Ias;GenEnviromentID", 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_massTVsGenID = dir.make<TH2F>("PostPreS_massTVsGenID", "PostPreS_massTVsGenID", 50, 0.0, 250.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoChgVsGenID = dir.make<TH2F>("PostPreS_miniIsoChgVsGenID", "PostPreS_miniIsoChgVsGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoAllVsGenID = dir.make<TH2F>("PostPreS_miniIsoAllVsGenID", "PostPreS_miniIsoAllVsGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_MassVsGenID = dir.make<TH2F>("PostPreS_MassVsGenID",";Mass (GeV);GenID",80,0.,4000.,4000, 0.0, 4000.0);

  tuple->PostPreS_EtaVsMomGenID = dir.make<TH2F>("PostPreS_EtaVsMomGenID", "PostPreS_EtaVsMomGenID",  50, -2.6, 2.6, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbQVsMomGenID = dir.make<TH2F>("PostPreS_ProbQVsMomGenID", "PostPreS_ProbQVsMomGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbXYVsMomGenID = dir.make<TH2F>("PostPreS_ProbXYVsMomGenID", "PostPreS_ProbXYVsMomGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_PtVsMomGenID = dir.make<TH2F>("PostPreS_PtVsMomGenID", "PostPreS_PtVsMomGenID", 50, 0, PtHistoUpperBound, 4000, 0.0, 4000.0);
  tuple->PostPreS_EoPVsMomGenID = dir.make<TH2F>("PostPreS_EoPVsMomGenID", "PostPreS_EoPVsMomGenID", 25, 0, 1.5, 4000, 0.0, 4000.0);
  tuple->PostPreS_IhVsMomGenID = dir.make<TH2F>("PostPreS_IhVsMomGenID", "PostPreS_IhVsMomGenID", 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsMomGenID = dir.make<TH2F>("PostPreS_IasVsMomGenID", "PostPreS_IasVsMomGenID", 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_massTVsMomGenID = dir.make<TH2F>("PostPreS_massTVsMomGenID", "PostPreS_massTVsMomGenID", 50, 0.0, 250.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoChgVsMomGenID = dir.make<TH2F>("PostPreS_miniIsoChgVsMomGenID", "PostPreS_miniIsoChgVsMomGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_miniIsoAllVsMomGenID = dir.make<TH2F>("PostPreS_miniIsoAllVsMomGenID", "PostPreS_miniIsoAllVsMomGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_MassVsMomGenID = dir.make<TH2F>("PostPreS_MassVsMomGenID","PostPreS_MassVsMomGenID;Mass;MomGenID",80,0.,4000.,4000, 0.0, 4000.0);
  
  tuple->PostPreS_EtaVsSiblingGenID = dir.make<TH2F>("PostPreS_EtaVsSiblingGenID", "PostPreS_EtaVsSiblingGenID",  50, -2.6, 2.6, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbQVsSiblingGenID = dir.make<TH2F>("PostPreS_ProbQVsSiblingGenID", "PostPreS_ProbQVsSiblingGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_ProbXYVsSiblingGenID = dir.make<TH2F>("PostPreS_ProbXYVsSiblingGenID", "PostPreS_ProbXYVsSiblingGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_PtVsSiblingGenID = dir.make<TH2F>("PostPreS_PtVsSiblingGenID", "PostPreS_PtVsSiblingGenID", 50, 0, PtHistoUpperBound, 4000, 0.0, 4000.0);
  tuple->PostPreS_EoPVsSiblingGenID = dir.make<TH2F>("PostPreS_EoPVsSiblingGenID", "PostPreS_EoPVsSiblingGenID", 25, 0, 1.5, 4000, 0.0, 4000.0);
  tuple->PostPreS_IhVsSiblingGenID = dir.make<TH2F>("PostPreS_IhVsSiblingGenID", "PostPreS_IhVsSiblingGenID", 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
  tuple->PostPreS_IasVsSiblingGenID = dir.make<TH2F>("PostPreS_IasVsSiblingGenID", "PostPreS_IasVsSiblingGenID", 10, 0., 1., 4000, 0.0, 4000.0);
  tuple->PostPreS_massTVsSiblingGenID = dir.make<TH2F>("PostPreS_massTVsSiblingGenID", ";massT;SiblingGenID", 50, 0.0, 250.0, 4000, 0.0, 4000.0);
  tuple->PostPreS_MassVsSiblingGenID = dir.make<TH2F>("PostPreS_MassVsSiblingGenID", ";Mass;SiblingGenID",80, 0.0, 4000.0, 4000, 0.0, 4000.0);
  
  tuple->PostPreS_EtaVsGenAngle = dir.make<TH2F>("PostPreS_EtaVsGenAngle", "PostPreS_EtaVsGenAngle",  50, -2.6, 2.6, 100, 0.0, 1.0);
  tuple->PostPreS_ProbQVsGenAngle = dir.make<TH2F>("PostPreS_ProbQVsGenAngle", "PostPreS_ProbQVsGenAngle", 20, 0.0, 1.0, 100, 0.0,1.0);
  tuple->PostPreS_ProbXYVsGenAngle = dir.make<TH2F>("PostPreS_ProbXYVsGenAngle", "PostPreS_ProbXYVsGenAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_PtVsGenAngle = dir.make<TH2F>("PostPreS_PtVsGenAngle", "PostPreS_PtVsGenAngle", 50, 0, PtHistoUpperBound, 100, 0.0, 1.0);
  tuple->PostPreS_EoPVsGenAngle = dir.make<TH2F>("PostPreS_EoPVsGenAngle", "PostPreS_EoPVsGenAngle", 25, 0, 1.5, 100, 0.0, 1.0);
  tuple->PostPreS_IhVsGenAngle = dir.make<TH2F>("PostPreS_IhVsGenAngle", "PostPreS_IhVsGenAngle", 200, 0, dEdxM_UpLim, 100, 0.0, 1.0);
  tuple->PostPreS_IasVsGenAngle = dir.make<TH2F>("PostPreS_IasVsGenAngle", "PostPreS_IasVsGenAngle", 10, 0., 1., 100, 0.0, 1.0);
  tuple->PostPreS_massTVsGenAngle = dir.make<TH2F>("PostPreS_massTVsGenAngle", "PostPreS_massTVsGenAngle", 50, 0.0, 250.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoChgVsGenAngle = dir.make<TH2F>("PostPreS_miniIsoChgVsGenAngle", "PostPreS_miniIsoChgVsGenAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoAllVsGenAngle = dir.make<TH2F>("PostPreS_miniIsoAllVsGenAngle", "PostPreS_miniIsoAllVsGenAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_MassVsGenAngle = dir.make<TH2F>("PostPreS_MassVsGenAngle",";Mass;GenAngle",80, 0.0,4000.0,  100, 0.0, 1.0);
  
  tuple->PostPreS_EtaVsGenMomAngle = dir.make<TH2F>("PostPreS_EtaVsGenMomAngle", "PostPreS_EtaVsGenMomAngle",  50, -2.6, 2.6, 100, 0.0, 1.0);
  tuple->PostPreS_ProbQVsGenMomAngle = dir.make<TH2F>("PostPreS_ProbQVsGenMomAngle", "PostPreS_ProbQVsGenMomAngle", 20, 0.0, 1.0, 100, 0.0,1.0);
  tuple->PostPreS_ProbXYVsGenMomAngle = dir.make<TH2F>("PostPreS_ProbXYVsGenMomAngle", "PostPreS_ProbXYVsGenMomAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_PtVsGenMomAngle = dir.make<TH2F>("PostPreS_PtVsGenMomAngle", "PostPreS_PtVsGenMomAngle", 50, 0, PtHistoUpperBound, 100, 0.0, 1.0);
  tuple->PostPreS_EoPVsGenMomAngle = dir.make<TH2F>("PostPreS_EoPVsGenMomAngle", "PostPreS_EoPVsGenMomAngle", 25, 0, 1.5, 100, 0.0, 1.0);
  tuple->PostPreS_IhVsGenMomAngle = dir.make<TH2F>("PostPreS_IhVsGenMomAngle", "PostPreS_IhVsGenMomAngle", 100, 0, dEdxM_UpLim, 100, 0.0, 1.0);
  tuple->PostPreS_IasVsGenMomAngle = dir.make<TH2F>("PostPreS_IasVsGenMomAngle", "PostPreS_IasVsGenMomAngle", 10, 0., 1., 100, 0.0, 1.0);
  tuple->PostPreS_massTVsGenMomAngle = dir.make<TH2F>("PostPreS_massTVsGenMomAngle", "PostPreS_massTVsGenMomAngle", 50, 0.0, 250.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoChgVsGenMomAngle = dir.make<TH2F>("PostPreS_miniIsoChgVsGenMomAngle", "PostPreS_miniIsoChgVsGenMomAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_miniIsoAllVsGenMomAngle = dir.make<TH2F>("PostPreS_miniIsoAllVsGenMomAngle", "PostPreS_miniIsoAllVsGenMomAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->PostPreS_MassVsGenMomAngle = dir.make<TH2F>("PostPreS_MassVsGenMomAngle",";Mass;GenMomAngle",80,0.,4000., 100, 0.0, 1.0);

  tuple->PostPreS_ProbQVsIas = dir.make<TH2F>("PostPreS_ProbQVsIas", "PostPreS_ProbQVsIas", 20, 0., 1., 20, 0., 1.);

  tuple->PostPreS_EtaVsGenNumSibling = dir.make<TH2F>("PostPreS_EtaVsGenNumSibling", "PostPreS_EtaVsGenNumSibling",  50, -2.6, 2.6, 100, 0.0, 10.0);
  tuple->PostPreS_ProbQVsGenNumSibling = dir.make<TH2F>("PostPreS_ProbQVsGenNumSibling", "PostPreS_ProbQVsGenNumSibling", 20, 0., 1., 10, 0.0,10.);
  tuple->PostPreS_ProbXYVsGenNumSibling = dir.make<TH2F>("PostPreS_ProbXYVsGenNumSibling", "PostPreS_ProbXYVsGenNumSibling", 20, 0.0, 1.0, 10, 0.0, 10.0);
  tuple->PostPreS_PtVsGenNumSibling = dir.make<TH2F>("PostPreS_PtVsGenNumSibling", "PostPreS_PtVsGenNumSibling", 50, 0, PtHistoUpperBound, 100, 0.0, 10.0);
  tuple->PostPreS_EoPVsGenNumSibling = dir.make<TH2F>("PostPreS_EoPVsGenNumSibling", "PostPreS_EoPVsGenNumSibling", 25, 0, 1.5, 100, 0.0, 10.0);
  tuple->PostPreS_IhVsGenNumSibling = dir.make<TH2F>("PostPreS_IhVsGenNumSibling", "PostPreS_IhVsGenNumSibling", 200, 0, dEdxM_UpLim, 100, 0.0, 10.0);
  tuple->PostPreS_IasVsGenNumSibling = dir.make<TH2F>("PostPreS_IasVsGenNumSibling", "PostPreS_IasVsGenNumSibling", 10, 0., 1., 100, 0.0, 10.0);
  tuple->PostPreS_massTVsGenNumSibling = dir.make<TH2F>("PostPreS_massTVsGenNumSibling", "PostPreS_massTVsGenNumSibling", 50, 0.0, 250.0, 100, 0.0, 10.0);
  tuple->PostPreS_miniIsoChgVsGenNumSibling = dir.make<TH2F>("PostPreS_miniIsoChgVsGenNumSibling", "PostPreS_miniIsoChgVsGenNumSibling", 20, 0.0, 1.0, 100, 0.0, 10.0);
  tuple->PostPreS_miniIsoAllVsGenNumSibling = dir.make<TH2F>("PostPreS_miniIsoAllVsGenNumSibling", "PostPreS_miniIsoAllVsGenNumSibling", 20, 0.0, 1.0, 100, 0.0, 10.0);

  tuple->PostPreS_EoPVsPfType = dir.make<TH2F>("PostPreS_EoPVsPfType", "PostPreS_EoPVsPfType", 25, 0.0, 1.5, 9, 0.0, 9.0);
  tuple->PostPreS_Mass = dir.make<TH1F>("PostPreS_Mass","PostPreS_Mass;Mass (GeV);Tracks / 50 GeV", 80,0.,4000.);
  tuple->PostPreS_MassVsPfType = dir.make<TH2F>("PostPreS_MassVsPfType","PostPreS_MassVsPfType;Mass (GeV);PF ID",80,0.,4000.,9, 0.0, 9.0);
  tuple->PostPreS_MassVsPt = dir.make<TH2F>("PostPreS_MassVsPt", ";Mass (GeV);", 80,0.,4000.,80,0.,4000.);
  tuple->PostPreS_MassVsP = dir.make<TH2F>("PostPreS_MassVsP", ";Mass (GeV);", 80,0.,4000.,80,0.,4000.);
  tuple->PostPreS_MassVsTNOHFraction = dir.make<TH2F>("PostPreS_MassVsTNOHFraction", ";Mass (GeV);", 80,0.,4000.,50, 0, 1);
  tuple->PostPreS_MassVsTNOPH = dir.make<TH2F>("PostPreS_MassVsTNOPH", ";Mass (GeV);", 80,0.,4000.,8, -0.5, 7.5);
  tuple->PostPreS_MassVsTNOM = dir.make<TH2F>("PostPreS_MassVsTNOM", ";Mass (GeV);", 80,0.,4000.,40, -0.5, 39.5);
  tuple->PostPreS_MassVsProbQNoL1 = dir.make<TH2F>("PostPreS_MassVsProbQNoL1", ";Mass (GeV);", 80,0.,4000.,20,0.,1.);
  tuple->PostPreS_MassVsProbXYNoL1 = dir.make<TH2F>("PostPreS_MassVsProbXYNoL1", ";Mass (GeV);", 80,0.,4000.,20,0.,1.);
  tuple->PostPreS_MassVsEoP = dir.make<TH2F>("PostPreS_MassVsEoP", ";Mass (GeV);", 80,0.,4000.,30, 0, 0.3);
  tuple->PostPreS_MassVsSumpTOverpT = dir.make<TH2F>("PostPreS_MassVsSumpTOverpT", ";Mass (GeV);", 80,0.,4000.,80, 0, 2);
  tuple->PostPreS_MassVsPtErrOverPt = dir.make<TH2F>("PostPreS_MassVsPtErrOverPt", ";Mass (GeV);", 80,0.,4000.,40, 0, 1);
  tuple->PostPreS_MassVsTIsol = dir.make<TH2F>("PostPreS_MassVsTIsol", ";Mass (GeV);", 80,0.,4000., 25, 0, 100);
  tuple->PostPreS_MassVsIh = dir.make<TH2F>("PostPreS_MassVsIh", ";Mass (GeV);", 80,0.,4000.,200, 0, dEdxM_UpLim);
  tuple->PostPreS_MassVsMassT = dir.make<TH2F>("PostPreS_MassVsMassT", ";Mass (GeV);", 80,0.,4000.,50, 0.0, 250.0);
  tuple->PostPreS_MassVsMiniRelIsoAll = dir.make<TH2F>("PostPreS_MassVsMiniRelIsoAll", ";Mass (GeV);", 80,0.,4000.,20, 0., 0.2);
  tuple->PostPreS_MassVsMassErr = dir.make<TH2F>("PostPreS_MassVsMassErr", ";Mass (GeV);", 80,0.,4000.,50, 0., 10.);
  
  // Maybe we dont need these anymore
  tuple->PostPreS_IasAllIhVsLayer = dir.make<TH3F>("PostPreS_IasAllIhVsLayer", "PostPreS_IasAllIhVsLayer", 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 35, 0.,35.);
  tuple->PostPreS_IasPixelIhVsLayer = dir.make<TH3F>("PostPreS_IasPixelIhVsLayer", "PostPreS_IasPixelIhVsLayer", 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 10, 0.,10.);
  tuple->PostPreS_IasStripIhVsLayer = dir.make<TH3F>("PostPreS_IasStripIhVsLayer", "PostPreS_IasStripIhVsLayer", 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 25, 0.,25.);

  tuple->PostPreS_CluProbQVsPixelLayer = dir.make<TH2F>("PostPreS_CluProbQVsPixelLayer",";CluProbQ;Layer",20,0.,1.,4,0.,4.);
  tuple->PostPreS_CluProbXYVsPixelLayer = dir.make<TH2F>("PostPreS_CluProbXYVsPixelLayer",";CluProbXY;Layer",100,0.,1.,4,0.,4.);
  tuple->PostPreS_CluSizeVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeVsPixelLayer",";CluSize;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeXVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeXVsPixelLayer",";CluSizeX;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeYVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeYVsPixelLayer",";CluSizeY;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSpecInCPEVsPixelLayer = dir.make<TH2F>("PostPreS_CluSpecInCPEVsPixelLayer",";CluSpecInCPE;Layer",4,0.,4.,4,0.,4.);

  tuple->PostPreS_CluProbQVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluProbQVsPixelLayer_highIas",";CluProbQ;Layer",20,0.,1.,4,0.,4.);
  tuple->PostPreS_CluProbXYVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluProbXYVsPixelLayer_highIas",";CluProbXY;Layer",100,0.,1.,4,0.,4.);
  tuple->PostPreS_CluSizeVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeVsPixelLayer_highIas",";CluSize;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeXVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeXVsPixelLayer_highIas",";CluSizeX;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSizeYVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeYVsPixelLayer_highIas",";CluSizeY;Layer",10,0.,10.,4,0.,4.);
  tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSpecInCPEVsPixelLayer_highIas",";CluSpecInCPE;Layer",4,0.,4.,4,0.,4.);

  tuple->PostPreS_CluCotBetaVsPixelLayer_lowProbXY = dir.make<TH2F>("PostPreS_CluCotBetaVsPixelLayer_lowProbXY",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->PostPreS_CluCotAlphaVsPixelLayer_lowProbXY = dir.make<TH2F>("PostPreS_CluCotAlphaVsPixelLayer_lowProbXY",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
  tuple->PostPreS_CluCotBetaVsPixelLayer = dir.make<TH2F>("PostPreS_CluCotBetaVsPixelLayer",";CotBeta;Layer",200,-10.,10.,4,0.,4.);
  tuple->PostPreS_CluCotAlphaVsPixelLayer = dir.make<TH2F>("PostPreS_CluCotAlphaVsPixelLayer",";CotAlpha;Layer",100,-1.,1.,4,0.,4.);

  tuple->PostPreS_CluNormChargeVsStripLayer_lowBetaGamma = dir.make<TH2F>("PostPreS_CluNormChargeVsStripLayer_lowBetaGamma",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->PostPreS_CluNormChargeVsStripLayer_higherBetaGamma = dir.make<TH2F>("PostPreS_CluNormChargeVsStripLayer_higherBetaGamma",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->PostPreS_CluNormChargeVsStripLayer_higherBetaGamma_Stat91 = dir.make<TH2F>("PostPreS_CluNormChargeVsStripLayer_higherBetaGamma_Stat91",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->PostPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatNot91 = dir.make<TH2F>("PostPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatNot91",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
  tuple->PostPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatHigherThan2 = dir.make<TH2F>("PostPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatHigherThan2",";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);

  tuple->PostPreS_dRMinPfJet = dir.make<TH1F>("PostPreS_dRMinPfJet",";dRMinPfJet",100,0.,1.5);
  tuple->PostPreS_closestPfJetMuonFraction = dir.make<TH1F>("PostPreS_closestPfJetMuonFraction",":closestPfJetMuonFraction",20,0.,1.);
  tuple->PostPreS_closestPfJetElectronFraction = dir.make<TH1F>("PostPreS_closestPfJetElectronFraction",";closestPfJetElectronFraction",20,0.,1.);
  tuple->PostPreS_closestPfJetPhotonFraction = dir.make<TH1F>("PostPreS_closestPfJetPhotonFraction",";closestPfJetPhotonFraction",20,0.,1.);

  tuple->PostPreS_closestPfJetMuonFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetMuonFractionVsIas",":closestPfJetMuonFraction;Ias",20,0.,1.,20,0.,1.);
  tuple->PostPreS_closestPfJetElectronFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetElectronFractionVsIas",";closestPfJetElectronFraction;Ias",20,0.,1.,20,0.,1.);
  tuple->PostPreS_closestPfJetPhotonFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetPhotonFractionVsIas",";closestPfJetPhotonFraction;Ias",20,0.,1.,20,0.,1.);

  tuple->PostPreS_dRMinPfJetVsIas = dir.make<TH2F>("PostPreS_dRMinPfJetVsIas", ";dRMinPfJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->PostPreS_dRMinCaloJet = dir.make<TH1F>("PostPreS_dRMinCaloJet",";dRMinCaloJet",100,0.,1.5);
  tuple->PostPreS_dPhiMinPfMet = dir.make<TH1F>("PostPreS_dPhiMinPfMet",";dPhiMinPfMet",100,0.,3.2);

  tuple->PostPreS_dRMinCaloJetVsIas =  dir.make<TH2F>("PostPreS_dRMinCaloJetVsIas",";dRMinCaloJet;Ias",100,0.,1.5,10,0.,1.);
  tuple->PostPreS_dPhiMinPfMetVsIas =  dir.make<TH2F>("PostPreS_dPhiMinPfMetVsIas",";dPhiMinPfMet;Ias",100,0.,3.2,10,0.,1.);

  tuple->PostPreS_PfMet = dir.make<TH1F>("PostPreS_PfMet",";PfMet",200,0.,2000.);
  tuple->PostPreS_PfMetPhi = dir.make<TH1F>("PostPreS_PfMetPhi",";PfMetPhi",30,0.,3.2);
  
  tuple->GenLevelBinning = dir.make<TH1F>("GenLevelBinning","GenLevelBinning",1200,0.,1200.);
  tuple->GenLevelpT = dir.make<TH1F>("GenLevelpT", "GenLevelpT;Generator p_{T} (GeV)", 50, 0, PtHistoUpperBound);
  tuple->GenLevelEta = dir.make<TH1F>("GenLevelEta",";Generator #eta", 60, -3, 3);
  tuple->GenLevelBeta = dir.make<TH1F>("GenLevelBeta",";Generator #beta", 20, 0, 1);
  tuple->GenLevelBetaGamma = dir.make<TH1F>("GenLevelBetaGamma",";Generator #beta #gamma;Gen canidate",4500,0.,450.);

  //Initialize histograms for number of bins.  For everything but muon only PredBins=0 so no histograms created
  for (int i = 0; i < PredBins; i++) {
    char Suffix[1024];
    sprintf(Suffix, "_%i", i);
      Name.append(Suffix);
    tuple->BefPreS_Pt_Binned[std::to_string(i)] = dir.make<TH1F>("BefPreS_Pt_Binned", "BefPreS_Pt_Binned", 50, 0, PtHistoUpperBound);
      Name.append(Suffix);
    tuple->BefPreS_TOF_Binned[std::to_string(i)] = dir.make<TH1F>("BefPreS_TOF_Binned", "BefPreS_TOF_Binned", 150, -1, 5);
  }

  tuple->PostS_CutIdVsEta_RegionA = dir.make<TH2F>("PostS_CutIdVsEta_RegionA", ";NCuts;#eta (RegionA)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->PostS_CutIdVsEta_RegionB = dir.make<TH2F>("PostS_CutIdVsEta_RegionB", ";NCuts;#eta (RegionB)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->PostS_CutIdVsEta_RegionC = dir.make<TH2F>("PostS_CutIdVsEta_RegionC", ";NCuts;#eta (RegionC)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->PostS_CutIdVsEta_RegionD = dir.make<TH2F>("PostS_CutIdVsEta_RegionD", ";NCuts;#eta (RegionD)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->PostS_CutIdVsEta_RegionE = dir.make<TH2F>("PostS_CutIdVsEta_RegionE", ";NCuts;#eta (RegionE)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->PostS_CutIdVsEta_RegionF = dir.make<TH2F>("PostS_CutIdVsEta_RegionF", ";NCuts;#eta (RegionF)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->PostS_CutIdVsEta_RegionG = dir.make<TH2F>("PostS_CutIdVsEta_RegionG", ";NCuts;#eta (RegionG)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  tuple->PostS_CutIdVsEta_RegionH = dir.make<TH2F>("PostS_CutIdVsEta_RegionH", ";NCuts;#eta (RegionH)", NCuts, 0, NCuts, 52, -2.6, 2.6);
  
  tuple->PostS_CutIdVsBeta_postPt = dir.make<TH2F>("PostS_CutIdVsBeta_postPt", ";NCuts;#beta (p_{T} > p_{T,cut})", NCuts, 0, NCuts, 20, 0, 1);
  tuple->PostS_CutIdVsBeta_postPtAndIas = dir.make<TH2F>("PostS_CutIdVsBeta_postPtAndIas", ";NCuts;#beta (p_{T} > p_{T,cut} and I_{as} > I_{as,cut} )", NCuts, 0, NCuts, 20, 0, 1);
  tuple->PostS_CutIdVsBeta_postPtAndIasAndTOF = dir.make<TH2F>("PostS_CutIdVsBeta_postPtAndIasAndTOF", ";NCuts;#beta (p_{T} > p_{T,cut} and I_{as} > I_{as,cut} and TOF > TOF_{cut} ", NCuts, 0, NCuts, 20, 0, 1);

  tuple->PostS_CutIdVsP = dir.make<TH2F>("PostS_CutIdVsP", ";NCuts;PostS_P", NCuts, 0, NCuts, 50, 0, PtHistoUpperBound);
  tuple->PostS_CutIdVsPt = dir.make<TH2F>("PostS_CutIdVsPt", ";NCuts;PostS_Pt", NCuts, 0, NCuts, 50, 0, PtHistoUpperBound);
  tuple->PostS_CutIdVsIas = dir.make<TH2F>("PostS_CutIdVsIas", ";NCuts;PostS_Ias", NCuts, 0, NCuts, 10, 0., 1.);
  tuple->PostS_CutIdVsIh = dir.make<TH2F>("PostS_CutIdVsIh", ";NCuts;PostS_Ih", NCuts, 0, NCuts, 100, 0, dEdxM_UpLim);
  tuple->PostS_CutIdVsTOF = dir.make<TH2F>("PostS_CutIdVsTOF", ";NCuts;PostS_TOF", NCuts, 0, NCuts, 50, 1, 5);
//tuple->PostS_CutIdVsEtaVsIas = dir.make<TH3F>("PostS_CutIdVsEtaVsIas", ";NCuts;PostS_EtaIs", NCuts, 0,  NCuts, 50,-3, 3, 10, 0., 1.);
//tuple->PostS_CutIdVsEtaVsIm = dir.make<TH3F>("PostS_CutIdVsEtaVsIm", ";NCuts;PostS_EtaIh", NCuts, 0,  NCuts, 50,-3, 3,100, 0, dEdxM_UpLim);
//tuple->PostS_CutIdVsEtaVsP  = dir.make<TH3F>("PostS_CutIdVsEtaVsP", ";NCuts;PostS_EtaP", NCuts, 0,  NCuts, 50,-3, 3, 50, 0, PtHistoUpperBound);
//tuple->PostS_CutIdVsEtaVsPt = dir.make<TH3F>("PostS_CutIdVsEtaVsPt", ";NCuts;PostS_EtaPt", NCuts, 0,  NCuts, 50,-3, 3, 50, 0, PtHistoUpperBound);
//tuple->PostS_CutIdVsEtaVsTOF = dir.make<TH3F>("PostS_CutIdVsEtaVsTOF", ";NCuts;PostS_EtaTOF", NCuts, 0,  NCuts, 50,-3, 3, 50, 0, 3);
  tuple->PostS_CutIdVsPVsIas = dir.make<TH3F>("PostS_CutIdVsPVsIas", ";NCuts;P;Ias", NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 10, 0., 1.);
  tuple->PostS_CutIdVsPVsIh = dir.make<TH3F>("PostS_CutIdVsPVsIh", ";NCuts;P;Ih", NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  tuple->PostS_CutIdVsPtVsIas = dir.make<TH3F>("PostS_CutIdVsPtVsIas", ";NCuts;Pt;Ias", NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 10, 0., 1.);
  tuple->PostS_CutIdVsPtVsIh = dir.make<TH3F>("PostS_CutIdVsPtVsIh", ";NCuts;Pt;Ih", NCuts, 0, NCuts, 50, 0, PtHistoUpperBound, 100, 0, dEdxM_UpLim);
  tuple->PostS_CutIdVsTOFVsIas = dir.make<TH3F>("PostS_CutIdVsTOFVsIas", ";NCuts;TOF;Ias", NCuts, 0, NCuts, 50, 0, 5, 10, 0., 1.);
  tuple->PostS_CutIdVsTOFVsIh = dir.make<TH3F>("PostS_CutIdVsTOFVsIh", ";NCuts;TOF;Ih", NCuts, 0, NCuts, 50, 0, 5, 100, 0, dEdxM_UpLim);

  tuple->H_D_DzSidebands = dir.make<TH2F>("H_D_DzSidebands", ";NCuts;H_D_DzSidebands", NCuts, 0, NCuts, DzRegions, 0, DzRegions);

  // Background prediction histograms don't need to be made for signal or individual MC samples
  // if (!isSignal) {
  // Although not needed let's still do it for every input
  if (true) {
    tuple->H_A = dir.make<TH1D>("H_A", ";NCuts;H_A", NCuts, 0, NCuts);
    tuple->H_B = dir.make<TH1D>("H_B", ";NCuts;H_B", NCuts, 0, NCuts);
    tuple->H_C = dir.make<TH1D>("H_C", ";NCuts;H_C", NCuts, 0, NCuts);
    tuple->H_D = dir.make<TH1D>("H_D", ";NCuts;H_D", NCuts, 0, NCuts);
    tuple->H_E = dir.make<TH1D>("H_E", ";NCuts;H_E", NCuts, 0, NCuts);
    tuple->H_F = dir.make<TH1D>("H_F", ";NCuts;H_F", NCuts, 0, NCuts);
    tuple->H_G = dir.make<TH1D>("H_G", ";NCuts;H_G", NCuts, 0, NCuts);
    tuple->H_H = dir.make<TH1D>("H_H", ";NCuts;H_H", NCuts, 0, NCuts);

    //Initialize histograms for number of bins.  For everything but muon only PredBins=0 so no histograms created
    for (int i = 0; i < PredBins; i++) {
      char Suffix[1024];
      sprintf(Suffix, "_%i", i);
      Name = "H_B_Binned";
      Name.append(Suffix);
      tuple->H_B_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;H_B_Binned", NCuts, 0, NCuts);
      Name = "H_D_Binned";
      Name.append(Suffix);
      tuple->H_D_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;H_D_Binned", NCuts, 0, NCuts);
      Name = "H_F_Binned";
      Name.append(Suffix);
      tuple->H_F_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;H_F_Binned", NCuts, 0, NCuts);
      Name = "H_H_Binned";
      Name.append(Suffix);
      tuple->H_H_Binned[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;H_H_Binned", NCuts, 0, NCuts);
    }

    // Where are these used?
      tuple->Hist_Is = dir.make<TH1D>("Hist_Is", "Hist_Is", 200, 0, dEdxS_UpLim);
      tuple->Hist_Pt = dir.make<TH1D>("Hist_Pt", "Hist_Pt", 200, 0, PtHistoUpperBound);
      tuple->Hist_TOF = dir.make<TH1D>("Hist_TOF", "Hist_TOF", 200, -10, 20);
    //The following are only used to create the predicted mass spectrum.  Memory intensive so don't initialize for analyses not doing mass fits
    if (TypeMode < 3) {
      tuple->Pred_I = dir.make<TH2F>("Pred_I", ";NCuts;Pred_I", NCuts, 0, NCuts, 400, 0, dEdxM_UpLim);
      tuple->Pred_EtaI = dir.make<TH3F>("Pred_EtaI", ";NCuts;Pred_EtaI", NCuts, 0, NCuts, EtaBins, -3., 3., 400, 0, dEdxM_UpLim);
      tuple->Pred_EtaB = dir.make<TH2F>("Pred_EtaB", ";NCuts;Pred_EtaB", NCuts, 0, NCuts, EtaBins, -3., 3.);
      tuple->Pred_EtaS = dir.make<TH2F>("Pred_EtaS", ";NCuts;Pred_EtaS", NCuts, 0, NCuts, EtaBins, -3., 3.);
      tuple->Pred_EtaS2 = dir.make<TH2F>("Pred_EtaS2", ";NCuts;Pred_EtaS2", NCuts, 0, NCuts, EtaBins, -3., 3.);
      tuple->Pred_EtaP = dir.make<TH3F>("Pred_EtaP", ";NCuts;#eta;P", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);
      tuple->Pred_TOF = dir.make<TH2F>("Pred_TOF", ";NCuts;Pred_TOF", NCuts, 0, NCuts, 200, GlobalMinTOF, 5);
      //pz

      tuple->PDF_G_EtaP = dir.make<TH3F>("PDF_G_EtaP", ";NCuts;PDF_G_Eta;P", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);
      tuple->PDF_C_EtaP = dir.make<TH3F>("PDF_C_EtaP", ";NCuts;PDF_C_Eta;P", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);

      tuple->PDF_A_Eta = dir.make<TH2F>("PDF_A_Eta", ";NCuts;PDF_A_Eta", NCuts, 0, NCuts, EtaBins, -3., 3.);
      tuple->PDF_E_Eta = dir.make<TH2F>("PDF_E_Eta", ";NCuts;PDF_E_Eta", NCuts, 0, NCuts, EtaBins, -3., 3.);

      tuple->PDF_B_EtaICK = dir.make<TH3F>("PDF_B_EtaICK", ";NCuts;PDF_B_EtaICK", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);
      tuple->PDF_F_EtaICK = dir.make<TH3F>("PDF_F_EtaICK", ";NCuts;PDF_F_EtaICK", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);

      tuple->PDF_H_EtaMass = dir.make<TH3F>("PDF_H_EtaMass", ";NCuts;PDF_H_Eta;Mass", NCuts, 0, NCuts, EtaBins, -3., 3., MassNBins, 0, MassHistoUpperBound);

      //pz FLIP
      tuple->PDF_G_EtaP_Flip = dir.make<TH3F>("PDF_G_EtaP_Flip", ";NCuts;PDF_G_EtaP_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);
      tuple->PDF_C_EtaP_Flip = dir.make<TH3F>("PDF_C_EtaP_Flip", ";NCuts;PDF_C_EtaP_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);

      tuple->PDF_A_Eta_Flip = dir.make<TH2F>("PDF_A_Eta_Flip", ";NCuts;PDF_A_Eta_Flip", NCuts, 0, NCuts, EtaBins, -3., 3.);
      tuple->PDF_E_Eta_Flip = dir.make<TH2F>("PDF_E_Eta_Flip", ";NCuts;PDF_E_Eta_Flip", NCuts, 0, NCuts, EtaBins, -3., 3.);

      tuple->PDF_B_EtaICK_Flip = dir.make<TH3F>("PDF_B_EtaICK_Flip", ";NCuts;PDF_B_EtaICK_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);
      tuple->PDF_F_EtaICK_Flip = dir.make<TH3F>("PDF_F_EtaICK_Flip", ";NCuts;PDF_F_EtaICK_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);

      tuple->PDF_H_EtaMass_Flip = dir.make<TH3F>("PDF_H_EtaMass_Flip", ";NCuts;PDF_H_EtaMass_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., MassNBins, 0, MassHistoUpperBound);
    }

      tuple->RegionD_I = dir.make<TH2F>("RegionD_I", ";NCuts;RegionD_I", NCuts, 0, NCuts, 400, 0, dEdxM_UpLim);
      tuple->RegionD_Ias = dir.make<TH2F>("RegionD_Ias", ";NCuts;RegionD_Ias", NCuts, 0, NCuts, 100, 0, dEdxS_UpLim);
      tuple->RegionD_P = dir.make<TH2F>("RegionD_P", ";NCuts;RegionD_P", NCuts, 0, NCuts, 200, GlobalMinPt, PtHistoUpperBound);
      tuple->RegionD_TOF = dir.make<TH2F>("RegionD_TOF", ";NCuts;RegionD_TOF", NCuts, 0, NCuts, 200, GlobalMinTOF, 5);

      tuple->RegionH_Ias = dir.make<TH2F>("RegionH_Ias", ";NCuts;RegionH_Ias", NCuts, 0, NCuts, 100, 0, dEdxS_UpLim);

      tuple->H_A_Flip = dir.make<TH1D>("H_A_Flip", ";NCuts_Flip;H_A_Flip", NCuts_Flip, 0, NCuts_Flip);
      tuple->H_B_Flip = dir.make<TH1D>("H_B_Flip", ";NCuts_Flip;H_B_Flip", NCuts_Flip, 0, NCuts_Flip);
      tuple->H_C_Flip = dir.make<TH1D>("H_C_Flip", ";NCuts_Flip;H_C_Flip", NCuts_Flip, 0, NCuts_Flip);
      tuple->H_D_Flip = dir.make<TH1D>("H_D_Flip", ";NCuts_Flip;H_D_Flip", NCuts_Flip, 0, NCuts_Flip);
      tuple->H_E_Flip = dir.make<TH1D>("H_E_Flip", ";NCuts_Flip;H_E_Flip", NCuts_Flip, 0, NCuts_Flip);
      tuple->H_F_Flip = dir.make<TH1D>("H_F_Flip", ";NCuts_Flip;H_F_Flip", NCuts_Flip, 0, NCuts_Flip);
      tuple->H_G_Flip = dir.make<TH1D>("H_G_Flip", ";NCuts_Flip;H_G_Flip", NCuts_Flip, 0, NCuts_Flip);
      tuple->H_H_Flip = dir.make<TH1D>("H_H_Flip", ";NCuts_Flip;H_H_Flip", NCuts_Flip, 0, NCuts_Flip);

    for (int i = 0; i < PredBins; i++) {
      char Suffix[1024];
      sprintf(Suffix, "_%i", i);
      Name = "H_B_Binned_Flip";
      Name.append(Suffix);
      tuple->H_B_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;", NCuts, 0, NCuts);
      Name = "H_D_Binned_Flip";
      Name.append(Suffix);
      tuple->H_D_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;", NCuts, 0, NCuts);
      Name = "H_F_Binned_Flip";
      Name.append(Suffix);
      tuple->H_F_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;", NCuts, 0, NCuts);
      Name = "H_H_Binned_Flip";
      Name.append(Suffix);
      tuple->H_H_Binned_Flip[std::to_string(i)] = dir.make<TH1D>(Name.c_str(), ";NCuts;", NCuts, 0, NCuts);
    }

    //The following are only used to create the predicted mass spectrum.  Memory intensive so don't initialize for analyses not doing mass fits
    if (TypeMode < 3) {

      tuple->Pred_I_Flip = dir.make<TH2F>("Pred_I_Flip", ";NCuts_Flip;Pred_I_Flip", NCuts_Flip, 0, NCuts_Flip, 400, 0, dEdxM_UpLim);

      tuple->Pred_EtaB_Flip = dir.make<TH2F>("Pred_EtaB_Flip", ";NCuts_Flip;Pred_EtaB_Flip", NCuts_Flip, 0, NCuts_Flip, EtaBins, -3., 3.);

      tuple->Pred_EtaS_Flip = dir.make<TH2F>("Pred_EtaS_Flip", ";NCuts_Flip;Pred_EtaS_Flip", NCuts_Flip, 0, NCuts_Flip, EtaBins, -3., 3.);

      tuple->Pred_EtaS2_Flip = dir.make<TH2F>("Pred_EtaS2_Flip", ";NCuts_Flip;Pred_EtaS2_Flip", NCuts_Flip, 0, NCuts_Flip, EtaBins, -3., 3.);

      tuple->Pred_EtaP_Flip = dir.make<TH3F>("Pred_EtaP_Flip", ";NCuts_Flip;Pred_EtaP_Flip", NCuts_Flip, 0, NCuts_Flip, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);

      tuple->Pred_TOF_Flip = dir.make<TH2F>("Pred_TOF_Flip", ";NCuts_Flip;Pred_TOF_Flip", NCuts_Flip, 0, NCuts_Flip, 200, GlobalMinTOF, 5);
    }


    tuple->RegionD_I_Flip = dir.make<TH2F>("RegionD_I_Flip", ";NCuts_Flip;RegionD_I_Flip", NCuts_Flip, 0, NCuts_Flip, 400, 0, dEdxM_UpLim);

    tuple->RegionD_Ias_Flip =
        dir.make<TH2F>("RegionD_Ias_Flip", ";NCuts_Flip;RegionD_Ias_Flip", NCuts_Flip, 0, NCuts_Flip, 100, 0, dEdxS_UpLim);

    tuple->RegionD_P_Flip =
        dir.make<TH2F>("RegionD_P_Flip", ";NCuts_Flip;RegionD_P_Flip", NCuts_Flip, 0, NCuts_Flip, 200, GlobalMinPt, PtHistoUpperBound);

    tuple->RegionD_TOF_Flip = dir.make<TH2F>("RegionD_TOF_Flip", ";NCuts_Flip;RegionD_TOF_Flip", NCuts_Flip, 0, NCuts_Flip, 200, -3, 1);


    tuple->RegionH_Ias_Flip =
        dir.make<TH2F>("RegionH_Ias_Flip", ";NCuts_Flip;RegionH_Ias_Flip", NCuts_Flip, 0, NCuts_Flip, 100, 0, dEdxS_UpLim);


    tuple->CtrlPt_S1_Is = dir.make<TH1D>("CtrlPt_S1_Is", "CtrlPt_S1_Is", 200, 0, dEdxS_UpLim);

    tuple->CtrlPt_S2_Is = dir.make<TH1D>("CtrlPt_S2_Is", "CtrlPt_S2_Is", 200, 0, dEdxS_UpLim);

    tuple->CtrlPt_S3_Is = dir.make<TH1D>("CtrlPt_S3_Is", "CtrlPt_S3_Is", 200, 0, dEdxS_UpLim);

    tuple->CtrlPt_S4_Is = dir.make<TH1D>("CtrlPt_S4_Is", "CtrlPt_S4_Is", 200, 0, dEdxS_UpLim);


    tuple->CtrlPt_S1_Ih = dir.make<TH1D>("CtrlPt_S1_Ih", "CtrlPt_S1_Ih", 400, 0, dEdxM_UpLim);

    tuple->CtrlPt_S2_Ih = dir.make<TH1D>("CtrlPt_S2_Ih", "CtrlPt_S2_Ih", 400, 0, dEdxM_UpLim);

    tuple->CtrlPt_S3_Ih = dir.make<TH1D>("CtrlPt_S3_Ih", "CtrlPt_S3_Ih", 400, 0, dEdxM_UpLim);

    tuple->CtrlPt_S4_Ih = dir.make<TH1D>("CtrlPt_S4_Ih", "CtrlPt_S4_Ih", 400, 0, dEdxM_UpLim);


    tuple->CtrlIs_S1_TOF = dir.make<TH1D>("CtrlIs_S1_TOF", "CtrlIs_S1_TOF", 200, 0, 5);

    tuple->CtrlIs_S2_TOF = dir.make<TH1D>("CtrlIs_S2_TOF", "CtrlIs_S2_TOF", 200, 0, 5);

    tuple->CtrlIs_S3_TOF = dir.make<TH1D>("CtrlIs_S3_TOF", "CtrlIs_S3_TOF", 200, 0, 5);

    tuple->CtrlIs_S4_TOF = dir.make<TH1D>("CtrlIs_S4_TOF", "CtrlIs_S4_TOF", 200, 0, 5);


    tuple->CtrlIh_S1_TOF = dir.make<TH1D>("CtrlIh_S1_TOF", "CtrlIh_S1_TOF", 200, 0, 5);

    tuple->CtrlIh_S2_TOF = dir.make<TH1D>("CtrlIh_S2_TOF", "CtrlIh_S2_TOF", 200, 0, 5);

    tuple->CtrlIh_S3_TOF = dir.make<TH1D>("CtrlIh_S3_TOF", "CtrlIh_S3_TOF", 200, 0, 5);

    tuple->CtrlIh_S4_TOF = dir.make<TH1D>("CtrlIh_S4_TOF", "CtrlIh_S4_TOF", 200, 0, 5);


    tuple->CtrlPt_S1_TOF = dir.make<TH1D>("CtrlPt_S1_TOF", "CtrlPt_S1_TOF", 200, -2, 7);

    tuple->CtrlPt_S2_TOF = dir.make<TH1D>("CtrlPt_S2_TOF", "CtrlPt_S2_TOF", 200, -2, 7);

    tuple->CtrlPt_S3_TOF = dir.make<TH1D>("CtrlPt_S3_TOF", "CtrlPt_S3_TOF", 200, -2, 7);

    tuple->CtrlPt_S4_TOF = dir.make<TH1D>("CtrlPt_S4_TOF", "CtrlPt_S4_TOF", 200, -2, 7);

    for (int i = 0; i < PredBins; i++) {
      char Suffix[1024];
      sprintf(Suffix, "_%i", i);

      Name.append(Suffix);
      tuple->CtrlPt_S1_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S1_TOF_Binned", "CtrlPt_S1_TOF_Binned", 200, -2, 7);

      Name.append(Suffix);
      tuple->CtrlPt_S2_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S2_TOF_Binned", "CtrlPt_S2_TOF_Binned", 200, -2, 7);

      Name.append(Suffix);
      tuple->CtrlPt_S3_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S3_TOF_Binned", "CtrlPt_S3_TOF_Binned", 200, -2, 7);

      Name.append(Suffix);
      tuple->CtrlPt_S4_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S4_TOF_Binned", "CtrlPt_S4_TOF_Binned", 200, -2, 7);
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
    tuple->Tree->Branch("HLT_PFMET120_PFMHT120_IDTight",
                        &tuple->Tree_HLT_PFMET120_PFMHT120_IDTight,
                        "HLT_PFMET120_PFMHT120_IDTight/O");
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
      tuple->PostS_CutIdVsEta_RegionD->Fill(CutIndex, track->eta());
    } else if (PassTOFCut && PassPtCut && !PassICut) {  //Region C
      tuple->H_C->Fill(CutIndex, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
      tuple->PDF_C_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz
      //Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
      tuple->PostS_CutIdVsEta_RegionC->Fill(CutIndex, track->eta());
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
      tuple->PostS_CutIdVsEta_RegionB->Fill(CutIndex, track->eta());
    } else if (PassTOFCut && !PassPtCut && !PassICut) {  //Region A
      tuple->H_A->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_TOF->Fill(CutIndex, MuonTOF, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaB->Fill(CutIndex, track->eta(), Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS2->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->PostS_CutIdVsEta_RegionA->Fill(CutIndex, track->eta());
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
        tuple->PostS_CutIdVsEta_RegionH->Fill(CutIndex, track->eta());
    } else if (!PassTOFCut && PassPtCut && !PassICut) {  //Region G
      tuple->H_G->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
      tuple->PostS_CutIdVsEta_RegionG->Fill(CutIndex, track->eta());
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
      tuple->PostS_CutIdVsEta_RegionF->Fill(CutIndex, track->eta());
      if (TypeMode == 2)
        tuple->PDF_F_EtaICK->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz

    } else if (!PassTOFCut && !PassPtCut && !PassICut) {  //Region E
      tuple->H_E->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaB->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->PostS_CutIdVsEta_RegionE->Fill(CutIndex, track->eta());
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
