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
                       bool calcSyst_,
                       int TypeMode,
                       bool isSignal,
                       bool doBefTrigPlots_,
                       bool doBefPreSplots_,
                       bool doPostPreSplots_,
                       bool createGiTemplates_,
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
                       float GlobalMinTOF,
                       bool tapeRecallOnly_);

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
                        const std::vector<int>  &BunchXing,
                        const std::vector<int>  &nPU,
                        const std::vector<float>  &nPUmean,
                        const unsigned int &nofVertices,
                        const int &npv,
                        const std::vector<float>  &pvX,
                        const std::vector<float>  &pvY,
                        const std::vector<float>  &pvZ,
                        const std::vector<float>  &pvRho,
                        const std::vector<int>  &pvNdof,
                        const std::vector<float>  &pvChi2,
                        const std::vector<float>  &pvSumPt2,

                        const unsigned int &Hscp,
                        const unsigned int &nMuons,
                        const unsigned int &njets,
                        const float &weight,
                        const float &generator_weight,
                        const float &generator_binning_values,
                        const std::vector<bool> &triggerDecision,
                        const std::vector<int> &triggerHLTPrescale,
                        const std::vector<std::vector<float>> &triggerObjectE,
                        const std::vector<std::vector<float>> &triggerObjectPt,
                        const std::vector<std::vector<float>> &triggerObjectEta,
                        const std::vector<std::vector<float>> &triggerObjectPhi,
                        const bool &HLT_Mu50,
                        const bool &HLT_PFMET120_PFMHT120_IDTight,
                        const bool &HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                        const bool &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                        const bool &HLT_MET105_IsoTrk50,
                        const float &RecoCaloMET,
                        const float &RecoCaloMET_phi,
                        const float &RecoCaloMET_sigf,
                        const float &RecoPFMET,
                        const float &RecoPFMET_phi,
                        const float &RecoPFMET_sigf,
                        const float &RecoPFMHT,
                        const float &HLTCaloMET,
                        const float &HLTCaloMET_phi,
                        const float &HLTCaloMET_sigf,
                        const float &HLTCaloMETClean,
                        const float &HLTCaloMETClean_phi,
                        const float &HLTCaloMETClean_sigf,
                        const float &HLTCaloMHT,
                        const float &HLTCaloMHT_phi,
                        const float &HLTCaloMHT_sigf,
                        const float &HLTPFMET,
                        const float &HLTPFMET_phi,
                        const float &HLTPFMET_sigf,
                        const float &HLTPFMHT,
                        const float &HLTPFMHT_phi,
                        const float &HLTPFMHT_sigf,
                        const bool &matchedMuonWasFound,
                        const std::vector<int> &gParticleId,
                        const std::vector<int> &gParticleStatus,
                        const std::vector<float> &gParticleE,
                        const std::vector<float> &gParticlePt,
                        const std::vector<float> &gParticlePz,
                        const std::vector<float> &gParticleEta,
                        const std::vector<float> &gParticlePhi,
                        const std::vector<float> &gParticleBeta,
                        const std::vector<int> &gParticleCharge,
                        const std::vector<float> &gParticleProdVertexX,
                        const std::vector<float> &gParticleProdVertexY,
                        const std::vector<float> &gParticleProdVertexZ,
                        const std::vector<int> &gParticleMotherId,
                        const std::vector<int> &gParticleMotherIndex,
                        const std::vector<float> &eleE,
                        const std::vector<float> &elePt,
                        const std::vector<float> &eleEta,
                        const std::vector<float> &elePhi,
                        const std::vector<float> &eleCharge,
                        const std::vector<float> &eleE_SC,
                        const std::vector<float> &eleEta_SC,
                        const std::vector<float> &elePhi_SC,
                        const std::vector<float> &eleSigmaIetaIeta,
                        const std::vector<float> &eleFull5x5SigmaIetaIeta,
                        const std::vector<float> &eleR9,
                        const std::vector<float> &ele_dEta,
                        const std::vector<float> &ele_dPhi,
                        const std::vector<float> &ele_HoverE,
                        const std::vector<float> &ele_d0,
                        const std::vector<float> &ele_dZ,
                        const std::vector<float> &ele_pileupIso,
                        const std::vector<float> &ele_chargedIso,
                        const std::vector<float> &ele_photonIso,
                        const std::vector<float> &ele_neutralHadIso,
                        const std::vector<int> &ele_MissHits,
                        const std::vector<bool> &ele_passCutBasedIDVeto,
                        const std::vector<bool> &ele_passCutBasedIDLoose,
                        const std::vector<bool> &ele_passCutBasedIDMedium,
                        const std::vector<bool> &ele_passCutBasedIDTight,
                        const std::vector<bool> &ele_passMVAIsoIDWP80,
                        const std::vector<bool> &ele_passMVAIsoIDWP90,
                        const std::vector<bool> &ele_passMVAIsoIDWPHZZ,
                        const std::vector<bool> &ele_passMVAIsoIDWPLoose,
                        const std::vector<bool> &ele_passMVANoIsoIDWP80,
                        const std::vector<bool> &ele_passMVANoIsoIDWP90,
                        const std::vector<bool> &ele_passMVANoIsoIDWPLoose,
                        const std::vector<bool> &ele_PassConvVeto,
                        const std::vector<float> &ele_OneOverEminusOneOverP,
                        const std::vector<float> &muonE,
                        const std::vector<float> &muonPt,
                        const std::vector<float> &muonEta,
                        const std::vector<float> &muonPhi,
                        const std::vector<int> &muonCharge,
                        const std::vector<bool> &muonIsLoose,
                        const std::vector<bool> &muonIsMedium,
                        const std::vector<bool> &muonIsTight,
                        const std::vector<float> &muon_d0,
                        const std::vector<float> &muon_d0Err,
                        const std::vector<float> &muon_dZ,
                        const std::vector<float> &muon_ip3d,
                        const std::vector<float> &muon_ip3dSignificance,
                        const std::vector<unsigned int> &muonType,
                        const std::vector<unsigned int> &muonQuality,
                        const std::vector<float> &muon_pileupIso,
                        const std::vector<float> &muon_chargedIso,
                        const std::vector<float> &muon_photonIso,
                        const std::vector<float> &muon_neutralHadIso,
                        const std::vector<float> &muon_validFractionTrackerHits,
                        const std::vector<float> &muon_normChi2,
                        const std::vector<float> &muon_chi2LocalPosition,
                        const std::vector<float> &muon_kinkFinder,
                        const std::vector<float> &muon_segmentCompatability,
                        const std::vector<float> &muon_trkIso,
                        const std::vector<float> &muon_tuneP_Pt,
                        const std::vector<float> &muon_tuneP_PtErr,
                        const std::vector<float> &muon_tuneP_Eta,
                        const std::vector<float> &muon_tuneP_Phi,
                        const std::vector<int> &muon_tuneP_MuonBestTrackType,
                        const std::vector<bool> &muon_isHighPtMuon,
                        const std::vector<bool> &muon_isTrackerHighPtMuon,
                        const std::vector<float> &vect_jet_pt,
                        const std::vector<float> &vect_jet_eta,
                        const std::vector<float> &vect_jet_phi,
                        const std::vector<float> &vect_jet_mass,
                        const std::vector<float> &vect_jet_energy,
                        const std::vector<float> &vect_jet_pdgId,
                        const std::vector<float> &vect_jet_et,
                        const std::vector<float> &vect_jet_chargedEmEnergyFraction,
                        const std::vector<float> &vect_jet_neutralEmEnergyFraction,
                        const std::vector<float> &vect_jet_chargedHadronEnergyFraction,
                        const std::vector<float> &vect_jet_neutralHadronEnergyFraction,
                        const std::vector<float> &vect_jet_muonEnergyFraction,
                        const std::vector<int> &vect_jet_chargedMultiplicity,
                        const std::vector<int> &vect_jet_neutralMultiplicity,
                        const std::vector<float> &vect_jet_jetArea,
                        const std::vector<float> &vect_jet_pileupE,
                        const std::vector<float> &vect_mT,
                        const std::vector<bool> &passCutPt55,
                        const std::vector<bool> &passPreselection,
                        const std::vector<bool> &passPreselectionSept8,
                        const std::vector<bool> &passSelection,
                        const std::vector<bool> &isPFMuon,
                        const std::vector<bool> &PFMuonPt,
                        const std::vector<float> &Charge,
                        const std::vector<float> &Pt,
                        const std::vector<float> &PtErr,
                        const std::vector<float> &Is_StripOnly,
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
                        const std::vector<bool>  &isPhoton,
                        const std::vector<bool>  &isElectron,
                        const std::vector<bool>  &isChHadron,
                        const std::vector<bool>  &isNeutHadron,
                        const std::vector<bool>  &isPfTrack,
                        const std::vector<bool>  &isUndefined,
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
                        const std::vector<float> &p,
                        const std::vector<float> &eta,
                        const std::vector<float> &phi,
                        const std::vector<unsigned int> &noh,
                        const std::vector<unsigned int> &noph,
                        const std::vector<float> &fovh,
                        const std::vector<unsigned int> &nomh,
                        const std::vector<float> &fovhd,
                        const std::vector<unsigned int> &nom,
                        const std::vector<float> &matchTrigMuon_minDeltaR,
                        const std::vector<float> &matchTrigMuon_pT,
                        const std::vector<float> &iso_TK,
                        const std::vector<float> &iso_ECAL,
                        const std::vector<float> &iso_HCAL,
                        const std::vector<float> &track_genTrackMiniIsoSumPt,
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
                        const std::vector<float> &genphi,
                        const std::vector<float> &HSCP_tuneP_Pt,
                        const std::vector<float> &HSCP_tuneP_PtErr,
                        const std::vector<float> &HSCP_tuneP_Eta,
                        const std::vector<float> &HSCP_tuneP_Phi,
                        const std::vector<int> &HSCP_tuneP_MuonBestTrackType,
                        const std::vector<int> &HSCP_ErrorHisto_bin,
                        const std::vector<int> &HSCP_type);


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
                   float Ias_quantiles[5],
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
                                 bool calcSyst_,
                                 int TypeMode,
                                 bool isSignal,
                                 bool doBefTrigPlots_,
                                 bool doBefPreSplots_,
                                 bool doPostPreSplots_,
                                 bool createGiTemplates_,
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
                                 float GlobalMinTOF,
                                 bool tapeRecallOnly_) {
  std::string Name;

  TH1::SetDefaultSumw2(kTRUE);
  // if the only purpose is to trick CRAB to do a TAPERECALL
  if (tapeRecallOnly_) return;

  tuple->IntLumi = dir.make<TProfile>("IntLumi", ";IntLumi", 1, 0, 1);
  tuple->XSection = dir.make<TProfile>("XSection", ";XSection", 1, 0, 1);

  tuple->NumEvents = dir.make<TH1F>("NumEvents", ";;Number of events / category", 4, -0.5, 3.5);
  tuple->NumEvents->GetXaxis()->SetBinLabel(1,"All events");
  tuple->NumEvents->GetXaxis()->SetBinLabel(2,"Events w/ PU syst");
  tuple->NumEvents->GetXaxis()->SetBinLabel(3,"After trigger");
  tuple->NumEvents->GetXaxis()->SetBinLabel(4,"After HLT obj to evt matching");

  tuple->dRMinHLTMuon = dir.make<TH1F>("dRMinHLTMuon", ";#Delta R_{min,mu,HLT};Events / 0.04",100,0.,4.);
  tuple->dRMinHLTMuon_lowDeltaR = dir.make<TH1F>("dRMinHLTMuon_lowDeltaR", ";#Delta R_{min,mu,HLT};Events / 0.01",40,0.,0.4);

  tuple->ErrorHisto = dir.make<TH1F>("ErrorHisto", ";;", 11, -0.5, 10.5);
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(1,"All tracks");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(2,"Not tracker / global muon");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(3,"Track is null");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(4,"No PV associated");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(5,"No gen match found");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(6,"Gen match found, but too far");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(7,"No dEdx associated");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(8,"Not a collision track");
  tuple->ErrorHisto->GetXaxis()->SetBinLabel(9,"Has status 91 around it");
    
  tuple->Gen_HSCPCandidateType = dir.make<TH1F>("Gen_HSCPCandidateType", ";;Number of generator candidate / category", 6, -0.5, 5.5);
  tuple->Gen_HSCPCandidateType->GetXaxis()->SetBinLabel(1,"All candidates");
  tuple->Gen_HSCPCandidateType->GetXaxis()->SetBinLabel(2,"Neutral HSCP");
  tuple->Gen_HSCPCandidateType->GetXaxis()->SetBinLabel(3,"Single-charged");
  tuple->Gen_HSCPCandidateType->GetXaxis()->SetBinLabel(4,"Double-charged R-hadrons");
  tuple->Gen_HSCPCandidateType->GetXaxis()->SetBinLabel(5,"Tau-prime (1e or 2e)");
  tuple->Gen_HSCPCandidateType->GetXaxis()->SetBinLabel(6,"Else");
  
  tuple->CutFlow = dir.make<TH1F>("CutFlow", ";;Tracks / category", 21, -0.5, 20.5);
  tuple->CutFlow->GetXaxis()->SetBinLabel(1,"All tracks");
  tuple->CutFlow->GetXaxis()->SetBinLabel(2,"Technical");
  tuple->CutFlow->GetXaxis()->SetBinLabel(3,"Trigger");
  tuple->CutFlow->GetXaxis()->SetBinLabel(4,"p_{T}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(5,"#eta");
  tuple->CutFlow->GetXaxis()->SetBinLabel(6,"N_{no-L1 pixel hits}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(7,"f_{valid/all hits}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(8,"N_{dEdx hits}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(9,"HighPurity");
  tuple->CutFlow->GetXaxis()->SetBinLabel(10,"#chi^{2} / N_{dof}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(11,"d_{z}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(12,"d_{xy}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(13,"MiniRelIsoAll");
  tuple->CutFlow->GetXaxis()->SetBinLabel(14,"MiniRelTkIso");
  tuple->CutFlow->GetXaxis()->SetBinLabel(15,"E/p");
  tuple->CutFlow->GetXaxis()->SetBinLabel(16,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(17,"F_{i}");
  tuple->CutFlow->GetXaxis()->SetBinLabel(18,"SR0");
  tuple->CutFlow->GetXaxis()->SetBinLabel(19,"SR1");
  tuple->CutFlow->GetXaxis()->SetBinLabel(20,"SR2");
  tuple->CutFlow->GetXaxis()->SetBinLabel(21,"SR3");
  
  tuple->EventCutFlow = dir.make<TH1F>("EventCutFlow", ";;Events / category", 21, -0.5, 20.5);
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(1,"All tracks");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(2,"Technical");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(3,"Trigger");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(4,"p_{T}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(5,"#eta");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(6,"N_{no-L1 pixel hits}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(7,"f_{valid/all hits}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(8,"N_{dEdx hits}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(9,"HighPurity");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(10,"#chi^{2} / N_{dof}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(11,"d_{z}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(12,"d_{xy}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(13,"MiniRelIsoAll");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(14,"MiniRelTkIso");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(15,"E/p");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(16,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(17,"F_{i}");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(18,"SR0");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(19,"SR1");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(20,"SR2");
  tuple->EventCutFlow->GetXaxis()->SetBinLabel(21,"SR3");
  
  tuple->CutFlowReverse = dir.make<TH1F>("CutFlowReverse", ";;Tracks / category", 17, -0.5, 16.5);
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(1,"F_{i}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(2,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(3,"E/p");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(4,"MiniRelTkIso");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(5,"MiniRelIsoAll");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(6,"d_{xy}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(7,"d_{z}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(8,"#chi^{2} / N_{dof}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(9,"HighPurity");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(10,"N_{dEdx hits}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(11,"f_{valid/all hits}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(12,"N_{no-L1 pixel hits}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(13,"#eta");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(14,"p_{T}");
  tuple->CutFlowReverse->GetXaxis()->SetBinLabel(15,"Trigger");
  
  tuple->CutFlowProbQ =  dir.make<TH2F>("CutFlowProbQ", ";F_{i}^{pixels};",10, 0., 1.,17, -0.5, 16.5);
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(5,"f_{valid/all hits}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(8,"#chi^{2} / N_{dof}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(13,"E/p");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(14,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->CutFlowProbQ->GetYaxis()->SetBinLabel(15,"F_{i}");

  tuple->CutFlowEta = dir.make<TH2F>("CutFlowEta", ";#eta;", 50, -2.6, 2.6, 17, -0.5, 16.5);
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(5,"f_{valid/all hits}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(8,"#chi^{2} / N_{dof}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(13,"E/p");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(14,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->CutFlowEta->GetYaxis()->SetBinLabel(15,"F_{i}");
  
  tuple->CutFlowEoP = dir.make<TH2F>("CutFlowEoP", ";PF calo energy / momentum;", 25, 0, 1.5, 17, -0.5, 16.5);
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(5,"f_{valid/all hits}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(8,"#chi^{2} / N_{dof}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(13,"E/p");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(14,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->CutFlowEoP->GetYaxis()->SetBinLabel(15,"F_{i}");

  tuple->CutFlowPfType = dir.make<TH2F>("CutFlowPfType", ";;", 9, -0.5, 8.5, 17, -0.5, 16.5);
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(5,"f_{valid/all hits}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(8,"#chi^{2} / N_{dof}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(13,"E/p");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(14,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->CutFlowPfType->GetYaxis()->SetBinLabel(15,"F_{i}");
  // now the x-axis titles
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(1,"AllTracks");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(2,"PFtracks");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(3,"isElectron");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(4,"isMuon");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(5,"isPhoton");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(6,"isChHadron");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(7,"isNeutHadron");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(8,"isUndefined");
  tuple->CutFlowPfType->GetXaxis()->SetBinLabel(9,"notPFtrack");
  
  tuple->N1_FlowEta = dir.make<TH2F>("N1_FlowEta", ";#eta;", 50, -2.6, 2.6, 17, -0.5, 16.5);
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(1,"Trigger");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(2,"p_{T}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(3,"#eta");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(4,"N_{no-L1 pixel hits}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(5,"f_{valid/all hits}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(6,"N_{dEdx hits}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(7,"HighPurity");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(8,"#chi^{2} / N_{dof}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(9,"d_{z}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(10,"d_{xy}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(11,"MiniRelIsoAll");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(12,"MiniRelTkIso");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(13,"E/p");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(14,"#sigma_{p_{T}} / p_{T}^{2}");
  tuple->N1_FlowEta->GetYaxis()->SetBinLabel(15,"F_{i}");

  tuple->N1_Eta = dir.make<TH1F>("N1_Eta", ";#eta;Tracks / 0.05", 52, -2.6, 2.6);
  tuple->N1_Pt = dir.make<TH1F>("N1_Pt", ";p_{T} (GeV);Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
  tuple->N1_Pt_lowPt = dir.make<TH1F>("N1_Pt_lowPt", ";p_{T} (GeV);Tracks / 10 GeV", 50, 0., 500);
  tuple->N1_Chi2oNdof = dir.make<TH1F>("N1_Chi2oNdof", ";#chi^{2} / N_{dof};Tracks / 1", 20, 0, 20);

  tuple->N1_Qual = dir.make<TH1F>("N1_Qual", ";;Tracks / category", 2, -0.5, 1.5);
  tuple->N1_Qual->GetXaxis()->SetBinLabel(1,"Not-HighPurity");
  tuple->N1_Qual->GetXaxis()->SetBinLabel(2,"HighPurity");
  tuple->N1_TNOM = dir.make<TH1F>("N1_TNOM", ";Number of measurement;Tracks / 1", 40, -0.5, 39.5);
  tuple->N1_TNOPH = dir.make<TH1F>("N1_TNOPH", ";Number of pixel hits;Tracks / 1", 8, -0.5, 7.5);
  tuple->N1_TNOHFraction = dir.make<TH1F>("N1_TNOHFraction", ";Number of valid hit fraction;Tracks / 0.02", 50, 0., 1);
//  tuple->N1_nDof = dir.make<TH1F>("nDof", ";N_{dof}", 40, -0.5, 39.5);
//  tuple->N1_tofError = dir.make<TH1F>("tofError", ";tofError", 25, 0, 0.25);
  tuple->N1_TIsol = dir.make<TH1F>("TIsol", ";#Sigma_{R<0.3} p_{T} - p_{T,cand} (GeV) / 4 GeV", 25, 0, 100);
  tuple->N1_EoP = dir.make<TH1F>("N1_EoP", ";PF calo energy / momentum; Tracks / 0.06", 25, 0, 1.5);
  tuple->N1_ECalEoP = dir.make<TH1F>("N1_ECalEoP", ";PF ECal energy / momentum; Tracks / 0.06", 25, 0, 1.5);
  tuple->N1_HCalEoP = dir.make<TH1F>("N1_HCalEoP", ";PF HCal energy / momentum; Tracks / 0.06", 25, 0, 1.5);
  tuple->N1_DrMinPfJet = dir.make<TH1F>("N1_DrMinPfJet", ";dRMinPfJet",100,0.,5.0);
  tuple->N1_SumpTOverpT = dir.make<TH1F>("N1_SumpTOverpT", ";#Sigma p_{T} / p_{T}; Tracks / 0.025", 80, 0, 2);
  tuple->N1_Ih = dir.make<TH1F>("N1_Ih", ";I_{h} (MeV/cm)", 200, 0, dEdxM_UpLim);
  tuple->N1_MTOF = dir.make<TH1F>("N1_MTOF", ";TOF", 50, -2, 5);
  tuple->N1_I = dir.make<TH1F>("N1_I", ";I", NCuts, 0, NCuts);
  tuple->N1_TOF = dir.make<TH1F>("TOF", ";TOF", NCuts, 0, NCuts);

  tuple->NVTrack = dir.make<TH1F>("NVTrack", ";No Vertex Track; Tracks / 1", 1, 0, 1);
  tuple->N1_Stations = dir.make<TH1F>("N1_Stations", ";Stations; Tracks / 1", 1, 0, 1);
  tuple->N1_Dxy = dir.make<TH1F>("N1_Dxy", ";d_{xy} (cm); Tracks / 0.001 cm", 200, -0.1, 0.1);
  tuple->N1_Dz = dir.make<TH1F>("N1_Dz", ";d_{z} (cm); Tracks / 0.003 cm", 200, -0.3, 0.3);

  tuple->N1_PtErrOverPt = dir.make<TH1F>("N1_PtErrOverPt", ";#sigma_{p_{T}}/p_{T}", 40, 0, 1);
  tuple->N1_PtErrOverPt2 = dir.make<TH1F>("N1_PtErrOverPt2", ";#sigma_{p_{T}}/p_{T}^{2};Tracks / 7.5e-5 ", 40, 0, 0.003);
  tuple->N1_PtErrOverPtVsPt = dir.make<TH2F>("N1_PtErrOverPtVsPt", ";#sigma_{p_{T}}/p_{T};p_{T} (GeV)",  40, 0., 1., 40, 0., 4000.);
  tuple->N1_PtErrOverPtVsPt_lowPt = dir.make<TH2F>("N1_PtErrOverPtVsPt_lowPt", ";#sigma_{p_{T}}/p_{T};p_{T} (GeV)",  40, 0., 1., 50, 0., 1000.);
  tuple->N1_PtErrOverPtVsGenBeta = dir.make<TH2F>("N1_PtErrOverPtVsGenBeta", ";#sigma_{p_{T}}/p_{T};Gen #beta",  40, 0., 1., 100, 0., 1.);
  tuple->N1_PtErrOverPt2VsIas =  dir.make<TH2F>("N1_PtErrOverPt2VsIas", ";#sigma_{p_{T}}/p_{T}^{2};G_{i}^{strips};Tracks", 40, 0., 0.003, 20, 0.,1.);
  tuple->N1_PtErrOverPt2VsProbQNoL1 =  dir.make<TH2F>("N1_PtErrOverPt2VsProbQNoL1", ";#sigma_{p_{T}}/p_{T}^{2};prob_{Q,pixelAV} (pixels);Tracks", 40, 0., 0.003, 20, 0.,1.);

  tuple->N1_ProbQNoL1 = dir.make<TH1F>("N1_ProbQNoL1", ";F_{i}^{pixels};Tracks / 0.01", 40, 0., 1.);
  tuple->N1_ProbQNoL1VsIas = dir.make<TH2F>("N1_ProbQNoL1VsIas", ";F_{i}^{pixels};G_{i}^{strips}", 100, 0.0, 1.0, 100, 0.0, 1.0);
  tuple->N1_ProbXY = dir.make<TH1F>("N1_ProbXY", ";ProbXY;Tracks / 0.01", 100, 0, 1);
  tuple->N1_MiniRelIsoAll = dir.make<TH1F>("N1_MiniRelIsoAll", ";MiniRelIsoAll;Tracks / 0.01",  150, 0.0, 1.5);
  tuple->N1_MiniRelIsoAll_lowMiniRelIso = dir.make<TH1F>("N1_MiniRelIsoAll_lowMiniRelIso", ";MiniRelIsoAll;Tracks / 0.0001",  100, 0.0, 0.1);
  tuple->N1_MiniRelTkIso = dir.make<TH1F>("N1_MiniRelTkIso", ";MiniRelTkIso;Tracks / 0.01",  150, 0.0, 1.5);
  tuple->N1_MiniRelTkIso_lowMiniRelIso = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso", ";MiniRelTkIso;Tracks / 0.001",  100, 0.0, 0.1);
  tuple->N1_MiniTkIso = dir.make<TH1F>("N1_MiniTkIso", ";MiniTkIso (GeV);Tracks / 0.5",  100, 0.0, 50.);
  tuple->N1_MiniRelTkIso_lowMiniRelIso_PUA = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso_PUA", ";MiniRelTkIso (PU < 15);Tracks / 0.01",  100, 0.0, 0.1);
  tuple->N1_MiniRelTkIso_lowMiniRelIso_PUB = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso_PUB", ";MiniRelTkIso (15 =< PU < 30);Tracks / 0.01",  100, 0.0, 0.1);
  tuple->N1_MiniRelTkIso_lowMiniRelIso_PUC = dir.make<TH1F>("N1_MiniRelTkIso_lowMiniRelIso_PUC", ";MiniRelTkIso (PU >= 30);Tracks / 0.01",  100, 0.0, 0.1);
  tuple->N1_MiniTkIso_PUA = dir.make<TH1F>("N1_MiniTkIso_PUA", ";MiniTkIso (PU < 15) (GeV);Tracks / 0.3",  150, 0.0, 50.);
  tuple->N1_MiniTkIso_PUB = dir.make<TH1F>("N1_MiniTkIso_PUB", ";MiniTkIso (15 =< PU < 30) (GeV);Tracks / 0.3",  150, 0.0, 50.);
  tuple->N1_MiniTkIso_PUC = dir.make<TH1F>("N1_MiniTkIso_PUC", ";MiniTkIso (PU >= 30) (GeV);Tracks / 0.3",  150, 0.0, 50.);

  tuple->N1_PfType = dir.make<TH1F>("N1_PfType", ";;Tracks / category", 9, -0.5, 8.5);
  tuple->N1_PfType->GetXaxis()->SetBinLabel(1,"AllTracks");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(2,"PFtracks");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(3,"isElectron");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(4,"isMuon");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(5,"isPhoton");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(6,"isChHadron");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(7,"isNeutHadron");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(8,"isUndefined");
  tuple->N1_PfType->GetXaxis()->SetBinLabel(9,"notPFtrack");

  tuple->Gen_DecayLength = dir.make<TH1F>("Gen_DecayLength", "DecayLength (maybe cm);Gen candidate / 1 cm", 1000, 0., 1000.);
  tuple->Gen_Beta_Charged = dir.make<TH1F>("Beta_GenCharged", ";#beta (GenCharged);Gen candidate / 0.05", 20, 0, 1);
  tuple->Gen_Beta_Triggered = dir.make<TH1F>("Beta_Triggered", ";#beta (Triggered);Gen candidate / 0.05", 20, 0, 1);

  tuple->Gen_Binning = dir.make<TH1F>("Gen_Binning", ";Gen_Binning",1200,0.,1200.);
  tuple->Gen_pT = dir.make<TH1F>("Gen_pT", ";Generator candidate p_{T} (GeV);Gen candidate / 100 GeV", 40, 0, PtHistoUpperBound);
  tuple->Gen_Eta = dir.make<TH1F>("Gen_Eta", ";Generator candidate #eta;Gen candidate / 0.05", 52, -2.6, 2.6);
  tuple->Gen_Beta = dir.make<TH1F>("Gen_Beta", ";Generator candidate #beta;Gen candidate / 0.05", 20, 0., 1.);
  tuple->Gen_BetaGamma = dir.make<TH1F>("Gen_BetaGamma", ";Generator candidate #beta #gamma;Gen candidate/ 0.1",4500,0.,450.);
  tuple->Gen_BetaGamma_lowBetaGamma = dir.make<TH1F>("Gen_BetaGamma_lowBetaGamma", ";Generator candidate #beta #gamma;Gen candidate/ 0.1",100,0.,10.);
  
  if (doBefTrigPlots_) {
    tuple->BefTrig_ProbQNoL1 = dir.make<TH1F>("BefTrig_ProbQNoL1", ";F_{i}^{pixels};Tracks / 0.1", 10, 0., 1.);
    tuple->BefTrig_Ih = dir.make<TH1F>("BefTrig_Ih", ";I_{h} (MeV/cm)", 200, 0, dEdxM_UpLim);
    tuple->BefTrig_Ias = dir.make<TH1F>("BefTrig_Ias", ";G_{i}^{strips};Tracks / 0.05", 20, 0., 1.);
    tuple->BefTrig_TriggerMuon50VsPt_lowPt = dir.make<TH2F>("BefTrig_TriggerMuon50VsPt_lowPt", ";Muon50 triggered;Track p_{T};Tracks / bin",2,-.5,1.5,40,-0.05,200.05);
    tuple->BefTrig_TriggerMuonAllVsPt_lowPt = dir.make<TH2F>("BefTrig_TriggerMuonAllVsPt_lowPt", ";Muon50 triggered;Track p_{T};Tracks / bin",2,-.5,1.5,40,-0.05,200.05);
  }
  if (doBefPreSplots_) {
    tuple->BefPreS_MatchedMuonPt25Pt = dir.make<TH1F>("BefPreS_MatchedMuonPt25Pt", ";Matched muon ( p_{T} > 25)  p_{T} (GeV);Events / 5 GeV", 40,-0.05,200.05);
    tuple->BefPreS_RelDiffMatchedMuonPtAndTrigObjPt = dir.make<TH1F>("BefPreS_RelDiffMatchedMuonPtAndTrigObjPt", ";(Matched offline muon p_{T} - trigger object p_{T}) / trigger object p_{T};Tracks / bin", 60,-1.0,2.0);
    tuple->BefPreS_RelDiffTrigObjPtAndMatchedMuonPt = dir.make<TH1F>("BefPreS_RelDiffTrigObjPtAndMatchedMuonPt", ";(trigger object p_{T} - matched offline muon p_{T} ) / matched offline muon p_{T};Tracks / bin", 60,-1.0,2.0);
    tuple->BefPreS_NumPassedMatchingTrigObj = dir.make<TH1F>("BefPreS_NumPassedMatchingTrigObj", ";Num of trig objects passing matching;Trig objects / 1", 5, -0.5, 4.5);
    tuple->BefPreS_NumPassedMatchingTrigObjEtaCut = dir.make<TH1F>("BefPreS_NumPassedMatchingTrigObjEtaCut", ";Num of trig objects passing matching and #eta < 1;Trig objects / 1", 5, -0.5, 4.5);
    
    
    tuple->BefPreS_RelDiffMuonPtAndTrackPt = dir.make<TH1F>("BefPreS_RelDiffMuonPtAndTrackPt", ";(TuneP muon p_{T} - general track p_{T}) / general track p_{T};Tracks / bin", 60,-1.0,2.0);
    tuple->BefPreS_MuonPtVsTrackPt = dir.make<TH2F>("BefPreS_MuonPtVsTrackPt", ";TuneP muon p_{T};general track p_{T};", 100, 0.0, 4000.0, 100, 0.0, 4000.0);
    tuple->BefPreS_MuonPtOverGenPtVsTrackPtOverGenPt = dir.make<TH2F>("BefPreS_MuonPtOverGenPtVsTrackPtOverGenPt", ";TuneP muon p_{T} / gen  p_{T};general track p_{T} / gen  p_{T};", 20, 0.0, 3.0, 20, 0.0, 3.0);
    tuple->BefPreS_RelDiffMuonPtAndTruthPt = dir.make<TH1F>("BefPreS_RelDiffMuonPtAndTruthPt", ";(TuneP muon p_{T} - gen  p_{T}) / gen  p_{T};", 60,-1.0,2.0);
    tuple->BefPreS_RelDiffTrackPtAndTruthPt = dir.make<TH1F>("BefPreS_RelDiffTrackPtAndTruthPt", ";(general track p_{T}- gen  p_{T}) / gen  p_{T};", 60,-1.0,2.0);
    
    
    tuple->BefPreS_HltMatchTrackLevel = dir.make<TH1F>("BefPreS_HltMatchTrackLevel", ";;Tracks / category", 5, -0.5, 4.5);
    tuple->BefPreS_HltMatchTrackLevel->GetXaxis()->SetBinLabel(1,"All tracks");
    tuple->BefPreS_HltMatchTrackLevel->GetXaxis()->SetBinLabel(2,"Tracks matched to HLT muon");
    
    tuple->BefPreS_TriggerGenMatch = dir.make<TH1F>("BefPreS_TriggerGenMatch", ";;Events / category", 7, 0.5, 7.5);
    tuple->BefPreS_TriggerGenMatch->GetXaxis()->SetBinLabel(1,"Triggered w/ muon match");
    tuple->BefPreS_TriggerGenMatch->GetXaxis()->SetBinLabel(2,"Gen match was found");
    tuple->BefPreS_TriggerGenMatch->GetXaxis()->SetBinLabel(3,"Gen match: HSCP");
    tuple->BefPreS_TriggerGenMatch->GetXaxis()->SetBinLabel(4,"Gen match: muon");
    tuple->BefPreS_TriggerGenMatch->GetXaxis()->SetBinLabel(5,"Gen match: kaon or pion");
    tuple->BefPreS_TriggerGenMatch->GetXaxis()->SetBinLabel(6,"Gen match: else");
    tuple->BefPreS_TriggerGenMatch->GetXaxis()->SetBinLabel(7,"Gen eta < 1.0");
    
    tuple->BefPreS_TriggerType = dir.make<TH1F>("BefPreS_TriggerType", ";;Events / category", 5, -0.5, 4.5);
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(1,"Neither Muon nor MET triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(2,"Muon triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(3,"MET triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(4,"Muon OR MET triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(5,"Muon AND MET triggered");
    tuple->BefPreS_HSCPCandidateType = dir.make<TH1F>("BefPreS_HSCPCandidateType", ";;Number of generator candidate / category", 6, -0.5, 5.5);
    tuple->BefPreS_HSCPCandidateType->GetXaxis()->SetBinLabel(1,"All candidates");
    tuple->BefPreS_HSCPCandidateType->GetXaxis()->SetBinLabel(2,"Neutral HSCP");
    tuple->BefPreS_HSCPCandidateType->GetXaxis()->SetBinLabel(3,"Single-charged");
    tuple->BefPreS_HSCPCandidateType->GetXaxis()->SetBinLabel(4,"Double-charged R-hadrons");
    tuple->BefPreS_HSCPCandidateType->GetXaxis()->SetBinLabel(5,"Tau-prime (1e or 2e)");
    tuple->BefPreS_HSCPCandidateType->GetXaxis()->SetBinLabel(6,"Else");
    
    tuple->BefPreS_TriggerType = dir.make<TH1F>("BefPreS_TriggerType", ";;Events / category", 5, -0.5, 4.5);
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(1,"Neither Muon nor MET triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(2,"Muon triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(3,"MET triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(4,"Muon OR MET triggered");
    tuple->BefPreS_TriggerType->GetXaxis()->SetBinLabel(5,"Muon AND MET triggered");
    
    tuple->BefPreS_RecoHSCParticleType = dir.make<TH1F>("BefPreS_RecoHSCParticleType", ";;Tracks / category", 6, -0.5, 5.5);
    tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(1,"globalMuon");
    tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(2,"trackerMuon");
    tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(3,"matchedStandAloneMuon");
    tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(4,"standAloneMuon");
    tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(5,"innerTrack");
    tuple->BefPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(6,"unknown");

    tuple->BefPreS_GenPtVsdRMinGen = dir.make<TH2F>("BefPreS_GenPtVsdRMinGen", ";GenPt;dRMinGen", 40, 0., PtHistoUpperBound, 100, 0., 1.);
    tuple->BefPreS_GendRMin = dir.make<TH1F>("BefPreS_GendRMin", ";dR_min;Gen candidate / 0.032",100,0.,3.2);
    tuple->BefPreS_GenPtVsdRMinGenPostCut = dir.make<TH2F>("BefPreS_GenPtVsdRMinGenPostCut", ";GenPt (GeV);dRMinGen (after cut)", 40, 0., PtHistoUpperBound, 50, 0., 0.05);
    tuple->BefPreS_GenPtVsGenMinPt = dir.make<TH2F>("BefPreS_GenPtVsGenMinPt", ";GenPtVsGenMinPt", 40, 0., PtHistoUpperBound, 100, 0, 1.);
    tuple->BefPreS_GenPtVsRecoPt = dir.make<TH2F>("BefPreS_GenPtVsRecoPt", ";Generator p_{T} (GeV);Reco p_{T} (GeV)", 40, 0., PtHistoUpperBound, 40, 0, PtHistoUpperBound);
    
    tuple->BefPreS_RatioCleanAndAllStripsClu = dir.make<TH1F>("BefPreS_RatioCleanAndAllStripsClu",";Clean / all strips clu;Tracks / 0.05",20,-0.05,1.05);
    tuple->BefPreS_RatioCleanAndAllPixelClu = dir.make<TH1F>("BefPreS_RatioCleanAndAllPixelClu",";Clean / all pixel clu;Tracks / 0.05",20,-0.05,1.05);

    tuple->BefPreS_PfType = dir.make<TH1F>("BefPreS_PfType", ";;Tracks / category", 9, -0.5, 8.5);
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(1,"AllTracks");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(2,"PFtracks");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(3,"isElectron");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(4,"isMuon");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(5,"isPhoton");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(6,"isChHadron");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(7,"isNeutHadron");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(8,"isUndefined");
    tuple->BefPreS_PfType->GetXaxis()->SetBinLabel(9,"notPFtrack");
    
    tuple->BefPreS_MassT = dir.make<TH1F>("BefPreS_MassT", ";m_{T} (GeV);Tracks / 5 GeV", 50, 0.0, 250.0);
    tuple->BefPreS_MassT_highMassT = dir.make<TH1F>("BefPreS_MassT_highMassT", ";m_{T} (GeV);Tracks / 10 GeV", 250, 0.0, 2500.0);

    tuple->BefPreS_MiniRelIsoAll = dir.make<TH1F>("BefPreS_MiniRelIsoAll", ";MiniRelIsoAll;Tracks / 0.01", 150, 0.0, 1.5);
    tuple->BefPreS_MiniRelIsoChg = dir.make<TH1F>("BefPreS_MiniRelIsoChg", ";MiniRelIsoChg;Tracks / 0.01",  150, 0.0, 1.5);
    tuple->BefPreS_MiniRelTkIso = dir.make<TH1F>("BefPreS_MiniRelTkIso", ";MiniRelTkIso;Tracks / 0.01",  150, 0.0, 1.5);
    tuple->BefPreS_MiniTkIso = dir.make<TH1F>("BefPreS_MiniTkIso", ";MiniTkIso;Tracks / 0.05",  100, 0.0, 50);
    
    tuple->BefPreS_RecoPFMET = dir.make<TH1F>("BefPreS_RecoPFMET", ";RecoPFMET (GeV);Events / 10 GeV",  200, 0.0, 2000.0);
    tuple->BefPreS_RecoPfHT = dir.make<TH1F>("BefPreS_RecoPfHT", ";RecoPFHT (GeV);Events / 10 GeV",  200, 0.0, 2000.0);
    tuple->BefPreS_RecoPfJetsNum = dir.make<TH1F>("BefPreS_RecoPfJetsNum", ";Number of PF jets;Tracks / 1",  16, -0.5, 15.5);
    tuple->BefPreS_CaloJetsNum = dir.make<TH1F>("BefPreS_CaloJetsNum", ";Number of calo jets;Tracks / 1",  16, -0.5, 15.5);

    tuple->BefPreS_Chi2oNdof = dir.make<TH1F>("BefPreS_Chi2oNdof", ";Chi2oNdof;Tracks / bin", 20, 0, 20);
    // This should just be a 2-bin plot where high-purity is not present = 0 or present = 1
    tuple->BefPreS_Qual = dir.make<TH1F>("BefPreS_Qual", ";;Tracks / category", 2, -0.5, 1.5);
    tuple->BefPreS_Qual->GetXaxis()->SetBinLabel(1,"Not-HighPurity");
    tuple->BefPreS_Qual->GetXaxis()->SetBinLabel(2,"HighPurity");

    tuple->BefPreS_TNOH_PUA = dir.make<TH1F>("BefPreS_TNOH_PUA", ";TNOH (PUA);Tracks / bin",  40, -0.5, 39.5);
    tuple->BefPreS_TNOH_PUB = dir.make<TH1F>("BefPreS_TNOH_PUB", ";TNOH (PUB);Tracks / bin", 40, -0.5, 39.5);
    tuple->BefPreS_TNOHFraction = dir.make<TH1F>("BefPreS_TNOHFraction", ";Number of valid hit fraction;Tracks / 0.05", 20, 0., 1.);
    tuple->BefPreS_TNOPH = dir.make<TH1F>("BefPreS_TNOPH", ";Number of pixel hits;Tracks / bin", 8, -0.5, 7.5);
    tuple->BefPreS_TNOHFractionTillLast = dir.make<TH1F>("BefPreS_TNOHFractionTillLast",";TNOHFractionTillLastTracks / bin", 50, 0., 1);
    tuple->BefPreS_TNOMHTillLast = dir.make<TH1F>("BefPreS_TNOMHTillLast", ";TNOMHTillLast;Tracks / bin", 20, 0, 20);
    tuple->BefPreS_Eta = dir.make<TH1F>("BefPreS_Eta", ";#eta;Tracks / 0.05", 52, -2.6, 2.6);
    tuple->BefPreS_TNOM = dir.make<TH1F>("BefPreS_TNOM", ";Number of measurement;Tracks / 1", 40, -0.5, 39.5);
    tuple->BefPreS_TNOM_PUA = dir.make<TH1F>("BefPreS_TNOM_PUA", ";Number of measurements (low PU); Tracks / 1 ", 40, -0.5, 39.5);
    tuple->BefPreS_TNOM_PUB = dir.make<TH1F>("BefPreS_TNOM_PUB", ";Number of measurements (high PU); Tracks / 1 ", 40, -0.5, 39.5);
    tuple->BefPreS_nDof = dir.make<TH1F>("BefPreS_nDof", ";Number of DF;Tracks / bin", 40, -0.5, 39.5);
    tuple->BefPreS_TOFError = dir.make<TH1F>("BefPreS_TOFError", ";TOFError;Tracks / bin", 25, 0, 0.25);
    tuple->BefPreS_PtErrOverPt = dir.make<TH1F>("BefPreS_PtErrOverPt", ";#sigma_{p_{T}}/p_{T};Tracks / 0.025", 40, 0, 1);
    tuple->BefPreS_PtErrOverPt2 = dir.make<TH1F>("BefPreS_PtErrOverPt2", ";#sigma_{p_{T}}/p_{T}^{2};Tracks / 7.5e-5 ", 40, 0, 0.003);
    tuple->BefPreS_Pt = dir.make<TH1F>("BefPreS_Pt", ";p_{T} (GeV);Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_lowPt = dir.make<TH1F>("BefPreS_Pt_lowPt", ";p_{T} (GeV);Tracks / 10 GeV", 50, 0., 500);
    tuple->BefPreS_Ias = dir.make<TH1F>("BefPreS_Ias", ";G_{i}^{strips};Tracks / 0.1", 10, 0., 1.);
    tuple->BefPreS_IasForStatus91 = dir.make<TH1F>("BefPreS_IasForStatus91", ";G_{i}^{strips};Tracks / 0.1", 10, 0., 1.);
    tuple->BefPreS_IasForStatusNot91 = dir.make<TH1F>("BefPreS_IasForStatusNot91", ";G_{i}^{strips};Tracks / 0.1", 10, 0., 1.);
    tuple->BefPreS_Ih = dir.make<TH1F>("BefPreS_Ih", ";I_{h} (MeV/cm)", 200, 0, dEdxM_UpLim);
    tuple->BefPreS_MTOF = dir.make<TH1F>("BefPreS_MTOF", ";MTOF;Tracks / bin", 50, -2, 5);
    tuple->BefPreS_TIsol = dir.make<TH1F>("BefPreS_TIsol", ";#Sigma_{R<0.3} p_{T} - p_{T,cand} (GeV);Tracks / 4 GeV", 25, 0, 100);
    tuple->BefPreS_EoP = dir.make<TH1F>("BefPreS_EoP", ";PF calo energy / momentum; Tracks / 0.06", 25, 0, 1.5);
    tuple->BefPreS_ECalEoP = dir.make<TH1F>("BefPreS_ECalEoP", ";PF ECal energy / momentum; Tracks / 0.06", 25, 0, 1.5);
    tuple->BefPreS_HCalEoP = dir.make<TH1F>("BefPreS_HCalEoP", ";PF HCal energy / momentum; Tracks / 0.06", 25, 0, 1.5);
    tuple->BefPreS_SumpTOverpT = dir.make<TH1F>("BefPreS_SumpTOverpT", ";#Sigma p_{T} / p_{T};Tracks / bin", 80, 0.0, 2.0);
    tuple->BefPreS_LastHitDXY = dir.make<TH1F>("BefPreS_LastHitDXY", ";LastHitDXY;Tracks / bin", 75, 0, 150);
    tuple->BefPreS_LastHitD3D = dir.make<TH1F>("BefPreS_LastHitD3D", ";LastHitD3D;Tracks / bin", 175, 0, 350);
    tuple->BefPreS_P = dir.make<TH1F>("BefPreS_P", ";P;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt = dir.make<TH1F>("BefPreS_Pt", ";Pt;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_PUA = dir.make<TH1F>("BefPreS_Pt_PUA", ";p_{T} (PUA);Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_PUB = dir.make<TH1F>("BefPreS_Pt_PUB", ";p_{T} (PUB);Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_Cosmic = dir.make<TH1F>("BefPreS_Pt_Cosmic", ";p_{T} (Cosmic);Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_DT = dir.make<TH1F>("BefPreS_Pt_DT", ";Pt_DT;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_CSC = dir.make<TH1F>("BefPreS_Pt_CSC", ";Pt_CSC;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Ias = dir.make<TH1F>("BefPreS_Ias", ";G_{i}^{strips};Tracks / bin", 100, 0, dEdxS_UpLim);
    tuple->BefPreS_Ias_PUA = dir.make<TH1F>("BefPreS_Ias_PUA", ";G_{i}^{strips} (PUA);Tracks / bin", 100, 0, dEdxS_UpLim);
    tuple->BefPreS_Ias_PUB = dir.make<TH1F>("BefPreS_Ias_PUB", ";G_{i}^{strips} (PUB);Tracks / bin", 100, 0, dEdxS_UpLim);
    tuple->BefPreS_Ias_Cosmic = dir.make<TH1F>("BefPreS_Ias_Cosmic", ";G_{i}^{strips}_Cosmic;Tracks / bin", 100, 0, dEdxS_UpLim);
    tuple->BefPreS_Ih_Cosmic = dir.make<TH1F>("BefPreS_Ih_Cosmic", ";I_{h} (Cosmic);Tracks / bin", 200, 0, dEdxM_UpLim);
    tuple->BefPreS_Ih = dir.make<TH1F>("BefPreS_Ih", ";I_{h} (MeV/cm);Tracks / bin", 200, 0, dEdxM_UpLim);
    tuple->BefPreS_Ih_PUA = dir.make<TH1F>("BefPreS_Ih_PUA", ";I_{h} (PUA);Tracks / bin", 200, 0, dEdxM_UpLim);
    tuple->BefPreS_Ih_PUB = dir.make<TH1F>("BefPreS_Ih_PUB", ";I_{h} (PUB);Tracks / bin", 200, 0, dEdxM_UpLim);
    tuple->BefPreS_TOF = dir.make<TH1F>("BefPreS_TOF", ";TOF;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_TOF_PUA = dir.make<TH1F>("BefPreS_TOF_PUA", ";TOF_PUA;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_TOF_PUB = dir.make<TH1F>("BefPreS_TOF_PUB", ";TOF_PUB;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_TOF_DT = dir.make<TH1F>("BefPreS_TOF_DT", ";TOF_DT;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_TOF_CSC = dir.make<TH1F>("BefPreS_TOF_CSC", ";TOF_CSC;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_dR_NVTrack = dir.make<TH1F>("BefPreS_dR_NVTrack", ";dR_NVTrack;Tracks / 0.025", 40, 0, 1);
    tuple->BefPreS_MatchedStations = dir.make<TH1F>("BefPreS_MatchedStations", ";MatchedStations;Tracks / bin", 8, -0.5, 7.5);
    tuple->BefPreS_InnerInvPtDiff = dir.make<TH1F>("BefPreS_InnerInvPtDiff", ";InnerInvPtDiff;Tracks / bin", 120, -4, 4);
    tuple->BefPreS_Phi = dir.make<TH1F>("BefPreS_Phi", ";Phi", 50, -3.14, 3.14);
    tuple->BefPreS_TimeAtIP = dir.make<TH1F>("BefPreS_TimeAtIP", ";TimeAtIP;Tracks / bin", 50, -100, 100);
    tuple->BefPreS_OpenAngle = dir.make<TH1F>("BefPreS_OpenAngle", ";OpenAngle;Tracks / bin", 50, -0.3, 3.15);
    tuple->BefPreS_OpenAngle_Cosmic = dir.make<TH1F>("BefPreS_OpenAngle_Cosmic", ";OpenAngle_Cosmic;Tracks / bin", 50, -0.3, 3.15);

    tuple->BefPreS_NVertex = dir.make<TH1F>("BefPreS_NVertex", ";NVertex;Tracks / 1", 50, -0.5, 49.5);
    tuple->BefPreS_NVertex_NoEventWeight = dir.make<TH1F>("BefPreS_NVertex_NoEventWeight", ";N_{Vertex} (NoEventWeight);Tracks / 1", 50, -0.5, 49.5);
    tuple->BefPreS_PV = dir.make<TH1F>("BefPreS_PV", ";N_{PV};Tracks / 1", 60, 0, 60);
    tuple->BefPreS_PV_NoEventWeight = dir.make<TH1F>("BefPreS_PV_NoEventWeight", ";N_{PV} (NoEventWeight);Tracks / 1", 60, -0.5, 59.5);
    tuple->BefPreS_NOMoNOH = dir.make<TH1F>("BefPreS_NOMoNOH", ";Num of measurment / num of hits;Tracks / bin",10,0.,1.0);
    tuple->BefPreS_NOMoNOHvsPV = dir.make<TProfile>("BefPreS_NOMoNOHvsPV", ";NOMoNOHvsPV;Tracks / bin", 60, 0, 60);
    tuple->BefPreS_Dz = dir.make<TH1F>("BefPreS_Dz",";d_{z} (cm);Tracks / 0.003 cm", 200, -0.3, 0.3);
    tuple->BefPreS_Dxy = dir.make<TH1F>("BefPreS_Dxy","d_{xy} (cm);Tracks / 0.001 cm", 200, -0.1, 0.1);

    tuple->BefPreS_SegSep = dir.make<TH1F>("BefPreS_SegSep", ";SegSep;Tracks / bin", 50, 0., 2.5);
    tuple->BefPreS_SegMinEtaSep = dir.make<TH1F>("BefPreS_SegMinEtaSep", ";SegMinEtaSep;Tracks / bin", 50, -1., 1.);
    tuple->BefPreS_SegMinPhiSep = dir.make<TH1F>("BefPreS_SegMinPhiSep", ";SegMinPhiSep;Tracks / bin", 50, -3.3, 3.3);
    tuple->BefPreS_SegMinEtaSep_FailDz = dir.make<TH1F>("BefPreS_SegMinEtaSep_FailDz", ";SegMinEtaSep_FailDz;Tracks / bin", 50, -1., 1.);
    tuple->BefPreS_SegMinEtaSep_PassDz = dir.make<TH1F>("BefPreS_SegMinEtaSep_PassDz", ";SegMinEtaSep_PassDz;Tracks / bin", 50, -1., 1.);
    tuple->BefPreS_Dz_FailSep = dir.make<TH1F>("BefPreS_Dz_FailSep", ";Dz_FailSep;Tracks / bin", 50, -150, 150);

    tuple->BefPreS_Dxy_Cosmic = dir.make<TH1F>("BefPreS_Dxy_Cosmic", ";Dxy_Cosmic;Tracks / bin", 150, -IPbound, IPbound);
    tuple->BefPreS_Dz_Cosmic = dir.make<TH1F>("BefPreS_Dz_Cosmic", ";Dz_Cosmic;Tracks / bin", 150, -IPbound, IPbound);
    tuple->BefPreS_Dz_CSC = dir.make<TH1F>("BefPreS_Dz_CSC", ";Dz_CSC;Tracks / bin", 150, -IPbound, IPbound);
    tuple->BefPreS_Dz_DT = dir.make<TH1F>("BefPreS_Dz_DT", ";Dz_DT;Tracks / bin", 150, -IPbound, IPbound);
    tuple->BefPreS_Pt_FailDz = dir.make<TH1F>("BefPreS_Pt_FailDz", ";Pt_FailDz;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_FailDz_DT = dir.make<TH1F>("BefPreS_Pt_FailDz_DT", ";Pt_FailDz_DT;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_Pt_FailDz_CSC = dir.make<TH1F>("BefPreS_Pt_FailDz_CSC", ";Pt_FailDz_CSC;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->BefPreS_TOF_FailDz = dir.make<TH1F>("BefPreS_TOF_FailDz", ";TOF_FailDz;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_TOF_FailDz_DT = dir.make<TH1F>("BefPreS_TOF_FailDz_DT", ";TOF_FailDz_DT;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_TOF_FailDz_CSC = dir.make<TH1F>("BefPreS_TOF_FailDz_CSC", ";TOF_FailDz_CSC;Tracks / 0.05", 120, -1, 5);
    tuple->BefPreS_GenPtVsRecoPt = dir.make<TH2F>("BefPreS_GenPtVsRecoPt", ";GenPt;RecoPt;Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 40, 0, PtHistoUpperBound);
    tuple->BefPreS_PtErrOverPtVsPt = dir.make<TH2F>("BefPreS_PtErrOverPtVsPt",  ";#sigma_{p_{T}}/p_{T};p_{T};Tracks / bin",  40, 0., 1., 40, 0., 4000);
    tuple->BefPreS_PtErrOverPtVsPtErrOverPt2 = dir.make<TH2F>("BefPreS_PtErrOverPtVsPtErrOverPt2",  ";#sigma_{p_{T}}/p_{T};p_{T}^{2};Tracks / bin",  40, 0., 1., 40, 0., 0.003);

    tuple->BefPreS_ProbQ = dir.make<TH1F>("BefPreS_ProbQ", ";F_{i}^{pixels};Tracks / bin", 100, 0, 1);
    tuple->BefPreS_ProbXY = dir.make<TH1F>("BefPreS_ProbXY", ";Prob_{XY,pixelAV} (pixels);Tracks / bin", 100, 0, 1);
    tuple->BefPreS_ProbQNoL1 = dir.make<TH1F>("BefPreS_ProbQNoL1", ";F_{i}^{pixels};Tracks / bin", 100, 0, 1);
    tuple->BefPreS_ProbXYNoL1 = dir.make<TH1F>("BefPreS_ProbXYNoL1", ";Prob_{XY,pixelAV} (pixels, no-L1);Tracks / bin", 100, 0, 1);
    tuple->BefPreS_MassErr = dir.make<TH1F>("BefPreS_MassErr", ";MassErr;Tracks / bin", 50, 0., 10.);
    tuple->BefPreS_ProbQVsIas = dir.make<TH2F>("BefPreS_ProbQVsIas", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin", 100, 0.0, 1.0, 100, 0.0, 1.0);

    tuple->BefPreS_EtaVsIas = dir.make<TH2F>("BefPreS_EtaVsIas", ";#eta;G_{i}^{strips};Tracks / bin", 50, -3, 3, 10, 0., 1.);
    tuple->BefPreS_EtaVsIh = dir.make<TH2F>("BefPreS_EtaVsIh", ";#eta;I_{h} (MeV/cm);Tracks / bin", 50, -3, 3, 100, 0, dEdxM_UpLim);
    tuple->BefPreS_EtaVsP = dir.make<TH2F>("BefPreS_EtaVsP", ";#eta;P (GeV);Tracks / bin", 50, -3, 3, 40, 0, PtHistoUpperBound);
    tuple->BefPreS_EtaVsPt = dir.make<TH2F>("BefPreS_EtaVsPt", ";#eta;p_{T} (GeV);Tracks / bin", 50, -3, 3, 40, 0, PtHistoUpperBound);
    tuple->BefPreS_EtaVsTOF = dir.make<TH2F>("BefPreS_EtaVsTOF", ";#eta;TOF;Tracks / bin", 50, -3, 3, 50, 0., 3);
    tuple->BefPreS_EtaVsNBH = dir.make<TH2F>("BefPreS_EtaVsNBH", ";#eta;Number of bad hits;Tracks / bin", 60, -3., 3., 24, 0., 24.);
    tuple->BefPreS_EtaVsDz = dir.make<TH2F>("BefPreS_EtaVsDz", ";#eta;D_{z} (cm);Tracks / bin", 50, -3, 3, 50, -IPbound, IPbound);
    tuple->BefPreS_PVsIas = dir.make<TH2F>("BefPreS_PVsIas", ";P;G_{i}^{strips};Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 100, 0, dEdxS_UpLim);
    tuple->BefPreS_IhVsIas = dir.make<TH2F>("BefPreS_IhVsIas", ";I_{h} (MeV/cm);G_{i}^{strips};Tracks / bin", 50, 0., dEdxM_UpLim, 20, 0, dEdxS_UpLim);
    tuple->BefPreS_PVsIh = dir.make<TH2F>("BefPreS_PVsIh", ";P;I_{h} (MeV/cm);Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 50, 0., dEdxM_UpLim);
    tuple->BefPreS_PtVsIas = dir.make<TH2F>("BefPreS_PtVsIas", ";p_{T} (GeV);G_{i}^{strips};Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 10, 0., 1.);
    tuple->BefPreS_PtVsIh = dir.make<TH2F>("BefPreS_PtVsIh", ";p_{T} (GeV);I_{h} (MeV/cm);Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 100, 0, dEdxM_UpLim);
    tuple->BefPreS_PtVsTOF = dir.make<TH2F>("BefPreS_PtVsTOF", ";Pt;TOF;Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 50, 0., 5);
      //tuple->BefPreS_TOFIs = dir.make<TH2F>("BefPreS_TOFIs", ";TOFIs", 100, 1, 5, 100, 0, dEdxS_UpLim);
    tuple->BefPreS_TOFVsIs = dir.make<TH2F>("BefPreS_TOFVsIs", ";TOF;G_{i}^{strips};Tracks / bin", 50, 0., 5, 10, 0., 1.);
      //tuple->BefPreS_TOFIm = dir.make<TH2F>("BefPreS_TOFIh", ";TOFIh", 100, 1, 5, 200, 0, dEdxM_UpLim);
    tuple->BefPreS_TOFVsIh = dir.make<TH2F>("BefPreS_TOFVsIh", ";TOF;I_{h} (MeV/cm);Tracks / bin", 50, 0., 5, 100, 0, dEdxM_UpLim);

    tuple->BefPreS_CluProbHasFilled = dir.make<TH1F>("BefPreS_CluProbHasFilled", ";;Clusters / category", 2, -0.5, 1.5);
    tuple->BefPreS_CluProbHasFilled->GetXaxis()->SetBinLabel(1,"ProbHasFailed");
    tuple->BefPreS_CluProbHasFilled->GetXaxis()->SetBinLabel(2,"ProbHasFilled");

    tuple->BefPreS_CluProbQVsPixelLayer = dir.make<TH2F>("BefPreS_CluProbQVsPixelLayer", ";Cluster Prob_{Q} (pixels);Layer",20,0.,1.,4,0.,4.);
    tuple->BefPreS_CluProbXYVsPixelLayer = dir.make<TH2F>("BefPreS_CluProbXYVsPixelLayer", ";Cluster Prob_{XZ,pixelAV} (pixels);Layer",100,0.,1.,4,0.,4.);
    tuple->BefPreS_CluNormChargeVsPixelLayer = dir.make<TH2F>("BefPreS_CluNormChargeVsPixelLayer", ";CluNormCharge (e/um);Layer",100,0.,600.,4,0.,4.);
    tuple->BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma", ";CluNormCharge (e/um);Layer",100,0.,600.,4,0.,4.);
    tuple->BefPreS_CluSizeVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeVsPixelLayer", ";CluSize;Layer",10,0.,10.,4,0.,4.);
    tuple->BefPreS_CluSizeXVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeXVsPixelLayer", ";CluSizeX;Layer",10,0.,10.,4,0.,4.);
    tuple->BefPreS_CluSizeYVsPixelLayer = dir.make<TH2F>("BefPreS_CluSizeYVsPixelLayer", ";CluSizeY;Layer",10,0.,10.,4,0.,4.);
    tuple->BefPreS_CluSpecInCPEVsPixelLayer = dir.make<TH2F>("BefPreS_CluSpecInCPEVsPixelLayer", ";CluSpecInCPE;Layer",4,0.,4.,4,0.,4.);

    tuple->BefPreS_CluCotBetaVsPixelLayer_lowProbXY = dir.make<TH2F>("BefPreS_CluCotBetaVsPixelLayer_lowProbXY", ";CotBeta;Layer",200,-10.,10.,4,0.,4.);
    tuple->BefPreS_CluCotAlphaVsPixelLayer_lowProbXY = dir.make<TH2F>("BefPreS_CluCotAlphaVsPixelLayer_lowProbXY", ";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
    tuple->BefPreS_CluCotBetaVsPixelLayer = dir.make<TH2F>("BefPreS_CluCotBetaVsPixelLayer", ";CotBeta;Layer",200,-10.,10.,4,0.,4.);
    tuple->BefPreS_CluCotAlphaVsPixelLayer = dir.make<TH2F>("BefPreS_CluCotAlphaVsPixelLayer", ";CotAlpha;Layer",100,-1.,1.,4,0.,4.);

    tuple->BefPreS_CluNormChargeVsStripLayer_lowBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_lowBetaGamma", ";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
    tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma", ";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
    tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_Stat91 = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_Stat91", ";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
    tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatNot91 = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatNot91", ";CluNormCharge;Layer",600,0.,600.,24,-0.5,23.5);
    tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatHigherThan2 = dir.make<TH2F>("BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatHigherThan2", ";CluNormCharge;Layer",600,0.,600.,20,0.,20.);

    tuple->BefPreS_dRMinPfJet = dir.make<TH1F>("BefPreS_dRMinPfJet", ";dRMinPfJet;Tracks / bin",100,0.,5.0);
    tuple->BefPreS_dRMinPfJetVsIas =  dir.make<TH2F>("BefPreS_dRMinPfJetVsIas", ";dRMinPfJet;G_{i}^{strips};Tracks / bin",100,0.,5.0,10,0.,1.);
    tuple->BefPreS_dRMinCaloJet= dir.make<TH1F>("BefPreS_dRMinCaloJet", ";dRMinCaloJet;Tracks / bin",100,0.,5.0);
    tuple->BefPreS_dRMinCaloJetVsIas =  dir.make<TH2F>("BefPreS_dRMinCaloJetVsIas", ";dRMinCaloJet;G_{i}^{strips};Tracks / bin",100,0.,5.0,10,0.,1.);
    tuple->BefPreS_genGammaBetaVsProbXYNoL1 =  dir.make<TH2F>("BefPreS_genGammaBetaVsProbXYNoL1", ";#gamma #beta;ProbXYNoL1",10,0.,1.3,20,0.,1.);
    tuple->BefPreS_dRVsPtPfJet = dir.make<TH2F>("BefPreS_dRVsPtPfJet", ";dR(cand,jet);p_{T} (GeV)",100,0.,1.5,100,0.,1000.);
    tuple->BefPreS_dRVsdPtPfCaloJet = dir.make<TH2F>("BefPreS_dRVsdPtPfCaloJet", ";dRmin;dPtPfCaloJet",100,0.,1.5,20,0.,100.);
    tuple->BefPreS_GenBeta = dir.make<TH1F>("BefPreS_GenBeta", ";#beta;Tracks / bin", 20, 0., 1.);
    tuple->BefPreS_NumCandidates = dir.make<TH1F>("BefPreS_NumCandidates", ";Number of HSCP candidates;Events / bin", 11, -0.5, 10.5);
    
    tuple->BefPreS_TriggerMuon50VsBeta = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaA = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaA", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaA_BetaUp = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaA_BetaUp", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaA_BetaDown = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaA_BetaDown", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaB = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaB", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaB_BetaUp = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaB_BetaUp", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaB_BetaDown = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaB_BetaDown", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaC = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaC", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaC_BetaUp = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaC_BetaUp", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuon50VsBeta_EtaC_BetaDown = dir.make<TH2F>("BefPreS_TriggerMuon50VsBeta_EtaC_BetaDown", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    
    tuple->BefPreS_TriggerMuonAllVsBeta = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaA = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaA", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaA_BetaUp = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaA_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaA_BetaDown = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaA_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaB = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaB", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaB_BetaUp = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaB_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaB_BetaDown = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaB_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaC = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaC", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaC_BetaUp = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaC_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMuonAllVsBeta_EtaC_BetaDown = dir.make<TH2F>("BefPreS_TriggerMuonAllVsBeta_EtaC_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    
    tuple->BefPreS_TriggerMuon50VsPt = dir.make<TH2F>("BefPreS_TriggerMuon50VsPt", ";Muon50 triggered;Track p_{T};Tracks / bin",2,-.5,1.5,50,-0.05,1000.05);
    tuple->BefPreS_TriggerMuonAllVsPt = dir.make<TH2F>("BefPreS_TriggerMuonAllVsPt", ";MuonAll triggered;Track p_{T};Tracks / bin",2,-.5,1.5,50,-0.05,1000.05);
    
    tuple->BefPreS_TriggerMETallVsBeta = dir.make<TH2F>("BefPreS_TriggerMETallVsBeta", ";OR of MET triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->BefPreS_TriggerMETallVsMet = dir.make<TH2F>("BefPreS_TriggerMETallVsMet", ";OR of MET triggered;MET (GeV);Tracks / bin",2,-.5,1.5,50,-0.05,2000.05);
    tuple->BefPreS_TriggerMETallVsHT = dir.make<TH2F>("BefPreS_TriggerMETallVsHT", ";OR of MET triggered;H_{T} (GeV);Tracks / bin",2,-.5,1.5,50,-0.05,2000.05);
    tuple->BefPreS_TriggerMETallVsMetOverHt = dir.make<TH2F>("BefPreS_TriggerMETallVsMetOverHt", ";OR of MET triggered;MET / H_{T};Tracks / bin",2,-.5,1.5,30,-0.05,2.95);
    tuple->BefPreS_TriggerMETallVsMetVsHT = dir.make<TH3F>("BefPreS_TriggerMETallVsMetVsHT", ";OR of MET triggered;MET (GeV);H_{T} (GeV);Tracks / bin",2,-.5,1.5,50,-0.05,2000.05,50,-0.05,2000.05);
  }
  
  if (createGiTemplates_) {
    tuple->Calibration_GiTemplate =  dir.make<TH3F>("Calibration_GiTemplate", ";Module geometry;Path lenght (mm);Path normalised charge (ke / cm)", 15, 1.0, 16.0, 42, 0.2, 1.6, 500, 0.0, 5000.0);
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(1,"IB1");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(2,"IB2");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(3,"OB1");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(4,"OB2");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(5,"W1A");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(6,"W2A");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(7,"W3A");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(8,"W1B");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(9,"W2B");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(10,"W3B");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(11,"W4");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(12,"W5");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(13,"W6");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(14,"W7");
    tuple->Calibration_GiTemplate->GetXaxis()->SetBinLabel(15,"Pixels");

    tuple->Calibration_GiTemplate_PU_1 =  dir.make<TH3F>("Calibration_GiTemplate_PU_1", ";Module geometry;Path lenght (mm);Path normalised charge (ke / cm)", 15, 1.0, 16.0, 42, 0.2, 1.6, 500, 0.0, 5000.0);
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(1,"IB1");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(2,"IB2");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(3,"OB1");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(4,"OB2");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(5,"W1A");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(6,"W2A");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(7,"W3A");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(8,"W1B");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(9,"W2B");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(10,"W3B");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(11,"W4");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(12,"W5");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(13,"W6");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(14,"W7");
    tuple->Calibration_GiTemplate_PU_1->GetXaxis()->SetBinLabel(15,"Pixels");
    tuple->Calibration_GiTemplate_PU_2 =  dir.make<TH3F>("Calibration_GiTemplate_PU_2", ";Module geometry;Path lenght (mm);Path normalised charge (ke / cm)", 15, 1.0, 16.0, 42, 0.2, 1.6, 500, 0.0, 5000.0);
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(1,"IB1");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(2,"IB2");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(3,"OB1");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(4,"OB2");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(5,"W1A");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(6,"W2A");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(7,"W3A");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(8,"W1B");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(9,"W2B");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(10,"W3B");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(11,"W4");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(12,"W5");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(13,"W6");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(14,"W7");
    tuple->Calibration_GiTemplate_PU_2->GetXaxis()->SetBinLabel(15,"Pixels");
    tuple->Calibration_GiTemplate_PU_3 =  dir.make<TH3F>("Calibration_GiTemplate_PU_3", ";Module geometry;Path lenght (mm);Path normalised charge (ke / cm)", 15, 1.0, 16.0, 42, 0.2, 1.6, 500, 0.0, 5000.0);
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(1,"IB1");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(2,"IB2");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(3,"OB1");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(4,"OB2");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(5,"W1A");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(6,"W2A");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(7,"W3A");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(8,"W1B");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(9,"W2B");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(10,"W3B");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(11,"W4");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(12,"W5");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(13,"W6");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(14,"W7");
    tuple->Calibration_GiTemplate_PU_3->GetXaxis()->SetBinLabel(15,"Pixels");
    tuple->Calibration_GiTemplate_PU_4 =  dir.make<TH3F>("Calibration_GiTemplate_PU_4", ";Module geometry;Path lenght (mm);Path normalised charge (ke / cm)", 15, 1.0, 16.0, 42, 0.2, 1.6, 500, 0.0, 5000.0);
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(1,"IB1");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(2,"IB2");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(3,"OB1");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(4,"OB2");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(5,"W1A");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(6,"W2A");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(7,"W3A");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(8,"W1B");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(9,"W2B");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(10,"W3B");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(11,"W4");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(12,"W5");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(13,"W6");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(14,"W7");
    tuple->Calibration_GiTemplate_PU_4->GetXaxis()->SetBinLabel(15,"Pixels");
    tuple->Calibration_GiTemplate_PU_5 =  dir.make<TH3F>("Calibration_GiTemplate_PU_5", ";Module geometry;Path lenght (mm);Path normalised charge (ke / cm)", 15, 1.0, 16.0, 42, 0.2, 1.6, 500, 0.0, 5000.0);
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(1,"IB1");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(2,"IB2");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(3,"OB1");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(4,"OB2");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(5,"W1A");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(6,"W2A");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(7,"W3A");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(8,"W1B");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(9,"W2B");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(10,"W3B");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(11,"W4");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(12,"W5");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(13,"W6");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(14,"W7");
    tuple->Calibration_GiTemplate_PU_5->GetXaxis()->SetBinLabel(15,"Pixels");
  }

  if (doPostPreSplots_) {
    tuple->PostPreS_NumCandidates = dir.make<TH1F>("PostPreS_NumCandidates", ";Number of HSCP candidates;Events / bin", 11, -0.5, 10.5);
    tuple->PostPreS_RelDiffMuonPtAndTrackPt = dir.make<TH1F>("PostPreS_RelDiffMuonPtAndTrackPt", ";(TuneP muon p_{T} - general track p_{T}) / general track p_{T};Tracks / bin", 60,-1.0,2.0);
    tuple->PostPreS_MuonPtVsTrackPt = dir.make<TH2F>("PostPreS_MuonPtVsTrackPt", ";TuneP muon p_{T};general track p_{T};", 100, 0.0, 4000.0, 100, 0.0, 4000.0);
    tuple->PostPreS_MuonPtOverGenPtVsTrackPtOverGenPt = dir.make<TH2F>("PostPreS_MuonPtOverGenPtVsTrackPtOverGenPt", ";TuneP muon p_{T} / gen  p_{T};general track p_{T} / gen  p_{T};", 20, 0.0, 3.0, 20, 0.0, 3.0);
    
    tuple->PostPreS_RelDiffMuonPtAndTruthPt = dir.make<TH1F>("PostPreS_RelDiffMuonPtAndTruthPt", ";(TuneP muon p_{T} - gen  p_{T}) / gen  p_{T};Tracks / bin", 60,-1.0,2.0);
    tuple->PostPreS_RelDiffTrackPtAndTruthPt = dir.make<TH1F>("PostPreS_RelDiffTrackPtAndTruthPt", ";(general track p_{T} - gen  p_{T}) / gen  p_{T};Tracks / bin", 60,-1.0,2.0);
    
    tuple->PostPreS_TriggerType = dir.make<TH1F>("PostPreS_TriggerType", ";;Events / category", 5, -0.5, 4.5);
    tuple->PostPreS_TriggerType->GetXaxis()->SetBinLabel(1,"Neither Muon nor MET triggered");
    tuple->PostPreS_TriggerType->GetXaxis()->SetBinLabel(2,"Muon triggered");
    tuple->PostPreS_TriggerType->GetXaxis()->SetBinLabel(3,"MET triggered");
    tuple->PostPreS_TriggerType->GetXaxis()->SetBinLabel(4,"Muon OR MET triggered");
    tuple->PostPreS_TriggerType->GetXaxis()->SetBinLabel(5,"Muon AND MET triggered");

    tuple->PostPreS_RecoHSCParticleType = dir.make<TH1F>("PostPreS_RecoHSCParticleType", ";;Tracks / category", 6, -0.5, 5.5);
    tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(1,"globalMuon");
    tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(2,"trackerMuon");
    tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(3,"matchedStandAloneMuon");
    tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(4,"standAloneMuon");
    tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(5,"innerTrack");
    tuple->PostPreS_RecoHSCParticleType->GetXaxis()->SetBinLabel(6,"unknown");

    tuple->PostPreS_PfType = dir.make<TH1F>("PostPreS_PfType", ";;Tracks / category", 9, -0.5, 8.5);
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(1,"AllTracks");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(2,"PFtracks");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(3,"isElectron");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(4,"isMuon");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(5,"isPhoton");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(6,"isChHadron");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(7,"isNeutHadron");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(8,"isUndefined");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(9,"notPFtrack");

    tuple->PostPreS_PfTypeVsIas = dir.make<TH2F>("PostPreS_PfTypeVsIas", ";;G_{i}^{strips}", 9, -0.5, 8.5,20,0.,1.);
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(1,"AllTracks");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(2,"PFtracks");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(3,"isElectron");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(4,"isMuon");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(5,"isPhoton");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(6,"isChHadron");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(7,"isNeutHadron");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(8,"isUndefined");
    tuple->PostPreS_PfType->GetXaxis()->SetBinLabel(9,"notPFtrack");

    tuple->PostPreS_MassT = dir.make<TH1F>("PostPreS_MassT", ";m_{T} (GeV);Tracks / 5 GeV", 50, 0.0, 250.0);
    tuple->PostPreS_MassT_highMassT = dir.make<TH1F>("PostPreS_MassT_highMassT", ";m_{T} (GeV);Tracks / 10 GeV", 250, 0.0, 2500.0);
    tuple->PostPreS_MassTVsIas = dir.make<TH2F>("PostPreS_MassTVsIas", ";m_{T} (GeV);G_{i}^{strips}",50, 0.0, 250.0, 20, 0., 1.);

    tuple->PostPreS_MiniRelIsoAll = dir.make<TH1F>("PostPreS_MiniRelIsoAll", ";MiniRelIsoAll;Tracks / 0.001", 100, 0.0, 0.1);
    tuple->PostPreS_MiniRelIsoAllVsIas =  dir.make<TH2F>("PostPreS_MiniRelIsoAllVsIas",";MiniRelIsoAll;G_{i}^{strips}", 100, 0.0, 0.1, 10, 0.,1.);
    tuple->PostPreS_MiniRelIsoChg = dir.make<TH1F>("PostPreS_MiniRelIsoChg", ";MiniRelIsoChg;Tracks / 0.01",  150, 0.0, 1.5);
    tuple->PostPreS_MiniTkIso = dir.make<TH1F>("PostPreS_MiniTkIso", ";MiniTkIso (GeV);Tracks / 0.05", 100, 0.0, 50.);
    tuple->PostPreS_MiniRelTkIso = dir.make<TH1F>("PostPreS_MiniRelTkIso", ";MiniRelTkIso;Tracks / 0.01", 150, 0.0, 1.5);

    tuple->PostPreS_RecoPFMET = dir.make<TH1F>("PostPreS_RecoPFMET", ";RecoPFMET (GeV);Tracks / 10 GeV",  200, 0.0, 2000.0);
    tuple->PostPreS_RecoPFHT = dir.make<TH1F>("PostPreS_RecoPFHT", ";RecoPFHT (GeV);Tracks / 10 GeV",  200, 0.0, 2000.0);
    tuple->PostPreS_CaloJetsNum = dir.make<TH1F>("PostPreS_CaloJetsNum", ";Number of calo jets;Tracks / 1",  16, -0.5, 15.5);
    tuple->PostPreS_Chi2oNdof = dir.make<TH1F>("PostPreS_Chi2oNdof", ";#chi^{2}/N_{dof};Tracks / 1", 20, 0, 20);
    tuple->PostPreS_Chi2oNdofVsIas = dir.make<TH2F>("PostPreS_Chi2oNdofVsIas", ";#chi^{2}/Ndof;G_{i}^{strips}",20, 0, 20,10,0.,1.);
    tuple->PostPreS_Qual = dir.make<TH1F>("PostPreS_Qual", ";;Tracks / category", 2, -0.5, 1.5);
    tuple->PostPreS_Qual->GetXaxis()->SetBinLabel(1,"Not-HighPurity");
    tuple->PostPreS_Qual->GetXaxis()->SetBinLabel(2,"HighPurity");

    tuple->PostPreS_TNOH_PUA = dir.make<TH1F>("PostPreS_TNOH_PUA", ";Number of hits (low PU);Tracks / 1", 40, -0.5, 39.5);
    tuple->PostPreS_TNOH_PUB = dir.make<TH1F>("PostPreS_TNOH_PUB", ";Number of hits (mid PU);Tracks / 1", 40, -0.5, 39.5);
    tuple->PostPreS_TNOH_PUC = dir.make<TH1F>("PostPreS_TNOH_PUC", ";Number of hits (high PU);Tracks / 1", 40, -0.5, 39.5);
    tuple->PostPreS_TNOHFraction = dir.make<TH1F>("PostPreS_TNOHFraction", ";Number of valid hit fraction;Tracks / 0.05", 20, 0, 1);
    tuple->PostPreS_TNOHFractionVsIas = dir.make<TH2F>("PostPreS_TNOHFractionVsIas",";TNOHFraction;Ias;Tracks",20, 0., 1.,10,0.,1.);
    tuple->PostPreS_TNOPH = dir.make<TH1F>( "PostPreS_TNOPH", ";Number of pixel hits;Tracks / 1", 8, -0.5, 7.5);
    tuple->PostPreS_RatioCleanAndAllStripsClu = dir.make<TH1F>("PostPreS_RatioCleanAndAllStripsClu",";Clean / all strips clu;Tracks / 0.05",20,-0.05,1.05);
    tuple->PostPreS_RatioCleanAndAllStripsCluVsIas = dir.make<TH2F>("PostPreS_RatioCleanAndAllStripsCluVsIas",";Clean / all strips clu;G_{i}^{strips};Tracks / 0.05",20,-0.05,1.05,20,0.,1.);
    tuple->PostPreS_RatioCleanAndAllPixelClu = dir.make<TH1F>("PostPreS_RatioCleanAndAllPixelClu",";Clean / all pixel clu;Tracks / 0.05",20,-0.05,1.05);
    
    tuple->PostPreS_TNOPHVsIas = dir.make<TH2F>("PostPreS_TNOPHVsIas", ";_TNOPH;G_{i}^{strips}", 8,-0.5, 7.5, 20, 0., 1.);
    tuple->PostPreS_TNOHFractionTillLast = dir.make<TH1F>("PostPreS_TNOHFractionTillLast", ";TNOHFractionTillLast;Tracks / 0.05", 20, 0., 1.);
    tuple->PostPreS_TNOMHTillLast = dir.make<TH1F>("PostPreS_TNOMHTillLast", ";TNOMHTillLast;Tracks / 1", 20, -0.5, 19.5);
    tuple->PostPreS_Eta = dir.make<TH1F>("PostPreS_Eta", ";#eta;Tracks / 0.05", 52, -2.6, 2.6);
    tuple->PostPreS_EtaVsIas =  dir.make<TH2F>("PostPreS_EtaVsIas", ";#eta;G_{i}^{strips};Tracks / 0.05", 52, -2.6, 2.6, 20,0.,1.);
    tuple->PostPreS_TNOM = dir.make<TH1F>("PostPreS_TNOM", ";Number of measurement;Tracks / 1", 40, -0.5, 39.5);
    tuple->PostPreS_TNOMVsIas = dir.make<TH2F>("PostPreS_TNOMVsIas", ";Number of measurement;G_{i}^{strips}",  40, -0.5, 39.5, 20, 0., 1.);
    tuple->PostPreS_EtaVsNBH = dir.make<TH2F>("PostPreS_EtaVsNBH", ";#eta;Number of bad hits;Tracks / bin", 60, -3., 3., 24, 0., 24.);
    tuple->PostPreS_TNOM_PUA = dir.make<TH1F>("PostPreS_TNOM_PUA",  ";Number of measurement (low PU);Tracks / 1", 40, -0.5, 39.5);
    tuple->PostPreS_TNOM_PUB = dir.make<TH1F>("PostPreS_TNOM_PUB",  ";Number of measurement (mid PU);Tracks / 1",  40, -0.5, 39.5);
    tuple->PostPreS_TNOM_PUC = dir.make<TH1F>("PostPreS_TNOM_PUB",  ";Number of measurement (high PU);Tracks / 1",  40, -0.5, 39.5);
    tuple->PostPreS_nDof = dir.make<TH1F>("PostPreS_nDof", ";nDof;Tracks / 1",  40, -0.5, 39.5);
    tuple->PostPreS_TOFError = dir.make<TH1F>("PostPreS_TOFError", ";TOFError;Tracks / 0.01", 25, 0, 0.25);
    tuple->PostPreS_PtErrOverPt = dir.make<TH1F>("PostPreS_PtErrOverPt", ";#sigma_{p_{T}}/p_{T};Tracks / 0.025", 40, 0, 1);
    tuple->PostPreS_PtErrOverPtVsIas =  dir.make<TH2F>("PostPreS_PtErrOverPtVsIas", ";#sigma_{p_{T}}/p_{T};G_{i}^{strips};Tracks / bin", 40, 0, 1, 20, 0.,1.);
    tuple->PostPreS_PtErrOverPt2VsIas =  dir.make<TH2F>("PostPreS_PtErrOverPt2VsIas", ";#sigma_{p_{T}}/p_{T}^{2};G_{i}^{strips};Tracks / bin", 40, 0., 0.003, 20, 0.,1.);
    tuple->PostPreS_PtErrOverPt2 = dir.make<TH1F>("PostPreS_PtErrOverPt2", ";#sigma_{p_{T}}/p_{T}^{2};Tracks / bin", 40, 0, 0.003);
    tuple->PostPreS_Pt = dir.make<TH1F>("PostPreS_Pt", ";p_{T} (GeV);Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->PostPreS_Pt_lowPt = dir.make<TH1F>("PostPreS_Pt_lowPt", ";p_{T} (GeV);Tracks / 10 GeV", 50, 0., 500.);
    tuple->PostPreS_PtVsIas = dir.make<TH2F>("PostPreS_PtVsIas",";p_{T};G_{i}^{strips};Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 20, 0., 1.);
    tuple->PostPreS_Ias_CR = dir.make<TH1F>("PostPreS_Ias_CR", ";G_{i}^{strips};Tracks / 0.1", 10, 0, dEdxS_UpLim);
    tuple->PostPreS_Pt_lowPt_CR = dir.make<TH1F>("PostPreS_Pt_lowPt_CR", ";p_{T} (GeV);Tracks / 10 GeV", 50, 0., 500.);
    tuple->PostPreS_Ias_CR_veryLowPt = dir.make<TH1F>("PostPreS_Ias_CR_veryLowPt", ";G_{i}^{strips};Tracks / 0.1", 10, 0, dEdxS_UpLim);
    tuple->PostPreS_P_CR_veryLowPt = dir.make<TH1F>("PostPreS_P_CR_veryLowPt", ";p (GeV);Tracks / 10 GeV", 50, 0., 50);
    tuple->PostPreS_Pt_lowPt = dir.make<TH1F>("PostPreS_Pt_lowPt", ";p_{T} (GeV);Tracks / 10 GeV", 50, 0., 500.);
    tuple->PostPreS_Ias_NoEventWeight = dir.make<TH1F>("PostPreS_Ias_NoEventWeight", ";G_{i}^{strips} (NoEventWeight);Tracks / 0.1", 10, 0, dEdxS_UpLim);
    tuple->PostPreS_FiStrips_NoEventWeight = dir.make<TH1F>("PostPreS_FiStrips_NoEventWeight", ";F_{i}^{strips} (NoEventWeight);Tracks / 0.05", 20, 0., 1.);
    tuple->PostPreS_Ih = dir.make<TH1F>("PostPreS_Ih", ";I_{h} (MeV/cm)", 200, 0, dEdxM_UpLim);
    tuple->PostPreS_IhVsIas = dir.make<TH2F>("PostPreS_IhVsIas",";I_{h} (MeV/cm);G_{i}^{strips}",200, 0, dEdxM_UpLim, 20, 0.,1.);
    tuple->PostPreS_Ih_NoEventWeight = dir.make<TH1F>("PostPreS_Ih_NoEventWeight", ";I_{h} (NoEventWeight) (MeV/cm);Tracks / bin", 200, 0, dEdxM_UpLim);
    tuple->PostPreS_MTOF = dir.make<TH1F>("PostPreS_MTOF", ";MTOF;Tracks / bin", 50, -2, 5);
    tuple->PostPreS_TIsol = dir.make<TH1F>("PostPreS_TIsol", ";#Sigma_{R<0.3} p_{T} - p_{T,cand} (GeV);Tracks / 4 GeV", 25, 0, 100);
    tuple->PostPreS_TIsolVsIas = dir.make<TH2F>("PostPreS_TIsolVsIas","TIsol;Cluster G_{i}^{strips};Tracks / bin",25, 0, 100, 10, 0., 1.);
    tuple->PostPreS_EoP = dir.make<TH1F>("PostPreS_EoP", ";PF calo energy / momentum; Tracks / 0.06", 25, 0, 1.5);
    tuple->PostPreS_ECalEoP = dir.make<TH1F>("PostPreS_ECalEoP", ";PF ECal energy / momentum; Tracks / 0.06", 25, 0, 1.5);
    tuple->PostPreS_HCalEoP = dir.make<TH1F>("PostPreS_HCalEoP", ";PF HCal energy / momentum; Tracks / 0.06", 25, 0, 1.5);
    tuple->PostPreS_EoPVsIas = dir.make<TH2F>("PostPreS_EoPVsIas","PostPreS_EoPVsIas;PF E/p;G_{i}^{strips}",25, 0, 1.5, 20, 0.,1.);
    tuple->PostPreS_SumpTOverpT = dir.make<TH1F>("PostPreS_SumpTOverpT", ";SumpTOverpT", 80, 0.0, 2.0);
    tuple->PostPreS_SumpTOverpTVsIas = dir.make<TH2F>("PostPreS_SumpTOverpTVsIas",";SumpTOverpT;G_{i}^{strips}", 80, 0.0, 2.0, 20, 0.,1.);
    tuple->PostPreS_LastHitDXY = dir.make<TH1F>("PostPreS_LastHitDXY", ";LastHitDXY;Tracks / bin", 75, 0, 150);
    tuple->PostPreS_LastHitDXYVsEta  = dir.make<TH2F>("PostPreS_LastHitDXYVsEta","PostPreS_LastHitDXYVsEta;LastHitDXY;#eta", 75, 0, 150, 20, 0.,1.);
    tuple->PostPreS_LastHitD3D = dir.make<TH1F>("PostPreS_LastHitD3D", ";LastHitD3D", 175, 0, 350);
    tuple->PostPreS_LastHitD3DVsEta = dir.make<TH2F>("PostPreS_LastHitD3DVsEta","PostPreS_LastHitD3DVsEta;LastHitD3D;#eta", 175, 0, 350, 20, 0.,1.);
    tuple->PostPreS_P = dir.make<TH1F>("PostPreS_P", ";Momentum (GeV);Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->PostPreS_dR_NVTrack = dir.make<TH1F>("PostPreS_dR_NVTrack", ";dR_NVTrack;Tracks / 0.025", 40, 0, 1);
    tuple->PostPreS_MatchedStations = dir.make<TH1F>("PostPreS_MatchedStations", ";MatchedStations",8, -0.5, 7.5);
    tuple->PostPreS_InnerInvPtDiff = dir.make<TH1F>("PostPreS_InnerInvPtDiff", ";InnerInvPtDiff", 120, -4, 4);
    tuple->PostPreS_Phi = dir.make<TH1F>("PostPreS_Phi", ";#phi;Tracks / bin", 50, -3.14, 3.14);
    tuple->PostPreS_TimeAtIP = dir.make<TH1F>("PostPreS_TimeAtIP", ";TimeAtIP;Tracks / bin", 50, -100, 100);
    tuple->PostPreS_OpenAngle = dir.make<TH1F>("PostPreS_OpenAngle", ";OpenAngle;Tracks / bin", 50, -0.3, 3.15);
    tuple->PostPreS_OpenAngle_Cosmic = dir.make<TH1F>("PostPreS_OpenAngle_Cosmic", ";OpenAngle_Cosmic;Tracks / bin", 50, -0.3, 3.15);

    tuple->PostPreS_NVertex = dir.make<TH1F>("PostPreS_NVertex", ";N_{vertex};Tracks / bin", 50, -0.5, 49.5);
    tuple->PostPreS_NVertex_NoEventWeight = dir.make<TH1F>("PostPreS_NVertex_NoEventWeight", ";N_{vertex} (NoEventWeight);Tracks / bin", 50, -0.5, 49.5);
    tuple->PostPreS_PV = dir.make<TH1F>("PostPreS_PV", ";PV;Tracks / bin", 60, -0.5, 59.5);
    tuple->PostPreS_PV_NoEventWeight = dir.make<TH1F>("PostPreS_PV_NoEventWeight", ";PV_NoEventWeight;Tracks / bin", 60, 0, 60);
    tuple->PostPreS_NOMoNOH = dir.make<TH1F>("PostPreS_NOMoNOH", ";Num of measurment / num of hits;Tracks / bin",10,0.,1.0);
    tuple->PostPreS_NOMoNOHvsPV = dir.make<TProfile>("PostPreS_NOMoNOHvsPV", ";NOMoNOHvsPV", 60, 0, 60);
    tuple->PostPreS_Dz = dir.make<TH1F>( "PostPreS_Dz", ";d_{z} (cm);Tracks / 0.003 cm", 200, -0.3, 0.3);
    tuple->PostPreS_DzVsIas = dir.make<TH2F>( "PostPreS_DzVsIas", ";d_{z} (cm);G_{i}^{strips};Tracks", 200, -0.3, 0.3, 20, 0.0, 1.0);
    tuple->PostPreS_DzVsGenID = dir.make<TH2F>( "PostPreS_DzVsGenID", ";d_{z} (cm);GenID;Tracks", 200, -0.3, 0.3, 4000, 0.0, 4000.0);
    tuple->PostPreS_Dxy = dir.make<TH1F>("PostPreS_Dxy", ";d_{xy} (cm);Tracks / 0.001 cm", 200, -0.1, 0.1);
    tuple->PostPreS_DxyVsIas = dir.make<TH2F>("PostPreS_DxyVsIas", ";d_{xy} (cm);G_{i}^{strips};Tracks", 200, -0.1, 0.1, 20, 0.0, 1.0);
    tuple->PostPreS_DxyVsGenID = dir.make<TH2F>("PostPreS_DxyVsGenID", ";d_{xy} (cm);GenID;Tracks", 200, -0.1, 0.1, 4000, 0.0, 4000.0);

    tuple->PostPreS_SegSep = dir.make<TH1F>("PostPreS_SegSep", ";SegSep;Tracks / bin", 50, 0., 2.5);
    tuple->PostPreS_SegMinEtaSep = dir.make<TH1F>("PostPreS_SegMinEtaSep", ";SegMinEtaSep;Tracks / bin", 50, -1., 1.);
    tuple->PostPreS_SegMinPhiSep = dir.make<TH1F>("PostPreS_SegMinPhiSep", ";SegMinPhiSep;Tracks / bin", 50, -3.3, 3.3);
    tuple->PostPreS_SegMinEtaSep_FailDz = dir.make<TH1F>("PostPreS_SegMinEtaSep_FailDz", ";SegMinEtaSep_FailDz;Tracks / bin", 50, -1., 1.);
    tuple->PostPreS_SegMinEtaSep_PassDz = dir.make<TH1F>("PostPreS_SegMinEtaSep_PassDz", ";SegMinEtaSep_PassDz;Tracks / bin", 50, -1., 1.);
    tuple->PostPreS_Dz_FailSep = dir.make<TH1F>("PostPreS_Dz_FailSep", ";Dz_FailSep;Tracks / bin", 50, -150, 150);

    tuple->PostPreS_Dxy_Cosmic = dir.make<TH1F>("PostPreS_Dxy_Cosmic", ";Dxy_Cosmic;Tracks / bin", 150, -IPbound, IPbound);
    tuple->PostPreS_Dz_Cosmic = dir.make<TH1F>("PostPreS_Dz_Cosmic", ";Dz_Cosmic;Tracks / bin", 150, -IPbound, IPbound);
    tuple->PostPreS_Dz_CSC = dir.make<TH1F>("PostPreS_Dz_CSC", ";Dz_CSC;Tracks / bin", 150, -IPbound, IPbound);
    tuple->PostPreS_Dz_DT = dir.make<TH1F>("PostPreS_Dz_DT", ";Dz_DT;Tracks / bin", 150, -IPbound, IPbound);
    tuple->PostPreS_Pt_FailDz = dir.make<TH1F>("PostPreS_Pt_FailDz", ";Pt_FailDz;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->PostPreS_Pt_FailDz_DT = dir.make<TH1F>("PostPreS_Pt_FailDz_DT", ";Pt_FailDz_DT;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->PostPreS_Pt_FailDz_CSC = dir.make<TH1F>("PostPreS_Pt_FailDz_CSC", ";Pt_FailDz_CSC;Tracks / 100 GeV", 40, 0, PtHistoUpperBound);
    tuple->PostPreS_TOF_FailDz = dir.make<TH1F>("PostPreS_TOF_FailDz", ";TOF_FailDz;Tracks / 0.05", 120, -1, 5);
    tuple->PostPreS_TOF_FailDz_DT = dir.make<TH1F>("PostPreS_TOF_FailDz_DT", ";TOF_FailDz_DT;Tracks / 0.05", 120, -1, 5);
    tuple->PostPreS_TOF_FailDz_CSC = dir.make<TH1F>("PostPreS_TOF_FailDz_CSC", ";TOF_FailDz_CSC;Tracks / 0.05", 120, -1, 5);
    tuple->PostPreS_PtErrOverPtVsPtErrOverPt2 = dir.make<TH2F>("PostPreS_PtErrOverPtVsPtErrOverPt2",  ";#sigma_{p_{T}}/p_{T};p_{T}^{2};Tracks / bin",  40, 0., 1., 40, 0., 0.003);
    tuple->PostPreS_PtErrOverPtVsPt = dir.make<TH2F>("PostPreS_PtErrOverPtVsPt",  ";#sigma_{p_{T}}/p_{T};p_{T};Tracks / bin",  40, 0., 1., 40, 0., 4000);
    tuple->PostPreS_GenPtVsRecoPt = dir.make<TH2F>("PostPreS_GenPtVsRecoPt", ";GenPt;RecoPt;Tracks / 80 GeV", 40, 0., PtHistoUpperBound, 40, 0, PtHistoUpperBound);

    tuple->PostPreS_ProbQ = dir.make<TH1F>("PostPreS_ProbQ", ";F_{i}^{pixels}", 100, 0., 1.);
    tuple->PostPreS_ProbQVsIas = dir.make<TH2F>("PostPreS_ProbQVsIas", ";F_{i}^{pixels};G_{i}^{strips}", 100, 0., 1., 100, 0., 1.);
    tuple->PostPreS_IhVsProbQNoL1VsIas = dir.make<TH3F>("PostPreS_IhVsProbQVsIas", ";I_{h} (MeV/cm);F_{i}^{pixels};G_{i}^{strips}",200, 0, dEdxM_UpLim, 100, 0.0, 1.0, 100, 0.0, 1.0);
    tuple->PostPreS_MomentumVsProbQNoL1VsIas = dir.make<TH3F>("PostPreS_MomentumVsProbQVsIas", ";p (GeV);F_{i}^{pixels};G_{i}^{strips}", 40, 0., PtHistoUpperBound, 100, 0.0, 1.0, 100, 0.0, 1.0);
    tuple->PostPreS_ProbXY = dir.make<TH1F>("PostPreS_ProbXY", ";Prob_{XY,pixelAV} (pixels);Tracks / bin", 100, 0, 1);
    tuple->PostPreS_ProbXY_highIas = dir.make<TH1F>("PostPreS_ProbXY_highIas", ";ProbXY (G_{i}^{strips} > 0.6);Tracks / bin", 100, 0, 1);
    tuple->PostPreS_ProbXYVsIas = dir.make<TH2F>("PostPreS_ProbXYVsIas", ";ProbXY;G_{i}^{strips};Tracks / bin",  100, 0, 1, 10, 0., 1.);
    tuple->PostPreS_ProbXYVsIas_highIas = dir.make<TH2F>("PostPreS_ProbXYVsIas_highIas", ";ProbXY (G_{i}^{strips} > 0.6);G_{i}^{strips} (G_{i}^{strips} > 0.6);Tracks / bin",  100, 0, 1, 10, 0., 1.);
    tuple->PostPreS_ProbXYVsProbQ = dir.make<TH2F>("PostPreS_ProbXYVsProbQ", ";Prob_{XY,pixelAV} (pixels);F_{i}^{pixels};Tracks / bin",  100, 0., 1., 10, 0., 1.);
    tuple->PostPreS_ProbXYVsProbQ_highIas = dir.make<TH2F>("PostPreS_ProbXYVsProbQ_highIas", ";ProbXY (G_{i}^{strips} > 0.6);ProbQ (G_{i}^{strips} > 0.6);Tracks / bin",  100, 0., 1., 10, 0., 1.);
    
    tuple->PostPreS_ProbQNoL1_CR = dir.make<TH1F>("PostPreS_ProbQNoL1_CR", ";F_{i}^{pixels};Tracks / bin", 20, 0., 1.);
    
    tuple->PostPreS_ProbQNoL1VsIas_CR = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_Pileup_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_Pileup_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_Pileup_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_Pileup_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_ProbQNoL1_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_ProbQNoL1_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_ProbQNoL1_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_ProbQNoL1_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_Ias_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_Ias_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_Ias_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_Ias_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_Pt_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_Pt_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_CR_Pt_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_CR_Pt_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    
    tuple->PostPreS_ProbQNoL1 = dir.make<TH1F>("PostPreS_ProbQNoL1", ";F_{i}^{pixels};Tracks / bin", 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsFiStrips = dir.make<TH2F>("PostPreS_ProbQNoL1VsFiStrips", ";F_{i}^{pixels};F_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_Pileup_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_Pileup_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_Pileup_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_Pileup_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_ProbQNoL1_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_ProbQNoL1_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_ProbQNoL1_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_ProbQNoL1_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_Pt_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_Pt_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_Pt_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_Pt_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_Ias_up = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_Ias_up", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_ProbQNoL1VsIas_Ias_down = dir.make<TH2F>("PostPreS_ProbQNoL1VsIas_Ias_down", ";F_{i}^{pixels};G_{i}^{strips};Tracks / bin",20, 0., 1., 20, 0., 1.);
    tuple->PostPreS_TriggerMuon50VsBeta = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaA = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaA", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaA_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaA_BetaUp", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaA_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaA_BetaDown", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaB = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaB", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaB_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaB_BetaUp", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaB_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaB_BetaDown", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaC = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaC", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaC_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaC_BetaUp", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuon50VsBeta_EtaC_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuon50VsBeta_EtaC_BetaDown", ";Muon50 triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    
    tuple->PostPreS_TriggerMuonAllVsBeta = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaA = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaA", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaA_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaA_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaA_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaA_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaB = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaB", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaB_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaB_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaB_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaB_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaC = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaC", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaC_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaC_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaC_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaC_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);

    tuple->PostPreS_TriggerMuon50VsPt = dir.make<TH2F>("PostPreS_TriggerMuon50VsPt", ";Muon50 triggered;Track p_{T};Tracks / bin",2,-.5,1.5,50,-0.05,1000.05);
    tuple->PostPreS_TriggerMuonAllVsPt = dir.make<TH2F>("PostPreS_TriggerMuonAllVsPt", ";MuonAll triggered;Track p_{T};Tracks / bin",2,-.5,1.5,50,-0.05,1000.05);
    
    tuple->PostPreS_TriggerMuonAllVsBeta = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaA = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaA", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaA_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaA_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaA_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaA_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaB = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaB", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaB_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaB_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaB_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaB_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaC = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaC", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaC_BetaUp = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaC_BetaUp", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMuonAllVsBeta_EtaC_BetaDown = dir.make<TH2F>("PostPreS_TriggerMuonAllVsBeta_EtaC_BetaDown", ";MuonAll triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    
    tuple->PostPreS_TriggerMuonAllVsPt = dir.make<TH2F>("PostPreS_TriggerMuonAllVsPt", ";MuonAll triggered;Track p_{T};Tracks / bin",2,-.5,1.5,50,-0.05,1000.05);
    
    tuple->PostPreS_TriggerMETallVsBeta = dir.make<TH2F>("PostPreS_TriggerMETallVsBeta", ";OR of MET triggered;Gen #beta;Tracks / bin",2,-.5,1.5,20,0.,1.);
    tuple->PostPreS_TriggerMETallVsMet = dir.make<TH2F>("PostPreS_TriggerMETallVsMet", ";OR of MET triggered;MET (GeV);Tracks / bin",2,-.5,1.5,50,-0.05,2000.05);
    tuple->PostPreS_TriggerMETallVsHT = dir.make<TH2F>("PostPreS_TriggerMETallVsHT", ";OR of MET triggered;H_{T} (GeV);Tracks / bin",2,-.5,1.5,50,-0.05,2000.05);
    tuple->PostPreS_TriggerMETallVsMetOverHt = dir.make<TH2F>("PostPreS_TriggerMETallVsMetOverHt", ";OR of MET triggered;MET / H_{T};Tracks / bin",2,-.5,1.5,30,-0.05,2.95);
    
    tuple->PostS_TriggerMETallVsMetOverHt_Cand0 = dir.make<TH2F>("PostS_TriggerMETallVsMetOverHt_Cand0", ";OR of MET triggered;MET / H_{T};Tracks / bin",2,-.5,1.5,30,-0.05,2.95);
    tuple->PostS_TriggerMETallVsMetOverHt_Cand1 = dir.make<TH2F>("PostS_TriggerMETallVsMetOverHt_Cand1", ";OR of MET triggered;MET / H_{T};Tracks / bin",2,-.5,1.5,30,-0.05,2.95);
    tuple->PostS_TriggerMETallVsMetOverHt_Cand2 = dir.make<TH2F>("PostS_TriggerMETallVsMetOverHt_Cand2", ";OR of MET triggered;MET / H_{T};Tracks / bin",2,-.5,1.5,30,-0.05,2.95);
    
    
    tuple->PostPreS_MetOverHt = dir.make<TH1F>("PostPreS_MetOverHt", ";MET / H_{T};Tracks / bin",30,-0.05,2.95);
    tuple->PostS_MetOverHt_Cand0 = dir.make<TH1F>("PostS_MetOverHt_Cand0", ";MET / H_{T};Events / bin",30,-0.05,2.95);
    tuple->PostS_MetOverHt_Cand1 = dir.make<TH1F>("PostS_MetOverHt_Cand1", ";MET / H_{T};Events / bin",30,-0.05,2.95);
    tuple->PostS_MetOverHt_Cand2 = dir.make<TH1F>("PostS_MetOverHt_Cand2", ";MET / H_{T};Events / bin",30,-0.05,2.95);
    
    
    tuple->PostPreS_MetVsHT = dir.make<TH2F>("PostPreS_MetVsHT", ";MET (GeV);H_{T} (GeV);Tracks / bin",25,-0.05,2000.05,25,-0.05,2000.05);
    tuple->PostPreS_TriggerMETallVsMetVsHT = dir.make<TH3F>("PostPreS_TriggerMETallVsMetVsHT", ";OR of MET triggered;MET (GeV);H_{T} (GeV);Tracks / bin",2,-.5,1.5,50,-0.05,2000.05,50,-0.05,2000.05);

    tuple->PostPreS_ProbXYNoL1 = dir.make<TH1F>("PostPreS_ProbXYNoL1", ";ProbXYNoL1;Tracks / bin", 100, 0, 1);
    tuple->PostPreS_ProbXYNoL1_highIas = dir.make<TH1F>("PostPreS_ProbXYNoL1_highIas", ";ProbXYNoL1 for (G_{i}^{strips} > 0.6);Tracks / bin", 100, 0, 1);
    tuple->PostPreS_ProbXYNoL1VsIas = dir.make<TH2F>("PostPreS_ProbXYNoL1VsIas", ";ProbXYNoL1;G_{i}^{strips};Tracks / bin", 100, 0., 1., 20, 0.,1.);
    tuple->PostPreS_ProbXYNoL1VsIas_highIas = dir.make<TH2F>("PostPreS_ProbXYNoL1VsIas_highIas", ";ProbXYNoL1;G_{i}^{strips} (G_{i}^{strips} > 0.6);Tracks / bin", 100, 0., 1., 20, 0.,1.);
    tuple->PostPreS_ProbXYNoL1VsProbQNoL1  = dir.make<TH2F>("PostPreS_ProbXYNoL1VsProbQNoL1", ";ProbXYNoL1;F_{i}^{pixels};Tracks / bin", 100, 0., 1., 20,0.,1.);
    tuple->PostPreS_ProbXYNoL1VsProbQNoL1_highIas  = dir.make<TH2F>("PostPreS_ProbXYNoL1VsProbQNoL1_highIas", ";ProbXYNoL1 (G_{i}^{strips} > 0.6);ProbQNoL1 (G_{i}^{strips} > 0.6)", 100, 0., 1., 20,0.,1.);

    tuple->PostPreS_MassErr = dir.make<TH1F>("PostPreS_MassErr", ";MassErr/Mass;Tracks / bin", 50, 0., 10.);
    tuple->PostPreS_MassErrVsIas = dir.make<TH2F>("PostPreS_MassErrVsIas", ";MassErr/Mass;G_{i}^{strips};Tracks / bin",50, 0., 10.,10,0.,1.);

    tuple->PostPreS_EtaVsGenID = dir.make<TH2F>("PostPreS_EtaVsGenID", ";#eta;GenID",  50, -2.6, 2.6, 4000, 0.0, 4000.0);
    tuple->PostPreS_ProbQVsGenID = dir.make<TH2F>("PostPreS_ProbQVsGenID", ";F_{i}^{pixels};GenID", 100, 0.0, 1.0, 4000, 0.0, 4000.0);

    tuple->PostPreS_IasForStatus91 = dir.make<TH1F>("PostPreS_IasForStatus91", ";G_{i}^{strips} when status=91;Tracks / bin", 10, 0., 1.);
    tuple->PostPreS_IasForStatusNot91 = dir.make<TH1F>("PostPreS_IasForStatusNot91", ";G_{i}^{strips} when status!=91;Tracks / bin", 10, 0., 1.);

    tuple->PostPreS_ProbQVsGenEnviromentID = dir.make<TH2F>("PostPreS_ProbQVsGenEnviromentID", ";F_{i}^{pixels};GenEnviromentID",20, 0.0, 1.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_ProbXYVsGenID = dir.make<TH2F>("PostPreS_ProbXYVsGenID", ";ProbXY;GenID", 100, 0.0, 1.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_PtVsGenID = dir.make<TH2F>("PostPreS_PtVsGenID", ";p_{T} (GeV);GenID", 40, 0., PtHistoUpperBound, 4000, 0.0, 4000.0);
    tuple->PostPreS_EoPVsGenID = dir.make<TH2F>("PostPreS_EoPVsGenID", ";EoPVsGenID", 25, 0, 1.5, 4000, 0.0, 4000.0);
    tuple->PostPreS_IhVsGenID = dir.make<TH2F>("PostPreS_IhVsGenID", ";I_{h} (MeV/cm);GenID", 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
    tuple->PostPreS_IasVsGenID = dir.make<TH2F>("PostPreS_IasVsGenID", ";G_{i}^{strips};GenID", 10, 0., 1., 4000, 0.0, 4000.0);
    tuple->PostPreS_IasVsGenEnviromentID = dir.make<TH2F>("PostPreS_IasVsGenEnviromentID", ";G_{i}^{strips};GenEnviromentID", 10, 0., 1., 4000, 0.0, 4000.0);
    tuple->PostPreS_MassTVsGenID = dir.make<TH2F>("PostPreS_massTVsGenID", ";m_{T} (GeV);GenID", 50, 0.0, 250.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_MiniIsoChgVsGenID = dir.make<TH2F>("PostPreS_MiniIsoChgVsGenID", ";miniIsoChg;GenID", 20, 0.0, 0.1, 4000, 0.0, 4000.0);
    tuple->PostPreS_MiniIsoAllVsGenID = dir.make<TH2F>("PostPreS_MiniIsoAllVsGenID", ";miniIsoAll;GenID", 20, 0.0, 0.1, 4000, 0.0, 4000.0);
    tuple->PostPreS_MassVsGenID = dir.make<TH2F>("PostPreS_MassVsGenID", ";Mass (GeV);GenID",80,0.,4000.,4000, 0.0, 4000.0);

    tuple->PostPreS_EtaVsMomGenID = dir.make<TH2F>("PostPreS_EtaVsMomGenID", ";#eta;MomGenID",  50, -2.6, 2.6, 4000, 0.0, 4000.0);
    tuple->PostPreS_ProbQVsMomGenID = dir.make<TH2F>("PostPreS_ProbQVsMomGenID", ";F_{i}^{pixels};MomGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_ProbXYVsMomGenID = dir.make<TH2F>("PostPreS_ProbXYVsMomGenID", ";ProbXY;MomGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_PtVsMomGenID = dir.make<TH2F>("PostPreS_PtVsMomGenID", ";p_{T} (GeV);MomGenID", 40, 0., PtHistoUpperBound, 4000, 0.0, 4000.0);
    tuple->PostPreS_EoPVsMomGenID = dir.make<TH2F>("PostPreS_EoPVsMomGenID", ";EoPVsMomGenID", 25, 0, 1.5, 4000, 0.0, 4000.0);
    tuple->PostPreS_IhVsMomGenID = dir.make<TH2F>("PostPreS_IhVsMomGenID", ";I_{h} (MeV/cm);MomGenID", 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
    tuple->PostPreS_IasVsMomGenID = dir.make<TH2F>("PostPreS_IasVsMomGenID", ";G_{i}^{strips};MomGenID", 10, 0., 1., 4000, 0.0, 4000.0);
    tuple->PostPreS_MassTVsMomGenID = dir.make<TH2F>("PostPreS_massTVsMomGenID", ";m_{T} (GeV);MomGenID", 50, 0.0, 250.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_MiniIsoChgVsMomGenID = dir.make<TH2F>("PostPreS_MiniIsoChgVsMomGenID", ";miniIsoChg;MomGenID", 20, 0.0, 0.1, 4000, 0.0, 4000.0);
    tuple->PostPreS_MiniIsoAllVsMomGenID = dir.make<TH2F>("PostPreS_MiniIsoAllVsMomGenID", ";miniIsoAll;MomGenID", 20, 0.0, 0.1, 4000, 0.0, 4000.0);
    tuple->PostPreS_MassVsMomGenID = dir.make<TH2F>("PostPreS_MassVsMomGenID","PostPreS_MassVsMomGenID;Mass;MomGenID",80,0.,4000.,4000, 0.0, 4000.0);

    tuple->PostPreS_EtaVsSiblingGenID = dir.make<TH2F>("PostPreS_EtaVsSiblingGenID", ";#eta;SiblingGenID",  50, -2.6, 2.6, 4000, 0.0, 4000.0);
    tuple->PostPreS_ProbQVsSiblingGenID = dir.make<TH2F>("PostPreS_ProbQVsSiblingGenID", ";F_{i}^{pixels};SiblingGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_ProbXYVsSiblingGenID = dir.make<TH2F>("PostPreS_ProbXYVsSiblingGenID", ";ProbXY;SiblingGenID", 20, 0.0, 1.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_PtVsSiblingGenID = dir.make<TH2F>("PostPreS_PtVsSiblingGenID", ";p_{T} (GeV);SiblingGenID", 40, 0., PtHistoUpperBound, 4000, 0.0, 4000.0);
    tuple->PostPreS_EoPVsSiblingGenID = dir.make<TH2F>("PostPreS_EoPVsSiblingGenID", ";EoPVsSiblingGenID", 25, 0, 1.5, 4000, 0.0, 4000.0);
    tuple->PostPreS_IhVsSiblingGenID = dir.make<TH2F>("PostPreS_IhVsSiblingGenID", ";I_{h} (MeV/cm);SiblingGenID", 200, 0, dEdxM_UpLim, 4000, 0.0, 4000.0);
    tuple->PostPreS_IasVsSiblingGenID = dir.make<TH2F>("PostPreS_IasVsSiblingGenID", ";G_{i}^{strips};SiblingGenID", 10, 0., 1., 4000, 0.0, 4000.0);
    tuple->PostPreS_MassTVsSiblingGenID = dir.make<TH2F>("PostPreS_massTVsSiblingGenID", ";massT;SiblingGenID", 50, 0.0, 250.0, 4000, 0.0, 4000.0);
    tuple->PostPreS_MassVsSiblingGenID = dir.make<TH2F>("PostPreS_MassVsSiblingGenID", ";Mass;SiblingGenID",80, 0.0, 4000.0, 4000, 0.0, 4000.0);

    tuple->PostPreS_EtaVsGenAngle = dir.make<TH2F>("PostPreS_EtaVsGenAngle", ";#eta;GenAngle",  50, -2.6, 2.6, 100, 0.0, 1.0);
    tuple->PostPreS_ProbQVsGenAngle = dir.make<TH2F>("PostPreS_ProbQVsGenAngle", ";F_{i}^{pixels};GenAngle", 20, 0.0, 1.0, 100, 0.0,1.0);
    tuple->PostPreS_ProbXYVsGenAngle = dir.make<TH2F>("PostPreS_ProbXYVsGenAngle", ";ProbXY;GenAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
    tuple->PostPreS_PtVsGenAngle = dir.make<TH2F>("PostPreS_PtVsGenAngle", ";p_{T} (GeV);GenAngle", 40, 0., PtHistoUpperBound, 100, 0.0, 1.0);
    tuple->PostPreS_EoPVsGenAngle = dir.make<TH2F>("PostPreS_EoPVsGenAngle", ";EoPVsGenAngle", 25, 0, 1.5, 100, 0.0, 1.0);
    tuple->PostPreS_IhVsGenAngle = dir.make<TH2F>("PostPreS_IhVsGenAngle", ";I_{h} (MeV/cm);GenAngle", 200, 0, dEdxM_UpLim, 100, 0.0, 1.0);
    tuple->PostPreS_IasVsGenAngle = dir.make<TH2F>("PostPreS_IasVsGenAngle", ";G_{i}^{strips};GenAngle", 10, 0., 1., 100, 0.0, 1.0);
    tuple->PostPreS_MassTVsGenAngle = dir.make<TH2F>("PostPreS_massTVsGenAngle", ";m_{T} (GeV);GenAngle", 50, 0.0, 250.0, 100, 0.0, 1.0);
    tuple->PostPreS_MiniIsoChgVsGenAngle = dir.make<TH2F>("PostPreS_MiniIsoChgVsGenAngle", ";miniIsoChg;GenAngle", 20, 0.0, 0.1, 100, 0.0, 1.0);
    tuple->PostPreS_MiniIsoAllVsGenAngle = dir.make<TH2F>("PostPreS_MiniIsoAllVsGenAngle", ";miniIsoAll;GenAngle", 20, 0.0, 0.1, 100, 0.0, 1.0);
    tuple->PostPreS_MassVsGenAngle = dir.make<TH2F>("PostPreS_MassVsGenAngle", ";Mass;GenAngle",80, 0.0,4000.0,  100, 0.0, 1.0);

    tuple->PostPreS_EtaVsGenMomAngle = dir.make<TH2F>("PostPreS_EtaVsGenMomAngle", ";#eta;GenMomAngle",  50, -2.6, 2.6, 100, 0.0, 1.0);
    tuple->PostPreS_ProbQVsGenMomAngle = dir.make<TH2F>("PostPreS_ProbQVsGenMomAngle", ";F_{i}^{pixels};GenMomAngle", 20, 0.0, 1.0, 100, 0.0,1.0);
    tuple->PostPreS_ProbXYVsGenMomAngle = dir.make<TH2F>("PostPreS_ProbXYVsGenMomAngle", ";ProbXY;GenMomAngle", 20, 0.0, 1.0, 100, 0.0, 1.0);
    tuple->PostPreS_PtVsGenMomAngle = dir.make<TH2F>("PostPreS_PtVsGenMomAngle", ";p_{T} (GeV);GenMomAngle", 40, 0., PtHistoUpperBound, 100, 0.0, 1.0);
    tuple->PostPreS_EoPVsGenMomAngle = dir.make<TH2F>("PostPreS_EoPVsGenMomAngle", ";EoPVsGenMomAngle", 25, 0, 1.5, 100, 0.0, 1.0);
    tuple->PostPreS_IhVsGenMomAngle = dir.make<TH2F>("PostPreS_IhVsGenMomAngle", ";I_{h} (MeV/cm);GenMomAngle", 100, 0, dEdxM_UpLim, 100, 0.0, 1.0);
    tuple->PostPreS_IasVsGenMomAngle = dir.make<TH2F>("PostPreS_IasVsGenMomAngle", ";G_{i}^{strips};GenMomAngle", 10, 0., 1., 100, 0.0, 1.0);
    tuple->PostPreS_MassTVsGenMomAngle = dir.make<TH2F>("PostPreS_massTVsGenMomAngle", ";m_{T} (GeV);GenMomAngle", 50, 0.0, 250.0, 100, 0.0, 1.0);
    tuple->PostPreS_MiniIsoChgVsGenMomAngle = dir.make<TH2F>("PostPreS_MiniIsoChgVsGenMomAngle", ";miniIsoChg;GenMomAngle", 20, 0.0, 0.1, 100, 0.0, 1.0);
    tuple->PostPreS_MiniIsoAllVsGenMomAngle = dir.make<TH2F>("PostPreS_MiniIsoAllVsGenMomAngle", ";miniIsoAll;GenMomAngle", 20, 0.0, 0.1, 100, 0.0, 1.0);
    tuple->PostPreS_MassVsGenMomAngle = dir.make<TH2F>("PostPreS_MassVsGenMomAngle", ";Mass;GenMomAngle",80,0.,4000., 100, 0.0, 1.0);

    tuple->PostPreS_ProbQVsIas = dir.make<TH2F>("PostPreS_ProbQVsIas", ";F_{i}^{pixels};G_{i}^{strips}", 20, 0., 1., 20, 0., 1.);

    tuple->PostPreS_EtaVsGenNumSibling = dir.make<TH2F>("PostPreS_EtaVsGenNumSibling", ";#eta;GenNumSibling",  50, -2.6, 2.6, 100, 0.0, 10.0);
    tuple->PostPreS_ProbQVsGenNumSibling = dir.make<TH2F>("PostPreS_ProbQVsGenNumSibling", ";F_{i}^{pixels};GenNumSibling", 20, 0., 1., 10, 0.0,10.);
    tuple->PostPreS_ProbXYVsGenNumSibling = dir.make<TH2F>("PostPreS_ProbXYVsGenNumSibling", ";ProbXY;GenNumSibling", 20, 0.0, 1.0, 10, 0.0, 10.0);
    tuple->PostPreS_PtVsGenNumSibling = dir.make<TH2F>("PostPreS_PtVsGenNumSibling", ";p_{T} (GeV);GenNumSibling", 40, 0., PtHistoUpperBound, 100, 0.0, 10.0);
    tuple->PostPreS_EoPVsGenNumSibling = dir.make<TH2F>("PostPreS_EoPVsGenNumSibling", ";EoPVsGenNumSibling", 25, 0, 1.5, 100, 0.0, 10.0);
    tuple->PostPreS_IhVsGenNumSibling = dir.make<TH2F>("PostPreS_IhVsGenNumSibling", ";I_{h} (MeV/cm);GenNumSibling", 200, 0, dEdxM_UpLim, 100, 0.0, 10.0);
    tuple->PostPreS_IasVsGenNumSibling = dir.make<TH2F>("PostPreS_IasVsGenNumSibling", ";G_{i}^{strips};GenNumSibling", 10, 0., 1., 100, 0.0, 10.0);
    tuple->PostPreS_MassTVsGenNumSibling = dir.make<TH2F>("PostPreS_massTVsGenNumSibling", ";m_{T} (GeV);GenNumSibling", 50, 0.0, 250.0, 100, 0.0, 10.0);
    tuple->PostPreS_MiniIsoChgVsGenNumSibling = dir.make<TH2F>("PostPreS_MiniIsoChgVsGenNumSibling", ";miniIsoChg;GenNumSibling", 20, 0.0, 0.1, 100, 0.0, 10.0);
    tuple->PostPreS_MiniIsoAllVsGenNumSibling = dir.make<TH2F>("PostPreS_MiniIsoAllVsGenNumSibling", ";miniIsoAll;GenNumSibling", 20, 0.0, 1.0, 100, 0.0, 10.0);

    tuple->PostPreS_EoPVsPfType = dir.make<TH2F>("PostPreS_EoPVsPfType", ";EoPVsPfType", 25, 0.0, 1.5, 9, 0.0, 9.0);
    tuple->PostPreS_Mass = dir.make<TH1F>("PostPreS_Mass","PostPreS_Mass;Mass (GeV);Tracks / 50 GeV", 80,0.,4000.);
    tuple->PostPreS_MassVsPfType = dir.make<TH2F>("PostPreS_MassVsPfType","PostPreS_MassVsPfType;Mass (GeV);PF ID",80,0.,4000.,9, 0.0, 9.0);
    tuple->PostPreS_MassVsPt = dir.make<TH2F>("PostPreS_MassVsPt", ";Mass (GeV);", 80,0.,4000.,80,0.,4000.);
    tuple->PostPreS_MassVsP = dir.make<TH2F>("PostPreS_MassVsP", ";Mass (GeV);", 80,0.,4000.,80,0.,4000.);
    tuple->PostPreS_MassVsTNOHFraction = dir.make<TH2F>("PostPreS_MassVsTNOHFraction", ";Mass (GeV);Number of valid hit fraction", 80,0.,4000.,20, 0., 1.);
    tuple->PostPreS_MassVsTNOPH = dir.make<TH2F>("PostPreS_MassVsTNOPH", ";Mass (GeV);TNOPH", 80,0.,4000.,8, -0.5, 7.5);
    tuple->PostPreS_MassVsTNOM = dir.make<TH2F>("PostPreS_MassVsTNOM", ";Mass (GeV);Number of measurements;", 80,0.,4000.,40, -0.5, 39.5);
    tuple->PostPreS_MassVsProbQNoL1 = dir.make<TH2F>("PostPreS_MassVsProbQNoL1", ";Mass (GeV);F_{i}^{pixels}", 80,0.,4000.,20,0.,1.);
    tuple->PostPreS_MassVsProbXYNoL1 = dir.make<TH2F>("PostPreS_MassVsProbXYNoL1", ";Mass (GeV);ProbXY", 80,0.,4000.,20,0.,1.);
    tuple->PostPreS_MassVsEoP = dir.make<TH2F>("PostPreS_MassVsEoP", ";Mass (GeV);E/p", 80,0.,4000.,30, 0, 0.3);
    tuple->PostPreS_MassVsSumpTOverpT = dir.make<TH2F>("PostPreS_MassVsSumpTOverpT", ";Mass (GeV);SumpTOverpT", 80,0.,4000.,80, 0, 2);
    tuple->PostPreS_MassVsPtErrOverPt = dir.make<TH2F>("PostPreS_MassVsPtErrOverPt", ";Mass (GeV);PtErrOverPt", 80,0.,4000.,40, 0, 1);
    tuple->PostPreS_MassVsTIsol = dir.make<TH2F>("PostPreS_MassVsTIsol", ";Mass (GeV);TIsol", 80,0.,4000., 25, 0, 100);
    tuple->PostPreS_MassVsIh = dir.make<TH2F>("PostPreS_MassVsIh", ";Mass (GeV);I_{h} (MeV/cm)", 80,0.,4000.,200, 0, dEdxM_UpLim);
    tuple->PostPreS_MassVsMassT = dir.make<TH2F>("PostPreS_MassVsMassT", ";Mass (GeV);m_{T} (GeV)", 80,0.,4000.,50, 0.0, 250.0);
    tuple->PostPreS_MassVsMiniRelIsoAll = dir.make<TH2F>("PostPreS_MassVsMiniRelIsoAll", ";Mass (GeV);MiniRelIsoAll", 80,0.,4000.,20, 0., 0.2);
    tuple->PostPreS_MassVsMassErr = dir.make<TH2F>("PostPreS_MassVsMassErr", ";Mass (GeV);MassErr", 80,0.,4000.,50, 0., 10.);

    // Maybe we dont need these anymore
    // Have to deal with this later on, should a boolean to have them or not
  //  tuple->PostPreS_IasAllIhVsLayer = dir.make<TH3F>("PostPreS_IasAllIhVsLayer", ";G_{i}^{strips};I_{h} (MeV/cm);LayerIndex (full tracker)", 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 35, 0.,35.);
  //  tuple->PostPreS_IasPixelIhVsLayer = dir.make<TH3F>("PostPreS_IasPixelIhVsLayer", ";G_{i}^{strips};I_{h} (MeV/cm);LayerIndex (pixels)", 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 10, 0.,10.);
  //  tuple->PostPreS_IasStripIhVsLayer = dir.make<TH3F>("PostPreS_IasStripIhVsLayer", ";G_{i}^{strips};I_{h} (MeV/cm);LayerIndex (strips)", 50, 0., dEdxS_UpLim, 200, 0., dEdxM_UpLim, 25, 0.,25.);

    tuple->PostPreS_CluProbQVsPixelLayer = dir.make<TH2F>("PostPreS_CluProbQVsPixelLayer", ";Cluster Prob_{Q} (pixels);Pixel Layer",20,0.,1.,4,0.,4.);
    tuple->PostPreS_CluProbXYVsPixelLayer = dir.make<TH2F>("PostPreS_CluProbXYVsPixelLayer", ";Cluster Prob_{XY,pixelAV} (pixels);Pixel Layer",100,0.,1.,4,0.,4.);
    tuple->PostPreS_CluSizeVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeVsPixelLayer", ";CluSize;Pixel Layer",10,0.,10.,4,0.,4.);
    tuple->PostPreS_CluSizeXVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeXVsPixelLayer", ";CluSizeX;Pixel Layer",10,0.,10.,4,0.,4.);
    tuple->PostPreS_CluSizeYVsPixelLayer = dir.make<TH2F>("PostPreS_CluSizeYVsPixelLayer", ";CluSizeY;Pixel Layer",10,0.,10.,4,0.,4.);
    tuple->PostPreS_CluSpecInCPEVsPixelLayer = dir.make<TH2F>("PostPreS_CluSpecInCPEVsPixelLayer", ";;Pixel Layer",4,-0.5,3.5,4,0.,4.);
    tuple->PostPreS_CluSpecInCPEVsPixelLayer->GetXaxis()->SetBinLabel(1,"isOnEdge");
    tuple->PostPreS_CluSpecInCPEVsPixelLayer->GetXaxis()->SetBinLabel(2,"hasBadPixels");
    tuple->PostPreS_CluSpecInCPEVsPixelLayer->GetXaxis()->SetBinLabel(3,"spansTwoROCs");

    tuple->PostPreS_CluProbQVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluProbQVsPixelLayer_highIas", ";Cluster Prob_{Q} (pixels);Layer",20,0.,1.,4,0.,4.);
    tuple->PostPreS_CluProbXYVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluProbXYVsPixelLayer_highIas", ";Cluster Prob_{XY,pixelAV} (pixels);Layer",100,0.,1.,4,0.,4.);
    tuple->PostPreS_CluSizeVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeVsPixelLayer_highIas", ";CluSize;Layer",10,0.,10.,4,0.,4.);
    tuple->PostPreS_CluSizeXVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeXVsPixelLayer_highIas", ";CluSizeX;Layer",10,0.,10.,4,0.,4.);
    tuple->PostPreS_CluSizeYVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSizeYVsPixelLayer_highIas", ";CluSizeY;Layer",10,0.,10.,4,0.,4.);
    tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas = dir.make<TH2F>("PostPreS_CluSpecInCPEVsPixelLayer_highIas", ";;Layer",4,-0.5,3.5,4,0.,4.);
    tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->GetXaxis()->SetBinLabel(1,"isOnEdge");
    tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->GetXaxis()->SetBinLabel(2,"hasBadPixels");
    tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->GetXaxis()->SetBinLabel(3,"spansTwoROCs");

    tuple->PostPreS_CluCotBetaVsPixelLayer_lowProbXY = dir.make<TH2F>("PostPreS_CluCotBetaVsPixelLayer_lowProbXY", ";CotBeta;Layer",200,-10.,10.,4,0.,4.);
    tuple->PostPreS_CluCotAlphaVsPixelLayer_lowProbXY = dir.make<TH2F>("PostPreS_CluCotAlphaVsPixelLayer_lowProbXY", ";CotAlpha;Layer",100,-1.,1.,4,0.,4.);
    tuple->PostPreS_CluCotBetaVsPixelLayer = dir.make<TH2F>("PostPreS_CluCotBetaVsPixelLayer", ";CotBeta;Layer",200,-10.,10.,4,0.,4.);
    tuple->PostPreS_CluCotAlphaVsPixelLayer = dir.make<TH2F>("PostPreS_CluCotAlphaVsPixelLayer", ";CotAlpha;Layer",100,-1.,1.,4,0.,4.);

    tuple->PostPreS_CluNormChargeVsStripLayer_lowBetaGamma = dir.make<TH2F>("PostPreS_CluNormChargeVsStripLayer_lowBetaGamma", ";CluNormCharge (e/um);Layer",600,0.,600.,24,-0.5,23.5);
    tuple->PostPreS_CluNormChargeVsStripLayer_higherBetaGamma = dir.make<TH2F>("PostPreS_CluNormChargeVsStripLayer_higherBetaGamma", ";CluNormCharge (e/um);Layer",600,0.,600.,24,-0.5,23.5);

    tuple->PostPreS_dRMinPfJet = dir.make<TH1F>("PostPreS_dRMinPfJet", ";dRMinPfJet",100,0.,5.0);
    tuple->PostPreS_closestPfJetMuonFraction = dir.make<TH1F>("PostPreS_closestPfJetMuonFraction",";closestPfJetMuonFraction; Tracks / 0.05",20,0.,1.);
    tuple->PostPreS_closestPfJetElectronFraction = dir.make<TH1F>("PostPreS_closestPfJetElectronFraction", ";closestPfJetElectronFraction; Tracks / 0.05",20,0.,1.);
    tuple->PostPreS_closestPfJetPhotonFraction = dir.make<TH1F>("PostPreS_closestPfJetPhotonFraction", ";closestPfJetPhotonFraction; Tracks / 0.05",20,0.,1.);

    tuple->PostPreS_closestPfJetMuonFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetMuonFractionVsIas",";closestPfJetMuonFraction;G_{i}^{strips}",20,0.,1.,20,0.,1.);
    tuple->PostPreS_closestPfJetElectronFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetElectronFractionVsIas", ";closestPfJetElectronFraction;G_{i}^{strips}",20,0.,1.,20,0.,1.);
    tuple->PostPreS_closestPfJetPhotonFractionVsIas = dir.make<TH2F>("PostPreS_closestPfJetPhotonFractionVsIas", ";closestPfJetPhotonFraction;G_{i}^{strips}",20,0.,1.,20,0.,1.);

    tuple->PostPreS_dRMinPfJetVsIas = dir.make<TH2F>("PostPreS_dRMinPfJetVsIas", ";dRMinPfJet;G_{i}^{strips}",100,0.,5.0,10,0.,1.);
    tuple->PostPreS_dRMinCaloJet = dir.make<TH1F>("PostPreS_dRMinCaloJet", ";dRMinCaloJet",100,0.,5.0);
    tuple->PostPreS_dPhiMinPfMet = dir.make<TH1F>("PostPreS_dPhiMinPfMet", ";dPhiMinPfMet",100,0.,3.2);

    tuple->PostPreS_dRMinCaloJetVsIas =  dir.make<TH2F>("PostPreS_dRMinCaloJetVsIas", ";dRMinCaloJet;G_{i}^{strips}",100,0.,5.0,10,0.,1.);
    tuple->PostPreS_dPhiMinPfMetVsIas =  dir.make<TH2F>("PostPreS_dPhiMinPfMetVsIas", ";dPhiMinPfMet;G_{i}^{strips}",100,0.,3.2,10,0.,1.);

    tuple->PostPreS_RecoPfJetsNum = dir.make<TH1F>("PostPreS_RecoPfJetsNum", ";Number of PF jets;Tracks / 1",  15, -0.5, 15.5);
    tuple->PostPreS_RecoPfHT = dir.make<TH1F>("PostPreS_RecoPfHT", ";PfHT",100,0.,2000.);
    tuple->PostPreS_GenBeta = dir.make<TH1F>("PostPreS_GenBeta", ";#beta;Gen candidate / 0.05", 20, 0., 1.);
  }
  tuple->PostS_HltMatchTrackLevel = dir.make<TH1F>("PostS_HltMatchTrackLevel", ";;Events / category", 4, 0.5, 4.5);
  tuple->PostS_HltMatchTrackLevel->GetXaxis()->SetBinLabel(1,"HLT + any muon match");
  tuple->PostS_HltMatchTrackLevel->GetXaxis()->SetBinLabel(2,"Bin 1 + best HSCP cand matches muon");
  tuple->PostS_HltMatchTrackLevel->GetXaxis()->SetBinLabel(3,"Bin 2 + Tight ID");
  tuple->PostS_HltMatchTrackLevel->GetXaxis()->SetBinLabel(4,"Bin 1 + best HSCP cand matches HLT muon");
  
  tuple->PostPreS_RecoPfMet = dir.make<TH1F>("PostPreS_RecoPfMet", ";PfMet",200,0.,2000.);
  tuple->PostPreS_RecoPfMetPhi = dir.make<TH1F>("PostPreS_RecoPfMetPhi", ";PfMetPhi",30,0.,3.2);

  //Initialize histograms for number of bins.  For everything but muon only PredBins=0 so no histograms created
  for (int i = 0; i < PredBins; i++) {
    char Suffix[1024];
    sprintf(Suffix, "_%i", i);
      Name.append(Suffix);
    tuple->BefPreS_Pt_Binned[std::to_string(i)] = dir.make<TH1F>("BefPreS_Pt_Binned", ";Pt_Binned", 40, 0, PtHistoUpperBound);
      Name.append(Suffix);
    tuple->BefPreS_TOF_Binned[std::to_string(i)] = dir.make<TH1F>("BefPreS_TOF_Binned", ";TOF_Binned", 150, -1, 5);
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
  tuple->PostS_CutIdVsBeta_postPtAndIas = dir.make<TH2F>("PostS_CutIdVsBeta_postPtAndIas", ";NCuts;#beta (p_{T} > p_{T,cut} and G_{i}^{strips} > I_{as,cut} )", NCuts, 0, NCuts, 20, 0, 1);
  tuple->PostS_CutIdVsBeta_postPtAndIasAndTOF = dir.make<TH2F>("PostS_CutIdVsBeta_postPtAndIasAndTOF", ";NCuts;#beta (p_{T} > p_{T,cut} and G_{i}^{strips} > I_{as,cut} and TOF > TOF_{cut} ", NCuts, 0, NCuts, 20, 0, 1);

  tuple->PostS_CutIdVsP = dir.make<TH2F>("PostS_CutIdVsP", ";NCuts;p (GeV)", NCuts, 0, NCuts, 40, 0, PtHistoUpperBound);
  tuple->PostS_CutIdVsPt = dir.make<TH2F>("PostS_CutIdVsPt", ";NCuts;p_{T} (GeV)", NCuts, 0, NCuts, 40, 0, PtHistoUpperBound);
  tuple->PostS_CutIdVsIas = dir.make<TH2F>("PostS_CutIdVsIas", ";NCuts;G_{i}^{strips}", NCuts, 0, NCuts, 10, 0., 1.);
  tuple->PostS_CutIdVsIh = dir.make<TH2F>("PostS_CutIdVsIh", ";NCuts;I_{h} (MeV/cm)", NCuts, 0, NCuts, 100, 0, dEdxM_UpLim);
    // not running this for a bit, they are not used currently, and being 3D histos they are quite big
//    tuple->PostS_CutIdVsPVsIas = dir.make<TH3F>("PostS_CutIdVsPVsIas", ";NCuts;p (GeV);G_{i}^{strips}", NCuts, 0, NCuts, 40, 0., PtHistoUpperBound, 10, 0., 1.);
//    tuple->PostS_CutIdVsPVsIh = dir.make<TH3F>("PostS_CutIdVsPVsIh", ";NCuts;P;I_{h} (MeV/cm)", NCuts, 0, NCuts, 40, 0., PtHistoUpperBound, 100, 0, dEdxM_UpLim);
//  tuple->PostS_CutIdVsPtVsIas = dir.make<TH3F>("PostS_CutIdVsPtVsIas", ";NCuts;p_{T} (GeV);G_{i}^{strips}", NCuts, 0, NCuts, 40, 0., PtHistoUpperBound, 10, 0., 1.);
//  tuple->PostS_CutIdVsPtVsIh = dir.make<TH3F>("PostS_CutIdVsPtVsIh", ";NCuts;p_{T} (GeV);I_{h} (MeV/cm)", NCuts, 0, NCuts, 40, 0., PtHistoUpperBound, 100, 0, dEdxM_UpLim);
//  if (TypeMode > 1) {
//    tuple->PostS_CutIdVsTOF = dir.make<TH2F>("PostS_CutIdVsTOF", ";NCuts;TOF", NCuts, 0, NCuts, 50, 1, 5);
//    tuple->PostS_CutIdVsTOFVsIas = dir.make<TH3F>("PostS_CutIdVsTOFVsIas", ";NCuts;TOF;G_{i}^{strips}", NCuts, 0, NCuts, 50, 0., 5, 10, 0., 1.);
//    tuple->PostS_CutIdVsTOFVsIh = dir.make<TH3F>("PostS_CutIdVsTOFVsIh", ";NCuts;TOF;I_{h} (MeV/cm)", NCuts, 0, NCuts, 50, 0., 5, 100, 0, dEdxM_UpLim);
//  }
  
  tuple->PostS_RelativePtShift = dir.make<TH1F>("PostS_RelativePtShift", ";#Delta p_{T} / p_{T}; Events / bin", 20, 0., 0.1);
  
  tuple->PostS_Ias = dir.make<TH1F>("PostS_Ias", ";G_{i}^{strips};Events / 0.1", 10, 0, dEdxS_UpLim);
  tuple->PostS_FiStrips = dir.make<TH1F>("PostS_FiStrips", ";F_{i}^{strips};Events / 0.1", 10, 0., 1.);
  tuple->PostS_FiStripsLog = dir.make<TH1F>("PostS_FiStripsLog", ";-log(1-F_{i}^{strips});Events / 0.05", 120, 0., 6.);
  tuple->PostS_IasVsFiStrips = dir.make<TH2F>("PostS_IasVsFiStrips", ";G_{i}^{strips};F_{i}^{strips}", 20, 0., 1.,20, 0., 1.);
  tuple->PostS_ProbQNoL1 = dir.make<TH1F>("PostS_ProbQNoL1", ";F_{i}^{pixels};Events / bin", 20, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas = dir.make<TH2F>("PostS_ProbQNoL1VsIas", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsFiStrips = dir.make<TH2F>("PostS_ProbQNoL1VsFiStrips", ";F_{i}^{pixels};F_{i}^{strips};Events",20, 0., 1., 20, 0., 1.);
  
  tuple->PostS_ProbQNoL1VsIasVsPt = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Pileup_up = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Pileup_up", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Pileup_down = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Pileup_down", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_ProbQNoL1_up = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_ProbQNoL1_up", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_ProbQNoL1_down = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_ProbQNoL1_down", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Ias_up = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Ias_up", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Ias_down = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Ias_down", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Pt_up = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Pt_up", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Pt_down = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Pt_down", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Trigger_up = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Trigger_up", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsIasVsPt_Trigger_down = dir.make<TH3F>("PostS_ProbQNoL1VsIasVsPt_Trigger_down", ";F_{i}^{pixels};G_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  
  tuple->PostS_ProbQNoL1VsFiStripsVsPt = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Pileup_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Pileup_up", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Pileup_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Pileup_down", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_ProbQNoL1_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_ProbQNoL1_up", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_ProbQNoL1_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_ProbQNoL1_down", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Ias_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Ias_up", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Ias_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Ias_down", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Pt_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Pt_up", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Pt_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Pt_down", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Trigger_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Trigger_up", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsVsPt_Trigger_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsVsPt_Trigger_down", ";F_{i}^{pixels};F_{i}^{strips};p_{T} (GeV)",20, 0., 1., 100, 0., 1.,160, 0., PtHistoUpperBound);

  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Pileup_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Pileup_up", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Pileup_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Pileup_down", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_ProbQNoL1_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_ProbQNoL1_up", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_ProbQNoL1_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_ProbQNoL1_down", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Ias_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Ias_up", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Ias_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Ias_down", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Pt_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Pt_up", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Pt_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Pt_down", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Trigger_up = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Trigger_up", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);
  tuple->PostS_ProbQNoL1VsFiStripsLogVsPt_Trigger_down = dir.make<TH3F>("PostS_ProbQNoL1VsFiStripsLogVsPt_Trigger_down", ";F_{i}^{pixels};-log(1-F_{i}^{strips});p_{T} (GeV)",20, 0., 1., 120, 0., 6.,160, 0., PtHistoUpperBound);

  // Inclusive 2D plots for Alphabet
  tuple->PostS_ProbQNoL1VsIas_Pileup_up = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Pileup_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_Pileup_down = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Pileup_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_ProbQNoL1_up = dir.make<TH2F>("PostS_ProbQNoL1VsIas_ProbQNoL1_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_ProbQNoL1_down = dir.make<TH2F>("PostS_ProbQNoL1VsIas_ProbQNoL1_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_Pt_up = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Pt_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_Pt_down = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Pt_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_Ias_up = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Ias_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_Ias_down = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Ias_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_Trigger_up = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Trigger_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_ProbQNoL1VsIas_Trigger_down = dir.make<TH2F>("PostS_ProbQNoL1VsIas_Trigger_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  
  tuple->PostS_SR1_Ias = dir.make<TH1F>("PostS_SR1_Ias", ";G_{i}^{strips};Events / 0.1", 10, 0, dEdxS_UpLim);
  tuple->PostS_SR1_ProbQNoL1 = dir.make<TH1F>("PostS_SR1_ProbQNoL1", ";F_{i}^{pixels};Events / bin", 20, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Pileup_up = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Pileup_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Pileup_down = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Pileup_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_ProbQNoL1_up = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_ProbQNoL1_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_ProbQNoL1_down = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_ProbQNoL1_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Pt_up = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Pt_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Pt_down = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Pt_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Ias_up = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Ias_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Ias_down = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Ias_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Trigger_up = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Trigger_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR1_ProbQNoL1VsIas_Trigger_down = dir.make<TH2F>("PostS_SR1_ProbQNoL1VsIas_Trigger_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  
  
  tuple->PostS_SR2_Ias = dir.make<TH1F>("PostS_SR2_Ias", ";G_{i}^{strips};Events / 0.1", 10, 0, dEdxS_UpLim);
  tuple->PostS_SR2_ProbQNoL1 = dir.make<TH1F>("PostS_SR2_ProbQNoL1", ";F_{i}^{pixels};Events / bin", 20, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Pileup_up = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Pileup_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Pileup_down = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Pileup_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_ProbQNoL1_up = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_ProbQNoL1_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_ProbQNoL1_down = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_ProbQNoL1_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Pt_up = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Pt_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Pt_down = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Pt_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Ias_up = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Ias_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Ias_down = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Ias_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Trigger_up = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Trigger_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR2_ProbQNoL1VsIas_Trigger_down = dir.make<TH2F>("PostS_SR2_ProbQNoL1VsIas_Trigger_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  
  tuple->PostS_SR2_ProbQNoL1VsIasVsMass = dir.make<TH3F>("PostS_SR2_ProbQNoL1VsIasVsMass", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.,80,0.,4000.);
  
  tuple->PostS_SR3_Ias = dir.make<TH1F>("PostS_SR3_Ias", ";G_{i}^{strips};Events / 0.1", 10, 0, dEdxS_UpLim);
  tuple->PostS_SR3_ProbQNoL1 = dir.make<TH1F>("PostS_SR3_ProbQNoL1", ";F_{i}^{pixels};Events / bin", 20, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Pileup_up = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Pileup_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Pileup_down = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Pileup_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_ProbQNoL1_up = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_ProbQNoL1_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_ProbQNoL1_down = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_ProbQNoL1_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Pt_up = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Pt_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Pt_down = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Pt_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Ias_up = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Ias_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Ias_down = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Ias_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Trigger_up = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Trigger_up", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);
  tuple->PostS_SR3_ProbQNoL1VsIas_Trigger_down = dir.make<TH2F>("PostS_SR3_ProbQNoL1VsIas_Trigger_down", ";F_{i}^{pixels};G_{i}^{strips};Events",20, 0., 1., 50, 0., 1.);

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
      
//      if (false) {
        // I think these are not used
//        tuple->PDF_G_EtaP = dir.make<TH3F>("PDF_G_EtaP", ";NCuts;PDF_G_Eta;P", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);
//        tuple->PDF_C_EtaP = dir.make<TH3F>("PDF_C_EtaP", ";NCuts;PDF_C_Eta;P", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);

//        tuple->PDF_A_Eta = dir.make<TH2F>("PDF_A_Eta", ";NCuts;PDF_A_Eta", NCuts, 0, NCuts, EtaBins, -3., 3.);
//        tuple->PDF_E_Eta = dir.make<TH2F>("PDF_E_Eta", ";NCuts;PDF_E_Eta", NCuts, 0, NCuts, EtaBins, -3., 3.);

//        tuple->PDF_B_EtaICK = dir.make<TH3F>("PDF_B_EtaICK", ";NCuts;PDF_B_EtaICK", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);
//        tuple->PDF_F_EtaICK = dir.make<TH3F>("PDF_F_EtaICK", ";NCuts;PDF_F_EtaICK", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);

//        tuple->PDF_H_EtaMass = dir.make<TH3F>("PDF_H_EtaMass", ";NCuts;PDF_H_Eta;Mass", NCuts, 0, NCuts, EtaBins, -3., 3., MassNBins, 0, MassHistoUpperBound);

        //pz FLIP
//        tuple->PDF_G_EtaP_Flip = dir.make<TH3F>("PDF_G_EtaP_Flip", ";NCuts;PDF_G_EtaP_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);
//        tuple->PDF_C_EtaP_Flip = dir.make<TH3F>("PDF_C_EtaP_Flip", ";NCuts;PDF_C_EtaP_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 200, GlobalMinPt, PtHistoUpperBound);

//        tuple->PDF_A_Eta_Flip = dir.make<TH2F>("PDF_A_Eta_Flip", ";NCuts;PDF_A_Eta_Flip", NCuts, 0, NCuts, EtaBins, -3., 3.);
//        tuple->PDF_E_Eta_Flip = dir.make<TH2F>("PDF_E_Eta_Flip", ";NCuts;PDF_E_Eta_Flip", NCuts, 0, NCuts, EtaBins, -3., 3.);

//        tuple->PDF_B_EtaICK_Flip = dir.make<TH3F>("PDF_B_EtaICK_Flip", ";NCuts;PDF_B_EtaICK_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);
//        tuple->PDF_F_EtaICK_Flip = dir.make<TH3F>("PDF_F_EtaICK_Flip", ";NCuts;PDF_F_EtaICK_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., 60, -2., 3.);

//        tuple->PDF_H_EtaMass_Flip = dir.make<TH3F>("PDF_H_EtaMass_Flip", ";NCuts;PDF_H_EtaMass_Flip", NCuts, 0, NCuts, EtaBins, -3., 3., MassNBins, 0, MassHistoUpperBound);
//      }
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
    
    tuple->PostS_RecoHSCParticleType = dir.make<TH1F>("PostS_RecoHSCParticleType", ";;Tracks / category", 6, -0.5, 5.5);
    tuple->PostS_RecoHSCParticleType->GetXaxis()->SetBinLabel(1,"globalMuon");
    tuple->PostS_RecoHSCParticleType->GetXaxis()->SetBinLabel(2,"trackerMuon");
    tuple->PostS_RecoHSCParticleType->GetXaxis()->SetBinLabel(3,"matchedStandAloneMuon");
    tuple->PostS_RecoHSCParticleType->GetXaxis()->SetBinLabel(4,"standAloneMuon");
    tuple->PostS_RecoHSCParticleType->GetXaxis()->SetBinLabel(5,"innerTrack");
    tuple->PostS_RecoHSCParticleType->GetXaxis()->SetBinLabel(6,"unknown");
    
    tuple->HSCPE = dir.make<TH1F>("HSCPE", ";NCuts;HSCPE", NCuts, 0, NCuts);
    tuple->Mass = dir.make<TH2F>("Mass", ";NCuts;Mass", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
    tuple->MaxEventMass = dir.make<TH2F>("MaxEventMass", ";NCuts;MaxEventMass", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
    if (calcSyst_) {
      tuple->HSCPE_SystP = dir.make<TH1F>("HSCPE_SystP", ";NCuts;HSCPE_SystP", NCuts, 0, NCuts);
      tuple->HSCPE_SystI = dir.make<TH1F>("HSCPE_SystI", ";NCuts;HSCPE_SystI", NCuts, 0, NCuts);
      tuple->HSCPE_SystM = dir.make<TH1F>("HSCPE_SystM", ";NCuts;HSCPE_SystM", NCuts, 0, NCuts);
      tuple->HSCPE_SystPU = dir.make<TH1F>("HSCPE_SystPU", ";NCuts;HSCPE_SystPU", NCuts, 0, NCuts);
      tuple->HSCPE_SystHUp = dir.make<TH1F>("HSCPE_SystHUp", ";NCuts;HSCPE_SystHUp", NCuts, 0, NCuts);
      tuple->HSCPE_SystHDown = dir.make<TH1F>("HSCPE_SystHDown", ";NCuts;HSCPE_SystHDown", NCuts, 0, NCuts);

      tuple->Mass_SystP = dir.make<TH2F>("Mass_SystP", ";NCuts;Mass_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MaxEventMass_SystP = dir.make<TH2F>("MaxEventMass_SystP", ";NCuts;MaxEventMass_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->Mass_SystI = dir.make<TH2F>("Mass_SystI", ";NCuts;Mass_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MaxEventMass_SystI = dir.make<TH2F>("MaxEventMass_SystI", ";NCuts;MaxEventMass_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->Mass_SystM = dir.make<TH2F>("Mass_SystM", ";NCuts;Mass_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MaxEventMass_SystM = dir.make<TH2F>("MaxEventMass_SystM", ";NCuts;MaxEventMass_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->Mass_SystPU = dir.make<TH2F>("Mass_SystPU", ";NCuts;Mass_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MaxEventMass_SystPU = dir.make<TH2F>("MaxEventMass_SystPU", ";NCuts;MaxEventMass_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->Mass_SystHUp = dir.make<TH2F>("Mass_SystHUp", ";NCuts;Mass_SystHUp", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MaxEventMass_SystHUp = dir.make<TH2F>("MaxEventMass_SystHUp", ";NCuts;MaxEventMass_SystHUp", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->Mass_SystHDown = dir.make<TH2F>("Mass_SystHDown", ";NCuts;Mass_SystHDown", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MaxEventMass_SystHDown = dir.make<TH2F>("MaxEventMass_SystHDown", ";NCuts;MaxEventMass_SystHDown", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
    }
    tuple->Mass_Flip = dir.make<TH2F>("Mass_Flip", ";NCuts;Mass_Flip", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
    
    if (TypeMode > 1) {
      tuple->HSCPE_SystT = dir.make<TH1F>("HSCPE_SystT", ";NCuts;HSCPE_SystT", NCuts, 0, NCuts);
      
      tuple->MassComb_Flip = dir.make<TH2F>("MassComb_Flip", ";NCuts;MassComb_Flip", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF_SystH = dir.make<TH2F>("MassTOF_SystH", ";NCuts;MassTOF_SystH", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb_SystHUp = dir.make<TH2F>("MassComb_SystHUp", ";NCuts;MassComb_SystHUp", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb_SystPU = dir.make<TH2F>("MassComb_SystPU", ";NCuts;MassComb_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->Mass_SystT = dir.make<TH2F>("Mass_SystT", ";NCuts;Mass_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF_SystT = dir.make<TH2F>("MassTOF_SystT", ";NCuts;MassTOF_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb_SystT = dir.make<TH2F>("MassComb_SystT", ";NCuts;MassComb_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MaxEventMass_SystT = dir.make<TH2F>("MaxEventMass_SystT", ";NCuts;MaxEventMass_SystT", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb_SystI = dir.make<TH2F>("MassComb_SystI", ";NCuts;MassComb_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF_SystM = dir.make<TH2F>("MassTOF_SystM", ";NCuts;MassTOF_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb_SystM = dir.make<TH2F>("MassComb_SystM", ";NCuts;MassComb_SystM", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb = dir.make<TH2F>("MassComb", ";NCuts;MassComb", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF_SystP = dir.make<TH2F>("MassTOF_SystP", ";NCuts;MassTOF_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb_SystP = dir.make<TH2F>("MassComb_SystP", ";NCuts;MassComb_SystP", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF = dir.make<TH2F>("MassTOF", ";NCuts;MassTOF", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassComb_SystHDown = dir.make<TH2F>("MassComb_SystHDown", ";NCuts;MassComb_SystHDown", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF_Flip = dir.make<TH2F>("MassTOF_Flip", ";NCuts;MassTOF_Flip", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF_SystI = dir.make<TH2F>("MassTOF_SystI", ";NCuts;MassTOF_SystI", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
      tuple->MassTOF_SystPU = dir.make<TH2F>("MassTOF_SystPU", ";NCuts;MassTOF_SystPU", NCuts, 0, NCuts, MassNBins, 0, MassHistoUpperBound);
    }


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

      //  I think this is not use, let's remove it if nobody complains
//    tuple->CtrlPt_S1_Is = dir.make<TH1D>("CtrlPt_S1_Is", "CtrlPt_S1_Is", 200, 0, dEdxS_UpLim);
//    tuple->CtrlPt_S2_Is = dir.make<TH1D>("CtrlPt_S2_Is", "CtrlPt_S2_Is", 200, 0, dEdxS_UpLim);
//    tuple->CtrlPt_S3_Is = dir.make<TH1D>("CtrlPt_S3_Is", "CtrlPt_S3_Is", 200, 0, dEdxS_UpLim);
//    tuple->CtrlPt_S4_Is = dir.make<TH1D>("CtrlPt_S4_Is", "CtrlPt_S4_Is", 200, 0, dEdxS_UpLim);
//
//    tuple->CtrlPt_S1_Ih = dir.make<TH1D>("CtrlPt_S1_Ih", "CtrlPt_S1_Ih", 400, 0, dEdxM_UpLim);
//    tuple->CtrlPt_S2_Ih = dir.make<TH1D>("CtrlPt_S2_Ih", "CtrlPt_S2_Ih", 400, 0, dEdxM_UpLim);
//    tuple->CtrlPt_S3_Ih = dir.make<TH1D>("CtrlPt_S3_Ih", "CtrlPt_S3_Ih", 400, 0, dEdxM_UpLim);
//    tuple->CtrlPt_S4_Ih = dir.make<TH1D>("CtrlPt_S4_Ih", "CtrlPt_S4_Ih", 400, 0, dEdxM_UpLim);
//
//    tuple->CtrlIs_S1_TOF = dir.make<TH1D>("CtrlIs_S1_TOF", "CtrlIs_S1_TOF", 200, 0, 5);
//    tuple->CtrlIs_S2_TOF = dir.make<TH1D>("CtrlIs_S2_TOF", "CtrlIs_S2_TOF", 200, 0, 5);
//    tuple->CtrlIs_S3_TOF = dir.make<TH1D>("CtrlIs_S3_TOF", "CtrlIs_S3_TOF", 200, 0, 5);
//    tuple->CtrlIs_S4_TOF = dir.make<TH1D>("CtrlIs_S4_TOF", "CtrlIs_S4_TOF", 200, 0, 5);
//
//    tuple->CtrlIh_S1_TOF = dir.make<TH1D>("CtrlIh_S1_TOF", "CtrlIh_S1_TOF", 200, 0, 5);
//    tuple->CtrlIh_S2_TOF = dir.make<TH1D>("CtrlIh_S2_TOF", "CtrlIh_S2_TOF", 200, 0, 5);
//    tuple->CtrlIh_S3_TOF = dir.make<TH1D>("CtrlIh_S3_TOF", "CtrlIh_S3_TOF", 200, 0, 5);
//    tuple->CtrlIh_S4_TOF = dir.make<TH1D>("CtrlIh_S4_TOF", "CtrlIh_S4_TOF", 200, 0, 5);
//
//    tuple->CtrlPt_S1_TOF = dir.make<TH1D>("CtrlPt_S1_TOF", "CtrlPt_S1_TOF", 200, -2, 7);
//    tuple->CtrlPt_S2_TOF = dir.make<TH1D>("CtrlPt_S2_TOF", "CtrlPt_S2_TOF", 200, -2, 7);
//    tuple->CtrlPt_S3_TOF = dir.make<TH1D>("CtrlPt_S3_TOF", "CtrlPt_S3_TOF", 200, -2, 7);
//    tuple->CtrlPt_S4_TOF = dir.make<TH1D>("CtrlPt_S4_TOF", "CtrlPt_S4_TOF", 200, -2, 7);
//
//    for (int i = 0; i < PredBins; i++) {
//      char Suffix[1024];
//      sprintf(Suffix, "_%i", i);
//
//      Name.append(Suffix);
//      tuple->CtrlPt_S1_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S1_TOF_Binned", "CtrlPt_S1_TOF_Binned", 200, -2, 7);
//
//      Name.append(Suffix);
//      tuple->CtrlPt_S2_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S2_TOF_Binned", "CtrlPt_S2_TOF_Binned", 200, -2, 7);
//
//      Name.append(Suffix);
//      tuple->CtrlPt_S3_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S3_TOF_Binned", "CtrlPt_S3_TOF_Binned", 200, -2, 7);
//
//      Name.append(Suffix);
//      tuple->CtrlPt_S4_TOF_Binned[std::to_string(i)] = dir.make<TH1D>("CtrlPt_S4_TOF_Binned", "CtrlPt_S4_TOF_Binned", 200, -2, 7);
//    }

 // Calibration
 // Scale Factors
    tuple->SF_HHit2DPix_loose   =  dir.make<TH2F>("SF_HHit2DPix_loose", "SF_HHit2DPix_loose",  50, 0., 100., 200, 0., 20.);
    tuple->SF_HHit2DStrip_loose =  dir.make<TH2F>("SF_HHit2DStrip_loose", "SF_HHit2DStrip_loose",  50, 0., 100., 200, 0., 20.);
    tuple->SF_HHit2DPix   =  dir.make<TH2F>("SF_HHit2DPix", "SF_HHit2DPix",  50, 0., 100., 200, 0., 20.);
    tuple->SF_HHit2DStrip =  dir.make<TH2F>("SF_HHit2DStrip", "SF_HHit2DStrip",  50, 0., 100., 200, 0., 20.);

  // K and C
    tuple->K_and_C_Ih_noL1_VsP_loose1 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_loose1","K_and_C_Ih_noL1_VsP_loose1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_loose2 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_loose2","K_and_C_Ih_noL1_VsP_loose2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta1_loose1 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta1_loose1","K_and_C_Ih_noL1_VsP_eta1_loose1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta1_loose2 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta1_loose2","K_and_C_Ih_noL1_VsP_eta1_loose2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta2_loose1 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta2_loose1","K_and_C_Ih_noL1_VsP_eta2_loose1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta2_loose2 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta2_loose2","K_and_C_Ih_noL1_VsP_eta2_loose2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_strip_VsP_loose1 = dir.make<TH2F>("K_and_C_Ih_strip_VsP_loose1","K_and_C_Ih_strip_VsP_loose1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_strip_VsP_loose2 = dir.make<TH2F>("K_and_C_Ih_strip_VsP_loose2","K_and_C_Ih_strip_VsP_loose2", 250,0,50, 80, 2.,10.);

    tuple->K_and_C_Ih_noL1_VsP_1 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_1","K_and_C_Ih_noL1_VsP_1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_2 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_2","K_and_C_Ih_noL1_VsP_2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta1_1 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta1_1","K_and_C_Ih_noL1_VsP_eta1_1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta1_2 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta1_2","K_and_C_Ih_noL1_VsP_eta1_2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta2_1 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta2_1","K_and_C_Ih_noL1_VsP_eta2_1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_eta2_2 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_eta2_2","K_and_C_Ih_noL1_VsP_eta2_2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_strip_VsP_1 = dir.make<TH2F>("K_and_C_Ih_strip_VsP_1","K_and_C_Ih_strip_VsP_1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_strip_VsP_2 = dir.make<TH2F>("K_and_C_Ih_strip_VsP_2","K_and_C_Ih_strip_VsP_2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_noFcut1 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_noFcut1","K_and_C_Ih_noL1_VsP_noFcut1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_noL1_VsP_noFcut2 = dir.make<TH2F>("K_and_C_Ih_noL1_VsP_noFcut2","K_and_C_Ih_noL1_VsP_noFcut2", 250,0,50, 80, 2.,10.);
    tuple->K_and_C_Ih_strip_VsP_noFcut1 = dir.make<TH2F>("K_and_C_Ih_strip_VsP_noFcut1","K_and_C_Ih_strip_VsP_noFcut1", 50,0,5, 80, 2.,10.);
    tuple->K_and_C_Ih_strip_VsP_noFcut2 = dir.make<TH2F>("K_and_C_Ih_strip_VsP_noFcut2","K_and_C_Ih_strip_VsP_noFcut2", 250,0,50, 80, 2.,10.);

    tuple->K_and_C_Kin_Mass = dir.make<TH1F>("K_and_C_Kin_Mass",";Mass (GeV)",100,0.,5.);
    tuple->K_and_C_Kin_p = dir.make<TH1F>("K_and_C_Kin_p",";p (GeV)",50,0.,5.);
    tuple->K_and_C_Kin_phi = dir.make<TH1F>("K_and_C_Kin_phi",";#phi",  24, -1.*acos(-1),acos(-1));
    tuple->K_and_C_Kin_eta = dir.make<TH1F>("K_and_C_Kin_eta",";#eta", 18, -2.25, 2.25);

  // Stability
    tuple->Stab_Ih_NoL1_VsRun    = dir.make<TH2F>("Stab_Ih_NoL1_VsRun",";Run;I_{h,NoL1}", 545, 271000,325500, 60, 0.,15.);
    tuple->Stab_Ih_pixNoL1_VsRun = dir.make<TH2F>("Stab_Ih_pixNoL1_VsRun",";Run;I_{h,NoPixL1}", 545, 271000,325500, 60, 0.,15.);
    tuple->Stab_Ih_strip_VsRun   = dir.make<TH2F>("Stab_Ih_strip_VsRun",";Run;dEdX(strip)", 545, 271000,325500, 60, 0.,15.);
    tuple->Stab_Gi_strip_VsRun   = dir.make<TH2F>("Stab_Gi_strip_VsRun",";Run;Gi(strip)", 545, 271000,325500, 80, 0.,1.);
    tuple->Stab_Gi_NoL1_VsRun    = dir.make<TH2F>("Stab_Gi_NoL1_VsRun", ";Run;G_{i}^{NoL1}", 545, 271000,325500, 80, 0.,1.);
    tuple->Stab_Fi_pixNoL1_VsRun = dir.make<TH2F>("Stab_Fi_pixNoL1_VsRun", ";Run;F_{i}^{Pixels}", 545, 271000,325500, 81, 0.,1.0125);
    tuple->Stab_invB_VsRun       = dir.make<TH2F>("Stab_invB_VsRun",";Run;invBeta", 545, 271000,325500, 90,-1,2);
    tuple->Stab_invB_DT_VsRun    = dir.make<TH2F>("Stab_invB_DT_VsRun",";Run;invBeta(DT)", 545, 271000,325500, 90,-1,2);
    tuple->Stab_invB_CSC_VsRun   = dir.make<TH2F>("Stab_invB_CSC_VsRun",";Run;invBeta(CSC)", 545, 271000,325500, 90,-1,2);



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
    tuple->Tree->Branch("BunchXing", &tuple->Tree_BunchXing);
    tuple->Tree->Branch("nPU", &tuple->Tree_nPU);
    tuple->Tree->Branch("nPUmean", &tuple->Tree_nPUmean);
    tuple->Tree->Branch("nofVtx", &tuple->Tree_nofVertices, "nofVtx/i");
    tuple->Tree->Branch("npv", &tuple->Tree_npv, "npv/I");

    tuple->Tree->Branch("pvX", &tuple->Tree_pvX);
    tuple->Tree->Branch("pvY", &tuple->Tree_pvY);
    tuple->Tree->Branch("pvZ", &tuple->Tree_pvZ);
    tuple->Tree->Branch("pvRho", &tuple->Tree_pvRho);
    tuple->Tree->Branch("pvNdof", &tuple->Tree_pvNdof);
    tuple->Tree->Branch("pvChi2", &tuple->Tree_pvChi2);
    tuple->Tree->Branch("pvSumPt2", &tuple->Tree_pvSumPt2);

    tuple->Tree->Branch("Hscp", &tuple->Tree_Hscp, "Hscp/i");
    tuple->Tree->Branch("nMuons", &tuple->Tree_nMuons, "nMuons/i");
    tuple->Tree->Branch("njets", &tuple->Tree_njets, "njets/i");
    tuple->Tree->Branch("Weight", &tuple->Tree_Weight, "Weight/F");
    tuple->Tree->Branch("GeneratorWeight", &tuple->Tree_GeneratorWeight, "GeneratorWeight/F");
    tuple->Tree->Branch("GeneratorBinningValues", &tuple->Tree_GeneratorBinningValues, "GeneratorBinningValues/F");
    tuple->Tree->Branch("triggerDecision", &tuple->Tree_triggerDecision);
    tuple->Tree->Branch("triggerHLTPrescale", &tuple->Tree_triggerHLTPrescale);

    tuple->Tree->Branch("triggerObjectE", &tuple->Tree_triggerObjectE);
    tuple->Tree->Branch("triggerObjectPt", &tuple->Tree_triggerObjectPt);
    tuple->Tree->Branch("triggerObjectEta", &tuple->Tree_triggerObjectEta);
    tuple->Tree->Branch("triggerObjectPhi", &tuple->Tree_triggerObjectPhi);

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
    tuple->Tree->Branch("RecoCaloMET", &tuple->Tree_RecoCaloMET, "RecoCaloMET/F");
    tuple->Tree->Branch("RecoCaloMET_phi", &tuple->Tree_RecoCaloMET_phi, "RecoCaloMET_phi/F");
    tuple->Tree->Branch("RecoCaloMET_sigf", &tuple->Tree_RecoCaloMET_sigf, "RecoCaloMET_sigf/F");
    tuple->Tree->Branch("RecoPFMET", &tuple->Tree_RecoPFMET, "RecoPFMET/F");
    tuple->Tree->Branch("RecoPFMET_phi", &tuple->Tree_RecoPFMET_phi, "RecoPFMET_phi/F");
    tuple->Tree->Branch("RecoPFMET_sigf", &tuple->Tree_RecoPFMET_sigf, "RecoPFMET_sigf/F");
    tuple->Tree->Branch("RecoPFMHT", &tuple->Tree_RecoPFMHT, "RecoPFMHT/F");
    tuple->Tree->Branch("HLTCaloMET", &tuple->Tree_HLTCaloMET, "HLTCaloMET/F");
    tuple->Tree->Branch("HLTCaloMET_phi", &tuple->Tree_HLTCaloMET_phi, "HLTCaloMET_phi/F");
    tuple->Tree->Branch("HLTCaloMET_sigf", &tuple->Tree_HLTCaloMET_sigf, "HLTCaloMET_sigf/F");
    tuple->Tree->Branch("HLTCaloMETClean", &tuple->Tree_HLTCaloMETClean, "HLTCaloMETClean/F");
    tuple->Tree->Branch("HLTCaloMETClean_phi", &tuple->Tree_HLTCaloMETClean_phi, "HLTCaloMETClean_phi/F");
    tuple->Tree->Branch("HLTCaloMETClean_sigf", &tuple->Tree_HLTCaloMETClean_sigf, "HLTCaloMETClean_sigf/F");
    tuple->Tree->Branch("HLTCaloMHT", &tuple->Tree_HLTCaloMHT, "HLTCaloMHT/F");
    tuple->Tree->Branch("HLTCaloMHT_phi", &tuple->Tree_HLTCaloMHT_phi, "HLTCaloMHT_phi/F");
    tuple->Tree->Branch("HLTCaloMHT_sigf", &tuple->Tree_HLTCaloMHT_sigf, "HLTCaloMHT_sigf/F");
    tuple->Tree->Branch("HLTPFMET", &tuple->Tree_HLTPFMET, "HLTPFMET/F");
    tuple->Tree->Branch("HLTPFMET_phi", &tuple->Tree_HLTPFMET_phi, "HLTPFMET_phi/F");
    tuple->Tree->Branch("HLTPFMET_sigf", &tuple->Tree_HLTPFMET_sigf, "HLTPFMET_sigf/F");
    tuple->Tree->Branch("HLTPFMHT", &tuple->Tree_HLTPFMHT, "HLTPFMHT/F");
    tuple->Tree->Branch("HLTPFMHT_phi", &tuple->Tree_HLTPFMHT_phi, "HLTPFMHT_phi/F");
    tuple->Tree->Branch("HLTPFMHT_sigf", &tuple->Tree_HLTPFMHT_sigf, "HLTPFMHT_sigf/F");
    tuple->Tree->Branch("matchedMuonWasFound", &tuple->Tree_matchedMuonWasFound, "matchedMuonWasFound/O");
    tuple->Tree->Branch("gParticleId", &tuple->Tree_gParticleId);
    tuple->Tree->Branch("gParticleStatus", &tuple->Tree_gParticleStatus);
    tuple->Tree->Branch("gParticleE", &tuple->Tree_gParticleE);
    tuple->Tree->Branch("gParticlePt", &tuple->Tree_gParticlePt);
    tuple->Tree->Branch("gParticlePz", &tuple->Tree_gParticlePz);
    tuple->Tree->Branch("gParticleEta", &tuple->Tree_gParticleEta);
    tuple->Tree->Branch("gParticlePhi", &tuple->Tree_gParticlePhi);
    tuple->Tree->Branch("gParticleBeta", &tuple->Tree_gParticleBeta);
    tuple->Tree->Branch("gParticleCharge", &tuple->Tree_gParticleCharge);
    tuple->Tree->Branch("gParticleProdVertexX", &tuple->Tree_gParticleProdVertexX);
    tuple->Tree->Branch("gParticleProdVertexY", &tuple->Tree_gParticleProdVertexY);
    tuple->Tree->Branch("gParticleProdVertexZ", &tuple->Tree_gParticleProdVertexZ);
    tuple->Tree->Branch("gParticleMotherId", &tuple->Tree_gParticleMotherId);
    tuple->Tree->Branch("gParticleMotherIndex", &tuple->Tree_gParticleMotherIndex);

    tuple->Tree->Branch("eleE", &tuple->Tree_eleE);
    tuple->Tree->Branch("elePt", &tuple->Tree_elePt);
    tuple->Tree->Branch("eleEta", &tuple->Tree_eleEta);
    tuple->Tree->Branch("elePhi", &tuple->Tree_elePhi);
    tuple->Tree->Branch("eleCharge", &tuple->Tree_eleCharge);
    tuple->Tree->Branch("eleE_SC", &tuple->Tree_eleE_SC);
    tuple->Tree->Branch("eleEta_SC", &tuple->Tree_eleEta_SC);
    tuple->Tree->Branch("elePhi_SC", &tuple->Tree_elePhi_SC);
    tuple->Tree->Branch("eleSigmaIetaIeta", &tuple->Tree_eleSigmaIetaIeta);
    tuple->Tree->Branch("eleFull5x5SigmaIetaIeta", &tuple->Tree_eleFull5x5SigmaIetaIeta);
    tuple->Tree->Branch("eleR9", &tuple->Tree_eleR9);
    tuple->Tree->Branch("ele_dEta", &tuple->Tree_ele_dEta);
    tuple->Tree->Branch("ele_dPhi", &tuple->Tree_ele_dPhi);
    tuple->Tree->Branch("ele_HoverE", &tuple->Tree_ele_HoverE);
    tuple->Tree->Branch("ele_d0", &tuple->Tree_ele_d0);
    tuple->Tree->Branch("ele_dZ", &tuple->Tree_ele_dZ);
    tuple->Tree->Branch("ele_pileupIso", &tuple->Tree_ele_pileupIso);
    tuple->Tree->Branch("ele_chargedIso", &tuple->Tree_ele_chargedIso);
    tuple->Tree->Branch("ele_photonIso", &tuple->Tree_ele_photonIso);
    tuple->Tree->Branch("ele_neutralHadIso", &tuple->Tree_ele_neutralHadIso);
    tuple->Tree->Branch("ele_MissHits", &tuple->Tree_ele_MissHits);
    tuple->Tree->Branch("ele_passCutBasedIDVeto", &tuple->Tree_ele_passCutBasedIDVeto);
    tuple->Tree->Branch("ele_passCutBasedIDLoose", &tuple->Tree_ele_passCutBasedIDLoose);
    tuple->Tree->Branch("ele_passCutBasedIDMedium", &tuple->Tree_ele_passCutBasedIDMedium);
    tuple->Tree->Branch("ele_passCutBasedIDTight", &tuple->Tree_ele_passCutBasedIDTight);
    tuple->Tree->Branch("ele_passMVAIsoIDWP80", &tuple->Tree_ele_passMVAIsoIDWP80);
    tuple->Tree->Branch("ele_passMVAIsoIDWP90", &tuple->Tree_ele_passMVAIsoIDWP90);
    tuple->Tree->Branch("ele_passMVAIsoIDWPHZZ", &tuple->Tree_ele_passMVAIsoIDWPHZZ);
    tuple->Tree->Branch("ele_passMVAIsoIDWPLoose", &tuple->Tree_ele_passMVAIsoIDWPLoose);
    tuple->Tree->Branch("ele_passMVANoIsoIDWP80", &tuple->Tree_ele_passMVANoIsoIDWP80);
    tuple->Tree->Branch("ele_passMVANoIsoIDWP90", &tuple->Tree_ele_passMVANoIsoIDWP90);
    tuple->Tree->Branch("ele_passMVANoIsoIDWPLoose", &tuple->Tree_ele_passMVANoIsoIDWPLoose);
    tuple->Tree->Branch("ele_PassConvVeto", &tuple->Tree_ele_PassConvVeto);
    tuple->Tree->Branch("ele_OneOverEminusOneOverP", &tuple->Tree_ele_OneOverEminusOneOverP);


    tuple->Tree->Branch("muonE", &tuple->Tree_muonE);
    tuple->Tree->Branch("muonPt", &tuple->Tree_muonPt);
    tuple->Tree->Branch("muonEta", &tuple->Tree_muonEta);
    tuple->Tree->Branch("muonPhi", &tuple->Tree_muonPhi);
    tuple->Tree->Branch("muonCharge", &tuple->Tree_muonCharge);
    tuple->Tree->Branch("muonIsLoose", &tuple->Tree_muonIsLoose);
    tuple->Tree->Branch("muonIsMedium", &tuple->Tree_muonIsMedium);
    tuple->Tree->Branch("muonIsTight", &tuple->Tree_muonIsTight);
    tuple->Tree->Branch("muon_d0", &tuple->Tree_muon_d0);
    tuple->Tree->Branch("muon_d0Err", &tuple->Tree_muon_d0Err);
    tuple->Tree->Branch("muon_dZ", &tuple->Tree_muon_dZ);
    tuple->Tree->Branch("muon_ip3d", &tuple->Tree_muon_ip3d);
    tuple->Tree->Branch("muon_ip3dSignificance", &tuple->Tree_muon_ip3dSignificance);
    tuple->Tree->Branch("muonType", &tuple->Tree_muonType);
    tuple->Tree->Branch("muonQuality", &tuple->Tree_muonQuality);
    tuple->Tree->Branch("muon_pileupIso", &tuple->Tree_muon_pileupIso);
    tuple->Tree->Branch("muon_chargedIso", &tuple->Tree_muon_chargedIso);
    tuple->Tree->Branch("muon_photonIso", &tuple->Tree_muon_photonIso);
    tuple->Tree->Branch("muon_neutralHadIso", &tuple->Tree_muon_neutralHadIso);
    tuple->Tree->Branch("muon_validFractionTrackerHits", &tuple->Tree_muon_validFractionTrackerHits);
    tuple->Tree->Branch("muTree_muon_normChi2onE", &tuple->Tree_muon_normChi2);
    tuple->Tree->Branch("muon_chi2LocalPosition", &tuple->Tree_muon_chi2LocalPosition);
    tuple->Tree->Branch("muon_kinkFinder", &tuple->Tree_muon_kinkFinder);
    tuple->Tree->Branch("muon_segmentCompatability", &tuple->Tree_muon_segmentCompatability);

    tuple->Tree->Branch("muon_trkIso", &tuple->Tree_muon_trkIso);
    tuple->Tree->Branch("muon_tuneP_Pt", &tuple->Tree_muon_tuneP_Pt);
    tuple->Tree->Branch("muon_tuneP_PtErr", &tuple->Tree_muon_tuneP_PtErr);
    tuple->Tree->Branch("muon_tuneP_Eta", &tuple->Tree_muon_tuneP_Eta);
    tuple->Tree->Branch("muon_tuneP_Phi", &tuple->Tree_muon_tuneP_Phi);
    tuple->Tree->Branch("muon_tuneP_MuonBestTrackType", &tuple->Tree_muon_tuneP_MuonBestTrackType);
    tuple->Tree->Branch("muon_isHighPtMuon", &tuple->Tree_muon_isHighPtMuon);
    tuple->Tree->Branch("muon_isTrackerHighPtMuon", &tuple->Tree_muon_isTrackerHighPtMuon);


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

      tuple->Tree->Branch("Jet_chargedHadronEnergyFraction", &tuple->Tree_jet_chargedHadronEnergyFraction);
      tuple->Tree->Branch("Jet_neutralHadronEnergyFraction", &tuple->Tree_jet_neutralHadronEnergyFraction);
      tuple->Tree->Branch("Jet_muonEnergyFraction", &tuple->Tree_jet_muonEnergyFraction);
      tuple->Tree->Branch("Jet_chargedMultiplicity", &tuple->Tree_jet_chargedMultiplicity);
      tuple->Tree->Branch("Jet_neutralMultiplicity", &tuple->Tree_jet_neutralMultiplicity);
      tuple->Tree->Branch("Jet_jetArea", &tuple->Tree_jet_jetArea);
      tuple->Tree->Branch("Jet_pileupE", &tuple->Tree_jet_pileupE);

    }
    tuple->Tree->Branch("mT", &tuple->Tree_vect_mT);
    if (saveTree > 1) {
      tuple->Tree->Branch("passCutPt55", &tuple->Tree_passCutPt55);
      tuple->Tree->Branch("passPreselection", &tuple->Tree_passPreselection);
      tuple->Tree->Branch("passPreselectionSept8", &tuple->Tree_passPreselectionSept8);
      tuple->Tree->Branch("passSelection", &tuple->Tree_passSelection);
    }
    tuple->Tree->Branch("isPFMuon", &tuple->Tree_isPFMuon);
    tuple->Tree->Branch("PFMuonPt", &tuple->Tree_PFMuonPt);
    tuple->Tree->Branch("Charge", &tuple->Tree_Charge);
    tuple->Tree->Branch("Pt", &tuple->Tree_Pt);
    tuple->Tree->Branch("PtErr", &tuple->Tree_PtErr);
    tuple->Tree->Branch("Is_StripOnly", &tuple->Tree_Is_StripOnly);
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
    tuple->Tree->Branch("isPhoton", &tuple->Tree_isPhoton);
    tuple->Tree->Branch("isElectron", &tuple->Tree_isElectron);
    tuple->Tree->Branch("isChHadron", &tuple->Tree_isChHadron);
    tuple->Tree->Branch("isNeutHadron", &tuple->Tree_isNeutHadron);
    tuple->Tree->Branch("isPfTrack", &tuple->Tree_isPfTrack);
    tuple->Tree->Branch("isUndefined", &tuple->Tree_isUndefined);
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
    tuple->Tree->Branch("p", &tuple->Tree_p);
    tuple->Tree->Branch("eta", &tuple->Tree_eta);
    tuple->Tree->Branch("phi", &tuple->Tree_phi);
    tuple->Tree->Branch("NOH", &tuple->Tree_NOH);
    tuple->Tree->Branch("NOPH", &tuple->Tree_NOPH);
    tuple->Tree->Branch("FOVH", &tuple->Tree_FOVH);
    tuple->Tree->Branch("NOMH", &tuple->Tree_NOMH);
    tuple->Tree->Branch("FOVHD", &tuple->Tree_FOVHD);
    tuple->Tree->Branch("NOM", &tuple->Tree_NOM);
    tuple->Tree->Branch("matchTrigMuon_minDeltaR", &tuple->Tree_matchTrigMuon_minDeltaR);
    tuple->Tree->Branch("matchTrigMuon_pT", &tuple->Tree_matchTrigMuon_pT);
    tuple->Tree->Branch("iso_TK", &tuple->Tree_iso_TK);
    tuple->Tree->Branch("iso_ECAL", &tuple->Tree_iso_ECAL);
    tuple->Tree->Branch("iso_HCAL", &tuple->Tree_iso_HCAL);
    tuple->Tree->Branch("track_genTrackMiniIsoSumPt", &tuple->Tree_track_genTrackMiniIsoSumPt);
    tuple->Tree->Branch("HSCP_tuneP_Pt", &tuple->Tree_HSCP_tuneP_Pt);
    tuple->Tree->Branch("HSCP_tuneP_PtErr", &tuple->Tree_HSCP_tuneP_PtErr);
    tuple->Tree->Branch("HSCP_tuneP_Eta", &tuple->Tree_HSCP_tuneP_Eta);
    tuple->Tree->Branch("HSCP_tuneP_Phi", &tuple->Tree_HSCP_tuneP_Phi);
    tuple->Tree->Branch("HSCP_tuneP_MuonBestTrackType", &tuple->Tree_HSCP_tuneP_MuonBestTrackType);
    tuple->Tree->Branch("HSCP_ErrorHisto_bin", &tuple->Tree_HSCP_ErrorHisto_bin);
    tuple->Tree->Branch("HSCP_type", &tuple->Tree_HSCP_type);


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
                                  const std::vector<int>  &BunchXing,
                                  const std::vector<int>  &nPU,
                                  const std::vector<float>  &nPUmean,
                                  const unsigned int &nofVertices,
                                  const int &npv,
                                  const std::vector<float>  &pvX,
                                  const std::vector<float>  &pvY,
                                  const std::vector<float>  &pvZ,
                                  const std::vector<float>  &pvRho,
                                  const std::vector<int>  &pvNdof,
                                  const std::vector<float>  &pvChi2,
                                  const std::vector<float>  &pvSumPt2,
                                  const unsigned int &Hscp,
                                  const unsigned int &nMuons,
                                  const unsigned int &njets,
                                  const float &weight,
                                  const float &generator_weight,
                                  const float &generator_binning_values,
                                  const std::vector<bool> &triggerDecision,
                                  const std::vector<int> &triggerHLTPrescale,
                                  const std::vector<std::vector<float>> &triggerObjectE,
                                  const std::vector<std::vector<float>> &triggerObjectPt,
                                  const std::vector<std::vector<float>> &triggerObjectEta,
                                  const std::vector<std::vector<float>> &triggerObjectPhi,
                                  const bool &HLT_Mu50,
                                  const bool &HLT_PFMET120_PFMHT120_IDTight,
                                  const bool &HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                                  const bool &HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                                  const bool &HLT_MET105_IsoTrk50,
                                  const float &RecoCaloMET,
                                  const float &RecoCaloMET_phi,
                                  const float &RecoCaloMET_sigf,
                                  const float &RecoPFMET,
                                  const float &RecoPFMET_phi,
                                  const float &RecoPFMET_sigf,
                                  const float &RecoPFMHT,
                                  const float &HLTCaloMET,
                                  const float &HLTCaloMET_phi,
                                  const float &HLTCaloMET_sigf,
                                  const float &HLTCaloMETClean,
                                  const float &HLTCaloMETClean_phi,
                                  const float &HLTCaloMETClean_sigf,
                                  const float &HLTCaloMHT,
                                  const float &HLTCaloMHT_phi,
                                  const float &HLTCaloMHT_sigf,
                                  const float &HLTPFMET,
                                  const float &HLTPFMET_phi,
                                  const float &HLTPFMET_sigf,
                                  const float &HLTPFMHT,
                                  const float &HLTPFMHT_phi,
                                  const float &HLTPFMHT_sigf,
                                  const bool &matchedMuonWasFound,
                                  const std::vector<int> &gParticleId,
                                  const std::vector<int> &gParticleStatus,
                                  const std::vector<float> &gParticleE,
                                  const std::vector<float> &gParticlePt,
                                  const std::vector<float> &gParticlePz,
                                  const std::vector<float> &gParticleEta,
                                  const std::vector<float> &gParticlePhi,
                                  const std::vector<float> &gParticleBeta,
                                  const std::vector<int> &gParticleCharge,
                                  const std::vector<float> &gParticleProdVertexX,
                                  const std::vector<float> &gParticleProdVertexY,
                                  const std::vector<float> &gParticleProdVertexZ,
                                  const std::vector<int> &gParticleMotherId,
                                  const std::vector<int> &gParticleMotherIndex,
                                  const std::vector<float> &eleE,
                                  const std::vector<float> &elePt,
                                  const std::vector<float> &eleEta,
                                  const std::vector<float> &elePhi,
                                  const std::vector<float> &eleCharge,
                                  const std::vector<float> &eleE_SC,
                                  const std::vector<float> &eleEta_SC,
                                  const std::vector<float> &elePhi_SC,
                                  const std::vector<float> &eleSigmaIetaIeta,
                                  const std::vector<float> &eleFull5x5SigmaIetaIeta,
                                  const std::vector<float> &eleR9,
                                  const std::vector<float> &ele_dEta,
                                  const std::vector<float> &ele_dPhi,
                                  const std::vector<float> &ele_HoverE,
                                  const std::vector<float> &ele_d0,
                                  const std::vector<float> &ele_dZ,
                                  const std::vector<float> &ele_pileupIso,
                                  const std::vector<float> &ele_chargedIso,
                                  const std::vector<float> &ele_photonIso,
                                  const std::vector<float> &ele_neutralHadIso,
                                  const std::vector<int> &ele_MissHits,
                                  const std::vector<bool> &ele_passCutBasedIDVeto,
                                  const std::vector<bool> &ele_passCutBasedIDLoose,
                                  const std::vector<bool> &ele_passCutBasedIDMedium,
                                  const std::vector<bool> &ele_passCutBasedIDTight,
                                  const std::vector<bool> &ele_passMVAIsoIDWP80,
                                  const std::vector<bool> &ele_passMVAIsoIDWP90,
                                  const std::vector<bool> &ele_passMVAIsoIDWPHZZ,
                                  const std::vector<bool> &ele_passMVAIsoIDWPLoose,
                                  const std::vector<bool> &ele_passMVANoIsoIDWP80,
                                  const std::vector<bool> &ele_passMVANoIsoIDWP90,
                                  const std::vector<bool> &ele_passMVANoIsoIDWPLoose,
                                  const std::vector<bool> &ele_PassConvVeto,
                                  const std::vector<float> &ele_OneOverEminusOneOverP,
                                  const std::vector<float> &muonE,
                                  const std::vector<float> &muonPt,
                                  const std::vector<float> &muonEta,
                                  const std::vector<float> &muonPhi,
                                  const std::vector<int> &muonCharge,
                                  const std::vector<bool> &muonIsLoose,
                                  const std::vector<bool> &muonIsMedium,
                                  const std::vector<bool> &muonIsTight,
                                  const std::vector<float> &muon_d0,
                                  const std::vector<float> &muon_d0Err,
                                  const std::vector<float> &muon_dZ,
                                  const std::vector<float> &muon_ip3d,
                                  const std::vector<float> &muon_ip3dSignificance,
                                  const std::vector<unsigned int> &muonType,
                                  const std::vector<unsigned int> &muonQuality,
                                  const std::vector<float> &muon_pileupIso,
                                  const std::vector<float> &muon_chargedIso,
                                  const std::vector<float> &muon_photonIso,
                                  const std::vector<float> &muon_neutralHadIso,
                                  const std::vector<float> &muon_validFractionTrackerHits,
                                  const std::vector<float> &muon_normChi2,
                                  const std::vector<float> &muon_chi2LocalPosition,
                                  const std::vector<float> &muon_kinkFinder,
                                  const std::vector<float> &muon_segmentCompatability,
                                  const std::vector<float> &muon_trkIso,
                                  const std::vector<float> &muon_tuneP_Pt,
                                  const std::vector<float> &muon_tuneP_PtErr,
                                  const std::vector<float> &muon_tuneP_Eta,
                                  const std::vector<float> &muon_tuneP_Phi,
                                  const std::vector<int> &muon_tuneP_MuonBestTrackType,
                                  const std::vector<bool> &muon_isHighPtMuon,
                                  const std::vector<bool> &muon_isTrackerHighPtMuon,
                                  const std::vector<float> &Jet_pt,
                                  const std::vector<float> &Jet_eta,
                                  const std::vector<float> &Jet_phi,
                                  const std::vector<float> &Jet_mass,
                                  const std::vector<float> &Jet_energy,
                                  const std::vector<float> &Jet_pdgId,
                                  const std::vector<float> &Jet_et,
                                  const std::vector<float> &Jet_chargedEmEnergyFraction,
                                  const std::vector<float> &Jet_neutralEmEnergyFraction,
                                  const std::vector<float> &Jet_chargedHadronEnergyFraction,
                                  const std::vector<float> &Jet_neutralHadronEnergyFraction,
                                  const std::vector<float> &Jet_muonEnergyFraction,
                                  const std::vector<int> &Jet_chargedMultiplicity,
                                  const std::vector<int> &Jet_neutralMultiplicity,
                                  const std::vector<float> &Jet_jetArea,
                                  const std::vector<float> &Jet_pileupE,
                                  const std::vector<float> &vect_mT,
                                  const std::vector<bool> &passCutPt55,
                                  const std::vector<bool> &passPreselection,
                                  const std::vector<bool> &passPreselectionSept8,
                                  const std::vector<bool> &passSelection,
                                  const std::vector<bool> &isPFMuon,
                                  const std::vector<bool> &PFMuonPt,
                                  const std::vector<float> &Charge,
                                  const std::vector<float> &Pt,
                                  const std::vector<float> &PtErr,
                                  const std::vector<float> &Is_StripOnly,
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
                                  const std::vector<bool>  &isPhoton,
                                  const std::vector<bool>  &isElectron,
                                  const std::vector<bool>  &isChHadron,
                                  const std::vector<bool>  &isNeutHadron,
                                  const std::vector<bool>  &isPfTrack,
                                  const std::vector<bool>  &isUndefined,
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
                                  const std::vector<float> &p,
                                  const std::vector<float> &eta,
                                  const std::vector<float> &phi,
                                  const std::vector<unsigned int> &noh,
                                  const std::vector<unsigned int> &noph,
                                  const std::vector<float> &fovh,
                                  const std::vector<unsigned int> &nomh,
                                  const std::vector<float> &fovhd,
                                  const std::vector<unsigned int> &nom,
                                  const std::vector<float> &matchTrigMuon_minDeltaR,
                                  const std::vector<float> &matchTrigMuon_pT,
                                  const std::vector<float> &iso_TK,
                                  const std::vector<float> &iso_ECAL,
                                  const std::vector<float> &iso_HCAL,
                                  const std::vector<float> &track_genTrackMiniIsoSumPt,


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
                                  const std::vector<float> &genphi,
                                  const std::vector<float> &HSCP_tuneP_Pt,
                                  const std::vector<float> &HSCP_tuneP_PtErr,
                                  const std::vector<float> &HSCP_tuneP_Eta,
                                  const std::vector<float> &HSCP_tuneP_Phi,
                                  const std::vector<int> &HSCP_tuneP_MuonBestTrackType,
                                  const std::vector<int> &HSCP_ErrorHisto_bin,
                                  const std::vector<int> &HSCP_type
                                  ) {
  tuple->Tree_Trig = Trig;
  tuple->Tree_Run = Run;
  tuple->Tree_Event = Event;
  tuple->Tree_Lumi = Lumi;
  tuple->Tree_PileUp = PileUp;
  tuple->Tree_BunchXing = BunchXing;
  tuple->Tree_nPU = nPU;
  tuple->Tree_nPUmean = nPUmean;
  tuple->Tree_nofVertices = nofVertices;
  tuple->Tree_npv = npv;

  tuple->Tree_pvX = pvX;
  tuple->Tree_pvY = pvY;
  tuple->Tree_pvZ = pvZ;
  tuple->Tree_pvRho = pvRho;
  tuple->Tree_pvNdof = pvNdof;
  tuple->Tree_pvChi2 = pvChi2;
  tuple->Tree_pvSumPt2 = pvSumPt2;

  tuple->Tree_Hscp = Hscp;
  tuple->Tree_nMuons = nMuons;
  tuple->Tree_njets = njets;
  tuple->Tree_Weight = weight;
  tuple->Tree_GeneratorWeight = generator_weight;
  tuple->Tree_GeneratorBinningValues = generator_binning_values;
  tuple->Tree_triggerDecision = triggerDecision;
  tuple->Tree_triggerHLTPrescale = triggerHLTPrescale;
  tuple->Tree_triggerObjectE = triggerObjectE;
  tuple->Tree_triggerObjectPt = triggerObjectPt;
  tuple->Tree_triggerObjectEta = triggerObjectEta;
  tuple->Tree_triggerObjectPhi = triggerObjectPhi;
  tuple->Tree_HLT_Mu50 = HLT_Mu50;
  tuple->Tree_HLT_PFMET120_PFMHT120_IDTight = HLT_PFMET120_PFMHT120_IDTight;
  tuple->Tree_HLT_PFHT500_PFMET100_PFMHT100_IDTight = HLT_PFHT500_PFMET100_PFMHT100_IDTight;
  tuple->Tree_HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60;
  tuple->Tree_HLT_MET105_IsoTrk50 = HLT_MET105_IsoTrk50;
  tuple->Tree_RecoCaloMET = RecoCaloMET;
  tuple->Tree_RecoCaloMET_phi = RecoCaloMET_phi;
  tuple->Tree_RecoCaloMET_sigf = RecoCaloMET_sigf;
  tuple->Tree_RecoPFMET = RecoPFMET;
  tuple->Tree_RecoPFMET_phi = RecoPFMET_phi;
  tuple->Tree_RecoPFMET_sigf = RecoPFMET_sigf;
  tuple->Tree_RecoPFMHT = RecoPFMHT;
  tuple->Tree_HLTCaloMET = HLTCaloMET;
  tuple->Tree_HLTCaloMET_phi = HLTCaloMET_phi;
  tuple->Tree_HLTCaloMET_sigf = HLTCaloMET_sigf;
  tuple->Tree_HLTCaloMETClean = HLTCaloMETClean;
  tuple->Tree_HLTCaloMETClean_phi = HLTCaloMETClean_phi;
  tuple->Tree_HLTCaloMETClean_sigf = HLTCaloMETClean_sigf;
  tuple->Tree_HLTCaloMHT = HLTCaloMHT;
  tuple->Tree_HLTCaloMHT_phi = HLTCaloMHT_phi;
  tuple->Tree_HLTCaloMHT_sigf = HLTCaloMHT_sigf;
  tuple->Tree_HLTPFMET = HLTPFMET;
  tuple->Tree_HLTPFMET_phi = HLTPFMET_phi;
  tuple->Tree_HLTPFMET_sigf = HLTPFMET_sigf;
  tuple->Tree_HLTPFMHT = HLTPFMHT;
  tuple->Tree_HLTPFMHT_phi = HLTPFMHT_phi;
  tuple->Tree_HLTPFMHT_sigf = HLTPFMHT_sigf;
  tuple->Tree_matchedMuonWasFound = matchedMuonWasFound;
  tuple->Tree_gParticleId = gParticleId;
  tuple->Tree_gParticleStatus = gParticleStatus;
  tuple->Tree_gParticleE = gParticleE;
  tuple->Tree_gParticlePt = gParticlePt;
  tuple->Tree_gParticlePz = gParticlePz;
  tuple->Tree_gParticleEta = gParticleEta;
  tuple->Tree_gParticlePhi = gParticlePhi;
  tuple->Tree_gParticleBeta = gParticleBeta;
  tuple->Tree_gParticleCharge = gParticleCharge;
  tuple->Tree_gParticleProdVertexX = gParticleProdVertexX;
  tuple->Tree_gParticleProdVertexY = gParticleProdVertexY;
  tuple->Tree_gParticleProdVertexZ = gParticleProdVertexZ;
  tuple->Tree_gParticleMotherId = gParticleMotherId;
  tuple->Tree_gParticleMotherIndex = gParticleMotherIndex;


  tuple->Tree_eleE = eleE;
  tuple->Tree_elePt = elePt;
  tuple->Tree_eleEta = eleEta;
  tuple->Tree_elePhi = elePhi;
  tuple->Tree_eleCharge = eleCharge;
  tuple->Tree_eleE_SC = eleE_SC;
  tuple->Tree_eleEta_SC = eleEta_SC;
  tuple->Tree_elePhi_SC = elePhi_SC;
  tuple->Tree_eleSigmaIetaIeta = eleSigmaIetaIeta;
  tuple->Tree_eleFull5x5SigmaIetaIeta = eleFull5x5SigmaIetaIeta;
  tuple->Tree_eleR9 = eleR9;
  tuple->Tree_ele_dEta = ele_dEta;
  tuple->Tree_ele_dPhi = ele_dPhi;
  tuple->Tree_ele_HoverE = ele_HoverE;
  tuple->Tree_ele_d0 = ele_d0;
  tuple->Tree_ele_dZ = ele_dZ;
  tuple->Tree_ele_pileupIso = ele_pileupIso;
  tuple->Tree_ele_chargedIso = ele_chargedIso;
  tuple->Tree_ele_photonIso = ele_photonIso;
  tuple->Tree_ele_neutralHadIso = ele_neutralHadIso;
  tuple->Tree_ele_MissHits = ele_MissHits;
  tuple->Tree_ele_passCutBasedIDVeto = ele_passCutBasedIDVeto;
  tuple->Tree_ele_passCutBasedIDLoose = ele_passCutBasedIDLoose;
  tuple->Tree_ele_passCutBasedIDMedium = ele_passCutBasedIDMedium;
  tuple->Tree_ele_passCutBasedIDTight = ele_passCutBasedIDTight;
  tuple->Tree_ele_passMVAIsoIDWP80 = ele_passMVAIsoIDWP80;
  tuple->Tree_ele_passMVAIsoIDWP90 = ele_passMVAIsoIDWP90;
  tuple->Tree_ele_passMVAIsoIDWPHZZ = ele_passMVAIsoIDWPHZZ;
  tuple->Tree_ele_passMVAIsoIDWPLoose = ele_passMVAIsoIDWPLoose;
  tuple->Tree_ele_passMVANoIsoIDWP80 = ele_passMVANoIsoIDWP80;
  tuple->Tree_ele_passMVANoIsoIDWP90 = ele_passMVANoIsoIDWP90;
  tuple->Tree_ele_passMVANoIsoIDWPLoose = ele_passMVANoIsoIDWPLoose;
  tuple->Tree_ele_PassConvVeto = ele_PassConvVeto;
  tuple->Tree_ele_OneOverEminusOneOverP = ele_OneOverEminusOneOverP;






  tuple->Tree_muonE = muonE;
  tuple->Tree_muonPt = muonPt;
  tuple->Tree_muonEta = muonEta;
  tuple->Tree_muonPhi = muonPhi;
  tuple->Tree_muonCharge = muonCharge;
  tuple->Tree_muonIsLoose = muonIsLoose;
  tuple->Tree_muonIsMedium = muonIsMedium;
  tuple->Tree_muonIsTight = muonIsTight;
  tuple->Tree_muon_d0 = muon_d0;
  tuple->Tree_muon_d0Err = muon_d0Err;
  tuple->Tree_muon_dZ = muon_dZ;
  tuple->Tree_muon_ip3d = muon_ip3d;
  tuple->Tree_muon_ip3dSignificance = muon_ip3dSignificance;
  tuple->Tree_muonType = muonType;
  tuple->Tree_muonQuality = muonQuality;
  tuple->Tree_muon_pileupIso = muon_pileupIso;
  tuple->Tree_muon_chargedIso = muon_chargedIso;
  tuple->Tree_muon_photonIso = muon_photonIso;
  tuple->Tree_muon_neutralHadIso = muon_neutralHadIso;
  tuple->Tree_muon_validFractionTrackerHits = muon_validFractionTrackerHits;
  tuple->Tree_muon_normChi2 = muon_normChi2;
  tuple->Tree_muon_chi2LocalPosition = muon_chi2LocalPosition;
  tuple->Tree_muon_kinkFinder = muon_kinkFinder;
  tuple->Tree_muon_segmentCompatability = muon_segmentCompatability;

  tuple->Tree_muon_trkIso = muon_trkIso;
  tuple->Tree_muon_tuneP_Pt = muon_tuneP_Pt;
  tuple->Tree_muon_tuneP_PtErr = muon_tuneP_PtErr;
  tuple->Tree_muon_tuneP_Eta = muon_tuneP_Eta;
  tuple->Tree_muon_tuneP_Phi = muon_tuneP_Phi;
  tuple->Tree_muon_tuneP_MuonBestTrackType = muon_tuneP_MuonBestTrackType;
  tuple->Tree_muon_isHighPtMuon = muon_isHighPtMuon;
  tuple->Tree_muon_isTrackerHighPtMuon = muon_isTrackerHighPtMuon;

  tuple->Tree_jet_pt = Jet_pt;
  tuple->Tree_jet_eta = Jet_eta;
  tuple->Tree_jet_phi = Jet_phi;
  tuple->Tree_jet_mass = Jet_mass;
  tuple->Tree_jet_energy = Jet_energy;
  tuple->Tree_jet_pdgId = Jet_pdgId;
  tuple->Tree_jet_et = Jet_et;
  tuple->Tree_jet_chargedEmEnergyFraction = Jet_chargedEmEnergyFraction;
  tuple->Tree_jet_neutralEmEnergyFraction = Jet_neutralEmEnergyFraction;
  tuple->Tree_jet_chargedHadronEnergyFraction = Jet_chargedHadronEnergyFraction;
  tuple->Tree_jet_neutralHadronEnergyFraction = Jet_neutralHadronEnergyFraction;
  tuple->Tree_jet_muonEnergyFraction = Jet_muonEnergyFraction;
  tuple->Tree_jet_chargedMultiplicity = Jet_chargedMultiplicity;
  tuple->Tree_jet_neutralMultiplicity = Jet_neutralMultiplicity;
  tuple->Tree_jet_jetArea = Jet_jetArea;
  tuple->Tree_jet_pileupE = Jet_pileupE;
  tuple->Tree_vect_mT = vect_mT;
  tuple->Tree_passCutPt55 = passCutPt55;
  tuple->Tree_passPreselection = passPreselection;
  tuple->Tree_passPreselectionSept8 = passPreselectionSept8;
  tuple->Tree_passSelection = passSelection;
  tuple->Tree_isPFMuon = isPFMuon;
  tuple->Tree_PFMuonPt = PFMuonPt;
  tuple->Tree_Charge = Charge;
  tuple->Tree_Pt = Pt;
  tuple->Tree_PtErr = PtErr;
  tuple->Tree_Is_StripOnly = Is_StripOnly;
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
  tuple->Tree_isPhoton = isPhoton;
  tuple->Tree_isElectron = isElectron;
  tuple->Tree_isChHadron = isChHadron;
  tuple->Tree_isNeutHadron = isNeutHadron;
  tuple->Tree_isPfTrack = isPfTrack;
  tuple->Tree_isUndefined = isUndefined;
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
  tuple->Tree_p = p;
  tuple->Tree_eta = eta;
  tuple->Tree_phi = phi;
  tuple->Tree_NOH = noh;
  tuple->Tree_NOPH = noph;
  tuple->Tree_FOVH = fovh;
  tuple->Tree_NOMH = nomh;
  tuple->Tree_FOVHD = fovhd;
  tuple->Tree_NOM = nom;
  tuple->Tree_matchTrigMuon_minDeltaR = matchTrigMuon_minDeltaR;
  tuple->Tree_matchTrigMuon_pT = matchTrigMuon_pT;
  tuple->Tree_iso_TK = iso_TK;
  tuple->Tree_iso_ECAL = iso_ECAL;
  tuple->Tree_iso_HCAL = iso_HCAL;
  tuple->Tree_track_genTrackMiniIsoSumPt = track_genTrackMiniIsoSumPt;
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

  tuple->Tree_HSCP_tuneP_Pt = HSCP_tuneP_Pt;
  tuple->Tree_HSCP_tuneP_PtErr = HSCP_tuneP_PtErr;
  tuple->Tree_HSCP_tuneP_Eta = HSCP_tuneP_Eta;
  tuple->Tree_HSCP_tuneP_Phi = HSCP_tuneP_Phi;
  tuple->Tree_HSCP_tuneP_MuonBestTrackType = HSCP_tuneP_MuonBestTrackType;
  tuple->Tree_HSCP_ErrorHisto_bin = HSCP_ErrorHisto_bin;
  tuple->Tree_HSCP_type = HSCP_type;


  // Save in the tree
  tuple->Tree->Fill();
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

  // Save in the gen tree
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

//  I think this is not use, let's remove it if nobody complains
//  if (!isCosmicSB) {
//    if (track->pt() > PtLimits[0]) {
//      tuple->CtrlPt_S4_Is->Fill(Is, Event_Weight);
//      tuple->CtrlPt_S4_Ih->Fill(Ih, Event_Weight);
//      if (tof)
//        tuple->CtrlPt_S4_TOF->Fill(MuonTOF, Event_Weight);
//      if (tof && bin >= 0 && bin < MaxPredBins)
//        tuple->CtrlPt_S4_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
//    } else if (track->pt() > PtLimits[1]) {
//      tuple->CtrlPt_S3_Is->Fill(Is, Event_Weight);
//      tuple->CtrlPt_S3_Ih->Fill(Ih, Event_Weight);
//      if (tof)
//        tuple->CtrlPt_S3_TOF->Fill(MuonTOF, Event_Weight);
//      if (tof && bin >= 0 && bin < MaxPredBins)
//        tuple->CtrlPt_S3_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
//    } else if (track->pt() > PtLimits[2]) {
//      tuple->CtrlPt_S2_Is->Fill(Is, Event_Weight);
//      tuple->CtrlPt_S2_Ih->Fill(Ih, Event_Weight);
//      if (tof)
//        tuple->CtrlPt_S2_TOF->Fill(MuonTOF, Event_Weight);
//      if (tof && bin >= 0 && bin < MaxPredBins)
//        tuple->CtrlPt_S2_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
//    } else {
//      tuple->CtrlPt_S1_Is->Fill(Is, Event_Weight);
//      tuple->CtrlPt_S1_Ih->Fill(Ih, Event_Weight);
//      if (tof)
//        tuple->CtrlPt_S1_TOF->Fill(MuonTOF, Event_Weight);
//      if (tof && bin >= 0 && bin < MaxPredBins)
//        tuple->CtrlPt_S1_TOF_Binned[to_string(bin)]->Fill(MuonTOF, Event_Weight);
//    }
//
//    if (Is > 0.2) {
//      if (tof)
//        tuple->CtrlIs_S4_TOF->Fill(MuonTOF, Event_Weight);
//    } else if (Is > 0.1) {
//      if (tof)
//        tuple->CtrlIs_S3_TOF->Fill(MuonTOF, Event_Weight);
//    } else if (Is > 0.05) {
//      if (tof)
//        tuple->CtrlIs_S2_TOF->Fill(MuonTOF, Event_Weight);
//    } else {
//      if (tof)
//        tuple->CtrlIs_S1_TOF->Fill(MuonTOF, Event_Weight);
//    }
//
//    if (Ih > 4.4) {
//      if (tof)
//        tuple->CtrlIh_S4_TOF->Fill(MuonTOF, Event_Weight);
//    } else if (Ih > 4.1) {
//      if (tof)
//        tuple->CtrlIh_S3_TOF->Fill(MuonTOF, Event_Weight);
//    } else if (Ih > 3.8) {
//      if (tof)
//        tuple->CtrlIh_S2_TOF->Fill(MuonTOF, Event_Weight);
//    } else {
//      if (tof)
//        tuple->CtrlIh_S1_TOF->Fill(MuonTOF, Event_Weight);
//    }
//  }

  float Ick = 0;
  if (dedxMObj)
    Ick = GetIck(Ih, DeDxK, DeDxC);
  if (false) {
    // just to make the compiler happy
    std::cout << "Ick: " << Ick << std::endl;
  }

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
//      tuple->PDF_C_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz
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
//      tuple->PDF_B_EtaICK->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz
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
//      tuple->PDF_A_Eta->Fill(CutIndex, track->eta(), Event_Weight);  //pz

    } else if (!PassTOFCut && PassPtCut && PassICut) {  //Region H
      tuple->H_H->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_H_Binned[to_string(bin)]->Fill(CutIndex, Event_Weight);
      tuple->RegionH_Ias->Fill(CutIndex, Is, Event_Weight);
//      if (TypeMode == 2 && Ick > 0)
//        tuple->PDF_H_EtaMass->Fill(CutIndex, track->eta(), track->p() * sqrt(Ick), Event_Weight);  //pz
      //Pred_P->Fill(CutIndex,track->p(),        Event_Weight);
      //Pred_I->Fill(CutIndex,Ih,   Event_Weight);
      if (TypeMode == 2)
        tuple->PostS_CutIdVsEta_RegionH->Fill(CutIndex, track->eta());
    } else if (!PassTOFCut && PassPtCut && !PassICut) {  //Region G
      tuple->H_G->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
      tuple->PostS_CutIdVsEta_RegionG->Fill(CutIndex, track->eta());
//      if (TypeMode == 2) tuple->PDF_G_EtaP->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz
    } else if (!PassTOFCut && !PassPtCut && PassICut) {                             //Region F
      tuple->H_F->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_F_Binned[to_string(bin)]->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_I->Fill(CutIndex, Ih, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->PostS_CutIdVsEta_RegionF->Fill(CutIndex, track->eta());
//      if (TypeMode == 2)
//        tuple->PDF_F_EtaICK->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz

    } else if (!PassTOFCut && !PassPtCut && !PassICut) {  //Region E
      tuple->H_E->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaB->Fill(CutIndex, track->eta(), Event_Weight);
      tuple->PostS_CutIdVsEta_RegionE->Fill(CutIndex, track->eta());
//      if (TypeMode == 2)
//        tuple->PDF_E_Eta->Fill(CutIndex, track->eta(), Event_Weight);  //pz
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
//      tuple->PDF_C_EtaP_Flip->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz
      //Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
    } else if (PassTOFCut && !PassPtCut && PassICut) {  //Region B
      tuple->H_B_Flip->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_B_Binned_Flip[to_string(bin)]->Fill(CutIndex, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_I_Flip->Fill(CutIndex, Ih, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaS_Flip->Fill(CutIndex, track->eta(), Event_Weight);
//      tuple->PDF_B_EtaICK_Flip->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz
      //Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
    } else if (PassTOFCut && !PassPtCut && !PassICut) {  //Region A
      tuple->H_A_Flip->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_TOF_Flip->Fill(CutIndex, MuonTOF, Event_Weight);
      if (TypeMode < 2)
        tuple->Pred_EtaB_Flip->Fill(CutIndex, track->eta(), Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS2_Flip->Fill(CutIndex, track->eta(), Event_Weight);
//      tuple->PDF_A_Eta_Flip->Fill(CutIndex, track->eta(), Event_Weight);  //pz
    } else if (!PassTOFCut && PassPtCut && PassICut) {                    //Region H
      tuple->H_H_Flip->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_H_Binned_Flip[to_string(bin)]->Fill(CutIndex, Event_Weight);
      tuple->RegionH_Ias_Flip->Fill(CutIndex, Is, Event_Weight);
//      if (TypeMode == 2 && Ick > 0)
//        tuple->PDF_H_EtaMass_Flip->Fill(CutIndex, track->eta(), track->p() * sqrt(Ick), Event_Weight);  //pz

      //Pred_P_Flip->Fill(CutIndex,track->p(),        Event_Weight);
      //Pred_I_Flip->Fill(CutIndex,Ih,   Event_Weight);
    } else if (!PassTOFCut && PassPtCut && !PassICut) {  //Region G
      tuple->H_G_Flip->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaP_Flip->Fill(CutIndex, track->eta(), track->p(), Event_Weight);
//      if (TypeMode == 2)
//        tuple->PDF_G_EtaP_Flip->Fill(CutIndex, track->eta(), track->p(), Event_Weight);  //pz

    } else if (!PassTOFCut && !PassPtCut && PassICut) {  //Region F
      tuple->H_F_Flip->Fill(CutIndex, Event_Weight);
      if (bin > -1 && bin < MaxPredBins)
        tuple->H_F_Binned_Flip[to_string(bin)]->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_I_Flip->Fill(CutIndex, Ih, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaS_Flip->Fill(CutIndex, track->eta(), Event_Weight);
//      if (TypeMode == 2)
//        tuple->PDF_F_EtaICK_Flip->Fill(CutIndex, track->eta(), Ick, Event_Weight);  //pz
    } else if (!PassTOFCut && !PassPtCut && !PassICut) {                            //Region E
      tuple->H_E_Flip->Fill(CutIndex, Event_Weight);
      if (TypeMode == 2)
        tuple->Pred_EtaB_Flip->Fill(CutIndex, track->eta(), Event_Weight);
//      if (TypeMode == 2)
//        tuple->PDF_E_Eta_Flip->Fill(CutIndex, track->eta(), Event_Weight);  //pz
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
                             float Ias_quantiles[5],
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
    if(ias>=Ias_quantiles[4])                         tuple->rB_90ias100.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
  }else{
    if(ias<Ias_quantiles[0]) tuple->rC_ias50.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
    if(ias>=Ias_quantiles[0] && ias<Ias_quantiles[1]) tuple->rD_50ias60.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
    if(ias>=Ias_quantiles[1] && ias<Ias_quantiles[2]) tuple->rD_60ias70.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
    if(ias>=Ias_quantiles[2] && ias<Ias_quantiles[3]) tuple->rD_70ias80.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
    if(ias>=Ias_quantiles[3] && ias<Ias_quantiles[4]) tuple->rD_80ias90.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
    if(ias>=Ias_quantiles[0] && ias<Ias_quantiles[4]) tuple->rD_50ias90.fill(eta,p,pt,pterr,ih,ias,m,tof,w);
    if(ias>=Ias_quantiles[4])        {m<500?m=m:m=-1; tuple->rD_90ias100.fill(eta,p,pt,pterr,ih,ias,m,tof,w);} //blind in the last quantile
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
