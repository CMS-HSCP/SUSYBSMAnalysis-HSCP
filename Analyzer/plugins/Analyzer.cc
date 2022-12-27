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
// Modifications by Dylan Angie Frank Apparu
//                  and Tamas Almos Vami

// - 41p0: - Refactor so no tuple is needed in the preslection function
// - 41p1: - Further code cleaning
// - 41p2: - 1D plots for CR, include syst on probQ, add passPreSept8, dz/dxy for non-vertex tracks, plots for cosmics checks
// - 41p3: - Add Pt,Ias systematics,trigger syst measurement. Now we dont return if trigger is not passed, but it's part of the cutflow, check num tracks before after HLT mu matching
// - 41p4: - Fix syst plots
// - 41p5: Do trigger syst plots for diff etas, Ias syst in low p region, Gi template calibration, add met trig eff, POG version of track rescaling
// - 41p6: CreateGiTemplates is false
// - 41p7: Plot for number of candidates, extend dRMinPfJet to 3.2, add 3 SRs, add HLT match, relative pt shift, tuneP muon checks
// - 41p8: Cutflow fixes, add fraction of cleaned clusters, getting rid of pixelProbs[], make trigger syst more general
// - 41p9: Increase pT in SR, dont load histos that are not needed for this typeMode
// - 42p0: Add trigger syst plots, add PostS_ProbQNoL1VsIasVsPt, extend dRMinPfJet to 5.0
// - 42p1: Let in MET triggers, for some MET trigger studies, tuneP muon checks extended axis, fix CutFlow 1st bin
// - 42p2: Changes so the code is compatible with MET triggers
// - 42p3: Dont cut on dR of gen vs candidate, add MET vs HT plot, go back to Mu only triggers, Add EventCutFlow
// - 43p4: Have more bins in Gi for the PostS histos, fix EventCutFlow's technical bin, letting in innerTracks
// - 42p4: Really 43p4 was a typo, and I also now when back with the Gi bins for PostS histos
// - 42p5: 3 SRs with exclusive pT bins in them, MetOverHt plots, use genParticles not genParticlesSkimmed
// - 31p0: same as 42p5 but with MET triggers
// - 31p2: candidate dependent MetOverHt plot
// - 42p5: Back to 42p5:
// - 32p0: Technical run, after PR75 is merged
// - 42p6: Adding Calibration plot, adding templates, adding electrons collection to the ntuple
// - 42p7: Add Gen_BetaGamma_lowBetaGamma, fix BefPreS HSCP type axis, add BefTrig plots, add per eta, per beta trigger systs,
//         unify vertex handling (use highest sum pt one),

// v25 Dylan
// - add EoP in the ntuple
// - add jets info in the ntuple

#include "SUSYBSMAnalysis/Analyzer/plugins/Analyzer.h"

Analyzer::Analyzer(const edm::ParameterSet& iConfig)
    // Read config file
    : hscpToken_(consumes<vector<susybsm::HSCParticle>>(iConfig.getParameter<edm::InputTag>("HscpCollection"))),
      genTrackToken_(consumes<reco::TrackCollection>(edm::InputTag("generalTracks"))),
      hscpIsoToken_(
          consumes<edm::ValueMap<susybsm::HSCPIsolation>>(iConfig.getParameter<edm::InputTag>("HscpIsoCollection"))),
      muonSegmentToken_(
          consumes<susybsm::MuonSegmentCollection>(iConfig.getParameter<edm::InputTag>("MuonSegmentCollection"))),
      dedxToken_(consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("DedxCollection"))),
      dedxPrescaleToken_(consumes<edm::ValueMap<int>>(iConfig.getParameter<edm::InputTag>("DedxCollectionPrescale"))),
      muonTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("MuonTimeCollection"))),
      muonDtTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("MuonDtTimeCollection"))),
      muonCscTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("MuonCscTimeCollection"))),
      muonDtSegmentToken_(
          consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("MuonDtSegmentCollection"))),
      muonCscSegmentToken_(
          consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("MuonCscSegmentCollection"))),
      offlinePrimaryVerticesToken_(
          consumes<vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("OfflinePrimaryVerticesCollection"))),
      lumiScalersToken_(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("LumiScalers"))),
      refittedStandAloneMuonsToken_(
          consumes<vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("RefittedStandAloneMuonsCollection"))),
      offlineBeamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("OfflineBeamSpotCollection"))),
      muonToken_(consumes<vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("MuonCollection"))),
      conversionsToken_(consumes<vector<reco::Conversion> >(iConfig.getParameter<edm::InputTag>("conversions"))),
      electronToken_(consumes<reco::GsfElectronCollection>(iConfig.getParameter<edm::InputTag>("ElectronCollection"))),
      electron_cutbasedID_decisions_veto_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_cutbasedID_decisions_veto"))),
      electron_cutbasedID_decisions_loose_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_cutbasedID_decisions_loose"))),
      electron_cutbasedID_decisions_medium_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_cutbasedID_decisions_medium"))),
      electron_cutbasedID_decisions_tight_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_cutbasedID_decisions_tight"))),
      electron_mvaIsoID_decisions_wp80_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_mvaIsoID_decisions_wp80"))),
      electron_mvaIsoID_decisions_wp90_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_mvaIsoID_decisions_wp90"))),
      electron_mvaIsoID_decisions_wpHZZ_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_mvaIsoID_decisions_wpHZZ"))),
      electron_mvaIsoID_decisions_wpLoose_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_mvaIsoID_decisions_wpLoose"))),
      electron_mvaNoIsoID_decisions_wp80_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_mvaNoIsoID_decisions_wp80"))),
      electron_mvaNoIsoID_decisions_wp90_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_mvaNoIsoID_decisions_wp90"))),
      electron_mvaNoIsoID_decisions_wpLoose_Token_(consumes<edm::ValueMap<bool> >(iConfig.getParameter<edm::InputTag>("electron_mvaNoIsoID_decisions_wpLoose"))),
      triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
      triggerPrescalesToken_(consumes<pat::PackedTriggerPrescales>(iConfig.getParameter<edm::InputTag>("triggerPrescales"))),
      trigEventToken_(consumes<trigger::TriggerEvent>(iConfig.getParameter<edm::InputTag>("TriggerSummary"))),
      filterName_(iConfig.getParameter<std::string>("FilterName")),
      triggerPathNamesFile_(iConfig.getParameter<string>("TriggerPathNamesFile")),
      muonHLTFilterNamesFile_(iConfig.getParameter<string>("MuonHLTFilterNamesFile")),
      pfMETToken_(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("PfMET"))),
      pfJetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PfJet"))),
      caloMETToken_(consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("CaloMET"))),
      caloJetToken_(consumes<std::vector<reco::CaloJet>>(iConfig.getParameter<edm::InputTag>("CaloJet"))),
      pileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("PileupInfo"))),
      genParticleToken_(
          consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("GenParticleCollection"))),
      trackToGenToken_(consumes<edm::Association<reco::GenParticleCollection>>(
          iConfig.getParameter<edm::InputTag>("TrackToGenAssoc"))),
      pfCandToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("PfCand"))),
      genEventToken_(consumes<GenEventInfoProduct>(iConfig.getParameter<edm::InputTag>("GenCollection"))),
      // HLT triggers
      trigger_met_(iConfig.getUntrackedParameter<vector<string>>("Trigger_MET")),
      trigger_mu_(iConfig.getUntrackedParameter<vector<string>>("Trigger_Mu")),
      // =========Analysis parameters===============
      typeMode_(iConfig.getUntrackedParameter<int>("TypeMode")),
      sampleType_(iConfig.getUntrackedParameter<int>("SampleType")),
      sampleName_(iConfig.getUntrackedParameter<string>("SampleName")),
      period_(iConfig.getUntrackedParameter<string>("Period")),
      skipSelectionPlot_(iConfig.getUntrackedParameter<bool>("SkipSelectionPlot")),
      ptHistoUpperBound_(iConfig.getUntrackedParameter<double>("PtHistoUpperBound")),
      pHistoUpperBound_(iConfig.getUntrackedParameter<double>("PHistoUpperBound")),
      massHistoUpperBound_(iConfig.getUntrackedParameter<double>("MassHistoUpperBound")),
      massNBins_(iConfig.getUntrackedParameter<int>("MassNBins")),
      cutOnIPbound_(iConfig.getUntrackedParameter<double>("IPbound")),
      predBins_(iConfig.getUntrackedParameter<int>("PredBins")),
      etaBins_(iConfig.getUntrackedParameter<int>("EtaBins")),
      reg_etabins_(iConfig.getUntrackedParameter<int>("RegEtaBins")),
      reg_ihbins_(iConfig.getUntrackedParameter<int>("RegIhBins")),
      reg_pbins_(iConfig.getUntrackedParameter<int>("RegPBins")),
      reg_massbins_(iConfig.getUntrackedParameter<int>("RegMassBins")),
      dEdxS_UpLim_(iConfig.getUntrackedParameter<double>("DeDxS_UpLim")),
      dEdxM_UpLim_(iConfig.getUntrackedParameter<double>("DeDxM_UpLim")),
      numDzRegions_(iConfig.getUntrackedParameter<int>("DzRegions")),
      globalMaxEta_(iConfig.getUntrackedParameter<double>("GlobalMaxEta")),
      globalMinPt_(iConfig.getUntrackedParameter<double>("GlobalMinPt")),
      globalMinNOPH_(iConfig.getUntrackedParameter<int>("GlobalMinNOPH")),
      globalMinFOVH_(iConfig.getUntrackedParameter<double>("GlobalMinFOVH")),
      globalMinNOM_(iConfig.getUntrackedParameter<int>("GlobalMinNOM")),
      globalMaxChi2_(iConfig.getUntrackedParameter<double>("GlobalMaxChi2")),
      globalMaxEoP_(iConfig.getUntrackedParameter<double>("GlobalMaxEoP")),
      globalMaxDZ_(iConfig.getUntrackedParameter<double>("GlobalMaxDZ")),
      globalMaxDXY_(iConfig.getUntrackedParameter<double>("GlobalMaxDXY")),
      globalMaxTIsol_(iConfig.getUntrackedParameter<double>("GlobalMaxTIsol")),
      globalMinDeltaRminJet_(iConfig.getUntrackedParameter<double>("GlobalMinDeltaRminJet")),
      globalMaxMiniRelIsoAll_(iConfig.getUntrackedParameter<double>("GlobalMaxMiniRelIsoAll")),
      globalMinIh_(iConfig.getUntrackedParameter<double>("GlobalMinIh")),
      globalMinTrackProbQCut_(iConfig.getUntrackedParameter<double>("GlobalMinTrackProbQCut")),
      globalMaxTrackProbQCut_(iConfig.getUntrackedParameter<double>("GlobalMaxTrackProbQCut")),
      globalMinTrackProbXYCut_(iConfig.getUntrackedParameter<double>("GlobalMinTrackProbXYCut")),
      globalMaxTrackProbXYCut_(iConfig.getUntrackedParameter<double>("GlobalMaxTrackProbXYCut")),
      minMuStations_(iConfig.getUntrackedParameter<int>("MinMuStations")),
      globalMinIs_(iConfig.getUntrackedParameter<double>("GlobalMinIs")),
      globalMinTOF_(iConfig.getUntrackedParameter<double>("GlobalMinTOF")),
      puTreatment_(iConfig.getUntrackedParameter<bool>("PileUpTreatment")),
      createGiTemplates_(iConfig.getUntrackedParameter<bool>("CreateGiTemplates")),
      CreateAndRunGitemplates_(iConfig.getUntrackedParameter<bool>("CreateAndRunGitemplates")),
      NbPuBins_(iConfig.getUntrackedParameter<int>("NbPileUpBins")),
      PuBins_(iConfig.getUntrackedParameter<vector<int>>("PileUpBins")),

      useTemplateLayer_(iConfig.getUntrackedParameter<bool>("UseTemplateLayer")),
      dEdxSF_0_(iConfig.getUntrackedParameter<double>("DeDxSF_0")),
      dEdxSF_1_(iConfig.getUntrackedParameter<double>("DeDxSF_1")),
      dEdxK_(iConfig.getUntrackedParameter<double>("DeDxK")),
      dEdxC_(iConfig.getUntrackedParameter<double>("DeDxC")),
      dEdxTemplate_(iConfig.getUntrackedParameter<string>("DeDxTemplate")),
      timeOffset_(iConfig.getUntrackedParameter<string>("TimeOffset")),
      theFMIPX_(iConfig.getUntrackedParameter<double>("FMIPX")),
      saveTree_(iConfig.getUntrackedParameter<int>("SaveTree")),
      saveGenTree_(iConfig.getUntrackedParameter<int>("SaveGenTree")),
      pixelCPE_(iConfig.getParameter<std::string>("PixelCPE")),
      debug_(iConfig.getUntrackedParameter<int>("DebugLevel")),
      hasMCMatch_(iConfig.getUntrackedParameter<bool>("HasMCMatch")),
      calcSyst_(iConfig.getUntrackedParameter<bool>("CalcSystematics"))
 {
  //now do what ever initialization is needed
  // define the selection to be considered later for the optimization
  // WARNING: recall that this has a huge impact on the analysis time AND on the output file size --> be carefull with your choice

  useClusterCleaning = true;
  if (typeMode_ == 4) {
    useClusterCleaning = false;  //switch off cluster cleaning for mCHAMPs
  }

  isData = (sampleType_ == 0);
  isBckg = (sampleType_ == 1);
  isSignal = (sampleType_ >= 2);

  //dEdxSF [0] = dEdxSF_0_;
  //dEdxSF [1] = dEdxSF_1_;

  bool splitByModuleType = true;
  //Get bool if we want PU or not : PU_Treat, get nb of PU bins (all in cfg)

  dEdxTemplatesPU.resize(NbPuBins_, NULL);

  // Option for Gi to have PU dependence
  if(!createGiTemplates_){
    if(puTreatment_){
      for (int i = 0; i < NbPuBins_ ; i++){
        dEdxTemplatesPU[i] = loadDeDxTemplate(dEdxTemplate_, splitByModuleType,true,(i+1));
      }
    } else {
        dEdxTemplates = loadDeDxTemplate(dEdxTemplate_, splitByModuleType,false,0);
    }
  }
  if(CreateAndRunGitemplates_){
    dEdxTemplates = loadDeDxTemplate(dEdxTemplate_, splitByModuleType,false,0);
  }

  tofCalculator.loadTimeOffset(timeOffset_);
}

Analyzer::~Analyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob() {
  static constexpr const char* const MOD = "Analyzer";

  // Book histograms using TFileService
  edm::Service<TFileService> fs;
  TFileDirectory dir = fs->mkdir(sampleName_.c_str(), sampleName_.c_str());

  // create histograms & trees
  tuple = new Tuple();
  initializeCuts(fs, CutPt_, CutI_, CutTOF_, CutPt_Flip_, CutI_Flip_, CutTOF_Flip_);

  tuple_maker->initializeTuple(tuple,
                               dir,
                               saveTree_,
                               saveGenTree_,
                               calcSyst_,
                               typeMode_,
                               isSignal,
                               CutPt_.size(),
                               CutPt_Flip_.size(),
                               ptHistoUpperBound_,
                               massHistoUpperBound_,
                               massNBins_,
                               cutOnIPbound_,
                               predBins_,
                               etaBins_,
                               dEdxS_UpLim_,
                               dEdxM_UpLim_,
                               numDzRegions_,
                               globalMinPt_,
                               globalMinTOF_);

  tuple_maker->initializeRegions(tuple,
                                 dir,
                                 reg_etabins_,
                                 reg_ihbins_,
                                 reg_pbins_,
                                 reg_massbins_);

  // Re-weighting
  // Functions defined in Analyzer/interface/MCWeight.h
  if (!isData) {
    mcWeight = new MCWeight();
    mcWeight->loadPileupWeights(period_);
    mcWeight->getSampleWeights(period_, sampleName_.c_str(), IntegratedLuminosity_, CrossSection_);
  }

  loadSFPixel();

  // Set in Analyzer/interface/MCWeight.h
  // 58970.47 for 2018
  // 41809.45 for 2017
  // 35552.24 for 2016
  // TODO: this should be revised, currently 2016 is used, it should be made year dependent,
  //       prob best to be able to control it from the config too
  //       or if not, then it should not be a global variable, I think
  tuple->IntLumi->Fill(0.5, IntegratedLuminosity_);

  // Get cross section from Analyzer/interface/MCWeight.h file
  // The SampleName in the config has to contain the HSCP flavor and mass
  tuple->XSection->Fill(0.5, CrossSection_);

  tof = nullptr;
  dttof = nullptr;
  csctof = nullptr;

  CurrentRun_ = 0;
  RNG = new TRandom3();

  // TODO: This is needed when there is no PU reweighting, i.e. data
  // Should be revised
  //PUSystFactor_.clear();
  PUSystFactor_.resize(2, 1.);
  PUSystFactor_[0] = PUSystFactor_[1] = 0.;

  HSCPTk = new bool[CutPt_.size()];
  HSCPTk_SystP = new bool[CutPt_.size()];
  HSCPTk_SystI = new bool[CutPt_.size()];
  HSCPTk_SystT = new bool[CutPt_.size()];
  HSCPTk_SystM = new bool[CutPt_.size()];
  HSCPTk_SystPU = new bool[CutPt_.size()];
  HSCPTk_SystHUp = new bool[CutPt_.size()];
  HSCPTk_SystHDown = new bool[CutPt_.size()];
  MaxMass = new float[CutPt_.size()];
  MaxMass_SystP = new float[CutPt_.size()];
  MaxMass_SystI = new float[CutPt_.size()];
  MaxMass_SystT = new float[CutPt_.size()];
  MaxMass_SystM = new float[CutPt_.size()];
  MaxMass_SystPU = new float[CutPt_.size()];
  MaxMass_SystHUp = new float[CutPt_.size()];
  MaxMass_SystHDown = new float[CutPt_.size()];

    // Check if we are dealing with data or MC
  if (sampleType_ == 0 ) {
    if (debug_> 0) edm::LogPrint(MOD) << "This is data processing";
  } else if (sampleType_ == 1) {
    if (debug_> 0) edm::LogPrint(MOD) << "This is background MC processing";
  } else if (sampleType_ == 2) {
    if (debug_> 0) edm::LogPrint(MOD) << "This is signal MC processing";
  } else {
    if (debug_> 0) edm::LogPrint(MOD) << "This is syst studies";
  }
}

// ------------ method called for each event  ------------
void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  static constexpr const char* const MOD = "Analyzer";
  static constexpr const float cm2umUnit = 0.0001;
  using namespace edm;
  using namespace trigger;

  //if run change, update conditions
  if (CurrentRun_ != iEvent.id().run()) {
    CurrentRun_ = iEvent.id().run();
    tofCalculator.setRun(CurrentRun_);
    dEdxSF[0] = dEdxSF_0_;
    dEdxSF[1] = dEdxSF_1_;
  }

  // Compute event weight
  if (!isData) {
    float PUWeight = mcWeight->getEventPUWeight(iEvent, pileupInfoToken_, PUSystFactor_);
    EventWeight_ = PUWeight;  // 1. : unweighted w.r.t pileup
  } else {
    EventWeight_ = 1.;
  }

    // useful at the event level -- because used both for the HSCP and the track loops
    string year = period_;
    if(!isData) year="";
    int run_number=iEvent.id().run();
    bool usePixel = true;
    bool useStrip = true;
    bool useTruncated = false;
    bool mustBeInside = true;
    size_t MaxStripNOM = 99;
    bool correctFEDSat = false;
    int crossTalkInvAlgo = 1;
    float dropLowerDeDxValue = 0.0;
    bool skipPixelL1 = false;
    int  skip_templates_ias = 0;
    TH3F* localdEdxTemplates;
    if(!createGiTemplates_ && !puTreatment_){
      localdEdxTemplates = dEdxTemplates;
    }

//  if (debug_ > 1 ) LogPrint(MOD) << "\nThis is a new event. Weight factor applied: " << EventWeight_;
  float HSCPgenBeta1 = -1., HSCPgenBeta2 = -1.;

  //get generator weight and pthat
  if (!isData){
      const edm::Handle<GenEventInfoProduct> genEvt = iEvent.getHandle(genEventToken_);
      if (genEvt.isValid()){
          GeneratorWeight_ = genEvt->weight();
        if(genEvt->binningValues().size()>0) {
          GeneratorBinningValues_ = genEvt->binningValues()[0];
          tuple->Gen_Binning->Fill(GeneratorBinningValues_);
        }
      }
  }

  //get the collection of generated Particles
  vector<reco::GenParticle> genColl;
  if (!isData) {
    const edm::Handle<vector<reco::GenParticle>> genCollH = iEvent.getHandle(genParticleToken_);
    if (!genCollH.isValid()) {
      LogPrint(MOD) << "Invalid GenParticle collection, this event will be ignored!";
      return;
    } else {
      genColl = *genCollH;
    }
  }

  float SignalEventWeight = 1.0;
  if (isSignal) {
    int NChargedHSCP = HowManyChargedHSCP(genColl);
    float HSCPDLength1 = -1, HSCPDLength2 = -1;

    SignalEventWeight = mcWeight->getFGluinoWeight(NChargedHSCP, typeMode_);

    GetGenHSCPDecayLength(genColl, HSCPDLength1, HSCPDLength2, true);

    tuple->Gen_DecayLength->Fill(HSCPDLength1, SignalEventWeight);
    tuple->Gen_DecayLength->Fill(HSCPDLength2, SignalEventWeight);

    // Get gen level beta for charged R-hadrons only
    GetGenHSCPBeta(genColl, HSCPgenBeta1, HSCPgenBeta2, true);
    if (HSCPgenBeta1 >= 0)
      tuple->Gen_Beta_Charged->Fill(HSCPgenBeta1, SignalEventWeight);
    if (HSCPgenBeta2 >= 0)
      tuple->Gen_Beta_Charged->Fill(HSCPgenBeta2, SignalEventWeight);

    // R-hadron weights needed due to wrong GenId
    // Wa is additional weight for single other, Wad for other+double_charged,
    // Waa for the event with 2 other R-hadron, Wan for other+neutral
    float Wa = 1.0, Wad = 1.0, Waa = 1.0, Wan = 1.0;
    bool Rhadron = false;  // default value - not R-hadron (not need to weight)
    string sample_name = "";
    mcWeight->getRHadronWeights(sample_name, Rhadron, Wa, Wad, Waa, Wan);
    // R-hadron wights needed due to wrong GenId
    // TODO: this prob is not true anymore in UL samples
  }  //End of isSignal
  else if (isBckg) {
    float notHSCPDLength1 = -1, notHSCPDLength2 = -1;

    // This returns a lot of zeros, I think we should not stop with the fist 2 on the list, maybe the first 2 non-zero?
    GetGenBcgDecayLength(genColl, notHSCPDLength1, notHSCPDLength2, true);
    tuple->Gen_DecayLength->Fill(notHSCPDLength1, EventWeight_);
    tuple->Gen_DecayLength->Fill(notHSCPDLength2, EventWeight_);
  }


  //initialize counters: nw - wrong, na - other, nd - double charged, nn - neutral
  unsigned int nw = 0, na = 0, nd = 0, nn = 0;
  vector<float> genid;
  vector<float> gencharge;
  vector<float> genmass;
  vector<float> genpt;
  vector<float> geneta;
  vector<float> genphi;


  if (!isData) {
    for (auto const& gen : genColl) {
      int GenId = abs(gen.pdgId());
      if (isSignal && isHSCPgenID(gen)) {
        // Categorise event with R-hadrons for additional weighting
        if (GenId == 1000612 || GenId == 1092214) {
          nw += 1;  // count wrong
        } else if (GenId == 1006223 || GenId == 1092224) {
           nd += 1;  // count doble charged
        } else if ((GenId) == 1006113 || (GenId) == 1006333 || (GenId) == 1006313 || (GenId) == 1000622 ||
                   (GenId) == 1092114 || (GenId) == 1093324 || (GenId) == 1093214 || (GenId) == 1009333 ||
                   (GenId) == 1009223 || (GenId) == 1009113 || (GenId) == 1009313 || (GenId) == 1000993) {
          nn += 1;  // count neutral
        } else if ((GenId) > 1000000) {
          na += 1;
        }  // count other R-hadrons

        // Fill up pT, eta, and beta plots for gen-level HSCP particles
        tuple->Gen_pT->Fill(gen.pt(), SignalEventWeight);
        tuple->Gen_Eta->Fill(gen.eta(), SignalEventWeight);
        tuple->Gen_Beta->Fill(gen.p() / gen.energy(), SignalEventWeight);
        tuple->Gen_BetaGamma->Fill(gen.p() / gen.mass(), SignalEventWeight);
        tuple->Gen_BetaGamma_lowBetaGamma->Fill(gen.p() / gen.mass(), SignalEventWeight);

        // Variables for the tuple gen tree branch
        genid.push_back(gen.pdgId());
        gencharge.push_back(gen.charge());
        genmass.push_back(gen.mass());
        genpt.push_back(gen.pt());
        geneta.push_back(gen.eta());
        genphi.push_back(gen.phi());
      } else if (isBckg) {
        // Fill up pT, eta, and beta plots for gen-level background particles
        tuple->Gen_pT->Fill(gen.pt(), EventWeight_);
        tuple->Gen_Eta->Fill(gen.eta(), EventWeight_);
        tuple->Gen_Beta->Fill(gen.p() / gen.energy(), EventWeight_);
        tuple->Gen_BetaGamma->Fill(gen.p() / gen.mass(), EventWeight_);
        tuple->Gen_BetaGamma_lowBetaGamma->Fill(gen.p() / gen.mass(), EventWeight_);

        // Variables for the tuple gen tree branch
        genid.push_back(gen.pdgId());
        gencharge.push_back(gen.charge());
        genmass.push_back(gen.mass());
        genpt.push_back(gen.pt());
        geneta.push_back(gen.eta());
        genphi.push_back(gen.phi());
      }
    }

    if (saveGenTree_ > 0) {
      if (debug_ > 4 ) LogPrint(MOD) << "Fill GenTree with basics gen info";
      tuple_maker->fillGenTreeBranches(tuple,
                                       iEvent.id().run(),
                                       iEvent.id().event(),
                                       iEvent.id().luminosityBlock(),
                                       EventWeight_,
                                       GeneratorWeight_,
                                       GeneratorBinningValues_,
                                       genid,
                                       gencharge,
                                       genmass,
                                       genpt,
                                       geneta,
                                       genphi);
    }
  }

  //------------------------------------------------------------------
  // GenParticles to be saved in the HSCP candidate ntuples (Christina Wang)
  //------------------------------------------------------------------
  std::vector<int> gParticleId;
  std::vector<int> gParticleStatus;
  std::vector<float> gParticleE;
  std::vector<float> gParticlePt;
  std::vector<float> gParticlePz;
  std::vector<float> gParticleEta;
  std::vector<float> gParticlePhi;
  std::vector<float> gParticleBeta;
  std::vector<int> gParticleCharge;
  std::vector<float> gParticleProdVertexX;
  std::vector<float> gParticleProdVertexY;
  std::vector<float> gParticleProdVertexZ;
  std::vector<int> gParticleMotherId;
  std::vector<int> gParticleMotherIndex;

  std::vector<const reco::Candidate*> prunedV;//Allows easier comparison for mother finding

    for(size_t i=0; i<genColl.size();i++)
    {
      if(
         (abs(genColl[i].pdgId()) >= 1 && abs(genColl[i].pdgId()) <= 6 && ( genColl[i].status() < 30 ))
         || (abs(genColl[i].pdgId()) >= 11 && abs(genColl[i].pdgId()) <= 18)
         || (abs(genColl[i].pdgId()) == 21 && genColl[i].status() < 30)
         || (abs(genColl[i].pdgId()) == 22 && genColl[i].pt() > 10.0 )
         || (abs(genColl[i].pdgId()) >= 23 && abs(genColl[i].pdgId()) <= 25)
         || (abs(genColl[i].pdgId()) >= 32 && abs(genColl[i].pdgId()) <= 42)
         || (abs(genColl[i].pdgId()) == 1023)
         || (abs(genColl[i].pdgId()) >= 1000000)
         )
        {
      	  prunedV.push_back(&genColl[i]);
        }
    } //loop over all gen particles

    // //Look for mother particle and Fill gen variables
    for(unsigned int i = 0; i < prunedV.size(); i++) {
      gParticleId.push_back(prunedV[i]->pdgId());
      gParticleStatus.push_back(prunedV[i]->status());
      gParticleE.push_back(prunedV[i]->energy());
      gParticlePt.push_back(prunedV[i]->pt());
      gParticlePz.push_back(prunedV[i]->pz());
      gParticleEta.push_back(prunedV[i]->eta());
      gParticlePhi.push_back(prunedV[i]->phi());
      gParticleProdVertexX.push_back(prunedV[i]->vx());
      gParticleProdVertexY.push_back(prunedV[i]->vy());
      gParticleProdVertexZ.push_back(prunedV[i]->vz());
      gParticleBeta.push_back(prunedV[i]->p()/prunedV[i]->energy());
      gParticleCharge.push_back(prunedV[i]->charge());
      gParticleMotherId.push_back(0);
      gParticleMotherIndex.push_back(-1);

      if(prunedV[i]->numberOfMothers() > 0) {
        //find the ID of the first mother that has a different ID than the particle itself
        const reco::Candidate* firstMotherWithDifferentID = findFirstMotherWithDifferentID(prunedV[i]);
        if (firstMotherWithDifferentID)gParticleMotherId[i] = firstMotherWithDifferentID->pdgId();

        //find the mother and keep going up the mother chain if the ID's are the same
        const reco::Candidate* originalMotherWithSameID = findOriginalMotherWithSameID(prunedV[i]);
        for(unsigned int j = 0; j < prunedV.size(); j++) {
          if(prunedV[j] == originalMotherWithSameID) {
            gParticleMotherIndex[i] = j;
            break;
          }
        }
      }
      else {
        gParticleMotherIndex[i] = -1;
      }
    }

  //------------------------------------------------------------------
  // Get trigger results for this event
  //------------------------------------------------------------------
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  const auto triggerNames = iEvent.triggerNames(*triggerH);
  // if (isData) iEvent.getByToken(triggerPrescalesToken_, triggerPrescales);

  //edm::Handle<pat::PackedTriggerPrescales> triggerPrescales = iEvent.getHandle(triggerPrescalesToken_);

  //------------------------------------------------------------------
  //Read in HLT Trigger Path List from config file
  //------------------------------------------------------------------
  for (int i = 0; i<NTriggersMAX; ++i) triggerPathNames[i] = "";
  ifstream myfile (edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str()) ;
  if (myfile.is_open()) {
    std::string line;
    int index;
    std::string hltpathname;
    while(myfile>>index>>hltpathname)
    {
      if (index < NTriggersMAX) triggerPathNames[index] = hltpathname;
    }
    myfile.close();
  }
  else {
    std::cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(triggerPathNamesFile_.c_str()).fullPath().c_str() << "\n";
  }


  //------------------------------------------------------------------
  // Save trigger decisions in array of booleans
  //------------------------------------------------------------------
  std::vector<bool> triggerDecision;
  std::vector<int> triggerHLTPrescale;
  for (unsigned int j = 0; j < NTriggersMAX; ++j) {
    triggerDecision.push_back(false);
    triggerHLTPrescale.push_back(1);

  }

  for (unsigned int i = 0, n = triggerH->size(); i < n; ++i) {
    string hltPathNameReq = "HLT_";
    if ((triggerNames.triggerName(i)).find(hltPathNameReq) == string::npos) continue;
    if ((triggerNames.triggerName(i)).find_last_of("_") == string::npos) continue;
    int lastUnderscorePos = (triggerNames.triggerName(i)).find_last_of("_");
    string hltPathNameWithoutVersionNumber = (triggerNames.triggerName(i)).substr(0,lastUnderscorePos);

    for (unsigned int j = 0; j < NTriggersMAX; ++j) {
      if (triggerPathNames[j] == "") continue;
      if (hltPathNameWithoutVersionNumber == triggerPathNames[j]) {
      	triggerDecision[j] = triggerH->accept(i);
      	//if (isData) triggerHLTPrescale[j] = triggerPrescales->getPrescaleForIndex(i);
      }
    }
  }

  //0: neither mu nor met, 1: mu only, 2: met only, 3: mu and met
  trigInfo_ = 0;

  // These are used in the tree alone, otherwise we use passTriggerPatterns to check the triggers
  bool HLT_Mu50 = false;
  bool HLT_PFMET120_PFMHT120_IDTight = false;
  bool HLT_PFHT500_PFMET100_PFMHT100_IDTight = false;
  bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = false;
  bool HLT_MET105_IsoTrk50 = false;

  for (unsigned int i = 0; i < triggerH->size(); i++) {
    if (TString(triggerNames.triggerName(i)).Contains("HLT_Mu50_v") && triggerH->accept(i))
      HLT_Mu50 = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_PFMET120_PFMHT120_IDTight_v") && triggerH->accept(i))
      HLT_PFMET120_PFMHT120_IDTight = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_PFHT500_PFMET100_PFMHT100_IDTight_v") && triggerH->accept(i))
      HLT_PFHT500_PFMET100_PFMHT100_IDTight = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v") &&
        triggerH->accept(i))
      HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_MET105_IsoTrk50_v") && triggerH->accept(i))
      HLT_MET105_IsoTrk50 = true;
  }

  // Should this be a bin in error histo?
  if (HLT_Mu50 != triggerDecision.at(196)) cout<<"TRIGGER DECISION DOESN'T AGREE!!"<<endl;

  // Number of (re-weighted) events
  tuple->NumEvents->Fill(0., EventWeight_);
  // Number of (re-weighted with PU syst fact) events
  tuple->NumEvents->Fill(1., EventWeight_ * PUSystFactor_[0]);

  // Get handle for trigEvent
  edm::Handle<trigger::TriggerEvent> trigEvent = iEvent.getHandle(trigEventToken_);
  // Get muon collections
  vector<reco::Muon> muonColl = iEvent.get(muonToken_);

  // Match muon track to HLT muon track
  std::vector<TLorentzVector> trigObjP4s;
  trigtools::getP4sOfObsPassingFilter(trigObjP4s,*trigEvent,filterName_,"HLT");
  
  bool matchedMuonWasFound = false;
  
  float dr_min_hlt_muon_inEvent = 9999.0;
  for(size_t objNr=0; objNr<trigObjP4s.size(); objNr++) {
    if (trigObjP4s[objNr].Pt() < 50) continue;
    for (unsigned int i = 0; i < muonColl.size(); i++) {
      const reco::Muon* mu = &(muonColl)[i];
      if (mu->isStandAloneMuon() && !mu->isGlobalMuon())  continue;
      
      float dr_hltmu_muon = deltaR(trigObjP4s[objNr].Eta(),trigObjP4s[objNr].Phi(),mu->eta(),mu->phi());
      if (dr_hltmu_muon < dr_min_hlt_muon_inEvent){
        dr_min_hlt_muon_inEvent = dr_hltmu_muon;
      }
    }
  }
  if (HLT_Mu50) {
    tuple->dRMinHLTMuon->Fill(dr_min_hlt_muon_inEvent, EventWeight_);
  }
  
    // Number of events that pass the matching
  if (HLT_Mu50 && dr_min_hlt_muon_inEvent < 0.15) {
    matchedMuonWasFound = true;
    tuple->NumEvents->Fill(3., EventWeight_);
  }

  // Check if the event is passing trigger
  if (debug_ > 1) LogPrint(MOD) << "\nChecking if the LS " << iEvent.id().luminosityBlock() <<  " / event " << iEvent.id().event() << "  is passing trigger...";
  bool metTrig = passTriggerPatterns(triggerH, triggerNames, trigger_met_);
  bool muTrig = passTriggerPatterns(triggerH, triggerNames, trigger_mu_);
  
  if (muTrig && dr_min_hlt_muon_inEvent < 0.15) {
    trigInfo_ = 1;
    tuple->BefPreS_TriggerType->Fill(1., EventWeight_);
  }
  if (metTrig) {
    trigInfo_ = 2;
    tuple->BefPreS_TriggerType->Fill(2., EventWeight_);
  }
  if (metTrig || muTrig) { trigInfo_ = 3; }
  if (metTrig && muTrig) { trigInfo_ = 4; }

  // why not let the MET triggers into this plot? then make a new boolean for the trigger choices
//  bool triggerPassed = (trigInfo_ == 1);
  tuple->BefPreS_TriggerType->Fill(trigInfo_, EventWeight_);
  // If triggering is intended (not the case when we make ntuples)
  if (trigInfo_ > 0) {
    if (debug_ > 2 ) LogPrint(MOD) << " > This event passeed the needed triggers! trigInfo_ = " << trigInfo_;
    // Number of events that pass the trigger
    tuple->NumEvents->Fill(2., EventWeight_);
  }
  else {
    if (debug_ > 2 ) LogPrint(MOD) << " > This event did not pass the needed triggers";
  }
  // For TOF only analysis if the event doesn't pass the signal triggers check if it was triggered by the no BPTX cosmic trigger

  //keep beta distribution for signal after the trigger
  if (trigInfo_ > 0) {
    if (isSignal) {
      if (HSCPgenBeta1 >= 0)
        tuple->Gen_Beta_Triggered->Fill(HSCPgenBeta1, EventWeight_);
      if (HSCPgenBeta2 >= 0)
        tuple->Gen_Beta_Triggered->Fill(HSCPgenBeta2, EventWeight_);
    }
  }

  // Define handles for DeDx Hits, Muon TOF Combined, Muon TOF DT, Muon TOF CSC
  const edm::Handle<reco::DeDxHitInfoAss> dedxCollH = iEvent.getHandle(dedxToken_);
  edm::Handle<edm::ValueMap<int>> dedxPrescCollH; //  = iEvent.getHandle(dedxPrescaleToken_)
  iEvent.getByToken(dedxPrescaleToken_, dedxPrescCollH);

  const edm::Handle<reco::MuonTimeExtraMap> tofMap = iEvent.getHandle(muonTimeToken_);
  const edm::Handle<reco::MuonTimeExtraMap> tofDtMap = iEvent.getHandle(muonDtTimeToken_);
  const edm::Handle<reco::MuonTimeExtraMap> tofCscMap = iEvent.getHandle(muonCscTimeToken_);

  edm::Handle<CSCSegmentCollection> CSCSegmentCollH;
  edm::Handle<DTRecSegment4DCollection> DTSegmentCollH;

  // Retrieve tracker topology from the event setup
  edm::ESHandle<TrackerTopology> TopoHandle;
  iSetup.get<TrackerTopologyRcd>().get(TopoHandle);
  const TrackerTopology* tTopo = TopoHandle.product();
  
  // Retrieve tracker geometry from the event setup
  edm::ESHandle<TrackerGeometry> tkGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(tkGeometry);
  
  // Retrieve CPE from the event setup

  edm::ESHandle<PixelClusterParameterEstimator> pixelCPE;
  iSetup.get<TkPixelCPERecord>().get(pixelCPE_, pixelCPE);

  // Handles for track collection, PF candidates, PF MET and PF jets, and Calo jets
  const edm::Handle<reco::TrackCollection> trackCollectionHandle = iEvent.getHandle(genTrackToken_);;
  const edm::Handle<reco::PFCandidateCollection> pfCandHandle = iEvent.getHandle(pfCandToken_);
  const edm::Handle<std::vector<reco::PFMET>> recoPFMETHandle = iEvent.getHandle(pfMETToken_);
  const edm::Handle<reco::PFJetCollection> pfJetHandle = iEvent.getHandle(pfJetToken_);
  const edm::Handle<std::vector<reco::CaloJet>> caloJetHandle = iEvent.getHandle(caloJetToken_);
  const edm::Handle<std::vector<reco::CaloMET>> recoCaloMETHandle = iEvent.getHandle(caloMETToken_);

  //================= Handle For Muon DT/CSC Segment ===============
  if (!isBckg) {  //do not recompute TOF on MC background
    iEvent.getByToken(muonCscSegmentToken_, CSCSegmentCollH);
    if (!CSCSegmentCollH.isValid()) {
      LogError("Analyzer") << "CSC Segment Collection not found!";
      return;
    }

    iEvent.getByToken(muonDtSegmentToken_, DTSegmentCollH);
    if (!DTSegmentCollH.isValid()) {
      LogError("Analyzer") << "DT Segment Collection not found!";
      return;
    }
  }

  //===================== Handle For PileUp ================
  unsigned int pileup_fromLumi = 0;
  const edm::Handle<LumiScalersCollection> lumiScalers = iEvent.getHandle(lumiScalersToken_);
  if (lumiScalers.isValid() && !lumiScalers->empty()) {
    LumiScalersCollection::const_iterator scalit = lumiScalers->begin();
    pileup_fromLumi = scalit->pileup();
  }
  
  std::vector<int> bunchXing;
  std::vector<int> nPU;
  std::vector<float> nPUmean;
  edm::Handle<std::vector<PileupSummaryInfo> > puInfo;

  if (!isData) {
    iEvent.getByToken(pileupInfoToken_,puInfo);
    for(const PileupSummaryInfo &pu : *puInfo) {
         bunchXing.push_back(pu.getBunchCrossing());
         nPU.push_back(pu.getPU_NumInteractions());
         nPUmean.push_back(pu.getTrueNumInteractions());
    }
  }

  // Collection for vertices
  vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken_);
  int numGoodVerts = 0;
  std::vector<float> pvX;
  std::vector<float> pvY;
  std::vector<float> pvZ;
  std::vector<float> pvRho;
  std::vector<int> pvNdof;
  std::vector<float> pvChi2;
  std::vector<float> pvSumPt2;

  // Loop on the vertices in the event
  bool foundPV = false;
  int highestSumPt2VertexIndex = -1;
  reco::Vertex highestSumPt2Vertex = vertexColl.front();

  for (unsigned int i = 0; i < vertexColl.size(); i++) {
    if (vertexColl[i].isFake())continue;
    if (!vertexColl[i].isValid())continue;
    numGoodVerts++;
    
    pvX.push_back(vertexColl[i].x());
    pvY.push_back(vertexColl[i].y());
    pvZ.push_back(vertexColl[i].z());
    pvRho.push_back(vertexColl[i].position().rho());
    pvNdof.push_back(vertexColl[i].ndof());
    pvChi2.push_back(vertexColl[i].chi2());

    // calculate sum pT squared
    double sum = 0.;
    double pT;
     for (reco::Vertex::trackRef_iterator it = vertexColl[i].tracks_begin(); it != vertexColl[i].tracks_end(); it++) {
       pT = (**it).pt();
       double epT = (**it).ptError();
       pT = pT > epT ? pT - epT : 0;

       sum += pT * pT;
     }
     pvSumPt2.push_back(sum);

     if (!foundPV) {
        highestSumPt2Vertex = vertexColl[i];
        foundPV = true;
        highestSumPt2VertexIndex = i;
      }
  } // end loop on vertexColl

  // Needs to be defined at event level, used in several functions
  float RecoCaloMET = -10, RecoCaloMET_phi = -10, RecoCaloMET_sigf = -10;

  // Add all electrons
  edm::Handle<reco::GsfElectronCollection> electrons = iEvent.getHandle(electronToken_);

  std::vector<float> eleE;
  std::vector<float> elePt;
  std::vector<float> eleEta;
  std::vector<float> elePhi;
  std::vector<float> eleCharge;
  std::vector<float> eleE_SC;
  std::vector<float> eleEta_SC;
  std::vector<float> elePhi_SC;
  std::vector<float> eleSigmaIetaIeta;
  std::vector<float> eleFull5x5SigmaIetaIeta;
  std::vector<float> eleR9;
  std::vector<float> ele_dEta;
  std::vector<float> ele_dPhi;
  std::vector<float> ele_HoverE;
  std::vector<float> ele_d0;
  std::vector<float> ele_dZ;
  std::vector<float> ele_pileupIso;
  std::vector<float> ele_chargedIso;
  std::vector<float> ele_photonIso;
  std::vector<float> ele_neutralHadIso;
  std::vector<int> ele_MissHits;
  std::vector<bool> ele_passCutBasedIDVeto;
  std::vector<bool> ele_passCutBasedIDLoose;
  std::vector<bool> ele_passCutBasedIDMedium;
  std::vector<bool> ele_passCutBasedIDTight;
  std::vector<bool> ele_passMVAIsoIDWP80;
  std::vector<bool> ele_passMVAIsoIDWP90;
  std::vector<bool> ele_passMVAIsoIDWPHZZ;
  std::vector<bool> ele_passMVAIsoIDWPLoose;
  std::vector<bool> ele_passMVANoIsoIDWP80;
  std::vector<bool> ele_passMVANoIsoIDWP90;
  std::vector<bool> ele_passMVANoIsoIDWPLoose;
  std::vector<bool> ele_PassConvVeto;
  std::vector<float> ele_OneOverEminusOneOverP;


  iEvent.getByToken(electron_cutbasedID_decisions_veto_Token_, electron_cutbasedID_decisions_veto);
  iEvent.getByToken(electron_cutbasedID_decisions_loose_Token_, electron_cutbasedID_decisions_loose);
  iEvent.getByToken(electron_cutbasedID_decisions_medium_Token_, electron_cutbasedID_decisions_medium);
  iEvent.getByToken(electron_cutbasedID_decisions_tight_Token_, electron_cutbasedID_decisions_tight);
  iEvent.getByToken(electron_mvaIsoID_decisions_wp80_Token_, electron_mvaIsoID_decisions_wp80);
  iEvent.getByToken(electron_mvaIsoID_decisions_wp90_Token_, electron_mvaIsoID_decisions_wp90);
  iEvent.getByToken(electron_mvaIsoID_decisions_wpHZZ_Token_, electron_mvaIsoID_decisions_wpHZZ);
  iEvent.getByToken(electron_mvaIsoID_decisions_wpLoose_Token_, electron_mvaIsoID_decisions_wpLoose);
  iEvent.getByToken(electron_mvaNoIsoID_decisions_wp80_Token_, electron_mvaNoIsoID_decisions_wp80);
  iEvent.getByToken(electron_mvaNoIsoID_decisions_wp90_Token_, electron_mvaNoIsoID_decisions_wp90);
  iEvent.getByToken(electron_mvaNoIsoID_decisions_wpLoose_Token_, electron_mvaNoIsoID_decisions_wpLoose);
  edm::Handle<vector<reco::Conversion> > conversions;
  iEvent.getByToken(conversionsToken_,conversions);
  edm::Handle<reco::BeamSpot> beamSpot;
  iEvent.getByToken(offlineBeamSpotToken_,beamSpot);


  for (uint i = 0; i < electrons->size(); ++i){
    const reco::GsfElectron ele = (*electrons)[i];
    reco::GsfElectronRef eleRef(electrons, i);
    if(ele.pt() < 5) continue;

    eleE.push_back(ele.energy());
    elePt.push_back(ele.pt());
    eleEta.push_back(ele.eta());
    elePhi.push_back(ele.phi());
    eleCharge.push_back(ele.charge());
    ele.hadronicOverEm();
    eleE_SC.push_back(ele.superCluster()->energy());
    eleEta_SC.push_back(ele.superCluster()->eta());
    elePhi_SC.push_back(ele.superCluster()->phi());

    eleSigmaIetaIeta.push_back(ele.sigmaIetaIeta());
    eleFull5x5SigmaIetaIeta.push_back(ele.full5x5_sigmaIetaIeta());
    eleR9.push_back(ele.r9());
    ele_dEta.push_back(ele.deltaEtaSuperClusterTrackAtVtx() - ele.superCluster()->eta() + ele.superCluster()->seed()->eta());

    ele_dPhi.push_back(ele.deltaPhiSuperClusterTrackAtVtx());
    ele_HoverE.push_back(ele.hcalOverEcal());
    ele_d0.push_back(ele.gsfTrack().get()->dxy(highestSumPt2Vertex.position()));
    ele_dZ.push_back(ele.gsfTrack().get()->dz(highestSumPt2Vertex.position()));

    ele_pileupIso.push_back(ele.pfIsolationVariables().sumPUPt);
    ele_chargedIso.push_back(ele.pfIsolationVariables().sumChargedHadronPt);
    ele_photonIso.push_back(ele.pfIsolationVariables().sumPhotonEt);
    ele_neutralHadIso.push_back(ele.pfIsolationVariables().sumNeutralHadronEt);
    ele_MissHits.push_back(ele.gsfTrack()->hitPattern().numberOfAllHits(reco::HitPattern::MISSING_INNER_HITS));
    ele_passCutBasedIDVeto.push_back((*electron_cutbasedID_decisions_veto)[eleRef]);
    ele_passCutBasedIDLoose.push_back((*electron_cutbasedID_decisions_loose)[eleRef]);
    ele_passCutBasedIDMedium.push_back((*electron_cutbasedID_decisions_medium)[eleRef]);
    ele_passCutBasedIDTight.push_back((*electron_cutbasedID_decisions_tight)[eleRef]);
    ele_passMVAIsoIDWP80.push_back((*electron_mvaIsoID_decisions_wp80)[eleRef]);
    ele_passMVAIsoIDWP90.push_back((*electron_mvaIsoID_decisions_wp90)[eleRef]);
    ele_passMVAIsoIDWPHZZ.push_back((*electron_mvaIsoID_decisions_wpHZZ)[eleRef]);
    ele_passMVAIsoIDWPLoose.push_back((*electron_mvaIsoID_decisions_wpLoose)[eleRef]);
    ele_passMVANoIsoIDWP80.push_back((*electron_mvaNoIsoID_decisions_wp80)[eleRef]);
    ele_passMVANoIsoIDWP90.push_back((*electron_mvaNoIsoID_decisions_wp90)[eleRef]);
    ele_passMVANoIsoIDWPLoose.push_back((*electron_mvaNoIsoID_decisions_wpLoose)[eleRef]);

    //---------------
    //Conversion Veto
    //---------------

    !ConversionTools::hasMatchedConversion(ele,(*conversions),beamSpot->position());

    if( beamSpot.isValid() && conversions.isValid() )
    {
      ele_PassConvVeto.push_back(!ConversionTools::hasMatchedConversion(ele,(*conversions),beamSpot->position()));
    } else {
      ele_PassConvVeto.push_back(false);
    }
    // // 1/E - 1/P
    if( ele.ecalEnergy() == 0 ){
      ele_OneOverEminusOneOverP.push_back(1e30);
    } else if( !std::isfinite(ele.ecalEnergy())){
      ele_OneOverEminusOneOverP.push_back(1e30);
    } else {
    ele_OneOverEminusOneOverP.push_back(1./ele.ecalEnergy()  -  ele.eSuperClusterOverP()/ele.ecalEnergy());
    }
  }
/////////////////////////
// add all muons
////////////////////////


unsigned int nMuons = 0;
std::vector<float> muonE;
std::vector<float> muonPt;
std::vector<float> muonEta;
std::vector<float> muonPhi;
std::vector<int> muonCharge;
std::vector<bool> muonIsLoose;
std::vector<bool> muonIsMedium;
std::vector<bool> muonIsTight;
std::vector<float> muon_d0;
std::vector<float> muon_d0Err;
std::vector<float> muon_dZ;
std::vector<float> muon_ip3d;
std::vector<float> muon_ip3dSignificance;
std::vector<unsigned int> muonType;
std::vector<unsigned int> muonQuality;
std::vector<float> muon_pileupIso;
std::vector<float> muon_chargedIso;
std::vector<float> muon_photonIso;
std::vector<float> muon_neutralHadIso;
std::vector<float> muon_validFractionTrackerHits;
std::vector<float> muon_normChi2;
std::vector<float> muon_chi2LocalPosition;
std::vector<float> muon_kinkFinder;
std::vector<float> muon_segmentCompatability;
std::vector<float> muon_trkIso;
std::vector<float> muon_tuneP_Pt;
std::vector<float> muon_tuneP_PtErr;
std::vector<float> muon_tuneP_Eta;
std::vector<float> muon_tuneP_Phi;
std::vector<int> muon_tuneP_MuonBestTrackType;
std::vector<bool> muon_isHighPtMuon;
std::vector<bool> muon_isTrackerHighPtMuon;

// read muon HLT filter names
for (int i = 0; i<MAX_MuonHLTFilters; ++i) muonHLTFilterNames[i] = "";
 ifstream myMuonHLTFilterFile (edm::FileInPath(muonHLTFilterNamesFile_.c_str()).fullPath().c_str()) ;
 if (myMuonHLTFilterFile.is_open()) {
   char tmp[1024];
   string line;
   int index;
   string hltfiltername;

   while(myMuonHLTFilterFile>>line) {

     if ( line.empty() || line.substr(0,1) == "#") {
       myMuonHLTFilterFile.getline(tmp,1024);
       continue;
     }
     index = atoi(line.c_str());
     myMuonHLTFilterFile >> hltfiltername;
     if (index < MAX_MuonHLTFilters) {
      muonHLTFilterNames[index] = hltfiltername;
     }
   }
   myMuonHLTFilterFile.close();
 } else cout << "ERROR!!! Could not open trigger path name file : " << edm::FileInPath(muonHLTFilterNamesFile_.c_str()).fullPath().c_str() << "\n";

for (unsigned int i = 0; i < muonColl.size(); i++) {
  const reco::Muon* mu = &(muonColl)[i];

    muonE.push_back(mu->energy());
    muonPt.push_back(mu->pt());
    muonEta.push_back(mu->eta());
    muonPhi.push_back(mu->phi());
    muonCharge.push_back(mu->charge());
    muonIsLoose.push_back(muon::isLooseMuon(*mu));
    muonIsMedium.push_back(muon::isMediumMuon(*mu));
    muonIsTight.push_back(muon::isTightMuon(*mu, highestSumPt2Vertex));
    muon_d0.push_back(-mu->muonBestTrack()->dxy(highestSumPt2Vertex.position()));
    muon_d0Err.push_back(-mu->muonBestTrack()->dxyError());
    muon_dZ.push_back(mu->muonBestTrack()->dz(highestSumPt2Vertex.position()));
    // muon_ip3d.push_back(mu->dB(pat::Muon::PV3D));
    // muon_ip3dSignificance.push_back(mu->dB(pat::Muon::PV3D)/mu->edB(pat::Muon::PV3D));
    muonType.push_back(mu->isMuon() + 2*mu->isGlobalMuon() + 4*mu->isTrackerMuon() + 8*mu->isStandAloneMuon()
      + 16*mu->isCaloMuon() + 32*mu->isPFMuon() + 64*mu->isRPCMuon());
    muonQuality.push_back(
      muon::isGoodMuon(*mu,muon::All)
      + pow(2,1)*muon::isGoodMuon(*mu,muon::AllGlobalMuons)
    + pow(2,2)*muon::isGoodMuon(*mu,muon::AllStandAloneMuons)
    + pow(2,3)*muon::isGoodMuon(*mu,muon::AllTrackerMuons)
    + pow(2,4)*muon::isGoodMuon(*mu,muon::TrackerMuonArbitrated)
    + pow(2,5)*muon::isGoodMuon(*mu,muon::AllArbitrated)
    + pow(2,6)*muon::isGoodMuon(*mu,muon::GlobalMuonPromptTight)
    + pow(2,7)*muon::isGoodMuon(*mu,muon::TMLastStationLoose)
    + pow(2,8)*muon::isGoodMuon(*mu,muon::TMLastStationTight)
    + pow(2,9)*muon::isGoodMuon(*mu,muon::TM2DCompatibilityLoose)
    + pow(2,10)*muon::isGoodMuon(*mu,muon::TM2DCompatibilityTight)
    + pow(2,11)*muon::isGoodMuon(*mu,muon::TMOneStationLoose)
    + pow(2,12)*muon::isGoodMuon(*mu,muon::TMOneStationTight)
    + pow(2,13)*muon::isGoodMuon(*mu,muon::TMLastStationOptimizedLowPtLoose)
    + pow(2,14)*muon::isGoodMuon(*mu,muon::TMLastStationOptimizedLowPtTight)
    + pow(2,15)*muon::isGoodMuon(*mu,muon::GMTkChiCompatibility)
    + pow(2,16)*muon::isGoodMuon(*mu,muon::GMStaChiCompatibility)
    + pow(2,17)*muon::isGoodMuon(*mu,muon::GMTkKinkTight)
    + pow(2,18)*muon::isGoodMuon(*mu,muon::TMLastStationAngLoose)
    + pow(2,19)*muon::isGoodMuon(*mu,muon::TMLastStationAngTight)
    + pow(2,20)*muon::isGoodMuon(*mu,muon::TMOneStationAngLoose)
    + pow(2,21)*muon::isGoodMuon(*mu,muon::TMOneStationAngTight)
    + pow(2,22)*muon::isGoodMuon(*mu,muon::TMLastStationOptimizedBarrelLowPtLoose)
    + pow(2,23)*muon::isGoodMuon(*mu,muon::TMLastStationOptimizedBarrelLowPtTight)
    + pow(2,24)*muon::isGoodMuon(*mu,muon::RPCMuLoose)
    //   //This is the soft muon ID
    + pow(2,25)*( muon::isGoodMuon(*mu,muon::TMOneStationTight)
		  && mu->innerTrack()->hitPattern().trackerLayersWithMeasurement() > 5
		  && mu->innerTrack()->hitPattern().pixelLayersWithMeasurement() > 0
		  && mu->innerTrack()->quality(reco::TrackBase::highPurity)
		  && fabs(mu->innerTrack()->dxy(highestSumPt2Vertex.position())) < 0.3
		  && fabs(mu->innerTrack()->dz(highestSumPt2Vertex.position())) < 20.
    ));

    muon_pileupIso.push_back(mu->pfIsolationR04().sumPUPt);
    muon_chargedIso.push_back(mu->pfIsolationR04().sumChargedHadronPt);
    muon_photonIso.push_back(mu->pfIsolationR04().sumPhotonEt);
    muon_neutralHadIso.push_back(mu->pfIsolationR04().sumNeutralHadronEt);
    muon_validFractionTrackerHits.push_back((mu->innerTrack().isNonnull() ? mu->track()->validFraction() : -99.0));
    muon_normChi2.push_back(( muon::isGoodMuon(*mu,muon::AllGlobalMuons) ? mu->globalTrack()->normalizedChi2() : -99.0));
    muon_chi2LocalPosition.push_back(mu->combinedQuality().chi2LocalPosition);
    muon_kinkFinder.push_back(mu->combinedQuality().trkKink);
    muon_segmentCompatability.push_back(muon::segmentCompatibility(*mu));

    muon_trkIso.push_back(mu->isolationR03().sumPt);


    // TuneP muon
    muon_tuneP_Pt.push_back(mu->tunePMuonBestTrack()->pt());
    muon_tuneP_PtErr.push_back(mu->tunePMuonBestTrack()->ptError());
    muon_tuneP_Eta.push_back(mu->tunePMuonBestTrack()->eta());
    muon_tuneP_Phi.push_back(mu->tunePMuonBestTrack()->phi());
    muon_tuneP_MuonBestTrackType.push_back(mu->tunePMuonBestTrackType());
    // case Muon::InnerTrack:
    //   return 0;
    // case Muon::OuterTrack:
    //   return 1;
    // case Muon::CombinedTrack:
    //   return 2;
    // case Muon::TPFMS:
    //   return 3;
    // case Muon::Picky:
    //   return 4;
    // case Muon::DYT:
    //   return 5;
    // case Muon::None:
    //   return 6;
    muon_isHighPtMuon.push_back(muon::isHighPtMuon(*mu, highestSumPt2Vertex));
    muon_isTrackerHighPtMuon.push_back(muon::isTrackerHighPtMuon(*mu, highestSumPt2Vertex));


    nMuons++;

}
std::vector<std::vector<float>> triggerObjectE;
std::vector<std::vector<float>> triggerObjectPt;
std::vector<std::vector<float>> triggerObjectEta;
std::vector<std::vector<float>> triggerObjectPhi;


// std::vector<TLorentzVector> trigObjP4s;
for ( int q=0; q<MAX_MuonHLTFilters;q++) {
  trigtools::getP4sOfObsPassingFilter(trigObjP4s,*trigEvent,muonHLTFilterNames[q].c_str(),"HLT");
  std::vector<float> triggerObjectE_temp;
  std::vector<float> triggerObjectPt_temp;
  std::vector<float> triggerObjectEta_temp;
  std::vector<float> triggerObjectPhi_temp;
  for(size_t objNr=0; objNr<trigObjP4s.size(); objNr++) {
    triggerObjectE_temp.push_back(trigObjP4s[objNr].Energy());
    triggerObjectPt_temp.push_back(trigObjP4s[objNr].Pt());
    triggerObjectEta_temp.push_back(trigObjP4s[objNr].Eta());
    triggerObjectPhi_temp.push_back(trigObjP4s[objNr].Phi());
  }
  triggerObjectE.push_back(triggerObjectE_temp);
  triggerObjectPt.push_back(triggerObjectPt_temp);
  triggerObjectEta.push_back(triggerObjectEta_temp);
  triggerObjectPhi.push_back(triggerObjectPhi_temp);
}

  float RecoPFMET = -10, RecoPFMET_phi = -10, RecoPFMET_sigf = -10, RecoPFMHT = -10;
  float HLTCaloMET = -10, HLTCaloMET_phi = -10, HLTCaloMET_sigf = -10;
  float HLTCaloMETClean = -10, HLTCaloMETClean_phi = -10, HLTCaloMETClean_sigf = -10;
  float HLTCaloMHT = -10, HLTCaloMHT_phi = -10, HLTCaloMHT_sigf = -10;
  float HLTPFMET = -10, HLTPFMET_phi = -10, HLTPFMET_sigf = -10;
  float HLTPFMHT = -10, HLTPFMHT_phi = -10, HLTPFMHT_sigf = -10;

  //===================== Handle For RecoCaloMET ===================
  if (recoCaloMETHandle.isValid() && !recoCaloMETHandle->empty()) {
    for (unsigned int i = 0; i < recoCaloMETHandle->size(); i++) {
      const reco::CaloMET* recoCaloMet = &(*recoCaloMETHandle)[i];
      RecoCaloMET = recoCaloMet->et();
      RecoCaloMET_phi = recoCaloMet->phi();
      RecoCaloMET_sigf = recoCaloMet->significance();
    }
  }

  //===================== Handle For RecoPFMET ===================
  if (recoPFMETHandle.isValid() && !recoPFMETHandle->empty()) {
    for (unsigned int i = 0; i < recoPFMETHandle->size(); i++) {
      const reco::PFMET* recoPFMet = &(*recoPFMETHandle)[i];
      RecoPFMET = recoPFMet->et();
      RecoPFMET_phi = recoPFMet->phi();
      RecoPFMET_sigf = recoPFMet->significance();
    }
  }

  //===================== Handle For HLT Trigger Summary ===================
  const edm::Handle<trigger::TriggerEvent> hltTriggerSummaryHandle = iEvent.getHandle(trigEventToken_);
  if (hltTriggerSummaryHandle.isValid()) {

    int caloMETKey = 0, caloMETCleanKey = 0, caloMHTKey = 0, pfMHTKey = 0, pfMETKey = 0;
    // loop over trigger object collections to find HLT CaloMET, CaloMETClean, CaloMHT, PFMHT, PFMET collections
    for (int iC = 0; iC < hltTriggerSummaryHandle->sizeCollections(); iC++) {
      if(hltTriggerSummaryHandle->collectionTag(iC).encode()=="hltMet::HLT") {
        // collectionKey(iC) gives trigger object key ONE PAST the object collection of interest
        caloMETKey = hltTriggerSummaryHandle->collectionKey(iC);
        // HLT MET object collections ALWAYS have four objects {MET, TET, MET significance, ELongitudinal}, hence -4 for MET value
        HLTCaloMET = hltTriggerSummaryHandle->getObjects()[caloMETKey-4].pt();
        HLTCaloMET_phi = hltTriggerSummaryHandle->getObjects()[caloMETKey-4].phi();
        // and -2 for MET significance
        // significance  saved as .pt() but obviously pt holds no meaning here
        HLTCaloMET_sigf = hltTriggerSummaryHandle->getObjects()[caloMETKey-2].pt();
      } if(hltTriggerSummaryHandle->collectionTag(iC).encode()=="hltMetClean::HLT") {
          caloMETCleanKey = hltTriggerSummaryHandle->collectionKey(iC);
          HLTCaloMETClean = hltTriggerSummaryHandle->getObjects()[caloMETCleanKey-4].pt();
          HLTCaloMETClean_phi = hltTriggerSummaryHandle->getObjects()[caloMETCleanKey-4].phi();
          HLTCaloMETClean_sigf = hltTriggerSummaryHandle->getObjects()[caloMETCleanKey-2].pt();
      } if(hltTriggerSummaryHandle->collectionTag(iC).encode()=="hltMht::HLT") {
          caloMHTKey = hltTriggerSummaryHandle->collectionKey(iC);
          // HLT MHT object collections ALWAYS have four objects {MHT, THT, MHT significance, HLongitudinal}, hence -4 for MHT value
          HLTCaloMHT = hltTriggerSummaryHandle->getObjects()[caloMHTKey-4].pt();
          HLTCaloMHT_phi = hltTriggerSummaryHandle->getObjects()[caloMHTKey-4].phi();
          // and -2 for MHT significance
          // significance  saved as .pt() but obviously pt holds no meaning here
          HLTCaloMHT_sigf = hltTriggerSummaryHandle->getObjects()[caloMHTKey-2].pt();
      } if(hltTriggerSummaryHandle->collectionTag(iC).encode()=="hltPFMHTTightID::HLT") {
        pfMHTKey = hltTriggerSummaryHandle->collectionKey(iC);
        HLTPFMHT = hltTriggerSummaryHandle->getObjects()[pfMHTKey-4].pt();
        HLTPFMHT_phi = hltTriggerSummaryHandle->getObjects()[pfMHTKey-4].phi();
        HLTPFMHT_sigf = hltTriggerSummaryHandle->getObjects()[pfMHTKey-2].pt();
      } if(hltTriggerSummaryHandle->collectionTag(iC).encode()=="hltPFMETProducer::HLT") {
        pfMETKey = hltTriggerSummaryHandle->collectionKey(iC);
        HLTPFMET = hltTriggerSummaryHandle->getObjects()[pfMETKey-4].pt();
        HLTPFMET_phi = hltTriggerSummaryHandle->getObjects()[pfMETKey-4].phi();
        HLTPFMET_sigf = hltTriggerSummaryHandle->getObjects()[pfMETKey-2].pt();
      }
    }
  }

  // PF jet info for the ntuple
  unsigned int pfJetsNum = 0;
  float pfJetHT = 0.;

  std::vector<float> Jets_pt;
  std::vector<float> Jets_eta;
  std::vector<float> Jets_phi;
  std::vector<float> Jets_mass;
  std::vector<float> Jets_E;
  std::vector<float> Jets_pdgId;
  std::vector<float> Jets_et;
  std::vector<float> Jets_chargedEmEnergyFraction;
  std::vector<float> Jets_neutralEmEnergyFraction;
  std::vector<float> Jets_chargedHadronEnergyFraction;
  std::vector<float> Jets_neutralHadronEnergyFraction;
  std::vector<float> Jets_muonEnergyFraction;
  std::vector<int> Jets_chargedMultiplicity;
  std::vector<int> Jets_neutralMultiplicity;
  std::vector<float> Jets_jetArea;
  std::vector<float> Jets_pileupE;
  // std::vector<float> Jets_pileupId;
  // std::vector<int> Jets_pileupIdFlag;



  // Loop on pfJetColl for the ntuple, for histos there is a matching to the HSCP candidates
  if (pfJetHandle.isValid() && !pfJetHandle->empty()) {
    const reco::PFJetCollection* pfJetColl = pfJetHandle.product();
    TLorentzVector pMHT;
    for (unsigned int i = 0; i < pfJetColl->size(); i++) {
      const reco::PFJet* jet = &(*pfJetColl)[i];
      if (jet->pt() < 30) {
        continue;
      }
      pfJetsNum++;
      Jets_pt.push_back(jet->pt());
      Jets_eta.push_back(jet->eta());
      Jets_phi.push_back(jet->phi());
      Jets_mass.push_back(jet->mass());
      Jets_E.push_back(jet->energy());
      Jets_pdgId.push_back(jet->pdgId());
      Jets_et.push_back(jet->et());
      Jets_chargedEmEnergyFraction.push_back(jet->chargedEmEnergyFraction());
      Jets_neutralEmEnergyFraction.push_back(jet->neutralEmEnergyFraction());
      Jets_chargedHadronEnergyFraction.push_back(jet->chargedHadronEnergyFraction());
      Jets_neutralHadronEnergyFraction.push_back(jet->neutralHadronEnergyFraction());
      Jets_muonEnergyFraction.push_back(jet->muonEnergyFraction());
      Jets_chargedMultiplicity.push_back(jet->chargedMultiplicity());
      Jets_neutralMultiplicity.push_back(jet->neutralMultiplicity());

      Jets_jetArea.push_back(jet->jetArea());
      Jets_pileupE.push_back(jet->pileup());
      // Jets_pileupId.push_back(jet->userFloat("pileupJetId:fullDiscriminant"));
      // Jets_pileupIdFlag.push_back(jet->userInt("pileupJetId:fullId")); //A bit map for loose, medium, and tight working points

      pfJetHT += jet->pt();
      TLorentzVector p4(jet->pt() * cos(jet->phi()), jet->pt() * sin(jet->phi()), 0, jet->et());
      pMHT += p4;
    }
    RecoPFMHT = pMHT.Pt();
  }



  //reinitialize the bookeeping array for each event
  for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
    HSCPTk[CutIndex] = false;
    HSCPTk_SystP[CutIndex] = false;
    HSCPTk_SystI[CutIndex] = false;
    HSCPTk_SystT[CutIndex] = false;
    HSCPTk_SystM[CutIndex] = false;
    HSCPTk_SystPU[CutIndex] = false;
    HSCPTk_SystHUp[CutIndex] = false;
    HSCPTk_SystHDown[CutIndex] = false;
    MaxMass[CutIndex] = -1;
    MaxMass_SystP[CutIndex] = -1;
    MaxMass_SystI[CutIndex] = -1;
    MaxMass_SystT[CutIndex] = -1;
    MaxMass_SystM[CutIndex] = -1;
    MaxMass_SystPU[CutIndex] = -1;
    MaxMass_SystHUp[CutIndex] = -1;
    MaxMass_SystHDown[CutIndex] = -1;
  }

  //load all event collection that will be used later on (HSCP, dEdx and TOF)
  unsigned int HSCP_count = 0;

  std::vector<float> HSCP_mT;
  std::vector<bool> HSCP_passCutPt55;
  std::vector<bool> HSCP_passPreselection;
  std::vector<bool> HSCP_passPreselectionSept8;
  std::vector<bool> HSCP_passSelection;
  std::vector<bool> HSCP_isPFMuon;
  std::vector<bool> HSCP_PFMuonPt;
  std::vector<float> HSCP_Charge;
  std::vector<float> HSCP_Pt;
  std::vector<float> HSCP_PtErr;
  std::vector<float> HSCP_Is_StripOnly;
  std::vector<float> HSCP_Ias;
  std::vector<float> HSCP_Ias_noPix_noTIB_noTID_no3TEC;
  std::vector<float> HSCP_Ias_PixelOnly;
  std::vector<float> HSCP_Ias_StripOnly;
  std::vector<float> HSCP_Ias_PixelOnly_noL1;
  std::vector<float> HSCP_Ih;
  std::vector<float> HSCP_Ick;  //return (Ih-C)/K
  std::vector<float> HSCP_Fmip;
  std::vector<float> HSCP_ProbXY;
  std::vector<float> HSCP_ProbXY_noL1;
  std::vector<float> HSCP_ProbQ;
  std::vector<float> HSCP_ProbQ_noL1;
  std::vector<float> HSCP_ProbQ_dEdx;
  std::vector<float> HSCP_Ndof;
  std::vector<float> HSCP_Chi2;
  std::vector<int>   HSCP_QualityMask;
  std::vector<bool>  HSCP_isHighPurity;
  std::vector<float> HSCP_EoverP;
  std::vector<bool>  HSCP_isMuon;
  std::vector<bool>  HSCP_isPhoton;
  std::vector<bool>  HSCP_isElectron;
  std::vector<bool>  HSCP_isChHadron;
  std::vector<bool>  HSCP_isNeutHadron;
  std::vector<bool>  HSCP_isPfTrack;
  std::vector<bool>  HSCP_isUndefined;
  std::vector<float> HSCP_ECAL_energy;
  std::vector<float> HSCP_HCAL_energy;
  std::vector<float> HSCP_TOF;
  std::vector<float> HSCP_TOFErr;
  std::vector<unsigned int> HSCP_TOF_ndof;
  std::vector<float> HSCP_DTTOF;
  std::vector<float> HSCP_DTTOFErr;
  std::vector<unsigned int> HSCP_DTTOF_ndof;
  std::vector<float> HSCP_CSCTOF;
  std::vector<float> HSCP_CSCTOFErr;
  std::vector<unsigned int> HSCP_CSCTOF_ndof;
  std::vector<float> HSCP_Mass;
  std::vector<float> HSCP_MassErr;
  std::vector<float> HSCP_dZ;
  std::vector<float> HSCP_dXY;
  std::vector<float> HSCP_dR;
  std::vector<float> HSCP_p;
  std::vector<float> HSCP_eta;
  std::vector<float> HSCP_phi;
    // Number of (valid) track pixel+strip hits
  std::vector<unsigned int> HSCP_NOH;
    // Number of (valid) track pixel hits
  std::vector<unsigned int> HSCP_NOPH;
    // Fraction of valid track hits
  std::vector<float> HSCP_FOVH;
    // Number of missing hits from IP till last hit (excluding hits behind the last hit)
  std::vector<unsigned int> HSCP_NOMH;
    // Fraction of valid hits divided by total expected hits until the last one
  std::vector<float> HSCP_FOVHD;
    // Number of dEdx hits (= #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used)
  std::vector<unsigned int> HSCP_NOM;
  std::vector<float> HSCP_matchTrigMuon_minDeltaR;
  std::vector<float> HSCP_matchTrigMuon_pT;

  std::vector<float> HSCP_iso_TK;
  std::vector<float> HSCP_iso_ECAL;
  std::vector<float> HSCP_iso_HCAL;
  std::vector<float> HSCP_track_genTrackMiniIsoSumPt;
  std::vector<float> HSCP_PFMiniIso_relative;
  std::vector<float> HSCP_PFMiniIso_wMuon_relative;
  std::vector<float> HSCP_track_PFIsolationR005_sumChargedHadronPt;
  std::vector<float> HSCP_track_PFIsolationR005_sumNeutralHadronPt;
  std::vector<float> HSCP_track_PFIsolationR005_sumPhotonPt;
  std::vector<float> HSCP_track_PFIsolationR005_sumPUPt;
  std::vector<float> HSCP_track_PFIsolationR01_sumChargedHadronPt;
  std::vector<float> HSCP_track_PFIsolationR01_sumNeutralHadronPt;
  std::vector<float> HSCP_track_PFIsolationR01_sumPhotonPt;
  std::vector<float> HSCP_track_PFIsolationR01_sumPUPt;
  std::vector<float> HSCP_track_PFIsolationR03_sumChargedHadronPt;
  std::vector<float> HSCP_track_PFIsolationR03_sumNeutralHadronPt;
  std::vector<float> HSCP_track_PFIsolationR03_sumPhotonPt;
  std::vector<float> HSCP_track_PFIsolationR03_sumPUPt;
  std::vector<float> HSCP_track_PFIsolationR05_sumChargedHadronPt;
  std::vector<float> HSCP_track_PFIsolationR05_sumNeutralHadronPt;
  std::vector<float> HSCP_track_PFIsolationR05_sumPhotonPt;
  std::vector<float> HSCP_track_PFIsolationR05_sumPUPt;
  std::vector<float> HSCP_muon_PFIsolationR03_sumChargedHadronPt;
  std::vector<float> HSCP_muon_PFIsolationR03_sumNeutralHadronPt;
  std::vector<float> HSCP_muon_PFIsolationR03_sumPhotonPt;
  std::vector<float> HSCP_muon_PFIsolationR03_sumPUPt;
  std::vector<float> HSCP_Ih_noL1;
  std::vector<float> HSCP_Ih_15drop;
  std::vector<float> HSCP_Ih_StripOnly;
  std::vector<float> HSCP_Ih_StripOnly_15drop;
  std::vector<float> HSCP_Ih_PixelOnly_noL1;
  std::vector<float> HSCP_Ih_SaturationCorrectionFromFits;
  std::vector<std::vector<float>> HSCP_clust_charge;  //dedx charge -> either strip or pixel
  std::vector<std::vector<float>> HSCP_clust_pathlength;
  std::vector<std::vector<unsigned int>> HSCP_clust_nstrip;
  std::vector<std::vector<bool>> HSCP_clust_sat254;
  std::vector<std::vector<bool>> HSCP_clust_sat255;
  std::vector<std::vector<uint32_t>> HSCP_clust_detid;
  std::vector<std::vector<bool>> HSCP_clust_isStrip;  //is it a SiStrip cluster?
  std::vector<std::vector<bool>> HSCP_clust_isPixel;  //is it a Pixel hit?
  std::vector<float> HSCP_GenId;
  std::vector<float> HSCP_GenCharge;
  std::vector<float> HSCP_GenMass;
  std::vector<float> HSCP_GenPt;
  std::vector<float> HSCP_GenEta;
  std::vector<float> HSCP_GenPhi;
  std::vector<float> HSCP_tuneP_Pt;
  std::vector<float> HSCP_tuneP_PtErr;
  std::vector<float> HSCP_tuneP_Eta;
  std::vector<float> HSCP_tuneP_Phi;
  std::vector<int> HSCP_tuneP_MuonBestTrackType;

  std::vector<int> HSCP_ErrorHisto_bin;
  std::vector<int> HSCP_type;

  float dr_min_hlt_muon = 9999.0;

  //====================loop over HSCP candidates===================
  if (debug_ > 0 && trigInfo_ > 0) LogPrint(MOD) << "Loop over HSCP candidates:";
  unsigned int candidate_count = 0;
  unsigned int postPreS_candidate_count = 0;
  int bestCandidateIndex = -1;
  float maxIhSoFar = -1.;
  float bestCandidateIas = -1.;
  float bestCandidateProbQNoL1 = -1.;
  float bestCandidateDrMinHltMuon = 9999.;
  float bestCandidateGenBeta = -1.;
  unsigned int bestSoFarCandCutInd = 0;
  bool passTechnicalChecks = false;
  
  tuple->EventCutFlow->Fill(0.0, EventWeight_);
  
  for (const auto& hscp : iEvent.get(hscpToken_)) {
    // Number of tracks before any trigger or preselection
    tuple->CutFlow->Fill(0., EventWeight_);
    if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "  --------------------------------------------";
    
    candidate_count++;
    int ErrorHisto_bin = 0;
    if (trigInfo_ > 0) {
      // First bin of the error histo is all tracks (after trigger)
      tuple->ErrorHisto->Fill(0.);
      
      if ( hscp.type() == susybsm::HSCParticleType::globalMuon) {
        tuple->BefPreS_RecoHSCParticleType->Fill(0.);
      } else if ( hscp.type() == susybsm::HSCParticleType::trackerMuon) {
        tuple->BefPreS_RecoHSCParticleType->Fill(1.);
      } else if ( hscp.type() == susybsm::HSCParticleType::matchedStandAloneMuon) {
        tuple->BefPreS_RecoHSCParticleType->Fill(2.);
      } else if ( hscp.type() == susybsm::HSCParticleType::standAloneMuon) {
        tuple->BefPreS_RecoHSCParticleType->Fill(3.);
      } else if ( hscp.type() == susybsm::HSCParticleType::innerTrack) {
        tuple->BefPreS_RecoHSCParticleType->Fill(4.);
      } else if ( hscp.type() == susybsm::HSCParticleType::unknown) {
        tuple->BefPreS_RecoHSCParticleType->Fill(5.);
      }
    }
    
    if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "  >> This is HSCP candidate track " << candidate_count;
    // Tracker only analysis must have either a track in the tracking system

    if (typeMode_ == 0)  {
      if ((hscp.type() != susybsm::HSCParticleType::trackerMuon) && (hscp.type() != susybsm::HSCParticleType::globalMuon) && (hscp.type() != susybsm::HSCParticleType::innerTrack)) {
        if (debug_ > 0 && trigInfo_ > 0) LogPrint(MOD) << "  >> Tracker only analysis w/o a tracker track";
        // Second bin of the error histo, num tracks that fail the track existence checks
        if (trigInfo_ > 0) {
          tuple->ErrorHisto->Fill(1.);
        }
        continue;
      }
    }
    
    // Tracker + muon analysis must have a muon, either a tracker muon or a global muon
    if (typeMode_ == 1)  {
      if ((hscp.type() != susybsm::HSCParticleType::trackerMuon) && (hscp.type() != susybsm::HSCParticleType::globalMuon)) {
        if (debug_ > 0 && trigInfo_ > 0) LogPrint(MOD) << "  >> Tracker + Muon analysis w/o a tracker muon or a global muon";
          // Second bin of the error histo, num tracks that fail the track existence checks
        if (trigInfo_ > 0) {
          tuple->ErrorHisto->Fill(1.);
        }
        continue;
      }
    }

    // Define muon reference
    // For TOF only analysis we must have a muon connected to the HSCP candidate
    reco::MuonRef muon = hscp.muonRef();
    if (typeMode_ == 3 && muon.isNull()) {
      if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "  >> TOF only mode but no muon connected to the candidate -- skipping it";
      // Second bin of the error histo, num tracks that fail the track existence checks
      if (trigInfo_ > 0) tuple->ErrorHisto->Fill(1.);
      continue;
    }

    // Define track reference
    // For TOF only analysis use updated stand alone muon track, otherwise use inner tracker track
    reco::TrackRef track = (typeMode_ != 3) ? hscp.trackRef() : muon->standAloneMuon();
    if (track.isNull()) {
      if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "  >> Event has no track associated to this HSCP, skipping it";
      // Third bin of the error histo, no tracks
      if (trigInfo_ > 0) tuple->ErrorHisto->Fill(2.);
      continue;
    }

    float trackerPt = track->pt();
    float muonPt = 0.f;
    if (!muon.isNull()) {
      reco::TrackRef trackMu;
      if (muon->tunePMuonBestTrack().isNonnull()) {
        trackMu = muon->tunePMuonBestTrack();
        muonPt = trackMu->pt();

      }
    }
    if (trigInfo_ > 0) {
      tuple->BefPreS_MuonPtVsTrackPt->Fill(muonPt, trackerPt, EventWeight_);
      tuple->BefPreS_RelDiffMuonPtAndTrackPt->Fill((muonPt-trackerPt)/trackerPt, EventWeight_);
    }
    
    // Require a track segment in the muon system
    if (typeMode_ > 1 && typeMode_ != 5 && (muon.isNull() || !muon->isStandAloneMuon())) {
      if (debug_> 0) LogPrint(MOD) << "  >> typeMode_ > 1 && typeMode_ != 5 && (muon.isNull() || !muon->isStandAloneMuon()), skipping it";
      // Second bin of the error histo, num tracks that fail the track existence checks
      if (trigInfo_ > 0) tuple->ErrorHisto->Fill(1.);
      continue;
    }

    // Apply a scale factor to muon only analysis to account for differences seen in data/MC preselection efficiency
    // For eta regions where Data > MC no correction to be conservative
    if (!isData && typeMode_ == 3 && scaleFactor(track->eta()) < RNG->Uniform(0, 1)) {
      if (debug_> 0) LogPrint(MOD) << "  >> This is a non-data but TOF onlypwd mode where the eta scale factor is non-uniform, skipping it";
      continue;
    }

    if (vertexColl.size() < 1) {
        if (debug_> 0) LogPrint(MOD) << "  >> Event has no primary vertices, skipping it";
        // 4-th bin of the error histo, no PV
      if (trigInfo_ > 0) tuple->ErrorHisto->Fill(3.);
        continue;
    }


    // tune P muon
    float tuneP_Pt = -999;
    float tuneP_PtErr = -999;
    float tuneP_Eta = -999;
    float tuneP_Phi = -999;
    int tunePMuonBestTrackType = -999;
    if (!muon.isNull()) {
      tuneP_Pt = muon->tunePMuonBestTrack()->pt();
      tuneP_PtErr = muon->tunePMuonBestTrack()->ptError();
      tuneP_Eta = muon->tunePMuonBestTrack()->eta();
      tuneP_Phi = muon->tunePMuonBestTrack()->phi();
      tunePMuonBestTrackType = muon->tunePMuonBestTrackType();
    }


    // Reco - GEN track matching
    // For signal only, make sure that the candidate is associated to a true HSCP
    int closestGenIndex = -1;
    float dRMinGen = 9999.0;
    float dPtMinGen = 9999.0;
    unsigned int closestHSCPsPDGsID = 0;
    if (!isData) {
      if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "  >> Background MC, Reco - GEN track matching";
      if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "  >> GEN signal PDG IDs in the event: ";
      if (debug_> 0 && trigInfo_ > 0) std::cout << "     ";
      for (unsigned int g = 0; g < genColl.size(); g++) {
        if (isSignal && !isHSCPgenID(genColl[g])) {
          continue;
        }

        if (genColl[g].pt() < 5) {
          continue;
        }
        if (genColl[g].status() != 1) {
          continue;
        }
        
        if (isSignal && trigInfo_ > 0) {
          if (debug_> 0) std::cout << genColl[g].pdgId() << " , " ;
        }

        float dr = deltaR(genColl[g].eta(),genColl[g].phi(),track->eta(),track->phi());
        float dPt = (fabs(genColl[g].pt() - track->pt()))/track->pt();
        // We need to make this gen charge dependent for TauPrime2e
//        cout << "genColl[g].charge(): " << genColl[g].charge() << endl;
        if (dr < dRMinGen) {
          dRMinGen = dr;
          closestGenIndex = g;
        }
        if (dPt < dPtMinGen) {
          dPtMinGen = dPt;
        }
      }
      if (debug_> 0 && trigInfo_ > 0) cout << endl;
    }
    
    if (!isData && closestGenIndex < 0 ) {
      if (debug_ > 4 && trigInfo_ > 0) {
        LogPrint(MOD) << "min dr: " << dRMinGen << endl;
        LogPrint(MOD) << "min dPt: " << dPtMinGen << endl;
        LogPrint(MOD) << "  >> Event where we didnt find the gen candidate";
      }
      // 5-th bin of the error histo, didnt find the gen canidate
      if (trigInfo_ > 0) {
        tuple->ErrorHisto->Fill(4.);
      }
    }
    
    int genPdgId =  (closestGenIndex > 0) ?  genColl[closestGenIndex].pdgId() : 0;
    float genPt = (closestGenIndex > 0) ? genColl[closestGenIndex].pt() : -1.;
    float genEta = (closestGenIndex > 0) ? genColl[closestGenIndex].eta() : -9999.;
    float genPhi = (closestGenIndex > 0) ? genColl[closestGenIndex].phi() : -9999.;
    float genGammaBeta = (closestGenIndex > 0) ? genColl[closestGenIndex].p() /  genColl[closestGenIndex].mass() : -1.;
    float genBeta = (closestGenIndex > 0) ? genColl[closestGenIndex].p() / genColl[closestGenIndex].energy() : -1.f;
  

    if (trigInfo_ > 0 && !isData) {
      if (debug_ > 5) {
        LogPrint(MOD) << "  >> The min Gen candidate distance is " << dRMinGen << " for PDG ID " << genPdgId << " with pT " << genPt << " and eta " << genEta ;
      }
      tuple->BefPreS_GendRMin->Fill(dRMinGen);
      tuple->BefPreS_GenPtVsdRMinGen->Fill(genPt, dRMinGen);
      tuple->BefPreS_MuonPtOverGenPtVsTrackPtOverGenPt->Fill(muonPt/genPt,trackerPt/genPt, EventWeight_);
      // There used to be a cut on dRMinGen 0.015
      tuple->BefPreS_GenPtVsdRMinGenPostCut->Fill(genPt, dRMinGen);
      tuple->BefPreS_GenPtVsGenMinPt->Fill(genPt, dPtMinGen);
      // 2D plot to compare gen pt vs reco pt
      tuple->BefPreS_GenPtVsRecoPt->Fill(genPt, track->pt());
    }

    // ID for the candidate, it's mother, and it's nearest sibling, and their angle
    // the pt of the candidate and the number of siblings
    float closestBackgroundPDGsIDs[8] = {0.,0.,0.,9999.,9999.,0.,9999.,0.};
    // Look at the properties of the closes gen candidate
    if (isSignal) {
      closestHSCPsPDGsID = abs(genPdgId);
      // All HSCP candidates
      tuple->Gen_HSCPCandidateType->Fill(0., EventWeight_);
      if (trigInfo_ > 0) {
        tuple->BefPreS_HSCPCandidateType->Fill(0., EventWeight_);
      }
      // Neutral HSCP candidates
      if (   closestHSCPsPDGsID == 1000993 || closestHSCPsPDGsID == 1009113
          || closestHSCPsPDGsID == 1009223 || closestHSCPsPDGsID == 1009313
          || closestHSCPsPDGsID == 1009333 || closestHSCPsPDGsID == 1092114
          || closestHSCPsPDGsID == 1093214 || closestHSCPsPDGsID == 1093324
          || closestHSCPsPDGsID == 1000622 || closestHSCPsPDGsID == 1000642
          || closestHSCPsPDGsID == 1006113 || closestHSCPsPDGsID == 1006311
          || closestHSCPsPDGsID == 1006313 || closestHSCPsPDGsID == 1006333) {
        tuple->Gen_HSCPCandidateType->Fill(1., EventWeight_);
        if (trigInfo_ > 0) {
          tuple->BefPreS_HSCPCandidateType->Fill(1., EventWeight_);
        }
      }
      // Single-charged HSCP
      else if (   closestHSCPsPDGsID == 1009213 || closestHSCPsPDGsID == 1009323
               || closestHSCPsPDGsID == 1091114 || closestHSCPsPDGsID == 1092214
               || closestHSCPsPDGsID == 1093114 || closestHSCPsPDGsID == 1093224
               || closestHSCPsPDGsID == 1093314 || closestHSCPsPDGsID == 1093334
               || closestHSCPsPDGsID == 1000612 || closestHSCPsPDGsID == 1000632
               || closestHSCPsPDGsID == 1000652 || closestHSCPsPDGsID == 1006211
               || closestHSCPsPDGsID == 1006213 || closestHSCPsPDGsID == 1006321
               || closestHSCPsPDGsID == 1006323 || closestHSCPsPDGsID == 1000015) {
        tuple->Gen_HSCPCandidateType->Fill(2., EventWeight_);
        if (trigInfo_ > 0) {
          tuple->BefPreS_HSCPCandidateType->Fill(2., EventWeight_);
        }
      }
      // Double-charged R-hadrons
      else if (closestHSCPsPDGsID == 1092224 || closestHSCPsPDGsID == 1006223) {
        tuple->Gen_HSCPCandidateType->Fill(3., EventWeight_);
        if (trigInfo_ > 0) {
          tuple->BefPreS_HSCPCandidateType->Fill(3., EventWeight_);
        }
        // Dont mix double charged R-hadrons with the rest
        // The reco pt of them is 1/2 the pt of the gen track
//        continue;
      }
      // tau prime, could be single or multiple charged
      else if (closestHSCPsPDGsID == 17) {
        tuple->Gen_HSCPCandidateType->Fill(4., EventWeight_);
        if (trigInfo_ > 0) {
          tuple->BefPreS_HSCPCandidateType->Fill(4., EventWeight_);
        }
      }
      else {
        tuple->Gen_HSCPCandidateType->Fill(5., EventWeight_);
        if (trigInfo_ > 0) {
          tuple->Gen_HSCPCandidateType->Fill(5., EventWeight_);
        }
      }
    }

    // Match candidate track to HLT muon
    float hlt_match_pt = -9999.0;
    float temp_dr = 9999.;
    for(size_t objNr=0; objNr<trigObjP4s.size(); objNr++) {
      temp_dr = deltaR(trigObjP4s[objNr].Eta(),trigObjP4s[objNr].Phi(), track->eta(), track->phi());
      if (temp_dr < dr_min_hlt_muon)
      {
        dr_min_hlt_muon = temp_dr;
        hlt_match_pt = trigObjP4s[objNr].Pt();
      }
    }
    
    // Compute transverse mass mT between HSCP with and MET
    float massT = -10.;
    if (RecoPFMET > 0) massT = sqrt(2*track->pt()*RecoPFMET*(1-cos(track->phi()-RecoPFMET_phi)));
    
    if (highestSumPt2VertexIndex < 0) {
      if (debug_> 0) LogPrint(MOD) << "  >> PV associated to this track has no **good** primary vertex match, skipping it";
      // 4-th bin of the error histo, no PV
      if (trigInfo_ > 0) tuple->ErrorHisto->Fill(3.);
      continue;
    }

    // Impact paramters dz and dxy
    float dz = track->dz(highestSumPt2Vertex.position());
    float dxy = track->dxy(highestSumPt2Vertex.position());
    // Impact parameters if no vertex
    float dzFromBeamSpot = 9999.;
    float dxyFromBeamSpot = 9999.;

    //For TOF only analysis match to a SA track without vertex constraint for IP cuts
    if (typeMode_ == 3) {
      //Find closest NV track
      const std::vector<reco::Track> noVertexTrackColl = iEvent.get(refittedStandAloneMuonsToken_);
      reco::Track NVTrack;
      float minDr = 15;
      for (unsigned int i = 0; i < noVertexTrackColl.size(); i++) {
        auto dR = deltaR(track->eta(), track->phi(), noVertexTrackColl[i].eta(), noVertexTrackColl[i].phi());
        if (dR < minDr) {
          minDr = dR;
          NVTrack = noVertexTrackColl[i];
        }
      }
      if (trigInfo_ > 0) {
      tuple->BefPreS_dR_NVTrack->Fill(minDr, EventWeight_);
      }
      if (minDr < 0.4) {
        // Find displacement of tracks with respect to beam spot
        const reco::BeamSpot beamSpotColl = iEvent.get(offlineBeamSpotToken_);
        dzFromBeamSpot = NVTrack.dz(beamSpotColl.position());
        dxyFromBeamSpot = NVTrack.dxy(beamSpotColl.position());
      }

//      if (muonStations(NVTrack.hitPattern()) < minMuStations_) return false;
    } // End condition for TOF only analysis

    // Save PF informations and isolation
    float track_PFIso005_sumCharHadPt = 0, track_PFIso005_sumNeutHadPt = 0, track_PFIso005_sumPhotonPt = 0, track_PFIso005_sumPUPt = 0;
    float track_PFIso01_sumCharHadPt = 0, track_PFIso01_sumNeutHadPt = 0, track_PFIso01_sumPhotonPt = 0, track_PFIso01_sumPUPt = 0;
    float track_PFIso03_sumCharHadPt = 0, track_PFIso03_sumNeutHadPt = 0, track_PFIso03_sumPhotonPt = 0, track_PFIso03_sumPUPt = 0;
    float track_PFIso05_sumCharHadPt = 0, track_PFIso05_sumNeutHadPt = 0, track_PFIso05_sumPhotonPt = 0, track_PFIso05_sumPUPt = 0;
    float track_PFMiniIso_sumCharHadPt = 0, track_PFMiniIso_sumNeutHadPt = 0, track_PFMiniIso_sumPhotonPt = 0, track_PFMiniIso_sumPUPt = 0, track_PFMiniIso_sumMuonPt = 0;
    float pf_energy=0, pf_ecal_energy = 0, pf_hcal_energy = 0;

    // loop on PF Jets for the histograms
    float dRMinPfJetWithCuts = 9999.0;
    float dRMinPfJetWithNoCuts = 9999.0;
    float closestPfJetMuonFraction = 0.0;
    float closestPfJetElectronFraction = 0.0;
    float closestPfJetPhotonFraction = 0.0;

    // This is again repeated in the preselection
    bool pf_isMuon = false, pf_isElectron = false, pf_isChHadron = false, pf_isNeutHadron = false;

    bool pf_isPfTrack = false,  pf_isPhoton = false, pf_isUndefined = false;
    float track_PFMiniIso_sumLeptonPt = 0;
    float track_PFMiniIso_otherPt = 0;
    if(pfCandHandle.isValid() && !pfCandHandle->empty()) {
      const reco::PFCandidateCollection* pf = pfCandHandle.product();
      for (unsigned int i = 0; i < pf->size(); i++){
          // https://github.com/cms-sw/cmssw/blob/72d0fc00976da53d1fb745eb7f37b2a4ad965d7e/
          // PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc#L555
        const reco::PFCandidate* pfCand = &(*pf)[i];

        bool pf_isElectronForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::e;
        bool pf_isMuonForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;
        bool pf_isPhotonForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::gamma;
        bool pf_isChHadronForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h;
        bool pf_isNeutHadronForIdx = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h0;

        if (pfCand->trackRef().isNonnull() && pfCand->trackRef().key() == track.key()) {
          pf_isElectron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::e;
          pf_isMuon = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;
          pf_isPhoton = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::gamma;

          pf_isChHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h;
          pf_isNeutHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h0;
          pf_isUndefined = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::X;
          pf_isPfTrack = true;
          pf_energy = pfCand->ecalEnergy() + pfCand->hcalEnergy();
          pf_ecal_energy = pfCand->ecalEnergy();
          pf_hcal_energy = pfCand->hcalEnergy();

          if (trigInfo_ > 0) {
              // Number of PF tracks matched to HSCP candidate track
            tuple->BefPreS_PfType->Fill(1., EventWeight_);
            if (pf_isElectron) {
              tuple->BefPreS_PfType->Fill(2., EventWeight_);
            } else if (pf_isMuon) {
              tuple->BefPreS_PfType->Fill(3., EventWeight_);
            } else if (pf_isPhoton) {
              tuple->BefPreS_PfType->Fill(4., EventWeight_);
            } else if (pf_isChHadron) {
              tuple->BefPreS_PfType->Fill(5., EventWeight_);
            } else if (pf_isNeutHadron) {
              tuple->BefPreS_PfType->Fill(6., EventWeight_);
            } else if (pf_isUndefined) {
              tuple->BefPreS_PfType->Fill(7., EventWeight_);
            } else {
              tuple->BefPreS_PfType->Fill(8., EventWeight_);
            }
          }
          if (debug_ > 4 && trigInfo_ > 0) LogPrint(MOD) << "      >> HSCP candidate track has ID " << pfCand->pdgId() << " categoriezed by PF as " << pfCand->translatePdgIdToType(pfCand->pdgId());
            // The sum of the pt in the cone does not contain the pt of the track
            // just the pt of the surrounding tracks in the cone
          continue;
        }

        float dr = deltaR(pfCand->eta(),pfCand->phi(),track->eta(),track->phi());
        bool fromPV = (fabs(dz) < 0.1);

        float pt = pfCand->p4().pt();
        float drForMiniIso = 0.0;
        if (track->pt() < 50 ) {
          drForMiniIso = 0.2;
        } else if (track->pt() < 200) {
          drForMiniIso = 10/track->pt();
        } else {
          drForMiniIso = 0.05;
        }
        if (dr<drForMiniIso) {
            // Leptons get added to trackIso (this is not in the official definition)
          if (pf_isElectronForIdx || pf_isMuonForIdx) track_PFMiniIso_sumLeptonPt+=pt;
            // charged cands from PV get added to trackIso
          if(pf_isChHadronForIdx && fromPV) track_PFMiniIso_sumCharHadPt+=pt;
            // charged cands not from PV get added to pileup iso
          else if(pf_isChHadronForIdx) track_PFMiniIso_sumPUPt+=pt;
            // neutral hadron iso
          if(pf_isNeutHadronForIdx) track_PFMiniIso_sumNeutHadPt+=pt;
            // photon iso
          if(pf_isPhotonForIdx) track_PFMiniIso_sumPhotonPt+=pt;
            // muon iso
          if(pf_isMuonForIdx) track_PFMiniIso_sumMuonPt+=pt;
          if (!pf_isElectronForIdx && !pf_isMuonForIdx && !pf_isChHadronForIdx && !pf_isNeutHadronForIdx && !pf_isPhotonForIdx) {
            track_PFMiniIso_otherPt+=pt;
          }
        }
        if (dr<0.05) {
          // charged cands from PV get added to trackIso
          if(pf_isChHadronForIdx && fromPV) track_PFIso005_sumCharHadPt+=pt;
          // charged cands not from PV get added to pileup iso
          else if(pf_isChHadronForIdx) track_PFIso005_sumPUPt+=pt;
          // neutral hadron iso
          if(pf_isNeutHadronForIdx) track_PFIso005_sumNeutHadPt+=pt;
          // photon iso
          if(pf_isPhotonForIdx) track_PFIso005_sumPhotonPt+=pt;
        } if(dr<0.1){
          if(pf_isChHadronForIdx && fromPV) track_PFIso01_sumCharHadPt+=pt;
          else if(pf_isChHadronForIdx) track_PFIso01_sumPUPt+=pt;
          if(pf_isNeutHadronForIdx) track_PFIso01_sumNeutHadPt+=pt;
            if(pf_isPhotonForIdx) track_PFIso01_sumPhotonPt+=pt;
        } if(dr<0.3){
          if(pf_isChHadronForIdx && fromPV) track_PFIso03_sumCharHadPt+=pt;
          else if(pf_isChHadronForIdx) track_PFIso03_sumPUPt+=pt;
          else if(pf_isNeutHadronForIdx) track_PFIso03_sumNeutHadPt+=pt;
          else if(pf_isPhotonForIdx) track_PFIso03_sumPhotonPt+=pt;
        } if(dr<0.5){
          if(pf_isChHadronForIdx && fromPV) track_PFIso05_sumCharHadPt+=pt;
          else if(pf_isChHadronForIdx) track_PFIso05_sumPUPt+=pt;
          else if(pf_isNeutHadronForIdx) track_PFIso05_sumNeutHadPt+=pt;
          else if(pf_isPhotonForIdx) track_PFIso05_sumPhotonPt+=pt;
        }
      } // end loop PFCandidates
    } // close condition on the validity of PFCandidates

    // Calculate PF mini relative isolation
    // float miniRelIsoOfficial = (track_PFMiniIso_sumCharHadPt + std::max(0.0, track_PFMiniIso_sumNeutHadPt + track_PFMiniIso_sumPhotonPt - 0.5* track_PFMiniIso_sumPUPt))/track->pt();
    float miniRelIsoAll = (track_PFMiniIso_sumLeptonPt + track_PFMiniIso_otherPt + track_PFMiniIso_sumCharHadPt + std::max(0.0, track_PFMiniIso_sumNeutHadPt + track_PFMiniIso_sumPhotonPt - 0.5* track_PFMiniIso_sumPUPt))/track->pt();
    float miniRelIsoChg = track_PFMiniIso_sumCharHadPt/track->pt();
//  float miniRelIsoAll = (track_PFMiniIso_sumCharHadPt + std::max(0.0, track_PFMiniIso_sumNeutHadPt + track_PFMiniIso_sumPhotonPt - 0.5* track_PFMiniIso_sumPUPt))/track->pt();
    float miniRelIsoAll_wMuon = (track_PFMiniIso_sumMuonPt + track_PFMiniIso_sumCharHadPt + std::max(0.0, track_PFMiniIso_sumNeutHadPt + track_PFMiniIso_sumPhotonPt - 0.5* track_PFMiniIso_sumPUPt))/track->pt();

    HSCP_count++;

    std::vector<float> clust_charge;
    std::vector<float> clust_pathlength;
    std::vector<unsigned int> clust_nstrip;
    std::vector<bool> clust_sat254;
    std::vector<bool> clust_sat255;
    std::vector<uint32_t> clust_detid;
    std::vector<bool> clust_isStrip;
    std::vector<bool> clust_isPixel;

    // Associate gen track to reco track
    // If not associate gen track exists that means this is a fake track!
    // Let's skip the track in that case
    if (!isData && hasMCMatch_) {
      // Handle for the gen association of the track
      edm::Handle<edm::Association<reco::GenParticleCollection>>  trackToGenAssocHandle = iEvent.getHandle(trackToGenToken_);
      if (!trackToGenAssocHandle.isValid()) {
        // This became default from 12_0_X, in 10_6_X it's gated behind the bParking modifier
        LogPrint(MOD) << "trackToGenAssocHandle is invalid -- this should never happen in the latest AODSIM"
                      << "Please set hasMCMatch_ to false or move to a newer campaign";
        continue;
      }
      const auto& trackToGenAssoc = *trackToGenAssocHandle;
      reco::GenParticleRef genCollForTrack = trackToGenAssoc[track]; //.key()];
      if (genCollForTrack.isNull()) {
        LogPrint(MOD) << "  >> No associated gen track to this candidate -- this is a fake track, skipping it";
        continue;
      }
    }

    //load quantity associated to this track (TOF and dEdx)
    const reco::DeDxHitInfo* dedxHits = nullptr;
    if (typeMode_ != 3 && !track.isNull()) {
      reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
      if (!dedxHitsRef.isNull())
        dedxHits = &(*dedxHitsRef);
    }

    reco::DeDxHitInfoRef dedxHitsPrescRef;
    dedxHitsPrescRef = dedxCollH->get(track.key());
    int preScaleForDeDx = 1;
    if (!dedxHitsPrescRef.isNull()) {
      preScaleForDeDx = (*dedxPrescCollH)[dedxHitsPrescRef];
    }

    if (typeMode_ > 1 && typeMode_ != 5 && !hscp.muonRef().isNull()) {
      if (!isData) {
        tof = &(*tofMap)[hscp.muonRef()];
        dttof = &(*tofDtMap)[hscp.muonRef()];
        csctof = &(*tofCscMap)[hscp.muonRef()];
      } else {
        const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollH;
        const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollH;
        // Apply T0 correction on data but not on signal MC
        tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, 1);
        tof = &tofCalculator.combinedTOF;
        dttof = &tofCalculator.dtTOF;
        csctof = &tofCalculator.cscTOF;
      }
    } // end conditions for TOF including analysis variables

    // skip tracks without hits otherwise there will be a crash
    if (!dedxHits) {
      if (debug_> 3 && trigInfo_ > 0) LogPrint(MOD) << "No dedxHits associated to this track, skipping it";
      // 7-th bin of the error histo, No dedxHits associated to this track
      if (trigInfo_ > 0) tuple->ErrorHisto->Fill(6.);
      continue;
    }

    // Loop on the gen particles (again) to find if the 0.001 enviroment of the candidate has 91 or >2
    bool candidateEnvHasStatus91 = false;
    bool candidateEnvHasStatusHigherThan2 = false;
    if (!isData) {
      unsigned int usignedIntclosestGenIndex = 0;
      if (closestGenIndex>0) usignedIntclosestGenIndex = closestGenIndex;

      for (unsigned int g = 0; g < genColl.size(); g++) {
        // Exclude the canidate when looking at its envirment
        if (g == usignedIntclosestGenIndex) continue;
        // Look only at the R=0.001 enviroment of the candidate
        if (deltaR(genColl[g].eta(),genColl[g].phi(),genColl[usignedIntclosestGenIndex].eta(),genColl[usignedIntclosestGenIndex].phi()) > 0.001) continue;
        if (genColl[g].status() == 91) {
          candidateEnvHasStatus91 = true;
        }
        if (genColl[g].status() > 2) {
          candidateEnvHasStatusHigherThan2 = true;
        }
      }
    }

    int nofClust_dEdxLowerThan = 0;


    if (!isData && candidateEnvHasStatus91) {
      if (trigInfo_ > 0) {
        tuple->BefPreS_IasForStatus91->Fill(globalIas_, EventWeight_);
        tuple->ErrorHisto->Fill(8.);
      }
      ErrorHisto_bin = 8;
      continue;
    }

    passTechnicalChecks = true;
    
    // Include probQonTrack, probXYonTrack, probQonTrackNoL1, probXYonTrackNoL1 into one array
    float pixelProbs[4] = {0.0,0.0,0.0,0.0};
    int numRecHitsQ = 0, numRecHitsXY = 0;
    int numRecHitsQNoL1 = 0, numRecHitsXYNoL1 = 0;
    float probQonTrackWMulti = 1;
    float probXYonTrackWMulti = 1;
    float probQonTrackWMultiNoL1 = 1;
    float probXYonTrackWMultiNoL1 = 1;
    
    float ratioCleanAndAllStripsClu = 0.0;
    float ratioCleanAndAllPixelClu = 0.0;
    
    // Loop through the rechits on the given track **before** preselection
    unsigned int nonL1PixHits = 0;
    unsigned int allStripsHits = 0;
    unsigned int cleanStripsHits = 0;
    for (unsigned int i = 0; i < dedxHits->size(); i++) {
      clust_charge.push_back(dedxHits->charge(i));
      clust_pathlength.push_back(dedxHits->pathlength(i));
      clust_isStrip.push_back(dedxHits->detId(i) >= 3 ? true : false);
      clust_isPixel.push_back(dedxHits->detId(i) >= 3 ? false : true);
      clust_detid.push_back(dedxHits->detId(i));
      DetId detid(dedxHits->detId(i));
      float factorChargeToE = (detid.subdetId() < 3) ? 3.61e-06 : 3.61e-06 * 265;

      if (detid.subdetId() < 3) {
        // Calculate probQ and probXY for this pixel rechit
        // Taking the pixel cluster
        auto const* pixelCluster =  dedxHits->pixelCluster(i);
        if (pixelCluster == nullptr) {
          if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "    >> No dedxHits associated to this pixel cluster, skipping it";
          if (debug_> 0 && trigInfo_ > 0) LogPrint(MOD) << "    >> At this point this should never happen";
          continue;
        }
        // Check on which geometry unit the hit is
        const GeomDetUnit& geomDet = *tkGeometry->idToDetUnit(detid);
        // Get the local vector for the track direction
        LocalVector lv = geomDet.toLocal(GlobalVector(track->px(), track->py(), track->pz()));
        // Re-run the CPE on this cluster with the lv above
        // getParameters will return std::tuple<LocalPoint, LocalError, SiPixelRecHitQuality::QualWordType>;
        // from this we pick the 2nd, the QualWordType
        auto reCPE = std::get<2>(pixelCPE->getParameters(*pixelCluster, geomDet, LocalTrajectoryParameters(dedxHits->pos(i), lv, track->charge())));
        // extract probQ and probXY from this
        float probQ = SiPixelRecHitQuality::thePacking.probabilityQ(reCPE);
        float probXY = SiPixelRecHitQuality::thePacking.probabilityXY(reCPE);

        // To measure how often the CPE fails
        bool cpeHasFailed = false;
        if (!SiPixelRecHitQuality::thePacking.hasFilledProb(reCPE)) {
          cpeHasFailed = true;
          if (trigInfo_ > 0) {
            tuple->BefPreS_CluProbHasFilled->Fill(0., EventWeight_);
          }
        } else {
          if (trigInfo_ > 0) {
            tuple->BefPreS_CluProbHasFilled->Fill(1., EventWeight_);
          }
        }
        if (cpeHasFailed) continue;

        if (probQ <= 0.0 || probQ >= 1.f) probQ = 1.f;
        if (probXY <= 0.0 || probXY >= 1.f) probXY = 0.f;

        bool isOnEdge = SiPixelRecHitQuality::thePacking.isOnEdge(reCPE);
        bool hasBadPixels = SiPixelRecHitQuality::thePacking.hasBadPixels(reCPE);
        bool spansTwoROCs = SiPixelRecHitQuality::thePacking.spansTwoROCs(reCPE);
        bool specInCPE = (isOnEdge || hasBadPixels || spansTwoROCs) ? true : false;
        auto cotAlpha = lv.x()/lv.z();
        auto cotBeta = lv.y()/lv.z();
        auto clustSize = pixelCluster->size();
        auto clustSizeX = pixelCluster->sizeX();
        auto clustSizeY = pixelCluster->sizeY();
        auto clustCharge = pixelCluster->charge();

        auto pixelNormCharge = cm2umUnit * dedxHits->charge(i) / dedxHits->pathlength(i);

        if (clustCharge != dedxHits->charge(i)) {
          LogPrint(MOD) << "clustCharge != dedxHits->charge(i) -- this shouldnt happen";
        }
        if (trigInfo_ > 0) {
          if ( detid.subdetId() == PixelSubdetector::PixelBarrel) {
            auto pixLayerIndex = abs(int(tTopo->pxbLayer(detid)));
            tuple->BefPreS_CluProbQVsPixelLayer->Fill(probQ, pixLayerIndex, EventWeight_);
            tuple->BefPreS_CluProbXYVsPixelLayer->Fill(probXY, pixLayerIndex, EventWeight_);
            tuple->BefPreS_CluNormChargeVsPixelLayer->Fill(pixelNormCharge, pixLayerIndex, EventWeight_);
            tuple->BefPreS_CluSizeVsPixelLayer->Fill(clustSize-0.5, pixLayerIndex, EventWeight_);
            tuple->BefPreS_CluSizeXVsPixelLayer->Fill(clustSizeX-0.5, pixLayerIndex, EventWeight_);
            
            tuple->BefPreS_CluSizeYVsPixelLayer->Fill(clustSizeY-0.5, pixLayerIndex, EventWeight_);
            if (isOnEdge) {
              tuple->BefPreS_CluSpecInCPEVsPixelLayer->Fill(0.5, pixLayerIndex, EventWeight_);
            } else if (hasBadPixels) {
              tuple->BefPreS_CluSpecInCPEVsPixelLayer->Fill(1.5, pixLayerIndex, EventWeight_);
            } else if (spansTwoROCs) {
              tuple->BefPreS_CluSpecInCPEVsPixelLayer->Fill(2.5, pixLayerIndex, EventWeight_);
            }
            tuple->BefPreS_CluSpecInCPEVsPixelLayer->Fill(3.5, pixLayerIndex, EventWeight_);
            
            if (probXY < globalMinTrackProbXYCut_ && !specInCPE) {
              tuple->BefPreS_CluCotBetaVsPixelLayer_lowProbXY->Fill(cotBeta, pixLayerIndex, EventWeight_);
              tuple->BefPreS_CluCotAlphaVsPixelLayer_lowProbXY->Fill(cotAlpha, pixLayerIndex, EventWeight_);
            } else if (probXY > globalMinTrackProbXYCut_ && !specInCPE) {
              tuple->BefPreS_CluCotBetaVsPixelLayer->Fill(cotBeta, pixLayerIndex, EventWeight_);
              tuple->BefPreS_CluCotAlphaVsPixelLayer->Fill(cotAlpha, pixLayerIndex, EventWeight_);
            }
            if (!isData && genGammaBeta > 0.31623 && genGammaBeta < 0.6 ) {
              tuple->BefPreS_CluNormChargeVsPixelLayer_lowBetaGamma->Fill(pixelNormCharge, pixLayerIndex, EventWeight_);
            }
          } // end of IF on barrel pixel
        }

//      if (probQ > 0.f && probXY > 0.1) {
        if (!specInCPE && probQ < 0.8) {
          numRecHitsQ++;
          // Calculate alpha term needed for the combination
          probQonTrackWMulti *= probQ;
        }

        if (!specInCPE && probQ < 0.8 && probXY > 0.f) {
          numRecHitsXY++;
          // Calculate alpha term needed for the combination
          probXYonTrackWMulti *= probXY;
        }

        // Have a separate variable that excludes Layer 1
        // Layer 1 was very noisy in 2017/2018
        if (( detid.subdetId() == PixelSubdetector::PixelEndcap) || (detid.subdetId() == PixelSubdetector::PixelBarrel && tTopo->pxbLayer(detid) != 1)) {
          nonL1PixHits++;
          float probQNoL1 = SiPixelRecHitQuality::thePacking.probabilityQ(reCPE);
          float probXYNoL1 = SiPixelRecHitQuality::thePacking.probabilityXY(reCPE);

          if (probQNoL1 <= 0.0 || probQNoL1 >= 1.f) probQNoL1 = 1.f;
          if (probXYNoL1 <= 0.0 || probXYNoL1 >= 1.f) probXYNoL1 = 0.f;

          if (!specInCPE && probQ < 0.8) {
            numRecHitsQNoL1++;
            // Calculate alpha term needed for the combination
            probQonTrackWMultiNoL1 *= probQNoL1;
          }
          if (!specInCPE && probQ < 0.8 && probXYNoL1 > 0.f) {
            numRecHitsXYNoL1++;
            // Calculate alpha term needed for the combination
            probXYonTrackWMultiNoL1 *= probXYNoL1;
          }
        } // end if on the noL1 pixel side
      } // end if on the pixel side
      else if (detid.subdetId() >= 3) {
        // Taking the strips cluster
        auto const* stripsCluster = dedxHits->stripCluster(i);
        if (stripsCluster== nullptr) {
           if (debug_> 0) LogPrint(MOD) << "    >> No dedxHits associated to this strips cluster, skipping it";
           if (debug_> 0) LogPrint(MOD) << "    >> At this point this should never happen";
           continue;
        }
        std::vector<int> ampl = convert(stripsCluster->amplitudes());
        bool sat254 = false, sat255 = false;
        for (unsigned int s = 0; s < ampl.size(); s++) {
          if (ampl[s] >= 254)
            sat254 = true;
          if (ampl[s] == 255)
            sat255 = true;
        }
        ampl = CrossTalkInv(ampl, 0.10, 0.04, true);
        clust_nstrip.push_back(ampl.size());
        clust_sat254.push_back(sat254);
        clust_sat255.push_back(sat255);
        
        allStripsHits++;
        if ((clusterCleaning(ampl, 1))) {
          cleanStripsHits++;
        }
        
        if (trigInfo_ > 0) {
          float stripNormCharge = cm2umUnit * dedxHits->charge(i) * 265 / dedxHits->pathlength(i);
          unsigned int stripLayerIndex = 0;
          if (detid.subdetId() == StripSubdetector::TIB) stripLayerIndex = abs(int(tTopo->tibLayer(detid)));
          if (detid.subdetId() == StripSubdetector::TOB) stripLayerIndex = abs(int(tTopo->tobLayer(detid))) + 4;
          if (detid.subdetId() == StripSubdetector::TID) stripLayerIndex = abs(int(tTopo->tidWheel(detid))) + 10;
          if (detid.subdetId() == StripSubdetector::TEC) stripLayerIndex = abs(int(tTopo->tecWheel(detid))) + 13;
          if (!isData && genGammaBeta > 0.31623 && genGammaBeta < 0.6 ) {
            tuple->BefPreS_CluNormChargeVsStripLayer_lowBetaGamma->Fill(stripNormCharge, stripLayerIndex, EventWeight_);
          } else if (!isData && genGammaBeta > 0.6 ) {
            tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma->Fill(stripNormCharge, stripLayerIndex, EventWeight_);
            if (!candidateEnvHasStatus91) {
              tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatNot91->Fill(stripNormCharge, stripLayerIndex, EventWeight_);
            }
            if (candidateEnvHasStatusHigherThan2) {
              tuple->BefPreS_CluNormChargeVsStripLayer_higherBetaGamma_StatHigherThan2->Fill(stripNormCharge, stripLayerIndex, EventWeight_);
            }
          }
        }

        if (dedxHits->charge(i) * factorChargeToE / dedxHits->pathlength(i) < theFMIPX_)
          nofClust_dEdxLowerThan++;
        }
      }// end loop on rechits on the given track

    // Combine probQ-s into HSCP candidate (track) level quantity
    float probQonTrack = -1.f;
    float probXYonTrack = -1.f;
    float probQonTrackNoL1 = -1.f;
    float probXYonTrackNoL1 = -1.f;
    probQonTrack = combineProbs(probQonTrackWMulti, numRecHitsQ);
    probXYonTrack = combineProbs(probXYonTrackWMulti, numRecHitsXY);
    probQonTrackNoL1 = combineProbs(probQonTrackWMultiNoL1, numRecHitsQNoL1);
    probXYonTrackNoL1 = combineProbs(probXYonTrackWMultiNoL1, numRecHitsXYNoL1);
    
    // Ratios of cleaned hits and all hits
    if (allStripsHits > 0) {
      ratioCleanAndAllStripsClu = cleanStripsHits / float(allStripsHits);
    }
    if (nonL1PixHits > 0) {
      ratioCleanAndAllPixelClu = numRecHitsQNoL1 / float(nonL1PixHits);
    }
    
    if (probQonTrackNoL1 < 0.0) cout << "probQonTrackNoL1 < 0.0, but it shouldnt...";
  
    // Cleaning of tracks that had failed the template CPE (prob <= 0.0 and prob >= 1.0 cases)
    if (probQonTrack < 0.0 || probXYonTrack < 0.0 || probQonTrack > 1.f || probXYonTrack > 1.f) {
      if (debug_> 2) LogPrint(MOD) << "    >> Probs out of bound: " <<
        " ProbQ = " << probQonTrack << " ProbXY = " << probXYonTrack <<  " ProbQNoL1 = "<< probQonTrackNoL1 << " ProbXYNoL1 = " << probXYonTrackNoL1;
    }
    
    if (trigInfo_ > 0) {
      tuple->BefPreS_genGammaBetaVsProbXYNoL1->Fill(genGammaBeta, probXYonTrackNoL1, EventWeight_);
    }
  
    float Fmip = (float)nofClust_dEdxLowerThan / (float)dedxHits->size();

    //computedEdx: hits, SF, templates, usePixel, useStrips,, useClusterCleaning, uneTrunc,
    //             mustBeInside, MaxStripNOM, correctFEDSat, XtalkInv, lowDeDxDrop, dedxErr, useTemplateLayer_, skip_templates_ias
    //
    //correction inverseXtalk = 0 --> take the raw amplitudes of the cluster
    //correction inverseXtalk = 1 --> modify the amplitudes based on xtalk for non-saturated cluster + correct for saturation
    //
    //skip_templates_ias = 0 --> no skip
    //skip_templates_ias = 1 --> no Pix, no TIB, no TID, no 3 first layers TEC
    //skip_templates_ias = 2 --> Pixel Only


    float dEdxErr = 0;
    bool symmetricSmirnov = false;

    // Ih
    auto dedxMObj_FullTrackerTmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = true,  useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, &dEdxErr, useTemplateLayer_);

    reco::DeDxData* dedxMObj_FullTracker = dedxMObj_FullTrackerTmp.numberOfMeasurements() > 0 ? &dedxMObj_FullTrackerTmp : nullptr;

    // Ih Up
    auto dedxMUpObjTmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.15, 0, useTemplateLayer_);

    reco::DeDxData* dedxMUpObj = dedxMUpObjTmp.numberOfMeasurements() > 0 ? &dedxMUpObjTmp : nullptr;

    // Ih Down
    // For now it's a copy of Ih Up, I doubt that's what it should be...
    // Also I think this should be done on the top of Ih no pixel L1 not the full tracker version
    auto dedxMDownObjTmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.15, 0, useTemplateLayer_);

    reco::DeDxData* dedxMDownObj = dedxMDownObjTmp.numberOfMeasurements() > 0 ? &dedxMDownObjTmp : nullptr;

    // Ih no pixel L1
    auto dedxIh_noL1_Tmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, &dEdxErr, useTemplateLayer_, skipPixelL1 = true);

    reco::DeDxData* dedxIh_noL1 = dedxIh_noL1_Tmp.numberOfMeasurements() > 0 ? &dedxIh_noL1_Tmp : nullptr;

    // Ih 0.15 low values drop
    // Should useTruncated be true ?
    auto dedxIh_15drop_Tmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = true,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.15, &dEdxErr, useTemplateLayer_);

    reco::DeDxData* dedxIh_15drop = dedxIh_15drop_Tmp.numberOfMeasurements() > 0 ? &dedxIh_15drop_Tmp : nullptr;

    // Ih Strip only
    auto dedxIh_StripOnly_Tmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = false, useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, &dEdxErr, useTemplateLayer_);

    reco::DeDxData* dedxIh_StripOnly = dedxIh_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIh_StripOnly_Tmp : nullptr;

    // Ih Strip only and 0.15 low values drop
    auto dedxIh_StripOnly_15drop_Tmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = false, useStrip = true, useClusterCleaning, useTruncated = true,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.15, &dEdxErr, useTemplateLayer_, skipPixelL1 = true);

    reco::DeDxData* dedxIh_StripOnly_15drop = dedxIh_StripOnly_15drop_Tmp.numberOfMeasurements() > 0 ? &dedxIh_StripOnly_15drop_Tmp : nullptr;

    // Ih Pixel only no BPIXL1
    auto dedxIh_PixelOnly_noL1_Tmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = true, useStrip = false, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, &dEdxErr, useTemplateLayer_, skipPixelL1 = true);

    reco::DeDxData* dedxIh_PixelOnlyh_noL1 = dedxIh_PixelOnly_noL1_Tmp.numberOfMeasurements() > 0 ? &dedxIh_PixelOnly_noL1_Tmp : nullptr;

    // Ih correct saturation from fits
    auto dedxIh_SaturationCorrectionFromFits_Tmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = false, useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 2, dropLowerDeDxValue = 0.0, &dEdxErr, useTemplateLayer_, skipPixelL1 = true);

    reco::DeDxData* dedxIh_SaturationCorrectionFromFits = dedxIh_SaturationCorrectionFromFits_Tmp.numberOfMeasurements() > 0 ? &dedxIh_SaturationCorrectionFromFits_Tmp : nullptr;



    // Temporary dEdx info
    reco::DeDxData dedxIas_FullTrackerTmp,dedxIas_noL1Tmp,dedxIas_noTIBnoTIDno3TEC_Tmp,dedxIas_PixelOnly_Tmp,dedxIas_StripOnly_Tmp,dedxIas_PixelOnly_noL1_Tmp,dedxIs_StripOnly_Tmp;

    //Pointers that will have the dEdx info laterr
    reco::DeDxData* dedxIas_FullTracker = nullptr; //globalIas_
    reco::DeDxData* dedxIas_noL1 = nullptr; //globalIas_ no BPIXL1
    reco::DeDxData* dedxIas_noTIBnoTIDno3TEC = nullptr; //globalIas_ without TIB, TID, and 3 first TEC layers

    reco::DeDxData* dedxIas_PixelOnly = nullptr; //globalIas_ Pixel only
    reco::DeDxData* dedxIas_StripOnly = nullptr; //globalIas_ Strip only
    reco::DeDxData* dedxIas_PixelOnly_noL1 = nullptr; //globalIas_ Pixel only no BPIXL1
    reco::DeDxData* dedxIs_StripOnly = nullptr; //symmetric Smirnov discriminator - Is

    int NPV = vertexColl.size();
    //if(!createGiTemplates_){
        if(!puTreatment_ || CreateAndRunGitemplates_){
            dedxIas_FullTrackerTmp =
            computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplates, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                       mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_);

            dedxIas_FullTracker = dedxIas_FullTrackerTmp.numberOfMeasurements() > 0 ? &dedxIas_FullTrackerTmp : nullptr;

            dedxIas_noL1Tmp =
            computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplates, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                   mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_,skipPixelL1 = true, skip_templates_ias = 2);

            dedxIas_noL1 = dedxIas_noL1Tmp.numberOfMeasurements() > 0 ? &dedxIas_noL1Tmp : nullptr;

            dedxIas_noTIBnoTIDno3TEC_Tmp =
            computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplates, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                        mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = true, skip_templates_ias = 1);


            dedxIas_noTIBnoTIDno3TEC = dedxIas_noTIBnoTIDno3TEC_Tmp.numberOfMeasurements() > 0 ? &dedxIas_noTIBnoTIDno3TEC_Tmp : nullptr;

            dedxIas_PixelOnly_Tmp =
            computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplates, usePixel = true, useStrip = false, useClusterCleaning, useTruncated = false,
                        mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = false, skip_templates_ias = 2);

            dedxIas_PixelOnly = dedxIas_PixelOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIas_PixelOnly_Tmp : nullptr;

            dedxIas_StripOnly_Tmp =
            computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplates, usePixel = false, useStrip = true, useClusterCleaning, useTruncated = false,
                        mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = false, skip_templates_ias = 0);

            dedxIas_StripOnly = dedxIas_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIas_StripOnly_Tmp : nullptr;

            dedxIas_PixelOnly_noL1_Tmp =
            computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplates, usePixel = true, useStrip = false, useClusterCleaning, useTruncated = false,
                        mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = true, skip_templates_ias = 2);

            dedxIas_PixelOnly_noL1 = dedxIas_PixelOnly_noL1_Tmp.numberOfMeasurements() > 0 ? &dedxIas_PixelOnly_noL1_Tmp : nullptr;

            dedxIs_StripOnly_Tmp =
        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplates, usePixel = true, useStrip = false, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = true, skip_templates_ias = 2, symmetricSmirnov = true);

            dedxIs_StripOnly = dedxIs_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIs_StripOnly_Tmp : nullptr;
        }

        else{
            for(int i = 0 ; i < NbPuBins_ ; i++){
                if ( NPV > PuBins_[i] && NPV <= PuBins_[i+1] ){
                    //globalIas_
                    dedxIas_FullTrackerTmp =
                    computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplatesPU[i], usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                                mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_);

                    dedxIas_FullTracker = dedxIas_FullTrackerTmp.numberOfMeasurements() > 0 ? &dedxIas_FullTrackerTmp : nullptr;

                    //globalIas_ no BPIXL1
                    dedxIas_noL1Tmp =
                    computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplatesPU[i], usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                           mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_,skipPixelL1 = true, skip_templates_ias = 2);

                    dedxIas_noL1 = dedxIas_noL1Tmp.numberOfMeasurements() > 0 ? &dedxIas_noL1Tmp : nullptr;

                    //globalIas_ without TIB, TID, and 3 first TEC layers
                    dedxIas_noTIBnoTIDno3TEC_Tmp =
                        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplatesPU[i], usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = true, skip_templates_ias = 1);

                    dedxIas_noTIBnoTIDno3TEC = dedxIas_noTIBnoTIDno3TEC_Tmp.numberOfMeasurements() > 0 ? &dedxIas_noTIBnoTIDno3TEC_Tmp : nullptr;

                    //globalIas_ Pixel only
                    dedxIas_PixelOnly_Tmp =
                        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplatesPU[i], usePixel = true, useStrip = false, useClusterCleaning, useTruncated = false,
                                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = false, skip_templates_ias = 2);

                    dedxIas_PixelOnly = dedxIas_PixelOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIas_PixelOnly_Tmp : nullptr;

                    //globalIas_ Strip only
                    dedxIas_StripOnly_Tmp =
                        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplatesPU[i], usePixel = false, useStrip = true, useClusterCleaning, useTruncated = false,
                                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = false, skip_templates_ias = 0);

                    dedxIas_StripOnly = dedxIas_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIas_StripOnly_Tmp : nullptr;

                    //globalIas_ Pixel only no BPIXL1
                    dedxIas_PixelOnly_noL1_Tmp =
                        computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplatesPU[i], usePixel = true, useStrip = false, useClusterCleaning, useTruncated = false,
                                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = true, skip_templates_ias = 2);

                    dedxIas_PixelOnly_noL1 = dedxIas_PixelOnly_noL1_Tmp.numberOfMeasurements() > 0 ? &dedxIas_PixelOnly_noL1_Tmp : nullptr;

                    //symmetric Smirnov discriminator - Is
                    dedxIs_StripOnly_Tmp =
                    computedEdx(track->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = dEdxTemplatesPU[i], usePixel = true, useStrip = false, useClusterCleaning, useTruncated = false,
                                mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, 0, useTemplateLayer_, skipPixelL1 = true, skip_templates_ias = 2, symmetricSmirnov = true);

                    dedxIs_StripOnly = dedxIs_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIs_StripOnly_Tmp : nullptr;

                }//end vertex coll size for PU treatment
            }
        }
    //}
    //Choose of Ih definition - Ih_nodrop_noPixL1
    auto dedxMObj = dedxIh_noL1;
    globalIh_ = (dedxMObj) ?  dedxMObj->dEdx() : -1.f;
    
    //Choose of Ias definition - strips only
    auto dedxSObj = dedxIas_StripOnly;
    globalIas_ = (dedxSObj) ? dedxSObj->dEdx() : -1.f;

    float MassErr = GetMassErr(track->p(),
                                track->ptError(),
                                dedxMObj ? dedxMObj->dEdx() : -1,
                                dEdxErr,
                                GetMass(track->p(), dedxMObj ? dedxMObj->dEdx() : -1, dEdxK_, dEdxC_),
                                dEdxK_,
                                dEdxC_);

    // ------------------------------------------------------------------------------------
    //compute systematic uncertainties on signal
    if (isSignal && calcSyst_) {
      if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> Compute systematic uncertainties on signal";
      calculateSyst(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, iSetup, pixelProbs, tuple, -1, MassErr, closestBackgroundPDGsIDs);
    }  //End of systematic computation for signal
    // ------------------------------------------------------------------------------------
    // Fill up the closestBackgroundPDGsIDs array (has to be done before preselection function)
    if (!isData && closestGenIndex > 0) {
        //      if (debug_> 0) LogPrint(MOD) << "  >> Background MC, set gen IDs, mother IDs, sibling IDs";
      closestBackgroundPDGsIDs[0] = (float)abs(genPdgId);
      float genEta = genColl[closestGenIndex].eta();
      float genPhi = genColl[closestGenIndex].phi();
      float dRMinGenAndSibling = 9999.0;
      float dRMinGenAndMom = 9999.0;
      float numSiblingsF = 9999.0;
      bool motherFound = false;
      float dRMinGenAndAunt = 9999.0;
      reco::GenParticle& genCandidateUnderStudy = genColl[closestGenIndex];

        // Loop through all the mothers of the gen particle
      for (unsigned int numMomIndx = 0; numMomIndx < genCandidateUnderStudy.numberOfMothers(); numMomIndx++) {
        if (abs(genCandidateUnderStudy.mother(numMomIndx)->pdgId())  != abs(genCandidateUnderStudy.pdgId())) {
          closestBackgroundPDGsIDs[1] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->pdgId());
          closestBackgroundPDGsIDs[7] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->pt());
          unsigned int numSiblings = genCandidateUnderStudy.mother(numMomIndx)->numberOfDaughters() -1;
          numSiblingsF  = float(numSiblings);
          for (unsigned int daughterIndx = 0; daughterIndx < numSiblings+1; daughterIndx++) {
            if (globalIas_ > 0.3 && debug_ > 4 && trigInfo_ > 0)  std::cout << "      >> " << genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->pdgId() ;
            float siblingEta = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->eta();
            float siblingPhi = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->phi();
            float siblingDr = deltaR(genEta, genPhi, siblingEta, siblingPhi);
            if( (siblingDr != 0.0) && (siblingDr < dRMinGenAndSibling)) {
              dRMinGenAndSibling = siblingDr;
              closestBackgroundPDGsIDs[2] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->pdgId());
            }
          }
          float momEta = genCandidateUnderStudy.mother(numMomIndx)->eta();
          float momPhi = genCandidateUnderStudy.mother(numMomIndx)->phi();
          dRMinGenAndMom = deltaR(genEta, genPhi, momEta, momPhi);
          motherFound = true;
          break;
        }
      }
      
      // If the loop on the mothers didnt find the mother (e.g. all moms had the same ID), let's look at the grandmas
      if (!motherFound) {

        for (unsigned int numMomIndx = 0; numMomIndx < genCandidateUnderStudy.numberOfMothers(); numMomIndx++) {
          for (unsigned int numGramMomIndx = 0; numGramMomIndx < genCandidateUnderStudy.mother(numMomIndx)->numberOfMothers(); numGramMomIndx++) {
            if (abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->pdgId())  != abs(genCandidateUnderStudy.pdgId())) {
              closestBackgroundPDGsIDs[1] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->pdgId());
              closestBackgroundPDGsIDs[7] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->pt());
              unsigned int numSiblings = genCandidateUnderStudy.mother(numMomIndx)->numberOfDaughters() -1;

              unsigned int numAunts = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->numberOfDaughters() -1;
              numSiblingsF  = float(numSiblings);
              for (unsigned int daughterIndx = 0; daughterIndx < numAunts+1; daughterIndx++) {
                float auntEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->eta();
                float auntPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->phi();
                float auntDr = deltaR(genEta, genPhi, auntEta, auntPhi);
                if( (auntDr != 0.0) && (auntDr < dRMinGenAndAunt)) {
                  dRMinGenAndAunt = auntDr;
                  closestBackgroundPDGsIDs[2] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->pdgId());
                }
              }
              float momEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->eta();
              float momPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->phi();
              dRMinGenAndMom = deltaR(genEta, genPhi, momEta, momPhi);
              motherFound = true;
              break;
            }
          }
          if (motherFound) break;
        }
      }

        // If none of the mothers' mother's is the real mother (e.g. all moms'moms had the same ID as the candidate), let's look at the grand-grandmas
      if (!motherFound) {
        for (unsigned int numMomIndx = 0; numMomIndx < genCandidateUnderStudy.numberOfMothers(); numMomIndx++) {
          for (unsigned int numGramMomIndx = 0; numGramMomIndx < genCandidateUnderStudy.mother(numMomIndx)->numberOfMothers(); numGramMomIndx++) {
            for (unsigned int numGrandGramMomIndx = 0; numGrandGramMomIndx < genCandidateUnderStudy.mother(numGramMomIndx)->mother(numGramMomIndx)->numberOfMothers(); numGrandGramMomIndx++) {
              if (abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->pdgId())  != abs(genCandidateUnderStudy.pdgId())) {
                closestBackgroundPDGsIDs[1] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->pdgId());
                closestBackgroundPDGsIDs[7] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->pt());
                unsigned int numSiblings = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->numberOfDaughters() -1;
                numSiblingsF  = float(numSiblings);
                for (unsigned int daughterIndx = 0; daughterIndx < numSiblings+1; daughterIndx++) {
                  if (globalIas_ > 0.3 && debug_ > 4 && trigInfo_ > 0) std::cout << "      >> "  << genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->pdgId() ;
                  float siblingEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->eta();
                  float siblingPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->phi();
                  float siblingDr = deltaR(genEta, genPhi, siblingEta, siblingPhi);
                  if( (siblingDr != 0.0) && (siblingDr < dRMinGenAndSibling)) {
                    dRMinGenAndSibling = siblingDr;
                    closestBackgroundPDGsIDs[2] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->pdgId());
                  }
                }
                float momEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->eta();
                float momPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->phi();
                dRMinGenAndMom = deltaR(genEta, genPhi, momEta, momPhi);
                motherFound = true;
                break;
              }
            }
            if (motherFound) break;
          }
          if (motherFound) break;
        }
      }

      // I'm sure this could be done better, if you agree and feel like it, please fix it
      // issue with a while loop and a recursive I faced is tha that mom doesnt have the same type as the genParticle
      // This is also a dupplicate code, prob should be a function

      closestBackgroundPDGsIDs[3] = dRMinGenAndSibling;
      closestBackgroundPDGsIDs[4] = dRMinGenAndMom;
      closestBackgroundPDGsIDs[5] = fabs(genPt);
      closestBackgroundPDGsIDs[6] = numSiblingsF;
    }
    // -- end TODO Sept 25

      //  // Before preselection print-outs
      //  if (debug_ > 7 ) {
      //    LogPrint(MOD) << "        >> Before preselection print-outs:";
      //    LogPrint(MOD) << "        >> Trigger passed!";
      //    LogPrint(MOD) << "        >>   track->eta()  " <<   track->eta() ;
      //    LogPrint(MOD) << "        >>   track->pt()  " <<   track->pt() ;
      //    LogPrint(MOD) << "        >>   track->found()  " <<   track->found() ;
      //    LogPrint(MOD) << "        >>   track->hitPattern().numberOfValidPixelHits()  " <<   track->hitPattern().numberOfValidPixelHits() ;
      //    LogPrint(MOD) << "        >>   track->validFraction()  " <<   track->validFraction() ;
      //    LogPrint(MOD) << "        >>   numDeDxHits  " <<   numDeDxHits ;
      //    LogPrint(MOD) << "        >>   track->chi2() / track->ndof()   " <<   track->chi2() / track->ndof()  ;
      //    LogPrint(MOD) << "        >>   EoP   " <<   EoP  ;
      //    LogPrint(MOD) << "        >>   PF E = " << pf_energy <<  " Cone based (0.3) E = " << hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy() << " p = " << track->p() ;
      //    LogPrint(MOD) << "        >>   dz  " <<   dz ;
      //    LogPrint(MOD) << "        >>   dxy  " <<   dxy ;
      //    LogPrint(MOD) << "        >>   track->ptError() / track->pt()  " <<   track->ptError() / track->pt() ;
      //    LogPrint(MOD) << "        >>   pTerr_over_pT_etaBin(track->pt(), track->eta())  " <<   pTerr_over_pT_etaBin(track->pt(), track->eta()) ;
      //    LogPrint(MOD) << "        >>   IsoTK_SumEt   " <<   IsoTK_SumEt  ;
      //    LogPrint(MOD) << "        >>   miniRelIsoAll   " <<   miniRelIsoAll  ;
      //    LogPrint(MOD) << "        >>   globalIh_  " <<   globalIh_ ;
      //    LogPrint(MOD) << "        >>   globalIas_  " << globalIas_;
      //    LogPrint(MOD) << "        >>   probQonTrack   " <<   probQonTrack  ;
      //    LogPrint(MOD) << "        >>   probXYonTrack   " <<  probXYonTrack  ;
      //
      //  }

      // Loop on generalTracks
    float track_genTrackMiniIsoSumPt = 0;
    for(unsigned int c=0;c<trackCollectionHandle->size();c++){
      reco::TrackRef genTrackRef = reco::TrackRef( trackCollectionHandle.product(), c );
        // Dont count the HSCP candidate in
      if (genTrackRef.isNonnull() && genTrackRef.key() != track.key()) {
        float drForMiniIso = 0.0;
        if (track->pt() < 50 ) {
          drForMiniIso = 0.2;
        } else if (track->pt() < 200) {
          drForMiniIso = 10/track->pt();
        } else {
          drForMiniIso = 0.05;
        }
        float pt = genTrackRef->pt();
        float dr = deltaR(genTrackRef->eta(),genTrackRef->phi(),track->eta(),track->phi());
        if (dr<drForMiniIso) {
          track_genTrackMiniIsoSumPt+=pt;
        }
      }
    }

    // number of tracks as the first bin
    if (trigInfo_ > 0) {
      tuple->BefPreS_PfType->Fill(0., EventWeight_);
    }

    int nearestJetIndex = -1;
    if (pfJetHandle.isValid() && !pfJetHandle->empty()) {
      const reco::PFJetCollection* pfJetColl = pfJetHandle.product();
      for (unsigned int i = 0; i < pfJetColl->size(); i++) {
        const reco::PFJet* jet = &(*pfJetColl)[i];
        float dr = deltaR(jet->eta(), jet->phi(), track->eta(), track->phi());

        float dRMinPfCaloJet = 9999.0;
        float dPtPfCaloJet = 9999.0;
        if (caloJetHandle.isValid() && !caloJetHandle->empty()) {
          for (unsigned int iCalo = 0; iCalo < caloJetHandle->size(); iCalo++) {
            const reco::CaloJet* jetCalo = &(*caloJetHandle)[iCalo];
            float drCalo = deltaR(jet->eta(), jet->phi(), jetCalo->eta(), jetCalo->phi());

            if (drCalo < dRMinPfCaloJet) {
              dRMinPfCaloJet = drCalo;
              dPtPfCaloJet = fabs(jet->pt()-jetCalo->pt());
            }
          }
        }
        if (trigInfo_ > 0) {
          tuple->BefPreS_dRVsdPtPfCaloJet->Fill(dRMinPfCaloJet,dPtPfCaloJet, EventWeight_);
        }
        
        if (dr < dRMinPfJetWithNoCuts) {
          dRMinPfJetWithNoCuts = dr;
          nearestJetIndex = i;
        }

        if (jet->pt() > 50) {
          if (dr < dRMinPfJetWithCuts) {
            dRMinPfJetWithCuts = dr;
          }
        }
      }

      const reco::PFJet* nearestJetWithNoCuts = &(*pfJetColl)[nearestJetIndex];
      if (trigInfo_ > 0) {
        tuple->BefPreS_dRVsPtPfJet->Fill(dRMinPfJetWithNoCuts, nearestJetWithNoCuts->pt(), EventWeight_);
      }

      closestPfJetMuonFraction = nearestJetWithNoCuts->muonEnergyFraction();
      closestPfJetElectronFraction = nearestJetWithNoCuts->electronEnergyFraction();
      closestPfJetPhotonFraction = nearestJetWithNoCuts->photonEnergyFraction();
    }

      // loop on Calo jets
    float dRMinCaloJet = 9999.0;
    int caloJetsNum = 0;
    if (caloJetHandle.isValid() && !caloJetHandle->empty()) {
      for (unsigned int i = 0; i < caloJetHandle->size(); i++) {
        const reco::CaloJet* jet = &(*caloJetHandle)[i];
          if (jet->pt() < 50) {
          continue;
        }
        caloJetsNum++;
        float dr = deltaR(jet->eta(), jet->phi(), track->eta(), track->phi());
        if (dr < dRMinCaloJet) {
          dRMinCaloJet = dr;
        }
      }
    }
    
    // Comput the phi angle between MET and the candidate
    float dPhiMinPfMet = 9999.0;
    dPhiMinPfMet = fabs(reco::deltaPhi(RecoPFMET_phi,track->phi()));

      // Number of DeDx hits
    unsigned int numDeDxHitsOnStrips = dedxSObj ? dedxSObj->numberOfMeasurements() : 0;
    unsigned int numDeDxHits = numDeDxHitsOnStrips+nonL1PixHits;
    unsigned int missingHitsTillLast =
      track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
      track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
    float validFractionTillLast = track->found() <= 0 ? -1 : track->found() / float(track->found() + missingHitsTillLast);
    
    // A,B,C for 3 categories of PU
    bool PUA = (numGoodVerts < 15);
    bool PUB = (numGoodVerts >= 15 && numGoodVerts < 30);
    bool PUC = (numGoodVerts >= 30 );

    const edm::ValueMap<susybsm::HSCPIsolation> IsolationMap = iEvent.get(hscpIsoToken_);
    susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());

    //  float EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy()) / track->p();
    float EoP = pf_energy / track->p();
    float IsoTK_SumEt = hscpIso.Get_TK_SumEt();

    float Mass =  dedxMObj ?  GetMass(track->p(), globalIh_, dEdxK_, dEdxC_) : -1;

    // Find distance to nearest segment on opposite side of detector
    float minPhi = 0.0, minEta = 0.0;
    // This will return stuff if typeMode_ = 3
    float segSep = SegSep(track, iEvent, minPhi, minEta);


    // maximal deltaR angle between the candidate and any other canidate
    float maxOpenAngle = -1;
    maxOpenAngle = deltaROpositeTrack(iEvent.get(hscpToken_), hscp);

      // Get the location of the outmost hit
    const GlobalPoint outerHit = getOuterHitPos(iSetup, dedxHits);
    const float furthersHitDxy = sqrt(outerHit.x()*outerHit.x()+outerHit.y()*outerHit.y());
    const float furthersHitDistance = sqrt(outerHit.x()*outerHit.x()+outerHit.y()*outerHit.y()+outerHit.z()*outerHit.z());
    
    bool doBefTrigPlots = true;
    if (doBefTrigPlots) {
      tuple->BefTrig_Ih->Fill(globalIh_, EventWeight_);
      tuple->BefTrig_ProbQNoL1->Fill(1 - probQonTrackNoL1, EventWeight_);
      tuple->BefTrig_Ias->Fill(globalIas_, EventWeight_);
    }

    bool doBefPreSplots = true;
    // Before pre-selection, after trigger plots
    if (doBefPreSplots && trigInfo_ > 0) {
      tuple->BefPreS_RecoPFMET->Fill(RecoPFMET, EventWeight_);
      tuple->BefPreS_Eta->Fill(track->eta(), EventWeight_);
      tuple->BefPreS_MatchedStations->Fill(muonStations(track->hitPattern()), EventWeight_);
      tuple->BefPreS_NVertex->Fill(vertexColl.size(), EventWeight_);
      tuple->BefPreS_NVertex_NoEventWeight->Fill(vertexColl.size());
      if (PUA) {
        tuple->BefPreS_TNOH_PUA->Fill(track->found(), EventWeight_);
        tuple->BefPreS_TNOM_PUA->Fill(numDeDxHits, EventWeight_);
        tuple->BefPreS_Ias_PUA->Fill(globalIas_, EventWeight_);
        tuple->BefPreS_Ih_PUA->Fill(globalIh_, EventWeight_);
        tuple->BefPreS_Pt_PUA->Fill(track->pt(), EventWeight_);
      }
      if (PUB) {
        tuple->BefPreS_TNOH_PUB->Fill(track->found(), EventWeight_);
        tuple->BefPreS_TNOM_PUB->Fill(numDeDxHits, EventWeight_);
        tuple->BefPreS_Ias_PUB->Fill(globalIas_, EventWeight_);
        tuple->BefPreS_Ih_PUB->Fill(globalIh_, EventWeight_);
        tuple->BefPreS_Pt_PUB->Fill(track->pt(), EventWeight_);
      }
      tuple->BefPreS_TNOHFraction->Fill(track->validFraction(), EventWeight_);
      tuple->BefPreS_TNOPH->Fill(nonL1PixHits, EventWeight_);
      tuple->BefPreS_RatioCleanAndAllStripsClu->Fill(ratioCleanAndAllStripsClu, EventWeight_);
      tuple->BefPreS_RatioCleanAndAllPixelClu->Fill(ratioCleanAndAllPixelClu, EventWeight_);
      tuple->BefPreS_TNOHFractionTillLast->Fill(validFractionTillLast, EventWeight_);
      tuple->BefPreS_TNOMHTillLast->Fill(missingHitsTillLast, EventWeight_);
      tuple->BefPreS_TNOM->Fill(numDeDxHits, EventWeight_);
      if (track->found() - numDeDxHits) {
        tuple->BefPreS_EtaVsNBH->Fill(track->eta(), track->found() - numDeDxHits, EventWeight_);
      }
      tuple->BefPreS_ProbQ->Fill(1 - probQonTrack, EventWeight_);
      tuple->BefPreS_ProbXY->Fill(probXYonTrack, EventWeight_);
      tuple->BefPreS_ProbQNoL1->Fill(1 - probQonTrackNoL1, EventWeight_);
      tuple->BefPreS_ProbXYNoL1->Fill(probXYonTrackNoL1, EventWeight_);
      if (tof) {
        tuple->BefPreS_nDof->Fill(tof->nDof(), EventWeight_);
        tuple->BefPreS_MTOF->Fill(tof->inverseBeta(), EventWeight_);
        tuple->BefPreS_TOFError->Fill(tof->inverseBetaErr(), EventWeight_);
        tuple->BefPreS_TimeAtIP->Fill(tof->timeAtIpInOut(), EventWeight_);
      }
      if (track->quality(reco::TrackBase::highPurity)) {
        tuple->BefPreS_Qual->Fill(1., EventWeight_);
      } else {
        tuple->BefPreS_Qual->Fill(0., EventWeight_);
      }

      tuple->BefPreS_Chi2oNdof->Fill(track->chi2() / track->ndof(), EventWeight_);
      tuple->BefPreS_Pt->Fill(track->pt(), EventWeight_);
      tuple->BefPreS_Pt_lowPt->Fill(track->pt(), EventWeight_);
      tuple->BefPreS_P->Fill(track->p(), EventWeight_);
      tuple->BefPreS_NOMoNOH->Fill(numDeDxHits / (float)track->found(), EventWeight_);
      tuple->BefPreS_NOMoNOHvsPV->Fill(numGoodVerts, numDeDxHits / (float)track->found(), EventWeight_);
      tuple->BefPreS_Dxy->Fill(dxy, EventWeight_);
      tuple->BefPreS_Dz->Fill(dz, EventWeight_);
      tuple->BefPreS_EtaVsDz->Fill(track->eta(), dz, EventWeight_);
      tuple->BefPreS_PV->Fill(numGoodVerts, EventWeight_);
      tuple->BefPreS_PV_NoEventWeight->Fill(numGoodVerts);
      tuple->BefPreS_EoP->Fill(EoP, EventWeight_);
      tuple->BefPreS_SumpTOverpT->Fill(IsoTK_SumEt / track->pt(), EventWeight_);
      tuple->BefPreS_PtErrOverPt->Fill(track->ptError() / track->pt(), EventWeight_);
      tuple->BefPreS_PtErrOverPt2->Fill(track->ptError() / (track->pt()*track->pt()), EventWeight_);
      tuple->BefPreS_PtErrOverPtVsPtErrOverPt2->Fill(track->ptError() / track->pt(),track->ptError() / (track->pt()*track->pt()), EventWeight_);
      tuple->BefPreS_PtErrOverPtVsPt->Fill(track->ptError() / track->pt(), track->pt(), EventWeight_);
      tuple->BefPreS_TIsol->Fill(IsoTK_SumEt, EventWeight_);
      tuple->BefPreS_Ih->Fill(globalIh_, EventWeight_);
      tuple->BefPreS_Ias->Fill(globalIas_, EventWeight_);
      tuple->BefPreS_MassT->Fill(massT, EventWeight_);
      tuple->BefPreS_MassT_highMassT->Fill(massT, EventWeight_);
        // Add PFCadidate based isolation info to the tuple
        // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/
        // PhysicsTools/NanoAOD/plugins/IsoValueMapProducer.cc#L157
      tuple->BefPreS_MiniRelIsoAll->Fill(miniRelIsoAll, EventWeight_);
      tuple->BefPreS_MiniRelIsoChg->Fill(miniRelIsoChg, EventWeight_);
      tuple->BefPreS_MiniRelTkIso->Fill(track_genTrackMiniIsoSumPt / track->pt(), EventWeight_);
      tuple->BefPreS_MiniTkIso->Fill(track_genTrackMiniIsoSumPt, EventWeight_);
      tuple->BefPreS_SegSep->Fill(segSep, EventWeight_);
      tuple->BefPreS_SegMinPhiSep->Fill(minPhi, EventWeight_);
      tuple->BefPreS_SegMinEtaSep->Fill(minEta, EventWeight_);
      tuple->BefPreS_OpenAngle->Fill(maxOpenAngle, EventWeight_);
      tuple->BefPreS_MassErr->Fill(MassErr, EventWeight_);
      tuple->BefPreS_ProbQVsIas->Fill(1 - probQonTrack, globalIas_, EventWeight_);
      tuple->BefPreS_EtaVsIas->Fill(track->eta(), globalIas_, EventWeight_);
      tuple->BefPreS_EtaVsIh->Fill(track->eta(), globalIh_, EventWeight_);
      tuple->BefPreS_EtaVsP->Fill(track->eta(), track->p(), EventWeight_);
      tuple->BefPreS_EtaVsPt->Fill(track->eta(), track->pt(), EventWeight_);
      tuple->BefPreS_PVsIas->Fill(track->p(), globalIas_, EventWeight_);
      tuple->BefPreS_IhVsIas->Fill(globalIh_, globalIas_, EventWeight_);
      tuple->BefPreS_PVsIh->Fill(track->p(), globalIh_, EventWeight_);
      tuple->BefPreS_PtVsIas->Fill(track->pt(), globalIas_, EventWeight_);
      tuple->BefPreS_PtVsIh->Fill(track->pt(), globalIh_, EventWeight_);
      tuple->BefPreS_CaloJetsNum->Fill(caloJetsNum,EventWeight_);
      tuple->BefPreS_dRMinCaloJet->Fill(dRMinCaloJet, EventWeight_);
      tuple->BefPreS_dRMinCaloJetVsIas->Fill(dRMinCaloJet, globalIas_, EventWeight_);
      tuple->BefPreS_dRMinPfJet->Fill(dRMinPfJetWithCuts, EventWeight_);
      tuple->BefPreS_dRMinPfJetVsIas->Fill(dRMinPfJetWithCuts, globalIas_, EventWeight_);
      tuple->BefPreS_RecoPfHT->Fill(pfJetHT, EventWeight_);
      tuple->BefPreS_RecoPfJetsNum->Fill(pfJetsNum, EventWeight_);
      if (genBeta >= 0) {
        tuple->BefPreS_GenBeta->Fill(genBeta, EventWeight_);
      }
    }

    //Plotting segment separation depending on whether track passed dz cut
    if (fabs(dz) > globalMaxDZ_) {
      tuple->BefPreS_SegMinEtaSep_FailDz->Fill(minEta, EventWeight_);
    } else {
      tuple->BefPreS_SegMinEtaSep_PassDz->Fill(minEta, EventWeight_);
    }
      //Plots for tracking failing Eta Sep cut
    if (fabs(minEta) < minSegEtaSep) {
        //Needed to compare dz distribution of cosmics in pure cosmic and main sample
      tuple->BefPreS_Dz_FailSep->Fill(dz);
    }

    if (tof) {
        //Plots for tracks in dz control region
      if (fabs(dz) > CosmicMinDz && fabs(dz) < CosmicMaxDz) {
        tuple->BefPreS_Pt_FailDz->Fill(track->pt(), EventWeight_);
        tuple->BefPreS_TOF_FailDz->Fill(tof->inverseBeta(), EventWeight_);
        if (fabs(track->eta()) > CSCRegion) {
          tuple->BefPreS_TOF_FailDz_CSC->Fill(tof->inverseBeta(), EventWeight_);
          tuple->BefPreS_Pt_FailDz_CSC->Fill(track->pt(), EventWeight_);
        } else if (fabs(track->eta()) < DTRegion) {
          tuple->BefPreS_TOF_FailDz_DT->Fill(tof->inverseBeta(), EventWeight_);
          tuple->BefPreS_Pt_FailDz_DT->Fill(track->pt(), EventWeight_);
        }
      }
        //Plots of dz
      if (fabs(track->eta()) > CSCRegion) {
        tuple->BefPreS_Dz_CSC->Fill(dz, EventWeight_);
      } else if (fabs(track->eta()) < DTRegion) {
        tuple->BefPreS_Dz_DT->Fill(dz, EventWeight_);
      }
      tuple->BefPreS_EtaVsTOF->Fill(track->eta(), tof->inverseBeta(), EventWeight_);
    }
    
    // Check if HSCP is compatible with cosmics
    bool isDxyHigh = ((typeMode_ != 3 && fabs(dxy) > globalMaxDXY_) || (typeMode_ == 3 && fabs(dxyFromBeamSpot) > globalMaxDXY_)) ? true : false;
    bool isDzHigh = ((typeMode_ != 3 && fabs(dz) > globalMaxDZ_) || (typeMode_ == 3 && fabs(dzFromBeamSpot) > globalMaxDZ_)) ? true : false;
    bool isOpenAngleHigh =  (maxOpenAngle >= 2.8) ? true : false;

    isCosmicSB = (isDxyHigh && isDzHigh && isOpenAngleHigh);
    isSemiCosmicSB = (!isCosmicSB && (isDxyHigh || isDzHigh || isOpenAngleHigh));

    // TODO: have histos that shows if we have these cosmics ones (BefPreS should show non-zero
    // while PostPreS should not contain any, preselection cuts this
    if (isDzHigh && isOpenAngleHigh) {
      tuple->BefPreS_Dxy_Cosmic->Fill(dxy, EventWeight_);
    }
    if (isDxyHigh && isOpenAngleHigh) {
      tuple->BefPreS_Dz_Cosmic->Fill(dz, EventWeight_);
    }
    if (isDxyHigh && isDzHigh) {
      tuple->BefPreS_OpenAngle_Cosmic->Fill(maxOpenAngle, EventWeight_);
    }

    if (isDxyHigh && isDzHigh && isOpenAngleHigh) {
      tuple->BefPreS_Pt_Cosmic->Fill(track->pt(), EventWeight_);
      tuple->BefPreS_Ias_Cosmic->Fill(globalIas_, EventWeight_);
      tuple->BefPreS_Ih_Cosmic->Fill(globalIh_, EventWeight_);
    }
      // Get the location of the outmost hit

    tuple->BefPreS_LastHitDXY->Fill(furthersHitDxy, EventWeight_);
    tuple->BefPreS_LastHitD3D->Fill(furthersHitDistance, EventWeight_);

    if (fabs(track->eta()) < DTRegion) {
      tuple->BefPreS_Pt_DT->Fill(track->pt(), EventWeight_);
    } else {
      tuple->BefPreS_Pt_CSC->Fill(track->pt(), EventWeight_);
    }

    if (tof) {
      tuple->BefPreS_TOF->Fill(tof->inverseBeta(), EventWeight_);
      if (PUA)
        tuple->BefPreS_TOF_PUA->Fill(tof->inverseBeta(), EventWeight_);
      if (PUB)
        tuple->BefPreS_TOF_PUB->Fill(tof->inverseBeta(), EventWeight_);
      if (dttof->nDof() > 6)
        tuple->BefPreS_TOF_DT->Fill(dttof->inverseBeta(), EventWeight_);
      if (csctof->nDof() > 6)
        tuple->BefPreS_TOF_CSC->Fill(csctof->inverseBeta(), EventWeight_);
      tuple->BefPreS_PtVsTOF->Fill(track->pt(), tof->inverseBeta(), EventWeight_);
      tuple->BefPreS_TOFVsIs->Fill(tof->inverseBeta(), globalIas_, EventWeight_);
      tuple->BefPreS_TOFVsIh->Fill(tof->inverseBeta(), globalIas_, EventWeight_);
    }

    // Muon only prediction binned depending on where in the detector the track is and how many muon stations it has
    // Binning not used for other analyses
    int bin = -1;
    if (typeMode_ == 3) {
      if (fabs(track->eta()) < DTRegion) {
        bin = muonStations(track->hitPattern()) - 2;
      } else {

        bin = muonStations(track->hitPattern()) + 1;
      }
      tuple->BefPreS_Pt_Binned[to_string(bin)]->Fill(track->pt(), EventWeight_);
    }
    // ----------------------------------------------------------------------------|
    // Define preselection cuts                                                    |
    // ----------------------------------------------------------------------------|
    // int sizeOfpassedCutsArrays = 15;
    // bool passedCutsArray[sizeOfpassedCutsArrays]; --> this didnt work out, TODO come back to this?
    bool passedCutsArray[15];
    std::fill(std::begin(passedCutsArray), std::end(passedCutsArray),false);

    // No cut, i.e. events after trigger
    passedCutsArray[0]  = (trigInfo_ > 0) ? true : false;
    // Cut on transverse momentum
    // Single muon trigger threshold is 50 GeV
    passedCutsArray[1]  = (track->pt() > globalMinPt_) ? true : false;
    // Check if eta is inside the max eta cut for detector homogeneity
    passedCutsArray[2]  = (fabs(track->eta()) < globalMaxEta_) ? true : false;
    // Check the number of non-layer-1 pixel hits to ensure good stats on the hits
    passedCutsArray[3]  = (typeMode_ != 3 && nonL1PixHits >= globalMinNOPH_) ? true : false;
    // Check the min fraction of valid hits to ensure good stats on the hits
    passedCutsArray[4]  = (typeMode_ != 3 && track->validFraction() > globalMinFOVH_) ? true : false;
    // Cut for the number of dEdx hits to ensure good stats on the hits
    passedCutsArray[5]  = (numDeDxHits >= globalMinNOM_)  ? true : false;
    // Select only high purity tracks to ensure good quality tracks
    passedCutsArray[6]  = (typeMode_ != 3 && track->quality(reco::TrackBase::highPurity)) ? true : false;
    // Cut on the chi2 / ndof to ensure good quality tracks
    passedCutsArray[7] = (typeMode_ != 3 && (track->chi2() / track->ndof()) < globalMaxChi2_) ? true : false;
    // Cut on the impact parameter to ensure the track is coming from the PV
    // for typeMode_ 3 (TOF only) dz and dxy is supposed to come from the beamspot
    passedCutsArray[8] = (  (typeMode_ != 3 && fabs(dz) < globalMaxDZ_)
                          || (typeMode_ == 3 && fabs(dzFromBeamSpot) < 4)) ? true : false;
    passedCutsArray[9] = (  (typeMode_ != 3 && fabs(dxy) < globalMaxDXY_)
                          || (typeMode_ == 3 && fabs(dxyFromBeamSpot) < 4)) ? true : false;
    // Cut on the PF based mini-isolation
    passedCutsArray[10] = ( miniRelIsoAll < globalMaxMiniRelIsoAll_ ) ? true : false;
    // Cut on the absolute pT-dependent cone size TkIsolation
    passedCutsArray[11] = ( track_genTrackMiniIsoSumPt < globalMaxTIsol_ ) ? true : false;
    // Cut on the energy over momenta
    passedCutsArray[12] = (EoP < globalMaxEoP_) ? true : false;
    // Cut on the uncertainty of the pt measurement, TODO the cut value should be made a global variable
    passedCutsArray[13] = (typeMode_ != 3 && (track->ptError() / (track->pt()*track->pt()) < 0.0008)) ? true : false;
    // Cut away background events based on the probQ
    passedCutsArray[14] = (probQonTrackNoL1 < globalMaxTrackProbQCut_ && probQonTrackNoL1 > globalMinTrackProbQCut_) ? true : false;
    //  passedCutsArray[13] = (typeMode_ != 3 && (track->ptError() / track->pt()) < pTerr_over_pT_etaBin(track->pt(), track->eta())) ? true : false;
    // Cut on the tracker based isolation
    //  passedCutsArray[12] = ( IsoTK_SumEt < globalMaxTIsol_) ? true : false;
  
    // Cut on the PF electron ID
    //  passedCutsArray[14] = ( !pf_isElectron  && !pf_isPhoton) ? true : false;
    // Cut on min Ih (or max for fractionally charged)
    //  passedCutsArray[15] = (  (typeMode_ != 5 &&  globalIh_ > globalMinIh_)
    //                        || (typeMode_ == 5 && globalIh_ < globalMinIh_)) ? true : false;
    //passedCutsArray[16] = ( MassErr < 3 ) ? true : false;
    // Cut away background events based on the probXY
    //  passedCutsArray[16] = ((probXYonTrackNoL1 > globalMinTrackProbXYCut_) && (probXYonTrackNoL1 < globalMaxTrackProbXYCut_))  ? true : false;

    //  // TOF only cuts
    //  passedCutsArray[18] = (typeMode_ != 3 || (typeMode_ == 3 && muonStations(track->hitPattern()) > minMuStations_)) ? true : false;
    //  passedCutsArray[19] = (typeMode_ != 3 || (typeMode_ == 3 && fabs(track->phi()) > 1.2 && fabs(track->phi()) < 1.9)) ? true : false;
    //  passedCutsArray[20] = (typeMode_ != 3 || (typeMode_ == 3 && fabs(minEta) > minSegEtaSep)) ? true : false;
    //
    // Not used cuts TODO: revise
    // cut on the number of missing hits from IP till last hit (excluding hits behind the last hit)
    // bool cutMinNumOfMissingHits = (typeMode_ != 3 && missingHitsTillLast > GlobalMaxNOMHTillLast) ? true : false;
    // cut on the fraction of valid hits divided by total expected hits until the last one
    // bool cutMinFractOfValidHitsTillLast = (typeMode_ != 3 && validFractionTillLast < GlobalMinFOVHTillLast) ? true : false;
    // cut on relative tracker isolation (SumPt/Pt)
    // bool cutRelTKIso = ( IsoTK_SumEt / track->pt() > GlobalMaxRelTIsol)  ? true : false;
    // Cut for number of DOF in TOF ana
    
    // Define preselection cuts for Sept 8 preselection
    bool passedCutsArraySept8[10];
    std::fill(std::begin(passedCutsArraySept8), std::end(passedCutsArraySept8),false);

    // To reduce mistakes, we just copy the first 10 cuts from earlier, those are the frozen cuts
    passedCutsArraySept8[0]  = passedCutsArray[0];
    passedCutsArraySept8[1]  = passedCutsArray[1];
    passedCutsArraySept8[2]  = passedCutsArray[2];
    passedCutsArraySept8[3]  = passedCutsArray[3];
    passedCutsArraySept8[4]  = passedCutsArray[4];
    passedCutsArraySept8[5]  = passedCutsArray[5];
    passedCutsArraySept8[6]  = passedCutsArray[6];
    passedCutsArraySept8[7] = passedCutsArray[7];
    passedCutsArraySept8[8] = passedCutsArray[8];
    passedCutsArraySept8[9] = passedCutsArray[9];

    // N-1 plots
    for (size_t i=0;i<sizeof(passedCutsArray);i++) {
      bool allOtherCutsPassed = true;
      if (passedCutsArray[i] && trigInfo_ > 0) {
        tuple->N1_FlowEta->Fill(track->eta(), i, EventWeight_);
      }
      for (size_t j=0;j<sizeof(passedCutsArray);j++) {
        if (i==j) continue;
        if (!passedCutsArray[j]) {
          allOtherCutsPassed = false;
          // We found a cut that's not passed, no point in looking into the rest of them
          break;
        }
      }
      if (allOtherCutsPassed && trigInfo_ > 0) {
          // Put the not used variables to the i==0, this will be always true
        if (i==0)  {
          tuple->N1_PfType->Fill(0., EventWeight_);
          if (pf_isPfTrack) {
            tuple->N1_PfType->Fill(1., EventWeight_);
          } else {
            tuple->N1_PfType->Fill(8., EventWeight_);
          }
          if (pf_isElectron) {
            tuple->N1_PfType->Fill(2., EventWeight_);
          } else if (pf_isMuon) {
            tuple->N1_PfType->Fill(3., EventWeight_);
          } else if (pf_isPhoton) {
            tuple->N1_PfType->Fill(4., EventWeight_);
          } else if (pf_isChHadron) {
            tuple->N1_PfType->Fill(5., EventWeight_);
          } else if (pf_isNeutHadron) {
            tuple->N1_PfType->Fill(6., EventWeight_);
          } else if (pf_isUndefined) {
            tuple->N1_PfType->Fill(7., EventWeight_);
          }
          tuple->N1_Ih->Fill(globalIh_, EventWeight_);
          tuple->N1_ProbXY->Fill(probXYonTrack, EventWeight_);
          tuple->N1_Stations->Fill(muonStations(track->hitPattern()), EventWeight_);
          tuple->N1_DrMinPfJet->Fill(dRMinPfJetWithCuts, EventWeight_);
          tuple->N1_SumpTOverpT->Fill(IsoTK_SumEt / track->pt(), EventWeight_);
        };

        if (i==1)  {
          tuple->N1_Pt->Fill(track->pt(), EventWeight_);
          tuple->N1_Pt_lowPt->Fill(track->pt(), EventWeight_);
        };
        if (i==2)  { tuple->N1_Eta->Fill(track->eta(), EventWeight_); };
        if (i==3)  { tuple->N1_TNOPH->Fill(nonL1PixHits, EventWeight_); };
        if (i==4)  { tuple->N1_TNOHFraction->Fill(track->validFraction(), EventWeight_); };
        if (i==5)  { tuple->N1_TNOM->Fill(numDeDxHits, EventWeight_); };
        if (i==6)  {
          if (track->quality(reco::TrackBase::highPurity)) {
            tuple->N1_Qual->Fill(1., EventWeight_);
          } else {
            tuple->N1_Qual->Fill(0., EventWeight_);
          }
        };
        if (i==7) { tuple->N1_Chi2oNdof->Fill(track->chi2() / track->ndof(), EventWeight_); };
        if (i==8) { tuple->N1_Dz->Fill(dz, EventWeight_); };
        if (i==9) { tuple->N1_Dxy->Fill(dxy, EventWeight_); };
        if (i==10) {
          tuple->N1_MiniRelIsoAll->Fill(miniRelIsoAll, EventWeight_);
          tuple->N1_MiniRelIsoAll_lowMiniRelIso->Fill(miniRelIsoAll, EventWeight_);
        }
        if (i==11) {
          tuple->N1_MiniTkIso->Fill(track_genTrackMiniIsoSumPt, EventWeight_);

          tuple->N1_MiniRelTkIso->Fill(track_genTrackMiniIsoSumPt / track->pt(), EventWeight_);
          tuple->N1_MiniRelTkIso_lowMiniRelIso->Fill(track_genTrackMiniIsoSumPt / track->pt(), EventWeight_);
          if (PUA) {
            tuple->N1_MiniTkIso_PUA->Fill(track_genTrackMiniIsoSumPt, EventWeight_);
            tuple->N1_MiniRelTkIso_lowMiniRelIso_PUA->Fill(track_genTrackMiniIsoSumPt / track->pt(), EventWeight_);
          } else if (PUB) {
            tuple->N1_MiniTkIso_PUB->Fill(track_genTrackMiniIsoSumPt, EventWeight_);
            tuple->N1_MiniRelTkIso_lowMiniRelIso_PUB->Fill(track_genTrackMiniIsoSumPt / track->pt(), EventWeight_);
          } else if (PUC) {
            tuple->N1_MiniTkIso_PUC->Fill(track_genTrackMiniIsoSumPt, EventWeight_);
            tuple->N1_MiniRelTkIso_lowMiniRelIso_PUC->Fill(track_genTrackMiniIsoSumPt / track->pt(), EventWeight_);
          }
        };
        if (i==12) {
          tuple->N1_EoP->Fill(EoP, EventWeight_);
        };
        if (i==13) {
          tuple->N1_PtErrOverPt->Fill(track->ptError() / track->pt(), EventWeight_);
          tuple->N1_PtErrOverPt2->Fill(track->ptError() / (track->pt()*track->pt()), EventWeight_);
          tuple->N1_PtErrOverPtVsPt->Fill(track->ptError() / track->pt(), track->pt(), EventWeight_);
          tuple->N1_PtErrOverPtVsPt_lowPt->Fill(track->ptError() / track->pt(), track->pt(), EventWeight_);
          tuple->N1_PtErrOverPtVsGenBeta->Fill(track->ptError() / track->pt(), genBeta, EventWeight_);
          tuple->N1_PtErrOverPt2VsIas->Fill(track->ptError() / (track->pt()*track->pt()), globalIas_, EventWeight_);
          tuple->N1_PtErrOverPt2VsProbQNoL1->Fill(track->ptError() / (track->pt()*track->pt()), 1 - probQonTrackNoL1, EventWeight_);
        };
        if (i==14) {
          tuple->N1_ProbQNoL1->Fill(1 - probQonTrackNoL1, EventWeight_);
          tuple->N1_ProbQNoL1VsIas->Fill(1 - probQonTrackNoL1, globalIas_, EventWeight_);
        };
      }
    }

    // All track after technical checks
    tuple->CutFlow->Fill(1., EventWeight_);
    // CutFlow in a single plot
    if (trigInfo_ > 0) {
      for (size_t i=0;i<sizeof(passedCutsArray);i++) {
        bool allCutsPassedSoFar = true;
        for (size_t j=0;j<=i;j++) {
          if (!passedCutsArray[j]) {
            allCutsPassedSoFar = false;
            break;
          }
        }
        if (allCutsPassedSoFar) {
          if (bestSoFarCandCutInd < i) {
            bestSoFarCandCutInd = i;
          }
          tuple->CutFlow->Fill((i+2), EventWeight_);
        }
      }
    }
    
    // Plot Eta after each cut
    for (size_t i=0;i<sizeof(passedCutsArray);i++) {
      bool allCutsPassedSoFar = true;
      for (size_t j=0;j<=i;j++) {
        if (!passedCutsArray[j]) {
          allCutsPassedSoFar = false;
          break;
        }
      }

      if (allCutsPassedSoFar) {
        tuple->CutFlowEta->Fill(track->eta(), i, EventWeight_);
        tuple->CutFlowProbQ->Fill(1 - probQonTrackNoL1, i, EventWeight_);
        tuple->CutFlowPfType->Fill(0., i, EventWeight_);
        if (pf_isPfTrack) {
          tuple->CutFlowPfType->Fill(1., i, EventWeight_);
        } else {
          tuple->CutFlowPfType->Fill(8., i, EventWeight_);
        }
        if (pf_isElectron) {
          tuple->CutFlowPfType->Fill(2., i, EventWeight_);
        } else if (pf_isMuon) {
          tuple->CutFlowPfType->Fill(3., i, EventWeight_);
        } else if (pf_isPhoton) {
          tuple->CutFlowPfType->Fill(4., i, EventWeight_);
        } else if (pf_isChHadron) {
          tuple->CutFlowPfType->Fill(5., i, EventWeight_);
        } else if (pf_isNeutHadron) {
          tuple->CutFlowPfType->Fill(6., i, EventWeight_);
        } else if (pf_isUndefined) {
          tuple->CutFlowPfType->Fill(7., i, EventWeight_);
        }
      }
    }

    // Reverse cutflow, i.e. start with the last cut from the original cutflow
    bool passedCutsArrayReverse[15];
    std::reverse_copy(std::begin(passedCutsArray), std::end(passedCutsArray), std::begin(passedCutsArrayReverse));
    for (size_t i=0;i<sizeof(passedCutsArrayReverse);i++) {
      bool allCutsPassedSoFar = true;
      for (size_t j=0;j<=i;j++) {
        if (!passedCutsArrayReverse[j]) {
          allCutsPassedSoFar = false;
        }
      }
      if (allCutsPassedSoFar) {
        tuple->CutFlowReverse->Fill((i), EventWeight_);
      }
    }

    // Preselection cuts for a CR where the pT cut is flipped
    bool passedCutsArrayForCR[15];
    std::copy(std::begin(passedCutsArray), std::end(passedCutsArray), std::begin(passedCutsArrayForCR));
    passedCutsArrayForCR[1] = (track->pt() > 50 && track->pt() < 55) ? true : false;
    
    if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "\n      >> Check if we pass Preselection for CR";
    if (passPreselection(passedCutsArrayForCR, false)) {
      tuple->PostPreS_Ias_CR->Fill(globalIas_, EventWeight_);
      tuple->PostPreS_Pt_lowPt_CR->Fill(track->pt(), EventWeight_);
      tuple->PostPreS_ProbQNoL1_CR->Fill(1 - probQonTrackNoL1, EventWeight_);
      tuple->PostPreS_ProbQNoL1VsIas_CR->Fill(1 - probQonTrackNoL1, globalIas_, EventWeight_);
      tuple->PostPreS_ProbQNoL1VsIas_CR_Pileup_up->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_ * PUSystFactor_[0]);
      tuple->PostPreS_ProbQNoL1VsIas_CR_Pileup_down->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_ * PUSystFactor_[1]);
      tuple->PostPreS_ProbQNoL1VsIas_CR_ProbQNoL1_up->Fill(std::min(1.,(1 - probQonTrackNoL1)*1.005), globalIas_,  EventWeight_);
      tuple->PostPreS_ProbQNoL1VsIas_CR_ProbQNoL1_down->Fill((1 - probQonTrackNoL1)*0.995, globalIas_,  EventWeight_);
      tuple->PostPreS_ProbQNoL1VsIas_CR_Ias_up->Fill(1 - probQonTrackNoL1, std::min(1.0,globalIas_*1.05),  EventWeight_);
      tuple->PostPreS_ProbQNoL1VsIas_CR_Ias_down->Fill(1 - probQonTrackNoL1, globalIas_*0.95,  EventWeight_);
      tuple->PostPreS_ProbQNoL1VsIas_CR_Pt_up->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_);
      tuple->PostPreS_ProbQNoL1VsIas_CR_Pt_down->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_);
    }

    bool passedCutsArrayForGiTemplates[15];
    std::copy(std::begin(passedCutsArray), std::end(passedCutsArray), std::begin(passedCutsArrayForGiTemplates));

    //Preselection steps that can be changed for template generation
    passedCutsArrayForGiTemplates[0] = true;
    passedCutsArrayForGiTemplates[1] = ( (track->p() > 20) && (track->p() < 48) );
    passedCutsArrayForGiTemplates[2]  = (fabs(track->eta()) < 2.1);

    passedCutsArrayForGiTemplates[10] = true; //miniIso
    //[11] is TrackerIso > 15
    passedCutsArrayForGiTemplates[12] = true; //EoP
    passedCutsArrayForGiTemplates[13] = true; //sigma pt / pt2 cut
    passedCutsArrayForGiTemplates[14] = true; //ProbQ cut

    // PAY ATTENTION : PU == NPV for the following loop !
    if (passPreselection(passedCutsArrayForGiTemplates, false)) {
      tuple->PostPreS_Ias_CR_veryLowPt->Fill(globalIas_,preScaleForDeDx*EventWeight_);
      for(unsigned int h=0; h< dedxHits->size(); h++){
        DetId detid(dedxHits->detId(h));
        int modulgeomForIndxH = 0;
        float pathlenghtForIndxH = 0.0;
        float chargeForIndxH = 0.0;
        float normMult = 265;
        float scaleFactor = dEdxSF_0_;
        if (detid.subdetId() < 3) {
          modulgeomForIndxH = 15;
        } else {
          SiStripDetId SSdetId(detid);
          modulgeomForIndxH = SSdetId.moduleGeometry();
        }
        pathlenghtForIndxH = dedxHits->pathlength(h) * 10;
        chargeForIndxH = dedxHits->charge(h);
        if(createGiTemplates_){
            if(puTreatment_){
               int npv = vertexColl.size();
               for(int i = 0 ; i < NbPuBins_ ; i++){
                   if ( npv > PuBins_[i] && npv <= PuBins_[i+1] ){
                       //std::cout << "We are then in PU bin #" << (i+1) << std::endl;
                       if(detid.subdetId() >= 3){ // For strips only 1 question : charge is already multiplied by sclae factor ?
                           if(i==0) tuple->Calibration_GiTemplate_PU_1->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/pathlenghtForIndxH, preScaleForDeDx);
                           else if(i==1) tuple->Calibration_GiTemplate_PU_2->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/pathlenghtForIndxH, preScaleForDeDx);
                           else if(i==2) tuple->Calibration_GiTemplate_PU_3->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/pathlenghtForIndxH, preScaleForDeDx);
                           else if(i==3) tuple->Calibration_GiTemplate_PU_4->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/pathlenghtForIndxH, preScaleForDeDx);
                           else if(i==4) tuple->Calibration_GiTemplate_PU_5->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/pathlenghtForIndxH, preScaleForDeDx);
                       }
                       else{ // for pixels , we multiply pathlenght by normMult
                           scaleFactor *= dEdxSF_1_;
                           if(i==0) tuple->Calibration_GiTemplate_PU_1->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/(pathlenghtForIndxH*normMult), preScaleForDeDx);
                           else if(i==1) tuple->Calibration_GiTemplate_PU_2->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/(pathlenghtForIndxH*normMult), preScaleForDeDx);
                           else if(i==2) tuple->Calibration_GiTemplate_PU_3->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/(pathlenghtForIndxH*normMult), preScaleForDeDx);
                           else if(i==3) tuple->Calibration_GiTemplate_PU_4->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/(pathlenghtForIndxH*normMult), preScaleForDeDx);
                           else if(i==4) tuple->Calibration_GiTemplate_PU_5->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/(pathlenghtForIndxH*normMult), preScaleForDeDx);

                       }
                   }
               }
            } else {
              // Standard approach, no PU dependence
              //We should have pixel and strips in different loops : in pixels, we divide charge * SF / pathlenght * normMult
              if(detid.subdetId() >= 3){
                  tuple->Calibration_GiTemplate->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/pathlenghtForIndxH, preScaleForDeDx);
                  tuple->Calibration_GiTemplate_noL1->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/pathlenghtForIndxH, preScaleForDeDx);
              }
              else{
                  scaleFactor *= dEdxSF_1_;
                  tuple->Calibration_GiTemplate->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/(pathlenghtForIndxH*normMult), preScaleForDeDx);
                  tuple->Calibration_GiTemplate_noL1->Fill(modulgeomForIndxH, pathlenghtForIndxH, scaleFactor*chargeForIndxH/(pathlenghtForIndxH*normMult), preScaleForDeDx);
              }
            }

        }//end condition createGitemplates
      } // end loop on dEdx hits
    } // end condition on passing preselection cuts for the Gi templates


    if (createGiTemplates_ && !CreateAndRunGitemplates_) {
      // If the purpose of this is to run the Gi templates only, exit here
      return;
    }
    
    if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> Check if we pass Preselection";
    bool passPre = passPreselection(passedCutsArray, true);
    if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> Check if we pass Preselection with Sept8 cuts";
    bool passPreSept8 = passPreselection(passedCutsArraySept8, true);
    
    // Few more bins in CutFlow for SRs
    unsigned int passedCutsArraySize = sizeof(passedCutsArray);
    if (globalIas_ > 0.25 && probQonTrackNoL1 < 0.1 && passPre) {
      tuple->CutFlow->Fill(passedCutsArraySize + 2, EventWeight_);
      if (track->pt() > 100) {
        tuple->CutFlow->Fill(passedCutsArraySize+ 3, EventWeight_);
      } if (track->pt() > 200) {
        tuple->CutFlow->Fill(passedCutsArraySize + 4, EventWeight_);
      }
    }
    
    // This should be moved more up in the code
    // Dont do TOF only is isCosmicSB is true
    if (typeMode_ == 3 && isCosmicSB) {
      if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> This is a cosmic track, skipping it";
      // 8-th bin of the error histo, not a collision track
      tuple->ErrorHisto->Fill(7.);
      continue;
    } else if (isCosmicSB) {
      if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> This is a cosmic track, please check what's up";
    }

    // Dont do TOF only is isSemiCosmicSB is true
    if (typeMode_ == 3 && isSemiCosmicSB) {
      if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> This is a semi-cosmic track, skipping it";
      // 8-th bin of the error histo, not a collision track
      tuple->ErrorHisto->Fill(7.);
      continue;
    } else if (isSemiCosmicSB) {
      if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> This is a semi-cosmic track, please check what's up";
    }

    // Systematics plots for pT rescaling
    
    float shiftForPtValue = shiftForPt(track->pt(),track->eta(),track->phi(),track->charge());
    float rescaledPtUp = track->pt()*(1+shiftForPtValue);
    bool passedCutsArrayForPtSyst[15];
    std::copy(std::begin(passedCutsArray), std::end(passedCutsArray), std::begin(passedCutsArrayForPtSyst));
    passedCutsArrayForPtSyst[1] = (rescaledPtUp > globalMinPt_) ? true : false;
    if (debug_ > 2) LogPrint(MOD) << "      >> Check if we pass Preselection for pT systematics";
    if (passPreselection(passedCutsArrayForPtSyst, false)) {
      tuple->PostPreS_ProbQNoL1VsIas_Pt_up->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_);
    }

    float rescaledPtDown = track->pt()*(1-shiftForPtValue);
    passedCutsArrayForPtSyst[1] = (rescaledPtDown > globalMinPt_) ? true : false;
    if (passPreselection(passedCutsArrayForPtSyst, false)) {
      tuple->PostPreS_ProbQNoL1VsIas_Pt_down->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_);
    }

    // Systematics plots for Ias rescaling
    // TODODec20 --> Ias smearing isntead of a shift
    tuple->PostPreS_ProbQNoL1VsIas_Ias_up->Fill(1 - probQonTrackNoL1, std::min(1.0,globalIas_*1.05),  EventWeight_);
    tuple->PostPreS_ProbQNoL1VsIas_Ias_down->Fill(1 - probQonTrackNoL1, globalIas_*0.95,  EventWeight_);

    // Systematics plots for trigger
    bool passedCutsArrayForTriggerSyst[15];
    std::copy(std::begin(passedCutsArray), std::end(passedCutsArray), std::begin(passedCutsArrayForTriggerSyst));
    passedCutsArrayForTriggerSyst[0] = (true) ? true : false;
    if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> Check if we pass Preselection for trigger systematics";
    if (passPreselection(passedCutsArrayForTriggerSyst, false)) {
      float distanceAtThisEta = 500.0/sin(track->phi());
      float speedOfLightInCmPerNs = 29.97;
      float timing = distanceAtThisEta / (genBeta*speedOfLightInCmPerNs);
      
      float genBetaPrimeUp =  -1.f;
      float genBetaPrimeDown = -1.f;
      genBetaPrimeUp = distanceAtThisEta / ((timing+1.5)*speedOfLightInCmPerNs);
      if (((timing-1.5)*speedOfLightInCmPerNs) > 0.0) {
        genBetaPrimeDown = distanceAtThisEta / ((timing-1.5)*speedOfLightInCmPerNs);
      }
      if (!HLT_Mu50) {
        tuple->PostPreS_TriggerMuon50VsBeta->Fill(0., genBeta);
        tuple->PostPreS_TriggerMuon50VsPt->Fill(0., track->pt());
        if (fabs(track->eta()) < 0.3) {
          tuple->PostPreS_TriggerMuon50VsBeta_EtaA->Fill(0., genBeta);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaA_BetaUp->Fill(0., genBetaPrimeUp);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaA_BetaDown->Fill(0., genBetaPrimeDown);
        } else  if (fabs(track->eta()) > 0.3 && fabs(track->eta()) < 0.6) {
          tuple->PostPreS_TriggerMuon50VsBeta_EtaB->Fill(0., genBeta);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaB_BetaUp->Fill(0., genBetaPrimeUp);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaB_BetaDown->Fill(0., genBetaPrimeDown);
        } else  if (fabs(track->eta()) > 0.6 && fabs(track->eta()) < 1.0) {
          tuple->PostPreS_TriggerMuon50VsBeta_EtaC->Fill(0., genBeta);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaC_BetaUp->Fill(0., genBetaPrimeUp);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaC_BetaDown->Fill(0., genBetaPrimeDown);
        }

      } else if (HLT_Mu50) {
        tuple->PostPreS_TriggerMuon50VsBeta->Fill(1., genBeta);
        tuple->PostPreS_TriggerMuon50VsPt->Fill(1., track->pt());
        if (fabs(track->eta()) < 0.3) {
          tuple->PostPreS_TriggerMuon50VsBeta_EtaA->Fill(1., genBeta);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaA_BetaUp->Fill(1., genBetaPrimeUp);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaA_BetaDown->Fill(1., genBetaPrimeDown);
        } else  if (fabs(track->eta()) > 0.3 && fabs(track->eta()) < 0.6) {
          tuple->PostPreS_TriggerMuon50VsBeta_EtaB->Fill(1., genBeta);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaB_BetaUp->Fill(1., genBetaPrimeUp);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaB_BetaDown->Fill(1., genBetaPrimeDown);
        } else  if (fabs(track->eta()) > 0.6 && fabs(track->eta()) < 1.0) {
          tuple->PostPreS_TriggerMuon50VsBeta_EtaC->Fill(1., genBeta);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaC_BetaUp->Fill(1., genBetaPrimeUp);
          tuple->PostPreS_TriggerMuon50VsBeta_EtaC_BetaDown->Fill(1., genBetaPrimeDown);
        }
      }
      if (!HLT_PFMET120_PFMHT120_IDTight && !HLT_PFHT500_PFMET100_PFMHT100_IDTight && !HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && !HLT_MET105_IsoTrk50) {
        tuple->PostPreS_TriggerMETallVsBeta->Fill(0., genBeta);
        tuple->PostPreS_TriggerMETallVsMet->Fill(0., RecoPFMET);
        tuple->PostPreS_TriggerMETallVsHT->Fill(0., pfJetHT);
        tuple->PostPreS_TriggerMETallVsMetVsHT->Fill(0., RecoPFMET, pfJetHT);
        tuple->PostPreS_TriggerMETallVsMetOverHt->Fill(0.,RecoPFMET/pfJetHT);
      } else if (HLT_PFMET120_PFMHT120_IDTight || HLT_PFHT500_PFMET100_PFMHT100_IDTight || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || HLT_MET105_IsoTrk50) {
        tuple->PostPreS_TriggerMETallVsBeta->Fill(1., genBeta);
        tuple->PostPreS_TriggerMETallVsMet->Fill(1., RecoPFMET);
        tuple->PostPreS_TriggerMETallVsHT->Fill(1., pfJetHT);
        tuple->PostPreS_TriggerMETallVsMetVsHT->Fill(1., RecoPFMET, pfJetHT);
        tuple->PostPreS_TriggerMETallVsMetOverHt->Fill(1., RecoPFMET/pfJetHT);
      }
    } // end of passedCutsArrayForTriggerSyst loop
    

    //fill the ABCD histograms and a few other control plots
    if (passPre) {
      if (trigInfo_ < 1) {
        cout << "we are pass preS but the trigger is not passed -- this is wrong!!!" << endl << endl << endl << endl ;
      }
      tuple->PostPreS_MetVsHT->Fill(RecoPFMET, pfJetHT, EventWeight_);
      tuple->PostPreS_MetOverHt->Fill(RecoPFMET/pfJetHT, EventWeight_);
      if (trigInfo_ > 0) postPreS_candidate_count++;
      if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> Passed pre-selection";
      if (debug_ > 2 && trigInfo_ > 0) LogPrint(MOD) << "      >> Fill control and prediction histos";
      
      tuple_maker->fillControlAndPredictionHist(hscp,
                                                dedxSObj,
                                                dedxMObj,
                                                tof,
                                                tuple,
                                                typeMode_,
                                                globalMinTOF_ ,
                                                EventWeight_,
                                                isCosmicSB,
                                                DTRegion,
                                                MaxPredBins,
                                                dEdxK_,
                                                dEdxC_,
                                                CutPt_,
                                                CutI_,
                                                CutTOF_,
                                                CutPt_Flip_,
                                                CutI_Flip_,
                                                CutTOF_Flip_);
      // After pre-selection plots
      // TODO: controll this from config
      bool doPostPreSplots = true;
      if (doPostPreSplots && trigInfo_ > 0) {
        if (debug_ > 3) LogPrint(MOD) << "      >> Fill post preselection histos";
        tuple->PostPreS_MuonPtVsTrackPt->Fill(muonPt, trackerPt, EventWeight_);
        tuple->PostPreS_MuonPtOverGenPtVsTrackPtOverGenPt->Fill(muonPt/genPt, trackerPt/genPt, EventWeight_);
        tuple->PostPreS_RelDiffMuonPtAndTrackPt->Fill((muonPt-trackerPt)/trackerPt, EventWeight_);
        tuple->PostPreS_PfType->Fill(0., EventWeight_);
        tuple->PostPreS_PfTypeVsIas->Fill(0., globalIas_, EventWeight_);
        if (pf_isPfTrack) {
          tuple->PostPreS_PfType->Fill(1., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(1., globalIas_, EventWeight_);
        } else {
          tuple->PostPreS_PfType->Fill(8., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(8., globalIas_, EventWeight_);
        }
        if (pf_isElectron) {
          tuple->PostPreS_PfType->Fill(2., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(2., globalIas_, EventWeight_);
        } else if (pf_isMuon) {
          tuple->PostPreS_PfType->Fill(3., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(3., globalIas_, EventWeight_);
        } else if (pf_isPhoton) {
          tuple->PostPreS_PfType->Fill(4., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(4., globalIas_, EventWeight_);
        } else if (pf_isChHadron) {
          tuple->PostPreS_PfType->Fill(5., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(5., globalIas_, EventWeight_);
        } else if (pf_isNeutHadron) {
          tuple->PostPreS_PfType->Fill(6., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(6., globalIas_, EventWeight_);
        } else if (pf_isUndefined) {
          tuple->PostPreS_PfType->Fill(7., EventWeight_);
          tuple->PostPreS_PfTypeVsIas->Fill(7., globalIas_, EventWeight_);
        }
        tuple->PostPreS_Eta->Fill(track->eta(), EventWeight_);
        tuple->PostPreS_EtaVsIas->Fill(track->eta(), globalIas_, EventWeight_);
        tuple->PostPreS_MatchedStations->Fill(muonStations(track->hitPattern()), EventWeight_);
        tuple->PostPreS_NVertex->Fill(vertexColl.size(), EventWeight_);
        tuple->PostPreS_NVertex_NoEventWeight->Fill(vertexColl.size());
        if (PUA) {
          tuple->PostPreS_TNOH_PUA->Fill(track->found(), EventWeight_);
          tuple->PostPreS_TNOM_PUA->Fill(numDeDxHits, EventWeight_);
        }
        else if (PUB) {
          tuple->PostPreS_TNOH_PUB->Fill(track->found(), EventWeight_);
          tuple->PostPreS_TNOM_PUB->Fill(numDeDxHits, EventWeight_);
        }
        else if (PUC) {
          tuple->PostPreS_TNOH_PUC->Fill(track->found(), EventWeight_);
          tuple->PostPreS_TNOM_PUC->Fill(numDeDxHits, EventWeight_);
        }
        tuple->PostPreS_TNOHFraction->Fill(track->validFraction(), EventWeight_);
        tuple->PostPreS_TNOHFractionVsIas->Fill(track->validFraction(), globalIas_, EventWeight_);
        tuple->PostPreS_TNOPH->Fill(nonL1PixHits, EventWeight_);
        tuple->PostPreS_RatioCleanAndAllStripsClu->Fill(ratioCleanAndAllStripsClu, EventWeight_);
        tuple->PostPreS_RatioCleanAndAllStripsCluVsIas->Fill(ratioCleanAndAllStripsClu, globalIas_, EventWeight_);
        tuple->PostPreS_RatioCleanAndAllPixelClu->Fill(ratioCleanAndAllPixelClu, EventWeight_);
        tuple->PostPreS_TNOPHVsIas->Fill(nonL1PixHits, globalIas_, EventWeight_);
        tuple->PostPreS_TNOHFractionTillLast->Fill(validFractionTillLast, EventWeight_);
        tuple->PostPreS_TNOMHTillLast->Fill(missingHitsTillLast, EventWeight_);
        tuple->PostPreS_TNOM->Fill(numDeDxHits, EventWeight_);
        tuple->PostPreS_TNOMVsIas->Fill(numDeDxHits, globalIas_, EventWeight_);
        tuple->PostPreS_ProbQ->Fill(1 - probQonTrack, EventWeight_);
        tuple->PostPreS_ProbQVsIas->Fill(1 - probQonTrack, globalIas_, EventWeight_);
        tuple->PostPreS_IhVsProbQNoL1VsIas->Fill(globalIh_, 1 - probQonTrackNoL1, globalIas_, EventWeight_);
        tuple->PostPreS_MomentumVsProbQNoL1VsIas->Fill(track->p(), 1 - probQonTrackNoL1, globalIas_, EventWeight_);
        tuple->PostPreS_ProbXY->Fill(probXYonTrack, EventWeight_);
        tuple->PostPreS_ProbXYVsIas->Fill(probXYonTrack, globalIas_, EventWeight_);
        tuple->PostPreS_ProbXYVsProbQ->Fill(probXYonTrack, 1 - probQonTrack, EventWeight_);
        tuple->PostPreS_ProbQNoL1->Fill(1 - probQonTrackNoL1, EventWeight_);
        tuple->PostPreS_ProbQNoL1VsIas->Fill(1 - probQonTrackNoL1, globalIas_, EventWeight_);
        tuple->PostPreS_ProbQNoL1VsIas_Pileup_up->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_ * PUSystFactor_[0]);
        tuple->PostPreS_ProbQNoL1VsIas_Pileup_down->Fill(1 - probQonTrackNoL1, globalIas_,  EventWeight_ * PUSystFactor_[1]);
        tuple->PostPreS_ProbQNoL1VsIas_ProbQNoL1_up->Fill(std::min(1.0,(1 - probQonTrackNoL1)*1.005), globalIas_,  EventWeight_);
        tuple->PostPreS_ProbQNoL1VsIas_ProbQNoL1_down->Fill((1 - probQonTrackNoL1)*0.995, globalIas_,  EventWeight_);

        tuple->PostPreS_ProbXYNoL1->Fill(probXYonTrackNoL1, EventWeight_);
        tuple->PostPreS_ProbXYNoL1VsIas->Fill(probXYonTrackNoL1, globalIas_, EventWeight_);
        tuple->PostPreS_ProbXYNoL1VsProbQNoL1->Fill(probXYonTrackNoL1, 1 - probQonTrackNoL1, EventWeight_);

        if (globalIas_ > 0.3) {
          tuple->PostPreS_ProbXY_highIas->Fill(probXYonTrack, EventWeight_);
          tuple->PostPreS_ProbXYVsIas_highIas->Fill(probXYonTrack, globalIas_, EventWeight_);
          tuple->PostPreS_ProbXYVsProbQ_highIas->Fill(probXYonTrack, 1 - probQonTrack, EventWeight_);
          tuple->PostPreS_ProbXYNoL1_highIas->Fill(probXYonTrackNoL1, EventWeight_);
          tuple->PostPreS_ProbXYNoL1VsIas_highIas->Fill(probXYonTrackNoL1, globalIas_, EventWeight_);
          tuple->PostPreS_ProbXYNoL1VsProbQNoL1_highIas->Fill(probXYonTrackNoL1, 1 - probQonTrackNoL1, EventWeight_);
        }
        if (tof) {
          tuple->PostPreS_nDof->Fill(tof->nDof(), EventWeight_);
          tuple->PostPreS_MTOF->Fill(tof->inverseBeta(), EventWeight_);
          tuple->PostPreS_TOFError->Fill(tof->inverseBetaErr(), EventWeight_);
          tuple->PostPreS_TimeAtIP->Fill(tof->timeAtIpInOut(), EventWeight_);
        }
        if (track->quality(reco::TrackBase::highPurity)) {
          tuple->PostPreS_Qual->Fill(1., EventWeight_);
        } else {
          tuple->PostPreS_Qual->Fill(0., EventWeight_);
        }
        tuple->PostPreS_Chi2oNdof->Fill(track->chi2() / track->ndof(), EventWeight_);
        tuple->PostPreS_Chi2oNdofVsIas->Fill(track->chi2() / track->ndof(), globalIas_, EventWeight_);
        tuple->PostPreS_Pt->Fill(track->pt(), EventWeight_);
        tuple->PostPreS_Pt_lowPt->Fill(track->pt(), EventWeight_);
        tuple->PostPreS_PtVsIas->Fill(track->pt(), globalIas_, EventWeight_);
        tuple->PostPreS_P->Fill(track->p(), EventWeight_);
        tuple->PostPreS_NOMoNOH->Fill(numDeDxHits / (float)track->found(), EventWeight_);
        tuple->PostPreS_NOMoNOHvsPV->Fill(numGoodVerts, numDeDxHits / (float)track->found(), EventWeight_);
        tuple->PostPreS_Dz->Fill(dz, EventWeight_);
        tuple->PostPreS_DzVsIas->Fill(dz, globalIas_, EventWeight_);
        tuple->PostPreS_DzVsGenID->Fill(dz, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_Dxy->Fill(dxy, EventWeight_);
        tuple->PostPreS_DxyVsIas->Fill(dxy, globalIas_, EventWeight_);
        tuple->PostPreS_DxyVsGenID->Fill(dxy, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_PV->Fill(numGoodVerts, EventWeight_);
        tuple->PostPreS_PV_NoEventWeight->Fill(numGoodVerts);
        tuple->PostPreS_EoP->Fill(EoP, EventWeight_);
        tuple->PostPreS_EoPVsIas->Fill(EoP, globalIas_, EventWeight_);
        tuple->PostPreS_SumpTOverpT->Fill(IsoTK_SumEt / track->pt(), EventWeight_);
        tuple->PostPreS_SumpTOverpTVsIas->Fill(IsoTK_SumEt / track->pt(), globalIas_, EventWeight_);
        tuple->PostPreS_PtErrOverPt->Fill(track->ptError() / track->pt(), EventWeight_);
        tuple->PostPreS_PtErrOverPtVsIas->Fill(track->ptError() / track->pt(), globalIas_, EventWeight_);
        tuple->PostPreS_PtErrOverPt2VsIas->Fill(track->ptError() / (track->pt()*track->pt()), globalIas_, EventWeight_);
        tuple->PostPreS_PtErrOverPt2->Fill(track->ptError() / (track->pt()*track->pt()), EventWeight_);
        tuple->PostPreS_PtErrOverPtVsPtErrOverPt2->Fill(track->ptError() / track->pt(),track->ptError() / (track->pt()*track->pt()), EventWeight_);
        tuple->PostPreS_PtErrOverPtVsPt->Fill(track->ptError() / track->pt(), track->pt(), EventWeight_);
        tuple->PostPreS_TIsol->Fill(IsoTK_SumEt, EventWeight_);
        tuple->PostPreS_TIsolVsIas->Fill(IsoTK_SumEt, globalIas_,EventWeight_);
        tuple->PostPreS_Ih->Fill(globalIh_, EventWeight_);
        tuple->PostPreS_IhVsIas->Fill(globalIh_, globalIas_, EventWeight_);
        tuple->PostPreS_Ih_NoEventWeight->Fill(globalIh_);
        tuple->PostPreS_Ias_NoEventWeight->Fill(globalIas_);
        tuple->PostPreS_MassT->Fill(massT, EventWeight_);
        tuple->PostPreS_MassT_highMassT->Fill(massT, EventWeight_);
        tuple->PostPreS_MassTVsIas->Fill(massT, globalIas_, EventWeight_);
          // Add PFCadidate based isolation info to the tuple
          // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/
          // PhysicsTools/NanoAOD/plugins/IsoValueMapProducer.cc#L157
        tuple->PostPreS_MiniRelIsoAll->Fill(miniRelIsoAll, EventWeight_);
        tuple->PostPreS_MiniRelIsoAllVsIas->Fill(miniRelIsoAll, globalIas_, EventWeight_);
        tuple->PostPreS_MiniRelIsoChg->Fill(miniRelIsoChg, EventWeight_);
        tuple->PostPreS_MiniTkIso->Fill(track_genTrackMiniIsoSumPt, EventWeight_);
        tuple->PostPreS_MiniRelTkIso->Fill(track_genTrackMiniIsoSumPt, EventWeight_);

        tuple->PostPreS_MassErr->Fill(MassErr, EventWeight_);
        tuple->PostPreS_MassErrVsIas->Fill(MassErr, globalIas_, EventWeight_);

        tuple->PostPreS_EtaVsGenID->Fill(track->eta(), closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_ProbQVsGenID->Fill(1 - probQonTrack, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_ProbXYVsGenID->Fill(probXYonTrack, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_PtVsGenID->Fill(track->pt(), closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_EoPVsGenID->Fill(EoP, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_IhVsGenID->Fill(globalIh_, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_IasVsGenID->Fill(globalIas_, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_MassTVsGenID->Fill(massT, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_MiniIsoChgVsGenID->Fill(miniRelIsoChg, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_MiniIsoChgVsGenID->Fill(miniRelIsoAll, closestBackgroundPDGsIDs[0], EventWeight_);
        tuple->PostPreS_MassVsGenID->Fill(Mass, closestBackgroundPDGsIDs[0], EventWeight_);

        tuple->PostPreS_EtaVsMomGenID->Fill(track->eta(), closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_ProbQVsMomGenID->Fill(1 - probQonTrack, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_ProbXYVsMomGenID->Fill(probXYonTrack, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_PtVsMomGenID->Fill(track->pt(), closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_EoPVsMomGenID->Fill(EoP, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_IhVsMomGenID->Fill(globalIh_, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_IasVsMomGenID->Fill(globalIas_, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_MassTVsMomGenID->Fill(massT, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_MiniIsoChgVsMomGenID->Fill(miniRelIsoChg, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_MiniIsoAllVsMomGenID->Fill(miniRelIsoAll, closestBackgroundPDGsIDs[1], EventWeight_);
        tuple->PostPreS_MassVsMomGenID->Fill(Mass, closestBackgroundPDGsIDs[1], EventWeight_);

        tuple->PostPreS_EtaVsSiblingGenID->Fill(track->eta(), closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_ProbQVsSiblingGenID->Fill(1 - probQonTrack, closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_ProbXYVsSiblingGenID->Fill(probXYonTrack, closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_PtVsSiblingGenID->Fill(track->pt(), closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_EoPVsSiblingGenID->Fill(EoP, closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_IhVsSiblingGenID->Fill(globalIh_, closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_IasVsSiblingGenID->Fill(globalIas_, closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_MassTVsSiblingGenID->Fill(massT, closestBackgroundPDGsIDs[2], EventWeight_);
        tuple->PostPreS_MassVsSiblingGenID->Fill(Mass, closestBackgroundPDGsIDs[2], EventWeight_);

        tuple->PostPreS_EtaVsGenAngle->Fill(track->eta(), closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_ProbQVsGenAngle->Fill(1 - probQonTrack, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_ProbXYVsGenAngle->Fill(probXYonTrack, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_PtVsGenAngle->Fill(track->pt(), closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_EoPVsGenAngle->Fill(EoP, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_IhVsGenAngle->Fill(globalIh_, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_IasVsGenAngle->Fill(globalIas_, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_MassTVsGenAngle->Fill(massT, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_MiniIsoChgVsGenAngle->Fill(miniRelIsoChg, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_MiniIsoAllVsGenAngle->Fill(miniRelIsoAll, closestBackgroundPDGsIDs[3], EventWeight_);
        tuple->PostPreS_MassVsGenAngle->Fill(Mass, closestBackgroundPDGsIDs[3], EventWeight_);

        tuple->PostPreS_EtaVsGenMomAngle->Fill(track->eta(), closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_ProbQVsGenMomAngle->Fill(1 - probQonTrack, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_ProbXYVsGenMomAngle->Fill(probXYonTrack, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_PtVsGenMomAngle->Fill(track->pt(), closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_EoPVsGenMomAngle->Fill(EoP, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_IhVsGenMomAngle->Fill(globalIh_, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_IasVsGenMomAngle->Fill(globalIas_, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_MassTVsGenMomAngle->Fill(massT, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_MiniIsoChgVsGenMomAngle->Fill(miniRelIsoChg, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_MiniIsoAllVsGenMomAngle->Fill(miniRelIsoAll, closestBackgroundPDGsIDs[4], EventWeight_);
        tuple->PostPreS_MassVsGenMomAngle->Fill(Mass, closestBackgroundPDGsIDs[4], EventWeight_);

        tuple->PostPreS_GenPtVsRecoPt->Fill(closestBackgroundPDGsIDs[5], track->pt());

        tuple->PostPreS_EtaVsGenNumSibling->Fill(track->eta(), closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_ProbQVsGenNumSibling->Fill(1 - probQonTrack, closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_ProbXYVsGenNumSibling->Fill(probXYonTrack, closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_PtVsGenNumSibling->Fill(track->pt(), closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_EoPVsGenNumSibling->Fill(EoP, closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_IhVsGenNumSibling->Fill(globalIh_, closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_IasVsGenNumSibling->Fill(globalIas_, closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_MassTVsGenNumSibling->Fill(massT, closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_MiniIsoChgVsGenNumSibling->Fill(miniRelIsoChg, closestBackgroundPDGsIDs[6], EventWeight_);
        tuple->PostPreS_MiniIsoAllVsGenNumSibling->Fill(miniRelIsoAll, closestBackgroundPDGsIDs[6], EventWeight_);

        tuple->PostPreS_LastHitDXY->Fill(furthersHitDxy, EventWeight_);
        tuple->PostPreS_LastHitDXYVsEta->Fill(furthersHitDxy, track->eta(), EventWeight_);
        tuple->PostPreS_LastHitD3D->Fill(furthersHitDistance, EventWeight_);
        tuple->PostPreS_LastHitD3DVsEta->Fill(furthersHitDistance, track->eta(), EventWeight_);

        tuple->PostPreS_EoPVsPfType->Fill(EoP, 0., EventWeight_);
        tuple->PostPreS_MassVsPfType->Fill(Mass, 0., EventWeight_);
        if (pf_isPfTrack) {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 1., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 1., EventWeight_);
        } else {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 8., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 8., EventWeight_);
        }
        if (pf_isElectron) {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 2., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 2., EventWeight_);
        } else if (pf_isMuon) {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 3., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 3., EventWeight_);
        } else if (pf_isPhoton) {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 4., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 4., EventWeight_);
        } else if (pf_isChHadron) {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 5., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 5., EventWeight_);
        } else if (pf_isNeutHadron) {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 6., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 6., EventWeight_);
        } else if (pf_isUndefined) {
          tuple->PostPreS_EoPVsPfType->Fill(EoP, 7., EventWeight_);
          tuple->PostPreS_MassVsPfType->Fill(Mass, 7., EventWeight_);
        }

        tuple->PostPreS_Mass->Fill(Mass, EventWeight_);
        tuple->PostPreS_MassVsPt->Fill(Mass, track->pt(), EventWeight_);
        tuple->PostPreS_MassVsP->Fill(Mass, track->p(), EventWeight_);
        tuple->PostPreS_MassVsTNOHFraction->Fill(Mass, track->validFraction(), EventWeight_);
        tuple->PostPreS_MassVsTNOPH->Fill(Mass, nonL1PixHits, EventWeight_);
        tuple->PostPreS_MassVsTNOM->Fill(Mass, numDeDxHits, EventWeight_);
        tuple->PostPreS_MassVsProbQNoL1->Fill(Mass,1 - probQonTrackNoL1, EventWeight_);
        tuple->PostPreS_MassVsProbXYNoL1->Fill(Mass,probXYonTrackNoL1, EventWeight_);
        tuple->PostPreS_MassVsEoP->Fill(Mass, EoP, EventWeight_);
        tuple->PostPreS_MassVsSumpTOverpT->Fill(Mass, IsoTK_SumEt / track->pt(), EventWeight_);
        tuple->PostPreS_MassVsPtErrOverPt->Fill(Mass, track->ptError() / track->pt(), EventWeight_);
        tuple->PostPreS_MassVsTIsol->Fill(Mass, IsoTK_SumEt,EventWeight_);
        tuple->PostPreS_MassVsIh->Fill(Mass, globalIh_, EventWeight_);
        tuple->PostPreS_MassVsMassT->Fill(Mass, massT, EventWeight_);
        tuple->PostPreS_MassVsMiniRelIsoAll->Fill(Mass, miniRelIsoAll, EventWeight_);
        tuple->PostPreS_MassVsMassErr->Fill(Mass, MassErr, EventWeight_);
        tuple->PostPreS_dRMinPfJet->Fill(dRMinPfJetWithCuts, EventWeight_);
        tuple->PostPreS_closestPfJetMuonFraction->Fill(closestPfJetMuonFraction, EventWeight_);
        tuple->PostPreS_closestPfJetElectronFraction->Fill(closestPfJetElectronFraction, EventWeight_);
        tuple->PostPreS_closestPfJetPhotonFraction->Fill(closestPfJetPhotonFraction, EventWeight_);
        tuple->PostPreS_dRMinPfJetVsIas->Fill(dRMinPfJetWithCuts, globalIas_, EventWeight_);
        tuple->PostPreS_closestPfJetMuonFractionVsIas->Fill(closestPfJetMuonFraction, globalIas_, EventWeight_);
        tuple->PostPreS_closestPfJetElectronFractionVsIas->Fill(closestPfJetElectronFraction, globalIas_, EventWeight_);
        tuple->PostPreS_closestPfJetPhotonFractionVsIas->Fill(closestPfJetPhotonFraction, globalIas_, EventWeight_);
        tuple->PostPreS_dRMinCaloJet->Fill(dRMinCaloJet, EventWeight_);
        tuple->PostPreS_dPhiMinPfMet->Fill(dPhiMinPfMet, EventWeight_);
        tuple->PostPreS_CaloJetsNum->Fill(caloJetsNum, EventWeight_);
        tuple->PostPreS_dRMinCaloJetVsIas->Fill(dRMinCaloJet, globalIas_, EventWeight_);
        tuple->PostPreS_dPhiMinPfMetVsIas->Fill(dPhiMinPfMet, globalIas_, EventWeight_);
        tuple->PostPreS_RecoPfMet->Fill(RecoPFMET, EventWeight_);
        tuple->PostPreS_RecoPfMetPhi->Fill(RecoPFMET_phi, EventWeight_);
        tuple->PostPreS_RecoPfJetsNum->Fill(pfJetsNum, EventWeight_);
        tuple->PostPreS_RecoPfHT->Fill(pfJetHT, EventWeight_);
        if (genBeta >= 0) {
          tuple->PostPreS_GenBeta->Fill(genBeta, EventWeight_);
        }

        // stability plots
        tuple->Stab_Ih_NoL1_VsRun->Fill(CurrentRun_, globalIh_ , EventWeight_);
        tuple->Stab_Ih_pixNoL1_VsRun->Fill(CurrentRun_, (dedxIh_PixelOnlyh_noL1 ? dedxIh_PixelOnlyh_noL1->dEdx() : -1) , EventWeight_);
        tuple->Stab_Ih_strip_VsRun->Fill(CurrentRun_, (dedxIh_StripOnly ? dedxIh_StripOnly->dEdx() : -1) , EventWeight_);
        tuple->Stab_Gi_strip_VsRun->Fill(CurrentRun_, globalIas_ , EventWeight_);
        tuple->Stab_Gi_NoL1_VsRun->Fill(CurrentRun_, (dedxIas_noL1 ? dedxIas_noL1->dEdx() : -1), EventWeight_);
        tuple->Stab_Fi_pixNoL1_VsRun->Fill(CurrentRun_, 1 - probQonTrackNoL1, EventWeight_);
        if (tof) {
          tuple->Stab_invB_VsRun->Fill(CurrentRun_, tof->inverseBeta(), EventWeight_);
          tuple->Stab_invB_DT_VsRun->Fill(CurrentRun_, dttof->inverseBeta(), EventWeight_);
          tuple->Stab_invB_CSC_VsRun->Fill(CurrentRun_, csctof->inverseBeta(), EventWeight_);
        }

      }

      if ((((globalIas_ > 0.3 || Mass > 1000) && !isSignal) || (debug_ > 7))  && trigInfo_ > 0) {
        if (globalIas_ > 0.3)    { LogPrint(MOD) << "\n        >> After passing preselection, the globalIas_ > 0.3";}
        if (Mass > 1000) { LogPrint(MOD) << "\n        >> After passing preselection, the Mass > 1000";}
        LogPrint(MOD) << "        >> LS: " << iEvent.luminosityBlock() << " Event number: " << iEvent.id().event();
        LogPrint(MOD) << "        >> -----------------------------------------------";
        LogPrint(MOD) << "        >> Trigger passed!" ;
        LogPrint(MOD) << "        >> track->eta()  " <<   track->eta() ;
        LogPrint(MOD) << "        >> track->pt()  " <<   track->pt() ;
        LogPrint(MOD) << "        >> track->found()  " <<   track->found() ;
        LogPrint(MOD) << "        >> track->hitPattern().numberOfValidPixelHits()  " <<   track->hitPattern().numberOfValidPixelHits() ;
        LogPrint(MOD) << "        >> track->validFraction()  " <<   track->validFraction() ;
        LogPrint(MOD) << "        >> numDeDxHits  " <<   numDeDxHits ;
        LogPrint(MOD) << "        >> track->chi2() / track->ndof()   " <<   track->chi2() / track->ndof() ;
        LogPrint(MOD) << "        >> EoP   " <<   EoP << "     --> | PF E = " << pf_energy <<  " | Cone based (0.3) E = " << hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy() << " | p = " << track->p() << " | " ;
        LogPrint(MOD) << "        >> dz  " <<   dz ;
        LogPrint(MOD) << "        >> dxy  " <<   dxy ;
        LogPrint(MOD) << "        >> track->ptError() / track->pt()  " <<   track->ptError() / track->pt() ;
        LogPrint(MOD) << "        >> pTerr_over_pT_etaBin(track->pt(), track->eta())  " <<   pTerr_over_pT_etaBin(track->pt(), track->eta()) ;
        LogPrint(MOD) << "        >> IsoTK_SumEt   " <<   IsoTK_SumEt  ;
        LogPrint(MOD) << "        >> miniRelIsoAll   " <<   miniRelIsoAll  ;
        LogPrint(MOD) << "        >> globalIh_  " <<   globalIh_ ;
        LogPrint(MOD) << "        >> globalIas_  " << globalIas_ ;
        LogPrint(MOD) << "        >> probQonTrack   " <<   probQonTrack << " | probQonTrackNoL1 " << probQonTrackNoL1 ;
        LogPrint(MOD) << "        >> probXYonTrack  " <<   probXYonTrack << " | probXYonTrackNoL1 " << probXYonTrackNoL1 ;
        LogPrint(MOD) << "        >> dRMinCaloJet   " <<   dRMinCaloJet ;
        LogPrint(MOD) << "        >> dRMinPfJet   " <<   dRMinPfJetWithCuts ;
        LogPrint(MOD) << "        >> closestPfJetMuonFraction   " <<   closestPfJetMuonFraction ;
        LogPrint(MOD) << "        >> closestPfJetElectronFraction   " <<   closestPfJetElectronFraction ;
        LogPrint(MOD) << "        >> closestPfJetPhotonFraction   " <<   closestPfJetPhotonFraction ;
      }
      if (globalIh_ > maxIhSoFar) {
        maxIhSoFar = globalIh_;
        bestCandidateIndex = candidate_count;
        bestCandidateIas = globalIas_;
        bestCandidateProbQNoL1 = probQonTrackNoL1;
        bestCandidateDrMinHltMuon = dr_min_hlt_muon;
        bestCandidateGenBeta = genBeta;
      }
    }

    // Let's do some printouts after preselections for gen particles
    if (passPre && trigInfo_ > 0 && closestGenIndex > 0) {
      if (!isData) {
        if (debug_> 3) LogPrint(MOD) << "  >> MC, set gen IDs, mother IDs, sibling IDs";
        closestBackgroundPDGsIDs[0] = (float)abs(genPdgId);
        float dRMinGenAndSibling = 9999.0;
        float dRMinGenAndMom = 9999.0;
        float numSiblingsF = 9999.0;
        bool motherFound = false;
        float dRMinGenAndAunt = 9999.0;
//        float dRMinGenAndGrandAunt = 9999.0;
        reco::GenParticle& genCandidateUnderStudy = genColl[closestGenIndex];

        if (genCandidateUnderStudy.numberOfMothers() == 0) {
          LogPrint(MOD) << "There are zero mothers, track ID" << abs(genCandidateUnderStudy.pdgId()) <<
          " Eta: " << genEta << " Phi: " << genPhi ;
        }
      if (debug_ > 3) {
        // HSCP muon
        std::cout << "           | Relation | ID | $p_{T}$ | $v_{x}$ |  $v_{y}$ |  $v_{z}$ |  $R_{xy}$ | " << endl;
        std::cout << "           |--- | ---| " << std::endl;
        std::cout << "           | Me | " << genCandidateUnderStudy.pdgId() << " | " << genCandidateUnderStudy.pt()
        << " | " << genCandidateUnderStudy.vx() << " | " << genCandidateUnderStudy.vy() << " | " << genCandidateUnderStudy.vz()
        << " | " << sqrt(genCandidateUnderStudy.vx()*genCandidateUnderStudy.vx()+genCandidateUnderStudy.vy()*genCandidateUnderStudy.vy()) << " | "  << endl;
        // photon
        if (genCandidateUnderStudy.mother(0)->numberOfDaughters() > 1) {
          std::cout << "           | Sibling (1) | " << genCandidateUnderStudy.mother(0)->daughter(1)->pdgId() << " | " << genCandidateUnderStudy.mother(0)->daughter(1)->pt()
          << " | " << genCandidateUnderStudy.mother(0)->daughter(1)->vx() << " | " << genCandidateUnderStudy.mother(0)->daughter(1)->vy() << " | " << genCandidateUnderStudy.mother(0)->daughter(1)->vz()
          << " | " << sqrt(genCandidateUnderStudy.mother(0)->daughter(1)->vx()*genCandidateUnderStudy.mother(0)->daughter(1)->vx() + genCandidateUnderStudy.mother(0)->daughter(1)->vy()*genCandidateUnderStudy.mother(0)->daughter(1)->vy()) << " | " << endl;
        }
        if ( genCandidateUnderStudy.numberOfMothers() > 0) {
          // mother muon
          std::cout << "           | Mom | " << genCandidateUnderStudy.mother(0)->pdgId() << " | " << genCandidateUnderStudy.mother(0)->pt()
          << " | " << genCandidateUnderStudy.mother()->vx() << " | " << genCandidateUnderStudy.mother()->vy() << " | " << genCandidateUnderStudy.mother()->vz()
          << " | " << sqrt(genCandidateUnderStudy.mother()->vx()*genCandidateUnderStudy.mother()->vx()+genCandidateUnderStudy.mother()->vy()*genCandidateUnderStudy.mother()->vy()) << " | " << endl;
        }
        if (genCandidateUnderStudy.mother(0)->mother(0)->numberOfDaughters() > 1) {
          // neutrino
          std::cout << "           | Aunt (1) | " << genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->pdgId() << " | " << genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->pt()
          << " | " << genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->vx() << " | " << genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->vy() << " | " << genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->vz()
          << " | " << sqrt(genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->vx() * genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->vx() + genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->vy() * genCandidateUnderStudy.mother(0)->mother(0)->daughter(1)->vy()) << " | " << endl;
        }
        // kaon / D+
        if (genCandidateUnderStudy.mother()->numberOfMothers() > 0) {
          std::cout << "           | Gramma |  " << genCandidateUnderStudy.mother()->mother()->pdgId()   << " | " << genCandidateUnderStudy.mother()->mother()->pt()
          << " | " << genCandidateUnderStudy.mother()->mother()->vx() << " | " << genCandidateUnderStudy.mother()->mother()->vy() << " | " << genCandidateUnderStudy.mother()->mother()->vz()
          << " | " << sqrt(genCandidateUnderStudy.mother()->mother()->vx()*genCandidateUnderStudy.mother()->mother()->vx() + genCandidateUnderStudy.mother()->mother()->vy()*genCandidateUnderStudy.mother()->mother()->vy()) << " | "  << endl;
        }
        if (genCandidateUnderStudy.mother()->mother()->numberOfMothers() > 0) {
          std::cout << "           | GrandGramma |  " << genCandidateUnderStudy.mother()->mother()->mother()->pdgId()   << " | " << genCandidateUnderStudy.mother()->mother()->mother()->pt()
        << " | " << genCandidateUnderStudy.mother()->mother()->mother()->vx() << " | " << genCandidateUnderStudy.mother()->mother()->mother()->vy() << " | " << genCandidateUnderStudy.mother()->mother()->mother()->vz()
        << " | " << sqrt(genCandidateUnderStudy.mother()->mother()->mother()->vx()*genCandidateUnderStudy.mother()->mother()->mother()->vx() + genCandidateUnderStudy.mother()->mother()->mother()->vy()*genCandidateUnderStudy.mother()->mother()->mother()->vy()) << " | "  << endl;
        }
        if (genCandidateUnderStudy.mother()->mother()->mother()->numberOfMothers() > 0) {
          std::cout << "           | GreatGrandGramma |  " << genCandidateUnderStudy.mother()->mother()->mother()->mother()->pdgId()   << " | " << genCandidateUnderStudy.mother()->mother()->mother()->mother()->pt()
          << " | " << genCandidateUnderStudy.mother()->mother()->mother()->mother()->vx() << " | " << genCandidateUnderStudy.mother()->mother()->mother()->mother()->vy() << " | " << genCandidateUnderStudy.mother()->mother()->mother()->mother()->vz()
          << " | " << sqrt(genCandidateUnderStudy.mother()->mother()->mother()->mother()->vx()*genCandidateUnderStudy.mother()->mother()->mother()->mother()->vx() + genCandidateUnderStudy.mother()->mother()->mother()->mother()->vy()*genCandidateUnderStudy.mother()->mother()->mother()->mother()->vy()) << " | "  << endl;
        }
      }
          // Loop through all the mothers of the gen particle
        for (unsigned int numMomIndx = 0; numMomIndx < genCandidateUnderStudy.numberOfMothers(); numMomIndx++) {
          if (abs(genCandidateUnderStudy.mother(numMomIndx)->pdgId())  != abs(genCandidateUnderStudy.pdgId())) {
            closestBackgroundPDGsIDs[1] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->pdgId());
            closestBackgroundPDGsIDs[7] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->pt());
            unsigned int numSiblings = genCandidateUnderStudy.mother(numMomIndx)->numberOfDaughters() -1;
            numSiblingsF  = float(numSiblings);
            if (globalIas_ > 0.3 && trigInfo_ > 0 && debug_ > 4) LogPrint(MOD) << "      >> Number of siblings: " << numSiblings << ". Me and my syblings: ";
            for (unsigned int daughterIndx = 0; daughterIndx < numSiblings+1; daughterIndx++) {
              if (globalIas_ > 0.3 && debug_ > 4 && trigInfo_ > 0) std::cout << "      >> " << genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->pdgId() ;
              float siblingEta = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->eta();
              float siblingPhi = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->phi();
              float siblingPt  = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->pt();
              float siblingDr = deltaR(genEta, genPhi, siblingEta, siblingPhi);
              if (globalIas_ > 0.3 && debug_ > 4) std::cout << " (dR = " << siblingDr << ", pt = " << siblingPt <<  ") , ";
              if( (siblingDr != 0.0) && (siblingDr < dRMinGenAndSibling)) {
                dRMinGenAndSibling = siblingDr;
                closestBackgroundPDGsIDs[2] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->pdgId());
              }
            }
            float momEta = genCandidateUnderStudy.mother(numMomIndx)->eta();
            float momPhi = genCandidateUnderStudy.mother(numMomIndx)->phi();
            dRMinGenAndMom = deltaR(genEta, genPhi, momEta, momPhi);
            motherFound = true;
            break;
          }
        }
        if (globalIas_ > 0.3 && trigInfo_ > 0 && debug_ > 4) std::cout << std::endl;
        
        // If the loop on the mothers didnt find the mother (e.g. all moms had the same ID), let's look at the grandmas
        if (!motherFound) {
          if (globalIas_ > 0.3 && trigInfo_ > 0 && debug_ > 4) LogPrint(MOD) << "      >> All moms had the same ID as the candidate, let's look at the grammas";
          for (unsigned int numMomIndx = 0; numMomIndx < genCandidateUnderStudy.numberOfMothers(); numMomIndx++) {
            for (unsigned int numGramMomIndx = 0; numGramMomIndx < genCandidateUnderStudy.mother(numMomIndx)->numberOfMothers(); numGramMomIndx++) {
              if (abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->pdgId())  != abs(genCandidateUnderStudy.pdgId())) {
                closestBackgroundPDGsIDs[1] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->pdgId());
                closestBackgroundPDGsIDs[7] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->pt());
                unsigned int numSiblings = genCandidateUnderStudy.mother(numMomIndx)->numberOfDaughters() -1;
                if (globalIas_ > 0.3 && trigInfo_ > 0 && debug_ > 4) LogPrint(MOD) << "      >> Number of siblings: " << numSiblings << ". Me and my syblings: ";
                for (unsigned int daughterIndx = 0; daughterIndx < numSiblings+1; daughterIndx++) {
                  if (globalIas_ > 0.3 && debug_ > 4 && trigInfo_ > 0) std::cout << "      >> " << genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->pdgId() ;
                  float siblingEta = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->eta();
                  float siblingPhi = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->phi();
                  float siblingPt  = genCandidateUnderStudy.mother(numMomIndx)->daughter(daughterIndx)->pt();
                  float siblingDr = deltaR(genEta, genPhi, siblingEta, siblingPhi);
                  if (globalIas_ > 0.3) std::cout << " (dR = " << siblingDr << ", pt = " << siblingPt <<  ") , ";
                }
                if (globalIas_ > 0.3) std::cout << std::endl;
                unsigned int numAunts = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->numberOfDaughters() -1;
                numSiblingsF  = float(numSiblings);
                if (globalIas_ > 0.3 && trigInfo_ > 0 && debug_ > 4) LogPrint(MOD) << "      >> Number of aunts: " << numAunts << ". Mom with same ID as the candidate and her syblings: ";
                for (unsigned int daughterIndx = 0; daughterIndx < numAunts+1; daughterIndx++) {
                  if (globalIas_ > 0.3 && debug_ > 4 && trigInfo_ > 0) std::cout << "      >> "  << genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->pdgId() ;
                  float auntEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->eta();
                  float auntPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->phi();
                  float auntPt = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->pt();
                  float auntDr = deltaR(genEta, genPhi, auntEta, auntPhi);
                  if (globalIas_ > 0.3 && debug_ > 4) std::cout << " (dR = " << auntDr << ", pt =  " << auntPt <<  ") , ";
                  if( (auntDr != 0.0) && (auntDr < dRMinGenAndAunt)) {
                    dRMinGenAndAunt = auntDr;
                    closestBackgroundPDGsIDs[2] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->daughter(daughterIndx)->pdgId());
                  }
                }
                float momEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->eta();
                float momPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->phi();
                dRMinGenAndMom = deltaR(genEta, genPhi, momEta, momPhi);
                motherFound = true;
                break;
              }
            }
            if (motherFound) break;
          }
          if (globalIas_ > 0.3 && trigInfo_ > 0 && debug_ > 4) std::cout << std::endl;
        }

        // If none of the mothers' mother's is the real mother (e.g. all moms'moms had the same ID as the candidate), let's look at the grand-grandmas
        if (!motherFound) {
          if (debug_ > 4 && trigInfo_ > 0 && debug_ > 4) LogPrint(MOD) << "      >> All moms' moms had the same ID as the candidate, let's look at the grand-grammas";
          for (unsigned int numMomIndx = 0; numMomIndx < genCandidateUnderStudy.numberOfMothers(); numMomIndx++) {
            for (unsigned int numGramMomIndx = 0; numGramMomIndx < genCandidateUnderStudy.mother(numMomIndx)->numberOfMothers(); numGramMomIndx++) {
              for (unsigned int numGrandGramMomIndx = 0; numGrandGramMomIndx < genCandidateUnderStudy.mother(numGramMomIndx)->mother(numGramMomIndx)->numberOfMothers(); numGrandGramMomIndx++) {
                if (abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->pdgId())  != abs(genCandidateUnderStudy.pdgId())) {
                  closestBackgroundPDGsIDs[1] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->pdgId());
                  closestBackgroundPDGsIDs[7] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->pt());
                  unsigned int numSiblings = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->numberOfDaughters() -1;
                  numSiblingsF  = float(numSiblings);
                  if (globalIas_ > 0.3 && trigInfo_ > 0 && debug_ > 4) LogPrint(MOD) << "      >> Number of great-aunts: " << numSiblings << ". Gramma with same ID the candidate and her syblings: ";
                  for (unsigned int daughterIndx = 0; daughterIndx < numSiblings+1; daughterIndx++) {
                    if (globalIas_ > 0.3 && debug_ > 4 && trigInfo_ > 0) std::cout << "      >> "  << genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->pdgId() ;
                    float siblingEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->eta();
                    float siblingPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->phi();
                    float siblingPt = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->pt();
                    float siblingDr = deltaR(genEta, genPhi, siblingEta, siblingPhi);
                    if (globalIas_ > 0.3 && debug_ > 4) std::cout << " (dR = " << siblingDr << ", pt =  " << siblingPt <<  ") , ";
                    if( (siblingDr != 0.0) && (siblingDr < dRMinGenAndSibling)) {
                      dRMinGenAndSibling = siblingDr;
                      closestBackgroundPDGsIDs[2] = (float)abs(genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->daughter(daughterIndx)->pdgId());
                    }
                  }
                  float momEta = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->eta();
                  float momPhi = genCandidateUnderStudy.mother(numMomIndx)->mother(numGramMomIndx)->mother(numGrandGramMomIndx)->phi();
                  dRMinGenAndMom = deltaR(genEta, genPhi, momEta, momPhi);
                  motherFound = true;
                  break;
                }
              }
              if (motherFound) break;
            }
            if (motherFound) break;
          }
        }
        if (!motherFound) {
          if (debug_ > 4 && trigInfo_ > 0) LogPrint(MOD) << "      >> All moms' mom's moms had the same ID as the candidate -- is this realy possible at this point???";
        }
          // I'm sure this could be done better, if you agree and feel like it, please fix it
          // issue with a while loop and a recursive I faced is tha that mom doesnt have the same type as the genParticle

        closestBackgroundPDGsIDs[3] = dRMinGenAndSibling;
        closestBackgroundPDGsIDs[4] = dRMinGenAndMom;
        closestBackgroundPDGsIDs[5] = fabs(genColl[closestGenIndex].pt());
        closestBackgroundPDGsIDs[6] = numSiblingsF;

        if ((((debug_> 2) || (!isSignal && globalIas_ > 0.3))) && (trigInfo_ > 0)) {
          LogPrint(MOD) << "      >> Track's gen ID: " << closestBackgroundPDGsIDs[0];
          LogPrint(MOD) << "      >> Track's gen pt: " << closestBackgroundPDGsIDs[5];

          LogPrint(MOD) << "      >> Track's mom/gramma/grand-gramma ID: " << closestBackgroundPDGsIDs[1];
          LogPrint(MOD) << "      >> Track's mom/gramma/grand-gramma pt: " << closestBackgroundPDGsIDs[7];
          LogPrint(MOD) << "      >> Track's gen angle wrt to mom/gramma/grand-gramma: " << closestBackgroundPDGsIDs[4];

          LogPrint(MOD) << "      >> Track's num siblings: " << closestBackgroundPDGsIDs[6];
          LogPrint(MOD) << "      >> Track's closest sibling gen ID: " << closestBackgroundPDGsIDs[2];
          LogPrint(MOD) << "      >> Track's closest sibling gen angle: " << closestBackgroundPDGsIDs[3];


          float candidateEta = genColl[closestGenIndex].eta();
          float candidatePhi = genColl[closestGenIndex].phi();
          cout << "          | ID  | distance | pt | status | " << endl;
          std::cout << "          |--- | ---| " << std::endl;
          for (unsigned int g = 0; g < genColl.size(); g++) {
            float status = genColl[g].status();
            float pt = genColl[g].pt();
            float ID = genColl[g].pdgId();

            float muonDr = deltaR(genColl[g].eta(), genColl[g].phi(), candidateEta, candidatePhi);
            if (muonDr > 0.1) continue;

            cout << "          | " << ID;
            std::cout << " | " << muonDr << " | " << pt <<  " | " <<  status << " |  " << endl;
          }
        }
      }
    }
    
    if (!isData && passPre) {
      if ( hscp.type() == susybsm::HSCParticleType::globalMuon) {
        tuple->PostPreS_RecoHSCParticleType->Fill(0.);
      } else if ( hscp.type() == susybsm::HSCParticleType::trackerMuon) {
        tuple->PostPreS_RecoHSCParticleType->Fill(1.);
      } else if ( hscp.type() == susybsm::HSCParticleType::matchedStandAloneMuon) {
        tuple->PostPreS_RecoHSCParticleType->Fill(2.);
      } else if ( hscp.type() == susybsm::HSCParticleType::standAloneMuon) {
        tuple->PostPreS_RecoHSCParticleType->Fill(3.);
      } else if ( hscp.type() == susybsm::HSCParticleType::innerTrack) {
        tuple->PostPreS_RecoHSCParticleType->Fill(4.);
      } else if ( hscp.type() == susybsm::HSCParticleType::unknown) {
        tuple->PostPreS_RecoHSCParticleType->Fill(5.);
      }
    }

    if (!isData && passPre) {
      bool hasStatus91Around = false;
      unsigned int usignedIntclosestGenIndex = 0;
      if (closestGenIndex>0) usignedIntclosestGenIndex = closestGenIndex;

      for (unsigned int g = 0; g < genColl.size(); g++) {
        // Exclude the canidate when looking at its envirment
        if (g == usignedIntclosestGenIndex) continue;
        // Look only at the R=0.1 enviroment of the candidate
        if (deltaR(genColl[g].eta(),genColl[g].phi(),track->eta(),track->phi()) > 0.001) continue;

        if (genColl[g].status() == 91) hasStatus91Around = true;
        // Consider non-status 1 particles
        if (genColl[g].status() != 1) continue;
        tuple->PostPreS_ProbQVsGenEnviromentID->Fill(probQonTrack, abs(genColl[g].pdgId()), EventWeight_);
        tuple->PostPreS_IasVsGenEnviromentID->Fill(globalIas_, abs(genColl[g].pdgId()), EventWeight_);
      }
      if (hasStatus91Around && passPre && trigInfo_ > 0) {
        tuple->PostPreS_IasForStatus91->Fill(globalIas_, EventWeight_);
      } else if (passPre && trigInfo_ > 0) {
        tuple->PostPreS_IasForStatusNot91->Fill(globalIas_, EventWeight_);
      }
    }
    
    if (passPre) {
    // Loop through the deDx hits after the preselection
    bool headerStripsPrintedAlready = false;
    bool headerPixPrintedAlready = false;
    for (unsigned int i = 0; i < dedxHits->size(); i++) {
      DetId detid(dedxHits->detId(i));

      // The pixel part
      if (detid.subdetId() < 3 && passPre && trigInfo_ > 0) {
        // Taking the pixel cluster
        auto const* pixelCluster =  dedxHits->pixelCluster(i);
        // Get the local angles (axproximate from global)
        const GeomDetUnit& geomDet = *tkGeometry->idToDetUnit(detid);
        LocalVector lv = geomDet.toLocal(GlobalVector(track->px(), track->py(), track->pz()));
        // Let's redo CPE too
        auto reCPE = std::get<2>(pixelCPE->getParameters(*pixelCluster, geomDet, LocalTrajectoryParameters(dedxHits->pos(i), lv, track->charge())));
          // extract probQ and probXY from this
        float probQ = SiPixelRecHitQuality::thePacking.probabilityQ(reCPE);
        float probXY = SiPixelRecHitQuality::thePacking.probabilityXY(reCPE);
        if (probQ <= 0.0 || probQ >= 1.f) probQ = 1.f;
        if (probXY <= 0.0 || probXY >= 1.f) probXY = 0.f;

        bool isOnEdge = SiPixelRecHitQuality::thePacking.isOnEdge(reCPE);
        bool hasBadPixels = SiPixelRecHitQuality::thePacking.hasBadPixels(reCPE);
        bool spansTwoROCs = SiPixelRecHitQuality::thePacking.spansTwoROCs(reCPE);
        auto momentum = track->p();
        auto cotAlpha = lv.x()/lv.z();
        auto cotBeta = lv.y()/lv.z();
        auto clustSize = pixelCluster->size();
        auto clustSizeX = pixelCluster->sizeX();
        auto clustSizeY = pixelCluster->sizeY();
        auto pixelNormCharge = cm2umUnit * dedxHits->charge(i) / dedxHits->pathlength(i);

        float tmp1 = geomDet.surface().toGlobal(Local3DPoint(0.,0.,0.)).perp();
        float tmp2 = geomDet.surface().toGlobal(Local3DPoint(0.,0.,1.)).perp();
        int isFlippedModule = 0;
        if (tmp2 < tmp1) isFlippedModule = 1;

        bool specInCPE = false;
        (isOnEdge || hasBadPixels || spansTwoROCs) ? specInCPE = true : specInCPE = false;


        if ( detid.subdetId() == PixelSubdetector::PixelBarrel) {

          auto pixLayerIndex = abs(int(tTopo->pxbLayer(detid)));
          tuple->PostPreS_CluProbQVsPixelLayer->Fill(probQ, pixLayerIndex, EventWeight_);
          tuple->PostPreS_CluProbXYVsPixelLayer->Fill(probXY, pixLayerIndex, EventWeight_);
          tuple->PostPreS_CluSizeVsPixelLayer->Fill(clustSize, pixLayerIndex, EventWeight_);
          tuple->PostPreS_CluSizeXVsPixelLayer->Fill(clustSizeX, pixLayerIndex, EventWeight_);
          tuple->PostPreS_CluSizeYVsPixelLayer->Fill(clustSizeY, pixLayerIndex, EventWeight_);
          if (isOnEdge) {
            tuple->PostPreS_CluSpecInCPEVsPixelLayer->Fill(0., pixLayerIndex, EventWeight_);
          } else if (hasBadPixels) {
            tuple->PostPreS_CluSpecInCPEVsPixelLayer->Fill(1., pixLayerIndex, EventWeight_);
          } else if (spansTwoROCs) {
            tuple->PostPreS_CluSpecInCPEVsPixelLayer->Fill(2., pixLayerIndex, EventWeight_);
          }
          tuple->PostPreS_CluSpecInCPEVsPixelLayer->Fill(3., pixLayerIndex, EventWeight_);
          if (globalIas_ > 0.3) {
            tuple->PostPreS_CluProbQVsPixelLayer_highIas->Fill(probQ, pixLayerIndex, EventWeight_);
            tuple->PostPreS_CluProbXYVsPixelLayer_highIas->Fill(probXY, pixLayerIndex, EventWeight_);
            tuple->PostPreS_CluSizeVsPixelLayer_highIas->Fill(clustSize, pixLayerIndex, EventWeight_);
            tuple->PostPreS_CluSizeXVsPixelLayer_highIas->Fill(clustSizeX, pixLayerIndex, EventWeight_);
            tuple->PostPreS_CluSizeYVsPixelLayer_highIas->Fill(clustSizeY, pixLayerIndex, EventWeight_);
            if (isOnEdge) {
              tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->Fill(0., pixLayerIndex, EventWeight_);
            } else if (hasBadPixels) {
              tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->Fill(1., pixLayerIndex, EventWeight_);
            } else if (spansTwoROCs) {
              tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->Fill(2., pixLayerIndex, EventWeight_);
            }
            tuple->PostPreS_CluSpecInCPEVsPixelLayer_highIas->Fill(3., pixLayerIndex, EventWeight_);
          }

          if (probXY < globalMinTrackProbXYCut_ && !specInCPE) {
            tuple->PostPreS_CluCotBetaVsPixelLayer_lowProbXY->Fill(cotBeta, pixLayerIndex, EventWeight_);
            tuple->PostPreS_CluCotAlphaVsPixelLayer_lowProbXY->Fill(cotAlpha, pixLayerIndex, EventWeight_);
          } else if (probXY > globalMinTrackProbXYCut_ && !specInCPE) {
            tuple->PostPreS_CluCotBetaVsPixelLayer->Fill(cotBeta, pixLayerIndex, EventWeight_);
            tuple->PostPreS_CluCotAlphaVsPixelLayer->Fill(cotAlpha, pixLayerIndex, EventWeight_);
          }
          // 0.31623 [Bichsel's smallest entry]  && genGammaBeta > 0.31623
          // (globalIas_ > 0.02 && globalIas_ < 0.03)
          if ((!isSignal && globalIas_ > 0.3) || debug_ > 4) {
            if (!headerPixPrintedAlready) {
              std::cout << std::endl << "        | $G_{i}^{strips}$ | Layer | gammaBeta | flipped | cotAlpha | cotBeta | momentum | sizeX | sizeY";
              std::cout << " | Norm. Charge | edge | bad | double | cProbXY | cProbQ | " << std::endl;
              std::cout << "        |--- | ---| " << std::endl;

              headerPixPrintedAlready = true;
            }

            std::cout  << "        | " <<  globalIas_ << " | L" << tTopo->pxbLayer(detid) << " | " << genGammaBeta << " | " << isFlippedModule << " | ";
            std::cout << cotAlpha << " | " << cotBeta << " | " << momentum<< " | " << clustSizeX << " | " << clustSizeY << " | ";
            std::cout << pixelNormCharge << " e/um | " << isOnEdge  << " | " << hasBadPixels  << " | " << spansTwoROCs << " | " << probXY << " | " << probQ <<  " | " << std::endl;
          } else if (!isData && genGammaBeta <= 0.31623 && trigInfo_ > 0 && debug_ > 4)  {
            LogPrint(MOD) << "BetaGamma is too low for Bichsel";
          }
        }

      // the strip part
      } else if (passPre && trigInfo_ > 0 && detid.subdetId() >= 3) {
          // Taking the strips cluster
        auto const* stripsCluster = dedxHits->stripCluster(i);
        std::vector<int> amplitudes = convert(stripsCluster->amplitudes());
        std::vector<int> amplitudesPrim = CrossTalkInv(amplitudes,0.10,0.04,true);
        unsigned int clusterCleaned = (clusterCleaning(amplitudesPrim, 1)) ? 0 : 1;


        float stripNormCharge = cm2umUnit * dedxHits->charge(i) * 265 / dedxHits->pathlength(i);
        float stripSize = stripsCluster->amplitudes().size();
        
        unsigned int stripLayerIndex = 0;
        if (detid.subdetId() == StripSubdetector::TIB) stripLayerIndex = abs(int(tTopo->tibLayer(detid)));
        if (detid.subdetId() == StripSubdetector::TOB) stripLayerIndex = abs(int(tTopo->tobLayer(detid))) + 4;
        if (detid.subdetId() == StripSubdetector::TID) stripLayerIndex = abs(int(tTopo->tidWheel(detid))) + 10;
        if (detid.subdetId() == StripSubdetector::TEC) stripLayerIndex = abs(int(tTopo->tecWheel(detid))) + 13;
        
        if (!isData && genGammaBeta > 0.31623 && genGammaBeta < 0.6 ) {
          tuple->PostPreS_CluNormChargeVsStripLayer_lowBetaGamma->Fill(stripNormCharge, stripLayerIndex, EventWeight_);
        } else if (!isData && genGammaBeta > 0.6 ) {
          tuple->PostPreS_CluNormChargeVsStripLayer_higherBetaGamma->Fill(stripNormCharge, stripLayerIndex, EventWeight_);
        }
        // || (globalIas_ > 0.025 && globalIas_ < 0.03)
        if  (globalIas_ > 0.3 && debug_ > 4) {
          unsigned int isGlued = (tTopo->glued(detid) > 0) ? 1 : 0;
          if (!headerStripsPrintedAlready) {
            std::cout << std::endl <<  "        | $G_{i}^{strips}$  | Layer | gammaBeta | eta | Norm. Charge | size | stereo | glued | cleaned | " << std::endl;
            std::cout << "        |--- | ---| " << std::endl;
            headerStripsPrintedAlready = true;
          }
          std::cout << "        | " <<  globalIas_;
          if (detid.subdetId() == StripSubdetector::TIB) {
            std::cout << "        | TIB L" << abs(int(tTopo->tibLayer(detid)));
          }
          if (detid.subdetId() == StripSubdetector::TOB) {
            std::cout << "        | TOB L" << abs(int(tTopo->tobLayer(detid)));
          }
          else if (detid.subdetId() == StripSubdetector::TID) {
            std::cout << "        | TID D" << abs(int(tTopo->tidWheel(detid)));
          }
          else if (detid.subdetId() == StripSubdetector::TEC) {
            std::cout << "        | TEC D" << abs(int(tTopo->tidWheel(detid)));
          }
          
          std::cout << "        | " << genGammaBeta << " | " << track->eta() << " | " << stripNormCharge << " e/um | " << stripSize << " | " << tTopo->isStereo(detid) << " | " << isGlued << " | " << clusterCleaned << " | " << std::endl;
        }

      } // end of the strip part
    } // end the loop on the rechits
    }
    
    //Find the number of tracks passing selection for TOF<1 that will be used to check the background prediction
    //float Mass = -1;
    if (passPre && tof) {
      //compute the mass of the candidate, for TOF mass flip the TOF over 1 to get the mass, so 0.8->1.2
      float MassTOF = (tof) ? GetTOFMass(track->p(), (2 - tof->inverseBeta())) : -1;
      float InvBeta = GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_);
      float MassComb = -1;
      if (tof && dedxMObj)
        MassComb = GetMassFromBeta(track->p(),(InvBeta + (1 / (2 - tof->inverseBeta()))) * 0.5);
      if (dedxMObj)
        MassComb = Mass;
      if (tof)
        MassComb = GetMassFromBeta(track->p(), (1 / (2 - tof->inverseBeta())));
      //Background check looking at region with TOF<1
      for (unsigned int CutIndex = 0; CutIndex < CutPt_Flip_.size(); CutIndex++) {
        tuple->Mass_Flip->Fill(CutIndex, Mass, EventWeight_);
        if (tof && typeMode_ > 1) {
          tuple->MassTOF_Flip->Fill(CutIndex, MassTOF, EventWeight_);
          tuple->MassComb_Flip->Fill(CutIndex, MassComb, EventWeight_);
        }
      }
    }

    //compute the mass of the candidate
    float MassTOF = tof ?  GetTOFMass(track->p(), tof->inverseBeta()) : -1;
    float MassComb = tof   ?  GetMassFromBeta(track->p(), (1 / tof->inverseBeta())) : Mass;
    MassComb  = (tof && dedxMObj) ? GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5) : Mass;

    float MassUp = -1;
    if (dedxMUpObj)
      MassUp = GetMass(track->p(), dedxMUpObj->dEdx(), dEdxK_, dEdxC_);
    float MassUpComb = -1;
    if (tof && dedxMUpObj)
      MassUpComb =
          GetMassFromBeta(track->p(), (GetIBeta(dedxMUpObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
    if (dedxMUpObj)
      MassUpComb = MassUp;
    if (tof)
      MassUpComb = GetMassFromBeta(track->p(), (1 / tof->inverseBeta()));

    float MassDown = -1;
    if (dedxMDownObj)
      MassDown = GetMass(track->p(), dedxMDownObj->dEdx(), dEdxK_, dEdxC_);
    float MassDownComb = -1;
    if (tof && dedxMDownObj)
      MassDownComb =
          GetMassFromBeta(track->p(), (GetIBeta(dedxMDownObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
    if (dedxMDownObj)
      MassDownComb = MassDown;
    if (tof)
      MassDownComb = GetMassFromBeta(track->p(), (1 / tof->inverseBeta()));

    bool PassNonTrivialSelection = false;
  
    if (passPre) {
        tuple_maker->fillRegions(tuple,
                                 pT_cut,
                                 Ias_quantiles,
                                 track->eta(),
                                 track->p(),
                                 track->pt(),
                                 track->ptError(),
                                 dedxMObj ? dedxMObj->dEdx() : -1,
                                 dedxSObj ? dedxSObj->dEdx() : -1,
                                 Mass,
                                 tof ? tof->inverseBeta() : -99,
                                 EventWeight_);

      if (debug_ > 3  && trigInfo_ > 0) LogPrint(MOD) << "      >> We enter the selection cut loop now";
      //==========================================================
      // Cut loop: over all possible selection (one of them, the optimal one, will be used later)
      for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
        //Full Selection
        if (!passSelection(track, dedxSObj, dedxMObj, tof, iEvent, CutIndex, tuple, false, genBeta, false, 0, 0)) {
          continue;
        }

        if (CutIndex != 0) {
          PassNonTrivialSelection = true;
        }
        HSCPTk[CutIndex] = true;
        HSCPTk_SystHUp[CutIndex] = true;
        HSCPTk_SystHDown[CutIndex] = true;

        if (Mass > MaxMass[CutIndex])
          MaxMass[CutIndex] = Mass;
        if (MassUp > MaxMass_SystHUp[CutIndex])
          MaxMass_SystHUp[CutIndex] = Mass;
        if (MassDown > MaxMass_SystHDown[CutIndex])
          MaxMass_SystHDown[CutIndex] = Mass;

        //Fill Mass Histograms
        tuple->Mass->Fill(CutIndex, Mass, EventWeight_);
        if (tof && typeMode_ > 1) {
          tuple->MassTOF->Fill(CutIndex, MassTOF, EventWeight_);
          tuple->MassComb->Fill(CutIndex, MassComb, EventWeight_);
        }

        //Fill Mass Histograms for different Ih syst
        if (calcSyst_) {
          tuple->Mass_SystHUp->Fill(CutIndex, MassUp, EventWeight_);
          tuple->Mass_SystHDown->Fill(CutIndex, MassDown, EventWeight_);
        }
        if (tof && typeMode_ > 1) {
          tuple->MassTOF_SystH->Fill(CutIndex, MassTOF, EventWeight_);
          tuple->MassComb_SystHUp->Fill(CutIndex, MassUpComb, EventWeight_);
          tuple->MassComb_SystHDown->Fill(CutIndex, MassDownComb, EventWeight_);
        }

      }  //end of Cut loop
    } // end of condition for passPre

    float Ick2 = (dedxMObj) ? GetIck(dedxMObj->dEdx(), dEdxK_, dEdxC_) : 0.f;
    float genid = 0, gencharge = -99, genmass = -99, genpt = -99, geneta = -99, genphi = -99;

    if (isSignal) {
      genid = genColl[closestGenIndex].pdgId();
      gencharge = genColl[closestGenIndex].charge();
      genmass = genColl[closestGenIndex].mass();
      genpt = genColl[closestGenIndex].pt();
      geneta = genColl[closestGenIndex].eta();
      genphi = genColl[closestGenIndex].phi();
    }

    float iso_TK = -1;
    float iso_ECAL = -1;
    float iso_HCAL = -1;

    if (typeMode_ != 3) {
      const edm::ValueMap<susybsm::HSCPIsolation> IsolationMap = iEvent.get(hscpIsoToken_);
      susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
      iso_TK = hscpIso.Get_TK_SumEt();
      iso_ECAL = hscpIso.Get_ECAL_Energy();
      iso_HCAL = hscpIso.Get_HCAL_Energy();
    }

    float muon_PFIso03_sumCharHadPt = -1;
    float muon_PFIso03_sumNeutHadPt = -1;
    float muon_PFIso03_sumPhotonPt = -1;
    float muon_PFIso03_sumPUPt = -1;

    if (typeMode_ == 2 && !muon.isNull()) {
      muon_PFIso03_sumCharHadPt = muon->pfIsolationR03().sumChargedHadronPt;
      muon_PFIso03_sumNeutHadPt = muon->pfIsolationR03().sumNeutralHadronEt;
      muon_PFIso03_sumPhotonPt = muon->pfIsolationR03().sumPhotonEt;
      muon_PFIso03_sumPUPt = muon->pfIsolationR03().sumPUPt;
    }

    // HSCP type
    if ( hscp.type() == susybsm::HSCParticleType::globalMuon) {
      HSCP_type.push_back(0);
    } else if ( hscp.type() == susybsm::HSCParticleType::trackerMuon) {
      HSCP_type.push_back(1);
    } else if ( hscp.type() == susybsm::HSCParticleType::matchedStandAloneMuon) {
      HSCP_type.push_back(2);
    } else if ( hscp.type() == susybsm::HSCParticleType::standAloneMuon) {
      HSCP_type.push_back(3);
    } else if ( hscp.type() == susybsm::HSCParticleType::innerTrack) {
      HSCP_type.push_back(4);
    } else if ( hscp.type() == susybsm::HSCParticleType::unknown) {
      HSCP_type.push_back(5);
    }
    HSCP_passCutPt55.push_back(track->pt() > 55 ? true : false);
    HSCP_passPreselection.push_back(passPre);
    HSCP_passPreselectionSept8.push_back(passPreSept8);
    HSCP_passSelection.push_back(PassNonTrivialSelection);
    HSCP_Charge.push_back(track->charge());
    HSCP_Pt.push_back(track->pt());
    HSCP_PtErr.push_back(track->ptError());
    HSCP_Is_StripOnly.push_back(dedxIs_StripOnly ? dedxIs_StripOnly->dEdx() : -1);
    HSCP_Ias.push_back(dedxIas_FullTracker ? dedxIas_FullTracker->dEdx() : -1);
    HSCP_Ias_noPix_noTIB_noTID_no3TEC.push_back(dedxIas_noTIBnoTIDno3TEC ? dedxIas_noTIBnoTIDno3TEC->dEdx() : -1);
    HSCP_Ias_PixelOnly.push_back(dedxIas_PixelOnly ? dedxIas_PixelOnly->dEdx() : -1);
    HSCP_Ias_StripOnly.push_back(dedxIas_StripOnly ? dedxIas_StripOnly->dEdx() : -1);
    HSCP_Ias_PixelOnly_noL1.push_back(dedxIas_PixelOnly_noL1 ? dedxIas_PixelOnly_noL1->dEdx() : -1);
    HSCP_Ih.push_back(dedxMObj_FullTracker ? dedxMObj_FullTracker->dEdx() : -1);
    HSCP_Ick.push_back(dedxMObj ? Ick2 : -99);
    HSCP_Fmip.push_back(Fmip);
    HSCP_ProbXY.push_back(probXYonTrack);
    HSCP_ProbXY_noL1.push_back(probXYonTrackNoL1);
    HSCP_ProbQ.push_back(probQonTrack);
    HSCP_ProbQ_noL1.push_back(probQonTrackNoL1);
    HSCP_Ndof.push_back(track->ndof());
    HSCP_Chi2.push_back(track->chi2());
    HSCP_QualityMask.push_back(track->qualityMask());
    HSCP_isHighPurity.push_back(track->quality(reco::TrackBase::highPurity));
    HSCP_EoverP.push_back(pf_energy/track->p());
    HSCP_isMuon.push_back(pf_isMuon);
    HSCP_isPhoton.push_back(pf_isPhoton);
    HSCP_isElectron.push_back(pf_isElectron);
    HSCP_isChHadron.push_back(pf_isChHadron);
    HSCP_isNeutHadron.push_back(pf_isNeutHadron);
    HSCP_isPfTrack.push_back(pf_isPfTrack);
    HSCP_isUndefined.push_back(pf_isUndefined);
    HSCP_ECAL_energy.push_back(pf_ecal_energy);
    HSCP_HCAL_energy.push_back(pf_hcal_energy);
    HSCP_TOF.push_back(tof ? tof->inverseBeta() : -99);
    HSCP_TOFErr.push_back(tof ? tof->inverseBetaErr() : -99);
    HSCP_TOF_ndof.push_back(tof ? tof->nDof() : -99);
    HSCP_DTTOF.push_back(dttof ? dttof->inverseBeta() : -99);
    HSCP_DTTOFErr.push_back(dttof ? dttof->inverseBetaErr() : -99);
    HSCP_DTTOF_ndof.push_back(dttof ? dttof->nDof() : -99);
    HSCP_CSCTOF.push_back(csctof ? csctof->inverseBeta() : -99);
    HSCP_CSCTOFErr.push_back(csctof ? csctof->inverseBetaErr() : -99);
    HSCP_CSCTOF_ndof.push_back(csctof ? csctof->nDof() : -99);
    HSCP_Mass.push_back(Mass);
    HSCP_MassErr.push_back(MassErr);
    HSCP_dZ.push_back(dz);
    HSCP_dXY.push_back(dxy);
    HSCP_mT.push_back(massT);
    HSCP_dR.push_back(maxOpenAngle);
    HSCP_p.push_back(track->p());
    HSCP_eta.push_back(track->eta());
    HSCP_phi.push_back(track->phi());
    HSCP_NOH.push_back(track->found());
    HSCP_NOPH.push_back(nonL1PixHits);
    HSCP_FOVH.push_back(track->validFraction());
    HSCP_NOMH.push_back(missingHitsTillLast);
    HSCP_FOVHD.push_back(validFractionTillLast);
    HSCP_NOM.push_back(numDeDxHits);
    HSCP_matchTrigMuon_minDeltaR.push_back(dr_min_hlt_muon);
    HSCP_matchTrigMuon_pT.push_back(hlt_match_pt);
    HSCP_iso_TK.push_back(iso_TK);
    HSCP_iso_ECAL.push_back(iso_ECAL);
    HSCP_iso_HCAL.push_back(iso_HCAL);
    HSCP_track_genTrackMiniIsoSumPt.push_back(track_genTrackMiniIsoSumPt);
    HSCP_PFMiniIso_relative.push_back(miniRelIsoAll);
    HSCP_PFMiniIso_wMuon_relative.push_back(miniRelIsoAll_wMuon);
    HSCP_track_PFIsolationR005_sumChargedHadronPt.push_back(track_PFIso005_sumCharHadPt);
    HSCP_track_PFIsolationR005_sumNeutralHadronPt.push_back(track_PFIso005_sumNeutHadPt);
    HSCP_track_PFIsolationR005_sumPhotonPt.push_back(track_PFIso005_sumPhotonPt);
    HSCP_track_PFIsolationR005_sumPUPt.push_back(track_PFIso005_sumPUPt);
    HSCP_track_PFIsolationR01_sumChargedHadronPt.push_back(track_PFIso01_sumCharHadPt);
    HSCP_track_PFIsolationR01_sumNeutralHadronPt.push_back(track_PFIso01_sumNeutHadPt);
    HSCP_track_PFIsolationR01_sumPhotonPt.push_back(track_PFIso01_sumPhotonPt);
    HSCP_track_PFIsolationR01_sumPUPt.push_back(track_PFIso01_sumPUPt);
    HSCP_track_PFIsolationR03_sumChargedHadronPt.push_back(track_PFIso03_sumCharHadPt);
    HSCP_track_PFIsolationR03_sumNeutralHadronPt.push_back(track_PFIso03_sumNeutHadPt);
    HSCP_track_PFIsolationR03_sumPhotonPt.push_back(track_PFIso03_sumPhotonPt);
    HSCP_track_PFIsolationR03_sumPUPt.push_back(track_PFIso03_sumPUPt);
    HSCP_track_PFIsolationR05_sumChargedHadronPt.push_back(track_PFIso05_sumCharHadPt);
    HSCP_track_PFIsolationR05_sumNeutralHadronPt.push_back(track_PFIso05_sumNeutHadPt);
    HSCP_track_PFIsolationR05_sumPhotonPt.push_back(track_PFIso05_sumPhotonPt);
    HSCP_track_PFIsolationR05_sumPUPt.push_back(track_PFIso05_sumPUPt);
    HSCP_muon_PFIsolationR03_sumChargedHadronPt.push_back(muon_PFIso03_sumCharHadPt);
    HSCP_muon_PFIsolationR03_sumNeutralHadronPt.push_back(muon_PFIso03_sumNeutHadPt);
    HSCP_muon_PFIsolationR03_sumPhotonPt.push_back(muon_PFIso03_sumPhotonPt);
    HSCP_muon_PFIsolationR03_sumPUPt.push_back(muon_PFIso03_sumPUPt);
    HSCP_Ih_noL1.push_back(dedxIh_noL1 ? dedxIh_noL1->dEdx() : -1);
    HSCP_Ih_15drop.push_back(dedxIh_15drop ? dedxIh_15drop->dEdx() : -1);
    HSCP_Ih_StripOnly.push_back(dedxIh_StripOnly ? dedxIh_StripOnly->dEdx() : -1);
    HSCP_Ih_StripOnly_15drop.push_back(dedxIh_StripOnly_15drop ? dedxIh_StripOnly_15drop->dEdx() : -1);
    HSCP_Ih_PixelOnly_noL1.push_back(dedxIh_PixelOnlyh_noL1 ? dedxIh_PixelOnlyh_noL1->dEdx() : -1);
    HSCP_Ih_SaturationCorrectionFromFits.push_back(
        dedxIh_SaturationCorrectionFromFits ? dedxIh_SaturationCorrectionFromFits->dEdx() : -1);
    HSCP_clust_charge.push_back(clust_charge);
    HSCP_clust_pathlength.push_back(clust_pathlength);
    HSCP_clust_nstrip.push_back(clust_nstrip);
    HSCP_clust_sat254.push_back(clust_sat254);
    HSCP_clust_sat255.push_back(clust_sat255);
    HSCP_clust_detid.push_back(clust_detid);
    HSCP_clust_isStrip.push_back(clust_isStrip);
    HSCP_clust_isPixel.push_back(clust_isPixel);
    HSCP_GenId.push_back(genid);
    HSCP_GenCharge.push_back(gencharge);
    HSCP_GenMass.push_back(genmass);
    HSCP_GenPt.push_back(genpt);
    HSCP_GenEta.push_back(geneta);
    HSCP_GenPhi.push_back(genphi);

    HSCP_tuneP_Pt.push_back(tuneP_Pt);
    HSCP_tuneP_PtErr.push_back(tuneP_PtErr);
    HSCP_tuneP_Eta.push_back(tuneP_Eta);
    HSCP_tuneP_Phi.push_back(tuneP_Phi);
    HSCP_tuneP_MuonBestTrackType.push_back(tunePMuonBestTrackType);

    HSCP_ErrorHisto_bin.push_back(ErrorHisto_bin);
  }  //END loop over HSCP candidates
  
  if (passTechnicalChecks) {
    tuple->EventCutFlow->Fill(1.0, EventWeight_);
  }
  
  // Event level information, after choosing the best candidate track
  if (trigInfo_ > 0) tuple->BefPreS_NumCandidates->Fill(candidate_count, EventWeight_);
  if (trigInfo_ > 0) {
    tuple->PostPreS_NumCandidates->Fill(postPreS_candidate_count, EventWeight_);
    tuple->EventCutFlow->Fill(2.0, EventWeight_);
  }
  
  if (trigInfo_ > 0) {
    if (postPreS_candidate_count == 0) {
      tuple->PostS_MetOverHt_Cand0->Fill(RecoPFMET/pfJetHT, EventWeight_);
    } else if (postPreS_candidate_count == 1) {
      tuple->PostS_MetOverHt_Cand1->Fill(RecoPFMET/pfJetHT, EventWeight_);
    } else if (postPreS_candidate_count == 2) {
      tuple->PostS_MetOverHt_Cand2->Fill(RecoPFMET/pfJetHT, EventWeight_);
    }
  }
  
  if (!HLT_PFMET120_PFMHT120_IDTight && !HLT_PFHT500_PFMET100_PFMHT100_IDTight && !HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 && !HLT_MET105_IsoTrk50) {
    if (postPreS_candidate_count == 0) {
      tuple->PostS_TriggerMETallVsMetOverHt_Cand0->Fill(0.,RecoPFMET/pfJetHT);
    } else if (postPreS_candidate_count == 1) {
      tuple->PostS_TriggerMETallVsMetOverHt_Cand1->Fill(0.,RecoPFMET/pfJetHT);
    } else if (postPreS_candidate_count == 2) {
      tuple->PostS_TriggerMETallVsMetOverHt_Cand2->Fill(0.,RecoPFMET/pfJetHT);
    }
  } else if (HLT_PFMET120_PFMHT120_IDTight || HLT_PFHT500_PFMET100_PFMHT100_IDTight || HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 || HLT_MET105_IsoTrk50) {
    if (postPreS_candidate_count == 0) {
      tuple->PostS_TriggerMETallVsMetOverHt_Cand0->Fill(1.,RecoPFMET/pfJetHT);
    } else if (postPreS_candidate_count == 1) {
      tuple->PostS_TriggerMETallVsMetOverHt_Cand1->Fill(1.,RecoPFMET/pfJetHT);
    } else if (postPreS_candidate_count == 2) {
      tuple->PostS_TriggerMETallVsMetOverHt_Cand2->Fill(1.,RecoPFMET/pfJetHT);
    }
  }
  
  for (uint32_t cutFlowIndex = 1; cutFlowIndex<18; cutFlowIndex++) {
    if (cutFlowIndex <= bestSoFarCandCutInd) {
      tuple->EventCutFlow->Fill(cutFlowIndex+2, EventWeight_);
    }
  }
  
  if (trigInfo_ > 0 && postPreS_candidate_count == 0) {
    if (debug_ > 2) LogPrint(MOD) << "Trigger passed, but number of postPreSelected candidates is zero";
  }
  
  
  int hscpIndex = 0;
  for (const auto& hscp : iEvent.get(hscpToken_)) {
    hscpIndex++;
    if (bestCandidateIndex != hscpIndex) {
      continue;
    }
    
    if (debug_ > 3 && trigInfo_ > 0) LogPrint(MOD) << "After choosing the best candidate track";
    
    if ( hscp.type() == susybsm::HSCParticleType::globalMuon) {
      tuple->PostS_RecoHSCParticleType->Fill(0.);
    } else if ( hscp.type() == susybsm::HSCParticleType::trackerMuon) {
      tuple->PostS_RecoHSCParticleType->Fill(1.);
    } else if ( hscp.type() == susybsm::HSCParticleType::matchedStandAloneMuon) {
      tuple->PostS_RecoHSCParticleType->Fill(2.);
    } else if ( hscp.type() == susybsm::HSCParticleType::standAloneMuon) {
      tuple->PostS_RecoHSCParticleType->Fill(3.);
    } else if ( hscp.type() == susybsm::HSCParticleType::innerTrack) {
      tuple->PostS_RecoHSCParticleType->Fill(4.);
    } else if ( hscp.type() == susybsm::HSCParticleType::unknown) {
      tuple->PostS_RecoHSCParticleType->Fill(5.);
    }

    reco::TrackRef bestCandidateTrack =  hscp.trackRef();
    reco::MuonRef bestCandidateMuon = hscp.muonRef();
    
    if (trigInfo_ > 0) {
      tuple->PostS_HltMatchTrackLevel->Fill(1.0, EventWeight_);
      if (bestCandidateMuon.isNonnull()) {
        tuple->PostS_HltMatchTrackLevel->Fill(2.0, EventWeight_);
      }
      if (bestCandidateDrMinHltMuon < 0.15) {
        if (trigInfo_ > 0) tuple->PostS_HltMatchTrackLevel->Fill(3.0, EventWeight_);
      }
    }
    
    if (bestCandidateIas > 0.25 && bestCandidateProbQNoL1 < 0.1) {
        tuple->EventCutFlow->Fill(17.0, EventWeight_);
      if (bestCandidateTrack->pt() > 100) {
        tuple->EventCutFlow->Fill(18.0, EventWeight_);
      } if (bestCandidateTrack->pt() > 200) {
        tuple->EventCutFlow->Fill(19.0, EventWeight_);
      }
    }
    
    float shiftForPtValue = shiftForPt(bestCandidateTrack->pt(),bestCandidateTrack->eta(),bestCandidateTrack->phi(),bestCandidateTrack->charge());
    tuple->PostS_RelativePtShift->Fill(shiftForPtValue,  EventWeight_);
    
    float rescaledPtUp = bestCandidateTrack->pt()*(1+shiftForPtValue);
    float rescaledPtDown = bestCandidateTrack->pt()*(1-shiftForPtValue);
    
    tuple->PostS_Ias->Fill(bestCandidateIas, EventWeight_);
    tuple->PostS_ProbQNoL1->Fill(1 - bestCandidateProbQNoL1, EventWeight_);
    tuple->PostS_ProbQNoL1VsIas->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas, EventWeight_);
    tuple->PostS_ProbQNoL1VsIasVsPt->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas, bestCandidateTrack->pt(), EventWeight_);
    
    // Systematics plots for pT rescaling
    if (rescaledPtUp > globalMinPt_) {
      tuple->PostS_ProbQNoL1VsIas_Pt_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
    }
    if (rescaledPtDown > globalMinPt_) {
      tuple->PostS_ProbQNoL1VsIas_Pt_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
    }
    
    // Systematics plots for Gi rescaling
    tuple->PostS_ProbQNoL1VsIas_Ias_up->Fill(1 - bestCandidateProbQNoL1, std::min(1.0,bestCandidateIas*1.05),  EventWeight_);
    tuple->PostS_ProbQNoL1VsIas_Ias_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas*0.95,  EventWeight_);
    
    // Systematics plots for PU rescaling
    tuple->PostS_ProbQNoL1VsIas_Pileup_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[0]);
    tuple->PostS_ProbQNoL1VsIas_Pileup_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[1]);
    
    // Systematics plots for Fi rescaling
    tuple->PostS_ProbQNoL1VsIas_ProbQNoL1_up->Fill(std::min(1.,(1 - bestCandidateProbQNoL1)*1.005), bestCandidateIas,  EventWeight_);
    tuple->PostS_ProbQNoL1VsIas_ProbQNoL1_down->Fill((1 - bestCandidateProbQNoL1)*0.995, bestCandidateIas,  EventWeight_);
    
    // Systematics plots for trigger rescaling
    float triggerSystFactorUp = triggerSystFactor(bestCandidateTrack->eta(),bestCandidateGenBeta,+1);
    float triggerSystFactorDown = triggerSystFactor(bestCandidateTrack->eta(),bestCandidateGenBeta,-1);

    tuple->PostS_ProbQNoL1VsIas_Trigger_up->Fill(std::min(1.f,(1 - bestCandidateProbQNoL1)), bestCandidateIas,  EventWeight_ * triggerSystFactorUp);
    tuple->PostS_ProbQNoL1VsIas_Trigger_down->Fill((1 - bestCandidateProbQNoL1), bestCandidateIas,  EventWeight_  * triggerSystFactorDown);
    
    if (bestCandidateTrack->pt() >= pT_cut && bestCandidateTrack->pt() < 100) {
      tuple->PostS_SR1_Ias->Fill(bestCandidateIas, EventWeight_);
      tuple->PostS_SR1_ProbQNoL1->Fill(1 - bestCandidateProbQNoL1, EventWeight_);
      tuple->PostS_SR1_ProbQNoL1VsIas->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas, EventWeight_);
      
        // Systematics plots for pT rescaling
      if (rescaledPtUp >= pT_cut && bestCandidateTrack->pt() < 100) {
        tuple->PostS_SR1_ProbQNoL1VsIas_Pt_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
      }
      if (rescaledPtDown >= pT_cut && bestCandidateTrack->pt() < 100) {
        tuple->PostS_SR1_ProbQNoL1VsIas_Pt_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
      }
      
        // Systematics plots for Gi rescaling
      tuple->PostS_SR1_ProbQNoL1VsIas_Ias_up->Fill(1 - bestCandidateProbQNoL1, std::min(1.0,bestCandidateIas*1.05),  EventWeight_);
      tuple->PostS_SR1_ProbQNoL1VsIas_Ias_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas*0.95,  EventWeight_);
      
        // Systematics plots for PU rescaling
      tuple->PostS_SR1_ProbQNoL1VsIas_Pileup_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[0]);
      tuple->PostS_SR1_ProbQNoL1VsIas_Pileup_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[1]);
      
        // Systematics plots for Fi rescaling
      tuple->PostS_SR1_ProbQNoL1VsIas_ProbQNoL1_up->Fill(std::min(1.,(1 - bestCandidateProbQNoL1)*1.005), bestCandidateIas,  EventWeight_);
      tuple->PostS_SR1_ProbQNoL1VsIas_ProbQNoL1_down->Fill((1 - bestCandidateProbQNoL1)*0.995, bestCandidateIas,  EventWeight_);
      
        // Systematics plots for trigger rescaling
      tuple->PostS_SR1_ProbQNoL1VsIas_Trigger_up->Fill(std::min(1.f,(1 - bestCandidateProbQNoL1)), bestCandidateIas,  EventWeight_ * triggerSystFactorUp);
      tuple->PostS_SR1_ProbQNoL1VsIas_Trigger_down->Fill((1 - bestCandidateProbQNoL1), bestCandidateIas,  EventWeight_  * triggerSystFactorDown);
    }
    // Repeat for several SRs
    if (bestCandidateTrack->pt() >= 100 && bestCandidateTrack->pt() < 200) {
      tuple->PostS_SR2_Ias->Fill(bestCandidateIas, EventWeight_);
      tuple->PostS_SR2_ProbQNoL1->Fill(1 - bestCandidateProbQNoL1, EventWeight_);
      tuple->PostS_SR2_ProbQNoL1VsIas->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas, EventWeight_);
      
        // Systematics plots for pT rescaling
      if (rescaledPtUp >= 100 && rescaledPtUp < 200) {
        tuple->PostS_SR2_ProbQNoL1VsIas_Pt_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
      }
      if (rescaledPtUp >= 100 && rescaledPtUp < 200) {
        tuple->PostS_SR2_ProbQNoL1VsIas_Pt_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
      }
      
      // Systematics plots for Gi rescaling
      tuple->PostS_SR2_ProbQNoL1VsIas_Ias_up->Fill(1 - bestCandidateProbQNoL1, std::min(1.0,bestCandidateIas*1.05),  EventWeight_);
      tuple->PostS_SR2_ProbQNoL1VsIas_Ias_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas*0.95,  EventWeight_);
      
      // Systematics plots for PU rescaling
      tuple->PostS_SR2_ProbQNoL1VsIas_Pileup_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[0]);
      tuple->PostS_SR2_ProbQNoL1VsIas_Pileup_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[1]);
      
      // Systematics plots for Fi rescaling
      tuple->PostS_SR2_ProbQNoL1VsIas_ProbQNoL1_up->Fill(std::min(1.,(1 - bestCandidateProbQNoL1)*1.005), bestCandidateIas,  EventWeight_);
      tuple->PostS_SR2_ProbQNoL1VsIas_ProbQNoL1_down->Fill((1 - bestCandidateProbQNoL1)*0.995, bestCandidateIas,  EventWeight_);
      
      // Systematics plots for trigger rescaling
      tuple->PostS_SR2_ProbQNoL1VsIas_Trigger_up->Fill(std::min(1.f,(1 - bestCandidateProbQNoL1)), bestCandidateIas,  EventWeight_ * triggerSystFactorUp);
      tuple->PostS_SR2_ProbQNoL1VsIas_Trigger_down->Fill((1 - bestCandidateProbQNoL1), bestCandidateIas,  EventWeight_  * triggerSystFactorDown);
    }
    if (bestCandidateTrack->pt() >= 200) {
      tuple->PostS_SR3_Ias->Fill(bestCandidateIas, EventWeight_);
      tuple->PostS_SR3_ProbQNoL1->Fill(1 - bestCandidateProbQNoL1, EventWeight_);
      tuple->PostS_SR3_ProbQNoL1VsIas->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas, EventWeight_);
      
        // Systematics plots for pT rescaling
      if (rescaledPtUp >= 200) {
        tuple->PostS_SR3_ProbQNoL1VsIas_Pt_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
      }
      if (rescaledPtUp >= 200) {
        tuple->PostS_SR3_ProbQNoL1VsIas_Pt_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_);
      }
      
      // Systematics plots for Gi rescaling
      tuple->PostS_SR3_ProbQNoL1VsIas_Ias_up->Fill(1 - bestCandidateProbQNoL1, std::min(1.0,bestCandidateIas*1.05),  EventWeight_);
      tuple->PostS_SR3_ProbQNoL1VsIas_Ias_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas*0.95,  EventWeight_);
      
      // Systematics plots for PU rescaling
      tuple->PostS_SR3_ProbQNoL1VsIas_Pileup_up->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[0]);
      tuple->PostS_SR3_ProbQNoL1VsIas_Pileup_down->Fill(1 - bestCandidateProbQNoL1, bestCandidateIas,  EventWeight_ * PUSystFactor_[1]);
      
      // Systematics plots for Fi rescaling
      tuple->PostS_SR3_ProbQNoL1VsIas_ProbQNoL1_up->Fill(std::min(1.,(1 - bestCandidateProbQNoL1)*1.005), bestCandidateIas,  EventWeight_);
      tuple->PostS_SR3_ProbQNoL1VsIas_ProbQNoL1_down->Fill((1 - bestCandidateProbQNoL1)*0.995, bestCandidateIas,  EventWeight_);
      
      // Systematics plots for trigger rescaling
      tuple->PostS_SR3_ProbQNoL1VsIas_Trigger_up->Fill(std::min(1.f,(1 - bestCandidateProbQNoL1)), bestCandidateIas,  EventWeight_ * triggerSystFactorUp);
      tuple->PostS_SR3_ProbQNoL1VsIas_Trigger_down->Fill((1 - bestCandidateProbQNoL1), bestCandidateIas,  EventWeight_  * triggerSystFactorDown);
    }
  } // end of "loop" on the single best HSCP candidate
  
  // Trigger type after preSelection at the event level
  // I need these two lines because the OR of MET/Muon will be overwriting the trigInfo_
  // so we never get 1 or 2
  if (postPreS_candidate_count > 0) {
    if (muTrig) {
      tuple->PostPreS_TriggerType->Fill(1., EventWeight_);
    }
    if (metTrig) {
      tuple->PostPreS_TriggerType->Fill(2., EventWeight_);
    }
    tuple->PostPreS_TriggerType->Fill(trigInfo_, EventWeight_);
  }

  // Trigger type after preSelection at the event level
  tuple->PostPreS_TriggerType->Fill(trigInfo_, EventWeight_);


  // Calibration information
  // loop on the tracks
  bool apply_trigger_sel_for_calib = true;
  if (typeMode_ != 3 && apply_trigger_sel_for_calib && trigInfo_ > 0) {
    for(unsigned int c=0;c<trackCollectionHandle->size();c++){
      reco::TrackRef generalTrack = reco::TrackRef( trackCollectionHandle.product(), c );

      // selection of the tracks

      bool passedTrackCutsArray[8];
      std::fill(std::begin(passedTrackCutsArray), std::end(passedTrackCutsArray),true);

      // Check if eta is inside the max eta cut for detector homogeneity
      passedTrackCutsArray[0]  = (fabs(generalTrack->eta()) < globalMaxEta_) ? true : false;
      // Select only high purity tracks to ensure good quality tracks
      passedTrackCutsArray[1]  = (generalTrack->quality(reco::TrackBase::highPurity)) ? true : false;
      // Cut on the chi2 / ndof to ensure good quality tracks
      passedTrackCutsArray[2] = ((generalTrack->chi2() / generalTrack->ndof()) < globalMaxChi2_) ? true : false;
      // Check the min fraction of valid hits to ensure good stats on the hits
      passedTrackCutsArray[3]  = (generalTrack->validFraction() > globalMinFOVH_) ? true : false;
      // Cut on the impact parameter to ensure the track is coming from the PV
      // for typeMode_ 3 (TOF only) dz and dxy is supposed to come from the beamspot
      float dzForCalib = generalTrack->dz(highestSumPt2Vertex.position());
      float dxyForCalib = generalTrack->dxy(highestSumPt2Vertex.position());
      passedTrackCutsArray[4] = (fabs(dzForCalib) < 0.5) ? true : false;
      passedTrackCutsArray[5] = (fabs(dxyForCalib) < 0.5)  ? true : false;

      if (!passPreselection(passedTrackCutsArray, false)) continue;


      //load dedx quantity associated to this track
      const reco::DeDxHitInfo* dedxHits = nullptr;
      reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(generalTrack.key());
      if (!dedxHitsRef.isNull())
        dedxHits = &(*dedxHitsRef);

//      if(!dedxHits)  cout << "track " << c << " pt : " <<  generalTrack->pt()  << ", eta : " << generalTrack->eta()  << " without dedxHit info " << endl;
      if(!dedxHits)continue;
//      cout << "track " << c << " pt : " <<  generalTrack->pt()  << ", eta : " << generalTrack->eta()  << " with dedxHit info " << endl;

      reco::DeDxHitInfoRef dedxHitsPrescRef;
      dedxHitsPrescRef = dedxCollH->get(generalTrack.key());
      int preScaleForDeDx = 1;
      if (!dedxHitsPrescRef.isNull()) {
        preScaleForDeDx = (*dedxPrescCollH)[dedxHitsPrescRef];
      }


      float dEdxErr = 0;

      // Ih no pixel L1
      auto dedxIh_noL1_Tmp =
        computedEdx(generalTrack->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = true, useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, &dEdxErr, useTemplateLayer_, skipPixelL1 = true);

      reco::DeDxData* dedxIh_noL1 = dedxIh_noL1_Tmp.numberOfMeasurements() > 0 ? &dedxIh_noL1_Tmp : nullptr;

      // Ih Strip only
      auto dedxIh_StripOnly_Tmp =
        computedEdx(generalTrack->eta(),run_number, year, dedxHits, dEdxSF, localdEdxTemplates = nullptr, usePixel = false, useStrip = true, useClusterCleaning, useTruncated = false,
                    mustBeInside, MaxStripNOM, correctFEDSat, crossTalkInvAlgo = 1, dropLowerDeDxValue = 0.0, &dEdxErr, useTemplateLayer_);

      reco::DeDxData* dedxIh_StripOnly = dedxIh_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIh_StripOnly_Tmp : nullptr;


      unsigned int nonL1PixHits = dedxIh_noL1_Tmp.numberOfMeasurements() - dedxIh_StripOnly_Tmp.numberOfMeasurements() ;
      unsigned int numDeDxHits = dedxIh_noL1_Tmp.numberOfMeasurements() ;

      // Check the number of non-layer-1 pixel hits to ensure good stats on the hits
      passedTrackCutsArray[6]  = (typeMode_ != 3 && nonL1PixHits >= globalMinNOPH_) ? true : false;
      // Cut for the number of dEdx hits to ensure good stats on the hits
      passedTrackCutsArray[7]  = (numDeDxHits >= globalMinNOM_)  ? true : false;
      if (!passPreselection(passedTrackCutsArray, false)) continue;

      float ih0_noL1 = (dedxIh_noL1) ?  dedxIh_noL1->dEdx() : 0.0;
      float ih0_strip = (dedxIh_StripOnly) ?  dedxIh_StripOnly->dEdx() : 0.0;

      if (generalTrack->p() < 50) {
        tuple->K_and_C_Ih_noL1_VsP_loose1->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        tuple->K_and_C_Ih_noL1_VsP_loose2->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        if (fabs(generalTrack->eta()) < 1) {
             tuple->K_and_C_Ih_noL1_VsP_eta1_loose1->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
             tuple->K_and_C_Ih_noL1_VsP_eta1_loose2->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        }
        else {
             tuple->K_and_C_Ih_noL1_VsP_eta2_loose1->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
             tuple->K_and_C_Ih_noL1_VsP_eta2_loose2->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        }
        tuple->K_and_C_Ih_strip_VsP_loose1->Fill(generalTrack->p(),ih0_strip, preScaleForDeDx*EventWeight_);
        tuple->K_and_C_Ih_strip_VsP_loose2->Fill(generalTrack->p(),ih0_strip, preScaleForDeDx*EventWeight_);

       if ( fabs(dzForCalib) < globalMaxDZ_ && fabs(dxyForCalib) < globalMaxDXY_) {

        tuple->K_and_C_Ih_noL1_VsP_1->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        tuple->K_and_C_Ih_noL1_VsP_2->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        if (fabs(generalTrack->eta())<1) {
             tuple->K_and_C_Ih_noL1_VsP_eta1_1->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
             tuple->K_and_C_Ih_noL1_VsP_eta1_2->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        }
        else {
             tuple->K_and_C_Ih_noL1_VsP_eta2_1->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
             tuple->K_and_C_Ih_noL1_VsP_eta2_2->Fill(generalTrack->p(),ih0_noL1, preScaleForDeDx*EventWeight_);
        }
        tuple->K_and_C_Ih_strip_VsP_1->Fill(generalTrack->p(),ih0_strip, preScaleForDeDx*EventWeight_);
        tuple->K_and_C_Ih_strip_VsP_2->Fill(generalTrack->p(),ih0_strip, preScaleForDeDx*EventWeight_);
       }
      }

      // Loop on the dedxHits belonging to the generalTrack under study
      for(unsigned int h=0;h< dedxHits->size();h++){
           DetId detid(dedxHits->detId(h));
           float dedx_charge=dedxHits->charge(h);
           float dedx_pathlength=dedxHits->pathlength(h);
           float factorChargeToE = (detid.subdetId() < 3) ? 3.61e-06 : 3.61e-06 * 265;
           float scaleFactor =dEdxSF[0];
           if (detid.subdetId() >=3 ) {
              // strip

              // saturation correction of the charge
              const SiStripCluster* cluster = dedxHits->stripCluster(h);
              std::vector<int> amplitudes = convert(cluster->amplitudes());
              amplitudes = SaturationCorrection(amplitudes,0.10,0.04,true,20,25);
              dedx_charge = 0;
              for (unsigned int s = 0; s < amplitudes.size(); s++) {
                 dedx_charge+=amplitudes[s];
              }
              float charge_over_pathlength = dedx_charge * scaleFactor * factorChargeToE / dedx_pathlength;

              // cleaning cuts
              std::vector<int> amplitudes2 = convert(cluster->amplitudes());
              std::vector<int> amplitudesPrim = CrossTalkInv(amplitudes2, 0.10, 0.04, true);
              bool cleaning = clusterCleaning(amplitudesPrim, 1);
              bool dedx_inside = isHitInsideTkModule(dedxHits->pos(h), dedxHits->detId(h), cluster);

              if (cleaning && dedx_inside)  {
                 if (fabs(generalTrack->eta()<0.4)){
                   tuple->SF_HHit2DStrip_loose->Fill(generalTrack->p(), charge_over_pathlength, preScaleForDeDx*EventWeight_);
                   if ( fabs(dzForCalib) < globalMaxDZ_ && fabs(dxyForCalib) < globalMaxDXY_)
                     tuple->SF_HHit2DStrip->Fill(generalTrack->p(), charge_over_pathlength, preScaleForDeDx*EventWeight_);
                 }
              }
           }
           else {
              // pixel
              scaleFactor *= dEdxSF[1];
              //correction de tamas
              float pixelScaling = GetSFPixel(detid.subdetId(), detid, year, run_number);
              scaleFactor *= pixelScaling;
              float charge_over_pathlength = dedx_charge * scaleFactor * factorChargeToE / dedx_pathlength;
              // BPIXL1
              bool isBPIXL1=false;
              if (detid.subdetId() == 1 && ((detid >> 20) & 0xF) == 1) isBPIXL1=true;
              if (!isBPIXL1) {
                 if (fabs(generalTrack->eta()<0.4)){
                   tuple->SF_HHit2DPix_loose->Fill(generalTrack->p(), charge_over_pathlength, preScaleForDeDx*EventWeight_);
                   if ( fabs(dzForCalib) < globalMaxDZ_ && fabs(dxyForCalib) < globalMaxDXY_)
                      tuple->SF_HHit2DPix->Fill(generalTrack->p(), charge_over_pathlength, preScaleForDeDx*EventWeight_);
                 }
              }
           }

      } // end dEdX loop
    } // end general track loop
   } // end if selection before looping on track


  // end of the Calibration part

  tuple_maker->fillTreeBranches(tuple,
                                trigInfo_,
                                iEvent.id().run(),
                                iEvent.id().event(),
                                iEvent.id().luminosityBlock(),
                                pileup_fromLumi,
                                bunchXing,
                                nPU,
                                nPUmean,
                                vertexColl.size(),
                                numGoodVerts,
                                pvX,
                                pvY,
                                pvZ,
                                pvRho,
                                pvNdof,
                                pvChi2,
                                pvSumPt2,
                                HSCP_count,
                                nMuons,
                                pfJetsNum,
                                EventWeight_,
                                GeneratorWeight_,
                                GeneratorBinningValues_,
                                triggerDecision,
                                triggerHLTPrescale,
                                triggerObjectE,
                                triggerObjectPt,
                                triggerObjectEta,
                                triggerObjectPhi,
                                HLT_Mu50,
                                HLT_PFMET120_PFMHT120_IDTight,
                                HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                                HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                                HLT_MET105_IsoTrk50,
                                RecoCaloMET,
                                RecoCaloMET_phi,
                                RecoCaloMET_sigf,
                                RecoPFMET,
                                RecoPFMET_phi,
                                RecoPFMET_sigf,
                                RecoPFMHT,
                                HLTCaloMET,
                                HLTCaloMET_phi,
                                HLTCaloMET_sigf,
                                HLTCaloMETClean,
                                HLTCaloMETClean_phi,
                                HLTCaloMETClean_sigf,
                                HLTCaloMHT,
                                HLTCaloMHT_phi,
                                HLTCaloMHT_sigf,
                                HLTPFMET,
                                HLTPFMET_phi,
                                HLTPFMET_sigf,
                                HLTPFMHT,
                                HLTPFMHT_phi,
                                HLTPFMHT_sigf,
                                matchedMuonWasFound,
                                gParticleId,
                                gParticleStatus,
                                gParticleE,
                                gParticlePt,
                                gParticlePz,
                                gParticleEta,
                                gParticlePhi,
                                gParticleBeta,
                                gParticleCharge,
                                gParticleProdVertexX,
                                gParticleProdVertexY,
                                gParticleProdVertexZ,
                                gParticleMotherId,
                                gParticleMotherIndex,
                                eleE,
                                elePt,
                                eleEta,
                                elePhi,
                                eleCharge,
                                eleE_SC,
                                eleEta_SC,
                                elePhi_SC,
                                eleSigmaIetaIeta,
                                eleFull5x5SigmaIetaIeta,
                                eleR9,
                                ele_dEta,
                                ele_dPhi,
                                ele_HoverE,
                                ele_d0,
                                ele_dZ,
                                ele_pileupIso,
                                ele_chargedIso,
                                ele_photonIso,
                                ele_neutralHadIso,
                                ele_MissHits,
                                ele_passCutBasedIDVeto,
                                ele_passCutBasedIDLoose,
                                ele_passCutBasedIDMedium,
                                ele_passCutBasedIDTight,
                                ele_passMVAIsoIDWP80,
                                ele_passMVAIsoIDWP90,
                                ele_passMVAIsoIDWPHZZ,
                                ele_passMVAIsoIDWPLoose,
                                ele_passMVANoIsoIDWP80,
                                ele_passMVANoIsoIDWP90,
                                ele_passMVANoIsoIDWPLoose,
                                ele_PassConvVeto,
                                ele_OneOverEminusOneOverP,
                                muonE,
                                muonPt,
                                muonEta,
                                muonPhi,
                                muonCharge,
                                muonIsLoose,
                                muonIsMedium,
                                muonIsTight,
                                muon_d0,
                                muon_d0Err,
                                muon_dZ,
                                muon_ip3d,
                                muon_ip3dSignificance,
                                muonType,
                                muonQuality,
                                muon_pileupIso,
                                muon_chargedIso,
                                muon_photonIso,
                                muon_neutralHadIso,
                                muon_validFractionTrackerHits,
                                muon_normChi2,
                                muon_chi2LocalPosition,
                                muon_kinkFinder,
                                muon_segmentCompatability,
                                muon_trkIso,
                                muon_tuneP_Pt,
                                muon_tuneP_PtErr,
                                muon_tuneP_Eta,
                                muon_tuneP_Phi,
                                muon_tuneP_MuonBestTrackType,
                                muon_isHighPtMuon,
                                muon_isTrackerHighPtMuon,
                                Jets_pt,
                                Jets_eta,
                                Jets_phi,
                                Jets_mass,
                                Jets_E,
                                Jets_pdgId,
                                Jets_et,
                                Jets_chargedEmEnergyFraction,
                                Jets_neutralEmEnergyFraction,
                                Jets_chargedHadronEnergyFraction,
                                Jets_neutralHadronEnergyFraction,
                                Jets_muonEnergyFraction,
                                Jets_chargedMultiplicity,
                                Jets_neutralMultiplicity,
                                Jets_jetArea,
                                Jets_pileupE,
                                HSCP_mT,
                                HSCP_passCutPt55,
                                HSCP_passPreselection,
                                HSCP_passPreselectionSept8,
                                HSCP_passSelection,
                                HSCP_isPFMuon,
                                HSCP_PFMuonPt,
                                HSCP_Charge,
                                HSCP_Pt,
                                HSCP_PtErr,
                                HSCP_Is_StripOnly,
                                HSCP_Ias,
                                HSCP_Ias_noPix_noTIB_noTID_no3TEC,
                                HSCP_Ias_PixelOnly,
                                HSCP_Ias_StripOnly,
                                HSCP_Ias_PixelOnly_noL1,
                                HSCP_Ih,
                                HSCP_Ick,
                                HSCP_Fmip,
                                HSCP_ProbXY,
                                HSCP_ProbXY_noL1,
                                HSCP_ProbQ,
                                HSCP_ProbQ_noL1,
                                HSCP_Ndof,
                                HSCP_Chi2,
                                HSCP_QualityMask,
                                HSCP_isHighPurity,
                                HSCP_EoverP,
                                HSCP_isMuon,
                                HSCP_isPhoton,
                                HSCP_isElectron,
                                HSCP_isChHadron,
                                HSCP_isNeutHadron,
                                HSCP_isPfTrack,
                                HSCP_isUndefined,
                                HSCP_ECAL_energy,
                                HSCP_HCAL_energy,
                                HSCP_TOF,
                                HSCP_TOFErr,
                                HSCP_TOF_ndof,
                                HSCP_DTTOF,
                                HSCP_DTTOFErr,
                                HSCP_DTTOF_ndof,
                                HSCP_CSCTOF,
                                HSCP_CSCTOFErr,
                                HSCP_CSCTOF_ndof,
                                HSCP_Mass,
                                HSCP_MassErr,
                                HSCP_dZ,
                                HSCP_dXY,
                                HSCP_dR,
                                HSCP_p,
                                HSCP_eta,
                                HSCP_phi,
                                HSCP_NOH,
                                HSCP_NOPH,
                                HSCP_FOVH,
                                HSCP_NOMH,
                                HSCP_FOVHD,
                                HSCP_NOM,
                                HSCP_matchTrigMuon_minDeltaR,
                                HSCP_matchTrigMuon_pT,
                                HSCP_iso_TK,
                                HSCP_iso_ECAL,
                                HSCP_iso_HCAL,
                                HSCP_track_genTrackMiniIsoSumPt,
                                HSCP_PFMiniIso_relative,
                                HSCP_PFMiniIso_wMuon_relative,
                                HSCP_track_PFIsolationR005_sumChargedHadronPt,
                                HSCP_track_PFIsolationR005_sumNeutralHadronPt,
                                HSCP_track_PFIsolationR005_sumPhotonPt,
                                HSCP_track_PFIsolationR005_sumPUPt,
                                HSCP_track_PFIsolationR01_sumChargedHadronPt,
                                HSCP_track_PFIsolationR01_sumNeutralHadronPt,
                                HSCP_track_PFIsolationR01_sumPhotonPt,
                                HSCP_track_PFIsolationR01_sumPUPt,
                                HSCP_track_PFIsolationR03_sumChargedHadronPt,
                                HSCP_track_PFIsolationR03_sumNeutralHadronPt,
                                HSCP_track_PFIsolationR03_sumPhotonPt,
                                HSCP_track_PFIsolationR03_sumPUPt,
                                HSCP_track_PFIsolationR05_sumChargedHadronPt,
                                HSCP_track_PFIsolationR05_sumNeutralHadronPt,
                                HSCP_track_PFIsolationR05_sumPhotonPt,
                                HSCP_track_PFIsolationR05_sumPUPt,
                                HSCP_muon_PFIsolationR03_sumChargedHadronPt,
                                HSCP_muon_PFIsolationR03_sumNeutralHadronPt,
                                HSCP_muon_PFIsolationR03_sumPhotonPt,
                                HSCP_muon_PFIsolationR03_sumPUPt,
                                HSCP_Ih_noL1,
                                HSCP_Ih_15drop,
                                HSCP_Ih_StripOnly,
                                HSCP_Ih_StripOnly_15drop,
                                HSCP_Ih_PixelOnly_noL1,
                                HSCP_Ih_SaturationCorrectionFromFits,
                                HSCP_clust_charge,
                                HSCP_clust_pathlength,
                                HSCP_clust_nstrip,
                                HSCP_clust_sat254,
                                HSCP_clust_sat255,
                                HSCP_clust_detid,
                                HSCP_clust_isStrip,
                                HSCP_clust_isPixel,
                                HSCP_GenId,
                                HSCP_GenCharge,
                                HSCP_GenMass,
                                HSCP_GenPt,
                                HSCP_GenEta,
                                HSCP_GenPhi,
                                HSCP_tuneP_Pt,
                                HSCP_tuneP_PtErr,
                                HSCP_tuneP_Eta,
                                HSCP_tuneP_Phi,
                                HSCP_tuneP_MuonBestTrackType,
                                HSCP_ErrorHisto_bin,
                                HSCP_type
                              );

  //save event dependent information thanks to the bookkeeping
  for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
    if (HSCPTk[CutIndex]) {
      tuple->HSCPE->Fill(CutIndex, EventWeight_);
      tuple->MaxEventMass->Fill(CutIndex, MaxMass[CutIndex], EventWeight_);
      if (isBckg) {
        tuple->HSCPE->Fill(CutIndex, EventWeight_);
        tuple->MaxEventMass->Fill(CutIndex, MaxMass[CutIndex], EventWeight_);
      }
    }
    if (calcSyst_) {
      if (HSCPTk_SystP[CutIndex]) {
        tuple->HSCPE_SystP->Fill(CutIndex, EventWeight_);
        tuple->MaxEventMass_SystP->Fill(CutIndex, MaxMass_SystP[CutIndex], EventWeight_);
      }
      if (HSCPTk_SystI[CutIndex]) {
        tuple->HSCPE_SystI->Fill(CutIndex, EventWeight_);
        tuple->MaxEventMass_SystI->Fill(CutIndex, MaxMass_SystI[CutIndex], EventWeight_);
      }
      if (HSCPTk_SystM[CutIndex]) {
        tuple->HSCPE_SystM->Fill(CutIndex, EventWeight_);
        tuple->MaxEventMass_SystM->Fill(CutIndex, MaxMass_SystM[CutIndex], EventWeight_);
      }
      if (HSCPTk_SystT[CutIndex] && typeMode_ > 1) {
        tuple->HSCPE_SystT->Fill(CutIndex, EventWeight_);
        tuple->MaxEventMass_SystT->Fill(CutIndex, MaxMass_SystT[CutIndex], EventWeight_);
      }
      if (HSCPTk_SystPU[CutIndex]) {
        tuple->HSCPE_SystPU->Fill(CutIndex, EventWeight_ * PUSystFactor_[0]);
        tuple->MaxEventMass_SystPU->Fill(CutIndex, MaxMass_SystPU[CutIndex], EventWeight_ * PUSystFactor_[0]);
      }
      if (HSCPTk_SystHUp[CutIndex]) {
        tuple->HSCPE_SystHUp->Fill(CutIndex, EventWeight_);
        tuple->MaxEventMass_SystHUp->Fill(CutIndex, MaxMass_SystHUp[CutIndex], EventWeight_);
      }
      if (HSCPTk_SystHDown[CutIndex]) {
        tuple->HSCPE_SystHDown->Fill(CutIndex, EventWeight_);
        tuple->MaxEventMass_SystHDown->Fill(CutIndex, MaxMass_SystHDown[CutIndex], EventWeight_);
      }
    } // end of calcSyst_
  }
} // end of analyze()

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() {
  delete RNG;
  delete RNG2;
  delete tuple;
  if (!isData) {
    delete mcWeight;
  }
  delete[] HSCPTk;
  delete[] HSCPTk_SystP;
  delete[] HSCPTk_SystI;
  delete[] HSCPTk_SystT;
  delete[] HSCPTk_SystM;
  delete[] HSCPTk_SystPU;
  delete[] HSCPTk_SystHUp;
  delete[] HSCPTk_SystHDown;
  delete[] MaxMass;
  delete[] MaxMass_SystP;
  delete[] MaxMass_SystI;
  delete[] MaxMass_SystT;
  delete[] MaxMass_SystM;
  delete[] MaxMass_SystPU;
  delete[] MaxMass_SystHUp;
  delete[] MaxMass_SystHDown;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setComment("Analyzer for HSCP search");
  desc.add("HscpCollection", edm::InputTag("HSCParticleProducer"))
    ->setComment("Input collection for HSCP candidate");
  desc.add("HscpIsoCollection", edm::InputTag("HSCPIsolation", "R03"))
    ->setComment("Input collection for HSCP isolation");
  desc.add("DedxCollection", edm::InputTag("dedxHitInfo"))
    ->setComment("Input collection for dEdx hit information");
  desc.add("DedxCollectionPrescale", edm::InputTag("dedxHitInfo","prescale"))
  ->setComment("Input collection for prescaled dEdx hit information");
  desc.add("MuonTimeCollection", edm::InputTag("muons", "combined"))
    ->setComment("Input collection for combined muon timing information");
  desc.add("MuonDtTimeCollection", edm::InputTag("muons", "dt"))
    ->setComment("A");
  desc.add("MuonCscTimeCollection", edm::InputTag("muons", "csc"))
    ->setComment("A");
  desc.add("MuonDtSegmentCollection", edm::InputTag("dt4DSegments"))
    ->setComment("A");
  desc.add("MuonCscSegmentCollection", edm::InputTag("cscSegments"))
    ->setComment("A");
  desc.add("OfflinePrimaryVerticesCollection", edm::InputTag("offlinePrimaryVertices"))
    ->setComment("A");
  desc.add("LumiScalers", edm::InputTag("scalersRawToDigi"))
    ->setComment("A");
  desc.add("RefittedStandAloneMuonsCollection", edm::InputTag("refittedStandAloneMuons"))
    ->setComment("A");
  desc.add("OfflineBeamSpotCollection", edm::InputTag("offlineBeamSpot"))
    ->setComment("A");
  desc.add("MuonSegmentCollection", edm::InputTag("MuonSegmentProducer"))
    ->setComment("A");
  desc.add("MuonCollection", edm::InputTag("muons"))
    ->setComment("A");
    desc.add("ElectronCollection", edm::InputTag("gedGsfElectrons"))
      ->setComment("A");

      desc.add("conversions", edm::InputTag("allConversions","","RECO"))
        ->setComment("A");


  desc.add("electron_cutbasedID_decisions_veto", edm::InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-veto", ""))
    ->setComment("A");
  desc.add("electron_cutbasedID_decisions_loose", edm::InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-loose", ""))
    ->setComment("A");
  desc.add("electron_cutbasedID_decisions_medium", edm::InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-medium", ""))
    ->setComment("A");
  desc.add("electron_cutbasedID_decisions_tight", edm::InputTag("egmGsfElectronIDs", "cutBasedElectronID-Fall17-94X-V2-tight", ""))
      ->setComment("A");
  desc.add("electron_mvaIsoID_decisions_wp80", edm::InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wp80", ""))
    ->setComment("A");
  desc.add("electron_mvaIsoID_decisions_wp90", edm::InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wp90", ""))
    ->setComment("A");
  desc.add("electron_mvaIsoID_decisions_wpHZZ", edm::InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wpHZZ", ""))
    ->setComment("A");
  desc.add("electron_mvaIsoID_decisions_wpLoose", edm::InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-iso-V2-wpLoose", ""))
    ->setComment("A");
  desc.add("electron_mvaNoIsoID_decisions_wp80", edm::InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-noIso-V2-wp80", ""))
    ->setComment("A");
  desc.add("electron_mvaNoIsoID_decisions_wp90", edm::InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-noIso-V2-wp90", ""))
    ->setComment("A");
  desc.add("electron_mvaNoIsoID_decisions_wpLoose", edm::InputTag("egmGsfElectronIDs", "mvaEleID-Fall17-noIso-V2-wpLoose", ""))
    ->setComment("A");
  desc.add("photon_cutbasedID_decisions_loose", edm::InputTag("egmPhotonIDs", "cutBasedPhotonID-Fall17-94X-V2-loose", ""))
    ->setComment("A");
  desc.add("photon_cutbasedID_decisions_medium", edm::InputTag("egmPhotonIDs", "cutBasedPhotonID-Fall17-94X-V2-medium", ""))
    ->setComment("A");
  desc.add("photon_cutbasedID_decisions_tight", edm::InputTag("egmPhotonIDs", "cutBasedPhotonID-Fall17-94X-V2-tight", ""))
    ->setComment("A");
  desc.add("photon_mvaID_decisions_wp80", edm::InputTag("egmPhotonIDs", "mvaPhoID-RunIIFall17-v2-wp80", ""))
    ->setComment("A");
  desc.add("photon_mvaID_decisions_wp90", edm::InputTag("egmPhotonIDs", "mvaPhoID-RunIIFall17-v2-wp90", ""))
    ->setComment("A");
  desc.add("TriggerResults", edm::InputTag("TriggerResults","","HLT"))
    ->setComment("A");
    desc.add("triggerPrescales",  edm::InputTag("patTrigger"))
      ->setComment("A");



  desc.add<std::string>("FilterName",std::string("hltL3fL1sMu22Or25L1f0L2f10QL3Filtered50Q"))
  ->setComment("Meaning of this is: L1f0  = L1 filtered with threshold 0, L2f10Q = L2 filtered at 10 GeV with quality cuts, L3Filtered50Q = L3 filtered at 50 GeV with quality cuts");
  desc.add("TriggerPathNamesFile", std::string("SUSYBSMAnalysis/Analyzer/data/trigger_names_llp_v3.dat"))
    ->setComment("A");
  desc.add("MuonHLTFilterNamesFile", std::string("SUSYBSMAnalysis/Analyzer/data/MuonHLTFilterNames.dat"))
    ->setComment("A");
  desc.add("PfMET", edm::InputTag("pfMet"))
    ->setComment("A");
  desc.add("PfJet", edm::InputTag("ak4PFJetsCHS"))
    ->setComment("A");
  desc.add("CaloMET", edm::InputTag("caloMet"))
    ->setComment("Take MET from the calorimeters");
  desc.add("CaloJet", edm::InputTag("ak4CaloJets"))
    ->setComment("Take jets from the calorimeters");
  desc.add("TriggerSummary", edm::InputTag("hltTriggerSummaryAOD"))
    ->setComment("A");



  desc.add("PileupInfo", edm::InputTag("addPileupInfo"))
    ->setComment("A");
  desc.add("GenParticleCollection", edm::InputTag("genParticles"))
//  desc.add("GenParticleCollection", edm::InputTag("genParticlesSkimmed"))
    ->setComment("A");
  desc.add("TrackToGenAssoc", edm::InputTag("allTrackMCMatch"))
    ->setComment("Collection used to match to gen thruth");
  desc.add("PfCand", edm::InputTag("particleFlow"))
    ->setComment("Input collection for particleFlow algorithm");
  desc.add("GenCollection", edm::InputTag("generator","","GEN"))
    ->setComment("A");
  desc.addUntracked("TypeMode", 0)
    ->setComment("0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1");
  desc.addUntracked("SampleType", 0)
    ->setComment("0:Data, 1:Background, 2:Signal, 3:Signal Systematics");
  desc.addUntracked<std::string>("SampleName","BaseName")->setComment("This can be used to distinguish different signal models");
  desc.addUntracked<std::string>("Period","2017")->setComment("A");
  desc.addUntracked("SkipSelectionPlot",false)->setComment("A");
  desc.addUntracked("PtHistoUpperBound",4000.0)->setComment("A");
  desc.addUntracked("PHistoUpperBound",10000.0)->setComment("A");
  desc.addUntracked("MassHistoUpperBound",4000.0)->setComment("A");
  desc.addUntracked("MassNBins",400)->setComment("Number of bins in the mass plot");
  desc.addUntracked("IPbound",1.0)
    ->setComment("Number of different Dz side regions used to make cosmic background prediction");
  desc.addUntracked("PredBins",0)
  ->setComment("How many different bins the prediction is split in for analysis being run, sets how many histograms are actually initialized.");
  desc.addUntracked("EtaBins",60)
  ->setComment("How many bins we use for the background prediction method in Eta -- impacts background prediction method -- histograms with the name of the form Pred_Eta in Analysis_PlotStructure.h");
  desc.addUntracked("RegEtaBins",120)->setComment("How many bins we use for the validation of the background estimate method");
  desc.addUntracked("RegIhBins",200)->setComment("How many bins we use for the validation of the background estimate method");
  desc.addUntracked("RegPBins",200)->setComment("How many bins we use for the validation of the background estimate method");
  desc.addUntracked("RegMassBins",50)->setComment("How many bins we use for the validation of the background estimate method");
  desc.addUntracked("DeDxS_UpLim",1.0)->setComment("A");
  desc.addUntracked("DeDxM_UpLim",30.0)->setComment("A");
  desc.addUntracked("DzRegions",6)->setComment("A");
  desc.addUntracked("UseTemplateLayer",false)->setComment("A");
  desc.addUntracked("DeDxSF_0",1.0)->setComment("A");
  desc.addUntracked("DeDxSF_1",1.0325)->setComment("A");
  desc.addUntracked("DeDxK",2.3)->setComment("A");
  desc.addUntracked("DeDxC",3.17)->setComment("A");
  desc.addUntracked("FMIPX",4.0)->setComment("A");
  desc.addUntracked("SaveTree",0)->setComment("0: do not save tree, 6: everything is saved");
  desc.addUntracked("SaveGenTree",0)->setComment("A");
  desc.addUntracked<std::string>("DeDxTemplate","SUSYBSMAnalysis/HSCP/data/template_2017B.root")
    ->setComment("globalIas_ vs Pt templates in eta binning");



  desc.addUntracked<std::string>("TimeOffset","SUSYBSMAnalysis/HSCP/data/MuonTimeOffset.txt")
    ->setComment("MuonTimeOffset info"); // I'm not sure we need this
  desc.add<std::string>("PixelCPE","PixelCPETemplateReco")
    ->setComment("CPE used in the pixel reco, PixelCPEClusterRepair is the best available so far, template only is PixelCPETemplateReco ");
  desc.addUntracked("DebugLevel",0)->setComment("Level of the debugging print statements ");
  desc.addUntracked("HasMCMatch",false)
    ->setComment("Boolean for having the TrackToGenAssoc collection, only new sample have it");
  desc.addUntracked("CalcSystematics",false)->setComment("Boolean to decide whether we want to calculate the systematics");

  // Trigger choice
  // Choice of HLT_Mu50_v is to simplify analysis
  desc.addUntracked("Trigger_Mu", std::vector<std::string>{"HLT_Mu50_v"})
  ->setComment("Add the list of muon triggers");
//  desc.addUntracked("Trigger_MET",  std::vector<std::string>{"HLT_PFMET120_PFMHT120_IDTight_v","HLT_PFHT500_PFMET100_PFMHT100_IDTight_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v","HLT_MET105_IsoTrk50_v"})
    // Possibly used in a next version of the analysis
     desc.addUntracked("Trigger_MET",  std::vector<std::string>{""})
    ->setComment("Add the list of MET triggers");
  // Choice of >55.0 is motivated by the fact that Single muon trigger threshold is 50 GeV
  desc.addUntracked("GlobalMinPt",55.0)->setComment("Cut on pT at PRE-SELECTION");
  // Choice of <1.0 is for detector homogeneity - use only barrel for now - not use disks
  desc.addUntracked("GlobalMaxEta",1.0)->setComment("Cut on inner tracker track eta");
  // Excluding the L1 in BPix because of hardware problems, require >=2 hits for track-probQ
  desc.addUntracked("GlobalMinNOPH",2)->setComment("Cut on number of (valid) track pixel hits");
  // Choice >0.8 is motivated by looking at N1 plot N1_TNOHFraction
  desc.addUntracked("GlobalMinFOVH",0.8)->setComment("Cut on fraction of valid track hits for track cleaning");
  // Choice of >=10 is motivated by N1_TNOM
  desc.addUntracked("GlobalMinNOM",10)->setComment("Cut on number of dEdx hits (#strip+#pixel-#ClusterCleaned-#BPixL1 hits)");
  // Use high purity tracks to ensure good quality tracks
  desc.addUntracked("GlobalUseHighPurity",true)->setComment("Cut on the quality of the track");
  // Choice of <5 is motivated by N1_Chi2oNdof
  desc.addUntracked("GlobalMaxChi2",5.0)->setComment("Cut on Track maximal Chi2/NDF");
  // Choice of <0.1 motivated by looking at N1 plot N1 plot
  desc.addUntracked("GlobalMaxDZ",0.1)->setComment("Cut on 1D distance (cm) to closest vertex in Z direction");
  // Choice of <0.02 motivated by looking at N1 plot N1 plot
  desc.addUntracked("GlobalMaxDXY",0.02)->setComment("Cut on 2D distance (cm) to closest vertex in R direction");

  desc.addUntracked("GlobalMaxMiniRelIsoAll",0.02)->setComment("Cut on the PF based mini-isolation");
  desc.addUntracked("GlobalMaxTIsol",15.0)->setComment("Cut on tracker isolation (SumPt of genTracks with variable cone)");
  desc.addUntracked("GlobalMaxEoP",0.3)->setComment("Cut on calorimeter isolation (E/P) using PF");

  desc.addUntracked("GlobalMinTrackProbQCut",0.0)->setComment("Min cut for probQ, 0.0 means no cuts applied");
  desc.addUntracked("GlobalMaxTrackProbQCut",0.7)->setComment("Max cut for probQ, 1.0 means no cuts applied");
  desc.addUntracked("GlobalMinIh",3.47)->setComment("Cut on dEdx estimator (Im,Ih,etc)");
  desc.addUntracked("GlobalMinTrackProbXYCut",0.01)->setComment("Min cut for probXY, -0.01 means no cuts applied");
  desc.addUntracked("GlobalMaxTrackProbXYCut",1.0)->setComment("Max cut for probXY, 1.0 means no cuts applied");
  desc.addUntracked("GlobalMinIs",0.0)->setComment("Cut on dEdx discriminator (Ias,Is,etc)");
  desc.addUntracked("GlobalMinDeltaRminJet",0.4)->setComment("Min distance in dR to the nearest jet");
  desc.addUntracked("MinMuStations",2)->setComment("Minimum number of muon stations");
  desc.addUntracked("GlobalMinNDOF",8.0)->setComment("Cut on number of DegreeOfFreedom used for muon TOF measurement");
  desc.addUntracked("GlobalMinNDOFDT",6.0)->setComment("Cut on number of DT DegreeOfFreedom used for muon TOF measurement");
  desc.addUntracked("GlobalMinNDOFCSC",6.0)->setComment("Cut on number of CSC DegreeOfFreedom used for muon TOF measurement");
  desc.addUntracked("GlobalMaxTOFErr",0.15)->setComment("Cut on error on muon TOF measurement");
  desc.addUntracked("GlobalMinTOF",1.0)->setComment("Cut on the min time-of-flight");

  //Templates related parameters

  desc.addUntracked("PileUpTreatment",true)->setComment("Boolean to decide whether we want to have pile up dependent templates or not");
  desc.addUntracked("CreateGiTemplates",true)->setComment("Boolean to decide whether we create templates or not, true means we generate");

  desc.addUntracked("CreateAndRunGitemplates",true)->setComment("Boolean to decide if we create PU dependent templates while using other template for analysis");


  desc.addUntracked("NbPileUpBins",5)->setComment("Number of Pile-Up bins for IAS templates");
  desc.addUntracked("PileUpBins",  std::vector<int>{0,20,25,30,35,200})->setComment("Choice of Pile-Up Bins");

 descriptions.add("HSCParticleAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);

//=============================================================
//
//     Method for initializing pT and globalIas_ cuts
//
//=============================================================
void Analyzer::initializeCuts(edm::Service<TFileService>& fs,
                              vector<float>& CutPt,
                              vector<float>& CutI,
                              vector<float>& CutTOF,
                              vector<float>& CutPt_Flip,
                              vector<float>& CutI_Flip,
                              vector<float>& CutTOF_Flip) {
  CutPt.clear();
  CutI.clear();
  CutTOF.clear();
  CutPt_Flip.clear();
  CutI_Flip.clear();
  CutTOF_Flip.clear();

  CutPt.push_back(globalMinPt_);
  CutI.push_back(globalMinIs_);
  CutTOF.push_back(globalMinTOF_);
  CutPt_Flip.push_back(globalMinPt_);
  CutI_Flip.push_back(globalMinIs_);
  CutTOF_Flip.push_back(globalMinTOF_ );

  if (typeMode_ < 2) {
    for (float Pt = globalMinPt_ + 5; Pt < 200; Pt += 5) {
      for (float I = globalMinIs_ + 0.025; I < 0.45; I += 0.025) {
        CutPt.push_back(Pt);
        CutI.push_back(I);
        CutTOF.push_back(-1);
      }
    }
  } else if (typeMode_ == 2) {
    for (float Pt = globalMinPt_ + 5; Pt < 120; Pt += 5) {
      if (Pt > 80 && ((int)Pt) % 10 != 0)
        continue;
      for (float I = globalMinIs_ + 0.025; I < 0.40; I += 0.025) {
        for (float TOF = globalMinTOF_ + 0.025; TOF < 1.35; TOF += 0.025) {
          CutPt.push_back(Pt);
          CutI.push_back(I);
          CutTOF.push_back(TOF);
        }
      }
    }
    for (float Pt = globalMinPt_ + 10; Pt < 90; Pt += 30) {
      for (float I = globalMinIs_ + 0.1; I < 0.30; I += 0.1) {
        for (float TOF = globalMinTOF_ - 0.05; TOF > 0.65; TOF -= 0.05) {
          CutPt_Flip.push_back(Pt);
          CutI_Flip.push_back(I);
          CutTOF_Flip.push_back(TOF);
        }
      }
    }
  } else if (typeMode_ == 3) {
    for (float Pt = globalMinPt_ + 30; Pt < 450; Pt += 30) {
      for (float TOF = globalMinTOF_ + 0.025; TOF < 1.5; TOF += 0.025) {
        CutPt.push_back(Pt);
        CutI.push_back(-1);
        CutTOF.push_back(TOF);
      }
    }
    for (float Pt = globalMinPt_ + 30; Pt < 450; Pt += 60) {
      for (float TOF = globalMinTOF_ - 0.025; TOF > 0.5; TOF -= 0.025) {
        CutPt_Flip.push_back(Pt);
        CutI_Flip.push_back(-1);
        CutTOF_Flip.push_back(TOF);
      }
    }
  } else if (typeMode_ == 4) {
    for (float I = globalMinIs_ + 0.025; I < 0.55; I += 0.025) {
      for (float TOF = globalMinTOF_ + 0.025; TOF < 1.46; TOF += 0.025) {
        CutPt.push_back(-1);
        CutI.push_back(I);
        CutTOF.push_back(TOF);
      }
    }
    for (float I = globalMinIs_ + 0.025; I < 0.55; I += 0.025) {
      for (float TOF = globalMinTOF_ - 0.025; TOF > 0.54; TOF -= 0.025) {
        CutPt_Flip.push_back(-1);
        CutI_Flip.push_back(I);
        CutTOF_Flip.push_back(TOF);
      }
    }
  } else if (typeMode_ == 5) {
    for (float Pt = 75; Pt <= 150; Pt += 25) {
      for (float I = 0.0; I <= 0.45; I += 0.025) {
        CutPt.push_back(Pt);
        CutI.push_back(I);
        CutTOF.push_back(-1);
        CutPt_Flip.push_back(Pt);
        CutI_Flip.push_back(I);
        CutTOF_Flip.push_back(-1);
      }
    }
  }

  //printf("%i Different Final Selection will be tested\n",(int)CutPt.size());
  //printf("%i Different Final Selection will be tested for background uncertainty\n",(int)CutPt_Flip.size());
  edm::LogInfo("Analyzer") << CutPt.size() << " Different Final Selection will be tested\n"
                           << CutPt_Flip.size()
                           << " Different Final Selection will be tested for background uncertainty";

  //Initialization of variables that are common to all samples
  HCuts["Pt"] = fs->make<TProfile>("HCuts_Pt", "HCuts_Pt", CutPt.size(), 0, CutPt.size());
  HCuts["I"] = fs->make<TProfile>("HCuts_I", "HCuts_I", CutPt.size(), 0, CutPt.size());
  HCuts["TOF"] = fs->make<TProfile>("HCuts_TOF", "HCuts_TOF", CutPt.size(), 0, CutPt.size());
  for (unsigned int i = 0; i < CutPt.size(); i++) {
    HCuts["Pt"]->Fill(i, CutPt[i]);
    HCuts["I"]->Fill(i, CutI[i]);
    HCuts["TOF"]->Fill(i, CutTOF[i]);
  }

  HCuts["Pt_Flip"] = fs->make<TProfile>("HCuts_Pt_Flip", "HCuts_Pt_Flip", CutPt_Flip.size(), 0, CutPt_Flip.size());
  HCuts["I_Flip"] = fs->make<TProfile>("HCuts_I_Flip", "HCuts_I_Flip", CutPt_Flip.size(), 0, CutPt_Flip.size());
  HCuts["TOF_Flip"] = fs->make<TProfile>("HCuts_TOF_Flip", "HCuts_TOF_Flip", CutPt_Flip.size(), 0, CutPt_Flip.size());
  for (unsigned int i = 0; i < CutPt_Flip.size(); i++) {
    HCuts["Pt_Flip"]->Fill(i, CutPt_Flip[i]);
    HCuts["I_Flip"]->Fill(i, CutI_Flip[i]);
    HCuts["TOF_Flip"]->Fill(i, CutTOF_Flip[i]);
  }
}

//=============================================================
//     Method for scaling eta Vs bin
//=============================================================
float Analyzer::scaleFactor(float eta) {
  float etaBins[15] = {-2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
  float scaleBins[15] = {0, 0.97, 1.06, 1.00, 0.89, 0.91, 0.93, 0.93, 0.92, 0.92, 0.91, 0.89, 1.00, 1.06, 0.99};
  for (int i = 0; i < 15; i++)
    if (eta < etaBins[i])
      return scaleBins[i];
  return 0;
}

//======================================================================
//     Method for returning eta and beta dependent trigger syst factors
//======================================================================
float Analyzer::triggerSystFactor(float eta, float beta, int syst) {
  float betaBins[17] = {0.06, 0.12, 0.17, 0.23, 0.28, 0.34, 0.39, 0.45, 0.5, 0.56, 0.61, 0.67, 0.72, 0.78, 0.83, 0.89, 1};
  if (syst > 0) {
    // Up systematics
    if (fabs(eta) < 0.3) {
      float scaleBins[17] = {0.0,0.98,1.04,0.98,1.02,0.98,1.01,0.97,1.11,1.05,1.0,0.98,0.99,0.96,1.0,1.05,0.94};
      for (int i = 0; i < 17; i++) {
        if (beta < betaBins[i]) {
          return scaleBins[i];
        }
      }
    } else if(fabs(eta) < 0.6) {
      float scaleBins[17] = {0.0,0.99,1.0,0.99,1.0,1.0,1.01,1.01,1.04,1.06,1.07,0.98,0.97,0.99,1.0,0.96,1.01};
      for (int i = 0; i < 17; i++) {
        if (beta < betaBins[i]) {
          return scaleBins[i];
        }
      }
    } else {
      float scaleBins[17] = {0.0,0.99,0.9,1.01,1.04,1.01,1.02,0.99,0.97,1.13,1.08,1.1,0.96,0.95,0.98,0.98,0.98};
      for (int i = 0; i < 17; i++) {
        if (beta < betaBins[i]) {
          return scaleBins[i];
        }
      }
    }
  } else {
    // Down systematics
    if (fabs(eta) < 0.3) {
      float scaleBins[17] = {0,0.49,0.98,0.97,0.93,1.01,1.02,0.97,0.96,0.73,0.78,0.86,0.93,0.95,0.97,1.06,};
      for (int i = 0; i < 17; i++) {
        if (beta < betaBins[i]) {
          return scaleBins[i];
        }
      }
    } else if(fabs(eta) < 0.6) {
      float scaleBins[17] = {0.0,1.73,0.9,1.18,1.1,0.86,1.07,0.95,0.98,0.79,0.7,0.77,0.88,0.98,1.0,0.99,1.02};
      for (int i = 0; i < 17; i++) {
        if (beta < betaBins[i]) {
          return scaleBins[i];
        }
      }
    } else {
      float scaleBins[17] = {0.0,0.0,0.95,1.1,0.81,1.11,1.14,0.95,0.99,0.95,0.79,0.7,0.73,0.86,0.96,1.0,0.99};
      for (int i = 0; i < 17; i++) {
        if (beta < betaBins[i]) {
          return scaleBins[i];
        }
      }
    }
  }
  return 0;
}


//=============================================================
//     Method for rescaling pT
//=============================================================
// with the new method, no need to pass phi and charge
float Analyzer::shiftForPt(const float& pt, const float& eta, const float& phi, const int& charge) {

// based on https://twiki.cern.ch/twiki/bin/viewauth/CMS/MuonUL2018#Momentum_Resolution
// AN-2018/008
  // factor 0.46 corresponds to the 10% shift uncertainty
  RNG2 = new TRandom3();
  float sigma = 0.0;
  if (fabs(eta) < 1.2) {
    sigma = 0.0141926 + 4.23456e-05*pt - 9.91644e-09*pt*pt;
  } else {
    sigma = 0.016883 + 6.21438e-05*pt - 9.79418e-09*pt*pt;
  }
  float shift = RNG2->Gaus(0,sigma*0.46);
  return shift;
// In the past the following was used:
//    float newInvPt = 1 / pt + 0.000236 - 0.000135 * pow(eta, 2) + charge * 0.000282 * TMath::Sin(phi - 1.337);
//    return 1 / newInvPt;
}

//=============================================================
//
//     Method to get hit position
//
//=============================================================
GlobalPoint Analyzer::getOuterHitPos(const edm::EventSetup& iSetup, const reco::DeDxHitInfo* dedxHits) {
  // Retrieve tracker geometry from the event setup
  edm::ESHandle<TrackerGeometry> tkGeometry;
  iSetup.get<TrackerDigiGeometryRecord>().get(tkGeometry);

  GlobalPoint point(0, 0, 0);
  if (!dedxHits) {
    return point;
  }
  float outerDistance=-1;
  for (unsigned int h = 0; h < dedxHits->size(); h++) {
    const auto detid = dedxHits->detId(h);
    const auto &surface = (*tkGeometry->idToDetUnit(detid)).surface();
    const GlobalPoint globalPoint = surface.toGlobal(Local3DPoint(dedxHits->pos(h).x(), dedxHits->pos(h).y(), dedxHits->pos(h).z()));
    float distanceForIdxH = sqrt(globalPoint.x()*globalPoint.x() + globalPoint.y()*globalPoint.y() + globalPoint.z()*globalPoint.z());
    if(distanceForIdxH>outerDistance) {
      outerDistance=distanceForIdxH;
      point=globalPoint;
    }
  }
  return point;
}

//=============================================================
//
//     Method for muon segment separation
//
//=============================================================
float Analyzer::SegSep(const reco::TrackRef track, const edm::Event& iEvent, float& minPhi, float& minEta) {
  if (typeMode_ != 3)
    return -1;

  if (track.isNull())
    return false;

  float minDr = 10;
  minPhi = 10;
  minEta = 10;

  //Look for segment on opposite side of detector from track
  for (const auto& segment : iEvent.get(muonSegmentToken_)) {
    GlobalPoint gp = segment.getGP();

    //Flip HSCP to opposite side of detector
    float eta_hscp = -1 * track->eta();
    float phi_hscp = track->phi() + M_PI;

    float deta = gp.eta() - eta_hscp;
    float dphi = gp.phi() - phi_hscp;
    while (dphi > M_PI)
      dphi -= 2 * M_PI;
    while (dphi <= -M_PI)
      dphi += 2 * M_PI;

    //Find segment most opposite in eta
    //Require phi difference of 0.5 so it doesn't match to own segment
    if (fabs(deta) < fabs(minEta) && fabs(dphi) < (M_PI - 0.5)) {
      minEta = deta;
    }
    //Find segment most opposite in phi
    if (fabs(dphi) < fabs(minPhi)) {
      minPhi = dphi;
    }
    //Find segment most opposite in Eta-Phi
    float dR = sqrt(deta * deta + dphi * dphi);
    if (dR < minDr) {
      minDr = dR;
    }
  }
  return minDr;
}

//=============================================================
//
//     Pre-Selection
//
//=============================================================
template <typename T, size_t n>
bool Analyzer::passPreselection(T (&passedCutsArray)[n], bool verbose) {
  using namespace edm;
  std::map<int, std::string> namingMap;
  namingMap[0] = "Trigger";
  namingMap[1] = "p_{T}";
  namingMap[2] = "#eta";
  namingMap[3] = "N_{no-L1 pixel hits}";
  namingMap[4] = "f_{valid/all hits}";
  namingMap[5] = "N_{dEdx hits}";
  namingMap[6] = "HighPurity";
  namingMap[7] = "#chi^{2} / N_{dof}";
  namingMap[8] = "d_{z}";
  namingMap[9] = "d_{xy}";
  namingMap[10] = "MiniRelIsoAll";
  namingMap[11] = "MiniRelTkIso";
  namingMap[12] = "E/p";
  namingMap[13] = "#sigma_{p_{T}} / p_{T}^{2}";
  namingMap[14] = "F_{i}";
  
  // Return false in the function if a given cut is not passed
  for (size_t i=0;i<sizeof(T) * n;i++) {
    if (passedCutsArray[i]) {
    } else {
      if (debug_ > 2 && trigInfo_ > 0 && verbose) LogPrint(MOD) << "        >> Preselection not passed for the " <<  namingMap[i];
      return false;
    }
  }


/*
 // will put this back with the new way of doing preselection


  if (tof) {
    if ((typeMode_ > 1 && typeMode_ != 5) && tof->inverseBetaErr() > GlobalMaxTOFErr) {
      return false;
    }

    if (typeMode_ == 3 && min(min(fabs(tof->timeAtIpInOut() - 100), fabs(tof->timeAtIpInOut() - 50)),
                              min(fabs(tof->timeAtIpInOut() + 100), fabs(tof->timeAtIpInOut() + 50))) < 5)
      return false;
  } // End condition on tof existence or not
*/
//  if (cutEtaTOFOnly) {
//    if (debug_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: for TOF only analysis, eta is too low";
//    return false;
//  } else if (typeMode_ == 3 && fabs(minEta) > minSegEtaSep) {
//    if (debug_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for TOF eta cut";
//  }
//  if (tuple)
//    tuple->BefPreS_Phi->Fill(track->phi(), EventWeight_);

//  if (cutPhiTOFOnly) {
//    if (debug_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: for TOF only analysis, 1.2 < phi < 1.9";
//    return false;
//  }

//    float RecoQoPt = track->charge() / track->pt();
//    if (!hscp.trackRef().isNull() && hscp.trackRef()->pt() > 200) {
//      float InnerRecoQoPt = hscp.trackRef()->charge() / hscp.trackRef()->pt();
//      tuple->BefPreS_InnerInvPtDiff->Fill((RecoQoPt - InnerRecoQoPt) / InnerRecoQoPt, EventWeight_);
//    }

  return true;
}
//=============================================================
//
//     Selection
//
//=============================================================
// TAV: iEvent is not needed to be passed
bool Analyzer::passSelection(const reco::TrackRef track,
                             const reco::DeDxData* dedxSObj,
                             const reco::DeDxData* dedxMObj,
                             const reco::MuonTimeExtra* tof,
                             const edm::Event& iEvent,
                             const int& CutIndex,
                             Tuple* tuple,
                             const bool isFlip,
                             const float genBeta,
                             const bool RescaleP,
                             const float RescaleI,
                             const float RescaleT) {
  static constexpr const char* const MOD = "Analyzer";
  using namespace edm;
  
  if (track.isNull()) {
    // Should be a LogError I believe...
    LogPrint(MOD) << "@passSelection: track.isNull() -- this should never happen!!!";
    return false;
  }

  float MuonTOF = (tof) ? tof->inverseBeta() : globalMinTOF_ ;
  float PtCut = CutPt_[CutIndex];
  float IasCut = CutI_[CutIndex];
  float TOFCut = -1.;
  if (typeMode_ > 1) TOFCut = CutTOF_[CutIndex];
  if (isFlip) {
    PtCut = CutPt_Flip_[CutIndex];
    IasCut = CutI_Flip_[CutIndex];
    if (typeMode_ > 1) TOFCut = CutTOF_Flip_[CutIndex];
  }

  // Check if we pass the momentum selection
  if (RescaleP) {
    if (track->pt()*(1+shiftForPt(track->pt(),track->eta(),track->phi(),track->charge())) < PtCut)
      return false;
  } else if (track->pt() < PtCut) {
      return false;
  }
  
  // Distribtution of genBeta after Pt selection is passed
  if (tuple && genBeta >= 0) {
    tuple->PostS_CutIdVsBeta_postPt->Fill(CutIndex, genBeta, EventWeight_);
  }

  // Check if we pass the (rescalled) Ias selection
  if (typeMode_ != 3 && globalIas_ + RescaleI < IasCut) {
    return false;
  }

  // Distribtution of genBeta after Pt and Ias selection is passed
  if (tuple && genBeta >= 0) {
      tuple->PostS_CutIdVsBeta_postPtAndIas->Fill(CutIndex, genBeta, EventWeight_);
  }

  // Check if we pass the TOF selection
  if ((typeMode_ > 1 && typeMode_ != 5) && !isFlip && MuonTOF + RescaleT < TOFCut) {
    return false;
  }
  if ((typeMode_ > 1 && typeMode_ != 5) && isFlip && MuonTOF + RescaleT > TOFCut) {
    return false;
  }

  if (tuple) {
    if (genBeta >= 0) {
      tuple->PostS_CutIdVsBeta_postPtAndIasAndTOF->Fill(CutIndex, genBeta, EventWeight_);
    }
    tuple->PostS_CutIdVsP->Fill(CutIndex, track->p(), EventWeight_);
    tuple->PostS_CutIdVsPt->Fill(CutIndex, track->pt(), EventWeight_);
    tuple->PostS_CutIdVsIas->Fill(CutIndex, globalIas_, EventWeight_);
    tuple->PostS_CutIdVsIh->Fill(CutIndex, globalIh_, EventWeight_);

    tuple->PostS_CutIdVsPVsIas->Fill(CutIndex, track->p(), globalIas_, EventWeight_);
    tuple->PostS_CutIdVsPVsIh->Fill(CutIndex, track->p(), globalIh_, EventWeight_);
    tuple->PostS_CutIdVsPtVsIas->Fill(CutIndex, track->pt(), globalIas_, EventWeight_);
    tuple->PostS_CutIdVsPtVsIh->Fill(CutIndex, track->pt(), globalIh_, EventWeight_);
    if (typeMode_ > 1) {
      tuple->PostS_CutIdVsTOF->Fill(CutIndex, MuonTOF, EventWeight_);
      tuple->PostS_CutIdVsTOFVsIas->Fill(CutIndex, MuonTOF, globalIas_, EventWeight_);
      tuple->PostS_CutIdVsTOFVsIh->Fill(CutIndex, MuonTOF, globalIh_, EventWeight_);
    }
  }
  return true;
}

//=============================================================
//
//     Combine individual probs into a track level one
//
//=============================================================
float Analyzer::combineProbs(float probOnTrackWMulti, int numRecHits) const {
  float logprobOnTrackWMulti = (probOnTrackWMulti > 0) ? log(probOnTrackWMulti) : 0;
  float factQ = -logprobOnTrackWMulti;
  float probOnTrackTerm = 0.f;

  if (numRecHits == 1) {
    probOnTrackTerm = 1.f;
  } else if (numRecHits > 1) {
    probOnTrackTerm = 1.f + factQ;
    for (int iTkRh = 2; iTkRh < numRecHits; ++iTkRh) {
      factQ *= -logprobOnTrackWMulti / float(iTkRh);
      probOnTrackTerm += factQ;
    }
  }
  float probOnTrack = probOnTrackWMulti * probOnTrackTerm;

  return probOnTrack;
}


//=============================================================
//     Check if the GenIDs are for a HSCP
//=============================================================
bool Analyzer::isHSCPgenID(const reco::GenParticle& gen) {
  int thePDGidForCandidate = abs(gen.pdgId());
  // R-hadrons
  if (   thePDGidForCandidate == 1000993 || thePDGidForCandidate == 1009113
      || thePDGidForCandidate == 1009223 || thePDGidForCandidate == 1009313
      || thePDGidForCandidate == 1009333 || thePDGidForCandidate == 1092114
      || thePDGidForCandidate == 1093214 || thePDGidForCandidate == 1093324
      || thePDGidForCandidate == 1000622 || thePDGidForCandidate == 1000642
      || thePDGidForCandidate == 1006113 || thePDGidForCandidate == 1006311
      || thePDGidForCandidate == 1006313 || thePDGidForCandidate == 1006333) {
    return true;
  }
  // Single-charged HSCP
  else if (   thePDGidForCandidate == 1009213 || thePDGidForCandidate == 1009323
       || thePDGidForCandidate == 1091114 || thePDGidForCandidate == 1092214
       || thePDGidForCandidate == 1093114 || thePDGidForCandidate == 1093224
       || thePDGidForCandidate == 1093314 || thePDGidForCandidate == 1093334
       || thePDGidForCandidate == 1000612 || thePDGidForCandidate == 1000632
       || thePDGidForCandidate == 1000652 || thePDGidForCandidate == 1006211
       || thePDGidForCandidate == 1006213 || thePDGidForCandidate == 1006321
       || thePDGidForCandidate == 1006323 || thePDGidForCandidate == 1000015) {
    return true;
  }
  // Double-charged R-hadrons
  else if (thePDGidForCandidate == 1092224 || thePDGidForCandidate == 1006223) {
    return true;
  }
  // tau prime, could be single or multiple charged
  else if (thePDGidForCandidate == 17) {
    return true;
  } else {
    return false;
  }
}

//=============================================================
//
//     Calculate systematics on signal
//
//=============================================================
void Analyzer::calculateSyst(reco::TrackRef track,
                             const reco::DeDxHitInfo* dedxHits,
                             const reco::DeDxData* dedxSObj,
                             const reco::DeDxData* dedxMObj,
                             const reco::MuonTimeExtra* tof,
                             const edm::Event& iEvent,
                             const edm::EventSetup& iSetup,
                             const float pixelProbs[],
                             Tuple* tuple,
                             const float genBeta,
                             float MassErr,
                             const float closestBackgroundPDGsIDs[]) {
  // Will put this back with the new way of preselections
//  //FIXME to be measured on 2015 data, currently assume 2012
//  bool PRescale = true;

//  float MRescale = 0.95;
//  float TRescale = -0.015;  //-0.005 (used in 2012); // added to the 1/beta value
//
//  // compute systematic due to momentum scale
//  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, iSetup, pixelProbs, nullptr, -1, PRescale, 0, 0, 0, closestBackgroundPDGsIDs)) {
//    float RescalingFactor = RescaledPt(track->pt(), track->eta(), track->phi(), track->charge()) / track->pt();
//
//    float Mass = -1;
//    if (dedxMObj)
//      Mass = GetMass(track->p() * RescalingFactor, dedxMObj->dEdx(), dEdxK_, dEdxC_);
//    float MassTOF = -1;
//    if (tof)
//      MassTOF = GetTOFMass(track->p() * RescalingFactor, tof->inverseBeta());
//    float MassComb = -1;
//    if (tof && dedxMObj)
//      MassComb = GetMassFromBeta(track->p() * RescalingFactor,
//                                 (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
//    else if (dedxMObj)
//      MassComb = Mass;
//    if (tof)
//      MassComb = MassTOF;
//
//    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
//      if (passSelection(track,
//                        dedxSObj,
//                        dedxMObj,
//                        tof,
//                        iEvent,
//                        CutIndex,
//                        nullptr,
//                        false,
//                        -1,
//                        PRescale,
//                        0,
//                        0)) {  //WAIT//
//        HSCPTk_SystP[CutIndex] = true;
//        if (Mass > MaxMass_SystP[CutIndex])
//          MaxMass_SystP[CutIndex] = Mass;
//        tuple->Mass_SystP->Fill(CutIndex, Mass, EventWeight_);
//        if (tof) {
//          tuple->MassTOF_SystP->Fill(CutIndex, MassTOF, EventWeight_);
//        }
//        tuple->MassComb_SystP->Fill(CutIndex, MassComb, EventWeight_);
//      }
//    } // end loop on cut index
//  } // end compute systematic due to momentum scale
//  // compute systematic due to dEdx (both globalIas_ and Ih)
//  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, iSetup, pixelProbs, nullptr, -1, false, 0.0, IRescale, 0.0, closestBackgroundPDGsIDs)) {
//      //if(TypeMode==5 && isSemiCosmicSB)continue;
//    float Mass = -1;
//    if (dedxMObj)
//      Mass = GetMass(track->p(), dedxMObj->dEdx() * MRescale, dEdxK_, dEdxC_);
//    float MassTOF = -1;
//    if (tof)
//      MassTOF = GetTOFMass(track->p(), tof->inverseBeta());
//    float MassComb = -1;
//    if (tof && dedxMObj)
//      MassComb =
//      GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
//    else if (dedxMObj)
//      MassComb = Mass;
//    if (tof)
//      MassComb = MassTOF;
//    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
//      if (passSelection(track, dedxSObj, dedxMObj, tof, iEvent, CutIndex, nullptr, false, -1, 0, IRescale, 0)) {
//        HSCPTk_SystI[CutIndex] = true;
//        if (Mass > MaxMass_SystI[CutIndex])
//          MaxMass_SystI[CutIndex] = Mass;
//        tuple->Mass_SystI->Fill(CutIndex, Mass, EventWeight_);
//        if (tof)
//          tuple->MassTOF_SystI->Fill(CutIndex, MassTOF, EventWeight_);
//        tuple->MassComb_SystI->Fill(CutIndex, MassComb, EventWeight_);
//      }
//    }
//  } // End compute systematic due to dEdx
//  // compute systematic due to Mass shift ??????????
//  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, iSetup, pixelProbs, nullptr, -1, 0, 0, 0, 0, closestBackgroundPDGsIDs)) {
//    /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
//    float Mass = -1;
//    if (dedxMObj)
//      Mass = GetMass(track->p(), dedxMObj->dEdx() * MRescale, dEdxK_, dEdxC_);
//    float MassTOF = -1;
//    if (tof)
//      MassTOF = GetTOFMass(track->p(), tof->inverseBeta());
//    float MassComb = -1;
//    if (tof && dedxMObj)
//      MassComb = GetMassFromBeta(
//                                 track->p(), (GetIBeta(dedxMObj->dEdx() * MRescale, dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
//    else if (dedxMObj)
//      MassComb = Mass;
//    if (tof)
//      MassComb = MassTOF;
//
//    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
//      if (passSelection(track, dedxSObj, dedxMObj, tof, iEvent, CutIndex, nullptr, false, -1, 0, 0, 0)) {
//        HSCPTk_SystM[CutIndex] = true;
//        if (Mass > MaxMass_SystM[CutIndex])
//          MaxMass_SystM[CutIndex] = Mass;
//        tuple->Mass_SystM->Fill(CutIndex, Mass, EventWeight_);
//        if (tof)
//          tuple->MassTOF_SystM->Fill(CutIndex, MassTOF, EventWeight_);
//        tuple->MassComb_SystM->Fill(CutIndex, MassComb, EventWeight_);
//      }
//    }
//  } // End compute systematic due to Mass shift
//  // compute systematic due to TOF
//  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, iSetup, pixelProbs, nullptr, -1, 0, 0, TRescale, 0, closestBackgroundPDGsIDs)) {
//    /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
//    float Mass = -1;
//    if (dedxMObj)
//      Mass = GetMass(track->p(), dedxMObj->dEdx(), dEdxK_, dEdxC_);
//    float MassTOF = -1;
//    if (tof)
//      MassTOF = GetTOFMass(track->p(), (tof->inverseBeta() + TRescale));
//    float MassComb = -1;
//    if (tof && dedxMObj)
//      MassComb = GetMassFromBeta(
//                                 track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / (tof->inverseBeta() + TRescale))) * 0.5);
//    else if (dedxMObj)
//      MassComb = Mass;
//    if (tof)
//      MassComb = MassTOF;
//
//    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
//      if (passSelection(track, dedxSObj, dedxMObj, tof, iEvent, CutIndex, nullptr, false, -1, 0, 0, TRescale)) {
//        HSCPTk_SystT[CutIndex] = true;
//        if (Mass > MaxMass_SystT[CutIndex])
//          MaxMass_SystT[CutIndex] = Mass;
//        tuple->Mass_SystT->Fill(CutIndex, Mass, EventWeight_);
//        if (tof)
//          tuple->MassTOF_SystT->Fill(CutIndex, MassTOF, EventWeight_);
//        tuple->MassComb_SystT->Fill(CutIndex, MassComb, EventWeight_);
//      }
//    }
//  } // End condition for compute systematic due to TOF
//  // compute systematics due to PU
//  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, iSetup, pixelProbs, nullptr, -1, 0, 0, 0, 0, closestBackgroundPDGsIDs)) {
//    /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
//    float Mass = -1;
//    if (dedxMObj)
//      Mass = GetMass(track->p(), dedxMObj->dEdx(), dEdxK_, dEdxC_);
//    float MassTOF = -1;
//    if (tof)
//      MassTOF = GetTOFMass(track->p(), tof->inverseBeta());
//    float MassComb = -1;
//    if (tof && dedxMObj)
//      MassComb =
//      GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
//    else if (dedxMObj)
//      MassComb = Mass;
//    if (tof)
//      MassComb = MassTOF;
//
//    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
//      if (passSelection(track, dedxSObj, dedxMObj, tof, iEvent, CutIndex, nullptr, false, -1, 0, 0, 0)) {
//        HSCPTk_SystPU[CutIndex] = true;
//        if (Mass > MaxMass_SystPU[CutIndex])
//          MaxMass_SystPU[CutIndex] = Mass;
//        tuple->Mass_SystPU->Fill(CutIndex, Mass, EventWeight_ * PUSystFactor_[0]);
//        if (tof)
//          tuple->MassTOF_SystPU->Fill(CutIndex, MassTOF, EventWeight_ * PUSystFactor_[0]);
//        tuple->MassComb_SystPU->Fill(CutIndex, MassComb, EventWeight_ * PUSystFactor_[0]);
//      }
//    }
//  }  // End compute systematics due to PU
}


const reco::Candidate* Analyzer::findFirstMotherWithDifferentID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is this the first parent with a different ID? If yes, return, otherwise
  // go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->pdgId() == particle->mother(0)->pdgId()
	&& particle->mother(0)->status() != 11  // prevent infinite loop for sherpa documentation gluons
	) {
      return findFirstMotherWithDifferentID(particle->mother(0));
    } else {
      return particle->mother(0);
    }
  }

  return 0;
};

const reco::Candidate* Analyzer::findOriginalMotherWithSameID(const reco::Candidate *particle){

  if( particle == 0 ){
    printf("ERROR! null candidate pointer, this should never happen\n");
    return 0;
  }

  // Is there another parent with the same ID? If yes, go deeper into recursion
  if (particle->numberOfMothers() > 0 && particle->pdgId() != 0) {
    if (particle->mother(0)->numberOfMothers() == 0 ||
	particle->mother(0)->status() == 11 ||  // prevent infinite loop for sherpa documentation gluons
	(particle->mother(0)->numberOfMothers() > 0 && particle->mother(0)->mother(0)->pdgId() != particle->mother(0)->pdgId())
	) {
      return particle->mother(0);
    } else {
      return findOriginalMotherWithSameID(particle->mother(0));
    }
  }

  return 0;
}
