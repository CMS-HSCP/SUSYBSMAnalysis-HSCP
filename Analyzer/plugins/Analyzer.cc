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

// v19p0
// - change double to float
// - create fillDescription
// - intro ptErrOverPt vs ptErrOverPt2
// - change the order of preselection cuts
// - N-1 plots
// - Add two more cutflow histos, change boundary for ptErrOverPt2
// - Fix logic for new cutflow, fix the  change boundary for ptErrOverPt2
// - Make cuts into an array
// - Fix logic with not used variales
// - Change the cut flow order
// - Add Ih vs Is plot in preselection, change boundary for dxy/dz plots
// - Change dxy/dz cut default
// - Add plots for MiniIsol, MET, mT
// - Change MiniIsol definition, and plot range, move it to preselection
// - Change EoP to 0.8, then to 2.0 (essentially no cut)
// - Change to allTrackMCMatch
// - 18p3: PF matching to gentracks, change the binning of MiniIso histo
// - 18p4: fix for cutflowProbQfirst index, get rid of EoP cut
// - 18p5 change to new templates
// - 18p5: remove TK iso
// - 18p8: Add postPreselection plots
// - 19p0: One try with TOF
// - 19p1: Change mass binning, remove massT cut
// - 19p3: Simplify probQ cut, change mini-iso def
// - 19p4: Change mini-iso binning
// - 19p5: use charged iso in cutflow, dont cut away out of bound probs, only in preselection
// - 19p6: intro CutFlowEta and PerGenID
// - 19p7: intro NumEvents and HSCPCandidateType, for comparrison, put back EoP cut and TkIso cut (will remove in 19p8)

#include "SUSYBSMAnalysis/Analyzer/plugins/Analyzer.h"

Analyzer::Analyzer(const edm::ParameterSet& iConfig)
    // Read config file
    : hscpToken_(consumes<vector<susybsm::HSCParticle>>(iConfig.getParameter<edm::InputTag>("HscpCollection"))),
      hscpIsoToken_(
          consumes<edm::ValueMap<susybsm::HSCPIsolation>>(iConfig.getParameter<edm::InputTag>("HscpIsoCollection"))),
      muonSegmentToken_(
          consumes<susybsm::MuonSegmentCollection>(iConfig.getParameter<edm::InputTag>("MuonSegmentCollection"))),
      dedxToken_(consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("DedxCollection"))),
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
      triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("TriggerResults"))),
      pfMETToken_(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("PfMET"))),
      pfJetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("PfJet"))),
      CaloMETToken_(consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("CaloMET"))),
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
      massHistoUpperBound_(iConfig.getUntrackedParameter<double>("MassHistoUpperBound")),
      massNBins_(iConfig.getUntrackedParameter<int>("MassNBins")),
      cutOnIPbound_(iConfig.getUntrackedParameter<double>("IPbound")),
      predBins_(iConfig.getUntrackedParameter<int>("PredBins")),
      etaBins_(iConfig.getUntrackedParameter<int>("EtaBins")),
      dEdxS_UpLim_(iConfig.getUntrackedParameter<double>("DeDxS_UpLim")),
      dEdxM_UpLim_(iConfig.getUntrackedParameter<double>("DeDxM_UpLim")),
      numDzRegions_(iConfig.getUntrackedParameter<int>("DzRegions")),
      globalMaxEta_(iConfig.getUntrackedParameter<double>("GlobalMaxEta")),
      globalMinPt_(iConfig.getUntrackedParameter<double>("GlobalMinPt")),
      globalMinNOH_(iConfig.getUntrackedParameter<int>("GlobalMinNOH")),
      globalMinNOPH_(iConfig.getUntrackedParameter<int>("GlobalMinNOPH")),
      globalMinFOVH_(iConfig.getUntrackedParameter<double>("GlobalMinFOVH")),
      globalMinNOM_(iConfig.getUntrackedParameter<int>("GlobalMinNOM")),
      globalMaxChi2_(iConfig.getUntrackedParameter<double>("GlobalMaxChi2")),
      globalMaxEIsol_(iConfig.getUntrackedParameter<double>("GlobalMaxEIsol")),
      globalMaxDZ_(iConfig.getUntrackedParameter<double>("GlobalMaxDZ")),
      globalMaxDXY_(iConfig.getUntrackedParameter<double>("GlobalMaxDXY")),
      globalMaxPtErr_(iConfig.getUntrackedParameter<double>("GlobalMaxPtErr")),
      globalMaxTIsol_(iConfig.getUntrackedParameter<double>("GlobalMaxTIsol")),
      globalMiniRelIsoChg_(iConfig.getUntrackedParameter<double>("GlobalMiniRelIsoChg")),
      globalMassT_(iConfig.getUntrackedParameter<double>("GlobalMassT")),
      globalMinIh_(iConfig.getUntrackedParameter<double>("GlobalMinIh")),
      trackProbQCut_(iConfig.getUntrackedParameter<double>("TrackProbQCut")),
      minMuStations_(iConfig.getUntrackedParameter<int>("MinMuStations")),
      globalMinIs_(iConfig.getUntrackedParameter<double>("GlobalMinIs")),
      globalMinTOF_(iConfig.getUntrackedParameter<double>("GlobalMinTOF")),
      skipPixel_(iConfig.getUntrackedParameter<bool>("SkipPixel")),
      useTemplateLayer_(iConfig.getUntrackedParameter<bool>("UseTemplateLayer")),
      dEdxSF_0_(iConfig.getUntrackedParameter<double>("DeDxSF_0")),
      dEdxSF_1_(iConfig.getUntrackedParameter<double>("DeDxSF_1")),
      dEdxK_(iConfig.getUntrackedParameter<double>("DeDxK")),
      dEdxC_(iConfig.getUntrackedParameter<double>("DeDxC")),
      dEdxTemplate_(iConfig.getUntrackedParameter<string>("DeDxTemplate")),
      enableDeDxCalibration_(iConfig.getUntrackedParameter<bool>("EnableDeDxCalibration")),
      dEdxCalibration_(iConfig.getUntrackedParameter<string>("DeDxCalibration")),
      geometry_(iConfig.getUntrackedParameter<string>("Geometry")),
      timeOffset_(iConfig.getUntrackedParameter<string>("TimeOffset")),
      theFMIPX_(iConfig.getUntrackedParameter<double>("FMIPX")),
      saveTree_(iConfig.getUntrackedParameter<int>("SaveTree")),
      saveGenTree_(iConfig.getUntrackedParameter<int>("SaveGenTree")),
      pixelCPE_(iConfig.getParameter<std::string>("PixelCPE")),
      debug_(iConfig.getUntrackedParameter<int>("DebugLevel")),
      hasMCMatch_(iConfig.getUntrackedParameter<bool>("HasMCMatch")),
      doTriggering_(iConfig.getUntrackedParameter<bool>("DoTriggering")),
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
  dEdxTemplates = loadDeDxTemplate(dEdxTemplate_, splitByModuleType);
  if (enableDeDxCalibration_)
    trackerCorrector.LoadDeDxCalibration(dEdxCalibration_);
  else
    trackerCorrector.TrackerGains = nullptr;

  moduleGeom::loadGeometry(geometry_);
  tofCalculator.loadTimeOffset(timeOffset_);
}

Analyzer::~Analyzer() {
  // do anything here that needs to be done at desctruction time
  // (e.g. close files, deallocate resources etc.)
}

// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob() {
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
                               skipSelectionPlot_,
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
  
  cout << "After initializeTuple" << endl;

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

  TrigInfo_ = 0;

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
}

// ------------ method called for each event  ------------
void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  static constexpr const char* const MOD = "Analyzer";
  using namespace edm;

  //if run change, update conditions
  if (CurrentRun_ != iEvent.id().run()) {
    CurrentRun_ = iEvent.id().run();
    tofCalculator.setRun(CurrentRun_);
    trackerCorrector.setRun(CurrentRun_);

    /* FIXME is it still relevant to use this function ?
    loadDeDxParameters(CurrentRun_, sampleType_, dEdxSF_0_, dEdxSF_1_, dEdxK_, dEdxC_);
    */
    dEdxSF[0] = dEdxSF_0_;
    dEdxSF[1] = dEdxSF_1_;
    //LogInfo(MOD) <<"------> dEdx parameters SF for Run "<<CurrentRun_<< ": "<< dEdxSF[1];
  }
    

  // Compute event weight
  if (!isData) {
    float PUWeight = mcWeight->getEventPUWeight(iEvent, pileupInfoToken_, PUSystFactor_);
    EventWeight_ = PUWeight;  // 1. : unweighted w.r.t pileup
  } else {
    EventWeight_ = 1.;
  }

  if (debug_ > 0 ) LogPrint(MOD) << "\nEvent weight factor applied: " << EventWeight_;

  float HSCPGenBeta1 = -1, HSCPGenBeta2 = -1;

  //get generator weight and pthat
  if (!isData){
      const edm::Handle<GenEventInfoProduct> genEvt = iEvent.getHandle(genEventToken_);
      if (genEvt.isValid()){
          GeneratorWeight_ = genEvt->weight();
          if(genEvt->binningValues().size()>0) GeneratorBinningValues_ = genEvt->binningValues()[0];
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

    GetGenHSCPBeta(genColl, HSCPGenBeta1, HSCPGenBeta2, false);
    if (HSCPGenBeta1 >= 0)
      tuple->Beta_Gen->Fill(HSCPGenBeta1, SignalEventWeight);
    if (HSCPGenBeta2 >= 0)
      tuple->Beta_Gen->Fill(HSCPGenBeta2, SignalEventWeight);

    GetGenHSCPBeta(genColl, HSCPGenBeta1, HSCPGenBeta2, true);
    if (HSCPGenBeta1 >= 0)
      tuple->Beta_GenCharged->Fill(HSCPGenBeta1, SignalEventWeight);
    if (HSCPGenBeta2 >= 0)
      tuple->Beta_GenCharged->Fill(HSCPGenBeta2, SignalEventWeight);

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
      if (isSignal && isGoodGenHSCP(gen,false)) {
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
        tuple->genlevelpT->Fill(gen.pt(), SignalEventWeight);
        tuple->genleveleta->Fill(gen.eta(), SignalEventWeight);
        tuple->genlevelbeta->Fill(gen.p() / gen.energy(), SignalEventWeight);
      
        // Variables for the tuple gen tree branch
        genid.push_back(gen.pdgId());
        gencharge.push_back(gen.charge());
        genmass.push_back(gen.mass());
        genpt.push_back(gen.pt());
        geneta.push_back(gen.eta());
        genphi.push_back(gen.phi());
      } else if (isBckg) {
        // Fill up pT, eta, and beta plots for gen-level background particles
        tuple->genlevelpT->Fill(gen.pt(), EventWeight_);
        tuple->genleveleta->Fill(gen.eta(), EventWeight_);
        tuple->genlevelbeta->Fill(gen.p() / gen.energy(), EventWeight_);
        // TODO: I'm not sure if this needs to be weighted

        // Variables for the tuple gen tree branch
        genid.push_back(gen.pdgId());
        gencharge.push_back(gen.charge());
        genmass.push_back(gen.mass());
        genpt.push_back(gen.pt());
        geneta.push_back(gen.eta());
        genphi.push_back(gen.phi());
      }
    }

    if (debug_ > 3 ) LogPrint(MOD) << "Fill GenTree with basics gen info";
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

  // Get trigger results for this event
  const edm::Handle<edm::TriggerResults> triggerH = iEvent.getHandle(triggerResultsToken_);
  const auto triggerNames = iEvent.triggerNames(*triggerH);

  // These are used in the tree alone, otherwise we use passTriggerPatterns to check the triggers
  bool HLT_Mu50 = false;
  bool HLT_PFMET120_PFMHT120_IDTight = false;
  bool HLT_PFHT500_PFMET100_PFMHT100_IDTight = false;
  bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = false;
  bool HLT_MET105_IsoTrk50 = false;

  for (unsigned int i = 0; i < triggerH->size(); i++) {
    if (TString(triggerNames.triggerName(i)).Contains("HLT_Mu50") && triggerH->accept(i))
      HLT_Mu50 = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_PFMET120_PFMHT120_IDTight") && triggerH->accept(i))
      HLT_PFMET120_PFMHT120_IDTight = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_PFHT500_PFMET100_PFMHT100_IDTight") && triggerH->accept(i))
      HLT_PFHT500_PFMET100_PFMHT100_IDTight = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60") &&
        triggerH->accept(i))
      HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60 = true;
    if (TString(triggerNames.triggerName(i)).Contains("HLT_MET105_IsoTrk50") && triggerH->accept(i))
      HLT_MET105_IsoTrk50 = true;
  }
  // Number of (re-weighted) events
  tuple->NumEvents->Fill(0.5, EventWeight_);
  // Number of (re-weighted with PU syst fact) events
  tuple->NumEvents->Fill(1.5, EventWeight_ * PUSystFactor_[0]);

  // Check if the event is passing trigger
  if (debug_ > 0) LogPrint(MOD) << "Checking if the event is passing trigger...";
  bool metTrig = passTriggerPatterns(triggerH, triggerNames, trigger_met_);
  bool muTrig = passTriggerPatterns(triggerH, triggerNames, trigger_mu_);

  if (!metTrig && muTrig) {
    // mu only
    TrigInfo_ = 1;
  } else if (metTrig && !muTrig) {
    // met only
    TrigInfo_ = 2;  // met only
  } else if (metTrig && muTrig) {
    // mu and met
    TrigInfo_ = 3;
  }

  // If triggering is intended (might not be for some studies and one of the triggers is passing let's analyze the event
  if (doTriggering_ && TrigInfo_ > 0) {
      if (debug_ > 0 ) LogPrint(MOD) << "This event passeed the needed triggers! TrigInfo_ = " << TrigInfo_;
  } else {
      if (debug_ > 0 ) LogPrint(MOD) << "This event did not pass the needed triggers, skipping it";
      return;
     //For TOF only analysis if the event doesn't pass the signal triggers check if it was triggered by the no BPTX cosmic trigger
  }

  // Number of events that pass the trigger
  tuple->NumEvents->Fill(2.5, EventWeight_);

  //keep beta distribution for signal after the trigger
  if (isSignal) {
    if (HSCPGenBeta1 >= 0)
      tuple->Beta_Triggered->Fill(HSCPGenBeta1, EventWeight_);
    if (HSCPGenBeta2 >= 0)
      tuple->Beta_Triggered->Fill(HSCPGenBeta2, EventWeight_);
  }

  // Define handles for DeDx Hits, Muon TOF Combined, Muon TOF DT, Muon TOF CSC
  const edm::Handle<reco::DeDxHitInfoAss> dedxCollH = iEvent.getHandle(dedxToken_);
  const edm::Handle<reco::MuonTimeExtraMap> tofMap = iEvent.getHandle(muonTimeToken_);
  const edm::Handle<reco::MuonTimeExtraMap> tofDtMap = iEvent.getHandle(muonDtTimeToken_);
  const edm::Handle<reco::MuonTimeExtraMap> tofCscMap = iEvent.getHandle(muonCscTimeToken_);

  //================= Handle For Muon DT/CSC Segment ===============
  edm::Handle<CSCSegmentCollection> CSCSegmentCollH;
  edm::Handle<DTRecSegment4DCollection> DTSegmentCollH;
  if (!isBckg) {  //do not recompute TOF on MC background
    iEvent.getByToken(muonCscSegmentToken_, CSCSegmentCollH);
    if (!CSCSegmentCollH.isValid()) {
      LogError(MOD) << "CSC Segment Collection not found!";
      return;
    }

    iEvent.getByToken(muonDtSegmentToken_, DTSegmentCollH);
    if (!DTSegmentCollH.isValid()) {
      LogError(MOD) << "DT Segment Collection not found!";
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
  
  // Collection for vertices
  vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken_);

  float CaloMET = -1, RecoPFMET_et = -1, RecoPFMHT = -1, HLTPFMET = -1, HLTPFMHT = -1;
  float RecoPFMET_eta = -1, RecoPFMET_phi = -1, RecoPFMET_significance = -1;

  //===================== Handle For PFMET ===================
  const edm::Handle<std::vector<reco::PFMET>> pfMETHandle = iEvent.getHandle(pfMETToken_);
  if (pfMETHandle.isValid() && !pfMETHandle->empty()) {
    for (unsigned int i = 0; i < pfMETHandle->size(); i++) {
      const reco::PFMET* pfMet = &(*pfMETHandle)[i];
      RecoPFMET_et = pfMet->et();
      RecoPFMET_eta = pfMet->eta();
      RecoPFMET_phi = pfMet->phi();
      RecoPFMET_significance = pfMet->significance();
    }
  }

  tuple->BS_RecoPFMET->Fill(RecoPFMET_et);

  //===================== Handle For CaloMET ===================
  const edm::Handle<std::vector<reco::CaloMET>> CaloMETHandle = iEvent.getHandle(CaloMETToken_);
  if (CaloMETHandle.isValid() && !CaloMETHandle->empty()) {
    for (unsigned int i = 0; i < CaloMETHandle->size(); i++) {
      const reco::CaloMET* calomet = &(*CaloMETHandle)[i];
      CaloMET = calomet->et();
    }
  }

  //===================== Handle For PFJet ===================
  unsigned int Jets_count = 0;
  const edm::Handle<reco::PFJetCollection> pfJetHandle = iEvent.getHandle(pfJetToken_);
  if (pfJetHandle.isValid() && !pfJetHandle->empty()) {
    const reco::PFJetCollection* pfJetColl = pfJetHandle.product();
    TLorentzVector pMHT;
    for (unsigned int i = 0; i < pfJetColl->size(); i++) {
      const reco::PFJet* jet = &(*pfJetColl)[i];
      Jets_count++;
      if (jet->pt() < 20 || abs(jet->eta()) > 5 ||
          jet->chargedEmEnergyFraction() + jet->neutralEmEnergyFraction() > 0.9)
        continue;
      TLorentzVector p4(jet->pt() * cos(jet->phi()), jet->pt() * sin(jet->phi()), 0, jet->et());
      pMHT += p4;
    }
    RecoPFMHT = pMHT.Pt();
  }
  //===================== Handle For PFCandidate ===================
  const edm::Handle<reco::PFCandidateCollection> pfCandHandle = iEvent.getHandle(pfCandToken_);

  //===================== Collection For Muons ===================
  // Get muon collections
  vector<reco::Muon> muonColl = iEvent.get(muonToken_);
  unsigned int Muons_count = 0; 
  float maxPtMuon1 = 0, maxPtMuon2 = 0;
  float etaMuon1 = 0, phiMuon1 = 0;
  float etaMuon2 = 0, phiMuon2 = 0;
  unsigned int muon1 = 0;
  for (unsigned int i = 0; i < muonColl.size(); i++) {
    const reco::Muon* mu = &(muonColl)[i];
    Muons_count++;
    if (mu->pt() > maxPtMuon1) {
      maxPtMuon1 = mu->pt();
      etaMuon1 = mu->eta();
      phiMuon1 = mu->phi();
      muon1 = i;
    }
  }
  for (unsigned int i = 0; i < muonColl.size(); i++) {
    if (i == muon1) continue;
    const reco::Muon* mu = &(muonColl)[i];
    if (mu->pt() > maxPtMuon2) {
      maxPtMuon2 = mu->pt();
      etaMuon2 = mu->eta();
      phiMuon2 = mu->phi();
    }
  }
  
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
  std::vector<bool> HSCP_passPreselection_noIsolation_noIh;
  std::vector<bool> HSCP_passPreselection;
  std::vector<bool> HSCP_passSelection;
  std::vector<float> HSCP_Charge;
  std::vector<float> HSCP_Pt;
  std::vector<float> HSCP_PtErr;
  std::vector<float> HSCP_Ias;
  std::vector<float> HSCP_Ias_noPix_noTIB_noTID_no3TEC;
  std::vector<float> HSCP_Ias_PixelOnly;
  std::vector<float> HSCP_Ias_StripOnly;
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
  std::vector<bool>  HSCP_isMuon;
  std::vector<int>   HSCP_MuonSelector;
  std::vector<bool>  HSCP_isElectron;
  std::vector<bool>  HSCP_isChHadron;
  std::vector<bool>  HSCP_isNeutHadron;
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
  std::vector<float> HSCP_iso_TK;
  std::vector<float> HSCP_iso_ECAL;
  std::vector<float> HSCP_iso_HCAL;
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

  //====================loop over HSCP candidates===================
  if (debug_ > 0 ) LogPrint(MOD) << "Loop over HSCP candidates:";
  unsigned int genTrack_count = 0;
  for (const auto& hscp : iEvent.get(hscpToken_)) {
    if (debug_> 0) LogPrint(MOD) << "  --------------------------------------------";
    genTrack_count++;
    if (debug_> 0) LogPrint(MOD) << "  >> This is general track " << genTrack_count;
    
    // Tracker only analysis must have either a tracker muon or a global muon
    if (typeMode_ == 1 &&
        !(hscp.type() == susybsm::HSCParticleType::trackerMuon || hscp.type() == susybsm::HSCParticleType::globalMuon)) {
      if (debug_ > 0 ) LogPrint(MOD) << "  >> Tracker only analysis  w/o a tracker muon or a global muon";
      continue;
    }
    
    // Tracker + Muon analysis  must have either a global muon
    if ((typeMode_ == 2 || typeMode_ == 4) && hscp.type() != susybsm::HSCParticleType::globalMuon) {
      if (debug_ > 0 ) LogPrint(MOD) << "  >> Tracker + Muon analysis w/o a global muon";
      continue;
    }
    
    // Define muon reference
    // For TOF only analysis we must have a muon connected to the HSCP candidate
    reco::MuonRef muon = hscp.muonRef();
    if (typeMode_ == 3 && muon.isNull()) {
      if (debug_> 0) LogPrint(MOD) << "  >> TOF only mode but no muon connected to the candidate -- skipping it";
      continue;
    }
    
    // Define track reference
    // For TOF only analysis use updated stand alone muon track, otherwise use inner tracker track
    reco::TrackRef track = (typeMode_ != 3) ? hscp.trackRef() : track = muon->standAloneMuon();
    
    // Skip events without track
    if (track.isNull()) {
      if (debug_> 0) LogPrint(MOD) << "  >> Event has no track associated to this HSCP, skipping it";
      continue;
    }

    // Require a track segment in the muon system
    if (typeMode_ > 1 && typeMode_ != 5 && (muon.isNull() || !muon->isStandAloneMuon())) {
      if (debug_> 0) LogPrint(MOD) << "  >> typeMode_ > 1 && typeMode_ != 5 && (muon.isNull() || !muon->isStandAloneMuon()), skipping it";
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
        continue;
    }
        
    // Reco - GEN track matching
    // For signal only, make sure that the candidate is associated to a true HSCP
    int closestGenIndex = -1;
    unsigned int closestHSCPsPDGsID = 0;
    unsigned int closestBackgroundPDGsIDs[2] = {0,0};
    if (isSignal && DistToHSCP(hscp, genColl, closestGenIndex, typeMode_) > 0.3) {
      if (debug_> 0) LogPrint(MOD) << "  >> Signal MC HSCP distance from gen to candidate is too big (" <<
        DistToHSCP(hscp, genColl, closestGenIndex, typeMode_) << "), skipping it";
      continue;
    } else if (isBckg) {
      float dRMinBckg = 9999.0;
      for (unsigned int g = 0; g < genColl.size(); g++) {
        float dr = deltaR(genColl[g].eta(),genColl[g].phi(),track->eta(),track->phi());
        if(dr < dRMinBckg){
          dRMinBckg = dr;
          closestBackgroundPDGsIDs[0] = abs(genColl[g].pdgId());
          for (unsigned int numMom = 0; numMom < genColl[g].numberOfMothers(); numMom++) {
            if (abs(genColl[g].mother(numMom)->pdgId())  != abs(genColl[g].pdgId()))
            {
              closestBackgroundPDGsIDs[1] = abs(genColl[g].mother(numMom)->pdgId());
              break;
            }
          }
          closestGenIndex = g;
        }
      }
    }
      
    if (isSignal) {
      closestHSCPsPDGsID = abs(genColl[closestGenIndex].pdgId());
      // All HSCP candidates
      tuple->HSCPCandidateType->Fill(0.5, EventWeight_);
      // Neutral HSCP candidates
      if (   closestHSCPsPDGsID == 1000993 || closestHSCPsPDGsID == 1009113
          || closestHSCPsPDGsID == 1009223 || closestHSCPsPDGsID == 1009313
          || closestHSCPsPDGsID == 1009333 || closestHSCPsPDGsID == 1092114
          || closestHSCPsPDGsID == 1093214 || closestHSCPsPDGsID == 1093324
          || closestHSCPsPDGsID == 1000622 || closestHSCPsPDGsID == 1000642
          || closestHSCPsPDGsID == 1006113 || closestHSCPsPDGsID == 1006311
          || closestHSCPsPDGsID == 1006313 || closestHSCPsPDGsID == 1006333) {
        tuple->HSCPCandidateType->Fill(1.5, EventWeight_);
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
        tuple->HSCPCandidateType->Fill(2.5, EventWeight_);
      }
      // Double-charged R-hadrons
      else if (closestHSCPsPDGsID == 1092224 || closestHSCPsPDGsID == 1006223) {
        tuple->HSCPCandidateType->Fill(3.5, EventWeight_);
        // Dont mix double charged R-hadrons with the rest
        // The reco pt of them is 1/2 the pt of the gen track
        continue;
      }
      // tau prime, could be single or multiple charged
      else if (closestHSCPsPDGsID == 17) {
        tuple->HSCPCandidateType->Fill(4.5, EventWeight_);
      }
      else {
        tuple->HSCPCandidateType->Fill(5.5, EventWeight_);
      }
    }
    
    // TODO this is repeated in the pre-selection
    int highestPtGoodVertex = -1;
    int goodVerts = 0;
    float dzMin = 10000;
    // Loop on the vertices in the event
    for (unsigned int i = 0; i < vertexColl.size(); i++) {
      if (vertexColl[i].isFake() || fabs(vertexColl[i].z()) > 24 || vertexColl[i].position().rho() > 2 ||
          vertexColl[i].ndof() <= 4)
        continue;  //only consider good vertex
      goodVerts++;
      if (fabs(track->dz(vertexColl[i].position())) < fabs(dzMin)) {
        dzMin = fabs(track->dz(vertexColl[i].position()));
        highestPtGoodVertex = i;
      }
    } // End loop on the vertices in the event
    if (highestPtGoodVertex < 0)
      highestPtGoodVertex = 0;

    // Impact paramters dz and dxy
    float dz = track->dz(vertexColl[highestPtGoodVertex].position());
    float dxy = track->dxy(vertexColl[highestPtGoodVertex].position());

    // Compute transverse mass mT between HSCP with and MET
    float massT = sqrt(2*track->pt()*RecoPFMET_et*(1-cos(track->phi()-RecoPFMET_phi))); 
    HSCP_mT.push_back(massT);
  
    // Save PF informations and isolation
    float pfIsolation_DZ_ = 0.1;
   
    float track_PFIso005_sumCharHadPt = 0, track_PFIso005_sumNeutHadPt = 0, track_PFIso005_sumPhotonPt = 0, track_PFIso005_sumPUPt = 0;
    float track_PFIso01_sumCharHadPt = 0, track_PFIso01_sumNeutHadPt = 0, track_PFIso01_sumPhotonPt = 0, track_PFIso01_sumPUPt = 0;
    float track_PFIso03_sumCharHadPt = 0, track_PFIso03_sumNeutHadPt = 0, track_PFIso03_sumPhotonPt = 0, track_PFIso03_sumPUPt = 0;
    float track_PFIso05_sumCharHadPt = 0, track_PFIso05_sumNeutHadPt = 0, track_PFIso05_sumPhotonPt = 0, track_PFIso05_sumPUPt = 0;

    float RMin = 9999.;
    unsigned int idx_pf_RMin = 9999;

    bool pf_isMuon = false, pf_isElectron = false, pf_isChHadron = false, pf_isNeutHadron = false;
    int pf_muon_selector = -1;
    float pf_ecal_energy = 0, pf_hcal_energy = 0;

    if(pfCandHandle.isValid() && !pfCandHandle->empty()) {
      const reco::PFCandidateCollection* pf = pfCandHandle.product();
      for (unsigned int i = 0; i < pf->size(); i++){
          const reco::PFCandidate* pfCand = &(*pf)[i];
          float dr = deltaR(pfCand->eta(),pfCand->phi(),track->eta(),track->phi());
          if(dr < RMin){
              RMin = dr;
              idx_pf_RMin = i;
           }
       }//end loop PFCandidates

      // https://github.com/cms-sw/cmssw/blob/72d0fc00976da53d1fb745eb7f37b2a4ad965d7e/
      // PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc#L555
      for(unsigned int i=0;i<pf->size();i++){
        const reco::PFCandidate* pfCand = &(*pf)[i];
        if(i == idx_pf_RMin) {
            pf_isMuon = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;
            pf_isElectron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::e;
            pf_isChHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h;
            pf_isNeutHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h0;
            pf_ecal_energy = pfCand->ecalEnergy();
            pf_hcal_energy = pfCand->hcalEnergy();
        }
        if(i == idx_pf_RMin) continue; //don't count itself
        float dr = deltaR(pfCand->eta(),pfCand->phi(),track->eta(),track->phi());
        bool fromPV = (fabs(dz) < pfIsolation_DZ_);
        int id = std::abs(pfCand->pdgId());
        float pt = pfCand->p4().pt();
        if(dr<0.05){
            // charged cands from PV get added to trackIso
            if(id == 211 && fromPV) track_PFIso005_sumCharHadPt+=pt;
            // charged cands not from PV get added to pileup iso
            else if(id == 211) track_PFIso005_sumPUPt+=pt;
            // neutral hadron iso
            if(id == 130) track_PFIso005_sumNeutHadPt+=pt;
            // photon iso
            if(id == 22) track_PFIso005_sumPhotonPt+=pt;
        }if(dr<0.1){
            if(id == 211 && fromPV) track_PFIso01_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso01_sumPUPt+=pt;
            if(id == 130) track_PFIso01_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso01_sumPhotonPt+=pt;
        }if(dr<0.3){
            if(id == 211 && fromPV) track_PFIso03_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso03_sumPUPt+=pt;
            if(id == 130) track_PFIso03_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso03_sumPhotonPt+=pt;
        }if(dr<0.5){
            if(id == 211 && fromPV) track_PFIso05_sumCharHadPt+=pt;
            else if(id == 211) track_PFIso05_sumPUPt+=pt;
            if(id == 130) track_PFIso05_sumNeutHadPt+=pt;
            if(id == 22) track_PFIso05_sumPhotonPt+=pt;
        }
      }//end loop PFCandidates
    }

    HSCP_count++;
    if (debug_> 0) LogPrint(MOD) << "  >> This is HSCP candidate track " << HSCP_count ;
    std::vector<float> clust_charge;
    std::vector<float> clust_pathlength;
    std::vector<unsigned int> clust_nstrip;
    std::vector<bool> clust_sat254;
    std::vector<bool> clust_sat255;
    std::vector<uint32_t> clust_detid;
    std::vector<bool> clust_isStrip;
    std::vector<bool> clust_isPixel;
    
    // Include probQonTrack, probXYonTrack, probQonTrackNoLayer1, probXYonTrackNoLayer1 into one array
    float pixelProbs[4] = {0.0,0.0,0.0,0.0};
    
    int numRecHits = 0;
    int numRecHitsNoLayer1 = 0;
    float probQonTrackWMulti = 1;
    float probXYonTrackWMulti = 1;
    float probQonTrackWMultiNoLayer1 = 1;
    float probXYonTrackWMultiNoLayer1 = 1;

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
        LogPrint(MOD) << "  >> No associated gen track to this candidatei -- this is a fake track, skipping it"; 
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

    if (typeMode_ > 1 && typeMode_ != 5 && !hscp.muonRef().isNull()) {
      if (isBckg) {
        tof = &(*tofMap)[hscp.muonRef()];
        dttof = &(*tofDtMap)[hscp.muonRef()];
        csctof = &(*tofCscMap)[hscp.muonRef()];
      } else {
        const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollH;
        const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollH;
        //tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, isData?1:0 ); //apply T0 correction on data but not on signal MC
        tofCalculator.computeTOF(
            muon, CSCSegmentColl, DTSegmentColl, 1);  //apply T0 correction on data but not on signal MC
        tof = &tofCalculator.combinedTOF;
        dttof = &tofCalculator.dtTOF;
        csctof = &tofCalculator.cscTOF;
      }
    } // end conditions for TOF including analysis variables

    // skip tracks without hits otherwise there will be a crash
    if (!dedxHits) {
      if (debug_> 3) LogPrint(MOD) << "No dedxHits associated to this track, skipping it";
      continue;
    }

    int nofClust_dEdxLowerThan = 0;
    float factorChargeToE = 3.61 * pow(10, -6) * 247;

    // Loop through the rechits on the given track
    for (unsigned int i = 0; i < dedxHits->size(); i++) {
      clust_charge.push_back(dedxHits->charge(i));
      clust_pathlength.push_back(dedxHits->pathlength(i));
      clust_isStrip.push_back(dedxHits->detId(i) >= 3 ? true : false);
      clust_isPixel.push_back(dedxHits->detId(i) >= 3 ? false : true);
      clust_detid.push_back(dedxHits->detId(i));
      DetId detid(dedxHits->detId(i));
  
      if (detid.subdetId() < 3) {
        // Calculate probQ and probXY for this pixel rechit
        // Taking the pixel cluster
        auto const* pixelCluster =  dedxHits->pixelCluster(i);
        if (pixelCluster == nullptr) {
          if (debug_> 0) LogPrint(MOD) << "    >> No dedxHits associated to this pixel cluster, skipping it";
          if (debug_> 0) LogPrint(MOD) << "    >> At this point this should never happen";
          continue;
        }
        // Check on which geometry unit the hit is
        const GeomDetUnit& geomDet = *tkGeometry->idToDetUnit(detid);
        // Get the local vector for the track direction
        LocalVector lv = geomDet.toLocal(GlobalVector(track->px(), track->py(), track->pz()));
        // Re-run the CPE on this cluster with the lv above
        auto reCPE = std::get<2>(pixelCPE->getParameters(
                                                         *pixelCluster, geomDet, LocalTrajectoryParameters(dedxHits->pos(i), lv, track->charge())));
        // extract probQ and probXY from this
        float probQ = SiPixelRecHitQuality::thePacking.probabilityQ(reCPE);
        float probXY = SiPixelRecHitQuality::thePacking.probabilityXY(reCPE);
        if (probQ > 0.f) {
          numRecHits++;
          // Calculate alpha term needed for the combination
          probQonTrackWMulti *= probQ;
          probXYonTrackWMulti *= probXY;
        }
        
        // Have a separate variable that excludes Layer 1
        // Layer 1 was very noisy in 2017/2018
        if (( detid.subdetId() == PixelSubdetector::PixelEndcap) || (detid.subdetId() == PixelSubdetector::PixelBarrel &&
                                                                     tTopo->pxbLayer(detid) != 1)) {
          float probQNoLayer1 = SiPixelRecHitQuality::thePacking.probabilityQ(reCPE);
          float probXYNoLayer1 = SiPixelRecHitQuality::thePacking.probabilityXY(reCPE);
          if (probQNoLayer1 > 0.f) {  // only save the non-zero rechits
            numRecHitsNoLayer1++;
            // Calculate alpha term needed for the combination
            probQonTrackWMultiNoLayer1 *= probQNoLayer1;
            probXYonTrackWMultiNoLayer1 *= probXYNoLayer1;
          }
        }
      } else if (detid.subdetId() >= 3) {
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
        if (dedxHits->charge(i) * factorChargeToE / dedxHits->pathlength(i) < theFMIPX_)
          nofClust_dEdxLowerThan++;
      }
    } // end loop on rechits on the given track

    // Combine probQ-s into HSCP candidate (track) level quantity
    pixelProbs[0] = combineProbs(probQonTrackWMulti, numRecHits);
    pixelProbs[1] = combineProbs(probXYonTrackWMulti, numRecHits);
    pixelProbs[2] = combineProbs(probQonTrackWMultiNoLayer1, numRecHitsNoLayer1);
    pixelProbs[3] = combineProbs(probXYonTrackWMultiNoLayer1, numRecHitsNoLayer1);
      
    // Cleaning of tracks that had failed the template CPE (prob <= 0.0 and prob >= 1.0 cases)
    if (pixelProbs[0] <= 0.0 || pixelProbs[1] <= 0.0 || pixelProbs[0] >= 1.00000001 || pixelProbs[1] >= 1.000000001) {
      if (debug_> 2) LogPrint(MOD) << "    >> Probs out of bound: " <<
        " ProbQ = " << pixelProbs[0] << " ProbXY = " << pixelProbs[1] <<  " ProbQNoL1 = "<< pixelProbs[2] << " ProbXYNoL1 = " << pixelProbs[3];
   //   continue;
    }
    
    TreeprobQonTrack = pixelProbs[0];
    TreeprobQonTracknoL1 = pixelProbs[1];
    TreeprobXYonTrack = pixelProbs[2];
    TreeprobXYonTracknoL1 = pixelProbs[3];

    float Fmip = (float)nofClust_dEdxLowerThan / (float)dedxHits->size();


    //computedEdx: hits, SF, templates, usePixel, useClusterCleaning, reverseProb, uneTrunc, TrackerGains,
    //             useStrips, mustBeInside, MaxStripNOM, correctFEDSat, XtalkInv, lowDeDxDrop, hipEmul, dedxErr, closestHSCPsPDGsID, skipPix, useTemplateLayer_, skipPixel_L1, DeDxprobQ, skip_templ_Ias
    //
    //correction inverseXtalk = 0 --> take the raw amplitudes of the cluster
    //correction inverseXtalk = 1 --> modify the amplitudes based on xtalk for non-saturated cluster + correct for saturation
    //
    //skip templates_Ias = 0 --> no skip
    //skip templates_Ias = 1 --> no Pix, no TIB, no TID, no 3 first layers TEC
    //skip templates_Ias = 2 --> Pixel Only

    string year = period_;
    if(!isData) year="";
    int run_number=iEvent.id().run();
    //Ias
    float dEdxErr = 0;
    auto dedxSObjTmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, typeMode_ == 5, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.00, nullptr, 0, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_);
    
    reco::DeDxData* dedxSObj = dedxSObjTmp.numberOfMeasurements() > 0 ? &dedxSObjTmp : nullptr;

    //Ih
    auto dedxMObjTmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, true, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.0, nullptr, &dEdxErr, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_);
    
    reco::DeDxData* dedxMObj = dedxMObjTmp.numberOfMeasurements() > 0 ? &dedxMObjTmp : nullptr;

    // Ih Up
    auto dedxMUpObjTmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, true, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.15, nullptr, 0, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_);
    
    reco::DeDxData* dedxMUpObj = dedxMUpObjTmp.numberOfMeasurements() > 0 ? &dedxMUpObjTmp : nullptr;

    // Ih Down
    auto dedxMDownObjTmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, true, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.15, nullptr, 0, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_);
    
    reco::DeDxData* dedxMDownObj = dedxMDownObjTmp.numberOfMeasurements() > 0 ? &dedxMDownObjTmp : nullptr;

    // Ih no pixel L1
    auto dedxIh_noL1_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, true, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                      true, true, 99, false, 1, 0.0, nullptr, &dEdxErr, closestHSCPsPDGsID, false, useTemplateLayer_, true);
    
    reco::DeDxData* dedxIh_noL1 = dedxIh_noL1_Tmp.numberOfMeasurements() > 0 ? &dedxIh_noL1_Tmp : nullptr;

    // Ih 0.15 low values drop
    auto dedxIh_15drop_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, true, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                      true, true, 99, false, 1, 0.15, nullptr, &dEdxErr, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_);
    
    reco::DeDxData* dedxIh_15drop = dedxIh_15drop_Tmp.numberOfMeasurements() > 0 ? &dedxIh_15drop_Tmp : nullptr;

    // Ih Strip only
    auto dedxIh_StripOnly_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, false, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.0, nullptr, &dEdxErr, closestHSCPsPDGsID, true, useTemplateLayer_);

    reco::DeDxData* dedxIh_StripOnly = dedxIh_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIh_StripOnly_Tmp : nullptr;

    // Ih Strip only and 0.15 low values drop
    auto dedxIh_StripOnly_15drop_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, false, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.15, nullptr, &dEdxErr, closestHSCPsPDGsID, true, useTemplateLayer_, true);
    
    reco::DeDxData* dedxIh_StripOnly_15drop = dedxIh_StripOnly_15drop_Tmp.numberOfMeasurements() > 0 ? &dedxIh_StripOnly_15drop_Tmp : nullptr;

    // Ih correct saturation from fits
    auto dedxIh_SaturationCorrectionFromFits_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, nullptr, false, useClusterCleaning, false, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 2, 0.0, nullptr, &dEdxErr, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_, true);

    reco::DeDxData* dedxIh_SaturationCorrectionFromFits = dedxIh_SaturationCorrectionFromFits_Tmp.numberOfMeasurements() > 0 ? &dedxIh_SaturationCorrectionFromFits_Tmp : nullptr;

    //dEdx probQ discriminator based on templates (same than Ias)
    auto dedx_probQ_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, true, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.00, nullptr, 0, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_, true, true);

    reco::DeDxData* dedx_probQ = dedx_probQ_Tmp.numberOfMeasurements() > 0 ? &dedx_probQ_Tmp : nullptr;


    //Ias without TIB, TID, and 3 first TEC layers
    auto dedxIas_noTIBnoTIDno3TEC_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, true, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.00, nullptr, 0, closestHSCPsPDGsID, skipPixel_, useTemplateLayer_, true, false, 1);

    reco::DeDxData* dedxIas_noTIBnoTIDno3TEC = dedxIas_noTIBnoTIDno3TEC_Tmp.numberOfMeasurements() > 0 ? &dedxIas_noTIBnoTIDno3TEC_Tmp : nullptr;

    //Ias Pixel only

    auto dedxIas_PixelOnly_Tmp =
        computedEdx(run_number, year, dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, true, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.00, nullptr, 0, closestHSCPsPDGsID, false, useTemplateLayer_, false, false, 2);

    reco::DeDxData* dedxIas_PixelOnly = dedxIas_PixelOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIas_PixelOnly_Tmp : nullptr;

    //Ias Strip only
    auto dedxIas_StripOnly_Tmp = 

        computedEdx(run_number, year, dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, true, false, trackerCorrector.TrackerGains,
                    true, true, 99, false, 1, 0.00, nullptr, 0, closestHSCPsPDGsID, true, useTemplateLayer_);

    reco::DeDxData* dedxIas_StripOnly = dedxIas_StripOnly_Tmp.numberOfMeasurements() > 0 ? &dedxIas_StripOnly_Tmp : nullptr;

    //Choose of Ih definition - Ih_nodrop_noPixL1
    dedxMObj = dedxIh_noL1;

    OpenAngle = deltaROpositeTrack(iEvent.get(hscpToken_), hscp);
    
    float MassErr = GetMassErr(track->p(),
                                track->ptError(),
                                dedxMObj ? dedxMObj->dEdx() : -1,
                                dEdxErr,
                                GetMass(track->p(), dedxMObj ? dedxMObj->dEdx() : -1, dEdxK_, dEdxC_),
                                dEdxK_,
                                dEdxC_);
    
    // ------------------------------------------------------------------------------------
    //compute systematic uncertainties on signal
    if (isSignal) {
      float genpT = -1.0;
        // Loop through the gen collection
      for (auto const gen : genColl) {
        if (!isGoodGenHSCP(gen,false)) continue;
        float separation = deltaR(track->eta(), track->phi(), gen.eta(), gen.phi());
        if (separation > 0.3) {
          continue;
        } else {
          genpT = gen.pt();
          break;
        }
      }
      if (genpT > 0) {
        // 2D plot to compare gen pt vs reco pt
        tuple->genrecopT->Fill(genpT, track->pt());
      }
      if (calcSyst_) {
        calculateSyst(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, pixelProbs, EventWeight_, tuple, -1, MassErr, true, closestBackgroundPDGsIDs);
      }
    }  //End of systematic computation for signal
    // ------------------------------------------------------------------------------------
    
    if (debug_ > 5 ) LogPrint(MOD)  << "        >> dEdxK_: " << dEdxK_ << " dEdxC_: " << dEdxC_;
    // Check if we pass the preselection
    if (debug_ > 2) LogPrint(MOD)  << "      >> Check if we pass Preselection";
    bool Ih_Iso_cut = true;

    bool passPre = passPreselection(
                          track,
                          dedxHits,
                          dedxSObj,
                          dedxMObj,
                          tof,
                          iEvent,
                          pixelProbs,
                          EventWeight_,
                          tuple,
                          isSignal ? genColl[closestGenIndex].p() / genColl[closestGenIndex].energy() : -1,
                          false,
                          0,
                          0,
                          MassErr,
                          Ih_Iso_cut,
                          closestBackgroundPDGsIDs);
    
    Ih_Iso_cut = false;
    bool passPre_noIh_noIso = passPreselection(
                          track,
                          dedxHits,
                          NULL,
                          NULL,
                          tof,
                          iEvent,
                          pixelProbs,
                          EventWeight_,
                          NULL,
                          isSignal ? genColl[closestGenIndex].p() / genColl[closestGenIndex].energy() : -1,
                          false,
                          0,
                          0,
                          MassErr,
                          Ih_Iso_cut,
                          closestBackgroundPDGsIDs);

    
      // Dont do TOF only is isCosmicSB is true
    if (typeMode_ == 5 && isCosmicSB) {
      if (debug_ > 2) LogPrint(MOD) << "      >> This is a cosmic track, skipping it";
      continue;
    } else if (isCosmicSB) {
      if (debug_ > 2) LogPrint(MOD) << "      >> This is a cosmic track, please check what's up";
    }
    
    // Dont do TOF only is isSemiCosmicSB is true
    if (typeMode_ == 5 && isSemiCosmicSB) {
      if (debug_ > 2) LogPrint(MOD) << "      >> This is a semi-cosmic track, skipping it";
      continue;
    } else if (isSemiCosmicSB) {
      if (debug_ > 2) LogPrint(MOD) << "      >> This is a semi-cosmic track, please check what's up";
    }
    
    //fill the ABCD histograms and a few other control plots
    //WAIT//else if(isBckg) Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, MCTrPlots);

    
    if (passPre) {
      if (debug_ > 2) LogPrint(MOD) << "      >> Passed pre-selection";
      if (debug_ > 2) LogPrint(MOD) << "      >> Fill control and prediction histos";
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
    } else {
      if (debug_ > 2) LogPrint(MOD) << "      >> Preselection not passed";
      continue;
    }
    
    //Find the number of tracks passing selection for TOF<1 that will be used to check the background prediction
    //float Mass = -1;
    if (isBckg || isData) {
      //compute the mass of the candidate, for TOF mass flip the TOF over 1 to get the mass, so 0.8->1.2
      float Mass = -1;
      if (dedxMObj)
        Mass = GetMass(track->p(), dedxMObj->dEdx(), dEdxK_, dEdxC_);
      float MassTOF = -1;
      if (tof)
        MassTOF = GetTOFMass(track->p(), (2 - tof->inverseBeta()));
      float MassComb = -1;
      if (tof && dedxMObj)
        MassComb = GetMassFromBeta(track->p(),
                                   (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / (2 - tof->inverseBeta()))) * 0.5);
      if (dedxMObj)
        MassComb = Mass;
      if (tof)
        MassComb = GetMassFromBeta(track->p(), (1 / (2 - tof->inverseBeta())));
      //Background check looking at region with TOF<1
      for (unsigned int CutIndex = 0; CutIndex < CutPt_Flip_.size(); CutIndex++) {
        //Fill Mass Histograms
        tuple->Mass_Flip->Fill(CutIndex, Mass, EventWeight_);
        if (tof)
          tuple->MassTOF_Flip->Fill(CutIndex, MassTOF, EventWeight_);
        tuple->MassComb_Flip->Fill(CutIndex, MassComb, EventWeight_);
      }
    }

    //compute the mass of the candidate
    float Mass =  dedxMObj ?  GetMass(track->p(), dedxMObj->dEdx(), dEdxK_, dEdxC_) : -1;
    float MassTOF = tof    ?  GetTOFMass(track->p(), tof->inverseBeta()) : -1;
    float MassComb  = (tof && dedxMObj) ? GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5) : -1;
    MassComb = tof   ?  GetMassFromBeta(track->p(), (1 / tof->inverseBeta())) : Mass;

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
      if (debug_ > 3 ) LogPrint(MOD) << "      >> We enter the selection cut loop now";
      //==========================================================
      // Cut loop: over all possible selection (one of them, the optimal one, will be used later)
      for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
        //Full Selection
        if (!passSelection(track,
                           dedxSObj,
                           dedxMObj,
                           tof,
                           iEvent,
                           EventWeight_,
                           CutIndex,
                           tuple,
                           false,
                           isSignal ? genColl[closestGenIndex].p() / genColl[closestGenIndex].energy() : -1,
                           false,
                           0,
                           0)) {
          if (debug_ > 6 ) LogPrint(MOD) << "        >> Selection failed, skipping this CutIndex = " << CutIndex;
          continue;
        } else {
          if (debug_ > 6 ) LogPrint(MOD) << "        >> Selection passed with CutIndex = " << CutIndex;
        }
        if (CutIndex != 0)
          PassNonTrivialSelection = true;
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
        if (tof)
          tuple->MassTOF->Fill(CutIndex, MassTOF, EventWeight_);
        if (isBckg)
          tuple->MassComb->Fill(CutIndex, MassComb, EventWeight_);

        //Fill Mass Histograms for different Ih syst
        tuple->Mass_SystHUp->Fill(CutIndex, MassUp, EventWeight_);
        tuple->Mass_SystHDown->Fill(CutIndex, MassDown, EventWeight_);
        if (tof)
          tuple->MassTOF_SystH->Fill(CutIndex, MassTOF, EventWeight_);
        tuple->MassComb_SystHUp->Fill(CutIndex, MassUpComb, EventWeight_);
        tuple->MassComb_SystHDown->Fill(CutIndex, MassDownComb, EventWeight_);

      }  //end of Cut loop
    } // end of condition for passPre

    float Ick2 = 0;
    if (dedxMObj)
      Ick2 = GetIck(dedxMObj->dEdx(), dEdxK_, dEdxC_);
    int nomh = 0;
    nomh = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
           track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
    float fovhd = track->found() <= 0 ? -1 : track->found() / float(track->found() + nomh);
    unsigned int nom = 0;
    if (dedxSObj)
      nom = dedxSObj->numberOfMeasurements();

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

    HSCP_passCutPt55.push_back(track->pt() > 55 ? true : false);
    HSCP_passPreselection_noIsolation_noIh.push_back(passPre_noIh_noIso);
    HSCP_passPreselection.push_back(passPre);
    HSCP_passSelection.push_back(PassNonTrivialSelection);
    HSCP_Charge.push_back(track->charge());
    HSCP_Pt.push_back(track->pt());
    HSCP_PtErr.push_back(track->ptError());
    HSCP_Ias.push_back(dedxSObj ? dedxSObj->dEdx() : -1);
    HSCP_Ias_noPix_noTIB_noTID_no3TEC.push_back(dedxIas_noTIBnoTIDno3TEC ? dedxIas_noTIBnoTIDno3TEC->dEdx() : -1);
    HSCP_Ias_PixelOnly.push_back(dedxIas_PixelOnly ? dedxIas_PixelOnly->dEdx() : -1);
    HSCP_Ias_StripOnly.push_back(dedxIas_StripOnly ? dedxIas_StripOnly->dEdx() : -1);
    HSCP_Ih.push_back(dedxMObj ? dedxMObj->dEdx() : -1);
    HSCP_Ick.push_back(dedxMObj ? Ick2 : -99);
    HSCP_Fmip.push_back(Fmip);
    HSCP_ProbXY.push_back(TreeprobXYonTrack);
    HSCP_ProbXY_noL1.push_back(TreeprobXYonTracknoL1);
    HSCP_ProbQ.push_back(TreeprobQonTrack);
    HSCP_ProbQ_noL1.push_back(TreeprobQonTracknoL1);
    HSCP_ProbQ_dEdx.push_back(dedx_probQ ? dedx_probQ->dEdx() : -1);
    HSCP_Ndof.push_back(track->ndof());
    HSCP_Chi2.push_back(track->chi2());
    HSCP_QualityMask.push_back(track->qualityMask());
    HSCP_isHighPurity.push_back(track->quality(reco::TrackBase::highPurity));
    HSCP_isMuon.push_back(pf_isMuon);
    HSCP_MuonSelector.push_back(pf_muon_selector);
    HSCP_isElectron.push_back(pf_isElectron);
    HSCP_isChHadron.push_back(pf_isChHadron);
    HSCP_isNeutHadron.push_back(pf_isNeutHadron);
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
    HSCP_dR.push_back(OpenAngle);
    HSCP_eta.push_back(track->eta());
    HSCP_phi.push_back(track->phi());
    HSCP_NOH.push_back(track->found());
    HSCP_NOPH.push_back(track->hitPattern().numberOfValidPixelHits());
    HSCP_FOVH.push_back(track->validFraction());
    HSCP_NOMH.push_back(nomh);
    HSCP_FOVHD.push_back(fovhd);
    HSCP_NOM.push_back(nom);
    HSCP_iso_TK.push_back(iso_TK);
    HSCP_iso_ECAL.push_back(iso_ECAL);
    HSCP_iso_HCAL.push_back(iso_HCAL);
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

  }  //END loop over HSCP candidates

  tuple_maker->fillTreeBranches(tuple,
                                TrigInfo_,
                                iEvent.id().run(),
                                iEvent.id().event(),
                                iEvent.id().luminosityBlock(),
                                pileup_fromLumi,
                                vertexColl.size(),
                                HSCP_count,
                                Muons_count,
                                Jets_count,
                                EventWeight_,
                                GeneratorWeight_,
                                GeneratorBinningValues_,
                                HLT_Mu50,
                                HLT_PFMET120_PFMHT120_IDTight,
                                HLT_PFHT500_PFMET100_PFMHT100_IDTight,
                                HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
                                HLT_MET105_IsoTrk50,
                                CaloMET,
                                RecoPFMET_et,
                                RecoPFMHT,
                                HLTPFMET,
                                HLTPFMHT,
                                RecoPFMET_eta,
                                RecoPFMET_phi,
                                RecoPFMET_significance,
                                maxPtMuon1,
                                etaMuon1,
                                phiMuon1,
                                maxPtMuon2,
                                etaMuon2,
                                phiMuon2,
                                HSCP_mT,
                                HSCP_passCutPt55,
                                HSCP_passPreselection_noIsolation_noIh,
                                HSCP_passPreselection,
                                HSCP_passSelection,
                                HSCP_Charge,
                                HSCP_Pt,
                                HSCP_PtErr,
                                HSCP_Ias,
                                HSCP_Ias_PixelOnly,
                                HSCP_Ias_StripOnly,
                                HSCP_Ias_noPix_noTIB_noTID_no3TEC,
                                HSCP_Ih,
                                HSCP_Ick,
                                HSCP_Fmip,
                                HSCP_ProbXY,
                                HSCP_ProbXY_noL1,
                                HSCP_ProbQ,
                                HSCP_ProbQ_noL1,
                                HSCP_ProbQ_dEdx,
                                HSCP_Ndof,
                                HSCP_Chi2,
                                HSCP_QualityMask,
                                HSCP_isHighPurity,
                                HSCP_isMuon,
                                HSCP_MuonSelector,
                                HSCP_isElectron,
                                HSCP_isChHadron,
                                HSCP_isNeutHadron,
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
                                HSCP_eta,
                                HSCP_phi,
                                HSCP_NOH,
                                HSCP_NOPH,
                                HSCP_FOVH,
                                HSCP_NOMH,
                                HSCP_FOVHD,
                                HSCP_NOM,
                                HSCP_iso_TK,
                                HSCP_iso_ECAL,
                                HSCP_iso_HCAL,
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
                                HSCP_GenPhi);

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
    if (HSCPTk_SystT[CutIndex]) {
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
  }
}

// ------------ method called once each job just after ending the event loop  ------------
void Analyzer::endJob() {
  delete RNG;
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
  desc.add("TriggerResults", edm::InputTag("TriggerResults","","HLT"))
    ->setComment("A");
  desc.add("PfMET", edm::InputTag("pfMet"))
    ->setComment("A");
  desc.add("PfJet", edm::InputTag("ak4PFJetsCHS"))
    ->setComment("A");
  desc.add("CaloMET", edm::InputTag("caloMet"))
    ->setComment("A");
  desc.add("PileupInfo", edm::InputTag("addPileupInfo"))
    ->setComment("A");
  desc.add("GenParticleCollection", edm::InputTag("genParticlesSkimmed"))
    ->setComment("A");
  desc.add("TrackToGenAssoc", edm::InputTag("allTrackMCMatch"))
    ->setComment("Collection used to match to gen thruth");
  desc.add("PfCand", edm::InputTag("particleFlow"))
    ->setComment("Input collection for particleFlow algorithm");
  desc.add("GenCollection", edm::InputTag("GenParticleCollection"))
    ->setComment("A");
  desc.addUntracked("Trigger_Mu", std::vector<std::string>{"HLT_Mu50_v"})
    ->setComment("Add the list of muon triggers");
  desc.addUntracked("Trigger_MET",  std::vector<std::string>{"HLT_PFMET120_PFMHT120_IDTight_v","HLT_PFHT500_PFMET100_PFMHT100_IDTight_v","HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60_v","HLT_MET105_IsoTrk50_v"})
    ->setComment("Add the list of MET triggers");
  desc.addUntracked("TypeMode", 0)
    ->setComment("0:Tk only, 1:Tk+Muon, 2:Tk+TOF, 3:TOF onlypwd, 4:Q>1, 5:Q<1");
  desc.addUntracked("SampleType", 0)
    ->setComment("0:Data, 1:Background, 2:Signal, 3:Signal Systematics");
  desc.addUntracked<std::string>("SampleName","BaseName")->setComment("A");
  desc.addUntracked<std::string>("Period","2017")->setComment("A");
  desc.addUntracked("SkipSelectionPlot",false)->setComment("A");
  desc.addUntracked("PtHistoUpperBound",4000.0)->setComment("A");
  desc.addUntracked("MassHistoUpperBound",4000.0)->setComment("A");
  desc.addUntracked("MassNBins",40)->setComment("Number of bins in the mass plot");
  desc.addUntracked("IPbound",1.0)
    ->setComment("Number of different Dz side regions used to make cosmic background prediction");
  desc.addUntracked("PredBins",0)
  ->setComment("How many different bins the prediction is split in for analysis being run, sets how many histograms are actually initialized.");
  desc.addUntracked("EtaBins",60)
  ->setComment("How many bins we use for the background prediction method in Eta -- impacts background prediction method -- histograms with the name of the form Pred_Eta in Analysis_PlotStructure.h");
  desc.addUntracked("DeDxS_UpLim",1.0)->setComment("A");
  desc.addUntracked("DeDxM_UpLim",30.0)->setComment("A");
  desc.addUntracked("DzRegions",6)->setComment("A");
  desc.addUntracked("GlobalMinTOF",1.0)->setComment("A");
  desc.addUntracked("SkipPixel",true)->setComment("A");
  desc.addUntracked("UseTemplateLayer",false)->setComment("A");
  desc.addUntracked("DeDxSF_0",1.0)->setComment("A");
  desc.addUntracked("DeDxSF_1",1.0325)->setComment("A");
  desc.addUntracked("DeDxK",2.3)->setComment("A");
  desc.addUntracked("DeDxC",3.17)->setComment("A");
  desc.addUntracked("FMIPX",4.0)->setComment("A");
  desc.addUntracked("SaveTree",0)->setComment("0: do not save tree, 6: everything is saved");
  desc.addUntracked("SaveGenTree",0)->setComment("A");
  desc.addUntracked("EnableDeDxCalibration",false)->setComment("A");
  desc.addUntracked<std::string>("DeDxCalibration","SUSYBSMAnalysis/HSCP/data/Data13TeVGains_v2.root")
    ->setComment("Second gain calibration for strips");
  desc.addUntracked<std::string>("DeDxTemplate","SUSYBSMAnalysis/HSCP/data/template_2017B.root")
    ->setComment("Ias vs Pt templates in eta binning");
  desc.addUntracked<std::string>("Geometry","SUSYBSMAnalysis/HSCP/data/CMS_GeomTree.txt")
  ->setComment("MuonTimeOffset info"); // I'm not sure we need this
  desc.addUntracked<std::string>("TimeOffset","SUSYBSMAnalysis/HSCP/data/MuonTimeOffset.txt")
    ->setComment("MuonTimeOffset info"); // I'm not sure we need this
  desc.add<std::string>("PixelCPE","PixelCPEClusterRepair")
    ->setComment("CPE used in the pixel reco, cluster repair is the best available so far ");
  desc.addUntracked("DebugLevel",0)->setComment("Level of the debugging print statements ");
  desc.addUntracked("HasMCMatch",false)
    ->setComment("Boolean for having the TrackToGenAssoc collection, only new sample have it");
  desc.addUntracked("DoTriggering",true)->setComment("Boolean to eecide whether we want to use triggers");
  desc.addUntracked("CalcSystematics",true)->setComment("Boolean to decide  whether we want to calculate the systematics");
  desc.addUntracked("GlobalMaxEta",2.1)->setComment("Cut on inner tracker track eta");
  desc.addUntracked("GlobalMinPt",55.0)->setComment("Cut on pT    at PRE-SELECTION");
  desc.addUntracked("GlobalMinNOH",8)->setComment("Cut on number of (valid) track pixel+strip hits");
  desc.addUntracked("GlobalMinNOPH",2)->setComment("Cut on number of (valid) track pixel hits");
  desc.addUntracked("GlobalMinFOVH",0.8)->setComment("Cut on fraction of valid track hits");
  desc.addUntracked("GlobalMinNOM",6)->setComment("Cut on number of dEdx hits (generally equal to #strip+#pixel-#ClusterCleaned hits)");
  desc.addUntracked("GlobalMaxChi2",5.0)->setComment("Cut on Track maximal Chi2/NDF");
  desc.addUntracked("GlobalMaxEIsol",2.0)->setComment("Cut on calorimeter isolation (E/P)");
  desc.addUntracked("GlobalMaxDZ",0.1)->setComment("Cut on 1D distance (cm) to closest vertex in Z direction");
  desc.addUntracked("GlobalMaxDXY",0.02)->setComment("Cut on 2D distance (cm) to closest vertex in R direction");
  desc.addUntracked("GlobalMaxPtErr",0.25)->setComment("Cut on error on track pT measurement");
  desc.addUntracked("GlobalMaxTIsol",50.0)->setComment("Cut on tracker isolation (SumPt)");
  desc.addUntracked("GlobalMiniRelIsoChg",0.1)->setComment("Cut on the charged PF based mini-isolation");
  desc.addUntracked("GlobalMassT",50.0)->setComment("Cut on the transverse mass");
  desc.addUntracked("GlobalMinIh",0.0)->setComment("Cut on dEdx estimator (Im,Ih,etc)");
  desc.addUntracked("TrackProbQCut",1.0)->setComment("Cut for probQ, 1.0 means no cuts applied");
  desc.addUntracked("GlobalMinIs",0.0)->setComment("Cut on dEdx discriminator (Is,Ias,etc)");
  desc.addUntracked("MinMuStations",2)->setComment("Minimum number of muon stations");
//  desc.addUntracked("GlobalMinNDOF",8.0)->setComment("Cut on number of DegreeOfFreedom used for muon TOF measurement");
//  desc.addUntracked("GlobalMinNDOFDT",6.0)->setComment("Cut on number of DT DegreeOfFreedom used for muon TOF measurement");
//  desc.addUntracked("GlobalMinNDOFCSC",6.0)->setComment("Cut on number of CSC DegreeOfFreedom used for muon TOF measurement");
//  desc.addUntracked("GlobalMaxTOFErr",0.15)->setComment("Cut on error on muon TOF measurement");
//  desc.addUntracked("globalMinTOF_",1.0)->setComment("Cut on minimal TOF");
  /*
   float  = 50;          // c
   */

 descriptions.add("HSCParticleAnalyzer",desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);

//=============================================================
//
//     Method for initializing pT and Is cuts
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
  edm::LogInfo(MOD) << CutPt.size() << " Different Final Selection will be tested\n"
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
//
//     Method for scaling eta per bin
//
//=============================================================
float Analyzer::scaleFactor(float eta) {
  float etaBins[15] = {-2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0, 0.3, 0.6, 0.9, 1.2, 1.5, 1.8, 2.1};
  float scaleBins[15] = {0, 0.97, 1.06, 1.00, 0.89, 0.91, 0.93, 0.93, 0.92, 0.92, 0.91, 0.89, 1.00, 1.06, 0.99};
  for (int i = 0; i < 15; i++)
    if (eta < etaBins[i])
      return scaleBins[i];
  return 0;
}

//=============================================================
//
//     Method for rescaling pT
//
//=============================================================
float Analyzer::RescaledPt(const float& pt, const float& eta, const float& phi, const int& charge) {
  if (typeMode_ != 3) {
    float newInvPt = 1 / pt + 0.000236 - 0.000135 * pow(eta, 2) + charge * 0.000282 * TMath::Sin(phi - 1.337);
    return 1 / newInvPt;
  } else {
    float newInvPt = (1. / pt) * 1.1;
    return 1 / newInvPt;
  }
}

//=============================================================
//
//     Method to get hit position
//
//=============================================================
TVector3 Analyzer::getOuterHitPos(const reco::DeDxHitInfo* dedxHits) {
  TVector3 point(0, 0, 0);
  if (!dedxHits)
    return point;
  //WAIT//float outerDistance=-1;
  for (unsigned int h = 0; h < dedxHits->size(); h++) {
    DetId detid(dedxHits->detId(h));
    //WAIT//moduleGeom* geomDet = moduleGeom::get(detid.rawId());
    //WAIT//TVector3 hitPos = geomDet->toGlobal(TVector3(dedxHits->pos(h).x(), dedxHits->pos(h).y(), dedxHits->pos(h).z()));
    //WAIT//if(hitPos.Mag()>outerDistance){outerDistance=hitPos.Mag();  point=hitPos;}
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
bool Analyzer::passPreselection(const reco::TrackRef track,
                                const reco::DeDxHitInfo* dedxHits,
                                const reco::DeDxData* dedxSObj,
                                const reco::DeDxData* dedxMObj,
                                const reco::MuonTimeExtra* tof,
                                const edm::Event& iEvent,
                                const float pixelProbs[],
                                const float Event_Weight,
                                Tuple* tuple,
                                const float GenBeta,
                                const bool RescaleP,
                                const float RescaleI,
                                const float RescaleT,
                                float MassErr,
                                const bool Ih_Iso_cut,
                                const unsigned int closestBackgroundPDGsIDs[]) {
  using namespace edm;
    
  //===================== Handle For vertex ===============
  vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken_);
  //===================== Handle For PFCandidate ===================
  const edm::Handle<reco::PFCandidateCollection> pfCandHandle = iEvent.getHandle(pfCandToken_);
    //===================== Handle For PFMET ===================
  const edm::Handle<std::vector<reco::PFMET>> pfMETHandle = iEvent.getHandle(pfMETToken_);
  
  
  if (vertexColl.size() < 1) {
    LogError(MOD) << "        >> Preselection not passed: there is no vertex"
                  << " -- this should never happen as there was a check before";
    return false;
  }
  
  // This is a repeated code here
  int highestPtGoodVertex = -1;
  int goodVerts = 0;
  float dzMin = 10000;
    // Loop on the vertices in the event
  for (unsigned int i = 0; i < vertexColl.size(); i++) {
    if (vertexColl[i].isFake() || fabs(vertexColl[i].z()) > 24 || vertexColl[i].position().rho() > 2 ||
        vertexColl[i].ndof() <= 4)
      continue;  //only consider good vertex
    goodVerts++;
    if (tuple) {
      tuple->BS_dzAll->Fill(track->dz(vertexColl[i].position()), Event_Weight);
      tuple->BS_dxyAll->Fill(track->dxy(vertexColl[i].position()), Event_Weight);
    }
    
    if (fabs(track->dz(vertexColl[i].position())) < fabs(dzMin)) {
      dzMin = fabs(track->dz(vertexColl[i].position()));
      highestPtGoodVertex = i;
    }
  } // End loop on the vertices in the event
  
  if (highestPtGoodVertex < 0) {
    highestPtGoodVertex = 0;
  }
  
  // Impact paramters dz and dxy
  float dz = track->dz(vertexColl[highestPtGoodVertex].position());
  float dxy = track->dxy(vertexColl[highestPtGoodVertex].position());
  
  // Loop on PF candidates
  bool pf_isPhoton = false, pf_isElectron = false, pf_isMuon = false;
  bool pf_isChHadron = false, pf_isNeutHadron = false, pf_isUndefined = false;
  float track_PFMiniIso_sumCharHadPt = 0, track_PFMiniIso_sumNeutHadPt = 0, track_PFMiniIso_sumPhotonPt = 0, track_PFMiniIso_sumPUPt = 0;
    
  // number of tracks as the first bin
  if (tuple) {
    tuple->pfType->Fill(0.5, EventWeight_);
  }
  
  
  if(pfCandHandle.isValid() && !pfCandHandle->empty()) {
    const reco::PFCandidateCollection* pf = pfCandHandle.product();
    for (unsigned int i = 0; i < pf->size(); i++){
    // https://github.com/cms-sw/cmssw/blob/72d0fc00976da53d1fb745eb7f37b2a4ad965d7e/
    // PhysicsTools/PatAlgos/plugins/PATIsolatedTrackProducer.cc#L555
      const reco::PFCandidate* pfCand = &(*pf)[i];

      pf_isElectron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::e;
      pf_isMuon = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;
      pf_isPhoton = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::gamma;

      pf_isChHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h;
      pf_isNeutHadron = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::h0;
      pf_isUndefined = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::X;

      if (pfCand->trackRef().isNonnull() && pfCand->trackRef().key() == track.key()) {
        if (tuple) {
          // Number of PF tracks matched to general track
            tuple->pfType->Fill(1.5, EventWeight_);
          if (pf_isElectron) {
            tuple->pfType->Fill(2.5, EventWeight_);
          } else if (pf_isMuon) {
            tuple->pfType->Fill(3.5, EventWeight_);
          } else if (pf_isPhoton) {
            tuple->pfType->Fill(4.5, EventWeight_);
          } else if (pf_isChHadron) {
           tuple->pfType->Fill(5.5, EventWeight_);
          } else if (pf_isNeutHadron) {
            tuple->pfType->Fill(6.5, EventWeight_);
          } else if (pf_isUndefined) {
            tuple->pfType->Fill(7.5, EventWeight_);
          } else {
           tuple->pfType->Fill(8.5, EventWeight_);
          }
        }
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
        // charged cands from PV get added to trackIso
        if(pf_isChHadron && fromPV) track_PFMiniIso_sumCharHadPt+=pt;
        // charged cands not from PV get added to pileup iso
        else if(pf_isChHadron) track_PFMiniIso_sumPUPt+=pt;
        // neutral hadron iso
        if(pf_isNeutHadron) track_PFMiniIso_sumNeutHadPt+=pt;
        // photon iso
        if(pf_isPhoton) track_PFMiniIso_sumPhotonPt+=pt;
      }
    }
  }//end loop PFCandidates
  
  // Calculate PF mini relative isolation
  float miniRelIsoAll = (track_PFMiniIso_sumCharHadPt + track_PFMiniIso_sumPUPt + track_PFMiniIso_sumNeutHadPt)/track->pt();
  float miniRelIsoChg = track_PFMiniIso_sumCharHadPt/track->pt();


  if (miniRelIsoChg > 2.0 && track->pt() > 55.0) {
     cout << "miniRelIsoChg > 2, it's " << miniRelIsoChg << " and pt is " << track->pt() << endl;
     cout << "miniRelIsoAll " << miniRelIsoAll << " and track_PFMiniIso_sumPUPt " << track_PFMiniIso_sumPUPt << endl;
  }
  
  // Calculate transverse mass
  float RecoPFMET_et = -1, RecoPFMET_phi = -1;
//  float RecoPFMET_eta = -1, RecoPFMET_significance = -1;

  if (pfMETHandle.isValid() && !pfMETHandle->empty()) {
    for (unsigned int i = 0; i < pfMETHandle->size(); i++) {
      const reco::PFMET* pfMet = &(*pfMETHandle)[i];
      RecoPFMET_et = pfMet->et();
//      RecoPFMET_eta = pfMet->eta();
      RecoPFMET_phi = pfMet->phi();
//      RecoPFMET_significance = pfMet->significance();
    }
  }
  float massT = sqrt(2*track->pt()*RecoPFMET_et*(1-cos(track->phi()-RecoPFMET_phi)));

  // Number of DeDx hits
  unsigned int numDeDxHits = (dedxSObj) ? dedxSObj->numberOfMeasurements() : 0;
  unsigned int missingHitsTillLast =
    track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
    track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  float validFractionTillLast =
    track->found() <= 0 ? -1 : track->found() / float(track->found() + missingHitsTillLast);
  
  float probQonTrack = pixelProbs[0];
  float probQonTrackNoLayer1 = pixelProbs[1];
  float probXYonTrack = pixelProbs[2];
  float probXYonTrackNoLayer1 = pixelProbs[3];
  
  // TODO: what do PUA and PUB stand for??
  bool PUA = (vertexColl.size() < 15);
  bool PUB = (vertexColl.size() >= 15);
  
  const edm::ValueMap<susybsm::HSCPIsolation> IsolationMap = iEvent.get(hscpIsoToken_);
  susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
  
  float EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy()) / track->p();
  float IsoTK_SumEt = (Ih_Iso_cut) ? hscpIso.Get_TK_SumEt() : 0.0;
  
  float Ih = (dedxMObj) ?  dedxMObj->dEdx() : 0.0;
  float Is = (dedxSObj) ? dedxSObj->dEdx() : 0.0;
  
  //Find distance to nearest segment on opposite side of detector
  float minPhi = 0.0, minEta = 0.0;
  float segSep = SegSep(track, iEvent, minPhi, minEta);

  // Preselection cuts
  bool passedCutsArray[22];
  std::fill(std::begin(passedCutsArray), std::end(passedCutsArray),false);
  
  // No cut, i.e. events after trigger
  passedCutsArray[0]  = true;
  // Check if eta is inside the max eta cut
  passedCutsArray[1]  = (fabs(track->eta()) < globalMaxEta_) ? true : false;
  // Cut on number of matched muon stations
  passedCutsArray[2]  = (track->pt() > globalMinPt_) ? true : false;
  // Check the number of found hits (measurements)
  passedCutsArray[3]  = (typeMode_ != 3 && track->found() > globalMinNOH_) ? true : false;
  // Check the number of pixel hits
  passedCutsArray[4]  = (typeMode_ != 3 && fabs(track->hitPattern().numberOfValidPixelHits()) > globalMinNOPH_) ? true : false;
  // Check the min fraction of valid hits
  passedCutsArray[5]  = (typeMode_ != 3 && track->validFraction() > globalMinFOVH_) ? true : false;
  // Cut for the number of dEdx hits
  passedCutsArray[6]  = (numDeDxHits > globalMinNOM_)  ? true : false;
  // This should be revised, for now switching it off
  passedCutsArray[7]  = (probXYonTrack >= 0.0 || probXYonTrack <= 1.0)  ? true : false;
  // Select only high purity tracks
  passedCutsArray[8]  = (typeMode_ != 3 && track->quality(reco::TrackBase::highPurity)) ? true : false;
  // Cut on the chi2 / ndof
  passedCutsArray[9] = (typeMode_ != 3 && track->chi2() / track->ndof() < globalMaxChi2_) ? true : false;
  // Cut on the energy over momenta
  passedCutsArray[10] = (EoP < globalMaxEIsol_) ? true : false;
//  passedCutsArray[10] = true;
  // Cut on the impact parameter
    // for typeMode_ 5 dz is supposed to come from the beamspot, TODO
  passedCutsArray[11] = (  (typeMode_ != 5 && fabs(dz) < globalMaxDZ_)
                        || (typeMode_ == 5 && fabs(dz) < 4)) ? true : false;
  // for typeMode_ 5 dxy is supposed to come from the beamspot, TODO
  passedCutsArray[12] = (  (typeMode_ != 5 && fabs(dxy) < globalMaxDXY_)
                        || (typeMode_ == 5 && fabs(dxy) < 4)) ? true : false;
  // Cut on the uncertainty of the pt measurement
  passedCutsArray[13] = (typeMode_ != 3 && (track->ptError() / track->pt()) < globalMaxPtErr_) ? true : false;
  // Cut on the tracker based isolation
//  passedCutsArray[14] = true;
  passedCutsArray[14] = ( IsoTK_SumEt < globalMaxTIsol_) ? true : false;
  // Cut on the PF based mini-isolation
  passedCutsArray[15] = ( miniRelIsoChg < globalMiniRelIsoChg_) ? true : false;
  // Cut on the transverse mass
  passedCutsArray[16] = true; //( massT < globalMassT_) ? true : false;
  // Cut on min Ih (or max for fractionally charged)
  passedCutsArray[17] = (  (typeMode_ != 5 &&  Ih > globalMinIh_)
                        || (typeMode_ == 5 && Ih < globalMinIh_)) ? true : false;
    // Cut away background events based on the probQ
  passedCutsArray[18]  = (probQonTrack < trackProbQCut_) ? true : false;
// passedCutsArray[18]  = (probQonTrack < trackProbQCut_ || probQonTrackNoLayer1 < trackProbQCut_) ? true : false;
  // TOF only cuts
  passedCutsArray[19] = (typeMode_ != 3 || (typeMode_ == 3 && muonStations(track->hitPattern()) > minMuStations_)) ? true : false;
  passedCutsArray[20] = (typeMode_ != 3 || (typeMode_ == 3 && fabs(track->phi()) > 1.2 && fabs(track->phi()) < 1.9)) ? true : false;
  passedCutsArray[21] = (typeMode_ != 3 || (typeMode_ == 3 && fabs(minEta) > minSegEtaSep)) ? true : false;
  
  // Not used cuts TODO: revise
  // cut on the number of missing hits from IP till last hit (excluding hits behind the last hit)
  // bool cutMinNumOfMissingHits = (typeMode_ != 3 && missingHitsTillLast > GlobalMaxNOMHTillLast) ? true : false;
  // cut on the fraction of valid hits divided by total expected hits until the last one
  // bool cutMinFractOfValidHitsTillLast = (typeMode_ != 3 && validFractionTillLast < GlobalMinFOVHTillLast) ? true : false;
  // cut on relative tracker isolation (SumPt/Pt)
  // bool cutRelTKIso = ( IsoTK_SumEt / track->pt() > GlobalMaxRelTIsol)  ? true : false;
  // Cut for number of DOF in TOF ana
  
  // CutFlow in a single plot
  if (tuple) {
    for (size_t i=0;i<sizeof(passedCutsArray);i++) {
      bool allCutsPassedSoFar = true;
      for (size_t j=0;j<=i;j++) {
        if (!passedCutsArray[j]) {
          allCutsPassedSoFar = false;
        }
      }
      if (allCutsPassedSoFar) {
        tuple->CutFlow->Fill((i+0.5), Event_Weight);
      }
    }
  }
    
    // Preselection cuts when probQ is one of the first cuts
    bool passedCutsArray2[22];
    std::fill(std::begin(passedCutsArray2), std::end(passedCutsArray2),false);
    passedCutsArray2[0]  = true; // passed trigger
    passedCutsArray2[1]  = (fabs(track->eta()) < globalMaxEta_) ? true : false;
    passedCutsArray2[2]  = (track->pt() > globalMinPt_) ? true : false;
    passedCutsArray2[3]  = (probQonTrack < trackProbQCut_) ? true : false;
//  passedCutsArray2[3]  = (probQonTrack < trackProbQCut_ || probQonTrackNoLayer1 < trackProbQCut_) ? true : false;
    passedCutsArray2[4]  = (typeMode_ != 3 && track->found() > globalMinNOH_) ? true : false;
    passedCutsArray2[5]  = (typeMode_ != 3 && fabs(track->hitPattern().numberOfValidPixelHits()) > globalMinNOPH_) ? true : false;
    passedCutsArray2[6]  = (typeMode_ != 3 && track->validFraction() > globalMinFOVH_) ? true : false;
    passedCutsArray2[7]  = (numDeDxHits > globalMinNOM_)  ? true : false;
    passedCutsArray2[8]  = (probXYonTrack >= 0.0 || probXYonTrack <= 1.0)  ? true : false;
    passedCutsArray2[9]  = (typeMode_ != 3 && track->quality(reco::TrackBase::highPurity)) ? true : false;
    passedCutsArray2[10] = (typeMode_ != 3 && track->chi2() / track->ndof() < globalMaxChi2_) ? true : false;
    passedCutsArray2[11] = (EoP < globalMaxEIsol_) ? true : false;
//    passedCutsArray2[11] = true;
    // for typeMode_ 5 dz is supposed to come from the beamspot, TODO
    passedCutsArray2[12] = (   (typeMode_ != 5 && fabs(dz) < globalMaxDZ_)
                            || (typeMode_ == 5 && fabs(dz) < 4.0)) ? true : false;
    // for typeMode_ 5 dxy is supposed to come from the beamspot, TODO
    passedCutsArray2[13] = (  (typeMode_ != 5 && fabs(dxy) < globalMaxDXY_)
                           || (typeMode_ == 5 && fabs(dxy)) < 4.0) ? true : false;
    passedCutsArray2[14] = (typeMode_ != 3 && (track->ptError() / track->pt()) < globalMaxPtErr_) ? true : false;
//    passedCutsArray2[15] = true;
    passedCutsArray2[15] = ( IsoTK_SumEt < globalMaxTIsol_) ? true : false;
    // Cut on the PF based mini-isolation
    passedCutsArray2[16] = ( miniRelIsoChg < globalMiniRelIsoChg_) ? true : false;
    // Cut on the transverse mass
    passedCutsArray2[17] = true;// ( massT < globalMassT_) ? true : false;
    // Cut on Ih
    passedCutsArray2[18] = (  (typeMode_ != 5 && Ih > globalMinIh_)
                           || (typeMode_ == 5 && Ih < globalMinIh_)) ? true : false;
    // TOF only cuts
    passedCutsArray2[19] = (typeMode_ != 3 || (typeMode_ == 3 && muonStations(track->hitPattern()) > minMuStations_)) ? true : false;
    passedCutsArray2[20] = (typeMode_ != 3 || (typeMode_ == 3 && fabs(track->phi()) > 1.2 && fabs(track->phi()) < 1.9)) ? true : false;
    passedCutsArray2[21] = (typeMode_ != 3 || (typeMode_ == 3 && fabs(minEta) > minSegEtaSep)) ? true : false;
    
    // CutFlow in a single plot when probQ is one of the first cuts
    if (tuple) {
      for (size_t i=0;i<sizeof(passedCutsArray2);i++) {
        bool allCutsPassedSoFar = true;
        for (size_t j=0;j<=i;j++) {
          if (!passedCutsArray2[j]) {
            allCutsPassedSoFar = false;
          }
        }
        if (allCutsPassedSoFar) {
          tuple->CutFlowProbQFirst->Fill((i+0.5), Event_Weight);
        }
      }
    }
  
  // Before (pre)selection plots
  if (tuple) {
    if (GenBeta >= 0) {
      tuple->Beta_Matched->Fill(GenBeta, Event_Weight);
    }
    tuple->BS_Eta->Fill(track->eta(), Event_Weight);
    tuple->BS_MatchedStations->Fill(muonStations(track->hitPattern()), Event_Weight);
    tuple->BS_NVertex->Fill(vertexColl.size(), Event_Weight);
    tuple->BS_NVertex_NoEventWeight->Fill(vertexColl.size());
    tuple->BS_TNOH->Fill(track->found(), Event_Weight);
    if (PUA) {
      tuple->BS_TNOH_PUA->Fill(track->found(), Event_Weight);
      tuple->BS_TNOM_PUA->Fill(numDeDxHits, Event_Weight);
    }
    if (PUB) {
      tuple->BS_TNOH_PUB->Fill(track->found(), Event_Weight);
      tuple->BS_TNOM_PUB->Fill(numDeDxHits, Event_Weight);
    }
    tuple->BS_TNOHFraction->Fill(track->validFraction(), Event_Weight);
    tuple->BS_TNOPH->Fill(track->hitPattern().numberOfValidPixelHits(), Event_Weight);
    tuple->BS_TNOHFractionTillLast->Fill(validFractionTillLast, Event_Weight);
    tuple->BS_TNOMHTillLast->Fill(missingHitsTillLast, Event_Weight);
    tuple->BS_TNOM->Fill(numDeDxHits, Event_Weight);
    if (track->found() - numDeDxHits) {
      tuple->BS_EtaNBH->Fill(track->eta(), track->found() - numDeDxHits, Event_Weight);
    }
    tuple->BS_ProbQ->Fill(probQonTrack, EventWeight_);
    tuple->BS_ProbXY->Fill(probXYonTrack, EventWeight_);
    tuple->BS_ProbQNoL1->Fill(probQonTrackNoLayer1, EventWeight_);
    tuple->BS_ProbXYNoL1->Fill(probXYonTrackNoLayer1, EventWeight_);
    if (tof) {
      tuple->BS_nDof->Fill(tof->nDof(), Event_Weight);
      tuple->BS_MTOF->Fill(tof->inverseBeta(), Event_Weight);
      tuple->BS_TOFError->Fill(tof->inverseBetaErr(), Event_Weight);
      tuple->BS_TimeAtIP->Fill(tof->timeAtIpInOut(), Event_Weight);
    }
    tuple->BS_Qual->Fill(track->qualityMask(), Event_Weight);
    tuple->BS_Chi2PerNdof->Fill(track->chi2() / track->ndof(), Event_Weight);
    tuple->BS_MPt->Fill(track->pt(), Event_Weight);
    tuple->BS_NOMoNOHvsPV->Fill(goodVerts, numDeDxHits / (float)track->found(), Event_Weight);
    tuple->BS_dzMinv3d->Fill(dz, Event_Weight);
    tuple->BS_dxyMinv3d->Fill(dxy, Event_Weight);
    tuple->BS_Dxy->Fill(dxy, Event_Weight);
    tuple->BS_Dz->Fill(dz, Event_Weight);
    tuple->BS_EtaDz->Fill(track->eta(), dz, Event_Weight);
    tuple->BS_PV->Fill(goodVerts, Event_Weight);
    tuple->BS_PV_NoEventWeight->Fill(goodVerts);
    tuple->BS_EIsol->Fill(EoP, Event_Weight);
    tuple->BS_SumpTOverpT->Fill(IsoTK_SumEt / track->pt(), Event_Weight);
    tuple->BS_PtErrOverPt->Fill(track->ptError() / track->pt(), Event_Weight);
    tuple->BS_PtErrOverPt2->Fill(track->ptError() / (track->pt()*track->pt()), Event_Weight);
    tuple->BS_PtErrOverPtVsPtErrOverPt2->Fill(track->ptError() / track->pt(),track->ptError() / (track->pt()*track->pt()), Event_Weight);
    tuple->BS_PtErrOverPtVsPt->Fill(track->ptError() / track->pt(), track->pt(), Event_Weight);
    tuple->BS_TIsol->Fill(IsoTK_SumEt, Event_Weight);
    tuple->BS_MIh->Fill(Ih, Event_Weight);
    tuple->BS_MIs->Fill(Is, Event_Weight);
    tuple->BS_massT->Fill(massT, Event_Weight);
    // Add PFCadidate based isolation info to the tuple
    // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/
    // PhysicsTools/NanoAOD/plugins/IsoValueMapProducer.cc#L157
    tuple->BS_MiniRelIsoAll->Fill(miniRelIsoAll, Event_Weight);
    tuple->BS_MiniRelIsoChg->Fill(miniRelIsoChg, Event_Weight);
    tuple->BS_SegSep->Fill(segSep, Event_Weight);
    tuple->BS_SegMinPhiSep->Fill(minPhi, Event_Weight);
    tuple->BS_SegMinEtaSep->Fill(minEta, Event_Weight);
    tuple->BS_OpenAngle->Fill(OpenAngle, Event_Weight);
    tuple->BS_MassErr->Fill(MassErr, Event_Weight);
  }
  
  // N-1 plots
  if (tuple) {
    for (size_t i=0;i<sizeof(passedCutsArray);i++) {
      bool allOtherCutsPassed = true;
      for (size_t j=0;j<sizeof(passedCutsArray);j++) {
        if (i==j) continue;
        if (!passedCutsArray[j]) {
          allOtherCutsPassed = false;
          // We found a cut that's not passed, no point in looking into the rest of them
          break;
        }
      }
      if (allOtherCutsPassed) {
        if (i==1)  { tuple->N1Eta->Fill(track->eta(), Event_Weight); };
        if (i==2)  { tuple->N1MPt->Fill(track->pt(), Event_Weight); };
        if (i==3)  { tuple->N1TNOH->Fill(track->found(), Event_Weight); };
        if (i==4)  { tuple->N1TNOPH->Fill(track->hitPattern().numberOfValidPixelHits(), Event_Weight); };
        if (i==5)  { tuple->N1TNOHFraction->Fill(track->validFraction(), Event_Weight); };
        if (i==6)  { tuple->N1TNOM->Fill(numDeDxHits, Event_Weight); };
        if (i==7)  { tuple->N1ProbQ->Fill(probQonTrack, EventWeight_); };
        if (i==8)  { tuple->N1ProbXY->Fill(probXYonTrack, EventWeight_); };
        if (i==9)  { tuple->N1Qual->Fill(track->qualityMask(), Event_Weight); };
        if (i==10) { tuple->N1Chi2PerNdof->Fill(track->chi2() / track->ndof(), Event_Weight); };
        if (i==11) { tuple->N1EIsol->Fill(EoP, Event_Weight); };
        if (i==12) { tuple->N1Dz->Fill(dz, Event_Weight); };
        if (i==13) { tuple->N1Dxy->Fill(dxy, Event_Weight); };
        if (i==14) { tuple->N1PtErrOverPt->Fill(track->ptError() / track->pt(), Event_Weight); };
        if (i==15) { tuple->N1SumpTOverpT->Fill(IsoTK_SumEt / track->pt(), Event_Weight); };
        if (i==16) { tuple->N1MIh->Fill(Ih, Event_Weight); };
        if (i==17) { tuple->N1Stations->Fill(muonStations(track->hitPattern()), Event_Weight); };
        if (i==18) { LogDebug(MOD) << "cutPhiTOFOnly"; };
        if (i==19) { LogDebug(MOD) << "cutEtaTOFOnly"; };
      }
    }
  }

  // Return false in the function is a given cut is not passed
  for (size_t i=0;i<sizeof(passedCutsArray);i++) {
    if (!passedCutsArray[i]) {
      if (debug_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed for the " <<  std::to_string(i) << "-th cut, please check the code what that corresponds to";
      // TODO: when the preselection list finalizes I might be more verbose than this
      // Plot Eta after each cut
      if (tuple) {
        tuple->CutFlowEta->Fill(track->eta(),i+0.5, Event_Weight);
      }
      return false;
    }
  }
  
  // After (pre)selection plots
  if (tuple) {
    tuple->PostPreS_Eta->Fill(track->eta(), Event_Weight);
    tuple->PostPreS_MatchedStations->Fill(muonStations(track->hitPattern()), Event_Weight);
    tuple->PostPreS_NVertex->Fill(vertexColl.size(), Event_Weight);
    tuple->PostPreS_NVertex_NoEventWeight->Fill(vertexColl.size());
    tuple->PostPreS_TNOH->Fill(track->found(), Event_Weight);
    if (PUA) {
      tuple->PostPreS_TNOH_PUA->Fill(track->found(), Event_Weight);
      tuple->PostPreS_TNOM_PUA->Fill(numDeDxHits, Event_Weight);
    }
    if (PUB) {
      tuple->PostPreS_TNOH_PUB->Fill(track->found(), Event_Weight);
      tuple->PostPreS_TNOM_PUB->Fill(numDeDxHits, Event_Weight);
    }
    tuple->PostPreS_TNOHFraction->Fill(track->validFraction(), Event_Weight);
    tuple->PostPreS_TNOPH->Fill(track->hitPattern().numberOfValidPixelHits(), Event_Weight);
    tuple->PostPreS_TNOHFractionTillLast->Fill(validFractionTillLast, Event_Weight);
    tuple->PostPreS_TNOMHTillLast->Fill(missingHitsTillLast, Event_Weight);
    tuple->PostPreS_TNOM->Fill(numDeDxHits, Event_Weight);
    tuple->PostPreS_ProbQ->Fill(probQonTrack, EventWeight_);
    tuple->PostPreS_ProbXY->Fill(probXYonTrack, EventWeight_);
    tuple->PostPreS_ProbQNoL1->Fill(probQonTrackNoLayer1, EventWeight_);
    tuple->PostPreS_ProbXYNoL1->Fill(probXYonTrackNoLayer1, EventWeight_);
    if (tof) {
      tuple->PostPreS_nDof->Fill(tof->nDof(), Event_Weight);
      tuple->PostPreS_MTOF->Fill(tof->inverseBeta(), Event_Weight);
      tuple->PostPreS_TOFError->Fill(tof->inverseBetaErr(), Event_Weight);
      tuple->PostPreS_TimeAtIP->Fill(tof->timeAtIpInOut(), Event_Weight);
    }
    tuple->PostPreS_Qual->Fill(track->qualityMask(), Event_Weight);
    tuple->PostPreS_Chi2PerNdof->Fill(track->chi2() / track->ndof(), Event_Weight);
    tuple->PostPreS_Pt->Fill(track->pt(), Event_Weight);
    tuple->PostPreS_NOMoNOHvsPV->Fill(goodVerts, numDeDxHits / (float)track->found(), Event_Weight);
    tuple->PostPreS_dzMinv3d->Fill(dz, Event_Weight);
    tuple->PostPreS_dxyMinv3d->Fill(dxy, Event_Weight);
    tuple->PostPreS_PV->Fill(goodVerts, Event_Weight);
    tuple->PostPreS_PV_NoEventWeight->Fill(goodVerts);
    
    tuple->PostPreS_EIsol->Fill(EoP, Event_Weight);
    tuple->PostPreS_SumpTOverpT->Fill(IsoTK_SumEt / track->pt(), Event_Weight);
    tuple->PostPreS_PtErrOverPt->Fill(track->ptError() / track->pt(), Event_Weight);
    tuple->PostPreS_PtErrOverPt2->Fill(track->ptError() / (track->pt()*track->pt()), Event_Weight);
    tuple->PostPreS_PtErrOverPtVsPtErrOverPt2->Fill(track->ptError() / track->pt(),track->ptError() / (track->pt()*track->pt()), Event_Weight);
    tuple->PostPreS_PtErrOverPtVsPt->Fill(track->ptError() / track->pt(), track->pt(), Event_Weight);
    tuple->PostPreS_TIsol->Fill(IsoTK_SumEt, Event_Weight);
    tuple->PostPreS_MIh->Fill(Ih, Event_Weight);
    tuple->PostPreS_MIh_NoEventWeight->Fill(Ih);
    tuple->PostPreS_MIs->Fill(Is, Event_Weight);
    tuple->PostPreS_MIs_NoEventWeight->Fill(Is);
    tuple->PostPreS_massT->Fill(massT, Event_Weight);
      // Add PFCadidate based isolation info to the tuple
      // https://github.com/cms-sw/cmssw/blob/6d2f66057131baacc2fcbdd203588c41c885b42c/
      // PhysicsTools/NanoAOD/plugins/IsoValueMapProducer.cc#L157
    tuple->PostPreS_MiniRelIsoAll->Fill(miniRelIsoAll, Event_Weight);
    tuple->PostPreS_MiniRelIsoChg->Fill(miniRelIsoChg, Event_Weight);
    
    tuple->PostPreS_EtaPerGenID->Fill(track->eta(), closestBackgroundPDGsIDs[0], Event_Weight);
    tuple->PostPreS_ProbQPerGenID->Fill(probQonTrack, closestBackgroundPDGsIDs[0], EventWeight_);
    tuple->PostPreS_ProbXYPerGenID->Fill(probXYonTrack, closestBackgroundPDGsIDs[0], EventWeight_);
    tuple->PostPreS_PtPerGenID->Fill(track->pt(), closestBackgroundPDGsIDs[0], Event_Weight);
    tuple->PostPreS_EIsolPerGenID->Fill(EoP, closestBackgroundPDGsIDs[0], Event_Weight);
    tuple->PostPreS_MIhPerGenID->Fill(Ih, closestBackgroundPDGsIDs[0], Event_Weight);
    tuple->PostPreS_MIsPerGenID->Fill(Is, closestBackgroundPDGsIDs[0], Event_Weight);
    tuple->PostPreS_massTPerGenID->Fill(massT, closestBackgroundPDGsIDs[0], Event_Weight);
    tuple->PostPreS_EtaPerMomGenID->Fill(track->eta(), closestBackgroundPDGsIDs[1], Event_Weight);
    tuple->PostPreS_ProbQPerMomGenID->Fill(probQonTrack, closestBackgroundPDGsIDs[1], EventWeight_);
    tuple->PostPreS_ProbXYPerMomGenID->Fill(probXYonTrack, closestBackgroundPDGsIDs[1], EventWeight_);
    tuple->PostPreS_PtPerMomGenID->Fill(track->pt(), closestBackgroundPDGsIDs[1], Event_Weight);
    tuple->PostPreS_EIsolPerMomGenID->Fill(EoP, closestBackgroundPDGsIDs[1], Event_Weight);
    tuple->PostPreS_MIhPerMomGenID->Fill(Ih, closestBackgroundPDGsIDs[1], Event_Weight);
    tuple->PostPreS_MIsPerMomGenID->Fill(Is, closestBackgroundPDGsIDs[1], Event_Weight);
    tuple->PostPreS_massTPerMomGenID->Fill(massT, closestBackgroundPDGsIDs[1], Event_Weight);
    tuple->PostPreS_MassErr->Fill(MassErr, Event_Weight);
  }

  // Fill up gen based beta histo after preselection
  if (tuple && GenBeta >= 0) {
    tuple->Beta_PreselectedA->Fill(GenBeta, Event_Weight);
  }

  // Cut on  Rescaled P
  if (RescaleP && RescaledPt(track->pt(), track->eta(), track->phi(), track->charge()) < globalMinPt_) {
      return false;
  }

  // Cut on  Rescaled Is
  if (dedxSObj && RescaleI != 0.0) {
    if (dedxSObj->dEdx() + RescaleI < globalMinIs_) {
      if (debug_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Rescaled Ias is too low for fractionally charged";
      return false;
    } else {
    if (debug_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for rescaled Ias cut";
    }
  }

  if (tof) {
    if ((typeMode_ > 1 && typeMode_ != 5) && tof->inverseBetaErr() > GlobalMaxTOFErr) {
      return false;
    }

    if (typeMode_ == 3 && min(min(fabs(tof->timeAtIpInOut() - 100), fabs(tof->timeAtIpInOut() - 50)),
                              min(fabs(tof->timeAtIpInOut() + 100), fabs(tof->timeAtIpInOut() + 50))) < 5)
      return false;
  } // End condition on tof existence or not

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
    if (tuple) {
      tuple->BS_dR_NVTrack->Fill(minDr, Event_Weight);
    }
    if (minDr > 0.4) {
      return false;
    }
    if (tuple) {
      tuple->NVTrack->Fill(0.0, Event_Weight);
    }

    // Find displacement of tracks with respect to beam spot
    const reco::BeamSpot beamSpotColl = iEvent.get(offlineBeamSpotToken_);
    float dzFromBeamSpot = NVTrack.dz(beamSpotColl.position());
    float dxyFromBeamSpot = NVTrack.dxy(beamSpotColl.position());
    if (debug_ > 8 ) LogPrint(MOD) << dzFromBeamSpot << " and " << dxyFromBeamSpot;
    // TODO use this for TOF only analysis, instead of dxy and dz

    if (muonStations(NVTrack.hitPattern()) < minMuStations_)
      return false;
  
    if (tuple) {
      tuple->MTOF->Fill(0.0, Event_Weight);
      if (GenBeta >= 0)
        tuple->Beta_PreselectedB->Fill(GenBeta, Event_Weight);
    }
  } // End condition for TOF only analysis

  

//  //mk if(MassErr > 0 && MassErr > 2.2)return false; //FIXME jozze -- cut on relative mass error in units of 8*MassErr/Mass

  if (tuple) {
    //Plotting segment separation depending on whether track passed dz cut
    if (fabs(dz) > globalMaxDZ_) {
      tuple->BS_SegMinEtaSep_FailDz->Fill(minEta, Event_Weight);
    } else {
      tuple->BS_SegMinEtaSep_PassDz->Fill(minEta, Event_Weight);
    }
    //Plots for tracking failing Eta Sep cut
    if (fabs(minEta) < minSegEtaSep) {
      //Needed to compare dz distribution of cosmics in pure cosmic and main sample
      tuple->BS_Dz_FailSep->Fill(dz);
    }

    if (tof) {
      //Plots for tracks in dz control region
      if (fabs(dz) > CosmicMinDz && fabs(dz) < CosmicMaxDz) {
        tuple->BS_Pt_FailDz->Fill(track->pt(), Event_Weight);
        tuple->BS_TOF_FailDz->Fill(tof->inverseBeta(), Event_Weight);
        if (fabs(track->eta()) > CSCRegion) {
          tuple->BS_TOF_FailDz_CSC->Fill(tof->inverseBeta(), Event_Weight);
          tuple->BS_Pt_FailDz_CSC->Fill(track->pt(), Event_Weight);
        } else if (fabs(track->eta()) < DTRegion) {
          tuple->BS_TOF_FailDz_DT->Fill(tof->inverseBeta(), Event_Weight);
          tuple->BS_Pt_FailDz_DT->Fill(track->pt(), Event_Weight);
        }
      }
    //Plots of dz
      if (fabs(track->eta()) > CSCRegion) {
        tuple->BS_Dz_CSC->Fill(dz, Event_Weight);
      } else if (fabs(track->eta()) < DTRegion) {
        tuple->BS_Dz_DT->Fill(dz, Event_Weight);
      }
    }
  }

  TreeDXY = dxy;
  TreeDZ = dz;
  
  bool DXYSB = (typeMode_ == 5 && fabs(dxy) > globalMaxDXY_) ? true : false;
  bool DZSB = (typeMode_ == 5 && fabs(dz) > globalMaxDZ_) ? true : false;
  
//  if (cutEtaTOFOnly) {
//    if (debug_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: for TOF only analysis, eta is too low";
//    return false;
//  } else if (typeMode_ == 3 && fabs(minEta) > minSegEtaSep) {
//    if (debug_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for TOF eta cut";
//  }
//  if (tuple)
//    tuple->BS_Phi->Fill(track->phi(), Event_Weight);

//  if (cutPhiTOFOnly) {
//    if (debug_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: for TOF only analysis, 1.2 < phi < 1.9";
//    return false;
//  }
  
  //check if HSCP is compatible with cosmics.
  bool OASB = (typeMode_ == 5 && OpenAngle >= 2.8) ? true : false;

  isCosmicSB = DXYSB && DZSB && OASB;
  isSemiCosmicSB = (!isCosmicSB && (DXYSB || DZSB || OASB));

  if (tuple) {
    if (dedxSObj) {
      tuple->BS_EtaIs->Fill(track->eta(), dedxSObj->dEdx(), Event_Weight);
    } 
    if (dedxMObj) {
      tuple->BS_EtaIh->Fill(track->eta(), dedxMObj->dEdx(), Event_Weight);
    }
    tuple->BS_EtaP->Fill(track->eta(), track->p(), Event_Weight);
    tuple->BS_EtaPt->Fill(track->eta(), track->pt(), Event_Weight);
    if (tof)
      tuple->BS_EtaTOF->Fill(track->eta(), tof->inverseBeta(), Event_Weight);
  }

  if (tuple) {
    if (GenBeta >= 0)
      tuple->Beta_PreselectedC->Fill(GenBeta, Event_Weight);
    if (DZSB && OASB)
      tuple->BS_Dxy_Cosmic->Fill(dxy, Event_Weight);
    if (DXYSB && OASB)
      tuple->BS_Dz_Cosmic->Fill(dz, Event_Weight);
    if (DXYSB && DZSB)
      tuple->BS_OpenAngle_Cosmic->Fill(OpenAngle, Event_Weight);

    //WAIT//
    TVector3 outerHit = getOuterHitPos(dedxHits);
    TVector3 vertex(vertexColl[highestPtGoodVertex].position().x(),
                    vertexColl[highestPtGoodVertex].position().y(),
                    vertexColl[highestPtGoodVertex].position().z());
    tuple->BS_LastHitDXY->Fill((outerHit).Perp(), Event_Weight);
    tuple->BS_LastHitD3D->Fill((outerHit).Mag(), Event_Weight);

    tuple->BS_P->Fill(track->p(), Event_Weight);
    tuple->BS_Pt->Fill(track->pt(), Event_Weight);
    if (PUA)
      tuple->BS_Pt_PUA->Fill(track->pt(), Event_Weight);
    if (PUB)
      tuple->BS_Pt_PUB->Fill(track->pt(), Event_Weight);
    if (DXYSB && DZSB && OASB)
      tuple->BS_Pt_Cosmic->Fill(track->pt(), Event_Weight);

    if (fabs(track->eta()) < DTRegion)
      tuple->BS_Pt_DT->Fill(track->pt(), Event_Weight);
    else
      tuple->BS_Pt_CSC->Fill(track->pt(), Event_Weight);

    
//    float RecoQoPt = track->charge() / track->pt();
//    if (!hscp.trackRef().isNull() && hscp.trackRef()->pt() > 200) {
//      float InnerRecoQoPt = hscp.trackRef()->charge() / hscp.trackRef()->pt();
//      tuple->BS_InnerInvPtDiff->Fill((RecoQoPt - InnerRecoQoPt) / InnerRecoQoPt, Event_Weight);
//    }

    if (dedxSObj)
      tuple->BS_Is->Fill(dedxSObj->dEdx(), Event_Weight);
    if (dedxSObj && PUA)
      tuple->BS_Is_PUA->Fill(dedxSObj->dEdx(), Event_Weight);
    if (dedxSObj && PUB)
      tuple->BS_Is_PUB->Fill(dedxSObj->dEdx(), Event_Weight);
    if (dedxSObj && DXYSB && DZSB && OASB)
      tuple->BS_Is_Cosmic->Fill(dedxSObj->dEdx(), Event_Weight);
    if (dedxMObj)
      tuple->BS_Ih->Fill(dedxMObj->dEdx(), Event_Weight);
    if (dedxMObj && PUA)
      tuple->BS_Ih_PUA->Fill(dedxMObj->dEdx(), Event_Weight);
    if (dedxMObj && PUB)
      tuple->BS_Ih_PUB->Fill(dedxMObj->dEdx(), Event_Weight);

    if (tof) {
      tuple->BS_TOF->Fill(tof->inverseBeta(), Event_Weight);
      if (PUA)
        tuple->BS_TOF_PUA->Fill(tof->inverseBeta(), Event_Weight);
      if (PUB)
        tuple->BS_TOF_PUB->Fill(tof->inverseBeta(), Event_Weight);
      if (dttof->nDof() > 6)
        tuple->BS_TOF_DT->Fill(dttof->inverseBeta(), Event_Weight);
      if (csctof->nDof() > 6)
        tuple->BS_TOF_CSC->Fill(csctof->inverseBeta(), Event_Weight);
      tuple->BS_PtTOF->Fill(track->pt(), tof->inverseBeta(), Event_Weight);
    }
    if (dedxSObj && dedxMObj) {
      tuple->BS_PIs->Fill(track->p(), dedxSObj->dEdx(), Event_Weight);
      tuple->BS_IhIs->Fill(dedxMObj->dEdx(), dedxSObj->dEdx(), Event_Weight);
      tuple->BS_PIh->Fill(track->p(), dedxMObj->dEdx(), Event_Weight);
      tuple->BS_PtIs->Fill(track->pt(), dedxSObj->dEdx(), Event_Weight);
      tuple->BS_PtIh->Fill(track->pt(), dedxMObj->dEdx(), Event_Weight);
    }
    if (dedxSObj && tof)
      tuple->BS_TOFIs->Fill(tof->inverseBeta(), dedxSObj->dEdx(), Event_Weight);
    if (dedxMObj && tof)
      tuple->BS_TOFIh->Fill(tof->inverseBeta(), dedxMObj->dEdx(), Event_Weight);

    //Muon only prediction binned depending on where in the detector the track is and how many muon stations it has
    //Binning not used for other analyses
    int bin = -1;
    if (typeMode_ == 3) {
      if (fabs(track->eta()) < DTRegion)
        bin = muonStations(track->hitPattern()) - 2;
      else
        bin = muonStations(track->hitPattern()) + 1;
      tuple->BS_Pt_Binned[to_string(bin)]->Fill(track->pt(), Event_Weight);
    }
  }

  return true;
}

//=============================================================
//
//     Selection
//
//=============================================================
bool Analyzer::passSelection(const reco::TrackRef track,
                             const reco::DeDxData* dedxSObj,
                             const reco::DeDxData* dedxMObj,
                             const reco::MuonTimeExtra* tof,
                             const edm::Event& iEvent,
                             float Event_Weight,
                             const int& CutIndex,
                             Tuple*& tuple,
                             const bool isFlip,
                             const float GenBeta,
                             const bool RescaleP,
                             const float RescaleI,
                             const float RescaleT) {
  static constexpr const char* const MOD = "Analyzer";
  using namespace edm;
  float MuonTOF, Is, Ih;
  
  if (track.isNull()) {
    LogPrint(MOD) << "@passSelection: track.isNull() -- this should never happen!!!";
    return false;
  }
  
  tof ? MuonTOF = tof->inverseBeta(): MuonTOF= globalMinTOF_ ;
  dedxSObj ? Is = dedxSObj->dEdx() : Is = 0;
  dedxMObj ? Ih = dedxMObj->dEdx() : Ih = 0;
  //WAIT//float Ick=0; // if(dedxMObj) Ick=GetIck(Ih,isBckg);

  float PtCut = CutPt_[CutIndex];
  float ICut = CutI_[CutIndex];
  float TOFCut = CutTOF_[CutIndex];
  if (isFlip) {
    PtCut = CutPt_Flip_[CutIndex];
    ICut = CutI_Flip_[CutIndex];
    TOFCut = CutTOF_Flip_[CutIndex];
  }

  if (RescaleP) {
    if (RescaledPt(track->pt(), track->eta(), track->phi(), track->charge()) < PtCut)
      return false;
  } else {
    if (track->pt() < PtCut) {
      if (debug_ > 6) LogPrint(MOD) << "        >> @passSelection: p_T less than p_T cut (" << PtCut << ")";
      return false;
    }
  }
  
  if (tuple) {
    tuple->Pt->Fill(CutIndex, Event_Weight);
    if (GenBeta >= 0)
      tuple->Beta_SelectedP->Fill(CutIndex, GenBeta, Event_Weight);
  }

  if (typeMode_ != 3 && Is + RescaleI < ICut) {
    if (debug_ > 6) LogPrint(MOD) << "        >> @passSelection: I_s less than I_s cut (" << ICut << ")";
    return false;
  }

  if (tuple) {
    tuple->I->Fill(CutIndex, Event_Weight);
    if (GenBeta >= 0)
      tuple->Beta_SelectedI->Fill(CutIndex, GenBeta, Event_Weight);
  }

  if ((typeMode_ > 1 && typeMode_ != 5) && !isFlip && MuonTOF + RescaleT < TOFCut)
    return false;
  if ((typeMode_ > 1 && typeMode_ != 5) && isFlip && MuonTOF + RescaleT > TOFCut)
    return false;

  if (tuple) {
    tuple->TOF->Fill(CutIndex, Event_Weight);
    if (GenBeta >= 0)
      tuple->Beta_SelectedT->Fill(CutIndex, GenBeta, Event_Weight);
    tuple->AS_P->Fill(CutIndex, track->p(), Event_Weight);
    tuple->AS_Pt->Fill(CutIndex, track->pt(), Event_Weight);
    tuple->AS_Is->Fill(CutIndex, Is, Event_Weight);
    tuple->AS_Ih->Fill(CutIndex, Ih, Event_Weight);
    tuple->AS_TOF->Fill(CutIndex, MuonTOF, Event_Weight);
    //tuple->AS_EtaIs->Fill(CutIndex,track->eta(),Is,Event_Weight);
    //tuple->AS_EtaIh->Fill(CutIndex,track->eta(),Ih,Event_Weight);
    //tuple->AS_EtaP ->Fill(CutIndex,track->eta(),track->p(),Event_Weight);
    //tuple->AS_EtaPt->Fill(CutIndex,track->eta(),track->pt(),Event_Weight);
    tuple->AS_PIs->Fill(CutIndex, track->p(), Is, Event_Weight);
    tuple->AS_PIh->Fill(CutIndex, track->p(), Ih, Event_Weight);
    tuple->AS_PtIs->Fill(CutIndex, track->pt(), Is, Event_Weight);
    tuple->AS_PtIh->Fill(CutIndex, track->pt(), Ih, Event_Weight);
    tuple->AS_TOFIs->Fill(CutIndex, MuonTOF, Is, Event_Weight);
    tuple->AS_TOFIh->Fill(CutIndex, MuonTOF, Ih, Event_Weight);
  }
  return true;
}

//=============================================================
//
//     Check if track is from pixel
//
//=============================================================
void Analyzer::isPixelTrack(const edm::Ref<std::vector<Trajectory>>& refTraj, bool& isBpixtrack, bool& isFpixtrack) {
  // Used in analyze() to see if it is pixel track
  std::vector<TrajectoryMeasurement> tmeasColl = refTraj->measurements();
  std::vector<TrajectoryMeasurement>::const_iterator tmeasIt;
  for (tmeasIt = tmeasColl.begin(); tmeasIt != tmeasColl.end(); tmeasIt++) {
    if (!tmeasIt->updatedState().isValid())
      continue;
    TransientTrackingRecHit::ConstRecHitPointer testhit = tmeasIt->recHit();
    if (!testhit->isValid() || testhit->geographicalId().det() != DetId::Tracker)
      continue;
    uint testSubDetID = (testhit->geographicalId().subdetId());
    if (testSubDetID == PixelSubdetector::PixelBarrel)
      isBpixtrack = true;
    if (testSubDetID == PixelSubdetector::PixelEndcap)
      isFpixtrack = true;
    if (isBpixtrack && isFpixtrack)
      break;
  }
}

//=============================================================
//
//     Combine individual probs into a track level one
//
//=============================================================
float Analyzer::combineProbs(float probOnTrackWMulti, int numRecHits) const {
  float logprobOnTrackWMulti = probOnTrackWMulti > 0 ? log(probOnTrackWMulti) : 0;
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
                             const float pixelProbs[],
                             const float Event_Weight,
                             Tuple* tuple,
                             const float GenBeta,
                             float MassErr,
                             const bool Ih_Iso_cut,
                             const unsigned int closestBackgroundPDGsIDs[]) {
  //FIXME to be measured on 2015 data, currently assume 2012
  bool PRescale = true;
  float IRescale = -0.05;  // added to the Ias value
  float MRescale = 0.95;
  float TRescale = -0.015;  //-0.005 (used in 2012); // added to the 1/beta value

  // compute systematic due to momentum scale
  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, pixelProbs, EventWeight_, nullptr, -1, PRescale, 0, 0, 0, Ih_Iso_cut, closestBackgroundPDGsIDs)) {
    float RescalingFactor = RescaledPt(track->pt(), track->eta(), track->phi(), track->charge()) / track->pt();
    
    float Mass = -1;
    if (dedxMObj)
      Mass = GetMass(track->p() * RescalingFactor, dedxMObj->dEdx(), dEdxK_, dEdxC_);
    float MassTOF = -1;
    if (tof)
      MassTOF = GetTOFMass(track->p() * RescalingFactor, tof->inverseBeta());
    float MassComb = -1;
    if (tof && dedxMObj)
      MassComb = GetMassFromBeta(track->p() * RescalingFactor,
                                 (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
    else if (dedxMObj)
      MassComb = Mass;
    if (tof)
      MassComb = MassTOF;
    
    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
      if (passSelection(track,
                        dedxSObj,
                        dedxMObj,
                        tof,
                        iEvent,
                        EventWeight_,
                        CutIndex,
                        tuple,
                        false,
                        -1,
                        PRescale,
                        0,
                        0)) {  //WAIT//
        HSCPTk_SystP[CutIndex] = true;
        if (Mass > MaxMass_SystP[CutIndex])
          MaxMass_SystP[CutIndex] = Mass;
        tuple->Mass_SystP->Fill(CutIndex, Mass, EventWeight_);
        if (tof) {
          tuple->MassTOF_SystP->Fill(CutIndex, MassTOF, EventWeight_);
        }
        tuple->MassComb_SystP->Fill(CutIndex, MassComb, EventWeight_);
      }
    } // end loop on cut index
  } // end compute systematic due to momentum scale
  // compute systematic due to dEdx (both Ias and Ih)
  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, pixelProbs, EventWeight_, nullptr, -1, false, 0.0, IRescale, 0.0, false, closestBackgroundPDGsIDs)) {
      //if(TypeMode==5 && isSemiCosmicSB)continue;
    float Mass = -1;
    if (dedxMObj)
      Mass = GetMass(track->p(), dedxMObj->dEdx() * MRescale, dEdxK_, dEdxC_);
    float MassTOF = -1;
    if (tof)
      MassTOF = GetTOFMass(track->p(), tof->inverseBeta());
    float MassComb = -1;
    if (tof && dedxMObj)
      MassComb =
      GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
    else if (dedxMObj)
      MassComb = Mass;
    if (tof)
      MassComb = MassTOF;
    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
      if (passSelection(
                        track, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1, 0, IRescale, 0)) {
        HSCPTk_SystI[CutIndex] = true;
        if (Mass > MaxMass_SystI[CutIndex])
          MaxMass_SystI[CutIndex] = Mass;
        tuple->Mass_SystI->Fill(CutIndex, Mass, EventWeight_);
        if (tof)
          tuple->MassTOF_SystI->Fill(CutIndex, MassTOF, EventWeight_);
        tuple->MassComb_SystI->Fill(CutIndex, MassComb, EventWeight_);
      }
    }
  } // End compute systematic due to dEdx
  // compute systematic due to Mass shift ??????????
  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, pixelProbs, EventWeight_, nullptr, -1, 0, 0, 0, 0, Ih_Iso_cut, closestBackgroundPDGsIDs)) {
    /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
    float Mass = -1;
    if (dedxMObj)
      Mass = GetMass(track->p(), dedxMObj->dEdx() * MRescale, dEdxK_, dEdxC_);
    float MassTOF = -1;
    if (tof)
      MassTOF = GetTOFMass(track->p(), tof->inverseBeta());
    float MassComb = -1;
    if (tof && dedxMObj)
      MassComb = GetMassFromBeta(
                                 track->p(), (GetIBeta(dedxMObj->dEdx() * MRescale, dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
    else if (dedxMObj)
      MassComb = Mass;
    if (tof)
      MassComb = MassTOF;
    
    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
      if (passSelection(track, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1, 0, 0, 0)) {
        HSCPTk_SystM[CutIndex] = true;
        if (Mass > MaxMass_SystM[CutIndex])
          MaxMass_SystM[CutIndex] = Mass;
        tuple->Mass_SystM->Fill(CutIndex, Mass, EventWeight_);
        if (tof)
          tuple->MassTOF_SystM->Fill(CutIndex, MassTOF, EventWeight_);
        tuple->MassComb_SystM->Fill(CutIndex, MassComb, EventWeight_);
      }
    }
  } // End compute systematic due to Mass shift
  // compute systematic due to TOF
  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, pixelProbs, EventWeight_, nullptr, -1, 0, 0, TRescale, 0, Ih_Iso_cut, closestBackgroundPDGsIDs)) {
    /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
    float Mass = -1;
    if (dedxMObj)
      Mass = GetMass(track->p(), dedxMObj->dEdx(), dEdxK_, dEdxC_);
    float MassTOF = -1;
    if (tof)
      MassTOF = GetTOFMass(track->p(), (tof->inverseBeta() + TRescale));
    float MassComb = -1;
    if (tof && dedxMObj)
      MassComb = GetMassFromBeta(
                                 track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / (tof->inverseBeta() + TRescale))) * 0.5);
    else if (dedxMObj)
      MassComb = Mass;
    if (tof)
      MassComb = MassTOF;
    
    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
      if (passSelection(
                        track, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1, 0, 0, TRescale)) {
        HSCPTk_SystT[CutIndex] = true;
        if (Mass > MaxMass_SystT[CutIndex])
          MaxMass_SystT[CutIndex] = Mass;
        tuple->Mass_SystT->Fill(CutIndex, Mass, EventWeight_);
        if (tof)
          tuple->MassTOF_SystT->Fill(CutIndex, MassTOF, EventWeight_);
        tuple->MassComb_SystT->Fill(CutIndex, MassComb, EventWeight_);
      }
    }
  } // End condition for compute systematic due to TOF
  // compute systematics due to PU
  if (passPreselection(track, dedxHits, dedxSObj, dedxMObj, tof, iEvent, pixelProbs, EventWeight_, nullptr, -1, 0, 0, 0, 0, Ih_Iso_cut, closestBackgroundPDGsIDs)) {
    /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
    float Mass = -1;
    if (dedxMObj)
      Mass = GetMass(track->p(), dedxMObj->dEdx(), dEdxK_, dEdxC_);
    float MassTOF = -1;
    if (tof)
      MassTOF = GetTOFMass(track->p(), tof->inverseBeta());
    float MassComb = -1;
    if (tof && dedxMObj)
      MassComb =
      GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(), dEdxK_, dEdxC_) + (1 / tof->inverseBeta())) * 0.5);
    else if (dedxMObj)
      MassComb = Mass;
    if (tof)
      MassComb = MassTOF;
    
    for (unsigned int CutIndex = 0; CutIndex < CutPt_.size(); CutIndex++) {
      if (passSelection(track, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1, 0, 0, 0)) {
        HSCPTk_SystPU[CutIndex] = true;
        if (Mass > MaxMass_SystPU[CutIndex])
          MaxMass_SystPU[CutIndex] = Mass;
        tuple->Mass_SystPU->Fill(CutIndex, Mass, EventWeight_ * PUSystFactor_[0]);
        if (tof)
          tuple->MassTOF_SystPU->Fill(CutIndex, MassTOF, EventWeight_ * PUSystFactor_[0]);
        tuple->MassComb_SystPU->Fill(CutIndex, MassComb, EventWeight_ * PUSystFactor_[0]);
      }
    }
  }  // End compute systematics due to PU
}


