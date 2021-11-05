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
//


#include "SUSYBSMAnalysis/Analyzer/plugins/Analyzer.h"


Analyzer::Analyzer(const edm::ParameterSet& iConfig)
   // Read config file
   :hscpToken_(consumes<vector<susybsm::HSCParticle>>(iConfig.getParameter<edm::InputTag>("hscpCollection")))
   ,hscpIsoToken_(consumes<edm::ValueMap<susybsm::HSCPIsolation>>(iConfig.getParameter<edm::InputTag>("hscpIsoCollection")))
   ,muonSegmentToken_(consumes<susybsm::MuonSegmentCollection>(iConfig.getParameter<edm::InputTag>("muonSegmentCollection")))
   ,dedxToken_(consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dedxCollection")))
   ,muonTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonTimeCollection")))
   ,muonDtTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonDtTimeCollection")))
   ,muonCscTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonCscTimeCollection")))
   ,muonDtSegmentToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("muonDtSegmentCollection")))
   ,muonCscSegmentToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("muonCscSegmentCollection")))
   ,offlinePrimaryVerticesToken_(consumes<vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVerticesCollection")))
   ,lumiScalersToken_(consumes<LumiScalersCollection>(iConfig.getParameter<edm::InputTag>("lumiScalers")))
   ,refittedStandAloneMuonsToken_(consumes<vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("refittedStandAloneMuonsCollection")))
   ,offlineBeamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpotCollection")))
   ,muonToken_(consumes<vector<reco::Muon>>(iConfig.getParameter<edm::InputTag>("muonCollection")))
   ,triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults")))
   ,pfMETToken_(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("pfMET")))
   ,pfJetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJet")))
   ,CaloMETToken_(consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("CaloMET")))
   ,pileupInfoToken_(consumes<std::vector<PileupSummaryInfo>>(iConfig.getParameter<edm::InputTag>("pileupInfo")))
   ,genParticleToken_(consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticleCollection")))
   ,m_trajTag        ( consumes<TrajTrackAssociationCollection> (iConfig.getUntrackedParameter<edm::InputTag>("trajInputLabel")))
   // HLT triggers
   ,trigger_met_(iConfig.getUntrackedParameter<vector<string>>("Trigger_MET"))
   ,trigger_mu_(iConfig.getUntrackedParameter<vector<string>>("Trigger_Mu"))
   // =========Analysis parameters================
   ,TypeMode_(iConfig.getUntrackedParameter<unsigned int>("TypeMode"))
   ,SampleType_(iConfig.getUntrackedParameter<unsigned int>("SampleType"))
   ,SampleName_(iConfig.getUntrackedParameter<string>("SampleName"))
   ,Period_(iConfig.getUntrackedParameter<string>("Period"))
   ,SkipSelectionPlot_(iConfig.getUntrackedParameter<bool>("SkipSelectionPlot"))
   ,PtHistoUpperBound(iConfig.getUntrackedParameter<double>("PtHistoUpperBound"))
   ,MassHistoUpperBound(iConfig.getUntrackedParameter<double>("MassHistoUpperBound"))
   ,MassNBins(iConfig.getUntrackedParameter<unsigned int>("MassNBins"))
   ,IPbound(iConfig.getUntrackedParameter<double>("IPbound"))
   ,PredBins(iConfig.getUntrackedParameter<unsigned int>("PredBins"))
   ,EtaBins(iConfig.getUntrackedParameter<unsigned int>("EtaBins"))
   ,dEdxS_UpLim(iConfig.getUntrackedParameter<double>("dEdxS_UpLim"))
   ,dEdxM_UpLim(iConfig.getUntrackedParameter<double>("dEdxM_UpLim"))
   ,DzRegions(iConfig.getUntrackedParameter<unsigned int>("DzRegions"))
   ,GlobalMaxPterr(iConfig.getUntrackedParameter<double>("GlobalMaxPterr"))
   ,GlobalMinPt(iConfig.getUntrackedParameter<double>("GlobalMinPt"))
   ,GlobalMinTOF(iConfig.getUntrackedParameter<double>("GlobalMinTOF"))
   ,skipPixel(iConfig.getUntrackedParameter<bool>("skipPixel")) 
   ,useTemplateLayer(iConfig.getUntrackedParameter<bool>("useTemplateLayer"))
   //,DeDxSF_0(iConfig.getUntrackedParameter<double>("DeDxSF_0"))
   //,DeDxSF_1(iConfig.getUntrackedParameter<double>("DeDxSF_1"))
   //,DeDxK(iConfig.getUntrackedParameter<double>("DeDxK"))
   //,DeDxC(iConfig.getUntrackedParameter<double>("DeDxC"))
   ,DeDxTemplate(iConfig.getUntrackedParameter<string>("DeDxTemplate"))
   ,enableDeDxCalibration(iConfig.getUntrackedParameter<bool>("enableDeDxCalibration"))
   ,DeDxCalibration(iConfig.getUntrackedParameter<string>("DeDxCalibration"))
   ,Geometry(iConfig.getUntrackedParameter<string>("Geometry"))
   ,TimeOffset(iConfig.getUntrackedParameter<string>("TimeOffset"))
   ,FMIPX(iConfig.getUntrackedParameter<double>("FMIPX"))
   ,STree(iConfig.getUntrackedParameter<unsigned int>("saveTree"))
   ,SGTree(iConfig.getUntrackedParameter<unsigned int>("saveGenTree"))
{
   //now do what ever initialization is needed
   // define the selection to be considered later for the optimization
   // WARNING: recall that this has a huge impact on the analysis time AND on the output file size --> be carefull with your choice
  
   useClusterCleaning = true;
   if(TypeMode_==4) {
      useClusterCleaning = false; //switch off cluster cleaning for mCHAMPs
   }
   
   isData   = (SampleType_==0);
   isBckg   = (SampleType_==1);
   isSignal = (SampleType_>=2);

   //dEdxSF [0] = DeDxSF_0;
   //dEdxSF [1] = DeDxSF_1;

   bool splitByModuleType = true;
   dEdxTemplates = loadDeDxTemplate(DeDxTemplate, splitByModuleType);
   if(enableDeDxCalibration)   trackerCorrector.LoadDeDxCalibration(DeDxCalibration); 
   else                        trackerCorrector.TrackerGains = nullptr; //FIXME check gain for MC

   moduleGeom::loadGeometry(Geometry);
   tofCalculator.loadTimeOffset(TimeOffset);

}


Analyzer::~Analyzer()
{

   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}

// ------------ method called once each job just before starting event loop  ------------
void
Analyzer::beginJob()
{
   
   // Book histograms
   edm::Service<TFileService> fs;
   tuple = new Tuple();

   TFileDirectory dir = fs->mkdir( SampleName_.c_str(), SampleName_.c_str() );

   // create histograms & trees
   initializeCuts(fs, CutPt_, CutI_, CutTOF_, CutPt_Flip_, CutI_Flip_, CutTOF_Flip_);
   tuple_maker->initializeTuple(tuple, dir, STree, SGTree, SkipSelectionPlot_, TypeMode_, isSignal, CutPt_.size(), CutPt_Flip_.size(), PtHistoUpperBound, MassHistoUpperBound, MassNBins, IPbound, PredBins, EtaBins, dEdxS_UpLim, dEdxM_UpLim, DzRegions, GlobalMinPt, GlobalMinTOF);

   // Re-weighting
   mcWeight = new MCWeight();
   if(!isData){
      mcWeight->loadPileupWeights(Period_);
   }
   mcWeight->getSampleWeights(Period_, SampleName_.c_str(), IntegratedLuminosity_, CrossSection_);
   tuple->IntLumi->Fill(0.0,IntegratedLuminosity_);
   tuple->XSection->Fill(0.0,CrossSection_);

   tof    = nullptr;
   dttof  = nullptr;
   csctof = nullptr;

   TrigInfo_ = 0;

   CurrentRun_ = 0;
   RNG = new TRandom3();
   is2016 = false;
   is2016G = false;

   //PUSystFactor_.clear();
   PUSystFactor_.resize(2, 1.);
   PUSystFactor_[0]=PUSystFactor_[1]=0.;

   HSCPTk              = new bool[CutPt_.size()];
   HSCPTk_SystP        = new bool[CutPt_.size()];
   HSCPTk_SystI        = new bool[CutPt_.size()];
   HSCPTk_SystT        = new bool[CutPt_.size()];
   HSCPTk_SystM        = new bool[CutPt_.size()];
   HSCPTk_SystPU       = new bool[CutPt_.size()];
   HSCPTk_SystHUp      = new bool[CutPt_.size()];
   HSCPTk_SystHDown    = new bool[CutPt_.size()];
   MaxMass           = new double[CutPt_.size()];
   MaxMass_SystP     = new double[CutPt_.size()];
   MaxMass_SystI     = new double[CutPt_.size()];
   MaxMass_SystT     = new double[CutPt_.size()];
   MaxMass_SystM     = new double[CutPt_.size()];
   MaxMass_SystPU    = new double[CutPt_.size()];
   MaxMass_SystHUp   = new double[CutPt_.size()];
   MaxMass_SystHDown = new double[CutPt_.size()];

   /*HIPemulator.    setPeriodHIPRate(is2016G);
   HIPemulatorUp.  setPeriodHIPRate(is2016G, "ratePdfPixel_Up", "ratePdfStrip_Up");
   HIPemulatorDown.setPeriodHIPRate(is2016G, "ratePdfPixel_Up", "ratePdfStrip_Up");*/

   //HIPemulatorUp(false, "ratePdfPixel_Up", "ratePdfStrip_Up");
   //HIPemulatorDown(false, "ratePdfPixel_Down", "ratePdfStrip_Down");
}


//
// member functions
//

// ------------ method called for each event  ------------
void
Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
   using namespace edm;

   tuple->EventsTotal->Fill(0.0,EventWeight_);

   //if run change, update conditions
   if(CurrentRun_ != iEvent.id().run()){
      CurrentRun_  = iEvent.id().run();
      tofCalculator.setRun(CurrentRun_);
      trackerCorrector.setRun(CurrentRun_);

      loadDeDxParameters(CurrentRun_, SampleType_, DeDxSF_0, DeDxSF_1, DeDxK, DeDxC);
      dEdxSF [0] = DeDxSF_0;
      dEdxSF [1] = DeDxSF_1;

      //LogInfo("Analyzer") <<"------> dEdx parameters SF for Run "<<CurrentRun_<< ": "<< dEdxSF[1];
   }

   //WAIT////compute event weight
   //vector<PileupSummaryInfo> pileupInfo;
   if(!isData){
      /*double PUWeight = 1.;
      Handle<vector<PileupSummaryInfo> >  pileupInfoH; 
      iEvent.getByToken(pileupInfoToken_, pileupInfoH);
      if(pileupInfoH.isValid()){
         PUWeight = mcWeight->getEventPUWeight(pileupInfoH, PUSystFactor_);
      }
      else {LogWarning("Analyzer") << "PileupSummaryInfo Collection NotFound";}*/
      double PUWeight = mcWeight->getEventPUWeight(iEvent, pileupInfoToken_,PUSystFactor_);
      EventWeight_ = PUWeight; // 1. : unweighted w.r.t pileup
   }
   else
      EventWeight_ = 1.;
   
   vector<reco::GenParticle> genColl;
   double HSCPGenBeta1=-1, HSCPGenBeta2=-1;
   double HSCPDLength1=-1, HSCPDLength2=-1;
   if(isSignal){
      //get the collection of generated Particles
      Handle< vector<reco::GenParticle> >  genCollH;
      iEvent.getByToken(genParticleToken_, genCollH);
      if(!genCollH.isValid()) {LogWarning("Analyzer") << "Invalid GenParticle!!, this event will be ignored"; return;}

      genColl = *genCollH;

      int NChargedHSCP=HowManyChargedHSCP(genColl);
      float SignalEventWeight = mcWeight->getFGluinoWeight(NChargedHSCP, TypeMode_);

      GetGenHSCPDecayLength(genColl,HSCPDLength1,HSCPDLength2,true);
      tuple->Gen_DecayLength->Fill(HSCPDLength1, SignalEventWeight); //????
      tuple->Gen_DecayLength->Fill(HSCPDLength2, SignalEventWeight);

      GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,false);
      if(HSCPGenBeta1>=0)tuple->Beta_Gen->Fill(HSCPGenBeta1, SignalEventWeight);  
      if(HSCPGenBeta2>=0)tuple->Beta_Gen->Fill(HSCPGenBeta2, SignalEventWeight);

      GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,true);
      if(HSCPGenBeta1>=0)tuple->Beta_GenCharged->Fill(HSCPGenBeta1, SignalEventWeight); 
      if(HSCPGenBeta2>=0)tuple->Beta_GenCharged->Fill(HSCPGenBeta2, SignalEventWeight);

      // R-hadron wights needed due to wrong GenId---------------------------------BEGIN
      double  Wa=1.0, Wad=1.0, Waa=1.0, Wan=1.0; // Wa is additional weight for single other, Wad for other+double_charged,
                                                 // Waa for the event with 2 other R-hadron, Wan for other+neutral
      bool Rhadron=0; // default value - not R-hadron (not need to weight)
      string sample_name = "";
      mcWeight->getRHadronWeights(sample_name, Rhadron, Wa, Wad, Waa, Wan);
      // R-hadron wights needed due to wrong GenId---------------------------------END

      unsigned int nw=0, na=0, nd=0, nn=0; //initialize counters: nw - wrong, na - other, nd - double charged, nn - neutral

      for(auto const &gen : genColl){
         if(gen.pt()<5)continue;
         if(gen.status()!=1)continue;
         int AbsPdg=abs(gen.pdgId());
         if(AbsPdg<1000000 && AbsPdg!=17)continue;

         // categorise event with R-hadrons for additional weighting-----------------------BEGIN
         int GenId=gen.pdgId();
         if        (    GenId ==1000612||    GenId ==1092214)  { nw+=1;  // count wrong
         } else if (abs(GenId)==1006223||abs(GenId)==1092224)  { nd+=1;  // count doble charged
         } else if (abs(GenId)==1006113.||
                    abs(GenId)==1006333.||
                    abs(GenId)==1006313.||
                    abs(GenId)==1000622.||
                    abs(GenId)==1092114.||
                    abs(GenId)==1093324.||
                    abs(GenId)==1093214.||
                    abs(GenId)==1009333.||
                    abs(GenId)==1009223.||
                    abs(GenId)==1009113.||
                    abs(GenId)==1009313.||
                    abs(GenId)==1000993.){ nn+=1;  // count neutral
         } else if (AbsPdg>1000000)                            { na+=1;} // count other R-hadrons
         // categorise event with R-hadrons for additional weighting-----------------------BEGIN

         tuple->genlevelpT->Fill(gen.pt(), SignalEventWeight);
         tuple->genleveleta->Fill(gen.eta(), SignalEventWeight);
         tuple->genlevelbeta->Fill(gen.p()/gen.energy(), SignalEventWeight);
      }

   } //End of isSignal

   // new genHSCP ntuple after correcting weights
   int nrha=0;
   vector<float> genid; 
   vector<float> gencharge; 
   vector<float> genmass;
   vector<float> genpt; 
   vector<float> geneta; 
   vector<float> genphi;
   for(auto const &gen : genColl) {
      if(gen.pt()<5)continue;
      if(gen.status()!=1)continue;
      int AbsPdg=abs(gen.pdgId());
      if(AbsPdg<1000000 && AbsPdg!=17)continue;
   
      nrha++;
      //mk rhadron ntuple
      if(isSignal){
        genid.push_back(gen.pdgId());
        gencharge.push_back(gen.charge());
        genmass.push_back(gen.mass());
        genpt.push_back(gen.pt());
        geneta.push_back(gen.eta());
        genphi.push_back(gen.phi());
      }

   }
	nrha=0;


	if (isSignal)
     tuple_maker->fillGenTreeBranches(tuple, iEvent.id().run(),iEvent.id().event(),iEvent.id().luminosityBlock(), EventWeight_,genid,gencharge,genmass,genpt,geneta,genphi);

   edm::Handle<edm::TriggerResults> triggerH;
   iEvent.getByToken(triggerResultsToken_,triggerH); 
   const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerH);

   bool HLT_Mu50=false;
   bool HLT_PFMET120_PFMHT120_IDTight=false;
   bool HLT_PFHT500_PFMET100_PFMHT100_IDTight=false;
   bool HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60=false;
   bool HLT_MET105_IsoTrk50=false;

   for(unsigned int i=0; i<triggerH->size(); i++)
   {
        if(TString(triggerNames.triggerName(i)).Contains("HLT_Mu50") && triggerH->accept(i)) HLT_Mu50=true;
        if(TString(triggerNames.triggerName(i)).Contains("HLT_PFMET120_PFMHT120_IDTight") && triggerH->accept(i)) HLT_PFMET120_PFMHT120_IDTight=true;
        if(TString(triggerNames.triggerName(i)).Contains("HLT_PFHT500_PFMET100_PFMHT100_IDTight") && triggerH->accept(i)) HLT_PFHT500_PFMET100_PFMHT100_IDTight=true;
        if(TString(triggerNames.triggerName(i)).Contains("HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60") && triggerH->accept(i)) HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60=true;
        if(TString(triggerNames.triggerName(i)).Contains("HLT_MET105_IsoTrk50") && triggerH->accept(i)) HLT_MET105_IsoTrk50=true;
   }


   //check if the event is passing trigger
   tuple->TotalE  ->Fill(0.0,EventWeight_);
   tuple->TotalEPU->Fill(0.0,EventWeight_*PUSystFactor_[0]);
   //See if event passed signal triggers
   //WAIT//if(!PassTrigger(iEvent, isData, false, (is2016&&!is2016G)?&L1Emul:nullptr) ) {
   if(!passTrigger(iEvent, isData)){ return;
      //For TOF only analysis if the event doesn't pass the signal triggers check if it was triggered by the no BPTX cosmic trigger
      //If not TOF only then move to next event
      /*if(TypeMode_!=3) continue;
      if(!passTrigger(iEvent, isData, true, (is2016&&!is2016G)?&L1Emul:NULL)) continue;*/

      //If is cosmic event then switch plots to use to the ones for cosmics
	   //WAIT//SamplePlots=&plotsMap[CosmicName];
	}
   //WAIT//else if(TypeMode==3) {
	   //WAIT//SamplePlots = &plotsMap[samples[s].Name];
   //WAIT//}

   tuple->TotalTE->Fill(0.0,EventWeight_);

   //keep beta distribution for signal
   if(isSignal){
      if(HSCPGenBeta1>=0) tuple->Beta_Triggered->Fill(HSCPGenBeta1, EventWeight_); 
      if(HSCPGenBeta2>=0) tuple->Beta_Triggered->Fill(HSCPGenBeta2, EventWeight_);
   }
   

   //===================== Handle For DeDx Hits ==============
   Handle<reco::DeDxHitInfoAss> dedxCollH;
   iEvent.getByToken(dedxToken_,dedxCollH);

   //================= Handle For Muon TOF Combined ===============
   Handle<reco::MuonTimeExtraMap>     tofMap;
   iEvent.getByToken(muonTimeToken_,  tofMap);

   //================= Handle For Muon TOF DT ===============
   Handle<reco::MuonTimeExtraMap>       tofDtMap;
   iEvent.getByToken(muonDtTimeToken_,  tofDtMap);

   //================= Handle For Muon TOF CSC ===============
   Handle<reco::MuonTimeExtraMap>        tofCscMap;
   iEvent.getByToken(muonCscTimeToken_,  tofCscMap);

   //================= Handle For Muon DT/CSC Segment ===============
   Handle<CSCSegmentCollection> CSCSegmentCollH;
   Handle<DTRecSegment4DCollection> DTSegmentCollH;
   if(!isBckg){ //do not recompute TOF on MC background
      iEvent.getByToken(muonCscSegmentToken_, CSCSegmentCollH);
      if(!CSCSegmentCollH.isValid()){LogError("Analyzer") << "CSC Segment Collection not found!"; return;}

      iEvent.getByToken(muonDtSegmentToken_, DTSegmentCollH);
      if(!DTSegmentCollH.isValid()){LogError("Analyzer") << "DT Segment Collection not found!"; return;}
   }

   //reinitialize the bookeeping array for each event
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk        [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk_SystP  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk_SystI  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk_SystT  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk_SystM  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk_SystPU [CutIndex] = false; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk_SystHUp[CutIndex] = false; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  HSCPTk_SystHDown[CutIndex] = false; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass       [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass_SystP [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass_SystI [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass_SystT [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass_SystM [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass_SystPU[CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass_SystHUp [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){  MaxMass_SystHDown[CutIndex] = -1; }

   //WAIT//HIPemulator.setEventRate(); //take it from a pdf
   //WAIT//HIPemulatorUp.setEventRate(HIPemulator.getEventRatePixel()*1.25, HIPemulator.getEventRateStrip()*1.80);  // deltaPixel = 3.653981e+02, basePixel = 1.332625e+03; deltaStrip = 4.662832e+02, baseStrip = 5.958308e+02, from Run257805
   //WAIT//HIPemulatorDown.setEventRate(HIPemulator.getEventRatePixel()*0.75, HIPemulator.getEventRateStrip()*0.20); 

	//WAIT//HIPTrackLossEmul.SetHIPTrackLossRate(iEvent);

    //===================== Handle For PileUp ================
   unsigned int pileup_fromLumi = 0;
   edm::Handle<LumiScalersCollection> lumiScalers;
   iEvent.getByToken(lumiScalersToken_, lumiScalers);
   if(lumiScalers.isValid() && !lumiScalers->empty())
   {
       LumiScalersCollection::const_iterator scalit = lumiScalers->begin();
       pileup_fromLumi = scalit->pileup(); 
   }
    //===================== Handle For vertices ================ 
   vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken_);

   
   float CaloMET=-1, RecoPFMET=-1, RecoPFMHT=-1, HLTPFMET=-1, HLTPFMHT=-1;

    //===================== Handle For PFMET ===================
   edm::Handle<std::vector<reco::PFMET>> pfMETHandle;
   iEvent.getByToken(pfMETToken_, pfMETHandle);
   if(pfMETHandle.isValid() && !pfMETHandle->empty())
   {
   for(unsigned int i=0; i<pfMETHandle->size(); i++)
   {
       const reco::PFMET* pfMet = &(*pfMETHandle)[i];
       RecoPFMET = pfMet->et();
   }
   }
   
    //===================== Handle For CaloMET ===================
   edm::Handle<std::vector<reco::CaloMET>> CaloMETHandle;
   iEvent.getByToken(CaloMETToken_, CaloMETHandle);
   if(CaloMETHandle.isValid() && !CaloMETHandle->empty())
   {
   for(unsigned int i=0; i<CaloMETHandle->size(); i++)
   {
       const reco::CaloMET* calomet = &(*CaloMETHandle)[i];
       CaloMET = calomet->et();
   }
   }
    
    //===================== Handle For PFJet ===================
   edm::Handle<reco::PFJetCollection> pfJetHandle;
   iEvent.getByToken(pfJetToken_, pfJetHandle);
   if(pfJetHandle.isValid() && !pfJetHandle->empty())
   {
   const reco::PFJetCollection* pfJetColl = pfJetHandle.product();
   TLorentzVector pMHT;
   for(unsigned int i=0; i<pfJetColl->size(); i++)
   {
       const reco::PFJet* jet = &(*pfJetColl)[i];
       if(jet->pt() < 20 || abs(jet->eta()) > 5 || jet->chargedEmEnergyFraction()+jet->neutralEmEnergyFraction() > 0.9) continue;
       TLorentzVector p4(jet->pt()*cos(jet->phi()),jet->pt()*sin(jet->phi()),0,jet->et());
       pMHT += p4;
   }
   RecoPFMHT = pMHT.Pt();
   }
    
    

   //load all event collection that will be used later on (HSCP, dEdx and TOF)
   unsigned int HSCP_count = 0;

   std::vector<bool>         HSCP_passCutPt55;
   std::vector<bool>         HSCP_passPreselection_noIsolation_noIh;
   std::vector<bool>         HSCP_passPreselection;
   std::vector<bool>         HSCP_passSelection;
   std::vector<float>        HSCP_Charge;
   std::vector<float>        HSCP_Pt;
   std::vector<float>        HSCP_PtErr;
   std::vector<float>        HSCP_Ias;
   std::vector<float>        HSCP_Ih;
   std::vector<float>        HSCP_Ick; //return (Ih-C)/K
   std::vector<float>        HSCP_Fmip;
   std::vector<float>        HSCP_ProbQ;
   std::vector<float>        HSCP_TOF;
   std::vector<float>        HSCP_TOFErr;
   std::vector<unsigned int> HSCP_TOF_ndof;
   std::vector<float>        HSCP_DTTOF;
   std::vector<float>        HSCP_DTTOFErr;
   std::vector<unsigned int> HSCP_DTTOF_ndof;
   std::vector<float>        HSCP_CSCTOF;
   std::vector<float>        HSCP_CSCTOFErr;
   std::vector<unsigned int> HSCP_CSCTOF_ndof;
   std::vector<float>        HSCP_Mass;
   std::vector<float>        HSCP_MassErr;
   std::vector<float>        HSCP_dZ;
   std::vector<float>        HSCP_dXY;
   std::vector<float>        HSCP_dR;
   std::vector<float>        HSCP_eta;
   std::vector<float>        HSCP_phi;
   std::vector<unsigned int> HSCP_NOH; //number of (valid) track pixel+strip hits 
   std::vector<unsigned int> HSCP_NOPH;//number of (valid) track pixel hits
   std::vector<float>        HSCP_FOVH;//fraction of valid track hits
   std::vector<unsigned int> HSCP_NOMH;//number of missing hits from IP till last hit (excluding hits behind the last hit)
   std::vector<float>        HSCP_FOVHD;//fraction of valid hits divided by total expected hits until the last one
   std::vector<unsigned int> HSCP_NOM;//number of dEdx hits (= #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used)
   std::vector<float>        HSCP_iso_TK;
   std::vector<float>        HSCP_iso_ECAL;
   std::vector<float>        HSCP_iso_HCAL;   
   std::vector<float>        HSCP_PFIsolationR03_sumChargedHadronPt;
   std::vector<float>        HSCP_PFIsolationR03_sumNeutralHadronPt;
   std::vector<float>        HSCP_PFIsolationR03_sumPhotonPt;
   std::vector<float>        HSCP_PFIsolationR03_sumPUPt;
   std::vector<float>        HSCP_Ih_noL1;
   std::vector<float>        HSCP_Ih_15drop;
   std::vector<float>        HSCP_Ih_StripOnly;
   std::vector<float>        HSCP_Ih_StripOnly_15drop;
   std::vector<float>        HSCP_Ih_SaturationCorrectionFromFits;
   std::vector<std::vector<float>>      HSCP_clust_charge; //dedx charge -> either strip or pixel 
   std::vector<std::vector<float>>      HSCP_clust_pathlength;
   std::vector<std::vector<bool>>       HSCP_clust_ClusterCleaning;
   std::vector<std::vector<unsigned int>>        HSCP_clust_nstrip;
   std::vector<std::vector<bool>>       HSCP_clust_sat254;
   std::vector<std::vector<bool>>       HSCP_clust_sat255;
   std::vector<std::vector<uint32_t>>   HSCP_clust_detid;
   std::vector<std::vector<bool>>       HSCP_clust_isStrip; //is it a SiStrip cluster?
   std::vector<std::vector<bool>>       HSCP_clust_isPixel; //is it a Pixel hit?
   std::vector<float>        HSCP_GenId;
   std::vector<float>        HSCP_GenCharge;
   std::vector<float>        HSCP_GenMass;
   std::vector<float>        HSCP_GenPt;
   std::vector<float>        HSCP_GenEta;
   std::vector<float>        HSCP_GenPhi;

   //====================loop over HSCP candidates===================
   for(const auto& hscp : iEvent.get(hscpToken_)){
      reco::MuonRef  muon  = hscp.muonRef();//const reco::MuonRef& muon = hscp.muonRef();

      //For TOF only analysis use updated stand alone muon track.
	   //Otherwise use inner tracker track
	   reco::TrackRef track;
      if(TypeMode_!=3) track = hscp.trackRef();
      else {
		   if(muon.isNull()) continue;
		   track = muon->standAloneMuon();
      }
      //skip events without track
	   if(track.isNull())continue;
      // FIXME jozze skip events with |Eta| > 0.9 (out of the barrel)
	   //if(track->eta()>0.9 || track->eta() < -0.9) continue;

      //require a track segment in the muon system
      if(TypeMode_>1 && TypeMode_!=5 && (muon.isNull() || !muon->isStandAloneMuon()))continue;

      //Apply a scale factor to muon only analysis to account for differences seen in data/MC preselection efficiency
      //For eta regions where Data > MC no correction to be conservative
      if(!isData && TypeMode_==3 && scaleFactor(track->eta())<RNG->Uniform(0, 1)) continue;
   

      HSCP_count++;
      std::vector<float> clust_charge;
      std::vector<float> clust_pathlength;
      std::vector<bool> clust_ClusterCleaning;
      std::vector<unsigned int> clust_nstrip;
      std::vector<bool> clust_sat254;
      std::vector<bool> clust_sat255;
      std::vector<uint32_t> clust_detid;
      std::vector<bool> clust_isStrip;
      std::vector<bool> clust_isPixel;


      //for signal only, make sure that the candidate is associated to a true HSCP
      int ClosestGen;
      if(isSignal && DistToHSCP(hscp, genColl, ClosestGen, TypeMode_)>0.03)continue;

      // we are losing some tracks due to HIP
	   //WAIT//if(!isData && is2016 && !HIPTrackLossEmul.TrackSurvivesHIPInefficiency()) continue;

      //load quantity associated to this track (TOF and dEdx)
      const reco::DeDxHitInfo* dedxHits = nullptr;
      if(TypeMode_!=3 && !track.isNull()) {
         reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
         if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
      }
      
      if(TypeMode_>1 && TypeMode_!=5 && !hscp.muonRef().isNull()){
         if(isBckg){
            tof    = &(*tofMap)[hscp.muonRef()];
            dttof  = &(*tofDtMap)[hscp.muonRef()];
            csctof = &(*tofCscMap)[hscp.muonRef()];
         }else{
            const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollH;
            const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollH;
            //tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, isData?1:0 ); //apply T0 correction on data but not on signal MC
            tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, 1 ); //apply T0 correction on data but not on signal MC
            tof    = &tofCalculator.combinedTOF; 
            dttof  = &tofCalculator.dtTOF;  
            csctof = &tofCalculator.cscTOF;

         }
      }

      if(!dedxHits) continue; // skip tracks without hits otherwise there will be a crash


      int nofClust_dEdxLowerThan=0;
      float factorChargeToE = 3.61*pow(10,-6)*247;
      for(unsigned int i=0; i<dedxHits->size(); i++)
      {          
          clust_charge.push_back(dedxHits->charge(i));
          clust_pathlength.push_back(dedxHits->pathlength(i));
          clust_isStrip.push_back(dedxHits->detId(i)>=3?true:false);
          clust_isPixel.push_back(dedxHits->detId(i)>=3?false:true);
          clust_detid.push_back(dedxHits->detId(i));
          DetId detid(dedxHits->detId(i));
          if(detid.subdetId()<3) continue;
          const SiStripCluster* cluster = dedxHits->stripCluster(i);
          std::vector<int> ampl = convert(cluster->amplitudes());
          bool sat254=false, sat255=false;
          for(unsigned int s=0; s<ampl.size(); s++)
          {
              if(ampl[s]>=254) sat254=true;
              if(ampl[s]==255) sat255=true;
          }
          ampl = CrossTalkInv(ampl,0.10,0.04,true);
          clust_ClusterCleaning.push_back(clusterCleaning(dedxHits->stripCluster(i),1));
          clust_nstrip.push_back(ampl.size());
          clust_sat254.push_back(sat254);
          clust_sat255.push_back(sat255);
          if(dedxHits->charge(i)*factorChargeToE/dedxHits->pathlength(i)<FMIPX) nofClust_dEdxLowerThan++; 
      }
      float Fmip = (float)nofClust_dEdxLowerThan/(float)dedxHits->size(); 


      HitDeDxCollection hitDeDx = getHitDeDx(dedxHits, dEdxSF, trackerCorrector.TrackerGains, false, 1);

      unsigned int pdgId = 0;
      if(isSignal){ 
         pdgId = genColl[ClosestGen].pdgId();
         LogDebug("Analyzer") << "GenId  " << pdgId;
      }

      double dEdxErr = 0;
      reco::DeDxData dedxSObjTmp  = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, TypeMode_==5, false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.00, nullptr,0,pdgId,skipPixel,useTemplateLayer);
      reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr, &dEdxErr,pdgId,skipPixel,useTemplateLayer);
      reco::DeDxData dedxMUpObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr,0,pdgId,skipPixel,useTemplateLayer);
      reco::DeDxData dedxMDownObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr,0,pdgId,skipPixel,useTemplateLayer);
      /*reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulator:nullptr, &dEdxErr,pdgId);
      reco::DeDxData dedxMUpObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulatorUp:nullptr,0,pdgId);
      reco::DeDxData dedxMDownObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulatorDown:nullptr,0,pdgId);*/
      reco::DeDxData* dedxSObj  = dedxSObjTmp.numberOfMeasurements()>0?&dedxSObjTmp:nullptr;
      reco::DeDxData* dedxMObj  = dedxMObjTmp.numberOfMeasurements()>0?&dedxMObjTmp:nullptr;
      reco::DeDxData* dedxMUpObj = dedxMUpObjTmp.numberOfMeasurements()>0?&dedxMUpObjTmp:nullptr;
      reco::DeDxData* dedxMDownObj = dedxMDownObjTmp.numberOfMeasurements()>0?&dedxMDownObjTmp:nullptr;

    
      reco::DeDxData dedxIh_noL1_Tmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.0, nullptr, &dEdxErr,pdgId,skipPixel,useTemplateLayer,true);
      reco::DeDxData dedxIh_15drop_Tmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr, &dEdxErr,pdgId,skipPixel,useTemplateLayer);
      reco::DeDxData dedxIh_StripOnly_Tmp = computedEdx(dedxHits, dEdxSF, nullptr,          false, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.0, nullptr, &dEdxErr,pdgId,skipPixel,useTemplateLayer);
      reco::DeDxData dedxIh_StripOnly_15drop_Tmp = computedEdx(dedxHits, dEdxSF, nullptr,          false, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr, &dEdxErr,pdgId,skipPixel,useTemplateLayer,true);
      reco::DeDxData dedxIh_SaturationCorrectionFromFits_Tmp = computedEdx(dedxHits, dEdxSF, nullptr, false, useClusterCleaning, false, false, trackerCorrector.TrackerGains, true, true, 99, false, 2, 0.0, nullptr, &dEdxErr, pdgId, skipPixel, useTemplateLayer, true);


      reco::DeDxData* dedxIh_noL1               = dedxIh_noL1_Tmp.numberOfMeasurements()>0?&dedxIh_noL1_Tmp:nullptr;
      reco::DeDxData* dedxIh_15drop             = dedxIh_15drop_Tmp.numberOfMeasurements()>0?&dedxIh_15drop_Tmp:nullptr;
      reco::DeDxData* dedxIh_StripOnly          = dedxIh_StripOnly_Tmp.numberOfMeasurements()>0?&dedxIh_StripOnly_Tmp:nullptr;
      reco::DeDxData* dedxIh_StripOnly_15drop   = dedxIh_StripOnly_15drop_Tmp.numberOfMeasurements()>0?&dedxIh_StripOnly_15drop_Tmp:nullptr;
      reco::DeDxData* dedxIh_SaturationCorrectionFromFits   = dedxIh_SaturationCorrectionFromFits_Tmp.numberOfMeasurements()>0?&dedxIh_SaturationCorrectionFromFits_Tmp:nullptr;

      reco::DeDxData dedx_probQ_Tmp  = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, true, false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.00, nullptr,0,pdgId,skipPixel,useTemplateLayer,true);

      reco::DeDxData* dedx_probQ                = dedx_probQ_Tmp.numberOfMeasurements()>0?&dedx_probQ_Tmp:nullptr;


      if(TypeMode_==5)OpenAngle = deltaROpositeTrack(iEvent.get(hscpToken_), hscp); //OpenAngle is a global variable... that's uggly C++, but that's the best I found so far

      //compute systematic uncertainties on signal
      if(isSignal){
         //FIXME to be measured on 2015 data, currently assume 2012
         bool   PRescale = true;
         double IRescale =-0.05; // added to the Ias value
         double MRescale = 0.95;
		   double TRescale =-0.015; //-0.005 (used in 2012); // added to the 1/beta value
		  
		   double genpT = -1.0;
		   for(auto const gen : genColl) {
            if(gen.pt()<5)continue;
            if(gen.status()!=1)continue;
            int AbsPdg=abs(gen.pdgId());
            if(AbsPdg!=17)continue;
            
            double separation = deltaR(track->eta(), track->phi(), gen.eta(), gen.phi());
            if (separation > 0.03) continue;
            genpT = gen.pt();
            break;
         }
         if (genpT>0) {  tuple->genrecopT->Fill(genpT, track->pt()); }

         // compute systematic due to momentum scale
         //WAIT//if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   PRescale, 0, 0)){..}
         if(passPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, nullptr, -1,   PRescale, 0, 0, 0)){//WAIT//
            double RescalingFactor = RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())/track->pt();

            if(TypeMode_==5 && isSemiCosmicSB)continue;
            double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p()*RescalingFactor,dedxMObj->dEdx(), DeDxK, DeDxC);
            double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p()*RescalingFactor,tof->inverseBeta());
            double MassComb = -1;
		      if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p()*RescalingFactor, (GetIBeta(dedxMObj->dEdx(), DeDxK, DeDxC) + (1/tof->inverseBeta()))*0.5);
		      else if(dedxMObj) MassComb = Mass;
		      if(tof) MassComb=MassTOF;

            for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
		         if(passSelection(hscp,  dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1,   PRescale, 0, 0)){//WAIT//
                  HSCPTk_SystP[CutIndex] = true;
                  if(Mass>MaxMass_SystP[CutIndex]) MaxMass_SystP[CutIndex]=Mass;
                  tuple->Mass_SystP->Fill(CutIndex, Mass,EventWeight_);
                  if(tof)
                     {tuple->MassTOF_SystP ->Fill(CutIndex, MassTOF , EventWeight_);}
                     tuple->MassComb_SystP->Fill(CutIndex, MassComb, EventWeight_);
               }
            }
         }


         // compute systematic due to dEdx (both Ias and Ih)
         //WAIT//if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, IRescale, 0))
         if(passPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_,  nullptr, -1,   0, 0, IRescale, 0)){
            //if(TypeMode==5 && isSemiCosmicSB)continue;
            double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx()*MRescale, DeDxK, DeDxC);
            double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
            double MassComb = -1;
            if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(), DeDxK, DeDxC) + (1/tof->inverseBeta()))*0.5);
            else if(dedxMObj) MassComb = Mass;
            if(tof) MassComb=MassTOF;
            for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
               if(passSelection(hscp,  dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1,   0, IRescale, 0)){
                  HSCPTk_SystI[CutIndex] = true;
                  if(Mass>MaxMass_SystI[CutIndex]) MaxMass_SystI[CutIndex]=Mass;
                  tuple->Mass_SystI->Fill(CutIndex, Mass, EventWeight_);
                  if(tof) tuple->MassTOF_SystI ->Fill(CutIndex, MassTOF , EventWeight_);
                  tuple->MassComb_SystI->Fill(CutIndex, MassComb, EventWeight_);
               }
            }
         }

         // compute systematic due to Mass shift ??????????
         if(passPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, nullptr, -1,   0, 0, 0, 0)){
            /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
            double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx()*MRescale,DeDxK, DeDxC);
            double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
            double MassComb = -1;
            if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx()*MRescale,DeDxK, DeDxC) + (1/tof->inverseBeta()))*0.5);
            else if(dedxMObj) MassComb = Mass;
            if(tof) MassComb=MassTOF;

            for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
               if(passSelection(hscp,  dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1,   0, 0, 0)){
                  HSCPTk_SystM[CutIndex] = true;
                  if(Mass>MaxMass_SystM[CutIndex]) MaxMass_SystM[CutIndex]=Mass;
                  tuple->Mass_SystM->Fill(CutIndex, Mass,EventWeight_);
                  if(tof) tuple->MassTOF_SystM ->Fill(CutIndex, MassTOF , EventWeight_);
                  tuple->MassComb_SystM->Fill(CutIndex, MassComb, EventWeight_);
               }
            }
         }

         // compute systematic due to TOF
         if(passPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, nullptr, -1,   0, 0, TRescale, 0)){
            /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
            double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx(), DeDxK, DeDxC);
            double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),(tof->inverseBeta()+TRescale));
            double MassComb = -1;
            if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),DeDxK, DeDxC) + (1/(tof->inverseBeta()+TRescale)))*0.5);
            else if(dedxMObj) MassComb = Mass;
            if(tof) MassComb=MassTOF;

            for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
               if(passSelection(hscp,  dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1,   0, 0, TRescale)){
                  HSCPTk_SystT[CutIndex] = true;
                  if(Mass>MaxMass_SystT[CutIndex]) MaxMass_SystT[CutIndex]=Mass;
                  tuple->Mass_SystT->Fill(CutIndex, Mass,EventWeight_);
                  if(tof) tuple->MassTOF_SystT ->Fill(CutIndex, MassTOF , EventWeight_);
                  tuple->MassComb_SystT->Fill(CutIndex, MassComb, EventWeight_);
               }
            }
         }

         // compute systematics due to PU
         if(passPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, nullptr, -1,   0, 0, 0, 0)){
            /*if(TypeMode==5 && isSemiCosmicSB)continue;*/
            double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx(),DeDxK, DeDxC);
            double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
            double MassComb = -1;
            if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),DeDxK, DeDxC) + (1/tof->inverseBeta()))*0.5);
            else if(dedxMObj) MassComb = Mass;
            if(tof) MassComb=MassTOF;

            for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
               if(passSelection(hscp,  dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1,   0, 0, 0)){
                  HSCPTk_SystPU[CutIndex] = true;
                  if(Mass>MaxMass_SystPU[CutIndex]) MaxMass_SystPU[CutIndex]=Mass;
                  tuple->Mass_SystPU->Fill(CutIndex, Mass,EventWeight_*PUSystFactor_[0]);
                  if(tof) tuple->MassTOF_SystPU ->Fill(CutIndex, MassTOF , EventWeight_*PUSystFactor_[0]);
                  tuple->MassComb_SystPU->Fill(CutIndex, MassComb, EventWeight_*PUSystFactor_[0]);
               }
            }
         }

      }//End of systematic computation for signal

      //check if the canddiate pass the preselection cuts
      /*const susybsm::HSCParticle& hscp, const DeDxHitInfo* dedxHits,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const ChainEvent& ev, stPlots* st, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT, double MassErr*/
      double MassErr = GetMassErr(track->p(), track->ptError(), dedxMObj?dedxMObj->dEdx():-1, dEdxErr, GetMass(track->p(), dedxMObj?dedxMObj->dEdx():-1, DeDxK,DeDxC), DeDxK,DeDxC);
      if(isBckg)passPreselection( hscp, dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, tuple  , -1, false, 0, 0, MassErr);
      if(    !passPreselection( hscp, dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, tuple, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1, false, 0, 0, MassErr)) continue;
      /*if(TypeMode==5 && isSemiCosmicSB)continue;*/

      bool passPre = true;
      bool passPre_noIh_noIso = true;

      if(isBckg){
         passPreselection(  hscp, dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, tuple, -1, false, 0, 0, MassErr );
      }
      if(!passPreselection( hscp, dedxHits, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, tuple, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1, false, 0, 0, MassErr) ) passPre=false;
      
      if(!passPreselection( hscp, dedxHits, NULL, NULL, tof, iEvent, EventWeight_, NULL, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1, false, 0, 0, MassErr, false) ) passPre_noIh_noIso=false;

      if(TypeMode_==5 && isSemiCosmicSB)continue;//WAIT//

// Loop on track trajectory association map // Tav
    edm::Handle<TrajTrackAssociationCollection> hTTAC;
    iEvent.getByToken(m_trajTag, hTTAC);
    if (hTTAC.isValid())  {
      const TrajTrackAssociationCollection ttac = *(hTTAC.product());
        for (TrajTrackAssociationCollection::const_iterator it = ttac.begin(); it !=  ttac.end(); ++it){
            const edm::Ref<std::vector<Trajectory> > refTraj = it->key;
            const reco::TrackRef trackReference = it->val;

	    // track is defined in line546 as track = hscp.trackRef();
	    // so if this track doesnt belong to the HSCP in question let's skip it
            if (track != refTraj) continue;

            // -- Check whether it is a pixel track
            bool isBpixTrack(false), isFpixTrack(false);
            isPixelTrack(refTraj, isBpixTrack, isFpixTrack);
            if (!isBpixTrack && !isFpixTrack) { continue; }
            
            // -- Clusters associated with a track
            float probQonTrackWMulti = 1;
            float probXYonTrackWMulti = 1;
            std::vector<TrajectoryMeasurement> tmeasColl = refTraj->measurements();
            int numRecHits = 0;
            for (auto const& tmeasIt : tmeasColl) {
                if (!tmeasIt.updatedState().isValid()) continue;
                const TrackingRecHit* hit = tmeasIt.recHit()->hit();
                const SiPixelRecHit* pixhit = dynamic_cast<const SiPixelRecHit*>(hit);
                if (hit->geographicalId().det() != DetId::Tracker) continue;
                if (pixhit == nullptr) continue;
                if (!pixhit->isValid()) continue;
                float probQ         = pixhit->probabilityQ();
                float probXY        = pixhit->probabilityXY();
                numRecHits++;
                probQonTrackWMulti *= probQ; // \alpha_n in formula
                probXYonTrackWMulti *= probXY; // \alpha_n in formula
//              LocalPoint lp = pixhit->localPosition();
//              float rechit_x = lp.x();
//              float rechit_y = lp.y();
//              std::cout << "rechit_x: " << rechit_x << " and " << "rechit_y" << rechit_y << std::endl;
//              std::cout << "probQ: " << probQ << " and " << "probXY " << probXY << std::endl;
                

            } // end loop on on-track-clusters
            float logprobQonTrackWMulti = log(probQonTrackWMulti);
            float logprobXYonTrackWMulti = log(probXYonTrackWMulti);
//            std::cout << "numRecHits: " << numRecHits << std::endl;
            float probQonTrackTerm = 0;
            float probXYonTrackTerm = 0;
            for(int iTkCl = 0; iTkCl < numRecHits; ++iTkCl) {
                        probQonTrackTerm += ((pow(-logprobQonTrackWMulti,iTkCl))/(factorial(iTkCl)));
                        probXYonTrackTerm += ((pow(-logprobXYonTrackWMulti,iTkCl))/(factorial(iTkCl)));
//                        cout << "For cluster " << iTkCl << " the probQonTrackTerm is " << probQonTrackTerm << " the probXYonTrackTerm is " << probXYonTrackTerm <<  endl;
            }
            
            float probQonTrack = probQonTrackWMulti*probQonTrackTerm;
            float probXYonTrack = probXYonTrackWMulti*probXYonTrackTerm;
            std::cout << "For this track probQonTrack is " << probQonTrack << " and probXYonTrack is  " << probXYonTrack << endl;
            
        } // end loop TrajTrackAssociationCollection
    } else {
        std::cout << "hTTAC is invalid" << std::endl;
    }

//std::cout << "track is " << std::endl << std::endl;
//std::cout << track << std::endl << std::endl;

      //fill the ABCD histograms and a few other control plots
      //WAIT//if(isData)Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, SamplePlots);
      //WAIT//else if(isBckg) Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, MCTrPlots);

      if(passPre && !isSignal) tuple_maker->fillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, tuple, TypeMode_, GlobalMinTOF, EventWeight_,isCosmicSB, DTRegion, MaxPredBins, DeDxK, DeDxC, CutPt_, CutI_, CutTOF_, CutPt_Flip_, CutI_Flip_, CutTOF_Flip_);

      if(TypeMode_==5 && isCosmicSB)continue; 

      //Find the number of tracks passing selection for TOF<1 that will be used to check the background prediction
      //double Mass = -1;
      if(isBckg || isData) {
         //compute the mass of the candidate, for TOF mass flip the TOF over 1 to get the mass, so 0.8->1.2
		   double Mass = -1; if(dedxMObj) Mass = GetMass(track->p(),dedxMObj->dEdx(),DeDxK,DeDxC);
		   double MassTOF  = -1; if(tof) MassTOF = GetTOFMass(track->p(),(2-tof->inverseBeta()));
		   double MassComb = -1;
		   if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),DeDxK,DeDxC) + (1/(2-tof->inverseBeta())))*0.5 ) ;
		   if(dedxMObj) MassComb = Mass;
		   if(tof) MassComb=GetMassFromBeta(track->p(),(1/(2-tof->inverseBeta())));

           if(passPre)
           {
         for(unsigned int CutIndex=0;CutIndex<CutPt_Flip_.size();CutIndex++){
            //Background check looking at region with TOF<1
            //WAIT//if(!PassSelection   (hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, true)) continue;

            //Fill Mass Histograms
            tuple->Mass_Flip->Fill(CutIndex, Mass,EventWeight_);
            if(tof) tuple->MassTOF_Flip->Fill(CutIndex, MassTOF, EventWeight_);
            tuple->MassComb_Flip->Fill(CutIndex, MassComb, EventWeight_);
		   }
           }

      }

      //compute the mass of the candidate
      double Mass     = -1; if(dedxMObj) Mass = GetMass(track->p(),dedxMObj->dEdx(),DeDxK,DeDxC);
      double MassTOF  = -1; if(tof) MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
      double MassComb = -1;
      if(tof && dedxMObj) MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),DeDxK,DeDxC) + (1/tof->inverseBeta()))*0.5 ) ;
      if(dedxMObj) MassComb = Mass;
      if(tof) MassComb=GetMassFromBeta(track->p(),(1/tof->inverseBeta()));

      double MassUp    = -1; if(dedxMUpObj) MassUp=GetMass(track->p(),dedxMUpObj->dEdx(),DeDxK,DeDxC);
      double MassUpComb = -1;
      if(tof && dedxMUpObj) MassUpComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMUpObj->dEdx(),DeDxK,DeDxC) + (1/tof->inverseBeta()))*0.5 ) ;
      if(dedxMUpObj) MassUpComb = MassUp;
      if(tof) MassUpComb=GetMassFromBeta(track->p(),(1/tof->inverseBeta()));

      double MassDown    = -1; if(dedxMDownObj) MassDown=GetMass(track->p(),dedxMDownObj->dEdx(),DeDxK,DeDxC);
      double MassDownComb = -1;
      if(tof && dedxMDownObj) MassDownComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMDownObj->dEdx(),DeDxK,DeDxC) + (1/tof->inverseBeta()))*0.5 ) ;
      if(dedxMDownObj) MassDownComb = MassDown;
      if(tof) MassDownComb=GetMassFromBeta(track->p(),(1/tof->inverseBeta()));

      bool PassNonTrivialSelection=false;

      if(passPre)
      {

      //==========================================================
      // Cut loop: over all possible selection (one of them, the optimal one, will be used later)
      for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
         //Full Selection
         //if(isBckg)passSelection   (hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, MCTrPlots);
         if(isBckg) passSelection(hscp, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, -1, false, 0, 0);
         if(     !passSelection(hscp, dedxSObj, dedxMObj, tof, iEvent, EventWeight_, CutIndex, tuple, false, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1, false, 0, 0) ) continue;

         if(CutIndex!=0)PassNonTrivialSelection=true;
         HSCPTk[CutIndex] = true;
         HSCPTk_SystHUp[CutIndex] = true;
         HSCPTk_SystHDown[CutIndex] = true;

         if(Mass>MaxMass[CutIndex]) MaxMass[CutIndex]=Mass;
         if(MassUp>MaxMass_SystHUp[CutIndex]) MaxMass_SystHUp[CutIndex]=Mass;
         if(MassDown>MaxMass_SystHDown[CutIndex]) MaxMass_SystHDown[CutIndex]=Mass;

         //Fill Mass Histograms
         tuple->Mass->Fill(CutIndex, Mass,EventWeight_);
         if(tof) tuple->MassTOF->Fill(CutIndex, MassTOF, EventWeight_);
         if(isBckg) tuple->MassComb->Fill(CutIndex, MassComb, EventWeight_);

         //Fill Mass Histograms for different Ih syst
         tuple->Mass_SystHUp  ->Fill(CutIndex, MassUp,EventWeight_);
         tuple->Mass_SystHDown->Fill(CutIndex, MassDown,EventWeight_);
         if(tof) tuple->MassTOF_SystH ->Fill(CutIndex, MassTOF, EventWeight_);
         tuple->MassComb_SystHUp  ->Fill(CutIndex, MassUpComb, EventWeight_);
         tuple->MassComb_SystHDown->Fill(CutIndex, MassDownComb, EventWeight_);

      }//end of Cut loop
      }

      double Ick2=0;  if(dedxMObj) Ick2=GetIck(dedxMObj->dEdx(),DeDxK,DeDxC);
      int nomh= 0;nomh = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
      double fovhd = track->found()<=0?-1:track->found() / float(track->found() + nomh);
      unsigned int nom=0; if(dedxSObj) nom=dedxSObj->numberOfMeasurements();

      double genid=0,gencharge=-99,genmass=-99,genpt=-99,geneta=-99,genphi=-99;
  
      if(isSignal){
         genid = genColl[ClosestGen].pdgId();
         gencharge = genColl[ClosestGen].charge();
         genmass = genColl[ClosestGen].mass();
         genpt = genColl[ClosestGen].pt();
         geneta = genColl[ClosestGen].eta();
         genphi = genColl[ClosestGen].phi();
      }
    
      
      
      float iso_TK = -1;
      float iso_ECAL = -1;
      float iso_HCAL = -1;

      if(TypeMode_!=3) {
        const edm::ValueMap<susybsm::HSCPIsolation> IsolationMap = iEvent.get(hscpIsoToken_);
        susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
        iso_TK = hscpIso.Get_TK_SumEt();
        iso_ECAL = hscpIso.Get_ECAL_Energy();
        iso_HCAL = hscpIso.Get_HCAL_Energy();
      }

      float PFIso03_sumCharHadPt = -1;
      float PFIso03_sumNeutHadPt = -1;
      float PFIso03_sumPhotonPt = -1;
      float PFIso03_sumPUPt = -1;

      if(TypeMode_==2 && !muon.isNull())
      {
          PFIso03_sumCharHadPt = muon->pfIsolationR03().sumChargedHadronPt;
          PFIso03_sumNeutHadPt = muon->pfIsolationR03().sumNeutralHadronEt;
          PFIso03_sumPhotonPt = muon->pfIsolationR03().sumPhotonEt;
          PFIso03_sumPUPt = muon->pfIsolationR03().sumPUPt;
      }

      HSCP_passCutPt55.push_back(track->pt()>55?true:false);
      HSCP_passPreselection_noIsolation_noIh.push_back(passPre_noIh_noIso);
      HSCP_passPreselection.push_back(passPre);
      HSCP_passSelection.push_back(PassNonTrivialSelection);
      HSCP_Charge.push_back(track->charge());
      HSCP_Pt.push_back(track->pt());
      HSCP_PtErr.push_back(track->ptError());
      HSCP_Ias.push_back(dedxSObj?dedxSObj->dEdx():-1);
      HSCP_Ih.push_back(dedxMObj?dedxMObj->dEdx():-1);
      HSCP_Ick.push_back(dedxMObj?Ick2:-99);
      HSCP_Fmip.push_back(Fmip);
      HSCP_ProbQ.push_back(dedx_probQ?dedx_probQ->dEdx():-1);
      HSCP_TOF.push_back(tof?tof->inverseBeta():-99);
      HSCP_TOFErr.push_back(tof?tof->inverseBetaErr():-99);
      HSCP_TOF_ndof.push_back(tof?tof->nDof():-99);
      HSCP_DTTOF.push_back(dttof?dttof->inverseBeta():-99);
      HSCP_DTTOFErr.push_back(dttof?dttof->inverseBetaErr():-99);
      HSCP_DTTOF_ndof.push_back(dttof?dttof->nDof():-99);
      HSCP_CSCTOF.push_back(csctof?csctof->inverseBeta():-99);
      HSCP_CSCTOFErr.push_back(csctof?csctof->inverseBetaErr():-99);
      HSCP_CSCTOF_ndof.push_back(csctof?csctof->nDof():-99);
      HSCP_Mass.push_back(Mass);
      HSCP_MassErr.push_back(MassErr);
      HSCP_dZ.push_back(TreeDZ);
      HSCP_dXY.push_back(TreeDXY);
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
      HSCP_PFIsolationR03_sumChargedHadronPt.push_back(PFIso03_sumCharHadPt);
      HSCP_PFIsolationR03_sumNeutralHadronPt.push_back(PFIso03_sumNeutHadPt);
      HSCP_PFIsolationR03_sumPhotonPt.push_back(PFIso03_sumPhotonPt);
      HSCP_PFIsolationR03_sumPUPt.push_back(PFIso03_sumPUPt);
      HSCP_Ih_noL1.push_back(dedxIh_noL1 ? dedxIh_noL1->dEdx() : -1);
      HSCP_Ih_15drop.push_back(dedxIh_15drop ? dedxIh_15drop->dEdx() : -1);
      HSCP_Ih_StripOnly.push_back(dedxIh_StripOnly ? dedxIh_StripOnly->dEdx() : -1);
      HSCP_Ih_StripOnly_15drop.push_back(dedxIh_StripOnly_15drop ? dedxIh_StripOnly_15drop->dEdx() : -1);
      HSCP_Ih_SaturationCorrectionFromFits.push_back(dedxIh_SaturationCorrectionFromFits ? dedxIh_SaturationCorrectionFromFits->dEdx() : -1);
      HSCP_clust_charge.push_back(clust_charge);
      HSCP_clust_pathlength.push_back(clust_pathlength);
      HSCP_clust_ClusterCleaning.push_back(clust_ClusterCleaning);
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

      
   } //END loop over HSCP candidates

   tuple_maker->fillTreeBranches(tuple,
           TrigInfo_,
           iEvent.id().run(),
           iEvent.id().event(),
           iEvent.id().luminosityBlock(),
           pileup_fromLumi,
           vertexColl.size(),
           HSCP_count,
           EventWeight_,
           HLT_Mu50,
           HLT_PFMET120_PFMHT120_IDTight,
           HLT_PFHT500_PFMET100_PFMHT100_IDTight,
           HLT_PFMETNoMu120_PFMHTNoMu120_IDTight_PFHT60,
           HLT_MET105_IsoTrk50,
           CaloMET,
           RecoPFMET,
           RecoPFMHT,
           HLTPFMET,
           HLTPFMHT,
           HSCP_passCutPt55,
           HSCP_passPreselection_noIsolation_noIh,
           HSCP_passPreselection,
           HSCP_passSelection,
           HSCP_Charge,
           HSCP_Pt,
           HSCP_PtErr,
           HSCP_Ias,
           HSCP_Ih,
           HSCP_Ick,
           HSCP_Fmip,
           HSCP_ProbQ,
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
           HSCP_PFIsolationR03_sumChargedHadronPt,
           HSCP_PFIsolationR03_sumNeutralHadronPt,
           HSCP_PFIsolationR03_sumPhotonPt,
           HSCP_PFIsolationR03_sumPUPt,
           HSCP_Ih_noL1,
           HSCP_Ih_15drop,
           HSCP_Ih_StripOnly,
           HSCP_Ih_StripOnly_15drop,
           HSCP_Ih_SaturationCorrectionFromFits,
           HSCP_clust_charge,
           HSCP_clust_pathlength,
           HSCP_clust_ClusterCleaning,
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
           HSCP_GenPhi
    );


   //save event dependent information thanks to the bookkeeping
   for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
      if(HSCPTk[CutIndex]){
         tuple->HSCPE->Fill(CutIndex,EventWeight_);
         tuple->MaxEventMass->Fill(CutIndex,MaxMass[CutIndex], EventWeight_);
         if(isBckg){
            tuple->HSCPE->Fill(CutIndex,EventWeight_);
            tuple->MaxEventMass->Fill(CutIndex,MaxMass[CutIndex], EventWeight_);
         }
      }
      if(HSCPTk_SystP[CutIndex]){
         tuple->HSCPE_SystP       ->Fill(CutIndex,EventWeight_);
         tuple->MaxEventMass_SystP->Fill(CutIndex,MaxMass_SystP[CutIndex], EventWeight_);
      }
      if(HSCPTk_SystI[CutIndex]){
         tuple->HSCPE_SystI       ->Fill(CutIndex,EventWeight_);
         tuple->MaxEventMass_SystI->Fill(CutIndex,MaxMass_SystI[CutIndex], EventWeight_);
      }
      if(HSCPTk_SystM[CutIndex]){
         tuple->HSCPE_SystM       ->Fill(CutIndex,EventWeight_);
         tuple->MaxEventMass_SystM->Fill(CutIndex,MaxMass_SystM[CutIndex], EventWeight_);
      }
      if(HSCPTk_SystT[CutIndex]){
         tuple->HSCPE_SystT       ->Fill(CutIndex,EventWeight_);
         tuple->MaxEventMass_SystT->Fill(CutIndex,MaxMass_SystT[CutIndex], EventWeight_);
      }
      if(HSCPTk_SystPU[CutIndex]){
         tuple->HSCPE_SystPU       ->Fill(CutIndex,EventWeight_*PUSystFactor_[0]);
         tuple->MaxEventMass_SystPU->Fill(CutIndex,MaxMass_SystPU[CutIndex], EventWeight_*PUSystFactor_[0]);
      }
      if(HSCPTk_SystHUp[CutIndex]){
         tuple->HSCPE_SystHUp     ->Fill(CutIndex,EventWeight_);
         tuple->MaxEventMass_SystHUp   ->Fill(CutIndex,MaxMass_SystHUp   [CutIndex], EventWeight_);
      }
      if(HSCPTk_SystHDown[CutIndex]){
         tuple->HSCPE_SystHDown   ->Fill(CutIndex,EventWeight_);
         tuple->MaxEventMass_SystHDown ->Fill(CutIndex,MaxMass_SystHDown [CutIndex], EventWeight_);
      }
   }

#ifdef THIS_IS_AN_EVENTSETUP_EXAMPLE
   ESHandle<SetupData> pSetup;
   iSetup.get<SetupRecord>().get(pSetup);
#endif
}

// ------------ method called once each job just after ending the event loop  ------------
void
Analyzer::endJob()
{
   delete RNG;
   delete [] HSCPTk;
   delete [] HSCPTk_SystP;
   delete [] HSCPTk_SystI;
   delete [] HSCPTk_SystT;
   delete [] HSCPTk_SystM;
   delete [] HSCPTk_SystPU;
   delete [] HSCPTk_SystHUp;
   delete [] HSCPTk_SystHDown;
   delete [] MaxMass;
   delete [] MaxMass_SystP;
   delete [] MaxMass_SystI;
   delete [] MaxMass_SystT;
   delete [] MaxMass_SystM;
   delete [] MaxMass_SystPU;
   delete [] MaxMass_SystHUp;
   delete [] MaxMass_SystHDown;
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);

  //Specify that only 'tracks' is allowed
  //To use, remove the default given above and uncomment below
  //ParameterSetDescription desc;
  //desc.addUntracked<edm::InputTag>("tracks","ctfWithMaterialTracks");
  //descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(Analyzer);

//=============================================================
//
//     Method for initializing pT and Is cuts
//
//=============================================================
void Analyzer::initializeCuts(edm::Service<TFileService> &fs, vector<double>  &CutPt, vector<double>  &CutI, vector<double>  &CutTOF, vector<double>  &CutPt_Flip, vector<double>  &CutI_Flip, vector<double>  &CutTOF_Flip){
   CutPt.clear();       CutI.clear();       CutTOF.clear();      
   CutPt_Flip.clear();  CutI_Flip.clear();  CutTOF_Flip.clear();  
   
   CutPt     .push_back(GlobalMinPt);   CutI       .push_back(GlobalMinIs);  CutTOF     .push_back(GlobalMinTOF);
   CutPt_Flip.push_back(GlobalMinPt);   CutI_Flip  .push_back(GlobalMinIs);  CutTOF_Flip.push_back(GlobalMinTOF);

   if(TypeMode_<2){   
      for(double Pt =GlobalMinPt+5 ; Pt <200;Pt+=5){
         for(double I  =GlobalMinIs+0.025  ; I  <0.45 ;I+=0.025){
            CutPt .push_back(Pt);   CutI  .push_back(I);  CutTOF.push_back(-1);
         }
      }
   }else if(TypeMode_==2){
      for(double Pt =GlobalMinPt+5 ; Pt <120;  Pt+=5){
      if(Pt>80 && ((int)Pt)%10!=0)continue;
      for(double I  =GlobalMinIs +0.025; I  <0.40;  I+=0.025){
      for(double TOF=GlobalMinTOF+0.025; TOF<1.35;TOF+=0.025){
	   CutPt .push_back(Pt);   CutI  .push_back(I);  CutTOF.push_back(TOF);
      }}}
      for(double Pt =GlobalMinPt+10 ; Pt <90;  Pt+=30){
      for(double I  =GlobalMinIs +0.1; I  <0.30;  I+=0.1){
      for(double TOF=GlobalMinTOF-0.05; TOF>0.65;TOF-=0.05){
	   CutPt_Flip .push_back(Pt);   CutI_Flip  .push_back(I);  CutTOF_Flip.push_back(TOF);
      }}}
   }else if(TypeMode_==3){
      for(double Pt =GlobalMinPt+30 ; Pt <450;  Pt+=30){
      for(double TOF=GlobalMinTOF+0.025; TOF<1.5;TOF+=0.025){
         CutPt .push_back(Pt);   CutI  .push_back(-1);  CutTOF.push_back(TOF);
      }}
      for(double Pt =GlobalMinPt+30 ; Pt <450;  Pt+=60){
      for(double TOF=GlobalMinTOF-0.025; TOF>0.5;TOF-=0.025){
      CutPt_Flip .push_back(Pt);   CutI_Flip  .push_back(-1);  CutTOF_Flip.push_back(TOF);
      }}
   }else if(TypeMode_==4){
      for(double I  =GlobalMinIs +0.025; I  <0.55;  I+=0.025){
      for(double TOF=GlobalMinTOF+0.025; TOF<1.46;TOF+=0.025){
 	   CutPt .push_back(-1);   CutI  .push_back(I);  CutTOF.push_back(TOF);
       }}
      for(double I  =GlobalMinIs +0.025; I  <0.55;  I+=0.025){
      for(double TOF=GlobalMinTOF-0.025; TOF>0.54;TOF-=0.025){
	 CutPt_Flip .push_back(-1);   CutI_Flip  .push_back(I);  CutTOF_Flip.push_back(TOF);
       }}
   }else if(TypeMode_==5){   
      for(double Pt =75 ; Pt <=150;Pt+=25){
      for(double I  =0.0; I  <=0.45 ;I+=0.025){
         CutPt     .push_back(Pt);   CutI     .push_back(I);  CutTOF     .push_back(-1);
         CutPt_Flip.push_back(Pt);   CutI_Flip.push_back(I);  CutTOF_Flip.push_back(-1);
     }}
   }

   //printf("%i Different Final Selection will be tested\n",(int)CutPt.size());
   //printf("%i Different Final Selection will be tested for background uncertainty\n",(int)CutPt_Flip.size());
   edm::LogInfo("Analyzer") << CutPt.size() << " Different Final Selection will be tested\n"
                       << CutPt_Flip.size() << " Different Final Selection will be tested for background uncertainty";

   //Initialization of variables that are common to all samples
   HCuts["Pt"]  = fs->make<TProfile>("HCuts_Pt" ,"HCuts_Pt" ,CutPt.size(),0,CutPt.size());
   HCuts["I"]   = fs->make<TProfile>("HCuts_I"  ,"HCuts_I"  ,CutPt.size(),0,CutPt.size());
   HCuts["TOF"] = fs->make<TProfile>("HCuts_TOF","HCuts_TOF",CutPt.size(),0,CutPt.size());
   for(unsigned int i=0;i<CutPt.size();i++){  
      HCuts["Pt"]->Fill(i,CutPt[i]);     
      HCuts["I"]->Fill(i,CutI[i]);    
      HCuts["TOF"]->Fill(i,CutTOF[i]); 
   }

   HCuts["Pt_Flip"]  = fs->make<TProfile>("HCuts_Pt_Flip" ,"HCuts_Pt_Flip" ,CutPt_Flip.size(),0,CutPt_Flip.size());
   HCuts["I_Flip"]   = fs->make<TProfile>("HCuts_I_Flip"  ,"HCuts_I_Flip"  ,CutPt_Flip.size(),0,CutPt_Flip.size());
   HCuts["TOF_Flip"] = fs->make<TProfile>("HCuts_TOF_Flip","HCuts_TOF_Flip",CutPt_Flip.size(),0,CutPt_Flip.size());
   for(unsigned int i=0;i<CutPt_Flip.size();i++){  
      HCuts["Pt_Flip"]->Fill(i,CutPt_Flip[i]);     
      HCuts["I_Flip"]->Fill(i,CutI_Flip[i]);    
      HCuts["TOF_Flip"]->Fill(i,CutTOF_Flip[i]); 
   }
}

//=============================================================
//
//     Method for scaling eta per bin
//
//=============================================================
double Analyzer::scaleFactor(double eta) {
  double etaBins[15]   = {-2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0 , 0.3 , 0.6 , 0.9 , 1.2 ,1.5 , 1.8 , 2.1 };
  double scaleBins[15] = {0,    0.97, 1.06, 1.00, 0.89, 0.91, 0.93, 0.93, 0.92, 0.92, 0.91, 0.89,1.00, 1.06, 0.99};
  for (int i=0; i<15; i++) if(eta<etaBins[i]) return scaleBins[i];
  return 0;
}

//=============================================================
//
//     Method for rescaling pT
//
//=============================================================
double Analyzer::RescaledPt(const double& pt, const double& eta, const double& phi, const int& charge){
  if(TypeMode_!=3) {
    double newInvPt = 1/pt+0.000236-0.000135*pow(eta,2)+charge*0.000282*TMath::Sin(phi-1.337);
    return 1/newInvPt;
  }
  else {
    double newInvPt = (1./pt)*1.1;
    return 1/newInvPt;
  }
}

//=============================================================
//
//     Method to get hit position
//
//=============================================================
TVector3 Analyzer::getOuterHitPos(const reco::DeDxHitInfo* dedxHits){
     TVector3 point(0,0,0);
     if(!dedxHits)return point;
     //WAIT//double outerDistance=-1;
     for(unsigned int h=0;h<dedxHits->size();h++){
        DetId detid(dedxHits->detId(h));  
        //WAIT//moduleGeom* geomDet = moduleGeom::get(detid.rawId());
        //WAIT//TVector3 hitPos = geomDet->toGlobal(TVector3(dedxHits->pos(h).x(), dedxHits->pos(h).y(), dedxHits->pos(h).z())); 
        //WAIT//if(hitPos.Mag()>outerDistance){outerDistance=hitPos.Mag();  point=hitPos;}
     }
     return point;
}

//=============================================================
//
//     Method for ...
//
//=============================================================
double Analyzer::SegSep(const susybsm::HSCParticle& hscp, const edm::Event& iEvent, double& minPhi, double& minEta){
  if(TypeMode_!=3)return -1;

  reco::MuonRef muon = hscp.muonRef();
  if(muon.isNull()) return false;
  reco::TrackRef  track = muon->standAloneMuon();
  if(track.isNull())return false;


  double minDr=10;
  minPhi=10;
  minEta=10;

  //Look for segment on opposite side of detector from track
  //susybsm::MuonSegmentCollection SegCollection = iEvent.get(muonSegmentToken_);
  //for (susybsm::MuonSegmentCollection::const_iterator segment = SegCollection.begin(); segment!=SegCollection.end();++segment) {  
  for( const auto& segment : iEvent.get(muonSegmentToken_) ){
    GlobalPoint gp = segment.getGP();

    //Flip HSCP to opposite side of detector
    double eta_hscp = -1*track->eta();
    double phi_hscp= track->phi()+M_PI;

    double deta = gp.eta() - eta_hscp;
    double dphi = gp.phi() - phi_hscp;
    while (dphi >   M_PI) dphi -= 2*M_PI;
    while (dphi <= -M_PI) dphi += 2*M_PI;

    //Find segment most opposite in eta
    //Require phi difference of 0.5 so it doesn't match to own segment
    if(fabs(deta)<fabs(minEta) && fabs(dphi)<(M_PI-0.5)) {
      minEta=deta;
    }
    //Find segment most opposite in phi
    if(fabs(dphi)<fabs(minPhi)) {
      minPhi=dphi;
    }
    //Find segment most opposite in Eta-Phi
    double dR=sqrt(deta*deta+dphi*dphi);
    if(dR<minDr) minDr=dR;
  }
  return minDr;
}

//=============================================================
//
//     Pre-Selection
//
//=============================================================
bool Analyzer::passPreselection(
   const susybsm::HSCParticle& hscp, 
   const reco::DeDxHitInfo* dedxHits,  
   const reco::DeDxData* dedxSObj, 
   const reco::DeDxData* dedxMObj,
   const reco::MuonTimeExtra* tof,
   const edm::Event& iEvent, 
   float Event_Weight,
   Tuple* tuple, 
   const double& GenBeta, 
   bool RescaleP, 
   const double& RescaleI, 
   const double& RescaleT, 
   double MassErr,
   bool Ih_Iso_cut)
{

   if(TypeMode_==1 && !(hscp.type() == susybsm::HSCParticleType::trackerMuon || hscp.type() == susybsm::HSCParticleType::globalMuon))return false;
   if( (TypeMode_==2 || TypeMode_==4) && hscp.type() != susybsm::HSCParticleType::globalMuon)return false;
   reco::TrackRef   track;
   reco::MuonRef muon = hscp.muonRef();

   if(TypeMode_!=3) track = hscp.trackRef();
   else {
     if(muon.isNull()) return false;
     track = muon->standAloneMuon();
   }
   if(track.isNull())return false;

   if(tuple){tuple->Total->Fill(0.0,Event_Weight);
     if(GenBeta>=0)tuple->Beta_Matched->Fill(GenBeta, Event_Weight);
     tuple->BS_Eta->Fill(track->eta(),Event_Weight);
   }

   if(fabs(track->eta())>GlobalMaxEta) return false;

   //Cut on number of matched muon stations
   int count = muonStations(track->hitPattern());
   if(tuple) {
     tuple->BS_MatchedStations->Fill(count, Event_Weight);
   }
   if(TypeMode_==3 && count<minMuStations) return false;
   if(tuple) tuple->Stations->Fill(0.0, Event_Weight);

   

   //===================== Handle For vertex ================
   vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken_);
   if(vertexColl.size()<1){edm::LogError("Analyzer") << "NO VERTEX"; return false;}

   int highestPtGoodVertex = -1;
   int goodVerts=0;
   double dzMin=10000;
   for(unsigned int i=0;i<vertexColl.size();i++){
      if(vertexColl[i].isFake() || fabs(vertexColl[i].z())>24 || vertexColl[i].position().rho()>2 || vertexColl[i].ndof()<=4)continue; //only consider good vertex
      goodVerts++;
      if(tuple) tuple->BS_dzAll->Fill( track->dz (vertexColl[i].position()),Event_Weight);
      if(tuple) tuple->BS_dxyAll->Fill(track->dxy(vertexColl[i].position()),Event_Weight);

      if(fabs(track->dz (vertexColl[i].position())) < fabs(dzMin) ){
         dzMin = fabs(track->dz (vertexColl[i].position()));
         highestPtGoodVertex = i;
      }
   }
   if(highestPtGoodVertex<0)highestPtGoodVertex=0;

   if(tuple){tuple->BS_NVertex->Fill(vertexColl.size(), Event_Weight);
     tuple->BS_NVertex_NoEventWeight->Fill(vertexColl.size());
   }
   double dz  = track->dz (vertexColl[highestPtGoodVertex].position());
   double dxy = track->dxy(vertexColl[highestPtGoodVertex].position());

   bool PUA = (vertexColl.size()<15);
   bool PUB = (vertexColl.size()>=15);

   if(tuple){tuple->BS_TNOH->Fill(track->found(),Event_Weight);
          if(PUA)tuple->BS_TNOH_PUA->Fill(track->found(),Event_Weight);
          if(PUB)tuple->BS_TNOH_PUB->Fill(track->found(),Event_Weight);
          tuple->BS_TNOHFraction->Fill(track->validFraction(),Event_Weight);
	  tuple->BS_TNOPH->Fill(track->hitPattern().numberOfValidPixelHits(),Event_Weight);
   }

   if(TypeMode_!=3 && track->found()<GlobalMinNOH)return false;

   if(TypeMode_!=3 && track->hitPattern().numberOfValidPixelHits()<GlobalMinNOPH)return false;
   if(TypeMode_!=3 && track->validFraction()<GlobalMinFOVH)return false;

   unsigned int missingHitsTillLast = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);;
   double validFractionTillLast = track->found()<=0?-1:track->found() / float(track->found() + missingHitsTillLast);
  
   if(tuple){
      tuple->BS_TNOHFractionTillLast->Fill(validFractionTillLast,Event_Weight);
	   tuple->BS_TNOMHTillLast->Fill(missingHitsTillLast,Event_Weight);
   }

   if(TypeMode_!=3 && missingHitsTillLast>GlobalMaxNOMHTillLast)return false;
   if(TypeMode_!=3 && validFractionTillLast<GlobalMinFOVHTillLast)return false;

   if(tuple){
      tuple->TNOH  ->Fill(0.0,Event_Weight);
      if(dedxSObj){
         tuple->BS_TNOM->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);
         if(track->found() - dedxSObj->numberOfMeasurements()) 
             tuple->BS_EtaNBH->Fill(track->eta(), track->found() - dedxSObj->numberOfMeasurements(), Event_Weight);
         if(PUA)tuple->BS_TNOM_PUA->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);
         if(PUB)tuple->BS_TNOM_PUB->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);
      }
   }
   if(dedxSObj) if(dedxSObj->numberOfMeasurements()<GlobalMinNOM)return false;
   if(tuple){tuple->TNOM  ->Fill(0.0,Event_Weight);}

   if(tof){
      if(tuple){tuple->BS_nDof->Fill(tof->nDof(),Event_Weight);}
      if((TypeMode_>1  && TypeMode_!=5) && tof->nDof()<GlobalMinNDOF && (dttof->nDof()<GlobalMinNDOFDT || csctof->nDof()<GlobalMinNDOFCSC) )return false;
   }
   if(tuple){
      tuple->nDof  ->Fill(0.0,Event_Weight);
      tuple->BS_Qual->Fill(track->qualityMask(),Event_Weight);
   }

   if(TypeMode_!=3 && track->qualityMask()<GlobalMinQual )return false; // FIXME Tracks with quality > 2 are bad also!
//   if(TypeMode_!=3 && track->qualityMask() != FixedQual)return false; // FIXME if this is true, no tracks pass eventually ... so what now?
   if(tuple){tuple->Qual  ->Fill(0.0,Event_Weight);
          tuple->BS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight);
   }
   if(TypeMode_!=3 && track->chi2()/track->ndof()>GlobalMaxChi2 )return false;
   if(tuple){tuple->Chi2  ->Fill(0.0,Event_Weight);}

   if(tuple && GenBeta>=0)tuple->Beta_PreselectedA->Fill(GenBeta, Event_Weight);

   if(tuple){tuple->BS_MPt ->Fill(track->pt(),Event_Weight);}
   if(RescaleP){ if(RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())<GlobalMinPt)return false;
   }else{        if(track->pt()<GlobalMinPt)return false;   }

   if(tuple){tuple->MPt   ->Fill(0.0,Event_Weight);
     if(dedxSObj) tuple->BS_MIs->Fill(dedxSObj->dEdx(),Event_Weight);
     if(dedxMObj) tuple->BS_MIm->Fill(dedxMObj->dEdx(),Event_Weight);
   }

   if(dedxSObj) {if(dedxSObj->dEdx()+RescaleI<GlobalMinIs)return false;}
   if(dedxMObj){if((TypeMode_!=5 && dedxMObj->dEdx()<GlobalMinIm) || (TypeMode_==5 && dedxMObj->dEdx()>GlobalMinIm)) return false;}
   if(tuple){tuple->MI   ->Fill(0.0,Event_Weight);}

   if(tof){
      if(tuple){tuple->BS_MTOF ->Fill(tof->inverseBeta(),Event_Weight);}
      //This cut is no longer applied here but rather in the PassSelection part to use the region
      //with TOF<GlobalMinTOF as a background check
      //if(TypeMode_>1 && tof->inverseBeta()+RescaleT<GlobalMinTOF)return false;

      if(tuple)tuple->BS_TOFError->Fill(tof->inverseBetaErr(),Event_Weight);
      if((TypeMode_>1  && TypeMode_!=5) && tof->inverseBetaErr()>GlobalMaxTOFErr)return false;

      if(tuple) tuple->BS_TimeAtIP->Fill(tof->timeAtIpInOut(),Event_Weight);
      if(TypeMode_==3 && min(min(fabs(tof->timeAtIpInOut()-100), fabs(tof->timeAtIpInOut()-50)), min(fabs(tof->timeAtIpInOut()+100), fabs(tof->timeAtIpInOut()+50)))<5) return false;
   }

   if(tuple) tuple->BS_dzMinv3d->Fill(dz,Event_Weight);
   if(tuple) tuple->BS_dxyMinv3d->Fill(dxy,Event_Weight);
   if(tuple) tuple->BS_PV->Fill(goodVerts,Event_Weight);   
   if(tuple) tuple->BS_PV_NoEventWeight->Fill(goodVerts);
   if(tuple && dedxSObj) tuple->BS_NOMoNOHvsPV->Fill(goodVerts,dedxSObj->numberOfMeasurements()/(double)track->found(),Event_Weight);

   //Require at least one good vertex except if cosmic event
   //WAIT//if(TypeMode_==3 && goodVerts<1 && (!tuple || tuple->Name.find("Cosmic")==string::npos)) return false;

   //For TOF only analysis match to a SA track without vertex constraint for IP cuts
   if(TypeMode_==3) {

      //Find closest NV track
      const std::vector<reco::Track> noVertexTrackColl = iEvent.get(refittedStandAloneMuonsToken_);
      reco::Track NVTrack;
      double minDr=15;
      for(unsigned int i=0;i<noVertexTrackColl.size();i++){
         double dR = deltaR(track->eta(), track->phi(), noVertexTrackColl[i].eta(), noVertexTrackColl[i].phi());
         if(dR<minDr) {
            minDr=dR;
	         NVTrack=noVertexTrackColl[i];
         }
      }
      if(tuple) tuple->BS_dR_NVTrack->Fill(minDr,Event_Weight);
      if(minDr>0.4) return false;
      if(tuple)tuple->NVTrack->Fill(0.0,Event_Weight);

      //Find displacement of tracks with respect to beam spot
      const reco::BeamSpot beamSpotColl = iEvent.get(offlineBeamSpotToken_);

      dz  = NVTrack.dz (beamSpotColl.position());
      dxy = NVTrack.dxy(beamSpotColl.position());
      if(muonStations(NVTrack.hitPattern())<minMuStations) return false;
   }

   if(tuple){tuple->MTOF ->Fill(0.0,Event_Weight);
     if(GenBeta>=0)tuple->Beta_PreselectedB->Fill(GenBeta, Event_Weight);
   }

   double v3d = sqrt(dz*dz+dxy*dxy);

   if(tuple){tuple->BS_V3D->Fill(v3d,Event_Weight);}
   if(v3d>GlobalMaxV3D )return false;
   if(tuple){tuple->V3D  ->Fill(0.0,Event_Weight);}

   if(tuple)tuple->BS_Dxy->Fill(dxy, Event_Weight);

   TreeDXY = dxy;   
   bool DXYSB = false;
   if(TypeMode_!=5 && fabs(dxy)>GlobalMaxDXY)return false;
   if(TypeMode_==5 && fabs(dxy)>4)return false;
   if(TypeMode_==5 && fabs(dxy)>GlobalMaxDXY) DXYSB = true;

   if(tuple){tuple->Dxy  ->Fill(0.0,Event_Weight);}

   if(TypeMode_!=3 && Ih_Iso_cut) {
     const edm::ValueMap<susybsm::HSCPIsolation> IsolationMap = iEvent.get(hscpIsoToken_);

     susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
     if(tuple){tuple->BS_TIsol ->Fill(hscpIso.Get_TK_SumEt(),Event_Weight);}
      //     if(TypeMode_!=4){       if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;     }
      if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;
     if(tuple){tuple->TIsol   ->Fill(0.0,Event_Weight);}

     double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
     if(tuple){tuple->BS_EIsol ->Fill(EoP,Event_Weight);}
      //     if(TypeMode_!=4){       if(EoP>GlobalMaxEIsol)return false;     }
     if(EoP>GlobalMaxEIsol)return false;
     if(tuple){tuple->EIsol   ->Fill(0.0,Event_Weight);}
     
     // relative tracker isolation
     if (tuple) {  tuple->BS_SumpTOverpT->Fill(hscpIso.Get_TK_SumEt()/track->pt(), Event_Weight); }
      //     if(TypeMode_==4) { if(hscpIso.Get_TK_SumEt()/track->pt()>GlobalMaxRelTIsol)return false;   }
     if(hscpIso.Get_TK_SumEt()/track->pt()>GlobalMaxRelTIsol)return false;
     if (tuple) {  tuple->SumpTOverpT   ->Fill(0.0,Event_Weight);} 
   }

   if(tuple){tuple->BS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);}
   if(TypeMode_!=3 && (track->ptError()/track->pt())>GlobalMaxPterr)return false;
   //mk if(MassErr > 0 && MassErr > 2.2)return false; //FIXME jozze -- cut on relative mass error in units of 8*MassErr/Mass

   if(std::max(0.0,track->pt())<GlobalMinPt)return false;
   if(tuple){tuple->Pterr   ->Fill(0.0,Event_Weight);}

   //Find distance to nearest segment on opposite side of detector
   double minPhi=0.0, minEta=0.0;
   double segSep=SegSep(hscp, iEvent, minPhi, minEta);

   if(tuple){
     tuple->BS_SegSep->Fill(segSep, Event_Weight);
     tuple->BS_SegMinPhiSep->Fill(minPhi, Event_Weight);
     tuple->BS_SegMinEtaSep->Fill(minEta, Event_Weight);
     //Plotting segment separation depending on whether track passed dz cut
     if(fabs(dz)>GlobalMaxDZ) {
       tuple->BS_SegMinEtaSep_FailDz->Fill(minEta, Event_Weight);
     }
     else {
       tuple->BS_SegMinEtaSep_PassDz->Fill(minEta, Event_Weight);
     }
     //Plots for tracking failing Eta Sep cut
     if(fabs(minEta)<minSegEtaSep) {
       //Needed to compare dz distribution of cosmics in pure cosmic and main sample
       tuple->BS_Dz_FailSep->Fill(dz);
     }
   }



   //Now cut Eta separation
   //if(TypeMode_==3 && fabs(minEta)<minSegEtaSep) return false;
   //WAIT//if(tuple){tuple->SegSep->Fill(0.0,Event_Weight);}

   if(tuple && tof) {
     //Plots for tracks in dz control region
     if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz && !muon->isGlobalMuon()) {
       tuple->BS_Pt_FailDz->Fill(track->pt(), Event_Weight);
       tuple->BS_TOF_FailDz->Fill(tof->inverseBeta(), Event_Weight);
       if(fabs(track->eta())>CSCRegion) {
	 tuple->BS_TOF_FailDz_CSC->Fill(tof->inverseBeta(), Event_Weight);
	 tuple->BS_Pt_FailDz_CSC->Fill(track->pt(), Event_Weight);
       }
       else if(fabs(track->eta())<DTRegion) {
	 tuple->BS_TOF_FailDz_DT->Fill(tof->inverseBeta(), Event_Weight);
	 tuple->BS_Pt_FailDz_DT->Fill(track->pt(), Event_Weight);
       }
     }
     //Plots of dz
     tuple->BS_Dz->Fill(dz, Event_Weight);
     if(fabs(track->eta())>CSCRegion) tuple->BS_Dz_CSC->Fill(dz,Event_Weight);
     else if(fabs(track->eta())<DTRegion) tuple->BS_Dz_DT->Fill(dz,Event_Weight);
     tuple->BS_EtaDz->Fill(track->eta(),dz,Event_Weight);
   }


   //Split into different dz regions, each different region used to predict cosmic background and find systematic
   if(TypeMode_==3 && !muon->isGlobalMuon() && tuple) {
     //WAIT//int DzType=-1;
     //WAIT//if(fabs(dz)<GlobalMaxDZ) DzType=0;
     //WAIT//else if(fabs(dz)<30) DzType=1;
     //WAIT//else if(fabs(dz)<50) DzType=2;
     //WAIT//else if(fabs(dz)<70) DzType=3;
     //WAIT//if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz) DzType=4;
     //WAIT//if(fabs(dz)>CosmicMaxDz) DzType=5;

     //Count number of tracks in dz sidebands passing the TOF cut
     //The pt cut is not applied to increase statistics
     for(unsigned int CutIndex=0;CutIndex<CutPt_.size();CutIndex++){
       //WAIT//if(tof->inverseBeta()>=CutTOF_[CutIndex]) {
	    //WAIT//tuple->H_D_DzSidebands->Fill(CutIndex, DzType);
       //WAIT//}
     }
   }

   TreeDZ = dz;
   bool DZSB = false;
   if(TypeMode_!=5 && fabs(dz)>GlobalMaxDZ) return false;
   if(TypeMode_==5 && fabs(dz)>4) return false;
   if(TypeMode_==5 && fabs(dz)>GlobalMaxDZ) DZSB = true;
   if(tuple){tuple->Dz  ->Fill(0.0,Event_Weight);}

   if(TypeMode_==3 && fabs(minEta)<minSegEtaSep) return false;
   if(tuple)tuple->BS_Phi->Fill(track->phi(),Event_Weight);
   if(TypeMode_==3 && fabs(track->phi())>1.2 && fabs(track->phi())<1.9) return false;

    //skip HSCP that are compatible with cosmics.
    if(tuple)tuple->BS_OpenAngle->Fill(OpenAngle,Event_Weight);

    bool OASB = false;
    if(TypeMode_==5 && OpenAngle>=2.8)OASB = true;

   isCosmicSB = DXYSB && DZSB && OASB;
   isSemiCosmicSB = (!isCosmicSB && (DXYSB || DZSB || OASB));
 
   if(tuple){if(dedxSObj) tuple->BS_EtaIs->Fill(track->eta(),dedxSObj->dEdx(),Event_Weight);
          if(dedxMObj) tuple->BS_EtaIm->Fill(track->eta(),dedxMObj->dEdx(),Event_Weight);
          tuple->BS_EtaP ->Fill(track->eta(),track->p(),Event_Weight);
          tuple->BS_EtaPt->Fill(track->eta(),track->pt(),Event_Weight);
          if(tof) tuple->BS_EtaTOF->Fill(track->eta(),tof->inverseBeta(),Event_Weight);
   }

   if(tuple){if(GenBeta>=0)tuple->Beta_PreselectedC->Fill(GenBeta, Event_Weight);
          if(DZSB  && OASB)tuple->BS_Dxy_Cosmic->Fill(dxy, Event_Weight);
          if(DXYSB && OASB)tuple->BS_Dz_Cosmic->Fill(dz, Event_Weight);
          if(DXYSB && DZSB)tuple->BS_OpenAngle_Cosmic->Fill(OpenAngle,Event_Weight);


          //WAIT// 
          TVector3 outerHit = getOuterHitPos(dedxHits);
          TVector3 vertex(vertexColl[highestPtGoodVertex].position().x(), vertexColl[highestPtGoodVertex].position().y(), vertexColl[highestPtGoodVertex].position().z());
          tuple->BS_LastHitDXY  ->Fill((outerHit).Perp(),Event_Weight);
          tuple->BS_LastHitD3D  ->Fill((outerHit).Mag(),Event_Weight);

          tuple->BS_P  ->Fill(track->p(),Event_Weight);
          tuple->BS_Pt ->Fill(track->pt(),Event_Weight);
          if(PUA)tuple->BS_Pt_PUA ->Fill(track->pt(),Event_Weight);
          if(PUB)tuple->BS_Pt_PUB ->Fill(track->pt(),Event_Weight);
          if(DXYSB && DZSB && OASB) tuple->BS_Pt_Cosmic->Fill(track->pt(),Event_Weight);

	  if(fabs(track->eta())<DTRegion) tuple->BS_Pt_DT->Fill(track->pt(),Event_Weight);
	  else tuple->BS_Pt_CSC->Fill(track->pt(),Event_Weight);

          double RecoQoPt = track->charge()/track->pt();
          if(!hscp.trackRef().isNull() && hscp.trackRef()->pt()>200) {
            double InnerRecoQoPt = hscp.trackRef()->charge()/hscp.trackRef()->pt();
            tuple->BS_InnerInvPtDiff->Fill((RecoQoPt-InnerRecoQoPt)/InnerRecoQoPt,Event_Weight);
          }

          if(dedxSObj) tuple->BS_Is ->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj && PUA) tuple->BS_Is_PUA ->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj && PUB) tuple->BS_Is_PUB ->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj && DXYSB && DZSB && OASB) tuple->BS_Is_Cosmic->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxMObj) tuple->BS_Im ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(dedxMObj && PUA) tuple->BS_Im_PUA ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(dedxMObj && PUB) tuple->BS_Im_PUB ->Fill(dedxMObj->dEdx(),Event_Weight);

          if(tof)
          {
	        tuple->BS_TOF->Fill(tof->inverseBeta(),Event_Weight);
            if(PUA)tuple->BS_TOF_PUA->Fill(tof->inverseBeta(),Event_Weight);
            if(PUB)tuple->BS_TOF_PUB->Fill(tof->inverseBeta(),Event_Weight);
	        if(dttof->nDof()>6) tuple->BS_TOF_DT->Fill(dttof->inverseBeta(),Event_Weight);
            if(csctof->nDof()>6) tuple->BS_TOF_CSC->Fill(csctof->inverseBeta(),Event_Weight);
            tuple->BS_PtTOF->Fill(track->pt() ,tof->inverseBeta(),Event_Weight);
          }
          if(dedxSObj && dedxMObj) {
	        tuple->BS_PIs  ->Fill(track->p()  ,dedxSObj->dEdx(),Event_Weight);
            tuple->BS_PImHD->Fill(track->p()  ,dedxMObj->dEdx(),Event_Weight);
            tuple->BS_PIm  ->Fill(track->p()  ,dedxMObj->dEdx(),Event_Weight);
            tuple->BS_PtIs ->Fill(track->pt() ,dedxSObj->dEdx(),Event_Weight);
            tuple->BS_PtIm ->Fill(track->pt() ,dedxMObj->dEdx(),Event_Weight);
	  }
          if(dedxSObj && tof)tuple->BS_TOFIs->Fill(tof->inverseBeta(),dedxSObj->dEdx(),Event_Weight);
          if(dedxMObj && tof)tuple->BS_TOFIm->Fill(tof->inverseBeta(),dedxMObj->dEdx(),Event_Weight);

	  //Muon only prediction binned depending on where in the detector the track is and how many muon stations it has
	  //Binning not used for other analyses
	  int bin=-1;
	  if(TypeMode_==3) {
	    if(fabs(track->eta())<DTRegion) bin=muonStations(track->hitPattern())-2;
	    else bin=muonStations(track->hitPattern())+1;
	    tuple->BS_Pt_Binned[to_string(bin)] ->Fill(track->pt(),Event_Weight);
	  }
   }
   if(tuple){tuple->Basic  ->Fill(0.0,Event_Weight);}

   return true;
}

//=============================================================
//
//     Selection
//
//=============================================================
bool Analyzer::passSelection(
         const susybsm::HSCParticle& hscp,  
         const reco::DeDxData* dedxSObj, 
         const reco::DeDxData* dedxMObj, 
         const reco::MuonTimeExtra* tof, 
         const edm::Event& iEvent,
         float Event_Weight,
         const int& CutIndex, 
         Tuple* &tuple, 
         const bool isFlip, 
         const double& GenBeta, 
         bool RescaleP, 
         const double& RescaleI, 
         const double& RescaleT)
{
   reco::TrackRef   track;
   if(TypeMode_!=3) track = hscp.trackRef();
   else {
     reco::MuonRef muon = hscp.muonRef();
     if(muon.isNull()) return false;
     track = muon->standAloneMuon();
   }
   if(track.isNull())return false;

   double MuonTOF = GlobalMinTOF;
   if(tof){
      MuonTOF = tof->inverseBeta();
   }

   double Is=0;   if(dedxSObj) Is=dedxSObj->dEdx();
   double Ih=0;   if(dedxMObj) Ih=dedxMObj->dEdx();
   //WAIT//double Ick=0; // if(dedxMObj) Ick=GetIck(Ih,isBckg);

   double PtCut=CutPt_[CutIndex];
   double ICut=CutI_[CutIndex];
   double TOFCut=CutTOF_[CutIndex];
   if(isFlip) {
     PtCut=CutPt_Flip_[CutIndex];
     ICut=CutI_Flip_[CutIndex];
     TOFCut=CutTOF_Flip_[CutIndex];
   }

   if(RescaleP){
     if(RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())<PtCut)return false;
     //if(std::max(0.0,RescaledPt(track->pt() - track->ptError(),track->eta(),track->phi(),track->charge()))<CutPt_[CutIndex])return false; 
   }else{
     if(track->pt()<PtCut)return false;
     //if(std::max(0.0,(track->pt() - track->ptError()))<CutPt_[CutIndex])return false;
   } 
   if(tuple){tuple->Pt    ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)tuple->Beta_SelectedP->Fill(CutIndex,GenBeta, Event_Weight);
   }

   if(TypeMode_!=3 && Is+RescaleI<ICut)return false;

   if(tuple){tuple->I    ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)tuple->Beta_SelectedI->Fill(CutIndex, GenBeta, Event_Weight);
   }

   if((TypeMode_>1  && TypeMode_!=5) && !isFlip && MuonTOF+RescaleT<TOFCut)return false;
   if((TypeMode_>1  && TypeMode_!=5) && isFlip && MuonTOF+RescaleT>TOFCut)return false;

   if(tuple){tuple->TOF  ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)tuple->Beta_SelectedT->Fill(CutIndex, GenBeta, Event_Weight);
          tuple->AS_P  ->Fill(CutIndex,track->p(),Event_Weight);
          tuple->AS_Pt ->Fill(CutIndex,track->pt(),Event_Weight);
          tuple->AS_Is ->Fill(CutIndex,Is,Event_Weight);
          tuple->AS_Im ->Fill(CutIndex,Ih,Event_Weight);
          tuple->AS_TOF->Fill(CutIndex,MuonTOF,Event_Weight);
          //tuple->AS_EtaIs->Fill(CutIndex,track->eta(),Is,Event_Weight);
          //tuple->AS_EtaIm->Fill(CutIndex,track->eta(),Ih,Event_Weight);
          //tuple->AS_EtaP ->Fill(CutIndex,track->eta(),track->p(),Event_Weight);
          //tuple->AS_EtaPt->Fill(CutIndex,track->eta(),track->pt(),Event_Weight);
          tuple->AS_PIs  ->Fill(CutIndex,track->p()  ,Is,Event_Weight);
          tuple->AS_PIm  ->Fill(CutIndex,track->p()  ,Ih,Event_Weight);
          tuple->AS_PtIs ->Fill(CutIndex,track->pt() ,Is,Event_Weight);
          tuple->AS_PtIm ->Fill(CutIndex,track->pt() ,Ih,Event_Weight);
          tuple->AS_TOFIs->Fill(CutIndex,MuonTOF     ,Is,Event_Weight);
          tuple->AS_TOFIm->Fill(CutIndex,MuonTOF     ,Ih,Event_Weight);
   }
   return true;      
}

//=============================================================
//
//     Trigger-Selection
//
//=============================================================
bool Analyzer::passTrigger(const edm::Event& iEvent, bool isData, bool isCosmic, L1BugEmulator* emul){

   edm::Handle<edm::TriggerResults> triggerH;
   iEvent.getByToken(triggerResultsToken_,triggerH); 
   bool valid = triggerH.isValid(); 
   if (not valid){
      edm::LogError("Analyzer") << "HLT TriggerResults not found!";
      return false;
   }

   const edm::TriggerNames& triggerNames = iEvent.triggerNames(*triggerH);

   bool metTrig = PassTriggerPatterns(triggerH, triggerNames, trigger_met_);
   bool muTrig  = PassTriggerPatterns(triggerH, triggerNames, trigger_mu_);

   if (!metTrig && muTrig) TrigInfo_ = 1; // mu only
   if (metTrig && !muTrig) TrigInfo_ = 2; // met only
   if (metTrig && muTrig)  TrigInfo_ = 3; // mu and met*/

   if (metTrig) return true;
   if (muTrig){
      if (!isData && emul){
         edm::Handle < vector<reco::Muon> > muonCollH;
         iEvent.getByToken(muonToken_,muonCollH);
         if(!muonCollH.isValid()) return false;
         bool KeepEvent=false;
         for (unsigned int c=0;c<muonCollH->size();c++){
            reco::MuonRef muon = reco::MuonRef(muonCollH.product(), c);
            if (muon.isNull()) continue;
            if (muon->track().isNull()) continue; 
            if (emul->PassesL1Inefficiency(muon->track()->pt(), std::fabs(muon->track()->eta()))){
               KeepEvent=true;
               break;
            }
         }
         return KeepEvent;
      }
      else return true;
   }

   return false; //FIXME triggers bellow will need to be adapted based on Run2 trigger menu
   /*
   //for(unsigned int i=0;i<tr.size();i++){
   //printf("Path %3i %50s --> %1i\n",i, tr.triggerName(i).c_str(),tr.accept(i));
   //}fflush(stdout);

   //if(tr.accept("HSCPHLTTriggerMetDeDxFilter"))return true;
   //if(tr.accept("HSCPHLTTriggerMuDeDxFilter"))return true;
   if(tr.accept("HSCPHLTTriggerMuFilter"))return true;
   if(tr.accept("HSCPHLTTriggerPFMetFilter"))return true;

   //Could probably use this trigger for the other analyses as well
   if(TypeMode_==3){
      if(tr.size()== tr.triggerIndex("HSCPHLTTriggerL2MuFilter")) return false;
      if(tr.accept(tr.triggerIndex("HSCPHLTTriggerL2MuFilter")))  return true;

      //Only accepted if looking for cosmic events
      if(isCosmic) {
         if(tr.size()== tr.triggerIndex("HSCPHLTTriggerCosmicFilter")) return false;
         if(tr.accept(tr.triggerIndex("HSCPHLTTriggerCosmicFilter"))) return true;
      }
   }*/
   return false;
}

//=============================================================
//
//     Check if track is from pixel
//
//=============================================================
void 
Analyzer::isPixelTrack(const edm::Ref<std::vector<Trajectory> > &refTraj, bool &isBpixtrack, bool &isFpixtrack) {
    // Used in analyze() to see if it is pixel track
    std::vector<TrajectoryMeasurement> tmeasColl = refTraj->measurements();
    std::vector<TrajectoryMeasurement>::const_iterator tmeasIt;
    for (tmeasIt = tmeasColl.begin(); tmeasIt != tmeasColl.end(); tmeasIt++) {
        if (!tmeasIt->updatedState().isValid()) continue;
        TransientTrackingRecHit::ConstRecHitPointer testhit = tmeasIt->recHit();
        if (!testhit->isValid() || testhit->geographicalId().det() != DetId::Tracker) continue;
        uint testSubDetID = (testhit->geographicalId().subdetId());
        if (testSubDetID==PixelSubdetector::PixelBarrel) isBpixtrack = true;
        if (testSubDetID==PixelSubdetector::PixelEndcap) isFpixtrack = true;
        if (isBpixtrack && isFpixtrack) break;
    }
}
