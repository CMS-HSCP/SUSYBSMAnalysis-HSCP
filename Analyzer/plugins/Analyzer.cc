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
   ,refittedStandAloneMuonsToken_(consumes<vector<reco::Track>>(iConfig.getParameter<edm::InputTag>("refittedStandAloneMuonsCollection")))
   ,offlineBeamSpotToken_(consumes<reco::BeamSpot>(iConfig.getParameter<edm::InputTag>("offlineBeamSpotCollection")))
   // Parameters
   ,Debug(iConfig.getUntrackedParameter<bool>("Debug"))
   ,TypeMode(iConfig.getUntrackedParameter<unsigned int>("TypeMode"))
   ,SampleType(iConfig.getUntrackedParameter<unsigned int>("SampleType"))
   ,SkipSelectionPlot(iConfig.getUntrackedParameter<bool>("SkipSelectionPlot"))
   ,PtHistoUpperBound(iConfig.getUntrackedParameter<double>("PtHistoUpperBound"))
   ,MassHistoUpperBound(iConfig.getUntrackedParameter<double>("MassHistoUpperBound"))
   ,MassNBins(iConfig.getUntrackedParameter<unsigned int>("MassNBins"))
   ,IPbound(iConfig.getUntrackedParameter<double>("IPbound"))
   ,PredBins(iConfig.getUntrackedParameter<unsigned int>("PredBins"))
   ,EtaBins(iConfig.getUntrackedParameter<unsigned int>("EtaBins"))
   ,dEdxS_UpLim(iConfig.getUntrackedParameter<double>("dEdxS_UpLim"))
   ,dEdxM_UpLim(iConfig.getUntrackedParameter<double>("dEdxM_UpLim"))
   ,DzRegions(iConfig.getUntrackedParameter<unsigned int>("DzRegions"))
   ,GlobalMinPt(iConfig.getUntrackedParameter<double>("GlobalMinPt"))
   ,GlobalMinTOF(iConfig.getUntrackedParameter<double>("GlobalMinTOF"))
   ,DeDxSF_0(iConfig.getUntrackedParameter<double>("DeDxSF_0"))
   ,DeDxSF_1(iConfig.getUntrackedParameter<double>("DeDxSF_1"))
   ,DeDxK(iConfig.getUntrackedParameter<double>("DeDxK"))
   ,DeDxC(iConfig.getUntrackedParameter<double>("DeDxC"))
   ,DeDxTemplate(iConfig.getUntrackedParameter<string>("DeDxTemplate"))
   ,DeDxCalibration(iConfig.getUntrackedParameter<string>("DeDxCalibration"))
   ,Geometry(iConfig.getUntrackedParameter<string>("Geometry"))
   ,TimeOffset(iConfig.getUntrackedParameter<string>("TimeOffset"))

{
   //now do what ever initialization is needed
   // define the selection to be considered later for the optimization
   // WARNING: recall that this has a huge impact on the analysis time AND on the output file size --> be carefull with your choice
   
   useClusterCleaning = true;
   if(TypeMode==4) {
      useClusterCleaning = false; //switch off cluster cleaning for mCHAMPs
   }
   
   isData   = (SampleType==0);
   isMC     = (SampleType==1);
   isSignal = (SampleType>=2);

   dEdxSF [0] = DeDxSF_0;
   dEdxSF [1] = DeDxSF_1;

   dEdxK_Data    = DeDxK;
   dEdxC_Data    = DeDxC;
   dEdxK_MC      = DeDxK;
   dEdxC_MC      = DeDxC;

   dEdxTemplates = loadDeDxTemplate(DeDxTemplate, true);
   if(isData)   trackerCorrector.LoadDeDxCalibration(DeDxCalibration); 
   else         trackerCorrector.TrackerGains = nullptr; //FIXME check gain for MC

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
   
   string BaseName;
   if(isData)
        BaseName = "Data";
   else BaseName = "MC";

   TFileDirectory dir = fs->mkdir( BaseName.c_str(), BaseName.c_str() );

   // create histograms & trees
   initializeCuts(fs);
   //initializeTuple(tuple_saver, dir);
   tuple_saver->initializeTuple(tuple, dir, SkipSelectionPlot, TypeMode, isSignal, CutPt.size(), CutPt_Flip.size(), PtHistoUpperBound, MassHistoUpperBound, MassNBins, IPbound, PredBins, EtaBins, dEdxS_UpLim, dEdxM_UpLim, DzRegions, GlobalMinPt, GlobalMinTOF);

   CurrentRun = 0;
   RNG = new TRandom3();
   is2016 = false;
   is2016G = false;

   


   HSCPTk              = new bool[CutPt.size()];
   HSCPTk_SystP        = new bool[CutPt.size()];
   HSCPTk_SystI        = new bool[CutPt.size()];
   HSCPTk_SystT        = new bool[CutPt.size()];
   HSCPTk_SystM        = new bool[CutPt.size()];
   HSCPTk_SystPU       = new bool[CutPt.size()];
   HSCPTk_SystHUp      = new bool[CutPt.size()];
   HSCPTk_SystHDown    = new bool[CutPt.size()];
   MaxMass           = new double[CutPt.size()];
   MaxMass_SystP     = new double[CutPt.size()];
   MaxMass_SystI     = new double[CutPt.size()];
   MaxMass_SystT     = new double[CutPt.size()];
   MaxMass_SystM     = new double[CutPt.size()];
   MaxMass_SystPU    = new double[CutPt.size()];
   MaxMass_SystHUp   = new double[CutPt.size()];
   MaxMass_SystHDown = new double[CutPt.size()];

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

   //TO BE TESTED// isData = iEvent.isRealData();

   if(CurrentRun != iEvent.id().run()){ //removed iEvent.eventAuxiliary().run() call
      CurrentRun  = iEvent.id().run();
      tofCalculator.setRun(CurrentRun);
      //trackerCorrector.setRun(CurrentRun);
   }

   //===================== Handle For DeDx Hits ==============
   Handle<reco::DeDxHitInfoAss> dedxCollH;
   iEvent.getByToken(dedxToken_,dedxCollH);
   if(!dedxCollH.isValid()){printf("Invalid dedxCollH\n");return;}

   //================= Handle For Muon TOF Combined ===============
   //Handle<ValueMap<reco::MuonTimeExtra>> TOFCollH; Handle<reco::MuonTimeExtraMap> tofMap
   Handle<reco::MuonTimeExtraMap>     TOFCollH;
   iEvent.getByToken(muonTimeToken_,  TOFCollH);
   const reco::MuonTimeExtraMap & tofMap = *TOFCollH;

   //================= Handle For Muon TOF DT ===============
   Handle<reco::MuonTimeExtraMap>      TOFDTCollH;
   iEvent.getByToken(muonDtTimeToken_, TOFDTCollH);
   const reco::MuonTimeExtraMap & tofDtMap = *TOFDTCollH;

   //================= Handle For Muon TOF CSC ===============
   Handle<reco::MuonTimeExtraMap>        TOFCSCCollH;
   iEvent.getByToken(muonCscTimeToken_,  TOFCSCCollH);
   const reco::MuonTimeExtraMap & tofCscMap = *TOFCSCCollH;

   //================= Handle For Muon DT/CSC Segment ===============
   Handle<CSCSegmentCollection> CSCSegmentCollH;
   Handle<DTRecSegment4DCollection> DTSegmentCollH;
   if(!isMC){ //do not recompute TOF on MC background
      iEvent.getByToken(muonCscSegmentToken_, CSCSegmentCollH);
      if(!CSCSegmentCollH.isValid()){printf("CSC Segment Collection not found!\n"); return;}

      iEvent.getByToken(muonDtSegmentToken_, DTSegmentCollH);
      if(!DTSegmentCollH.isValid()){printf("DT Segment Collection not found!\n"); return;}
   }

   //reinitialize the bookeeping array for each event
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk        [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystP  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystI  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystT  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystM  [CutIndex] = false;   }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystPU [CutIndex] = false; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystHUp[CutIndex] = false; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  HSCPTk_SystHDown[CutIndex] = false; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass       [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystP [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystI [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystT [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystM [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystPU[CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystHUp [CutIndex] = -1; }
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){  MaxMass_SystHDown[CutIndex] = -1; }

   //WAIT//HIPemulator.setEventRate(); //take it from a pdf
   //WAIT//HIPemulatorUp.setEventRate(HIPemulator.getEventRatePixel()*1.25, HIPemulator.getEventRateStrip()*1.80);  // deltaPixel = 3.653981e+02, basePixel = 1.332625e+03; deltaStrip = 4.662832e+02, baseStrip = 5.958308e+02, from Run257805
   //WAIT//HIPemulatorDown.setEventRate(HIPemulator.getEventRatePixel()*0.75, HIPemulator.getEventRateStrip()*0.20); 

	//WAIT//HIPTrackLossEmul.SetHIPTrackLossRate(iEvent);

   vector<reco::GenParticle> genColl;
   /*//double HSCPGenBeta1=-1, HSCPGenBeta2=-1;
   //double HSCPDLength1=-1, HSCPDLength2=-1;

   //if(isSignal){}*/

   //load all event collection that will be used later on (HSCP, dEdx and TOF)
   //====================loop over HSCP candidates===================
   unsigned int count = 0;
   for(const auto& hscp : iEvent.get(hscpToken_)){
      reco::MuonRef  muon  = hscp.muonRef();//const reco::MuonRef& muon = hscp.muonRef();

      //For TOF only analysis use updated stand alone muon track.
	   //Otherwise use inner tracker track
	   reco::TrackRef track;
      if(TypeMode!=3) track = hscp.trackRef();
      else {
		   if(muon.isNull()) continue;
		   track = muon->standAloneMuon();
      }
      //skip events without track
	   if(track.isNull())continue;
      // FIXME jozze skip events with |Eta| > 0.9 (out of the barrel)
	   //if(track->eta()>0.9 || track->eta() < -0.9) continue;

      //require a track segment in the muon system
      if(TypeMode>1 && TypeMode!=5 && (muon.isNull() || !muon->isStandAloneMuon()))continue;

      //Apply a scale factor to muon only analysis to account for differences seen in data/MC preselection efficiency
      //For eta regions where Data > MC no correction to be conservative
      if(!isData && TypeMode==3 && scaleFactor(track->eta())<RNG->Uniform(0, 1)) continue;

      //for signal only, make sure that the candidate is associated to a true HSCP
      int ClosestGen;
      if(isSignal && DistToHSCP(hscp, genColl, ClosestGen, TypeMode)>0.03)continue;

      // we are losing some tracks due to HIP
	   //WAIT//if(!isData && is2016 && !HIPTrackLossEmul.TrackSurvivesHIPInefficiency()) continue;

      //load quantity associated to this track (TOF and dEdx)
      const reco::DeDxHitInfo* dedxHits = nullptr;//const reco::DeDxHitInfo* dedxHits = nullptr;
      if(TypeMode!=3 && !track.isNull()) {
         reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
         if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
      }
      //bool useTOF = false;
      if(TypeMode>1 && TypeMode!=5 && !hscp.muonRef().isNull()){
         //useTOF = true;
         if(isMC){
            /*tof    = &TOFCollH->get(hscp.muonRef().key());
            dttof  = &_tofDtMap->get(hscp.muonRef().key());  
            csctof = &_tofCscMap->get(hscp.muonRef().key());*/
            tof    = tofMap[hscp.muonRef()];
            dttof  = tofDtMap[hscp.muonRef()];
            csctof = tofCscMap[hscp.muonRef()];
         }else{
            const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollH;
            const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollH;
            //std::cout<<"TESTA\n";
            tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, isData?1:0 ); //apply T0 correction on data but not on signal MC
            //std::cout<<"TESTB\n";
            tof    = tofCalculator.combinedTOF; 
            dttof  = tofCalculator.dtTOF;  
            csctof = tofCalculator.cscTOF;
            //std::cout<<"TESTC\n";
         }
      }

      if(!dedxHits) continue; // skip tracks without hits otherwise there will be a crash

      HitDeDxCollection hitDeDx = getHitDeDx(dedxHits, dEdxSF, trackerCorrector.TrackerGains, false, 1);

      unsigned int pdgId = 0;
      if(isSignal){ 
         pdgId = genColl[ClosestGen].pdgId();
         if (Debug) cout << "GenId  " << pdgId << endl;
      }

      double dEdxErr = 0;
      reco::DeDxData dedxSObjTmp  = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, TypeMode==5, false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.00, nullptr,0,pdgId);
      reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr, &dEdxErr,pdgId);
      reco::DeDxData dedxMUpObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr,0,pdgId);
      reco::DeDxData dedxMDownObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, nullptr,0,pdgId);
      /*reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulator:nullptr, &dEdxErr,pdgId);
      reco::DeDxData dedxMUpObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulatorUp:nullptr,0,pdgId);
      reco::DeDxData dedxMDownObjTmp = computedEdx(dedxHits, dEdxSF, nullptr,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulatorDown:nullptr,0,pdgId);*/
      reco::DeDxData* dedxSObj  = dedxSObjTmp.numberOfMeasurements()>0?&dedxSObjTmp:nullptr;
      reco::DeDxData* dedxMObj  = dedxMObjTmp.numberOfMeasurements()>0?&dedxMObjTmp:nullptr;
      reco::DeDxData* dedxMUpObj = dedxMUpObjTmp.numberOfMeasurements()>0?&dedxMUpObjTmp:nullptr;
      reco::DeDxData* dedxMDownObj = dedxMDownObjTmp.numberOfMeasurements()>0?&dedxMDownObjTmp:nullptr;
      if(TypeMode==5)OpenAngle = deltaROpositeTrack(iEvent.get(hscpToken_), hscp); //OpenAngle is a global variable... that's uggly C++, but that's the best I found so far

      //compute systematic uncertainties on signal
      if(isSignal){
         //FIXME to be measured on 2015 data, currently assume 2012
         /*bool   PRescale = true;
         double IRescale =-0.05; // added to the Ias value
         double MRescale = 0.95;
		   double TRescale =-0.015; //-0.005 (used in 2012); // added to the 1/beta value*/
		  
		   double genpT = -1.0;
		   for(unsigned int g=0;g<genColl.size();g++) {
            if(genColl[g].pt()<5)continue;
            if(genColl[g].status()!=1)continue;
            int AbsPdg=abs(genColl[g].pdgId());
            if(AbsPdg!=17)continue;
            
            double separation = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
            if (separation > 0.03) continue;
            genpT = genColl[g].pt();
            break;
         }
         if (genpT>0) {  tuple->genrecopT->Fill(genpT, track->pt()); }

         // compute systematic due to momentum scale
         //WAIT//if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   PRescale, 0, 0)){..}

         // compute systematic due to dEdx (both Ias and Ih)
         //WAIT//if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, IRescale, 0)){...}

         // compute systematic due to Mass shift
         //WAIT//if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, 0, 0)){...}

         // compute systematic due to TOF
         //WAIT//if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, 0, TRescale)){...}

         // compute systematics due to PU

      }//End of systematic computation for signal

      //check if the canddiate pass the preselection cuts
      /*const susybsm::HSCParticle& hscp, const DeDxHitInfo* dedxHits,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const ChainEvent& ev, stPlots* st, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT, double MassErr*/
      double MassErr=GetMassErr (track->p(), track->ptError(), dedxMObj?dedxMObj->dEdx():-1, dEdxErr, GetMass(track->p(), dedxMObj?dedxMObj->dEdx():-1, DeDxK,DeDxC), DeDxK,DeDxC);
      if(isMC){
         passPreselection(  hscp, dedxHits, dedxSObj, dedxMObj, iEvent, tuple, -1, false, 0, 0, MassErr );
      }
      if(!passPreselection( hscp, dedxHits, dedxSObj, dedxMObj, iEvent, tuple, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1, false, 0, 0, MassErr) ) continue;
      if(TypeMode==5 && isSemiCosmicSB)continue;

      //fill the ABCD histograms and a few other control plots
      //WAIT//if(isData)Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, SamplePlots);
      //WAIT//else if(isMC) Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, MCTrPlots);
      if(TypeMode==5 && isCosmicSB)continue; 

      //Find the number of tracks passing selection for TOF<1 that will be used to check the background prediction
      double Mass = -1;
      if(isMC || isData) {
         //compute the mass of the candidate, for TOF mass flip the TOF over 1 to get the mass, so 0.8->1.2
		   Mass = GetMass(track->p(),dedxMObj->dEdx(),DeDxK,DeDxC);
		   //double MassTOF  = -1; if(useTOF)MassTOF = GetTOFMass(track->p(),(2-1));//GetTOFMass(track->p(),(2-tof.inverseBeta()));
		   //double MassComb = -1;
		   ///if(useTOF && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),DeDxK,DeDxC) + (1/(2-tof.inverseBeta())))*0.5 ) ;
		   ///if(dedxMObj) MassComb = Mass;
		   ///if(useTOF)MassComb=GetMassFromBeta(track->p(),(1/(2-tof.inverseBeta())));

         for(unsigned int CutIndex=0;CutIndex<CutPt_Flip.size();CutIndex++){
            //Background check looking at region with TOF<1
            //WAIT//if(!PassSelection   (hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, true)) continue;

            //Fill Mass Histograms
            tuple->Mass_Flip->Fill(CutIndex, Mass,Event_Weight);
            ///if(tof){
            ///   MassTOF_Flip->Fill(CutIndex, MassTOF, Event_Weight);
            ///}
            ///MassComb_Flip->Fill(CutIndex, MassComb, Event_Weight);
		   }

      }

      double Ick2=0;  if(dedxMObj) Ick2=GetIck(dedxMObj->dEdx(),isMC,DeDxK,DeDxC);
      int nomh= 0;nomh = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
      double fovhd = track->found()<=0?-1:track->found() / float(track->found() + nomh);
      unsigned int nom=0; if(dedxSObj) nom=dedxSObj->numberOfMeasurements();

      double weight=0,genid=0,gencharge=-99,genmass=-99,genpt=-99,geneta=-99,genphi=-99;
      weight = Event_Weight;
  
      if(isSignal){
         genid = genColl[ClosestGen].pdgId();
         gencharge = genColl[ClosestGen].charge();
         genmass = genColl[ClosestGen].mass();
         genpt = genColl[ClosestGen].pt();
         geneta = genColl[ClosestGen].eta();
         genphi = genColl[ClosestGen].phi();
      }

      tuple_saver->fillTreeBranches(tuple,
         TrigInfo, iEvent.id().run(),iEvent.id().event(),iEvent.id().luminosityBlock(), 
         count, track->charge(), track->pt(),track->ptError(), 
         dedxSObj ? dedxSObj->dEdx() : -1,
         dedxSObj ? dedxMObj->dEdx() : -1,
         dedxMObj ? Ick2 : -99, 
         tof.inverseBeta(), //tof ? tof.inverseBeta() : -1, 
         Mass, TreeDZ, TreeDXY, OpenAngle, 
         track->eta(), track->phi(), track->found(), track->hitPattern().numberOfValidPixelHits(), track->validFraction(), 
         nomh,fovhd, nom, weight,genid,gencharge,genmass,genpt,geneta,genphi
      );
      // Save in the tree
      if (!SkipSelectionPlot) tuple->Tree->Fill();
   count++;
   } //END loop over HSCP candidates

   //save event dependent information thanks to the bookkeeping
   for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
      if(HSCPTk[CutIndex]){
         tuple->HSCPE             ->Fill(CutIndex,Event_Weight);
         tuple->MaxEventMass      ->Fill(CutIndex,MaxMass[CutIndex], Event_Weight);
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
{}

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
void Analyzer::initializeCuts(edm::Service<TFileService> &fs){
   CutPt.clear();       CutI.clear();       CutTOF.clear();      
   CutPt_Flip.clear();  CutI_Flip.clear();  CutTOF_Flip.clear();  
   
   CutPt     .push_back(GlobalMinPt);   CutI       .push_back(GlobalMinIs);  CutTOF     .push_back(GlobalMinTOF);
   CutPt_Flip.push_back(GlobalMinPt);   CutI_Flip  .push_back(GlobalMinIs);  CutTOF_Flip.push_back(GlobalMinTOF);

   if(TypeMode<2){   
      for(double Pt =GlobalMinPt+5 ; Pt <200;Pt+=5){
         for(double I  =GlobalMinIs+0.025  ; I  <0.45 ;I+=0.025){
            CutPt .push_back(Pt);   CutI  .push_back(I);  CutTOF.push_back(-1);
         }
      }
   }else if(TypeMode==2){
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
   }else if(TypeMode==3){
      for(double Pt =GlobalMinPt+30 ; Pt <450;  Pt+=30){
      for(double TOF=GlobalMinTOF+0.025; TOF<1.5;TOF+=0.025){
         CutPt .push_back(Pt);   CutI  .push_back(-1);  CutTOF.push_back(TOF);
      }}
      for(double Pt =GlobalMinPt+30 ; Pt <450;  Pt+=60){
      for(double TOF=GlobalMinTOF-0.025; TOF>0.5;TOF-=0.025){
         CutPt_Flip .push_back(Pt);   CutI_Flip  .push_back(-1);  CutTOF_Flip.push_back(TOF);
      }}
   }else if(TypeMode==4){
      for(double I  =GlobalMinIs +0.025; I  <0.55;  I+=0.025){
      for(double TOF=GlobalMinTOF+0.025; TOF<1.46;TOF+=0.025){
 	 CutPt .push_back(-1);   CutI  .push_back(I);  CutTOF.push_back(TOF);
       }}
      for(double I  =GlobalMinIs +0.025; I  <0.55;  I+=0.025){
      for(double TOF=GlobalMinTOF-0.025; TOF>0.54;TOF-=0.025){
	 CutPt_Flip .push_back(-1);   CutI_Flip  .push_back(I);  CutTOF_Flip.push_back(TOF);
       }}
   }else if(TypeMode==5){   
      for(double Pt =75 ; Pt <=150;Pt+=25){
      for(double I  =0.0; I  <=0.45 ;I+=0.025){
         CutPt     .push_back(Pt);   CutI     .push_back(I);  CutTOF     .push_back(-1);
         CutPt_Flip.push_back(Pt);   CutI_Flip.push_back(I);  CutTOF_Flip.push_back(-1);
     }}
   }

   printf("%i Different Final Selection will be tested\n",(int)CutPt.size());
   printf("%i Different Final Selection will be tested for background uncertainty\n",(int)CutPt_Flip.size());

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

//Counts the number of muon stations used in track fit only counting DT and CSC stations.
int  Analyzer::muonStations(const reco::HitPattern& hitPattern) {
  int stations[4] = { 0,0,0,0 };
  for (int i=0; i<hitPattern.numberOfAllHits(reco::HitPattern::HitCategory::TRACK_HITS); i++) {
    uint32_t pattern = hitPattern.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS, i );
    if(pattern == 0) break;
    if(hitPattern.muonHitFilter(pattern) && (int(hitPattern.getSubStructure(pattern)) == 1 || int(hitPattern.getSubStructure(pattern)) == 2) && hitPattern.getHitType(pattern) == 0){
      stations[hitPattern.getMuonStation(pattern)-1] = 1;
    }
  }
  return stations[0]+stations[1]+stations[2]+stations[3];

}

double Analyzer::RescaledPt(const double& pt, const double& eta, const double& phi, const int& charge){
  if(TypeMode!=3) {
    double newInvPt = 1/pt+0.000236-0.000135*pow(eta,2)+charge*0.000282*TMath::Sin(phi-1.337);
    return 1/newInvPt;
  }
  else {
    double newInvPt = (1./pt)*1.1;
    return 1/newInvPt;
  }
}

TVector3 Analyzer::getOuterHitPos(const reco::DeDxHitInfo* dedxHits){
     TVector3 point(0,0,0);
     if(!dedxHits)return point;
     double outerDistance=-1;
     for(unsigned int h=0;h<dedxHits->size();h++){
        DetId detid(dedxHits->detId(h));  
        moduleGeom* geomDet = moduleGeom::get(detid.rawId());
        TVector3 hitPos = geomDet->toGlobal(TVector3(dedxHits->pos(h).x(), dedxHits->pos(h).y(), dedxHits->pos(h).z())); 
        if(hitPos.Mag()>outerDistance){outerDistance=hitPos.Mag();  point=hitPos;}
     }
     return point;
}

double Analyzer::SegSep(const susybsm::HSCParticle& hscp, const edm::Event& iEvent, double& minPhi, double& minEta){
  if(TypeMode!=3)return -1;

  reco::MuonRef muon = hscp.muonRef();
  if(muon.isNull()) return false;
  reco::TrackRef  track = muon->standAloneMuon();
  if(track.isNull())return false;

  /*edm::Handle<MuonSegmentCollection> SegCollHandle;
  SegCollHandle.getByLabel(ev, "MuonSegmentProducer");
  if(!SegCollHandle.isValid()){printf("Segment Collection Not Found\n"); return -1;}
  MuonSegmentCollection SegCollection = *SegCollHandle;*/

  susybsm::MuonSegmentCollection SegCollection = iEvent.get(muonSegmentToken_);

  double minDr=10;
  minPhi=10;
  minEta=10;

  //Look for segment on opposite side of detector from track
  for (susybsm::MuonSegmentCollection::const_iterator segment = SegCollection.begin(); segment!=SegCollection.end();++segment) {  
    GlobalPoint gp = segment->getGP();

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
   const edm::Event& iEvent, 
   Tuple* &tuple, 
   const double& GenBeta, 
   bool RescaleP, 
   const double& RescaleI, 
   const double& RescaleT, 
   double MassErr)
{

   if(TypeMode==1 && !(hscp.type() == susybsm::HSCParticleType::trackerMuon || hscp.type() == susybsm::HSCParticleType::globalMuon))return false;
   if( (TypeMode==2 || TypeMode==4) && hscp.type() != susybsm::HSCParticleType::globalMuon)return false;

   reco::TrackRef   track;
   reco::MuonRef muon = hscp.muonRef();

   if(TypeMode!=3) track = hscp.trackRef();
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
   if(TypeMode==3 && count<minMuStations) return false;
   if(tuple) tuple->Stations->Fill(0.0, Event_Weight);

   //===================== Handle For vertex ================
   //Handle<vector<reco::Vertex>> vertexCollH;
   //iEvent.getByToken(offlinePrimaryVerticesToken_,vertexCollH);

   vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken_);
   if(vertexColl.size()<1){printf("NO VERTEX\n"); return false;}

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

   if(TypeMode!=3 && track->found()<GlobalMinNOH)return false;

   if(TypeMode!=3 && track->hitPattern().numberOfValidPixelHits()<GlobalMinNOPH)return false;
   if(TypeMode!=3 && track->validFraction()<GlobalMinFOVH)return false;

   unsigned int missingHitsTillLast = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);;
   double validFractionTillLast = track->found()<=0?-1:track->found() / float(track->found() + missingHitsTillLast);
  
   if(tuple){
      tuple->BS_TNOHFractionTillLast->Fill(validFractionTillLast,Event_Weight);
	   tuple->BS_TNOMHTillLast->Fill(missingHitsTillLast,Event_Weight);
   }

   if(TypeMode!=3 && missingHitsTillLast>GlobalMaxNOMHTillLast)return false;
   if(TypeMode!=3 && validFractionTillLast<GlobalMinFOVHTillLast)return false;

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

   if(tuple){tuple->BS_nDof->Fill(tof.nDof(),Event_Weight);}
   if((TypeMode>1  && TypeMode!=5) && tof.nDof()<GlobalMinNDOF && (dttof.nDof()<GlobalMinNDOFDT || csctof.nDof()<GlobalMinNDOFCSC) )return false;

   if(tuple){
      tuple->nDof  ->Fill(0.0,Event_Weight);
      tuple->BS_Qual->Fill(track->qualityMask(),Event_Weight);
   }

   if(TypeMode!=3 && track->qualityMask()<GlobalMinQual )return false; // FIXME Tracks with quality > 2 are bad also!
//   if(TypeMode!=3 && track->qualityMask() != FixedQual)return false; // FIXME if this is true, no tracks pass eventually ... so what now?
   if(tuple){tuple->Qual  ->Fill(0.0,Event_Weight);
          tuple->BS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight);
   }
   if(TypeMode!=3 && track->chi2()/track->ndof()>GlobalMaxChi2 )return false;
   if(tuple){tuple->Chi2  ->Fill(0.0,Event_Weight);}

   if(tuple && GenBeta>=0)tuple->Beta_PreselectedA->Fill(GenBeta, Event_Weight);

   if(tuple){tuple->BS_MPt ->Fill(track->pt(),Event_Weight);}
   if(RescaleP){ if(RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())<GlobalMinPt)return false;
   }else{        if(track->pt()<GlobalMinPt)return false;   }

   if(tuple){tuple->MPt   ->Fill(0.0,Event_Weight);
     if(dedxSObj) tuple->BS_MIs->Fill(dedxSObj->dEdx(),Event_Weight);
     if(dedxMObj) tuple->BS_MIm->Fill(dedxMObj->dEdx(),Event_Weight);
   }

   if(dedxSObj && dedxSObj->dEdx()+RescaleI<GlobalMinIs)return false;
   if(dedxMObj && ((TypeMode!=5 && dedxMObj->dEdx()<GlobalMinIm) || (TypeMode==5 && dedxMObj->dEdx()>GlobalMinIm)) )return false;
   if(tuple){tuple->MI   ->Fill(0.0,Event_Weight);}

   if(tuple){tuple->BS_MTOF ->Fill(tof.inverseBeta(),Event_Weight);}
   //This cut is no longer applied here but rather in the PassSelection part to use the region
   //with TOF<GlobalMinTOF as a background check
   //if(TypeMode>1 && tof.inverseBeta()+RescaleT<GlobalMinTOF)return false;

   if(tuple)tuple->BS_TOFError->Fill(tof.inverseBetaErr(),Event_Weight);
   if((TypeMode>1  && TypeMode!=5) && tof.inverseBetaErr()>GlobalMaxTOFErr)return false;

   if(tuple) tuple->BS_TimeAtIP->Fill(tof.timeAtIpInOut(),Event_Weight);
   if(TypeMode==3 && min(min(fabs(tof.timeAtIpInOut()-100), fabs(tof.timeAtIpInOut()-50)), min(fabs(tof.timeAtIpInOut()+100), fabs(tof.timeAtIpInOut()+50)))<5) return false;

   if(tuple) tuple->BS_dzMinv3d->Fill(dz,Event_Weight);
   if(tuple) tuple->BS_dxyMinv3d->Fill(dxy,Event_Weight);
   if(tuple) tuple->BS_PV->Fill(goodVerts,Event_Weight);   
   if(tuple) tuple->BS_PV_NoEventWeight->Fill(goodVerts);
   if(tuple && dedxSObj) tuple->BS_NOMoNOHvsPV->Fill(goodVerts,dedxSObj->numberOfMeasurements()/(double)track->found(),Event_Weight);

   //Require at least one good vertex except if cosmic event
   //WAIT//if(TypeMode==3 && goodVerts<1 && (!tuple || tuple->Name.find("Cosmic")==string::npos)) return false;

   //For TOF only analysis match to a SA track without vertex constraint for IP cuts
   if(TypeMode==3) {

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
      /*edm::Handle<reco::BeamSpot> beamSpotCollHandle;
      beamSpotCollHandle.getByLabel(iEvent,"offlineBeamSpot");
      if(!beamSpotCollHandle.isValid()){printf("Beam Spot Collection NotFound\n");return false;}
      const reco::BeamSpot& beamSpotColl = *beamSpotCollHandle;*/

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
   if(TypeMode!=5 && fabs(dxy)>GlobalMaxDXY)return false;
   if(TypeMode==5 && fabs(dxy)>4)return false;
   if(TypeMode==5 && fabs(dxy)>GlobalMaxDXY) DXYSB = true;

   if(tuple){tuple->Dxy  ->Fill(0.0,Event_Weight);}

   if(TypeMode!=3) {
     /*edm::Handle<HSCPIsolationValueMap> IsolationH;
     IsolationH.getByLabel(ev, "HSCPIsolation", "R03"); //New format used for data since 17-07-2015
     if(!IsolationH.isValid()){
        IsolationH.getByLabel(ev, "HSCPIsolation03");//Old format used for first 2015B data, Signal and MC Backgrounds
        if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
     }
     const ValueMap<HSCPIsolation>& IsolationMap = *IsolationH.product();*/

     const edm::ValueMap<susybsm::HSCPIsolation> IsolationMap = iEvent.get(hscpIsoToken_);

     susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
     if(tuple){tuple->BS_TIsol ->Fill(hscpIso.Get_TK_SumEt(),Event_Weight);}
//     if(TypeMode!=4){       if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;     }
      if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;
     if(tuple){tuple->TIsol   ->Fill(0.0,Event_Weight);}

     double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
     if(tuple){tuple->BS_EIsol ->Fill(EoP,Event_Weight);}
//     if(TypeMode!=4){       if(EoP>GlobalMaxEIsol)return false;     }
     if(EoP>GlobalMaxEIsol)return false;
     if(tuple){tuple->EIsol   ->Fill(0.0,Event_Weight);}
     
     // relative tracker isolation
     if (tuple) {  tuple->BS_SumpTOverpT->Fill(hscpIso.Get_TK_SumEt()/track->pt(), Event_Weight); }
//     if(TypeMode==4) { if(hscpIso.Get_TK_SumEt()/track->pt()>GlobalMaxRelTIsol)return false;   }
     if(hscpIso.Get_TK_SumEt()/track->pt()>GlobalMaxRelTIsol)return false;
     if (tuple) {  tuple->SumpTOverpT   ->Fill(0.0,Event_Weight);} 
   }

   if(tuple){tuple->BS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);}
   if(TypeMode!=3 && (track->ptError()/track->pt())>GlobalMaxPterr)return false;
   //mk if(MassErr > 0 && MassErr > 2.2)return false; //FIXME jozze -- cut on relative mass error in units of 8*MassErr/Mass

   if(std::max(0.0,track->pt())<GlobalMinPt)return false;
   if(tuple){tuple->Pterr   ->Fill(0.0,Event_Weight);}

   //Find distance to nearest segment on opposite side of detector
   double minPhi, minEta;
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
   //if(TypeMode==3 && fabs(minEta)<minSegEtaSep) return false;
   if(tuple){tuple->SegSep->Fill(0.0,Event_Weight);}

   if(tuple) {
     //Plots for tracks in dz control region
     if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz && !muon->isGlobalMuon()) {
       tuple->BS_Pt_FailDz->Fill(track->pt(), Event_Weight);
       tuple->BS_TOF_FailDz->Fill(tof.inverseBeta(), Event_Weight);
       if(fabs(track->eta())>CSCRegion) {
	 tuple->BS_TOF_FailDz_CSC->Fill(tof.inverseBeta(), Event_Weight);
	 tuple->BS_Pt_FailDz_CSC->Fill(track->pt(), Event_Weight);
       }
       else if(fabs(track->eta())<DTRegion) {
	 tuple->BS_TOF_FailDz_DT->Fill(tof.inverseBeta(), Event_Weight);
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
   if(TypeMode==3 && !muon->isGlobalMuon() && tuple) {
     int DzType=-1;
     if(fabs(dz)<GlobalMaxDZ) DzType=0;
     else if(fabs(dz)<30) DzType=1;
     else if(fabs(dz)<50) DzType=2;
     else if(fabs(dz)<70) DzType=3;
     if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz) DzType=4;
     if(fabs(dz)>CosmicMaxDz) DzType=5;

     //Count number of tracks in dz sidebands passing the TOF cut
     //The pt cut is not applied to increase statistics
     for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
       if(tof.inverseBeta()>=CutTOF[CutIndex]) {
	 tuple->H_D_DzSidebands->Fill(CutIndex, DzType);
       }
     }
   }

   TreeDZ = dz;
   bool DZSB = false;
   if(TypeMode!=5 && fabs(dz)>GlobalMaxDZ) return false;
   if(TypeMode==5 && fabs(dz)>4) return false;
   if(TypeMode==5 && fabs(dz)>GlobalMaxDZ) DZSB = true;
   if(tuple){tuple->Dz  ->Fill(0.0,Event_Weight);}

   if(TypeMode==3 && fabs(minEta)<minSegEtaSep) return false;
   if(tuple)tuple->BS_Phi->Fill(track->phi(),Event_Weight);
   if(TypeMode==3 && fabs(track->phi())>1.2 && fabs(track->phi())<1.9) return false;

    //skip HSCP that are compatible with cosmics.
    if(tuple)tuple->BS_OpenAngle->Fill(OpenAngle,Event_Weight);

    bool OASB = false;
    if(TypeMode==5 && OpenAngle>=2.8)OASB = true;

   isCosmicSB = DXYSB && DZSB && OASB;
   isSemiCosmicSB = (!isCosmicSB && (DXYSB || DZSB || OASB));
 
   if(tuple){if(dedxSObj) tuple->BS_EtaIs->Fill(track->eta(),dedxSObj->dEdx(),Event_Weight);
          if(dedxMObj) tuple->BS_EtaIm->Fill(track->eta(),dedxMObj->dEdx(),Event_Weight);
          tuple->BS_EtaP ->Fill(track->eta(),track->p(),Event_Weight);
          tuple->BS_EtaPt->Fill(track->eta(),track->pt(),Event_Weight);
          tuple->BS_EtaTOF->Fill(track->eta(),tof.inverseBeta(),Event_Weight);
   }

   if(tuple){if(GenBeta>=0)tuple->Beta_PreselectedC->Fill(GenBeta, Event_Weight);
          if(DZSB  && OASB)tuple->BS_Dxy_Cosmic->Fill(dxy, Event_Weight);
          if(DXYSB && OASB)tuple->BS_Dz_Cosmic->Fill(dz, Event_Weight);
          if(DXYSB && DZSB)tuple->BS_OpenAngle_Cosmic->Fill(OpenAngle,Event_Weight);


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
          if(dedxSObj) tuple->BS_Im ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(dedxSObj && PUA) tuple->BS_Im_PUA ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(dedxSObj && PUB) tuple->BS_Im_PUB ->Fill(dedxMObj->dEdx(),Event_Weight);

	    tuple->BS_TOF->Fill(tof.inverseBeta(),Event_Weight);
            if(PUA)tuple->BS_TOF_PUA->Fill(tof.inverseBeta(),Event_Weight);
            if(PUB)tuple->BS_TOF_PUB->Fill(tof.inverseBeta(),Event_Weight);
	    if(dttof.nDof()>6) tuple->BS_TOF_DT->Fill(dttof.inverseBeta(),Event_Weight);
            if(csctof.nDof()>6) tuple->BS_TOF_CSC->Fill(csctof.inverseBeta(),Event_Weight);
            tuple->BS_PtTOF->Fill(track->pt() ,tof.inverseBeta(),Event_Weight);

          if(dedxSObj) {
	    tuple->BS_PIs  ->Fill(track->p()  ,dedxSObj->dEdx(),Event_Weight);
            tuple->BS_PImHD->Fill(track->p()  ,dedxMObj->dEdx(),Event_Weight);
            tuple->BS_PIm  ->Fill(track->p()  ,dedxMObj->dEdx(),Event_Weight);
            tuple->BS_PtIs ->Fill(track->pt() ,dedxSObj->dEdx(),Event_Weight);
            tuple->BS_PtIm ->Fill(track->pt() ,dedxMObj->dEdx(),Event_Weight);
	  }
          if(dedxSObj)tuple->BS_TOFIs->Fill(tof.inverseBeta(),dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj)tuple->BS_TOFIm->Fill(tof.inverseBeta(),dedxMObj->dEdx(),Event_Weight);

	  //Muon only prediction binned depending on where in the detector the track is and how many muon stations it has
	  //Binning not used for other analyses
	  int bin=-1;
	  if(TypeMode==3) {
	    if(fabs(track->eta())<DTRegion) bin=muonStations(track->hitPattern())-2;
	    else bin=muonStations(track->hitPattern())+1;
	    tuple->BS_Pt_Binned[to_string(bin)] ->Fill(track->pt(),Event_Weight);
	    tuple->BS_TOF_Binned[to_string(bin)]->Fill(tof.inverseBeta(),Event_Weight);
	  }
   }
   if(tuple){tuple->Basic  ->Fill(0.0,Event_Weight);}

   return true;
}