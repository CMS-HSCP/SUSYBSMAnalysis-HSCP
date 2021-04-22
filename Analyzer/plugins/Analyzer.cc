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

// Collaborating Class Header
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


Analyzer::Analyzer(const edm::ParameterSet& iConfig)
 :_hscpToken(consumes<vector<susybsm::HSCParticle>>(iConfig.getParameter<edm::InputTag>("hscpCollection")))
 //,_dedxToken(consumes<vector<reco::DeDxHitInfo>>(iConfig.getParameter<edm::InputTag>("dedxCollection")))
 //,_tofToken(consumes<vector<edm::ValueMap<reco::MuonTimeExtra>>>(iConfig.getParameter<edm::InputTag>("tofCollection")))
 ,_dedxToken(consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dedxCollection")))
 /*,_tofToken(consumes<edm::ValueMap<reco::MuonTimeExtra>>(iConfig.getParameter<edm::InputTag>("tofCollection")))
 ,_tofDtToken(consumes<edm::ValueMap<reco::MuonTimeExtra>>(iConfig.getParameter<edm::InputTag>("tofDtCollection")))
 ,_tofCscToken(consumes<edm::ValueMap<reco::MuonTimeExtra>>(iConfig.getParameter<edm::InputTag>("tofCscCollection")))*/
 ,muonTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonTimeCollection")))
 ,muonDtTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonDtTimeCollection")))
 ,muonCscTimeToken_(consumes<reco::MuonTimeExtraMap>(iConfig.getParameter<edm::InputTag>("muonCscTimeCollection")))
 ,muonDtSegmentToken_(consumes<DTRecSegment4DCollection>(iConfig.getParameter<edm::InputTag>("muonDtSegmentCollection")))
 ,muonCscSegmentToken_(consumes<CSCSegmentCollection>(iConfig.getParameter<edm::InputTag>("muonCscSegmentCollection")))
 ///consumes<edm::ValueMap<reco::MuonTimeExtra>>(iConfig.getParameter<edm::InputTag>("sourceMuonTimeExtra"));
 //,_tracksToken(consumes<reco::Track>(iConfig.getParameter<edm::InputTag>("tracks")))
 ,Debug(iConfig.getUntrackedParameter<bool>("Debug"))
 ,AddTree(iConfig.getUntrackedParameter<bool>("AddTree"))
 ,SampleType(iConfig.getUntrackedParameter<unsigned int>("SampleType"))
 ,TypeMode(iConfig.getUntrackedParameter<unsigned int>("TypeMode"))
 ,SampleTxtFile(iConfig.getUntrackedParameter<string>("SampleTxtFile"))
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

   CutPt.clear();       CutI.clear();       CutTOF.clear();      
   CutPt_Flip.clear();  CutI_Flip.clear();  CutTOF_Flip.clear();  
   
   CutPt     .push_back(GlobalMinPt);   CutI       .push_back(GlobalMinIs);  CutTOF     .push_back(GlobalMinTOF);
   CutPt_Flip.push_back(GlobalMinPt);   CutI_Flip  .push_back(GlobalMinIs);  CutTOF_Flip.push_back(GlobalMinTOF);

   if(TypeMode<2){   
      for(double Pt =GlobalMinPt+5 ; Pt <200;Pt+=5){
      for(double I  =GlobalMinIs+0.025  ; I  <0.45 ;I+=0.025){
         CutPt .push_back(Pt);   CutI  .push_back(I);  CutTOF.push_back(-1);
      }}
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

   //InitBaseDirectory(): directory containing the EDM files
   vector<stSample> samples;
   vector<stSample> samplesFull;
   samplesFull = samples;
   GetSampleDefinition(samples, SampleTxtFile);// input "Analysis_Samples.txt" and return "samples"
   if(true){ //if(MODE.find("ANALYSE_")==0){
      //int sampleIdStart, sampleIdEnd; 
      //string MODE="ANALYSE_0_to_1";
      //sscanf(MODE.c_str(),"ANALYSE_%d_to_%d",&sampleIdStart, &sampleIdEnd);
      //keepOnlyTheXtoYSamples(samples,sampleIdStart,sampleIdEnd);
      keepOnlyValidSamples(samples);
      printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n");
      printf("Run on the following samples:\n");
      for(unsigned int s=0;s<samples.size();s++){samples[s].print();}
      printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
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
   //_hists["tracks_charge"] = fs->make<TH1F>("tracks charge","charge",4,-2,2);

   InitHist(fs);
   
   if (AddTree) InitTree(Tree, fs);//// Create tree variables and branches

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
      trackerCorrector.setRun(CurrentRun);
   }

   //===================== Handle For DeDx Hits ==============
   Handle<reco::DeDxHitInfoAss> dedxCollH;
   iEvent.getByToken(_dedxToken,dedxCollH);
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
   for(const auto& hscp : iEvent.get(_hscpToken)){
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
      if(TypeMode==5)OpenAngle = deltaROpositeTrack(iEvent.get(_hscpToken), hscp); //OpenAngle is a global variable... that's uggly C++, but that's the best I found so far

      //compute systematic uncertainties on signal
      if(isSignal){
         //FIXME to be measured on 2015 data, currently assume 2012
         /*bool   PRescale = true;
         double IRescale =-0.05; // added to the Ias value
         double MRescale = 0.95;
		   double TRescale =-0.015; //-0.005 (used in 2012); // added to the 1/beta value*/
		  
		   //double genpT = -1.0;
		   for(unsigned int g=0;g<genColl.size();g++) {
            if(genColl[g].pt()<5)continue;
            if(genColl[g].status()!=1)continue;
            int AbsPdg=abs(genColl[g].pdgId());
            if(AbsPdg!=17)continue;
            
            double separation = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
            if (separation > 0.03) continue;
            //genpT = genColl[g].pt();
            break;
         }
         //WAIT//if (genpT>0) {  SamplePlots->genrecopT->Fill(genpT, track->pt()); }

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
      //WAIT//if(isMC)PassPreselection( hscp, dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev, MCTrPlots  , -1, false, 0, 0, GetMassErr (track->p(), track->ptError(), dedxMObj?dedxMObj->dEdx():-1, dEdxErr, GetMass(track->p(), dedxMObj?dedxMObj->dEdx():-1, !isData)));
      //WAIT//if(    !PassPreselection( hscp, dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev, SamplePlots, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1, false, 0, 0, GetMassErr (track->p(), track->ptError(), dedxMObj?dedxMObj->dEdx():-1, dEdxErr, GetMass(track->p(), dedxMObj?dedxMObj->dEdx():-1, !isData)))) continue;
      if(TypeMode==5 && isSemiCosmicSB)continue;

      //fill the ABCD histograms and a few other control plots
      //WAIT//if(isData)Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, SamplePlots);
      //WAIT//else if(isMC) Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, MCTrPlots);
      if(TypeMode==5 && isCosmicSB)continue; 

      //Find the number of tracks passing selection for TOF<1 that will be used to check the background prediction
      //cout << ">>> P=" << track->p() << ", I=" << dedxMObj->dEdx() << ", K=" << DeDxK << ", C=" << DeDxC << " ==> Mass=" << GetMass(track->p(),dedxMObj->dEdx(),DeDxK,DeDxC) << endl;
      
      if(isMC || isData) {
         //compute the mass of the candidate, for TOF mass flip the TOF over 1 to get the mass, so 0.8->1.2
		   double Mass = GetMass(track->p(),dedxMObj->dEdx(),DeDxK,DeDxC);
		   //double MassTOF  = -1; if(useTOF)MassTOF = GetTOFMass(track->p(),(2-1));//GetTOFMass(track->p(),(2-tof.inverseBeta()));
		   //double MassComb = -1;
		   /*if(useTOF && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),DeDxK,DeDxC) + (1/(2-tof.inverseBeta())))*0.5 ) ;
		   if(dedxMObj) MassComb = Mass;
		   if(useTOF)MassComb=GetMassFromBeta(track->p(),(1/(2-tof.inverseBeta())));*/

         for(unsigned int CutIndex=0;CutIndex<CutPt_Flip.size();CutIndex++){
            //Background check looking at region with TOF<1
            //WAIT//if(!PassSelection   (hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, true)) continue;

            //Fill Mass Histograms
            HMass["Mass_Flip"]->Fill(CutIndex, Mass,Event_Weight);
            /*if(tof){
               MassTOF_Flip->Fill(CutIndex, MassTOF, Event_Weight);
            }
            MassComb_Flip->Fill(CutIndex, MassComb, Event_Weight);*/
		   }

      }

   } //END loop over HSCP candidates



   /*for(const auto& track : iEvent.get(_tracksToken) ) {
      // do something with track parameters, e.g, plot the charge.
      int charge = track.charge();
      _hists["tracks_charge"]->Fill(charge);
   }*/

   if (AddTree) {//FillTree(Tree, iEvent);
      Tree_Trig = 99;
      Tree_Run   = iEvent.id().run();
      Tree_Event = iEvent.id().event();
      Tree_Lumi  = iEvent.luminosityBlock();
   }

   if (AddTree) Tree->Fill();

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
//     Selection
//
//=============================================================
bool Analyzer::PassPreselection(
   const susybsm::HSCParticle& hscp, 
   const reco::DeDxHitInfo* dedxHits,  
   const reco::DeDxData* dedxSObj, 
   const reco::DeDxData* dedxMObj,
   const edm::Event& iEvent, 
   TTree* &tree,//stPlots* st, 
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

   Total->Fill(0.0,Event_Weight);
   //if(GenBeta>=0)st->Beta_Matched->Fill(GenBeta, Event_Weight);
   //st->BS_Eta->Fill(track->eta(),Event_Weight);

   return true;
}

//=============================================================
//
//     Create tree variables and branches
//
//=============================================================
void Analyzer::InitTree(TTree * &tree, edm::Service<TFileService> &fs){
   tree = fs->make<TTree>("HscpCandidates", "HscpCandidates");//tree  = new TTree("HscpCandidates", "HscpCandidates");

   tree->Branch("Trig"       ,&Tree_Trig         ,"Trig/I");
   tree->Branch("Run"        ,&Tree_Run          ,"Run/I");
   tree->Branch("Event"      ,&Tree_Event        ,"Event/I");
   tree->Branch("Lumi"       ,&Tree_Lumi         ,"Lumi/I");
}

void Analyzer::FillTree(TTree* &tree, const edm::Event& iEvent){
   Tree_Trig = 99;
   Tree_Run   = iEvent.id().run();
   Tree_Event = iEvent.id().event();
   Tree_Lumi  = iEvent.luminosityBlock();

   tree->Fill();
}

void Analyzer::InitHist(edm::Service<TFileService> &fs){
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

   //Mass
   HMass["Mass_Flip"] = fs->make<TH2F>("Mass_Flip" ,"Mass_Flip" ,CutPt_Flip.size(),0,CutPt_Flip.size(), MassNBins, 0, MassHistoUpperBound);

   //std::string Name;
   //Name = "Total";    Total   = fs->make<TH1F>(Name.c_str(), Name.c_str(),  1    , 0,  1);
   //Name = "Beta_Matched"     ; Beta_Matched     = fs->make<TH1F>(Name.c_str(), Name.c_str(),                 20, 0,  1);  Beta_Matched     ->Sumw2();
   Total   = fs->make<TH1F>("Total", "Total",  1    , 0,  1);
   Beta_Matched     = fs->make<TH1F>("Beta_Matched", "Beta_Matched",                 20, 0,  1);  Beta_Matched     ->Sumw2();
}