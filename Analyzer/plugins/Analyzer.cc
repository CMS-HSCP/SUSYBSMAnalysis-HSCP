// -*- C++ -*-
//
//
// Class : AnalyzerCustom
//
//
// Inspired by class Analyzer SUSYBSMAnalysis/Analyzer/plugins/Analyzer.cc
//
// Author : Raphael Haeberle
// 	 Created : Wed, 23 February 2022 16:12:48 GMT+1
//
//




#include "SUSYBSMAnalysis/Analyzer/plugins/Analyzer.h"



Analyzer::Analyzer(const edm::ParameterSet& iConfig)
    : hscpToken_(consumes<vector<susybsm::HSCParticle>>(iConfig.getParameter<edm::InputTag>("hscpCollection"))),
      hscpIsoToken_(consumes<edm::ValueMap<susybsm::HSCPIsolation>>(iConfig.getParameter<edm::InputTag>("hscpIsoCollection"))),
      l1resToken_(consumes<BXVector<GlobalAlgBlk>>(iConfig.getParameter<edm::InputTag>("l1results"))),
      dedxToken_(consumes<reco::DeDxHitInfoAss>(iConfig.getParameter<edm::InputTag>("dedxCollection"))),
      offlinePrimaryVerticesToken_(consumes<vector<reco::Vertex>>(iConfig.getParameter<edm::InputTag>("offlinePrimaryVerticesCollection"))),
      triggerResultsToken_(consumes<edm::TriggerResults>(iConfig.getParameter<edm::InputTag>("triggerResults"))),
      pfMETToken_(consumes<std::vector<reco::PFMET>>(iConfig.getParameter<edm::InputTag>("pfMET"))),
      pfJetToken_(consumes<reco::PFJetCollection>(iConfig.getParameter<edm::InputTag>("pfJet"))),
      CaloMETToken_(consumes<std::vector<reco::CaloMET>>(iConfig.getParameter<edm::InputTag>("CaloMET"))),
      genParticleToken_(
          consumes<std::vector<reco::GenParticle>>(iConfig.getParameter<edm::InputTag>("genParticleCollection"))),
      pfCandToken_(consumes<reco::PFCandidateCollection>(iConfig.getParameter<edm::InputTag>("PfCand"))),
      // HLT triggers
      trigger_met_(iConfig.getUntrackedParameter<vector<string>>("Trigger_MET")),    
      trigger_mu_(iConfig.getUntrackedParameter<vector<string>>("Trigger_Mu")),
      // =========Analysis parameters===============
      TypeMode_(iConfig.getUntrackedParameter<unsigned int>("TypeMode")),
      SampleType_(iConfig.getUntrackedParameter<unsigned int>("SampleType")),
      SampleName_(iConfig.getUntrackedParameter<string>("SampleName")),
      Period_(iConfig.getUntrackedParameter<string>("Period")),      
      DeDxTemplate(iConfig.getUntrackedParameter<string>("DeDxTemplate")),
      enableDeDxCalibration(iConfig.getUntrackedParameter<bool>("enableDeDxCalibration")),
      DeDxCalibration(iConfig.getUntrackedParameter<string>("DeDxCalibration")),
      Geometry(iConfig.getUntrackedParameter<string>("Geometry")),
      trackProbQCut_(iConfig.getUntrackedParameter<double>("trackProbQCut")),
      debugLevel_(iConfig.getUntrackedParameter<unsigned int>("debugLevel")),
      hltPSProv_(iConfig,consumesCollector(),*this),
      hltProcess_(iConfig.getParameter<std::string>("hltProcess")),
      scenarios_(iConfig.getUntrackedParameter<std::string>("Scenario")) 
 {
  //now do what ever initialization is needed
  // define the selection to be considered later for the optimization
  // WARNING: recall that this has a huge impact on the analysis time AND on the output file size --> be carefull with your choice


  useClusterCleaning = true;
  if (TypeMode_ == 4) {
    useClusterCleaning = false;  //switch off cluster cleaning for mCHAMPs
  }

  isData = (SampleType_ == 0);
  isBckg = (SampleType_ == 1);
  isSignal = (SampleType_ >= 2);

  EffL1Seeds.resize(nbl1names,0.0);
  L1Num.resize(nbl1names,0);
  L1Denom.resize(nbl1names,0);
  L1NumFail.resize(nbl1names,0);
  L1DenomFail.resize(nbl1names,0);
  L1Dec.resize(nbl1names,false);
  L1NoDec.resize(nbl1names,false);
  bool splitByModuleType = true;
  dEdxTemplates = loadDeDxTemplate(DeDxTemplate, splitByModuleType);
  if (enableDeDxCalibration)
    trackerCorrector.LoadDeDxCalibration(DeDxCalibration);
  else
    trackerCorrector.TrackerGains = nullptr;

  //moduleGeom::loadGeometry(Geometry);
  //tofCalculator.loadTimeOffset(TimeOffset);
}
      
Analyzer::~Analyzer() {
  // do anything here that needs to be done at desctruction time
  //   // (e.g. close files, deallocate resources etc.
}



// ------------ method called once each job just before starting event loop  ------------
void Analyzer::beginJob() {
  // Book histograms
  edm::Service<TFileService> fs;
  //tuple = new Tuple();

  TFileDirectory dir = fs->mkdir(SampleName_.c_str(), SampleName_.c_str());
  // create histograms & trees
  initializeCuts(fs, CutPt_, CutI_, CutTOF_, CutPt_Flip_, CutI_Flip_, CutTOF_Flip_);

  // Re-weighting

  tof = nullptr;
  dttof = nullptr;
  csctof = nullptr;

  TrigInfo_ = 0;

  CurrentRun_ = 0;

  is2016 = false;
  is2016G = false;
  

  /*
 *move out for migration part 2
 *
 *
 *
 *
 *
 *
 *
   */
}


void Analyzer::beginRun(const edm::Run& run,const edm::EventSetup& setup)
{
  bool changed=false;
  hltPSProv_.init(run,setup,hltProcess_,changed);
  const l1t::L1TGlobalUtil& l1GtUtils = hltPSProv_.l1tGlobalUtil();
  std::cout <<"l1 menu "<<l1GtUtils.gtTriggerMenuName()<<" version "<<l1GtUtils.gtTriggerMenuVersion()<<" comment "<<std::endl;
  std::cout << "Scenario for study : " << scenarios_ <<std::endl;
  std::cout <<"hlt name "<<hltPSProv_.hltConfigProvider().tableName()<<std::endl;
}


void Analyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup) {
  static constexpr const char* const MOD = "Analyzer";
  using namespace edm;

  vector<reco::GenParticle> genColl;
  double HSCPGenBeta1 = -1, HSCPGenBeta2 = -1;
  double HSCPDLength1 = -1, HSCPDLength2 = -1;
  int psColumn = hltPSProv_.prescaleSet(iEvent,iSetup);
  //std::cout <<"PS column "<<psColumn<<std::endl;
  if(psColumn==0 && iEvent.isRealData()){
    std::cout <<"PS column zero detected for data, this is unlikely (almost all triggers are disabled in normal menus here) and its more likely that you've not loaded the correct global tag in "<<std::endl;
  }
 
  l1t::L1TGlobalUtil& l1GtUtils = const_cast<l1t::L1TGlobalUtil&> (hltPSProv_.l1tGlobalUtil());
  //std::cout <<"l1 menu: name decisions prescale "<<std::endl;

  int nbl1=0;
  for(int j=0;j<nbl1names;j++){
    for(size_t bitNr=0;bitNr<l1GtUtils.decisionsFinal().size();bitNr++){
      const std::string& bitName = l1GtUtils.decisionsFinal()[bitNr].first;
      bool passInitial = l1GtUtils.decisionsInitial()[bitNr].second;
      bool passInterm = l1GtUtils.decisionsInterm()[bitNr].second;
      bool passFinal = l1GtUtils.decisionsFinal()[bitNr].second;
      int prescale = l1GtUtils.prescales()[bitNr].second;
      if(bitName == ListL1Names[j]){
        nbl1+=1;
        L1Dec[j] = passFinal;
      } 
    }
  }
  if(nbl1!=nbl1names){
    cout << "Didn't find every L1 seed for this event" << endl;
    return;
  }

  if (isSignal) {
    ntot +=1;
    Handle<vector<reco::GenParticle>> genCollH;
    iEvent.getByToken(genParticleToken_, genCollH);
    if (!genCollH.isValid()) {
      cout << "Warning, gen particle collection is invalid" << endl;
      LogWarning("Analyzer") << "Invalid  GenParticle!!, this event will be ignored";
      nbwrongcoll+=1;
      return;
    }

    genColl = *genCollH;

    int NChargedHSCP = HowManyChargedHSCP(genColl);
    int NNeutralHSCP = HowManyNeutralHSCP(genColl);
    int CheckNbHSCP = NChargedHSCP + NNeutralHSCP;

    if(CheckNbHSCP != 2 ){
      cout << " WARNING ----- Event does not contain 2 HSCPs, skipping it ----- " <<endl;
    }

    if(NChargedHSCP == 2){
        nchch+=1;
    }
    else if(NNeutralHSCP == 2){
        nnn+=1;

    }
    else if(NChargedHSCP == 1 && NNeutralHSCP == 1){
        nchn +=1;
    }

    //std::cout <<"# Charged : " << NChargedHSCP << " ,# Neutral : " << NNeutralHSCP << endl; 
    double Wa = 1.0, Wad = 1.0, Waa = 1.0,
           Wan = 1.0;  // Wa is additional weight for single other, Wad for other+double_charged,
                       // Waa for the event with 2 other R-hadron, Wan for other+neutral
                        
    bool Rhadron = 0;  // default value - not R-hadron (not need to weight)
    string sample_name = "";

    double calo_met;
    unsigned int nw = 0, na = 0, nd = 0;  //initialize counters: nw - wrong, na - other, nd - double charged, nn - neutral
    //===================== Handle For CaloMET =================
    const edm::Handle<std::vector<reco::CaloMET>> CaloMETHandle = iEvent.getHandle(CaloMETToken_);
    if (CaloMETHandle.isValid() && !CaloMETHandle->empty()) {
      //cout << "Size reco::calo MET collection : " << CaloMETHandle->size() << endl;
      for (unsigned int i = 0; i < CaloMETHandle->size(); i++) {
        const reco::CaloMET* calomet = &(*CaloMETHandle)[i];
        //cout << "calo MET[" << i << "] = " << calomet->et() << endl;
        calo_met = calomet->et();
      }
    }
    
    //===================== Handle For PFCandidate  =================
    const edm::Handle<reco::PFCandidateCollection> pfCandHandle = iEvent.getHandle(pfCandToken_);
    

   
    //int discr = 0;
    bool dcr = true,dcrno = true;
    //===================== Handle For DeDx Hits ==============

    Handle<reco::DeDxHitInfoAss> dedxCollH;
    iEvent.getByToken(dedxToken_, dedxCollH);

    //====================loop over HSCP candidates===================
    if (debugLevel_ > 0 ) LogPrint(MOD) << "Loop over HSCP candidates:";
    for (const auto& hscp : iEvent.get(hscpToken_)) {
      //if (debugLevel_> 0) LogPrint(MOD) << "  --------------------------------------------";
      reco::TrackRef track;
      track = hscp.trackRef();  
      if (track.isNull()) {
        //cout << "  >> Event has no track associated to this HSCP, skipping it";
        if (debugLevel_> 0) LogPrint(MOD) << "  >> Event has no track associated to this HSCP, skipping it";
        continue;
      }

      int ClosestGen;
      if (isSignal && DistToHSCP(hscp, genColl, ClosestGen, TypeMode_) > 0.03) {
        if (debugLevel_> 0) LogPrint(MOD) << "  >> Signal MC HSCP distance from gen to candidate is too big (" <<
      DistToHSCP(hscp, genColl, ClosestGen, TypeMode_) << "), skipping it";
      //cout << " Signal MC HSCP distance from gen to candidate is too big "<<endl;
      continue;
      }


      const reco::DeDxHitInfo* dedxHits = nullptr;
      if (TypeMode_ != 3 && !track.isNull()) {
        reco::DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());
        if (!dedxHitsRef.isNull())
          dedxHits = &(*dedxHitsRef);
      }

      int pdgId = 0;

      if (isSignal) {
        pdgId = genColl[ClosestGen].pdgId();

        //std::cout << "  >> GenId  " << pdgId << " >> gen status :" << status;
        //LogPrint(MOD) << "  >> GenId  " << pdgId << " >> gen status :" << status;
      }

      float RMin = 9999.;
      unsigned int idx_pf_RMin = 9999;
      bool pf_isMuon = false;



      if(pfCandHandle.isValid() && !pfCandHandle->empty()) {
        const reco::PFCandidateCollection* pf = pfCandHandle.product();
        for (unsigned int i = 0; i < pf->size(); i++){
          const reco::PFCandidate* pfCand = &(*pf)[i];
          float dr = deltaR(pfCand->eta(),pfCand->phi(),track->eta(),track->phi());
          //cout << dr; 
          if(dr < RMin){
            RMin = dr;
            idx_pf_RMin = i;
          }
        }
        //pf candidate ref to track 
        for(unsigned int i=0;i<pf->size();i++){
          const reco::PFCandidate* pfCand = &(*pf)[i];
          if(i == idx_pf_RMin) {
            pf_isMuon = pfCand->translatePdgIdToType(pfCand->pdgId()) == reco::PFCandidate::ParticleType::mu;
            //cout <<" pf_isMuon : " << pf_isMuon << endl;
            if(pf_isMuon)
              nbpfmuon+=1;
          }
        }
      }
      bool globlMuon = hscp.type() == susybsm::HSCParticleType::globalMuon;
      //cout <<"globalMuon : " << globlMuon << endl;
      if(globlMuon)
          nbglobalmuon+=1;
      bool TrkrMuon = hscp.type() == susybsm::HSCParticleType::trackerMuon;
      //cout <<"trackerMuon : " << TrkrMuon << endl;
      if(TrkrMuon)
          nbtrkmuon+=1;
      bool StdAloneMuon = hscp.type() == susybsm::HSCParticleType::standAloneMuon;
      //cout <<"StandAloneMatchedMuon : " << StdAloneMuon << endl;
      if(StdAloneMuon)
          nbstdalnmuon+=1;
      bool MtchStdAloneMuon = hscp.type() == susybsm::HSCParticleType::matchedStandAloneMuon;
      //cout <<"StandAloneMatchedMuon : " << MtchStdAloneMuon << endl;
      //frac muon or no muon : compare with old sample, or with sample from dylan 

      double dEdxErr = 0;
      reco::DeDxData dedxSObjTmp = computedEdx(dedxHits,
                                             dEdxSF,
                                             dEdxTemplates,
                                             true,
                                             useClusterCleaning,
                                             TypeMode_ == 5,
                                             false,
                                             trackerCorrector.TrackerGains,
                                             true,
                                             true,
                                             99,
                                             false,
                                             1,
                                             0.00,
                                             nullptr,
                                             0,
                                             pdgId,
                                             skipPixel,
                                             useTemplateLayer);
     
      reco::DeDxData dedxMObjTmp = computedEdx(dedxHits,
                                             dEdxSF,
                                             nullptr,
                                             true,
                                             useClusterCleaning,
                                             false,
                                             false,
                                             trackerCorrector.TrackerGains,
                                             true,
                                             true,
                                             99,
                                             false,
                                             1,
                                             0.15,
                                             nullptr,
                                             &dEdxErr,
                                             pdgId,
                                             skipPixel,
                                             useTemplateLayer);


      reco::DeDxData* dedxSObj = dedxSObjTmp.numberOfMeasurements() > 0 ? &dedxSObjTmp : nullptr;
      reco::DeDxData* dedxMObj = dedxMObjTmp.numberOfMeasurements() > 0 ? &dedxMObjTmp : nullptr;
      bool PRescale = false;
      nbtot+=1;
      /*
      Splitting scenarios r-hadrons 
      + infos L1 after HLT filters aswell
      ++ distrib cut variables, efficiencies step by step filters 
      */
      // comprendre les chooses differemment muon et MET"
      

      if(globlMuon){
        if(passPreselection(hscp, dedxHits, dedxSObj, dedxMObj, nullptr, iEvent, 1, nullptr, -1, PRescale, 0, 0, 0)){
          if(dcr){
          nbpasspresel+=1;
            if(scenarios_ == "both" ){
              AssoGenID(L1Dec,"pass");
            }      
            else if(scenarios_ == "chch" && NChargedHSCP == 2){
                AssoGenID(L1Dec,"pass");
            } 
            else if(scenarios_ == "chn" ){
              if(NChargedHSCP == 1 && NNeutralHSCP == 1 && calo_met>90 ){
                AssoGenID(L1Dec,"pass");
              }        
            }
          }
          dcr = false;
        }
        else{
          if(dcrno){
            nbnopasspresel+=1;
            if(scenarios_ == "both"){
              AssoGenID(L1Dec,"fail");
            }
            else if(scenarios_ == "chch" && NChargedHSCP == 2){
              AssoGenID(L1Dec,"fail");
            }
            else if(scenarios_ == "chn"){
              if(NChargedHSCP == 1 && NNeutralHSCP == 1 && calo_met>90 ){
                AssoGenID(L1Dec,"fail");
              }
            }
          }
          dcrno = false;
        }
      }      
    } // HSCP loop end
  } // end of IsSignal

}

void Analyzer::endJob() {

  //cout << "# HSCP passing preSelection : " << nbpasspresel << ", # total HSCPs : " << nbtot << endl;      cout << "--> " << (nbpasspresel*1.0/nbtot)*100 << " % " << endl;
  const string path = "MergeEff/PRESEL+MUON/";
  const string number = "1";
  const string passnm = "PassNum",passdnm = "PassDenom",failnm = "FailNum", faildnm = "FailDenom";
  const string exttxt = ".txt";
  const string passnumfile = path + passnm + number + exttxt,  passdenomfile = path +  passdnm + number + exttxt, failnumfile = path +  failnm + number + exttxt, faildenomfile = path +  faildnm + number + exttxt;
  ofstream numl1pass,denoml1pass,numl1fail,denoml1fail;
  numl1pass.open(passnumfile.c_str());
  denoml1pass.open(passdenomfile.c_str());
  numl1fail.open(failnumfile.c_str());
  denoml1fail.open(faildenomfile.c_str());
  
  
  cout << "Splitting scenarios on r-hadrons" << endl;
  cout << "--------------------------------" << endl;
  cout << "There was " << ntot << " - " << nbwrongcoll <<"(invalid Gen Col) = " <<ntot - nbwrongcoll << " events, and  " << nbtot  << " HSCPs"  <<endl;
  cout << "There was " << nbpasspresel << " unique (per event) HSCP passing presel, and " << nbnopasspresel << " who failed presel" << endl; 
  cout << nchch << " Charged-Charged, " << nchn << " Charged-Neutral, " << nnn << " Neutral-Neutral"<<endl;
  cout << "Charged-Charged proportion : " << (nchch*1.0/ntot)*100 << " %, Charged-Neutral : " << (nchn*1.0/ntot)*100 << " %, Neutral-Neutral : " << (nnn*1.0/ntot)*100 << " %"<<endl;

  cout << " Out of the " << nbtot << " HSCPS : " << endl;
  cout << "             " << nbglobalmuon << " were global muons" <<endl;
  cout << "             " << nbpfmuon << " were pf muons" <<endl;
  

  for (int i = 0; i < nbl1names; i++){
    if(L1Denom[i] == 0 || L1Num[i] == 0)
      EffL1Seeds.at(i) = 0.0;
    else{
      EffL1Seeds.at(i) = (L1Num[i]*1.0/L1Denom[i]);
    }
  }
 //here efficiencies of l1seed
  cout << "Efficiencies of L1 seeds after preselection, for scenario "<< scenarios_.c_str() <<endl;
  cout << "-------------------------------------------" << endl;

  cout << "HSCP PASSING PRESELECTION" << endl;

  for (int i = 0; i < nbl1names; i++){
    numl1pass << L1Num[i] << endl;
    denoml1pass << L1Denom[i] << endl;
    cout << ListL1Names[i] << " = " << L1Num[i] << "/" << L1Denom[i] << " = " <<EffL1Seeds.at(i)*100 << " ±" << sqrt((EffL1Seeds.at(i)*(1-EffL1Seeds.at(i)))/L1Denom[i])*100  << " %"<<endl;
  }

  cout << "HSCP FAILING PRESELECTION" << endl;
  for (int i = 0; i < nbl1names; i++){
    numl1fail << L1NumFail[i] << endl;
    denoml1fail << L1DenomFail[i] << endl;
    cout << ListL1Names[i] << " = " << L1NumFail[i] << "/" << L1DenomFail[i] << " = " << (L1NumFail[i]*1.0/L1DenomFail[i])*100 << " ±" <<  sqrt(( (L1NumFail[i]*1.0/L1DenomFail[i]) * ( 1- (L1NumFail[i]*1.0/L1DenomFail[i]) ) )/L1DenomFail[i])*100  << " %"<<endl;
  }


  cout << "There was " <<ntot << " valid events and " << nbwrongcoll << " wrong events (gen collection invalid)";
  numl1pass.close();
  denoml1pass.close();
  numl1fail.close();
  denoml1fail.close();
  EffL1Seeds.clear();
  
 //end eff L1
 
  
}



void Analyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

DEFINE_FWK_MODULE(Analyzer);



//=============================================================
//
//     Method for initializing pT and Is cuts
//
//=============================================================
void Analyzer::initializeCuts(edm::Service<TFileService>& fs,
                              vector<double>& CutPt,
                              vector<double>& CutI,
                              vector<double>& CutTOF,
                              vector<double>& CutPt_Flip,
                              vector<double>& CutI_Flip,
                              vector<double>& CutTOF_Flip) {
  CutPt.clear();
  CutI.clear();
  CutTOF.clear();
  CutPt_Flip.clear();
  CutI_Flip.clear();
  CutTOF_Flip.clear();

  CutPt.push_back(GlobalMinPt);
  CutI.push_back(GlobalMinIs);
  CutTOF.push_back(GlobalMinTOF);
  CutPt_Flip.push_back(GlobalMinPt);
  CutI_Flip.push_back(GlobalMinIs);
  CutTOF_Flip.push_back(GlobalMinTOF);

  if (TypeMode_ < 2) {
    for (double Pt = GlobalMinPt + 5; Pt < 200; Pt += 5) {
      for (double I = GlobalMinIs + 0.025; I < 0.45; I += 0.025) {
        CutPt.push_back(Pt);
        CutI.push_back(I);
        CutTOF.push_back(-1);
      }
    }
  } else if (TypeMode_ == 2) {
    for (double Pt = GlobalMinPt + 5; Pt < 120; Pt += 5) {
      if (Pt > 80 && ((int)Pt) % 10 != 0)
        continue;
      for (double I = GlobalMinIs + 0.025; I < 0.40; I += 0.025) {
        for (double TOF = GlobalMinTOF + 0.025; TOF < 1.35; TOF += 0.025) {
          CutPt.push_back(Pt);
          CutI.push_back(I);
          CutTOF.push_back(TOF);
        }
      }
    }
    for (double Pt = GlobalMinPt + 10; Pt < 90; Pt += 30) {
      for (double I = GlobalMinIs + 0.1; I < 0.30; I += 0.1) {
        for (double TOF = GlobalMinTOF - 0.05; TOF > 0.65; TOF -= 0.05) {
          CutPt_Flip.push_back(Pt);
          CutI_Flip.push_back(I);
          CutTOF_Flip.push_back(TOF);
        }
      }
    }
  } else if (TypeMode_ == 3) {
    for (double Pt = GlobalMinPt + 30; Pt < 450; Pt += 30) {
      for (double TOF = GlobalMinTOF + 0.025; TOF < 1.5; TOF += 0.025) {
        CutPt.push_back(Pt);
        CutI.push_back(-1);
        CutTOF.push_back(TOF);
      }
    }
    for (double Pt = GlobalMinPt + 30; Pt < 450; Pt += 60) {
      for (double TOF = GlobalMinTOF - 0.025; TOF > 0.5; TOF -= 0.025) {
        CutPt_Flip.push_back(Pt);
        CutI_Flip.push_back(-1);
        CutTOF_Flip.push_back(TOF);
      }
    }
  } else if (TypeMode_ == 4) {
    for (double I = GlobalMinIs + 0.025; I < 0.55; I += 0.025) {
      for (double TOF = GlobalMinTOF + 0.025; TOF < 1.46; TOF += 0.025) {
        CutPt.push_back(-1);
        CutI.push_back(I);
        CutTOF.push_back(TOF);
      }
    }
    for (double I = GlobalMinIs + 0.025; I < 0.55; I += 0.025) {
      for (double TOF = GlobalMinTOF - 0.025; TOF > 0.54; TOF -= 0.025) {
        CutPt_Flip.push_back(-1);
        CutI_Flip.push_back(I);
        CutTOF_Flip.push_back(TOF);
      }
    }
  } else if (TypeMode_ == 5) {
    for (double Pt = 75; Pt <= 150; Pt += 25) {
      for (double I = 0.0; I <= 0.45; I += 0.025) {
        CutPt.push_back(Pt);
        CutI.push_back(I);
        CutTOF.push_back(-1);
        CutPt_Flip.push_back(Pt);
        CutI_Flip.push_back(I);
        CutTOF_Flip.push_back(-1);
      }
    }
  }

  edm::LogInfo("Analyzer") << CutPt.size() << " Different Final Selection will be tested\n"
                           << CutPt_Flip.size()
                           << " Different Final Selection will be tested for background uncertainty";


  //Initialization of variables that are common to all samples

}





//=============================================================
//
//     Pre-Selection
//
//=============================================================

bool Analyzer::passPreselection(const susybsm::HSCParticle& hscp,
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
                                bool Ih_Iso_cut) {
  static constexpr const char* const MOD = "Analyzer";
  using namespace edm;
  // Tracker only analysis must have either a tracker muon or a global muon
  
  if (TypeMode_ == 1 &&
      !(hscp.type() == susybsm::HSCParticleType::trackerMuon || hscp.type() == susybsm::HSCParticleType::globalMuon)) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Tracker only analysis  w/o a tracker muon or a global muon";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for having a tracker muon or a global muon";
  }
  // Tracker + Muon analysis  must have either a global muon
  if ((TypeMode_ == 2 || TypeMode_ == 4) && hscp.type() != susybsm::HSCParticleType::globalMuon) {
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for having a global muon ";
  }

  // Define track reference and the muon reference
  reco::TrackRef track;
  reco::MuonRef muon = hscp.muonRef();
  // Assign track value, in case of TOF only track must come from the muons
  // otherwise take the track from the HSCP candidate track
  if (TypeMode_ != 3)
    track = hscp.trackRef();
  else {
    if (muon.isNull()) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: TOF only track w/o a muon";
      return false;
    } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for having a muon track";
    }
    track = muon->standAloneMuon();
  }

  // Check if the track collection exists

  if (track.isNull()) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: track collection is null";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for an existing track collecion";
  }
  
  // Check if eta is inside the max eta cut
  if (fabs(track->eta()) > GlobalMaxEta) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: eta is outside the max eta cut";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for eta cut";
  }

  // Cut on number of matched muon stations
  int count = muonStations(track->hitPattern());

  if (TypeMode_ == 3 && count < minMuStations)
    return false;

  //===================== Handle For vertex ================
  vector<reco::Vertex> vertexColl = iEvent.get(offlinePrimaryVerticesToken_);
  if (vertexColl.size() < 1) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: there is no vertex";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for vertex cut";
  }

  int highestPtGoodVertex = -1;
  int goodVerts = 0;
  double dzMin = 10000;
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
  double dz = track->dz(vertexColl[highestPtGoodVertex].position());
  double dxy = track->dxy(vertexColl[highestPtGoodVertex].position());
  bool PUA = (vertexColl.size() < 15);
  bool PUB = (vertexColl.size() >= 15);


  // Check the number of found hits
  if (TypeMode_ != 3 && track->found() < GlobalMinNOH) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Number of hits too low";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for number of hits cut";
  }

  // Check the number of pixel hits
  if (TypeMode_ != 3 && track->hitPattern().numberOfValidPixelHits() < GlobalMinNOPH) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Number of pixel hits too low";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for number of pixel hits cut ";
  }
  if (TypeMode_ != 3 && track->validFraction() < GlobalMinFOVH) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Valid hit fraction is too low";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for valid hit fraction";        
  }

  unsigned int missingHitsTillLast =
      track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) +
      track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
  double validFractionTillLast =
      track->found() <= 0 ? -1 : track->found() / float(track->found() + missingHitsTillLast);
  

  
  if (TypeMode_ != 3 && missingHitsTillLast > GlobalMaxNOMHTillLast) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Number of missing hits is too high";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for number of missing hits cut";
  }
  if (TypeMode_ != 3 && validFractionTillLast < GlobalMinFOVHTillLast) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Valid hit fraction w/ missing hits is too low";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for valid hit fraction w/ missing hits cut";
  }


  if (dedxSObj) {
    if (dedxSObj->numberOfMeasurements() < GlobalMinNOM) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Number of dEdx hits is too low";
      return false;
    } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for number of dEdx hits cut";
    }
  }

  // Select only high purity tracks 
  if (TypeMode_ != 3 && !track->quality(reco::TrackBase::highPurity)) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Not a high purity track";
    return false;
  }

  if (TypeMode_ != 3 && track->chi2() / track->ndof() > GlobalMaxChi2) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Chi2 / ndof is too high";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for chi2 / ndof cut";
  }

  if (track->pt() < GlobalMinPt)
    return false;
  
  if (dedxSObj) {
    if (dedxSObj->dEdx() + RescaleI < GlobalMinIs) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Rescaled Ih is too low for fractionally charged";
      return false;
    } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for rescaled Ih cut";
    }
  }
  if (dedxMObj) {
    if ((TypeMode_ != 5 && dedxMObj->dEdx() < GlobalMinIm) || (TypeMode_ == 5 && dedxMObj->dEdx() > GlobalMinIm)) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Ih is too low OR Ih is too high for fractionally charged";
      return false;
    } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for Ih cut ";
    }
  }

 

    //

    
    double v3d = sqrt(dz * dz + dxy * dxy);
    

    if (v3d > GlobalMaxV3D) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: 3D distance from vertex is too high";
      return false;
    } else {
      if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for 3D distance cut";
    }

  TreeDXY = dxy;
  //bool DXYSB = false;
  if (TypeMode_ != 5 && fabs(dxy) > GlobalMaxDXY) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: dxy is too high";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for dxy cut";
  }
  if (TypeMode_ == 5 && fabs(dxy) > 4) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: dxy is too high";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for dxy cut";
  }
  if (TypeMode_ == 5 && fabs(dxy) > GlobalMaxDXY) {
    //DXYSB = true;
  }



  if (TypeMode_ != 3 && Ih_Iso_cut) {
    const edm::ValueMap<susybsm::HSCPIsolation> IsolationMap = iEvent.get(hscpIsoToken_);

    susybsm::HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());

    if (hscpIso.Get_TK_SumEt() > GlobalMaxTIsol) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Tracker based isolation is too high";
      return false;
    } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for tracker based isolation cut";
    }


    double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy()) / track->p();

    if (EoP > GlobalMaxEIsol) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Calo based isolation is too high";
      return false;
    } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for calo based isolation cut";
    }

    // relative tracker isolation

    if (hscpIso.Get_TK_SumEt() / track->pt() > GlobalMaxRelTIsol) {
      if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: Tracker based relative isolation is too high";
      return false;
    } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for tracker based relative isolation cut";
    }

  }


  if (TypeMode_ != 3 && (track->ptError() / track->pt()) > GlobalMaxPterr) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: pt error is too high";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for pt error cut";
  }
  if (std::max(0.0, track->pt()) < GlobalMinPt) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: pt too is low";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for pt cut";
  }

  //Find distance to nearest segment on opposite side of detector
  double minPhi = 0.0, minEta = 0.0;


  // Missing split different dz region to predict cosmic background
  TreeDZ = dz;
  //bool DZSB = false;
  if (TypeMode_ != 5 && fabs(dz) > GlobalMaxDZ) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: dz it too high";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for dz cut";
  }
  if (TypeMode_ == 5 && fabs(dz) > 4) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: dz it too high";
    return false;
  } else {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for dz cut";
  }
  if (TypeMode_ == 5 && fabs(dz) > GlobalMaxDZ)
    //DZSB = true;


  if (TypeMode_ == 3 && fabs(minEta) < minSegEtaSep) {
    if (debugLevel_ > 4 ) LogPrint(MOD) << "        >> Preselection not passed: for TOF only analysis, eta is too low";
    return false;
  } else if (TypeMode_ == 3 && fabs(minEta) > minSegEtaSep) {
    if (debugLevel_ > 5 ) LogPrint(MOD) << "        >> Preselection criteria passed for TOF eta cut";
  }

  if (TypeMode_ == 3 && fabs(track->phi()) > 1.2 && fabs(track->phi()) < 1.9) {
    return false;
  }
  //skip HSCP that are compatible with cosmics.


  return true;
}


//=============================================================
//
//     Selection
//
//=============================================================
bool Analyzer::passSelection(const susybsm::HSCParticle& hscp,
                             const reco::DeDxData* dedxSObj,
                             const reco::DeDxData* dedxMObj,
                             const reco::MuonTimeExtra* tof,
                             const edm::Event& iEvent,
                             float Event_Weight,
                             const int& CutIndex,
                             Tuple*& tuple,
                             const bool isFlip,
                             const double& GenBeta,
                             bool RescaleP,
                             const double& RescaleI,
                             const double& RescaleT) {
  reco::TrackRef track;
  if (TypeMode_ != 3)
    track = hscp.trackRef();
  else {
    reco::MuonRef muon = hscp.muonRef();
    if (muon.isNull())
      return false;
    track = muon->standAloneMuon();
  }
  if (track.isNull())
    return false;


  double Is = 0;
  if (dedxSObj)
    Is = dedxSObj->dEdx();
/*
  double Ih = 0;
  if (dedxMObj)
    Ih = dedxMObj->dEdx();
*/
  double PtCut = CutPt_[CutIndex];
  double ICut = CutI_[CutIndex];
  //double TOFCut = CutTOF_[CutIndex];
  if (isFlip) {
    PtCut = CutPt_Flip_[CutIndex];
    ICut = CutI_Flip_[CutIndex];
    //TOFCut = CutTOF_Flip_[CutIndex];
  }

  if (track->pt() < PtCut)
    return false;
  

  if (TypeMode_ != 3 && Is + RescaleI < ICut)
    return false;

  return true;


}


void Analyzer::AssoGenID(const vector<bool>& decisionl1, const string& mode){
  
  if(mode == "pass"){
    for(int i =0; i< nbl1names; i++){
      L1Denom[i]+=1;
      if(decisionl1[i]==true){
        L1Num[i]+=1;
      }
    }
  }
  else if(mode == "fail"){
    for(int i =0; i< nbl1names; i++){
      L1DenomFail[i]+=1;
      if(decisionl1[i]==true){
        L1NumFail[i]+=1;
      }
    }
  }
}















