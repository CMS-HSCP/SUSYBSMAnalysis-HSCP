#include <exception>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphAsymmErrors.h"
#include "TProfile.h"
#include "TPaveText.h"


namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     {class TriggerResults; class TriggerResultsByName; class InputTag;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"
#include "FWCore/Common/interface/TriggerResultsByName.h"

#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;


#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"

#endif


struct plotSt{
   TH1D* NVert;
   TH1D* PreselEff;
   TH1D* TPsVNEvts;
//   TH1D* PtBPS;
   TH1D* Pt;
   TH2D* PtEta;
//   TH1D* PtErr;
   TH1D* dEdxHitStrip;
   TH1D* dEdxHitPixel;
   TH1D* dEdxMin1;
   TH1D* dEdxMin2;
   TH1D* dEdxMin3;
   TH1D* dEdxMin4;
   TH1D* dEdx;
   TH1D* dEdxOld;
   TH1D* dEdxMT;
   TH1D* dEdxM;
   TH1D* dEdxMS;
   TH1D* dEdxMP;
   TH1D* dEdxMSC;
   TH1D* dEdxMPC;
   TH1D* dEdxMSF;
   TH1D* dEdxMPF;
   TH1D* TOF;
   TH1D* TOFDT;
   TH1D* TOFCSC;
   TH1D* Vertex;
   TH1D* VertexDT;
   TH1D* VertexCSC;
   std::map<unsigned int, TH1D** > dEdxHitPerLumi;

   plotSt(string prefix, string sufix){
      string histoName;              
      histoName=prefix + "PreselEff"    + sufix ; PreselEff    = new TH1D(histoName.c_str(), histoName.c_str(),    1, 0.0, 1.0);
      histoName=prefix + "TPsVNEvts"    + sufix ; TPsVNEvts    = new TH1D(histoName.c_str(), histoName.c_str(),    1, 0.0, 1.0);
      histoName=prefix + "NVert"        + sufix ; NVert        = new TH1D(histoName.c_str(), histoName.c_str(),  100, 0.0, 100);
//      histoName=prefix + "PtBPS"        + sufix ; PtBPS        = new TH1D(histoName.c_str(), histoName.c_str(), 1000, 0.0,1000);
      histoName=prefix + "Pt"           + sufix ; Pt           = new TH1D(histoName.c_str(), histoName.c_str(), 1000, 0.0,1000);
      double EtaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4};
      histoName=prefix + "PtEta"        + sufix ; PtEta        = new TH2D(histoName.c_str(), histoName.c_str(), sizeof(EtaBins)/sizeof(double)-1, EtaBins, 1000, 0.0, 1000);
//      histoName=prefix + "PtErr"        + sufix ; PtErr        = new TH1D(histoName.c_str(), histoName.c_str(), 1000, 0.0,32.0);
      histoName=prefix + "dEdxHitStrip" + sufix ; dEdxHitStrip = new TH1D(histoName.c_str(), histoName.c_str(),  400, 0.0,20.0);
      histoName=prefix + "dEdxHitPixel" + sufix ; dEdxHitPixel = new TH1D(histoName.c_str(), histoName.c_str(),  400, 0.0,20.0);
      histoName=prefix + "dEdxMin1"     + sufix ; dEdxMin1     = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMin2"     + sufix ; dEdxMin2     = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMin3"     + sufix ; dEdxMin3     = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMin4"     + sufix ; dEdxMin4     = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdx"         + sufix ; dEdx         = new TH1D(histoName.c_str(), histoName.c_str(),  100, 0.0, 1.0);
      histoName=prefix + "dEdxOld"      + sufix ; dEdxOld      = new TH1D(histoName.c_str(), histoName.c_str(),  100, 0.0, 1.0);
      histoName=prefix + "dEdxMT"       + sufix ; dEdxMT       = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxM"        + sufix ; dEdxM        = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMS"       + sufix ; dEdxMS       = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMP"       + sufix ; dEdxMP       = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMSC"      + sufix ; dEdxMSC      = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMPC"      + sufix ; dEdxMPC      = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMSF"      + sufix ; dEdxMSF      = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "dEdxMPF"      + sufix ; dEdxMPF      = new TH1D(histoName.c_str(), histoName.c_str(),  200, 0.0,10.0);
      histoName=prefix + "TOF"          + sufix ; TOF          = new TH1D(histoName.c_str(), histoName.c_str(),  100, -1.0, 3.0);
      histoName=prefix + "TOFDT"        + sufix ; TOFDT        = new TH1D(histoName.c_str(), histoName.c_str(),  100, -1.0, 3.0);
      histoName=prefix + "TOFCSC"       + sufix ; TOFCSC       = new TH1D(histoName.c_str(), histoName.c_str(),  100, -1.0, 3.0);
      histoName=prefix + "Vertex"       + sufix ; Vertex       = new TH1D(histoName.c_str(), histoName.c_str(),  100, -10.0, 10.0);
      histoName=prefix + "VertexDT"     + sufix ; VertexDT     = new TH1D(histoName.c_str(), histoName.c_str(),  100, -10.0, 10.0);
      histoName=prefix + "VertexCSC"    + sufix ; VertexCSC    = new TH1D(histoName.c_str(), histoName.c_str(),  100, -10.0, 10.0);
   };
};
std::map<string, std::map<string, plotSt*> > MapRunTriggerPlots;

bool PassPreselection(const susybsm::HSCParticle& hscp,  const reco::DeDxData& dedxSObj, const reco::DeDxData& dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev);
bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold=false);



bool PassPreselection(const susybsm::HSCParticle& hscp,  const reco::DeDxData& dedxSObj, const reco::DeDxData& dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev)
{
   if(TypeMode==1 && !(hscp.type() == HSCParticleType::trackerMuon || hscp.type() == HSCParticleType::globalMuon))return false;
   if(TypeMode==2 && hscp.type() != HSCParticleType::globalMuon)return false;
   reco::TrackRef   track = hscp.trackRef(); if(track.isNull())return false;

   if(fabs(track->eta())>GlobalMaxEta) return false;
   if(track->found()<GlobalMinNOH)return false;
   if(track->hitPattern().numberOfValidPixelHits()<2)return false; 
   if(dedxSObj.numberOfMeasurements()<GlobalMinNOM)return false;
//   if(tof && tof->nDof()<GlobalMinNDOF && (dttof->nDof()<GlobalMinNDOFDT || csctof->nDof()<GlobalMinNDOFCSC) )return false;

   if(track->qualityMask()<GlobalMinQual )return false;
   if(track->chi2()/track->ndof()>GlobalMaxChi2 )return false;
   if(track->pt()<GlobalMinPt)return false;
   if(dedxSObj.dEdx()<GlobalMinIs)return false;
   if(dedxMObj.dEdx()<GlobalMinIm)return false;
//   if(tof && tof->inverseBeta()<GlobalMinTOF)return false;
//   if(tof && tof->inverseBetaErr()>GlobalMaxTOFErr)return false;



   fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
   vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
   if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");return false;}
   const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
   if(vertexColl.size()<1){printf("NO VERTEX\n"); return false;}

   int highestPtGoodVertex = -1;
   int goodVerts=0;
   double dzMin=10000;
   for(unsigned int i=0;i<vertexColl.size();i++){
     if(vertexColl[i].isFake() || fabs(vertexColl[i].z())>24 || vertexColl[i].position().rho()>2 || vertexColl[i].ndof()<=4)continue; //only consider good vertex
     goodVerts++;

//     if(highestPtGoodVertex<0)highestPtGoodVertex = i;
     if(fabs(track->dz (vertexColl[i].position())) < fabs(dzMin) ){
         dzMin = fabs(track->dz (vertexColl[i].position()));
         highestPtGoodVertex = i;
//       dz  = track->dz (vertexColl[i].position());
//       dxy = track->dxy(vertexColl[i].position());
     }
   }if(highestPtGoodVertex<0)highestPtGoodVertex=0;

   double dz  = track->dz (vertexColl[highestPtGoodVertex].position());
   double dxy = track->dxy(vertexColl[highestPtGoodVertex].position());

   double v3d = sqrt(dz*dz+dxy*dxy);
   if(v3d>GlobalMaxV3D )return false;

   fwlite::Handle<HSCPIsolationValueMap> IsolationH;
   IsolationH.getByLabel(ev, "HSCPIsolation", "R03"); //New format used for data since 17-07-2015
   if(!IsolationH.isValid()){
     IsolationH.getByLabel(ev, "HSCPIsolation03");//Old format used for first 2015B data, Signal and MC Backgrounds
     if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
   }
  const ValueMap<HSCPIsolation>& IsolationMap = *IsolationH.product();


   HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
    if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;

   double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
   if(EoP>GlobalMaxEIsol)return false;

   if((track->ptError()/track->pt())>GlobalMaxPterr)return false;
   if(std::max(0.0,track->pt() - track->ptError())<GlobalMinPt)return false;
   return true;
}


bool PassingTrigger(const fwlite::ChainEvent& ev, const std::string& TriggerName){
   edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
   if(!tr.isValid())         tr = ev.triggerResultsByName("MergeHLT");
   if(!tr.isValid())return false;

   if(TriggerName=="Any"){
      if(passTriggerPatterns(tr, "HLT_PFMET170_NoiseCleaned_v*"))return true;
      if(passTriggerPatterns(tr, "HLT_Mu45_eta2p1_v*"))return true;
      if(passTriggerPatterns(tr, "HLT_Mu50_v*"))return true;
   }else{
      if(passTriggerPatterns(tr, (TriggerName + "_v*").c_str()))return true;
   }
   return false;
}

void StabilityCheck(string DIRNAME="COMPILE", string OUTDIRNAME="pictures", string JobIndexStr="0", string NJobsStr="1")
{
  printf("DIRNAME = %s\n", DIRNAME.c_str());
  if(DIRNAME=="COMPILE") return;
  OUTDIRNAME+="/";

   std::vector<string> triggers;
   triggers.push_back("Any");
//   triggers.push_back("HLT_Mu45_eta2p1");
//   triggers.push_back("HLT_Mu50");
//   triggers.push_back("HLT_PFMET170_NoiseCleaned");
//   triggers.push_back("HLT_PFMET170_HBECleaned");

   std::vector<string> versions;
   versions.push_back("");
//   versions.push_back("AOD");
//   versions.push_back("FAKE");




  int JobIndex;  sscanf(JobIndexStr.c_str(),"%d",&JobIndex);
  int NJobs;     sscanf(NJobsStr   .c_str(),"%d",&NJobs);
   char OutputFileName[1024];  sprintf(OutputFileName, "%s/Histos_%i.root", OUTDIRNAME.c_str(), JobIndex);
   TFile* OutputHisto = new TFile(OutputFileName,"RECREATE");
   TypeMode      = 0;


   Event_Weight = 1;
   MaxEntry = -1;

   system((string("mkdir -p ") + OUTDIRNAME).c_str());

   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.07);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505,"X");
   TH1::AddDirectory(kTRUE);

   std::map<unsigned int, unsigned int> RunBinIndex;
   unsigned int NextIndex=0;

   InitBaseDirectory();
   GetSampleDefinition(samples , DIRNAME+"/../../AnalysisCode/Analysis_Samples.txt");
//   GetSampleDefinition(samples , DIRNAME+"/Analysis_Samples_tmp.txt");
   int sampleIdStart, sampleIdEnd; sscanf(JobIndexStr.c_str(),"%d",&sampleIdStart); sampleIdEnd=sampleIdStart;
   keepOnlyTheXtoYSamples(samples,sampleIdStart,sampleIdEnd);
   keepOnlyValidSamples(samples);
   printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n");
   printf("Run on the following samples:\n");
   for(unsigned int s=0;s<samples.size();s++){samples[s].print();}
   printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n\n");

   vector<string> DataFileName;
   TH3F* dEdxTemplatesOld = NULL;
   for(unsigned int s=0;s<samples.size();s++){
      GetInputFiles(samples[s], BaseDirectory, DataFileName, 0);

      for(unsigned int f=0;f<DataFileName.size();f++){printf("file %i : %s\n", f, DataFileName[f].c_str());}

      bool isData   = (samples[s].Type==0);
      bool isMC     = (samples[s].Type==1);
      bool isSignal = (samples[s].Type>=2);
      bool is2016   = (samples[s].Name.find("13TeV16")==std::string::npos)?false:true;
      bool is2016G  = (samples[s].Name.find("13TeV16G")==std::string::npos)?false:true;

      dEdxTemplatesOld = loadDeDxTemplate(DIRNAME+(isData?"/../../../data/Data13TeV16_dEdxTemplate.root":"/../../../data/MC13TeV16_dEdxTemplate.root"), true);
      if(isData){  // 2016 values
         dEdxSF [0] = 1.00000;
         dEdxSF [1] = 1.464;
	 if (is2016G) dEdxSF[1] = 1.611;  //PostG
         //dEdxTemplates = loadDeDxTemplate(DIRNAME+"/../../../data/Data13TeV16_dEdxTemplate.root", true);
	 dEdxTemplates = loadDeDxTemplate((!is2016G)?(DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_RunPreG.root"):(DIRNAME+"/../../../data/dEdxTemplate_hit_SP_in_noC_CCC_RunPostG.root"), true);
      }else{  
         dEdxSF [0] = 1.09711;
         dEdxSF [1] = 1.09256;
         dEdxTemplates = loadDeDxTemplate(DIRNAME+"/../../../data/MC13TeV16_dEdxTemplate.root", true);
      }

      if(isData){    trackerCorrector.LoadDeDxCalibration(DIRNAME+"/../../../data/Data13TeVGains_v2.root"); 
      }else{ trackerCorrector.TrackerGains = NULL; //FIXME check gain for MC
      }

   moduleGeom::loadGeometry(DIRNAME+"/../../../data/CMS_GeomTree.root");
   muonTimingCalculator tofCalculator;
   tofCalculator.loadTimeOffset(DIRNAME+"/../../../data/MuonTimeOffset.txt");
   unsigned int CurrentRun = 0;

//   dedxHIPEmulator HIPemulator;

   fwlite::ChainEvent ev(DataFileName);
   printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
   printf("Looping on Tree              :");

   std::map<string, plotSt*>* MapTriggerPlots=NULL;;
   int NEvents = ev.size();// / NJobs;
   int FirstEvent = 0;//JobIndex * NEvents;
   int TreeStep = NEvents/50;if(TreeStep==0)TreeStep=1;
   double totalNEvents      = 0.0;
   double totalPreselEvents = 0.0;
   double totalPreselTracks = 0.0;
   for(Long64_t e=FirstEvent;e<FirstEvent+NEvents;e++){
      ev.to(e); 
      if(e%TreeStep==0){printf(".");fflush(stdout);}

      //if run change, update conditions
      if(CurrentRun != ev.eventAuxiliary().run()){
         totalNEvents     =0;
         totalPreselEvents=0;
	 totalPreselTracks=0;
         CurrentRun = ev.eventAuxiliary().run();
         tofCalculator.setRun(CurrentRun);
         trackerCorrector.setRun(CurrentRun);


         char DIRECTORY[2048];
         if(isData){
            sprintf(DIRECTORY,"%6i", ev.eventAuxiliary().run());
         }else{
            sprintf(DIRECTORY,"%s", samples[s].Name.c_str() );            
         }
         if(MapRunTriggerPlots.find(DIRECTORY)==MapRunTriggerPlots.end()){
            TDirectory* dir = OutputHisto;
            TDirectory::AddDirectory(kTRUE);
            TH1::AddDirectory(kTRUE);
            dir = (TDirectory*)OutputHisto->Get(DIRECTORY);
            if(dir==NULL){
               dir = OutputHisto->mkdir(DIRECTORY, DIRECTORY);
               dir->cd();

               for(unsigned int i=0;i<triggers.size();i++){
               for(unsigned int v=0;v<versions.size();v++){
                  MapRunTriggerPlots[DIRECTORY][triggers[i]+versions[v]] = new plotSt(triggers[i], versions[v]);
               }}
            }else{printf("BUG\n");}
         }
         MapTriggerPlots = &MapRunTriggerPlots[DIRECTORY];
      }
      if(!PassingTrigger(ev,"Any")){continue;} //need to pass at least one of the trigger, otherwise save time


      totalNEvents += 1;


      fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
      hscpCollHandle.getByLabel(ev,"HSCParticleProducer");
      if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
      susybsm::HSCParticleCollection hscpColl = *hscpCollHandle;

      fwlite::Handle<DeDxHitInfoAss> dedxCollH;
      dedxCollH.getByLabel(ev, "dedxHitInfo");
      if(!dedxCollH.isValid()){printf("Invalid dedxCollH\n");continue;}

      fwlite::Handle<MuonTimeExtraMap> TOFCollH;
      TOFCollH.getByLabel(ev, "muons",TOF_Label.c_str());
      if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}

      fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
      TOFDTCollH.getByLabel(ev, "muons",TOFdt_Label.c_str());
      if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");continue;}

      fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
      TOFCSCCollH.getByLabel(ev, "muons",TOFcsc_Label.c_str());
      if(!TOFCSCCollH.isValid()){printf("Invalid CSCTOF collection\n");continue;}

      fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
      vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
      if(!vertexCollHandle.isValid()){printf("Vertex Collection NotFound\n");continue;}
      const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;

      fwlite::Handle<CSCSegmentCollection> CSCSegmentCollHandle;
      fwlite::Handle<DTRecSegment4DCollection> DTSegmentCollHandle;            
      if(!isMC){ //do not reocmpute TOF on MC background
         CSCSegmentCollHandle.getByLabel(ev, "cscSegments");
         if(!CSCSegmentCollHandle.isValid()){printf("CSC Segment Collection not found!\n"); continue;}

         DTSegmentCollHandle.getByLabel(ev, "dt4DSegments");
         if(!DTSegmentCollHandle.isValid()){printf("DT Segment Collection not found!\n"); continue;}
      }

      HIPemulator.setEventRate(); //take it from a pdf

      bool plotPerEvent=true;
      bool alreadyPassed=false;
      for(unsigned int c=0;c<hscpColl.size();c++){
         susybsm::HSCParticle hscp  = hscpColl[c];
         reco::TrackRef track = hscp.trackRef();
         if(track.isNull())continue;
         reco::MuonRef muon = hscp.muonRef();         
         if(muon.isNull())continue; 
         if(!hscp.muonRef()->isStandAloneMuon())continue;

         const DeDxHitInfo* dedxHits = NULL;
         if(TypeMode!=3 && !track.isNull()) {
            DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());		 
            if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
         }

         bool useClusterCleaning = true;
         DeDxData PSdedxSObj = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, TypeMode==5, false, trackerCorrector.TrackerGains, true, true, 99, false, 1);
         DeDxData PSdedxMObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1);
         const reco::MuonTimeExtra* PStof = NULL;
         const reco::MuonTimeExtra* PSdttof = NULL;
         const reco::MuonTimeExtra* PScsctof = NULL;        
         if(!hscp.muonRef().isNull() && hscp.muonRef()->isStandAloneMuon() ){
            if(isMC){
               PStof  = &TOFCollH->get(hscp.muonRef().key()); PSdttof = &TOFDTCollH->get(hscp.muonRef().key());  PScsctof = &TOFCSCCollH->get(hscp.muonRef().key());
            }else{
               const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
               const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;
               tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, 1 ); //apply T0 correction on data but not on signal MC
               PStof  = &tofCalculator.combinedTOF; PSdttof = &tofCalculator.dtTOF;  PScsctof = &tofCalculator.cscTOF;
            }
         }
         
//         for(unsigned int i=0;i<triggers.size();i++){ //PtDistribution before preselection - i.e. out of the box
//            if(!PassingTrigger(ev,triggers[i])){continue;}
//            for(unsigned int v=0;v<versions.size();v++){
//               plotSt* plots = (*MapTriggerPlots)[triggers[i]+versions[v]];
//               plots->PtBPS->Fill(hscp.trackRef()->ptError());
//            }
//         }

         if(!PassPreselection(hscp, PSdedxSObj, PSdedxMObj, PStof, PSdttof, PScsctof, ev)){continue;}
         if(!alreadyPassed) {totalPreselEvents += 1; alreadyPassed=true;}
	 totalPreselTracks += 1;
         for(unsigned int i=0;i<triggers.size();i++){
            if(!PassingTrigger(ev,triggers[i])){continue;}
            for(unsigned int v=0;v<versions.size();v++){

            plotSt* plots = (*MapTriggerPlots)[triggers[i]+versions[v]];
            std::map<unsigned int, TH1D**>::iterator dEdxHitPerLumiIt = plots->dEdxHitPerLumi.find(ev.eventAuxiliary().luminosityBlock());
            if(dEdxHitPerLumiIt==plots->dEdxHitPerLumi.end()){
               char LUMI[256];sprintf(LUMI,"_ls%i", ev.eventAuxiliary().luminosityBlock());
               plots->dEdxHitPerLumi[ev.eventAuxiliary().luminosityBlock()] = new TH1D*[3];
               dEdxHitPerLumiIt = plots->dEdxHitPerLumi.find(ev.eventAuxiliary().luminosityBlock());
               (dEdxHitPerLumiIt->second)[0] = (TH1D*)plots->NVert       ->Clone((string(plots->NVert       ->GetName())+LUMI).c_str()); 
               (dEdxHitPerLumiIt->second)[1] = (TH1D*)plots->dEdxHitPixel->Clone((string(plots->dEdxHitPixel->GetName())+LUMI).c_str());
               (dEdxHitPerLumiIt->second)[2] = (TH1D*)plots->dEdxHitStrip->Clone((string(plots->dEdxHitStrip->GetName())+LUMI).c_str());               
//               (dEdxHitPerLumiIt->second)[3] = (TH1D*)plots->dEdxHitPixel->Clone((string(plots->dEdxHitPixel->GetName())+"PIB"+LUMI).c_str());
//               (dEdxHitPerLumiIt->second)[4] = (TH1D*)plots->dEdxHitPixel->Clone((string(plots->dEdxHitPixel->GetName())+"PIE"+LUMI).c_str());
//               (dEdxHitPerLumiIt->second)[5] = (TH1D*)plots->dEdxHitStrip->Clone((string(plots->dEdxHitStrip->GetName())+"TIB"+LUMI).c_str());
//               (dEdxHitPerLumiIt->second)[6] = (TH1D*)plots->dEdxHitStrip->Clone((string(plots->dEdxHitStrip->GetName())+"TID"+LUMI).c_str());
//               (dEdxHitPerLumiIt->second)[7] = (TH1D*)plots->dEdxHitStrip->Clone((string(plots->dEdxHitStrip->GetName())+"TOB"+LUMI).c_str());
//               (dEdxHitPerLumiIt->second)[8] = (TH1D*)plots->dEdxHitStrip->Clone((string(plots->dEdxHitStrip->GetName())+"TEC"+LUMI).c_str());
            }


            bool useClusterCleaning = true;

            auto tkGains = trackerCorrector.TrackerGains;
            if(versions[v]=="AOD")tkGains=NULL;
 
            bool fake=false;
            if(versions[v]=="FAKE" && !isData)fake=true;

            HitDeDxCollection hitDeDx = getHitDeDx(dedxHits, dEdxSF, tkGains, false, 1);
            if(fake)HIPemulator.fakeHIP(hitDeDx);
            for(unsigned int h=0;h<hitDeDx.size();h++){
               if((useClusterCleaning && !hitDeDx[h].passClusterCleaning) || !hitDeDx[h].isInside)continue;
               if(hitDeDx[h].subDet< 3){plots->dEdxHitPixel->Fill(hitDeDx[h].dedx); (dEdxHitPerLumiIt->second)[1]->Fill(hitDeDx[h].dedx);}
               if(hitDeDx[h].subDet>=3){plots->dEdxHitStrip->Fill(hitDeDx[h].dedx); (dEdxHitPerLumiIt->second)[2]->Fill(hitDeDx[h].dedx);}
               //(dEdxHitPerLumiIt->second)[2+hitDeDx[h].subDet]->Fill(hitDeDx[h].dedx); 
            }

            DeDxData dedxMin1Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, tkGains, true, true, 99, false, 1, 0.15, fake?&HIPemulator:NULL);
            DeDxData dedxMin2Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, tkGains, true, true, 99, false, 1, 0.2, fake?&HIPemulator:NULL);
            DeDxData dedxMin3Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, tkGains, true, true, 99, false, 1, 0.3, fake?&HIPemulator:NULL);
            DeDxData dedxMin4Obj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, tkGains, true, true, 99, false, 1, 0.4, fake?&HIPemulator:NULL);
            DeDxData dedxSObj    = computedEdx(dedxHits, dEdxSF, dEdxTemplates   , true, useClusterCleaning, TypeMode==5, false, tkGains, true, true, 99, false, 1, 0.0, fake?&HIPemulator:NULL);
            DeDxData dedxSObjOld = computedEdx(dedxHits, dEdxSF, dEdxTemplatesOld, true, useClusterCleaning, TypeMode==5, false, tkGains, true, true, 99, false, 1, 0.0, fake?&HIPemulator:NULL);
            DeDxData dedxMObj  = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, tkGains, true, true, 99, false, 1, 0.0, fake?&HIPemulator:NULL);
            DeDxData dedxMTObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , true, tkGains, true, true, 99, false, 1, 0.0, fake?&HIPemulator:NULL);
            DeDxData dedxMSObj = computedEdx(dedxHits, dEdxSF, NULL,          false,useClusterCleaning, false      , false, tkGains, true, true, 99, false, 1, 0.0, fake?&HIPemulator:NULL);
            DeDxData dedxMPObj = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, tkGains, false, true, 99, false, 1, 0.0, fake?&HIPemulator:NULL);

            const reco::MuonTimeExtra* tof = NULL;
            const reco::MuonTimeExtra* dttof = NULL;
            const reco::MuonTimeExtra* csctof = NULL;
            if(versions[v]=="AOD" || isMC){
               if(!hscp.muonRef().isNull()){ tof  = &TOFCollH->get(hscp.muonRef().key()); dttof  = &TOFDTCollH->get(hscp.muonRef().key()); csctof  = &TOFCSCCollH->get(hscp.muonRef().key());}
            }else{
               if(!hscp.muonRef().isNull() && hscp.muonRef()->isStandAloneMuon() ){
                  const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
                  const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;
                  tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, 1 ); //apply T0 correction on data but not on signal MC
                  tof  = &tofCalculator.combinedTOF; dttof = &tofCalculator.dtTOF;  csctof = &tofCalculator.cscTOF;
               }          
            }

            if(plotPerEvent){plots->NVert->Fill(vertexColl.size()); (dEdxHitPerLumiIt->second)[0]->Fill(vertexColl.size());}
            plots->Pt   ->Fill(hscp.trackRef()->pt());
            plots->PtEta->Fill(std::fabs(hscp.trackRef()->eta()), hscp.trackRef()->pt());
//            plots->PtErr->Fill(hscp.trackRef()->ptError());

            plots->dEdxMin1->Fill(dedxMin1Obj.dEdx());
            plots->dEdxMin2->Fill(dedxMin2Obj.dEdx());
            plots->dEdxMin3->Fill(dedxMin3Obj.dEdx());
            plots->dEdxMin4->Fill(dedxMin4Obj.dEdx());
            plots->dEdx    ->Fill(dedxSObj.dEdx());
            plots->dEdxOld ->Fill(dedxSObjOld.dEdx());
            plots->dEdxMT  ->Fill(dedxMTObj.dEdx());
            plots->dEdxM   ->Fill(dedxMObj.dEdx());
            plots->dEdxMS  ->Fill(dedxMSObj.dEdx());
            plots->dEdxMP  ->Fill(dedxMPObj.dEdx());
            if(fabs(track->eta())<0.5){
            plots->dEdxMSC->Fill(dedxMSObj.dEdx());
            plots->dEdxMPC->Fill(dedxMPObj.dEdx());
            }
            if(fabs(track->eta())>1.5){
            plots->dEdxMSF->Fill(dedxMSObj.dEdx());
            plots->dEdxMPF->Fill(dedxMPObj.dEdx());
            }

            if(tof && tof->nDof()>=GlobalMinNDOF && (dttof->nDof()>=GlobalMinNDOFDT || csctof->nDof()>=GlobalMinNDOFCSC) && tof->inverseBetaErr()<=GlobalMaxTOFErr && fabs(dttof->inverseBeta()-1)<50){
               plots->TOF->Fill(tof->inverseBeta());
               if(dttof->nDof()>=GlobalMinNDOFDT) plots->TOFDT->Fill(dttof->inverseBeta());
               if(csctof->nDof()>=GlobalMinNDOFCSC) plots->TOFCSC->Fill(csctof->inverseBeta());
               plots->Vertex->Fill(tof->timeAtIpInOut());
               if(dttof->nDof()>=GlobalMinNDOFDT) plots->VertexDT->Fill(dttof->timeAtIpInOut());
               if(csctof->nDof()>=GlobalMinNDOFCSC) plots->VertexCSC->Fill(csctof->timeAtIpInOut());
            } 
           plots->PreselEff->SetBinContent(1, totalPreselEvents/totalNEvents);
	   plots->TPsVNEvts->SetBinContent(1, totalPreselTracks/totalNEvents);
           plots->PreselEff->SetBinError  (1, sqrt(totalPreselEvents + totalNEvents));
	   plots->TPsVNEvts->SetBinError  (1, sqrt(totalPreselTracks + totalNEvents));
           } 
         }
         plotPerEvent = false;
      }
   }printf("\n");
   }

   OutputHisto->Write();
   OutputHisto->Close();  
}



bool IncreasedTreshold(const trigger::TriggerEvent& trEv, const edm::InputTag& InputPath, double NewThreshold, int NObjectAboveThreshold, bool averageThreshold)
{
   unsigned int filterIndex = trEv.filterIndex(InputPath);
   //if(filterIndex<trEv.sizeFilters())printf("SELECTED INDEX =%i --> %s    XXX   %s\n",filterIndex,trEv.filterTag(filterIndex).label().c_str(), trEv.filterTag(filterIndex).process().c_str());
         
   if (filterIndex<trEv.sizeFilters()){
      const trigger::Vids& VIDS(trEv.filterIds(filterIndex));
      const trigger::Keys& KEYS(trEv.filterKeys(filterIndex));
      const size_type nI(VIDS.size());
      const size_type nK(KEYS.size());
      assert(nI==nK);
      const size_type n(max(nI,nK));
      const trigger::TriggerObjectCollection& TOC(trEv.getObjects());


      if(!averageThreshold){
         int NObjectAboveThresholdObserved = 0;
         for (size_type i=0; i!=n; ++i) {
            const TriggerObject& TO(TOC[KEYS[i]]);
            if(TO.pt()> NewThreshold) NObjectAboveThresholdObserved++;
   	    //cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "<< TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()<< endl;
         }          
         if(NObjectAboveThresholdObserved>=NObjectAboveThreshold)return true;

      }else{
         std::vector<double> ObjPt;

         for (size_type i=0; i!=n; ++i) {
            const TriggerObject& TO(TOC[KEYS[i]]);
            ObjPt.push_back(TO.pt());
            //cout << "   " << i << " " << VIDS[i] << "/" << KEYS[i] << ": "<< TO.id() << " " << TO.pt() << " " << TO.eta() << " " << TO.phi() << " " << TO.mass()<< endl;
         }  
         if((int)(ObjPt.size())<NObjectAboveThreshold)return false;
         std::sort(ObjPt.begin(), ObjPt.end());
         
         double Average = 0;
         for(int i=0; i<NObjectAboveThreshold;i++){
            Average+= ObjPt[ObjPt.size()-1-i];            
         }Average/=NObjectAboveThreshold;
	 //cout << "AVERAGE = " << Average << endl;
         
         if(Average>NewThreshold)return true;                  
      }
   }
   return false;
}


