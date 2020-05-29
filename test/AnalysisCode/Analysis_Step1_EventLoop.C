// Original Author:  Loic Quertenmont


namespace reco    { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra; class PFMET; class HitPattern;}
namespace susybsm { class HSCParticle; class HSCPIsolation; class MuonSegment; class HSCPDeDxInfo;}
namespace fwlite  { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm     { class TriggerResults; class TriggerResultsByName; class InputTag; class LumiReWeighting;}
namespace reweight{ class PoissonMeanShifter;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h" 

#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPDeDxInfo.h"
#include "DataFormats/TrackReco/interface/DeDxHitInfo.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/MuonSegment.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/METReco/interface/PFMETCollection.h"
#include "DataFormats/METReco/interface/PFMET.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;
using namespace reweight;

#endif


//the define here is simply need to load FWLITE code from the include
#define FWLITE
#include "Analysis_Global.h"
#include "Analysis_CommonFunction.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_PlotStructure.h"
#include "Analysis_Samples.h"
#include "Analysis_TOFUtility.h"
#include "tdrstyle.C"

/////////////////////////// FUNCTION DECLARATION /////////////////////////////

void InitHistos(stPlots* st=NULL);
void Analysis_Step1_EventLoop(char* SavePath);


bool PassTrigger(const fwlite::ChainEvent& ev, bool isData, bool isCosmic=false, L1BugEmulator* emul=NULL);
bool   PassPreselection(const susybsm::HSCParticle& hscp, const DeDxHitInfo* dedxHits, const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st=NULL, const double& GenBeta=-1, bool RescaleP=false, const double& RescaleI=0.0, const double& RescaleT=0.0, double MassErr=-1);
bool PassSelection(const susybsm::HSCParticle& hscp,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const fwlite::ChainEvent& ev, const int& CutIndex=0, stPlots* st=NULL, const bool isFlip=false, const double& GenBeta=-1, bool RescaleP=false, const double& RescaleI=0.0, const double& RescaleT=0.0);
void Analysis_FillControlAndPredictionHist(const susybsm::HSCParticle& hscp, const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, stPlots* st=NULL);
double SegSep(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, double& minPhi, double& minEta);
double RescaledPt(const double& pt, const double& eta, const double& phi, const int& charge);
int  muonStations(const reco::HitPattern& hitPattern);
double scaleFactor(double eta);
/////////////////////////// VARIABLE DECLARATION /////////////////////////////

bool isMCglobal = false;
float Event_Weight = 1;
int   MaxEntry = -1;

TFile* HistoFile;

std::vector<double>  CutPt ;
std::vector<double>  CutI  ;
std::vector<double>  CutTOF;

TProfile*  HCuts_Pt;
TProfile*  HCuts_I;
TProfile*  HCuts_TOF;

//The cuts used for the check on the background prediction by tracks with TOF<1
std::vector<double>  CutPt_Flip ;
std::vector<double>  CutI_Flip  ;
std::vector<double>  CutTOF_Flip;

TProfile*  HCuts_Pt_Flip;
TProfile*  HCuts_I_Flip;
TProfile*  HCuts_TOF_Flip;

std::vector<stSample> samples;
std::vector<stSample> samplesFull;
std::map<std::string, stPlots> plotsMap;

std::vector< float > BgLumiMC; //MC                                           
std::vector< float > TrueDist;
std::vector< float > TrueDistSyst; 
edm::LumiReWeighting LumiWeightsMC;
edm::LumiReWeighting LumiWeightsMCSyst;
//reweight::PoissonMeanShifter PShift(0.6);//0.6 for upshift, -0.6 for downshift

TH3F* dEdxTemplates = NULL;
dedxGainCorrector trackerCorrector;
double dEdxSF [2] = {
   1.00000,   // [0]  unchanged
   1.464//1.21836    // [1]  Pixel data to SiStrip data
};
dedxHIPEmulator HIPemulator;
dedxHIPEmulator HIPemulatorUp (false, "ratePdfPixel_Up", "ratePdfStrip_Up");
dedxHIPEmulator HIPemulatorDown (false, "ratePdfPixel_Down", "ratePdfStrip_Down");
L1BugEmulator        L1Emul;
HIPTrackLossEmulator HIPTrackLossEmul;

bool useClusterCleaning = true;
/////////////////////////// CODE PARAMETERS /////////////////////////////

void Analysis_Step1_EventLoop(string MODE="COMPILE", int TypeMode_=0, string InputSampleName="")
{
   if(MODE=="COMPILE")return;

   //setup ROOT global variables (mostly cosmetic and histo in file treatment)
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.05);
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadRightMargin (0.18);
   gStyle->SetPadLeftMargin  (0.13);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);
   TH1::AddDirectory(kTRUE);

   // redefine global variable dependent on the arguments given to the function
   TypeMode       = TypeMode_;

   //AnalysisType dependent cuts
   if(TypeMode<2){  
     GlobalMinNDOF      = 0; 
     GlobalMinTOF       = 0;
   }else if(TypeMode==3){
     GlobalMinPt        =      80;
     GlobalMaxDZ        =    15.0;
     GlobalMaxDXY       =    15.0;
     GlobalMaxV3D       =  999999;
     GlobalMinIs        =      -1;
     IPbound            =     150;
     PredBins           =       6;
   }else if(TypeMode==4){
     GlobalMaxEIsol     =  999999;   // cut on calorimeter isolation (E/P)
     useClusterCleaning = false; //switch off cluster cleaning for mCHAMPs
   } else if(TypeMode==5){
     IPbound            = 4.5;
     //GlobalMinIm        = 2.8; //is actually dEdx max at skim level (reverse logic for type5)
     GlobalMinIm        = 999999;
     GlobalMinNDOF      = 0; //tkOnly analysis --> comment these 2 lines to use only global muon tracks
     GlobalMinTOF       = 0;
   }
   
   // define the selection to be considered later for the optimization
   // WARNING: recall that this has a huge impact on the analysis time AND on the output file size --> be carefull with your choice
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
   //   for (int CutIndex = 0; CutIndex < CutPt.size(); ++CutIndex)      printf("%4.0i  %3.0f   %3.3f   %3.3f\n", CutIndex+1, CutPt[CutIndex], CutI[CutIndex],  CutTOF[CutIndex]);
   //   for (int CutIndex = 0; CutIndex < CutPt_Flip.size(); ++CutIndex) printf("%4.0i  %3.0f   %3.3f   %3.3f\n", CutIndex+1, CutPt_Flip[CutIndex], CutI_Flip[CutIndex],  CutTOF_Flip[CutIndex]);

   //make the directory structure corresponding to this analysis (depends on dEdx/TOF estimator being used, Eta/Pt cuts and Mode of the analysis)
   char Buffer[2048], Command[2048];
   sprintf(Buffer,"Results/Type%i/", TypeMode);
   sprintf(Command,"mkdir -p %s",Buffer); system(Command);

   // get all the samples and clean the list to keep only the one we want to run on... Also initialize the BaseDirectory
   printf(" InitBaseDirectory() <-------------\n");
   InitBaseDirectory();
   GetSampleDefinition(samples);
   samplesFull = samples;
   if(MODE.find("ANALYSE_")==0){
      int sampleIdStart, sampleIdEnd; sscanf(MODE.c_str(),"ANALYSE_%d_to_%d",&sampleIdStart, &sampleIdEnd);
      keepOnlyTheXtoYSamples(samples,sampleIdStart,sampleIdEnd);
      keepOnlyValidSamples(samples);
      printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n");
      printf("Run on the following samples:\n");
      for(unsigned int s=0;s<samples.size();s++){samples[s].print();}
      printf("----------------------------------------------------------------------------------------------------------------------------------------------------\n\n");
   }else{
      printf("You must select a MODE:\n");
      printf("MODE='ANALYSE_X_to_Y'   : Will run the analysis on the samples with index in the range [X,Y]\n"); 
      return;
   }

std::cout<<"A\n";

   //initialize LumiReWeighting
   //FIXME  pileup scenario must be updated based on data/mc
   bool is2016  = (samples[0].Name.find("13TeV16") !=std::string::npos)?true:false;
   bool is2016G = (samples[0].Name.find("13TeV16G")!=std::string::npos)?true:false;
   HIPemulator.    setPeriodHIPRate(is2016G);
   HIPemulatorUp.  setPeriodHIPRate(is2016G);
   HIPemulatorDown.setPeriodHIPRate(is2016G);
   if(samples[0].Pileup=="S15"){        for(int i=0; i<100; ++i) BgLumiMC.push_back(Pileup_MC_Startup2015_25ns[i]);
   }else if(samples[0].Pileup=="NoPU" && !is2016 && !is2016G){ for(int i=0; i<100; ++i) BgLumiMC.push_back(TrueDist2015_f[i]); //Push same as 2015 data to garantee no PU reweighting
   }else if(samples[0].Pileup=="NoPU" && is2016 && !is2016G) { for(int i=0; i<100; ++i) BgLumiMC.push_back(TrueDist2016_f[i]); //Push same as 2016 data to garantee no PU reweighting
   }else if(samples[0].Pileup=="NoPU" && is2016 && is2016G) { for(int i=0; i<100; ++i) BgLumiMC.push_back(TrueDist2016G_f[i]); //Push same as 2016 data to garantee no PU reweighting
   }else if (samples[0].Pileup=="S10"){ for(int i=0; i<100; ++i) BgLumiMC.push_back(Pileup_MC_Summer2012[i]);
   }else{                               for(int i=0; i<100; ++i) BgLumiMC.push_back(Pileup_MC_Fall11[i]);
   }
std::cout<<"A1\n";

   if (!is2016 && !is2016G){
      for(int i=0; i<100; ++i) TrueDist    .push_back(TrueDist2015_f[i]);
      for(int i=0; i<100; ++i) TrueDistSyst.push_back(TrueDist2015_XSecShiftUp_f[i]);
   } else if (is2016 && !is2016G){
      for(int i=0; i<100; ++i) TrueDist    .push_back(TrueDist2016_f[i]);
      for(int i=0; i<100; ++i) TrueDistSyst.push_back(TrueDist2016_XSecShiftUp_f[i]);
   }  else if (is2016 && is2016G){
      for(int i=0; i<100; ++i) TrueDist    .push_back(TrueDist2016G_f[i]);
      for(int i=0; i<100; ++i) TrueDistSyst.push_back(TrueDist2016G_XSecShiftUp_f[i]);
   }
   
std::cout<<"A2\n";

   LumiWeightsMC     = edm::LumiReWeighting(BgLumiMC, TrueDist);
std::cout<<"A3\n";

   LumiWeightsMCSyst = edm::LumiReWeighting(BgLumiMC, TrueDistSyst);

std::cout<<"B\n";


   //create histogram file and run the analyis
   HistoFile = new TFile((string(Buffer)+"/Histos_"+samples[0].Name+"_"+ReplacePartOfString(samples[0].FileName, "/", "_")+".root").c_str(),"RECREATE");
std::cout<<"C\n";

   Analysis_Step1_EventLoop(Buffer);
std::cout<<"Z\n";

   HistoFile->Write();
   HistoFile->Close();
   return;
}

TVector3 getOuterHitPos(const DeDxHitInfo* dedxHits){
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

// check if the event is passing trigger or not --> note that the function has two part (one for 2011 analysis and the other one for 2012)
bool PassTrigger(const fwlite::ChainEvent& ev, bool isData, bool isCosmic, L1BugEmulator* emul)
{
   edm::TriggerResultsByName tr = ev.triggerResultsByName("HLT");
   if(!tr.isValid())         tr = ev.triggerResultsByName("MergeHLT");
   if(!tr.isValid())return false;

   //mk TrigInfo - as global var
   TrigInfo=0;
   
   bool metTrig = false;
   bool muTrig = false;
   

   if(passTriggerPatterns(tr, "HLT_PFMET170_NoiseCleaned_v*") || passTriggerPatterns(tr, "HLT_PFMET170_HBHECleaned_v*")) metTrig = true;
   if(passTriggerPatterns(tr, "HLT_Mu45_eta2p1_v*") || passTriggerPatterns(tr, "HLT_Mu50_v*")) muTrig = true;

   if (!metTrig && muTrig) TrigInfo = 1;
   if (metTrig && !muTrig) TrigInfo = 2;
   if (metTrig && muTrig)  TrigInfo = 3;
   //unsigned int TrigInfo =0; //1 -mu only, 2- met only, 3 mu and met 



   if(passTriggerPatterns(tr, "HLT_PFMET170_NoiseCleaned_v*") || passTriggerPatterns(tr, "HLT_PFMET170_HBHECleaned_v*"))return true;
   if(passTriggerPatterns(tr, "HLT_Mu45_eta2p1_v*") || passTriggerPatterns(tr, "HLT_Mu50_v*")){
      if (!isData && emul){
         fwlite::Handle < std::vector<reco::Muon> > muonCollHandle;
         muonCollHandle.getByLabel(ev, "muons");
         if (!muonCollHandle.isValid()) return false;
         else{
            bool KeepEvent=false;
            for (unsigned int c=0;c<muonCollHandle->size();c++){
               reco::MuonRef muon = reco::MuonRef(muonCollHandle.product(), c);
               if (muon.isNull()) continue;
               if (muon->track().isNull()) continue;
               if (emul->PassesL1Inefficiency(muon->track()->pt(), std::fabs(muon->track()->eta()))){
                  KeepEvent=true;
                  break;
               }
            }
            return KeepEvent;
         }
      }
      else return true;
   }
//   if(passTriggerPatterns(tr, "HLT_Mu50_v*"))return true;

   return false; //FIXME triggers bellow will need to be adapted based on Run2 trigger menu

   //for(unsigned int i=0;i<tr.size();i++){
   //printf("Path %3i %50s --> %1i\n",i, tr.triggerName(i).c_str(),tr.accept(i));
   //}fflush(stdout);

   //if(tr.accept("HSCPHLTTriggerMetDeDxFilter"))return true;
   //if(tr.accept("HSCPHLTTriggerMuDeDxFilter"))return true;
   if(tr.accept("HSCPHLTTriggerMuFilter"))return true;
   if(tr.accept("HSCPHLTTriggerPFMetFilter"))return true;

   //Could probably use this trigger for the other analyses as well
   if(TypeMode==3){
      if(tr.size()== tr.triggerIndex("HSCPHLTTriggerL2MuFilter")) return false;
      if(tr.accept(tr.triggerIndex("HSCPHLTTriggerL2MuFilter")))  return true;

      //Only accepted if looking for cosmic events
      if(isCosmic) {
         if(tr.size()== tr.triggerIndex("HSCPHLTTriggerCosmicFilter")) return false;
         if(tr.accept(tr.triggerIndex("HSCPHLTTriggerCosmicFilter"))) return true;
      }
   }
   return false;
}


// check if one HSCP candidate is passing the preselection (the function also has many more arguments because it is used to fill some histograms AND to evaluate the systematics
double OpenAngle = -1; //global variable needed by PassPreselection... Ugly isn't it?!
double TreeDXY = -1;
double TreeDZ = -1;
bool isCosmicSB = false;
bool isSemiCosmicSB = false;
bool PassPreselection(const susybsm::HSCParticle& hscp, const DeDxHitInfo* dedxHits,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const reco::MuonTimeExtra* dttof, const reco::MuonTimeExtra* csctof, const fwlite::ChainEvent& ev, stPlots* st, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT, double MassErr)
{
   if(TypeMode==1 && !(hscp.type() == HSCParticleType::trackerMuon || hscp.type() == HSCParticleType::globalMuon))return false;
   if( (TypeMode==2 || TypeMode==4) && hscp.type() != HSCParticleType::globalMuon)return false;

   reco::TrackRef   track;
   reco::MuonRef muon = hscp.muonRef();

   if(TypeMode!=3) track = hscp.trackRef();
   else {
     if(muon.isNull()) return false;
     track = muon->standAloneMuon();
   }
   if(track.isNull())return false;

   if(st){st->Total->Fill(0.0,Event_Weight);
     if(GenBeta>=0)st->Beta_Matched->Fill(GenBeta, Event_Weight);
     st->BS_Eta->Fill(track->eta(),Event_Weight);
   }

   if(fabs(track->eta())>GlobalMaxEta) return false;
 //mk_ to byla jakas proba   if(fabs(track->eta())<1.2) return false;   
 //mk_ a new cut on eta
 
//max PtCut

//  if(TypeMode==0) if( track->pt*cosh(track->eta() ) < 1700 ) return false;
//  if(TypeMode==2) if( track->pt*cosh(track->eta() ) <  850 ) return false;

// no max PCut
//  if(TypeMode==0) if( track->p() > 1700 ) return false;
//  if(TypeMode==2) if( track->p() >  850 ) return false;




   //Cut on number of matched muon stations
   int count = muonStations(track->hitPattern());
   if(st) {
     st->BS_MatchedStations->Fill(count, Event_Weight);
   }

   if(TypeMode==3 && count<minMuStations) return false;
   if(st) st->Stations->Fill(0.0, Event_Weight);


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
     if(st) st->BS_dzAll->Fill( track->dz (vertexColl[i].position()),Event_Weight);
     if(st) st->BS_dxyAll->Fill(track->dxy(vertexColl[i].position()),Event_Weight);
//     if(highestPtGoodVertex<0){
//        printf("debug Nvert=%3i vertIndex=%3i dxy=%+6.2f dz=%+6.2f\n", (int)vertexColl.size(), (int) i, track->dxy (vertexColl[i].position()), track->dz (vertexColl[i].position()));
//     }


//     if(highestPtGoodVertex<0)highestPtGoodVertex = i;
     if(fabs(track->dz (vertexColl[i].position())) < fabs(dzMin) ){
         dzMin = fabs(track->dz (vertexColl[i].position()));
         highestPtGoodVertex = i;
//       dz  = track->dz (vertexColl[i].position());
//       dxy = track->dxy(vertexColl[i].position());
     }
   }if(highestPtGoodVertex<0)highestPtGoodVertex=0;

   if(st){st->BS_NVertex->Fill(vertexColl.size(), Event_Weight);
     st->BS_NVertex_NoEventWeight->Fill(vertexColl.size());
   }
   double dz  = track->dz (vertexColl[highestPtGoodVertex].position());
   double dxy = track->dxy(vertexColl[highestPtGoodVertex].position());
 
   bool PUA = (vertexColl.size()<15);
   bool PUB = (vertexColl.size()>=15);

   if(st){st->BS_TNOH->Fill(track->found(),Event_Weight);
          if(PUA)st->BS_TNOH_PUA->Fill(track->found(),Event_Weight);
          if(PUB)st->BS_TNOH_PUB->Fill(track->found(),Event_Weight);
          st->BS_TNOHFraction->Fill(track->validFraction(),Event_Weight);
	  st->BS_TNOPH->Fill(track->hitPattern().numberOfValidPixelHits(),Event_Weight);
   }

   if(TypeMode!=3 && track->found()<GlobalMinNOH)return false;

   if(TypeMode!=3 && track->hitPattern().numberOfValidPixelHits()<GlobalMinNOPH)return false;
   if(TypeMode!=3 && track->validFraction()<GlobalMinFOVH)return false;

   int missingHitsTillLast = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);;
   double validFractionTillLast = track->found()<=0?-1:track->found() / float(track->found() + missingHitsTillLast);
  
   if(st){st->BS_TNOHFractionTillLast->Fill(validFractionTillLast,Event_Weight);
	  st->BS_TNOMHTillLast->Fill(missingHitsTillLast,Event_Weight);
   }

   if(TypeMode!=3 && missingHitsTillLast>GlobalMaxNOMHTillLast)return false;
   if(TypeMode!=3 && validFractionTillLast<GlobalMinFOVHTillLast)return false;

   if(st){st->TNOH  ->Fill(0.0,Event_Weight);
     if(dedxSObj){
         st->BS_TNOM->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);
         if(track->found() - dedxSObj->numberOfMeasurements()) 
             st->BS_EtaNBH->Fill(track->eta(), track->found() - dedxSObj->numberOfMeasurements(), Event_Weight);
         if(PUA)st->BS_TNOM_PUA->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);
         if(PUB)st->BS_TNOM_PUB->Fill(dedxSObj->numberOfMeasurements(),Event_Weight);
     }
   }
   if(dedxSObj) if(dedxSObj->numberOfMeasurements()<GlobalMinNOM)return false;
   if(st){st->TNOM  ->Fill(0.0,Event_Weight);}

   if(tof){
   if(st){st->BS_nDof->Fill(tof->nDof(),Event_Weight);}
   if((TypeMode>1  && TypeMode!=5) && tof->nDof()<GlobalMinNDOF && (dttof->nDof()<GlobalMinNDOFDT || csctof->nDof()<GlobalMinNDOFCSC) )return false;
   }

   if(st){st->nDof  ->Fill(0.0,Event_Weight);
          st->BS_Qual->Fill(track->qualityMask(),Event_Weight);
   }

   if(TypeMode!=3 && track->qualityMask()<GlobalMinQual )return false; // FIXME Tracks with quality > 2 are bad also!
//   if(TypeMode!=3 && track->qualityMask() != FixedQual)return false; // FIXME if this is true, no tracks pass eventually ... so what now?
   if(st){st->Qual  ->Fill(0.0,Event_Weight);
          st->BS_Chi2->Fill(track->chi2()/track->ndof(),Event_Weight);
   }
   if(TypeMode!=3 && track->chi2()/track->ndof()>GlobalMaxChi2 )return false;
   if(st){st->Chi2  ->Fill(0.0,Event_Weight);}

   if(st && GenBeta>=0)st->Beta_PreselectedA->Fill(GenBeta, Event_Weight);

   if(st){st->BS_MPt ->Fill(track->pt(),Event_Weight);}
   if(RescaleP){ if(RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())<GlobalMinPt)return false;
   }else{        if(track->pt()<GlobalMinPt)return false;   }

   if(st){st->MPt   ->Fill(0.0,Event_Weight);
     if(dedxSObj) st->BS_MIs->Fill(dedxSObj->dEdx(),Event_Weight);
     if(dedxMObj) st->BS_MIm->Fill(dedxMObj->dEdx(),Event_Weight);
   }

   if(dedxSObj && dedxSObj->dEdx()+RescaleI<GlobalMinIs)return false;
   if(dedxMObj && ((TypeMode!=5 && dedxMObj->dEdx()<GlobalMinIm) || (TypeMode==5 && dedxMObj->dEdx()>GlobalMinIm)) )return false;
   if(st){st->MI   ->Fill(0.0,Event_Weight);}

   if(tof){
   if(st){st->BS_MTOF ->Fill(tof->inverseBeta(),Event_Weight);}
   //This cut is no longer applied here but rather in the PassSelection part to use the region
   //with TOF<GlobalMinTOF as a background check
   //if(TypeMode>1 && tof->inverseBeta()+RescaleT<GlobalMinTOF)return false;

   if(st)st->BS_TOFError->Fill(tof->inverseBetaErr(),Event_Weight);
   if((TypeMode>1  && TypeMode!=5) && tof->inverseBetaErr()>GlobalMaxTOFErr)return false;

   if(st) st->BS_TimeAtIP->Fill(tof->timeAtIpInOut(),Event_Weight);
   if(TypeMode==3 && min(min(fabs(tof->timeAtIpInOut()-100), fabs(tof->timeAtIpInOut()-50)), min(fabs(tof->timeAtIpInOut()+100), fabs(tof->timeAtIpInOut()+50)))<5) return false;
   }

   if(st) st->BS_dzMinv3d->Fill(dz,Event_Weight);
   if(st) st->BS_dxyMinv3d->Fill(dxy,Event_Weight);
   if(st) st->BS_PV->Fill(goodVerts,Event_Weight);   
   if(st) st->BS_PV_NoEventWeight->Fill(goodVerts);
   if(st && dedxSObj) st->BS_NOMoNOHvsPV->Fill(goodVerts,dedxSObj->numberOfMeasurements()/(double)track->found(),Event_Weight);

   //Require at least one good vertex except if cosmic event
   if(TypeMode==3 && goodVerts<1 && (!st || st->Name.find("Cosmic")==string::npos)) return false;

   //For TOF only analysis match to a SA track without vertex constraint for IP cuts
   if(TypeMode==3) {
     fwlite::Handle< std::vector<reco::Track> > noVertexTrackCollHandle;
     noVertexTrackCollHandle.getByLabel(ev,"refittedStandAloneMuons", "");

     //Find closest NV track
     const std::vector<reco::Track>& noVertexTrackColl = *noVertexTrackCollHandle;
     reco::Track NVTrack;
     double minDr=15;
     for(unsigned int i=0;i<noVertexTrackColl.size();i++){
       double dR = deltaR(track->eta(), track->phi(), noVertexTrackColl[i].eta(), noVertexTrackColl[i].phi());
       if(dR<minDr) {minDr=dR;
	 NVTrack=noVertexTrackColl[i];}
     }
     if(st) st->BS_dR_NVTrack->Fill(minDr,Event_Weight);
     if(minDr>0.4) return false;
     if(st)st->NVTrack->Fill(0.0,Event_Weight);

     //Find displacement of tracks with respect to beam spot
     fwlite::Handle<reco::BeamSpot> beamSpotCollHandle;
     beamSpotCollHandle.getByLabel(ev,"offlineBeamSpot");
     if(!beamSpotCollHandle.isValid()){printf("Beam Spot Collection NotFound\n");return false;}
     const reco::BeamSpot& beamSpotColl = *beamSpotCollHandle;

     dz  = NVTrack.dz (beamSpotColl.position());
     dxy = NVTrack.dxy(beamSpotColl.position());
     if(muonStations(NVTrack.hitPattern())<minMuStations) return false;
   }

   if(st){st->MTOF ->Fill(0.0,Event_Weight);
     if(GenBeta>=0)st->Beta_PreselectedB->Fill(GenBeta, Event_Weight);
   }

   double v3d = sqrt(dz*dz+dxy*dxy);

   if(st){st->BS_V3D->Fill(v3d,Event_Weight);}
   if(v3d>GlobalMaxV3D )return false;
   if(st){st->V3D  ->Fill(0.0,Event_Weight);}

   if(st)st->BS_Dxy->Fill(dxy, Event_Weight);

   TreeDXY = dxy;   
   bool DXYSB = false;
   if(TypeMode!=5 && fabs(dxy)>GlobalMaxDXY)return false;
   if(TypeMode==5 && fabs(dxy)>4)return false;
   if(TypeMode==5 && fabs(dxy)>GlobalMaxDXY) DXYSB = true;

   if(st){st->Dxy  ->Fill(0.0,Event_Weight);}

   if(TypeMode!=3) {
     fwlite::Handle<HSCPIsolationValueMap> IsolationH;
     IsolationH.getByLabel(ev, "HSCPIsolation", "R03"); //New format used for data since 17-07-2015
     if(!IsolationH.isValid()){
        IsolationH.getByLabel(ev, "HSCPIsolation03");//Old format used for first 2015B data, Signal and MC Backgrounds
        if(!IsolationH.isValid()){printf("Invalid IsolationH\n");return false;}
     }
     const ValueMap<HSCPIsolation>& IsolationMap = *IsolationH.product();

     HSCPIsolation hscpIso = IsolationMap.get((size_t)track.key());
     if(st){st->BS_TIsol ->Fill(hscpIso.Get_TK_SumEt(),Event_Weight);}
//     if(TypeMode!=4){       if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;     }
      if(hscpIso.Get_TK_SumEt()>GlobalMaxTIsol)return false;
     if(st){st->TIsol   ->Fill(0.0,Event_Weight);}

     double EoP = (hscpIso.Get_ECAL_Energy() + hscpIso.Get_HCAL_Energy())/track->p();
     if(st){st->BS_EIsol ->Fill(EoP,Event_Weight);}
//     if(TypeMode!=4){       if(EoP>GlobalMaxEIsol)return false;     }
     if(EoP>GlobalMaxEIsol)return false;
     if(st){st->EIsol   ->Fill(0.0,Event_Weight);}
     
     // relative tracker isolation
     if (st) {  st->BS_SumpTOverpT->Fill(hscpIso.Get_TK_SumEt()/track->pt(), Event_Weight); }
//     if(TypeMode==4) { if(hscpIso.Get_TK_SumEt()/track->pt()>GlobalMaxRelTIsol)return false;   }
     if(hscpIso.Get_TK_SumEt()/track->pt()>GlobalMaxRelTIsol)return false;
     if (st) {  st->SumpTOverpT   ->Fill(0.0,Event_Weight);} 
   }

   if(st){st->BS_Pterr ->Fill(track->ptError()/track->pt(),Event_Weight);}
   if(TypeMode!=3 && (track->ptError()/track->pt())>GlobalMaxPterr)return false;
   //mk if(MassErr > 0 && MassErr > 2.2)return false; //FIXME jozze -- cut on relative mass error in units of 8*MassErr/Mass

   if(std::max(0.0,track->pt())<GlobalMinPt)return false;
   if(st){st->Pterr   ->Fill(0.0,Event_Weight);}

   //Find distance to nearest segment on opposite side of detector
   double minPhi, minEta;
   double segSep=SegSep(hscp, ev, minPhi, minEta);

   if(st){
     st->BS_SegSep->Fill(segSep, Event_Weight);
     st->BS_SegMinPhiSep->Fill(minPhi, Event_Weight);
     st->BS_SegMinEtaSep->Fill(minEta, Event_Weight);
     //Plotting segment separation depending on whether track passed dz cut
     if(fabs(dz)>GlobalMaxDZ) {
       st->BS_SegMinEtaSep_FailDz->Fill(minEta, Event_Weight);
     }
     else {
       st->BS_SegMinEtaSep_PassDz->Fill(minEta, Event_Weight);
     }
     //Plots for tracking failing Eta Sep cut
     if(fabs(minEta)<minSegEtaSep) {
       //Needed to compare dz distribution of cosmics in pure cosmic and main sample
       st->BS_Dz_FailSep->Fill(dz);
     }
   }



   //Now cut Eta separation
   //if(TypeMode==3 && fabs(minEta)<minSegEtaSep) return false;
   if(st){st->SegSep->Fill(0.0,Event_Weight);}

   if(st) {
     //Plots for tracks in dz control region
     if(fabs(dz)>CosmicMinDz && fabs(dz)<CosmicMaxDz && !muon->isGlobalMuon()) {
       st->BS_Pt_FailDz->Fill(track->pt(), Event_Weight);
       st->BS_TOF_FailDz->Fill(tof->inverseBeta(), Event_Weight);
       if(fabs(track->eta())>CSCRegion) {
	 st->BS_TOF_FailDz_CSC->Fill(tof->inverseBeta(), Event_Weight);
	 st->BS_Pt_FailDz_CSC->Fill(track->pt(), Event_Weight);
       }
       else if(fabs(track->eta())<DTRegion) {
	 st->BS_TOF_FailDz_DT->Fill(tof->inverseBeta(), Event_Weight);
	 st->BS_Pt_FailDz_DT->Fill(track->pt(), Event_Weight);
       }
     }
     //Plots of dz
     st->BS_Dz->Fill(dz, Event_Weight);
     if(fabs(track->eta())>CSCRegion) st->BS_Dz_CSC->Fill(dz,Event_Weight);
     else if(fabs(track->eta())<DTRegion) st->BS_Dz_DT->Fill(dz,Event_Weight);
     st->BS_EtaDz->Fill(track->eta(),dz,Event_Weight);
   }


   //Split into different dz regions, each different region used to predict cosmic background and find systematic
   if(TypeMode==3 && !muon->isGlobalMuon() && st) {
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
       if(tof->inverseBeta()>=CutTOF[CutIndex]) {
	 st->H_D_DzSidebands->Fill(CutIndex, DzType);
       }
     }
   }

   TreeDZ = dz;
   bool DZSB = false;
   if(TypeMode!=5 && fabs(dz)>GlobalMaxDZ) return false;
   if(TypeMode==5 && fabs(dz)>4) return false;
   if(TypeMode==5 && fabs(dz)>GlobalMaxDZ) DZSB = true;
   if(st){st->Dz  ->Fill(0.0,Event_Weight);}

   if(TypeMode==3 && fabs(minEta)<minSegEtaSep) return false;
   if(st)st->BS_Phi->Fill(track->phi(),Event_Weight);
   if(TypeMode==3 && fabs(track->phi())>1.2 && fabs(track->phi())<1.9) return false;

    //skip HSCP that are compatible with cosmics.
    if(st)st->BS_OpenAngle->Fill(OpenAngle,Event_Weight);

    bool OASB = false;
    if(TypeMode==5 && OpenAngle>=2.8)OASB = true;

   isCosmicSB = DXYSB && DZSB && OASB;
   isSemiCosmicSB = (!isCosmicSB && (DXYSB || DZSB || OASB));
 
   if(st){if(dedxSObj) st->BS_EtaIs->Fill(track->eta(),dedxSObj->dEdx(),Event_Weight);
          if(dedxMObj) st->BS_EtaIm->Fill(track->eta(),dedxMObj->dEdx(),Event_Weight);
          st->BS_EtaP ->Fill(track->eta(),track->p(),Event_Weight);
          st->BS_EtaPt->Fill(track->eta(),track->pt(),Event_Weight);
          if(tof)st->BS_EtaTOF->Fill(track->eta(),tof->inverseBeta(),Event_Weight);
   }

   if(st){if(GenBeta>=0)st->Beta_PreselectedC->Fill(GenBeta, Event_Weight);
          if(DZSB  && OASB)st->BS_Dxy_Cosmic->Fill(dxy, Event_Weight);
          if(DXYSB && OASB)st->BS_Dz_Cosmic->Fill(dz, Event_Weight);
          if(DXYSB && DZSB)st->BS_OpenAngle_Cosmic->Fill(OpenAngle,Event_Weight);


          TVector3 outerHit = getOuterHitPos(dedxHits);
          TVector3 vertex(vertexColl[highestPtGoodVertex].position().x(), vertexColl[highestPtGoodVertex].position().y(), vertexColl[highestPtGoodVertex].position().z());
          st->BS_LastHitDXY  ->Fill((outerHit).Perp(),Event_Weight);
          st->BS_LastHitD3D  ->Fill((outerHit).Mag(),Event_Weight);

          st->BS_P  ->Fill(track->p(),Event_Weight);
          st->BS_Pt ->Fill(track->pt(),Event_Weight);
          if(PUA)st->BS_Pt_PUA ->Fill(track->pt(),Event_Weight);
          if(PUB)st->BS_Pt_PUB ->Fill(track->pt(),Event_Weight);
          if(DXYSB && DZSB && OASB) st->BS_Pt_Cosmic->Fill(track->pt(),Event_Weight);

	  if(fabs(track->eta())<DTRegion) st->BS_Pt_DT->Fill(track->pt(),Event_Weight);
	  else st->BS_Pt_CSC->Fill(track->pt(),Event_Weight);

          double RecoQoPt = track->charge()/track->pt();
          if(!hscp.trackRef().isNull() && hscp.trackRef()->pt()>200) {
            double InnerRecoQoPt = hscp.trackRef()->charge()/hscp.trackRef()->pt();
            st->BS_InnerInvPtDiff->Fill((RecoQoPt-InnerRecoQoPt)/InnerRecoQoPt,Event_Weight);
          }

          if(dedxSObj) st->BS_Is ->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj && PUA) st->BS_Is_PUA ->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj && PUB) st->BS_Is_PUB ->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj && DXYSB && DZSB && OASB) st->BS_Is_Cosmic->Fill(dedxSObj->dEdx(),Event_Weight);
          if(dedxSObj) st->BS_Im ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(dedxSObj && PUA) st->BS_Im_PUA ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(dedxSObj && PUB) st->BS_Im_PUB ->Fill(dedxMObj->dEdx(),Event_Weight);
          if(tof) {
	    st->BS_TOF->Fill(tof->inverseBeta(),Event_Weight);
            if(PUA)st->BS_TOF_PUA->Fill(tof->inverseBeta(),Event_Weight);
            if(PUB)st->BS_TOF_PUB->Fill(tof->inverseBeta(),Event_Weight);
	    if(dttof->nDof()>6) st->BS_TOF_DT->Fill(dttof->inverseBeta(),Event_Weight);
            if(csctof->nDof()>6) st->BS_TOF_CSC->Fill(csctof->inverseBeta(),Event_Weight);
            st->BS_PtTOF->Fill(track->pt() ,tof->inverseBeta(),Event_Weight);
	  }
          if(dedxSObj) {
	    st->BS_PIs  ->Fill(track->p()  ,dedxSObj->dEdx(),Event_Weight);
            st->BS_PImHD->Fill(track->p()  ,dedxMObj->dEdx(),Event_Weight);
            st->BS_PIm  ->Fill(track->p()  ,dedxMObj->dEdx(),Event_Weight);
            st->BS_PtIs ->Fill(track->pt() ,dedxSObj->dEdx(),Event_Weight);
            st->BS_PtIm ->Fill(track->pt() ,dedxMObj->dEdx(),Event_Weight);
	  }
          if(tof && dedxSObj)st->BS_TOFIs->Fill(tof->inverseBeta(),dedxSObj->dEdx(),Event_Weight);
          if(tof && dedxSObj)st->BS_TOFIm->Fill(tof->inverseBeta(),dedxMObj->dEdx(),Event_Weight);

	  //Muon only prediction binned depending on where in the detector the track is and how many muon stations it has
	  //Binning not used for other analyses
	  int bin=-1;
	  if(TypeMode==3) {
	    if(fabs(track->eta())<DTRegion) bin=muonStations(track->hitPattern())-2;
	    else bin=muonStations(track->hitPattern())+1;
	    st->BS_Pt_Binned[bin] ->Fill(track->pt(),Event_Weight);
	    if(tof) st->BS_TOF_Binned[bin]->Fill(tof->inverseBeta(),Event_Weight);
	  }
   }
   if(st){st->Basic  ->Fill(0.0,Event_Weight);}

   return true;
}

// check if one HSCP candidate is passing the selection (the function also has many more arguments because it is used to fill some histograms AND to evaluate the systematics
bool PassSelection(const susybsm::HSCParticle& hscp,  const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, const fwlite::ChainEvent& ev, const int& CutIndex, stPlots* st, const bool isFlip, const double& GenBeta, bool RescaleP, const double& RescaleI, const double& RescaleT){
   reco::TrackRef   track;
   if(TypeMode!=3) track = hscp.trackRef();
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
   double Ick=0; // if(dedxMObj) Ick=GetIck(Ih,isMC);

   double PtCut=CutPt[CutIndex];
   double ICut=CutI[CutIndex];
   double TOFCut=CutTOF[CutIndex];
   if(isFlip) {
     PtCut=CutPt_Flip[CutIndex];
     ICut=CutI_Flip[CutIndex];
     TOFCut=CutTOF_Flip[CutIndex];
   }

   if(RescaleP){
     if(RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())<PtCut)return false;
     //if(std::max(0.0,RescaledPt(track->pt() - track->ptError(),track->eta(),track->phi(),track->charge()))<CutPt[CutIndex])return false;
   }else{
     if(track->pt()<PtCut)return false;
     //if(std::max(0.0,(track->pt() - track->ptError()))<CutPt[CutIndex])return false;
   } 
   if(st){st->Pt    ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)st->Beta_SelectedP->Fill(CutIndex,GenBeta, Event_Weight);
   }

   if(TypeMode!=3 && Is+RescaleI<ICut)return false;

   if(st){st->I    ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)st->Beta_SelectedI->Fill(CutIndex, GenBeta, Event_Weight);
   }

   if((TypeMode>1  && TypeMode!=5) && !isFlip && MuonTOF+RescaleT<TOFCut)return false;
   if((TypeMode>1  && TypeMode!=5) && isFlip && MuonTOF+RescaleT>TOFCut)return false;

   if(st){st->TOF  ->Fill(CutIndex,Event_Weight);
          if(GenBeta>=0)st->Beta_SelectedT->Fill(CutIndex, GenBeta, Event_Weight);
          st->AS_P  ->Fill(CutIndex,track->p(),Event_Weight);
          st->AS_Pt ->Fill(CutIndex,track->pt(),Event_Weight);
          st->AS_Is ->Fill(CutIndex,Is,Event_Weight);
          st->AS_Im ->Fill(CutIndex,Ih,Event_Weight);
          st->AS_TOF->Fill(CutIndex,MuonTOF,Event_Weight);
//        st->AS_EtaIs->Fill(CutIndex,track->eta(),Is,Event_Weight);
//        st->AS_EtaIm->Fill(CutIndex,track->eta(),Ih,Event_Weight);
//        st->AS_EtaP ->Fill(CutIndex,track->eta(),track->p(),Event_Weight);
//        st->AS_EtaPt->Fill(CutIndex,track->eta(),track->pt(),Event_Weight);
          st->AS_PIs  ->Fill(CutIndex,track->p()  ,Is,Event_Weight);
          st->AS_PIm  ->Fill(CutIndex,track->p()  ,Ih,Event_Weight);
          st->AS_PtIs ->Fill(CutIndex,track->pt() ,Is,Event_Weight);
          st->AS_PtIm ->Fill(CutIndex,track->pt() ,Ih,Event_Weight);
          st->AS_TOFIs->Fill(CutIndex,MuonTOF     ,Is,Event_Weight);
          st->AS_TOFIm->Fill(CutIndex,MuonTOF     ,Ih,Event_Weight);
   }
   return true;
}


// all code for the filling of the ABCD related histograms --> this information will be used later in Step4 for the actual datadriven prediction
void Analysis_FillControlAndPredictionHist(const susybsm::HSCParticle& hscp, const reco::DeDxData* dedxSObj, const reco::DeDxData* dedxMObj, const reco::MuonTimeExtra* tof, stPlots* st){
	 reco::TrackRef   track;
         if(TypeMode!=3) track = hscp.trackRef();
         else {
	   reco::MuonRef muon = hscp.muonRef();
           if(muon.isNull()) return;
           track = muon->standAloneMuon();
         }

         double MuonTOF = GlobalMinTOF;
         if(tof){MuonTOF = tof->inverseBeta(); }

	 double Is=0; 	 if(dedxSObj) Is=dedxSObj->dEdx();
	 double Ih=0;	 if(dedxMObj) Ih=dedxMObj->dEdx();

         if(!isCosmicSB){
	 st->Hist_Pt->Fill(track->pt(),Event_Weight);
         st->Hist_Is->Fill(Is,Event_Weight);
         st->Hist_TOF->Fill(MuonTOF,Event_Weight);
         }

//          /\ I
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
	 std::vector<double> PtLimits;
	 if(TypeMode!=3) {
	   PtLimits.push_back(100);
           PtLimits.push_back(80);
           PtLimits.push_back(60);
	 }
	 else {
           PtLimits.push_back(240);
           PtLimits.push_back(170);
           PtLimits.push_back(120);
	 }

	    //Muon only prediction binned depending on where in the detector the track is and how many muon stations it has
	    //Binning not used for other analyses
	    int bin=-1;
	    if(TypeMode==3) {
	      if(fabs(track->eta())<DTRegion) bin=muonStations(track->hitPattern())-2;
	      else bin=muonStations(track->hitPattern())+1;
	    }

         if(!isCosmicSB){
            if(track->pt()>PtLimits[0]){
               st->CtrlPt_S4_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S4_Im->Fill(Ih, Event_Weight);
               if(tof)st->CtrlPt_S4_TOF->Fill(MuonTOF, Event_Weight);
               if(tof && bin>=0 && bin<MaxPredBins)st->CtrlPt_S4_TOF_Binned[bin]->Fill(MuonTOF, Event_Weight);
            }else if(track->pt()>PtLimits[1]){
               st->CtrlPt_S3_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S3_Im->Fill(Ih, Event_Weight);
               if(tof)st->CtrlPt_S3_TOF->Fill(MuonTOF, Event_Weight);
               if(tof && bin>=0 && bin<MaxPredBins)st->CtrlPt_S3_TOF_Binned[bin]->Fill(MuonTOF, Event_Weight);
            }else if(track->pt()>PtLimits[2]){
               st->CtrlPt_S2_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S2_Im->Fill(Ih, Event_Weight);
               if(tof)st->CtrlPt_S2_TOF->Fill(MuonTOF, Event_Weight);
               if(tof && bin>=0 && bin<MaxPredBins)st->CtrlPt_S2_TOF_Binned[bin]->Fill(MuonTOF, Event_Weight);
            }else{
               st->CtrlPt_S1_Is->Fill(Is, Event_Weight);
               st->CtrlPt_S1_Im->Fill(Ih, Event_Weight);
               if(tof)st->CtrlPt_S1_TOF->Fill(MuonTOF, Event_Weight);
               if(tof && bin>=0 && bin<MaxPredBins)st->CtrlPt_S1_TOF_Binned[bin]->Fill(MuonTOF, Event_Weight);
            }

            if(Is>0.2){           if(tof)st->CtrlIs_S4_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Is>0.1){     if(tof)st->CtrlIs_S3_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Is>0.05){    if(tof)st->CtrlIs_S2_TOF->Fill(MuonTOF, Event_Weight);
            }else{                if(tof)st->CtrlIs_S1_TOF->Fill(MuonTOF, Event_Weight);
            }

            if(Ih>4.4){           if(tof)st->CtrlIm_S4_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Ih>4.1){     if(tof)st->CtrlIm_S3_TOF->Fill(MuonTOF, Event_Weight);
            }else if(Ih>3.8){     if(tof)st->CtrlIm_S2_TOF->Fill(MuonTOF, Event_Weight);
            }else{                if(tof)st->CtrlIm_S1_TOF->Fill(MuonTOF, Event_Weight);
            }
         }



	 //	 if(dedxMObj) Ih=dedxMObj->dEdx();
	 double Ick=0;  if(dedxMObj) Ick=GetIck(Ih,isMCglobal);


         for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
 	    if(MuonTOF<GlobalMinTOF) continue;
            if(TypeMode==5 && isCosmicSB)continue;
            bool PassPtCut  = track->pt()>=CutPt[CutIndex];
            bool PassICut   = (Is>=CutI[CutIndex]);
            bool PassTOFCut = MuonTOF>=CutTOF[CutIndex];

            if(       PassTOFCut &&  PassPtCut &&  PassICut){   //Region D
               st->H_D      ->Fill(CutIndex,                Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_D_Binned[bin]->Fill(CutIndex,                Event_Weight);
               st->RegionD_P  ->Fill(CutIndex,track->p(),     Event_Weight);
               st->RegionD_I  ->Fill(CutIndex,Ih,Event_Weight);
	       st->RegionD_Ias->Fill(CutIndex,Is,Event_Weight);
               st->RegionD_TOF->Fill(CutIndex,MuonTOF,        Event_Weight);
	       st->AS_Eta_RegionD->Fill(CutIndex,track->eta());
            }else if( PassTOFCut &&  PassPtCut && !PassICut){   //Region C
               st->H_C     ->Fill(CutIndex,                 Event_Weight);
               if(TypeMode<2)st->Pred_EtaP  ->Fill(CutIndex,track->eta(), track->p(),     Event_Weight);
               st->PDF_C_EtaP ->Fill(CutIndex,track->eta(), track->p(),     Event_Weight); //pz
               //Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
               st->AS_Eta_RegionC->Fill(CutIndex,track->eta());
            }else if( PassTOFCut && !PassPtCut &&  PassICut){   //Region B
               st->H_B     ->Fill(CutIndex,                 Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_B_Binned[bin]->Fill(CutIndex,                Event_Weight);
               if(TypeMode<2)st->Pred_I  ->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode<2)st->Pred_EtaS->Fill(CutIndex,track->eta(),         Event_Weight);
	       st->PDF_B_EtaICK ->Fill(CutIndex,track->eta(),Ick, Event_Weight); //pz
               //Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
               st->AS_Eta_RegionB->Fill(CutIndex,track->eta());
            }else if( PassTOFCut && !PassPtCut && !PassICut){   //Region A
               st->H_A     ->Fill(CutIndex,                 Event_Weight);
               if(TypeMode==2)st->Pred_TOF->Fill(CutIndex,MuonTOF,         Event_Weight);
               if(TypeMode<2)st->Pred_EtaB->Fill(CutIndex,track->eta(),         Event_Weight);
               if(TypeMode==2)st->Pred_EtaS2->Fill(CutIndex,track->eta(),        Event_Weight);
               st->AS_Eta_RegionA->Fill(CutIndex,track->eta());
	       st->PDF_A_Eta->Fill(CutIndex,track->eta(),        Event_Weight);//pz

            }else if(!PassTOFCut &&  PassPtCut &&  PassICut){   //Region H
               st->H_H   ->Fill(CutIndex,          Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_H_Binned[bin]->Fill(CutIndex,                Event_Weight);
	       st->RegionH_Ias->Fill(CutIndex,Is,Event_Weight);
	       if(TypeMode==2 && Ick>0)st->PDF_H_EtaMass ->Fill(CutIndex,track->eta(),track->p()*sqrt(Ick), Event_Weight); //pz
               //Pred_P->Fill(CutIndex,track->p(),        Event_Weight);
               //Pred_I->Fill(CutIndex,Ih,   Event_Weight);
               if(TypeMode==2)st->AS_Eta_RegionH->Fill(CutIndex,track->eta());
            }else if(!PassTOFCut &&  PassPtCut && !PassICut){   //Region G
               st->H_G     ->Fill(CutIndex,                 Event_Weight);
               if(TypeMode==2)st->Pred_EtaP  ->Fill(CutIndex,track->eta(),track->p(),     Event_Weight);
               st->AS_Eta_RegionG->Fill(CutIndex,track->eta());
               if(TypeMode==2)st->PDF_G_EtaP ->Fill(CutIndex,track->eta(), track->p(),     Event_Weight); //pz
            }else if(!PassTOFCut && !PassPtCut &&  PassICut){   //Region F
               st->H_F     ->Fill(CutIndex,                 Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_F_Binned[bin]->Fill(CutIndex,                Event_Weight);
               if(TypeMode==2)st->Pred_I  ->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode==2)st->Pred_EtaS->Fill(CutIndex,track->eta(),         Event_Weight);
               st->AS_Eta_RegionF->Fill(CutIndex,track->eta());
	       if(TypeMode==2)st->PDF_F_EtaICK ->Fill(CutIndex,track->eta(),Ick, Event_Weight); //pz

            }else if(!PassTOFCut && !PassPtCut && !PassICut){   //Region E
               st->H_E     ->Fill(CutIndex,                 Event_Weight);
               if(TypeMode==2)st->Pred_EtaB->Fill(CutIndex,track->eta(),         Event_Weight);
               st->AS_Eta_RegionE->Fill(CutIndex,track->eta());
	       if(TypeMode==2)st->PDF_E_Eta->Fill(CutIndex,track->eta(),        Event_Weight);//pz

            }
         }

	 //Use events with low TOF to check accuracy of background prediction
         for(unsigned int CutIndex=0;CutIndex<CutPt_Flip.size();CutIndex++){
            if(TypeMode!=5 && MuonTOF>=GlobalMinTOF) continue;
            if(TypeMode==5 && !isCosmicSB)continue;

            bool PassPtCut  = track->pt()>=CutPt_Flip[CutIndex];
            bool PassICut   = (Is>=CutI_Flip[CutIndex]);
            bool PassTOFCut = MuonTOF<=CutTOF_Flip[CutIndex]; 


            if(TypeMode==5)PassTOFCut=true;

            if(       PassTOFCut &&  PassPtCut &&  PassICut){   //Region D
	      st->RegionD_P_Flip  ->Fill(CutIndex,track->p(),     Event_Weight);
	      st->RegionD_I_Flip  ->Fill(CutIndex,Ih,Event_Weight);
	      st->RegionD_Ias_Flip  ->Fill(CutIndex,Is,Event_Weight);
	      st->RegionD_TOF_Flip->Fill(CutIndex,MuonTOF,        Event_Weight);
               st->H_D_Flip->Fill(CutIndex,                Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_D_Binned_Flip[bin]->Fill(CutIndex,                Event_Weight);
            }else if( PassTOFCut &&  PassPtCut && !PassICut){   //Region C
               st->H_C_Flip->Fill(CutIndex,                 Event_Weight);
               if(TypeMode<2)st->Pred_EtaP_Flip->Fill(CutIndex,track->eta(), track->p(),     Event_Weight);
               st->PDF_C_EtaP_Flip ->Fill(CutIndex,track->eta(), track->p(),     Event_Weight); //pz
               //Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
            }else if( PassTOFCut && !PassPtCut &&  PassICut){   //Region B
               st->H_B_Flip->Fill(CutIndex,                 Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_B_Binned_Flip[bin]->Fill(CutIndex,                Event_Weight);
               if(TypeMode<2)st->Pred_I_Flip->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode<2)st->Pred_EtaS_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
	       st->PDF_B_EtaICK_Flip ->Fill(CutIndex,track->eta(),Ick, Event_Weight); //pz
               //Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
            }else if( PassTOFCut && !PassPtCut && !PassICut){   //Region A
               st->H_A_Flip->Fill(CutIndex,                 Event_Weight);
               if(TypeMode==2)st->Pred_TOF_Flip->Fill(CutIndex,MuonTOF,         Event_Weight);
               if(TypeMode<2)st->Pred_EtaB_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
               if(TypeMode==2)st->Pred_EtaS2_Flip->Fill(CutIndex,track->eta(),        Event_Weight);
	       st->PDF_A_Eta_Flip->Fill(CutIndex,track->eta(),        Event_Weight);//pz
            }else if(!PassTOFCut &&  PassPtCut &&  PassICut){   //Region H
               st->H_H_Flip->Fill(CutIndex,          Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_H_Binned_Flip[bin]->Fill(CutIndex,                Event_Weight);
	       st->RegionH_Ias_Flip  ->Fill(CutIndex,Is,Event_Weight);
	       if(TypeMode==2 && Ick>0)st->PDF_H_EtaMass_Flip ->Fill(CutIndex,track->eta(),track->p()*sqrt(Ick), Event_Weight); //pz

	       //Pred_P_Flip->Fill(CutIndex,track->p(),        Event_Weight);
	       //Pred_I_Flip->Fill(CutIndex,Ih,   Event_Weight);
            }else if(!PassTOFCut &&  PassPtCut && !PassICut){   //Region G
               st->H_G_Flip->Fill(CutIndex,                 Event_Weight);
               if(TypeMode==2)st->Pred_EtaP_Flip->Fill(CutIndex,track->eta(),track->p(),     Event_Weight);
               if(TypeMode==2)st->PDF_G_EtaP_Flip ->Fill(CutIndex,track->eta(), track->p(),     Event_Weight); //pz

            }else if(!PassTOFCut && !PassPtCut &&  PassICut){   //Region F
               st->H_F_Flip->Fill(CutIndex,                 Event_Weight);
               if(bin>-1 && bin<MaxPredBins) st->H_F_Binned_Flip[bin]->Fill(CutIndex,                Event_Weight);
               if(TypeMode==2)st->Pred_I_Flip->Fill(CutIndex,Ih, Event_Weight);
               if(TypeMode==2)st->Pred_EtaS_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
	       if(TypeMode==2)st->PDF_F_EtaICK_Flip ->Fill(CutIndex,track->eta(),Ick, Event_Weight); //pz
            }else if(!PassTOFCut && !PassPtCut && !PassICut){   //Region E
               st->H_E_Flip->Fill(CutIndex,                 Event_Weight);
               if(TypeMode==2)st->Pred_EtaB_Flip->Fill(CutIndex,track->eta(),         Event_Weight);
	       if(TypeMode==2)st->PDF_E_Eta_Flip->Fill(CutIndex,track->eta(),        Event_Weight);//pz
            }
         }
}


// Looping on all events, tracks, selection and check how many events are entering the mass distribution
void Analysis_Step1_EventLoop(char* SavePath)
{
   //Initialize a RandomNumberGenerator
   TRandom3* RNG = new TRandom3();

   //Initialize histo common to all samples
   InitHistos(NULL);

   for(unsigned int s=0;s<samples.size();s++){
      bool isData   = (samples[s].Type==0);
      bool isMC     = (samples[s].Type==1);
      bool isSignal = (samples[s].Type>=2);
      bool is2016   = (samples[s].Name.find("13TeV16")==std::string::npos)?false:true;
      bool is2016G  = (samples[s].Name.find("13TeV16G")==std::string::npos)?false:true;
      isMCglobal = isMC;

      dEdxK_Data = is2016?dEdxK_Data16:dEdxK_Data15;
      dEdxC_Data = is2016?dEdxC_Data16:dEdxC_Data15;
      dEdxK_MC   = is2016?dEdxK_MC16:dEdxK_MC15;
      dEdxC_MC   = is2016?dEdxC_MC16:dEdxC_MC15;

std::cout<<"D\n";

      char basepath [200]; sprintf (basepath, "%s/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/", getenv("CMSSW_BASE"));
      string analysis_path (basepath);
      if(isData){ 
         dEdxSF [0] = 1.00000;
	 //mk_dEdxSF [1] = 1.6107 *0.91345;  //PreG
         //mk_ if (is2016G) dEdxSF[1] = 1.6107 * 1.06665;  //PostG - first period  -- change for the other two periods later
         dEdxTemplates = loadDeDxTemplate(analysis_path+"../../data/Data13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true); //fix if you want to run on 2015
      }else{  
         dEdxSF [0] = 1.09711;
         dEdxSF [1] = 1.09256;
         dEdxTemplates = loadDeDxTemplate(analysis_path+"../../data/MC13TeV_Deco_SiStripDeDxMip_3D_Rcd_v2_CCwCI.root", true); // this will need to be checked if we rerun on MC
      }

std::cout<<"E\n";

 if(isData){    trackerCorrector.LoadDeDxCalibration(analysis_path+"../../data/Data13TeVGains_v2.root");  //Je: those are the correct calib tree to use. I don't know why the default was null, but it was wrong
	// if(isData){    trackerCorrector.TrackerGains = NULL;
      }else{ trackerCorrector.TrackerGains = NULL; //FIXME check gain for MC
      }

std::cout<<"F\n";

      //check that the plot container exist for this sample, otherwise create it
      if(plotsMap.find(samples[s].Name)==plotsMap.end()){plotsMap[samples[s].Name] = stPlots();}
      //For data and MCTr only initialize prediction histograms
      if(isData) stPlots_Init(HistoFile,plotsMap[samples[s].Name],samples[s].Name, CutPt.size(), false, false, CutPt_Flip.size());
      else stPlots_Init(HistoFile,plotsMap[samples[s].Name],samples[s].Name, CutPt.size());
      stPlots* SamplePlots = &plotsMap[samples[s].Name];
      if (is2016G)                 SamplePlots->IntLumi->Fill(0.0,IntegratedLuminosity13TeV16G);
      else if (!is2016G && is2016) SamplePlots->IntLumi->Fill(0.0,IntegratedLuminosity13TeV16PreG);
      else                         SamplePlots->IntLumi->Fill(0.0,IntegratedLuminosity13TeV15);

      string MCTrDirName = "MCTr_13TeV";
      if(isMC){
         if(samples[s].Name.find("7TeV")!=string::npos)          MCTrDirName = "MCTr_7TeV";
	 else if(samples[s].Name.find("8TeV")!=string::npos)     MCTrDirName = "MCTr_8TeV";
	 else if(samples[s].Name.find("13TeV16G")!=string::npos) MCTrDirName = "MCTr_13TeV16G";
	 else if(samples[s].Name.find("13TeV16")!=string::npos)  MCTrDirName = "MCTr_13TeV16";
	 else if(plotsMap.find(MCTrDirName)==plotsMap.end()){plotsMap[MCTrDirName] = stPlots();}
         stPlots_Init(HistoFile,plotsMap[MCTrDirName],MCTrDirName, CutPt.size(), false, false, CutPt_Flip.size());
      }stPlots* MCTrPlots = &plotsMap[MCTrDirName];

      //Initialize plot container for pure cosmic sample
      //Cosmic sample is contained in data file so for TOF-Only search
      //need a new set of plots for these events
      string CosmicName="";
      if(isData && TypeMode==3) {
         if(samples[s].Name.find("8TeV")!=string::npos) CosmicName="Cosmic8TeV";
         else CosmicName="Cosmic7TeV";
         if(plotsMap.find(CosmicName)==plotsMap.end()){plotsMap[CosmicName] = stPlots();}
         stPlots_Init(HistoFile,plotsMap[CosmicName],CosmicName, CutPt.size(), false, false, CutPt_Flip.size());
      }

      //needed for bookeeping
      bool* HSCPTk              = new bool[CutPt.size()];
      bool* HSCPTk_SystP        = new bool[CutPt.size()];
      bool* HSCPTk_SystI        = new bool[CutPt.size()];
      bool* HSCPTk_SystT        = new bool[CutPt.size()];
      bool* HSCPTk_SystM        = new bool[CutPt.size()];
      bool* HSCPTk_SystPU       = new bool[CutPt.size()];
      bool* HSCPTk_SystHUp      = new bool[CutPt.size()];
      bool* HSCPTk_SystHDown    = new bool[CutPt.size()];
      double* MaxMass           = new double[CutPt.size()];
      double* MaxMass_SystP     = new double[CutPt.size()];
      double* MaxMass_SystI     = new double[CutPt.size()];
      double* MaxMass_SystT     = new double[CutPt.size()];
      double* MaxMass_SystM     = new double[CutPt.size()];
      double* MaxMass_SystPU    = new double[CutPt.size()];
      double* MaxMass_SystHUp   = new double[CutPt.size()];
      double* MaxMass_SystHDown = new double[CutPt.size()];

      moduleGeom::loadGeometry(analysis_path+"../../data/CMS_GeomTree.root");
      muonTimingCalculator tofCalculator;
      tofCalculator.loadTimeOffset(analysis_path+"../../data/MuonTimeOffset.txt");
      unsigned int CurrentRun = 0;


      //do two loops through signal for samples with and without trigger changes.
      for (int period=0; period<(samples[s].Type>=2?RunningPeriods:1); period++){
         //load the files corresponding to this sample
         std::vector<string> FileName;
	 GetInputFiles(samples[s], BaseDirectory, FileName, period);
         fwlite::ChainEvent ev(FileName);

         DuplicatesClass duplicateChecker; 
         duplicateChecker.Clear();
         bool checkDuplicates = isData && FileName.size()>1;
	 checkDuplicates = true; // FIXME JOZE
         if(checkDuplicates){printf("Duplicated events will be removed\n");}

         //compute sample global weight
         Event_Weight = 1.0;
         double SampleWeight = 1.0;
         double PUSystFactor;

         if(samples[s].Type>0){           
            //get PU reweighted total # MC events.
            double NMCevents=0;
            for(Long64_t ientry=0;ientry<ev.size();ientry++){
              ev.to(ientry);
              if(MaxEntry>0 && ientry>MaxEntry)break;
              NMCevents += GetPUWeight(ev, samples[s].Pileup, PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
            }
            if(samples[s].Type==1){
              if      (is2016G) SampleWeight = GetSampleWeightMC (IntegratedLuminosity13TeV16G,    FileName, samples[s].XSec, ev.size(), NMCevents, numberOfMatchingSamples(samples[s].Name, samplesFull));
	      else if (is2016)  SampleWeight = GetSampleWeightMC (IntegratedLuminosity13TeV16PreG, FileName, samples[s].XSec, ev.size(), NMCevents, numberOfMatchingSamples(samples[s].Name, samplesFull));
	      else              SampleWeight = GetSampleWeightMC (IntegratedLuminosity13TeV15,     FileName, samples[s].XSec, ev.size(), NMCevents, numberOfMatchingSamples(samples[s].Name, samplesFull));
	    }
            else {
              if      (is2016G) SampleWeight = GetSampleWeight (IntegratedLuminosity13TeV16G,    IntegratedLuminosityBeforeTriggerChange,samples[s].XSec,NMCevents, period);
	      else if (is2016)  SampleWeight = GetSampleWeight (IntegratedLuminosity13TeV16PreG, IntegratedLuminosityBeforeTriggerChange,samples[s].XSec,NMCevents, period);
	      else              SampleWeight = GetSampleWeight (IntegratedLuminosity13TeV15,     IntegratedLuminosityBeforeTriggerChange,samples[s].XSec,NMCevents, period);
	    }
         }

	 if(SampleWeight==0) continue; //If sample weight 0 don't run, happens Int Lumi before change = 0

std::cout<<"G\n";

         //Loop on the events
         printf("Progressing Bar                   :0%%       20%%       40%%       60%%       80%%       100%%\n");
         printf("Building Mass for %10s (%1i/%1i) :",samples[s].Name.c_str(),period+1,(samples[s].Type>=2?RunningPeriods:1));
         int TreeStep = ev.size()/50;if(TreeStep==0)TreeStep=1;


	 //for(Long64_t ientry=0;ientry<20;ientry++){
		    //mkpz 
	   for(Long64_t ientry=0;ientry<ev.size();ientry++){
            ev.to(ientry);
            if(MaxEntry>0 && ientry>MaxEntry)break;
            if(ientry%TreeStep==0){printf(".");fflush(stdout);}
            if(checkDuplicates && duplicateChecker.isDuplicate(ev.eventAuxiliary().run(), ev.eventAuxiliary().event()))continue;

	    //Je: once we have the templates we have to update those lines below as well.. 
            //if run change, update conditions
            if(CurrentRun != ev.eventAuxiliary().run()){
               CurrentRun = ev.eventAuxiliary().run();
               tofCalculator.setRun(CurrentRun);
               trackerCorrector.setRun(CurrentRun);
              
               if(isData){//fix if you want to use 2015
                 //preG
               dEdxSF [0] = 1.00000;
               dEdxSF [1] = 1.6107*0.91345;
               dEdxK_Data = 2.062;
               dEdxC_Data = 3.430;

               if (278769 <= CurrentRun && CurrentRun <= 278822){
                dEdxSF [0] = 1.00000;
                dEdxSF [1] = 1.6107*1.06665;
                dEdxK_Data = 2.300;
                dEdxC_Data = 3.825;
               }

               if (278873 <= CurrentRun && CurrentRun <= 279116){
                dEdxSF [0] = 1.00000;
                dEdxSF [1] = 1.6107*1.07695;
                dEdxK_Data = 2.300;
                dEdxC_Data = 3.799;

               }

               if (279479 <= CurrentRun){
                dEdxSF [0] = 1.00000;
                dEdxSF [1] = 1.6107*1.0448500;
                dEdxK_Data = 2.275;
                dEdxC_Data = 3.675;
               }
              } 
              if(isMC){
//                  dEdxSF [0] = 1.09711;
//                  dEdxSF [1] = 1.16;
//old Joze parameters to be used
         dEdxSF [0] = 1.09711;
         dEdxSF [1] = 1.09256;
//KC takie jest od Joze z Global


         dEdxK_MC   = 2.935;
         dEdxC_MC   = 3.197;
                  }


              std::cout<<"\n------> dEdx parameters SF for run = "<<CurrentRun<< "  "<< dEdxSF[1]<<std::endl;
          } 



            //compute event weight
            if(samples[s].Type>0){Event_Weight = SampleWeight * GetPUWeight(ev, samples[s].Pileup, PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);}else{Event_Weight = 1;}
            std::vector<reco::GenParticle> genColl;
            double HSCPGenBeta1=-1, HSCPGenBeta2=-1;
            double HSCPDLength1=-1, HSCPDLength2=-1;
            if(isSignal){
               //get the collection of generated Particles
               fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
               genCollHandle.getByLabel(ev, "genParticlePlusGeant");
               if(!genCollHandle.isValid()){
                  genCollHandle.getByLabel(ev, "genParticlesSkimmed");
                  if(!genCollHandle.isValid()){
                     genCollHandle.getByLabel(ev, "genParticles");
                     if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound BUG\n");continue;}
                  }
               }

               genColl = *genCollHandle;
               int NChargedHSCP=HowManyChargedHSCP(genColl);
               Event_Weight*=samples[s].GetFGluinoWeight(NChargedHSCP);

               GetGenHSCPDecayLength(genColl,HSCPDLength1,HSCPDLength2,true);
               SamplePlots->Gen_DecayLength->Fill(HSCPDLength1, Event_Weight); SamplePlots->Gen_DecayLength->Fill(HSCPDLength2, Event_Weight);

               GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,false);
               if(HSCPGenBeta1>=0)SamplePlots->Beta_Gen      ->Fill(HSCPGenBeta1, Event_Weight);  if(HSCPGenBeta2>=0)SamplePlots->Beta_Gen       ->Fill(HSCPGenBeta2, Event_Weight);
               GetGenHSCPBeta(genColl,HSCPGenBeta1,HSCPGenBeta2,true);
               if(HSCPGenBeta1>=0)SamplePlots->Beta_GenCharged->Fill(HSCPGenBeta1, Event_Weight); if(HSCPGenBeta2>=0)SamplePlots->Beta_GenCharged->Fill(HSCPGenBeta2, Event_Weight);

// R-hadron wights needed due to wrong GenId---------------------------------BEGIN

   double  Wa=1.0, Wad=1.0, Waa=1.0, Wan=1.0;  // Wa is additional weight for single other, Wad for other+double_charged, 
                                              // Waa for the event with 2 other R-hadron, Wan for other+neutral
   bool Rhadron=0; // default value - not R-hadron (not need to weight)
     
   string sample_name = samples[s].Name;

   //------------ weights for R-hadron samples ----------Start

   if(sample_name== "Gluino_13TeV16_M100N_f10"){Rhadron=1; Wa=1.02544e+00; Wad=1.07061e+00; Waa=1.13625e+00; Wan=1.06637e+00;}
   if(sample_name== "Gluino_13TeV16G_M100N_f10"){Rhadron=1; Wa=1.01333e+00; Wad=1.05776e+00; Waa=1.11787e+00; Wan=1.06393e+00;}
   if(sample_name== "Gluino_13TeV16_M100N_f50"){Rhadron=1; Wa=1.02544e+00; Wad=1.07061e+00; Waa=1.13625e+00; Wan=1.06637e+00;}
   if(sample_name== "Gluino_13TeV16G_M100N_f50"){Rhadron=1; Wa=1.01333e+00; Wad=1.05776e+00; Waa=1.11788e+00; Wan=1.06393e+00;}
   if(sample_name== "Gluino_13TeV16_M200_f10"){Rhadron=1; Wa=1.01060e+00; Wad=1.08585e+00; Waa=1.14509e+00; Wan=1.07439e+00;}
   if(sample_name== "Gluino_13TeV16G_M200_f10"){Rhadron=1; Wa=1.00530e+00; Wad=1.08020e+00; Waa=1.13846e+00; Wan=1.07613e+00;}
   if(sample_name== "Gluino_13TeV16_M200N_f10"){Rhadron=1; Wa=1.00138e+00; Wad=1.05406e+00; Waa=1.13549e+00; Wan=1.06934e+00;}
   if(sample_name== "Gluino_13TeV16G_M200N_f10"){Rhadron=1; Wa=1.00012e+00; Wad=1.04085e+00; Waa=1.12277e+00; Wan=1.07311e+00;}
   if(sample_name== "Gluino_13TeV16_M200_f50"){Rhadron=1; Wa=1.01060e+00; Wad=1.08585e+00; Waa=1.14510e+00; Wan=1.07439e+00;}
   if(sample_name== "Gluino_13TeV16G_M200_f50"){Rhadron=1; Wa=1.00530e+00; Wad=1.08020e+00; Waa=1.13846e+00; Wan=1.07612e+00;}
   if(sample_name== "Gluino_13TeV16_M200N_f50"){Rhadron=1; Wa=1.00138e+00; Wad=1.05406e+00; Waa=1.13549e+00; Wan=1.06935e+00;}
   if(sample_name== "Gluino_13TeV16G_M200N_f50"){Rhadron=1; Wa=1.00012e+00; Wad=1.04085e+00; Waa=1.12277e+00; Wan=1.07311e+00;}
   if(sample_name== "Gluino_13TeV16_M400_f10"){Rhadron=1; Wa=1.01706e+00; Wad=1.11794e+00; Waa=1.14045e+00; Wan=1.06800e+00;}
   if(sample_name== "Gluino_13TeV16G_M400_f10"){Rhadron=1; Wa=1.00483e+00; Wad=1.24734e+00; Waa=1.12817e+00; Wan=1.07278e+00;}
   if(sample_name== "Gluino_13TeV16_M400N_f10"){Rhadron=1; Wa=1.08215e+00; Wad=1.06579e+00; Waa=1.14654e+00; Wan=1.05864e+00;}
   if(sample_name== "Gluino_13TeV16G_M400N_f10"){Rhadron=1; Wa=1.06555e+00; Wad=1.06825e+00; Waa=1.15421e+00; Wan=1.05185e+00;}
   if(sample_name== "Gluino_13TeV16_M400_f50"){Rhadron=1; Wa=1.01706e+00; Wad=1.11794e+00; Waa=1.14045e+00; Wan=1.06800e+00;}
   if(sample_name== "Gluino_13TeV16G_M400_f50"){Rhadron=1; Wa=1.00483e+00; Wad=1.24734e+00; Waa=1.12816e+00; Wan=1.07277e+00;}
   if(sample_name== "Gluino_13TeV16_M400N_f50"){Rhadron=1; Wa=1.08215e+00; Wad=1.06579e+00; Waa=1.14654e+00; Wan=1.05864e+00;}
   if(sample_name== "Gluino_13TeV16G_M400N_f50"){Rhadron=1; Wa=1.06555e+00; Wad=1.06825e+00; Waa=1.15421e+00; Wan=1.05185e+00;}
   if(sample_name== "Gluino_13TeV16_M600_f10"){Rhadron=1; Wa=1.07582e+00; Wad=1.05595e+00; Waa=1.14064e+00; Wan=1.06812e+00;}
   if(sample_name== "Gluino_13TeV16G_M600_f10"){Rhadron=1; Wa=1.05737e+00; Wad=1.04830e+00; Waa=1.13847e+00; Wan=1.07254e+00;}
   if(sample_name== "Gluino_13TeV16_M600N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.06769e+00; Waa=1.12687e+00; Wan=1.05991e+00;}
   if(sample_name== "Gluino_13TeV16G_M600N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.04890e+00; Waa=1.11368e+00; Wan=1.05742e+00;}
   if(sample_name== "Gluino_13TeV16_M600_f50"){Rhadron=1; Wa=1.07582e+00; Wad=1.05595e+00; Waa=1.14065e+00; Wan=1.06813e+00;}
   if(sample_name== "Gluino_13TeV16G_M600_f50"){Rhadron=1; Wa=1.05737e+00; Wad=1.04830e+00; Waa=1.13845e+00; Wan=1.07254e+00;}
   if(sample_name== "Gluino_13TeV16_M600N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.06769e+00; Waa=1.12687e+00; Wan=1.05991e+00;}
   if(sample_name== "Gluino_13TeV16G_M600N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.04890e+00; Waa=1.11368e+00; Wan=1.05742e+00;}
   if(sample_name== "Gluino_13TeV16_M800_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.03583e+00; Waa=1.13499e+00; Wan=1.06721e+00;}
   if(sample_name== "Gluino_13TeV16G_M800_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.02294e+00; Waa=1.13032e+00; Wan=1.07048e+00;}
   if(sample_name== "Gluino_13TeV16_M800N_f10"){Rhadron=1; Wa=1.17722e+00; Wad=1.03942e+00; Waa=1.13085e+00; Wan=1.06297e+00;}
   if(sample_name== "Gluino_13TeV16G_M800N_f10"){Rhadron=1; Wa=1.18090e+00; Wad=1.02265e+00; Waa=1.11868e+00; Wan=1.06410e+00;}
   if(sample_name== "Gluino_13TeV16_M800_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.03583e+00; Waa=1.13499e+00; Wan=1.06721e+00;}
   if(sample_name== "Gluino_13TeV16G_M800_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.02294e+00; Waa=1.13031e+00; Wan=1.07049e+00;}
   if(sample_name== "Gluino_13TeV16_M800N_f50"){Rhadron=1; Wa=1.17722e+00; Wad=1.03942e+00; Waa=1.13085e+00; Wan=1.06296e+00;}
   if(sample_name== "Gluino_13TeV16G_M800N_f50"){Rhadron=1; Wa=1.18090e+00; Wad=1.02265e+00; Waa=1.11868e+00; Wan=1.06410e+00;}
   if(sample_name== "Gluino_13TeV16_M1000_f10"){Rhadron=1; Wa=1.09463e+01; Wad=1.03897e+00; Waa=1.12760e+00; Wan=1.06332e+00;}
   if(sample_name== "Gluino_13TeV16G_M1000_f10"){Rhadron=1; Wa=3.78760e+01; Wad=1.03466e+00; Waa=1.12804e+00; Wan=1.06625e+00;}
   if(sample_name== "Gluino_13TeV16_M1000N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.07928e+00; Waa=1.13146e+00; Wan=1.05836e+00;}
   if(sample_name== "Gluino_13TeV16G_M1000N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.08561e+00; Waa=1.12628e+00; Wan=1.05502e+00;}
   if(sample_name== "Gluino_13TeV16_M1000_f50"){Rhadron=1; Wa=1.09463e+01; Wad=1.03897e+00; Waa=1.12760e+00; Wan=1.06331e+00;}
   if(sample_name== "Gluino_13TeV16G_M1000_f50"){Rhadron=1; Wa=3.78760e+01; Wad=1.03466e+00; Waa=1.12803e+00; Wan=1.06624e+00;}
   if(sample_name== "Gluino_13TeV16_M1000N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.07928e+00; Waa=1.13146e+00; Wan=1.05836e+00;}
   if(sample_name== "Gluino_13TeV16G_M1000N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.08561e+00; Waa=1.12627e+00; Wan=1.05502e+00;}
   if(sample_name== "Gluino_13TeV16_M1200_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.04179e+00; Waa=1.13407e+00; Wan=1.06200e+00;}
   if(sample_name== "Gluino_13TeV16G_M1200_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.02489e+00; Waa=1.12876e+00; Wan=1.06085e+00;}
   if(sample_name== "Gluino_13TeV16_M1200N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.08557e+00; Waa=1.14153e+00; Wan=1.05728e+00;}
   if(sample_name== "Gluino_13TeV16G_M1200N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.08035e+00; Waa=1.14312e+00; Wan=1.05405e+00;}
   if(sample_name== "Gluino_13TeV16_M1200_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.04179e+00; Waa=1.13407e+00; Wan=1.06199e+00;}
   if(sample_name== "Gluino_13TeV16G_M1200_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.02489e+00; Waa=1.12877e+00; Wan=1.06085e+00;}
   if(sample_name== "Gluino_13TeV16_M1200N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.08557e+00; Waa=1.14154e+00; Wan=1.05728e+00;}
   if(sample_name== "Gluino_13TeV16G_M1200N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.08035e+00; Waa=1.14313e+00; Wan=1.05405e+00;}
   if(sample_name== "Gluino_13TeV16_M1400_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.06238e+00; Waa=1.12582e+00; Wan=1.06014e+00;}
   if(sample_name== "Gluino_13TeV16G_M1400_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.04640e+00; Waa=1.12155e+00; Wan=1.06117e+00;}
   if(sample_name== "Gluino_13TeV16_M1400N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.06739e+00; Waa=1.13203e+00; Wan=1.05467e+00;}
   if(sample_name== "Gluino_13TeV16G_M1400N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.06340e+00; Waa=1.14317e+00; Wan=1.05446e+00;}
   if(sample_name== "Gluino_13TeV16_M1400_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.06238e+00; Waa=1.12582e+00; Wan=1.06014e+00;}
   if(sample_name== "Gluino_13TeV16G_M1400_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.04640e+00; Waa=1.12154e+00; Wan=1.06117e+00;}
   if(sample_name== "Gluino_13TeV16_M1400N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.06739e+00; Waa=1.13204e+00; Wan=1.05466e+00;}
   if(sample_name== "Gluino_13TeV16G_M1400N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.06340e+00; Waa=1.14317e+00; Wan=1.05446e+00;}
   if(sample_name== "Gluino_13TeV16_M1600_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.13721e+00; Waa=1.12260e+00; Wan=1.05958e+00;}
   if(sample_name== "Gluino_13TeV16G_M1600_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.15968e+00; Waa=1.11867e+00; Wan=1.05941e+00;}
   if(sample_name== "Gluino_13TeV16_M1600N_f10"){Rhadron=1; Wa=1.40489e+00; Wad=1.05869e+00; Waa=1.12437e+00; Wan=1.05206e+00;}
   if(sample_name== "Gluino_13TeV16G_M1600N_f10"){Rhadron=1; Wa=1.36716e+00; Wad=1.03935e+00; Waa=1.12433e+00; Wan=1.04814e+00;}
   if(sample_name== "Gluino_13TeV16_M1600_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.13721e+00; Waa=1.12260e+00; Wan=1.05957e+00;}
   if(sample_name== "Gluino_13TeV16G_M1600_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.15968e+00; Waa=1.11867e+00; Wan=1.05943e+00;}
   if(sample_name== "Gluino_13TeV16_M1600N_f50"){Rhadron=1; Wa=1.40489e+00; Wad=1.05869e+00; Waa=1.12438e+00; Wan=1.05205e+00;}
   if(sample_name== "Gluino_13TeV16G_M1600N_f50"){Rhadron=1; Wa=1.36716e+00; Wad=1.03935e+00; Waa=1.12433e+00; Wan=1.04815e+00;}
   if(sample_name== "Gluino_13TeV16_M1800_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.04399e+00; Waa=1.13701e+00; Wan=1.05993e+00;}
   if(sample_name== "Gluino_13TeV16G_M1800_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.02927e+00; Waa=1.15483e+00; Wan=1.06635e+00;}
   if(sample_name== "Gluino_13TeV16_M1800N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.03688e+00; Waa=1.11329e+00; Wan=1.05435e+00;}
   if(sample_name== "Gluino_13TeV16G_M1800N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.02108e+00; Waa=1.10542e+00; Wan=1.05980e+00;}
   if(sample_name== "Gluino_13TeV16_M1800_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.04399e+00; Waa=1.13700e+00; Wan=1.05994e+00;}
   if(sample_name== "Gluino_13TeV16G_M1800_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.02927e+00; Waa=1.15483e+00; Wan=1.06634e+00;}
   if(sample_name== "Gluino_13TeV16_M1800N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.03688e+00; Waa=1.11329e+00; Wan=1.05435e+00;}
   if(sample_name== "Gluino_13TeV16G_M1800N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.02108e+00; Waa=1.10541e+00; Wan=1.05979e+00;}
   if(sample_name== "Gluino_13TeV16_M2000_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.05497e+00; Waa=1.11383e+00; Wan=1.05472e+00;}
   if(sample_name== "Gluino_13TeV16G_M2000_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.05321e+00; Waa=1.10173e+00; Wan=1.05007e+00;}
   if(sample_name== "Gluino_13TeV16_M2000N_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.05937e+00; Waa=1.11855e+00; Wan=1.04864e+00;}
   if(sample_name== "Gluino_13TeV16G_M2000N_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.05173e+00; Waa=1.11570e+00; Wan=1.04457e+00;}
   if(sample_name== "Gluino_13TeV16_M2000_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.05497e+00; Waa=1.11383e+00; Wan=1.05471e+00;}
   if(sample_name== "Gluino_13TeV16G_M2000_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.05321e+00; Waa=1.10173e+00; Wan=1.05007e+00;}
   if(sample_name== "Gluino_13TeV16_M2000N_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.05937e+00; Waa=1.11856e+00; Wan=1.04862e+00;}
   if(sample_name== "Gluino_13TeV16G_M2000N_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.05173e+00; Waa=1.11569e+00; Wan=1.04457e+00;}
   if(sample_name== "Gluino_13TeV16_M2200_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.06535e+00; Waa=1.13193e+00; Wan=1.05624e+00;}
   if(sample_name== "Gluino_13TeV16G_M2200_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.06233e+00; Waa=1.15471e+00; Wan=1.05806e+00;}
   if(sample_name== "Gluino_13TeV16_M2200N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.05069e+00; Waa=1.11169e+00; Wan=1.04728e+00;}
   if(sample_name== "Gluino_13TeV16G_M2200N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.03741e+00; Waa=1.10548e+00; Wan=1.04138e+00;}
   if(sample_name== "Gluino_13TeV16_M2200_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.06535e+00; Waa=1.13192e+00; Wan=1.05624e+00;}
   if(sample_name== "Gluino_13TeV16G_M2200_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.06233e+00; Waa=1.15470e+00; Wan=1.05806e+00;}
   if(sample_name== "Gluino_13TeV16_M2200N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.05069e+00; Waa=1.11168e+00; Wan=1.04728e+00;}
   if(sample_name== "Gluino_13TeV16G_M2200N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.03741e+00; Waa=1.10548e+00; Wan=1.04138e+00;}
   if(sample_name== "Gluino_13TeV16_M2400_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.09729e+00; Waa=1.11488e+00; Wan=1.05083e+00;}
   if(sample_name== "Gluino_13TeV16G_M2400_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.11968e+00; Waa=1.09763e+00; Wan=1.04302e+00;}
   if(sample_name== "Gluino_13TeV16_M2400N_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.02856e+00; Waa=1.11822e+00; Wan=1.05504e+00;}
   if(sample_name== "Gluino_13TeV16G_M2400N_f10"){Rhadron=1; Wa=      1.0 ; Wad=1.01861e+00; Waa=1.10373e+00; Wan=1.06239e+00;}
   if(sample_name== "Gluino_13TeV16_M2400_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.09729e+00; Waa=1.11488e+00; Wan=1.05083e+00;}
   if(sample_name== "Gluino_13TeV16G_M2400_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.11968e+00; Waa=1.09763e+00; Wan=1.04302e+00;}
   if(sample_name== "Gluino_13TeV16_M2400N_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.02856e+00; Waa=1.11822e+00; Wan=1.05505e+00;}
   if(sample_name== "Gluino_13TeV16G_M2400N_f50"){Rhadron=1; Wa=      1.0 ; Wad=1.01861e+00; Waa=1.10372e+00; Wan=1.06238e+00;}
   if(sample_name== "Gluino_13TeV16_M2600_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.03106e+00; Waa=1.11380e+00; Wan=1.04804e+00;}
   if(sample_name== "Gluino_13TeV16G_M2600_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.01947e+00; Waa=1.10963e+00; Wan=1.04575e+00;}
   if(sample_name== "Gluino_13TeV16_M2600N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.09505e+00; Waa=1.11331e+00; Wan=1.05058e+00;}
   if(sample_name== "Gluino_13TeV16G_M2600N_f10"){Rhadron=1; Wa=1.00000e+00; Wad=1.10710e+00; Waa=1.10905e+00; Wan=1.05031e+00;}
   if(sample_name== "Gluino_13TeV16_M2600_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.03106e+00; Waa=1.11379e+00; Wan=1.04805e+00;}
   if(sample_name== "Gluino_13TeV16G_M2600_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.01947e+00; Waa=1.10963e+00; Wan=1.04575e+00;}
   if(sample_name== "Gluino_13TeV16_M2600N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.09505e+00; Waa=1.11331e+00; Wan=1.05059e+00;}
   if(sample_name== "Gluino_13TeV16G_M2600N_f50"){Rhadron=1; Wa=1.00000e+00; Wad=1.10710e+00; Waa=1.10905e+00; Wan=1.05031e+00;}
   if(sample_name== "Stop_13TeV16_M100"){Rhadron=1; Wa=1.48895e+00; Wad=1.63773e+00; Waa=3.72100e+00; Wan=1.58986e+00;}
   if(sample_name== "Stop_13TeV16G_M100"){Rhadron=1; Wa=1.85460e+00; Wad=1.79354e+00; Waa=3.57042e+00; Wan=1.59750e+00;}
   if(sample_name== "Stop_13TeV16_M100N"){Rhadron=1; Wa=1.52743e+00; Wad=1.39185e+00; Waa=3.89317e+00; Wan=1.58947e+00;}
   if(sample_name== "Stop_13TeV16G_M100N"){Rhadron=1; Wa=1.47291e+00; Wad=1.27392e+00; Waa=3.88127e+00; Wan=1.56703e+00;}
   if(sample_name== "Stop_13TeV16_M200"){Rhadron=1; Wa=1.41378e+00; Wad=1.53085e+00; Waa=3.69431e+00; Wan=1.58899e+00;}
   if(sample_name== "Stop_13TeV16G_M200"){Rhadron=1; Wa=1.31263e+00; Wad=1.38737e+00; Waa=3.69745e+00; Wan=1.57685e+00;}
   if(sample_name== "Stop_13TeV16_M200N"){Rhadron=1; Wa=1.94239e+00; Wad=1.65615e+00; Waa=3.78488e+00; Wan=1.57685e+00;}
   if(sample_name== "Stop_13TeV16G_M200N"){Rhadron=1; Wa=2.06572e+00; Wad=1.84068e+00; Waa=3.59161e+00; Wan=1.56411e+00;}
   if(sample_name== "Stop_13TeV16_M400"){Rhadron=1; Wa=1.48638e+00; Wad=1.77125e+00; Waa=3.91371e+00; Wan=1.57659e+00;}
   if(sample_name== "Stop_13TeV16G_M400"){Rhadron=1; Wa=1.34605e+00; Wad=1.70643e+00; Waa=4.05354e+00; Wan=1.59001e+00;}
   if(sample_name== "Stop_13TeV16_M400N"){Rhadron=1; Wa=3.33978e+00; Wad=1.60340e+00; Waa=3.65370e+00; Wan=1.57817e+00;}
   if(sample_name== "Stop_13TeV16G_M400N"){Rhadron=1; Wa=5.29832e+00; Wad=1.54759e+00; Waa=3.50364e+00; Wan=1.57003e+00;}
   if(sample_name== "Stop_13TeV16_M600"){Rhadron=1; Wa=1.64779e+00; Wad=1.73679e+00; Waa=3.71421e+00; Wan=1.57603e+00;}
   if(sample_name== "Stop_13TeV16G_M600"){Rhadron=1; Wa=1.62139e+00; Wad=2.00693e+00; Waa=3.59236e+00; Wan=1.56359e+00;}
   if(sample_name== "Stop_13TeV16_M600N"){Rhadron=1; Wa=1.29048e+00; Wad=1.64931e+00; Waa=4.03799e+00; Wan=1.60249e+00;}
   if(sample_name== "Stop_13TeV16G_M600N"){Rhadron=1; Wa=1.25517e+00; Wad=1.67057e+00; Waa=4.32971e+00; Wan=1.62918e+00;}
   if(sample_name== "Stop_13TeV16_M800"){Rhadron=1; Wa=1.33778e+00; Wad=1.60026e+00; Waa=3.83396e+00; Wan=1.58746e+00;}
   if(sample_name== "Stop_13TeV16G_M800"){Rhadron=1; Wa=1.29931e+00; Wad=1.59535e+00; Waa=3.83467e+00; Wan=1.58766e+00;}
   if(sample_name== "Stop_13TeV16_M800N"){Rhadron=1; Wa=2.02103e+00; Wad=1.59819e+00; Waa=3.86165e+00; Wan=1.59136e+00;}
   if(sample_name== "Stop_13TeV16G_M800N"){Rhadron=1; Wa=2.13045e+00; Wad=1.82254e+00; Waa=3.82256e+00; Wan=1.58237e+00;}
   if(sample_name== "Stop_13TeV16_M1000"){Rhadron=1; Wa=2.43539e+00; Wad=1.61244e+00; Waa=3.79207e+00; Wan=1.57390e+00;}
   if(sample_name== "Stop_13TeV16G_M1000"){Rhadron=1; Wa=2.29364e+00; Wad=1.68773e+00; Waa=3.62179e+00; Wan=1.57279e+00;}
   if(sample_name== "Stop_13TeV16_M1000N"){Rhadron=1; Wa=1.14838e+01; Wad=1.48860e+00; Waa=3.88234e+00; Wan=1.55547e+00;}
   if(sample_name== "Stop_13TeV16G_M1000N"){Rhadron=1; Wa=1.91021e+01; Wad=1.41199e+00; Waa=4.06618e+00; Wan=1.54431e+00;}
   if(sample_name== "Stop_13TeV16_M1200"){Rhadron=1; Wa=1.76112e+00; Wad=1.74799e+00; Waa=3.79217e+00; Wan=1.58862e+00;}
   if(sample_name== "Stop_13TeV16G_M1200"){Rhadron=1; Wa=1.59545e+00; Wad=1.84435e+00; Waa=3.71432e+00; Wan=1.59049e+00;}
   if(sample_name== "Stop_13TeV16_M1200N"){Rhadron=1; Wa=3.58410e+00; Wad=1.47529e+00; Waa=3.80086e+00; Wan=1.59485e+00;}
   if(sample_name== "Stop_13TeV16G_M1200N"){Rhadron=1; Wa=3.17943e+00; Wad=1.47822e+00; Waa=3.70784e+00; Wan=1.59335e+00;}
   if(sample_name== "Stop_13TeV16_M1400"){Rhadron=1; Wa=1.57971e+00; Wad=1.50758e+00; Waa=3.83374e+00; Wan=1.56841e+00;}
   if(sample_name== "Stop_13TeV16G_M1400"){Rhadron=1; Wa=1.66959e+00; Wad=1.36557e+00; Waa=3.99336e+00; Wan=1.52871e+00;}
   if(sample_name== "Stop_13TeV16_M1600"){Rhadron=1; Wa=1.00000e+00; Wad=1.63048e+00; Waa=3.90841e+00; Wan=1.57668e+00;}
   if(sample_name== "Stop_13TeV16G_M1600"){Rhadron=1; Wa=1.00000e+00; Wad=1.55628e+00; Waa=3.97900e+00; Wan=1.56856e+00;}
   if(sample_name== "Stop_13TeV16_M1600N"){Rhadron=1; Wa=1.00000e+00; Wad=1.74601e+00; Waa=3.88697e+00; Wan=1.61307e+00;}
   if(sample_name== "Stop_13TeV16G_M1600N"){Rhadron=1; Wa=1.00000e+00; Wad=1.67093e+00; Waa=4.04171e+00; Wan=1.64330e+00;}
   if(sample_name== "Stop_13TeV16_M1800"){Rhadron=1; Wa=      1.0 ; Wad=1.71460e+00; Waa=3.75776e+00; Wan=1.58930e+00;}
   if(sample_name== "Stop_13TeV16G_M1800"){Rhadron=1; Wa=      1.0 ; Wad=1.76000e+00; Waa=3.76333e+00; Wan=1.62603e+00;}
   if(sample_name== "Stop_13TeV16_M1800N"){Rhadron=1; Wa=1.00000e+00; Wad=1.67225e+00; Waa=3.68948e+00; Wan=1.58159e+00;}
   if(sample_name== "Stop_13TeV16G_M1800N"){Rhadron=1; Wa=1.00000e+00; Wad=1.67238e+00; Waa=3.84747e+00; Wan=1.59714e+00;}
   if(sample_name== "Stop_13TeV16_M2000"){Rhadron=1; Wa=1.00000e+00; Wad=1.55812e+00; Waa=3.94397e+00; Wan=1.57747e+00;}
   if(sample_name== "Stop_13TeV16G_M2000"){Rhadron=1; Wa=1.00000e+00; Wad=1.48680e+00; Waa=3.99095e+00; Wan=1.57254e+00;}
   if(sample_name== "Stop_13TeV16_M2000N"){Rhadron=1; Wa=1.00000e+00; Wad=1.50401e+00; Waa=3.83539e+00; Wan=1.56599e+00;}
   if(sample_name== "Stop_13TeV16G_M2000N"){Rhadron=1; Wa=1.00000e+00; Wad=1.45216e+00; Waa=3.79238e+00; Wan=1.53687e+00;}
   if(sample_name== "Stop_13TeV16_M2200"){Rhadron=1; Wa=      1.0 ; Wad=1.59142e+00; Waa=3.86218e+00; Wan=1.54876e+00;}
   if(sample_name== "Stop_13TeV16G_M2200"){Rhadron=1; Wa=      1.0 ; Wad=1.57193e+00; Waa=3.93874e+00; Wan=1.52037e+00;}
   if(sample_name== "Stop_13TeV16_M2200N"){Rhadron=1; Wa=1.00000e+00; Wad=1.47373e+00; Waa=3.69419e+00; Wan=1.58923e+00;}
   if(sample_name== "Stop_13TeV16G_M2200N"){Rhadron=1; Wa=1.00000e+00; Wad=1.48432e+00; Waa=3.66894e+00; Wan=1.61450e+00;}
   if(sample_name== "Stop_13TeV16_M2400"){Rhadron=1; Wa=      1.0 ; Wad=1.58353e+00; Waa=4.21290e+00; Wan=1.56575e+00;}
   if(sample_name== "Stop_13TeV16G_M2400"){Rhadron=1; Wa=      1.0 ; Wad=1.53402e+00; Waa=4.37498e+00; Wan=1.51360e+00;}
   if(sample_name== "Stop_13TeV16_M2400N"){Rhadron=1; Wa=1.52334e+01; Wad=1.64971e+00; Waa=3.80043e+00; Wan=1.58422e+00;}
   if(sample_name== "Stop_13TeV16G_M2400N"){Rhadron=1; Wa=1.33230e+02; Wad=1.74531e+00; Waa=3.85981e+00; Wan=1.57615e+00;}
   if(sample_name== "Stop_13TeV16_M2600"){Rhadron=1; Wa=1.40679e+00; Wad=1.62079e+00; Waa=3.80874e+00; Wan=1.56877e+00;}
   if(sample_name== "Stop_13TeV16G_M2600"){Rhadron=1; Wa=1.36173e+00; Wad=1.72353e+00; Waa=3.91675e+00; Wan=1.54487e+00;}
   if(sample_name== "Stop_13TeV16_M2600N"){Rhadron=1; Wa=3.45642e+00; Wad=1.59418e+00; Waa=3.88138e+00; Wan=1.60157e+00;}
   if(sample_name== "Stop_13TeV16G_M2600N"){Rhadron=1; Wa=4.42824e+00; Wad=1.58612e+00; Waa=3.94492e+00; Wan=1.61766e+00;}


   //------------ weights for R-hadron samples ----------Stop

   unsigned int nw=0, na=0, nd=0, nn=0; //initialize counters: nw - wrong, na - other, nd - double charged, nn - neutral

// R-hadron wights needed due to wrong GenId----------------------------------END


   // int nrha=0; 
     for(unsigned int g=0;g<genColl.size();g++) {
       if(genColl[g].pt()<5)continue;
       if(genColl[g].status()!=1)continue;
       int AbsPdg=abs(genColl[g].pdgId());
       if(AbsPdg<1000000 && AbsPdg!=17)continue;

// categorise event with R-hadrons for additional weighting-----------------------BEGIN

     int GenId=genColl[g].pdgId();
     
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
     
// categorise event with R-hadrons for additional weighting-----------------------END

     SamplePlots->genlevelpT->Fill(genColl[g].pt(), Event_Weight);
     SamplePlots->genleveleta->Fill(genColl[g].eta(), Event_Weight);
     SamplePlots->genlevelbeta->Fill(genColl[g].p()/genColl[g].energy(), Event_Weight);

     /*
     nrha++;
//mk rhadron ntuple
    if(isSignal)
       stPlots_GenFillTree(SamplePlots, ev.eventAuxiliary().run(),ev.eventAuxiliary().event(),ev.eventAuxiliary().luminosityBlock(), nrha, Event_Weight,genColl[g].pdgId(),genColl[g].charge(),genColl[g].mass(),genColl[g].pt(),genColl[g].eta(),genColl[g].phi(),-1);

     } // remove } below if uncommenting
     nrha=0;
     */

     }

// additional weighting for R-hadrons---------------------------------------------BEGIN

/*
     if(nw>0) continue; // to terminate the event
     if(nd==1&&na==1) Event_Weight*=Wad;       // double_charged + other R-hadron
     if(nn==1&&na==1) Event_Weight*=Wan;       // neutral        + other R-hadron
     if(nd==0&&nn==0&&na==1) Event_Weight*=Wa ;       // single other R-hadron (the second lost in action) //nn==0 was missing !!! 
     if(       na==2) Event_Weight*=Waa;       // 2 others in the event
// additional weighting for R-hadrons---------------------------------------------BEGIN

*///mkrh test of not removing events
            }

	    // new genHSCP ntuple after correcting weights

	    int nrha=0; 
	    for(unsigned int g=0;g<genColl.size();g++) {
	      if(genColl[g].pt()<5)continue;
	      if(genColl[g].status()!=1)continue;
	      int AbsPdg=abs(genColl[g].pdgId());
	      if(AbsPdg<1000000 && AbsPdg!=17)continue;

	      nrha++;
	      //mk rhadron ntuple
	      if(isSignal)
		stPlots_GenFillTree(SamplePlots, ev.eventAuxiliary().run(),ev.eventAuxiliary().event(),ev.eventAuxiliary().luminosityBlock(), nrha, Event_Weight,genColl[g].pdgId(),genColl[g].charge(),genColl[g].mass(),genColl[g].pt(),genColl[g].eta(),genColl[g].phi(),-1);

	    }
	    nrha=0;



            //check if the event is passing trigger
            SamplePlots      ->TotalE  ->Fill(0.0,Event_Weight);  
            if(isMC)MCTrPlots->TotalE  ->Fill(0.0,Event_Weight);
            SamplePlots      ->TotalEPU->Fill(0.0,Event_Weight*PUSystFactor);
            if(isMC)MCTrPlots->TotalEPU->Fill(0.0,Event_Weight*PUSystFactor);
	    //See if event passed signal triggers
            if(!PassTrigger(ev, isData, false, (is2016&&!is2016G)?&L1Emul:NULL) ) {
	      //For TOF only analysis if the event doesn't pass the signal triggers check if it was triggered by the no BPTX cosmic trigger
	      //If not TOF only then move to next event
	      if(TypeMode!=3) continue;
	      if(!PassTrigger(ev, isData, true, (is2016&&!is2016G)?&L1Emul:NULL)) continue;

	      //If is cosmic event then switch plots to use to the ones for cosmics
	      SamplePlots=&plotsMap[CosmicName];
	    }else if(TypeMode==3) {
	      SamplePlots = &plotsMap[samples[s].Name];
	    }

            SamplePlots       ->TotalTE->Fill(0.0,Event_Weight);
            if(isMC)MCTrPlots ->TotalTE->Fill(0.0,Event_Weight);

            //keep beta distribution for signal
            if(isSignal){if(HSCPGenBeta1>=0)SamplePlots->Beta_Triggered->Fill(HSCPGenBeta1, Event_Weight); if(HSCPGenBeta2>=0)SamplePlots->Beta_Triggered->Fill(HSCPGenBeta2, Event_Weight);}

            //load all event collection that will be used later on (HSCP COll, dEdx and TOF)
            fwlite::Handle<susybsm::HSCParticleCollection> hscpCollH;
            hscpCollH.getByLabel(ev,"HSCParticleProducer");
            if(!hscpCollH.isValid()){printf("HSCParticle Collection NotFound\n");continue;}
            const susybsm::HSCParticleCollection& hscpColl = *hscpCollH;

            fwlite::Handle<DeDxHitInfoAss> dedxCollH;
            dedxCollH.getByLabel(ev, "dedxHitInfo");
            if(!dedxCollH.isValid()){printf("Invalid dedxCollH\n");continue;}

            fwlite::Handle<MuonTimeExtraMap> TOFCollH;
            TOFCollH.getByLabel(ev, "muons",TOF_Label.c_str());
            if(!TOFCollH.isValid()){printf("Invalid TOF collection\n");return;}

            fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
            TOFDTCollH.getByLabel(ev, "muons",TOFdt_Label.c_str());
            if(!TOFDTCollH.isValid()){printf("Invalid DT TOF collection\n");return;}

            fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
            TOFCSCCollH.getByLabel(ev, "muons",TOFcsc_Label.c_str());
            if(!TOFCSCCollH.isValid()){printf("Invalid CSC TOF collection\n");return;}

            fwlite::Handle<CSCSegmentCollection> CSCSegmentCollHandle;
            fwlite::Handle<DTRecSegment4DCollection> DTSegmentCollHandle;            
            if(!isMC){ //do not reocmpute TOF on MC background
               CSCSegmentCollHandle.getByLabel(ev, "cscSegments");
               if(!CSCSegmentCollHandle.isValid()){printf("CSC Segment Collection not found!\n"); continue;}

               DTSegmentCollHandle.getByLabel(ev, "dt4DSegments");
               if(!DTSegmentCollHandle.isValid()){printf("DT Segment Collection not found!\n"); continue;}
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

            HIPemulator.setEventRate(); //take it from a pdf
            HIPemulatorUp.setEventRate(HIPemulator.getEventRatePixel()*1.25, HIPemulator.getEventRateStrip()*1.80);  // deltaPixel = 3.653981e+02, basePixel = 1.332625e+03; deltaStrip = 4.662832e+02, baseStrip = 5.958308e+02, from Run257805
            HIPemulatorDown.setEventRate(HIPemulator.getEventRatePixel()*0.75, HIPemulator.getEventRateStrip()*0.20); 

	    HIPTrackLossEmul.SetHIPTrackLossRate(ev);
//           if (HIPemulator.getEventRatePixel()>0 && HIPemulator.getEventRateStrip()>0)
//              fprintf(stderr, "HIPs: %lf\t%lf\t%lf\t%lf\t%lf\t%lf\n", HIPemulator.getEventRatePixel(), HIPemulatorUp.getEventRatePixel(), HIPemulatorDown.getEventRatePixel(), HIPemulator.getEventRateStrip(), HIPemulatorUp.getEventRateStrip(), HIPemulatorDown.getEventRateStrip());

            //loop on HSCP candidates
            for(unsigned int c=0;c<hscpColl.size();c++){
               //define alias for important variable
               susybsm::HSCParticle hscp  = hscpColl[c];
               reco::MuonRef  muon  = hscp.muonRef();

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
               if(isSignal && DistToHSCP(hscp, genColl, ClosestGen)>0.03)continue;

	       // we are losing some tracks due to HIP
	       if(!isData && is2016 && !is2016G && !HIPTrackLossEmul.TrackSurvivesHIPInefficiency()) continue;

               //load quantity associated to this track (TOF and dEdx)
               const DeDxHitInfo* dedxHits = NULL;
	       if(TypeMode!=3 && !track.isNull()) {
                  DeDxHitInfoRef dedxHitsRef = dedxCollH->get(track.key());		 
                  if(!dedxHitsRef.isNull())dedxHits = &(*dedxHitsRef);
	       }
               const reco::MuonTimeExtra* tof = NULL;
               const reco::MuonTimeExtra* dttof = NULL;
               const reco::MuonTimeExtra* csctof = NULL;
               if(TypeMode>1 && TypeMode!=5 && !hscp.muonRef().isNull()){
                  if(isMC){
                     tof  = &TOFCollH->get(hscp.muonRef().key()); dttof = &TOFDTCollH->get(hscp.muonRef().key());  csctof = &TOFCSCCollH->get(hscp.muonRef().key());
                  }else{
                     const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
                     const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;
//std::cout<<"TESTA\n";
                     tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, isData?1:0 ); //apply T0 correction on data but not on signal MC
//std::cout<<"TESTB\n";
                     tof  = &tofCalculator.combinedTOF; dttof = &tofCalculator.dtTOF;  csctof = &tofCalculator.cscTOF;
//std::cout<<"TESTC\n";
                  }
               }
if(!dedxHits) continue; // skip tracks without hits otherwise there will be a crash
            HitDeDxCollection hitDeDx = getHitDeDx(dedxHits, dEdxSF, trackerCorrector.TrackerGains, false, 1);

          
unsigned int codep = 0;

if(isSignal){
 codep = genColl[ClosestGen].pdgId();
//mk  std::cout<<"GenId  " <<codep<<std::endl;
 }

               
	       double dEdxErr = 0;
               DeDxData dedxSObjTmp  = computedEdx(dedxHits, dEdxSF, dEdxTemplates, true, useClusterCleaning, TypeMode==5, false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.00, NULL,0,codep);
               DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulator:NULL, &dEdxErr,codep);
               DeDxData dedxMUpObjTmp = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulatorUp:NULL,0,codep);
               DeDxData dedxMDownObjTmp = computedEdx(dedxHits, dEdxSF, NULL,          true, useClusterCleaning, false      , false, trackerCorrector.TrackerGains, true, true, 99, false, 1, 0.15, (!isData && !is2016G)?&HIPemulatorDown:NULL,0,codep);
               DeDxData* dedxSObj  = dedxSObjTmp.numberOfMeasurements()>0?&dedxSObjTmp:NULL;
               DeDxData* dedxMObj  = dedxMObjTmp.numberOfMeasurements()>0?&dedxMObjTmp:NULL;
               DeDxData* dedxMUpObj = dedxMUpObjTmp.numberOfMeasurements()>0?&dedxMUpObjTmp:NULL;
               DeDxData* dedxMDownObj = dedxMDownObjTmp.numberOfMeasurements()>0?&dedxMDownObjTmp:NULL;
               if(TypeMode==5)OpenAngle = deltaROpositeTrack(hscpColl, hscp); //OpenAngle is a global variable... that's uggly C++, but that's the best I found so far

               //compute systematic uncertainties on signal
               if(isSignal){
                  //FIXME to be measured on 2015 data, currently assume 2012
                  bool   PRescale = true;
                  double IRescale =-0.05; // added to the Ias value
                  double MRescale = 0.95;
		  double TRescale =-0.015; //-0.005 (used in 2012); // added to the 1/beta value
		  
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
                  if (genpT>0) {  SamplePlots->genrecopT->Fill(genpT, track->pt()); }
		  

                  //DEBUG TOF
/*                  if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   PRescale, 0, 0)){
                     const reco::MuonTimeExtra* dttofaod = &TOFDTCollH->get(hscp.muonRef().key());
                     if((dttof->inverseBeta()>0 && dttof->inverseBeta()<0.8 && fabs( dttof->inverseBeta()-dttofaod->inverseBeta())>0.2) ){
                        printf("event = %i\n", (int) ientry);
                        printf("OTF %f +- %f with %i NDOF\n", dttof->inverseBeta(), dttof->inverseBetaErr(), dttof->nDof());
                        printf("AOD %f +- %f with %i NDOF\n", dttofaod->inverseBeta(), dttofaod->inverseBetaErr(), dttofaod->nDof());



                        const CSCSegmentCollection& CSCSegmentColl = *CSCSegmentCollHandle;
                        const DTRecSegment4DCollection& DTSegmentColl = *DTSegmentCollHandle;
                        tofCalculator.computeTOF(muon, CSCSegmentColl, DTSegmentColl, isData?1:0, true ); //apply T0 correction on data but not on signal MC

                        exit(0);
                     }
                  }
*/


                  // compute systematic due to momentum scale
                  if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   PRescale, 0, 0)){
                     double RescalingFactor = RescaledPt(track->pt(),track->eta(),track->phi(),track->charge())/track->pt();

                     if(TypeMode==5 && isSemiCosmicSB)continue;
 		     double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p()*RescalingFactor,dedxMObj->dEdx(),!isData);
		     double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p()*RescalingFactor,tof->inverseBeta());
		     double MassComb = -1;
		     if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p()*RescalingFactor, (GetIBeta(dedxMObj->dEdx(),!isData) + (1/tof->inverseBeta()))*0.5);
		     else if(dedxMObj) MassComb = Mass;
		     if(tof) MassComb=MassTOF;

                     for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
		       if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, false, -1,   PRescale, 0, 0)){
                           HSCPTk_SystP[CutIndex] = true;
                           if(Mass>MaxMass_SystP[CutIndex]) MaxMass_SystP[CutIndex]=Mass;
                           SamplePlots->Mass_SystP->Fill(CutIndex, Mass,Event_Weight);
                           if(tof){
                              SamplePlots->MassTOF_SystP ->Fill(CutIndex, MassTOF , Event_Weight);
                           }
                           SamplePlots->MassComb_SystP->Fill(CutIndex, MassComb, Event_Weight);
                        }
                     }
                  }

                  // compute systematic due to dEdx (both Ias and Ih)
                  if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, IRescale, 0)){
                     if(TypeMode==5 && isSemiCosmicSB)continue;
		     double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx()*MRescale,!isData);
		     double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
		     double MassComb = -1;
		     if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),!isData) + (1/tof->inverseBeta()))*0.5);
		     else if(dedxMObj) MassComb = Mass;
		     if(tof) MassComb=MassTOF;

                     for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
		       if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, false, -1,   0, IRescale, 0)){
                           HSCPTk_SystI[CutIndex] = true;
                           if(Mass>MaxMass_SystI[CutIndex]) MaxMass_SystI[CutIndex]=Mass;
                           SamplePlots->Mass_SystI->Fill(CutIndex, Mass,Event_Weight);
                           if(tof){
                              SamplePlots->MassTOF_SystI ->Fill(CutIndex, MassTOF , Event_Weight);
                           }
                           SamplePlots->MassComb_SystI->Fill(CutIndex, MassComb, Event_Weight);
                        }
                     }
                  }




                  // compute systematic due to Mass shift
                  if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, 0, 0)){
                     if(TypeMode==5 && isSemiCosmicSB)continue;
		     double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx()*MRescale,!isData);
		     double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
		     double MassComb = -1;
		     if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx()*MRescale,!isData) + (1/tof->inverseBeta()))*0.5);
		     else if(dedxMObj) MassComb = Mass;
		     if(tof) MassComb=MassTOF;

                     for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
		       if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, false, -1,   0, 0, 0)){
                           HSCPTk_SystM[CutIndex] = true;
                           if(Mass>MaxMass_SystM[CutIndex]) MaxMass_SystM[CutIndex]=Mass;
                           SamplePlots->Mass_SystM->Fill(CutIndex, Mass,Event_Weight);
                           if(tof){
                              SamplePlots->MassTOF_SystM ->Fill(CutIndex, MassTOF , Event_Weight);
                           }
                           SamplePlots->MassComb_SystM->Fill(CutIndex, MassComb, Event_Weight);
                        }
                     }
                  }

                  // compute systematic due to TOF
                  if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, 0, TRescale)){
                     if(TypeMode==5 && isSemiCosmicSB)continue;
 		     double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx(),!isData);
		     double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),(tof->inverseBeta()+TRescale));
		     double MassComb = -1;
		     if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),!isData) + (1/(tof->inverseBeta()+TRescale)))*0.5);
		     else if(dedxMObj) MassComb = Mass;
		     if(tof) MassComb=MassTOF;

                     for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
		       if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, false, -1,   0, 0, TRescale)){
                           HSCPTk_SystT[CutIndex] = true;
                           if(Mass>MaxMass_SystT[CutIndex]) MaxMass_SystT[CutIndex]=Mass;
                           SamplePlots->Mass_SystT->Fill(CutIndex, Mass,Event_Weight);
                           if(tof){
                              SamplePlots->MassTOF_SystT ->Fill(CutIndex, MassTOF , Event_Weight);
                           }
                           SamplePlots->MassComb_SystT->Fill(CutIndex, MassComb, Event_Weight);
                        }
                     }
                  }

                  // compute systematics due to PU
                  if(PassPreselection( hscp,  dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev,  NULL, -1,   0, 0, 0)){
                     if(TypeMode==5 && isSemiCosmicSB)continue;
		     double Mass     = -1; if(dedxMObj) Mass=GetMass(track->p(),dedxMObj->dEdx(),!isData);
		     double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
		     double MassComb = -1;
		     if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),!isData) + (1/tof->inverseBeta()))*0.5);
		     else if(dedxMObj) MassComb = Mass;
		     if(tof) MassComb=MassTOF;

                     for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
		       if(PassSelection(hscp,  dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, false, -1,   0, 0, 0)){
                           HSCPTk_SystPU[CutIndex] = true;
                           if(Mass>MaxMass_SystPU[CutIndex]) MaxMass_SystPU[CutIndex]=Mass;
                           SamplePlots->Mass_SystPU->Fill(CutIndex, Mass,Event_Weight*PUSystFactor);
                           if(tof){
                              SamplePlots->MassTOF_SystPU ->Fill(CutIndex, MassTOF , Event_Weight*PUSystFactor);
                           }
                           SamplePlots->MassComb_SystPU->Fill(CutIndex, MassComb, Event_Weight*PUSystFactor);
                        }
                     }
                  }
               }//End of systematic computation for signal

               //check if the canddiate pass the preselection cuts
               if(isMC)PassPreselection( hscp, dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev, MCTrPlots  , -1, false, 0, 0, GetMassErr (track->p(), track->ptError(), dedxMObj?dedxMObj->dEdx():-1, dEdxErr, GetMass(track->p(), dedxMObj?dedxMObj->dEdx():-1, !isData)));
               if(    !PassPreselection( hscp, dedxHits, dedxSObj, dedxMObj, tof, dttof, csctof, ev, SamplePlots, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1, false, 0, 0, GetMassErr (track->p(), track->ptError(), dedxMObj?dedxMObj->dEdx():-1, dEdxErr, GetMass(track->p(), dedxMObj?dedxMObj->dEdx():-1, !isData)))) continue;
               if(TypeMode==5 && isSemiCosmicSB)continue;

               //fill the ABCD histograms and a few other control plots
               if(isData)Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, SamplePlots);
	       else if(isMC) Analysis_FillControlAndPredictionHist(hscp, dedxSObj, dedxMObj, tof, MCTrPlots);
               if(TypeMode==5 && isCosmicSB)continue; 

	       //Find the number of tracks passing selection for TOF<1 that will be used to check the background prediction
	       if(isMC || isData) {
               //compute the mass of the candidate, for TOF mass flip the TOF over 1 to get the mass, so 0.8->1.2
		 double Mass     = -1; if(dedxMObj) Mass = GetMass(track->p(),dedxMObj->dEdx(),!isData);
		 double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),(2-tof->inverseBeta()));
		 double MassComb = -1;
		 if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),!isData) + (1/(2-tof->inverseBeta())))*0.5 ) ;
		 if(dedxMObj) MassComb = Mass;
		 if(tof)MassComb=GetMassFromBeta(track->p(),(1/(2-tof->inverseBeta())));

		 for(unsigned int CutIndex=0;CutIndex<CutPt_Flip.size();CutIndex++){
		 //Background check looking at region with TOF<1
		   if(!PassSelection   (hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, true)) continue;
                  //Fill Mass Histograms

                  if(isMC)MCTrPlots->Mass_Flip->Fill(CutIndex, Mass,Event_Weight);
                  SamplePlots      ->Mass_Flip->Fill(CutIndex, Mass,Event_Weight);
                  if(tof){
                  if(isMC)MCTrPlots->MassTOF_Flip->Fill(CutIndex, MassTOF, Event_Weight);
                     SamplePlots   ->MassTOF_Flip->Fill(CutIndex, MassTOF, Event_Weight);
                  }
                  if(isMC)MCTrPlots->MassComb_Flip->Fill(CutIndex, MassComb, Event_Weight);
                  SamplePlots      ->MassComb_Flip->Fill(CutIndex, MassComb, Event_Weight);
		 }
	       }

               //compute the mass of the candidate
	       double Mass     = -1; if(dedxMObj) Mass = GetMass(track->p(),dedxMObj->dEdx(),!isData);
	       double MassTOF  = -1; if(tof)MassTOF = GetTOFMass(track->p(),tof->inverseBeta());
	       double MassComb = -1;
	       if(tof && dedxMObj)MassComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMObj->dEdx(),!isData) + (1/tof->inverseBeta()))*0.5 ) ;
	       if(dedxMObj) MassComb = Mass;
	       if(tof)MassComb=GetMassFromBeta(track->p(),(1/tof->inverseBeta()));


               double MassUp    = -1; if(dedxMUpObj) MassUp=GetMass(track->p(),dedxMUpObj->dEdx(),!isData);
               double MassUpComb = -1;
               if(tof && dedxMUpObj)MassUpComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMUpObj->dEdx(),!isData) + (1/tof->inverseBeta()))*0.5 ) ;
               if(dedxMUpObj) MassUpComb = MassUp;
               if(tof)MassUpComb=GetMassFromBeta(track->p(),(1/tof->inverseBeta()));

               double MassDown    = -1; if(dedxMDownObj) MassDown=GetMass(track->p(),dedxMDownObj->dEdx(),!isData);
               double MassDownComb = -1;
               if(tof && dedxMDownObj)MassDownComb=GetMassFromBeta(track->p(), (GetIBeta(dedxMDownObj->dEdx(),!isData) + (1/tof->inverseBeta()))*0.5 ) ;
               if(dedxMDownObj) MassDownComb = MassDown;
               if(tof)MassDownComb=GetMassFromBeta(track->p(),(1/tof->inverseBeta()));

               bool PassNonTrivialSelection=false;

               //loop on all possible selection (one of them, the optimal one, will be used later)
               for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
                  //Full Selection
		 if(isMC)PassSelection   (hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, MCTrPlots);
		 if(    !PassSelection   (hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, SamplePlots, false, isSignal?genColl[ClosestGen].p()/genColl[ClosestGen].energy():-1))continue;

                  if(CutIndex!=0)PassNonTrivialSelection=true;
                  HSCPTk[CutIndex] = true;
                  HSCPTk_SystHUp[CutIndex] = true;
                  HSCPTk_SystHDown[CutIndex] = true;

                  if(Mass>MaxMass[CutIndex]) MaxMass[CutIndex]=Mass;
                  if(MassUp>MaxMass_SystHUp[CutIndex]) MaxMass_SystHUp[CutIndex]=Mass;
                  if(MassDown>MaxMass_SystHDown[CutIndex]) MaxMass_SystHDown[CutIndex]=Mass;

                  //Fill Mass Histograms
                  if(isMC)MCTrPlots->Mass->Fill(CutIndex, Mass,Event_Weight);
                  SamplePlots      ->Mass->Fill(CutIndex, Mass,Event_Weight);
                  if(tof){
                  if(isMC)MCTrPlots->MassTOF->Fill(CutIndex, MassTOF, Event_Weight);
                     SamplePlots   ->MassTOF->Fill(CutIndex, MassTOF, Event_Weight);
                  }
                  if(isMC)MCTrPlots->MassComb->Fill(CutIndex, MassComb, Event_Weight);
                  SamplePlots      ->MassComb->Fill(CutIndex, MassComb, Event_Weight);


                  //Fill Mass Histograms for different Ih syst

                           SamplePlots->Mass_SystHUp  ->Fill(CutIndex, MassUp,Event_Weight);
                           SamplePlots->Mass_SystHDown->Fill(CutIndex, MassDown,Event_Weight);
                           if(tof){
                              SamplePlots->MassTOF_SystH ->Fill(CutIndex, MassTOF , Event_Weight);
                           }
                           SamplePlots->MassComb_SystHUp  ->Fill(CutIndex, MassUpComb, Event_Weight);
                           SamplePlots->MassComb_SystHDown->Fill(CutIndex, MassDownComb, Event_Weight);
 

               } //end of Cut loop
	       //              if(PassNonTrivialSelection) stPlots_FillTree(SamplePlots, ev.eventAuxiliary().run(),ev.eventAuxiliary().event(), c, track->pt(), dedxSObj ? dedxSObj->dEdx() : -1, tof ? tof->inverseBeta() : -1, Mass, TreeDZ, TreeDXY, OpenAngle, track->eta(), track->phi(), -1);
	       double Ick2=0;  if(dedxMObj) Ick2=GetIck(dedxMObj->dEdx(),isMC);
	       int nomh= 0;nomh = track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::MISSING_INNER_HITS) + track->hitPattern().trackerLayersWithoutMeasurement(reco::HitPattern::TRACK_HITS);
	       double fovhd = track->found()<=0?-1:track->found() / float(track->found() + nomh);
	       unsigned int nom=0; if(dedxSObj) nom=dedxSObj->numberOfMeasurements();

/*
if(isSignal)std::cout<<"PD "<<genColl[ClosestGen].pdgId()<< "  charge "<<genColl[ClosestGen].charge()<<"  p" <<genColl[ClosestGen].p() <<"  mass "<<genColl[ClosestGen].mass()
<<"  eta "<< genColl[ClosestGen].eta()
<<"  phi "<< genColl[ClosestGen].phi() 
<<"  pt "<< genColl[ClosestGen].pt() 
<<" charge "<<track->charge()
<<" nom "<<nom
<<std::endl;
*/
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


	      // if(PassNonTrivialSelection)
//fill ntuple without any preselection
//
//if(PassNonTrivialSelection)  

if(PassNonTrivialSelection||(dedxSObj && dedxSObj->dEdx()> 0. && track->pt()>60.))
stPlots_FillTree(SamplePlots, TrigInfo, ev.eventAuxiliary().run(),ev.eventAuxiliary().event(),ev.eventAuxiliary().luminosityBlock(), c, track->charge(), track->pt(),track->ptError(), dedxSObj ? dedxSObj->dEdx() : -1,dedxSObj ? dedxMObj->dEdx() : -1,dedxMObj ? Ick2 : -99, tof ? tof->inverseBeta() : -1, Mass, TreeDZ, TreeDXY, OpenAngle, track->eta(), track->phi(), track->found(), track->hitPattern().numberOfValidPixelHits(), track->validFraction(), nomh,fovhd, nom, weight,genid,gencharge,genmass,genpt,geneta,genphi,-1);
            }// end of Track Loop

            //save event dependent information thanks to the bookkeeping
            for(unsigned int CutIndex=0;CutIndex<CutPt.size();CutIndex++){
              if(HSCPTk[CutIndex]){
                 SamplePlots->HSCPE             ->Fill(CutIndex,Event_Weight);
                 SamplePlots->MaxEventMass      ->Fill(CutIndex,MaxMass[CutIndex], Event_Weight);
                 if(isMC){
                 MCTrPlots->HSCPE               ->Fill(CutIndex,Event_Weight);
                 MCTrPlots->MaxEventMass        ->Fill(CutIndex,MaxMass[CutIndex], Event_Weight);
                 }
              }
              if(HSCPTk_SystP[CutIndex]){
                 SamplePlots->HSCPE_SystP       ->Fill(CutIndex,Event_Weight);
                 SamplePlots->MaxEventMass_SystP->Fill(CutIndex,MaxMass_SystP[CutIndex], Event_Weight);
              }
              if(HSCPTk_SystI[CutIndex]){
                 SamplePlots->HSCPE_SystI       ->Fill(CutIndex,Event_Weight);
                 SamplePlots->MaxEventMass_SystI->Fill(CutIndex,MaxMass_SystI[CutIndex], Event_Weight);
              }
              if(HSCPTk_SystM[CutIndex]){
                 SamplePlots->HSCPE_SystM       ->Fill(CutIndex,Event_Weight);
                 SamplePlots->MaxEventMass_SystM->Fill(CutIndex,MaxMass_SystM[CutIndex], Event_Weight);
              }
              if(HSCPTk_SystT[CutIndex]){
                 SamplePlots->HSCPE_SystT       ->Fill(CutIndex,Event_Weight);
                 SamplePlots->MaxEventMass_SystT->Fill(CutIndex,MaxMass_SystT[CutIndex], Event_Weight);
              }
              if(HSCPTk_SystPU[CutIndex]){
                 SamplePlots->HSCPE_SystPU       ->Fill(CutIndex,Event_Weight*PUSystFactor);
                 SamplePlots->MaxEventMass_SystPU->Fill(CutIndex,MaxMass_SystPU[CutIndex], Event_Weight*PUSystFactor);
              }
              if(HSCPTk_SystHUp[CutIndex]){
                 SamplePlots->HSCPE_SystHUp     ->Fill(CutIndex,Event_Weight);
                 SamplePlots->MaxEventMass_SystHUp   ->Fill(CutIndex,MaxMass_SystHUp   [CutIndex], Event_Weight);
              }
              if(HSCPTk_SystHDown[CutIndex]){
                 SamplePlots->HSCPE_SystHDown   ->Fill(CutIndex,Event_Weight);
                 SamplePlots->MaxEventMass_SystHDown ->Fill(CutIndex,MaxMass_SystHDown [CutIndex], Event_Weight);
              }

           }
         }printf("\n");// end of Event Loop
      }//end of period loop
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

      stPlots_Clear(SamplePlots, true);
      if(isMC)stPlots_Clear(MCTrPlots, true);
   }// end of sample loop
   delete RNG;
}

void InitHistos(stPlots* st){   
   //Initialization of variables that are common to all samples
   HistoFile->cd();
   HCuts_Pt  = new TProfile("HCuts_Pt" ,"HCuts_Pt" ,CutPt.size(),0,CutPt.size());
   HCuts_I   = new TProfile("HCuts_I"  ,"HCuts_I"  ,CutPt.size(),0,CutPt.size());
   HCuts_TOF = new TProfile("HCuts_TOF","HCuts_TOF",CutPt.size(),0,CutPt.size());
   for(unsigned int i=0;i<CutPt.size();i++){  HCuts_Pt->Fill(i,CutPt[i]);     HCuts_I->Fill(i,CutI[i]);    HCuts_TOF->Fill(i,CutTOF[i]); }

   HCuts_Pt_Flip  = new TProfile("HCuts_Pt_Flip" ,"HCuts_Pt_Flip" ,CutPt_Flip.size(),0,CutPt_Flip.size());
   HCuts_I_Flip   = new TProfile("HCuts_I_Flip"  ,"HCuts_I_Flip"  ,CutPt_Flip.size(),0,CutPt_Flip.size());
   HCuts_TOF_Flip = new TProfile("HCuts_TOF_Flip","HCuts_TOF_Flip",CutPt_Flip.size(),0,CutPt_Flip.size());
   for(unsigned int i=0;i<CutPt_Flip.size();i++){  HCuts_Pt_Flip->Fill(i,CutPt_Flip[i]);     HCuts_I_Flip->Fill(i,CutI_Flip[i]);    HCuts_TOF_Flip->Fill(i,CutTOF_Flip[i]); }
}

// code needed for the evaluation of the systematics related to pt measurement
double RescaledPt(const double& pt, const double& eta, const double& phi, const int& charge)
{
  if(TypeMode!=3) {
    double newInvPt = 1/pt+0.000236-0.000135*pow(eta,2)+charge*0.000282*TMath::Sin(phi-1.337);
    return 1/newInvPt;
  }
  else {
    double newInvPt = (1./pt)*1.1;
    return 1/newInvPt;
  }
}

double SegSep(const susybsm::HSCParticle& hscp, const fwlite::ChainEvent& ev, double& minPhi, double& minEta) {
  if(TypeMode!=3)return -1;

  reco::MuonRef muon = hscp.muonRef();
  if(muon.isNull()) return false;
  reco::TrackRef  track = muon->standAloneMuon();
  if(track.isNull())return false;

  fwlite::Handle<MuonSegmentCollection> SegCollHandle;
  SegCollHandle.getByLabel(ev, "MuonSegmentProducer");
  if(!SegCollHandle.isValid()){printf("Segment Collection Not Found\n"); return -1;}
  MuonSegmentCollection SegCollection = *SegCollHandle;

  double minDr=10;
  minPhi=10;
  minEta=10;

  //Look for segment on opposite side of detector from track
  for (MuonSegmentCollection::const_iterator segment = SegCollection.begin(); segment!=SegCollection.end();++segment) {  
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

//Counts the number of muon stations used in track fit only counting DT and CSC stations.
int  muonStations(const reco::HitPattern& hitPattern) {
  int stations[4] = { 0,0,0,0 };
  for (int i=0; i<hitPattern.numberOfHits(reco::HitPattern::HitCategory::TRACK_HITS); i++) {
    uint32_t pattern = hitPattern.getHitPattern(reco::HitPattern::HitCategory::TRACK_HITS, i );
    if(pattern == 0) break;
    if(hitPattern.muonHitFilter(pattern) && (int(hitPattern.getSubStructure(pattern)) == 1 || int(hitPattern.getSubStructure(pattern)) == 2) && hitPattern.getHitType(pattern) == 0){
      stations[hitPattern.getMuonStation(pattern)-1] = 1;
    }
  }
  return stations[0]+stations[1]+stations[2]+stations[3];

}

double scaleFactor(double eta) {
  double etaBins[15]   = {-2.1, -1.8, -1.5, -1.2, -0.9, -0.6, -0.3, 0.0 , 0.3 , 0.6 , 0.9 , 1.2 ,1.5 , 1.8 , 2.1 };
  double scaleBins[15] = {0,    0.97, 1.06, 1.00, 0.89, 0.91, 0.93, 0.93, 0.92, 0.92, 0.91, 0.89,1.00, 1.06, 0.99};
  for (int i=0; i<15; i++) if(eta<etaBins[i]) return scaleBins[i];
  return 0;
}
