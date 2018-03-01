#include <exception>
#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDCacheFile.h"
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
#include "TPaveText.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TDirectory.h"
#include "TCanvas.h"
#include "TProfile.h"
#include "TPaveText.h"
#include "TRandom.h"

namespace reco { class Vertex; class Track; class GenParticle; class DeDxData; class MuonTimeExtra;}
namespace susybsm { class HSCParticle; class HSCPIsolation;}
namespace fwlite { class ChainEvent;}
namespace trigger { class TriggerEvent;}
namespace edm {class TriggerResults; class TriggerResultsByName; class InputTag; class LumiReWeighting;}
namespace reweight{class PoissonMeanShifter;}

#if !defined(__CINT__) && !defined(__MAKECINT__)
#include "FWCore/FWLite/interface/FWLiteEnabler.h"
#include "DataFormats/FWLite/interface/Handle.h"
#include "DataFormats/FWLite/interface/Event.h"
#include "DataFormats/FWLite/interface/Run.h"
#include "DataFormats/FWLite/interface/ChainEvent.h"
#include "DataFormats/Common/interface/MergeableCounter.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/HepMCCandidate/interface/GenParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCParticle.h"
#include "AnalysisDataFormats/SUSYBSMObjects/interface/HSCPIsolation.h"
#include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h"
#include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"
#include "DataFormats/Common/interface/TriggerResults.h"
#include "PhysicsTools/Utilities/interface/LumiReWeighting.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"
#include "SimDataFormats/GeneratorProducts/interface/GenRunInfoProduct.h"

using namespace fwlite;
using namespace reco;
using namespace susybsm;
using namespace std;
using namespace edm;
using namespace trigger;

#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"

#endif

/*
void computeSimpleLimits(string outpath, string ChannelName, string SignalName, double Obs, double Pred, double PredRelErr, double SignEff, double SignEffStat, double SignalEffUnc, double Lumi, double LumiUnc=1.044, double* limits=NULL){

   double NSign = SignEff*Lumi;
   double Scale = 1.0;
   while(NSign>100){NSign/=10.0; Scale*=10.0;}

   FILE* pFile = fopen(outpath.c_str(), "w");
   fprintf(pFile, "imax 1\n");
   fprintf(pFile, "jmax *\n");
   fprintf(pFile, "kmax *\n");
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "bin %s\n",ChannelName.c_str());
   fprintf(pFile, "Observation %f\n",Obs);
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "bin      %s %s\n",ChannelName.c_str(), ChannelName.c_str());
   fprintf(pFile, "process  %s pred\n",SignalName.c_str());
   fprintf(pFile, "process  0 1\n");
   fprintf(pFile, "rate    %f %f\n",NSign,std::max(1E-4, Pred) );  //if Pred<1E-4 we have troubles when merging datacards
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "%35s    %6s %5.3f     1.0  \n","Lumi" , "lnN", LumiUnc);
   fprintf(pFile, "%35s    %6s -         %5.3f\n",(ChannelName+"systP").c_str(), "lnN", PredRelErr);
   fprintf(pFile, "%35s    %6s %5.3f     -    \n",(ChannelName+"systS").c_str(), "lnN", SignalEffUnc);
   fprintf(pFile, "%35s    %6s %5.3f     -    \n",(ChannelName+"statS").c_str(), "lnN", std::min(SignEffStat,2.0));
   fclose(pFile);

   string JobName = SignalName;
   string massStr = "0";

   //prepare and run the script that will run the external "combine" tool from the Higgs group
   //If very low background range too small, set limit at 0.001.  Only affects scanning range not final limit
   if(Pred<0.001) Pred=0.001;
   char rangeStr[255];sprintf(rangeStr," --rMin %f --rMax %f ", 0.0f, 2*(3*sqrt(Pred)/NSign) );
   string CodeToExecute = "cd /tmp/;";
   CodeToExecute += "combine -M Asymptotic        -n " + JobName + " -m " + massStr + rangeStr + " " + outpath + " &> " + outpath + ".log;";

   system(CodeToExecute.c_str());

   //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
   //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
   TFile* file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+".Asymptotic.mH"+massStr+".root").c_str());
   if(!file || file->IsZombie())return;
   TTree* tree = (TTree*)file->Get("limit");
   if(!tree)return;
   double Tmass, Tlimit, TlimitErr; float TquantExp;
   tree->GetBranch("mh"              )->SetAddress(&Tmass    );
   tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
   tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
   tree->GetBranch("limitErr"        )->SetAddress(&TlimitErr);
   tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
   for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
     tree->GetEntry(ientry);
           if(TquantExp==0.025f){ limits[0] = Tlimit/Scale;
     }else if(TquantExp==0.160f){ limits[1] = Tlimit/Scale;
     }else if(TquantExp==0.500f){ limits[2] = Tlimit/Scale;
     }else if(TquantExp==0.840f){ limits[3] = Tlimit/Scale;
     }else if(TquantExp==0.975f){ limits[4] = Tlimit/Scale;
     }else if(TquantExp==-1    ){ limits[5] = Tlimit/Scale; //will be overwritten afterward
     }else{printf("Quantil %f unused by the analysis --> check the code\n", TquantExp);
     }
   }
   file->Close();
}

*/


void computeSimpleLimits(string outpath, string ChannelName, string SignalName, double Obs, double Pred, double PredRelErr, double SignEff, double SignEffStat, double SignalEffUnc, double Lumi, double LumiUnc=1.044, double* limits=NULL){

   double NSign = SignEff*Lumi;
   double Scale = 1.0;
   while(NSign>100){NSign/=10.0; Scale*=10.0;}

   FILE* pFile = fopen(outpath.c_str(), "w");
   fprintf(pFile, "imax 1\n");
   fprintf(pFile, "jmax *\n");
   fprintf(pFile, "kmax *\n");
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "bin %s\n",ChannelName.c_str());
   fprintf(pFile, "Observation %f\n",Obs);
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "bin      %s %s\n",ChannelName.c_str(), ChannelName.c_str());
   fprintf(pFile, "process  %s pred\n",SignalName.c_str());
   fprintf(pFile, "process  0 1\n");
   fprintf(pFile, "rate    %f %f\n",NSign,std::max(1E-4, Pred) );  //if Pred<1E-4 we have troubles when merging datacards
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "%35s    %6s %5.3f     1.0  \n","Lumi" , "lnN", LumiUnc);
   fprintf(pFile, "%35s    %6s -         %5.3f\n",(ChannelName+"systP").c_str(), "lnN", PredRelErr);
   fprintf(pFile, "%35s    %6s %5.3f     -    \n",(ChannelName+"systS").c_str(), "lnN", SignalEffUnc);
   fprintf(pFile, "%35s    %6s %5.3f     -    \n",(ChannelName+"statS").c_str(), "lnN", std::min(SignEffStat,2.0));
   fclose(pFile);

   string JobName = SignalName;
   string massStr = "0";

   //prepare and run the script that will run the external "combine" tool from the Higgs group
   //If very low background range too small, set limit at 0.001.  Only affects scanning range not final limit
   if(Pred<0.001) Pred=0.001;
   char rangeStr[255];sprintf(rangeStr," --rMin %f --rMax %f ", 0.0f, 2*(3*sqrt(Pred)/NSign) );
   string CodeToExecute = "cd /tmp/;";
   CodeToExecute += "combine -M Asymptotic        -n " + JobName + " -m " + massStr + rangeStr + " " + outpath + " &> " + outpath + ".log;";

   printf("Execute %s\n", CodeToExecute.c_str());fflush(stdout);
   system(CodeToExecute.c_str());
   printf("Done\n");

   //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
   //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
   TFile* file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+".Asymptotic.mH"+massStr+".root").c_str());
   if(!file || file->IsZombie())return;
   TTree* tree = (TTree*)file->Get("limit");
   if(!tree)return;
   double Tmass, Tlimit, TlimitErr; float TquantExp;
   tree->GetBranch("mh"              )->SetAddress(&Tmass    );
   tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
   tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
   tree->GetBranch("limitErr"        )->SetAddress(&TlimitErr);
   tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
   for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
     tree->GetEntry(ientry);
           if(TquantExp==0.025f){ limits[0] = Tlimit/Scale;
     }else if(TquantExp==0.160f){ limits[1] = Tlimit/Scale;
     }else if(TquantExp==0.500f){ limits[2] = Tlimit/Scale;
     }else if(TquantExp==0.840f){ limits[3] = Tlimit/Scale;
     }else if(TquantExp==0.975f){ limits[4] = Tlimit/Scale;
     }else if(TquantExp==-1    ){ limits[5] = Tlimit/Scale; //will be overwritten afterward
     }else{printf("Quantil %f unused by the analysis --> check the code\n", TquantExp);
     }
   }
   file->Close();

   return;


   sprintf(rangeStr," ");//--rMax 1E20 ");

   printf("RunCLS for observedLimit\n");fflush(stdout);
   double FullHybridSF = limits[5];

      //RUN FULL HYBRID CLS LIMIT (just for observed limit so far, because it is very slow for expected limits --> should be updated --> FIXME)
      CodeToExecute = "cd /tmp/;";
      CodeToExecute += "combine -M HybridNew -n " + JobName + " -m " + massStr + rangeStr + " " + outpath + ";";
//      CodeToExecute += "cd $OLDPWD; cp /tmp/shape_" + JobName + ".* " + InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/." + ";";
      system(CodeToExecute.c_str());

      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+".HybridNew.mH"+massStr+".root").c_str());
      if(!file || file->IsZombie())return;
      tree = (TTree*)file->Get("limit");
      if(!tree)return;
      tree->GetBranch("mh"              )->SetAddress(&Tmass    );
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("limitErr"        )->SetAddress(&TlimitErr);
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        if(TquantExp==-1    ){ limits[5] = Tlimit/Scale;
        }else{printf("Quantil %f unused by the analysis --> check the code\n", TquantExp);
        }
      }
      file->Close();

      printf("Done\n");fflush(stdout);

      FullHybridSF = limits[5]/FullHybridSF;
      limits[0] *= FullHybridSF;
      limits[1] *= FullHybridSF;
      limits[2] *= FullHybridSF;
      limits[3] *= FullHybridSF;
      limits[4] *= FullHybridSF;

      //uncomment to return before computing the full hybrid limits
      return;

      //Create grid to find expected limits in Cls
      //Number of different signal strengths to try
      int gridPoints=30;
      //Normalize to 10/fb
      double Down2 = 0.5*limits[0]*Scale;
      double Up2 = 1.3*limits[4]*Scale;
      double Step=(Up2-Down2)/gridPoints;

      CodeToExecute = "cd /tmp/;";
      for (int i=0; i<gridPoints+1; i++) {
        printf("RunCLS for expectedLimit %i/%i\n", i, gridPoints);fflush(stdout);

	char Seed[1024];
	sprintf(Seed,"%i",i);
        char PointStr[1024];

	double Point = Down2 + i*Step;
        sprintf(PointStr,"%6.8f",Point);
	//Don't include mass string here or else it won't work
//	CodeToExecute += "combine " + outpath + " -M HybridNew --freq --fork 1 -T 500 --clsAcc 0 -n " + JobName +              " --saveHybridResult --saveToys -s " + Seed + " -i 8 --rMax 1E20 --singlePoint " + PointStr + " >> " + outpath + "Exp.log;";
      CodeToExecute += "combine " + outpath + " -M HybridNew --freq --fork 1 -T 500 --clsAcc 0 -n " + JobName +              " --saveHybridResult --saveToys -s " + Seed + " -i 8 " + rangeStr + " --singlePoint " + PointStr + ";";// >> " + outpath + "Exp.log;";
      }

      CodeToExecute += "hadd -f higgsCombine"+JobName+".HybridNew.mH"+massStr+"grid.root higgsCombine"+JobName+".HybridNew.mH120.*.root >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine " + outpath + " -M HybridNew --grid=higgsCombine"+JobName+".HybridNew.mH"+massStr+"grid.root -n " + JobName + "Expected --expectedFromGrid 0.5 >>  " + outpath + "Exp.log;";
      CodeToExecute += "combine " + outpath + " -M HybridNew --grid=higgsCombine"+JobName+".HybridNew.mH"+massStr+"grid.root -n " + JobName + "Expected --expectedFromGrid 0.16 >>  " + outpath +  "Exp.log;";
      CodeToExecute += "combine " + outpath + " -M HybridNew --grid=higgsCombine"+JobName+".HybridNew.mH"+massStr+"grid.root -n " + JobName + "Expected --expectedFromGrid 0.84 >>  " + outpath +  "Exp.log;";
      CodeToExecute += "combine " + outpath + " -M HybridNew --grid=higgsCombine"+JobName+".HybridNew.mH"+massStr+"grid.root -n " + JobName + "Expected --expectedFromGrid 0.025 >>  " + outpath +  "Exp.log;";
      CodeToExecute += "combine " + outpath + " -M HybridNew --grid=higgsCombine"+JobName+".HybridNew.mH"+massStr+"grid.root -n " + JobName + "Expected --expectedFromGrid 0.975 >>  " + outpath +  "Exp.log;";
//      CodeToExecute += "cd $OLDPWD; cp /tmp/shape_" + JobName + "Exp.* " + InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/." + ";";

      printf("RunCLS for expectedLimit execute\n");fflush(stdout);
      system(CodeToExecute.c_str());
      printf("Done\n");fflush(stdout);


      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.500.root").c_str());
      if(!file || file->IsZombie()){printf("Can't fine file %s", (string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.500.root").c_str()); return;}
      tree = (TTree*)file->Get("limit");
      if(!tree){printf("Can't find tree named limit"); return;}
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        if(TquantExp==0.500f){ limits[2] = Tlimit/Scale;;
        }else{printf("Quantil %f should be 0.5 --> check the code\n", TquantExp);
        }
      }
      file->Close();
      
      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.160.root").c_str());
      if(!file || file->IsZombie()){printf("Can't fine file %s", (string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.160.root").c_str()); return;}
      tree = (TTree*)file->Get("limit");
      if(!tree){printf("Can't find tree named limit"); return;}
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        if(TquantExp==0.160f){ limits[1] = Tlimit/Scale;
        }else{printf("Quantil %f should be 0.16 --> check the code\n", TquantExp);
        }
      }
      file->Close();

      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.840.root").c_str());
      if(!file || file->IsZombie()){printf("Can't fine file %s", (string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.840.root").c_str()); return;}
      tree = (TTree*)file->Get("limit");
      if(!tree){printf("Can't find tree named limit"); return;}
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        if(TquantExp==0.840f){ limits[3] = Tlimit/Scale;
        }else{printf("Quantil %f should be 0.84 --> check the code\n", TquantExp);
        }
      }
      file->Close();

      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.025.root").c_str());
      if(!file || file->IsZombie()){printf("Can't fine file %s", (string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.025.root").c_str()); return;}
      tree = (TTree*)file->Get("limit");
      if(!tree){printf("Can't find tree named limit"); return;}
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        if(TquantExp==0.025f){ limits[0] = Tlimit/Scale;
        }else{printf("Quantil %f should be 0.025 --> check the code\n", TquantExp);
        }
      }
      file->Close();

      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.975.root").c_str());
      if(!file || file->IsZombie()){printf("Can't fine file %s", (string("/tmp/")+"higgsCombine"+JobName+"Expected.HybridNew.mH120.quant0.975.root").c_str()); return;}
      tree = (TTree*)file->Get("limit");
      if(!tree){printf("Can't find tree named limit"); return;}
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        if(TquantExp==0.975f){ limits[4] = Tlimit/Scale;
        }else{printf("Quantil %f should be 0.975 --> check the code\n", TquantExp);
        }
      }
      file->Close();
}



double getXSecFromFile(string url){
   TFile* file = TFile::Open(url.c_str());      
   fwlite::Run run( file );
   for(run.toBegin(); !run.atEnd(); ++run){
      fwlite::Handle<GenRunInfoProduct> genInfo;
      genInfo.getByLabel(run,"generator");
      if(!genInfo.isValid()){printf("Invalid genInfoH to get XSec\n");continue;}
      return genInfo->crossSection();
   }
   return 0;
}



bool isIsolated(const std::vector<reco::GenParticle>& genColl, const reco::GenParticle& cand, unsigned int currentCanddiate){
      double caloIso=0.0;
      double tkIso=0.0;
      for(unsigned int g=0;g<genColl.size();g++){
         if(genColl[g].status()!=1)continue;
//         if(genColl[g].pt()<0.5)continue;
         int AbsPdg=abs(genColl[g].pdgId());
         if(g==currentCanddiate)continue; //skip HSCPs
//         if(AbsPdg>1000000 || AbsPdg==17)continue; //skip HSCPs
         if(AbsPdg==12 || AbsPdg==14 || AbsPdg==16)continue; //skip neutrinos
         if(deltaR(cand.eta(), cand.phi(), genColl[g].eta(), genColl[g].phi())>0.3)continue;
//         printf("pdg = %i status=%i  charge=%i p=%f pt=%f\n", genColl[g].pdgId(), genColl[g].status(), genColl[g].charge(), genColl[g].p(), genColl[g].pt());
         if(abs(genColl[g].charge())!=0)tkIso+=genColl[g].pt();
         caloIso+=genColl[g].energy();
      }
//      printf("CaloIso = %f -> %f tkIso = %f\n", caloIso, caloIso/cand.p(), tkIso);
      if(caloIso/cand.p()>0.3)return false;
      if(tkIso>50)return false;
      return true;
}



void ModelIndependent_Acceptance(string MODE="COMPILE", string fileurl="", double ctau=-1)
{
  if(MODE=="COMPILE") return;

   Event_Weight = 1;
   MaxEntry = -1;


   system("mkdir -p pictures");

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

   TRandom randomGenerator;

   TFile* InputFile = new TFile("pictures/Histos.root");
   unsigned int nM = 6;

   string histoNames = "Norm_Beta_Triggered";
   TH3F* histoTrigger = (TH3F*)GetObjectFromPath(InputFile, histoNames);
   histoNames = "Norm_Beta_Preselected";
   TH3F* histoPreselected = (TH3F*)GetObjectFromPath(InputFile, histoNames);
   TH3F** histoOffline = new TH3F*[nM];
   for(unsigned int Mi=0;Mi<nM;Mi++){
      char tmpHistoName[256];
      sprintf(tmpHistoName, "Norm_Beta_SelectedM%i", Mi);
      histoNames = tmpHistoName;
      histoOffline[Mi] = (TH3F*)GetObjectFromPath(InputFile, histoNames);
   }

   vector<string> InputFiles;
   InputFiles.push_back(fileurl);

   for(unsigned int s=0;s<InputFiles.size();s++){
      int CutIndex=0;
      TypeMode = 0;
      GlobalMaxEta   = 2.1;
      GlobalMinPt    = 70;
      system("mkdir pictures/");

      double NEvents = 0;
      double NTEvents = 0, NTEventsErr = 0;
      double NPSEvents = 0, NPSEventsErr = 0;
      double* NSEvents = new double[nM];
      double* NSEventsErr = new double[nM];
      for(unsigned int Mi=0;Mi<nM;Mi++){NSEvents[Mi]=0; NSEventsErr[Mi]=0;}

      printf("%s\n",InputFiles[s].c_str());
      double ThXSec = getXSecFromFile(InputFiles[s]);



      vector<string> DataFileName;
      DataFileName.push_back(InputFiles[s]);

      fwlite::ChainEvent treeS(DataFileName);
      printf("Progressing Bar              :0%%       20%%       40%%       60%%       80%%       100%%\n");
      printf("Looping on Tree              :");
      int TreeStep = treeS.size()/50;if(TreeStep==0)TreeStep=1;
      for(Long64_t e=0;e<treeS.size();e++){
         treeS.to(e); 
         if(e%TreeStep==0){printf(".");fflush(stdout);}
         fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
         genCollHandle.getByLabel(treeS, "genParticles");
         if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound\n");continue;}
         std::vector<reco::GenParticle> genColl = *genCollHandle;

         NEvents++;
         bool onlyCharged = true;
         double ProbS1 = 1;
         double ProbS2 = 1;
         double ProbT1 = 0,  ProbErrT1=0;
         double ProbT2 = 0,  ProbErrT2=0;
         double ProbPS1 = 0,  ProbErrPS1=0;
         double ProbPS2 = 0,  ProbErrPS2=0;
         double* ProbO1 = new double[nM];
         double* ProbErrO1 = new double[nM];
         double* ProbO2 = new double[nM];
         double* ProbErrO2 = new double[nM];
         int H=0;
         for(unsigned int Mi=0;Mi<nM;Mi++){ProbO1[Mi]=0; ProbErrO1[Mi]=0; ProbO2[Mi]=0; ProbErrO2[Mi]=0;}
         int NHSCP = 0;
         for(unsigned int g=0;g<genColl.size();g++){
            if(genColl[g].pt()<5)continue;
            if(genColl[g].status()!=1)continue;
            int AbsPdg=abs(genColl[g].pdgId());
            if(AbsPdg<1000000 && AbsPdg!=17)continue;
            if(onlyCharged && (AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 || AbsPdg==1093214 || AbsPdg==1093324))continue; //Skip neutral gluino RHadrons
            if(onlyCharged && (AbsPdg==1000622 || AbsPdg==1000642 || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 || AbsPdg==1006333))continue;  //skip neutral stop RHadrons

            if(genColl[g].pt()<40 || fabs(genColl[g].eta())>2.1)continue;    //why?!

            if(!isIsolated(genColl, genColl[g], g))continue;

            NHSCP++;
            
            //printf("%f > %E=exp(%E)\n", randomGenerator.Uniform(1.0), exp(-genColl[g].mass()*10.0*(1.0/ctau)/genColl[g].p()), -genColl[g].mass()/genColl[g].p());

            //does it leaves long enough (10m)? --> skip randomly particles depending on their survival probability
//            double x = 10.0;
            double x = 11.0;
            if(fabs(genColl[g].eta())<0.3) x = 7.0;
            if(fabs(genColl[g].eta())<0.5) x = 8.0;
//            if(fabs(genColl[g].eta())<0.8) x = 9.5;
            if(fabs(genColl[g].eta())<1.1) x = 10.0;
//            if(ctau>0 && randomGenerator.Uniform(1.0)>exp(-genColl[g].mass()*x*(1.0/ctau)/genColl[g].p()))continue;
            double P = ctau>0?exp(-genColl[g].mass()*x*(1.0/ctau)/genColl[g].p()):1.0;
            
            int BinX = histoTrigger->GetXaxis()->FindBin(genColl[g].pt());
            int BinY = histoTrigger->GetYaxis()->FindBin(genColl[g].p()/genColl[g].energy()); 
            int BinZ = histoTrigger->GetZaxis()->FindBin(fabs(genColl[g].eta()) );

//            if(H==0){ ProbS1 = P;
//            }else{    ProbS2 = P;
//            }

            if(H==0){       ProbT1=P*histoTrigger->GetBinContent(BinX, BinY, BinZ);  ProbErrT1=P*histoTrigger->GetBinError(BinX, BinY, BinZ);
            }else{          ProbT2=P*histoTrigger->GetBinContent(BinX, BinY, BinZ);  ProbErrT2=P*histoTrigger->GetBinError(BinX, BinY, BinZ);
            }

            if(H==0){       ProbPS1=histoPreselected->GetBinContent(BinX, BinY, BinZ);  ProbErrPS1=histoPreselected->GetBinError(BinX, BinY, BinZ);
            }else{          ProbPS2=histoPreselected->GetBinContent(BinX, BinY, BinZ);  ProbErrPS2=histoPreselected->GetBinError(BinX, BinY, BinZ);
            }

            for(unsigned int Mi=0;Mi<nM;Mi++){            
               if(genColl[g].mass()<Mi*100)continue;
               if(H==0){   ProbO1[Mi]=histoOffline[Mi]->GetBinContent(BinX, BinY, BinZ);  ProbErrO1[Mi]=histoOffline[Mi]->GetBinError(BinX, BinY, BinZ);
               }else{      ProbO2[Mi]=histoOffline[Mi]->GetBinContent(BinX, BinY, BinZ);  ProbErrO2[Mi]=histoOffline[Mi]->GetBinError(BinX, BinY, BinZ);
               }
            }
            H++;
         }
//         printf("NHSCP=%i\n", NHSCP);
         double EventTProbErr  = 0;
         double EventTProb     = 0;      
         double EventPSProbErr = 0;
         double EventPSProb    = 0;      

         EventTProb    = ProbT1 + ProbT2 - ProbT1*ProbT2;
         EventTProbErr = pow(ProbErrT1*(1-ProbT2),2) + pow(ProbErrT2*(1-ProbT1),2);
         NTEvents     += EventTProb;
         NTEventsErr  += EventTProbErr;

         EventPSProb    = ProbPS1 + ProbPS2 - ProbPS1*ProbPS2;
         EventPSProbErr = pow(ProbErrPS1*(1-ProbPS2),2) + pow(ProbErrPS2*(1-ProbPS1),2);
         NPSEvents     += EventPSProb;
         NPSEventsErr  += EventPSProbErr;

         for(unsigned int Mi=0;Mi<nM;Mi++){
            double EventProb    = EventTProb;
            double EventProbErr = EventTProbErr;
            EventProb        *= (ProbO1[Mi] + ProbO2[Mi] - ProbO1[Mi]*ProbO2[Mi]);
            EventProbErr     += pow(ProbErrO1[Mi]*(1-ProbO2[Mi]),2) + pow(ProbErrO2[Mi]*(1-ProbO1[Mi]),2);
            NSEvents[Mi]     += EventProb;
            NSEventsErr[Mi]  += EventProbErr;
         }
        
      }printf("\n");// end of Event Loop

      char TxtOutput[20000]; sprintf(TxtOutput,"");
      for(unsigned int Mi=0;Mi<nM;Mi++){
         double NPred=-1, NPredErr=-1, NData=-1;
              if(Mi==0){NPred=44.0; NPredErr=9.0; NData=42;}
         else if(Mi==1){NPred= 5.6; NPredErr=1.1; NData= 7;}
         else if(Mi==2){NPred= 0.56;NPredErr=0.11;NData= 0;}
         else          {NPred= 0.02;NPredErr=0.004;NData= 0;}

         double limits[] = {-1, -1, -1, -1, -1, -1};
         printf("Compute the limits for M=%i\n", Mi);fflush(stdout);
         computeSimpleLimits(string("/tmp/combine_")+MODE+".dat", "TkTOF", MODE, NData, NPred, 1.0 + NPredErr/NPred, NSEvents[Mi]/NEvents, 1.0 + sqrt(pow(sqrt(NSEventsErr[Mi])/NEvents,2) + pow(NSEvents[Mi]*sqrt(NEvents)/pow(NEvents,2),2)), 1.40, IntegratedLuminosity8TeV, 1.044, limits);
         sprintf(TxtOutput, "%s%30s M>%3i Efficiencies: Trigger=%6.2f%%+-%6.2f%%  Presel=%6.2f%%+-%6.2f%% Offline=%6.2f%%+-%6.2f%%  Limits=%E %E %E %E %E %E  ThXSec=%E\n",TxtOutput, MODE.c_str(), Mi*100, 100.0*NTEvents/NEvents, 100.0*sqrt(pow(sqrt(NTEventsErr)/NEvents,2) + pow(NTEvents*sqrt(NEvents)/pow(NEvents,2),2)), 100.0*NPSEvents/NEvents, 100.0*sqrt(pow(sqrt(NPSEventsErr)/NEvents,2) + pow(NPSEvents*sqrt(NEvents)/pow(NEvents,2),2)), 100.0*NSEvents[Mi]/NEvents, 100.0*sqrt(pow(sqrt(NSEventsErr[Mi])/NEvents,2) + pow(NSEvents[Mi]*sqrt(NEvents)/pow(NEvents,2),2)), limits[0], limits[1], limits[2], limits[3], limits[4], limits[5], getXSecFromFile(InputFiles[s]));      
      }

      FILE* pFile = fopen((string("pictures/")+MODE+".txt").c_str(), "w");
      //FILE* pFile = fopen("testtt.tx", "w");
      fprintf(pFile , "%s", TxtOutput);
      fprintf(stdout, "%s", TxtOutput);      
      fclose(pFile);
   }
}


