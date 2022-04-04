// Original Author:  Loic Quertenmont

namespace reco {
  class Vertex;
  class Track;
  class GenParticle;
  class DeDxData;
  class MuonTimeExtra;
  class PFMET;
  class HitPattern;
}  // namespace reco
namespace susybsm {
  class HSCParticle;
  class HSCPIsolation;
  class MuonSegment;
  class HSCPDeDxInfo;
}  // namespace susybsm
namespace fwlite {
  class ChainEvent;
}
namespace trigger {
  class TriggerEvent;
}
namespace edm {
  class TriggerResults;
  class TriggerResultsByName;
  class InputTag;
  class LumiReWeighting;
}  // namespace edm
namespace reweight {
  class PoissonMeanShifter;
}

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
#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"
/////////////////////////// FUNCTION DECLARATION /////////////////////////////

// check if the event is passing trigger or not --> note that the function has two part (one for 2011 analysis and the other one for 2012)
bool PassTriggerCustom(const fwlite::ChainEvent& ev) {
  edm::TriggerResultsByName tr = ev.triggerResultsByName("MergeHLT");
  if (!tr.isValid())
    return false;

  if (tr.accept("HSCPHLTTriggerMuFilter"))
    return true;
  //	if(tr.accept("HSCPHLTTriggerPFMetFilter"))return true;

  return false;
}

void computeSimpleLimits(string outpath,
                         string ChannelName,
                         string SignalName,
                         double Obs,
                         double Pred,
                         double PredRelErr,
                         double SignEff,
                         double SignEffStat,
                         double SignalEffUnc,
                         double Lumi,
                         double LumiUnc = 1.044,
                         double* limits = NULL) {
  double NSign = SignEff * Lumi;
  double Scale = 1.0;
  while (NSign > 100) {
    NSign /= 10.0;
    Scale *= 10.0;
  }

  //    double SF = 1.0;
  //         if(Pred==44.0 ){SF = 0.95; } // SF = Hybrid/Asymptotic  M>0
  //    else if(Pred== 5.6 ){SF = 1.07; } // SF = Hybrid/Asymptotic  M>100
  //    else if(Pred== 0.56){SF = 1.38; } // SF = Hybrid/Asymptotic  M>200
  //    else if(Pred== 0.02){SF = 2.00; } // SF = Hybrid/Asymptotic  M>300
  //    else{printf("ERROR SF CANNOT BE COMPUTED.... YOU SHOULD USE FULL HYBRID LIMITS\n");}

  FILE* pFile = fopen(outpath.c_str(), "w");
  fprintf(pFile, "imax 1\n");
  fprintf(pFile, "jmax *\n");
  fprintf(pFile, "kmax *\n");
  fprintf(pFile, "-------------------------------\n");
  fprintf(pFile, "bin %s\n", ChannelName.c_str());
  fprintf(pFile, "Observation %f\n", Obs);
  fprintf(pFile, "-------------------------------\n");
  fprintf(pFile, "bin      %s %s\n", ChannelName.c_str(), ChannelName.c_str());
  fprintf(pFile, "process  %s pred\n", SignalName.c_str());
  fprintf(pFile, "process  0 1\n");
  fprintf(
      pFile, "rate    %f %f\n", NSign, std::max(1E-4, Pred));  //if Pred<1E-4 we have troubles when merging datacards
  fprintf(pFile, "-------------------------------\n");
  fprintf(pFile, "%35s    %6s %5.3f     1.0  \n", "Lumi", "lnN", LumiUnc);
  fprintf(pFile, "%35s    %6s -         %5.3f\n", (ChannelName + "systP").c_str(), "lnN", PredRelErr);
  fprintf(pFile, "%35s    %6s %5.3f     -    \n", (ChannelName + "systS").c_str(), "lnN", SignalEffUnc);
  fprintf(pFile, "%35s    %6s %5.3f     -    \n", (ChannelName + "statS").c_str(), "lnN", std::min(SignEffStat, 2.0));
  fclose(pFile);

  string JobName = SignalName;
  string massStr = "0";

  //prepare and run the script that will run the external "combine" tool from the Higgs group
  //If very low background range too small, set limit at 0.001.  Only affects scanning range not final limit
  if (Pred < 0.001)
    Pred = 0.001;
  char rangeStr[255];
  sprintf(rangeStr, " --rMin %f --rMax %f ", 0.0f, 2 * (3 * sqrt(Pred) / NSign));
  string CodeToExecute = "cd /tmp/;";
  CodeToExecute += "combine -M Asymptotic        -n " + JobName + " -m " + massStr + rangeStr + " " + outpath + " &> " +
                   outpath + ".log;";

  printf("Execute %s\n", CodeToExecute.c_str());
  fflush(stdout);
  system(CodeToExecute.c_str());
  printf("Done\n");

  //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
  //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
  TFile* file =
      TFile::Open((string("/tmp/") + "higgsCombine" + JobName + ".Asymptotic.mH" + massStr + ".root").c_str());
  if (!file || file->IsZombie())
    return;
  TTree* tree = (TTree*)file->Get("limit");
  if (!tree)
    return;
  double Tmass, Tlimit, TlimitErr;
  float TquantExp;
  tree->GetBranch("mh")->SetAddress(&Tmass);
  tree->GetBranch("limit")->SetAddress(&Tlimit);
  tree->GetBranch("limit")->SetAddress(&Tlimit);
  tree->GetBranch("limitErr")->SetAddress(&TlimitErr);
  tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
  for (int ientry = 0; ientry < tree->GetEntriesFast(); ientry++) {
    tree->GetEntry(ientry);
    if (TquantExp == 0.025f) {
      limits[0] = Tlimit / Scale;
    } else if (TquantExp == 0.160f) {
      limits[1] = Tlimit / Scale;
    } else if (TquantExp == 0.500f) {
      limits[2] = Tlimit / Scale;
    } else if (TquantExp == 0.840f) {
      limits[3] = Tlimit / Scale;
    } else if (TquantExp == 0.975f) {
      limits[4] = Tlimit / Scale;
    } else if (TquantExp == -1) {
      limits[5] = Tlimit / Scale;  //will be overwritten afterward
    } else {
      printf("Quantil %f unused by the analysis --> check the code\n", TquantExp);
    }
  }
  file->Close();

  return;

  /*
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
*/
}

void StandardAnalysis_Acceptance(string MODE = "COMPILE", int TypeMode_ = 0, double paperMassCut = 0) {
  if (MODE == "COMPILE")
    return;
  string sampleName = MODE;

  // redefine global variable dependent on the arguments given to the function
  InitdEdx(dEdxS_Label);
  TypeMode = TypeMode_;
  GlobalMaxEta = 2.1;
  GlobalMinPt = 45;
  if (TypeMode < 2) {
    GlobalMinNDOF = 0;
    GlobalMinTOF = 0;
  } else if (TypeMode == 2) {
  } else {
    printf("Code only suited for TkOnly and TkTOF Analysis, exit now\n");
    exit(0);
  }

  if (TypeMode < 2) {
    CutPt.push_back(70.0);
    CutI.push_back(0.4);
    CutTOF.push_back(-1);
  } else if (TypeMode == 2) {
    CutPt.push_back(70.0);
    CutI.push_back(0.125);
    CutTOF.push_back(1.225);
  }
  int CutIndex = 0;

  //determine the list of models that are considered
  std::vector<stSample> samples;
  GetSampleDefinition(samples, "../../AnalysisCode/Analysis_Samples.txt");
  InitBaseDirectory();
  int s = JobIdToIndex(sampleName, samples);
  if (s < 0) {
    printf("Current Sample (%s) Not Found!  Exit now\n", sampleName.c_str());
    exit(0);
  }
  printf("Analyze %s\n", sampleName.c_str());

  //	string SignalName="PPStau";
  //	string SignalPath="root://eoscms//eos/cms//store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/"+SignalName+"_8TeV_M";

  //initialize LumiReWeighting
  BgLumiMC.clear();
  TrueDist.clear();
  TrueDistSyst.clear();
  if (samples[s].Pileup == "S10") {
    for (int i = 0; i < 60; ++i)
      BgLumiMC.push_back(Pileup_MC_Summer2012[i]);
  } else {
    for (int i = 0; i < 60; ++i)
      BgLumiMC.push_back(Pileup_MC_Fall11[i]);
  }

  for (int i = 0; i < 60; ++i)
    TrueDist.push_back(TrueDist2012_f[i]);
  for (int i = 0; i < 60; ++i)
    TrueDistSyst.push_back(TrueDist2012_XSecShiftUp_f[i]);
  LumiWeightsMC = edm::LumiReWeighting(BgLumiMC, TrueDist);
  LumiWeightsMCSyst = edm::LumiReWeighting(BgLumiMC, TrueDistSyst);

  //Initialize a RandomNumberGenerator
  TRandom3* RNG = new TRandom3();

  bool isData = (samples[s].Type == 0);
  bool isMC = (samples[s].Type == 1);
  bool isSignal = (samples[s].Type >= 2);

  dEdxTemplates = loadDeDxTemplate("../../../data/Discrim_Templates_MC_2012.root");
  dEdxSF = 1.05;

  //do two loops through signal for samples with and without trigger changes.
  //load the files corresponding to this sample
  std::vector<string> FileName;
  GetInputFiles(samples[s], BaseDirectory, FileName);
  fwlite::ChainEvent ev(FileName);

  double NEvents = 0, NEventsErr = 0;
  double NTEvents = 0, NTEventsErr = 0;
  double NPSEvents = 0, NPSEventsErr = 0;
  double NSEvents = 0, NSEventsErr = 0;
  double NSEventsM[6] = {0, 0, 0, 0, 0, 0};
  double NSEventsMErr[6] = {0, 0, 0, 0, 0, 0};

  //Loop on the events
  printf("Progressing Bar                   :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Looping on Tree                   :");
  int TreeStep = ev.size() / 50;
  if (TreeStep == 0)
    TreeStep = 1;

  //loop once on the MC events to get the PU normalization
  double PUSystFactor = 1.0;  //not used
  double SampleWeight = ev.size();
  if (samples[s].Type > 0) {
    double NMCevents = 0;
    for (Long64_t ientry = 0; ientry < ev.size(); ientry++) {
      ev.to(ientry);
      NMCevents += GetPUWeight(ev, samples[s].Pileup, PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
    }
    SampleWeight /= NMCevents;
  }

  //real loop on events
  for (Long64_t ientry = 0; ientry < ev.size(); ientry++) {
    ev.to(ientry);
    if (MaxEntry > 0 && ientry > MaxEntry)
      break;
    if (ientry % TreeStep == 0) {
      printf(".");
      fflush(stdout);
    }
    Event_Weight = SampleWeight * GetPUWeight(ev, samples[s].Pileup, PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
    NEvents += Event_Weight;
    NEventsErr += Event_Weight * Event_Weight;

    //get the collection of generated Particles
    std::vector<reco::GenParticle> genColl;
    fwlite::Handle<std::vector<reco::GenParticle> > genCollHandle;
    genCollHandle.getByLabel(ev, "genParticles");
    if (!genCollHandle.isValid()) {
      printf("GenParticle Collection NotFound\n");
      continue;
    }
    genColl = *genCollHandle;
    int NChargedHSCP = HowManyChargedHSCP(genColl);

    //check if the event is passing trigger
    if (!PassTriggerCustom(ev))
      continue;
    NTEvents += Event_Weight;
    NTEventsErr += Event_Weight * Event_Weight;

    //load all event collection that will be used later on (HSCP COll, dEdx and TOF)
    fwlite::Handle<susybsm::HSCParticleCollection> hscpCollHandle;
    hscpCollHandle.getByLabel(ev, "HSCParticleProducer");
    //if(!hscpCollHandle.isValid()){printf("HSCP Collection NotFound\n");continue;}
    if (!hscpCollHandle.isValid())
      continue;
    const susybsm::HSCParticleCollection& hscpColl = *hscpCollHandle;

    fwlite::Handle<DeDxDataValueMap> dEdxSCollH;
    dEdxSCollH.getByLabel(ev, dEdxS_Label.c_str());
    if (!dEdxSCollH.isValid()) {
      printf("Invalid dEdx Selection collection\n");
      continue;
    }

    fwlite::Handle<DeDxDataValueMap> dEdxMCollH;
    dEdxMCollH.getByLabel(ev, dEdxM_Label.c_str());
    if (!dEdxMCollH.isValid()) {
      printf("Invalid dEdx Mass collection\n");
      continue;
    }

    fwlite::Handle<MuonTimeExtraMap> TOFCollH;
    TOFCollH.getByLabel(ev, "muontiming", TOF_Label.c_str());
    if (!TOFCollH.isValid()) {
      printf("Invalid TOF collection\n");
      return;
    }

    fwlite::Handle<MuonTimeExtraMap> TOFDTCollH;
    TOFDTCollH.getByLabel(ev, "muontiming", TOFdt_Label.c_str());
    if (!TOFDTCollH.isValid()) {
      printf("Invalid DT TOF collection\n");
      return;
    }

    fwlite::Handle<MuonTimeExtraMap> TOFCSCCollH;
    TOFCSCCollH.getByLabel(ev, "muontiming", TOFcsc_Label.c_str());
    if (!TOFCSCCollH.isValid()) {
      printf("Invalid CSC TOF collection\n");
      return;
    }

    bool Passed = false;
    bool PSPassed = false;
    bool PassedM[6] = {false, false, false, false, false, false};

    //loop on HSCP candidates
    for (unsigned int c = 0; c < hscpColl.size(); c++) {
      //define alias for important variable
      susybsm::HSCParticle hscp = hscpColl[c];
      reco::MuonRef muon = hscp.muonRef();
      reco::TrackRef track = hscp.trackRef();
      if (track.isNull())
        continue;

      //for signal only, make sure that the candidate is associated to a true HSCP
      int ClosestGen;
      if (isSignal && DistToHSCP(hscp, genColl, ClosestGen) > 0.03)
        continue;

      //load quantity associated to this track (TOF and dEdx)
      const DeDxData* dedxSObj = NULL;
      const DeDxData* dedxMObj = NULL;
      if (!track.isNull()) {
        dedxSObj = &dEdxSCollH->get(track.key());
        dedxMObj = &dEdxMCollH->get(track.key());
      }

      const reco::MuonTimeExtra* tof = NULL;
      const reco::MuonTimeExtra* dttof = NULL;
      const reco::MuonTimeExtra* csctof = NULL;
      if (TypeMode > 1 && !hscp.muonRef().isNull()) {
        tof = &TOFCollH->get(hscp.muonRef().key());
        dttof = &TOFDTCollH->get(hscp.muonRef().key());
        csctof = &TOFCSCCollH->get(hscp.muonRef().key());
      }

      //Recompute dE/dx on the fly
      if (dedxSObj) {
        dedxMObj = dEdxEstimOnTheFly(ev, track, dedxMObj, dEdxSF, false, useClusterCleaning);
        dedxSObj = dEdxOnTheFly(ev, track, dedxSObj, dEdxSF, dEdxTemplates, false, useClusterCleaning);
      }

      //                       printf("Run=%6i Event=%6i --> %6.2f   %6.2f\n",ev.eventAuxiliary().run(),ev.eventAuxiliary().event(), track->pt(), dedxMObj->dEdx());

      //check if the candiate pass the preselection cuts
      if (!PassPreselection(hscp, dedxSObj, dedxMObj, tof, dttof, csctof, ev, NULL, 0))
        continue;
      PSPassed = true;

      //check if the candiate pass the selection cuts
      if (!PassSelection(hscp, dedxSObj, dedxMObj, tof, ev, CutIndex, NULL, false, 0))
        continue;

      //compute the mass of the candidate
      double Mass = -1;
      if (dedxMObj)
        Mass = GetMass(track->p(), dedxMObj->dEdx(), !isData);
      for (int M = 0; M < 6; M++) {
        if (Mass > 100.0 * M)
          PassedM[M] = true;
      }
      if (Mass < paperMassCut)
        continue;
      Passed = true;
    }  // end of Track Loop

    if (PSPassed) {
      NPSEvents += Event_Weight;
      NPSEventsErr += Event_Weight * Event_Weight;
    }
    if (Passed) {
      NSEvents += Event_Weight;
      NSEventsErr += Event_Weight * Event_Weight;
    }
    for (int M = 0; M < 6; M++) {
      if (PassedM[M]) {
        NSEventsM[M] += Event_Weight;
        NSEventsMErr[M] += Event_Weight * Event_Weight;
      }
    }
  }
  printf("\n");  // end of Event Loop

  system("mkdir -p pictures");
  string outputname = "pictures/Std_" + samples[s].Name + ".txt";
  char OutTxt[20000];
  sprintf(OutTxt, "");
  printf("%30s M>%3.0f Efficiencies: Trigger=%6.2f%%+-%6.2f%%  Offline=%6.2f%%+-%6.2f%%\n",
         MODE.c_str(),
         paperMassCut,
         100.0 * NTEvents / NEvents,
         100.0 * sqrt(pow(sqrt(NTEventsErr) / NEvents, 2) + pow(NTEvents * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
         100.0 * NSEvents / NEvents,
         100.0 * sqrt(pow(sqrt(NSEventsErr) / NEvents, 2) + pow(NSEvents * sqrt(NEventsErr) / pow(NEvents, 2), 2)));
  for (unsigned int M = 0; M < 6; M++) {
    double NPred = -1, NPredErr = -1, NData = -1;
    if (M == 0) {
      NPred = 44.0;
      NPredErr = 9.0;
      NData = 42;
    } else if (M == 1) {
      NPred = 5.6;
      NPredErr = 1.1;
      NData = 7;
    } else if (M == 2) {
      NPred = 0.56;
      NPredErr = 0.11;
      NData = 0;
    } else {
      NPred = 0.02;
      NPredErr = 0.004;
      NData = 0;
    }

    double limits[] = {-1, -1, -1, -1, -1, -1};
    printf("Compute the limits for M=%i\n", M);
    fflush(stdout);
    computeSimpleLimits(
        string("/tmp/combine_") + MODE + ".dat",
        "TkTOF",
        MODE,
        NData,
        NPred,
        1.0 + NPredErr / NPred,
        NSEventsM[M] / NEvents,
        1.0 + sqrt(pow(sqrt(NSEventsMErr[M]) / NEvents, 2) + pow(NSEventsM[M] * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
        1.32,
        IntegratedLuminosity8TeV,
        1.044,
        limits);

    sprintf(OutTxt,
            "%s%30s M>%3.0f Efficiencies: Trigger=%6.2f%%+-%6.2f%%  Presel=%6.2f%%+-%6.2f%% Offline=%6.2f%%+-%6.2f%%  "
            "Limits=%E %E %E %E %E %E\n",
            OutTxt,
            MODE.c_str(),
            M * 100.0,
            100.0 * NTEvents / NEvents,
            100.0 * sqrt(pow(sqrt(NTEventsErr) / NEvents, 2) + pow(NTEvents * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
            100.0 * NPSEvents / NEvents,
            100.0 * sqrt(pow(sqrt(NPSEventsErr) / NEvents, 2) + pow(NPSEvents * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
            100.0 * NSEventsM[M] / NEvents,
            100.0 * sqrt(pow(sqrt(NSEventsMErr[M]) / NEvents, 2) +
                         pow(NSEventsM[M] * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
            limits[0],
            limits[1],
            limits[2],
            limits[3],
            limits[4],
            limits[5]);
  }

  double limits[] = {-1, -1, -1, -1, -1, -1};
  sprintf(OutTxt,
          "%s%30s M>%3.0f Efficiencies: Trigger=%6.2f%%+-%6.2f%%  Presel=%6.2f%%+-%6.2f%% Offline=%6.2f%%+-%6.2f%%  "
          "Limits=%E %E %E %E %E %E\n",
          OutTxt,
          MODE.c_str(),
          paperMassCut,
          100.0 * NTEvents / NEvents,
          100.0 * sqrt(pow(sqrt(NTEventsErr) / NEvents, 2) + pow(NTEvents * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
          100.0 * NPSEvents / NEvents,
          100.0 * sqrt(pow(sqrt(NPSEventsErr) / NEvents, 2) + pow(NPSEvents * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
          100.0 * NSEvents / NEvents,
          100.0 * sqrt(pow(sqrt(NSEventsErr) / NEvents, 2) + pow(NSEvents * sqrt(NEventsErr) / pow(NEvents, 2), 2)),
          limits[0],
          limits[1],
          limits[2],
          limits[3],
          limits[4],
          limits[5]);

  FILE* pFile = fopen(outputname.c_str(), "w");
  fprintf(pFile, "%s", OutTxt);
  fclose(pFile);

  delete RNG;
}
