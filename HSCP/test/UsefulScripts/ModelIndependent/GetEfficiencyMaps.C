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

void GetGenHSCPBetaCustom(const std::vector<reco::GenParticle>& genColl, double& beta1, bool onlyCharged, int& index) {
  beta1 = -1;
  for (unsigned int g = 0; g < genColl.size(); g++) {
    if (genColl[g].pt() < 5)
      continue;
    if (genColl[g].status() != 1)
      continue;
    int AbsPdg = abs(genColl[g].pdgId());
    if (AbsPdg < 1000000)
      continue;
    if (onlyCharged && (AbsPdg == 1000993 || AbsPdg == 1009313 || AbsPdg == 1009113 || AbsPdg == 1009223 ||
                        AbsPdg == 1009333 || AbsPdg == 1092114 || AbsPdg == 1093214 || AbsPdg == 1093324))
      continue;  //Skip neutral gluino RHadrons
    if (onlyCharged && (AbsPdg == 1000622 || AbsPdg == 1000642 || AbsPdg == 1006113 || AbsPdg == 1006311 ||
                        AbsPdg == 1006313 || AbsPdg == 1006333))
      continue;  //skip neutral stop RHadrons
    if (beta1 < 0) {
      beta1 = genColl[g].p() / genColl[g].energy();
      index = g;
      return;
    }
  }
}

void GetEfficiencyMaps(string MODE = "COMPILE", int TypeMode_ = 0, string inputURL = "") {
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
    CutI.push_back(0.0);
    CutTOF.push_back(-1);
    CutPt.push_back(70.0);
    CutI.push_back(0.4);
    CutTOF.push_back(-1);
    CutPt.push_back(70.0);
    CutI.push_back(0.4);
    CutTOF.push_back(-1);
  } else if (TypeMode == 2) {
    CutPt.push_back(70.0);
    CutI.push_back(0.0);
    CutTOF.push_back(-1);
    CutPt.push_back(70.0);
    CutI.push_back(0.125);
    CutTOF.push_back(-1);
    CutPt.push_back(70.0);
    CutI.push_back(0.125);
    CutTOF.push_back(1.225);
  }
  int CutIndex = 0;

  //initialize LumiReWeighting
  BgLumiMC.clear();
  TrueDist.clear();
  TrueDistSyst.clear();
  for (int i = 0; i < 60; ++i)
    BgLumiMC.push_back(Pileup_MC_Summer2012[i]);

  for (int i = 0; i < 60; ++i)
    TrueDist.push_back(TrueDist2012_f[i]);
  for (int i = 0; i < 60; ++i)
    TrueDistSyst.push_back(TrueDist2012_XSecShiftUp_f[i]);
  LumiWeightsMC = edm::LumiReWeighting(BgLumiMC, TrueDist);
  LumiWeightsMCSyst = edm::LumiReWeighting(BgLumiMC, TrueDistSyst);

  //////////////////////////////////////////////
  //prepare to save all maps to histos
  TFile* OutputHisto = new TFile((string("pictures/") + "/Histos_" + sampleName + ".root").c_str(), "RECREATE");

  //       double Bins_Pt  [] = {45, 70, 100, 150, 300, 500, 1000, 9999};
  //       double Bins_Pt  [] = {45, 75, 100, 120, 140, 160, 180, 200, 225, 250, 275, 300, 325, 350, 375, 400, 450, 500, 550, 600, 650, 700, 750, 800, 850, 900, 950, 1000, 1100, 1200, 1300, 1400, 1500, 1700, 1900, 2100, 2500, 9999};
  double Bins_Pt[] = {45,   50,   60,   70,   80,   90,   100,  110,  120,  130,  140, 150,  160,  170,
                      180,  190,  200,  220,  240,  260,  280,  300,  325,  350,  375, 400,  425,  450,
                      475,  500,  550,  600,  650,  700,  750,  800,  850,  900,  950, 1000, 1100, 1200,
                      1300, 1400, 1500, 1600, 1700, 1800, 1900, 2000, 2250, 2500, 9999};
  double Bins_Beta[] = {0.00, 0.05, 0.10, 0.15, 0.20, 0.25, 0.30, 0.35, 0.40, 0.45, 0.50,
                        0.55, 0.60, 0.65, 0.70, 0.75, 0.80, 0.85, 0.90, 0.95, 1.00};
  double Bins_Eta[] = {0.0, 0.25, 0.50, 0.75, 1.00, 1.10, 1.125, 1.50, 1.75, 2.00, 2.10};

  TString Name;
  Name = "Beta_Gen";
  TH3F* Beta_Gen = new TH3F(Name,
                            Name + ";p_{T};#beta;|#eta|",
                            sizeof(Bins_Pt) / sizeof(double) - 1,
                            Bins_Pt,
                            sizeof(Bins_Beta) / sizeof(double) - 1,
                            Bins_Beta,
                            sizeof(Bins_Eta) / sizeof(double) - 1,
                            Bins_Eta);
  Beta_Gen->Sumw2();
  Name = "Beta_GenChaged";
  TH3F* Beta_GenCharged = new TH3F(Name,
                                   Name + ";p_{T};#beta;|#eta|",
                                   sizeof(Bins_Pt) / sizeof(double) - 1,
                                   Bins_Pt,
                                   sizeof(Bins_Beta) / sizeof(double) - 1,
                                   Bins_Beta,
                                   sizeof(Bins_Eta) / sizeof(double) - 1,
                                   Bins_Eta);
  Beta_GenCharged->Sumw2();
  Name = "Beta_Triggered";
  TH3F* Beta_Triggered = new TH3F(Name,
                                  Name + ";p_{T};#beta;|#eta|",
                                  sizeof(Bins_Pt) / sizeof(double) - 1,
                                  Bins_Pt,
                                  sizeof(Bins_Beta) / sizeof(double) - 1,
                                  Bins_Beta,
                                  sizeof(Bins_Eta) / sizeof(double) - 1,
                                  Bins_Eta);
  Beta_Triggered->Sumw2();
  Name = "Beta_Skimmed";
  TH3F* Beta_Skimmed = new TH3F(Name,
                                Name + ";p_{T};#beta;|#eta|",
                                sizeof(Bins_Pt) / sizeof(double) - 1,
                                Bins_Pt,
                                sizeof(Bins_Beta) / sizeof(double) - 1,
                                Bins_Beta,
                                sizeof(Bins_Eta) / sizeof(double) - 1,
                                Bins_Eta);
  Beta_Skimmed->Sumw2();
  Name = "Beta_Matched";
  TH3F* Beta_Matched = new TH3F(Name,
                                Name + ";p_{T};#beta;|#eta|",
                                sizeof(Bins_Pt) / sizeof(double) - 1,
                                Bins_Pt,
                                sizeof(Bins_Beta) / sizeof(double) - 1,
                                Bins_Beta,
                                sizeof(Bins_Eta) / sizeof(double) - 1,
                                Bins_Eta);
  Beta_Matched->Sumw2();
  Name = "Beta_Preselected";
  TH3F* Beta_Preselected = new TH3F(Name,
                                    Name + ";p_{T};#beta;|#eta|",
                                    sizeof(Bins_Pt) / sizeof(double) - 1,
                                    Bins_Pt,
                                    sizeof(Bins_Beta) / sizeof(double) - 1,
                                    Bins_Beta,
                                    sizeof(Bins_Eta) / sizeof(double) - 1,
                                    Bins_Eta);
  Beta_Preselected->Sumw2();
  Name = "Beta_SelectedP";
  TH3F* Beta_SelectedP = new TH3F(Name,
                                  Name + ";p_{T};#beta;|#eta|",
                                  sizeof(Bins_Pt) / sizeof(double) - 1,
                                  Bins_Pt,
                                  sizeof(Bins_Beta) / sizeof(double) - 1,
                                  Bins_Beta,
                                  sizeof(Bins_Eta) / sizeof(double) - 1,
                                  Bins_Eta);
  Beta_SelectedP->Sumw2();
  Name = "Beta_SelectedI";
  TH3F* Beta_SelectedI = new TH3F(Name,
                                  Name + ";p_{T};#beta;|#eta|",
                                  sizeof(Bins_Pt) / sizeof(double) - 1,
                                  Bins_Pt,
                                  sizeof(Bins_Beta) / sizeof(double) - 1,
                                  Bins_Beta,
                                  sizeof(Bins_Eta) / sizeof(double) - 1,
                                  Bins_Eta);
  Beta_SelectedI->Sumw2();
  Name = "Beta_SelectedT";
  TH3F* Beta_SelectedT = new TH3F(Name,
                                  Name + ";p_{T};#beta;|#eta|",
                                  sizeof(Bins_Pt) / sizeof(double) - 1,
                                  Bins_Pt,
                                  sizeof(Bins_Beta) / sizeof(double) - 1,
                                  Bins_Beta,
                                  sizeof(Bins_Eta) / sizeof(double) - 1,
                                  Bins_Eta);
  Beta_SelectedT->Sumw2();
  Name = "Beta_SelectedM0";
  TH3F* Beta_SelectedM0 = new TH3F(Name,
                                   Name + ";p_{T};#beta;|#eta|",
                                   sizeof(Bins_Pt) / sizeof(double) - 1,
                                   Bins_Pt,
                                   sizeof(Bins_Beta) / sizeof(double) - 1,
                                   Bins_Beta,
                                   sizeof(Bins_Eta) / sizeof(double) - 1,
                                   Bins_Eta);
  Beta_SelectedM0->Sumw2();
  Name = "Beta_SelectedM1";
  TH3F* Beta_SelectedM1 = new TH3F(Name,
                                   Name + ";p_{T};#beta;|#eta|",
                                   sizeof(Bins_Pt) / sizeof(double) - 1,
                                   Bins_Pt,
                                   sizeof(Bins_Beta) / sizeof(double) - 1,
                                   Bins_Beta,
                                   sizeof(Bins_Eta) / sizeof(double) - 1,
                                   Bins_Eta);
  Beta_SelectedM1->Sumw2();
  Name = "Beta_SelectedM2";
  TH3F* Beta_SelectedM2 = new TH3F(Name,
                                   Name + ";p_{T};#beta;|#eta|",
                                   sizeof(Bins_Pt) / sizeof(double) - 1,
                                   Bins_Pt,
                                   sizeof(Bins_Beta) / sizeof(double) - 1,
                                   Bins_Beta,
                                   sizeof(Bins_Eta) / sizeof(double) - 1,
                                   Bins_Eta);
  Beta_SelectedM2->Sumw2();
  Name = "Beta_SelectedM3";
  TH3F* Beta_SelectedM3 = new TH3F(Name,
                                   Name + ";p_{T};#beta;|#eta|",
                                   sizeof(Bins_Pt) / sizeof(double) - 1,
                                   Bins_Pt,
                                   sizeof(Bins_Beta) / sizeof(double) - 1,
                                   Bins_Beta,
                                   sizeof(Bins_Eta) / sizeof(double) - 1,
                                   Bins_Eta);
  Beta_SelectedM3->Sumw2();
  Name = "Beta_SelectedM4";
  TH3F* Beta_SelectedM4 = new TH3F(Name,
                                   Name + ";p_{T};#beta;|#eta|",
                                   sizeof(Bins_Pt) / sizeof(double) - 1,
                                   Bins_Pt,
                                   sizeof(Bins_Beta) / sizeof(double) - 1,
                                   Bins_Beta,
                                   sizeof(Bins_Eta) / sizeof(double) - 1,
                                   Bins_Eta);
  Beta_SelectedM4->Sumw2();
  Name = "Beta_SelectedM5";
  TH3F* Beta_SelectedM5 = new TH3F(Name,
                                   Name + ";p_{T};#beta;|#eta|",
                                   sizeof(Bins_Pt) / sizeof(double) - 1,
                                   Bins_Pt,
                                   sizeof(Bins_Beta) / sizeof(double) - 1,
                                   Bins_Beta,
                                   sizeof(Bins_Eta) / sizeof(double) - 1,
                                   Bins_Eta);
  Beta_SelectedM5->Sumw2();

  bool isData = false;
  bool isMC = false;
  bool isSignal = true;

  dEdxTemplates = loadDeDxTemplate("../../../data/Discrim_Templates_MC_2012.root");
  dEdxSF = 1.05;

  //do two loops through signal for samples with and without trigger changes.
  //load the files corresponding to this sample
  std::vector<string> FileName;
  FileName.push_back(inputURL);
  fwlite::ChainEvent ev(FileName);

  //loop once on the MC events to get the PU normalization
  double PUSystFactor = 1.0;  //not used
  double SampleWeight = ev.size();
  double NMCevents = 0;
  for (Long64_t ientry = 0; ientry < ev.size(); ientry++) {
    ev.to(ientry);
    NMCevents += GetPUWeight(ev, "S10", PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);
  }
  SampleWeight /= NMCevents;

  //Loop on the events
  printf("Progressing Bar                   :0%%       20%%       40%%       60%%       80%%       100%%\n");
  printf("Looping on Tree                   :");
  int TreeStep = ev.size() / 50;
  if (TreeStep == 0)
    TreeStep = 1;

  //real loop on events
  for (Long64_t ientry = 0; ientry < ev.size(); ientry++) {
    ev.to(ientry);
    if (MaxEntry > 0 && ientry > MaxEntry)
      break;
    if (ientry % TreeStep == 0) {
      printf(".");
      fflush(stdout);
    }
    Event_Weight = SampleWeight * GetPUWeight(ev, "S10", PUSystFactor, LumiWeightsMC, LumiWeightsMCSyst);

    //get the collection of generated Particles
    std::vector<reco::GenParticle> genColl;
    fwlite::Handle<std::vector<reco::GenParticle> > genCollHandle;
    genCollHandle.getByLabel(ev, "genParticles");
    if (!genCollHandle.isValid()) {
      printf("GenParticle Collection NotFound\n");
      continue;
    }
    genColl = *genCollHandle;

    double BETA;
    int genIndex = -1;
    GetGenHSCPBetaCustom(genColl, BETA, false, genIndex);
    if (BETA < 0)
      continue;
    double PT = genColl[genIndex].pt();
    double ETA = fabs(genColl[genIndex].eta());
    if (genIndex < 0 || genColl[genIndex].pt() < 45 || fabs(genColl[genIndex].eta()) > 2.1)
      continue;

    Beta_Gen->Fill(PT, BETA, ETA, Event_Weight);
    GetGenHSCPBetaCustom(genColl, BETA, true, genIndex);
    if (BETA < 0)
      continue;  //make sure there is a charged HSCP
    Beta_GenCharged->Fill(PT, BETA, ETA, Event_Weight);

    //check if the event is passing trigger
    if (!PassTriggerCustom(ev))
      continue;
    Beta_Triggered->Fill(PT, BETA, ETA, Event_Weight);

    //if(!PassSkimCustom(ev) )continue;
    Beta_Skimmed->Fill(PT, BETA, ETA, Event_Weight);

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
      Beta_Matched->Fill(PT, BETA, ETA, Event_Weight);

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

      //check if the candiate pass the preselection cuts
      if (!PassPreselection(hscp, dedxSObj, dedxMObj, tof, dttof, csctof, ev, NULL, 0))
        continue;
      Beta_Preselected->Fill(PT, BETA, ETA, Event_Weight);

      //check if the candiate pass the selection cuts
      if (genColl[ClosestGen].pt() < 45 || fabs(genColl[ClosestGen].eta()) > 2.1)
        continue;

      if (!PassSelection(hscp, dedxSObj, dedxMObj, tof, ev, CutIndex + 0, NULL, false, 0))
        continue;
      Beta_SelectedP->Fill(PT, BETA, ETA, Event_Weight);
      if (!PassSelection(hscp, dedxSObj, dedxMObj, tof, ev, CutIndex + 1, NULL, false, 0))
        continue;
      Beta_SelectedI->Fill(PT, BETA, ETA, Event_Weight);
      if (!PassSelection(hscp, dedxSObj, dedxMObj, tof, ev, CutIndex + 2, NULL, false, 0))
        continue;
      Beta_SelectedT->Fill(PT, BETA, ETA, Event_Weight);

      //compute the mass of the candidate
      double Mass = -1;
      if (dedxMObj)
        Mass = GetMass(track->p(), dedxMObj->dEdx(), !isData);
      if (Mass < 0)
        continue;
      Beta_SelectedM0->Fill(PT, BETA, ETA, Event_Weight);
      if (Mass < 100)
        continue;
      Beta_SelectedM1->Fill(PT, BETA, ETA, Event_Weight);
      if (Mass < 200)
        continue;
      Beta_SelectedM2->Fill(PT, BETA, ETA, Event_Weight);
      if (Mass < 300)
        continue;
      Beta_SelectedM3->Fill(PT, BETA, ETA, Event_Weight);
      if (Mass < 400)
        continue;
      Beta_SelectedM4->Fill(PT, BETA, ETA, Event_Weight);
      if (Mass < 500)
        continue;
      Beta_SelectedM5->Fill(PT, BETA, ETA, Event_Weight);

    }  // end of Track Loop
  }
  printf("\n");  // end of Event Loop
  OutputHisto->Write();
  OutputHisto->Close();
}
