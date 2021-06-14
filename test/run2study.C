#define run2study_cxx
#include "run2study.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

void run2study::Loop(int year, TString Letter)
{
//   In a ROOT session, you can do:
//      root> .L run2study.C
//      root> run2study t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
   if (fChain == 0) return;

//   Long64_t nentries = fChain->GetEntriesFast();
   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;

   double pi=acos(-1);
   double dEdxSF [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF

// histograms

    
    TH1D* GenHSCP_probQ_wCuts = new TH1D("GenHSCP_probQ_wCuts", "Combined ProbQ on tracks for gen HSCP w/ cuts;Combined on-track charge probability;Entries (1/bin)",100,0.,1.);
    TH1D* GenHSCP_probXY_wCuts = new TH1D("GenHSCP_probXY_wCuts", "Combined ProbXY on tracks for gen HSCP w/ cuts;Combined on-track shape probability;Entries (1/bin)",100,0.,1.);
    TH1D* GenHSCP_probQ_wNoCuts = new TH1D("GenHSCP_probQ_wNoCuts", "Combined ProbQ on tracks  for gen HSCP w/ no cuts;Combined on-track charge probability;Entries (1/bin)",100,0.,1.);
    TH1D* GenHSCP_probXY_wNoCuts = new TH1D("GenHSCP_probXY_wNoCuts", "Combined ProbXY on tracks for gen HSCP w/ no cuts;Combined on-track shape probability;Entries (1/bin)",100,0.,1.);
    TH1D* GenBckg_probQ_wCuts = new TH1D("GenBckg_probQ_wCuts", "Combined ProbQ on tracks for gen background w/ cuts;Combined on-track charge probability;Entries (1/bin)",100,0.,1.);
    TH1D* GenBckg_probXY_wCuts = new TH1D("GenBckg_probXY_wCuts", "Combined ProbXY on tracks for gen background w/ cuts;Combined on-track shape probability;Entries (1/bin)",100,0.,1.);
    TH1D* GenBckg_probQ_wNoCuts = new TH1D("GenBckg_probQ_wNoCuts", "Combined ProbQ on tracks for gen background w/ no cuts;Combined on-track charge probability;Entries (1/bin)",100,0.,1.);
    TH1D* GenBckg_probXY_wNoCuts = new TH1D("GenBckg_probXY_wNoCuts", "Combined ProbXY on tracks  for gen background w/ no cuts;Combined on-track shape probability;Entries (1/bin)",100,0.,1.);
    
    
   TH1D* HNtracks = new TH1D("HNtracks", "HNtracks", 50, -0.5,49.5);
   TH1D* HNtracks1 = new TH1D("HNtracks1", "Ntracks with pT>1", 40, -0.5,39.5);
   TH1D* HNtracks20 = new TH1D("HNtracks20", "Ntracks with pT>20", 30, -0.5,29.5);
   TH1D* HNtracks50 = new TH1D("HNtracks50", "Ntracks with pT>50", 20, -0.5,19.5);
   TH1D* HNtracks100 = new TH1D("HNtracks100", "Ntracks with pT>100", 15, -0.5,14.5);

   TH1D* Htrackpt = new TH1D("Htrackpt", "track pT", 500, 0.,500);
   TH1D* Htracketa = new TH1D("Htracketa", "track eta", 24, -3.,3.);
   TH1D* Htracketa_lowp = new TH1D("Htracketa_lowp", "track eta", 24, -3.,3.);
   TH1D* Htrackphi = new TH1D("Htrackphi", "track phi", 24, -1.*pi,pi);
   TH1D* Htracknhit = new TH1D("Htracknhit", "track nhit", 50, -0.5,49.5);

   TH1D* Htrackih = new TH1D("Htrackih", "track ih", 50, 0.,20);
   TH1D* Htrackih_reco = new TH1D("Htrackih_reco", "track ih (reco)", 50, 0.,20);
   TH1D* Htrackih_pix = new TH1D("Htrackih_pix", "track ih in pix", 50, 0.,20);
   TH1D* Htrackih_strip = new TH1D("Htrackih_strip", "track ih in strip", 50, 0.,20);
   TH1D* Htrackdedx_pix = new TH1D("Htrackdedx_pix", "track dedx in pix", 100, 0.,500000);
   TH1D* Htrackdedx_strip = new TH1D("Htrackdedx_strip", "track dedx in strip", 100, 0.,2000);

   TH1D* Htrackih_lowp = new TH1D("Htrackih_lowp", "track ih", 50, 0.,20);
   TH1D* Htrackih_reco_lowp = new TH1D("Htrackih_reco_lowp", "track ih (reco)", 50, 0.,20);
   TH1D* Htrackih_pix_lowp = new TH1D("Htrackih_pix_lowp", "track ih in pix", 50, 0.,20);
   TH1D* Htrackih_strip_lowp = new TH1D("Htrackih_strip_lowp", "track ih in strip", 50, 0.,20);
   TH1D* Htrackdedx_pix_lowp = new TH1D("Htrackdedx_pix_lowp", "track dedx in pix", 100, 0.,500000);
   TH1D* Htrackdedx_strip_lowp = new TH1D("Htrackdedx_strip_lowp", "track dedx in strip", 100, 0.,2000);
   TH1D* Htrackdedx_strip_lowp1 = new TH1D("Htrackdedx_strip_lowp1", "track dedx in strip", 100, 0.,20);
   TH1D* Htrackdedx_strip_lowp2 = new TH1D("Htrackdedx_strip_lowp2", "track dedx in strip", 100, 0.,20);

   TH1D* Nsat = new TH1D("Nsat", "Nsat", 20, -0.5,19.5);
   TH1D* NPix = new TH1D("NPix", "NPix", 10, -0.5,9.5);
   TH1D* NStrip = new TH1D("NStrip", "NStrip", 30, -0.5,29.5);
   TH2D* dEdXVsP = new TH2D("dEdXVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXVsP_lowp = new TH2D("dEdXVsP_lowp", "dEdX:P", 50,0,5, 50, 2.,7.);
   TH2D* dEdXVsP_lowp2 = new TH2D("dEdXVsP_lowp2", "dEdX:P", 125,0,25, 50, 2.,7.);
   TH2D* dEdXVsP_lowp3 = new TH2D("dEdXVsP_lowp3", "dEdX:P", 100,0,20, 50, 2.,7.);
   TH2D* dEdXpixVsP = new TH2D("dEdXpixVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXpixVsP_lowp = new TH2D("dEdXpixVsP_lowp", "dEdX:P", 50,0,5, 50, 0.,7.);
   TH2D* dEdXstripVsP = new TH2D("dEdXstripVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXstripVsP_lowp = new TH2D("dEdXstripVsP_lowp", "dEdX:P", 50,0,5, 50, 2.,7.);
   TH2D* dEstrVsdE_lowp = new TH2D("dEstrVsdE_lowp", "dEdX(all):dEdX(strip)", 50,2.,7., 50, 2.,7.);
   TH2D* dEdXstripVsEta_lowp = new TH2D("dEdXstripVsEta_lowp", "dEdX:P", 12, -3.,3., 50, 2.,7.);
   TH2D* EtaVsPhi_nhit8 = new TH2D("EtaVsPhi_nhit8", "Eta:Phi", 12, -1.*pi,pi,12, -3.,3.);
   TH2D* EtaVsPhi_nhit = new TH2D("EtaVsPhi_nhit", "Eta:Phi", 12, -1.*pi,pi,12, -3.,3.);

   TH2D* NPixVsNStrip_lowIhp = new TH2D("NPixVsNStrip_lowIhp", "NPix:NStrip", 20, -0.5,19.5, 10, -0.5,9.5);
   TH2D* NPixVsNStrip_highIhp = new TH2D("NPixVsNStrip_highIhp", "NPix:NStrip", 20, -0.5,19.5, 10, -0.5,9.5);

   TH2D* dEdXstripVsNhit_lowp  = new TH2D("dEdXstripVsNhit_lowp", "dEdX(strip):Nhit", 20, -0.5,19.5,50, 2.,7.);
   TH2D* dEdXstripVsNhittrunc_lowp  = new TH2D("dEdXstripVsNhittrunc_lowp", "dEdX(strip):Nhit", 20, -0.5,19.5,50, 2.,7.);
   TH2D* dEdXstripVsCharge_lowp = new TH2D("dEdXstripVsCharge_lowp", "dEdX(strip):Charge", 100, 0.,20, 50, 2.,7.);

   TH1D* Charge_pixl1 = new TH1D("Charge_pixl1","Charge per layer",100, 0.,20);
   TH1D* Charge_pixl2 = new TH1D("Charge_pixl2","Charge per layer",100, 0.,20);
   TH1D* Charge_pixl3 = new TH1D("Charge_pixl3","Charge per layer",100, 0.,20);
   TH1D* Charge_pixl4 = new TH1D("Charge_pixl4","Charge per layer",100, 0.,20);
   TH1D* Charge_pixd1 = new TH1D("Charge_pixd1","Charge per layer",100, 0.,20);
   TH1D* Charge_pixd2 = new TH1D("Charge_pixd2","Charge per layer",100, 0.,20);
   TH1D* Charge_pixd3 = new TH1D("Charge_pixd3","Charge per layer",100, 0.,20);
   TH1D* Charge_tib1 = new TH1D("Charge_tib1","Charge per layer",100, 0.,20);
   TH1D* Charge_tib2 = new TH1D("Charge_tib2","Charge per layer",100, 0.,20);
   TH1D* Charge_tib3 = new TH1D("Charge_tib3","Charge per layer",100, 0.,20);
   TH1D* Charge_tib4 = new TH1D("Charge_tib4","Charge per layer",100, 0.,20);
   TH1D* Charge_tob1 = new TH1D("Charge_tob1","Charge per layer",100, 0.,20);
   TH1D* Charge_tob2 = new TH1D("Charge_tob2","Charge per layer",100, 0.,20);
   TH1D* Charge_tob3 = new TH1D("Charge_tob3","Charge per layer",100, 0.,20);
   TH1D* Charge_tob4 = new TH1D("Charge_tob4","Charge per layer",100, 0.,20);
   TH1D* Charge_tob5 = new TH1D("Charge_tob5","Charge per layer",100, 0.,20);
   TH1D* Charge_tob6 = new TH1D("Charge_tob6","Charge per layer",100, 0.,20);
   TH1D* Charge_tid1 = new TH1D("Charge_tid1","Charge per layer",100, 0.,20);
   TH1D* Charge_tid2 = new TH1D("Charge_tid2","Charge per layer",100, 0.,20);
   TH1D* Charge_tid3 = new TH1D("Charge_tid3","Charge per layer",100, 0.,20);
   TH1D* Charge_tec1 = new TH1D("Charge_tec1","Charge per layer",100, 0.,20);
   TH1D* Charge_tec2 = new TH1D("Charge_tec2","Charge per layer",100, 0.,20);
   TH1D* Charge_tec3 = new TH1D("Charge_tec3","Charge per layer",100, 0.,20);
   TH1D* Charge_tec4 = new TH1D("Charge_tec4","Charge per layer",100, 0.,20);
   TH1D* Charge_tec5 = new TH1D("Charge_tec5","Charge per layer",100, 0.,20);
   TH1D* Charge_tec6 = new TH1D("Charge_tec6","Charge per layer",100, 0.,20);
   TH1D* Charge_tec7 = new TH1D("Charge_tec7","Charge per layer",100, 0.,20);
   TH1D* Charge_tec8 = new TH1D("Charge_tec8","Charge per layer",100, 0.,20);
   TH1D* Charge_tec9 = new TH1D("Charge_tec9","Charge per layer",100, 0.,20);

   TH1D* LowCharge_tib1 = new TH1D("LowCharge_tib1","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tib2 = new TH1D("LowCharge_tib2","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tib3 = new TH1D("LowCharge_tib3","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tib4 = new TH1D("LowCharge_tib4","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tob1 = new TH1D("LowCharge_tob1","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tob2 = new TH1D("LowCharge_tob2","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tob3 = new TH1D("LowCharge_tob3","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tob4 = new TH1D("LowCharge_tob4","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tob5 = new TH1D("LowCharge_tob5","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tob6 = new TH1D("LowCharge_tob6","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tid1 = new TH1D("LowCharge_tid1","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tid2 = new TH1D("LowCharge_tid2","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tid3 = new TH1D("LowCharge_tid3","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec1 = new TH1D("LowCharge_tec1","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec2 = new TH1D("LowCharge_tec2","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec3 = new TH1D("LowCharge_tec3","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec4 = new TH1D("LowCharge_tec4","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec5 = new TH1D("LowCharge_tec5","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec6 = new TH1D("LowCharge_tec6","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec7 = new TH1D("LowCharge_tec7","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec8 = new TH1D("LowCharge_tec8","Charge per layer",100, 0.,20);
   TH1D* LowCharge_tec9 = new TH1D("LowCharge_tec9","Charge per layer",100, 0.,20);

   TH2D* ChargeVsRun_pixl1 = new TH2D("ChargeVsRun_pixl1","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_pixl2 = new TH2D("ChargeVsRun_pixl2","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_pixl3 = new TH2D("ChargeVsRun_pixl3","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_pixl4 = new TH2D("ChargeVsRun_pixl4","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_pixd1 = new TH2D("ChargeVsRun_pixd1","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_pixd2 = new TH2D("ChargeVsRun_pixd2","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_pixd3 = new TH2D("ChargeVsRun_pixd3","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tib1 = new TH2D("ChargeVsRun_tib1","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tib2 = new TH2D("ChargeVsRun_tib2","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tib3 = new TH2D("ChargeVsRun_tib3","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tib4 = new TH2D("ChargeVsRun_tib4","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tob1 = new TH2D("ChargeVsRun_tob1","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tob2 = new TH2D("ChargeVsRun_tob2","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tob3 = new TH2D("ChargeVsRun_tob3","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tob4 = new TH2D("ChargeVsRun_tob4","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tob5 = new TH2D("ChargeVsRun_tob5","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tob6 = new TH2D("ChargeVsRun_tob6","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tid1 = new TH2D("ChargeVsRun_tid1","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tid2 = new TH2D("ChargeVsRun_tid2","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tid3 = new TH2D("ChargeVsRun_tid3","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec1 = new TH2D("ChargeVsRun_tec1","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec2 = new TH2D("ChargeVsRun_tec2","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec3 = new TH2D("ChargeVsRun_tec3","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec4 = new TH2D("ChargeVsRun_tec4","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec5 = new TH2D("ChargeVsRun_tec5","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec6 = new TH2D("ChargeVsRun_tec6","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec7 = new TH2D("ChargeVsRun_tec7","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec8 = new TH2D("ChargeVsRun_tec8","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);
   TH2D* ChargeVsRun_tec9 = new TH2D("ChargeVsRun_tec9","ChargeVsRun per layer",545, 271000,325500,100, 0.,20);

   TH2D* ZooChargeVsRun_pixl1 = new TH2D("ZooChargeVsRun_pixl1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl2 = new TH2D("ZooChargeVsRun_pixl2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl3 = new TH2D("ZooChargeVsRun_pixl3","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl4 = new TH2D("ZooChargeVsRun_pixl4","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd1 = new TH2D("ZooChargeVsRun_pixd1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd2 = new TH2D("ZooChargeVsRun_pixd2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd3 = new TH2D("ZooChargeVsRun_pixd3","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_tib1 = new TH2D("ZooChargeVsRun_tib1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tib2 = new TH2D("ZooChargeVsRun_tib2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tib3 = new TH2D("ZooChargeVsRun_tib3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tib4 = new TH2D("ZooChargeVsRun_tib4","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob1 = new TH2D("ZooChargeVsRun_tob1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob2 = new TH2D("ZooChargeVsRun_tob2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob3 = new TH2D("ZooChargeVsRun_tob3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob4 = new TH2D("ZooChargeVsRun_tob4","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob5 = new TH2D("ZooChargeVsRun_tob5","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tob6 = new TH2D("ZooChargeVsRun_tob6","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tid1 = new TH2D("ZooChargeVsRun_tid1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tid2 = new TH2D("ZooChargeVsRun_tid2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tid3 = new TH2D("ZooChargeVsRun_tid3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec1 = new TH2D("ZooChargeVsRun_tec1","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec2 = new TH2D("ZooChargeVsRun_tec2","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec3 = new TH2D("ZooChargeVsRun_tec3","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec4 = new TH2D("ZooChargeVsRun_tec4","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec5 = new TH2D("ZooChargeVsRun_tec5","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec6 = new TH2D("ZooChargeVsRun_tec6","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec7 = new TH2D("ZooChargeVsRun_tec7","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec8 = new TH2D("ZooChargeVsRun_tec8","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);
   TH2D* ZooChargeVsRun_tec9 = new TH2D("ZooChargeVsRun_tec9","ChargeVsRun per layer",545, 271000,325500,50, 2., 4.5);

   TH1D* HSCP_dEdX = new TH1D("HSCP_dEdX", "dEdX", 60, 0.,15.);
   TH1D* HSCP_dEdXpix = new TH1D("HSCP_dEdXpix", "dEdX(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdXstrip = new TH1D("HSCP_dEdXstrip", "dEdX(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdX0 = new TH1D("HSCP_dEdX0", "dEdX0", 60, 0.,15.);
   TH1D* HSCP_dEdX0pix = new TH1D("HSCP_dEdX0pix", "dEdX0(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdX0strip = new TH1D("HSCP_dEdX0strip", "dEdX0(strip)", 60, 0.,15.);

   TH1D* HSCP_FMIP4 = new TH1D("HSCP_FMIP4", "FMIP(4)", 50, 0.,1.);
   TH1D* HSCP_FMIP3p5 = new TH1D("HSCP_FMIP3p5", "FMIP(3.5)", 50, 0.,1.);
   TH1D* HSCP_FMIP3p2 = new TH1D("HSCP_FMIP3p2", "FMIP(3.2)", 50, 0.,1.);

   TH1D* HSCP_pt = new TH1D("HSCP_pt", "pT",  50, 55.,1550);
   TH1D* HSCP_iso_eop = new TH1D("HSCP_iso_eop", "Isolation (ECAL+HCAL)/p",  50, 0., 5.);
   TH1D* nPV = new TH1D("nPV", "nPV",  40,0,80);
   TH1D* HSCP_invB = new TH1D("HSCP_invB", "invBeta",  300,-1,2);
   TH1D* HSCP_errinvB = new TH1D("HSCP_errinvB", "err_invBeta",  50,0,0.5);
   TH1D* HSCP_invBDT = new TH1D("HSCP_invBDT", "invBeta(DT)",  90,-1,2);
   TH1D* HSCP_invBCSC = new TH1D("HSCP_invBCSC", "invBeta(CSC)",  90,-1,2);
   TH1D* HSCP_time = new TH1D("HSCP_time", "VertexTiming",  200,-100,100);
   TH1D* HSCP_probQ_stdAna = new TH1D("HSCP_probQ_stdAna", "Combined ProbQ on tracks in std analysis;Combined on-track charge probability;Entries (1/bin)",100,0.,1.);
   TH1D* HSCP_probXY_stdAna = new TH1D("HSCP_probXY_stdAna", "Combined ProbXY on tracks in std analysis;Combined on-track shape probability;Entries (1/bin)",100,0.,1.);
   TH1D* HSCPCandidateMatchedGenID_probQ = new TH1D("HSCPCandidateMatchedGenID_probQ", "Combined ProbQ on tracks for matched gen HSCPs;Combined on-track charge probability;Entries (1/bin)",100,0.,1.);
   TH1D* HSCPCandidateMatchedGenID_probXY = new TH1D("HSCPCandidateMatchedGenID_probXY", "Combined ProbXY on tracks for matched gen HSCPs;Combined on-track shape probability;Entries (1/bin)",100,0.,1.);


   TH2D* dEdXVsRun = new TH2D("dEdXVsRun", "dEdX:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXpixVsRun = new TH2D("dEdXpixVsRun", "dEdX(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXstripVsRun = new TH2D("dEdXstripVsRun", "dEdX(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0VsRun = new TH2D("dEdX0VsRun", "dEdX0:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0pixVsRun = new TH2D("dEdX0pixVsRun", "dEdX0(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0stripVsRun = new TH2D("dEdX0stripVsRun", "dEdX0(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4VsRun = new TH2D("dEdXV4sRun", "dEdX:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4pixVsRun = new TH2D("dEdX4pixVsRun", "dEdX(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4stripVsRun = new TH2D("dEdX4stripVsRun", "dEdX(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40VsRun = new TH2D("dEdX40VsRun", "dEdX0:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40pixVsRun = new TH2D("dEdX40pixVsRun", "dEdX0(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40stripVsRun = new TH2D("dEdX40stripVsRun", "dEdX0(strip):Run", 545, 271000,325500, 60, 0.,15.);

   TH2D* FMIP4VsRun = new TH2D("FMIP4VsRun", "FMIP(4):Run", 545, 271000,325500, 50, 0.,1.);
   TH2D* FMIP3p5VsRun = new TH2D("FMIP3p5VsRun", "FMIP(3.5):Run", 545, 271000,325500, 50, 0.,1.);
   TH2D* FMIP3p2VsRun = new TH2D("FMIP3p2VsRun", "FMIP(3.2):Run", 545, 271000,325500, 50, 0.,1.);
   TH2D* FMIP4VsEta = new TH2D("FMIP4VsEta", "FMIP(4):P", 12, -3.,3., 50, 0.,1.);

   TH2D* NmeasVsRun = new TH2D("NmeasVsRun", "#meas:Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasPixVsRun = new TH2D("NmeasPixVsRun", "#meas(pix):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasStrVsRun = new TH2D("NmeasStrVsRun", "#meas(strip):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* Nmeas0VsRun = new TH2D("Nmeas0VsRun", "#meas:Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasPix0VsRun = new TH2D("NmeasPix0VsRun", "#meas(pix):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NmeasStr0VsRun = new TH2D("NmeasStr0VsRun", "#meas(strip):Run", 545, 271000,325500, 50, -0.5,49.5);
   TH2D* NsatVsRun = new TH2D("NsatVsRun", "#sat:Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatPixVsRun = new TH2D("NsatPixVsRun", "#sat(pix):Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatStrVsRun = new TH2D("NsatStrVsRun", "#sat(strip):Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* Nsat0VsRun = new TH2D("Nsat0VsRun", "#sat:Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatPix0VsRun = new TH2D("NsatPix0VsRun", "#sat(pix):Run", 545, 271000,325500, 10, -0.5,9.5);
   TH2D* NsatStr0VsRun = new TH2D("NsatStr0VsRun", "#sat(strip):Run", 545, 271000,325500, 10, -0.5,9.5);
/*
   TH2D* dEdXVsIL = new TH2D("dEdXVsIL", "dEdX:InstLumi", 28, 0,14000, 60, 0.,15.);
   TH2D* dEdXpixVsIL = new TH2D("dEdXpixVsIL", "dEdX(strip):InstLumi", 28, 0,14000, 60, 0.,15.);
   TH2D* dEdXstripVsIL = new TH2D("dEdXstripVsIL", "dEdX(strip):InstLumi", 28, 0,14000, 60, 0.,15.);
*/

   TH2D* ptVsRun = new TH2D("ptVsRun", "pT:Run", 545, 271000,325500, 50, 55.,1550);
   TH2D* nPVVsRun = new TH2D("nPVVsRun", "nPV:Run", 545, 271000,325500, 40,0,80);
   TH2D* invBVsRun = new TH2D("invBVsRun", "invBeta:Run", 545, 271000,325500, 90,-1,2);
   TH2D* errinvBVsRun = new TH2D("errinvBVsRun", "err_invBeta:Run", 545, 271000,325500, 50,0,0.5);
   TH2D* invBDTVsRun = new TH2D("invBDTVsRun", "invBeta(DT):Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBCSCVsRun = new TH2D("invBCSCVsRun", "invBeta(CSC):Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBnewVsRun = new TH2D("invnewBVsRun", "invBeta:Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBnewDTVsRun = new TH2D("invnewBDTVsRun", "invBeta(DT):Run", 545, 271000,325500, 90,-1,2);
   TH2D* invBnewCSCVsRun = new TH2D("invnewBCSCVsRun", "invBeta(CSC):Run", 545, 271000,325500, 90,-1,2);
   TH2D* timeVsRun = new TH2D("timeVsRun", "VertexTimung:Run", 545, 271000,325500, 100,-100,100);
   TH2D* lumiVsRun = new TH2D("lumiVsRun", "Lumi:Run", 545, 271000,325500, 56, 0,14000);


   TH2D* R1_StdEdXVsEvent = new TH2D("R1_StdEdXVsEvent", "dEdX(strip):Event",500,0.,2000000000,60, 0.,15.); 
   TH2D* R1_StdEdXVsLumi = new TH2D("R1_StdEdXVsLumi", "dEdX(strip):Lumi",80,0.,16000,60, 0.,15.); 
   TH2D* R1_LumiVsEvent = new TH2D("R1_LumiVsEvent", "Lumi:Event",500,0.,2000000000,80,0.,16000); 
   TH2D* R1_nPVVsEvent = new TH2D("R1_nPVVsEvent", "nPV:Event",500,0.,2000000000,40, 0.,80.); 
   TH1D* R1_CandVsEvent = new TH1D("R1_CandVsEvent", "Cand per 4M events vs Event", 500,0.,2000000000);

   TH2D* R2_StdEdXVsEvent = new TH2D("R2_StdEdXVsEvent", "dEdX(strip):Event",500,0.,2000000000,60, 0.,15.); 
   TH2D* R2_StdEdXVsLumi = new TH2D("R2_StdEdXVsLumi", "dEdX(strip):Lumi",80,0.,16000,60, 0.,15.); 
   TH2D* R2_LumiVsEvent = new TH2D("R2_LumiVsEvent", "Lumi:Event",500,0.,2000000000,80,0.,16000); 
   TH2D* R2_nPVVsEvent = new TH2D("R2_nPVVsEvent", "nPV:Event",500,0.,2000000000,40, 0.,80.); 
   TH1D* R2_CandVsEvent = new TH1D("R2_CandVsEvent", "Cand per 4M events vs Event", 500,0.,2000000000);

   TH2D* R3_StdEdXVsEvent = new TH2D("R3_StdEdXVsEvent", "dEdX(strip):Event",500,0.,2000000000,60, 0.,15.); 
   TH2D* R3_StdEdXVsLumi = new TH2D("R3_StdEdXVsLumi", "dEdX(strip):Lumi",80,0.,16000,60, 0.,15.); 
   TH2D* R3_LumiVsEvent = new TH2D("R3_LumiVsEvent", "Lumi:Event",500,0.,2000000000,80,0.,16000); 
   TH2D* R3_nPVVsEvent = new TH2D("R3_nPVVsEvent", "nPV:Event",500,0.,2000000000,40, 0.,80.); 
   TH1D* R3_CandVsEvent = new TH1D("R3_CandVsEvent", "Cand per 4M events vs Event", 500,0.,2000000000);
    
    
   GenHSCP_probQ_wCuts->Sumw2();
   GenHSCP_probXY_wCuts->Sumw2();
   GenHSCP_probQ_wNoCuts->Sumw2();
   GenHSCP_probXY_wNoCuts->Sumw2();
   GenBckg_probQ_wCuts->Sumw2();
   GenBckg_probXY_wCuts->Sumw2();
   GenBckg_probQ_wNoCuts->Sumw2();
   GenBckg_probXY_wNoCuts->Sumw2();

   HNtracks->Sumw2();
   HNtracks1->Sumw2();
   HNtracks20->Sumw2();
   HNtracks50->Sumw2();
   HNtracks100->Sumw2();

   Htrackpt->Sumw2();
   Htracketa->Sumw2();
   Htracketa_lowp->Sumw2();
   Htrackphi->Sumw2();
   Htracknhit->Sumw2();

   Htrackih->Sumw2();
   Htrackih_reco->Sumw2();
   Htrackih_pix->Sumw2();
   Htrackih_strip->Sumw2();
   Htrackdedx_pix->Sumw2();
   Htrackdedx_strip->Sumw2();

   Htrackih_lowp->Sumw2();
   Htrackih_reco_lowp->Sumw2();
   Htrackih_pix_lowp->Sumw2();
   Htrackih_strip_lowp->Sumw2();
   Htrackdedx_pix_lowp->Sumw2();
   Htrackdedx_strip_lowp->Sumw2();
   Htrackdedx_strip_lowp1->Sumw2();
   Htrackdedx_strip_lowp2->Sumw2();


   Nsat->Sumw2();
   NPix->Sumw2();
   NStrip->Sumw2();
   dEdXVsP->Sumw2();
   dEdXVsP_lowp->Sumw2();
   dEdXVsP_lowp2->Sumw2();
   dEdXVsP_lowp3->Sumw2();
   dEdXpixVsP->Sumw2();
   dEdXstripVsP->Sumw2();
   dEdXpixVsP_lowp->Sumw2();
   dEdXstripVsP_lowp->Sumw2();
   dEdXstripVsEta_lowp->Sumw2();
   dEstrVsdE_lowp->Sumw2();
   NPixVsNStrip_highIhp->Sumw2();
   NPixVsNStrip_lowIhp->Sumw2();
   dEdXstripVsNhit_lowp->Sumw2();
   dEdXstripVsNhittrunc_lowp->Sumw2();
   dEdXstripVsCharge_lowp->Sumw2();
   EtaVsPhi_nhit8->Sumw2();
   EtaVsPhi_nhit->Sumw2();

   Charge_pixl1->Sumw2();
   Charge_pixl2->Sumw2();
   Charge_pixl3->Sumw2();
   Charge_pixl4->Sumw2();
   Charge_pixd1->Sumw2();
   Charge_pixd2->Sumw2();
   Charge_pixd3->Sumw2();
   Charge_tib1->Sumw2();
   Charge_tib2->Sumw2();
   Charge_tib3->Sumw2();
   Charge_tib4->Sumw2();
   Charge_tob1->Sumw2();
   Charge_tob2->Sumw2();
   Charge_tob3->Sumw2();
   Charge_tob4->Sumw2();
   Charge_tob5->Sumw2();
   Charge_tob6->Sumw2();
   Charge_tid1->Sumw2();
   Charge_tid2->Sumw2();
   Charge_tid3->Sumw2();
   Charge_tec1->Sumw2();
   Charge_tec2->Sumw2();
   Charge_tec3->Sumw2();
   Charge_tec4->Sumw2();
   Charge_tec5->Sumw2();
   Charge_tec6->Sumw2();
   Charge_tec7->Sumw2();
   Charge_tec8->Sumw2();
   Charge_tec9->Sumw2();

   LowCharge_tib1->Sumw2();
   LowCharge_tib2->Sumw2();
   LowCharge_tib3->Sumw2();
   LowCharge_tib4->Sumw2();
   LowCharge_tob1->Sumw2();
   LowCharge_tob2->Sumw2();
   LowCharge_tob3->Sumw2();
   LowCharge_tob4->Sumw2();
   LowCharge_tob5->Sumw2();
   LowCharge_tob6->Sumw2();
   LowCharge_tid1->Sumw2();
   LowCharge_tid2->Sumw2();
   LowCharge_tid3->Sumw2();
   LowCharge_tec1->Sumw2();
   LowCharge_tec2->Sumw2();
   LowCharge_tec3->Sumw2();
   LowCharge_tec4->Sumw2();
   LowCharge_tec5->Sumw2();
   LowCharge_tec6->Sumw2();
   LowCharge_tec7->Sumw2();
   LowCharge_tec8->Sumw2();
   LowCharge_tec9->Sumw2();

   ChargeVsRun_pixl1->Sumw2();
   ChargeVsRun_pixl2->Sumw2();
   ChargeVsRun_pixl3->Sumw2();
   ChargeVsRun_pixl4->Sumw2();
   ChargeVsRun_pixd1->Sumw2();
   ChargeVsRun_pixd2->Sumw2();
   ChargeVsRun_pixd3->Sumw2();
   ChargeVsRun_tib1->Sumw2();
   ChargeVsRun_tib2->Sumw2();
   ChargeVsRun_tib3->Sumw2();
   ChargeVsRun_tib4->Sumw2();
   ChargeVsRun_tob1->Sumw2();
   ChargeVsRun_tob2->Sumw2();
   ChargeVsRun_tob3->Sumw2();
   ChargeVsRun_tob4->Sumw2();
   ChargeVsRun_tob5->Sumw2();
   ChargeVsRun_tob6->Sumw2();
   ChargeVsRun_tid1->Sumw2();
   ChargeVsRun_tid2->Sumw2();
   ChargeVsRun_tid3->Sumw2();
   ChargeVsRun_tec1->Sumw2();
   ChargeVsRun_tec2->Sumw2();
   ChargeVsRun_tec3->Sumw2();
   ChargeVsRun_tec4->Sumw2();
   ChargeVsRun_tec5->Sumw2();
   ChargeVsRun_tec6->Sumw2();
   ChargeVsRun_tec7->Sumw2();
   ChargeVsRun_tec8->Sumw2();
   ChargeVsRun_tec9->Sumw2();

   ZooChargeVsRun_pixl1->Sumw2();
   ZooChargeVsRun_pixl2->Sumw2();
   ZooChargeVsRun_pixl3->Sumw2();
   ZooChargeVsRun_pixl4->Sumw2();
   ZooChargeVsRun_pixd1->Sumw2();
   ZooChargeVsRun_pixd2->Sumw2();
   ZooChargeVsRun_pixd3->Sumw2();
   ZooChargeVsRun_tib1->Sumw2();
   ZooChargeVsRun_tib2->Sumw2();
   ZooChargeVsRun_tib3->Sumw2();
   ZooChargeVsRun_tib4->Sumw2();
   ZooChargeVsRun_tob1->Sumw2();
   ZooChargeVsRun_tob2->Sumw2();
   ZooChargeVsRun_tob3->Sumw2();
   ZooChargeVsRun_tob4->Sumw2();
   ZooChargeVsRun_tob5->Sumw2();
   ZooChargeVsRun_tob6->Sumw2();
   ZooChargeVsRun_tid1->Sumw2();
   ZooChargeVsRun_tid2->Sumw2();
   ZooChargeVsRun_tid3->Sumw2();
   ZooChargeVsRun_tec1->Sumw2();
   ZooChargeVsRun_tec2->Sumw2();
   ZooChargeVsRun_tec3->Sumw2();
   ZooChargeVsRun_tec4->Sumw2();
   ZooChargeVsRun_tec5->Sumw2();
   ZooChargeVsRun_tec6->Sumw2();
   ZooChargeVsRun_tec7->Sumw2();
   ZooChargeVsRun_tec8->Sumw2();
   ZooChargeVsRun_tec9->Sumw2();


   dEdXVsRun->Sumw2();
   dEdXpixVsRun->Sumw2();
   dEdXstripVsRun->Sumw2();
   dEdX0VsRun->Sumw2();
   dEdX0pixVsRun->Sumw2();
   dEdX0stripVsRun->Sumw2();
   dEdX4VsRun->Sumw2();
   dEdX4pixVsRun->Sumw2();
   dEdX4stripVsRun->Sumw2();
   dEdX40VsRun->Sumw2();
   dEdX40pixVsRun->Sumw2();
   dEdX40stripVsRun->Sumw2();

   FMIP4VsRun->Sumw2();
   FMIP3p5VsRun->Sumw2();
   FMIP3p2VsRun->Sumw2();
   FMIP4VsEta->Sumw2();

   NmeasVsRun->Sumw2();
   NmeasPixVsRun->Sumw2();
   NmeasStrVsRun->Sumw2();
   Nmeas0VsRun->Sumw2();
   NmeasPix0VsRun->Sumw2();
   NmeasStr0VsRun->Sumw2();
   NsatVsRun->Sumw2();
   NsatPixVsRun->Sumw2();
   NsatStrVsRun->Sumw2();
   Nsat0VsRun->Sumw2();
   NsatPix0VsRun->Sumw2();
   NsatStr0VsRun->Sumw2();

   ptVsRun->Sumw2();
   nPVVsRun->Sumw2();
   invBVsRun->Sumw2();
   timeVsRun->Sumw2();
   lumiVsRun->Sumw2();
   errinvBVsRun->Sumw2();
   invBDTVsRun->Sumw2();
   invBCSCVsRun->Sumw2();
   invBnewVsRun->Sumw2();
   invBnewDTVsRun->Sumw2();
   invBnewCSCVsRun->Sumw2();

   HSCP_dEdX->Sumw2(); 
   HSCP_dEdXpix->Sumw2(); 
   HSCP_dEdXstrip->Sumw2(); 

   HSCP_FMIP4->Sumw2(); 
   HSCP_FMIP3p5->Sumw2(); 
   HSCP_FMIP3p2->Sumw2(); 
   HSCP_pt->Sumw2();
   HSCP_iso_eop->Sumw2(); 
   nPV->Sumw2(); 
   HSCP_invB->Sumw2(); 
   HSCP_errinvB->Sumw2(); 
   HSCP_invBDT->Sumw2(); 
   HSCP_invBCSC->Sumw2(); 
   HSCP_time->Sumw2();
   HSCP_probQ_stdAna->Sumw2();
   HSCP_probXY_stdAna->Sumw2();
   HSCPCandidateMatchedGenID_probQ->Sumw2();
   HSCPCandidateMatchedGenID_probXY->Sumw2();

   R1_StdEdXVsEvent->Sumw2();
   R1_StdEdXVsLumi->Sumw2();
   R1_LumiVsEvent->Sumw2();
   R1_nPVVsEvent->Sumw2();
   R1_CandVsEvent->Sumw2();
   R2_StdEdXVsEvent->Sumw2();
   R2_StdEdXVsLumi->Sumw2();
   R2_LumiVsEvent->Sumw2();
   R2_nPVVsEvent->Sumw2();
   R2_CandVsEvent->Sumw2();
   R3_StdEdXVsEvent->Sumw2();
   R3_StdEdXVsLumi->Sumw2();
   R3_LumiVsEvent->Sumw2();
   R3_nPVVsEvent->Sumw2();
   R3_CandVsEvent->Sumw2();

   TString outputfilename="analysis_rereco";
   if (year==2016) outputfilename+="_2016";
   else if (year==2017) outputfilename+="_2017";
   else if (year==2018) outputfilename+="_2018";
   outputfilename+=Letter;
   outputfilename+=".root";
   TFile* OutputHisto = new TFile(outputfilename,"RECREATE");


   loadSFPixelCalib();
    
   

//   nentries = 100;
   cout << "run on  " << nentries << " entries " << endl;

//    for (Long64_t jentry=89; jentry<90;jentry++) { // this is good for testing
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry%10000 ==0 && jentry!=0) cout << " number of processed events is " << jentry <<  " = " << (100.*jentry)/(1.*nentries) << "%" <<endl;

//      if (runNumber<305150 || runNumber>305300) continue;

     
       
//       for (int igenpart=0; igenpart<ngenpart; igenpart++) {
//           bool selection=true;
//           if (gen_pt[igenpart]<55) selection=false;
//           if (abs(gen_eta[ngenpart])>2.1) selection=false;
//           if (!hlt_mu50 && !hlt_tkmu100 && !hlt_oldmu100) selection=false;
//
//           bool isHSCPBasedOnGenPDGIDs=false;
//           int PDGid = gen_pdg[ngenpart];
//
//           isHSCPBasedOnGenPDGIDs = (PDGid==17||PDGid==1000993||PDGid==1009213||PDGid==1009313||PDGid==1009323||PDGid==1009113||PDGid==1009223||PDGid==1009333||PDGid==1091114||PDGid==1092114||PDGid==1092214||PDGid==1092224||PDGid==1093114||PDGid==1093214||PDGid==1093224||PDGid==1093314||PDGid==1093324||PDGid==1093334);
//
//           cout << "Filling HSCP gen plots at the index " << track_index_hit[igenpart] << " where igenpart is " << igenpart << endl;
//           if (isHSCPBasedOnGenPDGIDs) {
//               if (selection) {
//                   GenHSCP_probQ_wCuts->Fill(track_probQ[track_index_hit[igenpart]]);
//                   GenHSCP_probXY_wCuts->Fill(track_probXY[track_index_hit[igenpart]]);
//               } else {
//                   GenHSCP_probQ_wNoCuts->Fill(track_probQ[track_index_hit[igenpart]]);
//                   GenHSCP_probXY_wNoCuts->Fill(track_probXY[track_index_hit[igenpart]]);
//               }
//           } else {
//               if (selection) {
//                   GenBckg_probQ_wCuts->Fill(track_probQ[track_index_hit[igenpart]]);
//                   GenBckg_probXY_wCuts->Fill(track_probXY[track_index_hit[igenpart]]);
//               } else {
//                   GenBckg_probQ_wNoCuts->Fill(track_probQ[track_index_hit[igenpart]]);
//                   GenBckg_probXY_wNoCuts->Fill(track_probXY[track_index_hit[igenpart]]);
//               }
//           }
//       } // end loop on gen part
       
      if (!hlt_mu50 && !hlt_tkmu100 && !hlt_oldmu100) continue;
       
      nPVVsRun->Fill(runNumber,npv);
      nPV->Fill(npv);
       
       // start loop on tracks
       for (int iTrack=0; iTrack<ntracks; iTrack++) {
           bool selection=true;
           bool isHSCPBasedOnGenPDGIDs=false;

           if (track_pt[iTrack]<55) selection=false;
           if (track_nvalidhits[iTrack]<8) selection=false;
           if (track_chi2[iTrack]>5) selection=false;
           if (abs(track_eta[iTrack])>2.1) selection=false;
           if (track_npixhits[iTrack]<2) selection=false;
           if (track_qual[iTrack]<2) selection=false;
           if (track_pterr[iTrack]/track_pt[iTrack]>0.25) selection=false;

//           cout << "----------------------------------------" << endl;
//           cout << "Lookin at track index iTrack=" << iTrack << endl;
//           cout << "track_probQ[iTrack] before selection: " << track_probQ[1] << track_probQ[2] << track_probQ[3] << track_probQ[4] << endl;
           
           int numOfClusOnTrack = track_nclus[iTrack];

           for (int ihs=0; ihs<numOfClusOnTrack; ihs++) {
               int PDGid = track_clus_PID[iTrack][ihs];
               isHSCPBasedOnGenPDGIDs = (PDGid==17||PDGid==1000993||PDGid==1009213||PDGid==1009313||PDGid==1009323||PDGid==1009113||PDGid==1009223||PDGid==1009333||PDGid==1091114||PDGid==1092114||PDGid==1092214||PDGid==1092224||PDGid==1093114||PDGid==1093214||PDGid==1093224||PDGid==1093314||PDGid==1093324||PDGid==1093334);

//               cout << "Looking at HSCP candidate ihs=" << ihs << endl;
//               cout << "on the track index hscp_track_idx=" << hscp_track_idx << endl;
//               cout << "The HSCP PDGid is " << PDGid << endl;
//               cout << "out of the total " << nhscp << " HSCP candidates on this track" << endl;
               if (isHSCPBasedOnGenPDGIDs) {
                   if (selection) {
                       GenHSCP_probQ_wCuts->Fill(track_probQ[iTrack]);
                       GenHSCP_probXY_wCuts->Fill(track_probXY[iTrack]);
                   } else {
                       GenHSCP_probQ_wNoCuts->Fill(track_probQ[iTrack]);
                       GenHSCP_probXY_wNoCuts->Fill(track_probXY[iTrack]);
                   }
               } else {
                   if (selection) {
                       GenBckg_probQ_wCuts->Fill(track_probQ[iTrack]);
                       GenBckg_probXY_wCuts->Fill(track_probXY[iTrack]);
                   } else {
                       GenBckg_probQ_wNoCuts->Fill(track_probQ[iTrack]);
                       GenBckg_probXY_wNoCuts->Fill(track_probXY[iTrack]);
                   }
               }
           } // end loop on rechtis on each track
       } // end loop on tracks
       
      // start loop on HSCP candidates
      for (int ihs=0; ihs<nhscp; ihs++) {
          if (hscp_track_idx[ihs]>-1) { 
             // une trace existe pour le candidat HSCP
             bool selection=true;
             bool isHSCPBasedOnGenPDGIDs=false;
             int PDGid = hscp_gen_id[ihs];
             isHSCPBasedOnGenPDGIDs = (PDGid==17||PDGid==1000993||PDGid==1009213||PDGid==1009313||PDGid==1009323||PDGid==1009113||PDGid==1009223||PDGid==1009333||PDGid==1091114||PDGid==1092114||PDGid==1092214||PDGid==1092224||PDGid==1093114||PDGid==1093214||PDGid==1093224||PDGid==1093314||PDGid==1093324||PDGid==1093334);
              
             if (track_pt[hscp_track_idx[ihs]]<55) selection=false;
             if (track_nvalidhits[hscp_track_idx[ihs]]<8) selection=false;
             if (track_chi2[hscp_track_idx[ihs]]>5) selection=false;
             if (abs(track_eta[hscp_track_idx[ihs]])>2.1) selection=false;
             if (track_npixhits[hscp_track_idx[ihs]]<2) selection=false;
             if (track_qual[hscp_track_idx[ihs]]<2) selection=false;
             float eop=(hscp_iso2_ecal[ihs] + hscp_iso2_hcal[ihs])/track_p[hscp_track_idx[ihs]];
             if (track_pterr[hscp_track_idx[ihs]]/track_pt[hscp_track_idx[ihs]]>0.25) selection=false;
              
             
//             cout << "[hscp_track_idx[ihs]] before selection: " << hscp_track_idx[ihs] << endl;
//             cout << "track_probQ[hscp_track_idx[ihs]] before selection: " << track_probQ[hscp_track_idx[ihs]] << endl;
            
             if (selection) {
                ptVsRun->Fill(runNumber,track_pt[hscp_track_idx[ihs]]);
                HSCP_pt->Fill(track_pt[hscp_track_idx[ihs]]);
                HSCP_iso_eop->Fill(eop);
                HSCP_probQ_stdAna->Fill(track_probQ[hscp_track_idx[ihs]]);
                HSCP_probXY_stdAna->Fill(track_probXY[hscp_track_idx[ihs]]);
                 
//                 cout << "[hscp_track_idx[ihs]] after selection: " << hscp_track_idx[ihs] << endl;
//                 cout << "track_probQ[hscp_track_idx[ihs]] after selection: " << track_probQ[hscp_track_idx[ihs]] << endl;
                 
                 if (isHSCPBasedOnGenPDGIDs) {
                     HSCPCandidateMatchedGenID_probQ->Fill(track_probQ[hscp_track_idx[ihs]]);
                     HSCPCandidateMatchedGenID_probXY->Fill(track_probXY[hscp_track_idx[ihs]]);
                 }
                 

                if (boolILumi) lumiVsRun->Fill(runNumber,InstLumi);

                std::vector <float>  charge_uncorr;
                std::vector <float>  pathlength;
                std::vector <int>    subdetId;
                std::vector <int>    moduleGeometry;
                std::vector <bool>   bool_cleaning;
                std::vector <bool>   mustBeInside;

                std::vector <float>  charge_uncorr1;
                std::vector <float>  pathlength1;
                std::vector <int>    subdetId1;
                std::vector <UInt_t> detId1;
                std::vector <int>    moduleGeometry1;
                std::vector <bool>   bool_cleaning1;
                std::vector <bool>   mustBeInside1;

                std::vector <float> charge_uncorr2;
                std::vector <float> pathlength2;
                std::vector <int> subdetId2;
                std::vector <int> moduleGeometry2;
                std::vector <bool> bool_cleaning2;
                std::vector <bool> mustBeInside2;
                for (int iclu=track_index_hit[hscp_track_idx[ihs]]; iclu<track_index_hit[hscp_track_idx[ihs]]+track_nhits[hscp_track_idx[ihs]]; iclu++) {
                     float ch1=dedx_charge[iclu];
                     bool clean1=true;
                     if (dedx_subdetid[iclu]>=3) {
                        // strip
                        ch1=sclus_charge[iclu];
                        clean1=sclus_clusclean[iclu];
                     }
                     else {
                        // pixel
                        float scaling =GetSFPixel(dedx_subdetid[iclu], dedx_detid[iclu], year, abs(track_eta[hscp_track_idx[ihs]]), runNumber);
                        ch1*=scaling;
                     }
                     charge_uncorr.push_back(ch1);
                     pathlength.push_back(dedx_pathlength[iclu]);
                     subdetId.push_back(dedx_subdetid[iclu]);
                     moduleGeometry.push_back(dedx_modulgeom[iclu]);
                     mustBeInside.push_back(dedx_insideTkMod[iclu]);
                     bool_cleaning.push_back(clean1);

                     if (dedx_isstrip[iclu]) {
                     charge_uncorr1.push_back(ch1);
                     pathlength1.push_back(dedx_pathlength[iclu]);
                     subdetId1.push_back(dedx_subdetid[iclu]);
                     detId1.push_back(dedx_detid[iclu]);
                     moduleGeometry1.push_back(dedx_modulgeom[iclu]);
                     mustBeInside1.push_back(dedx_insideTkMod[iclu]);
                     bool_cleaning1.push_back(clean1);
                     }
                     else {
                     charge_uncorr2.push_back(ch1);
                     pathlength2.push_back(dedx_pathlength[iclu]);
                     subdetId2.push_back(dedx_subdetid[iclu]);
                     moduleGeometry2.push_back(dedx_modulgeom[iclu]);
                     mustBeInside2.push_back(dedx_insideTkMod[iclu]);
                     bool_cleaning2.push_back(clean1);
                     }
                    }

                    int nval1=0;
                    int nval2=0;
                    int nval3=0;
                    int nval1_0=0;
                    int nval2_0=0;
                    int nval3_0=0;
                    int nsatv1=0;
                    int nsatv2=0;
                    int nsatv3=0;
                    int nsatv1_0=0;
                    int nsatv2_0=0;
                    int nsatv3_0=0;
                    double ih_uncor = getdEdX(charge_uncorr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,2, 0.15,  nval1, nsatv1);
                    double ih_strip = getdEdX(charge_uncorr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0.15, nval2, nsatv2);
                    double ih_pix = getdEdX(charge_uncorr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0.15, nval3, nsatv3);

                    float fmip_strip4 = FMIP(charge_uncorr1, pathlength1, 4);
                    float fmip_strip3p5 = FMIP(charge_uncorr1, pathlength1, 3.5);
                    float fmip_strip3p2 = FMIP(charge_uncorr1, pathlength1, 3.2);

                    double ih0_uncor = getdEdX(charge_uncorr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,2, 0., nval1_0, nsatv1_0);
                    double ih0_strip = getdEdX(charge_uncorr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0., nval2_0, nsatv2_0);
                    double ih0_pix = getdEdX(charge_uncorr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0.,nval3_0, nsatv3_0);
                    dEdXVsRun->Fill(runNumber,ih_uncor);
                    dEdXpixVsRun->Fill(runNumber,ih_pix);
                    dEdXstripVsRun->Fill(runNumber,ih_strip);
                    FMIP4VsRun->Fill(runNumber,fmip_strip4);
                    FMIP3p5VsRun->Fill(runNumber,fmip_strip3p5);
                    FMIP3p2VsRun->Fill(runNumber,fmip_strip3p2);
                    FMIP4VsEta->Fill(track_eta[hscp_track_idx[ihs]],fmip_strip4);

                    NmeasVsRun->Fill(runNumber,nval1);
                    NmeasPixVsRun->Fill(runNumber,nval3);
                    NmeasStrVsRun->Fill(runNumber,nval2);
                    Nmeas0VsRun->Fill(runNumber,nval1_0);
                    NmeasPix0VsRun->Fill(runNumber,nval3_0);
                    NmeasStr0VsRun->Fill(runNumber,nval2_0);
                    if (nsatv1>9) nsatv1=9;
                    if (nsatv2>9) nsatv2=9;
                    if (nsatv3>9) nsatv3=9;
                    if (nsatv1_0>9) nsatv1_0=9;
                    if (nsatv2_0>9) nsatv2_0=9;
                    if (nsatv3_0>9) nsatv3_0=9;
                    NsatVsRun->Fill(runNumber,nsatv1);
                    NsatPixVsRun->Fill(runNumber,nsatv3);
                    NsatStrVsRun->Fill(runNumber,nsatv2);
                    Nsat0VsRun->Fill(runNumber,nsatv1_0);
                    NsatPix0VsRun->Fill(runNumber,nsatv3_0);
                    NsatStr0VsRun->Fill(runNumber,nsatv2_0);

/*
                    dEdXVsIL->Fill(InstLumi,ih_uncor);
                    dEdXpixVsIL->Fill(InstLumi,ih_pix);
                    dEdXstripVsIL->Fill(InstLumi,ih_strip);
*/
                    HSCP_dEdX->Fill(ih_uncor);
                    HSCP_FMIP4->Fill(fmip_strip4);
                    HSCP_FMIP3p5->Fill(fmip_strip3p5);
                    HSCP_FMIP3p2->Fill(fmip_strip3p2);
                    HSCP_dEdXpix->Fill(ih_pix);
                    HSCP_dEdXstrip->Fill(ih_strip);
                    dEdX0VsRun->Fill(runNumber,ih0_uncor);
                    dEdX0pixVsRun->Fill(runNumber,ih0_pix);
                    dEdX0stripVsRun->Fill(runNumber,ih0_strip);
                    HSCP_dEdX0->Fill(ih0_uncor);
                    HSCP_dEdX0pix->Fill(ih0_pix);
                    HSCP_dEdX0strip->Fill(ih0_strip);

                    if (runNumber==305186) {
                        R1_StdEdXVsEvent->Fill(event,ih_strip);
                        if (boolILumi) {
                         R1_StdEdXVsLumi->Fill(InstLumi,ih_strip);
                         R1_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R1_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R1_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R1_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
                    }
                    else if (runNumber==305188) {
                        R2_StdEdXVsEvent->Fill(event,ih_strip);
                        if (boolILumi) {
                         R2_StdEdXVsLumi->Fill(InstLumi,ih_strip);
                         R2_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R2_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R2_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R2_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
                    }
                    else if (runNumber==305204) {
                        R3_StdEdXVsEvent->Fill(event,ih_strip);
                        if (boolILumi) {
                         R3_StdEdXVsLumi->Fill(InstLumi,ih_strip);
                         R3_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R3_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R3_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R3_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
		    }
                    

                    double ih4_uncor = getdEdX(charge_uncorr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,4, 0.15,  nval1, nsatv1);
                    double ih4_strip = getdEdX(charge_uncorr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,4, 0.15, nval2, nsatv2);
                    double ih4_pix = getdEdX(charge_uncorr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,4, 0.15, nval3, nsatv3);

                    double ih40_uncor = getdEdX(charge_uncorr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,4, 0., nval1_0, nsatv1_0);
                    double ih40_strip = getdEdX(charge_uncorr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,4, 0., nval2_0, nsatv2_0);
                    double ih40_pix = getdEdX(charge_uncorr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,4, 0.,nval3_0, nsatv3_0);

                    dEdX4VsRun->Fill(runNumber,ih4_uncor);
                    dEdX4pixVsRun->Fill(runNumber,ih4_pix);
                    dEdX4stripVsRun->Fill(runNumber,ih4_strip);
                    dEdX40VsRun->Fill(runNumber,ih40_uncor);
                    dEdX40pixVsRun->Fill(runNumber,ih40_pix);
                    dEdX40stripVsRun->Fill(runNumber,ih40_strip);

                    if (hscp_muon_idx[ihs]>-1) { 
                       // un muon existe pour le candidat HSCP
                       if (muon_comb_tofndof[hscp_muon_idx[ihs]]>=8
                             && (muon_dt_tofndof[hscp_muon_idx[ihs]]>=6 || muon_csc_tofndof[hscp_muon_idx[ihs]]>=6) 
                             && muon_comb_inversebetaerr[hscp_muon_idx[ihs]]<=0.15 
                             && fabs(muon_comb_inversebeta[hscp_muon_idx[ihs]]-1)<50) {

                         invBVsRun->Fill(runNumber,muon_comb_inversebeta[hscp_muon_idx[ihs]]);
                         errinvBVsRun->Fill(runNumber,muon_comb_inversebetaerr[hscp_muon_idx[ihs]]);
                         timeVsRun->Fill(runNumber,muon_comb_vertextime[hscp_muon_idx[ihs]]);
                         HSCP_invB->Fill(muon_comb_inversebeta[hscp_muon_idx[ihs]]);
                         HSCP_errinvB->Fill(muon_comb_inversebetaerr[hscp_muon_idx[ihs]]);
                         HSCP_time->Fill(muon_comb_vertextime[hscp_muon_idx[ihs]]);
                         if (muon_dt_tofndof[hscp_muon_idx[ihs]]>=6) {
                         invBDTVsRun->Fill(runNumber,muon_dt_inversebeta[hscp_muon_idx[ihs]]);
                         HSCP_invBDT->Fill(muon_dt_inversebeta[hscp_muon_idx[ihs]]);
                         }
                         if (muon_csc_tofndof[hscp_muon_idx[ihs]]>=6) {
                         invBCSCVsRun->Fill(runNumber,muon_csc_inversebeta[hscp_muon_idx[ihs]]);
                         HSCP_invBCSC->Fill(muon_csc_inversebeta[hscp_muon_idx[ihs]]);
                         }
                       }
                       if (year==2016) {
                        if (muon_newcomb_tofndof[hscp_muon_idx[ihs]]>=8
                             && (muon_newdt_tofndof[hscp_muon_idx[ihs]]>=6 || muon_newcsc_tofndof[hscp_muon_idx[ihs]]>=6) 
                             && muon_newcomb_inversebetaerr[hscp_muon_idx[ihs]]<=0.15 
                             && fabs(muon_newcomb_inversebeta[hscp_muon_idx[ihs]]-1)<50) {

                         invBnewVsRun->Fill(runNumber,muon_newcomb_inversebeta[hscp_muon_idx[ihs]]);
                         if (muon_newdt_tofndof[hscp_muon_idx[ihs]]>=6) {
                         invBnewDTVsRun->Fill(runNumber,muon_newdt_inversebeta[hscp_muon_idx[ihs]]);
                         }
                         if (muon_newcsc_tofndof[hscp_muon_idx[ihs]]>=6) {
                         invBnewCSCVsRun->Fill(runNumber,muon_newcsc_inversebeta[hscp_muon_idx[ihs]]);
                         }
                        }
                       }
                         
                    }
             }
          }

      } // end loop HSCP candidates

      int ntracks1=0;
      int ntracks20=0;
      int ntracks50=0;
      int ntracks100=0;
      for (int itr=0; itr<ntracks; itr++) {
         int presk= 1;
         if (year!=2016) presk=track_prescale[itr];
         if (track_chi2[itr]>5) continue;
         if (track_pt[itr]>100) Htracknhit->Fill(track_nvalidhits[itr]);

         if (track_nvalidhits[itr]<8) continue;
         Htrackpt->Fill(track_pt[itr]);

         if (track_pt[itr]>1) ntracks1+=presk;
         if (track_pt[itr]>20) ntracks20+=presk;
         if (track_pt[itr]>50) ntracks50+=presk;
         if (track_pt[itr]>100) ntracks100+=presk;


         std::vector <float> charge_uncorr;
         std::vector <float> pathlength;
         std::vector <int> subdetId;
         std::vector <UInt_t> detId;
         std::vector <int> moduleGeometry;
         std::vector <bool> bool_cleaning;
         std::vector <bool> mustBeInside;

         std::vector <float> charge_uncorr1;
         std::vector <float> pathlength1;
         std::vector <int> subdetId1;
         std::vector <UInt_t> detId1;
         std::vector <int> moduleGeometry1;
         std::vector <bool> bool_cleaning1;
         std::vector <bool> mustBeInside1;

         std::vector <float> charge_uncorr2;
         std::vector <float> pathlength2;
         std::vector <int> subdetId2;
         std::vector <int> moduleGeometry2;
         std::vector <bool> bool_cleaning2;
         std::vector <bool> mustBeInside2;
         for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
              float ch1=dedx_charge[iclu];
              bool clean1=true;
              if (dedx_subdetid[iclu]>=3) {
                 // strip
                 ch1=sclus_charge[iclu];
                 clean1=sclus_clusclean[iclu];
              }
              else {
                 // pixel
                 float scaling =GetSFPixel(dedx_subdetid[iclu], dedx_detid[iclu], year, abs(track_eta[itr]), runNumber);
                 ch1*=scaling;
              }
              charge_uncorr.push_back(ch1);
              pathlength.push_back(dedx_pathlength[iclu]);
              subdetId.push_back(dedx_subdetid[iclu]);
              detId.push_back(dedx_detid[iclu]);
              moduleGeometry.push_back(dedx_modulgeom[iclu]);
              mustBeInside.push_back(dedx_insideTkMod[iclu]);
              bool_cleaning.push_back(clean1);

              if (dedx_isstrip[iclu]) {
              charge_uncorr1.push_back(ch1);
              pathlength1.push_back(dedx_pathlength[iclu]);
              subdetId1.push_back(dedx_subdetid[iclu]);
              detId1.push_back(dedx_detid[iclu]);
              moduleGeometry1.push_back(dedx_modulgeom[iclu]);
              mustBeInside1.push_back(dedx_insideTkMod[iclu]);
              bool_cleaning1.push_back(clean1);
              if (track_p[itr]<5) Htrackdedx_strip_lowp->Fill(ch1,presk);
              if (track_pt[itr]>100) Htrackdedx_strip->Fill(ch1, presk);
              }
              else {
              charge_uncorr2.push_back(ch1);
              pathlength2.push_back(dedx_pathlength[iclu]);
              subdetId2.push_back(dedx_subdetid[iclu]);
              moduleGeometry2.push_back(dedx_modulgeom[iclu]);
              mustBeInside2.push_back(dedx_insideTkMod[iclu]);
              bool_cleaning2.push_back(clean1);
              if (track_p[itr]<5) Htrackdedx_pix_lowp->Fill(ch1,presk);
              if (track_pt[itr]>100) Htrackdedx_pix->Fill(ch1,presk);
              }
         }

         int nv=0;
         int ns=0;
         double ih_uncor = getdEdX(charge_uncorr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF,  NULL,2, 0.15,  nv, ns);
         double ih_strip = getdEdX(charge_uncorr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF,  NULL,2, 0.15,  nv, ns);
         double ih_pix = getdEdX(charge_uncorr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF,  NULL,2, 0.15,  nv, ns);
         if (track_p[itr]<20) {
//            dEdXVsP_lowp2->Fill(track_p[itr],track_ih_ampl[itr],presk);
//            dEdXVsP_lowp3->Fill(track_p[itr],track_ih_ampl[itr],presk);
            dEdXVsP_lowp2->Fill(track_p[itr],ih_uncor,presk);
            dEdXVsP_lowp3->Fill(track_p[itr],ih_uncor,presk);
         }
         if (track_p[itr]<5) {
//          Htrackih_lowp->Fill(track_ih_ampl[itr], presk);
//          dEdXVsP_lowp->Fill(track_p[itr],track_ih_ampl[itr],presk);
          Htrackih_lowp->Fill(ih_uncor, presk);
          dEdXVsP_lowp->Fill(track_p[itr],ih_uncor,presk);
          Htrackih_reco_lowp->Fill(ih_uncor,presk);
          Htrackih_pix_lowp->Fill(ih_pix,presk);
          Htrackih_strip_lowp->Fill(ih_strip,presk);
          dEdXpixVsP_lowp->Fill(track_p[itr],ih_pix,presk);
          dEdXstripVsP_lowp->Fill(track_p[itr],ih_strip,presk);
          dEstrVsdE_lowp->Fill(ih_uncor,ih_strip,presk);
          dEdXstripVsEta_lowp->Fill(track_eta[itr],ih_strip,presk);
         } 
         else {
          std::vector<double> vect;
          std::vector <int> subdetvec;
          std::vector <UInt_t> detvec;
//          for (int h=0; h<charge_uncorr1.size(); h++) {   //only strip
          for (int h=0; h<charge_uncorr.size(); h++) {      //all cluster
            if(!bool_cleaning[h])continue;
            if(!mustBeInside[h])continue;
            double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
            double ChargeOverPathlength = Norm*charge_uncorr[h]/pathlength[h];
            vect.push_back(ChargeOverPathlength); //save charge
            subdetvec.push_back(subdetId[h]);
            detvec.push_back(detId[h]);
          }
          std::vector <double> tmp (vect.size());
          std::copy (vect.begin(), vect.end(), tmp.begin());
//          std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
//          int nTrunc = tmp.size()*0.15;
//          dEdXstripVsNhit_lowp->Fill(vect.size(),ih_strip,presk);
//          for(unsigned int t=0;t+nTrunc<tmp.size();t++){
          for(unsigned int t=0;t<tmp.size();t++){
/*
              if (subdetvec[t]>2) {
                dEdXstripVsCharge_lowp->Fill(tmp[t],ih_strip,presk);
                if (ih_strip<3.4) Htrackdedx_strip_lowp1->Fill(tmp[t],presk);
                else Htrackdedx_strip_lowp2->Fill(tmp[t],presk);
              }
*/
              int lab=0;
/*
              int indexx=-1;
              for (int h=0; h<vect.size(); h++) {
                float dif = tmp[t]-vect[h];
                if (dif<0) dif*=-1;
                if (dif<0.000001) indexx=h;
              }
//              if (indexx>-1) lab=GetLayerLabel(subdetvec[t], detvec[t]);
              if (indexx>-1) lab=GetLayerLabel(subdetvec[indexx], detvec[indexx]);
*/

              lab=GetLayerLabel(subdetvec[t], detvec[t],year);
              if (lab==1)      Charge_tib1->Fill(tmp[t],presk);
              else if (lab==2) Charge_tib2->Fill(tmp[t],presk);
              else if (lab==3) Charge_tib3->Fill(tmp[t],presk);
              else if (lab==4) Charge_tib4->Fill(tmp[t],presk);
              else if (lab==5) Charge_tob1->Fill(tmp[t],presk);
              else if (lab==6) Charge_tob2->Fill(tmp[t],presk);
              else if (lab==7) Charge_tob3->Fill(tmp[t],presk);
              else if (lab==8) Charge_tob4->Fill(tmp[t],presk);
              else if (lab==9) Charge_tob5->Fill(tmp[t],presk);
              else if (lab==10) Charge_tob6->Fill(tmp[t],presk);
              else if (lab==11) Charge_tid1->Fill(tmp[t],presk);
              else if (lab==12) Charge_tid2->Fill(tmp[t],presk);
              else if (lab==13) Charge_tid3->Fill(tmp[t],presk);
              else if (lab==14) Charge_tec1->Fill(tmp[t],presk);
              else if (lab==15) Charge_tec2->Fill(tmp[t],presk);
              else if (lab==16) Charge_tec3->Fill(tmp[t],presk);
              else if (lab==17) Charge_tec4->Fill(tmp[t],presk);
              else if (lab==18) Charge_tec5->Fill(tmp[t],presk);
              else if (lab==19) Charge_tec6->Fill(tmp[t],presk);
              else if (lab==20) Charge_tec7->Fill(tmp[t],presk);
              else if (lab==21) Charge_tec8->Fill(tmp[t],presk);
              else if (lab==22) Charge_tec9->Fill(tmp[t],presk);
              else if (lab==23) Charge_pixl1->Fill(tmp[t],presk);
              else if (lab==24) Charge_pixl2->Fill(tmp[t],presk);
              else if (lab==25) Charge_pixl3->Fill(tmp[t],presk);
              else if (lab==26) Charge_pixl4->Fill(tmp[t],presk);
              else if (lab==27) Charge_pixd1->Fill(tmp[t],presk);
              else if (lab==28) Charge_pixd2->Fill(tmp[t],presk);
              else if (lab==29) Charge_pixd3->Fill(tmp[t],presk);
/*
              else {
                cout << " lab " << subdetvec[t] << "  "  << detvec[t] ;
                cout << "   "  << lab << endl;
              }
*/

              if (lab==1)      ChargeVsRun_tib1->Fill(runNumber,tmp[t],presk);
              else if (lab==2) ChargeVsRun_tib2->Fill(runNumber,tmp[t],presk);
              else if (lab==3) ChargeVsRun_tib3->Fill(runNumber,tmp[t],presk);
              else if (lab==4) ChargeVsRun_tib4->Fill(runNumber,tmp[t],presk);
              else if (lab==5) ChargeVsRun_tob1->Fill(runNumber,tmp[t],presk);
              else if (lab==6) ChargeVsRun_tob2->Fill(runNumber,tmp[t],presk);
              else if (lab==7) ChargeVsRun_tob3->Fill(runNumber,tmp[t],presk);
              else if (lab==8) ChargeVsRun_tob4->Fill(runNumber,tmp[t],presk);
              else if (lab==9) ChargeVsRun_tob5->Fill(runNumber,tmp[t],presk);
              else if (lab==10) ChargeVsRun_tob6->Fill(runNumber,tmp[t],presk);
              else if (lab==11) ChargeVsRun_tid1->Fill(runNumber,tmp[t],presk);
              else if (lab==12) ChargeVsRun_tid2->Fill(runNumber,tmp[t],presk);
              else if (lab==13) ChargeVsRun_tid3->Fill(runNumber,tmp[t],presk);
              else if (lab==14) ChargeVsRun_tec1->Fill(runNumber,tmp[t],presk);
              else if (lab==15) ChargeVsRun_tec2->Fill(runNumber,tmp[t],presk);
              else if (lab==16) ChargeVsRun_tec3->Fill(runNumber,tmp[t],presk);
              else if (lab==17) ChargeVsRun_tec4->Fill(runNumber,tmp[t],presk);
              else if (lab==18) ChargeVsRun_tec5->Fill(runNumber,tmp[t],presk);
              else if (lab==19) ChargeVsRun_tec6->Fill(runNumber,tmp[t],presk);
              else if (lab==20) ChargeVsRun_tec7->Fill(runNumber,tmp[t],presk);
              else if (lab==21) ChargeVsRun_tec8->Fill(runNumber,tmp[t],presk);
              else if (lab==22) ChargeVsRun_tec9->Fill(runNumber,tmp[t],presk);
              else if (lab==23) ChargeVsRun_pixl1->Fill(runNumber,tmp[t],presk);
              else if (lab==24) ChargeVsRun_pixl2->Fill(runNumber,tmp[t],presk);
              else if (lab==25) ChargeVsRun_pixl3->Fill(runNumber,tmp[t],presk);
              else if (lab==26) ChargeVsRun_pixl4->Fill(runNumber,tmp[t],presk);
              else if (lab==27) ChargeVsRun_pixd1->Fill(runNumber,tmp[t],presk);
              else if (lab==28) ChargeVsRun_pixd2->Fill(runNumber,tmp[t],presk);
              else if (lab==29) ChargeVsRun_pixd3->Fill(runNumber,tmp[t],presk);

              if (lab==1)      ZooChargeVsRun_tib1->Fill(runNumber,tmp[t],presk);
              else if (lab==2) ZooChargeVsRun_tib2->Fill(runNumber,tmp[t],presk);
              else if (lab==3) ZooChargeVsRun_tib3->Fill(runNumber,tmp[t],presk);
              else if (lab==4) ZooChargeVsRun_tib4->Fill(runNumber,tmp[t],presk);
              else if (lab==5) ZooChargeVsRun_tob1->Fill(runNumber,tmp[t],presk);
              else if (lab==6) ZooChargeVsRun_tob2->Fill(runNumber,tmp[t],presk);
              else if (lab==7) ZooChargeVsRun_tob3->Fill(runNumber,tmp[t],presk);
              else if (lab==8) ZooChargeVsRun_tob4->Fill(runNumber,tmp[t],presk);
              else if (lab==9) ZooChargeVsRun_tob5->Fill(runNumber,tmp[t],presk);
              else if (lab==10) ZooChargeVsRun_tob6->Fill(runNumber,tmp[t],presk);
              else if (lab==11) ZooChargeVsRun_tid1->Fill(runNumber,tmp[t],presk);
              else if (lab==12) ZooChargeVsRun_tid2->Fill(runNumber,tmp[t],presk);
              else if (lab==13) ZooChargeVsRun_tid3->Fill(runNumber,tmp[t],presk);
              else if (lab==14) ZooChargeVsRun_tec1->Fill(runNumber,tmp[t],presk);
              else if (lab==15) ZooChargeVsRun_tec2->Fill(runNumber,tmp[t],presk);
              else if (lab==16) ZooChargeVsRun_tec3->Fill(runNumber,tmp[t],presk);
              else if (lab==17) ZooChargeVsRun_tec4->Fill(runNumber,tmp[t],presk);
              else if (lab==18) ZooChargeVsRun_tec5->Fill(runNumber,tmp[t],presk);
              else if (lab==19) ZooChargeVsRun_tec6->Fill(runNumber,tmp[t],presk);
              else if (lab==20) ZooChargeVsRun_tec7->Fill(runNumber,tmp[t],presk);
              else if (lab==21) ZooChargeVsRun_tec8->Fill(runNumber,tmp[t],presk);
              else if (lab==22) ZooChargeVsRun_tec9->Fill(runNumber,tmp[t],presk);
              else if (lab==23) ZooChargeVsRun_pixl1->Fill(runNumber,tmp[t],presk);
              else if (lab==24) ZooChargeVsRun_pixl2->Fill(runNumber,tmp[t],presk);
              else if (lab==25) ZooChargeVsRun_pixl3->Fill(runNumber,tmp[t],presk);
              else if (lab==26) ZooChargeVsRun_pixl4->Fill(runNumber,tmp[t],presk);
              else if (lab==27) ZooChargeVsRun_pixd1->Fill(runNumber,tmp[t],presk);
              else if (lab==28) ZooChargeVsRun_pixd2->Fill(runNumber,tmp[t],presk);
              else if (lab==29) ZooChargeVsRun_pixd3->Fill(runNumber,tmp[t],presk);

              if (ih_strip<3.4) {
              if (lab==1)      LowCharge_tib1->Fill(tmp[t],presk);
              else if (lab==2) LowCharge_tib2->Fill(tmp[t],presk);
              else if (lab==3) LowCharge_tib3->Fill(tmp[t],presk);
              else if (lab==4) LowCharge_tib4->Fill(tmp[t],presk);
              else if (lab==5) LowCharge_tob1->Fill(tmp[t],presk);
              else if (lab==6) LowCharge_tob2->Fill(tmp[t],presk);
              else if (lab==7) LowCharge_tob3->Fill(tmp[t],presk);
              else if (lab==8) LowCharge_tob4->Fill(tmp[t],presk);
              else if (lab==9) LowCharge_tob5->Fill(tmp[t],presk);
              else if (lab==10) LowCharge_tob6->Fill(tmp[t],presk);
              else if (lab==11) LowCharge_tid1->Fill(tmp[t],presk);
              else if (lab==12) LowCharge_tid2->Fill(tmp[t],presk);
              else if (lab==13) LowCharge_tid3->Fill(tmp[t],presk);
              else if (lab==14) LowCharge_tec1->Fill(tmp[t],presk);
              else if (lab==15) LowCharge_tec2->Fill(tmp[t],presk);
              else if (lab==16) LowCharge_tec3->Fill(tmp[t],presk);
              else if (lab==17) LowCharge_tec4->Fill(tmp[t],presk);
              else if (lab==18) LowCharge_tec5->Fill(tmp[t],presk);
              else if (lab==19) LowCharge_tec6->Fill(tmp[t],presk);
              else if (lab==20) LowCharge_tec7->Fill(tmp[t],presk);
              else if (lab==21) LowCharge_tec8->Fill(tmp[t],presk);
              else if (lab==22) LowCharge_tec9->Fill(tmp[t],presk);
              }

          }
//          dEdXstripVsNhittrunc_lowp->Fill(vect.size()-nTrunc,ih_strip,presk);
         }


         int nsatclu=0;
         int nstrip=0;
         int npix=0;
         for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
            if (!dedx_insideTkMod[iclu]) continue;
            if (dedx_isstrip[iclu] && !sclus_clusclean[iclu]) continue;
            if (dedx_isstrip[iclu] && (sclus_sat254[iclu] || sclus_sat255[iclu])) nsatclu++;
            if (dedx_isstrip[iclu]) nstrip++;
            if (dedx_ispixel[iclu]) npix++;
         }

         if (track_p[itr]<5) {
           if (ih_uncor<3.4) NPixVsNStrip_lowIhp->Fill(nstrip,npix,presk);
           else if (ih_uncor>3.5) NPixVsNStrip_highIhp->Fill(nstrip,npix,presk);
          Htracketa_lowp->Fill(track_eta[itr],presk);
         }

         if (track_pt[itr]<100) continue;
         Htracketa->Fill(track_eta[itr]);
         Htrackphi->Fill(track_phi[itr]);
//         Htrackih->Fill(track_ih_ampl[itr]);
         Htrackih->Fill(ih_uncor);
//         dEdXVsP->Fill(track_p[itr],track_ih_ampl[itr]);
         dEdXVsP->Fill(track_p[itr],ih_uncor);
         dEdXpixVsP->Fill(track_p[itr],ih_pix);
         dEdXstripVsP->Fill(track_p[itr],ih_strip);
         if(track_nvalidhits[itr]==8 || track_nvalidhits[itr]==9) EtaVsPhi_nhit8->Fill(track_phi[itr],track_eta[itr]);
         EtaVsPhi_nhit->Fill(track_phi[itr],track_eta[itr]);

         Nsat->Fill(nsatclu);
         NPix->Fill(npix);
         NStrip->Fill(nstrip);

         Htrackih_reco->Fill(ih_uncor);
         Htrackih_pix->Fill(ih_pix);
         Htrackih_strip->Fill(ih_strip);
      } // end loop tracks
      HNtracks->Fill(ntracks);
      HNtracks1->Fill(ntracks1);
      HNtracks20->Fill(ntracks20);
      HNtracks50->Fill(ntracks50);
      HNtracks100->Fill(ntracks100);

   } // end loop num entries
   OutputHisto->cd();
    
    
   GenHSCP_probQ_wCuts->Write();
   GenHSCP_probXY_wCuts->Write();
   GenHSCP_probQ_wNoCuts->Write();
   GenHSCP_probXY_wNoCuts->Write();
   GenBckg_probQ_wCuts->Write();
   GenBckg_probXY_wCuts->Write();
   GenBckg_probQ_wNoCuts->Write();
   GenBckg_probXY_wNoCuts->Write();
    
    
   HNtracks->Write();
   HNtracks1->Write();
   HNtracks20->Write();
   HNtracks50->Write();
   HNtracks100->Write();

   Htrackpt->Write();
   Htracketa->Write();
   Htracketa_lowp->Write();
   Htrackphi->Write();
   Htracknhit->Write();

   Htrackih->Write();
   Htrackih_reco->Write();
   Htrackih_pix->Write();
   Htrackih_strip->Write();
   Htrackdedx_pix->Write();
   Htrackdedx_strip->Write();

   Htrackih_lowp->Write();
   Htrackih_reco_lowp->Write();
   Htrackih_pix_lowp->Write();
   Htrackih_strip_lowp->Write();
   Htrackdedx_pix_lowp->Write();
   Htrackdedx_strip_lowp->Write();
   Htrackdedx_strip_lowp1->Write();
   Htrackdedx_strip_lowp2->Write();


   Nsat->Write();
   NPix->Write();
   NStrip->Write();
   dEdXVsP->Write();
   dEdXVsP_lowp->Write();
   dEdXVsP_lowp2->Write();
   dEdXVsP_lowp3->Write();
   dEdXpixVsP->Write();
   dEdXstripVsP->Write();
   dEdXpixVsP_lowp->Write();
   dEdXstripVsP_lowp->Write();
   dEdXstripVsEta_lowp->Write();
   dEstrVsdE_lowp->Write();
   NPixVsNStrip_highIhp->Write();
   NPixVsNStrip_lowIhp->Write();
   dEdXstripVsNhit_lowp->Write();
   dEdXstripVsNhittrunc_lowp->Write();
   dEdXstripVsCharge_lowp->Write();
   EtaVsPhi_nhit8->Write();
   EtaVsPhi_nhit->Write();

   Charge_pixl1->Write();
   Charge_pixl2->Write();
   Charge_pixl3->Write();
   Charge_pixl4->Write();
   Charge_pixd1->Write();
   Charge_pixd2->Write();
   Charge_pixd3->Write();
   Charge_tib1->Write();
   Charge_tib2->Write();
   Charge_tib3->Write();
   Charge_tib4->Write();
   Charge_tob1->Write();
   Charge_tob2->Write();
   Charge_tob3->Write();
   Charge_tob4->Write();
   Charge_tob5->Write();
   Charge_tob6->Write();
   Charge_tid1->Write();
   Charge_tid2->Write();
   Charge_tid3->Write();
   Charge_tec1->Write();
   Charge_tec2->Write();
   Charge_tec3->Write();
   Charge_tec4->Write();
   Charge_tec5->Write();
   Charge_tec6->Write();
   Charge_tec7->Write();
   Charge_tec8->Write();
   Charge_tec9->Write();

   LowCharge_tib1->Write();
   LowCharge_tib2->Write();
   LowCharge_tib3->Write();
   LowCharge_tib4->Write();
   LowCharge_tob1->Write();
   LowCharge_tob2->Write();
   LowCharge_tob3->Write();
   LowCharge_tob4->Write();
   LowCharge_tob5->Write();
   LowCharge_tob6->Write();
   LowCharge_tid1->Write();
   LowCharge_tid2->Write();
   LowCharge_tid3->Write();
   LowCharge_tec1->Write();
   LowCharge_tec2->Write();
   LowCharge_tec3->Write();
   LowCharge_tec4->Write();
   LowCharge_tec5->Write();
   LowCharge_tec6->Write();
   LowCharge_tec7->Write();
   LowCharge_tec8->Write();
   LowCharge_tec9->Write();

   ChargeVsRun_pixl1->Write();
   ChargeVsRun_pixl2->Write();
   ChargeVsRun_pixl3->Write();
   ChargeVsRun_pixl4->Write();
   ChargeVsRun_pixd1->Write();
   ChargeVsRun_pixd2->Write();
   ChargeVsRun_pixd3->Write();
   ChargeVsRun_tib1->Write();
   ChargeVsRun_tib2->Write();
   ChargeVsRun_tib3->Write();
   ChargeVsRun_tib4->Write();
   ChargeVsRun_tob1->Write();
   ChargeVsRun_tob2->Write();
   ChargeVsRun_tob3->Write();
   ChargeVsRun_tob4->Write();
   ChargeVsRun_tob5->Write();
   ChargeVsRun_tob6->Write();
   ChargeVsRun_tid1->Write();
   ChargeVsRun_tid2->Write();
   ChargeVsRun_tid3->Write();
   ChargeVsRun_tec1->Write();
   ChargeVsRun_tec2->Write();
   ChargeVsRun_tec3->Write();
   ChargeVsRun_tec4->Write();
   ChargeVsRun_tec5->Write();
   ChargeVsRun_tec6->Write();
   ChargeVsRun_tec7->Write();
   ChargeVsRun_tec8->Write();
   ChargeVsRun_tec9->Write();

   ZooChargeVsRun_pixl1->Write();
   ZooChargeVsRun_pixl2->Write();
   ZooChargeVsRun_pixl3->Write();
   ZooChargeVsRun_pixl4->Write();
   ZooChargeVsRun_pixd1->Write();
   ZooChargeVsRun_pixd2->Write();
   ZooChargeVsRun_pixd3->Write();
   ZooChargeVsRun_tib1->Write();
   ZooChargeVsRun_tib2->Write();
   ZooChargeVsRun_tib3->Write();
   ZooChargeVsRun_tib4->Write();
   ZooChargeVsRun_tob1->Write();
   ZooChargeVsRun_tob2->Write();
   ZooChargeVsRun_tob3->Write();
   ZooChargeVsRun_tob4->Write();
   ZooChargeVsRun_tob5->Write();
   ZooChargeVsRun_tob6->Write();
   ZooChargeVsRun_tid1->Write();
   ZooChargeVsRun_tid2->Write();
   ZooChargeVsRun_tid3->Write();
   ZooChargeVsRun_tec1->Write();
   ZooChargeVsRun_tec2->Write();
   ZooChargeVsRun_tec3->Write();
   ZooChargeVsRun_tec4->Write();
   ZooChargeVsRun_tec5->Write();
   ZooChargeVsRun_tec6->Write();
   ZooChargeVsRun_tec7->Write();
   ZooChargeVsRun_tec8->Write();
   ZooChargeVsRun_tec9->Write();


   dEdXVsRun->Write();
   dEdXpixVsRun->Write();
   dEdXstripVsRun->Write();
/*
   dEdXVsIL->Write();
   dEdXpixVsIL->Write();
   dEdXstripVsIL->Write();
*/
   dEdX0VsRun->Write();
   dEdX0pixVsRun->Write();
   dEdX0stripVsRun->Write();
   dEdX4VsRun->Write();
   dEdX4pixVsRun->Write();
   dEdX4stripVsRun->Write();
   dEdX40VsRun->Write();
   dEdX40pixVsRun->Write();
   dEdX40stripVsRun->Write();
   NmeasVsRun->Write();
   NmeasPixVsRun->Write();
   NmeasStrVsRun->Write();
   Nmeas0VsRun->Write();
   NmeasPix0VsRun->Write();
   NmeasStr0VsRun->Write();
   NsatVsRun->Write();
   NsatPixVsRun->Write();
   NsatStrVsRun->Write();
   Nsat0VsRun->Write();
   NsatPix0VsRun->Write();
   NsatStr0VsRun->Write();
   ptVsRun->Write();
   nPVVsRun->Write();
   invBVsRun->Write();
   errinvBVsRun->Write();
   invBDTVsRun->Write();
   invBCSCVsRun->Write();
   invBnewVsRun->Write();
   invBnewDTVsRun->Write();
   invBnewCSCVsRun->Write();
   timeVsRun->Write();
   lumiVsRun->Write();
   HSCP_dEdX->Write(); 
   HSCP_dEdXpix->Write(); 
   HSCP_dEdXstrip->Write(); 
   HSCP_dEdX0->Write(); 
   HSCP_dEdX0pix->Write(); 
   HSCP_dEdX0strip->Write(); 

   HSCP_FMIP4->Write(); 
   HSCP_FMIP3p5->Write(); 
   HSCP_FMIP3p2->Write(); 
   FMIP4VsRun->Write();
   FMIP3p5VsRun->Write();
   FMIP3p5VsRun->Write();
   FMIP4VsEta->Write();

   HSCP_pt->Write(); 
   HSCP_iso_eop->Write(); 
   nPV->Write(); 
   HSCP_invB->Write(); 
   HSCP_errinvB->Write(); 
   HSCP_invBDT->Write(); 
   HSCP_invBCSC->Write(); 
   HSCP_time->Write();
   HSCP_probQ_stdAna->Write();
   HSCP_probXY_stdAna->Write();
   HSCPCandidateMatchedGenID_probQ->Write();
   HSCPCandidateMatchedGenID_probXY->Write();

   R1_StdEdXVsEvent->Write();
   R1_StdEdXVsLumi->Write();
   R1_LumiVsEvent->Write();
   R1_nPVVsEvent->Write();
   R1_CandVsEvent->Write();
   R2_StdEdXVsEvent->Write();
   R2_StdEdXVsLumi->Write();
   R2_LumiVsEvent->Write();
   R2_nPVVsEvent->Write();
   R2_CandVsEvent->Write();
   R3_StdEdXVsEvent->Write();
   R3_StdEdXVsLumi->Write();
   R3_LumiVsEvent->Write();
   R3_nPVVsEvent->Write();
   R3_CandVsEvent->Write();

   OutputHisto->Close();

}



double run2study::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, int & nv, int & ns) {
  double result=-1;
//     double dropLowerDeDxValue=0.15;
     size_t MaxStripNOM=99;
     bool usePixel=true;
     bool useStrip=true;

     std::vector<double> vect;

     bool debugprint=false;
     unsigned int SiStripNOM = 0;
     ns=0;

     for(unsigned int h=0;h<charge.size();h++){
        if (debugprint) std::cout << "look on dedxHits in computedEdx " << h << std::endl;
        if(!usePixel && subdetId[h]<3)continue; // skip pixels
        if(!useStrip && subdetId[h]>=3)continue; // skip strips        
        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;

        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;
        if(useStrip && subdetId[h]>=3 && !mustBeInside[h])continue;
        if(useStrip && subdetId[h]>=3 && ++SiStripNOM > MaxStripNOM) continue; // skip remaining strips, but not pixel

        int ClusterCharge = charge[h];
        if (subdetId[h]>=3 && charge[h]>=254) ns++;

        double scaleFactor = scaleFactors[0];
        if (subdetId[h]<3) scaleFactor *= scaleFactors[1]; // add pixel scaling
        if (debugprint) std::cout << " after SF " << std::endl;

        if(templateHisto){  //save discriminator probability
           double ChargeOverPathlength = scaleFactor*ClusterCharge/(pathlength[h]*10.0*(subdetId[h]<3?265:1));
           int    BinX   = templateHisto->GetXaxis()->FindBin(moduleGeometry[h]);
           int    BinY   = templateHisto->GetYaxis()->FindBin(pathlength[h]*10.0); //*10 because of cm-->mm
           int    BinZ   = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
           double Prob   = templateHisto->GetBinContent(BinX,BinY,BinZ);
           vect.push_back(Prob); //save probability
           if (debugprint) std::cout << " after Prob vect.push_back " << std::endl;
        }else{
           double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
           double ChargeOverPathlength = scaleFactor*Norm*ClusterCharge/pathlength[h];
           vect.push_back(ChargeOverPathlength); //save charge
           if (debugprint) std::cout << " after ChargeOverPathlength vect.push_back " << std::endl;
        }
     }

     if(dropLowerDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
         int nTrunc = tmp.size()*dropLowerDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropLowerDeDxValue " << std::endl;

     int size = vect.size();
     nv = size;

     if(size>0){
        if(templateHisto){  //dEdx discriminator
          //Ias discriminator
          result = 1.0/(12*size);
           std::sort(vect.begin(), vect.end(), std::less<double>() );
           for(int i=1;i<=size;i++){
              result += vect[i-1] * pow(vect[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
           }
           result *= (3.0/size);
           if (debugprint) std::cout << " Ias discriminator " << result << std::endl;
        }else{  //dEdx estimator
           //harmonic2 estimator
           result=0;
//           double expo = -2;
           double expo = -1* n_estim;
           for(int i = 0; i< size; i ++){
              result+=pow(vect[i],expo);
           }
           result = pow(result/size,1./expo);
           if (debugprint) std::cout << " harmonic discriminator " << result << " with expo " << expo << std::endl;
        }
     }else{
        result = -1;
     }
     if (debugprint) std::cout << " ok finished computeDeDx " << std::endl;


  return result;
}

//float run2study::FMIP(const vector<float>& charge, const<vector>& path, float thre = 4){
float run2study::FMIP(std::vector <float> charge, std::vector <float> path, float thre = 4){
   if(charge.size()!=path.size()) return -999; // error

   int nclusters = charge.size();
   if(nclusters==0) return -888;

   int nlow = 0; 
   for(int i=0;i<charge.size();i++){
       //hard-coded conversion factor to go for dEdx in MeV.cm2/g
       // ATTENTION : HERE THE CONVERSION FACTOR IS 247 WHILE ABOVE WE USE 265 !!!!!
       float dEdx = charge[i]*(3.61*pow(10,-9)*247)*1000/path[i];
       if(dEdx<thre) nlow++;
   }
   return nlow*1./nclusters;
}


int run2study::GetLayerLabel(int subdetid_, UInt_t detid_, int year)
{
// from https://github.com/dapparu/HSCP/blob/c69cf1c71dd99289f72ab6d03077915776c85690/src/Cluster.cc
// and https://cmssdt.cern.ch/lxr/source/DataFormats/SiPixelDetId/interface/PXBDetId.h
// https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder
        if(subdetid_==1)
        {
             if (year==2016) {
                if(((detid_>>16)&0xF)==1) return 23;
                else if(((detid_>>16)&0xF)==2) return 24;
                else if(((detid_>>16)&0xF)==3) return 25;
                else if(((detid_>>16)&0xF)==4) return 26;  // do not exist in 2016
             }
             else {
                if(((detid_>>20)&0xF)==1) return 23;
                else if(((detid_>>20)&0xF)==2) return 24;
                else if(((detid_>>20)&0xF)==3) return 25;
                else if(((detid_>>20)&0xF)==4) return 26;
             }

        }
        else if(subdetid_==2)
        {
             if (year==2016) {
                if(((detid_>>16)&0xF)==1) return 27;
                else if(((detid_>>16)&0xF)==2) return 28;
                else if(((detid_>>16)&0xF)==3) return 29; // do not exist in 2016
             }
             else {
                if(((detid_>>18)&0xF)==1) return 27;
                else if(((detid_>>18)&0xF)==2) return 28;
                else if(((detid_>>18)&0xF)==3) return 29;
             }
/*  https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder for 2016
                cout << "  side " << int((detid_>>23)&0x3) <<
                        "  disk " << int((detid_>>16)&0xF) << 
                        "  blade "  << int((detid_>>10)&0x3F) <<
                        "  panel "  << int((detid_>>8)&0x3) <<
                        "  mod "  << int((detid_>>2)&0x3F) << endl;
*/
//                if (((detid_>>16)&0xF)==0) cout << " disk 0 ? " << std::endl;
        }
        else if(subdetid_==3)  // TIB
        {
                if(((detid_>>14)&0x7)==1) return 1;
                else if(((detid_>>14)&0x7)==2) return 2;
                else if(((detid_>>14)&0x7)==3) return 3;
                else if(((detid_>>14)&0x7)==4) return 4;
        }
        else if(subdetid_==5) // TOB
        {
                if(((detid_>>14)&0x7)==1) return 5;
                else if(((detid_>>14)&0x7)==2) return 6;
                else if(((detid_>>14)&0x7)==3) return 7;
                else if(((detid_>>14)&0x7)==4) return 8;
                else if(((detid_>>14)&0x7)==5) return 9;
                else if(((detid_>>14)&0x7)==6) return 10;
        }
        else if(subdetid_==4)  //TID
        {
                if(((detid_>>11)&0x3)==1) return 11;
                else if(((detid_>>11)&0x3)==2) return 12;
                else if(((detid_>>11)&0x3)==3) return 13;
        }
        else if(subdetid_==6) // TEC
        {
                if(((detid_>>14)&0xF)==1) return 14;
                else if(((detid_>>14)&0xF)==2) return 15;
                else if(((detid_>>14)&0xF)==3) return 16;
                else if(((detid_>>14)&0xF)==4) return 17;
                else if(((detid_>>14)&0xF)==5) return 18;
                else if(((detid_>>14)&0xF)==6) return 19;
                else if(((detid_>>14)&0xF)==7) return 20;
                else if(((detid_>>14)&0xF)==8) return 21;
                else if(((detid_>>14)&0xF)==9) return 22;
        }
        return -1;
}


void run2study::loadSFPixelCalib() {

   std::ifstream file_calib;
   file_calib.open("scale_for_cmssw2017.txt");
   icalib=0;
   if (!file_calib) { std::cerr << "cannot open file scale_for_cmssw2017.txt " << std::endl; }
   else
   {
      while (!file_calib.eof () && icalib<calmax) {
       // pix", "layerorside", "ladderorblade", "etaMin", "etaMax", "irunMin", "irunMax", "value
       file_calib >> pixVal[icalib] >> layerSideVal[icalib] >> ladderBladeVal[icalib] >> etaMinVal[icalib] >> etaMaxVal[icalib]
                  >> irunMinVal[icalib] >> irunMaxVal[icalib] >> scaleVal[icalib] ;
       if (icalib<10) { std::cout << pixVal[icalib] << " " << layerSideVal[icalib]  << " "  <<  ladderBladeVal[icalib] ; 
           std::cout << " " << etaMinVal[icalib] << " "  << etaMaxVal[icalib] << " " << irunMinVal[icalib] << " " << irunMaxVal[icalib] ;
           std::cout << " " << scaleVal[icalib] << std::endl ;
       }
       icalib++;
      }
   }
   std::cout << " file_calib : " << icalib << " entries in 2017" << std::endl;
   if (icalib>0) std::cout << " example runMin "<< irunMinVal[0] << " runMax "<< irunMaxVal[0] << " pix " << pixVal[0] << 
                " scale " << scaleVal[0] <<  std::endl;
   file_calib.close ();

   icalib2017=icalib;
   std::ifstream file_calib2;
   file_calib2.open("scale_for_cmssw2018.txt");
   if (!file_calib2) { std::cerr << "cannot open file scale_for_cmssw2018.txt " << std::endl; }
   else
   {  
      while (!file_calib2.eof () && icalib<calmax) {
       // pix", "layerorside", "ladderorblade", "etaMin", "etaMax", "irunMin", "irunMax", "value
       file_calib2 >> pixVal[icalib] >> layerSideVal[icalib] >> ladderBladeVal[icalib] >> etaMinVal[icalib] >> etaMaxVal[icalib]
                  >> irunMinVal[icalib] >> irunMaxVal[icalib] >> scaleVal[icalib] ;
       icalib++;
      }
   }
   icalib2018=icalib-icalib2017;
   std::cout << " file_calib : " << icalib2018 << " entries in 2018" << std::endl;
   if (icalib>icalib2017) std::cout << " example runMin "<< irunMinVal[icalib2017] << " runMax "<< irunMaxVal[icalib2017] 
                 << " pix " << pixVal[icalib2017] << " scale " << scaleVal[icalib2017] <<  std::endl;
   file_calib2.close ();
 



}

float run2study::GetSFPixel(int subdetid_, UInt_t detid_, int year, float eta, int run) {
  // https://github.com/cms-sw/cmssw/tree/master/Geometry/TrackerNumberingBuilder

   int pix = 0;
   int layerorside = 0;
   int ladderorblade = 0;

   if (subdetid_==1) {
        pix =1 ;
        if (year==2016) {
                layerorside = int((detid_>>16)&0xF);
                ladderorblade = int((detid_>>8)&0xFF);
        }
        else {
                layerorside = int((detid_>>20)&0xF);
                ladderorblade = int((detid_>>12)&0xFF);
        }
   }
   if (subdetid_==2) {
        pix =2 ;
        if (year==2016) {
                layerorside = int((detid_>>16)&0xF);
                ladderorblade = int((detid_>>10)&0x3F);
        }
        else {
                layerorside = int((detid_>>18)&0xF);
                ladderorblade = int((detid_>>12)&0x3F);
        }
   }

   float scale =1 ;

   if (year==2017) {
    for (int i=0; i<icalib2017; i++) {
     if (run >= irunMinVal[i] && run<irunMaxVal[i]) {
       if (eta >= etaMinVal[i] && eta < etaMaxVal[i]) { 
          if (pix == pixVal[i] ) {
             if (layerorside == layerSideVal[i]) {
               if (ladderorblade == ladderBladeVal[i]) {
                  scale = scaleVal[i];
                  return scale;
               }
             }
          }
       }
     }
    }
   }
   else if (year==2018) {
    for (int i=icalib2017; i<icalib; i++) {
     if (run >= irunMinVal[i] && run<irunMaxVal[i]) {
       if (eta >= etaMinVal[i] && eta < etaMaxVal[i]) { 
          if (pix == pixVal[i] ) {
             if (layerorside == layerSideVal[i]) {
               if (ladderorblade == ladderBladeVal[i]) {
                  scale = scaleVal[i];
                  return scale;
               }
             }
          }
       }
     }
    }
   }

   return scale;
   
}
