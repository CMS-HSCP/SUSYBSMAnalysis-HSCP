#define run2analysis_cxx
#include "run2analysis.h"
#include <TMatrix.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>
#include <fstream>

enum TrackQuality {
      undefQuality = -1,
      loose = 0,
      tight = 1,
      highPurity = 2,
      confirmed = 3,      // means found by more than one iteration
      goodIterative = 4,  // meaningless
      looseSetWithPV = 5,
      highPuritySetWithPV = 6,
      discarded = 7,  // because a better track found. kept in the collection for reference....
      qualitySize = 8
    };

bool writeTemplateOnDisk=true;
bool computeSpecial= true;
bool boolDeDxTemp= true;

// Modification of the code in September 2021 to 
// align it with the use of xtalk inversion (only for cluster cleaning) and saturation 
//
void run2analysis::Loop(int year, TString Letter, bool dataFlag=true)
{
//   In a ROOT session, you can do:
//      root> .L run2analysis.C
//      root> run2analysis t
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

   Long64_t nentries = fChain->GetEntries();

   Long64_t nbytes = 0, nb = 0;


   double pi=acos(-1);
   double dEdxSF [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF


   // values for 2017 UL data & MC
   if (dataFlag) {
      if (year==2017) {
      dEdxSF[0]=1.;
      dEdxSF[1]=1.0325;
      }
      else if (year==2018) { 
      dEdxSF[0]=1.;
      dEdxSF[1]=1.0817;
      }
      else {
         cout << " AIE AIE AIE : Year not well specified for dEdXSF in data !!!! " << endl;
         return;
      }
   }
   else {
      if (year==2017) {
        dEdxSF[0]=1.0079;
        dEdxSF[1]=1.0875;
      } 
      else if (year==2018) {
        dEdxSF[0]=1.0047;
        dEdxSF[1]=1.1429;
      }
      else {
         cout << " AIE AIE AIE : Year not well specified for dEdXSF in MC !!!! " << endl;
         return;
      }
   }
   cout << "  dEdXSF : " << dEdxSF[0] << "   &   " << dEdxSF[1] << endl;

   // update all on feb 25 (on jan21 rootfile)
   // with Cval extracted on 3-5 
   // MC 2017
   float Kval_ld = 2.33;
   float Cval_ld =3.37 ;
   float Kval_ldstrip = 2.45;
   float Cval_ldstrip = 3.41;
   float Kval_all = 2.19;
   float Cval_all = 3.18;

   float Cval_nol1 = 3.31; // fit C in 5-25
   float Kval_nol1 = 2.20; 
   float Cval_nol1_2 = 3.22; // fit C in 3-5
   float Kval_nol1_2 = 2.26;
   float Cval_nol1_3 = 2.88; // with log term
   float Kval_nol1_3 = 2.57;
   float Nval_nol1_3 = 0.099;

   float Kval_strip = 2.43;
   float Cval_strip = 3.23;
   float Kval_hdnol1 = 2.20; 
   float Cval_hdnol1 = 3.07; 

   float Kval_pix = 1.75;
   float Cval_pix = 3.03;
   float Kval_pixnol1 = 1.62;
   float Cval_pixnol1 = 3.13;

   if (dataFlag) {
    if (year==2017) {
     // extracted on analysis_ul_2017_21jan.root
     Kval_ld = 2.36 ;
     Cval_ld = 3.37 ;
     Kval_ldstrip = 2.58 ;
     Cval_ldstrip = 3.37 ;
     Kval_all = 2.19 ;
     Cval_all = 3.10 ;

     Cval_nol1 = 3.19 ; // fit C in 5-20
     Kval_nol1 = 2.28 ;
     Cval_nol1_2 = 3.17 ; // fit C in 3-5
     Kval_nol1_2 = 2.30 ;
     Cval_nol1_3 = 2.80 ;
     Kval_nol1_3 = 2.83 ;
     Nval_nol1_3 = 0.090 ;

     Kval_strip = 2.50 ;
     Cval_strip = 3.19 ;
     Kval_hdnol1 = 2.23 ;
     Cval_hdnol1 = 3.03 ;

/* not possible to extract K and C values on the pixel only samples  ==> do not update and keep the MC values
     Kval_pix = 2.173 ;
     Cval_pix = 3.242 ;
     Kval_pixnol1 = 2.173 ;
     Cval_pixnol1 = 3.242 ;
*/
    }
    else if (year==2018) {
     //extracted on analysis_ul_2018_28feb.root
     Kval_ld = 2.35 ;
     Cval_ld = 3.36 ;
     Kval_ldstrip = 2.60 ;
     Cval_ldstrip = 3.38 ;
     Kval_all = 2.17 ;
     Cval_all = 3.09 ;

     Cval_nol1 = 3.18 ; // fit C in 5-20
     Kval_nol1 = 2.25 ;
     Cval_nol1_2 = 3.16 ; // fit C in 3-5
     Kval_nol1_2 = 2.27 ;
     Cval_nol1_3 = 2.90 ;
     Kval_nol1_3 = 2.52 ;
     Nval_nol1_3 = 0.066 ;

     Kval_strip = 2.50 ;
     Cval_strip = 3.19 ;
     Kval_hdnol1 = 2.21 ;
     Cval_hdnol1 = 3.02 ;

// not possible to extract K and C values on the pixel only samples  ==> do not update and keep the MC values
     Kval_pix = 1.72 ;
     Cval_pix = 2.96 ;
// MC 2018 not available 
/*
     Kval_pixnol1 = .. ;
     Cval_pixnol1 = .. ;
*/
    }
   }
   else if (year==2018) {
     // MC : extracted on analysis_ul_2018MC_w18_MC_28feb.root
     Kval_ld = 2.34 ;
     Cval_ld = 3.35 ;
     Kval_ldstrip = 2.48 ;
     Cval_ldstrip = 3.38 ;
     Kval_all = 2.20 ;
     Cval_all = 3.14 ;

     Cval_nol1 = 3.29 ; // fit C in 5-20
     Kval_nol1 = 2.22 ;
     Cval_nol1_2 = 3.16 ; // fit C in 3-5
     Kval_nol1_2 = 2.27 ;
     Cval_nol1_3 = 2.86 ;
     Kval_nol1_3 = 2.61 ;
     Nval_nol1_3 = 0.103 ;

     Kval_strip = 2.40 ;
     Cval_strip = 3.21 ;
     Kval_hdnol1 = 2.23 ;
     Cval_hdnol1 = 3.05 ;

// MC values : no convergence... oscillation between 1.64 and 1.72 for Kval_pix :(
     Kval_pix = 1.72 ;
     Cval_pix = 2.96 ;
// MC values : not possible to fit --> K =1 at limit :(
/*
     Kval_pixnol1 = .. ;
     Cval_pixnol1 = .. ;
*/
   }
   bool blind_data= false;
   if (dataFlag)  blind_data=true;


// histograms

   TH1D* HNtracks = new TH1D("HNtracks", "HNtracks", 50, -0.5,49.5);
   TH1D* HNtracks1 = new TH1D("HNtracks1", "Ntracks with pT>1", 40, -0.5,39.5);
   TH1D* HNtracks20 = new TH1D("HNtracks20", "Ntracks with pT>20", 30, -0.5,29.5);

   TH1D* Htrackpt = new TH1D("Htrackpt", "track pT", 500, 0.,500);
   TH1D* Htracketa = new TH1D("Htracketa", "track eta", 24, -3.,3.);
   TH1D* Htracketa_lowp = new TH1D("Htracketa_lowp", "track eta", 24, -3.,3.);
   TH1D* Htrackphi = new TH1D("Htrackphi", "track phi", 24, -1.*pi,pi);
   TH1D* Htracknhit = new TH1D("Htracknhit", "track nhit", 50, -0.5,49.5);

   TH1D* Htrackih_reco = new TH1D("Htrackih_reco", "track ih (reco)", 50, 0.,20);
   TH1D* Htrackih_pix = new TH1D("Htrackih_pix", "track ih in pix", 50, 0.,20);
   TH1D* Htrackih_strip = new TH1D("Htrackih_strip", "track ih in strip", 50, 0.,20);
   TH1D* Htrackdedx_pix = new TH1D("Htrackdedx_pix", "track dedx in pix", 100, 0.,500000);
   TH1D* Htrackdedx_strip = new TH1D("Htrackdedx_strip", "track dedx in strip", 100, 0.,2000);
   TH1D* Htrackias = new TH1D("Htrackias", "track ias (noL1)", 80, 0.,1);
   TH1D* Htrackiasall = new TH1D("Htrackiasall", "track ias (all)", 80, 0.,1);

   TH1D* Htrackih_lowp = new TH1D("Htrackih_lowp", "track ih", 50, 0.,20);
   TH1D* Htrackih_pix_lowp = new TH1D("Htrackih_pix_lowp", "track ih in pix", 50, 0.,20);
   TH1D* Htrackih_strip_lowp = new TH1D("Htrackih_strip_lowp", "track ih in strip", 50, 0.,20);
   TH1D* Htrackih0_lowp = new TH1D("Htrackih0_lowp", "track ih0 ", 50, 0.,20);
   TH1D* Htrackih0noL1_lowp = new TH1D("Htrackih0noL1_lowp", "track ih0 noL1", 50, 0.,20);
   TH1D* Htrackias_lowp = new TH1D("Htrackias_lowp", "track ias (noL1)", 80, 0.,1);
   TH1D* Htrackiasall_lowp = new TH1D("Htrackiasall_lowp", "track ias (all)", 80, 0.,1);

   TH1D* Htrackdedx_pix_lowp = new TH1D("Htrackdedx_pix_lowp", "track dedx in pix", 100, 0.,500000);
   TH1D* Htrackdedx_strip_lowp = new TH1D("Htrackdedx_strip_lowp", "track dedx in strip", 100, 0.,2000);
   TH1D* Htrackdedx_strip_lowp1 = new TH1D("Htrackdedx_strip_lowp1", "track dedx in strip", 100, 0.,20);
   TH1D* Htrackdedx_strip_lowp2 = new TH1D("Htrackdedx_strip_lowp2", "track dedx in strip", 100, 0.,20);

   TH1D* Nsat = new TH1D("Nsat", "Nsat", 20, -0.5,19.5);
   TH1D* NPix = new TH1D("NPix", "NPix", 10, -0.5,9.5);
   TH1D* NStrip = new TH1D("NStrip", "NStrip", 30, -0.5,29.5);
   TH2D* dEdXVsP = new TH2D("dEdXVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXVsP_lowp = new TH2D("dEdXVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXVsP_lowp2 = new TH2D("dEdXVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdXpixVsP = new TH2D("dEdXpixVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXpixVsP_lowp = new TH2D("dEdXpixVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXpixVsP_lowp2 = new TH2D("dEdXpixVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdXstripVsP = new TH2D("dEdXstripVsP", "dEdX:P", 20,0,5000, 40, 0.,20);
   TH2D* dEdXstripVsP_lowp = new TH2D("dEdXstripVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXstripVsP_lowp2 = new TH2D("dEdXstripVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0VsP_lowp = new TH2D("dEdX0VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0VsP_lowp2 = new TH2D("dEdX0VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0stripVsP_lowp = new TH2D("dEdX0stripVsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0stripVsP_lowp2 = new TH2D("dEdX0stripVsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_lowp = new TH2D("dEdX0noL1VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_lowp2 = new TH2D("dEdX0noL1VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta1_lowp = new TH2D("dEdX0noL1VsP_eta1_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta1_lowp2 = new TH2D("dEdX0noL1VsP_eta1_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta2_lowp = new TH2D("dEdX0noL1VsP_eta2_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta2_lowp2 = new TH2D("dEdX0noL1VsP_eta2_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta3_lowp = new TH2D("dEdX0noL1VsP_eta3_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_eta3_lowp2 = new TH2D("dEdX0noL1VsP_eta3_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu1_lowp = new TH2D("dEdX0noL1VsP_pu1_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu1_lowp2 = new TH2D("dEdX0noL1VsP_pu1_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu2_lowp = new TH2D("dEdX0noL1VsP_pu2_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu2_lowp2 = new TH2D("dEdX0noL1VsP_pu2_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu3_lowp = new TH2D("dEdX0noL1VsP_pu3_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu3_lowp2 = new TH2D("dEdX0noL1VsP_pu3_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu4_lowp = new TH2D("dEdX0noL1VsP_pu4_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0noL1VsP_pu4_lowp2 = new TH2D("dEdX0noL1VsP_pu4_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdXHDnoL1VsP_lowp = new TH2D("dEdXHDnoL1VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdXHDnoL1VsP_lowp2 = new TH2D("dEdXHDnoL1VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEdX0pixnoL1VsP_lowp = new TH2D("dEdX0pixnoL1VsP_lowp", "dEdX:P", 50,0,5, 80, 2.,10.);
   TH2D* dEdX0pixnoL1VsP_lowp2 = new TH2D("dEdX0pixnoL1VsP_lowp2", "dEdX:P", 125,0,25, 80, 2.,10.);
   TH2D* dEstrVsdE_lowp = new TH2D("dEstrVsdE_lowp", "dEdX(all):dEdX(strip)", 50,2.,7., 80, 2.,10.);
   TH2D* dEdXstripVsEta_lowp = new TH2D("dEdXstripVsEta_lowp", "dEdX:P", 12, -3.,3., 80, 2.,10.);
   TH2D* EtaVsPhi_nhit = new TH2D("EtaVsPhi_nhit", "Eta:Phi", 12, -1.*pi,pi,12, -3.,3.);


   TH2D* dEdXstripVsNhit_lowp  = new TH2D("dEdXstripVsNhit_lowp", "dEdX(strip):Nhit", 20, -0.5,19.5,80, 2.,10.);
   TH2D* dEdXstripVsNhittrunc_lowp  = new TH2D("dEdXstripVsNhittrunc_lowp", "dEdX(strip):Nhit", 20, -0.5,19.5,80, 2.,10.);
   TH2D* dEdXstripVsCharge_lowp = new TH2D("dEdXstripVsCharge_lowp", "dEdX(strip):Charge", 100, 0.,20, 80, 2.,10.);

   TH1D* Charge_pixl1 = new TH1D("Charge_pixl1","Charge per layer",400, 0.,20);
   TH1D* Charge_pixl2 = new TH1D("Charge_pixl2","Charge per layer",400, 0.,20);
   TH1D* Charge_pixl3 = new TH1D("Charge_pixl3","Charge per layer",400, 0.,20);
   TH1D* Charge_pixl4 = new TH1D("Charge_pixl4","Charge per layer",400, 0.,20);
   TH1D* Charge_pixd1 = new TH1D("Charge_pixd1","Charge per layer",400, 0.,20);
   TH1D* Charge_pixd2 = new TH1D("Charge_pixd2","Charge per layer",400, 0.,20);
   TH1D* Charge_pixd3 = new TH1D("Charge_pixd3","Charge per layer",400, 0.,20);
   TH1D* Charge_pixr1 = new TH1D("Charge_pixr1","Charge per layer",400, 0.,20);
   TH1D* Charge_pixr2 = new TH1D("Charge_pixr2","Charge per layer",400, 0.,20);
   TH1D* Charge_tib1 = new TH1D("Charge_tib1","Charge per layer",400, 0.,20);
   TH1D* Charge_tib2 = new TH1D("Charge_tib2","Charge per layer",400, 0.,20);
   TH1D* Charge_tib3 = new TH1D("Charge_tib3","Charge per layer",400, 0.,20);
   TH1D* Charge_tib4 = new TH1D("Charge_tib4","Charge per layer",400, 0.,20);
   TH1D* Charge_tob1 = new TH1D("Charge_tob1","Charge per layer",400, 0.,20);
   TH1D* Charge_tob2 = new TH1D("Charge_tob2","Charge per layer",400, 0.,20);
   TH1D* Charge_tob3 = new TH1D("Charge_tob3","Charge per layer",400, 0.,20);
   TH1D* Charge_tob4 = new TH1D("Charge_tob4","Charge per layer",400, 0.,20);
   TH1D* Charge_tob5 = new TH1D("Charge_tob5","Charge per layer",400, 0.,20);
   TH1D* Charge_tob6 = new TH1D("Charge_tob6","Charge per layer",400, 0.,20);
   TH1D* Charge_tid1 = new TH1D("Charge_tid1","Charge per layer",400, 0.,20);
   TH1D* Charge_tid2 = new TH1D("Charge_tid2","Charge per layer",400, 0.,20);
   TH1D* Charge_tid3 = new TH1D("Charge_tid3","Charge per layer",400, 0.,20);
   TH1D* Charge_tec1 = new TH1D("Charge_tec1","Charge per layer",400, 0.,20);
   TH1D* Charge_tec2 = new TH1D("Charge_tec2","Charge per layer",400, 0.,20);
   TH1D* Charge_tec3 = new TH1D("Charge_tec3","Charge per layer",400, 0.,20);
   TH1D* Charge_tec4 = new TH1D("Charge_tec4","Charge per layer",400, 0.,20);
   TH1D* Charge_tec5 = new TH1D("Charge_tec5","Charge per layer",400, 0.,20);
   TH1D* Charge_tec6 = new TH1D("Charge_tec6","Charge per layer",400, 0.,20);
   TH1D* Charge_tec7 = new TH1D("Charge_tec7","Charge per layer",400, 0.,20);
   TH1D* Charge_tec8 = new TH1D("Charge_tec8","Charge per layer",400, 0.,20);
   TH1D* Charge_tec9 = new TH1D("Charge_tec9","Charge per layer",400, 0.,20);

   TH1D* LowpCharge_tib1 = new TH1D("LowpCharge_tib1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tib2 = new TH1D("LowpCharge_tib2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tib3 = new TH1D("LowpCharge_tib3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tib4 = new TH1D("LowpCharge_tib4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob1 = new TH1D("LowpCharge_tob1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob2 = new TH1D("LowpCharge_tob2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob3 = new TH1D("LowpCharge_tob3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob4 = new TH1D("LowpCharge_tob4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob5 = new TH1D("LowpCharge_tob5","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tob6 = new TH1D("LowpCharge_tob6","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tid1 = new TH1D("LowpCharge_tid1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tid2 = new TH1D("LowpCharge_tid2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tid3 = new TH1D("LowpCharge_tid3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec1 = new TH1D("LowpCharge_tec1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec2 = new TH1D("LowpCharge_tec2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec3 = new TH1D("LowpCharge_tec3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec4 = new TH1D("LowpCharge_tec4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec5 = new TH1D("LowpCharge_tec5","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec6 = new TH1D("LowpCharge_tec6","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec7 = new TH1D("LowpCharge_tec7","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec8 = new TH1D("LowpCharge_tec8","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_tec9 = new TH1D("LowpCharge_tec9","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl1 = new TH1D("LowpCharge_pixl1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl2 = new TH1D("LowpCharge_pixl2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl3 = new TH1D("LowpCharge_pixl3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixl4 = new TH1D("LowpCharge_pixl4","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixd1 = new TH1D("LowpCharge_pixd1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixd2 = new TH1D("LowpCharge_pixd2","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixd3 = new TH1D("LowpCharge_pixd3","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixr1 = new TH1D("LowpCharge_pixr1","Charge per layer",400, 0.,20);
   TH1D* LowpCharge_pixr2 = new TH1D("LowpCharge_pixr2","Charge per layer",400, 0.,20);
   TH2D* LowpCharge_Eta_pix = new TH2D("LowpCharge_Eta_pix", "Charge/pathlength:eta", 12, -3.,3., 200, 0.,10);
   TH2D* LowpCharge_Eta_strip = new TH2D("LowpCharge_Eta_strip", "Charge/pathlength:eta", 12, -3.,3., 200, 0.,10);

   TH2D* ChargeVsRun_pixl1 = new TH2D("ChargeVsRun_pixl1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixl2 = new TH2D("ChargeVsRun_pixl2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixl3 = new TH2D("ChargeVsRun_pixl3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixl4 = new TH2D("ChargeVsRun_pixl4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixd1 = new TH2D("ChargeVsRun_pixd1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixd2 = new TH2D("ChargeVsRun_pixd2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixd3 = new TH2D("ChargeVsRun_pixd3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixr1 = new TH2D("ChargeVsRun_pixr1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_pixr2 = new TH2D("ChargeVsRun_pixr2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib1 = new TH2D("ChargeVsRun_tib1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib2 = new TH2D("ChargeVsRun_tib2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib3 = new TH2D("ChargeVsRun_tib3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tib4 = new TH2D("ChargeVsRun_tib4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob1 = new TH2D("ChargeVsRun_tob1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob2 = new TH2D("ChargeVsRun_tob2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob3 = new TH2D("ChargeVsRun_tob3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob4 = new TH2D("ChargeVsRun_tob4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob5 = new TH2D("ChargeVsRun_tob5","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tob6 = new TH2D("ChargeVsRun_tob6","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tid1 = new TH2D("ChargeVsRun_tid1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tid2 = new TH2D("ChargeVsRun_tid2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tid3 = new TH2D("ChargeVsRun_tid3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec1 = new TH2D("ChargeVsRun_tec1","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec2 = new TH2D("ChargeVsRun_tec2","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec3 = new TH2D("ChargeVsRun_tec3","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec4 = new TH2D("ChargeVsRun_tec4","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec5 = new TH2D("ChargeVsRun_tec5","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec6 = new TH2D("ChargeVsRun_tec6","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec7 = new TH2D("ChargeVsRun_tec7","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec8 = new TH2D("ChargeVsRun_tec8","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);
   TH2D* ChargeVsRun_tec9 = new TH2D("ChargeVsRun_tec9","ChargeVsRun per layer",545, 271000,325500,400, 0.,20);

   TH2D* ZooChargeVsRun_pixl1 = new TH2D("ZooChargeVsRun_pixl1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl2 = new TH2D("ZooChargeVsRun_pixl2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl3 = new TH2D("ZooChargeVsRun_pixl3","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixl4 = new TH2D("ZooChargeVsRun_pixl4","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd1 = new TH2D("ZooChargeVsRun_pixd1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd2 = new TH2D("ZooChargeVsRun_pixd2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixd3 = new TH2D("ZooChargeVsRun_pixd3","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixr1 = new TH2D("ZooChargeVsRun_pixr1","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
   TH2D* ZooChargeVsRun_pixr2 = new TH2D("ZooChargeVsRun_pixr2","ChargeVsRun per layer",545, 271000,325500,90, 0., 4.5);
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
   TH1D* HSCP_MassIh = new TH1D("HSCP_MassIh", "Mass via Ih", 200, 0.,4000.);
   TH1D* HSCP_MassIh0 = new TH1D("HSCP_MassIh0", "Mass via Ih(no drop)", 200, 0.,4000.);
   TH1D* HSCP_MassIhstrip = new TH1D("HSCP_MassIhstrip", "Mass via Ih(strip)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1 = new TH1D("HSCP_MassIh0noL1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2 = new TH1D("HSCP_MassIh0noL1_2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3 = new TH1D("HSCP_MassIh0noL1_3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_11 = new TH1D("HSCP_MassIh0noL1_11", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_12 = new TH1D("HSCP_MassIh0noL1_12", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_13 = new TH1D("HSCP_MassIh0noL1_13", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_1s1 = new TH1D("HSCP_MassIh0noL1_1s1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_1s2 = new TH1D("HSCP_MassIh0noL1_1s2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_1s3 = new TH1D("HSCP_MassIh0noL1_1s3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2s1 = new TH1D("HSCP_MassIh0noL1_2s1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2s2 = new TH1D("HSCP_MassIh0noL1_2s2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_2s3 = new TH1D("HSCP_MassIh0noL1_2s3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3s1 = new TH1D("HSCP_MassIh0noL1_3s1", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3s2 = new TH1D("HSCP_MassIh0noL1_3s2", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0noL1_3s3 = new TH1D("HSCP_MassIh0noL1_3s3", "Mass via Ih(no drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassIh0strip = new TH1D("HSCP_MassIh0strip", "Mass via Ih(no drop strip)", 200, 0.,4000.);
   TH1D* HSCP_MassIhHDnoL1 = new TH1D("HSCP_MassIhHDnoL1", "Mass via Ih(High drop noL1)", 200, 0.,4000.);
   TH1D* HSCP_MassTOF = new TH1D("HSCP_MassTOF", "Mass via 1/beta", 200, 0.,4000.);
   TH2D* HSCP2d_MassTOFvsIh = new TH2D("HSCP2d_MassTOFvsIh", "Mass Ih vs 1/beta",50, 0.,4000., 50, 0.,4000.);
   TH2D* HSCP2d_MassIh = new TH2D("HSCP2d_MassIh", "Mass via Ih",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIh0 = new TH2D("HSCP2d_MassIh0", "Mass via Ih(no drop)",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIhstrip = new TH2D("HSCP2d_MassIhstrip", "Mass via Ih(strip)",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIh0noL1 = new TH2D("HSCP2d_MassIh0noL1", "Mass via Ih(no drop noL1)",30,0.,3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIh0strip = new TH2D("HSCP2d_MassIh0strip", "Mass via Ih(no drop strip)",30, 0, 3000, 100, 0.,4000.);
   TH2D* HSCP2d_MassIhHDnoL1 = new TH2D("HSCP2d_MassIhHDnoL1", "Mass via Ih(High drop noL1)",30,0,3000, 100, 0.,4000.);

   TH2D* HSCP2d_Mass_pix_strip15 = new TH2D("HSCP2d_Mass_pix_strip15", "Mass Pixel vs Mass Strip (15p drop)",100,0,4000, 100, 0.,4000.);
   TH2D* HSCP2d_Mass_pix_strip0 = new TH2D("HSCP2d_Mass_pix_strip0", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,4000, 100, 0.,4000.);
   TH2D* HSCP2d_Mass_pix_strip = new TH2D("HSCP2d_Mass_pix_strip", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,500, 100, 0.,500.);
   TH1D* HSCP_MassDiff_pix_strip0 = new TH1D("HSCP_MassDiff_pix_strip0", "MassDiff Pixel-Strip (no drop noL1)",200,-4000,4000);
   TH1D* HSCP_MassDiff_pix_strip15 = new TH1D("HSCP_MassDiff_pix_strip15", "MassDiff Pixel-Strip (15p drop)",200,-4000,4000);
   TH1D* HSCP_MassResol_pix_strip0 = new TH1D("HSCP_MassResol_pix_strip0", "MassResol Pixel-Strip (no drop noL1)",200,-2,2);
   TH1D* HSCP_MassResol_pix_strip15 = new TH1D("HSCP_MassResol_pix_strip15", "MassResol Pixel-Strip (15p drop)",200,-2,2);

   TH1D* lowp_MassIh = new TH1D("lowp_MassIh", "Mass via Ih", 100, 0.,5.);
   TH1D* lowp_MassIh0 = new TH1D("lowp_MassIh0", "Mass via Ih(no drop)", 100, 0.,5.);
   TH1D* lowp_MassIhstrip = new TH1D("lowp_MassIhstrip", "Mass via Ih(strip)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1 = new TH1D("lowp_MassIh0noL1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2 = new TH1D("lowp_MassIh0noL1_2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3 = new TH1D("lowp_MassIh0noL1_3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_11 = new TH1D("lowp_MassIh0noL1_11", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_12 = new TH1D("lowp_MassIh0noL1_12", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_13 = new TH1D("lowp_MassIh0noL1_13", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_1s1 = new TH1D("lowp_MassIh0noL1_1s1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_1s2 = new TH1D("lowp_MassIh0noL1_1s2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_1s3 = new TH1D("lowp_MassIh0noL1_1s3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2s1 = new TH1D("lowp_MassIh0noL1_2s1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2s2 = new TH1D("lowp_MassIh0noL1_2s2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_2s3 = new TH1D("lowp_MassIh0noL1_2s3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3s1 = new TH1D("lowp_MassIh0noL1_3s1", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3s2 = new TH1D("lowp_MassIh0noL1_3s2", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0noL1_3s3 = new TH1D("lowp_MassIh0noL1_3s3", "Mass via Ih(no drop noL1)", 100, 0.,5.);
   TH1D* lowp_MassIh0strip = new TH1D("lowp_MassIh0strip", "Mass via Ih(no drop strip)", 100, 0.,5.);
   TH1D* lowp_MassIhHDnoL1 = new TH1D("lowp_MassIhHDnoL1", "Mass via Ih(High drop noL1)", 100, 0.,5.);
   TH2D* lowp2d_MassIh = new TH2D("lowp2d_MassIh", "Mass via Ih", 30,0,3, 100, 0.,2.);
   TH2D* lowp2d_MassIh0 = new TH2D("lowp2d_MassIh0", "Mass via Ih(no drop)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIhstrip = new TH2D("lowp2d_MassIhstrip", "Mass via Ih(strip)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIh0noL1 = new TH2D("lowp2d_MassIh0noL1", "Mass via Ih(no drop noL1)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIh0strip = new TH2D("lowp2d_MassIh0strip", "Mass via Ih(no drop strip)",30,0,3, 100, 0.,5.);
   TH2D* lowp2d_MassIhHDnoL1 = new TH2D("lowp2d_MassIhHDnoL1", "Mass via Ih(High drop noL1)",30,0,3, 100, 0.,5.);
   TH2D* lowp_dEdXpixVsstrip = new TH2D("lowp_dEdXpixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* lowp_dEdX0pixVsstrip = new TH2D("lowp_dEdX0pixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* lowp2d_Mass_pix_strip15 = new TH2D("lowp2d_Mass_pix_strip15", "Mass Pixel vs Mass Strip (15p drop)",100,0,5, 100, 0.,5.);
   TH2D* lowp2d_Mass_pix_strip0 = new TH2D("lowp2d_Mass_pix_strip0", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,5, 100, 0.,5.);
   TH2D* bg_lowp2d_Mass_pix_strip0 = new TH2D("bg_lowp2d_Mass_pix_strip0", "Mass Pixel vs Mass Strip (no drop noL1)",100,0,100, 100, 0.,100.);
   TH2D* bg_transf_Mass = new TH2D("bg_transf_Mass", "Diff Mass Pixel-Strip  vs Mean (no drop noL1)",100,0,100, 100, -2.,2.);
   TH1D* lowp_MassDiff_pix_strip0 = new TH1D("lowp_MassDiff_pix_strip0", "MassDiff Pixel-Strip (no drop noL1)",200,-5,5);
   TH1D* lowp_MassDiff_pix_strip15 = new TH1D("lowp_MassDiff_pix_strip15", "MassDiff Pixel-Strip (15p drop)",200,-5,5);
   TH2D* bg_dEdXVsIas = new TH2D("bg_dEdXVsIas","Ih:Ias ", 100, 0.,1., 100, 0.,10.);

   TH1D* HSCP_dEdXHiDrop  = new TH1D("HSCP_dEdXHiDrop", "dEdX HighDrop", 60, 0.,15.);
   TH1D* HSCP_dEdXstripHiDrop  = new TH1D("HSCP_dEdXstripHiDrop", "dEdX HighDrop(strip)", 60, 0.,15.);
   TH1D* HSCP_dEdXpixHiDrop = new TH1D("HSCP_dEdXpixHiDrop", "dEdX HighDrop(pix)", 60, 0.,15.);
   TH1D* HSCP_dEdXHiDropNoL1  = new TH1D("HSCP_dEdXHiDropNoL1", "dEdX HighDrop No L1", 60, 0.,15.);
   TH1D* HSCP_dEdX0NoL1  = new TH1D("HSCP_dEdX0NoL1", "dEdX No L1", 60, 0.,15.);

   TH1D* HSCP_FMIP4 = new TH1D("HSCP_FMIP4", "FMIP(4)", 50, 0.,1.);
   TH1D* HSCP_FMIP3p5 = new TH1D("HSCP_FMIP3p5", "FMIP(3.5)", 50, 0.,1.);
   TH1D* HSCP_FMIP3p2 = new TH1D("HSCP_FMIP3p2", "FMIP(3.2)", 50, 0.,1.);

   TH1D* HSCP_iasnol1 = new TH1D("HSCP_iasnol1", "Ias (noL1)", 80, 0.,1.);
   TH1D* HSCP_iasall = new TH1D("HSCP_iasall", "Ias (all)", 80, 0.,1.);
   TH1D* HSCP_probQ = new TH1D("HSCP_probQ", "ProbQ", 80, 0.,1.);
   TH1D* HSCP_probQNoL1 = new TH1D("HSCP_probQNoL1", "ProbQNoL1", 80, 0.,1.);
   TH1D* HSCP_probXY = new TH1D("HSCP_probXY", "ProbXY", 80, 0.,1.);
   TH1D* HSCP_probXYNoL1 = new TH1D("HSCP_probXYNoL1", "ProbXYNoL1", 80, 0.,1.);

   TH1D* HSCP_pt = new TH1D("HSCP_pt", "pT",  50, 55.,1550);
   TH1D* HSCP_eta = new TH1D("HSCP_eta", "eta",  24, -3.,3.);
   TH1D* HSCP_iso_eop = new TH1D("HSCP_iso_eop", "Isolation (ECAL+HCAL)/p",  50, 0., 5.);
   TH1D* nPV = new TH1D("nPV", "nPV",  40,0,80);
   TH1D* HSCP_invB = new TH1D("HSCP_invB", "invBeta",  300,-1,2);
   TH1D* HSCP_errinvB = new TH1D("HSCP_errinvB", "err_invBeta",  50,0,0.5);
   TH1D* HSCP_invBDT = new TH1D("HSCP_invBDT", "invBeta(DT)",  90,-1,2);
   TH1D* HSCP_invBCSC = new TH1D("HSCP_invBCSC", "invBeta(CSC)",  90,-1,2);
   TH1D* HSCP_time = new TH1D("HSCP_time", "VertexTiming",  200,-100,100);
   TH1D* HSCP_npix = new TH1D("HSCP_npix", "#(pix)", 15, -0.5,14.5);
   TH1D* HSCP_nstrip = new TH1D("HSCP_nstrip", "#(strip)", 40, -0.5,39.5);
   TH1D* HSCP_nmpix = new TH1D("HSCP_nmpix", "#m(pix)", 15, -0.5,14.5);
   TH1D* HSCP_nmstrip = new TH1D("HSCP_nmstrip", "#m(strip)", 40, -0.5,39.5);
   TH1D* HSCP_nratio = new TH1D("HSCP_nratio", "#(ratio)", 100, 0.,1.);
   TH1D* HSCP_nmratio = new TH1D("HSCP_nmratio", "#m(ratio)", 100, 0.,1.);

   TH2D* HSCP_dEdXpixVsstrip = new TH2D("HSCP_dEdXpixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* HSCP_dEdX0pixVsstrip = new TH2D("HSCP_dEdX0pixVsstrip","Pix dE/dX:Strip dE/dX", 100, 0.,10., 100, 0.,10.);
   TH2D* HSCP_dEdXstripVsall  = new TH2D("HSCP_dEdXstripVsall","Strip dE/dX: Ih(all)", 100, 0.,10., 100, 0.,10.);
   TH2D* HSCP_dEdXpixVsall   = new TH2D("HSCP_dEdXpixVsall","Pix dE/dX: Ih(all)", 100, 0.,10., 100, 0.,10.);

   TH2D* dEdXVsRun = new TH2D("dEdXVsRun", "dEdX:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXpixVsRun = new TH2D("dEdXpixVsRun", "dEdX(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXstripVsRun = new TH2D("dEdXstripVsRun", "dEdX(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXNoL1VsRun = new TH2D("dEdXNoL1VsRun", "dEdX(noL1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXNoL1pixVsRun = new TH2D("dEdXNoL1pixVsRun", "dEdX(Pix only with no L1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXHiDropVsRun = new TH2D("dEdXHiDropVsRun","dEdX HighDrop:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXpixHiDropVsRun = new TH2D("dEdXpixHiDropVsRun","dEdX HighDrop(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXstripHiDropVsRun = new TH2D("dEdXstripHiDropVsRun","dEdX HighDrop(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdXHiDropNoL1VsRun = new TH2D("dEdXHiDropNoL1VsRun","dEdX HighDrop:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0VsRun = new TH2D("dEdX0VsRun", "dEdX0:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0pixVsRun = new TH2D("dEdX0pixVsRun", "dEdX0(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0stripVsRun = new TH2D("dEdX0stripVsRun", "dEdX0(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* MassStripVsRun = new TH2D("MassStripVsRun", "Mass:Run", 545, 271000,325500, 80, 0.,4000.);
   TH2D* MassNoL1VsRun = new TH2D("MassNoL1VsRun", "Mass:Run", 545, 271000,325500, 80, 0.,4000.);
   TH2D* dEdX0NoL1pixVsRun = new TH2D("dEdX0NoL1pixVsRun", "dEdX(Pix only with no L1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX0NoL1VsRun = new TH2D("dEdX0NoL1VsRun", "dEdX(noL1Pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4VsRun = new TH2D("dEdXV4sRun", "dEdX:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4pixVsRun = new TH2D("dEdX4pixVsRun", "dEdX(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX4stripVsRun = new TH2D("dEdX4stripVsRun", "dEdX(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40VsRun = new TH2D("dEdX40VsRun", "dEdX0:Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40pixVsRun = new TH2D("dEdX40pixVsRun", "dEdX0(pix):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* dEdX40stripVsRun = new TH2D("dEdX40stripVsRun", "dEdX0(strip):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* bg_dEdX0NoL1VsRun = new TH2D("bg_dEdX0NoL1VsRun", "dEdX(noL1):Run", 545, 271000,325500, 60, 0.,15.);
   TH2D* iasNoL1VsRun = new TH2D("iasNoL1VsRun", "Ias(noL1):Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* iasAllVsRun = new TH2D("iasAllVsRun", "Ias(all):Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probQVsRun = new TH2D("probQVsRun", "ProbQ:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probQNoL1VsRun = new TH2D("probQNoL1VsRun", "ProbQ:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probXYVsRun = new TH2D("probXYVsRun", "ProbXY:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probXYNoL1VsRun = new TH2D("probXYNoL1VsRun", "ProbXY:Run", 545, 271000,325500, 80, 0.,1.);
   TH2D* probQVsIas = new TH2D("probQVsIas", "ProbQ:Ias", 80, 0.,1., 80, 0.,1.);

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


   TH1D* HHitPix = new TH1D( "HHitPix", "HHitPix", 200, 0, 20);
   TProfile* HHitProfilePix = new TProfile( "HHitProfilePix", "HHitProfilePix", 50, 0, 100);
   TH2D* HHit2DPix = new TH2D( "HHit2DPix", "HHit2DPix", 50, 0, 100,200, 0, 20);
   TH2D* HHit2DPix_NoM = new TH2D( "HHit2DPix_NoM", "HHit2DPix", 50, 0, 100,200, 0, 20);
   TH1D* HHitStrip = new TH1D( "HHitStrip", "HHitStrip", 200, 0, 20);
   TProfile* HHitProfileStrip = new TProfile( "HHitProfileStrip", "HHitProfileStrip", 50, 0, 100);
   TH2D* HHit2DStrip = new TH2D( "HHit2DStrip", "HHit2DStrip", 50, 0, 100,200, 0, 20);
   TH2D* HHit2DStrip_NoM = new TH2D( "HHit2DStrip_NoM", "HHit2DStrip", 50, 0, 100,200, 0, 20);

   const double P_Min               = 1   ;
   const double P_Max               = 16  ; // 1 + 14 + 1; final one is for pixel!
   const int    P_NBins             = 15  ; // 15th bin = pixel; 0 is underflow
   const double Path_Min            = 0.2 ;
   const double Path_Max            = 1.6 ;
   const int    Path_NBins          = 42  ;
   const double Charge_Min          = 0   ;
   const double Charge_Max          = 5000;
   const int Charge_NBins = 500 ;
   TH3D* Charge_Vs_Path = new TH3D( "Charge_Vs_Path", "Charge_Vs_Path", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH3D* Charge_Vs_Path_noL1 = new TH3D( "Charge_Vs_Path_noL1", "Charge_Vs_Path_noL1", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH3D* Charge_Vs_Path_NoM = new TH3D( "Charge_Vs_Path_NoM", "Charge_Vs_Path", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH3D* Charge_Vs_Path_noL1_NoM = new TH3D( "Charge_Vs_Path_noL1_NoM", "Charge_Vs_Path_noL1", P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);


   HHitPix->Sumw2();
   HHitProfilePix->Sumw2();
   HHit2DPix->Sumw2();
   HHit2DPix_NoM->Sumw2();
   HHitStrip->Sumw2();
   HHitProfileStrip->Sumw2();
   HHit2DStrip->Sumw2();
   HHit2DStrip_NoM->Sumw2();
   Charge_Vs_Path->Sumw2();
   Charge_Vs_Path_noL1->Sumw2();
   Charge_Vs_Path_NoM->Sumw2();
   Charge_Vs_Path_noL1_NoM->Sumw2();

   HNtracks->Sumw2();
   HNtracks1->Sumw2();
   HNtracks20->Sumw2();

   Htrackpt->Sumw2();
   Htracketa->Sumw2();
   Htracketa_lowp->Sumw2();
   Htrackphi->Sumw2();
   Htracknhit->Sumw2();

   Htrackih_reco->Sumw2();
   Htrackih_pix->Sumw2();
   Htrackih_strip->Sumw2();
   Htrackdedx_pix->Sumw2();
   Htrackdedx_strip->Sumw2();
   Htrackias->Sumw2();
   Htrackiasall->Sumw2();

   Htrackih_lowp->Sumw2();
   Htrackih0noL1_lowp->Sumw2();
   Htrackih0_lowp->Sumw2();
   Htrackih_pix_lowp->Sumw2();
   Htrackih_strip_lowp->Sumw2();
   Htrackdedx_pix_lowp->Sumw2();
   Htrackdedx_strip_lowp->Sumw2();
   Htrackdedx_strip_lowp1->Sumw2();
   Htrackdedx_strip_lowp2->Sumw2();
   Htrackias_lowp->Sumw2();
   Htrackiasall_lowp->Sumw2();


   Nsat->Sumw2();
   NPix->Sumw2();
   NStrip->Sumw2();
   dEdXVsP->Sumw2();
   dEdXVsP_lowp->Sumw2();
   dEdXVsP_lowp2->Sumw2();
   dEdXpixVsP->Sumw2();
   dEdXpixVsP_lowp->Sumw2();
   dEdXpixVsP_lowp2->Sumw2();
   dEdXstripVsP->Sumw2();
   dEdXstripVsP_lowp->Sumw2();
   dEdXstripVsP_lowp2->Sumw2();
   dEdX0VsP_lowp->Sumw2();
   dEdX0VsP_lowp2->Sumw2();
   dEdX0stripVsP_lowp->Sumw2();
   dEdX0stripVsP_lowp2->Sumw2();
   dEdX0noL1VsP_lowp->Sumw2();
   dEdX0noL1VsP_lowp2->Sumw2();
   dEdX0noL1VsP_eta1_lowp->Sumw2();
   dEdX0noL1VsP_eta1_lowp2->Sumw2();
   dEdX0noL1VsP_eta2_lowp->Sumw2();
   dEdX0noL1VsP_eta2_lowp->Sumw2();
   dEdX0noL1VsP_eta3_lowp->Sumw2();
   dEdX0noL1VsP_eta3_lowp2->Sumw2();
   dEdX0noL1VsP_pu1_lowp->Sumw2();
   dEdX0noL1VsP_pu1_lowp2->Sumw2();
   dEdX0noL1VsP_pu2_lowp->Sumw2();
   dEdX0noL1VsP_pu2_lowp2->Sumw2();
   dEdX0noL1VsP_pu3_lowp->Sumw2();
   dEdX0noL1VsP_pu3_lowp2->Sumw2();
   dEdX0noL1VsP_pu4_lowp->Sumw2();
   dEdX0noL1VsP_pu4_lowp2->Sumw2();
   dEdXHDnoL1VsP_lowp->Sumw2();
   dEdXHDnoL1VsP_lowp2->Sumw2();
   dEdX0pixnoL1VsP_lowp->Sumw2();
   dEdX0pixnoL1VsP_lowp2->Sumw2();
   dEdXstripVsEta_lowp->Sumw2();
   dEstrVsdE_lowp->Sumw2();
   dEdXstripVsNhit_lowp->Sumw2();
   dEdXstripVsNhittrunc_lowp->Sumw2();
   dEdXstripVsCharge_lowp->Sumw2();
   EtaVsPhi_nhit->Sumw2();

   HSCP_dEdXpixVsstrip->Sumw2();
   HSCP_dEdX0pixVsstrip->Sumw2();
   HSCP_dEdXstripVsall->Sumw2();
   HSCP_dEdXpixVsall->Sumw2();

   Charge_pixl1->Sumw2();
   Charge_pixl2->Sumw2();
   Charge_pixl3->Sumw2();
   Charge_pixl4->Sumw2();
   Charge_pixd1->Sumw2();
   Charge_pixd2->Sumw2();
   Charge_pixd3->Sumw2();
   Charge_pixr1->Sumw2();
   Charge_pixr2->Sumw2();
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

   LowpCharge_tib1->Sumw2();
   LowpCharge_tib2->Sumw2();
   LowpCharge_tib3->Sumw2();
   LowpCharge_tib4->Sumw2();
   LowpCharge_tob1->Sumw2();
   LowpCharge_tob2->Sumw2();
   LowpCharge_tob3->Sumw2();
   LowpCharge_tob4->Sumw2();
   LowpCharge_tob5->Sumw2();
   LowpCharge_tob6->Sumw2();
   LowpCharge_tid1->Sumw2();
   LowpCharge_tid2->Sumw2();
   LowpCharge_tid3->Sumw2();
   LowpCharge_tec1->Sumw2();
   LowpCharge_tec2->Sumw2();
   LowpCharge_tec3->Sumw2();
   LowpCharge_tec4->Sumw2();
   LowpCharge_tec5->Sumw2();
   LowpCharge_tec6->Sumw2();
   LowpCharge_tec7->Sumw2();
   LowpCharge_tec8->Sumw2();
   LowpCharge_tec9->Sumw2();
   LowpCharge_pixl1->Sumw2();
   LowpCharge_pixl2->Sumw2();
   LowpCharge_pixl3->Sumw2();
   LowpCharge_pixl4->Sumw2();
   LowpCharge_pixd1->Sumw2();
   LowpCharge_pixd2->Sumw2();
   LowpCharge_pixd3->Sumw2();
   LowpCharge_pixr1->Sumw2();
   LowpCharge_pixr2->Sumw2();
   LowpCharge_Eta_pix->Sumw2();
   LowpCharge_Eta_strip->Sumw2();

   ChargeVsRun_pixl1->Sumw2();
   ChargeVsRun_pixl2->Sumw2();
   ChargeVsRun_pixl3->Sumw2();
   ChargeVsRun_pixl4->Sumw2();
   ChargeVsRun_pixd1->Sumw2();
   ChargeVsRun_pixd2->Sumw2();
   ChargeVsRun_pixd3->Sumw2();
   ChargeVsRun_pixr1->Sumw2();
   ChargeVsRun_pixr2->Sumw2();
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
   ZooChargeVsRun_pixr1->Sumw2();
   ZooChargeVsRun_pixr2->Sumw2();
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
   dEdXNoL1pixVsRun->Sumw2();
   dEdXNoL1VsRun->Sumw2();
   dEdXHiDropVsRun->Sumw2();
   dEdXpixHiDropVsRun->Sumw2();
   dEdXstripHiDropVsRun->Sumw2();
   dEdXHiDropNoL1VsRun->Sumw2();
   dEdX0VsRun->Sumw2();
   dEdX0pixVsRun->Sumw2();
   dEdX0stripVsRun->Sumw2();
   dEdX0NoL1pixVsRun->Sumw2();
   dEdX0NoL1VsRun->Sumw2();
   MassStripVsRun->Sumw2();
   MassNoL1VsRun->Sumw2();
   dEdX4VsRun->Sumw2();
   dEdX4pixVsRun->Sumw2();
   dEdX4stripVsRun->Sumw2();
   dEdX40VsRun->Sumw2();
   dEdX40pixVsRun->Sumw2();
   dEdX40stripVsRun->Sumw2();
   bg_dEdX0NoL1VsRun->Sumw2();
   iasNoL1VsRun->Sumw2();
   iasAllVsRun->Sumw2();

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
   HSCP_dEdX0->Sumw2(); 
   HSCP_dEdX0pix->Sumw2(); 
   HSCP_dEdX0strip->Sumw2(); 
   HSCP_MassIh->Sumw2(); 
   HSCP_MassIh0->Sumw2(); 
   HSCP_MassIhstrip->Sumw2(); 
   HSCP_MassIh0noL1->Sumw2(); 
   HSCP_MassIh0noL1_2->Sumw2(); 
   HSCP_MassIh0noL1_3->Sumw2(); 
   HSCP_MassIh0noL1_11->Sumw2(); 
   HSCP_MassIh0noL1_12->Sumw2(); 
   HSCP_MassIh0noL1_13->Sumw2(); 
   HSCP_MassIh0noL1_1s1->Sumw2(); 
   HSCP_MassIh0noL1_1s2->Sumw2(); 
   HSCP_MassIh0noL1_1s3->Sumw2(); 
   HSCP_MassIh0noL1_2s1->Sumw2(); 
   HSCP_MassIh0noL1_2s2->Sumw2(); 
   HSCP_MassIh0noL1_2s3->Sumw2(); 
   HSCP_MassIh0noL1_3s1->Sumw2(); 
   HSCP_MassIh0noL1_3s2->Sumw2(); 
   HSCP_MassIh0noL1_3s3->Sumw2(); 
   HSCP_MassIh0strip->Sumw2(); 
   HSCP_MassIhHDnoL1->Sumw2(); 
   HSCP_MassTOF->Sumw2(); 
   HSCP2d_MassTOFvsIh->Sumw2(); 
   HSCP2d_MassIh->Sumw2(); 
   HSCP2d_MassIh0->Sumw2(); 
   HSCP2d_MassIhstrip->Sumw2(); 
   HSCP2d_MassIh0noL1->Sumw2(); 
   HSCP2d_MassIh0strip->Sumw2(); 
   HSCP2d_MassIhHDnoL1->Sumw2(); 
   HSCP2d_Mass_pix_strip15->Sumw2(); 
   HSCP2d_Mass_pix_strip0->Sumw2(); 
   HSCP2d_Mass_pix_strip->Sumw2(); 
   HSCP_MassDiff_pix_strip0->Sumw2(); 
   HSCP_MassDiff_pix_strip15->Sumw2(); 
   HSCP_MassResol_pix_strip0->Sumw2(); 
   HSCP_MassResol_pix_strip15->Sumw2(); 
   HSCP_dEdXHiDrop->Sumw2();
   HSCP_dEdXstripHiDrop->Sumw2();
   HSCP_dEdXpixHiDrop->Sumw2();
   HSCP_dEdXHiDropNoL1->Sumw2();
   HSCP_dEdX0NoL1->Sumw2();
   lowp_MassIh->Sumw2(); 
   lowp_MassIh0->Sumw2(); 
   lowp_MassIhstrip->Sumw2(); 
   lowp_MassIh0noL1->Sumw2(); 
   lowp_MassIh0noL1_2->Sumw2(); 
   lowp_MassIh0noL1_3->Sumw2(); 
   lowp_MassIh0noL1_11->Sumw2(); 
   lowp_MassIh0noL1_12->Sumw2(); 
   lowp_MassIh0noL1_13->Sumw2(); 
   lowp_MassIh0noL1_1s1->Sumw2(); 
   lowp_MassIh0noL1_1s2->Sumw2(); 
   lowp_MassIh0noL1_1s3->Sumw2(); 
   lowp_MassIh0noL1_2s1->Sumw2(); 
   lowp_MassIh0noL1_2s2->Sumw2(); 
   lowp_MassIh0noL1_2s3->Sumw2(); 
   lowp_MassIh0noL1_3s1->Sumw2(); 
   lowp_MassIh0noL1_3s2->Sumw2(); 
   lowp_MassIh0noL1_3s3->Sumw2(); 
   lowp_MassIh0strip->Sumw2(); 
   lowp_MassIhHDnoL1->Sumw2(); 
   lowp2d_MassIh->Sumw2(); 
   lowp2d_MassIh0->Sumw2(); 
   lowp2d_MassIhstrip->Sumw2(); 
   lowp2d_MassIh0noL1->Sumw2(); 
   lowp2d_MassIh0strip->Sumw2(); 
   lowp2d_MassIhHDnoL1->Sumw2(); 
   lowp_dEdXpixVsstrip->Sumw2();
   lowp_dEdX0pixVsstrip->Sumw2();
   lowp2d_Mass_pix_strip15->Sumw2(); 
   lowp2d_Mass_pix_strip0->Sumw2(); 
   bg_lowp2d_Mass_pix_strip0->Sumw2(); 
   bg_transf_Mass->Sumw2(); 
   lowp_MassDiff_pix_strip0->Sumw2(); 
   lowp_MassDiff_pix_strip15->Sumw2(); 
   bg_dEdXVsIas->Sumw2(); 

   HSCP_FMIP4->Sumw2(); 
   HSCP_FMIP3p5->Sumw2(); 
   HSCP_FMIP3p2->Sumw2(); 
   HSCP_iasnol1->Sumw2(); 
   HSCP_iasall->Sumw2(); 
   HSCP_probQ->Sumw2();
   HSCP_probQNoL1->Sumw2();
   HSCP_probXY->Sumw2();
   HSCP_probXYNoL1->Sumw2();
   probQVsRun->Sumw2();
   probQNoL1VsRun->Sumw2();
   probXYVsRun->Sumw2();
   probXYNoL1VsRun->Sumw2();
   probQVsIas->Sumw2();

   HSCP_pt->Sumw2(); 
   HSCP_eta->Sumw2(); 
   HSCP_iso_eop->Sumw2(); 
   nPV->Sumw2(); 
   HSCP_invB->Sumw2(); 
   HSCP_errinvB->Sumw2(); 
   HSCP_invBDT->Sumw2(); 
   HSCP_invBCSC->Sumw2(); 
   HSCP_time->Sumw2(); 
   HSCP_npix->Sumw2(); 
   HSCP_nstrip->Sumw2(); 
   HSCP_nmpix->Sumw2(); 
   HSCP_nmstrip->Sumw2(); 
   HSCP_nratio->Sumw2(); 
   HSCP_nmratio->Sumw2(); 

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

   TString outputfilename="analysis_ul";
   if (year==2016) outputfilename+="_2016";
   else if (year==2017) outputfilename+="_2017";
   else if (year==2018) outputfilename+="_2018";
   outputfilename+=Letter;
   if (!dataFlag) outputfilename+="_MC";
   outputfilename+="_2mars.root";
   TFile* OutputHisto = new TFile(outputfilename,"RECREATE");
   TString templateFileName="template";
   if (year==2016) templateFileName+="_2016";
   else if (year==2017) templateFileName+="_2017";
   else if (year==2018) templateFileName+="_2018";
   templateFileName+=Letter;
   if (!dataFlag) templateFileName+="_MC";
   templateFileName+="_2mars.root";
   TFile* OutputTemplate;
   if (writeTemplateOnDisk) OutputTemplate = new TFile(templateFileName,"RECREATE");


//   loadSFPixelCalib();
   loadSFPixelTamas();
   bool PixelCorr2Apply = true;
   if (!dataFlag) PixelCorr2Apply = false;  // No Pixel correction for MC !


//   loadDeDxTemplates

    if (boolDeDxTemp) {
     if (!dataFlag) {
//     dEdxTemplatesNoL1 = loadDeDxTemplate("templateMC_w17_MC_21jan.root", "Charge_Vs_Path_noL1_NoM", true);
      if (year==2017) {
       dEdxTemplatesNoL1 = loadDeDxTemplate("templateMC_w17_MC_21jan.root", "Charge_Vs_Path_noL1", true);
       dEdxTemplatesAll = loadDeDxTemplate("templateMC_w17_MC_21jan.root", "Charge_Vs_Path", true);
      }
      else if (year==2018){
       dEdxTemplatesNoL1 = loadDeDxTemplate("template_2018MC_w18_MC_28feb.root", "Charge_Vs_Path_noL1", true);
       dEdxTemplatesAll = loadDeDxTemplate("template_2018MC_w18_MC_28feb.root", "Charge_Vs_Path", true);
      }
     }
     else {
     // note : for data, we should write something smart to read template per era
      std::string TemplateDataName = "template";
      if (year>2016) { 
//      if (year==2016) TemplateDataName+="_2016";
       if (year==2017) TemplateDataName+="_2017";
       else if (year==2018) TemplateDataName+="_2018";
       TemplateDataName+=Letter;
       if (year==2017) TemplateDataName+="_21jan.root";
       else if (year==2018) TemplateDataName+="_28feb.root";
      }
      else { 
        TemplateDataName+="_2017B_21jan.root";
      }

//      dEdxTemplatesNoL1 = loadDeDxTemplate(TemplateDataName, "Charge_Vs_Path_noL1_NoM", true);
      dEdxTemplatesNoL1 = loadDeDxTemplate(TemplateDataName, "Charge_Vs_Path_noL1", true);
      dEdxTemplatesAll = loadDeDxTemplate(TemplateDataName, "Charge_Vs_Path", true);
     }
    }


//   nentries = 200000;
   cout << "run on  " << nentries << " entries " << endl;


   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      if(jentry%1000000 ==0 && jentry!=0) cout << " number of processed events is " << jentry <<  " = " << (100.*jentry)/(1.*nentries) << "%" <<endl;
//      if(jentry%10000 ==0 && jentry!=0) cout << " number of processed events is " << jentry <<  " = " << (100.*jentry)/(1.*nentries) << "%" <<endl;

      if (!dataFlag) runNumber = 300000;

//      if (runNumber<305150 || runNumber>305300) continue;
//      if (runNumber>323232) continue;  // because 2018 is rereco at the moment --> remove this if when running on UL

      // only Muon Trigger at the moment (because I am running on SingleMu data)
      if (!hlt_mu50 && !hlt_tkmu100 && !hlt_oldmu100) continue;

      nPVVsRun->Fill(runNumber,npv);
      nPV->Fill(npv);

      if (npv<1) continue;

      // Loop on HSCP candidate
      for (int ihs=0; ihs<nhscp; ihs++) {
          int index_of_the_track=hscp_track_idx[ihs];
          if (index_of_the_track>-1) { 
             bool selection=true;
             if (track_pt[index_of_the_track]<55) selection=false;
             if (abs(track_eta[index_of_the_track])>2.1) selection=false;

             if (track_nvalidhits[index_of_the_track]<8) selection=false;
             if (track_npixhits[index_of_the_track]<2) selection=false;
             if (track_validfraction[index_of_the_track]<0.8) selection=false;
             if (track_missing[index_of_the_track]>99999) selection=false;
             if (track_validlast[index_of_the_track]<-99999) selection=false;

             
             bool is_high_qual =  (track_qual[index_of_the_track] & (1 << TrackQuality::highPurity)) >> TrackQuality::highPurity ;
             if (!is_high_qual) selection=false;
             if (track_chi2[index_of_the_track]>5) selection=false;

             // which selection for dz and dyx w/r to the primary vertex ? here 1st PV to pass the selection
             if (abs(track_dxy[index_of_the_track])>0.5) selection=false;
             if (abs(track_dz[index_of_the_track])>0.5) selection=false;

             // which DR used for iso ? here 0.3
             if (hscp_iso2_tk[ihs]>50) selection=false;
             float eop=(hscp_iso2_ecal[ihs] + hscp_iso2_hcal[ihs])/track_p[index_of_the_track];
             if (eop>0.3) selection=false;
             // no cut on relative iso : false if hscpIso.Get_TK_SumEt() / track->pt() > 9999999

             // no cut on sigma pT/pT for signal yet :
             //if (track_pterr[index_of_the_track]/track_pt[index_of_the_track]>0.25) selection=false; // commented on September 22, 2021
             
             // no cut in TOF yet : false if tof->nDof() < 8
             //     &&  (dttof->nDof() < 6 || csctof->nDof() < 6)
             // no cut in TOF yet : false if tof->inverseBetaErr() > GlobalMaxTOFErr
            
             

             // no cut yet :  false if  dedxSObj->numberOfMeasurements() < 6 ;
             
             if (selection) {
                ptVsRun->Fill(runNumber,track_pt[index_of_the_track]);

                std::vector <float> charge_corr;
                std::vector <float> pathlength;
                std::vector <int> subdetId;
                std::vector <int> moduleGeometry;
                std::vector <bool> bool_cleaning;
                std::vector <bool> mustBeInside;
             
                std::vector <float> charge_corr1;
                std::vector <float> pathlength1;
                std::vector <int> subdetId1;
                std::vector <UInt_t> detId1;
                std::vector <int> moduleGeometry1;
                std::vector <bool> bool_cleaning1;
                std::vector <bool> mustBeInside1;

                std::vector <float> charge_corr2;
                std::vector <float> pathlength2;
                std::vector <int> subdetId2;
                std::vector <int> moduleGeometry2;
                std::vector <bool> bool_cleaning2;
                std::vector <bool> mustBeInside2;

                std::vector <float> charge_corr3;
                std::vector <float> pathlength3;
                std::vector <int> subdetId3;
                std::vector <UInt_t> detId3;
                std::vector <int> moduleGeometry3;
                std::vector <bool> bool_cleaning3;
                std::vector <bool> mustBeInside3;

                std::vector <float> charge_corr5;
                std::vector <float> pathlength5;
                std::vector <int> subdetId5;
                std::vector <UInt_t> detId5;
                std::vector <int> moduleGeometry5;
                std::vector <bool> bool_cleaning5;
                std::vector <bool> mustBeInside5;

                int nstip_=0;
                int npix_=0;
                for (int iclu=track_index_hit[hscp_track_idx[ihs]]; iclu<track_index_hit[hscp_track_idx[ihs]]+track_nhits[hscp_track_idx[ihs]]; iclu++) {
                     float ch1=dedx_charge[iclu];
                     bool clean1=true;
                     if (dedx_subdetid[iclu]>=3) {
                        // strip
                        // without any correction :
//                        ch1=sclus_charge[iclu];  // charge without any correction
//                        clean1=sclus_clusclean[iclu];  
//                      // with Saturation only (but no Xtalk inversion) :
                          nstip_++;
                          float check_charge=0;
                          vector<int> Quncor;
                          for (int istrip=sclus_index_strip[iclu]; istrip<sclus_index_strip[iclu]+sclus_nstrip[iclu]; istrip++) {
                            check_charge+=strip_ampl[istrip];
                            Quncor.push_back(strip_ampl[istrip]);
                          }
                          float deltaDif=check_charge-ch1;
                          if (deltaDif<0) deltaDif*=-1;
                          if (deltaDif>0.001) std::cout << "difference dans le calcul de la charge " << ch1 << " " << check_charge << " --> probleme acces ampl ???? " << std::endl; 
                          vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25);
                          float newcharge =0;
                          for (unsigned int inwc=0; inwc<Qcor.size(); inwc++) { newcharge+=Qcor[inwc]; }
                          ch1=newcharge;
                          clean1=sclus_clusclean2[iclu]; // clusterCleaning with Xtalk inversion and Saturation (September 22, 2021)
                          // for the record : there is a possibility to recompute the cluster cleaning with the Saturation only : clusterCleaning(Qcor, 1, &exitCode)
                          // but it is not what we want here
                     }
                     else {
                        // pixel
                          // float scaling =GetSFPixel(dedx_subdetid[iclu], dedx_detid[iclu], year, abs(track_eta[hscp_track_idx[ihs]]), runNumber);
                          npix_++;
                          if (PixelCorr2Apply) {
                            float scaling =GetSFPixelTamas(dedx_subdetid[iclu], dedx_detid[iclu], year, runNumber);
                            ch1*=scaling;
                          }
                     }
                     if (clean1 && dedx_insideTkMod[iclu]) {
                       // fill both Strip and Pixel
                       charge_corr.push_back(ch1);
                       pathlength.push_back(dedx_pathlength[iclu]);
                       subdetId.push_back(dedx_subdetid[iclu]);
                       moduleGeometry.push_back(dedx_modulgeom[iclu]);
                       mustBeInside.push_back(dedx_insideTkMod[iclu]);
                       bool_cleaning.push_back(clean1);

                       if (dedx_isstrip[iclu]) {
                        // fill only Strip
                        charge_corr1.push_back(ch1);
                        pathlength1.push_back(dedx_pathlength[iclu]);
                        subdetId1.push_back(dedx_subdetid[iclu]);
                        detId1.push_back(dedx_detid[iclu]);
                        moduleGeometry1.push_back(dedx_modulgeom[iclu]);
                        mustBeInside1.push_back(dedx_insideTkMod[iclu]);
                        bool_cleaning1.push_back(clean1);
                       }
                       else {
                        // fill only Pixel
                        charge_corr2.push_back(ch1);
                        pathlength2.push_back(dedx_pathlength[iclu]);
                        subdetId2.push_back(dedx_subdetid[iclu]);
                        moduleGeometry2.push_back(dedx_modulgeom[iclu]);
                        mustBeInside2.push_back(dedx_insideTkMod[iclu]);
                        bool_cleaning2.push_back(clean1);
                       }
                     
                       bool no_in_L1_pixel =true;
                       int info_layr=GetLayerLabel(dedx_subdetid[iclu], dedx_detid[iclu],year);
                       if (info_layr==23) no_in_L1_pixel=false;
                       if (no_in_L1_pixel) {
                         // fill both Strip and Pixel but without L1 Pixel
                         charge_corr3.push_back(ch1);
                         pathlength3.push_back(dedx_pathlength[iclu]);
                         subdetId3.push_back(dedx_subdetid[iclu]);
                         detId3.push_back(dedx_detid[iclu]);
                         moduleGeometry3.push_back(dedx_modulgeom[iclu]);
                         mustBeInside3.push_back(dedx_insideTkMod[iclu]);
                         bool_cleaning3.push_back(clean1);
                         if (dedx_subdetid[iclu]<3) {  // pixel only
                            // fill only Pixel, without L1 Pixel
                            charge_corr5.push_back(ch1);
                            pathlength5.push_back(dedx_pathlength[iclu]);
                            subdetId5.push_back(dedx_subdetid[iclu]);
                            detId5.push_back(dedx_detid[iclu]);
                            moduleGeometry5.push_back(dedx_modulgeom[iclu]);
                            mustBeInside5.push_back(dedx_insideTkMod[iclu]);
                            bool_cleaning5.push_back(clean1);
                         }
                       }
                    }


                } // end loop iclu


                if (charge_corr.size()>6) {  // the one for which we don't remove any hits
                  int nval1=0;
                  int nsatv1=0;
                  int nval1_0=0;
                  int nsatv1_0=0;
                  // Ih 15% drop of lowest values
                  double ih_LD        = getdEdX(charge_corr,    pathlength,  subdetId,  moduleGeometry,  bool_cleaning,  mustBeInside,  dEdxSF, NULL,2, 0.15, nval1, nsatv1);
                  dEdXVsRun->Fill(       runNumber,ih_LD);
                  NmeasVsRun->Fill(runNumber,nval1);
                  if (nsatv1>9) nsatv1=9;
                  // Ih no drop
                  double ih0_cor     = getdEdX(charge_corr,  pathlength,  subdetId,  moduleGeometry,  bool_cleaning,  mustBeInside,  dEdxSF, NULL,2, 0., nval1_0, nsatv1_0);
                  dEdX0VsRun->Fill(runNumber,ih0_cor);
                  int nval1HD=0;
                  int nsatv1HD=0;
                  // Ih  15% drop of highest values
                  double ih_corHiDrop = getdEdX(charge_corr,  pathlength,  subdetId,  moduleGeometry,  bool_cleaning,  mustBeInside,  dEdxSF, NULL,2, 0., 0.15, nval1HD, nsatv1HD);
                  dEdXHiDropVsRun->Fill(     runNumber,ih_corHiDrop);

                  int nval2=0;
                  int nval2_0=0;
                  int nsatv2=0;
                  int nsatv2_0=0;
                  int nval2HD=0;
                  int nsatv2HD=0;
                  // Ih 15% drop of lowest values
                  double ih_LDstrip   = getdEdX(charge_corr1,   pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0.15, nval2, nsatv2);
                  dEdXstripVsRun->Fill(  runNumber,ih_LDstrip);
                  NmeasStrVsRun->Fill(runNumber,nval2);
                  if (nsatv2>9) nsatv2=9;
                  NsatVsRun->Fill(runNumber,nsatv1);
                  NsatStrVsRun->Fill(runNumber,nsatv2);
                  // Ih no drop
                  double ih0_strip   = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0., nval2_0, nsatv2_0);
                  dEdX0stripVsRun->Fill(runNumber,ih0_strip);


                  // Ih  15% drop of highest values
                  double ih_stripHiDrop = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,2, 0., 0.15, nval2HD, nsatv2HD);
                  dEdXstripHiDropVsRun->Fill(runNumber,ih_stripHiDrop );

                  int nval3=0;
                  int nval3_0=0;
                  int nsatv3=0;
                  int nsatv3_0=0;
                  int nval3HD=0;
                  int nsatv3HD=0;
                  // Ih 15% drop of lowest values
                  double ih_LDpix     = getdEdX(charge_corr2,   pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0.15, nval3, nsatv3);
                  dEdXpixVsRun->Fill(    runNumber,ih_LDpix);
                  NmeasPixVsRun->Fill(runNumber,nval3);
                  if (nsatv3>9) nsatv3=9;
                  NsatPixVsRun->Fill(runNumber,nsatv3);
                  // Ih no drop
                   double ih0_pix     = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0., nval3_0, nsatv3_0);
                  dEdX0pixVsRun->Fill(runNumber,ih0_pix);



                  // Ih  15% drop of highest values
                  double ih_pixHiDrop   = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,2, 0., 0.15, nval3HD, nsatv3HD);
                  dEdXpixHiDropVsRun->Fill(  runNumber,ih_pixHiDrop);


                  int nval6=0;
                  int nval6_0=0;
                  int nsatv6=0;
                  int nsatv6_0=0;
                  // Ih 15% drop of lowest values
                  double ih_LDnoL1pix = getdEdX(charge_corr5, pathlength5, subdetId5, moduleGeometry5, bool_cleaning5, mustBeInside5, dEdxSF, NULL,2, 0.15, nval6, nsatv6);
                  dEdXNoL1pixVsRun->Fill(runNumber,ih_LDnoL1pix);
                  // Ih no drop
                  double ih0_noL1pix = getdEdX(charge_corr5, pathlength5, subdetId5, moduleGeometry5, bool_cleaning5, mustBeInside5, dEdxSF, NULL,2, 0., nval6_0, nsatv6_0);
                  dEdX0NoL1pixVsRun->Fill(runNumber,ih0_noL1pix);
                 

                  float mass_ih0noL1=-1;
                  float mass_ihHDnoL1=-1;
                  float mass_ih0noL1_2=-1;
                  float mass_ih0noL1_3=-1;


                  // Ias
                  int nval20_0=0;
                  int nsat20_0=0;
                  double ias_all =-1;
                  if (boolDeDxTemp)  ias_all=getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);
                  iasAllVsRun->Fill(   runNumber,ias_all);


                  nval20_0=0;
                  nsat20_0=0;
                  double ias_noL1 = -1;
                  if (boolDeDxTemp) ias_noL1= getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, dEdxTemplatesNoL1,2, 0., nval20_0, nsat20_0);
                  iasNoL1VsRun->Fill(   runNumber,ias_noL1);

                  float mass_strip=-1;
                  if (ih_LDstrip - Cval_ldstrip>0) mass_strip=sqrt((ih_LDstrip - Cval_ldstrip)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_ldstrip);

                  float mass_ih=-1;
                  if (ih_LD - Cval_ld>0) mass_ih=sqrt((ih_LD - Cval_ld)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_ld);

                  float mass_ih0=-1;
                  if (ih0_cor - Cval_all>0) mass_ih0= sqrt((ih0_cor - Cval_all)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_all);


                  float mass_0strip=-1;
                  if (ih0_strip - Cval_strip>0) mass_0strip= sqrt((ih0_strip - Cval_strip)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_strip);


                  float mass_pixel=-1;
                  if (ih_LDpix - Cval_pix>0) mass_pixel= sqrt((ih_LDpix - Cval_pix)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_pix);

                  float mass_pixnL1=-1;
                  if (ih0_noL1pix - Cval_pixnol1>0) mass_pixnL1= sqrt((ih0_noL1pix - Cval_pixnol1)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_pixnol1);


                  // additional cut on the number of measurement points
                  if (charge_corr3.size() >=6) {

                     int nval4=0;
                     int nval4_0=0;
                     int nsatv4=0;
                     int nsatv4_0=0;


                     // Ih 15% drop of lowest values
                     double ih_LDnoL1    = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, NULL,2, 0.15, nval4, nsatv4);
                     dEdXNoL1VsRun->Fill(   runNumber,ih_LDnoL1);
 
                     // Ih no drop
                     double ih0_noL1 = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, NULL,2, 0., nval4_0, nsatv4_0);

                     dEdX0NoL1VsRun->Fill(runNumber,ih0_noL1);  //  <==== our best estimate for Ih

                     Nmeas0VsRun->Fill(runNumber,nval4_0);    // no L1
                     NmeasPix0VsRun->Fill(runNumber,nval6_0); // no L1
                     NmeasStr0VsRun->Fill(runNumber,nval2_0);  

                     // Ih  15% drop of highest values
                     nval1HD=0;
                     nsatv1HD=0;
                     double ih_corHiDropNoL1 = getdEdX(charge_corr3,  pathlength3,  subdetId3,  moduleGeometry3,  bool_cleaning3,  mustBeInside3,  dEdxSF, NULL,2, 0., 0.15, nval1HD, nsatv1HD);
                     dEdXHiDropNoL1VsRun->Fill(     runNumber,ih_corHiDropNoL1);

                     if (ih0_noL1 - Cval_nol1>0) mass_ih0noL1= sqrt((ih0_noL1 - Cval_nol1)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1);
                     if (ih0_noL1 - Cval_nol1_2>0) mass_ih0noL1_2= sqrt((ih0_noL1 - Cval_nol1_2)*track_p[hscp_track_idx[ihs]]*track_p[hscp_track_idx[ihs]]/Kval_nol1_2);
                     float val_nol1_3 = Cval_nol1_3 + (0.92*0.92*Kval_nol1_3)/(6.5*6.5) +Nval_nol1_3*log(6.5/0.92);
                     if (ih0_noL1 - val_nol1_3 >0 && computeSpecial) 
                                          mass_ih0noL1_3= getMassSpecial(ih0_noL1,track_p[hscp_track_idx[ihs]], Kval_nol1_3, Cval_nol1_3,Nval_nol1_3);
                     if (ih_corHiDropNoL1 - Cval_hdnol1>0) mass_ihHDnoL1= sqrt((ih_corHiDropNoL1 - Cval_hdnol1)*track_p[index_of_the_track]*track_p[index_of_the_track]/Kval_hdnol1);

                     float fmip_strip4 =   FMIP(charge_corr1, pathlength1,dEdxSF[0],4);
                     float fmip_strip3p5 = FMIP(charge_corr1, pathlength1,dEdxSF[0], 3.5);
                     float fmip_strip3p2 = FMIP(charge_corr1, pathlength1,dEdxSF[0], 3.2);

                     FMIP4VsRun->Fill(   runNumber,fmip_strip4);
                     FMIP3p5VsRun->Fill( runNumber,fmip_strip3p5);
                     FMIP3p2VsRun->Fill( runNumber,fmip_strip3p2);

                     if (nsatv4_0>9) nsatv4_0=9;
                     if (nsatv2_0>9) nsatv2_0=9;
                     if (nsatv6_0>9) nsatv6_0=9;
                     Nsat0VsRun->Fill(runNumber,nsatv4_0);
                     NsatPix0VsRun->Fill(runNumber,nsatv6_0);
                     NsatStr0VsRun->Fill(runNumber,nsatv2_0);


                     if ((!blind_data) || (mass_ih0noL1_2<500)) {

                       HSCP_pt->Fill(track_pt[index_of_the_track]);
                       HSCP_eta->Fill(track_eta[index_of_the_track]);
                       HSCP_iso_eop->Fill(eop);
                       if (boolILumi) lumiVsRun->Fill(runNumber,InstLumi);
                       HSCP_dEdX->Fill(ih_LD);
                       HSCP_dEdX0->Fill(ih0_cor);
                       HSCP_dEdXHiDrop->Fill(ih_corHiDrop);
                       HSCP_dEdXstrip->Fill(ih_LDstrip);
                       HSCP_dEdX0strip->Fill(ih0_strip);
                       HSCP_dEdXHiDropNoL1->Fill(ih_corHiDropNoL1);
                       // comparison
                       HSCP_dEdXstripHiDrop->Fill(ih_stripHiDrop);
                       HSCP_dEdX0pix->Fill(ih0_pix);
                       HSCP_dEdXpix->Fill(ih_LDpix);
                       HSCP_dEdXpixHiDrop->Fill(ih_pixHiDrop);
                       HSCP_dEdXpixVsstrip->Fill(ih_LDstrip,ih_LDpix);
                       HSCP_dEdXstripVsall->Fill(ih_LD,ih_LDstrip);
                       HSCP_dEdXpixVsall->Fill(ih_LD,ih_LDpix);
  
                       // our best Ih estimate :
                       HSCP_dEdX0NoL1->Fill(ih0_noL1);
                       HSCP_dEdX0pixVsstrip->Fill(ih0_strip,ih0_noL1pix);


                       HSCP_iasall->Fill(ias_all);
                       HSCP_iasnol1->Fill(ias_noL1);
                       if (boolProbQ) {
                        HSCP_probQ->Fill(track_probQ[index_of_the_track]);
                        HSCP_probQNoL1->Fill(track_probQNoL1[index_of_the_track]);
                        HSCP_probXY->Fill(track_probXY[index_of_the_track]);
                        HSCP_probXYNoL1->Fill(track_probXYNoL1[index_of_the_track]);
                        probQVsRun->Fill(   runNumber,track_probQ[index_of_the_track]);
                        probQNoL1VsRun->Fill(   runNumber,track_probQNoL1[index_of_the_track]);
                        probXYVsRun->Fill(   runNumber,track_probXY[index_of_the_track]);
                        probXYNoL1VsRun->Fill(   runNumber,track_probXYNoL1[index_of_the_track]);
                        probQVsIas->Fill(ias_all,track_probQ[index_of_the_track]);
                       }

                       HSCP_nstrip->Fill(nstip_);
                       HSCP_npix->Fill(npix_);
                       if (nstip_>0) HSCP_nratio->Fill(npix_*1./nstip_);
                       HSCP_nmstrip->Fill(charge_corr1.size());
                       HSCP_nmpix->Fill(charge_corr5.size());
                       if (charge_corr1.size()>0) HSCP_nmratio->Fill(charge_corr5.size()*1./charge_corr1.size());

                       HSCP_MassIh0noL1->Fill(mass_ih0noL1);
                       HSCP_MassIh0noL1_2->Fill(mass_ih0noL1_2);
                       HSCP_MassIh0noL1_3->Fill(mass_ih0noL1_3);

                       if (ih0_noL1> Cval_nol1) {
                       HSCP_MassIh0noL1_11->Fill(mass_ih0noL1);
                       HSCP_MassIh0noL1_12->Fill(mass_ih0noL1_2);
                       HSCP_MassIh0noL1_13->Fill(mass_ih0noL1_3);
                       MassNoL1VsRun->Fill(runNumber,mass_ih0noL1_2);
                       }
                       if (ih0_noL1>3.27 + 0.21) {
                       HSCP_MassIh0noL1_1s1->Fill(mass_ih0noL1);
                       HSCP_MassIh0noL1_1s2->Fill(mass_ih0noL1_2);
                       HSCP_MassIh0noL1_1s3->Fill(mass_ih0noL1_3);
                         if (ih0_noL1>3.27 + 2*0.21) {
                         HSCP_MassIh0noL1_2s1->Fill(mass_ih0noL1);
                         HSCP_MassIh0noL1_2s2->Fill(mass_ih0noL1_2);
                         HSCP_MassIh0noL1_2s3->Fill(mass_ih0noL1_3);
                           if (ih0_noL1>3.27 + 3*0.21) {
                           HSCP_MassIh0noL1_3s1->Fill(mass_ih0noL1);
                           HSCP_MassIh0noL1_3s2->Fill(mass_ih0noL1_2);
                           HSCP_MassIh0noL1_3s3->Fill(mass_ih0noL1_3);
                           }
                         }
                       }

                       HSCP_MassIhHDnoL1->Fill(mass_ihHDnoL1);
                       HSCP2d_MassIh0noL1->Fill(track_p[index_of_the_track],mass_ih0noL1);
                       HSCP2d_MassIhHDnoL1->Fill(track_p[index_of_the_track],mass_ihHDnoL1);
                       HSCP_MassIhstrip->Fill(mass_strip);
                       HSCP_MassIh->Fill(mass_ih);
                       HSCP_MassIh0->Fill(mass_ih0);
                       HSCP_MassIh0strip->Fill(mass_0strip);
                       HSCP2d_MassIhstrip->Fill(track_p[index_of_the_track],mass_strip);
                       HSCP2d_MassIh->Fill(track_p[index_of_the_track],mass_ih);
                       HSCP2d_MassIh0->Fill(track_p[index_of_the_track],mass_ih0);
                       HSCP2d_MassIh0strip->Fill(track_p[index_of_the_track],mass_0strip);

                       HSCP2d_Mass_pix_strip15->Fill(mass_strip,mass_pixel);
                       HSCP2d_Mass_pix_strip0->Fill(mass_0strip,mass_pixnL1);
                       HSCP2d_Mass_pix_strip->Fill(mass_0strip,mass_pixnL1);
                       if (mass_0strip>0 && mass_pixnL1>0 ) {
                        HSCP_MassDiff_pix_strip0->Fill(mass_0strip - mass_pixnL1);
                        HSCP_MassResol_pix_strip0->Fill((mass_0strip - mass_pixnL1)/mass_0strip);
                       }
                       if (mass_strip>0 && mass_pixel>0) {
                        HSCP_MassDiff_pix_strip15->Fill(mass_strip - mass_pixel);
                        HSCP_MassResol_pix_strip15->Fill((mass_strip - mass_pixel)/mass_strip);
                       }

                       MassStripVsRun->Fill(runNumber,mass_strip);

                       FMIP4VsEta->Fill(track_eta[index_of_the_track],fmip_strip4);

                       HSCP_FMIP4->Fill(fmip_strip4);
                       HSCP_FMIP3p5->Fill(fmip_strip3p5);
                       HSCP_FMIP3p2->Fill(fmip_strip3p2);
                     }




/*
                     if (boolILumi) {
                      dEdXVsIL->Fill(InstLumi,ih_LD);
                      dEdXpixVsIL->Fill(InstLumi,ih_pix);
                      dEdXstripVsIL->Fill(InstLumi,ih_strip);
                     }
*/

                  } // end cut on #meas charge_corr3.size


                  if (runNumber==305186) {
                        R1_StdEdXVsEvent->Fill(event,ih_LDstrip);
                        if (boolILumi) {
                         R1_StdEdXVsLumi->Fill(InstLumi,ih_LDstrip);
                         R1_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R1_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R1_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R1_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
                  }
                  else if (runNumber==305188) {
                        R2_StdEdXVsEvent->Fill(event,ih_LDstrip);
                        if (boolILumi) {
                         R2_StdEdXVsLumi->Fill(InstLumi,ih_LDstrip);
                         R2_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R2_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R2_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R2_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
                  }
                  else if (runNumber==305204) {
                        R3_StdEdXVsEvent->Fill(event,ih_LDstrip);
                        if (boolILumi) {
                         R3_StdEdXVsLumi->Fill(InstLumi,ih_LDstrip);
                         R3_LumiVsEvent->Fill(event,InstLumi);
                        }
                        R3_nPVVsEvent->Fill(event,npv);
                        int binEv = (int) event/4000000 ;
                        int nCandperEv = R3_CandVsEvent->GetBinContent(binEv+1) + 1;
                        R3_CandVsEvent->SetBinContent(binEv+1, nCandperEv);
		  }
                    
                  // Test with a power 4 in place of 2 in the formula for Ih, 15% drop of the lowest values
                  double ih4_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,4, 0.15,  nval1, nsatv1);
                  double ih4_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,4, 0.15, nval2, nsatv2);
                  double ih4_pix = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,4, 0.15, nval3, nsatv3);

                  double ih40_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, NULL,4, 0., nval1_0, nsatv1_0);
                  double ih40_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF, NULL,4, 0., nval2_0, nsatv2_0);
                  double ih40_pix = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF, NULL,4, 0.,nval3_0, nsatv3_0);

                  dEdX4VsRun->Fill(runNumber,ih4_cor);
                  dEdX4pixVsRun->Fill(runNumber,ih4_pix);
                  dEdX4stripVsRun->Fill(runNumber,ih4_strip);
                  dEdX40VsRun->Fill(runNumber,ih40_cor);
                  dEdX40pixVsRun->Fill(runNumber,ih40_pix);
                  dEdX40stripVsRun->Fill(runNumber,ih40_strip);

                  if (hscp_muon_idx[ihs]>-1 && 
                       charge_corr3.size() >=6) {
                       // un muon existe pour le candidat HSCP
             // no cut in TOF yet : false if tof->nDof() < 8
             //     &&  (dttof->nDof() < 6 || csctof->nDof() < 6)
             // no cut in TOF yet : false if tof->inverseBetaErr() > GlobalMaxTOFErr
                       bool tof_sel= true;
                       if (muon_comb_tofndof[hscp_muon_idx[ihs]] <8 &&
                          (muon_dt_tofndof[hscp_muon_idx[ihs]]<6 || muon_csc_tofndof[hscp_muon_idx[ihs]]<6)) tof_sel = false;
                       if (muon_comb_inversebetaerr[hscp_muon_idx[ihs]]>0.15 ) tof_sel = false;
                       if (fabs(muon_comb_inversebeta[hscp_muon_idx[ihs]]-1)>=50) tof_sel = false;  // not in the official selection in 
                       if (tof_sel) {
                            // https://github.com/enibigir/SUSYBSMAnalysis-HSCP/blob/dev/Analyzer/plugins/Analyzer.cc  ???
/*
 // what I did before
                       if (muon_comb_tofndof[hscp_muon_idx[ihs]]>=8
                             && (muon_dt_tofndof[hscp_muon_idx[ihs]]>=6 || muon_csc_tofndof[hscp_muon_idx[ihs]]>=6) 
                             && muon_comb_inversebetaerr[hscp_muon_idx[ihs]]<=0.15 
                             && fabs(muon_comb_inversebeta[hscp_muon_idx[ihs]]-1)<50) {
*/

                         if (muon_comb_inversebeta[hscp_muon_idx[ihs]]!=0.) {
                          float beta_for_mass = 1./muon_comb_inversebeta[hscp_muon_idx[ihs]];
                          float gamma_for_mass = 1/sqrt(1-beta_for_mass*beta_for_mass);
                          float mass_tof = track_p[index_of_the_track]/(beta_for_mass*gamma_for_mass) ;
//                          if (!blind_data) {
                          if ((!blind_data) || (mass_ih0noL1_2<500)) {
                              HSCP_MassTOF->Fill(mass_tof); 
                              HSCP2d_MassTOFvsIh->Fill(mass_tof,mass_ih0noL1); 
                          }
                         }

                         invBVsRun->Fill(runNumber,muon_comb_inversebeta[hscp_muon_idx[ihs]]);
                         errinvBVsRun->Fill(runNumber,muon_comb_inversebetaerr[hscp_muon_idx[ihs]]);
                         timeVsRun->Fill(runNumber,muon_comb_vertextime[hscp_muon_idx[ihs]]);
                         if ((!blind_data) || (mass_ih0noL1_2<500)) {
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
                       }
                       if (year==2016) {
                        bool tof_sel= true;
                        if (muon_newcomb_tofndof[hscp_muon_idx[ihs]] <8 &&
                          (muon_newdt_tofndof[hscp_muon_idx[ihs]]<6 || muon_newcsc_tofndof[hscp_muon_idx[ihs]]<6)) tof_sel = false;
                        if (muon_newcomb_inversebetaerr[hscp_muon_idx[ihs]]>0.15 ) tof_sel = false;
                        if (fabs(muon_newcomb_inversebeta[hscp_muon_idx[ihs]]-1)>=50) tof_sel = false;  // not in the official selection in 
                        if (tof_sel) {

                         if ((!blind_data) || (mass_ih0noL1_2<500)) {
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
                         
                 } // muon
               } // #measurements
             } // selection
          } // index
      } // hscp

      int ntracks1=0;
      int ntracks20=0;
      for (int itr=0; itr<ntracks; itr++) {
//          cout << " debug loop tracks "  << itr << endl;
         int presk= 1;
         if (year!=2016) presk=track_prescale[itr];
         int index_of_the_track=itr;
         bool selection=true;
         if (track_pt[index_of_the_track]>50) selection=false;
         if (abs(track_eta[index_of_the_track])>2.1) selection=false;

         if (track_nvalidhits[index_of_the_track]<8) selection=false;
         if (track_npixhits[index_of_the_track]<2) selection=false;
         if (track_validfraction[index_of_the_track]<0.8) selection=false;
         if (track_missing[index_of_the_track]>99999) selection=false;
         if (track_validlast[index_of_the_track]<-99999) selection=false;

             
         bool is_high_qual =  (track_qual[index_of_the_track] & (1 << TrackQuality::highPurity)) >> TrackQuality::highPurity ;
         if (!is_high_qual) selection=false;
         if (track_chi2[index_of_the_track]>5) selection=false;

         // which selection for dz and dyx w/r to the primary vertex ? here 1st PV to pass the selection
         if (abs(track_dxy[index_of_the_track])>0.5) selection=false;
         if (abs(track_dz[index_of_the_track])>0.5) selection=false;

         // no cut on the isolation

         // cut on sigma pT/pT for signal  :
         if (track_pterr[index_of_the_track]/track_pt[index_of_the_track]>0.25) selection=false;
         if (!selection) continue;

         Htrackpt->Fill(track_pt[itr],presk);

         if (track_pt[itr]>1) ntracks1+=presk;
         if (track_pt[itr]>20) ntracks20+=presk;


         std::vector <float> charge_corr;
         std::vector <float> pathlength;
         std::vector <int> subdetId;
         std::vector <UInt_t> detId;
         std::vector <int> moduleGeometry;
         std::vector <bool> bool_cleaning;
         std::vector <bool> mustBeInside;

         std::vector <float> charge_corr1;
         std::vector <float> pathlength1;
         std::vector <int> subdetId1;
         std::vector <UInt_t> detId1;
         std::vector <int> moduleGeometry1;
         std::vector <bool> bool_cleaning1;
         std::vector <bool> mustBeInside1;

         std::vector <float> charge_corr2;
         std::vector <float> pathlength2;
         std::vector <int> subdetId2;
         std::vector <int> moduleGeometry2;
         std::vector <bool> bool_cleaning2;
         std::vector <bool> mustBeInside2;

         std::vector <float> charge_corr3;
         std::vector <float> pathlength3;
         std::vector <int> subdetId3;
         std::vector <UInt_t> detId3;
         std::vector <int> moduleGeometry3;
         std::vector <bool> bool_cleaning3;
         std::vector <bool> mustBeInside3;

         std::vector <float> charge_corr4;
         std::vector <float> pathlength4;
         std::vector <int> subdetId4;
         std::vector <UInt_t> detId4;
         std::vector <int> moduleGeometry4;
         std::vector <bool> bool_cleaning4;
         std::vector <bool> mustBeInside4;


         int nsatclu=0;
         for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
               float ch1=dedx_charge[iclu];
               bool clean1=true;
               if (dedx_subdetid[iclu]>=3) {
                        // strip
                        // without any correction :
//                        ch1=sclus_charge[iclu];  // charge without any correction
//                        clean1=sclus_clusclean[iclu];  
//                      // with Saturation only (but no Xtalk inversion) :
                          float check_charge=0;
                          vector<int> Quncor;
                          for (int istrip=sclus_index_strip[iclu]; istrip<sclus_index_strip[iclu]+sclus_nstrip[iclu]; istrip++) {
                            check_charge+=strip_ampl[istrip];
                            Quncor.push_back(strip_ampl[istrip]);
                          }
                          float deltaDif=check_charge-ch1;
                          if (deltaDif<0) deltaDif*=-1;
                          if (deltaDif>0.001) std::cout << "difference dans le calcul de la charge " << ch1 << " " << check_charge << " --> probleme acces ampl ???? " << std::endl; 
                          vector<int> Qcor = SaturationCorrection(Quncor,0.10,0.04, true,20,25);
                          float newcharge =0;
                          for (unsigned int inwc=0; inwc<Qcor.size(); inwc++) { newcharge+=Qcor[inwc]; }
                          ch1=newcharge;
                          clean1=sclus_clusclean2[iclu]; // clusterCleaning with Xtalk inversion and Saturation (September 22, 2021)
                          // for the record : there is a possibility to recompute the cluster cleaning with the Saturation only : clusterCleaning(Qcor, 1, &exitCode)
                          // but it is not what we want here
                          if (dedx_isstrip[iclu] && clean1 && dedx_insideTkMod[iclu] && (sclus_sat254[iclu] || sclus_sat255[iclu])) nsatclu++;
               }
               else {
                        // pixel
                          // float scaling =GetSFPixel(dedx_subdetid[iclu], dedx_detid[iclu], year, abs(track_eta[hscp_track_idx[ihs]]), runNumber);
                          if (PixelCorr2Apply) {
                            float scaling =GetSFPixelTamas(dedx_subdetid[iclu], dedx_detid[iclu], year, runNumber);
                            ch1*=scaling;
                          }
               }
               if (clean1 && dedx_insideTkMod[iclu]) {
                 // fill both Strip and Pixel
                 charge_corr.push_back(ch1);
                 pathlength.push_back(dedx_pathlength[iclu]);
                 subdetId.push_back(dedx_subdetid[iclu]);
                 detId.push_back(dedx_detid[iclu]);
                 moduleGeometry.push_back(dedx_modulgeom[iclu]);
                 mustBeInside.push_back(dedx_insideTkMod[iclu]);
                 bool_cleaning.push_back(clean1);

                 if (dedx_isstrip[iclu]) {
                 // fill only Strip
                     charge_corr1.push_back(ch1);
                     pathlength1.push_back(dedx_pathlength[iclu]);
                     subdetId1.push_back(dedx_subdetid[iclu]);
                     detId1.push_back(dedx_detid[iclu]);
                     moduleGeometry1.push_back(dedx_modulgeom[iclu]);
                     mustBeInside1.push_back(dedx_insideTkMod[iclu]);
                     bool_cleaning1.push_back(clean1);
                     if (track_p[itr]<5) Htrackdedx_strip_lowp->Fill(ch1,presk);
                     if (track_pt[itr]>55) Htrackdedx_strip->Fill(ch1, presk);
                 }
                 else {
                 // fill only Pixel
                     charge_corr2.push_back(ch1);
                     pathlength2.push_back(dedx_pathlength[iclu]);
                     subdetId2.push_back(dedx_subdetid[iclu]);
                     moduleGeometry2.push_back(dedx_modulgeom[iclu]);
                     mustBeInside2.push_back(dedx_insideTkMod[iclu]);
                     bool_cleaning2.push_back(clean1);
                     if (track_p[itr]<5) Htrackdedx_pix_lowp->Fill(ch1,presk);
                     if (track_pt[itr]>55) Htrackdedx_pix->Fill(ch1,presk);
                 }
                 bool no_in_L1_pixel =true;
                 int info_layr=GetLayerLabel(dedx_subdetid[iclu], dedx_detid[iclu],year);
                 if (info_layr==23) no_in_L1_pixel=false;
                 if (no_in_L1_pixel) {
                       // fill both Strip and Pixel but without L1 Pixel
                       charge_corr3.push_back(ch1);
                       pathlength3.push_back(dedx_pathlength[iclu]);
                       subdetId3.push_back(dedx_subdetid[iclu]);
                       detId3.push_back(dedx_detid[iclu]);
                       moduleGeometry3.push_back(dedx_modulgeom[iclu]);
                       mustBeInside3.push_back(dedx_insideTkMod[iclu]);
                       bool_cleaning3.push_back(clean1);
                       if (!dedx_isstrip[iclu]) {
                       // fill only Pixel without L1
                         charge_corr4.push_back(ch1);
                         pathlength4.push_back(dedx_pathlength[iclu]);
                         subdetId4.push_back(dedx_subdetid[iclu]);
                         detId4.push_back(dedx_detid[iclu]);
                         moduleGeometry4.push_back(dedx_modulgeom[iclu]);
                         mustBeInside4.push_back(dedx_insideTkMod[iclu]);
                         bool_cleaning4.push_back(clean1);
                       }
                 }


              float norm_mult = 265; // 247 or 265?
              double Norm = 3.61e-06*norm_mult;
              double scaleFactor = dEdxSF[0];

              if (dedx_ispixel[iclu]) {
                  Norm = 3.61e-06;
                  scaleFactor *=dEdxSF[1];
                  double ChargeOverPathlength = scaleFactor*Norm*ch1;
                  if (dedx_pathlength[iclu]>0) ChargeOverPathlength/=dedx_pathlength[iclu];
                  else ChargeOverPathlength=0;
                  HHitPix->Fill(ChargeOverPathlength, presk);
                  if (no_in_L1_pixel) {
                   if (fabs(track_eta[itr])<0.4){
                    HHitProfilePix->Fill(track_p[itr], ChargeOverPathlength, presk);
                    HHit2DPix->Fill(track_p[itr], ChargeOverPathlength, presk);
		   }
                   if (track_p[itr]>5 && track_p[itr]<45) {
                      Charge_Vs_Path_noL1->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                   }
                  }
                  if (track_p[itr]>5 && track_p[itr]<45) {
                      Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10*norm_mult),presk);
                     
                  }
              }
              else {
                  double ChargeOverPathlength = scaleFactor*Norm*ch1;
		  if (dedx_pathlength[iclu]>0) ChargeOverPathlength/=dedx_pathlength[iclu];
	  	  else ChargeOverPathlength=0;
		  HHitStrip->Fill(ChargeOverPathlength, presk);
		   if (fabs(track_eta[itr])<0.4 && sclus_clusclean2[iclu] && dedx_insideTkMod[iclu] ){
			    HHitProfileStrip->Fill(track_p[itr], ChargeOverPathlength, presk);
			    HHit2DStrip->Fill(track_p[itr], ChargeOverPathlength, presk);
		  }
                  if (track_p[itr]>5 && track_p[itr]<45) {
                   Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10),presk);
                   Charge_Vs_Path_noL1->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor*ch1/(dedx_pathlength[iclu]*10), presk);
                  }
	      }

            }  // end if cleaning & inside

			     
         }
  
         // APPLY Selection on the Number of Measurements :
         if (charge_corr3.size()>6) {
            // Pixel no L1;
            float norm_mult = 265; // 247 or 265?
            for (unsigned int jch=0;jch<charge_corr4.size();jch++) {
		   float Norm = 3.61e-06;
		   float scaleFactor =dEdxSF[0]*dEdxSF[1];
		   double ChargeOverPathlength = scaleFactor*Norm*charge_corr4[jch];
		   if (pathlength4[jch]>0) ChargeOverPathlength/=pathlength4[jch];
		   else ChargeOverPathlength=0;
                   if (fabs(track_eta[itr])<0.4) HHit2DPix_NoM->Fill(track_p[itr], ChargeOverPathlength, presk);
                   if (track_p[itr]>5 && track_p[itr]<45) {
                      Charge_Vs_Path_noL1_NoM->Fill (moduleGeometry4[jch], pathlength4[jch]*10, scaleFactor*charge_corr4[jch]/(pathlength4[jch]*10*norm_mult),presk);
                      Charge_Vs_Path_NoM->Fill (moduleGeometry4[jch], pathlength4[jch]*10, scaleFactor*charge_corr4[jch]/(pathlength4[jch]*10*norm_mult),presk);
                   }
            }
            // Strip
            for (unsigned int jch=0;jch<charge_corr1.size();jch++) {
                  double Norm = 3.61e-06*norm_mult;
                  double scaleFactor = dEdxSF[0];
                  double ChargeOverPathlength = scaleFactor*Norm*charge_corr1[jch];
		  if (pathlength1[jch]>0) ChargeOverPathlength/=pathlength1[jch];
	  	  else ChargeOverPathlength=0;
		  HHitStrip->Fill(ChargeOverPathlength, presk);
		  if (fabs(track_eta[itr])<0.4 && bool_cleaning1[jch] && mustBeInside1[jch] )
                             HHit2DStrip_NoM->Fill(track_p[itr], ChargeOverPathlength, presk);
                  if (track_p[itr]>5 && track_p[itr]<45) {
                    Charge_Vs_Path_NoM->Fill (moduleGeometry1[jch], pathlength1[jch]*10, scaleFactor*charge_corr1[jch]/(pathlength1[jch]*10),presk);
                    Charge_Vs_Path_noL1_NoM->Fill (moduleGeometry1[jch], pathlength1[jch]*10, scaleFactor*charge_corr1[jch]/(pathlength1[jch]*10),presk);
                  }
             }
         }
         // endAPPLY
                     

         int nv=0;
         int ns=0;
         double ih_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF,  NULL,2, 0.15,  nv, ns);
         double ih_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF,  NULL,2, 0.15,  nv, ns);
         double ih_pix = getdEdX(charge_corr2, pathlength2, subdetId2, moduleGeometry2, bool_cleaning2, mustBeInside2, dEdxSF,  NULL,2, 0.15,  nv, ns);

         double ih0_cor = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF,  NULL,2, 0.,  nv, ns);
         double ih0_strip = getdEdX(charge_corr1, pathlength1, subdetId1, moduleGeometry1, bool_cleaning1, mustBeInside1, dEdxSF,  NULL,2, 0.,  nv, ns);
         double ih0_noL1 = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, NULL,2, 0., nv, ns);  // <======= important one
         double ih0_pixnoL1 = getdEdX(charge_corr4, pathlength4, subdetId4, moduleGeometry4, bool_cleaning4, mustBeInside4, dEdxSF, NULL,2, 0., nv, ns);

         double ih_corHiDropNoL1 = getdEdX(charge_corr3,  pathlength3,  subdetId3,  moduleGeometry3,  bool_cleaning3,  mustBeInside3,  dEdxSF, NULL,2, 0., 0.15, nv, ns);
         // Ias
         int nval20_0=0;
         int nsat20_0=0;
         double ias_noL1 = getdEdX(charge_corr3, pathlength3, subdetId3, moduleGeometry3, bool_cleaning3, mustBeInside3, dEdxSF, dEdxTemplatesNoL1,2, 0., nval20_0, nsat20_0);
         double ias_all = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, dEdxSF, dEdxTemplatesAll,2, 0., nval20_0, nsat20_0);

         float mass_strip=-1;
         if (ih_strip - Cval_ldstrip>0.) mass_strip=sqrt((ih_strip - Cval_ldstrip)*track_p[itr]*track_p[itr]/Kval_ldstrip);

         float mass_ih=-1;
         if (ih_cor - Cval_ld>0.) mass_ih=sqrt((ih_cor - Cval_ld)*track_p[itr]*track_p[itr]/Kval_ld);

         float mass_ih0=-1;
         if (ih0_cor - Cval_all>0.) mass_ih0= sqrt((ih0_cor - Cval_all)*track_p[itr]*track_p[itr]/Kval_all);

         float mass_ih0noL1=-1;
         float mass_ih0noL1_2=-1;
         float mass_ih0noL1_3=-1;
         if (ih0_noL1 - Cval_nol1>0.) mass_ih0noL1= sqrt((ih0_noL1 - Cval_nol1)*track_p[itr]*track_p[itr]/Kval_nol1);
         if (ih0_noL1 - Cval_nol1_2>0.) mass_ih0noL1_2= sqrt((ih0_noL1 - Cval_nol1_2)*track_p[itr]*track_p[itr]/Kval_nol1_2);
         float val_nol1_3 = Cval_nol1_3 + (0.92*0.92*Kval_nol1_3)/(6.5*6.5) +Nval_nol1_3*log(6.5/0.92);
         if (track_p[itr]<3) {
           if (ih0_noL1 - val_nol1_3  >0. && computeSpecial) mass_ih0noL1_3= getMassSpecial(ih0_noL1,track_p[itr], Kval_nol1_3, Cval_nol1_3,Nval_nol1_3);
         }

         float mass_0strip=-1;
         if (ih0_strip - Cval_strip>0.) mass_0strip= sqrt((ih0_strip - Cval_strip)*track_p[itr]*track_p[itr]/Kval_strip);

         float mass_ihHDnoL1=-1;
         if (ih_corHiDropNoL1 - Cval_hdnol1>0.) mass_ihHDnoL1= sqrt((ih_corHiDropNoL1 - Cval_hdnol1)*track_p[itr]*track_p[itr]/Kval_hdnol1);

         float mass_pixel=-1;
         if (ih_pix - Cval_pix>0) mass_pixel= sqrt((ih_pix - Cval_pix)*track_p[itr]*track_p[itr]/Kval_pix);

         float mass_pixnL1=-1;
         if (ih0_pixnoL1 - Cval_pixnol1>0) mass_pixnL1= sqrt((ih0_pixnoL1 - Cval_pixnol1)*track_p[itr]*track_p[itr]/Kval_pixnol1);


         if (track_p[itr]<20) {
            dEdXVsP_lowp2->Fill(track_p[itr],ih_cor,presk);

            dEdX0VsP_lowp2->Fill(track_p[itr],ih0_cor,presk);
            dEdX0noL1VsP_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            // eta bins motivated by Dylan's presentation on March 12, 2021
            if (fabs(track_eta[itr])<0.91) {
                 dEdX0noL1VsP_eta1_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (fabs(track_eta[itr])<1.74) {
                 dEdX0noL1VsP_eta2_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_eta3_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            if (npv<20.) {
                 dEdX0noL1VsP_pu1_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<30.) {
                 dEdX0noL1VsP_pu2_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<40.) {
                 dEdX0noL1VsP_pu3_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_pu4_lowp2->Fill(track_p[itr],ih0_noL1,presk);
            }
         
            dEdXHDnoL1VsP_lowp2->Fill(track_p[itr],ih_corHiDropNoL1,presk);
            dEdXstripVsP_lowp2->Fill(track_p[itr],ih_strip,presk);
            dEdXpixVsP_lowp2->Fill(track_p[itr],ih_pix,presk);
            dEdX0stripVsP_lowp2->Fill(track_p[itr],ih0_strip,presk);
            dEdX0pixnoL1VsP_lowp2->Fill(track_p[itr],ih0_pixnoL1,presk);
         }
         if (track_p[itr]<5) {
            Htrackih_lowp->Fill(ih_cor, presk);
            Htrackih_pix_lowp->Fill(ih_pix,presk);
            Htrackih_strip_lowp->Fill(ih_strip,presk);
            Htrackih0_lowp->Fill(ih0_cor,presk);
            Htrackih0noL1_lowp->Fill(ih0_noL1,presk);
            Htracketa_lowp->Fill(track_eta[itr],presk);
            Htrackias_lowp->Fill(ias_noL1,presk);
            Htrackiasall_lowp->Fill(ias_all,presk);

            dEdXVsP_lowp->Fill(track_p[itr],ih_cor,presk);
            dEdXpixVsP_lowp->Fill(track_p[itr],ih_pix,presk);
            dEdXstripVsP_lowp->Fill(track_p[itr],ih_strip,presk);
            dEdX0stripVsP_lowp->Fill(track_p[itr],ih0_strip,presk);
            dEdX0VsP_lowp->Fill(track_p[itr],ih0_cor,presk);
            dEdX0noL1VsP_lowp->Fill(track_p[itr],ih0_noL1,presk);
            // eta bins motivated by Dylan's presentation on March 12, 2021
            if (fabs(track_eta[itr])<0.91) {
                 dEdX0noL1VsP_eta1_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (fabs(track_eta[itr])<1.74) {
                 dEdX0noL1VsP_eta2_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_eta3_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            if (npv<20.) {
                 dEdX0noL1VsP_pu1_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<30.) {
                 dEdX0noL1VsP_pu2_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else if (npv<40.) {
                 dEdX0noL1VsP_pu3_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            else {
                 dEdX0noL1VsP_pu4_lowp->Fill(track_p[itr],ih0_noL1,presk);
            }
            dEdXHDnoL1VsP_lowp->Fill(track_p[itr],ih_corHiDropNoL1,presk);
            dEdX0pixnoL1VsP_lowp->Fill(track_p[itr],ih0_pixnoL1,presk);

            dEstrVsdE_lowp->Fill(ih_cor,ih_strip,presk);
            dEdXstripVsEta_lowp->Fill(track_eta[itr],ih_strip,presk);

            // test proton Mass 

            if (track_p[itr]<3) {
                    lowp_MassIhstrip->Fill(mass_strip,presk);
                    lowp_MassIh->Fill(mass_ih,presk);
                    lowp_MassIh0->Fill(mass_ih0,presk);
                    lowp_MassIh0noL1->Fill(mass_ih0noL1,presk);
                    lowp_MassIh0noL1_2->Fill(mass_ih0noL1_2,presk);
                    lowp_MassIh0noL1_3->Fill(mass_ih0noL1_3,presk);
                       if (ih0_noL1> Cval_nol1) {
                       lowp_MassIh0noL1_11->Fill(mass_ih0noL1);
                       lowp_MassIh0noL1_12->Fill(mass_ih0noL1_2);
                       lowp_MassIh0noL1_13->Fill(mass_ih0noL1_3);
                       }
                       if (ih0_noL1>3.27 + 0.21) {
                       lowp_MassIh0noL1_1s1->Fill(mass_ih0noL1);
                       lowp_MassIh0noL1_1s2->Fill(mass_ih0noL1_2);
                       lowp_MassIh0noL1_1s3->Fill(mass_ih0noL1_3);
                         if (ih0_noL1>3.27 + 2*0.21) {
                         lowp_MassIh0noL1_2s1->Fill(mass_ih0noL1);
                         lowp_MassIh0noL1_2s2->Fill(mass_ih0noL1_2);
                         lowp_MassIh0noL1_2s3->Fill(mass_ih0noL1_3);
                           if (ih0_noL1>3.27 + 3*0.21) {
                           lowp_MassIh0noL1_3s1->Fill(mass_ih0noL1);
                           lowp_MassIh0noL1_3s2->Fill(mass_ih0noL1_2);
                           lowp_MassIh0noL1_3s3->Fill(mass_ih0noL1_3);
                           }
                         }
                       }
                    lowp_MassIh0strip->Fill(mass_0strip,presk);
                    lowp_MassIhHDnoL1->Fill(mass_ihHDnoL1,presk);
                    lowp2d_MassIhstrip->Fill(track_p[itr],mass_strip,presk);
                    lowp2d_MassIh->Fill(track_p[itr],mass_ih,presk);
                    lowp2d_MassIh0->Fill(track_p[itr],mass_ih0,presk);
                    lowp2d_MassIh0noL1->Fill(track_p[itr],mass_ih0noL1,presk);
                    lowp2d_MassIh0strip->Fill(track_p[itr],mass_0strip,presk);
                    lowp2d_MassIhHDnoL1->Fill(track_p[itr],mass_ihHDnoL1,presk);
                    lowp_dEdXpixVsstrip->Fill(ih_strip,ih_pix,presk);
                    lowp_dEdX0pixVsstrip->Fill(ih0_strip,ih0_pixnoL1,presk);
                    lowp2d_Mass_pix_strip15->Fill(mass_strip,mass_pixel,presk);
                    lowp2d_Mass_pix_strip0->Fill(mass_0strip,mass_pixnL1,presk);
                    if (mass_0strip>0 && mass_pixnL1>0 ) {
                      lowp_MassDiff_pix_strip0->Fill(mass_0strip - mass_pixnL1,presk);
                    }
                    if (mass_strip>0 && mass_pixel>0) {
                      lowp_MassDiff_pix_strip15->Fill(mass_strip - mass_pixel,presk);
                    }
              }

         } 
         else {
          // track_p>5
          //
          // test for the background
                bg_lowp2d_Mass_pix_strip0->Fill(mass_0strip,mass_pixnL1,presk);
                if (mass_0strip>0 && mass_pixnL1>0) {
                     float meanval=(mass_0strip+mass_pixnL1)/2;
                     float alpha= (mass_0strip-mass_pixnL1)/meanval;
                     bg_transf_Mass->Fill(meanval,alpha,presk);
                }
          // stability
                bg_dEdX0NoL1VsRun->Fill(runNumber,ih0_noL1,presk);  //  <==== our best estimate for Ih
                bg_dEdXVsIas->Fill(ias_noL1,ih0_noL1,presk);

          std::vector<double> vect;
          std::vector <int> subdetvec;
          std::vector <UInt_t> detvec;
//          for (int h=0; h<charge_corr1.size(); h++) {   //only strip
          for (unsigned int h=0; h<charge_corr3.size(); h++) {      //all cluster
            if(!bool_cleaning[h])continue;
            if(!mustBeInside[h])continue;
            double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
            double scaleFactor = dEdxSF[0];
            if (subdetId[h]<3) scaleFactor*=dEdxSF[1];
            double ChargeOverPathlength = 0;
            if (pathlength[h]>0) ChargeOverPathlength = scaleFactor*Norm*charge_corr[h]/pathlength[h];
            vect.push_back(ChargeOverPathlength); //save charge
            subdetvec.push_back(subdetId[h]);
            detvec.push_back(detId[h]);
            if (10<track_pt[itr] && track_pt[itr]<30) {
              if (subdetId[h]<3) {
                 LowpCharge_Eta_pix->Fill(track_eta[itr],ChargeOverPathlength,presk);
              }
              else { 
                 LowpCharge_Eta_strip->Fill(track_eta[itr],ChargeOverPathlength,presk);
              }
            }
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

   
              if (subdetvec[t]==2) {
                 int layerorside = int((detvec[t]>>23)&0x3); // 1=FPIX- 2=FPIX+
                 if (layerorside==1) {  
                     Charge_pixr1->Fill(tmp[t],presk);
                     ChargeVsRun_pixr1->Fill(runNumber,tmp[t],presk);
                     ZooChargeVsRun_pixr1->Fill(runNumber,tmp[t],presk);
                 }
                 else if (layerorside==2) {
                     Charge_pixr2->Fill(tmp[t],presk);
                     ChargeVsRun_pixr2->Fill(runNumber,tmp[t],presk);
                     ZooChargeVsRun_pixr2->Fill(runNumber,tmp[t],presk);
                 }
              }
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

            if (10<track_pt[itr] && track_pt[itr]<30) {
              if (lab==1)      LowpCharge_tib1->Fill(tmp[t],presk);
              else if (lab==2) LowpCharge_tib2->Fill(tmp[t],presk);
              else if (lab==3) LowpCharge_tib3->Fill(tmp[t],presk);
              else if (lab==4) LowpCharge_tib4->Fill(tmp[t],presk);
              else if (lab==5) LowpCharge_tob1->Fill(tmp[t],presk);
              else if (lab==6) LowpCharge_tob2->Fill(tmp[t],presk);
              else if (lab==7) LowpCharge_tob3->Fill(tmp[t],presk);
              else if (lab==8) LowpCharge_tob4->Fill(tmp[t],presk);
              else if (lab==9) LowpCharge_tob5->Fill(tmp[t],presk);
              else if (lab==10) LowpCharge_tob6->Fill(tmp[t],presk);
              else if (lab==11) LowpCharge_tid1->Fill(tmp[t],presk);
              else if (lab==12) LowpCharge_tid2->Fill(tmp[t],presk);
              else if (lab==13) LowpCharge_tid3->Fill(tmp[t],presk);
              else if (lab==14) LowpCharge_tec1->Fill(tmp[t],presk);
              else if (lab==15) LowpCharge_tec2->Fill(tmp[t],presk);
              else if (lab==16) LowpCharge_tec3->Fill(tmp[t],presk);
              else if (lab==17) LowpCharge_tec4->Fill(tmp[t],presk);
              else if (lab==18) LowpCharge_tec5->Fill(tmp[t],presk);
              else if (lab==19) LowpCharge_tec6->Fill(tmp[t],presk);
              else if (lab==20) LowpCharge_tec7->Fill(tmp[t],presk);
              else if (lab==21) LowpCharge_tec8->Fill(tmp[t],presk);
              else if (lab==22) LowpCharge_tec9->Fill(tmp[t],presk);
              else if (lab==23) LowpCharge_pixl1->Fill(tmp[t],presk);
              else if (lab==24) LowpCharge_pixl2->Fill(tmp[t],presk);
              else if (lab==25) LowpCharge_pixl3->Fill(tmp[t],presk);
              else if (lab==26) LowpCharge_pixl4->Fill(tmp[t],presk);
              else if (lab==27) LowpCharge_pixd1->Fill(tmp[t],presk);
              else if (lab==28) LowpCharge_pixd2->Fill(tmp[t],presk);
              else if (lab==29) LowpCharge_pixd3->Fill(tmp[t],presk);

   
              if (subdetvec[t]==2) {
                 int layerorside = int((detvec[t]>>23)&0x3); // 1=FPIX- 2=FPIX+
                 if (layerorside==1) {  
                     LowpCharge_pixr1->Fill(tmp[t],presk);
                 }
                 else if (layerorside==2) {
                     LowpCharge_pixr2->Fill(tmp[t],presk);
                 }
              }
            }



          } // end loop on tmp
//          dEdXstripVsNhittrunc_lowp->Fill(vect.size()-nTrunc,ih_strip,presk);
         } // end cut on p



         Htracketa->Fill(track_eta[itr],presk);
         Htrackphi->Fill(track_phi[itr],presk);
         dEdXVsP->Fill(track_p[itr],ih0_noL1,presk);
         dEdXpixVsP->Fill(track_p[itr],ih0_pixnoL1,presk);
         dEdXstripVsP->Fill(track_p[itr],ih0_strip,presk);
         EtaVsPhi_nhit->Fill(track_phi[itr],track_eta[itr],presk);

         Nsat->Fill(nsatclu,presk);
         NPix->Fill(charge_corr4.size(),presk);
         NStrip->Fill(charge_corr1.size(),presk);

         Htrackih_reco->Fill(ih0_noL1,presk);
         Htrackih_pix->Fill(ih0_pixnoL1,presk);
         Htrackih_strip->Fill(ih0_strip,presk);
         Htrackias->Fill(ias_noL1,presk);
         Htrackiasall->Fill(ias_all,presk);
      }
      HNtracks->Fill(ntracks);
      HNtracks1->Fill(ntracks1);
      HNtracks20->Fill(ntracks20);

   }
   OutputHisto->cd();
   HNtracks->Write();
   HNtracks1->Write();
   HNtracks20->Write();

   Htrackpt->Write();
   Htracketa->Write();
   Htracketa_lowp->Write();
   Htrackphi->Write();
   Htracknhit->Write();

   Htrackih_reco->Write();
   Htrackih_pix->Write();
   Htrackih_strip->Write();
   Htrackdedx_pix->Write();
   Htrackdedx_strip->Write();
   Htrackias->Write();
   Htrackiasall->Write();

   Htrackih_lowp->Write();
   Htrackih_pix_lowp->Write();
   Htrackih_strip_lowp->Write();
   Htrackih0_lowp->Write();
   Htrackih0noL1_lowp->Write();
   Htrackdedx_pix_lowp->Write();
   Htrackdedx_strip_lowp->Write();
   Htrackdedx_strip_lowp1->Write();
   Htrackdedx_strip_lowp2->Write();
   Htrackias_lowp->Write();
   Htrackiasall_lowp->Write();


   Nsat->Write();
   NPix->Write();
   NStrip->Write();
   dEdXVsP->Write();
   dEdXpixVsP->Write();
   dEdXstripVsP->Write();

   dEdXVsP_lowp->Write();
   dEdXVsP_lowp2->Write();

   dEdX0VsP_lowp->Write();
   dEdX0VsP_lowp2->Write();
   dEdX0noL1VsP_lowp->Write();
   dEdX0noL1VsP_lowp2->Write();
   dEdX0noL1VsP_eta1_lowp->Write();
   dEdX0noL1VsP_eta1_lowp2->Write();
   dEdX0noL1VsP_eta2_lowp->Write();
   dEdX0noL1VsP_eta2_lowp2->Write();
   dEdX0noL1VsP_eta3_lowp->Write();
   dEdX0noL1VsP_eta3_lowp2->Write();
   dEdX0noL1VsP_pu1_lowp->Write();
   dEdX0noL1VsP_pu1_lowp2->Write();
   dEdX0noL1VsP_pu2_lowp->Write();
   dEdX0noL1VsP_pu2_lowp2->Write();
   dEdX0noL1VsP_pu3_lowp->Write();
   dEdX0noL1VsP_pu3_lowp2->Write();
   dEdX0noL1VsP_pu4_lowp->Write();
   dEdX0noL1VsP_pu4_lowp2->Write();
   dEdXHDnoL1VsP_lowp->Write();
   dEdXHDnoL1VsP_lowp2->Write();
   dEdX0pixnoL1VsP_lowp->Write();
   dEdX0pixnoL1VsP_lowp2->Write();

   dEdXpixVsP_lowp->Write();
   dEdXpixVsP_lowp2->Write();
   dEdXstripVsP_lowp->Write();
   dEdXstripVsP_lowp2->Write();
   dEdX0stripVsP_lowp->Write();
   dEdX0stripVsP_lowp2->Write();

   dEdXstripVsEta_lowp->Write();
   dEstrVsdE_lowp->Write();
   dEdXstripVsNhit_lowp->Write();
   dEdXstripVsNhittrunc_lowp->Write();
   dEdXstripVsCharge_lowp->Write();
   EtaVsPhi_nhit->Write();

   Charge_pixl1->Write();
   Charge_pixl2->Write();
   Charge_pixl3->Write();
   Charge_pixl4->Write();
   Charge_pixd1->Write();
   Charge_pixd2->Write();
   Charge_pixd3->Write();
   Charge_pixr1->Write();
   Charge_pixr2->Write();
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

   LowpCharge_tib1->Write();
   LowpCharge_tib2->Write();
   LowpCharge_tib3->Write();
   LowpCharge_tib4->Write();
   LowpCharge_tob1->Write();
   LowpCharge_tob2->Write();
   LowpCharge_tob3->Write();
   LowpCharge_tob4->Write();
   LowpCharge_tob5->Write();
   LowpCharge_tob6->Write();
   LowpCharge_tid1->Write();
   LowpCharge_tid2->Write();
   LowpCharge_tid3->Write();
   LowpCharge_tec1->Write();
   LowpCharge_tec2->Write();
   LowpCharge_tec3->Write();
   LowpCharge_tec4->Write();
   LowpCharge_tec5->Write();
   LowpCharge_tec6->Write();
   LowpCharge_tec7->Write();
   LowpCharge_tec8->Write();
   LowpCharge_tec9->Write();
   LowpCharge_pixl1->Write();
   LowpCharge_pixl2->Write();
   LowpCharge_pixl3->Write();
   LowpCharge_pixl4->Write();
   LowpCharge_pixd1->Write();
   LowpCharge_pixd2->Write();
   LowpCharge_pixd3->Write();
   LowpCharge_pixr1->Write();
   LowpCharge_pixr2->Write();
   LowpCharge_Eta_pix->Write();
   LowpCharge_Eta_strip->Write();

   ChargeVsRun_pixl1->Write();
   ChargeVsRun_pixl2->Write();
   ChargeVsRun_pixl3->Write();
   ChargeVsRun_pixl4->Write();
   ChargeVsRun_pixd1->Write();
   ChargeVsRun_pixd2->Write();
   ChargeVsRun_pixd3->Write();
   ChargeVsRun_pixr1->Write();
   ChargeVsRun_pixr2->Write();
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
   ZooChargeVsRun_pixr1->Write();
   ZooChargeVsRun_pixr2->Write();
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
   dEdXNoL1pixVsRun->Write();
   dEdXNoL1VsRun->Write();
   dEdXHiDropVsRun->Write();
   dEdXpixHiDropVsRun->Write();
   dEdXstripHiDropVsRun->Write();
   dEdXHiDropNoL1VsRun->Write();
/*
   dEdXVsIL->Write();
   dEdXpixVsIL->Write();
   dEdXstripVsIL->Write();
*/
   dEdX0VsRun->Write();
   dEdX0pixVsRun->Write();
   dEdX0stripVsRun->Write();
   MassStripVsRun->Write();
   MassNoL1VsRun->Write();
   dEdX0NoL1pixVsRun->Write();
   dEdX0NoL1VsRun->Write();
   dEdX4VsRun->Write();
   dEdX4pixVsRun->Write();
   dEdX4stripVsRun->Write();
   dEdX40VsRun->Write();
   dEdX40pixVsRun->Write();
   dEdX40stripVsRun->Write();
   bg_dEdX0NoL1VsRun->Write();
   iasNoL1VsRun->Write();
   iasAllVsRun->Write();
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
   HSCP_MassIh->Write(); 
   HSCP_MassIh0->Write(); 
   HSCP_MassIhstrip->Write(); 
   HSCP_MassIh0noL1->Write(); 
   HSCP_MassIh0noL1_2->Write(); 
   HSCP_MassIh0noL1_3->Write(); 
   HSCP_MassIh0noL1_11->Write(); 
   HSCP_MassIh0noL1_12->Write(); 
   HSCP_MassIh0noL1_13->Write(); 
   HSCP_MassIh0noL1_1s1->Write(); 
   HSCP_MassIh0noL1_1s2->Write(); 
   HSCP_MassIh0noL1_1s3->Write(); 
   HSCP_MassIh0noL1_2s1->Write(); 
   HSCP_MassIh0noL1_2s2->Write(); 
   HSCP_MassIh0noL1_2s3->Write(); 
   HSCP_MassIh0noL1_3s1->Write(); 
   HSCP_MassIh0noL1_3s2->Write(); 
   HSCP_MassIh0noL1_3s3->Write(); 
   HSCP_MassIhHDnoL1->Write(); 
   HSCP_MassTOF->Write(); 
   HSCP2d_MassTOFvsIh->Write(); 
   HSCP_MassIh0strip->Write(); 
   HSCP2d_MassIh->Write(); 
   HSCP2d_MassIh0->Write(); 
   HSCP2d_MassIhstrip->Write(); 
   HSCP2d_MassIh0noL1->Write(); 
   HSCP2d_MassIhHDnoL1->Write(); 
   HSCP2d_MassIh0strip->Write(); 
   HSCP2d_Mass_pix_strip15->Write(); 
   HSCP2d_Mass_pix_strip0->Write(); 
   HSCP2d_Mass_pix_strip->Write(); 
   HSCP_MassDiff_pix_strip0->Write(); 
   HSCP_MassDiff_pix_strip15->Write(); 
   HSCP_MassResol_pix_strip0->Write(); 
   HSCP_MassResol_pix_strip15->Write(); 
   lowp_MassIh->Write(); 
   lowp_MassIh0->Write(); 
   lowp_MassIhstrip->Write(); 
   lowp_MassIh0noL1->Write(); 
   lowp_MassIh0noL1_2->Write(); 
   lowp_MassIh0noL1_3->Write(); 
   lowp_MassIh0noL1_11->Write(); 
   lowp_MassIh0noL1_12->Write(); 
   lowp_MassIh0noL1_13->Write(); 
   lowp_MassIh0noL1_1s1->Write(); 
   lowp_MassIh0noL1_1s2->Write(); 
   lowp_MassIh0noL1_1s3->Write(); 
   lowp_MassIh0noL1_2s1->Write(); 
   lowp_MassIh0noL1_2s2->Write(); 
   lowp_MassIh0noL1_2s3->Write(); 
   lowp_MassIh0noL1_3s1->Write(); 
   lowp_MassIh0noL1_3s2->Write(); 
   lowp_MassIh0noL1_3s3->Write(); 
   lowp_MassIhHDnoL1->Write(); 
   lowp_MassIh0strip->Write(); 
   lowp2d_MassIh->Write(); 
   lowp2d_MassIh0->Write(); 
   lowp2d_MassIhstrip->Write(); 
   lowp2d_MassIh0noL1->Write(); 
   lowp2d_MassIhHDnoL1->Write(); 
   lowp2d_MassIh0strip->Write(); 
   lowp_dEdXpixVsstrip->Write();
   lowp_dEdX0pixVsstrip->Write();
   lowp2d_Mass_pix_strip15->Write(); 
   lowp2d_Mass_pix_strip0->Write(); 
   bg_lowp2d_Mass_pix_strip0->Write(); 
   bg_transf_Mass->Write(); 
   lowp_MassDiff_pix_strip0->Write(); 
   lowp_MassDiff_pix_strip15->Write(); 
   bg_dEdXVsIas->Write(); 

   HSCP_dEdXpixVsstrip->Write();
   HSCP_dEdX0pixVsstrip->Write();
   HSCP_dEdXstripVsall->Write();
   HSCP_dEdXpixVsall->Write();
   HSCP_dEdXHiDrop->Write();
   HSCP_dEdXstripHiDrop->Write();
   HSCP_dEdXpixHiDrop->Write();
   HSCP_dEdXHiDropNoL1->Write();
   HSCP_dEdX0NoL1->Write();

   HSCP_FMIP4->Write(); 
   HSCP_FMIP3p5->Write(); 
   HSCP_FMIP3p2->Write(); 
   FMIP4VsRun->Write();
   FMIP3p5VsRun->Write();
   FMIP3p5VsRun->Write();
   FMIP4VsEta->Write();
   HSCP_iasnol1->Write(); 
   HSCP_iasall->Write(); 
   HSCP_probQ->Write();
   HSCP_probQNoL1->Write();
   HSCP_probXY->Write();
   HSCP_probXYNoL1->Write();
   probQVsRun->Write();
   probQNoL1VsRun->Write();
   probXYVsRun->Write();
   probXYNoL1VsRun->Write();
   probQVsIas->Write();

   HSCP_pt->Write(); 
   HSCP_eta->Write(); 
   HSCP_iso_eop->Write(); 
   nPV->Write(); 
   HSCP_invB->Write(); 
   HSCP_errinvB->Write(); 
   HSCP_invBDT->Write(); 
   HSCP_invBCSC->Write(); 
   HSCP_time->Write(); 
   HSCP_npix->Write(); 
   HSCP_nstrip->Write(); 
   HSCP_nmpix->Write(); 
   HSCP_nmstrip->Write(); 
   HSCP_nratio->Write(); 
   HSCP_nmratio->Write(); 

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


   HHitPix->Write();
   HHitProfilePix->Write();
   HHit2DPix->Write();
   HHit2DPix_NoM->Write();
   HHitStrip->Write();
   HHitProfileStrip->Write();
   HHit2DStrip->Write();
   HHit2DStrip_NoM->Write();

   OutputHisto->Close();

   if (writeTemplateOnDisk) {
   OutputTemplate->cd();
   Charge_Vs_Path->Write();
   Charge_Vs_Path_noL1->Write();
   Charge_Vs_Path_NoM->Write();
   Charge_Vs_Path_noL1_NoM->Write();
   OutputTemplate->Close();
   }

}


double run2analysis::getMassSpecial(float ih, float p, float K, float C, float N){

 float m1=sqrt((ih-C)*p*p/K);

 char fitfunc[2048];
 sprintf(fitfunc,"%6.4f*pow(x/%6.4f,2) + %6.4f + %6.4f*log(%6.4f/x) -%6.4f",K,p,C,N,p,ih);

 TF1 *fa1 = new TF1("fa1",fitfunc,0.1,m1*3.);

 float evalval = fa1->Eval(m1);
 float mbest=0;
 float ebest=1000;

 for (int i=0;i<2000;i++) {
   if (evalval<0) {
     float m2=m1+i*0.0001*m1;
     float evalval2 = fa1->Eval(m2);
     if (fabs(evalval2)<0.001) {
          if (fabs(evalval2)<ebest) {
                      mbest=m2;
                      ebest=fabs(evalval2);
          }
     }
     else if (evalval2>0.001) i+=2000; // this should allow to stop the search of solution
   }
   else {
     float m2=m1-i*0.0001*m1;
     if (m2>0) {
           float evalval2 = fa1->Eval(m2);
           if (fabs(evalval2)<0.001)  {
                   if (fabs(evalval2)<ebest) {
                         mbest=m2;
                         ebest=fabs(evalval2);
                    }
           }
           else if (evalval2<-0.001) i+=2000; // this should allow to stop the search of solution
      }
   }
  }

  return mbest;
}

double run2analysis::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, int & nv, int & ns) {
  double result= getdEdX(charge, pathlength, subdetId, moduleGeometry, bool_cleaning, mustBeInside, scaleFactors, templateHisto, n_estim, dropLowerDeDxValue, 0., nv, ns);
  return result;
}


double run2analysis::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, double dropHigherDeDxValue, int & nv, int & ns) {
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

     if(dropHigherDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::less<double>() );
         int nTrunc = tmp.size()*dropHigherDeDxValue;

         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }
     if (debugprint) std::cout << " after dropHigherDeDxValue " << std::endl;






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

//float run2analysis::FMIP(const vector<float>& charge, const<vector>& path, float thre = 4){
float run2analysis::FMIP(std::vector <float> charge, std::vector <float> path, float SFactor, float thre = 4){
   if(charge.size()!=path.size()) return -999; // error

   int nclusters = charge.size();
   if(nclusters==0) return -888;

   int nlow = 0; 
   for(unsigned int i=0;i<charge.size();i++){
       //hard-coded conversion factor to go for dEdx in MeV.cm2/g
       // ATTENTION : HERE THE CONVERSION FACTOR IS 247 WHILE ABOVE WE USE 265 !!!!!
       float dEdx = SFactor*charge[i]*(3.61*pow(10,-9)*247)*1000/path[i];
       if(dEdx<thre) nlow++;
   }
   return nlow*1./nclusters;
}


int run2analysis::GetLayerLabel(int subdetid_, UInt_t detid_, int year)
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


void run2analysis::loadSFPixelCalib() {

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

float run2analysis::GetSFPixel(int subdetid_, UInt_t detid_, int year, float eta, int run) {
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

void run2analysis::loadSFPixelTamas() {

   std::ifstream file_calib1;
   file_calib1.open("Tamas/CorrFactHistory2017L1.txt");
   icalibL1=0;
   if (!file_calib1) { std::cerr << "cannot open file CorrFactHistory2017L1.txt " << std::endl; }
   else
   {
      while (!file_calib1.eof () && icalibL1<calmax) {
       // from_run_number, corr_val, error_val
       file_calib1 >> irunMinValL1[icalibL1] >> scaleValL1[icalibL1] >> errorScaleValL1[icalibL1] ;
//       if (icalibL1<10) { std::cout << irunMinValL1[icalibL1] << " " <<  scaleValL1[icalibL1] << " " <<  errorScaleValL1[icalibL1] << std::endl ; }
       icalibL1++;
      }
   }
   std::cout << " file_calib1 : " << icalibL1 << " entries in L1 2017" << std::endl;
//   if (icalibL1>1) std::cout << " example runMin "<< irunMinValL1[0] << " runMax "<< irunMinValL1[1] << " scale " << scaleValL1[0] << " +/- " << errorScaleValL1[0] << std::endl;
   file_calib1.close ();

   icalibL1_2017=icalibL1;
   std::ifstream file_calib2;
   file_calib2.open("Tamas/CorrFactHistory2018L1.txt");
   if (!file_calib2) { std::cerr << "cannot open file CorrFactHistory2018L1.txt " << std::endl; }
   else
   {  
      while (!file_calib2.eof () && icalibL1<calmax) {
       file_calib2 >> irunMinValL1[icalibL1] >> scaleValL1[icalibL1] >> errorScaleValL1[icalibL1] ;
       icalibL1++;
      }
   }
   if (irunMinValL1[icalibL1-1]==0) icalibL1--;
   icalibL1_2018=icalibL1-icalibL1_2017;
   std::cout << " file_calib : " << icalibL1_2018 << " entries in L1 2018" << std::endl;
//   if (icalibL1>icalibL1_2017+1) std::cout << " example runMin "<< irunMinValL1[icalibL1_2017] << " runMax "<< irunMinValL1[icalibL1_2017+1] 
//                                           << " scale " << scaleValL1[icalibL1_2017] << " +/- " << errorScaleValL1[icalibL1_2017] << std::endl;
   file_calib2.close ();
 

   std::ifstream file_calib3;
   file_calib3.open("Tamas/CorrFactHistory2017L2.txt");
   icalibL2=0;
   if (!file_calib3) { std::cerr << "cannot open file CorrFactHistory2017L2.txt " << std::endl; }
   else
   {
      while (!file_calib3.eof () && icalibL2<calmax) {
       // from_run_number, corr_val, error_val
       file_calib3 >> irunMinValL2[icalibL2] >> scaleValL2[icalibL2] >> errorScaleValL2[icalibL2] ;
//       if (icalibL2<10) { std::cout << irunMinValL2[icalibL2] << " " <<  scaleValL2[icalibL2] << " " <<  errorScaleValL2[icalibL2] << std::endl ; }
       icalibL2++;
      }
   }
//   std::cout << " file_calib3 : " << icalibL2 << " entries in L2 2017" << std::endl;
//   if (icalibL2>1) std::cout << " example runMin "<< irunMinValL2[0] << " runMax "<< irunMinValL2[1] << " scale " << scaleValL2[0] << " +/- " << errorScaleValL2[0] << std::endl;
   file_calib3.close ();

   icalibL2_2017=icalibL2;
   std::ifstream file_calib4;
   file_calib4.open("Tamas/CorrFactHistory2018L2.txt");
   if (!file_calib4) { std::cerr << "cannot open file CorrFactHistory2018L2.txt " << std::endl; }
   else
   {  
      while (!file_calib4.eof () && icalibL2<calmax) {
       file_calib4 >> irunMinValL2[icalibL2] >> scaleValL2[icalibL2] >> errorScaleValL2[icalibL2] ;
       icalibL2++;
      }
   }
   if (irunMinValL2[icalibL2-1]==0) icalibL2--;
   icalibL2_2018=icalibL2-icalibL2_2017;
//   std::cout << " file_calib : " << icalibL2_2018 << " entries in L2 2018" << std::endl;
//   if (icalibL2>icalibL2_2017+1) std::cout << " example runMin "<< irunMinValL2[icalibL2_2017] << " runMax "<< irunMinValL2[icalibL2_2017+1] 
//                                           << " scale " << scaleValL2[icalibL2_2017] << " +/- " << errorScaleValL2[icalibL2_2017] << std::endl;
   file_calib4.close ();
 
   std::ifstream file_calib5;
   file_calib5.open("Tamas/CorrFactHistory2017L3.txt");
   icalibL3=0;
   if (!file_calib5) { std::cerr << "cannot open file CorrFactHistory2017L3.txt " << std::endl; }
   else
   {
      while (!file_calib5.eof () && icalibL3<calmax) {
       // from_run_number, corr_val, error_val
       file_calib5 >> irunMinValL3[icalibL3] >> scaleValL3[icalibL3] >> errorScaleValL3[icalibL3] ;
//       if (icalibL3<10) { std::cout << irunMinValL3[icalibL3] << " " <<  scaleValL3[icalibL3] << " " <<  errorScaleValL3[icalibL3] << std::endl ; }
       icalibL3++;
      }
   }
//   std::cout << " file_calib5 : " << icalibL3 << " entries in L3 2017" << std::endl;
//   if (icalibL3>1) std::cout << " example runMin "<< irunMinValL3[0] << " runMax "<< irunMinValL3[1] << " scale " << scaleValL3[0] << " +/- " << errorScaleValL3[0] << std::endl;
   file_calib5.close ();

   icalibL3_2017=icalibL3;
   std::ifstream file_calib6;
   file_calib6.open("Tamas/CorrFactHistory2018L3.txt");
   if (!file_calib6) { std::cerr << "cannot open file CorrFactHistory2018L3.txt " << std::endl; }
   else
   {  
      while (!file_calib6.eof () && icalibL3<calmax) {
       file_calib6 >> irunMinValL3[icalibL3] >> scaleValL3[icalibL3] >> errorScaleValL3[icalibL3] ;
       icalibL3++;
      }
   }
   if (irunMinValL3[icalibL3-1]==0) icalibL3--;
   icalibL3_2018=icalibL3-icalibL3_2017;
//   std::cout << " file_calib : " << icalibL3_2018 << " entries in L3 2018" << std::endl;
//   if (icalibL3>icalibL3_2017+1) std::cout << " example runMin "<< irunMinValL3[icalibL3_2017] << " runMax "<< irunMinValL3[icalibL3_2017+1] 
//                                           << " scale " << scaleValL3[icalibL3_2017] << " +/- " << errorScaleValL3[icalibL3_2017] << std::endl;
   file_calib6.close ();
 
   std::ifstream file_calib7;
   file_calib7.open("Tamas/CorrFactHistory2017L4.txt");
   icalibL4=0;
   if (!file_calib7) { std::cerr << "cannot open file CorrFactHistory2017L4.txt " << std::endl; }
   else
   {
      while (!file_calib7.eof () && icalibL4<calmax) {
       // from_run_number, corr_val, error_val
       file_calib7 >> irunMinValL4[icalibL4] >> scaleValL4[icalibL4] >> errorScaleValL4[icalibL4] ;
//       if (icalibL4<10) { std::cout << irunMinValL4[icalibL4] << " " <<  scaleValL4[icalibL4] << " " <<  errorScaleValL4[icalibL4] << std::endl ; }
       icalibL4++;
      }
   }
//   std::cout << " file_calib7 : " << icalibL4 << " entries in L4 2017" << std::endl;
//   if (icalibL4>1) std::cout << " example runMin "<< irunMinValL4[0] << " runMax "<< irunMinValL4[1] << " scale " << scaleValL4[0] << " +/- " << errorScaleValL4[0] << std::endl;
   file_calib7.close ();

   icalibL4_2017=icalibL4;
   std::ifstream file_calib8;
   file_calib8.open("Tamas/CorrFactHistory2018L4.txt");
   if (!file_calib8) { std::cerr << "cannot open file CorrFactHistory2018L4.txt " << std::endl; }
   else
   {  
      while (!file_calib8.eof () && icalibL4<calmax) {
       file_calib8 >> irunMinValL4[icalibL4] >> scaleValL4[icalibL4] >> errorScaleValL4[icalibL4] ;
       icalibL4++;
      }
   }
   if (irunMinValL4[icalibL4-1]==0) icalibL4--;
   icalibL4_2018=icalibL4-icalibL4_2017;
//   std::cout << " file_calib : " << icalibL4_2018 << " entries in L4 2018" << std::endl;
//   if (icalibL4>icalibL4_2017+1) std::cout << " example runMin "<< irunMinValL4[icalibL4_2017] << " runMax "<< irunMinValL4[icalibL4_2017+1] 
//                                           << " scale " << scaleValL4[icalibL4_2017] << " +/- " << errorScaleValL4[icalibL4_2017] << std::endl;
   file_calib8.close ();
 
   std::ifstream file_calib9;
   file_calib9.open("Tamas/CorrFactHistory2017R1.txt");
   icalibR1=0;
   if (!file_calib9) { std::cerr << "cannot open file CorrFactHistory2017R1.txt " << std::endl; }
   else
   {
      while (!file_calib9.eof () && icalibR1<calmax) {
       // from_run_number, corr_val, error_val
       file_calib9 >> irunMinValR1[icalibR1] >> scaleValR1[icalibR1] >> errorScaleValR1[icalibR1] ;
//       if (icalibR1<10) { std::cout << irunMinValR1[icalibR1] << " " <<  scaleValR1[icalibR1] << " " <<  errorScaleValR1[icalibR1] << std::endl ; }
       icalibR1++;
      }
   }
//   std::cout << " file_calib9 : " << icalibR1 << " entries in R1 2017" << std::endl;
//   if (icalibR1>1) std::cout << " example runMin "<< irunMinValR1[0] << " runMax "<< irunMinValR1[1] << " scale " << scaleValR1[0] << " +/- " << errorScaleValR1[0] << std::endl;
   file_calib9.close ();

   icalibR1_2017=icalibR1;
   std::ifstream file_calib10;
   file_calib10.open("Tamas/CorrFactHistory2018R1.txt");
   if (!file_calib10) { std::cerr << "cannot open file CorrFactHistory2018R1.txt " << std::endl; }
   else
   {  
      while (!file_calib10.eof () && icalibR1<calmax) {
       file_calib10 >> irunMinValR1[icalibR1] >> scaleValR1[icalibR1] >> errorScaleValR1[icalibR1] ;
       icalibR1++;
      }
   }
   if (irunMinValR1[icalibR1-1]==0) icalibR1--;
   icalibR1_2018=icalibR1-icalibR1_2017;
//   std::cout << " file_calib : " << icalibR1_2018 << " entries in R1 2018" << std::endl;
//   if (icalibR1>icalibR1_2017+1) std::cout << " example runMin "<< irunMinValR1[icalibR1_2017] << " runMax "<< irunMinValR1[icalibR1_2017+1] 
//                                           << " scale " << scaleValR1[icalibR1_2017] << " +/- " << errorScaleValR1[icalibR1_2017] << std::endl;
   file_calib10.close ();
 
   std::ifstream file_calib11;
   file_calib11.open("Tamas/CorrFactHistory2017R2.txt");
   icalibR2=0;
   if (!file_calib11) { std::cerr << "cannot open file CorrFactHistory2017R2.txt " << std::endl; }
   else
   {
      while (!file_calib11.eof () && icalibR2<calmax) {
       // from_run_number, corr_val, error_val
       file_calib11 >> irunMinValR2[icalibR2] >> scaleValR2[icalibR2] >> errorScaleValR2[icalibR2] ;
//       if (icalibR2<10) { std::cout << irunMinValR2[icalibR2] << " " <<  scaleValR2[icalibR2] << " " <<  errorScaleValR2[icalibR2] << std::endl ; }
       icalibR2++;
      }
   }
//   std::cout << " file_calib11 : " << icalibR2 << " entries in R2 2017" << std::endl;
//   if (icalibR2>1) std::cout << " example runMin "<< irunMinValR2[0] << " runMax "<< irunMinValR2[1] << " scale " << scaleValR2[0] << " +/- " << errorScaleValR2[0] << std::endl;
   file_calib11.close ();

   if (irunMinValR2[icalibR2-1]==0) icalibR2--;
   icalibR2_2017=icalibR2;
   std::ifstream file_calib12;
   file_calib12.open("Tamas/CorrFactHistory2018R2.txt");
   if (!file_calib12) { std::cerr << "cannot open file CorrFactHistory2018R2.txt " << std::endl; }
   else
   {  
      while (!file_calib12.eof () && icalibR2<calmax) {
       file_calib12 >> irunMinValR2[icalibR2] >> scaleValR2[icalibR2] >> errorScaleValR2[icalibR2] ;
       icalibR2++;
      }
   }
   icalibR2_2018=icalibR2-icalibR2_2017;
//   std::cout << " file_calib : " << icalibR2_2018 << " entries in R2 2018" << std::endl;
//   if (icalibR2>icalibR2_2017+1) std::cout << " example runMin "<< irunMinValR2[icalibR2_2017] << " runMax "<< irunMinValR2[icalibR2_2017+1] 
//                                           << " scale " << scaleValR2[icalibR2_2017] << " +/- " << errorScaleValR2[icalibR2_2017] << std::endl;
   file_calib12.close ();
 


}
float run2analysis::GetSFPixelTamas(int subdetid_, UInt_t detid_, int year, int run) {

   int pix = 0;
   int layerorside = 0;

   if (subdetid_==1) {
        pix =1 ;
        if (year==2016) {
                layerorside = int((detid_>>16)&0xF);
        }
        else {
                layerorside = int((detid_>>20)&0xF);
        }
   }
   if (subdetid_==2) {
        pix =2 ;
        layerorside = int((detid_>>23)&0x3); // 1=FPIX- 2=FPIX+
   }

   float scale =1 ;

   if (pix==1) {
     if (layerorside==1) {

       if (year==2017) {
        for (int i=0; i<icalibL1_2017; i++) {
         if (run >= irunMinValL1[i] && run<irunMinValL1[i+1]) {
                scale = scaleValL1[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL1_2017; i<icalibL1-1; i++) {
         if (run >= irunMinValL1[i] && run<irunMinValL1[i+1]) {
                scale = scaleValL1[i];
                return scale;
         }
        }
        if (run >= irunMinValL1[icalibL1-1]) {
                scale = scaleValL1[icalibL1-1];
                return scale;
        }
       }

     }
     else if (layerorside==2) {

       if (year==2017) {
        for (int i=0; i<icalibL2_2017; i++) {
         if (run >= irunMinValL2[i] && run<irunMinValL2[i+1]) {
                scale = scaleValL2[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL2_2017; i<icalibL2-1; i++) {
         if (run >= irunMinValL2[i] && run<irunMinValL2[i+1]) {
                scale = scaleValL2[i];
                return scale;
         }
        }
        if (run >= irunMinValL2[icalibL2-1]) {
                scale = scaleValL2[icalibL2-1];
                return scale;
        }
       }

     }
     else if (layerorside==3) {
     }

       if (year==2017) {
        for (int i=0; i<icalibL3_2017; i++) {
         if (run >= irunMinValL3[i] && run<irunMinValL3[i+1]) {
                scale = scaleValL3[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL3_2017; i<icalibL3-1; i++) {
         if (run >= irunMinValL3[i] && run<irunMinValL3[i+1]) {
                scale = scaleValL3[i];
                return scale;
         }
        }
        if (run >= irunMinValL3[icalibL3-1]) {
                scale = scaleValL3[icalibL3-1];
                return scale;
        }
       }

     else if (layerorside==4) {

       if (year==2017) {
        for (int i=0; i<icalibL4_2017; i++) {
         if (run >= irunMinValL4[i] && run<irunMinValL4[i+1]) {
                scale = scaleValL4[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibL4_2017; i<icalibL4-1; i++) {
         if (run >= irunMinValL4[i] && run<irunMinValL4[i+1]) {
                scale = scaleValL4[i];
                return scale;
         }
        }
        if (run >= irunMinValL4[icalibL4-1]) {
                scale = scaleValL4[icalibL4-1];
                return scale;
        }
       }

     }
     else { 
        std::cout << "unknowm layer number in Barrel Pixel " << layerorside << std::endl; 
        return 0; 
     } 
     
   }
   else if (pix==2) {
     if (layerorside==1) {

       if (year==2017) {
        for (int i=0; i<icalibR1_2017; i++) {
         if (run >= irunMinValR1[i] && run<irunMinValR1[i+1]) {
                scale = scaleValR1[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibR1_2017; i<icalibR1-1; i++) {
         if (run >= irunMinValR1[i] && run<irunMinValR1[i+1]) {
                scale = scaleValR1[i];
                return scale;
         }
        }
        if (run >= irunMinValR1[icalibR1-1]) {
                scale = scaleValR1[icalibR1-1];
                return scale;
        }
       }

     }
     else if (layerorside==2) {

       if (year==2017) {
        for (int i=0; i<icalibR2_2017; i++) {
         if (run >= irunMinValR2[i] && run<irunMinValR2[i+1]) {
                scale = scaleValR2[i];
                return scale;
         }
        }
       }
       else if (year==2018) {
        for (int i=icalibR2_2017; i<icalibR2-1; i++) {
         if (run >= irunMinValR2[i] && run<irunMinValR2[i+1]) {
                scale = scaleValR2[i];
                return scale;
         }
        }
        if (run >= irunMinValR2[icalibR2-1]) {
                scale = scaleValR2[icalibR2-1];
                return scale;
        }
       }

     }
     else { 
        std::cout << "unknowm ring number in EndCap Pixel " << layerorside << std::endl; 
        return 0; 
     } 
   }

   return scale;
   
}




///
// Code from Dylan

std::vector<int> run2analysis::SaturationCorrection(const std::vector<int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat) {
  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);
  Double_t a=1-2*x1-2*x2;
//  TMatrix A(N,N);

//---  que pour 1 max bien net
 if(Q.size()<2 || Q.size()>8){
        for (unsigned int i=0;i<Q.size();i++){
                QII.push_back((int) Q[i]);
        }
        return QII;
  }
 if(way){
          vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())      ;
          if(*mQ>253){
                 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
                 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
                     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
          }
      else{
          return Q; // no saturation --> no x-talk inversion
      }
  }
//---
 // do nothing else
 return Q; 
}

///

bool run2analysis::clusterCleaning(const std::vector<int>&  Q, int crosstalkInv, uint8_t * exitCode)
{
     vector<int>  ampls = Q;
//   if(crosstalkInv==1)ampls = CrossTalkInv(ampls,0.10,0.04, true,20,25);


  // ----------------  COMPTAGE DU NOMBRE DE MAXIMA   --------------------------
  //----------------------------------------------------------------------------
  //
         Int_t NofMax=0; Int_t recur255=1; Int_t recur254=1;
         bool MaxOnStart=false;bool MaxInMiddle=false, MaxOnEnd =false;
         Int_t MaxPos=0;
        // D?but avec max
        if(ampls.size()!=1 && ((ampls[0]>ampls[1])
            || (ampls.size()>2 && ampls[0]==ampls[1] && ampls[1]>ampls[2] && ampls[0]!=254 && ampls[0]!=255)
            || (ampls.size()==2 && ampls[0]==ampls[1] && ampls[0]!=254 && ampls[0]!=255)) ){
          NofMax=NofMax+1;  MaxOnStart=true;  }

        // Maximum entour?
        if(ampls.size()>2){
          for (unsigned int i =1; i < ampls.size()-1; i++) {
                if( (ampls[i]>ampls[i-1] && ampls[i]>ampls[i+1])
                    || (ampls.size()>3 && i>0 && i<ampls.size()-2 && ampls[i]==ampls[i+1] && ampls[i]>ampls[i-1] && ampls[i]>ampls[i+2] && ampls[i]!=254 && ampls[i]!=255) ){
                 NofMax=NofMax+1; MaxInMiddle=true;  MaxPos=i;
                }
                if(ampls[i]==255 && ampls[i]==ampls[i-1]) {
                        recur255=recur255+1;
                        MaxPos=i-(recur255/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
                if(ampls[i]==254 && ampls[i]==ampls[i-1]) {
                        recur254=recur254+1;
                        MaxPos=i-(recur254/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
            }
        }
        // Fin avec un max
       if(ampls.size()>1){
          if(ampls[ampls.size()-1]>ampls[ampls.size()-2]
             || (ampls.size()>2 && ampls[ampls.size()-1]==ampls[ampls.size()-2] && ampls[ampls.size()-2]>ampls[ampls.size()-3] )
             ||  ampls[ampls.size()-1]==255){
           NofMax=NofMax+1;  MaxOnEnd=true;   }
         }
        // Si une seule strip touch?e
        if(ampls.size()==1){    NofMax=1;}

  // ---  SELECTION EN FONCTION DE LA FORME POUR LES MAXIMA UNIQUES ---------
  //------------------------------------------------------------------------
  //
  //               ____
  //              |    |____
  //          ____|    |    |
  //         |    |    |    |____
  //     ____|    |    |    |    |
  //    |    |    |    |    |    |____
  //  __|____|____|____|____|____|____|__
  //    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
  //

   bool shapecdtn=false;
   if (exitCode) *exitCode = 255;

      if(crosstalkInv==1){
        if(NofMax==1){shapecdtn=true; if (exitCode) *exitCode=0;}
        return shapecdtn;
      }

        Float_t C_M=0.0;        Float_t C_D=0.0;        Float_t C_Mn=10000;     Float_t C_Dn=10000;     Float_t C_Mnn=10000;    Float_t C_Dnn=10000;
        Int_t CDPos;
        Float_t coeff1=1.7;     Float_t coeff2=2.0;
        Float_t coeffn=0.10;    Float_t coeffnn=0.02; Float_t noise=4.0;

        if(NofMax==1){

                if(MaxOnStart==true){
                        C_M=(Float_t)ampls[0]; C_D=(Float_t)ampls[1];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[2] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=2;}
                                else if(ampls.size()>3){ C_Dn=(Float_t)ampls[2];  C_Dnn=(Float_t)ampls[3] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=3;
                                }
                }

                if(MaxOnEnd==true){
                        C_M=(Float_t)ampls[ampls.size()-1]; C_D=(Float_t)ampls[ampls.size()-2];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[0] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=4;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[ampls.size()-3] ; C_Dnn=(Float_t)ampls[ampls.size()-4] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=5;
                                }
                }
                if(MaxInMiddle==true){
                        C_M=(Float_t)ampls[MaxPos];
                        int LeftOfMaxPos=MaxPos-1;if(LeftOfMaxPos<=0)LeftOfMaxPos=0;
                        int RightOfMaxPos=MaxPos+1;if(RightOfMaxPos>=(int)ampls.size())RightOfMaxPos=ampls.size()-1;
                        if(ampls[LeftOfMaxPos]<ampls[RightOfMaxPos]){ C_D=(Float_t)ampls[RightOfMaxPos]; C_Mn=(Float_t)ampls[LeftOfMaxPos];CDPos=RightOfMaxPos;} else{ C_D=(Float_t)ampls[LeftOfMaxPos]; C_Mn=(Float_t)ampls[RightOfMaxPos];CDPos=LeftOfMaxPos;}
                        if(C_Mn<coeff1*coeffn*C_M+coeff2*coeffnn*C_D+2*noise || C_M==255){
                                if(ampls.size()==3) shapecdtn=true ;
                                else if(ampls.size()>3){
                                        if(CDPos>MaxPos){
                                                if(ampls.size()-CDPos-1==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1==1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1>1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=(Float_t)ampls[CDPos+2];
                                                }
                                                if(MaxPos>=2){
                                                        C_Mnn=(Float_t)ampls[MaxPos-2];
                                                }
                                                else if(MaxPos<2) C_Mnn=0;
                                        }
                                        if(CDPos<MaxPos){
                                                if(CDPos==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(CDPos==1){
                                                        C_Dn=(Float_t)ampls[0];
                                                        C_Dnn=0;
                                                }
                                                if(CDPos>1){
                                                        C_Dn=(Float_t)ampls[CDPos-1];
                                                        C_Dnn=(Float_t)ampls[CDPos-2];
                                                }
                                                if(ampls.size()-LeftOfMaxPos>1 && MaxPos+2<(int)(ampls.size())-1){
                                                        C_Mnn=(Float_t)ampls[MaxPos+2];
                                                }else C_Mnn=0;
                                        }
                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                           && C_Mnn<=coeff1*coeffn*C_Mn+coeff2*coeffnn*C_M+2*noise
                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise) {
                                                shapecdtn=true;
                                        }

                                }
                        } else if (exitCode) *exitCode=6;
                }
        }
        else if (NofMax>1 && exitCode) *exitCode = 1; // more than one maximum
        if(ampls.size()==1){shapecdtn=true;}
        if(shapecdtn && exitCode) *exitCode=0;

   return shapecdtn;
}


TObject* run2analysis::GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
   size_t pos = Path.find("/");
   if(pos < 256){
      std::string firstPart = Path.substr(0,pos);
      std::string endPart   = Path.substr(pos+1,Path.length());
      TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
      if(TMP!=NULL)return GetObjectFromPath(TMP,endPart,GetACopy);

      printf("ObjectNotFound: %s::%s\n",File->GetName(), Path.c_str());
      return NULL;
   }else{
      if(GetACopy){
         return (File->Get(Path.c_str()))->Clone();
      }else{
         return File->Get(Path.c_str());
      }
   }
}


TH3F* run2analysis::loadDeDxTemplate(std::string path, std::string nameHisto, bool splitByModuleType){
   cout << "  load DeDx template --> root file : " << path.c_str() << ", histo :  " << nameHisto.c_str() << endl;
   TFile* InputFile = new TFile(path.c_str());
//   TH3F* DeDxMap_ = (TH3F*)GetObjectFromPath(InputFile, "Charge_Vs_Path");
   TH3F* DeDxMap_ = (TH3F*)GetObjectFromPath(InputFile, nameHisto.c_str());
   cout << "  load DeDx template : " << nameHisto.c_str() << endl;
   if(!DeDxMap_){printf("dEdx templates in file %s can't be open\n", path.c_str()); exit(0);}

   TH3F* Prob_ChargePath  = (TH3F*)(DeDxMap_->Clone("Prob_ChargePath"));
   Prob_ChargePath->Reset();
   Prob_ChargePath->SetDirectory(0);

   if(!splitByModuleType){
      Prob_ChargePath->RebinX(Prob_ChargePath->GetNbinsX()-1); // <-- do not include pixel in the inclusive
   }

   for(int i=0;i<=Prob_ChargePath->GetXaxis()->GetNbins()+1;i++){
      for(int j=0;j<=Prob_ChargePath->GetYaxis()->GetNbins()+1;j++){
         double Ni = 0;
         for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){Ni+=DeDxMap_->GetBinContent(i,j,k);}

         for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){
            double tmp = 0;
            for(int l=0;l<=k;l++){ tmp+=DeDxMap_->GetBinContent(i,j,l);}

            if(Ni>0){
               Prob_ChargePath->SetBinContent (i, j, k, tmp/Ni);
            }else{
               Prob_ChargePath->SetBinContent (i, j, k, 0);
            }
         }
      }
   }
   InputFile->Close();
   return Prob_ChargePath;
}

