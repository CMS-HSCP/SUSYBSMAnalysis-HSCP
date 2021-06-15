#define dEdXtemplate_cxx
#include "dEdXtemplate.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


float norm_mult = 247; // 247 or 265?

void dEdXtemplate::Loop(TString name, bool uncorr)
{
//   In a ROOT session, you can do:
//      root> .L dEdXtemplate.C
//      root> dEdXtemplate t
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

   Long64_t nentries = fChain->GetEntriesFast();

   // histo definition
   // HistoName = Name + "_ChargeVsPath"; Charge_Vs_Path = new TH3D( HistoName.c_str(), HistoName.c_str(), P_NBins, P_Min, P_Max, Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   //
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
   TH2D* Charge_Vs_Path0 = new TH2D( "Charge_Vs_Path0", "Charge_Vs_Path0", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path1 = new TH2D( "Charge_Vs_Path1", "Charge_Vs_Path1", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path2 = new TH2D( "Charge_Vs_Path2", "Charge_Vs_Path2", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path3 = new TH2D( "Charge_Vs_Path3", "Charge_Vs_Path3", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path4 = new TH2D( "Charge_Vs_Path4", "Charge_Vs_Path4", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path5 = new TH2D( "Charge_Vs_Path5", "Charge_Vs_Path5", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path6 = new TH2D( "Charge_Vs_Path6", "Charge_Vs_Path6", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path7 = new TH2D( "Charge_Vs_Path7", "Charge_Vs_Path7", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path8 = new TH2D( "Charge_Vs_Path8", "Charge_Vs_Path8", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path9 = new TH2D( "Charge_Vs_Path9", "Charge_Vs_Path9", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path10 = new TH2D( "Charge_Vs_Path10", "Charge_Vs_Path10", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path11 = new TH2D( "Charge_Vs_Path11", "Charge_Vs_Path11", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path12 = new TH2D( "Charge_Vs_Path12", "Charge_Vs_Path12", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path13 = new TH2D( "Charge_Vs_Path13", "Charge_Vs_Path13", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path14 = new TH2D( "Charge_Vs_Path14", "Charge_Vs_Path14", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* Charge_Vs_Path15 = new TH2D( "Charge_Vs_Path15", "Charge_Vs_Path15", Path_NBins, Path_Min, Path_Max, Charge_NBins, Charge_Min, Charge_Max);
   TH2D* HdedxVsP       = new TH2D( "HdedxVsP", "HdedxVsP", 500, 0, 5,1000,0, 15);
   TH2D* HdedxVsP2       = new TH2D( "HdedxVsP2", "HdedxVsP2", 500, 0, 100,1000,0, 15);

   TH1D* HHitPix = new TH1D( "HHitPix", "HHitPix", 200, 0, 20);
   TProfile* HHitProfilePix = new TProfile( "HHitProfilePix", "HHitProfilePix", 50, 0, 100); 
   TH2D* HHit2DPix = new TH2D( "HHit2DPix", "HHit2DPix", 50, 0, 100,200, 0, 20); 
   TH1D* HHitStrip = new TH1D( "HHitStrip", "HHitStrip", 200, 0, 20);
   TProfile* HHitProfileStrip = new TProfile( "HHitProfileStrip", "HHitProfileStrip", 50, 0, 100); 
   TH2D* HHit2DStrip = new TH2D( "HHit2DStrip", "HHit2DStrip", 50, 0, 100,200, 0, 20); 
   TH1D* HHitAll = new TH1D( "HHitAll", "HHitAll", 200, 0, 20);
   TProfile* HHitProfileAll = new TProfile( "HHitProfileAll", "HHitProfileAll", 50, 0, 100); 
   TH2D* HHit2DAll = new TH2D( "HHit2DAll", "HHit2DAll", 50, 0, 100,200, 0, 20); 

   TH2D* HHitPixVsP = new TH2D( "HHitPixVsP", "HHitPixVsP", 50, 0, 100, 200, 0, 20);
   TH2D* HHitStripVsP = new TH2D( "HHitStripVsP", "HHitStripVsP", 50, 0, 100, 200, 0, 20);
   TH2D* HHitAllVsP = new TH2D( "HHitAllVsP", "HHitAllVsP", 50, 0, 100, 200, 0, 20);

   TH2D* HHitPixVsPw = new TH2D( "HHitPixVsPw", "HHitPixVsPw", 10, 0, 100, 200, 0, 20);
   TH2D* HHitStripVsPw = new TH2D( "HHitStripVsPw", "HHitStripVsPw", 10, 0, 100, 200, 0, 20);
   TH2D* HHitAllVsPw = new TH2D( "HHitAllVsPw", "HHitAllVsPw", 10, 0, 100, 200, 0, 20);

   TH2D* HHitPixVsPExt = new TH2D( "HHitPixVsPExt", "HHitPixVsPExt", 50, 0, 100, 200, 0, 20);
   TH2D* HHitStripVsPExt = new TH2D( "HHitStripVsPExt", "HHitStripVsPExt", 50, 0, 100, 200, 0, 20);
   TH2D* HHitAllVsPExt = new TH2D( "HHitAllVsPExt", "HHitAllVsPExt", 50, 0, 100, 200, 0, 20);

   TH2D* HHitPixVsPExtw = new TH2D( "HHitPixVsPExtw", "HHitPixVsPExtw", 10, 0, 100, 200, 0, 20);
   TH2D* HHitStripVsPExtw = new TH2D( "HHitStripVsPExtw", "HHitStripVsPExtw", 10, 0, 100, 200, 0, 20);
   TH2D* HHitAllVsPExtw = new TH2D( "HHitAllVsPExtw", "HHitAllVsPExtw", 10, 0, 100, 200, 0, 20);

   TH2D* HHitPixVsEta = new TH2D( "HHitPixVsEta", "HHitPixVsEta", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitStripVsEta = new TH2D( "HHitStripVsEta", "HHitStripVsEta", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitAllVsEta = new TH2D( "HHitAllVsEta", "HHitAllVsEta", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitPixVsEtaP5 = new TH2D( "HHitPixVsEtaP5", "HHitPixVsEtaP5", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitStripVsEtaP5 = new TH2D( "HHitStripVsEtaP5", "HHitStripVsEtaP5", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitAllVsEtaP5 = new TH2D( "HHitAllVsEtaP5", "HHitAllVsEtaP5", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitPixVsEtapL5 = new TH2D( "HHitPixVsEtapL5", "HHitPixVsEtapL5", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitStripVsEtapL5 = new TH2D( "HHitStripVsEtapL5", "HHitStripVsEtapL5", 50, -2.5, 2.5, 200, 0, 20);
   TH2D* HHitAllVsEtapL5 = new TH2D( "HHitAllVsEtapL5", "HHitAllVsEtapL5", 50, -2.5, 2.5, 200, 0, 20);

   TH1D* Mass = new TH1D("Mass", "Mass", 150, 0.,1.5);
   TH1D* Mass2 = new TH1D("Mass2", "Mass2", 150, 0.,15);
   TH1D* Mass_0p5M0p7 = new TH1D("Mass_0p5M0p7", "Mass_0p5M0p7", 150, 0.,1.5);
   TH1D* Mass_0p7M1 = new TH1D("Mass_0p7M1", "Mass_0p7M1", 150, 0.,1.5);
   TH1D* Mass_1M2 = new TH1D("Mass_1M2", "Mass_1M2", 150, 0.,1.5);
   TH1D* Mass_2M5 = new TH1D("Mass_2M5", "Mass_2M5", 150, 0.,1.5);

   double OLDscaleFactors[2] =  { 1.09711, 1.09256 };
//   double dEdxSF [2] = { 1.09711, 1.09256 };  // 0 : Strip SF, 1 : Pixel to Strip SF
//   double dEdxSF [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF
//   MIP
//   double dEdxSF [2] = { 1., 1.1140203 };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr
//   double dEdxSF [2] = { 1., 1.2367051 };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr
//   Profile
//   double dEdxSF_uncorr [2] = { 1., 1.3201000 };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr
//   double dEdxSF_corr [2] = { 1.,  1.4327 };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr

//   double dEdxSF_uncorr [2] = { 1., 1.2025 };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr avec Profile
//   double dEdxSF_corr [2] = { 1.,  1.5583 };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr avec Profile
   double dEdxSF_uncorr [2] = { 1., 1. };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr avec Fit gaussien
   double dEdxSF_corr [2] = { 1.,  1. };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr avec Fit gaussien
//   double dEdxSF_uncorr [2] = { 1., 1.0841*1.0057  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr avec Fit gaussien in MB
//   double dEdxSF_corr [2] = { 1.,  1.1529*1.0117  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr avec Fit gaussien in mB
//   double dEdxSF_uncorr [2] = { 1., 1.0675*1.0051  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr avec Fit gaussien in tt
//   double dEdxSF_corr [2] = { 1.,  1.1339*1.0131  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr avec Fit gaussien in tt
//   double dEdxSF_uncorr [2] = { 1., 1.0141*0.9993  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr avec Fit gaussien in MB and 247
//   double dEdxSF_corr [2] = { 1.,  1.0849*1.0053  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr avec Fit gaussien in mB and 247
   double KC_uncorr[2] = { 1.988, 3.773 } ;
   double KC_corr[2] = { 2.403, 4.259 } ;
 

   TFile* OutputHisto;
   if (uncorr) OutputHisto = new TFile(name+"_template_uncorr.root","RECREATE");
   else OutputHisto = new TFile(name+"_template_corr.root","RECREATE");

   int idebug=0;
   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;

      for (int itr=0; itr<ntracks; itr++) {
      // if (Cut(ientry) < 0) continue;
          if (track_nvalidhits[itr]<8) continue;
          if (track_chi2[itr]>5) continue;


          double ih_uncorr=0;
          ih_uncorr=track_ih_ampl[itr];
          double ih_corr=0;
          ih_corr=track_ih_ampl_corr[itr];

          std::vector <float> charge_uncorr;
          std::vector <float> charge_corr;
          std::vector <float> pathlength;
          std::vector <int> subdetId;
          std::vector <int> moduleGeometry;
          std::vector <bool> bool_cleaning1;
          std::vector <bool> bool_cleaning2;
          std::vector <bool> mustBeInside;

          int nstatinfo=0;
          for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
              float ch1=dedx_charge[iclu];
              float ch2=dedx_charge[iclu];
              bool clean1=true;
              bool clean2=true;
              if (dedx_subdetid[iclu]>=3) {
                 ch1=sclus_charge[iclu];
                 ch2=sclus_charge_corr[iclu];
                 clean1=sclus_clusclean[iclu];
                 clean2=sclus_clusclean2[iclu];
                 if (sclus_sat254[iclu] || sclus_sat255[iclu]) nstatinfo++; 
              }
              charge_uncorr.push_back(ch1);
              charge_corr.push_back(ch2);
              pathlength.push_back(dedx_pathlength[iclu]);
              subdetId.push_back(dedx_subdetid[iclu]);
              moduleGeometry.push_back(dedx_modulgeom[iclu]);
              mustBeInside.push_back(dedx_insideTkMod[iclu]);
              bool_cleaning1.push_back(clean1);
              bool_cleaning2.push_back(clean2);
          }

          ih_uncorr = getdEdX(charge_uncorr, pathlength, subdetId, moduleGeometry, bool_cleaning1, mustBeInside, dEdxSF_uncorr, NULL);
          ih_corr = getdEdX(charge_corr, pathlength, subdetId, moduleGeometry, bool_cleaning2, mustBeInside, dEdxSF_corr, NULL);

          float delta2 = ih_corr - track_ih_ampl_corr[itr];
          float delta = ih_uncorr - track_ih_ampl[itr];
          if (delta<0) delta*=-1.;
          if (delta2<0) delta2*=-1.;
          if (delta>0.01 && ih_corr!=-1) {
             if (idebug<50) {
//                 cout << "probleme Ih " << ih_corr << "  " << track_ih_ampl_corr[itr] << " for entry " << jentry << endl;
/*
                 cout <<   "   nstatinfo " << nstatinfo << " nclu " << track_nhits[itr] << endl;
                 for (int ic=0; ic<charge_corr.size() ; ic++) {
                   cout << "   ch " << charge_corr[ic] << " ph " << pathlength[ic] << " sub " << subdetId[ic] << " mod " << moduleGeometry[ic] << " in "<< mustBeInside[ic] << " clean " << bool_cleaning2[ic] << endl;
*/
             }
             
            idebug++;
          }
          
          double massval=0;
          if (uncorr) massval = sqrt((ih_uncorr-KC_uncorr[1])/KC_uncorr[0])*track_p[itr];
          else  massval = sqrt((ih_corr-KC_corr[1])/KC_corr[0])*track_p[itr];

          if (track_p[itr]<5) { 
           if (uncorr) HdedxVsP ->Fill(track_p[itr], ih_uncorr);
           else HdedxVsP ->Fill(track_p[itr], ih_corr);
           Mass->Fill(massval);
           Mass2->Fill(massval);
           if (0.5<= track_p[itr] && track_p[itr]< 0.7) {
                Mass_0p5M0p7->Fill(massval);
           }
           else if (0.7 <= track_p[itr] && track_p[itr]< 1) {
                Mass_0p7M1->Fill(massval);
           }
           else if (1.0 <= track_p[itr] && track_p[itr]< 2.0) {
                Mass_1M2->Fill(massval);
           }
           else if (2.0 <= track_p[itr] && track_p[itr]< 5.0) {
                Mass_2M5->Fill(massval);
           }
          }
          if (uncorr) HdedxVsP2 ->Fill(track_p[itr], ih_uncorr);
          else HdedxVsP2 ->Fill(track_p[itr], ih_corr);

//
          for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
              if (!dedx_insideTkMod[iclu]) continue;
              if (uncorr && dedx_isstrip[iclu] && !sclus_clusclean[iclu]) continue;
              if (!uncorr && dedx_isstrip[iclu] && !sclus_clusclean2[iclu]) continue;

//              double Norm = 3.61e-06*265;
              double Norm = 3.61e-06*norm_mult;
              double scaleFactor = 1;
              int charge = 1;

              if (uncorr) {
                  scaleFactor = dEdxSF_uncorr[0];
                  if (dedx_ispixel[iclu]) scaleFactor *= dEdxSF_uncorr[1];
                  charge = sclus_charge[iclu];
              } 
              else {
                scaleFactor = dEdxSF_corr[0];
                if (dedx_ispixel[iclu]) scaleFactor *= dEdxSF_corr[1];
                charge = sclus_charge_corr[iclu];
              }
              if (dedx_ispixel[iclu]) {
                 Norm = 3.61e-06;
                 charge = dedx_charge[iclu];
              }
              
              double ChargeOverPathlength = scaleFactor*Norm*charge/dedx_pathlength[iclu];
              if (dedx_ispixel[iclu]) HHitPixVsEta->Fill(track_eta[itr],ChargeOverPathlength);
              else  HHitStripVsEta->Fill(track_eta[itr],ChargeOverPathlength);
              HHitAllVsEta->Fill(track_eta[itr],ChargeOverPathlength);
              if (track_p[itr]<5) {
               if (dedx_ispixel[iclu]) HHitPixVsEtapL5->Fill(track_eta[itr],ChargeOverPathlength);
               else  HHitStripVsEtapL5->Fill(track_eta[itr],ChargeOverPathlength);
               HHitAllVsEtapL5->Fill(track_eta[itr],ChargeOverPathlength);
              }
              else {
               if (dedx_ispixel[iclu]) HHitPixVsEtaP5->Fill(track_eta[itr],ChargeOverPathlength);
               else  HHitStripVsEtaP5->Fill(track_eta[itr],ChargeOverPathlength);
               HHitAllVsEtaP5->Fill(track_eta[itr],ChargeOverPathlength);
              }
              if (fabs(track_eta[itr])<0.4){
	          HHitAllVsP->Fill(track_p[itr],ChargeOverPathlength);
	          HHitAllVsPw->Fill(track_p[itr],ChargeOverPathlength);
                  if (dedx_ispixel[iclu]){
                     HHitPixVsP->Fill(track_p[itr],ChargeOverPathlength);
                     HHitPixVsPw->Fill(track_p[itr],ChargeOverPathlength);
                  }
                  else {
                     HHitStripVsP->Fill(track_p[itr], ChargeOverPathlength);
                     HHitStripVsPw->Fill(track_p[itr], ChargeOverPathlength);
                  }
                  
              }
              else {
	          HHitAllVsPExt->Fill(track_p[itr],ChargeOverPathlength);
	          HHitAllVsPExtw->Fill(track_p[itr],ChargeOverPathlength);
                  if (dedx_ispixel[iclu]){
                     HHitPixVsPExt->Fill(track_p[itr],ChargeOverPathlength);
                     HHitPixVsPExtw->Fill(track_p[itr],ChargeOverPathlength);
                  }
                  else {
                     HHitStripVsPExt->Fill(track_p[itr], ChargeOverPathlength);
                     HHitStripVsPExtw->Fill(track_p[itr], ChargeOverPathlength);
                  }
              }
          }
//

          if (track_pt[itr]<5) continue;

          for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
              //condition isInsideTkMod
              if (!dedx_insideTkMod[iclu]) continue;
              //condition ClusterCleaning
              if (uncorr && dedx_isstrip[iclu] && !sclus_clusclean[iclu]) continue;
              if (!uncorr && dedx_isstrip[iclu] && !sclus_clusclean2[iclu]) continue;

              float norma=1;
//              double Norm = 3.61e-06*265;
              double Norm = 3.61e-06*norm_mult;
              double scaleFactor = 1;
              int charge = 1;

              if (uncorr) {
                  scaleFactor = dEdxSF_uncorr[0];
                  if (dedx_ispixel[iclu]) scaleFactor *= dEdxSF_uncorr[1];
                  charge = sclus_charge[iclu];
              } 
              else {
                scaleFactor = dEdxSF_corr[0];
                if (dedx_ispixel[iclu]) scaleFactor *= dEdxSF_corr[1];
                charge = sclus_charge_corr[iclu];
              }
              if (dedx_ispixel[iclu]) {
//                 norma=265;
                 norma=norm_mult;
                 Norm = 3.61e-06;
                 charge = dedx_charge[iclu];
              }
              
              // results[R]->Charge_Vs_Path->Fill (moduleGeometry, dedxHits->pathlength(h)*10, scaleFactor * charge/(dedxHits->pathlength(h)*10*(detid.subdetId()<3?265:1))); 
              Charge_Vs_Path->Fill (dedx_modulgeom[iclu], dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              Charge_Vs_Path0->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              if (dedx_modulgeom[iclu]==1) Charge_Vs_Path1->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==2) Charge_Vs_Path2->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==3) Charge_Vs_Path3->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==4) Charge_Vs_Path4->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==5) Charge_Vs_Path5->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==6) Charge_Vs_Path6->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==7) Charge_Vs_Path7->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==8) Charge_Vs_Path8->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==9) Charge_Vs_Path9->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==10) Charge_Vs_Path10->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==11) Charge_Vs_Path11->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==12) Charge_Vs_Path12->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==13) Charge_Vs_Path13->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==14) Charge_Vs_Path14->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));
              else if (dedx_modulgeom[iclu]==15) Charge_Vs_Path15->Fill (dedx_pathlength[iclu]*10, scaleFactor* charge/(dedx_pathlength[iclu]*10*norma));

              double ChargeOverPathlength = scaleFactor*Norm*charge/dedx_pathlength[iclu];
              if (dedx_ispixel[iclu]) { 
                  HHitPix->Fill(ChargeOverPathlength);
                  if (fabs(track_eta[itr])<0.4){
                    HHitProfilePix->Fill(track_p[itr], ChargeOverPathlength);
                    HHit2DPix->Fill(track_p[itr], ChargeOverPathlength);
                  }
              }
              else {
                  HHitStrip->Fill(ChargeOverPathlength);
                  if (fabs(track_eta[itr])<0.4){
                    HHitProfileStrip->Fill(track_p[itr], ChargeOverPathlength);
                    HHit2DStrip->Fill(track_p[itr], ChargeOverPathlength);
                  }
              }
              HHitAll->Fill(ChargeOverPathlength);
              if (fabs(track_eta[itr])<0.4){
                 HHitProfileAll->Fill(track_p[itr], ChargeOverPathlength);
                 HHit2DAll->Fill(track_p[itr], ChargeOverPathlength);
              }
          }
      }
   }
   cout << " nombre de fois ou ih_corr n est pas exactement ce qui est sauvegarde " << idebug << endl;
   OutputHisto->cd();
   Charge_Vs_Path->Write();
   Charge_Vs_Path0->Write();
   Charge_Vs_Path1->Write();
   Charge_Vs_Path2->Write();
   Charge_Vs_Path3->Write();
   Charge_Vs_Path4->Write();
   Charge_Vs_Path5->Write();
   Charge_Vs_Path6->Write();
   Charge_Vs_Path7->Write();
   Charge_Vs_Path8->Write();
   Charge_Vs_Path9->Write();
   Charge_Vs_Path10->Write();
   Charge_Vs_Path11->Write();
   Charge_Vs_Path12->Write();
   Charge_Vs_Path13->Write();
   Charge_Vs_Path14->Write();
   Charge_Vs_Path15->Write();
   HdedxVsP->Write();
   HdedxVsP2->Write();
   HHitPix->Write();
   HHitProfilePix->Write();
   HHit2DPix->Write();
   HHitStrip->Write();
   HHitProfileStrip->Write();
   HHit2DStrip->Write();
   HHitAll->Write();
   HHitProfileAll->Write();
   HHit2DAll->Write();
   HHitPixVsP->Write();
   HHitStripVsP->Write();
   HHitAllVsP->Write();
   HHitPixVsPw->Write();
   HHitStripVsPw->Write();
   HHitAllVsPw->Write();
   HHitPixVsPExt->Write();
   HHitStripVsPExt->Write();
   HHitAllVsPExt->Write();
   HHitPixVsPExtw->Write();
   HHitStripVsPExtw->Write();
   HHitAllVsPExtw->Write();
   HHitPixVsEta->Write();
   HHitStripVsEta->Write();
   HHitAllVsEta->Write();
   HHitPixVsEtapL5->Write();
   HHitStripVsEtapL5->Write();
   HHitAllVsEtapL5->Write();
   HHitPixVsEtaP5->Write();
   HHitStripVsEtaP5->Write();
   HHitAllVsEtaP5->Write();
   Mass->Write();
   Mass2->Write();
   Mass_0p5M0p7->Write();
   Mass_0p7M1->Write();
   Mass_1M2->Write();
   Mass_2M5->Write();
   OutputHisto->Close();
}



double dEdXtemplate::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto) {
  double result=-1;

//reco::DeDxData computedEdx(const reco::DeDxHitInfo* dedxHits, double* scaleFactors, TH3* templateHisto, bool usePixel, bool useClusterCleaning, bool reverseProb, bool useTruncated, std::unordered_map<unsigned int,double>* TrackerGains, bool useStrip, bool mustBeInside, size_t MaxStripNOM, bool correctFEDSat, int crossTalkInvAlgo, double dropLowerDeDxValue,  double* dEdxErr){

//reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, NULL, true, true, false , false, NULL, true, true, 99, false, 0, 0.15,  NULL);
//reco::DeDxData dedxMObjTmp2 = computedEdx(dedxHits, dEdxSF_corr, NULL, true, true, false , false, NULL, true, true, 99, false, 1, 0.15,  NULL);


     double dropLowerDeDxValue=0.15;
     size_t MaxStripNOM=99;
     bool usePixel=true;
     bool useStrip=true;
     
     std::vector<double> vect;

     bool debugprint=false;
     unsigned int SiStripNOM = 0;

     for(unsigned int h=0;h<charge.size();h++){
        if (debugprint) std::cout << "look on dedxHits in computedEdx " << h << std::endl;
        if(!usePixel && subdetId[h]<3)continue; // skip pixels
        if(!useStrip && subdetId[h]>=3)continue; // skip strips        
        if(useStrip && subdetId[h]>=3 && !bool_cleaning[h])continue;
        if(useStrip && subdetId[h]>=3 && !mustBeInside[h])continue;
        if(useStrip && subdetId[h]>=3 && ++SiStripNOM > MaxStripNOM) continue; // skip remaining strips, but not pixel

        int ClusterCharge = charge[h];

        double scaleFactor = scaleFactors[0];
        if (subdetId[h]<3) scaleFactor *= scaleFactors[1]; // add pixel scaling
        if (debugprint) std::cout << " after SF " << std::endl;

        if(templateHisto){  //save discriminator probability
//           double ChargeOverPathlength = scaleFactor*ClusterCharge/(pathlength[h]*10.0*(subdetId[h]<3?265:1));
           double ChargeOverPathlength = scaleFactor*ClusterCharge/(pathlength[h]*10.0*(subdetId[h]<3?norm_mult:1));
           int    BinX   = templateHisto->GetXaxis()->FindBin(moduleGeometry[h]);
           int    BinY   = templateHisto->GetYaxis()->FindBin(pathlength[h]*10.0); //*10 because of cm-->mm
           int    BinZ   = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
           double Prob   = templateHisto->GetBinContent(BinX,BinY,BinZ);
           vect.push_back(Prob); //save probability
           if (debugprint) std::cout << " after Prob vect.push_back " << std::endl;
        }else{
//           double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*265;
           double Norm = (subdetId[h]<3)?3.61e-06:3.61e-06*norm_mult;
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
           double expo = -2;
           for(int i = 0; i< size; i ++){
              result+=pow(vect[i],expo);
           }
           result = pow(result/size,1./expo);
           if (debugprint) std::cout << " harmonic2 discriminator " << result << std::endl;
        }
     }else{
        result = -1;
     }
     if (debugprint) std::cout << " ok finished computeDeDx " << std::endl;


  return result;
}
