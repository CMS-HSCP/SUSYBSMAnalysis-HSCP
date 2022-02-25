#define GluinoAnalysis_cxx
#include "GluinoAnalysis.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>

void GluinoAnalysis::Loop()
{
//   In a ROOT session, you can do:
//      root> .L GluinoAnalysis.C
//      root> GluinoAnalysis t
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

   double KC_uncorr[2] = { 1.988, 3.773 } ;
   double KC_corr[2] = { 2.403, 4.259 } ;
   double dEdxSF_uncorr [2] = { 1., 1.0841*1.0057  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr avec Fit gaussien in MB
   double dEdxSF_corr [2] = { 1.,  1.1529*1.0117  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr avec Fit gaussien in mB
//    double dEdxSF_uncorr [2] = { 1., 1.0675*1.0051  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 uncorr avec Fit gaussien in tt
//    double dEdxSF_corr [2] = { 1.,  1.1339*1.0131  };  // 0 : Strip SF, 1 : Pixel to Strip SF  // iteration 1 corr avec Fit gaussien in tt


   TH1D* Mass = new TH1D("Mass", "Mass", 200, 0.,5000);
   TH1D* Mass2 = new TH1D("Mass2", "Mass2", 200, 0.,5000);
   TH1D* Nsat = new TH1D("Nsat", "Nsat", 20, 0.,20);
   TH2D* MassVsNsat = new TH2D("MassVsNsat", "MassVsNsat", 20, 0.,5000,20,0,20);
   TH2D* MassVsNsat2 = new TH2D("MassVsNsat2", "MassVsNsat2", 20, 0.,5000,20,0,20);
   TH2D* MassVsP = new TH2D("MassVsP", "MassVsP", 20, 0.,5000,50,0,5000);
   TH2D* MassVsP2 = new TH2D("MassVsP2", "MassVsP2", 20, 0.,5000,50,0,5000);
   TH2D* dEdXVsP = new TH2D("dEdXVsP", "dEdXVsP", 20,0,5000, 40, 0.,20);
   TH2D* dEdXVsP2 = new TH2D("dEdXVsP2", "dEdXVsP2", 20,0,5000, 40, 0.,20);
   TH2D* dEdXVsP3 = new TH2D("dEdXVsP3", "dEdXVsP3", 20,0,5000, 40, 0.,20);


   TFile* OutputHisto = new TFile("testMass.root","RECREATE");

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      //
      for (int itr=0; itr<ntracks; itr++) {
         if (track_nvalidhits[itr]<8) continue;
         if (track_chi2[itr]>5) continue;
         if (track_pt[itr]<60) continue;
         int nsatclu=0;
         for (int iclu=track_index_hit[itr]; iclu<track_index_hit[itr]+track_nhits[itr]; iclu++) {
            if (!dedx_insideTkMod[iclu]) continue;
            if (dedx_isstrip[iclu] && !sclus_clusclean[iclu]) continue;
            if (dedx_isstrip[iclu] && (sclus_sat254[iclu] || sclus_sat255[iclu])) nsatclu++;
         }
         double massval= sqrt((track_ih_ampl[itr]-KC_uncorr[1])/KC_uncorr[0])*track_p[itr];
         double masscor= sqrt((track_ih_ampl_corr[itr]-KC_corr[1])/KC_corr[0])*track_p[itr];
         if (track_p[itr]>200) Mass->Fill(massval);
         if (track_p[itr]>200) Mass2->Fill(masscor);
         Nsat->Fill(nsatclu);
         MassVsNsat->Fill(massval,nsatclu);
         MassVsNsat2->Fill(masscor,nsatclu);
         MassVsP->Fill(track_p[itr],massval);
         MassVsP2->Fill(track_p[itr],masscor);
         dEdXVsP->Fill(track_p[itr],track_ih_ampl[itr]);
         dEdXVsP2->Fill(track_p[itr],track_ih_ampl_corr[itr]);



         //test caro
          std::vector <float> charge_uncorr;
          std::vector <float> charge_corr;
          std::vector <float> charge_eloss;
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
          // test caro



      }

   }
   OutputHisto->cd();
   Mass->Write();  
   Mass2->Write();  
   Nsat->Write();  
   MassVsNsat2->Write();
   MassVsNsat->Write();
   MassVsP->Write();
   MassVsP2->Write();
   dEdXVsP->Write();
   dEdXVsP2->Write();
   dEdXVsP3->Write();
   OutputHisto->Close();
}


double dEdXtemplate::getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto) {
  double result=-1;

//reco::DeDxData computedEdx(const reco::DeDxHitInfo* dedxHits, double* scaleFactors, TH3* templateHisto, bool usePixel, bool useClusterCleaning, bool reverseProb, bool useTruncated, std::unordered_map<unsigned int,double>* TrackerGains, bool useStrip, bool mustBeInside, size_t MaxStripNOM, bool correctFEDSat, int crossTalkInvAlgo, double dropLowerDeDxValue,  double* dEdxErr){

//reco::DeDxData dedxMObjTmp = computedEdx(dedxHits, dEdxSF, NULL, true, true, false , false, NULL, true, true, 99, false, 0, 0.15,  NULL);
//reco::DeDxData dedxMObjTmp2 = computedEdx(dedxHits, dEdxSF_corr, NULL, true, true, false , false, NULL, true, true, 99, false, 1, 0.15,  NULL);
//
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

