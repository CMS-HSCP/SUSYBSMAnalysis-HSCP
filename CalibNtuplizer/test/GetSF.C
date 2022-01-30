// inspired by https://github.com/CMS-HSCP/SUSYBSMAnalysis-HSCP/blob/Run2_2016/test/UsefulScripts/DeDxStudy/MakePlot.C#L436:L572

#include <exception>
#include <vector>
#include <fstream>

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
#include "TGraphErrors.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TCutG.h"


void SaveCanvas(TCanvas* c, std::string path, std::string name, bool OnlyPPNG=false);

void GetSF() {

   TProfile*   HdedxVsPProfile1;
   TProfile*   HdedxVsPProfile2;
   TH1D*       HdedxMIP1;
   TH1D*       HdedxMIP2;


   // ouvrir le fichier  et recuperer les objets
   TFile* myfile = new TFile ("minbias_template_uncorr_aod4_ch.root");
//   TFile* myfile = new TFile ("minbias_template_corr_iter1.root");
   myfile->cd();

   HdedxVsPProfile2 = (TProfile*) gROOT->FindObject("HHitProfilePix");
   HdedxVsPProfile1 = (TProfile*) gROOT->FindObject("HHitProfileStrip");

   HdedxMIP2        = (TH1D*)     gROOT->FindObject("HHitPix");
   HdedxMIP1        = (TH1D*)     gROOT->FindObject("HHitStrip");


   //

   double SFMip, SFProfile;
   std::string SaveDir = "SFdir";
   FILE* fout = fopen ((SaveDir+"ScaleFactors.txt").c_str(), "w");
   fprintf (fout, "=================\n >>> Old CCC >>>\n=================\n");

   TF1* mygausMIP = new TF1("mygausMIP","gaus", 1, 5);
   HdedxMIP1->Fit("mygausMIP","Q0","");
   double peakMIP1  = mygausMIP->GetParameter(1);
   HdedxMIP2->Fit("mygausMIP","Q0","");
   double peakMIP2  = mygausMIP->GetParameter(1);

   std::cout << "SCALE FACTOR WITH MIP     = " << peakMIP1/peakMIP2 << endl;
   SFMip = peakMIP1/peakMIP2;

   TH1D* Chi2Dist = new TH1D("Chi2Dist","Chi2Dist",5000, 0.8 ,1.8);

   double Minimum = 999999;
   double AbsGain = -1;

   for(int i=1;i<=Chi2Dist->GetNbinsX();i++){
      double ScaleFactor = Chi2Dist->GetXaxis()->GetBinCenter(i);
      TProfile* Rescaled = (TProfile*)HdedxVsPProfile2->Clone("Cloned");
      Rescaled->Scale(ScaleFactor);
      double Dist = 0;
      double Error = 0;
      for(int x=1;x<=HdedxVsPProfile1->GetNbinsX();x++){
         double Momentum = HdedxVsPProfile1->GetXaxis()->GetBinCenter(x);
         if(Momentum<5)continue;  //|| Momentum>20)continue;
         if(HdedxVsPProfile1->GetBinError(x)<=0)continue;
         Dist += pow(HdedxVsPProfile1->GetBinContent(x) - Rescaled->GetBinContent(x),2) / std::max(1E-8,pow(HdedxVsPProfile1->GetBinError(x),2));
         Error += pow(HdedxVsPProfile1->GetBinError(x),2);
      }
      Dist *= Error;

      if(Dist<Minimum){Minimum=Dist;AbsGain=ScaleFactor;}

      //std::cout << "Rescale = " << ScaleFactor << " --> SquareDist = " << Dist << endl;
      Chi2Dist->Fill(ScaleFactor,Dist);
      delete Rescaled;
   }

   std::cout << "SCALE FACTOR WITH PROFILE = " << AbsGain << endl;
   SFProfile = AbsGain;

   TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
   TLegend* leg = new TLegend (0.50, 0.75, 0.80, 0.90);
   leg->SetHeader ("Fitting the MIP");
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   HdedxMIP2->SetStats(kFALSE);
   HdedxMIP2->SetAxisRange(0,10,"X");
   HdedxMIP2->GetXaxis()->SetNdivisions(516);
   HdedxMIP2->GetXaxis()->SetTitle("dE/dx (MeV/cm)");
   HdedxMIP2->GetYaxis()->SetTitle("Tracker hits");
   HdedxMIP2->SetLineColor(1);
   if (HdedxMIP2->GetMaximum() < HdedxMIP1->GetMaximum()) HdedxMIP2->SetMaximum(HdedxMIP1->GetMaximum()*3);
   HdedxMIP2->Draw("");
   TH1D* HdedxMIP3 = (TH1D*)HdedxMIP2->Clone("aaa");
   HdedxMIP3->SetLineColor(8);
   HdedxMIP3->GetXaxis()->Set(HdedxMIP3->GetXaxis()->GetNbins(), HdedxMIP3->GetXaxis()->GetXmin()*2.0, HdedxMIP3->GetXaxis()->GetXmax()*(peakMIP1/peakMIP2) );
//   HdedxMIP3->Draw("same");
   TH1D* HdedxMIP4 = (TH1D*)HdedxMIP2->Clone("bbb");
   HdedxMIP4->SetLineColor(4);
   HdedxMIP4->GetXaxis()->Set(HdedxMIP4->GetXaxis()->GetNbins(), HdedxMIP4->GetXaxis()->GetXmin()*2.0, HdedxMIP4->GetXaxis()->GetXmax()*(AbsGain) );
   HdedxMIP4->Draw("same");
   HdedxMIP1->SetLineColor(2);
   HdedxMIP1->Draw("same");
   c1->SetLogy(true);
   c1->SetGridx(true); 
   leg->AddEntry (HdedxMIP1, "Strip", "L");
   leg->AddEntry (HdedxMIP2, "unscaled Pixel", "L");
   leg->AddEntry (HdedxMIP4, "scaled   Pixel", "L");
   leg->Draw();
   SaveCanvas(c1, SaveDir, "Rescale_HdedxMIP");
   delete c1;


   c1 = new TCanvas("c1", "c1", 600,600);
   Chi2Dist->SetStats(kFALSE);
   Chi2Dist->GetXaxis()->SetNdivisions(504);
   Chi2Dist->GetXaxis()->SetTitle("Rescale Factor");
   Chi2Dist->GetYaxis()->SetTitle("Weighted Square Distance");
   Chi2Dist->Draw("");
   c1->SetLogy(true);
   c1->SetGridx(true); 
   SaveCanvas(c1, SaveDir, "Rescale_Chi2Dist");
   delete c1;

   c1 = new TCanvas("c1", "c1", 600,600);
//   leg = new TLegend (0.30, 0.25, 0.80, 0.55);
   leg = new TLegend (0.25, 0.70, 0.65, 0.90);
   leg->SetHeader ("Fitting the Profile");
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   HdedxVsPProfile1->SetStats(kFALSE);
   HdedxVsPProfile1->SetAxisRange(5,50,"X");
//   HdedxVsPProfile1->SetAxisRange(3.0,5.5,"Y");
   HdedxVsPProfile1->SetAxisRange(0.0,10.,"Y");
   HdedxVsPProfile1->GetXaxis()->SetTitle("track momentum (GeV)");
   HdedxVsPProfile1->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
   HdedxVsPProfile1->SetMarkerColor(2);
   HdedxVsPProfile1->SetLineColor(2);
   HdedxVsPProfile1->Draw("");

   HdedxVsPProfile2->SetMarkerColor(1);
   HdedxVsPProfile2->Draw("same");
   TProfile* HdedxVsPProfile3 = (TProfile*)HdedxVsPProfile2->Clone("abc");
   HdedxVsPProfile3->SetMarkerColor(8);
   HdedxVsPProfile3->SetLineColor(8);
   HdedxVsPProfile3->Scale(peakMIP1/peakMIP2);
//   HdedxVsPProfile3->Draw("same");
   TProfile* HdedxVsPProfile4 = (TProfile*)HdedxVsPProfile2->Clone("afs");
   HdedxVsPProfile4->SetMarkerColor(4);
   HdedxVsPProfile4->SetLineColor(4);
   HdedxVsPProfile4->Scale(AbsGain);
   HdedxVsPProfile4->Draw("same");
   leg->AddEntry (HdedxVsPProfile1, "Strip", "LP");
   leg->AddEntry (HdedxVsPProfile2, "unscaled Pixel", "LP");
   leg->AddEntry (HdedxVsPProfile4, "scaled Pixel",   "LP");
   leg->Draw();

//   DrawPreliminary("", 13, "");
   SaveCanvas(c1, SaveDir, "Rescale_HdedxVsPProfile");
   delete c1;

   fprintf (fout, "HHit Pixel to Strip :: MIP: %.7lf \t Profile: %.7lf \n", SFMip, SFProfile);
   fclose (fout);

}


void SaveCanvas(TCanvas* c, std::string path, std::string name, bool OnlyPPNG){
   std::string tmppath = path;
   if(tmppath[tmppath.length()-1]!='/')tmppath += "_";
   tmppath += name;

   std::string filepath;
   filepath = tmppath + ".png"; c->SaveAs(filepath.c_str()); if(OnlyPPNG)return;
   filepath = tmppath +  ".eps"; c->SaveAs(filepath.c_str());
   filepath = tmppath + ".C"  ; c->SaveAs(filepath.c_str());
   filepath = tmppath +  ".pdf"; c->SaveAs(filepath.c_str());
}
