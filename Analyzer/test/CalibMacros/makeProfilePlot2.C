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

using namespace std;
void MakeMyProfile(TH2D* input, TH1D* output);
void SaveCanvas(TCanvas* c, std::string path, std::string name, bool OnlyPPNG=false);
void plot(bool mb);
void GetSF(bool mb, bool corr);
void GetSF();


void plot( bool mb) {
TFile *_file0 ;
TFile *_file1;

if (mb)  {
_file0= TFile::Open("minbias_template_corr.root");
_file1 = TFile::Open("minbias_template_uncorr.root");
}
else {
_file0 = TFile::Open("ttbar_template_corr.root");
_file1 = TFile::Open("ttbar_template_uncorr.root");
}

TH2D* CorrHHitPixVsEtaP5;
TH2D* CorrHHitPixVsEtapL5;

TH2D* UnCorrHHitPixVsEtaP5;
TH2D* UnCorrHHitPixVsEtapL5;

TH2D* CorrHHitStripVsEtaP5;
TH2D* CorrHHitStripVsEtapL5;

TH2D* UnCorrHHitStripVsEtaP5;
TH2D* UnCorrHHitStripVsEtapL5;


TH2D *CorrHHitPixVsP;
TH2D *CorrHHitStripVsP;
TH2D *UnCorrHHitPixVsP;
TH2D *UnCorrHHitStripVsP;
TH2D *CorrHHitPixVsPExt;
TH2D *CorrHHitStripVsPExt;
TH2D *UnCorrHHitPixVsPExt;
TH2D *UnCorrHHitStripVsPExt;

_file0->cd();


CorrHHitStripVsEtaP5 = (TH2D*)gROOT->FindObject("HHitStripVsEtaP5");
CorrHHitStripVsEtapL5 = (TH2D*) gROOT->FindObject("HHitStripVsEtapL5");
CorrHHitPixVsEtaP5 = (TH2D*) gROOT->FindObject("HHitPixVsEtaP5");
CorrHHitPixVsEtapL5 = (TH2D*) gROOT->FindObject("HHitPixVsEtapL5");

TH1D* CorrPixVsEtaP5 = new TH1D("CorrPixVsEtaP5"       , "CorrPixVsEtaP5"      ,CorrHHitPixVsEtaP5->GetXaxis()->GetNbins(),CorrHHitPixVsEtaP5->GetXaxis()->GetXmin(),CorrHHitPixVsEtaP5->GetXaxis()->GetXmax());
TH1D* CorrPixVsEtapL5 = new TH1D("CorrPixVsEtapL5"       , "CorrPixVsEtapL5"      ,CorrHHitPixVsEtapL5->GetXaxis()->GetNbins(),CorrHHitPixVsEtapL5->GetXaxis()->GetXmin(),CorrHHitPixVsEtapL5->GetXaxis()->GetXmax());
TH1D* CorrStripVsEtaP5 = new TH1D("CorrStripVsEtaP5"       , "CorrStripVsEtaP5"      ,CorrHHitStripVsEtaP5->GetXaxis()->GetNbins(),CorrHHitStripVsEtaP5->GetXaxis()->GetXmin(),CorrHHitStripVsEtaP5->GetXaxis()->GetXmax());
TH1D* CorrStripVsEtapL5 = new TH1D("CorrStripVsEtapL5"       , "CorrStripVsEtapL5"      ,CorrHHitStripVsEtapL5->GetXaxis()->GetNbins(),CorrHHitStripVsEtapL5->GetXaxis()->GetXmin(),CorrHHitStripVsEtapL5->GetXaxis()->GetXmax());

MakeMyProfile(CorrHHitPixVsEtaP5, CorrPixVsEtaP5);
MakeMyProfile(CorrHHitPixVsEtapL5, CorrPixVsEtapL5);
MakeMyProfile(CorrHHitStripVsEtaP5, CorrStripVsEtaP5);
MakeMyProfile(CorrHHitStripVsEtapL5, CorrStripVsEtapL5);

CorrHHitPixVsP= (TH2D*)gROOT->FindObject("HHitPixVsPw");
CorrHHitStripVsP= (TH2D*)gROOT->FindObject("HHitStripVsPw");
TH1D* CorrPixVsP = new TH1D("CorrPixVsP"       , "CorrPixVsP"      ,CorrHHitPixVsP->GetXaxis()->GetNbins(),CorrHHitPixVsP->GetXaxis()->GetXmin(),CorrHHitPixVsP->GetXaxis()->GetXmax());
TH1D* CorrStripVsP = new TH1D("CorrStripVsP"       , "CorrStripVsP"      ,CorrHHitStripVsP->GetXaxis()->GetNbins(),CorrHHitStripVsP->GetXaxis()->GetXmin(),CorrHHitStripVsP->GetXaxis()->GetXmax());
MakeMyProfile(CorrHHitPixVsP,CorrPixVsP);
MakeMyProfile(CorrHHitStripVsP,CorrStripVsP);
CorrHHitPixVsPExt= (TH2D*)gROOT->FindObject("HHitPixVsPExtw");
CorrHHitStripVsPExt= (TH2D*)gROOT->FindObject("HHitStripVsPExtw");
TH1D* CorrPixVsPExt = new TH1D("CorrPixVsPExt"       , "CorrPixVsPExt"      ,CorrHHitPixVsPExt->GetXaxis()->GetNbins(),CorrHHitPixVsPExt->GetXaxis()->GetXmin(),CorrHHitPixVsPExt->GetXaxis()->GetXmax());
TH1D* CorrStripVsPExt = new TH1D("CorrStripVsPExt"       , "CorrStripVsPExt"      ,CorrHHitStripVsPExt->GetXaxis()->GetNbins(),CorrHHitStripVsPExt->GetXaxis()->GetXmin(),CorrHHitStripVsPExt->GetXaxis()->GetXmax());
MakeMyProfile(CorrHHitPixVsPExt,CorrPixVsPExt);
MakeMyProfile(CorrHHitStripVsPExt,CorrStripVsPExt);

cout <<  "fill from file0"<< endl;
_file1->cd();


UnCorrHHitStripVsEtaP5 = (TH2D*)gROOT->FindObject("HHitStripVsEtaP5");
UnCorrHHitStripVsEtapL5 = (TH2D*) gROOT->FindObject("HHitStripVsEtapL5");
UnCorrHHitPixVsEtaP5 = (TH2D*) gROOT->FindObject("HHitPixVsEtaP5");
UnCorrHHitPixVsEtapL5 = (TH2D*) gROOT->FindObject("HHitPixVsEtapL5");

TH1D* UnCorrPixVsEtaP5 = new TH1D("UnCorrPixVsEtaP5"       , "UnCorrPixVsEtaP5"      ,UnCorrHHitPixVsEtaP5->GetXaxis()->GetNbins(),UnCorrHHitPixVsEtaP5->GetXaxis()->GetXmin(),UnCorrHHitPixVsEtaP5->GetXaxis()->GetXmax());
TH1D* UnCorrPixVsEtapL5 = new TH1D("UnCorrPixVsEtapL5"       , "UnCorrPixVsEtapL5"      ,UnCorrHHitPixVsEtapL5->GetXaxis()->GetNbins(),UnCorrHHitPixVsEtapL5->GetXaxis()->GetXmin(),UnCorrHHitPixVsEtapL5->GetXaxis()->GetXmax());
TH1D* UnCorrStripVsEtaP5 = new TH1D("UnCorrStripVsEtaP5"       , "UnCorrStripVsEtaP5"      ,UnCorrHHitStripVsEtaP5->GetXaxis()->GetNbins(),UnCorrHHitStripVsEtaP5->GetXaxis()->GetXmin(),UnCorrHHitStripVsEtaP5->GetXaxis()->GetXmax());
TH1D* UnCorrStripVsEtapL5 = new TH1D("UnCorrStripVsEtapL5"       , "UnCorrStripVsEtapL5"      ,UnCorrHHitStripVsEtapL5->GetXaxis()->GetNbins(),UnCorrHHitStripVsEtapL5->GetXaxis()->GetXmin(),UnCorrHHitStripVsEtapL5->GetXaxis()->GetXmax());

MakeMyProfile(UnCorrHHitPixVsEtaP5,UnCorrPixVsEtaP5);
MakeMyProfile(UnCorrHHitPixVsEtapL5, UnCorrPixVsEtapL5);
MakeMyProfile(UnCorrHHitStripVsEtaP5, UnCorrStripVsEtaP5);
MakeMyProfile(UnCorrHHitStripVsEtapL5, UnCorrStripVsEtapL5);

UnCorrHHitPixVsP= (TH2D*)gROOT->FindObject("HHitPixVsPw");
UnCorrHHitStripVsP= (TH2D*)gROOT->FindObject("HHitStripVsPw");
TH1D* UnCorrPixVsP = new TH1D("UnCorrPixVsP5"       , "UnCorrPixVsP5"      ,UnCorrHHitPixVsP->GetXaxis()->GetNbins(),UnCorrHHitPixVsP->GetXaxis()->GetXmin(),UnCorrHHitPixVsP->GetXaxis()->GetXmax());
TH1D* UnCorrStripVsP = new TH1D("UnCorrStripVsP5"       , "UnCorrStripVsP5"      ,UnCorrHHitStripVsP->GetXaxis()->GetNbins(),UnCorrHHitStripVsP->GetXaxis()->GetXmin(),UnCorrHHitStripVsP->GetXaxis()->GetXmax());
MakeMyProfile(UnCorrHHitPixVsP,UnCorrPixVsP);
MakeMyProfile(UnCorrHHitStripVsP,UnCorrStripVsP);
UnCorrHHitPixVsPExt= (TH2D*)gROOT->FindObject("HHitPixVsPExtw");
UnCorrHHitStripVsPExt= (TH2D*)gROOT->FindObject("HHitStripVsPExtw");
TH1D* UnCorrPixVsPExt = new TH1D("UnCorrPixVsPExt"       , "UnCorrPixVsPExt"      ,UnCorrHHitPixVsPExt->GetXaxis()->GetNbins(),UnCorrHHitPixVsPExt->GetXaxis()->GetXmin(),UnCorrHHitPixVsPExt->GetXaxis()->GetXmax());
TH1D* UnCorrStripVsPExt = new TH1D("UnCorrStripVsPExt"       , "UnCorrStripVsPExt"      ,UnCorrHHitStripVsPExt->GetXaxis()->GetNbins(),UnCorrHHitStripVsPExt->GetXaxis()->GetXmin(),UnCorrHHitStripVsPExt->GetXaxis()->GetXmax());
MakeMyProfile(UnCorrHHitPixVsPExt,UnCorrPixVsPExt);
MakeMyProfile(UnCorrHHitStripVsPExt,UnCorrStripVsPExt);

cout <<  "fill from file1"<< endl;


TCanvas *c5 = new TCanvas("c5", "EtaDep Corr", 800,600);
c5->Divide(2,2);                       
c5->cd(2);                       
CorrHHitPixVsEtaP5->SetStats(kFALSE);
CorrHHitPixVsEtaP5 ->Draw("colz");
CorrPixVsEtaP5->Draw("same");
TProfile* pr1 = CorrHHitPixVsEtaP5->ProfileX();
pr1->SetLineColor(2);
pr1->Draw("same");

c5->cd(1);                       
CorrHHitPixVsEtapL5->SetStats(kFALSE);
CorrHHitPixVsEtapL5 ->Draw("colz");
CorrPixVsEtapL5->Draw("same");
TProfile* pr2 = CorrHHitPixVsEtapL5->ProfileX();
pr2->SetLineColor(2);
pr2->Draw("same");

c5->cd(4);                       
CorrHHitStripVsEtaP5->SetStats(kFALSE);
CorrHHitStripVsEtaP5 ->Draw("colz");
CorrStripVsEtaP5->Draw("same");
TProfile* pr3 = CorrHHitStripVsEtaP5->ProfileX();
pr3->SetLineColor(2);
pr3->Draw("same");

c5->cd(3);                       
CorrHHitStripVsEtapL5->SetStats(kFALSE);
CorrHHitStripVsEtapL5 ->Draw("colz");
CorrStripVsEtapL5->Draw("same");
TProfile* pr4 = CorrHHitStripVsEtapL5->ProfileX();
pr4->SetLineColor(2);
pr4->Draw("same");

c5->cd();
if (mb) c5->SaveAs("colorProfile_mb_corr.png");
else c5->SaveAs("colorProfile_tt_corr.png");
cout <<  "c5 ok"<< endl;

TCanvas *c5b = new TCanvas("c5b", "EtaDep UnCorr", 800,600);
c5b->Divide(2,2);                       
c5b->cd(2);                       
UnCorrHHitPixVsEtaP5->SetStats(kFALSE);
UnCorrHHitPixVsEtaP5 ->Draw("colz");
UnCorrPixVsEtaP5->Draw("same");
pr1 = UnCorrHHitPixVsEtaP5->ProfileX();
pr1->SetLineColor(2);
pr1->Draw("same");

c5b->cd(1);                       
UnCorrHHitPixVsEtapL5->SetStats(kFALSE);
UnCorrHHitPixVsEtapL5 ->Draw("colz");
UnCorrPixVsEtapL5->Draw("same");
pr2 = UnCorrHHitPixVsEtapL5->ProfileX();
pr2->SetLineColor(2);
pr2->Draw("same");

c5b->cd(4);                       
UnCorrHHitStripVsEtaP5->SetStats(kFALSE);
UnCorrHHitStripVsEtaP5 ->Draw("colz");
UnCorrStripVsEtaP5->Draw("same");
pr3 = UnCorrHHitStripVsEtaP5->ProfileX();
pr3->SetLineColor(2);
pr3->Draw("same");

c5b->cd(3);                       
UnCorrHHitStripVsEtapL5->SetStats(kFALSE);
UnCorrHHitStripVsEtapL5 ->Draw("colz");
UnCorrStripVsEtapL5->Draw("same");
pr4 = UnCorrHHitStripVsEtapL5->ProfileX();
pr4->SetLineColor(2);
pr4->Draw("same");

c5b->cd();
if (mb) c5b->SaveAs("colorProfile_mb_uncorr.png");
else c5b->SaveAs("colorProfile_tt_uncorr.png");
cout <<  "c5b ok"<< endl;


TCanvas *c6 = new TCanvas("c6", "Corr EtaDep", 600,600);
c6->cd();

CorrPixVsEtaP5->SetLineColor(2);
CorrPixVsEtapL5->SetLineColor(2);
CorrStripVsEtaP5->SetLineColor(4);
CorrStripVsEtapL5->SetLineColor(4);

CorrPixVsEtaP5->SetMarkerColor(2);
CorrPixVsEtapL5->SetMarkerColor(2);
CorrStripVsEtaP5->SetMarkerColor(4);
CorrStripVsEtapL5->SetMarkerColor(4);

CorrPixVsEtaP5->SetMarkerStyle(21);
CorrPixVsEtapL5->SetMarkerStyle(26);
CorrStripVsEtaP5->SetMarkerStyle(21);
CorrStripVsEtapL5->SetMarkerStyle(26);


CorrPixVsEtaP5->SetStats(kFALSE);
CorrPixVsEtaP5->SetMaximum(8);
CorrPixVsEtaP5->SetMinimum(1);
CorrPixVsEtaP5->Draw();
CorrPixVsEtapL5->Draw("same");
CorrStripVsEtaP5->Draw("same");
CorrStripVsEtapL5->Draw("same");

TLegend* qw3 =  new TLegend(0.15,0.12,0.40,0.28);
TString type;
if (mb) type=" (MB)";
else type=" (tt)";

qw3->AddEntry(CorrPixVsEtaP5,    "Pixel p>5"+type,                       "p");
qw3->AddEntry(CorrPixVsEtapL5,   "Pixel p<5"+type,                       "p");
qw3->AddEntry(CorrStripVsEtaP5,      "Strip p>5"+type,                       "p");
qw3->AddEntry(CorrStripVsEtapL5,     "Strip p<5"+type,                       "p");
qw3->Draw("same");

if (mb) c6->SaveAs("corr_etadep_mb.png");
else c6->SaveAs("corr_etadep_tt.png");

TCanvas *c6b = new TCanvas("c6b", "UnCorr EtaDep", 600,600);
c6b->cd();

UnCorrPixVsEtaP5->SetLineColor(2);
UnCorrPixVsEtapL5->SetLineColor(2);
UnCorrStripVsEtaP5->SetLineColor(4);
UnCorrStripVsEtapL5->SetLineColor(4);

UnCorrPixVsEtaP5->SetMarkerColor(2);
UnCorrPixVsEtapL5->SetMarkerColor(2);
UnCorrStripVsEtaP5->SetMarkerColor(4);
UnCorrStripVsEtapL5->SetMarkerColor(4);

UnCorrPixVsEtaP5->SetMarkerStyle(21);
UnCorrPixVsEtapL5->SetMarkerStyle(26);
UnCorrStripVsEtaP5->SetMarkerStyle(21);
UnCorrStripVsEtapL5->SetMarkerStyle(26);



UnCorrPixVsEtaP5->SetStats(kFALSE);
UnCorrPixVsEtaP5->SetMaximum(6);
UnCorrPixVsEtaP5->SetMinimum(1);
UnCorrPixVsEtaP5->Draw();
UnCorrPixVsEtapL5->Draw("same");
UnCorrStripVsEtaP5->Draw("same");
UnCorrStripVsEtapL5->Draw("same");

qw3->Draw("same");
if (mb) c6b->SaveAs("uncorr_etadep_mb.png");
else c6b->SaveAs("uncorr_etadep_tt.png");

TCanvas *c7 = new TCanvas("c7", "Corr PDep", 600,600);
c7->cd();

CorrPixVsP->SetLineColor(2);
CorrPixVsPExt->SetLineColor(2);
CorrStripVsP->SetLineColor(4);
CorrStripVsPExt->SetLineColor(4);

CorrPixVsP->SetMarkerColor(2);
CorrPixVsPExt->SetMarkerColor(2);
CorrStripVsP->SetMarkerColor(4);
CorrStripVsPExt->SetMarkerColor(4);

CorrPixVsP->SetMarkerStyle(21);
CorrPixVsPExt->SetMarkerStyle(26);
CorrStripVsP->SetMarkerStyle(21);
CorrStripVsPExt->SetMarkerStyle(26);

CorrPixVsP->SetStats(kFALSE);
CorrPixVsP->SetMaximum(6);
CorrPixVsP->SetMinimum(1);
CorrPixVsP->Draw();
CorrPixVsPExt->Draw("same");
CorrStripVsP->Draw("same");
CorrStripVsPExt->Draw("same");
TLegend* qw4 =  new TLegend(0.15,0.12,0.40,0.28);
qw4->AddEntry(CorrPixVsP,    "Pixel |#eta|<0.4"+type,                       "p");
qw4->AddEntry(CorrPixVsPExt,   "Pixel |#eta|>0.4"+type,                       "p");
qw4->AddEntry(CorrStripVsP,      "Strip |#eta|<0.4"+type,                       "p");
qw4->AddEntry(CorrStripVsPExt,     "Strip |#eta|>0.4"+type,                       "p");
qw4->Draw("same");
if (mb) c7->SaveAs("corr_pdepw_mb.png");
else c7->SaveAs("corr_pdepw_tt.png");

TCanvas *c7b = new TCanvas("c7b", "UnCorr PDep", 600,600);
c7b->cd();

UnCorrPixVsP->SetLineColor(2);
UnCorrPixVsPExt->SetLineColor(2);
UnCorrStripVsP->SetLineColor(4);
UnCorrStripVsPExt->SetLineColor(4);

UnCorrPixVsP->SetMarkerColor(2);
UnCorrPixVsPExt->SetMarkerColor(2);
UnCorrStripVsP->SetMarkerColor(4);
UnCorrStripVsPExt->SetMarkerColor(4);

UnCorrPixVsP->SetMarkerStyle(21);
UnCorrPixVsPExt->SetMarkerStyle(26);
UnCorrStripVsP->SetMarkerStyle(21);
UnCorrStripVsPExt->SetMarkerStyle(26);

UnCorrPixVsP->SetStats(kFALSE);
UnCorrPixVsP->SetMaximum(6);
UnCorrPixVsP->SetMinimum(1);
UnCorrPixVsP->Draw();
UnCorrPixVsPExt->Draw("same");
UnCorrStripVsP->Draw("same");
UnCorrStripVsPExt->Draw("same");
qw4->Draw("same");
if (mb) c7b->SaveAs("uncorr_pdepw_mb.png");
else c7b->SaveAs("uncorr_pdepw_tt.png");



}

void MakeMyProfile(TH2D* input, TH1D* FitResult) 
{
	       TH2D* inputnew = (TH2D*)input->Clone("tempTH2D");


	       for(int x=1;x<inputnew->GetXaxis()->GetNbins()+1;x++){
		  double eta = inputnew->GetXaxis()->GetBinCenter(x);
		  TH1D* Projection = (TH1D*)(inputnew->ProjectionY("proj",x,x))->Clone();
                  cout << " x " << x << " = "  << eta << " integral " << Projection->Integral() << endl;
		  if(Projection->Integral()<100) {
		   FitResult->SetBinContent(x, 0.);
		   FitResult->SetBinError  (x, 0.);
                  } 
                  else {
   		   Projection->SetAxisRange(1.,8.,"X");
		   Projection->Sumw2();
		   Projection->Scale(1.0/Projection->Integral());


		   TF1* mygaus = new TF1("mygaus","gaus", 5., 1.);
		   Projection->Fit("mygaus","Q0 RME");
		   double chiFromFit  = (mygaus->GetChisquare())/(mygaus->GetNDF());
		   FitResult->SetBinContent(x, mygaus->GetParameter(1));
		   FitResult->SetBinError  (x, mygaus->GetParError (1));
                   cout << "      fit " << mygaus->GetParameter(1) << "  chi2 " << chiFromFit << endl;
                   delete mygaus;

/*

                   float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
                   float liminf = maxval -2;
                   if (liminf<0) liminf=0;
                   float limsup = maxval + 3;
                   TF1* mylandau = new TF1("mylandau","[0]*TMath::Landau(x,[1],[2])",liminf, limsup);
                   mylandau->SetParameters(1,maxval,0.3);
                   mylandau->SetParLimits(0, 0.1, 1000000000.0);
                   mylandau->SetParLimits(1, 0.0001, 100.0);
                   mylandau->SetParLimits(2, 0.0001, 100.0);
                   Projection->Fit("mylandau","Q 0 RME");
                   FitResult->SetBinContent(x, mylandau->GetParameter(1));
                   FitResult->SetBinError  (x, mylandau->GetParError (1));
                   delete mylandau;

*/


                  }

                  delete Projection;
	       }

               delete inputnew;
               return FitResult;
}


//void GetSF(bool mb, bool corr) {
void GetSF() {
  TFile *_file0 ;
  TFile *_file1 ;

/*
  if (mb)  {
//   if(corr) _file0= TFile::Open("minbias_template_corr.root");
//   else _file0 = TFile::Open("minbias_template_uncorr.root");
//   if(corr) _file0= TFile::Open("minb3_template_corr.root");
//   else _file0 = TFile::Open("minb3_template_uncorr.root");
   if(corr) _file0= TFile::Open("mb_247_3_template_corr.root");
   else _file0 = TFile::Open("mb_247_3_template_uncorr.root");
  }
  else {
//   if(corr) _file0 = TFile::Open("ttbar_template_corr.root");
//   else _file0 = TFile::Open("ttbar_template_uncorr.root");
   if(corr) _file0 = TFile::Open("tt3_template_corr.root");
   else _file0 = TFile::Open("tt3_template_uncorr.root");
  }
*/
   TH2D*   HdedxVsP1;
   TH2D*   HdedxVsP2;

   bool mc_to_data=false;
   bool mc_alone=false;
   bool data_alone=false;
   bool test_year=true;

  if (mc_alone)  {
     _file0 = TFile::Open("crab_Analysis_2018_AllBackground_CodeV43p3_v1.root");
     _file0->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DPix",HdedxVsP2);
     _file0->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DStrip",HdedxVsP1);
  }
  else if (data_alone) {
//     _file0 = new TFile ("crab_Analysis_SingleMuon_Run2018_CodeV43p3_v1.root");  // 
     _file0 = new TFile ("crab_Analysis_SingleMuon_Run2017_CodeV43p3_v1.root");  // 
     _file0->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DPix",HdedxVsP2);
     _file0->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DStrip",HdedxVsP1);
  }
  else if (mc_to_data)  {
     _file0 = new TFile ("crab_Analysis_SingleMuon_Run2018_CodeV43p3_v1.root");
     _file0->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DStrip",HdedxVsP1);
     _file1 = TFile::Open("crab_Analysis_2018_AllBackground_CodeV43p3_v1.root");
     _file1->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DStrip",HdedxVsP2);
  }
  else if (test_year)  {
     _file0 = new TFile ("crab_Analysis_SingleMuon_Run2018_CodeV43p3_v1.root");
     _file0->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DStrip",HdedxVsP1);
     _file1 = TFile::Open("crab_Analysis_SingleMuon_Run2017_CodeV43p3_v1.root");
     _file1->GetObject("HSCParticleAnalyzer/BaseName/SF_HHit2DStrip",HdedxVsP2);
  }


   HdedxVsP2->Rebin2D(2,1);
   HdedxVsP1->Rebin2D(2,1);
    TCanvas* c0 = new TCanvas("c0", "c0", 600,600);
   c0->Divide(2,1);
   c0->cd(1);
   HdedxVsP2->Draw("colz");
   c0->cd(2);
   HdedxVsP1->Draw("colz");
   c0->cd();
   cout << " 2Dhisto ok " << endl;

   TH1D* HdedxVsPProfile1 = new TH1D("HdedxVsPProfile1"       , "HdedxVsPProfile1"      ,HdedxVsP1->GetXaxis()->GetNbins(),HdedxVsP1->GetXaxis()->GetXmin(),HdedxVsP1->GetXaxis()->GetXmax());
   TH1D* HdedxVsPProfile2 = new TH1D("HdedxVsPProfile2"       , "HdedxVsPProfile2"      ,HdedxVsP2->GetXaxis()->GetNbins(),HdedxVsP2->GetXaxis()->GetXmin(),HdedxVsP2->GetXaxis()->GetXmax());

   HdedxVsPProfile1->Rebin(2);
   HdedxVsPProfile2->Rebin(2);


   MakeMyProfile(HdedxVsP1, HdedxVsPProfile1);
   MakeMyProfile(HdedxVsP2, HdedxVsPProfile2);
   cout << " Profile ok " << endl;

   //

   double SFProfile;
   std::string SaveDir = "SFdir";
   FILE* fout = fopen ((SaveDir+"ScaleFactors.txt").c_str(), "w");
   fprintf (fout, "=================\n >>> Old CCC >>>\n=================\n");

   TH1D* Chi2Dist = new TH1D("Chi2Dist","Chi2Dist",5000, 0.8 ,1.8);

   double Minimum = 99999999999;
   double AbsGain = -1;

   for(int i=1;i<=Chi2Dist->GetNbinsX();i++){
      double ScaleFactor = Chi2Dist->GetXaxis()->GetBinCenter(i);
      TH1D* Rescaled = (TH1D*)HdedxVsPProfile2->Clone("Cloned");
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
   Chi2Dist->SetStats(kFALSE);
   Chi2Dist->GetXaxis()->SetNdivisions(504);
   Chi2Dist->GetXaxis()->SetTitle("Rescale Factor");
   Chi2Dist->GetYaxis()->SetTitle("Weighted Square Distance");
   Chi2Dist->Draw("");
   c1->SetLogy(true);
   c1->SetGridx(true); 
   SaveCanvas(c1, SaveDir, "Rescale_Chi2Dist");

   TCanvas* c2 = new TCanvas("c2", "c2", 600,600);
   TLegend* leg = new TLegend (0.25, 0.70, 0.65, 0.90);
   leg->SetHeader ("Fitting the Gaussian means");
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   HdedxVsPProfile1->SetStats(kFALSE);
   HdedxVsPProfile1->SetAxisRange(5,50,"X");
   HdedxVsPProfile1->SetAxisRange(0.0,10.,"Y");
   HdedxVsPProfile1->GetXaxis()->SetTitle("track momentum (GeV)");
   HdedxVsPProfile1->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
   HdedxVsPProfile1->SetMarkerColor(2);
   HdedxVsPProfile1->SetLineColor(2);
   HdedxVsPProfile1->SetLineWidth(2);
   HdedxVsPProfile1->Draw("");

   HdedxVsPProfile2->SetMarkerColor(1);
   HdedxVsPProfile2->Draw("same");
   TH1D* HdedxVsPProfile4 = (TH1D*)HdedxVsPProfile2->Clone("afs");
   HdedxVsPProfile4->SetMarkerColor(4);
   HdedxVsPProfile4->SetLineColor(4);
   HdedxVsPProfile4->Scale(AbsGain);
   HdedxVsPProfile4->Draw("same");
   if (data_alone || mc_alone) { 
   leg->AddEntry (HdedxVsPProfile1, "Strip", "LP");
   leg->AddEntry (HdedxVsPProfile2, "unscaled Pixel", "LP");
   leg->AddEntry (HdedxVsPProfile4, "scaled Pixel",   "LP");
   }
   else if (mc_to_data) {
   leg->AddEntry (HdedxVsPProfile1, "Data Strip ", "LP");
   leg->AddEntry (HdedxVsPProfile2, "unscaled MC Pixel", "LP");
   leg->AddEntry (HdedxVsPProfile4, "scaled MC Pixel",   "LP");
   }
   leg->Draw();

//   DrawPreliminary("", 13, "");
   SaveCanvas(c2, SaveDir, "Rescale_HdedxVsPProfile");

   if (data_alone || mc_alone) fprintf (fout, "HHit Pixel to Strip ::  Profile: %.7lf \n", SFProfile);
   else if (mc_to_data)  fprintf (fout, "HHit MC to Data Strip ::  Profile: %.7lf \n", SFProfile);
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
