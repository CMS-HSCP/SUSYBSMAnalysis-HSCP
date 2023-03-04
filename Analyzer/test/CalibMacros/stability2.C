#include "../tdrstyle.C"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TH1D.h"
#include "TMultiGraph.h"
#include "TGraph.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <iomanip>
#include "THStack.h"
#include "TFile.h"
#include "TROOT.h"
#include "TColor.h"
#include "TF1.h"
#include "TMath.h"
#include <TLegend.h>
#include <TProfile.h>
#include <TLatex.h>
#include <TAxis.h>
#include <TProfile.h>
#include <TString.h>
#include "Math/Vavilov.h"
#include "Math/VavilovAccurate.h"

bool check=false;
int which_check=1;

using namespace std;

struct Vavilov_Func {
   Vavilov_Func() {}

   double operator() (const double *x, const double *p) {
      double kappa = p[0];
      double beta2 = p[1];
      // float xvav = (qtotal-mpv)/sigmaQ;
      return p[4]*( pdf.Pdf( (x[0]-p[2])/p[3], kappa,beta2) );
   }
   
   ROOT::Math::VavilovAccurate pdf;
};

TH1D* extractMaximum(TH2D* histo2d);
TH1D* extractMaximum(TH2D* histo2d, TString namehist);
TH1D* extractRMS(TH2D* histo2d, TString namehist);
TH1D* extractLandau(TH2D* histo2d, TString namehist);
TH1D* extractGaus(TH2D* histo2d, TString namehist);
TH1D* extractVavilov(TH2D* histo2d, TString namehist);
TH1D* extractQuantile(TH2D* histo2d, TString namehist, int caseQ);
TH1D* extractVariationGaus(TH2D* histo2d, TString namehist, int nbin, float x1, float x2);
TH1D* extractVariationMaximum(TH2D* histo2d, TString namehist , int nbin, float x1, float x2);

TH1D* extractDiff(TFile* file1, TFile* file2, TString namehistin, TString namehistout);
TH2D* modif2DHisto(TH2D* histo, TString namehist);
TH2D* correl2DHisto(TH1D* histo1, TH1D* histo2, TString namehist, TString namecap);

void stability(){

int caseFit2=2;
bool full=false;   // Charge Pixel in 0-20 or in 0-4.5 
//TFile *myfile = TFile::Open("crab_Analysis_SingleMuon_CodeV42p6_v1.root"); // 12jan23
TFile *myfile = TFile::Open("crab_Analysis_SingleMuon_CodeV43p3_v1.root"); //16feb23

TCanvas *c1 = new TCanvas("Ih", "Ih",10,32,782,552);
    c1->SetFillColor(10);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

TH1D* LdEdXVsRun;
TH1D* LdEdXVsRunSel;
TH1D* LdEdXVsRun2;
TH1D* LdEdXVsRunSel2;

// GAUS 

TH2D* dEdX0NoL1VsRun;
myfile->GetObject("HSCParticleAnalyzer/BaseName/Stab_Ih_NoL1_VsRun",dEdX0NoL1VsRun);

 LdEdXVsRun = extractGaus(dEdX0NoL1VsRun,"LdEdXVsRun");

 LdEdXVsRun->SetMarkerColor(2);
 LdEdXVsRun->SetLineColor(2);

 LdEdXVsRun->SetMinimum(3.);
 LdEdXVsRun->SetMaximum(3.7);
 LdEdXVsRun->GetXaxis()->SetRangeUser(295000, 325500);
 LdEdXVsRun->GetXaxis()->SetTitle("Run number");
 LdEdXVsRun->GetYaxis()->SetTitle("#mu(Ih per bin)");
 LdEdXVsRun->Draw();

 c1->SaveAs("Stab_IhNoL1_Gauss.pdf");


TCanvas *c1_b = new TCanvas("Ias", "Ias",10,32,782,552);
    c1_b->SetFillColor(10);
    c1_b->SetBorderMode(0);
    c1_b->SetBorderSize(2);
    c1_b->SetFrameFillColor(0);
    c1_b->SetFrameBorderMode(0);
    c1_b->SetLeftMargin(0.15);
    c1_b->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);


TH2D* iasStripVsRun;
myfile->GetObject("HSCParticleAnalyzer/BaseName/Stab_Gi_strip_VsRun",iasStripVsRun);
TH1D* Lias = extractMaximum(iasStripVsRun,"Lias");

Lias->SetMarkerColor(2);
Lias->SetLineColor(2);

Lias->SetMinimum(0.02);
Lias->SetMaximum(0.025);
Lias->GetXaxis()->SetRangeUser(295000, 325500);
 Lias->GetXaxis()->SetTitle("Run number");
 
 Lias->GetYaxis()->SetTitleOffset(1.5);
 Lias->GetYaxis()->SetTitle("<G per bin>");
Lias->Draw();
c1_b->SaveAs("Stab_Ias_Mean.pdf");


//----
//



TCanvas *c1_c = new TCanvas("ProbQ", "ProbQ",10,32,782,552);
    c1_c->SetFillColor(10);
    c1_c->SetBorderMode(0);
    c1_c->SetBorderSize(2);
    c1_c->SetFrameFillColor(0);
    c1_c->SetFrameBorderMode(0);
    c1_c->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

TH2D* InvprobQNoL1VsRunCut;
myfile->GetObject("HSCParticleAnalyzer/BaseName/Stab_Fi_pixNoL1_VsRun",InvprobQNoL1VsRunCut);
TH1D* Lprobq = extractMaximum(InvprobQNoL1VsRunCut,"Lprobq");


Lprobq->SetMarkerColor(2);
Lprobq->SetLineColor(2);



Lprobq->SetMinimum(0.6);
Lprobq->SetMaximum(0.7);
Lprobq->GetXaxis()->SetRangeUser(295000, 325500);
 Lprobq->GetXaxis()->SetTitle("Run number");
 Lprobq->GetYaxis()->SetTitleOffset(1.2);
 Lprobq->GetYaxis()->SetTitle("<F per bin>");
Lprobq->Draw();
c1_c->SaveAs("Stab_ProbQ_Mean.pdf");

TCanvas *c2 = new TCanvas("InvB", "InvB",10,32,782,552);
c2->SetFillColor(10);
c2->SetBorderMode(0);
c2->SetBorderSize(2);
c2->SetFrameFillColor(0);
c2->SetFrameBorderMode(0);
c2->cd();
gStyle->SetOptStat(0);


TH2D* invBVsRun;
myfile->GetObject("HSCParticleAnalyzer/BaseName/Stab_invB_VsRun",invBVsRun);
TH1D* LinvBVsRun = extractGaus(invBVsRun,"LinvBVsRun");

LinvBVsRun->SetMarkerColor(2);
LinvBVsRun->SetLineColor(2);

LinvBVsRun->SetMinimum(0.5);
LinvBVsRun->SetMaximum(1.5);
LinvBVsRun->GetXaxis()->SetRangeUser(295000, 325500);
LinvBVsRun->Draw();

//c2->SaveAs("Stab_invBeta_Gauss.pdf");


TCanvas *c3 = new TCanvas("Ih1D", "Ih",10,32,782,552);
    c3->SetFillColor(10);
    c3->SetBorderMode(0);
    c3->SetBorderSize(2);
    c3->SetFrameFillColor(0);
    c3->SetFrameBorderMode(0);
    c3->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1);

TH1D* L1DdEdXVsRun;
TH1D* L1DdEdXVsRunSel;

// GAUS 

 L1DdEdXVsRun = extractVariationGaus(dEdX0NoL1VsRun,"L1DdEdXVsRun",100,3.1,3.4);

 L1DdEdXVsRun->SetMarkerColor(2);
 L1DdEdXVsRun->SetLineColor(2);
 L1DdEdXVsRun->GetXaxis()->SetTitle("weigthed #mu(Ih)");

 L1DdEdXVsRun->Draw("hist");

 c3->SaveAs("Stab1D_IhNoL1_Gauss.pdf");




//
//---


TCanvas *c3_b = new TCanvas("Ias1D", "Ias",10,32,782,552);
    c3_b->SetFillColor(10);
    c3_b->SetBorderMode(0);
    c3_b->SetBorderSize(2);
    c3_b->SetFrameFillColor(0);
    c3_b->SetFrameBorderMode(0);
    c3_b->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1);

TH1D* L1Dias = extractVariationMaximum(iasStripVsRun,"L1Dias",100,0.02,0.025);


L1Dias->SetMarkerColor(2);
L1Dias->SetLineColor(2);



L1Dias->GetXaxis()->SetTitle("weigthed <G>");
L1Dias->Draw("hist");
c3_b->SaveAs("Stab1D_Ias_Mean.pdf");


//----
//



TCanvas *c3_c = new TCanvas("ProbQ1D", "ProbQ",10,32,782,552);
    c3_c->SetFillColor(10);
    c3_c->SetBorderMode(0);
    c3_c->SetBorderSize(2);
    c3_c->SetFrameFillColor(0);
    c3_c->SetFrameBorderMode(0);
    c3_c->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(1);

TH1D* L1Dprobq = extractVariationMaximum(InvprobQNoL1VsRunCut,"L1Dprobq",100,0.6, 0.7);


L1Dprobq->SetMarkerColor(2);
L1Dprobq->SetLineColor(2);



L1Dprobq->GetXaxis()->SetTitle("weigthed <F>");
L1Dprobq->Draw("hist" );
c3_c->SaveAs("Stab1D_ProbQ_Mean.pdf");

}

void plotDiff(){

TFile * file1 = TFile::Open("/opt/sbg/cms/safe1/cms/ccollard/HSCP/CMSSW_10_6_2/src/stage/ntuple/test/analysisRun2/analysis_full.root");
TFile * file2 = TFile::Open("analysis_ul.root");

TCanvas *c1 = new TCanvas("Ih", "Ih",10,32,782,552);
    c1->SetFillColor(10);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

TH1D* LdEdXVsRun = extractDiff(file1, file2,  "dEdXVsRun", "LdEdXVsRun");
TH1D* LdEdXpixVsRun = extractDiff(file1, file2,  "dEdXpixVsRun","LdEdXpixVsRun");
TH1D* LdEdXstripVsRun = extractDiff(file1, file2,  "dEdXstripVsRun","LdEdXstripVsRun");

LdEdXpixVsRun->SetMarkerColor(2);
LdEdXpixVsRun->SetLineColor(2);
LdEdXstripVsRun->SetMarkerColor(8);
LdEdXstripVsRun->SetLineColor(8);

LdEdXVsRun->SetMinimum(-1.);
LdEdXVsRun->SetMaximum(1.);
LdEdXVsRun->Draw();
LdEdXpixVsRun->Draw("same");
LdEdXstripVsRun->Draw("same");

c1->SaveAs("Ih.pdf");
//
//----
//

TCanvas *c0 = new TCanvas("Ih_0", "Ih_0",10,32,782,552);
c0->SetFillColor(10);
c0->SetBorderMode(0);
c0->SetBorderSize(2);
c0->SetFrameFillColor(0);
c0->SetFrameBorderMode(0);
c0->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LdEdX0VsRun = extractDiff(file1, file2,  "dEdX0VsRun", "LdEdX0VsRun");
TH1D* LdEdX0pixVsRun = extractDiff(file1, file2,  "dEdX0pixVsRun","LdEdX0pixVsRun");
TH1D* LdEdX0stripVsRun = extractDiff(file1, file2,  "dEdX0stripVsRun","LdEdX0stripVsRun");

LdEdX0pixVsRun->SetMarkerColor(2);
LdEdX0pixVsRun->SetLineColor(2);
LdEdX0stripVsRun->SetMarkerColor(8);
LdEdX0stripVsRun->SetLineColor(8);

LdEdX0VsRun->SetMinimum(-1.);
LdEdX0VsRun->SetMaximum(1.);
LdEdX0VsRun->Draw();
LdEdX0pixVsRun->Draw("same");
LdEdX0stripVsRun->Draw("same");
c0->SaveAs("Ih_0.pdf");

TCanvas *c11 = new TCanvas("Ih_4", "Ih_4",10,32,782,552);
    c11->SetFillColor(10);
    c11->SetBorderMode(0);
    c11->SetBorderSize(2);
    c11->SetFrameFillColor(0);
    c11->SetFrameBorderMode(0);
    c11->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

TH1D* LdEdX4VsRun = extractDiff(file1, file2,"dEdXV4sRun","LdEdX4VsRun");
TH1D* LdEdX4pixVsRun = extractDiff(file1, file2,"dEdX4pixVsRun","LdEdX4pixVsRun");
TH1D* LdEdX4stripVsRun = extractDiff(file1, file2,"dEdX4stripVsRun","LdEdX4stripVsRun");

LdEdX4pixVsRun->SetMarkerColor(2);
LdEdX4pixVsRun->SetLineColor(2);
LdEdX4stripVsRun->SetMarkerColor(8);
LdEdX4stripVsRun->SetLineColor(8);

LdEdX4VsRun->SetMinimum(-1.);
LdEdX4VsRun->SetMaximum(1.);
LdEdX4VsRun->Draw();
LdEdX4pixVsRun->Draw("same");
LdEdX4stripVsRun->Draw("same");

//
//----
//

TCanvas *c01 = new TCanvas("Ih_40", "Ih_40",10,32,782,552);
c01->SetFillColor(10);
c01->SetBorderMode(0);
c01->SetBorderSize(2);
c01->SetFrameFillColor(0);
c01->SetFrameBorderMode(0);
c01->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LdEdX40VsRun = extractDiff(file1, file2,"dEdX40VsRun","LdEdX40VsRun");
TH1D* LdEdX40pixVsRun = extractDiff(file1, file2,"dEdX40pixVsRun","LdEdX40pixVsRun");
TH1D* LdEdX40stripVsRun = extractDiff(file1, file2,"dEdX40stripVsRun","LdEdX40stripVsRun");

LdEdX40pixVsRun->SetMarkerColor(2);
LdEdX40pixVsRun->SetLineColor(2);
LdEdX40stripVsRun->SetMarkerColor(8);
LdEdX40stripVsRun->SetLineColor(8);

LdEdX40VsRun->SetMinimum(-1.);
LdEdX40VsRun->SetMaximum(1.);
LdEdX40VsRun->Draw();
LdEdX40pixVsRun->Draw("same");
LdEdX40stripVsRun->Draw("same");

//
//----
//

TCanvas *c2 = new TCanvas("InvB", "InvB",10,32,782,552);
c2->SetFillColor(10);
c2->SetBorderMode(0);
c2->SetBorderSize(2);
c2->SetFrameFillColor(0);
c2->SetFrameBorderMode(0);
c2->cd();
gStyle->SetOptStat(0);


TH1D* LinvBVsRun = extractDiff(file1, file2,"invBVsRun","LinvBVsRun");
TH1D* LinvBDTVsRun = extractDiff(file1, file2,"invBDTVsRun","LinvBDTVsRun");
TH1D* LinvBCSCVsRun = extractDiff(file1, file2,"invBCSCVsRun","LinvBCSCVsRun");

LinvBDTVsRun->SetMarkerColor(2);
LinvBDTVsRun->SetLineColor(2);
LinvBCSCVsRun->SetMarkerColor(8);
LinvBCSCVsRun->SetLineColor(8);

LinvBVsRun->SetMinimum(-1.);
LinvBVsRun->SetMaximum(1.);
LinvBVsRun->Draw();
LinvBDTVsRun->Draw("same");
LinvBCSCVsRun->Draw("same");


//
//----
//

//
//----
//

TCanvas *cChargPix = new TCanvas("ChargPix", "ChargPix",10,32,782,552);
cChargPix->SetFillColor(10);
cChargPix->SetBorderMode(0);
cChargPix->SetBorderSize(2);
cChargPix->SetFrameFillColor(0);
cChargPix->SetFrameBorderMode(0);
cChargPix->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LZooChargeVsRun_pixl1 = extractDiff(file1, file2,"ZooChargeVsRun_pixl1","LZooChargeVsRun_pixl1");
TH1D* LZooChargeVsRun_pixl2 = extractDiff(file1, file2,"ZooChargeVsRun_pixl2","LZooChargeVsRun_pixl2");
TH1D* LZooChargeVsRun_pixl3 = extractDiff(file1, file2,"ZooChargeVsRun_pixl3","LZooChargeVsRun_pixl3");
TH1D* LZooChargeVsRun_pixl4 = extractDiff(file1, file2,"ZooChargeVsRun_pixl4","LZooChargeVsRun_pixl4");

LZooChargeVsRun_pixl2->SetMarkerColor(2);
LZooChargeVsRun_pixl2->SetLineColor(2);
LZooChargeVsRun_pixl3->SetMarkerColor(8);
LZooChargeVsRun_pixl3->SetLineColor(8);
LZooChargeVsRun_pixl4->SetMarkerColor(4);
LZooChargeVsRun_pixl4->SetLineColor(4);

LZooChargeVsRun_pixl1->SetMinimum(-1.);
LZooChargeVsRun_pixl1->SetMaximum(1.);
LZooChargeVsRun_pixl1->Draw();
LZooChargeVsRun_pixl2->Draw("same");
LZooChargeVsRun_pixl3->Draw("same");
LZooChargeVsRun_pixl4->Draw("same");

cChargPix->SaveAs("ChargPix.pdf");
//
//--
//

TCanvas *cChargDisk = new TCanvas("ChargDisk", "ChargDisk",10,32,782,552);
cChargDisk->SetFillColor(10);
cChargDisk->SetBorderMode(0);
cChargDisk->SetBorderSize(2);
cChargDisk->SetFrameFillColor(0);
cChargDisk->SetFrameBorderMode(0);
cChargDisk->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LZooChargeVsRun_pixd1 = extractDiff(file1, file2,"ZooChargeVsRun_pixd1","LZooChargeVsRun_pixd1");
TH1D* LZooChargeVsRun_pixd2 = extractDiff(file1, file2,"ZooChargeVsRun_pixd2","LZooChargeVsRun_pixd2");
TH1D* LZooChargeVsRun_pixd3 = extractDiff(file1, file2,"ZooChargeVsRun_pixd3","LZooChargeVsRun_pixd3");

LZooChargeVsRun_pixd2->SetMarkerColor(2);
LZooChargeVsRun_pixd2->SetLineColor(2);
LZooChargeVsRun_pixd3->SetMarkerColor(8);
LZooChargeVsRun_pixd3->SetLineColor(8);

LZooChargeVsRun_pixd1->SetMinimum(-1.);
LZooChargeVsRun_pixd1->SetMaximum(1.);
LZooChargeVsRun_pixd1->Draw();
LZooChargeVsRun_pixd2->Draw("same");
LZooChargeVsRun_pixd3->Draw("same");
cChargDisk->SaveAs("ChargDisk.pdf");

//
//--
//

TCanvas *cChargTib = new TCanvas("ChargTib", "ChargTib",10,32,782,552);
cChargTib->SetFillColor(10);
cChargTib->SetBorderMode(0);
cChargTib->SetBorderSize(2);
cChargTib->SetFrameFillColor(0);
cChargTib->SetFrameBorderMode(0);
cChargTib->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LZooChargeVsRun_tib1 = extractDiff(file1, file2,"ZooChargeVsRun_tib1","LZooChargeVsRun_tib1");
TH1D* LZooChargeVsRun_tib2 = extractDiff(file1, file2,"ZooChargeVsRun_tib2","LZooChargeVsRun_tib2");
TH1D* LZooChargeVsRun_tib3 = extractDiff(file1, file2,"ZooChargeVsRun_tib3","LZooChargeVsRun_tib3");
TH1D* LZooChargeVsRun_tib4 = extractDiff(file1, file2,"ZooChargeVsRun_tib4","LZooChargeVsRun_tib4");

LZooChargeVsRun_tib2->SetMarkerColor(2);
LZooChargeVsRun_tib2->SetLineColor(2);
LZooChargeVsRun_tib3->SetMarkerColor(8);
LZooChargeVsRun_tib3->SetLineColor(8);
LZooChargeVsRun_tib4->SetMarkerColor(4);
LZooChargeVsRun_tib4->SetLineColor(4);

LZooChargeVsRun_tib1->SetMinimum(-1.);
LZooChargeVsRun_tib1->SetMaximum(1.);
LZooChargeVsRun_tib1->Draw();
LZooChargeVsRun_tib2->Draw("same");
LZooChargeVsRun_tib3->Draw("same");
LZooChargeVsRun_tib4->Draw("same");

cChargTib->SaveAs("ChargTib.pdf");
//
//--
//

TCanvas *cChargTob = new TCanvas("ChargTob", "ChargTob",10,32,782,552);
cChargTob->SetFillColor(10);
cChargTob->SetBorderMode(0);
cChargTob->SetBorderSize(2);
cChargTob->SetFrameFillColor(0);
cChargTob->SetFrameBorderMode(0);
cChargTob->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LZooChargeVsRun_tob1 = extractDiff(file1, file2,"ZooChargeVsRun_tob1","LZooChargeVsRun_tob1");
TH1D* LZooChargeVsRun_tob2 = extractDiff(file1, file2,"ZooChargeVsRun_tob2","LZooChargeVsRun_tob2");

LZooChargeVsRun_tob2->SetMarkerColor(2);
LZooChargeVsRun_tob2->SetLineColor(2);

LZooChargeVsRun_tob1->SetMinimum(-1.);
LZooChargeVsRun_tob1->SetMaximum(1.);
LZooChargeVsRun_tob1->Draw();
LZooChargeVsRun_tob2->Draw("same");
cChargTob->SaveAs("ChargTob.pdf");

//
//--
//

TCanvas *cChargTid = new TCanvas("ChargTid", "ChargTid",10,32,782,552);
cChargTid->SetFillColor(10);
cChargTid->SetBorderMode(0);
cChargTid->SetBorderSize(2);
cChargTid->SetFrameFillColor(0);
cChargTid->SetFrameBorderMode(0);
cChargTid->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LZooChargeVsRun_tid1 = extractDiff(file1, file2,"ZooChargeVsRun_tid1","LZooChargeVsRun_tid1");
TH1D* LZooChargeVsRun_tid2 = extractDiff(file1, file2,"ZooChargeVsRun_tid2","LZooChargeVsRun_tid2");
TH1D* LZooChargeVsRun_tid3 = extractDiff(file1, file2,"ZooChargeVsRun_tid3","LZooChargeVsRun_tid3");
TH1D* LZooChargeVsRun_tid4 = extractDiff(file1, file2,"ZooChargeVsRun_tid4","LZooChargeVsRun_tid4");
TH1D* LZooChargeVsRun_tid5 = extractDiff(file1, file2,"ZooChargeVsRun_tid5","LZooChargeVsRun_tid5");
TH1D* LZooChargeVsRun_tid6 = extractDiff(file1, file2,"ZooChargeVsRun_tid6","LZooChargeVsRun_tid6");

LZooChargeVsRun_tid2->SetMarkerColor(2);
LZooChargeVsRun_tid2->SetLineColor(2);
LZooChargeVsRun_tid3->SetMarkerColor(8);
LZooChargeVsRun_tid3->SetLineColor(8);
LZooChargeVsRun_tid4->SetMarkerColor(4);
LZooChargeVsRun_tid4->SetLineColor(4);
LZooChargeVsRun_tid5->SetMarkerColor(6);
LZooChargeVsRun_tid5->SetLineColor(6);
LZooChargeVsRun_tid6->SetMarkerColor(13);
LZooChargeVsRun_tid6->SetLineColor(13);

LZooChargeVsRun_tid1->SetMinimum(-1.);
LZooChargeVsRun_tid1->SetMaximum(1.);
LZooChargeVsRun_tid1->Draw();
LZooChargeVsRun_tid2->Draw("same");
LZooChargeVsRun_tid3->Draw("same");
LZooChargeVsRun_tid4->Draw("same");
LZooChargeVsRun_tid5->Draw("same");
LZooChargeVsRun_tid6->Draw("same");
cChargTid->SaveAs("ChargTid.pdf");


//
//--
//

TCanvas *cChargTec = new TCanvas("ChargTec", "ChargTec",10,32,782,552);
cChargTec->SetFillColor(10);
cChargTec->SetBorderMode(0);
cChargTec->SetBorderSize(2);
cChargTec->SetFrameFillColor(0);
cChargTec->SetFrameBorderMode(0);
cChargTec->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LZooChargeVsRun_tec1 = extractDiff(file1, file2,"ZooChargeVsRun_tec1","LZooChargeVsRun_tec1");
TH1D* LZooChargeVsRun_tec2 = extractDiff(file1, file2,"ZooChargeVsRun_tec2","LZooChargeVsRun_tec2");
TH1D* LZooChargeVsRun_tec3 = extractDiff(file1, file2,"ZooChargeVsRun_tec3","LZooChargeVsRun_tec3");
TH1D* LZooChargeVsRun_tec4 = extractDiff(file1, file2,"ZooChargeVsRun_tec4","LZooChargeVsRun_tec4");
TH1D* LZooChargeVsRun_tec5 = extractDiff(file1, file2,"ZooChargeVsRun_tec5","LZooChargeVsRun_tec5");
TH1D* LZooChargeVsRun_tec6 = extractDiff(file1, file2,"ZooChargeVsRun_tec6","LZooChargeVsRun_tec6");
TH1D* LZooChargeVsRun_tec7 = extractDiff(file1, file2,"ZooChargeVsRun_tec7","LZooChargeVsRun_tec7");
TH1D* LZooChargeVsRun_tec8 = extractDiff(file1, file2,"ZooChargeVsRun_tec8","LZooChargeVsRun_tec8");
TH1D* LZooChargeVsRun_tec9 = extractDiff(file1, file2,"ZooChargeVsRun_tec9","LZooChargeVsRun_tec9");

LZooChargeVsRun_tec2->SetMarkerColor(2);
LZooChargeVsRun_tec2->SetLineColor(2);
LZooChargeVsRun_tec3->SetMarkerColor(8);
LZooChargeVsRun_tec3->SetLineColor(8);
LZooChargeVsRun_tec4->SetMarkerColor(4);
LZooChargeVsRun_tec4->SetLineColor(4);
LZooChargeVsRun_tec5->SetMarkerColor(6);
LZooChargeVsRun_tec5->SetLineColor(6);
LZooChargeVsRun_tec6->SetMarkerColor(13);
LZooChargeVsRun_tec6->SetLineColor(13);
LZooChargeVsRun_tec7->SetMarkerColor(7);
LZooChargeVsRun_tec7->SetLineColor(7);
LZooChargeVsRun_tec8->SetMarkerColor(kOrange);
LZooChargeVsRun_tec8->SetLineColor(kOrange);
LZooChargeVsRun_tec9->SetMarkerColor(kYellow+2);
LZooChargeVsRun_tec9->SetLineColor(kYellow+2);

LZooChargeVsRun_tec1->SetMinimum(-1);
LZooChargeVsRun_tec1->SetMaximum(1.);
LZooChargeVsRun_tec1->Draw();
LZooChargeVsRun_tec2->Draw("same");
LZooChargeVsRun_tec3->Draw("same");
LZooChargeVsRun_tec4->Draw("same");
LZooChargeVsRun_tec5->Draw("same");
LZooChargeVsRun_tec6->Draw("same");
LZooChargeVsRun_tec7->Draw("same");
LZooChargeVsRun_tec8->Draw("same");
LZooChargeVsRun_tec9->Draw("same");
cChargTec->SaveAs("ChargTec.pdf");

//
//----
//

}

TH1D* extractMaximum(TH2D* histo2d,  TString namehist){

   TH1D* histo_modif = new TH1D(namehist,namehist, histo2d->GetXaxis()->GetNbins(),histo2d->GetXaxis()->GetXmin(),histo2d->GetXaxis()->GetXmax());
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      
      if (Projection->Integral()>100) {
//       float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
//       float rmsval=Projection->GetXaxis()->GetBinWidth(Projection->GetMaximumBin());
       float maxval=Projection->GetMean();
       float rmsval=Projection->GetMeanError();
       histo_modif->SetBinContent(x,maxval);
       histo_modif->SetBinError(x,rmsval);
      }
   }
   return histo_modif;
}
TH1D* extractRMS(TH2D* histo2d,  TString namehist){

   TH1D* histo_modif = new TH1D(namehist,namehist, histo2d->GetXaxis()->GetNbins(),histo2d->GetXaxis()->GetXmin(),histo2d->GetXaxis()->GetXmax());
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      
      if (Projection->Integral()>100) {
       float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
//       float rmsval=Projection->GetXaxis()->GetBinWidth(Projection->GetMaximumBin());
       float rmsval=Projection->GetRMS();
       float errrmsval=Projection->GetRMSError();
       histo_modif->SetBinContent(x,rmsval);
       histo_modif->SetBinError(x,errrmsval);
      }
   }
   return histo_modif;
}

TH1D* extractLandau(TH2D* histo2d, TString namehist){

   int debug=0;
   int debug2=0;
   TH1D* histo_modif = new TH1D(namehist,namehist, histo2d->GetXaxis()->GetNbins(),histo2d->GetXaxis()->GetXmin(),histo2d->GetXaxis()->GetXmax());
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      
      if (Projection->Integral()>1000) {
       debug++;
       float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
       float liminf = maxval -2;
       if (liminf<0) liminf=0;
       float limsup = maxval + 3;
       TF1* mylandau = new TF1("mylandau","[0]*TMath::Landau(x,[1],[2])",liminf, limsup);
       mylandau->SetParameters(1,maxval,0.3); 
       //int fstatus = Projection->Fit("mylandau","Q 0 RME");
       mylandau->SetParLimits(0, 0.1, 1000000000.0);
       mylandau->SetParLimits(1, 0.0001, 100.0);
       mylandau->SetParLimits(2, 0.0001, 100.0);

       TFitResultPtr r;
       if (debug<-1) {
          r = Projection->Fit("mylandau","S");
          cout << " status :  " << r->Status() << "  validity " << r->IsValid() << endl;
       }
       else {r = Projection->Fit("mylandau","Q 0 RME S"); }

//       if (!check && mylandau->GetParError(1)<0.1) {
       if (!check) {
         if (r->IsValid() && (r->Status()==0 || r->Status()==4000)) {
          histo_modif->SetBinContent(x,mylandau->GetParameter(1));
          histo_modif->SetBinError(x,mylandau->GetParError(1));
         }
       }
       else if (which_check==1) {
        int ndf=mylandau->GetNDF();
        if (ndf<1) ndf=1;
        histo_modif->SetBinContent(x,mylandau->GetChisquare()/ndf);
        histo_modif->SetBinError(x,0.05*mylandau->GetChisquare()/ndf);
       }
       else if (which_check==2) {
        if (r->IsValid()) {
         histo_modif->SetBinContent(x,r->Status());
//        histo_modif->SetBinContent(x,r->IsValid());
         histo_modif->SetBinError(x,0.1);
         if (r->Status()==4000 && debug2<10) {
            debug2++;
            cout << " ***** " << endl;
            Projection->Fit("mylandau","S");
            cout << " status :  " << r->Status() << "  validity " << r->IsValid() << endl;
            cout << " ***** " << endl;
         }
        }
      }
     }
   }
   return histo_modif;
}
TH1D* extractGaus(TH2D* histo2d, TString namehist){

   TH1D* histo_modif = new TH1D(namehist,namehist, histo2d->GetXaxis()->GetNbins(),histo2d->GetXaxis()->GetXmin(),histo2d->GetXaxis()->GetXmax());
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      
      if (Projection->Integral()>1000) {
       TF1* mygaus = new TF1("mygaus","[0]*TMath::Gaus(x,[1],[2])", 1, 6);
       float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
       mygaus->SetParameters(1,maxval,0.3); 
       mygaus->SetParLimits(1, 0., 10.0);
       mygaus->SetParLimits(2, 0., 10.0);
//       int fstatus = Projection->Fit("mygaus","Q 0 RME");
       TFitResultPtr r = Projection->Fit("mygaus","Q 0 RME S"); 
//       if (!check && mygaus->GetParError(1)<0.3) {
       if (!check) {
        if(r->IsValid() && (r->Status()==0 || r->Status()==4000) && mygaus->GetParError(1)<0.1) {
         histo_modif->SetBinContent(x,mygaus->GetParameter(1));
         histo_modif->SetBinError(x,mygaus->GetParError(1));
/*
         if( mygaus->GetParameter(1)>0) {
           histo_modif->SetBinContent(x,mygaus->GetParameter(2)/mygaus->GetParameter(1));
           histo_modif->SetBinError(x,mygaus->GetParError(2)/mygaus->GetParameter(1));
         }
*/
        }

       }
       else if (which_check==1) {
        int ndf=mygaus->GetNDF();
        if (ndf<1) ndf=1;
        histo_modif->SetBinContent(x,mygaus->GetChisquare()/ndf);
        histo_modif->SetBinError(x,0.05*mygaus->GetChisquare()/ndf);
       }
       else if (which_check==2) {
        histo_modif->SetBinContent(x,r->Status());
        histo_modif->SetBinError(x,0.1);
       }
      }
   }
   return histo_modif;
}

TH1D* extractVavilov(TH2D* histo2d, TString namehist){

   TH1D* histo_modif = new TH1D(namehist,namehist, histo2d->GetXaxis()->GetNbins(),histo2d->GetXaxis()->GetXmin(),histo2d->GetXaxis()->GetXmax());
   int debug=0;
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();
      float maxval2=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
      if (Projection->Integral()>1000) {
       debug++;
       Vavilov_Func* VVf=new Vavilov_Func();
       TF1* TamasFunc  = new TF1("TamasFunc",VVf, maxval2-0.5,8,5,"Vavilov_Func");
       TamasFunc->SetParameters(
            0.01*maxval2,   // p[0]: kappa [0.01:10]
            1.0,      
            0.8*maxval2,        // p[2]: mpv
            0.2*Projection->GetRMS(),         // p[3]: sigmaQ
            0.5*Projection->Integral()      // p[4]: norm
       );
       TamasFunc->SetParLimits(0, 0.001, 10.0);
       TamasFunc->SetParLimits(1, 0.9, 1.);
       TamasFunc->SetParLimits(2, 0.5, 10.0);
       TamasFunc->SetParLimits(3, 0.0, 100.0);
       TamasFunc->SetParLimits(4, 0.0, 200000000.0);
//       int fstatus = Projection->Fit(TamasFunc,"Q 0 RME");
//       if (!check && TamasFunc->GetParError(2)<0.3) {
//       TFitResultPtr r = Projection->Fit("TamasFunc","Q 0 RME S"); 
       TFitResultPtr r;
       if (debug<10) {
          r = Projection->Fit("TamasFunc","S");
          cout << " status :  " << r->Status() << "  validity " << r->IsValid() << endl;
       }
       else { r = Projection->Fit("TamasFunc","Q 0 RME S"); }

       if (!check) {
        if (r->IsValid() && (r->Status()==0 || r->Status()==4000)) {
          histo_modif->SetBinContent(x,TamasFunc->GetParameter(2));
          histo_modif->SetBinError(x,TamasFunc->GetParError(2));
        }
       }
       else if (which_check==1) {
        int ndf=TamasFunc->GetNDF();
        if (ndf<1) ndf=1;
        histo_modif->SetBinContent(x,TamasFunc->GetChisquare()/ndf);
        histo_modif->SetBinError(x,0.05*TamasFunc->GetChisquare()/ndf);
       }
       else if (which_check==2) {
//        histo_modif->SetBinContent(x,fstatus);
        histo_modif->SetBinContent(x,r->Status());
        histo_modif->SetBinError(x,0.1);
       }  
      }
   }
   return histo_modif;
}


TH1D* extractDiff(TFile* file1, TFile* file2, TString namehistin, TString namehistout){
TH1D* L1 ;
TH1D* L2 ;
TH1D* Lout ;
file1->cd();
TString L1name = namehistin+"L1";
L1 = extractGaus((TH2D*)gROOT->FindObject(namehistin),L1name);
file2->cd();
TString L2name = namehistin+"L2";
L2 = extractGaus((TH2D*)gROOT->FindObject(namehistin),L2name);

Lout = (TH1D*) L1->Clone();
Lout->SetName(namehistout);
Lout->Add(L2,-1.);

return Lout;
}


void IntegratedResolution(){

 TFile *_file0 = TFile::Open("analysis_ul_2017_tamas_v2.root");
 TCanvas *c1 = new TCanvas("I", "I",10,32,782,552);
    c1->SetFillColor(10);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

    TH2D* IhNoL1 = (TH2D*)gROOT->FindObject("dEdX0NoL1VsRun");
    TH2D* Ih15drop = (TH2D*)gROOT->FindObject("dEdXVsRun");
    TH2D* IhStrip = (TH2D*)gROOT->FindObject("dEdX0stripVsRun");
    TH2D* Ih15Strip = (TH2D*)gROOT->FindObject("dEdXstripVsRun");

    TH1D* ProjectionIhNoL1 = (TH1D*)(IhNoL1->ProjectionY("proj"))->Clone();   
    TH1D* ProjectionIh15drop = (TH1D*)(Ih15drop->ProjectionY("proj"))->Clone();   
    TH1D* ProjectionIhStrip = (TH1D*)(IhStrip->ProjectionY("proj"))->Clone();   
    TH1D* ProjectionIh15Strip = (TH1D*)(Ih15Strip->ProjectionY("proj"))->Clone();   
      
    TF1* mygaus1 = new TF1("mygaus1","[0]*TMath::Gaus(x,[1],[2])", 1, 6);
    float maxval=ProjectionIhNoL1->GetXaxis()->GetBinCenter(ProjectionIhNoL1->GetMaximumBin());
    mygaus1->SetParameters(1,maxval,0.3); 
    ProjectionIhNoL1->Fit("mygaus1"); 
    if ((mygaus1->GetParameter(1))>0) cout << " Resolution ProjectionIhNoL1 " << (mygaus1->GetParameter(2))/(mygaus1->GetParameter(1)) << endl;
    cout << endl;

    TF1* mygaus2 = new TF1("mygaus2","[0]*TMath::Gaus(x,[1],[2])", 1, 6);
    float maxval2=ProjectionIh15drop->GetXaxis()->GetBinCenter(ProjectionIh15drop->GetMaximumBin());
    mygaus2->SetParameters(1,maxval2,0.3); 
    ProjectionIh15drop->Fit("mygaus2"); 
    if ((mygaus2->GetParameter(1))>0) cout << " Resolution ProjectionIh15drop " << (mygaus2->GetParameter(2))/ (mygaus2->GetParameter(1)) << endl;

    TF1* mygaus3 = new TF1("mygaus3","[0]*TMath::Gaus(x,[1],[2])", 1, 6);
    float maxval3=ProjectionIhStrip->GetXaxis()->GetBinCenter(ProjectionIhStrip->GetMaximumBin());
    mygaus3->SetParameters(1,maxval3,0.3); 
    ProjectionIhStrip->Fit("mygaus3"); 
    if ((mygaus3->GetParameter(1))>0) cout << " Resolution ProjectionIhStrip " << (mygaus3->GetParameter(2))/ (mygaus3->GetParameter(1)) << endl;

    TF1* mygaus4 = new TF1("mygaus4","[0]*TMath::Gaus(x,[1],[2])", 1, 6);
    float maxval4=ProjectionIh15Strip->GetXaxis()->GetBinCenter(ProjectionIh15Strip->GetMaximumBin());
    mygaus4->SetParameters(1,maxval4,0.3); 
    ProjectionIh15Strip->Fit("mygaus4"); 
    if ((mygaus4->GetParameter(1))>0) cout << " Resolution ProjectionIh15Strip " << (mygaus4->GetParameter(2))/ (mygaus4->GetParameter(1)) << endl;


    c1->cd();
    ProjectionIhNoL1->SetLineColor(2);
    ProjectionIhNoL1->SetMarkerColor(2);
    mygaus1->SetLineColor(2);
    ProjectionIh15drop->SetLineColor(4);
    ProjectionIh15drop->SetMarkerColor(4);
    mygaus2->SetLineColor(4);
    ProjectionIhStrip->SetLineColor(8);
    ProjectionIhStrip->SetMarkerColor(8);
    mygaus3->SetLineColor(8);
    ProjectionIh15Strip->SetLineColor(7);
    ProjectionIh15Strip->SetMarkerColor(7);
    mygaus4->SetLineColor(7);

    if (ProjectionIhNoL1->GetMaximum() < ProjectionIh15drop->GetMaximum()) ProjectionIhNoL1->SetMaximum( ProjectionIh15drop->GetMaximum()*1.2);
    if (ProjectionIhNoL1->GetMaximum() < ProjectionIhStrip->GetMaximum()) ProjectionIhNoL1->SetMaximum( ProjectionIhStrip->GetMaximum()*1.2);
   
    ProjectionIhNoL1->Draw();
    mygaus1->Draw("same");
    ProjectionIh15drop->Draw("same");
    mygaus2->Draw("same");
/*
    ProjectionIhStrip->Draw("same");
    mygaus3->Draw("same");
*/
    ProjectionIh15Strip->Draw("same");
    mygaus4->Draw("same");



}
void extractQuantile(int caseQ=0){

//TFile *_file0 = TFile::Open("analysis_orig_v2.root");
TFile *_file0 = TFile::Open("analysis_tamas_17v2_18v1.root");

TCanvas *c1 = new TCanvas("Ih", "Ih",10,32,782,552);
    c1->SetFillColor(10);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

TH1D* LdEdXVsRun;
TH1D* LdEdXpixVsRun;
TH1D* LdEdXstripVsRun;

 
 LdEdXVsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdXVsRun"),"LdEdXVsRun",caseQ);
 LdEdXpixVsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdXpixVsRun"),"LdEdXpixVsRun",caseQ);
 LdEdXstripVsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdXstripVsRun"),"LdEdXstripVsRun",caseQ);

 LdEdXpixVsRun->SetMarkerColor(2);
 LdEdXpixVsRun->SetLineColor(2);
 LdEdXstripVsRun->SetMarkerColor(8);
 LdEdXstripVsRun->SetLineColor(8);

  LdEdXVsRun->SetMinimum(2.);
  LdEdXVsRun->SetMaximum(4.5);
//  LdEdXVsRun->GetXaxis()->SetRangeUser(295000, 325000);

 LdEdXVsRun->Draw();
 LdEdXpixVsRun->Draw("same");
 LdEdXstripVsRun->Draw("same");

 if (caseQ==0) c1->SaveAs("Ih_quantile90.pdf");
 else if (caseQ==1) c1->SaveAs("Ih_quantile95.pdf");
 else if (caseQ==2) c1->SaveAs("Ih_quantile99.pdf");
 else if (caseQ==2) c1->SaveAs("Ih_quantile75.pdf");


//
//---



TCanvas *c0 = new TCanvas("Ih_0", "Ih_0",10,32,782,552);
c0->SetFillColor(10);
c0->SetBorderMode(0);
c0->SetBorderSize(2);
c0->SetFrameFillColor(0);
c0->SetFrameBorderMode(0);
c0->cd();
gStyle->SetOptTitle(0);
gStyle->SetOptStat(0);

TH1D* LdEdX0VsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdX0VsRun"),"LdEdX0VsRun",caseQ);
TH1D* LdEdX0pixVsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdX0pixVsRun"),"LdEdX0pixVsRun",caseQ);
TH1D* LdEdX0stripVsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdX0stripVsRun"),"LdEdX0stripVsRun",caseQ);

LdEdX0pixVsRun->SetMarkerColor(2);
LdEdX0pixVsRun->SetLineColor(2);
LdEdX0stripVsRun->SetMarkerColor(8);
LdEdX0stripVsRun->SetLineColor(8);

  LdEdX0VsRun->SetMinimum(2.);
  LdEdX0VsRun->SetMaximum(4.5);
//  LdEdX0VsRun->GetXaxis()->SetRangeUser(295000, 325000);

LdEdX0VsRun->Draw();
LdEdX0pixVsRun->Draw("same");
LdEdX0stripVsRun->Draw("same");
 if (caseQ==0) c0->SaveAs("Ih0_quantile90.pdf");
 else if (caseQ==1) c0->SaveAs("Ih0_quantile95.pdf");
 else if (caseQ==2) c0->SaveAs("Ih0_quantile99.pdf");
 else if (caseQ==3) c0->SaveAs("Ih0_quantile75.pdf");

TCanvas *c2_c = new TCanvas("Ih0 NoL1", "Ih0 NoL1",10,32,782,552);
    c2_c->SetFillColor(10);
    c2_c->SetBorderMode(0);
    c2_c->SetBorderSize(2);
    c2_c->SetFrameFillColor(0);
    c2_c->SetFrameBorderMode(0);
    c2_c->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

TH1D* LadEdX0VsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdX0NoL1VsRun"),"LadEdX0VsRun",caseQ);
TH1D* LadEdX0pixVsRun = extractQuantile((TH2D*)gROOT->FindObject("dEdX0NoL1pixVsRun"),"LadEdX0pixVsRun",caseQ);

LadEdX0pixVsRun->SetMarkerColor(2);
LadEdX0pixVsRun->SetLineColor(2);

LadEdX0VsRun->SetMinimum(2.);
LadEdX0VsRun->SetMaximum(4.5);
//  LadEdX0VsRun->GetXaxis()->SetRangeUser(295000, 325000);

LadEdX0VsRun->Draw();
LadEdX0pixVsRun->Draw("same");
LdEdX0stripVsRun->Draw("same");
 if (caseQ==0) c2_c->SaveAs("Ih0NoL1_quantile90.pdf");
 else if (caseQ==1) c2_c->SaveAs("Ih0NoL1_quantile95.pdf");
 else if (caseQ==2) c2_c->SaveAs("Ih0NoL1_quantile99.pdf");
 else if (caseQ==3) c2_c->SaveAs("Ih0NoL1_quantile75.pdf");



TCanvas *c2_a = new TCanvas("SuperPos","SuperPos",10,32,782,552);
    c2_a->SetFillColor(10);
    c2_a->SetBorderMode(0);
    c2_a->SetBorderSize(2);
    c2_a->SetFrameFillColor(0);
    c2_a->SetFrameBorderMode(0);
    c2_a->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

LadEdX0VsRun->SetMarkerColor(2);
LadEdX0VsRun->SetLineColor(2);
LdEdXVsRun->SetMarkerColor(4);
LdEdXVsRun->SetLineColor(4);
LdEdX0stripVsRun->SetMarkerColor(8);
LdEdX0stripVsRun->SetLineColor(8);
LdEdXstripVsRun->SetMarkerColor(7);
LdEdXstripVsRun->SetLineColor(7);

LadEdX0VsRun->SetMinimum(2.);
LadEdX0VsRun->SetMaximum(4.5);
//LadEdX0VsRun->GetXaxis()->SetRangeUser(295000, 325000);

LadEdX0VsRun->Draw(); // Ih with no L1
LdEdXVsRun->Draw("same"); // Ih with 15% drop
LdEdX0stripVsRun->Draw("same"); // Ih, strip only
LdEdXstripVsRun->Draw("same"); // Ih with 15% drop, strip only

 if (caseQ==0) c2_a->SaveAs("IhSup_quantile90.pdf");
 else if (caseQ==1) c2_a->SaveAs("IhSup_quantile95.pdf");
 else if (caseQ==2) c2_a->SaveAs("IhSup_quantile99.pdf");
 else if (caseQ==3) c2_a->SaveAs("IhSup_quantile75.pdf");


}

TH1D* extractQuantile(TH2D* histo2d, TString namehist, int caseQ){

   TH1D* histo_modif = new TH1D(namehist,namehist, histo2d->GetXaxis()->GetNbins(),histo2d->GetXaxis()->GetXmin(),histo2d->GetXaxis()->GetXmax());
   Double_t xq[4] = {0.9,0.95,0.99,0.75};
   Double_t yq[4];
   if (caseQ<0 || caseQ>=4) {
       cout << "wrong caseQ value for the quantile determination --> between 0 and 3 ! " << endl;
       return histo_modif;
   }
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      if (Projection->Integral()>1000) {
       yq[0]=0;
       yq[1]=0;
       yq[2]=0;
       Projection->GetQuantiles(4,yq,xq);
       histo_modif->SetBinContent(x,yq[caseQ]);
       histo_modif->SetBinError(x,yq[caseQ]*0.001);
      }
   }
   return histo_modif;
}


void testProj(){

TFile *_file0 = TFile::Open("analysis_orig_v2.root");
//TFile *_file0 = TFile::Open("analysis_tamas_17v2_18v1.root");

TCanvas *c1 = new TCanvas("Ih", "Ih",10,32,782,552);
    c1->SetFillColor(10);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->cd();
    gStyle->SetOptTitle(0);

 TH2D* histo2d= (TH2D*) gROOT->FindObject("dEdX0stripVsRun");
 TH1D* histo_modif = new TH1D("histo_proj","histo_proj", 200, 3.1, 3.5);
 TH1D* histo_modif2 = new TH1D("run","run", histo2d->GetXaxis()->GetNbins(),histo2d->GetXaxis()->GetXmin(),histo2d->GetXaxis()->GetXmax());

 float minx=10;
 float maxx=0;
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      
      if (Projection->Integral()>1000) {
       TF1* mygaus = new TF1("mygaus","[0]*TMath::Gaus(x,[1],[2])", 1, 6);
       float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
       mygaus->SetParameters(1,maxval,0.3); 
       mygaus->SetParLimits(1, 0., 10.0);
       mygaus->SetParLimits(2, 0., 10.0);
       TFitResultPtr r = Projection->Fit("mygaus","Q 0 RME S"); 
        if(r->IsValid() && (r->Status()==0 || r->Status()==4000) && mygaus->GetParError(1)<0.1) {
           histo_modif->Fill(mygaus->GetParameter(1));
           histo_modif2->SetBinContent(x,mygaus->GetParameter(1));
           histo_modif2->SetBinError(x,mygaus->GetParError(1));
           if (mygaus->GetParameter(1)<minx) minx=mygaus->GetParameter(1);
           if (mygaus->GetParameter(1)>maxx) maxx=mygaus->GetParameter(1);
        }
      }
   }

  histo_modif->Draw();
  cout << "minx = "<< minx << endl;
  cout << "maxx = "<< maxx << endl;
  cout << "diff = "<< maxx-minx << endl;

TCanvas *c2 = new TCanvas("Ih2", "Ih2",10,32,782,552);
    c2->SetFillColor(10);
    c2->SetBorderMode(0);
    c2->SetBorderSize(2);
    c2->SetFrameFillColor(0);
    c2->SetFrameBorderMode(0);
    c2->cd();
    gStyle->SetOptTitle(0);
   histo_modif2->Draw();

}

void testIh(){

TFile *_file0 = TFile::Open("analysis_orig_v2.root");
//TH2D* LdEdXVsRun0 = (TH2D*)gROOT->FindObject("dEdXVsRun");
//TH2D* LdEdXVsRun0 = (TH2D*)gROOT->FindObject("ChargeVsRun_pixl2");
TH2D* LdEdXVsRun0 = (TH2D*)gROOT->FindObject("dEdXpixVsRun");
TFile *_file1 = TFile::Open("analysis_tamas_v2_n3.root");
//TH2D* LdEdXVsRun1 = (TH2D*)gROOT->FindObject("dEdXVsRun");
//TH2D* LdEdXVsRun1 = (TH2D*)gROOT->FindObject("ChargeVsRun_pixl2");
TH2D* LdEdXVsRun1 = (TH2D*)gROOT->FindObject("dEdXpixVsRun");


TCanvas *c1 = new TCanvas("Ih", "Ih",10,32,782,552);
    c1->SetFillColor(10);
    c1->SetBorderMode(0);
    c1->SetBorderSize(2);
    c1->SetFrameFillColor(0);
    c1->SetFrameBorderMode(0);
    c1->cd();
    gStyle->SetOptTitle(0);
    gStyle->SetOptStat(0);

 TH1D* Projection0 = (TH1D*)(LdEdXVsRun0->ProjectionY("proj",511,511))->Clone();   
 TH1D* Projection1 = (TH1D*)(LdEdXVsRun1->ProjectionY("proj",511,511))->Clone();   

 Projection0->Draw("hist");
 Projection1->SetLineColor(2);
 Projection1->Draw("hist,same");
 c1->SaveAs("test_Ihpix.pdf");

 float a0 = Projection0->GetMean();
 float a1 = Projection1->GetMean();
 cout <<" a1/a0 " << a1*1./a0 << endl;


}

TH2D* modif2DHisto(TH2D* histo, TString namehist) {
     TH2D* modif= (TH2D*) histo->Clone();
     modif->SetName(namehist);
     for(int x=1;x<histo->GetXaxis()->GetNbins()+1;x++){
          TH1D* Projection = (TH1D*)(histo->ProjectionY("proj",x,x))->Clone();
          for (int y=1;y<histo->GetYaxis()->GetNbins()+1;y++) {
            if (Projection->Integral()>0) {
             modif->SetBinContent(x,y,histo->GetBinContent(x,y)/Projection->Integral());
             modif->SetBinError(x,y,histo->GetBinError(x,y)/Projection->Integral());
            }
          }
     }
     return modif;
}


TH2D* correl2DHisto(TH1D* histo1, TH1D* histo2, TString namehist, TString namecap){
//     TH2D* modif = new TH2D(namehist,namecap,50,histo1->GetMinimum()*0.90, histo1->GetMaximum()*1.1, 50, histo2->GetMinimum()*0.90, histo2->GetMaximum()*1.1); 
     TH2D* modif = new TH2D(namehist,namecap,50,histo1->GetMinimum(), histo1->GetMaximum(), 50, histo2->GetMinimum(), histo2->GetMaximum()); 
     for (int i=0; i<histo1->GetXaxis()->GetNbins()+1; i++) {
        float valx1=histo1->GetBinContent(i);
        float valx2=histo2->GetBinContent(i);
        if (valx1>0 && valx2>0) {
           modif->Fill(valx1,valx2);
        }
  
     }
     modif->SetMarkerStyle(21);
     return modif;

}

TH1D* extractVariationGaus(TH2D* histo2d, TString namehist, int nbin, float x1, float x2){

   TH1D* histo_modif = new TH1D(namehist,namehist, nbin, x1, x2);
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      
      if (Projection->Integral()>1000) {
       TF1* mygaus = new TF1("mygaus","[0]*TMath::Gaus(x,[1],[2])", 1, 6);
       float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
       mygaus->SetParameters(1,maxval,0.3); 
       mygaus->SetParLimits(1, 0., 10.0);
       mygaus->SetParLimits(2, 0., 10.0);
//       int fstatus = Projection->Fit("mygaus","Q 0 RME");
       TFitResultPtr r = Projection->Fit("mygaus","Q 0 RME S"); 
//       if (!check && mygaus->GetParError(1)<0.3) {
       if (!check) {
        if(r->IsValid() && (r->Status()==0 || r->Status()==4000) && mygaus->GetParError(1)<0.1) {
         histo_modif->Fill(mygaus->GetParameter(1),Projection->Integral());
/*
         if( mygaus->GetParameter(1)>0) {
           histo_modif->SetBinContent(x,mygaus->GetParameter(2)/mygaus->GetParameter(1));
           histo_modif->SetBinError(x,mygaus->GetParError(2)/mygaus->GetParameter(1));
         }
*/
        }

       }
      }
   }
   return histo_modif;
}
TH1D* extractVariationMaximum(TH2D* histo2d,  TString namehist, int nbin, float x1, float x2){

   TH1D* histo_modif = new TH1D(namehist,namehist, nbin, x1, x2);
   for(int x=1;x<=histo2d->GetNbinsX();x++){
      TH1D* Projection = (TH1D*)(histo2d->ProjectionY("proj",x,x))->Clone();   
      
      if (Projection->Integral()>100) {
//       float maxval=Projection->GetXaxis()->GetBinCenter(Projection->GetMaximumBin());
//       float rmsval=Projection->GetXaxis()->GetBinWidth(Projection->GetMaximumBin());
       float maxval=Projection->GetMean();
       float rmsval=Projection->GetMeanError();
       histo_modif->Fill(maxval,Projection->Integral());
      }
   }
   return histo_modif;
}
