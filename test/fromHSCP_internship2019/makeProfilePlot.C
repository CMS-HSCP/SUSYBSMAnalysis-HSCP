{
TFile *_file0 = TFile::Open("root_13mai/mbias_template_corr.root");
TFile *_file1 = TFile::Open("root_13mai/mbias_template_uncorr.root");
TFile *_file2 = TFile::Open("root_13mai/ttbar_template_corr.root");
TFile *_file3 = TFile::Open("root_13mai/ttbar_template_uncorr.root");

TProfile *CorrPixel;
TProfile *CorrStrip;
TProfile *UnCorrPixel;
TProfile *UnCorrStrip;

TProfile *CorrPixelExt;
TProfile *CorrStripExt;
TProfile *UnCorrPixelExt;
TProfile *UnCorrStripExt;

TProfile *CorrPixelttbar;
TProfile *CorrStripttbar;
TProfile *UnCorrPixelttbar;
TProfile *UnCorrStripttbar;

TProfile *CorrPixelExttt;
TProfile *CorrStripExttt;
TProfile *UnCorrPixelExttt;
TProfile *UnCorrStripExttt;

TProfile *CorrPixelw;
TProfile *CorrStripw;
TProfile *UnCorrPixelw;
TProfile *UnCorrStripw;

TProfile *CorrPixelExtw;
TProfile *CorrStripExtw;
TProfile *UnCorrPixelExtw;
TProfile *UnCorrStripExtw;

TProfile *CorrPixelttbarw;
TProfile *CorrStripttbarw;
TProfile *UnCorrPixelttbarw;
TProfile *UnCorrStripttbarw;
TProfile *CorrPixelExtttw;
TProfile *CorrStripExtttw;
TProfile *UnCorrPixelExtttw;
TProfile *UnCorrStripExtttw;

TH2D* CorrHHitPixVsEtaP5;
TProfile *CorrPixVsEtaP5;
TH2D* CorrHHitPixVsEtapL5;
TProfile *CorrPixVsEtapL5;

TH2D* UnCorrHHitPixVsEtaP5;
TProfile *UnCorrPixVsEtaP5;
TH2D* UnCorrHHitPixVsEtapL5;
TProfile *UnCorrPixVsEtapL5;

TH2D* CorrHHitStripVsEtaP5;
TProfile *CorrStripVsEtaP5;
TH2D* CorrHHitStripVsEtapL5;
TProfile *CorrStripVsEtapL5;

TH2D* UnCorrHHitStripVsEtaP5;
TProfile *UnCorrStripVsEtaP5;
TH2D* UnCorrHHitStripVsEtapL5;
TProfile *UnCorrStripVsEtapL5;

TProfile *cCorrPixVsEtaP5;
TProfile *cCorrPixVsEtapL5;
TProfile *cUnCorrPixVsEtaP5;
TProfile *cUnCorrPixVsEtapL5;
TProfile *cCorrStripVsEtaP5;
TProfile *cCorrStripVsEtapL5;
TProfile *cUnCorrStripVsEtaP5;
TProfile *cUnCorrStripVsEtapL5;
TProfile *tCorrPixVsEtaP5;
TProfile *tCorrPixVsEtapL5;
TProfile *tUnCorrPixVsEtaP5;
TProfile *tUnCorrPixVsEtapL5;
TProfile *tCorrStripVsEtaP5;
TProfile *tCorrStripVsEtapL5;
TProfile *tUnCorrStripVsEtaP5;
TProfile *tUnCorrStripVsEtapL5;

_file0->cd();


CorrPixel= HHitPixVsP->ProfileX();
CorrPixelExt= HHitPixVsPExt->ProfileX();

CorrStrip= HHitStripVsP->ProfileX();
CorrStripExt= HHitStripVsPExt->ProfileX();

CorrPixelw= HHitPixVsPw->ProfileX();
CorrPixelExtw= HHitPixVsPExtw->ProfileX();

CorrStripw= HHitStripVsPw->ProfileX();
CorrStripExtw= HHitStripVsPExtw->ProfileX();

CorrHHitStripVsEtaP5 = (TH2D*)gROOT->FindObject("HHitStripVsEtaP5");
CorrStripVsEtaP5 = HHitStripVsEtaP5->ProfileX(); 
CorrHHitStripVsEtapL5 = (TH2D*) gROOT->FindObject("HHitStripVsEtapL5");
CorrStripVsEtapL5 = HHitStripVsEtapL5->ProfileX(); 
CorrHHitPixVsEtaP5 = (TH2D*) gROOT->FindObject("HHitPixVsEtaP5");
CorrPixVsEtaP5 = HHitPixVsEtaP5->ProfileX(); 
CorrHHitPixVsEtapL5 = (TH2D*) gROOT->FindObject("HHitPixVsEtapL5");
CorrPixVsEtapL5 = HHitPixVsEtapL5->ProfileX(); 

cCorrPixVsEtaP5 = HHitPixVsEtaP5->ProfileX(); 
cCorrPixVsEtapL5 = HHitPixVsEtapL5->ProfileX(); 
cCorrStripVsEtaP5 = HHitStripVsEtaP5->ProfileX(); 
cCorrStripVsEtapL5 = HHitStripVsEtapL5->ProfileX(); 
/*
cCorrPixVsEtaP5 = (TProfile*) CorrPixVsEtaP5->Clone();
cCorrPixVsEtapL5 = (TProfile*) CorrPixVsEtapL5->Clone();
cCorrStripVsEtaP5 = (TProfile*) CorrStripVsEtaP5->Clone();
CorrStripVsEtapL5 = (TProfile*) CorrStripVsEtapL5->Clone();
*/

cout <<  "fill from file0"<< endl;
_file1->cd();

UnCorrPixel= HHitPixVsP->ProfileX();
UnCorrPixelExt= HHitPixVsPExt->ProfileX();

UnCorrStrip= HHitStripVsP->ProfileX();
UnCorrStripExt= HHitStripVsPExt->ProfileX();

UnCorrPixelw= HHitPixVsPw->ProfileX();
UnCorrPixelExtw= HHitPixVsPExtw->ProfileX();

UnCorrStripw= HHitStripVsPw->ProfileX();
UnCorrStripExtw= HHitStripVsPExtw->ProfileX();

UnCorrHHitStripVsEtaP5 = (TH2D*)gROOT->FindObject("HHitStripVsEtaP5");
UnCorrStripVsEtaP5 = HHitStripVsEtaP5->ProfileX(); 
UnCorrHHitStripVsEtapL5 = (TH2D*) gROOT->FindObject("HHitStripVsEtapL5");
UnCorrStripVsEtapL5 = HHitStripVsEtapL5->ProfileX(); 
UnCorrHHitPixVsEtaP5 = (TH2D*) gROOT->FindObject("HHitPixVsEtaP5");
UnCorrPixVsEtaP5 = HHitPixVsEtaP5->ProfileX(); 
UnCorrHHitPixVsEtapL5 = (TH2D*) gROOT->FindObject("HHitPixVsEtapL5");
UnCorrPixVsEtapL5 = HHitPixVsEtapL5->ProfileX(); 

cUnCorrPixVsEtaP5 = HHitPixVsEtaP5->ProfileX(); 
cUnCorrPixVsEtapL5 = HHitPixVsEtapL5->ProfileX(); 
cUnCorrStripVsEtaP5 = HHitStripVsEtaP5->ProfileX(); 
cUnCorrStripVsEtapL5 = HHitStripVsEtapL5->ProfileX(); 

cout <<  "fill from file1"<< endl;
_file2->cd();
CorrPixelttbar= HHitPixVsP->ProfileX();
CorrPixelExttt= HHitPixVsPExt->ProfileX();
CorrStripttbar= HHitStripVsP->ProfileX();
CorrStripExttt= HHitStripVsPExt->ProfileX();

CorrPixelttbarw= HHitPixVsPw->ProfileX();
CorrStripttbarw= HHitStripVsPw->ProfileX();
CorrPixelExtttw= HHitPixVsPExtw->ProfileX();
CorrStripExtttw= HHitStripVsPExtw->ProfileX();

tCorrPixVsEtaP5 = HHitPixVsEtaP5->ProfileX(); 
tCorrPixVsEtapL5 = HHitPixVsEtapL5->ProfileX(); 
tCorrStripVsEtaP5 = HHitStripVsEtaP5->ProfileX(); 
tCorrStripVsEtapL5 = HHitStripVsEtapL5->ProfileX(); 

cout <<  "fill from file2"<< endl;
_file3->cd();
UnCorrPixelttbar= HHitPixVsP->ProfileX();
UnCorrPixelExttt= HHitPixVsPExt->ProfileX();
UnCorrStripttbar= HHitStripVsP->ProfileX();
UnCorrStripExttt= HHitStripVsPExt->ProfileX();
UnCorrPixelttbarw= HHitPixVsPw->ProfileX();
UnCorrPixelExtttw= HHitPixVsPExtw->ProfileX();
UnCorrStripttbarw= HHitStripVsPw->ProfileX();
UnCorrStripExtttw= HHitStripVsPExtw->ProfileX();
tUnCorrPixVsEtaP5 = HHitPixVsEtaP5->ProfileX(); 
tUnCorrPixVsEtapL5 = HHitPixVsEtapL5->ProfileX(); 
tUnCorrStripVsEtaP5 = HHitStripVsEtaP5->ProfileX(); 
tUnCorrStripVsEtapL5 = HHitStripVsEtapL5->ProfileX(); 

cout <<  "fill from file3"<< endl;

//000

CorrPixel->SetLineColor(4);
CorrPixel->SetMarkerColor(4);
CorrPixel->SetMarkerStyle(21);
UnCorrPixel->SetLineColor(4);
UnCorrPixel->SetMarkerColor(4);
UnCorrPixel->SetMarkerStyle(21);

CorrPixelExt->SetLineColor(4);
CorrPixelExt->SetMarkerColor(4);
CorrPixelExt->SetMarkerStyle(26);
UnCorrPixelExt->SetLineColor(4);
UnCorrPixelExt->SetMarkerColor(4);
UnCorrPixelExt->SetMarkerStyle(26);


CorrPixelttbar->SetLineColor(6);
CorrPixelttbar->SetMarkerColor(6);
CorrPixelttbar->SetMarkerStyle(20);
UnCorrPixelttbar->SetLineColor(6);
UnCorrPixelttbar->SetMarkerColor(6);
UnCorrPixelttbar->SetMarkerStyle(20);

CorrPixelExttt->SetLineColor(6);
CorrPixelExttt->SetMarkerColor(6);
CorrPixelExttt->SetMarkerStyle(25);
UnCorrPixelExttt->SetLineColor(6);
UnCorrPixelExttt->SetMarkerColor(6);
UnCorrPixelExttt->SetMarkerStyle(25);

//000

CorrStrip->SetLineColor(2);
CorrStrip->SetMarkerColor(2);
CorrStrip->SetMarkerStyle(21);
UnCorrStrip->SetLineColor(2);
UnCorrStrip->SetMarkerColor(2);
UnCorrStrip->SetMarkerStyle(21);

CorrStripExt->SetLineColor(2);
CorrStripExt->SetMarkerColor(2);
CorrStripExt->SetMarkerStyle(26);
UnCorrStripExt->SetLineColor(2);
UnCorrStripExt->SetMarkerColor(2);
UnCorrStripExt->SetMarkerStyle(26);


CorrStripttbar->SetLineColor(8);
CorrStripttbar->SetMarkerColor(8);
CorrStripttbar->SetMarkerStyle(20);
UnCorrStripttbar->SetLineColor(8);
UnCorrStripttbar->SetMarkerColor(8);
UnCorrStripttbar->SetMarkerStyle(20);

CorrStripExttt->SetLineColor(8);
CorrStripExttt->SetMarkerColor(8);
CorrStripExttt->SetMarkerStyle(25);
UnCorrStripExttt->SetLineColor(8);
UnCorrStripExttt->SetMarkerColor(8);
UnCorrStripExttt->SetMarkerStyle(25);

//000
//
CorrPixelw->SetLineColor(4);
CorrPixelw->SetMarkerColor(4);
CorrPixelw->SetMarkerStyle(21);
UnCorrPixelw->SetLineColor(4);
UnCorrPixelw->SetMarkerColor(4);
UnCorrPixelw->SetMarkerStyle(21);

CorrPixelExtw->SetLineColor(4);
CorrPixelExtw->SetMarkerColor(4);
CorrPixelExtw->SetMarkerStyle(26);
UnCorrPixelExtw->SetLineColor(4);
UnCorrPixelExtw->SetMarkerColor(4);
UnCorrPixelExtw->SetMarkerStyle(26);


CorrPixelttbarw->SetLineColor(6);
CorrPixelttbarw->SetMarkerColor(6);
CorrPixelttbarw->SetMarkerStyle(20);
UnCorrPixelttbarw->SetLineColor(6);
UnCorrPixelttbarw->SetMarkerColor(6);
UnCorrPixelttbarw->SetMarkerStyle(20);

CorrPixelExtttw->SetLineColor(6);
CorrPixelExtttw->SetMarkerColor(6);
CorrPixelExtttw->SetMarkerStyle(25);
UnCorrPixelExtttw->SetLineColor(6);
UnCorrPixelExtttw->SetMarkerColor(6);
UnCorrPixelExtttw->SetMarkerStyle(25);

//0000

CorrStripw->SetLineColor(2);
CorrStripw->SetMarkerColor(2);
CorrStripw->SetMarkerStyle(21);
UnCorrStripw->SetLineColor(2);
UnCorrStripw->SetMarkerColor(2);
UnCorrStripw->SetMarkerStyle(21);

CorrStripExtw->SetLineColor(2);
CorrStripExtw->SetMarkerColor(2);
CorrStripExtw->SetMarkerStyle(26);
UnCorrStripExtw->SetLineColor(2);
UnCorrStripExtw->SetMarkerColor(2);
UnCorrStripExtw->SetMarkerStyle(26);


CorrStripttbarw->SetLineColor(8);
CorrStripttbarw->SetMarkerColor(8);
CorrStripttbarw->SetMarkerStyle(20);
UnCorrStripttbarw->SetLineColor(8);
UnCorrStripttbarw->SetMarkerColor(8);
UnCorrStripttbarw->SetMarkerStyle(20);

CorrStripExtttw->SetLineColor(8);
CorrStripExtttw->SetMarkerColor(8);
CorrStripExtttw->SetMarkerStyle(25);
UnCorrStripExtttw->SetLineColor(8);
UnCorrStripExtttw->SetMarkerColor(8);
UnCorrStripExtttw->SetMarkerStyle(25);

TCanvas *c1 = new TCanvas("c1", "CorrPixel", 600,600);
c1->cd();

CorrPixel->SetStats(kFALSE);
CorrPixel->SetMinimum(0);
CorrPixel->SetMaximum(10);
CorrPixel->Draw();
CorrPixelExt->Draw("same");
//CorrPixelttbar->Draw("same");
CorrStrip->Draw("same");
CorrStripExt->Draw("same");
//CorrStripttbar->Draw("same");


TLegend* qw =  new TLegend(0.15,0.12,0.40,0.28);
qw->AddEntry(CorrPixel,     "Pixel |#eta|<0.4 (MB)",                       "p");
qw->AddEntry(CorrPixelExt,     "Pixel |#eta|>0.4 (MB)",                       "p");
qw->AddEntry(CorrPixelttbar,     "Pixel |#eta|<0.4 (tt)",                       "p");
qw->AddEntry(CorrPixelExttt,     "Pixel |#eta|>0.4 (tt)",                       "p");

TLegend* qw2 =  new TLegend(0.4,0.12,0.65,0.28);
qw2->AddEntry(CorrStrip,     "Strip |#eta|<0.4 (MB)",                       "p");
qw2->AddEntry(CorrStripExt,     "Strip |#eta|>0.4 (MB)",                       "p");
qw2->AddEntry(CorrStripttbar,     "Strip |#eta|<0.4 (tt)",                       "p");
qw2->AddEntry(CorrStripExttt,     "Strip |#eta|>0.4 (tt)",                       "p");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c1 ok"<< endl;

TCanvas *c1b = new TCanvas("c1b", "CorrPixelb", 600,600);
c1b->cd();

CorrPixel->Draw();
CorrPixelttbar->Draw("same");
CorrStrip->Draw("same");
CorrStripttbar->Draw("same");


qw->Draw("same");
qw2->Draw("same");

cout <<  "c1b ok"<< endl;

TCanvas *c1c = new TCanvas("c1c", "CorrPixelc", 600,600);
c1c->cd();

CorrPixelttbar->SetStats(kFALSE);
CorrPixelttbar->SetMinimum(0);
CorrPixelttbar->SetMaximum(10);
CorrPixelttbar->Draw();
CorrPixelExttt->Draw("same");
CorrStripttbar->Draw("same");
CorrStripExttt->Draw("same");


qw->Draw("same");
qw2->Draw("same");

cout <<  "c1c ok"<< endl;




TCanvas *c2 = new TCanvas("c2", "UnCorrPixel", 600,600);
c2->cd();                       
                                
UnCorrPixel->SetStats(kFALSE);
UnCorrPixel->SetMinimum(0);
UnCorrPixel->SetMaximum(10);
UnCorrPixel->Draw();              
UnCorrPixelExt->Draw("same");
UnCorrStrip->Draw("same");              
UnCorrStripExt->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c2 ok"<< endl;

TCanvas *c2b = new TCanvas("c2b", "UnCorrPixelb", 600,600);
c2b->cd();

UnCorrPixel->SetStats(kFALSE);
UnCorrPixel->SetMinimum(0);
UnCorrPixel->SetMaximum(10);
UnCorrPixel->Draw();
UnCorrPixelttbar->Draw("same");
UnCorrStrip->Draw("same");
UnCorrStripttbar->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c2b ok"<< endl;

TCanvas *c2c = new TCanvas("c2c", "UnCorrPixelc", 600,600);
c2c->cd();

UnCorrPixelttbar->SetStats(kFALSE);
UnCorrPixelttbar->SetMinimum(0);
UnCorrPixelttbar->SetMaximum(10);
UnCorrPixelttbar->Draw();
UnCorrPixelExttt->Draw("same");
UnCorrStripttbar->Draw("same");
UnCorrStripExttt->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c2c ok"<< endl;


TCanvas *c3 = new TCanvas("c3", "CorrPixelw", 600,600);
c3->cd();

CorrPixelw->SetStats(kFALSE);
CorrPixelw->SetMinimum(0);
CorrPixelw->SetMaximum(10);
CorrPixelw->Draw();
CorrPixelExtw->Draw("same");
CorrStripw->Draw("same");
CorrStripExtw->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c3 ok"<< endl;
TCanvas *c3b = new TCanvas("c3b", "CorrPixelw", 600,600);
c3b->cd();

CorrPixelw->Draw();
CorrPixelttbarw->Draw("same");
CorrStripw->Draw("same");
CorrStripttbarw->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c3b ok"<< endl;
TCanvas *c3c = new TCanvas("c3c", "CorrPixelw", 600,600);
c3c->cd();

CorrPixelttbarw->SetStats(kFALSE);
CorrPixelttbarw->SetMinimum(0);
CorrPixelttbarw->SetMaximum(10);
CorrPixelttbarw->Draw();
CorrPixelttbarw->Draw();
CorrPixelExtttw->Draw("same");
CorrStripttbarw->Draw("same");
CorrStripExtttw->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c3c ok"<< endl;
TCanvas *c4 = new TCanvas("c4", "UnCorrPixelw", 600,600);
c4->cd();                       
                                
UnCorrPixelw->SetStats(kFALSE);
UnCorrPixelw->SetMinimum(0);
UnCorrPixelw->SetMaximum(10);
UnCorrPixelw->Draw();              
UnCorrPixelExtw->Draw("same");
UnCorrStripw->Draw("same");              
UnCorrStripExtw->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c4 ok"<< endl;

TCanvas *c4b = new TCanvas("c4b", "UnCorrPixelw", 600,600);
c4b->cd();                       
                                
UnCorrPixelw->Draw();              
UnCorrPixelttbarw->Draw("same");
UnCorrStripw->Draw("same");              
UnCorrStripttbarw->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c4b ok"<< endl;
TCanvas *c4c = new TCanvas("c4c", "UnCorrPixelw", 600,600);
c4c->cd();                       
                                
UnCorrPixelttbarw->SetStats(kFALSE);
UnCorrPixelttbarw->SetMinimum(0);
UnCorrPixelttbarw->SetMaximum(10);
UnCorrPixelttbarw->Draw();              
UnCorrPixelExtttw->Draw("same");
UnCorrStripttbarw->Draw("same");
UnCorrStripExtttw->Draw("same");
qw->Draw("same");
qw2->Draw("same");

cout <<  "c4c ok"<< endl;


TCanvas *c5 = new TCanvas("c5", "EtaDep1", 800,600);
c5->Divide(2,2);                       
c5->cd(2);                       
CorrHHitPixVsEtaP5->SetStats(kFALSE);
CorrHHitPixVsEtaP5 ->Draw("colz");
CorrPixVsEtaP5->Draw("same");

c5->cd(1);                       
CorrHHitPixVsEtapL5->SetStats(kFALSE);
CorrHHitPixVsEtapL5 ->Draw("colz");
CorrPixVsEtapL5->Draw("same");

c5->cd(4);                       
CorrHHitStripVsEtaP5->SetStats(kFALSE);
CorrHHitStripVsEtaP5 ->Draw("colz");
CorrStripVsEtaP5->Draw("same");

c5->cd(3);                       
CorrHHitStripVsEtapL5->SetStats(kFALSE);
CorrHHitStripVsEtapL5 ->Draw("colz");
CorrStripVsEtapL5->Draw("same");

c5->cd();
c5->SaveAs("colorProfile.png");
cout <<  "c5 ok"<< endl;


TCanvas *c6 = new TCanvas("c6", "Corr EtaDep", 600,600);
c6->cd();

cCorrPixVsEtaP5->SetLineColor(2);
cCorrPixVsEtaP5->SetMarkerColor(2);
cCorrPixVsEtaP5->SetMarkerStyle(21);
cCorrPixVsEtapL5->SetLineColor(2);
cCorrPixVsEtapL5->SetMarkerColor(2);
cCorrPixVsEtapL5->SetMarkerStyle(26);

cUnCorrPixVsEtaP5->SetLineColor(2);
cUnCorrPixVsEtaP5->SetMarkerColor(2);
cUnCorrPixVsEtaP5->SetMarkerStyle(21);
cUnCorrPixVsEtapL5->SetLineColor(2);
cUnCorrPixVsEtapL5->SetMarkerColor(2);
cUnCorrPixVsEtapL5->SetMarkerStyle(26);

cCorrStripVsEtaP5->SetLineColor(4);
cCorrStripVsEtaP5->SetMarkerColor(4);
cCorrStripVsEtaP5->SetMarkerStyle(21);
cCorrStripVsEtapL5->SetLineColor(4);
cCorrStripVsEtapL5->SetMarkerColor(4);
cCorrStripVsEtapL5->SetMarkerStyle(26);

cUnCorrStripVsEtaP5->SetLineColor(4);
cUnCorrStripVsEtaP5->SetMarkerColor(4);
cUnCorrStripVsEtaP5->SetMarkerStyle(21);
cUnCorrStripVsEtapL5->SetLineColor(4);
cUnCorrStripVsEtapL5->SetMarkerColor(4);
cUnCorrStripVsEtapL5->SetMarkerStyle(26);

cCorrPixVsEtaP5->SetStats(kFALSE);
cCorrPixVsEtaP5->SetMaximum(10);
cCorrPixVsEtaP5->SetMinimum(1);
cCorrPixVsEtaP5->Draw();
cCorrPixVsEtapL5->Draw("same");
cCorrStripVsEtaP5->Draw("same");
cCorrStripVsEtapL5->Draw("same");

TLegend* qw3 =  new TLegend(0.15,0.12,0.40,0.28);
qw3->AddEntry(cCorrPixVsEtaP5,    "Pixel p>5 (MB)",                       "p");
qw3->AddEntry(cCorrPixVsEtapL5,   "Pixel p<5 (MB)",                       "p");
qw3->AddEntry(cCorrStripVsEtaP5,      "Strip p>5 (MB)",                       "p");
qw3->AddEntry(cCorrStripVsEtapL5,     "Strip p<5 (MB)",                       "p");
qw3->Draw("same");

TCanvas *c6b = new TCanvas("c6b", "Corr EtaDep", 600,600);
c6b->cd();

tCorrPixVsEtaP5->SetLineColor(6);
tCorrPixVsEtaP5->SetMarkerColor(6);
tCorrPixVsEtaP5->SetMarkerStyle(20);
tCorrPixVsEtapL5->SetLineColor(6);
tCorrPixVsEtapL5->SetMarkerColor(6);
tCorrPixVsEtapL5->SetMarkerStyle(25);

tUnCorrPixVsEtaP5->SetLineColor(6);
tUnCorrPixVsEtaP5->SetMarkerColor(6);
tUnCorrPixVsEtaP5->SetMarkerStyle(20);
tUnCorrPixVsEtapL5->SetLineColor(6);
tUnCorrPixVsEtapL5->SetMarkerColor(6);
tUnCorrPixVsEtapL5->SetMarkerStyle(25);

tCorrStripVsEtaP5->SetLineColor(8);
tCorrStripVsEtaP5->SetMarkerColor(8);
tCorrStripVsEtaP5->SetMarkerStyle(20);
tCorrStripVsEtapL5->SetLineColor(8);
tCorrStripVsEtapL5->SetMarkerColor(8);
tCorrStripVsEtapL5->SetMarkerStyle(25);

tUnCorrStripVsEtaP5->SetLineColor(8);
tUnCorrStripVsEtaP5->SetMarkerColor(8);
tUnCorrStripVsEtaP5->SetMarkerStyle(20);
tUnCorrStripVsEtapL5->SetLineColor(8);
tUnCorrStripVsEtapL5->SetMarkerColor(8);
tUnCorrStripVsEtapL5->SetMarkerStyle(25);

tCorrPixVsEtaP5->SetStats(kFALSE);
tCorrPixVsEtaP5->SetMaximum(10);
tCorrPixVsEtaP5->SetMinimum(1);
tCorrPixVsEtaP5->Draw();
tCorrPixVsEtapL5->Draw("same");
tCorrStripVsEtaP5->Draw("same");
tCorrStripVsEtapL5->Draw("same");

TLegend* qw4 =  new TLegend(0.15,0.12,0.40,0.28);
qw4->AddEntry(tCorrPixVsEtaP5,    "Pixel p>5 (tt)",                       "p");
qw4->AddEntry(tCorrPixVsEtapL5,   "Pixel p<5 (tt)",                       "p");
qw4->AddEntry(tCorrStripVsEtaP5,      "Strip p>5 (tt)",                       "p");
qw4->AddEntry(tCorrStripVsEtapL5,     "Strip p<5 (tt)",                       "p");
qw4->Draw("same");

TCanvas *c7 = new TCanvas("c7", "UnCorr EtaDep", 600,600);
c7->cd();

cUnCorrPixVsEtaP5->SetStats(kFALSE);
cUnCorrPixVsEtaP5->SetMaximum(10);
cUnCorrPixVsEtaP5->SetMinimum(1);
cUnCorrPixVsEtaP5->Draw();
cUnCorrPixVsEtapL5->Draw("same");
cUnCorrStripVsEtaP5->Draw("same");
cUnCorrStripVsEtapL5->Draw("same");
qw4->Draw("same");

TCanvas *c7b = new TCanvas("c7", "UnCorr EtaDep", 600,600);
c7b->cd();

tUnCorrPixVsEtaP5->SetStats(kFALSE);
tUnCorrPixVsEtaP5->SetMaximum(10);
tUnCorrPixVsEtaP5->SetMinimum(1);
tUnCorrPixVsEtaP5->Draw();
tUnCorrPixVsEtapL5->Draw("same");
tUnCorrStripVsEtaP5->Draw("same");
tUnCorrStripVsEtapL5->Draw("same");
qw4->Draw("same");

}

