#include "Analysis_Global.h"
//#include "Analysis_CommonFunction.h"
#include "Analysis_PlotFunction.h"
//#include "Analysis_PlotStructure.h"
//#include "Analysis_Samples.h"
#include "tdrstyle.C"
#include "TMarker.h"

void MassLimitSummaryPlot(){
  setTDRStyle();
  gStyle->SetPadTopMargin   (0.06);
  gStyle->SetPadBottomMargin(0.16);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadLeftMargin  (0.14);
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.7);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505);
  //gStyle->SetTextFont(43);

   TCanvas *c1 = new TCanvas("c1","c1",600,600);
   Double_t mchx[12] = {0.0, 0.3333, 0.6666, 1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
   Double_t mchy[12] = {-100.0, 200.0, 480.0, 574.0, 685.0, 752.0, 793.0, 796.0, 781.0, 757.0, 715.0};
   TGraph* gr1 = new TGraph(11,mchx,mchy);
   gr1->SetName("gr1");
   gr1->SetTitle("");
   gr1->SetMaximum(1700.0);
   gr1->SetMinimum(0.0);
   gr1->SetMarkerStyle(20);
   gr1->SetMarkerSize(1.3);
   gr1->SetMarkerColor(2);
   gr1->GetXaxis()->SetTitle("Charge (e)");
   gr1->GetYaxis()->SetTitle("95% C.L. lower mass limit (GeV/#font[12]{c}^{2})");
   gr1->GetYaxis()->SetTitleOffset(1.9);
   gr1->Draw("AP");

   //fitting highly charged points
   Double_t highQmchx[10] = {1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0};
   Double_t highQmchy[10] = {517, 685.0, 752.0, 793.0, 796.0, 781.0, 757.0, 715.0};
   TGraph* highQgr1 = new TGraph(8,highQmchx,highQmchy);
   TF1* highQFit = new TF1("highQFit","pol2", 1, 8);
   highQgr1->Fit("highQFit","QR0","");
   printf("FitParameter for Q>=1= %f + %fQ + %fQ^2\n", highQFit->GetParameter(0), highQFit->GetParameter(1), highQFit->GetParameter(2));
   highQFit->SetLineColor(2);
   highQFit->SetLineStyle(2);
   highQFit->SetLineWidth(2);
   //highQFit->Draw("same"); //was in the 2012 AN
   //add tje mCHAMP Q=1 point with a smaller size
   TMarker* mCHAMPQ1 = new TMarker(1.0, 517.0, 20);
   mCHAMPQ1->SetMarkerSize(0.8);
   mCHAMPQ1->SetMarkerColor(2);
   //mCHAMPQ1->Draw(); //was in the 2012 AN

   //previous CMS limits on frac. charge
   TMarker* PCDY1e3 = new TMarker(0.3333, 145, 21);
   PCDY1e3->SetMarkerSize(1.3);
   PCDY1e3->SetMarkerColor(4);
   PCDY1e3->Draw();
   TMarker* PCDY2e3 = new TMarker(0.6666, 310, 21);
   PCDY2e3->SetMarkerSize(1.3);
   PCDY2e3->SetMarkerColor(4);
   PCDY2e3->Draw();

   //previous CMS limits on frac. charge
   TMarker* PCDY1e3New = new TMarker(0.3333, 172, 24);
   PCDY1e3New->SetMarkerSize(1.3);
   PCDY1e3New->SetMarkerColor(2);
   PCDY1e3New->Draw();
   TMarker* PCDY2e3New = new TMarker(0.6666, 378, 24);
   PCDY2e3New->SetMarkerSize(1.3);
   PCDY2e3New->SetMarkerColor(2);
   PCDY2e3New->Draw();


   //Previous ATLAS result on multiply charged
   TMarker* PADY2e = new TMarker(2, 430, 23);
   PADY2e->SetMarkerSize(1.3);
   PADY2e->SetMarkerColor(1);
   PADY2e->Draw();
   TMarker* PADY3e = new TMarker(3, 480, 23);
   PADY3e->SetMarkerSize(1.3);
   PADY3e->SetMarkerColor(1);
   PADY3e->Draw();
   TMarker* PADY4e = new TMarker(4, 490, 23);
   PADY4e->SetMarkerSize(1.3);
   PADY4e->SetMarkerColor(1);
   PADY4e->Draw();
   TMarker* PADY5e = new TMarker(5, 470, 23);
   PADY5e->SetMarkerSize(1.3);
   PADY5e->SetMarkerColor(1);
   PADY5e->Draw();
   TMarker* PADY6e = new TMarker(6, 420, 23);
   PADY6e->SetMarkerSize(1.3);
   PADY6e->SetMarkerColor(1);
   PADY6e->Draw();

   //Previous CMS result on multiply charged at 7TeV
   TMarker* PCDY1e = new TMarker(1, 427, 34);
   PCDY1e->SetMarkerSize(1.3);
   PCDY1e->SetMarkerColor(8);
   //PCDY1e->Draw();  //was in the 2012AN
   TMarker* PCDY2e = new TMarker(2, 489, 34);
   PCDY2e->SetMarkerSize(1.3);
   PCDY2e->SetMarkerColor(8);
   //PCDY2e->Draw(); //was in the 2012AN
   TMarker* PCDY3e = new TMarker(3, 545, 34);
   PCDY3e->SetMarkerSize(1.3);
   PCDY3e->SetMarkerColor(8);
   //PCDY3e->Draw(); was in the 2012AN
   TMarker* PCDY4e = new TMarker(4, 572, 34);
   PCDY4e->SetMarkerSize(1.3);
   PCDY4e->SetMarkerColor(8);
   //PCDY4e->Draw(); was in the 2012AN
   TMarker* PCDY5e = new TMarker(5, 504, 34);
   PCDY5e->SetMarkerSize(1.3);
   PCDY5e->SetMarkerColor(8);
   //PCDY5e->Draw(); //was in the 2012AN


   //CMS result on multiply charged at 7TeV
   TMarker* PCNDY1e = new TMarker(1, 411, 24);
   PCNDY1e->SetMarkerSize(1.3);
   PCNDY1e->SetMarkerColor(2);
   PCNDY1e->Draw();
   TMarker* PCNDY2e = new TMarker(2, 543, 24);
   PCNDY2e->SetMarkerSize(1.3);
   PCNDY2e->SetMarkerColor(2);
   PCNDY2e->Draw();
   TMarker* PCNDY3e = new TMarker(3, 591, 24);
   PCNDY3e->SetMarkerSize(1.3);
   PCNDY3e->SetMarkerColor(2);
   PCNDY3e->Draw();
   TMarker* PCNDY4e = new TMarker(4, 613, 24);
   PCNDY4e->SetMarkerSize(1.3); 
   PCNDY4e->SetMarkerColor(2);
   PCNDY4e->Draw();
   TMarker* PCNDY5e = new TMarker(5, 597, 24);
   PCNDY5e->SetMarkerSize(1.3);
   PCNDY5e->SetMarkerColor(2);
   PCNDY5e->Draw();
   TMarker* PCNDY6e = new TMarker(6, 566, 24);
   PCNDY6e->SetMarkerSize(1.3);
   PCNDY6e->SetMarkerColor(2);
   PCNDY6e->Draw();
   TMarker* PCNDY7e = new TMarker(7, 477, 24);
   PCNDY7e->SetMarkerSize(1.3);
   PCNDY7e->SetMarkerColor(2);
   PCNDY7e->Draw();


   TLegend* leg = new TLegend(0.15,0.91-(4)*0.040,0.45,0.91);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextFont(43); //give the font size in pixel (instead of fraction)
   leg->SetTextSize(18); //font size
   leg->SetTextAlign(12);
   leg->AddEntry(gr1,"#bf{CMS}      #sqrt{s}=7 TeV, L=5.0fb^{-1}  #sqrt{s}=8 TeV, L=18.8fb^{-1}","p");
   leg->AddEntry(PCDY2e3New,"#bf{CMS}      #sqrt{s}=7 TeV, L=5.0fb^{-1} (2012)","p");
   leg->AddEntry(PCDY2e3,"#bf{CMS}      #sqrt{s}=7 TeV, L=5.0fb^{-1} (2011)","p");
   //leg->AddEntry(PCDY1e,"#bf{CMS}      #sqrt{s}=7 TeV, L=5.0fb^{-1} (2011 PAS)","p");  //was in the 2012AN
   
   leg->AddEntry(PADY2e,"#bf{ATLAS} #sqrt{s}=7 TeV, L=4.4fb^{-1}","p");
   leg->Draw();

//   SQRTS=78.0;
//   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS), false, 0.15, 0.995, 0.95, 0.945);
   c1->Print("massLimitSummaryPlot_Q.pdf");
   c1->Print("massLimitSummaryPlot_Q.png");
   delete c1;
   delete leg;

   c1 = new TCanvas("c1","c1",600,600);
   TH1D* frame = new TH1D("frame", "frame", 8,0.5,8.5);
   frame->SetTitle("");
   frame->SetStats(kFALSE);
   frame->GetXaxis()->SetTitle("");
   frame->GetYaxis()->SetTitle("95% C.L. lower mass limit (GeV/#font[12]{c}^{2})");
   frame->GetYaxis()->SetTitleOffset(1.9);
   frame->SetMaximum(1700);
   frame->SetMinimum(0);
   frame->Draw();

   // pair prod. stau
   frame->GetXaxis()->SetBinLabel(1, "#splitline{  stau}{#splitline{  Pair}{  prod.}}");
   TMarker* ppstau    = new TMarker(1.0, 339, 20);    ppstau->SetMarkerSize(1.3);    ppstau->SetMarkerColor(2);    ppstau->Draw();
// TMarker* PCppstau  = new TMarker(1.0, 223, 21);  PCppstau->SetMarkerSize(1.3);  PCppstau->SetMarkerColor(4);  PCppstau->Draw();
   TMarker* PCNppstau = new TMarker(1.0, 190, 24); PCNppstau->SetMarkerSize(1.3); PCNppstau->SetMarkerColor(2); PCNppstau->Draw();

   // gmsb stau
   frame->GetXaxis()->SetBinLabel(2, "#splitline{  stau}{GMSB}");
   TMarker* stau    = new TMarker(2.0, 435, 20);     stau->SetMarkerSize(1.3);     stau->SetMarkerColor(2);     stau->Draw();
// TMarker* PCstau  = new TMarker(2.0, 314, 21);   PCstau->SetMarkerSize(1.3);   PCstau->SetMarkerColor(4);   PCstau->Draw();
   TMarker* PCNstau = new TMarker(2.0, 285, 24);  PCNstau->SetMarkerSize(1.3);  PCNstau->SetMarkerColor(2);  PCNstau->Draw();
   TMarker* PAstau  = new TMarker(2.0, 300, 23);   PAstau->SetMarkerSize(1.3);   PAstau->SetMarkerColor(1);   PAstau->Draw();

   //stop
   frame->GetXaxis()->SetBinLabel(3, "#splitline{ stop}{#splitline{  ch.}{suppr.}}");
   TMarker* stopN    = new TMarker(3.0, 818, 20);      stopN->SetMarkerSize(1.3);      stopN->SetMarkerColor(2);      stopN->Draw();
// TMarker* PCstopN  = new TMarker(3.0, 626, 21);    PCstopN->SetMarkerSize(1.3);    PCstopN->SetMarkerColor(4);    PCstopN->Draw();
   TMarker* PCNstopN = new TMarker(3.0, 641, 24);   PCNstopN->SetMarkerSize(1.3);   PCNstopN->SetMarkerColor(2);   PCNstopN->Draw();
   TMarker* PAstopN  = new TMarker(3.0, 604, 23);    PAstopN->SetMarkerSize(1.3);    PAstopN->SetMarkerColor(1);    PAstopN->Draw();

   frame->GetXaxis()->SetBinLabel(4, "stop");   
   TMarker* stop    = new TMarker(4.0, 935, 20);      stop->SetMarkerSize(1.3);      stop->SetMarkerColor(2);      stop->Draw();
// TMarker* PCstop  = new TMarker(4.0, 737, 21);    PCstop->SetMarkerSize(1.3);    PCstop->SetMarkerColor(4);    PCstop->Draw();
   TMarker* PCNstop = new TMarker(4.0, 749, 24);   PCNstop->SetMarkerSize(1.3);   PCNstop->SetMarkerColor(2);   PCNstop->Draw();
   TMarker* PAstop  = new TMarker(4.0, 683, 23);    PAstop->SetMarkerSize(1.3);    PAstop->SetMarkerColor(1);    PAstop->Draw();

   //gluino
   frame->GetXaxis()->SetBinLabel(5, "#splitline{gluino}{#splitline{  ch.}{#splitline{suppr.}{(f=0.1)}}}");
   TMarker* gluinoN    = new TMarker(5.0, 1233, 20);    gluinoN->SetMarkerSize(1.3);    gluinoN->SetMarkerColor(2);    gluinoN->Draw();
// TMarker* PCgluinoN  = new TMarker(5.0, 940, 21);   PCgluinoN->SetMarkerSize(1.3);  PCgluinoN->SetMarkerColor(4);  PCgluinoN->Draw();
   TMarker* PCNgluinoN = new TMarker(5.0, 963, 24);  PCNgluinoN->SetMarkerSize(1.3); PCNgluinoN->SetMarkerColor(2); PCNgluinoN->Draw();
   TMarker* PAgluinoN  = new TMarker(5.0, 885, 23);   PAgluinoN->SetMarkerSize(1.3);  PAgluinoN->SetMarkerColor(1);  PAgluinoN->Draw();

   frame->GetXaxis()->SetBinLabel(6, "#splitline{gluino}{(f=0.1)}");
   TMarker* gluino    = new TMarker(6.0, 1322, 20);       gluino->SetMarkerSize(1.3);       gluino->SetMarkerColor(2);       gluino->Draw();
// TMarker* PCgluino  = new TMarker(6.0, 1098, 21);     PCgluino->SetMarkerSize(1.3);     PCgluino->SetMarkerColor(4);     PCgluino->Draw();
   TMarker* PCNgluino = new TMarker(6.0, 1099, 24);    PCNgluino->SetMarkerSize(1.3);    PCNgluino->SetMarkerColor(2);    PCNgluino->Draw();
   TMarker* PAgluino  = new TMarker(6.0, 985, 23);      PAgluino->SetMarkerSize(1.3);     PAgluino->SetMarkerColor(1);     PAgluino->Draw();

   frame->GetXaxis()->SetBinLabel(7, "#splitline{gluino}{(f=0.5)}");
   TMarker* gluino_f50    = new TMarker(7.0, 1276, 20);     gluino_f50->SetMarkerSize(1.3);      gluino_f50->SetMarkerColor(2);      gluino_f50->Draw();
// TMarker* PCgluino_f50  = new TMarker(7.0, 1046, 21);   PCgluino_f50->SetMarkerSize(1.3);    PCgluino_f50->SetMarkerColor(4);    PCgluino_f50->Draw();
   TMarker* PCNgluino_f50 = new TMarker(7.0, 1049, 24);  PCNgluino_f50->SetMarkerSize(1.3);   PCNgluino_f50->SetMarkerColor(2);   PCNgluino_f50->Draw();
   TMarker* PAgluino_f50  = new TMarker(7.0, 537, 22);    PAgluino_f50->SetMarkerSize(1.3);    PAgluino_f50->SetMarkerColor(1);    PAgluino_f50->Draw();


   frame->GetXaxis()->SetBinLabel(8, "#splitline{gluino}{(f=1.0)}");
   TMarker*   gluino_f100 = new TMarker(8.0, 1250, 20);    gluino_f100->SetMarkerSize(1.3);     gluino_f100->SetMarkerColor(2);     gluino_f100->Draw();
   TMarker* PAgluino_f100 = new TMarker(8.0, 530, 22);   PAgluino_f100->SetMarkerSize(1.3);   PAgluino_f100->SetMarkerColor(1);   PAgluino_f100->Draw();


   frame->GetXaxis()->LabelsOption("h");
   leg = new TLegend(0.15,0.91-5*0.040,0.45,0.91);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->SetTextFont(43); //give the font size in pixel (instead of fraction)
   leg->SetTextSize(18); //font size
   leg->SetTextAlign(12);
   leg->AddEntry(  stau ,"#bf{CMS}     5.0+18.8 fb^{-1} (7+8 TeV)","p");
   leg->AddEntry(PCNstau,"#bf{CMS}     5.0 fb^{-1} (7 TeV)","p");
//   leg->AddEntry(PCstau,"#bf{CMS}   5.0 fb^{-1} (7 TeV)","p");
   leg->AddEntry(PAstau,"#bf{ATLAS} 4.7 fb^{-1} (7 TeV)","p");
   leg->AddEntry(PAgluino_f50,"#bf{ATLAS} 37 pb^{-1} (7 TeV)","p");
   leg->Draw();

//   SQRTS=78.0;
//   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS), false, 0.15, 0.995, 0.95, 0.945);
   c1->Print("massLimitSummaryPlot_S.pdf");
   c1->Print("massLimitSummaryPlot_S.png");
   //return c1;
}
