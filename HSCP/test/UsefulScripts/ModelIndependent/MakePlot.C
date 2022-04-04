

//TEST

#include <vector>

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
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "tdrstyle.C"
#include "TCutG.h"
#include "TProfile.h"

#include "../../AnalysisCode/Analysis_PlotFunction.h"
#include "../../AnalysisCode/Analysis_Samples.h"


void DumpHistoToASCIIFile(string fileName, TH3* Histo){
   FILE* pFile = fopen(fileName.c_str(), "w");

   for(int x=1;x<=Histo->GetNbinsX(); x++){
   for(int y=1;y<=Histo->GetNbinsY(); y++){
   for(int z=1;z<=Histo->GetNbinsZ(); z++){
      fprintf(pFile, "%f<$p_T$<%f\t%f<$\\beta$<%f\t%f<$|\\eta|$<%f --> Bin=$%g\\pm%g$\n", Histo->GetXaxis()->GetBinLowEdge(x), Histo->GetXaxis()->GetBinUpEdge(x), Histo->GetYaxis()->GetBinLowEdge(y), Histo->GetYaxis()->GetBinUpEdge(y), Histo->GetZaxis()->GetBinLowEdge(z), Histo->GetZaxis()->GetBinUpEdge(z), Histo->GetBinContent(x, y, z), Histo->GetBinError(x, y, z));
   }}}
   fclose(pFile);   
}

void DumpManyHistoToASCIIFile(string fileName, TH3* Histo1, TH3* Histo2,TH3* Histo3,TH3* Histo4,TH3* Histo5, TH3* Histo6){
   FILE* pFile = fopen(fileName.c_str(), "w");

//\begin{table}
// \begin{center}
//  \topcaption{TabulatedValues\FIXME{write a legend}.
//     \label{tab:TabulatedValues}}
//  \begin{tabular}{|c|c|c|l|l|l|l|l|l|} \hline


   for(int x=1;x<=Histo1->GetNbinsX(); x++){
   for(int y=1;y<=Histo1->GetNbinsY(); y++){
   for(int z=1;z<=Histo1->GetNbinsZ(); z++){
      if(Histo1->GetBinContent(x, y, z)<=0)continue;
      if(Histo2->GetBinContent(x, y, z)<=0)continue;
      fprintf(pFile, "$%4.0f<p_T<%4.0f$ & $%4.2f<\\beta<%4.2f$ & $%5.2f<|\\eta|<%5.2f$ & $%g$ & $%g\\pm%g$ & $%g\\pm%g$ & $%g\\pm%g$ & $%g\\pm%g$ & $%g\\pm%g$ \\\\ \n", Histo1->GetXaxis()->GetBinLowEdge(x), Histo1->GetXaxis()->GetBinUpEdge(x), Histo1->GetYaxis()->GetBinLowEdge(y), Histo1->GetYaxis()->GetBinUpEdge(y), Histo1->GetZaxis()->GetBinLowEdge(z), Histo1->GetZaxis()->GetBinUpEdge(z), Histo1->GetBinContent(x, y, z), Histo2->GetBinContent(x, y, z), Histo2->GetBinError(x, y, z), Histo3->GetBinContent(x, y, z), Histo3->GetBinError(x, y, z), Histo4->GetBinContent(x, y, z), Histo4->GetBinError(x, y, z), Histo5->GetBinContent(x, y, z), Histo5->GetBinError(x, y, z), Histo6->GetBinContent(x, y, z), Histo6->GetBinError(x, y, z));
   }}}

//  \end{tabular}
// \end{center}
//\end{table}

   fclose(pFile);   
}


void MakePlot()
{
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.05);//0.06
   gStyle->SetPadBottomMargin(0.08);
   gStyle->SetPadRightMargin (0.10);
   gStyle->SetPadLeftMargin  (0.10);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);
   gStyle->SetPaintTextFormat("4.2f");
   gStyle->SetOptStat("");

   TCanvas* c1;
   TObject** Histos = new TObject*[10];
   std::vector<string> legend;



   TFile* InputFile = new TFile("pictures/Histos.root", "UPDATE");
   std::vector<string> histoNames;
   histoNames.push_back("Beta_Gen");
   histoNames.push_back("Beta_GenChaged");
   histoNames.push_back("Beta_Triggered");
   histoNames.push_back("Beta_Skimmed");
   histoNames.push_back("Beta_Matched");
   histoNames.push_back("Beta_Preselected");
   histoNames.push_back("Beta_SelectedP");
   histoNames.push_back("Beta_SelectedI");
   histoNames.push_back("Beta_SelectedT");
   histoNames.push_back("Beta_SelectedM0");
   histoNames.push_back("Beta_SelectedM1");
   histoNames.push_back("Beta_SelectedM2");
   histoNames.push_back("Beta_SelectedM3");
   histoNames.push_back("Beta_SelectedM4");
   histoNames.push_back("Beta_SelectedM5");

   TH3F** Beta_Map         = new TH3F*[histoNames.size()];
   for(unsigned int i=0;i<histoNames.size();i++){
      Beta_Map[i] = (TH3F*)GetObjectFromPath(InputFile, histoNames[i]);

      if(i>0){//if not GEN MAP
        Beta_Map[i]->Divide(Beta_Map[0]); //normalize to NGen
      }
   }

   for(unsigned int i=0;i<histoNames.size();i++){
      if(i>2){//normalized to efficiency of passing cuts after having passed the trigger
         Beta_Map[i]->Divide(Beta_Map[2]); 
      }  
      if(i!=0){//efficiencies
         Beta_Map[i]->SetMaximum(1.0); 
      }  
      Beta_Map[i]->Write(TString("Norm_")+Beta_Map[i]->GetName());
   }


//Rebin to accomodate the large number of pT bin
   for(unsigned int i=0;i<histoNames.size();i++){
      Beta_Map[i]->Rebin3D(8,1,1, "");
      Beta_Map[i]->Scale(1.0/8);
   }


  int NBins = std::min(Beta_Map[0]->GetNbinsX(), 12);
  int NCol = 2;
  int NRow = NBins/NCol + (NBins%NCol==0?0:1);
  printf("Bins = %i NCol=%i NRow=%i\n",NBins, NCol, NRow);
  TH2F* frame = new TH2F("frame", "frame;#beta;|#eta|", 10, 0.0, 1.0, 10, 0.0, 2.1);

  for(unsigned int i=0;i<histoNames.size();i++){
      TCanvas* c1 = new TCanvas("c1","c1",600*NCol,600*NRow);
/*
      TPaveText* T = new TPaveText(0.0, 0.95, 1.0, 1.0);
      T->SetTextFont(43);  //give the font size in pixel (instead of fraction)
      T->SetTextSize(41);  //font size
      T->SetFillColor(0);
      T->SetTextAlign(22);
      T->AddText("CMS Preliminary  -  #sqrt{s} = 8 TeV  -  Simulation");
      T->Draw("same");
*/

/*
      TPaveText* T1 = new TPaveText(0.50, 0.95, 0.95, 1.0);
      T1->SetTextFont(43); T1->SetTextSize(45);   T1->SetTextAlign(32);
      T1->SetFillColor(0); T1->SetFillStyle(0);   T1->SetBorderSize(0);
      T1->AddText("Simulation (8 TeV)");
      T1->Draw("same");

      TPaveText* T2 = new TPaveText(0.05, 0.95, 0.50, 1.0);
      T2->SetTextFont(63); T2->SetTextSize(54);   T2->SetTextAlign(12);
      T2->SetFillColor(0); T2->SetFillStyle(0);   T2->SetBorderSize(0);
      T2->AddText("CMS");
      T2->Draw("same");
*/
//      TPaveText* T3 = new TPaveText(0.15, 0.95, 0.50, 1.0);
//      T3->SetTextFont(53); T3->SetTextSize(45);   T3->SetTextAlign(12);
//      T3->SetFillColor(0); T3->SetFillStyle(0);   T3->SetBorderSize(0);
//      T3->AddText("simulation");
//      T3->Draw("same");

      c1->cd(0);
      TPad* p1 = new TPad("p1", "p1", 0.0, 0.0, 1.0, 1.00); p1->Draw();
      p1->Divide(NCol,NRow);
      std::vector<TObject*> toBeDeleted;             
      for(int x=1;x<=NBins;x++){
         (p1->cd(x))->SetLogz(true);
         (p1->cd(x))->SetRightMargin(0.18);//0.14);
         (p1->cd(x))->SetTopMargin(0.06);
         (p1->cd(x))->SetBottomMargin(0.10);

         char Buffer[256]; sprintf(Buffer, "p_{T}=[%i, %i] GeV/c", (int)Beta_Map[i]->GetXaxis()->GetBinLowEdge(x), (int)Beta_Map[i]->GetXaxis()->GetBinUpEdge(x));
         TPaveText* T = new TPaveText(0.05, 0.96, 0.90, 1.0, "NDC");
         T->SetTextFont(43); T->SetTextSize(23);   T->SetTextAlign(12);
         T->SetFillColor(0); T->SetFillStyle(0);   T->SetBorderSize(0);
         T->AddText(Buffer);
         Beta_Map[i]->GetXaxis()->SetRange(x,x);
         TH2F* proj = (TH2F*)Beta_Map[i]->Project3D((TString(histoNames[i].c_str())+Buffer)+"_zy");

         frame->Draw("Axis");
         if(i==0){proj->GetZaxis()->SetTitle("#Events");  proj->SetMaximum(1.5E4); proj->SetMinimum(9E1); }
         if(i==2){proj->GetZaxis()->SetTitle("P^{on}(k)");  proj->SetMaximum(1); proj->SetMinimum(1E-2);}
         if(i>=9){
                 proj->SetMaximum(1); proj->SetMinimum(1E-4);
                 sprintf(Buffer, "P^{off}(M#kern[-0.1]{_{req}} = %i GeV/c^{2}, k)", (i-9)*100);
                 proj->GetZaxis()->SetTitle(Buffer);
         }

         proj->GetXaxis()->SetTitle("#beta");
         proj->GetYaxis()->SetTitle("|#eta|");
         proj->GetXaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
         proj->GetXaxis()->SetLabelSize(22); //font size
         proj->GetXaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
         proj->GetXaxis()->SetTitleSize(22); //font size
         proj->GetYaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
         proj->GetYaxis()->SetLabelSize(22); //font size
         proj->GetYaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
         proj->GetYaxis()->SetTitleSize(22); //font size
         proj->GetZaxis()->SetLabelFont(43); //give the font size in pixel (instead of fraction)
         proj->GetZaxis()->SetLabelSize(22); //font size
         proj->GetZaxis()->SetTitleFont(43); //give the font size in pixel (instead of fraction)
         proj->GetZaxis()->SetTitleSize(22); //font size
         proj->GetXaxis()->SetTitleOffset(1.5);
         proj->GetYaxis()->SetTitleOffset(1.5);
         proj->GetZaxis()->SetTitleOffset(5.0);

         proj->Draw("same COLZ text E");
         proj->Rebin2D(2,1);  if(i>0)proj->Scale(1.0/2.0);
         T->Draw("same");
         toBeDeleted.push_back(proj);
         toBeDeleted.push_back(T);
         frame->Draw("Axis same");

         TPaveText* T2 = new TPaveText(0.12, 0.865, 0.24, 0.92, "NDC");
         T2->SetTextFont(63); T2->SetTextSize(30);   T2->SetTextAlign(22);
         T2->SetFillColor(0); T2->SetFillStyle(1001);   T2->SetBorderSize(0);
         T2->AddText("CMS");
         T2->Draw("same");

         TPaveText* T4 = new TPaveText(0.50, 0.96, 0.90, 1.0, "NDC");
         T4->SetTextFont(43); T4->SetTextSize(23);   T4->SetTextAlign(32);
         T4->SetFillColor(0); T4->SetFillStyle(0);   T4->SetBorderSize(0);
         T4->AddText("Simulation (8 TeV)");
         T4->Draw("same");

      }
      c1->cd();


      char plotindex[255];sprintf(plotindex,"%02i_",i);
      c1->SaveAs((TString("pictures/")+plotindex)+histoNames[i]+".png");
      c1->SaveAs((TString("pictures/")+plotindex)+histoNames[i]+".pdf");
      c1->SaveAs((TString("pictures/")+plotindex)+histoNames[i]+".C");
//      for(unsigned int d=0;d<toBeDeleted.size();d++){delete toBeDeleted[d];}
      delete c1;
   }
   InputFile->Close();


   //CLEANUP THE ROOT FILE FOR PUBLIC USAGE
   TFile* InputFileIn  = new TFile("pictures/Histos.root", "READ");
   TFile* InputFileOut = new TFile("pictures/EXO-13-006-ProbabilityMaps.root", "RECREATE");

   TH3F* histo = NULL;

   histo = (TH3F*)GetObjectFromPath(InputFileIn, "Beta_Gen");
   histo->SetTitle("Number of single particle event generated");
   histo->SetName("NGenEvents");
   DumpHistoToASCIIFile("EXO13006_ASCII_TABLES_NGenEvents.txt", histo);
   histo->Write();

   histo = (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_Triggered");
   histo->SetTitle("Probability to pass single muon trigger requirements");
   histo->SetName("P_On(k)");
   DumpHistoToASCIIFile("EXO13006_ASCII_TABLES_POn.txt", histo);
   histo->Write();

   histo = (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM0");
   histo->SetTitle("Probability to pass offline requirements with M_{req} = 0 GeV/c^{2}");
   histo->SetName("P_Off(k,0GeV)");
   DumpHistoToASCIIFile("EXO13006_ASCII_TABLES_POff_MReq000.txt", histo);
   histo->Write();

   histo = (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM1");
   histo->SetTitle("Probability to pass offline requirements with M_{req} = 100 GeV/c^{2}");
   histo->SetName("P_Off(k,100GeV)");
   DumpHistoToASCIIFile("EXO13006_ASCII_TABLES_POff_MReq100.txt", histo);
   histo->Write();

   histo = (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM2");
   histo->SetTitle("Probability to pass offline requirements with M_{req} = 200 GeV/c^{2}");
   histo->SetName("P_Off(k,200GeV)");
   DumpHistoToASCIIFile("EXO13006_ASCII_TABLES_POff_MReq200.txt", histo);
   histo->Write();

   histo = (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM3");
   histo->SetTitle("Probability to pass offline requirements with M_{req} = 300 GeV/c^{2}");
   histo->SetName("P_Off(K,300GeV)");
   DumpHistoToASCIIFile("EXO13006_ASCII_TABLES_POff_MReq300.txt", histo);
   histo->Write();



   DumpManyHistoToASCIIFile("EXO13006_ASCII_TABLES_ALL.txt", (TH3F*)GetObjectFromPath(InputFileIn, "Beta_Gen"), (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_Triggered"), (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM0"), (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM1"), (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM2"), (TH3F*)GetObjectFromPath(InputFileIn, "Norm_Beta_SelectedM3") );

   InputFileOut->Close();

}


