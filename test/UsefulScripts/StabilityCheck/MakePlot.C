
#include <vector>
#include <algorithm>

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
#include "TProfile.h"
#include "TCutG.h"


#include "../../AnalysisCode/Analysis_Step1_EventLoop.C"


std::unordered_map<unsigned int, std::unordered_map<unsigned int, std::pair<double,double> > > LumiMap;
std::unordered_map<unsigned int, double> RunIntLumi;

bool LoadInstLumi(string inputFile="parsedLumi.txt")
{
   FILE* pFile = fopen(inputFile.c_str(),"r");
   if(!pFile){
      printf("Not Found: %s\n",inputFile.c_str());
      return false;
   }

   unsigned int Run, Lumi;
   double InstLumi; double InteLumi;
   while ( ! feof (pFile) ){
     fscanf(pFile,"%d:%d --> %lf / %lf\n",&Run, &Lumi, &InstLumi, &InteLumi);
     LumiMap[Run][Lumi] = std::make_pair(InstLumi, InteLumi);
     RunIntLumi[Run] += InteLumi;
   }
   fclose(pFile);
   return true;
}


std::map<unsigned int, double> RunToIntLumi;

bool LoadLumiToRun()
{
   float TotalIntLuminosity = 0;

   FILE* pFile = fopen("out.txt","r");
   if(!pFile){
      printf("Not Found: %s\n","out.txt");
      return false;
   }

   unsigned int Run; float IntLumi;
   unsigned int DeliveredLs; double DeliveredLumi;
   char Line[2048], Tmp1[2048], Tmp2[2048], Tmp3[2048];
   while ( ! feof (pFile) ){
     fscanf(pFile,"%s\n",Line);
//     printf("%s\n",Line);
     for(unsigned int i=0;Line[i]!='\0';i++){if(Line[i]==',')Line[i]=' ';} 
     sscanf(Line,"%d %s %s %s %f\n",&Run,Tmp1,Tmp2,Tmp3,&IntLumi);
     TotalIntLuminosity+= IntLumi/1000000.0;
//     printf("%6i --> %f/pb   (%s | %s | %s)\n",Run,TotalIntLuminosity,Tmp1,Tmp2,Tmp3);

     RunToIntLumi[Run] = TotalIntLuminosity;
   }
   fclose(pFile);
   return true;
}


TGraph* ConvertFromRunToIntLumi(TProfile* Object, const char* DrawOption, string YLabel, double YRange_Min=3.1, double YRange_Max=3.7){
   TGraphErrors* graph = new TGraphErrors(Object->GetXaxis()->GetNbins());
   for(int i=1;i<Object->GetXaxis()->GetNbins()+1;i++){
      int RunNumber;
      sscanf(Object->GetXaxis()->GetBinLabel(i),"%d",&RunNumber);
      graph->SetPoint(i-1, RunToIntLumi[RunNumber], Object->GetBinContent(i));
      graph->SetPointError(i-1, 0.0*RunToIntLumi[RunNumber], Object->GetBinError(i));
   }
   graph->Draw(DrawOption);
   graph->SetTitle("");
   graph->GetYaxis()->SetTitle(Object->GetYaxis()->GetTitle());
   graph->GetYaxis()->SetTitleOffset(1.10);
   graph->GetXaxis()->SetTitle("Int. Luminosity (/pb)");
   graph->GetYaxis()->SetTitle(YLabel.c_str());
   graph->SetMarkerColor(Object->GetMarkerColor());
   graph->SetMarkerStyle(Object->GetMarkerStyle());
   graph->GetXaxis()->SetNdivisions(510);
   if(YRange_Min!=YRange_Max)graph->GetYaxis()->SetRangeUser(YRange_Min,YRange_Max);
   return graph;
}

void MakedEdxPlot()
{
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.09);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);

   TCanvas* c1;
   TObject** Histos = new TObject*[10];
   std::vector<string> legend;

   TFile* InputFile = new TFile("pictures/Histos.root");

   std::vector<string> runList;
   TList* ObjList = InputFile->GetListOfKeys();
   for(int i=0;i<ObjList->GetSize();i++){
      TObject* tmp = GetObjectFromPath(InputFile,ObjList->At(i)->GetName(),false);
      if(tmp->InheritsFrom("TDirectory")){
         runList.push_back(ObjList->At(i)->GetName());
         printf("Add a new run: %s\n", ObjList->At(i)->GetName() );
      }
      delete tmp;
   }

   TProfile* SingleMu_PtProf           = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50PtProf");      
   TProfile* SingleMu_dEdxProf         = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxProf");   
   TProfile* SingleMu_dEdxMProf        = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMSProf       = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMSProf");
   TProfile* SingleMu_dEdxMPProf       = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMPProf");
   TProfile* SingleMu_dEdxMSCProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMSCProf");
   TProfile* SingleMu_dEdxMPCProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMPCProf");
   TProfile* SingleMu_dEdxMSFProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMSFProf");
   TProfile* SingleMu_dEdxMPFProf      = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50dEdxMPFProf");

   TProfile* SingleMu_NVertProf        = (TProfile*)GetObjectFromPath(InputFile, "HLT_Mu50NVertProf");

   SingleMu_NVertProf->LabelsDeflate("X");
   SingleMu_NVertProf->LabelsOption("av","X");

/*
   TFile* InputFileLumi166380 = new TFile("pictures/HistosLumi166380.root");
   TFile* InputFileLumi166512 = new TFile("pictures/HistosLumi166512.root");
   TFile* InputFileLumi167807 = new TFile("pictures/HistosLumi167807.root");
   TFile* InputFileLumi167898 = new TFile("pictures/HistosLumi167898.root");

   TProfile* SingleMu_dEdxMProfLumi166380         = (TProfile*)GetObjectFromPath(InputFileLumi166380, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMProfLumi166512         = (TProfile*)GetObjectFromPath(InputFileLumi166512, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMProfLumi167807         = (TProfile*)GetObjectFromPath(InputFileLumi167807, "HLT_Mu50dEdxMProf");
   TProfile* SingleMu_dEdxMProfLumi167898         = (TProfile*)GetObjectFromPath(InputFileLumi167898, "HLT_Mu50dEdxMProf");
*/

   if(LoadLumiToRun()){
      TLegend* leg;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graph =  ConvertFromRunToIntLumi(SingleMu_dEdxMProf  , "A*", "I_{h} (MeV/cm)");
      TGraph* graphS = ConvertFromRunToIntLumi(SingleMu_dEdxMSProf, "*" , "I_{h} (MeV/cm)");
      TGraph* graphP = ConvertFromRunToIntLumi(SingleMu_dEdxMPProf, "*" , "I_{h} (MeV/cm)");
      graphS->SetMarkerColor(2);    graphS->SetMarkerStyle(26);
      graphP->SetMarkerColor(4);    graphP->SetMarkerStyle(32);


      TF1* myfunc = new TF1("Fitgraph" ,"pol1",250,5000);  graph ->Fit(myfunc ,"QN","",250,5000); myfunc ->SetLineWidth(2); myfunc ->SetLineColor(graph ->GetMarkerColor()); myfunc ->Draw("same");
      TF1* myfuncS= new TF1("FitgraphS","pol1",250,5000);  graphS->Fit(myfuncS,"QN","",250,5000); myfuncS->SetLineWidth(2); myfuncS->SetLineColor(graphS->GetMarkerColor()); myfuncS->Draw("same");
      TF1* myfuncP= new TF1("FitgraphP","pol1",250,5000);  graphP->Fit(myfuncP,"QN","",250,5000); myfuncP->SetLineWidth(2); myfuncP->SetLineColor(graphP->GetMarkerColor()); myfuncP->Draw("same");
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip+Pixel)", myfunc ->GetChisquare()/ myfunc ->GetNDF(), myfunc ->GetParameter(0),myfunc ->GetParError(0),myfunc ->GetParameter(1),myfunc ->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip)"      , myfuncS->GetChisquare()/ myfuncS->GetNDF(), myfuncS->GetParameter(0),myfuncS->GetParError(0),myfuncS->GetParameter(1),myfuncS->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Pixel)"      , myfuncP->GetChisquare()/ myfuncP->GetNDF(), myfuncP->GetParameter(0),myfuncP->GetParError(0),myfuncP->GetParameter(1),myfuncP->GetParError(1));
      leg = new TLegend(0.79,0.92,0.79-0.20,0.92 - 3*0.05);     leg->SetFillColor(0);     leg->SetBorderSize(0);
      leg->AddEntry(graph, "dE/dx (Strip+Pixel)" ,"P");
      leg->AddEntry(graphS, "dE/dx (Strip)" ,"P");
      leg->AddEntry(graphP, "dE/dx (Pixel)" ,"P");
      leg->Draw();
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_dEdxM");
      delete c1;  delete leg;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphSC = ConvertFromRunToIntLumi(SingleMu_dEdxMSCProf, "A*", "I_{h} (MeV/cm)");
      TGraph* graphSF = ConvertFromRunToIntLumi(SingleMu_dEdxMSFProf, "*" , "I_{h} (MeV/cm)");
      graphSC->SetMarkerColor(2);    graphSC->SetMarkerStyle(26);
      graphSF->SetMarkerColor(4);    graphSF->SetMarkerStyle(32);
      TF1* myfuncSC= new TF1("FitgraphSC","pol1",250,5000);  graphSC->Fit(myfuncSC,"QN","",250,5000); myfuncSC->SetLineWidth(2); myfuncSC->SetLineColor(graphSC->GetMarkerColor()); myfuncSC->Draw("same");
      TF1* myfuncSF= new TF1("FitgraphSF","pol1",250,5000);  graphSF->Fit(myfuncSF,"QN","",250,5000); myfuncSF->SetLineWidth(2); myfuncSF->SetLineColor(graphSF->GetMarkerColor()); myfuncSF->Draw("same");
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip) |eta|<0.5", myfuncSC->GetChisquare()/ myfuncSC->GetNDF(), myfuncSC->GetParameter(0),myfuncSC->GetParError(0),myfuncSC->GetParameter(1),myfuncSC->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Strip) |eta|>1.5", myfuncSF->GetChisquare()/ myfuncSF->GetNDF(), myfuncSF->GetParameter(0),myfuncSF->GetParError(0),myfuncSF->GetParameter(1),myfuncSF->GetParError(1));
      leg = new TLegend(0.79,0.92,0.79-0.20,0.92 - 3*0.05);     leg->SetFillColor(0);     leg->SetBorderSize(0);
      leg->AddEntry(graphSC, "dE/dx (Strip) |#eta|<0.5" ,"P");
      leg->AddEntry(graphSF, "dE/dx (Strip) |#eta|>1.5"  ,"P");
      leg->Draw();
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_dEdxMS");
      delete c1; delete leg;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphPC = ConvertFromRunToIntLumi(SingleMu_dEdxMPCProf, "A*", "I_{h} (MeV/cm)");
      TGraph* graphPF = ConvertFromRunToIntLumi(SingleMu_dEdxMPFProf, "*" , "I_{h} (MeV/cm)");
      graphPC->SetMarkerColor(2);    graphPC->SetMarkerStyle(26);
      graphPF->SetMarkerColor(4);    graphPF->SetMarkerStyle(32);
      TF1* myfuncPC= new TF1("FitgraphPC","pol1",250,5000);  graphPC->Fit(myfuncPC,"QN","",250,5000); myfuncPC->SetLineWidth(2); myfuncPC->SetLineColor(graphPC->GetMarkerColor()); myfuncPC->Draw("same");
      TF1* myfuncPF= new TF1("FitgraphPF","pol1",250,5000);  graphPF->Fit(myfuncPF,"QN","",250,5000); myfuncPF->SetLineWidth(2); myfuncPF->SetLineColor(graphPF->GetMarkerColor()); myfuncPF->Draw("same");
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Pixel) |eta|<0.5", myfuncPC->GetChisquare()/ myfuncPC->GetNDF(), myfuncPC->GetParameter(0),myfuncPC->GetParError(0),myfuncPC->GetParameter(1),myfuncPC->GetParError(1));
      printf("%25s --> Chi2/ndf = %6.2f --> a=%6.2E+-%6.2E   b=%6.2E+-%6.2E\n","dE/dx (Pixel) |eta|>1.5", myfuncPF->GetChisquare()/ myfuncPF->GetNDF(), myfuncPF->GetParameter(0),myfuncPF->GetParError(0),myfuncPF->GetParameter(1),myfuncPF->GetParError(1));
      leg = new TLegend(0.79,0.92,0.79-0.20,0.92 - 3*0.05);     leg->SetFillColor(0);     leg->SetBorderSize(0);
      leg->AddEntry(graphPC, "dE/dx (Pixel) |#eta|<0.5" ,"P");
      leg->AddEntry(graphPF, "dE/dx (Pixel) |#eta|>1.5" ,"P");
      leg->Draw();
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_dEdxMP");
      delete c1; delete leg;




      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphNV = ConvertFromRunToIntLumi(SingleMu_NVertProf, "A*" , "<#Reco Vertices>",0,0);
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_Vert");
      delete c1;

      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      TGraph* graphpT = ConvertFromRunToIntLumi(SingleMu_PtProf, "A*" , "<p_{T}> (GeV)",0,0);
      SaveCanvas(c1,"pictures/","GraphdEdx_Profile_pT");
      delete c1;


   }else{
      printf("TEST TEST TEST\n");
   }


   for(int i=0;i<SingleMu_PtProf->GetXaxis()->GetNbins();i++){
//      if((i+3)%4==0)continue;
      SingleMu_PtProf->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxProf->GetXaxis()->SetBinLabel(i,"");
      SingleMu_dEdxMProf->GetXaxis()->SetBinLabel(i,"");
      SingleMu_NVertProf->GetXaxis()->SetBinLabel(i,"");
   }  


   SQRTS = 1316;
   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_NVertProf;                 legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "<#Reco Vertices>", 0,0, 0,0);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_NVert");
   delete c1;

 

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_PtProf;                    legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "p_{T} (GeV)", 0,0, 0,0);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_Pt");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_dEdxProf;                  legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "I_{as}", 0,0, 0.02,0.06);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_dEdx");
   delete c1;

   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_dEdxMProf;                  legend.push_back("SingleMu50");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "", "I_{h}", 0,0, 3.2,3.4);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_dEdxM");
   delete c1;

/*
   c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
   Histos[0] = SingleMu_dEdxMProfLumi166380;       legend.push_back("SingleMu50 - Run166380");
   Histos[1] = SingleMu_dEdxMProfLumi166512;       legend.push_back("SingleMu50 - Run166512");
   Histos[2] = SingleMu_dEdxMProfLumi167807;       legend.push_back("SingleMu50 - Run167807");
   Histos[3] = SingleMu_dEdxMProfLumi167898;       legend.push_back("SingleMu50 - Run167898");
   DrawSuperposedHistos((TH1**)Histos, legend, "E1",  "Lumi", "I_{h}", 0,0, 3.2,3.4);
   for(unsigned int i=0;i<legend.size();i++){((TProfile*)Histos[i])->SetMarkerSize(0.5);           ((TProfile*)Histos[i])->GetYaxis()->SetTitleOffset(0.9);}
   DrawLegend(Histos,legend,"","P");
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"pictures/","dEdx_Profile_dEdxMRun");
   delete c1;
*/


}

double GetMeanAndRMS(TFile* InputFile, string path, double& mean, double& rms, bool gaussianFit=false){
   double toReturn=-1;
   TH1D* histo = (TH1D*)GetObjectFromPath(InputFile, path);
   if(!histo){
      mean = -999;
      rms  = -999;
   }else{
      toReturn = histo->GetEntries();
      mean = histo->GetMean();
      rms  = histo->GetMeanError();

      if(histo->GetEntries()>10 && gaussianFit){
         rms  = histo->GetRMS(); 
         TF1* fit = new TF1("MyFit","gaus", mean-2*rms, mean+2*rms);
         fit->SetParameters(histo->GetEntries(),  histo->GetMean(), histo->GetRMS());
         histo->Fit(fit, "0QR"); 
         mean = fit->GetParameter(1);
         rms  = fit->GetParError(1);
         delete fit;
      }
      delete histo;
   }
   return toReturn;
}

void printAverageValue(std::vector<string>& runList, TFile* InputFile, string HistoName){
   printf("Average value for %s\n", HistoName.c_str());
   double mean, rms;
   for(unsigned int r=0;r<runList.size();r++){
      char name[1024]; sprintf(name, "%s/%s", runList[r].c_str(), HistoName.c_str());
      double NEntries=GetMeanAndRMS(InputFile, name, mean, rms);
      printf("   Run=%s --> average=%8.4E  NEntries=%8.4E\n", runList[r].c_str(), mean, NEntries);
   }
}

TH2D* getHIPRateGraph(std::vector<string>& runList, TFile* InputFile, string HistoName, float threshold, string HistoWeightName){
   printf("HIP value for %s\n", HistoName.c_str());

   //TGraph* toReturn = new TGraph(65000);
   TH2D* toReturn = new TH2D("h2d", "h2d", 20, -5.5, -3.5, 100, 0, 10);
   int INDEX=0;
   for(unsigned int r=0;r<runList.size();r++){
      unsigned int Run;
      sscanf(runList[r].c_str(),"%d",&Run);
      std::unordered_map<unsigned int, std::pair<double,double> >& ThisRunLumiMap = LumiMap[Run];

      //perLumi Numbers
      TDirectory* dir = (TDirectory*)InputFile->Get(runList[r].c_str());
      TList* ObjList = dir->GetListOfKeys();
      for(int i=0;i<ObjList->GetSize();i++){         
         string name = string(ObjList->At(i)->GetName());
         if(name.find(HistoName+"_ls")==0 && name.find("_ls")!=string::npos){
            unsigned int Lumi;
            sscanf(name.c_str()+name.find("_ls"),"_ls%d",&Lumi);
            double instLumi = ThisRunLumiMap[Lumi].first;
            if(instLumi<=0)continue;


            TH1D* histo = (TH1D*)GetObjectFromPath(dir,ObjList->At(i)->GetName(),false);
            if(!histo)continue;
            if(histo->GetEntries()<100)continue;

            double rate = histo->Integral(0, histo->GetXaxis()->FindBin(threshold)) / histo->Integral(0, histo->GetNbinsX()+2);
//            toReturn->SetPoint(INDEX, log10(instLumi), 100.0*rate);
            toReturn->Fill(log10(instLumi), 100.0*rate);

            //printf("%i - %i --> %f vs %f\n", Run, Lumi, log10(instLumi), 100.0*rate);
            INDEX++;
            //if(INDEX>1000)break;
         }
      }
      //if(INDEX>1000)break;
   }
//   toReturn->Set(INDEX);
   return toReturn;
}



TH1D* getHIPRate(std::vector<string>& runList, TFile* InputFile, string HistoName, float threshold, string HistoWeightName){
   printf("HIP value for %s\n", HistoName.c_str());

   TH1D* toReturn = new TH1D((HistoName+"HIPRate").c_str(), "HIP Rate", 26, 0.0, 13.0);
   toReturn->SetTitle("");
   toReturn->SetStats(kFALSE);
   toReturn->GetXaxis()->SetNdivisions(505);
   toReturn->GetXaxis()->SetTitle("");
   toReturn->GetYaxis()->SetTitleOffset(0.95);
   toReturn->GetYaxis()->SetTitle("Entries (a.u.)");
   toReturn->GetXaxis()->SetTitle("Fraction of low dEdx hits (%)");


   double mean, rms;
   for(unsigned int r=0;r<runList.size();r++){
      unsigned int Run=0; double runIntLumi = -1;  
      if(runList[r].find("MC")==string::npos){
          sscanf(runList[r].c_str(),"%d",&Run);
          runIntLumi = RunIntLumi[Run];
      }
      std::unordered_map<unsigned int, std::pair<double,double> >& ThisRunLumiMap = LumiMap[Run];


      char nameStr[1024]; sprintf(nameStr, "%s/%s", runList[r].c_str(), HistoName.c_str());
      TH1D* histo = (TH1D*)GetObjectFromPath(InputFile, nameStr);
      if(histo){
         double rate = histo->Integral(0, histo->GetXaxis()->FindBin(threshold)) / histo->Integral(0, histo->GetNbinsX()+2);
//         sprintf(nameStr, "%s/%s", runList[r].c_str(), HistoWeightName.c_str());
//         double NEntries=GetMeanAndRMS(InputFile, nameStr, mean, rms);
         printf("   Run=%s --> average=%6.2f  IntLumi=%8.4E --> %6.2f%%\n", runList[r].c_str(), mean, runIntLumi, 100.0*rate);
         //toReturn->Fill(100.0*rate, NEntries);
      }

      if(runList[r].find("MC")!=string::npos)continue; 

      //perLumi Numbers
      TDirectory* dir = (TDirectory*)InputFile->Get(runList[r].c_str());
      TList* ObjList = dir->GetListOfKeys();
      for(int i=0;i<ObjList->GetSize();i++){         
         string name = string(ObjList->At(i)->GetName());
         if(name.find(HistoName+"_ls")==0 && name.find("_ls")!=string::npos){
            unsigned int Lumi;
            sscanf(name.c_str()+name.find("_ls"),"_ls%d",&Lumi);
            double intLumi = ThisRunLumiMap[Lumi].second;
            if(intLumi<=0)continue;

            TH1D* histo = (TH1D*)GetObjectFromPath(dir,ObjList->At(i)->GetName(),false);
            if(!histo)continue;

            double rate = histo->Integral(0, histo->GetXaxis()->FindBin(threshold)) / histo->Integral(0, histo->GetNbinsX()+2);

//            sprintf(nameStr, "%s/%s%s", runList[r].c_str(), HistoWeightName.c_str(),name.c_str()+name.find("_ls"));
//            double NEntries=GetMeanAndRMS(InputFile, nameStr, mean, rms);
//            if(NEntries<100)continue;
            //printf("   Run=%s - %s --> average=%6.2f  NEntries=%8.4E --> %6.2f%%\n", runList[r].c_str(), name.c_str()+name.find("_ls"), mean, NEntries, 100.0*rate);            
            toReturn->Fill(100.0*rate, intLumi );
         }
      }
   }
   toReturn->SetBinContent(1, toReturn->GetBinContent(0) + toReturn->GetBinContent(1));
   toReturn->SetBinContent(toReturn->GetNbinsX()+1, toReturn->GetBinContent(toReturn->GetNbinsX()+2) + toReturn->GetBinContent(toReturn->GetNbinsX()+1));
   toReturn->SetBinContent(0, 0.0);
   toReturn->SetBinContent(toReturn->GetNbinsX()+2, 0.0);
   return toReturn;
}




TGraphErrors* getStabilityGraph(std::vector<string>& runList, TFile* InputFile, string HistoName, bool gaussianFit=false, unsigned int color=4, unsigned int marker=20){
   TGraphErrors* graph = new TGraphErrors(runList.size()); 
   graph->SetName(HistoName.c_str());

   double mean, rms;
   for(unsigned int r=0;r<runList.size();r++){
      char name[1024]; sprintf(name, "%s/%s", runList[r].c_str(), HistoName.c_str());
      GetMeanAndRMS(InputFile, name, mean, rms); graph->SetPoint(r, r+0.5, mean);   graph->SetPointError(r, 0, rms);
   }
   graph->SetLineWidth(2);
   graph->SetLineColor(color);
   graph->SetFillColor(1);
   graph->SetMarkerSize(0.5);
   graph->SetMarkerStyle(marker);
   graph->SetMarkerColor(color);
   return graph;
}

TGraphErrors* getEfficiencyGraph(std::vector<string>& runList, TFile* InputFile, string HistoName, bool gaussianFit=false, unsigned int color=4, unsigned int marker=20){
   TGraphErrors* graph = new TGraphErrors(runList.size()); 
   graph->SetName(HistoName.c_str());

   double mean, rms;
   for(unsigned int r=0;r<runList.size();r++){
      char name[1024]; sprintf(name, "%s/%s", runList[r].c_str(), HistoName.c_str());
      TH1D* histo = (TH1D*) GetObjectFromPath(InputFile, name);
      graph->SetPoint(r, r+0.5, histo->GetBinContent(1));   graph->SetPointError(r, 0, 0);
   }
   graph->SetLineWidth(2);
   graph->SetLineColor(color);
   graph->SetFillColor(1);
   graph->SetMarkerSize(0.5);
   graph->SetMarkerStyle(marker);
   graph->SetMarkerColor(color);
   return graph;
}

double GetMeanAndRMS2D(TFile* InputFile, string path, int XBin, double& mean, double& rms, bool gaussianFit=false){
   double toReturn=-1;
   char BinNum [1024]; sprintf(BinNum, "%d", XBin);
   TH2D* tmp = (TH2D*)GetObjectFromPath(InputFile, path);
   TH1D* histo = tmp->ProjectionY((string(tmp->GetName())+"_Y"+string(BinNum)).c_str(), XBin, XBin);

   if(!histo){
      mean = -999;
      rms  = -999;
   }else{
      toReturn = histo->GetEntries();
      mean = histo->GetMean();
      rms  = histo->GetMeanError();

      if(histo->GetEntries()>10 && gaussianFit){
         rms  = histo->GetRMS(); 
         TF1* fit = new TF1("MyFit","gaus", mean-2*rms, mean+2*rms);
         fit->SetParameters(histo->GetEntries(),  histo->GetMean(), histo->GetRMS());
         histo->Fit(fit, "0QR"); 
         mean = fit->GetParameter(1);
         rms  = fit->GetParError(1);
         delete fit;
      }
      delete tmp;
      delete histo;
   }
   return toReturn;
}

TGraphErrors* getStabilityGraphFrom2D(std::vector<string>& runList, TFile* InputFile, string HistoName, int XBin, bool gaussianFit=false, unsigned int color=4, unsigned int marker=20){
   TGraphErrors* graph = new TGraphErrors(runList.size()); 
   graph->SetName(HistoName.c_str());

   double mean, rms;
   for(unsigned int r=0;r<runList.size();r++){
      char name[1024]; sprintf(name, "%s/%s", runList[r].c_str(), HistoName.c_str());
      GetMeanAndRMS2D(InputFile, name, XBin, mean, rms); graph->SetPoint(r, r+0.5, mean);   graph->SetPointError(r, 0, rms);
   }
   graph->SetLineWidth(2);
   graph->SetLineColor(color);
   graph->SetFillColor(1);
   graph->SetMarkerSize(0.5);
   graph->SetMarkerStyle(marker);
   graph->SetMarkerColor(color);
   return graph;
}

void overlay(std::vector<string>& runList, TFile* InputFile, string HistoName, double xMin, double xMax, string savePath, double YMIN=1E-3, string YLegend="", std::vector<string>* runLegend=NULL, std::vector<string>* selRuns=NULL){
   TCanvas* c = new TCanvas("c","c,",600,600);     
   c->SetLeftMargin(0.12);
   c->SetRightMargin(0.05);

   TH1D* frame = new TH1D("frame", "frame", 1, xMin, xMax);
   frame->SetTitle("");
   frame->SetStats(kFALSE);
   frame->GetXaxis()->SetNdivisions(505);
   frame->GetXaxis()->SetTitle("");
   frame->GetYaxis()->SetTitleOffset(0.95);
   frame->GetYaxis()->SetTitle("Entries (a.u.)");
   frame->GetXaxis()->SetTitle(YLegend.c_str());
   frame->Draw("AXIS");

   double yMax=-1E100;
   double yMin= 1E100;

   std::vector<TH1D*> histoVec;
   std::vector<double> entriesVec;

   TH1D* Average=NULL;

   for(unsigned int r=0;r<runList.size();r++){
      if(runLegend==NULL && runList[r].find("MC")!=string::npos)continue;

      char name[1024]; sprintf(name, "%s/%s", runList[r].c_str(), HistoName.c_str());
      TH1D* histo = (TH1D*)GetObjectFromPath(InputFile, name);  
      if(!histo)continue;

      if(runList[r].find("MC")==string::npos){ //only for data
          if(Average==NULL){Average = (TH1D*)histo->Clone("Data Average");
          }else{ Average->Add(histo);
          }
      }

      if(selRuns==NULL || find(selRuns->begin(), selRuns->end(), runList[r])!=selRuns->end()){
         histo->SetName((runList[r]).c_str());
         entriesVec.push_back(histo->GetEntries());
         histoVec.push_back(histo);
      }
   }

   unsigned int NPlotToShow = 10;
   if(!runLegend)NPlotToShow=25;
   std::sort(entriesVec.begin(), entriesVec.end(), std::greater<double>());
   double MinEntryCut = entriesVec.size()<NPlotToShow?0:entriesVec[NPlotToShow-1];
  

   TLegend* LEG = new TLegend(0.40,0.70,0.85,0.95);      LEG->SetFillColor(0);      LEG->SetFillStyle(0);      LEG->SetBorderSize(0);    LEG->SetNColumns(2);
   for(unsigned int r=0;r<histoVec.size();r++){
      histoVec[r]->Scale(1.0/histoVec[r]->Integral());
      if(string(histoVec[r]->GetName()).find("MC")!=string::npos){
         histoVec[r]->SetLineColor(2);
         histoVec[r]->SetLineWidth(4);
         histoVec[r]->Draw("L same");
      }else{
         if(histoVec[r]->GetEntries()<MinEntryCut)continue;
         histoVec[r]->SetLineColor(gStyle->GetColorPalette(int(r*(255.0/NPlotToShow))));
         histoVec[r]->SetMarkerColor(gStyle->GetColorPalette(int(r*(255.0/NPlotToShow))));
         histoVec[r]->SetMarkerStyle(20);
         histoVec[r]->SetLineWidth(1);
         histoVec[r]->Draw("P same");
      }
      if(runLegend!=NULL){ 
         if(string(histoVec[r]->GetName()).find("MC")!=string::npos)LEG->AddEntry(histoVec[r],(*runLegend)[r].c_str()  , "L");
         else                                                       LEG->AddEntry(histoVec[r],(*runLegend)[r].c_str()  , "P");    
      }else{     LEG->AddEntry(histoVec[r], histoVec[r]->GetName(), "P");     }

      yMax = std::max(yMax, histoVec[r]->GetMaximum());
      yMin = std::min(yMin, histoVec[r]->GetMinimum());
   }

   if(Average!=NULL){
      Average->Scale(1.0/Average->Integral());
      Average->SetLineColor(1);
      Average->SetLineWidth(2);     
      Average->Draw("L same");
      LEG->AddEntry(Average, "2016 - Overall", "L"); 
   }

   frame->SetMaximum(yMax*10.0);
   frame->SetMinimum(std::max(YMIN, yMin));
   c->SetLogy(true);
   LEG->Draw();
   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c,"pictures/",savePath);
 
   delete LEG;
   delete c;
}

void MakePlot()
{
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.15);
   gStyle->SetPadRightMargin (0.03);
   gStyle->SetPadLeftMargin  (0.07);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
   gStyle->SetPalette(55);  //Rainbow --> deifne the color of all plots
   gStyle->SetNdivisions(505);

   system("mkdir -p pictures");

   // LOAD LUMI ONCE AND FOR ALL
   LoadInstLumi();

   TCanvas* c1;
   TObject** Histos = new TObject*[10];
   std::vector<string> legend;

   TFile* InputFile = new TFile("pictures/Histos.root");
   std::vector<string> runList;
   TList* ObjList = InputFile->GetListOfKeys();
   for(int i=0;i<ObjList->GetSize();i++){
      TObject* tmp = GetObjectFromPath(InputFile,ObjList->At(i)->GetName(),false);
      if(tmp->InheritsFrom("TDirectory")){
         runList.push_back(ObjList->At(i)->GetName());
         printf("Add a new run: %s\n", ObjList->At(i)->GetName() );
      }
      delete tmp;
   }
   unsigned int N= runList.size();
   std::sort(runList.begin(), runList.end());

   std::vector<string> selectedRuns;   std::vector<string> selectedLegs;
   selectedRuns.push_back("273446");   selectedLegs.push_back("R273446 <#vtx>=17.1");
   selectedRuns.push_back("273447");   selectedLegs.push_back("R273447 <#vtx>=16.5");
   selectedRuns.push_back("273448");   selectedLegs.push_back("R273448 <#vtx>=15.4");
   selectedRuns.push_back("273449");   selectedLegs.push_back("R273449 <#vtx>=14.7");
   selectedRuns.push_back("273450");   selectedLegs.push_back("R273450 <#vtx>=13.5");
   selectedRuns.push_back("274998");   selectedLegs.push_back("R274998 <#vtx>=17.7");
   selectedRuns.push_back("274999");   selectedLegs.push_back("R274999 <#vtx>=15.0");
   selectedRuns.push_back("275000");   selectedLegs.push_back("R275000 <#vtx>=13.3");
   selectedRuns.push_back("275001");   selectedLegs.push_back("R275001 <#vtx>=11.4");
   selectedRuns.push_back("MC_13TeV_DYToMuMu");  selectedLegs.push_back("MC: DY to #mu#mu");

   std::vector<string> selectedRuns275001;   std::vector<string> selectedLegs275001;
   selectedRuns275001.push_back("275001");   selectedLegs275001.push_back("R275001");
   selectedRuns275001.push_back("MC_13TeV_DYToMuMu");  selectedLegs275001.push_back("MC: DY to #mu#mu");



   TH1D* frameR = new TH1D("frameR", "frameR", N, 0, N);
   frameR->SetTitle("");
   frameR->SetStats(kFALSE);
   frameR->GetXaxis()->SetNdivisions(505);
   frameR->GetXaxis()->SetTitle("");
   frameR->GetYaxis()->SetTitleOffset(0.95);
   for(unsigned int r=0;r<N;r++){frameR->GetXaxis()->SetBinLabel(r+1, ((r+0)%2==0)?runList[r].c_str():"");}  //plot only a label every 2


   std::vector<string> triggers;
   triggers.push_back("Any");
//   triggers.push_back("HLT_Mu45_eta2p1");
//   triggers.push_back("HLT_Mu50");
//   triggers.push_back("HLT_PFMET170_NoiseCleaned");
//   triggers.push_back("HLT_PFMET170_HBHECleaned");

   std::vector<string> versions;
   versions.push_back("");
//   versions.push_back("AOD");
//   versions.push_back("FAKE");

   for(unsigned int T=0;T<triggers.size();T++){
   for(unsigned int V=0;V<versions.size();V++){
      string trigger = triggers[T];
      string version = versions[V];
      TGraphErrors *g1, *g2, *g3;
      TLegend* LEG;

      SQRTS=1316.0;

      if(V==0 && T==0){
/*         c1 = new TCanvas("c1","c1,",600,600);                      
         TH1D* frame = new TH1D("frameR", "frameR", 1, -5.5, -3.5);
         frame->SetTitle("");
         frame->SetStats(kFALSE);
         frame->GetXaxis()->SetNdivisions(505);
         frame->GetXaxis()->SetTitle("Instantaneous Luminosity Log10 (#mu b / Xing)");
         frame->GetYaxis()->SetTitleOffset(0.95);
         frame->GetYaxis()->SetTitle("HIP Rate (%)");   
         frame->SetMinimum(0.0);   
         frame->SetMaximum(10.0);
         frame->Draw("AXIS");

         c1->SetLogy(true);
         TH2D* PixelHIPG = getHIPRateGraph(runList, InputFile, trigger+"dEdxHitPixel", 2.0, trigger+"NVert");   PixelHIPG->SetMarkerColor(2);  PixelHIPG->Draw("P same");
         TH2D* StripHIPG = getHIPRateGraph(runList, InputFile, trigger+"dEdxHitStrip", 2.0, trigger+"NVert");   StripHIPG->SetMarkerColor(4);  StripHIPG->Draw("COLZ");
         DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
         SaveCanvas(c1,"pictures/","HIP_RateGraph");
         delete c1;

         c1 = new TCanvas("c1","c1,",600,600);
         c1->SetLogy(true);
         TH1D* PixelHIP = getHIPRate(runList, InputFile, trigger+"dEdxHitPixel", 2.0, trigger+"NVert");   PixelHIP->SetLineColor(2); PixelHIP->SetLineWidth(2); PixelHIP->Draw("HIST E1");
         TH1D* StripHIP = getHIPRate(runList, InputFile, trigger+"dEdxHitStrip", 2.0, trigger+"NVert");   StripHIP->SetLineColor(4); StripHIP->SetLineWidth(2); StripHIP->Draw("HIST E1 same");
         DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
         SaveCanvas(c1,"pictures/","HIP_Rate");
         delete c1;


         std::vector<string> subDets = {"PIB", "PIE", "TIB", "TID", "TOB", "TEC"};
         c1 = new TCanvas("c1","c1,",600,600);
         c1->SetLogy(true);
         for(unsigned int d=0;d<subDets.size();d++){
            if((d+1)<3){
               TH1D* PixelHIP = getHIPRate(runList, InputFile, trigger+"dEdxHitPixel"+subDets[d], 2.0, trigger+"NVert");   PixelHIP->SetLineColor(1+d); PixelHIP->SetLineWidth(2); PixelHIP->Draw(d==0?"HIST E1":"HIST E1 same");
            }else{
               TH1D* StripHIP = getHIPRate(runList, InputFile, trigger+"dEdxHitStrip"+subDets[d], 2.0, trigger+"NVert");   StripHIP->SetLineColor(1+d); StripHIP->SetLineWidth(2); StripHIP->Draw("HIST E1 same");
            }
         }
         DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
         SaveCanvas(c1,"pictures/","HIP_RateSubDet");
         delete c1;*/
      }
     



      SQRTS = 1316;
      c1 = new TCanvas("c1","c1,",1200,600);        
      frameR->GetYaxis()->SetTitle("<#Vertices>");       frameR->SetMinimum(7.0);   frameR->SetMaximum(27);  frameR->Draw("AXIS");
      g1 = getStabilityGraph(runList, InputFile, trigger+"NVert");  g1->Draw("0 P same");
      DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1,"pictures/","Summary_Profile_NVert_"+trigger+"_"+version);
      delete c1;

      SQRTS = 1316;
      c1 = new TCanvas("c1","c1,",1200,600);        
      frameR->GetYaxis()->SetTitle("<Preselection Efficiency>");       frameR->SetMinimum(0.0);   frameR->SetMaximum(0.4);  frameR->Draw("AXIS");
      g1 = getEfficiencyGraph(runList, InputFile, trigger+"PreselEff");  g1->Draw("0 P same");
      DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1,"pictures/","Summary_Profile_PreselEff_"+trigger+"_"+version);
      delete c1;

      SQRTS = 1316;
      c1 = new TCanvas("c1","c1,",1200,600);        
      frameR->GetYaxis()->SetTitle("<Presel. tracks / Num. of Events>");       frameR->SetMinimum(0.0);   frameR->SetMaximum(0.4);  frameR->Draw("AXIS");
      g1 = getEfficiencyGraph(runList, InputFile, trigger+"TPsVNEvts");  g1->Draw("0 P same");
      DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1,"pictures/","Summary_Profile_TPsVNEvts_"+trigger+"_"+version);
      delete c1;

      SQRTS = 1316;
      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      frameR->GetYaxis()->SetTitle("p_{T} (GeV)");   frameR->SetMinimum(0.0);   frameR->SetMaximum(150.0);  frameR->Draw("AXIS");
      g1 = getStabilityGraph(runList, InputFile, trigger+"Pt");  g1->Draw("0 P same");
      DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1,"pictures/","Summary_Profile_Pt_"+trigger+"_"+version);
      delete c1;

      SQRTS = 1316;
      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
      frameR->GetYaxis()->SetTitle("p_{T} (GeV)");   frameR->SetMinimum(0.0);   frameR->SetMaximum(150.0);  frameR->Draw("AXIS");
      LEG = new TLegend(0.70,0.80,0.90,0.90);      LEG->SetFillColor(0);      LEG->SetFillStyle(0);      LEG->SetBorderSize(0);
      TGraphErrors* gr1 = getStabilityGraphFrom2D(runList, InputFile, trigger+"PtEta", 1, false, kBlack);  gr1->Draw("0 P same"); LEG->AddEntry(gr1, "|#eta| < 0.9", "PE");
      TGraphErrors* gr2 = getStabilityGraphFrom2D(runList, InputFile, trigger+"PtEta", 2, false, kBlue);   gr2->Draw("0 P same"); LEG->AddEntry(gr2, "0.9 < |#eta| < 1.2", "PE");
      TGraphErrors* gr3 = getStabilityGraphFrom2D(runList, InputFile, trigger+"PtEta", 3, false, kRed);    gr3->Draw("0 P same"); LEG->AddEntry(gr3, "1.2 < |#eta| < 2.1", "PE");
//      TGraphErrors* gr4 = getStabilityGraphFrom2D(runList, InputFile, trigger+"PtEta", 4, false, kGreen);  gr4->Draw("0 P same"); LEG->AddEntry(gr4, "2.1 < |#eta| < 2.4", "PE");
      LEG->Draw();
      DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
      SaveCanvas(c1,"pictures/","Summary_Profile_PtEta_"+trigger+"_"+version);
      delete c1;


//
//      SQRTS = 1316;
//      c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
//      frameR->GetYaxis()->SetTitle("p_{T_{BPS}} (GeV)");   frameR->SetMinimum(0.0);   frameR->SetMaximum(5.0);  frameR->Draw("AXIS");
//      g1 = getStabilityGraph(runList, InputFile, trigger+"PtBPS");  g1->Draw("0 P same");
//      DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
//      SaveCanvas(c1,"pictures/","Summary_Profile_PtBPS_"+trigger+"_"+version);
//      delete c1;


      std::string dEdxVariables[] = {"dEdx", "dEdxM", "dEdxMS", "dEdxMP", "dEdxMSC", "dEdxMPC", "dEdxMSF", "dEdxMPF", "dEdxMT", "dEdxMin1", "dEdxMin2", "dEdxMin3", "dEdxMin4", "dEdxHitStrip", "dEdxHitPixel"};
      std::string dEdxLegends[]   = {"I_{as}", "I_{h}", "I_{h} StripOnly", "I_{h} PixelOnly", "I_{h} StripOnly Barrel", "I_{h} PixelOnly Barrel", "I_{h} StripOnly Endcap", "I_{h} PixelOnly Endcap", "I_{T40}", "I_{h} (Hybrid-2-15)", "I_{h} drop lowHit (20%)", "I_{h} drop lowHit (30%)", "I_{h} drop lowHit (40%)", "Strip Hit dEdx", "Pixel Hit dEdx"};
      for(unsigned int S=0;S<sizeof(dEdxVariables)/sizeof(string);S++){
         if(dEdxLegends[S].find("I_{as}")!=std::string::npos){
            //overlay(runList     , InputFile,  trigger+dEdxVariables[S]+version, 0.0, 0.5, "overlay_"+dEdxVariables[S]+version+"All", 1E-5, dEdxLegends[S].c_str()); 
            overlay(runList, InputFile,  trigger+dEdxVariables[S]+version, 0.0, 0.5, "overlay_"+dEdxVariables[S]+"_"+trigger+"_"+version      , 1E-5, dEdxLegends[S].c_str(), &selectedLegs, &selectedRuns);

            overlay(runList, InputFile,  trigger+dEdxVariables[S]+version, 0.0, 0.5, "275001overlay_"+dEdxVariables[S]+"_"+trigger+"_"+version      , 1E-5, dEdxLegends[S].c_str(), &selectedLegs275001, &selectedRuns275001);
         }else{
            //overlay(runList     , InputFile,  trigger+dEdxVariables[S]+version, 0.0, 10.0, "overlay_"+dEdxVariables[S]+version+"All", 1E-5, dEdxLegends[S].c_str());
            overlay(runList, InputFile,  trigger+dEdxVariables[S]+version, 0.0, 10.0, "overlay_"+dEdxVariables[S]+"_"+trigger+"_"+version      , 1E-5, dEdxLegends[S].c_str(), &selectedLegs, &selectedRuns);

            overlay(runList, InputFile,  trigger+dEdxVariables[S]+version, 0.0, 10.0, "275001overlay_"+dEdxVariables[S]+"_"+trigger+"_"+version      , 1E-5, dEdxLegends[S].c_str(), &selectedLegs275001, &selectedRuns275001);




         }
     }

//      overlay(runList, InputFile,  trigger+"dEdxHitPixel"+version, 0.0, 10, "overlay_dEdxHitPixel"+version+"All", 1E-4, "Pixel Hit dEdx (MeV/cm)");
      overlay(runList, InputFile,  trigger+"dEdxHitPixel"+version, 0.0, 6, string("overlay_dEdxHitPixel")+"_"+trigger+"_"+version, 1E-4, "Pixel Hit dEdx (MeV/cm)", &selectedLegs, &selectedRuns);

//      overlay(runList, InputFile,  trigger+"dEdxHitStrip"+version, 0.0, 10, "overlay_dEdxHitStrip"+version+"All", 1E-4, "Strip Hit dEdx (MeV/cm)");
      overlay(runList, InputFile,  trigger+"dEdxHitStrip"+version, 0.0, 6, string("overlay_dEdxHitStrip")+"_"+trigger+"_"+version, 1E-4, "Strip Hit dEdx (MeV/cm)", &selectedLegs, &selectedRuns);

      overlay(runList, InputFile,  trigger+"dEdxHitStrip"+version, 0.0, 6, string("275001overlay_dEdxHitStrip")+"_"+trigger+"_"+version, 1E-4, "Strip Hit dEdx (MeV/cm)", &selectedLegs275001, &selectedRuns275001);


      for(unsigned int S=0;S<sizeof(dEdxVariables)/sizeof(string);S++){
         c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
         if(dEdxLegends[S].find("I_{as}")!=std::string::npos){
            frameR->GetYaxis()->SetTitle("I_{as}");   frameR->SetMinimum(0.0);   frameR->SetMaximum(0.05);  frameR->Draw("AXIS");
         }else{
           frameR->GetYaxis()->SetTitle(dEdxLegends[S].c_str());   frameR->SetMinimum(2.5);   frameR->SetMaximum(5.0);  frameR->Draw("AXIS");
         }
         LEG = new TLegend(0.70,0.80,0.90,0.90);      LEG->SetFillColor(0);      LEG->SetFillStyle(0);      LEG->SetBorderSize(0);
         g1 = getStabilityGraph(runList, InputFile, trigger+dEdxVariables[S]      , false);         g1->Draw("0 P same");  LEG->AddEntry(g1, "Calibrated" ,"P");
         g2 = getStabilityGraph(runList, InputFile, trigger+dEdxVariables[S]+"AOD", false, 1, 24);  g2->Draw("0 P same");  LEG->AddEntry(g2, "Prompt" ,"P");
         LEG->Draw();
         DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
         SaveCanvas(c1,"pictures/","Summary_Profile_"+dEdxVariables[S]+"_"+trigger+"_"+version);
         delete c1;
     }


         c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
         frameR->GetYaxis()->SetTitle("I_{h}");   frameR->SetMinimum(2.5);   frameR->SetMaximum(5.0);  frameR->Draw("AXIS");
         LEG = new TLegend(0.70,0.80,0.90,0.90);      LEG->SetFillColor(0);      LEG->SetFillStyle(0);      LEG->SetBorderSize(0);
         g1 = getStabilityGraph(runList, InputFile, trigger+"dEdxM", false, 1, 20);            g1->Draw("0 P same");     LEG->AddEntry(g1, "drop  0%" ,"P");
         g1 = getStabilityGraph(runList, InputFile, trigger+"dEdxMin1", false, 2, 20);            g1->Draw("0 P same");     LEG->AddEntry(g1, "drop 15%" ,"P");
         g1 = getStabilityGraph(runList, InputFile, trigger+"dEdxMin2", false, 4, 20);            g1->Draw("0 P same");     LEG->AddEntry(g1, "drop 20%" ,"P");
         g1 = getStabilityGraph(runList, InputFile, trigger+"dEdxMin3", false, 8, 20);            g1->Draw("0 P same");     LEG->AddEntry(g1, "drop 30%" ,"P");
         g1 = getStabilityGraph(runList, InputFile, trigger+"dEdxMin4", false, 6, 20);            g1->Draw("0 P same");     LEG->AddEntry(g1, "drop 40%" ,"P");
         LEG->Draw();
         DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
         SaveCanvas(c1,"pictures/","Summary_Profile_AllMin_"+trigger+"_"+version);
         delete c1;


         c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
         frameR->GetYaxis()->SetTitle("I_{as}");   frameR->SetMinimum(0.0);   frameR->SetMaximum(0.05);  frameR->Draw("AXIS");
         LEG = new TLegend(0.70,0.80,0.90,0.90);      LEG->SetFillColor(0);      LEG->SetFillStyle(0);      LEG->SetBorderSize(0);
         g1 = getStabilityGraph(runList, InputFile, trigger+"dEdx", false, 1, 20);            g1->Draw("0 P same");     LEG->AddEntry(g1, "I_{as} (new)" ,"P");
         g1 = getStabilityGraph(runList, InputFile, trigger+"dEdxOld", false, 2, 20);            g1->Draw("0 P same");     LEG->AddEntry(g1, "I_{as} (old)" ,"P");
         LEG->Draw();
         DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
         SaveCanvas(c1,"pictures/","Summary_Profile_TemplateComp_"+trigger+"_"+version);
         delete c1;


      std::string TOFVariables[] = {"TOF", "TOFDT", "TOFCSC", "Vertex", "VertexDT", "VertexCSC"};
      std::string TOFLegends[]   = {"1/#beta_{TOF}", "1/#beta_{TOF DT}", "1/#beta_{TOF CSC}", "Vertex time [ns]", "Vertex time from DT [ns]", "Vertex time from CSC  [ns]" };
      for(unsigned int T=0;T<sizeof(TOFVariables)/sizeof(string);T++){
         c1 = new TCanvas("c1","c1,",1200,600);          legend.clear();
         if(TOFVariables[T].find("Vertex")!=std::string::npos){
         frameR->GetYaxis()->SetTitle(TOFLegends[T].c_str());   frameR->SetMinimum(-5.0);   frameR->SetMaximum(5.0);  frameR->Draw("AXIS");
         }else{
         frameR->GetYaxis()->SetTitle(TOFLegends[T].c_str());   frameR->SetMinimum(0.85);   frameR->SetMaximum(1.15);  frameR->Draw("AXIS");
         }
         LEG = new TLegend(0.70,0.80,0.90,0.90);      LEG->SetFillColor(0);      LEG->SetFillStyle(0);      LEG->SetBorderSize(0);
         g1 = getStabilityGraph(runList, InputFile, trigger+TOFVariables[T], false);               g1->Draw("0 P same");  LEG->AddEntry(g1, "Calibrated" ,"P");
         g2 = getStabilityGraph(runList, InputFile, trigger+TOFVariables[T]+"AOD", false, 1, 24);  g2->Draw("0 P same");  LEG->AddEntry(g2, "Prompt" ,"P");
         LEG->Draw();
         DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
         SaveCanvas(c1,"pictures/","Summary_Profile_"+TOFVariables[T]+"_"+trigger+"_"+version);
         delete c1;
      }

   }
   }//end of trigger loop

//   MakedEdxPlot();
}


