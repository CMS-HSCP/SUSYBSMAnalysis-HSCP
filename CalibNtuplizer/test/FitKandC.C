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
//#include "../../AnalysisCode/tdrstyle.C"
//#include "../../AnalysisCode/Analysis_PlotFunction.h"

// inspiré de https://github.com/CMS-HSCP/SUSYBSMAnalysis-HSCP/blob/Run2_2016/test/UsefulScripts/DeDxStudy/MakePlot.C


using namespace std;

void ExtractConstants(TH2D* input1, TH2D* input2, double* K, double* C, double* Kerr, double* Cerr, double MinRange = 1.0, double MaxRange = 1.6, double MassCenter = 1.875, double LeftMassMargin = 0.2, double RightMassMargin = 0.8, double yPionMax=4.2); // by default use protons
double GetMass(double P, double I, double* K, double* C);
void SaveCanvas(TCanvas* c, std::string path, std::string name, bool OnlyPPNG=false);
TF1* GetMassLine(double M, double K, double C, bool left=false);
void ExtracttheCConstants(TH2D* inputlong, double* C,  double* Cerr) ;
void FitCalone(TString filename);
void FitKandC_sig();
void ExtractConstants_sig(TH2D* input, double* K, double* C, double* Kerr, double* Cerr, double MinRange = 1.0, double MaxRange = 1.6, double MassCenter = 1.875); // by default use protons
void ExtractProj(TH2D* inputlong, TH2D* inputshort,double* K, double* C, 
      double MinRange, double MaxRange, double MassCenter, double LeftMassMargin, double RightMassMargin, double yPionMax,
      TH1D* FitResultpion,  TH1D* FitResult);
void compaFits(TString filename1, TString filename2, double Kval, double Cval);
void TestExtractConstants(TH2D* input1, TH2D* input2, double* K, double* C, double* N, double* Kerr, double* Cerr, double* Nerr, double MinRange = 1.0, double MaxRange = 1.6, double MassCenter = 1.875, double LeftMassMargin = 0.2, double RightMassMargin = 0.8, double yPionMax=4.2); // by default use protons
TF1* TestGetMassLine(double M, double K, double C, double N, bool left=false);
Double_t TestfitSpecial(Double_t *x, Double_t *par);
Double_t TestfitSpecial2(Double_t *x, Double_t *par);

void FitKandC(TString filename){
   double K_ = 3.0; double Kerr_ = 0.001; double C_ = 3.; double Cerr_ = 0.001;
   double Ktmp = K_ ; double Ctmp = C_; double KerrTmp = Kerr_; double CerrTmp = Cerr_;
   bool uncorr = true;
   TFile* myfile;
   myfile = new TFile (filename);

   TH2D* HdedxVsP_1fit;
   TH2D* HdedxVsP_2fit;
   bool pixelOnly = false;
   myfile->cd();
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("HdedxVsP");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("HdedxVsP2");
//    // Ih 15% drop
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdXVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXVsP_lowp2");
//    // Ih 15% drop, strip only
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdXstripVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXstripVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXstripVsP_lowp2");
//    // Ih no drop
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0VsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0VsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0VsP_lowp2");
//    // Ih no drop, no L1




   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp");
   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp2");
//   HdedxVsP_2fit->Rebin2D(2,1);





//    // Ih no drop, strip only
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0stripVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0stripVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0stripVsP_lowp2");
//    // Ih 15% high drop, no L1 pix
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdXHDnoL1VsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXHDnoL1VsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXHDnoL1VsP_lowp2");
//    // Ih 15% drop, pix only
//    pixelOnly=true;
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdXpixVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXpixVsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXpixVsP_lowp2");
//    // Ih no drop, pix only, no L1
//    pixelOnly=true;
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0pixnoL1VsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0pixnoL1VsP_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0pixnoL1VsP_lowp2");
//    // Ih no drop, no L1 eta bins
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_eta1_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_eta1_lowp2");
//    // Ih no drop, no L1 PU bins
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_pu2_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_pu2_lowp2");

//    // Ih no drop, strip only, coupure 
//   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0stripVsP_c5_lowp");
//   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0stripVsP_c5_lowp2");

/*
*/

//   ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.5, 3., 0.938, 0.2, 0.8, 4.2);
//   ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.5, 1.2, 0.938, 0.2, 0.8, 4.2);
//   ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.8, 1.2, 0.938, 0.2, 0.8, 4.2);
//   ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.75, 1.5, 0.938, 0.2, 0.8, 4.2);
//     if (!pixelOnly) ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.6, 1.5, 0.938, 0.2, 0.8, 4.2);
     if (!pixelOnly) ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.5, 1.5, 0.938, 0.2, 0.5, 4.2);
//     if (!pixelOnly) ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.5, 1.3, 0.938, 0.2, 0.5, 4.2);
  // for Pixel only 
     else ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.6, 1., 0.938, 0.1, 0.8, 4.2);

     // for eta3
//     ExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 1.5, 2.5, 0.938, 0.1, 0.8, 4.2);

   

   TCanvas* c2 = new TCanvas("c2", "c2", 800,600);
   c2->cd();
   TF1* PionLine = GetMassLine(0.140, Ktmp, Ctmp, false);
   PionLine->SetLineColor(1);
   PionLine->SetLineWidth(2);
   TF1* KaonLine = GetMassLine(0.494, Ktmp, Ctmp, false);
   KaonLine->SetLineColor(1);
   KaonLine->SetLineWidth(2);
   TF1* ProtonLine = GetMassLine(0.938,Ktmp, Ctmp, false);
   ProtonLine->SetLineColor(1);
   ProtonLine->SetLineWidth(2);

   TLine* linepion = new TLine(2.,Ctmp,25., Ctmp);
   linepion->SetLineColor(1);
   linepion->SetLineWidth(2);

//   HdedxVsP_2fit->Rebin2D(5,10);
   HdedxVsP_1fit->Draw("colz");
   PionLine->Draw("same");
   KaonLine->Draw("same");
   ProtonLine->Draw("same");
   linepion->Draw("same");
   SaveCanvas(c2, "dirtest/", "dedxVsP_lines");
   

  TCanvas* c3 = new TCanvas("c3", "c3", 800,600);
   c3->cd();
   HdedxVsP_2fit->Draw("colz");
   PionLine->Draw("same");
   KaonLine->Draw("same");
   ProtonLine->Draw("same");
   linepion->Draw("same");
   SaveCanvas(c3, "dirtest/", "dedxVsP_long");
   
  gStyle->SetOptStat(0);
  c2->cd();
  c2->SetLogz(1);
  SaveCanvas(c2, "dirtest/", "dedxVsP_lines_logZ"); 
   c3->cd();
  c3->SetLogz(1);
   SaveCanvas(c3, "dirtest/", "dedxVsP_long_logZ");

   

}
void FitKandC_sig(){
   double K_ = 1.6; double Kerr_ = 0.001; double C_ = 4.3; double Cerr_ = 0.001;
   double Ktmp = K_ ; double Ctmp = C_; double KerrTmp = Kerr_; double CerrTmp = Cerr_;
   TFile* myfile;
        myfile = new TFile ("testMass2400.root");
        cout << " opening testMass2400.root" << endl;

   TH2D* HdedxVsP_2fit;
   myfile->cd();
   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdXVsP");
   ExtractConstants_sig (HdedxVsP_2fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 500., 3000., 2400);

   TCanvas* c2 = new TCanvas("c2", "c2", 800,600);
   c2->cd();
   TF1* SigLine = GetMassLine(2.4, Ktmp, Ctmp, false);
   SigLine->SetLineColor(1);
   SigLine->SetLineWidth(2);

   HdedxVsP_2fit->Draw("colz");
   SigLine->Draw("same");
   SaveCanvas(c2, "dirtest/", "SigdedxVsP_lines");
   


   

}
void FitCalone(TString filename){
   double C_ = 3.; double Cerr_ = 0.001;
   double Ctmp = C_;  double CerrTmp = Cerr_;
   TFile* myfile;
   myfile = new TFile (filename);

   TH2D* HdedxVsP_2fit;
   myfile->cd();
   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp2");
   HdedxVsP_2fit->Rebin2D(2,1);

   ExtracttheCConstants (HdedxVsP_2fit, &Ctmp, &CerrTmp);


   TCanvas* c2 = new TCanvas("c2", "c2", 800,600);
   c2->cd();

   TLine* linepion = new TLine(2.,Ctmp,25., Ctmp);
   linepion->SetLineColor(1);
   linepion->SetLineWidth(2);

   HdedxVsP_2fit->Draw("colz");
   linepion->Draw("same");
   SaveCanvas(c2, "dirtest/", "test_LineC");

}


void ExtractConstants(TH2D* inputlong, TH2D* inputshort,double* K, double* C, double* Kerr, double* Cerr,
      double MinRange, double MaxRange, double MassCenter, double LeftMassMargin, double RightMassMargin, double yPionMax)
{
       char buffer[2048];
       bool hasConverged = false;

//       for(unsigned int loop=0;loop<5 and !hasConverged; loop++){
//       for(unsigned int loop=0;loop<8 and !hasConverged; loop++){
       for(unsigned int loop=0;loop<20 and !hasConverged; loop++){
	      TH2D* inputnew = (TH2D*)inputshort->Clone("tempTH2D");
	      TH2D* inputnewPion = (TH2D*)inputlong->Clone("tempTH2D2");
//	      inputnew->Rebin2D(5,10);
	      for(int x=0;x<=inputnew->GetNbinsX()+1;x++){
	      for(int y=0;y<=inputnew->GetNbinsY()+1;y++){
		double Mass = GetMass(inputnew->GetXaxis()->GetBinCenter(x),inputnew->GetYaxis()->GetBinCenter(y), K, C);
		if(isnan (float(Mass)) || Mass<MassCenter-(LeftMassMargin) || Mass>MassCenter+RightMassMargin){
		  inputnew->SetBinContent(x,y,0);        
		  inputnew->SetBinError(x,y,0);        
		  //cout<<x<<"   "<<y<<endl;
		}
	      }}
	      for(int x=0;x<=inputnewPion->GetNbinsX()+1;x++){
   	        for(int y=0;y<=inputnewPion->GetNbinsY()+1;y++){
//                    if (inputnewPion->GetXaxis()->GetBinCenter(x)<10) inputnewPion->SetBinContent(x,y,0);
                    if (inputnewPion->GetXaxis()->GetBinCenter(x)<2.5) {
                          inputnewPion->SetBinContent(x,y,0);
                          inputnewPion->SetBinError(x,y,0);
                    }
//                    if (inputnewPion->GetXaxis()->GetBinCenter(x)>10) inputnewPion->SetBinContent(x,y,0);
	  	    if (inputnewPion->GetYaxis()->GetBinCenter(y)<2 || inputnewPion->GetYaxis()->GetBinCenter(y)>yPionMax) {
                          inputnewPion->SetBinContent(x,y,0);
                          inputnewPion->SetBinError(x,y,0);
                    }
	      }}
//	      inputnewPion->Rebin2D(10,10);

	      

	      TCanvas* c1 = new TCanvas("c1", "c1", 800,600);
	      c1->SetLogz(true);
	      inputnew->SetStats(kFALSE);
	      inputnew->GetXaxis()->SetTitle("track momentum (GeV)");
	      inputnew->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
	      inputnew->SetAxisRange(0,25,"X");
	      inputnew->SetAxisRange(0,15,"Y");
	      inputnew->Draw("COLZ");

	//      KaonLine->Draw("same");
	//      ProtonLine->Draw("same");
	//      DeuteronLine->Draw("same");
	//      TritonLine->Draw("same");

	      SaveCanvas(c1, "dirtest/", "dedxVsP");


	      inputnewPion->SetStats(kFALSE);
	      inputnewPion->GetXaxis()->SetTitle("track momentum (GeV)");
	      inputnewPion->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
	      inputnewPion->SetAxisRange(0,25,"X");
	      inputnewPion->SetAxisRange(0,15,"Y");
              inputnewPion->Draw("COLZ");

	      SaveCanvas(c1, "dirtest/", "dedxVsP_pion");

	      delete c1;
               TH1D* FitResult = new TH1D("FitResult"       , "FitResult"      ,inputnew->GetXaxis()->GetNbins(),inputnew->GetXaxis()->GetXmin(),inputnew->GetXaxis()->GetXmax());
	       FitResult->SetTitle("");
	       FitResult->SetStats(kFALSE);  
	       FitResult->GetXaxis()->SetTitle("P [GeV]");
	       FitResult->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
	       FitResult->GetYaxis()->SetTitleOffset(1.20);
	       FitResult->Reset();     

	       TH1D* FitResultPion = new TH1D("FitResultPion", "FitResultPion" ,inputnewPion->GetXaxis()->GetNbins(),inputnewPion->GetXaxis()->GetXmin(),inputnewPion->GetXaxis()->GetXmax());
               FitResultPion->SetTitle("");
               FitResultPion->SetStats(kFALSE);
               FitResultPion->GetXaxis()->SetTitle("P [GeV]");
               FitResultPion->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
               FitResultPion->GetYaxis()->SetTitleOffset(1.20);
               FitResultPion->Reset();

               for(int x=1;x<inputnewPion->GetNbinsX();x++){
		  TH1D* ProjectionPion = (TH1D*)(inputnewPion->ProjectionY("proj",x,x))->Clone();
                  if(ProjectionPion->Integral()<100)continue;
                  ProjectionPion->SetAxisRange(0.1,25,"X");
//                  ProjectionPion->Sumw2();
                  ProjectionPion->Scale(1.0/ProjectionPion->Integral());

                  TF1* mygausPion = new TF1("mygausPion","gaus", 2., 15);
                  ProjectionPion->Fit("mygausPion","Q0 RME");
                  FitResultPion->SetBinContent(x, mygausPion->GetParameter(1));
		  cout<< " mygausPion     " <<x<<"  "<<mygausPion->GetParameter(1)<<endl;
		  FitResultPion->SetBinError  (x, mygausPion->GetParError (1));
               }

	       for(int x=1;x<inputnew->GetXaxis()->FindBin(1.5);x++){
		  double P       = inputnew->GetXaxis()->GetBinCenter(x);
	    
		  TH1D* Projection = (TH1D*)(inputnew->ProjectionY("proj",x,x))->Clone();
		  if(Projection->Integral()<100)continue;
		  Projection->SetAxisRange(0.1,25,"X");
//		  Projection->Sumw2();
		  Projection->Scale(1.0/Projection->Integral());

		  TF1* mygaus = new TF1("mygaus","gaus", 2.5, 15);
		  Projection->Fit("mygaus","Q0 RME");
		  double chiFromFit  = (mygaus->GetChisquare())/(mygaus->GetNDF());
		  FitResult->SetBinContent(x, mygaus->GetParameter(1));
		  FitResult->SetBinError  (x, mygaus->GetParError (1));
		  mygaus->SetLineColor(2);
		  mygaus->SetLineWidth(2);
		  cout<<x<<"  "<<mygaus->GetParameter(1)<<endl;

		  /*
		  c1  = new TCanvas("canvas", "canvas", 600,600);
		  Projection->Draw();
		  Projection->SetTitle("");
		  Projection->SetStats(kFALSE);
		  Projection->GetXaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
		  Projection->GetYaxis()->SetTitle("#Entries");
		  Projection->GetYaxis()->SetTitleOffset(1.30);
		  Projection->SetAxisRange(1E-5,1.0,"Y");

		  mygaus->Draw("same");

		  TPaveText* stt = new TPaveText(0.55,0.82,0.79,0.92, "NDC");
		  stt->SetFillColor(0);
		  stt->SetTextAlign(31);
		  sprintf(buffer,"Proton  #mu:%5.1fMeV/cm",mygaus->GetParameter(1));      stt->AddText(buffer);
		  sprintf(buffer,"Proton  #sigma:%5.1fMeV/cm",mygaus->GetParameter(2));      stt->AddText(buffer);
		  stt->Draw("same");

		  //std::cout << "P = " << P << "  --> Proton dE/dx = " << mygaus->GetParameter(1) << endl;
		  c1->SetLogy(true);
		  sprintf(buffer,"%sProjectionFit_P%03i_%03i","fit/",(int)(100*FitResult->GetXaxis()->GetBinLowEdge(x)),(int)(100*FitResult->GetXaxis()->GetBinUpEdge(x)) );
		  if(P>=MinRange && P<=MaxRange){SaveCanvas(c1,"dirtest/",buffer);}
		  delete c1;
                  delete stt;
                  */
                  delete Projection;
                  delete mygaus;
	       }

	       c1  = new TCanvas("canvas", "canvas", 800,600);
	       FitResult->SetAxisRange(0,25,"X");
	       FitResult->SetAxisRange(0,15,"Y");
	       FitResult->Draw("");

	       TLine* line1 = new TLine(MinRange, FitResult->GetMinimum(), MinRange, FitResult->GetMaximum());
	       line1->SetLineWidth(2);
	       line1->SetLineStyle(2);
	       line1->Draw();

	       TLine* line2 = new TLine(MaxRange, FitResult->GetMinimum(), MaxRange, FitResult->GetMaximum());
	       line2->SetLineWidth(2);
	       line2->SetLineStyle(2);
	       line2->Draw();


	       double prevConstants [] = {*K, *Kerr, *C, *Cerr};

	       TF1* fitC =  new TF1("fitC","[0]", 1,25);
	       fitC->SetParName(0,"C");
	       fitC->SetParameter(0, 3.);
	       fitC->SetParLimits(0, 2,5);
	       fitC->SetLineWidth(2);
               fitC->SetLineColor(2);
               FitResultPion->Fit("fitC", "M R E I 0");
//               fitC->SetRange(5.,25);
               fitC->SetRange(3.,5);
//               fitC->SetRange(3.5,5);
               fitC->Draw("same");
	       *C    = fitC->GetParameter(0);
	       *Cerr = fitC->GetParError(0);
	       cout<< "FitResultPion : "<<*C<<"   " <<*Cerr<< endl;
               
               inputnewPion->Draw("COLZ");
               FitResultPion->Draw("same");
               fitC->Draw("same");
	       sprintf(buffer,"%sFitPion","fit/");
	       SaveCanvas(c1,"dirtest/",buffer);              



               char fitfunc[2048];
//   si on fixe C
//               *C    = 2.93; 
               sprintf(fitfunc,"[0]*pow(%6.4f/x,2) +%4.3f",MassCenter, *C);
               cout << "fitfunc " << fitfunc << endl;
	       TF1* myfit = new TF1("myfit",fitfunc, MinRange, MaxRange);
	       // TF1* myfit = new TF1("myfit","[0]*pow(1.8756/x,2) + [1]", MinRange, MaxRange); //1875.6 MeV  deuteron mass
	       myfit->SetParName  (0,"K");
	       //myfit->SetParName  (1,"C");
//	       myfit->SetParameter(0, 1.8);
	       myfit->SetParameter(0, 2.);
	       //myfit->SetParameter(1, *C);
//	       myfit->SetParLimits(0, 1.3,4.0); //
	       myfit->SetParLimits(0, 1.0,4.0); //
	       // myfit->SetParLimits(1, 3.8,3.8); // 
	       myfit->SetLineWidth(2);
	       myfit->SetLineColor(2);
	       FitResult->Fit("myfit", "M R E I 0");
	       myfit->SetRange(MinRange,MaxRange);
               inputnew->Draw("COLZ");
               FitResult->Draw("same");
	       myfit->Draw("same");

	       //double prevConstants [] = {*K, *Kerr, *C, *Cerr};
	       *K    = myfit->GetParameter(0);
               //*C    = myfit->GetParameter(1);
               *Kerr = myfit->GetParError(0);
               //*Cerr = myfit->GetParError(1);
	       cout<< "FitResult : "<<*K<<"   " <<*Kerr<< endl;
	       //cout<< "FitResult : "<<myfit->GetParameter(1)<<"   " <<myfit->GetParError(1)<< endl;
               printf("K Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n",
                prevConstants[0], prevConstants[1], *K, *Kerr, 100.0*((*K)-prevConstants[0])/(*K));
	       printf("C Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n",
                prevConstants[2], prevConstants[3], *C, *Cerr, 100.0*((*C)-prevConstants[2])/(*C));

          if(std::max(fabs(100.0*((*K)-prevConstants[0])/(*K)), fabs(100.0*((*C)-prevConstants[2])/(*C)))<1.0)
             hasConverged=true;  //<1% variation of the constant --> converged

	       TPaveText* st = new TPaveText(0.40,0.78,0.79,0.89, "NDC");
	       st->SetFillColor(0);
	       sprintf(buffer,"K = %4.3f +- %6.4f",myfit->GetParameter(0), myfit->GetParError(0));
	       st->AddText(buffer);
	       //sprintf(buffer,"C = %4.3f +- %6.4f",myfit->GetParameter(1), myfit->GetParError(1));
	       sprintf(buffer,"C = %4.3f +- %6.4f",*C, *Cerr);
	       st->AddText(buffer);
	       st->Draw("same");
	       sprintf(buffer,"%sFit","fit/");
	       SaveCanvas(c1,"dirtest/",buffer);              
	       delete c1;

          delete line1;
          delete line2;
          delete myfit;
          delete FitResult;
          delete inputnew;
       }
}

void ExtracttheCConstants(TH2D* inputlong, double* C,  double* Cerr) 
{
       char buffer[2048];

	      TH2D* inputnewPion = (TH2D*)inputlong->Clone("tempTH2D2");
	      for(int x=0;x<=inputnewPion->GetNbinsX()+1;x++){
   	        for(int y=0;y<=inputnewPion->GetNbinsY()+1;y++){
                    if (inputnewPion->GetXaxis()->GetBinCenter(x)<10) {
                       inputnewPion->SetBinContent(x,y,0); 
                       inputnewPion->SetBinError(x,y,0); 
                    }
	  	    if (inputnewPion->GetYaxis()->GetBinCenter(y)<2 || inputnewPion->GetYaxis()->GetBinCenter(y)>4.2) {
                       inputnewPion->SetBinContent(x,y,0);
                       inputnewPion->SetBinError(x,y,0);
                    }
	      }}


	      TCanvas* c1 = new TCanvas("c1", "c1", 800,600);
	      c1->SetLogz(true);

	      inputnewPion->SetStats(kFALSE);
	      inputnewPion->GetXaxis()->SetTitle("track momentum (GeV)");
	      inputnewPion->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
	      inputnewPion->SetAxisRange(0,25,"X");
	      inputnewPion->SetAxisRange(0,15,"Y");
              inputnewPion->Draw("COLZ");

	      SaveCanvas(c1, "dirtest/", "testCalone_region");


	      TH1D* ProjectionPion = (TH1D*)inputnewPion->ProjectionY("proj");
              cout << " mean value " << ProjectionPion->GetMean() << " RMS " << ProjectionPion->GetRMS() << endl;
              float xminval = ProjectionPion->GetMean() - 0.5;
              float xmaxval = ProjectionPion->GetMean() + 0.5;
//              TF1* mygausPion = new TF1("mygausPion","gaus", 2., 5.);
              TF1* mygausPion = new TF1("mygausPion","gaus", xminval, xmaxval);
              ProjectionPion->Fit("mygausPion","Q0 RME");
	      cout<< " mygausPion  mean   " <<mygausPion->GetParameter(1) << "  "<<mygausPion->GetParError(1)<<endl;
	      cout<< " mygausPion  sigma   " <<mygausPion->GetParameter(2)<<"  "<<mygausPion->GetParError(2)<<endl;

	       *C    = mygausPion->GetParameter(1);
	       *Cerr = mygausPion->GetParameter(2);
	       cout<< "FitResultPion : "<<*C<<"   " <<*Cerr<< endl;
               
               ProjectionPion->Draw();
               mygausPion->Draw("same");
	       sprintf(buffer,"%sFitC_highp","fit/");
	       SaveCanvas(c1,"dirtest/",buffer);              

}

void ExtractConstants_sig(TH2D* input, double* K, double* C, double* Kerr, double* Cerr,
      double MinRange, double MaxRange, double MassCenter)
{
       char buffer[2048];
       bool hasConverged = false;

       for(unsigned int loop=0;loop<5 and !hasConverged; loop++){
	      TH2D* inputnew = (TH2D*)input->Clone("tempTH2D");
	      for(int x=1;x<=inputnew->GetNbinsX();x++){
	       for(int y=1;y<=inputnew->GetNbinsY();y++){
                if (inputnew->GetYaxis()->GetBinCenter(y)<5. && inputnew->GetXaxis()->GetBinCenter(x)<1000) {
		  inputnew->SetBinContent(x,y,0);        
		}
                else if (inputnew->GetYaxis()->GetBinCenter(y)<2.) {
		  inputnew->SetBinContent(x,y,0);        
		}
	       }
              }  

	      

	      TCanvas* c1 = new TCanvas("c1", "c1", 800,600);
	      c1->SetLogz(true);
	      inputnew->SetStats(kFALSE);
	      inputnew->GetXaxis()->SetTitle("track momentum (GeV)");
	      inputnew->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
	      inputnew->Draw("COLZ");

	      SaveCanvas(c1, "dirtest/", "Sig_dedxVsP");

	      delete c1;
               TH1D* FitResult = new TH1D("FitResult"       , "FitResult"      ,inputnew->GetXaxis()->GetNbins(),inputnew->GetXaxis()->GetXmin(),inputnew->GetXaxis()->GetXmax());
	       FitResult->SetTitle("");
	       FitResult->SetStats(kFALSE);  
	       FitResult->GetXaxis()->SetTitle("P [GeV]");
	       FitResult->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
	       FitResult->GetYaxis()->SetTitleOffset(1.20);
	       FitResult->Reset();     

	       for(int x=inputnew->GetXaxis()->FindBin(MinRange);x<inputnew->GetXaxis()->FindBin(MaxRange);x++){
		  double P       = inputnew->GetXaxis()->GetBinCenter(x);
	    
		  TH1D* Projection = (TH1D*)(inputnew->ProjectionY("proj",x,x))->Clone();
		  if(Projection->Integral()<100)continue;
		  Projection->SetAxisRange(0.1,25,"X");
		  Projection->Sumw2();
		  Projection->Scale(1.0/Projection->Integral());


		  TF1* mygaus = new TF1("mygaus","gaus", 2.5, 20);
		  Projection->Fit("mygaus","Q0 RME");
		  double chiFromFit  = (mygaus->GetChisquare())/(mygaus->GetNDF());
		  FitResult->SetBinContent(x, mygaus->GetParameter(1));
		  FitResult->SetBinError  (x, mygaus->GetParError (1));
		  mygaus->SetLineColor(2);
		  mygaus->SetLineWidth(2);

		  c1  = new TCanvas("canvas", "canvas", 800,600);
		  Projection->Draw();
		  Projection->SetTitle("");
		  Projection->SetStats(kFALSE);
		  Projection->GetXaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
		  Projection->GetYaxis()->SetTitle("#Entries");
		  Projection->GetYaxis()->SetTitleOffset(1.30);
		  Projection->SetAxisRange(1E-5,1.0,"Y");

		  mygaus->Draw("same");


		  TPaveText* stt = new TPaveText(0.55,0.82,0.79,0.92, "NDC");
		  stt->SetFillColor(0);
		  stt->SetTextAlign(31);
		  sprintf(buffer,"Proton  #mu:%5.1fMeV/cm",mygaus->GetParameter(1));      stt->AddText(buffer);
		  sprintf(buffer,"Proton  #sigma:%5.1fMeV/cm",mygaus->GetParameter(2));      stt->AddText(buffer);
		  stt->Draw("same");

		  c1->SetLogy(true);
		  sprintf(buffer,"%sProjectionFit_P%03i_%03i","fit/",(int)(100*FitResult->GetXaxis()->GetBinLowEdge(x)),(int)(100*FitResult->GetXaxis()->GetBinUpEdge(x)) );
		  if(P>=MinRange && P<=MaxRange){SaveCanvas(c1,"dirtest/",buffer);}
		  delete c1;
                  delete Projection;
                  delete mygaus;
                  delete stt;
	       }

	       c1  = new TCanvas("canvas", "canvas", 800,600);
	       FitResult->Draw("");

	       TLine* line1 = new TLine(MinRange, FitResult->GetMinimum(), MinRange, FitResult->GetMaximum());
	       line1->SetLineWidth(2);
	       line1->SetLineStyle(2);
	       line1->Draw();

	       TLine* line2 = new TLine(MaxRange, FitResult->GetMinimum(), MaxRange, FitResult->GetMaximum());
	       line2->SetLineWidth(2);
	       line2->SetLineStyle(2);
	       line2->Draw();


	       double prevConstants [] = {*K, *Kerr, *C, *Cerr};

               char fitfunc[2048];
               sprintf(fitfunc,"[0]*pow(%6.4f/x,2) +[1]",MassCenter);
               cout << "fitfunc " << fitfunc << endl;
	       TF1* myfit = new TF1("myfit",fitfunc, MinRange, MaxRange);
	       myfit->SetParName  (0,"K");
	       myfit->SetParName  (1,"C");
	       myfit->SetParameter(0, 1.8);
	       myfit->SetParameter(1, 4.0);
	       myfit->SetParLimits(0, 1.3,4.0); //
	       myfit->SetParLimits(1, 3.0,7.0); // 
	       myfit->SetLineWidth(2);
	       myfit->SetLineColor(2);
	       FitResult->Fit("myfit", "M R E I 0");
	       myfit->SetRange(MinRange,MaxRange);
               inputnew->Draw("COLZ");
               FitResult->Draw("same");
	       myfit->Draw("same");

	       //double prevConstants [] = {*K, *Kerr, *C, *Cerr};
	       *K    = myfit->GetParameter(0);
               *C    = myfit->GetParameter(1);
               *Kerr = myfit->GetParError(0);
               *Cerr = myfit->GetParError(1);
	       cout<< "FitResult : "<<*K<<"   " <<*Kerr<< endl;
	       cout<< "FitResult : "<<*C<<"   " <<*Cerr<< endl;
               printf("K Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n",
                prevConstants[0], prevConstants[1], *K, *Kerr, 100.0*((*K)-prevConstants[0])/(*K));
	       printf("C Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n",
                prevConstants[2], prevConstants[3], *C, *Cerr, 100.0*((*C)-prevConstants[2])/(*C));

          if(std::max(fabs(100.0*((*K)-prevConstants[0])/(*K)), fabs(100.0*((*C)-prevConstants[2])/(*C)))<1.0)
             hasConverged=true;  //<1% variation of the constant --> converged

	       TPaveText* st = new TPaveText(0.40,0.78,0.79,0.89, "NDC");
	       st->SetFillColor(0);
	       sprintf(buffer,"K = %4.3f +- %6.4f",myfit->GetParameter(0), myfit->GetParError(0));
	       st->AddText(buffer);
	       //sprintf(buffer,"C = %4.3f +- %6.4f",myfit->GetParameter(1), myfit->GetParError(1));
	       sprintf(buffer,"C = %4.3f +- %6.4f",*C, *Cerr);
	       st->AddText(buffer);
	       st->Draw("same");
	       sprintf(buffer,"%sFit","fit/");
	       SaveCanvas(c1,"dirtest/",buffer);              
	       delete c1;

          delete line1;
          delete line2;
          delete myfit;
          delete FitResult;
          delete inputnew;
       }
}



double GetMass(double P, double I, double* K, double* C){
  //  *C=4.591; // to remove
  //  *K = 2.002; //to remove
         
      return sqrt((I-(*C))/(*K))*P;
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


TF1* GetMassLine(double M, double K, double C, bool left=false)
{  
   double BetaMax = 0.99;
   double PMax = sqrt((BetaMax*BetaMax*M*M)/(1-BetaMax*BetaMax));

   double BetaMin = 0.1;
   double PMin = sqrt((BetaMin*BetaMin*M*M)/(1-BetaMin*BetaMin));
   
   if(left){PMax*=-1; PMin*=-1;}

   TF1* MassLine = new TF1("MassLine","[2] + ([0]*[0]*[1])/(x*x)", PMin, PMax);
   MassLine->SetParName  (0,"M");
   MassLine->SetParName  (1,"K");
   MassLine->SetParName  (2,"C");
   MassLine->SetParameter(0, M);
   MassLine->SetParameter(1, K);
   MassLine->SetParameter(2, C);
   MassLine->SetLineWidth(2);
   return MassLine;
}
void compaFits(TString filename1, TString filename2, double Kval, double Cval){

   
   TFile* myfile1;
   myfile1 = new TFile (filename1);
   TFile* myfile2;
   myfile2 = new TFile (filename2);

   TH2D* HdedxVsP_1fit1;
   TH2D* HdedxVsP_2fit1;
   TH2D* HdedxVsP_1fit2;
   TH2D* HdedxVsP_2fit2;
   TH1D* resuHdedxVsP_1fit1;
   TH1D* resuHdedxVsP_2fit1 ;
   TH1D* resuHdedxVsP_1fit2;
   TH1D* resuHdedxVsP_2fit2 ;


   myfile1->cd();
//    // Ih no drop, no L1
   HdedxVsP_1fit1 = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp");
   HdedxVsP_2fit1 = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp2");
   HdedxVsP_2fit1->Rebin2D(2,1);

   resuHdedxVsP_1fit1 = new TH1D("resuHdedxVsP_1fit1", "resuHdedxVsP_1fit1", HdedxVsP_1fit1->GetXaxis()->GetNbins(),HdedxVsP_1fit1->GetXaxis()->GetXmin(),HdedxVsP_1fit1->GetXaxis()->GetXmax());
	       resuHdedxVsP_1fit1->SetTitle("");
	       resuHdedxVsP_1fit1->SetStats(kFALSE);  
	       resuHdedxVsP_1fit1->GetXaxis()->SetTitle("P [GeV]");
	       resuHdedxVsP_1fit1->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
	       resuHdedxVsP_1fit1->GetYaxis()->SetTitleOffset(1.20);
	       resuHdedxVsP_1fit1->Reset();     

  resuHdedxVsP_2fit1 = new TH1D("resuHdedxVsP_2fit1", "resuHdedxVsP_2fit1" ,HdedxVsP_2fit1->GetXaxis()->GetNbins(),HdedxVsP_2fit1->GetXaxis()->GetXmin(),HdedxVsP_2fit1->GetXaxis()->GetXmax());
               resuHdedxVsP_2fit1->SetTitle("");
               resuHdedxVsP_2fit1->SetStats(kFALSE);
               resuHdedxVsP_2fit1->GetXaxis()->SetTitle("P [GeV]");
               resuHdedxVsP_2fit1->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
               resuHdedxVsP_2fit1->GetYaxis()->SetTitleOffset(1.20);
               resuHdedxVsP_2fit1->Reset();


   ExtractProj(HdedxVsP_2fit1, HdedxVsP_1fit1, &Kval, &Cval, 0.5, 1.5, 0.938, 0.2, 0.5, 4.2,resuHdedxVsP_2fit1,resuHdedxVsP_1fit1);

   myfile2->cd();
   HdedxVsP_1fit2 = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp");
   HdedxVsP_2fit2 = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp2");
   HdedxVsP_2fit2->Rebin2D(2,1);
   resuHdedxVsP_1fit2 = new TH1D("resuHdedxVsP_1fit2"       , "resuHdedxVsP_1fit2"      ,HdedxVsP_1fit2->GetXaxis()->GetNbins(),HdedxVsP_1fit2->GetXaxis()->GetXmin(),HdedxVsP_1fit2->GetXaxis()->GetXmax());
	       resuHdedxVsP_1fit2->SetTitle("");
	       resuHdedxVsP_1fit2->SetStats(kFALSE);  
	       resuHdedxVsP_1fit2->GetXaxis()->SetTitle("P [GeV]");
	       resuHdedxVsP_1fit2->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
	       resuHdedxVsP_1fit2->GetYaxis()->SetTitleOffset(1.20);
	       resuHdedxVsP_1fit2->Reset();     

   resuHdedxVsP_2fit2 = new TH1D("resuHdedxVsP_2fit2", "resuHdedxVsP_2fit2" ,HdedxVsP_2fit2->GetXaxis()->GetNbins(),HdedxVsP_2fit2->GetXaxis()->GetXmin(),HdedxVsP_2fit2->GetXaxis()->GetXmax());
               resuHdedxVsP_2fit2->SetTitle("");
               resuHdedxVsP_2fit2->SetStats(kFALSE);
               resuHdedxVsP_2fit2->GetXaxis()->SetTitle("P [GeV]");
               resuHdedxVsP_2fit2->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
               resuHdedxVsP_2fit2->GetYaxis()->SetTitleOffset(1.20);
               resuHdedxVsP_2fit2->Reset();

   ExtractProj(HdedxVsP_2fit2, HdedxVsP_1fit2, &Kval, &Cval, 0.5, 1.5, 0.938, 0.2, 0.5, 4.2,resuHdedxVsP_2fit2,resuHdedxVsP_1fit2);


   TLegend* leg = new TLegend (0.25, 0.70, 0.65, 0.90);
//   leg->SetHeader ("Fitting the Profile");
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   leg->AddEntry (resuHdedxVsP_1fit1, filename1, "LP");
   leg->AddEntry (resuHdedxVsP_1fit2, filename2, "LP");


   TCanvas* c2 = new TCanvas("c2", "c2", 800,600);
   c2->cd();
   TF1* ProtonLine = GetMassLine(0.938,Kval, Cval, false);
   ProtonLine->SetLineColor(1);
   ProtonLine->SetLineWidth(2);

   TLine* linepion = new TLine(2.,Cval,25., Cval);
   linepion->SetLineColor(6);
   linepion->SetLineWidth(2);

   resuHdedxVsP_1fit1->SetLineColor(2);
   resuHdedxVsP_1fit1->SetMarkerColor(2);
   resuHdedxVsP_2fit1->SetLineColor(2);
   resuHdedxVsP_2fit1->SetMarkerColor(2);
   resuHdedxVsP_1fit2->SetLineColor(4);
   resuHdedxVsP_1fit2->SetMarkerColor(4);
   resuHdedxVsP_2fit2->SetLineColor(4);
   resuHdedxVsP_2fit2->SetMarkerColor(4);

   resuHdedxVsP_1fit1->SetMaximum(13.);
   resuHdedxVsP_1fit1->Draw();
   resuHdedxVsP_1fit2->Draw("same");
   ProtonLine->Draw("same");
   linepion->Draw("same");
   leg->Draw("same");
   SaveCanvas(c2, "dirtest/", "compa_lines");
   

  TCanvas* c3 = new TCanvas("c3", "c3", 800,600);
   c3->cd();
   resuHdedxVsP_2fit1->SetMinimum(2.5);
   resuHdedxVsP_2fit1->SetMaximum(5.);
   resuHdedxVsP_2fit1->Draw();
   resuHdedxVsP_2fit2->Draw("same");
//   ProtonLine->Draw("same");
   linepion->Draw("same");
   leg->Draw("same");
   SaveCanvas(c3, "dirtest/", "compa_long");
   
  gStyle->SetOptStat(0);
  c2->cd();
  c2->SetLogz(1);
  SaveCanvas(c2, "dirtest/", "compa_lines_logZ"); 
   c3->cd();
  c3->SetLogz(1);
   SaveCanvas(c3, "dirtest/", "compa_long_logZ");


   

}

void ExtractProj(TH2D* inputlong, TH2D* inputshort,double* K, double* C, 
      double MinRange, double MaxRange, double MassCenter, double LeftMassMargin, double RightMassMargin, double yPionMax,
      TH1D* FitResultPion,  TH1D* FitResult)
{
	      TH2D* inputnew = (TH2D*)inputshort->Clone("tempTH2D");
	      TH2D* inputnewPion = (TH2D*)inputlong->Clone("tempTH2D2");
	      for(int x=0;x<=inputnew->GetNbinsX()+1;x++){
	      for(int y=0;y<=inputnew->GetNbinsY()+1;y++){
		double Mass = GetMass(inputnew->GetXaxis()->GetBinCenter(x),inputnew->GetYaxis()->GetBinCenter(y), K, C);
		if(isnan (float(Mass)) || Mass<MassCenter-(LeftMassMargin) || Mass>MassCenter+RightMassMargin){
		  inputnew->SetBinContent(x,y,0);        
		  inputnew->SetBinError(x,y,0);        
		}
	      }}
	      for(int x=0;x<=inputnewPion->GetNbinsX()+1;x++){
   	        for(int y=0;y<=inputnewPion->GetNbinsY()+1;y++){
                    if (inputnewPion->GetXaxis()->GetBinCenter(x)<2.5) {
                          inputnewPion->SetBinContent(x,y,0);
                          inputnewPion->SetBinError(x,y,0);
                    }
	  	    if (inputnewPion->GetYaxis()->GetBinCenter(y)<2 || inputnewPion->GetYaxis()->GetBinCenter(y)>yPionMax) {
                          inputnewPion->SetBinContent(x,y,0);
                          inputnewPion->SetBinError(x,y,0);
                    }
	      }}

	      

               for(int x=1;x<inputnewPion->GetNbinsX();x++){
		  TH1D* ProjectionPion = (TH1D*)(inputnewPion->ProjectionY("proj",x,x))->Clone();
                  if(ProjectionPion->Integral()<100)continue;
                  ProjectionPion->SetAxisRange(0.1,25,"X");
//                  ProjectionPion->Sumw2();
                  ProjectionPion->Scale(1.0/ProjectionPion->Integral());

                  TF1* mygausPion = new TF1("mygausPion","gaus", 2., 15);
                  ProjectionPion->Fit("mygausPion","Q0 RME");
                  FitResultPion->SetBinContent(x, mygausPion->GetParameter(1));
//		  cout<< " mygausPion     " <<x<<"  "<<mygausPion->GetParameter(1)<<endl;
		  FitResultPion->SetBinError  (x, mygausPion->GetParError (1));
                  delete ProjectionPion;
                  delete mygausPion;
               }

	       for(int x=1;x<inputnew->GetXaxis()->FindBin(1.5);x++){
		  double P       = inputnew->GetXaxis()->GetBinCenter(x);
	    
		  TH1D* Projection = (TH1D*)(inputnew->ProjectionY("proj",x,x))->Clone();
		  if(Projection->Integral()<100)continue;
		  Projection->SetAxisRange(0.1,25,"X");
//		  Projection->Sumw2();
		  Projection->Scale(1.0/Projection->Integral());

		  TF1* mygaus = new TF1("mygaus","gaus", 2.5, 15);
		  Projection->Fit("mygaus","Q0 RME");
		  double chiFromFit  = (mygaus->GetChisquare())/(mygaus->GetNDF());
		  FitResult->SetBinContent(x, mygaus->GetParameter(1));
		  FitResult->SetBinError  (x, mygaus->GetParError (1));
		  mygaus->SetLineColor(2);
		  mygaus->SetLineWidth(2);

                  delete Projection;
                  delete mygaus;
	       }

          delete inputnewPion;
          delete inputnew;
}
void TestFitKandC(TString filename){
   double K_ = 3.0; double Kerr_ = 0.001; double C_ = 3.; double Cerr_ = 0.001;
   double Ktmp = K_ ; double Ctmp = C_; double KerrTmp = Kerr_; double CerrTmp = Cerr_;
   double Ntmp = 0.3 ;  double NerrTmp = 0.001;
   TFile* myfile;
   myfile = new TFile (filename);

   TH2D* HdedxVsP_1fit;
   TH2D* HdedxVsP_2fit;
   myfile->cd();
//    // Ih no drop, no L1
   HdedxVsP_1fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp");
   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("dEdX0noL1VsP_lowp2");
   HdedxVsP_2fit->Rebin2D(2,1);
   TestExtractConstants (HdedxVsP_2fit, HdedxVsP_1fit, &Ktmp, &Ctmp, &Ntmp, &KerrTmp, &CerrTmp, &NerrTmp, 0.5, 20., 0.938, 0.2, 0.5, 4.2);
   

   TCanvas* c2 = new TCanvas("c2", "c2", 800,600);
   c2->cd();
   gStyle->SetOptStat(0);
   TF1* PionLine = TestGetMassLine(0.140, Ktmp, Ctmp, Ntmp, false);
   PionLine->SetLineColor(1);
   PionLine->SetLineWidth(2);
   TF1* KaonLine = TestGetMassLine(0.494, Ktmp, Ctmp, Ntmp, false);
   KaonLine->SetLineColor(1);
   KaonLine->SetLineWidth(2);
   TF1* ProtonLine = TestGetMassLine(0.938,Ktmp, Ctmp, Ntmp, false);
   ProtonLine->SetLineColor(1);
   ProtonLine->SetLineWidth(2);
   TF1* DeuteronLine = TestGetMassLine(1.8756,Ktmp, Ctmp, Ntmp, false);
   DeuteronLine->SetLineColor(1);
   DeuteronLine->SetLineWidth(2);

   HdedxVsP_2fit->Draw("colz");
   PionLine->Draw("same");
   KaonLine->Draw("same");
   ProtonLine->Draw("same");
   if (filename=="analysis_ul_2017_21jan.root") DeuteronLine->Draw("same");
   SaveCanvas(c2, "dirtest/", "TestdedxVsP_long");
   
   c2->SetLogz(1);
    SaveCanvas(c2, "dirtest/", "TestdedxVsP_long_logZ"); 
   
   TCanvas* c3 = new TCanvas("c3", "c3", 800,600);
   c3->cd();
   gStyle->SetOptStat(0);
   HdedxVsP_1fit->Draw("colz");
   PionLine->Draw("same");
   KaonLine->Draw("same");
   ProtonLine->Draw("same");
   if (filename=="analysis_ul_2017_21jan.root") DeuteronLine->Draw("same");
   SaveCanvas(c3, "dirtest/", "TestdedxVsP_lines");
   c3->SetLogz(1);
    SaveCanvas(c3, "dirtest/", "TestdedxVsP_lines_logZ"); 

}
void TestExtractConstants(TH2D* inputlong, TH2D* inputshort, double* K, double* C, double* N, double* Kerr, double* Cerr, double* Nerr, double MinRange , double MaxRange , double MassCenter , double LeftMassMargin , double RightMassMargin , double yPionMax) {
       char buffer[2048];
       bool hasConverged = false;

       for(unsigned int loop=0;loop<8 and !hasConverged; loop++){
              TH2D* Myinputnew = new TH2D("Myinputnew","Myinputnew", 250,0,25, 80, 2.,10.);
              Myinputnew->Sumw2();
	      TH2D* inputnew = (TH2D*)inputshort->Clone("tempTH2D");
	      TH2D* inputnewPion = (TH2D*)inputlong->Clone("tempTH2D2");
	      for(int x=0;x<=inputnew->GetNbinsX()+1;x++){
	      for(int y=0;y<=inputnew->GetNbinsY()+1;y++){
		double Mass = GetMass(inputnew->GetXaxis()->GetBinCenter(x),inputnew->GetYaxis()->GetBinCenter(y), K, C);
		if(isnan (float(Mass)) || Mass<MassCenter-(LeftMassMargin) || Mass>MassCenter+RightMassMargin){
		  inputnew->SetBinContent(x,y,0);        
		  inputnew->SetBinError(x,y,0);        
		}
                else {
                  if (Myinputnew->GetXaxis()->GetBinCenter(x)<1.5) {
                         Myinputnew->SetBinContent(x,y,inputnew->GetBinContent(x,y)); 
                         Myinputnew->SetBinError(x,y,inputnew->GetBinError(x,y)); 
                  }
             
                }
	      }}
	      for(int x=0;x<=inputnewPion->GetNbinsX()+1;x++){
   	        for(int y=0;y<=inputnewPion->GetNbinsY()+1;y++){
                    if (inputnewPion->GetXaxis()->GetBinCenter(x)<3.0) {
                        inputnewPion->SetBinContent(x,y,0);
                        inputnewPion->SetBinError(x,y,0);
                    }
	  	    if (inputnewPion->GetYaxis()->GetBinCenter(y)<2 || inputnewPion->GetYaxis()->GetBinCenter(y)>yPionMax) {
                       inputnewPion->SetBinContent(x,y,0);
                       inputnewPion->SetBinError(x,y,0);
                    }
                    Myinputnew->Fill(inputnewPion->GetXaxis()->GetBinCenter(x),inputnewPion->GetYaxis()->GetBinCenter(y),inputnewPion->GetBinContent(x,y));
	      }}

	      

	      TCanvas* c1 = new TCanvas("c1", "c1", 800,600);
	      c1->SetLogz(true);
	      Myinputnew->SetStats(kFALSE);
	      Myinputnew->GetXaxis()->SetTitle("track momentum (GeV)");
	      Myinputnew->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
	      Myinputnew->SetAxisRange(0,25,"X");
	      Myinputnew->SetAxisRange(0,15,"Y");
	      Myinputnew->Draw("COLZ");

	//      KaonLine->Draw("same");
	//      ProtonLine->Draw("same");
	//      DeuteronLine->Draw("same");
	//      TritonLine->Draw("same");

	      SaveCanvas(c1, "dirtest/", "testMyinputnew");

	      delete c1;
               TH1D* FitResult = new TH1D("FitResult"       , "FitResult"      ,Myinputnew->GetXaxis()->GetNbins(),Myinputnew->GetXaxis()->GetXmin(),Myinputnew->GetXaxis()->GetXmax());
	       FitResult->SetTitle("");
	       FitResult->SetStats(kFALSE);  
	       FitResult->GetXaxis()->SetTitle("P [GeV]");
	       FitResult->GetYaxis()->SetTitle("dE/dx Estimator [MeV/cm]");
	       FitResult->GetYaxis()->SetTitleOffset(1.20);
	       FitResult->Reset();     

               for(int x=1;x<Myinputnew->GetNbinsX();x++){
		  TH1D* ProjectionPion = (TH1D*)(Myinputnew->ProjectionY("proj",x,x))->Clone();
                  if(ProjectionPion->Integral()<100)continue;
                  ProjectionPion->SetAxisRange(0.1,25,"X");
                  ProjectionPion->Scale(1.0/ProjectionPion->Integral());

                  TF1* mygausPion = new TF1("mygausPion","gaus", 2., 15);
                  ProjectionPion->Fit("mygausPion","Q0 RME");
                  FitResult->SetBinContent(x, mygausPion->GetParameter(1));
		  FitResult->SetBinError  (x, mygausPion->GetParError (1));
               }

	       c1  = new TCanvas("canvas", "canvas", 800,600);
	       FitResult->SetAxisRange(0,25,"X");
	       FitResult->SetAxisRange(0,15,"Y");
	       FitResult->Draw("");
	      SaveCanvas(c1, "dirtest/", "testFitResult");

	       double prevConstants [] = {*K, *Kerr, *C, *Cerr, *N, *Nerr};


               char fitfunc1[2048];
               sprintf(fitfunc1,"[0]*pow(%6.4f/x,2) +[1] + [2]*log(x/%6.4f)",MassCenter,MassCenter);
//	       TF1* myfit = new TF1("myfit",fitfunc1, MinRange, MaxRange);
	       TF1* myfit = new TF1("myfit",TestfitSpecial, MinRange, MaxRange,3);
//	       TF1* myfit = new TF1("myfit",TestfitSpecial2, MinRange, MaxRange,3);
	       myfit->SetParName  (0,"K");
	       myfit->SetParName  (1,"C");
	       myfit->SetParName  (2,"N");
	       myfit->SetParameter(0, 2.);
	       myfit->SetParameter(1, 3.3);
	       myfit->SetParameter(2, 0.2);
	       myfit->SetParLimits(0, 1.0,4.0); //
	       myfit->SetParLimits(1, 2.,5); // 
	       myfit->SetParLimits(2, 0.0001,1000); // 

	       myfit->SetLineWidth(2);
	       myfit->SetLineColor(2);
	       FitResult->Fit("myfit", "M R E I 0");
	       myfit->SetRange(MinRange,MaxRange);
               Myinputnew->Draw("COLZ");
               FitResult->SetMarkerSize(21);
               FitResult->SetLineWidth(2);
               FitResult->Draw("same");
	       myfit->Draw("same");
               FitResult->Draw("same");
               c1->SetLogz(1);

	       //double prevConstants [] = {*K, *Kerr, *C, *Cerr};
	       *K    = myfit->GetParameter(0);
               *C    = myfit->GetParameter(1);
               *N    = myfit->GetParameter(2);
               *Kerr = myfit->GetParError(0);
               *Cerr = myfit->GetParError(1);
               *Nerr = myfit->GetParError(2);
	       cout<< "FitResult : "<<*K<<"   " <<*Kerr<< "   " << *C <<"   " <<*Cerr<< "   " << *N<<"   " <<*Nerr<<endl;
	       //cout<< "FitResult : "<<myfit->GetParameter(1)<<"   " <<myfit->GetParError(1)<< endl;
               printf("K Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n",
                prevConstants[0], prevConstants[1], *K, *Kerr, 100.0*((*K)-prevConstants[0])/(*K));
	       printf("C Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n",
                prevConstants[2], prevConstants[3], *C, *Cerr, 100.0*((*C)-prevConstants[2])/(*C));
	       if (*N>0) printf("N Constant changed from %6.4f+-%6.4f to %6.4f+-%6.4f    (diff = %6.3f%%)\n",
                prevConstants[4], prevConstants[5], *N, *Nerr, 100.0*((*N)-prevConstants[4])/(*N));

          if(std::max(fabs(100.0*((*K)-prevConstants[0])/(*K)), fabs(100.0*((*C)-prevConstants[2])/(*C)))<1.0) {
             if (*N>0) {
                if ( fabs(100.0*((*N)-prevConstants[4])/(*N))<1.0) hasConverged=true; 
             }
             else {
               hasConverged=true;  //<1% variation of the constant --> converged
             }
          }

	       TPaveText* st = new TPaveText(0.40,0.78,0.79,0.89, "NDC");
	       st->SetFillColor(0);
	       sprintf(buffer,"K = %4.3f +- %6.4f",myfit->GetParameter(0), myfit->GetParError(0));
	       st->AddText(buffer);
	       //sprintf(buffer,"C = %4.3f +- %6.4f",myfit->GetParameter(1), myfit->GetParError(1));
	       sprintf(buffer,"C = %4.3f +- %6.4f",*C, *Cerr);
	       st->AddText(buffer);
	       sprintf(buffer,"N = %4.3f +- %6.4f",*N, *Nerr);
	       st->AddText(buffer);
	       st->Draw("same");
	       sprintf(buffer,"%sFit","fit/");
	       SaveCanvas(c1,"dirtest/",buffer);              
	       delete c1;

          delete myfit;
          delete FitResult;
          delete Myinputnew;
       }
}

Double_t TestfitSpecial(Double_t *x, Double_t *par) {
  if (x[0]<2.5) {
    return par[0]*pow(0.938/x[0],2) + par[1] + par[2]*log(x[0]/0.938);
  }
  else {
    return par[0]*pow(0.140/x[0],2) + par[1] + par[2]*log(x[0]/0.140);
  }
}

Double_t TestfitSpecial2(Double_t *x, Double_t *par) {
  if (x[0]<2.5) {
    return par[0]*pow(0.938/x[0],2)*log(par[2]*x[0]/0.938) + par[1] ;
  }
  else {
    return par[0]*pow(0.140/x[0],2)*log(par[2]*x[0]/0.140) + par[1] ;
  }
}

TF1* TestGetMassLine(double M, double K, double C, double N, bool left=false)
{  

   TF1* MassLine = new TF1("MassLine","[2] + ([0]*[0]*[1])/(x*x) +[3]*log(x/[0])", 0.5, 20.);
   MassLine->SetParName  (0,"M");
   MassLine->SetParName  (1,"K");
   MassLine->SetParName  (2,"C");
   MassLine->SetParName  (3,"N");
   MassLine->SetParameter(0, M);
   MassLine->SetParameter(1, K);
   MassLine->SetParameter(2, C);
   MassLine->SetParameter(3, N);
   MassLine->SetLineWidth(2);
   return MassLine;
}
