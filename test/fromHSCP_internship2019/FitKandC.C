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

void ExtractConstants(TH2D* input, double* K, double* C, double* Kerr, double* Cerr, double MinRange = 1.0, double MaxRange = 1.6, double MassCenter = 1.875, double LeftMassMargin = 0.2, double RightMassMargin = 0.8, double yPionMax=4.2); // by default use protons
double GetMass(double P, double I, double* K, double* C);
void SaveCanvas(TCanvas* c, std::string path, std::string name, bool OnlyPPNG=false);
TF1* GetMassLine(double M, double K, double C, bool left=false);
void FitKandC_sig();
void ExtractConstants_sig(TH2D* input, double* K, double* C, double* Kerr, double* Cerr, double MinRange = 1.0, double MaxRange = 1.6, double MassCenter = 1.875); // by default use protons

void FitKandC(){
   double K_ = 1.6; double Kerr_ = 0.001; double C_ = 4.3; double Cerr_ = 0.001;
   double Ktmp = K_ ; double Ctmp = C_; double KerrTmp = Kerr_; double CerrTmp = Cerr_;
   bool uncorr = false;
   TFile* myfile;
   if (uncorr) {
        myfile = new TFile ("minbias_template_uncorr_iter1.root");
        cout << " opening minbias_template_uncorr_iter1.root" << endl;
   }
   else {
        myfile = new TFile ("minbias_template_corr_iter1.root");
        cout << " opening minbias_template_corr_iter1.root" << endl;
   }
   
//   TFile* myfile = new TFile ("minbias_template_uncorr_iter1.root");
//   TFile* myfile = new TFile ("minbias_template_uncorr.root");
//   TFile* myfile = new TFile ("minbias_template_corr.root");
   TH2D* HdedxVsP_2fit;
   myfile->cd();
   HdedxVsP_2fit = (TH2D*)gROOT->FindObject("HdedxVsP");
//   ExtractConstants (HdedxVsP_2fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.7, 1.4, 0.938, 0.4, 0.8);
   if (uncorr) {
      ExtractConstants (HdedxVsP_2fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.5, 1.2, 0.938, 0.35, 0.6, 4.2);
   }
   else {
      ExtractConstants (HdedxVsP_2fit, &Ktmp, &Ctmp, &KerrTmp, &CerrTmp, 0.5, 1.2, 0.938, 0.35, 0.6, 4.6);
   }

   TCanvas* c2 = new TCanvas("c2", "c2", 600,600);
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

   HdedxVsP_2fit->Rebin2D(5,10);
   HdedxVsP_2fit->Draw("colz");
   PionLine->Draw("same");
   KaonLine->Draw("same");
   ProtonLine->Draw("same");
   SaveCanvas(c2, "dirtest/", "dedxVsP_lines");
   


   

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

   TCanvas* c2 = new TCanvas("c2", "c2", 600,600);
   c2->cd();
   TF1* SigLine = GetMassLine(2.4, Ktmp, Ctmp, false);
   SigLine->SetLineColor(1);
   SigLine->SetLineWidth(2);

   HdedxVsP_2fit->Draw("colz");
   SigLine->Draw("same");
   SaveCanvas(c2, "dirtest/", "SigdedxVsP_lines");
   


   

}


void ExtractConstants(TH2D* input, double* K, double* C, double* Kerr, double* Cerr,
      double MinRange, double MaxRange, double MassCenter, double LeftMassMargin, double RightMassMargin, double yPionMax)
{
       char buffer[2048];
       bool hasConverged = false;

       for(unsigned int loop=0;loop<5 and !hasConverged; loop++){
	      TH2D* inputnew = (TH2D*)input->Clone("tempTH2D");
	      TH2D* inputnewPion = (TH2D*)input->Clone("tempTH2D");
	      inputnew->Rebin2D(5,10);
	      inputnewPion->Rebin2D(5,10);
	      for(int x=1;x<=inputnew->GetNbinsX();x++){
	      for(int y=1;y<=inputnew->GetNbinsY();y++){
		double Mass = GetMass(inputnew->GetXaxis()->GetBinCenter(x),inputnew->GetYaxis()->GetBinCenter(y), K, C);
		if(isnan (float(Mass)) || Mass<MassCenter-(LeftMassMargin) || Mass>MassCenter+RightMassMargin){
		  inputnew->SetBinContent(x,y,0);        
		  //cout<<x<<"   "<<y<<endl;
		}
		if (inputnew->GetYaxis()->GetBinCenter(y)<2 || inputnew->GetYaxis()->GetBinCenter(y)>yPionMax) inputnewPion->SetBinContent(x,y,0);
//		if (inputnew->GetYaxis()->GetBinCenter(y)<2 || inputnew->GetYaxis()->GetBinCenter(y)>4.6) inputnewPion->SetBinContent(x,y,0);
		//cout<< inputnewPion->GetBinContent(x,y)<<endl;
	      }}

	      

	      TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
	      c1->SetLogz(true);
	      inputnew->SetStats(kFALSE);
	      inputnew->GetXaxis()->SetTitle("track momentum (GeV)");
	      inputnew->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
	      inputnew->SetAxisRange(0,5,"X");
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
	      inputnewPion->SetAxisRange(0,5,"X");
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

	       for(int x=1;x<inputnew->GetXaxis()->FindBin(5);x++){
		  double P       = inputnew->GetXaxis()->GetBinCenter(x);
	    
		  TH1D* Projection = (TH1D*)(inputnew->ProjectionY("proj",x,x))->Clone();
		  if(Projection->Integral()<100)continue;
		  Projection->SetAxisRange(0.1,25,"X");
		  Projection->Sumw2();
		  Projection->Scale(1.0/Projection->Integral());


		  TF1* mygaus = new TF1("mygaus","gaus", 2.5, 15);
		  Projection->Fit("mygaus","Q0 RME");
		  double chiFromFit  = (mygaus->GetChisquare())/(mygaus->GetNDF());
		  FitResult->SetBinContent(x, mygaus->GetParameter(1));
		  //cout<<x<<"  "<<mygaus->GetParameter(1)<<endl;
		  FitResult->SetBinError  (x, mygaus->GetParError (1));
		  mygaus->SetLineColor(2);
		  mygaus->SetLineWidth(2);

		  TH1D* ProjectionPion = (TH1D*)(inputnewPion->ProjectionY("proj",x,x))->Clone();
                  if(ProjectionPion->Integral()<100)continue;
                  ProjectionPion->SetAxisRange(0.1,25,"X");
                  ProjectionPion->Sumw2();
                  ProjectionPion->Scale(1.0/ProjectionPion->Integral());
		  

                  TF1* mygausPion = new TF1("mygausPion","gaus", 2.5, 15);
                  ProjectionPion->Fit("mygausPion","Q0 RME");
                  FitResultPion->SetBinContent(x, mygausPion->GetParameter(1));
		  cout<< " mygausPion     " <<x<<"  "<<mygausPion->GetParameter(1)<<endl;
		  FitResultPion->SetBinError  (x, mygausPion->GetParError (1));

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
                  delete Projection;
                  delete mygaus;
                  delete stt;
	       }

	       c1  = new TCanvas("canvas", "canvas", 600,600);
	       FitResult->SetAxisRange(0,2.5,"X");
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

	       TF1* fitC =  new TF1("fitC","[0]", 1,4);
	       fitC->SetParName(0,"C");
	       fitC->SetParameter(0, 3.);
	       fitC->SetParLimits(0, 2,5);
	       fitC->SetLineWidth(2);
               fitC->SetLineColor(2);
               FitResultPion->Fit("fitC", "M R E I 0");
//               fitC->SetRange(0.5,4);
               fitC->SetRange(1.,4);
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
               sprintf(fitfunc,"[0]*pow(%6.4f/x,2) +%4.3f",MassCenter, *C);
               cout << "fitfunc " << fitfunc << endl;
	       TF1* myfit = new TF1("myfit",fitfunc, MinRange, MaxRange);
	       // TF1* myfit = new TF1("myfit","[0]*pow(1.8756/x,2) + [1]", MinRange, MaxRange); //1875.6 MeV  deuteron mass
	       myfit->SetParName  (0,"K");
	       //myfit->SetParName  (1,"C");
	       myfit->SetParameter(0, 1.8);
	       //myfit->SetParameter(1, *C);
	       myfit->SetParLimits(0, 1.3,4.0); //
               //myfit->SetParLimits(1, *C,*C);
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

	      

	      TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
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

		  c1->SetLogy(true);
		  sprintf(buffer,"%sProjectionFit_P%03i_%03i","fit/",(int)(100*FitResult->GetXaxis()->GetBinLowEdge(x)),(int)(100*FitResult->GetXaxis()->GetBinUpEdge(x)) );
		  if(P>=MinRange && P<=MaxRange){SaveCanvas(c1,"dirtest/",buffer);}
		  delete c1;
                  delete Projection;
                  delete mygaus;
                  delete stt;
	       }

	       c1  = new TCanvas("canvas", "canvas", 600,600);
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
