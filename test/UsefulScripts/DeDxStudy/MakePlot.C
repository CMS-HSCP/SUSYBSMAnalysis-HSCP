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
#include "../../AnalysisCode/tdrstyle.C"
#include "../../AnalysisCode/Analysis_PlotFunction.h"
#include "FWCore/FWLite/interface/FWLiteEnabler.h"

using namespace std;

void ExtractConstants(TH2D* input, double* K, double* C, double* Kerr, double* Cerr, double MinRange = 1.0, double MaxRange = 1.6, double MassCenter = 1.875, double LeftMassMargin = 0.2, double RightMassMargin = 0.8); // by default use protons
double GetMass(double P, double I, double* K, double* C);

typedef struct dEdxPlotObj
{ // one such object per file
   TFile* InputFile;

   string FileName;
   string LegEntry;
   string SavePrefix;
   std::vector <string> StdObjName;
   std::vector <string> HitObjName;
   std::vector <string> StdObjLegend;
   std::vector <string> HitObjLegend;

   unsigned short type; // 0 - data, 1 - SMMC, 2 - signal MC
   map <string, double> K;   map <string, double> Kerr;
   map <string, double> C;   map <string, double> Cerr;
   
   // List of histograms to plot -- one for each *ObjName vector
   TH3F**       dEdxTemplate;
   TH1D**       hit_MIP;
   TH2D***      Charge_Vs_XYH;
   TH2D***      Charge_Vs_XYHN;
   TH2D***      Charge_Vs_XYLN;

   TH1D**       HdedxMIP;
   TH1D**       HMass;
   TH1D**       HMassHSCP;
   TH1D**       HProtonHitSO;
   TH1D**       HProtonHitPO;
   TH2D**       HdedxVsP;
   TH2D**       HdedxVsPSyst;
   TProfile**   HdedxVsPProfile;
   TH2D**       HdedxVsEta;
   TProfile**   HdedxVsEtaProfile;

   // constructor
   dEdxPlotObj (string FileName_, string LegEntry_, string SavePrefix_,
         vector<string> HitObjName_, vector<string> StdObjName_, vector<string> HitObjLegend_, vector<string> StdObjLegend_,
         unsigned short type_,
         double K_ = 2.26949, double Kerr_ = 0.001, double C_ = 3.32082, double Cerr_ = 0.001){

      // Initialize all the variables
      FileName    = FileName_;
      LegEntry    = LegEntry_;
      SavePrefix  = SavePrefix_;
      StdObjName  = StdObjName_;
      HitObjName  = HitObjName_;
      StdObjLegend= StdObjLegend_;
      HitObjLegend= HitObjLegend_;
      type        = type_;

      // open the root file
      InputFile = new TFile (FileName.c_str());

      // initialize all the histograms
      dEdxTemplate       = new TH3F*     [HitObjName.size()];
      hit_MIP            = new TH1D*     [HitObjName.size()];
      Charge_Vs_XYH      = new TH2D**    [HitObjName.size()];
      Charge_Vs_XYHN     = new TH2D**    [HitObjName.size()];
      Charge_Vs_XYLN     = new TH2D**    [HitObjName.size()];
      HdedxMIP           = new TH1D*     [StdObjName.size()];
      HMass              = new TH1D*     [StdObjName.size()];
      HMassHSCP          = new TH1D*     [StdObjName.size()];
      HProtonHitSO       = new TH1D*     [StdObjName.size()];
      HProtonHitPO       = new TH1D*     [StdObjName.size()];
      HdedxVsP           = new TH2D*     [StdObjName.size()];
      HdedxVsPSyst       = new TH2D*     [StdObjName.size()];
      HdedxVsPProfile    = new TProfile* [StdObjName.size()];
      HdedxVsEta         = new TH2D*     [StdObjName.size()];
      HdedxVsEtaProfile  = new TProfile* [StdObjName.size()];

      for (unsigned int i = 0; i < HitObjName.size(); i++){
         dEdxTemplate   [i] = (TH3F*) GetObjectFromPath (InputFile, (HitObjName[i] + "_ChargeVsPath").c_str());
         hit_MIP        [i] = (TH1D*) GetObjectFromPath (InputFile, (HitObjName[i] + "_Hit"         ).c_str());
	 //cout << (HitObjName[i] + "_ChargeVsPath").c_str()<<endl;
	 if (dEdxTemplate[i]){
	   cout<<"object found "<< (HitObjName[i] + "_ChargeVsPath").c_str()<<endl;
	   cout<<"name is "<<dEdxTemplate[i]->GetName()<<endl;
	 }
	 if (!dEdxTemplate[i]) cout<<"object not  found "<< (HitObjName[i] + "_ChargeVsPath").c_str()<<endl;

         Charge_Vs_XYH  [i] = new TH2D* [15];
         Charge_Vs_XYHN [i] = new TH2D* [15];
         Charge_Vs_XYLN [i] = new TH2D* [15];
         for (unsigned int g = 1; g < 16; g++){
            char Id[255]; sprintf (Id, "%02i", g);
            Charge_Vs_XYH  [i][g-1] = (TH2D*) GetObjectFromPath (InputFile, (HitObjName[i]+"_ChargeVsXYH"  + Id).c_str());
            Charge_Vs_XYHN [i][g-1] = (TH2D*) GetObjectFromPath (InputFile, (HitObjName[i]+"_ChargeVsXYHN" + Id).c_str());
            Charge_Vs_XYLN [i][g-1] = (TH2D*) GetObjectFromPath (InputFile, (HitObjName[i]+"_ChargeVsXYLN" + Id).c_str());
         }
      }


      for (unsigned int i = 0; i < StdObjName.size(); i++){
         HdedxMIP          [i] = (TH1D*)     GetObjectFromPath (InputFile, (StdObjName[i] + "_MIP"         ).c_str() );
         HdedxVsP          [i] = (TH2D*)     GetObjectFromPath (InputFile, (StdObjName[i] + "_dedxVsP"     ).c_str() );
         HdedxVsPSyst      [i] = (TH2D*)     GetObjectFromPath (InputFile, (StdObjName[i] + "_dedxVsPSyst" ).c_str() );
         HdedxVsPProfile   [i] = (TProfile*) GetObjectFromPath (InputFile, (StdObjName[i] + "_Profile"     ).c_str() );
         HdedxVsEta        [i] = (TH2D*)     GetObjectFromPath (InputFile, (StdObjName[i] + "_Eta2D"       ).c_str() );
         HdedxVsEtaProfile [i] = (TProfile*) GetObjectFromPath (InputFile, (StdObjName[i] + "_Eta"         ).c_str() );

         if (StdObjName[i].find("Ias")==string::npos){

            HMass        [i]  = (TH1D*) GetObjectFromPath(InputFile, (StdObjName[i] + "_Mass").c_str() );
            HProtonHitSO [i]  = (TH1D*) GetObjectFromPath(InputFile, (StdObjName[i] + "_ProtonHitSO").c_str() );
            HProtonHitPO [i]  = (TH1D*) GetObjectFromPath(InputFile, (StdObjName[i] + "_ProtonHitPO").c_str() );
            K [StdObjName[i]] = K_;  Kerr[StdObjName[i]] = Kerr_;
            C [StdObjName[i]] = C_;  Cerr[StdObjName[i]] = Cerr_;

            //while we're at it, get the constants
            if (type != 2){

               double Ktmp = K_; double Ctmp = C_; double KerrTmp = Kerr_; double CerrTmp = Cerr_;
               ExtractConstants (HdedxVsP[i], &Ktmp, &Ctmp, &KerrTmp, &CerrTmp);
               K [StdObjName[i]] = Ktmp;  Kerr[StdObjName[i]] = KerrTmp;
               C [StdObjName[i]] = Ctmp;  Cerr[StdObjName[i]] = CerrTmp;
               printf ("FINAL %s :: %s :: K = %.3lf\tC = %.3lf\n", LegEntry.c_str(), StdObjLegend[i].c_str(), Ktmp, Ctmp);
            }

            if (type == 2)
               HMassHSCP[i] = (TH1D*) GetObjectFromPath(InputFile, (StdObjName[i] + "_MassHSCP").c_str() );
            else
               HMassHSCP[i] = NULL;
         } else {

            HMass [i] = NULL;
            K[StdObjName[i]] = -1; Kerr[StdObjName[i]] = -1;
            C[StdObjName[i]] = -1; Cerr[StdObjName[i]] = -1;
         }
      }
   } // end of constructor
} dEdxPlotObj;

void getScaleFactor(TFile* InputFile1, TFile* InputFile2, string ObjName1, string ObjName2, string SaveDir, string Prefix, double* SFMip=NULL, double* SFProfile=NULL);
double GetMass(double P, double I, dEdxPlotObj* plotObj, string ObjName);
void MakeMapPlots(TH3F* Charge_Vs_Path3D, string ObjName, string SaveDir, string Prefix);
void MakeROCGeneral (TFile* InputFile1, TFile* InputFile2, vector<string> HistoNames, vector<string> LegendLabels, vector<Color_t> Colors, string SaveDir, string suffix, bool WithErrorBarse=false, bool Every2ndIsDashed=false);

void Draw2D (string SaveDir, vector <dEdxPlotObj*> plotObj);
void SuperposeFilesOnDeDxObj (string SaveDir, vector<dEdxPlotObj*> plotObj);
void CrossCompareAndControlPlots (string SaveDir, vector <dEdxPlotObj*> plotObj, string Reject, string Select);
void SaveKC (vector<dEdxPlotObj*> plotObj, string filename = "Report.txt");
void SystStudy(string SaveDir, vector<dEdxPlotObj*> plotObj, bool createTable=true, bool showChi2=false);
void ExtractSlope (string SaveDir, vector <dEdxPlotObj*> plotObj);
void HitPlots (string SaveDir, vector <dEdxPlotObj*> plotObj);

double GetMass(double P, double I, dEdxPlotObj* plotObj, string ObjName){
   return sqrt((I-plotObj->C[ObjName])/plotObj->K[ObjName])*P;
}


double GetMass(double P, double I, double* K, double* C){
  //  *C=4.591; // to remove
  // *K = 2.002; //to remove
   
   return sqrt((I-(*C))/(*K))*P;
}

TF1* GetMassLine(double M, dEdxPlotObj* plotObj, string ObjName, bool left=false)
{  
   double BetaMax = 0.9;
   double PMax = sqrt((BetaMax*BetaMax*M*M)/(1-BetaMax*BetaMax));

   double BetaMin = 0.2;
   double PMin = sqrt((BetaMin*BetaMin*M*M)/(1-BetaMin*BetaMin));
   
   if(left){PMax*=-1; PMin*=-1;}

   TF1* MassLine = new TF1("MassLine","[2] + ([0]*[0]*[1])/(x*x)", PMin, PMax);
   MassLine->SetParName  (0,"M");
   MassLine->SetParName  (1,"K");
   MassLine->SetParName  (2,"C");
   MassLine->SetParameter(0, M);
   MassLine->SetParameter(1, plotObj->K[ObjName]);
   MassLine->SetParameter(2, plotObj->C[ObjName]);
   MassLine->SetLineWidth(2);
   return MassLine;
}

TF1* GetMassLine(double M, double K, double C, bool left=false)
{  
   double BetaMax = 0.9;
   double PMax = sqrt((BetaMax*BetaMax*M*M)/(1-BetaMax*BetaMax));

   double BetaMin = 0.2;
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



void MakePlot()
{
   setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.10);
   gStyle->SetPadRightMargin (0.18);
   gStyle->SetPadLeftMargin  (0.125);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.35);
//   gStyle->SetPalette(51); 
   gStyle->SetPalette(1); 
   gStyle->SetNdivisions(510,"X");


   vector<string> HitObjName;                         vector<string> HitObjLegend;
   HitObjName.push_back("hit_PO");                    HitObjLegend.push_back("Pixel");
   HitObjName.push_back("hit_SO_in_noC_CCC");         HitObjLegend.push_back("Strip");
   HitObjName.push_back("hit_SP_in_noC_CCC");         HitObjLegend.push_back("Strip+Pixel charges"); //questo forse serve..

   vector<string> StdObjName;                         vector<string> StdObjLegend;
   StdObjName.push_back("Hybr2015_SP_in_noC_CCC");    StdObjLegend.push_back("hybrid2-15, SP");  // i think that is the one to keep

   StdObjName.push_back("Ias_SP_in_noC_CCC16");      StdObjLegend.push_back("Ias Strip+Pixel");  // to uncomment second step
   StdObjName.push_back("Ias_SO_in_noC_CCC16");      StdObjLegend.push_back("Ias Strip-Only");  // to uncomment second step 
   StdObjName.push_back("Ias_PO_in_noC_CCC16");      StdObjLegend.push_back("Ias Pixel-Only");   // to uncomment second step 


   vector <dEdxPlotObj*> plotObj;
 //  plotObj.push_back(new dEdxPlotObj("Histos_MCMinBias.root", "MC (MinBias)", "MCMinBias", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 1));
//   plotObj.push_back(new dEdxPlotObj("Histos_Data2015.root", "Data 2015",   "Data", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
//   plotObj.push_back(new dEdxPlotObj("Histos_Data2016.root", "Data 2016",   "Data", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
//   plotObj.push_back(new dEdxPlotObj("Histos_Run273158.root", "Run 273158",   "Run273158", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
//   plotObj.push_back(new dEdxPlotObj("Histos_Run273502.root", "Run 273502",   "Run273502", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
//////   plotObj.push_back(new dEdxPlotObj("Histos_Run278018.root", "Run 278018",   "Run278018", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
//////   plotObj.push_back(new dEdxPlotObj("Histos_Run278308.root", "Run 278308",   "Run278308", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
/////   plotObj.push_back(new dEdxPlotObj("Histos_Run279931.root", "Run 279931",   "Run279931", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
////   plotObj.push_back(new dEdxPlotObj("Histos_RunPostG.root", "Run PostG",   "RunPostG", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
/////   plotObj.push_back(new dEdxPlotObj("Histos_RunPostG.root", "Runs PostG",   "RunPostG", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
plotObj.push_back(new dEdxPlotObj("Histos_Region2.root", "Region2",   "Region2", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
plotObj.push_back(new dEdxPlotObj("Histos_Region3.root", "Region3",   "Region3", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
plotObj.push_back(new dEdxPlotObj("Histos_Region4.root", "Region4",   "Region4", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
//   plotObj.push_back(new dEdxPlotObj("Histos_Run280385.root", "Run 280385",   "Run280385", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 0));
//   plotObj.push_back(new dEdxPlotObj("Histos_MCDYM2600Q2.root",        "DY, Q = 2, M = 2.6TeV",    "DYM2600Q2",        HitObjName, StdObjName, HitObjLegend, StdObjLegend, 2));
//   plotObj.push_back(new dEdxPlotObj("Histos_MCGluino_M1000_f10.root", "Gluino, f=10, M = 1TeV",   "Gluino_M1000_f10", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 2));
//   plotObj.push_back(new dEdxPlotObj("Histos_MCGluino_M1800_f10.root", "Gluino, f=10, M = 1.8TeV", "Gluino_M1800_f10", HitObjName, StdObjName, HitObjLegend, StdObjLegend, 2));
//   plotObj.push_back(new dEdxPlotObj("Histos_MCStop_M1000.root",       "Stop, M = 1TeV",           "Stop_M1000",       HitObjName, StdObjName, HitObjLegend, StdObjLegend, 2));
//   plotObj.push_back(new dEdxPlotObj("Histos_MCGMStau_M494.root",      "GMStau, M = 494GeV",       "GMStau_M494",      HitObjName, StdObjName, HitObjLegend, StdObjLegend, 2));

   string SaveDir = "pictures_FAST/";
   system (string("rm -rf "+SaveDir+" && mkdir "+SaveDir).c_str());
   system ("rm -rf systematics && mkdir systematics");
   SaveKC (plotObj, SaveDir+"ConstantReport.txt");
   // copy K and C from SMMC to signal MC
   size_t MCIndex = 0;
   for (size_t m = 0; m < plotObj.size()-1; m++)
      if (plotObj[m]->type == 1){ MCIndex = m; break; }
   for (size_t i = 0; i < plotObj.size()-1; i++){
      if (plotObj[i]->type != 2) continue;
      for (size_t j = 0; j < plotObj[i]->StdObjName.size(); j++){
         if (plotObj[i]->StdObjName[j].find("Ias")!=string::npos) continue;
         plotObj[i]->K[plotObj[i]->StdObjName[j]] = plotObj[MCIndex]->K[plotObj[i]->StdObjName[j]];
         plotObj[i]->C[plotObj[i]->StdObjName[j]] = plotObj[MCIndex]->C[plotObj[i]->StdObjName[j]];
      }
   }

 
   cerr << "====== TESTA :: 2D Histograms ======" << endl;
   Draw2D (SaveDir, plotObj);

   cerr << "====== TESTB :: Standard Histos ======" << endl;
   SuperposeFilesOnDeDxObj (SaveDir, plotObj);
   // HitPlots (SaveDir, plotObj);

   cerr << "====== TESTC :: Cross-compare estimators ======" << endl;
   //   CrossCompareAndControlPlots (SaveDir, plotObj, "SO", "Hybr");
   
   cerr << "====== TESTD :: Systematics study ======" << endl;
   //SystStudy (SaveDir, plotObj);

   cerr << "====== TESTE :: Scale Factors ======" << endl;
   double SFMip, SFProfile;
   FILE* fout = fopen ((SaveDir+"ScaleFactors.txt").c_str(), "w");
   fprintf (fout, "=================\n >>> Old CCC >>>\n=================\n");
   for (size_t i = 0; i < plotObj.size(); i++){
      if(plotObj[i]->type==2) continue;
      getScaleFactor(plotObj[i]->InputFile, NULL, "hit_SO_in_noC_CCC", "hit_PO", SaveDir, plotObj[i]->SavePrefix, &SFMip, &SFProfile);   // here is where it computes the scale factor!
      fprintf (fout, (plotObj[i]->LegEntry + "\tPO to SO :: MIP:%.7lf\tProfile:%.7lf\n").c_str(), SFMip, SFProfile);
      if(i == MCIndex) continue;
      getScaleFactor (plotObj[i]->InputFile, plotObj[MCIndex]->InputFile, "hit_SO_in_noC_CCC", "", SaveDir, plotObj[i]->SavePrefix+"hit_SO", &SFMip, &SFProfile);
      fprintf (fout, (plotObj[i]->LegEntry + "\tSO MC to Data :: MIP:%.7lf\tProfile:%.7lf\n").c_str(), SFMip, SFProfile);
   }
//   fprintf (fout, "=================\n >>> New CCC >>>\n=================\n");
//   for (size_t i = 0; i < plotObj.size(); i++){
//      if(plotObj[i]->type==2) continue;
//      getScaleFactor(plotObj[i]->InputFile, NULL, "hit_SO_in_noC_newCCC", "hit_PO", SaveDir, plotObj[i]->SavePrefix, &SFMip, &SFProfile);
//      fprintf (fout, (plotObj[i]->LegEntry + "\tPO to SO :: MIP:%.7lf\tProfile:%.7lf\n").c_str(), SFMip, SFProfile);
//      if(i == MCIndex) continue;
//      getScaleFactor (plotObj[i]->InputFile, plotObj[MCIndex]->InputFile, "hit_SO_in_noC_newCCC", "", SaveDir, plotObj[i]->SavePrefix+"hit_SO", &SFMip, &SFProfile);
//      fprintf (fout, (plotObj[i]->LegEntry + "\tSO MC to Data :: MIP:%.7lf\tProfile:%.7lf\n").c_str(), SFMip, SFProfile);
//   }
   fclose (fout);
/*
   vector <string> ObjNames; vector <string> LegendLabels; vector <Color_t> Colors;
   ObjNames.push_back("harm2_SP_in_noC_CCC_MIP");    LegendLabels.push_back("I_{h}");  Colors.push_back(kRed);
   ObjNames.push_back("Ias_SP_in_noC_CCC_MIP");      LegendLabels.push_back("I_{as}"); Colors.push_back(kBlue);
   MakeROCGeneral (plotObj[2]->InputFile, plotObj[5]->InputFile, ObjNames, LegendLabels, Colors, SaveDir, "Estimators_LOIC");

   ObjNames.clear(); LegendLabels.clear(); Colors.clear();
   ObjNames.push_back("harm2_SP_in_noC_CCC_noF_MIP");     LegendLabels.push_back("harmonic-2");  Colors.push_back(kRed);
   ObjNames.push_back("Hybr201_SP_in_noC_CCC_noF_MIP");   LegendLabels.push_back("hybrid-2-10"); Colors.push_back(kRed);
   ObjNames.push_back("Hybr2015_SP_in_noC_CCC_noF_MIP");  LegendLabels.push_back("hybrid-2-15"); Colors.push_back(kRed);
   ObjNames.push_back("Hybr202_SP_in_noC_CCC_noF_MIP");   LegendLabels.push_back("hybrid-2-20"); Colors.push_back(kRed);
   ObjNames.push_back("Hybr2025_SP_in_noC_CCC_noF_MIP");  LegendLabels.push_back("hybrid-2-25"); Colors.push_back(kRed);
*/

//   vector <string> ObjNames; vector <string> LegendLabels; vector <Color_t> Colors;
//   ObjNames.push_back("Hybr2015_SP_in_noC_CCC_MIP");  LegendLabels.push_back("Hybrid-2-15"); Colors.push_back(kBlack);
//   ObjNames.push_back("Hybr2020_SP_in_noC_CCC_MIP");  LegendLabels.push_back("Hybrid-2-20"); Colors.push_back(kRed);
//   ObjNames.push_back("Hybr2025_SP_in_noC_CCC_MIP");  LegendLabels.push_back("Hybrid-2-25"); Colors.push_back(kBlue);
//   ObjNames.push_back("Hybr2030_SP_in_noC_CCC_MIP");  LegendLabels.push_back("Hybrid-2-30"); Colors.push_back(kGreen);
//   ObjNames.push_back("Hybr2035_SP_in_noC_CCC_MIP");  LegendLabels.push_back("Hybrid-2-35"); Colors.push_back(kCyan);
//   ObjNames.push_back("Ias_SP_in_noC_CCC_MIP");           LegendLabels.push_back("Ias (old)"); Colors.push_back(kMagenta);   //to uncomment later
//   ObjNames.push_back("Ias_SP_in_noC_CCC16_MIP");         LegendLabels.push_back("Ias (new)"); Colors.push_back(kAzure+2);   // to uncomment later
//   MakeROCGeneral (plotObj[8]->InputFile, plotObj[12]->InputFile, ObjNames, LegendLabels, Colors, SaveDir, "Estimators");  //to uncomment later
}

   /* 
   std::cout << "TESTD\n";

   getScaleFactor(InputFile, NULL, "harm2_SO_in_noC_CCC", "harm2_PO_raw", SaveDir, SaveName); // shift PO_raw to SO for File1
   if (InputFile2) {
      getScaleFactor(InputFile, InputFile2, "harm2_SO_in_noC_CCC", "", SaveDir, SaveName+SaveName2); // shift File2 to File1
      getScaleFactor(InputFile2, NULL, "harm2_SO_in_noC_CCC", "harm2_PO_raw", SaveDir, SaveName2);   // shift PO_raw to SO for File2
   }
}
*/

void getScaleFactor(TFile* InputFile1, TFile* InputFile2, string ObjName1, string ObjName2, string SaveDir, string Prefix,
      double* SFMip, double* SFProfile){
   TProfile*   HdedxVsPProfile1;
   TProfile*   HdedxVsPProfile2;
   TH1D*       HdedxMIP1;
   TH1D*       HdedxMIP2;

   if (InputFile2!=NULL){
      HdedxVsPProfile1 = (TProfile*)GetObjectFromPath(InputFile1, (ObjName1 + "_HitProfile"  ).c_str() );
      HdedxVsPProfile2 = (TProfile*)GetObjectFromPath(InputFile2, (ObjName1 + "_HitProfile"  ).c_str() );

      HdedxMIP1        = (TH1D*)    GetObjectFromPath(InputFile1, (ObjName1 + "_Hit"  ).c_str() );
      HdedxMIP2        = (TH1D*)    GetObjectFromPath(InputFile2, (ObjName1 + "_Hit"  ).c_str() );
   } else if (ObjName2!=""){
      HdedxVsPProfile1 = (TProfile*)GetObjectFromPath(InputFile1, (ObjName1 + "_HitProfile"  ).c_str() );
      HdedxVsPProfile2 = (TProfile*)GetObjectFromPath(InputFile1, (ObjName2 + "_HitProfile"  ).c_str() );

      HdedxMIP1        = (TH1D*)    GetObjectFromPath(InputFile1, (ObjName1 + "_Hit"  ).c_str() );
      HdedxMIP2        = (TH1D*)    GetObjectFromPath(InputFile1, (ObjName2 + "_Hit"  ).c_str() );
   } else return;

   TF1* mygausMIP = new TF1("mygausMIP","gaus", 1, 5);
   HdedxMIP1->Fit("mygausMIP","Q0","");
   double peakMIP1  = mygausMIP->GetParameter(1);
   HdedxMIP2->Fit("mygausMIP","Q0","");
   double peakMIP2  = mygausMIP->GetParameter(1);

   std::cout << "SCALE FACTOR WITH MIP     = " << peakMIP1/peakMIP2 << endl;
   *SFMip = peakMIP1/peakMIP2;

   TH1D* Chi2Dist = new TH1D("Chi2Dist","Chi2Dist",5000, 0.9 ,1.4);

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
         if(Momentum<5)continue;//|| Momentum>20)continue;
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
   *SFProfile = AbsGain;

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
   leg->AddEntry (HdedxMIP1, "Data", "L");
   leg->AddEntry (HdedxMIP2, "unscaled MC", "L");
   leg->AddEntry (HdedxMIP4, "scaled   MC", "L");
   leg->Draw();
   SaveCanvas(c1, SaveDir, "Rescale_"+Prefix+"_"+ObjName1+ObjName2 + "_Hit");
   delete c1;


   c1 = new TCanvas("c1", "c1", 600,600);
   Chi2Dist->SetStats(kFALSE);
   Chi2Dist->GetXaxis()->SetNdivisions(504);
   Chi2Dist->GetXaxis()->SetTitle("Rescale Factor");
   Chi2Dist->GetYaxis()->SetTitle("Weighted Square Distance");
   Chi2Dist->Draw("");
   c1->SetLogy(true);
   c1->SetGridx(true); 
   SaveCanvas(c1, SaveDir, "Rescale_"+Prefix+"_"+ObjName1+ObjName2 + "_Dist");
   delete c1;

   c1 = new TCanvas("c1", "c1", 600,600);
   leg = new TLegend (0.30, 0.25, 0.80, 0.55);
   leg->SetHeader ("Fitting the Profile");
   leg->SetFillColor(0);
   leg->SetFillStyle(0);
   leg->SetBorderSize(0);
   HdedxVsPProfile1->SetStats(kFALSE);
   HdedxVsPProfile1->SetAxisRange(5,50,"X");
   HdedxVsPProfile1->SetAxisRange(3.0,5.5,"Y");
   HdedxVsPProfile1->GetXaxis()->SetTitle("track momentum (GeV)");
   HdedxVsPProfile1->GetYaxis()->SetTitle("dE/dx (MeV/cm)");
   HdedxVsPProfile1->SetMarkerColor(2);
   HdedxVsPProfile1->Draw("");

   HdedxVsPProfile2->SetMarkerColor(1);
   HdedxVsPProfile2->Draw("same");
   TProfile* HdedxVsPProfile3 = (TProfile*)HdedxVsPProfile2->Clone("abc");
   HdedxVsPProfile3->SetMarkerColor(8);
   HdedxVsPProfile3->Scale(peakMIP1/peakMIP2);
//   HdedxVsPProfile3->Draw("same");
   TProfile* HdedxVsPProfile4 = (TProfile*)HdedxVsPProfile2->Clone("afs");
   HdedxVsPProfile4->SetMarkerColor(4);
   HdedxVsPProfile4->Scale(AbsGain);
   HdedxVsPProfile4->Draw("same");
   leg->AddEntry (HdedxVsPProfile1, "Data", "LP");
   leg->AddEntry (HdedxVsPProfile2, "unscaled MC", "LP");
   leg->AddEntry (HdedxVsPProfile4, "scaled MC",   "LP");
   leg->Draw();

   DrawPreliminary("", 13, "");
   SaveCanvas(c1, SaveDir, "Rescale_"+Prefix+"_"+ObjName1+ObjName2 + "_Profile");
   delete c1;
}


void getScaleFactor_new(TFile* InputFile1, TFile* InputFile2, string ObjName1, string ObjName2, string SaveDir, string Prefix,
      double* SFMip, double* SFProfile){
    return;
}


void ExtractConstants(TH2D* input, double* K, double* C, double* Kerr, double* Cerr,
      double MinRange, double MaxRange, double MassCenter, double LeftMassMargin, double RightMassMargin)
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
		if (inputnew->GetYaxis()->GetBinCenter(y)<2 || inputnew->GetYaxis()->GetBinCenter(y)>4.2) inputnewPion->SetBinContent(x,y,0);
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
	      SaveCanvas(c1, "fit/", "dedxVsP");
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
		  //cout<<x<<"  "<<mygausPion->GetParameter(1)<<endl;
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
		  if(P>=MinRange && P<=MaxRange){SaveCanvas(c1,"../DeDxStudy/",buffer);}
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
	       //fitC->SetParameter(0, 2.5);
	       fitC->SetParLimits(0, 2,4);
	       fitC->SetLineWidth(2);
               fitC->SetLineColor(2);
	       cout<<"prima del fit"<<endl;
               FitResultPion->Fit("fitC", "M R E I 0");
	       cout<<"dopo il fit"<<endl;
               fitC->SetRange(1,4);
               fitC->Draw("same");
	       cout<<"ha fatto il fit"<<endl;
	       *C    = fitC->GetParameter(0);
	       *Cerr = fitC->GetParError(0);
	       cout<< "il valore e`: "<<*C<<"   " <<*Cerr<< endl;




	       //   TF1* myfit = new TF1("myfit","[1]+(pow(0.93827,2) + x*x)/([0]*x*x)", MinRange, MaxRange);
	       TF1* myfit = new TF1("myfit","[0]*pow(1.8756/x,2) + [1]", MinRange, MaxRange); //1875.6 MeV  deuteron mass
	       myfit->SetParName  (0,"K");
	       myfit->SetParName  (1,"C");
	       myfit->SetParameter(0, 1.8);
	       myfit->SetParameter(1, *C);
	       myfit->SetParLimits(0, 1.3,3.0); //
	       myfit->SetParLimits(1, *C,*C);
	       // myfit->SetParLimits(1, 3.8,3.8); // 
	       myfit->SetLineWidth(2);
	       myfit->SetLineColor(2);
	       FitResult->Fit("myfit", "M R E I 0");
	       myfit->SetRange(MinRange,MaxRange);
	       myfit->Draw("same");

	       //double prevConstants [] = {*K, *Kerr, *C, *Cerr};
	       *K    = myfit->GetParameter(0);
	       //*C    = myfit->GetParameter(1);
	       *Kerr = myfit->GetParError(0);
	       //*Cerr = myfit->GetParError(1);

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
	       sprintf(buffer,"C = %4.3f +- %6.4f",myfit->GetParameter(1), myfit->GetParError(1));
	       st->AddText(buffer);
	       st->Draw("same");
	       sprintf(buffer,"%sFit","fit/");
	       SaveCanvas(c1,"../DeDxStudy/",buffer);              
	       delete c1;

          delete line1;
          delete line2;
          delete myfit;
          delete FitResult;
          delete inputnew;
       }
}


void SystStudy(string SaveDir, vector<dEdxPlotObj*> plotObj, bool createTable, bool showChi2){
   double MinRange = 0.8;
   double MaxRange = 1.8;
   char buffer[2048];
   for (size_t j = 0; j < plotObj[0]->StdObjName.size(); j++){
      bool isEstim = (plotObj[0]->StdObjName[j].find("Ias")==string::npos);
   
      vector <TH2D*> inputnew;
      vector <Color_t> colors;
      vector <string> SaveNames;
      vector <string> SavePrefix;
      colors.push_back(kBlack);
      colors.push_back(kBlue);
      colors.push_back(kGreen);
      colors.push_back(kRed);
      colors.push_back(kCyan);
      colors.push_back(kSpring-7);
      colors.push_back(kPink+9);
      colors.push_back(kOrange+8);
      colors.push_back(kYellow);
      colors.push_back(kAzure+5);
      size_t MCIndex = 0;
      for (size_t i = 0; i < plotObj.size(); i++){
         if (plotObj[i]->type == 2) continue;
         if (plotObj[i]->type == 1 && MCIndex != i) MCIndex = i;
         char buffer2[20]; sprintf (buffer2, "tempTH2D_%lu", i);
         inputnew.push_back((TH2D*) plotObj[i]->HdedxVsPSyst[j]->Clone(buffer2));
         inputnew[inputnew.size()-1]->Rebin2D(5,10);
         if (isEstim){
            for(int x=1;x<=inputnew[inputnew.size()-1]->GetNbinsX();x++){
            for(int y=1;y<=inputnew[inputnew.size()-1]->GetNbinsY();y++){
               double Mass = GetMass(inputnew[inputnew.size()-1]->GetXaxis()->GetBinCenter(x),inputnew[inputnew.size()-1]->GetYaxis()->GetBinCenter(y), plotObj[i], plotObj[i]->StdObjName[j]);
               if (Mass < 1.88-0.40 || Mass > 1.88+0.40) inputnew[inputnew.size()-1]->SetBinContent (x,y,0);  //Mass was at 0.938
            }}
         }
         SaveNames.push_back (plotObj[i]->LegEntry);
         SavePrefix.push_back (plotObj[i]->SavePrefix);
      }

      double** means            = new double* [inputnew.size()];
      double** sigmas           = new double* [inputnew.size()];
      double** momenta          = new double* [inputnew.size()];

      double** means_err        = new double* [inputnew.size()];
      double** sigmas_err       = new double* [inputnew.size()];
      double** momenta_err      = new double* [inputnew.size()];

      double** ratio_means      = new double* [inputnew.size()];
      double** ratio_means_err  = new double* [inputnew.size()];
      double** ratio_sigmas     = new double* [inputnew.size()];
      double** ratio_sigmas_err = new double* [inputnew.size()];

      for (size_t i = 0; i < inputnew.size(); i++){
         TCanvas* c1 = new TCanvas("canvas", "canvas", 600,600);
         c1->SetLogz(true);
	 TH2D htmp ("2Dtmp", "2Dtmp", 1, 0, 1.5, 1, isEstim?4:0, isEstim?16:1);
	 htmp.Draw("");
         inputnew[i]->Draw("same COLZ");
         SaveCanvas (c1, "systematics/", "Syst_"+SavePrefix[i]+"_"+plotObj[0]->StdObjName[j]+"_inputnew");
         delete c1;

         means[i]             = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];
         sigmas[i]            = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];
         momenta[i]           = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];

         means_err[i]         = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];
         sigmas_err[i]        = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];
         momenta_err[i]       = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];

         ratio_means[i]       = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];
         ratio_means_err[i]   = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];

         ratio_sigmas[i]      = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];
         ratio_sigmas_err[i]  = new double [inputnew[0]->GetXaxis()->FindBin(MaxRange) - inputnew[0]->GetXaxis()->FindBin(MinRange)];   
      }

      TF1**  mygaus     = new TF1*  [inputnew.size()];
      TH1D** Projection = new TH1D* [inputnew.size()];

      for(int x=inputnew[0]->GetXaxis()->FindBin(MinRange); x<inputnew[0]->GetXaxis()->FindBin(MaxRange); x++){
         double P = inputnew[0]->GetXaxis()->GetBinCenter(x);
         char momentumP    [20]; sprintf (momentumP,    "%d",    (int) (100*P));
         char momentumLow  [20]; sprintf (momentumLow,  "%.2lf", P-0.05);
         char momentumHigh [20]; sprintf (momentumHigh, "%.2lf", P+0.05);
         TCanvas* c1  = new TCanvas("canvas", "canvas", 600,600);
         TLegend* leg = new TLegend(0.40, 0.65, 0.95, 0.89);
         leg->SetFillColor(0);
         leg->SetFillStyle(0);
         leg->SetBorderSize(0);
         TH1D h ("temp", "temp", 1, 0.1, isEstim?15.0:1.0);
         h.GetXaxis()->SetTitle(isEstim?"dE/dx estimator (MeV/cm)":"I_{as} discriminator");
         h.GetYaxis()->SetTitle("arbitrary units");
         h.SetTitle("");
         h.SetStats(0);
         h.SetTitleOffset(1.30);
         double maxY = 0;
         vector <string> chi2;
         for(int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){

            Projection[fileIndex] = (TH1D*)(inputnew[fileIndex]->ProjectionY("proj",x,x))->Clone();
//            if(Projection[fileIndex]->Integral()<100)continue;
            Projection[fileIndex]->Sumw2();
            Projection[fileIndex]->Scale(1.0/Projection[fileIndex]->Integral());
            Projection[fileIndex]->SetMarkerColor(colors[fileIndex]);
            Projection[fileIndex]->SetLineColor(colors[fileIndex]);

            double maximum = Projection[fileIndex]->GetMaximumBin()*((isEstim?15.0:1.0)- 0.0)/Projection[fileIndex]->GetNbinsX();
            mygaus[fileIndex] = new TF1("mygaus","gaus",
                  isEstim?(maximum*(1 - 0.10)):(maximum - 0.03 - 0.09*sin(maximum*M_PI)*sin(maximum*M_PI)),
                  isEstim?(maximum*(1 + 0.10)):(maximum + 0.03 + 0.09*sin(maximum*M_PI)*sin(maximum*M_PI)));
            mygaus[fileIndex]->SetParameter(1, maximum);
            mygaus[fileIndex]->SetParameter(2, 0.05);
            mygaus[fileIndex]->SetLineColor(colors[fileIndex]);
            Projection[fileIndex]->Fit("mygaus","Q0 RME");
            double chiFromFit  = (mygaus[fileIndex]->GetChisquare())/(mygaus[fileIndex]->GetNDF());

            means  [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] = mygaus[fileIndex]->GetParameter(1);
            sigmas [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] = mygaus[fileIndex]->GetParameter(2);
            momenta[fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] = P;

            means_err  [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] = mygaus[fileIndex]->GetParError(1);
            sigmas_err [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] = mygaus[fileIndex]->GetParError(2);
            momenta_err[fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] = 0.05;

            cerr << "means      = " << means      [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] << endl;
            cerr << "sigmas     = " << sigmas     [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] << endl;
            cerr << "momenta    = " << momenta    [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] << endl;
            cerr << "means_err  = " << means_err  [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] << endl;
            cerr << "sigmas_err = " << sigmas_err [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] << endl;
            cerr << "r-chi2     = " << chiFromFit << endl;
            cerr << "---------------------------" << endl;

//         TText text1, text2;

            if (Projection[fileIndex]->GetMaximum()*1.4 > maxY)
               maxY = Projection[fileIndex]->GetMaximum()*1.4;
            char drawStuff [2048]; sprintf (drawStuff, "#chi^{2} = %.1lf", chiFromFit);
            chi2.push_back(string(drawStuff));
//         char draw1 [2048];     sprintf (draw1, "#mu_{1} = %lf", maximum);
//         text1.DrawText (0.2, Projection->GetMaximum(), drawStuff);
//         text2.DrawText (0.2, Projection->GetMaximum()*0.5, draw1);
         }
         h.GetYaxis()->SetRangeUser(1E-5, maxY*1.15);
         h.GetXaxis()->SetRangeUser (0.1,isEstim?15:1.0);
         h.Draw();
         leg->SetHeader(((string)momentumLow + "GeV < p < " + (string)momentumHigh+"GeV").c_str());
         for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
            Projection[fileIndex]->Draw("same EP");
            mygaus[fileIndex]->Draw("same L");
            leg->AddEntry (Projection[fileIndex], (SaveNames[fileIndex]+(showChi2?", "+chi2[fileIndex]:"")).c_str(), "LP");
         }
         leg->Draw();
         DrawPreliminary("", 13.0, plotObj[0]->StdObjLegend[j]);
         SaveCanvas (c1, "systematics/", "Syst_"+plotObj[0]->StdObjName[j]+"_"+string(momentumP));
         delete c1;
         delete leg;
         for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
            delete Projection[fileIndex];
            delete mygaus[fileIndex];
         }
      }

      delete [] mygaus;
      delete [] Projection;

      for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
      for(int x=inputnew[fileIndex]->GetXaxis()->FindBin(MinRange); x<inputnew[fileIndex]->GetXaxis()->FindBin(MaxRange); x++){
            ratio_means     [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] =
               means[MCIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] / means [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)];

            ratio_means_err [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] =
               ratio_means  [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] *
               sqrt(
                     pow (means_err[MCIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] / means[MCIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)],2) +
                     pow (means_err[fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] / means [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)], 2)
                  );

            ratio_sigmas     [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] =
               sigmas[MCIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] / sigmas [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)];

            ratio_sigmas_err [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] =
               ratio_sigmas  [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] *
               sqrt(
                     pow (sigmas_err[MCIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] / sigmas[MCIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)],2) +
                     pow (sigmas_err[fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] / sigmas [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)], 2)
                  );
      }}

      if (createTable){
         FILE * fmeans  = fopen (("systematics/reports_"+plotObj[0]->StdObjName[j]+"_means.txt" ).c_str(), "w");
         FILE * fsigmas = fopen (("systematics/reports_"+plotObj[0]->StdObjName[j]+"_sigmas.txt").c_str(), "w");
         fprintf (fmeans,  "\\begin{table}\\centering\n  \\caption{}\n  \\begin{tabular}{| c || c | c | c | c |}\\hline\n     $p$ (GeV)");
         fprintf (fsigmas, "\\begin{table}\\centering\n  \\caption{}\n  \\begin{tabular}{| c || c | c | c | c |}\\hline\n     $p$ (GeV)");
         for (int fileIndex = SaveNames.size()-1; fileIndex >=0; fileIndex--){
            fprintf (fmeans,  " & %s", SaveNames[fileIndex].c_str());
            fprintf (fsigmas, " & %s", SaveNames[fileIndex].c_str());
         }
         fprintf (fmeans,  " \\\\ \\hline\n");
         fprintf (fsigmas, " \\\\ \\hline\n");
         for(int x=inputnew[0]->GetXaxis()->FindBin(MinRange); x<inputnew[0]->GetXaxis()->FindBin(MaxRange); x++){
            fprintf (fmeans,  "     %.2lf -- %.2lf", momenta[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)]-momenta_err[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)], momenta[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)]+momenta_err[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)]);
            fprintf (fsigmas, "     %.2lf -- %.2lf", momenta[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)]-momenta_err[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)], momenta[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)]+momenta_err[0][x - inputnew[0]->GetXaxis()->FindBin(MinRange)]);
            for (int fileIndex = SaveNames.size()-1; fileIndex >=0; fileIndex--){
               fprintf (fmeans,  "   &   % .2lf", means   [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)]);
               fprintf (fsigmas, "   &   % .2lf", sigmas  [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)]);
            }
            fprintf (fmeans,  " \\\\%%");
            fprintf (fsigmas, " \\\\%%");
            for (int fileIndex = SaveNames.size()-2; fileIndex >=0; fileIndex--){
               fprintf (fmeans,  "   &   % .4lf %%", 100*(ratio_means   [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] - 1.0));
               fprintf (fsigmas, "   &   % .4lf %%", 100*(ratio_sigmas  [fileIndex][x - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange)] - 1.0));
            }
            fprintf (fmeans,  "\n");
            fprintf (fsigmas, "\n");
         }
         fprintf (fmeans,  "  \\hline\n");
         fprintf (fsigmas, "  \\hline\n");
         fprintf (fmeans,  "  \\end{tabular}\n\\end{table}");
         fprintf (fsigmas, "  \\end{tabular}\n\\end{table}");
         fclose(fmeans);
         fclose(fsigmas);
      }

      TGraphErrors** gMeans = new TGraphErrors* [SaveNames.size()];
      TGraphErrors** gSigmas = new TGraphErrors* [SaveNames.size()];
      TGraphErrors** gRatio = new TGraphErrors* [SaveNames.size()];

      // Means
      TCanvas* c1  = new TCanvas("canvas", "canvas", 600,600);
      TLegend* leg = new TLegend(0.50, 0.80, 0.80, 0.90);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      TH1D h ("temp", "temp", 1, MinRange, MaxRange);
      h.GetXaxis()->SetTitle("p (GeV)");
      h.GetYaxis()->SetTitle(isEstim?"#mu (I_{h}) (MeV/cm)":"#mu (I_{as})");
      h.GetXaxis()->SetRangeUser (MinRange, MaxRange);
      h.SetAxisRange (0, isEstim?15:1.0, "Y");
      h.GetYaxis()->SetTitleOffset(1.35);
      h.SetStats(0);
      h.Draw();
      for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
         gMeans[fileIndex] = new TGraphErrors (inputnew[fileIndex]->GetXaxis()->FindBin(MaxRange) - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange), momenta[fileIndex], means[fileIndex], momenta_err[fileIndex], means_err[fileIndex]);
         gMeans[fileIndex]->SetMarkerColor (colors[fileIndex]);
         gMeans[fileIndex]->SetLineColor (colors[fileIndex]);
         gMeans[fileIndex]->SetMarkerStyle (20);
         gMeans[fileIndex]->Draw("same ELP");
      }
      leg->AddEntry (gMeans[0], SaveNames[0].c_str(), "LP");
      leg->AddEntry (gMeans[1], SaveNames[1].c_str(), "LP");
      leg->AddEntry (gMeans[2], SaveNames[2].c_str(), "LP");
      leg->AddEntry (gMeans[3], SaveNames[3].c_str(), "LP");
      leg->Draw();
      DrawPreliminary("", 13.0, plotObj[0]->StdObjLegend[j]);
      SaveCanvas (c1, SaveDir, "Syst_"+plotObj[0]->StdObjName[j]+"_mean");
      delete c1;
      delete leg;

      // Sigmas
      c1  = new TCanvas("canvas", "canvas", 600,600);
      leg = new TLegend(0.50, 0.80, 0.80, 0.90);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      h.GetYaxis()->SetTitle (isEstim?"#sigma (I_{h}) (MeV/cm)":"#sigma (I_{as})");
      h.SetAxisRange (0, isEstim?1.5:0.2, "Y");
      h.GetXaxis()->SetRangeUser (MinRange, MaxRange);
      h.Draw();
      for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
         gSigmas[fileIndex] = new TGraphErrors (inputnew[fileIndex]->GetXaxis()->FindBin(MaxRange) - inputnew[fileIndex]->GetXaxis()->FindBin(MinRange), momenta[fileIndex], sigmas[fileIndex], momenta_err[fileIndex], sigmas_err[fileIndex]);
         gSigmas[fileIndex]->SetMarkerColor (colors[fileIndex]);
         gSigmas[fileIndex]->SetLineColor (colors[fileIndex]);
         gSigmas[fileIndex]->SetMarkerStyle (20);
         gSigmas[fileIndex]->Draw("same ELP");
      }
      leg->AddEntry (gSigmas[0], SaveNames[0].c_str(), "LP");
      leg->AddEntry (gSigmas[1], SaveNames[1].c_str(), "LP");
      leg->AddEntry (gSigmas[2], SaveNames[2].c_str(), "LP");
      leg->AddEntry (gSigmas[3], SaveNames[3].c_str(), "LP");
      leg->Draw();
      DrawPreliminary("", 13.0, plotObj[0]->StdObjLegend[j]);
      SaveCanvas (c1, SaveDir, "Syst_"+plotObj[0]->StdObjName[j]+"_mean");
      delete c1;
      delete leg;

      TPad* t1 = NULL;
      TPad* t2 = NULL;
      // first pad
      double upEstim, upNoEstim;
      string labelEstim, labelNoEstim;
      for (uint8_t plot = 0; plot < 2; plot++){
         c1  = new TCanvas ("canvas", "canvas", 600, 600);
         leg = new TLegend (0.50, 0.60, 0.80, 0.90);
         leg->SetFillColor(0);
         leg->SetFillStyle(0);
         leg->SetBorderSize(0);
         t1 = new TPad ("t1", "t1", 0.0, 0.20, 1.0, 1.0);
         t1->Draw();
         t1->cd();
         t1->SetTopMargin(0.06);
         h.Reset();
         h.SetTitle("");
         h.SetStats(0);
         h.GetXaxis()->SetRangeUser (MinRange, MaxRange);
         h.GetXaxis()->SetTitle ("p (GeV)");
         if (plot == 0){
            upEstim = 15;
            upNoEstim = 1.0;
            labelEstim = "#mu (I_{h}) (MeV/cm)";
            labelNoEstim = "#mu (I_{as})";
         } else {
            upEstim = 1.0;
            upNoEstim = 0.2;
            labelEstim = "#sigma (I_{h}) (MeV/cm)";
            labelNoEstim = "#sigma (I_{as})";
         }
         h.GetYaxis()->SetRangeUser (0.0, isEstim?upEstim:upNoEstim);
         h.GetYaxis()->SetTitle (isEstim?labelEstim.c_str():labelNoEstim.c_str());
         h.Draw();

         if (plot == 0){
            for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
            gMeans[fileIndex]->Draw("same ELP");
            leg->AddEntry (gMeans[fileIndex],  SaveNames[fileIndex] .c_str(), "LP");
            }
         } else {
            for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
            gSigmas[fileIndex]->Draw("same ELP");
            leg->AddEntry (gSigmas[fileIndex],  SaveNames[fileIndex] .c_str(), "LP");
            }
         }
         leg->Draw();
         DrawPreliminary("", 13.0, plotObj[0]->StdObjLegend[j]);

         // now the 2nd pad
         c1->cd();
         t2 = new TPad ("t2", "t2", 0.0, 0.0, 1.0, 0.2);
         t2->Draw();
         t2->cd();
         t2->SetGridy(true);
         t2->SetPad(0.0, 0.0, 1.0, 0.2);
         t2->SetTopMargin(0);
         t2->SetBottomMargin(0.5);

         TH1D h2 ("temp2", "temp2", 1, MinRange, MaxRange);
         h2.SetTitle("");
         h2.SetStats(0);
         h2.GetXaxis()->SetRangeUser (0.8, 2.0);
         h2.GetXaxis()->SetTitle ("");
         h2.GetXaxis()->SetLabelFont (43);
         h2.GetXaxis()->SetLabelSize (15);
         h2.GetXaxis()->SetTitleFont (43);
         h2.GetXaxis()->SetTitleSize (15);
         h2.GetYaxis()->SetRangeUser (0.85, 1.15);
         if (plot == 1){
            h2.GetYaxis()->SetRangeUser (0.40, 1.60);
	      }
         h2.GetYaxis()->SetTitle ("");
         h2.GetYaxis()->SetLabelFont (43);
         h2.GetYaxis()->SetLabelSize (15);
         h2.GetYaxis()->SetTitleFont (43);
         h2.GetYaxis()->SetTitleSize (15);
         h2.Draw("same");
         TLine LineAtOne (0.6, 1, 1.2, 1);
         LineAtOne.SetLineStyle (3);
         LineAtOne.SetLineWidth (1);
         LineAtOne.SetLineColor (kBlack);
         LineAtOne.Draw("same");

         // finally the ratio graph
         if (plot == 0){
            for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
               gRatio[fileIndex] = new TGraphErrors (inputnew[fileIndex]->GetXaxis()->FindBin(MaxRange)-inputnew[fileIndex]->GetXaxis()->FindBin(MinRange), momenta[fileIndex], ratio_means[fileIndex], momenta_err[fileIndex], ratio_means_err[fileIndex]);
               gRatio[fileIndex]->SetMarkerColor (colors[fileIndex]);
               gRatio[fileIndex]->SetLineColor (colors[fileIndex]);
               gRatio[fileIndex]->SetMarkerStyle (23);
               if(fileIndex == MCIndex) continue;
               gRatio[fileIndex]->Draw("same EP");
            }
         } else{
            for (int fileIndex = 0; fileIndex < SaveNames.size(); fileIndex++){
               gRatio[fileIndex] = new TGraphErrors (inputnew[fileIndex]->GetXaxis()->FindBin(MaxRange)-inputnew[fileIndex]->GetXaxis()->FindBin(MinRange), momenta[fileIndex], ratio_sigmas[fileIndex], momenta_err[fileIndex], ratio_sigmas_err[fileIndex]);
               gRatio[fileIndex]->SetMarkerColor (colors[fileIndex]);
               gRatio[fileIndex]->SetLineColor (colors[fileIndex]);
               gRatio[fileIndex]->SetMarkerStyle (23);
               if(fileIndex == MCIndex) continue;
               gRatio[fileIndex]->Draw("same EP");
            }
         }

         SaveCanvas (c1, SaveDir, "Syst_"+plotObj[0]->StdObjName[j]+(!plot?"_means_WRatio":"_sigma_WRatio"));

         delete t1;
         delete t2;
         delete c1;
         delete leg;
         for (size_t i = 0; i < SaveNames.size(); i++)
            delete gRatio[i];
      }

      for (size_t i = 0; i < SaveNames.size(); i++){
         delete inputnew[i];
         delete gMeans[i];
         delete gSigmas[i];
         delete [] means[i];
         delete [] sigmas[i];
         delete [] momenta[i];
         delete [] means_err[i];
         delete [] sigmas_err[i];
         delete [] momenta_err[i];
         delete [] ratio_means[i];
         delete [] ratio_means_err[i];
         delete [] ratio_sigmas[i];
         delete [] ratio_sigmas_err[i];
      }
      delete [] gMeans;
      delete [] gSigmas;
      delete [] gRatio;
      delete [] means;
      delete [] sigmas;
      delete [] momenta;
      delete [] means_err;
      delete [] sigmas_err;
      delete [] momenta_err;
      delete [] ratio_means;
      delete [] ratio_means_err;
      delete [] ratio_sigmas;
      delete [] ratio_sigmas_err;
   }
}


void CompareDeDx (TFile* InputFile, string SaveDir, string SaveName, string ObjName1, string ObjName2){
   if (ObjName1.find("hit")==string::npos && ObjName2.find("hit")==string::npos){
      TProfile*   HdedxVsEtaProfile1  = (TProfile*)  GetObjectFromPath(InputFile, (ObjName1 + "_Eta" ).c_str() );
      TProfile*   HdedxVsEtaProfile2  = (TProfile*)  GetObjectFromPath(InputFile, (ObjName2 + "_Eta" ).c_str() );
      TH1D*       HdedxMIP1           = (TH1D*)      GetObjectFromPath(InputFile, (ObjName1 + "_MIP" ).c_str() );
      TH1D*       HdedxMIP2           = (TH1D*)      GetObjectFromPath(InputFile, (ObjName2 + "_MIP" ).c_str() );
	
   	TCanvas* c1  = new TCanvas("c1", "c1", 600,600);
   	TLegend* leg = new TLegend(0.50, 0.80, 0.80, 0.90);
   	leg->SetHeader (SaveName.c_str());
   	leg->SetFillColor(0);
   	leg->SetFillStyle(0);
   	leg->SetBorderSize(0);
   	leg->AddEntry (HdedxVsEtaProfile1, ObjName1.c_str(), "P");
   	leg->AddEntry (HdedxVsEtaProfile2, ObjName2.c_str(), "P");
   	HdedxVsEtaProfile1->SetStats(kFALSE);
   	HdedxVsEtaProfile2->SetMarkerStyle(23);
   	HdedxVsEtaProfile1->SetMarkerColor(kBlack);
   	HdedxVsEtaProfile2->SetMarkerColor(kBlue);
   	HdedxVsEtaProfile1->GetXaxis()->SetTitle("#eta");
   	HdedxVsEtaProfile1->GetYaxis()->SetTitle(ObjName1.find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
   	HdedxVsEtaProfile1->Draw("");
   	HdedxVsEtaProfile2->Draw("same");
   	leg->Draw();
      DrawPreliminary("", 13, "");
   	SaveCanvas(c1, SaveDir, "Comparison"+SaveName+"_"+ObjName1+"_"+ObjName2+"_HdedxVsEtaProfile");
   	delete leg;
   	delete c1;

   	c1 = new TCanvas("c1", "c1", 600,600);
   	leg = new TLegend (0.50, 0.80, 0.80, 0.90);
   	c1->SetLogy(true);
   	c1->SetGridx(true);
   	leg->SetHeader (SaveName.c_str());
   	leg->SetFillColor(0);
   	leg->SetFillStyle(0);
   	leg->SetBorderSize(0);
   	HdedxMIP1->SetStats(kFALSE);
   	HdedxMIP1->GetXaxis()->SetTitle(ObjName1.find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
   	HdedxMIP1->GetYaxis()->SetTitle("fraction of tracks");
   	HdedxMIP2->Scale(1.0/HdedxMIP2->Integral());
   	HdedxMIP1->Scale(1.0/HdedxMIP1->Integral());
   	HdedxMIP1->GetXaxis()->SetRangeUser(0,8);
   	HdedxMIP1->GetYaxis()->SetRangeUser(2e-5,1);
   	HdedxMIP1->SetLineColor (kBlack);
   	HdedxMIP2->SetLineColor (kBlue);
   	leg->AddEntry (HdedxMIP1, ObjName1.c_str(), "L");
   	leg->AddEntry (HdedxMIP2, ObjName2.c_str(), "L");
   	HdedxMIP1->Draw("");
   	HdedxMIP2->Draw("same");
   	leg->Draw();
      DrawPreliminary("", 13, "");
      SaveCanvas(c1, SaveDir, "Comparison"+SaveName+"_"+ObjName1+"_"+ObjName2+"_MIP", true);
      delete leg;
      delete c1;

   } else if (ObjName1.find("hit")!=string::npos && ObjName2.find("hit")!=string::npos){
      for (unsigned int g=0;g<16;g++){
         char Id[255]; sprintf (Id, "%02i", g);
         TH2D* Charge_Vs_XYLN1 = (TH2D*) GetObjectFromPath (InputFile, (ObjName1 + "_ChargeVsXYLN" + Id).c_str());
         TH2D* Charge_Vs_XYLN2 = (TH2D*) GetObjectFromPath (InputFile, (ObjName2 + "_ChargeVsXYLN" + Id).c_str());
         TH1D* ProjX1          = Charge_Vs_XYLN1->ProjectionX (("X1_"+string(Id)).c_str());
         TH1D* ProjY1          = Charge_Vs_XYLN1->ProjectionY (("Y1_"+string(Id)).c_str());
         TH1D* ProjX2          = Charge_Vs_XYLN2->ProjectionX (("X2_"+string(Id)).c_str());
         TH1D* ProjY2          = Charge_Vs_XYLN2->ProjectionY (("Y2_"+string(Id)).c_str());

         TCanvas* c1  = new TCanvas ("c1", "c1", 600, 600);
         TLegend* leg = new TLegend (0.50, 0.75, 0.80, 0.90);
         c1->SetLogy (true);
         leg->SetHeader (("Module No. " + string(Id)).c_str());
         leg->SetHeader (SaveName.c_str());
         leg->SetFillColor(0);
         leg->SetFillStyle(0);
         leg->SetBorderSize(0);
         ProjX1->SetStats(kFALSE);
         ProjX1->SetLineColor (kBlack);
         ProjX2->SetLineColor (kBlue);
         ProjX1->GetXaxis()->SetTitle("normalized x coordinate");
         ProjX1->GetYaxis()->SetTitle("number of hits");
         ProjX1->SetAxisRange (-1.5, 1.5, "X");
         ProjX1->Draw("L");
         ProjX2->Draw("same");
         leg->AddEntry(ProjX1, ObjName1.c_str(), "L");
         leg->AddEntry(ProjX2, ObjName2.c_str(), "L");
         DrawPreliminary("", 13, "");
         SaveCanvas(c1, SaveDir, "Comparison"+SaveName+"_"+ObjName1+"_"+ObjName2+"_ProjX"+string(Id), true);
         delete leg;
         delete c1;
 
         c1  = new TCanvas ("c1", "c1", 600, 600);
         leg = new TLegend (0.50, 0.75, 0.80, 0.90);
         c1->SetLogy (true);
         leg->SetHeader (("Module No. " + string(Id)).c_str());
         leg->SetHeader (SaveName.c_str());
         leg->SetFillColor(0);
         leg->SetFillStyle(0);
         leg->SetBorderSize(0);
         ProjY1->SetStats(kFALSE);
         ProjY1->SetLineColor (kBlack);
         ProjY2->SetLineColor (kBlue);
         ProjY1->GetXaxis()->SetTitle("normalized y coordinate");
         ProjY1->GetYaxis()->SetTitle("number of hits");
         ProjY1->SetAxisRange (-1.5, 1.5, "X");
         ProjY1->Draw("L");
         ProjY2->Draw("same");
         leg->AddEntry(ProjY1, ObjName1.c_str(), "L");
         leg->AddEntry(ProjY2, ObjName2.c_str(), "L");
         DrawPreliminary("", 13, "");
         SaveCanvas(c1, SaveDir, "Comparison"+SaveName+"_"+ObjName1+"_"+ObjName2+"_ProjY"+string(Id), true);
         delete leg;
         delete c1;
         delete ProjX1;
         delete ProjX2;
         delete ProjY1;
         delete ProjY2;
      }
   }
}

void MakeMapPlots(TH3F* Charge_Vs_Path3D, string ObjName, string SaveDir, string Prefix)
{
   for(int x=0;x<17;x++){
      char xProjName[255];
      if(x==0){
         sprintf(xProjName,"%s","SO_inc");
         Charge_Vs_Path3D->GetXaxis()->SetRange(1,14);
      }else if (x==16){
         sprintf(xProjName,"%s","SP_inc");
         Charge_Vs_Path3D->GetXaxis()->SetRange(1,15);
      }else if (x==15){
         sprintf(xProjName,"%s", "PO");
         Charge_Vs_Path3D->GetXaxis()->SetRange(x,x);
      }else{
         sprintf(xProjName,"%02i",x);
         Charge_Vs_Path3D->GetXaxis()->SetRange(x,x);
      }
      printf("---------------\n%s------------\n",xProjName);
      string xProjNameStr(xProjName);


      TH2D*  Charge_Vs_Path2D = (TH2D*)Charge_Vs_Path3D->Project3D("zy");
      char legEntry[128];
      double binMinA = Charge_Vs_Path2D->GetXaxis()->GetBinLowEdge(4);
      double binMaxA = Charge_Vs_Path2D->GetXaxis()->GetBinUpEdge(6);

      TH1D*  Charge_Vs_PathA  = (TH1D*)Charge_Vs_Path2D->ProjectionY("projA",4,6);
      Charge_Vs_PathA->Rebin(2);
      sprintf(legEntry,"[%5.2f,%5.2f]",binMinA,binMaxA); string ALegend (legEntry);

      double binMinB = Charge_Vs_Path2D->GetXaxis()->GetBinLowEdge(7);
      double binMaxB = Charge_Vs_Path2D->GetXaxis()->GetBinUpEdge(9);
      TH1D*  Charge_Vs_PathB  = (TH1D*)Charge_Vs_Path2D->ProjectionY("projB",7,9);
      Charge_Vs_PathB->Rebin(2);
      sprintf(legEntry,"[%5.2f,%5.2f]",binMinB,binMaxB); string BLegend (legEntry);

      double binMinC = Charge_Vs_Path2D->GetXaxis()->GetBinLowEdge(10);
      double binMaxC = Charge_Vs_Path2D->GetXaxis()->GetBinUpEdge(12);
      TH1D*  Charge_Vs_PathC  = (TH1D*)Charge_Vs_Path2D->ProjectionY("projC",10,12);
      Charge_Vs_PathC->Rebin(2);
      sprintf(legEntry,"[%5.2f,%5.2f]",binMinC,binMaxC); string CLegend (legEntry);

      double binMinD = Charge_Vs_Path2D->GetXaxis()->GetBinLowEdge(13);
      double binMaxD = Charge_Vs_Path2D->GetXaxis()->GetBinUpEdge(15);
      TH1D*  Charge_Vs_PathD  = (TH1D*)Charge_Vs_Path2D->ProjectionY("projD",13,15);
      Charge_Vs_PathD->Rebin(2);
      sprintf(legEntry,"[%5.2f,%5.2f]",binMinD,binMaxD); string DLegend (legEntry);

      printf("%f to %f\n",binMinA,binMaxA);
      printf("%f to %f\n",binMinB,binMaxB);
      printf("%f to %f\n",binMinC,binMaxC);
      printf("%f to %f\n",binMinD,binMaxD);


      TCanvas* c0;
      TH1D** Histos = new TH1D* [4]; 
      vector <string> legend;

      c0  = new TCanvas("c0", "c0", 600,600);
      Charge_Vs_Path2D->SetTitle("");
      Charge_Vs_Path2D->SetStats(kFALSE);
      Charge_Vs_Path2D->GetXaxis()->SetTitle("pathlength (mm)");
      Charge_Vs_Path2D->GetYaxis()->SetTitle("#Delta E/#Delta x (ADC/mm)");
      Charge_Vs_Path2D->GetYaxis()->SetTitleOffset(1.80);
      Charge_Vs_Path2D->GetZaxis()->SetRangeUser(1.0,4e+4);
      Charge_Vs_Path2D->Draw("COLZ");

      c0->SetLogz(true);
      SaveCanvas(c0, SaveDir, Prefix + "_" + ObjName+xProjName + "_TH2", true);
      delete c0;


      //Compute Probability Map.
      TH2D* Prob_ChargePath  = new TH2D ("Prob_ChargePath"     , "Prob_ChargePath" , Charge_Vs_Path2D->GetXaxis()->GetNbins(), Charge_Vs_Path2D->GetXaxis()->GetXmin(), Charge_Vs_Path2D->GetXaxis()->GetXmax(), Charge_Vs_Path2D->GetYaxis()->GetNbins(), Charge_Vs_Path2D->GetYaxis()->GetXmin(), Charge_Vs_Path2D->GetYaxis()->GetXmax());
      for(int j=0;j<=Prob_ChargePath->GetXaxis()->GetNbins()+1;j++){
         double Ni = 0;
         for(int k=0;k<=Prob_ChargePath->GetYaxis()->GetNbins()+1;k++){ Ni+=Charge_Vs_Path2D->GetBinContent(j,k);} 

         for(int k=0;k<=Prob_ChargePath->GetYaxis()->GetNbins()+1;k++){
            double tmp = 1E-10;
	         for(int l=0;l<=k;l++){ tmp+=Charge_Vs_Path2D->GetBinContent(j,l);}

            if(Ni>0){
               Prob_ChargePath->SetBinContent (j, k, tmp/Ni);
            }else{
               Prob_ChargePath->SetBinContent (j, k, 0);
            }
         }
      }

      c0  = new TCanvas("c0", "c0", 600,600);
      Prob_ChargePath->SetTitle("Probability MIP(#DeltaE/#DeltaX) < Obs (#DeltaE/#DeltaX)");
      Prob_ChargePath->SetStats(kFALSE);
      Prob_ChargePath->GetXaxis()->SetTitle("pathlength (mm)");
      Prob_ChargePath->GetYaxis()->SetTitle("Observed #DeltaE/#DeltaX (ADC/mm)");
      Prob_ChargePath->GetYaxis()->SetTitleOffset(1.80);
      Prob_ChargePath->GetXaxis()->SetRangeUser(0.28,1.2);
      Prob_ChargePath->GetYaxis()->SetRangeUser(0.0, 1000);
//      Prob_ChargePath->GetYaxis()->SetRangeUser(0,1000);
      Prob_ChargePath->Draw("COLZ");

      //c0->SetLogz(true);
      SaveCanvas(c0, SaveDir, Prefix + "_" + ObjName+xProjName + "_TH2Proba", true);
      delete c0;

      c0 = new TCanvas("c1","c1,",600,600);          legend.clear();
      Histos[0] = Charge_Vs_PathA;                   legend.push_back(ALegend);
      Histos[1] = Charge_Vs_PathB;                   legend.push_back(BLegend);
      Histos[2] = Charge_Vs_PathC;                   legend.push_back(CLegend);
      Histos[3] = Charge_Vs_PathD;                   legend.push_back(DLegend);
      if(Histos[0]->Integral()>=1)Histos[0]->Scale(1/Histos[0]->Integral());
      if(Histos[1]->Integral()>=1)Histos[1]->Scale(1/Histos[1]->Integral());
      if(Histos[2]->Integral()>=1)Histos[2]->Scale(1/Histos[2]->Integral());
      if(Histos[3]->Integral()>=1)Histos[3]->Scale(1/Histos[3]->Integral());
   //   DrawSuperposedHistos((TH1D**)Histos, legend, "",  "Normalized Cluster Charge (ADC/mm)", "u.a.", 0,1200, 0,0);
      DrawSuperposedHistos((TH1**)Histos, legend, "",  "Normalized Cluster Charge (ADC/mm)", "u.a.", 0,600, 0,0);
      DrawLegend((TObject**)Histos,legend,"PathLength (mm):","L");
      c0->SetGridx(true);
      Charge_Vs_PathA->GetXaxis()->SetNdivisions(520);
      SaveCanvas(c0, SaveDir, Prefix+"_"+ObjName+xProjName+"_TH1Linear");
      delete c0;


      c0 = new TCanvas("c1","c1,",600,600);          legend.clear();
      Histos[0] = Charge_Vs_PathA;                   legend.push_back(ALegend);
      Histos[1] = Charge_Vs_PathB;                   legend.push_back(BLegend);
      Histos[2] = Charge_Vs_PathC;                   legend.push_back(CLegend);
      Histos[3] = Charge_Vs_PathD;                   legend.push_back(DLegend);
      if(Histos[0]->Integral()>=1)Histos[0]->Scale(1.0/Histos[0]->Integral());
      if(Histos[1]->Integral()>=1)Histos[1]->Scale(1.0/Histos[1]->Integral());
      if(Histos[2]->Integral()>=1)Histos[2]->Scale(1.0/Histos[2]->Integral());
      if(Histos[3]->Integral()>=1)Histos[3]->Scale(1.0/Histos[3]->Integral());
      DrawSuperposedHistos((TH1**)Histos, legend, "",  "Normalized Cluster Charge (ADC/mm)", "u.a.", 0,3000, 0,0);
//      DrawLegend((TObject**)Histos,legend,"PathLength (mm):","L");
      c0->SetLogy(true);
      SaveCanvas(c0, SaveDir, Prefix + "_"+ObjName+xProjName+"_TH1");
      delete c0;
      delete Charge_Vs_Path2D;
      delete Charge_Vs_PathA;
      delete Charge_Vs_PathB;
      delete Charge_Vs_PathC;
      delete Charge_Vs_PathD;
      delete Prob_ChargePath;
   }
}

void MakeROCGeneral (TFile* InputFile1, TFile* InputFile2, vector<string> HistoNames, vector<string> LegendLabels, vector<Color_t> Colors, string SaveDir, string suffix, bool withErrorBars, bool Every2ndIsDashed) {

      TCanvas* c1   = new TCanvas ("c1", "c1", 600,600); 
      TLegend* leg  = new TLegend (0.30, 0.15, 0.80, 0.15+0.05*HistoNames.size()/(Every2ndIsDashed?2:1));
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      c1->SetLogx(true);
      TH1D h ("tmp", "tmp", 1, 8E-7, 1);
      h.GetXaxis()->SetTitle("background efficiency");
      h.GetXaxis()->SetNdivisions(5);
      h.GetYaxis()->SetTitle("signal efficiency");
      h.GetYaxis()->SetNdivisions(5);
      h.SetAxisRange (0.65,1.0,"Y");
      h.SetStats(0);
      h.Draw();
      TGraphErrors** ROC = new TGraphErrors* [HistoNames.size()];
      for (size_t NameIndex = 0; NameIndex < HistoNames.size(); NameIndex++)
      {
         int divide = 1;
         TH1D* HdedxMIP1 = (TH1D*) GetObjectFromPath(InputFile1, (HistoNames[NameIndex]).c_str() );
         TH1D* HdedxMIP2 = (TH1D*) GetObjectFromPath(InputFile2, (HistoNames[NameIndex]).c_str() );
         ROC[NameIndex]  = new TGraphErrors(HdedxMIP1->GetNbinsX()/divide + 1);

         double fullBkg  = HdedxMIP1->Integral(0, HdedxMIP1->GetNbinsX()+1),
                fullSig  = HdedxMIP2->Integral(0, HdedxMIP2->GetNbinsX()+1);
         for (unsigned int cut_i = 1; cut_i <= HdedxMIP1->GetNbinsX()/divide; cut_i++){
            double a = HdedxMIP2->Integral(0, cut_i*divide),
                   e = (fullSig-a)/fullSig;
            ROC[NameIndex]->SetPoint (cut_i-1, 1 - HdedxMIP1->Integral(0, cut_i*divide)/fullBkg, 1 - a/fullSig);
            ROC[NameIndex]->SetPointError (cut_i-1, 0, withErrorBars?sqrt(e*(1-e)/(fullSig-a)):0);
         }
         double a = HdedxMIP2->Integral(0, HdedxMIP2->GetNbinsX()+1),
                e = (fullSig-a)/fullSig;
         ROC[NameIndex]->SetPoint (HdedxMIP1->GetNbinsX(), 1 - HdedxMIP1->Integral(0, HdedxMIP1->GetNbinsX()+1)/fullBkg, 1 - a/fullSig);
         ROC[NameIndex]->SetPointError (HdedxMIP1->GetNbinsX(), 0, withErrorBars?sqrt(e*(1-e)/(fullSig-a)):0);

         ROC[NameIndex]->SetLineColor(Colors[NameIndex]);
         if (Every2ndIsDashed && NameIndex&1) {ROC[NameIndex]->SetLineStyle(3); ROC[NameIndex]->SetLineWidth(1);}
         else  {leg->AddEntry (ROC[NameIndex], LegendLabels[NameIndex].c_str(), "L"); ROC[NameIndex]->SetLineWidth(2);}

         ROC[NameIndex]->Draw("same");
      }
      leg->Draw();
      DrawPreliminary("", 13, "");
      SaveCanvas(c1, SaveDir, "Comparison_ROC_" + string(withErrorBars?"WError_":"") + suffix);
      for (size_t NameIndex = 0; NameIndex < HistoNames.size(); NameIndex++)
         delete ROC[NameIndex];
      delete ROC;
      delete leg;
      delete c1;
}

void CrossCompareAndControlPlots (string SaveDir, vector <dEdxPlotObj*> plotObj, string Reject, string Select){
   // some ROC plots
/*
  for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type != 2) continue;
      vector <string> ObjNames; vector <string> LegendLabels; vector <Color_t> Colors;

      ObjNames.push_back("harm2_SP_in_noC_CCC_MIP");       LegendLabels.push_back("harmonic2, old CC");  Colors.push_back(kRed);
      ObjNames.push_back("harm2_SP_in_noC_newCCC_MIP");    LegendLabels.push_back("harmonic2, new CC");  Colors.push_back(kBlue);
      ObjNames.push_back("Hybr2015_SP_in_noC_CCC_MIP");    LegendLabels.push_back("hybrid2-15, old CC"); Colors.push_back(kOrange);
      ObjNames.push_back("Hybr2015_SP_in_noC_newCCC_MIP"); LegendLabels.push_back("hybrid2-15, new CC"); Colors.push_back(kGreen);
      ObjNames.push_back("Ias_SP_in_noC_CCC_MIP");         LegendLabels.push_back("I_{as}, old CC");     Colors.push_back(kBlack);
      ObjNames.push_back("Ias_SP_in_noC_newACCC_MIP");      LegendLabels.push_back("I_{as}, new CC");     Colors.push_back(kOrange+3);
      MakeROCGeneral (plotObj[2]->InputFile, plotObj[i]->InputFile, ObjNames, LegendLabels, Colors, SaveDir, "Estimators_hybrid-"+plotObj[i]->SavePrefix);
   }

   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type != 2) continue;
      vector <string> ObjNames; vector <string> LegendLabels; vector <Color_t> Colors;
      ObjNames.push_back("harm2_SP_in_noC_CCC_MIP");    LegendLabels.push_back("harmonic2");  Colors.push_back(kRed);
      ObjNames.push_back("Hybr201_SP_in_noC_CCC_MIP");  LegendLabels.push_back("Hybrid2-10"); Colors.push_back(kBlue);
      ObjNames.push_back("Hybr2015_SP_in_noC_CCC_MIP"); LegendLabels.push_back("Hybrid2-15"); Colors.push_back(kOrange);
      ObjNames.push_back("Hybr202_SP_in_noC_CCC_MIP");  LegendLabels.push_back("Hybrid2-20"); Colors.push_back(kGreen);
      ObjNames.push_back("Hybr2025_SP_in_noC_CCC_MIP"); LegendLabels.push_back("Hybrid2-25"); Colors.push_back(kMagenta);
      MakeROCGeneral (plotObj[2]->InputFile, plotObj[i]->InputFile, ObjNames, LegendLabels, Colors, SaveDir, "Estimators_Hybrid"+plotObj[i]->SavePrefix);
   }

   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type != 2) continue;
      vector <string> ObjNames; vector <string> LegendLabels; vector <Color_t> Colors;
      ObjNames.push_back("harm2_SP_in_noC_CCC_noF_MIP");     LegendLabels.push_back("harmonic2");  Colors.push_back(kRed);
      ObjNames.push_back("Hybr201_SP_in_noC_CCC_noF_MIP");   LegendLabels.push_back("Hybrid2-10"); Colors.push_back(kBlue);
      ObjNames.push_back("Hybr2015_SP_in_noC_CCC_noF_MIP");  LegendLabels.push_back("Hybrid2-15"); Colors.push_back(kOrange);
      ObjNames.push_back("Hybr202_SP_in_noC_CCC_noF_MIP");   LegendLabels.push_back("Hybrid2-20"); Colors.push_back(kGreen);
      ObjNames.push_back("Hybr2025_SP_in_noC_CCC_noF_MIP");  LegendLabels.push_back("Hybrid2-25"); Colors.push_back(kMagenta);
      MakeROCGeneral (plotObj[2]->InputFile, plotObj[i]->InputFile, ObjNames, LegendLabels, Colors, SaveDir, "Estimators_Hybrid_noF"+plotObj[i]->SavePrefix);
   }
*/
   for (size_t i = 0; i < plotObj.size(); i++){
      TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
      c1->SetLogy(true);

      vector<TH1D*> histos;
      vector <string> legend;
      double max = 0;
      for (size_t j = 0; j < plotObj[i]->StdObjName.size()-1; j++){
         if (plotObj[i]->StdObjName[j].find(Reject)!=string::npos) continue;
         if (plotObj[i]->StdObjName[j].find(Select)==string::npos && plotObj[i]->StdObjName[j].find("harm")==string::npos) continue;
         histos.push_back(plotObj[i]->HdedxMIP[j]);
          histos[histos.size()-1]->Scale(1.0/histos[histos.size()-1]->Integral());
         legend.push_back(plotObj[i]->StdObjLegend[j]);

         max = (max<histos[histos.size()-1]->GetMaximum())?histos[histos.size()-1]->GetMaximum():max;
      }
      DrawSuperposedHistos((TH1**) (&histos[0]), legend, "L", "dE/dx (MeV/cm)", "arbitrary units", 0, (plotObj[i]->type==2)?15:8, 5e-7, max*12, true);
      DrawLegend ((TObject**) (&histos[0]), legend, "Estimators", "L", 0.8, 0.9, 0.3, 0.05);
      DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
      SaveCanvas (c1, SaveDir, "Comp_Estim_" + Select + "_" + (Reject!=""?string("No_")+Reject+"_":"") + plotObj[i]->SavePrefix + "_MIP");
      delete c1;
   }


   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type == 2) continue;
      TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
      c1->SetLogy(true);

      vector<TH1D*> histos;
      vector <string> legend;
      for (size_t j = 0; j < plotObj[i]->StdObjName.size()-1; j++){
         if (plotObj[i]->StdObjName[j].find(Reject)!=string::npos) continue;
         if (plotObj[i]->StdObjName[j].find(Select)==string::npos && plotObj[i]->StdObjName[j].find("harm")==string::npos) continue;
         histos.push_back(plotObj[i]->HProtonHitSO[j]);
         legend.push_back(plotObj[i]->StdObjLegend[j]);
      }
      histos[0]->Scale(1.0/histos[0]->Integral());
      DrawSuperposedHistos((TH1**) (&histos[0]), legend, "L", "dE/dx (MeV/cm)", "arbitrary units", 0, 20, 5e-7, histos[0]->GetMaximum()*1.5, true);
      DrawLegend ((TObject**) (&histos[0]), legend, "Estimators", "L", 0.8, 0.9, 0.3, 0.05);
      DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
      SaveCanvas (c1, SaveDir, "Comp_Estim_" + Select + "_" + (Reject!=""?string("No_")+Reject+"_":"") + plotObj[i]->SavePrefix + "_ProtonHitSO");
      delete c1;
   }

   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type == 2) continue;
      TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
      c1->SetLogy(true);

      vector<TH1D*> histos;
      vector <string> legend;
      for (size_t j = 0; j < plotObj[i]->StdObjName.size()-1; j++){
         if (plotObj[i]->StdObjName[j].find(Reject)!=string::npos) continue;
         if (plotObj[i]->StdObjName[j].find(Select)==string::npos && plotObj[i]->StdObjName[j].find("harm")==string::npos) continue;
         histos.push_back(plotObj[i]->HProtonHitPO[j]);
         legend.push_back(plotObj[i]->StdObjLegend[j]);
      }
      histos[0]->Scale(1.0/histos[0]->Integral());
      DrawSuperposedHistos((TH1**) (&histos[0]), legend, "L", "dE/dx (MeV/cm)", "arbitrary units", 0, 20, 5e-7, histos[0]->GetMaximum()*1.5, true);
      DrawLegend ((TObject**) (&histos[0]), legend, "Estimators", "L", 0.8, 0.9, 0.3, 0.05);
      DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
      SaveCanvas (c1, SaveDir, "Comp_Estim_" + Select + "_" + (Reject!=""?string("No_")+Reject+"_":"") + plotObj[i]->SavePrefix + "_ProtonHitPO");
      delete c1;
   }

   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type != 2) continue;
      TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
      c1->SetLogy(true);
      c1->SetGridx(true);

      vector<TH1D*> histos;
      vector <string> legend;
      double max = 0, min = 9999;
      for (size_t j = 0; j < plotObj[i]->StdObjName.size()-1; j++){
         if (plotObj[i]->StdObjName[j].find(Reject)!=string::npos) continue;
         if (plotObj[i]->StdObjName[j].find(Select)==string::npos && plotObj[i]->StdObjName[j].find("harm")==string::npos) continue;
         histos.push_back(plotObj[i]->HMassHSCP[j]);
         legend.push_back(plotObj[i]->StdObjLegend[j]);
         histos[histos.size()-1]->Scale(1.0/histos[histos.size()-1]->Integral());
         max = histos[histos.size()-1]->GetMaximum()>max?histos[histos.size()-1]->GetMaximum():max;
         min = histos[histos.size()-1]->GetMinimum()>min?histos[histos.size()-1]->GetMinimum():min;
      }
      DrawSuperposedHistos((TH1**) (&histos[0]), legend, "hist", "Mass (GeV)", "arbitrary units", 0, 3000, min, max*12, true);
      DrawLegend ((TObject**) (&histos[0]), legend, "Estimators", "L", 0.8, 0.9, 0.3, 0.05);
      DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
      SaveCanvas (c1, SaveDir, "Comp_Estim_" + Select + "_" + (Reject!=""?string("No_")+Reject+"_":"") + plotObj[i]->SavePrefix + "_HSCPMassOld");
      delete c1;
   }

   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type==2) continue;
      TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
      TLegend* leg = new TLegend(0.50, 0.80, 0.80, 0.90);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);

      vector <string> legend;
      vector <Color_t> colors;
      TGraphErrors** g = new TGraphErrors* [2];
      vector <double> x;
      vector <double> xErr;
      vector <double> y0;
      vector <double> y0Err;
      vector <double> y1;
      vector <double> y1Err;
      for (size_t j = 0; j < plotObj[i]->StdObjName.size()-1; j++){
         if (plotObj[i]->StdObjName[j].find(Reject)==string::npos && plotObj[i]->StdObjName[j].find(Select)!=string::npos){
            if (plotObj[i]->StdObjName[j].find("harm")!=string::npos) x.push_back(0.0);
            else{
               size_t zeroPos = plotObj[i]->StdObjName[j].find("0") + 1;
               size_t endPos  = plotObj[i]->StdObjName[j].find("_") - zeroPos;
               x.push_back(atof(("0."+plotObj[i]->StdObjName[j].substr(zeroPos, endPos)).c_str()));
            }
            y0.push_back(plotObj[i]->HProtonHitSO[j]->Integral(0, plotObj[i]->HProtonHitSO[j]->FindBin(2.0))/
                  plotObj[i]->HProtonHitSO[j]->Integral(0, plotObj[i]->HProtonHitSO[j]->GetNbinsX()+1));
            y1.push_back(plotObj[i]->HProtonHitPO[j]->Integral(0, plotObj[i]->HProtonHitPO[j]->FindBin(2.0))/
                  plotObj[i]->HProtonHitPO[j]->Integral(0, plotObj[i]->HProtonHitPO[j]->GetNbinsX()+1));

            xErr.push_back(0);
         }
      }

      double MaxY = 0;
      double MinY = 1e+300;
      for (size_t j = 0; j < x.size(); j++){
         MaxY = (MaxY>std::max(y0[j],y1[j]))?MaxY:std::max(y0[j],y1[j]);
         MinY = (MinY<std::min(y0[j],y1[j]))?MinY:std::min(y0[j],y1[j]);

         y0Err.push_back(0/*sqrt(y0[j])*/);
         y1Err.push_back(0/*sqrt(y1[j])*/);
            
         cerr << "y0[" << j << "] = " << y0[j] << "\t" << "y1[" << j << "] = " << y1[j] << endl;
      }

      g[0] = new TGraphErrors (x.size(), &x[0], &y0[0], &xErr[0], &y0Err[0]);
      g[1] = new TGraphErrors (x.size(), &x[0], &y1[0], &xErr[0], &y1Err[0]);

      colors.push_back (kRed);
      colors.push_back (kBlue);

      legend.push_back ("Strip Charges");
      legend.push_back ("Pixel Charges");

      if (MinY == 0) MinY = 5e-7;

      TH1D h ("temp", "temp", 1, x[0], x[x.size()-1]*1.1);
      h.GetXaxis()->SetTitle("hybrid index");
      h.GetYaxis()->SetTitle("fraction of clusters");
      h.GetXaxis()->SetRangeUser (x[0], x[x.size()-1]*1.1);
      h.SetAxisRange (MinY, MaxY*1.2, "Y");
      h.GetYaxis()->SetTitleOffset(1.35);
      h.SetStats(0);
      h.Draw();

      for (unsigned int k = 0; k < 2; k++){
         g[k]->SetMarkerColor (colors[k]);
         g[k]->SetLineColor (colors[k]);
         g[k]->SetMarkerStyle (20);
         g[k]->Draw("same EP");
         leg->AddEntry (g[k], legend[k].c_str(), "LP");
      }
      leg->Draw();
      DrawPreliminary ("", 13.0, plotObj[i]->SavePrefix);
      SaveCanvas (c1, SaveDir, "Comp_Estim_" + Select + "_"+ (Reject!=""?string("No_")+Reject+"_":"") + plotObj[i]->SavePrefix + "_LowHits");
      delete c1;
      delete leg;
   }

   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type != 2) continue;
      TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
      c1->SetLogy(true);
      c1->SetGridx(true);

      vector<TH1D*> histos;
      vector <string> legend;
      double max = 0;
      for (size_t j = 0; j < plotObj[i]->StdObjName.size()-1; j++){
         if (plotObj[i]->StdObjName[j].find(Reject)!=string::npos) continue;
         if (plotObj[i]->StdObjName[j].find(Select)==string::npos && plotObj[i]->StdObjName[j].find("harm")==string::npos) continue;
         histos.push_back(plotObj[i]->HMassHSCP[j]);
         histos[histos.size()-1]->Reset();
         for (int x = 1; x <= plotObj[i]->HdedxVsP[j]->GetNbinsX(); x++){
            for (int y = 1; y <= plotObj[i]->HdedxVsP[j]->GetNbinsY(); y++){
               if(plotObj[i]->HdedxVsP[j]->GetYaxis()->GetBinCenter(y)<5.0)continue;
               histos[histos.size()-1]->Fill(GetMass (plotObj[i]->HdedxVsP[j]->GetXaxis()->GetBinCenter(x),
                                        plotObj[i]->HdedxVsP[j]->GetYaxis()->GetBinCenter(y),
                                        plotObj[i], plotObj[i]->StdObjName[j]),
                                        plotObj[i]->HdedxVsP[j]->GetBinContent(x,y));
            }
         }
         histos[histos.size()-1]->Scale(1.0/histos[histos.size()-1]->Integral());
         if (max < histos[histos.size()-1]->GetMaximum()) max = histos[histos.size()-1]->GetMaximum();
         legend.push_back(plotObj[i]->StdObjLegend[j]);
      }

      DrawSuperposedHistos((TH1**) (&histos[0]), legend, "hist", "HSCP Mass (GeV)", "arbitrary units", 0, 2500, 5e-5, max*12.0, true);
      DrawLegend ((TObject**) (&histos[0]), legend, "Estimators", "L", 0.8, 0.9, 0.3, 0.05);
      DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
      SaveCanvas (c1, SaveDir, "Comp_Estim_" + Select + "_" + (Reject!=""?string("No_")+Reject+"_":"") + plotObj[i]->SavePrefix + "_HSCPMass");
      delete c1;
   }
}

void SuperposeFilesOnDeDxObj (string SaveDir, vector<dEdxPlotObj*> plotObj){
   for (size_t j=0; j < plotObj[0]->StdObjName.size(); j++){
      TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
      vector<string> legend;
      vector<TH1D*> histos;

      bool isEstim = (plotObj[0]->StdObjName[j].find("Ias")==string::npos);
      for (size_t i=0; i < plotObj.size(); i++){
         histos.push_back(plotObj[i]->HdedxVsPProfile[j]);
         legend.push_back (plotObj[i]->LegEntry);
      }
      DrawSuperposedHistos((TH1**) &(histos[0]), legend, "E1", "track momentum (GeV)", isEstim?"dE/dx (MeV/cm)":"I_{as}", 0, 5,
            isEstim?2.5:0, isEstim?5:1);
      DrawLegend ((TObject**) &(histos[0]), legend, plotObj[0]->StdObjLegend[j], "LP", 0.80, 0.90, 0.40, 0.05);
      DrawPreliminary("", 13.0, "");
      SaveCanvas(c1, SaveDir, plotObj[0]->StdObjName[j] + "_Comp_Profile");
      histos.clear();
      delete c1;


      c1 = new TCanvas("c1", "c1", 600,600);
      legend.clear();
      for (size_t i=0; i < plotObj.size(); i++){
         histos.push_back(plotObj[i]->HdedxVsEtaProfile[j]);
         legend.push_back (plotObj[i]->LegEntry);
      }
      DrawSuperposedHistos((TH1**) &(histos[0]), legend, "E1", "#eta", isEstim?"dE/dx (MeV/cm)":"I_{as}", -2.1, 2.1,
            isEstim?2.5:0, isEstim?5:1);
      DrawLegend ((TObject**) &(histos[0]), legend, plotObj[0]->StdObjLegend[j], "LP", 0.80, 0.90, 0.40, 0.05);
      DrawPreliminary("", 13.0, "");
      SaveCanvas(c1, SaveDir, plotObj[0]->StdObjName[j] + "_Comp_HdedxVsEtaProfile");
      histos.clear();
      delete c1;


      c1 = new TCanvas("c1", "c1", 600,600);
      legend.clear();
      c1->SetLogy(true);
      for (size_t i=0; i < plotObj.size(); i++){
	 if (plotObj[i]->type == 2) continue;
         histos.push_back(plotObj[i]->HdedxMIP[j]);
         legend.push_back (plotObj[i]->LegEntry);
      }
      DrawSuperposedHistos((TH1**) &(histos[0]), legend, "L", isEstim?"dE/dx (MeV/cm)":"I_{as}", "arbitrary units",
            0, isEstim?7:1, 5e-7, 6, true);
      DrawLegend ((TObject**) &(histos[0]), legend, plotObj[0]->StdObjLegend[j], "LP", 0.80, 0.90, 0.40, 0.05);
      DrawPreliminary("", 13.0, "");
      SaveCanvas(c1, SaveDir, plotObj[0]->StdObjName[j] + "_Comp_MIP");
      histos.clear();
      delete c1;
/*
      c1 = new TCanvas("c1", "c1", 600,600);
      legend.clear();
      c1->SetLogy(true);
      for (size_t i=0; i < plotObj.size(); i++){
	 if (plotObj[i]->type != 2) continue;
         histos.push_back(plotObj[i]->HdedxMIP[j]);
         legend.push_back (plotObj[i]->LegEntry);
      }
      DrawSuperposedHistos((TH1**) &(histos[0]), legend, "L", isEstim?"dE/dx (MeV/cm)":"I_{as}", "arbitrary units",
            0, isEstim?7:1, 5e-7, 6, true);
      DrawLegend ((TObject**) &(histos[0]), legend, plotObj[0]->StdObjLegend[j], "LP", 0.80, 0.90, 0.40, 0.05);
      DrawPreliminary("", 13.0, "");
      SaveCanvas(c1, SaveDir, plotObj[0]->StdObjName[j] + "_Comp_HSCPMIP");
      histos.clear();
      delete c1;
*/
      if (isEstim){
         c1 = new TCanvas("c1", "c1", 600,600);
         legend.clear();
         c1->SetLogy(true);
         double min = 1005, max = 0;
         for (size_t i=0; i < plotObj.size(); i++){
            if (plotObj[i]->type==2) continue;
            histos.push_back(plotObj[i]->HMass[j]);
            histos[histos.size()-1]->Reset();
            for (int x = 1; x <= plotObj[i]->HdedxVsP[j]->GetNbinsX(); x++){
	      if (plotObj[i]->HdedxVsP[j]->GetXaxis()->GetBinCenter(x) > 3.0) continue; //3.0 was the default
               for (int y = 1; y <= plotObj[i]->HdedxVsP[j]->GetNbinsY(); y++){
		 if(plotObj[i]->HdedxVsP[j]->GetYaxis()->GetBinCenter(y)<plotObj[i]->C[plotObj[i]->StdObjName[j]]+1.5)continue; //3.0 was the default
                  histos[histos.size()-1]->Fill(GetMass (plotObj[i]->HdedxVsP[j]->GetXaxis()->GetBinCenter(x),
                                           plotObj[i]->HdedxVsP[j]->GetYaxis()->GetBinCenter(y),
                                           plotObj[i], plotObj[i]->StdObjName[j]),
                                           plotObj[i]->HdedxVsP[j]->GetBinContent(x,y));
               }
            }
            histos[histos.size()-1]->Scale(1.0/histos[histos.size()-1]->Integral());
            if (max < histos[histos.size()-1]->GetMaximum()) max = histos[histos.size()-1]->GetMaximum();
            if (min > histos[histos.size()-1]->GetMinimum()) min = histos[histos.size()-1]->GetMinimum();
            legend.push_back (plotObj[i]->LegEntry);
         }

         if (min == 0) min = 1e-7;
         DrawSuperposedHistos((TH1**) &(histos[0]), legend, "L", "Mass (GeV)", "arbitrary units",
               0, 5, min, max*12.0);
         DrawLegend ((TObject**) &(histos[0]), legend, plotObj[0]->StdObjLegend[j], "P", 0.80, 0.90, 0.30, 0.05);

         TLine* lineKaon = new TLine(0.493667, min, 0.493667, max);
         lineKaon->SetLineWidth(2);
         lineKaon->SetLineStyle(2);
         lineKaon->SetLineColor(9);
         TLine* lineProton = new TLine(0.938272, min, 0.938272, max);
         lineProton->SetLineWidth(2);
         lineProton->SetLineStyle(2);
         lineProton->SetLineColor(9);
         TLine* lineDeuteron = new TLine(1.88, min, 1.88, max);
         lineDeuteron->SetLineWidth(2);
         lineDeuteron->SetLineStyle(2);
         lineDeuteron->SetLineColor(9);

         lineKaon->Draw("same");
         lineProton->Draw("same");
         lineDeuteron->Draw("same");

         DrawPreliminary("", 13.0, "");
         SaveCanvas(c1, SaveDir, plotObj[0]->StdObjName[j] + "_Comp_Mass");
         histos.clear();
         delete c1;
         delete lineKaon;
         delete lineProton;
         delete lineDeuteron;
      }


      // SPECIAL PLOT REQUESTED BY LOIC -- CAN DELETE AFTERWARDS
      if (isEstim){
         c1 = new TCanvas("c1", "c1", 600,600);
         legend.clear();
         c1->SetLogy(true);
         double min = 1005, max = 0;
         double dEdxK_Data2015 = 2.684; // 2015 Data&MC K
         double dEdxC_Data2015 = 3.375; // 2015 Data&MC C
         double dEdxK_Data1 = 2.580; // 2016 Data K
         double dEdxC_Data1 = 3.922; // 2016 Data C
         double dEdxK_MC1   = 2.935; // 2016 MC K
         double dEdxC_MC1   = 3.197; // 2016 MC C
         for (size_t i=0; i < plotObj.size(); i++){
	    bool is2016 = (plotObj[i]->FileName.find("2016")!=std::string::npos);
            if (plotObj[i]->type==2) continue;
            double dEdxK = plotObj[i]->type==0?dEdxK_Data1:dEdxK_MC1;
            double dEdxC = plotObj[i]->type==0?dEdxC_Data1:dEdxC_MC1;
	    if (plotObj[i]->FileName.find("Data")!=std::string::npos){
		   dEdxK = is2016?dEdxK_Data1:dEdxK_Data2015;
		   dEdxC = is2016?dEdxC_Data1:dEdxC_Data2015;
	    }
            histos.push_back(plotObj[i]->HMass[j]);
            histos[histos.size()-1]->Reset();
            for (int x = 1; x <= plotObj[i]->HdedxVsP[j]->GetNbinsX(); x++){
               if (plotObj[i]->HdedxVsP[j]->GetXaxis()->GetBinCenter(x) > 3.0) continue;
               for (int y = 1; y <= plotObj[i]->HdedxVsP[j]->GetNbinsY(); y++){
                  if(plotObj[i]->HdedxVsP[j]->GetYaxis()->GetBinCenter(y)<plotObj[i]->C[plotObj[i]->StdObjName[j]]+(is2016?5.0:3.0))continue;
                  histos[histos.size()-1]->Fill(GetMass (plotObj[i]->HdedxVsP[j]->GetXaxis()->GetBinCenter(x),
                                           plotObj[i]->HdedxVsP[j]->GetYaxis()->GetBinCenter(y),
                                           &dEdxK, &dEdxC),
                                           plotObj[i]->HdedxVsP[j]->GetBinContent(x,y));
               }
            }
            histos[histos.size()-1]->Scale(1.0/histos[histos.size()-1]->Integral());
            if (max < histos[histos.size()-1]->GetMaximum()) max = histos[histos.size()-1]->GetMaximum();
            if (min > histos[histos.size()-1]->GetMinimum()) min = histos[histos.size()-1]->GetMinimum();
            legend.push_back (plotObj[i]->LegEntry);
         }

         if (min == 0) min = 1e-7;
         DrawSuperposedHistos((TH1**) &(histos[0]), legend, "L", "Mass (GeV)", "arbitrary units",
               0, 5, min, max*12.0);
         DrawLegend ((TObject**) &(histos[0]), legend, plotObj[0]->StdObjLegend[j], "P", 0.80, 0.90, 0.30, 0.05);

         TLine* lineKaon = new TLine(0.493667, min, 0.493667, max);
         lineKaon->SetLineWidth(2);
         lineKaon->SetLineStyle(2);
         lineKaon->SetLineColor(9);
         TLine* lineProton = new TLine(0.938272, min, 0.938272, max);
         lineProton->SetLineWidth(2);
         lineProton->SetLineStyle(2);
         lineProton->SetLineColor(9);
         TLine* lineDeuteron = new TLine(1.88, min, 1.88, max);
         lineDeuteron->SetLineWidth(2);
         lineDeuteron->SetLineStyle(2);
         lineDeuteron->SetLineColor(9);

         lineKaon->Draw("same");
         lineProton->Draw("same");
         lineDeuteron->Draw("same");

         DrawPreliminary("", 13.0, "");
         SaveCanvas(c1, SaveDir, plotObj[0]->StdObjName[j] + "_Comp_Mass_KCData");
         histos.clear();
         delete c1;
         delete lineKaon;
         delete lineProton;
         delete lineDeuteron;
      }
   }
}

void Draw2D (string SaveDir, vector<dEdxPlotObj*> plotObj){
   for (size_t i=0; i < plotObj.size(); i++){
      for (size_t j=0; j < plotObj[i]->StdObjName.size(); j++){
         TCanvas* c1 = new TCanvas("c1", "c1", 600,600);
         c1->SetLogz(true);
         TH2D h2 ("tmp", "tmp", 1, 0, 5, 1, 0, (plotObj[i]->StdObjName[j].find("Ias")==string::npos)?17:1.2);
         h2.SetStats(kFALSE);
         h2.SetStats(kFALSE);
         h2.GetXaxis()->SetTitle("p (GeV)");
         h2.GetYaxis()->SetTitle(plotObj[i]->StdObjName[j].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
         if (plotObj[i]->type != 2){
            h2.SetAxisRange(0,5,"X");
            h2.SetAxisRange(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?17:1.2,"Y");
            h2.GetYaxis()->SetRangeUser(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?17:1.2);
            h2.GetXaxis()->SetRangeUser(0, 5.0);
         } else {
            h2.SetBins(1,0,1200,1,0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?30:1.2);
            h2.SetAxisRange(0,1200,"X");
            h2.SetAxisRange(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?30:1.2,"Y");
            h2.GetYaxis()->SetRangeUser(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?30:1.2);
            h2.GetXaxis()->SetRangeUser(0, 1200.0);
         }
         h2.Draw("COLZ");
         plotObj[i]->HdedxVsP[j]->Draw("same COLZ");

         TF1* PionLine      = NULL;
	 TF1* PionLineFit      = NULL;
         TF1* KaonLine      = NULL;
         TF1* ProtonLine    = NULL;
         TF1* DeuteronLine  = NULL;
         TF1* DeuteronLineFit = NULL;

/*
            plotObj[i]->HdedxVsP[j]->SetStats(kFALSE);
            plotObj[i]->HdedxVsP[j]->GetXaxis()->SetTitle("p (GeV)");
            plotObj[i]->HdedxVsP[j]->GetYaxis()->SetTitle(plotObj[i]->StdObjName[j].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
            plotObj[i]->HdedxVsP[j]->SetAxisRange(0,5,"X");
            plotObj[i]->HdedxVsP[j]->SetAxisRange(0,17,"Y");
            plotObj[i]->HdedxVsP[j]->GetYaxis()->SetRangeUser(0,17);
            plotObj[i]->HdedxVsP[j]->Draw("COLZ");
*/
         if (plotObj[i]->type != 2 && plotObj[i]->StdObjName[j].find("Ias")==string::npos){
            TF1* PionLine = GetMassLine(0.140, plotObj[i], plotObj[i]->StdObjName[j]);
            PionLine->SetLineColor(1);
            PionLine->SetLineWidth(2);
            PionLine->SetRange(PionLine->GetX(15),1 );//PionLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

	    TF1* PionLineFit = GetMassLine(0.140, plotObj[i], plotObj[i]->StdObjName[j]);
            PionLineFit->SetLineColor(8);
	    PionLineFit->SetLineWidth(2);
            //PionLine->SetRange(PionLine->GetX(15), PionLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));
	    PionLineFit->SetRange(1,4);

            TF1* KaonLine = GetMassLine(0.494, plotObj[i], plotObj[i]->StdObjName[j]);
            KaonLine->SetLineColor(1);
            KaonLine->SetLineWidth(2);
            KaonLine->SetRange(KaonLine->GetX(15), KaonLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

            TF1* ProtonLine = GetMassLine(0.938, plotObj[i], plotObj[i]->StdObjName[j]);
            ProtonLine->SetLineColor(1);
            ProtonLine->SetLineWidth(2);
            ProtonLine->SetRange(ProtonLine->GetX(15), ProtonLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

            TF1* DeuteronLine = GetMassLine(1.88, plotObj[i], plotObj[i]->StdObjName[j]);
            DeuteronLine->SetLineColor(1);
            DeuteronLine->SetLineWidth(2);
            DeuteronLine->SetRange(DeuteronLine->GetX(15), DeuteronLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

            TF1* DeuteronLineFit = GetMassLine(1.88, plotObj[i], plotObj[i]->StdObjName[j]);
            DeuteronLineFit->SetLineColor(2);
            DeuteronLineFit->SetLineWidth(2);
            DeuteronLineFit->SetRange(1.0,1.6); //range to fix

            PionLine->Draw("same");
	    PionLineFit->Draw("same");
            KaonLine->Draw("same");
            ProtonLine->Draw("same");
            DeuteronLine->Draw("same");
            DeuteronLineFit->Draw("same");
         }
         DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
         SaveCanvas(c1, SaveDir, plotObj[i]->StdObjName[j] + "_" + plotObj[i]->SavePrefix + "_dedxVsP", true);
         if (PionLine){ delete PionLine; delete KaonLine; delete ProtonLine; delete DeuteronLine; delete DeuteronLineFit; }
         delete c1;

         // REAL K and C :: FOR AN -- LOIC
         c1 = new TCanvas("c1", "c1", 600,600);
         c1->SetLogz(true);
         TH2D h3 ("tmp", "tmp", 1, 0, 5, 1, 0, (plotObj[i]->StdObjName[j].find("Ias")==string::npos)?17:1.2);
         h3.SetStats(kFALSE);
         h3.SetStats(kFALSE);
         h3.GetXaxis()->SetTitle("p (GeV)");
         h3.GetYaxis()->SetTitle(plotObj[i]->StdObjName[j].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
         if (plotObj[i]->type != 2){
            h3.SetAxisRange(0,5,"X");
            h3.SetAxisRange(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?17:1.2,"Y");
            h3.GetYaxis()->SetRangeUser(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?17:1.2);
            h3.GetXaxis()->SetRangeUser(0, 5.0);
         } else {
            h3.SetBins(1,0,1200,1,0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?30:1.2);
            h3.SetAxisRange(0,1200,"X");
            h3.SetAxisRange(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?30:1.2,"Y");
            h3.GetYaxis()->SetRangeUser(0,(plotObj[i]->StdObjName[j].find("Ias")==string::npos)?30:1.2);
            h3.GetXaxis()->SetRangeUser(0, 1200.0);
         }
         h3.Draw("COLZ");
         plotObj[i]->HdedxVsP[j]->Draw("same COLZ");

         PionLine      = NULL;
         KaonLine      = NULL;
         ProtonLine    = NULL;
         DeuteronLine  = NULL;
         DeuteronLineFit = NULL;

/*
            plotObj[i]->HdedxVsP[j]->SetStats(kFALSE);
            plotObj[i]->HdedxVsP[j]->GetXaxis()->SetTitle("p (GeV)");
            plotObj[i]->HdedxVsP[j]->GetYaxis()->SetTitle(plotObj[i]->StdObjName[j].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
            plotObj[i]->HdedxVsP[j]->SetAxisRange(0,5,"X");
            plotObj[i]->HdedxVsP[j]->SetAxisRange(0,17,"Y");
            plotObj[i]->HdedxVsP[j]->GetYaxis()->SetRangeUser(0,17);
            plotObj[i]->HdedxVsP[j]->Draw("COLZ");
*/
         if (plotObj[i]->type != 2 && plotObj[i]->StdObjName[j].find("Ias")==string::npos){
            double dEdxK_Data1 = 2.580; // 2016 Data K
            double dEdxC_Data1 = 3.922; // 2016 Data C
            double dEdxK_MC1   = 2.935; // 2016 MC K
            double dEdxC_MC1   = 3.197; // 2016 MC C
            TF1* PionLine = GetMassLine(0.140, (plotObj[i]->type == 0)?dEdxK_Data1:dEdxK_MC1, (plotObj[i]->type == 0)?dEdxC_Data1:dEdxC_MC1);
            PionLine->SetLineColor(1);
            PionLine->SetLineWidth(2);
            PionLine->SetRange(PionLine->GetX(15), PionLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

            TF1* KaonLine = GetMassLine(0.494, (plotObj[i]->type == 0)?dEdxK_Data1:dEdxK_MC1, (plotObj[i]->type == 0)?dEdxC_Data1:dEdxC_MC1);
            KaonLine->SetLineColor(1);
            KaonLine->SetLineWidth(2);
            KaonLine->SetRange(KaonLine->GetX(15), KaonLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

            TF1* ProtonLine = GetMassLine(0.938, (plotObj[i]->type == 0)?dEdxK_Data1:dEdxK_MC1, (plotObj[i]->type == 0)?dEdxC_Data1:dEdxC_MC1);
            ProtonLine->SetLineColor(1);
            ProtonLine->SetLineWidth(2);
            ProtonLine->SetRange(ProtonLine->GetX(15), ProtonLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

            TF1* DeuteronLine = GetMassLine(1.88, (plotObj[i]->type == 0)?dEdxK_Data1:dEdxK_MC1, (plotObj[i]->type == 0)?dEdxC_Data1:dEdxC_MC1);
            DeuteronLine->SetLineColor(1);
            DeuteronLine->SetLineWidth(2);
            DeuteronLine->SetRange(DeuteronLine->GetX(15), DeuteronLine->GetX(plotObj[i]->C[plotObj[i]->StdObjName[j]]+0.1));

            TF1* DeuteronLineFit = GetMassLine(1.88, (plotObj[i]->type == 0)?dEdxK_Data1:dEdxK_MC1, (plotObj[i]->type == 0)?dEdxC_Data1:dEdxC_MC1);
            DeuteronLineFit->SetLineColor(2);
            DeuteronLineFit->SetLineWidth(2);
            DeuteronLineFit->SetRange(0.8,2.0);

            PionLine->Draw("same");
            KaonLine->Draw("same");
            ProtonLine->Draw("same");
            DeuteronLine->Draw("same");
            DeuteronLineFit->Draw("same");
         }
         DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
         SaveCanvas(c1, SaveDir, plotObj[i]->StdObjName[j] + "_" + plotObj[i]->SavePrefix + "_dedxVsP_KCData", true);
         if (PionLine){ delete PionLine; delete KaonLine; delete ProtonLine; delete DeuteronLine; delete DeuteronLineFit; }
         delete c1;

         c1 = new TCanvas("c1", "c1", 600,600);
         c1->SetLogz(true);
         plotObj[i]->HdedxVsEta[j]->SetStats(kFALSE);
         plotObj[i]->HdedxVsEta[j]->GetXaxis()->SetTitle("Eta");
         plotObj[i]->HdedxVsEta[j]->GetYaxis()->SetTitle(plotObj[i]->StdObjName[j].find("Ias")!=std::string::npos?"I_{as}":"dE/dx (MeV/cm)");
         plotObj[i]->HdedxVsEta[j]->SetAxisRange(-2.1,2.1,"X");
         plotObj[i]->HdedxVsEta[j]->Draw("COLZ");
         DrawPreliminary ("", 13.0, plotObj[i]->LegEntry);
         SaveCanvas(c1, SaveDir, plotObj[i]->StdObjName[j] + "_" + plotObj[i]->SavePrefix + "_Eta2D", true);
         delete c1;
      }

      cout<<"qui arriva"<<endl;
      for (size_t j=0; j < plotObj[i]->HitObjName.size(); j++){
	cout <<(plotObj[i]->HitObjName[j]).c_str()<<endl;
         if (plotObj[i]->HitObjName[j].find("SP")==string::npos) continue;
	 cout<<"e` qui il problema"<<endl;
	 cout<<plotObj[i]<<endl;
	 cout<<plotObj[i]->dEdxTemplate[j]<<endl;
	 if (plotObj[i]->dEdxTemplate[j]){
	   cout<<"abbiamo il template"<<endl;
	   plotObj[i]->dEdxTemplate[j]->SetName("Charge_Vs_Path");
	   plotObj[i]->dEdxTemplate[j]->SaveAs (("dEdxTemplate_" + plotObj[i]->HitObjName[j] + "_" + plotObj[i]->SavePrefix + ".root").c_str());
	   cout<<"qui anche"<<endl;
	   MakeMapPlots (plotObj[i]->dEdxTemplate[j], plotObj[i]->HitObjName[j], SaveDir, "Map" + plotObj[i]->SavePrefix);
	   cout <<"qui pure"<<endl;}
/*         
         TH1D* hit_MIP = (TH1D*) GetObjectFromPath (InputFiles[fileIndex], (ObjName[i]+"_Hit").c_str());

         TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
         c1->SetLogy(true);
         hit_MIP->SetStats(kFALSE);
         hit_MIP->GetXaxis()->SetTitle("cluster dE/dx");
         hit_MIP->GetYaxis()->SetTitle("number of hits");
         hit_MIP->Draw("hist");
         DrawPreliminary(extraStringMode1);
         SaveCanvas (c1, SaveDir, ObjName[i]+"_"+SaveNames[fileIndex]+"_Hit");
         delete c1;
*/
         // all the other graphs -- Charge_Vs_XYNLetc.


	 
	 cout<<"qui pure"<<endl;
         string topologies [] = {"IB1", "IB2", "OB1", "OB2", "W1A", "W2A", "W3A", "W1B", "W2B", "W3B", "W4", "W5", "W6", "W7", "Pixel"};

         for (unsigned int g=0;g<15;g++){
            char Id [255]; sprintf (Id, "%02i", g+1);

	    TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
	    if (plotObj[i]->Charge_Vs_XYH[j][g]){
	      plotObj[i]->Charge_Vs_XYH[j][g]->SetStats(kFALSE);
	      plotObj[i]->Charge_Vs_XYH[j][g]->GetXaxis()->SetTitle("local x coordinate");
	      plotObj[i]->Charge_Vs_XYH[j][g]->GetYaxis()->SetTitle("local y coordinate");
	      plotObj[i]->Charge_Vs_XYH[j][g]->SetAxisRange (-7,7,"X");
	      plotObj[i]->Charge_Vs_XYH[j][g]->SetAxisRange (-15,15,"Y");
	      plotObj[i]->Charge_Vs_XYH[j][g]->Draw("COLZ");
	      DrawPreliminary (topologies[g], 13.0, plotObj[i]->LegEntry);
	      SaveCanvas (c1, SaveDir, plotObj[i]->HitObjName[j]+"_"+plotObj[i]->SavePrefix +"_ChargeVsXYH"+string(Id), true);
	    }
	    delete c1;

	    if (plotObj[i]->Charge_Vs_XYHN[j][g]){
	      c1 = new TCanvas ("c1", "c1", 600, 600);
	      plotObj[i]->Charge_Vs_XYHN[j][g]->SetStats(kFALSE);
	      plotObj[i]->Charge_Vs_XYHN[j][g]->GetXaxis()->SetTitle("normalized x coordinate");
	      plotObj[i]->Charge_Vs_XYHN[j][g]->GetYaxis()->SetTitle("normalized y coordinate");
	      plotObj[i]->Charge_Vs_XYHN[j][g]->SetAxisRange (-1.5,1.5,"X");
	      plotObj[i]->Charge_Vs_XYHN[j][g]->SetAxisRange (-1.5,1.5,"Y");
	      plotObj[i]->Charge_Vs_XYHN[j][g]->Draw("COLZ");
	      DrawPreliminary (topologies[g], 13.0, plotObj[i]->LegEntry);
	      SaveCanvas (c1, SaveDir, plotObj[i]->HitObjName[j]+"_"+plotObj[i]->SavePrefix +"_ChargeVsXYHN"+string(Id), true);
	      delete c1;
	    }

	    if(plotObj[i]->Charge_Vs_XYLN[j][g]){
	      c1 = new TCanvas ("c1", "c1", 600, 600);
	      plotObj[i]->Charge_Vs_XYLN[j][g]->SetStats(kFALSE);
	      plotObj[i]->Charge_Vs_XYLN[j][g]->GetXaxis()->SetTitle("normalized x coordinate");
	      plotObj[i]->Charge_Vs_XYLN[j][g]->GetYaxis()->SetTitle("normalized y coordinate");
	      plotObj[i]->Charge_Vs_XYLN[j][g]->SetAxisRange (-1.5,1.5,"X");
	      plotObj[i]->Charge_Vs_XYLN[j][g]->SetAxisRange (-1.5,1.5,"Y");
	      plotObj[i]->Charge_Vs_XYLN[j][g]->Draw("COLZ");
	      DrawPreliminary (topologies[g], 13.0, plotObj[i]->LegEntry);
	      SaveCanvas (c1, SaveDir, plotObj[i]->HitObjName[j]+"_"+plotObj[i]->SavePrefix +"_ChargeVsXYLN"+string(Id), true);
	      delete c1;
	    }
         }
      }
   }
}

void HitPlots (string SaveDir, vector <dEdxPlotObj*> plotObj){
   for (size_t j=0; j < plotObj[0]->HitObjName.size()-1; j++){
      vector <TH1D*> histos;
      vector <string> legend;
      TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
      c1->SetLogy(true);

      double min = 1e-5,
	     max = 2e-1;
      for (size_t i=0; i < plotObj.size(); i++){
         if (plotObj[i]->type==2) continue;
         histos.push_back(plotObj[i]->hit_MIP[j]);
         histos[histos.size()-1]->Scale(1.0/histos[histos.size()-1]->Integral());
         legend.push_back(plotObj[i]->LegEntry);
//         max = (max<histos[histos.size()-1]->GetMaximum())?histos[histos.size()-1]->GetMaximum():max;
//         min = (min>histos[histos.size()-1]->GetMinimum())?histos[histos.size()-1]->GetMinimum():min;
      }

      DrawSuperposedHistos ((TH1**) &(histos[0]), legend, "HIST", "dE/dx (MeV/cm)", "arbitrary units", 0, 5, min, max*12);
      DrawLegend ((TObject**) &(histos[0]), legend, plotObj[0]->HitObjLegend[j], "LP", 0.80, 0.90, 0.30, 0.05);
      DrawPreliminary ("", 13.0, "");
      SaveCanvas (c1, SaveDir, plotObj[0]->HitObjName[j] + "_HHit");
      legend.clear();
      histos.clear();
      delete c1;
   }
}

void ExtractSlope (string SaveDir, vector <dEdxPlotObj*> plotObj){
   vector<TH1D*> histos;

   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type != 0) continue;
      for (size_t j = 0; j < plotObj[i]->HitObjName.size(); j++)
         histos.push_back (plotObj[i]->hit_MIP[j]);
   }

   TF1** myLin = new TF1* [histos.size()];
   vector <double> deltas;
   for (size_t i = 0; i < histos.size(); i++){
      char tmp [100]; sprintf (tmp, "linear_%0d", (int) i);
      myLin[i] = new TF1 (tmp, "[0]*x + [1]", 0.7, 1.3);

	   myLin[i]->SetParName  (0,"k");
	   myLin[i]->SetParName  (1,"n");
	   myLin[i]->SetParLimits(0,  1000, 6000);
	   myLin[i]->SetParLimits(1, -2000, 6000);
      myLin[i]->SetParameter (0, 3000);
      myLin[i]->SetParameter (1, 3000);
	   myLin[i]->SetLineWidth(2);
	   myLin[i]->SetLineColor(2);
      histos[i]->Fit (myLin[i], "Q M R E 0");

      deltas.push_back (myLin[i]->Eval(1.5) - myLin[i]->Eval(1.0));
   }

   //TCanvas* c1 = new TCanvas ("c1", "c1", 600, 600);
   //for Run257805 delta = 9.270155e+02 -- interactively
   //or            delta = 1.112419e+03 -- boundaries between 0.6 and 1.5
}

void SaveKC (vector<dEdxPlotObj*> plotObj, string filename){
   FILE* pFile = fopen (filename.c_str(), "w");
   for (size_t i = 0; i < plotObj.size(); i++){
      if (plotObj[i]->type == 2) continue;
      fprintf (pFile, "%s\n--------------------------\n", plotObj[i]->LegEntry.c_str());
      for (size_t j = 0; j < plotObj[i]->StdObjName.size(); j++){
         if (plotObj[i]->StdObjName[j].find("Ias")!=string::npos) continue;
         fprintf (pFile, "%s\tK = %.3lf +/- %.3lf\tC = %.3lf +/- %.3lf\n",
               plotObj[i]->StdObjLegend[j].c_str(),
               plotObj[i]->K[plotObj[i]->StdObjName[j]],
               plotObj[i]->Kerr[plotObj[i]->StdObjName[j]],
               plotObj[i]->C[plotObj[i]->StdObjName[j]],
               plotObj[i]->Cerr[plotObj[i]->StdObjName[j]]);
      }
      fprintf (pFile, "\n");
   }
   fclose (pFile);
}


