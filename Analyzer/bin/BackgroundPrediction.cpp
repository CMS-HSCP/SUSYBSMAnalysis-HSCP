// Original Author:  Loic Quertenmont

// ~~~~~~~~~c++ include files ~~~~~~~~~
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <exception>
#include <unordered_map>

// ~~~~~~~~~ ROOT include files ~~~~~~~~~
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TObject.h"
#include "TVector3.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TProfile.h"
//#include "TCanvas.h"

#include <TBenchmark.h>
#include <TSystemDirectory.h>
#include <TSystem.h>
//#include "tdrstyle.C"

#include "../interface/CommonFunction.h"
#include "../interface/SaturationCorrection.h"
#include "tdrstyle.h"
#include "HistoTools.h"
#include "ArgumentParser.h"

class ArgumentParser;

using namespace std;

void Analysis_Step2_BackgroundPrediction(TFile* InputFile, int TypeMode = 0)
{
   bool symmetrizeHistos = false;
   ////if(InputPattern=="COMPILE")return;

   /*//WAIT//setTDRStyle();
   gStyle->SetPadTopMargin   (0.06);
   gStyle->SetPadBottomMargin(0.12);
   gStyle->SetPadRightMargin (0.16);
   gStyle->SetPadLeftMargin  (0.14);
   gStyle->SetTitleSize(0.04, "XYZ");
   gStyle->SetTitleXOffset(1.1);
   gStyle->SetTitleYOffset(1.45);
   gStyle->SetPalette(1);
   gStyle->SetNdivisions(505);*/

   unsigned int NPseudoExp = 100; //Number of PseudoExperiment to run

   //7TeV DXY/DZ/ANGLE 85/326  86/327   10/251
   double CosmicVetoInEfficiency7TeV    = 0.26 * 0.26 * 0.04 ; 
   double CosmicVetoInEfficiency7TeVErr = sqrt(pow(0.03*0.26*0.04,2) + pow(0.26*0.03*0.04,2) +pow(0.26*0.04*0.01,2) + pow(0.50*CosmicVetoInEfficiency7TeV,2));  //add 50% syst uncertainty

   //8TeV DXY/DZ/ANGLE 19/22    0/3      0/3
   double CosmicVetoInEfficiency8TeV    = 0.86 * 0.26 * 0.04 ; 
   double CosmicVetoInEfficiency8TeVErr = sqrt(pow(0.20*0.26*0.04,2) + pow(0.86*0.03*0.04,2) + pow(0.86*0.04*0.01,2) + pow(0.50*CosmicVetoInEfficiency8TeV,2)); //add 50% syst uncertainty

   //13TeV DXY/DZ/ANGLE 19/22    0/3      0/3 #FIXME TO BE UPDATED... COPY FROM 8TeV
   double CosmicVetoInEfficiency13TeV    = 0.86 * 0.26 * 0.04 ; 
   double CosmicVetoInEfficiency13TeVErr = sqrt(pow(0.20*0.26*0.04,2) + pow(0.86*0.03*0.04,2) + pow(0.86*0.04*0.01,2) + pow(0.50*CosmicVetoInEfficiency13TeV,2)); //add 50% syst uncertainty


   double CosmicVetoInEfficiency    = CosmicVetoInEfficiency7TeV;
   double CosmicVetoInEfficiencyErr = CosmicVetoInEfficiency7TeVErr;
   

    //string Input     = InputPattern + "Histos.root";
    //TFile* InputFile = new TFile(Input.c_str(), "UPDATE");
    //WAIT//TypeMode = TypeFromPattern(InputPattern);
    
    int PredBins=0;
    const int MaxPredBins =   6;
    int MassNBins           = 400;
    double             MassHistoUpperBound = 4000;

    double dEdxK_MC(1.0), dEdxC_MC(1.0);

    TDirectory* dir = (TDirectory*)InputFile->FindObjectAny("analyzer");
    if (!dir){
        dir = InputFile;
    }
    TList* list = dir->GetListOfKeys();
    //Do two loops, one for the actual background prediction and one for the 
    //region with TOF<1
    for(unsigned int S=0; S<2; S++) {
	    string Suffix="";
        if(S==1) Suffix="_Flip";

        TH1D*  HCuts_Pt       = (TH1D*)GetObjectFromPath(dir, ("HCuts_Pt" + Suffix).c_str());
	    TH1D*  HCuts_I        = (TH1D*)GetObjectFromPath(dir, ("HCuts_I" + Suffix).c_str());
	    TH1D*  HCuts_TOF      = (TH1D*)GetObjectFromPath(dir, ("HCuts_TOF" + Suffix).c_str());
  
        for(int d=0;d<list->GetEntries();d++){
            TObject *key = list->At(d);
            if(!key->IsFolder()) continue;
	        string DirName;
	        DirName = key->GetName();
            //if (DirName.find("Data13TeV16")==string::npos) continue;
            //if (DirName.find("Data13TeV16G")!=string::npos) continue;
	        if(DirName.find("Cosmic")!=string::npos) continue;

            if(DirName.find("7TeV")!=string::npos){
                CosmicVetoInEfficiency    = CosmicVetoInEfficiency7TeV;
                CosmicVetoInEfficiencyErr = CosmicVetoInEfficiency7TeVErr;
            }else if(DirName.find("8TeV")!=string::npos){
                CosmicVetoInEfficiency    = CosmicVetoInEfficiency8TeV;
                CosmicVetoInEfficiencyErr = CosmicVetoInEfficiency8TeVErr;
            }else{
                CosmicVetoInEfficiency    = CosmicVetoInEfficiency13TeV;
                CosmicVetoInEfficiencyErr = CosmicVetoInEfficiency13TeVErr;
            }

            bool is2016 = (DirName.find("13TeV16")!=string::npos);

            /*// to chyba i tak nigddzie nie jest uzywane
            dEdxK_MC   = is2016?dEdxK_MC16:dEdxK_MC15;
            dEdxC_MC   = is2016?dEdxC_MC16:dEdxC_MC15;


            // parametry K C and postG

             dEdxK_Data = 2.062;
             dEdxC_Data = 3.430;
 
            // Usrenione wazone parametry dla KC and preG
            if(DirName.find("13TeV16G")!=string::npos){
                dEdxK_Data = 2.278;
                dEdxC_Data = 3.694;
            }*/



            TDirectory* directory = dir->GetDirectory(key->GetName()); //TDirectory* directory = InputFile->GetDirectory(file->GetName());
            
	        directory->cd();

            TH1D*  H_A            = (TH1D*)GetObjectFromPath(directory, ("H_A" + Suffix).c_str());      if(!H_A)continue; //ABCD INFO NOT SAVED IN THIS DIRECTORY --> Skip it
            TH1D*  H_B            = (TH1D*)GetObjectFromPath(directory, ("H_B" + Suffix).c_str());
            TH1D*  H_C            = (TH1D*)GetObjectFromPath(directory, ("H_C" + Suffix).c_str());
            TH1D*  H_D            = (TH1D*)GetObjectFromPath(directory, ("H_D" + Suffix).c_str());
            TH1D*  H_E            = (TH1D*)GetObjectFromPath(directory, ("H_E" + Suffix).c_str());
            TH1D*  H_F            = (TH1D*)GetObjectFromPath(directory, ("H_F" + Suffix).c_str());
            TH1D*  H_G            = (TH1D*)GetObjectFromPath(directory, ("H_G" + Suffix).c_str());
            TH1D*  H_H            = (TH1D*)GetObjectFromPath(directory, ("H_H" + Suffix).c_str());
        
            TH1D*  H_A_Cosmic     = (S==0) ? (TH1D*)GetObjectFromPath(directory, ("H_A_Flip" + Suffix).c_str()) : NULL;
            TH1D*  H_B_Cosmic     = (S==0) ? (TH1D*)GetObjectFromPath(directory, ("H_B_Flip" + Suffix).c_str()) : NULL;
            TH1D*  H_C_Cosmic     = (S==0) ? (TH1D*)GetObjectFromPath(directory, ("H_C_Flip" + Suffix).c_str()) : NULL;
            TH1D*  H_D_Cosmic     = (S==0) ? (TH1D*)GetObjectFromPath(directory, ("H_D_Flip" + Suffix).c_str()) : NULL;

            TH1D*  H_B_Binned[MaxPredBins];
            TH1D*  H_F_Binned[MaxPredBins];
            TH1D*  H_H_Binned[MaxPredBins];
            TH1D*  H_P_Binned[MaxPredBins];
            if(TypeMode==3) PredBins=6;
	        for(int i=0; i<PredBins; i++) {
	            string Version=Suffix;
	            char Bin[1024];
	            sprintf(Bin,"_%i", i);
	            Version.append(Bin);
	            H_B_Binned[i]            = (TH1D*)GetObjectFromPath(directory, ("H_B_Binned" + Version).c_str());
	            H_F_Binned[i]            = (TH1D*)GetObjectFromPath(directory, ("H_F_Binned" + Version).c_str());
	            H_H_Binned[i]            = (TH1D*)GetObjectFromPath(directory, ("H_H_Binned" + Version).c_str());
	        }

            TH3D*  Pred_EtaP      = (TH3D*)GetObjectFromPath(directory, ("Pred_EtaP" + Suffix).c_str());
	        TH2D*  Pred_I         = (TH2D*)GetObjectFromPath(directory, ("Pred_I" + Suffix).c_str());
	        TH2D*  Pred_TOF       = (TH2D*)GetObjectFromPath(directory, ("Pred_TOF" + Suffix).c_str());
	        TH2D*  Pred_EtaB      = (TH2D*)GetObjectFromPath(directory, ("Pred_EtaB" + Suffix).c_str());
	        TH2D*  Pred_EtaS      = (TH2D*)GetObjectFromPath(directory, ("Pred_EtaS" + Suffix).c_str());
	        TH2D*  Pred_EtaS2     = (TH2D*)GetObjectFromPath(directory, ("Pred_EtaS2" + Suffix).c_str());
            //----------------------------------------------------------------------------------------------------Prediction in eta bins  PZ1
            //                                                                                                                      pT I ibeta
            TH3D*  PDF_G_EtaP      = (TH3D*)GetObjectFromPath(directory, ("PDF_G_EtaP"   + Suffix).c_str());       //G <-> x o o
            TH3D*  PDF_C_EtaP      = (TH3D*)GetObjectFromPath(directory, ("PDF_C_EtaP"   + Suffix).c_str());       //C <-> x o x
            TH3D*  PDF_F_EtaICK    = (TH3D*)GetObjectFromPath(directory, ("PDF_F_EtaICK" + Suffix).c_str());       //F <-> o x o
            TH3D*  PDF_B_EtaICK    = (TH3D*)GetObjectFromPath(directory, ("PDF_B_EtaICK" + Suffix).c_str());       //B <-> o x x
            TH2D*  PDF_E_Eta       = (TH2D*)GetObjectFromPath(directory, ("PDF_E_Eta"    + Suffix).c_str());       //E <-> o o o
            TH2D*  PDF_A_Eta       = (TH2D*)GetObjectFromPath(directory, ("PDF_A_Eta"    + Suffix).c_str());       //A <-> o o x
            TH3D*  PDF_H_EtaMass   = (TH3D*)GetObjectFromPath(directory, ("PDF_H_EtaMass"+ Suffix).c_str());       //H <-> x x o
            //----------------------------------------------------------------------------------------------------Prediction in eta bins


            TH2D*  H_D_DzSidebands= (TH2D*)GetObjectFromPath(directory, "H_D_DzSidebands");
	        TH1D*  H_D_CosmicMO=NULL;
            TH2D*  H_D_DzSidebands_Cosmic=NULL;

	        TH1D*  H_B_Cosmic_Binned[MaxPredBins];
	        TH1D*  H_F_Cosmic_Binned[MaxPredBins];
            TH1D*  H_H_Cosmic_Binned[MaxPredBins];

	        if(TypeMode==3 && DirName.find("Data")!=string::npos) {
	            //Only 2012 sample has pure cosmic sample, as only ratio used can use 2012 sample to make 2011 cosmic prediction	    
                string CosmicDir = "Cosmic8TeV";
                if(DirName.find("13TeV")) CosmicDir = "Cosmic13TeV";
	            //string CosmicDir = DirName.replace(0, 4, "Cosmic");
	            H_D_DzSidebands_Cosmic = (TH2D*)GetObjectFromPath(InputFile, (CosmicDir + "/H_D_DzSidebands").c_str());
	            H_D_CosmicMO           = (TH1D*)GetObjectFromPath(InputFile, (CosmicDir + "/H_D" + Suffix).c_str());

	            for(int i=0; i<PredBins; i++) {
	                string Version=Suffix;
	                char Bin[1024];
	                sprintf(Bin,"_%i", i);
	                Version.append(Bin);

	                H_B_Cosmic_Binned[i]            = (TH1D*)GetObjectFromPath(InputFile, (CosmicDir + "/H_B_Binned" + Version).c_str());
	                H_F_Cosmic_Binned[i]            = (TH1D*)GetObjectFromPath(InputFile, (CosmicDir + "/H_F_Binned" + Version).c_str());
	                H_H_Cosmic_Binned[i]            = (TH1D*)GetObjectFromPath(InputFile, (CosmicDir + "/H_H_Binned" + Version).c_str());
	            }
	        }

            //erase histogram created at previous iteration
	        //directory->Delete(("Pred_P" + Suffix + ";*").c_str());
            //directory->Delete(("Pred_Mass" + Suffix + ";*").c_str());
            //directory->Delete(("Pred_MassTOF" + Suffix + ";*").c_str());
            //directory->Delete(("Pred_MassComb" + Suffix + ";*").c_str());
            //directory->Delete(("H_P" + Suffix + ";*").c_str());
            //directory->Delete(("H_P_Coll" + Suffix + ";*").c_str());
            //directory->Delete(("H_P_Cosmic" + Suffix + ";*").c_str());

            //take data histogram to save the resulting momentum distribution
            TH1D*  H_P            = (TH1D*)GetObjectFromPath(directory, ("H_D" + Suffix).c_str())->Clone(("H_P" + Suffix).c_str());                   H_P->Reset();
	        TH1D*  H_P_Coll       = (TH1D*)H_P->Clone(("H_P_Coll" + Suffix).c_str());
            TH1D*  H_P_Cosmic     = (TH1D*)H_P->Clone(("H_P_Cosmic" + Suffix).c_str());
            TH2D*  Pred_Mass      = (TH2D*)GetObjectFromPath(directory, ("Mass" + Suffix).c_str())->Clone(("Pred_Mass" + Suffix).c_str());         Pred_Mass->Reset();
	        TH2D*  Pred_MassTOF   = (TH2D*)GetObjectFromPath(directory, ("MassTOF" + Suffix).c_str())->Clone(("Pred_MassTOF" + Suffix).c_str());  Pred_MassTOF->Reset();
	        TH2D*  Pred_MassComb  = (TH2D*)GetObjectFromPath(directory, ("MassComb" + Suffix).c_str())->Clone(("Pred_MassComb" + Suffix).c_str()); Pred_MassComb->Reset();
	        TH2D*  Pred_P         = (TH2D*)GetObjectFromPath(directory, ("RegionD_P" + Suffix).c_str())->Clone(("Pred_P" + Suffix).c_str());           Pred_P->Reset();
            //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ2
            TH2D*  Pred_Mass_CB      = (TH2D*)GetObjectFromPath(directory, ("Mass" + Suffix).c_str())->Clone(("Pred_Mass_CB" + Suffix).c_str());         
            Pred_Mass_CB ->Reset();
            TH2D*  Pred_Mass_GFA     = (TH2D*)GetObjectFromPath(directory, ("Mass" + Suffix).c_str())->Clone(("Pred_Mass_GFA"+ Suffix).c_str());         
            Pred_Mass_GFA->Reset();
            TH2D*  Pred_Mass_S       = (TH2D*)GetObjectFromPath(directory, ("Mass" + Suffix).c_str())->Clone(("Pred_Mass_S"+   Suffix).c_str());         
            Pred_Mass_S->  Reset();
            TH2D*  Pred_Mass_GB      = (TH2D*)GetObjectFromPath(directory, ("Mass" + Suffix).c_str())->Clone(("Pred_Mass_GB" + Suffix).c_str());         
            Pred_Mass_GB ->Reset();
            TH2D*  Pred_Mass_HA      = (TH2D*)GetObjectFromPath(directory, ("Mass" + Suffix).c_str())->Clone(("Pred_Mass_HA" + Suffix).c_str());         
            Pred_Mass_HA ->Reset();
            TH2D*  Pred_Mass_CF      = (TH2D*)GetObjectFromPath(directory, ("Mass" + Suffix).c_str())->Clone(("Pred_Mass_CF" + Suffix).c_str());         
            Pred_Mass_CF ->Reset();
            //----------------------------------------------------------------------------------------------------Prediction in eta bins

            for(int i=0; i<PredBins; i++) {
                string Version=Suffix;
                char Bin[1024];
                sprintf(Bin,"_%i", i);
                Version.append(Bin);
                H_P_Binned[i]       = (TH1D*)H_P->Clone(("H_P_Binned" + Version).c_str());
            }

            printf("Making prediction for %s\n",directory->GetName());
            //////////////////////////////////////////////////      MAKING THE PREDICTION
            for(unsigned int CutIndex=0;CutIndex<(unsigned int)HCuts_Pt->GetXaxis()->GetNbins();CutIndex++){
                //if(CutIndex<86 || CutIndex>87)continue;

                double A=H_A->GetBinContent(CutIndex+1);  double AErr = sqrt(A);
                double B=H_B->GetBinContent(CutIndex+1);  double BErr = sqrt(B);
                double C=H_C->GetBinContent(CutIndex+1);  double CErr = sqrt(C);
                double D=H_D->GetBinContent(CutIndex+1); // double DErr = sqrt(D); //UNUSED
                double E=H_E->GetBinContent(CutIndex+1); // double EErr = sqrt(E);
                double F=H_F->GetBinContent(CutIndex+1); // double FErr = sqrt(F);
                double G=H_G->GetBinContent(CutIndex+1); // double GErr = sqrt(G);
                double H=H_H->GetBinContent(CutIndex+1); // double HErr = sqrt(H);

                double A_Cosmic=1, B_Cosmic=0, C_Cosmic=0, D_Cosmic=0;
                double AErr_Cosmic=0, BErr_Cosmic=0, CErr_Cosmic=0, DErr_Cosmic=0;
                if(S==0 && TypeMode==5){
                    A_Cosmic=H_A_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiency;
                    B_Cosmic=H_B_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiency;
                    C_Cosmic=H_C_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiency;
                    D_Cosmic=H_D_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiency;

                    AErr_Cosmic=sqrt( pow(sqrt(H_A_Cosmic->GetBinContent(CutIndex+1)) * CosmicVetoInEfficiency,2) + pow(H_A_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiencyErr,2) );
                    BErr_Cosmic=sqrt( pow(sqrt(H_B_Cosmic->GetBinContent(CutIndex+1)) * CosmicVetoInEfficiency,2) + pow(H_B_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiencyErr,2) );
                    CErr_Cosmic=sqrt( pow(sqrt(H_C_Cosmic->GetBinContent(CutIndex+1)) * CosmicVetoInEfficiency,2) + pow(H_C_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiencyErr,2) );
                    DErr_Cosmic=sqrt( pow(sqrt(H_D_Cosmic->GetBinContent(CutIndex+1)) * CosmicVetoInEfficiency,2) + pow(H_D_Cosmic->GetBinContent(CutIndex+1) * CosmicVetoInEfficiencyErr,2) );

                    if(CutIndex==44){
                        printf("scale factor = %f+-%f\n", CosmicVetoInEfficiency, CosmicVetoInEfficiencyErr);
                        printf("%E+-%E   %E+-%E   %E+-%E\n", A_Cosmic, AErr_Cosmic, B_Cosmic, BErr_Cosmic, C_Cosmic, CErr_Cosmic);
                    }

                    A = A - A_Cosmic;    AErr = sqrt(AErr*AErr + AErr_Cosmic*AErr_Cosmic);
                    B = B - B_Cosmic;    BErr = sqrt(BErr*BErr + BErr_Cosmic*BErr_Cosmic);
                    C = C - C_Cosmic;    CErr = sqrt(CErr*CErr + CErr_Cosmic*CErr_Cosmic);
                    //D = D - D_Cosmic;    DErr = sqrt(DErr*DErr + DErr_Cosmic*DErr_Cosmic);
                }

                double B_Binned[MaxPredBins];
                double F_Binned[MaxPredBins];
                double H_Binned[MaxPredBins];
                //double P_Binned[MaxPredBins];
                //double Perr_Binned[MaxPredBins];

                double B_Cosmic_Binned[MaxPredBins];
                double F_Cosmic_Binned[MaxPredBins];
                double H_Cosmic_Binned[MaxPredBins];

                for(int i=0; i<PredBins; i++) {
                    B_Binned[i]=H_B_Binned[i]->GetBinContent(CutIndex+1);
                    F_Binned[i]=H_F_Binned[i]->GetBinContent(CutIndex+1);
                    H_Binned[i]=H_H_Binned[i]->GetBinContent(CutIndex+1);
                    //P_Binned[i]=0;
                    //Perr_Binned[i]=0;

                    B_Cosmic_Binned[i]=H_B_Cosmic_Binned[i]->GetBinContent(CutIndex+1);
                    F_Cosmic_Binned[i]=H_F_Cosmic_Binned[i]->GetBinContent(CutIndex+1);
                    H_Cosmic_Binned[i]=H_H_Cosmic_Binned[i]->GetBinContent(CutIndex+1);
                }

                double P=0;
                double Perr=0;
	            double P_Coll=0;
	            double Perr_Coll=0;
	            double P_Cosmic=0;
	            double Perr_Cosmic=0;

                printf("%4i --> Pt>%7.2f  I>%6.2f  TOF>%+5.2f --> A=%6.2E B=%6.2E C=%6.2E D=%6.2E E=%6.2E F=%6.2E G=%6.2E H=%6.2E",CutIndex,HCuts_Pt->GetBinContent(CutIndex+1), HCuts_I->GetBinContent(CutIndex+1), HCuts_TOF->GetBinContent(CutIndex+1),A, B, C, D, E, F, G, H );

                if(E>0){
	                //Prediction in Pt-Is-TOF plane
                    P    = (A*F*G)/(E*E);
                    Perr = sqrt( ((pow(F*G,2)* A + pow(A*G,2)*F + pow(A*F,2)*G)/pow(E,4)) + (pow((2*A*F*G)/pow(E,3),2)*E));
                }else if(A>0){
	                //Prediction in Pt-Is plane
                    P    = ((C*B)/A);
                    Perr = sqrt( pow(CErr*B/A,2) + pow(BErr*C/A,2) + pow((AErr*B*C/(A*A)),2) );

                    if(S==0 && TypeMode==5){
                        P_Coll      = P;
                        Perr_Coll   = Perr;
                        P_Cosmic    = -1;
                        if(A_Cosmic>0 && B_Cosmic>0 && C_Cosmic>0){
                            P_Cosmic    = ((C_Cosmic*B_Cosmic)/A_Cosmic);
                            Perr_Cosmic = sqrt( pow(CErr_Cosmic*B_Cosmic/A_Cosmic,2) + pow(BErr_Cosmic*C_Cosmic/A_Cosmic,2) + pow(AErr_Cosmic*B_Cosmic*C_Cosmic/(A_Cosmic*A_Cosmic),2) );
                        }else if(D_Cosmic>0){
                            P_Cosmic = D_Cosmic * CosmicVetoInEfficiency;
                            Perr_Cosmic = sqrt( pow(DErr_Cosmic * CosmicVetoInEfficiency,2) + pow(D_Cosmic*CosmicVetoInEfficiencyErr,2));
                        }else{
                            P_Cosmic = 3/2.0 * CosmicVetoInEfficiency;
                            Perr_Cosmic = P_Cosmic; //100% uncertainty
                        }
                        P           = P_Coll + P_Cosmic;
                        Perr        = sqrt( Perr_Coll*Perr_Coll + Perr_Cosmic*Perr_Cosmic);
                    }
	            }else if(F>0){
                    //Predict the number of cosmics passing all cuts as number passing in dz sideband times the ratio of tracks in the sideband
                    //vs number in central region as determined by pure cosmic sample
                    //Multile sidebands are made to check for background consistency, the fifth one is used for the actual prediction
                    double D_Sideband = 0;
                    double D_Sideband_Cosmic = 0;
                    if(DirName.find("Data")!=string::npos) {
                        D_Sideband = H_D_DzSidebands->GetBinContent(CutIndex+1, 5);
                        double D_Cosmic = H_D_CosmicMO->GetBinContent(CutIndex+1);
                        D_Sideband_Cosmic = H_D_DzSidebands_Cosmic->GetBinContent(CutIndex+1, 5);
                        if(D_Sideband_Cosmic>0) {
                            P_Cosmic = D_Sideband * D_Cosmic / D_Sideband_Cosmic;
                            Perr_Cosmic = sqrt( (pow(D_Cosmic/D_Sideband_Cosmic,2)*D_Sideband) + (pow(D_Sideband/D_Sideband_Cosmic,2)*D_Cosmic) + (pow((D_Cosmic*(D_Sideband)/(D_Sideband_Cosmic*D_Sideband_Cosmic)),2)*D_Sideband_Cosmic) );
                        }
	                }

                    //Prediction in Pt-TOF plane
                    for(int i=0; i<PredBins; i++) {
                        //Subtract the expected cosmic tracks from each region
                        double B_Bin = B_Binned[i] - B_Cosmic_Binned[i]*D_Sideband/D_Sideband_Cosmic;
                        double F_Bin = F_Binned[i] - F_Cosmic_Binned[i]*D_Sideband/D_Sideband_Cosmic;
                        double H_Bin = H_Binned[i] - H_Cosmic_Binned[i]*D_Sideband/D_Sideband_Cosmic;

                        double Berr = sqrt(B_Binned[i] + (pow(B_Cosmic_Binned[i]/D_Sideband_Cosmic,2)*D_Sideband) + (pow(D_Sideband/D_Sideband_Cosmic,2)*B_Cosmic_Binned[i]) + (pow((B_Cosmic_Binned[i]*(D_Sideband)/(D_Sideband_Cosmic*D_Sideband_Cosmic)),2)*D_Sideband_Cosmic) );
                        double Ferr = sqrt(F_Binned[i] + (pow(F_Cosmic_Binned[i]/D_Sideband_Cosmic,2)*D_Sideband) + (pow(D_Sideband/D_Sideband_Cosmic,2)*F_Cosmic_Binned[i]) + (pow((F_Cosmic_Binned[i]*(D_Sideband)/(D_Sideband_Cosmic*D_Sideband_Cosmic)),2)*D_Sideband_Cosmic) );
                        double Herr = sqrt(H_Binned[i] + (pow(H_Cosmic_Binned[i]/D_Sideband_Cosmic,2)*D_Sideband) + (pow(D_Sideband/D_Sideband_Cosmic,2)*H_Cosmic_Binned[i]) + (pow((H_Cosmic_Binned[i]*(D_Sideband)/(D_Sideband_Cosmic*D_Sideband_Cosmic)),2)*D_Sideband_Cosmic) );


                        double P_Binned = ((H_Bin*B_Bin)/F_Bin);
                        double Perr_Binned = (pow(Berr/Ferr,2)*Herr) + (pow(Herr/Ferr,2)*Berr) + (pow((Berr*(Herr)/(Ferr*Ferr)),2)*Ferr);

                        H_P_Binned[i]->SetBinContent(CutIndex+1, P_Binned);
                        H_P_Binned[i]->SetBinError(CutIndex+1, sqrt(Perr_Binned));
                        P_Coll    += P_Binned;
                        Perr_Coll += Perr_Binned;
                    }
                    Perr_Coll = sqrt(Perr_Coll);

                    P    = P_Coll + P_Cosmic;
                    Perr = sqrt(Perr_Coll*Perr_Coll + Perr_Cosmic*Perr_Cosmic);
                    //Add in systematic contribution, now done in step 6
                    //Perr = sqrt(Perr*Perr + P_Coll*P_Coll*0.2*0.2 + P_Cosmic*P_Cosmic*0.8*0.8);
                }else if(G>0){
                    //Prediction in Ias-TOF plane
                    P    = ((C*H)/G);
                    Perr = sqrt( (pow(H/G,2)*C) + (pow(C/G,2)*H) + (pow((H*(C)/(G*G)),2)*G) );
                }

                H_P_Coll->SetBinContent(CutIndex+1, P_Coll);
                H_P_Coll->SetBinError  (CutIndex+1, Perr_Coll);

                H_P_Cosmic->SetBinContent(CutIndex+1,P_Cosmic);
                H_P_Cosmic->SetBinError  (CutIndex+1,Perr_Cosmic);

                H_P->SetBinContent(CutIndex+1,P);
                H_P->SetBinError  (CutIndex+1,Perr);

                if(P==0 || isnan((float)P)) {printf("\n"); continue;} //Skip this CutIndex --> No Prediction possible

                printf(" --> D=%6.2E vs Pred = %6.2E +- %6.2E (%6.2E%%)\n", D, P,  Perr, 100.0*Perr/P );
                if(TypeMode>2)continue; //Need to compute mass predicted distribution ONLY for TkOnly and TkTOF

                TH1D* Pred_EtaB_Proj     = Pred_EtaB ->ProjectionY("ProjEtaB" ,CutIndex+1,CutIndex+1);  
                TH1D* Pred_EtaS_Proj     = Pred_EtaS ->ProjectionY("ProjEtaS" ,CutIndex+1,CutIndex+1); 
                TH1D* Pred_EtaS2_Proj    = Pred_EtaS2->ProjectionY("ProjEtaS2",CutIndex+1,CutIndex+1);
                // here

                //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ3
                TH1D* PDF_A_Eta_2D    = PDF_A_Eta ->ProjectionY("ProjAEta" ,CutIndex+1,CutIndex+1);  
                TH1D* PDF_E_Eta_2D    = PDF_E_Eta ->ProjectionY("ProjEEta" ,CutIndex+1,CutIndex+1);  
                //----------------------------------------------------------------------------------------------------Prediction in eta bins

                if (CutIndex==4) compareForwardToBackwardWeights (Pred_EtaS_Proj, Pred_EtaB_Proj);
                if (symmetrizeHistos) symmetrizeHisto(Pred_EtaB_Proj,TypeMode);//WAIT//
                if (symmetrizeHistos) symmetrizeHisto(Pred_EtaS_Proj,TypeMode);
                if (symmetrizeHistos) symmetrizeHisto(Pred_EtaS2_Proj,TypeMode);

                TH1D* Pred_EtaB_Proj_PE  = (TH1D*)Pred_EtaB_Proj ->Clone("Pred_EtaB_Proj_PE");  Pred_EtaB_Proj_PE ->Reset();
                TH1D* Pred_EtaS_Proj_PE  = (TH1D*)Pred_EtaS_Proj ->Clone("Pred_EtaS_Proj_PE");  Pred_EtaS_Proj_PE ->Reset();
                TH1D* Pred_EtaS2_Proj_PE = (TH1D*)Pred_EtaS2_Proj->Clone("Pred_EtaS2_Proj_PE"); Pred_EtaS2_Proj_PE->Reset();

                Pred_EtaP->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
                TH2D* Pred_EtaPWeighted    = (TH2D*)Pred_EtaP->Project3D("zy");
                //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ4
                PDF_G_EtaP->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
                TH2D* PDF_G_EtaP_2D    =  (TH2D*)PDF_G_EtaP->Project3D("zy"); PDF_C_EtaP->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
                TH2D* PDF_C_EtaP_2D    =  (TH2D*)PDF_C_EtaP->Project3D("zy");

                PDF_F_EtaICK->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
                TH2D* PDF_F_EtaICK_2D    =  (TH2D*)PDF_F_EtaICK->Project3D("zy"); PDF_B_EtaICK->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
                TH2D* PDF_B_EtaICK_2D    =  (TH2D*)PDF_B_EtaICK->Project3D("zy");

                PDF_H_EtaMass->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
                TH2D* PDF_H_EtaMass_2D    =  (TH2D*)PDF_H_EtaMass->Project3D("zy");
                //----------------------------------------------------------------------------------------------------Prediction in eta bins
                if (symmetrizeHistos) symmetrizeHisto(Pred_EtaPWeighted,TypeMode);
                TH2D* Pred_EtaPWeighted_PE = (TH2D*)Pred_EtaPWeighted->Clone("Pred_EtaPWeightedPE");   Pred_EtaPWeighted_PE->Reset();

                TH1D* Pred_I_Proj = Pred_I->ProjectionY("ProjI",CutIndex+1,CutIndex+1);
                TH1D* Pred_T_Proj = Pred_TOF->ProjectionY("ProjT",CutIndex+1,CutIndex+1);
                TH1D* Pred_I_ProjPE = (TH1D*) Pred_I_Proj->Clone("Pred_I_ProjPE"); Pred_I_ProjPE->Reset();
                TH1D* Pred_T_ProjPE = (TH1D*) Pred_T_Proj->Clone("Pred_T_ProjPE"); Pred_T_ProjPE->Reset();

                TH2D* Pred_Prof_Mass     =  new TH2D("Pred_Prof_Mass"    ,"Pred_Prof_Mass"    ,MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp); 
                TH2D* Pred_Prof_MassTOF  =  new TH2D("Pred_Prof_MassTOF" ,"Pred_Prof_MassTOF" ,MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp);  
                TH2D* Pred_Prof_MassComb =  new TH2D("Pred_Prof_MassComb","Pred_Prof_MassComb",MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp);
                //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ5
                TH2D* Pred_Prof_Mass_CB     =  new TH2D("Pred_Prof_Mass_CB"    ,"Pred_Prof_Mass_CB"    ,MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp); 
                TH2D* Pred_Prof_Mass_GFA    =  new TH2D("Pred_Prof_Mass_GFA"   ,"Pred_Prof_Mass_GFA"   ,MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp); 
                TH2D* Pred_Prof_Mass_GB     =  new TH2D("Pred_Prof_Mass_GB"    ,"Pred_Prof_Mass_GB"    ,MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp); 
                TH2D* Pred_Prof_Mass_HA     =  new TH2D("Pred_Prof_Mass_HA"    ,"Pred_Prof_Mass_HA"    ,MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp); 
                TH2D* Pred_Prof_Mass_CF     =  new TH2D("Pred_Prof_Mass_CF"    ,"Pred_Prof_Mass_CF"    ,MassNBins,0,MassHistoUpperBound, NPseudoExp, 0, NPseudoExp); 
                //----------------------------------------------------------------------------------------------------Prediction in eta bins

                for(int x=0;x<=Pred_Mass->GetNbinsY()+1;x++){ // taking the overflow as well
                    for(unsigned int pe=0;pe<NPseudoExp;pe++){
                        Pred_Prof_Mass    ->SetBinContent(x, pe, 0);
                        Pred_Prof_MassTOF ->SetBinContent(x, pe, 0);
                        Pred_Prof_MassComb->SetBinContent(x, pe, 0);
                        //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ6
                        Pred_Prof_Mass_CB ->SetBinContent(x, pe, 0);
                        Pred_Prof_Mass_GFA->SetBinContent(x, pe, 0);
                        Pred_Prof_Mass_GB ->SetBinContent(x, pe, 0);
                        Pred_Prof_Mass_HA ->SetBinContent(x, pe, 0);
                        Pred_Prof_Mass_CF ->SetBinContent(x, pe, 0);
                        //----------------------------------------------------------------------------------------------------Prediction in eta bins
                    }
                }

                TRandom3* RNG = new TRandom3();
                printf("Predicting (%4i / %4i)     :",CutIndex+1,HCuts_Pt->GetXaxis()->GetNbins());
                int TreeStep = NPseudoExp/50;if(TreeStep==0)TreeStep=1;
                // loop over pseudo-experiments
                for(unsigned int pe=0;pe<NPseudoExp;pe++){    
                    if(pe%TreeStep==0){printf(".");fflush(stdout);}

                    TH1D* tmpH_Mass     =  new TH1D("tmpH_Mass"    ,"tmpH_Mass"    ,MassNBins,0,MassHistoUpperBound);
                    TH1D* tmpH_MassTOF  =  new TH1D("tmpH_MassTOF" ,"tmpH_MassTOF" ,MassNBins,0,MassHistoUpperBound);
                    TH1D* tmpH_MassComb =  new TH1D("tmpH_MassComb","tmpH_MassComb",MassNBins,0,MassHistoUpperBound);


                    double PE_A=RNG->Poisson(A);
                    double PE_B=RNG->Poisson(B);
                    double PE_C=RNG->Poisson(C);
                    double PE_E=RNG->Poisson(E);
                    double PE_F=RNG->Poisson(F);
                    double PE_G=RNG->Poisson(G);
                    double PE_P = 0;

                    if(E>0){        PE_P    = (PE_E>0 ? (PE_A*PE_F*PE_G)/(PE_E*PE_E) : 0);}
                    else if(A>0){   PE_P    = (PE_A>0 ? ((PE_C*PE_B)/PE_A) : 0);}

                    for(int i=0;i<Pred_EtaB_Proj_PE->GetNbinsX()+1;i++){Pred_EtaB_Proj_PE->SetBinContent(i,RNG->Poisson(Pred_EtaB_Proj->GetBinContent(i)) );}    Pred_EtaB_Proj_PE->Scale(1.0/Pred_EtaB_Proj_PE->Integral());
                    for(int i=0;i<Pred_EtaS_Proj_PE->GetNbinsX()+1;i++){Pred_EtaS_Proj_PE->SetBinContent(i,RNG->Poisson(Pred_EtaS_Proj->GetBinContent(i)) );}    Pred_EtaS_Proj_PE->Scale(1.0/Pred_EtaS_Proj_PE->Integral());
                    for(int i=0;i<Pred_EtaS2_Proj_PE->GetNbinsX()+1;i++){Pred_EtaS2_Proj_PE->SetBinContent(i,RNG->Poisson(Pred_EtaS2_Proj->GetBinContent(i)) );} Pred_EtaS2_Proj_PE->Scale(1.0/Pred_EtaS2_Proj_PE->Integral());

                    for(int i=0;i<Pred_EtaPWeighted_PE->GetNbinsX()+1;i++){
                        for(int j=0;j<=Pred_EtaPWeighted_PE->GetNbinsY()+1;j++){ // take the overflow as well
                            Pred_EtaPWeighted_PE->SetBinContent(i,j,RNG->Poisson(Pred_EtaPWeighted->GetBinContent(i,j)));
                        }
                    }
                    /*
                    // this part does the reweighing in terms of removing either high or
                    // low momentum tracks (currently set on removing low momenta first) to match the Eta distribution
                    // to remove high momenta, change this loop a few lines below
                    // for (int y = 0; y <= EtaPC->GetNbinsY()+1; y++)
                    // loop into
                    // for (int y = EtaPC->GetNbinsY()+1; --y;)
                    // if you opt for this reweighing, be sure to use EtaPC in the end for your momentum distribution
                    // instead of Pred_EtaPWeighted_PE
                    //
                        TH1D* AbsoluteWeight = Pred_EtaPWeighted_PE->ProjectionX("PWeighted");
                        Pred_EtaB_Proj_PE->Scale(AbsoluteWeight->GetMaximum()/Pred_EtaB_Proj_PE->GetMaximum());
                        Pred_EtaS_Proj_PE->Scale(AbsoluteWeight->GetMaximum()/Pred_EtaS_Proj_PE->GetMaximum());
                        for (int x = 0; x < AbsoluteWeight->GetNbinsX()+1; x++) {
                                AbsoluteWeight->SetBinContent(x, Pred_EtaB_Proj_PE->GetBinContent(x) - Pred_EtaS_Proj_PE->GetBinContent(x));
                                if (AbsoluteWeight->GetBinContent(x)<0) AbsoluteWeight->SetBinContent(x, 0);
                        }
                        TH2D* EtaPC = (TH2D*) Pred_EtaPWeighted_PE->Clone("EtaPC");
                        for (int x = 0; x < EtaPC->GetNbinsX()+1; x++){
                                if (AbsoluteWeight->GetBinContent(x)<=0) continue;
                                for (int y = 0; y <= EtaPC->GetNbinsY()+1; y++){
                                if (EtaPC->GetBinContent(x,y)>0){
                                    if (AbsoluteWeight->GetBinContent(x) >= EtaPC->GetBinContent(x,y)){
                                        AbsoluteWeight->SetBinContent(x, AbsoluteWeight->GetBinContent(x) - EtaPC->GetBinContent(x,y));
                                        EtaPC->SetBinContent(x, y, 0);
                                        if (AbsoluteWeight->GetBinContent(x)<=0) break;
                                    } else {
                                        EtaPC->SetBinContent(x, y, EtaPC->GetBinContent(x, y) - AbsoluteWeight->GetBinContent(x));
                                        break;
                                    }
                                }
                                }
                        }
                        delete AbsoluteWeight;
                    */
                    if (symmetrizeHistos) symmetrizeHisto(Pred_EtaB_Proj_PE,TypeMode);
                    if (symmetrizeHistos) symmetrizeHisto(Pred_EtaS_Proj_PE,TypeMode);
                    if (symmetrizeHistos) symmetrizeHisto(Pred_EtaS2_Proj_PE,TypeMode);
                    if (symmetrizeHistos) symmetrizeHisto(Pred_EtaPWeighted_PE,TypeMode);
                    Pred_EtaB_Proj_PE->Scale(1.0/Pred_EtaB_Proj_PE->Integral());
                    Pred_EtaS_Proj_PE->Scale(1.0/Pred_EtaS_Proj_PE->Integral());
                    if(TypeMode==2) Pred_EtaS2_Proj_PE->Scale(1.0/Pred_EtaS2_Proj_PE->Integral());// jesli normalizowac, to wszystkie


                    double WeightP = 0.0;
                    for(int x=0;x<=Pred_EtaPWeighted_PE->GetXaxis()->GetNbins();x++){
                        WeightP = 0.0;
                        //          reweigh C_Eta on B_Eta/A_Eta
                        if(Pred_EtaB_Proj_PE->GetBinContent(x)>0){
                            //mk pz    WeightP = 1.;
                            WeightP = Pred_EtaS_Proj_PE ->GetBinContent(x)/Pred_EtaB_Proj_PE->GetBinContent(x);
                            if(TypeMode==2)WeightP*= Pred_EtaS2_Proj_PE->GetBinContent(x)/Pred_EtaB_Proj_PE->GetBinContent(x);
                        }

                        /*
                        //          reweigh C_Eta on B_Eta/C_Eta -- not out of the box, in this case we don't need the    
                        //	    TH1D* Pred_EtaCRegion = (TH1D*) Pred_EtaPWeighted->ProjectionX ("Pred_EtaCRegion");
                                TH1D* Pred_EtaCRegion = (TH1D*) EtaPC->ProjectionX ("Pred_EtaCRegion", 0, EtaPC->GetNbinsY()+1);
                                Pred_EtaCRegion->Scale(1.0/Pred_EtaCRegion->Integral());
                                    if(Pred_EtaCRegion->GetBinContent(x)>0){
                                                    WeightP = Pred_EtaS_Proj_PE ->GetBinContent(x)/Pred_EtaCRegion->GetBinContent(x);
                                    if(TypeMode==2)WeightP*= Pred_EtaS2_Proj_PE->GetBinContent(x)/Pred_EtaCRegion->GetBinContent(x);
                                    }
                                delete Pred_EtaCRegion;
                        */
                        for(int y=0;y<=Pred_EtaPWeighted_PE->GetYaxis()->GetNbins()+1;y++){ // take the overflow as well
                            Pred_EtaPWeighted_PE->SetBinContent(x,y,Pred_EtaPWeighted_PE->GetBinContent(x,y)*WeightP);
                                //EtaPC->SetBinContent(x,y,EtaPC->GetBinContent(x,y)*WeightP);
                        }
                    }

                    // Plus side

                    //         TH1D* Pred_P_ProjPE = Pred_EtaPWeighted_PE /*EtaPC*/->ProjectionY("Pred_P_ProjPE", 0, 30);                                                        Pred_P_ProjPE->Scale(1.0/Pred_P_ProjPE->Integral(0, Pred_P_ProjPE->GetNbinsX()+1)); // include overflow


                    TH1D* Pred_P_ProjPE = Pred_EtaPWeighted_PE /*EtaPC*/->ProjectionY("Pred_P_ProjPE", 1, Pred_EtaPWeighted_PE->GetXaxis()->GetNbins());                                                        Pred_P_ProjPE->Scale(1.0/Pred_P_ProjPE->Integral(0, Pred_P_ProjPE->GetNbinsX()+1)); // include overflow
                    // Minus side
                    //         TH1D* Pred_P_ProjPE_M = Pred_EtaPWeighted_PE /*EtaPC*/->ProjectionY("Pred_P_ProjPE", 0, 30);                                                        Pred_P_ProjPE->Scale(1.0/Pred_P_ProjPE->Integral(0, Pred_P_ProjPE->GetNbinsX()+1)); // include overflow
                    for(int i=0;i<Pred_I_ProjPE->GetNbinsX()+1;i++){Pred_I_ProjPE->SetBinContent(i,RNG->Poisson(Pred_I_Proj->GetBinContent(i)) );}   Pred_I_ProjPE->Scale(1.0/Pred_I_ProjPE->Integral());
                    for(int i=0;i<Pred_T_ProjPE->GetNbinsX()+1;i++){Pred_T_ProjPE->SetBinContent(i,RNG->Poisson(Pred_T_Proj->GetBinContent(i)) );}   Pred_T_ProjPE->Scale(1.0/Pred_T_ProjPE->Integral());

                    //save the predP distribution
                    for(int x=0;x<=Pred_P_ProjPE->GetNbinsX()+1;x++){Pred_P->SetBinContent(CutIndex+1, x, Pred_P->GetBinContent(CutIndex+1, x) + Pred_P_ProjPE->GetBinContent(x) * PE_P);}; // overflow

                    double Proba, MI, MComb;//, MT=0, ProbaT=0;
                    for(int x=0;x<=Pred_P_ProjPE->GetNbinsX()+1;x++){    
                        if(Pred_P_ProjPE->GetBinContent(x)<=0.0){continue;}  
                        const double& p = (x<Pred_P_ProjPE->GetNbinsX()+1)?Pred_P_ProjPE->GetBinCenter(x):Pred_P_ProjPE->GetBinCenter(x+3); //overflow bin
                        for(int y=0;y<Pred_I_ProjPE->GetNbinsX()+1;y++){    if(Pred_I_ProjPE->GetBinContent(y)<=0.0){continue;}  const double& i = Pred_I_ProjPE->GetBinCenter(y);
                            Proba = Pred_P_ProjPE->GetBinContent(x) * Pred_I_ProjPE->GetBinContent(y);  if(Proba<=0 || isnan((float)Proba))continue;
                            //MI = GetMass(p,i, false);
                            MI = GetMass(p,i, dEdxK_MC, dEdxC_MC);
                            MComb = MI;
                            tmpH_Mass->Fill(MI,Proba);

                            //if(MI>500 && MI<800 && pe==0 && (CutIndex==3 || CutIndex==4)){printf("Index=%i p=%7.2f i=%7.2f M=%7.2f P=%7.5E( = %7.5E x %7.5E) --> %7.5E | %7.5E\n",CutIndex,p,i,MI,Proba,Pred_P_ProjPE->GetBinContent(x),Pred_I_ProjPE->GetBinContent(y),Pred_P_ProjPE->GetBinContent(x)*1,Pred_I_ProjPE->GetBinContent(y)*1);}

                            //commented part there is related to the prediction of the mass reconstructed from TOF.
                            //if(TypeMode==2){
                            //for(int z=0;z<Pred_T_ProjPE->GetNbinsX()+1;z++){   if(Pred_T_ProjPE->GetBinContent(z)<=0.0){continue;}   const double& t = Pred_T_ProjPE->GetBinCenter(z);
                            //   ProbaT = Proba * Pred_T_ProjPE->GetBinContent(z);  if(ProbaT<=0 || isnan(ProbaT))continue;
                            //   MT = GetTOFMass(p,t);
                            //   tmpH_MassTOF->Fill(MT,ProbaT);
                            //   MComb = GetMassFromBeta(p, (GetIBeta(i, false) + (1/t))*0.5 );        
                            //   tmpH_MassComb->Fill(MComb,ProbaT);
                            //}}else{
                            tmpH_MassComb->Fill(MComb,Proba);
                            //}
                        }
                    }

                    for(int x=0;x<tmpH_Mass->GetNbinsX()+1;x++){
                        Pred_Prof_Mass    ->SetBinContent(x, pe, tmpH_Mass    ->GetBinContent(x) * PE_P);
                        Pred_Prof_MassTOF ->SetBinContent(x, pe, tmpH_MassTOF ->GetBinContent(x) * PE_P);
                        Pred_Prof_MassComb->SetBinContent(x, pe, tmpH_MassComb->GetBinContent(x) * PE_P);
                        if(isnan((float)(tmpH_Mass    ->GetBinContent(x) * PE_P))){printf("%f x %f\n",tmpH_Mass    ->GetBinContent(x), P /*PE_P*/); fflush(stdout);exit(0);}
                    }
            
                    if (DirName.find("13TeV16G")==string::npos && is2016 && CutIndex == 4) {Pred_EtaPWeighted_PE->SaveAs("Pred_EtaPWeighted_PE.root");/* if (EtaPC) EtaPC->SaveAs("EtaPC.root");*/}
                    //         delete EtaPC;
                    delete Pred_P_ProjPE;
                    delete tmpH_Mass;
                    delete tmpH_MassTOF;
                    delete tmpH_MassComb;
                    //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ7
                    /*
                    if(E>0){continue;          //Prediction in Pt-Is-TOF space
                    }else if(A>0){continue;    //Prediction in Pt-Is     space
                    }else {printf("A==0&&E==0 no eta bins prediction possible, CutIndex= %i5\n",CutIndex);}
                    */

                    double EE=0, AA=0, HH=0, GG=0, CC=0, FF=0, BB=0, MOM=0, ICK=0, MassH=0;

                    //loop over eta bins (without overflows)
                    for (int ee=1;ee<PDF_A_Eta_2D->GetNbinsX()+1;ee++){ 
                        EE=RNG->Poisson(PDF_E_Eta_2D->GetBinContent(ee));
                        AA=RNG->Poisson(PDF_A_Eta_2D->GetBinContent(ee));
                        if(EE>0){         //Prediction in Pt-Is-TOF space for HA only where one has a Mass PDF already
                            //loop over p bins (including overflow)
                            for(int mm=1;    mm<PDF_H_EtaMass_2D->GetNbinsY()+2;mm++){ 
                                HH=RNG->Poisson(  PDF_H_EtaMass_2D->GetBinContent(ee,mm));  
                                MassH=((TAxis*)   PDF_H_EtaMass_2D->GetYaxis())->GetBinCenter(mm); // Mass
                                Pred_Prof_Mass_HA ->Fill(MassH,pe,HH*AA/EE);
                            }
                        }
                        //loop over p bins (without overflows)
                        for (int pp=1;pp<PDF_C_EtaP_2D->GetNbinsY()+1;pp++){ 
                            MOM=((TAxis*)   PDF_C_EtaP_2D->GetYaxis())->GetBinCenter(pp); // momentum
                            GG=RNG->Poisson(PDF_G_EtaP_2D->GetBinContent(ee,pp));
                            CC=RNG->Poisson(PDF_C_EtaP_2D->GetBinContent(ee,pp));             
                            for (int ii=0;ii<PDF_B_EtaICK_2D->GetNbinsY()+1;ii++){ //loop over Ias bins (without overflows)
                                ICK=((TAxis*)   PDF_B_EtaICK_2D->GetYaxis())->GetBinCenter(ii); // ICK
                                FF=RNG->Poisson(PDF_F_EtaICK_2D->GetBinContent(ee,ii));
                                BB=RNG->Poisson(PDF_B_EtaICK_2D->GetBinContent(ee,ii));             
                                if(E>0){         //Prediction in Pt-Is-TOF space
                                    if(EE>0&&ICK>0){
                                        Pred_Prof_Mass_GFA->Fill(MOM*sqrt(ICK),pe,GG*FF/EE*AA/EE);
                                        Pred_Prof_Mass_GB ->Fill(MOM*sqrt(ICK),pe,GG*BB/EE);
                                        Pred_Prof_Mass_CF ->Fill(MOM*sqrt(ICK),pe,CC*FF/EE);
                                    }
                                }else if(A>0){   //Prediction in Pt-Is     space
                                    if(AA>0&&ICK>0){
                                        Pred_Prof_Mass_CB ->Fill(MOM*sqrt(ICK),pe,CC*BB/AA);
                                    }
                                }
                            }
                        }// END loop over p bins
                    }// END loop over eta bins
                    //----------------------------------------------------------------------------------------------------Prediction in eta bins
	            }// END // loop over pseudo-experiments
                printf("\n");

                TH2D* MassesLooseCut = NULL;
                if (CutIndex == 4) MassesLooseCut = new TH2D ("MassesLooseCut", "MassesLooseCut", NPseudoExp, 0, NPseudoExp, MassNBins,0,MassHistoUpperBound);
                //       vector <double> masses;
                for(int x=0;x<Pred_Mass->GetNbinsY()+1;x++){
                    double Mean=0, MeanTOF=0, MeanComb=0, Median=0;
                    for(unsigned int pe=0;pe<NPseudoExp;pe++){
                        Mean     += Pred_Prof_Mass    ->GetBinContent(x, pe);
                        MeanTOF  += Pred_Prof_MassTOF ->GetBinContent(x, pe);
                        MeanComb += Pred_Prof_MassComb->GetBinContent(x, pe);
                        //	     masses.push_back(Pred_Prof_Mass    ->GetBinContent(x, pe));
                        if (MassesLooseCut) MassesLooseCut->SetBinContent (pe, x, Pred_Prof_Mass->GetBinContent (x, pe));
                    }
                    Mean/=NPseudoExp; MeanTOF/=NPseudoExp;  MeanComb/=NPseudoExp;
                    //	  sort (masses.begin(), masses.end());
                    //	  if (NPseudoExp%2==0) Median = masses[NPseudoExp/2];
                    //	  if (NPseudoExp%2==1) Median = masses[(NPseudoExp+1)/2];

                    double Err=0, ErrTOF=0, ErrComb=0, MAD=0;
                    for(unsigned int pe=0;pe<NPseudoExp;pe++){
                        Err     += pow(Mean     - Pred_Prof_Mass    ->GetBinContent(x, pe),2);
                        ErrTOF  += pow(MeanTOF  - Pred_Prof_MassTOF ->GetBinContent(x, pe),2);
                        ErrComb += pow(MeanComb - Pred_Prof_MassComb->GetBinContent(x, pe),2);
                        MAD     += fabs(Median - Pred_Prof_Mass->GetBinContent(x, pe));
                    }
                    Err=sqrt(Err/(NPseudoExp-1)); ErrTOF=sqrt(ErrTOF/(NPseudoExp-1));  ErrComb=sqrt(ErrComb/(NPseudoExp-1));
                    //	  MAD /= NPseudoExp;

                    Pred_Mass    ->SetBinContent(CutIndex+1,x,Mean    ); Pred_Mass      ->SetBinError(CutIndex+1,x,Err    );
                    //          Pred_Mass    ->SetBinContent(CutIndex+1,x,Median  ); Pred_Mass      ->SetBinError(CutIndex+1,x,MAD    );
                    Pred_MassTOF ->SetBinContent(CutIndex+1,x,MeanTOF ); Pred_MassTOF   ->SetBinError(CutIndex+1,x,ErrTOF );
                    Pred_MassComb->SetBinContent(CutIndex+1,x,MeanComb); Pred_MassComb  ->SetBinError(CutIndex+1,x,ErrComb);
                }
                if (MassesLooseCut && is2016 && DirName.find("13TeV16G")==string::npos){
                    MassesLooseCut->SaveAs("MassesLooseCut.root");Pred_EtaPWeighted_PE->SaveAs("Pred_EtaP_Weighted_PE.root"); delete MassesLooseCut;
                }
                delete Pred_EtaB_Proj_PE;
                delete Pred_EtaS_Proj_PE;
                delete Pred_EtaS2_Proj_PE;

                delete Pred_Prof_Mass;
                delete Pred_Prof_MassTOF;
                delete Pred_Prof_MassComb;
                delete Pred_EtaPWeighted_PE;
                delete Pred_I_ProjPE;
                delete Pred_T_ProjPE;

                delete Pred_I_Proj;
                delete Pred_T_Proj;
                delete Pred_EtaB_Proj;
                delete Pred_EtaS_Proj;
                delete Pred_EtaS2_Proj;
                delete Pred_EtaPWeighted;
                //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ8
                delete PDF_E_Eta_2D;
                delete PDF_A_Eta_2D;
                delete PDF_H_EtaMass_2D;
                delete PDF_G_EtaP_2D;
                delete PDF_C_EtaP_2D;
                delete PDF_F_EtaICK_2D;
                delete PDF_B_EtaICK_2D;


                for(int x=0;x<Pred_Mass->GetNbinsY()+1;x++){
                    double MeanGFA=0, MeanGB=0, MeanCF=0, MeanHA=0, MeanCB=0;
                    for(unsigned int pe=0;pe<NPseudoExp;pe++){
                        MeanGFA     += Pred_Prof_Mass_GFA    ->GetBinContent(x, pe);
                        MeanGB      += Pred_Prof_Mass_GB     ->GetBinContent(x, pe);
                        MeanCF      += Pred_Prof_Mass_CF     ->GetBinContent(x, pe);
                        MeanHA      += Pred_Prof_Mass_HA     ->GetBinContent(x, pe);
                        MeanCB      += Pred_Prof_Mass_CB     ->GetBinContent(x, pe);
                    }
                    MeanGFA/=NPseudoExp; MeanGB/=NPseudoExp; MeanCF/=NPseudoExp; MeanHA/=NPseudoExp; MeanCB/=NPseudoExp;

                    double ERR=0, ErrGFA=0, ErrGB=0, ErrCF=0, ErrHA=0, ErrCB=0;
                    for(unsigned int pe=0;pe<NPseudoExp;pe++){
                        ErrGFA += pow(MeanGFA - Pred_Prof_Mass_GFA->GetBinContent(x, pe),2);
                        ErrGB  += pow(MeanGB  - Pred_Prof_Mass_GB ->GetBinContent(x, pe),2);
                        ERR    += pow(MeanGB  - Pred_Prof_Mass_GFA->GetBinContent(x, pe),2);
                        ErrCF  += pow(MeanCF  - Pred_Prof_Mass_CF ->GetBinContent(x, pe),2);
                        ERR    += pow(MeanCF  - Pred_Prof_Mass_GFA->GetBinContent(x, pe),2);
                        ErrHA  += pow(MeanHA  - Pred_Prof_Mass_HA ->GetBinContent(x, pe),2);
                        ERR    += pow(MeanHA  - Pred_Prof_Mass_GFA->GetBinContent(x, pe),2);
                        ErrCB  += pow(MeanCB  - Pred_Prof_Mass_CB ->GetBinContent(x, pe),2);
                    }
                    ErrGFA=sqrt(ErrGFA/(NPseudoExp-1)); 
                    ErrGB =sqrt(ErrGB /(NPseudoExp-1)); 
                    ErrCF =sqrt(ErrCF /(NPseudoExp-1));
                    ErrHA =sqrt(ErrHA /(NPseudoExp-1));
                    ERR   =sqrt(ERR   /(3*NPseudoExp)  +ErrGFA*ErrGFA);
                    ErrCB =sqrt(ErrCB /(NPseudoExp-1));
            
                    Pred_Mass_GFA    ->SetBinContent(CutIndex+1,x,MeanGFA); Pred_Mass_GFA ->SetBinError(CutIndex+1,x,ErrGFA);
                    Pred_Mass_S      ->SetBinContent(CutIndex+1,x,MeanGFA); Pred_Mass_S   ->SetBinError(CutIndex+1,x,ERR   );
                    Pred_Mass_GB     ->SetBinContent(CutIndex+1,x,MeanGB ); Pred_Mass_GB  ->SetBinError(CutIndex+1,x,ErrGB );
                    Pred_Mass_CF     ->SetBinContent(CutIndex+1,x,MeanCF ); Pred_Mass_CF  ->SetBinError(CutIndex+1,x,ErrCF );
                    Pred_Mass_HA     ->SetBinContent(CutIndex+1,x,MeanHA ); Pred_Mass_HA  ->SetBinError(CutIndex+1,x,ErrHA );
                    Pred_Mass_CB     ->SetBinContent(CutIndex+1,x,MeanCB ); Pred_Mass_CB  ->SetBinError(CutIndex+1,x,ErrCB );
                }
                //----------------------------------------------------------------------------------------------------Prediction in eta bins
                delete Pred_Prof_Mass_GFA;
                delete Pred_Prof_Mass_GB;
                delete Pred_Prof_Mass_CF;
                delete Pred_Prof_Mass_HA;
                delete Pred_Prof_Mass_CB; 
    
        
            }// END loop over CutIndex
            //scale it down by the number of PseudoExperiment to get right normalization
            Pred_P->Scale(1.0/NPseudoExp);


            //save histogram to file
            Pred_P       ->Write();
            H_P          ->Write();
            Pred_Mass    ->Write();
            Pred_MassTOF ->Write();
            Pred_MassComb->Write();
            H_P_Coll->Write();
            H_P_Cosmic->Write();
            //----------------------------------------------------------------------------------------------------Prediction in eta bins PZ9
            Pred_Mass_GFA->Write();
            Pred_Mass_S  ->Write();
            Pred_Mass_GB ->Write();
            Pred_Mass_CF ->Write();
            Pred_Mass_HA ->Write();
            Pred_Mass_CB ->Write();
            //----------------------------------------------------------------------------------------------------Prediction in eta bins

            if(TypeMode==3) {
                for(int i=0; i<PredBins; i++) {
                    H_P_Binned[i]->Write();
                }
            }
            //directory->Delete("H_P;1");

            //////////////////////////////////////////////////     DUMP USEFUL INFORMATION
            //WAIT//FILE* pFile = fopen((InputPattern+"/Info_"+directory->GetName()+Suffix+".txt").c_str(),"w");
            FILE* pFile = fopen("Info_test.txt","w");
            for(unsigned int CutIndex=0;CutIndex<(unsigned int)HCuts_Pt->GetXaxis()->GetNbins();CutIndex++){
                const double& A=H_A->GetBinContent(CutIndex+1);
                const double& B=H_B->GetBinContent(CutIndex+1);
                const double& C=H_C->GetBinContent(CutIndex+1);
                const double& D=H_D->GetBinContent(CutIndex+1);
                const double& E=H_E->GetBinContent(CutIndex+1);
                const double& F=H_F->GetBinContent(CutIndex+1);
                const double& G=H_G->GetBinContent(CutIndex+1);
                const double& H=H_H->GetBinContent(CutIndex+1);
                fprintf(pFile  ,"CutIndex=%4i --> (Pt>%6.2f I>%6.3f TOF>%6.3f) Ndata=%+6.2E  NPred=%6.3E+-%6.3E (=%6.3E+-%6.3E + %6.3E+-%6.3E) <--> A=%6.2E B=%6.2E C=%6.2E D=%6.2E E=%6.2E F=%6.2E G=%6.2E H=%6.2E\n",CutIndex,HCuts_Pt ->GetBinContent(CutIndex+1), HCuts_I  ->GetBinContent(CutIndex+1), HCuts_TOF->GetBinContent(CutIndex+1), D,H_P->GetBinContent(CutIndex+1),H_P->GetBinError(CutIndex+1), H_P_Coll->GetBinContent(CutIndex+1),H_P_Coll->GetBinError(CutIndex+1), H_P_Cosmic->GetBinContent(CutIndex+1),H_P_Cosmic->GetBinError(CutIndex+1), A, B, C, D, E, F, G, H);
            }
            fprintf(pFile,"--------------------\n");
            fclose(pFile);      

        }//End loop on sub directories
    }//End of loop on two predictions
}

int main(int argc, char* argv[]) {

    string usage = "Usage: ./RunBackgroundPrediction --inputFiles <file.txt> [--mode <value>]\n";
    usage       += "Or   : ./RunBackgroundPrediction --inputFiles <file1 file2 ...> [--mode <value>]";

    vector<string> input;
    int typeMode = 0;

    ArgumentParser parser(argc,argv);

    if( parser.findOption("-h") ){
        cout << usage << endl;
        return 0;
    }
    if( parser.findOption("--inputFiles") ) parser.getArgument("--inputFiles", input);
    else {
        cout << usage << endl;
        return 0;
    }
    if( parser.findOption("--mode") ) parser.getArgument("--mode", typeMode);


    cout << "======================" << endl;
    cout << " BackgroundPrediction " << endl;
    cout << "======================\n" << endl;

    setTDRStyle();
    /*gStyle->SetPadTopMargin   (0.06);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadRightMargin (0.16);
    gStyle->SetPadLeftMargin  (0.14);
    gStyle->SetTitleSize(0.04, "XYZ");
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.45);
    gStyle->SetPalette(1);
    gStyle->SetNdivisions(505);*/

    TBenchmark clock;
    clock.Start("BackgroundPrediction");

    vector<string> inputFiles;
    if (endsWith(input[0],".txt")){
        ifstream fin(input[0], std::ios::binary);
        if (not fin.is_open()) {
            cout << "Failed to open " << input[0] << endl;
            return 0;
        }
        string rootfile;
        while (fin >> rootfile)
            inputFiles.push_back(rootfile);
    }
    else
        inputFiles = input;

    for(auto const &inputFile : inputFiles){
    
        cout << "opening file " << inputFile << endl;
        //TFile* inFile = TFile::Open(inputFile.c_str(),"UPDATE");
        TFile* tfile = TFile::Open(inputFile.c_str(),"UPDATE");
        //TFile* tfile = TFile::Open(inputFile.c_str());
        if(not tfile){
            cout << "Failed to open " << inputFile << endl;
            return 0;
        }
        Analysis_Step2_BackgroundPrediction(tfile, typeMode);
    }

    cout << "" << endl;
    clock.Show("BackgroundPrediction");
    cout << "" << endl;

    return 0;
}
