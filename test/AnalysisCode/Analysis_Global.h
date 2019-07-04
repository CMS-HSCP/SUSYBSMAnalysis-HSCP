// Original Author:  Loic Quertenmont

#ifndef HSCP_ANALYSIS_GLOBAL
#define HSCP_ANALYSIS_GLOBAL

//Include widely used in all the codes
#include <string>
#include <vector>
#include <unordered_map>
#include <fstream>
#include "TROOT.h"
#include "TCanvas.h"
#include "TChain.h"
#include "TCutG.h" 
#include "TDCacheFile.h"
#include "TDirectory.h"
#include "TF1.h"
#include "TFile.h"
#include "TGraph.h"
#include "TGraphAsymmErrors.h"
#include "TGraphErrors.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TLegend.h"
#include "TLegendEntry.h"
#include "TMath.h"
#include "TMultiGraph.h"
#include "TObject.h"
#include "TPaveText.h"
#include "TProfile.h"
#include "TProfile2D.h"
#include "TRandom3.h"
#include "TTree.h"

//double IntegratedLuminosity13TeV               = 72.63; //pb
//double IntegratedLuminosity13TeV               = 84.557; //pb
double IntegratedLuminosity13TeV15             = 2490.518; //2439.264; //pb
//double IntegratedLuminosity13TeV16             = 4002.946; // pb
double preTrackingChangeL1IntLumi              = 29679.982; // pb
double IntegratedLuminosity13TeV16             = 34228.911;  // pb -> not from brilcalc, but from the email, brilcalc numpy error
double IntegratedLuminosity13TeV16G            = 16554.454;  // pb from brilcalc
//double IntegratedLuminosity13TeV16G            = 18148.982;  // pb from brilcalc
double IntegratedLuminosity13TeV16PreG         = 17674.457;

double               SQRTS                     = 13;
double               SQRTS15                   = 1315;
double               SQRTS16                   = 1316;
double               SQRTS1615                 = 131615;
int                  RunningPeriods            = 1;
double IntegratedLuminosity                    = IntegratedLuminosity13TeV15; //pb
double IntegratedLuminosityBeforeTriggerChange =    0; //pb

bool correctMassCut                            = true;


string IntegratedLuminosityFromE(double SQRTS_){
  char LumiText[1024];

  if(SQRTS_==13 || SQRTS_==1315)
     sprintf(LumiText,"%1.1f fb^{-1} (%1.0f TeV)", 0.001*IntegratedLuminosity13TeV15, 13.0);
  else if(SQRTS_==1316.0)
     sprintf(LumiText,"%1.1f fb^{-1} (%1.0f TeV)", 0.001*IntegratedLuminosity13TeV16PreG, 13.0);
  else if(SQRTS_==13167.0)
     sprintf(LumiText,"%1.1f fb^{-1} (%1.0f TeV)", 0.001*IntegratedLuminosity13TeV16G, 13.0);
  else if(SQRTS_==131615 || SQRTS_==131516)
     sprintf(LumiText,"2016: %1.1f fb^{-1}   2015: %1.1f fb^{-1} (%1.0f TeV)", 0.001*IntegratedLuminosity13TeV16, 0.001*IntegratedLuminosity13TeV15, 13.0);
  else if(SQRTS_==131677)
     sprintf(LumiText,"2016 - Pre-G: %1.1f fb^{-1}, Post-G: %1.1f fb^{-1} (%1.0f TeV)", 0.001*IntegratedLuminosity13TeV16PreG, 0.001*IntegratedLuminosity13TeV16G, 13.0);
  else sprintf(LumiText, "unknown energy and int. lumi");
//  if(SQRTS_==78 || SQRTS_==87)sprintf(LumiText,"#sqrt{s} = %1.0f TeV, L = %1.1f fb^{-1}   #sqrt{s} = %1.0f TeV, L = %1.1f fb^{-1}", 7.0, 0.001*IntegratedLuminosity7TeV,8.0, 0.001*IntegratedLuminosity8TeV);
  return LumiText;
}

// Type of the analysis
int		   TypeMode         = 0; //0 = Tracker-Only analysis (used in 2010 and 2011 papers)
					 //1 = Tracker+Muon analysis (used in 2010 paper)
                                         //2 = Tracker+TOF  analysis (used in 2011 paper)
                                         //3 = TOF Only     analysis (to be used in 2012 paper)
                                         //4 = Q>1          analysis (to be used in 2012 paper)
                                         //5 = Q<1          analysis (to be used in 2012 paper)
                                         //? do not hesitate to define your own --> TOF-Only, mCHAMPs, fractional charge

// directory where to find the EDM files --> check the function at the end of this file, to see how it is defined interactively
std::string BaseDirectory = "undefined... Did you call InitBaseDirectory() ? --> ";


// binning for the pT, mass, and IP distributions
double             PtHistoUpperBound   = 1200;
double             MassHistoUpperBound = 3000;
int                MassNBins           = 300;
double             IPbound             = 1.0;

// Thresholds for candidate preselection --> note that some of the followings can be replaced at the beginning of Analysis_Step1_EventLoop function
double             GlobalMaxEta       =   2.1;    // cut on inner tracker track eta
double             GlobalMaxV3D       =   99999;  //0.50;   // cut on 3D distance (cm) to closest vertex
double             GlobalMaxDZ        =   0.50;   // cut on 1D distance (cm) to closest vertex in "Z" direction
double             GlobalMaxDXY       =   0.50;   // cut on 2D distance (cm) to closest vertex in "R" direction
double             GlobalMaxChi2      =   5.0;    // cut on Track maximal Chi2/NDF
int                GlobalMinQual      =   2;      // cut on track quality (2 meaning HighPurity tracks)
int                FixedQual          =   2;      // cut on track quality (2 meaning HighPurity tracks) -- select ONLY high purity tracks, other cuts are irrelevant
unsigned int       GlobalMinNOH       =   8;//7AMSB;      // cut on number of (valid) track pixel+strip hits 
int                GlobalMinNOPH      =   2;      // cut on number of (valid) track pixel hits 
double             GlobalMinFOVH      =   0.8;//0.0AMSB;    // cut on fraction of valid track hits
unsigned int       GlobalMaxNOMHTillLast = 99999;//1AMSB;     // cut on the number of missing hits from IP till last hit (excluding hits behind the last hit)
double             GlobalMinFOVHTillLast =-99999;//0.85AMSB;   // cut on the fraction of valid hits divided by total expected hits until the last one
unsigned int       GlobalMinNOM       =   6;//7AMSB;      // cut on number of dEdx hits (generally equal to #strip+#pixel-#ClusterCleaned hits, but this depend on estimator used) -- should be more stringent -- 10 should be used
double             GlobalMinNDOF      =   8;      // cut on number of     DegreeOfFreedom used for muon TOF measurement
double             GlobalMinNDOFDT    =   6;      // cut on number of DT  DegreeOfFreedom used for muon TOF measurement
double             GlobalMinNDOFCSC   =   6;      // cut on number of CSC DegreeOfFreedom used for muon TOF measurement
double             GlobalMaxTOFErr    =   0.15;//0.07;   // cut on error on muon TOF measurement
double             GlobalMaxPterr     = 1;//0.50;//0.25;   // cut on error on track pT measurement 
double             GlobalMaxTIsol     =  50;      // cut on tracker isolation (SumPt)
double             GlobalMaxRelTIsol  =  9999999; // cut on relative tracker isolation (SumPt/Pt)
double             GlobalMaxEIsol     =  0.30;    // cut on calorimeter isolation (E/P)
double             GlobalMinPt        =  55.00;   // cut on pT    at PRE-SELECTION
double             GlobalMinIs        =   0.0;    // cut on dEdxS at PRE-SELECTION (dEdxS is generally a  discriminator)
double             GlobalMinIm        =   0.0;    // cut on dEdxM at PRE-SELECTION (dEdxM is generally an estimator    )
double             GlobalMinTOF       =   1.0;    // cut on TOF   at PRE-SELECTION
const int          MaxPredBins        =   6;      // The maximum number of different bins prediction is done in for any of the analyses (defines array size)
int                PredBins           =   0;      //  How many different bins the prediction is split in for analysis being run, sets how many histograms are actually initialized.
int                EtaBins            =  60;      //  How many bins we use for the background prediction method in Eta -- impacts background prediction method -- histograms with the name of the form "Pred_Eta*" in Analysis_PlotStructure.h

// dEdx related variables, Name of dEdx estimator/discriminator to be used for selection (dEdxS) and for mass reconstruction (dEdxM)
// as well as the range for the dEdx variable and K/C constant for mass reconstruction
std::string        dEdxS_Label     = "dedxASmi";
double             dEdxS_UpLim     = 1.0;
std::string        dEdxS_Legend    = "I_{as}";
std::string        dEdxM_Label     = "dedxHarm2";
double             dEdxM_UpLim     = 30.0;
std::string        dEdxM_Legend    = "I_{h} (MeV/cm)";
//double             dEdxK_Data      = 2.779;
//double             dEdxC_Data      = 2.879;
//double             dEdxK_MC        = 2.683;
//double             dEdxC_MC        = 2.453;
// below are the new constants computed with harm2_SP_in for both MC and Data USED FOR THE PAS
//double             dEdxK_Data      = 2.535;
//double             dEdxC_Data      = 3.339;
//double             dEdxK_MC        = 2.535;
//double             dEdxC_MC        = 3.339;
// below are the new constants computed with hhybrid2_SP_in for both MC and Data USED FOR THE PAPER
double             dEdxK_Data15    = 2.535;
double             dEdxC_Data15    = 3.339;
double             dEdxK_MC15      = 2.535;
double             dEdxC_MC15      = 3.339;

double             dEdxK_Data16    = 2.270; //updated with pion fit (march 2019)
double             dEdxC_Data16    = 3.475; //updated  with pion fit (march 2019)
double             dEdxK_MC16      = 2.935;
double             dEdxC_MC16      = 3.197;

double             dEdxK_Data      = dEdxK_Data16;
double             dEdxC_Data      = dEdxC_Data16;
double             dEdxK_MC        = dEdxK_MC16;
double             dEdxC_MC        = dEdxC_MC16;


//below is for harm2_SP where 20% of the lower dEdx hits are drop
//double             dEdxK_Data      = 2.40;
//double             dEdxC_Data      = 3.45;
//double             dEdxK_MC        = 2.78;
//double             dEdxC_MC        = 3.35;

// TOF object to be used for combined, DT and CSC TOF measurement
std::string        TOF_Label       = "combined";
std::string        TOFdt_Label     = "dt";
std::string        TOFcsc_Label    = "csc";

//Variables used in the TOF only HSCP search
float              DTRegion        =   0.9;  //Define the dividing line between DT and 
float              CSCRegion       =   0.9;  //CSC regions of CMS
float              CosmicMinDz     =   70.;  //Min dz displacement to be tagged as cosmic muon
float              CosmicMaxDz     =   120.; //Max dz displacement for cosmic tagged tracks
double             minSegEtaSep    =  0.1;   //Minimum eta separation between SA track and muon segment on opposite side of detector
const int          DzRegions       =  6;     //Number of different Dz side regions used to make cosmic background prediction
int                minMuStations   =  2;


//for initializing PileupReweighting utility.
const   float TrueDist2012_f[100] = {6.53749e-07 ,1.73877e-06 ,4.7972e-06 ,1.57721e-05 ,2.97761e-05 ,0.000162201 ,0.000931952 ,0.00272619 ,0.0063166 ,0.0128901 ,0.0229009 ,0.0355021 ,0.045888 ,0.051916 ,0.0555598 ,0.0580188 ,0.059286 ,0.0596022 ,0.059318 ,0.0584214 ,0.0570249 ,0.0553875 ,0.0535731 ,0.0512788 ,0.0480472 ,0.0436582 ,0.0382936 ,0.0323507 ,0.0262419 ,0.0203719 ,0.0151159 ,0.0107239 ,0.00727108 ,0.00470101 ,0.00288906 ,0.00168398 ,0.000931041 ,0.000489695 ,0.000246416 ,0.00011959 ,5.65558e-05 ,2.63977e-05 ,1.23499e-05 ,5.89242e-06 ,2.91502e-06 ,1.51247e-06 ,8.25545e-07 ,4.71584e-07 ,2.79203e-07 ,1.69571e-07 ,1.04727e-07 ,6.53264e-08 ,4.09387e-08 ,2.56621e-08 ,1.60305e-08 ,9.94739e-09 ,6.11516e-09 ,3.71611e-09 ,2.22842e-09 ,1.3169e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};  // MB xsec = 69.3mb
const   float TrueDist2012_XSecShiftUp_f[100] = {6.53749e-07 ,1.73877e-06 ,4.7972e-06 ,1.57721e-05 ,2.97761e-05 ,0.000162201 ,0.000931952 ,0.00272619 ,0.0063166 ,0.0128901 ,0.0229009 ,0.0355021 ,0.045888 ,0.051916 ,0.0555598 ,0.0580188 ,0.059286 ,0.0596022 ,0.059318 ,0.0584214 ,0.0570249 ,0.0553875 ,0.0535731 ,0.0512788 ,0.0480472 ,0.0436582 ,0.0382936 ,0.0323507 ,0.0262419 ,0.0203719 ,0.0151159 ,0.0107239 ,0.00727108 ,0.00470101 ,0.00288906 ,0.00168398 ,0.000931041 ,0.000489695 ,0.000246416 ,0.00011959 ,5.65558e-05 ,2.63977e-05 ,1.23499e-05 ,5.89242e-06 ,2.91502e-06 ,1.51247e-06 ,8.25545e-07 ,4.71584e-07 ,2.79203e-07 ,1.69571e-07 ,1.04727e-07 ,6.53264e-08 ,4.09387e-08 ,2.56621e-08 ,1.60305e-08 ,9.94739e-09 ,6.11516e-09 ,3.71611e-09 ,2.22842e-09 ,1.3169e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // MB xsec = 73.5mb; observed in Z-->MuMu see https://twiki.cern.ch/twiki/bin/viewauth/CMS/PileupJSONFileforData#Calculating_Your_Pileup_Distribu
const   float TrueDist2012_XSecShiftDown_f[100] = {6.53749e-07 ,1.73877e-06 ,4.7972e-06 ,1.57721e-05 ,2.97761e-05 ,0.000162201 ,0.000931952 ,0.00272619 ,0.0063166 ,0.0128901 ,0.0229009 ,0.0355021 ,0.045888 ,0.051916 ,0.0555598 ,0.0580188 ,0.059286 ,0.0596022 ,0.059318 ,0.0584214 ,0.0570249 ,0.0553875 ,0.0535731 ,0.0512788 ,0.0480472 ,0.0436582 ,0.0382936 ,0.0323507 ,0.0262419 ,0.0203719 ,0.0151159 ,0.0107239 ,0.00727108 ,0.00470101 ,0.00288906 ,0.00168398 ,0.000931041 ,0.000489695 ,0.000246416 ,0.00011959 ,5.65558e-05 ,2.63977e-05 ,1.23499e-05 ,5.89242e-06 ,2.91502e-06 ,1.51247e-06 ,8.25545e-07 ,4.71584e-07 ,2.79203e-07 ,1.69571e-07 ,1.04727e-07 ,6.53264e-08 ,4.09387e-08 ,2.56621e-08 ,1.60305e-08 ,9.94739e-09 ,6.11516e-09 ,3.71611e-09 ,2.22842e-09 ,1.3169e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // MB xsec = 65.835mb

const float TrueDist2015_f[100] = {141899.9,    622421.5,    818662.3,    1253138,    2114164,    4937264,    1.627961e+07,    6.475337e+07,    1.817159e+08,    3.220406e+08,    4.16444e+08,    4.291912e+08,    3.663307e+08,    2.557795e+08,    1.447021e+08,    6.727821e+07,    2.736262e+07,    1.139197e+07,    5668287,    3025054,    1401679,    511811,    146735.4,    35294.85,    8269.747,    2235.301,    721.3224,    258.8397,    97.26937,    36.8714,    13.72746,    4.931709,    1.692407,    0.5518936,    0.1706013,    0.0499358,    0.01383391,    0.003626617,    0.0008996063,    0.0002111484,    4.689276e-05,    9.85392e-06,    1.959294e-06,    3.686198e-07,    6.562482e-08,    1.105342e-08,    1.762478e-09,    2.614969e-10,    4.768003e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const   float TrueDist2015_XSecShiftUp_f[100] = {   119011.7,    575187.4,    742505.6,    1119922,    1752064,    3603643,    1.025708e+07,    3.766355e+07,    1.209274e+08,    2.47339e+08,    3.605043e+08,    4.119143e+08,    3.926434e+08,    3.147836e+08,    2.092462e+08,    1.150021e+08,    5.334319e+07,    2.238585e+07,    9895584,    5174317,    2843878,    1368343,    530329.1,    164616.8,    43068.57,    10705.86,    2954.413,    964.7733,    354.4973,    138.0159,    54.82292,    21.6373,    8.325977,    3.086652,    1.095163,    0.3706287,    0.1194443,    0.03663022,    0.0106862,    0.002965238,    0.0007825775,    0.0001964351,    4.689582e-05,    1.064815e-05,    2.29954e-06,    4.723236e-07,    9.227752e-08,    1.714345e-08,    3.031843e-09,    5.050254e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const   float TrueDist2015_XSecShiftDown_f[100] = {166787.7,    674049,    913830.3,    1425367,    2638444,    7155207,    2.778629e+07,    1.094527e+08,    2.587231e+08,    3.982389e+08,    4.560547e+08,    4.197129e+08,    3.116313e+08,    1.836509e+08,    8.643619e+07,    3.42708e+07,    1.339142e+07,    6282790,    3245101,    1445026,    494664.3,    130033.8,    28603.11,    6332.219,    1680.587,    533.0118,    185.3784,    66.66396,    23.86613,    8.286534,    2.746725,    0.8620616,    0.2551907,    0.07113181,    0.01865672,    0.004603238,    0.001068332,    0.0002332118,    4.788436e-05,    9.247787e-06,    1.679908e-06,    2.870404e-07,    4.613098e-08,    6.975087e-09,    9.940954e-10,    1.331527e-10,    1.59604e-11, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

const float TrueDist2016_f[100] = {2822.19, 250452.4, 937691, 1664697, 2619751, 3565840, 4627719, 9533987, 3.106812e+07, 7.442141e+07, 1.54961e+08, 2.738765e+08, 4.205078e+08, 5.877119e+08, 7.642573e+08, 9.326678e+08, 1.066947e+09, 1.159318e+09, 1.213587e+09, 1.23433e+09, 1.227176e+09, 1.194125e+09, 1.135918e+09, 1.057543e+09, 9.64454e+08, 8.594009e+08, 7.458684e+08, 6.304864e+08, 5.210056e+08, 4.228937e+08, 3.382319e+08, 2.668304e+08, 2.076012e+08, 1.591548e+08, 1.199555e+08, 8.852149e+07, 6.36179e+07, 4.427985e+07, 2.970336e+07, 1.912839e+07, 1.179261e+07, 6947340, 3907544, 2097994, 1075998, 528117.8, 249112, 114029.5, 51803.16, 24509.65, 13094.47, 8515.386, 6722.897, 5996.101, 5643.242, 5402.081, 5179.803, 4946.978, 4696.771, 4430.086, 4150.377, 3861.862, 3568.872, 3275.582, 2985.862, 2703.188, 2430.58, 2170.564, 1925.156, 1695.865, 1483.716, 1289.274, 1112.694, 953.7694, 811.9872, 686.584, 576.6037, 480.951, 398.4413, 327.8447, 267.9245, 217.4689, 175.3163, 140.3744, 111.6333, 88.17391, 69.17134, 53.89545, 41.70782, 32.05694, 24.47178, 18.55443, 13.97231, 10.45025, 7.762874, 5.727378, 4.196876, 3.05445, 2.207883, 1.585094};
const float TrueDist2016_XSecShiftUp_f[100] = {1634.393, 206470.1, 816751.6, 1397443, 2371275, 3112905, 4069668, 6158024, 1.883994e+07, 4.87039e+07, 1.034972e+08, 1.94571e+08, 3.135822e+08, 4.548117e+08, 6.102477e+08, 7.706301e+08, 9.180892e+08, 1.032709e+09, 1.11146e+09, 1.157717e+09, 1.175342e+09, 1.169009e+09, 1.140347e+09, 1.089751e+09, 1.021046e+09, 9.390486e+08, 8.462858e+08, 7.451003e+08, 6.402245e+08, 5.380952e+08, 4.442198e+08, 3.614785e+08, 2.904354e+08, 2.304752e+08, 1.80571e+08, 1.395211e+08, 1.060552e+08, 7.899706e+07, 5.738697e+07, 4.046513e+07, 2.758198e+07, 1.811387e+07, 1.143421e+07, 6926789, 4023701, 2240806, 1196924, 614035, 303433, 145372.2, 68504.14, 32757.41, 16848.39, 10054.66, 7247.247, 6095.602, 5589.093, 5311.823, 5101.729, 4899.81, 4686.984, 4458.986, 4216.868, 3963.516, 3702.395, 3437.072, 3171, 2907.407, 2649.224, 2399.033, 2159.037, 1931.043, 1716.457, 1516.298, 1331.215, 1161.514, 1007.198, 867.9995, 743.4292, 632.8133, 535.3378, 450.0875, 376.0823, 312.3101, 257.7546, 211.4191, 172.3453, 139.6278, 112.4246, 89.96389, 71.54708, 56.54997, 44.4211, 34.67869, 26.90618, 20.74709, 15.89929, 12.10915, 9.165679, 6.894935};
const float TrueDist2016_XSecShiftDown_f[100] = {4788.162,    308025.7,    1065583,    1992484,    2934566,    4077637,    5602066,    1.653722e+07,    4.975603e+07,    1.151307e+08,    2.284721e+08,    3.785913e+08,    5.56859e+08,    7.504048e+08,    9.420001e+08,    1.100246e+09,    1.209735e+09,    1.274151e+09,    1.298987e+09,    1.291077e+09,    1.252941e+09,    1.185818e+09,    1.096138e+09,    9.900185e+08,    8.706727e+08,    7.433731e+08,    6.170774e+08,    5.006005e+08,    3.988648e+08,    3.128956e+08,    2.417499e+08,    1.83846e+08,    1.373459e+08,    1.003848e+08,    7.135496e+07,    4.90081e+07,    3.233369e+07,    2.039632e+07,    1.226051e+07,    7008285,    3805522,    1962832,    962574.4,    449980.6,    201758.3,    88065.66,    38768.71,    18513.46,    10601.05,    7627.918,    6508.273,    6029.02,    5742.798,    5495.859,    5240.517,    4965.349,    4670.767,    4360.884,    4040.863,    3716.013,    3391.426,    3071.78,    2761.23,    2463.325,    2180.97,    1916.41,    1671.243,    1446.456,    1242.472,    1059.217,    896.1916,    752.5507,    627.177,    518.7576,    425.8526,    346.9564,    280.5506,    225.1483,    179.3275,    141.7576,    111.2157,    86.59779,    66.92181,    51.32733,    39.07053,    29.51681,    22.13139,    16.46897,    12.16305,    8.915316,    6.485573,    4.682493,    3.355232,    2.386077,    1.684078,    1.179656,    0.8200945,    0.5658315,    0.3874581,    0.2633155};

const float TrueDist2016G_f[100] = {236340, 582615.3, 1253111, 1339639, 1592782, 2230405, 2269792, 3492157, 4744094, 5830225, 2.516408e+07, 8.630452e+07, 1.72147e+08, 2.686555e+08, 3.78458e+08, 5.078014e+08, 6.281231e+08, 7.121934e+08, 7.592498e+08, 7.874632e+08, 8.231413e+08, 8.690824e+08, 9.016346e+08, 9.082807e+08, 8.991826e+08, 8.846946e+08, 8.652907e+08, 8.378568e+08, 7.996478e+08, 7.488523e+08, 6.855186e+08, 6.121396e+08, 5.331206e+08, 4.535068e+08, 3.777277e+08, 3.087596e+08, 2.479366e+08, 1.953358e+08, 1.504502e+08, 1.127188e+08, 8.171348e+07, 5.704874e+07, 3.821431e+07, 2.449266e+07, 1.499169e+07, 8752220, 4869676, 2581174, 1303232, 626888.5, 287437.7, 125747.1, 52571, 21055.85, 8109.3, 3018.746, 1093.606, 388.7934, 136.9462, 48.28288, 17.21885, 6.277593, 2.36426, 0.9280961, 0.3815795, 0.1640383, 0.0731353, 0.03341637, 0.01545884, 0.007168489, 0.003308073, 0.001512026, 0.0006825101, 0.0003037189, 0.0001331095, 5.742019e-05, 2.437142e-05, 1.017559e-05, 4.178426e-06, 1.687318e-06, 6.69966e-07, 2.61531e-07, 1.003673e-07, 3.786621e-08, 1.403427e-08, 5.114176e-09, 1.832959e-09, 6.43311e-10, 2.231323e-10, 7.38758e-11, 2.543465e-11, 7.922507e-12, 3.568623e-12, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const float TrueDist2016G_XSecShiftUp_f[100] = {230877.5, 434415.7, 1252878, 1212087, 1438685, 2002748, 2146839, 2654803, 4396017, 4411504, 1.044585e+07, 4.524728e+07, 1.132154e+08, 1.953775e+08, 2.856046e+08, 3.905565e+08, 5.094412e+08, 6.136256e+08, 6.848468e+08, 7.249238e+08, 7.501299e+08, 7.823475e+08, 8.241488e+08, 8.559718e+08, 8.651781e+08, 8.586469e+08, 8.463808e+08, 8.303999e+08, 8.08359e+08, 7.778129e+08, 7.368952e+08, 6.848652e+08, 6.227833e+08, 5.535541e+08, 4.811713e+08, 4.097129e+08, 3.424965e+08, 2.815883e+08, 2.278039e+08, 1.810706e+08, 1.409402e+08, 1.069573e+08, 7.877462e+07, 5.607867e+07, 3.846059e+07, 2.534948e+07, 1.602883e+07, 9711714, 5634017, 3128076, 1661860, 844876.3, 411154.4, 191645.8, 85651.46, 36765.16, 15194.71, 6068.155, 2353.076, 891.5058, 332.47, 123.0692, 45.62064, 17.08917, 6.528392, 2.566539, 1.046674, 0.4449771, 0.1971711, 0.09053342, 0.04266113, 0.02041044, 0.009820438, 0.00471729, 0.002250773, 0.001063207, 0.0004962109, 0.0002285305, 0.0001037849, 4.645635e-05, 2.049045e-05, 8.903809e-06, 3.811072e-06, 1.606658e-06, 6.670418e-07, 2.727128e-07, 1.097888e-07, 4.350653e-08, 1.697828e-08, 6.522799e-09, 2.462457e-09, 9.200519e-10, 3.368865e-10, 1.193577e-10, 4.027529e-11, 1.625383e-11, 6.760731e-12, 9.151568e-14, 0.0, 0.0};
const float TrueDist2016G_XSecShiftDown_f[100] = {243824.2, 775789.8, 1242178, 1465862, 1819384, 2397521, 2580611, 4527904, 4882343, 1.192443e+07, 5.723073e+07, 1.438728e+08, 2.470935e+08, 3.627031e+08, 5.015209e+08, 6.401359e+08, 7.402594e+08, 7.962849e+08, 8.281701e+08, 8.677372e+08, 9.183114e+08, 9.513652e+08, 9.550462e+08, 9.432325e+08, 9.259881e+08, 9.020212e+08, 8.675139e+08, 8.194436e+08, 7.564222e+08, 6.799177e+08, 5.943686e+08, 5.058126e+08, 4.200957e+08, 3.414938e+08, 2.72159e+08, 2.124611e+08, 1.618521e+08, 1.19636e+08, 8.527745e+07, 5.829893e+07, 3.806016e+07, 2.365462e+07, 1.396653e+07, 7823502, 4154428, 2090545, 996877.5, 450622.8, 193254.1, 78745.93, 30560.74, 11338.69, 4043.815, 1396.595, 471.4881, 157.3049, 52.48465, 17.72755, 6.137182, 2.204245, 0.8296434, 0.3287436, 0.1366453, 0.05893023, 0.02599405, 0.01157097, 0.005144383, 0.002268424, 0.0009877694, 0.0004236625, 0.0001787262, 7.409739e-05, 3.017569e-05, 1.206752e-05, 4.737972e-06, 1.826052e-06, 6.907386e-07, 2.564187e-07, 9.340797e-08, 3.337595e-08, 1.170396e-08, 4.030145e-09, 1.354385e-09, 4.487284e-10, 1.454326e-10, 4.567288e-11, 1.371745e-11, 5.755874e-12, 9.151568e-14, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};

const   float Pileup_MC_Fall11[100]= {1.45346E-01, 6.42802E-02, 6.95255E-02, 6.96747E-02, 6.92955E-02, 6.84997E-02, 6.69528E-02, 6.45515E-02, 6.09865E-02, 5.63323E-02, 5.07322E-02, 4.44681E-02, 3.79205E-02, 3.15131E-02, 2.54220E-02, 2.00184E-02, 1.53776E-02, 1.15387E-02, 8.47608E-03, 6.08715E-03, 4.28255E-03, 2.97185E-03, 2.01918E-03, 1.34490E-03, 8.81587E-04, 5.69954E-04, 3.61493E-04, 2.28692E-04, 1.40791E-04, 8.44606E-05, 5.10204E-05, 3.07802E-05, 1.81401E-05, 1.00201E-05, 5.80004E-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const   float Pileup_MC_Summer2012[100] = { 2.560E-06, 5.239E-06, 1.420E-05, 5.005E-05, 1.001E-04, 2.705E-04, 1.999E-03, 6.097E-03, 1.046E-02, 1.383E-02, 1.685E-02, 2.055E-02, 2.572E-02, 3.262E-02, 4.121E-02, 4.977E-02, 5.539E-02, 5.725E-02, 5.607E-02, 5.312E-02, 5.008E-02, 4.763E-02, 4.558E-02, 4.363E-02, 4.159E-02, 3.933E-02, 3.681E-02, 3.406E-02, 3.116E-02, 2.818E-02, 2.519E-02, 2.226E-02, 1.946E-02, 1.682E-02, 1.437E-02, 1.215E-02, 1.016E-02, 8.400E-03, 6.873E-03, 5.564E-03, 4.457E-03, 3.533E-03, 2.772E-03, 2.154E-03, 1.656E-03, 1.261E-03, 9.513E-04, 7.107E-04, 5.259E-04, 3.856E-04, 2.801E-04, 2.017E-04, 1.439E-04, 1.017E-04, 7.126E-05, 4.948E-05, 3.405E-05, 2.322E-05, 1.570E-05, 5.005E-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
const   float Pileup_MC_Startup2015_25ns[100] = {4.8551E-07, 1.74806E-06, 3.30868E-06, 1.62972E-05, 4.95667E-05, 0.000606966, 0.003307249, 0.010340741, 0.022852296, 0.041948781, 0.058609363, 0.067475755, 0.072817826, 0.075931405, 0.076782504, 0.076202319, 0.074502547, 0.072355135, 0.069642102, 0.064920999, 0.05725576, 0.047289348, 0.036528446, 0.026376131, 0.017806872, 0.011249422, 0.006643385, 0.003662904, 0.001899681, 0.00095614, 0.00050028, 0.000297353, 0.000208717, 0.000165856, 0.000139974, 0.000120481, 0.000103826, 8.88868E-05, 7.53323E-05, 6.30863E-05, 5.21356E-05, 4.24754E-05, 3.40876E-05, 2.69282E-05, 2.09267E-05, 1.5989E-05, 4.8551E-06, 2.42755E-06, 4.8551E-07, 2.42755E-07, 1.21378E-07, 4.8551E-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};


// function used to define interactively the directory containing the EDM files
// you are please to add the line for your case and not touch the line of the other users
void InitBaseDirectory(){  
   char* analystTmp=getenv("USER");
   char* hostTmp   =getenv("HOSTNAME");
   char* remoteStoragePath =getenv("REMOTESTORAGEPATH");
   if(!hostTmp||!analystTmp)return;
   string analyst(analystTmp);
   string host   (hostTmp);
  
   if(getenv("PWD")!=NULL)host+=string(" PWD=") + getenv("PWD");

   // BaseDirectory is first set for the AAA protocol ...
   if(remoteStoragePath!=NULL){
      string RemoteStorageServer = "cms-xrd-global.cern.ch"; // default value
      if (getenv("REMOTESTORAGESERVER")!=NULL) RemoteStorageServer = string(getenv("REMOTESTORAGESERVER"));
      BaseDirectory = string("root://")+RemoteStorageServer+"/"+string(remoteStoragePath);
      printf("Accessing remote files using BaseDirectory = %s\n.", BaseDirectory.c_str());
   // if we give no path to remote storage site, BaseDirectory is then set via hostname
   }else if(host.find("ucl.ac.be")!=std::string::npos){
      //BaseDirectory = "/storage/data/cms/users/quertenmont/HSCP/CMSSW_4_2_8/12_08_16/"; //for run1
      BaseDirectory = "/storage/data/cms/store/user/jozobec/HSCP2016/"; //for run2
   }else if(host.find("cern.ch")!=std::string::npos){
      //BaseDirectory = "rfio:/castor/cern.ch/user/r/rybinska/HSCPEDMFiles/";
      //BaseDirectory = "root://eoscms//eos/cms/store/cmst3/user/querten/12_08_30_HSCP_EDMFiles/";    //for run1
      BaseDirectory = "root://eoscms//eos/cms/store/cmst3/user/querten/15_03_25_HSCP_Run2EDMFiles/";  //for run2
   }else if(host.find("fnal.gov")!=std::string::npos){
     BaseDirectory = "dcache:/pnfs/cms/WAX/11/store/user/lpchscp/2012HSCPEDMFiles/"; //for run1
   }else if(host.find("ingrid-ui")!=std::string::npos){
      BaseDirectory = "/storage/data/cms/store/user/jozobec/HSCP2016/"; // run2 ingrid
   }else if(host.find(".cis.gov.pl")!=std::string::npos){
      BaseDirectory = "root://se.cis.gov.pl//cms/store/user/fruboes/HSCP/15_03_25_HSCP_Run2EDMFiles/"; // run2 Swierk
   }else{
      BaseDirectory = "/storage/data/cms/store/user/jozobec/HSCP2016/";  //for run1
      printf("YOUR MACHINE (%s) IS NOT KNOWN --> please add your machine to the 'InitBaseDirectory' function of 'Analysis_Global.h'\n", host.c_str());
      printf("HOST=%s  USER=%s\n",host.c_str(), analyst.c_str());
      printf("In the mean time, the directory containing the HSCP EDM file is assumed to be %s\n",BaseDirectory.c_str());
   }

   // BaseDirectory is defined a function of the username
//   if(analyst.find("querten")!=std::string::npos && host.find("ucl.ac.be")!=std::string::npos){
//      BaseDirectory = "/storage/data/cms/users/quertenmont/HSCP/CMSSW_4_2_3/11_11_01/";
//   }   
}


#endif
