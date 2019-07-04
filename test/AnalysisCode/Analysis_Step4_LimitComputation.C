// Original Author:  Loic Quertenmont

#include "Analysis_Global.h"
#include "Analysis_CommonFunction.h"
#include "Analysis_PlotFunction.h"
#include "Analysis_Samples.h"
#include "tdrstyle.C"

#include "TGraphAsymmErrors.h"

using namespace std;

class stAllInfo{
   public:
   double Mass, MassMean, MassSigma, MassCut;
   double XSec_Th, XSec_Err, XSec_Exp, XSec_ExpUp, XSec_ExpDown, XSec_Exp2Up, XSec_Exp2Down, XSec_Obs;
   double  Eff, Eff_SYSTP, Eff_SYSTI, Eff_SYSTM, Eff_SYSTHUp, Eff_SYSTHDown, Eff_SYSTT, Eff_SYSTPU, TotalUnc;
   double  EffE, EffE_SYSTP, EffE_SYSTI, EffE_SYSTM, EffE_SYSTHUp, EffE_SYSTHDown, EffE_SYSTT, EffE_SYSTPU;
   double Significance; double XSec_5Sigma;
   double Index, WP_Pt, WP_I, WP_TOF;
   float  NData, NPred, NPredErr, NSign;
   double LInt;

   stAllInfo(string path=""){
      //Default Values
      Mass          = 0;      MassMean      = 0;      MassSigma     = 0;      MassCut       = 0;
      Index         = 0;      WP_Pt         = 0;      WP_I          = 0;      WP_TOF        = 0;
      XSec_Th       = 0;      XSec_Err      = 0;      XSec_Exp      = 1E50;   XSec_ExpUp    = 1E50;   XSec_ExpDown  = 1E50;    XSec_Exp2Up   = 1E50;    XSec_Exp2Down = 1E50;    XSec_Obs    = 1E50;
      Significance  = 1E50;   XSec_5Sigma   = 1E50;
      Eff           = 0;      Eff_SYSTP     = 0;      Eff_SYSTI     = 0;      Eff_SYSTM     = 0;      Eff_SYSTHUp     = 0;   Eff_SYSTHDown  = 0;  Eff_SYSTT     = 0;   Eff_SYSTPU  = 0;
      EffE          = 0;      EffE_SYSTP    = 0;      EffE_SYSTI    = 0;      EffE_SYSTM    = 0;      EffE_SYSTHUp    = 0;   EffE_SYSTHDown = 0;  EffE_SYSTT    = 0;   EffE_SYSTPU = 0;
      NData         = 0;      NPred         = 0;      NPredErr      = 0;      NSign         = 0;      LInt          = 0;
      if(path=="")return;
      FILE* pFile = fopen(path.c_str(),"r");
      if(!pFile){return; printf("Can't open %s\n",path.c_str()); return;}
      fscanf(pFile,"Mass         : %lf\n",&Mass);
      fscanf(pFile,"MassMean     : %lf\n",&MassMean);
      fscanf(pFile,"MassSigma    : %lf\n",&MassSigma);
      fscanf(pFile,"MassCut      : %lf\n",&MassCut);
      fscanf(pFile,"Index        : %lf\n",&Index);
      fscanf(pFile,"WP_Pt        : %lf\n",&WP_Pt);
      fscanf(pFile,"WP_I         : %lf\n",&WP_I);
      fscanf(pFile,"WP_TOF       : %lf\n",&WP_TOF);
      fscanf(pFile,"Eff          : %lf +- %lf\n",&Eff       , &EffE );
      fscanf(pFile,"Eff_SystP    : %lf +- %lf\n",&Eff_SYSTP , &EffE_SYSTP );
      fscanf(pFile,"Eff_SystI    : %lf +- %lf\n",&Eff_SYSTI , &EffE_SYSTI );
      fscanf(pFile,"Eff_SystM    : %lf +- %lf\n",&Eff_SYSTM , &EffE_SYSTM );
      fscanf(pFile,"Eff_SystHUp  : %lf +- %lf\n",&Eff_SYSTHUp , &EffE_SYSTHUp );
      fscanf(pFile,"Eff_SystHDown: %lf +- %lf\n",&Eff_SYSTHDown , &EffE_SYSTHDown );
      fscanf(pFile,"Eff_SystT    : %lf +- %lf\n",&Eff_SYSTT , &EffE_SYSTT );
      fscanf(pFile,"Eff_SystPU   : %lf +- %lf\n",&Eff_SYSTPU, &EffE_SYSTPU);
      fscanf(pFile,"TotalUnc     : %lf\n",&TotalUnc);
      fscanf(pFile,"Signif       : %lf\n",&Significance);
      fscanf(pFile,"XSec_Th      : %lE\n",&XSec_Th);
      fscanf(pFile,"XSec_Exp     : %lE\n",&XSec_Exp);
      fscanf(pFile,"XSec_ExpUp   : %lE\n",&XSec_ExpUp);
      fscanf(pFile,"XSec_ExpDown : %lE\n",&XSec_ExpDown);
      fscanf(pFile,"XSec_Exp2Up  : %lE\n",&XSec_Exp2Up);
      fscanf(pFile,"XSec_Exp2Down: %lE\n",&XSec_Exp2Down);
      fscanf(pFile,"XSec_Obs     : %lE\n",&XSec_Obs);
      fscanf(pFile,"NData        : %E\n" ,&NData);
      fscanf(pFile,"NPred        : %E\n" ,&NPred);
      fscanf(pFile,"NPredErr     : %E\n" ,&NPredErr);
      fscanf(pFile,"NSign        : %E\n" ,&NSign);
      fscanf(pFile,"LInt         : %lf\n",&LInt);
      fscanf(pFile,"XSec_5Sigma  : %lf\n",&XSec_5Sigma);
      fclose(pFile);
   }

   void Save(string path=""){
      FILE* pFile = fopen(path.c_str(),"w");
      if(!pFile)printf("Can't open file : %s\n",path.c_str());
      fprintf(pFile,"Mass         : %f\n",Mass);
      fprintf(pFile,"MassMean     : %f\n",MassMean);
      fprintf(pFile,"MassSigma    : %f\n",MassSigma);
      fprintf(pFile,"MassCut      : %f\n",MassCut);
      fprintf(pFile,"Index        : %f\n",Index);
      fprintf(pFile,"WP_Pt        : %f\n",WP_Pt);
      fprintf(pFile,"WP_I         : %f\n",WP_I);
      fprintf(pFile,"WP_TOF       : %f\n",WP_TOF);
      fprintf(pFile,"Eff          : %f +- %f\n",Eff       , EffE);
      fprintf(pFile,"Eff_SystP    : %f +- %f\n",Eff_SYSTP , EffE_SYSTP);
      fprintf(pFile,"Eff_SystI    : %f +- %f\n",Eff_SYSTI , EffE_SYSTI);
      fprintf(pFile,"Eff_SystM    : %f +- %f\n",Eff_SYSTM , EffE_SYSTM);
      fprintf(pFile,"Eff_SystHUp  : %f +- %f\n",Eff_SYSTHUp , EffE_SYSTHUp);
      fprintf(pFile,"Eff_SystHDown: %f +- %f\n",Eff_SYSTHDown , EffE_SYSTHDown );
      fprintf(pFile,"Eff_SystT    : %f +- %f\n",Eff_SYSTT , EffE_SYSTT);
      fprintf(pFile,"Eff_SystPU   : %f +- %f\n",Eff_SYSTPU, EffE_SYSTPU);
      fprintf(pFile,"TotalUnc     : %f\n",TotalUnc);
      fprintf(pFile,"Signif       : %f\n",Significance);
      fprintf(pFile,"XSec_Th      : %.12E\n",XSec_Th);
      fprintf(pFile,"XSec_Exp     : %.12E\n",XSec_Exp);
      fprintf(pFile,"XSec_ExpUp   : %.12E\n",XSec_ExpUp);
      fprintf(pFile,"XSec_ExpDown : %.12E\n",XSec_ExpDown);
      fprintf(pFile,"XSec_Exp2Up  : %.12E\n",XSec_Exp2Up);
      fprintf(pFile,"XSec_Exp2Down: %.12E\n",XSec_Exp2Down);
      fprintf(pFile,"XSec_Obs     : %.12E\n",XSec_Obs);     
      fprintf(pFile,"NData        : %+6.2E\n",NData);
      fprintf(pFile,"NPred        : %+6.2E\n",NPred);
      fprintf(pFile,"NPredErr     : %+6.2E\n",NPredErr);
      fprintf(pFile,"NSign        : %+6.2E\n",NSign);
      fprintf(pFile,"LInt         : %f\n",LInt);
      fprintf(pFile,"XSec_5Sigma  : %f\n",XSec_5Sigma);
      fclose(pFile); 
   }
};


string EXCLUSIONDIR = "EXCLUSION";

//Background prediction rescale and uncertainty
double RescaleFactor = 1.0;
double RescaleError  = 0.2;

//final Plot y-axis range
double PlotMinScale = 2.0e-7;
double PlotMaxScale = 1000;

//Easy flag to skip running time consuming Cls expected limits. True runs the limit, false does not
bool FullExpLimit=true;

void Optimize(string InputPattern, string Data, string signal, bool shape, bool cutFromFile, int* OptimCutIndex=nullptr);
double GetSignalMeanHSCPPerEvent(string InputPattern, unsigned int CutIndex, double MinRange, double MaxRange);
TGraph* MakePlot(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, int XSectionType, std::vector<stSample>& modelSamples, double& LInt);
TGraph* CheckSignalUncertainty(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, std::vector<stSample>& modelSample);
void DrawModelLimitWithBand(string InputPattern);
void DrawRatioBands(string InputPattern);
void printSummary(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, std::vector<stSample>& modelSamples);
void printSummaryPaper(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, std::vector<stSample>& modelSamples);
string toLatex(double value);

void makeDataCard(string outpath, string rootPath, string ChannelName, string SignalName, double Obs, double Pred, double PredRelErr, double Sign, double SignalStat, double SignalUnc, bool Shape);
void saveHistoForLimit(TH1* histo, string Name, string Id);
void saveVariationHistoForLimit(TH1* histo, TH1* vardown, string Name, string variationName);
void testShapeBasedAnalysis(string InputPattern, string signal);
double computeSignificance(string datacard, bool expected, string& signal, string massStr, float Strength);
bool runCombine(bool fastOptimization, bool getXsection, bool getSignificance, string& InputPattern, string& signal, unsigned int CutIndex, bool Shape, bool Temporary, stAllInfo& result, TH1* MassData, TH1* MassPred, TH1* MassSign, TH1* MassSignP, TH1* MassSignI, TH1* MassSignM, TH1* MassSignHUp, TH1* MassSignHDown, TH1* MassSignT, TH1* MassSignPU);
bool Combine(string InputPattern, string signal7, string signal8, int* OptimCutIndex=nullptr);
bool useSample(int TypeMode, string sample);

double MinRange = 0;
double MaxRange = 9999;
int    CurrentSampleIndex;

std::vector<stSample> samples;
std::vector<std::string> modelVector;
std::map<std::string, std::vector<stSample> > modelMap;

string SHAPESTRING="";

//void Analysis_Step4_LimitComputation(){
//   return;
//}

void Analysis_Step4_LimitComputation(string MODE="COMPILE", string InputPattern="", string signal=""){
  setTDRStyle();
  gStyle->SetPadTopMargin   (0.05);
  gStyle->SetPadBottomMargin(0.11);
  gStyle->SetPadRightMargin (0.05);
  gStyle->SetPadLeftMargin  (0.15);
  gStyle->SetTitleSize(0.05, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.4);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505,"X");
  gStyle->SetNdivisions(550,"Y");
  //gStyle->SetTextFont(43);

   printf("Test %s\n", MODE.c_str());
   if(MODE=="COMPILE")return;

   string Data;
   if(MODE.find("SHAPE")!=string::npos){SHAPESTRING="SHAPE";}else{SHAPESTRING="";}

/*   if(MODE.find("COMPUTELIMIT")!=string::npos || MODE.find("OPTIMIZE")!=string::npos){
      if(signal.find("7TeV")!=string::npos){ Data = "Data7TeV";  SQRTS= 7.0; EXCLUSIONDIR+= "7TeV"; }
      if(signal.find("8TeV")!=string::npos){ Data = "Data8TeV";  SQRTS= 8.0; EXCLUSIONDIR+= "8TeV"; }
      if(signal.find("13TeV")!=string::npos){Data = "Data13TeV"; SQRTS=13.0; EXCLUSIONDIR+="13TeV"; }
      printf("EXCLUSIONDIR = %s\nData = %s\n",EXCLUSIONDIR.c_str(), Data.c_str());  

      if(MODE.find("COMPUTELIMIT")!=string::npos){Optimize(InputPattern, Data, signal, SHAPESTRING!="", true);      return;}
      if(MODE.find("OPTIMIZE")!=string::npos){    Optimize(InputPattern, Data, signal, SHAPESTRING!="", false);     return;} //testShapeBasedAnalysis(InputPattern,signal);  //use the second part if you want to run shape based analyssi on optimal point form c&c      
   }
*/
   if(MODE.find("COMPUTELIMIT2015")!=string::npos || MODE.find("OPTIMIZE2015")!=string::npos){
      if(signal.find("13TeV")!=string::npos){Data = "Data13TeV"; SQRTS=1315.0; EXCLUSIONDIR+="13TeV15"; }
      printf("EXCLUSIONDIR = %s\nData = %s\n",EXCLUSIONDIR.c_str(), Data.c_str());  

      if(MODE.find("COMPUTELIMIT")!=string::npos){Optimize(InputPattern, Data, signal, SHAPESTRING!="", true);      return;}
      if(MODE.find("OPTIMIZE")!=string::npos){    Optimize(InputPattern, Data, signal, SHAPESTRING!="", false);     return;} //testShapeBasedAnalysis(InputPattern,signal);  //use the second part if you want to run shape based analyssi on optimal point form c&c      
   }

   if(MODE.find("COMPUTELIMIT2016G")!=string::npos || MODE.find("OPTIMIZE2016G")!=string::npos){
      if(signal.find("13TeV16G")==string::npos) return;
      Data = "Data13TeV16G"; SQRTS=13167.0; EXCLUSIONDIR+="13TeV16G";
      printf("EXCLUSIONDIR = %s\nData = %s\n",EXCLUSIONDIR.c_str(), Data.c_str());  

      if(MODE.find("COMPUTELIMIT")!=string::npos){Optimize(InputPattern, Data, signal, SHAPESTRING!="", true);      return;}
      if(MODE.find("OPTIMIZE")!=string::npos){    Optimize(InputPattern, Data, signal, SHAPESTRING!="", false);     return;} //testShapeBasedAnalysis(InputPattern,signal);  //use the second part if you want to run shape based analyssi on optimal point form c&c      
   }

   if(MODE.find("16G") == string::npos && (MODE.find("COMPUTELIMIT2016")!=string::npos || MODE.find("OPTIMIZE2016")!=string::npos)){
      if(signal.find("13TeV16")==string::npos || signal.find("13TeV16G")!=string::npos) return;
      Data = "Data13TeV16"; SQRTS=1316.0; EXCLUSIONDIR+="13TeV16";
      printf("EXCLUSIONDIR = %s\nData = %s\n",EXCLUSIONDIR.c_str(), Data.c_str());  

      if(MODE.find("COMPUTELIMIT")!=string::npos){Optimize(InputPattern, Data, signal, SHAPESTRING!="", true);      return;}
      if(MODE.find("OPTIMIZE")!=string::npos){    Optimize(InputPattern, Data, signal, SHAPESTRING!="", false);     return;} //testShapeBasedAnalysis(InputPattern,signal);  //use the second part if you want to run shape based analyssi on optimal point form c&c      
   }

   if(MODE.find("COMBINE_Run1")!=string::npos){
      printf("COMBINE!!!\n");

      string signal7TeV = signal; if(signal7TeV.find("_8TeV")!=string::npos) signal7TeV = signal7TeV.replace(signal7TeV.find("_8TeV"),5, "_7TeV");
      string signal8TeV = signal; if(signal8TeV.find("_7TeV")!=string::npos) signal8TeV = signal8TeV.replace(signal8TeV.find("_7TeV"),5, "_8TeV");

      string EXCLUSIONDIR_SAVE = EXCLUSIONDIR;
      //2011 Limits
      Data = "Data7TeV"; SQRTS=7.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"7TeV";
      Optimize(InputPattern, Data, signal7TeV, SHAPESTRING!="", true);

      //2012 Limits
      Data = "Data8TeV"; SQRTS=8.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"8TeV";
      Optimize(InputPattern, Data, signal8TeV, SHAPESTRING!="", true);

      //Combined Limits
      EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"COMB";  SQRTS=78.0;
      Combine(InputPattern, signal7TeV, signal8TeV);
      return;
   }

   if(MODE.find("COMBINE2016")!=string::npos){
      printf("COMBINE!!!\n");

      string signal13TeV16  = ReplacePartOfString(signal, "13TeV16", "13TeV") + "W13TeV16";
      string signal13TeV16G = ReplacePartOfString(signal, "13TeV16G", "13TeV") + "W13TeV16G";

      string EXCLUSIONDIR_SAVE = EXCLUSIONDIR;

      //2016 PreG Limits
      signed int OptCutIndex = -1;
      printf("2016 pre-G Data ...\n");
      Data = "Data13TeV16"; SQRTS=1316.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"13TeV16";
      Optimize(InputPattern, Data, signal, SHAPESTRING!="", true, &OptCutIndex);

      //2016 PostG Limits
      printf("2016G post-G Data ...\n");
      Data = "Data13TeV16G"; SQRTS=13167.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"13TeV16G";
      Optimize(InputPattern, Data, ReplacePartOfString(signal, "13TeV16", "13TeV16G"), SHAPESTRING!="", true, &OptCutIndex);

      //Combined Limits
      printf("Combining ...\n");
      EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"COMB2016";  SQRTS=131667.0;
      Combine(InputPattern, signal13TeV16, signal13TeV16G, &OptCutIndex);
      return;
   }

   if(MODE.find("COMBINE_Run2")!=string::npos){
      printf("COMBINE!!!\n");

      string signal13TeV15 = ReplacePartOfString(signal, "13TeV16", "13TeV") + "W13TeV15";
      string signal13TeV16 = signal + "W13TeV16";

      string EXCLUSIONDIR_SAVE = EXCLUSIONDIR;

      //2015 Limits
      printf("2015 Data ...\n");
      Data = "Data13TeV"; SQRTS=1315.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"13TeV15";
      Optimize(InputPattern, Data, ReplacePartOfString(signal, "13TeV16", "13TeV"), SHAPESTRING!="", true);

      //2016 Limits
      printf("2016 pre G Data ...\n");
      Data = "Data13TeV16"; SQRTS=1316.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"13TeV16";
      Optimize(InputPattern, Data, signal, SHAPESTRING!="", true);

      //2016G Limits
      printf("2016G Data ...\n");
      Data = "Data13TeV16G"; SQRTS=13167.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"13TeV16G";
      Optimize(InputPattern, Data, ReplacePartOfString(signal, "13TeV16", "13TeV"), SHAPESTRING!="", true);

      //Combined Limits
      printf("Combining ...\n");
      EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"COMB";  SQRTS=131615.0;
      Combine(InputPattern, signal13TeV15, signal13TeV16);
      return;
   }
   if     (MODE.find("7TeV"    )!=string::npos){Data = "Data7TeV"    ; SQRTS=7.0    ; EXCLUSIONDIR+="7TeV"    ; }
   else if(MODE.find("8TeV"    )!=string::npos){Data = "Data8TeV"    ; SQRTS=8.0    ; EXCLUSIONDIR+="8TeV"    ; }
   else if(MODE.find("13TeV15" )!=string::npos){Data = "Data13TeV"   ; SQRTS=1315.0 ; EXCLUSIONDIR+="13TeV15" ; }
   else if(MODE.find("13TeV16G")!=string::npos){Data = "Data13TeV16G"; SQRTS=13167.0; EXCLUSIONDIR+="13TeV16G"; }
   else if(MODE.find("13TeV16" )!=string::npos){Data = "Data13TeV16" ; SQRTS=1316.0 ; EXCLUSIONDIR+="13TeV16" ; }
   else if(MODE.find("13TeV"   )!=string::npos){Data = "Data13TeV"   ; SQRTS=13.0   ; EXCLUSIONDIR+="13TeV"   ; }

   // #FIXME JOZE
   string TkPattern  = "Results/Type0/";
   string MuPattern  = "Results/Type2/";
//   string MOPattern  = "Results/Type3/";  //Disabling this analyis
//   string HQPattern  = "Results/Type4/";  //Disabling this analyis
//   string LQPattern  = "Results/Type5/";  //Disabling this analyis

   bool Combine = (MODE.find("COMB")!=string::npos);
   if (MODE.find("Run1")!=string::npos){ EXCLUSIONDIR+="COMB"    ; SQRTS=78.0    ; }
   if (MODE.find("Run2")!=string::npos){ EXCLUSIONDIR+="COMB"    ; SQRTS=131615.0; }
   if (MODE.find("2016")!=string::npos){ EXCLUSIONDIR+="COMB2016"; SQRTS=131677.0; }
//   if(Combine) {PlotMinScale=1E-6; PlotMaxScale=300;} // FIXME JOZZE -- for 13TeV we don't need this

   string outpath = string("Results/"+SHAPESTRING+EXCLUSIONDIR+"/");
   MakeDirectories(outpath);

   //determine the list of models that are considered
   GetSampleDefinition(samples);

   if(SQRTS!=78.0 && SQRTS!=131615.0 && SQRTS!=131516.0 && SQRTS!=131677.0) keepOnlySamplesAtSQRTS(samples, SQRTS);

   for(unsigned int s=0;s<samples.size();s++){
      if(samples[s].Type!=2)continue;
      //printf("Name-->Model >>  %30s --> %s\n",samples[s].Name.c_str(), samples[s].ModelName().c_str());
      
      if(SQRTS== 7.0    && samples[s].Name.find( "_7TeV")==string::npos){continue;}
      if(SQRTS== 8.0    && samples[s].Name.find( "_8TeV")==string::npos){continue;}
      if(SQRTS==13.0    && samples[s].Name.find("_13TeV")==string::npos){continue;}
      if(SQRTS==1315.0  && samples[s].Name.find("_13TeV")==string::npos){continue;}
      if(SQRTS==1316.0  && (samples[s].Name.find("_13TeV16")==string::npos || samples[s].Name.find("13TeV16G")!=string::npos)){continue;}
      if(SQRTS==13167.0 && samples[s].Name.find("_13TeV16G")==string::npos){continue;}
//      if(SQRTS==78.0){if(samples[s].Name.find("_7TeV")==string::npos){continue;}else{samples[s].Name.replace(samples[s].Name.find("_7TeV"),5, ""); } }
      if(SQRTS==78.0){if(samples[s].Name.find("_8TeV")==string::npos){continue;}else{samples[s].Name.replace(samples[s].Name.find("_8TeV"),5, ""); } }
      if(SQRTS==131615.0 || SQRTS==131516.0){
         if  (samples[s].Name.find("_13TeV")==string::npos){continue;}
         else{
            samples[s].Name = ReplacePartOfString (samples[s].Name,"_13TeV16" , "");
            samples[s].Name = ReplacePartOfString (samples[s].Name,"_13TeV15" , "");
            samples[s].Name = ReplacePartOfString (samples[s].Name,"_13TeV"   , "");
         }
      } else if (SQRTS==131677){
//              if (samples[s].Name.find("13TeV16")==std::npos) continue;
         if (samples[s].Name.find("13TeV16G")==string::npos) continue;
         samples[s].Name = ReplacePartOfString (samples[s].Name,"_13TeV16G", "");
      }
        
      if (!samples[s].ModelName().empty()){
         modelMap[samples[s].ModelName()].push_back(samples[s]);
         if(modelMap[samples[s].ModelName()].size()==1) modelVector.push_back(samples[s].ModelName());
      }
   }
   printf("EXCLUSIONDIR = %s\nData = %s\n",EXCLUSIONDIR.c_str(), Data.c_str());  

   //unti we have all the samples at both 7 and 8TeV, add the 7TeV models
//   if(SQRTS== 8.0){
//      for(unsigned int s=0;s<samples.size();s++){
//       if(samples[s].Type!=2)continue;
// 
//        if(modelMap.find(samples[s].ModelName())==modelMap.end()){
//          modelMap[samples[s].ModelName()].push_back(samples[s]);
//          if(modelMap[samples[s].ModelName()].size()==1)modelVector.push_back(samples[s].ModelName());
//        }
//      }      
//   } 

   //based on the modelMap
   //DrawRatioBands(TkPattern); 
   //DrawRatioBands(MuPattern);
   //DrawRatioBands(MOPattern); 
   //DrawRatioBands(HQPattern);
   //DrawRatioBands(LQPattern);

   //draw the cross section limit for all model
   DrawModelLimitWithBand(TkPattern);
   DrawModelLimitWithBand(MuPattern);
  // DrawModelLimitWithBand(MOPattern);
   //DrawModelLimitWithBand(HQPattern);
   //DrawModelLimitWithBand(LQPattern);
   //make plots of the observed limit for all signal model (and mass point) and save the result in a latex table

   TCanvas* c1;
   TLegend* LEG;
   double LInt = 0;

   FILE* pFile    = fopen((outpath+string("Summary") + ".txt").c_str(),"w");
   FILE* talkFile = fopen((outpath + "TalkPlots" + ".txt").c_str(),"w");

   fprintf(pFile   , "\\documentclass{article}\n");
   fprintf(pFile   , "\\begin{document}\n\n");
   fprintf(talkFile, "\\documentclass{article}\n");
   fprintf(talkFile, "\\usepackage{rotating}\n");
   fprintf(talkFile, "\\begin{document}\n\n");
   fprintf(talkFile, "\\begin{tiny}\n\n");

   fprintf(pFile   , "%% %50s\n", "TkOnly");
   fprintf(pFile   , "\\begin{table}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(pFile   , "%-20s & %-4s & %-6s & %-5s & %-5s & %-4s & %-15s & %-2s & %-4s & %-6s & %-6s & %-6s & %-3s \\\\\n", "Model", "Mass", "pT>", "I>", "TOF>", "Mass>", "Bckg. Pred", "Obs", "Seff", "ThXsec", "ExpXsec", "ObsXsec", "Signif");
   fprintf(talkFile, "\\begin{sidewaystable}\n   \\centering\n      \\begin{tabular}{|l|cccccccc|}\n      \\hline\n");
   fprintf(talkFile,"Sample & Mass(GeV) & Pt(GeV) & $I_{as}$ & TOF & Mass Cut (GeV) & N pred & N observed & Eff & Signif \\\\\n");

   fprintf(talkFile, "\\hline\n");
   TGraph** TkGraphs    = new TGraph*[modelVector.size()];
   TGraph** TkExpGraphs = new TGraph*[modelVector.size()];
   for(unsigned int k=0; k<modelVector.size(); k++){
     TkGraphs[k]    = MakePlot(pFile,talkFile,TkPattern,modelVector[k], 2, modelMap[modelVector[k]], LInt);
     TkExpGraphs[k] = MakePlot(pFile,talkFile,TkPattern,modelVector[k], 1, modelMap[modelVector[k]], LInt);
   }
   fprintf(pFile   ,"      \\end{tabular}\n\\end{table}\n\n");
   fprintf(talkFile,"      \\end{tabular}\n\\end{sidewaystable}\n\n");

   fprintf(pFile   , "%% %50s\n", "TkMuon");
   fprintf(pFile   , "\\begin{table}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile, "\\begin{sidewaystable}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile,"Sample & Mass(GeV) & Pt(GeV) & $I_{as}$ & $#beta^{-1]$ & Mass Cut (GeV) & N pred & N observed & Eff \\\\\n");
   fprintf(talkFile, "\\hline\n");
   TGraph** MuGraphs    = new TGraph*[modelVector.size()];
   TGraph** MuExpGraphs = new TGraph*[modelVector.size()];
   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
      if(isNeutral) continue;//skip charged suppressed models
      MuGraphs[k]    = MakePlot(pFile,talkFile,MuPattern,modelVector[k], 2, modelMap[modelVector[k]], LInt);
      MuExpGraphs[k] = MakePlot(pFile,talkFile,MuPattern,modelVector[k], 1, modelMap[modelVector[k]], LInt);
   }
   fprintf(pFile   ,"      \\end{tabular}\n\\end{table}\n\n");
   fprintf(talkFile,"      \\end{tabular}\n\\end{sidewaystable}\n\n");

   //DISABLED FOR NOW ON the MuonOnly, mCHamps and Frac Charge plots
   /*
   fprintf(pFile   , "%% %50s\n", "MuOnly");
   fprintf(pFile   , "\\begin{table}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile, "\\begin{sidewaystable}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile,"Sample & Mass(GeV) & Pt(GeV) & $I_{as}$ & $#beta^{-1]$ & Mass Cut (GeV) & N pred & N observed & Eff \\\\\n");
   fprintf(talkFile, "\\hline\n");
   TGraph** MOGraphs = new TGraph*[modelVector.size()];
   for(unsigned int k=0; k<modelVector.size(); k++){
     bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
     if(isNeutral) continue;//skip charged suppressed models                                                                                                                      
     MOGraphs[k] = MakePlot(pFile,talkFile,MOPattern,modelVector[k], 2, modelMap[modelVector[k]], LInt);
   }
   fprintf(pFile   ,"      \\end{tabular}\n\\end{table}\n\n");
   fprintf(talkFile,"      \\end{tabular}\n\\end{sidewaystable}\n\n");


   fprintf(pFile   , "%% %50s\n", "multiple charge");
   fprintf(pFile   , "\\begin{table}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile, "\\begin{sidewaystable}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile,"Sample & Mass(GeV) & Pt(GeV) & $I_{as}$ & $#beta^{-1]$ & Mass Cut (GeV) & N pred & N observed & Eff \\\\\n");
   fprintf(talkFile, "\\hline\n");
   TGraph** HQGraphs  = new TGraph*[modelVector.size()];
   for(unsigned int k=0; k<modelVector.size(); k++){
     bool ismHSCP = false;
     if(modelVector[k].find("DY_Q")!=string::npos && modelVector[k].find("o3")==string::npos ) ismHSCP = true;
     if(!ismHSCP) continue;  //  only mHSCP models
     HQGraphs[k] = MakePlot(pFile,talkFile,HQPattern,modelVector[k], 2, modelMap[modelVector[k]], LInt);
   }
   fprintf(pFile   ,"      \\end{tabular}\n\\end{table}\n\n");
   fprintf(talkFile,"      \\end{tabular}\n\\end{sidewaystable}\n\n");

   fprintf(pFile   , "%% %50s\n", "fractionnally charge");
   fprintf(pFile   , "\\begin{table}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile, "\\begin{sidewaystable}\n   \\centering\n      \\begin{tabular}{|l|cccccc|}\n      \\hline\n");
   fprintf(talkFile,"Sample & Mass(GeV) & Pt(GeV) & $I_{as}$ & $#beta^{-1]$ & Mass Cut (GeV) & N pred & N observed & Eff \\\\\n");
   fprintf(talkFile, "\\hline\n");
   TGraph** LQGraphs = new TGraph*[modelVector.size()];
   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isFractional = false;if(modelVector[k].find("1o3")!=string::npos || modelVector[k].find("2o3")!=string::npos ||modelVector[k].find("Q1")!=string::npos)isFractional = true;
      if(!isFractional) continue;//skip q>=1 charged suppressed models
      LQGraphs[k] = MakePlot(pFile,talkFile,LQPattern,modelVector[k], 2, modelMap[modelVector[k]], LInt);
   }
   fprintf(pFile   ,"      \\end{tabular}\n\\end{table}\n\n");
   fprintf(talkFile,"      \\end{tabular}\n\\end{sidewaystable}\n\n");
   */
   fprintf(pFile   ,"\\end{document}\n\n");
   fprintf(talkFile,"\\end{document}\n\n");

   //DISABLED FOR NOW ON, WILL NEED TO BE ADAPTED FOR THE PAPER SUMMARY TABLE
   /*
   if(SQRTS==8.0){
      fprintf(pFile,"%%TKONLY\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummary(pFile, talkFile, TkPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%TKTOF\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummary(pFile, talkFile, MuPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%MUONLY\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummary(pFile, talkFile, MOPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%Q>1\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummary(pFile, talkFile, HQPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%Q<1\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummary(pFile, talkFile, LQPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"\n\n\n\n\n\n");
      fprintf(pFile,"%%TKONLY\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummaryPaper(pFile, talkFile, TkPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%TKTOF\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummaryPaper(pFile, talkFile, MuPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%MUONLY\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummaryPaper(pFile, talkFile, MOPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%Q>1\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummaryPaper(pFile, talkFile, HQPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");

      fprintf(pFile,"%%Q<1\n");
      fprintf(pFile,"Sample & Mass  & Cut   & \\multicolumn{4}{c|}{$\\sqrt{s}=7TeV$} & \\multicolumn{4}{c|}{$\\sqrt{s}=8TeV$} & \\multicolumn{2}{c|}{$\\sqrt{s}=7+8TeV$} \\\\\\hline\n");
      fprintf(pFile,"       & (GeV) & (GeV) & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & Eff & $\\sigma_{TH}$ & $\\sigma_{obs}$ & $\\sigma_{pred}$ & $\\mu_{obs}$ & $\\mu_{pred}$ \\\\\\hline\n");
      for(unsigned int k=0; k<modelVector.size(); k++){printSummaryPaper(pFile, talkFile, LQPattern , modelVector[k], modelMap[modelVector[k]]); }
      fprintf(pFile,"\\hline\n\n\n");
   }*/

   //print a table with all uncertainty on signal efficiency

   c1 = new TCanvas("c1", "c1",600,600);
   c1->SetLeftMargin(0.15);

   TMultiGraph* TkSystGraphs = new TMultiGraph();

   LEG = new TLegend(0.20,0.65,0.80,0.80);
   LEG->SetNColumns(2) ;
   LEG->SetFillColor(0);
   LEG->SetFillStyle(0);
   LEG->SetBorderSize(0);

   fprintf(pFile   ,"\n\n %20s \n\n", LegendFromType(TkPattern).c_str());
   fprintf(pFile   ,          "%20s    Eff   --> PScale |  DeDxScale | PUScale | TotalUncertainty     \n","Model");
   fprintf(talkFile, "\\hline\n%20s &  Eff     & PScale &  DeDxScale & PUScale & TotalUncertainty \\\\\n","Model");
   int Graphs=0;

   for(unsigned int k=0; k<modelVector.size(); k++){
     TGraph* Uncertainty = CheckSignalUncertainty(pFile,talkFile,TkPattern, modelVector[k], modelMap[modelVector[k]]);
     if(Uncertainty!=NULL && useSample(0, modelVector[k])) {
       Uncertainty->SetLineColor(Color[Graphs]);  Uncertainty->SetMarkerColor(Color[Graphs]);   Uncertainty->SetMarkerStyle(20); Uncertainty->SetLineWidth(2);
       TkSystGraphs->Add(Uncertainty,"LP");
       LEG->AddEntry(Uncertainty,  modelVector[k].c_str() ,"L");
       Graphs++;
     }
   }

   if(Graphs>0) {
     TkSystGraphs->Draw("A");
     TkSystGraphs->SetTitle("");
     TkSystGraphs->GetXaxis()->SetTitle("Mass (GeV)");
     TkSystGraphs->GetYaxis()->SetTitle("Relative Uncertainty");
     TkSystGraphs->GetYaxis()->SetTitleOffset(1.40);
     TkSystGraphs->GetYaxis()->SetRangeUser(0.0, 0.70);
     TkSystGraphs->GetYaxis()->SetNdivisions(505, "X");
     
     LEG->Draw();
     c1->SetLogy(false);
     c1->SetGridy(false);
     
     DrawPreliminary(LegendFromType(TkPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
     SaveCanvas(c1,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", "TkUncertainty");
     delete c1;
     delete TkSystGraphs;
     delete LEG;
   }

   c1 = new TCanvas("c1", "c1",600,600);
   c1->SetLeftMargin(0.15);
   TMultiGraph* MuSystGraphs = new TMultiGraph();

   LEG = new TLegend(0.20,0.65,0.80,0.80);
   LEG->SetNColumns(2) ;
   LEG->SetFillColor(0);
   LEG->SetFillStyle(0);
   LEG->SetBorderSize(0);

   fprintf(pFile   ,"\n\n %20s \n\n", LegendFromType(MuPattern).c_str());
   fprintf(pFile,             "%20s   Eff   --> PScale |  DeDxScale | MassScale | HIP | PUScale | TOFScale | TotalUncertainty     \n","Model");
   fprintf(talkFile, "\\hline\n%20s &  Eff    & PScale &  DeDxScale & MassScale & HIP & PUScale & TOFScale & TotalUncertainty \\\\\n","Model");
   Graphs=0;
   for(unsigned int k=0; k<modelVector.size(); k++){
     TGraph* Uncertainty = CheckSignalUncertainty(pFile,talkFile,MuPattern, modelVector[k], modelMap[modelVector[k]]);

     if(Uncertainty!=NULL && useSample(2, modelVector[k])) {
       Uncertainty->SetLineColor(Color[Graphs]);  Uncertainty->SetMarkerColor(Color[Graphs]);   Uncertainty->SetMarkerStyle(20); Uncertainty->SetLineWidth(2);
       MuSystGraphs->Add(Uncertainty,"LP");
       LEG->AddEntry(Uncertainty,  modelVector[k].c_str() ,"L");
       Graphs++;
     }
   }

   if(Graphs>0) {
     MuSystGraphs->Draw("A");
     MuSystGraphs->SetTitle("");
     MuSystGraphs->GetXaxis()->SetTitle("Mass (GeV)");
     MuSystGraphs->GetYaxis()->SetTitle("Relative Uncertainty");
     MuSystGraphs->GetYaxis()->SetTitleOffset(1.40);
     MuSystGraphs->GetYaxis()->SetRangeUser(0.0, 0.7);
     MuSystGraphs->GetYaxis()->SetNdivisions(505, "X");
     
     LEG->Draw();
     c1->SetLogy(false);
     c1->SetGridy(false);
     
     DrawPreliminary(LegendFromType(MuPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
     SaveCanvas(c1,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", "MuUncertainty");
     delete c1;
     delete MuSystGraphs;
     delete LEG;
   }

   /*
   fprintf(pFile   ,"\n\n %20s \n\n", LegendFromType(MOPattern).c_str());
   fprintf(pFile,             "%20s   Eff   --> PScale |  DeDxScale | PUScale | TOFScale | TotalUncertainty     \n","Model");
   fprintf(talkFile, "\\hline\n%20s &  Eff    & PScale &  DeDxScale & PUScale & TOFScale & TotalUncertainty \\\\\n","Model");

   c1 = new TCanvas("c1", "c1",600,600);
   c1->SetLeftMargin(0.15);
   TMultiGraph* MOSystGraphs = new TMultiGraph();

   LEG = new TLegend(0.20,0.75,0.80,0.90);
   LEG->SetNColumns(2) ;
   LEG->SetFillColor(0);
   LEG->SetFillStyle(0);
   LEG->SetBorderSize(0);

   Graphs=0;
   for(unsigned int k=0; k<modelVector.size(); k++){
     TGraph* Uncertainty = CheckSignalUncertainty(pFile,talkFile,MOPattern, modelVector[k], modelMap[modelVector[k]]);
     if(Uncertainty!=NULL && useSample(3, modelVector[k])) {
       Uncertainty->SetLineColor(Color[Graphs]);  Uncertainty->SetMarkerColor(Color[Graphs]);   Uncertainty->SetMarkerStyle(20); Uncertainty->SetLineWidth(2);
       MOSystGraphs->Add(Uncertainty,"LP");
       LEG->AddEntry(Uncertainty,  modelVector[k].c_str() ,"L");
       Graphs++;
     }
   }
   if(Graphs>0) {
   MOSystGraphs->Draw("A");
   MOSystGraphs->SetTitle("");
   MOSystGraphs->GetXaxis()->SetTitle("Mass (GeV)");
   MOSystGraphs->GetYaxis()->SetTitle("Relative Uncertainty");
   MOSystGraphs->GetYaxis()->SetTitleOffset(1.40);
   MOSystGraphs->GetYaxis()->SetRangeUser(0.0, 0.6);
   MOSystGraphs->GetYaxis()->SetNdivisions(520, "X");

   LEG->Draw();
   c1->SetLogy(false);
   c1->SetGridy(false);

   DrawPreliminary(LegendFromType(MOPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", "MOUncertainty");
   delete c1;
   delete MOSystGraphs;
   delete LEG;
   }

   c1 = new TCanvas("c1", "c1",600,600);
   c1->SetLeftMargin(0.15);
   TMultiGraph* LQSystGraphs = new TMultiGraph();

   LEG = new TLegend(0.20,0.75,0.80,0.90);
   LEG->SetNColumns(2) ;
   LEG->SetFillColor(0);
   LEG->SetFillStyle(0);
   LEG->SetBorderSize(0);

   fprintf(pFile   ,"\n\n %20s \n\n", LegendFromType(LQPattern).c_str());
   fprintf(pFile   ,          "%20s    Eff   --> PScale |  DeDxScale | PUScale | TotalUncertainty     \n","Model");
   fprintf(talkFile, "\\hline\n%20s &  Eff     & PScale &  DeDxScale & PUScale & TotalUncertainty \\\\\n","Model");

   Graphs=0;
   for(unsigned int k=0; k<modelVector.size(); k++){
     TGraph* Uncertainty = CheckSignalUncertainty(pFile,talkFile,LQPattern, modelVector[k], modelMap[modelVector[k]]);
     if(Uncertainty!=NULL && useSample(5, modelVector[k])) {
       Uncertainty->SetLineColor(Color[Graphs]);  Uncertainty->SetMarkerColor(Color[Graphs]);   Uncertainty->SetMarkerStyle(20); Uncertainty->SetLineWidth(2);
       LQSystGraphs->Add(Uncertainty,"LP");
       LEG->AddEntry(Uncertainty,  modelVector[k].c_str() ,"L");
       Graphs++;
     }
   }

   if(Graphs>0) {
   LQSystGraphs->Draw("A");
   LQSystGraphs->SetTitle("");
   LQSystGraphs->GetXaxis()->SetTitle("Mass (GeV)");
   LQSystGraphs->GetYaxis()->SetTitle("Relative Uncertainty");
   LQSystGraphs->GetYaxis()->SetTitleOffset(1.40);
   LQSystGraphs->GetYaxis()->SetRangeUser(0.0, 0.6);
   LQSystGraphs->GetYaxis()->SetNdivisions(520, "X");

   LEG->Draw();
   c1->SetLogy(false);
   c1->SetGridy(false);

   DrawPreliminary(LegendFromType(LQPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", "LQUncertainty");
   delete c1;
   delete LQSystGraphs;
   delete LEG;
   }


   c1 = new TCanvas("c1", "c1",600,600);
   c1->SetLeftMargin(0.15);
   TMultiGraph* HQSystGraphs = new TMultiGraph();

   LEG = new TLegend(0.20,0.75,0.80,0.90);
   LEG->SetNColumns(2) ;
   LEG->SetFillColor(0);
   LEG->SetFillStyle(0);
   LEG->SetBorderSize(0);

   fprintf(pFile   ,"\n\n %20s \n\n", LegendFromType(HQPattern).c_str());
   fprintf(pFile   ,          "%20s    Eff   --> PScale |  DeDxScale | PUScale | TotalUncertainty     \n","Model");
   fprintf(talkFile, "\\hline\n%20s &  Eff     & PScale &  DeDxScale & PUScale & TotalUncertainty \\\\\n","Model");

   Graphs=0;
   for(unsigned int k=0; k<modelVector.size(); k++){
     TGraph* Uncertainty = CheckSignalUncertainty(pFile,talkFile,HQPattern, modelVector[k], modelMap[modelVector[k]]);
     if(Uncertainty!=NULL && useSample(4, modelVector[k])) {
       Uncertainty->SetLineColor(Color[Graphs]);  Uncertainty->SetMarkerColor(Color[Graphs]);   Uncertainty->SetMarkerStyle(20); Uncertainty->SetLineWidth(2);
       HQSystGraphs->Add(Uncertainty,"LP");
       LEG->AddEntry(Uncertainty,  modelVector[k].c_str() ,"L");
       Graphs++;
     }
   }

   if(Graphs>0) {
   HQSystGraphs->Draw("A");
   HQSystGraphs->SetTitle("");
   HQSystGraphs->GetXaxis()->SetTitle("Mass (GeV)");
   HQSystGraphs->GetYaxis()->SetTitle("Relative Uncertainty");
   HQSystGraphs->GetYaxis()->SetTitleOffset(1.40);
   HQSystGraphs->GetYaxis()->SetRangeUser(0.0, 0.6);
   HQSystGraphs->GetYaxis()->SetNdivisions(520, "X");

   LEG->Draw();
   c1->SetLogy(false);
   c1->SetGridy(false);

   DrawPreliminary(LegendFromType(HQPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c1,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", "HQUncertainty");
   delete c1;
   delete HQSystGraphs;
   delete LEG;
   }
   */
std::cout<<"TESTA\n";

   //Get Theoretical xsection and error bands
   TGraph** ThXSec    = new TGraph*[modelVector.size()];
   TCutG ** ThXSecErr = new TCutG* [modelVector.size()];
   for(unsigned int k=0; k<modelVector.size(); k++){
     if(modelVector[k].find("Gluino")!=string::npos){
         if(SQRTS==7){
            ThXSec   [k] = new TGraph(sizeof(THXSEC7TeV_Gluino_Mass)/sizeof(double),THXSEC7TeV_Gluino_Mass,THXSEC7TeV_Gluino_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr",sizeof(THXSEC7TeV_Gluino_Mass)/sizeof(double),THXSEC7TeV_Gluino_Mass,THXSEC7TeV_Gluino_Low,THXSEC7TeV_Gluino_High, PlotMinScale, PlotMaxScale);
         }else if(SQRTS==8){
            ThXSec   [k] = new TGraph(sizeof(THXSEC8TeV_Gluino_Mass)/sizeof(double),THXSEC8TeV_Gluino_Mass,THXSEC8TeV_Gluino_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr",sizeof(THXSEC8TeV_Gluino_Mass)/sizeof(double),THXSEC8TeV_Gluino_Mass,THXSEC8TeV_Gluino_Low,THXSEC8TeV_Gluino_High, PlotMinScale, PlotMaxScale);
         }else if(SQRTS==13 || SQRTS==1315 || SQRTS==1316 || SQRTS==13167 || SQRTS==131677){ 
            ThXSec   [k] = new TGraph(sizeof(THXSEC13TeV_Gluino_Mass)/sizeof(double),THXSEC13TeV_Gluino_Mass,THXSEC13TeV_Gluino_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr",sizeof(THXSEC13TeV_Gluino_Mass)/sizeof(double),THXSEC13TeV_Gluino_Mass,THXSEC13TeV_Gluino_Low,THXSEC13TeV_Gluino_High, PlotMinScale, PlotMaxScale);
         }else{
 	    const int NMass=sizeof(THXSEC13TeV_Gluino_Mass)/sizeof(double);
	    double ones[NMass]; for(int i=0; i<NMass; i++) ones[i]=1;
	    ThXSec   [k] = new TGraph(NMass,THXSEC13TeV_Gluino_Mass,ones);
	 }
      }else if(modelVector[k].find("Stop"  )!=string::npos){
         if(SQRTS==7){
            ThXSec   [k] = new TGraph(sizeof(THXSEC7TeV_Stop_Mass)/sizeof(double),THXSEC7TeV_Stop_Mass,THXSEC7TeV_Stop_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr",sizeof(THXSEC7TeV_Stop_Mass)/sizeof(double),THXSEC7TeV_Stop_Mass,THXSEC7TeV_Stop_Low,THXSEC7TeV_Stop_High, PlotMinScale, PlotMaxScale);
         }else if(SQRTS==8){
            ThXSec   [k] = new TGraph(sizeof(THXSEC8TeV_Stop_Mass)/sizeof(double),THXSEC8TeV_Stop_Mass,THXSEC8TeV_Stop_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr",sizeof(THXSEC8TeV_Stop_Mass)/sizeof(double),THXSEC8TeV_Stop_Mass,THXSEC8TeV_Stop_Low,THXSEC8TeV_Stop_High, PlotMinScale, PlotMaxScale);
         }else if(SQRTS==13 || SQRTS==1315 || SQRTS==1316 || SQRTS==13167 || SQRTS==131677){
            ThXSec   [k] = new TGraph(sizeof(THXSEC13TeV_Stop_Mass)/sizeof(double),THXSEC13TeV_Stop_Mass,THXSEC13TeV_Stop_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr",sizeof(THXSEC13TeV_Stop_Mass)/sizeof(double),THXSEC13TeV_Stop_Mass,THXSEC13TeV_Stop_Low,THXSEC13TeV_Stop_High, PlotMinScale, PlotMaxScale);
         }else{
            const int NMass=sizeof(THXSEC13TeV_Stop_Mass)/sizeof(double);
            double ones[NMass]; for(int i=0; i<NMass; i++) ones[i]=1;
            ThXSec   [k] = new TGraph(NMass,THXSEC13TeV_Stop_Mass,ones);
         }
      }else if(modelVector[k].find("GMStau"  )!=string::npos){
         if(SQRTS==7){
            ThXSec   [k] = MakePlot(NULL, NULL, TkPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", sizeof(THXSEC7TeV_GMStau_Mass)/sizeof(double),THXSEC7TeV_GMStau_Mass,THXSEC7TeV_GMStau_Low,THXSEC7TeV_GMStau_High, PlotMinScale, PlotMaxScale);
         }else if(SQRTS==8){
            ThXSec   [k] = MakePlot(NULL, NULL, TkPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", sizeof(THXSEC8TeV_GMStau_Mass)/sizeof(double),THXSEC8TeV_GMStau_Mass,THXSEC8TeV_GMStau_Low,THXSEC8TeV_GMStau_High, PlotMinScale, PlotMaxScale);
         }else if(SQRTS==13 || SQRTS==1315 || SQRTS==1316 || SQRTS==13167 || SQRTS==131677){
//            #Prospino xsection that I get looks very weird, use pythia for the time being
            ThXSec   [k] = new TGraph(sizeof(THXSEC13TeV_GMStau_Mass)/sizeof(double),THXSEC13TeV_GMStau_Mass,THXSEC13TeV_GMStau_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", sizeof(THXSEC13TeV_GMStau_Mass)/sizeof(double),THXSEC13TeV_GMStau_Mass,THXSEC13TeV_GMStau_Low,THXSEC13TeV_GMStau_High, PlotMinScale, PlotMaxScale);
//            ThXSec   [k] = MakePlot(NULL, NULL, TkPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);
//            double* XSecErrLow  = new double[ThXSec[k]->GetN()];
//            double* XSecErrHigh = new double[ThXSec[k]->GetN()];
//            for(int i=0;i<ThXSec[k]->GetN();i++){ XSecErrLow[i] = ThXSec[k]->GetY()[i]*0.90; XSecErrHigh[i] = ThXSec[k]->GetY()[i]*1.10; }
//            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", ThXSec[k]->GetN(),ThXSec[k]->GetX(),XSecErrLow,XSecErrHigh, PlotMinScale, PlotMaxScale);            
         }else{
            const int NMass=sizeof(THXSEC13TeV_GMStau_Mass)/sizeof(double);
            double ones[NMass]; for(int i=0; i<NMass; i++) ones[i]=1;
            ThXSec   [k] = new TGraph(NMass,THXSEC13TeV_GMStau_Mass,ones);
         }
      }else if(modelVector[k].find("PPStau"  )!=string::npos){
         if(SQRTS==7){
            ThXSec   [k] = MakePlot(NULL, NULL, TkPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);   
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", sizeof(THXSEC7TeV_PPStau_Mass)/sizeof(double),THXSEC7TeV_PPStau_Mass,THXSEC7TeV_PPStau_Low,THXSEC7TeV_PPStau_High, PlotMinScale, PlotMaxScale); 
         }else if(SQRTS==8){
            ThXSec   [k] = MakePlot(NULL, NULL, TkPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", sizeof(THXSEC8TeV_PPStau_Mass)/sizeof(double),THXSEC8TeV_PPStau_Mass,THXSEC8TeV_PPStau_Low,THXSEC8TeV_PPStau_High, PlotMinScale, PlotMaxScale);
         }else if(SQRTS==13 || SQRTS==1315 || SQRTS==1316 || SQRTS==13167 || SQRTS==131677){
//            #Prospino xsection that I get looks very weird, use pythia for the time being
            ThXSec   [k] = new TGraph(sizeof(THXSEC13TeV_PPStau_Mass)/sizeof(double),THXSEC13TeV_PPStau_Mass,THXSEC13TeV_PPStau_Cen);
            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", sizeof(THXSEC13TeV_PPStau_Mass)/sizeof(double),THXSEC13TeV_PPStau_Mass,THXSEC13TeV_PPStau_Low,THXSEC13TeV_PPStau_High, PlotMinScale, PlotMaxScale);
//            ThXSec   [k] = MakePlot(NULL, NULL, TkPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);            
//            double* XSecErrLow  = new double[ThXSec[k]->GetN()];
//            double* XSecErrHigh = new double[ThXSec[k]->GetN()];
//            for(int i=0;i<ThXSec[k]->GetN();i++){ XSecErrLow[i] = ThXSec[k]->GetY()[i]*0.90; XSecErrHigh[i] = ThXSec[k]->GetY()[i]*1.10; }
//            ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", ThXSec[k]->GetN(),ThXSec[k]->GetX(),XSecErrLow,XSecErrHigh, PlotMinScale, PlotMaxScale);
         }else{
           const int NMass=sizeof(THXSEC13TeV_PPStau_Mass)/sizeof(double);
           double ones[NMass]; for(int i=0; i<NMass; i++) ones[i]=1;
           ThXSec   [k] = new TGraph(NMass,THXSEC13TeV_PPStau_Mass,ones);
         }
       }else{
         //if(modelVector[k].find("o3")!=string::npos){
         //   ThXSec   [k] = MakePlot(NULL, NULL, LQPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);
         //}else if(modelVector[k].find("DY_Q")!=string::npos){
         //   ThXSec   [k] = MakePlot(NULL, NULL, HQPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);
         //}else{
            ThXSec   [k] = MakePlot(NULL, NULL, TkPattern,modelVector[k], 0, modelMap[modelVector[k]], LInt);
         //}
         double* XSecErrLow  = new double[ThXSec[k]->GetN()];
         double* XSecErrHigh = new double[ThXSec[k]->GetN()];
         //assume 10% error on xsection
         for(int i=0;i<ThXSec[k]->GetN();i++){ XSecErrLow[i] = ThXSec[k]->GetY()[i]*0.90; XSecErrHigh[i] = ThXSec[k]->GetY()[i]*1.10; }
         ThXSecErr[k] = GetErrorBand(modelVector[k]+"ThErr", ThXSec[k]->GetN(),ThXSec[k]->GetX(),XSecErrLow,XSecErrHigh, 1E-5, PlotMaxScale); 
       }
   }

std::cout<<"TESTB\n";

   //Print the excluded mass range
   fprintf(pFile,"\n\n\n=======================\n Observed mass range excluded   \n=========================\n");
   fprintf(pFile,"-----------------------\n0%% TK-ONLY       \n-------------------------\n");
   for(unsigned int k=0; k<modelVector.size(); k++){
      if(TkGraphs[k]->GetN()==0) continue;
      if(TkGraphs[k]->GetX()[TkGraphs[k]->GetN()-1]<0) continue;
      fprintf(pFile,"%20s --> Excluded mass below %8.3fGeV\n", modelVector[k].c_str(), FindIntersectionBetweenTwoGraphs(TkGraphs[k],  ThXSec[k], TkGraphs[k]->GetX()[0], TkGraphs[k]->GetX()[TkGraphs[k]->GetN()-1], 1, 0.00));
   }

   fprintf(pFile,"-----------------------\n0%% MU+TOF        \n-------------------------\n");
   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
      if(isNeutral) continue;//skip charged suppressed models
      if(MuGraphs[k]->GetN()==0) continue;
      if(MuGraphs[k]->GetX()[MuGraphs[k]->GetN()-1]<0) continue;
      double minMass=-1, maxMass=-1;
      FindRangeBetweenTwoGraphs(MuGraphs[k],  ThXSec[k], MuGraphs[k]->GetX()[0], MuGraphs[k]->GetX()[MuGraphs[k]->GetN()-1], 1, 0.00, minMass, maxMass);
      fprintf(pFile,"%20s --> Excluded mass range %8.3f - %8.3fGeV\n", modelVector[k].c_str(), minMass, maxMass);
      //fprintf(pFile,"%20s --> Excluded mass below %8.3fGeV\n", modelVector[k].c_str(), FindIntersectionBetweenTwoGraphs(MuGraphs[k],  ThXSec[k], MuGraphs[k]->GetX()[0], MuGraphs[k]->GetX()[MuGraphs[k]->GetN()-1], 1, 0.00));
   }  

   // JOZZE :: EDIT EXPECTED MASS RANGE EXCLUDED
   fprintf(pFile,"\n\n\n=======================\n Expected mass range excluded   \n=========================\n");
   fprintf(pFile,"-----------------------\n0%% TK-ONLY       \n-------------------------\n");
   for(unsigned int k=0; k<modelVector.size(); k++){
      if(TkExpGraphs[k]->GetN()==0) continue;
      if(TkExpGraphs[k]->GetX()[TkExpGraphs[k]->GetN()-1]<0) continue;
      fprintf(pFile,"%20s --> Excluded mass below %8.3fGeV\n", modelVector[k].c_str(), FindIntersectionBetweenTwoGraphs(TkExpGraphs[k],  ThXSec[k], TkExpGraphs[k]->GetX()[0], TkExpGraphs[k]->GetX()[TkExpGraphs[k]->GetN()-1], 1, 0.00));
   }

   fprintf(pFile,"-----------------------\n0%% MU+TOF        \n-------------------------\n");
   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
      if(isNeutral) continue;//skip charged suppressed models
      if(MuExpGraphs[k]->GetN()==0) continue;
      if(MuExpGraphs[k]->GetX()[MuExpGraphs[k]->GetN()-1]<0) continue;
      double minMass=-1, maxMass=-1;
      FindRangeBetweenTwoGraphs(MuExpGraphs[k],  ThXSec[k], MuExpGraphs[k]->GetX()[0], MuExpGraphs[k]->GetX()[MuExpGraphs[k]->GetN()-1], 1, 0.00, minMass, maxMass);
      fprintf(pFile,"%20s --> Excluded mass range %8.3f - %8.3fGeV\n", modelVector[k].c_str(), minMass, maxMass);
      //fprintf(pFile,"%20s --> Excluded mass below %8.3fGeV\n", modelVector[k].c_str(), FindIntersectionBetweenTwoGraphs(MuGraphs[k],  ThXSec[k], MuGraphs[k]->GetX()[0], MuGraphs[k]->GetX()[MuGraphs[k]->GetN()-1], 1, 0.00));
   }
   /*
   fprintf(pFile,"-----------------------\n0%% MU+Only        \n-------------------------\n");
   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
      if(isNeutral) continue;//skip charged suppressed models
      if(MOGraphs[k]->GetN()==0) continue;
      if(MOGraphs[k]->GetX()[MOGraphs[k]->GetN()-1]<0) continue;
      fprintf(pFile,"%20s --> Excluded mass below %8.3fGeV\n", modelVector[k].c_str(), FindIntersectionBetweenTwoGraphs(MOGraphs[k],  ThXSec[k], MOGraphs[k]->GetX()[0], MOGraphs[k]->GetX()[MOGraphs[k]->GetN()-1], 1, 0.00));
   }

   fprintf(pFile,"-----------------------\n0%% Q>1            \n-------------------------\n");
   for(unsigned int k=0; k<modelVector.size(); k++){
     bool ismHSCP = false;
     if(modelVector[k].find("DY_Q")!=string::npos && modelVector[k].find("o3")==string::npos ) ismHSCP = true;
     if(!ismHSCP) continue; 

      if(HQGraphs[k]->GetN()==0) continue;
      if(HQGraphs[k]->GetX()[HQGraphs[k]->GetN()-1]<0) continue;
      double minMass=-1, maxMass=-1;
      FindRangeBetweenTwoGraphs(HQGraphs[k],  ThXSec[k], HQGraphs[k]->GetX()[0], HQGraphs[k]->GetX()[HQGraphs[k]->GetN()-1], 1, 0.00, minMass, maxMass);
      fprintf(pFile,"%20s --> Excluded mass range %8.3f - %8.3fGeV\n", modelVector[k].c_str(), minMass, maxMass);
      //fprintf(pFile,"%20s --> Excluded mass below %8.3fGeV\n", modelVector[k].c_str(), FindIntersectionBetweenTwoGraphs(HQGraphs[k],  ThXSec[k], HQGraphs[k]->GetX()[0], HQGraphs[k]->GetX()[HQGraphs[k]->GetN()-1], 1, 0.00));
   }

   fprintf(pFile,"-----------------------\n0%% Q<1             \n-------------------------\n");
   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isFractional = false;if(modelVector[k].find("1o3")!=string::npos || modelVector[k].find("2o3")!=string::npos || modelVector[k].find("DY_Q1")!=string::npos)isFractional = true;
      if(!isFractional) continue;//skip non fractional charge models
      if(LQGraphs[k]->GetN()==0) continue;
      if(LQGraphs[k]->GetX()[LQGraphs[k]->GetN()-1]<0) continue;
      fprintf(pFile,"%20s --> Excluded mass below %8.3fGeV\n", modelVector[k].c_str(), FindIntersectionBetweenTwoGraphs(LQGraphs[k],  ThXSec[k], LQGraphs[k]->GetX()[0], LQGraphs[k]->GetX()[LQGraphs[k]->GetN()-1], 1, 0.00));
   }
   */
   fclose(pFile);
  

std::cout<<"TESTC\n";

 std::cout<<"qui arriva\n";
   //Make the final plot with all curves in it
   // I don't like much this part because it is dependent of what is in Analysis_Samples.h in an hardcoded way   
   std::map<string, TGraph*> TkGraphMap;
   std::map<string, TGraph*> MuGraphMap;
   std::map<string, TGraph*> LQGraphMap;
   std::map<string, TGraph*> HQGraphMap;
   std::map<string, TGraph*> MOGraphMap;
   std::map<string, TGraph*> ThGraphMap;
   std::map<string, TCutG* > ThErrorMap;
   for(unsigned int k=0; k<modelVector.size(); k++){
     std::cout << modelVector[k] << std::endl;
      TkGraphMap[modelVector[k]] = TkGraphs [k];
      MuGraphMap[modelVector[k]] = MuGraphs [k];
      //LQGraphMap[modelVector[k]] = LQGraphs [k];
      //HQGraphMap[modelVector[k]] = HQGraphs [k];
      //MOGraphMap[modelVector[k]] = MOGraphs [k];
      ThGraphMap[modelVector[k]] = ThXSec   [k];
      ThErrorMap[modelVector[k]] = ThXSecErr[k];
   }
   std::cout<<"qui anche \n";

   // FIXME add GluinoN_f50 maybe as well?
   if (MODE.find("13TeV16")==std::string::npos){
     std::cout<<"entrare anche 1"<<std::endl;
     ThGraphMap["Gluino_f10"   ]->SetLineColor(4);  ThGraphMap["Gluino_f10"   ]->SetMarkerColor(4);   ThGraphMap["Gluino_f10"   ]->SetLineWidth(1);   ThGraphMap["Gluino_f10"   ]->SetLineStyle(1);  ThGraphMap["Gluino_f10"   ]->SetMarkerStyle(1);
     MuGraphMap["Gluino_f10"   ]->SetLineColor(4);  MuGraphMap["Gluino_f10"   ]->SetMarkerColor(4);   MuGraphMap["Gluino_f10"   ]->SetLineWidth(2);   MuGraphMap["Gluino_f10"   ]->SetLineStyle(1);  MuGraphMap["Gluino_f10"   ]->SetMarkerStyle(22);
     MuGraphMap["Gluino_f50"   ]->SetLineColor(4);  MuGraphMap["Gluino_f50"   ]->SetMarkerColor(4);   MuGraphMap["Gluino_f50"   ]->SetLineWidth(2);   MuGraphMap["Gluino_f50"   ]->SetLineStyle(1);  MuGraphMap["Gluino_f50"   ]->SetMarkerStyle(23);
     TkGraphMap["Gluino_f10"   ]->SetLineColor(4);  TkGraphMap["Gluino_f10"   ]->SetMarkerColor(4);   TkGraphMap["Gluino_f10"   ]->SetLineWidth(2);   TkGraphMap["Gluino_f10"   ]->SetLineStyle(1);  TkGraphMap["Gluino_f10"   ]->SetMarkerStyle(22);
     TkGraphMap["Gluino_f50"   ]->SetLineColor(4);  TkGraphMap["Gluino_f50"   ]->SetMarkerColor(4);   TkGraphMap["Gluino_f50"   ]->SetLineWidth(2);   TkGraphMap["Gluino_f50"   ]->SetLineStyle(1);  TkGraphMap["Gluino_f50"   ]->SetMarkerStyle(23);
     TkGraphMap["GluinoN_f10"  ]->SetLineColor(4);  TkGraphMap["GluinoN_f10"  ]->SetMarkerColor(4);   TkGraphMap["GluinoN_f10"  ]->SetLineWidth(2);   TkGraphMap["GluinoN_f10"  ]->SetLineStyle(1);  TkGraphMap["GluinoN_f10"  ]->SetMarkerStyle(26);
      ThGraphMap["Stop"         ]->SetLineColor(2);  ThGraphMap["Stop"         ]->SetMarkerColor(2);   ThGraphMap["Stop"         ]->SetLineWidth(1);   ThGraphMap["Stop"         ]->SetLineStyle(2);  ThGraphMap["Stop"         ]->SetMarkerStyle(1);
      MuGraphMap["Stop"         ]->SetLineColor(2);  MuGraphMap["Stop"         ]->SetMarkerColor(2);   MuGraphMap["Stop"         ]->SetLineWidth(2);   MuGraphMap["Stop"         ]->SetLineStyle(1);  MuGraphMap["Stop"         ]->SetMarkerStyle(21);
      TkGraphMap["Stop"         ]->SetLineColor(2);  TkGraphMap["Stop"         ]->SetMarkerColor(2);   TkGraphMap["Stop"         ]->SetLineWidth(2);   TkGraphMap["Stop"         ]->SetLineStyle(1);  TkGraphMap["Stop"         ]->SetMarkerStyle(21);
      TkGraphMap["StopN"        ]->SetLineColor(2);  TkGraphMap["StopN"        ]->SetMarkerColor(2);   TkGraphMap["StopN"        ]->SetLineWidth(2);   TkGraphMap["StopN"        ]->SetLineStyle(1);  TkGraphMap["StopN"        ]->SetMarkerStyle(25);
      ThGraphMap["GMStau"       ]->SetLineColor(1);  ThGraphMap["GMStau"       ]->SetMarkerColor(1);   ThGraphMap["GMStau"       ]->SetLineWidth(1);   ThGraphMap["GMStau"       ]->SetLineStyle(3);  ThGraphMap["GMStau"       ]->SetMarkerStyle(1);
      ThGraphMap["PPStau"       ]->SetLineColor(6);  ThGraphMap["PPStau"       ]->SetMarkerColor(6);   ThGraphMap["PPStau"       ]->SetLineWidth(1);   ThGraphMap["PPStau"       ]->SetLineStyle(4);  ThGraphMap["PPStau"       ]->SetMarkerStyle(1);
      MuGraphMap["GMStau"       ]->SetLineColor(1);  MuGraphMap["GMStau"       ]->SetMarkerColor(1);   MuGraphMap["GMStau"       ]->SetLineWidth(2);   MuGraphMap["GMStau"       ]->SetLineStyle(1);  MuGraphMap["GMStau"       ]->SetMarkerStyle(20);
      MuGraphMap["PPStau"       ]->SetLineColor(6);  MuGraphMap["PPStau"       ]->SetMarkerColor(6);   MuGraphMap["PPStau"       ]->SetLineWidth(2);   MuGraphMap["PPStau"       ]->SetLineStyle(1);  MuGraphMap["PPStau"       ]->SetMarkerStyle(20);
      TkGraphMap["GMStau"       ]->SetLineColor(1);  TkGraphMap["GMStau"       ]->SetMarkerColor(1);   TkGraphMap["GMStau"       ]->SetLineWidth(2);   TkGraphMap["GMStau"       ]->SetLineStyle(1);  TkGraphMap["GMStau"       ]->SetMarkerStyle(20);
      TkGraphMap["PPStau"       ]->SetLineColor(6);  TkGraphMap["PPStau"       ]->SetMarkerColor(6);   TkGraphMap["PPStau"       ]->SetLineWidth(2);   TkGraphMap["PPStau"       ]->SetLineStyle(1);  TkGraphMap["PPStau"       ]->SetMarkerStyle(20);
      ThGraphMap["DY_Q1"        ]->SetLineColor(46); ThGraphMap["DY_Q1"        ]->SetMarkerColor(46);  ThGraphMap["DY_Q1"        ]->SetLineWidth(1);   ThGraphMap["DY_Q1"        ]->SetLineStyle(8);  ThGraphMap["DY_Q1"        ]->SetMarkerStyle(1);
      MuGraphMap["DY_Q1"        ]->SetLineColor(46); MuGraphMap["DY_Q1"        ]->SetMarkerColor(46);  MuGraphMap["DY_Q1"        ]->SetLineWidth(2);   MuGraphMap["DY_Q1"        ]->SetLineStyle(1);  MuGraphMap["DY_Q1"      ]->SetMarkerStyle(20);
      TkGraphMap["DY_Q1"        ]->SetLineColor(46); TkGraphMap["DY_Q1"        ]->SetMarkerColor(46);  TkGraphMap["DY_Q1"        ]->SetLineWidth(2);   TkGraphMap["DY_Q1"        ]->SetLineStyle(1);  TkGraphMap["DY_Q1"      ]->SetMarkerStyle(20);
      ThGraphMap["DY_Q2"        ]->SetLineColor(43); ThGraphMap["DY_Q2"        ]->SetMarkerColor(43);  ThGraphMap["DY_Q2"        ]->SetLineWidth(1);   ThGraphMap["DY_Q2"        ]->SetLineStyle(10); ThGraphMap["DY_Q2"      ]->SetMarkerStyle(1);
      MuGraphMap["DY_Q2"        ]->SetLineColor(43); MuGraphMap["DY_Q2"        ]->SetMarkerColor(43);  MuGraphMap["DY_Q2"        ]->SetLineWidth(2);   MuGraphMap["DY_Q2"        ]->SetLineStyle(1);  MuGraphMap["DY_Q2"      ]->SetMarkerStyle(34);
      TkGraphMap["DY_Q2"        ]->SetLineColor(43); TkGraphMap["DY_Q2"        ]->SetMarkerColor(43);  TkGraphMap["DY_Q2"        ]->SetLineWidth(2);   TkGraphMap["DY_Q2"        ]->SetLineStyle(1);  TkGraphMap["DY_Q2"      ]->SetMarkerStyle(34);
      std::cout<<"esce anche"<<std::endl;
   
}


   else if (MODE.find("13TeV16")!=std::string::npos && MODE.find("13TeV16G")==std::string::npos){
     ThGraphMap["Gluino16_f10"   ]->SetLineColor(4);  ThGraphMap["Gluino16_f10"   ]->SetMarkerColor(4);   ThGraphMap["Gluino16_f10"   ]->SetLineWidth(1);   ThGraphMap["Gluino16_f10"   ]->SetLineStyle(1);  ThGraphMap["Gluino16_f10"   ]->SetMarkerStyle(1);
     MuGraphMap["Gluino16_f10"   ]->SetLineColor(4);  MuGraphMap["Gluino16_f10"   ]->SetMarkerColor(4);   MuGraphMap["Gluino16_f10"   ]->SetLineWidth(2);   MuGraphMap["Gluino16_f10"   ]->SetLineStyle(1);  MuGraphMap["Gluino16_f10"   ]->SetMarkerStyle(22);
     MuGraphMap["Gluino16_f50"   ]->SetLineColor(4);  MuGraphMap["Gluino16_f50"   ]->SetMarkerColor(4);   MuGraphMap["Gluino16_f50"   ]->SetLineWidth(2);   MuGraphMap["Gluino16_f50"   ]->SetLineStyle(1);  MuGraphMap["Gluino16_f50"   ]->SetMarkerStyle(23);
     TkGraphMap["Gluino16_f10"   ]->SetLineColor(4);  TkGraphMap["Gluino16_f10"   ]->SetMarkerColor(4);   TkGraphMap["Gluino16_f10"   ]->SetLineWidth(2);   TkGraphMap["Gluino16_f10"   ]->SetLineStyle(1);  TkGraphMap["Gluino16_f10"   ]->SetMarkerStyle(22);
     TkGraphMap["Gluino16_f50"   ]->SetLineColor(4);  TkGraphMap["Gluino16_f50"   ]->SetMarkerColor(4);   TkGraphMap["Gluino16_f50"   ]->SetLineWidth(2);   TkGraphMap["Gluino16_f50"   ]->SetLineStyle(1);  TkGraphMap["Gluino16_f50"   ]->SetMarkerStyle(23);
     TkGraphMap["Gluino16N_f10"  ]->SetLineColor(4);  TkGraphMap["Gluino16N_f10"  ]->SetMarkerColor(4);   TkGraphMap["Gluino16N_f10"  ]->SetLineWidth(2);   TkGraphMap["Gluino16N_f10"  ]->SetLineStyle(1);  TkGraphMap["Gluino16N_f10"  ]->SetMarkerStyle(26);
     ThGraphMap["Stop16"         ]->SetLineColor(2);  ThGraphMap["Stop16"         ]->SetMarkerColor(2);   ThGraphMap["Stop16"         ]->SetLineWidth(1);   ThGraphMap["Stop16"         ]->SetLineStyle(2);  ThGraphMap["Stop16"         ]->SetMarkerStyle(1);
     MuGraphMap["Stop16"         ]->SetLineColor(2);  MuGraphMap["Stop16"         ]->SetMarkerColor(2);   MuGraphMap["Stop16"         ]->SetLineWidth(2);   MuGraphMap["Stop16"         ]->SetLineStyle(1);  MuGraphMap["Stop16"         ]->SetMarkerStyle(21);
     TkGraphMap["Stop16"         ]->SetLineColor(2);  TkGraphMap["Stop16"         ]->SetMarkerColor(2);   TkGraphMap["Stop16"         ]->SetLineWidth(2);   TkGraphMap["Stop16"         ]->SetLineStyle(1);  TkGraphMap["Stop16"         ]->SetMarkerStyle(21);
     TkGraphMap["Stop16N"        ]->SetLineColor(2);  TkGraphMap["Stop16N"        ]->SetMarkerColor(2);   TkGraphMap["Stop16N"        ]->SetLineWidth(2);   TkGraphMap["Stop16N"        ]->SetLineStyle(1);  TkGraphMap["Stop16N"        ]->SetMarkerStyle(25);
     ThGraphMap["GMStau16"       ]->SetLineColor(1);  ThGraphMap["GMStau16"       ]->SetMarkerColor(1);   ThGraphMap["GMStau16"       ]->SetLineWidth(1);   ThGraphMap["GMStau16"       ]->SetLineStyle(3);  ThGraphMap["GMStau16"       ]->SetMarkerStyle(1);
     ThGraphMap["PPStau16"       ]->SetLineColor(6);  ThGraphMap["PPStau16"       ]->SetMarkerColor(6);   ThGraphMap["PPStau16"       ]->SetLineWidth(1);   ThGraphMap["PPStau16"       ]->SetLineStyle(4);  ThGraphMap["PPStau16"       ]->SetMarkerStyle(1);
     MuGraphMap["GMStau16"       ]->SetLineColor(1);  MuGraphMap["GMStau16"       ]->SetMarkerColor(1);   MuGraphMap["GMStau16"       ]->SetLineWidth(2);   MuGraphMap["GMStau16"       ]->SetLineStyle(1);  MuGraphMap["GMStau16"       ]->SetMarkerStyle(20);
     MuGraphMap["PPStau16"       ]->SetLineColor(6);  MuGraphMap["PPStau16"       ]->SetMarkerColor(6);   MuGraphMap["PPStau16"       ]->SetLineWidth(2);   MuGraphMap["PPStau16"       ]->SetLineStyle(1);  MuGraphMap["PPStau16"       ]->SetMarkerStyle(20);
     TkGraphMap["GMStau16"       ]->SetLineColor(1);  TkGraphMap["GMStau16"       ]->SetMarkerColor(1);   TkGraphMap["GMStau16"       ]->SetLineWidth(2);   TkGraphMap["GMStau16"       ]->SetLineStyle(1);  TkGraphMap["GMStau16"       ]->SetMarkerStyle(20);
     TkGraphMap["PPStau16"       ]->SetLineColor(6);  TkGraphMap["PPStau16"       ]->SetMarkerColor(6);   TkGraphMap["PPStau16"       ]->SetLineWidth(2);   TkGraphMap["PPStau16"       ]->SetLineStyle(1);  TkGraphMap["PPStau16"       ]->SetMarkerStyle(20);
     ThGraphMap["DY16_Q1"        ]->SetLineColor(46); ThGraphMap["DY16_Q1"        ]->SetMarkerColor(46);  ThGraphMap["DY16_Q1"        ]->SetLineWidth(1);   ThGraphMap["DY16_Q1"        ]->SetLineStyle(8);  ThGraphMap["DY16_Q1"        ]->SetMarkerStyle(1);
     MuGraphMap["DY16_Q1"        ]->SetLineColor(46); MuGraphMap["DY16_Q1"        ]->SetMarkerColor(46);  MuGraphMap["DY16_Q1"        ]->SetLineWidth(2);   MuGraphMap["DY16_Q1"        ]->SetLineStyle(1);  MuGraphMap["DY16_Q1"      ]->SetMarkerStyle(20);
     TkGraphMap["DY16_Q1"        ]->SetLineColor(46); TkGraphMap["DY16_Q1"        ]->SetMarkerColor(46);  TkGraphMap["DY16_Q1"        ]->SetLineWidth(2);   TkGraphMap["DY16_Q1"        ]->SetLineStyle(1);  TkGraphMap["DY16_Q1"      ]->SetMarkerStyle(20);
     ThGraphMap["DY16_Q2"        ]->SetLineColor(43); ThGraphMap["DY16_Q2"        ]->SetMarkerColor(43);  ThGraphMap["DY16_Q2"        ]->SetLineWidth(1);   ThGraphMap["DY16_Q2"        ]->SetLineStyle(10); ThGraphMap["DY16_Q2"      ]->SetMarkerStyle(1);
     MuGraphMap["DY16_Q2"        ]->SetLineColor(43); MuGraphMap["DY16_Q2"        ]->SetMarkerColor(43);  MuGraphMap["DY16_Q2"        ]->SetLineWidth(2);   MuGraphMap["DY16_Q2"        ]->SetLineStyle(1);  MuGraphMap["DY16_Q2"      ]->SetMarkerStyle(34);
     TkGraphMap["DY16_Q2"        ]->SetLineColor(43); TkGraphMap["DY16_Q2"        ]->SetMarkerColor(43);  TkGraphMap["DY16_Q2"        ]->SetLineWidth(2);   TkGraphMap["DY16_Q2"        ]->SetLineStyle(1);  TkGraphMap["DY16_Q2"      ]->SetMarkerStyle(34);
   }


   else if (MODE.find("13TeV16G")!=std::string::npos){
     ThGraphMap["Gluino16G_f10"   ]->SetLineColor(4);  ThGraphMap["Gluino16G_f10"   ]->SetMarkerColor(4);   ThGraphMap["Gluino16G_f10"   ]->SetLineWidth(1);   ThGraphMap["Gluino16G_f10"   ]->SetLineStyle(1);  ThGraphMap["Gluino16G_f10"   ]->SetMarkerStyle(1);
     MuGraphMap["Gluino16G_f10"   ]->SetLineColor(4);  MuGraphMap["Gluino16G_f10"   ]->SetMarkerColor(4);   MuGraphMap["Gluino16G_f10"   ]->SetLineWidth(2);   MuGraphMap["Gluino16G_f10"   ]->SetLineStyle(1);  MuGraphMap["Gluino16G_f10"   ]->SetMarkerStyle(22);
     MuGraphMap["Gluino16G_f50"   ]->SetLineColor(4);  MuGraphMap["Gluino16G_f50"   ]->SetMarkerColor(4);   MuGraphMap["Gluino16G_f50"   ]->SetLineWidth(2);   MuGraphMap["Gluino16G_f50"   ]->SetLineStyle(1);  MuGraphMap["Gluino16G_f50"   ]->SetMarkerStyle(23);
     TkGraphMap["Gluino16G_f10"   ]->SetLineColor(4);  TkGraphMap["Gluino16G_f10"   ]->SetMarkerColor(4);   TkGraphMap["Gluino16G_f10"   ]->SetLineWidth(2);   TkGraphMap["Gluino16G_f10"   ]->SetLineStyle(1);  TkGraphMap["Gluino16G_f10"   ]->SetMarkerStyle(22);
     TkGraphMap["Gluino16G_f50"   ]->SetLineColor(4);  TkGraphMap["Gluino16G_f50"   ]->SetMarkerColor(4);   TkGraphMap["Gluino16G_f50"   ]->SetLineWidth(2);   TkGraphMap["Gluino16G_f50"   ]->SetLineStyle(1);  TkGraphMap["Gluino16G_f50"   ]->SetMarkerStyle(23);
     TkGraphMap["Gluino16GN_f10"  ]->SetLineColor(4);  TkGraphMap["Gluino16GN_f10"  ]->SetMarkerColor(4);   TkGraphMap["Gluino16GN_f10"  ]->SetLineWidth(2);   TkGraphMap["Gluino16GN_f10"  ]->SetLineStyle(1);  TkGraphMap["Gluino16GN_f10"  ]->SetMarkerStyle(26);
     ThGraphMap["Stop16G"         ]->SetLineColor(2);  ThGraphMap["Stop16G"         ]->SetMarkerColor(2);   ThGraphMap["Stop16G"         ]->SetLineWidth(1);   ThGraphMap["Stop16G"         ]->SetLineStyle(2);  ThGraphMap["Stop16G"         ]->SetMarkerStyle(1);
     MuGraphMap["Stop16G"         ]->SetLineColor(2);  MuGraphMap["Stop16G"         ]->SetMarkerColor(2);   MuGraphMap["Stop16G"         ]->SetLineWidth(2);   MuGraphMap["Stop16G"         ]->SetLineStyle(1);  MuGraphMap["Stop16G"         ]->SetMarkerStyle(21);
     TkGraphMap["Stop16G"         ]->SetLineColor(2);  TkGraphMap["Stop16G"         ]->SetMarkerColor(2);   TkGraphMap["Stop16G"         ]->SetLineWidth(2);   TkGraphMap["Stop16G"         ]->SetLineStyle(1);  TkGraphMap["Stop16G"         ]->SetMarkerStyle(21);
     TkGraphMap["Stop16GN"        ]->SetLineColor(2);  TkGraphMap["Stop16GN"        ]->SetMarkerColor(2);   TkGraphMap["Stop16GN"        ]->SetLineWidth(2);   TkGraphMap["Stop16GN"        ]->SetLineStyle(1);  TkGraphMap["Stop16GN"        ]->SetMarkerStyle(25);
     ThGraphMap["GMStau16G"       ]->SetLineColor(1);  ThGraphMap["GMStau16G"       ]->SetMarkerColor(1);   ThGraphMap["GMStau16G"       ]->SetLineWidth(1);   ThGraphMap["GMStau16G"       ]->SetLineStyle(3);  ThGraphMap["GMStau16G"       ]->SetMarkerStyle(1);
     ThGraphMap["PPStau16G"       ]->SetLineColor(6);  ThGraphMap["PPStau16G"       ]->SetMarkerColor(6);   ThGraphMap["PPStau16G"       ]->SetLineWidth(1);   ThGraphMap["PPStau16G"       ]->SetLineStyle(4);  ThGraphMap["PPStau16G"       ]->SetMarkerStyle(1);
     MuGraphMap["GMStau16G"       ]->SetLineColor(1);  MuGraphMap["GMStau16G"       ]->SetMarkerColor(1);   MuGraphMap["GMStau16G"       ]->SetLineWidth(2);   MuGraphMap["GMStau16G"       ]->SetLineStyle(1);  MuGraphMap["GMStau16G"       ]->SetMarkerStyle(20);
     MuGraphMap["PPStau16G"       ]->SetLineColor(6);  MuGraphMap["PPStau16G"       ]->SetMarkerColor(6);   MuGraphMap["PPStau16G"       ]->SetLineWidth(2);   MuGraphMap["PPStau16G"       ]->SetLineStyle(1);  MuGraphMap["PPStau16G"       ]->SetMarkerStyle(20);
     TkGraphMap["GMStau16G"       ]->SetLineColor(1);  TkGraphMap["GMStau16G"       ]->SetMarkerColor(1);   TkGraphMap["GMStau16G"       ]->SetLineWidth(2);   TkGraphMap["GMStau16G"       ]->SetLineStyle(1);  TkGraphMap["GMStau16G"       ]->SetMarkerStyle(20);
     TkGraphMap["PPStau16G"       ]->SetLineColor(6);  TkGraphMap["PPStau16G"       ]->SetMarkerColor(6);   TkGraphMap["PPStau16G"       ]->SetLineWidth(2);   TkGraphMap["PPStau16G"       ]->SetLineStyle(1);  TkGraphMap["PPStau16G"       ]->SetMarkerStyle(20);
     ThGraphMap["DY16G_Q1"        ]->SetLineColor(46); ThGraphMap["DY16G_Q1"        ]->SetMarkerColor(46);  ThGraphMap["DY16G_Q1"        ]->SetLineWidth(1);   ThGraphMap["DY16G_Q1"        ]->SetLineStyle(8);  ThGraphMap["DY16G_Q1"        ]->SetMarkerStyle(1);
     MuGraphMap["DY16G_Q1"        ]->SetLineColor(46); MuGraphMap["DY16G_Q1"        ]->SetMarkerColor(46);  MuGraphMap["DY16G_Q1"        ]->SetLineWidth(2);   MuGraphMap["DY16G_Q1"        ]->SetLineStyle(1);  MuGraphMap["DY16G_Q1"      ]->SetMarkerStyle(20);
     TkGraphMap["DY16G_Q1"        ]->SetLineColor(46); TkGraphMap["DY16G_Q1"        ]->SetMarkerColor(46);  TkGraphMap["DY16G_Q1"        ]->SetLineWidth(2);   TkGraphMap["DY16G_Q1"        ]->SetLineStyle(1);  TkGraphMap["DY16G_Q1"      ]->SetMarkerStyle(20);
     ThGraphMap["DY16G_Q2"        ]->SetLineColor(43); ThGraphMap["DY16G_Q2"        ]->SetMarkerColor(43);  ThGraphMap["DY16G_Q2"        ]->SetLineWidth(1);   ThGraphMap["DY16G_Q2"        ]->SetLineStyle(10); ThGraphMap["DY16G_Q2"      ]->SetMarkerStyle(1);
     MuGraphMap["DY16G_Q2"        ]->SetLineColor(43); MuGraphMap["DY16G_Q2"        ]->SetMarkerColor(43);  MuGraphMap["DY16G_Q2"        ]->SetLineWidth(2);   MuGraphMap["DY16G_Q2"        ]->SetLineStyle(1);  MuGraphMap["DY16G_Q2"      ]->SetMarkerStyle(34);
     TkGraphMap["DY16G_Q2"        ]->SetLineColor(43); TkGraphMap["DY16G_Q2"        ]->SetMarkerColor(43);  TkGraphMap["DY16G_Q2"        ]->SetLineWidth(2);   TkGraphMap["DY16G_Q2"        ]->SetLineStyle(1);  TkGraphMap["DY16G_Q2"      ]->SetMarkerStyle(34);
   }


   //MOGraphMap["Gluino_f10"   ]->SetLineColor(4);  MOGraphMap["Gluino_f10"   ]->SetMarkerColor(4);   MOGraphMap["Gluino_f10"   ]->SetLineWidth(2);   MOGraphMap["Gluino_f10"   ]->SetLineStyle(1);  MOGraphMap["Gluino_f10"   ]->SetMarkerStyle(22);
   //MOGraphMap["Gluino_f50"   ]->SetLineColor(4);  MOGraphMap["Gluino_f50"   ]->SetMarkerColor(4);   MOGraphMap["Gluino_f50"   ]->SetLineWidth(2);   MOGraphMap["Gluino_f50"   ]->SetLineStyle(1);  MOGraphMap["Gluino_f50"   ]->SetMarkerStyle(23);
   //MOGraphMap["Gluino_f100"  ]->SetLineColor(4);  MOGraphMap["Gluino_f100"  ]->SetMarkerColor(4);   MOGraphMap["Gluino_f100"  ]->SetLineWidth(2);   MOGraphMap["Gluino_f100"  ]->SetLineStyle(1);  MOGraphMap["Gluino_f100"  ]->SetMarkerStyle(26);
   //MOGraphMap["Stop"         ]->SetLineColor(2);  MOGraphMap["Stop"         ]->SetMarkerColor(2);   MOGraphMap["Stop"         ]->SetLineWidth(2);   MOGraphMap["Stop"         ]->SetLineStyle(1);  MOGraphMap["Stop"         ]->SetMarkerStyle(21);
   //ThGraphMap["DY_Q1o3"      ]->SetLineColor(41); ThGraphMap["DY_Q1o3"      ]->SetMarkerColor(41);  ThGraphMap["DY_Q1o3"      ]->SetLineWidth(1);   ThGraphMap["DY_Q1o3"      ]->SetLineStyle(9);  ThGraphMap["DY_Q1o3"      ]->SetMarkerStyle(1);
   //TkGraphMap["DY_Q1o3"      ]->SetLineColor(41); TkGraphMap["DY_Q1o3"      ]->SetMarkerColor(41);  TkGraphMap["DY_Q1o3"      ]->SetLineWidth(2);   TkGraphMap["DY_Q1o3"      ]->SetLineStyle(1);  TkGraphMap["DY_Q1o3"      ]->SetMarkerStyle(33);
   //MOGraphMap["DY_Q1o3"      ]->SetLineColor(41); MOGraphMap["DY_Q1o3"      ]->SetMarkerColor(41);  MOGraphMap["DY_Q1o3"      ]->SetLineWidth(2);   MOGraphMap["DY_Q1o3"      ]->SetLineStyle(1);  MOGraphMap["DY_Q1o3"      ]->SetMarkerStyle(33);
   //LQGraphMap["DY_Q1o3"      ]->SetLineColor(41); LQGraphMap["DY_Q1o3"      ]->SetMarkerColor(41);  LQGraphMap["DY_Q1o3"      ]->SetLineWidth(2);   LQGraphMap["DY_Q1o3"      ]->SetLineStyle(1);  LQGraphMap["DY_Q1o3"      ]->SetMarkerStyle(33);
   //ThGraphMap["DY_Q2o3"      ]->SetLineColor(43); ThGraphMap["DY_Q2o3"      ]->SetMarkerColor(43);  ThGraphMap["DY_Q2o3"      ]->SetLineWidth(1);   ThGraphMap["DY_Q2o3"      ]->SetLineStyle(10); ThGraphMap["DY_Q2o3"      ]->SetMarkerStyle(1);
   //TkGraphMap["DY_Q2o3"      ]->SetLineColor(43); TkGraphMap["DY_Q2o3"      ]->SetMarkerColor(43);  TkGraphMap["DY_Q2o3"      ]->SetLineWidth(2);   TkGraphMap["DY_Q2o3"      ]->SetLineStyle(1);  TkGraphMap["DY_Q2o3"      ]->SetMarkerStyle(34);
   //MuGraphMap["DY_Q2o3"      ]->SetLineColor(43); MuGraphMap["DY_Q2o3"      ]->SetMarkerColor(43);  MuGraphMap["DY_Q2o3"      ]->SetLineWidth(2);   MuGraphMap["DY_Q2o3"      ]->SetLineStyle(1);  MuGraphMap["DY_Q2o3"      ]->SetMarkerStyle(34);
   //MOGraphMap["DY_Q2o3"      ]->SetLineColor(43); MOGraphMap["DY_Q2o3"      ]->SetMarkerColor(43);  MOGraphMap["DY_Q2o3"      ]->SetLineWidth(2);   MOGraphMap["DY_Q2o3"      ]->SetLineStyle(1);  MOGraphMap["DY_Q2o3"      ]->SetMarkerStyle(34);
   //LQGraphMap["DY_Q2o3"      ]->SetLineColor(43); LQGraphMap["DY_Q2o3"      ]->SetMarkerColor(43);  LQGraphMap["DY_Q2o3"      ]->SetLineWidth(2);   LQGraphMap["DY_Q2o3"      ]->SetLineStyle(1);  LQGraphMap["DY_Q2o3"      ]->SetMarkerStyle(34);
   //LQGraphMap["DY_Q1"        ]->SetLineColor(46); LQGraphMap["DY_Q1"        ]->SetMarkerColor(46);  LQGraphMap["DY_Q1"        ]->SetLineWidth(2);   LQGraphMap["DY_Q1"        ]->SetLineStyle(1);  LQGraphMap["DY_Q1"        ]->SetMarkerStyle(20);
   //HQGraphMap["DY_Q1"        ]->SetLineColor(46); HQGraphMap["DY_Q1"        ]->SetMarkerColor(46);  HQGraphMap["DY_Q1"        ]->SetLineWidth(2);   HQGraphMap["DY_Q1"        ]->SetLineStyle(1);  HQGraphMap["DY_Q1"        ]->SetMarkerStyle(20);
   //HQGraphMap["DY_Q2"        ]->SetLineColor(2 ); HQGraphMap["DY_Q2"        ]->SetMarkerColor(2 );  HQGraphMap["DY_Q2"        ]->SetLineWidth(2);   HQGraphMap["DY_Q2"        ]->SetLineStyle(1);  HQGraphMap["DY_Q2"        ]->SetMarkerStyle(21);
   //ThGraphMap["DY_Q3"        ]->SetLineColor(1 ); ThGraphMap["DY_Q3"        ]->SetMarkerColor(1 );  ThGraphMap["DY_Q3"        ]->SetLineWidth(1);   ThGraphMap["DY_Q3"        ]->SetLineStyle(9);  ThGraphMap["DY_Q3"        ]->SetMarkerStyle(1);
   //HQGraphMap["DY_Q3"        ]->SetLineColor(1 ); HQGraphMap["DY_Q3"        ]->SetMarkerColor(1 );  HQGraphMap["DY_Q3"        ]->SetLineWidth(2);   HQGraphMap["DY_Q3"        ]->SetLineStyle(1);  HQGraphMap["DY_Q3"        ]->SetMarkerStyle(22);
   //ThGraphMap["DY_Q4"        ]->SetLineColor(6 ); ThGraphMap["DY_Q4"        ]->SetMarkerColor(6 );  ThGraphMap["DY_Q4"        ]->SetLineWidth(1);   ThGraphMap["DY_Q4"        ]->SetLineStyle(3);  ThGraphMap["DY_Q4"        ]->SetMarkerStyle(1);
   //HQGraphMap["DY_Q4"        ]->SetLineColor(6 ); HQGraphMap["DY_Q4"        ]->SetMarkerColor(6 );  HQGraphMap["DY_Q4"        ]->SetLineWidth(2);   HQGraphMap["DY_Q4"        ]->SetLineStyle(1);  HQGraphMap["DY_Q4"        ]->SetMarkerStyle(23);
   //ThGraphMap["DY_Q5"        ]->SetLineColor(4 ); ThGraphMap["DY_Q5"        ]->SetMarkerColor(4 );  ThGraphMap["DY_Q5"        ]->SetLineWidth(1);   ThGraphMap["DY_Q5"        ]->SetLineStyle(8);  ThGraphMap["DY_Q5"        ]->SetMarkerStyle(1);
   //HQGraphMap["DY_Q5"        ]->SetLineColor(4 ); HQGraphMap["DY_Q5"        ]->SetMarkerColor(4 );  HQGraphMap["DY_Q5"        ]->SetLineWidth(2);   HQGraphMap["DY_Q5"        ]->SetLineStyle(1);  HQGraphMap["DY_Q5"        ]->SetMarkerStyle(29);
   //ThGraphMap["DY_Q6"        ]->SetLineColor(9 ); ThGraphMap["DY_Q6"        ]->SetMarkerColor(9 );  ThGraphMap["DY_Q6"        ]->SetLineWidth(1);   ThGraphMap["DY_Q6"        ]->SetLineStyle(6);  ThGraphMap["DY_Q6"        ]->SetMarkerStyle(1);
   //HQGraphMap["DY_Q6"        ]->SetLineColor(9 ); HQGraphMap["DY_Q6"        ]->SetMarkerColor(9 );  HQGraphMap["DY_Q6"        ]->SetLineWidth(2);   HQGraphMap["DY_Q6"        ]->SetLineStyle(1);  HQGraphMap["DY_Q6"        ]->SetMarkerStyle(33);
   //ThGraphMap["DY_Q7"        ]->SetLineColor(12); ThGraphMap["DY_Q7"        ]->SetMarkerColor(12);  ThGraphMap["DY_Q7"        ]->SetLineWidth(1);   ThGraphMap["DY_Q7"        ]->SetLineStyle(7);  ThGraphMap["DY_Q7"        ]->SetMarkerStyle(1);
   //HQGraphMap["DY_Q7"        ]->SetLineColor(12); HQGraphMap["DY_Q7"        ]->SetMarkerColor(12);  HQGraphMap["DY_Q7"        ]->SetLineWidth(2);   HQGraphMap["DY_Q7"        ]->SetLineStyle(1);  HQGraphMap["DY_Q7"        ]->SetMarkerStyle(34);
   //ThGraphMap["DY_Q8"        ]->SetLineColor(14); ThGraphMap["DY_Q8"        ]->SetMarkerColor(14);  ThGraphMap["DY_Q8"        ]->SetLineWidth(1);   ThGraphMap["DY_Q8"        ]->SetLineStyle(10); ThGraphMap["DY_Q8"        ]->SetMarkerStyle(1);
   //HQGraphMap["DY_Q8"        ]->SetLineColor(14); HQGraphMap["DY_Q8"        ]->SetMarkerColor(14);  HQGraphMap["DY_Q8"        ]->SetLineWidth(2);   HQGraphMap["DY_Q8"        ]->SetLineStyle(1);  HQGraphMap["DY_Q8"        ]->SetMarkerStyle(24);

   std::cout<<"qui anche 2\n";

std::cout<<"TESTD\n";
   c1 = new TCanvas("c1", "c1",600,600);
   c1->SetLogy(true);

   TH1D* frame = new TH1D("frame", "frame", 1,50, 2650);
   frame->GetXaxis()->SetNdivisions(505);
   frame->SetTitle("");
   frame->SetStats(kFALSE);
   frame->GetXaxis()->SetTitle("Mass (GeV)");
   frame->GetYaxis()->SetTitle(/*Combine?"95% CL limit on #sigma/#sigma_{th}":*/"95% CL limit on #sigma (pb)");
   frame->GetYaxis()->SetTitleOffset(1.40);
   frame->SetMaximum(PlotMaxScale);
   frame->SetMinimum(PlotMinScale);
   frame->GetYaxis()->SetRangeUser(!Combine?1e-4:2e-5, !Combine?1.5e1:2.5e1); // JOZZE EDIT
   frame->Draw("AXIS");


//   MGMu->Draw("A");

   std::string tmp = "";
   if      (MODE.find("13TeV16G")!=std::string::npos) tmp = "16G";
   else if (MODE.find("13TeV16")!=std::string::npos)  tmp = "16" ;

////   if(!Combine) {
      ThErrorMap["Gluino"+tmp+"_f10"]->Draw("F");
      ThGraphMap["Gluino"+tmp+"_f10" ]->Draw("L");

      ThErrorMap["Stop"+tmp      ]->Draw("F");
      ThGraphMap["Stop"+tmp      ]->Draw("L");

      ThErrorMap["GMStau"+tmp    ]->Draw("F");
      ThGraphMap["GMStau"+tmp     ]->Draw("L");

      ThErrorMap["PPStau"+tmp    ]->Draw("F");
      ThGraphMap["PPStau"+tmp     ]->Draw("L");

      //ThErrorMap["DY_Q2o3"   ]->Draw("F");
      //ThGraphMap["DY_Q2o3"    ]->Draw("L");

      ThErrorMap["DY"+tmp+"_Q1"   ]->Draw("F");
      ThGraphMap["DY"+tmp+"_Q1"    ]->Draw("L");
                         
      ThErrorMap["DY"+tmp+"_Q2"   ]->Draw("F");
      ThGraphMap["DY"+tmp+"_Q2"    ]->Draw("L");
////   }else{
////      TLine* LineAtOne = new TLine(50,1,1550,1);      LineAtOne->SetLineStyle(3);   LineAtOne->Draw(); // FIXME do we really need this for 2016? Since it's at 13TeV, we might not really need it finally
////   }

   MuGraphMap["Gluino"+tmp+"_f10"]->Draw("LP");
   MuGraphMap["Gluino"+tmp+"_f50"]->Draw("LP");
   MuGraphMap["Stop"+tmp      ]->Draw("LP");
   MuGraphMap["GMStau"+tmp    ]->Draw("LP");
   MuGraphMap["PPStau"+tmp    ]->Draw("LP");
   MuGraphMap["DY"+tmp+"_Q1"     ]->Draw("LP");
   MuGraphMap["DY"+tmp+"_Q2"     ]->Draw("LP");

   DrawPreliminary(LegendFromType(MuPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
   TLegend* LEGMu = /*!Combine ? */new TLegend(0.60,0.82-7*0.043,0.93,0.82)/* : new TLegend(0.60,0.15,0.93,0.15+7*0.043)*/;
   LEGMu->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGMu->SetTextSize(18); //font size
   LEGMu->SetFillColor(0); 
   LEGMu->SetFillStyle(0);
   LEGMu->SetBorderSize(0);
   LEGMu->AddEntry(MuGraphMap["Gluino"+tmp+"_f10" ] , "gluino; 50% #tilde{g}g"    ,"LP");
   LEGMu->AddEntry(MuGraphMap["Gluino"+tmp+"_f50"] , "gluino; 10% #tilde{g}g"    ,"LP");
   LEGMu->AddEntry(MuGraphMap["Stop"+tmp      ] , "stop"                      ,"LP");
   LEGMu->AddEntry(MuGraphMap["GMStau"+tmp    ] , "stau; dir. prod."           ,"LP");
   LEGMu->AddEntry(MuGraphMap["PPStau"+tmp    ] , "stau"                 ,"LP");
   LEGMu->AddEntry(MuGraphMap["DY"+tmp+"_Q1"     ], "|Q| = 1e"                ,"LP");
   LEGMu->AddEntry(MuGraphMap["DY"+tmp+"_Q2"     ], "|Q| = 2e"                ,"LP");

   TLegend* LEGTh = new TLegend(0.25,0.82-(1+6)*0.043,0.60,0.82);
   LEGTh->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGTh->SetTextSize(18); //font size
////   if(!Combine) {
   LEGTh->SetHeader("Theoretical Prediction");
   LEGTh->SetFillColor(0);
   //LEGTh->SetFillStyle(0);
   LEGTh->SetBorderSize(0);

   TGraph* GlThLeg = (TGraph*) ThGraphMap["Gluino"+tmp+"_f10"]->Clone("GluinoThLeg");
   GlThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGTh->AddEntry(GlThLeg, "gluino (NLO+NLL)" ,"LF");
   TGraph* StThLeg = (TGraph*) ThGraphMap["Stop"+tmp      ]->Clone("StopThLeg");
   StThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGTh->AddEntry(StThLeg   ,"stop (NLO+NLL)" ,"LF");
   TGraph* PPStauThLeg = (TGraph*) ThGraphMap["PPStau"+tmp        ]->Clone("PPStauThLeg");
   PPStauThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGTh->AddEntry(PPStauThLeg   ,"stau, dir. prod. (NLO)" ,"LF");
   TGraph* StauThLeg = (TGraph*) ThGraphMap["GMStau"+tmp        ]->Clone("StauThLeg");
   StauThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGTh->AddEntry(StauThLeg   ,"stau (NLO)" ,"LF");
   //TGraph* DYQ2o3ThLeg = (TGraph*) ThGraphMap["DY_Q2o3"        ]->Clone("DYQ2o3ThLeg");
   //DYQ2o3ThLeg->SetFillColor(ThErrorMap["DY_Q2o3"]->GetFillColor());
   //LEGTh->AddEntry(DYQ2o3ThLeg   ,"|Q| = 2e/3 (LO)" ,"LF");
   TGraph* DYQ1ThLeg = (TGraph*) ThGraphMap["DY"+tmp+"_Q1"        ]->Clone("DYQ1ThLeg");
   DYQ1ThLeg->SetFillColor(ThErrorMap["DY"+tmp+"_Q1"]->GetFillColor());
   LEGTh->AddEntry(DYQ1ThLeg   ,"|Q| = 1e (LO)" ,"LF");
   TGraph* DYQ2ThLeg = (TGraph*) ThGraphMap["DY"+tmp+"_Q2"        ]->Clone("DYQ2ThLeg");
   DYQ2ThLeg->SetFillColor(ThErrorMap["DY"+tmp+"_Q2"]->GetFillColor());
   LEGTh->AddEntry(DYQ2ThLeg   ,"|Q| = 2e (LO)" ,"LF");
   LEGTh->Draw();
////   }
   LEGMu->Draw();

   c1->SetLogy(true);
   SaveCanvas(c1, outpath, string("MuExclusionLog"));
   delete c1;

   c1 = new TCanvas("c1", "c1",600,600);
   c1->SetLogy(true);
//   frame = new TH1D("frame", "frame", 1,50, 2650);
//   frame->GetXaxis()->SetNdivisions(505);
//   frame->SetTitle("");
//   frame->SetStats(kFALSE);
//   frame->GetXaxis()->SetTitle("Mass (GeV)");
//   frame->GetYaxis()->SetTitle(Combine?"95% CL limit on #sigma/#sigma_{th}":"95% CL limit on #sigma (pb)");
//   frame->GetYaxis()->SetTitleOffset(1.40);
//   frame->SetMaximum(PlotMaxScale);
//   frame->SetMinimum(PlotMinScale);
   frame->GetYaxis()->SetRangeUser(!Combine?1e-4:2e-5, !Combine?1.5e1:2.5e1); // JOZZE EDIT
   frame->Draw("AXIS");

////   if(!Combine) {
      ThErrorMap["Gluino"+tmp+"_f10"]->Draw("F");
      ThGraphMap["Gluino"+tmp+"_f10"]->Draw("L");

      ThErrorMap["Stop"+tmp      ]->Draw("F");
      ThGraphMap["Stop"+tmp      ]->Draw("L");

      ThErrorMap["GMStau"+tmp    ]->Draw("F");
      ThGraphMap["GMStau"+tmp    ]->Draw("L");

      ThErrorMap["PPStau"+tmp    ]->Draw("F");
      ThGraphMap["PPStau"+tmp    ]->Draw("L");

      ThErrorMap["DY"+tmp+"_Q1"   ]->Draw("F");
      ThGraphMap["DY"+tmp+"_Q1"   ]->Draw("L");
                          
      ThErrorMap["DY"+tmp+"_Q2"   ]->Draw("F");
      ThGraphMap["DY"+tmp+"_Q2"   ]->Draw("L");

      //ThErrorMap["DY_Q2o3"   ]->Draw("F");
      //ThGraphMap["DY_Q2o3"   ]->Draw("L");
////   }else{
////      TLine* LineAtOne = new TLine(50,1,1550,1);      LineAtOne->SetLineStyle(3);   LineAtOne->Draw();
////   }

   TkGraphMap["Gluino"+tmp+"_f10" ]->Draw("LP");
   TkGraphMap["Gluino"+tmp+"_f50" ]->Draw("LP");
   TkGraphMap["Gluino"+tmp+"N_f10"]->Draw("LP");
   TkGraphMap["Stop"+tmp      ]->Draw("LP");
   TkGraphMap["Stop"+tmp+"N"     ]->Draw("LP");
   TkGraphMap["GMStau"+tmp    ]->Draw("LP");
   TkGraphMap["PPStau"+tmp    ]->Draw("LP");
   TkGraphMap["DY"+tmp+"_Q1"     ]->Draw("LP");
   TkGraphMap["DY"+tmp+"_Q2"     ]->Draw("LP");
   //TkGraphMap["DY_Q2o3"    ]->Draw("LP");

   DrawPreliminary(LegendFromType(TkPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));

   TLegend* LEGTk = /*!Combine ? */new TLegend(0.60,0.82-8*0.043,0.93,0.82)/* : new TLegend(0.60,0.15,0.93,0.15+8*0.043)*/;
   LEGTk->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGTk->SetTextSize(18); //font size
   LEGTk->SetFillColor(0); 
   LEGTk->SetFillStyle(0);
   LEGTk->SetBorderSize(0);
   LEGTk->AddEntry(TkGraphMap["Gluino"+tmp+"_f50" ], "gluino; 50% #tilde{g}g"            ,"LP");
   LEGTk->AddEntry(TkGraphMap["Gluino"+tmp+"_f10" ], "gluino; 10% #tilde{g}g"            ,"LP");
   LEGTk->AddEntry(TkGraphMap["Gluino"+tmp+"N_f10"], "gluino; 10% #tilde{g}g; CS"        ,"LP");
   LEGTk->AddEntry(TkGraphMap["Stop"+tmp      ], "stop"                              ,"LP");
   LEGTk->AddEntry(TkGraphMap["Stop"+tmp+"N"     ], "stop; CS"                          ,"LP");
   LEGTk->AddEntry(TkGraphMap["GMStau"+tmp    ], "stau; dir. prod."                ,"LP");
   LEGTk->AddEntry(TkGraphMap["PPStau"+tmp    ], "stau"                              ,"LP");
   LEGTk->AddEntry(TkGraphMap["DY"+tmp+"_Q1"     ], "|Q| = 1e"                            ,"LP");
   LEGTk->AddEntry(TkGraphMap["DY"+tmp+"_Q2"     ], "|Q| = 2e"                            ,"LP");
   //LEGTk->AddEntry(TkGraphMap["DY_Q2o3"    ], "|Q| = 2e/3"                            ,"LP");
   /* 2016 G JOZZE */

   TLegend* LEGThTk = new TLegend(0.25,0.82-(1+6)*0.043,0.60,0.82);
   LEGThTk->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGThTk->SetTextSize(18); //font size
////   if(!Combine) {
   LEGThTk->SetHeader("Theoretical Prediction");
   LEGThTk->SetFillColor(0);
   //LEGThTk->SetFillStyle(0);
   LEGThTk->SetBorderSize(0);
////   TGraph* GlThLeg = (TGraph*) ThGraphMap["Gluino"+tmp+"_f10"]->Clone("GluinoThLeg");
   GlThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGThTk->AddEntry(GlThLeg, "gluino (NLO+NLL)" ,"LF");
////   TGraph* StThLeg = (TGraph*) ThGraphMap["Stop"+tmp      ]->Clone("StopThLeg");
   StThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGThTk->AddEntry(StThLeg   ,"stop (NLO+NLL)" ,"LF");
////   TGraph* PPStauThLeg = (TGraph*) ThGraphMap["PPStau"+tmp        ]->Clone("PPStauThLeg");
   PPStauThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGThTk->AddEntry(PPStauThLeg   ,"stau; dir. prod. (NLO)" ,"LF");
////   TGraph* StauThLeg = (TGraph*) ThGraphMap["GMStau"+tmp        ]->Clone("StauThLeg");
   StauThLeg->SetFillColor(ThErrorMap["Gluino"+tmp+"_f10"]->GetFillColor());
   LEGThTk->AddEntry(StauThLeg   ,"stau (NLO)" ,"LF");
   //TGraph* DYQ2o3ThLeg = (TGraph*) ThGraphMap["DY_Q2o3"        ]->Clone("DYQ2o3ThLeg");
   //DYQ2o3ThLeg->SetFillColor(ThErrorMap["DY_Q2o3"]->GetFillColor());
   //LEGThTk->AddEntry(DYQ2o3ThLeg   ,"|Q| = 2e/3 (LO)" ,"LF");
////   TGraph* DYQ1ThLeg = (TGraph*) ThGraphMap["DY"+tmp+"_Q1"        ]->Clone("DYQ1ThLeg");
   DYQ1ThLeg->SetFillColor(ThErrorMap["DY"+tmp+"_Q1"]->GetFillColor());
   LEGThTk->AddEntry(DYQ1ThLeg   ,"|Q| = 1e (LO)" ,"LF");
////   TGraph* DYQ2ThLeg = (TGraph*) ThGraphMap["DY"+tmp+"_Q2"        ]->Clone("DYQ2ThLeg");
   DYQ2ThLeg->SetFillColor(ThErrorMap["DY"+tmp+"_Q2"]->GetFillColor());
   LEGThTk->AddEntry(DYQ2ThLeg   ,"|Q| = 2e (LO)" ,"LF");

   LEGThTk->Draw();
////   }

/*   if(!Combine)*/ LEGThTk->Draw();
   LEGTk->Draw();
   c1->SetLogy(true);
   SaveCanvas(c1, outpath, string("TkExclusionLog"));
   delete c1;

std::cout<<"F\n";

////////////////////////////////////// plot requested by G. Landsberg
if(Combine){

   c1 = new TCanvas("c1", "c1",600,600);
//   frame = new TH1D("frame", "frame", 1,90, 570);
//   frame->GetXaxis()->SetNdivisions(505);
//   frame->SetTitle("");
//   frame->SetStats(kFALSE);
//   frame->GetXaxis()->SetTitle("Mass (GeV)");
//   frame->GetYaxis()->SetTitle(Combine?"95% CL limit on #sigma/#sigma_{th}":"95% CL limit on #sigma (pb)");
//   frame->GetYaxis()->SetTitleOffset(1.40);
   frame->SetMaximum(20);
   frame->SetMinimum(1e-3);
   frame->Draw("AXIS");

   MuGraphMap["GMStau"       ]->SetLineColor(2);  MuGraphMap["GMStau"       ]->SetMarkerColor(2);   MuGraphMap["GMStau"       ]->SetLineWidth(2);   MuGraphMap["GMStau"       ]->SetLineStyle(1);  MuGraphMap["GMStau"       ]->SetMarkerStyle(22);
   TkGraphMap["GMStau"       ]->SetLineColor(1);  TkGraphMap["GMStau"       ]->SetMarkerColor(1);   TkGraphMap["GMStau"       ]->SetLineWidth(2);   TkGraphMap["GMStau"       ]->SetLineStyle(1);  TkGraphMap["GMStau"       ]->SetMarkerStyle(23);

   MuGraphMap["PPStau"       ]->SetLineColor(2);  MuGraphMap["PPStau"       ]->SetMarkerColor(2);   MuGraphMap["PPStau"       ]->SetLineWidth(2);   MuGraphMap["PPStau"       ]->SetLineStyle(2);  MuGraphMap["PPStau"       ]->SetMarkerStyle(26);
   TkGraphMap["PPStau"       ]->SetLineColor(1);  TkGraphMap["PPStau"       ]->SetMarkerColor(1);   TkGraphMap["PPStau"       ]->SetLineWidth(2);   TkGraphMap["PPStau"       ]->SetLineStyle(2);  TkGraphMap["PPStau"       ]->SetMarkerStyle(32);

   MuGraphMap["PPStau"     ]->Draw("LP");
   TkGraphMap["PPStau"     ]->Draw("LP");
   MuGraphMap["GMStau"     ]->Draw("LP");
   TkGraphMap["GMStau"     ]->Draw("LP");

   TLine* LineAtOne = new TLine(90,1,570,1);      LineAtOne->SetLineStyle(3);   LineAtOne->Draw();

   DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));

   LEGTk =/* !Combine ? */new TLegend(0.50,0.92-3*0.043,0.83,0.92)/* : new TLegend(0.45,0.15+4*0.043,0.80,0.15+7*0.043)*/;
   LEGTk->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGTk->SetTextSize(18); //font size
   LEGTk->SetFillColor(0); 
   LEGTk->SetFillStyle(0);
   LEGTk->SetBorderSize(0);
   LEGTk->SetHeader("Stau prod. (direct+indirect)");
   LEGTk->AddEntry(MuGraphMap["GMStau"     ], LegendFromType(MuPattern).c_str()                              ,"LP");
   LEGTk->AddEntry(TkGraphMap["GMStau"     ], LegendFromType(TkPattern).c_str()                              ,"LP");
   LEGTk->Draw();

   LEGTh =/* !Combine ? */new TLegend(0.50,0.92-3*0.043,0.83,0.92)/* : new TLegend(0.45,0.15+0*0.043,0.80,0.15+3*0.043)*/;
   LEGTh->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGTh->SetTextSize(18); //font size     
   LEGTh->SetFillColor(0); 
   LEGTh->SetFillStyle(0);
   LEGTh->SetBorderSize(0);
   LEGTh->SetHeader("Stau prod. (direct)");
   LEGTh->AddEntry(MuGraphMap["PPStau"     ], LegendFromType(MuPattern).c_str()                              ,"LP");
   LEGTh->AddEntry(TkGraphMap["PPStau"     ], LegendFromType(TkPattern).c_str()                              ,"LP");
   LEGTh->Draw();

   c1->SetLogy(true);
   SaveCanvas(c1, outpath, string("StauExclusionLog"));
   delete c1;
}
//////////////////////////////////////


   /*
   c1 = new TCanvas("c1", "c1",600,600);

   TMultiGraph* MGMO = new TMultiGraph();
   if(!Combine) {
   MGMO->Add(ThGraphMap["Gluino_f10" ]     ,"L");
   MGMO->Add(ThGraphMap["Stop"       ]     ,"L");
   //MGMO->Add(ThGraphMap["DY_Q1o3"    ]     ,"L");
   //MGMO->Add(ThGraphMap["DY_Q2o3"    ]     ,"L");
   }

   //TGraph* DYQ1o3ThLeg = (TGraph*) ThGraphMap["DY_Q1o3"        ]->Clone("DYQ1o3ThLeg");
   //DYQ1o3ThLeg->SetFillColor(ThErrorMap["DY_Q1o3"]->GetFillColor());
   //LQLEGTh->AddEntry(DYQ1o3ThLeg   ,"|Q| = e/3   (LO)" ,"LF");
   //TGraph* DYQ2o3ThLeg = (TGraph*) ThGraphMap["DY_Q2o3"        ]->Clone("DYQ2o3ThLeg");
   //DYQ2o3ThLeg->SetFillColor(ThErrorMap["DY_Q2o3"]->GetFillColor());
   //LQLEGTh->AddEntry(DYQ2o3ThLeg   ,"|Q| = 2e/3   (LO)" ,"LF");
   //LQLEGTh->Draw();


   MGMO->Add(MOGraphMap["Gluino_f10" ]     ,"LP");
   MGMO->Add(MOGraphMap["Gluino_f50" ]     ,"LP");
   MGMO->Add(MOGraphMap["Gluino_f100"]     ,"LP");
   MGMO->Add(MOGraphMap["Stop"       ]     ,"LP");
   //MGMO->Add(MOGraphMap["DY_Q1o3"    ]     ,"LP");
   //MGMO->Add(MOGraphMap["DY_Q2o3"    ]     ,"LP");

   MGMO->Draw("A");
   if(!Combine) {
     ThErrorMap["Gluino_f10"]->Draw("f");
     ThErrorMap["Stop"      ]->Draw("f");
     //ThErrorMap["DY_Q1o3"   ]->Draw("f");
     //ThErrorMap["DY_Q2o3"   ]->Draw("f");
   }else{
      TLine* LineAtOne = new TLine(50,1,1550,1);      LineAtOne->SetLineStyle(3);   LineAtOne->Draw();
   }

   MGMO->Draw("same");
   MGMO->SetTitle("");
   MGMO->GetXaxis()->SetTitle("Mass (GeV)");
   MGMO->GetYaxis()->SetTitle(Combine?"95% CL limit on #sigma/#sigma_{th}":"95% CL limit on #sigma (pb)");
   MGMO->GetYaxis()->SetTitleOffset(1.40);
   MGMO->GetYaxis()->SetRangeUser(PlotMinScale,PlotMaxScale);
   MGMO->GetXaxis()->SetRangeUser(50,1550);
   
   DrawPreliminary(LegendFromType(MOPattern).c_str(), 8.0, IntegratedLuminosityFromE(8.0), true);

   TLegend* LEGMO = !Combine ? new TLegend(0.50,0.92-4*0.043,0.83,0.92) : new TLegend(0.55,0.25,0.80,0.25+4*0.043);
   LEGMO->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGMO->SetTextSize(18); //font size
   LEGMO->SetFillColor(0); 
   LEGMO->SetFillStyle(0);
   LEGMO->SetBorderSize(0);
   LEGMO->AddEntry(MOGraphMap["Gluino_f100" ], "gluino; 100% #tilde{g}g"            ,"LP");
   LEGMO->AddEntry(MOGraphMap["Gluino_f50" ], "gluino; 50% #tilde{g}g"            ,"LP");
   LEGMO->AddEntry(MOGraphMap["Gluino_f10" ], "gluino; 10% #tilde{g}g"            ,"LP");
   LEGMO->AddEntry(MOGraphMap["Stop"       ], "stop"                              ,"LP");
   //LEGMO->AddEntry(TkGraphMap["DY_Q1o3"    ], "|Q| = e/3"            ,"LP");
   //LEGMO->AddEntry(TkGraphMap["DY_Q2o3"    ], "|Q| = 2e/3"            ,"LP");
   LEGMO->Draw();

   TLegend* MOLEGTh = new TLegend(0.15,0.92-(1+2)*0.043,0.50,0.92);
   MOLEGTh->SetTextFont(43); //give the font size in pixel (instead of fraction)
   MOLEGTh->SetTextSize(18); //font size
   if(!Combine) {
     MOLEGTh->SetHeader("Theoretical Prediction");
     MOLEGTh->SetFillColor(0);
     //MOLEGTh->SetFillStyle(0);
     MOLEGTh->SetBorderSize(0);

     TGraph* GlThLeg = (TGraph*) ThGraphMap["Gluino_f10"]->Clone("GluinoThLeg");
     GlThLeg->SetFillColor(ThErrorMap["Gluino_f10"]->GetFillColor());
     MOLEGTh->AddEntry(GlThLeg, "gluino (NLO+NLL)" ,"LF");
     TGraph* StThLeg = (TGraph*) ThGraphMap["Stop"      ]->Clone("StopThLeg");
     StThLeg->SetFillColor(ThErrorMap["Gluino_f10"]->GetFillColor());
     MOLEGTh->AddEntry(StThLeg   ,"stop (NLO+NLL)" ,"LF");
     //TGraph* DYQ1o3ThLeg = (TGraph*) ThGraphMap["DY_Q1o3"        ]->Clone("DYQ1o3ThLeg");
     //DYQ1o3ThLeg->SetFillColor(ThErrorMap["DY_Q1o3"]->GetFillColor());
     //MOLEGTh->AddEntry(DYQ1o3ThLeg   ,"|Q| = e/3 (LO)" ,"LF");
     //TGraph* DYQ2o3ThLeg = (TGraph*) ThGraphMap["DY_Q2o3"        ]->Clone("DYQ2o3ThLeg");
     //DYQ2o3ThLeg->SetFillColor(ThErrorMap["DY_Q2o3"]->GetFillColor());
     //MOLEGTh->AddEntry(DYQ2o3ThLeg   ,"|Q| = 2e/3 (LO)" ,"LF");
     MOLEGTh->Draw();
   }

   c1->SetLogy(true);
   SaveCanvas(c1, outpath, string("MOExclusionLog"));
   delete c1;








   /////////////////////////////// LQ Analysis
   TLegend* LQLEGTh = new TLegend(0.15,0.92-(1+2)*0.043,0.50,0.92);
   LQLEGTh->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LQLEGTh->SetTextSize(18); //font size
   c1 = new TCanvas("c1", "c1",600,600);
   if(!Combine) {
   LQLEGTh->SetHeader("Theoretical Prediction");
   LQLEGTh->SetFillColor(0);
   //LQLEGTh->SetFillStyle(0);
   LQLEGTh->SetBorderSize(0);

   TGraph* DYQ1o3ThLeg = (TGraph*) ThGraphMap["DY_Q1o3"        ]->Clone("DYQ1o3ThLeg");
   DYQ1o3ThLeg->SetFillColor(ThErrorMap["DY_Q1o3"]->GetFillColor());
   LQLEGTh->AddEntry(DYQ1o3ThLeg   ,"|Q| = e/3   (LO)" ,"LF");
   TGraph* DYQ2o3ThLeg = (TGraph*) ThGraphMap["DY_Q2o3"        ]->Clone("DYQ2o3ThLeg");
   DYQ2o3ThLeg->SetFillColor(ThErrorMap["DY_Q2o3"]->GetFillColor());
   LQLEGTh->AddEntry(DYQ2o3ThLeg   ,"|Q| = 2e/3   (LO)" ,"LF");
   LQLEGTh->Draw();
   }

   TMultiGraph* MGLQ = new TMultiGraph();
   if(!Combine) {
   MGLQ->Add(ThGraphMap["DY_Q1o3"    ]     ,"L");
   MGLQ->Add(ThGraphMap["DY_Q2o3"    ]     ,"L");
   }

   MGLQ->Add(LQGraphMap["DY_Q1o3"    ]     ,"LP");
   MGLQ->Add(LQGraphMap["DY_Q2o3"    ]     ,"LP");

   MGLQ->Draw("A");
   if(!Combine) {
   ThErrorMap["DY_Q1o3"   ]->Draw("f");
   ThErrorMap["DY_Q2o3"   ]->Draw("f");
   }else{
      TLine* LineAtOne = new TLine(75,1,625,1);      LineAtOne->SetLineStyle(3);   LineAtOne->Draw();
   }

   MGLQ->Draw("same");
   MGLQ->SetTitle("");
   MGLQ->GetXaxis()->SetTitle("Mass (GeV)");
   MGLQ->GetYaxis()->SetTitle(Combine?"95% CL limit on #sigma/#sigma_{th}":"95% CL limit on #sigma (pb)");
   MGLQ->GetYaxis()->SetTitleOffset(1.40);
   MGLQ->GetYaxis()->SetRangeUser(PlotMinScale,PlotMaxScale);
   MGLQ->GetXaxis()->SetRangeUser(75,625);

   DrawPreliminary(LegendFromType(LQPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS), true);

   TLegend* LEGLQ = !Combine ? new TLegend(0.50,0.92-2*0.043,0.83,0.92) : new TLegend(0.20,0.88-2*0.043,0.50,0.88);
   LEGLQ->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGLQ->SetTextSize(18); //font size
//   LEGLQ->SetHeader("Q<1");
   LEGLQ->SetFillColor(0); 
   LEGLQ->SetFillStyle(0);
   LEGLQ->SetBorderSize(0);
   LEGLQ->AddEntry(TkGraphMap["DY_Q1o3"    ], "|Q| = e/3"            ,"LP");
   LEGLQ->AddEntry(TkGraphMap["DY_Q2o3"    ], "|Q| = 2e/3"            ,"LP");
   if(!Combine) LQLEGTh->Draw();

   LEGLQ->Draw();
   c1->SetLogy(true);
   SaveCanvas(c1, outpath, string("LQExclusionLog"));
   delete c1;


   c1 = new TCanvas("c1", "c1",600,600);
   TMultiGraph* MGHQ = new TMultiGraph();

   if(!Combine) {
     MGHQ->Add(ThGraphMap["DY_Q1" ]      ,"L");
     MGHQ->Add(ThGraphMap["DY_Q2" ]      ,"L");
     MGHQ->Add(ThGraphMap["DY_Q3" ]      ,"L");
     MGHQ->Add(ThGraphMap["DY_Q4" ]      ,"L");
     MGHQ->Add(ThGraphMap["DY_Q5" ]      ,"L");
     MGHQ->Add(ThGraphMap["DY_Q6" ]      ,"L");
     MGHQ->Add(ThGraphMap["DY_Q7" ]      ,"L");
     MGHQ->Add(ThGraphMap["DY_Q8" ]      ,"L");
   }
   MGHQ->Add(HQGraphMap["DY_Q1" ]      ,"LP");
   MGHQ->Add(HQGraphMap["DY_Q2" ]      ,"LP");
   MGHQ->Add(HQGraphMap["DY_Q3" ]      ,"LP");
   MGHQ->Add(HQGraphMap["DY_Q4" ]      ,"LP");
   MGHQ->Add(HQGraphMap["DY_Q5" ]      ,"LP");
   MGHQ->Add(HQGraphMap["DY_Q6" ]      ,"LP");
   MGHQ->Add(HQGraphMap["DY_Q7" ]      ,"LP");
   MGHQ->Add(HQGraphMap["DY_Q8" ]      ,"LP");
   MGHQ->Draw("A");
   
   if(!Combine) {
     ThErrorMap["DY_Q1"]->Draw("f");
     ThErrorMap["DY_Q2"]->Draw("f");
     ThErrorMap["DY_Q3"]->Draw("f");
     ThErrorMap["DY_Q4"]->Draw("f");
     ThErrorMap["DY_Q5"]->Draw("f");
     ThErrorMap["DY_Q6"]->Draw("f");
     ThErrorMap["DY_Q7"]->Draw("f");
     ThErrorMap["DY_Q8"]->Draw("f");
   }else{
      TLine* LineAtOne = new TLine(50,1,1050,1);      LineAtOne->SetLineStyle(3);   LineAtOne->Draw();
   }

   
   MGHQ->Draw("same");
   MGHQ->SetTitle("");
   MGHQ->GetXaxis()->SetTitle("Mass (GeV)");
   MGHQ->GetYaxis()->SetTitle(Combine?"95% CL limit on #sigma/#sigma_{th}":"95% CL limit on #sigma (pb)");
   MGHQ->GetYaxis()->SetTitleOffset(1.40);
   MGHQ->GetYaxis()->SetRangeUser(PlotMinScale,PlotMaxScale);
   //MGHQ->GetYaxis()->SetRangeUser(PlotMinScale,100);
   MGHQ->GetXaxis()->SetRangeUser(50,1050);

   DrawPreliminary(LegendFromType(HQPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS), true);
   TLegend* LEGHQ = !Combine ? new TLegend(0.62,0.92-0.043-8*0.043,0.83,0.92-0.043) : new TLegend(0.55,0.35,0.80,0.35+6*0.043);
//   TLegend* LEGHQ = !Combine ? new TLegend(0.62,0.92-5*0.043,0.83,0.92) : new TLegend(0.55,0.35,0.80,0.35+6*0.043);
   LEGHQ->SetTextFont(43); //give the font size in pixel (instead of fraction)
   LEGHQ->SetTextSize(18); //font size

   LEGHQ->SetFillColor(0); 
   //LEGHQ->SetFillStyle(0);
   LEGHQ->SetBorderSize(0);
   LEGHQ->AddEntry(HQGraphMap["DY_Q1"] , "|Q| = 1e "    ,"LP");
   LEGHQ->AddEntry(HQGraphMap["DY_Q2"] , "|Q| = 2e "    ,"LP");
   LEGHQ->AddEntry(HQGraphMap["DY_Q3"] , "|Q| = 3e "    ,"LP");
   LEGHQ->AddEntry(HQGraphMap["DY_Q4"] , "|Q| = 4e "    ,"LP");
   LEGHQ->AddEntry(HQGraphMap["DY_Q5"] , "|Q| = 5e "    ,"LP");
   LEGHQ->AddEntry(HQGraphMap["DY_Q6"] , "|Q| = 6e "    ,"LP");
   LEGHQ->AddEntry(HQGraphMap["DY_Q7"] , "|Q| = 7e "    ,"LP");
   LEGHQ->AddEntry(HQGraphMap["DY_Q8"] , "|Q| = 8e "    ,"LP");

   TLegend* HQLEGTh = new TLegend(0.35,0.92-(1+8)*0.043,0.57,0.92);
//   TLegend* HQLEGTh = new TLegend(0.3,0.92-(1+5)*0.043,0.57,0.92);
   HQLEGTh->SetTextFont(43); //give the font size in pixel (instead of fraction)
   HQLEGTh->SetTextSize(18); //font size
   if(!Combine){
   HQLEGTh->SetHeader("Theoretical Prediction");
   HQLEGTh->SetFillColor(0);
   HQLEGTh->SetFillStyle(0);
   HQLEGTh->SetBorderSize(0);
   
   TGraph* Q1ThLeg = (TGraph*) ThGraphMap["DY_Q1"]->Clone("HSCPQ1ThLeg");
   Q1ThLeg->SetFillColor(ThErrorMap["DY_Q1"]->GetFillColor());
   HQLEGTh->AddEntry(Q1ThLeg, "|Q| = 1e (LO)" ,"LF");

   TGraph* Q2ThLeg = (TGraph*) ThGraphMap["DY_Q2"]->Clone("HSCPQ2ThLeg");
   Q2ThLeg->SetFillColor(ThErrorMap["DY_Q2"]->GetFillColor());
   HQLEGTh->AddEntry(Q2ThLeg, "|Q| = 2e (LO)" ,"LF");

   TGraph* Q3ThLeg = (TGraph*) ThGraphMap["DY_Q3"]->Clone("HSCPQ3ThLeg");
   Q3ThLeg->SetFillColor(ThErrorMap["DY_Q3"]->GetFillColor());
   HQLEGTh->AddEntry(Q3ThLeg, "|Q| = 3e (LO)" ,"LF");

   TGraph* Q4ThLeg = (TGraph*) ThGraphMap["DY_Q4"]->Clone("HSCPQ4ThLeg");
   Q4ThLeg->SetFillColor(ThErrorMap["DY_Q4"]->GetFillColor());
   HQLEGTh->AddEntry(Q4ThLeg, "|Q| = 4e (LO)" ,"LF");

   TGraph* Q5ThLeg = (TGraph*) ThGraphMap["DY_Q5"]->Clone("HSCPQ5ThLeg");
   Q5ThLeg->SetFillColor(ThErrorMap["DY_Q5"]->GetFillColor());
   HQLEGTh->AddEntry(Q5ThLeg, "|Q| = 5e (LO)" ,"LF");

   TGraph* Q6ThLeg = (TGraph*) ThGraphMap["DY_Q6"]->Clone("HSCPQ6ThLeg");
   Q6ThLeg->SetFillColor(ThErrorMap["DY_Q6"]->GetFillColor());
   HQLEGTh->AddEntry(Q6ThLeg, "|Q| = 6e (LO)" ,"LF");

   TGraph* Q7ThLeg = (TGraph*) ThGraphMap["DY_Q7"]->Clone("HSCPQ7ThLeg");
   Q7ThLeg->SetFillColor(ThErrorMap["DY_Q7"]->GetFillColor());
   HQLEGTh->AddEntry(Q7ThLeg, "|Q| = 7e (LO)" ,"LF");

   TGraph* Q8ThLeg = (TGraph*) ThGraphMap["DY_Q8"]->Clone("HSCPQ8ThLeg");
   Q8ThLeg->SetFillColor(ThErrorMap["DY_Q8"]->GetFillColor());
   HQLEGTh->AddEntry(Q8ThLeg, "|Q| = 8e (LO)" ,"LF");

   HQLEGTh->Draw();
   }

   LEGHQ->Draw();
   c1->SetLogy(true);
   SaveCanvas(c1, outpath, string("HQExclusionLog"));
   delete c1;
   */


   return;
}

TGraph* CheckSignalUncertainty(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, std::vector<stSample>& modelSample){
  int TypeMode = TypeFromPattern(InputPattern);
  string prefix = "BUG";
  switch(TypeMode){
  case 0: prefix   = "Tk"; break;
  case 2: prefix   = "Mu"; break;
  case 3: prefix   = "Mo"; break;
  case 4: prefix   = "HQ"; break;
  case 5: prefix   = "LQ"; break;
  }

   unsigned int N   = 0;

   double* Mass      = new double   [modelSample.size()];  double* MassErr      = new double   [modelSample.size()];
   double* SystP     = new double   [modelSample.size()];  double* SystErrP     = new double   [modelSample.size()];
   double* SystI     = new double   [modelSample.size()];  double* SystErrI     = new double   [modelSample.size()];
   double* SystM     = new double   [modelSample.size()];  double* SystErrM     = new double   [modelSample.size()];
   double* SystH     = new double   [modelSample.size()];  double* SystErrH     = new double   [modelSample.size()];
   double* SystPU    = new double   [modelSample.size()];  double* SystErrPU    = new double   [modelSample.size()];
   double* SystT     = new double   [modelSample.size()];  double* SystErrT     = new double   [modelSample.size()];
   double* SystTr    = new double   [modelSample.size()];  double* SystErrTr    = new double   [modelSample.size()];
   double* SystRe    = new double   [modelSample.size()];  double* SystErrRe    = new double   [modelSample.size()];
   double* SystMB    = new double   [modelSample.size()];  double* SystErrMB    = new double   [modelSample.size()];
   double* SystTotal = new double   [modelSample.size()];  double* SystErrTotal = new double   [modelSample.size()];
   double* SystTotal2 = new double   [modelSample.size()]; 


   for(unsigned int s=0;s<modelSample.size();s++){
      if(modelSample[s].Type!=2)continue;
      bool IsNeutral = (modelSample[s].ModelName().find("N")!=std::string::npos);
      if(TypeMode!=0 && IsNeutral)continue;
      stAllInfo tmp(InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR + "/"+modelSample[s].Name+".txt");
      if(tmp.Eff==0) continue;

      Mass[N]        = tmp.Mass;                            MassErr[N]        = 0.0;
      SystP[N]       = (tmp.Eff_SYSTP  - tmp.Eff)/tmp.Eff;  SystErrP[N]       = (tmp.EffE_SYSTP)/tmp.Eff;
      SystI[N]       = (tmp.Eff_SYSTI  - tmp.Eff)/tmp.Eff;  SystErrI[N]       = (tmp.EffE_SYSTI)/tmp.Eff;
      if(TypeMode==5){
         if(modelSample[s].ModelName().find("1o3")!=string::npos) SystI[N]=-0.25; SystErrP[N]       = 0;
         if(modelSample[s].ModelName().find("2o3")!=string::npos) SystI[N]=-0.10; SystErrI[N]       = 0;
      }
      SystM[N]       = (tmp.Eff_SYSTM  - tmp.Eff)/tmp.Eff;  SystErrM[N]       = (tmp.EffE_SYSTM)/tmp.Eff;
      SystH[N]       = (tmp.Eff_SYSTHUp - tmp.Eff)/tmp.Eff;  SystErrH[N]       = (tmp.EffE_SYSTHUp)/tmp.Eff;
      SystPU[N]      = (tmp.Eff_SYSTPU - tmp.Eff)/tmp.Eff;  SystErrPU[N]      = (tmp.EffE_SYSTPU)/tmp.Eff;
      SystT[N]       = (tmp.Eff_SYSTT  - tmp.Eff)/tmp.Eff;  SystErrT[N]       = (tmp.EffE_SYSTT )/tmp.Eff;
      SystRe[N]      = -0.02; SystErrRe[N]      = 0.0;
      SystMB[N]=0.;  SystErrMB[N]=0.0;
//      if((modelSample[s].ModelName().find("Q2")!=string::npos && modelSample[s].ModelName().find("Q2o3")==string::npos) || modelSample[s].ModelName().find("Q3")!=string::npos || modelSample[s].ModelName().find("Q4")!=string::npos || modelSample[s].ModelName().find("Q5")!=string::npos) SystMB[N]=-0.2;
      if(modelSample[s].ModelName().find("DY_Q")!=string::npos && modelSample[s].ModelName().find("o3")==string::npos) SystMB[N]=-0.2;

      
      if(SQRTS==7) {
	if(modelSample[s].ModelName().find("1o3")!=string::npos) SystTr[N] = -1*sqrt(0.15*0.15 + 0.02*0.02 + 0.05*0.05);
	else if(modelSample[s].ModelName().find("2o3")!=string::npos) SystTr[N] = -1*sqrt(0.03*0.03 + 0.02*0.02 + 0.05*0.05);
	else if(IsNeutral) SystTr[N] = -0.05;
	else SystTr[N] = -1*sqrt(0.05*0.05 + 0.02*0.02 + 0.02*0.02); 
        SystErrTr[N] = 0.0;
      }else if(SQRTS==8){
	if(IsNeutral) SystTr[N] = -0.01;
	else if(modelSample[s].ModelName().find("1o3")!=string::npos) SystTr[N] = -1*sqrt(0.15*0.15 + 0.04*0.04 + 0.05*0.05);
	else if(modelSample[s].ModelName().find("2o3")!=string::npos) SystTr[N] = -1*sqrt(0.03*0.03 + 0.04*0.04 + 0.05*0.05);
	else SystTr[N] = -1*sqrt(0.05*0.05 + 0.04*0.04 + 0.01*0.01);
        SystErrTr[N] = 0.0;
      }else{
	if(IsNeutral) SystTr[N] = -0.13;
//	else if(modelSample[s].ModelName().find("1o3")!=string::npos) SystTr[N] = -1*sqrt(0.15*0.15 + 0.04*0.04 + 0.05*0.05);
//	else if(modelSample[s].ModelName().find("2o3")!=string::npos) SystTr[N] = -1*sqrt(0.03*0.03 + 0.04*0.04 + 0.05*0.05);
	else SystTr[N] = -0.13;
        SystErrTr[N] = 0.0;
      }

//      double Ptemp=max(SystP[N], 0.0), Itemp=max(SystI[N], 0.0), PUtemp=max(SystPU[N], 0.0), Ttemp=max(SystT[N], 0.0);
      double Ptemp=SystP[N], Itemp=SystI[N], Mtemp=SystM[N], Htemp=SystH[N], PUtemp=SystPU[N], Ttemp=SystT[N];
      SystTotal[N] = -1*sqrt(Ptemp*Ptemp + Itemp*Itemp + Mtemp*Mtemp + Htemp*Htemp + PUtemp*PUtemp + Ttemp*Ttemp + SystTr[N]*SystTr[N] + SystRe[N]*SystRe[N] + SystMB[N]*SystMB[N]);
//      SystErrTotal[N] = sqrt(pow(SystErrP[N],2) + pow(SystErrI[N],2) + pow(SystErrPU[N],2) + pow(SystErrT[N],2) + pow(SystErrTr[N],2) + pow(SystErrRe[N],2) + pow(SystErrMB[N],2) );
      SystErrTotal[N] = (tmp.EffE)/tmp.Eff;
      SystTotal2[N] = -1*SystTotal[N]; 

      if(TypeMode==0 || TypeMode==5)fprintf(pFile, "%30s   %7.3f --> %7.3f  |  %7.3f  | %7.3f  | %7.3f | %7.3f  | %7.3f"        ,modelSample[N].Name.c_str(), tmp.Eff, SystP[N], SystI[N], SystM[N], SystH[N], SystPU[N]          , SystTotal[N]);  
      else                          fprintf(pFile, "%30s   %7.3f --> %7.3f  |  %7.3f  | %7.3f  | %7.3f | %7.3f  | %7.3f | %7.3f",modelSample[N].Name.c_str(), tmp.Eff, SystP[N], SystI[N], SystM[N], SystH[N], SystPU[N], SystT[N], SystTotal[N]);


      //Check if they are variated samples for this point
      for(unsigned int sv=0;sv<samples.size();sv++){
         if(samples[sv].Type!=3)continue;
         if(samples[sv].Name.find(modelSample[s].Name)!=0)continue; //we expect to have the beginning of the name being identical
         if(samples[sv].Name.length()!=modelSample[s].Name.length()+3) continue; // variated samples are exactly 3char longer
         stAllInfo tmpVaried(InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR + "/"+samples[sv].Name+".txt");
         if(tmp.Mass<=0) continue;
//         fprintf(pFile, "%20s%10s    %7.3f --> RelDiff=%7.3f\n", "   Variation --> ",samples[sv].Name.c_str()+modelSample[s].Name.length()+1, tmpVaried.Eff, (tmpVaried.Eff-tmp.Eff)/tmp.Eff);
         fprintf(pFile, " | %7.3f (%2s)", (tmpVaried.Eff - tmp.Eff)/tmp.Eff, samples[sv].Name.c_str()+modelSample[s].Name.length()+1);
      }

      fprintf(pFile, "\n");
      N++;
   }

   TGraph* graphSystP = NULL;
   TGraph* graphSystI = NULL;
   TGraph* graphSystM = NULL;
   TGraph* graphSystH = NULL;
   TGraph* graphSystPU = NULL;
   TGraph* graphSystT = NULL;
   TGraph* graphSystTr = NULL;
   TGraph* graphSystRe = NULL;
   TGraph* graphSystMB = NULL;
   TGraph* graphSystTotal = NULL;
   TGraph* graphSystTotal2 = NULL;

   if(N>0) {
     TCanvas* c2 = new TCanvas("c2", "c2",600,600);
     c2->SetLeftMargin(0.15);

//     graphSystP = new TGraphErrors(N,Mass,SystP, MassErr, SystErrP);
//     graphSystI = new TGraphErrors(N,Mass,SystI, MassErr,SystErrI);
//     graphSystPU = new TGraphErrors(N,Mass,SystPU, MassErr,SystErrPU);
//     graphSystT = new TGraphErrors(N,Mass,SystT, MassErr,SystErrT);
//     graphSystTr = new TGraphErrors(N,Mass,SystTr, MassErr,SystErrTr);
//     graphSystRe = new TGraphErrors(N,Mass,SystRe, MassErr,SystErrRe);
//     graphSystMB = new TGraphErrors(N,Mass,SystMB, MassErr,SystErrMB);
//     graphSystTotal = new TGraphErrors(N,Mass,SystTotal, MassErr, SystErrTotal);

     graphSystP = new TGraph(N,Mass,SystP);//, MassErr, SystErrP);
     graphSystI = new TGraph(N,Mass,SystI);//, MassErr,SystErrI);
     graphSystM = new TGraph(N,Mass,SystM);//, MassErr,SystErrI);
     graphSystH = new TGraph(N,Mass,SystH);//, MassErr,SystErrI);
     graphSystPU = new TGraph(N,Mass,SystPU);//, MassErr,SystErrPU);
     graphSystT = new TGraph(N,Mass,SystT);//, MassErr,SystErrT);
     graphSystTr = new TGraph(N,Mass,SystTr);//, MassErr,SystErrTr);
     graphSystRe = new TGraph(N,Mass,SystRe);//, MassErr,SystErrRe);
     graphSystMB = new TGraph(N,Mass,SystMB);//, MassErr,SystErrMB);
     graphSystTotal = new TGraphErrors(N,Mass,SystTotal, MassErr, SystErrTotal);
     graphSystTotal2 = new TGraph(N,Mass,SystTotal2);//, MassErr, SystErrTotal);
     TMultiGraph* SystGraphs = new TMultiGraph();

     graphSystTotal->GetYaxis()->SetTitle("CrossSection ( pb )");
     graphSystTotal->SetLineColor(Color[0]);  graphSystTotal->SetMarkerColor(Color[0]);   graphSystTotal->SetMarkerStyle(20);    graphSystTotal->SetLineWidth(2);
     graphSystP->SetLineColor(Color[1]);      graphSystP->SetMarkerColor(Color[1]);       graphSystP->SetMarkerStyle(Marker[1]); graphSystP->SetLineWidth(2);
     graphSystI->SetLineColor(Color[2]);      graphSystI->SetMarkerColor(Color[2]);       graphSystI->SetMarkerStyle(Marker[2]); graphSystI->SetLineWidth(2);
     graphSystM->SetLineColor(Color[8]);      graphSystM->SetMarkerColor(Color[8]);       graphSystM->SetMarkerStyle(Marker[8]); graphSystM->SetLineWidth(2);
     graphSystH->SetLineColor(Color[9]);      graphSystH->SetMarkerColor(Color[9]);       graphSystH->SetMarkerStyle(Marker[9]); graphSystH->SetLineWidth(2);
     graphSystPU->SetLineColor(Color[3]);     graphSystPU->SetMarkerColor(Color[3]);      graphSystPU->SetMarkerStyle(Marker[3]);graphSystPU->SetLineWidth(2);
     graphSystT->SetLineColor(Color[4]);      graphSystT->SetMarkerColor(Color[4]);       graphSystT->SetMarkerStyle(Marker[4]); graphSystT->SetLineWidth(2);
     graphSystTr->SetLineColor(Color[5]);     graphSystTr->SetMarkerColor(Color[5]);      graphSystTr->SetMarkerStyle(Marker[5]);graphSystTr->SetLineWidth(2);
     graphSystRe->SetLineColor(Color[6]);     graphSystRe->SetMarkerColor(Color[6]);      graphSystRe->SetMarkerStyle(Marker[6]);graphSystRe->SetLineWidth(2);
     graphSystMB->SetLineColor(Color[7]);     graphSystMB->SetMarkerColor(Color[7]);      graphSystMB->SetMarkerStyle(Marker[7]);graphSystMB->SetLineWidth(2);

     SystGraphs->Add(graphSystP,"LP0");
     SystGraphs->Add(graphSystTr,"LP");
     SystGraphs->Add(graphSystRe,"LP");
     if(TypeMode!=3)SystGraphs->Add(graphSystI,"LP0");
     if(TypeMode!=3)SystGraphs->Add(graphSystH,"LP0");
     if(TypeMode!=3)SystGraphs->Add(graphSystM,"LP0");
     SystGraphs->Add(graphSystPU,"LP0");
     if(TypeMode!=0 && TypeMode!=5)SystGraphs->Add(graphSystT,"LP0");
     if(TypeMode==4) SystGraphs->Add(graphSystMB,"LP");
     SystGraphs->Add(graphSystTotal,"P");

     SystGraphs->Draw("A");
     SystGraphs->SetTitle("");
     SystGraphs->GetXaxis()->SetTitle("Mass (GeV)");
     SystGraphs->GetYaxis()->SetTitle("Relative Change in Efficiency");
     SystGraphs->GetYaxis()->SetTitleOffset(1.40);
     SystGraphs->GetYaxis()->SetRangeUser(-0.45, 0.45);
     SystGraphs->GetYaxis()->SetNdivisions(510, "X");

     TLegend* LEG = new TLegend(0.30,0.80,0.80,0.65);
     LEG->SetFillColor(0);
     LEG->SetFillStyle(0);
     LEG->SetBorderSize(0);
     LEG->AddEntry(graphSystTr,  "Trigger" ,"LP");
     LEG->AddEntry(graphSystRe,  "Reconstruction" ,"LP");
     if(TypeMode==4)LEG->AddEntry(graphSystMB,  "MB" ,"LP");
     LEG->AddEntry(graphSystP,  "P" ,"LP");
     if(TypeMode!=3)LEG->AddEntry(graphSystI,  "dE/dx Ias" ,"LP");
     if(TypeMode!=3)LEG->AddEntry(graphSystM,  "dE/dx Mass" ,"LP");
     if(TypeMode!=3)LEG->AddEntry(graphSystH,  "dE/dx HIP" ,"LP");
     LEG->AddEntry(graphSystPU,  "Pile Up" ,"LP");
     if(TypeMode!=0 && TypeMode!=5)LEG->AddEntry(graphSystT,  "1/#beta" ,"LP");
     LEG->AddEntry(graphSystTotal,  "Total" ,"P");
     LEG->SetNColumns(2);
     LEG->Draw();
     c2->SetLogy(false);
     c2->SetGridy(false);

   DrawPreliminary(LegendFromType(InputPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
   SaveCanvas(c2,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", string(prefix+ ModelName + "Uncertainty"));
   delete c2;
   //delete SystGraphs;
   }

   return graphSystTotal2;
}


TGraph* MakePlot(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, int XSectionType, std::vector<stSample>& modelSamples, double& LInt){
   std::vector<int> signalPoints;
   for(unsigned int i=0;i<modelSamples.size();i++) if(XSectionType==0 || stAllInfo(InputPattern+""+SHAPESTRING+EXCLUSIONDIR+"/" + modelSamples[i].Name +".txt").XSec_Exp<1E10) {
     //Skip 100GeV for DY Q=7 and Q=8
     if(XSectionType>0 && (ModelName.find("DY_Q7")!=string::npos || ModelName.find("DY_Q8")!=string::npos) && stAllInfo(InputPattern+""+SHAPESTRING+EXCLUSIONDIR+"/" + modelSamples[i].Name +".txt").Mass==100.0 )continue;
     signalPoints.push_back(i);
   }
   unsigned int N   = signalPoints.size();

   double* Mass     = new double   [signalPoints.size()];
   double* XSecTh   = new double   [signalPoints.size()];
   double* XSecObs  = new double   [signalPoints.size()];
   double* XSecExp  = new double   [signalPoints.size()];
   stAllInfo* Infos = new stAllInfo[signalPoints.size()];

   bool FileFound=false;

//   int I=0;
   for(unsigned int i=0;i<signalPoints.size();i++){
     Infos     [i] = stAllInfo(InputPattern+""+SHAPESTRING+EXCLUSIONDIR+"/" + modelSamples[signalPoints[i]].Name +".txt");
     if(Infos[i].Mass!=0) FileFound=true;
//     if(XSectionType>0 && Infos[i].XSec_Exp>1E10)continue;
     Mass      [i] = Infos[i].Mass;
     XSecTh    [i] = Infos[i].XSec_Th;
     XSecObs   [i] = Infos[i].XSec_Obs;
     XSecExp   [i] = Infos[i].XSec_Exp;
     LInt          = std::max(LInt, Infos[i].LInt);

     //printf("%i %s\n", (int)FileFound, (InputPattern+""+SHAPESTRING+EXCLUSIONDIR+"/" + modelSamples[signalPoints[i]].Name +".txt").c_str());
//     I++;
   }
//   N=I;
 

   if(XSectionType>0 && FileFound){
      //for(unsigned int i=0;i<N;i++)printf("%-18s %5.0f --> Pt>%+6.1f & I>%+5.3f & TOF>%+4.3f & M>%3.0f--> NData=%2.0f  NPred=%6.1E+-%6.1E  NSign=%6.1E (Eff=%3.2f) Local Significance %3.2f\n",ModelName.c_str(),Infos[i].Mass,Infos[i].WP_Pt,Infos[i].WP_I,Infos[i].WP_TOF,Infos[i].MassCut, Infos[i].NData, Infos[i].NPred, Infos[i].NPredErr, Infos[i].NSign, Infos[i].Eff, Infos[i].Significance);

     for(unsigned int i=0;i<signalPoints.size();i++){
       //for(unsigned int i=0;i<N;i++){
        if(Infos[i].WP_TOF==-1){fprintf(pFile,"%-20s & %4.0f & %6.0f & %5.3f & / & %4.0f & %6.3f $\\pm$ %6.3f & %2.0f & %4.3f & %6.1E & %6.1E & %6.1E & %3.2f \\\\\n", ModelName.c_str(), Infos[i].Mass,  Infos[i].WP_Pt,Infos[i].WP_I,Infos[i].MassCut, Infos[i].NPred, Infos[i].NPredErr, Infos[i].NData, Infos[i].Eff, Infos[i].XSec_Th,Infos[i].XSec_Exp, Infos[i].XSec_Obs, Infos[i].Significance);
        }else{                  fprintf(pFile,"%-20s & %4.0f & %6.0f & %5.3f & %4.3f & %4.0f & %6.3f $\\pm$ %6.3f & %2.0f & %4.3f & %6.1E & %6.1E & %6.1E & %3.2f \\\\\n", ModelName.c_str(), Infos[i].Mass,  Infos[i].WP_Pt,Infos[i].WP_I,Infos[i].WP_TOF,Infos[i].MassCut, Infos[i].NPred, Infos[i].NPredErr, Infos[i].NData, Infos[i].Eff, Infos[i].XSec_Th,Infos[i].XSec_Exp, Infos[i].XSec_Obs, Infos[i].Significance);
        }
        bool IsNeutral = (ModelName.find("N",0)<std::string::npos);
        if(Infos[i].WP_TOF==-1 && (ModelName=="GMSB Stau" || (int)Infos[i].Mass%200==0)) {
          fprintf(talkFile,"%-20s & %4.0f & %6.0f & %5.3f & / & %4.0f & %6.3f $\\pm$ %6.3f & %2.0f & %4.3f & %3.2f \\\\\n", ModelName.c_str(), Infos[i].Mass,  Infos[i].WP_Pt,Infos[i].WP_I,Infos[i].MassCut, Infos[i].NPred, Infos[i].NPredErr, Infos[i].NData, Infos[i].Eff, Infos[i].Significance);
          fprintf(talkFile, "\\hline\n");
        }
        if(Infos[i].WP_TOF!=-1 && !IsNeutral) {
          fprintf(talkFile,"%-20s & %4.0f & %6.0f & %5.3f & %4.3f & %4.0f & %6.3f $\\pm$ %6.3f & %2.0f & %4.3f %3.2f \\\\\n", ModelName.c_str(), Infos[i].Mass,  Infos[i].WP_Pt,Infos[i].WP_I,Infos[i].WP_TOF,Infos[i].MassCut, Infos[i].NPred, Infos[i].NPredErr, Infos[i].NData, Infos[i].Eff, Infos[i].Significance);
          fprintf(talkFile, "\\hline\n");
        }
      }
   }

   TGraph* graph = NULL;
   if(XSectionType==0)graph = new TGraph(N,Mass,XSecTh);
   if(XSectionType==1)graph = new TGraph(N,Mass,XSecExp);
   if(XSectionType==2)graph = new TGraph(N,Mass,XSecObs);
   graph->SetTitle("");
   graph->GetYaxis()->SetTitle("CrossSection ( pb )");
   graph->GetYaxis()->SetTitleOffset(1.40);
   return graph;
}


void printSummary(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, std::vector<stSample>& modelSamples){
  TypeMode = TypeFromPattern(InputPattern);
   for(unsigned int i=0;i<modelSamples.size();i++){
      string signal7TeV = modelSamples[i].Name; if(signal7TeV.find("_8TeV")!=string::npos) signal7TeV = signal7TeV.replace(signal7TeV.find("_8TeV"),5, "_7TeV");
      string signal8TeV = modelSamples[i].Name; if(signal8TeV.find("_7TeV")!=string::npos) signal8TeV = signal8TeV.replace(signal8TeV.find("_7TeV"),5, "_8TeV");
      string signal     = signal8TeV;           if(signal    .find("_8TeV")!=string::npos) signal     = signal    .replace(signal    .find("_8TeV"),5, "");
      stAllInfo Infos7(InputPattern+""+SHAPESTRING+"EXCLUSION7TeV"+"/" + signal7TeV +".txt");
      stAllInfo Infos8(InputPattern+""+SHAPESTRING+"EXCLUSION8TeV"+"/" + signal8TeV +".txt");
      stAllInfo InfosC(InputPattern+""+SHAPESTRING+"EXCLUSIONCOMB"+"/" + signal     +".txt");
      if(Infos7.Mass<=0 && Infos8.Mass<=0 && InfosC.Mass<=0)continue;
      if(Infos7.Eff<=0 && Infos8.Eff<=0)continue;
      double Mass = std::max(Infos7.Mass, Infos8.Mass);
      TString ModelNameTS =  ModelName.c_str();  ModelNameTS.ReplaceAll("_"," ");  ModelNameTS.ReplaceAll("8TeV",""); ModelNameTS.ReplaceAll("7TeV","");

      if(ModelNameTS.Contains("Stop")   && ((int)(Mass)/100)%2!=0)continue;
      if(ModelNameTS.Contains("Gluino") && ((int)(Mass)/100)%2!=1)continue;
      if(ModelNameTS.Contains("DY")     && ((int)(Mass)/100)%2!=0)continue;
      if(ModelNameTS.Contains("DC")                              )continue;

      char massCut[255];  if(TypeMode<3){sprintf(massCut,"$>%.0f$",Infos8.MassCut);}else{sprintf(massCut," - ");}
      char Results7[255]; if(Infos7.Mass>0 && TypeMode!=3){sprintf(Results7, "%6.2f & %6.2E & %6.2E & %6.2E", Infos7.Eff, Infos7.XSec_Th,Infos7.XSec_Obs, Infos7.XSec_Exp);}else{sprintf(Results7, "   -   &    -     &   -      &   -     ");}
      char Results8[255]; if(Infos8.Mass>0){sprintf(Results8, "%6.2f & %6.2E & %6.2E & %6.2E", Infos8.Eff, Infos8.XSec_Th,Infos8.XSec_Obs, Infos8.XSec_Exp);}else{sprintf(Results8, "   -    &    -     &    -     &   -     ");}
      char ResultsC[255]; if(InfosC.Mass>0 && TypeMode!=3){sprintf(ResultsC, "%6.2E & %6.2E", InfosC.XSec_Obs, InfosC.XSec_Exp);}
      else if(InfosC.Mass>0){sprintf(ResultsC, "%6.2E & %6.2E", Infos8.XSec_Obs/Infos8.XSec_Th, Infos8.XSec_Exp/Infos8.XSec_Th);}
      else{sprintf(ResultsC, "   -     &    -    ");}

//    char Results7[255]; if(Infos7.Mass>0){sprintf(Results7, "%10s & %10s & %10s & %10s", toLatexRounded(Infos7.Eff).c_str(), toLatexRounded(Infos7.XSec_Th).c_str(),toLatexRounded(Infos7.XSec_Obs).c_str(), toLatexRounded(Infos7.XSec_Exp).c_str());}else{sprintf(Results7, "%10s & %10s & %10s & %10s", "", "", "", "");}
//    char Results8[255]; if(Infos8.Mass>0){sprintf(Results8, "%10s & %10s & %10s & %10s", toLatexRounded(Infos8.Eff).c_str(), toLatexRounded(Infos8.XSec_Th).c_str(),toLatexRounded(Infos8.XSec_Obs).c_str(), toLatexRounded(Infos8.XSec_Exp).c_str());}else{sprintf(Results8, "%10s & %10s & %10s & %10s", "", "", "", "");}
//    char ResultsC[255]; if(InfosC.Mass>0){sprintf(ResultsC, "%10s & %10s", toLatexRounded(InfosC.XSec_Obs).c_str(), toLatexRounded(InfosC.XSec_Exp).c_str());}else{sprintf(ResultsC, "%10s & %10s", "", "");}

      fprintf(pFile,"%-20s & %4.0f & %-7s & %s & %s & %s\\\\\n", ModelNameTS.Data(), Mass, massCut, Results7, Results8, ResultsC);
   }
}



void printSummaryPaper(FILE* pFile, FILE* talkFile, string InputPattern, string ModelName, std::vector<stSample>& modelSamples){
   for(unsigned int i=0;i<modelSamples.size();i++){
     TypeMode = TypeFromPattern(InputPattern);
      string signal7TeV = modelSamples[i].Name; if(signal7TeV.find("_8TeV")!=string::npos) signal7TeV = signal7TeV.replace(signal7TeV.find("_8TeV"),5, "_7TeV");
      string signal8TeV = modelSamples[i].Name; if(signal8TeV.find("_7TeV")!=string::npos) signal8TeV = signal8TeV.replace(signal8TeV.find("_7TeV"),5, "_8TeV");
      string signal     = signal8TeV;           if(signal    .find("_8TeV")!=string::npos) signal     = signal    .replace(signal    .find("_8TeV"),5, "");
      stAllInfo Infos7(InputPattern+""+SHAPESTRING+"EXCLUSION7TeV"+"/" + signal7TeV +".txt");
      stAllInfo Infos8(InputPattern+""+SHAPESTRING+"EXCLUSION8TeV"+"/" + signal8TeV +".txt");
      stAllInfo InfosC(InputPattern+""+SHAPESTRING+"EXCLUSIONCOMB"+"/" + signal     +".txt");
      if(Infos7.Mass<=0 && Infos8.Mass<=0 && InfosC.Mass<=0)continue;
      if(Infos7.Eff<=0 && Infos8.Eff<=0)continue;
      double Mass = std::max(Infos7.Mass, Infos8.Mass);
      TString ModelNameTS =  ModelName.c_str();  ModelNameTS.ReplaceAll("8TeV",""); ModelNameTS.ReplaceAll("7TeV","");

      if((ModelNameTS.Contains("Gluino_f10") && !ModelNameTS.Contains("f100") && ((int)(Mass)/100)%4==3 && TypeMode==0) ||
	 (ModelNameTS.Contains("GluinoN_f10") && !ModelNameTS.Contains("f100") && ((int)(Mass)/100)%4==3 && TypeMode==0) ||
	 (ModelNameTS.Contains("Gluino_f50") && ((int)(Mass)/100)%4==3 && TypeMode==3) ||
	 (ModelNameTS.Contains("Gluino_f100") && ((int)(Mass)/100)%4==3 && TypeMode==3) ||
	 (ModelNameTS.Contains("Stop") && ((int)(Mass)/100)%3==2 && TypeMode==0) ||
	 (ModelNameTS.Contains("GMStau") && (Mass==126 || Mass==308 || Mass==494) && TypeMode==2) ||
	 (ModelNameTS.Contains("PPStau") && (Mass==126 || Mass==308 || Mass==494) && TypeMode==2) ||
         (ModelNameTS.Contains("DY")     && ModelNameTS.Contains("Q1")     && !ModelNameTS.Contains("o3") && ((int)(Mass)/100)%3==2 && TypeMode==2) ||
	 (ModelNameTS.Contains("DY")     && ((int)(Mass)/100)%3==2 && TypeMode==4) ||
         (ModelNameTS.Contains("DY")     && ((int)(Mass)/100)%1==0 && TypeMode==5)) {

	fprintf(pFile,"%s\\\\\n", ModelName.c_str());

      char massCut[255];  if(TypeMode<3){sprintf(massCut,"$>%.0f$",Infos8.MassCut);}else{sprintf(massCut," - ");}
      char Results7[255]; if(Infos7.Mass>0 && TypeMode!=3){sprintf(Results7, "%s & %s & %6.2f", toLatex(Infos7.XSec_Exp).c_str(), toLatex(Infos7.XSec_Obs).c_str(), Infos7.Eff);}else{sprintf(Results7, "    -     &   -      &   -     ");}
      char Results8[255]; if(Infos8.Mass>0){sprintf(Results8, "%s & %s & %6.2f", toLatex(Infos8.XSec_Exp).c_str(), toLatex(Infos8.XSec_Obs).c_str(), Infos8.Eff);}else{sprintf(Results8, "  -    &   -      &    -     &   -     ");}
      char ResultsC[255]; if(InfosC.Mass>0 && TypeMode!=3){sprintf(ResultsC, "%s & %s", toLatex(InfosC.XSec_Exp).c_str(), toLatex(InfosC.XSec_Obs).c_str());}
      else if(InfosC.Mass>0){sprintf(ResultsC, "%s & %s", toLatex(Infos8.XSec_Exp/Infos8.XSec_Th).c_str(), toLatex(Infos8.XSec_Obs/Infos8.XSec_Th).c_str());}
      else{sprintf(ResultsC, "   -     &    -    ");}

//    char Results7[255]; if(Infos7.Mass>0){sprintf(Results7, "%10s & %10s & %10s & %10s", toLatexRounded(Infos7.Eff).c_str(), toLatexRounded(Infos7.XSec_Th).c_str(),toLatexRounded(Infos7.XSec_Obs).c_str(), toLatexRounded(Infos7.XSec_Exp).c_str());}else{sprintf(Results7, "%10s & %10s & %10s & %10s", "", "", "", "");}
//    char Results8[255]; if(Infos8.Mass>0){sprintf(Results8, "%10s & %10s & %10s & %10s", toLatexRounded(Infos8.Eff).c_str(), toLatexRounded(Infos8.XSec_Th).c_str(),toLatexRounded(Infos8.XSec_Obs).c_str(), toLatexRounded(Infos8.XSec_Exp).c_str());}else{sprintf(Results8, "%10s & %10s & %10s & %10s", "", "", "", "");}
//    char ResultsC[255]; if(InfosC.Mass>0){sprintf(ResultsC, "%10s & %10s", toLatexRounded(InfosC.XSec_Obs).c_str(), toLatexRounded(InfosC.XSec_Exp).c_str());}else{sprintf(ResultsC, "%10s & %10s", "", "");}

      fprintf(pFile," %4.0f & %-7s & %s & %s & %s\\\\\n", Mass, massCut, Results7, Results8, ResultsC);
   }
   }
}


double GetSignalMeanHSCPPerEvent(string InputPattern, unsigned int CutIndex, double MinRange_, double MaxRange_){
   TFile* InputFile     = new TFile((InputPattern + "Histos.root").c_str());

   TH2D*  Mass                     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name          + "/Mass");
   TH2D*  MaxEventMass             = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name          + "/MaxEventMass");
   TH1D*  NTracksPassingSelection  = Mass->ProjectionY("NTracksPassingSelection",CutIndex+1,CutIndex+1);
   TH1D*  NEventsPassingSelection  = MaxEventMass->ProjectionY("NEventsPassingSelection",CutIndex+1,CutIndex+1);

   double NTracks       = NTracksPassingSelection->Integral(NTracksPassingSelection->GetXaxis()->FindBin(MinRange_), NTracksPassingSelection->GetXaxis()->FindBin(MaxRange_)); // all the way to the overflow
   double NEvents       = NEventsPassingSelection->Integral(NEventsPassingSelection->GetXaxis()->FindBin(MinRange_), NEventsPassingSelection->GetXaxis()->FindBin(MaxRange_));
   double toReturn      = (float)std::max(1.0,NTracks/ NEvents);
   delete Mass;
   delete MaxEventMass;
   delete NTracksPassingSelection;
   delete NEventsPassingSelection;

   delete InputFile;
   return toReturn;
}



void DrawModelLimitWithBand(string InputPattern){
   int TypeMode = TypeFromPattern(InputPattern);
   string prefix = "BUG";
   switch(TypeMode){
      case 0: prefix   = "Tk"; break;
      case 2: prefix   = "Mu"; break;
      case 3: prefix   = "Mo"; break;
      case 4: prefix   = "HQ"; break;
      case 5: prefix   = "LQ"; break;
   }

   double LInt = 0;
   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
      bool skip = false;
      bool isComb = EXCLUSIONDIR.find("COMB2016")!=string::npos;
      
      printf ("Model = %s\n", (modelVector[k]).c_str());

      if (modelVector[k] == "") continue;
      if (modelVector[k].find("16")!=string::npos && isComb) continue;
      if(TypeMode!=0 && isNeutral) continue;

      unsigned int N = modelMap[modelVector[k]].size();
      stAllInfo Infos;
      vector<double> Mass, XSecTh, XSecExp,XSecObs, XSecExpUp,XSecExpDown,XSecExp2Up,XSecExp2Down;
      for(unsigned int i=0;i<N;i++){
         string samplePath = InputPattern+""+SHAPESTRING+EXCLUSIONDIR+"/" + modelMap[modelVector[k]][i].Name +".txt";
	 if (isComb){
            samplePath = cleanSampleName (samplePath);
         }
         Infos = stAllInfo(samplePath);
         if (Infos.Mass < 100){
            printf ("Point at %s not found ...\n", samplePath.c_str());
	    skip = true;
	    N--;
	    continue;
         }
	 if (Mass.size() > 0 && Infos.Mass < Mass[Mass.size()-1]) break;
         std::cout << samplePath << std::endl;
         Mass        .push_back(Infos.Mass);
         XSecTh      .push_back(Infos.XSec_Th);
         XSecObs     .push_back(Infos.XSec_Obs);
         XSecExp     .push_back(Infos.XSec_Exp);
         XSecExpUp   .push_back(Infos.XSec_ExpUp);
         XSecExpDown .push_back(Infos.XSec_ExpDown);
         XSecExp2Up  .push_back(Infos.XSec_Exp2Up);
         XSecExp2Down.push_back(Infos.XSec_Exp2Down);
         LInt           =std::max(LInt, Infos.LInt);
      }
      if (skip) continue;
      TGraph* graphtheory  = new TGraph(Mass.size(),&(Mass[0]),&(XSecTh [0]));
      TGraph* graphobs     = new TGraph(Mass.size(),&(Mass[0]),&(XSecObs[0]));
      TGraph* graphexp     = new TGraph(Mass.size(),&(Mass[0]),&(XSecExp[0]));
      TCutG*  ExpErr       = GetErrorBand("ExpErr"      ,Mass.size(),&(Mass[0]), &(XSecExpDown[0]),&(XSecExpUp[0]),PlotMinScale,PlotMaxScale);
      TCutG*  Exp2SigmaErr = GetErrorBand("Exp2SigmaErr",Mass.size(),&(Mass[0]),&(XSecExp2Down[0]),&(XSecExp2Up[0]),PlotMinScale,PlotMaxScale);

      graphtheory->SetLineStyle(3);
      graphtheory->SetFillColor(kBlue);
      graphexp->SetLineStyle(4); 
      graphexp->SetLineColor(kRed);
      graphexp->SetMarkerStyle(); 
      graphexp->SetMarkerSize(0.); 
      Exp2SigmaErr->SetFillColor(kYellow);
      Exp2SigmaErr->SetLineColor(kWhite);
      ExpErr->SetFillColor(kGreen);
      ExpErr->SetLineColor(kWhite);
      graphobs->SetLineColor(kBlack);
      graphobs->SetLineWidth(2);
      graphobs->SetMarkerColor(kBlack);
      graphobs->SetMarkerStyle(23);

      TCanvas* c1 = new TCanvas("c1", "c1",600,600);
      TH1D* frame = new TH1D("frame", "frame", 1,50, 2650);
      frame->GetXaxis()->SetNdivisions(505);
      frame->SetTitle("");
      frame->SetStats(kFALSE);
      frame->GetXaxis()->SetTitle("Mass (GeV)");
      frame->GetYaxis()->SetTitle("95% CL limit on #sigma (pb)");
      frame->GetYaxis()->SetTitleOffset(1.40);
      frame->SetMaximum(PlotMaxScale);
      frame->SetMinimum(PlotMinScale);
      frame->GetYaxis()->SetRangeUser(2e-5, 1.5e1); // JOZZE EDIT
      frame->Draw("AXIS");


      Exp2SigmaErr->Draw("F");
      ExpErr      ->Draw("F");
      graphexp->Draw("LP");
      graphobs->Draw("LP");
      graphtheory->Draw("L");
      DrawPreliminary(LegendFromType(InputPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));

      TLegend* LEG = new TLegend(0.60,0.82-8*0.043,0.93,0.82);
      //TLegend* LEG = new TLegend(0.40,0.65,0.8,0.90);
      string headerstr = "95% CL Limits (";
      headerstr += LegendFromType(InputPattern) + string(")");
      LEG->SetHeader(headerstr.c_str());
      LEG->SetFillColor(0); 
      LEG->SetBorderSize(0);
      fprintf (stderr, "k/kMax = (%u / %lu)\tN=%lu out of %lu\n", k+1, modelVector.size(), Mass.size(), modelMap[modelVector[k]].size());
      LEG->AddEntry(graphtheory,  modelMap[modelVector[k]][0].ModelLegend().c_str() ,"L");
      LEG->AddEntry(graphexp    , "Expected"             ,"L");
      LEG->AddEntry(ExpErr      , "Expected #pm 1#sigma" ,"F");
      LEG->AddEntry(Exp2SigmaErr, "Expected #pm 2#sigma ","F");
      LEG->AddEntry(graphobs    , "Observed"             ,"LP");
      LEG->Draw();
      c1->SetLogy(true);

      SaveCanvas(c1,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", string(prefix+ modelVector[k] + "ExclusionLog"));
      delete frame;
      delete c1;
   }
}

// This code make the Expected Limit error band divided by expected limit plot for all signal models
// I don't like much this function... I started to rewrite it, but more work is still needed to improve it.
// I don't think two loops are needed, neither all these arrays...
void DrawRatioBands(string InputPattern)
{
   int TypeMode = TypeFromPattern(InputPattern);
   string prefix = "BUG";
   switch(TypeMode){
      case 0: prefix   = "Tk"; break;
      case 2: prefix   = "Mu"; break;
      case 3: prefix   = "Mo"; break;
      case 4: prefix   = "HQ"; break;
      case 5: prefix   = "LQ"; break;
   }

   TCanvas* c2            = new TCanvas("c2", "c2",600,800);
   TGraph** graphAtheory  = new TGraph*[modelVector.size()];
   TGraph** graphAobs     = new TGraph*[modelVector.size()];
   TGraph** graphAexp     = new TGraph*[modelVector.size()];
   TCutG**  ExpAErr       = new TCutG* [modelVector.size()];
   TCutG**  Exp2SigmaAErr = new TCutG* [modelVector.size()];
   TPad** padA            = new TPad*  [modelVector.size()];
   double step, top;
   double LInt = 0;

   top= 1.0/(modelVector.size()+2);
   step=(1.0-2.*top)/(modelVector.size());

   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
      if(TypeMode!=0 && isNeutral) continue;
      TPad* pad;
      if(k<(modelVector.size()-1)){
         pad = new TPad(Form("pad%i",k),Form("ExpErr%i",k),0.1,1-top-(k+1)*step,0.9,1-top-step*k);//lower left x, y, topright x, y
         pad->SetBottomMargin(0.);
      }else {
         pad = new TPad(Form("pad%i",k),Form("ExpErr%i",k),0.1,0.0,0.9,1-top-step*(k));//lower left x, y, topright x, y
         pad->SetBottomMargin(top/(step+top));
      }
      pad->SetLeftMargin(0.1);
      pad->SetRightMargin(0.);
      pad->SetTopMargin(0.);
      padA[k] = pad;  
      padA[k]->Draw();
   }

   for(unsigned int k=0; k<modelVector.size(); k++){
      bool isNeutral = false;if(modelVector[k].find("GluinoN")!=string::npos || modelVector[k].find("StopN")!=string::npos)isNeutral = true;
      if(TypeMode>0 && isNeutral) continue;

      TMultiGraph* MG = new TMultiGraph();
      unsigned int N = modelMap[modelVector[k]].size();
      stAllInfo Infos;double Mass[N], XSecTh[N], XSecExp[N],XSecObs[N], XSecExpUp[N],XSecExpDown[N],XSecExp2Up[N],XSecExp2Down[N];
      for(unsigned int i=0;i<N;i++){
         Infos = stAllInfo(InputPattern+""+SHAPESTRING+EXCLUSIONDIR+"/" + modelMap[modelVector[k]][i].Name +".txt");
         Mass        [i]=Infos.Mass;
         XSecTh      [i]=Infos.XSec_Th;
         XSecObs     [i]=Infos.XSec_Obs     /Infos.XSec_Exp;
         XSecExp     [i]=Infos.XSec_Exp     /Infos.XSec_Exp;
         XSecExpUp   [i]=Infos.XSec_ExpUp   /Infos.XSec_Exp;
         XSecExpDown [i]=Infos.XSec_ExpDown /Infos.XSec_Exp;
         XSecExp2Up  [i]=Infos.XSec_Exp2Up  /Infos.XSec_Exp;
         XSecExp2Down[i]=Infos.XSec_Exp2Down/Infos.XSec_Exp;
         LInt           =std::max(LInt, Infos.LInt);
      }

      TGraph* graphtheory  = new TGraph(N,Mass,XSecTh);
      TGraph* graphobs     = new TGraph(N,Mass,XSecObs);
      TGraph* graphexp     = new TGraph(N,Mass,XSecExp);
      TCutG*  ExpErr       = GetErrorBand(Form("ExpErr%i",k)      ,N,Mass,XSecExpDown ,XSecExpUp,  PlotMinScale, PlotMaxScale);
      TCutG*  Exp2SigmaErr = GetErrorBand(Form("Exp2SigmaErr%i",k),N,Mass,XSecExp2Down,XSecExp2Up, PlotMinScale, PlotMaxScale);

      graphAtheory [k] = graphtheory;      
      graphAobs    [k] = graphobs;
      graphAexp    [k] = graphexp;
      ExpAErr      [k] = ExpErr;

      Exp2SigmaAErr[k] = Exp2SigmaErr;
      graphAtheory [k]->SetLineStyle(3);
      graphAexp    [k]->SetLineStyle(4); 
      graphAexp    [k]->SetLineColor(kRed);
      graphAexp    [k]->SetMarkerStyle(); 
      graphAexp    [k]->SetMarkerSize(0.); 
      Exp2SigmaAErr[k]->SetFillColor(kYellow);
      Exp2SigmaAErr[k]->SetLineColor(kWhite);
      ExpAErr      [k]->SetFillColor(kGreen);
      ExpAErr      [k]->SetLineColor(kWhite);
      graphAobs    [k]->SetLineColor(kBlack);
      graphAobs    [k]->SetLineWidth(2);
      graphAobs    [k]->SetMarkerColor(kBlack);
      graphAobs    [k]->SetMarkerStyle(23);
      padA[k]->cd();
   
      int masst[2] = {0,1250};
      int xsect[2] = {2, 1};
      TGraph* graph = new TGraph(2,masst,xsect); //fake graph to set xaxis right
      graph->SetMarkerSize(0.);
      MG->Add(graph      ,"P");
      MG->Add(graphAobs[k]      ,"LP");
      MG->Draw("A");

      if(k==0){
	 TLegend* LEG;
	 LEG = new TLegend(0.13,0.01,0.32,0.99);
         string headerstr;
         headerstr = LegendFromType(InputPattern);
         LEG->SetHeader(headerstr.c_str());
         LEG->SetFillStyle(0); 
         LEG->SetBorderSize(0);
         LEG->AddEntry(ExpAErr[0], "Expected #pm 1#sigma","F");
         LEG->SetMargin(0.1);
         LEG->Draw();
      }  

      if(k==1){
         TLegend* LEG;
         LEG = new TLegend(0.13,0.01,0.32,0.99);
	 string headerstr;
	 LEG->SetFillStyle(0);
	 LEG->SetBorderSize(0);
	 LEG->AddEntry(Exp2SigmaAErr[0], "Expected #pm 2#sigma","F");
	 LEG->AddEntry(graphAobs[0],"Observed" ,"LP");
	 LEG->SetMargin(0.1);
	 LEG->Draw();
      }

      Exp2SigmaAErr[k]->Draw("f");
      ExpAErr[k]  ->Draw("f");
      MG->Draw("same");
      MG->SetTitle("");
      if(k==modelVector.size()-1) {
         MG->GetXaxis()->SetTitle("Mass (GeV)");
         MG->GetXaxis()->SetTitleSize(0.2);
         MG->GetXaxis()->SetLabelSize(0.2);
      }

      TPaveText *pt;
      if(TypeMode==0) {
      if(k!=modelVector.size()-1) pt = new TPaveText(0.45, 0.6, 0.95, 0.87,"LBNDC");
      else pt = new TPaveText(0.45, 0.82, 0.95, 0.935,"LBNDC");
      }
      else {
	if(k!=modelVector.size()-1) pt = new TPaveText(0.55, 0.6, 0.95, 0.87,"LBNDC");
	else pt = new TPaveText(0.55, 0.82, 0.95, 0.935,"LBNDC");
      }

      pt->SetBorderSize(0);
      pt->SetLineWidth(0);
      pt->SetFillStyle(kWhite);
      TText *text = pt->AddText(modelMap[modelVector[k]][0].ModelLegend().c_str()); 
      text ->SetTextAlign(12);
      text ->SetTextSize(0.3);
      if(k==modelVector.size()-1) text ->SetTextSize(0.5*text ->GetTextSize());
      pt->Draw();
      
      MG->GetXaxis()->SetRangeUser(0,1250);    
      MG->GetXaxis()->SetNdivisions(506,"Z");

      MG->GetYaxis()->SetRangeUser(0.001,2.99);
      MG->GetYaxis()->SetNdivisions(303, "Z");
      MG->GetYaxis()->SetLabelSize(0.3);
      if(k==modelVector.size()-1){
	MG->GetYaxis()->SetLabelSize(0.15);
      }
   }
   c2->cd();
   DrawPreliminary(LegendFromType(InputPattern).c_str(), SQRTS, IntegratedLuminosityFromE(SQRTS));
   TPaveText *pt = new TPaveText(0.1, 0., 0.15, 0.7,"NDC");
   string tmp = "95% CL Limits (Relative to Expected Limit)";
   TText *text = pt->AddText(tmp.c_str()); 
   text ->SetTextAlign(12);
   text ->SetTextAngle(90);
   text ->SetTextSize(0.04);
   pt->SetBorderSize(0);
   pt->SetFillColor(0);
   pt->Draw();
   SaveCanvas(c2,"Results/"+SHAPESTRING+EXCLUSIONDIR+"/", string(prefix+"LimitsRatio"));
   delete c2;
}

//will run on all possible selection and try to identify which is the best one for this sample
void Optimize(string InputPattern, string Data, string signal, bool shape, bool cutFromFile, int* OptCutIndex){
   printf("Optimize selection for %s in %s\n",signal.c_str(), InputPattern.c_str());fflush(stdout);

   //get the typeMode from pattern
   TypeMode = TypeFromPattern(InputPattern); 

   if (TypeMode == 3)    RescaleError = 0.0; //Done in step 4
   if (TypeMode == 4)    RescaleError = 0.20;

   //Identify the signal sample
   GetSampleDefinition(samples);
   CurrentSampleIndex        = JobIdToIndex(signal,samples); 
   printf("CurrentSampleIndez = %d\n", CurrentSampleIndex);
   if(CurrentSampleIndex<0){  printf("There is no signal corresponding to the JobId Given\n");  return;  } 

   if      (Data.find("7TeV"    )!=string::npos){SQRTS=7.0;    } //IntegratedLuminosity = IntegratedLuminosityFromE(SQRTS);
   else if (Data.find("8TeV"    )!=string::npos){SQRTS=8.0;    } //IntegratedLuminosity = IntegratedLuminosityFromE(SQRTS);
   else if (Data.find("13TeV15" )!=string::npos){SQRTS=1315.0; } //IntegratedLuminosity = IntegratedLuminosityFromE(SQRTS);
   else if (Data.find("13TeV16G")!=string::npos){SQRTS=13167.0;} //IntegratedLuminosity = IntegratedLuminosityFromE(SQRTS);
   else if (Data.find("13TeV16" )!=string::npos){SQRTS=1316.0; } //IntegratedLuminosity = IntegratedLuminosityFromE(SQRTS);
   else if (Data.find("13TeV"   )!=string::npos){SQRTS=13.0;   } //IntegratedLuminosity = IntegratedLuminosityFromE(SQRTS);

   //For muon only don't run on neutral samples as near zero efficiency can make jobs take very long time
   if((signal.find("Gluino")!=string::npos || signal.find("Stop")!=string::npos) && signal.find("N")!=string::npos && TypeMode==3) return;

   //Load all input histograms
   TFile*InputFile     = new TFile((InputPattern + "Histos.root").c_str());
   TH1D* HCuts_Pt      = (TH1D*)GetObjectFromPath(InputFile, "HCuts_Pt");
   TH1D* HCuts_I       = (TH1D*)GetObjectFromPath(InputFile, "HCuts_I");
   TH1D* HCuts_TOF     = (TH1D*)GetObjectFromPath(InputFile, "HCuts_TOF");
   TH1D* H_Lumi        = (TH1D*)GetObjectFromPath(InputFile, Data+"/IntLumi");
   TH1D* H_A           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_A");
   TH1D* H_B           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_B");
   TH1D* H_C           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_C");
   TH1D* H_D           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_D");
   TH1D* H_E           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_E");
   TH1D* H_F           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_F");
   TH1D* H_G           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_G");
   TH1D* H_H           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_H");
   TH1D* H_P           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_P");
   TH1D* H_S           = (TH1D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/TOF");
   TH2D* MassData      = (TH2D*)GetObjectFromPath(InputFile, Data+"/Mass");
   TH2D* MassPred      = (TH2D*)GetObjectFromPath(InputFile, Data+"/Pred_Mass");
   TH2D* MassSign      = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass" );
   if(!MassSign){printf("The sample %s is not present in the root file, returns\n", signal.c_str());return;}
   TH2D* MassSignP     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystP");
   TH2D* MassSignI     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystI");
   TH2D* MassSignM     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystM");
   TH2D* MassSignHUp   = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystHUp");
   TH2D* MassSignHDown = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystHDown");
   TH2D* MassSignT     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystT");
   TH2D* MassSignPU    = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystPU" );
   TH1D* TotalE        = (TH1D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/TotalE" );
   TH1D* TotalEPU      = (TH1D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/TotalEPU" );

   if(TypeMode==3) {
     //Need to add in systematic uncertainty of NPred, can't do it later because it comes from two different sources that have different uncertainties
     TH1D* H_P_Coll           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_P_Coll");
     TH1D* H_P_Cosmic           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_P_Cosmic");
     for(int i=0; i<H_P->GetNbinsX()+2; i++) {
       H_P->SetBinError(i, sqrt(H_P->GetBinError(i)*H_P->GetBinError(i) + H_P_Coll->GetBinContent(i)*0.2*H_P_Coll->GetBinContent(i)*0.2 + H_P_Cosmic->GetBinContent(i)*0.8*H_P_Cosmic->GetBinContent(i)*0.8));
     }
   }

   //If Take the cuts From File --> Load the actual cut index
   int OptimCutIndex = -1;  //int OptimMassWindow;
   if(cutFromFile && OptimCutIndex < 0){ // if less than zero, read from the file, otherwise don't
      FILE* pFile = fopen("Analysis_Cuts.txt","r");
      if(!pFile){printf("Can't open %s\n","Analysis_Cuts.txt"); return;}

      while(true){
         char line[4096];  string Name_;  int TypeMode_; double cutPt_; double cutI_; double cutTOF_; int massWindow_;
         if(!fgets(line, 4096, pFile))break;
         char* pch=strtok(line,","); int Arg=0; string tmp;
         while (pch!=NULL){
            switch(Arg){
               case  0: tmp = pch;  Name_     = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
               case  1: sscanf(pch, "%d",  &TypeMode_); break;
               case  2: sscanf(pch, "%lf", &cutPt_ ); break;
               case  3: sscanf(pch, "%lf", &cutI_  ); break;
               case  4: sscanf(pch, "%lf", &cutTOF_); break;
               case  5: sscanf(pch, "%d",  &massWindow_);break;
               default:break;
            }
            pch=strtok(NULL,",");Arg++;
         }
         //printf("%s %i %f %f %f %i\n",Name_.c_str(), TypeMode_, cutPt_, cutI_, cutTOF_, massWindow_);
         if(TypeMode_!=TypeMode)continue; //Not reading the cut line for the right TypeMode 

         string signalNameWithoutEnergy = signal;
         signalNameWithoutEnergy = ReplacePartOfString (signalNameWithoutEnergy, "_7TeV"    , "");
         signalNameWithoutEnergy = ReplacePartOfString (signalNameWithoutEnergy, "_8TeV"    , "");
         signalNameWithoutEnergy = ReplacePartOfString (signalNameWithoutEnergy, "_13TeV15" , "");
         signalNameWithoutEnergy = ReplacePartOfString (signalNameWithoutEnergy, "_13TeV16G", "");
         signalNameWithoutEnergy = ReplacePartOfString (signalNameWithoutEnergy, "_13TeV16" , "");
         signalNameWithoutEnergy = ReplacePartOfString (signalNameWithoutEnergy, "_13TeV"   , "");

//         if(signalNameWithoutEnergy.find(str7TeV)!=string::npos)signalNameWithoutEnergy.erase(signalNameWithoutEnergy.find(str7TeV), str7TeV.length());
//         if(signalNameWithoutEnergy.find(str8TeV)!=string::npos)signalNameWithoutEnergy.erase(signalNameWithoutEnergy.find(str8TeV), string(str8TeV).length()); 
//         if(signalNameWithoutEnergy.find(str8TeV)!=string::npos)signalNameWithoutEnergy.erase(signalNameWithoutEnergy.find(str8TeV), string(str8TeV).length()); 
//         if(signalNameWithoutEnergy.find(str8TeV)!=string::npos)signalNameWithoutEnergy.erase(signalNameWithoutEnergy.find(str8TeV), string(str8TeV).length()); 
//         if(signalNameWithoutEnergy.find(str13TeV)!=string::npos)signalNameWithoutEnergy.erase(signalNameWithoutEnergy.find(str13TeV), string(str13TeV).length()); 


         //printf("%s vs %s\n",Name_.c_str(), signalNameWithoutEnergy.c_str());
         if(Name_!=signalNameWithoutEnergy    )continue; //Not reading the cut line for the right signal sample

         //We are looking at the right sample --> Now loop over all cuts and identify the cut index of the optimal cut
         double MinDistance = 10000;
         for(int CutIndex=0;CutIndex<HCuts_Pt->GetNbinsX();CutIndex++){
            double cutDistance = fabs(cutPt_ - HCuts_Pt ->GetBinContent(CutIndex+1)) + fabs(cutI_ - HCuts_I ->GetBinContent(CutIndex+1)) + fabs(cutTOF_ - HCuts_TOF ->GetBinContent(CutIndex+1));
            if(cutDistance<MinDistance){MinDistance=cutDistance; OptimCutIndex=CutIndex;  }//OptimMassWindow=massWindow_;
         }
         printf("Closest cut index to the cuts provided: %i\n",OptimCutIndex);
         break; 
      }
      fclose(pFile);
      if(OptimCutIndex<0){printf("DID NOT FIND THE CUT TO USE FOR THIS SAMPLE %s\n",signal.c_str());return;}
   }

   //normalise the signal samples to XSection * IntLuminosity
   double LInt  = H_Lumi->GetBinContent(1); // FIXME JOZE
   if (signal.find("13TeV16G")!=string::npos) LInt = IntegratedLuminosity13TeV16G; // FIXME quick and dirty patch, but shouldn't be needed if you run from step 1 - 5
   else if (signal.find("13TeV16")!=string::npos) LInt = IntegratedLuminosity13TeV16PreG;
//   double LInt  = IntegratedLuminosity13TeV16 - IntegratedLuminosity13TeV16PostHIP;
//   LInt = Data.find("13TeV16")!=string::npos?IntegratedLuminosity13TeV16:IntegratedLuminosity13TeV15; // from before, but a neat trick
   double norm  = samples[CurrentSampleIndex].XSec*LInt/TotalE  ->Integral(); //normalize the samples to the actual lumi used for limits
   double normPU= samples[CurrentSampleIndex].XSec*LInt/(TotalEPU->Integral()>0?TotalEPU->Integral():TotalE->Integral());

   MassSign      ->Scale(norm);
   MassSignP     ->Scale(norm);
   MassSignI     ->Scale(norm);
   MassSignM     ->Scale(norm);
   MassSignHUp   ->Scale(norm);
   MassSignHDown ->Scale(norm);
   MassSignT     ->Scale(norm);
   MassSignPU    ->Scale(normPU);

   //prepare output directory and log file
   string outpath = InputPattern + "/"+SHAPESTRING+EXCLUSIONDIR+"/";
   MakeDirectories(outpath);
   FILE* pFile=NULL;
   if(OptimCutIndex<0){ //create the info file only if we need to optimize the cuts
      pFile = fopen((outpath+"/"+signal+".info").c_str(),"w");
      if(!pFile)printf("Can't open file : %s\n",(outpath+"/"+signal+".info").c_str());
   }

   //Compute mass range for the cut&count search
   stAllInfo result;
   stAllInfo toReturn;
   double Mean=-1,Width=-1;
   if(!shape && TypeMode<=2){
      TH1D* tmpMassSignProj = MassSign->ProjectionY("MassSignProj0",1,1);
      Mean                  = tmpMassSignProj->GetMean();
      Width                 = tmpMassSignProj->GetRMS();
      MinRange              = std::max(0.0, Mean-2*Width);
      if (OptimCutIndex>=0) { // if we have an optimal cut index
         TH1D* MassPredTMP = ((TH2D*)MassPred)->ProjectionY("MassPredProj", OptimCutIndex, OptimCutIndex); // FIXME
         result.NPred = MassPredTMP->Integral(MassPredTMP->GetXaxis()->FindBin(MinRange), MassPredTMP->GetXaxis()->FindBin(MaxRange)); // FIXME make sure that the mass cut has NPred > 1E-2
         for (; result.NPred < 1e-3; MinRange -= 10)
            result.NPred = MassPredTMP->Integral(MassPredTMP->GetXaxis()->FindBin(MinRange), MassPredTMP->GetXaxis()->FindBin(MaxRange));
         MinRange = tmpMassSignProj->GetXaxis()->GetBinLowEdge(tmpMassSignProj->GetXaxis()->FindBin(MinRange)); //Round to a bin value to avoid counting prpoblem due to the binning. 
         delete tmpMassSignProj; delete MassPredTMP;
      }
   }else{
      MinRange = 0;
   }
   printf ("MassCut = %.0lf, NPred = % .2e\n", MinRange, result.NPred);

   //loop on all possible selections and determine which is the best one
   for(int CutIndex=0;CutIndex<MassData->GetNbinsX();CutIndex++){
      //if(CutIndex>25)break; //for debugging purposes

      //if we don't want to optimize but take instead the cuts from a file, we can skip all other cuts
      if(OptimCutIndex>=0 && CutIndex!=OptimCutIndex)continue;

      //make sure the pT cut is high enough to get some statistic for both ABCD and mass shape
      if(HCuts_Pt ->GetBinContent(CutIndex+1) < 64 && TypeMode!=4){printf("Skip cut=%i because of too lose pT cut\n", CutIndex); continue; }

      //make sure we have a reliable prediction of the number of events
      if(OptimCutIndex<0 && H_E->GetBinContent(CutIndex+1) >0 &&(H_A->GetBinContent(CutIndex+1)<25 ||  H_F->GetBinContent(CutIndex+1)<25 || H_G->GetBinContent(CutIndex+1)<25)){printf("Skip cut=%i because of unreliable ABCD prediction\n", CutIndex); continue;}  //Skip events where Prediction (AFG/EE) is not reliable
      if(OptimCutIndex<0 && H_E->GetBinContent(CutIndex+1)==0 && H_A->GetBinContent(CutIndex+1)<=0 && (H_C->GetBinContent(CutIndex+1)<25 || H_B->GetBinContent(CutIndex+1)<25)){printf("Skip cut=%i because of unreliable ABCD prediction\n", CutIndex); continue;}  //Skip events where Prediction (CB/A) is not reliable
      if(OptimCutIndex<0 && H_F->GetBinContent(CutIndex+1) >0 && H_A->GetBinContent(CutIndex+1)<=0 && (H_B->GetBinContent(CutIndex+1)<25 || H_H->GetBinContent(CutIndex+1)<25)){printf("Skip cut=%i because of unreliable ABCD prediction\n", CutIndex); continue;}  //Skip events where Prediction (CB/A) is not reliable
      if(OptimCutIndex<0 && H_G->GetBinContent(CutIndex+1) >0 && H_F->GetBinContent(CutIndex+1)==0 && (H_C->GetBinContent(CutIndex+1)<25 || H_H->GetBinContent(CutIndex+1)<25)){printf("Skip cut=%i because of unreliable ABCD prediction\n", CutIndex); continue;}  //Skip events where Prediction (CH/G) is not reliable
   
      //make sure we have a reliable prediction of the shape 
      if(TypeMode<=2){
         double N_P = H_P->GetBinContent(CutIndex+1);       
         if(H_E->GetBinContent(CutIndex+1) >0 && (H_A->GetBinContent(CutIndex+1)<0.25*N_P || H_F->GetBinContent(CutIndex+1)<0.25*N_P || H_G->GetBinContent(CutIndex+1)<0.25*N_P)){printf("Skip cut=%i because of unreliable mass prediction\n", CutIndex); continue;}  //Skip events where Mass Prediction is not reliable
         if(H_E->GetBinContent(CutIndex+1)==0 && (H_C->GetBinContent(CutIndex+1)<0.25*N_P || H_B->GetBinContent(CutIndex+1)<0.25*N_P)){printf("Skip cut=%i because of unreliable mass prediction\n", CutIndex); continue;}  //Skip events where Mass Prediction is not reliable
      }

      //prepare outputs result structure
      result = toReturn;
      result.MassMean  = Mean;
      result.MassSigma = Width;
//      if (correctMassCut && TypeMode<=2 && !shape){
//         TH1D* tmpMassSignProj = MassSign->ProjectionY("MassSignProj0",1,1);
//         MinRange = tmpMassSignProj->GetXaxis()->GetBinLowEdge(tmpMassSignProj->GetXaxis()->FindBin(samples[JobIdToIndex(signal,samples)].Mass*0.625));
//      }
      result.MassCut   = TypeMode<=2?MinRange:0;
      result.Mass      = samples[JobIdToIndex(signal,samples)].Mass;
//      if (signal.find("Q2")!=string::npos && SQRTS==1316.0 && result.Mass > 1000/* && TypeMode==0*/){
//         if      (result.Mass == 1400) result.MassCut = 300;
//         else if (result.Mass == 1800) result.MassCut = 190;
//         else if (result.Mass == 2200) result.MassCut =  70;
//         else if (result.Mass == 2600) result.MassCut =  50;
//         result.MassCut = result.Mass*0.5;
//	 MinRange = result.MassCut;
//      }
      result.XSec_Th   = samples[JobIdToIndex(signal,samples)].XSec;
      result.XSec_Err  = samples[JobIdToIndex(signal,samples)].XSec * 0.15;
      toReturn = result;
      result.Index     = CutIndex;
      result.WP_Pt     = HCuts_Pt ->GetBinContent(CutIndex+1);
      result.WP_I      = HCuts_I  ->GetBinContent(CutIndex+1);
      result.WP_TOF    = HCuts_TOF->GetBinContent(CutIndex+1);
      result.LInt      = LInt;

      //best significance --> is actually best reach
      if(TypeMode<=2){if(!runCombine(true, false, true, InputPattern, signal, CutIndex, shape, true, result, MassData, MassPred, MassSign, MassSignP, MassSignI, MassSignM, MassSignHUp, MassSignHDown, MassSignT, MassSignPU)){printf("runCombine did not converge\n"); continue;}
      }else          {if(!runCombine(true, false, true, InputPattern, signal, CutIndex, shape, true, result, H_D, H_P, MassSign, MassSignP, MassSignI, MassSignM, MassSignHUp, MassSignHDown, MassSignT, MassSignPU)){printf("runCombine did not converge\n"); continue;}
      }

      //report the result for this point in the log file
      if(pFile)fprintf(pFile  ,"%10s: Testing CutIndex=%4i (Pt>%6.2f I>%6.3f TOF>%6.3f) %3.0f<M<inf Ndata=%+6.2E NPred=%6.3E+-%6.3E SignalEff=%6.3f ExpLimit=%6.3E (%6.3E) Reach=%6.3E",signal.c_str(),CutIndex,HCuts_Pt ->GetBinContent(CutIndex+1), HCuts_I  ->GetBinContent(CutIndex+1), HCuts_TOF->GetBinContent(CutIndex+1), MinRange,result.NData,result.NPred, result.NPredErr,result.Eff,result.XSec_Exp, result.XSec_Obs, result.XSec_5Sigma);
      fprintf(stdout ,"%10s: Testing CutIndex=%4i (Pt>%6.2f I>%6.3f TOF>%6.3f) %3.0f<M<inf Ndata=%+6.2E NPred=%6.3E+-%6.3E SignalEff=%6.3f ExpLimit=%6.3E (%6.3E) Reach=%6.3E",signal.c_str(),CutIndex,HCuts_Pt ->GetBinContent(CutIndex+1), HCuts_I  ->GetBinContent(CutIndex+1), HCuts_TOF->GetBinContent(CutIndex+1), MinRange,result.NData,result.NPred, result.NPredErr,result.Eff,result.XSec_Exp, result.XSec_Obs, result.XSec_5Sigma);
      if(OptimCutIndex>=0 || (result.XSec_5Sigma>0 && result.XSec_5Sigma<toReturn.XSec_5Sigma)){
         toReturn=result;
         if(pFile)fprintf(pFile  ," BestSelection\n");
         fprintf(stdout ," BestSelection\n");
      }else{
         if(pFile)fprintf(pFile  ,"\n");
         fprintf(stdout ,"\n");
      }
 
      if(pFile)fflush(pFile);
      fflush(stdout);
   }//end of selection cut loop
   if(pFile)fclose(pFile);   
 
   //recompute the limit for the final point and save the output in the final directory (also save some plots for the shape based analysis)
   if(TypeMode<=2){runCombine(false, true, true, InputPattern, signal, toReturn.Index, shape, false, toReturn, MassData, MassPred, MassSign, MassSignP, MassSignI, MassSignM, MassSignHUp, MassSignHDown, MassSignT, MassSignPU);
   }else          {runCombine(false, true, true, InputPattern, signal, toReturn.Index, shape, false, toReturn, H_D, H_P, MassSign, MassSignP, MassSignI, MassSignM, MassSignHUp, MassSignHDown, MassSignT, MassSignPU);
   }
  
   //all done, save the result to file
   toReturn.Save(outpath+"/"+signal+".txt");
}

// produce the Higgs combine stat tool datacard
void makeDataCard(string outpath, string rootPath, string ChannelName, string SignalName, double Obs, double Pred, double PredRelErr, double Sign, double SignStat, double SignalUnc, bool Shape){

   double LumiUnc   = 1.0;
   if(SQRTS==7   ) LumiUnc=1.022;
   if(SQRTS==8   ) LumiUnc=1.044;
   if(SQRTS==13  ) LumiUnc=1.027;
   if(SQRTS==1315) LumiUnc=1.027;
   if(SQRTS==1316) LumiUnc=1.062; // noted from the experts, 6.2%

   if(isnan(float(PredRelErr)))PredRelErr= 1.2;


   FILE* pFile = fopen(outpath.c_str(), "w");
   fprintf(pFile, "imax 1\n");
   fprintf(pFile, "jmax *\n");
   fprintf(pFile, "kmax *\n");
   if(Shape){
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "shapes * * %s $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n",rootPath.c_str());
   }
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "bin %s\n",ChannelName.c_str());
   fprintf(pFile, "Observation %f\n",Obs);
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "bin      %s %s\n",ChannelName.c_str(), ChannelName.c_str());
   fprintf(pFile, "process  %s pred%i\n",SignalName.c_str(), TypeMode);
   fprintf(pFile, "process  0 1\n");
   fprintf(pFile, "rate    %f %f\n",Sign,std::max(1E-4, Pred) );  //if Pred<1E-4 we have troubles when merging datacards
   fprintf(pFile, "-------------------------------\n");
   fprintf(pFile, "%35s    %6s %5.3f     1.0  \n","Lumi" , "lnN", LumiUnc);
   fprintf(pFile, "%35s    %6s -         %5.3f\n",(ChannelName+"systP").c_str(), "lnN", PredRelErr);
   fprintf(pFile, "%35s    %6s %5.3f     -    \n",(ChannelName+"systS").c_str(), "lnN", SignalUnc);
   if(!Shape){
   fprintf(pFile, "%35s    %6s %5.3f     -    \n",(ChannelName+"statS").c_str(), "lnN", std::min(SignStat,2.0));
   }else{
   fprintf(pFile, "%35s    %6s 1.000     -    \n",(ChannelName+"statS").c_str(), "shapeN2");
   fprintf(pFile, "%35s    %6s -         1    \n",(ChannelName+"statP").c_str(), "shapeN2");
   fprintf(pFile, "%35s    %6s 1.000     -    \n","mom"  , "shapeN2");
   fprintf(pFile, "%35s    %6s 1.000     -    \n","ias"  , "shapeN2");
   fprintf(pFile, "%35s    %6s 1.000     -    \n","ih"   , "shapeN2");
   fprintf(pFile, "%35s    %6s 1.000     -    \n","hip"  , "shapeN2");
   fprintf(pFile, "%35s    %6s 1.000     -    \n","tof"  , "shapeN2");
   fprintf(pFile, "%35s    %6s 1.000     -    \n","pu"   , "shapeN2");
   }
   fclose(pFile);
}

//save histogram in root file (and it's statistical variation if it's not observed date histogram)
void saveHistoForLimit(TH1* histo, string Name, string Id){            
      histo   ->Write( Name                   .c_str());
      if(Name=="data_obs")return;
      
      TH1* statup   = (TH1*)histo->Clone((Name+"_stat"+Id+"Up").c_str());
      TH1* statdown = (TH1*)histo->Clone((Name+"_stat"+Id+"Down").c_str());       
      for(int ibin=1; ibin<=statup->GetXaxis()->GetNbins(); ibin++){
         statup  ->SetBinContent(ibin,statup  ->GetBinContent(ibin) + histo->GetBinError(ibin));
         statdown->SetBinContent(ibin,statdown->GetBinContent(ibin) - histo->GetBinError(ibin));
      }
      statup  ->Write((Name+"_stat"+Id+"Up"  ).c_str());
      statdown->Write((Name+"_stat"+Id+"Down").c_str());      

      delete statup;
      delete statdown;
}

//save histogram with variation in root file  (save it as symetrical up and down variation)
void saveVariationHistoForLimit(TH1* histo, TH1* vardown, string Name, string variationName){
      TH1* varup   = (TH1*)histo ->Clone((Name+"_"+variationName+"Up"  ).c_str());
      varup->Add(vardown,-1); //varup=x
      varup->Add(histo,1);
      varup  ->Write((Name+"_"+variationName+"Up"  ).c_str());
      vardown->Write((Name+"_"+variationName+"Down").c_str());
}

//function for debugging only, should be remove soon
//this function just compute the result of the shape based analysis using the optimal point from the cut&count analysis
void testShapeBasedAnalysis(string InputPattern, string signal){
   CurrentSampleIndex        = JobIdToIndex(signal, samples); if(CurrentSampleIndex<0){  printf("There is no signal corresponding to the JobId Given\n");  return;  }
   int s = CurrentSampleIndex;

   string outpath = InputPattern + "/"+SHAPESTRING+EXCLUSIONDIR+"/";
   MakeDirectories(outpath);

   //Get Optimal cut from cut&count optimization
   stAllInfo result =  stAllInfo(InputPattern+"/"+EXCLUSIONDIR+"/"+signal+".txt");

   //load all intput histograms
   TFile* InputFile  = new TFile((InputPattern+"/Histos.root").c_str());
   TH1D*  H_Lumi     = (TH1D*)GetObjectFromPath(InputFile, "Data7TeV/IntLumi");
   TH2D*  MassData   = (TH2D*)GetObjectFromPath(InputFile, "Data7TeV/Mass");
   TH2D*  MassPred   = (TH2D*)GetObjectFromPath(InputFile, "Data7TeV/Pred_Mass");
   TH2D*  MassSign   = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass");
   TH2D*  MassSignP  = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass_SystP");
   TH2D*  MassSignI  = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass_SystI");
   TH2D*  MassSignM  = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass_SystM");
   TH2D*  MassSignHUp   = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass_SystHUp");
   TH2D*  MassSignHDown = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass_SystHDown");
   TH2D*  MassSignT  = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass_SystT");
   TH2D*  MassSignPU = (TH2D*)GetObjectFromPath(InputFile, samples[s].Name+"/Mass_SystPU");
   TH1D*  TotalE     = (TH1D*)GetObjectFromPath(InputFile, samples[s].Name+"/TotalE" );
   TH1D*  TotalEPU   = (TH1D*)GetObjectFromPath(InputFile, samples[s].Name+"/TotalEPU" );

   //normalise the signal samples to XSection * IntLuminosity
   double LInt  = H_Lumi->GetBinContent(1); // FIXME JOZE
//   double LInt  = IntegratedLuminosity13TeV16 - IntegratedLuminosity13TeV16PostHIP;
   double norm  = samples[CurrentSampleIndex].XSec*LInt/TotalE  ->Integral(); //normalize the samples to the actual lumi used for limits
   double normPU= samples[CurrentSampleIndex].XSec*LInt/TotalEPU->Integral();
   MassSign      ->Scale(norm);
   MassSignP     ->Scale(norm);
   MassSignI     ->Scale(norm);
   MassSignM     ->Scale(norm);
   MassSignHUp   ->Scale(norm);
   MassSignHDown ->Scale(norm);
   MassSignT     ->Scale(norm);
   MassSignPU    ->Scale(normPU);

   bool Shape = true;

   //find range
   if(Shape){
      MinRange = 0;
   }else{
      MinRange = std::max(0.0, result.MassMean-2*result.MassSigma);
      MinRange = MassSign->GetYaxis()->GetBinLowEdge(MassSign->GetYaxis()->FindBin(MinRange)); //Round to a bin value to avoid counting prpoblem due to the binning. 
   }

   //compute shape based limits and save it's output
   runCombine(false, true, true, InputPattern, signal, result.Index, Shape, false, result, MassData, MassPred, MassSign, MassSignP, MassSignI, MassSignM, MassSignHUp, MassSignHDown, MassSignT, MassSignPU);

   //all done, save the results to file
   result.Save(InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/"+signal+".txt");
}

//compute the significance using a ProfileLikelihood (assuming datacard is already produced)
double computeSignificance(string datacard, bool expected, string& signal, string massStr, float Strength){
   double toReturn = -1;
   char strengthStr[255]; sprintf(strengthStr,"--expectSignal=%f",Strength);
   string CodeToExecute = "cd /tmp/;";
   if(expected)CodeToExecute += "combine -M ProfileLikelihood -n " + signal + " -m " + massStr + " --significance -t 100 " + strengthStr + " " + datacard + " &> shape_" + signal + ".log;";
   else        CodeToExecute += "combine -M ProfileLikelihood -n " + signal + " -m " + massStr + " --significance                            " + datacard + " &> shape_" + signal + ".log;";
   CodeToExecute += "cd $OLDPWD;";
   system(CodeToExecute.c_str());   

   char line[4096];
   FILE* sFile = fopen((string("/tmp/shape_")+signal + ".log").c_str(), "r");
   if(!sFile)std::cout<<"FILE NOT OPEN:"<< (string("/tmp/shape_")+signal + ".log").c_str() << endl;
   int LineIndex=0; int GarbageI; double GarbageD;     
   while(fgets(line, 4096, sFile)){LineIndex++;       
     if(!expected && LineIndex==3){sscanf(line,"Significance: %lf",&toReturn);     break;}
 //    if( expected && LineIndex==7){sscanf(line,"median expected limit: r < %lf @ 95%%CL (%i toyMC)",&toReturn,&GarbageI); break;}
     if( expected && LineIndex==6){sscanf(line,"mean   expected limit: r < %lf +/- %lf @ 95%%CL (%i toyMC)",&toReturn, &GarbageD, &GarbageI); break;}

     continue;
   }fclose(sFile);
   return toReturn;
}

//run the higgs combine stat tool using predicted mass shape distribution (possibly do shape based analysis and/or cut on mass) OR 1D histogram output from ABCD  (only do cut and count without mass cut)
bool runCombine(bool fastOptimization, bool getXsection, bool getSignificance, string& InputPattern, string& signal, unsigned int CutIndex, bool Shape, bool Temporary, stAllInfo& result, TH1* MassData, TH1* MassPred, TH1* MassSign, TH1* MassSignP, TH1* MassSignI, TH1* MassSignM, TH1* MassSignHUp, TH1* MassSignHDown, TH1* MassSignT, TH1* MassSignPU){
   TH1D *MassDataProj=NULL, *MassPredProj=NULL, *MassSignProj=NULL, *MassSignProjP=NULL, *MassSignProjI=NULL, *MassSignProjM=NULL, *MassSignProjHUp=NULL, *MassSignProjHDown=NULL, *MassSignProjT=NULL, *MassSignProjPU=NULL;
   double NData, NPredErr, NPred, NSign, NSignP, NSignI, NSignM, NSignHUp, NSignHDown, NSignT, NSignPU;
   double NSignErr, NSignPErr, NSignIErr, NSignMErr, NSignHErrUp, NSignHErrDown, NSignTErr, NSignPUErr;
   double signalsMeanHSCPPerEvent = GetSignalMeanHSCPPerEvent(InputPattern,CutIndex, MinRange, MaxRange);

   printf ("signalsMeanHSCPPerEvent = %lf\n", signalsMeanHSCPPerEvent);
   printf ("CurrentSampleIndex = %d, CutIndex = %d\n", CurrentSampleIndex, CutIndex);
   //IF 2D histograms --> we get all the information from there (and we can do shape based analysis AND/OR cut on mass)
   if(MassData->InheritsFrom("TH2")){
      //make the projection of all the 2D input histogram to get the shape for this single point
      MassDataProj       = ((TH2D*)MassData  )->ProjectionY("MassDataProj"  ,CutIndex+1,CutIndex+1);
      MassPredProj       = ((TH2D*)MassPred  )->ProjectionY("MassPredProj"  ,CutIndex+1,CutIndex+1);
      MassSignProj       = ((TH2D*)MassSign  )->ProjectionY("MassSignProj"  ,CutIndex+1,CutIndex+1);
      MassSignProjP      = ((TH2D*)MassSignP )->ProjectionY("MassSignProP"  ,CutIndex+1,CutIndex+1);
      MassSignProjI      = ((TH2D*)MassSignI )->ProjectionY("MassSignProI"  ,CutIndex+1,CutIndex+1);
      MassSignProjM      = ((TH2D*)MassSignM )->ProjectionY("MassSignProM"  ,CutIndex+1,CutIndex+1);
      MassSignProjHUp    = ((TH2D*)MassSignHUp )->ProjectionY("MassSignProH"  ,CutIndex+1,CutIndex+1);
      MassSignProjHDown  = ((TH2D*)MassSignHDown )->ProjectionY("MassSignProHDown"  ,CutIndex+1,CutIndex+1);
      MassSignProjT      = ((TH2D*)MassSignT )->ProjectionY("MassSignProT"  ,CutIndex+1,CutIndex+1);
      MassSignProjPU     = ((TH2D*)MassSignPU)->ProjectionY("MassSignProPU" ,CutIndex+1,CutIndex+1);

      //count events in the allowed range (infinite for shape based and restricted for cut&count)
      NData       = MassDataProj->Integral(MassDataProj->GetXaxis()->FindBin(MinRange), MassDataProj->GetXaxis()->FindBin(MaxRange));
      NPred       = MassPredProj->Integral(MassPredProj->GetXaxis()->FindBin(MinRange), MassPredProj->GetXaxis()->FindBin(MaxRange));
      NPredErr    = pow(NPred*RescaleError,2);
      for(int i=MassPredProj->GetXaxis()->FindBin(MinRange); i<=MassPredProj->GetXaxis()->FindBin(MaxRange) ;i++){NPredErr+=pow(MassPredProj->GetBinError(i),2);}NPredErr=sqrt(NPredErr);
      NSign       = (MassSignProj  ->IntegralAndError(MassSignProj  ->GetXaxis()->FindBin(MinRange), MassSignProj  ->GetXaxis()->FindBin(MaxRange), NSignErr  )) / signalsMeanHSCPPerEvent;  NSignErr  /= signalsMeanHSCPPerEvent;
      NSignP      = (MassSignProjP ->IntegralAndError(MassSignProjP ->GetXaxis()->FindBin(MinRange), MassSignProjP ->GetXaxis()->FindBin(MaxRange), NSignPErr )) / signalsMeanHSCPPerEvent;  NSignPErr /= signalsMeanHSCPPerEvent;
      NSignI      = (MassSignProjI ->IntegralAndError(MassSignProjI ->GetXaxis()->FindBin(MinRange), MassSignProjI ->GetXaxis()->FindBin(MaxRange), NSignIErr )) / signalsMeanHSCPPerEvent;  NSignIErr /= signalsMeanHSCPPerEvent;
      NSignM      = (MassSignProjM ->IntegralAndError(MassSignProjM ->GetXaxis()->FindBin(MinRange), MassSignProjM ->GetXaxis()->FindBin(MaxRange), NSignMErr )) / signalsMeanHSCPPerEvent;  NSignMErr /= signalsMeanHSCPPerEvent;
      NSignHUp    = (MassSignProjHUp   ->IntegralAndError(MassSignProjHUp   ->GetXaxis()->FindBin(MinRange), MassSignProjHUp   ->GetXaxis()->FindBin(MaxRange), NSignHErrUp ))   / signalsMeanHSCPPerEvent;  NSignHErrUp   /= signalsMeanHSCPPerEvent;
      NSignHDown  = (MassSignProjHDown ->IntegralAndError(MassSignProjHDown ->GetXaxis()->FindBin(MinRange), MassSignProjHDown ->GetXaxis()->FindBin(MaxRange), NSignHErrDown )) / signalsMeanHSCPPerEvent;  NSignHErrDown /= signalsMeanHSCPPerEvent;
      NSignT      = (MassSignProjT ->IntegralAndError(MassSignProjT ->GetXaxis()->FindBin(MinRange), MassSignProjT ->GetXaxis()->FindBin(MaxRange), NSignTErr )) / signalsMeanHSCPPerEvent;  NSignTErr /= signalsMeanHSCPPerEvent;
      NSignPU     = (MassSignProjPU->IntegralAndError(MassSignProjPU->GetXaxis()->FindBin(MinRange), MassSignProjPU->GetXaxis()->FindBin(MaxRange), NSignPUErr)) / signalsMeanHSCPPerEvent;  NSignPUErr/= signalsMeanHSCPPerEvent;

   //IF 1D histograms --> we get all the information from the ABCD method output 
   }else{
      Shape=false; //can not do shape based if we don't get the shapes
      NData       = MassData  ->GetBinContent(CutIndex+1);
      NPredErr    = MassPred  ->GetBinError  (CutIndex+1);
      NPred       = MassPred  ->GetBinContent(CutIndex+1);
      //NSign       = MassSign  ->GetBinContent(CutIndex+1) / signalsMeanHSCPPerEvent;
      //NSignErr    = MassSign  ->GetBinError  (CutIndex+1) / signalsMeanHSCPPerEvent;

      MassSignProj       = ((TH2D*)MassSign )->ProjectionY("MassSignPro"  ,CutIndex+1,CutIndex+1);
      MassSignProjP      = ((TH2D*)MassSignP )->ProjectionY("MassSignProP"  ,CutIndex+1,CutIndex+1);
      MassSignProjI      = ((TH2D*)MassSignI )->ProjectionY("MassSignProI"  ,CutIndex+1,CutIndex+1);
      MassSignProjM      = ((TH2D*)MassSignM )->ProjectionY("MassSignProM"  ,CutIndex+1,CutIndex+1);
      MassSignProjHUp    = ((TH2D*)MassSignHUp  )->ProjectionY("MassSignProHUp"    ,CutIndex+1,CutIndex+1);
      MassSignProjHDown  = ((TH2D*)MassSignHDown)->ProjectionY("MassSignProHDown"  ,CutIndex+1,CutIndex+1);
      MassSignProjT      = ((TH2D*)MassSignT )->ProjectionY("MassSignProT"  ,CutIndex+1,CutIndex+1);
      MassSignProjPU     = ((TH2D*)MassSignPU)->ProjectionY("MassSignProPU" ,CutIndex+1,CutIndex+1);

      NSign       = MassSignProj  ->IntegralAndError(0, MassSignProj  ->GetNbinsX()+1, NSignErr )  / signalsMeanHSCPPerEvent;  NSignErr  /= signalsMeanHSCPPerEvent;
      NSignP      = MassSignProjP ->IntegralAndError(0, MassSignProjP ->GetNbinsX()+1, NSignPErr)  / signalsMeanHSCPPerEvent;  NSignPErr /= signalsMeanHSCPPerEvent;
      NSignI      = MassSignProjI ->IntegralAndError(0, MassSignProjI ->GetNbinsX()+1, NSignIErr)  / signalsMeanHSCPPerEvent;  NSignIErr /= signalsMeanHSCPPerEvent;
      NSignM      = MassSignProjM ->IntegralAndError(0, MassSignProjM ->GetNbinsX()+1, NSignMErr)  / signalsMeanHSCPPerEvent;  NSignMErr /= signalsMeanHSCPPerEvent;
      NSignHUp    = MassSignProjHUp   ->IntegralAndError(0, MassSignProjHUp   ->GetNbinsX()+1, NSignHErrUp  )  / signalsMeanHSCPPerEvent;  NSignHErrUp   /= signalsMeanHSCPPerEvent;
      NSignHDown  = MassSignProjHDown ->IntegralAndError(0, MassSignProjHDown ->GetNbinsX()+1, NSignHErrDown)  / signalsMeanHSCPPerEvent;  NSignHErrDown /= signalsMeanHSCPPerEvent;
      NSignT      = MassSignProjT ->IntegralAndError(0, MassSignProjT ->GetNbinsX()+1, NSignTErr)  / signalsMeanHSCPPerEvent;  NSignTErr /= signalsMeanHSCPPerEvent;
      NSignPU     = MassSignProjPU->IntegralAndError(0, MassSignProjPU->GetNbinsX()+1, NSignPUErr) / signalsMeanHSCPPerEvent;  NSignPUErr/= signalsMeanHSCPPerEvent;

      NPredErr = sqrt(pow((NPred* RescaleError), 2) + pow(NPredErr,2));      // incorporate background uncertainty
   }

   //skip pathological selection point
   if(isnan((float)NPred)){printf("SkipThisPoint --> Nan Background prediction\n"); return false;}
//   if(NPred<=0){printf("SkipThisPoint --> NoBackgroundPrediction\n");return false;} //Is <=0 only when prediction failed or is not meaningful (i.e. WP=(0,0,0) ) //FIXME
   if(!Shape && NPred>1000){printf("SkipThisPoint --> NPred is way too big\n"); return false;}  //When NPred is too big, expected limits just take an infinite time! 

   //compute all efficiencies (not really needed anymore, but it's nice to look at these numbers afterward)
   double Eff         = NSign   / (result.XSec_Th*result.LInt);
   double EffP        = NSignP  / (result.XSec_Th*result.LInt);
   double EffI        = NSignI  / (result.XSec_Th*result.LInt);
   double EffM        = NSignM  / (result.XSec_Th*result.LInt);
   double EffHUp      = NSignHUp/ (result.XSec_Th*result.LInt);
   double EffHDown    = NSignHDown  / (result.XSec_Th*result.LInt);
   double EffT        = NSignT  / (result.XSec_Th*result.LInt);
   double EffPU       = NSignPU / (result.XSec_Th*result.LInt);

   double EffErr      = NSignErr   / (result.XSec_Th*result.LInt);
   double EffPErr     = NSignPErr  / (result.XSec_Th*result.LInt);
   double EffIErr     = NSignIErr  / (result.XSec_Th*result.LInt);
   double EffMErr     = NSignMErr  / (result.XSec_Th*result.LInt);
   double EffHErrUp   = NSignHErrUp  / (result.XSec_Th*result.LInt);
   double EffHErrDown = NSignHErrDown  / (result.XSec_Th*result.LInt);
   double EffTErr     = NSignTErr  / (result.XSec_Th*result.LInt);
   double EffPUErr    = NSignPUErr / (result.XSec_Th*result.LInt);


   if(Eff==0){printf("SkipThisPoint --> Signal acceptance is null\n"); return false;}
//   if(Eff<=1E-5)return false; // if Eff<0.001% -> limit will hardly converge and we are probably not interested by this point anyway

   //save these info to the result structure
   result.Eff       = Eff;
   result.Eff_SYSTP = EffP;
   result.Eff_SYSTI = EffI;
   result.Eff_SYSTM = EffM;
   result.Eff_SYSTHUp   = EffHUp;
   result.Eff_SYSTHDown = EffHDown;
   result.Eff_SYSTT = EffT;
   result.Eff_SYSTPU= EffPU;
   result.EffE       = EffErr;
   result.EffE_SYSTP = EffPErr;
   result.EffE_SYSTI = EffIErr;
   result.EffE_SYSTM = EffMErr;
   result.EffE_SYSTHUp   = EffHErrUp;
   result.EffE_SYSTHDown = EffHErrDown;
   result.EffE_SYSTT = EffTErr;
   result.EffE_SYSTPU= EffPUErr;
   result.NData     = NData;
   result.NPred     = NPred;
   result.NPredErr  = NPredErr;
   result.NSign     = NSign;
   printf ("NSign = %lf = %lf/(%.2e*%.4lf)\n", NSign/(result.XSec_Th*result.LInt), NSign, result.XSec_Th, result.LInt);
   NSign /= (result.XSec_Th*result.LInt); //normalize signal to 1pb
//   NPred /= (result.XSec_Th*result.LInt); //normalize signal to 1pb
   double SignalScaleFactor = 1.0;
   for(unsigned int i=0;i<20 && NSign<1e-1; i++){SignalScaleFactor*=10.0;  NSign*=10.0;}  
   if(NPred<0.001) {NPred=0.001; NPredErr=NPred;}


   //no way that this point is optimal
   bool pointMayBeOptimal = (fastOptimization && !getXsection && getSignificance && ((NPred-3*NPredErr)<=result.NPred || Eff>=result.Eff));


   //for shape based analysis we need to save all histograms into a root file
   char CutIndexStr[255];sprintf(CutIndexStr, "SQRTS%02.0fCut%03.0f",SQRTS, result.Index);
   if(Shape){
      //prepare the histograms and variation
      //scale to 10fb xsection and to observed events instead of observed tracks
      MassSignProj  ->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));
      MassSignProjP ->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));
      MassSignProjI ->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));
      MassSignProjM ->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));
      MassSignProjHUp   ->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));
      MassSignProjHDown ->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));
      MassSignProjT ->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));
      MassSignProjPU->Scale(1.0/(result.XSec_Th*signalsMeanHSCPPerEvent*result.LInt));

      //Rebin --> keep CPU time reasonable and error small
      MassDataProj  ->Rebin(2);
      MassPredProj  ->Rebin(2);
      MassSignProj  ->Rebin(2);
      MassSignProjP ->Rebin(2);
      MassSignProjI ->Rebin(2);
      MassSignProjM ->Rebin(2);
      MassSignProjHUp   ->Rebin(2);
      MassSignProjHDown ->Rebin(2);
      MassSignProjT ->Rebin(2);
      MassSignProjPU->Rebin(2);


      //make histo that will contains the shapes for limit
      string shapeFilePath = "/tmp/shape_"+signal+".root";
      TFile* out = new TFile(shapeFilePath.c_str(),"RECREATE");   
      out->cd();
      out->mkdir(CutIndexStr);
      out->cd(CutIndexStr);

      //save histo into the file   
      saveHistoForLimit(MassDataProj, "data_obs","");
      saveHistoForLimit(MassPredProj, "pred", "P");
      saveHistoForLimit(MassSignProj, signal, "S");
      saveVariationHistoForLimit(MassSignProj, MassSignProjP , signal, "mom");
      saveVariationHistoForLimit(MassSignProj, MassSignProjI , signal, "ias");
      saveVariationHistoForLimit(MassSignProj, MassSignProjM , signal, "ih");
      saveVariationHistoForLimit(MassSignProj, MassSignProjHUp, signal, "hip");
      saveVariationHistoForLimit(MassSignProj, MassSignProjT , signal, "tof");
      saveVariationHistoForLimit(MassSignProj, MassSignProjPU, signal, "pu");

      //close the output file
      out->Close();
   }

   //Need to set what the systematic uncertainty on signal is
   bool IsNeutral = (signal.find("N")!=std::string::npos);

   //Reconstruction Efficiency uncertainty
   double UncEffP=(EffP-Eff)/Eff;
   double UncEffI=(EffI-Eff)/Eff;
   double UncEffM=(EffM-Eff)/Eff;
   double UncEffH= (EffHUp-Eff)/Eff;
   double UncEffPU=(EffPU-Eff)/Eff;
   double UncEffT=(EffT-Eff)/Eff;
   double UncEffRe  = -0.02;
   double UncEffTr  = -0.13;
   double UncEffMB  = 0.0;

   //Reset Reco and dEdx uncertainty for fractional as dedicated samples for this
   if(TypeMode==5){
   if(signal.find("1o3")!=string::npos) {UncEffI= -0.25; UncEffRe=0.;}
   if(signal.find("2o3")!=string::npos) {UncEffI= -0.10; UncEffRe=0.;}
   }
   //Reset MB for mCHAMP
//   if((signal.find("Q2")!=string::npos && signal.find("Q2o3")==string::npos) || signal.find("Q3")!=string::npos || signal.find("Q4")!=string::npos || signal.find("Q5")!=string::npos) UncEffMB=-0.2;
//   if(signal.find("DY_Q")!=string::npos && signal.find("o3")==string::npos) UncEffMB=-0.2;
   if(signal.find("DY")!=string::npos && signal.find("o3")==string::npos) UncEffMB=-0.2;

   /*
   //Trigger efficiency uncertainty
   if(SQRTS==7) {
     //Numbers here are 0.05 for muon reconstruction uncertainty, 0.02 for muon trigger synchronization, 0.05 for MET charge suppresses
     //and 0.02 for MET non charge suppressed
     if(signal.find("1o3")!=string::npos) UncEffTr = -1*sqrt(0.15*0.15 + 0.02*0.02 + 0.05*0.05);
     else if(signal.find("2o3")!=string::npos) UncEffTr = -1*sqrt(0.03*0.03 + 0.02*0.02 + 0.05*0.05);
     else if(IsNeutral) UncEffTr = -0.05;
     else UncEffTr = -1*sqrt(0.05*0.05 + 0.02*0.02 + 0.02*0.02);
   }
   else {
     //Numbers here are 0.05 for muon reconstruction uncertainty, 0.04 for muon trigger synchronization, 0.01 for MET (charge suppresses or not)
     if(IsNeutral) UncEffTr = -0.01;
     else if(signal.find("1o3")!=string::npos) UncEffTr = -1*sqrt(0.15*0.15 + 0.04*0.04 + 0.05*0.05);
     else if(signal.find("2o3")!=string::npos) UncEffTr = -1*sqrt(0.03*0.03 + 0.04*0.04 + 0.05*0.05);
     else UncEffTr = -1*sqrt(0.05*0.05 + 0.04*0.04 + 0.01*0.01);
   }*/

//   printf("uncertainties = %f %f %f %f %f %f %f\n", UncEffP, UncEffI, UncEffPU, UncEffT, UncEffRe, UncEffTr, UncEffMB);
   if(isnan((float)UncEffPU))UncEffPU=0.0;
 
   double SignalUnc = 1 + sqrt(UncEffP*UncEffP + UncEffI*UncEffI + UncEffM*UncEffM + UncEffH*UncEffH + UncEffPU*UncEffPU + UncEffT*UncEffT + UncEffTr*UncEffTr + UncEffRe*UncEffRe + UncEffMB*UncEffMB);
   result.TotalUnc = SignalUnc-1;

//   printf("SIGNAL UNCERTAINTY = %f\n",SignalUnc);
   //build the combine datacard, the same code is used both for cut&count and shape base
   char TypeStr[255]; sprintf(TypeStr,"Type%i", TypeMode);
   string JobName = TypeStr+signal;
   string datacardPath = "/tmp/shape_"+JobName+".dat";
   makeDataCard(datacardPath,string("shape_")+JobName+".root", TypeStr,signal, NData, NPred, 1.0+(Shape?RescaleError:NPredErr/NPred), NSign, 1.0+fabs(EffErr/Eff), SignalUnc, Shape);

   char massStr[255]; sprintf(massStr,"%.0f",result.Mass);
   string test = massStr + signal;
   if(getSignificance && Temporary){
      if(NPred<0.001) NPred=0.001;
      double SignifValue=0.0; double PrevSignifValue=0; double Strength=0.1*(3*sqrt(NPred)/NSign);  if(result.XSec_5Sigma>0 && result.XSec_5Sigma<1E48)Strength=result.XSec_5Sigma / (SignalScaleFactor/result.LInt);
//      double SignifValue=0.0;double Strength=0.0005;  if(result.XSec_5Sigma>0 && result.XSec_5Sigma<1E50)Strength=result.XSec_5Sigma/result.XSec_Th;
      double previousXSec_5Sigma=result.XSec_5Sigma; result.XSec_5Sigma = -1;
      //find signal strength needed to get a 5sigma significance
      unsigned int l=0;
      double CountDecrease=0;
      for(l=0;l<10 && pointMayBeOptimal;l++){
         PrevSignifValue = SignifValue;
         SignifValue = computeSignificance(datacardPath, true, (JobName), massStr, Strength);
         printf("%i SIGNAL STRENGTH = %E --> SIGNIFICANCE=%E (countDecrease=%f)\n",l, Strength,SignifValue,CountDecrease);fflush(stdout);

         if(SignifValue<=PrevSignifValue || SignifValue<=0){CountDecrease++;}else{CountDecrease=0;}
         if(CountDecrease>=3){result.XSec_5Sigma  = 1E49; break;}

         //we found the signal strength that lead to a significance close enough to the 5sigma to stop the loop 
         //OR we know that this point is not going to be a good one --> can do a coarse approximation since the begining
//         if(fabs(SignifValue-5)<0.75 || (fastOptimization && Strength>=previousXSec_5Sigma && SignifValue<5)){
         if(fabs(SignifValue-5)<1.00 || (fastOptimization && previousXSec_5Sigma>=0 && Strength>=previousXSec_5Sigma*(SignalScaleFactor/result.LInt) && SignifValue<5 && SignifValue>=0)){
            if(fabs(SignifValue-5)<1.00){printf("Full Converge\n");
            }else{printf("Fast Converge\n");}
            result.XSec_5Sigma  = Strength * (5/SignifValue) * (SignalScaleFactor/result.LInt);//xsection in pb
            printf("XSection for 5sigma discovery = %f = %f * %f * (%f/%f)\n",result.XSec_5Sigma, Strength,  (5/SignifValue), SignalScaleFactor, result.LInt);
            break;
         }

         //Not yet at the right significance, change the strength to get close
         if(isinf((float)SignifValue)){Strength/=5;                 continue;} //strength is already way too high
         if(SignifValue<=0           ){Strength*=10;                continue;} //strength is already way too low
         if(SignifValue>5            ){Strength*=std::max( 0.1,(4.95/SignifValue)); continue;} //5/significance could be use but converge faster with 4.9
         if(SignifValue<5            ){Strength*=std::min(10.0,(5.05/SignifValue)); continue;} //5/significance could be use, but it converges faster with 5.1
         break;                    
      }
   }

   if(getXsection){
      //prepare and run the script that will run the external "combine" tool from the Higgs group
      //If very low background range too small, set limit at 0.001.  Only affects scanning range not final limit
      if(NPred<0.001) NPred=0.001;
      string CodeToExecute = "cd /tmp/;";
      CodeToExecute += "combine -M HybridNew  -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat > shape_" + JobName + ".log;";
      CodeToExecute += "cd $OLDPWD; cp /tmp/shape_" + JobName + ".* " + InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/." + ";";
      CodeToExecute += "cd /tmp/;";
      CodeToExecute += "combine -M HybridNew -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.500 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.160 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.840 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.025 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.975 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "hadd -f higgsCombine"+JobName+".HybridNew.mH"+massStr+".Merged.root higgsCombine"+JobName+".HybridNew.mH"+massStr+"*.root";
      system(CodeToExecute.c_str());

      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      TFile* file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+".HybridNew.mH"+massStr+".Merged.root").c_str());
      if(!file || file->IsZombie()){printf("SkipThisPoint --> combine output files do not exist\n"); return false;}
      TTree* tree = (TTree*)file->Get("limit");
      if(!tree){printf("SkipThisPoint --> limit file does not containing a valid TTree\n"); return false;}
      double Tmass, Tlimit, TlimitErr; float TquantExp;
      tree->GetBranch("mh"              )->SetAddress(&Tmass    );
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("limitErr"        )->SetAddress(&TlimitErr);
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        printf("Quantile=%f --> Limit = %f\n", TquantExp, Tlimit*(SignalScaleFactor/result.LInt));
        if(TquantExp==0.025f){ result.XSec_Exp2Down = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==0.160f){ result.XSec_ExpDown  = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==0.500f){ result.XSec_Exp      = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==0.840f){ result.XSec_ExpUp    = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==0.975f){ result.XSec_Exp2Up   = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==-1    ){ result.XSec_Obs      = Tlimit*(SignalScaleFactor/result.LInt);
        }else{printf("Quantil %f unused by the analysis --> check the code\n", TquantExp);
        }
      }
      file->Close();
   }

   if(!Temporary && getSignificance){
     result.Significance = computeSignificance(datacardPath, false, (JobName), massStr, 1.0);
   }


   //makePlots (only for shape based analysis)
   if(Shape && !Temporary){ 
      TCanvas* c1 = new TCanvas("c1", "c1",1200,600);
      c1->Divide(2,1);
      (c1->cd(1))->SetLogy(true);

      double Max = 2.0 * std::max(std::max(MassDataProj->GetMaximum(), MassPredProj->GetMaximum()), MassSignProj->GetMaximum());
      double Min = 0.01;
      MassSignProj->SetStats(kFALSE);
      MassSignProj->SetMaximum(Max);
      MassSignProj->SetMinimum(Min);
      MassSignProj->GetXaxis()->SetRangeUser(0,1400);
      MassSignProj->GetXaxis()->SetNdivisions(505,"X");
      MassSignProj->SetMarkerStyle(21);
      MassSignProj->SetMarkerColor(kBlue-9);
      MassSignProj->SetMarkerSize(1.5);
      MassSignProj->SetLineColor(kBlue-9);
      MassSignProj->SetFillColor(kBlue-9);
      MassSignProj->Draw("HIST");

      TH1D* MassPredSignProj = (TH1D*)MassPredProj->Clone("predWithSign");
      MassPredSignProj->Add(MassSignProj);
      MassPredSignProj->SetMarkerStyle(0);
      MassPredSignProj->SetLineColor(kBlue-10);
      MassPredSignProj->SetLineStyle(1);
      MassPredSignProj->SetLineWidth(4);
      MassPredSignProj->SetFillColor(0);
      MassPredSignProj->SetFillStyle(0);
      MassPredSignProj->Draw("same HIST C");

      MassPredProj->SetBinContent(MassPredProj->GetNbinsX(), MassPredProj->GetBinContent(MassPredProj->GetNbinsX()) + MassPredProj->GetBinContent(MassPredProj->GetNbinsX()+1));
      MassPredProj->SetMarkerStyle(22);
      MassPredProj->SetMarkerColor(2);
      MassPredProj->SetMarkerSize(1.5);
      MassPredProj->SetLineColor(2);
      MassPredProj->SetFillColor(8);
      MassPredProj->SetFillStyle(3001);
      MassPredProj->Draw("same E3");

      MassDataProj->SetBinContent(MassDataProj->GetNbinsX(), MassDataProj->GetBinContent(MassDataProj->GetNbinsX()) + MassDataProj->GetBinContent(MassDataProj->GetNbinsX()+1));
      MassDataProj->SetMarkerStyle(20);
      MassDataProj->SetMarkerColor(1);
      MassDataProj->SetMarkerSize(1.0);
      MassDataProj->SetLineColor(1);
      MassDataProj->SetFillColor(0);
      MassDataProj->Draw("same E1");

      TLegend* leg = new TLegend(0.40,0.75,0.80,0.93);
      leg->SetHeader(NULL);
      leg->SetFillColor(0);
      leg->SetFillStyle(0);
      leg->SetBorderSize(0);
      leg->AddEntry(MassDataProj,"Data", "P");
      leg->AddEntry(MassPredProj,"Prediction", "FP");
      leg->AddEntry(MassSignProj,signal.c_str(), "F");
      leg->AddEntry(MassPredSignProj,(string("Prediction + ")+signal).c_str(), "L");
      leg->Draw();

      (c1->cd(2))->SetLogy(true);
//      (c1->cd(2))->SetGridy(true);

      TH1* MassSignProjRatio  = (TH1D*)MassSignProj ->Clone(signal.c_str() );  MassSignProjRatio ->SetLineColor(1); MassSignProjRatio ->SetMarkerColor(1); MassSignProjRatio ->SetMarkerStyle(0); 
      TH1* MassSignProjPRatio = (TH1D*)MassSignProjP->Clone("mom");            MassSignProjPRatio->SetLineColor(2); MassSignProjPRatio->SetMarkerColor(2); MassSignProjPRatio->SetMarkerStyle(20);
      TH1* MassSignProjIRatio = (TH1D*)MassSignProjI->Clone("Ias");            MassSignProjIRatio->SetLineColor(4); MassSignProjIRatio->SetMarkerColor(4); MassSignProjIRatio->SetMarkerStyle(21);
      TH1* MassSignProjMRatio = (TH1D*)MassSignProjM->Clone("Ih");             MassSignProjMRatio->SetLineColor(3); MassSignProjMRatio->SetMarkerColor(3); MassSignProjMRatio->SetMarkerStyle(22);
      TH1* MassSignProjHRatioUp   = (TH1D*)MassSignProjHUp  ->Clone("hipUp");  MassSignProjHRatioUp  ->SetLineColor(5); MassSignProjHRatioUp  ->SetMarkerColor(5); MassSignProjHRatioUp  ->SetMarkerStyle(24);
      TH1* MassSignProjHRatioDown = (TH1D*)MassSignProjHDown->Clone("hipDown");MassSignProjHRatioDown->SetLineColor(5); MassSignProjHRatioDown->SetMarkerColor(5); MassSignProjHRatioDown->SetMarkerStyle(24);
      TH1* MassSignProjTRatio = (TH1D*)MassSignProjT->Clone("TOF");            MassSignProjTRatio->SetLineColor(8); MassSignProjTRatio->SetMarkerColor(8); MassSignProjTRatio->SetMarkerStyle(23);
      TH1* MassSignProjLRatio = (TH1D*)MassSignProjPU->Clone("pu");            MassSignProjLRatio->SetLineColor(6); MassSignProjLRatio->SetMarkerColor(6); MassSignProjLRatio->SetMarkerStyle(33);

      //MassSignProjPRatio->Divide(MassSignProjPRatio, MassSignProjRatio,1,1, "B");
      //MassSignProjIRatio->Divide(MassSignProjIRatio, MassSignProjRatio,1,1, "B");
      //MassSignProjMRatio->Divide(MassSignProjMRatio, MassSignProjRatio,1,1, "B");
      //MassSignProjTRatio->Divide(MassSignProjTRatio, MassSignProjRatio,1,1, "B");
      //MassSignProjLRatio->Divide(MassSignProjLRatio, MassSignProjRatio,1,1, "B");
      //MassSignProjRatio ->Divide(MassSignProjRatio , MassSignProjRatio,1,1, "B");

      MassSignProjRatio->SetStats(kFALSE);
      MassSignProjRatio->SetFillColor(0);
      MassSignProjRatio->SetLineWidth(2);
      MassSignProjRatio->GetXaxis()->SetRangeUser(0,1400);
      MassSignProjRatio->GetXaxis()->SetNdivisions(505,"X");
      MassSignProjRatio->GetYaxis()->SetNdivisions(505,"X");
      MassSignProjRatio->SetMaximum(MassSignProjPRatio->GetMaximum()*2.0);
      MassSignProjRatio->SetMinimum(Min);
      //MassSignProjRatio->SetMaximum(2);//Max);
      //MassSignProjRatio->SetMinimum(0);//Min);
      //MassSignProjRatio->Reset(); //use this histogram as a framework only
      MassSignProjRatio->Draw("HIST");
      //TLine l1(0,1,1400,1);l1.SetLineColor(1); l1.SetLineWidth(2); l1.Draw("same");
      MassSignProjPRatio->Draw("same E1");
      MassSignProjIRatio->Draw("same E1");
      MassSignProjMRatio->Draw("same E1");
      MassSignProjTRatio->Draw("same E1");
      MassSignProjLRatio->Draw("same E1");

      TLegend* leg2 = new TLegend(0.45,0.65,0.85,0.93);
      leg2->SetHeader("Variations");
      leg2->SetFillColor(0);
      leg2->SetFillStyle(0);
      leg2->SetBorderSize(0);
      leg2->AddEntry(MassSignProjRatio,signal.c_str(), "L");
      leg2->AddEntry(MassSignProjPRatio,"momentum", "P");
      leg2->AddEntry(MassSignProjIRatio,"Ias", "P");
      leg2->AddEntry(MassSignProjMRatio,"Ih", "P");
      leg2->AddEntry(MassSignProjTRatio,"tof", "P");
      leg2->AddEntry(MassSignProjLRatio,"pu", "P");
      leg2->Draw();

      SaveCanvas(c1, InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/shape", signal, true);
      delete leg2; delete leg; delete c1; delete MassPredSignProj;
      delete MassSignProjRatio; delete MassSignProjPRatio; delete MassSignProjIRatio; delete MassSignProjMRatio; delete MassSignProjTRatio; delete MassSignProjLRatio;
   }

   //all done, clean everything and return true
   delete MassDataProj; delete MassPredProj; delete MassSignProj; delete MassSignProjP; delete MassSignProjI; delete MassSignProjM; delete MassSignProjHUp; delete MassSignProjT; delete MassSignProjPU; delete MassSignProjHDown;
   return true;
}


bool Combine(string InputPattern, string signal1, string signal2, int* OptCutIndex){
//   CurrentSampleIndex        = JobIdToIndex(signal, samples); if(CurrentSampleIndex<0){  printf("There is no signal corresponding to the JobId Given\n");  return false;  }
//   int s = CurrentSampleIndex;
   GetSampleDefinition(samples);
   int SampleIndex1 = JobIdToIndex(signal1,samples),
       SampleIndex2 = JobIdToIndex(signal2,samples);

   int TypeMode = TypeFromPattern(InputPattern);
   
   string outpath = InputPattern + "/"+SHAPESTRING+EXCLUSIONDIR+"/";
   MakeDirectories(outpath);
   string JobName;
   stAllInfo result11;
   stAllInfo result12;
   stAllInfo result;
   string signal;
   string EXCLUSION1;
   string EXCLUSION2;
   string signal11;
   string signal12;
   string CodeToExecute = "combineCards.py ";
   char massStr[255];
   
   if (signal1.find("7TeV")!=string::npos && signal2.find("8TeV")!=string::npos){
      //Get Optimal cut from sample11
      EXCLUSION1 = "EXCLUSION7TeV";
      EXCLUSION2 = "EXCLUSION8TeV";
      signal11   = signal1;
      signal12   = signal2;
      result11   = stAllInfo(InputPattern+EXCLUSION1+signal1+".txt");
      //Get Optimal cut from sample12
      result12 =  stAllInfo(InputPattern+EXCLUSION2+signal2+".txt");

      result = stAllInfo(InputPattern+EXCLUSION2+signal2+".txt");
      sprintf(massStr,"%.0f",result.Mass);

      signal = signal1;
      if(signal.find("_7TeV")!=string::npos){signal.replace(signal.find("_7TeV"),5, "");}
      char TypeStr[100] ;sprintf(TypeStr,"Type%i", TypeMode);
      JobName = TypeStr+signal;

      FILE* pFileTmp = NULL;

      bool is7TeVPresent = true;
      pFileTmp = fopen((InputPattern+EXCLUSION1+"shape_"+(TypeStr+signal1)+".dat").c_str(), "r");
      if(!pFileTmp){is7TeVPresent=false;}else{fclose(pFileTmp);}
      if(TypeMode==3) is7TeVPresent=false;

      bool is8TeVPresent = true;
      pFileTmp = fopen((InputPattern+EXCLUSION2+"shape_"+(TypeStr+signal2)+".dat").c_str(), "r");
      if(!pFileTmp){is8TeVPresent=false;}else{fclose(pFileTmp);}


      if(is7TeVPresent)CodeToExecute+="   " + InputPattern+EXCLUSION1+"shape_"+(TypeStr+signal1)+".dat ";
      if(is8TeVPresent)CodeToExecute+="   " + InputPattern+EXCLUSION2+"shape_"+(TypeStr+signal2)+".dat ";
   }

   // COMBINE2016
   else if (signal1.find("W13TeV16")!=string::npos && signal2.find("W13TeV16G")!=string::npos){
      size_t toBreak1   = signal1.find("W13TeV");
      size_t toBreak2   = signal2.find("W13TeV");
      EXCLUSION1 = "/EXCLUSION"+signal1.substr(toBreak1+1, signal1.size()-toBreak1-1);
      EXCLUSION2 = "/EXCLUSION"+signal2.substr(toBreak2+1, signal2.size()-toBreak2-1);
      signal11   = signal1.substr(0, toBreak1);
      signal12   = signal2.substr(0, toBreak2);
      SampleIndex1 = JobIdToIndex(ReplacePartOfString(signal11, "13TeV", "13TeV16"), samples);
      SampleIndex2 = JobIdToIndex(ReplacePartOfString(signal12, "13TeV16", "13TeV16G"), samples);
      //Get Optimal cut from sample11
      std::cout << "Loading " << InputPattern+EXCLUSION1+"/"+ReplacePartOfString(signal11,"13TeV", "13TeV16")+".txt" << std::endl;
      result11 =  stAllInfo(InputPattern+EXCLUSION1+"/"+ReplacePartOfString(signal11,"13TeV", "13TeV16")+".txt");
      //Get Optimal cut from sample22
      std::cout << "Loading " << InputPattern+EXCLUSION2+"/"+ReplacePartOfString(signal12,"13TeV16", "13TeV16G")+".txt" << std::endl;
      result12 =  stAllInfo(InputPattern+EXCLUSION2+"/"+ReplacePartOfString(signal12, "13TeV16", "13TeV16G")+".txt");
      char TypeStr[100] ;sprintf(TypeStr,"Type%i", TypeMode);
      JobName = TypeStr+signal11;

      result =  stAllInfo(InputPattern+EXCLUSION2+"/"+ReplacePartOfString(signal12, "13TeV16", "13TeV16G")+".txt");
      sprintf(massStr,"%.0f",result.Mass);


      signal = signal11;
      signal = ReplacePartOfString(signal, "_13TeV16G", "");
      signal = ReplacePartOfString(signal, "_13TeV16" , "");
      signal = ReplacePartOfString(signal, "_13TeV15" , "");
      signal = ReplacePartOfString(signal, "_13TeV"   , "");

      FILE* pFileTmp = NULL;

      bool is2016Present = true;
      printf("Accessing the file %s\n", (InputPattern+EXCLUSION1+"/shape_"+(TypeStr+ReplacePartOfString(signal11, "13TeV", "13TeV16"))+".dat").c_str());
      pFileTmp = fopen((InputPattern+EXCLUSION1+"/shape_"+(TypeStr+ReplacePartOfString(signal11, "13TeV", "13TeV16"))+".dat").c_str(), "r");
      if(!pFileTmp){is2016Present=false;}else{fclose(pFileTmp);}
//      if(TypeMode!=0 && TypeMode!=2) is2016Present=false;

      bool is2016GPresent = true;
      printf("Accessing the file %s\n", (InputPattern+EXCLUSION2+"/shape_"+(TypeStr+ReplacePartOfString(signal12, "13TeV16", "13TeV16G"))+".dat").c_str());
      pFileTmp = fopen((InputPattern+EXCLUSION2+"/shape_"+(TypeStr+ReplacePartOfString(signal12, "13TeV16", "13TeV16G"))+".dat").c_str(), "r");
      if(!pFileTmp){is2016GPresent=false;}else{fclose(pFileTmp);}


      if(is2016Present)CodeToExecute+="   " + InputPattern+EXCLUSION1+"/shape_"+(TypeStr+ReplacePartOfString(signal11, "13TeV", "13TeV16"))+".dat ";
      if(is2016GPresent){
         CodeToExecute+="   " + InputPattern+EXCLUSION2+"/shape_"+(TypeStr+ReplacePartOfString(signal12,"13TeV16","13TeV16G"))+".dat ";
//	 string PreExecuteCode = string("sed -i \'s:13TeV16G:13TeV16:g\' ")+InputPattern+EXCLUSION2+"/shape_"+(TypeStr+ReplacePartOfString(signal12,"13TeV16","13TeV16G"))+".dat";
//	 printf ("Renaming the signal in the datacard:\n%s\n", PreExecuteCode.c_str());
//	 system(PreExecuteCode.c_str());
      }
   }
   if (result11.MassCut != result12.MassCut/* && OptCutIndex != nullptr*/ && (SampleIndex1 < 0 || SampleIndex2 < 0)){printf("These two signals are not listed!\n");return false;} 
   else if (result11.MassCut != result12.MassCut && SampleIndex1 >= 0 && SampleIndex2 >=0 && SampleIndex1 != SampleIndex2){
      string analysisPath = "/home/ucl/cp3/jzobec/Run2HSCP16/Run2HSCP16_v4/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode";
      FILE* fdebug = fopen("DifferentMassCuts.txt", "a");
      fprintf (fdebug, "%s\t%s\n", (ReplacePartOfString(signal12, "13TeV16", "13TeV16G")+".txt").c_str(), (ReplacePartOfString(signal11, "13TeV", "13TeV16")+".txt").c_str());
      fclose(fdebug);

      /* take the smaller of the two mass cuts */
      bool SecondIsGreater = (result11.MassCut<result12.MassCut);
      double MassCut2 = SecondIsGreater?result11.MassCut:result12.MassCut;

      // runCombine must be rerun for a different mass cut for the less permissive datacard
      // TOFIXCOMBINE
      string Data = SecondIsGreater?ReplacePartOfString(signal12,"13TeV16","13TeV16G"):ReplacePartOfString(signal11, "13TeV", "13TeV16");
      string tmp = Data;
      double SQRTS_OLD = SQRTS;
      if (Data.find("7TeV")!=string::npos) {tmp = "7TeV"; SQRTS = 7.0;}
      else if (Data.find("8TeV")!=string::npos) {tmp = "8TeV"; SQRTS = 8.0;}
      else if (Data.find("13TeV16G")!=string::npos) {tmp = "13TeV16G"; SQRTS = 13167.0;}
      else if (Data.find("13TeV16")!=string::npos) {tmp = "13TeV16"; SQRTS = 1316.0;}
      else if (Data.find("13TeV")!=string::npos) {tmp = "13TeV"; SQRTS = 13.0;}
      Data = "Data"+tmp;

      printf ("DATA = %s\n", Data.c_str());
      printf ("Rerunning %s (sample index = %d) to change its mass cut from %.0lf to %.0lf ...\n", SecondIsGreater?ReplacePartOfString(signal12,"13TeV16", "13TeV16G").c_str():ReplacePartOfString(signal11,"13TeV","13TeV16").c_str(), SecondIsGreater?SampleIndex2:SampleIndex1, SecondIsGreater?result12.MassCut:result11.MassCut, MassCut2);
      if (SecondIsGreater) result12.MassCut = result11.MassCut;
      else result11.MassCut = result12.MassCut;
      CurrentSampleIndex = SecondIsGreater?SampleIndex2:SampleIndex1;
      MinRange = MassCut2;
      result.MassCut      = MinRange;
/*
      printf("2016G Data ...\n");
      Data = "Data13TeV16G"; SQRTS=13167.0; EXCLUSIONDIR=EXCLUSIONDIR_SAVE+"13TeV16G";
      Optimize(InputPattern, Data, ReplacePartOfString(signal, "13TeV16", "13TeV"), SHAPESTRING!="", true);
*/
      TFile*InputFile     = new TFile((InputPattern + "Histos.root").c_str());
      TH1D* H_P           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_P");
      TH1D* H_D           = (TH1D*)GetObjectFromPath(InputFile, Data+"/H_D");
      TH2D* MassData      = (TH2D*)GetObjectFromPath(InputFile, Data+"/Mass");
      TH2D* MassPred      = (TH2D*)GetObjectFromPath(InputFile, Data+"/Pred_Mass");
      TH2D* MassSign      = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass" );
      if(!MassSign){printf("The sample %s is not present in the root file, returns\n", ReplacePartOfString(signal11, "13TeV", "13TeV16").c_str());return false;}
      TH2D* MassSignP     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystP"    );
      TH2D* MassSignI     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystI"    );
      TH2D* MassSignM     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystM"    );
      TH2D* MassSignHUp   = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystHUp"  );
      TH2D* MassSignHDown = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystHDown");
      TH2D* MassSignT     = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystT"    );
      TH2D* MassSignPU    = (TH2D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name + "/Mass_SystPU"   );
      TH1D* TotalE        = (TH1D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name+"/TotalE" );
      TH1D* TotalEPU      = (TH1D*)GetObjectFromPath(InputFile, samples[CurrentSampleIndex].Name+"/TotalEPU" );

      double norm  = samples[CurrentSampleIndex].XSec*(SecondIsGreater?result12.LInt:result11.LInt)/TotalE  ->Integral(); //normalize the samples to the actual lumi used for limits
      double normPU= samples[CurrentSampleIndex].XSec*(SecondIsGreater?result12.LInt:result11.LInt)/(TotalEPU->Integral()>0?TotalEPU->Integral():TotalE->Integral());

      MassSign      ->Scale(norm);
      MassSignP     ->Scale(norm);
      MassSignI     ->Scale(norm);
      MassSignM     ->Scale(norm);
      MassSignHUp   ->Scale(norm);
      MassSignHDown ->Scale(norm);
      MassSignT     ->Scale(norm);
      MassSignPU    ->Scale(normPU);



      // so far the shape analysis is disabled here! change the lines here to change it
      bool shape = false;
      string signalToReturn = SecondIsGreater?ReplacePartOfString(signal12,"13TeV16", "13TeV16G"):ReplacePartOfString(signal11,"13TeV","13TeV16");
      if(TypeMode<=2){runCombine(false, true, true, InputPattern, signalToReturn, SecondIsGreater?result12.Index:result11.Index, false, shape, SecondIsGreater?result12:result11, MassData, MassPred, MassSign, MassSignP, MassSignI, MassSignM, MassSignHUp, MassSignHDown, MassSignT, MassSignPU);
      }else          {runCombine(false, true, true, InputPattern, signalToReturn, SecondIsGreater?result12.Index:result11.Index , false, shape, SecondIsGreater?result12:result11, H_D, H_P, MassSign, MassSignP, MassSignI, MassSignM, MassSignHUp, MassSignHDown, MassSignT, MassSignPU);
      }
     
      //overwrite the last result and update it to the new mass cut
      printf ("New NPred = % .2e\n", SecondIsGreater?result12.NPred:result11.NPred);
      string OverwritePath = InputPattern+(SecondIsGreater?EXCLUSION2:EXCLUSION1)+"/"+(SecondIsGreater?ReplacePartOfString(signal12,"13TeV16","13TeV16G"):ReplacePartOfString(signal11,"13TeV", "13TeV16"))+".txt";
      printf ("Results written in %s.\n", OverwritePath.c_str());
      if (SecondIsGreater)
         result12.Save(OverwritePath);
      if (!SecondIsGreater)
         result11.Save(OverwritePath);
      if (SQRTS != SQRTS_OLD) SQRTS = SQRTS_OLD;

      delete H_P           ;
      delete H_D           ;
      delete MassData      ;
      delete MassPred      ;
      delete MassSign      ;
      delete MassSignP     ;
      delete MassSignI     ;
      delete MassSignM     ;
      delete MassSignHUp   ;
      delete MassSignHDown ;
      delete MassSignT     ;
      delete MassSignPU    ;
      delete TotalE        ;
      delete TotalEPU      ;
      delete InputFile;
   }

   double NSign1 = result11.NSign/(result11.XSec_Th*result11.LInt),
          NSign2 = result12.NSign/(result12.XSec_Th*result12.LInt),
          SignalScaleFactor1 = 1.0,
          SignalScaleFactor2 = 1.0,
	  SignalScaleFactor;

   for(unsigned int i=0;i<20 && NSign1<1e-1; i++) {SignalScaleFactor1*=10.0; NSign1*=10.0;}
   for(unsigned int i=0;i<20 && NSign2<1e-1; i++) {SignalScaleFactor2*=10.0; NSign2*=10.0;}
   if (SignalScaleFactor1 != SignalScaleFactor2){
      printf ("Warning! SignalScaleFactor1 = %.2e is different from SignalScaleFactor2 = %.2e\n", SignalScaleFactor1, SignalScaleFactor2);
      printf ("Rescaling the signal rates to the same level in both datacards before combining!\n");
      double NPred1        = (result11.NPred>=0.001)?result11.NPred:0.001,
	     NPredErrRatio = (result11.NPred>=0.001)?(result11.NPredErr/result11.NPred):0.001,
             EffErrRatio   = result11.EffE/result11.Eff,
             SignalUnc     = result11.TotalUnc+1,
             NData1        = result11.NData;

      string datacardPath1 = InputPattern+EXCLUSION1+string("/shape_Type")+((char) (48 + TypeMode))+(ReplacePartOfString(signal11, "13TeV", "13TeV16"))+".dat ",
             datacardPath2 = InputPattern+EXCLUSION2+string("/shape_Type")+((char) (48 + TypeMode))+(ReplacePartOfString(signal12,"13TeV16","13TeV16G"))+".dat";
      bool SecondIsGreater = (SignalScaleFactor1<SignalScaleFactor2);
      double ratio    = SignalScaleFactor2/SignalScaleFactor1;
      if (!SecondIsGreater) {
         datacardPath1         = datacardPath2;
         NSign1                = NSign2;
         NData1                = result12.NData;
         NPred1                = (result12.NPred>=0.001)?result12.NPred:0.001;
	 NPredErrRatio         = (result12.NPred>=0.001)?(result12.NPredErr/result12.NPred):0.001;
         EffErrRatio           = result12.EffE/result12.Eff,
         SignalUnc             = result12.TotalUnc+1;
	 SignalScaleFactor     = SignalScaleFactor1;
	 signal11              = signal12;
	 ratio                 = SignalScaleFactor1/SignalScaleFactor2;
      } else SignalScaleFactor = SignalScaleFactor2;

      std::cout << "Readjusting the rate in " + datacardPath1 + " from " << NSign1 << " to " << ratio*NSign1 << std::endl;
      // new datacard for this case -> the combined datacard emerging from these two datacards will be correct
      makeDataCard(datacardPath1,string("shape_")+JobName+".root", string("Type")+((char) (48 + TypeMode)), signal11, NData1, NPred1, 1.0+NPredErrRatio, NSign1*ratio, 1.0+fabs(EffErrRatio), SignalUnc, false);

      printf ("Datacard %s has been rewritten!\n", datacardPath1.c_str());
   } else SignalScaleFactor = SignalScaleFactor1;


   CodeToExecute+=" > " + outpath+"shape_"+JobName+".dat";

   system(CodeToExecute.c_str());   
   printf("%s \n",CodeToExecute.c_str());

   printf ("Signal scale factor = %.2e\n", SignalScaleFactor);

//   result.XSec_Th = 1.0; // FIXME NOT ADEQUATE FOR COMBINING SAMPLES AT THE SAME ENERGY! -- take whatever is in any 13TeV sample
   //Muon only uses just 2012
   if(TypeMode==3 && signal1.find("7TeV")!=string::npos && signal2.find("8TeV")!=string::npos) {
     result.XSec_Obs=result12.XSec_Obs/result12.XSec_Th;
     result.XSec_Exp=result12.XSec_Exp/result12.XSec_Th;
     result.XSec_ExpUp=result12.XSec_ExpUp/result12.XSec_Th;
     result.XSec_ExpDown=result12.XSec_ExpDown/result12.XSec_Th;
     result.XSec_Exp2Up=result12.XSec_Exp2Up/result12.XSec_Th;
     result.XSec_Exp2Down=result12.XSec_Exp2Down/result12.XSec_Th;
     result.Save(InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/"+signal+".txt");
     return true;
   }

   result.LInt  = result11.LInt  + result12.LInt ;
   result.NSign = result11.NSign + result12.NSign;
   result.NData = result11.NData + result12.NData;
   result.NPred = result11.NPred + result12.NPred;
   result.NPredErr = sqrt(pow(result11.NPredErr,2) + pow(result12.NPredErr,2));

   //compute combined significance
   double NPred = result.NPred;
   double NSign = result.NSign;
   if (SignalScaleFactor1 != SignalScaleFactor2){
      printf ("\nMaking a note in ListOfMismatches.log!\n");
      string analysisPath = "/home/ucl/cp3/jzobec/Run2HSCP16/Run2HSCP16_v4/CMSSW_8_0_30/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode";
      FILE * fdebug = fopen ((analysisPath+"/ListOfMismatches.log").c_str(), "a");
      fprintf (fdebug, "Type %d Signal %s does not match!\n", TypeMode, signal.c_str());
      fclose(fdebug);
   }

   CodeToExecute = "cp " + outpath+"shape_"+JobName+".dat /tmp/.;";   
   system(CodeToExecute.c_str());
   result.Significance = computeSignificance(outpath+"shape_"+JobName+".dat", false, signal, massStr, 1.0);
   printf("Combined Significance = %f (%s)\n", result.Significance, (outpath+"shape_"+JobName+".dat").c_str());
 
   //ALL CODE BELOW IS A BIT DIFFERENT THAN THE ONE USED IN runCombined, BECAUSE HERE WE KEEP THE RESULTS ON LIMIT IN TERMS OF SIGNAL STRENGTH (r=SigmaObs/SigmaTH)
   if(true){
      //prepare and run the script that will run the external "combine" tool from the Higgs group
      //If very low background range too small, set limit at 0.001.  Only affects scanning range not final limit
            if(NPred<0.001) NPred=0.001;
      string CodeToExecute = "cp " + outpath+"shape_"+JobName+".dat /tmp/.;";
      CodeToExecute += "cd /tmp/;";
      CodeToExecute += "combine -M HybridNew  -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat > shape_" + JobName + ".log;";
      CodeToExecute += "cd $OLDPWD; cp /tmp/shape_" + JobName + ".* " + InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/." + ";";
      CodeToExecute += "cd /tmp/;";
      CodeToExecute += "combine -M HybridNew  -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.500 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew  -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.160 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew  -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.840 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew  -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.025 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "combine -M HybridNew  -n " + JobName + " -m " + massStr + " shape_" + JobName+".dat --expectedFromGrid 0.975 >> shape_" + JobName + "Exp.log;";
      CodeToExecute += "hadd -f higgsCombine"+JobName+".HybridNew.mH"+massStr+".Merged.root higgsCombine"+JobName+".HybridNew.mH"+massStr+"*.root";
      system(CodeToExecute.c_str());

      //if all went well, the combine tool created a new file containing the result of the limit in the form of a TTree
      //we can open this TTree and access the values for the expected limit, uncertainty bands, and observed limits.
      TFile* file = TFile::Open((string("/tmp/")+"higgsCombine"+JobName+".HybridNew.mH"+massStr+".Merged.root").c_str());
      if(!file || file->IsZombie())return false;
      TTree* tree = (TTree*)file->Get("limit");
      if(!tree)return false;
      double Tmass, Tlimit, TlimitErr; float TquantExp;
      tree->GetBranch("mh"              )->SetAddress(&Tmass    );
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("limit"           )->SetAddress(&Tlimit   );
      tree->GetBranch("limitErr"        )->SetAddress(&TlimitErr);
      tree->GetBranch("quantileExpected")->SetAddress(&TquantExp);
      for(int ientry=0;ientry<tree->GetEntriesFast();ientry++){
        tree->GetEntry(ientry);
        printf("Quantile=%f --> Limit = %f\n", TquantExp, Tlimit*(SignalScaleFactor/result.LInt));
        if(TquantExp==0.025f){ result.XSec_Exp2Down = Tlimit*(SignalScaleFactor/result.LInt); // FIXME jozze -- had to rescale with LInt, or limits are wrong
        }else if(TquantExp==0.160f){ result.XSec_ExpDown  = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==0.500f){ result.XSec_Exp      = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==0.840f){ result.XSec_ExpUp    = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==0.975f){ result.XSec_Exp2Up   = Tlimit*(SignalScaleFactor/result.LInt);
        }else if(TquantExp==-1    ){ result.XSec_Obs      = Tlimit*(SignalScaleFactor/result.LInt);
        }else{printf("Quantil %f unused by the analysis --> check the code\n", TquantExp);
        }
      }
      file->Close();
   }

   //all done, save the results to file
   result.Save(InputPattern+"/"+SHAPESTRING+EXCLUSIONDIR+"/"+signal+".txt");
   return true;
}


bool useSample(int TypeMode, string sample) {
  if(TypeMode==0 && (sample=="Gluino16_f10" || sample=="Gluino16N_f10" || sample=="Stop16N" || sample=="Stop16" || sample=="DY16_Q2o3" || sample=="GMStau16" || sample=="PPStau16" ||
    sample=="Gluino16G_f10" || sample=="Gluino16GN_f10" || sample=="Stop16GN" || sample=="Stop16G" || sample=="DY16G_Q2o3" || sample=="GMStau16G" || sample=="PPStau16G"))
    return true;

  if(TypeMode==2 && (sample=="Gluino16_f10" || sample=="Gluino16_f50" || sample=="Stop16" || sample=="GMStau16" || sample=="PPStau16" || sample=="DY16_Q2o3" || sample=="DY16_Q1" ||
    sample=="Gluino16G_f10" || sample=="Gluino16G_f50" || sample=="Stop16G" || sample=="GMStau16G" || sample=="PPStau16G" || sample=="DY16G_Q2o3" || sample=="DY16G_Q1"))
    return true;

  if(TypeMode==3 && (sample=="Stop16" || sample=="Gluino16_f10" || sample=="Gluino16_f50" || sample=="Gluino16_f100" || 
    sample=="Stop16G" || sample=="Gluino16G_f10" || sample=="Gluino16G_f50" || sample=="Gluino16G_f100"))
    return true;

  if(TypeMode==4) return true;
  if(TypeMode==5) return true;
  return false;
}

string toLatex(double value) {
  char toReturn[255];
  if(value<0.00095) {
    int exponent=0;
    while (value<1) {
      exponent++;
      value*=10;
    }
    sprintf(toReturn,"$   %6.1f \\times 10^{-%i}$",value, exponent);
  }
  else sprintf(toReturn,"%.2g",value);
  string stringReturn = toReturn;
  return stringReturn;
}

