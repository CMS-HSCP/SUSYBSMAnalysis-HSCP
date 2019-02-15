// Original Author:  Loic Quertenmont

#ifndef HSCP_ANALYSIS_SAMPLE
#define HSCP_ANALYSIS_SAMPLE

class stSample{
   public:
   std::string CMSSW;
   int         Type;
   std::string Name;
   std::string FileName;
   std::string Legend;
   std::string Pileup;
   double      Mass;
   double      XSec;
   bool        MakePlot;
   float       WNC0;//weight for signal event with 0 Charged HSCP
   float       WNC1;//weight for signal event with 1 Charged HSCP
   float       WNC2;//weight for signal event with 2 Charged HSCP

   //These weights are used to reweight the gluino samples for different scenario of gluinobll fraction (f).  
   //The weights are independent on the mass and are computed as a reweighting factor w.r.t a f=10% gluino ball sample.
   //Weights have been computed in 2010 to be : 
   //for f=00% --> [0.2524 / 0.3029 , 0.4893 / 0.4955 , 0.2583 / 0.2015]
   //for f=10% --> [0.3029 / 0.3029 , 0.4955 / 0.4955 , 0.2015 / 0.2015]
   //for f=50% --> [0.5739 / 0.3029 , 0.3704 / 0.4955 , 0.0557 / 0.2015]

   stSample(){};
   stSample(std::string CMSSW_, int Type_, std::string Name_, std::string FileName_, std::string Legend_, std::string Pileup_, double Mass_, double XSec_, bool MakePlot_, float WNC0_, float WNC1_, float WNC2_){CMSSW=CMSSW_; Type=Type_; Name=Name_; FileName=FileName_; Legend=Legend_; Pileup=Pileup_; Mass=Mass_; XSec=XSec_; MakePlot=MakePlot_; WNC0=WNC0_; WNC1=WNC1_; WNC2=WNC2_;}

   int readFromFile(FILE* pFile){ 
      char line[4096];
      char* toReturn = fgets(line, 4096, pFile);
      if(!toReturn)return EOF;

      //delete whatever is found after the comment delimeter
      bool commenterFound = false;  bool inQuote=false;
      int i=0;while(toReturn[i]!='\0'){if(toReturn[i]=='"')inQuote=!inQuote; if(toReturn[i]=='#' && !inQuote){commenterFound=true; break;} i++; }
      if(commenterFound)while(toReturn[i]!='\0'){toReturn[i]='\0'; i++;}

      Name = ""; //just to make sure that the sample name is properly initialized
      char* pch=strtok(line,","); int Arg=0; string tmp; int temp;
      while (pch!=NULL){
         switch(Arg){
            case  0: tmp = pch;  CMSSW    = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  1: sscanf(pch, "%d", &Type); break;
            case  2: tmp = pch;  Name     = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  3: tmp = pch;  FileName = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  4: tmp = pch;  Legend   = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  5: tmp = pch;  Pileup   = tmp.substr(tmp.find("\"")+1,tmp.rfind("\"")-tmp.find("\"")-1); break;
            case  6: sscanf(pch, "%lf", &Mass); break;
            case  7: sscanf(pch, "%lf", &XSec); break;
            case  8: sscanf(pch, "%d", &temp); MakePlot=(temp>0); break;
            case  9: sscanf(pch, "%f", &WNC0); break;
            case 10: sscanf(pch, "%f", &WNC1); break;
            case 11: sscanf(pch, "%f", &WNC2); break;
            default:return EOF;break;
         }
         pch=strtok(NULL,",");Arg++;
      }
      return 0;
   }

   void print(FILE* pFile=stdout){
      fprintf(pFile, "%-9s, %1d, %-30s, %-40s, %-60s, %-8s, % 5g, %+12.10E, %1i, %5.3f, %5.3f, %5.3f\n", (string("\"")+CMSSW+"\"").c_str(), Type, (string("\"")+Name+"\"").c_str(), (string("\"")+FileName+"\"").c_str(), (string("\"")+Legend+"\"").c_str(), (string("\"")+Pileup+"\"").c_str(), Mass, XSec, MakePlot, WNC0, WNC1, WNC2);
   }

   std::string ModelName(){
      char strMass[255];sprintf(strMass, "_M%.0f",Mass);
      char str7TeV[]="_7TeV";
      char str8TeV[]="_8TeV";
      char str13TeV[]="_13TeV";
      char str13TeV16[]="_13TeV16";
      char str13TeV16G[]="_13TeV16G";
      string toReturn=Name;
      if(toReturn.find(strMass) !=string::npos)toReturn.erase(toReturn.find( strMass), string( strMass).length());
      if(toReturn.find(str7TeV) !=string::npos)toReturn.erase(toReturn.find( str7TeV), string( str7TeV).length());
      if(toReturn.find(str8TeV) !=string::npos)toReturn.erase(toReturn.find( str8TeV), string( str8TeV).length());
      if(toReturn.find(str13TeV16G)!=string::npos)toReturn.erase(toReturn.find(str13TeV16G), string(str13TeV16G).length());
      if(toReturn.find(str13TeV16)!=string::npos)toReturn.erase(toReturn.find(str13TeV16), string(str13TeV16).length());
      if(toReturn.find(str13TeV)!=string::npos)toReturn.erase(toReturn.find(str13TeV), string(str13TeV).length());
      return toReturn;
   }

   std::string ModelLegend(){
      char MassStr[255];sprintf(MassStr, "%.0f",Mass);
      string toReturn=Legend; toReturn.erase(toReturn.find(MassStr), string(MassStr).length());
      if(toReturn.find(" GeV/#font[12]{c}^{2}")!=string::npos)toReturn.erase(toReturn.find(" GeV/#font[12]{c}^{2}"), string(" GeV/#font[12]{c}^{2}").length());
      return toReturn;
   }

   double GetFGluinoWeight(int NChargedHSCP){ 
      if(Type!=2)return 1;
      double Weight;
      switch(NChargedHSCP){
         case 0: Weight=WNC0; break;
         case 1: Weight=WNC1; break;
         case 2: Weight=WNC2; break;
         default: return 1.0;
      }
      if(Weight<0)return 1.0;
      return Weight;
   }

};

void GetSampleDefinition(std::vector<stSample>& samples, std::string sampleTxtFile="Analysis_Samples.txt"){
      FILE* pFile = fopen(sampleTxtFile.c_str(),"r");
         if(!pFile){
		printf("Can't open %s\n","Analysis_Samples.txt,\nattempting to open from default explicit path ...\n");
		char* cmssw_base   = getenv ("CMSSW_BASE");
		char samples_path [200];
		sprintf (samples_path, "%s/src/SUSYBSMAnalysis/HSCP/test/AnalysisCode/Analysis_Samples.txt", cmssw_base);
		pFile = fopen (samples_path, "r");
		if (!pFile){printf("Still cannot open Analysis_Samples.txt ... giving up!\n"); return;}
	 }
         stSample newSample;      
         while(newSample.readFromFile(pFile)!=EOF){samples.push_back(newSample);}
      fclose(pFile);
}

void GetInputFiles(stSample sample, std::string BaseDirectory_, std::vector<std::string>& inputFiles, int period=0){
   std::vector<string> fileNames;
   char* tmp = (char*)sample.FileName.c_str();
   char* pch=strtok(tmp,";");
   while (pch!=NULL){
      fileNames.push_back(pch);
      pch=strtok(NULL,",");
   }

   for(unsigned int f=0;f<fileNames.size();f++){  
      if(fileNames[f].find('/')==0){//full file path is provided)
         if(fileNames[f].find("/store/")==0){
            inputFiles.push_back(string("root://eoscms//eos/cms") + fileNames[f]);
         }else{
            inputFiles.push_back(fileNames[f]);
         }
      }else{
         if(sample.Type>=2){ //MC Signal
           if(fileNames[f].find("7TeV")<string::npos){ //7TeV
              if (period==0) inputFiles.push_back(BaseDirectory_ + fileNames[f] + "_BX1.root");
              if (period==1) inputFiles.push_back(BaseDirectory_ + fileNames[f] + "_BX0.root");
            }else{ //8TeV
   	      inputFiles.push_back(BaseDirectory_ + fileNames[f] + ".root");
            }
         }else{ //Data or MC Background
           inputFiles.push_back(BaseDirectory_ + fileNames[f] + ".root");
         }
      }
   }
}

int JobIdToIndex(string JobId, const std::vector<stSample>& samples){
   for(unsigned int s=0;s<samples.size();s++){
      if(samples[s].Name==JobId)return s;
   }
   return -1;

   //if sample is not found, use the 7 or 8TeV sample instead...
   //replace 7TeV by 8TeV or vice versa
   if(     JobId.find("_7TeV")!=string::npos){JobId.replace(JobId.find("_7TeV"),5, "_8TeV");}
   else if(JobId.find("_8TeV")!=string::npos){JobId.replace(JobId.find("_8TeV"),5, "_7TeV");}

   for(unsigned int s=0;s<samples.size();s++){
      if(samples[s].Name==JobId)return s;
   }
   return -1;
}


void keepOnlyValidSamples(std::vector<stSample>& samples){
   for(unsigned int s=0;s<samples.size();s++){if(samples[s].Name==""){samples.erase(samples.begin()+s);s--;} }
}


void keepOnlySamplesOfTypeX(std::vector<stSample>& samples, int TypeX){
   for(unsigned int s=0;s<samples.size();s++){if(samples[s].Type!=TypeX){samples.erase(samples.begin()+s);s--;} }
}

void keepOnlyTheXtoYSamples(std::vector<stSample>& samples, unsigned int X, unsigned int Y){
   std::vector<stSample> tmp; for(unsigned int s=X;s<=Y && s<samples.size();s++){tmp.push_back(samples[s]);} samples.clear();
   for(unsigned int s=0;s<tmp.size();s++)samples.push_back(tmp[s]);
}

void keepOnlySamplesOfNameX(std::vector<stSample>& samples, string NameX){
   for(unsigned int s=0;s<samples.size();s++){if(samples[s].Name!=NameX){samples.erase(samples.begin()+s);s--;} }
}  

void keepOnlySamplesOfNamesXtoY(std::vector<stSample>& samples, std::vector<string> NamesXtoY){
    for(unsigned int s=0;s<samples.size();s++){
      bool keep=false;
      for(unsigned int i=0;i<NamesXtoY.size();i++){
	if(samples[s].Name==NamesXtoY[i]) {keep=true; break;}
      }
      if(!keep) {samples.erase(samples.begin()+s);s--;}
    }
}

void keepOnlySamplesAtSQRTS(std::vector<stSample>& samples, double SQRTS_){
   if(SQRTS_==78 || SQRTS_==87){
      std::vector<stSample> samples_tmp;
      for(unsigned int s=0;s<samples.size();s++){
         string tmpName = samples[s].Name;
         if(     tmpName.find("_7TeV")!=string::npos){tmpName.replace(tmpName.find("_7TeV"),5, "_8TeV");}
         else if(tmpName.find("_8TeV")!=string::npos){tmpName.replace(tmpName.find("_8TeV"),5, "_7TeV");}
         for(unsigned int t=s;t<samples.size();t++){
            if(samples[t].Name==tmpName){
               samples_tmp.push_back(samples[s]);
               samples_tmp.push_back(samples[t]);
            } 
         }
      }
      samples.clear();
      for(unsigned int s=0;s<samples_tmp.size();s++){samples.push_back(samples_tmp[s]);}
   }else{
      for(unsigned int s=0;s<samples.size();s++){
         string tmpName = samples[s].Name;
         if(     SQRTS_== 7 && tmpName.find("_7TeV" )==string::npos){samples.erase(samples.begin()+s);s--;}
         if(     SQRTS_== 8 && tmpName.find("_8TeV" )==string::npos){samples.erase(samples.begin()+s);s--;}
         if((SQRTS_==13 || SQRTS_==1315 || SQRTS_==1316 || SQRTS_==131516 || SQRTS_==131615) && tmpName.find("_13TeV")==string::npos){samples.erase(samples.begin()+s);s--;}
	 if(SQRTS_==13167 && tmpName.find("_13TeV16G")==string::npos){samples.erase(samples.begin()+s);s--;}
         if(SQRTS_==131677 && tmpName.find("_13TeV16")==string::npos){samples.erase(samples.begin()+s);s--;}
      }
   }
}

int numberOfMatchingSamples(string JobId, std::vector<stSample>& samples){ 
   int toReturn=0;
   for(unsigned int s=0;s<samples.size();s++){
      if(samples[s].Name==JobId)toReturn++;
   }
   return toReturn;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genertic code related to samples processing in FWLITE --> functions below will be loaded only if FWLITE compiler variable is defined

#ifdef FWLITE

unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames);
double GetSampleWeight  (const double& IntegratedLuminosityInPb=-1, const double& IntegratedLuminosityInPbBeforeTriggerChange=-1, const double& CrossSection=0, const double& MCEvents=0, int period=0);
double GetSampleWeightMC(const double& IntegratedLuminosityInPb, const std::vector<string> fileNames, const double& XSection, const double& SampleSize, double MaxEvent);
double GetPUWeight      (const fwlite::ChainEvent& ev, const std::string& pileup, double &PUSystFactor);

// loop on all the lumi blocks for an EDM file in order to count the number of events that are in a sample
// this is useful to determine how to normalize the events (compute weight)
unsigned long GetInitialNumberOfMCEvent(const vector<string>& fileNames)
{
   unsigned long Total = 0;
   for(unsigned int f=0;f<fileNames.size();f++){
      TFile *file;
      size_t place=fileNames[f].find("dcache");
      if(place!=string::npos){
 	 string name=fileNames[f];
         //name.replace(place, 7, "dcap://cmsgridftp.fnal.gov:24125");
         file = new TDCacheFile (name.c_str());
      }else{
         file = TFile::Open(fileNames[f].c_str());
      }

      fwlite::LuminosityBlock ls( file );
      for(ls.toBegin(); !ls.atEnd(); ++ls){
         fwlite::Handle<edm::MergeableCounter> nEventsTotalCounter;
         nEventsTotalCounter.getByLabel(ls,"nEventsBefSkim");
         if(!nEventsTotalCounter.isValid()){printf("Invalid nEventsTotalCounterH\n");continue;}
         Total+= nEventsTotalCounter->value;
      }
   }
   return Total;
}

// compute event weight for the signal samples based on number of events in the sample, cross section and intergated luminosity
double GetSampleWeight(const double& IntegratedLuminosityInPb, const double& IntegratedLuminosityInPbBeforeTriggerChange, const double& CrossSection, const double& MCEvents, int period){
  double Weight = 1.0;
  if(IntegratedLuminosityInPb>=IntegratedLuminosityInPbBeforeTriggerChange && IntegratedLuminosityInPb>0){
    double NMCEvents = MCEvents;
    //if(MaxEntry>0)NMCEvents=std::min(MCEvents,(double)MaxEntry);
    if     (period==0)Weight = (CrossSection * (IntegratedLuminosityInPb-IntegratedLuminosityInPbBeforeTriggerChange)) / NMCEvents;
    else if(period==1)Weight = (CrossSection * IntegratedLuminosityInPbBeforeTriggerChange) / NMCEvents;
  }
  return Weight;
}

// compute event weight for the MC background samples based on number of events in the sample, cross section and intergated luminosity
double GetSampleWeightMC(const double& IntegratedLuminosityInPb, const std::vector<string> fileNames, const double& XSection, const double& SampleSize, double MaxEvent, int NSubSamples=1){
  double Weight = 1.0;
   unsigned long InitNumberOfEvents = GetInitialNumberOfMCEvent(fileNames); 
   double SampleEquivalentLumi = InitNumberOfEvents / XSection;
   if(MaxEvent<0)MaxEvent=SampleSize;
   printf("GetSampleWeight MC: IntLumi = %6.2E  SampleLumi = %6.2E --> EventWeight = %6.2E --> ",IntegratedLuminosityInPb,SampleEquivalentLumi, IntegratedLuminosityInPb/SampleEquivalentLumi);
//   printf("Sample NEvent = %6.2E   SampleEventUsed = %6.2E --> Weight Rescale = %6.2E\n",SampleSize, MaxEvent, SampleSize/MaxEvent);
   Weight = (IntegratedLuminosityInPb/SampleEquivalentLumi) * (SampleSize/MaxEvent) * (1.0/NSubSamples);
   printf("FinalWeight = %6.2f\n",Weight);
   return Weight;
}

// compute a weight to make sure that the pileup distribution in Signal/Background MC is compatible to the pileup distribution in data
double GetPUWeight(const fwlite::ChainEvent& ev, const std::string& pileup, double &PUSystFactor, edm::LumiReWeighting& LumiWeightsMC, reweight::PoissonMeanShifter& PShift){
   fwlite::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   PupInfo.getByLabel(ev, "addPileupInfo");
   if(!PupInfo.isValid()){printf("PileupSummaryInfo Collection NotFound\n");return 1.0;}
   double PUWeight_thisevent=1;
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   int npv = -1; float Tnpv = -1;

   if(pileup=="NoPU"){
      PUWeight_thisevent = 1.0;
      PUSystFactor = 1.0;
   }else if(pileup=="S4"){
      float sum_nvtx = 0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         npv = PVI->getPU_NumInteractions();
         sum_nvtx += float(npv);
      }
      float ave_nvtx = sum_nvtx/3.;
      PUWeight_thisevent = LumiWeightsMC.weight( ave_nvtx );
      PUSystFactor = PShift.ShiftWeight( ave_nvtx );
   }else if(pileup=="S3"){
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         int BX = PVI->getBunchCrossing();
         if(BX == 0) {
            npv = PVI->getPU_NumInteractions();
            continue;
         }
      }
      PUWeight_thisevent = LumiWeightsMC.weight( npv );
      PUSystFactor = PShift.ShiftWeight( npv );
   }else if(pileup=="S10"){
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       int BX = PVI->getBunchCrossing();
       if(BX == 0) {
	 Tnpv = PVI->getTrueNumInteractions();
	 continue;
       }
     }
     PUWeight_thisevent = LumiWeightsMC.weight( Tnpv );
     PUSystFactor = PShift.ShiftWeight( Tnpv );   
   }else if(pileup=="S15"){
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       int BX = PVI->getBunchCrossing();
       if(BX == 0) {
         Tnpv = PVI->getTrueNumInteractions();
         continue;
       }
     }
     PUWeight_thisevent = LumiWeightsMC.weight( Tnpv );
     PUSystFactor = PShift.ShiftWeight( Tnpv );
   }else {
     printf("Can not find pile up scenario: %s", pileup.c_str());
   }
   return PUWeight_thisevent;
}



// compute a weight to make sure that the pileup distribution in Signal/Background MC is compatible to the pileup distribution in data
double GetPUWeight(const fwlite::ChainEvent& ev, const std::string& pileup, double &PUSystFactor, edm::LumiReWeighting& LumiWeightsMC, edm::LumiReWeighting& LumiWeightsMCSyst){
   fwlite::Handle<std::vector<PileupSummaryInfo> > PupInfo;
   PupInfo.getByLabel(ev, "addPileupInfo");
   if(!PupInfo.isValid()){printf("PileupSummaryInfo Collection NotFound\n");return 1.0;}
   double PUWeight_thisevent=1;
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   int npv = -1; float Tnpv = -1;

   if(pileup=="NoPU"){
      PUWeight_thisevent = 1.0;
      PUSystFactor = 1.0;
   }else if(pileup=="S4"){
      float sum_nvtx = 0;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         npv = PVI->getPU_NumInteractions();
         sum_nvtx += float(npv);
      }
      float ave_nvtx = sum_nvtx/3.;
      PUWeight_thisevent = LumiWeightsMC.weight( ave_nvtx );
      if(PUWeight_thisevent==0) PUSystFactor=1;
      else PUSystFactor = LumiWeightsMCSyst.weight( ave_nvtx ) / PUWeight_thisevent;
   }else if(pileup=="S3"){
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
         int BX = PVI->getBunchCrossing();
         if(BX == 0) {
            npv = PVI->getPU_NumInteractions();
            continue;
         }
      }
      PUWeight_thisevent = LumiWeightsMC.weight( npv );
      if(PUWeight_thisevent==0) PUSystFactor=1;
      else PUSystFactor = LumiWeightsMCSyst.weight( npv ) / PUWeight_thisevent;
   }else if(pileup=="S10"){
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       int BX = PVI->getBunchCrossing();
       if(BX == 0) {
	 Tnpv = PVI->getTrueNumInteractions();
	 continue;
       }
     }
     PUWeight_thisevent = LumiWeightsMC.weight( Tnpv );
     if(PUWeight_thisevent==0) PUSystFactor=1;
     else PUSystFactor = LumiWeightsMCSyst.weight( Tnpv ) / PUWeight_thisevent;
   }else if(pileup=="S15"){
     for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
       int BX = PVI->getBunchCrossing();
       if(BX == 0) {
         Tnpv = PVI->getTrueNumInteractions();
         continue;
       }
     }
     PUWeight_thisevent = LumiWeightsMC.weight( Tnpv );
     if(PUWeight_thisevent==0) PUSystFactor=1;
     else PUSystFactor = LumiWeightsMCSyst.weight( Tnpv ) / PUWeight_thisevent;
   }else {
     printf("Can not find pile up scenario");
   }
   return PUWeight_thisevent;
}

#endif //end FWLITE block



//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Define Theoretical XSection and error band for 7TeVStops and 7TeVGluino

//NLO Stop@7TeV
const double THXSEC7TeV_Stop_Mass [] = {100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0};
const double THXSEC7TeV_Stop_Cen  [] = {390.0,165.0,76.8,38.8,21.1,12.1,7.21,4.46,2.84,1.86,1.24,0.844,0.585,0.412,0.294,0.212,0.155,0.114,0.0847,0.0634,0.0479,0.0365,0.0279,0.0215,0.0166,0.0129,0.0101,0.00796,0.0063,0.00499,0.00397,0.00317,0.00253,0.00203,0.00164,0.00132,0.00107,0.000864,0.000702,0.000572,0.000467,0.000381,0.000312,0.000255,0.00021,0.000173 };
const double THXSEC7TeV_Stop_Low  [] = {336.687212239,142.161704016,65.8513147385,33.1253989627,17.918668408,10.2276057109,6.06740619265,3.73576507352,2.36659235722,1.54221509094,1.02378918186,0.693358458913,0.478780586501,0.335646071783,0.238501707841,0.171183802392,0.124608032654,0.0912510120893,0.0674833984023,0.0502792202566,0.037816978755,0.0286883112529,0.0218301431028,0.0167475917045,0.012877462018,0.00995875547206,0.00775938400598,0.00608185167256,0.00479071471604,0.00377411064367,0.00298741802917,0.00237365916816,0.00188189555314,0.00150193564851,0.00120631066392,0.000965374807573,0.000778619696638,0.000624245761673,0.000504320503585,0.00040804137369,0.000330700238198,0.000267790633317,0.000217758112762,0.000176730791724,0.000144344874883,0.000117997311653 };
const double THXSEC7TeV_Stop_High [] = {455.337270601,192.484749389,89.6237634573,45.38,24.7104320442,14.2154509779,8.49882260582,5.27246023469,3.3687269963,2.21436345753,1.48113888098,1.01267602748,0.704377399955,0.497929290742,0.356909806888,0.258306503513,0.189623679693,0.14001964183,0.104443008835,0.0785510793723,0.059587810898,0.0455979172605,0.035030746117,0.0271154457368,0.021046423859,0.0164444810585,0.0129238101605,0.0102301624614,0.00812681856605,0.00646997801663,0.00517242916474,0.004153874823,0.00333274574668,0.00268995484779,0.00218320499759,0.00176632187855,0.00143974932708,0.00116964830458,0.000955623561151,0.000783143384745,0.000643472508176,0.000528336426727,0.000435441662818,0.000358165927707,0.000296784169801,0.000246047049368 };


//NLO Stop@8TeV
const double THXSEC8TeV_Stop_Mass [] = {100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0};
const double THXSEC8TeV_Stop_Cen  [] = {526.0,226.0,107.0,55.0,30.3,17.7,10.7,6.69,4.32,2.85,1.93,1.33,0.935,0.666,0.48,0.351,0.259,0.194,0.146,0.111,0.0846,0.0652,0.0505,0.0394,0.0309,0.0244,0.0193,0.0153,0.0123,0.00985,0.00794,0.00642,0.0052,0.00423,0.00345,0.00282,0.00231,0.0019,0.00156,0.00129,0.00107,0.000888,0.000738,0.000614,0.000512,0.000427,0.000357,0.000298,0.00025,0.00021,0.000177,0.000149,0.000126,0.000106,8.95e-05,7.54e-05,6.37e-05,5.4e-05,4.58e-05,3.88e-05,3.29e-05,2.79e-05,2.37e-05,2.01e-05,1.71e-05,1.45e-05,1.24e-05,1.05e-05,8.98e-06,7.65e-06,6.51e-06,5.53e-06,4.7e-06,3.99e-06,3.38e-06,2.86e-06,2.4e-06,2.02e-06,1.73e-06,1.49e-06,1.29e-06,1.1e-06,9.35e-07,7.95e-07,6.77e-07,5.77e-07,4.92e-07,4.19e-07,3.57e-07,3.04e-07,2.58e-07,2.19e-07,1.87e-07,1.59e-07,1.35e-07};
const double THXSEC8TeV_Stop_Low [] = {456.40976717,195.648917903,92.2540349723,47.2865665445,25.9397832106,15.0956844963,9.08927760689,5.65777076356,3.6403946419,2.38881618792,1.61175270749,1.10604627471,0.775030543643,0.549718373166,0.394375877481,0.287159,0.211071455628,0.157614743299,0.118115158583,0.0896125182004,0.067981415327,0.0521458627606,0.0402437661268,0.0312470919086,0.0244182952548,0.019218893828,0.0151215196692,0.0119430303758,0.00956579804952,0.00763312081598,0.00613258287528,0.0049339106769,0.00397887414586,0.00321918065883,0.00261564964342,0.00212811380764,0.00173487125664,0.00142268794902,0.00116204134147,0.000956279372476,0.000788408712819,0.000650528932861,0.000538210968336,0.000445184820841,0.000368817901272,0.000305789996708,0.000254186271492,0.000210958143595,0.000175969661542,0.000146735307698,0.000123013224389,0.000102929880545,8.6402235565e-05,7.222995e-05,6.05164427274e-05,5.0580220783e-05,4.24247029546e-05,3.5715773314e-05,3.00635721872e-05,2.52451336458e-05,2.12193012116e-05,1.78446031904e-05,1.50354723785e-05,1.26276243982e-05,1.06536409697e-05,8.94732163468e-06,7.59026356243e-06,6.35376571526e-06,5.37945064143e-06,4.53447994822e-06,3.81776527832e-06,3.21539615667e-06,2.71447659044e-06,2.28387219518e-06,1.91493168116e-06,1.60622722428e-06,1.34502348655e-06,1.1253862558e-06,9.53375342241e-07,8.07186767383e-07,6.85400212545e-07,5.77264615054e-07,4.84520464268e-07,4.06967984563e-07,3.41287787427e-07,2.86338664102e-07,2.40809064183e-07,2.02382848192e-07,1.70268348502e-07,1.4267677199e-07,1.19072781243e-07,9.95434872399e-08,8.37672797081e-08,7.02461674403e-08,5.8546055207e-08};
const double THXSEC8TeV_Stop_High  [] = {612.151729265,262.612588103,124.288464259,63.9150906166,35.2825924691,20.6572041873,12.51422531,7.84901646422,5.08233985546,3.36513880209,2.28626927651,1.57937,1.11448052719,0.796049681797,0.575829786962,0.422291028968,0.312903555384,0.235042620074,0.177634138358,0.135439676246,0.103594494844,0.0801823918578,0.0623096661467,0.0488132526053,0.0384157751721,0.030457195299,0.02419270263,0.0192543023786,0.0155348683796,0.0124867180798,0.0101081955359,0.00820720230018,0.00667610315198,0.00545393409975,0.00446719444332,0.00366690218399,0.0030168136425,0.00249290452531,0.00205643855325,0.00170910492719,0.0014249051169,0.00118823222984,0.000992225516805,0.00082975230322,0.000695336835754,0.000582782847848,0.000489743672799,0.000410791069838,0.000346582593019,0.000292517320796,0.000247919343036,0.00020993472894,0.000178644320362,0.000151183052034,0.000128497606094,0.000108904879045,9.25617106816e-05,7.88942071601e-05,6.73641928454e-05,5.74235824351e-05,4.90287275026e-05,4.18712011298e-05,3.58244154887e-05,3.05926261591e-05,2.62141570875e-05,2.23830855768e-05,1.92778016416e-05,1.64463855993e-05,1.41720928606e-05,1.21588956729e-05,1.04173416472e-05,8.9066137467e-06,7.61913263458e-06,6.50912689733e-06,5.54723440685e-06,4.71830691087e-06,3.96587175902e-06,3.34925029375e-06,2.88869395402e-06,2.51958334975e-06,2.20968359995e-06,1.90103535042e-06,1.62823072952e-06,1.3949661943e-06,1.1975604273e-06,1.02928057305e-06,8.84816255638e-07,7.59765184201e-07,6.52949889306e-07,5.6101392013e-07,4.80393896187e-07,4.11598677578e-07,3.54793913074e-07,3.04265637883e-07,2.60687256456e-07};

//NLO Stop@13TeV
const double THXSEC13TeV_Stop_Mass [] = {100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0,2000.0,2020.0,2040.0,2060.0,2080.0,2100.0,2120.0,2140.0,2160.0,2180.0,2200.0,2220.0,2240.0,2260.0,2280.0,2300.0,2320.0,2340.0,2360.0,2380.0,2400.0,2420.0,2440.0,2460.0,2480.0,2500.0,2520.0,2540.0,2560.0,2580.0,2600.0,2620.0,2640.0,2660.0,2680.0,2700.0,2720.0,2740.0,2760.0,2780.0,2800.0,2820.0,2840.0,2860.0,2880.0,2900.0,2920.0,2940.0,2960.0,2980.0,3000.0};
const double THXSEC13TeV_Stop_Cen  [] = {1430.0,651.0,324.0,175.0,101.0,61.2,38.6,25.2,16.9,11.6,8.16,5.84,4.24,3.13,2.34,1.77,1.35,1.04,0.813,0.638,0.505,0.403,0.323,0.261,0.211,0.172,0.141,0.117,0.0964,0.0799,0.0664,0.0555,0.0466,0.0392,0.0331,0.028,0.0238,0.0202,0.0173,0.0148,0.0127,0.0109,0.00937,0.00807,0.00698,0.00604,0.00524,0.00455,0.00396,0.00345,0.00301,0.00263,0.00231,0.00202,0.00178,0.00156,0.00137,0.0012,0.00106,0.000937,0.000829,0.000733,0.000649,0.000575,0.00051,0.000452,0.000402,0.000357,0.000318,0.000282,0.000251,0.000223,0.000198,0.000176,0.000157,0.000139,0.000122,0.000107,9.6e-05,8.74e-05,7.96e-05,7.15e-05,6.4e-05,5.73e-05,5.14e-05,4.62e-05,4.15e-05,3.72e-05,3.35e-05,3.01e-05,2.71e-05,2.44e-05,2.19e-05,1.98e-05,1.78e-05,1.61e-05,1.45e-05,1.31e-05,1.18e-05,1.06e-05,9.62e-06,8.68e-06,7.83e-06,7.07e-06,6.39e-06,5.78e-06,5.22e-06,4.72e-06,4.27e-06,3.86e-06,3.49e-06,3.16e-06,2.86e-06,2.59e-06,2.34e-06,2.12e-06,1.92e-06,1.73e-06,1.57e-06,1.42e-06,1.29e-06,1.17e-06,1.05e-06,9.52e-07,8.63e-07,7.83e-07,7.09e-07,6.42e-07,5.81e-07,5.27e-07,4.77e-07,4.32e-07,3.91e-07,3.54e-07,3.2e-07,2.9e-07,2.62e-07,2.38e-07,2.15e-07,1.95e-07,1.77e-07,1.6e-07,1.45e-07,1.31e-07,1.19e-07,1.07e-07};
const double THXSEC13TeV_Stop_Low [] = {1256.5610395,571.261699469,284.062061983,153.302765158,88.5282049222,53.5883052288,33.705633428,21.982706154,14.6802081348,10.0481072625,7.0423554603,5.03255599779,3.64942918745,2.6858808265,2.0051985316,1.51209615168,1.14928564212,0.882418832497,0.688432011719,0.539252043663,0.426256201508,0.339578239856,0.271101052868,0.218642,0.176076469231,0.143063266297,0.117151068195,0.0971513412524,0.0798747177311,0.0659900866275,0.0546894596644,0.0456398104122,0.0382056586819,0.0320846830923,0.0270266853734,0.0228362505255,0.0193784352081,0.0163950359271,0.0140260565975,0.011954394223,0.010240120217,0.00876916472914,0.00752292464905,0.00646459672699,0.00557616477535,0.00481262413365,0.00416583956047,0.00360850010833,0.00313327357196,0.00272463149596,0.00237057006467,0.00206598998657,0.00181081503247,0.00157940743077,0.00138671148257,0.00121253055697,0.00106273613282,0.000927513309843,0.000817564771167,0.000720855458188,0.000636274788562,0.000560893691838,0.000494768520734,0.0004371666394,0.000386941446323,0.000342190557844,0.000303391930989,0.000268653105042,0.000238555208375,0.000210708310071,0.000187140862063,0.000166027760726,0.000147082530071,0.000130642243789,0.000116325338671,0.000103022269472,9.07878863427e-05,7.96969057042e-05,7.13411629362e-05,6.43598012648e-05,5.8117102697e-05,5.2006967459e-05,4.64648201519e-05,4.14888484066e-05,3.70669255445e-05,3.31761047918e-05,2.97396724727e-05,2.65578946232e-05,2.38700231066e-05,2.13811798886e-05,1.91773243071e-05,1.71952902967e-05,1.53691449581e-05,1.38376924142e-05,1.23800400704e-05,1.11470375047e-05,9.9901424048e-06,8.98781129874e-06,8.0604046832e-06,7.213337766e-06,6.52304022367e-06,5.85446630869e-06,5.25297544971e-06,4.72059330312e-06,4.24935081511e-06,3.82967608264e-06,3.4394174603e-06,3.09189737931e-06,2.78444279208e-06,2.50721546849e-06,2.25700967378e-06,2.03323130265e-06,1.82981539341e-06,1.64583242741e-06,1.48185138054e-06,1.33459447985e-06,1.20180420197e-06,1.07557225599e-06,9.69080738751e-07,8.72879502996e-07,7.88316804623e-07,7.11166174966e-07,6.33559127448e-07,5.70159078997e-07,5.13961464984e-07,4.63103042588e-07,4.16503566364e-07,3.75003279988e-07,3.36952929132e-07,3.03409902741e-07,2.72569568096e-07,2.45121334016e-07,2.20534771127e-07,1.98179581858e-07,1.77701704085e-07,1.59783974849e-07,1.43245273872e-07,1.29184259826e-07,1.15779971241e-07,1.04219701573e-07,9.38234888415e-08,8.40711501724e-08,7.55416381428e-08,6.76830159402e-08,6.0995990023e-08,5.43379171514e-08};
const double THXSEC13TeV_Stop_High  [] = {1655.4389605,749.659041405,371.747993042,200.408768416,115.671795078,70.1689992208,44.3130221333,28.9533186014,19.4504817919,13.3602363129,9.40189817771,6.74056929826,4.90423236034,3.62789120375,2.71900623471,2.06131651419,1.57441120363,1.21506736368,0.951055510673,0.74807684546,0.592882997822,0.47435432439,0.380973711128,0.308661358047,0.250067412312,0.20412935419,0.167678554169,0.139298003165,0.115091727636,0.0956080098584,0.0797001026777,0.0667537444146,0.0561350837241,0.0473143941983,0.0400177174555,0.0339217004101,0.0289205051724,0.0246120947729,0.0211334312799,0.0181141215504,0.0155689066196,0.0133939923644,0.0115415594333,0.00997132034622,0.00864386920142,0.00749968021293,0.00652039883583,0.00567199817812,0.00494678469058,0.00431841160476,0.00377771433589,0.00331087299637,0.00291435632145,0.00255663153112,0.00225736026771,0.00198324816028,0.00174554065756,0.00153341394573,0.00135852984494,0.00120460899617,0.00106911589897,0.000947943511605,0.000841066696334,0.000747438147953,0.000664207422307,0.000590541937873,0.000526857998883,0.00046975,0.000419610176018,0.000373123735911,0.00033278307296,0.000296263672267,0.00026342449443,0.000234454822267,0.000209372566885,0.000185336120032,0.000161917843469,0.000141652217832,0.000127453589089,0.000116951011678,0.000107519746269,9.69923818685e-05,8.70711329015e-05,7.81224374231e-05,7.02452362347e-05,6.32850459313e-05,5.70397800658e-05,5.12746886223e-05,4.63618226483e-05,4.1808861403e-05,3.77713857076e-05,3.41063492412e-05,3.07159457849e-05,2.78539908055e-05,2.51240462964e-05,2.27970791942e-05,2.06067576482e-05,1.86840746819e-05,1.68975546131e-05,1.52276722846e-05,1.38617682601e-05,1.25635467257e-05,1.13775015929e-05,1.03152070253e-05,9.35242737307e-06,8.49131704165e-06,7.70257722403e-06,6.99553118258e-06,6.35480160239e-06,5.77040051692e-06,5.23731075256e-06,4.76240870378e-06,4.32869232274e-06,3.93631002171e-06,3.57186680479e-06,3.25010581909e-06,2.95805618218e-06,2.67896548332e-06,2.44106448446e-06,2.21720513453e-06,2.02259030218e-06,1.84309361598e-06,1.66240700347e-06,1.51544883703e-06,1.37996527272e-06,1.25750631895e-06,1.14380637097e-06,1.04096434124e-06,9.46902488849e-07,8.63189855035e-07,7.85226645799e-07,7.14815762798e-07,6.5015815103e-07,5.91632879557e-07,5.37749483941e-07,4.89656752483e-07,4.44646214417e-07,4.0600225538e-07,3.68634396903e-07,3.36103041377e-07,3.06707385218e-07,2.78758127166e-07,2.53939392793e-07,2.30782842334e-07,2.10714626111e-07,1.9049299126e-07};





//NLO Gluino@7TeV
const double THXSEC7TeV_Gluino_Mass [] = {300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0 };
const double THXSEC7TeV_Gluino_Cen  [] = {65.8,44.8,31.0,21.7,15.5,11.2,8.15,6.01,4.47,3.36,2.54,1.93,1.48,1.14,0.888,0.693,0.543,0.428,0.339,0.269,0.214,0.172,0.137,0.111,0.0895,0.0725,0.0588,0.0479,0.0391,0.032,0.0262,0.0215,0.0176,0.0145,0.0119,0.00987,0.00815,0.00674,0.00559,0.00464,0.00386,0.0032,0.00266,0.00222,0.00185,0.00154,0.00128,0.00107,0.000893,0.000746,0.000623,0.00052,0.000435,0.000363,0.000304,0.000254,0.000212,0.000177,0.000148,0.000124,0.000103,8.64e-05,7.22e-05,6.02e-05,5.02e-05,4.18e-05,3.49e-05,2.9e-05,2.42e-05,2.02e-05,1.67e-05,1.39e-05,1.16e-05,9.59e-06,7.95e-06,6.59e-06,5.45e-06,4.51e-06,3.73e-06,3.08e-06,2.54e-06,2.09e-06,1.72e-06,1.42e-06,1.16e-06 };
const double THXSEC7TeV_Gluino_Low  [] = {53.9915787502,36.3679605186,25.04773106,17.5484531801,12.5029969909,8.9137230832,6.48814044136,4.71554753537,3.51079991399,2.63456118649,1.96607018967,1.50525499617,1.14834208312,0.886682415379,0.67892705475,0.528365261709,0.415476340911,0.322928423644,0.255343774802,0.202059219533,0.161164211822,0.127363153734,0.102023478594,0.0818182928751,0.065393518095,0.0529567854047,0.0430807023909,0.0345084584098,0.0280914623477,0.0229311955562,0.0184911518907,0.0152076608485,0.0124075251946,0.0102192110979,0.00832945628418,0.00686771000809,0.00568217384168,0.00462615538234,0.00382900127407,0.00313534410297,0.00259805136043,0.00215484050605,0.00180048617066,0.001467341696,0.00122061808,0.00100694617882,0.000838421855019,0.000684769467067,0.000571644050552,0.000476594263719,0.000392527368983,0.000326777849928,0.000268212922309,0.000224036696952,0.000184330609415,0.0001530784447,0.000126064184114,0.000105014721299,8.62123946791e-05,7.21482688124e-05,5.88983434237e-05,4.93030990383e-05,4.04422743416e-05,3.31101477016e-05,2.75456714875e-05,2.25762430633e-05,1.87694070612e-05,1.53040540155e-05,1.24350948334e-05,1.04092345492e-05,8.49269279945e-06,6.97519707277e-06,5.71827736139e-06,4.62718381962e-06,3.82358391349e-06,3.10935813704e-06,2.50561363539e-06,2.0258883204e-06,1.67348322755e-06,1.34494877209e-06,1.08865400522e-06,8.98317665716e-07,7.1490062492e-07,5.71185862116e-07,4.67630704264e-07 };
const double THXSEC7TeV_Gluino_High [] = {79.5668288977,54.2211352869,37.5922538006,26.5343616416,19.075,13.786,10.1121682798,7.45102358218,5.57973835763,4.22381192237,3.19868564961,2.44746309321,1.88138665258,1.45913982451,1.13478647663,0.892932349164,0.704893627141,0.556473798188,0.440065907012,0.352253724644,0.282471833163,0.226566443289,0.182846908818,0.147392902942,0.119981370949,0.0972614292036,0.0794933605454,0.0647760279165,0.0531432044804,0.0435131063759,0.0359897509829,0.0297064480657,0.0243363202323,0.0201945890666,0.0166888696477,0.013839213608,0.0115485526692,0.00955917250779,0.00797993717028,0.00667993235001,0.00555403311772,0.00464008175518,0.00388720169471,0.00326472742026,0.00273449708667,0.00230578041856,0.00191402087824,0.00161285546003,0.00135345358155,0.00113904856138,0.000958125352577,0.000805385272244,0.000678241538475,0.000570293065787,0.000481037075751,0.000404734213278,0.000340132555915,0.000285662391546,0.000240024886815,0.000203182671651,0.000169737973538,0.000143208791006,0.000120394572153,0.000101748910812,8.53401846951e-05,7.16178115223e-05,6.0214226498e-05,5.05737327547e-05,4.23835072134e-05,3.58341122051e-05,2.98807306886e-05,2.50681230531e-05,2.11198589408e-05,1.76419321064e-05,1.47213338727e-05,1.23406750758e-05,1.02636214127e-05,8.59048670095e-06,7.18944995696e-06,5.96042280205e-06,4.97224153032e-06,4.14424242035e-06,3.4282369471e-06,2.85290331058e-06,2.35628618585e-06 };


//NLO Gluino@8TeV
const double THXSEC8TeV_Gluino_Mass [] = {300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0};
const double THXSEC8TeV_Gluino_Cen  [] = {103.0,70.6,49.3,35.1,25.3,18.5,13.6,10.2,7.66,5.82,4.46,3.44,2.66,2.08,1.63,1.29,1.02,0.816,0.653,0.525,0.424,0.344,0.279,0.228,0.186,0.153,0.126,0.104,0.0856,0.0709,0.0588,0.0489,0.0408,0.0341,0.0285,0.0239,0.02,0.0168,0.0141,0.0119,0.01,0.00845,0.00714,0.00605,0.00513,0.00435,0.00369,0.00314,0.00267,0.00227,0.00193,0.00164,0.0014,0.00119,0.00102,0.000866,0.000739,0.000631,0.000539,0.00046,0.000393,0.000336,0.000287,0.000245,0.000209,0.000179,0.000153,0.000131,0.000112,9.52e-05,8.13e-05,6.94e-05,5.93e-05,5.06e-05,4.31e-05,3.68e-05,3.14e-05,2.67e-05,2.28e-05,1.94e-05,1.65e-05,1.4e-05,1.19e-05,1.02e-05,8.63e-06};
const double THXSEC8TeV_Gluino_Low  [] = {85.4522933503,58.3024212534,40.5404072161,28.7311959371,20.6110330718,14.9981486143,10.9817116552,8.2045980396,6.13599856717,4.64639599616,3.54679219469,2.7264626978,2.09915349508,1.63562861604,1.27699308705,1.00685323305,0.792356354919,0.632441971506,0.503862668979,0.403411581822,0.324940255218,0.262909295502,0.212319587364,0.173005424976,0.140570743272,0.115316873907,0.0947326706086,0.0778816775066,0.0638144532538,0.0526220125758,0.0434789865444,0.0359951306958,0.029944567911,0.0249059836578,0.020754521583,0.0173207451689,0.0144267991917,0.0120799394065,0.0100863120603,0.00847183872318,0.00708717142879,0.00595886341748,0.00501035858415,0.00422415500896,0.00356494485684,0.00300727877119,0.00253860291188,0.00214977087027,0.00181901909611,0.00153605196417,0.00129944484705,0.00109821151599,0.000931476154739,0.000787499214937,0.000671494422707,0.000566354149908,0.000479894886148,0.000407420453834,0.000345671282251,0.000292733637163,0.000248262306276,0.000210814299455,0.000178718470665,0.000151342195875,0.000127869581852,0.000108804207517,9.23062817067e-05,7.82498615939e-05,6.63400478633e-05,5.58615231572e-05,4.72597016726e-05,3.99537211035e-05,3.38756321908e-05,2.86208917333e-05,2.41297344895e-05,2.0403820875e-05,1.72341338273e-05,1.44941525431e-05,1.22590511119e-05,1.03016340378e-05,8.66483254755e-06,7.26134469144e-06,6.10042739754e-06,5.16703392541e-06,4.32379602699e-06};
const double THXSEC8TeV_Gluino_High [] = {123.042915069,84.6204861909,59.3621530586,42.4270342231,30.720245281,22.5429742994,16.6400543991,12.521,9.44012836998,7.19888313862,5.53783839898,4.28718343406,3.32951446035,2.61303391779,2.0546626767,1.63189165427,1.29573309725,1.03959612686,0.835270344823,0.673688637441,0.546222752456,0.444771492296,0.362012997158,0.297010039438,0.243115317971,0.200703894146,0.165936031728,0.137421289844,0.113579295087,0.094359434471,0.0785768709139,0.0656067780802,0.054914187159,0.0460700308168,0.0386624800886,0.0325504144943,0.0273588316481,0.0230737980716,0.019440405649,0.0164822382849,0.0139146663585,0.0118012857472,0.0100182943312,0.00852079529061,0.00725778032384,0.00618194507243,0.00526737821502,0.00450197140212,0.00384501422077,0.00328618755105,0.00280677982274,0.00239852281744,0.00205691112089,0.00175754634874,0.00151421701532,0.00129260450299,0.00110862278184,0.000951522205324,0.000817402023602,0.000701118190558,0.000602322798422,0.000517825498222,0.000444861065588,0.000382177349458,0.000327902228753,0.00028253307733,0.000243007105712,0.00020936089293,0.000180051621273,0.000154016576997,0.000132434560814,0.000113806344946,9.78639305045e-05,8.40542824159e-05,7.21176753055e-05,6.20064018569e-05,5.32790981711e-05,4.56654057754e-05,3.92622063114e-05,3.3654042914e-05,2.88430383571e-05,2.46542444561e-05,2.11238854891e-05,1.82236973231e-05,1.55451905173e-05};

//NLO Gluino@13TeV
const double THXSEC13TeV_Gluino_Mass [] = {100.0,120.0,140.0,160.0,180.0,200.0,220.0,240.0,260.0,280.0,300.0,320.0,340.0,360.0,380.0,400.0,420.0,440.0,460.0,480.0,500.0,520.0,540.0,560.0,580.0,600.0,620.0,640.0,660.0,680.0,700.0,720.0,740.0,760.0,780.0,800.0,820.0,840.0,860.0,880.0,900.0,920.0,940.0,960.0,980.0,1000.0,1020.0,1040.0,1060.0,1080.0,1100.0,1120.0,1140.0,1160.0,1180.0,1200.0,1220.0,1240.0,1260.0,1280.0,1300.0,1320.0,1340.0,1360.0,1380.0,1400.0,1420.0,1440.0,1460.0,1480.0,1500.0,1520.0,1540.0,1560.0,1580.0,1600.0,1620.0,1640.0,1660.0,1680.0,1700.0,1720.0,1740.0,1760.0,1780.0,1800.0,1820.0,1840.0,1860.0,1880.0,1900.0,1920.0,1940.0,1960.0,1980.0,2000.0,2020.0,2040.0,2060.0,2080.0,2100.0,2120.0,2140.0,2160.0,2180.0,2200.0,2220.0,2240.0,2260.0,2280.0,2300.0,2320.0,2340.0,2360.0,2380.0,2400.0,2420.0,2440.0,2460.0,2480.0,2500.0,2520.0,2540.0,2560.0,2580.0,2600.0,2620.0,2640.0,2660.0,2680.0,2700.0,2720.0,2740.0,2760.0,2780.0,2800.0,2820.0,2840.0,2860.0,2880.0,2900.0,2920.0,2940.0,2960.0,2980.0,3000.0};
const double THXSEC13TeV_Gluino_Cen  [] = {81700.0,38300.0,19200.0,10200.0,5740.0,3400.0,2110.0,1360.0,911.0,628.0,444.0,316.0,228.0,168.0,125.0,94.8,72.2,55.6,43.2,33.9,26.8,21.3,17.1,13.7,11.1,9.07,7.43,6.11,5.05,4.19,3.49,2.92,2.45,2.06,1.74,1.47,1.25,1.06,0.908,0.777,0.667,0.574,0.494,0.427,0.369,0.32,0.278,0.241,0.21,0.183,0.16,0.14,0.123,0.108,0.0949,0.0836,0.0737,0.065,0.0574,0.0507,0.0449,0.0398,0.0353,0.0313,0.0278,0.0247,0.022,0.0196,0.0175,0.0156,0.0139,0.0124,0.0111,0.00993,0.00889,0.00796,0.00713,0.0064,0.00574,0.00515,0.00463,0.00416,0.00374,0.00337,0.00303,0.00273,0.00246,0.00221,0.00199,0.0018,0.00162,0.00146,0.00132,0.00119,0.00108,0.000974,0.00088,0.000796,0.00072,0.000651,0.000589,0.000533,0.000482,0.000437,0.000395,0.000358,0.000324,0.000294,0.000266,0.000242,0.000219,0.000198,0.00018,0.000163,0.000148,0.000134,0.000121,0.00011,9.99e-05,9.05e-05,8.21e-05,7.44e-05,6.75e-05,6.12e-05,5.55e-05,5.03e-05,4.56e-05,4.14e-05,3.76e-05,3.41e-05,3.09e-05,2.8e-05,2.54e-05,2.3e-05,2.09e-05,1.89e-05,1.71e-05,1.56e-05,1.41e-05,1.28e-05,1.16e-05,1.05e-05,9.52e-06,8.62e-06,7.79e-06,7.04e-06};
const double THXSEC13TeV_Gluino_Low  [] = {72161.7063696,33837.5930717,16940.3444138,8983.2,5043.22871595,2977.20017271,1843.19484067,1183.40495545,789.384113457,542.307938127,381.915741715,270.902756156,194.943429131,143.305366433,106.196546837,80.4245811605,61.0334334009,46.8543434732,36.2773506503,28.3804293937,22.3782678927,17.730707219,14.205642051,11.3349673078,9.16241700515,7.47329203341,6.10744507272,5.0069653823,4.12858395283,3.41538829416,2.8374698137,2.37011075404,1.98380628989,1.66492968592,1.40306820284,1.18249991944,1.00301789065,0.849391956206,0.726099321247,0.620241134607,0.531075412297,0.455912517442,0.391079662199,0.337143503437,0.290891662193,0.251526297158,0.218235472419,0.188658002585,0.164238466538,0.142956567697,0.124800952403,0.108992444786,0.0955680495706,0.08375125149,0.0733934409601,0.064539664751,0.0568307653218,0.050000140099,0.0440998202223,0.0388565359835,0.034370216963,0.0304226612132,0.026902617145,0.0238061228493,0.0210854855025,0.0186935941961,0.0166242363827,0.0147837193785,0.0131811847709,0.0117164319207,0.0104204595025,0.00926054375268,0.00827779612419,0.0073902148408,0.00659652371963,0.00589015748147,0.00526827751583,0.00471589380334,0.00422122608767,0.0037765380786,0.00339029387912,0.00303726986541,0.00272214029258,0.00244896532521,0.0021951656103,0.00197161559956,0.00177045952293,0.00158672392247,0.00142392327118,0.0012859961877,0.00115254696007,0.00103538910071,0.000932937468043,0.000838746901827,0.000759405380401,0.000682132899996,0.000614302294493,0.000553877292717,0.00049982958146,0.000450446945555,0.00040615506283,0.000366192271946,0.000330020454143,0.000298344533118,0.000268303300034,0.000242351419144,0.000218668853146,0.000197914565634,0.000178198197728,0.000161723672562,0.000145856084349,0.000131366629909,0.000119047595613,0.000107420159415,9.71200269058e-05,8.76349672194e-05,7.87303163987e-05,7.13071401279e-05,6.44219161384e-05,5.81143706392e-05,5.24519797316e-05,4.73494658574e-05,4.27175063803e-05,3.85741644712e-05,3.4789379876e-05,3.14063745846e-05,2.83135685919e-05,2.56229482196e-05,2.31532770394e-05,2.09246959517e-05,1.88584261303e-05,1.70033234031e-05,1.53430949808e-05,1.38487298508e-05,1.25258658332e-05,1.12614454525e-05,1.01237097909e-05,9.1925206745e-06,8.26450453623e-06,7.46221409952e-06,6.71392516207e-06,6.03423810205e-06,5.44094198926e-06,4.89520648906e-06,4.39345902644e-06,3.94648018213e-06};
const double THXSEC13TeV_Gluino_High [] = {93344.2921956,43735.3294037,21913.7364067,11647.3495234,6556.99110339,3889.62460578,2417.9975287,1561.11753821,1047.53748855,723.25261799,512.583645667,365.690363525,264.611427602,195.410774817,145.838834765,110.882137406,84.6831995073,65.3155121257,50.8905039656,40.0223995803,31.7077221975,25.2795792728,20.3341546841,16.3424215238,13.2832967783,10.8762977262,8.92530250464,7.35624063003,6.09718564836,5.0686321648,4.23150949152,3.54753753207,2.98431349097,2.51583238461,2.12986506427,1.80460040321,1.53734744882,1.30768928941,1.12238390423,0.962034340895,0.82784614874,0.713401675895,0.615326292428,0.532736390486,0.461097175682,0.400828673991,0.349044973575,0.303459005222,0.265322992836,0.231854444034,0.203194054375,0.1782,0.156789157379,0.138027104033,0.121468081308,0.107238060707,0.094782227773,0.0837859421338,0.0741260166163,0.0656308768985,0.0582683956675,0.0517281772427,0.045990218397,0.0408434802558,0.0363650509399,0.0323627220477,0.0288888838493,0.0258131269014,0.0230964073626,0.0206485474245,0.0184421076064,0.0164817730536,0.0147873789048,0.0132522768791,0.0118825398021,0.0106667803967,0.0095763684134,0.00862144249975,0.0077495984105,0.00697553779103,0.00628579188165,0.00566109784231,0.00510225048898,0.00460854769233,0.0041540534974,0.00375249608915,0.00338892281926,0.00305307146342,0.00275582567926,0.00249881758212,0.00225514624647,0.00203696704864,0.00184750853474,0.00166952479669,0.00151901941197,0.00137360877083,0.00124476951123,0.00112873074795,0.00102424970704,0.000928254539262,0.000842482592478,0.000764587141632,0.000693198822557,0.000630261494851,0.000571506111076,0.000519474786663,0.00047169398167,0.000429200127609,0.000389558656885,0.000355688898633,0.000322884360291,0.000292896639198,0.000267150702403,0.000242648437319,0.000221095138235,0.000200789877818,0.000181867271721,0.0001659246771,0.00015106666297,0.000137303221447,0.000124989618821,0.000113705121407,0.000103553136886,9.43093583596e-05,8.58359430324e-05,7.8137939507e-05,7.10869098301e-05,6.48051499649e-05,5.90498489808e-05,5.37685213751e-05,4.88959924378e-05,4.44909508276e-05,4.05007182877e-05,3.68238586057e-05,3.35860595122e-05,3.05027589972e-05,2.77064010426e-05,2.53466050824e-05,2.30053023926e-05,2.09604668494e-05,1.90727887348e-05,1.73340905265e-05,1.57863435767e-05,1.43493148317e-05,1.30364353443e-05,1.18365175937e-05};


//LO GMSB Stau@7TeV
const double THXSEC7TeV_GMStau_Mass [] = {100,126,156,200,247,308,370,432,494, 557};
const double THXSEC7TeV_GMStau_Cen  [] = {1.85875, 0.437087, 0.124238, 0.0271784, 0.00848349, 0.00238113, 0.000797847, 0.000300147, 0.000123075, 5.39018e-05};
const double THXSEC7TeV_GMStau_CenO  [] = {1.3398,0.274591,0.0645953,0.0118093,0.00342512,0.00098447,0.000353388,0.000141817,6.17749e-05, -1};
const double THXSEC7TeV_GMStau_Low  [] = {1.18163*THXSEC7TeV_GMStau_Cen[0]/THXSEC7TeV_GMStau_CenO[0],0.242982*THXSEC7TeV_GMStau_Cen[1]/THXSEC7TeV_GMStau_CenO[1],0.0581651*THXSEC7TeV_GMStau_Cen[2]/THXSEC7TeV_GMStau_CenO[2],0.0109992*THXSEC7TeV_GMStau_Cen[3]/THXSEC7TeV_GMStau_CenO[3],0.00324853*THXSEC7TeV_GMStau_Cen[4]/THXSEC7TeV_GMStau_CenO[4],0.00093519*THXSEC7TeV_GMStau_Cen[5]/THXSEC7TeV_GMStau_CenO[5],0.000335826*THXSEC7TeV_GMStau_Cen[6]/THXSEC7TeV_GMStau_CenO[6],0.000134024*THXSEC7TeV_GMStau_Cen[7]/THXSEC7TeV_GMStau_CenO[7],5.83501e-05*THXSEC7TeV_GMStau_Cen[8]/THXSEC7TeV_GMStau_CenO[8], 5.83501e-05*THXSEC7TeV_GMStau_Cen[9]/THXSEC7TeV_GMStau_CenO[8]};
const double THXSEC7TeV_GMStau_High [] = {1.48684*THXSEC7TeV_GMStau_Cen[0]/THXSEC7TeV_GMStau_CenO[0],0.304386*THXSEC7TeV_GMStau_Cen[1]/THXSEC7TeV_GMStau_CenO[1],0.0709262*THXSEC7TeV_GMStau_Cen[2]/THXSEC7TeV_GMStau_CenO[2],0.012632*THXSEC7TeV_GMStau_Cen[3]/THXSEC7TeV_GMStau_CenO[3],0.00358232*THXSEC7TeV_GMStau_Cen[4]/THXSEC7TeV_GMStau_CenO[4],0.00102099*THXSEC7TeV_GMStau_Cen[5]/THXSEC7TeV_GMStau_CenO[5],0.000366819*THXSEC7TeV_GMStau_Cen[6]/THXSEC7TeV_GMStau_CenO[6],0.000147665*THXSEC7TeV_GMStau_Cen[7]/THXSEC7TeV_GMStau_CenO[7],6.45963e-05*THXSEC7TeV_GMStau_Cen[8]/THXSEC7TeV_GMStau_CenO[8], 6.45963e-05*THXSEC7TeV_GMStau_Cen[9]/THXSEC7TeV_GMStau_CenO[8]};

const double THXSEC8TeV_GMStau_Mass [] = {100, 126, 156, 200, 247, 308, 370, 432, 494, 557 };
const double THXSEC8TeV_GMStau_Cen  [] = {2.95785, 0.715665, 0.204696, 0.0446997, 0.01391,  0.00400497, 0.0014057,  0.000555594,  0.000238697, 0.000109033};
const double THXSEC8TeV_GMStau_CenO [] = {2.2834, 0.495872, 0.120584, 0.0214979, 0.00575159, 0.00156062, 0.00056426, 0.000234305, 0.000105751, 5.07802e-05};
const double THXSEC8TeV_GMStau_Low  [] = {2.02025*THXSEC8TeV_GMStau_Cen[0]/THXSEC8TeV_GMStau_CenO[0], 0.437789*THXSEC8TeV_GMStau_Cen[1]/THXSEC8TeV_GMStau_CenO[1], 0.107208*THXSEC8TeV_GMStau_Cen[2]/THXSEC8TeV_GMStau_CenO[2], 0.019625*THXSEC8TeV_GMStau_Cen[3]/THXSEC8TeV_GMStau_CenO[3], 0.00540193*THXSEC8TeV_GMStau_Cen[4]/THXSEC8TeV_GMStau_CenO[4], 0.00148344*THXSEC8TeV_GMStau_Cen[5]/THXSEC8TeV_GMStau_CenO[5], 0.000535717*THXSEC8TeV_GMStau_Cen[6]/THXSEC8TeV_GMStau_CenO[6], 0.000222165*THXSEC8TeV_GMStau_Cen[7]/THXSEC8TeV_GMStau_CenO[7], 9.99954e-05*THXSEC8TeV_GMStau_Cen[8]/THXSEC8TeV_GMStau_CenO[8], 4.79772e-05*THXSEC8TeV_GMStau_Cen[9]/THXSEC8TeV_GMStau_CenO[9]};
const double THXSEC8TeV_GMStau_High [] = {2.52494*THXSEC8TeV_GMStau_Cen[0]/THXSEC8TeV_GMStau_CenO[0], 0.548973*THXSEC8TeV_GMStau_Cen[1]/THXSEC8TeV_GMStau_CenO[1], 0.133042*THXSEC8TeV_GMStau_Cen[2]/THXSEC8TeV_GMStau_CenO[2], 0.0233234*THXSEC8TeV_GMStau_Cen[3]/THXSEC8TeV_GMStau_CenO[3], 0.00609887*THXSEC8TeV_GMStau_Cen[4]/THXSEC8TeV_GMStau_CenO[4], 0.00162316*THXSEC8TeV_GMStau_Cen[5]/THXSEC8TeV_GMStau_CenO[5], 0.000585186*THXSEC8TeV_GMStau_Cen[6]/THXSEC8TeV_GMStau_CenO[6], 0.000243149*THXSEC8TeV_GMStau_Cen[7]/THXSEC8TeV_GMStau_CenO[7], 0.000109999*THXSEC8TeV_GMStau_Cen[8]/THXSEC8TeV_GMStau_CenO[8], 5.30139e-05*THXSEC8TeV_GMStau_Cen[9]/THXSEC8TeV_GMStau_CenO[9]};

//FIXME error must be computed properly, currently it's just Cen+-5%
const double THXSEC13TeV_GMStau_Mass [] = {100, 126, 156, 200, 247, 308, 370, 432, 494, 557, 619, 651, 745, 871, 1029, 1218, 1409, 1599};
const double THXSEC13TeV_GMStau_Cen  [] = { 14.2831, 3.90949, 1.20093, 0.283988, 0.0882234, 0.0249859, 0.00907959, 0.00390754, 0.0018932, 0.000990061, 0.000545049, 0.000409813, 0.000184222, 6.89994e-05, 2.24488e-05, 6.41171e-06, 2.01349e-06, 5.26691e-07};
const double THXSEC13TeV_GMStau_Low  [] = {0.95*14.2831, 0.95*3.90949, 0.95*1.20093, 0.95*0.283988, 0.95*0.0882234, 0.95*0.0249859, 0.95*0.00907959, 0.95*0.00390754, 0.95*0.0018932, 0.95*0.000990061, 0.95*0.000545049, 0.95*0.000409813, 0.95*0.000184222, 0.95*6.89994e-05, 0.95*2.24488e-05, 0.95*6.41171e-06, 0.95*2.01349e-06, 0.95*5.26691e-07};
const double THXSEC13TeV_GMStau_High [] = { 1.05*14.2831, 1.05*3.90949, 1.05*1.20093, 1.05*0.283988, 1.05*0.0882234, 1.05*0.0249859, 1.05*0.00907959, 1.05*0.00390754, 1.05*0.0018932, 1.05*0.000990061, 1.05*0.000545049, 1.05*0.000409813, 1.05*0.000184222, 1.05*6.89994e-05, 1.05*2.24488e-05, 1.05*6.41171e-06, 1.05*2.01349e-06, 1.05*5.26691e-07};


//LO PairProduced Stau@7TeV
const double THXSEC7TeV_PPStau_Mass [] = {100,126,156,200,247,308, 370, 432, 494, 557};
const double THXSEC7TeV_PPStau_Cen  [] = {0.038200,0.0161,0.007040,0.002470,0.001010,0.000353, 0.00014, 6.01e-05, 2.75e-05, 1.31E-05};
const double THXSEC7TeV_PPStau_Low  [] = {0.037076,0.0155927,0.0067891,0.00237277,0.00096927,0.000335308, 0.000132238,  5.62644e-05, 2.56415e-05, 1.31E-05*(2.10179e-05/2.24e-05)};
const double THXSEC7TeV_PPStau_High [] = {0.0391443,0.016527,0.00723151,0.00253477,0.00103844,0.000363699, 0.00014685, 6.31201e-05, 2.90217e-05, 1.31E-05*(2.32576e-05/2.24e-05)};

const double THXSEC8TeV_PPStau_Mass [] = {100, 126, 156, 200, 247, 308, 370, 432, 494, 557};
const double THXSEC8TeV_PPStau_Cen  [] = {0.047, 0.0202, 0.00902, 0.00327, 0.00137, 0.000503, 0.000208, 9.36e-05, 4.47e-05, 2.24e-05};
const double THXSEC8TeV_PPStau_Low  [] = {0.0456616, 0.0196246, 0.00873146, 0.00315596, 0.00130957, 0.000479108, 0.000197003, 8.85171e-05, 4.21191e-05, 2.10179e-05};
const double THXSEC8TeV_PPStau_High [] = {0.0481431, 0.0206482, 0.00924512, 0.00335015, 0.00140238, 0.000516666, 0.000214196, 9.66572e-05, 4.63125e-05, 2.32576e-05};

//FIXME error must be computed properly, currently it's just Cen+-5%
const double THXSEC13TeV_PPStau_Mass [] = {100, 126, 156, 200, 247, 308, 370, 432, 494, 557, 619, 651, 745, 871, 1029, 1218, 1409, 1599};
const double THXSEC13TeV_PPStau_Cen  [] = {0.0944, 0.0426, 0.0202, 0.00797, 0.00365, 0.0015, 0.0007, 0.000354, 0.000191, 0.000108, 6.36e-05, 4.93e-05, 2.4e-05, 9.87e-06, 3.54e-06, 1.14e-06, 3.91e-07, 1.41e-07};
const double THXSEC13TeV_PPStau_Low  [] = {0.95*0.0944, 0.95*0.0426, 0.95*0.0202, 0.95*0.00797, 0.95*0.00365, 0.95*0.0015, 0.95*0.0007, 0.95*0.000354, 0.95*0.000191, 0.95*0.000108, 0.95*6.36e-05, 0.95*4.93e-05, 0.95*2.4e-05, 0.95*9.87e-06, 0.95*3.54e-06, 0.95*1.14e-06, 0.95*3.91e-07, 0.95*1.41e-07};
const double THXSEC13TeV_PPStau_High [] = {1.05*0.0944, 1.05*0.0426, 1.05*0.0202, 1.05*0.00797, 1.05*0.00365, 1.05*0.0015, 1.05*0.0007, 1.05*0.000354, 1.05*0.000191, 1.05*0.000108, 1.05*6.36e-05, 1.05*4.93e-05, 1.05*2.4e-05, 1.05*9.87e-06, 1.05*3.54e-06, 1.05*1.14e-06, 1.05*3.91e-07, 1.05*1.41e-07};

#endif


