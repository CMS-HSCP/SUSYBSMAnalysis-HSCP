// Original Author:  Loic Quertenmont

#include "Analysis_Global.h"
#include "Analysis_PlotFunction.h"
#include "TVector3.h"
#include <string>

////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// general purpose code 
std::string ReplacePartOfString (std::string s, std::string a, std::string b){
   size_t pos=0;
   string toReturn=s;
   while(1){
      pos = toReturn.find(a, pos);
      if (pos==std::string::npos) break;
      toReturn.replace(pos++,a.size(),b);
   }
   return toReturn;
}

std::string cleanSampleName (std::string s){
   string toReturn = ReplacePartOfString (s, "_13TeV16G", "");
   toReturn = ReplacePartOfString (toReturn, "_13TeV16", "");
   toReturn = ReplacePartOfString (toReturn, "_13TeV", "");
   return toReturn;
}

// return the TypeMode from a string inputPattern
int TypeFromPattern(const std::string& InputPattern){
   if(InputPattern.find("Type0",0)<std::string::npos){       return 0;
   }else if(InputPattern.find("Type1",0)<std::string::npos){ return 1;
   }else if(InputPattern.find("Type2",0)<std::string::npos){ return 2;
   }else if(InputPattern.find("Type3",0)<std::string::npos){ return 3;
   }else if(InputPattern.find("Type4",0)<std::string::npos){ return 4;
   }else if(InputPattern.find("Type5",0)<std::string::npos){ return 5;
   }else{                                                    return 6;
   }
}

// define the legend corresponding to a Type
std::string LegendFromType(const std::string& InputPattern){
   switch(TypeFromPattern(InputPattern)){
      case 0:  return std::string("Tracker - Only"); break;
      case 1:  return std::string("Tracker + Muon"); break;
      case 2:  return std::string("Tracker + TOF" ); break;
      case 3:  return std::string("Muon - Only"); break;
      case 4:  return std::string("|Q|>1e"); break;
      case 5:  return std::string("|Q|<1e"); break;
      default : std::string("unknown");
   }
   return std::string("unknown");
}

// compute deltaR between two point (eta,phi) (eta,phi)
double deltaR(double eta1, double phi1, double eta2, double phi2) {
   double deta = eta1 - eta2;
   double dphi = phi1 - phi2;
   while (dphi >   M_PI) dphi -= 2*M_PI;
   while (dphi <= -M_PI) dphi += 2*M_PI;
   return sqrt(deta*deta + dphi*dphi);
}

// function to go form a TH3 plot with cut index on the x axis to a  TH2
TH2D* GetCutIndexSliceFromTH3(TH3D* tmp, unsigned int CutIndex, string Name="zy"){
   tmp->GetXaxis()->SetRange(CutIndex+1,CutIndex+1);
   return (TH2D*)tmp->Project3D(Name.c_str());
}

// function to go form a TH2 plot with cut index on the x axis to a  TH1
TH1D* GetCutIndexSliceFromTH2(TH2D* tmp, unsigned int CutIndex, string Name="_py"){
   return tmp->ProjectionY(Name.c_str(),CutIndex+1,CutIndex+1);
}


////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// Genertic code for beta/mass reconstruction starting from p&dEdx or p&TOF

// compute mass out of a beta and momentum value
double GetMassFromBeta(double P, double beta){
   double gamma = 1/sqrt(1-beta*beta);
   return P/(beta*gamma);
} 

// compute mass out of a momentum and tof value
double GetTOFMass(double P, double TOF){
   return GetMassFromBeta(P, 1/TOF);
}

double GetMassErr (double P, double PErr, double dEdx, double dEdxErr, double M){
   if (M < 0) return -1;
   double KErr     = 0.2;
   double CErr     = 0.4;
   double cErr     = 0.01;
   double Criteria = dEdx - dEdxC_Data;
   double Fac1     = P*P/(2*M*dEdxK_Data);
   double Fac2     = pow(2*M*M*dEdxK_Data/(P*P), 2);
   double MassErr  = Fac1*sqrt(Fac2*pow(PErr/P, 2) + Criteria*Criteria*pow(KErr/dEdxK_Data,2) + dEdxErr*dEdxErr + dEdxC_Data*dEdxC_Data);

   if (std::isnan(MassErr) || std::isinf(MassErr)) MassErr = -1;

   return MassErr/M;
}

// estimate beta from a dEdx value, if dEdx is below the allowed threshold returns a negative beta value
double GetIBeta(double I, bool MC){
   double& K = dEdxK_Data;
   double& C = dEdxC_Data;
   if(MC){ K = dEdxK_MC;
           C = dEdxC_MC;  }

   double a = K / (I-C);
   double b2 = a / (a+1);
   if(b2<0)return -1*sqrt(b2);
   return sqrt(b2);
}

// compute mass out of a momentum and dEdx value
double GetMass(double P, double I, bool MC){
   double& K = dEdxK_Data;
   double& C = dEdxC_Data;
   if(MC){ K = dEdxK_MC;
           C = dEdxC_MC;  }

   if(I-C<0)return -1;
   return sqrt((I-C)/K)*P;
}

// pz compute Ick out of dEdx value
double GetIck(double I, bool MC){
   double& K = dEdxK_Data;
   double& C = dEdxC_Data;
   if(MC){ K = dEdxK_MC;
           C = dEdxC_MC;  }

   return (I-C)/K;
}


// return a TF1 corresponding to a mass line in the momentum vs dEdx 2D plane
TF1* GetMassLine(double M, bool MC){
   double& K = dEdxK_Data;
   double& C = dEdxC_Data;
   if(MC){ K = dEdxK_MC;
           C = dEdxC_MC;  }

   double BetaMax = 0.9;
   double PMax = sqrt((BetaMax*BetaMax*M*M)/(1-BetaMax*BetaMax));

   double BetaMin = 0.2;
   double PMin = sqrt((BetaMin*BetaMin*M*M)/(1-BetaMin*BetaMin));

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


TF1* GetMassLineQ(double M, double Charge=1, bool MC=false)
{
   double K;   double C;
   if(MC){
      K = dEdxK_MC;
      C = dEdxC_MC;
   }else{ 
      K = dEdxK_Data;
      C = dEdxC_Data;
   }

   double BetaMax = 0.999;
   double PMax = sqrt((BetaMax*BetaMax*M*M)/(1-BetaMax*BetaMax));

   double BetaMin = 0.01;
   double PMin = sqrt((BetaMin*BetaMin*M*M)/(1-BetaMin*BetaMin));

   TF1* MassLine = new TF1("MassLine","[3] * ([2] + ([0]*[0]*[1])/(x*x*[3]))", PMin, PMax);
   MassLine->SetParName  (0,"M");
   MassLine->SetParName  (1,"K");
   MassLine->SetParName  (2,"C");
   MassLine->SetParName  (3,"z2");
   MassLine->SetParameter(0, M);
   MassLine->SetParameter(1, K);
   MassLine->SetParameter(2, C);
   MassLine->SetParameter(3, Charge*Charge);
   MassLine->SetLineWidth(2);
   return MassLine;
}







////////////////////////////////////////////////////////////////////////////////////////////////////////// 
// Functions below were used for the 2009 and 2010 papers, but are probably not used anymore 

// return the selection efficiency for a given disribution (TH1) and for a given cut... Nothing but the integral above the cut divided by the number of events
double Efficiency(TH1* Histo, double CutX){
   double Entries  = Histo->Integral(0,Histo->GetNbinsX()+1);
   double Integral = Histo->Integral(Histo->GetXaxis()->FindBin(CutX),Histo->GetNbinsX()+1);
   return Integral/Entries;
}

// same as before but also compute the error in the efficiency
double EfficiencyAndError(TH1* Histo, double CutX, double& error){
   double Entries  = Histo->Integral(0,Histo->GetNbinsX()+1);
   double Integral = 0;
          error    = 0;
   for(Int_t binx = Histo->GetXaxis()->FindBin(CutX); binx<= Histo->GetNbinsX()+1; ++binx){
      Integral += Histo->GetBinContent(binx);
      error    += Histo->GetBinError(binx)*Histo->GetBinError(binx);
   }
   error = sqrt(error);
   error /= Entries;
   return Integral/Entries;
}

// return the selection efficiency for a given disribution (TH2) and for a given recangular signal region ... Nothing but the integral above the cut divided by the number of events
double Efficiency(TH2* Histo, double CutX, double CutY){
   double Entries  = Histo->Integral(0,Histo->GetNbinsX()+1, 0,Histo->GetNbinsY()+1);
   double Integral = Histo->Integral(Histo->GetXaxis()->FindBin(CutX),Histo->GetNbinsX()+1, Histo->GetYaxis()->FindBin(CutY),Histo->GetNbinsY()+1);
   return Integral/Entries;
}

// return the number of entry in an histogram (and it's error) in a window defined by two cuts
double GetEventInRange(double min, double max, TH1D* hist, double& error){
  int binMin = hist->GetXaxis()->FindBin(min);
  int binMax = hist->GetXaxis()->FindBin(max);
  error = 0; for(int i=binMin;i<binMax;i++){ error += pow(hist->GetBinError(i),2); }  error = sqrt(error);
  return hist->Integral(binMin,binMax);
}

// not used anymore, was computing a scale factor between datadriven prediction and observation using the M<75GeV events
void GetPredictionRescale(std::string InputPattern, double& Rescale, double& RMS, bool ForceRecompute=false)
{
   size_t CutIndex = InputPattern.find("/Type");
   InputPattern    = InputPattern.substr(0,CutIndex+7);
   std::string Input    = InputPattern + "PredictionRescale.txt";


   FILE* pFile = fopen(Input.c_str(),"r");
   if(pFile && !ForceRecompute){
      float tmp1, tmp2;
      fscanf(pFile,"Rescale=%f RMS=%f\n",&tmp1,&tmp2);
      Rescale = tmp1;
      RMS = tmp2;
      fclose(pFile);
   }else{
      Rescale = 0;
      RMS     = 0;
      int    NPoints = 0;

      std::vector<double> DValue;
      std::vector<double> PValue;
  
      for(float WP_Pt=0;WP_Pt>=-5;WP_Pt-=0.5f){
      for(float WP_I =0;WP_I >=-5;WP_I -=0.5f){
         char Buffer[2048];
         sprintf(Buffer,"%sWPPt%+03i/WPI%+03i/DumpHistos.root",InputPattern.c_str(),(int)(10*WP_Pt),(int)(10*WP_I));
         TFile* InputFile = new TFile(Buffer); 
         if(!InputFile || InputFile->IsZombie() || !InputFile->IsOpen() || InputFile->TestBit(TFile::kRecovered) )continue;

         double d=0, p=0;//, m=0;
         double error =0;
         TH1D* Hd = (TH1D*)GetObjectFromPath(InputFile, "Mass_Data");if(Hd){d=GetEventInRange(0,75,Hd,error);delete Hd;}
         TH1D* Hp = (TH1D*)GetObjectFromPath(InputFile, "Mass_Pred");if(Hp){p=GetEventInRange(0,75,Hp,error);delete Hp;}
//       TH1D* Hm = (TH1D*)GetObjectFromPath(InputFile, "Mass_MCTr");if(Hm){m=GetEventInRange(0,75,Hm);delete Hm;}

//       if(!(d!=d) && p>0 && d>10 && (WP_Pt+WP_I)<=-3){
//         if(!(d!=d) && p>0 && d>20 && (WP_Pt+WP_I)<=-3){
         if(!(d!=d) && p>0 && d>20 && (WP_Pt+WP_I)<=-2){
//         if(!(d!=d) && p>0 && d>500 && (WP_Pt+WP_I)<=-2){
            DValue.push_back(d);
            PValue.push_back(p);
            printf("%6.2f %6.2f (eff=%6.2E) --> %f  (d=%6.2E)\n",WP_Pt,WP_I, pow(10,WP_Pt+WP_I),d/p, d);
            Rescale += (d/p);
            NPoints++;
         }
         InputFile->Close();
      }}
      printf("----------------------------\n");
      Rescale /= NPoints;

      for(unsigned int i=0;i<DValue.size();i++){
          RMS += pow( (DValue[i]/(PValue[i]*Rescale)) - 1.0 ,2);
      }
      RMS /= NPoints;
      RMS = sqrt(RMS);

      pFile = fopen(Input.c_str(),"w");
      if(!pFile)return;
      fprintf(pFile,"Rescale=%6.2f RMS=%6.2f\n",Rescale,RMS);
      fclose(pFile);
   }
   printf("Mean Rescale = %f   RMS = %f\n",Rescale, RMS);
}

// find the intersection between two graphs... very useful to know what is the excluded mass range from an observed xsection limit and theoretical xsection
double FindIntersectionBetweenTwoGraphs(TGraph* obs, TGraph* th, double Min, double Max, double Step, double ThUncertainty=0, bool debug=false){

   double Intersection = -1;
   double ThShift = 1.0-ThUncertainty;
   double PreviousX = Min;
   double PreviousV = obs->Eval(PreviousX, 0, "") - (ThShift * th->Eval(PreviousX, 0, "")) ;
   if(PreviousV>0){if(debug){printf("FindIntersectionBetweenTwoGraphs returns -1 because observed xsection is above th xsection for the first mass already : %f vs %f\n", obs->Eval(PreviousX, 0, ""), th->Eval(PreviousX, 0, ""));}return -2;}
   for(double x=Min+=Step;x<Max;x+=Step){                 
      double V = obs->Eval(x, 0, "") - (ThShift * th->Eval(x, 0, "") );
      if(debug){
         printf("%7.2f --> Obs=%6.2E  Th=%6.2E",x,obs->Eval(x, 0, ""),ThShift * th->Eval(x, 0, ""));
         if(V>=0)printf("   X\n");
         else printf("\n");
      }
      if(V<0){
         PreviousX = x;
         PreviousV = V;
      }else{
         Intersection = PreviousX;
      }
   }
   if(Intersection!=-1)return Intersection;
   return PreviousV>0 ? Intersection : -1*Max;
}


// find the range of excluded masses. Necessary when low mass not excluded but higher masses are
void FindRangeBetweenTwoGraphs(TGraph* obs, TGraph* th, double Min, double Max, double Step, double ThUncertainty, double& minExclMass, double& maxExclMass){

   double ThShift = 1.0-ThUncertainty;
   double PreviousX = Min;
   double PreviousV = obs->Eval(PreviousX, 0, "") - (ThShift * th->Eval(PreviousX, 0, "")) ;
   bool LowMassNeeded=true;
   if(PreviousV<0) LowMassNeeded=false;
   for(double x=Min+=Step;x<Max;x+=Step){                 
      double V = obs->Eval(x, 0, "") - (ThShift * th->Eval(x, 0, "") );
      if(LowMassNeeded && V<0) {
	minExclMass=x;
	LowMassNeeded=false;
      }
      if(V<0){
         PreviousX = x;
         PreviousV = V;
      }else{
         maxExclMass = PreviousX;
      }
   }
   return;
}


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Genertic code related to samples processing in FWLITE --> functions below will be loaded only if FWLITE compiler variable is defined

#ifdef FWLITE
bool isGoodGenHSCP(const reco::GenParticle& gen, bool onlyCharged=false);   
double DistToHSCP        (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest);
int    HowManyChargedHSCP(const std::vector<reco::GenParticle>& genColl);
double FastestHSCP       (const fwlite::ChainEvent& ev);
void   GetGenHSCPBeta    (const std::vector<reco::GenParticle>& genColl, double& beta1, double& beta2, bool onlyCharged=true);
void   GetGenHSCPDecayLength (const std::vector<reco::GenParticle>& genColl, double& length1, double& length2, bool onlyCharged);    

// compute the distance between a "reconstructed" HSCP candidate and the closest generated HSCP
double DistToHSCP (const susybsm::HSCParticle& hscp, const std::vector<reco::GenParticle>& genColl, int& IndexOfClosest){
   reco::TrackRef   track;
   if(TypeMode!=3) track = hscp.trackRef();
   else {
     reco::MuonRef muon = hscp.muonRef();
     if(muon.isNull()) return false;
     track = muon->standAloneMuon();
   }
   if(track.isNull())return false;

   double RMin = 9999; IndexOfClosest=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(genColl[g].pt()<5)continue;
      if(genColl[g].status()!=1)continue;
      int AbsPdg=abs(genColl[g].pdgId());
      if(AbsPdg<1000000 && AbsPdg!=17)continue;    
      double dR = deltaR(track->eta(), track->phi(), genColl[g].eta(), genColl[g].phi());
      if(dR<RMin){RMin=dR;IndexOfClosest=g;}
   }
   return RMin;
}

//selection of gen HSCP
bool isGoodGenHSCP(const reco::GenParticle& gen, bool onlyCharged){
   if(gen.status()!=1)return false;
   if(gen.pt()<5)return false;
   int AbsPdg=abs(gen.pdgId());
   if(AbsPdg<1000000 && AbsPdg!=17)return false;
   if(onlyCharged && (AbsPdg==1000993 || AbsPdg==1009313 || AbsPdg==1009113 || AbsPdg==1009223 || AbsPdg==1009333 || AbsPdg==1092114 || AbsPdg==1093214 || AbsPdg==1093324))return false; //Skip neutral gluino RHadrons
   if(onlyCharged && (AbsPdg==1000622 || AbsPdg==1000642 || AbsPdg==1006113 || AbsPdg==1006311 || AbsPdg==1006313 || AbsPdg==1006333))return false;  //skip neutral stop RHadrons
   if(onlyCharged && AbsPdg==1000022)return false; //skip neutralino
   return true;
}

// count the number of charged generated HSCP in the event --> this is needed to reweights the events for different gluino ball fraction starting from f=10% samples
int HowManyChargedHSCP (const std::vector<reco::GenParticle>& genColl){
   int toReturn = 0;
   for(unsigned int g=0;g<genColl.size();g++){
      if(isGoodGenHSCP(genColl[g], true))toReturn++;
   }
   return toReturn;
}


// returns the generated beta of the two firsts HSCP in the events
void  GetGenHSCPBeta (const std::vector<reco::GenParticle>& genColl, double& beta1, double& beta2, bool onlyCharged){
   beta1=-1; beta2=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(!isGoodGenHSCP(genColl[g], onlyCharged))continue;
      if(beta1<0){beta1=genColl[g].p()/genColl[g].energy();}else if(beta2<0){beta2=genColl[g].p()/genColl[g].energy();return;}
   }
}


// returns the generated decay length of the two firsts HSCP in the events
 void  GetGenHSCPDecayLength (const std::vector<reco::GenParticle>& genColl, double& length1, double& length2, bool onlyCharged){
   length1=-1; length2=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(!isGoodGenHSCP(genColl[g], onlyCharged))continue;

      //code below copied from disapearing track source code.  FIXME we should make sure that the logic is ok.
      //https://raw.githubusercontent.com/OSU-CMS/DisappTrks/master/SignalMC/plugins/DecayAnalyzer.cc  
      const reco::Candidate* mother   = &genColl[g];
      const reco::Candidate* daughter = &genColl[g];

      // Descend the decay chain until no daughters have the same PDG ID as mcParticle.
      while (true) {
        bool foundDauSamePdgId = false;
        for (unsigned int i=0; i<daughter->numberOfDaughters(); i++) {
          if (daughter->daughter(i)->pdgId() == genColl[g].pdgId()) {
            foundDauSamePdgId = true;
            mother = daughter;
            daughter = daughter->daughter(i);
            break;
          }
        }
        if (!foundDauSamePdgId) break;
      }

      // Now daughter has no daughters with the same PDG ID as mcParticle.
      // Next choose the daughter with the outermost production vertex, in case there are multiple vertices
      // (e.g., an electron delta ray can produce a vertex before the decay vertex)
       double radiusLastVtx = -99;
       int idxDauLastVtx = -99;
       for(unsigned int i=0; i<daughter->numberOfDaughters(); i++) {
         double radius = daughter->daughter(i)->vertex().R();
         if (radius > radiusLastVtx) {
           radiusLastVtx = radius;
           idxDauLastVtx = i;
         }
       }
       if(idxDauLastVtx<0)continue;
 
       mother = daughter;
       daughter = daughter->daughter(idxDauLastVtx);

       TVector3 source (genColl[g].vx (), genColl[g].vy (), genColl[g].vz ());
       TVector3 decay  (daughter->vx (), daughter->vy (), daughter->vz ()); 
       double ctau = (decay - source).Mag () / (genColl[g].p4 ().Beta () * genColl[g].p4 ().Gamma ());

      if(length1<0){length1=ctau;}else if(length2<0){length2=ctau;return;}
   }
   if(length1<0){length1=9999;}
   if(length2<0){length2=9999;}
}

 

double FastestHSCP(const fwlite::ChainEvent& ev){
   //get the collection of generated Particles
   fwlite::Handle< std::vector<reco::GenParticle> > genCollHandle;
   genCollHandle.getByLabel(ev, "genParticlesSkimmed");
   if(!genCollHandle.isValid()){
      genCollHandle.getByLabel(ev, "genParticles");
      if(!genCollHandle.isValid()){printf("GenParticle Collection NotFound (from FastestHSCP)\n");return -1;}
   }

   std::vector<reco::GenParticle> genColl = *genCollHandle;

   double MaxBeta=-1;
   for(unsigned int g=0;g<genColl.size();g++){
      if(!isGoodGenHSCP(genColl[g]))continue;
      double beta=genColl[g].p()/genColl[g].energy();
      if(MaxBeta<beta)MaxBeta=beta;
   }
   return MaxBeta;
}

#include "FWCore/Utilities/interface/RegexMatch.h"
bool passTriggerPatterns(edm::TriggerResultsByName& tr, std::string pattern){
  if(edm::is_glob(pattern)){
     std::vector< std::vector<std::string>::const_iterator > matches = edm::regexMatch(tr.triggerNames(), pattern);
     for(size_t t=0;t<matches.size();t++){
        if(tr.accept( matches[t]->c_str() ) )return true;
     }
  }else{
     if(tr.accept( pattern.c_str() ) ) return true;
  }
  return false;
}

double deltaROpositeTrack(const susybsm::HSCParticleCollection& hscpColl, const susybsm::HSCParticle& hscp){
   reco::TrackRef track1=hscp.trackRef();

   double maxDr=-0.1;
   for(unsigned int c=0;c<hscpColl.size();c++){
      reco::TrackRef track2;
      if(!hscpColl[c].trackRef().isNull()){
         track2=hscpColl[c].trackRef();
      }else if(!hscpColl[c].muonRef().isNull() && hscpColl[c].muonRef()->combinedQuality().updatedSta){
         track2= hscpColl[c].muonRef()->standAloneMuon();
      }else{
         continue;
      }

      if(fabs(track1->pt()-track2->pt())<1 && deltaR(track1->eta(), track1->phi(), track2->eta(), track2->phi())<0.1)continue; //Skip same tracks
//      double dR = deltaR(-1*track1->eta(), M_PI+track1->phi(), track2->eta(), track2->phi());
      TVector3 v1 = TVector3(track1->momentum().x(), track1->momentum().y(), track1->momentum().z());
      TVector3 v2 = TVector3(track2->momentum().x(), track2->momentum().y(), track2->momentum().z());
      double dR = v1.Angle(v2);
      if(dR>maxDr)maxDr=dR;
   }
   return maxDr;
}

#endif


//////////////////////////////////////////////////////////////////////////////////////////////////////////
// Handfull class to check for duplicated events

class DuplicatesClass{
   private :
      typedef std::map<std::pair<unsigned int, unsigned int>, int > RunEventHashMap;
      RunEventHashMap map;
   public :
        DuplicatesClass(){}
        ~DuplicatesClass(){}
        void Clear(){map.clear();}
        bool isDuplicate(unsigned int Run, unsigned int Event){
	   RunEventHashMap::iterator it = map.find(std::make_pair(Run,Event));
           if(it==map.end()){
   	      map[std::make_pair(Run,Event)] = 1;
              return false;
           }else{
              map[std::make_pair(Run,Event)]++;
           }
           return true;
        }

        void printDuplicate(){
           printf("Duplicate event summary:\n##########################################");
           for(RunEventHashMap::iterator it = map.begin(); it != map.end(); it++){
              if(it->second>1)printf("Run %6i Event %10i is duplicated (%i times)\n",it->first.first, it->first.second, it->second);
           }          
           printf("##########################################");
        }
};
 

#ifdef FWLITE

struct HitDeDx{double dedx; double dx; bool isSat; bool passClusterCleaning; bool isInside; unsigned char subDet;};
typedef std::vector<HitDeDx> HitDeDxCollection;



class dedxHIPEmulator{
   private:
      TH1D* ratePdfPixel;
      TH1D* ratePdfStrip;

      double eventRatePixel;
      double eventRateStrip;

      bool is2016G, is2016;
   public:
     void setEventRate(double ratePixel=-1, double rateStrip=-1){
//           std::cout <<__LINE__<< " setEventRate  " << is2016 <<" "<< is2016G << std::endl;

        if(ratePixel<0){
           eventRatePixel = ratePdfPixel->GetRandom();
           eventRatePixel -= 3.2;//2.4; //subtract rate already present in the MC
           eventRatePixel *= 100; //for random generator usage
//	   std::cout <<__LINE__<< " pixelrate  " << eventRatePixel << std::endl;

        }else{eventRatePixel = ratePixel;}
        if(rateStrip<0){
           eventRateStrip = ratePdfStrip->GetRandom();
           eventRateStrip -= 1.1;//0.8; //subtract rate already present in the MC
           eventRateStrip *= 100; //for random generator usage
//           std::cout <<__LINE__<< " striprate  " << eventRateStrip << std::endl;

        }else{eventRateStrip = rateStrip;}  

        eventRatePixel = std::max(eventRatePixel, 0.0);
        eventRateStrip = std::max(eventRateStrip, 0.0);
//           std::cout <<__LINE__<< " MAX pixelrate  " << eventRatePixel << std::endl;
//           std::cout <<__LINE__<< " MAX striprate  " << eventRateStrip << std::endl;

     }

     void setPeriodHIPRate(bool is2016_, const char* ratePdfPixelName=NULL, const char* ratePdfStripName=NULL){
      //the meanig of "is2016_" changed (20190911): now "true" means that it is 2016_PostG, "false" : 2016_PreG
      // the procedure can not be used for 2015
        if (is2016G == is2016_ && ratePdfPixel && ratePdfStrip) return;
        else {
           is2016G = is2016_;
	   char tmp1 [1024];
	   char tmp2 [1024];
	   if (ratePdfPixel && !ratePdfPixelName) sprintf (tmp1, "%s_%s", ratePdfPixel->GetName(), is2016G?"2016postG":"2016pre");
	   else if (ratePdfPixelName)             sprintf (tmp1, "%s_%s", ratePdfPixelName, is2016G?"2016postG":"2016pre");
	   else                                   sprintf (tmp1, "ratePdfPixel_%s", is2016G?"2016postG":"2016pre");

	   if (ratePdfStrip && !ratePdfStripName) sprintf (tmp2, "%s_%s", ratePdfStrip->GetName(), is2016G?"2016postG":"2016pre");
	   else if (ratePdfStripName)             sprintf (tmp2, "%s_%s", ratePdfStripName, is2016G?"2016postG":"2016pre");
	   else                                   sprintf (tmp2, "ratePdfStrip_%s",is2016G?"2016postG":"2016pre");

	   std::cout << "In Histo1 name is " << tmp1 << std::endl;
	   std::cout << "In Histo2 name is " << tmp2 << std::endl;
           if (ratePdfPixel) {ratePdfPixel->SetBins(26,0,13); ratePdfPixel->SetNameTitle(tmp1, tmp1);} else ratePdfPixel = new TH1D(tmp1, tmp1,26,0,13);
           if (ratePdfStrip) {ratePdfStrip->SetBins(26,0,13); ratePdfStrip->SetNameTitle(tmp2, tmp2);} else ratePdfStrip = new TH1D(tmp2, tmp2,26,0,13);
           if (!is2016G){  // PreG

                    ratePdfPixel->SetBinContent(3,2435235);
                    ratePdfPixel->SetBinContent(5,80981);
                    ratePdfPixel->SetBinContent(6,4474);
                    ratePdfPixel->SetBinContent(7,881461);
                    ratePdfPixel->SetBinContent(8,1566102);
                    ratePdfPixel->SetBinContent(9,1398653);
                    ratePdfPixel->SetBinContent(10,1148772);
                    ratePdfPixel->SetBinContent(11,1410443);
                    ratePdfPixel->SetBinContent(12,1533880);
                    ratePdfPixel->SetBinContent(13,1983972);
                    ratePdfPixel->SetBinContent(14,1166108);
                    ratePdfPixel->SetBinContent(15,1356109);
                    ratePdfPixel->SetBinContent(16,176708);
                    ratePdfPixel->SetBinContent(17,784558);
                    ratePdfPixel->SetBinContent(18,267585);
                    ratePdfPixel->SetBinContent(19,327523);
                    ratePdfPixel->SetBinContent(20,378951);
                    ratePdfPixel->SetBinContent(21,123427);
                    ratePdfPixel->Scale(1.0/ratePdfPixel->Integral());

                    ratePdfStrip->SetBinContent(2,2435235);
                    ratePdfStrip->SetBinContent(3,667);
                    ratePdfStrip->SetBinContent(4,126408);
                    ratePdfStrip->SetBinContent(5,184591);
                    ratePdfStrip->SetBinContent(6,33254);
                    ratePdfStrip->SetBinContent(7,34475);
                    ratePdfStrip->SetBinContent(8,266368);
                    ratePdfStrip->SetBinContent(9,162917);
                    ratePdfStrip->SetBinContent(10,536416);
                    ratePdfStrip->SetBinContent(11,302761);
                    ratePdfStrip->SetBinContent(12,876949);
                    ratePdfStrip->SetBinContent(13,712690);
                    ratePdfStrip->SetBinContent(14,645763);
                    ratePdfStrip->SetBinContent(15,1083180);
                    ratePdfStrip->SetBinContent(16,711389);
                    ratePdfStrip->SetBinContent(17,1347154);
                    ratePdfStrip->SetBinContent(18,1996104);
                    ratePdfStrip->SetBinContent(19,1320364);
                    ratePdfStrip->SetBinContent(20,688952);
                    ratePdfStrip->SetBinContent(21,471939);
                    ratePdfStrip->SetBinContent(22,1180686);
                    ratePdfStrip->SetBinContent(23,676072);
                    ratePdfStrip->SetBinContent(24,651857);
                    ratePdfStrip->SetBinContent(25,118120);
                    ratePdfStrip->SetBinContent(26,210558);
                    ratePdfStrip->SetBinContent(27,500146);
                    ratePdfStrip->Scale(1.0/ratePdfStrip->Integral());



                 } else {            //PostG


                    ratePdfPixel->SetBinContent(3,2435235);
                    ratePdfPixel->SetBinContent(5,12465);
                    ratePdfPixel->SetBinContent(7,21122);
                    ratePdfPixel->SetBinContent(8,67380);
                    ratePdfPixel->SetBinContent(9,335441);
                    ratePdfPixel->SetBinContent(10,603035);
                    ratePdfPixel->SetBinContent(11,1420760);
                    ratePdfPixel->SetBinContent(12,1609384);
                    ratePdfPixel->SetBinContent(13,2606849);
                    ratePdfPixel->SetBinContent(14,1642034);
                    ratePdfPixel->SetBinContent(15,2021694);
                    ratePdfPixel->SetBinContent(16,1972932);
                    ratePdfPixel->SetBinContent(17,795034);
                    ratePdfPixel->SetBinContent(18,973985);
                    ratePdfPixel->SetBinContent(19,597100);
                    ratePdfPixel->SetBinContent(20,76477);
                    ratePdfPixel->SetBinContent(21,8446);
                    ratePdfPixel->Scale(1.0/ratePdfPixel->Integral());

                    ratePdfStrip->SetBinContent(2,3195357);
                    ratePdfStrip->SetBinContent(3,1676482);
                    ratePdfStrip->SetBinContent(4,8859200);
                    ratePdfStrip->SetBinContent(5,2549081);
                    ratePdfStrip->SetBinContent(15,80425);
                    ratePdfStrip->SetBinContent(16,53688);
                    ratePdfStrip->SetBinContent(19,77926);
                    ratePdfStrip->SetBinContent(21,223841);
                    ratePdfStrip->SetBinContent(22,276528);
                    ratePdfStrip->SetBinContent(24,108671);
                    ratePdfStrip->SetBinContent(27,196348);
                    ratePdfStrip->Scale(1.0/ratePdfStrip->Integral());



           }

        }
     }

     dedxHIPEmulator(bool is2016_=false, const char* ratePdfPixelName = NULL, const char* ratePdfStripName = NULL){
         setPeriodHIPRate(is2016_, ratePdfPixelName, ratePdfStripName);
     }
     ~dedxHIPEmulator(){}



     double getEventRatePixel() { return eventRatePixel; }
     double getEventRateStrip() { return eventRateStrip; }

     double fakeHIP(unsigned int subDet, double dedx){
 
        if(subDet< 3 && rand()%10000<eventRatePixel)dedx = ( 0.6 + ((rand()%15000)/10000.0) );
        if(subDet>=3 && rand()%10000<eventRateStrip)dedx = ( 0.3 + ((rand()%20000)/10000.0) );

        return dedx;
     }

     void fakeHIP(HitDeDxCollection& hitDeDx){
//           std::cout <<__LINE__<< " fakeHIP dla kolekcji hitow wolane  "  << std::endl;

        for(unsigned int h=0;h<hitDeDx.size();h++){
           hitDeDx[h].dedx = fakeHIP(hitDeDx[h].subDet, hitDeDx[h].dedx);
       }
    }

  
};






TH3F* loadDeDxTemplate(string path, bool splitByModuleType=false);
reco::DeDxData computedEdx(const DeDxHitInfo* dedxHits, double* scaleFactors, TH3* templateHisto=NULL, bool usePixel=false, bool useClusterCleaning=true, bool reverseProb=false, bool useTruncated=false, std::unordered_map<unsigned int,double>* TrackerGains=NULL, bool useStrip=true, bool mustBeInside=false, size_t MaxStripNOM=999, bool correctFEDSat=false, int crossTalkInvAlgo=0, double dropLowerDeDxValue=0.0, dedxHIPEmulator* hipEmulator=NULL, double* dEdxErr = NULL, unsigned int pdgId=0);
HitDeDxCollection getHitDeDx(const DeDxHitInfo* dedxHits, double* scaleFactors, std::unordered_map<unsigned int,double>* TrackerGains=NULL, bool correctFEDSat=false, int crossTalkInvAlgo=0);

bool clusterCleaning(const SiStripCluster*   cluster,  int crosstalkInv=0, uint8_t* exitCode=NULL);
void printStripCluster(FILE* pFile, const SiStripCluster*   cluster, const DetId& DetId, bool crossTalkInvAlgo);
void printClusterCleaningMessage (uint8_t exitCode);
std::vector<int> convert(const vector<unsigned char>& input);
std::vector<int> CrossTalkInv(const std::vector<int>&  Q, const float x1=0.10, const float x2=0.04, bool way=true,float threshold=20,float thresholdSat=25);


class dedxGainCorrector{
   private:
      std::map<unsigned int, std::unordered_map<unsigned int, double> > TrackerGainsPerRuns;

   public:
      std::unordered_map<unsigned int, double>* TrackerGains; 
      dedxGainCorrector(){TrackerGains=NULL;}
      ~dedxGainCorrector(){}

      void setRun(unsigned int currentRun){
         if(TrackerGainsPerRuns.size()<=0){TrackerGains=NULL; return;}
         std::map<unsigned int, std::unordered_map<unsigned int, double> >::iterator it, itPrev=TrackerGainsPerRuns.begin();
         for(it=TrackerGainsPerRuns.begin(); it!=TrackerGainsPerRuns.end(); it++){
            if(it->first>currentRun){TrackerGains = &(itPrev->second); return;}//runs are ordered, so the previous iterator correspond to our run
            itPrev=it;
         }
         TrackerGains = &(itPrev->second); //just in case we go beyond the list of run for which we have a correciton
      }



      void LoadDeDxCalibration(string path){
         TrackerGainsPerRuns.clear();
         TrackerGains=NULL;

         TFile* InputFile = new TFile(path.c_str(), "r");
         TList* ObjList = InputFile->GetListOfKeys();
         for(int i=0;i<ObjList->GetSize();i++){
            TObject* tmp = GetObjectFromPath(InputFile,ObjList->At(i)->GetName(),false);
            if(tmp->InheritsFrom("TTree")){
               string dirName = ObjList->At(i)->GetName();
               unsigned int FirstRun, LastRun;  sscanf(dirName.c_str(), "Gains_%d_to_%d", &FirstRun, &LastRun);
               printf("Add a new gain starting at run %d\n", FirstRun);
               
               TTree* t1 = (TTree*) tmp;
               unsigned int  tree_DetId;   t1->SetBranchAddress("DetId"             ,&tree_DetId      );
               unsigned char tree_APVId;   t1->SetBranchAddress("APVId"             ,&tree_APVId      );
               double        tree_Gain;    t1->SetBranchAddress("Gain"              ,&tree_Gain       );
//               double        tree_PrevGain;t1->SetBranchAddress("PrevGain"          ,&tree_PrevGain   );

               TrackerGains = &TrackerGainsPerRuns[FirstRun];
               for (unsigned int ientry = 0; ientry < t1->GetEntries(); ientry++) {
                  t1->GetEntry(ientry);
                  (*TrackerGains)[tree_DetId<<3 | tree_APVId] = tree_Gain;
               }
            }
         }
         InputFile->Close();
         delete InputFile;
      }
};



TH3F* loadDeDxTemplate(string path, bool splitByModuleType){
   TFile* InputFile = new TFile(path.c_str());
   TH3F* DeDxMap_ = (TH3F*)GetObjectFromPath(InputFile, "Charge_Vs_Path");
   if(!DeDxMap_){printf("dEdx templates in file %s can't be open\n", path.c_str()); exit(0);}

   TH3F* Prob_ChargePath  = (TH3F*)(DeDxMap_->Clone("Prob_ChargePath")); 
   Prob_ChargePath->Reset();
   Prob_ChargePath->SetDirectory(0); 

   if(!splitByModuleType){
      Prob_ChargePath->RebinX(Prob_ChargePath->GetNbinsX()-1); // <-- do not include pixel in the inclusive
   }

   for(int i=0;i<=Prob_ChargePath->GetXaxis()->GetNbins()+1;i++){
      for(int j=0;j<=Prob_ChargePath->GetYaxis()->GetNbins()+1;j++){
         double Ni = 0;
         for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){Ni+=DeDxMap_->GetBinContent(i,j,k);}

         for(int k=0;k<=Prob_ChargePath->GetZaxis()->GetNbins()+1;k++){
            double tmp = 0;
            for(int l=0;l<=k;l++){ tmp+=DeDxMap_->GetBinContent(i,j,l);}

            if(Ni>0){
               Prob_ChargePath->SetBinContent (i, j, k, tmp/Ni);
            }else{
               Prob_ChargePath->SetBinContent (i, j, k, 0);
            }
         }
      }
   }
   InputFile->Close();
   return Prob_ChargePath;
}

#include "DataFormats/SiStripDetId/interface/SiStripDetId.h"

const double TkModGeomThickness[] = {1, 0.029000, 0.029000, 0.047000, 0.047000, 0.029000, 0.029000, 0.029000, 0.029000, 0.029000, 0.029000, 0.029000, 0.047000, 0.047000, 0.047000};
const double TkModGeomLength   [] = {1, 5.844250, 5.844250, 9.306700, 9.306700,  5.542900, 4.408000, 5.533000, 4.258000, 4.408000, 5.533000, 5.758000, 7.363125, 9.204400, 10.235775};
const double TkModGeomWidthB   [] = {1, 3.072000, 3.072000, 4.684800, 4.684800, 4.579445, 5.502407, 3.158509, 4.286362, 5.502407, 4.049915, 3.561019, 6.002559, 5.235483, 3.574395};
const double TkModGeomWidthT   [] = {1, 3.072000, 3.072000, 4.684800,  4.684800, 3.095721, 4.322593, 4.049915, 3.146580, 4.322593, 3.158509, 2.898798, 4.824683, 4.177638,  4.398049};

bool isHitInsideTkModule(const LocalPoint hitPos, const DetId& detid, const SiStripCluster* cluster=NULL){
   if(detid.subdetId()<3){return true;} //do nothing for pixel modules
   SiStripDetId SSdetId(detid);
   int moduleGeometry = SSdetId.moduleGeometry();

   //clean along the apv lines
   if(cluster && (cluster->firstStrip()%128 == 0 || (cluster->firstStrip() + cluster->amplitudes().size()%128==127))) return false;

   double nx, ny;
   if(moduleGeometry<=4){
      ny = hitPos.y() /  TkModGeomLength[moduleGeometry];
      nx = hitPos.x() /  TkModGeomWidthT[moduleGeometry];
   }else{
      double  offset = TkModGeomLength[moduleGeometry] * (TkModGeomWidthT[moduleGeometry]+TkModGeomWidthB[moduleGeometry]) / (TkModGeomWidthT[moduleGeometry]-TkModGeomWidthB[moduleGeometry]);  // check sign if GeomWidthT[moduleGeometry] < TkModGeomWidthB[moduleGeometry] !!! 
      double  tan_a = TkModGeomWidthT[moduleGeometry] / std::abs(offset + TkModGeomLength[moduleGeometry]);
      ny = hitPos.y() /  TkModGeomLength[moduleGeometry];
      nx = hitPos.x() / (tan_a*std::abs(hitPos.y()+offset));
   }

   // "blacklists" for the gaps and edges
   switch (moduleGeometry){
      case  0: return true;
      case  1: if (fabs(ny) > 0.96 || fabs(nx) > 0.98) return false; break;
      case  2: if (fabs(ny) > 0.97 || fabs(nx) > 0.99) return false; break;
      case  3: if (fabs(ny) > 0.98 || fabs(nx) > 0.98 || fabs(ny) <  0.04) return false; break;
      case  4: if (fabs(ny) > 0.98 || fabs(nx) > 0.98 || fabs(ny) <  0.04) return false; break;
      case  5: if (fabs(ny) > 0.98 || fabs(nx) > 0.98) return false; break;
      case  6: if (fabs(ny) > 0.98 || fabs(nx) > 0.99) return false; break;
      case  7: if (fabs(ny) > 0.97 || fabs(nx) > 0.98) return false; break;
      case  8: if (fabs(ny) > 0.98 || fabs(nx) > 0.99) return false; break;
      case  9: if (fabs(ny) > 0.98 || fabs(nx) > 0.99) return false; break;
      case 10: if (fabs(ny) > 0.97 || fabs(nx) > 0.99) return false; break;
      case 11: if (fabs(ny) > 0.97 || fabs(nx) > 0.99) return false; break;
      case 12: if (fabs(ny) > 0.98 || fabs(nx) > 0.99 || (-0.17 < ny && ny < -0.07)) return false; break;
      case 13: if (fabs(ny) > 0.97 || fabs(nx) > 0.99 || (-0.10 < ny && ny < -0.01)) return false; break;
      case 14: if (fabs(ny) > 0.95 || fabs(nx) > 0.98 || ( 0.01 < ny && ny <  0.12)) return false; break;
      default: std::cerr << "Unknown module geometry! Exiting!" << std::endl; exit (EXIT_FAILURE);
   }

   return true;
}

HitDeDxCollection getHitDeDx(const DeDxHitInfo* dedxHits, double* scaleFactors, std::unordered_map<unsigned int,double>* TrackerGains, bool correctFEDSat, int crossTalkInvAlgo){
     HitDeDxCollection toReturn;
     for(unsigned int h=0;h<dedxHits->size();h++){
        DetId detid(dedxHits->detId(h));  

        HitDeDx hit; 
        hit.subDet              = detid.subdetId();
        hit.passClusterCleaning = clusterCleaning(dedxHits->stripCluster(h), crossTalkInvAlgo);
        hit.isInside            = isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId()>=3?dedxHits->stripCluster(h):NULL);
        hit.isSat               = false;

        int ClusterCharge = dedxHits->charge(h);
        if(detid.subdetId()>=3){//for strip only
           const SiStripCluster* cluster = dedxHits->stripCluster(h);
           vector<int> amplitudes = convert(cluster->amplitudes());
           if (crossTalkInvAlgo) amplitudes = CrossTalkInv(amplitudes, 0.10, 0.04, true);
           int firstStrip = cluster->firstStrip();
           int prevAPV = -1;
           double gain = 1.0;

           ClusterCharge = 0;
           for(unsigned int s=0;s<amplitudes.size();s++){
              if(TrackerGains!=NULL){ //don't reload the gain if unnecessary  since map access are slow
                 int APV = (firstStrip+s)/128;
                 if(APV != prevAPV){gain = TrackerGains->at(detid.rawId()<<3 |APV); prevAPV=APV; }
              }

              int StripCharge =  amplitudes[s];
              if(StripCharge<254){
                 StripCharge=(int)(StripCharge/gain);
                 if(StripCharge>=1024){         StripCharge = 255;
                 }else if(StripCharge>=254){    StripCharge = 254;
                 }
              }

              if(StripCharge>=254){hit.isSat=true;}
              if(StripCharge>=255 && correctFEDSat){StripCharge=512;}
              ClusterCharge += StripCharge;
            } 
        }

        double scaleFactor = scaleFactors[0];
        if (detid.subdetId()<3) scaleFactor *= scaleFactors[1]; // add pixel scaling

        double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;
        hit.dx   = dedxHits->pathlength(h);
        hit.dedx = (scaleFactor*Norm*ClusterCharge) / hit.dx;

        toReturn.push_back(hit);
     } 
   return toReturn;
}



DeDxData computedEdx(const DeDxHitInfo* dedxHits, double* scaleFactors, TH3* templateHisto, bool usePixel, bool useClusterCleaning, bool reverseProb, bool useTruncated, std::unordered_map<unsigned int,double>* TrackerGains, bool useStrip, bool mustBeInside, size_t MaxStripNOM, bool correctFEDSat, int crossTalkInvAlgo, double dropLowerDeDxValue, dedxHIPEmulator* hipEmulator, double* dEdxErr, unsigned int pdgId){


//GenId==1092214, (gluino_uud) and
//GenId==1000612, (stop_anti-d).

     bool isStrangePdgId = false;
     if(pdgId==1092214|| pdgId==1000612) isStrangePdgId= true;

     if(!dedxHits) return DeDxData(-1, -1, -1);

//     if(templateHisto)usePixel=false; //never use pixel for discriminator

     std::vector<double> vect;
     std::vector<double> vectStrip;
     std::vector<double> vectPixel;

     unsigned int NSat=0;
     unsigned int SiStripNOM = 0;
     double lowerStripDeDx=1000;
     int lowerStripDeDxIndex=-1;
     for(unsigned int h=0;h<dedxHits->size();h++){
        DetId detid(dedxHits->detId(h));  
        if(!usePixel && detid.subdetId()<3)continue; // skip pixels
        if(!useStrip && detid.subdetId()>=3)continue; // skip strips        
        if(useClusterCleaning && !clusterCleaning(dedxHits->stripCluster(h), crossTalkInvAlgo))continue;
         //printStripCluster(stdout, dedxHits->stripCluster(h), dedxHits->detId(h));

        if(mustBeInside && !isHitInsideTkModule(dedxHits->pos(h), detid, detid.subdetId()>=3?dedxHits->stripCluster(h):NULL))continue;
        if(detid.subdetId()>=3 && ++SiStripNOM > MaxStripNOM) continue; // skip remaining strips, but not pixel

        int ClusterCharge = dedxHits->charge(h);

        if(detid.subdetId()>=3){//for strip only
           const SiStripCluster* cluster = dedxHits->stripCluster(h);
           vector<int> amplitudes = convert(cluster->amplitudes());
           if (crossTalkInvAlgo) amplitudes = CrossTalkInv(amplitudes, 0.10, 0.04, true);
           int firstStrip = cluster->firstStrip();
           int prevAPV = -1;
           double gain = 1.0;

           bool isSatCluster = false;
           ClusterCharge = 0;
           for(unsigned int s=0;s<amplitudes.size();s++){
              if(TrackerGains!=NULL){ //don't reload the gain if unnecessary  since map access are slow
                 int APV = (firstStrip+s)/128;
                 if(APV != prevAPV){gain = TrackerGains->at(detid.rawId()<<3 |APV); prevAPV=APV; }
              }

              int StripCharge =  amplitudes[s];
              if(StripCharge<254){
                 StripCharge=(int)(StripCharge/gain);
                 if(StripCharge>=1024){         StripCharge = 255;
                 }else if(StripCharge>=254){    StripCharge = 254;
                 }
              }

              if(StripCharge>=254){isSatCluster=true;}
              if(StripCharge>=255 && correctFEDSat){StripCharge=512;}
              ClusterCharge += StripCharge;
            } 
            if(isSatCluster)NSat++;
        }

        double scaleFactor = scaleFactors[0];
        if (detid.subdetId()<3) scaleFactor *= scaleFactors[1]; // add pixel scaling

        if(templateHisto){  //save discriminator probability
           double ChargeOverPathlength = scaleFactor*ClusterCharge/(dedxHits->pathlength(h)*10.0*(detid.subdetId()<3?265:1));
           if(isStrangePdgId) ChargeOverPathlength /= 2; 
//           if(fakeHIP && detid.subdetId()>=3 && rand()%1000<35)ChargeOverPathlength = ( 0.5 + ((rand()%15000)/10000.0) ) / (3.61e-06*265*10);
//           if(fakeHIP && detid.subdetId() <3 && rand()%1000<20)ChargeOverPathlength = ( 0.3 + ((rand()%12000)/10000.0) ) / (3.61e-06*265*10*265);

           int moduleGeometry = 0; // underflow for debug
           if (detid.subdetId()<3) moduleGeometry = 15; // 15 == pixel
           else {SiStripDetId SSdetId(detid); moduleGeometry = SSdetId.moduleGeometry();}
           int    BinX   = templateHisto->GetXaxis()->FindBin(moduleGeometry);
           int    BinY   = templateHisto->GetYaxis()->FindBin(dedxHits->pathlength(h)*10.0); //*10 because of cm-->mm
           int    BinZ   = templateHisto->GetZaxis()->FindBin(ChargeOverPathlength);
           double Prob   = templateHisto->GetBinContent(BinX,BinY,BinZ);
           //printf("%i %i %i  %f\n", BinX, BinY, BinZ, Prob);
           if(reverseProb)Prob = 1.0 - Prob;
           vect.push_back(Prob); //save probability
        }else{
           double Norm = (detid.subdetId()<3)?3.61e-06:3.61e-06*265;
           double ChargeOverPathlength = scaleFactor*Norm*ClusterCharge/dedxHits->pathlength(h);
           if(hipEmulator)ChargeOverPathlength = hipEmulator->fakeHIP(detid.subdetId(), ChargeOverPathlength);
    // mk change 
           if(isStrangePdgId) ChargeOverPathlength /= 2;

           vect.push_back(ChargeOverPathlength); //save charge
           if(detid.subdetId()< 3)vectPixel.push_back(ChargeOverPathlength);
           if(detid.subdetId()>=3)vectStrip.push_back(ChargeOverPathlength);
//           printf("%i - %f / %f = %f\n", h, scaleFactor*Norm*dedxHits->charge(h), dedxHits->pathlength(h), ChargeOverPathlength);
        }
     }

     if(dropLowerDeDxValue>0){
         std::vector <double> tmp (vect.size());
         std::copy (vect.begin(), vect.end(), tmp.begin());
         std::sort(tmp.begin(), tmp.end(), std::greater<double>() );
         int nTrunc = tmp.size()*dropLowerDeDxValue;
         vect.clear();
         for(unsigned int t=0;t+nTrunc<tmp.size();t++){vect.push_back(tmp[t]);}
     }

     double result;
     int size = vect.size();

     if(size>0){
        if(templateHisto){  //dEdx discriminator
           //Prod discriminator
           //result = 1;
           //for(int i=0;i<size;i++){
           //   if(vect[i]<=0.0001){result *= pow(0.0001 , 1.0/size);}
           //   else               {result *= pow(vect[i], 1.0/size);}
           //}

           //Ias discriminator
           result = 1.0/(12*size);
           std::sort(vect.begin(), vect.end(), std::less<double>() );
           for(int i=1;i<=size;i++){
              result += vect[i-1] * pow(vect[i-1] - ((2.0*i-1.0)/(2.0*size)),2);
           }
           result *= (3.0/size);
//mk printf("------> In Ias = %f\n",result);
        }else{  //dEdx estimator
           if(useTruncated){
              //truncated40 estimator
              std::sort(vect.begin(), vect.end(), std::less<double>() );              
              result=0;
              int nTrunc = size*0.40;
              for(int i = 0; i+nTrunc<size; i ++){
                 result+=vect[i];
              }
              result /= (size-nTrunc);
           }else{
              //harmonic2 estimator           
              result=0;
              double expo = -2;
	      if (dEdxErr) *dEdxErr = 0;
              for(int i = 0; i< size; i ++){
                 result+=pow(vect[i],expo);
		 if (dEdxErr) *dEdxErr += pow(vect[i],2*(expo-1))*pow(0.01,2);
              }
              result = pow(result/size,1./expo);
	      if (dEdxErr) *dEdxErr = result*result*result*sqrt(*dEdxErr)/size;
           }
           //mk if(isStrangePdgId) result /= 2;
           //mk if(isStrangePdgId) printf("------> In pdgId = %u\n",pdgId);
           //mk printf("------> In Ih = %f\n",result);
        }
     }else{
        result = -1;
     }
      //mk printf("Results =%f , NSat = %d  , size = %d\n",result, NSat, size);
     return DeDxData(result, NSat, size);
}



void printStripCluster(FILE* pFile, const SiStripCluster*   cluster, const DetId& DetId, bool crossTalkInvAlgo)
{
        if(!cluster)return;
        const vector<unsigned char>&  ampls       = cluster->amplitudes();

        int Charge=0;
        for(unsigned int i=0;i<ampls.size();i++){Charge+=ampls[i];}
        char clusterCleaningOutput = clusterCleaning(cluster, crossTalkInvAlgo) ? 'V' : 'X';

        fprintf(pFile,"DetId = %7i --> %4i = %3i ",DetId.rawId(),Charge,ampls[0]);
        for(unsigned int i=1;i<ampls.size();i++){
           fprintf(pFile,"%3i ",ampls[i]);
        }
        fprintf(pFile,"   %c\n", clusterCleaningOutput);
}




std::vector<int> convert(const vector<unsigned char>& input)
{
  std::vector<int> output;
  for(unsigned int i=0;i<input.size();i++){
        output.push_back((int)input[i]);
  }
  return output;
}


std::vector<int> CrossTalkInv(const std::vector<int>&  Q, const float x1, const float x2, bool way,float threshold,float thresholdSat) {
  const unsigned N=Q.size();
  std::vector<int> QII;
  std::vector<float> QI(N,0);
  Double_t a=1-2*x1-2*x2;
//  bool debugbool=false;
  TMatrix A(N,N);

//---  que pour 1 max bien net 
 if(Q.size()<2 || Q.size()>8){
	for (unsigned int i=0;i<Q.size();i++){
		QII.push_back((int) Q[i]);
  	}
	return QII;
  }
 if(way){ 
	  vector<int>::const_iterator mQ = max_element(Q.begin(), Q.end())	;
	  if(*mQ>253){
	 	 if(*mQ==255 && *(mQ-1)>253 && *(mQ+1)>253 ) return Q ;
	 	 if(*(mQ-1)>thresholdSat && *(mQ+1)>thresholdSat && *(mQ-1)<254 && *(mQ+1)<254 &&  abs(*(mQ-1) - *(mQ+1)) < 40 ){
		     QII.push_back((10*(*(mQ-1))+10*(*(mQ+1)))/2); return QII;}
	  }
  }
//---

  for(unsigned int i=0; i<N; i++) {
        A(i,i) =a;
        if(i<N-1){ A(i+1,i)=x1;A(i,i+1)=x1;}
        else continue; 
        if(i<N-2){ A(i+2,i)=x2;A(i,i+2)=x2;}
  }

  if(N==1) A(0,0)=1/a;
  else  A.InvertFast();

  for(unsigned int i=0; i<N; i++) {
        for(unsigned int j=0; j<N; j++) {
        QI[i]+=A(i,j)*(float)Q[j];
        }
  }

 for (unsigned int i=0;i<QI.size();i++){
	if(QI[i]<threshold) QI[i]=0; 
	QII.push_back((int) QI[i]);
  }

return QII;
}


bool clusterCleaning(const SiStripCluster*   cluster,  int crosstalkInv, uint8_t * exitCode)
{
   if(!cluster) return true;
   vector<int>  ampls = convert(cluster->amplitudes());
   if(crosstalkInv==1)ampls = CrossTalkInv(ampls,0.10,0.04, true);
      

  // ----------------  COMPTAGE DU NOMBRE DE MAXIMA   --------------------------
  //----------------------------------------------------------------------------
         Int_t NofMax=0; Int_t recur255=1; Int_t recur254=1;
         bool MaxOnStart=false;bool MaxInMiddle=false, MaxOnEnd =false;
         Int_t MaxPos=0;
        // Début avec max
         if(ampls.size()!=1 && ((ampls[0]>ampls[1])
            || (ampls.size()>2 && ampls[0]==ampls[1] && ampls[1]>ampls[2] && ampls[0]!=254 && ampls[0]!=255) 
            || (ampls.size()==2 && ampls[0]==ampls[1] && ampls[0]!=254 && ampls[0]!=255)) ){
          NofMax=NofMax+1;  MaxOnStart=true;  }

        // Maximum entouré
         if(ampls.size()>2){
          for (unsigned int i =1; i < ampls.size()-1; i++) {
                if( (ampls[i]>ampls[i-1] && ampls[i]>ampls[i+1]) 
                    || (ampls.size()>3 && i>0 && i<ampls.size()-2 && ampls[i]==ampls[i+1] && ampls[i]>ampls[i-1] && ampls[i]>ampls[i+2] && ampls[i]!=254 && ampls[i]!=255) ){ 
                 NofMax=NofMax+1; MaxInMiddle=true;  MaxPos=i; 
                }
                if(ampls[i]==255 && ampls[i]==ampls[i-1]) {
                        recur255=recur255+1;
                        MaxPos=i-(recur255/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
                if(ampls[i]==254 && ampls[i]==ampls[i-1]) {
                        recur254=recur254+1;
                        MaxPos=i-(recur254/2);
                        if(ampls[i]>ampls[i+1]){NofMax=NofMax+1;MaxInMiddle=true;}
                }
            }
         }
        // Fin avec un max
         if(ampls.size()>1){
          if(ampls[ampls.size()-1]>ampls[ampls.size()-2]
             || (ampls.size()>2 && ampls[ampls.size()-1]==ampls[ampls.size()-2] && ampls[ampls.size()-2]>ampls[ampls.size()-3] ) 
             ||  ampls[ampls.size()-1]==255){
           NofMax=NofMax+1;  MaxOnEnd=true;   }
         }
        // Si une seule strip touchée
        if(ampls.size()==1){    NofMax=1;}



  // ---  SELECTION EN FONCTION DE LA FORME POUR LES MAXIMA UNIQUES ---------
  //------------------------------------------------------------------------
//  
//               ____
//              |    |____
//          ____|    |    |
//         |    |    |    |____
//     ____|    |    |    |    |
//    |    |    |    |    |    |____
//  __|____|____|____|____|____|____|__
//    C_Mnn C_Mn C_M  C_D  C_Dn C_Dnn
//  
//   bool shapetest=true;
   bool shapecdtn=false;
   if (exitCode) *exitCode = 255;

      if(crosstalkInv==1){
        if(NofMax==1){shapecdtn=true; if (exitCode) *exitCode=0;}
        return shapecdtn;
      }

//      Float_t C_M;    Float_t C_D;    Float_t C_Mn;   Float_t C_Dn;   Float_t C_Mnn;  Float_t C_Dnn;
        Float_t C_M=0.0;        Float_t C_D=0.0;        Float_t C_Mn=10000;     Float_t C_Dn=10000;     Float_t C_Mnn=10000;    Float_t C_Dnn=10000;
        Int_t CDPos;
        Float_t coeff1=1.7;     Float_t coeff2=2.0;
        Float_t coeffn=0.10;    Float_t coeffnn=0.02; Float_t noise=4.0;

        if(NofMax==1){

                if(MaxOnStart==true){
                        C_M=(Float_t)ampls[0]; C_D=(Float_t)ampls[1];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[2] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=2;}
                                else if(ampls.size()>3){ C_Dn=(Float_t)ampls[2];  C_Dnn=(Float_t)ampls[3] ;
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){
                                                         shapecdtn=true;} else if (exitCode) *exitCode=3;
                                }
                }

                if(MaxOnEnd==true){
                        C_M=(Float_t)ampls[ampls.size()-1]; C_D=(Float_t)ampls[ampls.size()-2];
                                if(ampls.size()<3) shapecdtn=true ;
                                else if(ampls.size()==3){C_Dn=(Float_t)ampls[0] ; if(C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255) shapecdtn=true; else if (exitCode) *exitCode=4;}
                                else if(ampls.size()>3){C_Dn=(Float_t)ampls[ampls.size()-3] ; C_Dnn=(Float_t)ampls[ampls.size()-4] ; 
                                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise){ 
                                                         shapecdtn=true;} else if (exitCode) *exitCode=5;
                                }
                }

                if(MaxInMiddle==true){
                        C_M=(Float_t)ampls[MaxPos];
                        int LeftOfMaxPos=MaxPos-1;if(LeftOfMaxPos<=0)LeftOfMaxPos=0;
                        int RightOfMaxPos=MaxPos+1;if(RightOfMaxPos>=(int)ampls.size())RightOfMaxPos=ampls.size()-1;
                        //int after = RightOfMaxPos; int before = LeftOfMaxPos; if (after>=(int)ampls.size() ||  before<0)  std::cout<<"invalid read MaxPos:"<<MaxPos <<"size:"<<ampls.size() <<std::endl; 
                        if(ampls[LeftOfMaxPos]<ampls[RightOfMaxPos]){ C_D=(Float_t)ampls[RightOfMaxPos]; C_Mn=(Float_t)ampls[LeftOfMaxPos];CDPos=RightOfMaxPos;} else{ C_D=(Float_t)ampls[LeftOfMaxPos]; C_Mn=(Float_t)ampls[RightOfMaxPos];CDPos=LeftOfMaxPos;}
                        if(C_Mn<coeff1*coeffn*C_M+coeff2*coeffnn*C_D+2*noise || C_M==255){ 
                                if(ampls.size()==3) shapecdtn=true ;
                                else if(ampls.size()>3){
                                        if(CDPos>MaxPos){
                                                if(ampls.size()-CDPos-1==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1==1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=0;
                                                }
                                                if(ampls.size()-CDPos-1>1){
                                                        C_Dn=(Float_t)ampls[CDPos+1];
                                                        C_Dnn=(Float_t)ampls[CDPos+2];
                                                }
                                                if(MaxPos>=2){
                                                        C_Mnn=(Float_t)ampls[MaxPos-2];
                                                }
                                                else if(MaxPos<2) C_Mnn=0;
                                        }
                                        if(CDPos<MaxPos){
                                                if(CDPos==0){
                                                        C_Dn=0; C_Dnn=0;
                                                }
                                                if(CDPos==1){
                                                        C_Dn=(Float_t)ampls[0];
                                                        C_Dnn=0;
                                                }
                                                if(CDPos>1){
                                                        C_Dn=(Float_t)ampls[CDPos-1];
                                                        C_Dnn=(Float_t)ampls[CDPos-2];
                                                }
                                                if(ampls.size()-LeftOfMaxPos>1 && MaxPos+2<(int)(ampls.size())-1){
                                                        C_Mnn=(Float_t)ampls[MaxPos+2];
                                                }else C_Mnn=0;                                                  
                                        }
                                        if((C_Dn<=coeff1*coeffn*C_D+coeff2*coeffnn*C_M+2*noise || C_D==255)
                                           && C_Mnn<=coeff1*coeffn*C_Mn+coeff2*coeffnn*C_M+2*noise
                                           && C_Dnn<=coeff1*coeffn*C_Dn+coeff2*coeffnn*C_D+2*noise) {
                                                shapecdtn=true;
                                        }

                                }
                        } else if (exitCode) *exitCode=6;
                }
        }
        else if (NofMax>1 && exitCode) *exitCode = 1; // more than one maximum
        if(ampls.size()==1){shapecdtn=true;}
        if(shapecdtn && exitCode) *exitCode=0;

   return shapecdtn;
}

class L1BugEmulator{
   public:
      L1BugEmulator(){
         double Weight    = preTrackingChangeL1IntLumi/IntegratedLuminosity13TeV16;
         double pTBins[]  = {0.0, 52, 80, 120, 200, 500, 999999};
         double EtaBins[] = {0.0, 0.9, 1.2, 2.1, 2.4, 999999};

         unsigned int NumPtBins  = static_cast<unsigned int> (sizeof(pTBins)/sizeof(double))  -1;
         unsigned int NumEtaBins = static_cast<unsigned int> (sizeof(EtaBins)/sizeof(double)) -1;

         // first period
         h1 = new TH2D ("Inefficiency_Run273158_274093", "Inefficiency_Run273158_274093", NumPtBins, pTBins, NumEtaBins, EtaBins);
         // eta<0.9 -> pass
         h1->SetBinContent(h1->FindBin( 25.0, 0.5), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 0.5), 0.92);
         h1->SetBinContent(h1->FindBin( 85.0, 0.5), 0.92);
         h1->SetBinContent(h1->FindBin(125.0, 0.5), 0.92);
         h1->SetBinContent(h1->FindBin(255.0, 0.5), 0.88);
         h1->SetBinContent(h1->FindBin(555.0, 0.5), 0.86);
         // 0.9<eta<1.2 -> depends on pT, pT>500 GeV get rejected
         h1->SetBinContent(h1->FindBin( 25.0, 1.0), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 1.0), 0.70);
         h1->SetBinContent(h1->FindBin( 85.0, 1.0), 0.68);
         h1->SetBinContent(h1->FindBin(125.0, 1.0), 0.65);
         h1->SetBinContent(h1->FindBin(255.0, 1.0), 0.68);
         h1->SetBinContent(h1->FindBin(555.0, 1.0), 0.00);
         // 1.2<eta<2.1
         h1->SetBinContent(h1->FindBin( 25.0, 1.5), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 1.5), 0.86);
         h1->SetBinContent(h1->FindBin( 85.0, 1.5), 0.87);
         h1->SetBinContent(h1->FindBin(125.0, 1.5), 0.78);
         h1->SetBinContent(h1->FindBin(255.0, 1.5), 0.69);
         h1->SetBinContent(h1->FindBin(555.0, 1.5), 0.00);
         // 2.1<eta<2.4
         h1->SetBinContent(h1->FindBin( 25.0, 2.2), 0.00); // we do not trigger on muons with pT<50
         h1->SetBinContent(h1->FindBin( 55.0, 2.2), 0.75);
         h1->SetBinContent(h1->FindBin( 85.0, 2.2), 0.80);
         h1->SetBinContent(h1->FindBin(125.0, 2.2), 0.68);
         h1->SetBinContent(h1->FindBin(255.0, 2.2), 0.00);
         h1->SetBinContent(h1->FindBin(555.0, 2.2), 0.00);
         // 2.4<eta is impossible -> remove muon
         h1->SetBinContent(h1->FindBin( 25.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin( 55.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin( 85.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin(125.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin(255.0, 3.2), 0.00);
         h1->SetBinContent(h1->FindBin(555.0, 3.2), 0.00);

         // second period
         h2 = new TH2D ("Inefficiency_274094_", "Inefficiency_Run_274094_", NumPtBins, pTBins, NumEtaBins, EtaBins);
         // eta<0.9 -> pass
         h2->SetBinContent(h2->FindBin( 25.0, 0.5), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin( 85.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin(125.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin(255.0, 0.5), 0.92);
         h2->SetBinContent(h2->FindBin(555.0, 0.5), 0.92);
         // 0.9<eta<1.2 -> depends on pT, pT>500 GeV get rejected
         h2->SetBinContent(h2->FindBin( 25.0, 1.0), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 1.0), 0.88);
         h2->SetBinContent(h2->FindBin( 85.0, 1.0), 0.88);
         h2->SetBinContent(h2->FindBin(125.0, 1.0), 0.89);
         h2->SetBinContent(h2->FindBin(255.0, 1.0), 0.84);
         h2->SetBinContent(h2->FindBin(555.0, 1.0), 0.00);
         // 1.2<eta<2.1
         h2->SetBinContent(h2->FindBin( 25.0, 1.5), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 1.5), 0.84);
         h2->SetBinContent(h2->FindBin( 85.0, 1.5), 0.84);
         h2->SetBinContent(h2->FindBin(125.0, 1.5), 0.76);
         h2->SetBinContent(h2->FindBin(255.0, 1.5), 0.63);
         h2->SetBinContent(h2->FindBin(555.0, 1.5), 0.00);
         // 2.1<eta<2.4
         h2->SetBinContent(h2->FindBin( 25.0, 2.2), 0.00); // we do not trigger on muons with pT<50
         h2->SetBinContent(h2->FindBin( 55.0, 2.2), 0.72);
         h2->SetBinContent(h2->FindBin( 85.0, 2.2), 0.76);
         h2->SetBinContent(h2->FindBin(125.0, 2.2), 0.78);
         h2->SetBinContent(h2->FindBin(255.0, 2.2), 0.66);
         h2->SetBinContent(h2->FindBin(555.0, 2.2), 0.00);
         // 2.4<eta is impossible -> remove muon
         h2->SetBinContent(h2->FindBin( 25.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin( 55.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin( 85.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin(125.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin(255.0, 3.2), 0.00);
         h2->SetBinContent(h2->FindBin(555.0, 3.2), 0.00);
      }

      ~L1BugEmulator(){
//         delete h1;
//         delete h2;
      }

      bool PassesL1Inefficiency(double pT, double Eta){
         if ((rand()%9999*1.0)/10000 < Weight)
            return (((rand()%9999)*1.0/10000) < h1->GetBinContent(h1->FindBin(pT, Eta)));
         else
            return (((rand()%9999)*1.0/10000) < h2->GetBinContent(h2->FindBin(pT, Eta)));
      }

   private:
      TH2D* h1;
      TH2D* h2;
      double Weight;
};

class HIPTrackLossEmulator{
  private:
     TH1D*  h;
     double lossRate;
  public:
     HIPTrackLossEmulator(){
        double NVertexBins[] = {0, 5, 10, 15, 22, 28, 35, 9999};
        h = new TH1D("TrackLoss_vs_NVertex", "TrackLoss_vs_NVertex", static_cast<unsigned int> (sizeof(NVertexBins)/sizeof(double)) -1, NVertexBins);
        h->SetBinContent(h->FindBin( 2.5), 1.000);
        h->SetBinContent(h->FindBin( 7.5), 0.998);
        h->SetBinContent(h->FindBin(12.5), 0.995);
        h->SetBinContent(h->FindBin(17.5), 0.992);
        h->SetBinContent(h->FindBin(24.5), 0.985);
        h->SetBinContent(h->FindBin(30.0), 0.980);
        h->SetBinContent(h->FindBin(37.0), 0.970);

        lossRate = 0.0;
     }

     ~HIPTrackLossEmulator(){
//        delete h;
     }

     void SetHIPTrackLossRate(fwlite::ChainEvent& ev){
        fwlite::Handle< std::vector<reco::Vertex> > vertexCollHandle;
        vertexCollHandle.getByLabel(ev,"offlinePrimaryVertices");
        if(!vertexCollHandle.isValid()){ lossRate = 1.0; return;}
        const std::vector<reco::Vertex>& vertexColl = *vertexCollHandle;
        if(vertexColl.size()<1){lossRate = 1.0; return;}
        else lossRate = 1.0 - h->GetBinContent(h->FindBin(vertexColl.size()));
     }

     bool TrackSurvivesHIPInefficiency(){
        return (((rand()%999999)*1.0/1000000) > lossRate);
     }
};

void printClusterCleaningMessage (uint8_t exitCode){
   switch (exitCode){
      case 0:  std::cout << "All went well"                     << std::endl; break;
      case 1:  std::cout << "More than one maximum"             << std::endl; break;
      case 2:  std::cout << "MFirst; CSize=3; CDn too big"      << std::endl; break;
      case 3:  std::cout << "MFirst; CSize>3; CDn||CDnn too big"<< std::endl; break;
      case 4:  std::cout << "MEnd; CSize=3; CDn too big"        << std::endl; break;
      case 5:  std::cout << "MEnd; CSize>3; CDn||CDnn too big"  << std::endl; break;
      case 6:  std::cout << "MMid; Sides are too big"           << std::endl; break;
      case 255:std::cout << "Failed all shape tests"            << std::endl; break;
      default: std::cout << "Unknown exit code!"<< std::endl;
   }
}
#endif
