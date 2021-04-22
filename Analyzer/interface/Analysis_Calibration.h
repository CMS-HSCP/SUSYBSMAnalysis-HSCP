#ifndef SUSYBSMAnalysis_Analyzer_Analysis_Calibration_h
#define SUSYBSMAnalysis_Analyzer_Analysis_Calibration_h

TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy=false)
{
   size_t pos = Path.find("/");
   if(pos < 256){
      std::string firstPart = Path.substr(0,pos);
      std::string endPart   = Path.substr(pos+1,Path.length());
      TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
      if(TMP!=NULL)return GetObjectFromPath(TMP,endPart,GetACopy);

      printf("ObjectNotFound: %s::%s\n",File->GetName(), Path.c_str());
      return NULL;
   }else{
      if(GetACopy){
         return (File->Get(Path.c_str()))->Clone();
      }else{
         return File->Get(Path.c_str());
      }
   }
}

// similar to the above code
TObject* GetObjectFromPath(TDirectory* Container, TDirectory* File, std::string Path, bool GetACopy=false){
   TObject* toreturn = GetObjectFromPath(File,Path,GetACopy);
   if(TH1* th1 = dynamic_cast<TH1*>(toreturn))th1->SetDirectory(Container);
   return toreturn;
}


TH1D* GetProjectionFromPath(TDirectory* File, std::string Path, int CutIndex, std::string Name){
      TH2D* tmp = (TH2D*)GetObjectFromPath(File, Path, false);
      if(!tmp)return NULL;
      return tmp->ProjectionY(Name.c_str()   ,CutIndex+1,CutIndex+1,"o");
}

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



      void LoadDeDxCalibration(std::string path){
         TrackerGainsPerRuns.clear();
         TrackerGains=NULL;

         TFile* InputFile = new TFile(path.c_str(), "r");
         TList* ObjList = InputFile->GetListOfKeys();
         for(int i=0;i<ObjList->GetSize();i++){
            TObject* tmp = GetObjectFromPath(InputFile,ObjList->At(i)->GetName(),false);
            if(tmp->InheritsFrom("TTree")){
               std::string dirName = ObjList->At(i)->GetName();
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



TH3F* loadDeDxTemplate(std::string path, bool splitByModuleType){
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
#endif