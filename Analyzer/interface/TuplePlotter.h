#include "../interface/Tuple.h"

struct Tuple;

class TuplePlotter{
   public:
      TuplePlotter(int TypeMode, std::string SavePath){TypeMode_ = TypeMode; SavePath_ = SavePath;};
      ~TuplePlotter();

      void getObjects(Tuple* &tuple, TDirectory* dir);
      void draw(Tuple* tuple, TDirectory* dir, std::string LegendTitle, unsigned int CutIndex);
   private:
      int TypeMode_;
      std::string SavePath_;
};

//=============================================================
//
//     get objects: histos, trees
//
//=============================================================


void TuplePlotter::getObjects(Tuple* &tuple, TDirectory* dir){
   
   
   std::string Name;

   Name = "IntLumi";  tuple->IntLumi = (TProfile*)GetObjectFromPath(dir, Name); 
   Name = "TotalE";   tuple->TotalE  = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "TotalEPU"; tuple->TotalEPU= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "TotalTE";  tuple->TotalTE = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "Total";    tuple->Total   = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "V3D";      tuple->V3D     = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "Chi2";     tuple->Chi2    = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "Qual";     tuple->Qual    = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "TNOH";     tuple->TNOH    = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "TNOM";     tuple->TNOM    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "nDof";     tuple->nDof    = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "tofError"; tuple->tofError= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Pterr";    tuple->Pterr   = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "TIsol";    tuple->TIsol   = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "EIsol";    tuple->EIsol   = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "SumpTOverpT";  tuple->SumpTOverpT = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "MPt";      tuple->MPt     = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "MI";       tuple->MI      = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "MTOF";     tuple->MTOF    = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "Pt";       tuple->Pt      = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "I";        tuple->I       = (TH1F*)GetObjectFromPath(dir, Name);     
   Name = "TOF";      tuple->TOF     = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "HSCPE";    tuple->HSCPE   = (TH1F*)GetObjectFromPath(dir, Name);

   Name = "NVTrack";  tuple->NVTrack = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Stations"; tuple->Stations= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Dxy";      tuple->Dxy     = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Dz";       tuple->Dz      = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "SegSep";   tuple->SegSep  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "FailDz";   tuple->FailDz  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Basic";    tuple->Basic   = (TH1F*)GetObjectFromPath(dir, Name);

   Name = "HSCPE_SystP";    tuple->HSCPE_SystP  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "HSCPE_SystI";    tuple->HSCPE_SystI  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "HSCPE_SystM";    tuple->HSCPE_SystM  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "HSCPE_SystT";    tuple->HSCPE_SystT  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "HSCPE_SystPU";   tuple->HSCPE_SystPU = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "HSCPE_SystHUp";  tuple->HSCPE_SystHUp  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "HSCPE_SystHDown";tuple->HSCPE_SystHDown  = (TH1F*)GetObjectFromPath(dir, Name);

   
   Name = "Mass";     tuple->Mass     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF";  tuple->MassTOF  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb"; tuple->MassComb = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass";     tuple->MaxEventMass     = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "Mass_SystP";     tuple->Mass_SystP     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF_SystP";  tuple->MassTOF_SystP  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_SystP"; tuple->MassComb_SystP = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass_SystP";     tuple->MaxEventMass_SystP = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "Mass_SystI";     tuple->Mass_SystI     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF_SystI";  tuple->MassTOF_SystI  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_SystI"; tuple->MassComb_SystI = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass_SystI";     tuple->MaxEventMass_SystI = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "Mass_SystM";     tuple->Mass_SystM     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF_SystM";  tuple->MassTOF_SystM  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_SystM"; tuple->MassComb_SystM = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass_SystM";     tuple->MaxEventMass_SystM = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "Mass_SystT";     tuple->Mass_SystT     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF_SystT";  tuple->MassTOF_SystT  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_SystT"; tuple->MassComb_SystT = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass_SystT";     (TH2F*)GetObjectFromPath(dir, Name);

   Name = "Mass_SystPU";    tuple->Mass_SystPU     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF_SystPU"; tuple->MassTOF_SystPU  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_SystPU";tuple->MassComb_SystPU = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass_SystPU";  tuple->MaxEventMass_SystPU = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "Mass_SystHUp";     tuple->Mass_SystHUp     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF_SystH";  tuple->MassTOF_SystH  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_SystHUp"; tuple->MassComb_SystHUp = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass_SystHUp";     tuple->MaxEventMass_SystHUp = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "Mass_SystHDown";     tuple->Mass_SystHDown     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_SystHDown"; tuple->MassComb_SystHDown = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MaxEventMass_SystHDown";     tuple->MaxEventMass_SystHDown = (TH2F*)GetObjectFromPath(dir, Name);


   Name = "Mass_Flip";     tuple->Mass_Flip     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassTOF_Flip";  tuple->MassTOF_Flip  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "MassComb_Flip"; tuple->MassComb_Flip = (TH2F*)GetObjectFromPath(dir, Name);

   //UNUSED//if(SkipSelectionPlot)return;
   
   Name = "Gen_DecayLength"  ; tuple->Gen_DecayLength  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_Gen"         ; tuple->Beta_Gen         = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_GenChaged"   ; tuple->Beta_GenCharged  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_Triggered"   ; tuple->Beta_Triggered   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_Matched"     ; tuple->Beta_Matched     = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_PreselectedA"; tuple->Beta_PreselectedA= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_PreselectedB"; tuple->Beta_PreselectedB= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_PreselectedC"; tuple->Beta_PreselectedC= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "Beta_SelectedP"   ; tuple->Beta_SelectedP   = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "Beta_SelectedI"   ; tuple->Beta_SelectedI   = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "Beta_SelectedT"   ; tuple->Beta_SelectedT   = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "BS_V3D"  ; tuple->BS_V3D   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Chi2" ; tuple->BS_Chi2  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Qual" ; tuple->BS_Qual  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOH" ; tuple->BS_TNOH  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOH_PUA" ; tuple->BS_TNOH_PUA  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOH_PUB" ; tuple->BS_TNOH_PUB  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOHFraction" ; tuple->BS_TNOHFraction  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOPH" ; tuple->BS_TNOPH  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOHFractionTillLast" ; tuple->BS_TNOHFractionTillLast  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOMHTillLast" ; tuple->BS_TNOMHTillLast  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Eta" ; tuple->BS_Eta  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOM" ; tuple->BS_TNOM  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOM_PUA" ; tuple->BS_TNOM_PUA  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TNOM_PUB" ; tuple->BS_TNOM_PUB  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_nDof" ; tuple->BS_nDof  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOFError" ; tuple->BS_TOFError  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_PtErr"; tuple->BS_Pterr = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_MPt"  ; tuple->BS_MPt   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_MIs"  ; tuple->BS_MIs   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_MIm"  ; tuple->BS_MIm   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_MTOF" ; tuple->BS_MTOF  = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TIsol"; tuple->BS_TIsol = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_EIsol"; tuple->BS_EIsol = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_SumpTOverpT";  tuple->BS_SumpTOverpT = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_LastHitDXY"    ; tuple->BS_LastHitDXY     = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_LastHitD3D"    ; tuple->BS_LastHitD3D     = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_P"    ; tuple->BS_P     = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt"   ; tuple->BS_Pt    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_PUA"   ; tuple->BS_Pt_PUA    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_PUB"   ; tuple->BS_Pt_PUB    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_Cosmic"   ; tuple->BS_Pt_Cosmic    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_DT"   ; tuple->BS_Pt_DT    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_CSC"   ; tuple->BS_Pt_CSC    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Is"   ; tuple->BS_Is    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Is_PUA"   ; tuple->BS_Is_PUA    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Is_PUB"   ; tuple->BS_Is_PUB    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Is_Cosmic"   ; tuple->BS_Is_Cosmic    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Im"   ; tuple->BS_Im    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Im_PUA"   ; tuple->BS_Im_PUA    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Im_PUB"   ; tuple->BS_Im_PUB    = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF"  ; tuple->BS_TOF   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF_PUA"  ; tuple->BS_TOF_PUA   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF_PUB"  ; tuple->BS_TOF_PUB   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF_DT"  ; tuple->BS_TOF_DT   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF_CSC"  ; tuple->BS_TOF_CSC   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_dR_NVTrack"  ; tuple->BS_dR_NVTrack = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_MatchedStations"  ; tuple->BS_MatchedStations= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_InnerInvPtDiff"  ; tuple->BS_InnerInvPtDiff = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Phi"  ; tuple->BS_Phi = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TimeAtIP"  ; tuple->BS_TimeAtIP = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_OpenAngle"  ; tuple->BS_OpenAngle = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_OpenAngle_Cosmic"  ; tuple->BS_OpenAngle_Cosmic = (TH1F*)GetObjectFromPath(dir, Name);

   Name = "BS_NVertex";  tuple->BS_NVertex = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_NVertex_NoEventWeight";    (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_PV"  ; tuple->BS_PV = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_PV_NoEventWeight"  ; tuple->BS_PV_NoEventWeight = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_NOMoNOHvsPV"  ; tuple->BS_NOMoNOHvsPV = (TProfile*)GetObjectFromPath(dir, Name);
   Name = "BS_dzAll";      tuple->BS_dzAll = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_dxyAll";     tuple->BS_dxyAll = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_dzMinv3d";   tuple->BS_dzMinv3d = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_dxyMinv3d";  tuple->BS_dxyMinv3d = (TH1F*)GetObjectFromPath(dir, Name);

   Name = "BS_SegSep"  ; tuple->BS_SegSep= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_SegMinEtaSep"  ; tuple->BS_SegMinEtaSep= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_SegMinPhiSep"  ; tuple->BS_SegMinPhiSep= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_SegMinEtaSep_FailDz"  ; tuple->BS_SegMinEtaSep_FailDz= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_SegMinEtaSep_PassDz"  ; tuple->BS_SegMinEtaSep_PassDz= (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Dz_FailSep"; tuple->BS_Dz_FailSep   = (TH1F*)GetObjectFromPath(dir, Name);

   Name = "BS_Dxy"; tuple->BS_Dxy   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Dxy_Cosmic"; tuple->BS_Dxy_Cosmic   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Dz"; tuple->BS_Dz   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Dz_Cosmic"; tuple->BS_Dz_Cosmic   = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Dz_CSC"; tuple->BS_Dz_CSC = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Dz_DT"; tuple->BS_Dz_DT=(TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_FailDz"; tuple->BS_Pt_FailDz = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_FailDz_DT"; tuple->BS_Pt_FailDz_DT = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_Pt_FailDz_CSC"; tuple->BS_Pt_FailDz_CSC = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF_FailDz"; tuple->BS_TOF_FailDz = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF_FailDz_DT"; tuple->BS_TOF_FailDz_DT = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOF_FailDz_CSC"; tuple->BS_TOF_FailDz_CSC = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "genrecopT"; tuple->genrecopT = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "genlevelpT";    tuple->genlevelpT = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "genleveleta";   tuple->genleveleta = (TH1F*)GetObjectFromPath(dir, Name);
   Name = "genlevelbeta";  tuple->genlevelbeta = (TH1F*)GetObjectFromPath(dir, Name);

   //Initialize histograms for number of bins.
   for(int i=0; i<6; i++) {
     char Suffix[1024];
     sprintf(Suffix,"_%i",i);
     Name = "BS_Pt_Binned"; Name.append(Suffix); tuple->BS_Pt_Binned[std::to_string(i)] = (TH1F*)GetObjectFromPath(dir, Name);
     Name = "BS_TOF_Binned"; Name.append(Suffix); tuple->BS_TOF_Binned[std::to_string(i)] = (TH1F*)GetObjectFromPath(dir, Name);
   }

   Name = "AS_Eta_RegionA" ; tuple->AS_Eta_RegionA  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Eta_RegionB" ; tuple->AS_Eta_RegionB  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Eta_RegionC" ; tuple->AS_Eta_RegionC  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Eta_RegionD" ; tuple->AS_Eta_RegionD  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Eta_RegionE" ; tuple->AS_Eta_RegionE  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Eta_RegionF" ; tuple->AS_Eta_RegionF  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Eta_RegionG" ; tuple->AS_Eta_RegionG  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Eta_RegionH" ; tuple->AS_Eta_RegionH  = (TH2F*)GetObjectFromPath(dir, Name);

   Name = "AS_P"    ; tuple->AS_P     = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Pt"   ; tuple->AS_Pt    = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Is"   ; tuple->AS_Is    = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_Im"   ; tuple->AS_Im    = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "AS_TOF"  ; tuple->AS_TOF   = (TH2F*)GetObjectFromPath(dir, Name);


   Name = "BS_EtaIs"; tuple->BS_EtaIs = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_EtaIm"; tuple->BS_EtaIm = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_EtaP" ; tuple->BS_EtaP  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_EtaPt"; tuple->BS_EtaPt = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_EtaTOF" ; tuple->BS_EtaTOF = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_EtaNBH" ; tuple->BS_EtaNBH = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_EtaDz"; tuple->BS_EtaDz  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_PIs"  ; tuple->BS_PIs   = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_PImHD"; tuple->BS_PImHD = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_PIm"  ; tuple->BS_PIm   = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_PtIs" ; tuple->BS_PtIs  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_PtIm" ; tuple->BS_PtIm  = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_PtTOF" ; tuple->BS_PtTOF= (TH2F*)GetObjectFromPath(dir, Name);
   //   Name = "BS_TOFIs"; tuple->BS_TOFIs = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOFIs"; tuple->BS_TOFIs = (TH2F*)GetObjectFromPath(dir, Name);
   //   Name = "BS_TOFIm"; BS_TOFIm = (TH2F*)GetObjectFromPath(dir, Name);
   Name = "BS_TOFIm"; tuple->BS_TOFIm = (TH2F*)GetObjectFromPath(dir, Name);

   //   Name = "AS_EtaIs"; AS_EtaIs = (TH3F*)GetObjectFromPath(dir, Name);
   //   Name = "AS_EtaIm"; AS_EtaIm = (TH3F*)GetObjectFromPath(dir, Name);
   //   Name = "AS_EtaP" ; AS_EtaP  = (TH3F*)GetObjectFromPath(dir, Name);
   //   Name = "AS_EtaPt"; AS_EtaPt = (TH3F*)GetObjectFromPath(dir, Name);
   //   Name = "AS_EtaTOF"; AS_EtaTOF = (TH3F*)GetObjectFromPath(dir, Name);
   Name = "AS_PIs"  ; tuple->AS_PIs   = (TH3F*)GetObjectFromPath(dir, Name);
   Name = "AS_PIm"  ; tuple->AS_PIm   = (TH3F*)GetObjectFromPath(dir, Name);
   Name = "AS_PtIs" ; tuple->AS_PtIs  = (TH3F*)GetObjectFromPath(dir, Name);
   Name = "AS_PtIm" ; tuple->AS_PtIm  = (TH3F*)GetObjectFromPath(dir, Name);
   Name = "AS_TOFIs"; tuple->AS_TOFIs = (TH3F*)GetObjectFromPath(dir, Name);
   Name = "AS_TOFIm"; tuple->AS_TOFIm = (TH3F*)GetObjectFromPath(dir, Name);

   
}

//=============================================================
//
//  Draw all plots that are not meant for comparison with other samples (mostly 2D plots that can't be superimposed)
//
//=============================================================
void TuplePlotter::draw(Tuple* tuple, TDirectory* dir, std::string LegendTitle, unsigned int CutIndex)
{
   TH1D* HCuts_Pt  = (TH1D*) GetObjectFromPath(dir, "HCuts_Pt");
   TH1D* HCuts_Is  = (TH1D*) GetObjectFromPath(dir, "HCuts_I");
   TH1D* HCuts_TOF = (TH1D*) GetObjectFromPath(dir, "HCuts_TOF");
   char PtCutStr [1024]; sprintf (PtCutStr, "%.0f GeV", HCuts_Pt ->GetBinContent(CutIndex+1));
   char ICutStr  [1024]; sprintf (ICutStr,  "%.2f",     HCuts_Is ->GetBinContent(CutIndex+1));
   char TOFCutStr[1024]; sprintf (TOFCutStr,"%.3f",     HCuts_TOF->GetBinContent(CutIndex+1));

   TObject** Histos = new TObject*[10];
   std::vector<std::string> legend;

   TCanvas* c1 = new TCanvas("c1","c1,",600,600);          legend.clear();
   Histos[0] = (TH1*)tuple->BS_EtaIs;             legend.push_back("Before Cut");
   std::string        dEdxS_Legend    = "I_{as}";
   DrawSuperposedHistos((TH1**)Histos, legend, "COLZ",  "#eta", dEdxS_Legend.c_str(), 0,0, 0,0, false);
   c1->SetLogz(true); c1->SetRightMargin(0.15);
   //DrawPreliminary(LegendTitle, SQRTS, IntegratedLuminosityFromE(SQRTS));
   c1->SaveAs("EtaIs_BS.png"); //SaveCanvas(c1,SavePath,"EtaIs_BS", true);
   delete c1;
}
