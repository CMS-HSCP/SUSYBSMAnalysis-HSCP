//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Mon Nov  4 12:08:02 2019 by ROOT version 6.14/09
// from TTree ttree/ttree
// found on file: /opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2017B/191023_141349/0000/nt_data_aod_44.root
//////////////////////////////////////////////////////////

#ifndef run2study_h
#define run2study_h
#define nmax_tr 1000
#define nmax_cl 10000
#define nmax_st 30000
#define nmax_mu 200
#define nmax_hscp 200
#define calmax 50000

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class run2study {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.
//
   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           event;
   bool            boolILumi;
   Float_t         InstLumi;   
   Int_t           npv;
   Int_t           ngoodpv;
   Bool_t          hlt_mu50;
   Bool_t          hlt_tkmu100;
   Bool_t          hlt_oldmu100;
   Bool_t          hlt_pfmet_mht;
   Int_t           ngenpart;
   Int_t           gen_pdg[1];   //[ngenpart]
   Float_t         gen_pt[1];   //[ngenpart]
   Float_t         gen_eta[1];   //[ngenpart]
   Float_t         gen_phi[1];   //[ngenpart]
   Float_t         gen_mass[1];   //[ngenpart]
   Bool_t          gen_isHardProcess[1];   //[ngenpart]
   Int_t           gen_status[1];   //[ngenpart]
   Int_t           gen_moth_pdg[1];   //[ngenpart]
   Int_t           gen_ndaughter[1];   //[ngenpart]
   Int_t           gen_daughter_pdg[1];   //[ngenpart]
   Int_t           ntracks;
   Float_t         track_pt[nmax_tr];   //[ntracks]
   Float_t         track_pterr[nmax_tr];   //[ntracks]
   Float_t         track_p[nmax_tr];   //[ntracks]
   Float_t         track_eta[nmax_tr];   //[ntracks]
   Float_t         track_phi[nmax_tr];   //[ntracks]
   Float_t         track_charge[nmax_tr];   //[ntracks]
   Float_t         track_chi2[nmax_tr];   //[ntracks]
   Int_t           track_nvalidhits[nmax_tr];   //[ntracks]
   Int_t           track_npixhits[nmax_tr];   //[ntracks]
   Int_t           track_missing[nmax_tr];   //[ntracks]
   Float_t         track_validfraction[nmax_tr];   //[ntracks]
   Float_t         track_validlast[nmax_tr];   //[ntracks]
   Int_t           track_qual[nmax_tr];   //[ntracks]
   Float_t         track_dz[nmax_tr];   //[ntracks]
   Float_t         track_dxy[nmax_tr];   //[ntracks]
   Int_t           track_index_hit[nmax_tr];   //[ntracks]
   Int_t           track_nhits[nmax_tr];   //[ntracks]
   Int_t           track_prescale[nmax_tr];   //[ntracks]
   Float_t         track_ih_ampl[nmax_tr];   //[ntracks]
   Float_t         track_ih_ampl_corr[nmax_tr];   //[ntracks]
   Float_t         track_ias_ampl[nmax_tr];   //[ntracks]
   Float_t         track_ias_ampl_corr[nmax_tr];   //[ntracks]
   Int_t           track_nclus[nmax_tr];   //[ntracks]
   Int_t           track_clus_PID[nmax_tr][10];   //[ntracks][10]
   Float_t         track_probQ[nmax_tr];   //[ntracks]
   Float_t         track_probXY[nmax_tr];   //[ntracks]
   Int_t           ndedxhits;
   UInt_t          dedx_detid[nmax_cl];   //[ndedxhits]
   Int_t           dedx_subdetid[nmax_cl];   //[ndedxhits]
   Int_t           dedx_modulgeom[nmax_cl];   //[ndedxhits]
   Float_t         dedx_charge[nmax_cl];   //[ndedxhits]
   Float_t         dedx_pathlength[nmax_cl];   //[ndedxhits]
   Float_t         dedx_posx[nmax_cl];   //[ndedxhits]
   Float_t         dedx_posy[nmax_cl];   //[ndedxhits]
   Float_t         dedx_posz[nmax_cl];   //[ndedxhits]
   Bool_t          dedx_isstrip[nmax_cl];   //[ndedxhits]
   Bool_t          dedx_ispixel[nmax_cl];   //[ndedxhits]
   Bool_t          dedx_insideTkMod[nmax_cl];   //[ndedxhits]
   Int_t           sclus_firstsclus[nmax_cl];   //[ndedxhits]
   Float_t         sclus_barycenter[nmax_cl];   //[ndedxhits]
   Float_t         sclus_charge[nmax_cl];   //[ndedxhits]
   Float_t         sclus_errorclus[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_ismerged[nmax_cl];   //[ndedxhits]
   Int_t           sclus_index_strip[nmax_cl];   //[ndedxhits]
   Int_t           sclus_nstrip[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_sat254[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_sat255[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_shape[nmax_cl];   //[ndedxhits]
   Int_t           sclus_index_strip_corr[nmax_cl];   //[ndedxhits]
   Int_t           sclus_nstrip_corr[nmax_cl];   //[ndedxhits]
   Float_t         sclus_charge_corr[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_clusclean[nmax_cl];   //[ndedxhits]
   Bool_t          sclus_clusclean2[nmax_cl];   //[ndedxhits]
   Int_t           sclus_index_simhit[nmax_cl];   //[ndedxhits]
   Int_t           sclus_nsimhit[nmax_cl];   //[ndedxhits]
   Float_t         sclus_eloss[nmax_cl];   //[ndedxhits]
   Int_t           nstrips;
   Int_t           strip_ampl[nmax_st];   //[nstrips]
   Int_t           nstrips_corr;
   Int_t           strip_ampl_corr[nmax_st];   //[nstrips_corr]
   Int_t           nsimhits;
   Int_t           simhit_pid[1];   //[nsimhits]
   Int_t           simhit_process[1];   //[nsimhits]
   Float_t         simhit_p[1];   //[nsimhits]
   Float_t         simhit_eloss[1];   //[nsimhits]
   Float_t         simhit_tof[1];   //[nsimhits]
   Float_t         simhit_segment[1];   //[nsimhits]
   Float_t         simhit_xentry[1];   //[nsimhits]
   Float_t         simhit_yentry[1];   //[nsimhits]
   Float_t         simhit_zentry[1];   //[nsimhits]
   Float_t         simhit_xexit[1];   //[nsimhits]
   Float_t         simhit_yexit[1];   //[nsimhits]
   Float_t         simhit_zexit[1];   //[nsimhits]
   Int_t           nmuons;
   Float_t         muon_pt[nmax_mu];   //[nmuons]
   Float_t         muon_ptSA[nmax_mu];   //[nmuons]
   Float_t         muon_ptIT[nmax_mu];   //[nmuons]
   Float_t         muon_p[nmax_mu];   //[nmuons]
   Float_t         muon_eta[nmax_mu];   //[nmuons]
   Float_t         muon_phi[nmax_mu];   //[nmuons]
   Float_t         muon_comb_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_comb_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_comb_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_comb_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_dt_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_dt_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_dt_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_dt_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_csc_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_csc_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_csc_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_csc_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_newcomb_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_newcomb_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_newcomb_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_newcomb_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_newdt_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_newdt_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_newdt_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_newdt_vertextime[nmax_mu];   //[nmuons]
   Float_t         muon_newcsc_inversebeta[nmax_mu];   //[nmuons]
   Float_t         muon_newcsc_inversebetaerr[nmax_mu];   //[nmuons]
   Int_t           muon_newcsc_tofndof[nmax_mu];   //[nmuons]
   Float_t         muon_newcsc_vertextime[nmax_mu];   //[nmuons]
   Int_t           nhscp;
   Int_t           hscp_gen_id[nmax_hscp];   //[nhscp]
   Float_t         hscp_gen_dr[nmax_hscp];   //[nhscp]
   Int_t           hscp_track_idx[nmax_hscp];   //[nhscp]
   Int_t           hscp_muon_idx[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso0_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso0_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso0_hcal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso1_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso1_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso1_hcal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso2_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso2_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso2_hcal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso3_tk[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso3_ecal[nmax_hscp];   //[nhscp]
   Float_t         hscp_iso3_hcal[nmax_hscp];   //[nhscp]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_event;   //!
   TBranch        *b_InstLumi;   //!
   TBranch        *b_npv;   //!
   TBranch        *b_ngoodpv;   //!
   TBranch        *b_hlt_mu50;   //!
   TBranch        *b_hlt_tkmu100;   //!
   TBranch        *b_hlt_oldmu100;   //!
   TBranch        *b_hlt_pfmet_mht;   //!
   TBranch        *b_ngenpart;   //!
   TBranch        *b_gen_pdg;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_mass;   //!
   TBranch        *b_gen_isHardProcess;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_moth_pdg;   //!
   TBranch        *b_gen_ndaughter;   //!
   TBranch        *b_gen_daughter_pdg;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_pterr;   //!
   TBranch        *b_track_p;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_charge;   //!
   TBranch        *b_track_chi2;   //!
   TBranch        *b_track_nvalidhits;   //!
   TBranch        *b_track_npixhits;   //!
   TBranch        *b_track_missing;   //!
   TBranch        *b_track_validfraction;   //!
   TBranch        *b_track_validlast;   //!
   TBranch        *b_track_qual;   //!
   TBranch        *b_track_dz;   //!
   TBranch        *b_track_dxy;   //!
   TBranch        *b_track_index_hit;   //!
   TBranch        *b_track_nhits;   //!
   TBranch        *b_track_prescale;   //!
   TBranch        *b_track_ih_ampl;   //!
   TBranch        *b_track_ih_ampl_corr;   //!
   TBranch        *b_track_ias_ampl;   //!
   TBranch        *b_track_ias_ampl_corr;   //!
   TBranch        *b_track_nclus; //!
   TBranch        *b_track_clus_PID; //!
   TBranch        *b_track_probQ;   //!
   TBranch        *b_track_probXY;   //!
   TBranch        *b_ndedxhits;   //!
   TBranch        *b_dedx_detid;   //!
   TBranch        *b_dedx_subdetid;   //!
   TBranch        *b_dedx_modulgeom;   //!
   TBranch        *b_dedx_charge;   //!
   TBranch        *b_dedx_pathlength;   //!
   TBranch        *b_dedx_posx;   //!
   TBranch        *b_dedx_posy;   //!
   TBranch        *b_dedx_posz;   //!
   TBranch        *b_dedx_isstrip;   //!
   TBranch        *b_dedx_ispixel;   //!
   TBranch        *b_dedx_insideTkMod;   //!
   TBranch        *b_sclus_firstsclus;   //!
   TBranch        *b_sclus_barycenter;   //!
   TBranch        *b_sclus_charge;   //!
   TBranch        *b_sclus_errorclus;   //!
   TBranch        *b_sclus_ismerged;   //!
   TBranch        *b_sclus_index_strip;   //!
   TBranch        *b_sclus_nstrip;   //!
   TBranch        *b_sclus_sat254;   //!
   TBranch        *b_sclus_sat255;   //!
   TBranch        *b_sclus_shape;   //!
   TBranch        *b_sclus_index_strip_corr;   //!
   TBranch        *b_sclus_nstrip_corr;   //!
   TBranch        *b_sclus_charge_corr;   //!
   TBranch        *b_sclus_clusclean;   //!
   TBranch        *b_sclus_clusclean2;   //!
   TBranch        *b_sclus_index_simhit;   //!
   TBranch        *b_sclus_nsimhit;   //!
   TBranch        *b_sclus_eloss;   //!
   TBranch        *b_nstrips;   //!
   TBranch        *b_strip_ampl;   //!
   TBranch        *b_nstrips_corr;   //!
   TBranch        *b_strip_ampl_corr;   //!
   TBranch        *b_nsimhits;   //!
   TBranch        *b_simhit_pid;   //!
   TBranch        *b_simhit_process;   //!
   TBranch        *b_simhit_p;   //!
   TBranch        *b_simhit_eloss;   //!
   TBranch        *b_simhit_tof;   //!
   TBranch        *b_simhit_segment;   //!
   TBranch        *b_simhit_xentry;   //!
   TBranch        *b_simhit_yentry;   //!
   TBranch        *b_simhit_zentry;   //!
   TBranch        *b_simhit_xexit;   //!
   TBranch        *b_simhit_yexit;   //!
   TBranch        *b_simhit_zexit;   //!
   TBranch        *b_nmuons;   //!
   TBranch        *b_muon_pt;   //!
   TBranch        *b_muon_ptSA;   //!
   TBranch        *b_muon_ptIT;   //!
   TBranch        *b_muon_p;   //!
   TBranch        *b_muon_eta;   //!
   TBranch        *b_muon_phi;   //!
   TBranch        *b_muon_comb_inversebeta;   //!
   TBranch        *b_muon_comb_inversebetaerr;   //!
   TBranch        *b_muon_comb_tofndof;   //!
   TBranch        *b_muon_comb_vertextime;   //!
   TBranch        *b_muon_dt_inversebeta;   //!
   TBranch        *b_muon_dt_inversebetaerr;   //!
   TBranch        *b_muon_dt_tofndof;   //!
   TBranch        *b_muon_dt_vertextime;   //!
   TBranch        *b_muon_csc_inversebeta;   //!
   TBranch        *b_muon_csc_inversebetaerr;   //!
   TBranch        *b_muon_csc_tofndof;   //!
   TBranch        *b_muon_csc_vertextime;   //!
   TBranch        *b_muon_newcomb_inversebeta;   //!
   TBranch        *b_muon_newcomb_inversebetaerr;   //!
   TBranch        *b_muon_newcomb_tofndof;   //!
   TBranch        *b_muon_newcomb_vertextime;   //!
   TBranch        *b_muon_newdt_inversebeta;   //!
   TBranch        *b_muon_newdt_inversebetaerr;   //!
   TBranch        *b_muon_newdt_tofndof;   //!
   TBranch        *b_muon_newdt_vertextime;   //!
   TBranch        *b_muon_newcsc_inversebeta;   //!
   TBranch        *b_muon_newcsc_inversebetaerr;   //!
   TBranch        *b_muon_newcsc_tofndof;   //!
   TBranch        *b_muon_newcsc_vertextime;   //!
   TBranch        *b_nhscp;   //!
   TBranch        *b_hscp_gen_id;   //!
   TBranch        *b_hscp_gen_dr;   //!
   TBranch        *b_hscp_track_idx;   //!
   TBranch        *b_hscp_muon_idx;   //!
   TBranch        *b_hscp_iso0_tk;   //!
   TBranch        *b_hscp_iso0_ecal;   //!
   TBranch        *b_hscp_iso0_hcal;   //!
   TBranch        *b_hscp_iso1_tk;   //!
   TBranch        *b_hscp_iso1_ecal;   //!
   TBranch        *b_hscp_iso1_hcal;   //!
   TBranch        *b_hscp_iso2_tk;   //!
   TBranch        *b_hscp_iso2_ecal;   //!
   TBranch        *b_hscp_iso2_hcal;   //!
   TBranch        *b_hscp_iso3_tk;   //!
   TBranch        *b_hscp_iso3_ecal;   //!
   TBranch        *b_hscp_iso3_hcal;   //!

   run2study(TTree *tree=0);
   virtual ~run2study();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(int year, TString Letter);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
   double getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto, int n_estim, double dropLowerDeDxValue, int &nv, int &ns);
//   float FMIP(const vector<float>& charge, const<vector>& path,  float thre);
   float FMIP(std::vector <float> charge, std::vector <float> path,  float thre);
   int GetLayerLabel(int subdetid_, UInt_t detid_, int year);
   void loadSFPixelCalib();
   float GetSFPixel(int subdetid_, UInt_t detid_, int year, float eta, int run);

   int icalib;
   int icalib2017;
   int icalib2018;
   int pixVal[calmax];
   int layerSideVal[calmax];
   int ladderBladeVal[calmax];
   float etaMinVal[calmax];
   float etaMaxVal[calmax];
   int irunMinVal[calmax];
   int irunMaxVal[calmax];
   float scaleVal[calmax]; 

};

#endif

#ifdef run2study_cxx
run2study::run2study(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2017B/191023_141349/0000/nt_data_aod_44.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2017B/191023_141349/0000/nt_data_aod_44.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("/opt/sbg/cms/ui3_data1/ccollard/HSCP_prod/prodOct2019_CMSSW_10_6_2/SingleMuon/run2017B/191023_141349/0000/nt_data_aod_44.root:/stage");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

run2study::~run2study()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t run2study::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t run2study::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void run2study::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("runNumber", &runNumber, &b_runNumber);
   fChain->SetBranchAddress("event", &event, &b_event);
   if (fChain->GetBranch("InstLumi")) {
    boolILumi = true;
    fChain->SetBranchAddress("InstLumi"   , &InstLumi, &b_InstLumi);
   }
   else {
    boolILumi=false;
   }
   fChain->SetBranchAddress("npv", &npv, &b_npv);
   fChain->SetBranchAddress("ngoodpv", &ngoodpv, &b_ngoodpv);
   fChain->SetBranchAddress("hlt_mu50", &hlt_mu50, &b_hlt_mu50);
   fChain->SetBranchAddress("hlt_tkmu100", &hlt_tkmu100, &b_hlt_tkmu100);
   fChain->SetBranchAddress("hlt_oldmu100", &hlt_oldmu100, &b_hlt_oldmu100);
   fChain->SetBranchAddress("hlt_pfmet_mht", &hlt_pfmet_mht, &b_hlt_pfmet_mht);
   fChain->SetBranchAddress("ngenpart", &ngenpart, &b_ngenpart);
   fChain->SetBranchAddress("gen_pdg", &gen_pdg, &b_gen_pdg);
   fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_eta", &gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_mass", &gen_mass, &b_gen_mass);
   fChain->SetBranchAddress("gen_isHardProcess", &gen_isHardProcess, &b_gen_isHardProcess);
   fChain->SetBranchAddress("gen_status", &gen_status, &b_gen_status);
   fChain->SetBranchAddress("gen_moth_pdg", &gen_moth_pdg, &b_gen_moth_pdg);
   fChain->SetBranchAddress("gen_ndaughter", &gen_ndaughter, &b_gen_ndaughter);
   fChain->SetBranchAddress("gen_daughter_pdg", &gen_daughter_pdg, &b_gen_daughter_pdg);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("track_pt", track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_pterr", track_pterr, &b_track_pterr);
   fChain->SetBranchAddress("track_p", track_p, &b_track_p);
   fChain->SetBranchAddress("track_eta", track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_charge", track_charge, &b_track_charge);
   fChain->SetBranchAddress("track_chi2", track_chi2, &b_track_chi2);
   fChain->SetBranchAddress("track_nvalidhits", track_nvalidhits, &b_track_nvalidhits);
   fChain->SetBranchAddress("track_npixhits", track_npixhits, &b_track_npixhits);
   fChain->SetBranchAddress("track_missing", track_missing, &b_track_missing);
   fChain->SetBranchAddress("track_validfraction", track_validfraction, &b_track_validfraction);
   fChain->SetBranchAddress("track_validlast", track_validlast, &b_track_validlast);
   fChain->SetBranchAddress("track_qual", track_qual, &b_track_qual);
   fChain->SetBranchAddress("track_dz", track_dz, &b_track_dz);
   fChain->SetBranchAddress("track_dxy", track_dxy, &b_track_dxy);
   fChain->SetBranchAddress("track_index_hit", track_index_hit, &b_track_index_hit);
   fChain->SetBranchAddress("track_nhits", track_nhits, &b_track_nhits);
   fChain->SetBranchAddress("track_prescale", track_prescale, &b_track_prescale);
   fChain->SetBranchAddress("track_ih_ampl", track_ih_ampl, &b_track_ih_ampl);
   fChain->SetBranchAddress("track_ih_ampl_corr", track_ih_ampl_corr, &b_track_ih_ampl_corr);
   fChain->SetBranchAddress("track_ias_ampl", track_ias_ampl, &b_track_ias_ampl);
   fChain->SetBranchAddress("track_ias_ampl_corr", track_ias_ampl_corr, &b_track_ias_ampl_corr);
   fChain->SetBranchAddress("track_nclus", track_nclus, &b_track_nclus);
   fChain->SetBranchAddress("track_clus_PID", track_clus_PID, &b_track_clus_PID);
   fChain->SetBranchAddress("track_probQ", track_probQ, &b_track_probQ);
   fChain->SetBranchAddress("track_probXY", track_probXY, &b_track_probXY);
   fChain->SetBranchAddress("ndedxhits", &ndedxhits, &b_ndedxhits);
   fChain->SetBranchAddress("dedx_detid", dedx_detid, &b_dedx_detid);
   fChain->SetBranchAddress("dedx_subdetid", dedx_subdetid, &b_dedx_subdetid);
   fChain->SetBranchAddress("dedx_modulgeom", dedx_modulgeom, &b_dedx_modulgeom);
   fChain->SetBranchAddress("dedx_charge", dedx_charge, &b_dedx_charge);
   fChain->SetBranchAddress("dedx_pathlength", dedx_pathlength, &b_dedx_pathlength);
   fChain->SetBranchAddress("dedx_posx", dedx_posx, &b_dedx_posx);
   fChain->SetBranchAddress("dedx_posy", dedx_posy, &b_dedx_posy);
   fChain->SetBranchAddress("dedx_posz", dedx_posz, &b_dedx_posz);
   fChain->SetBranchAddress("dedx_isstrip", dedx_isstrip, &b_dedx_isstrip);
   fChain->SetBranchAddress("dedx_ispixel", dedx_ispixel, &b_dedx_ispixel);
   fChain->SetBranchAddress("dedx_insideTkMod", dedx_insideTkMod, &b_dedx_insideTkMod);
   fChain->SetBranchAddress("sclus_firstsclus", sclus_firstsclus, &b_sclus_firstsclus);
   fChain->SetBranchAddress("sclus_barycenter", sclus_barycenter, &b_sclus_barycenter);
   fChain->SetBranchAddress("sclus_charge", sclus_charge, &b_sclus_charge);
   fChain->SetBranchAddress("sclus_errorclus", sclus_errorclus, &b_sclus_errorclus);
   fChain->SetBranchAddress("sclus_ismerged", sclus_ismerged, &b_sclus_ismerged);
   fChain->SetBranchAddress("sclus_index_strip", sclus_index_strip, &b_sclus_index_strip);
   fChain->SetBranchAddress("sclus_nstrip", sclus_nstrip, &b_sclus_nstrip);
   fChain->SetBranchAddress("sclus_sat254", sclus_sat254, &b_sclus_sat254);
   fChain->SetBranchAddress("sclus_sat255", sclus_sat255, &b_sclus_sat255);
   fChain->SetBranchAddress("sclus_shape", sclus_shape, &b_sclus_shape);
   fChain->SetBranchAddress("sclus_index_strip_corr", sclus_index_strip_corr, &b_sclus_index_strip_corr);
   fChain->SetBranchAddress("sclus_nstrip_corr", sclus_nstrip_corr, &b_sclus_nstrip_corr);
   fChain->SetBranchAddress("sclus_charge_corr", sclus_charge_corr, &b_sclus_charge_corr);
   fChain->SetBranchAddress("sclus_clusclean", sclus_clusclean, &b_sclus_clusclean);
   fChain->SetBranchAddress("sclus_clusclean2", sclus_clusclean2, &b_sclus_clusclean2);
   fChain->SetBranchAddress("sclus_index_simhit", sclus_index_simhit, &b_sclus_index_simhit);
   fChain->SetBranchAddress("sclus_nsimhit", sclus_nsimhit, &b_sclus_nsimhit);
   fChain->SetBranchAddress("sclus_eloss", sclus_eloss, &b_sclus_eloss);
   fChain->SetBranchAddress("nstrips", &nstrips, &b_nstrips);
   fChain->SetBranchAddress("strip_ampl", strip_ampl, &b_strip_ampl);
   fChain->SetBranchAddress("nstrips_corr", &nstrips_corr, &b_nstrips_corr);
   fChain->SetBranchAddress("strip_ampl_corr", strip_ampl_corr, &b_strip_ampl_corr);
   fChain->SetBranchAddress("nsimhits", &nsimhits, &b_nsimhits);
   fChain->SetBranchAddress("simhit_pid", &simhit_pid, &b_simhit_pid);
   fChain->SetBranchAddress("simhit_process", &simhit_process, &b_simhit_process);
   fChain->SetBranchAddress("simhit_p", &simhit_p, &b_simhit_p);
   fChain->SetBranchAddress("simhit_eloss", &simhit_eloss, &b_simhit_eloss);
   fChain->SetBranchAddress("simhit_tof", &simhit_tof, &b_simhit_tof);
   fChain->SetBranchAddress("simhit_segment", &simhit_segment, &b_simhit_segment);
   fChain->SetBranchAddress("simhit_xentry", &simhit_xentry, &b_simhit_xentry);
   fChain->SetBranchAddress("simhit_yentry", &simhit_yentry, &b_simhit_yentry);
   fChain->SetBranchAddress("simhit_zentry", &simhit_zentry, &b_simhit_zentry);
   fChain->SetBranchAddress("simhit_xexit", &simhit_xexit, &b_simhit_xexit);
   fChain->SetBranchAddress("simhit_yexit", &simhit_yexit, &b_simhit_yexit);
   fChain->SetBranchAddress("simhit_zexit", &simhit_zexit, &b_simhit_zexit);
   fChain->SetBranchAddress("nmuons", &nmuons, &b_nmuons);
   fChain->SetBranchAddress("muon_pt", muon_pt, &b_muon_pt);
   fChain->SetBranchAddress("muon_ptSA", muon_ptSA, &b_muon_ptSA);
   fChain->SetBranchAddress("muon_ptIT", muon_ptIT, &b_muon_ptIT);
   fChain->SetBranchAddress("muon_p", muon_p, &b_muon_p);
   fChain->SetBranchAddress("muon_eta", muon_eta, &b_muon_eta);
   fChain->SetBranchAddress("muon_phi", muon_phi, &b_muon_phi);
   fChain->SetBranchAddress("muon_comb_inversebeta", muon_comb_inversebeta, &b_muon_comb_inversebeta);
   fChain->SetBranchAddress("muon_comb_inversebetaerr", muon_comb_inversebetaerr, &b_muon_comb_inversebetaerr);
   fChain->SetBranchAddress("muon_comb_tofndof", muon_comb_tofndof, &b_muon_comb_tofndof);
   fChain->SetBranchAddress("muon_comb_vertextime", muon_comb_vertextime, &b_muon_comb_vertextime);
   fChain->SetBranchAddress("muon_dt_inversebeta", muon_dt_inversebeta, &b_muon_dt_inversebeta);
   fChain->SetBranchAddress("muon_dt_inversebetaerr", muon_dt_inversebetaerr, &b_muon_dt_inversebetaerr);
   fChain->SetBranchAddress("muon_dt_tofndof", muon_dt_tofndof, &b_muon_dt_tofndof);
   fChain->SetBranchAddress("muon_dt_vertextime", muon_dt_vertextime, &b_muon_dt_vertextime);
   fChain->SetBranchAddress("muon_csc_inversebeta", muon_csc_inversebeta, &b_muon_csc_inversebeta);
   fChain->SetBranchAddress("muon_csc_inversebetaerr", muon_csc_inversebetaerr, &b_muon_csc_inversebetaerr);
   fChain->SetBranchAddress("muon_csc_tofndof", muon_csc_tofndof, &b_muon_csc_tofndof);
   fChain->SetBranchAddress("muon_csc_vertextime", muon_csc_vertextime, &b_muon_csc_vertextime);
   if (fChain->GetBranch("muon_newcomb_inversebeta")) {
   fChain->SetBranchAddress("muon_newcomb_inversebeta", muon_newcomb_inversebeta, &b_muon_newcomb_inversebeta);
   fChain->SetBranchAddress("muon_newcomb_inversebetaerr", muon_newcomb_inversebetaerr, &b_muon_newcomb_inversebetaerr);
   fChain->SetBranchAddress("muon_newcomb_tofndof", muon_newcomb_tofndof, &b_muon_newcomb_tofndof);
   fChain->SetBranchAddress("muon_newcomb_vertextime", muon_newcomb_vertextime, &b_muon_newcomb_vertextime);
   fChain->SetBranchAddress("muon_newdt_inversebeta", muon_newdt_inversebeta, &b_muon_newdt_inversebeta);
   fChain->SetBranchAddress("muon_newdt_inversebetaerr", muon_newdt_inversebetaerr, &b_muon_newdt_inversebetaerr);
   fChain->SetBranchAddress("muon_newdt_tofndof", muon_newdt_tofndof, &b_muon_newdt_tofndof);
   fChain->SetBranchAddress("muon_newdt_vertextime", muon_newdt_vertextime, &b_muon_newdt_vertextime);
   fChain->SetBranchAddress("muon_newcsc_inversebeta", muon_newcsc_inversebeta, &b_muon_newcsc_inversebeta);
   fChain->SetBranchAddress("muon_newcsc_inversebetaerr", muon_newcsc_inversebetaerr, &b_muon_newcsc_inversebetaerr);
   fChain->SetBranchAddress("muon_newcsc_tofndof", muon_newcsc_tofndof, &b_muon_newcsc_tofndof);
   fChain->SetBranchAddress("muon_newcsc_vertextime", muon_newcsc_vertextime, &b_muon_newcsc_vertextime);
   }
   fChain->SetBranchAddress("nhscp", &nhscp, &b_nhscp);
   fChain->SetBranchAddress("hscp_gen_id", hscp_gen_id, &b_hscp_gen_id);
   fChain->SetBranchAddress("hscp_gen_dr", hscp_gen_dr, &b_hscp_gen_dr);
   fChain->SetBranchAddress("hscp_track_idx", hscp_track_idx, &b_hscp_track_idx);
   fChain->SetBranchAddress("hscp_muon_idx", hscp_muon_idx, &b_hscp_muon_idx);
   fChain->SetBranchAddress("hscp_iso0_tk", hscp_iso0_tk, &b_hscp_iso0_tk);
   fChain->SetBranchAddress("hscp_iso0_ecal", hscp_iso0_ecal, &b_hscp_iso0_ecal);
   fChain->SetBranchAddress("hscp_iso0_hcal", hscp_iso0_hcal, &b_hscp_iso0_hcal);
   fChain->SetBranchAddress("hscp_iso1_tk", hscp_iso1_tk, &b_hscp_iso1_tk);
   fChain->SetBranchAddress("hscp_iso1_ecal", hscp_iso1_ecal, &b_hscp_iso1_ecal);
   fChain->SetBranchAddress("hscp_iso1_hcal", hscp_iso1_hcal, &b_hscp_iso1_hcal);
   fChain->SetBranchAddress("hscp_iso2_tk", hscp_iso2_tk, &b_hscp_iso2_tk);
   fChain->SetBranchAddress("hscp_iso2_ecal", hscp_iso2_ecal, &b_hscp_iso2_ecal);
   fChain->SetBranchAddress("hscp_iso2_hcal", hscp_iso2_hcal, &b_hscp_iso2_hcal);
   fChain->SetBranchAddress("hscp_iso3_tk", hscp_iso3_tk, &b_hscp_iso3_tk);
   fChain->SetBranchAddress("hscp_iso3_ecal", hscp_iso3_ecal, &b_hscp_iso3_ecal);
   fChain->SetBranchAddress("hscp_iso3_hcal", hscp_iso3_hcal, &b_hscp_iso3_hcal);
   Notify();
}

Bool_t run2study::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void run2study::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t run2study::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef run2study_cxx
