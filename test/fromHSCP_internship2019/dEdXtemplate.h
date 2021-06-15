//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr 10 16:59:33 2019 by ROOT version 6.10/09
// from TTree ttree/ttree
// found on file: test_minbias_aod.root
//////////////////////////////////////////////////////////

#ifndef dEdXtemplate_h
#define dEdXtemplate_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

class dEdXtemplate {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

// Fixed size dimensions of array or collections stored in the TTree if any.

   // Declaration of leaf types
   Int_t           runNumber;
   Int_t           event;
   Int_t           ngenpart;
   Int_t           gen_pdg[10000];   //[ngenpart]
   Float_t         gen_pt[10000];   //[ngenpart]
   Float_t         gen_eta[10000];   //[ngenpart]
   Float_t         gen_phi[10000];   //[ngenpart]
   Bool_t          gen_isHardProcess[10000];   //[ngenpart]
   Int_t           gen_status[10000];   //[ngenpart]
   Int_t           gen_moth_pdg[10000];   //[ngenpart]
   Int_t           gen_ndaughter[10000];   //[ngenpart]
   Int_t           gen_daughter_pdg[10000];   //[ngenpart]
   Int_t           ntracks;
   Float_t         track_pt[10000];   //[ntracks]
   Float_t         track_p[10000];   //[ntracks]
   Float_t         track_eta[10000];   //[ntracks]
   Float_t         track_phi[10000];   //[ntracks]
   Float_t         track_chi2[10000];   //[ntracks]
   Int_t           track_nvalidhits[10000];   //[ntracks]
   Int_t           track_index_hit[10000];   //[ntracks]
   Int_t           track_nhits[10000];   //[ntracks]
   Float_t         track_ih_ampl[10000];   //[ntracks]
   Float_t         track_ih_ampl_corr[10000];   //[ntracks]
   Int_t           ndedxhits;
   UInt_t          dedx_detid[100000];   //[ndedxhits]
   Int_t           dedx_subdetid[100000];   //[ndedxhits]
   Int_t           dedx_modulgeom[100000];   //[ndedxhits]
   Float_t         dedx_charge[100000];   //[ndedxhits]
   Float_t         dedx_pathlength[100000];   //[ndedxhits]
   Float_t         dedx_posx[100000];   //[ndedxhits]
   Float_t         dedx_posy[100000];   //[ndedxhits]
   Float_t         dedx_posz[100000];   //[ndedxhits]
   Bool_t          dedx_isstrip[100000];   //[ndedxhits]
   Bool_t          dedx_ispixel[100000];   //[ndedxhits]
   Bool_t          dedx_insideTkMod[100000];   //[ndedxhits]
   Int_t           sclus_firstsclus[100000];   //[ndedxhits]
   Float_t         sclus_barycenter[100000];   //[ndedxhits]
   Float_t         sclus_charge[100000];   //[ndedxhits]
   Float_t         sclus_errorclus[100000];   //[ndedxhits]
   Bool_t          sclus_ismerged[100000];   //[ndedxhits]
   Int_t           sclus_index_strip[100000];   //[ndedxhits]
   Int_t           sclus_nstrip[100000];   //[ndedxhits]
   Bool_t          sclus_sat254[100000];   //[ndedxhits]
   Bool_t          sclus_sat255[100000];   //[ndedxhits]
   Bool_t          sclus_shape[100000];   //[ndedxhits]
   Int_t           sclus_index_strip_corr[100000];   //[ndedxhits]
   Int_t           sclus_nstrip_corr[100000];   //[ndedxhits]
   Float_t         sclus_charge_corr[100000];   //[ndedxhits]
   Bool_t          sclus_clusclean[100000];   //[ndedxhits]
   Bool_t          sclus_clusclean2[100000];   //[ndedxhits]
   Int_t           sclus_index_simhit[100000];   //[ndedxhits]
   Int_t           sclus_nsimhit[100000];   //[ndedxhits]
   Float_t         sclus_eloss[100000];   //[ndedxhits]
   Int_t           nstrips;
   Int_t           strip_ampl[50000];   //[nstrips]
   Int_t           nstrips_corr;
   Int_t           strip_ampl_corr[50000];   //[nstrips_corr]
   Int_t           nsimhits;
   Int_t           simhit_pid[10000];   //[nsimhits]
   Int_t           simhit_process[10000];   //[nsimhits]
   Float_t         simhit_p[10000];   //[nsimhits]
   Float_t         simhit_eloss[10000];   //[nsimhits]
   Float_t         simhit_tof[10000];   //[nsimhits]
   Float_t         simhit_segment[10000];   //[nsimhits]
   Float_t         simhit_xentry[10000];   //[nsimhits]
   Float_t         simhit_yentry[10000];   //[nsimhits]
   Float_t         simhit_zentry[10000];   //[nsimhits]
   Float_t         simhit_xexit[10000];   //[nsimhits]
   Float_t         simhit_yexit[10000];   //[nsimhits]
   Float_t         simhit_zexit[10000];   //[nsimhits]

   // List of branches
   TBranch        *b_runNumber;   //!
   TBranch        *b_event;   //!
   TBranch        *b_ngenpart;   //!
   TBranch        *b_gen_pdg;   //!
   TBranch        *b_gen_pt;   //!
   TBranch        *b_gen_eta;   //!
   TBranch        *b_gen_phi;   //!
   TBranch        *b_gen_isHardProcess;   //!
   TBranch        *b_gen_status;   //!
   TBranch        *b_gen_moth_pdg;   //!
   TBranch        *b_gen_ndaughter;   //!
   TBranch        *b_gen_daughter_pdg;   //!
   TBranch        *b_ntracks;   //!
   TBranch        *b_track_pt;   //!
   TBranch        *b_track_p;   //!
   TBranch        *b_track_eta;   //!
   TBranch        *b_track_phi;   //!
   TBranch        *b_track_chi2;   //!
   TBranch        *b_track_nvalidhits;   //!
   TBranch        *b_track_index_hit;   //!
   TBranch        *b_track_nhits;   //!
   TBranch        *b_track_ih_ampl;   //!
   TBranch        *b_track_ih_ampl_corr;   //!
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

   dEdXtemplate(TTree *tree=0);
   virtual ~dEdXtemplate();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(TString name, bool uncorr=false);
   double getdEdX(std::vector <float> charge, std::vector <float> pathlength, std::vector <int> subdetId, std::vector <int> moduleGeometry, std::vector <bool> bool_cleaning, std::vector <bool> mustBeInside, double* scaleFactors, TH3* templateHisto);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef dEdXtemplate_cxx
dEdXtemplate::dEdXtemplate(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("test_minbias_aod.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("test_minbias_aod.root");
      }
      TDirectory * dir = (TDirectory*)f->Get("test_minbias_aod.root:/stage");
      dir->GetObject("ttree",tree);

   }
   Init(tree);
}

dEdXtemplate::~dEdXtemplate()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t dEdXtemplate::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t dEdXtemplate::LoadTree(Long64_t entry)
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

void dEdXtemplate::Init(TTree *tree)
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
   fChain->SetBranchAddress("ngenpart", &ngenpart, &b_ngenpart);
   fChain->SetBranchAddress("gen_pdg", &gen_pdg, &b_gen_pdg);
   fChain->SetBranchAddress("gen_pt", &gen_pt, &b_gen_pt);
   fChain->SetBranchAddress("gen_eta", &gen_eta, &b_gen_eta);
   fChain->SetBranchAddress("gen_phi", &gen_phi, &b_gen_phi);
   fChain->SetBranchAddress("gen_isHardProcess", &gen_isHardProcess, &b_gen_isHardProcess);
   fChain->SetBranchAddress("gen_status", &gen_status, &b_gen_status);
   fChain->SetBranchAddress("gen_moth_pdg", &gen_moth_pdg, &b_gen_moth_pdg);
   fChain->SetBranchAddress("gen_ndaughter", &gen_ndaughter, &b_gen_ndaughter);
   fChain->SetBranchAddress("gen_daughter_pdg", &gen_daughter_pdg, &b_gen_daughter_pdg);
   fChain->SetBranchAddress("ntracks", &ntracks, &b_ntracks);
   fChain->SetBranchAddress("track_pt", track_pt, &b_track_pt);
   fChain->SetBranchAddress("track_p", track_p, &b_track_p);
   fChain->SetBranchAddress("track_eta", track_eta, &b_track_eta);
   fChain->SetBranchAddress("track_phi", track_phi, &b_track_phi);
   fChain->SetBranchAddress("track_chi2", track_chi2, &b_track_chi2);
   fChain->SetBranchAddress("track_nvalidhits", track_nvalidhits, &b_track_nvalidhits);
   fChain->SetBranchAddress("track_index_hit", track_index_hit, &b_track_index_hit);
   fChain->SetBranchAddress("track_nhits", track_nhits, &b_track_nhits);
   fChain->SetBranchAddress("track_ih_ampl", track_ih_ampl, &b_track_ih_ampl);
   fChain->SetBranchAddress("track_ih_ampl_corr", track_ih_ampl_corr, &b_track_ih_ampl_corr);
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
   Notify();
}

Bool_t dEdXtemplate::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void dEdXtemplate::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t dEdXtemplate::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef dEdXtemplate_cxx
