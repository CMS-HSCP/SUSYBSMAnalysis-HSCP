// Usage:
// root -l -q -b step2_backgroundPrediction.C

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include "TRandom3.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphErrors.h>

#include "SUSYBSMAnalysis/Analyzer/interface/Regions.h"

void step2_backgroundPrediction(){
    ifstream infile;
    infile.open("configFile_readHist.txt");
    std::string line;
    std::string filename;
    std::string st_sample;
    int nPE, cutIndex;
    int rebineta,rebinih,rebinp,rebinmass;
    bool rebin;
    while(std::getline(infile,line)){
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> st_sample >> nPE >> cutIndex >> rebin >> rebineta >> rebinih >> rebinp >> rebinmass;
    }

    std::string outfilename_;
    if(!rebin)outfilename_ = filename+"_cutIndex"+to_string(cutIndex)+"_analysed";
    else outfilename_ = filename+"_cutIndex"+to_string(cutIndex)+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_analysed";

    std::cout << outfilename_ << std::endl;

    bool bool_rebin=rebin;
    
    TFile* ifile = new TFile((filename+".root").c_str());

    // histograms used for the mass prediction
    //------------


    //std::string dir = "analyzer/BaseName/";
    std::string dir = "HSCParticleAnalyzer/BaseName/";
    TH2F* eta_cutIndex_regA = (TH2F*)ifile->Get((dir+"Pred_EtaB").c_str())->Clone(); 
    TH2F* eta_cutIndex_regB = (TH2F*)ifile->Get((dir+"Pred_EtaS").c_str())->Clone(); 
    TH3F* ih_eta_cutIndex_regB = (TH3F*)ifile->Get((dir+"Pred_EtaI").c_str())->Clone(); 
    TH3F* eta_p_cutIndex_regC = (TH3F*)ifile->Get((dir+"Pred_EtaP").c_str())->Clone(); 
    TH1F* H_A = (TH1F*)ifile->Get((dir+"H_A").c_str())->Clone();
    TH1F* H_B = (TH1F*)ifile->Get((dir+"H_B").c_str())->Clone();
    TH1F* H_C = (TH1F*)ifile->Get((dir+"H_C").c_str())->Clone();
    TH2F* mass_cutIndex = (TH2F*)ifile->Get((dir+"Mass").c_str())->Clone();

    //------------

    Region ra_ias50;
    Region rc_ias50;

    Region rb_50ias60;
    Region rb_60ias70;
    Region rb_70ias80;
    Region rb_80ias90;
    Region rb_50ias90;
    Region rb_90ias100;

    Region rd_50ias60;
    Region rd_60ias70;
    Region rd_70ias80;
    Region rd_80ias90;
    Region rd_50ias90;
    Region rd_90ias100;
    
    Region rbc_50ias60;
    Region rbc_60ias70;
    Region rbc_70ias80;
    Region rbc_80ias90;
    Region rbc_50ias90;
    Region rbc_90ias100;

    Region ra_sc1;
    Region rb_sc1;
    Region rc_sc1;
    Region rd_sc1;
    Region rbc_sc1;
    
   
    // loading histograms used to validate the background estimate method in data --> base on Ias slices 
    // ------------------------------------------------------------------------------------------------------
    
    /*loadHistograms(ra_ias50,ifile,"regionA_med",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_ias50,ifile,"regionC_med",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_50ias60,ifile,"regionB_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60ias70,ifile,"regionB_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70ias80,ifile,"regionB_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80ias90,ifile,"regionB_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50ias90,ifile,"regionB_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_90ias100,ifile,"regionB_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rd_50ias60,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60ias70,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70ias80,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80ias90,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50ias90,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90ias100,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
     
    loadHistograms(rbc_50ias60,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60ias70,ifile,"regionD_60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70ias80,ifile,"regionD_70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80ias90,ifile,"regionD_80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50ias90,ifile,"regionD_50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90ias100,ifile,"regionD_90",bool_rebin,rebineta,rebinp,rebinih,rebinmass);

    loadHistograms(ra_sc1,ifile,"regA_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rb_sc1,ifile,"regB_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rc_sc1,ifile,"regC_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rd_sc1,ifile,"regD_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
    loadHistograms(rbc_sc1,ifile,"regD_sc1",bool_rebin,rebineta,rebinp,rebinih,rebinmass);*/


    loadHistograms(ra_ias50,ifile,"regionA_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_ias50,ifile,"regionC_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_50ias60,ifile,"regionB_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60ias70,ifile,"regionB_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70ias80,ifile,"regionB_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80ias90,ifile,"regionB_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50ias90,ifile,"regionB_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_90ias100,ifile,"regionB_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rd_50ias60,ifile,"regionD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60ias70,ifile,"regionD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70ias80,ifile,"regionD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80ias90,ifile,"regionD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90ias100,ifile,"regionD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
     
    loadHistograms(rbc_50ias60,ifile,"regionD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60ias70,ifile,"regionD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70ias80,ifile,"regionD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80ias90,ifile,"regionD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50ias90,ifile,"regionD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90ias100,ifile,"regionD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);



/*    loadHistograms(ra_ias50,ifile,"regA_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rc_ias50,ifile,"regC_ias50",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    
    loadHistograms(rb_50ias60,ifile,"regB_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_60ias70,ifile,"regB_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_70ias80,ifile,"regB_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_80ias90,ifile,"regB_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_50ias90,ifile,"regB_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rb_90ias100,ifile,"regB_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 

    loadHistograms(rd_50ias60,ifile,"regD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_60ias70,ifile,"regD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_70ias80,ifile,"regD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_80ias90,ifile,"regD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_50ias90,ifile,"regD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rd_90ias100,ifile,"regD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);
     
    loadHistograms(rbc_50ias60,ifile,"regD_50ias60",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_60ias70,ifile,"regD_60ias70",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_70ias80,ifile,"regD_70ias80",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_80ias90,ifile,"regD_80ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_50ias90,ifile,"regD_50ias90",bool_rebin,rebineta,rebinp,rebinih,rebinmass); 
    loadHistograms(rbc_90ias100,ifile,"regD_90ias100",bool_rebin,rebineta,rebinp,rebinih,rebinmass);*/
    
    // ------------------------------------------------------------------------------------------------------
    
    std::cout << "Regions loaded" << std::endl;

    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::cout << "saving... " << std::endl;

    // estimate the background in different Ias slices, each containing 10% of the statistic 
    // ------------------------------------------------------------------------------------------------------
    
    bckgEstimate(st_sample, rb_50ias60, rc_ias50, rbc_50ias60, ra_ias50, rd_50ias60, "50ias60", nPE);
    bckgEstimate(st_sample, rb_60ias70, rc_ias50, rbc_60ias70, ra_ias50, rd_60ias70, "60ias70", nPE);
    bckgEstimate(st_sample, rb_70ias80, rc_ias50, rbc_70ias80, ra_ias50, rd_70ias80, "70ias80", nPE);
    bckgEstimate(st_sample, rb_80ias90, rc_ias50, rbc_80ias90, ra_ias50, rd_80ias90, "80ias90", nPE);
    //bckgEstimate(st_sample, rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", nPE);
    bckgEstimate(st_sample, rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", nPE);
    //bckgEstimate(st_sample, rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100", nPE);

    
    // ------------------------------------------------------------------------------------------------------
   
    // bkg estimate for a selected cut index 
    // cutIndex = 3 --> pT > 60 GeV & Ias > 0.05
   


    //bckgEstimate_fromHistos(st_sample, *mass_cutIndex, *eta_cutIndex_regA, *eta_cutIndex_regB, *ih_eta_cutIndex_regB, *eta_p_cutIndex_regC, *H_A, *H_B, *H_C, cutIndex, nPE);
    //bckgEstimate_fromHistos(st_sample, *mass_cutIndex, *eta_cutIndex_regA, *eta_cutIndex_regB, *ih_eta_cutIndex_regB, *eta_p_cutIndex_regC, *H_A, *H_B, *H_C, cutIndex, 20);
    //bckgEstimate_fromHistos(st_sample, *mass_cutIndex, *eta_cutIndex_regA, *eta_cutIndex_regB, *ih_eta_cutIndex_regB, *eta_p_cutIndex_regC, *H_A, *H_B, *H_C, cutIndex, 50);
    //bckgEstimate_fromHistos(st_sample, *mass_cutIndex, *eta_cutIndex_regA, *eta_cutIndex_regB, *ih_eta_cutIndex_regB, *eta_p_cutIndex_regC, *H_A, *H_B, *H_C, cutIndex, 100);
    //bckgEstimate_fromHistos(st_sample, *mass_cutIndex, *eta_cutIndex_regA, *eta_cutIndex_regB, *ih_eta_cutIndex_regB, *eta_p_cutIndex_regC, *H_A, *H_B, *H_C, cutIndex, 200);
    //bckgEstimate_fromHistos(st_sample, *mass_cutIndex, *eta_cutIndex_regA, *eta_cutIndex_regB, *ih_eta_cutIndex_regB, *eta_p_cutIndex_regC, *H_A, *H_B, *H_C, cutIndex=28, nPE);

    delete ofile;
    delete mass_cutIndex;

    return;
}
