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
    int nPE, cutIndex;
    int rebineta,rebinih,rebinp,rebinmass;
    bool rebin;
    while(std::getline(infile,line)){
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> nPE >> cutIndex >> rebin >> rebineta >> rebinih >> rebinp >> rebinmass;
    }

    std::string outfilename_;
    if(!rebin)outfilename_ = filename+"_cutIndex"+to_string(cutIndex)+"_analysed";
    else outfilename_ = filename+"_cutIndex"+to_string(cutIndex)+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_analysed";

    std::cout << outfilename_ << std::endl;

    bool bool_rebin=rebin;
    
    TFile* ifile = new TFile((filename+".root").c_str());

    // histograms used for the mass prediction
    //------------

    std::string dir = "analyzer/BaseName/";
    TH2F* eta_cutIndex_regA = (TH2F*)ifile->Get((dir+"Pred_EtaB").c_str())->Clone(); 
    TH2F* eta_cutIndex_regB = (TH2F*)ifile->Get((dir+"Pred_EtaS").c_str())->Clone(); 
    TH3F* ih_eta_cutIndex_regB = (TH3F*)ifile->Get((dir+"Pred_EtaI").c_str())->Clone(); 
    TH3F* eta_p_cutIndex_regC = (TH3F*)ifile->Get((dir+"Pred_EtaP").c_str())->Clone(); 

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
    
   
    // loading histograms used to validate the background estimate method in data --> base on Ias slices 
    // ------------------------------------------------------------------------------------------------------
    
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
    
    // ------------------------------------------------------------------------------------------------------
    
    std::cout << "Regions loaded" << std::endl;

    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::cout << "saving... " << std::endl;

    // estimate the background in different Ias slices, each containing 10% of the statistic 
    // ------------------------------------------------------------------------------------------------------
    
    bckgEstimate(rb_50ias60, rc_ias50, rbc_50ias60, ra_ias50, rd_50ias60, "50ias60", nPE);
    bckgEstimate(rb_60ias70, rc_ias50, rbc_60ias70, ra_ias50, rd_60ias70, "60ias70", nPE);
    bckgEstimate(rb_70ias80, rc_ias50, rbc_70ias80, ra_ias50, rd_70ias80, "70ias80", nPE);
    bckgEstimate(rb_80ias90, rc_ias50, rbc_80ias90, ra_ias50, rd_80ias90, "80ias90", nPE);
    bckgEstimate(rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", nPE);
    bckgEstimate(rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100", nPE);
    
    // ------------------------------------------------------------------------------------------------------
    
    // bkg estimate for a selected cut index 
    // cutIndex = 3 --> pT > 60 GeV & Ias > 0.05
    //bckgEstimate_fromHistos(eta_cutIndex_regA, eta_cutIndex_regB, ih_eta_cutIndex_regB, eta_p_cutIndex_regC, cutIndex, true, nPE)->Write();

    //TODO not yet finished
    // 2D-histograms PredMass_Vs_CutIndex used in later steps to compare to data 
    int n_cutIndex = 1;
    TH2F* mass_cutIndex = new TH2F();

    for(int cutI = 0; cutI < n_cutIndex; cutI++){
        //bckgEstimate_fromHistos(eta_cutIndex_regA, eta_cutIndex_regB, ih_eta_cutIndex_regB, eta_p_cutIndex_regC, cutI, true, nPE);

    }

    delete ofile;
    delete mass_cutIndex;

    return;
}
