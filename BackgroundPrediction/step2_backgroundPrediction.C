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
    int rebineta,rebinih,rebinp,rebinmass,thresholdP,thresholdMass;
    bool rebin,varBinsP,varBinsMass;
    while(std::getline(infile,line)){
        if(std::strncmp(line.c_str(),"#",1)==0) continue;
        std::cout << line << std::endl;
        std::stringstream ss(line);
        ss >> filename >> rebin >> rebineta >> rebinih >> rebinp >> rebinmass;
    }

    std::string outfilename_;
    outfilename_ = filename+"_rebinEta"+to_string(rebineta)+"_rebinIh"+to_string(rebinih)+"_rebinP"+to_string(rebinp)+"_rebinMass"+to_string(rebinmass)+"_analysed";

    std::cout << outfilename_ << std::endl;

    TFile* ifile = new TFile((filename+".root").c_str());

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
    
    bool bool_rebin=rebin;
    
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

    std::cout << "Regions loaded" << std::endl;

    TFile* ofile = new TFile((outfilename_+".root").c_str(),"RECREATE");

    std::cout << "saving... " << std::endl;

    bckgEstimate(rb_50ias60, rc_ias50, rbc_50ias60, ra_ias50, rd_50ias60, "50ias60", 100);
    bckgEstimate(rb_60ias70, rc_ias50, rbc_60ias70, ra_ias50, rd_60ias70, "60ias70", 100);
    bckgEstimate(rb_70ias80, rc_ias50, rbc_70ias80, ra_ias50, rd_70ias80, "70ias80", 100);
    bckgEstimate(rb_80ias90, rc_ias50, rbc_80ias90, ra_ias50, rd_80ias90, "80ias90", 100);
    bckgEstimate(rb_50ias90, rc_ias50, rbc_50ias90, ra_ias50, rd_50ias90, "50ias90", 100);
    bckgEstimate(rb_90ias100, rc_ias50, rbc_90ias100, ra_ias50, rd_90ias100, "90ias100", 100);

    return;
}
