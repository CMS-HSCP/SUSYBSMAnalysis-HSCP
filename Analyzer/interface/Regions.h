#ifndef SUSYBSMAnalysis_Analyzer_Regions_h
#define SUSYBSMAnalysis_Analyzer_Regions_h

#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TCanvas.h>
#include <TLegend.h>
#include "TFile.h"
#include "TDirectory.h"


TRandom3* RNG = new TRandom3();

float K = 1;
float C = 1;

// Scale the 1D-histogram given to the unit 
void scale(TH1F* h){
    h->Scale(1./h->Integral(0,h->GetNbinsX()+1));
}

// class using to definite signal and control regions. 
class Region{
    public:
        Region();
        Region(TFileDirectory &dir,std::string suffix,int& etabins,int& ihbins,int& pbins,int& massbins);
        ~Region();
        void setSuffix(std::string suffix);
        void initHisto(TFileDirectory &dir,int etabins,int ihbins,int pbins,int massbins);
        void fill(float& eta, float&p, float& pt, float& pterr, float& ih, float& ias, float& m, float& tof, float& w);
        void fillMassFrom1DTemplatesEtaBinning(float weight_);
        void write();

        int np;
        float plow;
        float pup;
        int npt;
        float ptlow;
        float ptup;
        int nih;
        float ihlow;
        float ihup;
        int nias;
        float iaslow;
        float iasup;
        int neta;
        float etalow;
        float etaup;
        int nmass;
        float masslow;
        float massup;
        std::vector<double> VectOfBins_P_;
        std::string suffix_;
        TH3F* ih_p_eta;
        TH2F* eta_p;
        TH2F* ih_eta;
        TH2F* ih_p;
        TH2F* ias_p;
        TH2F* ias_pt;
        TH1F* mass;
        TH1F* massFrom1DTemplatesEtaBinning;
        TH2F* eta_p_rebinned;
        TH2F* pt_pterroverpt;
        TH1F* hTOF;
};

Region::Region(){}

Region::Region(TFileDirectory &dir, std::string suffix,int& etabins,int& ihbins,int& pbins,int& massbins){
    suffix_ = suffix;
    initHisto(dir,etabins,ihbins,pbins,massbins);
} 

Region::~Region(){}

void Region::setSuffix(std::string suffix){
    suffix_ = suffix;
}

// Function which intializes the histograms with given binnings 
void Region::initHisto(TFileDirectory &dir,int etabins,int ihbins,int pbins,int massbins){
    np = pbins;
    plow = 0;
    pup = 10000;
    npt = pbins;
    ptlow = 0;
    ptup = 10000; 
    nih = ihbins;
    ihlow = 0;
    ihup = 20;
    nias = ihbins;
    iaslow = 0;
    iasup = 1;
    neta = etabins;
    etalow = -3;
    etaup = 3;
    nmass = massbins;
    masslow = 0;
    massup = 4000;
    std::string suffix = suffix_;
    ih_p_eta = dir.make<TH3F>(("ih_p_eta"+suffix).c_str(),";#eta;p [GeV];I_{h} [MeV/cm]",neta,etalow,etaup,np,plow,pup,nih,ihlow,ihup); ih_p_eta->Sumw2();
    eta_p = dir.make<TH2F>(("eta_p"+suffix).c_str(),";p [GeV];#eta",np,plow,pup,neta,etalow,etaup); eta_p->Sumw2();
    ih_eta = dir.make<TH2F>(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup); ih_eta->Sumw2();
    ih_p = dir.make<TH2F>(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup); ih_p->Sumw2();
    ias_p = dir.make<TH2F>(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",np,plow,pup,nias,iaslow,iasup); ias_p->Sumw2();
    ias_pt = dir.make<TH2F>(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup); ias_pt->Sumw2();
    mass = dir.make<TH1F>(("massFromTree"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); mass->Sumw2();
    massFrom1DTemplatesEtaBinning = dir.make<TH1F>(("massFrom1DTemplatesEtaBinning"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); massFrom1DTemplatesEtaBinning->Sumw2();
    mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    massFrom1DTemplatesEtaBinning->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    pt_pterroverpt = dir.make<TH2F>(("pt_pterroverpt"+suffix).c_str(),";p_{T} [GeV];#frac{#sigma_{pT}}{p_{T}}",npt,ptlow,ptup,100,0,1); pt_pterroverpt->Sumw2();
    hTOF    = dir.make<TH1F>(("hTOF_"+suffix).c_str(),";TOF",200,-10,10); hTOF->Sumw2();
}

// Function which fills histograms
void Region::fill(float& eta, float& p, float& pt, float& pterr, float& ih, float& ias, float& m, float& tof, float& w){
   ih_p_eta->Fill(eta,p,ih,w);
   eta_p->Fill(p,eta,w);
   ih_eta->Fill(eta,ih,w);
   ih_p->Fill(p,ih,w);
   ias_p->Fill(p,ias,w);
   ias_pt->Fill(pt,ias,w);
   mass->Fill(m,w);
   pt_pterroverpt->Fill(pt,pterr/pt,w);
   hTOF->Fill(tof,w);
}

// in order to compute properly the uncertainties we use the methods SetBinContent SetBinError instead of Fill
// as several couples of bins in (p,ih) can provide the same mass estimate we need to properly sum the entries and errors
// for a couple of bins in (p,ih) where the bin content were (N_p,N_ih) the associated quantities should be 
// content: (N_p * N_ih) / N_total, where N_total represents the total number of events in the region (integral of p, ih & mass distributions)
// error: content * sqrt( 1 / N_p + 1 / N_ih ) where we assume Poisson uncertainties in both distributions (independent distributions) and we neglect the uncertainty on N_total
// While combining the input for several couples leading to the same mass: 
// contents are added 
// errors: the sqrt of the squared uncertainties are added
void Region::fillMassFrom1DTemplatesEtaBinning(float weight_=-1) {
    TH1F* eta = (TH1F*) ih_eta->ProjectionX();
    for(int i=1;i<eta->GetNbinsX();i++)
    {
        TH1F* p = (TH1F*) eta_p->ProjectionX("proj_p",i,i);
        if(VectOfBins_P_.size()>1) p = (TH1F*)p->Rebin(VectOfBins_P_.size()-1,"",VectOfBins_P_.data());
        TH1F* ih = (TH1F*) ih_eta->ProjectionY("proj_ih",i,i);
        scale(p); //only scale one of the two distributions ih or p --> keep the information of the normalisation 
        for(int j=1;j<p->GetNbinsX();j++)
        {
            for(int k=1;k<ih->GetNbinsX();k++)
            {
                if(p->GetBinContent(j)<=0) continue;
                if(ih->GetBinContent(k)<=0) continue;
                float mom = p->GetBinCenter(j);
                float dedx = ih->GetBinCenter(k);
                float prob = p->GetBinContent(j) * ih->GetBinContent(k);
                float weight = prob;
                if(weight_>0) weight = weight_;
                float err_weight = weight*sqrt((1./(ih->GetBinContent(k)))+(1./(p->GetBinContent(j)*ih->Integral())));
                float mass = GetMass(mom,dedx,K,C);
                int bin_mass = massFrom1DTemplatesEtaBinning->FindBin(mass);
                if(prob>=0)
                {
                    massFrom1DTemplatesEtaBinning->SetBinContent(bin_mass,massFrom1DTemplatesEtaBinning->GetBinContent(bin_mass)+weight);
                    massFrom1DTemplatesEtaBinning->SetBinError(bin_mass,sqrt(pow(massFrom1DTemplatesEtaBinning->GetBinError(bin_mass),2)+pow(err_weight,2)));
                }
            }
        }
        delete p;
        delete ih;
    }
}

void Region::write(){
    ih_p_eta->Write();
    eta_p->Write();
    ih_eta->Write();
    ih_p->Write();
    ias_p->Write();
    ias_pt->Write();
    mass->Write();
    massFrom1DTemplatesEtaBinning->Write();
    pt_pterroverpt->Write();
    hTOF->Write();
}



void loadHistograms(Region& r, TFile* f, const std::string& regionName, bool bool_rebin=true, int rebineta=1, int rebinp=1, int rebinih=1, int rebinmass=1){
    r.ih_p_eta                          = (TH3F*)f->Get(("ih_p_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p_eta->Rebin3D(rebineta,rebinp,rebinih);
    r.eta_p                             = (TH2F*)f->Get(("eta_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.eta_p->Rebin2D(rebinp,rebineta);
    r.ih_eta                            = (TH2F*)f->Get(("ih_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_eta->Rebin2D(rebineta,rebinih);
    r.ih_p                              = (TH2F*)f->Get(("ih_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p->Rebin2D(rebinp,rebinih);
    r.ias_p                             = (TH2F*)f->Get(("ias_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_p->Rebin2D(rebinp,rebinih);
    r.ias_pt                            = (TH2F*)f->Get(("ias_pt_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_pt->Rebin2D(rebinp,rebinih);
    r.mass                              = (TH1F*)f->Get(("massFromTree_"+regionName).c_str())->Clone(); if(bool_rebin) r.mass->Rebin(rebinmass);
    r.massFrom1DTemplatesEtaBinning     = (TH1F*)f->Get(("massFrom1DTemplatesEtaBinning_"+regionName).c_str())->Clone(); r.massFrom1DTemplatesEtaBinning->Reset(); if(bool_rebin) r.massFrom1DTemplatesEtaBinning->Rebin(rebinmass);
}



// Return randomly select histo 
void poissonHisto(TH2F &h){
    for(int i=0;i<h.GetNbinsX()+1;i++){
        for(int j=0;j<h.GetNbinsY()+1;j++){
            h.SetBinContent(i,j,RNG->Poisson(h.GetBinContent(i,j)));
        }
    }
}
// Function doing the eta reweighing between two 2D-histograms as done in the Hscp background estimate method,
// because of the correlation between variables (momentum & transverse momentum). 
// The first given 2D-histogram is weighted in respect to the 1D-histogram 
void etaReweighingP(TH2F* eta_p_1, TH1F* eta2)
{
    TH1F* eta1 = (TH1F*) eta_p_1->ProjectionY(); eta1->Scale(1./eta1->Integral());
    eta2->Scale(1./eta2->Integral());
    eta2->Divide(eta1);
    for(int i=0;i<eta_p_1->GetNbinsX()+1;i++)
    {
        for(int j=0;j<eta_p_1->GetNbinsY()+1;j++)
        {
            float val_ij = eta_p_1->GetBinContent(i,j);
            float err_ij = eta_p_1->GetBinError(i,j);
            eta_p_1->SetBinContent(i,j,val_ij*eta2->GetBinContent(j));
            eta_p_1->SetBinError(i,j,err_ij*eta2->GetBinContent(j));
        }
    }
    eta_p_1->Sumw2();
}


// normalizes an histogram to a given norm 
void massNormalisation(TH1F* h, const float& normalisation){
    for(int k=0;k<h->GetNbinsX()+1;k++){
        h->SetBinContent(k,h->GetBinContent(k)*normalisation);
        h->SetBinError(k,h->GetBinError(k)*normalisation);
    }
}

// add the overflow bin to the last one
void overflowLastBin(TH1F* h){
    h->SetBinContent(h->GetNbinsX(),h->GetBinContent(h->GetNbinsX())+h->GetBinContent(h->GetNbinsX()+1));
    h->SetBinContent(h->GetNbinsX()+1,0);
}

// rebinning histogram according to an array of bins
TH1F* rebinHisto(TH1F* h){
    overflowLastBin(h);
    double xbins[17] = {0,50,100,150,200,250,300,350,400,450,500,600,700,800,1000,1500,2000};
    std::vector<double> xbins_v;
    for(double i=0.0;i<=1000.0;i+=50) xbins_v.push_back(i);
    std::string newname = h->GetName(); 
    newname += "_rebinned";
    TH1F* hres = (TH1F*) h->Rebin(16,newname.c_str(),xbins);
    overflowLastBin(hres);
    return hres;
}

// Function returning the ratio of right integer (from x to infty) for two 1D-histograms
// This function is used in the Hscp data-driven background estimate to test the mass shape prediction
// The argument to use this type of ratio is that we're in case of cut & count experiment 
TH1F* ratioIntegral(TH1F* h1, TH1F* h2){    
    float SystError = 0.2;
    TH1F* res = (TH1F*) h1->Clone(); res->Reset();
    for(int i=0;i<h1->GetNbinsX()+1;i++)
    {   
        double Perr=0, Derr=0;
        double P=h1->IntegralAndError(i,h1->GetNbinsX()+1,Perr); if(P<=0) continue;
        double D=h2->IntegralAndError(i,h2->GetNbinsX()+1,Derr);
        Perr = sqrt(Perr*Perr + pow(P*SystError,2));
        res->SetBinContent(i,D/P);
        res->SetBinError(i,sqrt(pow(Derr*P,2)+pow(Perr*D,2))/pow(P,2));
    }
    return res;
}

void saveHistoRatio(TH1F* h1,TH1F* h2,std::string st1,std::string st2,std::string st3,bool rebin=false){
    h1->SetName(st1.c_str());
    h2->SetName(st2.c_str());
    if(rebin){
        h1=rebinHisto(h1);
        h2=rebinHisto(h2);
    }
    h1->Write();
    h2->Write();
    TH1F* R = (TH1F*) ratioIntegral(h2,h1)->Clone();
    if(rebin) st3+="_rebinned";
    R->SetName(st3.c_str());
    R->Write();
}

TH1F meanHistoPE(std::vector<TH1F> vPE){
    TH1F h(vPE[0]);
    h.Reset();
    h.SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    for(int i=0;i<h.GetNbinsX()+1;i++){
        float mean=0, err=0;
        for(unsigned int pe=0;pe<vPE.size();pe++){
            mean += vPE[pe].GetBinContent(i);
        }
        mean /= vPE.size();
        for(unsigned int pe=0;pe<vPE.size();pe++){
            err += pow(mean - vPE[pe].GetBinContent(i),2);
        }
        err = sqrt(err/(vPE.size()-1));
        h.SetBinContent(i,mean);
        h.SetBinError(i,err);
    }
    return h;
}



// Function returning a canvas divide in two windows
// The first one contains the two 1D-histograms, given as arguments, superposed.
// There also is a legend associated to this window where the names are defined as arguments.
// The second window contains the ratio of these 1D-histograms or the ratio of right integers of them. 
// We define which kind of ratio we want with tha 'ratioSimple' boolean.
// The 'name' given corresponds to the name of the canvas 
TCanvas* plotting(TH1F* h1, TH1F* h2, bool ratioSimple=true, std::string name="", std::string leg1="", std::string leg2="", bool rebin=false){
    if(rebin) h1=rebinHisto(h1);
    if(rebin) h2=rebinHisto(h2);
    TCanvas* c1 = new TCanvas(("plotting_"+name).c_str(),"");
    c1->Divide(1,2);
    gStyle->SetOptStat(0);
    c1->cd(1);
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(h1,leg1.c_str(),"lep");
    leg->AddEntry(h2,leg2.c_str(),"lep");
    h1->Draw();
    h2->SetLineColor(2);
    h2->Draw("esame");
    leg->Draw("same");
    c1->SetLogy();
    c1->cd(2);
    TH1F* tmp = (TH1F*) h1->Clone(); tmp->Reset();
    if(ratioSimple){
        tmp = (TH1F*)h1->Clone();
        tmp->Divide(h2);
        tmp->GetYaxis()->SetTitle("#frac{N_{obs}}{N_{pred}}");
        tmp->GetYaxis()->SetTitleSize(0.06);
    }else{
        tmp=ratioIntegral(h2,h1);
        tmp->GetYaxis()->SetTitle("#int_{M}^{#infty} dm_{obs} / #int_{M}^{#infty} dm_{pred}");
        tmp->GetYaxis()->SetTitleSize(0.06);
    }
    tmp->GetYaxis()->SetRangeUser(0,2);
    tmp->Draw();
    return c1;
}

void bckgEstimate(Region& b, Region& c, Region& bc, Region& a, Region& d, std::string st, int nPE=100){
    std::vector<TH1F> vPE;
    std::cout << st << std::endl;
    for(int pe=0;pe<nPE;pe++){
        TH2F a_ih_eta(*a.ih_eta);
        TH2F b_ih_eta(*b.ih_eta);
        TH2F c_ih_eta(*c.ih_eta);
        TH2F b_eta_p(*b.eta_p);
        TH2F c_eta_p(*c.eta_p);
        poissonHisto(a_ih_eta);
        poissonHisto(b_ih_eta);
        poissonHisto(c_ih_eta);
        poissonHisto(b_eta_p);
        poissonHisto(c_eta_p);
        TH1F* b_eta = (TH1F*)b_eta_p.ProjectionY();
        etaReweighingP(&c_eta_p,b_eta);
        bc.eta_p = &c_eta_p;    bc.ih_eta = &b_ih_eta;
        float A = a_ih_eta.Integral();
        float B = b_ih_eta.Integral();
        float C = c_ih_eta.Integral();
        float normalisationABC = B*C/A;
        bc.fillMassFrom1DTemplatesEtaBinning();
        scale(bc.massFrom1DTemplatesEtaBinning);
        massNormalisation(bc.massFrom1DTemplatesEtaBinning,normalisationABC);
        vPE.push_back(*bc.massFrom1DTemplatesEtaBinning);
    }
    TH1F h_temp = meanHistoPE(vPE);
    if(nPE>1) bc.massFrom1DTemplatesEtaBinning = &h_temp;
    
    saveHistoRatio(d.mass,bc.massFrom1DTemplatesEtaBinning,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str());
    saveHistoRatio(d.mass,bc.massFrom1DTemplatesEtaBinning,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str(),true);
    
    overflowLastBin(d.mass);
    overflowLastBin(bc.massFrom1DTemplatesEtaBinning);
    
    plotting(d.mass,bc.massFrom1DTemplatesEtaBinning,false,("mass1D_regionBC_"+st).c_str(),"Observed","Prediction")->Write();
    plotting(d.mass,bc.massFrom1DTemplatesEtaBinning,false,("mass1D_regionBC_"+st).c_str(),"Observed","Prediction",true)->Write();
}

#endif
