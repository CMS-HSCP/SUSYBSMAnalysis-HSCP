#ifndef SUSYBSMAnalysis_Analyzer_Regions_h
#define SUSYBSMAnalysis_Analyzer_Regions_h

#include "SUSYBSMAnalysis/Analyzer/interface/CommonFunction.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include <TCanvas.h>
#include <TLegend.h>
#include "TFile.h"
#include "TDirectory.h"

//data 2017
float K_data2017 = 2.30;
float C_data2017 = 3.17;
//data 2018
float K_data2018 = 2.27;
float C_data2018 = 3.16;
//MC 2017
float K_mc2017 = 2.26;
float C_mc2017 = 3.22;
//MC 2018
float K_mc2018 = 2.27;
float C_mc2018 = 3.22;

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
        void fillPredMass(const std::string&, float weight_);
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
        TH1F* pred_mass;
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
    TH1::SetDefaultSumw2(kTRUE);
    TH2::SetDefaultSumw2(kTRUE);
    TH3::SetDefaultSumw2(kTRUE);
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
    ih_p_eta = dir.make<TH3F>(("ih_p_eta"+suffix).c_str(),";#eta;p [GeV];I_{h} [MeV/cm]",neta,etalow,etaup,np,plow,pup,nih,ihlow,ihup); 
    eta_p = dir.make<TH2F>(("eta_p"+suffix).c_str(),";p [GeV];#eta",np,plow,pup,neta,etalow,etaup); 
    ih_eta = dir.make<TH2F>(("ih_eta"+suffix).c_str(),";#eta;I_{h} [MeV/cm]",neta,etalow,etaup,nih,ihlow,ihup); 
    ih_p = dir.make<TH2F>(("ih_p"+suffix).c_str(),";p [GeV];I_{h} [MeV/cm]",np,plow,pup,nih,ihlow,ihup);
    ias_p = dir.make<TH2F>(("ias_p"+suffix).c_str(),";p [GeV];I_{as}",np,plow,pup,nias,iaslow,iasup); 
    ias_pt = dir.make<TH2F>(("ias_pt"+suffix).c_str(),";pt [GeV];I_{as}",npt,ptlow,ptup,nias,iaslow,iasup);
    mass = dir.make<TH1F>(("mass"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); 
    pred_mass = dir.make<TH1F>(("pred_mass"+suffix).c_str(),";Mass [GeV]",nmass,masslow,massup); 
    mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    pred_mass->SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    pt_pterroverpt = dir.make<TH2F>(("pt_pterroverpt"+suffix).c_str(),";p_{T} [GeV];#frac{#sigma_{pT}}{p_{T}}",npt,ptlow,ptup,100,0,1); 
    hTOF    = dir.make<TH1F>(("hTOF_"+suffix).c_str(),";TOF",200,-10,10); 
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
void Region::fillPredMass(const std::string& st_sample,float weight_=-1) {
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
                //if(weight_>0) weight = weight_;
                float err_weight = weight*sqrt((1./(ih->GetBinContent(k)))+(1./(p->GetBinContent(j)*ih->Integral())));
                float K=0, C=0;
                if(st_sample=="data2017"){K=K_data2017;C=C_data2017;}
                if(st_sample=="data2018"){K=K_data2018;C=C_data2018;}
                if(st_sample=="mc2017"){K=K_mc2017;C=C_mc2017;}
                if(st_sample=="mc2018"){K=K_mc2018;C=C_mc2018;}
                float mass = GetMass(mom,dedx,K,C);
                int bin_mass = pred_mass->FindBin(mass);
                if(prob>=0)
                {
                    pred_mass->SetBinContent(bin_mass,pred_mass->GetBinContent(bin_mass)+weight);
                    pred_mass->SetBinError(bin_mass,sqrt(pow(pred_mass->GetBinError(bin_mass),2)+pow(err_weight,2)));
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
    pred_mass->Write();
    pt_pterroverpt->Write();
    hTOF->Write();
}

void loadHistograms(Region& r, TFile* f, const std::string& regionName, bool bool_rebin=true, int rebineta=1, int rebinp=1, int rebinih=1, int rebinmass=1){
    //std::string dir = "analyzer/BaseName/";
    std::string dir = "HSCParticleAnalyzer/BaseName/";
    r.ih_p_eta                          = (TH3F*)f->Get((dir+"ih_p_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p_eta->Rebin3D(rebineta,rebinp,rebinih);
    r.eta_p                             = (TH2F*)f->Get((dir+"eta_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.eta_p->Rebin2D(rebinp,rebineta);
    r.ih_eta                            = (TH2F*)f->Get((dir+"ih_eta_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_eta->Rebin2D(rebineta,rebinih);
    r.ih_p                              = (TH2F*)f->Get((dir+"ih_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ih_p->Rebin2D(rebinp,rebinih);
    r.ias_p                             = (TH2F*)f->Get((dir+"ias_p_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_p->Rebin2D(rebinp,rebinih);
    r.ias_pt                            = (TH2F*)f->Get((dir+"ias_pt_"+regionName).c_str())->Clone(); if(bool_rebin) r.ias_pt->Rebin2D(rebinp,rebinih);
    r.mass                              = (TH1F*)f->Get((dir+"mass_"+regionName).c_str())->Clone(); if(bool_rebin) r.mass->Rebin(rebinmass);
    //r.mass                              = (TH1F*)f->Get((dir+"massFromTree_"+regionName).c_str())->Clone(); if(bool_rebin) r.mass->Rebin(rebinmass);
    r.pred_mass                         = (TH1F*)f->Get((dir+"pred_mass_"+regionName).c_str())->Clone(); r.pred_mass->Reset(); if(bool_rebin) r.pred_mass->Rebin(rebinmass);
    //r.pred_mass                         = (TH1F*)f->Get((dir+"massFrom1DTemplatesEtaBinning_"+regionName).c_str())->Clone(); r.pred_mass->Reset(); if(bool_rebin) r.pred_mass->Rebin(rebinmass);
}

// Return randomly select histo 
void poissonHisto(TH1F &h,TRandom3* RNG){
    for(int i=0;i<h.GetNbinsX()+1;i++){
        h.SetBinContent(i,RNG->Poisson(h.GetBinContent(i)));
    }
}
void poissonHisto(TH2F &h,TRandom3* RNG){
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
    TH1F* eta1 = (TH1F*) eta_p_1->ProjectionY(); 
    eta1->Scale(1./eta1->Integral());
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
    //double xbins[17] = {0,50,100,150,200,250,300,350,400,450,500,600,700,800,1000,1500,2000};
    double xbins[12] = {0,100,200,300,400,500,600,700,800,1000,1500,2000};
    std::vector<double> xbins_v;
    for(double i=0.0;i<=1000.0;i+=50) xbins_v.push_back(i);
    std::string newname = h->GetName(); 
    newname += "_rebinned";
    TH1F* hres = (TH1F*) h->Rebin(11,newname.c_str(),xbins);
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

TH1F* pull(TH1F* h1, TH1F* h2){
    float SystError = 0.2;
    TH1F* res = (TH1F*) h1->Clone(); res->Reset();
    for(int i=0;i<h1->GetNbinsX()+1;i++){
        double Perr = 0, Derr = 0;
        double P = h1->GetBinContent(i); if(P<=0) continue;
        double D = h2->GetBinContent(i); if(D<=0) continue;
        Perr = sqrt(P + pow(P*SystError,2));
        Derr = sqrt(D);
        res->SetBinContent(i,(D-P)/sqrt(pow(Derr,2)+pow(Perr,2)));
        res->SetBinError(i,res->GetBinContent(i)*((Derr/D)+(Perr/P)));
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
    float SystError = 0.2;
    TH1F h(vPE[0]);
    h.Reset();
    h.SetBinErrorOption(TH1::EBinErrorOpt::kPoisson);
    for(int i=0;i<h.GetNbinsX()+1;i++){
        float mean=0, err=0;
        TH1F* htemp = new TH1F("htemp","htemp",100,0,1e6);
        for(unsigned int pe=0;pe<vPE.size();pe++){
            mean += vPE[pe].GetBinContent(i);
            htemp->Fill(vPE[pe].GetBinContent(i));
        }
        mean /= vPE.size();
        for(unsigned int pe=0;pe<vPE.size();pe++){
            err += pow(mean - vPE[pe].GetBinContent(i),2);
        }
        float fact=1;
        if(vPE.size()>1) {err = sqrt(err/(vPE.size()-1));fact=vPE.size()/(vPE.size()-1);}
        else err = sqrt(err);
        mean = htemp->GetMean();
        err = fact*htemp->GetStdDev();
        err = sqrt(pow(err,2)+pow(SystError*mean,2));
        h.SetBinContent(i,mean);
        h.SetBinError(i,err);
        delete htemp;
    }
    return h;
}



// Function returning a canvas divide in two windows
// The first one contains the two 1D-histograms, given as arguments, superposed.
// There also is a legend associated to this window where the names are defined as arguments.
// The second window contains the ratio of these 1D-histograms or the ratio of right integers of them. 
// We define which kind of ratio we want with tha 'ratioSimple' boolean.
// The 'name' given corresponds to the name of the canvas 
TCanvas* plotting(TH1F* h1, TH1F* h2, bool ratioSimple=true, std::string dirname="", std::string name="", std::string leg1="", std::string leg2="", bool rebin=false){
    if(rebin) h1=rebinHisto(h1);
    if(rebin) h2=rebinHisto(h2);
    std::string canvName;
    if(rebin) {
      canvName = "plotting_"+name+"_rebin";
    } else {
      canvName = "plotting_"+name;
    }
    TCanvas* c1 = new TCanvas(canvName.c_str(),"", 800,800);
    c1->Divide(1,3);
    gStyle->SetOptStat(0);
    //c1->cd(1);
    //TPad* p1 = (TPad*)(c1->cd(1));
    TPad *p1 = new TPad("pad1", "pad1", 0, 0.3, 1, 1.0);
    p1->Draw();
    p1->cd();
    p1->SetLogy();
    p1->SetTopMargin(0.05);
    p1->SetBottomMargin(0.1);
    p1->SetLeftMargin(0.12);
    p1->SetRightMargin(0.05);
    TLegend* leg = new TLegend(0.7,0.7,0.9,0.9);
    leg->AddEntry(h1,leg1.c_str(),"lep");
    leg->AddEntry(h2,leg2.c_str(),"lep");
    h1->SetStats(0);
    h1->SetMarkerStyle(20);
    h1->Draw();
    h1->GetYaxis()->SetRangeUser(1e-4,1e6);
    h2->SetLineColor(2);
    h2->SetStats(0);
    h2->SetMarkerStyle(20);
    h2->SetMarkerColor(2);
    h2->Draw("esame");
    leg->Draw("same");
    h1->SetName((name+"_obs").c_str());
    h1->Write();
    h1->SetName((name+"_pred").c_str());
    h2->Write();
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
    tmp->SetLineColor(1);
    tmp->SetMarkerColor(1);
    tmp->SetName((name+"_ratioInt").c_str());
    tmp->Write();
    c1->cd(3);
    TH1F* tmp2 = (TH1F*) pull(h2,h1)->Clone();
    tmp2->Draw("E0");
    tmp2->GetYaxis()->SetTitle("pulls");
    tmp2->GetYaxis()->SetRangeUser(-3,3);
    tmp2->SetLineColor(1);
    tmp2->SetMarkerColor(1);
    tmp2->SetName((name+"_pull").c_str());
    tmp2->Write();
    c1->SaveAs((dirname+canvName+".pdf").c_str());
    c1->SaveAs((dirname+canvName+".root").c_str());
    return c1;
}

void bckgEstimate(const std::string& st_sample, const std::string& dirname, const Region& B, const Region& C, const Region& BC, const Region& A, const Region& D, const std::string& st, const int& nPE=100, const int& rebinMass=1){
    Region bc = BC;
    Region d = D;
    std::vector<TH1F> vPE;
    TRandom3* RNG = new TRandom3();
    for(int pe=0;pe<nPE;pe++){
        Region a = A;
        Region b = B;
        Region c = C;
        TH2F a_ih_eta(*a.ih_eta);
        TH2F b_ih_eta(*b.ih_eta);
        TH2F c_ih_eta(*c.ih_eta);
        TH2F b_eta_p(*b.eta_p);
        TH2F c_eta_p(*c.eta_p);
        poissonHisto(a_ih_eta,RNG);
        poissonHisto(b_ih_eta,RNG);
        poissonHisto(c_ih_eta,RNG);
        poissonHisto(b_eta_p,RNG);
        poissonHisto(c_eta_p,RNG);
        TH1F* b_eta = (TH1F*)b_eta_p.ProjectionY();
        etaReweighingP(&c_eta_p,b_eta);
        bc.eta_p = &c_eta_p;    bc.ih_eta = &b_ih_eta;
        float A = a_ih_eta.Integral();
        float B = b_ih_eta.Integral();
        float C = c_ih_eta.Integral();
        float normalisationABC = 1.;
        if(A>0) normalisationABC = B*C/A;
        bc.fillPredMass(st_sample);
        scale(bc.pred_mass);
        massNormalisation(bc.pred_mass,normalisationABC);
        vPE.push_back(*bc.pred_mass);
    }
    TH1F h_temp = meanHistoPE(vPE);
    if(nPE>1) bc.pred_mass = &h_temp;
    
    saveHistoRatio(d.mass,bc.pred_mass,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str());
    saveHistoRatio(d.mass,bc.pred_mass,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str(),true);
    
    overflowLastBin(d.mass);
    overflowLastBin(bc.pred_mass);

    d.mass->Rebin(rebinMass);
    bc.mass->Rebin(rebinMass);
    
    plotting(d.mass,bc.pred_mass,false,dirname,("mass1D_regionBC_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction")->Write();
    plotting(d.mass,bc.pred_mass,false,dirname,("mass1D_regionBC_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction",true)->Write();
    delete RNG;
}

void bckgEstimate_fromHistos(const std::string& st_sample, const std::string& dirname, const TH2F& mass_cutInd, const TH2F& eta_cutIndex_A, const TH2F& eta_cutIndex_B, const TH3F& ih_eta_cutIndex_B, const TH3F& eta_p_cutIndex_C, const TH1F& HA, const TH1F& HB, const TH1F& HC, int cutIndex=3, int nPE=100){
    TH2F* mass_cutIndex = (TH2F*) mass_cutInd.Clone();
    TH1F* mass_obs = (TH1F*)mass_cutIndex->ProjectionY("_projD",cutIndex+1,cutIndex+1);
    Region rBC;
    rBC.pred_mass = (TH1F*)mass_obs->Clone();
    rBC.pred_mass->Reset();
    std::vector<TH1F> vPE;
    TRandom3* RNG = new TRandom3();
    for(int pe=0;pe<nPE;pe++){
   
        TH2F* eta_cutIndex_regA = (TH2F*) eta_cutIndex_A.Clone();
        TH2F* eta_cutIndex_regB = (TH2F*) eta_cutIndex_B.Clone();
        TH3F* ih_eta_cutIndex_regB = (TH3F*) ih_eta_cutIndex_B.Clone();
        TH3F* eta_p_cutIndex_regC = (TH3F*) eta_p_cutIndex_C.Clone();
        TH1F* H_A = (TH1F*) HA.Clone();
        TH1F* H_B = (TH1F*) HB.Clone();
        TH1F* H_C = (TH1F*) HC.Clone();
    
        TH1F* eta_regA = (TH1F*)eta_cutIndex_regA->ProjectionY("_projA",cutIndex+1,cutIndex+1);
        TH1F* eta_regB = (TH1F*)eta_cutIndex_regB->ProjectionY("_projB",cutIndex+1,cutIndex+1);
        ih_eta_cutIndex_regB->GetXaxis()->SetRange(cutIndex+1,cutIndex+1);
        TH2F* ih_eta_regB =  (TH2F*)ih_eta_cutIndex_regB->Project3D("zyB");
        eta_p_cutIndex_regC->GetXaxis()->SetRange(cutIndex+1,cutIndex+1);
        TH2F* eta_p_regC = (TH2F*)eta_p_cutIndex_regC->Project3D("yzC"); 

        poissonHisto(*eta_regA,RNG);
        poissonHisto(*eta_regB,RNG);
        poissonHisto(*ih_eta_regB,RNG);
        poissonHisto(*eta_p_regC,RNG);
        poissonHisto(*H_A,RNG);
        poissonHisto(*H_B,RNG);
        poissonHisto(*H_C,RNG);
        etaReweighingP(eta_p_regC, eta_regB);
        rBC.eta_p = eta_p_regC; rBC.ih_eta = ih_eta_regB;
        float A = H_A->Integral(cutIndex+1,cutIndex+1);
        float B = H_B->Integral(cutIndex+1,cutIndex+1);
        float C = H_C->Integral(cutIndex+1,cutIndex+1);
        float norm = 1;
        if(A>0) norm = B*C/A;
        rBC.fillPredMass(st_sample);
        scale(rBC.pred_mass);
        massNormalisation(rBC.pred_mass,norm);
        vPE.push_back(*rBC.pred_mass);
    }
    TH1F h_tmp = meanHistoPE(vPE);
    if(nPE>1) rBC.pred_mass = &h_tmp;

    std::string st = "cutIndex"+to_string(cutIndex);
    
    saveHistoRatio(mass_obs,rBC.pred_mass,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str());
    saveHistoRatio(mass_obs,rBC.pred_mass,("mass_obs_"+st).c_str(),("mass_predBC_"+st).c_str(),("mass_predBCR_"+st).c_str(),true);
    
    overflowLastBin(mass_obs);
    overflowLastBin(rBC.pred_mass);
    
    plotting(mass_obs,rBC.pred_mass,false,dirname,("mass1D_regionBC_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction")->Write();
    plotting(mass_obs,rBC.pred_mass,false,dirname,("mass1D_regionBC_"+st+"_nPE-"+to_string(nPE)).c_str(),"Observed","Prediction",true)->Write();

    delete RNG;
}

#endif
