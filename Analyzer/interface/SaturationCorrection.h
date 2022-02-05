#ifndef SUSYBSMAnalysis_Analyzer_SaturationCorrection_h
#define SUSYBSMAnalysis_Analyzer_SaturationCorrection_h

#include <iostream>
#include <fstream>
#include <vector>
#include <stdlib.h>
#include "TFile.h"
#include "TTree.h"
#include "TBranch.h"
#include "TChain.h"
#include "TF1.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "THStack.h"
#include "TProfile.h"
#include "TFitResult.h"
#include "TCanvas.h"
#include "TLine.h"
#include "TStyle.h"
#include "TEfficiency.h"
#include "TText.h"
#include "TLegend.h"
#include "TMath.h"
#include "Rtypes.h"
#include "TLatex.h"
#include "TMathText.h"

using namespace std;

string LabelModul(int modulgeom) {
  if (modulgeom == 1)
    return "IB1";
  if (modulgeom == 2)
    return "IB2";
  if (modulgeom == 3)
    return "OB1";
  if (modulgeom == 4)
    return "OB2";
  if (modulgeom == 5)
    return "W1A";
  if (modulgeom == 6)
    return "W2A";
  if (modulgeom == 7)
    return "W3A";
  if (modulgeom == 8)
    return "W1B";
  if (modulgeom == 9)
    return "W2B";
  if (modulgeom == 10)
    return "W3B";
  if (modulgeom == 11)
    return "W4";
  if (modulgeom == 12)
    return "W5";
  if (modulgeom == 13)
    return "W6";
  if (modulgeom == 14)
    return "W7";
  else
    return "UNKNOWN";
  return "";
}

class SaturationCorrection {
private:
  TFile* file_;
  TTree* tree_;
  TTree* inputtree_;
  float p0_;
  float p1_;
  float p0err_;
  float p1err_;
  float chi2_;
  int label_;
  int nstrip_;
  int nsat_;
  TH2F* histo_;
  TH1F* historec_;

  float p0[15][4][2];
  float p1[15][4][2];
  float chi2[15][4][2];

  vector<vector<vector<TH2F*>>> vModulGeomvNStripvNSat;

public:
  SaturationCorrection();
  ~SaturationCorrection();
  void InitHisto();
  void FillHisto(int charge_sim, int charge_reco, int label, int nstrip, int nsat);
  void SetFileAndTreeName(string file_name, string tree_name);
  void SetHisto(TH2F& histo);
  void FillHistoRec();
  void SetTree(TTree& tree);
  void SetBranch();
  void Write(int label, int nstrip, int nsat);
  void WriteFile();
  void ReadParameters();
  void FillProfile();
  void WriteParameters();
  float ChargeCorrected(float charge, int label, int nstrip, int nsat);
};

SaturationCorrection::SaturationCorrection() {
  p0_ = 0.;
  p1_ = 1.;
  p0err_ = 0.;
  p1err_ = 0.;
  chi2_ = 0.;
  label_ = 0;
  nstrip_ = 0;
  nsat_ = 0;
  for (int label = 1; label < 15; label++) {
    for (int nstrip = 3; nstrip < 7; nstrip++) {
      for (int nsat = 0; nsat < 3; nsat++) {
        if (label_ == label && nstrip_ == nstrip && nsat_ == nsat) {
          p0[label - 1][nstrip - 3][nsat] = 0.;
          p1[label - 1][nstrip - 3][nsat] = 1.;
          chi2[label - 1][nstrip - 3][nsat] = 0.;
        }
      }
    }
  }
}

SaturationCorrection::~SaturationCorrection() {}

void SaturationCorrection::InitHisto() {
  for (int ModulGeom = 1; ModulGeom < 15; ModulGeom++) {
    vector<vector<TH2F*>> vNStripvNSat;
    for (int nstrip = 3; nstrip < 7; nstrip++) {
      vector<TH2F*> vNSat;
      for (int nsat = 1; nsat < 3; nsat++) {
        string title = LabelModul(ModulGeom) + " NStrip=" + to_string(nstrip) + " NSat=" + to_string(nsat);
        vNSat.push_back(new TH2F(title.c_str(), title.c_str(), 200, 0, 5000, 200, 0, 5000));
      }
      vNStripvNSat.push_back(vNSat);
    }
    vModulGeomvNStripvNSat.push_back(vNStripvNSat);
  }
}

void SaturationCorrection::FillHisto(int charge_sim, int charge_reco, int label, int nstrip, int nsat) {
  for (int counterlayer = 1; counterlayer < 15; counterlayer++) {
    for (int counterstrip = 3; counterstrip < 7; counterstrip++) {
      for (int countersat = 1; countersat < 3; countersat++) {
        if (counterlayer == label && counterstrip == nstrip && countersat == nsat) {
          vModulGeomvNStripvNSat[counterlayer - 1][counterstrip - 3][countersat - 1]->Fill(charge_sim, charge_reco);
        }
      }
    }
  }
}

void SaturationCorrection::SetFileAndTreeName(string file_name, string tree_name) {
  file_ = new TFile(file_name.c_str(), "RECREATE");
  tree_ = new TTree(tree_name.c_str(), tree_name.c_str());
}

void SaturationCorrection::SetHisto(TH2F& histo) {
  histo_ = &histo;
  histo_->Write();
}

void SaturationCorrection::FillHistoRec() {
  float lowedgexu = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX());
  //float lowedgexd = histo_->GetXaxis()->GetBinCenter(0); UNUSED
  float lowedgeyu = histo_->GetYaxis()->GetBinCenter(histo_->GetNbinsY());
  int integral = histo_->Integral();
  int countInt = 0;
  float ratioIntegral = 0.;
  TH2F* histo_clone = (TH2F*)histo_->Clone();
  histo_clone->Reset();
  int x_histo = 0;
  int x_low = 0;
  float integer = 0.9;
  while (ratioIntegral < integer) {
    x_histo++;
    for (int y_histo = 0; y_histo < histo_->GetNbinsY() + 2; y_histo++) {
      countInt += histo_->GetBinContent(x_histo, y_histo);
    }
    ratioIntegral = ((float)countInt / (float)integral);
    if (ratioIntegral < 0.1)
      x_low++;
  }
  x_low = 0;
  for (int i = x_low; i < x_histo; i++) {
    float LowEdgeX = histo_->GetXaxis()->GetBinLowEdge(i);
    //float WidthX = histo_->GetXaxis()->GetBinWidth(i); UNUSED
    for (int j = 0; j < histo_clone->GetNbinsY() + 2; j++) {
      float LowEdgeY = histo_->GetXaxis()->GetBinLowEdge(j);
      //float WidthY = histo_->GetXaxis()->GetBinWidth(j); UNUSED
      if (LowEdgeX > LowEdgeY && LowEdgeX < lowedgexu &&
          LowEdgeY < lowedgeyu)  // je remplis que la partie sous la diagonale Erec=Eloss
      {
        histo_clone->SetBinContent(i, j, histo_->GetBinContent(i, j));
        histo_clone->SetBinError(i, j, histo_->GetBinError(i, j));
      }
    }
  }
  histo_ = histo_clone;
  string first_title = histo_clone->GetTitle();
  first_title += "_NoDiag";
  histo_->SetTitle(first_title.c_str());
  histo_->SetName(first_title.c_str());
  histo_->Write();
}

void SaturationCorrection::FillProfile() {
  int firstdiv = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX()) / 2;
  int divider = histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX()) / firstdiv;
  historec_ = new TH1F(histo_->GetTitle(),
                       histo_->GetTitle(),
                       histo_->GetNbinsX() / divider,
                       0,
                       histo_->GetXaxis()->GetBinCenter(histo_->GetNbinsX()));
  TH2F* histoclone = (TH2F*)histo_->Clone();
  histoclone->Reset();
  int binrec = 0;

  for (int x = 1; x <= histo_->GetNbinsX(); x++) {
    for (int y = 1; y <= histo_->GetNbinsY(); y++) {
      if (histo_->GetBinContent(x, y) > 0)
        histoclone->SetBinContent(x, y, histo_->GetBinContent(x, y));
    }
    if (x % divider == 0) {
      binrec++;
      float mpv = 1.;
      float mpverror = 0.;
      TH1F* projY = (TH1F*)histoclone->ProjectionY();
      if (projY->GetEntries() >= 5) {
        int fitStatus = projY->Fit("gaus", "QL", "");
        TFitResultPtr func_fit = projY->Fit("gaus", "QLS", "");
        if (fitStatus == 0) {
          mpv = func_fit->Parameter(1);
          mpverror = func_fit->Error(1);
        }
        if (mpv < projY->GetMean() + 3 * projY->GetStdDev() && mpv > projY->GetMean() - 3 * projY->GetStdDev()) {
          historec_->SetBinContent(binrec, mpv);
          historec_->SetBinError(binrec, mpverror);
        }
        projY->Write();
      }
      histoclone->Reset();
      delete projY;
    }
  }
  delete histoclone;
  delete histo_;
}

void SaturationCorrection::SetTree(TTree& tree) { inputtree_ = &tree; }

void SaturationCorrection::SetBranch() {
  tree_->Branch("p0", &p0_, "p0/F");
  tree_->Branch("p1", &p1_, "p1/F");
  tree_->Branch("p0err", &p0err_, "p0err/F");
  tree_->Branch("p1err", &p1err_, "p1err/F");
  tree_->Branch("chi2", &chi2_, "chi2/F");
  tree_->Branch("label", &label_, "label/I");
  tree_->Branch("nstrip", &nstrip_, "nstrip/I");
  tree_->Branch("nsat", &nsat_, "nsat/I");
}

void SaturationCorrection::Write(int label, int nstrip, int nsat) {
  if (historec_->GetEntries() >= 3 &&
      (nsat >
       0))  //je ne recupere les resultats du fit que pour des profils avec plus de dix entrees et ou il y a eu saturation
  {
    TFitResultPtr FitRes_ = historec_->Fit("pol1", "SQ", "");
    int fitStat = historec_->Fit("pol1", "Q", "");
    if (fitStat == 0) {
      p0_ = FitRes_->Parameter(0);
      p1_ = FitRes_->Parameter(1);
      p0err_ = FitRes_->Error(0);
      p1err_ = FitRes_->Error(1);
      float chicalc = 0.;
      int nentries = 0;
      for (int bin = 0; bin < historec_->GetNbinsX(); bin++) {
        if (historec_->GetBinContent(bin) > 0) {
          nentries++;
          float valueCalc = (historec_->GetXaxis()->GetBinCenter(bin) * p1_ + p0_);
          chicalc += pow(valueCalc - historec_->GetBinContent(bin), 2);
        }
      }
      chicalc = (sqrt(chicalc) / nentries);
      chi2_ = chicalc;
    }
  } else  //si il n'y a pas eu de saturation ou s'il n'y a pas assez d'entrees dans le profil, je parametrise tel qu'il n'y ait pas de changement
  {
    p0_ = 0.;
    p1_ = 1.;
    p0err_ = 0.;
    p1err_ = 0.;
    chi2_ = 0.;
  }
  label_ = label;
  nstrip_ = nstrip;
  nsat_ = nsat;
  tree_->Fill();
  string title = historec_->GetTitle();
  title += "_rec";
  historec_->SetTitle(title.c_str());
  historec_->SetName(title.c_str());
  historec_->Write();
  delete historec_;
}

void SaturationCorrection::WriteFile() {
  file_->Write();
  file_->Close();
}

void SaturationCorrection::ReadParameters() {
  inputtree_->SetBranchAddress("p0", &p0_);
  inputtree_->SetBranchAddress("p1", &p1_);
  inputtree_->SetBranchAddress("chi2", &chi2_);
  inputtree_->SetBranchAddress("label", &label_);
  inputtree_->SetBranchAddress("nstrip", &nstrip_);
  inputtree_->SetBranchAddress("nsat", &nsat_);
  for (int i = 0; i < inputtree_->GetEntries(); i++) {
    inputtree_->GetEntry(i);
    for (int label = 1; label < 15; label++) {
      for (int nstrip = 3; nstrip < 7; nstrip++) {
        for (int nsat = 1; nsat < 3; nsat++) {
          if (label_ == label && nstrip_ == nstrip && nsat_ == nsat) {
            p0[label - 1][nstrip - 3][nsat - 1] = p0_;
            p1[label - 1][nstrip - 3][nsat - 1] = p1_;
            chi2[label - 1][nstrip - 3][nsat - 1] = chi2_;
          }
        }
      }
    }
  }
}

void SaturationCorrection::WriteParameters() {
  for (int countlayer = 1; countlayer < 15; countlayer++) {
    for (int countnstrip = 3; countnstrip < 7; countnstrip++) {
      for (int countnsat = 1; countnsat < 3; countnsat++) {
        SetHisto(*vModulGeomvNStripvNSat[countlayer - 1][countnstrip - 3][countnsat - 1]);
        FillHistoRec();
        FillProfile();
        Write(countlayer, countnstrip, countnsat);
      }
    }
  }
  file_->Write();
  file_->Close();
}

float SaturationCorrection::ChargeCorrected(float charge, int label, int nstrip, int nsat) {
  bool test = true;
  if (nstrip >= 6)
    nstrip = 6;  //inclusif pour le nombre de strips du cluster
  if (nsat >= 2)
    nsat = 2;  //inclusif pour le nombre de strips saturees

  float p0Calc = p0[label - 1][nstrip - 3][nsat - 1];
  float p1Calc = p1[label - 1][nstrip - 3][nsat - 1];
  float res = charge;
  if (test)
    res = (charge - p0Calc) / p1Calc;
  return res;
}

#endif
