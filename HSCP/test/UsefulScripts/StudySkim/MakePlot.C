

//TEST

#include <vector>

#include "TROOT.h"
#include "TFile.h"
#include "TDirectory.h"
#include "TChain.h"
#include "TObject.h"
#include "TCanvas.h"
#include "TMath.h"
#include "TLegend.h"
#include "TGraph.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TTree.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TGraphAsymmErrors.h"
#include "TPaveText.h"
#include "TCutG.h"
#include "TProfile.h"
#include "tdrstyle.C"

TGraphErrors* GetGraph(string SampleName, string PT, string NH, string DEDX) {
  string listPT[] = {"20.0", "40.0", "60.0", "100"};
  string listNH[] = {"3", "5", "7"};
  string listDEDX[] = {"2.7", "3.0", "3.3"};

  std::vector<string> fileList;
  if (PT == "") {
    for (size_t i = 0; i < sizeof(listPT) / sizeof(string); i++) {
      fileList.push_back(string("Results/SkimEff_") + SampleName + "_pT" + listPT[i] + "_nH" + NH + "_dEdx" + DEDX +
                         ".txt");
    }
  } else if (NH == "") {
    for (size_t i = 0; i < sizeof(listNH) / sizeof(string); i++) {
      fileList.push_back(string("Results/SkimEff_") + SampleName + "_pT" + PT + "_nH" + listNH[i] + "_dEdx" + DEDX +
                         ".txt");
    }
  } else if (DEDX == "") {
    for (size_t i = 0; i < sizeof(listDEDX) / sizeof(string); i++) {
      fileList.push_back(string("Results/SkimEff_") + SampleName + "_pT" + PT + "_nH" + NH + "_dEdx" + listDEDX[i] +
                         ".txt");
    }
  } else {
    fileList.push_back(string("Results/SkimEff_") + SampleName + "_pT" + PT + "_nH" + NH + "_dEdx" + DEDX + ".txt");
  }

  TGraphErrors* graph = new TGraphErrors(100);

  size_t N = 0;
  for (size_t f = 0; f < fileList.size(); f++) {
    FILE* pFile = fopen(fileList[f].c_str(), "r");
    if (!pFile) {
      printf("Skip file %50s (missing)\n", fileList[f].c_str());
      continue;
    }

    char name[1024];
    double Skimmed;
    double Total;
    double Eff;
    double Size;
    fscanf(pFile, "%s : Skimmed/Total = %lf/%lf = %lf%%  EventSize=%lfKB", name, &Skimmed, &Total, &Eff, &Size);
    Eff /= 100.0;

    double pT, nh, dedx;
    sscanf(name + 8 + SampleName.length(), "_pT%lf_nH%lf_dEdx%lf", &pT, &nh, &dedx);
    printf("%s --> %s -> %f %f %f --> %f+-%f\n",
           fileList[f].c_str(),
           name,
           pT,
           nh,
           dedx,
           Eff,
           sqrt(Eff * (1 - Eff) / Total));

    double val = pT;
    if (NH == "")
      val = nh;
    if (DEDX == "")
      val = dedx;

    graph->SetPoint(N, val, Eff);
    graph->SetPointError(N, 0, sqrt(Eff * (1 - Eff) / Total));
    N++;
    fclose(pFile);
  }
  graph->Set(N);
  return graph;
}

void MakePlot() {
  setTDRStyle();
  gStyle->SetPadTopMargin(0.05);
  gStyle->SetPadBottomMargin(0.13);
  gStyle->SetPadRightMargin(0.30);
  gStyle->SetPadLeftMargin(0.13);
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.35);
  gStyle->SetPalette(1);
  gStyle->SetNdivisions(505);
  gStyle->SetPaintTextFormat("4.2f");
  gStyle->SetOptStat("");

  int COLOR[] = {1, 12, 11, 7, 8, 3, 9, 4};
  string SAMPLES[] = {"MU", "MET", "ELEC", "Gluino300", "Gluino1000", "Gluino1500", "GMStau126", "GMStau494"};
  string LEGENDS[] = {
      "Mu PD", "MET PD", "EGamma PD", "Gluino (300)", "Gluino (1000)", "Gluino (1500)", "GMStau (126)", "GMStau (494)"};
  size_t NSamples = sizeof(SAMPLES) / sizeof(string);

  for (unsigned int SET = 1; SET <= 3; SET++) {
    string CUTS[] = {"40.0", "3", ""};
    if (SET == 1) {
      CUTS[0] = "40.0";
      CUTS[1] = "3";
      CUTS[2] = "";
    }
    if (SET == 2) {
      CUTS[0] = "40.0";
      CUTS[1] = "";
      CUTS[2] = "3.0";
    }
    if (SET == 3) {
      CUTS[0] = "";
      CUTS[1] = "3";
      CUTS[2] = "3.0";
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
    TH1D* frame = NULL;
    if (CUTS[0] == "")
      frame = new TH1D("frame", "frame;p_{T} [GeV];Efficiency", 1, 15, 105);
    if (CUTS[1] == "")
      frame = new TH1D("frame", "frame;nHits;Efficiency", 1, 2.5, 7.5);
    if (CUTS[2] == "")
      frame = new TH1D("frame", "frame;dEdx [MeV/cm];Efficiency", 1, 2.6, 3.4);
    frame->SetTitle("");
    frame->SetStats(kFALSE);
    frame->GetXaxis()->SetTitleOffset(1.50);
    frame->GetYaxis()->SetTitleOffset(1.50);
    frame->SetMaximum(1.05);
    frame->SetMinimum(0.00);
    frame->GetXaxis()->SetLabelFont(43);  //give the font size in pixel (instead of fraction)
    frame->GetXaxis()->SetLabelSize(22);  //font size
    frame->GetXaxis()->SetTitleFont(43);  //give the font size in pixel (instead of fraction)
    frame->GetXaxis()->SetTitleSize(22);  //font size
    frame->GetXaxis()->SetNdivisions(510);
    frame->GetYaxis()->SetLabelFont(43);  //give the font size in pixel (instead of fraction)
    frame->GetYaxis()->SetLabelSize(22);  //font size
    frame->GetYaxis()->SetTitleFont(43);  //give the font size in pixel (instead of fraction)
    frame->GetYaxis()->SetTitleSize(22);  //font size
    frame->GetYaxis()->SetNdivisions(510);
    frame->Draw("AXIS");

    TLegend* leg = new TLegend(0.70, 0.98, 0.98, 0.98 - (NSamples + 3) * 0.075);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);
    leg->SetTextSize(23);
    leg->SetHeader(NULL);  //"pMSSM LHC subspace");
    if (CUTS[0] != "")
      leg->AddEntry("", (string("p_{T}#geq") + CUTS[0]).c_str(), "");
    if (CUTS[1] != "")
      leg->AddEntry("", (string("nHits#geq") + CUTS[1]).c_str(), "");
    if (CUTS[2] != "")
      leg->AddEntry("", (string("dEdx#geq") + CUTS[2]).c_str(), "");
    leg->AddEntry("", "  ", "");

    TGraphErrors** graphs = new TGraphErrors*[NSamples];
    for (size_t s = 0; s < NSamples; s++) {
      graphs[s] = GetGraph(SAMPLES[s], CUTS[0], CUTS[1], CUTS[2]);
      graphs[s]->SetLineWidth(3);
      graphs[s]->SetLineColor(COLOR[s]);
      graphs[s]->SetMarkerColor(COLOR[s]);
      graphs[s]->Draw("same");
      leg->AddEntry(graphs[s], LEGENDS[s].c_str(), "L");
    }
    leg->Draw();

    char saveName[256];
    sprintf(saveName, "Results/PlotSet%i.png", SET);
    c1->SaveAs(saveName);
  }
}
