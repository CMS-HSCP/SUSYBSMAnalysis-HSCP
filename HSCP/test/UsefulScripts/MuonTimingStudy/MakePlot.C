// Original Author:  Loic Quertenmont

#include <vector>
#include <algorithm>

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
#include "TProfile.h"
#include "TCutG.h"
#include "TGaxis.h"

#include "../../AnalysisCode/tdrstyle.C"

//#include "DataFormats/MuonReco/interface/Muon.h"
//#include "DataFormats/MuonReco/interface/MuonFwd.h"
#include "DataFormats/MuonDetId/interface/CSCDetId.h"
#include "DataFormats/MuonDetId/interface/DTWireId.h"
//#include "DataFormats/CSCRecHit/interface/CSCSegmentCollection.h"
//#include "DataFormats/MuonDetId/interface/CSCIndexer.h"
//#include "DataFormats/DTRecHit/interface/DTRecSegment4DCollection.h"
//#include "DataFormats/DTRecHit/interface/DTRecHitCollection.h"

TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy = false);  //defined below

std::map<string, unsigned int> nameToDetId;
std::map<string, TGraphErrors*> graphMap;

string getStationName(unsigned int detId) {
  char stationName[255];
  DetId geomDetId(detId);
  if (geomDetId.subdetId() == 1) {
    DTChamberId dtId(detId);
    sprintf(stationName, "DT%i%i", dtId.wheel(), dtId.station());
  } else if (geomDetId.subdetId() == 2) {
    CSCDetId cscId(detId);
    CSCDetId cscchamberId = cscId.chamberId();
    sprintf(
        stationName, "ME%c%i%i", cscchamberId.zendcap() > 0 ? '+' : '-', cscchamberId.station(), cscchamberId.ring());
  }
  return stationName;
}

string getChamberName(unsigned int detId) {
  char chamberName[255];
  DetId geomDetId(detId);
  if (geomDetId.subdetId() == 1) {
    DTChamberId dtId(detId);
    sprintf(chamberName, "%s_%02i", getStationName(detId).c_str(), dtId.sector());
  } else if (geomDetId.subdetId() == 2) {
    CSCDetId cscId(detId);
    CSCDetId cscchamberId = cscId.chamberId();
    sprintf(chamberName, "%s_%02i", getStationName(detId).c_str(), cscchamberId.chamber());
  }
  return chamberName;
}

void GetMeanAndRMS(TFile* InputFile, string path, double& mean, double& rms) {
  TH1D* histo = (TH1D*)GetObjectFromPath(InputFile, path);
  if (!histo) {
    mean = -999;
    rms = -999;
  } else {
    mean = histo->GetMean();
    rms = histo->GetRMS();

    if (histo->GetEntries() > 10) {
      TF1* fit = new TF1("MyFit", "gaus", mean - 2 * rms, mean + 2 * rms);
      fit->SetParameters(histo->GetEntries(), histo->GetMean(), histo->GetRMS());
      histo->Fit(fit, "0QR");
      mean = fit->GetParameter(1);
      rms = fit->GetParError(1);
      delete fit;
    }
    delete histo;
  }
}

void makeFigure(string outfile, TH1* frame, string chambersList) {
  TCanvas* c1 = new TCanvas("c1", "c1,", 1200, 600);
  c1->SetGrid(1, 1);
  frame->Draw("AXIS");

  std::vector<string> chambersToDraw;
  char* pch = strtok((char*)chambersList.c_str(), ";");
  while (pch != NULL) {
    chambersToDraw.push_back(pch);
    pch = strtok(NULL, ";");
  }

  TLegend* LEG = new TLegend(0.15, 0.20, 0.90, 0.40);
  if (chambersToDraw.size() <= 10)
    LEG->SetNColumns(5);
  else if (chambersToDraw.size() <= 18)
    LEG->SetNColumns(9);
  else
    LEG->SetNColumns(9);

  LEG->SetFillColor(0);
  LEG->SetFillStyle(0);
  LEG->SetBorderSize(0);

  for (unsigned int c = 0; c < chambersToDraw.size(); c++) {
    TGraphErrors* graph = graphMap[chambersToDraw[c]];
    double* xpoints = graph->GetX();
    for (unsigned int x = 0; x < graph->GetN(); x++) {
      xpoints[x] = floor(xpoints[x]) + double(1 + c % chambersToDraw.size()) / (chambersToDraw.size() + 1);
    }
    if (outfile.find("Station") != string::npos) {
      double* yErr = graph->GetEY();
      double* yVal = graph->GetY();
      for (unsigned int x = 0; x < graph->GetN(); x++) {
        const char* runString = frame->GetXaxis()->GetBinLabel(frame->FindBin(xpoints[x]));
        if (yErr[x] > 10.0 || std::fabs(yVal[x]) > 25.0)
          printf("%s: %lf +- %lf\n", runString, yVal[x], yErr[x]);
      }
    }
    graph->SetMarkerStyle(20);
    graph->SetMarkerColor(gStyle->GetColorPalette(int(c * (255.0 / chambersToDraw.size()))));
    graph->SetLineColor(graph->GetMarkerColor());
    graph->Draw("P same");
    LEG->AddEntry(graph, graph->GetName(), "P");
  }

  TLine line(frame->GetXaxis()->GetXmin(), 0.0, frame->GetXaxis()->GetXmax(), 0.0);
  line.SetLineColor(1);
  line.SetLineWidth(2);
  line.SetLineStyle(2);
  line.Draw();

  LEG->Draw();
  //DrawPreliminary("", SQRTS, IntegratedLuminosityFromE(SQRTS));
  c1->SaveAs(outfile.c_str());
  delete c1;
  delete LEG;
}

void MakePlot() {
  setTDRStyle();
  gStyle->SetPadTopMargin(0.06);
  gStyle->SetPadBottomMargin(0.15);
  gStyle->SetPadRightMargin(0.03);
  gStyle->SetPadLeftMargin(0.07);
  gStyle->SetTitleSize(0.04, "XYZ");
  gStyle->SetTitleXOffset(1.1);
  gStyle->SetTitleYOffset(1.35);
  gStyle->SetPalette(55);  //Rainbow --> deifne the color of all plots
  gStyle->SetNdivisions(505);

  system("mkdir -p pictures");

  TFile* InputFile = new TFile("Histos.root");

  std::vector<string> runList;
  string chambers[] = {
      "ME+11",    "ME+12",    "ME+13",    "ME+14",    "ME+21",    "ME+22",    "ME+31",    "ME+32",    "ME+41",
      "ME+42",    "ME+11_01", "ME+11_02", "ME+11_03", "ME+11_04", "ME+11_05", "ME+11_06", "ME+11_07", "ME+11_08",
      "ME+11_09", "ME+11_10", "ME+11_11", "ME+11_12", "ME+11_13", "ME+11_14", "ME+11_15", "ME+11_16", "ME+11_17",
      "ME+11_18", "ME+11_19", "ME+11_20", "ME+11_21", "ME+11_22", "ME+11_23", "ME+11_24", "ME+11_25", "ME+11_26",
      "ME+11_27", "ME+11_28", "ME+11_29", "ME+11_30", "ME+11_31", "ME+11_32", "ME+11_33", "ME+11_34", "ME+11_35",
      "ME+11_36", "ME+12_01", "ME+12_02", "ME+12_03", "ME+12_04", "ME+12_05", "ME+12_06", "ME+12_07", "ME+12_08",
      "ME+12_09", "ME+12_10", "ME+12_11", "ME+12_12", "ME+12_13", "ME+12_14", "ME+12_15", "ME+12_16", "ME+12_17",
      "ME+12_18", "ME+12_19", "ME+12_20", "ME+12_21", "ME+12_22", "ME+12_23", "ME+12_24", "ME+12_25", "ME+12_26",
      "ME+12_27", "ME+12_28", "ME+12_29", "ME+12_30", "ME+12_31", "ME+12_32", "ME+12_33", "ME+12_34", "ME+12_35",
      "ME+12_36", "ME+13_01", "ME+13_02", "ME+13_03", "ME+13_04", "ME+13_05", "ME+13_06", "ME+13_07", "ME+13_08",
      "ME+13_09", "ME+13_10", "ME+13_11", "ME+13_12", "ME+13_13", "ME+13_14", "ME+13_15", "ME+13_16", "ME+13_17",
      "ME+13_18", "ME+13_19", "ME+13_20", "ME+13_21", "ME+13_22", "ME+13_23", "ME+13_24", "ME+13_25", "ME+13_26",
      "ME+13_27", "ME+13_28", "ME+13_29", "ME+13_30", "ME+13_31", "ME+13_32", "ME+13_33", "ME+13_34", "ME+13_35",
      "ME+13_36", "ME+14_01", "ME+14_02", "ME+14_03", "ME+14_04", "ME+14_05", "ME+14_06", "ME+14_07", "ME+14_08",
      "ME+14_09", "ME+14_10", "ME+14_11", "ME+14_12", "ME+14_13", "ME+14_14", "ME+14_15", "ME+14_16", "ME+14_17",
      "ME+14_18", "ME+14_19", "ME+14_20", "ME+14_21", "ME+14_22", "ME+14_23", "ME+14_24", "ME+14_25", "ME+14_26",
      "ME+14_27", "ME+14_28", "ME+14_29", "ME+14_30", "ME+14_31", "ME+14_32", "ME+14_33", "ME+14_34", "ME+14_35",
      "ME+14_36", "ME+21_01", "ME+21_02", "ME+21_03", "ME+21_04", "ME+21_05", "ME+21_06", "ME+21_07", "ME+21_08",
      "ME+21_09", "ME+21_10", "ME+21_11", "ME+21_12", "ME+21_13", "ME+21_14", "ME+21_15", "ME+21_16", "ME+21_17",
      "ME+21_18", "ME+22_01", "ME+22_02", "ME+22_03", "ME+22_04", "ME+22_05", "ME+22_06", "ME+22_07", "ME+22_08",
      "ME+22_09", "ME+22_10", "ME+22_11", "ME+22_12", "ME+22_13", "ME+22_14", "ME+22_15", "ME+22_16", "ME+22_17",
      "ME+22_18", "ME+22_19", "ME+22_20", "ME+22_21", "ME+22_22", "ME+22_23", "ME+22_24", "ME+22_25", "ME+22_26",
      "ME+22_27", "ME+22_28", "ME+22_29", "ME+22_30", "ME+22_31", "ME+22_32", "ME+22_33", "ME+22_34", "ME+22_35",
      "ME+22_36", "ME+31_01", "ME+31_02", "ME+31_03", "ME+31_04", "ME+31_05", "ME+31_06", "ME+31_07", "ME+31_08",
      "ME+31_09", "ME+31_10", "ME+31_11", "ME+31_12", "ME+31_13", "ME+31_14", "ME+31_15", "ME+31_16", "ME+31_17",
      "ME+31_18", "ME+32_01", "ME+32_02", "ME+32_03", "ME+32_04", "ME+32_05", "ME+32_06", "ME+32_07", "ME+32_08",
      "ME+32_09", "ME+32_10", "ME+32_11", "ME+32_12", "ME+32_13", "ME+32_14", "ME+32_15", "ME+32_16", "ME+32_17",
      "ME+32_18", "ME+32_19", "ME+32_20", "ME+32_21", "ME+32_22", "ME+32_23", "ME+32_24", "ME+32_25", "ME+32_26",
      "ME+32_27", "ME+32_28", "ME+32_29", "ME+32_30", "ME+32_31", "ME+32_32", "ME+32_33", "ME+32_34", "ME+32_35",
      "ME+32_36", "ME+41_01", "ME+41_02", "ME+41_03", "ME+41_04", "ME+41_05", "ME+41_06", "ME+41_07", "ME+41_08",
      "ME+41_09", "ME+41_10", "ME+41_11", "ME+41_12", "ME+41_13", "ME+41_14", "ME+41_15", "ME+41_16", "ME+41_17",
      "ME+41_18", "ME+42_01", "ME+42_02", "ME+42_03", "ME+42_04", "ME+42_05", "ME+42_06", "ME+42_07", "ME+42_08",
      "ME+42_09", "ME+42_10", "ME+42_11", "ME+42_12", "ME+42_13", "ME+42_14", "ME+42_15", "ME+42_16", "ME+42_17",
      "ME+42_18", "ME+42_19", "ME+42_20", "ME+42_21", "ME+42_22", "ME+42_23", "ME+42_24", "ME+42_25", "ME+42_26",
      "ME+42_27", "ME+42_28", "ME+42_29", "ME+42_30", "ME+42_31", "ME+42_32", "ME+42_33", "ME+42_34", "ME+42_35",
      "ME+42_36",

      "ME-11",    "ME-12",    "ME-13",    "ME-14",    "ME-21",    "ME-22",    "ME-31",    "ME-32",    "ME-41",
      "ME-42",    "ME-11_01", "ME-11_02", "ME-11_03", "ME-11_04", "ME-11_05", "ME-11_06", "ME-11_07", "ME-11_08",
      "ME-11_09", "ME-11_10", "ME-11_11", "ME-11_12", "ME-11_13", "ME-11_14", "ME-11_15", "ME-11_16", "ME-11_17",
      "ME-11_18", "ME-11_19", "ME-11_20", "ME-11_21", "ME-11_22", "ME-11_23", "ME-11_24", "ME-11_25", "ME-11_26",
      "ME-11_27", "ME-11_28", "ME-11_29", "ME-11_30", "ME-11_31", "ME-11_32", "ME-11_33", "ME-11_34", "ME-11_35",
      "ME-11_36", "ME-12_01", "ME-12_02", "ME-12_03", "ME-12_04", "ME-12_05", "ME-12_06", "ME-12_07", "ME-12_08",
      "ME-12_09", "ME-12_10", "ME-12_11", "ME-12_12", "ME-12_13", "ME-12_14", "ME-12_15", "ME-12_16", "ME-12_17",
      "ME-12_18", "ME-12_19", "ME-12_20", "ME-12_21", "ME-12_22", "ME-12_23", "ME-12_24", "ME-12_25", "ME-12_26",
      "ME-12_27", "ME-12_28", "ME-12_29", "ME-12_30", "ME-12_31", "ME-12_32", "ME-12_33", "ME-12_34", "ME-12_35",
      "ME-12_36", "ME-13_01", "ME-13_02", "ME-13_03", "ME-13_04", "ME-13_05", "ME-13_06", "ME-13_07", "ME-13_08",
      "ME-13_09", "ME-13_10", "ME-13_11", "ME-13_12", "ME-13_13", "ME-13_14", "ME-13_15", "ME-13_16", "ME-13_17",
      "ME-13_18", "ME-13_19", "ME-13_20", "ME-13_21", "ME-13_22", "ME-13_23", "ME-13_24", "ME-13_25", "ME-13_26",
      "ME-13_27", "ME-13_28", "ME-13_29", "ME-13_30", "ME-13_31", "ME-13_32", "ME-13_33", "ME-13_34", "ME-13_35",
      "ME-13_36", "ME-14_01", "ME-14_02", "ME-14_03", "ME-14_04", "ME-14_05", "ME-14_06", "ME-14_07", "ME-14_08",
      "ME-14_09", "ME-14_10", "ME-14_11", "ME-14_12", "ME-14_13", "ME-14_14", "ME-14_15", "ME-14_16", "ME-14_17",
      "ME-14_18", "ME-14_19", "ME-14_20", "ME-14_21", "ME-14_22", "ME-14_23", "ME-14_24", "ME-14_25", "ME-14_26",
      "ME-14_27", "ME-14_28", "ME-14_29", "ME-14_30", "ME-14_31", "ME-14_32", "ME-14_33", "ME-14_34", "ME-14_35",
      "ME-14_36", "ME-21_01", "ME-21_02", "ME-21_03", "ME-21_04", "ME-21_05", "ME-21_06", "ME-21_07", "ME-21_08",
      "ME-21_09", "ME-21_10", "ME-21_11", "ME-21_12", "ME-21_13", "ME-21_14", "ME-21_15", "ME-21_16", "ME-21_17",
      "ME-21_18", "ME-22_01", "ME-22_02", "ME-22_03", "ME-22_04", "ME-22_05", "ME-22_06", "ME-22_07", "ME-22_08",
      "ME-22_09", "ME-22_10", "ME-22_11", "ME-22_12", "ME-22_13", "ME-22_14", "ME-22_15", "ME-22_16", "ME-22_17",
      "ME-22_18", "ME-22_19", "ME-22_20", "ME-22_21", "ME-22_22", "ME-22_23", "ME-22_24", "ME-22_25", "ME-22_26",
      "ME-22_27", "ME-22_28", "ME-22_29", "ME-22_30", "ME-22_31", "ME-22_32", "ME-22_33", "ME-22_34", "ME-22_35",
      "ME-22_36", "ME-31_01", "ME-31_02", "ME-31_03", "ME-31_04", "ME-31_05", "ME-31_06", "ME-31_07", "ME-31_08",
      "ME-31_09", "ME-31_10", "ME-31_11", "ME-31_12", "ME-31_13", "ME-31_14", "ME-31_15", "ME-31_16", "ME-31_17",
      "ME-31_18", "ME-32_01", "ME-32_02", "ME-32_03", "ME-32_04", "ME-32_05", "ME-32_06", "ME-32_07", "ME-32_08",
      "ME-32_09", "ME-32_10", "ME-32_11", "ME-32_12", "ME-32_13", "ME-32_14", "ME-32_15", "ME-32_16", "ME-32_17",
      "ME-32_18", "ME-32_19", "ME-32_20", "ME-32_21", "ME-32_22", "ME-32_23", "ME-32_24", "ME-32_25", "ME-32_26",
      "ME-32_27", "ME-32_28", "ME-32_29", "ME-32_30", "ME-32_31", "ME-32_32", "ME-32_33", "ME-32_34", "ME-32_35",
      "ME-32_36", "ME-41_01", "ME-41_02", "ME-41_03", "ME-41_04", "ME-41_05", "ME-41_06", "ME-41_07", "ME-41_08",
      "ME-41_09", "ME-41_10", "ME-41_11", "ME-41_12", "ME-41_13", "ME-41_14", "ME-41_15", "ME-41_16", "ME-41_17",
      "ME-41_18", "ME-42_01", "ME-42_02", "ME-42_03", "ME-42_04", "ME-42_05", "ME-42_06", "ME-42_07", "ME-42_08",
      "ME-42_09", "ME-42_10", "ME-42_11", "ME-42_12", "ME-42_13", "ME-42_14", "ME-42_15", "ME-42_16", "ME-42_17",
      "ME-42_18", "ME-42_19", "ME-42_20", "ME-42_21", "ME-42_22", "ME-42_23", "ME-42_24", "ME-42_25", "ME-42_26",
      "ME-42_27", "ME-42_28", "ME-42_29", "ME-42_30", "ME-42_31", "ME-42_32", "ME-42_33", "ME-42_34", "ME-42_35",
      "ME-42_36",

      "DT-21",    "DT-22",    "DT-23",    "DT-24",    "DT-11",    "DT-12",    "DT-13",    "DT-14",    "DT01",
      "DT02",     "DT03",     "DT04",     "DT11",     "DT12",     "DT13",     "DT14",     "DT21",     "DT22",
      "DT23",     "DT24",     "DT-21_01", "DT-21_02", "DT-21_03", "DT-21_04", "DT-21_05", "DT-21_06", "DT-21_07",
      "DT-21_08", "DT-21_09", "DT-21_10", "DT-21_11", "DT-21_12", "DT-22_01", "DT-22_02", "DT-22_03", "DT-22_04",
      "DT-22_05", "DT-22_06", "DT-22_07", "DT-22_08", "DT-22_09", "DT-22_10", "DT-22_11", "DT-22_12", "DT-23_01",
      "DT-23_02", "DT-23_03", "DT-23_04", "DT-23_05", "DT-23_06", "DT-23_07", "DT-23_08", "DT-23_09", "DT-23_10",
      "DT-23_11", "DT-23_12", "DT-24_01", "DT-24_02", "DT-24_03", "DT-24_04", "DT-24_05", "DT-24_06", "DT-24_07",
      "DT-24_08", "DT-24_09", "DT-24_10", "DT-24_11", "DT-24_12", "DT-24_13", "DT-24_14", "DT-11_01", "DT-11_02",
      "DT-11_03", "DT-11_04", "DT-11_05", "DT-11_06", "DT-11_07", "DT-11_08", "DT-11_09", "DT-11_10", "DT-11_11",
      "DT-11_12", "DT-12_01", "DT-12_02", "DT-12_03", "DT-12_04", "DT-12_05", "DT-12_06", "DT-12_07", "DT-12_08",
      "DT-12_09", "DT-12_10", "DT-12_11", "DT-12_12", "DT-13_01", "DT-13_02", "DT-13_03", "DT-13_04", "DT-13_05",
      "DT-13_06", "DT-13_07", "DT-13_08", "DT-13_09", "DT-13_10", "DT-13_11", "DT-13_12", "DT-14_01", "DT-14_02",
      "DT-14_03", "DT-14_04", "DT-14_05", "DT-14_06", "DT-14_07", "DT-14_08", "DT-14_09", "DT-14_10", "DT-14_11",
      "DT-14_12", "DT-14_13", "DT-14_14", "DT01_01",  "DT01_02",  "DT01_03",  "DT01_04",  "DT01_05",  "DT01_06",
      "DT01_07",  "DT01_08",  "DT01_09",  "DT01_10",  "DT01_11",  "DT01_12",  "DT02_01",  "DT02_02",  "DT02_03",
      "DT02_04",  "DT02_05",  "DT02_06",  "DT02_07",  "DT02_08",  "DT02_09",  "DT02_10",  "DT02_11",  "DT02_12",
      "DT03_01",  "DT03_02",  "DT03_03",  "DT03_04",  "DT03_05",  "DT03_06",  "DT03_07",  "DT03_08",  "DT03_09",
      "DT03_10",  "DT03_11",  "DT03_12",  "DT04_01",  "DT04_02",  "DT04_03",  "DT04_04",  "DT04_05",  "DT04_06",
      "DT04_07",  "DT04_08",  "DT04_09",  "DT04_10",  "DT04_11",  "DT04_12",  "DT04_13",  "DT04_14",  "DT11_01",
      "DT11_02",  "DT11_03",  "DT11_04",  "DT11_05",  "DT11_06",  "DT11_07",  "DT11_08",  "DT11_09",  "DT11_10",
      "DT11_11",  "DT11_12",  "DT12_01",  "DT12_02",  "DT12_03",  "DT12_04",  "DT12_05",  "DT12_06",  "DT12_07",
      "DT12_08",  "DT12_09",  "DT12_10",  "DT12_11",  "DT12_12",  "DT13_01",  "DT13_02",  "DT13_03",  "DT13_04",
      "DT13_05",  "DT13_06",  "DT13_07",  "DT13_08",  "DT13_09",  "DT13_10",  "DT13_11",  "DT13_12",  "DT14_01",
      "DT14_02",  "DT14_03",  "DT14_04",  "DT14_05",  "DT14_06",  "DT14_07",  "DT14_08",  "DT14_09",  "DT14_10",
      "DT14_11",  "DT14_12",  "DT14_13",  "DT14_14",  "DT21_01",  "DT21_02",  "DT21_03",  "DT21_04",  "DT21_05",
      "DT21_06",  "DT21_07",  "DT21_08",  "DT21_09",  "DT21_10",  "DT21_11",  "DT21_12",  "DT22_01",  "DT22_02",
      "DT22_03",  "DT22_04",  "DT22_05",  "DT22_06",  "DT22_07",  "DT22_08",  "DT22_09",  "DT22_10",  "DT22_11",
      "DT22_12",  "DT23_01",  "DT23_02",  "DT23_03",  "DT23_04",  "DT23_05",  "DT23_06",  "DT23_07",  "DT23_08",
      "DT23_09",  "DT23_10",  "DT23_11",  "DT23_12",  "DT24_01",  "DT24_02",  "DT24_03",  "DT24_04",  "DT24_05",
      "DT24_06",  "DT24_07",  "DT24_08",  "DT24_09",  "DT24_10",  "DT24_11",  "DT24_12",  "DT24_13",  "DT24_14",
  };

  unsigned int NChambers = sizeof(chambers) / sizeof(string);

  TList* ObjList = InputFile->GetListOfKeys();
  for (int i = 0; i < ObjList->GetSize(); i++) {
    TObject* tmp = GetObjectFromPath(InputFile, ObjList->At(i)->GetName(), false);
    if (tmp->InheritsFrom("TDirectory")) {
      runList.push_back(ObjList->At(i)->GetName());
      printf("Add a new run: %s\n", ObjList->At(i)->GetName());

      //list all histo (done for each run to avoid empty=missing histograms)
      if (runList.size() == 1)
        nameToDetId.clear();
      TDirectory* runDir = (TDirectory*)InputFile->Get(ObjList->At(i)->GetName());
      TList* HistoList = runDir->GetListOfKeys();
      for (int j = 0; j < HistoList->GetSize(); j++) {
        unsigned int detId;
        if (sscanf(HistoList->At(j)->GetName(), "%u", &detId) > 0) {
          DetId geomDetId(detId);
          string detName;
          if (geomDetId.subdetId() == 1 && (detId & 0x003C0000) == 0) {
            detName = getStationName(detId);  //dt stations
          } else if (geomDetId.subdetId() == 2 && (detId & 0x000001F8) == 0) {
            detName = getStationName(detId);  //csc stations
          } else {
            detName = getChamberName(detId);  //dt or csc chambers
          }
          //printf("%i (--> %s --> 0x%x --> 0x%x\n", detId, detName.c_str(), detId, detId&0x003C0000);
          nameToDetId[detName] = detId;
        }
      }
    }
    delete tmp;
  }

  unsigned int N = runList.size();
  graphMap.clear();
  for (unsigned int c = 0; c < NChambers; c++) {
    graphMap[chambers[c]] = new TGraphErrors(N);
    graphMap[chambers[c]]->SetName(chambers[c].c_str());
  }

  TH1D* frame = new TH1D("frame", "frame", N, 0, N);
  frame->SetTitle("");
  frame->SetStats(kFALSE);
  frame->GetXaxis()->SetNdivisions(505);
  frame->GetXaxis()->SetTitle("Run Number");
  frame->GetYaxis()->SetTitle("Segment Time (ns)");
  frame->GetYaxis()->SetTitleOffset(0.95);
  frame->SetMaximum(30.0);
  frame->SetMinimum(-30.0);
  for (unsigned int r = 0; r < N; r++) {
    frame->GetXaxis()->SetBinLabel(r + 1, runList[r].c_str());
  }

  FILE* pFile = fopen("MuonTimeOffset.txt", "w");  //used to dump all the corrections per stations
  for (unsigned int c = 0; c < NChambers; c++) {
    if (nameToDetId[chambers[c]] == 0)
      printf("BUG with %s -->%i\n", chambers[c].c_str(), nameToDetId[chambers[c]]);
  }

  fprintf(pFile, "%10s", "runs");
  for (unsigned int r = 0; r < N; r++) {
    fprintf(pFile, ", %8s", runList[r].c_str());
  }
  fprintf(pFile, "\n");
  fprintf(pFile, "%10s", "chambers");
  for (unsigned int c = 0; c < NChambers; c++) {
    fprintf(pFile, ", %-6i", nameToDetId[chambers[c]]);
  }
  fprintf(pFile, "\n");

  double mean, rms;
  for (unsigned int r = 0; r < N; r++) {
    fprintf(pFile, "run %6s", runList[r].c_str());
    for (unsigned int c = 0; c < NChambers; c++) {
      TGraphErrors* graph = graphMap[chambers[c]];
      char histoName[256];
      sprintf(histoName, "%s/%i", runList[r].c_str(), nameToDetId[chambers[c]]);
      GetMeanAndRMS(InputFile, histoName, mean, rms);
      graph->SetPoint(r, r + 0.5, mean);
      graph->SetPointError(r, 0, rms);
      fprintf(pFile, ", %+08.4f", mean > -999 ? mean : 0.0);
    }
    fprintf(pFile, "\n");
  }
  fclose(pFile);

  makeFigure("pictures/CSC_StationsP.png", frame, "ME+11;ME+12;ME+13;ME+14;ME+21;ME+22;ME+31;ME+32;ME+41;ME+42");
  makeFigure("pictures/CSC_StationsM.png", frame, "ME-11;ME-12;ME-13;ME-14;ME-21;ME-22;ME-31;ME-32;ME-41;ME-42");

  makeFigure(
      "pictures/CSC_ChambersMEP11.png",
      frame,
      "ME+11_01;ME+11_02;ME+11_03;ME+11_04;ME+11_05;ME+11_06;ME+11_07;ME+11_08;ME+11_09;ME+11_10;ME+11_11;ME+11_12;ME+"
      "11_13;ME+11_14;ME+11_15;ME+11_16;ME+11_17;ME+11_18;ME+11_19;ME+11_20;ME+11_21;ME+11_22;ME+11_23;ME+11_24;ME+11_"
      "25;ME+11_26;ME+11_27;ME+11_28;ME+11_29;ME+11_30;ME+11_31;ME+11_32;ME+11_33;ME+11_34;ME+11_35;ME+11_36");
  makeFigure(
      "pictures/CSC_ChambersMEP12.png",
      frame,
      "ME+12_01;ME+12_02;ME+12_03;ME+12_04;ME+12_05;ME+12_06;ME+12_07;ME+12_08;ME+12_09;ME+12_10;ME+12_11;ME+12_12;ME+"
      "12_13;ME+12_14;ME+12_15;ME+12_16;ME+12_17;ME+12_18;ME+12_19;ME+12_20;ME+12_21;ME+12_22;ME+12_23;ME+12_24;ME+12_"
      "25;ME+12_26;ME+12_27;ME+12_28;ME+12_29;ME+12_30;ME+12_31;ME+12_32;ME+12_33;ME+12_34;ME+12_35;ME+12_36");
  makeFigure(
      "pictures/CSC_ChambersMEP13.png",
      frame,
      "ME+13_01;ME+13_02;ME+13_03;ME+13_04;ME+13_05;ME+13_06;ME+13_07;ME+13_08;ME+13_09;ME+13_10;ME+13_11;ME+13_12;ME+"
      "13_13;ME+13_14;ME+13_15;ME+13_16;ME+13_17;ME+13_18;ME+13_19;ME+13_20;ME+13_21;ME+13_22;ME+13_23;ME+13_24;ME+13_"
      "25;ME+13_26;ME+13_27;ME+13_28;ME+13_29;ME+13_30;ME+13_31;ME+13_32;ME+13_33;ME+13_34;ME+13_35;ME+13_36");
  makeFigure(
      "pictures/CSC_ChambersMEP14.png",
      frame,
      "ME+14_01;ME+14_02;ME+14_03;ME+14_04;ME+14_05;ME+14_06;ME+14_07;ME+14_08;ME+14_09;ME+14_10;ME+14_11;ME+14_12;ME+"
      "14_13;ME+14_14;ME+14_15;ME+14_16;ME+14_17;ME+14_18;ME+14_19;ME+14_20;ME+14_21;ME+14_22;ME+14_23;ME+14_24;ME+14_"
      "25;ME+14_26;ME+14_27;ME+14_28;ME+14_29;ME+14_30;ME+14_31;ME+14_32;ME+14_33;ME+14_34;ME+14_35;ME+14_36");
  makeFigure("pictures/CSC_ChambersMEP21.png",
             frame,
             "ME+21_01;ME+21_02;ME+21_03;ME+21_04;ME+21_05;ME+21_06;ME+21_07;ME+21_08;ME+21_09;ME+21_10;ME+21_11;ME+21_"
             "12;ME+21_13;ME+21_14;ME+21_15;ME+21_16;ME+21_17;ME+21_18");
  makeFigure(
      "pictures/CSC_ChambersMEP22.png",
      frame,
      "ME+22_01;ME+22_02;ME+22_03;ME+22_04;ME+22_05;ME+22_06;ME+22_07;ME+22_08;ME+22_09;ME+22_10;ME+22_11;ME+22_12;ME+"
      "22_13;ME+22_14;ME+22_15;ME+22_16;ME+22_17;ME+22_18;ME+22_19;ME+22_20;ME+22_21;ME+22_22;ME+22_23;ME+22_24;ME+22_"
      "25;ME+22_26;ME+22_27;ME+22_28;ME+22_29;ME+22_30;ME+22_31;ME+22_32;ME+22_33;ME+22_34;ME+22_35;ME+22_36");
  makeFigure("pictures/CSC_ChambersMEP31.png",
             frame,
             "ME+31_01;ME+31_02;ME+31_03;ME+31_04;ME+31_05;ME+31_06;ME+31_07;ME+31_08;ME+31_09;ME+31_10;ME+31_11;ME+31_"
             "12;ME+31_13;ME+31_14;ME+31_15;ME+31_16;ME+31_17;ME+31_18");
  makeFigure(
      "pictures/CSC_ChambersMEP32.png",
      frame,
      "ME+32_01;ME+32_02;ME+32_03;ME+32_04;ME+32_05;ME+32_06;ME+32_07;ME+32_08;ME+32_09;ME+32_10;ME+32_11;ME+32_12;ME+"
      "32_13;ME+32_14;ME+32_15;ME+32_16;ME+32_17;ME+32_18;ME+32_19;ME+32_20;ME+32_21;ME+32_22;ME+32_23;ME+32_24;ME+32_"
      "25;ME+32_26;ME+32_27;ME+32_28;ME+32_29;ME+32_30;ME+32_31;ME+32_32;ME+32_33;ME+32_34;ME+32_35;ME+32_36");
  makeFigure("pictures/CSC_ChambersMEP41.png",
             frame,
             "ME+41_01;ME+41_02;ME+41_03;ME+41_04;ME+41_05;ME+41_06;ME+41_07;ME+41_08;ME+41_09;ME+41_10;ME+41_11;ME+41_"
             "12;ME+41_13;ME+41_14;ME+41_15;ME+41_16;ME+41_17;ME+41_18");
  makeFigure(
      "pictures/CSC_ChambersMEP42.png",
      frame,
      "ME+42_01;ME+42_02;ME+42_03;ME+42_04;ME+42_05;ME+42_06;ME+42_07;ME+42_08;ME+42_09;ME+42_10;ME+42_11;ME+42_12;ME+"
      "42_13;ME+42_14;ME+42_15;ME+42_16;ME+42_17;ME+42_18;ME+42_19;ME+42_20;ME+42_21;ME+42_22;ME+42_23;ME+42_24;ME+42_"
      "25;ME+42_26;ME+42_27;ME+42_28;ME+42_29;ME+42_30;ME+42_31;ME+42_32;ME+42_33;ME+42_34;ME+42_35;ME+42_36");
  makeFigure(
      "pictures/CSC_ChambersMEM11.png",
      frame,
      "ME-11_01;ME-11_02;ME-11_03;ME-11_04;ME-11_05;ME-11_06;ME-11_07;ME-11_08;ME-11_09;ME-11_10;ME-11_11;ME-11_12;ME-"
      "11_13;ME-11_14;ME-11_15;ME-11_16;ME-11_17;ME-11_18;ME-11_19;ME-11_20;ME-11_21;ME-11_22;ME-11_23;ME-11_24;ME-11_"
      "25;ME-11_26;ME-11_27;ME-11_28;ME-11_29;ME-11_30;ME-11_31;ME-11_32;ME-11_33;ME-11_34;ME-11_35;ME-11_36");
  makeFigure(
      "pictures/CSC_ChambersMEM12.png",
      frame,
      "ME-12_01;ME-12_02;ME-12_03;ME-12_04;ME-12_05;ME-12_06;ME-12_07;ME-12_08;ME-12_09;ME-12_10;ME-12_11;ME-12_12;ME-"
      "12_13;ME-12_14;ME-12_15;ME-12_16;ME-12_17;ME-12_18;ME-12_19;ME-12_20;ME-12_21;ME-12_22;ME-12_23;ME-12_24;ME-12_"
      "25;ME-12_26;ME-12_27;ME-12_28;ME-12_29;ME-12_30;ME-12_31;ME-12_32;ME-12_33;ME-12_34;ME-12_35;ME-12_36");
  makeFigure(
      "pictures/CSC_ChambersMEM13.png",
      frame,
      "ME-13_01;ME-13_02;ME-13_03;ME-13_04;ME-13_05;ME-13_06;ME-13_07;ME-13_08;ME-13_09;ME-13_10;ME-13_11;ME-13_12;ME-"
      "13_13;ME-13_14;ME-13_15;ME-13_16;ME-13_17;ME-13_18;ME-13_19;ME-13_20;ME-13_21;ME-13_22;ME-13_23;ME-13_24;ME-13_"
      "25;ME-13_26;ME-13_27;ME-13_28;ME-13_29;ME-13_30;ME-13_31;ME-13_32;ME-13_33;ME-13_34;ME-13_35;ME-13_36");
  makeFigure(
      "pictures/CSC_ChambersMEM14.png",
      frame,
      "ME-14_01;ME-14_02;ME-14_03;ME-14_04;ME-14_05;ME-14_06;ME-14_07;ME-14_08;ME-14_09;ME-14_10;ME-14_11;ME-14_12;ME-"
      "14_13;ME-14_14;ME-14_15;ME-14_16;ME-14_17;ME-14_18;ME-14_19;ME-14_20;ME-14_21;ME-14_22;ME-14_23;ME-14_24;ME-14_"
      "25;ME-14_26;ME-14_27;ME-14_28;ME-14_29;ME-14_30;ME-14_31;ME-14_32;ME-14_33;ME-14_34;ME-14_35;ME-14_36");
  makeFigure("pictures/CSC_ChambersMEM21.png",
             frame,
             "ME-21_01;ME-21_02;ME-21_03;ME-21_04;ME-21_05;ME-21_06;ME-21_07;ME-21_08;ME-21_09;ME-21_10;ME-21_11;ME-21_"
             "12;ME-21_13;ME-21_14;ME-21_15;ME-21_16;ME-21_17;ME-21_18");
  makeFigure(
      "pictures/CSC_ChambersMEM22.png",
      frame,
      "ME-22_01;ME-22_02;ME-22_03;ME-22_04;ME-22_05;ME-22_06;ME-22_07;ME-22_08;ME-22_09;ME-22_10;ME-22_11;ME-22_12;ME-"
      "22_13;ME-22_14;ME-22_15;ME-22_16;ME-22_17;ME-22_18;ME-22_19;ME-22_20;ME-22_21;ME-22_22;ME-22_23;ME-22_24;ME-22_"
      "25;ME-22_26;ME-22_27;ME-22_28;ME-22_29;ME-22_30;ME-22_31;ME-22_32;ME-22_33;ME-22_34;ME-22_35;ME-22_36");
  makeFigure("pictures/CSC_ChambersMEM31.png",
             frame,
             "ME-31_01;ME-31_02;ME-31_03;ME-31_04;ME-31_05;ME-31_06;ME-31_07;ME-31_08;ME-31_09;ME-31_10;ME-31_11;ME-31_"
             "12;ME-31_13;ME-31_14;ME-31_15;ME-31_16;ME-31_17;ME-31_18");
  makeFigure(
      "pictures/CSC_ChambersMEM32.png",
      frame,
      "ME-32_01;ME-32_02;ME-32_03;ME-32_04;ME-32_05;ME-32_06;ME-32_07;ME-32_08;ME-32_09;ME-32_10;ME-32_11;ME-32_12;ME-"
      "32_13;ME-32_14;ME-32_15;ME-32_16;ME-32_17;ME-32_18;ME-32_19;ME-32_20;ME-32_21;ME-32_22;ME-32_23;ME-32_24;ME-32_"
      "25;ME-32_26;ME-32_27;ME-32_28;ME-32_29;ME-32_30;ME-32_31;ME-32_32;ME-32_33;ME-32_34;ME-32_35;ME-32_36");
  makeFigure("pictures/CSC_ChambersMEM41.png",
             frame,
             "ME-41_01;ME-41_02;ME-41_03;ME-41_04;ME-41_05;ME-41_06;ME-41_07;ME-41_08;ME-41_09;ME-41_10;ME-41_11;ME-41_"
             "12;ME-41_13;ME-41_14;ME-41_15;ME-41_16;ME-41_17;ME-41_18");
  makeFigure(
      "pictures/CSC_ChambersMEM42.png",
      frame,
      "ME-42_01;ME-42_02;ME-42_03;ME-42_04;ME-42_05;ME-42_06;ME-42_07;ME-42_08;ME-42_09;ME-42_10;ME-42_11;ME-42_12;ME-"
      "42_13;ME-42_14;ME-42_15;ME-42_16;ME-42_17;ME-42_18;ME-42_19;ME-42_20;ME-42_21;ME-42_22;ME-42_23;ME-42_24;ME-42_"
      "25;ME-42_26;ME-42_27;ME-42_28;ME-42_29;ME-42_30;ME-42_31;ME-42_32;ME-42_33;ME-42_34;ME-42_35;ME-42_36");

  makeFigure(
      "pictures/DT_Stations.png",
      frame,
      "DT-21;DT-22;DT-23;DT-24;DT-11;DT-12;DT-13;DT-14;DT01;DT02;DT03;DT04;DT11;DT12;DT13;DT14;DT21;DT22;DT23;DT24");
  makeFigure(
      "pictures/DT_ChambersM21.png",
      frame,
      "DT-21_01;DT-21_02;DT-21_03;DT-21_04;DT-21_05;DT-21_06;DT-21_07;DT-21_08;DT-21_09;DT-21_10;DT-21_11;DT-21_12");
  makeFigure(
      "pictures/DT_ChambersM22.png",
      frame,
      "DT-22_01;DT-22_02;DT-22_03;DT-22_04;DT-22_05;DT-22_06;DT-22_07;DT-22_08;DT-22_09;DT-22_10;DT-22_11;DT-22_12");
  makeFigure(
      "pictures/DT_ChambersM23.png",
      frame,
      "DT-23_01;DT-23_02;DT-23_03;DT-23_04;DT-23_05;DT-23_06;DT-23_07;DT-23_08;DT-23_09;DT-23_10;DT-23_11;DT-23_12");
  makeFigure("pictures/DT_ChambersM24.png",
             frame,
             "DT-24_01;DT-24_02;DT-24_03;DT-24_04;DT-24_05;DT-24_06;DT-24_07;DT-24_08;DT-24_09;DT-24_10;DT-24_11;DT-24_"
             "12;DT-24_13;DT-24_14");
  makeFigure(
      "pictures/DT_ChambersM11.png",
      frame,
      "DT-11_01;DT-11_02;DT-11_03;DT-11_04;DT-11_05;DT-11_06;DT-11_07;DT-11_08;DT-11_09;DT-11_10;DT-11_11;DT-11_12");
  makeFigure(
      "pictures/DT_ChambersM12.png",
      frame,
      "DT-12_01;DT-12_02;DT-12_03;DT-12_04;DT-12_05;DT-12_06;DT-12_07;DT-12_08;DT-12_09;DT-12_10;DT-12_11;DT-12_12");
  makeFigure(
      "pictures/DT_ChambersM13.png",
      frame,
      "DT-13_01;DT-13_02;DT-13_03;DT-13_04;DT-13_05;DT-13_06;DT-13_07;DT-13_08;DT-13_09;DT-13_10;DT-13_11;DT-13_12");
  makeFigure("pictures/DT_ChambersM14.png",
             frame,
             "DT-14_01;DT-14_02;DT-14_03;DT-14_04;DT-14_05;DT-14_06;DT-14_07;DT-14_08;DT-14_09;DT-14_10;DT-14_11;DT-14_"
             "12;DT-14_13;DT-14_14");
  makeFigure("pictures/DT_ChambersP01.png",
             frame,
             "DT01_01;DT01_02;DT01_03;DT01_04;DT01_05;DT01_06;DT01_07;DT01_08;DT01_09;DT01_10;DT01_11;DT01_12");
  makeFigure("pictures/DT_ChambersP02.png",
             frame,
             "DT02_01;DT02_02;DT02_03;DT02_04;DT02_05;DT02_06;DT02_07;DT02_08;DT02_09;DT02_10;DT02_11;DT02_12");
  makeFigure("pictures/DT_ChambersP03.png",
             frame,
             "DT03_01;DT03_02;DT03_03;DT03_04;DT03_05;DT03_06;DT03_07;DT03_08;DT03_09;DT03_10;DT03_11;DT03_12");
  makeFigure("pictures/DT_ChambersP04.png",
             frame,
             "DT04_01;DT04_02;DT04_03;DT04_04;DT04_05;DT04_06;DT04_07;DT04_08;DT04_09;DT04_10;DT04_11;DT04_12;DT04_13;"
             "DT04_14");
  makeFigure("pictures/DT_ChambersP11.png",
             frame,
             "DT11_01;DT11_02;DT11_03;DT11_04;DT11_05;DT11_06;DT11_07;DT11_08;DT11_09;DT11_10;DT11_11;DT11_12");
  makeFigure("pictures/DT_ChambersP12.png",
             frame,
             "DT12_01;DT12_02;DT12_03;DT12_04;DT12_05;DT12_06;DT12_07;DT12_08;DT12_09;DT12_10;DT12_11;DT12_12");
  makeFigure("pictures/DT_ChambersP13.png",
             frame,
             "DT13_01;DT13_02;DT13_03;DT13_04;DT13_05;DT13_06;DT13_07;DT13_08;DT13_09;DT13_10;DT13_11;DT13_12");
  makeFigure("pictures/DT_ChambersP14.png",
             frame,
             "DT14_01;DT14_02;DT14_03;DT14_04;DT14_05;DT14_06;DT14_07;DT14_08;DT14_09;DT14_10;DT14_11;DT14_12;DT14_13;"
             "DT14_14");
  makeFigure("pictures/DT_ChambersP21.png",
             frame,
             "DT21_01;DT21_02;DT21_03;DT21_04;DT21_05;DT21_06;DT21_07;DT21_08;DT21_09;DT21_10;DT21_11;DT21_12");
  makeFigure("pictures/DT_ChambersP22.png",
             frame,
             "DT22_01;DT22_02;DT22_03;DT22_04;DT22_05;DT22_06;DT22_07;DT22_08;DT22_09;DT22_10;DT22_11;DT22_12");
  makeFigure("pictures/DT_ChambersP23.png",
             frame,
             "DT23_01;DT23_02;DT23_03;DT23_04;DT23_05;DT23_06;DT23_07;DT23_08;DT23_09;DT23_10;DT23_11;DT23_12");
  makeFigure("pictures/DT_ChambersP24.png",
             frame,
             "DT24_01;DT24_02;DT24_03;DT24_04;DT24_05;DT24_06;DT24_07;DT24_08;DT24_09;DT24_10;DT24_11;DT24_12;DT24_13;"
             "DT24_14");

  for (unsigned int L = 0; L < 3; L++) {
    TGaxis::SetMaxDigits(2);
    TCanvas* c1 = new TCanvas("c1", "c1,", 600, 600);
    c1->SetLeftMargin(0.16);
    //c1->SetLogy(true);
    frame = new TH1D("frame", "frame", 1, 0.5, 1.5);
    frame->SetTitle("");
    frame->SetStats(kFALSE);
    frame->GetXaxis()->SetNdivisions(505);
    if (L == 0)
      frame->GetXaxis()->SetTitle("1/#beta_{DT}");
    if (L == 1)
      frame->GetXaxis()->SetTitle("1/#beta_{CSC}");
    if (L == 2)
      frame->GetXaxis()->SetTitle("1/#beta_{DT+CSC}");
    frame->GetXaxis()->SetTitleOffset(1.25);
    frame->GetYaxis()->SetTitle("Muons");
    frame->GetYaxis()->SetTitleOffset(1.75);
    frame->SetMaximum(0.2);
    frame->SetMinimum(0.1);
    frame->Draw("AXIS");

    TLegend* LEG = new TLegend(0.18, 0.60, 0.60, 0.90);
    LEG->SetNColumns(1);
    LEG->SetFillColor(0);
    LEG->SetFillStyle(0);
    LEG->SetBorderSize(0);

    if (L == 0) {
      TH1D* HiBetaAOD = (TH1D*)GetObjectFromPath(InputFile, "DT_iBeta_AOD", false);
      HiBetaAOD->SetLineColor(1);
      HiBetaAOD->SetLineWidth(2);
      HiBetaAOD->Draw("same");
      LEG->AddEntry(HiBetaAOD, "AOD default", "L");
      TH1D* HiBetaFLY0 = (TH1D*)GetObjectFromPath(InputFile, "DT_iBeta_FLY0", false);
      HiBetaFLY0->SetLineColor(4);
      HiBetaFLY0->SetLineWidth(2);
      HiBetaFLY0->Draw("same");
      LEG->AddEntry(HiBetaFLY0, "FWlite: no t0 corr.", "L");
      TH1D* HiBetaFLY1 = (TH1D*)GetObjectFromPath(InputFile, "DT_iBeta_FLY1", false);
      HiBetaFLY1->SetLineColor(2);
      HiBetaFLY1->SetLineWidth(2);
      HiBetaFLY1->Draw("same");
      LEG->AddEntry(HiBetaFLY1, "FWlite: station t0 corr.", "L");
      TH1D* HiBetaFLY2 = (TH1D*)GetObjectFromPath(InputFile, "DT_iBeta_FLY2", false);
      HiBetaFLY2->SetLineColor(8);
      HiBetaFLY2->SetLineWidth(2);
      HiBetaFLY2->Draw("same");
      LEG->AddEntry(HiBetaFLY2, "FWlite: chamber t0 corr.", "L");
      TH1D* HiBetaFLY3 = (TH1D*)GetObjectFromPath(InputFile, "DT_iBeta_FLY3", false);
      HiBetaFLY3->SetLineColor(9);
      HiBetaFLY3->SetLineWidth(2);
      HiBetaFLY3->Draw("same");
      LEG->AddEntry(HiBetaFLY3, "same with tighter pruning", "L");
      frame->SetMaximum(HiBetaFLY3->GetMaximum() * 1.2);
    } else if (L == 1) {
      TH1D* HiBetaAOD = (TH1D*)GetObjectFromPath(InputFile, "CSC_iBeta_AOD", false);
      HiBetaAOD->SetLineColor(1);
      HiBetaAOD->SetLineWidth(2);
      HiBetaAOD->Draw("same");
      LEG->AddEntry(HiBetaAOD, "AOD default", "L");
      TH1D* HiBetaFLY0 = (TH1D*)GetObjectFromPath(InputFile, "CSC_iBeta_FLY0", false);
      HiBetaFLY0->SetLineColor(4);
      HiBetaFLY0->SetLineWidth(2);
      HiBetaFLY0->Draw("same");
      LEG->AddEntry(HiBetaFLY0, "FWlite: no t0 corr.", "L");
      TH1D* HiBetaFLY1 = (TH1D*)GetObjectFromPath(InputFile, "CSC_iBeta_FLY1", false);
      HiBetaFLY1->SetLineColor(2);
      HiBetaFLY1->SetLineWidth(2);
      HiBetaFLY1->Draw("same");
      LEG->AddEntry(HiBetaFLY1, "FWlite: station t0 corr.", "L");
      TH1D* HiBetaFLY2 = (TH1D*)GetObjectFromPath(InputFile, "CSC_iBeta_FLY2", false);
      HiBetaFLY2->SetLineColor(8);
      HiBetaFLY2->SetLineWidth(2);
      HiBetaFLY2->Draw("same");
      LEG->AddEntry(HiBetaFLY2, "FWlite: chamber t0 corr.", "L");
      TH1D* HiBetaFLY3 = (TH1D*)GetObjectFromPath(InputFile, "CSC_iBeta_FLY3", false);
      HiBetaFLY3->SetLineColor(9);
      HiBetaFLY3->SetLineWidth(2);
      HiBetaFLY3->Draw("same");
      LEG->AddEntry(HiBetaFLY3, "same with tighter pruning", "L");
      frame->SetMaximum(HiBetaFLY3->GetMaximum() * 1.2);
    } else {
      TH1D* HiBetaAOD = (TH1D*)GetObjectFromPath(InputFile, "iBeta_AOD", false);
      HiBetaAOD->SetLineColor(1);
      HiBetaAOD->SetLineWidth(2);
      HiBetaAOD->Draw("same");
      LEG->AddEntry(HiBetaAOD, "AOD default", "L");
      TH1D* HiBetaFLY0 = (TH1D*)GetObjectFromPath(InputFile, "iBeta_FLY0", false);
      HiBetaFLY0->SetLineColor(4);
      HiBetaFLY0->SetLineWidth(2);
      HiBetaFLY0->Draw("same");
      LEG->AddEntry(HiBetaFLY0, "FWlite: no t0 corr.", "L");
      TH1D* HiBetaFLY1 = (TH1D*)GetObjectFromPath(InputFile, "iBeta_FLY1", false);
      HiBetaFLY1->SetLineColor(2);
      HiBetaFLY1->SetLineWidth(2);
      HiBetaFLY1->Draw("same");
      LEG->AddEntry(HiBetaFLY1, "FWlite: station t0 corr.", "L");
      TH1D* HiBetaFLY2 = (TH1D*)GetObjectFromPath(InputFile, "iBeta_FLY2", false);
      HiBetaFLY2->SetLineColor(8);
      HiBetaFLY2->SetLineWidth(2);
      HiBetaFLY2->Draw("same");
      LEG->AddEntry(HiBetaFLY2, "FWlite: chamber t0 corr.", "L");
      TH1D* HiBetaFLY3 = (TH1D*)GetObjectFromPath(InputFile, "iBeta_FLY3", false);
      HiBetaFLY3->SetLineColor(9);
      HiBetaFLY3->SetLineWidth(2);
      HiBetaFLY3->Draw("same");
      LEG->AddEntry(HiBetaFLY3, "same with tighter pruning", "L");
      frame->SetMaximum(HiBetaFLY3->GetMaximum() * 1.2);
    }

    TLine line(1.0, frame->GetMinimum(), 1.0, frame->GetMaximum());
    line.SetLineColor(1);
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    line.Draw();

    LEG->Draw();
    if (L == 0)
      c1->SaveAs("pictures/iBeta_DT.png");
    if (L == 1)
      c1->SaveAs("pictures/iBeta_CSC.png");
    if (L == 2)
      c1->SaveAs("pictures/iBeta_Combined.png");
    delete c1;
    delete LEG;
  }

  if (true) {
    TGaxis::SetMaxDigits(2);
    TCanvas* c1 = new TCanvas("c1", "c1,", 600, 600);
    c1->SetLeftMargin(0.16);
    c1->SetRightMargin(0.16);
    c1->SetLogz(true);
    frame = new TH1D("frame", "frame", 1, 0.5, 1.5);
    frame->SetTitle("");
    frame->SetStats(kFALSE);
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetXaxis()->SetTitle("AOD default");
    frame->GetXaxis()->SetTitleOffset(1.25);
    frame->GetYaxis()->SetTitle("FWlite: chamber t0 corr.");
    frame->GetYaxis()->SetTitleOffset(1.75);
    frame->SetMaximum(1.5);
    frame->SetMinimum(0.5);
    frame->Draw("AXIS");

    TH2D* H_iBetaAODvsFly = (TH2D*)GetObjectFromPath(InputFile, "iBeta_AODvsFly", false);
    H_iBetaAODvsFly->Draw("same COLZ");
    TLine line(0.5, 0.5, 1.5, 1.5);
    line.SetLineColor(1);
    line.SetLineWidth(2);
    line.SetLineStyle(2);
    line.Draw();

    c1->SaveAs("pictures/iBeta_Correlation.png");
    delete c1;
  }
}

/////////////////////////// /////////////////////////// ///////////////////////////
/////////////////////////// ADDITIONAL PLOTTING FUNCTIONS /////////////////////////
/////////////////////////// /////////////////////////// ///////////////////////////

// handfull function to get one TObject from a complex cirectory stucture in a file
TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy) {
  size_t pos = Path.find("/");
  if (pos < 256) {
    std::string firstPart = Path.substr(0, pos);
    std::string endPart = Path.substr(pos + 1, Path.length());
    TDirectory* TMP = (TDirectory*)File->Get(firstPart.c_str());
    if (TMP != NULL)
      return GetObjectFromPath(TMP, endPart, GetACopy);

    printf("ObjectNotFound: %s::%s\n", File->GetName(), Path.c_str());
    return NULL;
  } else {
    if (GetACopy) {
      return (File->Get(Path.c_str()))->Clone();
    } else {
      return File->Get(Path.c_str());
    }
  }
}
