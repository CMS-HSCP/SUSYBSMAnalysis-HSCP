// Original Author:  Loic Quertenmont

// ~~~~~~~~~c++ include files ~~~~~~~~~
#include <memory>
#include <vector>
#include <iostream>
#include <fstream>
#include <string>
#include <map>
#include <exception>
#include <unordered_map>

// ~~~~~~~~~ ROOT include files ~~~~~~~~~
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TDirectory.h"
#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TObject.h"
#include "TVector3.h"
#include "TChain.h"
#include "TRandom3.h"
#include "TProfile.h"
#include "TCanvas.h"

#include <TBenchmark.h>
#include <TSystemDirectory.h>
#include <TSystem.h>
//#include "tdrstyle.C"

#include "../interface/CommonFunction.h"
#include "HistoTools.h"
#include "../interface/TuplePlotter.h"
#include "tdrstyle.h"
#include "ArgumentParser.h"

class ArgumentParser;

using namespace std;

/////////////////////////// MAIN FUNCTION /////////////////////////////

int main(int argc, char* argv[]) {
  string usage = "Usage: LimitComputation --inputFiles <file.txt> [--mode <value>]\n";

  vector<string> input;
  int typeMode = 0;

  ArgumentParser parser(argc, argv);

  if (parser.findOption("-h")) {
    cout << usage << endl;
    return 0;
  }
  if (parser.findOption("--inputFiles"))
    parser.getArgument("--inputFiles", input);
  else {
    cout << usage << endl;
    return 0;
  }
  if (parser.findOption("--mode"))
    parser.getArgument("--mode", typeMode);

  cout << "==================" << endl;
  cout << " LimitComputation " << endl;
  cout << "==================\n" << endl;

  setTDRStyle();
  /*gStyle->SetPadTopMargin   (0.05);
    gStyle->SetPadBottomMargin(0.12);
    gStyle->SetPadRightMargin (0.05);
    gStyle->SetPadLeftMargin  (0.14);
    gStyle->SetTitleSize(0.05, "XYZ");
    gStyle->SetTitleXOffset(1.1);
    gStyle->SetTitleYOffset(1.4);
    gStyle->SetPalette(1);
    gStyle->SetNdivisions(505);
    //gStyle->SetTextFont(43);*/

  TBenchmark clock;
  clock.Start("LimitComputation");

  vector<string> inputFiles;
  if (endsWith(input[0], ".txt")) {
    ifstream fin(input[0], std::ios::binary);
    if (not fin.is_open()) {
      cout << "Failed to open " << input[0] << endl;
      return 0;
    }
    string rootfile;
    while (fin >> rootfile)
      inputFiles.push_back(rootfile);
  } else
    inputFiles = input;

  for (auto const& inputFile : inputFiles) {
    cout << "opening file " << inputFile << endl;
    //TFile* inFile = TFile::Open(inputFile.c_str(),"UPDATE");
    TFile* tfile = TFile::Open(inputFile.c_str());
    if (not tfile) {
      cout << "Failed to open " << inputFile << endl;
      return 0;
    }
    //Analysis_Step4_LimitComputation(string MODE="COMPILE", string InputPattern="", string signal="");//(TFile* InputFile, int TypeMode)
  }

  cout << "" << endl;
  clock.Show("LimitComputation");
  cout << "" << endl;

  return 0;
}