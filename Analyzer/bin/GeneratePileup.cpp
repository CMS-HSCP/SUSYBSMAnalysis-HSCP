#include <memory>
#include <iostream>
#include <string>
#include <vector>
#include <fstream>

#include "TFile.h"
#include "TBenchmark.h"
#include "TVector3.h"

#include "../interface/CommonFunction.h"
#include "ArgumentParser.h"

#include <boost/filesystem.hpp>
boost::filesystem::path BASE(__FILE__);
#define __FILENAME__ BASE.stem().c_str()

using namespace std;

void GeneratePileup(string infile, string key_word);

int main(int argc, char* argv[]) {
  string usage = "Usage: " + (string)__FILENAME__ + " -f file1 file2 fileN --keyword key_word\n";

  vector<string> input;
  string key_word = "mix.input.nbPileupEvents.probValue";

  ArgumentParser parser(argc, argv);

  if (parser.findOption("-h")) {
    cout << usage << endl;
    return 0;
  }
  if (parser.findOption("-f"))
    parser.getArgument("-f", input);
  else {
    cout << usage << endl;
    return 0;
  }
  if (parser.findOption("--keyword"))
    parser.getArgument("--keyword", key_word);

  cout << "======================" << endl;
  cout << " " << __FILENAME__ << " " << endl;
  cout << "======================\n" << endl;

  TBenchmark clock;
  clock.Start(__FILENAME__);

  for (auto const& file : input) {
    GeneratePileup(file, key_word);
  }
  cout << "" << endl;
  clock.Show(__FILENAME__);
  cout << "" << endl;

  return 0;
}

void GeneratePileup(string infile, string key_word) {
  ifstream ifs(infile, ios_base::in);

  if (!ifs) {
    cerr << "Can't open file '" << infile << "'\n";
    return;
  }
  cout << "Reading file " << infile << endl;

  vector<float> pu;

  bool found(false);
  for (string line; getline(ifs, line);) {
    if (line.find(key_word) != string::npos) {
      found = true;
    }
    if (found && line.find(",") != string::npos) {
      split(line, ',', pu);
    }
    if (found && line.find(")") != string::npos)
      break;
  }
  Int_t Nbins = pu.size();

  TString rootfile = infile, suffix = "-" + to_string(Nbins) + "bins.root";
  rootfile = rootfile.ReplaceAll(".py", suffix);

  TFile* fo = new TFile(rootfile, "recreate");

  TH1F* h = new TH1F("pileup", "pileup", Nbins, 0, float(Nbins));
  for (int ibin = 1; ibin < Nbins + 1; ibin++) {
    h->SetBinContent(ibin, pu[ibin - 1]);
  }

  h->Write();
  fo->Close();
}