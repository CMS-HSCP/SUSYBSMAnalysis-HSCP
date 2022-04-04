#include "../../AnalysisCode/Analysis_Global.h"
#include "../../AnalysisCode/Analysis_CommonFunction.h"
#include "../../AnalysisCode/Analysis_PlotFunction.h"
#include "../../AnalysisCode/Analysis_Samples.h"
#include "../../AnalysisCode/tdrstyle.C"

#include "TGraphAsymmErrors.h"

using namespace std;

size_t GetNumberOfLines(FILE* fin) {
  size_t N = 0;
  char c;
  do {
    c = fgetc(fin);
    if (c == '\n')
      N++;
  } while (c != EOF);
  rewind(fin);
  return N;
}

string PatternFromType(unsigned int TypeMode) {
  char tmp[128];
  sprintf(tmp, "Type%d/", TypeMode);
  return string(tmp);
}

class TxtInfo {
private:
  double Mass_, Eff_, XSec_Th_, XSec_Exp_, XSec_Obs_, NData_, NPred_, NPredErr_, NSign_, LInt_;

public:
  TxtInfo(string AnalysisPath, string Signal, unsigned int TypeMode, bool debug = false) {
    FILE* fin = fopen(
        (AnalysisPath + "Results/" + PatternFromType(TypeMode) + "EXCLUSION13TeV/" + Signal + ".txt").c_str(), "r");
    //	 if (!fin) fin = fopen ((AnalysisPath+"Results/"+PatternFromType(TypeMode)+"EXCLUSION13TeV15/"+Signal+".txt").c_str(), "r");
    if (!fin)
      fin = fopen(
          (AnalysisPath + "Results/" + PatternFromType(TypeMode) + "EXCLUSION13TeV16/" + Signal + ".txt").c_str(), "r");
    if (!fin)
      fprintf(stderr,
              "Unable to open EXCLUSION dir in the analysis path %s for signal %s.txt\n",
              AnalysisPath.c_str(),
              Signal.c_str());
    Mass_ = -1;
    Eff_ = -1;
    XSec_Th_ = -1;
    XSec_Exp_ = -1;
    XSec_Obs_ = -1;
    NData_ = -1;
    NPred_ = -1;
    NPredErr_ = -1;
    NSign_ = -1;
    LInt_ = -1;

    if (fin) {
      size_t N = GetNumberOfLines(fin);
      char VarName[128];
      char VarValue[128];
      for (size_t line = 0; line <= N; line++) {
        fscanf(fin, "%s : %s +- %*s\n", VarName, VarValue);

        if (debug)
          printf("%s\t%s\n", VarName, VarValue);

        if (strcmp(VarName, "Mass") == 0)
          Mass_ = atof(VarValue);
        else if (strcmp(VarName, "Eff") == 0)
          Eff_ = atof(VarValue);
        else if (strcmp(VarName, "XSec_Th") == 0)
          XSec_Th_ = atof(VarValue);
        else if (strcmp(VarName, "XSec_Exp") == 0)
          XSec_Exp_ = atof(VarValue);
        else if (strcmp(VarName, "XSec_Obs") == 0)
          XSec_Obs_ = atof(VarValue);
        else if (strcmp(VarName, "NData") == 0)
          NData_ = atof(VarValue);
        else if (strcmp(VarName, "NPred") == 0)
          NPred_ = atof(VarValue);
        else if (strcmp(VarName, "NPredErr") == 0)
          NPredErr_ = atof(VarValue);
        else if (strcmp(VarName, "NSign") == 0)
          NSign_ = atof(VarValue);
        else if (strcmp(VarName, "LInt") == 0)
          LInt_ = atof(VarValue);
      }
      fclose(fin);
    }
  }

  void PrintSummary() {
    cout << "Mass     :" << Mass_ << endl;
    cout << "Eff      :" << Eff_ << endl;
    cout << "XSec_Th  :" << XSec_Th_ << endl;
    cout << "XSec_Exp :" << XSec_Exp_ << endl;
    cout << "XSec_Obs :" << XSec_Obs_ << endl;
    cout << "NData    :" << NData_ << endl;
    cout << "NPred    :" << NPred_ << endl;
    cout << "NPredErr :" << NPredErr_ << endl;
    cout << "NSign    :" << NSign_ << endl;
    cout << "LInt     :" << LInt_ << endl;
  }

  double Mass() { return Mass_; }
  double Eff() { return Eff_; }
  double XSec_Th() { return XSec_Th_; }
  double XSec_Exp() { return XSec_Exp_; }
  double XSec_Obs() { return XSec_Obs_; }
  double NData() { return NData_; }
  double NPred() { return NPred_; }
  double NPredErr() { return NPredErr_; }
  double NSign() { return NSign_; }
  double LInt() { return LInt_; }
};

class SignalSummaryClass {
private:
  vector<TxtInfo*> Infos;
  size_t N;

public:
  SignalSummaryClass(string AnalysisPath, string ModelName, unsigned int TypeMode, bool debug = false) {
    vector<stSample> samples;
    InitBaseDirectory();
    GetSampleDefinition(samples, "../../AnalysisCode/Analysis_Samples.txt");
    for (size_t s = 0; s < samples.size(); s++) {
      if (samples[s].ModelName() == ModelName) {
        Infos.push_back(new TxtInfo(AnalysisPath, samples[s].Name, TypeMode, debug));
      }
    }

    N = Infos.size();
  }

  ~SignalSummaryClass() {
    for (size_t i = 0; i < Infos.size(); i++)
      delete Infos[i];
    Infos.clear();
  }

  vector<double> MassVector() {
    vector<double> toReturn;
    for (auto it = Infos.begin(); it != Infos.end(); it++)
      if ((*it)->Mass() > 0)
        toReturn.push_back((*it)->Mass());

    return toReturn;
  }

  vector<double> EffVector() {
    vector<double> toReturn;
    for (auto it = Infos.begin(); it != Infos.end(); it++)
      if ((*it)->Mass() > 0)
        toReturn.push_back((*it)->Eff());

    return toReturn;
  }

  vector<double> XSecThVector() {
    vector<double> toReturn;
    for (auto it = Infos.begin(); it != Infos.end(); it++)
      if ((*it)->Mass() > 0)
        toReturn.push_back((*it)->XSec_Th());

    return toReturn;
  }

  vector<double> XSecExpVector() {
    vector<double> toReturn;
    for (auto it = Infos.begin(); it != Infos.end(); it++)
      if ((*it)->Mass() > 0)
        toReturn.push_back((*it)->XSec_Exp());

    return toReturn;
  }

  vector<double> XSecObsVector() {
    vector<double> toReturn;
    for (auto it = Infos.begin(); it != Infos.end(); it++)
      if ((*it)->Mass() > 0)
        toReturn.push_back((*it)->XSec_Obs());

    return toReturn;
  }

  size_t size() { return N; }
};

void DrawFamilyOfGraphs(vector<TGraph*> graphs,
                        TH1D* tmp,
                        string style,
                        string xTitle,
                        string yTitle,
                        bool logy = true,
                        double xmin = -1,
                        double xmax = -1,
                        double ymin = -1,
                        double ymax = -1) {
  if (xmin == -1 || xmax == -1 || ymin == -1 || ymax == -1) {
    double XminOld = 999999, XmaxOld = 0, YminOld = 999999, YmaxOld = 0;
    double Xmin, Xmax, Ymin, Ymax;
    for (auto graph = graphs.begin(); graph != graphs.end(); graph++) {
      (*graph)->ComputeRange(Xmin, Ymin, Xmax, Ymax);
      Xmin = std::min(Xmin, XminOld);
      Xmax = std::max(Xmax, XmaxOld);
      Ymin = std::min(Ymin, YminOld);
      Ymax = std::max(Ymax, YmaxOld);

      XminOld = Xmin;
      XmaxOld = Xmax;
      YminOld = Ymin;
      YmaxOld = Ymax;
    }

    if (xmin == -1)
      xmin = Xmin;
    if (ymin == -1)
      ymin = Ymin;
    if (xmax == -1)
      xmax = Xmax;
    if (ymax == -1)
      ymax = Ymax;
  }

  if (ymin < 0)
    ymin = 1e-5;
  tmp->GetXaxis()->SetLimits(xmin, xmax);
  tmp->SetTitle("");
  tmp->SetStats(0);
  tmp->GetXaxis()->SetTitle(xTitle.c_str());
  tmp->GetYaxis()->SetTitle(yTitle.c_str());
  tmp->GetXaxis()->SetRangeUser(xmin, xmax);
  tmp->GetYaxis()->SetRangeUser(ymin * 0.1, ymax * (logy ? 12 : 1.5));
  tmp->Draw();

  for (size_t g = 0; g < graphs.size(); g++) {
    graphs[g]->SetLineWidth(2);
    graphs[g]->SetLineColor(Color[g]);
    graphs[g]->SetMarkerColor(Color[g]);
    graphs[g]->SetMarkerStyle(Marker[g]);
    graphs[g]->Draw((string("same ") + style).c_str());
  }
}

void Plotter(void) {
  setTDRStyle();
  vector<string> AnalysesPaths;
  vector<string> legend;
  vector<string> SignalsToProcess;
  vector<string> SignalsLegend;

  AnalysesPaths.push_back("../../AnalysisCode_Paper/");
  legend.push_back("AN 2015");
  AnalysesPaths.push_back("../../AnalysisCode/");
  legend.push_back("Data 2016");

  SignalsToProcess.push_back("Gluino_f10");
  SignalsLegend.push_back("Gluino, f = 10");
  SignalsToProcess.push_back("GluinoN_f10");
  SignalsLegend.push_back("Gluino, CS, f = 10");
  SignalsToProcess.push_back("Gluino_f50");
  SignalsLegend.push_back("Gluino, f = 50");
  SignalsToProcess.push_back("Stop");
  SignalsLegend.push_back("Stop");
  SignalsToProcess.push_back("StopN");
  SignalsLegend.push_back("Stop, CS");
  SignalsToProcess.push_back("GMStau");
  SignalsLegend.push_back("GMSBStau");
  SignalsToProcess.push_back("PPStau");
  SignalsLegend.push_back("PPStau");
  SignalsToProcess.push_back("DY_Q1");
  SignalsLegend.push_back("DY: |Q|=1");
  SignalsToProcess.push_back("DY_Q2");
  SignalsLegend.push_back("DY: |Q|=1");

  for (size_t s = 0; s < SignalsToProcess.size(); s++) {
    vector<TGraph*> graphsXSecObs;
    vector<TGraph*> graphsXSecExp;
    vector<TGraph*> graphsEff;
    for (size_t a = 0; a < AnalysesPaths.size(); a++) {
      SignalSummaryClass test(AnalysesPaths[a], SignalsToProcess[s], 0);
      graphsXSecObs.push_back(new TGraph(test.size(), &(test.MassVector()[0]), &(test.XSecObsVector()[0])));
      graphsXSecExp.push_back(new TGraph(test.size(), &(test.MassVector()[0]), &(test.XSecExpVector()[0])));
      graphsEff.push_back(new TGraph(test.size(), &(test.MassVector()[0]), &(test.EffVector()[0])));
    }

    TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
    TH1D tmp("tmp", "tmp", 1, 0, 99999);
    c1->SetLogy(true);
    DrawFamilyOfGraphs(graphsXSecObs, &tmp, "LP", "Mass (GeV)", "95% CL limit on #sigma_{obs} (pb)");
    DrawPreliminary("Tracker - Only", 13.0, SignalsLegend[s]);
    DrawLegend((TObject**)&(graphsXSecObs[0]), legend, "", "LP", 0.8, 0.9, 0.3, 0.05);
    SaveCanvas(c1, "TkOnly", SignalsToProcess[s] + "_XSecObs", true);
    delete c1;

    c1 = new TCanvas("c1", "c1", 600, 600);
    c1->SetLogy(true);
    DrawFamilyOfGraphs(graphsXSecExp, &tmp, "LP", "Mass (GeV)", "95% CL limit on #sigma_{exp} (pb)");
    DrawPreliminary("Tracker - Only", 13.0, SignalsLegend[s]);
    DrawLegend((TObject**)&(graphsXSecExp[0]), legend, "", "LP", 0.8, 0.9, 0.3, 0.05);
    SaveCanvas(c1, "TkOnly", SignalsToProcess[s] + "_XSecExp", true);
    delete c1;

    c1 = new TCanvas("c1", "c1", 600, 600);
    DrawFamilyOfGraphs(graphsEff, &tmp, "LP", "Mass (GeV)", "Efficiency", false);
    DrawPreliminary("Tracker - Only", 13.0, SignalsLegend[s]);
    DrawLegend((TObject**)&(graphsEff[0]), legend, "", "LP", 0.8, 0.9, 0.3, 0.05);
    SaveCanvas(c1, "TkOnly", SignalsToProcess[s] + "_Eff", true);
    delete c1;

    for (size_t g = 0; g < graphsXSecObs.size(); g++)
      delete graphsXSecObs[g];
    graphsXSecObs.clear();

    for (size_t g = 0; g < graphsXSecExp.size(); g++)
      delete graphsXSecExp[g];
    graphsXSecExp.clear();

    for (size_t g = 0; g < graphsEff.size(); g++)
      delete graphsEff[g];
    graphsEff.clear();
  }
}
