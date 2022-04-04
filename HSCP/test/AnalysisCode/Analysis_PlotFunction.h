// Original Author:  Loic Quertenmont

#ifndef PLOT_FUNCTION
#define PLOT_FUNCTION

int Color[] = {1, 4, 2, 8, 6, 9, 3, 7, 5, 46, 42};
int Marker[] = {20, 22, 21, 23, 29, 27, 2, 30, 24, 25};
int Style[] = {1, 2, 5, 7, 9, 10, 11, 12, 13, 14};
int GraphStyle[] = {20, 21, 22, 23, 24, 25, 26, 27};

// handfull function to get one TObject from a complex cirectory stucture in a file
TObject* GetObjectFromPath(TDirectory* File, std::string Path, bool GetACopy = false) {
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

// similar to the above code
TObject* GetObjectFromPath(TDirectory* Container, TDirectory* File, std::string Path, bool GetACopy = false) {
  TObject* toreturn = GetObjectFromPath(File, Path, GetACopy);
  if (TH1* th1 = dynamic_cast<TH1*>(toreturn))
    th1->SetDirectory(Container);
  return toreturn;
}

TH1D* GetProjectionFromPath(TDirectory* File, std::string Path, int CutIndex, string Name) {
  TH2D* tmp = (TH2D*)GetObjectFromPath(File, Path, false);
  if (!tmp)
    return NULL;
  return tmp->ProjectionY(Name.c_str(), CutIndex + 1, CutIndex + 1, "o");
}

// create a directory/subdirectory on disk
void MakeDirectories(std::string path) { system((std::string("mkdir -p ") + path).c_str()); }

// save a TCanvas on disk in a few different format (mind that 2D plots can be huge if saved in eps/C/pdf)
void SaveCanvas(TCanvas* c, std::string path, std::string name, bool OnlyPPNG = false) {
  std::string tmppath = path;
  if (tmppath[tmppath.length() - 1] != '/')
    tmppath += "_";
  tmppath += name;

  std::string filepath;
  filepath = tmppath + ".png";
  c->SaveAs(filepath.c_str());
  if (OnlyPPNG)
    return;
  filepath = tmppath + ".eps";
  c->SaveAs(filepath.c_str());
  filepath = tmppath + ".C";
  c->SaveAs(filepath.c_str());
  filepath = tmppath + ".pdf";
  c->SaveAs(filepath.c_str());
}

// function that add the TPaveText on the current canvas with the "CMS Preliminary...." on top of the Histograms. For split Lumi
void DrawPreliminary(string Text,
                     double SQRTS_,
                     string LumiText,
                     string MoreText,
                     double L,
                     double B,
                     double R,
                     double T,
                     bool preliminary = true) {
  //TOP RIGHT OUT-FRAME
  TPaveText* T1 = new TPaveText(1.0 - R - 0.50, 1.0 - T, 1.01 - R, 1.0, "NDC");
  T1->SetTextFont(43);
  T1->SetTextSize(23);
  T1->SetTextAlign(32);
  T1->SetFillColor(0);
  T1->SetFillStyle(0);
  T1->SetBorderSize(0);
  T1->AddText(LumiText.c_str());
  T1->Draw();

  //TOP LEFT IN-FRAME
  TPaveText* T2 = new TPaveText(L + 0.01, 1.0 - T - 0.06, L + 0.20, 1.0 - T - 0.01, "NDC");
  T2->SetTextFont(63);
  T2->SetTextSize(30);
  T2->SetTextAlign(12);
  T2->SetFillColor(0);
  T2->SetFillStyle(0);
  T2->SetBorderSize(0);
  T2->AddText("CMS");
  T2->Draw();

  if (preliminary) {  //Bellow CMS
    TPaveText* T3 = new TPaveText(L + 0.01, 1.0 - T - 0.11, L + 0.20, 1.0 - T - 0.05, "NDC");
    T3->SetTextFont(53);
    T3->SetTextSize(23);
    T3->SetTextAlign(12);
    T3->SetFillColor(0);
    T3->SetFillStyle(0);
    T3->SetBorderSize(0);
    T3->AddText("Preliminary");
    T3->Draw();
  }

  if (Text != "") {  //TOP right IN-FRAME
    TPaveText* T4 = new TPaveText(1.0 - R - 0.50, 1.0 - T - 0.06, 1.0 - R, 1.0 - T - 0.01, "NDC");
    T4->SetTextFont(43);
    T4->SetTextSize(23);
    T4->SetTextAlign(32);
    T4->SetFillColor(0);
    T4->SetFillStyle(0);
    T4->SetBorderSize(0);
    T4->AddText(Text.c_str());
    T4->Draw();
  }

  if (MoreText != "") {  //Right of CMS
    TPaveText* T5 = new TPaveText(L + 0.12, 1.0 - T - 0.06, L + 0.32, 1.0 - T - 0.02, "NDC");
    T5->SetTextFont(53);
    T5->SetTextSize(23);
    T5->SetTextAlign(12);
    T5->SetFillColor(0);
    T5->SetFillStyle(0);
    T5->SetBorderSize(0);
    T5->AddText(MoreText.c_str());
    T5->Draw();
  }
}

void DrawPreliminary(
    string Text, double SQRTS_, string LumiText, string MoreText = "", TAttPad* pad = NULL, bool preliminary = true) {
  if (pad)
    DrawPreliminary(Text,
                    SQRTS_,
                    LumiText,
                    MoreText,
                    pad->GetLeftMargin(),
                    pad->GetBottomMargin(),
                    pad->GetRightMargin(),
                    pad->GetTopMargin(),
                    preliminary);
  else if (gPad)
    DrawPreliminary(Text,
                    SQRTS_,
                    LumiText,
                    MoreText,
                    gPad->GetLeftMargin(),
                    gPad->GetBottomMargin(),
                    gPad->GetRightMargin(),
                    gPad->GetTopMargin(),
                    preliminary);
  else
    DrawPreliminary(Text, SQRTS_, LumiText, MoreText, 0.15, 0.15, 0.15, 0.15, preliminary);
}

// handfull function to draw the legend associated to a vector of histogram
void DrawLegend(TObject** Histos,
                std::vector<std::string> legend,
                std::string Title,
                std::string Style_,
                double X = 0.79,
                double Y = 0.92,
                double W = 0.20,
                double H = 0.05) {
  int N = legend.size();

  if (legend[0] != "") {
    TLegend* leg;
    leg = new TLegend(X, Y, X - W, Y - N * H);
    leg->SetFillStyle(0);
    leg->SetBorderSize(0);
    leg->SetTextFont(43);  //give the font size in pixel (instead of fraction)
    leg->SetTextSize(20);  //font size
                           //      leg->SetTextAlign(11);
    if (Title != "")
      leg->SetHeader(Title.c_str());

    if (Style_ == "DataMC") {
      for (int i = 0; i < N; i++) {
        TH2D* temp = (TH2D*)Histos[i]->Clone();
        temp->SetMarkerSize(1.3);
        if (i == 0) {
          leg->AddEntry(temp, legend[i].c_str(), "P");
        } else {
          leg->AddEntry(temp, legend[i].c_str(), "L");
        }
      }
    } else {
      for (int i = 0; i < N; i++) {
        TH2D* temp = (TH2D*)Histos[i]->Clone();
        temp->SetMarkerSize(1.3);
        leg->AddEntry(temp, legend[i].c_str(), Style_.c_str());
      }
    }
    leg->Draw();
  }
}

// draw the stat box
void DrawStatBox(TObject** Histos,
                 std::vector<std::string> legend,
                 bool Mean,
                 double X = 0.15,
                 double Y = 0.93,
                 double W = 0.15,
                 double H = 0.03) {
  int N = legend.size();
  char buffer[255];

  if (Mean)
    H *= 3;
  for (int i = 0; i < N; i++) {
    TPaveText* stat = new TPaveText(X, Y - (i * H), X + W, Y - (i + 1) * H, "NDC");
    TH1* Histo = (TH1*)Histos[i];
    sprintf(buffer, "Entries : %i\n", (int)Histo->GetEntries());
    stat->AddText(buffer);

    if (Mean) {
      sprintf(buffer, "Mean    : %6.2f\n", Histo->GetMean());
      stat->AddText(buffer);

      sprintf(buffer, "RMS     : %6.2f\n", Histo->GetRMS());
      stat->AddText(buffer);
    }

    stat->SetFillColor(0);
    stat->SetLineColor(Color[i]);
    stat->SetTextColor(Color[i]);
    stat->SetBorderSize(0);
    stat->SetMargin(0.05);
    stat->SetTextAlign(12);
    stat->Draw();
  }
}

// draw a TH2D histogram
void DrawTH2D(TH2D** Histos,
              std::vector<std::string> legend,
              std::string Style_,
              std::string Xlegend,
              std::string Ylegend,
              double xmin,
              double xmax,
              double ymin,
              double ymax) {
  int N = legend.size();

  for (int i = 0; i < N; i++) {
    if (!Histos[i])
      continue;
    Histos[i]->SetTitle("");
    Histos[i]->SetStats(kFALSE);
    Histos[i]->GetXaxis()->SetTitle(Xlegend.c_str());
    Histos[i]->GetYaxis()->SetTitle(Ylegend.c_str());
    Histos[i]->GetYaxis()->SetTitleOffset(1.60);
    if (xmin != xmax)
      Histos[i]->SetAxisRange(xmin, xmax, "X");
    if (ymin != ymax)
      Histos[i]->SetAxisRange(ymin, ymax, "Y");
    Histos[i]->SetMarkerStyle(Marker[i]);
    Histos[i]->SetMarkerColor(Color[i]);
    Histos[i]->SetMarkerSize(0.3);
  }

  char Buffer[256];
  Histos[0]->Draw(Style_.c_str());
  for (int i = 1; i < N; i++) {
    sprintf(Buffer, "%s same", Style_.c_str());
    Histos[i]->Draw(Buffer);
  }
}

// Draw a list of TH1 and superimposed them
void DrawSuperposedHistos(TH1** Histos,
                          std::vector<std::string> legend,
                          std::string Style_,
                          std::string Xlegend,
                          std::string Ylegend,
                          double xmin,
                          double xmax,
                          double ymin,
                          double ymax,
                          bool Normalize = false,
                          bool same = false,
                          bool lastBinOverflow = false,
                          bool firstBinOverflow = false) {
  int N = legend.size();

  double HistoMax = -1;
  int HistoHeighest = -1;

  for (int i = 0; i < N; i++) {
    if (!Histos[i])
      continue;
    if (Normalize && Histos[i]->Integral() != 0)
      Histos[i]->Scale(1.0 / Histos[i]->Integral());
    Histos[i]->SetTitle("");
    Histos[i]->SetStats(kFALSE);
    Histos[i]->GetXaxis()->SetTitle(Xlegend.c_str());
    Histos[i]->GetYaxis()->SetTitle(Ylegend.c_str());
    Histos[i]->GetXaxis()->SetTitleOffset(1.1);
    Histos[i]->GetYaxis()->SetTitleOffset(1.40);
    Histos[i]->GetXaxis()->SetNdivisions(505);
    Histos[i]->GetYaxis()->SetNdivisions(505);
    Histos[i]->GetXaxis()->SetTitleSize(0.05);
    if (xmin != xmax)
      Histos[i]->SetAxisRange(xmin, xmax, "X");
    if (ymin != ymax)
      Histos[i]->SetAxisRange(ymin, ymax, "Y");
    if (ymin == ymax && ymin < 0)
      Histos[i]->SetMaximum(Histos[i]->GetMaximum() * fabs(ymax));
    Histos[i]->SetFillColor(0);
    Histos[i]->SetMarkerStyle(Marker[i]);
    Histos[i]->SetMarkerColor(Color[i]);
    Histos[i]->SetMarkerSize(1.5);
    Histos[i]->SetLineColor(Color[i]);
    Histos[i]->SetLineWidth(2);
    if (lastBinOverflow) {
      if (xmin != xmax) {
        int lastBin = Histos[i]->GetXaxis()->FindBin(xmax);
        double sum = 0;
        double error = 0;
        for (int b = lastBin; b < Histos[i]->GetNbinsX() + 2; b++) {
          sum += Histos[i]->GetBinContent(b);
          error += Histos[i]->GetBinError(b) * Histos[i]->GetBinError(b);
        }
        Histos[i]->SetBinContent(lastBin, sum);
        Histos[i]->SetBinError(lastBin, sqrt(error));
      } else {
        Histos[i]->SetBinContent(
            Histos[i]->GetNbinsX(),
            Histos[i]->GetBinContent(Histos[i]->GetNbinsX()) + Histos[i]->GetBinContent(Histos[i]->GetNbinsX() + 1));
        double error = sqrt(pow(Histos[i]->GetBinError(Histos[i]->GetNbinsX()), 2) +
                            pow(Histos[i]->GetBinError(Histos[i]->GetNbinsX() + 1), 2));
        Histos[i]->SetBinError(Histos[i]->GetNbinsX(), error);
      }
    }
    if (firstBinOverflow) {
      if (xmin != xmax) {
        int firstBin = Histos[i]->GetXaxis()->FindBin(xmin);
        double sum = 0;
        double error = 0;
        for (int b = 0; b < firstBin; b++) {
          sum += Histos[i]->GetBinContent(b);
          error += Histos[i]->GetBinError(b) * Histos[i]->GetBinError(b);
        }
        Histos[i]->SetBinContent(firstBin, sum);
        Histos[i]->SetBinError(firstBin, sqrt(error));
      } else {
        Histos[i]->SetBinContent(1, Histos[i]->GetBinContent(1) + Histos[i]->GetBinContent(0));
        double error = sqrt(pow(Histos[i]->GetBinError(1), 2) + pow(Histos[i]->GetBinError(0), 2));
        Histos[i]->SetBinError(1, error);
      }
    }
    if (Style_ == "DataMC" && i == 0) {
      Histos[i]->SetFillColor(0);
      Histos[i]->SetMarkerStyle(20);
      Histos[i]->SetMarkerColor(1);
      Histos[i]->SetMarkerSize(1);
      Histos[i]->SetLineColor(1);
      Histos[i]->SetLineWidth(2);
    }

    if (Histos[i]->GetMaximum() >= HistoMax) {
      HistoMax = Histos[i]->GetMaximum();
      HistoHeighest = i;
    }
  }

  char Buffer[256];
  if (Style_ == "DataMC") {
    if (HistoHeighest == 0) {
      Histos[HistoHeighest]->Draw("E1");
    } else {
      Histos[HistoHeighest]->Draw("HIST");
    }
    for (int i = 0; i < N; i++) {
      if (i == HistoHeighest)
        continue;
      if (i == 0) {
        Histos[i]->Draw("same E1");
      } else {
        Histos[i]->Draw("same");
      }
    }
  } else {
    if (same) {
      sprintf(Buffer, "same %s", Style_.c_str());
      Histos[HistoHeighest]->Draw(Buffer);
    } else
      Histos[HistoHeighest]->Draw(Style_.c_str());
    for (int i = 0; i < N; i++) {
      if (i == HistoHeighest)
        continue;
      if (Style_ != "") {
        sprintf(Buffer, "same %s", Style_.c_str());
      } else {
        sprintf(Buffer, "same");
      }
      Histos[i]->Draw(Buffer);
    }
  }
}

// automatically determined what is the best axis ranges for a TH2D
void Smart_SetAxisRange(TH2D* histo) {
  double Min = 1E50;
  double Max = 1E-50;
  for (int x = 1; x <= histo->GetNbinsX(); x++) {
    for (int y = 1; y <= histo->GetNbinsY(); y++) {
      double c = histo->GetBinContent(x, y);
      if (c < Min && c > 0)
        Min = c;
      if (c > Max)
        Max = c;
    }
  }
  if (Max / Min < 10) {
    Max *= 5.0;
    Min /= 5.0;
  } else if (Max / Min < 100) {
    Max *= 10.0;
    Min /= 10.0;
  }
  histo->SetAxisRange(Min, Max, "Z");
}

// return a TCUTG corresponding to the uncertainty on a xsection
TCutG* GetErrorBand(
    string name, int N, const double* Mass, const double* Low, const double* High, double MinLow, double MaxHigh) {
  TCutG* cutg = new TCutG(name.c_str(), 2 * N);
  cutg->SetFillColor(kGreen - 7);
  double I = 0;
  for (int i = 0; i < N; i++) {
    if (High[i] < MinLow) {
      continue;
    }

    double Min = std::max(Low[i], MinLow);
    cutg->SetPoint(I, Mass[i], Min);
    I++;
  }
  for (int i = 0; i < N; i++) {
    if (Low[N - 1 - i] > MaxHigh) {
      continue;
    }

    double Max = std::min(High[N - 1 - i], MaxHigh);
    cutg->SetPoint(I, Mass[N - 1 - i], Max);
    I++;
  }
  cutg->Set(I);
  return cutg;
}

std::string toLatexRounded(double value, double error = -1, double systError = -1) {
  using namespace std;
  if (value == 0.0 && error == 0.0)
    return string("");

  double power = floor(log10(value));
  if (power <= -3) {
    power = power + 3;
  } else if (power >= 2) {
    power = power - 2;
  } else {
    power = 0;
  }

  value = value / pow(10, power);
  if (error >= 0)
    error = error / pow(10, power);
  if (systError >= 0)
    systError = systError / pow(10, power);
  int ValueFloating;
  if (systError < 0) {
    ValueFloating = 1 + std::max(-1 * log10(error), 0.0);
  } else {
    ValueFloating = 1 + std::max(-1 * log10(systError), std::max(-1 * log10(error), 0.0));
  }
  int ErrorFloating = ValueFloating;

  char tmpchar[255];
  if (power != 0) {
    if (systError < 0) {
      if (error < 0) {
        sprintf(tmpchar, "$(%.*f)\\times 10^{%g}$", ValueFloating, value, power);
      } else {
        sprintf(tmpchar, "$(%.*f\\pm%.*f)\\times 10^{%g}$", ValueFloating, value, ErrorFloating, error, power);
      }
    } else {
      sprintf(tmpchar,
              "$(%.*f\\pm%.*f\\pm%.*f)\\times 10^{%g}$",
              ValueFloating,
              value,
              ErrorFloating,
              error,
              ErrorFloating,
              systError,
              power);
    }

  } else {
    if (systError < 0) {
      if (error < 0) {
        sprintf(tmpchar, "$%.*f$", ValueFloating, value);
      } else {
        sprintf(tmpchar, "$%.*f\\pm%.*f$", ValueFloating, value, ErrorFloating, error);
      }
    } else {
      sprintf(tmpchar, "$%.*f\\pm%.*f\\pm%.*f$", ValueFloating, value, ErrorFloating, error, ErrorFloating, systError);
    }
  }
  return string(tmpchar);
}

//https://twiki.cern.ch/twiki/bin/view/CMS/PoissonErrorBars
#include "Math/QuantFuncMathCore.h"
#include "TMath.h"
TGraphAsymmErrors* getGarwoodErrorBars(TH1D* h1) {
  const double alpha = 1 - 0.6827;
  TGraphAsymmErrors* g = new TGraphAsymmErrors(h1->GetNbinsX());
  g->SetMarkerSize(h1->GetMarkerSize());
  g->SetMarkerStyle(h1->GetMarkerStyle());
  g->SetMarkerColor(h1->GetMarkerColor());
  g->SetLineWidth(h1->GetLineWidth());
  g->SetLineStyle(h1->GetLineStyle());
  g->SetLineColor(h1->GetLineColor());

  int lastPoint = -1;
  for (int i = 1; i <= h1->GetNbinsX(); ++i) {
    if (h1->GetBinContent(i) > 0)
      lastPoint = i;
  }

  int INDEX = 0;
  for (int i = 1; i <= lastPoint; ++i) {
    int N = h1->GetBinContent(i);
    double L = (N == 0) ? 0 : (ROOT::Math::gamma_quantile(alpha / 2, N, 1.));
    double U = ROOT::Math::gamma_quantile_c(alpha / 2, N + 1, 1);
    g->SetPoint(INDEX, h1->GetBinCenter(i), N);
    g->SetPointEYlow(INDEX, N - L);
    g->SetPointEYhigh(INDEX, U - N);
    INDEX++;
  }
  g->Set(INDEX);
  return g;
}

TH1D* getGarwoodErrorBarsNewROOT(TH1D* h1) {
  TH1D* grassCorrector = (TH1D*)h1->Clone((string(h1->GetName()) + "SomeHist").c_str());
  grassCorrector->SetBinErrorOption(TH1::kPoisson);
  int lastPoint = -1;
  for (int i = 1; i <= grassCorrector->GetNbinsX(); i++) {
    if (grassCorrector->GetBinContent(i) > 0)
      lastPoint = i;
  }
  for (int i = 1; i <= lastPoint; i++) {
    if (grassCorrector->GetBinContent(i) == 0)
      grassCorrector->SetBinError(i, 0.5);
    else
      grassCorrector->SetBinError(i, 0.0);
  }
  return grassCorrector;
}

#endif
