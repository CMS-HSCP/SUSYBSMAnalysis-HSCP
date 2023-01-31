import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
from tqdm import tqdm

codeVersion = sys.argv[1]
#just the number, like 18p2

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

BackgroundSamples = [
"crab_Analysis_2018_HSCPgluino_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2600_CodeV"+codeVersion+"_v1.root",
]

fileInArray = []
for sample in BackgroundSamples:
  fileInArray.append(ROOT.TFile.Open(sample,"UPDATE"))
  
  
cstackedSummedBackgroundStringN1 = "cstackedSummedBackgroundStringN1"
canvas = ROOT.TCanvas(cstackedSummedBackgroundStringN1, cstackedSummedBackgroundStringN1, 800,800)

EffForData = ROOT.TH1F("EffForDataN1",";;Integral to (from) max (min) cut value /  Integral to inf",10,0.,10.)

for fileIn in fileInArray:
  if not (fileIn.Get("HSCParticleAnalyzer/BaseName/NumEvents")):
    print("NumEvents not found, exit")
    continue
  
  obj = fileIn.Get("HSCParticleAnalyzer/BaseName/PostPreS_ProbQNoL1")

  firstBinInt = obj.Integral(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(0.1)+1)
  totalInt = obj.Integral()
  efficiency = firstBinInt/float(totalInt)

  
      
      EffForData = ROOT.TH1F("EffForDataN1",";;Integral to (from) max (min) cut value /  Integral to inf",10,0.,10.)
      EffForData.Draw("HISTTEXT00")
      EffForData.SetBinContent(1,N1effForMiniRelTkIso)
      EffForData.SetBinContent(2,N1effForMiniRelTkIsoPUA)
      EffForData.SetBinContent(3,N1effForMiniRelTkIsoPUB)
      EffForData.SetBinContent(4,N1effForMiniRelTkIsoPUC)
      EffForData.SetBinContent(5,N1effForMiniTkIso)
      EffForData.SetBinContent(6,N1effForMiniTkIsoPUA)
      EffForData.SetBinContent(7,N1effForMiniTkIsoPUB)
      EffForData.SetBinContent(8,N1effForMiniTkIsoPUC)
      EffForData.SetBinContent(9,N1effForMiniIso)
      EffForData.SetBinContent(10,N1effForTIsol)

      
      EffForData.SetLineColor(1)
      EffForData.SetMarkerColor(1)
#      EffForData.SetMarkerStyle(20)
      EffForData.SetStats(0)
      EffForData.GetXaxis().SetTitle("")
      EffForData.GetYaxis().SetRangeUser(0.,1.3)
      EffForData.GetXaxis().SetBinLabel(1,"MiniRelTkIso")
      EffForData.GetXaxis().SetBinLabel(2,"MiniRelTkIsoPUA")
      EffForData.GetXaxis().SetBinLabel(3,"MiniRelTkIsoPUB")
      EffForData.GetXaxis().SetBinLabel(4,"MiniRelTkIsoPUC")
      EffForData.GetXaxis().SetBinLabel(5,"MiniTkIso")
      EffForData.GetXaxis().SetBinLabel(6,"MiniTkIsoPUA")
      EffForData.GetXaxis().SetBinLabel(7,"MiniTkIsoPUB")
      EffForData.GetXaxis().SetBinLabel(8,"MiniTkIsoPUC")
      EffForData.GetXaxis().SetBinLabel(9,"MiniRelIso < 0.02")
      EffForData.GetXaxis().SetBinLabel(10,"TIsol < 15 GeV")
      EffForData.GetXaxis().SetTitle("")
      
      legend =  ROOT.TLegend(.55,.80,.80,.89,"","brNDC")
      legend.SetTextFont(42)
      legend.SetTextSize(0.02)
      legend.SetBorderSize(1);
      legend.SetBorderSize(0);
      legend.SetLineColor(1);
      legend.SetLineStyle(1);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(1001);
      legend.AddEntry(EffForData,"DATA (2018C)","LP")
