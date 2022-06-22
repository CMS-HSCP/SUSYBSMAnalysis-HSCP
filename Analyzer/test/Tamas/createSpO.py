import ROOT, sys, os, time, re
import numpy as np
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog fileNameSignal.root")
(opt,args) = parser.parse_args()

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.SetBatch(True)

fileNameSignal = sys.argv[1]

f1 = ROOT.TFile.Open(fileNameSignal)
f2 = ROOT.TFile.Open("crab_Analysis_2018_TTToSemiLeptonic_wProbQ_CodeV17p9_v1.root")
f3 = ROOT.TFile.Open("crab_Analysis_2018_TTTo2L2Nu_wProbQ_CodeV17p9_v1.root")
f4 = ROOT.TFile.Open("crab_Analysis_2018_TTToHadronic_wProbQ_CodeV17p9_v1.root")
f5 = ROOT.TFile.Open("crab_Analysis_2018_QCDwPt1000_wProbQ_wIhCut_CodeV17p9_v1.root")

obj1 = f1.Get("analyzer/BaseName/BS_MiniRelIsoAll")
#obj1 = f1.Get("analyzer/BaseName/BS_EIsol")
obj2 = f2.Get("analyzer/BaseName/BS_MiniRelIsoAll")
obj3 = f3.Get("analyzer/BaseName/BS_MiniRelIsoAll")
obj4 = f4.Get("analyzer/BaseName/BS_MiniRelIsoAll")
obj5 = f5.Get("analyzer/BaseName/BS_MiniRelIsoAll")

can = ROOT.TCanvas("CutEffs")
SignalOverNoise = ROOT.TH1F("SignalOverNoise","SignalOverNoise",50,0.0,1.0)
#SignalOverNoise = ROOT.TH1F("SignalOverNoise","SignalOverNoise",25,0.0,1.5)
SignalOverNoiseTTto2L2N = ROOT.TH1F("SignalOverNoiseTTto2L2N","SignalOverNoiseTTto2L2N",50,0.0,1.0)
SignalOverNoiseTTtHadronic = ROOT.TH1F("SignalOverNoiseTTto2L2N","SignalOverNoiseTTto2L2N",50,0.0,1.0)
SignalOverNoiseQCD = ROOT.TH1F("SignalOverNoiseQCD","SignalOverNoiseQCD",50,0.0,1.0)

#for x in np.arange(0,1.5,0.05):
for x in np.arange(0,1.0,0.02):
  Signal = obj1.Integral(obj1.GetXaxis().FindBin(0),obj1.GetXaxis().FindBin(x))
  BackgroundTTtoSemileptonic = obj2.Integral(obj1.GetXaxis().FindBin(0),obj2.GetXaxis().FindBin(x))
  BackgroundTTto2L2N = obj3.Integral(obj1.GetXaxis().FindBin(0),obj3.GetXaxis().FindBin(x))
  BackgroundTTtoHadronic = obj4.Integral(obj1.GetXaxis().FindBin(0),obj4.GetXaxis().FindBin(x))
  BackgroundQCD = obj5.Integral(obj1.GetXaxis().FindBin(0),obj5.GetXaxis().FindBin(x))
#  print("For ",x," the S = ",Signal," and B = ",Background, " and the S/sqrt(S+B)",Signal/np.sqrt(Signal+Background))
  SignalOverNoise.SetBinContent(obj1.GetXaxis().FindBin(x),(Signal/np.sqrt(Signal+BackgroundTTtoSemileptonic)))
  SignalOverNoiseTTto2L2N.SetBinContent(obj1.GetXaxis().FindBin(x),(Signal/np.sqrt(Signal+BackgroundTTto2L2N)))
  SignalOverNoiseTTtHadronic.SetBinContent(obj1.GetXaxis().FindBin(x),(Signal/np.sqrt(Signal+BackgroundTTtoHadronic)))
  SignalOverNoiseQCD.SetBinContent(obj1.GetXaxis().FindBin(x),(Signal/np.sqrt(Signal+BackgroundQCD)))
  # s / sqrt(s+b)

max1 = np.maximum(SignalOverNoise.GetMaximum(),SignalOverNoiseTTto2L2N.GetMaximum())
max2 = np.maximum(SignalOverNoiseTTtHadronic.GetMaximum(),SignalOverNoiseQCD.GetMaximum())
max = np.maximum(max1,max2)
SignalOverNoise.SetMaximum(max*1.4)
SignalOverNoise.SetLineColor(2)
SignalOverNoise.SetMarkerColor(2)
SignalOverNoise.SetMarkerStyle(20)
SignalOverNoise.GetYaxis().SetTitle("S/#sqrt{S+B}")
SignalOverNoise.GetYaxis().SetTitleOffset(1.1)
#SignalOverNoise.GetXaxis().SetTitle("(E_{HCAL} + E_{ECAL})/p")
SignalOverNoise.GetXaxis().SetTitle("BS_MiniRelIsoAll")
SignalOverNoise.SetStats(0)

SignalOverNoiseTTto2L2N.SetMarkerStyle(20)
SignalOverNoiseTTtHadronic.SetMarkerStyle(20)
SignalOverNoiseQCD.SetMarkerStyle(20)

SignalOverNoiseTTto2L2N.SetMarkerColor(3)
SignalOverNoiseTTtHadronic.SetMarkerColor(4)
SignalOverNoiseQCD.SetMarkerColor(5)

legend =  ROOT.TLegend(.11,.75,.35,.89,"","brNDC")
legend.SetTextFont(42)
#legend.SetTextSize(0.035)
legend.SetTextSize(0.025)
legend.SetBorderSize(1);
legend.SetLineColor(0);
legend.SetLineStyle(1);
legend.SetLineWidth(1);
legend.SetFillColor(0);
legend.SetFillStyle(1001);
legend.AddEntry(SignalOverNoise,"S="+fileNameSignal[14:-5]+" B=TTtoSemiLeptonic","LP")
legend.AddEntry(SignalOverNoiseTTto2L2N,"S="+fileNameSignal[14:-5]+" B=TTto2L2N","LP")
legend.AddEntry(SignalOverNoiseTTtHadronic,"S="+fileNameSignal[14:-5]+" B=TTtoHadronic","LP")
legend.AddEntry(SignalOverNoiseQCD,"S="+fileNameSignal[14:-5]+" B=QCD","LP")

SignalOverNoise.Draw("P")
SignalOverNoiseTTto2L2N.Draw("SAMEP")
SignalOverNoiseTTtHadronic.Draw("SAMEP")
SignalOverNoiseQCD.Draw("SAMEP")
legend.Draw("SAME")
can.SaveAs("SignalOverNoise_MiniRelIso_"+fileNameSignal[14:-5]+".png")
