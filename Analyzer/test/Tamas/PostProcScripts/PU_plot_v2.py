import ROOT, sys
import numpy as np
#import tdrstyle

from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog fileName.root BinNumber")
(opt,args) = parser.parse_args()

#ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.065);
ROOT.gStyle.SetPadBottomMargin(0.17);
ROOT.gStyle.SetPadLeftMargin(0.15);

fileName = sys.argv[1]
BinNumber = sys.argv[2] if (len(sys.argv)==1) else 3

StandardAnalysisInFile = ROOT.TFile.Open(fileName)

PUA = StandardAnalysisInFile.Get("/HSCParticleAnalyzer/BaseName/PostS_SR2PASS_Ias_PUA")
PUB = StandardAnalysisInFile.Get("/HSCParticleAnalyzer/BaseName/PostS_SR2PASS_Ias_PUB")
PUC = StandardAnalysisInFile.Get("/HSCParticleAnalyzer/BaseName/PostS_SR2PASS_Ias_PUC")


c = ROOT.TCanvas("c", "c", 800, 800)
c.Divide(1,2)

ratio1 = ROOT.TRatioPlot(PUB, PUA)
ratio2 = ROOT.TRatioPlot(PUB, PUC)



c.cd(1)
PUB.Draw("hist")
PUA.Draw("hist same")
PUB.SetLineColor(ROOT.kRed)
PUA.SetLineColor(ROOT.kBlue)

c.cd(2)
ratio1.Draw()
ratio1.SetLeftMargin(0.13);
ratio1.SetRightMargin(0.05);
ratio1.SetUpTopMargin(0.1);
ratio1.SetLowTopMargin(0.02);
ratio1.SetLowBottomMargin(0.35);
ratio2.Draw("same")

ratio1.GetLowerRefYaxis().SetTitle("#frac{pileup B}{pileup A}")
ratio2.GetLowerRefYaxis().SetTitle("#frac{pileup B}{pileup C}")
ratio1.GetLowerRefXaxis().SetTitle("vertex multiplicity")
ratio2.GetLowerRefXaxis().SetTitle("vertex multiplicity")

ratio1.GetUpperRefYaxis().SetTitle("events")
ratio2.GetUpperRefYaxis().SetTitle("events")
ratio1.GetUpperRefYaxis().SetTitleOffset(1.2)
ratio2.GetUpperRefYaxis().SetTitleOffset(1.2)

ratio1.GetLowerRefGraph().SetMinimum(0)
ratio2.GetLowerRefGraph().SetMinimum(0)

ratio1.GetLowerRefGraph().SetMaximum(2)
ratio2.GetLowerRefGraph().SetMaximum(2)

ratio1.GetUpperRefObject().SetLineColor(ROOT.kBlue)
ratio2.GetUpperRefObject().SetLineColor(ROOT.kGreen)

ratio1.GetLowerRefGraph().SetLineColor(ROOT.kBlue)
ratio2.GetLowerRefGraph().SetLineColor(ROOT.kGreen)

ratio1.GetLowerRefGraph().SetMarkerColor(ROOT.kBlue)
ratio2.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen)

ratio1.GetLowerRefGraph().SetMarkerStyle(20)
ratio2.GetLowerRefGraph().SetMarkerStyle(21)

ratio1.GetLowerRefGraph().SetLineWidth(2)
ratio2.GetLowerRefGraph().SetLineWidth(2)

c.Update()

c.cd(1)
leg = ROOT.TLegend(0.6, 0.7, 0.9, 0.9)
leg.AddEntry(PUB, "pileup B", "l")
leg.AddEntry(PUA, "pileup A", "l")
leg.Draw()

c.SaveAs("PU_plot_v2.png")
