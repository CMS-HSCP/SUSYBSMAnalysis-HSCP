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

PUA = StandardAnalysisInFile.Get("/HSCParticleAnalyzer/BaseName/PostS_SR2FAIL_PV")
PUB = StandardAnalysisInFile.Get("/HSCParticleAnalyzer/BaseName/PostS_SR2PASS_PV")

PUA.Scale(1/PUA.Integral())
PUB.Scale(1/PUB.Integral())

PUA.Rebin(5)
PUB.Rebin(5)

KSvalueBtoA = round(PUB.Chi2Test(PUA,"XD"),3) #KolmogorovTest
KSvalueAtoB = round(PUA.Chi2Test(PUB,"XD"),3) #KolmogorovTest


max1 = np.maximum(PUA.GetMaximum(),PUB.GetMaximum())
max2 = 0.01
max = np.maximum(max1,max2)

PUA.SetMarkerColor(2)
PUA.SetLineColor(2)
PUA.SetMarkerStyle(20)
PUA.GetXaxis().SetTitleSize(0.015)
PUA.GetXaxis().SetTitleOffset(1)
PUA.GetYaxis().SetTitleSize(0.05)
PUA.GetYaxis().SetTitleOffset(1.3)
PUA.GetYaxis().SetTitle("Normalized events / bin") # + PUA.GetYaxis().GetTitle())
PUA.GetXaxis().SetTitle("PV (p_{T} > 200 GeV)")
PUA.SetStats(0)
PUA.SetMaximum(max*10)
PUA.SetMinimum(0.0001)
#PUA.GetYaxis().SetRangeUser(0.1,10000)


PUB.SetMarkerColor(1)
PUB.SetLineColor(1)
PUB.SetMarkerStyle(20)
PUB.SetStats(0)
    
legPUA =  ROOT.TLegend(.45,.75,.80,.85,"","brNDC")
legPUA.SetTextFont(42)
legPUA.SetTextSize(0.035)
legPUA.SetBorderSize(0);
legPUA.SetLineColor(1);
legPUA.SetLineStyle(1);
legPUA.SetLineWidth(1);
legPUA.SetFillColor(0);
legPUA.SetFillStyle(1001);
legPUA.AddEntry(PUA,"PV in FAIL, #chi^{2} test wrt to PASS: "+str(KSvalueBtoA),"LP")
legPUA.AddEntry(PUB,"PV in PASS, #chi^{2} test wrt to FAIL: "+str(KSvalueAtoB),"LP")

tex2 = ROOT.TLatex(0.13,0.94,"CMS");
    #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
tex2.SetNDC();
tex2.SetTextFont(61);
tex2.SetTextSize(0.0675);
tex2.SetLineWidth(2);

#tex3 = ROOT.TLatex(0.27,0.96,"Simulation"); # for square plots
    #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
tex3 = ROOT.TLatex(0.24,0.94,"Internal");
tex3.SetNDC();
tex3.SetTextFont(52);
tex3.SetTextSize(0.0485);
tex3.SetLineWidth(2);

tex4 = ROOT.TLatex(0.7,0.93,"Post Selection");
tex4.SetNDC();
tex4.SetTextFont(52);
tex4.SetTextSize(0.0485);
tex4.SetLineWidth(2);



cPUA = ROOT.TCanvas('cPUA', 'cPUA',800,800)

cPUA.SetLogy()

rp = ROOT.TRatioPlot(PUA,PUB,"diffsig")
rp.SetDrawOption("pois")

rp.SetH1DrawOpt("P");
rp.SetH2DrawOpt("P");

rp.Draw()

#rp.GetUpperPad().BuildLegend()
rp.SetLeftMargin(0.13);
rp.SetRightMargin(0.05);
rp.SetUpTopMargin(0.1);
rp.SetLowTopMargin(0.02);
rp.SetLowBottomMargin(0.35);
    
rp.GetLowerRefGraph().SetMinimum(-5);
rp.GetLowerRefGraph().SetMaximum(5);
rp.GetLowerRefGraph().SetMarkerColor(2)
rp.GetLowerRefGraph().SetLineColor(2) #0
rp.GetLowerRefGraph().SetMarkerStyle(20)
rp.GetLowerRefGraph().SetMarkerSize(1);
rp.GetLowYaxis().SetNdivisions(505);
rp.GetLowerRefYaxis().SetTitle("Pull");
rp.GetLowerRefYaxis().SetTitleSize(0.05);
rp.GetLowerRefYaxis().SetTitleOffset(1);
rp.GetLowerRefYaxis().SetLabelSize(0.035);



    
rp.GetLowerRefXaxis().SetTitleSize(0.04);
rp.GetLowerRefXaxis().SetTitleOffset(1);
rp.GetLowerRefXaxis().SetLabelSize(0.035);
cPUA.Modified()
cPUA.Update()
#PUA.Draw()
#PUB.Draw("SAME")
#rp.Draw("X")

rp.GetUpperPad().cd();
legPUA.Draw("SAME")
tex2.Draw("SAME")
tex3.Draw("SAME")
tex4.Draw("SAME")
rp.GetLowerPad().cd();

#name = fileName[0:-5] + "_Bin" + str(ProjBin)+ "/cPUA_new.png"
#cPUA.SaveAs(name)
cPUA.SaveAs("name2.png")
