# Usage
# python3 compareTheProfiles.py

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
BinNumber = sys.argv[2]

#StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQ_newMethod_v2.root")
#StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wNoProbQ_NewMethod_v1.root")
#StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQNoL1_newMethod_v1.root")
#StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQNoL1NoProbXY_newMethod_v1.root")
#StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_QCDwPt470To600_woProbQ_CodeV2_v1.root")
#StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_QCDwPt50To80_woProbQ_CodeV2_v1.root")
#StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_CodeV2_v2.root")
StandardAnalysisInFile = ROOT.TFile.Open(fileName)


Mass = StandardAnalysisInFile.Get("/analyzer/BaseName/Mass")

#StandardAnalysisWpredInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQ_newMethod_v2.root")
#StandardAnalysisWpredInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wNoProbQ_NewMethod_v1.root")
#StandardAnalysisWpredInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQNoL1_newMethod_v1.root")
#StandardAnalysisWpredInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQNoL1NoProbXY_newMethod_v1.root")
#StandardAnalysisWpredInFile = ROOT.TFile.Open("crab_Analysis_QCDwPt470To600_woProbQ_CodeV2_v1.root")
#StandardAnalysisWpredInFile = ROOT.TFile.Open("crab_Analysis_QCDwPt50To80_woProbQ_CodeV2_v1.root")
#StandardAnalysisWpredInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_CodeV2_v2.root")
StandardAnalysisWpredInFile = ROOT.TFile.Open(fileName)


#Mass_wPred= StandardAnalysisWpredInFile.Get("/analyzer/BaseName/Pred_Mass") # Same as the MassComb
#Mass_wPred= StandardAnalysisWpredInFile.Get("/analyzer/BaseName/Pred_MassComb") # Combination with TOF
Mass_wPred = StandardAnalysisWpredInFile.Get("/analyzer/BaseName/Pred_Mass_CB") # Fill(MOM*sqrt(ICK),pe,CC*BB/AA);
#Mass_wPred= StandardAnalysisWpredInFile.Get("/analyzer/BaseName/Pred_Mass_CB_Flip") # no search result

ProjBin = int(BinNumber)
# https://cds.cern.ch/record/2205281/files/EXO-16-036-pas.pdf
# had cuts on pt>65 and I_as > 0.3 which corresponds to bins 28,35 --> it's just 28
# bins 25,35 mean  pt>65 and I_as > 0.175 --> maybe this is 25,-1 -->  just 25
# bins 18,35 mean no cut on I_as but a cut on the pt
# bins 3,3 mean  pt>60 and I_as > 0.05

massBins = [0.,50.,100.,150.,200.,250.,300.,350.,400.,500.,600.,700.,800.,900.,1000.,4000.]
#massBins = [0.,50.,100.,200.,300.,400.,500.,600.,700.,800.,900.,1000.,4000.]

massBinsArray = np.array(massBins)

Mass_projY_NotRebinned = Mass.ProjectionY("Mass_projY",ProjBin,ProjBin,"e")
Mass_projY = ROOT.TH1F("Mass_projY" , "Mass_projY" , len(massBinsArray)-1, massBinsArray)

Mass_wPred_projY_NotRebinned = Mass_wPred.ProjectionY("Mass_wPred_projY",ProjBin,ProjBin,"e")
Mass_wPred_projY = ROOT.TH1F("Mass_wPred_projY" , "Mass_wPred_projY" , len(massBinsArray)-1, massBinsArray)

print("Mass_projY_NotRebinned.Integral(): ",Mass_projY_NotRebinned.Integral())
print("Mass_wPred_projY_NotRebinned.Integral(): ",Mass_wPred_projY_NotRebinned.Integral())

KSvalue = Mass_projY_NotRebinned.KolmogorovTest(Mass_wPred_projY_NotRebinned,"XD")
print("KS-test: "+str(KSvalue))

for i, bin in enumerate(massBinsArray):
  Mass_projYCont = Mass_projY_NotRebinned.GetBinContent(i+1)
  Mass_projY.SetBinContent(i+1,Mass_projYCont)
  Mass_projYCont_err = Mass_projY_NotRebinned.GetBinError(i+1)
#  Mass_projYCont_err = Mass.GetBinErrorLow(ProjBin,i+1)
  Mass_projY.SetBinError(i+1,Mass_projYCont_err)
 
  Mass_wPred_projYCont = Mass_wPred_projY_NotRebinned.GetBinContent(i+1)
  Mass_wPred_projY.SetBinContent(i+1,Mass_wPred_projYCont)
  Mass_wPred_projYCont_err = Mass_wPred_projY_NotRebinned.GetBinError(i+1)
#  Mass_wPred_projYCont_err = Mass_wPred.GetBinErrorLow(ProjBin,i+1)
  Mass_wPred_projY.SetBinError(i+1,Mass_wPred_projYCont_err)
  
print("----------------------------------------------")
KSvalue2 = Mass_projY.KolmogorovTest(Mass_wPred_projY,"XD")
print("KS-test after rebinning: "+str(KSvalue2))

Mass_projY.SetMarkerColor(1)
Mass_projY.SetLineColor(1)
Mass_projY.SetMarkerStyle(20)
Mass_projY.SetTitle("")
Mass_projY.GetXaxis().SetTitleSize(0.05)
Mass_projY.GetXaxis().SetTitleOffset(1)
Mass_projY.GetXaxis().SetTitle("Mass [GeV]")
Mass_projY.GetYaxis().SetTitle("Tracks/bin")
Mass_projY.GetYaxis().SetTitleSize(0.05)
Mass_projY.GetYaxis().SetTitleOffset(1)
Mass_projY.SetStats(0)
Mass_projY.GetYaxis().SetRangeUser(0.1,10000)


Mass_wPred_projY.SetMarkerColor(2)
Mass_wPred_projY.SetLineColor(2)
Mass_wPred_projY.SetMarkerStyle(20)
Mass_wPred_projY.SetTitle("")
Mass_wPred_projY.GetXaxis().SetTitleSize(0.05)
Mass_wPred_projY.GetXaxis().SetTitleOffset(1)
Mass_wPred_projY.GetXaxis().SetTitle("Mass [GeV]")
Mass_wPred_projY.GetYaxis().SetTitle("Tracks/bin")
Mass_wPred_projY.GetYaxis().SetTitleSize(0.05)
Mass_wPred_projY.GetYaxis().SetTitleOffset(1)
Mass_wPred_projY.SetStats(0)


print("Mass_projY.Integral(): ",Mass_projY.Integral())
print("Mass_wPred_projY.Integral(): ",Mass_wPred_projY.Integral())
    
legMass =  ROOT.TLegend(.45,.75,.80,.9,"","brNDC")
legMass.SetTextFont(42)
legMass.SetTextSize(0.035)
legMass.SetBorderSize(1);
legMass.SetLineColor(1);
legMass.SetLineStyle(1);
legMass.SetLineWidth(1);
legMass.SetFillColor(0);
legMass.SetFillStyle(1001);
legMass.AddEntry(Mass_wPred_projY,"Prediction","LP")
legMass.AddEntry(Mass_projY,"Observation","LP")

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

tex4 = ROOT.TLatex(0.7,0.93,"K-S test v2: "+str(round(KSvalue2,4)));
tex4.SetNDC();
tex4.SetTextFont(52);
tex4.SetTextSize(0.0485);
tex4.SetLineWidth(2);





cMass_projY = ROOT.TCanvas('cMass_projY', 'cMass_projY',800,800)

cMass_projY.SetLogy()

rp = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY)

rp.SetH1DrawOpt("P");
rp.SetH2DrawOpt("P");
    
rp.Draw()
#rp.GetUpperPad().BuildLegend()
rp.SetLeftMargin(0.13);
rp.SetRightMargin(0.05);
rp.SetUpTopMargin(0.1);
rp.SetLowTopMargin(0.02);
rp.SetLowBottomMargin(0.35);
    
rp.GetLowerRefGraph().SetMinimum(0.01);
rp.GetLowerRefGraph().SetMaximum(2);
#rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
#rp.GetLowerRefGraph().SetLineColor(0) #0
rp.GetLowerRefGraph().SetMarkerStyle(20)
rp.GetLowerRefGraph().SetMarkerSize(1);
rp.GetLowYaxis().SetNdivisions(505);
rp.GetLowerRefYaxis().SetTitle("Ratio");
rp.GetLowerRefYaxis().SetTitleSize(0.05);
rp.GetLowerRefYaxis().SetTitleOffset(1);
rp.GetLowerRefYaxis().SetLabelSize(0.035);
    
    
rp.GetLowerRefXaxis().SetTitleSize(0.05);
rp.GetLowerRefXaxis().SetTitleOffset(0.8);
rp.GetLowerRefXaxis().SetLabelSize(0.035);
cMass_projY.Modified()
cMass_projY.Update()
#Mass_projY.Draw()
#Mass_wPred_projY.Draw("SAME")
#rp.Draw("X")

rp.GetUpperPad().cd();
legMass.Draw("SAME")
tex2.Draw("SAME")
tex3.Draw("SAME")
tex4.Draw("SAME")

name = fileName[0:-5] + "_Bin" + str(ProjBin)+ "/cMass_new.png"
cMass_projY.SaveAs(name)
