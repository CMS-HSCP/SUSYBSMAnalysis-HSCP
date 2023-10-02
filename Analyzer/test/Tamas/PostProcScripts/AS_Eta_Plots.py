# Usage
# python3 compareTheProfiles.py

import ROOT
#import tdrstyle

#ROOT.gROOT.SetBatch(True)
#tdrstyle.setTDRStyle()

StandardAnalysisInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQNoL1NoProbXY_newMethod_v2.root")

AS_Eta_RegionA = StandardAnalysisInFile.Get("/analyzer/BaseName/AS_Eta_RegionA")
AS_Eta_RegionB = StandardAnalysisInFile.Get("/analyzer/BaseName/AS_Eta_RegionB")
AS_Eta_RegionC = StandardAnalysisInFile.Get("/analyzer/BaseName/AS_Eta_RegionC")
AS_Eta_RegionD = StandardAnalysisInFile.Get("/analyzer/BaseName/AS_Eta_RegionD")

StandardAnalysisWprobQInFile = ROOT.TFile.Open("crab_Analysis_SingleMuon_UL2017CEra_wProbQNoL1NoProbXY_newMethod_v2.root")
AS_Eta_RegionA_wProbQ = StandardAnalysisWprobQInFile.Get("/analyzer/BaseName/AS_Eta_RegionA")
AS_Eta_RegionB_wProbQ = StandardAnalysisWprobQInFile.Get("/analyzer/BaseName/AS_Eta_RegionB")
AS_Eta_RegionC_wProbQ = StandardAnalysisWprobQInFile.Get("/analyzer/BaseName/AS_Eta_RegionC")
AS_Eta_RegionD_wProbQ = StandardAnalysisWprobQInFile.Get("/analyzer/BaseName/AS_Eta_RegionD")

# https://cds.cern.ch/record/2205281/files/EXO-16-036-pas.pdf
# had cuts on pt>65 and I_as > 0.3 which corresponds to bins 28,35
# bins 25,35 mean  pt>65 and I_as > 0.175
# bins 18,35 mean no cut on I_as but a cut on the pt
AS_Eta_RegionA_projY = AS_Eta_RegionA.ProjectionY("AS_Eta_RegionA_projY",25,35,"e")
AS_Eta_RegionB_projY = AS_Eta_RegionB.ProjectionY("AS_Eta_RegionB_projY",25,35,"e")
AS_Eta_RegionC_projY = AS_Eta_RegionC.ProjectionY("AS_Eta_RegionC_projY",25,35,"e")
AS_Eta_RegionD_projY = AS_Eta_RegionD.ProjectionY("AS_Eta_RegionD_projY",25,35,"e")

AS_Eta_RegionA_wProbQ_projY = AS_Eta_RegionA_wProbQ.ProjectionY("AS_Eta_RegionA_wProbQ_projY",25,35,"e")
AS_Eta_RegionB_wProbQ_projY = AS_Eta_RegionB_wProbQ.ProjectionY("AS_Eta_RegionB_wProbQ_projY",25,35,"e")
AS_Eta_RegionC_wProbQ_projY = AS_Eta_RegionC_wProbQ.ProjectionY("AS_Eta_RegionC_wProbQ_projY",25,35,"e")
AS_Eta_RegionD_wProbQ_projY = AS_Eta_RegionD_wProbQ.ProjectionY("AS_Eta_RegionD_wProbQ_projY",25,35,"e")


AS_Eta_RegionA_projY.SetMarkerColor(2)
AS_Eta_RegionB_projY.SetMarkerColor(2)
AS_Eta_RegionC_projY.SetMarkerColor(2)
AS_Eta_RegionD_projY.SetMarkerColor(2)

AS_Eta_RegionA_projY.SetLineColor(2)
AS_Eta_RegionB_projY.SetLineColor(2)
AS_Eta_RegionC_projY.SetLineColor(2)
AS_Eta_RegionD_projY.SetLineColor(2)

AS_Eta_RegionA_wProbQ_projY.SetMarkerColor(3)
AS_Eta_RegionB_wProbQ_projY.SetMarkerColor(3)
AS_Eta_RegionC_wProbQ_projY.SetMarkerColor(3)
AS_Eta_RegionD_wProbQ_projY.SetMarkerColor(3)

AS_Eta_RegionA_wProbQ_projY.SetLineColor(3)
AS_Eta_RegionB_wProbQ_projY.SetLineColor(3)
AS_Eta_RegionC_wProbQ_projY.SetLineColor(3)
AS_Eta_RegionD_wProbQ_projY.SetLineColor(3)



print("AS_Eta_RegionA_projY.Integral(): ",AS_Eta_RegionA_projY.Integral())
print("AS_Eta_RegionA_wProbQ_projY.Integral(): ",AS_Eta_RegionA_wProbQ_projY.Integral())

print("AS_Eta_RegionD_projY.Integral(): ",AS_Eta_RegionD_projY.Integral())
print("AS_Eta_RegionD_wProbQ_projY.Integral(): ",AS_Eta_RegionD_wProbQ_projY.Integral())

#    GenHSCP_probQ_wNoCuts.SetMarkerColor(2)
#    GenHSCP_probQ_wNoCuts.SetLineColor(2)
#    GenHSCP_probQ_wNoCuts.SetMarkerStyle(20)
#
#    legProbQ =  ROOT.TLegend(.25,.75,.60,.9)
#    legProbQ.SetTextFont(42)
#    legProbQ.SetTextSize(0.035)
#    legProbQ.AddEntry(GenHSCP_probQ_wNoCuts,"HSCP gen truth #tracks:","LP")
#
#
#    tex2 = ROOT.TLatex(0.18,0.96,"CMS");
#    #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
#    tex2.SetNDC();
#    tex2.SetTextFont(61);
#    tex2.SetTextSize(0.0375);
#    tex2.SetLineWidth(2);
#
#
##    tex3 = ROOT.TLatex(0.27,0.96,"Simulation"); # for square plots
#    #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
#    tex3 = ROOT.TLatex(0.24,0.96,"Simulation");
#    tex3.SetNDC();
#    tex3.SetTextFont(52);
#    tex3.SetTextSize(0.0285);
#    tex3.SetLineWidth(2);
#
cAS_Eta_RegionA_projY = ROOT.TCanvas('cAS_Eta_RegionA_projY', 'cAS_Eta_RegionA_projY',800,600)
AS_Eta_RegionA_projY.Draw()
AS_Eta_RegionA_wProbQ_projY.Draw("SAME")
#
#    tex2.Draw("SAME")
#    tex3.Draw("SAME")
#
cAS_Eta_RegionA_projY.SaveAs("cAS_Eta_RegionA_projY.png")

cAS_Eta_RegionB_projY = ROOT.TCanvas('cAS_Eta_RegionB_projY', 'cAS_Eta_RegionB_projY',800,600)
AS_Eta_RegionB_projY.Draw()
AS_Eta_RegionB_wProbQ_projY.Draw("SAME")
#
#    tex2.Draw("SAME")
#    tex3.Draw("SAME")
#
cAS_Eta_RegionB_projY.SaveAs("cAS_Eta_RegionB_projY.png")

cAS_Eta_RegionC_projY = ROOT.TCanvas('cAS_Eta_RegionC_projY', 'cAS_Eta_RegionC_projY',800,600)
AS_Eta_RegionC_projY.Draw()
AS_Eta_RegionC_wProbQ_projY.Draw("SAME")
#
#    tex2.Draw("SAME")
#    tex3.Draw("SAME")
#
cAS_Eta_RegionC_projY.SaveAs("cAS_Eta_RegionC_projY.png")



cAS_Eta_RegionD_projY = ROOT.TCanvas('cAS_Eta_RegionD_projY', 'cAS_Eta_RegionD_projY',800,600)
AS_Eta_RegionD_projY.Draw()
AS_Eta_RegionD_wProbQ_projY.Draw("SAME")
#
#    tex2.Draw("SAME")
#    tex3.Draw("SAME")
#
cAS_Eta_RegionD_projY.SaveAs("cAS_Eta_RegionD_projY.png")
#
