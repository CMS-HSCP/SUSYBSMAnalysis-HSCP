import ROOT, sys

ROOT.gStyle.SetPaintTextFormat(".2g");

fileName = sys.argv[1]
codeVersion = fileName[fileName.find("CodeV")+5:fileName.find("CodeV")+9]

can = ROOT.TCanvas("newname","newname",800,800)
f = ROOT.TFile.Open(fileName)
HalfSigma_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_Beta")
HalfSigma_UP = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_BetaUpHalfSigma")
HalfSigma_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_BetaDownHalfSigma")
#OneSigma_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_Beta")
OneSigma_UP = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_BetaUpOneSigma")
OneSigma_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_BetaDownOneSigma")
#TwoSigma_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_Beta")
TwoSigma_UP = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_BetaUpTwoSigma")
TwoSigma_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostS_SR2PASS_TriggerMuon50VsBeta_BetaDownTwoSigma")

HalfSigma_Nom_ProfY = HalfSigma_Nom.RebinY(2).ProfileY().ProjectionX()
HalfSigma_UP_ProfY = HalfSigma_UP.RebinY(2).ProfileY().ProjectionX()
HalfSigma_DOWN_ProfY = HalfSigma_DOWN.RebinY(2).ProfileY().ProjectionX()
#OneSigma_Nom_ProfY = OneSigma_Nom.RebinY(2).ProfileY().ProjectionX()
OneSigma_UP_ProfY = OneSigma_UP.RebinY(2).ProfileY().ProjectionX()
OneSigma_DOWN_ProfY = OneSigma_DOWN.RebinY(2).ProfileY().ProjectionX()
#TwoSigma_Nom_ProfY = TwoSigma_Nom.RebinY(2).ProfileY().ProjectionX()
TwoSigma_UP_ProfY = TwoSigma_UP.RebinY(2).ProfileY().ProjectionX()
TwoSigma_DOWN_ProfY = TwoSigma_DOWN.RebinY(2).ProfileY().ProjectionX()

# --------------------------------------------------
HalfSigma_UpRatio = HalfSigma_UP_ProfY.Clone()
HalfSigma_UpRatio.Divide(HalfSigma_Nom_ProfY)

HalfSigma_DownRatio = HalfSigma_DOWN_ProfY.Clone()
HalfSigma_DownRatio.Divide(HalfSigma_Nom_ProfY)

HalfSigma_UpRatio.Draw("COLZTEXT40")
HalfSigma_UpRatio.SetMarkerColor(2)
HalfSigma_UpRatio.SetMarkerStyle(20)
HalfSigma_UpRatio.SetStats(0)
HalfSigma_UpRatio.SetLineColor(2)
HalfSigma_UpRatio.GetYaxis().SetTitle("Ratio")
HalfSigma_DownRatio.Draw("SAMETEXT40")
HalfSigma_DownRatio.SetMarkerStyle(20)
HalfSigma_DownRatio.SetMarkerColor(3)
HalfSigma_DownRatio.SetLineColor(3)
legRatiosSyst_HalfSigma =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_HalfSigma.SetTextFont(42)
legRatiosSyst_HalfSigma.SetTextSize(0.035)
legRatiosSyst_HalfSigma.SetBorderSize(1);
legRatiosSyst_HalfSigma.SetLineColor(0);
legRatiosSyst_HalfSigma.SetLineStyle(1);
legRatiosSyst_HalfSigma.SetLineWidth(1);
legRatiosSyst_HalfSigma.SetFillColor(0);
legRatiosSyst_HalfSigma.SetFillStyle(1001);
legRatiosSyst_HalfSigma.AddEntry(HalfSigma_UpRatio,"Up (0.5#sigma) systematics / nominal","LP")
legRatiosSyst_HalfSigma.AddEntry(HalfSigma_DownRatio,"Down (0.5#sigma) systematics / nominal","LP")
legRatiosSyst_HalfSigma.Draw("SAME")
can.SaveAs(fileName[0:-5] + "_TriggerEff_RatiosSyst_HalfSigma.png")

HalfSigma_Nom_ProfY.Draw("COLZ")
HalfSigma_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
HalfSigma_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
HalfSigma_UP_ProfY.Draw("SAME")
HalfSigma_UP_ProfY.SetLineColor(2)
HalfSigma_DOWN_ProfY.Draw("SAME")
HalfSigma_DOWN_ProfY.SetLineColor(3)
HalfSigma_Nom_ProfY.SetStats(0)
legTriggerEff_HalfSigma =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_HalfSigma.SetTextFont(42)
legTriggerEff_HalfSigma.SetTextSize(0.035)
legTriggerEff_HalfSigma.SetBorderSize(1);
legTriggerEff_HalfSigma.SetLineColor(0);
legTriggerEff_HalfSigma.SetLineStyle(1);
legTriggerEff_HalfSigma.SetLineWidth(1);
legTriggerEff_HalfSigma.SetFillColor(0);
legTriggerEff_HalfSigma.SetFillStyle(1001);
legTriggerEff_HalfSigma.AddEntry(HalfSigma_Nom_ProfY,"Nominal","LP")
legTriggerEff_HalfSigma.AddEntry(HalfSigma_UP_ProfY,"Up (0.5#sigma) systematics","LP")
legTriggerEff_HalfSigma.AddEntry(HalfSigma_DOWN_ProfY,"Down (0.5#sigma) systematics","LP")
legTriggerEff_HalfSigma.Draw("SAME")
can.SaveAs(fileName[0:-5] + "_TriggerEff_Syst_HalfSigma.png")

can.Clear()

# --------------------------------------------------
OneSigma_UpRatio = OneSigma_UP_ProfY.Clone()
OneSigma_UpRatio.Divide(HalfSigma_Nom_ProfY)

OneSigma_DownRatio = OneSigma_DOWN_ProfY.Clone()
OneSigma_DownRatio.Divide(HalfSigma_Nom_ProfY)

OneSigma_UpRatio.Draw("COLZTEXT40")
OneSigma_UpRatio.SetMarkerColor(2)
OneSigma_UpRatio.SetMarkerStyle(20)
OneSigma_UpRatio.SetStats(0)
OneSigma_UpRatio.SetLineColor(2)
OneSigma_UpRatio.GetYaxis().SetTitle("Ratio")
OneSigma_DownRatio.Draw("SAMETEXT40")
OneSigma_DownRatio.SetMarkerStyle(20)
OneSigma_DownRatio.SetMarkerColor(3)
OneSigma_DownRatio.SetLineColor(3)
legRatiosSyst_OneSigma =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_OneSigma.SetTextFont(42)
legRatiosSyst_OneSigma.SetTextSize(0.035)
legRatiosSyst_OneSigma.SetBorderSize(1);
legRatiosSyst_OneSigma.SetLineColor(0);
legRatiosSyst_OneSigma.SetLineStyle(1);
legRatiosSyst_OneSigma.SetLineWidth(1);
legRatiosSyst_OneSigma.SetFillColor(0);
legRatiosSyst_OneSigma.SetFillStyle(1001);
legRatiosSyst_OneSigma.AddEntry(OneSigma_UpRatio,"Up (1#sigma) systematics / nominal","LP")
legRatiosSyst_OneSigma.AddEntry(OneSigma_DownRatio,"Down (1#sigma) systematics / nominal","LP")
legRatiosSyst_OneSigma.Draw("SAME")
can.SaveAs(fileName[0:-5] + "_TriggerEff_RatiosSyst_OneSigma.png")
#print("OneSigma_DownRatio content:")
#for xBin in range(1,OneSigma_DownRatio.GetNbinsX()+1) :
  #"" + str(round(OneSigma_UpRatio.GetXaxis().GetBinUpEdge(xBin),2)) +
#  print(str(round(OneSigma_DownRatio.GetBinContent(xBin),2))+",")

HalfSigma_Nom_ProfY.Draw("COLZ")
HalfSigma_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
HalfSigma_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
OneSigma_UP_ProfY.Draw("SAME")
OneSigma_UP_ProfY.SetLineColor(2)
OneSigma_DOWN_ProfY.Draw("SAME")
OneSigma_DOWN_ProfY.SetLineColor(3)
HalfSigma_Nom_ProfY.SetStats(0)
legTriggerEff_OneSigma =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_OneSigma.SetTextFont(42)
legTriggerEff_OneSigma.SetTextSize(0.035)
legTriggerEff_OneSigma.SetBorderSize(1);
legTriggerEff_OneSigma.SetLineColor(0);
legTriggerEff_OneSigma.SetLineStyle(1);
legTriggerEff_OneSigma.SetLineWidth(1);
legTriggerEff_OneSigma.SetFillColor(0);
legTriggerEff_OneSigma.SetFillStyle(1001);
legTriggerEff_OneSigma.AddEntry(HalfSigma_Nom_ProfY,"Nominal","LP")
legTriggerEff_OneSigma.AddEntry(OneSigma_UP_ProfY,"Up (1#sigma) systematics","LP")
legTriggerEff_OneSigma.AddEntry(OneSigma_DOWN_ProfY,"Down (1#sigma) systematics","LP")
legTriggerEff_OneSigma.Draw("SAME")
can.SaveAs(fileName[0:-5] + "_TriggerEff_Syst_OneSigma.png")

can.Clear()

# --------------------------------------------------
TwoSigma_UpRatio = TwoSigma_UP_ProfY.Clone()
TwoSigma_UpRatio.Divide(HalfSigma_Nom_ProfY)

TwoSigma_DownRatio = TwoSigma_DOWN_ProfY.Clone()
TwoSigma_DownRatio.Divide(HalfSigma_Nom_ProfY)

TwoSigma_UpRatio.Draw("COLZTEXT40")
TwoSigma_UpRatio.SetMarkerColor(2)
TwoSigma_UpRatio.SetMarkerStyle(20)
TwoSigma_UpRatio.SetStats(0)
TwoSigma_UpRatio.SetLineColor(2)
TwoSigma_UpRatio.GetYaxis().SetTitle("Ratio")
TwoSigma_DownRatio.Draw("SAMETEXT40")
TwoSigma_DownRatio.SetMarkerStyle(20)
TwoSigma_DownRatio.SetMarkerColor(3)
TwoSigma_DownRatio.SetLineColor(3)
legRatiosSyst_TwoSigma =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_TwoSigma.SetTextFont(42)
legRatiosSyst_TwoSigma.SetTextSize(0.035)
legRatiosSyst_TwoSigma.SetBorderSize(1);
legRatiosSyst_TwoSigma.SetLineColor(0);
legRatiosSyst_TwoSigma.SetLineStyle(1);
legRatiosSyst_TwoSigma.SetLineWidth(1);
legRatiosSyst_TwoSigma.SetFillColor(0);
legRatiosSyst_TwoSigma.SetFillStyle(1001);
legRatiosSyst_TwoSigma.AddEntry(TwoSigma_UpRatio,"Up (2#sigma) systematics / nominal","LP")
legRatiosSyst_TwoSigma.AddEntry(TwoSigma_DownRatio,"Down (2#sigma) systematics / nominal","LP")
legRatiosSyst_TwoSigma.Draw("SAME")
can.SaveAs(fileName[0:-5] + "_TriggerEff_RatiosSyst_TwoSigma.png")
#print("TwoSigma_DownRatio content:")
#for xBin in range(1,TwoSigma_DownRatio.GetNbinsX()+1) :
  #"" + str(round(TwoSigma_DownRatio.GetXaxis().GetBinUpEdge(xBin),2)) +
#  print(str(round(TwoSigma_DownRatio.GetBinContent(xBin),2))+",")

HalfSigma_Nom_ProfY.Draw("COLZ")
HalfSigma_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
HalfSigma_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
TwoSigma_UP_ProfY.Draw("SAME")
TwoSigma_UP_ProfY.SetLineColor(2)
TwoSigma_DOWN_ProfY.Draw("SAME")
TwoSigma_DOWN_ProfY.SetLineColor(3)
HalfSigma_Nom_ProfY.SetStats(0)
legTriggerEff_TwoSigma =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_TwoSigma.SetTextFont(42)
legTriggerEff_TwoSigma.SetTextSize(0.035)
legTriggerEff_TwoSigma.SetBorderSize(1);
legTriggerEff_TwoSigma.SetLineColor(0);
legTriggerEff_TwoSigma.SetLineStyle(1);
legTriggerEff_TwoSigma.SetLineWidth(1);
legTriggerEff_TwoSigma.SetFillColor(0);
legTriggerEff_TwoSigma.SetFillStyle(1001);
legTriggerEff_TwoSigma.AddEntry(HalfSigma_Nom_ProfY,"Nominal","LP")
legTriggerEff_TwoSigma.AddEntry(TwoSigma_UP_ProfY,"Up (2#sigma) systematics","LP")
legTriggerEff_TwoSigma.AddEntry(TwoSigma_DOWN_ProfY,"Down (2#sigma) systematics","LP")
legTriggerEff_TwoSigma.Draw("SAME")
can.SaveAs(fileName[0:-5] + "_TriggerEff_Syst_TwoSigma.png")

can.Clear()
print("scp -r "+ fileName[0:-5] + "*.png tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/.")
