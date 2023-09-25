import ROOT, sys

ROOT.gStyle.SetPaintTextFormat(".2g");

fileName = sys.argv[1]
codeVersion = fileName[fileName.find("CodeV")+5:fileName.find("CodeV")+9]

can = ROOT.TCanvas("newname","newname",800,800)
f = ROOT.TFile.Open(fileName)
EtaD_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaD")
EtaD_UP = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaD_BetaUpAtL4DT")
EtaD_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaD_BetaDownAtL4DT")
EtaE_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaE")
EtaE_UP = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaE_BetaUpAtL4DT")
EtaE_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaE_BetaDownAtL4DT")
EtaF_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaF")
EtaF_UP = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaF_BetaUpAtL4DT")
EtaF_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaF_BetaDownAtL4DT")

EtaD_Nom_ProfY = EtaD_Nom.RebinY(2).ProfileY().ProjectionX()
EtaD_UP_ProfY = EtaD_UP.RebinY(2).ProfileY().ProjectionX()
EtaD_DOWN_ProfY = EtaD_DOWN.RebinY(2).ProfileY().ProjectionX()
EtaE_Nom_ProfY = EtaE_Nom.RebinY(2).ProfileY().ProjectionX()
EtaE_UP_ProfY = EtaE_UP.RebinY(2).ProfileY().ProjectionX()
EtaE_DOWN_ProfY = EtaE_DOWN.RebinY(2).ProfileY().ProjectionX()
EtaF_Nom_ProfY = EtaF_Nom.RebinY(2).ProfileY().ProjectionX()
EtaF_UP_ProfY = EtaF_UP.RebinY(2).ProfileY().ProjectionX()
EtaF_DOWN_ProfY = EtaF_DOWN.RebinY(2).ProfileY().ProjectionX()

# --------------------------------------------------
EtaD_UpRatio = EtaD_UP_ProfY.Clone()
EtaD_UpRatio.Divide(EtaD_Nom_ProfY)

EtaD_DownRatio = EtaD_DOWN_ProfY.Clone()
EtaD_DownRatio.Divide(EtaD_Nom_ProfY)

EtaD_UpRatio.Draw("COLZTEXT40")
EtaD_UpRatio.SetMarkerColor(2)
EtaD_UpRatio.SetMarkerStyle(20)
EtaD_UpRatio.SetStats(0)
EtaD_UpRatio.SetLineColor(2)
EtaD_UpRatio.GetYaxis().SetTitle("Ratio")
EtaD_DownRatio.Draw("SAMETEXT40")
EtaD_DownRatio.SetMarkerStyle(20)
EtaD_DownRatio.SetMarkerColor(3)
EtaD_DownRatio.SetLineColor(3)
legRatiosSyst_EtaD =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_EtaD.SetTextFont(42)
legRatiosSyst_EtaD.SetTextSize(0.035)
legRatiosSyst_EtaD.SetBorderSize(1);
legRatiosSyst_EtaD.SetLineColor(0);
legRatiosSyst_EtaD.SetLineStyle(1);
legRatiosSyst_EtaD.SetLineWidth(1);
legRatiosSyst_EtaD.SetFillColor(0);
legRatiosSyst_EtaD.SetFillStyle(1001);
legRatiosSyst_EtaD.AddEntry(EtaD_UpRatio,"Up systematics / nominal","LP")
legRatiosSyst_EtaD.AddEntry(EtaD_DownRatio,"Down systematics / nominal","LP")
legRatiosSyst_EtaD.SetFillColorAlpha(ROOT.kWhite, 0.5)
legRatiosSyst_EtaD.Draw("SAME")
can.SaveAs("CodeV"+str(codeVersion)+"_TriggerEff_RatiosSyst_EtaDAtL4.png")

EtaD_Nom_ProfY.Draw("COLZ")
EtaD_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
EtaD_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
EtaD_UP_ProfY.Draw("SAME")
EtaD_UP_ProfY.SetLineColor(2)
EtaD_DOWN_ProfY.Draw("SAME")
EtaD_DOWN_ProfY.SetLineColor(3)
EtaD_Nom_ProfY.SetStats(0)
legTriggerEff_EtaD =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_EtaD.SetTextFont(42)
legTriggerEff_EtaD.SetTextSize(0.035)
legTriggerEff_EtaD.SetBorderSize(1);
legTriggerEff_EtaD.SetLineColor(0);
legTriggerEff_EtaD.SetLineStyle(1);
legTriggerEff_EtaD.SetLineWidth(1);
legTriggerEff_EtaD.SetFillColor(0);
legTriggerEff_EtaD.SetFillStyle(1001);
legTriggerEff_EtaD.AddEntry(EtaD_Nom_ProfY,"Nominal","LP")
legTriggerEff_EtaD.AddEntry(EtaD_UP_ProfY,"Up systematics","LP")
legTriggerEff_EtaD.AddEntry(EtaD_DOWN_ProfY,"Down systematics","LP")
legTriggerEff_EtaD.SetFillColorAlpha(ROOT.kWhite, 0.5)
legTriggerEff_EtaD.Draw("SAME")
can.SaveAs("CodeV"+str(codeVersion)+"_TriggerEff_Syst_EtaDAtL4.png")

can.Clear()

# --------------------------------------------------
EtaE_UpRatio = EtaE_UP_ProfY.Clone()
EtaE_UpRatio.Divide(EtaE_Nom_ProfY)

EtaE_DownRatio = EtaE_DOWN_ProfY.Clone()
EtaE_DownRatio.Divide(EtaE_Nom_ProfY)

EtaE_UpRatio.Draw("COLZTEXT40")
EtaE_UpRatio.SetMarkerColor(2)
EtaE_UpRatio.SetMarkerStyle(20)
EtaE_UpRatio.SetStats(0)
EtaE_UpRatio.SetLineColor(2)
EtaE_UpRatio.GetYaxis().SetTitle("Ratio")
EtaE_DownRatio.Draw("SAMETEXT40")
EtaE_DownRatio.SetMarkerStyle(20)
EtaE_DownRatio.SetMarkerColor(3)
EtaE_DownRatio.SetLineColor(3)
legRatiosSyst_EtaE =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_EtaE.SetTextFont(42)
legRatiosSyst_EtaE.SetTextSize(0.035)
legRatiosSyst_EtaE.SetBorderSize(1);
legRatiosSyst_EtaE.SetLineColor(0);
legRatiosSyst_EtaE.SetLineStyle(1);
legRatiosSyst_EtaE.SetLineWidth(1);
legRatiosSyst_EtaE.SetFillColor(0);
legRatiosSyst_EtaE.SetFillStyle(1001);
legRatiosSyst_EtaE.AddEntry(EtaE_UpRatio,"Up systematics / nominal","LP")
legRatiosSyst_EtaE.AddEntry(EtaE_DownRatio,"Down systematics / nominal","LP")
legRatiosSyst_EtaE.SetFillColorAlpha(ROOT.kWhite, 0.5)
legRatiosSyst_EtaE.Draw("SAME")
can.SaveAs("CodeV"+str(codeVersion)+"_TriggerEff_RatiosSyst_EtaEAtL4.png")
#print("EtaE_DownRatio content:")
#for xBin in range(1,EtaE_DownRatio.GetNbinsX()+1) :
  #"" + str(round(EtaE_UpRatio.GetXaxis().GetBinUpEdge(xBin),2)) +
#  print(str(round(EtaE_DownRatio.GetBinContent(xBin),2))+",")

EtaE_Nom_ProfY.Draw("COLZ")
EtaE_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
EtaE_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
EtaE_UP_ProfY.Draw("SAME")
EtaE_UP_ProfY.SetLineColor(2)
EtaE_DOWN_ProfY.Draw("SAME")
EtaE_DOWN_ProfY.SetLineColor(3)
EtaE_Nom_ProfY.SetStats(0)
legTriggerEff_EtaE =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_EtaE.SetTextFont(42)
legTriggerEff_EtaE.SetTextSize(0.035)
legTriggerEff_EtaE.SetBorderSize(1);
legTriggerEff_EtaE.SetLineColor(0);
legTriggerEff_EtaE.SetLineStyle(1);
legTriggerEff_EtaE.SetLineWidth(1);
legTriggerEff_EtaE.SetFillColor(0);
legTriggerEff_EtaE.SetFillStyle(1001);
legTriggerEff_EtaE.AddEntry(EtaE_Nom_ProfY,"Nominal","LP")
legTriggerEff_EtaE.AddEntry(EtaE_UP_ProfY,"Up systematics","LP")
legTriggerEff_EtaE.AddEntry(EtaE_DOWN_ProfY,"Down systematics","LP")
legTriggerEff_EtaE.SetFillColorAlpha(ROOT.kWhite, 0.5)
legTriggerEff_EtaE.Draw("SAME")
can.SaveAs("CodeV"+str(codeVersion)+"_TriggerEff_Syst_EtaEAtL4.png")

can.Clear()

# --------------------------------------------------
EtaF_UpRatio = EtaF_UP_ProfY.Clone()
EtaF_UpRatio.Divide(EtaF_Nom_ProfY)

EtaF_DownRatio = EtaF_DOWN_ProfY.Clone()
EtaF_DownRatio.Divide(EtaF_Nom_ProfY)

EtaF_UpRatio.Draw("COLZTEXT40")
EtaF_UpRatio.SetMarkerColor(2)
EtaF_UpRatio.SetMarkerStyle(20)
EtaF_UpRatio.SetStats(0)
EtaF_UpRatio.SetLineColor(2)
EtaF_UpRatio.GetYaxis().SetTitle("Ratio")
EtaF_DownRatio.Draw("SAMETEXT40")
EtaF_DownRatio.SetMarkerStyle(20)
EtaF_DownRatio.SetMarkerColor(3)
EtaF_DownRatio.SetLineColor(3)
legRatiosSyst_EtaF =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_EtaF.SetTextFont(42)
legRatiosSyst_EtaF.SetTextSize(0.035)
legRatiosSyst_EtaF.SetBorderSize(1);
legRatiosSyst_EtaF.SetLineColor(0);
legRatiosSyst_EtaF.SetLineStyle(1);
legRatiosSyst_EtaF.SetLineWidth(1);
legRatiosSyst_EtaF.SetFillColor(0);
legRatiosSyst_EtaF.SetFillStyle(1001);
legRatiosSyst_EtaF.AddEntry(EtaF_UpRatio,"Up systematics / nominal","LP")
legRatiosSyst_EtaF.AddEntry(EtaF_DownRatio,"Down systematics / nominal","LP")
legRatiosSyst_EtaF.SetFillColorAlpha(ROOT.kWhite, 0.5)
legRatiosSyst_EtaF.Draw("SAME")
can.SaveAs("CodeV"+str(codeVersion)+"_TriggerEff_RatiosSyst_EtaFAtL4.png")
#print("EtaF_DownRatio content:")
#for xBin in range(1,EtaF_DownRatio.GetNbinsX()+1) :
  #"" + str(round(EtaF_DownRatio.GetXaxis().GetBinUpEdge(xBin),2)) +
#  print(str(round(EtaF_DownRatio.GetBinContent(xBin),2))+",")

EtaF_Nom_ProfY.Draw("COLZ")
EtaF_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
EtaF_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
EtaF_UP_ProfY.Draw("SAME")
EtaF_UP_ProfY.SetLineColor(2)
EtaF_DOWN_ProfY.Draw("SAME")
EtaF_DOWN_ProfY.SetLineColor(3)
EtaF_Nom_ProfY.SetStats(0)
legTriggerEff_EtaF =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_EtaF.SetTextFont(42)
legTriggerEff_EtaF.SetTextSize(0.035)
legTriggerEff_EtaF.SetBorderSize(1);
legTriggerEff_EtaF.SetLineColor(0);
legTriggerEff_EtaF.SetLineStyle(1);
legTriggerEff_EtaF.SetLineWidth(1);
legTriggerEff_EtaF.SetFillColor(0);
legTriggerEff_EtaF.SetFillStyle(1001);
legTriggerEff_EtaF.AddEntry(EtaF_Nom_ProfY,"Nominal","LP")
legTriggerEff_EtaF.AddEntry(EtaF_UP_ProfY,"Up systematics","LP")
legTriggerEff_EtaF.AddEntry(EtaF_DOWN_ProfY,"Down systematics","LP")
legTriggerEff_EtaF.SetFillColorAlpha(ROOT.kWhite, 0.5)
legTriggerEff_EtaF.Draw("SAME")
can.SaveAs("CodeV"+str(codeVersion)+"_TriggerEff_Syst_EtaFAtL4.png")

can.Clear()
