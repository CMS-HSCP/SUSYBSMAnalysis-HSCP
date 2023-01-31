import ROOT

ROOT.gStyle.SetPaintTextFormat(".2g");

can = ROOT.TCanvas("newname","newname",800,800)
f = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-All_CodeV42p6_v1.root")
EtaA_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaA")
EtaA_UP = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaA_BetaUp")
EtaA_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaA_BetaDown")
EtaB_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaB")
EtaB_UP = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaB_BetaUp")
EtaB_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaB_BetaDown")
EtaC_Nom = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaC")
EtaC_UP = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaC_BetaUp")
EtaC_DOWN = f.Get("HSCParticleAnalyzer/BaseName/PostPreS_TriggerMuon50VsBeta_EtaC_BetaDown")

EtaA_Nom_ProfY = EtaA_Nom.RebinY(2).ProfileY().ProjectionX()
EtaA_UP_ProfY = EtaA_UP.RebinY(2).ProfileY().ProjectionX()
EtaA_DOWN_ProfY = EtaA_DOWN.RebinY(2).ProfileY().ProjectionX()
EtaB_Nom_ProfY = EtaB_Nom.RebinY(2).ProfileY().ProjectionX()
EtaB_UP_ProfY = EtaB_UP.RebinY(2).ProfileY().ProjectionX()
EtaB_DOWN_ProfY = EtaB_DOWN.RebinY(2).ProfileY().ProjectionX()
EtaC_Nom_ProfY = EtaC_Nom.RebinY(2).ProfileY().ProjectionX()
EtaC_UP_ProfY = EtaC_UP.RebinY(2).ProfileY().ProjectionX()
EtaC_DOWN_ProfY = EtaC_DOWN.RebinY(2).ProfileY().ProjectionX()

# --------------------------------------------------
EtaA_UpRatio = EtaA_UP_ProfY.Clone()
EtaA_UpRatio.Divide(EtaA_Nom_ProfY)

EtaA_DownRatio = EtaA_DOWN_ProfY.Clone()
EtaA_DownRatio.Divide(EtaA_Nom_ProfY)

EtaA_UpRatio.Draw("COLZTEXT40")
EtaA_UpRatio.SetMarkerColor(2)
EtaA_UpRatio.SetMarkerStyle(20)
EtaA_UpRatio.SetStats(0)
EtaA_UpRatio.SetLineColor(2)
EtaA_UpRatio.GetYaxis().SetTitle("Ratio")
EtaA_DownRatio.Draw("SAMETEXT40")
EtaA_DownRatio.SetMarkerStyle(20)
EtaA_DownRatio.SetMarkerColor(3)
EtaA_DownRatio.SetLineColor(3)
legRatiosSyst_EtaA =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_EtaA.SetTextFont(42)
legRatiosSyst_EtaA.SetTextSize(0.035)
legRatiosSyst_EtaA.SetBorderSize(1);
legRatiosSyst_EtaA.SetLineColor(0);
legRatiosSyst_EtaA.SetLineStyle(1);
legRatiosSyst_EtaA.SetLineWidth(1);
legRatiosSyst_EtaA.SetFillColor(0);
legRatiosSyst_EtaA.SetFillStyle(1001);
legRatiosSyst_EtaA.AddEntry(EtaA_UpRatio,"Up systematics / nominal","LP")
legRatiosSyst_EtaA.AddEntry(EtaA_DownRatio,"Down systematics / nominal","LP")
legRatiosSyst_EtaA.Draw("SAME")
can.SaveAs("TriggerEff_RatiosSyst_EtaA.png")

EtaA_Nom_ProfY.Draw("COLZ")
EtaA_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
EtaA_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
EtaA_UP_ProfY.Draw("SAME")
EtaA_UP_ProfY.SetLineColor(2)
EtaA_DOWN_ProfY.Draw("SAME")
EtaA_DOWN_ProfY.SetLineColor(3)
EtaA_Nom_ProfY.SetStats(0)
legTriggerEff_EtaA =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_EtaA.SetTextFont(42)
legTriggerEff_EtaA.SetTextSize(0.035)
legTriggerEff_EtaA.SetBorderSize(1);
legTriggerEff_EtaA.SetLineColor(0);
legTriggerEff_EtaA.SetLineStyle(1);
legTriggerEff_EtaA.SetLineWidth(1);
legTriggerEff_EtaA.SetFillColor(0);
legTriggerEff_EtaA.SetFillStyle(1001);
legTriggerEff_EtaA.AddEntry(EtaA_Nom_ProfY,"Nominal","LP")
legTriggerEff_EtaA.AddEntry(EtaA_UP_ProfY,"Up systematics","LP")
legTriggerEff_EtaA.AddEntry(EtaA_DOWN_ProfY,"Down systematics","LP")
legTriggerEff_EtaA.Draw("SAME")
can.SaveAs("TriggerEff_Syst_EtaA.png")

can.Clear()

# --------------------------------------------------
EtaB_UpRatio = EtaB_UP_ProfY.Clone()
EtaB_UpRatio.Divide(EtaB_Nom_ProfY)

EtaB_DownRatio = EtaB_DOWN_ProfY.Clone()
EtaB_DownRatio.Divide(EtaB_Nom_ProfY)

EtaB_UpRatio.Draw("COLZTEXT40")
EtaB_UpRatio.SetMarkerColor(2)
EtaB_UpRatio.SetMarkerStyle(20)
EtaB_UpRatio.SetStats(0)
EtaB_UpRatio.SetLineColor(2)
EtaB_UpRatio.GetYaxis().SetTitle("Ratio")
EtaB_DownRatio.Draw("SAMETEXT40")
EtaB_DownRatio.SetMarkerStyle(20)
EtaB_DownRatio.SetMarkerColor(3)
EtaB_DownRatio.SetLineColor(3)
legRatiosSyst_EtaB =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_EtaB.SetTextFont(42)
legRatiosSyst_EtaB.SetTextSize(0.035)
legRatiosSyst_EtaB.SetBorderSize(1);
legRatiosSyst_EtaB.SetLineColor(0);
legRatiosSyst_EtaB.SetLineStyle(1);
legRatiosSyst_EtaB.SetLineWidth(1);
legRatiosSyst_EtaB.SetFillColor(0);
legRatiosSyst_EtaB.SetFillStyle(1001);
legRatiosSyst_EtaB.AddEntry(EtaB_UpRatio,"Up systematics / nominal","LP")
legRatiosSyst_EtaB.AddEntry(EtaB_DownRatio,"Down systematics / nominal","LP")
legRatiosSyst_EtaB.Draw("SAME")
can.SaveAs("TriggerEff_RatiosSyst_EtaB.png")
#print("EtaB_DownRatio content:")
#for xBin in range(1,EtaB_DownRatio.GetNbinsX()+1) :
  #"" + str(round(EtaB_UpRatio.GetXaxis().GetBinUpEdge(xBin),2)) +
#  print(str(round(EtaB_DownRatio.GetBinContent(xBin),2))+",")

EtaB_Nom_ProfY.Draw("COLZ")
EtaB_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
EtaB_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
EtaB_UP_ProfY.Draw("SAME")
EtaB_UP_ProfY.SetLineColor(2)
EtaB_DOWN_ProfY.Draw("SAME")
EtaB_DOWN_ProfY.SetLineColor(3)
EtaB_Nom_ProfY.SetStats(0)
legTriggerEff_EtaB =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_EtaB.SetTextFont(42)
legTriggerEff_EtaB.SetTextSize(0.035)
legTriggerEff_EtaB.SetBorderSize(1);
legTriggerEff_EtaB.SetLineColor(0);
legTriggerEff_EtaB.SetLineStyle(1);
legTriggerEff_EtaB.SetLineWidth(1);
legTriggerEff_EtaB.SetFillColor(0);
legTriggerEff_EtaB.SetFillStyle(1001);
legTriggerEff_EtaB.AddEntry(EtaB_Nom_ProfY,"Nominal","LP")
legTriggerEff_EtaB.AddEntry(EtaB_UP_ProfY,"Up systematics","LP")
legTriggerEff_EtaB.AddEntry(EtaB_DOWN_ProfY,"Down systematics","LP")
legTriggerEff_EtaB.Draw("SAME")
can.SaveAs("TriggerEff_Syst_EtaB.png")

can.Clear()

# --------------------------------------------------
EtaC_UpRatio = EtaC_UP_ProfY.Clone()
EtaC_UpRatio.Divide(EtaC_Nom_ProfY)

EtaC_DownRatio = EtaC_DOWN_ProfY.Clone()
EtaC_DownRatio.Divide(EtaC_Nom_ProfY)

EtaC_UpRatio.Draw("COLZTEXT40")
EtaC_UpRatio.SetMarkerColor(2)
EtaC_UpRatio.SetMarkerStyle(20)
EtaC_UpRatio.SetStats(0)
EtaC_UpRatio.SetLineColor(2)
EtaC_UpRatio.GetYaxis().SetTitle("Ratio")
EtaC_DownRatio.Draw("SAMETEXT40")
EtaC_DownRatio.SetMarkerStyle(20)
EtaC_DownRatio.SetMarkerColor(3)
EtaC_DownRatio.SetLineColor(3)
legRatiosSyst_EtaC =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legRatiosSyst_EtaC.SetTextFont(42)
legRatiosSyst_EtaC.SetTextSize(0.035)
legRatiosSyst_EtaC.SetBorderSize(1);
legRatiosSyst_EtaC.SetLineColor(0);
legRatiosSyst_EtaC.SetLineStyle(1);
legRatiosSyst_EtaC.SetLineWidth(1);
legRatiosSyst_EtaC.SetFillColor(0);
legRatiosSyst_EtaC.SetFillStyle(1001);
legRatiosSyst_EtaC.AddEntry(EtaC_UpRatio,"Up systematics / nominal","LP")
legRatiosSyst_EtaC.AddEntry(EtaC_DownRatio,"Down systematics / nominal","LP")
legRatiosSyst_EtaC.Draw("SAME")
can.SaveAs("TriggerEff_RatiosSyst_EtaC.png")
#print("EtaC_DownRatio content:")
#for xBin in range(1,EtaC_DownRatio.GetNbinsX()+1) :
  #"" + str(round(EtaC_DownRatio.GetXaxis().GetBinUpEdge(xBin),2)) +
#  print(str(round(EtaC_DownRatio.GetBinContent(xBin),2))+",")

EtaC_Nom_ProfY.Draw("COLZ")
EtaC_Nom_ProfY.GetYaxis().SetRangeUser(0.,1.)
EtaC_Nom_ProfY.GetYaxis().SetTitle("Trigger efficiency")
EtaC_UP_ProfY.Draw("SAME")
EtaC_UP_ProfY.SetLineColor(2)
EtaC_DOWN_ProfY.Draw("SAME")
EtaC_DOWN_ProfY.SetLineColor(3)
EtaC_Nom_ProfY.SetStats(0)
legTriggerEff_EtaC =  ROOT.TLegend(.30,.80,.80,.885,"","brNDC")
legTriggerEff_EtaC.SetTextFont(42)
legTriggerEff_EtaC.SetTextSize(0.035)
legTriggerEff_EtaC.SetBorderSize(1);
legTriggerEff_EtaC.SetLineColor(0);
legTriggerEff_EtaC.SetLineStyle(1);
legTriggerEff_EtaC.SetLineWidth(1);
legTriggerEff_EtaC.SetFillColor(0);
legTriggerEff_EtaC.SetFillStyle(1001);
legTriggerEff_EtaC.AddEntry(EtaC_Nom_ProfY,"Nominal","LP")
legTriggerEff_EtaC.AddEntry(EtaC_UP_ProfY,"Up systematics","LP")
legTriggerEff_EtaC.AddEntry(EtaC_DOWN_ProfY,"Down systematics","LP")
legTriggerEff_EtaC.Draw("SAME")
can.SaveAs("TriggerEff_Syst_EtaC.png")

can.Clear()
