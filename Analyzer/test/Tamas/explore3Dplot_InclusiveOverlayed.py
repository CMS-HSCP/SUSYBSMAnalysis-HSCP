import ROOT

ROOT.gROOT.SetStyle("Plain")
#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.SetBatch(True)

ROOT.gStyle.SetPadTopMargin(0.07);
ROOT.gStyle.SetPadBottomMargin(0.1);
ROOT.gStyle.SetPadLeftMargin(0.15);
ROOT.gStyle.SetPadRightMargin(0.13);

can = ROOT.TCanvas("newname","newname",800,800)
f = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-1800_CodeV42p8_v1.root")

# Signal fail
Histo3Dfail = f.Get("HSCParticleAnalyzer/BaseName/PostS_ProbQNoL1VsIasVsPt")
Histo3Dfail.GetXaxis().SetRange(1,18) # makes it the FAIL region
projFail = Histo3Dfail.Project3D("YZ")
#proj.GetZaxis().SetLabelOffset(.7)
#proj.GetYaxis().SetRangeUser(0.3,1.0)
projFail.Draw("COLZ")
projFail.SetTitle("FAIL region (F < 0.9)")
projFail.SetStats(0)
can.SaveAs("HSCPgluino_GiVsPt_FAIL.png")

# Signal pass
Histo3Dpass = f.Get("HSCParticleAnalyzer/BaseName/PostS_ProbQNoL1VsIasVsPt")
Histo3Dpass.GetXaxis().SetRange(19,21) # makes it the PASS region
projPass = Histo3Dpass.Project3D("YZ")
#proj.GetXaxis().SetRangeUser(1,300)
#proj.GetYaxis().SetRangeUser(0.3,1.0)
#proj.Draw("COLZTEXT")
projPass.Draw("COLZ")
projPass.SetTitle("PASS region (F > 0.9)")
projPass.SetStats(0)
can.SaveAs("HSCPgluino_GiVsPt_PASS.png")

# Signal fail log Gi
canLog = ROOT.TCanvas("newname2","newname2",800,800)
canLog.SetLogy()
legProjXInRegions =  ROOT.TLegend(.5,.75,.80,.9,"","brNDC")
legProjXInRegions.SetTextFont(42)
legProjXInRegions.SetTextSize(0.035)
legProjXInRegions.SetBorderSize(1);
legProjXInRegions.SetLineColor(0);
legProjXInRegions.SetLineStyle(1);
legProjXInRegions.SetLineWidth(1);
legProjXInRegions.SetFillColor(0);
legProjXInRegions.SetFillStyle(1001);


Histo3Dfail.GetXaxis().SetRange(1,18) # makes it the FAIL region
projFail = Histo3Dfail.Project3D("YZ")
projFail.SetTitle("FAIL region (F < 0.9)")
projFail.GetXaxis().SetRangeUser(200,4000)
projFailPt1 = projFail.ProjectionY("name1").Rebin(2)
projFailPt1.SetStats(0)
projFailPt1.SetMarkerStyle(20)
projFailPt1.SetLineColor(1)
projFailPt1.SetMarkerColor(1)
projFailPt1.SetMaximum(2)
projFailPt1.DrawClone("SAMEP")

legProjXInRegions.AddEntry(projFailPt1, "p_{T} = 200 - inf GeV","LP")

projFail.GetXaxis().SetRangeUser(300,4000)
projFailPt2 = projFail.ProjectionY("name2").Rebin(2)
projFailPt2.SetStats(0)
projFailPt2.SetMarkerStyle(20)
projFailPt2.SetLineColor(2)
projFailPt2.SetMarkerColor(2)
projFailPt2.DrawClone("SAMEP")
legProjXInRegions.AddEntry(projFailPt2, "p_{T} = 300 - inf  GeV","LP")

projFail.GetXaxis().SetRangeUser(400,4000)
projFailPt3 = projFail.ProjectionY("name3").Rebin(2)
projFailPt3.SetStats(0)
projFailPt3.SetMarkerStyle(20)
projFailPt3.SetLineColor(3)
projFailPt3.SetMarkerColor(3)
projFailPt3.DrawClone("SAMEP")

legProjXInRegions.AddEntry(projFailPt3, "p_{T} = 400 - inf GeV","LP")
legProjXInRegions.Draw("SAMEP")
canLog.SaveAs("HSCPgluino_Gi_PtBins_FAIL.png")

# Signal PASS log Gi
canLog2 = ROOT.TCanvas("newname2-pass","newname2-pass",800,800)
canLog2.SetLogy()
Histo3Dpass.GetXaxis().SetRange(19,21) # makes it the PASS region
projPass = Histo3Dpass.Project3D("YZ")
projPass.SetTitle("PASS region (F > 0.9)")

projPass.GetXaxis().SetRangeUser(240,4000)
projPassPt1 = projPass.ProjectionY("name1-pass").Rebin(2)
projPassPt1.SetStats(0)
projPassPt1.SetMarkerStyle(20)
projPassPt1.SetLineColor(1)
projPassPt1.SetMarkerColor(1)
projPassPt1.SetMaximum(200)
projPassPt1.DrawClone("SAMEP")

projPass.GetXaxis().SetRangeUser(320,4000)
projPassPt2 = projPass.ProjectionY("name2-pass").Rebin(2)
projPassPt2.SetStats(0)
projPassPt2.SetMarkerStyle(20)
projPassPt2.SetLineColor(2)
projPassPt2.SetMarkerColor(2)
projPassPt2.DrawClone("SAMEP")

projPass.GetXaxis().SetRangeUser(400,4000)
projPassPt3 = projPass.ProjectionY("name3-pass").Rebin(2)
projPassPt3.SetStats(0)
projPassPt3.SetMarkerStyle(20)
projPassPt3.SetLineColor(3)
projPassPt3.SetMarkerColor(3)
projPassPt3.DrawClone("SAMEP")

legProjXInRegions.Draw("SAMEP")
canLog2.SaveAs("HSCPgluino_Gi_PtBins_PASS.png")



# Now repeat everything for the background
f2 = ROOT.TFile.Open("crab_Analysis_SingleMuon_RunPhase1_CodeV42p8_v1.root")
Histo3Dfail = f2.Get("HSCParticleAnalyzer/BaseName/PostS_ProbQNoL1VsIasVsPt")
# Signal fail log Gi
canLog3 = ROOT.TCanvas("canLog3anewname2","canLog3anewname2",800,800)
canLog3.SetLogy()
legProjXInRegions =  ROOT.TLegend(.5,.75,.80,.9,"","brNDC")
legProjXInRegions.SetTextFont(42)
legProjXInRegions.SetTextSize(0.035)
legProjXInRegions.SetBorderSize(1);
legProjXInRegions.SetLineColor(0);
legProjXInRegions.SetLineStyle(1);
legProjXInRegions.SetLineWidth(1);
legProjXInRegions.SetFillColor(0);
legProjXInRegions.SetFillStyle(1001);
canLog3.SetLogy()

Histo3Dfail.GetXaxis().SetRange(1,18) # makes it the FAIL region
projFail = Histo3Dfail.Project3D("YZ")
projFail.SetTitle("FAIL region (F < 0.9)")
projFail.GetXaxis().SetRangeUser(240,4000)
projFailPt1 = projFail.ProjectionY("name1").Rebin(2)
projFailPt1.SetStats(0)
projFailPt1.SetMarkerStyle(20)
projFailPt1.SetLineColor(1)
projFailPt1.SetMarkerColor(1)
projFailPt1.SetMaximum(20000000)
projFailPt1.DrawClone("SAMEP")

legProjXInRegions.AddEntry(projFailPt1, "p_{T} = 240 - inf GeV","LP")

projFail.GetXaxis().SetRangeUser(320,4000)
projFailPt2 = projFail.ProjectionY("name2").Rebin(2)
projFailPt2.SetStats(0)
projFailPt2.SetMarkerStyle(20)
projFailPt2.SetLineColor(2)
projFailPt2.SetMarkerColor(2)
projFailPt2.DrawClone("SAMEP")
legProjXInRegions.AddEntry(projFailPt2, "p_{T} = 320 - inf  GeV","LP")

projFail.GetXaxis().SetRangeUser(400,4000)
projFailPt3 = projFail.ProjectionY("name3").Rebin(2)
projFailPt3.SetStats(0)
projFailPt3.SetMarkerStyle(20)
projFailPt3.SetLineColor(3)
projFailPt3.SetMarkerColor(3)
projFailPt3.DrawClone("SAMEP")

legProjXInRegions.AddEntry(projFailPt3, "p_{T} = 400 - inf GeV","LP")
legProjXInRegions.Draw("SAMEP")
canLog3.SaveAs("AllBkg_data_Gi_PtBinsIncl_FAIL.png")

# Signal PASS log Gi
Histo3Dpass = f2.Get("HSCParticleAnalyzer/BaseName/PostS_ProbQNoL1VsIasVsPt")
canLog32 = ROOT.TCanvas("canLog3anewname2-pass","canLog3anewname2-pass",800,800)
canLog32.SetLogy()
Histo3Dpass.GetXaxis().SetRange(19,21) # makes it the PASS region
projPass = Histo3Dpass.Project3D("YZ")
projPass.SetTitle("PASS region (F > 0.9)")

projPass.GetXaxis().SetRangeUser(240,4000)
projPassPt1 = projPass.ProjectionY("name1-pass").Rebin(2)
projPassPt1.SetStats(0)
projPassPt1.SetMarkerStyle(20)
projPassPt1.SetLineColor(1)
projPassPt1.SetMarkerColor(1)
projPassPt1.SetMaximum(20000000)
projPassPt1.DrawClone("SAMEP")

projPass.GetXaxis().SetRangeUser(320,4000)
projPassPt2 = projPass.ProjectionY("name2-pass").Rebin(2)
projPassPt2.SetStats(0)
projPassPt2.SetMarkerStyle(20)
projPassPt2.SetLineColor(2)
projPassPt2.SetMarkerColor(2)
projPassPt2.DrawClone("SAMEP")

projPass.GetXaxis().SetRangeUser(400,4000)
projPassPt3 = projPass.ProjectionY("name3-pass").Rebin(2)
projPassPt3.SetStats(0)
projPassPt3.SetMarkerStyle(20)
projPassPt3.SetLineColor(3)
projPassPt3.SetMarkerColor(3)
projPassPt3.DrawClone("SAMEP")

legProjXInRegions.Draw("SAMEP")
canLog32.SaveAs("AllBkg_data_Gi_PtBinsIncl_PASS.png")
