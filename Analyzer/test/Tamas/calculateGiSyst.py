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
can.SetLogy()
fileName = "crab_Analysis_2018_AllBackground_CodeV42p6_v1.root"
mc = ROOT.TFile.Open(fileName)
data = ROOT.TFile.Open("crab_Analysis_SingleMuon_RunPhase1_CodeV42p6_v1.root")
Gi_MC = mc.Get("HSCParticleAnalyzer/BaseName/PostPreS_Ias_CR_lowPt")
Gi_Data = data.Get("HSCParticleAnalyzer/BaseName/PostPreS_Ias_CR_lowPt")

Gi_MC.SetStats(0)
Gi_MC.SetMarkerColor(1)
Gi_MC.SetLineColor(1)
Gi_MC.SetMarkerStyle(20)
Gi_MC.Scale(1/Gi_MC.Integral())
Gi_MC.GetYaxis().SetTitle("Normalized Tracks / 0.1")
#Gi_MC.GetYaxis().SetTitleOffset(1.5)
Gi_MC.Draw()

Gi_Data.SetStats(0)
Gi_Data.SetMarkerStyle(20)
Gi_Data.SetMarkerColor(2)
Gi_Data.SetLineColor(2)
Gi_Data.Scale(1/Gi_Data.Integral())
Gi_Data.Draw()

legIasForCR =  ROOT.TLegend(.60,.80,.80,.90,"","brNDC")
legIasForCR.SetTextFont(42)
legIasForCR.SetTextSize(0.035)
legIasForCR.SetBorderSize(1);
legIasForCR.SetLineColor(0);
legIasForCR.SetLineStyle(1);
legIasForCR.SetLineWidth(1);
legIasForCR.SetFillColor(0);
legIasForCR.SetFillStyle(1001);
legIasForCR.AddEntry(Gi_Data,"Data CR","LP")
legIasForCR.AddEntry(Gi_MC,"SM MC CR","LP")

tex2 = ROOT.TLatex(0.13,0.94,"CMS");
tex2.SetNDC();
tex2.SetTextFont(61);
tex2.SetTextSize(0.0675);
tex2.SetLineWidth(2);

tex3 = ROOT.TLatex(0.27,0.94,"Internal");
tex3.SetNDC();
tex3.SetTextFont(52);
tex3.SetTextSize(0.0485);
tex3.SetLineWidth(2);

tex4 = ROOT.TLatex()
tex4 = ROOT.TLatex(0.6,0.95,"After pre-selection")
tex4.SetNDC();
tex4.SetTextFont(52);
tex4.SetTextSize(0.045);
tex4.SetLineWidth(2);

codeVersion = fileName[fileName.find("CodeV")+5:fileName.find("CodeV")+9]
fileVersion = fileName[fileName.find("2018")+5:fileName.find("CodeV")+9]
tex5 = ROOT.TLatex(0.07,0.03,fileVersion);
tex5.SetNDC();
tex5.SetTextFont(52);
tex5.SetTextSize(0.0185);
tex5.SetLineWidth(2);


rp2 = ROOT.TRatioPlot(Gi_Data,Gi_MC,"divsym") #, "diffsigerrasym"
 
rp2.SetH1DrawOpt("P");
rp2.SetH2DrawOpt("P");

rp2.Draw("SAME")


rp2.SetLeftMargin(0.13);
rp2.SetRightMargin(0.05);
rp2.SetUpTopMargin(0.1);
rp2.SetLowTopMargin(0.02);
rp2.SetLowBottomMargin(0.35);


rp2.GetLowerRefGraph().SetMinimum(0)
rp2.GetLowerRefGraph().SetMaximum(1.05);
rp2.GetLowerRefGraph().SetMarkerStyle(20)
rp2.GetLowerRefGraph().SetMarkerSize(1);
rp2.GetLowYaxis().SetNdivisions(510);
rp2.GetLowerRefYaxis().SetTitle("Ratio");
rp2.GetLowerRefYaxis().SetTitleSize(0.05);
rp2.GetLowerRefYaxis().SetTitleOffset(1);
rp2.GetLowerRefYaxis().SetLabelSize(0.035);


rp2.GetLowerRefXaxis().SetTitleSize(0.05);
rp2.GetLowerRefXaxis().SetTitleOffset(0.8);
rp2.GetLowerRefXaxis().SetLabelSize(0.035);

can.Update()


legIasForCR.Draw("SAME")
tex2.Draw("SAME")
tex3.Draw("SAME")
tex4.Draw("SAME")
tex5.Draw("SAME")
can.SaveAs("GiSysts_DataOverMC_Norm.png")

