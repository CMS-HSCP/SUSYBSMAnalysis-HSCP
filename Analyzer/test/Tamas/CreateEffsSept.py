import ROOT, sys, os, time, re, numpy, random
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

codeVersion = sys.argv[1]

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

ROOT.gStyle.SetPaintTextFormat("4.2f");

InputListSamples = [
"BackgroundSamplesCode"+codeVersion+".txt",
#"HSCPgluinoOnlyNeutralSamples_Code"+codeVersion+".txt",
"HSCPgluinoSamples_Code"+codeVersion+".txt",
#"HSCPgmsbStauSamples_Code"+codeVersion+".txt",
#"HSCPpairStauSamples_Code"+codeVersion+".txt",
#"HSCPstopOnlyNeutralSamples_Code"+codeVersion+".txt",
#"HSCPstopSamples_Code"+codeVersion+".txt",
#"HSCPtauPrime1Samples_Code"+codeVersion+".txt",
#"HSCPtauPrime2Samples_Code"+codeVersion+".txt",
]

SingleMuonSample = ROOT.TFile.Open("crab_Analysis_2018_SingleMuon_Run2018C_woProbQ_CodeV"+codeVersion+"_v1.root")
AllBcgSample = ROOT.TFile.Open("crab_Analysis_2018_AllBackground_woProbQ_CodeV"+codeVersion+"_v1.root")
Rhadron1800GeV = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-1800_woProbQ_CodeV"+codeVersion+"_v1.root")

cutValues = {
'BefPreS_Eta' : 1.0,
'N1_Eta' : 1.0,
'BefPreS_Pt' : 55.,
'N1_Pt' : 55.,
'BefPreS_Dxy' : 0.02,
'N1_Dxy' : 0.02,
'BefPreS_Dz' : 0.02,
'N1_Dz' : 0.02,
'BefPreS_Chi2oNdof' : 5.,
'N1_Chi2oNdof' : 5.,
'BefPreS_TNOM' : 10.,
'N1_TNOM' : 10.,
'BefPreS_TNOPH' : 2.,
'N1_TNOPH' : 2.,
'BefPreS_TNOHFraction' : .9,
'N1_TNOHFraction' : .9,
'N1_EoP' : 0.3,
'BefPreS_EoP' : 0.3,
'BefPreS_TIsol' : 15,
'N1_TIsol' : 15,
'BefPreS_Ih' : 3.47,
'N1_Ih' : 3.47,
'BefPreS_ProbXY' : 0.001,
'N1_ProbXY' : 0.001,
'BefPreS_MiniRelIsoAll' : 0.02,
'N1_MiniRelIsoAll' : 0.02,
}


sigma = 3

for i in range(0, SingleMuonSample.GetListOfKeys().GetEntries()):
  dirname = SingleMuonSample.GetListOfKeys().At(i).GetName()
  curr_dir = SingleMuonSample.GetDirectory(dirname)
  if not (curr_dir) :
    continue
  N1eff = 0
  BefPreSeffForEta = BefPreSeffForpT =  BefPreSeffForNumPixHits = BefPreSeffForValidFract = 0
  BefPreSeffForNumDeDx = BefPreSeffForChi2oDOF = BefPreSeffForEoP = BefPreSeffFordz = 0
  BefPreSeffFordxy = BefPreSeffForMiniIso = BefPreSeffForIh = BefPreSeffForProbXY = 0
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = SingleMuonSample.GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
#          if not ("N1_Eta" in keyname2) : continue
          newname = dirname + "/" + keyname+ "/" + keyname2
          if ("N1_ProbQVsIas" in keyname2 or "N1_pfType" in keyname2) : continue
          histo = SingleMuonSample.Get(newname)
          histoAllBcg = AllBcgSample.Get(newname)
          histoSignal = Rhadron1800GeV.Get(newname)
          histoDenomAllBcgFixedN1 = AllBcgSample.Get(dirname + "/" + keyname + "/N1_Eta")
          histoDenomAllBcgFixedBefPreS = AllBcgSample.Get(dirname + "/" + keyname + "/BefPreS_Eta")

          if not (cutValues.get(keyname2)) : continue
          Num = histo.Integral(histo.GetXaxis().FindBin(0),histo.GetXaxis().FindBin(cutValues[keyname2]))
          NumAllBcg = histoAllBcg.Integral(histoAllBcg.GetXaxis().FindBin(0),histoAllBcg.GetXaxis().FindBin(cutValues[keyname2]))
          NumSignal = histoSignal.Integral(histoSignal.GetXaxis().FindBin(0),histoSignal.GetXaxis().FindBin(cutValues[keyname2]))
          Denom = histo.Integral(histo.GetXaxis().FindBin(0),histo.GetXaxis().FindBin(histo.GetXaxis().GetXmax()))
          DenomAllBcg = histoAllBcg.Integral(histoAllBcg.GetXaxis().FindBin(0),histoAllBcg.GetXaxis().FindBin(histoAllBcg.GetXaxis().GetXmax()))
          DenomAllBcgFixedBefPreS = histoDenomAllBcgFixedBefPreS.Integral(histoDenomAllBcgFixedBefPreS.GetXaxis().FindBin(0),histoDenomAllBcgFixedBefPreS.GetXaxis().FindBin(histoDenomAllBcgFixedBefPreS.GetXaxis().GetXmax()))
          DenomAllBcgFixedN1 = histoDenomAllBcgFixedN1.Integral(histoDenomAllBcgFixedN1.GetXaxis().FindBin(0),histoDenomAllBcgFixedN1.GetXaxis().FindBin(histoDenomAllBcgFixedN1.GetXaxis().GetXmax()))
          
#          print(keyname2 + ": " + str(DenomAllBcgFixedN1))
          DenomSignal = histoSignal.Integral(histoSignal.GetXaxis().FindBin(0),histoSignal.GetXaxis().FindBin(histoSignal.GetXaxis().GetXmax()))
          if (Denom>0) :  Eff = Num / Denom
          else : Eff = 0
          if (DenomAllBcg>0) : EffAllBcg = NumAllBcg / DenomAllBcg
          else : EffAllBcg = 0
          if (DenomSignal>0) : EffSignal = NumSignal / DenomSignal
          else : EffSignal = 0
         
          if ("BefPreS_Eta" in keyname2) :
            BefPreSeffForEta = Eff
            BefPreSeffForEtaBcg = EffAllBcg
            BefPreSeffForEtaSignal = EffSignal
            BefPreSeffForEtaPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_Pt" in keyname2) :
            BefPreSeffForpT = 1-Eff
            BefPreSeffForpTBcg = 1-EffAllBcg
            BefPreSeffForpTSignal = 1-EffSignal
            BefPreSeffForpTPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_TNOPH" in keyname2) :
            BefPreSeffForNumPixHits = 1-Eff
            BefPreSeffForNumPixHitsBcg = 1-EffAllBcg
            BefPreSeffForNumPixHitsSignal = 1-EffSignal
            BefPreSeffForNumPixHitsPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_TNOHFraction" in keyname2) :
            BefPreSeffForValidFract = 1-Eff
            BefPreSeffForValidFractBcg = 1-EffAllBcg
            BefPreSeffForValidFractSignal = 1-EffSignal
            BefPreSeffForValidFractPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_TNOM" in keyname2) :
            BefPreSeffForNumDeDx = 1-Eff
            BefPreSeffForNumDeDxBcg = 1-EffAllBcg
            BefPreSeffForNumDeDxSignal = 1-EffSignal
            BefPreSeffForNumDeDxPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_Chi2oNdof" in keyname2) :
            BefPreSeffForChi2oDOF = Eff
            BefPreSeffForChi2oDOFBcg = EffAllBcg
            BefPreSeffForChi2oDOFSignal = EffSignal
            BefPreSeffForChi2oDOFPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_EoP" in keyname2) :
            BefPreSeffForEoP = Eff
            BefPreSeffForEoPBcg = EffAllBcg
            BefPreSeffForEoPSignal = EffSignal
            BefPreSeffForEoPPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_Dz" in keyname2) :
            BefPreSeffFordz = Eff
            BefPreSeffFordzBcg = EffAllBcg
            BefPreSeffFordzSignal = EffSignal
            BefPreSeffFordzPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_Dxy" in keyname2) :
            BefPreSeffFordxy = Eff
            BefPreSeffFordxyBcg = EffAllBcg
            BefPreSeffFordxySignal = EffSignal
            BefPreSeffFordxyPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_MiniRelIsoAll" in keyname2) :
            BefPreSeffForMiniIso = Eff
            BefPreSeffForMiniIsoBcg = EffAllBcg
            BefPreSeffForMiniIsoSignal = EffSignal
            BefPreSeffForMiniIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_Ih" in keyname2) :
            BefPreSeffForIh = 1-Eff
            BefPreSeffForIhBcg = 1-EffAllBcg
            BefPreSeffForIhSignal = 1-EffSignal
            BefPreSeffForIhPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          if ("BefPreS_ProbXY" in keyname2) :
            BefPreSeffForProbXY = 1-Eff
            BefPreSeffForProbXYBcg = 1-EffAllBcg
            BefPreSeffForProbXYSignal = 1-EffSignal
            BefPreSeffForProbXYPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedBefPreS))
          
          
          if ("N1_Eta" in keyname2) :
            N1effForEta = Eff
            N1effForEtaBcg = EffAllBcg
            N1effForEtaSignal = EffSignal
            N1effForEtaPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_Pt" in keyname2) :
            N1effForpT = 1-Eff
            N1effForpTBcg = 1-EffAllBcg
            N1effForpTSignal = 1-EffSignal
            N1effForpTPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_TNOPH" in keyname2) :
            N1effForNumPixHits = 1-Eff
            N1effForNumPixHitsBcg = 1-EffAllBcg
            N1effForNumPixHitsSignal = 1-EffSignal
            N1effForNumPixHitsPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_TNOHFraction" in keyname2) :
            N1effForValidFract = 1-Eff
            N1effForValidFractBcg = 1-EffAllBcg
            N1effForValidFractSignal = 1-EffSignal
            N1effForValidFractPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_TNOM" in keyname2) :
            N1effForNumDeDx = 1-Eff
            N1effForNumDeDxBcg = 1-EffAllBcg
            N1effForNumDeDxSignal = 1-EffSignal
            N1effForNumDeDxPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_Chi2oNdof" in keyname2) :
            N1effForChi2oDOF = Eff
            N1effForChi2oDOFBcg = EffAllBcg
            N1effForChi2oDOFSignal = EffSignal
            N1effForChi2oDOFPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_EoP" in keyname2) :
            N1effForEoP = Eff
            N1effForEoPBcg = EffAllBcg
            N1effForEoPSignal = EffSignal
            N1effForEoPPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_Dz" in keyname2) :
            N1effFordz = Eff
            N1effFordzBcg = EffAllBcg
            N1effFordzSignal = EffSignal
            N1effFordzPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_Dxy" in keyname2) :
            N1effFordxy = Eff
            N1effFordxyBcg = EffAllBcg
            N1effFordxySignal = EffSignal
            N1effFordxyPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_MiniRelIsoAll" in keyname2) :
            N1effForMiniIso = Eff
            N1effForMiniIsoBcg = EffAllBcg
            N1effForMiniIsoSignal = EffSignal
            N1effForMiniIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_Ih" in keyname2) :
            N1effForIh = 1-Eff
            N1effForIhBcg = 1-EffAllBcg
            N1effForIhSignal = 1-EffSignal
            N1effForIhPunzi = (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          if ("N1_ProbXY" in keyname2) :
            N1effForProbXY = 1-Eff
            N1effForProbXYBcg = 1-EffAllBcg
            N1effForProbXYSignal = 1-EffSignal
            N1effForProbXYPunzi= (1-EffSignal) / (sigma + numpy.sqrt(DenomAllBcgFixedN1))
          
      tex2 = ROOT.TLatex(0.13,0.94,"CMS");
      tex2.SetNDC();
      tex2.SetTextFont(61);
      tex2.SetTextSize(0.0675);
      tex2.SetLineWidth(2);

      tex3 = ROOT.TLatex(0.27,0.94,"Internal"); # for square plots
      tex3.SetNDC();
      tex3.SetTextFont(52);
      tex3.SetTextSize(0.0485);
      tex3.SetLineWidth(2);
      
      tex4 = ROOT.TLatex()

      tex5 = ROOT.TLatex(0.07,0.04,"Code version: "+codeVersion);
      tex5.SetNDC();
      tex5.SetTextFont(52);
      tex5.SetTextSize(0.0185);
      tex5.SetLineWidth(2);
      
      cstackedSummedBackgroundString = "cstackedSummedBackgroundString"+str(keyname2)
      canvas = ROOT.TCanvas(cstackedSummedBackgroundString, cstackedSummedBackgroundString, 800,800)
      
      EffForData = ROOT.TH1F("EffForData",";;Integral to (from) max (min) cut value /  Integral to inf",12,0.,12.)
      EffForData.Draw("HISTTEXT00")
      EffForData.SetBinContent(1,BefPreSeffForEta)
      EffForData.SetBinContent(2,BefPreSeffForpT)
      EffForData.SetBinContent(3,BefPreSeffForNumPixHits)
      EffForData.SetBinContent(4,BefPreSeffForValidFract)
      EffForData.SetBinContent(5,BefPreSeffForNumDeDx)
      EffForData.SetBinContent(6,BefPreSeffForChi2oDOF)
      EffForData.SetBinContent(7,BefPreSeffForEoP)
      EffForData.SetBinContent(8,BefPreSeffFordz)
      EffForData.SetBinContent(9,BefPreSeffFordxy)
      EffForData.SetBinContent(10,BefPreSeffForMiniIso)
      EffForData.SetBinContent(11,BefPreSeffForIh)
      EffForData.SetBinContent(12,BefPreSeffForProbXY)

      
      EffForData.SetLineColor(1)
      EffForData.SetMarkerColor(1)
#      EffForData.SetMarkerStyle(20)
      EffForData.SetStats(0)
      EffForData.GetXaxis().SetTitle("")
      EffForData.GetYaxis().SetRangeUser(0.,1.3)
      EffForData.GetXaxis().SetBinLabel(1,"Eta < 1.0")
      EffForData.GetXaxis().SetBinLabel(2,"pT > 55")
      EffForData.GetXaxis().SetBinLabel(3,"NumPixHits > 2")
      EffForData.GetXaxis().SetBinLabel(4,"ValidFract > 0.9")
      EffForData.GetXaxis().SetBinLabel(5,"NumDeDx > 10")
      EffForData.GetXaxis().SetBinLabel(6,"Chi2oDOF < 5")
      EffForData.GetXaxis().SetBinLabel(7,"EoP < 0.3")
      EffForData.GetXaxis().SetBinLabel(8,"dz < 0.02")
      EffForData.GetXaxis().SetBinLabel(9,"dxy < 0.02")
      EffForData.GetXaxis().SetBinLabel(10,"MiniIso < 0.02")
      EffForData.GetXaxis().SetBinLabel(11,"Ih > 3.47")
      EffForData.GetXaxis().SetBinLabel(12,"ProbXY > 0.001")
      EffForData.GetXaxis().SetTitle("")
      
      legend =  ROOT.TLegend(.55,.80,.80,.89,"","brNDC")
      legend.SetTextFont(42)
      legend.SetTextSize(0.02)
      legend.SetBorderSize(1);
      legend.SetBorderSize(0);
      legend.SetLineColor(1);
      legend.SetLineStyle(1);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(1001);
      legend.AddEntry(EffForData,"DATA (2018C)","LP")
      

      EffForBcg = ROOT.TH1F("EffForBcg",";;Integral to (from) max (min) cut value /  Integral to inf",12,0.,12.)
      EffForBcg.Draw("SAMEHISTOTEXT00")
      EffForBcg.SetLineColor(3)
      EffForBcg.SetMarkerColor(3)
      EffForBcg.SetBinContent(1,BefPreSeffForEtaBcg)
      EffForBcg.SetBinContent(2,BefPreSeffForpTBcg)
      EffForBcg.SetBinContent(3,BefPreSeffForNumPixHitsBcg)
      EffForBcg.SetBinContent(4,BefPreSeffForValidFractBcg)
      EffForBcg.SetBinContent(5,BefPreSeffForNumDeDxBcg)
      EffForBcg.SetBinContent(6,BefPreSeffForChi2oDOFBcg)
      EffForBcg.SetBinContent(7,BefPreSeffForEoPBcg)
      EffForBcg.SetBinContent(8,BefPreSeffFordzBcg)
      EffForBcg.SetBinContent(9,BefPreSeffFordxyBcg)
      EffForBcg.SetBinContent(10,BefPreSeffForMiniIsoBcg)
      EffForBcg.SetBinContent(11,BefPreSeffForIhBcg)
      EffForBcg.SetBinContent(12,BefPreSeffForProbXY)
      legend.AddEntry(EffForBcg,"ALLBcg MC","LP")
      
      EffForSignal = ROOT.TH1F("EffForSignal",";;Integral to (from) max (min) cut value /  Integral to inf",12,0.,12.)
      EffForSignal.Draw("SAMEHISTOTEXT00")
      EffForSignal.SetLineColor(6)
      EffForSignal.SetMarkerColor(6)
      EffForSignal.SetBinContent(1,BefPreSeffForEtaSignal)
      EffForSignal.SetBinContent(2,BefPreSeffForpTSignal)
      EffForSignal.SetBinContent(3,BefPreSeffForNumPixHitsSignal)
      EffForSignal.SetBinContent(4,BefPreSeffForValidFractSignal)
      EffForSignal.SetBinContent(5,BefPreSeffForNumDeDxSignal)
      EffForSignal.SetBinContent(6,BefPreSeffForChi2oDOFSignal)
      EffForSignal.SetBinContent(7,BefPreSeffForEoPSignal)
      EffForSignal.SetBinContent(8,BefPreSeffFordzSignal)
      EffForSignal.SetBinContent(9,BefPreSeffFordxySignal)
      EffForSignal.SetBinContent(10,BefPreSeffForMiniIsoSignal)
      EffForSignal.SetBinContent(11,BefPreSeffForIhSignal)
      EffForSignal.SetBinContent(12,BefPreSeffForProbXY)
      legend.AddEntry(EffForSignal,"HSCP Gluino 1800 GeV","LP")
      
      tex4 = ROOT.TLatex(0.6,0.95,"Before pre-selection")
      tex4.SetNDC();
      tex4.SetTextFont(52);
      tex4.SetTextSize(0.045);
      tex4.SetLineWidth(2);
#      tex4 = ROOT.TLatex(0.6,0.95,"After N-1 selection")
      
      legend.Draw("SAME")
      tex2.Draw("SAME")
      tex3.Draw("SAME")
      tex4.Draw("SAME")
      tex5.Draw("SAME")
      
      canvas.SaveAs("Effs_CodeV"+codeVersion+"/"+"EffsInCutflowBefPreS.png")
      
      cstackedSummedBackgroundStringN1 = "cstackedSummedBackgroundStringN1"+str(keyname2)
      canvas = ROOT.TCanvas(cstackedSummedBackgroundStringN1, cstackedSummedBackgroundStringN1, 800,800)
      
      EffForData = ROOT.TH1F("EffForData",";;Integral to (from) max (min) cut value /  Integral to inf",12,0.,12.)
      EffForData.Draw("HISTTEXT00")
      EffForData.SetBinContent(1,N1effForEta)
      EffForData.SetBinContent(2,N1effForpT)
      EffForData.SetBinContent(3,N1effForNumPixHits)
      EffForData.SetBinContent(4,N1effForValidFract)
      EffForData.SetBinContent(5,N1effForNumDeDx)
      EffForData.SetBinContent(6,N1effForChi2oDOF)
      EffForData.SetBinContent(7,N1effForEoP)
      EffForData.SetBinContent(8,N1effFordz)
      EffForData.SetBinContent(9,N1effFordxy)
      EffForData.SetBinContent(10,N1effForMiniIso)
      EffForData.SetBinContent(11,N1effForIh)
      EffForData.SetBinContent(12,N1effForProbXY)

      
      EffForData.SetLineColor(1)
      EffForData.SetMarkerColor(1)
#      EffForData.SetMarkerStyle(20)
      EffForData.SetStats(0)
      EffForData.GetXaxis().SetTitle("")
      EffForData.GetYaxis().SetRangeUser(0.,1.3)
      EffForData.GetXaxis().SetBinLabel(1,"Eta < 1.0")
      EffForData.GetXaxis().SetBinLabel(2,"pT > 55")
      EffForData.GetXaxis().SetBinLabel(3,"NumPixHits > 2")
      EffForData.GetXaxis().SetBinLabel(4,"ValidFract > 0.9")
      EffForData.GetXaxis().SetBinLabel(5,"NumDeDx > 10")
      EffForData.GetXaxis().SetBinLabel(6,"Chi2oDOF < 5")
      EffForData.GetXaxis().SetBinLabel(7,"EoP < 0.3")
      EffForData.GetXaxis().SetBinLabel(8,"dz < 0.02")
      EffForData.GetXaxis().SetBinLabel(9,"dxy < 0.02")
      EffForData.GetXaxis().SetBinLabel(10,"MiniIso < 0.02")
      EffForData.GetXaxis().SetBinLabel(11,"Ih > 3.47")
      EffForData.GetXaxis().SetBinLabel(12,"ProbXY > 0.001")
      EffForData.GetXaxis().SetTitle("")
      
      legend =  ROOT.TLegend(.55,.80,.80,.89,"","brNDC")
      legend.SetTextFont(42)
      legend.SetTextSize(0.02)
      legend.SetBorderSize(1);
      legend.SetBorderSize(0);
      legend.SetLineColor(1);
      legend.SetLineStyle(1);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(1001);
      legend.AddEntry(EffForData,"DATA (2018C)","LP")
      

      EffForBcg = ROOT.TH1F("EffForBcg",";;Integral to (from) max (min) cut value /  Integral to inf",12,0.,12.)
      EffForBcg.Draw("SAMEHISTOTEXT00")
      EffForBcg.SetLineColor(3)
      EffForBcg.SetMarkerColor(3)
      EffForBcg.SetBinContent(1,N1effForEtaBcg)
      EffForBcg.SetBinContent(2,N1effForpTBcg)
      EffForBcg.SetBinContent(3,N1effForNumPixHitsBcg)
      EffForBcg.SetBinContent(4,N1effForValidFractBcg)
      EffForBcg.SetBinContent(5,N1effForNumDeDxBcg)
      EffForBcg.SetBinContent(6,N1effForChi2oDOFBcg)
      EffForBcg.SetBinContent(7,N1effForEoPBcg)
      EffForBcg.SetBinContent(8,N1effFordzBcg)
      EffForBcg.SetBinContent(9,N1effFordxyBcg)
      EffForBcg.SetBinContent(10,N1effForMiniIsoBcg)
      EffForBcg.SetBinContent(11,N1effForIhBcg)
      EffForBcg.SetBinContent(12,N1effForProbXY)
      legend.AddEntry(EffForBcg,"ALLBcg MC","LP")
      
      EffForSignal = ROOT.TH1F("EffForSignal",";;Integral to (from) max (min) cut value /  Integral to inf",12,0.,12.)
      EffForSignal.Draw("SAMEHISTOTEXT00")
      EffForSignal.SetLineColor(6)
      EffForSignal.SetMarkerColor(6)
      EffForSignal.SetBinContent(1,N1effForEtaSignal)
      EffForSignal.SetBinContent(2,N1effForpTSignal)
      EffForSignal.SetBinContent(3,N1effForNumPixHitsSignal)
      EffForSignal.SetBinContent(4,N1effForValidFractSignal)
      EffForSignal.SetBinContent(5,N1effForNumDeDxSignal)
      EffForSignal.SetBinContent(6,N1effForChi2oDOFSignal)
      EffForSignal.SetBinContent(7,N1effForEoPSignal)
      EffForSignal.SetBinContent(8,N1effFordzSignal)
      EffForSignal.SetBinContent(9,N1effFordxySignal)
      EffForSignal.SetBinContent(10,N1effForMiniIsoSignal)
      EffForSignal.SetBinContent(11,N1effForIhSignal)
      EffForSignal.SetBinContent(12,N1effForProbXY)
      legend.AddEntry(EffForSignal,"HSCP Gluino 1800 GeV","LP")
      
      tex4 = ROOT.TLatex(0.55,0.95,"After (N-1)+1 selection")
      tex4.SetNDC();
      tex4.SetTextFont(52);
      tex4.SetTextSize(0.045);
      tex4.SetLineWidth(2);
      
      legend.Draw("SAME")
      tex2.Draw("SAME")
      tex3.Draw("SAME")
      tex4.Draw("SAME")
      tex5.Draw("SAME")
      
      canvas.SaveAs("Effs_CodeV"+codeVersion+"/"+"EffsInCutflowN1.png")
      
####################################################################################################################
      # Punzi for BefPreS

      ROOT.gStyle.SetPaintTextFormat(".2g");
      cstackedSummedBackgroundStringBefPreSPunzi = "cstackedSummedBackgroundStringBefPreSPunzi"+str(keyname2)
      canvas = ROOT.TCanvas(cstackedSummedBackgroundStringBefPreSPunzi, cstackedSummedBackgroundStringBefPreSPunzi, 800,800)
      
      EffForData = ROOT.TH1F("EffForData",";;Punzi-significance",11,0.,11.)
      EffForData.Draw("HISTTEXT00")
      EffForData.SetBinContent(1,BefPreSeffForEtaPunzi)
      EffForData.SetBinContent(2,BefPreSeffForpTPunzi)
      EffForData.SetBinContent(3,BefPreSeffForNumPixHitsPunzi)
      EffForData.SetBinContent(4,BefPreSeffForValidFractPunzi)
      EffForData.SetBinContent(5,BefPreSeffForNumDeDxPunzi)
      EffForData.SetBinContent(6,BefPreSeffForChi2oDOFPunzi)
      EffForData.SetBinContent(7,BefPreSeffForEoPPunzi)
      EffForData.SetBinContent(8,BefPreSeffFordzPunzi)
      EffForData.SetBinContent(9,BefPreSeffFordxyPunzi)
      EffForData.SetBinContent(10,BefPreSeffForMiniIsoPunzi)
      EffForData.SetBinContent(11,BefPreSeffForIhPunzi)
#      EffForData.SetBinContent(12,BefPreSeffForProbXYPunzi)

      
      EffForData.SetLineColor(1)
      EffForData.SetMarkerColor(1)
      EffForData.SetStats(0)
      EffForData.GetXaxis().SetTitle("")
#      EffForData.GetYaxis().SetRangeUser(0.,1.3)
      EffForData.GetXaxis().SetBinLabel(1,"Eta < 1.0")
      EffForData.GetXaxis().SetBinLabel(2,"pT > 55")
      EffForData.GetXaxis().SetBinLabel(3,"NumPixHits > 2")
      EffForData.GetXaxis().SetBinLabel(4,"ValidFract > 0.9")
      EffForData.GetXaxis().SetBinLabel(5,"NumDeDx > 10")
      EffForData.GetXaxis().SetBinLabel(6,"Chi2oDOF < 5")
      EffForData.GetXaxis().SetBinLabel(7,"EoP < 0.3")
      EffForData.GetXaxis().SetBinLabel(8,"dz < 0.02")
      EffForData.GetXaxis().SetBinLabel(9,"dxy < 0.02")
      EffForData.GetXaxis().SetBinLabel(10,"MiniIso < 0.02")
      EffForData.GetXaxis().SetBinLabel(11,"Ih > 3.47")
#      EffForData.GetXaxis().SetBinLabel(12,"ProbXY > 0.001")
      EffForData.GetXaxis().SetTitle("")
      
      tex4 = ROOT.TLatex(0.55,0.95,"Before pre-selection")
      tex4.SetNDC();
      tex4.SetTextFont(52);
      tex4.SetTextSize(0.045);
      tex4.SetLineWidth(2);
      
      tex2.Draw("SAME")
      tex3.Draw("SAME")
      tex4.Draw("SAME")
      tex5.Draw("SAME")
      
      canvas.SaveAs("Effs_CodeV"+codeVersion+"/"+"EffsInCutflowBefPreSPunzi.png")
      
####################################################################################################################
      # Punzi for N-1
      
      ROOT.gStyle.SetPaintTextFormat(".4f");
      cstackedSummedBackgroundStringN1Punzi = "cstackedSummedBackgroundStringN1Punzi"+str(keyname2)
      canvas = ROOT.TCanvas(cstackedSummedBackgroundStringN1Punzi, cstackedSummedBackgroundStringN1Punzi, 800,800)
      
      EffForData = ROOT.TH1F("EffForData",";;Punzi-significance",11,0.,11.)
      EffForData.Draw("HISTTEXT00")
      EffForData.SetBinContent(1,N1effForEtaPunzi)
      EffForData.SetBinContent(2,N1effForpTPunzi)
      EffForData.SetBinContent(3,N1effForNumPixHitsPunzi)
      EffForData.SetBinContent(4,N1effForValidFractPunzi)
      EffForData.SetBinContent(5,N1effForNumDeDxPunzi)
      EffForData.SetBinContent(6,N1effForChi2oDOFPunzi)
      EffForData.SetBinContent(7,N1effForEoPPunzi)
      EffForData.SetBinContent(8,N1effFordzPunzi)
      EffForData.SetBinContent(9,N1effFordxyPunzi)
      EffForData.SetBinContent(10,N1effForMiniIsoPunzi)
      EffForData.SetBinContent(11,N1effForIhPunzi)
#      EffForData.SetBinContent(12,N1effForProbXYPunzi)

      
      EffForData.SetLineColor(1)
      EffForData.SetMarkerColor(1)
      EffForData.SetStats(0)
      EffForData.GetXaxis().SetTitle("")
#      EffForData.GetYaxis().SetRangeUser(0.,1.3)
      EffForData.GetXaxis().SetBinLabel(1,"Eta < 1.0")
      EffForData.GetXaxis().SetBinLabel(2,"pT > 55")
      EffForData.GetXaxis().SetBinLabel(3,"NumPixHits > 2")
      EffForData.GetXaxis().SetBinLabel(4,"ValidFract > 0.9")
      EffForData.GetXaxis().SetBinLabel(5,"NumDeDx > 10")
      EffForData.GetXaxis().SetBinLabel(6,"Chi2oDOF < 5")
      EffForData.GetXaxis().SetBinLabel(7,"EoP < 0.3")
      EffForData.GetXaxis().SetBinLabel(8,"dz < 0.02")
      EffForData.GetXaxis().SetBinLabel(9,"dxy < 0.02")
      EffForData.GetXaxis().SetBinLabel(10,"MiniIso < 0.02")
      EffForData.GetXaxis().SetBinLabel(11,"Ih > 3.47")
#      EffForData.GetXaxis().SetBinLabel(12,"ProbXY > 0.001")
      EffForData.GetXaxis().SetTitle("")
      
      tex4 = ROOT.TLatex(0.55,0.95,"After (N-1)+1 selection")
      tex4.SetNDC();
      tex4.SetTextFont(52);
      tex4.SetTextSize(0.045);
      tex4.SetLineWidth(2);
      
      tex2.Draw("SAME")
      tex3.Draw("SAME")
      tex4.Draw("SAME")
      tex5.Draw("SAME")
      
      canvas.SaveAs("Effs_CodeV"+codeVersion+"/"+"EffsInCutflowN1Punzi.png")

