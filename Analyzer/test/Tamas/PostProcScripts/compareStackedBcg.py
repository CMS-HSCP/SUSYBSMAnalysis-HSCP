import ROOT, sys, os, time, re, numpy
import numpy as np
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion BinNumber")
(opt,args) = parser.parse_args()

codeVersion = sys.argv[1]
BinNumber = sys.argv[2] if (len(sys.argv)==1) else 3
bin = int(BinNumber)
# bin 3: pt>60 and I_as > 0.05
# bin 25: pt>65 and I_as > 0.175
# bin 28: pt>65 and I_as > 0.3

blind = True

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

BackgroundSamples = [
"crab_Analysis_2018_AllTTbar_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_AllWJets_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_AllQCD_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_AllZToMuMu_CodeV"+codeVersion+"_v1.root",
]

SingleMuon = ROOT.TFile.Open("crab_Analysis_SingleMuon_RunPhase1_CodeV"+codeVersion+"_v1.root")
#SingleMuon = ROOT.TFile.Open("crab_Analysis_SingleMuon_Run2018C_CodeV"+codeVersion+"_v1.root")

SelectedSignalSamples1 = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root")
SelectedSignalSamples2 = ROOT.TFile.Open("crab_Analysis_2018_HSCPpairStau_M-557_CodeV"+codeVersion+"_v1.root")
#SelectedSignalSamples2 = ROOT.TFile.Open("crab_Analysis_2018_HSCPpairStau_M-1599_CodeV"+codeVersion+"_v1.root")
#SelectedSignalSamples2 = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-2200_CodeV"+codeVersion+"_v1.root")

bckArray = []
for sample in BackgroundSamples:
  bckArray.append(ROOT.TFile.Open(sample))

name = "StackedComparrison_CodeV"+codeVersion +"_Bin" + str(bin)+"/A.png"
print(os.path.dirname(name))
if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))

dirs = []
for i in range(0, bckArray[0].GetListOfKeys().GetEntries()):
  dirname = bckArray[0].GetListOfKeys().At(i).GetName()
  curr_dir = bckArray[0].GetDirectory(dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = bckArray[0].GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          newname = dirname + "/" + keyname+ "/" + keyname2
          obj = bckArray[0].Get(newname)
          if not (obj) : continue
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if not obj.InheritsFrom("TObject"): continue
          if any(substring in keyname2 for substring in ["_regionA_", "_regionB_", "_regionC_", "_regionD_", "VR1_Mass",  "VR2_Mass", "VR3_Mass", "SR1_Mass",  "SR2_Mass", "SR3_Mass"]) : continue
          
          if not ("MiniRelIsoAll" in keyname2) : continue

          if (obj.GetEntries() == 0 ) : continue
          stackedSummedBackground = ROOT.THStack("stackedSummedBackground","")
          stackedNormSummedBackground = ROOT.THStack("stackedNormSummedBackground","")
          
          SingleMuonHisto = SingleMuon.Get(newname)
          SelectedSignalSamples1Histo = SelectedSignalSamples1.Get(newname)
          SelectedSignalSamples2Histo = SelectedSignalSamples2.Get(newname)
          
          if (obj.ClassName() == "TH2F") :
            xAxisTitle = SingleMuon.Get(newname).GetYaxis().GetTitle()
            SingleMuonHisto = SingleMuon.Get(newname).ProjectionY(str(keyname2)+"Data_ProjY",bin,bin,"e")
            SingleMuonHisto.GetXaxis().SetTitle(xAxisTitle)
            SelectedSignalSamples1Histo = SelectedSignalSamples1.Get(newname).ProjectionY(str(keyname2)+"Signal1_ProjY",bin,bin,"e")
            SelectedSignalSamples2Histo = SelectedSignalSamples2.Get(newname).ProjectionY(str(keyname2)+"Signal2_ProjY",bin,bin,"e")
          
          SingleMuonHisto.SetMarkerStyle(20)
          SingleMuonHisto.SetMarkerColor(1)
          SingleMuonHisto.SetLineColor(1)
          
          SelectedSignalSamples1Histo.SetMarkerColor(6)
          SelectedSignalSamples1Histo.SetLineColor(6)
          SelectedSignalSamples1Histo.SetMarkerStyle(22)
          
          SelectedSignalSamples2Histo.SetMarkerColor(7)
          SelectedSignalSamples2Histo.SetLineColor(7)
          SelectedSignalSamples2Histo.SetMarkerStyle(22)
          
          legend =  ROOT.TLegend(.5,.70,.80,.89,"","brNDC")
          legend.SetTextFont(42)
          legend.SetTextSize(0.02)
          legend.SetBorderSize(1);
          legend.SetBorderSize(0);
          legend.SetLineColor(1);
          legend.SetLineStyle(1);
          legend.SetLineWidth(1);
          legend.SetFillColor(0);
          legend.SetFillStyle(1001);
#          if (keyname2.find("Per")==-1) :
#            axisXTitle = keyname2[keyname2.find("_")+1:]
#            axisYTitle = "Tracks/bin"
#          else :
#            axisXTitle = keyname2[keyname2.find("_")+1:keyname2.find("Per")]
#            axisYTitle = keyname2[keyname2.find("Per")+3:]

          if ((obj.ClassName() == "TH1F") or (obj.ClassName() == "TH2F" and (("PostS" in keyname2) or ("Mass" in keyname2)) and not "Vs" in keyname2) or ("K_and_C_Kin" in keyname2)) :
            # array to contain a specific (keyname2) histogram for all samples
#            print(keyname2)
            histoArray = []
            nEventsPostTrigArray = []
            for index,fileIn in enumerate(bckArray):
              if (obj.ClassName() == "TH1F" or obj.ClassName() == "TH1D") :
                histo = fileIn.Get(newname)
                if not histo : continue
                xAxisTitle = histo.GetXaxis().GetTitle()
                yAxisTitle = histo.GetYaxis().GetTitle()
                if (yAxisTitle == "") : yAxisTitle = "Tracks / bin"
                if ("Ias" in keyname2) : xAxisTitle = "G_{i}"
                if ("Pt" in keyname2) and not ("PtErr" in keyname2) : xAxisTitle = "p_{T} [GeV]"
                if ("Chi2oNdof" in keyname2) : xAxisTitle = "#chi^{2} / N_{dof}"
                if ("Dxy" in keyname2) : xAxisTitle = "d_{xy} [cm]"
              if (obj.ClassName() == "TH2F") :
                histo = obj.ProjectionY(str(newname)+"ProjY"+str(index),bin,bin,"e")
                if not histo : continue
                xAxisTitle = histo.GetYaxis().GetTitle()
                yAxisTitle = histo.GetZaxis().GetTitle()
              histo.GetXaxis().SetTitle(xAxisTitle)
              histoArray.append(histo)
              stackedSummedBackground.Add(histo)
              if ((index==2)) :
                histo.SetLineColor(2)
                histo.SetFillColor(2)
                histo.SetMarkerColor(2)
                legend.AddEntry(histo,"mu-QCD p_{T}={50,inf}","LP")
              elif ((index==1)) :
                histo.SetLineColor(3)
                histo.SetFillColor(3)
                histo.SetMarkerColor(4)
                legend.AddEntry(histo,"WJets","LP")
              elif ((index==0)) :
                histo.SetLineColor(4)
                histo.SetFillColor(4)
                histo.SetMarkerColor(4)
                legend.AddEntry(histo,"TTBar","LP")
              elif ((index==3)) :
                histo.SetLineColor(28)
                histo.SetFillColor(28)
                histo.SetMarkerColor(28)
                legend.AddEntry(histo,"ZToMuMu","LP")
          else :
            continue

          cstackedSummedBackgroundString = 'cstackedSummedBackground'+str(j)
          cstackedSummedBackground = ROOT.TCanvas(cstackedSummedBackgroundString, cstackedSummedBackgroundString, 800,800)
          #  convert stacks to (summed) histos
          stackedSummedBackground.Draw()
          OFbin = stackedSummedBackground.GetStack().Last().GetNbinsX()+1
          if any(substring in keyname2 for substring in ["Ias", "ProbQ", "Beta", "CutFlow", "Qual", "Eta", "MiniTkIso", "closestPfJet"]): #"MiniRelIsoAll"
            OFbin = stackedSummedBackground.GetStack().Last().GetNbinsX()
          UFbin = 1
          if any(substring in keyname2 for substring in ["Dxy", "Dz"]):
            UFbin = 0
          stackedSummedBackground.GetXaxis().SetRange(UFbin,OFbin)
#          stackedSummedBackground.GetXaxis().SetTitle(xAxisTitle)
          stackedSummedBackground.GetYaxis().SetTitle(yAxisTitle)
          stackedSummedBackground.GetXaxis().SetTitle(xAxisTitle)
          stackedSummedBackgroundTmp = stackedSummedBackground.GetStack().Last()
          
          for hist in stackedSummedBackground.GetStack():
            hist.GetXaxis().SetRange(UFbin,OFbin)
          
          SelectedSignalSamples1Histo.GetXaxis().SetRange(UFbin,OFbin)
          SelectedSignalSamples2Histo.GetXaxis().SetRange(UFbin,OFbin)
          SingleMuonHisto.GetXaxis().SetRange(UFbin,OFbin)
#          stackedSummedBackgroundTmp
          for index,fileIn in enumerate(bckArray):
            normHisto = fileIn.Get(newname)
            if not normHisto : continue
            max = stackedSummedBackgroundTmp.GetMaximum()
            if (max==0) : continue
            normHisto.Scale(1/max)
            stackedNormSummedBackground.Add(normHisto)
          stackedNormSummedBackground.Draw()
          
          if (keyname2== "EventCutFlow") :
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
            if (blind) :
              SingleMuonHisto.SetBinContent(18,0)
              SingleMuonHisto.SetBinContent(19,0)
              SingleMuonHisto.SetBinContent(20,0)
            if (SelectedSignalSamples1Histo.GetMaximum() == 0 ): continue
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
            obj.SetStats(0)
            ROOT.gStyle.SetPaintTextFormat(".2g");
            obj.Draw("HISTOTEXT00")
            tex2.Draw("SAME")
            tex3.Draw("SAME")
            tex4.Draw("SAME")
            tex5.Draw("SAME")
            overFlowText.Draw("SAME")
            overFlowLine.DrawLine(overFlowLocX,stackedSummedBackgroundTmp.GetMinimum(),overFlowLocX,max*3)
          elif (keyname2== "CutFlow") :
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
            if (blind) :
              SingleMuonHisto.SetBinContent(18,-1)
              SingleMuonHisto.SetBinContent(19,-1)
              SingleMuonHisto.SetBinContent(20,-1)
          elif (keyname2== "ProbQ") :
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
#            if (blind) :
#              SingleMuonHisto.SetBinContent(18,-1)
#              SingleMuonHisto.SetBinContent(19,-1)
#              SingleMuonHisto.SetBinContent(20,-1)
#              SingleMuonHisto.SetBinContent(21,-1)
            if (SelectedSignalSamples1Histo.GetMaximum() == 0 ): continue
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
          elif (keyname2== "CutFlowReverse") :
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
            if (SelectedSignalSamples1Histo.GetMaximum() == 0 ): continue
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
          elif (keyname2 == "ErrorHisto") :
            if (SelectedSignalSamples1Histo.GetMaximum() == 0 ): continue
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
            stackedNormSummedBackground.SetMaximum(1.4)
          elif ("pfType" in keyname2) :
            if (SelectedSignalSamples1Histo.GetMaximum() == 0 ): continue
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
          legend.AddEntry(SingleMuonHisto,"SingleMuon 2017-2018","LP")
#          legend.AddEntry(SingleMuonHisto,"SingleMuon 2018C","LP")
          legend.AddEntry(SelectedSignalSamples1Histo, "HSCPgluino M=1800 GeV", "LP")
          legend.AddEntry(SelectedSignalSamples2Histo, "HSCPstau M=557 GeV", "LP")
#          legend.AddEntry(METHisto,"MET-EraC","LP")
          
          tex2 = ROOT.TLatex(0.13,0.94,"CMS");
          #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

          tex3 = ROOT.TLatex(0.27,0.94,""); # for square plots
#          tex3 = ROOT.TLatex(0.27,0.94,"Internal"); # for square plots
          #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
          #tex3 = ROOT.TLatex(0.28,0.94,"Internal");
          tex3.SetNDC();
          tex3.SetTextFont(52);
          tex3.SetTextSize(0.0485);
          tex3.SetLineWidth(2);
          
          tex4 = ROOT.TLatex()

          if ("Calibration" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"Calibration selection")
          elif ("BefPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"Before preselection")
          elif ("N1" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After N-1 selection")
          elif ("PostPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After preselection")
          elif ("PostS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After selection")
            
          tex4.SetNDC();
          tex4.SetTextFont(52);
          tex4.SetTextSize(0.045);
          tex4.SetLineWidth(2);
          
          tex5 = ROOT.TLatex(0.07,0.04,"Code version: "+codeVersion);
          tex5.SetNDC();
          tex5.SetTextFont(52);
          tex5.SetTextSize(0.0185);
          tex5.SetLineWidth(2);
          max = numpy.maximum(stackedSummedBackground.GetMaximum(),SingleMuonHisto.GetMaximum())
          
          overFlowText = ROOT.TLatex(0.836,0.15,"OF");
#          overFlowText.SetTextAngle(-30)
          overFlowText.SetNDC();
          overFlowText.SetTextFont(52);
          overFlowText.SetTextSize(0.01);
          overFlowText.SetLineWidth(2);
          
          overFlowLine = ROOT.TLine();
          overFlowLine.SetLineWidth(2);
          overFlowLine.SetLineStyle(ROOT.kDashed);
          
          overFlowLine2 = ROOT.TLine();
          overFlowLine2.SetLineWidth(2);
          overFlowLine2.SetLineStyle(ROOT.kDashed);
          
          binWidth = histo.GetXaxis().GetBinWidth(1)
          firstBinLocXCent = histo.GetXaxis().GetBinCenter(1)
          overFlowLocXCent = histo.GetXaxis().GetBinCenter(histo.GetNbinsX()+ 1)
          overFlowLocX = (overFlowLocXCent - (binWidth)/2)

          if ("CutFlow" in keyname2 or "pfType" in keyname2):
            stackedNormSummedBackground.Draw("HISTO")
            stackedNormSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedNormSummedBackground.GetYaxis().SetRangeUser(0,1.5)
            stackedNormSummedBackground.SetMaximum(1.3)
          else :
            stackedSummedBackground.Draw("HISTO")
            stackedSummedBackground.GetXaxis().SetTitleOffset(1.1)
            stackedSummedBackground.GetXaxis().SetTitle(xAxisTitle)
            stackedSummedBackground.SetMaximum(max*1.4)
            stackedSummedBackground.SetMinimum(0.0)
          if not ("Mass" in keyname2 or "Ias" in keyname2 or "SR2PASS" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          if ("SR2FAIL" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          if ("CR" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          if ("K_and_C_Kin_Mass" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          SelectedSignalSamples1Histo.Draw("SAME")
          SelectedSignalSamples2Histo.Draw("SAME")
           
          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")
          overFlowText.Draw("SAME")
          overFlowLine.DrawLine(overFlowLocX,stackedSummedBackgroundTmp.GetMinimum(),overFlowLocX,max*1.4)
          cstackedSummedBackground.SaveAs("StackedComparrison_CodeV"+codeVersion+"_Bin"+str(bin)+"/"+keyname2+".png")
          
#---------------------------------------------------------------------------------------------
          # now let's do everything again but on a log Y scale
          cstackedSummedBackgroundLogString = 'cstackedSummedBackgroundLog'+str(j)
          cstackedSummedBackground = ROOT.TCanvas(cstackedSummedBackgroundLogString, cstackedSummedBackgroundLogString, 800,800)
          cstackedSummedBackground.SetLogy()
          
          stackedSummedBackground.SetMinimum(0.0001)
            
          if ("CutFlow" in keyname2 or "pfType" in keyname2):
            stackedNormSummedBackground.Draw("HISTO")
            stackedNormSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedNormSummedBackground.GetYaxis().SetRangeUser(0,1.5)
            stackedNormSummedBackground.SetMaximum(200)
            stackedNormSummedBackground.SetMinimum(0.0000001)
          else :
            stackedSummedBackground.Draw("HISTO")
#            stackedSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedSummedBackground.GetXaxis().SetTitleOffset(1.1)
            stackedSummedBackground.GetXaxis().SetTitle(xAxisTitle)
            stackedSummedBackground.SetMaximum(max*10000)
          if not ("Mass" in keyname2 or "Ias" in keyname2 or "SR2PASS" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          if ("CR" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          if ("K_and_C_Kin_Mass" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          if ("SR2FAIL" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
          SelectedSignalSamples1Histo.Draw("SAME")
          SelectedSignalSamples2Histo.Draw("SAME")
          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")
          overFlowText.Draw("SAME")
          overFlowLine.DrawLine(overFlowLocX,0.00001,overFlowLocX,max*40000)
           # max was max*40000 earlier, or stackedSummedBackground.GetMaximum()
          
          cstackedSummedBackground.SaveAs("StackedComparrison_CodeV"+codeVersion+"_Bin"+str(bin)+"/"+keyname2+"_log.png")
          
          
#---------------------------------------------------------------------------------------------
          # now let's do everything again but on a log Y scale normalized to the area
          cstackedSummedBackgroundNormString = 'cstackedSummedBackgroundNorm'+str(j)
          cstackedSummedBackgroundNorm = ROOT.TCanvas(cstackedSummedBackgroundNormString, cstackedSummedBackgroundNormString, 800,800)
          cstackedSummedBackgroundNorm.SetLogy(0)

          normFactToArea = 0 #
          normFactToArea2 = 0
          for hist in stackedSummedBackground.GetStack():
            normFactToArea += hist.Integral(1,hist.GetNbinsX()+1)
            
          if (normFactToArea == 0) :
            print("No data")
            normFactToArea = 1
            
          for hist in stackedSummedBackground.GetStack():
            normFactToArea2 += hist.Integral(1,hist.GetNbinsX()+1)
            
          stackedSummedBackgroundNormToArea = ROOT.THStack("stackedSummedBackgroundNormToArea","")
          for hist in stackedSummedBackground.GetStack():
            hist.GetXaxis().SetRange(UFbin,OFbin)
            normHistoToArea = hist.Clone()
            normHistoToArea.Scale(1/normFactToArea)
            stackedSummedBackgroundNormToArea.Add(normHistoToArea)
            
          if (SingleMuonHisto.Integral() > 0) :
            SingleMuonHisto.Scale(1/SingleMuonHisto.Integral(1,SingleMuonHisto.GetNbinsX()+1))
          if (SelectedSignalSamples1Histo.Integral() > 0) :
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.Integral(1,SelectedSignalSamples1Histo.GetNbinsX()+1))
          if (SelectedSignalSamples2Histo.Integral() > 0) :
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.Integral(1,SelectedSignalSamples2Histo.GetNbinsX()+1))

          
          if ("PrePreS" in keyname2 or "PostPreS" in keyname2 or "N1" in keyname2) :
            stackedSummedBackground.SetMinimum(0.0001)
          if ("CutFlow" in keyname2 or "pfType" in keyname2):
            stackedSummedBackgroundNormToArea.Draw("HISTO")
            stackedSummedBackgroundNormToArea.GetXaxis().SetTitleSize(0.05)
            stackedSummedBackgroundNormToArea.GetYaxis().SetRangeUser(0,1.5)
            stackedSummedBackgroundNormToArea.SetMaximum(1)
            stackedSummedBackgroundNormToArea.SetMinimum(0.00000000001)
          else :
            stackedSummedBackgroundNormToArea.Draw("HISTO")
            stackedSummedBackgroundNormToArea.GetYaxis().SetTitle("Normalized "+stackedSummedBackground.GetYaxis().GetTitle())
            stackedSummedBackgroundNormToArea.GetXaxis().SetTitleOffset(1.1)
            stackedSummedBackgroundNormToArea.GetXaxis().SetRange(UFbin,OFbin)
            stackedSummedBackgroundNormToArea.GetXaxis().SetTitle(xAxisTitle)
            max1 = np.maximum(stackedSummedBackgroundNormToArea.GetMaximum(),SelectedSignalSamples1Histo.GetMaximum())
            max2 = np.maximum(SelectedSignalSamples1Histo.GetMaximum(),SingleMuonHisto.GetMaximum())
            max = np.maximum(max1,max2)
            stackedSummedBackgroundNormToArea.SetMaximum(max*1.2)
            stackedSummedBackgroundNormToArea.SetMinimum(0.00000000001)
          SelectedSignalSamples1Histo.Draw("SAME")
          SelectedSignalSamples2Histo.Draw("SAME")
          if not ("Mass" in keyname2 or "Ias" in keyname2 or "SR2PASS" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
            SingleMuonHisto.SetMarkerStyle(24)
          if ("K_and_C_Kin_Mass" in keyname2 or "SR2FAIL" in keyname2 or "CR" in keyname2) :
            SingleMuonHisto.Draw("SAMEP")
            SingleMuonHisto.SetMarkerStyle(24)

          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")
          overFlowText.Draw("SAME")
          overFlowLine.DrawLine(overFlowLocX,stackedSummedBackgroundNormToArea.GetMinimum(),overFlowLocX,max*1.2)
          cstackedSummedBackgroundNorm.SaveAs("StackedComparrison_CodeV"+codeVersion+"_Bin"+str(bin)+"/"+keyname2+"_Norm.png")

#################################### same again but on log scale ####################################
          
          cstackedSummedBackgroundLogNormString = 'cstackedSummedBackgroundLogNorm'+str(j)
          cstackedSummedBackgroundLogNorm = ROOT.TCanvas(cstackedSummedBackgroundLogNormString, cstackedSummedBackgroundLogNormString, 800,800)
          cstackedSummedBackgroundLogNorm.SetLogy(1)
#          stackedSummedBackgroundNormToArea.SetMaximum(1000)

          if ("PrePreS" in keyname2 or "PostPreS" in keyname2 or "N1" in keyname2) :
            stackedSummedBackground.SetMinimum(0.0001)
          if ("CutFlow" in keyname2 or "pfType" in keyname2):
            stackedSummedBackgroundNormToArea.Draw("HISTO")
            stackedSummedBackgroundNormToArea.GetXaxis().SetTitleSize(0.05)
            stackedSummedBackgroundNormToArea.GetYaxis().SetRangeUser(0,1.5)
            stackedSummedBackgroundNormToArea.SetMaximum(1)
            stackedSummedBackgroundNormToArea.SetMinimum(0.00000000001)
          else :
            stackedSummedBackgroundNormToArea.Draw("HISTO")
#            stackedSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedSummedBackgroundNormToArea.GetYaxis().SetTitle("Normalized "+stackedSummedBackground.GetYaxis().GetTitle())
            stackedSummedBackgroundNormToArea.GetXaxis().SetTitleOffset(1.1)
            stackedSummedBackgroundNormToArea.GetXaxis().SetRange(UFbin,OFbin)
            stackedSummedBackgroundNormToArea.GetXaxis().SetTitle(xAxisTitle)
            stackedSummedBackgroundNormToArea.SetMaximum(max*10000)
            stackedSummedBackgroundNormToArea.SetMinimum(0.00000000001)
          SelectedSignalSamples1Histo.Draw("SAME")
          SelectedSignalSamples2Histo.Draw("SAME")
          if not any(substring in keyname2 for substring in ["Mass", "Ias", "SR2PASS"]):
            SingleMuonHisto.Draw("SAMEP")
            SingleMuonHisto.SetMarkerStyle(24)
          if any(substring in keyname2 for substring in ["K_and_C_Kin_Mass", "SR2FAIL", "CR"]):
            SingleMuonHisto.Draw("SAMEP")
            SingleMuonHisto.SetMarkerStyle(24)

          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")
          overFlowText.Draw("SAME")
          
          overFlowLine2.DrawLine(overFlowLocX,stackedSummedBackgroundNormToArea.GetMinimum(),overFlowLocX,max*40000)
          cstackedSummedBackgroundLogNorm.SaveAs("StackedComparrison_CodeV"+codeVersion+"_Bin"+str(bin)+"/"+keyname2+"_LogNorm.png")

os.system("cp forWebpage/* StackedComparrison_CodeV"+codeVersion+"_Bin"+str(bin)+"/.")
os.system("cp forWebpage/.htaccess StackedComparrison_CodeV"+codeVersion+"_Bin"+str(bin)+"/.")
print("scp -r StackedComparrison_CodeV"+ codeVersion+"_Bin"+str(bin) + " tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/.")
