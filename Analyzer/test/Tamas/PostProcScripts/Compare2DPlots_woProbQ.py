import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

codeVersion = sys.argv[1]
#BinNumber = sys.argv[2]
#bin = int(BinNumber)
bin = 3

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

SingleMuon = ROOT.TFile.Open("crab_Analysis_2018_SingleMuon_Run2018_woProbQ_CodeV"+codeVersion+"_v1.root")

AllTTbar = ROOT.TFile.Open("crab_Analysis_2018_AllTTbar_woProbQ_CodeV"+codeVersion+"_v1.root")
AllWJets = ROOT.TFile.Open("crab_Analysis_2018_AllWJets_woProbQ_CodeV"+codeVersion+"_v1.root")
AllQCD   = ROOT.TFile.Open("crab_Analysis_2018_AllQCD_woProbQ_CodeV"+codeVersion+"_v1.root")

SelectedSignalSamples1 = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-1800_woProbQ_CodeV"+codeVersion+"_v1.root")
SelectedSignalSamples2 = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-2400_woProbQ_CodeV"+codeVersion+"_v1.root")

  
if os.path.exists(os.path.dirname("2DComparrison_CodeV"+codeVersion)):
  os.system("mkdir 2DComparrison_woProbQ_CodeV"+codeVersion)

dirs = []
for i in range(0, AllTTbar.GetListOfKeys().GetEntries()):
  dirname = AllTTbar.GetListOfKeys().At(i).GetName()
  curr_dir = AllTTbar.GetDirectory(dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = AllTTbar.GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          newname = dirname + "/" + keyname+ "/" + keyname2
          obj = AllTTbar.Get(newname)
#          print(keyname2)
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if not obj.InheritsFrom("TObject"): continue

          if (obj.GetEntries() == 0 ) : continue
          
          if not (obj.ClassName() == "TH2F") : continue
          
          if ("CutFlow" in keyname2) : continue
          if ("Gen" in keyname2) : continue
          if not ("PostPre" in keyname2) : continue
          
          SingleMuonHisto = SingleMuon.Get(newname)
          SingleMuonHisto.SetMarkerStyle(20)
          SingleMuonHisto.SetMarkerColor(1)
          SingleMuonHisto.SetLineColor(1)
        
          AllTTbarHisto = AllTTbar.Get(newname)
          AllTTbarHisto.SetMarkerStyle(20)
          AllTTbarHisto.SetMarkerColor(2)
          AllTTbarHisto.SetLineColor(2)
          
          AllWJetsHisto = AllWJets.Get(newname)
          AllWJetsHisto.SetMarkerStyle(20)
          AllWJetsHisto.SetMarkerColor(3)
          AllWJetsHisto.SetLineColor(3)
          
          AllQCDHisto = AllQCD.Get(newname)
          AllQCDHisto.SetMarkerStyle(20)
          AllQCDHisto.SetMarkerColor(4)
          AllQCDHisto.SetLineColor(4)
          
          SelectedSignalSamples1Histo = SelectedSignalSamples1.Get(newname)
          SelectedSignalSamples1Histo.SetMarkerStyle(20)
          SelectedSignalSamples1Histo.SetMarkerColor(6)
          SelectedSignalSamples1Histo.SetLineColor(6)
          
          SelectedSignalSamples2Histo = SelectedSignalSamples2.Get(newname)
          SelectedSignalSamples2Histo.SetMarkerStyle(20)
          SelectedSignalSamples2Histo.SetMarkerColor(7)
          SelectedSignalSamples2Histo.SetLineColor(7)
          
          legend =  ROOT.TLegend(.55,.70,.80,.89,"","brNDC")
          legend.SetTextFont(42)
          legend.SetTextSize(0.02)
          legend.SetBorderSize(1);
          legend.SetBorderSize(0);
          legend.SetLineColor(1);
          legend.SetLineStyle(1);
          legend.SetLineWidth(1);
          legend.SetFillColor(0);
          legend.SetFillStyle(1001);
          if (keyname2.find("Vs")==-1) :
            axisXTitle = keyname2[keyname2.find("_")+1:]
            axisYTitle = "Tracks/bin"
          else :
            axisXTitle = keyname2[keyname2.find("_")+1:keyname2.find("Vs")]
            axisYTitle = keyname2[keyname2.find("Vs")+2:]
          
          legend.AddEntry(AllTTbarHisto,"TTBar","LP")
          legend.AddEntry(AllWJetsHisto,"WJets","LP")
          legend.AddEntry(AllQCDHisto,"mu-QCD p_{T}={50,inf}","LP")
#          legend.AddEntry(SingleMuonHisto,"SingleMuon 2018C","LP")
          legend.AddEntry(SelectedSignalSamples1Histo, "HSCPgluino M=1800 GeV", "LP")
          legend.AddEntry(SelectedSignalSamples2Histo, "HSCPgluino M=2400 GeV", "LP")


          
          tex2 = ROOT.TLatex(0.13,0.94,"CMS");
          #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

          tex3 = ROOT.TLatex(0.27,0.94,"Internal"); # for square plots
          #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
          #tex3 = ROOT.TLatex(0.28,0.94,"Internal");
          tex3.SetNDC();
          tex3.SetTextFont(52);
          tex3.SetTextSize(0.0485);
          tex3.SetLineWidth(2);
          
          tex4 = ROOT.TLatex()
          if ("PrePreS" in keyname2) :
            tex4 = ROOT.TLatex(0.75,0.94,"Before pre-selection")
          elif ("N1" in keyname2) :
            tex4 = ROOT.TLatex(0.75,0.94,"After N-1 selection")
          elif ("PostPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.75,0.94,"After pre-selection")
            
          tex4.SetNDC();
          tex4.SetTextFont(52);
          tex4.SetTextSize(0.0285);
          tex4.SetLineWidth(2);
          
          tex5 = ROOT.TLatex(0.07,0.04,"Code version: "+codeVersion);
          tex5.SetNDC();
          tex5.SetTextFont(52);
          tex5.SetTextSize(0.0185);
          tex5.SetLineWidth(2);
#          max = numpy.maximum(stackedSummedBackground.GetMaximum(),SingleMuonHisto.GetMaximum())

          cstackedSummedBackgroundString = 'cstackedSummedBackground'+str(j)
          cstackedSummedBackground = ROOT.TCanvas(cstackedSummedBackgroundString, cstackedSummedBackgroundString, 800,800)

          AllTTbarHisto.Draw("P")
          AllTTbarHisto.SetStats(0)
          AllTTbarHisto.SetTitle("")
          AllTTbarHisto.GetXaxis().SetTitleSize(0.05)
          AllTTbarHisto.GetXaxis().SetTitleOffset(1)
          AllTTbarHisto.GetXaxis().SetTitle(axisXTitle)
          AllTTbarHisto.GetYaxis().SetTitleSize(0.05)
#          AllTTbarHisto.GetYaxis().SetTitleOffset(1)
          AllTTbarHisto.GetYaxis().SetTitle(axisYTitle)
          AllWJetsHisto.Draw("SAMEP")
          AllQCDHisto.Draw("SAMEP")
#          SingleMuonHisto.Draw("SAMEP")
          SelectedSignalSamples1Histo.Draw("SAMEP")
          SelectedSignalSamples2Histo.Draw("SAMEP")

          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")

          cstackedSummedBackground.SaveAs("2DComparrison_CodeV"+codeVersion+"/"+keyname2+".png")
          
#---------------------------------------------------------------------------------------------
#          # now let's do everything again but on a log Y scale
#          cstackedSummedBackgroundLogString = 'cstackedSummedBackgroundLog'+str(j)
#          cstackedSummedBackground = ROOT.TCanvas(cstackedSummedBackgroundLogString, cstackedSummedBackgroundLogString, 800,800)
#          cstackedSummedBackground.SetLogy()
#
#          if ("PrePreS" in keyname2 or "PostPreS" in keyname2 or "N1" in keyname2) :
#            stackedSummedBackground.SetMinimum(0.0001)
#          else:
#            stackedSummedBackground.SetMinimum(0.000000000000000001)
#
#          if ("CutFlow" in keyname2 or "pfType" in keyname2):
#            stackedNormSummedBackground.Draw("HISTO")
#            stackedNormSummedBackground.SetTitle("")
#            stackedNormSummedBackground.GetXaxis().SetTitleSize(0.05)
#            stackedNormSummedBackground.GetXaxis().SetTitle("")
#            stackedNormSummedBackground.GetYaxis().SetRangeUser(0,1.5)
#            stackedNormSummedBackground.SetMaximum(12)
#          else :
#            stackedSummedBackground.Draw("HISTO")
#            stackedSummedBackground.GetXaxis().SetTitleSize(0.05)
#            stackedSummedBackground.GetXaxis().SetTitleOffset(1)
#            stackedSummedBackground.SetMaximum(max*1000)
#          SingleMuonHisto.Draw("SAMEP")
#          SelectedSignalSamples1Histo.Draw("SAME")
#          SelectedSignalSamples2Histo.Draw("SAME")
#          legend.Draw("SAME")
#          tex2.Draw("SAME")
#          tex3.Draw("SAME")
#          tex4.Draw("SAME")
#          tex5.Draw("SAME")
#
#          cstackedSummedBackground.SaveAs("2DComparrison_CodeV"+codeVersion+"/"+keyname2+"_log.png")
