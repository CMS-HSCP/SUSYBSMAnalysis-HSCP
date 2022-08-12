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

BackgroundSamples = [
"crab_Analysis_2018_AllQCD_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_AllWJets_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_AllTTbar_wProbQ_CodeV"+codeVersion+"_v1.root",
]

SingleMuon = ROOT.TFile.Open("crab_Analysis_SingleMuon_Run2018_wProbQ_CodeV"+codeVersion+"_v1.root")

SelectedSignalSamples1 = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-1800_wProbQ_CodeV"+codeVersion+"_v1.root")
SelectedSignalSamples2 = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-2400_wProbQ_CodeV"+codeVersion+"_v1.root")

bckArray = []
for sample in BackgroundSamples:
  bckArray.append(ROOT.TFile.Open(sample))
  
if not os.path.exists(os.path.dirname("StackedComparrison_wProbQ_CodeV"+codeVersion)) :
  os.system("mkdir StackedComparrison_wProbQ_CodeV"+codeVersion)

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
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if not obj.InheritsFrom("TObject"): continue

          if (obj.GetEntries() == 0 ) : continue
          stackedSummedBackground = ROOT.THStack("stackedSummedBackground","")
          stackedNormSummedBackground = ROOT.THStack("stackedNormSummedBackground","")
          
          SingleMuonHisto = SingleMuon.Get(newname)
          SingleMuonHisto.SetMarkerStyle(20)
          SingleMuonHisto.SetMarkerColor(1)
          SingleMuonHisto.SetLineColor(1)
          
          SelectedSignalSamples1Histo = SelectedSignalSamples1.Get(newname)
          SelectedSignalSamples1Histo.SetMarkerColor(6)
          SelectedSignalSamples1Histo.SetLineColor(6)
          
          SelectedSignalSamples2Histo = SelectedSignalSamples2.Get(newname)
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
          if (keyname2.find("Per")==-1) :
            axisXTitle = keyname2[keyname2.find("_")+1:]
            axisYTitle = "Tracks/bin"
          else :
            axisXTitle = keyname2[keyname2.find("_")+1:keyname2.find("Per")]
            axisYTitle = keyname2[keyname2.find("Per")+3:]
          
#              print("------------------------------------------------------------")
          if (obj.ClassName() == "TH1F") : # and "BS_" in keyname2):

            # array to contain a specific (keyname2) histogram for all samples
            histoArray = []
            nEventsPostTrigArray = []
            for index,fileIn in enumerate(bckArray):
              histo = fileIn.Get(newname)
              histoArray.append(histo)
              stackedSummedBackground.Add(histo)
              histo.SetLineColor(index+2)
              histo.SetFillColor(index+2)
              histo.SetMarkerColor(index+2)
              if ((index==0)) :
                legend.AddEntry(histo,"QCD p_{T}={50,inf}","LP")
              elif ((index==1)) :
                legend.AddEntry(histo,"WJets","LP")
              elif ((index==2)) :
                legend.AddEntry(histo,"TTBar","LP")
#              elif ((obj.ClassName() == "TH2F" or obj.ClassName() == "TH2D") and not (keyname2 == "genrecopT" or "BS_" in keyname2)):
#                histo2DArray = []
#                nEventsPostTrig2DArray = []
#                for fileIn in bckArray:
#                  if (fileIn.Get(newname).GetEntries() == 0 ) :
#                    continue
#                  fileIn.Get(newname).ProjectionY(newname,bin,bin,"e").Draw()
#                for index in range(0, len(histo2DArray)):
#                  if (index < 9) :
#                    stackedSummedQCD.Add(histo2DArray[index])
#                  elif (index>=9 and index<10) :
#                    stackedSummedW.Add(histo2DArray[index])
#                  elif (index>=10 and index<13) :
#                    stackedSummedTT.Add(histo2DArray[index])
            
#              elif (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
#                histo3DArray = []
#                nEventsPostTrig3DArray = []
#                for fileIn in bckArray:
#                  if (fileIn.Get(newname).GetEntries() == 0 ) :
#                      continue
#                  fileIn.Get(newname).GetXaxis().SetRange(bin,bin)
#                  fileIn.Get(newname).Project3DProfile("ZY").Draw()
#                  histo3DArray.append(fileIn.Get(newname).Project3DProfile("ZY"))
#                  nEvetsPostTrig3D = fileIn.Get("analyzer/BaseName/TotalTE").Integral()
#                  nEventsPostTrig3DArray.append(nEvetsPostTrig3D)
#                for index in range(0, len(histo3DArray)):
          else :
            continue
          #          convert stacks to (summed) histos
          stackedSummedBackground.Draw()
          stackedSummedBackgroundTmp = stackedSummedBackground.GetStack().Last()
          for index,fileIn in enumerate(bckArray):
            normHisto = fileIn.Get(newname)
            max = stackedSummedBackgroundTmp.GetMaximum()
            if (max==0) : continue
            normHisto.Scale(1/max)
            stackedNormSummedBackground.Add(normHisto)

          stackedNormSummedBackground.Draw()

          if (keyname2== "CutFlow") :
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
            stackedNormSummedBackground.GetXaxis().SetBinLabel(1,"Trigger")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(2,"Eta")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(3,"pT")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(4,"NumPixHits")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(5,"ValidFract")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(6,"NumDeDx")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(7,"ProbXY")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(8,"HighPurity")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(9,"Chi2oDOF")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(10,"EoP")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(11,"dz")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(12,"dxy")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(13,"pTerrOverpT")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(14,"N/A")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(15,"MiniIso")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(16,"PF ID")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(17,"Ih")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(18,"ProbQ")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(19,"MuStat")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(20,"PhiTOF")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(21,"EtaTOF")
          elif (keyname2== "CutFlowProbQFirst") :
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
            stackedNormSummedBackground.GetXaxis().SetBinLabel(1,"Trigger")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(2,"Eta")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(3,"pT")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(4,"ProbQ")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(5,"NumPixHits")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(6,"ValidFract")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(7,"NumDeDx")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(8,"ProbXY")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(9,"HighPurity")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(10,"Chi2oDOF")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(11,"EoP")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(12,"dz")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(13,"dxy")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(14,"pTerrOverpT")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(15,"N/A")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(16,"MiniIso")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(17,"PF ID")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(18,"Ih")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(19,"MuStat")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(20,"PhiTOF")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(21,"EtaTOF")
          elif ("pfType" in keyname2) :
            SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
            SelectedSignalSamples1Histo.Scale(1/SelectedSignalSamples1Histo.GetMaximum())
            SelectedSignalSamples2Histo.Scale(1/SelectedSignalSamples2Histo.GetMaximum())
            stackedNormSummedBackground.GetXaxis().SetBinLabel(1,"AllTracks")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(2,"PFtracks")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(3,"isElectron")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(4,"isMuon")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(5,"isPhoton")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(6,"isChHadron")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(7,"isNeutHadron")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(8,"isUndefined")
            stackedNormSummedBackground.GetXaxis().SetBinLabel(9,"else")
#          else :
#
#
#
          legend.AddEntry(SingleMuonHisto,"SingleMuon 2018C","LP")
          legend.AddEntry(SelectedSignalSamples1Histo, "HSCPgluino M=1800 GeV", "LP")
          legend.AddEntry(SelectedSignalSamples2Histo, "HSCPgluino M=2400 GeV", "LP")
#          legend.AddEntry(METHisto,"MET-EraC","LP")
          
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
          max = numpy.maximum(stackedSummedBackground.GetMaximum(),SingleMuonHisto.GetMaximum())

          cstackedSummedBackgroundString = 'cstackedSummedBackground'+str(j)
          cstackedSummedBackground = ROOT.TCanvas(cstackedSummedBackgroundString, cstackedSummedBackgroundString, 800,800)
          if ("CutFlow" in keyname2 or "pfType" in keyname2):
            stackedNormSummedBackground.Draw("HISTO")
            stackedNormSummedBackground.SetTitle("")
            stackedNormSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedNormSummedBackground.GetXaxis().SetTitle("")
            stackedNormSummedBackground.GetYaxis().SetRangeUser(0,1.5)
            stackedNormSummedBackground.SetMaximum(1.3)
          else :
            stackedSummedBackground.Draw("HISTO")
            stackedSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedSummedBackground.GetXaxis().SetTitleOffset(1)
            stackedSummedBackground.GetXaxis().SetTitle(axisXTitle)
            stackedSummedBackground.GetYaxis().SetTitle(axisYTitle)
            stackedSummedBackground.SetMaximum(max*1.4)
            stackedSummedBackground.SetMinimum(0.0)
          SingleMuonHisto.Draw("SAMEP")
          SelectedSignalSamples1Histo.Draw("SAME")
          SelectedSignalSamples2Histo.Draw("SAME")
          if ("Mass" in keyname2) :
            stackedSummedBackground.GetXaxis().SetRangeUser(0,1800)
#                stackedSummedBackground.GetXaxis().SetTitle("Mass [GeV]")
#                stackedSummedBackground.GetYaxis().SetTitle("Tracks/bin")
#                stackedSummedBackground.GetYaxis().SetTitleSize(0.05)
#                stackedSummedBackground.GetYaxis().SetTitleOffset(1)

 
          
          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")

          cstackedSummedBackground.SaveAs("StackedComparrison_wProbQ_CodeV"+codeVersion+"/"+keyname2+".png")
          
#---------------------------------------------------------------------------------------------
          # now let's do everything again but on a log Y scale
          cstackedSummedBackgroundLogString = 'cstackedSummedBackgroundLog'+str(j)
          cstackedSummedBackground = ROOT.TCanvas(cstackedSummedBackgroundLogString, cstackedSummedBackgroundLogString, 800,800)
          cstackedSummedBackground.SetLogy()
          
          if ("PrePreS" in keyname2 or "PostPreS" in keyname2 or "N1" in keyname2) :
            stackedSummedBackground.SetMinimum(0.0001)
          else:
            stackedSummedBackground.SetMinimum(0.000000000000000001)
            
          if ("CutFlow" in keyname2 or "pfType" in keyname2):
            stackedNormSummedBackground.Draw("HISTO")
            stackedNormSummedBackground.SetTitle("")
            stackedNormSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedNormSummedBackground.GetXaxis().SetTitle("")
            stackedNormSummedBackground.GetYaxis().SetRangeUser(0,1.5)
            stackedNormSummedBackground.SetMaximum(12)
          else :
            stackedSummedBackground.Draw("HISTO")
            stackedSummedBackground.GetXaxis().SetTitleSize(0.05)
            stackedSummedBackground.GetXaxis().SetTitleOffset(1)
            stackedSummedBackground.SetMaximum(max*1000)
          SingleMuonHisto.Draw("SAMEP")
          SelectedSignalSamples1Histo.Draw("SAME")
          SelectedSignalSamples2Histo.Draw("SAME")
          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")

          cstackedSummedBackground.SaveAs("StackedComparrison_wProbQ_CodeV"+codeVersion+"/"+keyname2+"_log.png")
