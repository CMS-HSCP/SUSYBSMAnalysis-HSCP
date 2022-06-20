import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog binNumber")
(opt,args) = parser.parse_args()

BinNumber = sys.argv[1]
bin = int(BinNumber)

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

BackgroundSamples = [
#"crab_Analysis_2018_QCD_Pt-15To20_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
#"crab_Analysis_2018_QCD_Pt-20To30_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
#"crab_Analysis_2018_QCD_Pt-30To50_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_WJetsToLNu_0J_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_TTToHadronic_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_TTToSemiLeptonic_wProbQ_CodeV19p1_v1.root",
"crab_Analysis_2018_TTTo2L2Nu_wProbQ_CodeV19p1_v1.root",
]

QCDFlat = ROOT.TFile.Open("crab_Analysis_2018_QCD_Pt-20_MuEnrichedPt15_wProbQ_CodeV18p9_v1.root")
SingleMuon = ROOT.TFile.Open("crab_Analysis_SingleMuon_Run2017_wProbQ_CodeV19p2_v1.root")
MET = ROOT.TFile.Open("crab_Analysis_MET_Run2017_wProbQ_CodeV19p2_v1.root")

intLumi = 10.0 # 137.0

crossSectionArray = [
# QCD_Pt-20_MuEnrichedPt15_TuneCP5_13TeV-pythia8, 239000.0 +-	755.8
# 2797000.0, #+-8800.0, 2018_QCD_Pt-15To20_MuEnrichedPt5
# 2518000.0, #+-7940.0, 2018_QCD_Pt-20To30_MuEnrichedPt5
# 1361000.0, #+-4263.0, 2018_QCD_Pt-30To50_MuEnrichedPt5
 377800.0, #+-1184.0, 2018_QCD_Pt-50To80_MuEnrichedPt5
 88620.0, #+-275.1, 2018_QCD_Pt-80To120_MuEnrichedPt5
 21070.0, #+-65.28, 2018_QCD_Pt-120To170_MuEnrichedPt5
 7019.0, #+-21.61, 2018_QCD_Pt-170To300_MuEnrichedPt5
 622.4, #+-1.891, 2018_QCD_Pt-300To470_MuEnrichedPt5
 58.86, #+-0.1776, 2018_QCD_Pt-470To600_MuEnrichedPt5
 18.22, #+-0.05471, 2018_QCD_Pt-600To800_MuEnrichedPt5
 3.25, #+-0.0148, 2018_QCD_Pt-600To800_MuEnrichedPt5
 1.613, #+-, 2018_QCD_Pt-1000_MuEnrichedPt5
 53330.0, #+-	90.89, WJetsToLNu
 687.1, # or 377.96, # or , TTtoHadronic
 687.1, # or 65.34, TTtoSemiLeptonic
 88.29, # or 687.1, TTto2L2N
]

crossSectionFlat = 239000.0	# +-755.8

fileInArray = []
for sample in BackgroundSamples:
  fileInArray.append(ROOT.TFile.Open(sample))
  

nEventsPostTrigArrayFlat = QCDFlat.Get("analyzer/BaseName/TotalTE").Integral()


dirs = []
for i in range(0, fileInArray[0].GetListOfKeys().GetEntries()):
  dirname = fileInArray[0].GetListOfKeys().At(i).GetName()
  curr_dir = fileInArray[0].GetDirectory(dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = fileInArray[0].GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          newname = dirname + "/" + keyname+ "/" + keyname2
          obj = fileInArray[0].Get(newname)
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if obj.InheritsFrom("TObject"):
              if not os.path.exists(os.path.dirname("StackedComparrison_"+(BackgroundSamples[0])[BackgroundSamples[0].find("Code"):BackgroundSamples[0].find("Code")+9]+"/")): os.makedirs(os.path.dirname("StackedComparrison_"+(BackgroundSamples[0])[BackgroundSamples[0].find("Code"):BackgroundSamples[0].find("Code")+9]+"/"))
              if (obj.GetEntries() == 0 ) :
                continue
              if ("Total" in keyname2) :
                continue
                
              stackPlotsTemp = ROOT.THStack("stackPlotsTemp","")
              stackPlots = ROOT.THStack("stackPlots","")
              stackedSummedQCD = ROOT.THStack("stackedSummedQCD","")
              stackedSummedTT = ROOT.THStack("stackedSummedTT","")
              stackedSummedW = ROOT.THStack("stackedSummedW","")
              histoSummedQCD = ROOT.TH1F() #"histoSummedQCD")
              histoSummedTT = ROOT.TH1F() #"histoSummedTT")
              histoSummedW = ROOT.TH1F() #"histoSummedW")
              
              legMass =  ROOT.TLegend(.55,.70,.80,.89,"","brNDC")
              legMass.SetTextFont(42)
              legMass.SetTextSize(0.02)
              legMass.SetBorderSize(1);
              legMass.SetBorderSize(0);
              legMass.SetLineColor(1);
              legMass.SetLineStyle(1);
              legMass.SetLineWidth(1);
              legMass.SetFillColor(0);
              legMass.SetFillStyle(1001);
              
#              print("------------------------------------------------------------")
              if (obj.ClassName() == "TH1F") : # and "BS_" in keyname2):
                # array to contain a specific (keyname2) histogram for all samples
                histoArray = []
                nEventsPostTrigArray = []
                for fileIn in fileInArray:
                  histoArray.append(fileIn.Get(newname))
                  nEvetsPostTrig = fileIn.Get("analyzer/BaseName/TotalTE").Integral()
                  nEventsPostTrigArray.append(nEvetsPostTrig)
                for index in range(0, len(histoArray)):
                  histoArray[index].Scale(intLumi*crossSectionArray[index]/nEventsPostTrigArray[index])
                  if (index < 9) :
                    stackedSummedQCD.Add(histoArray[index])
                  elif (index>=9 and index<10) :
                    stackedSummedW.Add(histoArray[index])
                  elif (index>=10 and index<13) :
                    stackedSummedTT.Add(histoArray[index])
                    
              elif ((obj.ClassName() == "TH2F" or obj.ClassName() == "TH2D") and not (keyname2 == "genrecopT" or "BS_" in keyname2)):
                histo2DArray = []
                nEventsPostTrig2DArray = []
                for fileIn in fileInArray:
                  if (fileIn.Get(newname).GetEntries() == 0 ) :
                    continue
                  fileIn.Get(newname).ProjectionY(newname,bin,bin,"e").Draw()
                  histo2DArray.append(fileIn.Get(newname).ProjectionY(newname,bin,bin,"e"))
                  nEvetsPostTrig2D = fileIn.Get("analyzer/BaseName/TotalTE").Integral()
                  nEventsPostTrig2DArray.append(nEvetsPostTrig2D)
                for index in range(0, len(histo2DArray)):
                  histo2DArray[index].Scale(intLumi*crossSectionArray[index]/nEventsPostTrig2DArray[index])
                  if (index < 9) :
                    stackedSummedQCD.Add(histo2DArray[index])
                  elif (index>=9 and index<10) :
                    stackedSummedW.Add(histo2DArray[index])
                  elif (index>=10 and index<13) :
                    stackedSummedTT.Add(histo2DArray[index])
                
                
#              elif (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
#                histo3DArray = []
#                nEventsPostTrig3DArray = []
#                for fileIn in fileInArray:
#                  if (fileIn.Get(newname).GetEntries() == 0 ) :
#                      continue
#                  fileIn.Get(newname).GetXaxis().SetRange(bin,bin)
#                  fileIn.Get(newname).Project3DProfile("ZY").Draw()
#                  histo3DArray.append(fileIn.Get(newname).Project3DProfile("ZY"))
#                  nEvetsPostTrig3D = fileIn.Get("analyzer/BaseName/TotalTE").Integral()
#                  nEventsPostTrig3DArray.append(nEvetsPostTrig3D)
#                for index in range(0, len(histo3DArray)):
#                  histo3DArray[index].Scale(intLumi*crossSectionArray[index]/nEventsPostTrig3DArray[index])
#                  if (index < 9) :
#                    stackedSummedQCD.Add(histo3DArray[index])
#                  elif (index>=9 and index<10) :
#                    stackedSummedW.Add(histo3DArray[index])
#                  elif (index>=10 and index<13) :
#                    stackedSummedTT.Add(histo3DArray[index])
              else :
                continue

#               convert stacks to (summed) histos
              histoSummedQCD = stackedSummedQCD.GetStack().Last()
              histoSummedTT = stackedSummedTT.GetStack().Last()
              histoSummedW = stackedSummedW.GetStack().Last()
              
              stackPlotsTemp.Add(histoSummedTT)
              histoSummedTT.SetMarkerStyle(20)
              histoSummedTT.SetLineColor(3)
              histoSummedTT.SetFillColor(3)
              histoSummedTT.SetMarkerColor(3)
              
              stackPlotsTemp.Add(histoSummedW)
              histoSummedW.SetMarkerStyle(20)
              histoSummedW.SetLineColor(4)
              histoSummedW.SetFillColor(4)
              histoSummedW.SetMarkerColor(4)
              
              
              stackPlotsTemp.Add(histoSummedQCD)
              histoSummedQCD.SetMarkerStyle(20)
              histoSummedQCD.SetLineColor(2)
              histoSummedQCD.SetFillColor(2)
              histoSummedQCD.SetMarkerColor(2)
              
              QCDFlatHisto = QCDFlat.Get(newname)
              SingleMuonHisto = SingleMuon.Get(newname)
              METHisto = MET.Get(newname)
              
              QCDFlatHisto.Draw()
              SingleMuonHisto.Draw()
              METHisto.Draw()
              QCDFlatHisto.SetMarkerStyle(20)
              QCDFlatHisto.SetMarkerColor(1)
              QCDFlatHisto.SetLineColor(1)
              
              QCDFlatHisto.Scale(intLumi*crossSectionFlat/nEventsPostTrigArrayFlat)
              
              SingleMuonHisto.SetMarkerStyle(20)
              SingleMuonHisto.SetMarkerColor(4)
              SingleMuonHisto.SetLineColor(4)
              
              METHisto.SetMarkerStyle(20)
              METHisto.SetMarkerColor(5)
              METHisto.SetLineColor(5)
              
              
              if (keyname2== "CutFlow") :
                histoSummedTT.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedW.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedQCD.Scale(1/stackPlotsTemp.GetMaximum())
                QCDFlatHisto.Scale(1/QCDFlatHisto.GetMaximum())
                SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
                METHisto.Scale(1/METHisto.GetMaximum())
                histoSummedTT.GetXaxis().SetBinLabel(1,"Trigger")
                histoSummedTT.GetXaxis().SetBinLabel(2,"Eta")
                histoSummedTT.GetXaxis().SetBinLabel(3,"pT")
                histoSummedTT.GetXaxis().SetBinLabel(4,"NumHits")
                histoSummedTT.GetXaxis().SetBinLabel(5,"NumPixHits")
                histoSummedTT.GetXaxis().SetBinLabel(6,"ValidFract")
                histoSummedTT.GetXaxis().SetBinLabel(7,"NumDeDx")
                histoSummedTT.GetXaxis().SetBinLabel(8,"ProbXY")
                histoSummedTT.GetXaxis().SetBinLabel(9,"HighPurity")
                histoSummedTT.GetXaxis().SetBinLabel(10,"Chi2oDOF")
                histoSummedTT.GetXaxis().SetBinLabel(11,"EoP")
                histoSummedTT.GetXaxis().SetBinLabel(12,"dz")
                histoSummedTT.GetXaxis().SetBinLabel(13,"dxy")
                histoSummedTT.GetXaxis().SetBinLabel(14,"pTerrOverpT")
                histoSummedTT.GetXaxis().SetBinLabel(15,"TKIso")
                histoSummedTT.GetXaxis().SetBinLabel(16,"MiniIso")
                histoSummedTT.GetXaxis().SetBinLabel(17,"MassT")
                histoSummedTT.GetXaxis().SetBinLabel(18,"Ih")
                histoSummedTT.GetXaxis().SetBinLabel(19,"ProbQ")
                histoSummedTT.GetXaxis().SetBinLabel(20,"MuStat")
                histoSummedTT.GetXaxis().SetBinLabel(21,"PhiTOF")
                histoSummedTT.GetXaxis().SetBinLabel(22,"EtaTOF")
              elif (keyname2== "CutFlowProbQFirst") :
                histoSummedTT.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedW.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedQCD.Scale(1/stackPlotsTemp.GetMaximum())
                QCDFlatHisto.Scale(1/QCDFlatHisto.GetMaximum())
                SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
                METHisto.Scale(1/METHisto.GetMaximum())
                histoSummedTT.GetXaxis().SetBinLabel(1,"Trigger")
                histoSummedTT.GetXaxis().SetBinLabel(2,"Eta")
                histoSummedTT.GetXaxis().SetBinLabel(3,"pT")
                histoSummedTT.GetXaxis().SetBinLabel(4,"ProbQ")
                histoSummedTT.GetXaxis().SetBinLabel(5,"NumHits")
                histoSummedTT.GetXaxis().SetBinLabel(6,"NumPixHits")
                histoSummedTT.GetXaxis().SetBinLabel(7,"ValidFract")
                histoSummedTT.GetXaxis().SetBinLabel(8,"NumDeDx")
                histoSummedTT.GetXaxis().SetBinLabel(9,"ProbXY")
                histoSummedTT.GetXaxis().SetBinLabel(10,"HighPurity")
                histoSummedTT.GetXaxis().SetBinLabel(11,"Chi2oDOF")
                histoSummedTT.GetXaxis().SetBinLabel(12,"EoP")
                histoSummedTT.GetXaxis().SetBinLabel(13,"dz")
                histoSummedTT.GetXaxis().SetBinLabel(14,"dxy")
                histoSummedTT.GetXaxis().SetBinLabel(15,"pTerrOverpT")
                histoSummedTT.GetXaxis().SetBinLabel(16,"TKIso")
                histoSummedTT.GetXaxis().SetBinLabel(17,"MiniIso")
                histoSummedTT.GetXaxis().SetBinLabel(18,"MassT")
                histoSummedTT.GetXaxis().SetBinLabel(19,"Ih")
                histoSummedTT.GetXaxis().SetBinLabel(20,"MuStat")
                histoSummedTT.GetXaxis().SetBinLabel(21,"PhiTOF")
                histoSummedTT.GetXaxis().SetBinLabel(22,"EtaTOF")
              elif (keyname2== "pfType") :
                histoSummedTT.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedW.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedQCD.Scale(1/stackPlotsTemp.GetMaximum())
                QCDFlatHisto.Scale(1/QCDFlatHisto.GetMaximum())
                SingleMuonHisto.Scale(1/SingleMuonHisto.GetMaximum())
                METHisto.Scale(1/METHisto.GetMaximum())
                histoSummedTT.GetXaxis().SetBinLabel(1,"AllTracks")
                histoSummedTT.GetXaxis().SetBinLabel(2,"PFtracks")
                histoSummedTT.GetXaxis().SetBinLabel(3,"isElectron")
                histoSummedTT.GetXaxis().SetBinLabel(4,"isMuon")
                histoSummedTT.GetXaxis().SetBinLabel(5,"isPhoton")
                histoSummedTT.GetXaxis().SetBinLabel(6,"isChHadron")
                histoSummedTT.GetXaxis().SetBinLabel(7,"isNeutHadron")
                histoSummedTT.GetXaxis().SetBinLabel(8,"isUndefined")
                histoSummedTT.GetXaxis().SetBinLabel(9,"else")
              else :
                histoSummedTT.GetYaxis().SetTitle("Tracks/bin")
                histoSummedTT.GetXaxis().SetTitle(keyname2)
                
              stackPlots.Add(histoSummedTT)
              stackPlots.Add(histoSummedW)
              stackPlots.Add(histoSummedQCD)
            
              legMass.AddEntry(histoSummedQCD,"#mu-enriched QCD (full pt)","LP")
              legMass.AddEntry(histoSummedW,"W process","LP")
              legMass.AddEntry(histoSummedTT,"ttbar process","LP")
              legMass.AddEntry(QCDFlatHisto,"2018_QCD_Pt-20ToInf_MuPt15","LP")
              legMass.AddEntry(SingleMuonHisto,"SingleMuon-EraC","LP")
              legMass.AddEntry(METHisto,"MET-EraC","LP")
              
              tex2 = ROOT.TLatex(0.13,0.94,"CMS");
              #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
              tex2.SetNDC();
              tex2.SetTextFont(61);
              tex2.SetTextSize(0.0675);
              tex2.SetLineWidth(2);

              tex3 = ROOT.TLatex(0.27,0.94,"Simulation"); # for square plots
              #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
              #tex3 = ROOT.TLatex(0.28,0.94,"Internal");
              tex3.SetNDC();
              tex3.SetTextFont(52);
              tex3.SetTextSize(0.0485);
              tex3.SetLineWidth(2);

              cstackPlotsString = 'cstackPlots'+str(j)
              cstackPlots = ROOT.TCanvas(cstackPlotsString, cstackPlotsString, 800,800)
              
              stackPlots.Draw("HISTO")
              QCDFlatHisto.Draw("SAME")
              SingleMuonHisto.Draw("SAME")
              METHisto.Draw("SAME")
#                stackPlots.SetTitle("")
              stackPlots.GetXaxis().SetTitleSize(0.05)
              stackPlots.GetXaxis().SetTitleOffset(1)
              if ("Mass" in keyname2) :
                stackPlots.GetXaxis().SetRangeUser(0,1800)
#                stackPlots.GetXaxis().SetTitle("Mass [GeV]")
#                stackPlots.GetYaxis().SetTitle("Tracks/bin")
#                stackPlots.GetYaxis().SetTitleSize(0.05)
#                stackPlots.GetYaxis().SetTitleOffset(1)

              max1 = numpy.maximum(stackPlots.GetMaximum(),QCDFlatHisto.GetMaximum())
              max2 = numpy.maximum(SingleMuonHisto.GetMaximum(),METHisto.GetMaximum())
              max = numpy.maximum(max1,max2)
              stackPlots.SetMaximum(max*1.4)
              stackPlots.SetMinimum(0.0)
              legMass.Draw("SAME")
              tex2.Draw("SAME")
              tex3.Draw("SAME")
              
#              QCDFlatHisto = QCDFlat.Get(newname)
#              QCDFlatHisto.SetMarkerStyle(20)
#              SingleMuonHisto = SingleMuon.Get(newname)
#              SingleMuonHisto.SetMarkerStyle(20)
#              SingleMuonHisto.SetMarkerColor(4)
#              SingleMuonHisto.SetLineColor(4)
#              METHisto = MET.Get(newname)
#              METHisto.SetMarkerStyle(20)
#              METHisto.SetMarkerColor(5)
#              METHisto.SetLineColor(5)

              cstackPlots.SaveAs("StackedComparrison_"+(BackgroundSamples[0])[BackgroundSamples[0].find("Code"):BackgroundSamples[0].find("Code")+9]+"/"+keyname2+".png")
              
              if ("Mass" in keyname2) :
                cstackPlotsZoomedString = 'cstackPlotsZoomed'+str(j)
                cstackPlotsZoomed = ROOT.TCanvas(cstackPlotsZoomedString, cstackPlotsZoomedString, 800,800)
                stackPlotsZoomed = stackPlots.Clone()
                stackPlotsZoomed.GetXaxis().SetRangeUser(500,1800)
                stackPlotsZoomed.GetYaxis().SetRangeUser(0,0.0000000000000025)
                stackPlotsZoomed.SetMaximum(0.0000000000000025)
#                stackPlots.SetMinimum(0.0)
                stackPlotsZoomed.Draw("HISTO")
                QCDFlatHisto.Draw("SAME")
                SingleMuonHisto.Draw("SAME")
                METHisto.Draw("SAME")
                legMass.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                cstackPlotsZoomed.SaveAs("StackedComparrison_"+(BackgroundSamples[0])[BackgroundSamples[0].find("Code"):BackgroundSamples[0].find("Code")+9]+"/"+keyname2+"_zoomed.png")
                
              
              cstackPlotsLogString = 'cstackPlotsLog'+str(j)
              cstackPlots = ROOT.TCanvas(cstackPlotsLogString, cstackPlotsLogString, 800,800)
              cstackPlots.SetLogy()
              stackPlots.Draw()

#                stackPlots.GetXaxis().SetTitle("Mass [GeV]")
#                stackPlots.GetYaxis().SetTitle("Tracks/bin")
#                stackPlots.GetYaxis().SetTitleSize(0.05)
#                stackPlots.GetYaxis().SetTitleOffset(1)

              max = numpy.maximum(stackPlots.GetMaximum(),stackPlots.GetMaximum())
              stackPlots.SetMaximum(max*100000)
              if ("BS" in keyname2 or "PostPreS" in keyname2) :
                stackPlots.SetMinimum(0.0001)
              else:
                stackPlots.SetMinimum(0.000000000000000001)
                
              stackPlots.Draw("HISTO")
              QCDFlatHisto.Draw("SAME")
              SingleMuonHisto.Draw("SAME")
              METHisto.Draw("SAME")
#                stackPlots.SetTitle("")
              stackPlots.GetXaxis().SetTitleSize(0.05)
              stackPlots.GetXaxis().SetTitleOffset(1)
              stackPlots.GetYaxis().SetRangeUser(0,1)
              
              legMass.Draw("SAME")
              tex2.Draw("SAME")
              tex3.Draw("SAME")

              cstackPlots.SaveAs("StackedComparrison_"+(BackgroundSamples[0])[BackgroundSamples[0].find("Code"):BackgroundSamples[0].find("Code")+9]+"/"+keyname2+"_log.png")
