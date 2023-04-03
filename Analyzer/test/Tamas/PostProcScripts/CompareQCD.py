import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog")
(opt,args) = parser.parse_args()

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

QCDBinnedSamples = [
#"crab_Analysis_2018_QCD_Pt-15To20_MuEnrichedPt5_wProbQ_CodeV18p2_v1.root",
#"crab_Analysis_2018_QCD_Pt-20To30_MuEnrichedPt5_wProbQ_CodeV18p2_v1.root",
#"crab_Analysis_2018_QCD_Pt-30To50_MuEnrichedPt5_wProbQ_CodeV18p2_v1.root",
"crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
"crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_wProbQ_CodeV18p9_v1.root",
]

QCDFlat = ROOT.TFile.Open("crab_Analysis_2018_QCD_Pt-20_MuEnrichedPt15_wProbQ_CodeV18p9_v1.root")

intLumi = 50.0 #137.0

SingleMuon = ROOT.TFile.Open("crab_Analysis_SingleMuon_Run2017_wProbQ_CodeV19p2_v1.root")

MET = ROOT.TFile.Open("crab_Analysis_MET_Run2017_wProbQ_CodeV19p2_v1.root")

crossSectionArray = [
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
]

fileInArray = []
for sample in QCDBinnedSamples:
  fileInArray.append(ROOT.TFile.Open(sample))
  

crossSectionFlat = 239000.0	# +-755.8
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
          if ("Total" in keyname2) :
            continue
          if obj.InheritsFrom("TObject"):
              if not os.path.exists(os.path.dirname("CompareQCDv4/")): os.makedirs(os.path.dirname("CompareQCDv4/"))
              if (obj.GetEntries() == 0 ) :
                continue
              if (obj.ClassName() == "TH1F"): # and "BS_" in keyname2):
                legMass =  ROOT.TLegend(.55,.60,.80,.89,"","brNDC")
                legMass.SetTextFont(42)
                legMass.SetTextSize(0.02)
                legMass.SetBorderSize(1);
                legMass.SetBorderSize(0);
                legMass.SetLineColor(1);
                legMass.SetLineStyle(1);
                legMass.SetLineWidth(1);
                legMass.SetFillColor(0);
                legMass.SetFillStyle(1001);
                
                stackPlots = ROOT.THStack("stackPlots","");
                histoArray = []
                nEventsPostTrigArray = []
                for fileIn in fileInArray:
                  histoArray.append(fileIn.Get(newname))
                  nEvetsPostTrig = fileIn.Get("analyzer/BaseName/TotalTE").Integral()
                  nEventsPostTrigArray.append(nEvetsPostTrig)
                for index in range(0, len(histoArray)):
                  
                  histoArray[index].Scale(intLumi*crossSectionArray[index]/nEventsPostTrigArray[index])
#                  finalHisto.Add(histoArray[index],intLumi*crossSectionArray[index]/nEventsPostTrigArray[index])
                  stackPlots.Add(histoArray[index])
                  legMass.AddEntry(histoArray[index],QCDBinnedSamples[index][14:QCDBinnedSamples[index].find("MuEnrichedPt5_")-1],"LP")
                  indexNew = -1
                  if (index>-1):
                    indexNew = index+2
                  if (indexNew==10) :
                    indexNew = 30
                  elif (indexNew==11) :
                    indexNew = 40
                  elif (indexNew==12) :
                    indexNew = 41
                  elif (indexNew==13) :
                    indexNew = 42
                  elif (indexNew==14) :
                    indexNew = 46
                  histoArray[index].SetLineColor(indexNew)
                  histoArray[index].SetFillColor(indexNew)

                QCDFlatHisto = QCDFlat.Get(newname)
                QCDFlatHisto.SetMarkerStyle(20)
                SingleMuonHisto = SingleMuon.Get(newname)
                SingleMuonHisto.SetMarkerStyle(20)
                SingleMuonHisto.SetMarkerColor(4)
                SingleMuonHisto.SetLineColor(4)
                METHisto = MET.Get(newname)
                METHisto.SetMarkerStyle(20)
                METHisto.SetMarkerColor(5)
                METHisto.SetLineColor(5)
                
                if (keyname2 == "BS_MPt") :
                  stackPlots.Draw()
                  stackPlots.GetXaxis().SetRangeUser(0,1000)
                  QCDFlatHisto.GetXaxis().SetRangeUser(0,1000)
                
                stackPlotsTemp = stackPlots.Clone()
                if (keyname2== "CutFlow") :
                histoSummedTT.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedW.Scale(1/stackPlotsTemp.GetMaximum())
                histoSummedQCD.Scale(1/stackPlotsTemp.GetMaximum())
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
                histoSummedTT.GetXaxis().SetBinLabel(1,"AllTracks")
                histoSummedTT.GetXaxis().SetBinLabel(2,"PFtracks")
                histoSummedTT.GetXaxis().SetBinLabel(3,"isElectron")
                histoSummedTT.GetXaxis().SetBinLabel(4,"isMuon")
                histoSummedTT.GetXaxis().SetBinLabel(5,"isPhoton")
                histoSummedTT.GetXaxis().SetBinLabel(6,"isChHadron")
                histoSummedTT.GetXaxis().SetBinLabel(7,"isNeutHadron")
                histoSummedTT.GetXaxis().SetBinLabel(8,"isUndefined")
                histoSummedTT.GetXaxis().SetBinLabel(9,"else")
                
                QCDFlatHisto.Scale(intLumi*crossSectionFlat/nEventsPostTrigArrayFlat)
                
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
                stackPlots.Draw("HIST")
                QCDFlatHisto.Draw("SAME")
                SingleMuonHisto.Draw("SAME")
                METHisto.Draw("SAME")
#                stackPlots.SetTitle("")
                stackPlots.GetXaxis().SetTitleSize(0.05)
                stackPlots.GetXaxis().SetTitleOffset(1)
                stackPlots.GetYaxis().SetRangeUser(0,1)
#                stackPlots.GetXaxis().SetTitle("Mass [GeV]")
#                stackPlots.GetYaxis().SetTitle("Tracks/bin")
#                stackPlots.GetYaxis().SetTitleSize(0.05)
#                stackPlots.GetYaxis().SetTitleOffset(1)

                max1 = numpy.maximum(stackPlots.GetMaximum(),QCDFlatHisto.GetMaximum())
                max2 = numpy.maximum(SingleMuonHisto.GetMaximum(),METHisto.GetMaximum())
                max = numpy.maximum(max1,max2)
                stackPlots.SetMaximum(max*1.5)
                stackPlots.SetMinimum(0.001)
                legMass.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")

                cstackPlots.SaveAs("CompareQCDv4/"+keyname2+".png")
