import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python3 %prog sample.txt")
(opt,args) = parser.parse_args()

sampleInFile = sys.argv[1]

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

SamplesArray = []

bin = 3

with open(sampleInFile, "r") as a_file:
  for line in a_file:
    stripped_line = line.strip()
    SamplesArray.append(stripped_line)

fileInArray = []
for sample in SamplesArray:
  fileInArray.append(ROOT.TFile.Open(sample))
  
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
          obj.SetMarkerStyle(20)
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if ("Total" in keyname2 or "IntLumi" in keyname2 or "XSection" in keyname2) :
            continue
          if obj.InheritsFrom("TObject"):
              if not os.path.exists(os.path.dirname("Compare"+sampleInFile[:-4]+"/")): os.makedirs(os.path.dirname("Compare"+sampleInFile[:-4]+"/"))
              if (obj.GetEntries() == 0 ) :
                continue
              if (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
                # this maybe should go from bin to bin+1 ?
                obj.SetMarkerStyle(20)
                if ("" in keyname2):
                  obj.GetXaxis().SetRange(bin,bin)
                  obj.Project3D("ZY").Draw("COLZ")
                else :
                  obj.GetXaxis().SetRange(bin,bin)
                  obj.Project3D("ZY").Draw("COLZ")
              if (obj.ClassName() == "TH1F" or obj.ClassName() == "TProfile"): # and "BS_" in keyname2):
                canvasString = 'canvas'+str(j)
                canvas = ROOT.TCanvas(canvasString, canvasString, 800,800)
                legend =  ROOT.TLegend(.55,.60,.80,.89,"","brNDC")
                legend.SetTextFont(42)
                legend.SetTextSize(0.017)
                legend.SetBorderSize(1);
                legend.SetBorderSize(0);
                legend.SetLineColor(1);
                legend.SetLineStyle(1);
                legend.SetLineWidth(1);
                legend.SetFillColor(0);
                legend.SetFillStyle(1001);
                
                histoArray = []
                for fileIn in fileInArray:
                  histoArray.append(fileIn.Get(newname))
                for index in range(0, len(histoArray)):
                  histoArray[index].SetStats(0)
                  histoArray[index].SetMarkerStyle(20)
                  
                  legend.AddEntry(histoArray[index],SamplesArray[index][19:SamplesArray[index].find("Code")+9],"LP")
                  indexNew = -1
                  if (index>-1):
                    indexNew = index+2
                  if (indexNew==10) :
                    indexNew = 40
                  elif (indexNew==11) :
                    indexNew = 46
                  elif (indexNew==12) :
                    indexNew = 41
                  elif (indexNew==13) :
                    indexNew = 30
                  elif (indexNew==14) :
                    indexNew = 42
                  histoArray[index].SetLineColor(indexNew)
                  histoArray[index].SetMarkerColor(indexNew)
#                  histoArray[index].SetFillColor(indexNew)
                  histoArray[index].SetTitle("")
                  max = 0.0
#                  print(keyname2)
                  for index2 in range(0, len(histoArray)):
                    max = numpy.maximum(max,histoArray[index2].GetMaximum())
                  histoArray[index].GetYaxis().SetTitle("Tracks/bin")
                  histoArray[index].GetXaxis().SetTitle(keyname2)
                  if (keyname2 == "pfType") :
                    histoArray[index].Scale(1/histoArray[index].GetMaximum())
                    histoArray[index].GetYaxis().SetRangeUser(0.,1.4)
                    histoArray[index].GetXaxis().SetBinLabel(1,"AllTracks")
                    histoArray[index].GetXaxis().SetBinLabel(2,"PFtracks")
                    histoArray[index].GetXaxis().SetBinLabel(3,"isElectron")
                    histoArray[index].GetXaxis().SetBinLabel(4,"isMuon")
                    histoArray[index].GetXaxis().SetBinLabel(5,"isPhoton")
                    histoArray[index].GetXaxis().SetBinLabel(6,"isChHadron")
                    histoArray[index].GetXaxis().SetBinLabel(7,"isNeutHadron")
                    histoArray[index].GetXaxis().SetBinLabel(8,"isUndefined")
                    histoArray[index].GetXaxis().SetBinLabel(9,"else")
                    histoArray[index].GetXaxis().SetTitle("")
                    histoArray[index].GetYaxis().SetTitle("")
                    histoArray[index].SetMaximum(1.4)
                  elif (keyname2== "CutFlow") :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    histoArray[index].Scale(1/histoArray[index].GetMaximum())
                    histoArray[index].GetYaxis().SetRangeUser(0.,1.4)
                    histoArray[index].GetXaxis().SetBinLabel(1,"Trigger")
                    histoArray[index].GetXaxis().SetBinLabel(2,"Eta")
                    histoArray[index].GetXaxis().SetBinLabel(3,"pT")
                    histoArray[index].GetXaxis().SetBinLabel(4,"NumHits")
                    histoArray[index].GetXaxis().SetBinLabel(5,"NumPixHits")
                    histoArray[index].GetXaxis().SetBinLabel(6,"ValidFract")
                    histoArray[index].GetXaxis().SetBinLabel(7,"NumDeDx")
                    histoArray[index].GetXaxis().SetBinLabel(8,"ProbXY")
                    histoArray[index].GetXaxis().SetBinLabel(9,"HighPurity")
                    histoArray[index].GetXaxis().SetBinLabel(10,"Chi2oDOF")
                    histoArray[index].GetXaxis().SetBinLabel(11,"EoP")
                    histoArray[index].GetXaxis().SetBinLabel(12,"dz")
                    histoArray[index].GetXaxis().SetBinLabel(13,"dxy")
                    histoArray[index].GetXaxis().SetBinLabel(14,"pTerrOverpT")
                    histoArray[index].GetXaxis().SetBinLabel(15,"SVfromNI")
                    histoArray[index].GetXaxis().SetBinLabel(16,"MiniIso")
                    histoArray[index].GetXaxis().SetBinLabel(17,"PFid")
                    histoArray[index].GetXaxis().SetBinLabel(18,"Ih")
                    histoArray[index].GetXaxis().SetBinLabel(19,"ProbQ")
                    histoArray[index].GetXaxis().SetBinLabel(20,"MuStat")
                    histoArray[index].GetXaxis().SetBinLabel(21,"PhiTOF")
                    histoArray[index].GetXaxis().SetBinLabel(22,"EtaTOF")
                    histoArray[index].GetXaxis().SetTitle("")
                    histoArray[index].GetYaxis().SetTitle("")
                    histoArray[index].SetMaximum(1.4)
                  elif (keyname2== "CutFlowProbQFirst") :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    histoArray[index].Scale(1/histoArray[index].GetMaximum())
                    histoArray[index].GetYaxis().SetRangeUser(0.,1.4)
                    histoArray[index].GetXaxis().SetBinLabel(1,"Trigger")
                    histoArray[index].GetXaxis().SetBinLabel(2,"Eta")
                    histoArray[index].GetXaxis().SetBinLabel(3,"pT")
                    histoArray[index].GetXaxis().SetBinLabel(4,"ProbQ")
                    histoArray[index].GetXaxis().SetBinLabel(5,"NumHits")
                    histoArray[index].GetXaxis().SetBinLabel(6,"NumPixHits")
                    histoArray[index].GetXaxis().SetBinLabel(7,"ValidFract")
                    histoArray[index].GetXaxis().SetBinLabel(8,"NumDeDx")
                    histoArray[index].GetXaxis().SetBinLabel(9,"ProbXY")
                    histoArray[index].GetXaxis().SetBinLabel(10,"HighPurity")
                    histoArray[index].GetXaxis().SetBinLabel(11,"Chi2oDOF")
                    histoArray[index].GetXaxis().SetBinLabel(12,"EoP")
                    histoArray[index].GetXaxis().SetBinLabel(13,"dz")
                    histoArray[index].GetXaxis().SetBinLabel(14,"dxy")
                    histoArray[index].GetXaxis().SetBinLabel(15,"pTerrOverpT")
                    histoArray[index].GetXaxis().SetBinLabel(16,"SVfromNI")
                    histoArray[index].GetXaxis().SetBinLabel(17,"MiniIso")
                    histoArray[index].GetXaxis().SetBinLabel(18,"PFid")
                    histoArray[index].GetXaxis().SetBinLabel(19,"Ih")
                    histoArray[index].GetXaxis().SetBinLabel(20,"MuStat")
                    histoArray[index].GetXaxis().SetBinLabel(21,"PhiTOF")
                    histoArray[index].GetXaxis().SetBinLabel(22,"EtaTOF")
                    histoArray[index].GetXaxis().SetTitle("")
                    histoArray[index].GetYaxis().SetTitle("")
                    histoArray[index].SetMaximum(1.4)
                  elif ("pfType" in keyname2) :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    histoArray[index].Scale(1/histoArray[index].GetMaximum())
                    histoArray[index].GetXaxis().SetBinLabel(1,"AllTracks")
                    histoArray[index].GetXaxis().SetBinLabel(2,"PFtracks")
                    histoArray[index].GetXaxis().SetBinLabel(3,"isElectron")
                    histoArray[index].GetXaxis().SetBinLabel(4,"isMuon")
                    histoArray[index].GetXaxis().SetBinLabel(5,"isPhoton")
                    histoArray[index].GetXaxis().SetBinLabel(6,"isChHadron")
                    histoArray[index].GetXaxis().SetBinLabel(7,"isNeutHadron")
                    histoArray[index].GetXaxis().SetBinLabel(8,"isUndefined")
                    histoArray[index].GetXaxis().SetBinLabel(9,"else")
                    histoArray[index].Draw("COLZ L")
                  elif (keyname2== "CutFlowEta" or keyname2 == "CutFlowPfType") :
                    histoArray[index].SetStats(0)
                    histoArray[index].GetYaxis().SetBinLabel(1,"Trigger")
                    histoArray[index].GetYaxis().SetBinLabel(2,"Eta")
                    histoArray[index].GetYaxis().SetBinLabel(3,"pT")
                    histoArray[index].GetYaxis().SetBinLabel(4,"NumHits")
                    histoArray[index].GetYaxis().SetBinLabel(5,"NumPixHits")
                    histoArray[index].GetYaxis().SetBinLabel(6,"ValidFract")
                    histoArray[index].GetYaxis().SetBinLabel(7,"NumDeDx")
                    histoArray[index].GetYaxis().SetBinLabel(8,"ProbXY")
                    histoArray[index].GetYaxis().SetBinLabel(9,"HighPurity")
                    histoArray[index].GetYaxis().SetBinLabel(10,"Chi2oDOF")
                    histoArray[index].GetYaxis().SetBinLabel(11,"EoP")
                    histoArray[index].GetYaxis().SetBinLabel(12,"dz")
                    histoArray[index].GetYaxis().SetBinLabel(13,"dxy")
                    histoArray[index].GetYaxis().SetBinLabel(14,"pTerrOverpT")
                    histoArray[index].GetYaxis().SetBinLabel(15,"TKIso")
                    histoArray[index].GetYaxis().SetBinLabel(16,"MiniIso")
                    histoArray[index].GetYaxis().SetBinLabel(17,"MassT")
                    histoArray[index].GetYaxis().SetBinLabel(18,"Ih")
                    histoArray[index].GetYaxis().SetBinLabel(19,"ProbQ")
                    histoArray[index].GetYaxis().SetBinLabel(20,"MuStat")
                    histoArray[index].GetYaxis().SetBinLabel(21,"PhiTOF")
                    histoArray[index].GetYaxis().SetBinLabel(22,"EtaTOF")
                    histoArray[index].Draw("COLZ")
                  elif ("IsPer" in keyname2) :
                    histoArray[index].SetStats(0)
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].ProjectionY(newname,histoArray[index].GetXaxis().FindBin(0.7),-1,"e").Draw("COLZ")
                  elif ("EIsolPer" in keyname2) :
                    histoArray[index].SetStats(0)
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].ProjectionY(newname,histoArray[index].GetXaxis().FindBin(0.85),-1,"e").Draw("COLZ")
                  else :
                    histoArray[index].SetMaximum(max*1.5)
                  histoArray[index].Draw("SAME")
                
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
                
                tex4 = ROOT.TLatex()
                if ("PrePreS" in keyname2) :
                  tex4 = ROOT.TLatex(0.6,0.95,"Before pre-selection")
                elif ("N1" in keyname2) :
                  tex4 = ROOT.TLatex(0.6,0.95,"After N-1 selection")
                elif ("PostPreS" in keyname2) :
                  tex4 = ROOT.TLatex(0.6,0.95,"After pre-selection")
                tex4.SetNDC();
                tex4.SetTextFont(52);
                tex4.SetTextSize(0.045);
                tex4.SetLineWidth(2);


#                stackPlots.Draw("HIST")
#                QCDFlatHisto.Draw("SAME")
#                stackPlots.SetTitle("")
#                stackPlots.GetXaxis().SetTitleSize(0.05)
#                stackPlots.GetXaxis().SetTitleOffset(1)
#                stackPlots.GetXaxis().SetTitle("Mass [GeV]")
#                stackPlots.GetYaxis().SetTitle("Tracks/bin")
#                stackPlots.GetYaxis().SetTitleSize(0.05)
#                stackPlots.GetYaxis().SetTitleOffset(1)

#                stackPlots.SetMaximum(max*1.1)
#                stackPlots.SetMinimum(0.001)
                legend.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")

                canvas.SaveAs("Compare"+sampleInFile[:-4]+"/"+keyname2+".png")

                cstackPlotsLogString = 'cstackPlotsLog'+str(j)
                cstackPlots = ROOT.TCanvas(cstackPlotsLogString, cstackPlotsLogString, 800,800)
                cstackPlots.SetLogy()
                for index in range(0, len(histoArray)):
                  histoArray[index].Draw("SAME")
                  if ("CutFlow" in keyname2 or "pfType"  in keyname2 ) :
                     histoArray[index].SetMaximum(100)
                     histoArray[index].SetMinimum(0.0001)
                  else :
                    histoArray[index].SetMaximum(max*10000)
                    histoArray[index].SetMinimum(0.000001)
                legend.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")

                cstackPlots.SaveAs("Compare"+sampleInFile[:-4]+"/"+keyname2+"_log.png")
              
