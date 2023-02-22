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

codeVersion = sampleInFile[sampleInFile.find("Code")+5:sampleInFile.find("Code")+9]

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
          
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if ("Total" in keyname2 or "IntLumi" in keyname2 or "XSection" in keyname2) :
            continue
          if obj.InheritsFrom("TObject"):
              if not os.path.exists(os.path.dirname("Compare"+sampleInFile[:-4]+"/a.png")):
                print("Create dir")
                os.makedirs(os.path.dirname("Compare"+sampleInFile[:-4]+"/"))
#              print(keyname2)
#              if not ("EventCutFlow" in keyname2) : continue

              if (obj.GetEntries() == 0 ) : continue
              if (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
                obj.SetMarkerStyle(20)
                obj.GetXaxis().SetRange(bin,bin)
                obj.Project3D("ZY").Draw("COLZ")
              if (obj.ClassName() == "TH2F" and  (("Trigger" in keyname2) and ( "Vs" in keyname2))) :
                canvasString = 'canvas1'+str(j)
                canvas = ROOT.TCanvas(canvasString, canvasString, 800,800)
                legend =  ROOT.TLegend(.55,.65,.80,.89,"","brNDC")
                legend.SetTextFont(42)
                legend.SetTextSize(0.017)
                legend.SetBorderSize(1);
                legend.SetBorderSize(0);
                legend.SetLineColor(1);
                legend.SetLineStyle(1);
                legend.SetLineWidth(1);
                legend.SetFillColor(0);
                legend.SetFillStyle(1001);
                
                histoArrayProfY = []
                k = 0
                for fileIn in fileInArray:
                  k = k + 1
                  projXName = keyname2 + str(k)
                  if (fileIn.Get(newname)) : histoArrayProfY.append(fileIn.Get(newname).ProfileY().ProjectionX(projXName))
                indexNew = -1
                for index in range(0, len(histoArrayProfY)):
                  histoArrayProfY[index].SetStats(0)
                  histoArrayProfY[index].SetMarkerStyle(20)
                  legend.AddEntry(histoArrayProfY[index],SamplesArray[index][19:SamplesArray[index].find("Code")+9],"LP")
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
                  histoArrayProfY[index].SetLineColor(indexNew)
                  histoArrayProfY[index].SetMarkerColor(indexNew)
#                  histoArray[index].SetFillColor(indexNew)
                  histoArrayProfY[index].SetTitle("")
                  histoArrayProfY[index].GetYaxis().SetTitle("Efficiency")
                  histoArrayProfY[index].SetMaximum(1.6)
                  histoArrayProfY[index].DrawClone("SAME")
                legend.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                canvas.SaveAs("Compare"+sampleInFile[:-4]+"/"+keyname2+"_profileY.png")
              
              
              # Now do the 1D histos
              
              if (obj.ClassName() == "TH1F" or obj.ClassName() == "TProfile"): # and "BS_" in keyname2):
                canvasString = 'canvas'+str(j)
                canvas = ROOT.TCanvas(canvasString, canvasString, 800,800)
                legend =  ROOT.TLegend(.55,.65,.80,.89,"","brNDC")
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
                  if (fileIn.Get(newname)) : histoArray.append(fileIn.Get(newname))
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
                    if not (histoArray[index2]) : continue
                    max = numpy.maximum(max,histoArray[index2].GetMaximum())
#                  histoArray[index].GetYaxis().SetTitle("Tracks/bin")
#                  histoArray[index].GetXaxis().SetTitle(keyname2)
                  if (keyname2 == "pfType") :
                    if (histoArray[index].GetMaximum() > 0 ) :  histoArray[index].Scale(1/histoArray[index].GetMaximum())
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
                  elif (keyname2== "EventCutFlow") :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    SecondBinSize =  histoArray[0].GetBinContent(1) * 1.3
                    histoArray[index].GetYaxis().SetRangeUser(0.,SecondBinSize)
                    histoArray[index].SetMaximum(SecondBinSize)
                    ROOT.gStyle.SetPaintTextFormat(".2g");
                    histoArray[index].Draw("SAMEHISTOTEXT00")
                  elif (keyname2== "CutFlow") :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    if (histoArray[index].GetMaximum() > 0 ) :  histoArray[index].Scale(1/histoArray[index].GetMaximum())
                    SecondBinSize =  histoArray[0].GetBinContent(1) * 1.3
                    histoArray[index].GetYaxis().SetRangeUser(0.,SecondBinSize)
                    histoArray[index].SetMaximum(SecondBinSize)
                  elif (keyname2== "CutFlowReverse") :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    if (histoArray[index].GetMaximum() > 0 ) :  histoArray[index].Scale(1/histoArray[index].GetMaximum())
                    histoArray[index].GetYaxis().SetRangeUser(0.,1.6)
                    histoArray[index].SetMaximum(1.6)
                  elif ("TriggerType" in keyname2) :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    if (histoArray[index].Integral()>0) : histoArray[index].Scale(1/histoArray[index].Integral())
                    histoArray[index].GetXaxis().SetTitle("")
                    histoArray[index].GetXaxis().SetRangeUser(0.,5.)
                    histoArray[index].GetYaxis().SetRangeUser(0.01,1.)
                    histoArray[index].SetMaximum(1.)
                    ROOT.gStyle.SetPaintTextFormat(".2g");
                    histoArray[index].Draw("SAMEHISTOTEXT00")
                  elif ("HltMatchTrackLevel" in keyname2) :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    if (histoArray[index].GetBinContent(1)>0) : histoArray[index].Scale(1/histoArray[index].GetBinContent(1))
                    histoArray[index].SetBinError(1,0)
                    histoArray[index].SetBinError(2,0)
                    histoArray[index].GetYaxis().SetTitle("Norm events / category")
                    histoArray[index].GetYaxis().SetRangeUser(0.92,1.01)
                    ROOT.gStyle.SetPaintTextFormat(".2g");
                    histoArray[index].Draw("SAMEHISTOTEXT00")
                  elif ("pfType" in keyname2) :
                    histoArray[index].SetMarkerStyle(20)
                    histoArray[index].SetStats(0)
                    if (histoArray[index].GetMaximum() > 0 ) :  histoArray[index].Scale(1/histoArray[index].GetMaximum())
                    histoArray[index].Draw("COLZ L")
                  elif (keyname2== "CutFlowEta" or keyname2 == "CutFlowPfType") :
                    histoArray[index].SetStats(0)
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
                tex4.Draw("SAME")
                tex5.Draw("SAME")

                canvas.SaveAs("Compare"+sampleInFile[:-4]+"/"+keyname2+".png")

                cstackPlotsLogString = 'cstackPlotsLog'+str(j)
                cstackPlots = ROOT.TCanvas(cstackPlotsLogString, cstackPlotsLogString, 800,800)
                cstackPlots.SetLogy()
                for index in range(0, len(histoArray)):
                  histoArray[index].Draw("SAME")
                  if ("CutFlow" in keyname2 or "pfType"  in keyname2 ) :
                     histoArray[index].SetMaximum(100)
                     histoArray[index].SetMinimum(0.0001)
                  if (keyname2== "EventCutFlow") :
                     histoArray[index].SetStats(0)
                     if ("SingleMuon" in SamplesArray[index]) :
                       histoArray[index].SetMinimum(300000)
                       histoArray[index].SetMaximum(max*2)
                       print(str(SamplesArray[index][14:SamplesArray[index].find("Code")-1])+": "+str(histoArray[index].GetBinContent(17)))
                     else :
                      histoArray[index].SetMinimum(max*1000000)
                      histoArray[index].SetMinimum(0.000001)

                     
                     ROOT.gStyle.SetPaintTextFormat(".2g");
                     histoArray[index].Draw("SAMEHISTOTEXT25")
                  elif ("ProbQ" in keyname2) :
                    histoArray[index].SetMaximum(max*1000000)
                    histoArray[index].SetMinimum(0.000001)
                  else :
                    histoArray[index].SetMaximum(max*10000)
                    histoArray[index].SetMinimum(0.000001)
                legend.Draw("SAME")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")

                cstackPlots.SaveAs("Compare"+sampleInFile[:-4]+"/"+keyname2+"_log.png")
                
                cstackPlotsNormString = 'cstackPlotsNorm'+str(j)
                cstackPlotsNorm = ROOT.TCanvas(cstackPlotsNormString, cstackPlotsLogString, 800,800)
                maxNorm = 0.
                for index in range(0, len(histoArray)):
                  int = histoArray[index].Integral()
                  if ("CandidateType" in keyname2) :
                    int = histoArray[index].GetBinContent(1)
                  if (int>0) :
                    histoArray[index].Scale(1/int)
                  yAxisTitle = histoArray[index].GetYaxis().GetTitle()
                  histoArray[index].GetYaxis().SetTitle("Normalized " + yAxisTitle)
                  maxNorm = numpy.maximum(maxNorm,histoArray[index].GetMaximum())
                  histoArray[index].Draw("SAME")
                histoArray[0].SetMaximum(maxNorm*1.6)
                canvas.SaveAs("Compare"+sampleInFile[:-4]+"/"+keyname2+"_Norm.png")
              
os.system("cp forWebpage/* Compare"+sampleInFile[:-4]+"/.")
os.system("cp forWebpage/.htaccess Compare"+sampleInFile[:-4]+"/.")
print("scp -r Compare"+sampleInFile[:-4]+" tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/.")
