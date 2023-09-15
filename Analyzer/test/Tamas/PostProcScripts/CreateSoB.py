import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python3 %prog sample.txt")
(opt,args) = parser.parse_args()

sampleInFile = sys.argv[1]

backFile = ROOT.TFile.Open("crab_Analysis_2018_AllBackground_woProbQ_CodeV23p8_v1.root")

ROOT.gROOT.SetBatch(True)
#ROOT.gStyle.SetPalette(1)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

SamplesArray = []

with open(sampleInFile, "r") as a_file:
  for line in a_file:
    stripped_line = line.strip()
    SamplesArray.append(stripped_line)

fileInArray = []
for sample in SamplesArray:
  fileInArray.append(ROOT.TFile.Open(sample))
  
dirs = []
# loop on the outer directory
for i in range(0, fileInArray[0].GetListOfKeys().GetEntries()):
  dirname = fileInArray[0].GetListOfKeys().At(i).GetName()
  curr_dir = fileInArray[0].GetDirectory(dirname)
  if not (curr_dir) :
    continue
  # loop on the second level
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = fileInArray[0].GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      # loop on the third level, this is very the rel plots are
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          newname = dirname + "/" + keyname+ "/" + keyname2
          obj = fileInArray[0].Get(newname)
          if not obj.InheritsFrom("TObject"):
            continue
          if (obj.GetEntries() == 0 ) :
            continue
          if not ("N1" in keyname2):
            continue
          if ("N1_Stations" in keyname2) or ("N1_Qual" in keyname2) or ("pfType" in keyname2):
            continue
          if not os.path.exists(os.path.dirname("CompareSoB"+sampleInFile[:-4]+"/")): os.makedirs(os.path.dirname("CompareSoB"+sampleInFile[:-4]+"/"))
          
          if (obj.ClassName() == "TH1F" or obj.ClassName() == "TProfile"): # and "BS_" in keyname2):
            canvasString = 'canvas'+str(j)
            canvas = ROOT.TCanvas(canvasString, canvasString, 800,800)
            legend =  ROOT.TLegend(.45,.60,.80,.89,"","brNDC")
            legend.SetTextFont(42)
            legend.SetTextSize(0.017)
            legend.SetBorderSize(1);
            legend.SetBorderSize(0);
            legend.SetLineColor(1);
            legend.SetLineStyle(1);
            legend.SetLineWidth(1);
            legend.SetFillColor(0);
            legend.SetFillStyle(1001);
            
            backObj = backFile.Get(newname)
            histoArray = []
            SignalOverNoise = []
            for fileIn in fileInArray:
              histoArray.append(fileIn.Get(newname))
              SignalOverNoise.append(fileIn.Get(newname))
            for index in range(0, len(histoArray)):
              histoArray[index].GetNbinsX()
              SoBstring = 'SoB'+str(index)
              numBins = histoArray[index].GetNbinsX()
              maxXaxis = histoArray[index].GetXaxis().GetXmax()
              minXaxis = histoArray[index].GetXaxis().GetXmin()

              stepSize = (maxXaxis-minXaxis)/numBins
              SignalOverNoise[index] = ROOT.TH1F(SoBstring,SoBstring,numBins,minXaxis,maxXaxis)
              SignalOverNoise[index].SetStats(0)
              SignalOverNoise[index].SetMarkerStyle(20)
              
              legend.AddEntry(SignalOverNoise[index],SamplesArray[index][19:SamplesArray[index].find("Code")+9],"LP")
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
              SignalOverNoise[index].SetLineColor(indexNew)
              SignalOverNoise[index].SetMarkerColor(indexNew)
              SignalOverNoise[index].SetTitle("")
              SignalOverNoise[index].GetYaxis().SetTitle("Normalized S/#sqrt{S+B}")
              axisTitle = keyname2[3:]
              SignalOverNoise[index].GetXaxis().SetTitle(axisTitle)
              for x in numpy.arange(0,maxXaxis,stepSize):
                Signal = histoArray[index].Integral(histoArray[index].GetXaxis().FindBin(0),histoArray[index].GetXaxis().FindBin(x))
                Background = backObj.Integral(backObj.GetXaxis().FindBin(0),backObj.GetXaxis().FindBin(x))
                

#                else :
#                  Signal = numpy.absolute(histoArray[index].Integral(histoArray[index].GetXaxis().FindBin(minXaxis),histoArray[index].GetXaxis().FindBin(x)))
#                  Background = numpy.absolutebackObj.Integral(backObj.GetXaxis().FindBin(minXaxis),backObj.GetXaxis().FindBin(x)))
                if ((Signal+Background)!=0) :
                  SoB = Signal/numpy.sqrt(Signal+Background)
                else :
                  SoB = 0
                if ("Dxy" in keyname2) :
                  print("For ",x," the S = ",Signal," and B = ",Background, " and the S/sqrt(S+B)",Signal/numpy.sqrt(Signal+Background))
                
                SignalOverNoise[index].SetBinContent(histoArray[index].GetXaxis().FindBin(x),SoB)
              SignalOverNoise[index].Scale(1/SignalOverNoise[index].GetMaximum())
              SignalOverNoise[index].GetYaxis().SetRangeUser(0.,2)
              SignalOverNoise[index].Draw("SAME L")
                
#
#                  histoArray[index].SetMaximum(1.4)
            
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
            
            legend.Draw("SAME")
            tex2.Draw("SAME")
            tex3.Draw("SAME")

            canvas.SaveAs("CompareSoB"+sampleInFile[:-4]+"/"+keyname2+".png")
