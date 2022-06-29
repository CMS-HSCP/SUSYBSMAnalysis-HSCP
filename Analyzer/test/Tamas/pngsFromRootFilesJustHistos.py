import ROOT, sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog fileName.root BinNumber")
(opt,args) = parser.parse_args()

ROOT.gROOT.SetStyle("Plain")
#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.SetBatch(True)

ROOT.gStyle.SetPadTopMargin(0.07);
ROOT.gStyle.SetPadBottomMargin(0.1);
ROOT.gStyle.SetPadLeftMargin(0.15);
ROOT.gStyle.SetPadRightMargin(0.13);

fileName = sys.argv[1]
BinNumber = sys.argv[2]

bin = int(BinNumber)
# bin 3: pt>60 and I_as > 0.05
# bin 25: pt>65 and I_as > 0.175
# bin 28: pt>65 and I_as > 0.3

print("Filename: "+fileName)
input_file = fileName
output_dir = "Plots/"

ProjBin = int(BinNumber)
newFileDir = fileName[0:-5] + "_Bin" + str(ProjBin)

if not os.path.exists(output_dir): os.mkdir(output_dir)

f = ROOT.TFile.Open(input_file)
fileOut = open("SignalBackgroundEff.txt", "a")

dirs = []
for i in range(0, f.GetListOfKeys().GetEntries()):
  # Remove/modify unnecessary stuff from the name of the plot that was required by SmartHistos to ditinguish plots
  dirname = f.GetListOfKeys().At(i).GetName()
  curr_dir = f.GetDirectory(dirname)
# print("dirname: "+dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      # Match the plot of interest
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = f.GetDirectory(dirname+"/"+keyname)
#                    print("keyname: "+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          # The plot should be TCanvas
          newname = dirname + "/" + keyname+ "/" + keyname2
#         print("newname: "+newname)
          obj = f.Get(newname)
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
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
#         print(obj.InheritsFrom())
          if obj.InheritsFrom("TObject"):
              can = obj
              can = ROOT.TCanvas(newname,newname,800,800)
              # Save plot
              name = fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  ".png"
#              print(keyname2)
              if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
              if (obj.GetEntries() == 0 ) :
                continue
#                 print(obj.ClassName())
              if (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
                obj.SetStats(0)
                obj.SetTitle("")
                if ("ProbQPerProbXY" in keyname2) :
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(3.22),-1)
                  obj.GetYaxis().SetRange(obj.GetYaxis().FindBin(0.0),obj.GetYaxis().FindBin(0.1))
                  obj.Project3D("YZ").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIh.png")
                  obj.GetXaxis().UnZoom()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(3.22))
                  obj.GetYaxis().SetRange(obj.GetYaxis().FindBin(0.0),obj.GetYaxis().FindBin(0.1))
                  obj.Project3D("YZ").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIh.png")
                  obj.GetXaxis().UnZoom()
                  obj.Project3D("YZ").Draw("COLZ")
                  can.SaveAs(name)
                # this maybe should go from bin to bin+1 ?
                if ("IhPerLayer" in keyname2) :
                  obj.SetMarkerStyle(20)
                  if ("PostPreS_MIsPixelIhPerLayer" in keyname2) :
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(1,"BPix L1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(2,"BPix L2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(3,"BPix L3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(4,"BPix L4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(5,"FPix D1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(6,"FPix D2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(7,"FPix D3")
                      obj.Project3D("YZ").GetYaxis().SetTitleOffset(0.8)
                      obj.Project3D("YZ").GetYaxis().SetTitle("Ih")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(1,"BPix L1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(2,"BPix L2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(3,"BPix L3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(4,"BPix L4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(5,"FPix D1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(6,"FPix D2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(7,"FPix D3")
                      obj.Project3D("XZ").GetYaxis().SetTitle("Ias")
                      obj.Project3D("XZ").GetYaxis().SetTitleOffset(0.8)
                  elif ("PostPreS_MIsStripIhPerLayer" in keyname2) :
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(1,"TIB L1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(2,"TIB L2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(3,"TIB L3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(4,"TIB L4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(5,"TOB L1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(6,"TOB L2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(7,"TOB L3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(8,"TOB L4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(9,"TOB L5")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(10,"TOB L6")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(11,"TID D1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(12,"TID D2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(13,"TID D3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(14,"TEC D1")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(15,"TEC D2")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(16,"TEC D3")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(17,"TEC D4")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(18,"TEC D5")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(19,"TEC D6")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(20,"TEC D7")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(21,"TEC D8")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(22,"TEC D9")
                      obj.Project3D("YZ").GetXaxis().SetBinLabel(23,"TEC D10")
                      obj.Project3D("YZ").GetYaxis().SetTitleOffset(0.8)
                      obj.Project3D("YZ").GetYaxis().SetTitle("Ih")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(1,"TIB L1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(2,"TIB L2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(3,"TIB L3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(4,"TIB L4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(5,"TOB L1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(6,"TOB L2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(7,"TOB L3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(8,"TOB L4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(9,"TOB L5")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(10,"TOB L6")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(11,"TID D1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(12,"TID D2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(13,"TID D3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(14,"TEC D1")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(15,"TEC D2")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(16,"TEC D3")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(17,"TEC D4")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(18,"TEC D5")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(19,"TEC D6")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(20,"TEC D7")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(21,"TEC D8")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(22,"TEC D9")
                      obj.Project3D("XZ").GetXaxis().SetBinLabel(23,"TEC D10")
                      obj.Project3D("XZ").GetYaxis().SetTitle("Ias")
                      obj.Project3D("XZ").GetYaxis().SetTitleOffset(0.8)
#                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(1.0))
                  obj.Project3D("YZ").SetStats(0)
                  obj.Project3D("YZ").Draw("COLZ")
                  obj.Project3D("YZ").SetTitle("")
                  tex4.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_allIs.png")
                  obj.Project3D("XZ").SetStats(0)
                  obj.Project3D("XZ").Draw("COLZ")
                  obj.Project3D("XZ").SetTitle("")
                  tex4.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_allIh.png")
#                  obj.Project3D("YZ").GetYaxis().SetRange(1,obj.Project3D("YZ").GetYaxis().FindBin(5.0))
#                  lowIhPart.obj.Project3D("YZ").Draw()
#                  obj.Project3D("YZ").GetYaxis().SetRange(obj.Project3D("YZ").GetYaxis().FindBin(5.0),-1)
#                  highIhPart = obj.Project3D("YZ").Draw()
#                  ratio = lowIhPart.Multiply(1/highIhPart)
#                  ratio.Draw()
                  tex4.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_RatioOfLowIsHighIas.png")
                  
                  obj.Project3D("YZ").GetYaxis().UnZoom()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.7),obj.GetXaxis().FindBin(1.0))
                  obj.Project3D("YZ").Draw("COLZ")
                  tex4.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIs.png")
                  
                  obj.GetXaxis().UnZoom()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(0.7))
                  obj.Project3D("YZ").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIs.png")
                  if ("PostPreS_MIsPixelIhPerLayer" in keyname2):
                    for i in range(7) :
                      obj.Project3D("YZ").ProjectionY(newname,i,i+1,"e").Draw()
                      can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_PixLayer"+str(i)+".png")
                  elif ("PostPreS_MIsStripIhPerLayer" in keyname2):
                    for i in range(23) :
                      obj.Project3D("YZ").ProjectionY(newname,i,i+1,"e").Draw()
                      can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_StripLayer"+str(i)+".png")
                else :
                  obj.SetMarkerStyle(20)
                  obj.GetXaxis().SetRange(bin,bin)
                  obj.Project3D("ZY").Draw("COLZ") # was Project3DProfile
              elif ((obj.ClassName() == "TH2F" or obj.ClassName() == "TH2D") and not (keyname2 == "GenPtVsRecoPt" or "PreS_" in keyname2 or keyname2 == "CutFlowEta" or  keyname2 == "CutFlowPfType"  or "N1_" in keyname2)):
                obj.SetMarkerStyle(20)
                obj.ProjectionY(newname,bin,bin,"e").Draw("COLZ")
              elif ("GenPtVsRecoPt" in keyname2) :
                tex4.Draw("SAME")
                obj.Draw("COLZ")
                obj.SetMinimum(0.000001)
                can.SetLogz();
                print(str(keyname2) + ": " + str(round(obj.GetCorrelationFactor(),2)))
              elif ("IsPer" in keyname2) :
                obj.SetMarkerStyle(20)
                if ("Angle" in keyname2 or "NumSibling" in keyname2) :
                  tex4.Draw("SAME")
                  obj.ProjectionY(newname,obj.GetXaxis().FindBin(0.7),-1,"e").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
                  obj.ProjectionY(newname,obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(0.7),"e").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
                  obj.ProjectionY(newname,obj.GetXaxis().FindBin(0.0),-1,"e").Draw("COLZ")
                else :
                  myPie = ROOT.TPie(obj.ProjectionY(newname+"_lowIas",1,obj.GetXaxis().FindBin(0.7),"e"))
                  myPie.SetLabelFormat("%txt (%perc)")
                  myPie.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  myPie.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
                  
                  myPie2 = ROOT.TPie(obj.ProjectionY(newname+"_highIas",obj.GetXaxis().FindBin(0.7),obj.GetNbinsX(),"e"))
                  myPie2.SetLabelFormat("%txt (%perc)")
                  myPie2.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  myPie2.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
                  
                  myPie3 = ROOT.TPie(obj.ProjectionY(newname,1,-1,"e"))
                  myPie3.SetLabelFormat("%txt (%perc)")
                  myPie3.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  myPie3.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  ".png")
              elif ("EIsolPer" in keyname2) :
                obj.SetMarkerStyle(20)
                if ("Angle" in keyname2) :
                  obj.ProjectionY(newname+"_lowEoP",1,obj.GetXaxis().FindBin(0.85),"e").Draw("COLZ")
                  tex4.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowEoP.png")
#                  obj.ProjectionY(newname+"_highEoP",obj.GetXaxis().FindBin(0.85),obj.GetNbinsX(),"e").Draw("COLZ")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highEoP.png")
                elif ("EIsolPerPfType" in keyname2) :
                  obj.SetMarkerColor(1)
                  obj.SetLineColor(1)
                  obj.SetMarkerStyle(20)
                  obj.SetStats(0)
#                  obj.Scale(1/obj.GetMaximum())
                  obj.GetYaxis().SetBinLabel(1,"AllTracks")
                  obj.GetYaxis().SetBinLabel(2,"PFtracks")
                  obj.GetYaxis().SetBinLabel(3,"isElectron")
                  obj.GetYaxis().SetBinLabel(4,"isMuon")
                  obj.GetYaxis().SetBinLabel(5,"isPhoton")
                  obj.GetYaxis().SetBinLabel(6,"isChHadron")
                  obj.GetYaxis().SetBinLabel(7,"isNeutHadron")
                  obj.GetYaxis().SetBinLabel(8,"isUndefined")
                  obj.GetYaxis().SetBinLabel(9,"notPFtrack")
                  obj.GetXaxis().SetTitle("EoP")
                  obj.Draw("COLZ")
                else :
                  myPie = ROOT.TPie(obj.ProjectionY(newname+"_lowEoP",1,obj.GetXaxis().FindBin(0.85),"e"))
                  myPie.SetLabelFormat("%txt (%perc)")
                  myPie.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  myPie.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowEoP.png")
                  
#                  myPie2 = ROOT.TPie(obj.ProjectionY(newname+"_highEoP",obj.GetXaxis().FindBin(0.85),obj.GetNbinsX(),"e"))
#                  myPie2.SetLabelFormat("%txt (%perc)")
#                  myPie2.SetLabelsOffset(-.27)
#                  tex4.Draw("SAME")
#                  myPie2.Draw("R<")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highEoP.png")
              elif ("ProbQVsIas" in keyname2) :
                print(str(keyname2) + ": " + str(round(obj.GetCorrelationFactor(),2)))
              elif (keyname2== "CutFlow") :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetXaxis().SetBinLabel(1,"Trigger")
                obj.GetXaxis().SetBinLabel(2,"Eta")
                obj.GetXaxis().SetBinLabel(3,"pT")
                obj.GetXaxis().SetBinLabel(4,"NumHits")
                obj.GetXaxis().SetBinLabel(5,"NumPixHits")
                obj.GetXaxis().SetBinLabel(6,"ValidFract")
                obj.GetXaxis().SetBinLabel(7,"NumDeDx")
                obj.GetXaxis().SetBinLabel(8,"ProbXY")
                obj.GetXaxis().SetBinLabel(9,"HighPurity")
                obj.GetXaxis().SetBinLabel(10,"Chi2oDOF")
                obj.GetXaxis().SetBinLabel(11,"EoP")
                obj.GetXaxis().SetBinLabel(12,"dz")
                obj.GetXaxis().SetBinLabel(13,"dxy")
                obj.GetXaxis().SetBinLabel(14,"pTerrOverpT")
                obj.GetXaxis().SetBinLabel(15,"SVfromNI")
                obj.GetXaxis().SetBinLabel(16,"MiniIso")
                obj.GetXaxis().SetBinLabel(17,"PFid")
                obj.GetXaxis().SetBinLabel(18,"Ih")
                obj.GetXaxis().SetBinLabel(19,"ProbQ")
                obj.GetXaxis().SetBinLabel(20,"MuStat")
                obj.GetXaxis().SetBinLabel(21,"PhiTOF")
                obj.GetXaxis().SetBinLabel(22,"EtaTOF")
                tex4.Draw("SAME")
                obj.Draw("COLZ L")
              elif (keyname2== "CutFlowProbQFirst") :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetXaxis().SetBinLabel(1,"Trigger")
                obj.GetXaxis().SetBinLabel(2,"Eta")
                obj.GetXaxis().SetBinLabel(3,"pT")
                obj.GetXaxis().SetBinLabel(4,"ProbQ")
                obj.GetXaxis().SetBinLabel(5,"NumHits")
                obj.GetXaxis().SetBinLabel(6,"NumPixHits")
                obj.GetXaxis().SetBinLabel(7,"ValidFract")
                obj.GetXaxis().SetBinLabel(8,"NumDeDx")
                obj.GetXaxis().SetBinLabel(9,"ProbXY")
                obj.GetXaxis().SetBinLabel(10,"HighPurity")
                obj.GetXaxis().SetBinLabel(11,"Chi2oDOF")
                obj.GetXaxis().SetBinLabel(12,"EoP")
                obj.GetXaxis().SetBinLabel(13,"dz")
                obj.GetXaxis().SetBinLabel(14,"dxy")
                obj.GetXaxis().SetBinLabel(15,"pTerrOverpT")
                obj.GetXaxis().SetBinLabel(16,"SVfromNI")
                obj.GetXaxis().SetBinLabel(17,"MiniIso")
                obj.GetXaxis().SetBinLabel(18,"PFid")
                obj.GetXaxis().SetBinLabel(19,"Ih")
                obj.GetXaxis().SetBinLabel(20,"MuStat")
                obj.GetXaxis().SetBinLabel(21,"PhiTOF")
                obj.GetXaxis().SetBinLabel(22,"EtaTOF")
                tex4.Draw("SAME")
                obj.Draw("COLZ L")
              elif ("_pfType" in keyname2) :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetXaxis().SetBinLabel(1,"AllTracks")
                obj.GetXaxis().SetBinLabel(2,"PFtracks")
                obj.GetXaxis().SetBinLabel(3,"isElectron")
                obj.GetXaxis().SetBinLabel(4,"isMuon")
                obj.GetXaxis().SetBinLabel(5,"isPhoton")
                obj.GetXaxis().SetBinLabel(6,"isChHadron")
                obj.GetXaxis().SetBinLabel(7,"isNeutHadron")
                obj.GetXaxis().SetBinLabel(8,"isUndefined")
                obj.GetXaxis().SetBinLabel(9,"notPFtrack")
                tex4.Draw("SAME")
                obj.Draw("COLZ L")
              elif ((keyname2 == "CutFlowEta") or (keyname2 == "CutFlowPfType")) :
                obj.SetStats(0)
                obj.GetYaxis().SetBinLabel(1,"Trigger")
                obj.GetYaxis().SetBinLabel(2,"Eta")
                obj.GetYaxis().SetBinLabel(3,"pT")
                obj.GetYaxis().SetBinLabel(4,"NumHits")
                obj.GetYaxis().SetBinLabel(5,"NumPixHits")
                obj.GetYaxis().SetBinLabel(6,"ValidFract")
                obj.GetYaxis().SetBinLabel(7,"NumDeDx")
                obj.GetYaxis().SetBinLabel(8,"ProbXY")
                obj.GetYaxis().SetBinLabel(9,"HighPurity")
                obj.GetYaxis().SetBinLabel(10,"Chi2oDOF")
                obj.GetYaxis().SetBinLabel(11,"EoP")
                obj.GetYaxis().SetBinLabel(12,"dz")
                obj.GetYaxis().SetBinLabel(13,"dxy")
                obj.GetYaxis().SetBinLabel(14,"pTerrOverpT")
                obj.GetYaxis().SetBinLabel(15,"SVfromNI")
                obj.GetYaxis().SetBinLabel(16,"MiniIso")
                obj.GetYaxis().SetBinLabel(17,"PfID")
                obj.GetYaxis().SetBinLabel(18,"Ih")
                obj.GetYaxis().SetBinLabel(19,"ProbQ")
                obj.GetYaxis().SetBinLabel(20,"MuStat")
                obj.GetYaxis().SetBinLabel(21,"PhiTOF")
                obj.GetYaxis().SetBinLabel(22,"EtaTOF")
                if (keyname2 == "CutFlowPfType"):
                  obj.Scale(1/obj.GetMaximum())
                  obj.GetXaxis().SetBinLabel(1,"AllTracks")
                  obj.GetXaxis().SetBinLabel(2,"PFtracks")
                  obj.GetXaxis().SetBinLabel(3,"isElectron")
                  obj.GetXaxis().SetBinLabel(4,"isMuon")
                  obj.GetXaxis().SetBinLabel(5,"isPhoton")
                  obj.GetXaxis().SetBinLabel(6,"isChHadron")
                  obj.GetXaxis().SetBinLabel(7,"isNeutHadron")
                  obj.GetXaxis().SetBinLabel(8,"isUndefined")
                  obj.GetXaxis().SetBinLabel(9,"notPFtrack")
                tex4.Draw("SAME")
                obj.Draw("COLZ")
              else :
                obj.SetMarkerStyle(20)
#                tex4 = ROOT.TLatex()
#                if ("PrePreS" in keyname2) :
#                  tex4 = ROOT.TLatex(0.6,0.95,"Before pre-selection")
#                elif ("N1" in keyname2) :
#                  tex4 = ROOT.TLatex(0.6,0.95,"After N-1 selection")
#                elif ("PostPreS" in keyname2) :
#                  tex4 = ROOT.TLatex(0.6,0.95,"After pre-selection")
#                tex4.SetNDC();
#                tex4.SetTextFont(52);
#                tex4.SetTextSize(0.045);
#                tex4.SetLineWidth(2);
                
                if (keyname2.find("Per")==-1) :
                  axisXTitle = keyname2[keyname2.find("_")+1:]
                  axisYTitle = "Yield/bin"
                else :
                  axisXTitle = keyname2[keyname2.find("_")+1:keyname2.find("Per")]
                  axisYTitle = keyname2[keyname2.find("Per")+3:]
            
                obj.GetXaxis().SetTitle(axisXTitle)
                obj.GetYaxis().SetTitle(axisYTitle)
                obj.Draw("COLZ")
                tex4.Draw("SAME")
#                print(obj.GetListOfFunctions())
#                palette = ROOT.TPaletteAxis(obj.GetListOfFunctions().FindObject("palette"))
#                palette->SetX1NDC(0.9);
#                palette->SetX2NDC(0.95);
#                palette->SetY1NDC(0.2);
#                palette->SetY2NDC(0.8);
#                palette.SetX1NDC(0.6);
#                palette.SetX2NDC(0.65);
#                palette.SetY1NDC(0.2);
#                palette.SetY2NDC(0.8);
              can.SaveAs(name)
              
              if ("Angle" in keyname2 and obj.ClassName() == "TH2F" ) :
                obj.Draw("COLZ")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_2d.png")
              
              # now let's plot everything in a logy scale
              obj.SetMarkerStyle(20)
              obj.SetMinimum(0.000001)
#              obj.SetMaximum(10000)
              
              if (keyname2=="CutFlowPfType" or keyname2=="CutFlowEta" or "PostPreS_MIsPixelIhPerLayer" in keyname2 or "PostPreS_MIsStripIhPerLayer" in keyname2) :
                obj.SetMaximum(obj.GetMaximum())
                can.SetLogz()
              elif ("GenID" in keyname2) :
                continue
              else :
                obj.SetMaximum(obj.GetMaximum()*100)
                can.SetLogy()
                obj.SetStats(0)
                obj.SetTitle("")
              can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_logy.png")
#              if (keyname2 == "BS_ProbQ") :
#                firstBinInt = obj.Integral(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(0.1))
#                totalInt = obj.Integral(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(1.0))
#                efficiency = firstBinInt/float(totalInt)
#                fileOut.write(name+": "+str(efficiency)+"\n")
#                fileOut.close()
              #can.SaveAs(name.replace(".png",".pdf"))
              #can.SaveAs(name.replace(".png",".C"))
              can.Close()

          else:
              print(keyname+"   "+newname + " does not inherit from TObject" )

Mass = f.Get("/analyzer/BaseName/Mass")
Mass_wPred = f.Get("/analyzer/BaseName/Pred_Mass_CB")
if Mass_wPred :
  name = fileName[0:-5] + "_Bin" + str(bin)+ "/"
  if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
  massBins = [10.,50.,100.,150.,200.,250.,300.,350.,400.,500.,600.,700.,800.,900.,1000.,4000.]
  massBinsArray = np.array(massBins)
  Mass_projY_NotRebinned = Mass.ProjectionY("Mass_projY_NotRebinned",ProjBin,ProjBin,"e")
  Mass_wPred_projY_NotRebinned = Mass_wPred.ProjectionY("Mass_wPred_projY_NotRebinned",ProjBin,ProjBin,"e")

  Mass_projY = ROOT.TH1F("Mass_projY" , "Mass_projY" , len(massBinsArray)-1, massBinsArray)
  Mass_wPred_projY = ROOT.TH1F("Mass_wPred_projY" , "Mass_wPred_projY" , len(massBinsArray)-1, massBinsArray)

  print("Mass_projY_NotRebinned.Integral(): ",Mass_projY_NotRebinned.Integral())
  print("Mass_wPred_projY_NotRebinned.Integral(): ",Mass_wPred_projY_NotRebinned.Integral())

  KSvalue = Mass_projY_NotRebinned.KolmogorovTest(Mass_wPred_projY_NotRebinned,"XD")
  print("KS-test: "+str(KSvalue))

  for i in range(1,len(massBinsArray)) :
    voltBin = Mass_projY_NotRebinned.FindBin(massBinsArray[i-1])+1
    currentBin = Mass_projY_NotRebinned.FindBin(massBinsArray[i])
    Mass_projYCont = 0.0
    Mass_wPred_projYCont = 0.0
    Mass_projYCont_err2 = 0.0
    Mass_wPred_projYCont_err2 = 0.0
    for j in range(voltBin,currentBin) :
      Mass_projYCont += Mass_projY_NotRebinned.GetBinContent(j)
      Mass_projYCont_err2 += (Mass_projY_NotRebinned.GetBinError(j) * Mass_projY_NotRebinned.GetBinError(j))
      
      Mass_wPred_projYCont += Mass_wPred_projY_NotRebinned.GetBinContent(j)
      Mass_wPred_projYCont_err2 += (Mass_wPred_projY_NotRebinned.GetBinError(j)*Mass_wPred_projY_NotRebinned.GetBinError(j))
    Mass_projY.SetBinContent(i,Mass_projYCont)
    Mass_projY.SetBinError(i,np.sqrt(Mass_projYCont_err2))
    Mass_wPred_projY.SetBinContent(i,Mass_wPred_projYCont)
    Mass_wPred_projY.SetBinError(i,np.sqrt(Mass_wPred_projYCont_err2))
    

  print("----------------------------------------------")
  KSvalue2 = Mass_projY.KolmogorovTest(Mass_wPred_projY,"XD")
  print("KS-test after rebinning: "+str(KSvalue2))

  Mass_projY.SetMarkerColor(1)
  Mass_projY.SetLineColor(1)
  Mass_projY.SetMarkerStyle(20)
  Mass_projY.SetTitle("")
  Mass_projY.GetXaxis().SetTitleSize(0.05)
  Mass_projY.GetXaxis().SetTitleOffset(1)
  Mass_projY.GetXaxis().SetTitle("Mass [GeV]")
  Mass_projY.GetYaxis().SetTitle("Tracks/bin")
  Mass_projY.GetYaxis().SetTitleSize(0.05)
  Mass_projY.GetYaxis().SetTitleOffset(1)
  Mass_projY.SetStats(0)
#  Mass_projY.GetYaxis().SetRangeUser(0.001,Mass_projY.GetMaximum())


  Mass_wPred_projY.SetMarkerColor(2)
  Mass_wPred_projY.SetLineColor(2)
  Mass_wPred_projY.SetMarkerStyle(20)
  Mass_wPred_projY.SetTitle("")
  Mass_wPred_projY.GetXaxis().SetTitleSize(0.05)
  Mass_wPred_projY.GetXaxis().SetTitleOffset(1)
  Mass_wPred_projY.GetXaxis().SetTitle("Mass [GeV]")
  Mass_wPred_projY.GetYaxis().SetTitle("Tracks/bin")
  Mass_wPred_projY.GetYaxis().SetTitleSize(0.05)
  Mass_wPred_projY.GetYaxis().SetTitleOffset(1)
  Mass_wPred_projY.SetStats(0)
  
  Mass_projY_NotRebinned.SetMarkerColor(1)
  Mass_projY_NotRebinned.SetLineColor(1)
  Mass_projY_NotRebinned.SetMarkerStyle(20)
  Mass_projY_NotRebinned.SetTitle("")
  Mass_projY_NotRebinned.GetXaxis().SetTitleSize(0.05)
  Mass_projY_NotRebinned.GetXaxis().SetTitleOffset(1)
  Mass_projY_NotRebinned.GetXaxis().SetTitle("Mass [GeV]")
  Mass_projY_NotRebinned.GetYaxis().SetTitle("Tracks/bin")
  Mass_projY_NotRebinned.GetYaxis().SetTitleSize(0.05)
  Mass_projY_NotRebinned.GetYaxis().SetTitleOffset(1)
  Mass_projY_NotRebinned.SetStats(0)
#  Mass_projY_NotRebinned.GetYaxis().SetRangeUser(0.001,Mass_projY_NotRebinned.GetMaximum())


  Mass_wPred_projY_NotRebinned.SetMarkerColor(2)
  Mass_wPred_projY_NotRebinned.SetLineColor(2)
  Mass_wPred_projY_NotRebinned.SetMarkerStyle(20)
  Mass_wPred_projY_NotRebinned.SetTitle("")
  Mass_wPred_projY_NotRebinned.GetXaxis().SetTitleSize(0.05)
  Mass_wPred_projY_NotRebinned.GetXaxis().SetTitleOffset(1)
  Mass_wPred_projY_NotRebinned.GetXaxis().SetTitle("Mass [GeV]")
  Mass_wPred_projY_NotRebinned.GetYaxis().SetTitle("Tracks/bin")
  Mass_wPred_projY_NotRebinned.GetYaxis().SetTitleSize(0.05)
  Mass_wPred_projY_NotRebinned.GetYaxis().SetTitleOffset(1)
  Mass_wPred_projY_NotRebinned.SetStats(0)


  print("Mass_projY.Integral(): ",Mass_projY.Integral())
  print("Mass_wPred_projY.Integral(): ",Mass_wPred_projY.Integral())

  legMass =  ROOT.TLegend(.45,.75,.80,.9,"","brNDC")
  legMass.SetTextFont(42)
  legMass.SetTextSize(0.035)
  legMass.SetBorderSize(1);
  legMass.SetLineColor(1);
  legMass.SetLineStyle(1);
  legMass.SetLineWidth(1);
  legMass.SetFillColor(0);
  legMass.SetFillStyle(1001);
  legMass.AddEntry(Mass_wPred_projY,"Prediction","LP")
  legMass.AddEntry(Mass_projY,"Observation","LP")

  tex2 = ROOT.TLatex(0.13,0.94,"CMS");
  #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
  tex2.SetNDC();
  tex2.SetTextFont(61);
  tex2.SetTextSize(0.0675);
  tex2.SetLineWidth(2);

  #tex3 = ROOT.TLatex(0.27,0.96,"Simulation"); # for square plots
  #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
  tex3 = ROOT.TLatex(0.24,0.94,"Internal");
  tex3.SetNDC();
  tex3.SetTextFont(52);
  tex3.SetTextSize(0.0485);
  tex3.SetLineWidth(2);

  tex4 = ROOT.TLatex(0.7,0.93,"K-S test v2: "+str(round(KSvalue2,4)));
  tex4.SetNDC();
  tex4.SetTextFont(52);
  tex4.SetTextSize(0.0485);
  tex4.SetLineWidth(2);

  cMass_projY = ROOT.TCanvas('cMass_projY', 'cMass_projY',800,800)

  rp = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY)

  rp.SetH1DrawOpt("P");
  rp.SetH2DrawOpt("P");

  rp.Draw()
  #rp.GetUpperPad().BuildLegend()
  rp.SetLeftMargin(0.13);
  rp.SetRightMargin(0.05);
  rp.SetUpTopMargin(0.1);
  rp.SetLowTopMargin(0.02);
  rp.SetLowBottomMargin(0.35);

  max = Mass_projY.GetMaximum()*1.2
  Mass_projY.SetMaximum(max);
  rp.GetLowerRefGraph().SetMinimum(0.01);
  rp.GetLowerRefGraph().SetMaximum(2);
  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp.GetLowerRefGraph().SetLineColor(0) #0
  rp.GetLowerRefGraph().SetMarkerStyle(20)
  rp.GetLowerRefGraph().SetMarkerSize(1);
  rp.GetLowYaxis().SetNdivisions(505);
  rp.GetLowerRefYaxis().SetTitle("Ratio");
  rp.GetLowerRefYaxis().SetTitleSize(0.05);
  rp.GetLowerRefYaxis().SetTitleOffset(1);
  rp.GetLowerRefYaxis().SetLabelSize(0.035);


  rp.GetLowerRefXaxis().SetTitleSize(0.05);
  rp.GetLowerRefXaxis().SetTitleOffset(0.8);
  rp.GetLowerRefXaxis().SetLabelSize(0.035);
  cMass_projY.Modified()
  cMass_projY.Update()
  #Mass_projY.Draw()
  #Mass_wPred_projY.Draw("SAME")
  #rp.Draw("X")

  rp.GetUpperPad().cd();
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  
  name = newFileDir + "/cMass.png"
  cMass_projY.SaveAs(name)
  
  
  #############################################################################
  cMass_projY_log = ROOT.TCanvas('cMass_projY_log', 'cMass_projY_log',800,800)
  cMass_projY_log.SetLogy()

  rp2 = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY)

  rp2.SetH1DrawOpt("P");
  rp2.SetH2DrawOpt("P");

  rp2.Draw()
  #rp2.GetUpperPad().BuildLegend()
  rp2.SetLeftMargin(0.13);
  rp2.SetRightMargin(0.05);
  rp2.SetUpTopMargin(0.1);
  rp2.SetLowTopMargin(0.02);
  rp2.SetLowBottomMargin(0.35);

  max2 = np.maximum(Mass_projY.GetMaximum()*10,10000)
  Mass_projY.SetMaximum(max2);
#  Mass_projY.GetYaxis().SetRangeUser(0.1,100)
  rp2.GetLowerRefGraph().SetMinimum(0.01);
  rp2.GetLowerRefGraph().SetMaximum(2);
  #rp2.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp2.GetLowerRefGraph().SetLineColor(0) #0
  rp2.GetLowerRefGraph().SetMarkerStyle(20)
  rp2.GetLowerRefGraph().SetMarkerSize(1);
  rp2.GetLowYaxis().SetNdivisions(505);
  rp2.GetLowerRefYaxis().SetTitle("Ratio");
  rp2.GetLowerRefYaxis().SetTitleSize(0.05);
  rp2.GetLowerRefYaxis().SetTitleOffset(1);
  rp2.GetLowerRefYaxis().SetLabelSize(0.035);


  rp2.GetLowerRefXaxis().SetTitleSize(0.05);
  rp2.GetLowerRefXaxis().SetTitleOffset(0.8);
  rp2.GetLowerRefXaxis().SetLabelSize(0.035);
  cMass_projY_log.Modified()
  cMass_projY_log.Update()
  #Mass_projY.Draw()
  #Mass_wPred_projY.Draw("SAME")
  #rp2.Draw("X")

  rp2.GetUpperPad().cd();
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  
  name = newFileDir + "/cMass_log.png"
  cMass_projY_log.SaveAs(name)
  
  cMassOrig_projY = ROOT.TCanvas('cMassOrig_projY', 'cMassOrig_projY',800,800)
  rp0 = ROOT.TRatioPlot(Mass_projY_NotRebinned,Mass_wPred_projY_NotRebinned)
  rp0.Draw()
  rp0.SetH1DrawOpt("P");
  rp0.SetH2DrawOpt("P");
  rp0.SetLeftMargin(0.13);
  rp0.SetRightMargin(0.05);
  rp0.SetUpTopMargin(0.1);
  rp0.SetLowTopMargin(0.02);
  rp0.SetLowBottomMargin(0.35);
  rp0.GetLowerRefGraph().SetMinimum(0.01);
  rp0.GetLowerRefGraph().SetMaximum(2);
  #rp0.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp0.GetLowerRefGraph().SetLineColor(0) #0
  rp0.GetLowerRefGraph().SetMarkerStyle(20)
  rp0.GetLowerRefGraph().SetMarkerSize(1);
  rp0.GetLowYaxis().SetNdivisions(505);
  rp0.GetLowerRefYaxis().SetTitle("Ratio");
  rp0.GetLowerRefYaxis().SetTitleSize(0.05);
  rp0.GetLowerRefYaxis().SetTitleOffset(1);
  rp0.GetLowerRefYaxis().SetLabelSize(0.035);


  rp0.GetLowerRefXaxis().SetTitleSize(0.05);
  rp0.GetLowerRefXaxis().SetTitleOffset(0.8);
  rp0.GetLowerRefXaxis().SetLabelSize(0.035);
  cMassOrig_projY.Modified()
  cMassOrig_projY.Update()

  rp0.GetUpperPad().cd();
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  cMassOrig_projY.SaveAs(newFileDir + "/cMass_NotRebinned.png")


os.system("cp forWebpage/* "+newFileDir+"/.")
print("scp -r "+ newFileDir + " tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV18/.")
