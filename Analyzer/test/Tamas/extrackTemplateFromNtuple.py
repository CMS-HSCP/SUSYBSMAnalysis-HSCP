import ROOT, sys, os, time, re
import numpy as np
from ctypes import c_double as double
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

blind = True
#blind = False

print("Filename: "+fileName)
input_file = fileName

ProjBin = int(BinNumber)
newFileDir = fileName[0:-5] + "_Bin" + str(ProjBin)


f = ROOT.TFile.Open(input_file)
fileOut = open("SignalBackgroundEff.txt", "a")

isData = False
if ("SingleMuon" in fileName) : isData = True

iDontWannaRunPlots = False
#iDontWannaRunPlots = True


dirs = []
for i in range(0, f.GetListOfKeys().GetEntries()):
  # Remove/modify unnecessary stuff from the name of the plot that was required by SmartHistos to ditinguish plots
  dirname = f.GetListOfKeys().At(i).GetName()
  curr_dir = f.GetDirectory(dirname)
#  print("dirname: "+dirname)
  if (True):
      # Match the plot of interest
      keyname = f.GetListOfKeys().At(i).GetName()
      keyname2 = f.GetListOfKeys().At(i).GetName()
      curr_dir2 = f.GetDirectory(dirname)
#      print("keyname: "+keyname)
      if (True):
#          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          if ("__" in keyname2) : continue
          # The plot should be TCanvas
          newname = dirname
#          print("newname: "+newname)
          obj = f.Get(newname)
          
          
          tex2 = ROOT.TLatex(0.13,0.94,"CMS");
          #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

          #tex3 = ROOT.TLatex(0.27,0.96,"Simulation"); # for square plots
          #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
          tex3 = ROOT.TLatex(0.27,0.94,"Internal");
          tex3.SetNDC();
          tex3.SetTextFont(52);
          tex3.SetTextSize(0.0485);
          tex3.SetLineWidth(2);


          tex4 = ROOT.TLatex()
          if ("BefPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"Before pre-selection")
#            if ("BefPreS_Eta" in keyname2) :
#              print("BefPreS number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          elif ("N1" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After N-1 selection")
#            if ("N1_Eta" in keyname2) :
#              print("N-1 number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          elif ("PostPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After pre-selection")
          elif ("PostS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After selection")
#            if ("PostPreS_Eta" in keyname2) :
#              print("PostPreS number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          tex4.SetNDC();
          tex4.SetTextFont(52);
          tex4.SetTextSize(0.045);
          tex4.SetLineWidth(2);
          
          codeVersion = fileName[fileName.find("CodeV")+5:fileName.find("CodeV")+9]
          fileVersion = fileName[fileName.find("2018")+5:fileName.find("CodeV")+9]
          tex5 = ROOT.TLatex(0.07,0.03,fileVersion);
          tex5.SetNDC();
          tex5.SetTextFont(52);
          tex5.SetTextSize(0.0185);
          tex5.SetLineWidth(2);
          

          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if not (obj) : continue
          if obj.InheritsFrom("TObject"):
              can = obj
              obj.SetStats(0)
              can = ROOT.TCanvas(newname,newname,800,800)
              # Name of the png to be saved
              name = fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  ".png"
              if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
              if (obj.GetEntries() == 0 ) :
                continue
                
#             when I dont want to plot everything
              if (iDontWannaRunPlots) : continue
              if ("_region" in keyname2 or "CtrlPt_" in keyname2 or "Pred_" in keyname2 or "PDF" in keyname2 or "Hist_" in keyname2) : continue
#              if not ("Trigger" in keyname2 and obj.ClassName() == "TH3F") : continue
              if not ("GiTemplate" in keyname2) : continue

              if (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
                obj.SetTitle("")
                if ("VsProbQVsIas" in keyname2) :
                  can2 = ROOT.TCanvas(newname+"3",newname+"3",800,800)
                  can2.SetLogy()
                  
                  projA = obj.ProjectionX(newname+"_RegionA",obj.GetYaxis().FindBin(0.0),obj.GetYaxis().FindBin(0.1),obj.GetZaxis().FindBin(0.0),obj.GetZaxis().FindBin(0.1),"e")
                  projB = obj.ProjectionX(newname+"_RegionB",obj.GetYaxis().FindBin(0.0),obj.GetYaxis().FindBin(0.1),obj.GetZaxis().FindBin(0.1),obj.GetZaxis().FindBin(1.0),"e")
                  projC = obj.ProjectionX(newname+"_RegionC",obj.GetYaxis().FindBin(0.1),obj.GetYaxis().FindBin(1.0),obj.GetZaxis().FindBin(0.0),obj.GetZaxis().FindBin(0.1),"e")
                  projD = obj.ProjectionX(newname+"_RegionD",obj.GetYaxis().FindBin(0.1),obj.GetYaxis().FindBin(1.0),obj.GetZaxis().FindBin(0.1),obj.GetZaxis().FindBin(1.0),"e")
                  projA.SetMarkerColor(1)
                  projA.SetTitle("")
                  projB.SetMarkerColor(2)
                  projC.SetMarkerColor(3)
                  projD.SetMarkerColor(4)
                  projA.SetMarkerStyle(20)
                  projB.SetMarkerStyle(20)
                  projC.SetMarkerStyle(20)
                  projD.SetMarkerStyle(20)
                  projA.SetLineColor(1)
                  projB.SetLineColor(2)
                  projC.SetLineColor(3)
                  projD.SetLineColor(4)
                  
                  max1 = np.maximum(projA.GetMaximum(),projB.GetMaximum())
                  max2 = np.maximum(projC.GetMaximum(),projD.GetMaximum())
                  max = np.maximum(max1,max2)
                  
                  legProjXInRegions =  ROOT.TLegend(.4,.65,.80,.9,"","brNDC")
                  legProjXInRegions.SetTextFont(42)
                  legProjXInRegions.SetTextSize(0.035)
                  legProjXInRegions.SetBorderSize(1);
                  legProjXInRegions.SetLineColor(0);
                  legProjXInRegions.SetLineStyle(1);
                  legProjXInRegions.SetLineWidth(1);
                  legProjXInRegions.SetFillColor(0);
                  legProjXInRegions.SetFillStyle(1001);
                  
                  Err = double(0.0)
                  legProjXInRegions.AddEntry(projA,"Region A: " +str(int(projA.IntegralAndError(1,projA.GetNbinsX() + 1,Err,""))) + " #pm " + str(int(Err.value)),"LP")
                  legProjXInRegions.AddEntry(projB,"Region B: " +str(int(projB.IntegralAndError(1,projB.GetNbinsX() + 1,Err,""))) + " #pm " + str(int(Err.value)),"LP")
                  legProjXInRegions.AddEntry(projC,"Region C: " +str(int(projC.IntegralAndError(1,projC.GetNbinsX() + 1,Err,""))) + " #pm " + str(int(Err.value)),"LP")
                  legProjXInRegions.AddEntry(projD,"Region D: " +str(int(projD.IntegralAndError(1,projD.GetNbinsX() + 1,Err,""))) + " #pm " + str(int(Err.value)),"LP")
                  
                  projA.SetStats(0)
                  projA.SetMinimum(0.0001)
                  if (projA.Integral() == 0) : continue
                  if (projB.Integral() == 0) : continue
                  if (projC.Integral() == 0) : continue
                  if (projD.Integral() == 0) : continue
                  projA.Scale(1/projA.Integral())
                  newMax = projA.GetMaximum()
                  projB.Scale(1/projB.Integral())
                  projC.Scale(1/projC.Integral())
                  projD.Scale(1/projD.Integral())
                  projA.Draw("SAME")
                  projA.GetYaxis().SetRangeUser(0.00000001, newMax*10)
                  projA.GetYaxis().SetTitle("Norm tracks / bin")
                  projA.GetYaxis().SetTitleOffset(1.7)
                  projB.Draw("SAME")
                  projC.Draw("SAME")
                  projD.Draw("SAME")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  legProjXInRegions.Draw("SAME")

                  can2.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_ProjXInRegions.png")
                
                if ("ProbQVsProbXY" in keyname2) :
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
                if ("GiTemplate" in keyname2) :
                  projX = obj.ProjectionX()
                  projX.SetTitle("")
                  projX.SetStats(0)
                  projX.SetMarkerStyle(20)
                  projX.GetYaxis().SetTitle("Clusters")
                  projX.GetYaxis().SetTitleOffset(1.9)
                  projX.Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_ProjX.png")
                  
                  legGiCalib =  ROOT.TLegend(.6,.55,.80,.9,"","brNDC")
                  legGiCalib.SetTextFont(42)
                  legGiCalib.SetTextSize(0.03)
                  legGiCalib.SetBorderSize(1);
                  legGiCalib.SetLineColor(0);
                  legGiCalib.SetLineStyle(1);
                  legGiCalib.SetLineWidth(1);
                  legGiCalib.SetFillColor(0);
                  legGiCalib.SetFillStyle(1001);
                  
                  can3 = ROOT.TCanvas(newname+"3",newname+"3",800,800)
                  can3.SetLogy()
                  for x in range(1,obj.GetNbinsX()+1) :
                    projY = obj.ProjectionY(keyname2 +  "_ProjY"+str(x),x,x+1,1,obj.GetNbinsZ()+1,"e")
                    projY.SetTitle("")
                    projY.SetStats(0)
                    projY.SetMarkerStyle(20)
                    projY.SetLineColor(x)
                    projY.SetMarkerColor(x)
                    projY.GetYaxis().SetTitle("Clusters")
                    projY.GetYaxis().SetTitleOffset(1.4)
                    legGiCalib.AddEntry(projY, "Module index: " + str(x),"LP")
                    projY.Draw("SAME")
                  legGiCalib.Draw("SAME")
                  can3.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_ProjY.png")
                  can4 = ROOT.TCanvas(newname+"4",newname+"4",800,800)
                  can4.SetLogy()
                  for x in range(1,obj.GetNbinsX()+2) :
                    projZ = obj.ProjectionZ(keyname2 + "_ProjZ"+str(x),x,x+1,1,obj.GetNbinsY()+1,"e").Rebin(10)
                    projZ.SetTitle("")
                    projZ.SetStats(0)
                    projZ.SetMarkerStyle(20)
                    projZ.SetLineColor(x)
                    projZ.SetMarkerColor(x)
                    projZ.GetYaxis().SetTitle("Clusters")
                    projZ.GetYaxis().SetTitleOffset(1.4)
                    projZ.Draw("SAME")
                  legGiCalib.Draw("SAME")
                  can4.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_ProjZ.png")
                  continue
                else :
                  print("Following plot was skipped: "+str(keyname2))
                  continue
              if ((obj.ClassName() == "TH2F" or obj.ClassName() == "TH2D") and not (keyname2 == "GenPtVsRecoPt" or "PreS_" in keyname2 or "CutFlow" in keyname2 or "N1_" in keyname2 or "_p_" in keyname2 or "_pterr" in keyname2 )):
                obj.SetTitle("")
                obj.SetMarkerStyle(20)
                projOb = obj.ProjectionY(newname,bin,bin,"e")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                projOb.Draw("COLZ L")
                can.SaveAs(name)

              if ((obj.ClassName() == "TH2F") and "VsPixelLayer" in keyname2) :
                for i in range(4) :
                  obj.SetTitle("")
                  obj.SetMarkerStyle(20)
                  obj.SetStats(0)
                  obj.ProjectionX(newname+"_PixLayer"+str(i),obj.GetYaxis().FindBin(i),obj.GetYaxis().FindBin(i),"e").Draw()
                  obj.ProjectionX(newname+"_PixLayer"+str(i),obj.GetYaxis().FindBin(i),obj.GetYaxis().FindBin(i),"e").SetStats(0)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_PixLayer"+str(i+1)+".png")
              if ((obj.ClassName() == "TH2F") and "VsStripLayer" in keyname2) :
                obj.GetYaxis().SetBinLabel(1,"TIB L1")
                obj.GetYaxis().SetBinLabel(2,"TIB L2")
                obj.GetYaxis().SetBinLabel(3,"TIB L3")
                obj.GetYaxis().SetBinLabel(4,"TIB L4")
                obj.GetYaxis().SetBinLabel(5,"TOB L1")
                obj.GetYaxis().SetBinLabel(6,"TOB L2")
                obj.GetYaxis().SetBinLabel(7,"TOB L3")
                obj.GetYaxis().SetBinLabel(8,"TOB L4")
                obj.GetYaxis().SetBinLabel(9,"TOB L5")
                obj.GetYaxis().SetBinLabel(10,"TOB L6")
                obj.GetYaxis().SetBinLabel(11,"TID D1")
                obj.GetYaxis().SetBinLabel(12,"TID D2")
                obj.GetYaxis().SetBinLabel(13,"TID D3")
                obj.GetYaxis().SetBinLabel(14,"TEC D1")
                obj.GetYaxis().SetBinLabel(15,"TEC D2")
                obj.GetYaxis().SetBinLabel(16,"TEC D3")
                obj.GetYaxis().SetBinLabel(17,"TEC D4")
                obj.GetYaxis().SetBinLabel(18,"TEC D5")
                obj.GetYaxis().SetBinLabel(19,"TEC D6")
                obj.GetYaxis().SetBinLabel(20,"TEC D7")
                obj.GetYaxis().SetBinLabel(21,"TEC D8")
                obj.GetYaxis().SetBinLabel(22,"TEC D9")
                obj.GetYaxis().SetBinLabel(23,"TEC D10")
                for i in range(20) :
                  obj.SetTitle("")
                  obj.SetMarkerStyle(20)
                  obj.SetStats(0)
                  obj.ProjectionX(newname+"_StripLayer"+str(i),obj.GetYaxis().FindBin(i),obj.GetYaxis().FindBin(i),"e").Draw()
                  obj.ProjectionX(newname+"_StripLayer"+str(i),obj.GetYaxis().FindBin(i),obj.GetYaxis().FindBin(i),"e").SetStats(0)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_StripLayer"+str(i+1)+".png")
              if ((obj.ClassName() == "TH2F") and ("Clu" in keyname2)) :
                profYobj = obj.ProfileY()
#                profYobj.GetYaxis().SetTitle(axisXTitle)
                profYobj.GetYaxis().SetTitleOffset(1.5)
                profYobj.GetYaxis().SetLabelSize(0.03)
                profYobj.SetStats(0)
                if ("CluNormChargeVsStripLayer" in keyname2) :
                  profYobj.GetXaxis().SetBinLabel(1,"TIB L1")
                  profYobj.GetXaxis().SetBinLabel(2,"TIB L2")
                  profYobj.GetXaxis().SetBinLabel(3,"TIB L3")
                  profYobj.GetXaxis().SetBinLabel(4,"TIB L4")
                  profYobj.GetXaxis().SetBinLabel(5,"TOB L1")
                  profYobj.GetXaxis().SetBinLabel(6,"TOB L2")
                  profYobj.GetXaxis().SetBinLabel(7,"TOB L3")
                  profYobj.GetXaxis().SetBinLabel(8,"TOB L4")
                  profYobj.GetXaxis().SetBinLabel(9,"TOB L5")
                  profYobj.GetXaxis().SetBinLabel(10,"TOB L6")
                  profYobj.GetXaxis().SetBinLabel(11,"TID D1")
                  profYobj.GetXaxis().SetBinLabel(12,"TID D2")
                  profYobj.GetXaxis().SetBinLabel(13,"TID D3")
                  profYobj.GetXaxis().SetBinLabel(14,"TEC D1")
                  profYobj.GetXaxis().SetBinLabel(15,"TEC D2")
                  profYobj.GetXaxis().SetBinLabel(16,"TEC D3")
                  profYobj.GetXaxis().SetBinLabel(17,"TEC D4")
                  profYobj.GetXaxis().SetBinLabel(18,"TEC D5")
                  profYobj.GetXaxis().SetBinLabel(19,"TEC D6")
                  profYobj.GetXaxis().SetBinLabel(20,"TEC D7")
                  profYobj.GetXaxis().SetBinLabel(21,"TEC D8")
                  profYobj.GetXaxis().SetBinLabel(22,"TEC D9")
                  profYobj.GetXaxis().SetBinLabel(23,"TEC D10")
                profYobj.DrawClone("COLZ")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_profileY.png")
                
                if ("CluSpecInCPE" in keyname2) :
                  obj.SetStats(0)
                  obj.GetXaxis().SetTitle("")
                  obj.GetXaxis().SetBinLabel(1,"isOnEdge")
                  obj.GetXaxis().SetBinLabel(2,"hasBadPixels")
                  obj.GetXaxis().SetBinLabel(3,"spansTwoROCs")
                  obj.GetXaxis().SetBinLabel(4,"AllClusters")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  obj.Draw("COLZ L")
#                  obj.GetXaxis().SetTitle(axisXTitle)
                obj.DrawClone("COLZ L")
#                obj.GetYaxis().SetTitle(axisYTitle)
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(name)
              elif ("GenPtVsRecoPt" in keyname2) :
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ")
                print(str(keyname2) + ": " + str(round(obj.GetCorrelationFactor(),2)))
              elif ("Stab" in keyname2) :
                profXobj = ROOT.TH1F(keyname2+"profX" , keyname2+"profX" , 545, 271000,325500)
                profXobj = obj.ProfileX()
#                profYobj.GetYaxis().SetTitle(axisXTitle)
#                profYobj.GetYaxis().SetTitleOffset(1.5)
#                profYobj.GetYaxis().SetLabelSize(0.03)
                profXobj.SetStats(0)
                profXobj.DrawClone("COLZ")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_profileX.png")
              elif ("IasVs" in keyname2 and not obj.ClassName() == "TH3F") :
                obj.SetMarkerStyle(20)
                if ("Angle" in keyname2 or "NumSibling" in keyname2) :
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  obj.ProjectionY(newname,obj.GetXaxis().FindBin(0.7),obj.GetNbinsX()+1,"e").Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
                  projObjs = obj.ProjectionY(newname,1,obj.GetXaxis().FindBin(0.7),"e")
                  if (projObjs.GetEntries()==0) : continue
                  projObjs.Draw("COLZ")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
                  obj.ProjectionY(newname,1,obj.GetNbinsX()+1,"e").Draw("COLZ")
                else :
                  projObject = obj.ProjectionY(newname+"_lowIas",1,obj.GetXaxis().FindBin(0.7),"e")
                  if (projObject.GetEntries()==0) : continue
                  myPie = ROOT.TPie(projObject)
                  myPie.SetLabelFormat("%txt (%perc)")
                  myPie.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
                  
                  objProj = obj.ProjectionY(newname+"_highIas",obj.GetXaxis().FindBin(0.7),obj.GetNbinsX()+1,"e")
                  if (objProj.GetEntries()==0) : continue
                  myPie2 = ROOT.TPie(objProj)
                  myPie2.SetLabelFormat("%txt (%perc)")
                  myPie2.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie2.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
                  
                  myPie3 = ROOT.TPie(obj.ProjectionY(newname,1,obj.GetNbinsX()+1,"e"))
                  myPie3.SetLabelFormat("%txt (%perc)")
                  myPie3.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie3.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  ".png")
              elif ("Trigger" in keyname2 and obj.ClassName() == "TH3F") :
                profZobj = obj.Project3DProfile("zy")
                profZobj.SetStats(0)
                profZobj.GetZaxis().SetTitle("Efficiency")
                profZobj.SetTitle("")
                profZobj.GetXaxis().SetTitle("MET (GeV)")
                profZobj.GetYaxis().SetTitle("H_{T} (GeV)")
                profZobj.GetYaxis().SetTitleOffset(1.7)
                profZobj.Rebin2D(2)
                profZobj.DrawClone("COLZ")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_profileZ.png")
              elif ("EoPVs" in keyname2) :
                obj.SetMarkerStyle(20)
                if ("Angle" in keyname2) :
                  obj.ProjectionY(newname+"_lowEoP",1,obj.GetXaxis().FindBin(0.85),"e").Draw("COLZ")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowEoP.png")
#                  obj.ProjectionY(newname+"_highEoP",obj.GetXaxis().FindBin(0.85),obj.GetNbinsX()+1,"e").Draw("COLZ")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highEoP.png")
                elif ("EoPVsPfType" in keyname2) :
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
                  projObj = obj.ProjectionY(newname+"_lowEoP",1,obj.GetXaxis().FindBin(0.85),"e")
                  if (projObj.GetEntries() == 0) : continue
                  myPie = ROOT.TPie(projObj)
                  myPie.SetLabelFormat("%txt (%perc)")
                  myPie.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowEoP.png")
                  
                  projObj2 = obj.ProjectionY(newname+"_highEoP",obj.GetXaxis().FindBin(0.85),obj.GetNbinsX()+1,"e")
                  if (projObj2.GetEntries() == 0) : continue
                  myPie2 = ROOT.TPie(projObj2)
                  myPie2.SetLabelFormat("%txt (%perc)")
                  myPie2.SetLabelsOffset(-.27)
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  myPie2.Draw("R<")
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highEoP.png")
              elif (obj.ClassName() == "TH2F" and  (("ProbQVsIas" in keyname2) or ( "ProbQNoL1VsIas" in keyname2))) :
#                obj.GetXaxis().SetTitle(axisXTitle)
#                obj.GetYaxis().SetTitle(axisYTitle)
                obj.GetYaxis().SetTitleOffset(1.3)
                obj.GetYaxis().SetLabelSize(0.03)
                obj.GetYaxis().SetTitle("G_{i}^{Strips}")
                obj.GetXaxis().SetTitle("F_{i}^{Pixels}")
                obj.SetStats(0)
                obj.DrawClone("COLZ")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                print(str(keyname2) + ": " + str(round(obj.GetCorrelationFactor(),5)))
                can.SaveAs(name)
                if ("ProbQNoL1VsIas" in keyname2) :
                  can2 = ROOT.TCanvas(newname+"2",newname+"2",800,800)
                  can2.SetLogy()
                  legIasForProbQSlice =  ROOT.TLegend(.30,.80,.80,.90,"","brNDC")
                  projY1 = obj.ProjectionY("IasForProbQSlice_Slice1",obj.GetXaxis().FindBin(0.9),obj.GetXaxis().FindBin(1.0),"e")
                  projY1.SetStats(0)
                  projY1.SetMarkerColor(1)
                  projY1.SetLineColor(1)
                  projY1.SetMarkerStyle(20)
                  projY1.Draw("SAME")
                  projY1.GetYaxis().SetTitle("Normalized Tracks / 0.05")
                  projY1.GetYaxis().SetTitleOffset(1.5)
                  
                  projY2 = obj.ProjectionY("IasForProbQSlice_Slice2",obj.GetXaxis().FindBin(0.3),obj.GetXaxis().FindBin(0.9),"e")
                  projY2.SetMarkerStyle(20)
                  projY2.SetMarkerColor(2)
                  projY2.SetLineColor(2)
                  projY2.Draw("SAME")
                  
                  numTracks1 = projY1.Integral()
                  numTracks2 = projY2.Integral()
                  
                  if (numTracks1>0) : projY1.Scale(1/projY1.Integral())
                  if (numTracks2>0) : projY2.Scale(1/projY2.Integral())
                  
                  blind = True
#                  blind = False

                  if (blind and not "CR" in keyname2 and isData) :
                    projY1.SetBinContent(7,0)
                    projY1.SetBinContent(8,0)
                    projY1.SetBinContent(9,0)
                    projY1.SetBinContent(10,0)
                  projY1.SetMaximum(projY1.GetMaximum()*100)
                  
                  rp2 = ROOT.TRatioPlot(projY1,projY2,"divsym") #, "diffsigerrasym"

                  rp2.SetH1DrawOpt("P");
                  rp2.SetH2DrawOpt("P");

                  rp2.Draw()
                  
                  rp2.SetLeftMargin(0.13);
                  rp2.SetRightMargin(0.05);
                  rp2.SetUpTopMargin(0.1);
                  rp2.SetLowTopMargin(0.02);
                  rp2.SetLowBottomMargin(0.35);


                  rp2.GetLowerRefGraph().SetMinimum(0)
                  rp2.GetLowerRefGraph().SetMaximum(3.5);
                  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
                  #rp.GetLowerRefGraph().SetLineColor(0) #0
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
                  legIasForProbQSlice.SetTextFont(42)
                  legIasForProbQSlice.SetTextSize(0.035)
                  legIasForProbQSlice.SetBorderSize(1);
                  legIasForProbQSlice.SetLineColor(0);
                  legIasForProbQSlice.SetLineStyle(1);
                  legIasForProbQSlice.SetLineWidth(1);
                  legIasForProbQSlice.SetFillColor(0);
                  legIasForProbQSlice.SetFillStyle(1001);
                  legIasForProbQSlice.AddEntry(projY1,"F_{i}^{Pixels} (0.9-1.0), #Tracks: " +str(round(numTracks1)),"LP")
                  legIasForProbQSlice.AddEntry(projY2,"F_{i}^{Pixels}  (0.3-0.9), #Tracks: " +str(round(numTracks2)),"LP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  legIasForProbQSlice.Draw("SAME")
                  can2.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 + "_Norm_ProbQSlice.png")
                  
                  can3 = ROOT.TCanvas(newname+"2",newname+"2",800,800)
                  can3.SetLogy()
                  legIasForProbQSlice =  ROOT.TLegend(.30,.80,.80,.90,"","brNDC")
                  projY1 = obj.ProjectionX("IasForProbQSlice_Slice1",1,obj.GetYaxis().FindBin(0.1),"e")
                  projY1.SetStats(0)
                  projY1.SetMarkerColor(1)
                  projY1.SetLineColor(1)
                  projY1.SetMarkerStyle(20)
                  projY1.Draw("SAME")
                  projY1.GetYaxis().SetTitle("Normalized Tracks / 0.05")
                  projY1.GetYaxis().SetTitleOffset(1.5)

                  projY2 = obj.ProjectionX("IasForProbQSlice_Slice2",obj.GetYaxis().FindBin(0.1),obj.GetYaxis().FindBin(1.0),"e")
                  projY2.SetMarkerStyle(20)
                  projY2.SetMarkerColor(2)
                  projY2.SetLineColor(2)
                  projY2.Draw("SAME")
                  
                  numTracks3 = projY1.Integral()
                  numTracks4 = projY2.Integral()
                  if (numTracks3>0) : projY1.Scale(1/projY1.Integral())
                  if (numTracks4>0) : projY2.Scale(1/projY2.Integral())
                
                  if (blind and not "CR" in keyname2 and isData) :
                    projY1.SetBinContent(19,0)
                    projY1.SetBinContent(20,0)
                  
                  projY1.SetMaximum(projY1.GetMaximum()*100)
                  
                  rp3 = ROOT.TRatioPlot(projY1,projY2,"divsym") #, "diffsigerrasym"

                  rp3.SetH1DrawOpt("P");
                  rp3.SetH2DrawOpt("P");

                  rp3.Draw()
                  
                  rp3.SetLeftMargin(0.13);
                  rp3.SetRightMargin(0.05);
                  rp3.SetUpTopMargin(0.1);
                  rp3.SetLowTopMargin(0.02);
                  rp3.SetLowBottomMargin(0.35);


                  rp3.GetLowerRefGraph().SetMinimum(0)
                  rp3.GetLowerRefGraph().SetMaximum(3.5)
                  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
                  #rp.GetLowerRefGraph().SetLineColor(0) #0
                  rp3.GetLowerRefGraph().SetMarkerStyle(20)
                  rp3.GetLowerRefGraph().SetMarkerSize(1);
                  rp3.GetLowYaxis().SetNdivisions(505);
                  rp3.GetLowerRefYaxis().SetTitle("Ratio");
                  rp3.GetLowerRefYaxis().SetTitleSize(0.05);
                  rp3.GetLowerRefYaxis().SetTitleOffset(1);
                  rp3.GetLowerRefYaxis().SetLabelSize(0.035);


                  rp3.GetLowerRefXaxis().SetTitleSize(0.05);
                  rp3.GetLowerRefXaxis().SetTitleOffset(0.8);
                  rp3.GetLowerRefXaxis().SetLabelSize(0.035);
                  legIasForProbQSlice.SetTextFont(42)
                  legIasForProbQSlice.SetTextSize(0.035)
                  legIasForProbQSlice.SetBorderSize(1);
                  legIasForProbQSlice.SetLineColor(0);
                  legIasForProbQSlice.SetLineStyle(1);
                  legIasForProbQSlice.SetLineWidth(1);
                  legIasForProbQSlice.SetFillColor(0);
                  legIasForProbQSlice.SetFillStyle(1001);
                  legIasForProbQSlice.AddEntry(projY1,"G_{i}^{Strips} (0.0-0.1), #Tracks: "+str(round(numTracks3)),"LP")
                  legIasForProbQSlice.AddEntry(projY2,"G_{i}^{Strips} (0.1-1.0), #Tracks: "+str(round(numTracks4)),"LP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  legIasForProbQSlice.Draw("SAME")
                  can3.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 + "_Norm_IasSlice.png")

                if ("ProbQNoL1VsIas" in keyname2) :
                  can2 = ROOT.TCanvas(newname+"2",newname+"2",800,800)
                  can2.SetLogy()
                  legIasForProbQSlice =  ROOT.TLegend(.30,.80,.80,.90,"","brNDC")
                  projY1 = obj.ProjectionY("IasForProbQSlice_Slice1",obj.GetXaxis().FindBin(0.9),obj.GetXaxis().FindBin(1.0),"e")
                  projY1.SetStats(0)
                  projY1.SetMarkerColor(1)
                  projY1.SetLineColor(1)
                  projY1.SetMarkerStyle(20)
                  projY1.Draw("SAME")
                  projY1.GetYaxis().SetTitle("Tracks / 0.05")
                  projY1.GetYaxis().SetTitleOffset(1.5)
                  
                  projY2 = obj.ProjectionY("IasForProbQSlice_Slice2",obj.GetXaxis().FindBin(0.3),obj.GetXaxis().FindBin(0.9),"e")
                  projY2.SetMarkerStyle(20)
                  projY2.SetMarkerColor(2)
                  projY2.SetLineColor(2)
                  projY2.Draw("SAME")
                  
                  numTracks1 = projY1.Integral()
                  numTracks2 = projY2.Integral()

                  if (blind and not "CR" in keyname2 and isData) :
                    projY1.SetBinContent(7,0)
                    projY1.SetBinContent(8,0)
                    projY1.SetBinContent(9,0)
                    projY1.SetBinContent(10,0)
                  projY1.SetMaximum(projY2.GetMaximum()*100)
                  projY1.SetMinimum(0.1)
                  
                  rp2 = ROOT.TRatioPlot(projY1,projY2,"divsym") #, "diffsigerrasym"

                  rp2.SetH1DrawOpt("P");
                  rp2.SetH2DrawOpt("P");

                  rp2.Draw()
                  
                  rp2.SetLeftMargin(0.13);
                  rp2.SetRightMargin(0.05);
                  rp2.SetUpTopMargin(0.1);
                  rp2.SetLowTopMargin(0.02);
                  rp2.SetLowBottomMargin(0.35);


                  rp2.GetLowerRefGraph().SetMinimum(0)
                  rp2.GetLowerRefGraph().SetMaximum(3.5);
                  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
                  #rp.GetLowerRefGraph().SetLineColor(0) #0
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
                  legIasForProbQSlice.SetTextFont(42)
                  legIasForProbQSlice.SetTextSize(0.035)
                  legIasForProbQSlice.SetBorderSize(1);
                  legIasForProbQSlice.SetLineColor(0);
                  legIasForProbQSlice.SetLineStyle(1);
                  legIasForProbQSlice.SetLineWidth(1);
                  legIasForProbQSlice.SetFillColor(0);
                  legIasForProbQSlice.SetFillStyle(1001);
                  legIasForProbQSlice.AddEntry(projY1,"F_{i}^{Pixels} (0.9-1.0), #Tracks: " +str(round(numTracks1)),"LP")
                  legIasForProbQSlice.AddEntry(projY2,"F_{i}^{Pixels}  (0.3-0.9), #Tracks: " +str(round(numTracks2)),"LP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  legIasForProbQSlice.Draw("SAME")
                  can2.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 + "_NotNorm_ProbQSlice.png")
                  
                  can3 = ROOT.TCanvas(newname+"2",newname+"2",800,800)
                  can3.SetLogy()
                  legIasForProbQSlice =  ROOT.TLegend(.30,.80,.80,.90,"","brNDC")
                  projY1 = obj.ProjectionX("IasForProbQSlice_Slice1",1,obj.GetYaxis().FindBin(0.1),"e")
                  projY1.SetStats(0)
                  projY1.SetMarkerColor(1)
                  projY1.SetLineColor(1)
                  projY1.SetMarkerStyle(20)
                  projY1.Draw("SAME")
                  projY1.GetYaxis().SetTitle("Normalized Tracks / 0.05")
                  projY1.GetYaxis().SetTitleOffset(1.5)

                  projY2 = obj.ProjectionX("IasForProbQSlice_Slice2",obj.GetYaxis().FindBin(0.1),obj.GetYaxis().FindBin(1.0),"e")
                  projY2.SetMarkerStyle(20)
                  projY2.SetMarkerColor(2)
                  projY2.SetLineColor(2)
                  projY2.Draw("SAME")
                  
                  numTracks3 = projY1.Integral()
                  numTracks4 = projY2.Integral()
                
                  if (blind and not "CR" in keyname2 and isData) :
                    projY1.SetBinContent(19,0)
                    projY1.SetBinContent(20,0)
                  
                  projY1.SetMaximum(projY1.GetMaximum()*100)
                  
                  rp3 = ROOT.TRatioPlot(projY1,projY2,"divsym") #, "diffsigerrasym"

                  rp3.SetH1DrawOpt("P");
                  rp3.SetH2DrawOpt("P");

                  rp3.Draw()
                  
                  rp3.SetLeftMargin(0.13);
                  rp3.SetRightMargin(0.05);
                  rp3.SetUpTopMargin(0.1);
                  rp3.SetLowTopMargin(0.02);
                  rp3.SetLowBottomMargin(0.35);


                  rp3.GetLowerRefGraph().SetMinimum(0)
                  rp3.GetLowerRefGraph().SetMaximum(3.5)
                  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
                  #rp.GetLowerRefGraph().SetLineColor(0) #0
                  rp3.GetLowerRefGraph().SetMarkerStyle(20)
                  rp3.GetLowerRefGraph().SetMarkerSize(1);
                  rp3.GetLowYaxis().SetNdivisions(505);
                  rp3.GetLowerRefYaxis().SetTitle("Ratio");
                  rp3.GetLowerRefYaxis().SetTitleSize(0.05);
                  rp3.GetLowerRefYaxis().SetTitleOffset(1);
                  rp3.GetLowerRefYaxis().SetLabelSize(0.035);


                  rp3.GetLowerRefXaxis().SetTitleSize(0.05);
                  rp3.GetLowerRefXaxis().SetTitleOffset(0.8);
                  rp3.GetLowerRefXaxis().SetLabelSize(0.035);
                  legIasForProbQSlice.SetTextFont(42)
                  legIasForProbQSlice.SetTextSize(0.035)
                  legIasForProbQSlice.SetBorderSize(1);
                  legIasForProbQSlice.SetLineColor(0);
                  legIasForProbQSlice.SetLineStyle(1);
                  legIasForProbQSlice.SetLineWidth(1);
                  legIasForProbQSlice.SetFillColor(0);
                  legIasForProbQSlice.SetFillStyle(1001);
                  legIasForProbQSlice.AddEntry(projY1,"G_{i}^{Strips} (0.0-0.1), #Tracks: "+str(round(numTracks3)),"LP")
                  legIasForProbQSlice.AddEntry(projY2,"G_{i}^{Strips} (0.1-1.0), #Tracks: "+str(round(numTracks4)),"LP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  legIasForProbQSlice.Draw("SAME")
                  can3.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 + "_NotNorm_IasSlice.png")

              elif (obj.ClassName() == "TH2F" and  (("Trigger" in keyname2) and ( "Vs" in keyname2))) :
                profYobj = obj.ProfileY()
                profYobj.SetStats(0)
                profYobj.GetYaxis().SetTitle("Efficiency")
                profYobj.DrawClone("COLZ")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_profileY.png")
              
              if (keyname2== "CutFlow" or keyname2== "EventCutFlow") :
                if (blind and isData) :
                  obj.SetBinContent(18,0)
                  obj.SetBinContent(19,0)
                  obj.SetBinContent(20,0)
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.GetXaxis().SetTitle("")
                obj.GetXaxis().SetTitle("")
                ROOT.gStyle.SetPaintTextFormat(".2g");
                obj.Draw("HISTOTEXT00")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
#                obj.GetYaxis().SetRangeUser(0.,1.3)
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_withNumbersNotNorm.png")
              if (keyname2== "CutFlow" or keyname2== "EventCutFlow") :
                if (blind and isData) :
                  obj.SetBinContent(18,0)
                  obj.SetBinContent(19,0)
                  obj.SetBinContent(20,0)
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetXaxis().SetTitle("")
                if (keyname2== "CutFlow") :  obj.GetYaxis().SetTitle("Track efficiency")
                if (keyname2== "EventCutFlow") : obj.GetYaxis().SetTitle("Event efficiency")
                obj.GetXaxis().SetTitle("")
                ROOT.gStyle.SetPaintTextFormat(".2g");
                obj.Draw("HISTOTEXT00")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.GetYaxis().SetRangeUser(0.,1.3)
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_withNumbers.png")
              elif (keyname2== "CutFlowReverse") :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
#                obj.GetXaxis().SetBinLabel(1,"Trigger")
#                obj.GetXaxis().SetBinLabel(2,"pT")
#                obj.GetXaxis().SetBinLabel(3,"Eta")
#                obj.GetXaxis().SetBinLabel(4,"NumPixHits")
#                obj.GetXaxis().SetBinLabel(5,"ValidFract")
#                obj.GetXaxis().SetBinLabel(6,"NumDeDx")
#                obj.GetXaxis().SetBinLabel(7,"HighPurity")
#                obj.GetXaxis().SetBinLabel(8,"Chi2oDOF")
#                obj.GetXaxis().SetBinLabel(9,"dz")
#                obj.GetXaxis().SetBinLabel(10,"dxy")
#                obj.GetXaxis().SetBinLabel(11,"EoP")
#                obj.GetXaxis().SetBinLabel(12,"")
#                obj.GetXaxis().SetBinLabel(13,"dRminCaloJet")
#                obj.GetXaxis().SetBinLabel(14,"MiniIso")
#                obj.GetXaxis().SetBinLabel(15,"PFid")
#                obj.GetXaxis().SetBinLabel(16,"Ih")
#                obj.GetXaxis().SetBinLabel(17,"")
                obj.GetXaxis().SetTitle("")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ L")
              if (keyname2 == "ErrorHisto") :
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
#                obj.GetXaxis().SetBinLabel(1,"All tracks")
#                obj.GetXaxis().SetBinLabel(2,"No track/glob muon")
#                obj.GetXaxis().SetBinLabel(3,"Track is null")
#                obj.GetXaxis().SetBinLabel(4,"No PV")
#                obj.GetXaxis().SetBinLabel(5,"No gen match")
#                obj.GetXaxis().SetBinLabel(6,"Gen match too far")
#                obj.GetXaxis().SetBinLabel(7,"No dEdx")
#                obj.GetXaxis().SetBinLabel(8,"Cosmic track")
#                obj.GetXaxis().SetBinLabel(9,"Has91 status")
                obj.GetXaxis().SetTitle("")
                obj.GetYaxis().SetTitle("")
                obj.SetMaximum(1.4)
                obj.GetXaxis().SetTitle("")
                obj.SetTitle("")
                obj.Draw("COLZ L")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(name)
              if ("TriggerType" in keyname2) :
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.GetYaxis().SetTitle("Events / category")
                ROOT.gStyle.SetPaintTextFormat(".2g");
                obj.Draw("SAMEHISTOTEXT00")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(name)
              if ("TriggerType" in keyname2) :
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.Integral(1,3))
                obj.GetYaxis().SetTitle("Norm events / category")
                ROOT.gStyle.SetPaintTextFormat(".2g");
                obj.Draw("SAMEHISTOTEXT00")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_Normalized.png")
              elif ("_pfType" in keyname2) :
                obj.SetMarkerColor(1)
                obj.SetLineColor(1)
                obj.SetMarkerStyle(20)
                obj.SetStats(0)
                obj.Scale(1/obj.GetMaximum())
                obj.GetYaxis().SetTitle("Efficiency")
#                obj.GetXaxis().SetBinLabel(1,"AllTracks")
#                obj.GetXaxis().SetBinLabel(2,"PFtracks")
#                obj.GetXaxis().SetBinLabel(3,"isElectron")
#                obj.GetXaxis().SetBinLabel(4,"isMuon")
#                obj.GetXaxis().SetBinLabel(5,"isPhoton")
#                obj.GetXaxis().SetBinLabel(6,"isChHadron")
#                obj.GetXaxis().SetBinLabel(7,"isNeutHadron")
#                obj.GetXaxis().SetBinLabel(8,"isUndefined")
#                obj.GetXaxis().SetBinLabel(9,"notPFtrack")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ L")

              elif ((keyname2 == "CutFlowEta") or (keyname2 == "CutFlowProbQ") or (keyname2 == "CutFlowPfType") or (keyname2 == "CutFlowProbQ")) :
                obj.SetStats(0)
                obj.GetYaxis().SetTitle("Efficiency")
                obj.GetYaxis().SetTitle("")
                obj.GetYaxis().SetTitle("")
                if (keyname2 == "CutFlowPfType"):
                  obj.Scale(1/obj.GetMaximum())
                  for x in range(1,obj.GetNbinsX()+1) :
                    localMax = obj.GetBinContent(x,1)
                    for y in range(1,obj.GetNbinsY()+1) :
                      value = obj.GetBinContent(x,y)
                      if (value <= 0 or localMax <= 0) : continue
                      obj.SetBinContent(x,y, value/localMax)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                obj.Draw("COLZ")
                can.SaveAs(name)
              #and ("Mass" not in keyname2)
              #and ("PostS" not in keyname2)
              elif ("Max" not in keyname2) and ("Pred" not in keyname2) and ("PDF" not in keyname2):
                obj.SetMarkerStyle(20)
                obj.SetTitle("")
#                obj.GetXaxis().SetTitle(axisXTitle)
#                obj.GetYaxis().SetTitle(axisYTitle)
                obj.Draw("COLZ")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(name)
              
              if ("Angle" in keyname2 and obj.ClassName() == "TH2F" ) :
                obj.Draw("COLZ")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_2d.png")
              
              ## ----------------------------------------------------------------------------------------------
              # now let's plot everything in a logy scale
              obj.SetMarkerStyle(20)
              obj.SetMinimum(0.000001)
#              obj.SetMaximum(10000)
              if ("GenID" in keyname2 or "GenEnviromentID" in keyname2) :
                continue
                
              if (keyname2=="CutFlowProbQ" or keyname2=="CutFlowPfType" or keyname2=="CutFlowEta" or "Vs" in keyname2) :
#              or "PostPreS_IasPixelIhVsLayer" in keyname2 or "PostPreS_IasStripIhVsLayer" in keyname2
                obj.SetMaximum(obj.GetMaximum())
                obj.SetMinimum(0.000001)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                can.SetLogz()

              else :
                obj.SetMaximum(obj.GetMaximum()*100)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
#                obj.SetMinimum(0.0001)
                can.SetLogy()

                obj.SetTitle("")
              can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_logy.png")
              #can.SaveAs(name.replace(".png",".pdf"))
              #can.SaveAs(name.replace(".png",".C"))
              can.Close()

          else:
              print(keyname+"   "+newname + " does not inherit from TObject" )

##template_2018A_v2.root
#name = "2017A"
#version = "v4"
#fOut = ROOT.TFile.Open('template_{}_{}.root'.format(name,version),'UPDATE')
#fIn = ROOT.TFile.Open("Histos_numEvent2000.root")
#templ1 = fIn.Get('HSCParticleAnalyzer/BaseName/Calibration_GiTemplate_PU_1')
#
#templ1.Write()
#fOut.Close()
