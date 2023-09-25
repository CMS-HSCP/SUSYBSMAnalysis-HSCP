import ROOT, sys, os, time, re, array
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
BinNumber = sys.argv[2] if (len(sys.argv)==1) else 3

bin = int(BinNumber)
# bin 3: pt>60 and I_as > 0.05
# bin 25: pt>65 and I_as > 0.175
# bin 28: pt>65 and I_as > 0.3

#blind = True
blind = False

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

class Vavilov_Func:
    def __init__(self):
        self.pdf = ROOT.Math.VavilovAccurate()
    
    def __call__(self, x, p):
        kappa = p[0]
        beta2 = p[1]
        return p[4] * (self.pdf.Pdf((x[0] - p[2]) / p[3], kappa, beta2))

def GaussWithLandau(x, par):
    landau = ROOT.TMath.Landau(x[0], par[0], par[1])
    gauss = ROOT.TMath.Gaus(x[0], 0, par[2], 1)
    dx = x[1] - x[0]
    conv = np.convolve(landau, gauss, mode='same') * dx
    return par[3] * conv
    
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
          if ("__" in keyname2) : continue
          # The plot should be TCanvas
          newname = dirname + "/" + keyname+ "/" + keyname2
#          print("newname: "+newname)
          obj = f.Get(newname)
          
          
          tex2 = ROOT.TLatex(0.15,0.94,"CMS");
          #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

#          tex3 = ROOT.TLatex(0.27,0.94,"Simulation"); # for square plots
#          tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
          tex3 = ROOT.TLatex(0.29,0.94,"Internal");
#          tex3 = ROOT.TLatex(0.27,0.94,"");
          tex3.SetNDC();
          tex3.SetTextFont(52);
          tex3.SetTextSize(0.0485);
          tex3.SetLineWidth(2);


          tex4 = ROOT.TLatex()
          if ("BefPreS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.94,"Before preselection")
#            if ("BefPreS_Eta" in keyname2) :
#              print("BefPreS number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          elif ("N1" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.94,"After N-1 selection")
#            if ("N1_Eta" in keyname2) :
#              print("N-1 number of tracks in plot (" +keyname2 + ") : "+str(obj.Integral()))
          elif ("PostPreS" in keyname2 or "Stab" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.94,"After preselection")
          elif ("PostS" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.94,"After selection")
          elif ("K_and_C" in keyname2) :
            tex4 = ROOT.TLatex(0.5,0.94,"After calibration selection")

          tex4.SetNDC();
          tex4.SetTextFont(52);
          tex4.SetTextSize(0.045);
          tex4.SetLineWidth(2);
          
          codeVersion = fileName[fileName.find("CodeV")+5:fileName.find("CodeV")+9]
          beginFileVersion = fileName.find("2018")+5 if (fileName.find("2018") > 0) else fileName.find("SingleMuon")
          fileVersion = fileName[beginFileVersion:fileName.find("CodeV")+9]
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
#              if not ("RunNumVsPixCluChargeAfter" in keyname2) : continue
#              if not ("PostS_SR2PASS" in keyname2) : continue

              if ("PostS_SR2PASS_RunVsLs" in keyname2) and not blind :
                  for ix in range(1, obj.GetNbinsX() + 1):
                    for iy in range(1, obj.GetNbinsY() + 1):
                        content = obj.GetBinContent(ix, iy)
                        if content != 0.0:
                            x = obj.GetXaxis().GetBinCenter(ix)
                            y = obj.GetYaxis().GetBinCenter(iy)
                            print(f"Non-zero yield at ({x}, {y}): {content}")

#              if not ((obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D") and "VsProbQVsIas" in keyname2) : continue
              if ("Gen" in keyname2 and isData) : continue
#                 print(obj.ClassName())
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
                if ("ProbQNoL1VsIasVsPt" in keyname2) :
                  canLog = ROOT.TCanvas("newname"+keyname2,"newname"+keyname2,800,800)
                  canLog.SetLogy()
                  legProjXInRegions =  ROOT.TLegend(.5,.75,.80,.9,"","brNDC")
                  legProjXInRegions.SetHeader("FAIL region (F_{i}^{pixels} < 0.9)","C")
                  legProjXInRegions.SetTextFont(42)
                  legProjXInRegions.SetTextSize(0.035)
                  legProjXInRegions.SetBorderSize(1);
                  legProjXInRegions.SetLineColor(0);
                  legProjXInRegions.SetLineStyle(1);
                  legProjXInRegions.SetLineWidth(1);
                  legProjXInRegions.SetFillColor(0);
                  legProjXInRegions.SetFillStyle(1001);
                  obj.GetXaxis().SetRange(1,obj.GetXaxis().FindBin(0.9)-1)
                 
                  projFail = obj.Project3D("YZ")
                  projFail.SetTitle("")
                  projFail.GetXaxis().SetRangeUser(200,4000)
                  projFailPt1 = projFail.ProjectionY("name1").Rebin(2)
                  projFailPt1.SetStats(0)
                  projFailPt1.SetMarkerStyle(20)
                  projFailPt1.SetLineColor(1)
                  projFailPt1.SetMarkerColor(1)
                  projFailPt1.SetMaximum(projFailPt1.GetMaximum()*100)
                  projFailPt1.DrawClone("SAMEP")
                  legProjXInRegions.AddEntry(projFailPt1, "p_{T} = 200 - inf GeV","LP")

                  projFail.GetXaxis().SetRangeUser(300,4000)
                  projFailPt2 = projFail.ProjectionY("name2").Rebin(2)
                  projFailPt2.SetStats(0)
                  projFailPt2.SetMarkerStyle(20)
                  projFailPt2.SetLineColor(2)
                  projFailPt2.SetMarkerColor(2)
                  projFailPt2.DrawClone("SAMEP")
                  legProjXInRegions.AddEntry(projFailPt2, "p_{T} = 300 - inf  GeV","LP")

                  projFail.GetXaxis().SetRangeUser(400,4000)
                  projFailPt3 = projFail.ProjectionY("name3").Rebin(2)
                  projFailPt3.SetStats(0)
                  projFailPt3.SetMarkerStyle(20)
                  projFailPt3.SetLineColor(3)
                  projFailPt3.SetMarkerColor(3)
                  projFailPt3.DrawClone("SAMEP")

                  legProjXInRegions.AddEntry(projFailPt3, "p_{T} = 400 - inf GeV","LP")
                  legProjXInRegions.Draw("SAMEP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  canLog.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_FAIL_PtBins.png")
                  
                  # repeat the same for PASS now
                  canLog2 = ROOT.TCanvas("newname2-pass"+keyname2,"newname2-pass"+keyname2,800,800)
                  canLog2.SetLogy()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.9),obj.GetXaxis().FindBin(1.0)+1) # makes it the PASS region
                  projPass = obj.Project3D("YZ")
                  projPass.SetTitle("")
                  
                  legProjXInRegions.SetHeader("PASS region (F_{i}^{pixels} > 0.9)","C")

                  projPass.GetXaxis().SetRangeUser(240,4000)
                  projPassPt1 = projPass.ProjectionY("name1-pass").Rebin(2)
                  
                  projPassPt1.SetStats(0)
                  projPassPt1.SetMarkerStyle(20)
                  projPassPt1.SetLineColor(1)
                  projPassPt1.SetMarkerColor(1)
                  projPassPt1.SetMaximum(projPassPt1.GetMaximum()*100)
                  projPassPt1.DrawClone("SAMEP")

                  projPass.GetXaxis().SetRangeUser(320,4000)
                  projPassPt2 = projPass.ProjectionY("name2-pass").Rebin(2)
                  projPassPt2.SetStats(0)
                  projPassPt2.SetMarkerStyle(20)
                  projPassPt2.SetLineColor(2)
                  projPassPt2.SetMarkerColor(2)
                  projPassPt2.DrawClone("SAMEP")

                  projPass.GetXaxis().SetRangeUser(400,4000)
                  projPassPt3 = projPass.ProjectionY("name3-pass").Rebin(2)
                  projPassPt3.SetStats(0)
                  projPassPt3.SetMarkerStyle(20)
                  projPassPt3.SetLineColor(3)
                  projPassPt3.SetMarkerColor(3)
                  projPassPt3.DrawClone("SAMEP")

                  legProjXInRegions.Draw("SAMEP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  canLog2.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_PASS_PtBins.png")

                if ("PostS_ProbQNoL1VsFiStripsVsPt" in keyname2) :
                  canLog = ROOT.TCanvas("newnameF"+keyname2,"newnameF"+keyname2,800,800)
                  canLog.SetLogy()
                  legProjXInRegions =  ROOT.TLegend(.2,.75,.7,.9,"","brNDC")
                  legProjXInRegions.SetHeader("FAIL region (F_{i}^{pixels} < 0.9)","C")
                  legProjXInRegions.SetTextFont(42)
                  legProjXInRegions.SetTextSize(0.035)
                  legProjXInRegions.SetBorderSize(1);
                  legProjXInRegions.SetLineColor(0);
                  legProjXInRegions.SetLineStyle(1);
                  legProjXInRegions.SetLineWidth(1);
                  legProjXInRegions.SetFillColor(0);
                  legProjXInRegions.SetFillStyle(1001);
                  obj.GetXaxis().SetRange(1,obj.GetXaxis().FindBin(0.9)-1)
                 
                  projFail = obj.Project3D("YZ")
                  projFail.SetTitle("")
                  projFail.GetXaxis().SetRangeUser(100,4000)
                  projFailPt1 = projFail.ProjectionY("name1F").Rebin(2)
                  projFailPt1.GetYaxis().SetTitle("Normalized events")
                  numEvents1 = projFailPt1.Integral()
                  projFailPt1.SetStats(0)
                  projFailPt1.SetMarkerStyle(20)
                  projFailPt1.SetLineColor(1)
                  projFailPt1.SetMarkerColor(1)
                  if (numEvents1 > 0) : projFailPt1.Scale(1/numEvents1)
                  projFailPt1.SetMaximum(projFailPt1.GetMaximum()*50)
                  projFailPt1.DrawClone("SAMEP")
                  legProjXInRegions.AddEntry(projFailPt1, "p_{T} = 100 - inf GeV, #Events: " +str(round(numEvents1)),"LP")

                  projFail.GetXaxis().SetRangeUser(200,4000)
                  projFailPt2 = projFail.ProjectionY("name2F").Rebin(2)
                  numEvents2 = projFailPt2.Integral()
                  projFailPt2.SetStats(0)
                  projFailPt2.SetMarkerStyle(20)
                  projFailPt2.SetLineColor(2)
                  projFailPt2.SetMarkerColor(2)
                  if (numEvents2 > 0) : projFailPt2.Scale(1/numEvents2)
                  projFailPt2.DrawClone("SAMEP")
                  legProjXInRegions.AddEntry(projFailPt2, "p_{T} = 200 - inf GeV, #Events: " +str(round(numEvents2)),"LP")

                  projFail.GetXaxis().SetRangeUser(300,4000)
                  projFailPt3 = projFail.ProjectionY("name3F").Rebin(2)
                  numEvents3 = projFailPt3.Integral()
                  projFailPt3.SetStats(0)
                  projFailPt3.SetMarkerStyle(20)
                  projFailPt3.SetLineColor(3)
                  projFailPt3.SetMarkerColor(3)
                  if (numEvents3 > 0) : projFailPt3.Scale(1/numEvents3)
                  projFailPt3.DrawClone("SAMEP")

                  legProjXInRegions.AddEntry(projFailPt3, "p_{T} = 300 - inf GeV, #Events: " +str(round(numEvents3)),"LP")
                  legProjXInRegions.Draw("SAMEP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  canLog.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_FAIL_PtBins.png")
                  
                  # repeat the same for PASS now
                  legProjXInRegionsPass =  ROOT.TLegend(.2,.75,.7,.9,"","brNDC")
                  legProjXInRegionsPass.SetHeader("FAIL region (F_{i}^{pixels} < 0.9)","C")
                  legProjXInRegionsPass.SetTextFont(42)
                  legProjXInRegionsPass.SetTextSize(0.035)
                  legProjXInRegionsPass.SetBorderSize(1);
                  legProjXInRegionsPass.SetLineColor(0);
                  legProjXInRegionsPass.SetLineStyle(1);
                  legProjXInRegionsPass.SetLineWidth(1);
                  legProjXInRegionsPass.SetFillColor(0);
                  legProjXInRegionsPass.SetFillStyle(1001);
                  canLog2 = ROOT.TCanvas("newname2-passF"+keyname2,"newname2-passF"+keyname2,800,800)
                  canLog2.SetLogy()
                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.9),obj.GetXaxis().FindBin(1.0)+1) # makes it the PASS region
                  projPass = obj.Project3D("YZ")
                  projPass.SetTitle("")
                  
                  legProjXInRegionsPass.SetHeader("PASS region (F_{i}^{pixels} > 0.9)","C")

                  projPass.GetXaxis().SetRangeUser(100,4000)
                  projPassPt1 = projPass.ProjectionY("name1-passF").Rebin(2)
                  projPassPt1.GetYaxis().SetTitle("Normalized events")
                  numEvents1 = projPassPt1.Integral()
                  legProjXInRegionsPass.AddEntry(projPassPt1, "p_{T} = 100 - inf GeV, #Events: " +str(round(numEvents1)),"LP")
                  
                  projPassPt1.SetStats(0)
                  projPassPt1.SetMarkerStyle(20)
                  projPassPt1.SetLineColor(1)
                  projPassPt1.SetMarkerColor(1)
                  if (numEvents1 > 0) : projPassPt1.Scale(1/numEvents1)
                  projPassPt1.SetMaximum(projPassPt1.GetMaximum()*50)
                  projPassPt1.DrawClone("SAMEP")

                  projPass.GetXaxis().SetRangeUser(200,4000)
                  projPassPt2 = projPass.ProjectionY("name2-passF").Rebin(2)
                  numEvents2 = projPassPt2.Integral()
                  legProjXInRegionsPass.AddEntry(projPassPt2, "p_{T} = 200 - inf GeV, #Events: " +str(round(numEvents2)),"LP")
                  projPassPt2.SetStats(0)
                  projPassPt2.SetMarkerStyle(20)
                  projPassPt2.SetLineColor(2)
                  projPassPt2.SetMarkerColor(2)
                  if (numEvents2 > 0) : projPassPt2.Scale(1/numEvents2)
                  projPassPt2.DrawClone("SAMEP")

                  projPass.GetXaxis().SetRangeUser(300,4000)
                  projPassPt3 = projPass.ProjectionY("name3-passF").Rebin(2)
                  numEvents3 = projPassPt3.Integral()
                  legProjXInRegionsPass.AddEntry(projPassPt3, "p_{T} = 300 - inf GeV, #Events: " +str(round(numEvents3)),"LP")
                  projPassPt3.SetStats(0)
                  projPassPt3.SetMarkerStyle(20)
                  projPassPt3.SetLineColor(3)
                  projPassPt3.SetMarkerColor(3)
                  if (numEvents3 > 0) : projPassPt3.Scale(1/numEvents3)
                  projPassPt3.DrawClone("SAMEP")

                  legProjXInRegionsPass.Draw("SAMEP")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  canLog2.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_PASS_PtBins.png")


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
                if ("Calibration_GiTemplate" in keyname2) :
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
                # this maybe should go from bin to bin+1 ?
#                if ("IhVsLayer" in keyname2 or "IhVsLayer" in keyname2) :
#                  obj.SetMarkerStyle(20)
#                  if ("PostPreS_IasPixelIhVsLayer" in keyname2) :
#                      ratioOfhighIhPartOlowIh = ROOT.TH1F("RatioOfHighIhOverLowIh",";;Ratio of Ih>5 over Ih<5",7,0.,7.)
#                      for x in range(1,obj.GetNbinsX()) :
#                        lowIhPart = obj.Project3D("YZ").Integral(x,x,1,obj.Project3D("YZ").GetYaxis().FindBin(5.0))
#                        highIhPart = obj.Project3D("YZ").Integral(x,x,obj.Project3D("YZ").GetYaxis().FindBin(5.0),obj.Project3D("YZ").GetNbinsY())
#                        ratio = 0
#                        if (lowIhPart>0) :
#                          ratio = highIhPart/lowIhPart
#                        ratioOfhighIhPartOlowIh.SetBinContent(x,ratio)
#                      ratioOfhighIhPartOlowIh.Draw("COLZ")
#                      ratioOfhighIhPartOlowIh.SetStats(0)
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(1,"BPix L1")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(2,"BPix L2")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(3,"BPix L3")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(4,"BPix L4")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(5,"FPix D1")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(6,"FPix D2")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(7,"FPix D3")
#                      ratioOfhighIhPartOlowIh.GetYaxis().SetTitleOffset(1.5)
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(1,"BPix L1")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(2,"BPix L2")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(3,"BPix L3")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(4,"BPix L4")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(5,"FPix D1")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(6,"FPix D2")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(7,"FPix D3")
#                      obj.Project3D("YZ").GetYaxis().SetTitleOffset(0.8)
#                      obj.Project3D("YZ").GetYaxis().SetTitle("Ih")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(1,"BPix L1")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(2,"BPix L2")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(3,"BPix L3")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(4,"BPix L4")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(5,"FPix D1")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(6,"FPix D2")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(7,"FPix D3")
#                      obj.Project3D("XZ").GetYaxis().SetTitle("Ias")
#                      obj.Project3D("XZ").GetYaxis().SetTitleOffset(0.8)
#                  elif ("PostPreS_IasStripIhVsLayer" in keyname2) :
#                      ratioOfhighIhPartOlowIh = ROOT.TH1F("RatioOfHighIhOverLowIh",";;Ratio of Ih>5 over Ih<5",23,0.,23.)
#                      for x in range(1,obj.GetNbinsX()) :
#                        lowIhPart = obj.Project3D("YZ").Integral(x,x,1,obj.Project3D("YZ").GetYaxis().FindBin(5.0))
#                        highIhPart = obj.Project3D("YZ").Integral(x,x,obj.Project3D("YZ").GetYaxis().FindBin(5.0),obj.Project3D("YZ").GetNbinsY())
#                        ratio = 0
#                        if (lowIhPart>0) :
#                          ratio = highIhPart/lowIhPart
#                        ratioOfhighIhPartOlowIh.SetBinContent(x,ratio)
#                      ratioOfhighIhPartOlowIh.Draw("COLZ")
#                      ratioOfhighIhPartOlowIh.SetStats(0)
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(1,"TIB L1")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(2,"TIB L2")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(3,"TIB L3")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(4,"TIB L4")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(5,"TOB L1")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(6,"TOB L2")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(7,"TOB L3")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(8,"TOB L4")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(9,"TOB L5")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(10,"TOB L6")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(11,"TID D1")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(12,"TID D2")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(13,"TID D3")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(14,"TEC D1")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(15,"TEC D2")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(16,"TEC D3")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(17,"TEC D4")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(18,"TEC D5")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(19,"TEC D6")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(20,"TEC D7")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(21,"TEC D8")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(22,"TEC D9")
#                      ratioOfhighIhPartOlowIh.GetXaxis().SetBinLabel(23,"TEC D10")
#                      ratioOfhighIhPartOlowIh.GetYaxis().SetTitleOffset(1.5)
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(1,"TIB L1")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(2,"TIB L2")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(3,"TIB L3")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(4,"TIB L4")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(5,"TOB L1")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(6,"TOB L2")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(7,"TOB L3")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(8,"TOB L4")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(9,"TOB L5")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(10,"TOB L6")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(11,"TID D1")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(12,"TID D2")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(13,"TID D3")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(14,"TEC D1")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(15,"TEC D2")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(16,"TEC D3")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(17,"TEC D4")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(18,"TEC D5")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(19,"TEC D6")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(20,"TEC D7")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(21,"TEC D8")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(22,"TEC D9")
#                      obj.Project3D("YZ").GetXaxis().SetBinLabel(23,"TEC D10")
#                      obj.Project3D("YZ").GetYaxis().SetTitleOffset(0.8)
#                      obj.Project3D("YZ").GetYaxis().SetTitle("I_{h} (MeV/cm)")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(1,"TIB L1")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(2,"TIB L2")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(3,"TIB L3")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(4,"TIB L4")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(5,"TOB L1")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(6,"TOB L2")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(7,"TOB L3")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(8,"TOB L4")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(9,"TOB L5")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(10,"TOB L6")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(11,"TID D1")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(12,"TID D2")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(13,"TID D3")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(14,"TEC D1")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(15,"TEC D2")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(16,"TEC D3")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(17,"TEC D4")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(18,"TEC D5")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(19,"TEC D6")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(20,"TEC D7")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(21,"TEC D8")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(22,"TEC D9")
#                      obj.Project3D("XZ").GetXaxis().SetBinLabel(23,"TEC D10")
#                      obj.Project3D("XZ").GetYaxis().SetTitle("G_{i}^{Strips}")
#                      obj.Project3D("XZ").GetYaxis().SetTitleOffset(0.8)
##                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(1.0))
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_RatioOfLowIasHighIas.png")
#                  obj.Project3D("YZ").SetStats(0)
#                  obj.Project3D("YZ").Draw("COLZ")
#                  obj.Project3D("YZ").SetTitle("")
#                  tex4.Draw("SAME")
#                  tex5.Draw("SAME")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_allIas.png")
#                  obj.Project3D("XZ").SetStats(0)
#                  obj.Project3D("XZ").Draw("COLZ")
#                  obj.Project3D("XZ").SetTitle("")
#                  tex4.Draw("SAME")
#                  tex5.Draw("SAME")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_allIh.png")
#
#                  obj.Project3D("YZ").GetYaxis().UnZoom()
#                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.7),obj.GetXaxis().FindBin(1.0))
#                  obj.Project3D("YZ").Draw("COLZ")
#                  tex4.Draw("SAME")
#                  tex5.Draw("SAME")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_highIas.png")
#
#                  obj.GetXaxis().UnZoom()
#                  obj.GetXaxis().SetRange(obj.GetXaxis().FindBin(0.0),obj.GetXaxis().FindBin(0.7))
#                  projObj = obj.Project3D("YZ")
#                  if (projObj.GetEntries()==0) : continue
#                  projObj.Draw("COLZ")
#                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_lowIas.png")
#                  if ("PostPreS_IasPixelIhVsLayer" in keyname2):
#                    for i in range(7) :
#                      obj.Project3D("YZ").ProjectionY(newname,i,i+1,"e").Draw()
#                      can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_PixLayer"+str(i)+".png")
#                  elif ("PostPreS_IasStripIhVsLayer" in keyname2):
#                    for i in range(23) :
#                      obj.Project3D("YZ").ProjectionY(newname,i,i+1,"e").Draw()
#                      can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_StripLayer"+str(i)+".png")
#                else :
#                  obj.SetMarkerStyle(20)
#                  obj.GetXaxis().SetRange(bin,bin)
#                  obj.Project3D("ZY").Draw("COLZ")
              if ((obj.ClassName() == "TH2F" or obj.ClassName() == "TH2D") and not (keyname2 == "GenPtVsRecoPt" or "PreS_" in keyname2 or "CutFlow" in keyname2 or "N1_" in keyname2 or "_p_" in keyname2 or "_pterr" in keyname2 )):
                obj.SetTitle("")
                obj.SetMarkerStyle(20)
                projOb = obj.ProjectionY(newname,bin,bin,"e")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                projOb.Draw("COLZ L")
                can.SaveAs(name)

              elif ("Layer_CR_veryLowPt" in keyname2) :
                for x in range(1,obj.GetNbinsY()) :
                  CluDeDxForLayerX =  obj.ProjectionX(newname+"_Layer"+str(x),x,x,"e")
                  CluDeDxForLayerX.GetYaxis().SetTitle("Clusters")
                  CluDeDxForLayerX.Draw()
                  can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_Layer"+str(x)+".png")
                  #.Rebin(3)
              if (isData and (obj.ClassName() == "TH2F") and ("RunNumVsPixCluChargeAfter" in keyname2)) :
                x_values = []
                y_values = []
                for i in range(obj.GetNbinsX()+1) :
#                  if not (i==302): continue

                  chargeForRunI = obj.ProjectionY(newname+"_RunBin"+str(i),i,i,"e").Rebin(3)
                  if (chargeForRunI.Integral() > 10 and chargeForRunI.GetMean() > 10) :
                    
                    vf = Vavilov_Func()
                    fitForI = ROOT.TF1("fitForI"+str(i), vf, 0, 300, 5)
#                    fitForI = ROOT.TF1("fitForI", GaussWithLandau, 30, 300, 4)
#                    fitForI.SetParameters(100,20,20,100000) #,chargeForRunI.GetRMS(),chargeForRunI.GetRMS(),chargeForRunI.GetMaximum())
                    fitForI.SetParameters(
                        0.001*chargeForRunI.GetMean(),   # p[0]: kappa [0.01:10]
                        1.0,                                # p[1]: beta^2
                        0.9*chargeForRunI.GetMean(),      # p[2]: mpv
                        0.4*chargeForRunI.GetRMS(),       # p[3]: sigmaQ
                        0.5*chargeForRunI.GetEntries()    # p[4]: norm
                    )
                    fitForI.SetParLimits(0, 0.01, 10.0)
                    fitForI.SetParLimits(1, 0.9, 1)
                    fitForI.SetParLimits(2, 75, 100.);
                    fitForI.SetParLimits(3, 0, 10000.)
#                    fitForI.SetParLimits(4, 0.0, 100000.0)
#                    chargeForRunI.Fit("landau")
  
                    chargeForRunI.Fit(fitForI,"Q")
                    chargeForRunI.Draw()
                    fitForI.Draw("SAME")
                    fitForI.SetLineColor(ROOT.kRed)
                    fitForI.SetLineWidth(2)
                    if (fitForI.GetChisquare() / fitForI.GetNDF() < 100 ) :
                      x_values.append(obj.GetXaxis().GetBinCenter(i))
                      y_values.append(fitForI.GetParameter(2))
                    can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_RunNumber"+str(i)+".png")
                
                g_mpv = ROOT.TGraph(len(x_values), array.array('d', x_values), array.array('d', y_values))
                g_mpv.Draw("AP")
                g_mpv.SetMarkerStyle(20)
                g_mpv.SetMarkerColor(4)
                g_mpv.SetTitle("")
                g_mpv.SetMinimum(0.)
                g_mpv.GetYaxis().SetRangeUser(0.,200.)
                g_mpv.GetXaxis().SetTitleOffset(1)
                g_mpv.GetXaxis().SetTitle("Run number")
                g_mpv.GetXaxis().SetLabelSize(0.035)
                g_mpv.GetYaxis().SetTitle("MPV of norm cluster charge (e/um)")
                g_mpv.GetYaxis().SetTitleOffset(1.5)
                
                nPoints = g_mpv.GetN()
                yVals = g_mpv.GetY()
                mean = np.mean(yVals)
                variance = 0

                for i in range(nPoints):
                    diff = yVals[i] - mean
                    variance += diff * diff

                variance /= nPoints
                stdDev = ROOT.TMath.Sqrt(variance)
                stdDevOMean = round(stdDev / mean, 4)
                print("Standard deviation of y-axis values: ", stdDevOMean)
                
                yearSepLine = ROOT.TLine();
                yearSepLine.SetLineWidth(2);
                yearSepLine.SetLineStyle(ROOT.kDashed);
                yearSepLine.DrawLine(296500,g_mpv.GetMinimum(),296500,200)
                yearSepLine.DrawLine(307000,g_mpv.GetMinimum(),307000,200)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex5.Draw("SAME")
                tex2017 = ROOT.TLatex(0.225,0.8,"2017");
                tex2017.SetNDC();
                tex2017.SetTextFont(52)
                tex2017.SetTextSize(0.0485)
                tex2017.SetLineWidth(2)
                tex2017.Draw("SAME")
                yearSepLine.DrawLine(314500,g_mpv.GetMinimum(),314500,200)
                yearSepLine.DrawLine(325500,g_mpv.GetMinimum(),325500,200)
                tex2018 = ROOT.TLatex(0.65,0.8,"2018");
                tex2018.SetNDC();
                tex2018.SetTextFont(52)
                tex2018.SetTextSize(0.0485)
                tex2018.SetLineWidth(2)
                tex2018.Draw("SAME")
                meanLine = ROOT.TLine();
                meanLine.SetLineWidth(2);
                meanLine.SetLineColor(ROOT.kGreen)
#                meanLine.SetLineStyle(ROOT.kDashed);
                meanLine.DrawLine(295000,mean,327000,mean)
                texStdOverMean = ROOT.TLatex(0.5,0.94,"L"+str(keyname2[keyname2.find("SFsL")+4:keyname2.find("SFsL")+5]) + ": std / mean = " + str(round(100*stdDevOMean,2)) + "%");
                texStdOverMean.SetNDC();
                texStdOverMean.SetTextFont(52)
                texStdOverMean.SetTextSize(0.0485)
                texStdOverMean.SetLineWidth(2)
                texStdOverMean.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_MPVtrend.png")
                
                obj.GetXaxis().SetRangeUser(298000,324000)
                obj.ProfileX().Draw()
                obj.ProfileX().GetYaxis().SetRangeUser(0,200)
                obj.ProfileX().GetYaxis().SetTitle("Avg norm clu charge")
                obj.ProfileX().GetYaxis().SetTitleOffset(1.5)
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_profileX.png")
                
              if (isData and (obj.ClassName() == "TH2F") and ("Stab_CluDeDx" in keyname2)) :
                # ---------------------------- now the avg trend ------------------------------------
                obj.GetXaxis().SetRangeUser(298000,324000)
                obj.ProfileX().Draw()
                obj.ProfileX().SetStats(0)
                obj.ProfileX().GetYaxis().SetRangeUser(0,10)
                obj.ProfileX().GetYaxis().SetTitle("Avg dE/dx [MeV/cm]")
                obj.ProfileX().GetYaxis().SetTitleOffset(1.5)
                yearSepLine = ROOT.TLine();
                yearSepLine.SetLineWidth(2);
                yearSepLine.SetLineStyle(ROOT.kDashed);
                yearSepLine.DrawLine(307000,obj.ProfileX().GetMinimum(),307000,10)
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex5.Draw("SAME")
                tex2017 = ROOT.TLatex(0.225,0.8,"2017");
                tex2017.SetNDC();
                tex2017.SetTextFont(52)
                tex2017.SetTextSize(0.0485)
                tex2017.SetLineWidth(2)
                tex2017.Draw("SAME")
                yearSepLine.DrawLine(314500,obj.ProfileX().GetMinimum(),314500,10)
                tex2018 = ROOT.TLatex(0.65,0.8,"2018");
                tex2018.SetNDC();
                tex2018.SetTextFont(52)
                tex2018.SetTextSize(0.0485)
                tex2018.SetLineWidth(2)
                tex2018.Draw("SAME")
                meanLine = ROOT.TLine();
                meanLine.SetLineWidth(2);
                meanLine.SetLineColor(ROOT.kGreen)
                meanLine.DrawLine(298000,obj.GetMean(2),324000,obj.GetMean(2))
                axisYTitle = obj.GetYaxis().GetTitle()
                legendText = str(axisYTitle[0:axisYTitle.find(" Cluster")])
                if (keyname2 == "Stab_CluDeDxStripsLayer10_VsRun_CR_veryLowPt") : legendText = "TOB6"
                texStdOverMean = ROOT.TLatex(0.5,0.94,legendText + ": Mean = " + str(round(obj.GetMean(2),2)))
                texStdOverMean.SetNDC();
                texStdOverMean.SetTextFont(52)
                texStdOverMean.SetTextSize(0.0485)
                texStdOverMean.SetLineWidth(2)
                texStdOverMean.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_profileX.png")
                
              if (isData and (obj.ClassName() == "TH2F") and ("Stab_" in keyname2 and not "Clu" in keyname2)) :
                # ---------------------------- now the avg trend ------------------------------------
                axisYTitle = obj.GetYaxis().GetTitle()
                obj.GetXaxis().SetRangeUser(298000,324000)
                obj.ProfileX().Draw()
                obj.ProfileX().SetStats(0)
                maxi = 10 if "Ih" in keyname2 else 1.0
                if "Gi" in keyname2  : maxi = 0.1
                obj.ProfileX().GetYaxis().SetRangeUser(0,maxi)
                
                obj.ProfileX().GetYaxis().SetTitle("Avg " + axisYTitle)
                obj.ProfileX().GetYaxis().SetTitleOffset(1.5)
                yearSepLine = ROOT.TLine();
                yearSepLine.SetLineWidth(2);
                yearSepLine.SetLineStyle(ROOT.kDashed);
                yearSepLine.DrawLine(307000,obj.ProfileX().GetMinimum(),307000,obj.ProfileX().GetMaximum())
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex5.Draw("SAME")
                tex2017 = ROOT.TLatex(0.225,0.8,"2017");
                tex2017.SetNDC();
                tex2017.SetTextFont(52)
                tex2017.SetTextSize(0.0485)
                tex2017.SetLineWidth(2)
                tex2017.Draw("SAME")
                yearSepLine.DrawLine(314500,obj.ProfileX().GetMinimum(),314500,obj.ProfileX().GetMaximum())
                tex2018 = ROOT.TLatex(0.65,0.8,"2018");
                tex2018.SetNDC();
                tex2018.SetTextFont(52)
                tex2018.SetTextSize(0.0485)
                tex2018.SetLineWidth(2)
                tex2018.Draw("SAME")
                meanLine = ROOT.TLine();
                meanLine.SetLineWidth(2);
                meanLine.SetLineColor(ROOT.kGreen)
                meanLine.DrawLine(298000,obj.GetMean(2),324000,obj.GetMean(2))
                
                legendText = str(axisYTitle)
                texStdOverMean = ROOT.TLatex(0.5,0.94,legendText + ": Mean = " + str(round(obj.GetMean(2),2)))
                texStdOverMean.SetNDC();
                texStdOverMean.SetTextFont(52)
                texStdOverMean.SetTextSize(0.0485)
                texStdOverMean.SetLineWidth(2)
                texStdOverMean.Draw("SAME")
                tex2.Draw("SAME")
                tex5.Draw("SAME")
                can.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_avgX.png")
                  
                  
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
                  
                  if (blind and not "CR" in keyname2 and isData) :
                    for binIndex in range(16,50) :
                      projY1.SetBinContent(binIndex,0)
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


                  rp2.GetLowerRefXaxis().SetTitleSize(0.05)
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

#                if ("ProbQNoL1VsIas" in keyname2) :
#                  can2 = ROOT.TCanvas(newname+"2",newname+"2",800,800)
#                  can2.SetLogy()
#                  legIasForProbQSlice =  ROOT.TLegend(.30,.80,.80,.90,"","brNDC")
#                  projY1 = obj.ProjectionY("IasForProbQSlice_Slice1",obj.GetXaxis().FindBin(0.9),obj.GetXaxis().FindBin(1.0),"e")
#                  projY1.SetStats(0)
#                  projY1.SetMarkerColor(1)
#                  projY1.SetLineColor(1)
#                  projY1.SetMarkerStyle(20)
#                  projY1.Draw("SAME")
#                  projY1.GetYaxis().SetTitle("Tracks / 0.05")
#                  projY1.GetYaxis().SetTitleOffset(1.5)
#
#                  projY2 = obj.ProjectionY("IasForProbQSlice_Slice2",obj.GetXaxis().FindBin(0.3),obj.GetXaxis().FindBin(0.9),"e")
#                  projY2.SetMarkerStyle(20)
#                  projY2.SetMarkerColor(2)
#                  projY2.SetLineColor(2)
#                  projY2.Draw("SAME")
#
#                  numTracks1 = projY1.Integral()
#                  numTracks2 = projY2.Integral()
#
#                  if (blind and not "CR" in keyname2 and isData) :
#                    for binIndex in range(1,50) :
#                      projY1.SetBinContent(binIndex,0)
#                  projY1.SetMaximum(projY2.GetMaximum()*100)
#                  projY1.SetMinimum(0.1)
#
#                  rp2 = ROOT.TRatioPlot(projY1,projY2,"divsym") #, "diffsigerrasym"
#
#                  rp2.SetH1DrawOpt("P");
#                  rp2.SetH2DrawOpt("P");
#
#                  rp2.Draw()
#
#                  rp2.SetLeftMargin(0.13);
#                  rp2.SetRightMargin(0.05);
#                  rp2.SetUpTopMargin(0.1);
#                  rp2.SetLowTopMargin(0.02);
#                  rp2.SetLowBottomMargin(0.35);
#
#
#                  rp2.GetLowerRefGraph().SetMinimum(0)
#                  rp2.GetLowerRefGraph().SetMaximum(3.5);
#                  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
#                  #rp.GetLowerRefGraph().SetLineColor(0) #0
#                  rp2.GetLowerRefGraph().SetMarkerStyle(20)
#                  rp2.GetLowerRefGraph().SetMarkerSize(1);
#                  rp2.GetLowYaxis().SetNdivisions(505);
#                  rp2.GetLowerRefYaxis().SetTitle("Ratio");
#                  rp2.GetLowerRefYaxis().SetTitleSize(0.05);
#                  rp2.GetLowerRefYaxis().SetTitleOffset(1);
#                  rp2.GetLowerRefYaxis().SetLabelSize(0.035);
#
#
#                  rp2.GetLowerRefXaxis().SetTitleSize(0.05);
#                  rp2.GetLowerRefXaxis().SetTitleOffset(0.8);
#                  rp2.GetLowerRefXaxis().SetLabelSize(0.035);
#                  legIasForProbQSlice.SetTextFont(42)
#                  legIasForProbQSlice.SetTextSize(0.035)
#                  legIasForProbQSlice.SetBorderSize(1);
#                  legIasForProbQSlice.SetLineColor(0);
#                  legIasForProbQSlice.SetLineStyle(1);
#                  legIasForProbQSlice.SetLineWidth(1);
#                  legIasForProbQSlice.SetFillColor(0);
#                  legIasForProbQSlice.SetFillStyle(1001);
#                  legIasForProbQSlice.AddEntry(projY1,"F_{i}^{Pixels} (0.9-1.0), #Tracks: " +str(round(numTracks1)),"LP")
#                  legIasForProbQSlice.AddEntry(projY2,"F_{i}^{Pixels}  (0.3-0.9), #Tracks: " +str(round(numTracks2)),"LP")
#                  tex2.Draw("SAME")
#                  tex3.Draw("SAME")
#                  tex4.Draw("SAME")
#                  tex5.Draw("SAME")
#                  legIasForProbQSlice.Draw("SAME")
#                  can2.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 + "_NotNorm_ProbQSlice.png")
#
#                  can3 = ROOT.TCanvas(newname+"2",newname+"2",800,800)
#                  can3.SetLogy()
#                  legIasForProbQSlice =  ROOT.TLegend(.30,.80,.80,.90,"","brNDC")
#                  projY1 = obj.ProjectionX("IasForProbQSlice_Slice1",1,obj.GetYaxis().FindBin(0.1),"e")
#                  projY1.SetStats(0)
#                  projY1.SetMarkerColor(1)
#                  projY1.SetLineColor(1)
#                  projY1.SetMarkerStyle(20)
#                  projY1.Draw("SAME")
#                  projY1.GetYaxis().SetTitle("Normalized Tracks / 0.05")
#                  projY1.GetYaxis().SetTitleOffset(1.5)
#
#                  projY2 = obj.ProjectionX("IasForProbQSlice_Slice2",obj.GetYaxis().FindBin(0.1),obj.GetYaxis().FindBin(1.0),"e")
#                  projY2.SetMarkerStyle(20)
#                  projY2.SetMarkerColor(2)
#                  projY2.SetLineColor(2)
#                  projY2.Draw("SAME")
#
#                  numTracks3 = projY1.Integral()
#                  numTracks4 = projY2.Integral()
#
#                  if (blind and not "CR" in keyname2 and isData) :
#                    projY1.SetBinContent(19,0)
#                    projY1.SetBinContent(20,0)
#
#                  projY1.SetMaximum(projY1.GetMaximum()*100)
#
#                  rp3 = ROOT.TRatioPlot(projY1,projY2,"divsym") #, "diffsigerrasym"
#
#                  rp3.SetH1DrawOpt("P");
#                  rp3.SetH2DrawOpt("P");
#
#                  rp3.Draw()
#
#                  rp3.SetLeftMargin(0.13);
#                  rp3.SetRightMargin(0.05);
#                  rp3.SetUpTopMargin(0.1);
#                  rp3.SetLowTopMargin(0.02);
#                  rp3.SetLowBottomMargin(0.35);
#
#
#                  rp3.GetLowerRefGraph().SetMinimum(0)
#                  rp3.GetLowerRefGraph().SetMaximum(3.5)
#                  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
#                  #rp.GetLowerRefGraph().SetLineColor(0) #0
#                  rp3.GetLowerRefGraph().SetMarkerStyle(20)
#                  rp3.GetLowerRefGraph().SetMarkerSize(1);
#                  rp3.GetLowYaxis().SetNdivisions(505);
#                  rp3.GetLowerRefYaxis().SetTitle("Ratio");
#                  rp3.GetLowerRefYaxis().SetTitleSize(0.05);
#                  rp3.GetLowerRefYaxis().SetTitleOffset(1);
#                  rp3.GetLowerRefYaxis().SetLabelSize(0.035);
#
#
#                  rp3.GetLowerRefXaxis().SetTitleSize(0.05);
#                  rp3.GetLowerRefXaxis().SetTitleOffset(0.8);
#                  rp3.GetLowerRefXaxis().SetLabelSize(0.035);
#                  legIasForProbQSlice.SetTextFont(42)
#                  legIasForProbQSlice.SetTextSize(0.035)
#                  legIasForProbQSlice.SetBorderSize(1);
#                  legIasForProbQSlice.SetLineColor(0);
#                  legIasForProbQSlice.SetLineStyle(1);
#                  legIasForProbQSlice.SetLineWidth(1);
#                  legIasForProbQSlice.SetFillColor(0);
#                  legIasForProbQSlice.SetFillStyle(1001);
#                  legIasForProbQSlice.AddEntry(projY1,"G_{i}^{Strips} (0.0-0.1), #Tracks: "+str(round(numTracks3)),"LP")
#                  legIasForProbQSlice.AddEntry(projY2,"G_{i}^{Strips} (0.1-1.0), #Tracks: "+str(round(numTracks4)),"LP")
#                  tex2.Draw("SAME")
#                  tex3.Draw("SAME")
#                  tex4.Draw("SAME")
#                  tex5.Draw("SAME")
#                  legIasForProbQSlice.Draw("SAME")
#                  can3.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 + "_NotNorm_IasSlice.png")
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
              
              if ("FracSat" in keyname2) :
                print("Mean: ",obj.GetMean())
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
                normFactor = obj.GetBinContent(1) + obj.GetBinContent(4)
                obj.Scale(1/normFactor)
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
              elif ("K_and_C_Kin_Mass" in keyname2) :
                print("TADA")
                tex2.Draw("SAME")
                tex3.Draw("SAME")
                tex4.Draw("SAME")
                tex5.Draw("SAME")
                tex3.Draw("SAME")
                obj.GetYaxis().SetTitle("Track / 0.05 GeV")
                
                g1 = ROOT.TF1( 'g1', 'gaus',  0.3,  0.7 )
                g2 = ROOT.TF1( 'g2', 'gaus',  0.7, 1.5 )
                g3 = ROOT.TF1( 'g3', 'gaus', 1.5, 2.5 )
                g1.SetLineColor(ROOT.kBlue)
                g2.SetLineColor(ROOT.kRed)
                g3.SetLineColor(ROOT.kGreen)
                obj.Fit(g1, "R")
                obj.Fit(g2, "R+")
                obj.Fit(g3, "R+")
                
                legMass =  ROOT.TLegend(.45,.75,.80,.9,"","brNDC")
                legMass.SetTextFont(42)
                legMass.SetTextSize(0.035)
                legMass.SetBorderSize(1);
                legMass.SetLineColor(1);
                legMass.SetLineStyle(1);
                legMass.SetLineWidth(1);
                legMass.SetFillColor(0);
                legMass.SetFillStyle(1001);
                legMass.AddEntry(g1,"#mu = "+str(round(g1.GetParameter(1),2)) + " #pm "+str(round(g1.GetParameter(2),2)),"LP")
                legMass.AddEntry(g2,"#mu = "+str(round(g2.GetParameter(1),2)) + " #pm "+str(round(g2.GetParameter(2),2)),"LP")
                legMass.AddEntry(g3,"#mu = "+str(round(g3.GetParameter(1),2)) + " #pm "+str(round(g3.GetParameter(2),2)),"LP")
                legMass.Draw("SAME")
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
                
              if (keyname2=="CutFlowProbQ" or keyname2=="CutFlowPfType" or keyname2=="CutFlowEta" or keyname2=="CutFlowEoP" or "Vs" in keyname2) :
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

Mass = f.Get("/analyzer/BaseName/Mass")
Mass_wPred = f.Get("/analyzer/BaseName/Pred_Mass_CB")
if Mass_wPred :
  tex5m = ROOT.TLatex(0.07,0.01,fileVersion)
  tex5m.SetNDC();
  tex5m.SetTextFont(52);
  tex5m.SetTextSize(0.0185);
  tex5m.SetLineWidth(2);
  name = fileName[0:-5] + "_Bin" + str(bin)+ "/"
  if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
  massBins = [10.,50.,100.,200.,300.,500.,1000.,4000.]
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
  Mass_projY.GetYaxis().SetLabelSize(0.03)
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
  Mass_wPred_projY.GetYaxis().SetLabelSize(0.03)
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

  tex4 = ROOT.TLatex(0.7,0.93,"K-S test v2: "+str(round(KSvalue2,4)));
  tex4.SetNDC();
  tex4.SetTextFont(52);
  tex4.SetTextSize(0.0485);
  tex4.SetLineWidth(2);

  cMass_projY = ROOT.TCanvas('cMass_projY', 'cMass_projY',800,800)

  rp = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY, "diffsigerrasym")

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
  rp.GetLowerRefGraph().SetMinimum(-4);
  rp.GetLowerRefGraph().SetMaximum(4);
  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp.GetLowerRefGraph().SetLineColor(0) #0
  rp.GetLowerRefGraph().SetMarkerStyle(20)
  rp.GetLowerRefGraph().SetMarkerSize(1);
  rp.GetLowYaxis().SetNdivisions(505);
  rp.GetLowerRefYaxis().SetTitle("Pull");
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
  tex5m.Draw("SAME")
  
  name = newFileDir + "/cMass.png"
  cMass_projY.SaveAs(name)
  
  
  #############################################################################
  cMass_projY_log = ROOT.TCanvas('cMass_projY_log', 'cMass_projY_log',800,800)
  cMass_projY_log.SetLogy()

  rp2 = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY, "diffsigerrasym")

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
  Mass_projY.SetMinimum(0.000001);
#  Mass_projY.GetYaxis().SetRangeUser(0.1,100)
  rp2.GetLowerRefGraph().SetMinimum(-4);
  rp2.GetLowerRefGraph().SetMaximum(4);
  #rp2.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp2.GetLowerRefGraph().SetLineColor(0) #0
  rp2.GetLowerRefGraph().SetMarkerStyle(20)
  rp2.GetLowerRefGraph().SetMarkerSize(1);
  rp2.GetLowYaxis().SetNdivisions(505);
  rp2.GetLowerRefYaxis().SetTitle("Pull");
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
  tex5m.Draw("SAME")
  
  name = newFileDir + "/cMass_log.png"
  cMass_projY_log.SaveAs(name)
  
  cMassOrig_projY = ROOT.TCanvas('cMassOrig_projY', 'cMassOrig_projY',800,800)
  rp0 = ROOT.TRatioPlot(Mass_projY_NotRebinned,Mass_wPred_projY_NotRebinned, "diffsigerrasym")
  rp0.Draw()
  rp0.SetH1DrawOpt("P");
  rp0.SetH2DrawOpt("P");
  rp0.SetLeftMargin(0.13);
  rp0.SetRightMargin(0.05);
  rp0.SetUpTopMargin(0.1);
  rp0.SetLowTopMargin(0.02);
  rp0.SetLowBottomMargin(0.35);
  rp0.GetLowerRefGraph().SetMinimum(-4);
  rp0.GetLowerRefGraph().SetMaximum(4);
  #rp0.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
  #rp0.GetLowerRefGraph().SetLineColor(0) #0
  rp0.GetLowerRefGraph().SetMarkerStyle(20)
  rp0.GetLowerRefGraph().SetMarkerSize(1);
  rp0.GetLowYaxis().SetNdivisions(505);
  rp0.GetLowerRefYaxis().SetTitle("Pull");
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
  tex5m.Draw("SAME")
  cMassOrig_projY.SaveAs(newFileDir + "/cMass_NotRebinned.png")
  
  cMassOrig_projY_log = ROOT.TCanvas('cMassOrig_projY_log', 'cMassOrig_projY_log',800,800)
  cMassOrig_projY_log.SetLogy()
  rp0.Draw()
  legMass.Draw("SAME")
  tex2.Draw("SAME")
  tex3.Draw("SAME")
  tex4.Draw("SAME")
  tex5m.Draw("SAME")
  cMassOrig_projY_log.SaveAs(newFileDir + "/cMass_logy_NotRebinned.png")
  


os.system("cp forWebpage/* "+newFileDir+"/.")
os.system("cp forWebpage/.htaccess "+newFileDir+"/.")
print("scp -r "+ newFileDir + " tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/.")
