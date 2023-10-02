import ROOT, sys, os, time, re
import numpy as np
from ctypes import c_double as double
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog fileName.root")
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
bin = 3
# bin 3: pt>60 and I_as > 0.05
# bin 25: pt>65 and I_as > 0.175
# bin 28: pt>65 and I_as > 0.3

blind = True
#blind = False

print("Filename: "+fileName)
input_file = fileName

ProjBin = 3
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
          
          
#          tex2 = ROOT.TLatex(0.13,0.94,"CMS");
          tex2 = ROOT.TLatex(0.23,0.94,"CMS");#if there is 10^x
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

          #tex3 = ROOT.TLatex(0.27,0.96,"Simulation"); # for square plots
          tex3 = ROOT.TLatex(0.37,0.94,"Internal"); #if there is 10^x
#          tex3 = ROOT.TLatex(0.27,0.94,"Internal");
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
          elif ("Calibration" in keyname2) :
            tex4 = ROOT.TLatex(0.6,0.95,"After calibration cuts")
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
              if not ("GiTemplate" in keyname2) : continue
              if (obj.ClassName() == "TH3F" or obj.ClassName() == "TH3D"):
                obj.SetTitle("")
                if ("Calibration_GiTemplate" in keyname2) :
                  projX = obj.ProjectionX()
                  projX.SetTitle("")
                  projX.SetStats(0)
                  projX.SetMarkerStyle(20)
                  projX.GetYaxis().SetTitle("Clusters")
                  projX.GetYaxis().SetTitleOffset(1.9)
                  projX.GetXaxis().SetRange(1,15)
                  projX.Draw("COLZ")
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
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
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
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
                  tex2.Draw("SAME")
                  tex3.Draw("SAME")
                  tex4.Draw("SAME")
                  tex5.Draw("SAME")
                  legGiCalib.Draw("SAME")
                  can4.SaveAs(fileName[0:-5] + "_Bin" + str(bin)+ "/" + keyname2 +  "_ProjZ.png")
              
name = fileName[fileName.find("Run")+3:fileName.find("_CodeV")]
version = "v5"
outFileName = str('template_{}_{}.root'.format(name,version))
fOut = ROOT.TFile.Open(outFileName,'UPDATE')
templ0 = f.Get('HSCParticleAnalyzer/BaseName/Calibration_GiTemplate')
templ1 = f.Get('HSCParticleAnalyzer/BaseName/Calibration_GiTemplate_PU_1')
templ2 = f.Get('HSCParticleAnalyzer/BaseName/Calibration_GiTemplate_PU_2')
templ3 = f.Get('HSCParticleAnalyzer/BaseName/Calibration_GiTemplate_PU_3')
templ4 = f.Get('HSCParticleAnalyzer/BaseName/Calibration_GiTemplate_PU_4')
templ5 = f.Get('HSCParticleAnalyzer/BaseName/Calibration_GiTemplate_PU_5')

templ0.Write()
templ1.Write()
templ2.Write()
templ3.Write()
templ4.Write()
templ5.Write()
fOut.Close()

os.system("mv " + outFileName + " " + fileName[0:-5] + "_Bin" + str(bin))

print("scp "+fileName[0:-5] + "_Bin" + str(bin)+"/"+outFileName + " vami@ui3.kfki.hu:/data/vami/projects/HSCP/v1/CMSSW_10_6_30/src/SUSYBSMAnalysis/HSCP/data/.")
