import ROOT, sys, os, time, re
import numpy as np
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog")
(opt,args) = parser.parse_args()

ROOT.gROOT.SetStyle("Plain")
#ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.SetBatch(True)


f1 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_1800_wProbQ_CodeV16p8_v1.root")
f2 = ROOT.TFile.Open("crab_Analysis_2018_TTToHadronic_wProbQ_CodeV16p8_v1.root")
f3 = ROOT.TFile.Open("crab_Analysis_2018_TTToSemiLeptonic_wProbQ_CodeV16p8_v1.root")
f4 = ROOT.TFile.Open("crab_Analysis_2018_TTTo2L2Nu_wProbQ_CodeV16p8_v1.root")
f5 = ROOT.TFile.Open("crab_Analysis_2018_QCDwPt1000_wProbQ_CodeV16p8_v1.root")
f6 = ROOT.TFile.Open("crab_Analysis_2018_WJetsToLNu_wProbQ_CodeV16p8_v1.root")
f7 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_1400_CodeV8_v1.root")
f8 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_1800_CodeV8_v1.root")
f9 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2000_CodeV8_v1.root")
f10 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2200_CodeV8_v1.root")
f11 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2400_CodeV8_v1.root")
f12 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2600_CodeV8_v1.root")

fileOut = open("BetaGammaValuesMPV.txt", "a")

dirs = []
for i in range(0, f1.GetListOfKeys().GetEntries()):
  # Remove/modify unnecessary stuff from the name of the plot that was required by SmartHistos to ditinguish plots
  dirname = f1.GetListOfKeys().At(i).GetName()
  curr_dir = f1.GetDirectory(dirname)
# print("dirname: "+dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      # Match the plot of interest
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = f1.GetDirectory(dirname+"/"+keyname)
#                    print("keyname: "+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          # The plot should be TCanvas
          newname = dirname + "/" + keyname+ "/" + keyname2
#         print("newname: "+newname)
          obj1 = f1.Get(newname)
          obj2 = f2.Get(newname)
          obj3 = f3.Get(newname)
          obj4 = f4.Get(newname)
          obj5 = f5.Get(newname)
          obj6 = f6.Get(newname)
#          obj7 = f7.Get(newname)
#          obj8 = f8.Get(newname)
#          obj9 = f9.Get(newname)
#          obj10 = f10.Get(newname)
#          obj11 = f11.Get(newname)
#          obj12 = f12.Get(newname)
          if obj1.InheritsFrom("TObject"):
              can = obj1
              can = ROOT.TCanvas(newname)
              # bin 3: pt>60 and I_as > 0.05
              # bin 25: pt>65 and I_as > 0.175
              # bin 28: pt>65 and I_as > 0.3
              name = "CompareHSCPandBcg/" + keyname2 +  ".png"
#                 print(name)
              if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
              if not (obj1.ClassName() == "TH1F"):
                  continue
              if not (obj1.Integral()) :
                      continue
              if ("BS_" in keyname2 or keyname2 == "CutFlow" or keyname2 == "CutFlowProbQLast" or keyname2 == "CutFlowProbQFirst") :
#                  can.SetLogy()
                  obj1.SetMarkerColor(1)
                  obj2.SetMarkerColor(2)
                  obj3.SetMarkerColor(3)
                  obj4.SetMarkerColor(4)
                  obj5.SetMarkerColor(5)
                  obj6.SetMarkerColor(6)
#                  obj7.SetMarkerColor(7)
#                  obj8.SetMarkerColor(8)
#                  obj9.SetMarkerColor(9)
#                  obj10.SetMarkerColor(10)
#                  obj11.SetMarkerColor(11)
#                  obj12.SetMarkerColor(12)
                  
                  obj1.SetLineColor(1)
                  obj2.SetLineColor(2)
                  obj3.SetLineColor(3)
                  obj4.SetLineColor(4)
                  obj5.SetLineColor(5)
                  obj6.SetLineColor(6)
#                  obj7.SetLineColor(7)
#                  obj8.SetLineColor(8)
#                  obj9.SetLineColor(9)
#                  obj10.SetLineColor(10)
#                  obj11.SetLineColor(11)
#                  obj12.SetLineColor(12)
                  
                  
                  obj1.SetMarkerStyle(20)
                  obj2.SetMarkerStyle(20)
                  obj3.SetMarkerStyle(20)
                  obj4.SetMarkerStyle(20)
                  obj5.SetMarkerStyle(20)
                  obj6.SetMarkerStyle(20)
#                  obj7.SetMarkerStyle(20)
#                  obj8.SetMarkerStyle(20)
#                  obj9.SetMarkerStyle(20)
#                  obj10.SetMarkerStyle(20)
#                  obj11.SetMarkerStyle(20)
#                  obj12.SetMarkerStyle(20)
                  
#                  obj1.SetTitle("")
#                  obj1.GetXaxis().SetTitleSize(0.05)
#                  obj1.GetXaxis().SetTitleOffset(1)
#                  obj1.GetYaxis().SetTitleSize(0.05)
#                  obj1.GetYaxis().SetTitleOffset(1)
                  obj1.SetStats(0)


                  if (keyname2 == "CutFlow" or keyname2 == "CutFlowProbQLast" or keyname2 == "CutFlowProbQFirst") :
                    obj1.Scale(1/obj1.GetMaximum())
                    obj2.Scale(1/obj2.GetMaximum())
                    obj3.Scale(1/obj3.GetMaximum())
                    obj4.Scale(1/obj4.GetMaximum())
                    obj5.Scale(1/obj5.GetMaximum())
                    obj6.Scale(1/obj6.GetMaximum())
                  else:
                    obj1.Scale(1/obj1.Integral(), "width")
                    obj2.Scale(1/obj2.Integral(), "width")
                    obj3.Scale(1/obj3.Integral(), "width")
                    obj4.Scale(1/obj4.Integral(), "width")
                    obj5.Scale(1/obj5.Integral(), "width")
                    obj6.Scale(1/obj6.Integral(), "width")
                    
                  max = np.maximum(obj2.GetMaximum(), obj1.GetMaximum())
                  obj1.SetMaximum(max*1.5)
                  
                  legMass =  ROOT.TLegend(.55,.55,.80,.9,"","brNDC")
                  legMass.SetTextFont(42)
                  legMass.SetTextSize(0.035)
                  legMass.SetBorderSize(1);
                  legMass.SetLineColor(1);
                  legMass.SetLineStyle(1);
                  legMass.SetLineWidth(1);
                  legMass.SetFillColor(0);
                  legMass.SetFillStyle(1001);
#                  legMass.AddEntry(obj12,"Mass=2600 GeV","LP")
#                  legMass.AddEntry(obj11,"Mass=2400 GeV","LP")
#                  legMass.AddEntry(obj10,"Mass=2200 GeV","LP")
#                  legMass.AddEntry(obj9,"Mass=2000 GeV","LP")
#                  legMass.AddEntry(obj8,"Mass=1800 GeV","LP")
#                  legMass.AddEntry(obj7,"Mass=1400 GeV","LP")
                  legMass.AddEntry(obj1,"HSCP 1800 GeV","LP")
                  legMass.AddEntry(obj2,"TTToHadronic","LP")
                  legMass.AddEntry(obj3,"TTToSemiLeptonic","LP")
                  legMass.AddEntry(obj4,"TTTo2L2N","LP")
                  legMass.AddEntry(obj5,"QCDwPt1000","LP")
                  legMass.AddEntry(obj6,"WJetsToLNu","LP")
                  
                  if (keyname2 == "CutFlow") :
                    obj1.GetXaxis().SetBinLabel(1,"Trigger")
                    obj1.GetXaxis().SetBinLabel(2,"Eta")
                    obj1.GetXaxis().SetBinLabel(3,"pT")
                    obj1.GetXaxis().SetBinLabel(4,"NumHits")
                    obj1.GetXaxis().SetBinLabel(5,"NumPixHits")
                    obj1.GetXaxis().SetBinLabel(6,"ValidFract")
                    obj1.GetXaxis().SetBinLabel(7,"NumDeDx")
                    obj1.GetXaxis().SetBinLabel(8,"ProbQ")
                    obj1.GetXaxis().SetBinLabel(9,"ProbXY")
                    obj1.GetXaxis().SetBinLabel(10,"HighPurity")
                    obj1.GetXaxis().SetBinLabel(11,"Chi2oDOF")
                    obj1.GetXaxis().SetBinLabel(12,"EoP")
                    obj1.GetXaxis().SetBinLabel(13,"dz")
                    obj1.GetXaxis().SetBinLabel(14,"dxy")
                    obj1.GetXaxis().SetBinLabel(15,"pTerrOverpT")
                    obj1.GetXaxis().SetBinLabel(16,"TKIso")
                    obj1.GetXaxis().SetBinLabel(17,"Ih")
                    obj1.GetXaxis().SetBinLabel(18,"MuStat")
                    obj1.GetXaxis().SetBinLabel(19,"PhiTOF")
                    obj1.GetXaxis().SetBinLabel(20,"EtaTOF")
                  
                  if (keyname2 == "CutFlowProbQFirst") :
                    obj1.GetXaxis().SetBinLabel(1,"Trigger")
                    obj1.GetXaxis().SetBinLabel(2,"pT")
                    obj1.GetXaxis().SetBinLabel(3,"ProbQ")
                    obj1.GetXaxis().SetBinLabel(4,"Eta")
                    obj1.GetXaxis().SetBinLabel(5,"NumHits")
                    obj1.GetXaxis().SetBinLabel(6,"NumPixHits")
                    obj1.GetXaxis().SetBinLabel(7,"ValidFract")
                    obj1.GetXaxis().SetBinLabel(8,"NumDeDx")
                    obj1.GetXaxis().SetBinLabel(9,"ProbXY")
                    obj1.GetXaxis().SetBinLabel(10,"HighPurity")
                    obj1.GetXaxis().SetBinLabel(11,"Chi2oDOF")
                    obj1.GetXaxis().SetBinLabel(12,"EoP")
                    obj1.GetXaxis().SetBinLabel(13,"dz")
                    obj1.GetXaxis().SetBinLabel(14,"dxy")
                    obj1.GetXaxis().SetBinLabel(15,"pTerrOverpT")
                    obj1.GetXaxis().SetBinLabel(16,"TKIso")
                    obj1.GetXaxis().SetBinLabel(17,"Ih")
                    obj1.GetXaxis().SetBinLabel(18,"MuStat")
                    obj1.GetXaxis().SetBinLabel(19,"PhiTOF")
                    obj1.GetXaxis().SetBinLabel(20,"EtaTOF")
                    
                  if (keyname2 == "CutFlowProbQLast") :
                    obj1.GetXaxis().SetBinLabel(1,"Trigger")
                    obj1.GetXaxis().SetBinLabel(2,"Eta")
                    obj1.GetXaxis().SetBinLabel(3,"pT")
                    obj1.GetXaxis().SetBinLabel(4,"NumHits")
                    obj1.GetXaxis().SetBinLabel(5,"NumPixHits")
                    obj1.GetXaxis().SetBinLabel(6,"ValidFract")
                    obj1.GetXaxis().SetBinLabel(7,"NumDeDx")
                    obj1.GetXaxis().SetBinLabel(8,"ProbXY")
                    obj1.GetXaxis().SetBinLabel(9,"HighPurity")
                    obj1.GetXaxis().SetBinLabel(10,"Chi2oDOF")
                    obj1.GetXaxis().SetBinLabel(11,"EoP")
                    obj1.GetXaxis().SetBinLabel(12,"dz")
                    obj1.GetXaxis().SetBinLabel(13,"dxy")
                    obj1.GetXaxis().SetBinLabel(14,"pTerrOverpT")
                    obj1.GetXaxis().SetBinLabel(15,"TKIso")
                    obj1.GetXaxis().SetBinLabel(16,"Ih")
                    obj1.GetXaxis().SetBinLabel(17,"ProbQ")
                    obj1.GetXaxis().SetBinLabel(18,"MuStat")
                    obj1.GetXaxis().SetBinLabel(19,"PhiTOF")
                    obj1.GetXaxis().SetBinLabel(20,"EtaTOF")
                    
                  obj1.Draw("COLZ L")
                  obj2.Draw("SAME")
                  obj3.Draw("SAME")
                  obj4.Draw("SAME")
                  obj5.Draw("SAME")
                  obj6.Draw("SAME")
#                  obj10.Draw("SAME")
#                  obj9.Draw("SAME")
#                  obj8.Draw("SAME")
#                  obj7.Draw("SAME")
                  

                  
                  legMass.Draw("SAME")
              else :
                  continue
              can.SaveAs(name)
              #can.SaveAs(name.replace(".png",".pdf"))
              #can.SaveAs(name.replace(".png",".C"))
              can.Close()
