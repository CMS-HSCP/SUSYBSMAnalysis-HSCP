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


f1 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_100_CodeV17_v1.root")
f2 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_200_CodeV17_v1.root")
f3 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_400_CodeV17_v1.root")
f4 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_500_CodeV17_v1.root")
f5 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_800_CodeV17_v1.root")
f6 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_1000_CodeV17_v1.root")
f7 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_1400_CodeV17_v1.root")
f8 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_1800_CodeV17_v1.root")
f9 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2000_CodeV17_v1.root")
f10 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2200_CodeV17_v1.root")
f11 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2400_CodeV17_v1.root")
f12 = ROOT.TFile.Open("crab_Analysis_HSCPgluinoCentral_M_2600_CodeV17_v1.root")

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
          obj7 = f7.Get(newname)
          obj8 = f8.Get(newname)
          obj9 = f9.Get(newname)
          obj10 = f10.Get(newname)
          obj11 = f11.Get(newname)
          obj12 = f12.Get(newname)
          if obj1.InheritsFrom("TObject"):
              can = obj1
              can = ROOT.TCanvas(newname)
              # bin 3: pt>60 and I_as > 0.05
              # bin 25: pt>65 and I_as > 0.175
              # bin 28: pt>65 and I_as > 0.3
              name = "CompareHSCP/" + keyname2 +  ".png"
#                 print(name)
              if not os.path.exists(os.path.dirname(name)): os.makedirs(os.path.dirname(name))
              if not (obj1.ClassName() == "TH1F"):
                  continue
              if ("BS_" in keyname2 or keyname2 == "genlevelbeta") :
#                  can.SetLogy()
                  obj1.SetMarkerColor(1)
                  obj2.SetMarkerColor(2)
                  obj3.SetMarkerColor(3)
                  obj4.SetMarkerColor(4)
                  obj5.SetMarkerColor(5)
                  obj6.SetMarkerColor(6)
                  obj7.SetMarkerColor(7)
                  obj8.SetMarkerColor(8)
                  obj9.SetMarkerColor(9)
                  obj10.SetMarkerColor(10)
                  obj11.SetMarkerColor(11)
                  obj12.SetMarkerColor(12)
                  
                  obj1.SetLineColor(1)
                  obj2.SetLineColor(2)
                  obj3.SetLineColor(3)
                  obj4.SetLineColor(4)
                  obj5.SetLineColor(5)
                  obj6.SetLineColor(6)
                  obj7.SetLineColor(7)
                  obj8.SetLineColor(8)
                  obj9.SetLineColor(9)
                  obj10.SetLineColor(10)
                  obj11.SetLineColor(11)
                  obj12.SetLineColor(12)
                  
                  
                  obj1.SetMarkerStyle(20)
                  obj2.SetMarkerStyle(20)
                  obj3.SetMarkerStyle(20)
                  obj4.SetMarkerStyle(20)
                  obj5.SetMarkerStyle(20)
                  obj6.SetMarkerStyle(20)
                  obj7.SetMarkerStyle(20)
                  obj8.SetMarkerStyle(20)
                  obj9.SetMarkerStyle(20)
                  obj10.SetMarkerStyle(20)
                  obj11.SetMarkerStyle(20)
                  obj12.SetMarkerStyle(20)
                  
#                  obj1.SetTitle("")
#                  obj1.GetXaxis().SetTitleSize(0.05)
#                  obj1.GetXaxis().SetTitleOffset(1)
#                  obj1.GetYaxis().SetTitleSize(0.05)
#                  obj1.GetYaxis().SetTitleOffset(1)
                  obj12.SetStats(0)
                  max = np.maximum(obj12.GetMaximum(), obj1.GetMaximum())
                  obj12.SetMaximum(max*1.5)
                  
                  legMass =  ROOT.TLegend(.15,.55,.40,.9,"","brNDC")
                  legMass.SetTextFont(42)
                  legMass.SetTextSize(0.035)
                  legMass.SetBorderSize(1);
                  legMass.SetLineColor(1);
                  legMass.SetLineStyle(1);
                  legMass.SetLineWidth(1);
                  legMass.SetFillColor(0);
                  legMass.SetFillStyle(1001);
                  legMass.AddEntry(obj12,"Mass=2600 GeV","LP")
                  legMass.AddEntry(obj11,"Mass=2400 GeV","LP")
                  legMass.AddEntry(obj10,"Mass=2200 GeV","LP")
                  legMass.AddEntry(obj9,"Mass=2000 GeV","LP")
                  legMass.AddEntry(obj8,"Mass=1800 GeV","LP")
                  legMass.AddEntry(obj7,"Mass=1400 GeV","LP")
                  legMass.AddEntry(obj6,"Mass=1000 GeV","LP")
                  legMass.AddEntry(obj5,"Mass=800 GeV","LP")
                  legMass.AddEntry(obj4,"Mass=500 GeV","LP")
                  legMass.AddEntry(obj3,"Mass=400 GeV","LP")
                  legMass.AddEntry(obj2,"Mass=200 GeV","LP")
                  legMass.AddEntry(obj1,"Mass=100 GeV","LP")
                  
                  
                  obj12.Draw("COLZ L")
                  obj11.Draw("SAME")
                  obj10.Draw("SAME")
                  obj9.Draw("SAME")
                  obj8.Draw("SAME")
                  obj7.Draw("SAME")
                  obj6.Draw("SAME")
                  obj5.Draw("SAME")
                  obj4.Draw("SAME")
                  obj3.Draw("SAME")
                  obj2.Draw("SAME")
                  obj1.Draw("SAME")
                  
                  legMass.Draw("SAME")
              else :
                  continue
              can.SaveAs(name)
              if (keyname2 == "genlevelbeta") :
                  meanValue1 = obj1.GetMean()
                  meanValue2 = obj2.GetMean()
                  meanValue3 = obj3.GetMean()
                  meanValue4 = obj4.GetMean()
                  meanValue5 = obj5.GetMean()
                  meanValue6 = obj6.GetMean()
                  meanValue7 = obj7.GetMean()
                  meanValue8 = obj8.GetMean()
                  meanValue9 = obj9.GetMean()
                  meanValue10 = obj10.GetMean()
                  meanValue11 = obj11.GetMean()
                  meanValue12 = obj12.GetMean()
                  
#                  meanValue1 = obj1.GetBinCenter(obj1.GetMaximumBin())
#                  meanValue2 = obj2.GetBinCenter(obj2.GetMaximumBin())
#                  meanValue3 = obj3.GetBinCenter(obj3.GetMaximumBin())
#                  meanValue4 = obj4.GetBinCenter(obj4.GetMaximumBin())
#                  meanValue5 = obj5.GetBinCenter(obj5.GetMaximumBin())
#                  meanValue6 = obj6.GetBinCenter(obj6.GetMaximumBin())
#                  meanValue7 = obj7.GetBinCenter(obj7.GetMaximumBin())
#                  meanValue8 = obj8.GetBinCenter(obj8.GetMaximumBin())
#                  meanValue9 = obj9.GetBinCenter(obj9.GetMaximumBin())
#                  meanValue10 = obj10.GetBinCenter(obj10.GetMaximumBin())
#                  meanValue11 = obj11.GetBinCenter(obj11.GetMaximumBin())
#                  meanValue12 = obj12.GetBinCenter(obj12.GetMaximumBin())
                  
                  error1 = obj1.GetStdDev()
                  error2 = obj2.GetStdDev()
                  error3 = obj3.GetStdDev()
                  error4 = obj4.GetStdDev()
                  error5 = obj5.GetStdDev()
                  error6 = obj6.GetStdDev()
                  error7 = obj7.GetStdDev()
                  error8 = obj8.GetStdDev()
                  error9 = obj9.GetStdDev()
                  error10 = obj10.GetStdDev()
                  error11 = obj11.GetStdDev()
                  error12 = obj12.GetStdDev()
                  
                  betaGammaMPV1 = meanValue1 * (1/np.sqrt(1-meanValue1*meanValue1))
                  betaGammaMPV2 = meanValue2 * (1/np.sqrt(1-meanValue2*meanValue2))
                  betaGammaMPV3 = meanValue3 * (1/np.sqrt(1-meanValue3*meanValue3))
                  betaGammaMPV4 = meanValue4 * (1/np.sqrt(1-meanValue4*meanValue4))
                  betaGammaMPV5 = meanValue5 * (1/np.sqrt(1-meanValue5*meanValue5))
                  betaGammaMPV6 = meanValue6 * (1/np.sqrt(1-meanValue6*meanValue6))
                  betaGammaMPV7 = meanValue7 * (1/np.sqrt(1-meanValue7*meanValue7))
                  betaGammaMPV8 = meanValue8 * (1/np.sqrt(1-meanValue8*meanValue8))
                  betaGammaMPV9 = meanValue9 * (1/np.sqrt(1-meanValue9*meanValue9))
                  betaGammaMPV10 = meanValue10 * (1/np.sqrt(1-meanValue10*meanValue10))
                  betaGammaMPV11 = meanValue11 * (1/np.sqrt(1-meanValue11*meanValue11))
                  betaGammaMPV12 = meanValue12 * (1/np.sqrt(1-meanValue12*meanValue12))
                  
                  betaGammaMPV1Up = (meanValue1+error1) * (1/np.sqrt(1-(meanValue1+error1)*(meanValue1+error1)))
                  betaGammaMPV2Up = (meanValue2+error2) * (1/np.sqrt(1-(meanValue2+error2)*(meanValue2+error2)))
                  betaGammaMPV3Up = (meanValue3+error3) * (1/np.sqrt(1-(meanValue3+error3)*(meanValue3+error3)))
                  betaGammaMPV4Up = (meanValue4+error4) * (1/np.sqrt(1-(meanValue4+error4)*(meanValue4+error4)))
                  betaGammaMPV5Up = (meanValue5+error5) * (1/np.sqrt(1-(meanValue5+error5)*(meanValue5+error5)))
                  betaGammaMPV6Up = (meanValue6+error6) * (1/np.sqrt(1-(meanValue6+error6)*(meanValue6+error6)))
                  betaGammaMPV7Up = (meanValue7+error7) * (1/np.sqrt(1-(meanValue7+error7)*(meanValue7+error7)))
                  betaGammaMPV8Up = (meanValue8+error8) * (1/np.sqrt(1-(meanValue8+error8)*(meanValue8+error8)))
                  betaGammaMPV9Up = (meanValue9+error9) * (1/np.sqrt(1-(meanValue9+error9)*(meanValue9+error9)))
                  betaGammaMPV10Up = (meanValue10+error10) * (1/np.sqrt(1-(meanValue10+error10)*(meanValue10+error10)))
                  betaGammaMPV11Up = (meanValue11+error11) * (1/np.sqrt(1-(meanValue11+error11)*(meanValue11+error11)))
                  betaGammaMPV12Up = (meanValue12+error12) * (1/np.sqrt(1-(meanValue12+error12)*(meanValue12+error12)))
                  
#                  betaGammaErr1 = error1 * error1 / (np.sqrt((1-error1 * error1)**3))
                  betaGammaErr1 = betaGammaMPV1Up - betaGammaMPV1
                  betaGammaErr2 = betaGammaMPV2Up - betaGammaMPV2
                  betaGammaErr3 = betaGammaMPV3Up - betaGammaMPV3
                  betaGammaErr4 = betaGammaMPV4Up - betaGammaMPV4
                  betaGammaErr5 = betaGammaMPV5Up - betaGammaMPV5
                  betaGammaErr6 = betaGammaMPV6Up - betaGammaMPV6
                  betaGammaErr7 = betaGammaMPV7Up - betaGammaMPV7
                  betaGammaErr8 = betaGammaMPV8Up - betaGammaMPV8
                  betaGammaErr9 = betaGammaMPV9Up - betaGammaMPV9
                  betaGammaErr10 = betaGammaMPV10Up - betaGammaMPV10
                  betaGammaErr11 = betaGammaMPV11Up - betaGammaMPV11
                  betaGammaErr12 = betaGammaMPV12Up - betaGammaMPV12
                  
                  
                  fileOut.write("100\t"+str(betaGammaMPV1)+"\t"+str(betaGammaErr1)+"\n")
                  fileOut.write("200\t"+str(betaGammaMPV2)+"\t"+str(betaGammaErr2)+"\n")
                  fileOut.write("400\t"+str(betaGammaMPV3)+"\t"+str(betaGammaErr3)+"\n")
                  fileOut.write("500\t"+str(betaGammaMPV4)+"\t"+str(betaGammaErr4)+"\n")
                  fileOut.write("800\t"+str(betaGammaMPV5)+"\t"+str(betaGammaErr5)+"\n")
                  fileOut.write("1000\t"+str(betaGammaMPV6)+"\t"+str(betaGammaErr6)+"\n")
                  fileOut.write("1400\t"+str(betaGammaMPV7)+"\t"+str(betaGammaErr7)+"\n")
                  fileOut.write("1800\t"+str(betaGammaMPV8)+"\t"+str(betaGammaErr8)+"\n")
                  fileOut.write("2000\t"+str(betaGammaMPV9)+"\t"+str(betaGammaErr9)+"\n")
                  fileOut.write("2200\t"+str(betaGammaMPV10)+"\t"+str(betaGammaErr10)+"\n")
                  fileOut.write("2400\t"+str(betaGammaMPV11)+"\t"+str(betaGammaErr11)+"\n")
                  fileOut.write("2600\t"+str(betaGammaMPV12)+"\t"+str(betaGammaErr12)+"\n")
                  fileOut.close()
              #can.SaveAs(name.replace(".png",".pdf"))
              #can.SaveAs(name.replace(".png",".C"))
              can.Close()

#Mass = f.Get("/analyzer/BaseName/Mass")
#Mass_wPred = f.Get("/analyzer/BaseName/Pred_Mass_CB")
#if Mass_wPred :
#  massBins = [0.,50.,100.,150.,200.,250.,300.,350.,400.,500.,600.,700.,800.,900.,1000.,4000.]
#  massBinsArray = np.array(massBins)
#  Mass_projY_NotRebinned = Mass.ProjectionY("Mass_projY",ProjBin,ProjBin,"e")
#  Mass_projY = ROOT.TH1F("Mass_projY" , "Mass_projY" , len(massBinsArray)-1, massBinsArray)
#
#  Mass_wPred_projY_NotRebinned = Mass_wPred.ProjectionY("Mass_wPred_projY",ProjBin,ProjBin,"e")
#  Mass_wPred_projY = ROOT.TH1F("Mass_wPred_projY" , "Mass_wPred_projY" , len(massBinsArray)-1, massBinsArray)
#
#  print("Mass_projY_NotRebinned.Integral(): ",Mass_projY_NotRebinned.Integral())
#  print("Mass_wPred_projY_NotRebinned.Integral(): ",Mass_wPred_projY_NotRebinned.Integral())
#
#  KSvalue = Mass_projY_NotRebinned.KolmogorovTest(Mass_wPred_projY_NotRebinned,"XD")
#  print("KS-test: "+str(KSvalue))
#
#  for i, bin in enumerate(massBinsArray):
#    Mass_projYCont = Mass_projY_NotRebinned.GetBinContent(i+1)
#    Mass_projY.SetBinContent(i+1,Mass_projYCont)
#    Mass_projYCont_err = Mass_projY_NotRebinned.GetBinError(i+1)
#    #Mass_projYCont_err = Mass.GetBinErrorLow(ProjBin,i+1)
#    Mass_projY.SetBinError(i+1,Mass_projYCont_err)
#
#    Mass_wPred_projYCont = Mass_wPred_projY_NotRebinned.GetBinContent(i+1)
#    Mass_wPred_projY.SetBinContent(i+1,Mass_wPred_projYCont)
#    Mass_wPred_projYCont_err = Mass_wPred_projY_NotRebinned.GetBinError(i+1)
#    #Mass_wPred_projYCont_err = Mass_wPred.GetBinErrorLow(ProjBin,i+1)
#    Mass_wPred_projY.SetBinError(i+1,Mass_wPred_projYCont_err)
#
#  print("----------------------------------------------")
#  KSvalue2 = Mass_projY.KolmogorovTest(Mass_wPred_projY,"XD")
#  print("KS-test after rebinning: "+str(KSvalue2))
#
#  Mass_projY.SetMarkerColor(1)
#  Mass_projY.SetLineColor(1)
#  Mass_projY.SetMarkerStyle(20)
#  Mass_projY.SetTitle("")
#  Mass_projY.GetXaxis().SetTitleSize(0.05)
#  Mass_projY.GetXaxis().SetTitleOffset(1)
#  Mass_projY.GetXaxis().SetTitle("Mass [GeV]")
#  Mass_projY.GetYaxis().SetTitle("Tracks/bin")
#  Mass_projY.GetYaxis().SetTitleSize(0.05)
#  Mass_projY.GetYaxis().SetTitleOffset(1)
#  Mass_projY.SetStats(0)
#  Mass_projY.GetYaxis().SetRangeUser(0.1,10000)
#
#
#  Mass_wPred_projY.SetMarkerColor(2)
#  Mass_wPred_projY.SetLineColor(2)
#  Mass_wPred_projY.SetMarkerStyle(20)
#  Mass_wPred_projY.SetTitle("")
#  Mass_wPred_projY.GetXaxis().SetTitleSize(0.05)
#  Mass_wPred_projY.GetXaxis().SetTitleOffset(1)
#  Mass_wPred_projY.GetXaxis().SetTitle("Mass [GeV]")
#  Mass_wPred_projY.GetYaxis().SetTitle("Tracks/bin")
#  Mass_wPred_projY.GetYaxis().SetTitleSize(0.05)
#  Mass_wPred_projY.GetYaxis().SetTitleOffset(1)
#  Mass_wPred_projY.SetStats(0)
#
#
#  print("Mass_projY.Integral(): ",Mass_projY.Integral())
#  print("Mass_wPred_projY.Integral(): ",Mass_wPred_projY.Integral())
#
#  legMass =  ROOT.TLegend(.45,.75,.80,.9,"","brNDC")
#  legMass.SetTextFont(42)
#  legMass.SetTextSize(0.035)
#  legMass.SetBorderSize(1);
#  legMass.SetLineColor(1);
#  legMass.SetLineStyle(1);
#  legMass.SetLineWidth(1);
#  legMass.SetFillColor(0);
#  legMass.SetFillStyle(1001);
#  legMass.AddEntry(Mass_wPred_projY,"Prediction","LP")
#  legMass.AddEntry(Mass_projY,"Observation","LP")
#
#  tex2 = ROOT.TLatex(0.13,0.94,"CMS");
#  #tex2 = ROOT.TLatex(0.20,0.94,"CMS");#if there is 10^x
#  tex2.SetNDC();
#  tex2.SetTextFont(61);
#  tex2.SetTextSize(0.0675);
#  tex2.SetLineWidth(2);
#
#  #tex3 = ROOT.TLatex(0.27,0.96,"Simulation"); # for square plots
#  #tex3 = ROOT.TLatex(0.28,0.94,"Work in Progress 2018"); #if there is 10^x
#  tex3 = ROOT.TLatex(0.24,0.94,"Internal");
#  tex3.SetNDC();
#  tex3.SetTextFont(52);
#  tex3.SetTextSize(0.0485);
#  tex3.SetLineWidth(2);
#
#  tex4 = ROOT.TLatex(0.7,0.93,"K-S test v2: "+str(round(KSvalue2,4)));
#  tex4.SetNDC();
#  tex4.SetTextFont(52);
#  tex4.SetTextSize(0.0485);
#  tex4.SetLineWidth(2);
#
#
#
#
#
#  cMass_projY = ROOT.TCanvas('cMass_projY', 'cMass_projY',800,800)
#  cMass_projY.SetLogy()
#
#  rp = ROOT.TRatioPlot(Mass_projY,Mass_wPred_projY)
#
#  rp.SetH1DrawOpt("P");
#  rp.SetH2DrawOpt("P");
#
#  rp.Draw()
#  #rp.GetUpperPad().BuildLegend()
#  rp.SetLeftMargin(0.13);
#  rp.SetRightMargin(0.05);
#  rp.SetUpTopMargin(0.1);
#  rp.SetLowTopMargin(0.02);
#  rp.SetLowBottomMargin(0.35);
#
#  rp.GetLowerRefGraph().SetMinimum(0.01);
#  rp.GetLowerRefGraph().SetMaximum(2);
#  #rp.GetLowerRefGraph().SetMarkerColor(ROOT.kGreen+2)
#  #rp.GetLowerRefGraph().SetLineColor(0) #0
#  rp.GetLowerRefGraph().SetMarkerStyle(20)
#  rp.GetLowerRefGraph().SetMarkerSize(1);
#  rp.GetLowYaxis().SetNdivisions(505);
#  rp.GetLowerRefYaxis().SetTitle("Ratio");
#  rp.GetLowerRefYaxis().SetTitleSize(0.05);
#  rp.GetLowerRefYaxis().SetTitleOffset(1);
#  rp.GetLowerRefYaxis().SetLabelSize(0.035);
#
#
#  rp.GetLowerRefXaxis().SetTitleSize(0.05);
#  rp.GetLowerRefXaxis().SetTitleOffset(0.8);
#  rp.GetLowerRefXaxis().SetLabelSize(0.035);
#  cMass_projY.Modified()
#  cMass_projY.Update()
#  #Mass_projY.Draw()
#  #Mass_wPred_projY.Draw("SAME")
#  #rp.Draw("X")
#
#  rp.GetUpperPad().cd();
#  legMass.Draw("SAME")
#  tex2.Draw("SAME")
#  tex3.Draw("SAME")
#  tex4.Draw("SAME")
#
#  name = newFileDir + "/cMass_new.png"
#  cMass_projY.SaveAs(name)
#
#os.system("cp forWebpage/* "+newFileDir+"/.")
#print("scp -r "+ newFileDir + " tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV8/.")
