import ROOT, sys, os, time, re
import numpy as np
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog fileName.root")
(opt,args) = parser.parse_args()

ROOT.gROOT.SetStyle("Plain")
ROOT.gStyle.SetPalette(1)
ROOT.gROOT.SetBatch(True)

fileName = sys.argv[1]

print("Filename: "+fileName)

f1 = ROOT.TFile.Open(fileName)

for i in range(0, f1.GetListOfKeys().GetEntries()):
  dirname = f1.GetListOfKeys().At(i).GetName()
  curr_dir = f1.GetDirectory(dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = f1.GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          newname = dirname + "/" + keyname+ "/" + keyname2
          obj1 = f1.Get(newname)
          if obj1.InheritsFrom("TObject"):
              can = obj1
              can = ROOT.TCanvas(newname)
              name = keyname2 +  ".png"
              if (keyname2 == "CutFlow") :
                obj1.SetMarkerColor(1)
                obj1.SetLineColor(1)
                obj1.SetMarkerStyle(20)
                obj1.SetStats(0)

                obj1.Scale(1/obj1.GetMaximum())
              
                legend =  ROOT.TLegend(.35,.85,.90,.90,"","brNDC")
                legend.SetTextFont(42)
                legend.SetTextSize(0.03)
                legend.SetBorderSize(1);
                legend.SetLineColor(1);
                legend.SetLineStyle(1);
                legend.SetLineWidth(1);
                legend.SetFillColor(0);
                legend.SetFillStyle(1001);
                legend.AddEntry(obj1,fileName,"LP")
                    
                obj1.Draw("COLZ L")
                    
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
                obj1.GetXaxis().SetBinLabel(16,"MiniIso")
                obj1.GetXaxis().SetBinLabel(17,"MassT")
                obj1.GetXaxis().SetBinLabel(18,"Ih")
                obj1.GetXaxis().SetBinLabel(19,"ProbQ")
                obj1.GetXaxis().SetBinLabel(20,"MuStat")
                obj1.GetXaxis().SetBinLabel(21,"PhiTOF")
                obj1.GetXaxis().SetBinLabel(22,"EtaTOF")
                    
                legend.Draw("SAME")
                can.SaveAs("CutFlow_"+fileName[:-5]+".png")
                #can.SaveAs(name.replace(".png",".pdf"))
                #can.SaveAs(name.replace(".png",".C"))
                can.Close()
            
              if (keyname2 == "N1Eta"):
                N1EtaEff = obj1.Integral(obj1.GetXaxis().FindBin(-2.1),obj1.GetXaxis().FindBin(2.1))/(obj1.Integral(obj1.GetXaxis().FindBin(-3),obj1.GetXaxis().FindBin(3)))
              elif (keyname2 == "N1Chi2PerNdof"):
                N1Chi2PerNdofEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(5))/(obj1.Integral())
              elif (keyname2 == "N1Dxy"):
                N1DxyEff = obj1.Integral(obj1.GetXaxis().FindBin(-0.5),obj1.GetXaxis().FindBin(0.5))/(obj1.Integral()) #obj1.GetXaxis().FindBin(0),obj1.GetXaxis().FindBin(10))
              elif (keyname2 == "N1Dz"):
                N1DzEff = obj1.Integral(obj1.GetXaxis().FindBin(-0.5),obj1.GetXaxis().FindBin(0.5))/(obj1.Integral()) #obj1.GetXaxis().FindBin(0),obj1.GetXaxis().FindBin(10))
              elif (keyname2 == "N1EIsol"):
                N1EIsolEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(0.3))/(obj1.Integral())
              elif (keyname2 == "N1MIh"):
                N1MIhEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(40))/(obj1.Integral())
              elif (keyname2 == "N1MPt"):
                N1MPtEff = obj1.Integral(obj1.GetXaxis().FindBin(55.0),obj1.GetXaxis().FindBin(4001.0))/(obj1.Integral())
              elif (keyname2 == "N1ProbQ"):
                N1ProbQEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(0.1))/(obj1.Integral())
              elif (keyname2 == "N1ProbXY"):
                N1ProbXYEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(1.0))/(obj1.Integral())
              elif (keyname2 == "N1PtErrOverPt"):
                N1PtErrOverPtEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(0.25))/(obj1.Integral())
              elif (keyname2 == "N1Qual"):
                N1QualEff = obj1.Integral(obj1.GetXaxis().FindBin(1),obj1.GetXaxis().FindBin(2.1))/(obj1.Integral())
#              elif (keyname2 == "N1SegSep"):
#                N1SegSepEff = obj1.Integral(obj1.GetXaxis().FindBin(-2.1),obj1.GetXaxis().FindBin(2.1))/(obj1.Integral())
              elif (keyname2 == "N1Stations"):
                N1StationsEff = obj1.Integral(obj1.GetXaxis().FindBin(2),obj1.GetXaxis().FindBin(8))/(obj1.Integral())
              elif (keyname2 == "N1TNOH"):
                N1TNOHEff = obj1.Integral(obj1.GetXaxis().FindBin(9),obj1.GetXaxis().FindBin(40))/(obj1.Integral())
              elif (keyname2 == "N1TNOHFraction"):
                N1TNOHFractionEff = obj1.Integral(obj1.GetXaxis().FindBin(0.8),obj1.GetXaxis().FindBin(0.999))/(obj1.Integral())
              elif (keyname2 == "N1TNOPH"):
                N1TNOPHEff = obj1.Integral(obj1.GetXaxis().FindBin(3),obj1.GetXaxis().FindBin(8))/(obj1.Integral())
                
              elif (keyname2 == "BS_Eta"):
                BSEtaEff = obj1.Integral(obj1.GetXaxis().FindBin(-2.1),obj1.GetXaxis().FindBin(2.1))/(obj1.Integral(obj1.GetXaxis().FindBin(-3),obj1.GetXaxis().FindBin(3)))
              elif (keyname2 == "BS_Chi2PerNdof"):
                BS_Chi2PerNdofEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(5))/(obj1.Integral())
              elif (keyname2 == "BS_dxyMinv3d"):
                BS_DxyEff = obj1.Integral(obj1.GetXaxis().FindBin(-0.5),obj1.GetXaxis().FindBin(0.5))/(obj1.Integral()) #obj1.GetXaxis().FindBin(0),obj1.GetXaxis().FindBin(10))
              elif (keyname2 == "BS_dzMinv3d"):
                BS_DzEff = obj1.Integral(obj1.GetXaxis().FindBin(-0.5),obj1.GetXaxis().FindBin(0.5))/(obj1.Integral()) #obj1.GetXaxis().FindBin(0),obj1.GetXaxis().FindBin(10))
              elif (keyname2 == "BS_EIsol"):
                BS_EIsolEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(0.3))/(obj1.Integral())
              elif (keyname2 == "BS_MIh"):
                BS_MIhEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(40))/(obj1.Integral())
              elif (keyname2 == "BS_MPt"):
                BS_MPtEff = obj1.Integral(obj1.GetXaxis().FindBin(55),obj1.GetXaxis().FindBin(3999))/(obj1.Integral())
              elif (keyname2 == "BS_ProbQ"):
                BS_ProbQEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(0.1))/(obj1.Integral())
              elif (keyname2 == "BS_ProbXY"):
                BS_ProbXYEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(1.0))/(obj1.Integral())
              elif (keyname2 == "BS_PtErrOverPt"):
                BS_PtErrOverPtEff = obj1.Integral(obj1.GetXaxis().FindBin(0.0),obj1.GetXaxis().FindBin(0.25))/(obj1.Integral())
              elif (keyname2 == "BS_Qual"):
                BS_QualEff = obj1.Integral(obj1.GetXaxis().FindBin(1.1),obj1.GetXaxis().FindBin(2.1))/(obj1.Integral())
#              elif (keyname2 == "BS_SegSep"):
#                BS_SegSepEff = obj1.Integral(obj1.GetXaxis().FindBin(-2.1),obj1.GetXaxis().FindBin(2.1))/(obj1.Integral())
#              elif (keyname2 == "BS_Stations"):
#                BS_StationsEff = obj1.Integral(obj1.GetXaxis().FindBin(2),obj1.GetXaxis().FindBin(8))/(obj1.Integral())
              elif (keyname2 == "BS_TNOH"):
                BS_TNOHEff = obj1.Integral(obj1.GetXaxis().FindBin(9),obj1.GetXaxis().FindBin(40))/(obj1.Integral())
              elif (keyname2 == "BS_TNOHFraction"):
                BS_TNOHFractionEff = obj1.Integral(obj1.GetXaxis().FindBin(0.8),obj1.GetXaxis().FindBin(0.999))/(obj1.Integral())
              elif (keyname2 == "BS_TNOPH"):
                BS_TNOPHEff = obj1.Integral(obj1.GetXaxis().FindBin(3),obj1.GetXaxis().FindBin(8))/(obj1.Integral())



fileOutEtaEff = open("EtaEff.txt", "a")
fileOutMPtEff = open("MPtEff.txt", "a")
fileOutTNOHEff = open("TNOHEff.txt", "a")
fileOutTNOPHEff = open("TNOPHEff.txt", "a")
fileOutTNOHFractionEff = open("TNOHFractionEff.txt", "a")
fileOutProbQEff = open("ProbQEff.txt", "a")
fileOutProbXYEff = open("ProbXYEff.txt", "a")
fileOutChi2PerNdofEff = open("Chi2PerNdofEff.txt", "a")
fileOutEIsolEff = open("EIsolEff.txt", "a")
fileOutMIhEff = open("MIhEff.txt", "a")
fileOutPtErrOverPtEff = open("PtErrOverPtEff.txt", "a")
fileOutDzEff = open("DzEff.txt", "a")
fileOutDxyEff = open("DxyEff.txt", "a")

idx = fileName.find("_M_")
fileOutEtaEff.write(fileName[idx+3:idx+7]+"\t"+str(BSEtaEff)+"\n")
fileOutMPtEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_MPtEff)+"\n")
fileOutTNOHEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_TNOHEff)+"\n")
fileOutTNOPHEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_TNOPHEff)+"\n")
fileOutTNOHFractionEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_TNOHFractionEff)+"\n")
fileOutProbQEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_ProbQEff)+"\n")
fileOutProbXYEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_ProbXYEff)+"\n")
fileOutChi2PerNdofEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_Chi2PerNdofEff)+"\n")
fileOutEIsolEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_EIsolEff)+"\n")
fileOutMIhEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_MIhEff)+"\n")
fileOutPtErrOverPtEff .write(fileName[idx+3:idx+7]+"\t"+str(BS_PtErrOverPtEff)+"\n")
fileOutDzEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_DzEff)+"\n")
fileOutDxyEff.write(fileName[idx+3:idx+7]+"\t"+str(BS_DxyEff)+"\n")

can = ROOT.TCanvas("CutEffs")
obj1 = ROOT.TH1F("CutEffs","CutEffs",14,0.0,14.0)
obj2 = ROOT.TH1F("CutEffsN1","CutEffsN1",14,0.0,14.0)
obj1.SetBinContent(1,1)
obj1.SetBinContent(2,BSEtaEff)
obj1.SetBinContent(3,BS_MPtEff)
obj1.SetBinContent(4,BS_TNOHEff)
obj1.SetBinContent(5,BS_TNOPHEff)
obj1.SetBinContent(6,BS_TNOHFractionEff)
obj1.SetBinContent(7,BS_ProbXYEff)
obj1.SetBinContent(8,BS_Chi2PerNdofEff)
obj1.SetBinContent(9,BS_EIsolEff)
obj1.SetBinContent(10,BS_DzEff)
obj1.SetBinContent(11,BS_DxyEff)
obj1.SetBinContent(12,BS_PtErrOverPtEff)
obj1.SetBinContent(13,BS_MIhEff)
obj1.SetBinContent(14,BS_ProbQEff)
obj1.SetMaximum(1.1)
#obj1.SetBinContent(9,1) #N1QualEff/BS_QualEff
#obj1.SetBinContent(16,1)
#obj1.SetBinContent(17,1)
#obj1.SetBinContent(18,1)
#obj1.SetBinContent(19,1)
#obj1.SetBinContent(20,1)

obj2.SetBinContent(1,1)
obj2.SetBinContent(2,N1EtaEff)
obj2.SetBinContent(3,N1MPtEff)
obj2.SetBinContent(4,N1TNOHEff)
obj2.SetBinContent(5,N1TNOPHEff)
obj2.SetBinContent(6,N1TNOHFractionEff)
obj2.SetBinContent(7,N1ProbXYEff)
obj2.SetBinContent(8,N1Chi2PerNdofEff)
obj2.SetBinContent(9,N1EIsolEff)
obj2.SetBinContent(10,N1DzEff)
obj2.SetBinContent(11,N1DxyEff)
obj2.SetBinContent(12,N1PtErrOverPtEff)
obj2.SetBinContent(13,N1MIhEff)
obj2.SetBinContent(14,N1ProbQEff)
#obj2.SetBinContent(9,1) #N1QualEff/BS_QualEff
#obj2.SetBinContent(16,1)
#obj2.SetBinContent(17,1)
#obj2.SetBinContent(18,1)
#obj2.SetBinContent(19,1)
#obj2.SetBinContent(20,1)
obj2.SetLineColor(2)
obj2.SetMarkerColor(2)

obj1.GetXaxis().SetBinLabel(1,"Trigger")
obj1.GetXaxis().SetBinLabel(2,"Eta")
obj1.GetXaxis().SetBinLabel(3,"pT")
obj1.GetXaxis().SetBinLabel(4,"NumHits")
obj1.GetXaxis().SetBinLabel(5,"NumPixHits")
obj1.GetXaxis().SetBinLabel(6,"ValidFract")
obj1.GetXaxis().SetBinLabel(7,"ProbXY")
obj1.GetXaxis().SetBinLabel(8,"Chi2oDOF")
obj1.GetXaxis().SetBinLabel(9,"EoP")
obj1.GetXaxis().SetBinLabel(10,"dz")
obj1.GetXaxis().SetBinLabel(11,"dxy")
obj1.GetXaxis().SetBinLabel(12,"pTerrOverpT")
obj1.GetXaxis().SetBinLabel(13,"Ih")
obj1.GetXaxis().SetBinLabel(14,"ProbQ")
#obj1.GetXaxis().SetBinLabel(15,"TKIso")
#obj1.GetXaxis().SetBinLabel(9,"HighPurity")
#obj1.GetXaxis().SetBinLabel(17,"MuStat")
#obj1.GetXaxis().SetBinLabel(18,"EtaTOF")
#obj1.GetXaxis().SetBinLabel(19,"PhiTOF")
#obj1.GetXaxis().SetBinLabel(20,"N/A")
obj1.GetYaxis().SetTitle("Efficiency")
obj1.GetYaxis().SetTitleOffset(1.1)
obj1.SetStats(0)

legend =  ROOT.TLegend(.65,.75,.90,.9,"","brNDC")
legend.SetTextFont(42)
legend.SetTextSize(0.035)
legend.SetBorderSize(1);
legend.SetLineColor(1);
legend.SetLineStyle(1);
legend.SetLineWidth(1);
legend.SetFillColor(0);
legend.SetFillStyle(1001);
legend.AddEntry(obj1,"Cut alone","LP")
legend.AddEntry(obj2,"Cut with other cuts","LP")

obj1.Draw()
obj2.Draw("SAME")
legend.Draw("SAME")
can.SaveAs("CutEffs_"+fileName[:-5]+".png")
