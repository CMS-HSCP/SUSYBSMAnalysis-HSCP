import ROOT, sys, os, time, re, numpy, random
from optparse import OptionParser
from ctypes import c_double as double
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

codeVersion = sys.argv[1]

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

InputListSamples = [
"HSCPgluinoSamplesSyst_CodeV"+codeVersion+".txt"
]

OpenedSamples = []
for list in InputListSamples:
  MassPointSamples = []
  with open(list, "r") as a_file:
    for line in a_file:
      stripped_line = line.strip()
      MassPointSamples.append(ROOT.TFile.Open(stripped_line))
    OpenedSamples.append(MassPointSamples)

if os.path.exists(os.path.dirname("ROC_CodeV"+codeVersion)):
  os.makedirs("ROC_CodeV"+codeVersion)
  
cutValues = {
'N1_Eta' : 1.,
'N1_Pt' : 55.,
'N1_Pt_lowPt' : 55.,
'N1_Dxy' : 0.02,
'N1_Dz' : 0.1,
'N1_dRMinPfJet' : 1000.,
'N1_Chi2oNdof' : 5.,
'N1_Qual' : 2.,
'N1_TNOH' : 10.,
'N1_TNOM' : 6.,
'N1_Qual' : 2.,
'N1_TNOPH' : 2.,
'N1_TNOHFraction' : .8,
'N1_EoP' : .3,
'N1_SumpTOverpT' : 1000.,
'N1_Ih' : 3.2,
'N1_ProbQ' : 0.85,
'PostPreS_ProbQNoL1' : 0.85,
'N1_Stations' : 9.,
'N1_PtErrOverPt' : .25,
'N1_SegSep' : 25,
'N1_ProbXY' : 1.0,
'N1_MiniRelIsoAll' : 0.02,
'N1_MiniRelTkIso' : 0.02,
'N1_MiniRelIsoAll_lowMiniRelIso' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso_PUA' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso_PUB' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso_PUC' : 0.02,
'N1_MiniTkIso' : 15.,
'N1_MiniTkIso_PUA' : 15.,
'N1_MiniTkIso_PUB' : 15.,
'N1_MiniTkIso_PUC' : 15.,
'N1_MiniTkIso_lowMiniRelIso' : 15.,
'N1_MiniTkIso_lowMiniRelIso_PUA' : 15.,
'N1_MiniTkIso_lowMiniRelIso_PUB' : 15.,
'N1_MiniTkIso_lowMiniRelIso_PUC' : 15.,
'N1_MTOF' : 1000.,
'N1_I' : 1000.,
'N1_PtErrOverPt2' : 0.0007,
'PostPreS_Ias' : 0.3,
'PostPreS_Ias_NoEventWeight' : 0.3,
'PostPreS_IasForStatus91' : 0.3,
'PostPreS_IasForStatusNot91' : 0.3,
'N1_DrMinPfJet': 10000,
'N1_PfType': 10000,
'N1_ProbQNoLayer1': 0.85,

}


for i in range(0, MassPointSamples[0].GetListOfKeys().GetEntries()):
  dirname = MassPointSamples[0].GetListOfKeys().At(i).GetName()
  curr_dir = MassPointSamples[0].GetDirectory(dirname)
  if not (curr_dir) :
    continue
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = MassPointSamples[0].GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
#          if not ("N1_Eta" in keyname2) : continue
          if not ("N1" in keyname2 or "PostPreS_ProbQNoL1" in keyname2 or "PostPreS_Ias"  in keyname2): continue
          if ("Vs" in keyname2 or "N1_pfType" in keyname2) : continue
          cstackedSummedBackgroundString = "cstackedSummedBackgroundString"+str(j)
          canvas = ROOT.TCanvas(cstackedSummedBackgroundString, cstackedSummedBackgroundString, 800,800)
          newname = dirname + "/" + keyname+ "/" + keyname2
          
          
          legend =  ROOT.TLegend(.55,.70,.80,.89,"","brNDC")
          legend.SetTextFont(42)
          legend.SetTextSize(0.02)
          legend.SetBorderSize(1);
          legend.SetBorderSize(0);
          legend.SetLineColor(1);
          legend.SetLineStyle(1);
          legend.SetLineWidth(1);
          legend.SetFillColor(0);
          legend.SetFillStyle(1001);
        
          # array to contain a specific (keyname2) histogram for all samples
          histoArray = []
          maxNumBins = 0.0
          for i in range(0,len(OpenedSamples)) :
            maxNumBins = numpy.maximum(maxNumBins,len(OpenedSamples[i]))
 
          for indexOnSamples in range(0,len(OpenedSamples)) :
            singleROCcurveString = 'singleROCcurve'+str(indexOnSamples)+str(j)
            if ("N1" in keyname2) : axisXTitle = keyname2[keyname2.find("N1_")+3:]
            if ("PostPreS" in keyname2) : axisXTitle = keyname2[keyname2.find("PostPreS_")+9:]
            singleROCcurve = ROOT.TH1F(singleROCcurveString,";;Efficiency for "+axisXTitle,int(maxNumBins),0.,maxNumBins)
            singleROCcurve.SetStats(0)
            singleROCcurve.SetMarkerColor(indexOnSamples+1)
            singleROCcurve.SetMarkerStyle(20)
            singleROCcurve.SetLineColor(indexOnSamples+1)
#            singleROCcurve.SetMaximum(1.5)
            
            for indexOnBins in range(0,len(OpenedSamples[indexOnSamples])) :
              givenSampleWmassPoint = OpenedSamples[indexOnSamples][indexOnBins]
              givenSampleWmassPointStr = str(givenSampleWmassPoint)
              stringIndexWhere2018 = givenSampleWmassPointStr.find("2018_")+5
              stringIndexWhereM = givenSampleWmassPointStr.find("_M")
              stringIndexWhereW = givenSampleWmassPointStr.find("_w")
              if (givenSampleWmassPointStr[stringIndexWhere2018:stringIndexWhereW].find("_M")>0) :
                sampleName = (givenSampleWmassPointStr[stringIndexWhere2018:stringIndexWhereM])
              else:
                sampleName = (givenSampleWmassPointStr[stringIndexWhere2018:stringIndexWhereW])
              histo = givenSampleWmassPoint.Get(newname)
              
              
              if not (cutValues.get(keyname2)) :
                print("No cut value for " + str(keyname2))
                continue
              DenomBkgErr = NumErr = NumBkgErr = DenomErr = double(0.0)

              Num = histo.IntegralAndError(histo.GetXaxis().FindBin(0),histo.GetXaxis().FindBin(cutValues[keyname2]),NumErr,"width")
              Denom = histo.IntegralAndError(histo.GetXaxis().FindBin(0),histo.GetXaxis().FindBin(histo.GetXaxis().GetXmax()),DenomErr,"width")
              
              if (Denom>0) :
                Eff = Num / Denom
              else : Eff = -1
              
              if (Denom>0 and Num>0) :
                EffErr =  Eff * numpy.sqrt( ((NumErr.value / Num) * (NumErr.value / Num)) + ((DenomErr.value / Denom) * (DenomErr.value / Denom) ))
              else :
                EffErr = 0.
      
              
              if ("N1_Pt" in keyname2) :
                Eff = 1-Eff
              if ("N1_Pt_lowPt" in keyname2) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
              if ("N1_TNOPH" in keyname2) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
              if ("N1_TNOHFraction" in keyname2) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
              if ("N1_TNOM" in keyname2) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
              if ("N1_Ih" in keyname2) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
              if ("N1_ProbXY" in keyname2) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
              if ("PostPreS_Ias" in keyname2) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
              if ("ProbQNoL1" in keyname2 or "ProbQNoLayer1" in keyname2 ) :
                Eff = 1-Eff
#                EffBkg = 1-EffBkg
                
          
#              print("For j-th: "+str(indexOnBins)+" Eff: "+str(Eff))
              singleROCcurve.SetBinContent(indexOnBins+1,Eff)
              singleROCcurve.SetBinError(indexOnBins+1,EffErr)
#              singleROCcurve.SetBinContent(8,EffBkg)
            histoArray.append(singleROCcurve)
            legend.AddEntry(singleROCcurve,sampleName,"LP")
     
          tex2 = ROOT.TLatex(0.13,0.94,"CMS");
          tex2.SetNDC();
          tex2.SetTextFont(61);
          tex2.SetTextSize(0.0675);
          tex2.SetLineWidth(2);

          tex3 = ROOT.TLatex(0.27,0.94,"Internal"); # for square plots
          tex3.SetNDC();
          tex3.SetTextFont(52);
          tex3.SetTextSize(0.0485);
          tex3.SetLineWidth(2);
          
          tex5 = ROOT.TLatex(0.07,0.04,"Code version: "+codeVersion);
          tex5.SetNDC();
          tex5.SetTextFont(52);
          tex5.SetTextSize(0.0185);
          tex5.SetLineWidth(2);
          
          for finalHisto in histoArray :
            finalHisto.Draw("SAMEHISTOTEXT00")
            finalHisto.GetXaxis().SetBinLabel(1,"Syst down")
            finalHisto.GetXaxis().SetBinLabel(2,"Syst nominal")
            finalHisto.GetXaxis().SetBinLabel(3,"Syst up")
            finalHisto.GetYaxis().SetRangeUser(0.,1.4)
          legend.Draw("SAME")
          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex5.Draw("SAME")

          canvas.SaveAs("ROC_CodeV"+codeVersion+"/"+keyname2+".png")
