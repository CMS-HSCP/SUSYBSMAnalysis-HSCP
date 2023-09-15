import ROOT, sys, os, time, re, numpy, random
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
(opt,args) = parser.parse_args()

codeVersion = sys.argv[1]

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

ROOT.gStyle.SetPaintTextFormat("4.2f");

SingleMuonSample = ROOT.TFile.Open("crab_Analysis_SingleMuon_Run2018C_CodeV"+codeVersion+"_v1.root")
AllBcgSample = ROOT.TFile.Open("crab_Analysis_2018_AllBackground_CodeV"+codeVersion+"_v1.root")
Rhadron1800GeV = ROOT.TFile.Open("crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root")

cutValues = {
'BefPreS_Pt' : 55.,
'N1_Pt' : 55.,
'N1_Pt_lowPt' : 55.,
'BefPreS_Eta' : 1.0,
'N1_Eta' : 1.0,
'BefPreS_TNOPH' : 1., # after 29p2 should be 1
'N1_TNOPH' : 1., # after 29p2 should be 1
'BefPreS_TNOHFraction' : .8,
'N1_TNOHFraction' : .8,
'BefPreS_TNOM' : 9.,
'N1_TNOM' : 9.,
'BefPreS_Chi2oNdof' : 5.,
'N1_Chi2oNdof' : 5.,
'BefPreS_Dz' : 0.1,
'N1_Dz' : 0.1,
'BefPreS_Dxy' : 0.02,
'N1_Dxy' : 0.02,
'BefPreS_PtErrOverPt2' : 0.001,
'PostPreS_PtErrOverPt2' : 0.001,
'N1_EoP' : 0.3,
'BefPreS_EoP' : 0.3,
'BefPreS_Ih' : 3.47,
'N1_Ih' : 3.47,
'BefPreS_ProbXY' : 0.01,
'N1_ProbXY' : 0.01,
'BefPreS_MiniRelIsoAll_lowMiniRelIso' : 0.02,
'N1_MiniRelIsoAll_lowMiniRelIso' : 0.02,
'BefPreS_MiniRelTkIso_lowMiniRelIso' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso' : 0.02,
'BefPreS_MiniRelTkIso_lowMiniRelIso_PUA' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso_PUA' : 0.02,
'BefPreS_MiniRelTkIso_lowMiniRelIso_PUB' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso_PUB' : 0.02,
'BefPreS_MiniRelTkIso_lowMiniRelIso_PUC' : 0.02,
'N1_MiniRelTkIso_lowMiniRelIso_PUC' : 0.02,
'BefPreS_P' : 3000,
'PostPreS_P' : 3000,
'BefPreS_TIsol' : 15,
'PostPreS_TIsol' : 15,
'BefPreS_MiniTkIso' : 15,
'N1_MiniTkIso' : 15,
'BefPreS_MiniTkIso_PUA' : 15,
'N1_MiniTkIso_PUA' : 15,
'BefPreS_MiniTkIso_PUB' : 15,
'N1_MiniTkIso_PUB' : 15,
'BefPreS_MiniTkIso_PUC' : 15,
'N1_MiniTkIso_PUC' : 15,
}


sigma = 3

if not os.path.exists(os.path.dirname("ComparePunzi_CodeV"+codeVersion+"/a.png")):
  print("Create dir")
  os.makedirs(os.path.dirname("ComparePunzi_CodeV"+codeVersion+"/"))

for i in range(0, SingleMuonSample.GetListOfKeys().GetEntries()):
  dirname = SingleMuonSample.GetListOfKeys().At(i).GetName()
  curr_dir = SingleMuonSample.GetDirectory(dirname)
  if not (curr_dir) :
    continue
  N1eff = N1effForP = N1effForPPunzi = N1effForPBcg = N1effForPSignal = 0
  BefPreSeffForEta = BefPreSeffForpT =  BefPreSeffForNumPixHits = BefPreSeffForValidFract = 0
  BefPreSeffForNumDeDx = BefPreSeffForChi2oDOF = BefPreSeffForEoP = BefPreSeffFordz = 0
  BefPreSeffFordxy = BefPreSeffForMiniIso = BefPreSeffForIh = BefPreSeffForProbXY = 0
  for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
      keyname = curr_dir.GetListOfKeys().At(i).GetName()
      curr_dir2 = SingleMuonSample.GetDirectory(dirname+"/"+keyname)
      if not (curr_dir2) :
        continue
      for j in range(0, curr_dir2.GetListOfKeys().GetEntries()):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          newname = dirname + "/" + keyname+ "/" + keyname2
          if not ("N1" in keyname2) : continue
          if ("Qual" in keyname2) : continue
          
          cstackedSummedBackgroundString = "cstackedSummedBackgroundString"+str(keyname2)
          canvas = ROOT.TCanvas(cstackedSummedBackgroundString, cstackedSummedBackgroundString, 800,800)
          
          histo = SingleMuonSample.Get(newname)
          histoAllBcg = AllBcgSample.Get(newname)
          histoSignal = Rhadron1800GeV.Get(newname)
          
          if not (histo.ClassName() == "TH1F") : continue
          
          numBins = histo.GetNbinsX()
          maxXaxis = histo.GetXaxis().GetXmax()
          minXaxis = histo.GetXaxis().GetXmin()
          stepSize = (maxXaxis-minXaxis)/numBins
          
          SoBstring = "Punzi"+str(j)

          PunziHisto = ROOT.TH1F(SoBstring,SoBstring,numBins,minXaxis,maxXaxis)
          PunziHisto.SetStats(0)
          PunziHisto.SetMarkerStyle(20)
          
          for x in numpy.arange(0,maxXaxis,stepSize):
            NumSignal = histoSignal.Integral(histoSignal.GetXaxis().FindBin(0),histoSignal.GetXaxis().FindBin(x))
            DenomSignal = histoSignal.Integral(histoSignal.GetXaxis().FindBin(0),histoSignal.GetXaxis().FindBin(histoSignal.GetXaxis().GetXmax()+1))
            NumBcg = histoAllBcg.Integral(histoAllBcg.GetXaxis().FindBin(0),histoAllBcg.GetXaxis().FindBin(x))
            DenomBcg = histoAllBcg.Integral()

            
            if (DenomSignal>0) : EffSignal = NumSignal / DenomSignal
            else : EffSignal = 0
            
#            print("EffSignal" + ": " + str(EffSignal))
#            print("sqrt(B)" + ": " + str(numpy.sqrt(NumBcg)))
            Bcg = NumBcg
            if ((keyname2 == "N1_TNOM") or (keyname2 == "N1_Pt") or (keyname2 == "N1_Pt_lowPt") or (keyname2 == "N1_TNOPH") or (keyname2 == "N1_TNOHFraction")) :
              EffSignal = 1-EffSignal
              Bcg = DenomBcg - NumBcg
            PunziForX = (EffSignal) / (sigma + numpy.sqrt(Bcg))

            
            PunziHisto.SetBinContent(histo.GetXaxis().FindBin(x),PunziForX)
          PunziHisto.Draw() #HISTTEXT00
          PunziHisto.SetTitle("")
          PunziHisto.GetYaxis().SetTitle("Punzi-significance")
          axisTitle = keyname2[3:]
          PunziHisto.GetXaxis().SetTitle(axisTitle)
            
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
          
          tex4 = ROOT.TLatex()

          tex5 = ROOT.TLatex(0.07,0.04,"Code version: "+codeVersion);
          tex5.SetNDC();
          tex5.SetTextFont(52);
          tex5.SetTextSize(0.0185);
          tex5.SetLineWidth(2);
                
          tex4 = ROOT.TLatex(0.5,0.95,"After (N-1)+1 selection")
          tex4.SetNDC();
          tex4.SetTextFont(52);
          tex4.SetTextSize(0.045);
          tex4.SetLineWidth(2);

          tex2.Draw("SAME")
          tex3.Draw("SAME")
          tex4.Draw("SAME")
          tex5.Draw("SAME")
    
          canvas.SaveAs("ComparePunzi_CodeV"+codeVersion+"/"+keyname2+".png")

os.system("cp forWebpage/* ComparePunzi_CodeV"+codeVersion+"/.")
os.system("cp forWebpage/.htaccess ComparePunzi_CodeV"+codeVersion+"/.")
print("scp -r ComparePunzi_CodeV"+ codeVersion + " tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/.")
