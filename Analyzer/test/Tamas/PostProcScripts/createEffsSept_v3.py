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
'BefPreS_TIsol' : 10,
'PostPreS_TIsol' : 10,
'BefPreS_MiniTkIso' : 10,
'N1_MiniTkIso' : 10,
'BefPreS_MiniTkIso_PUA' : 10,
'N1_MiniTkIso_PUA' : 10,
'BefPreS_MiniTkIso_PUB' : 10,
'N1_MiniTkIso_PUB' : 10,
'BefPreS_MiniTkIso_PUC' : 10,
'N1_MiniTkIso_PUC' : 10,
}


sigma = 3

if not os.path.exists(os.path.dirname("Effs_CodeV"+codeVersion+"/a.png")):
  print("Create dir")
  os.makedirs(os.path.dirname("Effs_CodeV"+codeVersion+"/"))

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
#          if not ("N1_Eta" in keyname2) : continue
          newname = dirname + "/" + keyname+ "/" + keyname2
          if ("N1_ProbQVsIas" in keyname2 or "N1_pfType" in keyname2) : continue
          histo = SingleMuonSample.Get(newname)
          histoAllBcg = AllBcgSample.Get(newname)
          histoSignal = Rhadron1800GeV.Get(newname)
#          histoDenomAllBcg = AllBcgSample.Get(dirname + "/" + keyname + "/N1_Eta")
#          histoDenomAllBcg = AllBcgSample.Get(dirname + "/" + keyname + "/BefPreS_Eta")

          if not (cutValues.get(keyname2)) : continue
          Num = histo.Integral(histo.GetXaxis().FindBin(0),histo.GetXaxis().FindBin(cutValues[keyname2]))
          NumAllBcg = histoAllBcg.Integral(histoAllBcg.GetXaxis().FindBin(0),histoAllBcg.GetXaxis().FindBin(cutValues[keyname2]))
          NumSignal = histoSignal.Integral(histoSignal.GetXaxis().FindBin(0),histoSignal.GetXaxis().FindBin(cutValues[keyname2]))
          Denom = histo.Integral(histo.GetXaxis().FindBin(0),histo.GetXaxis().FindBin(histo.GetXaxis().GetXmax()+1))
          DenomAllBcg = histoAllBcg.Integral(histoAllBcg.GetXaxis().FindBin(0),histoAllBcg.GetXaxis().FindBin(histoAllBcg.GetXaxis().GetXmax()+1))
#          DenomAllBcg = histoDenomAllBcg.Integral(histoDenomAllBcg.GetXaxis().FindBin(0),histoDenomAllBcg.GetXaxis().FindBin(histoDenomAllBcg.GetXaxis().GetXmax()+1))
#          DenomAllBcg = histoDenomAllBcg.Integral(histoDenomAllBcg.GetXaxis().FindBin(0),histoDenomAllBcg.GetXaxis().FindBin(histoDenomAllBcg.GetXaxis().GetXmax()+1))
#
#          print(keyname2 + ": " + str(DenomAllBcg))
          DenomSignal = histoSignal.Integral(histoSignal.GetXaxis().FindBin(0),histoSignal.GetXaxis().FindBin(histoSignal.GetXaxis().GetXmax()+1))
          if (Denom>0) :  Eff = Num / Denom
          else : Eff = 0
          if (DenomAllBcg>0) : EffAllBcg = NumAllBcg / DenomAllBcg
          else : EffAllBcg = 0
          if (DenomSignal>0) : EffSignal = NumSignal / DenomSignal
          else : EffSignal = 0
         
          if ("BefPreS_Eta" in keyname2) :
            BefPreSeffForEta = Eff
            BefPreSeffForEtaBcg = EffAllBcg
            BefPreSeffForEtaSignal = EffSignal
            BefPreSeffForEtaPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_P") :
            BefPreSeffForP = Eff
            BefPreSeffForPBcg = EffAllBcg
            BefPreSeffForPSignal = EffSignal
            BefPreSeffForPPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_Pt" in keyname2) :
            BefPreSeffForpT = 1-Eff
            BefPreSeffForpTBcg = 1-EffAllBcg
            BefPreSeffForpTSignal = 1-EffSignal
            BefPreSeffForpTPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_PtErrOverPt2" in keyname2) :
            BefPreSeffForpTErrOverPt2 = Eff
            BefPreSeffForpTErrOverPt2Bcg = EffAllBcg
            BefPreSeffForpTErrOverPt2Signal = EffSignal
            BefPreSeffForpTErrOverPt2Punzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_TNOPH" in keyname2) :
            BefPreSeffForNumPixHits = 1-Eff
            BefPreSeffForNumPixHitsBcg = 1-EffAllBcg
            BefPreSeffForNumPixHitsSignal = 1-EffSignal
            BefPreSeffForNumPixHitsPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_TNOHFraction" in keyname2) :
            BefPreSeffForValidFract = 1-Eff
            BefPreSeffForValidFractBcg = 1-EffAllBcg
            BefPreSeffForValidFractSignal = 1-EffSignal
            BefPreSeffForValidFractPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_TNOM" in keyname2) :
            BefPreSeffForNumDeDx = 1-Eff
            BefPreSeffForNumDeDxBcg = 1-EffAllBcg
            BefPreSeffForNumDeDxSignal = 1-EffSignal
            BefPreSeffForNumDeDxPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_Chi2oNdof" in keyname2) :
            BefPreSeffForChi2oDOF = Eff
            BefPreSeffForChi2oDOFBcg = EffAllBcg
            BefPreSeffForChi2oDOFSignal = EffSignal
            BefPreSeffForChi2oDOFPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_EoP" in keyname2) :
            BefPreSeffForEoP = Eff
            BefPreSeffForEoPBcg = EffAllBcg
            BefPreSeffForEoPSignal = EffSignal
            BefPreSeffForEoPPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_Dz" in keyname2) :
            BefPreSeffFordz = Eff
            BefPreSeffFordzBcg = EffAllBcg
            BefPreSeffFordzSignal = EffSignal
            BefPreSeffFordzPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_Dxy" in keyname2) :
            BefPreSeffFordxy = Eff
            BefPreSeffFordxyBcg = EffAllBcg
            BefPreSeffFordxySignal = EffSignal
            BefPreSeffFordxyPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_MiniRelIsoAll_lowMiniRelIso" in keyname2) :
            BefPreSeffForMiniIso = Eff
            BefPreSeffForMiniIsoBcg = EffAllBcg
            BefPreSeffForMiniIsoSignal = EffSignal
            BefPreSeffForMiniIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_Ih" in keyname2) :
            BefPreSeffForIh = 1-Eff
            BefPreSeffForIhBcg = 1-EffAllBcg
            BefPreSeffForIhSignal = 1-EffSignal
            BefPreSeffForIhPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("BefPreS_ProbXY" in keyname2) :
            BefPreSeffForProbXY = 1-Eff
            BefPreSeffForProbXYBcg = 1-EffAllBcg
            BefPreSeffForProbXYSignal = 1-EffSignal
            BefPreSeffForProbXYPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_TIsol") :
            BefPreSeffForTIsol = Eff
            BefPreSeffForTIsolBcg = EffAllBcg
            BefPreSeffForTIsolSignal = EffSignal
            BefPreSeffForTIsolPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniRelTkIso_lowMiniRelIso") :
            BefPreSeffForMiniRelTkIso = Eff
            BefPreSeffForMiniRelTkIsoBcg = EffAllBcg
            BefPreSeffForMiniRelTkIsoSignal = EffSignal
            BefPreSeffForMiniRelTkIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniRelTkIso_lowMiniRelIso_PUA") :
            BefPreSeffForMiniRelTkIsoPUA = Eff
            BefPreSeffForMiniRelTkIsoPUABcg = EffAllBcg
            BefPreSeffForMiniRelTkIsoPUASignal = EffSignal
            BefPreSeffForMiniRelTkIsoPUAPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniRelTkIso_lowMiniRelIso_PUB") :
            BefPreSeffForMiniRelTkIsoPUB = Eff
            BefPreSeffForMiniRelTkIsoPUBBcg = EffAllBcg
            BefPreSeffForMiniRelTkIsoPUBSignal = EffSignal
            BefPreSeffForMiniRelTkIsoPUBPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniRelTkIso_lowMiniRelIso_PUC") :
            BefPreSeffForMiniRelTkIsoPUC = Eff
            BefPreSeffForMiniRelTkIsoPUCBcg = EffAllBcg
            BefPreSeffForMiniRelTkIsoPUCSignal = EffSignal
            BefPreSeffForMiniRelTkIsoPUCPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniTkIso") :
            BefPreSeffForMiniTkIso = Eff
            BefPreSeffForMiniTkIsoBcg = EffAllBcg
            BefPreSeffForMiniTkIsoSignal = EffSignal
            BefPreSeffForMiniTkIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniTkIso_PUA") :
            BefPreSeffForMiniTkIsoPUA = Eff
            BefPreSeffForMiniTkIsoPUABcg = EffAllBcg
            BefPreSeffForMiniTkIsoPUASignal = EffSignal
            BefPreSeffForMiniTkIsoPUAPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniTkIso_PUB") :
            BefPreSeffForMiniTkIsoPUB = Eff
            BefPreSeffForMiniTkIsoPUBBcg = EffAllBcg
            BefPreSeffForMiniTkIsoPUBSignal = EffSignal
            BefPreSeffForMiniTkIsoPUBPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "BefPreS_MiniTkIso_PUC") :
            BefPreSeffForMiniTkIsoPUC = Eff
            BefPreSeffForMiniTkIsoPUCBcg = EffAllBcg
            BefPreSeffForMiniTkIsoPUCSignal = EffSignal
            BefPreSeffForMiniTkIsoPUCPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
            
            
          if ("N1_Eta" in keyname2) :
            N1effForEta = Eff
            N1effForEtaBcg = EffAllBcg
            N1effForEtaSignal = EffSignal
            N1effForEtaPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_Pt" in keyname2) :
            N1effForpT = 1-Eff
            N1effForpTBcg = 1-EffAllBcg
            N1effForpTSignal = 1-EffSignal
            N1effForpTPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("PostPreS_PtErrOverPt2" in keyname2) :
            N1effForpTErrOverPt2 = Eff
            N1effForpTErrOverPt2Bcg = EffAllBcg
            N1effForpTErrOverPt2Signal = EffSignal
            N1effForpTErrOverPt2Punzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_TNOPH" in keyname2) :
            N1effForNumPixHits = 1-Eff
            N1effForNumPixHitsBcg = 1-EffAllBcg
            N1effForNumPixHitsSignal = 1-EffSignal
            N1effForNumPixHitsPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_TNOHFraction" in keyname2) :
            N1effForValidFract = 1-Eff
            N1effForValidFractBcg = 1-EffAllBcg
            N1effForValidFractSignal = 1-EffSignal
            N1effForValidFractPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_TNOM" in keyname2) :
            N1effForNumDeDx = 1-Eff
            N1effForNumDeDxBcg = 1-EffAllBcg
            N1effForNumDeDxSignal = 1-EffSignal
            N1effForNumDeDxPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_Chi2oNdof" in keyname2) :
            N1effForChi2oDOF = Eff
            N1effForChi2oDOFBcg = EffAllBcg
            N1effForChi2oDOFSignal = EffSignal
            N1effForChi2oDOFPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_EoP" in keyname2) :
            N1effForEoP = Eff
            N1effForEoPBcg = EffAllBcg
            N1effForEoPSignal = EffSignal
            N1effForEoPPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_Dz" in keyname2) :
            N1effFordz = Eff
            N1effFordzBcg = EffAllBcg
            N1effFordzSignal = EffSignal
            N1effFordzPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_Dxy" in keyname2) :
            N1effFordxy = Eff
            N1effFordxyBcg = EffAllBcg
            N1effFordxySignal = EffSignal
            N1effFordxyPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_MiniRelIsoAll_lowMiniRelIso" in keyname2) :
            N1effForMiniIso = Eff
            N1effForMiniIsoBcg = EffAllBcg
            N1effForMiniIsoSignal = EffSignal
            print("MiniRelIsoAll EffSignal: "+str(EffSignal))
            print("DenomAllBcg: "+str(numpy.sqrt(NumAllBcg)))
            N1effForMiniIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_Ih" in keyname2) :
            N1effForIh = 1-Eff
            N1effForIhBcg = 1-EffAllBcg
            N1effForIhSignal = 1-EffSignal
            N1effForIhPunzi = (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if ("N1_ProbXY" in keyname2) :
            N1effForProbXY = 1-Eff
            N1effForProbXYBcg = 1-EffAllBcg
            N1effForProbXYSignal = 1-EffSignal
            N1effForProbXYPunzi= (1-EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "PostPreS_P" ) :
            N1effForP = Eff
            N1effForPBcg = EffAllBcg
            N1effForPSignal = EffSignal
            N1effForPPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "PostPreS_TIsol") :
            N1effForTIsol = Eff
            N1effForTIsolBcg = EffAllBcg
            N1effForTIsolSignal = EffSignal
            N1effForTIsolPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniRelTkIso_lowMiniRelIso") :
            N1effForMiniRelTkIso = Eff
            N1effForMiniRelTkIsoBcg = EffAllBcg
            N1effForMiniRelTkIsoSignal = EffSignal
            print("MiniRelTkIso EffSignal: "+str(EffSignal))
            print("DenomAllBcg: "+str(numpy.sqrt(NumAllBcg)))
            N1effForMiniRelTkIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniRelTkIso_lowMiniRelIso_PUA") :
            N1effForMiniRelTkIsoPUA = Eff
            N1effForMiniRelTkIsoPUABcg = EffAllBcg
            N1effForMiniRelTkIsoPUASignal = EffSignal
#            print("PUA EffSignal: "+str(EffSignal))
#            print("DenomAllBcg: "+str(numpy.sqrt(NumAllBcg)))
            N1effForMiniRelTkIsoPUAPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniRelTkIso_lowMiniRelIso_PUB") :
            N1effForMiniRelTkIsoPUB = Eff
            N1effForMiniRelTkIsoPUBBcg = EffAllBcg
            N1effForMiniRelTkIsoPUBSignal = EffSignal
#            print("PUB EffSignal: "+str(EffSignal))
#            print("DenomAllBcg: "+str(numpy.sqrt(NumAllBcg)))
            N1effForMiniRelTkIsoPUBPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniRelTkIso_lowMiniRelIso_PUC") :
            N1effForMiniRelTkIsoPUC = Eff
            N1effForMiniRelTkIsoPUCBcg = EffAllBcg
            N1effForMiniRelTkIsoPUCSignal = EffSignal
#            print("PUC EffSignal: "+str(EffSignal))
#            print("DenomAllBcg: "+str(numpy.sqrt(NumAllBcg)))
            N1effForMiniRelTkIsoPUCPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniTkIso") :
            N1effForMiniTkIso = Eff
            N1effForMiniTkIsoBcg = EffAllBcg
            N1effForMiniTkIsoSignal = EffSignal
            N1effForMiniTkIsoPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniTkIso_PUA") :
            N1effForMiniTkIsoPUA = Eff
            N1effForMiniTkIsoPUABcg = EffAllBcg
            N1effForMiniTkIsoPUASignal = EffSignal
            N1effForMiniTkIsoPUAPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniTkIso_PUB") :
            N1effForMiniTkIsoPUB = Eff
            N1effForMiniTkIsoPUBBcg = EffAllBcg
            N1effForMiniTkIsoPUBSignal = EffSignal
            N1effForMiniTkIsoPUBPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          if (keyname2 == "N1_MiniTkIso_PUC") :
            N1effForMiniTkIsoPUC = Eff
            N1effForMiniTkIsoPUCBcg = EffAllBcg
            N1effForMiniTkIsoPUCSignal = EffSignal
            N1effForMiniTkIsoPUCPunzi = (EffSignal) / (sigma + numpy.sqrt(NumAllBcg))
          
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
      
####################################################################################
      # EffsInCutflowN1.png
      cstackedSummedBackgroundStringN1 = "cstackedSummedBackgroundStringN1"+str(keyname2)
      canvas = ROOT.TCanvas(cstackedSummedBackgroundStringN1, cstackedSummedBackgroundStringN1, 800,800)
      
      EffForData = ROOT.TH1F("EffForDataN1",";;Integral to (from) max (min) cut value /  Integral to inf",10,0.,10.)
      EffForData.Draw("HISTTEXT00")
      EffForData.SetBinContent(1,N1effForMiniRelTkIso)
      EffForData.SetBinContent(2,N1effForMiniRelTkIsoPUA)
      EffForData.SetBinContent(3,N1effForMiniRelTkIsoPUB)
      EffForData.SetBinContent(4,N1effForMiniRelTkIsoPUC)
      EffForData.SetBinContent(5,N1effForMiniTkIso)
      EffForData.SetBinContent(6,N1effForMiniTkIsoPUA)
      EffForData.SetBinContent(7,N1effForMiniTkIsoPUB)
      EffForData.SetBinContent(8,N1effForMiniTkIsoPUC)
      EffForData.SetBinContent(9,N1effForMiniIso)
      EffForData.SetBinContent(10,N1effForTIsol)

      
      EffForData.SetLineColor(1)
      EffForData.SetMarkerColor(1)
#      EffForData.SetMarkerStyle(20)
      EffForData.SetStats(0)
      EffForData.GetXaxis().SetTitle("")
      EffForData.GetYaxis().SetRangeUser(0.,1.3)
      EffForData.GetXaxis().SetBinLabel(1,"MiniRelTkIso")
      EffForData.GetXaxis().SetBinLabel(2,"MiniRelTkIsoPUA")
      EffForData.GetXaxis().SetBinLabel(3,"MiniRelTkIsoPUB")
      EffForData.GetXaxis().SetBinLabel(4,"MiniRelTkIsoPUC")
      EffForData.GetXaxis().SetBinLabel(5,"MiniTkIso")
      EffForData.GetXaxis().SetBinLabel(6,"MiniTkIsoPUA")
      EffForData.GetXaxis().SetBinLabel(7,"MiniTkIsoPUB")
      EffForData.GetXaxis().SetBinLabel(8,"MiniTkIsoPUC")
      EffForData.GetXaxis().SetBinLabel(9,"MiniRelIso < 0.02")
      EffForData.GetXaxis().SetBinLabel(10,"TIsol < 15 GeV")
      EffForData.GetXaxis().SetTitle("")
      
      legend =  ROOT.TLegend(.55,.80,.80,.89,"","brNDC")
      legend.SetTextFont(42)
      legend.SetTextSize(0.02)
      legend.SetBorderSize(1);
      legend.SetBorderSize(0);
      legend.SetLineColor(1);
      legend.SetLineStyle(1);
      legend.SetLineWidth(1);
      legend.SetFillColor(0);
      legend.SetFillStyle(1001);
      legend.AddEntry(EffForData,"DATA (2018C)","LP")
      

      EffForBcg = ROOT.TH1F("EffForBcgN1",";;Integral to (from) max (min) cut value /  Integral to inf",10,0.,10.)
      EffForBcg.Draw("SAMEHISTOTEXT00")
      EffForBcg.SetLineColor(3)
      EffForBcg.SetMarkerColor(3)
      EffForBcg.SetBinContent(1,N1effForMiniRelTkIsoBcg)
      EffForBcg.SetBinContent(2,N1effForMiniRelTkIsoPUABcg)
      EffForBcg.SetBinContent(3,N1effForMiniRelTkIsoPUBBcg)
      EffForBcg.SetBinContent(4,N1effForMiniRelTkIsoPUCBcg)
      EffForBcg.SetBinContent(5,N1effForMiniTkIsoBcg)
      EffForBcg.SetBinContent(6,N1effForMiniTkIsoPUABcg)
      EffForBcg.SetBinContent(7,N1effForMiniTkIsoPUBBcg)
      EffForBcg.SetBinContent(8,N1effForMiniTkIsoPUCBcg)
      EffForBcg.SetBinContent(9,N1effForMiniIsoBcg)
      EffForBcg.SetBinContent(10,N1effForTIsolBcg)
      EffForBcg.SetBinContent(14,N1effForpTErrOverPt2Bcg)
      legend.AddEntry(EffForBcg,"ALLBcg MC","LP")
      
      EffForSignal = ROOT.TH1F("EffForSignalN1",";;Integral to (from) max (min) cut value /  Integral to inf",10,0.,10.)
      EffForSignal.Draw("SAMEHISTOTEXT00")
      EffForSignal.SetLineColor(6)
      EffForSignal.SetMarkerColor(6)
      EffForSignal.SetBinContent(1,N1effForMiniRelTkIsoSignal)
      EffForSignal.SetBinContent(2,N1effForMiniRelTkIsoPUASignal)
      EffForSignal.SetBinContent(3,N1effForMiniRelTkIsoPUBSignal)
      EffForSignal.SetBinContent(4,N1effForMiniRelTkIsoPUCSignal)
      EffForSignal.SetBinContent(5,N1effForMiniTkIsoSignal)
      EffForSignal.SetBinContent(6,N1effForMiniTkIsoPUASignal)
      EffForSignal.SetBinContent(7,N1effForMiniTkIsoPUBSignal)
      EffForSignal.SetBinContent(8,N1effForMiniTkIsoPUCSignal)
      EffForSignal.SetBinContent(9,N1effForMiniIsoSignal)
      EffForSignal.SetBinContent(10,N1effForTIsolSignal)
      EffForSignal.SetBinContent(14,N1effForpTErrOverPt2Signal)
      legend.AddEntry(EffForSignal,"HSCP Gluino 1800 GeV","LP")
      
      tex4 = ROOT.TLatex(0.55,0.95,"After (N-1)+1 selection")
      tex4.SetNDC();
      tex4.SetTextFont(52);
      tex4.SetTextSize(0.045);
      tex4.SetLineWidth(2);
      
      legend.Draw("SAME")
      tex2.Draw("SAME")
      tex3.Draw("SAME")
      tex4.Draw("SAME")
      tex5.Draw("SAME")
      
      canvas.SaveAs("Effs_CodeV"+codeVersion+"/EffsInCutflowN1.png")

####################################################################################################################
      # Punzi for N-1 (EffsInCutflowN1Punzi.png)
      
      ROOT.gStyle.SetPaintTextFormat(".4f");
      cstackedSummedBackgroundStringN1Punzi = "cstackedSummedBackgroundStringN1Punzi"+str(keyname2)
      canvas = ROOT.TCanvas(cstackedSummedBackgroundStringN1Punzi, cstackedSummedBackgroundStringN1Punzi, 800,800)
      
      EffForData = ROOT.TH1F("EffForDataPunzi",";;Punzi-significance",10,0.,10.)
      EffForData.Draw("HISTTEXT00")
      EffForData.SetBinContent(1,N1effForMiniRelTkIsoPunzi)
      EffForData.SetBinContent(2,N1effForMiniRelTkIsoPUAPunzi)
      EffForData.SetBinContent(3,N1effForMiniRelTkIsoPUBPunzi)
      EffForData.SetBinContent(4,N1effForMiniRelTkIsoPUCPunzi)
      EffForData.SetBinContent(5,N1effForMiniTkIsoPunzi)
      EffForData.SetBinContent(6,N1effForMiniTkIsoPUAPunzi)
      EffForData.SetBinContent(7,N1effForMiniTkIsoPUBPunzi)
      EffForData.SetBinContent(8,N1effForMiniTkIsoPUCPunzi)
      EffForData.SetBinContent(9,N1effForMiniIsoPunzi)
      EffForData.SetBinContent(10,N1effForTIsolPunzi)

      
      EffForData.SetLineColor(1)
      EffForData.SetMarkerColor(1)
      EffForData.SetStats(0)
      EffForData.GetXaxis().SetTitle("")
      EffForData.GetXaxis().SetBinLabel(1,"MiniRelTkIso")
      EffForData.GetXaxis().SetBinLabel(2,"MiniRelTkIsoPUA")
      EffForData.GetXaxis().SetBinLabel(3,"MiniRelTkIsoPUB")
      EffForData.GetXaxis().SetBinLabel(4,"MiniRelTkIsoPUC")
      EffForData.GetXaxis().SetBinLabel(5,"MiniTkIso")
      EffForData.GetXaxis().SetBinLabel(6,"MiniTkIsoPUA")
      EffForData.GetXaxis().SetBinLabel(7,"MiniTkIsoPUB")
      EffForData.GetXaxis().SetBinLabel(8,"MiniTkIsoPUC")
      EffForData.GetXaxis().SetBinLabel(9,"MiniRelIso < 0.02")
      EffForData.GetXaxis().SetBinLabel(10,"TIsol < 15 GeV")
      EffForData.GetXaxis().SetTitle("")
      
      tex4 = ROOT.TLatex(0.55,0.95,"After (N-1)+1 selection")
      tex4.SetNDC();
      tex4.SetTextFont(52);
      tex4.SetTextSize(0.045);
      tex4.SetLineWidth(2);
      
      tex2.Draw("SAME")
      tex3.Draw("SAME")
      tex4.Draw("SAME")
      tex5.Draw("SAME")
      
      canvas.SaveAs("Effs_CodeV"+codeVersion+"/EffsInCutflowN1Punzi.png")

os.system("cp forWebpage/* Effs_CodeV"+codeVersion+"/.")
os.system("cp forWebpage/.htaccess Effs_CodeV"+codeVersion+"/.")
print("scp -r Effs_CodeV"+ codeVersion + " tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/.")
