import ROOT, sys, os, time, re, numpy
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python %prog codeVersion")
from tqdm import tqdm

codeVersion = sys.argv[1]
#just the number, like 18p2

ROOT.gROOT.SetBatch(True)
ROOT.gStyle.SetPadRightMargin(.15)
ROOT.gStyle.SetPadTopMargin(0.1);
ROOT.gStyle.SetPadBottomMargin(0.14);
ROOT.gStyle.SetPadLeftMargin(0.15);

BackgroundSamples = [
#"crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
##"crab_Analysis_2018_WJetsToLNu_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_WJetsToLNu_0J_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_WJetsToLNu_1J_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_WJetsToLNu_2J_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_TTToHadronic_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_TTToSemiLeptonic_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_TTTo2L2Nu_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-1000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-1400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-1600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-2000_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-2200_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-2400_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-2600_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-500_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_HSCPgluino_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-557_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-247_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-308_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-432_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-557_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-651_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-745_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-871_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1029_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1218_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1409_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPpairStau_M-1599_CodeV"+codeVersion+"_v1.root",
]

#intLumi = 4598.37 #10.0 #
intLumi = 2000.0
#intLumi = 30000.0

crossSectionArray = {
'crab_Analysis_2018_QCD_Pt-20_MuEnrichedPt15_CodeV'+codeVersion+'_v1.root' : 239000, #+-755.0,
'crab_Analysis_2018_QCD_Pt-15To20_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 2797000, #+-8800.0,
'crab_Analysis_2018_QCD_Pt-20To30_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 2518000, #+-7940.0,
'crab_Analysis_2018_QCD_Pt-30To50_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 1361000, #+-4263.0,
'crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 377800.0, #+-1184.0,
'crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 88620.0, #+-275.1
'crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 21070.0, #+-65.28
'crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 7019.0, #+-21.61
'crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  622.4, #+-1.891
'crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  58.86, #+-0.1776
'crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  18.22, #+-0.05471
'crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  3.25, #+-0.0148
'crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  1.61, #+-
# 61526.7, #+-	61526.7, WJetsToLNu, this is also what FCP used, (NNLO)
'crab_Analysis_2018_WJetsToLNu_0J_CodeV'+codeVersion+'_v1.root' :  53330.0,#	+-90.89
'crab_Analysis_2018_WJetsToLNu_1J_CodeV'+codeVersion+'_v1.root' : 8875.0, #	+-55.31
'crab_Analysis_2018_WJetsToLNu_2J_CodeV'+codeVersion+'_v1.root' :  3338.0, #	+-34.64
'crab_Analysis_2018_TTToHadronic_CodeV'+codeVersion+'_v1.root' : 377.96, # or 377.96,
'crab_Analysis_2018_TTToSemiLeptonic_CodeV'+codeVersion+'_v1.root' : 365.35, # or 65.34
'crab_Analysis_2018_TTTo2L2Nu_CodeV'+codeVersion+'_v1.root' :  88.29, # or 687.1 or 91.47600
'crab_Analysis_2018_HSCPgluino_M-500_CodeV'+codeVersion+'_v1.root' :   33800,
'crab_Analysis_2018_HSCPgluino_M-800_CodeV'+codeVersion+'_v1.root' :   1810,
'crab_Analysis_2018_HSCPgluino_M-1000_CodeV'+codeVersion+'_v1.root' :  385 ,
'crab_Analysis_2018_HSCPgluino_M-1400_CodeV'+codeVersion+'_v1.root' :  28.4 ,
'crab_Analysis_2018_HSCPgluino_M-1600_CodeV'+codeVersion+'_v1.root' :  8.87 ,
'crab_Analysis_2018_HSCPgluino_M-1800_CodeV'+codeVersion+'_v1.root' :  2.93 ,
'crab_Analysis_2018_HSCPgluino_M-2000_CodeV'+codeVersion+'_v1.root' :  1.01 ,
'crab_Analysis_2018_HSCPgluino_M-2200_CodeV'+codeVersion+'_v1.root' :  0.356 ,
'crab_Analysis_2018_HSCPgluino_M-2400_CodeV'+codeVersion+'_v1.root' :  0.128 ,
'crab_Analysis_2018_HSCPgluino_M-2600_CodeV'+codeVersion+'_v1.root' :  0.0462 ,
'crab_Analysis_2018_HSCPpairStau_M-200_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-247_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-308_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-432_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-557_CodeV'+codeVersion+'_v1.root' : 0.001  ,
'crab_Analysis_2018_HSCPpairStau_M-651_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-745_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-871_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-1029_CodeV'+codeVersion+'_v1.root' : 0.001,
'crab_Analysis_2018_HSCPpairStau_M-1218_CodeV'+codeVersion+'_v1.root' : 0.001, #fake
'crab_Analysis_2018_HSCPpairStau_M-1409_CodeV'+codeVersion+'_v1.root' : 0.001, #fake
'crab_Analysis_2018_HSCPpairStau_M-1599_CodeV'+codeVersion+'_v1.root' : 0.001  ,
}

## TTBar is 832 pb (NNLO)
# gluino x-sections from
# https://github.com/fuenfundachtzig/xsec/blob/master/json/pp13_gluino_NNLO%2BNNLL.json#L1820

fileInArray = []
for sample in BackgroundSamples:
  fileInArray.append(ROOT.TFile.Open(sample,"UPDATE"))

for fileIn in fileInArray:
  if not (fileIn.Get("HSCParticleAnalyzer/BaseName/NumEvents")):
    print("NumEvents not found, exit")
    continue
  
  nEvetsPreTrig = fileIn.Get("HSCParticleAnalyzer/BaseName/NumEvents").GetBinContent(1)
  nEvetsPostTrig = fileIn.Get("HSCParticleAnalyzer/BaseName/NumEvents").GetBinContent(2)
  if (nEvetsPreTrig == 0):
    print("nEvetsPreTrig is zero, exiting")
    continue
  nameFromTFile = str(fileIn)[str(fileIn).find("Name")+6:str(fileIn).find("Title")-1]
  weight = intLumi*crossSectionArray.get(nameFromTFile)/nEvetsPreTrig
  
  for i in range(0, fileIn.GetListOfKeys().GetEntries()):
    dirname = fileIn.GetListOfKeys().At(i).GetName()
    curr_dir = fileIn.GetDirectory(dirname)
#    print("dirname: "+dirname)
    if not (curr_dir) :
      continue
    for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
        keyname = curr_dir.GetListOfKeys().At(i).GetName()
        curr_dir2 = fileIn.GetDirectory(dirname+"/"+keyname)
        if not (curr_dir2) :
          continue
        for j in tqdm(range(0, curr_dir2.GetListOfKeys().GetEntries())):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
#          print("keyname2: "+keyname2)
          newname = dirname + "/" + keyname+ "/" + keyname2
          obj = fileIn.Get(newname)
          if not (obj) : continue
          if (obj.GetEntries() == 0 ) :
#            print("obj.GetEntries() == 0")
            continue
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if obj.InheritsFrom("TObject"):
            obj.Scale(weight)
  fileIn.Write("",ROOT.TObject.kOverwrite)
  fileIn.Close()

print("hadd crab_Analysis_2018_AllBackground_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_WJetsToLNu_*J_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToHadronic_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToSemiLeptonic_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTTo2L2Nu_CodeV"+codeVersion+"_v1.root")

#os.system
print("hadd crab_Analysis_2018_AllQCD_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root")

print("hadd crab_Analysis_2018_AllTTbar_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToHadronic_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToSemiLeptonic_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTTo2L2Nu_CodeV"+codeVersion+"_v1.root")

print("hadd crab_Analysis_2018_AllWJets_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_WJetsToLNu*J_CodeV"+codeVersion+"_v1.root ")

#if not os.path.exists("CodeV"+codeVersion): os.mkdir("CodeV"+codeVersion)
#
#os.system("mv *"+codeVersion+"*root "+ "CodeV"+codeVersion+"/.")
#os.system("cp CodeV"+codeVersion+"/crab_Analysis_2018_SingleMuon_Run2018C_CodeV*"+codeVersion+"_v1.root .")
#os.system("cp CodeV"+codeVersion+"/crab_Analysis_2018_HSCPgluino_M-*_CodeV*"+codeVersion+"_v1.root .")
#os.system("mv CodeV"+codeVersion+"/crab_Analysis_2018_All*"+codeVersion+"_v1.root .")
