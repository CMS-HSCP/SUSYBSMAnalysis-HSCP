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
"crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_WJetsToLNu_0J_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTToHadronic_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTToSemiLeptonic_wProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTTo2L2Nu_wProbQ_CodeV"+codeVersion+"_v1.root",
]

intLumi = 10.0 # 137.0

crossSectionArray = [
# QCD_Pt-20_MuEnrichedPt15_TuneCP5_13TeV-pythia8, 239000.0 +-	755.8
# 2797000.0, #+-8800.0, 2018_QCD_Pt-15To20_MuEnrichedPt5
# 2518000.0, #+-7940.0, 2018_QCD_Pt-20To30_MuEnrichedPt5
# 1361000.0, #+-4263.0, 2018_QCD_Pt-30To50_MuEnrichedPt5
 3778000.0, #+-1184.0, 2018_QCD_Pt-50To80_MuEnrichedPt5
 886200.0, #+-275.1, 2018_QCD_Pt-80To120_MuEnrichedPt5
 210700.0, #+-65.28, 2018_QCD_Pt-120To170_MuEnrichedPt5
 70190.0, #+-21.61, 2018_QCD_Pt-170To300_MuEnrichedPt5
 6224.0, #+-1.891, 2018_QCD_Pt-300To470_MuEnrichedPt5
 588.6, #+-0.1776, 2018_QCD_Pt-470To600_MuEnrichedPt5
 182.2, #+-0.05471, 2018_QCD_Pt-600To800_MuEnrichedPt5
 32.5, #+-0.0148, 2018_QCD_Pt-600To800_MuEnrichedPt5
 16.1, #+-, 2018_QCD_Pt-1000_MuEnrichedPt5
 53330.0, #+-	61526.7, WJetsToLNu
 687.1, # or 377.96, # or , TTtoHadronic or 382.53600
 719.1, # or 65.34, TTtoSemiLeptonic or 365.90400
 88.29, # or 687.1, TTto2L2N or 91.47600
]

fileInArray = []
for sample in BackgroundSamples:
  fileInArray.append(ROOT.TFile.Open(sample,"UPDATE"))

for index,fileIn in enumerate(fileInArray):
  if not (fileIn.Get("analyzer/BaseName/NumEvents")):
    continue
  
  nEvetsPostTrig = fileIn.Get("analyzer/BaseName/NumEvents").GetBinContent(2)
  if (nEvetsPostTrig == 0):
    continue
  weight = intLumi*crossSectionArray[index]/nEvetsPostTrig
  
  for i in range(0, fileIn.GetListOfKeys().GetEntries()):
    dirname = fileIn.GetListOfKeys().At(i).GetName()
    curr_dir = fileIn.GetDirectory(dirname)
  # print("dirname: "+dirname)
    if not (curr_dir) :
      continue
    for i in range(0, curr_dir.GetListOfKeys().GetEntries()):
        keyname = curr_dir.GetListOfKeys().At(i).GetName()
        curr_dir2 = fileIn.GetDirectory(dirname+"/"+keyname)
        if not (curr_dir2) :
          continue
        for j in tqdm(range(0, curr_dir2.GetListOfKeys().GetEntries())):
          keyname2 = curr_dir2.GetListOfKeys().At(j).GetName()
          newname = dirname + "/" + keyname+ "/" + keyname2
          obj = fileIn.Get(newname)
          if (obj.GetEntries() == 0 ) :
#            print("obj.GetEntries() == 0")
            continue
          if (keyname2=="HscpCandidates" or keyname2=="GenHscpCandidates"):
            continue
          if obj.InheritsFrom("TObject"):
            obj.Scale(weight)
  fileIn.Write("",ROOT.TObject.kOverwrite)
  fileIn.Close()

os.system("hadd crab_Analysis_2018_AllBackground_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_WJetsToLNu_0J_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToHadronic_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToSemiLeptonic_wProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTTo2L2Nu_wProbQ_CodeV"+codeVersion+"_v1.root")

if not os.path.exists("CodeV"+codeVersion): os.mkdir("CodeV"+codeVersion)

os.system("mv *"+codeVersion+"*root "+ "CodeV"+codeVersion+"/.")
os.system("mv CodeV"+codeVersion+"/crab_Analysis_2018_AllBackground_wProbQ_CodeV"+codeVersion+"_v1.root .")
