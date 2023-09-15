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
"crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_WJetsToLNu_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTToHadronic_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTToSemiLeptonic_woProbQ_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTTo2L2Nu_woProbQ_CodeV"+codeVersion+"_v1.root",
]


if not os.path.exists("CodeV"+codeVersion+"_old"): os.mkdir("CodeV"+codeVersion+"_old")
os.system("cp *"+codeVersion+"*root "+ "CodeV"+codeVersion+"_old/.")

#this needs to be measured in pb
intLumi = 4598.37 #137.0 #10.0 #

#crossSections are measured in 1/pb
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

for sample,xsec in zip(BackgroundSamples,crossSectionArray):
  fileIn = ROOT.TFile.Open(sample,"UPDATE")
  print("fileIn ",fileIn)
  fileIn.cd("analyzer/BaseName")
  h = ROOT.gDirectory.Get("NumEvents").Clone()
  nEventsPreTrig = h.GetBinContent(1)
   
  if (nEventsPreTrig == 0):
    print("nEvetsPreTrig is zero, exiting")
    continue
  weight = intLumi*xsec/nEventsPreTrig
  keys = ROOT.gDirectory.GetListOfKeys().Clone()
  for k in tqdm(keys):
#    print(k.GetName())
    if ROOT.gDirectory.Get(k.GetName()) != None:
      h = ROOT.gDirectory.Get(k.GetName()).Clone()
#      if h.GetEntries()==0 : 
##        print (k.GetName(), " has zero entries")
#        continue
      if k.GetName() != "HscpCandidates" and k.GetName() != "GenHscpCandidates":
        h.Scale(weight)
      h.Write("",ROOT.TObject.kOverwrite)
#    else: print(k.GetName(), " is either a null pointer or not a histogram (likely the latter)")

#os.system("hadd crab_Analysis_2018_AllBackground_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_WJetsToLNu_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToHadronic_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToSemiLeptonic_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTTo2L2Nu_woProbQ_CodeV"+codeVersion+"_v1.root")
#
#os.system("hadd crab_Analysis_2018_AllQCD_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_woProbQ_CodeV"+codeVersion+"_v1.root")

#os.system("hadd crab_Analysis_2018_AllTTbar_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToHadronic_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTToSemiLeptonic_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_TTTo2L2Nu_woProbQ_CodeV"+codeVersion+"_v1.root")
#
#os.system("mv crab_Analysis_2018_WJetsToLNu_woProbQ_CodeV"+codeVersion+"_v1.root crab_Analysis_2018_AllWJetsToLNu_woProbQ_CodeV"+codeVersion+"_v1.root")
#
#if not os.path.exists("CodeV"+codeVersion): os.mkdir("CodeV"+codeVersion)
#
#os.system("mv *"+codeVersion+"*root "+ "CodeV"+codeVersion+"/.")
#os.system("mv CodeV"+codeVersion+"/crab_Analysis_2018_All*_woProbQ_CodeV"+codeVersion+"_v1.root .")
#os.system("mv CodeV"+codeVersion+"/crab_Analysis_SingleMuon_Run2018_woProbQ_CodeV"+codeVersion+"_v1.root .")
