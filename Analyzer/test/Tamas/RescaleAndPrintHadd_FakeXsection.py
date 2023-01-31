import ROOT, sys, os, time, re, numpy
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
"crab_Analysis_2018_QCD_Pt-50To80_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-80To120_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-120To170_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-170To300_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-300To470_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-470To600_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-600To800_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-800To1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_QCD_Pt-1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2018_WJetsToLNu_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_WJetsToLNu_0J_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_WJetsToLNu_1J_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_WJetsToLNu_2J_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTToHadronic_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTToSemiLeptonic_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_TTTo2L2Nu_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluino_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root",
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
"crab_Analysis_2018_HSCPgmsbStau_M-1029_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-1218_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-1409_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-1599_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-247_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-308_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-432_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-557_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-651_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-745_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPgmsbStau_M-871_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstopOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPstop_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge1e_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2018_HSCPtauPrimeCharge2e_M-500_CodeV"+codeVersion+"_v1.root",
]

#intLumi = 4598.37 #10.0 #
intLumi = 100000.0 #100/fb
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
'crab_Analysis_2018_HSCPgluino_M-500_CodeV'+codeVersion+'_v1.root' :   10.,
'crab_Analysis_2018_HSCPgluino_M-800_CodeV'+codeVersion+'_v1.root' :   10.,
'crab_Analysis_2018_HSCPgluino_M-1000_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-1200_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-1400_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-1600_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-1800_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-2000_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-2200_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-2400_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluino_M-2600_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-500_CodeV'+codeVersion+'_v1.root' :   10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-800_CodeV'+codeVersion+'_v1.root' :   10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1000_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1200_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1400_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1600_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-1800_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2000_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2200_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2400_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPgluinoOnlyNeutral_M-2600_CodeV'+codeVersion+'_v1.root' :  10.,
'crab_Analysis_2018_HSCPpairStau_M-200_CodeV'+codeVersion+'_v1.root' : 0.0303143567,
'crab_Analysis_2018_HSCPpairStau_M-247_CodeV'+codeVersion+'_v1.root' : 0.0129202991,
'crab_Analysis_2018_HSCPpairStau_M-308_CodeV'+codeVersion+'_v1.root' : 0.0048020183,
'crab_Analysis_2018_HSCPpairStau_M-432_CodeV'+codeVersion+'_v1.root' : 0.0012159719,
'crab_Analysis_2018_HSCPpairStau_M-557_CodeV'+codeVersion+'_v1.root' : 0.000389904,
'crab_Analysis_2018_HSCPpairStau_M-651_CodeV'+codeVersion+'_v1.root' : 0.0001234716,
'crab_Analysis_2018_HSCPpairStau_M-745_CodeV'+codeVersion+'_v1.root' : 5.8632201e-5,
'crab_Analysis_2018_HSCPpairStau_M-871_CodeV'+codeVersion+'_v1.root' : 2.9178605e-5,
'crab_Analysis_2018_HSCPpairStau_M-1029_CodeV'+codeVersion+'_v1.root' : 1.50415955e-5,
'crab_Analysis_2018_HSCPpairStau_M-1218_CodeV'+codeVersion+'_v1.root' : 0.00057/1000., #fake
'crab_Analysis_2018_HSCPpairStau_M-1409_CodeV'+codeVersion+'_v1.root' : 0.00057/1000., #fake
'crab_Analysis_2018_HSCPpairStau_M-1599_CodeV'+codeVersion+'_v1.root' : 0.00057/1000.  ,
'crab_Analysis_2018_HSCPgmsbStau_M-200_CodeV'+codeVersion+'_v1.root' :  2.8E-01,
'crab_Analysis_2018_HSCPgmsbStau_M-247_CodeV'+codeVersion+'_v1.root' :  8.8E-02,
'crab_Analysis_2018_HSCPgmsbStau_M-308_CodeV'+codeVersion+'_v1.root' :  2.5E-02,
'crab_Analysis_2018_HSCPgmsbStau_M-432_CodeV'+codeVersion+'_v1.root' :  3.9E-03,
'crab_Analysis_2018_HSCPgmsbStau_M-557_CodeV'+codeVersion+'_v1.root' :  1.9E-03,
'crab_Analysis_2018_HSCPgmsbStau_M-651_CodeV'+codeVersion+'_v1.root' :  4.1E-04,
'crab_Analysis_2018_HSCPgmsbStau_M-745_CodeV'+codeVersion+'_v1.root' :  1.9E-04,
'crab_Analysis_2018_HSCPgmsbStau_M-871_CodeV'+codeVersion+'_v1.root' :  6.9E-05,
'crab_Analysis_2018_HSCPgmsbStau_M-1029_CodeV'+codeVersion+'_v1.root' :  2.2E-05,
'crab_Analysis_2018_HSCPgmsbStau_M-1218_CodeV'+codeVersion+'_v1.root' :  6.4E-06,
'crab_Analysis_2018_HSCPgmsbStau_M-1409_CodeV'+codeVersion+'_v1.root' :  2.0E-06,
'crab_Analysis_2018_HSCPgmsbStau_M-1599_CodeV'+codeVersion+'_v1.root' :  5.3E-07,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-500_CodeV'+codeVersion+'_v1.root' :  0.000257,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-800_CodeV'+codeVersion+'_v1.root' :  0.0326,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-1000_CodeV'+codeVersion+'_v1.root' :  0.00683,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-1200_CodeV'+codeVersion+'_v1.root' :  0.0017,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-1400_CodeV'+codeVersion+'_v1.root' :  0.000473,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-1600_CodeV'+codeVersion+'_v1.root' :  0.000142,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-1800_CodeV'+codeVersion+'_v1.root' :  4.51e-05,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-2000_CodeV'+codeVersion+'_v1.root' :  1.48e-05,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-2200_CodeV'+codeVersion+'_v1.root' :  5e-06,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-2400_CodeV'+codeVersion+'_v1.root' :  1.71e-06,
'crab_Analysis_2018_HSCPstopOnlyNeutral_M-2600_CodeV'+codeVersion+'_v1.root' :  5.9e-07,
'crab_Analysis_2018_HSCPstop_M-500_CodeV'+codeVersion+'_v1.root' :  0.000257,
'crab_Analysis_2018_HSCPstop_M-800_CodeV'+codeVersion+'_v1.root' :  0.0326,
'crab_Analysis_2018_HSCPstop_M-1000_CodeV'+codeVersion+'_v1.root' :  0.00683,
'crab_Analysis_2018_HSCPstop_M-1200_CodeV'+codeVersion+'_v1.root' :  0.0017,
'crab_Analysis_2018_HSCPstop_M-1400_CodeV'+codeVersion+'_v1.root' :  0.000473,
'crab_Analysis_2018_HSCPstop_M-1600_CodeV'+codeVersion+'_v1.root' :  0.000142,
'crab_Analysis_2018_HSCPstop_M-1800_CodeV'+codeVersion+'_v1.root' :  4.51e-05,
'crab_Analysis_2018_HSCPstop_M-2000_CodeV'+codeVersion+'_v1.root' :  1.48e-05,
'crab_Analysis_2018_HSCPstop_M-2200_CodeV'+codeVersion+'_v1.root' :  5e-06,
'crab_Analysis_2018_HSCPstop_M-2400_CodeV'+codeVersion+'_v1.root' :  1.71e-06,
'crab_Analysis_2018_HSCPstop_M-2600_CodeV'+codeVersion+'_v1.root' :  5.9e-07,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-200_CodeV'+codeVersion+'_v1.root' :  1.1E-01,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-400_CodeV'+codeVersion+'_v1.root' :  7.3E-03,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-500_CodeV'+codeVersion+'_v1.root' :  1.2E-03,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-800_CodeV'+codeVersion+'_v1.root' :  2.6E-04,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1000_CodeV'+codeVersion+'_v1.root' :  7.6E-05,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1400_CodeV'+codeVersion+'_v1.root' :  8.5E-06,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-1800_CodeV'+codeVersion+'_v1.root' :  1.2E-06,
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2200_CodeV'+codeVersion+'_v1.root' :  8.2E-07, #fake
'crab_Analysis_2018_HSCPtauPrimeCharge1e_M-2600_CodeV'+codeVersion+'_v1.root' :  1.2E-07,#fake
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-200_CodeV'+codeVersion+'_v1.root' :  3.0E-01,
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-400_CodeV'+codeVersion+'_v1.root' :  2.3E-02,
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-500_CodeV'+codeVersion+'_v1.root' :  3.5E-03,
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1000_CodeV'+codeVersion+'_v1.root' :  2.4E-04,
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1400_CodeV'+codeVersion+'_v1.root' :  2.7E-05,
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-1800_CodeV'+codeVersion+'_v1.root' :  3.9E-06,
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2200_CodeV'+codeVersion+'_v1.root' :  8.2E-07, #fake
'crab_Analysis_2018_HSCPtauPrimeCharge2e_M-2600_CodeV'+codeVersion+'_v1.root' :  1.2E-07,#fake
}

## TTBar is 832 pb (NNLO)
# gluino x-sections from
# https://github.com/fuenfundachtzig/xsec/blob/master/json/pp13_gluino_NNLO%2BNNLL.json#L1820

fileInArray = []
for sample in BackgroundSamples:
  if not os.path.exists(sample): continue
  fileInArray.append(ROOT.TFile.Open(sample,"UPDATE"))

for fileIn in fileInArray:
  if not (fileIn.Get("HSCParticleAnalyzer/BaseName/NumEvents")):
    print("NumEvents not found, exit for "+str(fileIn))
    continue
  
  nEvetsPreTrig = fileIn.Get("HSCParticleAnalyzer/BaseName/NumEvents").GetBinContent(1)
  nEvetsPostTrig = fileIn.Get("HSCParticleAnalyzer/BaseName/NumEvents").GetBinContent(2)
  if (nEvetsPreTrig == 0):
    print("nEvetsPreTrig is zero, exiting")
    continue
  nameFromTFile = str(fileIn)[str(fileIn).find("Name")+6:str(fileIn).find("Title")-1]
  if not (crossSectionArray.get(nameFromTFile)) :
    print("No crossSectionArray for "+str(nameFromTFile))
    continue
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

print("hadd crab_Analysis_SingleMuon_RunPhase1_CodeV"+codeVersion+"_v1.root *SingleMuon*_CodeV"+codeVersion+"_v1.root")

#if not os.path.exists("CodeV"+codeVersion): os.mkdir("CodeV"+codeVersion)
#
#os.system("mv *"+codeVersion+"*root "+ "CodeV"+codeVersion+"/.")
#os.system("cp CodeV"+codeVersion+"/crab_Analysis_2018_SingleMuon_Run2018C_CodeV*"+codeVersion+"_v1.root .")
#os.system("cp CodeV"+codeVersion+"/crab_Analysis_2018_HSCPgluino_M-*_CodeV*"+codeVersion+"_v1.root .")
#os.system("mv CodeV"+codeVersion+"/crab_Analysis_2018_All*"+codeVersion+"_v1.root .")
