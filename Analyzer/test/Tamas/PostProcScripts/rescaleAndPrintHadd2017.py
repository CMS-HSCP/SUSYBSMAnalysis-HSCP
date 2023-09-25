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
"crab_Analysis_2017_QCD_Pt-50To80_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-80To120_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-120To170_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-170To300_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-300To470_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-470To600_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-600To800_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-800To1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_QCD_Pt-1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root",
#"crab_Analysis_2017_WJetsToLNu_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_WJetsToLNu_0J_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_WJetsToLNu_1J_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_WJetsToLNu_2J_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_TTToHadronic_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_TTToSemiLeptonic_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_TTTo2L2Nu_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-50To120_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-120To200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-200To400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-400To800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-800To1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-1400To2300_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-2300To3500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-3500To4500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-4500To6000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_ZToMuMu_M-6000ToInf_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluino_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-247_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-308_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-432_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-557_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-651_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-745_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-871_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-1029_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-1218_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-1409_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPpairStau_M-1599_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-1029_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-1218_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-1409_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-1599_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-247_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-308_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-432_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-557_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-651_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-745_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPgmsbStau_M-871_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-1200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstopOnlyNeutral_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-1200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-1600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-2000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-2400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPstop_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-500_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge1e_M-800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-1000_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-1400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-1800_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-2200_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-2600_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-400_CodeV"+codeVersion+"_v1.root",
"crab_Analysis_2017_HSCPtauPrimeCharge2e_M-500_CodeV"+codeVersion+"_v1.root",
]

#intLumi = 4598.37 #10.0 #
intLumi = 100000.0 #100/fb
#intLumi = 30000.0

crossSectionArray = {
'crab_Analysis_2017_QCD_Pt-20_MuEnrichedPt15_CodeV'+codeVersion+'_v1.root' : 239000, #+-755.0,
'crab_Analysis_2017_QCD_Pt-15To20_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 2797000, #+-8800.0,
'crab_Analysis_2017_QCD_Pt-20To30_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 2518000, #+-7940.0,
'crab_Analysis_2017_QCD_Pt-30To50_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 1361000, #+-4263.0,
'crab_Analysis_2017_QCD_Pt-50To80_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 377800.0, #+-1184.0,
'crab_Analysis_2017_QCD_Pt-80To120_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 88620.0, #+-275.1
'crab_Analysis_2017_QCD_Pt-120To170_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 21070.0, #+-65.28
'crab_Analysis_2017_QCD_Pt-170To300_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' : 7019.0, #+-21.61
'crab_Analysis_2017_QCD_Pt-300To470_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  622.4, #+-1.891
'crab_Analysis_2017_QCD_Pt-470To600_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  58.86, #+-0.1776
'crab_Analysis_2017_QCD_Pt-600To800_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  18.22, #+-0.05471
'crab_Analysis_2017_QCD_Pt-800To1000_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  3.25, #+-0.0148
'crab_Analysis_2017_QCD_Pt-1000_MuEnrichedPt5_CodeV'+codeVersion+'_v1.root' :  1.61, #+-
# 61526.7, #+-	61526.7, WJetsToLNu, this is also what FCP used, (NNLO)
'crab_Analysis_2017_WJetsToLNu_0J_CodeV'+codeVersion+'_v1.root' :  53330.0,#	+-90.89
'crab_Analysis_2017_WJetsToLNu_1J_CodeV'+codeVersion+'_v1.root' : 8875.0, #	+-55.31
'crab_Analysis_2017_WJetsToLNu_2J_CodeV'+codeVersion+'_v1.root' :  3338.0, #	+-34.64
'crab_Analysis_2017_TTToHadronic_CodeV'+codeVersion+'_v1.root' : 377.96, # or 377.96,
'crab_Analysis_2017_TTToSemiLeptonic_CodeV'+codeVersion+'_v1.root' : 365.35, # or 65.34
'crab_Analysis_2017_TTTo2L2Nu_CodeV'+codeVersion+'_v1.root' :  88.29, # or 687.1 or 91.47600
'crab_Analysis_2017_ZToMuMu_M-50To120_CodeV'+codeVersion+'_v1.root' : 2.116e+03, #+-9.883e-01
'crab_Analysis_2017_ZToMuMu_M-120To200_CodeV'+codeVersion+'_v1.root' : 2.058e+01, # +/- 1.362e-02
'crab_Analysis_2017_ZToMuMu_M-200To400_CodeV'+codeVersion+'_v1.root' : 2.890e+00, # +/- 1.977e-03
'crab_Analysis_2017_ZToMuMu_M-400To800_CodeV'+codeVersion+'_v1.root' : 2.515e-01, # +/- 1.770e-04
'crab_Analysis_2017_ZToMuMu_M-800To1400_CodeV'+codeVersion+'_v1.root' : 1.709e-02, # +/- 1.232e-05
'crab_Analysis_2017_ZToMuMu_M-1400To2300_CodeV'+codeVersion+'_v1.root' : 1.370e-03, # +/- 9.867e-07
'crab_Analysis_2017_ZToMuMu_M-2300To3500_CodeV'+codeVersion+'_v1.root' : 8.282e-05, #  +/- 5.740e-08
'crab_Analysis_2017_ZToMuMu_M-3500To4500_CodeV'+codeVersion+'_v1.root' : 4.650e-06,
'crab_Analysis_2017_ZToMuMu_M-4500To6000_CodeV'+codeVersion+'_v1.root' : 3.650e-07,# +/- 1.874e-10
'crab_Analysis_2017_ZToMuMu_M-6000ToInf_CodeV'+codeVersion+'_v1.root' : 2.526e-08, # +/- 2.331e-11
'crab_Analysis_2017_HSCPgluino_M-500_CodeV'+codeVersion+'_v1.root' :   33.800,
'crab_Analysis_2017_HSCPgluino_M-800_CodeV'+codeVersion+'_v1.root' :   1.810,
'crab_Analysis_2017_HSCPgluino_M-1000_CodeV'+codeVersion+'_v1.root' :  0.385,
'crab_Analysis_2017_HSCPgluino_M-1200_CodeV'+codeVersion+'_v1.root' :  0.0985,
'crab_Analysis_2017_HSCPgluino_M-1400_CodeV'+codeVersion+'_v1.root' :  0.0284,
'crab_Analysis_2017_HSCPgluino_M-1600_CodeV'+codeVersion+'_v1.root' :  0.00887,
'crab_Analysis_2017_HSCPgluino_M-1800_CodeV'+codeVersion+'_v1.root' :  0.00293,
'crab_Analysis_2017_HSCPgluino_M-2000_CodeV'+codeVersion+'_v1.root' :  0.00101,
'crab_Analysis_2017_HSCPgluino_M-2200_CodeV'+codeVersion+'_v1.root' :  0.000356,
'crab_Analysis_2017_HSCPgluino_M-2400_CodeV'+codeVersion+'_v1.root' :  0.000128,
'crab_Analysis_2017_HSCPgluino_M-2600_CodeV'+codeVersion+'_v1.root' :  4.62e-5,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-500_CodeV'+codeVersion+'_v1.root' :   33.800,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-800_CodeV'+codeVersion+'_v1.root' :   1.810,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1000_CodeV'+codeVersion+'_v1.root' :  0.385,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1200_CodeV'+codeVersion+'_v1.root' :  0.0985,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1400_CodeV'+codeVersion+'_v1.root' :  0.0284,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1600_CodeV'+codeVersion+'_v1.root' :  0.00887,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-1800_CodeV'+codeVersion+'_v1.root' :  0.00293,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2000_CodeV'+codeVersion+'_v1.root' :  0.00101,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2200_CodeV'+codeVersion+'_v1.root' :  0.000356,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2400_CodeV'+codeVersion+'_v1.root' :  0.000128,
'crab_Analysis_2017_HSCPgluinoOnlyNeutral_M-2600_CodeV'+codeVersion+'_v1.root' :  4.62e-5,
'crab_Analysis_2017_HSCPpairStau_M-200_CodeV'+codeVersion+'_v1.root' : 0.0303143567,
'crab_Analysis_2017_HSCPpairStau_M-247_CodeV'+codeVersion+'_v1.root' : 0.01354210931669562,
'crab_Analysis_2017_HSCPpairStau_M-308_CodeV'+codeVersion+'_v1.root' : 0.005617969677491303,
'crab_Analysis_2017_HSCPpairStau_M-432_CodeV'+codeVersion+'_v1.root' : 0.0013205993021946,
'crab_Analysis_2017_HSCPpairStau_M-557_CodeV'+codeVersion+'_v1.root' : 0.0003982919391389629,
'crab_Analysis_2017_HSCPpairStau_M-651_CodeV'+codeVersion+'_v1.root' : 0.00018170273864904894,
'crab_Analysis_2017_HSCPpairStau_M-745_CodeV'+codeVersion+'_v1.root' : 8.760649473515753e-5,
'crab_Analysis_2017_HSCPpairStau_M-871_CodeV'+codeVersion+'_v1.root' : 3.548748280460506e-5,
'crab_Analysis_2017_HSCPpairStau_M-1029_CodeV'+codeVersion+'_v1.root' : 1.1684554271190975e-5,
'crab_Analysis_2017_HSCPpairStau_M-1218_CodeV'+codeVersion+'_v1.root' : 0.00057/1000., #fake
'crab_Analysis_2017_HSCPpairStau_M-1409_CodeV'+codeVersion+'_v1.root' : 0.00057/1000., #fake
'crab_Analysis_2017_HSCPpairStau_M-1599_CodeV'+codeVersion+'_v1.root' : 0.00057/1000.  ,
'crab_Analysis_2017_HSCPgmsbStau_M-200_CodeV'+codeVersion+'_v1.root' :  2.8E-01,
'crab_Analysis_2017_HSCPgmsbStau_M-247_CodeV'+codeVersion+'_v1.root' :  8.8E-02,
'crab_Analysis_2017_HSCPgmsbStau_M-308_CodeV'+codeVersion+'_v1.root' :  2.5E-02,
'crab_Analysis_2017_HSCPgmsbStau_M-432_CodeV'+codeVersion+'_v1.root' :  3.9E-03,
'crab_Analysis_2017_HSCPgmsbStau_M-557_CodeV'+codeVersion+'_v1.root' :  1.9E-03,
'crab_Analysis_2017_HSCPgmsbStau_M-651_CodeV'+codeVersion+'_v1.root' :  4.1E-04,
'crab_Analysis_2017_HSCPgmsbStau_M-745_CodeV'+codeVersion+'_v1.root' :  1.9E-04,
'crab_Analysis_2017_HSCPgmsbStau_M-871_CodeV'+codeVersion+'_v1.root' :  6.9E-05,
'crab_Analysis_2017_HSCPgmsbStau_M-1029_CodeV'+codeVersion+'_v1.root' :  2.2E-05,
'crab_Analysis_2017_HSCPgmsbStau_M-1218_CodeV'+codeVersion+'_v1.root' :  6.4E-06,
'crab_Analysis_2017_HSCPgmsbStau_M-1409_CodeV'+codeVersion+'_v1.root' :  2.0E-06,
'crab_Analysis_2017_HSCPgmsbStau_M-1599_CodeV'+codeVersion+'_v1.root' :  5.3E-07,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-500_CodeV'+codeVersion+'_v1.root' :  0.000257,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-800_CodeV'+codeVersion+'_v1.root' :  0.0326,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-1000_CodeV'+codeVersion+'_v1.root' :  0.00683,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-1200_CodeV'+codeVersion+'_v1.root' :  0.0017,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-1400_CodeV'+codeVersion+'_v1.root' :  0.000473,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-1600_CodeV'+codeVersion+'_v1.root' :  0.000142,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-1800_CodeV'+codeVersion+'_v1.root' :  4.51e-05,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-2000_CodeV'+codeVersion+'_v1.root' :  1.48e-05,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-2200_CodeV'+codeVersion+'_v1.root' :  5e-06,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-2400_CodeV'+codeVersion+'_v1.root' :  1.71e-06,
'crab_Analysis_2017_HSCPstopOnlyNeutral_M-2600_CodeV'+codeVersion+'_v1.root' :  5.9e-07,
'crab_Analysis_2017_HSCPstop_M-500_CodeV'+codeVersion+'_v1.root' :  0.000257,
'crab_Analysis_2017_HSCPstop_M-800_CodeV'+codeVersion+'_v1.root' :  0.0326,
'crab_Analysis_2017_HSCPstop_M-1000_CodeV'+codeVersion+'_v1.root' :  0.00683,
'crab_Analysis_2017_HSCPstop_M-1200_CodeV'+codeVersion+'_v1.root' :  0.0017,
'crab_Analysis_2017_HSCPstop_M-1400_CodeV'+codeVersion+'_v1.root' :  0.000473,
'crab_Analysis_2017_HSCPstop_M-1600_CodeV'+codeVersion+'_v1.root' :  0.000142,
'crab_Analysis_2017_HSCPstop_M-1800_CodeV'+codeVersion+'_v1.root' :  4.51e-05,
'crab_Analysis_2017_HSCPstop_M-2000_CodeV'+codeVersion+'_v1.root' :  1.48e-05,
'crab_Analysis_2017_HSCPstop_M-2200_CodeV'+codeVersion+'_v1.root' :  5e-06,
'crab_Analysis_2017_HSCPstop_M-2400_CodeV'+codeVersion+'_v1.root' :  1.71e-06,
'crab_Analysis_2017_HSCPstop_M-2600_CodeV'+codeVersion+'_v1.root' :  5.9e-07,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-200_CodeV'+codeVersion+'_v1.root' :  1.1E-01,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-400_CodeV'+codeVersion+'_v1.root' :  7.3E-02,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-500_CodeV'+codeVersion+'_v1.root' :  1.2E-03,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-800_CodeV'+codeVersion+'_v1.root' :  2.6E-04,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-1000_CodeV'+codeVersion+'_v1.root' :  7.6E-05,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-1400_CodeV'+codeVersion+'_v1.root' :  8.5E-06,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-1800_CodeV'+codeVersion+'_v1.root' :  1.2E-06,
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-2200_CodeV'+codeVersion+'_v1.root' :  8.2E-07, #fake
'crab_Analysis_2017_HSCPtauPrimeCharge1e_M-2600_CodeV'+codeVersion+'_v1.root' :  1.2E-07,#fake
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-200_CodeV'+codeVersion+'_v1.root' :  3.0E-01,
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-400_CodeV'+codeVersion+'_v1.root' :  2.3E-02,
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-500_CodeV'+codeVersion+'_v1.root' :  3.5E-03,
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-1000_CodeV'+codeVersion+'_v1.root' :  2.4E-04,
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-1400_CodeV'+codeVersion+'_v1.root' :  2.7E-05,
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-1800_CodeV'+codeVersion+'_v1.root' :  3.9E-06,
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-2200_CodeV'+codeVersion+'_v1.root' :  8.2E-07, #fake
'crab_Analysis_2017_HSCPtauPrimeCharge2e_M-2600_CodeV'+codeVersion+'_v1.root' :  1.2E-07,#fake
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

print("hadd crab_Analysis_2017_AllBackground_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-50To80_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-80To120_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-120To170_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-170To300_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-300To470_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-470To600_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-600To800_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-800To1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_WJetsToLNu_*J_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_TTToHadronic_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_TTToSemiLeptonic_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_TTTo2L2Nu_CodeV"+codeVersion+"_v1.root")

#os.system
print("hadd crab_Analysis_2017_AllQCD_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-50To80_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-80To120_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-120To170_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-170To300_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-300To470_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-470To600_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-600To800_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-800To1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_QCD_Pt-1000_MuEnrichedPt5_CodeV"+codeVersion+"_v1.root")

print("hadd crab_Analysis_2017_AllTTbar_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_TTToHadronic_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_TTToSemiLeptonic_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_TTTo2L2Nu_CodeV"+codeVersion+"_v1.root")

print("hadd crab_Analysis_2017_AllWJets_CodeV"+codeVersion+"_v1.root crab_Analysis_2017_WJetsToLNu*J_CodeV"+codeVersion+"_v1.root ")

print("hadd crab_Analysis_SingleMuon_RunPhase1_CodeV"+codeVersion+"_v1.root *SingleMuon*_CodeV"+codeVersion+"_v1.root")

#if not os.path.exists("CodeV"+codeVersion): os.mkdir("CodeV"+codeVersion)
#
#os.system("mv *"+codeVersion+"*root "+ "CodeV"+codeVersion+"/.")
#os.system("cp CodeV"+codeVersion+"/crab_Analysis_2017_SingleMuon_Run2017C_CodeV*"+codeVersion+"_v1.root .")
#os.system("cp CodeV"+codeVersion+"/crab_Analysis_2017_HSCPgluino_M-*_CodeV*"+codeVersion+"_v1.root .")
#os.system("mv CodeV"+codeVersion+"/crab_Analysis_2017_All*"+codeVersion+"_v1.root .")
