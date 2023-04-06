import sys, os
from optparse import OptionParser
from threading import Thread

parser = OptionParser(usage="Usage: python3 %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = []

codeVersion = sys.argv[1]
#just the number, like 18p2

for fname in os.listdir(".") :
  if (codeVersion in fname and "root" in fname) :
    if ("gluino_" in fname):
      if not os.path.exists("HSCPgluinoSamples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPgluinoSamples_CodeV40p9.txt HSCPgluinoSamples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPgluinoSamples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)
    elif ("gluinoOnlyNeutral_" in fname):
      if not os.path.exists("HSCPgluinoOnlyNeutralSamples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPgluinoOnlyNeutralSamples_CodeV40p9.txt HSCPgluinoOnlyNeutralSamples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPgluinoOnlyNeutralSamples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)
    elif ("stop_" in fname):
      if not os.path.exists("HSCPstopSamples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPstopSamples_CodeV40p9.txt HSCPstopSamples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPstopSamples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)
    elif ("stopOnlyNeutral_" in fname):
      if not os.path.exists("HSCPstopOnlyNeutralSamples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPstopOnlyNeutralSamples_CodeV40p9.txt HSCPstopOnlyNeutralSamples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPstopOnlyNeutralSamples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)
    elif ("pairStau_" in fname):
      if not os.path.exists("HSCPpairStauSamples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPpairStauSamples_CodeV40p9.txt HSCPpairStauSamples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPpairStauSamples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)
    elif ("gmsbStau_" in fname):
      if not os.path.exists("HSCPgmsbStauSamples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPgmsbStauSamples_CodeV40p9.txt HSCPgmsbStauSamples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPgmsbStauSamples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)
    elif ("tauPrimeCharge1e_" in fname):
      if not os.path.exists("HSCPtauPrime1Samples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPtauPrime1Samples_CodeV40p9.txt HSCPtauPrime1Samples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPtauPrime1Samples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)
    elif ("tauPrimeCharge2e_" in fname):
      if not os.path.exists("HSCPtauPrime2Samples_CodeV"+codeVersion+".txt"):
        os.system("cp HSCPtauPrime2Samples_CodeV40p9.txt HSCPtauPrime2Samples_CodeV"+codeVersion+".txt")
        replaceMINTA = "sed -i '' 's/40p9/"+codeVersion+"/g' HSCPtauPrime2Samples_CodeV"+codeVersion+".txt"
        os.system(replaceMINTA)

for fname in os.listdir(".") :
  if (codeVersion in fname and "txt" in fname and "Samples_" in fname) :
    datasetList.append(fname)

for dataset in datasetList:
  outTask = "python3 compareWithArguementList.py "+dataset
  print(outTask)
  
print("python3 pngsFromRootFilesJustHistos.py crab_Analysis_2018_AllBackground_CodeV"+codeVersion+"_v1.root 3")
print("python3 pngsFromRootFilesJustHistos.py crab_Analysis_2018_HSCPgluino_M-1800_CodeV"+codeVersion+".root 3")
print("python3 pngsFromRootFilesJustHistos.py crab_Analysis_SingleMuon_RunPhase1_CodeV"+codeVersion+".root 3")
