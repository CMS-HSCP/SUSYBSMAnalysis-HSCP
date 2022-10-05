import sys, os
#from common_functions import *
from optparse import OptionParser
parser = OptionParser(usage="Usage: python3 %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = []

codeVersion = sys.argv[1]
#just the number, like 18p2

#didVoms = input("Push enter if you alread did voms-proxy-init -rfc -voms cms -valid 192:00 otherwise say no and do it\n")
#if(didVoms):
# sys.exit()
 
for fname in os.listdir("crab_projects") :
  if (fname.find(codeVersion)>0) :
    datasetList.append("crab_projects/"+fname)

for i in datasetList:
  print("Download for sample "+i)
#  createTask = "crab remake --task="+i
#  os.system(createTask)
#  outTask = "crab out -d "+i[(i.find('crab')):]+" --checksum=no --jobids 12"
  outTask = "crab out -d "+i+" --checksum=no --jobids 2"
  os.system(outTask)
#  print(outTask)
#  haddTask = "hadd "+i[(i.find('crab_projects'))+14:]+".root "+i[(i.find('crab_projects')):]+"/results/*root"
#  os.system(haddTask)
#  print(haddTask)
#  backgroundPred = "BackgroundPrediction -f "+i[(i.find('crab_projects'))+14:]+".root"
#  os.system(backgroundPred)
#  print(backgroundPred)
