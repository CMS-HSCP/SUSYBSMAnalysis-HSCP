import sys, os
from optparse import OptionParser
from threading import Thread

parser = OptionParser(usage="Usage: python3 %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = []

codeVersion = sys.argv[1]
#just the number, like 18p2

didVoms = input("Push enter if you alread did voms-proxy-init -rfc -voms cms -valid 192:00 otherwise say no and do it\n")
if(didVoms):
 sys.exit()
 
for fname in os.listdir("crab_projects") :
  if (codeVersion in fname) :
    datasetList.append("crab_projects/"+fname)

def task(i):
  print("Download for sample "+i)
  outTask = "crab out -d "+i+" --checksum=no --jobids 12"
  os.system(outTask)
  haddTask = "hadd "+i[(i.find('crab_projects'))+14:]+".root "+i[(i.find('crab_projects')):]+"/results/*root"
#  os.system(haddTask)
  backgroundPred = "BackgroundPrediction -f "+i[(i.find('crab_projects'))+14:]+".root"
#  os.system(backgroundPred)
  print("Done for sample "+i)

for dataset in datasetList:
  t = Thread(target=task, args=(dataset,))
  t.start()
