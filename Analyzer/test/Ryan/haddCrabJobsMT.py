import sys, os
from optparse import OptionParser
from threading import Thread

parser = OptionParser(usage="Usage: python3 %prog codeVersion")
(opt,args) = parser.parse_args()

datasetList = []

codeVersion = sys.argv[1]
#just the number, like 18p2
 
for fname in os.listdir("crab_projects") :
  if (codeVersion in fname) :
    datasetList.append("crab_projects/"+fname)

def task(i):
  haddTask = "hadd "+i[(i.find('crab_projects'))+14:]+".root "+i[(i.find('crab_projects')):]+"/results/*root"
  os.system(haddTask)

for dataset in datasetList:
  t = Thread(target=task, args=(dataset,))
  t.start()

