import sys, os
from optparse import OptionParser
from threading import Thread

parser = OptionParser(usage="Usage: python3 %prog codeVersion")
(opt,args) = parser.parse_args()

didConda = input("Push enter if you alread did conda activate withRoot or say no and do it\n")
if(didConda):
 sys.exit()

datasetList = []

codeVersion = sys.argv[1]
#just the number, like 18p2
 
for fname in os.listdir(".") :
  if (codeVersion in fname) and (".root" in fname) :
    datasetList.append(fname)

def task(i):
  runPngFromRoot3 = "python3 pngsFromRootFilesJustHistos.py "+i+" 3"
  os.system(runPngFromRoot3)
  os.system("cp forWebpage/* "+i[0:-5]+"_Bin3/.")
  runPngFromRoot25 = "python3 pngsFromRootFilesJustHistos.py "+i+" 25"
#  os.system(runPngFromRoot25)
#  os.system("cp forWebpage/* "+i[0:-5]+"_Bin25/.")
  
  scpBin3 = "scp -r "+i[0:-5]+"_Bin3 tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/."
#  os.system(scpBin3)
  scpBin25 = "scp -r "+i[0:-5]+"_Bin25 tvami@lxplus.cern.ch:/eos/home-t/tvami/www/projects/HSCP/2022CodeV"+codeVersion+"/."
#  os.system(scpBin25)
  print("Done for sample "+i)

for dataset in datasetList:
  t = Thread(target=task, args=(dataset,))
  t.start()
