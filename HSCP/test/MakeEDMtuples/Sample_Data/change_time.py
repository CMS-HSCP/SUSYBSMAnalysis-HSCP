import shutil
import os
import glob
import time
import os.path
import fileinput
source = glob.glob('/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/Sample_Data/FARM/inputs/29*HSCPEdmProd.sh')
print source
for el in source:
    for line in fileinput.input(el, inplace=True): 
      print line.rstrip().replace('#SBATCH --time=24:00:00', '#SBATCH --time=12:00:00')
      #print line.rstrip().replace('voms-proxy-init ;','voms-proxy-init --noregen;')

