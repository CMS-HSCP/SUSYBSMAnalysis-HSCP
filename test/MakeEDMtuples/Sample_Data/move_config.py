import shutil
import os
import glob
import time
import os.path

source_sh = glob.glob('/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/Sample_Data/FARM/inputs/*.sh')
source_py = glob.glob('/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/Sample_Data/FARM/inputs/*.py')
source_err = glob.glob('/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/Sample_Data/FARM/logs/*.err')
destination = '/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/Sample_Data/FARM/outputs/'

for el in source_err:
    fileName = el.split('/')[-1]
    jobNumb = fileName.split('SCP')[0]

    for pyFile in source_py:
        if jobNumb in pyFile:
            print pyFile
            shutil.move(pyFile,destination)
    for shFile in source_sh:
        if jobNumb in shFile:
            print shFile
            shutil.move(shFile,destination)
                
