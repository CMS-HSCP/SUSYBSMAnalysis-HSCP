import shutil
import os
import glob
import time
import os.path

while True:
    source = glob.glob('/home/ucl/cp3/jpriscia/CMSSW_9_4_3/src/SUSYBSMAnalysis/HSCP/test/MakeEDMtuples/Sample_Data/out/*/*.txt')
    print source
    destination = "/storage/data/cms/store/user/jpriscia/2017ReReco/"
    for files in source:
        if files.endswith(".txt"):
            print "file", files.replace( ".txt", ".root" )
            if os.path.isfile(files.replace( ".txt", ".root" )):
                shutil.move(files.replace( ".txt", ".root" ),destination)
                os.rename(files, files.replace( ".txt",".done" ))
            else: 
                print  files, " does not have correspondent root file"
    print "I moved ", len(source) , "files" 
    time.sleep(10)
    #threading.Timer(600,moveOutputs).start()  #called every  10 minutes


