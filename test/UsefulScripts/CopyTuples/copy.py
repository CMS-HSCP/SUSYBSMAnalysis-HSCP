#! /usr/bin/env python

import subprocess, os, sys

# note: make sure you have valid proxy

def main():
    version = "15_03_25_HSCP_Run2EDMFiles"
    src = ("root://eoscms.cern.ch/", "/eos/cms/store/cmst3/user/querten/"+version+"/")
    trg = ("root://se.cis.gov.pl/",  "/cms/store/user/fruboes/HSCP/"+version+"/")
    FNULL = open(os.devnull, 'w')

    # TODO: lustre (???)

    print "Copying from", src[0]+src[1]
    print "Copying to", trg[0]+trg[1]
    val = subprocess.call(["xrdfs", trg[0], "mkdir", "-p", trg[1] ], stderr=FNULL)
    if val==54:
        print "Output dir exists, you can probably ignore (above/this) warnings"


    filelist = ([f.strip() 
                for f in subprocess.check_output(["xrdfs", src[0], "ls", src[1]]).split("\n") 
                if len(f) > 0])

    for i,f in enumerate(filelist):
        print f
        fname = os.path.basename(f)
        print "doing {} ({}/{})".format(fname,i+1,len(filelist))
        source = src[0]+f
        targetPartiall = os.path.join(trg[1],fname)
        targetFull = trg[0]+targetPartiall
        exists = subprocess.call(["xrdfs", trg[0], "stat", targetPartiall], stdout=FNULL, stderr=FNULL)

        if exists == 0:
            print "   file exists, skipping"
            continue

        val = subprocess.call(["xrdcp", source, targetFull])


if __name__ == "__main__":
    main()

