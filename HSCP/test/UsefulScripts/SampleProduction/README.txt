#Execute the followin commands:
#
#copy gen fragments in your project
cmsenv
mkdir $CMSSW_BASE/src/Configuration
cp -rd /afs/cern.ch/user/c/cnuttens/public/TriggerSamples/Configuration/* $CMSSW_BASE/src/Configuration/.

#make GEN-SIM samples  (list of samples is defined at the head of Launch.py
python Launch.py 1

#make the RECO samples (to be executed when the GEN-SIM step is finished)
#python Launch.py 2
