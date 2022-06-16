# BackgroundStudies

This repositery is dedicated to background estimate

## Run background estimate: read histograms & produce mass plots 

Configure configFile_readHist.txt respecting dedicated columns
```bash
root_file   rebin_boolean   rebin_eta   rebin_Ih    rebin_p     rebin_mass  
```
where rebin is a boolean in order to know if we want to rebin the different distributions. Then we give the values of rebinning for each distribution. 

then run,
```bash
root -l -q -b step2_backgroundPrediction.C
```
