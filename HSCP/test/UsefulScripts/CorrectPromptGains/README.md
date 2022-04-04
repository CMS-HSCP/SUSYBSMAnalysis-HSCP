# README

Simple `python` script which produces correct strip tracker gains file for prompt datasets.

## Introduction
As gains are computed over a run range with enough statistics, the correct gains (also payloads) are unknown during the time of the datataking. This is why some old payloads are applied on the prompt data. Payloads are applied on the `APV` level.

After enough statistics have been collected, correct payloads are computed. However, we cannot use them, because they are meant to be applied on data which has no payloads applied "online". So to get correct payloads, one must compute

```
correct gains = applied gains * (correct gains / applied gains)
```
meaning the correct gains for _prompt_ datasets are (correct gains / applied gains). This is not needed for _re-reco_ as those datasets have correct gains already applied.

## Step 1

Gains applied online on prompt data can be found using `conddb`, which gives you run which had those gains applied at first and the hash checksum of the payload file. You do not know which file it is. It could be any of the files listed in the directories `/afs/cern.ch/cms/tracker/sistrvalidation/WWW/CalibrationValidation/ParticleGain`. Names of the directories listed there tell us which run range was used to compute the payloads. Script checks the `Gains_Sqlite.db` files and searches for `PAYLOAD_HASH` entry in the `IOV` table to find a matching payload.

Now we know which payloads are applied on prompt data for which run range. Correct gains are any of the gains listed in the previously mentioned directory. Script then merely combines them, and creates a file `gains.txt` in the first step which can be edited. Script also supports whitelisting or blacklisting specific payloads which may be detrimental to use for the analysis.

File `gains.txt` has the following format:

```
<Gains_firstRun_to_LastRun> <Correct payloads> <Applied payloads>
```


To run the first step, run
```
python Launch.py 1
```

## Step 2

After the `gains.txt` file has been checked and everything is in order, `CombineGains.sh` is run, which is just a shell wrapper that runs `CombineGains.C`. This creates a new gains file that is valid for a run range specified in `gains.txt`. Lastly the script uses `hadd -f` to merge those files and combines it with `2015` gains file.

To run the second step, producing the final gains file, run
```
python Launch.py 2
```

## Control check

The gains file should then be placed in the `../../../data` directory for use in the entire analysis code and `../StabilityCheck` should be run on it to see its impact on the overall stability.

## TO-DO

`CombineGains.C` usually segfaults: it tries to produce the final file twice, and since that file exists already it crashes. All the files required to run are still produced and are good, so fixing is not necessary, but it would be nice to have it run without segfaults.
