#!/bin/bash

#=============================================================================
#
#   Ce script ne fonctionne que sur les ui a IPHC
#   Il sert a merger des root files en parallele
#   Il faut indiquer le path du type /my/path/to/0000/*.root
#   ou encore du type /my/path/to/0001/*.root sans ajouter les
#   subdirectories 0000/ 0001/ ... a la fin
#   On peut regler le nombre de jobs que l on veut en reglant 
#   le nombre de fichiers par jobs. On Peut aussi decider de merger
#   les ttrees ou non (opt. -T) 
#   L option -j X sert a paralleliser sur X processeurs
#
#=============================================================================


filesperjobs=100
odir='/opt/sbg/cms/ui3_data1/dapparu/HSCP/Production/'
idir='/dpm/in2p3.fr/home/cms/phedex/store/user/dapparu/HSCP/Analysis/v2-0/prodJuly2021_CMSSW_10_6_20/MET/Analysis_MET_UL2017B_firstprod/210703_091704'
#idir='/dpm/in2p3.fr/home/cms/phedex/store/user/dapparu/HSCP/Analysis/v2-1/prodJuly2021_CMSSW_10_6_20/MET/Analysis_MET_UL2017B_297047-297666/210706_125248'
pathin2p3='root://sbgse1.in2p3.fr:/'
#options='-j 8 -T -f'
options='-f'

listofsubdir=`rfdir $idir | awk '{print $9'}`

for subdir in ${listofsubdir}; do
    indice_job=1
    listoffiles=`rfdir ${idir}/${subdir} | awk '{print $9}'`
    #echo ${listoffiles}
    #tr '\n' ' ' <<< "$listoffiles"
    indice_files=1
    outputname="${odir}${subdir}_${indice_job}.root"
    commandHadd="${outputname} "

    for file in ${listoffiles};do
        commandHadd+="${pathin2p3}${idir}/${subdir}/${file} "
        eucli=$((indice_files%filesperjobs))
        if [ ${eucli} -eq 0 ]; then
            let indice_job++
            #echo "hadd ${options} ${commandHadd}&"
            hadd ${options} ${commandHadd}&
            outputname="${odir}${subdir}_${indice_job}.root"
            commandHadd="${outputname} "
        fi
        let indice_files++
        
    done
    #echo "hadd ${options} ${commandHadd}&"
    hadd ${options} ${commandHadd}&

done

while [ `ps -u ${USER} | grep hadd | wc -l | awk '{print $1}'` != "0" ] ; do sleep 1s ; done
mkdir result
#echo "hadd -f ${odir}result/histo.root ${odir}00*.root"
hadd -f ${odir}result/histo.root ${odir}00*.root
rm ${odir}00*.root 
