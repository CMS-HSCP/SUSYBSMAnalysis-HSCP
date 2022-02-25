# HSCP_internship2019

note: commit of 4th of October 2019
code working under 10_6_2

in SUSYBSMAnalysis/HSCP  (from https://github.com/CMS-HSCP/SUSYBSMAnalysis-HSCP.git a long long time ago)
modification of SUSYBSMAnalysis/HSCP/src/BetaCalculatorECAL.cc
by https://github.com/cms-sw/cmssw/blob/CMSSW_10_6_X/SUSYBSMAnalysis/HSCP/src/BetaCalculatorECAL.cc
in order to compile

in principle now the ntuple.cc code should run on miniAOD (without the HSCPCandidate part and with reduced info), on AOD, on data and MC.
or at least it worked sometimes during the development phase.



note on 25th of March 2020
/plugin/ntuple.cc is the code producing the ntuple in the standard conditions, on AOD or miniAOD
/plugin/calib_ntuple.cc is the code producing the ntuple for alca-reco samples (reduced information)

