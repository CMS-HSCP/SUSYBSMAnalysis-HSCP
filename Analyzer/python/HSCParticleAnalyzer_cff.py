import FWCore.ParameterSet.Config as cms
import os

from SUSYBSMAnalysis.Analyzer.HSCParticleAnalyzer_cfi import HSCParticleAnalyzer as analyzer 
PATH_TO_DATA = "{}/src/SUSYBSMAnalysis/HSCP/data".format(os.getenv('CMSSW_BASE'))

''''
analyzer_2016 = analyzer.clone(
    DeDxSF_0        = cms.untracked.double(1.00000), #=1 if data
    DeDxSF_1        = cms.untracked.double(1.0325), #2017 data
    DeDxK           = cms.untracked.double(2.30),
    DeDxC           = cms.untracked.double(3.17),
)

analyzer_2017 = analyzer.clone(
    DeDxSF_0        = cms.untracked.double(1.00000), #=1 if data
    DeDxSF_1        = cms.untracked.double(1.0325), #2017 data
    DeDxK           = cms.untracked.double(2.30),
    DeDxC           = cms.untracked.double(3.17),
    DeDxTemplate    = cms.untracked.string("{}/template_2017B.root".format(PATH_TO_DATA))
)

analyzer_2018 = analyzer.Analyzer.HSCParticleAnalyzer_cfi.HSCParticleAnalyzer.clone(
    DeDxSF_0        = cms.untracked.double(1.00000), #=1 if data
    DeDxSF_1        = cms.untracked.double(1.0817), #2017 data
    DeDxK           = cms.untracked.double(2.27),
    DeDxC           = cms.untracked.double(3.16),
    DeDxTemplate    = cms.untracked.string("{}/template_2017B.root".format(PATH_TO_DATA))
)

from Configuration.Eras.Era_Run2_2016_cff import Run2_2016
from Configuration.Eras.Era_Run2_2017_cff import Run2_2017
from Configuration.Eras.Era_Run2_2018_cff import Run2_2018

Run2_2016.toReplaceWith(analyzer, analyzer_2016)
Run2_2017.toReplaceWith(analyzer, analyzer_2017)
Run2_2018.toReplaceWith(analyzer, analyzer_2018)
'''

