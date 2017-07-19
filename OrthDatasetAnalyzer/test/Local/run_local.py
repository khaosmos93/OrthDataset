import FWCore.ParameterSet.Config as cms
process = cms.Process( "MSAnalyser" )

process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
                    #'file:/u/user/msoh/ZPrime/OrthDataset/CMSSW_8_0_28/src/OrthDataset/OrthDatasetAnalyzer/test/1CFA8097-8AEA-E611-980F-001E67E6F855.root', #/JetHT/Run2016H-03Feb2017_ver3-v1/MINIAOD
                    #'file:/u/user/msoh/ZPrime/OrthDataset/CMSSW_8_0_28/src/OrthDataset/OrthDatasetAnalyzer/test/16F28614-84EA-E611-8083-A0369F310374.root', #SingleMuon
                    #'file:/u/user/msoh/ZPrime/OrthDataset/CMSSW_8_0_28/src/OrthDataset/OrthDatasetAnalyzer/test/BAF2A4E0-5CEC-E611-BFBC-848F69FD2925.root',  #JetHT RunG
                    'file:/u/user/msoh/ZPrime/OrthDataset/CMSSW_8_0_28/src/OrthDataset/OrthDatasetAnalyzer/test/F4C1CC9E-C9EA-E611-93AA-001517FAD934.root',  #SingleMuon RunG
                    #'file:/u/user/msoh/ZPrime/OrthDataset/CMSSW_8_0_28/src/OrthDataset/OrthDatasetAnalyzer/test/6A4BD1CD-C7EB-E611-AADE-008CFA197CCC.root',   #exit code 139
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

import FWCore.PythonUtilities.LumiList as LumiList
#process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON_MuonPhys.txt').getVLuminosityBlockRange()
#process.source.lumisToProcess = LumiList.LumiList(filename = '/afs/cern.ch/cms/CAF/CMSCOMM/COMM_DQM/certification/Collisions16/13TeV/ReReco/Final/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()
process.source.lumisToProcess = LumiList.LumiList(filename = '/u/user/msoh/Samples/Data/Cert_271036-284044_13TeV_23Sep2016ReReco_Collisions16_JSON.txt').getVLuminosityBlockRange()


process.load("FWCore.MessageService.MessageLogger_cfi")
process.MessageLogger.cerr.INFO.limit = 0
process.MessageLogger.cout.threshold = cms.untracked.string('WARNING')
process.MessageLogger.cerr.FwkSummary = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1000),
    limit = cms.untracked.int32(10000000)
)
process.MessageLogger.cerr.FwkReport = cms.untracked.PSet(
    reportEvery = cms.untracked.int32(1000),
    limit = cms.untracked.int32(10000000)
)

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( -1 )
#   input = cms.untracked.int32( 20 )
#   input = cms.untracked.int32( 100 )
#   input = cms.untracked.int32( 2000 )
#   input = cms.untracked.int32( 50000 )
#   input = cms.untracked.int32( 100000 )
#   input = cms.untracked.int32( 300000 )
)

process.OrthDataset = cms.EDAnalyzer('OrthDatasetAnalyzer',
                            Verbose = cms.bool(True),
                            #Verbose = cms.bool(False),
                            MinMass = cms.double(0),
                        )


#### Standard Configurations
#process.load('Configuration.StandardSequences.Services_cff')
#process.load('Configuration.StandardSequences.Geometry_cff')
#process.load('Configuration.StandardSequences.Reconstruction_cff')
#process.load('Configuration.StandardSequences.MagneticField_cff')
#
#### conditions
GT='80X_dataRun2_2016SeptRepro_v7'
#GT='80X_dataRun2_Prompt_v16'
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_condDBv2_cff")
process.GlobalTag.globaltag = GT
#process.GlobalTag.globaltag = '80X_dataRun2_2016LegacyRepro_v3' #Data 80X
#process.GlobalTag.globaltag = '80X_dataRun2_Prompt_v16'  #Feb RunHv3
##from Configuration.AlCa.GlobalTag import GlobalTag
##process.GlobalTag = GlobalTag(process.GlobalTag, 'auto:run2_mc', '')

process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
process.load('Configuration.StandardSequences.MagneticField_AutoFromDBCurrent_cff')
process.load('TrackingTools.TransientTrack.TransientTrackBuilder_cfi')

PD = 'SingleMuon'
Period = 'MaybeG'
OUTPUT='OrthDatasetTree_'+ PD +'_' + Period +'.root'
process.TFileService = cms.Service("TFileService",
                                       fileName = cms.string(OUTPUT)
                                   )

process.p = cms.Path(process.OrthDataset)

process.options = cms.untracked.PSet(
    SkipEvent = cms.untracked.vstring('ProductNotFound'),
    numberOfThreads = cms.untracked.uint32(8)
)
