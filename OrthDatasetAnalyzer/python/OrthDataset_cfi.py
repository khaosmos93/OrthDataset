import FWCore.ParameterSet.Config as cms

OrthDataset = cms.EDAnalyzer('OrthDatasetAnalyzer',
                            Verbose = cms.bool(True),
                            MinMass = cms.double(900),
                        )
