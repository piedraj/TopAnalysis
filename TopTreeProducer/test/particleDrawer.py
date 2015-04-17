import FWCore.ParameterSet.Config as cms

process = cms.Process("particleDrawer")

process.load("FWCore.MessageService.MessageLogger_cfi")
process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )

process.source = cms.Source("PoolSource",
                            fileNames = cms.untracked.vstring('file:/afs/cern.ch/user/p/piedra/work/store/mc/Phys14DR/TTJets_MSDecaysCKM_central_Tune4C_13TeV-madgraph-tauola/MINIAODSIM/PU20bx25_PHYS14_25_V1-v1/00000/00C90EFC-3074-E411-A845-002590DB9262.root'),
                            skipEvents = cms.untracked.uint32(4)
                            )

process.printTree1 = cms.EDAnalyzer("ParticleListDrawer",
                                    src = cms.InputTag("prunedGenParticles"),
                                    maxEventsToPrint  = cms.untracked.int32(5)
                                    )

process.printTree2 = cms.EDAnalyzer("ParticleTreeDrawer",
                                    src           = cms.InputTag("prunedGenParticles"),
                                    printP4       = cms.untracked.bool(False),
                                    printPtEtaPhi = cms.untracked.bool(False),
                                    printVertex   = cms.untracked.bool(False),
                                    printStatus   = cms.untracked.bool(False),
                                    printIndex    = cms.untracked.bool(False)
                                    )

process.printEventNumber = cms.OutputModule("AsciiOutputModule")

process.p = cms.Path(process.printTree1*process.printTree2)
process.outpath = cms.EndPath(process.printEventNumber)
