import FWCore.ParameterSet.Config as cms

process = cms.Process("T2DigitizerTest")

process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring('file:simTotem.root'),
)

process.load("Configuration.TotemCommon.RandomNumbers_cfi")

process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.load("SimGeneral.MixingModule.mixNoPU_cfi")

process.load("SimG4CMS.Forward.test.testGeometryXML_cfi")

process.load("SimGeneral.HepPDTESSource.pythiapdt_cfi")

process.printf = cms.OutputModule("AsciiOutputModule")

process.load("SimTotem.T2Digitizer.T2Digis_cfi")

process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('t2digis.root')
)

process.p1 = cms.Path(process.mix*process.T2Digis)

process.outpath = cms.EndPath(process.o1)

