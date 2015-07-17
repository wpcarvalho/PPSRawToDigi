import FWCore.ParameterSet.Config as cms

process = cms.Process("NTUPLETOTEM")

process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring("file:MinbiasTest_RECO_TOTEM.root"),
        noEventSort = cms.untracked.bool(True)
   )

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")
process.load("Configuration.TotemCommon.geometryRPT1T2CMS_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")

process.load("TotemRPValidation.HitDistributions.HitDistributions_cfi")
process.RPHitDists.outputFile = cms.string('file:MinbiasTest_RPHITS_TOTEM.root')

process.load("TotemAnalysis.TotemNtuplizer.TotemNtuplizer_cfi")
process.TotemNtuplizer.outputFileName = 'file:MinbiasTest_NTUPLE_TOTEM.root'
process.TotemNtuplizer.verbosity = 2
process.TotemNtuplizer.ProductLabelSimu = 'RPCCRawBits'
process.TotemNtuplizer.ModulLabelSimu = 'RPDataCCProducer'
process.TotemNtuplizer.T1TrackLabel = 't1tracks2'

process.p1 = cms.Path(process.RPHitDists*process.TotemNtuplizer)

#print process.dumpPython()
