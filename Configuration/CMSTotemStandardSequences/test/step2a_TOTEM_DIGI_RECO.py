import FWCore.ParameterSet.Config as cms

from Configuration.TotemStandardSequences.prodT1T2Default_cfg import *

process.setName_("RECOTOTEM")
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring("file:MinbiasTest_SIM.root")
   )

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

# Specify the output filename
exec 'process.' + str(process.outpath) + '.fileName = cms.untracked.string("file:MinbiasTest_RECO_TOTEM.root")'


# process.TrackerGeometricDetESModule = cms.ESProducer("TrackerGeometricDetESModule",
#     fromDDD = cms.bool(True)
#)

process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")

#Overwrite module chain
process.p1 = cms.Path(process.mix*process.T1Digis*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks2*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3*process.T2CC)



print process.dumpPython()