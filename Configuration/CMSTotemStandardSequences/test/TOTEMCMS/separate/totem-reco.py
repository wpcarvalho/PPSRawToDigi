import FWCore.ParameterSet.Config as cms

from Configuration.TotemStandardSequences.prodRPT1T2Default_cfg import *

process.setName_("RECOTOTEM")
process.source = cms.Source("PoolSource",
        fileNames = cms.untracked.vstring("file:TOTEM-SIM.root")
   )

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
  input = cms.untracked.int32(-1)
)

# Specify the output filename
exec 'process.' + str(process.outpath) + '.fileName = cms.untracked.string("file:TOTEM-RECO.root")'

process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")
process.load("Configuration.TotemCommon.geometryRPT1T2CMS_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_6500GeV_90_cfi")


#Overwrite module chain
process.p1 = cms.Path(process.mix*process.T1Digis*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks2*process.T2Digis*process.T2MCl*process.T2Hits*process.T2RoadPadFinder*process.T2TrackColl3*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.RPCC*process.T2CC)



#print process.dumpPython()
