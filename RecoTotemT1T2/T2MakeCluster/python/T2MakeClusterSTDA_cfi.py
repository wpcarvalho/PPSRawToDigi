import FWCore.ParameterSet.Config as cms

T2MClSTDA = cms.EDProducer("T2MakeCluster",

TakeCleanEventOnly=cms.bool(False),
maskvect = cms.vuint32(),
Projectionthreshold=cms.double(0.1),
BlobMinSize = cms.int32(5),
SimuClusterEfficiency=cms.bool(False),
EffiGeoRootFileName= cms.string("")
)
