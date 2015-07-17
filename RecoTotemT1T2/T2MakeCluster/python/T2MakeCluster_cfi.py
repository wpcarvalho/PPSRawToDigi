import FWCore.ParameterSet.Config as cms

T2MCl = cms.EDProducer("T2MakeCluster",

TakeCleanEventOnly=cms.bool(False),
maskvect = cms.vuint32(),
Projectionthreshold=cms.double(0.1),
BlobMinSize = cms.int32(5),
SimuClusterEfficiency=cms.bool(False),
EffiGeoRootFileName= cms.string(""),
T2PadDigiCollectionLabel = cms.InputTag("T2Digis", "T2PadDigi"),
T2StripDigiCollectionLabel = cms.InputTag("T2Digis", "T2StripDigi")
)

#OLD METHOD
# mydet.arm()*20+mydet.plane()*4+mydet.halfTelescope()*2+mydet.planeSide()

#0 1 2 3 4 -->  0 4 8 12 16 (+0/2)
#	   -->  2 6 10 14 18 
#	   -->  1 5 9 13 17
#	   -->  3 7 11 15 19


# NEW METHOD
# pl*2+pls+  ht*10+  20*arm;
