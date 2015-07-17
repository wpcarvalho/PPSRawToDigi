import FWCore.ParameterSet.Config as cms

T2Hits = cms.EDProducer("T2RecHit",
    IncludeClass0Hits = cms.bool(True),
    InsertAlignmentbyCFG=cms.bool(False),
    InsertAlignmentbyDB=cms.bool(False), # osolete, this flag is no longer supported, because of removed TotemDatabaseService module                   
    Cl1MaxPad = cms.uint32(30),
    Cl1MaxStrip = cms.uint32(30),

    # #should be of size = 40
    DXdisp = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0),
    DYdisp = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0)
)


