import FWCore.ParameterSet.Config as cms

T2Hits = cms.EDProducer("T2RecHit",
    IncludeClass0Hits = cms.bool(True),
    InsertAlignmentbyCFG=cms.bool(False),
    InsertAlignmentbyDB=cms.bool(False), # osolete, this flag is no longer supported, because of removed TotemDatabaseService module                   
    Cl1MaxPad = cms.uint32(25),
    Cl1MaxStrip = cms.uint32(25),
    inputFileNameMisal= cms.untracked.string('/home/mirko/SL/WorkingArea/CMSSW311_TEST/CMSSW_3_1_1/src/TotemAlignment/T2TrkBasedInternalAlignment/test/MisalOut/RecoRZ_Merged_OptorRX1_4And6clk.dat'),
    useTXTfile=cms.bool(False),
    CorrectWithResolution=cms.bool(False),                    
    verbosity=cms.untracked.bool(False),                    
    # #should be of size = 40
    DXdisp = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0),
    DYdisp = cms.vdouble(0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0, 
        0.0, 0.0, 0.0, 0.0, 0.0),
 
    LabelproductInstanceName=cms.string("T2Hits"),
    ModuleLabelInput=cms.string("T2MCl")
)



