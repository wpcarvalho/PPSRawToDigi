import FWCore.ParameterSet.Config as cms
process = cms.Process("recoTB")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

process.load("Configuration.TotemCommon.LoggerMin_cfi")

# raw data source
process.source = cms.Source("RawDataSource",
        verbosity = cms.untracked.uint32(0),
        eventsToCheck = cms.uint32(10),
    skipCorruptedEvents = cms.bool(False),
                            performChecks = cms.bool(True),
    fileNames = cms.vstring(
      'TotemRawData/RawToDigi/python/RawData/runT1_091101_009.vme2')
)


process.DAQInformationSourceXML = cms.ESSource("DAQInformationSourceXML",
    xmlFileName = cms.string('TotemRawData/RawToDigi/python/T1_neg_half_2.xml')
)

process.t1raw2digi = cms.EDProducer("T1XMLDataDigiProducer",
                                    verbosity = cms.untracked.uint32(0),
                                    GET_NOISY_MAP = cms.untracked.uint32(1),
                                    NOISY_MAP = cms.untracked.string("NOISY_MAP.txt") 
)




process.load("RecoTotemT1T2.T1MakeCluster.T1MakeCluster_cfi")

process.load("RecoTotemT1T2.T1RecHit.T1RecHit_cfi")

process.t1recoanal = cms.EDAnalyzer("T1RecHitAnalyzer",
                                    OutputFile = cms.string('./hhh_mu_3.root') )

process.t1roads = cms.EDProducer("T1RoadProducer",
    Verbosity = cms.int32(0),
    DeltaEta = cms.double(.9),
    DeltaPhi = cms.double(.3),
    MinHits = cms.int32(3),
    MaxHits = cms.int32(50)
)

process.t1tracks = cms.EDProducer("T1TrackProducer",
    Verbosity = cms.int32(0),
    ChiRidThr = cms.double(40.0)
)

process.t1tracksanalTB = cms.EDAnalyzer("T1TrackAnalyzerTB",
    Verbosity = cms.int32(0),
    SeeTrack = cms.double(1.0),
    SeeHits = cms.double(0),
    ChiOverNCut = cms.double(20.0),
    ZRange = cms.double(1000.0),
    Zmin = cms.double(-10000.0),
    Zmax = cms.double(10000.0),
    realBeamPosX = cms.double(0.0),
    realBeamPosY = cms.double(0.0),
    realBeamAngle = cms.double(0.0)
)


process.t1valid = cms.EDAnalyzer("T1Validation",
    Verbosity = cms.int32(0),
    SIM = cms.double(0.0),
    DIGI = cms.double(1.0),
    RECO = cms.double(1.0),
    TrackLabel = cms.string('t1tracks'),
    OutputFile = cms.string('./T1Valid_091101_mux.root')
)


process.primaryvertex = cms.EDFilter("PrimaryVertexFinder",
SeeHits = cms.double(0.0),
    Verbosity = cms.int32(0),
    ChiOverNCut = cms.double(30.0),
    TrackRmin = cms.double(1000.0),
    Telescopes = cms.int32(1),
    SeeTrack = cms.double(0.0),
    ZRange = cms.double(10000.0)
                                     
                                      )
process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('./recoT1_091101_mux.root'),
                              outputCommands=cms.untracked.vstring('drop totemRawEvent_*_*_*')
                              )
process.p1 = cms.Path(process.t1raw2digi*process.t1cluster*process.t1rechit*process.t1roads*process.t1tracks*process.t1valid)

















                              
#process.p1 = cms.Path(process.RawToDigi)

process.outpath = cms.EndPath(process.o1)


