import FWCore.ParameterSet.Config as cms

t1analreco3 = cms.EDAnalyzer('T1Analysis_reco3',
                              Verbosity = cms.int32(0),
                             OutputFile = cms.string('t1analysis_reco3_4.3_4.5.root'),
                              DumpOnly0TrackEvents = cms.bool(False),
                              MaxTracks = cms.int32(1000),
                              ChiProbCut=cms.double(0),
                              TrackLabel = cms.string("t1tracks3"),
                              EtaMin = cms.double(4.3),
                              EtaMax = cms.double(4.5),
			      RawEventLabel = cms.InputTag("source")

)
