import FWCore.ParameterSet.Config as cms

T1Digis = cms.EDProducer("T1DigiProducer",
    # thresholds and noise in mV
    Verbosity = cms.int32(0),
    Electronics = cms.string('VFAT'),
    THR1 = cms.double(80.0),
    THR2 = cms.double(80.0),
    THR = cms.double(50.0),
    NOISE = cms.double(0.0),
    WIRETHR = cms.double(60.0),
    WIRENOISE = cms.double(0.0),
    simulateDeadChannels = cms.bool(False)

)


