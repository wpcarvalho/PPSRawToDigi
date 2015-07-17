import FWCore.ParameterSet.Config as cms

T1Digis = cms.EDProducer("T1DigiProducer",
    nScaBins = cms.int32(8),
    doNoise = cms.bool(False),
    doCrosstalk = cms.bool(False),
    stripSamplingTime = cms.double(25.0),
    stripSignalStartTime = cms.double(-250.0),
    stripSignalStopTime = cms.double(500.0),
    scaPeakBin = cms.int32(4),
    scaNoiseMode = cms.string('simple'),
    analogNoise = cms.double(2.7),
    pedestal = cms.double(600.0),
    pedestalWidth = cms.double(0.0),
    wireSignalStartTime = cms.double(-100.0),
    wireSignalStopTime = cms.double(150.0),
    wireSamplingTime = cms.double(5.0),
    wireTimingError = cms.double(0.0),    

    Electronics = cms.string('Santiard'),
    Verbosity = cms.int32(1),
    # thresholds and noise in mV
    THR1 = cms.double(500.0),
    THR2 = cms.double(1000.0),
    THR = cms.double(5.0),
    NOISE = cms.double(0.0),
    WIRETHR = cms.double(1000.0),
    WIRENOISE = cms.double(0.0)
)


