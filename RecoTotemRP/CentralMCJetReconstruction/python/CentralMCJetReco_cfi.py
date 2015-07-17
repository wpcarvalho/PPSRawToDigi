import FWCore.ParameterSet.Config as cms

CentralMCJetReco = cms.EDProducer("CentralMCJetReconstruction",
  hepmcInstance = cms.string(""),
  hepmcLabel = cms.string("generator"),
  jetOutputLabel = cms.string("kt6algorithm"),
  verbosity = cms.int32(5),
  R = cms.double(0.6),
  minEta = cms.double(-2.5),
  maxEta = cms.double(2.5),
  minPt = cms.double(36), #GeV
  jetAlgName = cms.string("kt_algorithm"),    #kt_algorithm, cambridge_algorithm, antikt_algorithm, genkt_algorithm, ee_kt_algorithm, ee_genkt_algorithm
  beamEnergy = cms.double(4000),  #GeV
  findJets = cms.bool(True)  #find jets or only gather information on diffractive system 
)

