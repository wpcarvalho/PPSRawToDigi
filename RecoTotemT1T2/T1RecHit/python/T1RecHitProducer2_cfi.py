import FWCore.ParameterSet.Config as cms

T1RecHit2 = cms.EDProducer('T1RecHitProducer2',
                           t1DigiModuleLabel = cms.string('Raw2DigiProducer'),
                           t1DigiProductLabel = cms.string('t1DataOutput'),
                           geometryFile = cms.string('Geometry/TotemGeometry/data/T1_data_geometry.dat'),
                           geomParsFile = cms.string(''),
                           verbosity = cms.untracked.uint32(1),
                           printProgressFrequency = cms.untracked.uint32(1000),
                           tripleTolerance = cms.double(1.)
                           )

