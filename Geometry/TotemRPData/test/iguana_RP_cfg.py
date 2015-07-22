import FWCore.ParameterSet.Config as cms

process = cms.Process("IGUANA")

#from FWCore.MessageLogger.MessageLogger_cfi import *

process.MessageLogger = cms.Service( "MessageLogger",
  debugModules = cms.untracked.vstring( '*' ),
  debug = cms.untracked.PSet(
    threshold = cms.untracked.string( 'DEBUG' )
  ),
  cerr = cms.untracked.PSet(
    threshold = cms.untracked.string( 'ERROR' )
  ),
  cout = cms.untracked.PSet(
    threshold = cms.untracked.string( 'INFO' )
  ),
  destinations = cms.untracked.vstring( 'cout', 'cerr', 'debug' )
)

process.load("Geometry.TotemRPData.Geom_RP_150_220_90_cfi")

process.VisConfigurationService = cms.Service("VisConfigurationService",
                                              Views = cms.untracked.vstring('3D Window'),
                                              ContentProxies = cms.untracked.vstring('Simulation/Core', 
                                                                                     'Simulation/Geometry', 
                                                                                     'Reco/CMS Magnetic Field')
                                              )

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1)
    )
process.source = cms.Source("EmptySource")

process.prod = cms.EDProducer("GeometryProducer",
                              MagneticField = cms.PSet(
    delta = cms.double(1.0)
    ),
                              UseMagneticField = cms.bool(False),
                              UseSensitiveDetectors = cms.bool(False)
                              )

process.p = cms.Path(process.prod)
