import FWCore.ParameterSet.Config as cms

process = cms.Process("prodRPinelasticBetaXXXEnergyYYYTeV")

process.load("Configuration.TotemCommon.LoggerMin_cfi")
process.load("Configuration.TotemCommon.RandomNumbers_cfi")
process.load("IOMC.FlatProtonLogKsiLogTGun.Beta90_cfi")
process.generator.Verbosity = 0
process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.verbosity = 1
# Magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.load("Configuration.TotemCommon.g4SimHits_cfi")
# optics
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")
################## STEP 4 process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit
process.load("Configuration.TotemStandardSequences/RP_Digi_and_TrackReconstruction_cfi")
process.RPHecoHitProd.Verbosity = 1
process.RPSinglTrackCandFind.Verbosity = 1
process.RPSingleTrackCandCollFit.Verbosity = 1
################## STEP 5 process.InelasticReconstruction

# reconstruction
process.load("RecoTotemRP.RPInelasticReconstruction.RPRec220_cfi")
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup

process.load("L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi")



process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *','drop *_*mix*_*_*', 'drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*'),
    fileName = cms.untracked.string('file:prodRPinelasticBetaXXXEnergyYYYTeV.root')
)
process.outpath = cms.EndPath(process.o1)
exec 'process.' + str(process.outpath) + '.fileName = cms.untracked.string("file:RPInelasticRecoTest_cfg.root")' 



process.source = cms.Source("EmptySource")

process.setName_("prodRPinelasticBeta11Energy1TeV")


process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90/RP_Dist_Beam_Cent.xml')


process.g4SimHits.Generator.HepMCProductLabel = 'generator'
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
process.g4SimHits.UseMagneticField = False




process.RPClustProd.DigiLabel = cms.InputTag("RPSiDetDigitizer")
process.RPHecoHitProd.ClusterLabel = cms.InputTag("RPClustProd")
process.RPSinglTrackCandFind.RPRecoHitLabel = cms.InputTag("RPHecoHitProd")
process.RPSingleTrackCandCollFit.RPTrackCandidateCollectionLabel = cms.InputTag("RPSinglTrackCandFind")
process.RP220Reconst.RPFittedTrackCollectionLabel = cms.InputTag("RPSingleTrackCandCollFit")

process.p1 = cms.Path(process.generator
                      *process.SmearingGenerator
                      *process.g4SimHits
                      *process.mix
                      *process.RPSiDetDigitizer
                      *process.RPClustProd
                      *process.RPHecoHitProd
                      *process.RPSinglTrackCandFind
                      *process.RPSingleTrackCandCollFit
                      *process.RP220Reconst
                      )



