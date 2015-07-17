import FWCore.ParameterSet.Config as cms

process = cms.Process("prodRPelasticBeta90Energy7TeV")

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including: errors.log, warnings.log, infos.log, debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
process.load("Configuration.TotemCommon.LoggerMin_cfi")

# Use random number generator service
process.load("Configuration.TotemCommon.RandomNumbers_cfi")

import IOMC.Elegent.ElegentSource_cfi

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")


# Magnetic Field, by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")

# optics
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")

process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")

# G4 geometry
process.load("Configuration.TotemCommon.geometryRP_cfi")

# G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")

# Use particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

# No pile up for the mixing module
process.load("Configuration.TotemCommon.mixNoPU_cfi")

# RP Strip digitization
process.load("SimTotem.RPDigiProducer.RPSiDetConf_cfi")
process.RPSiDetDigitizer.RPVerbosity = 0



process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
#process.RPClustProd.Verbosity = 1

process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")
#process.RPHecoHitProd.Verbosity = 1

#
process.load("RecoTotemRP.RPNonParallelTrackCandidateFinder.RPNonParallelTrackCandidateFinder_cfi")
process.NonParallelTrackFinder.verbosity = 10

#
process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")
#process.RPSinglTrackCandFind.Verbosity = 1
# 
process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")

#
process.load("RecoTotemRP.RPMulCandidateTrackFinder.RPMulTrackCandFindConf_cfi")
process.RPMulTrackCandFind.Verbosity = 1
#
process.load("RecoTotemRP.RPMulTrackCandidateCollectionFitter.RPMulTrackCandCollFitter_cfi")

# # Elastic reconstruction options
process.load("RecoTotemRP.RPElasticReconstruction.rper_7000GeV_90_cfi")


process.load("L1TriggerTotem.CoincidenceChip.RPCoincidenceProducer_cfi")

process.load("RecoTotemRP.CentralMCJetReconstruction.CentralMCJetReco_cfi")

# Specify the maximum events to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(10)
)

# Configure the output module (save the result in a file)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *','drop *_*mix*_*_*', 'drop *_*_*TrackerHits*_*', 'drop *_*_*Muon*_*', 'drop *_*_*Ecal*_*', 'drop *_*_*Hcal*_*', 'drop *_*_*Calo*_*', 'drop *_*_*Castor*_*', 'drop *_*_*FP420SI_*', 'drop *_*_*ZDCHITS_*', 'drop *_*_*BSCHits_*', 'drop *_*_*ChamberHits_*', 'drop *_*_*FibreHits_*', 'drop *_*_*WedgeHits_*'),
    fileName = cms.untracked.string('file:RPElasticRecoTest_cfg.root')
)

process.outpath = cms.EndPath(process.o1)




################## STEP 1
process.source = cms.Source("EmptySource")

################## STEP 2 - process.generator



# Monte Carlo gun - elastic specific
energy = "7000"
process.generator = IOMC.Elegent.ElegentSource_cfi.generator
process.generator.fileName = IOMC.Elegent.ElegentSource_cfi.ElegentDefaultFileName(energy)

################## STEP 3 process.SmearingGenerator

# declare optics parameters
# process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_90_cfi")



################## STEP 4 process.OptInfo

process.OptInfo = cms.EDAnalyzer("OpticsInformation")

################## STEP 5 process.*process.g4SimHits

# Geometry - beta* specific
# process.load("Configuration.TotemCommon.geometryRP_cfi")
# process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90_150_out/RP_Dist_Beam_Cent.xml')


  # The input source for G4 module is connected to "process.source".

################## STEP 6 #process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit


#process.RPSingleTrackCandCollFit.Verbosity = 1

################## STEP 6'
#process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPMulTrackCandFind*process.RPMulTrackCandCollFit

################## STEP 7 process.ElasticReconstruction

# reconstruction options
# process.load("RecoTotemRP.RPElasticReconstruction.rper_7000GeV_90_cfi")

################## STEP 8 process.RPCC







# particle generator paramteres
process.generator.t_min = '2.5E-2'  # beta* specific


process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_90_150_out/RP_Dist_Beam_Cent.xml')

process.g4SimHits.Generator.HepMCProductLabel = 'generator'  
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
process.g4SimHits.UseMagneticField = False

process.RPClustProd.DigiLabel = cms.InputTag("RPSiDetDigitizer")
process.RPHecoHitProd.ClusterLabel = cms.InputTag("RPClustProd")

process.NonParallelTrackFinder.DetSetVectorRPRecoHitLabel = cms.InputTag("RPHecoHitProd")

process.RPSinglTrackCandFind.RPRecoHitLabel = cms.InputTag("RPHecoHitProd")
process.RPSingleTrackCandCollFit.RPTrackCandidateCollectionLabel = cms.InputTag("RPSinglTrackCandFind")

process.RPMulTrackCandFind.RPRecoHitDetSetLabel = cms.InputTag("RPHecoHitProd")
process.RPMulTrackCandCollFit.RPMulTrackCandidateCollectionLabel = cms.InputTag("RPMulTrackCandFind")

process.ElasticReconstruction.RPFittedTrackCollLabel = cms.InputTag("RPSingleTrackCandCollFit")

process.RPCC.DetTriggerLabel = cms.InputTag("RPSiDetDigitizer")

process.p1 = cms.Path(process.generator
                      *process.SmearingGenerator
                      *process.OptInfo
                      *process.g4SimHits
                      *process.mix
                      *process.RPSiDetDigitizer
                      *process.RPClustProd
                      *process.RPHecoHitProd
		      *process.NonParallelTrackFinder
                      *process.RPSinglTrackCandFind
                      *process.RPSingleTrackCandCollFit
		      *process.RPMulTrackCandFind
		      *process.RPMulTrackCandCollFit
                      *process.ElasticReconstruction
                      *process.RPCC
                      *process.CentralMCJetReco
                      )




