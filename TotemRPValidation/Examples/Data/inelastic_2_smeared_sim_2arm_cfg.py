import FWCore.ParameterSet.Config as cms

process = cms.Process("GunVal220")

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
)

process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.load("Configuration.TotemCommon.RandomNumbers_cfi")

process.load("SimGeneral.HepPDTESSource.pdt_cfi")

# Geometry - beta* specific
process.load("Configuration.TotemCommon.geometryRP_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Beta_2_220_RPs_only/RP_Dist_Beam_Cent.xml')
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_2_cfi")

# Magnetic Field
### Full field map, static configuration for each field value
# process.load("Configuration.StandardSequences.MagneticField_20T_cff")
# process.load("Configuration.StandardSequences.MagneticField_30T_cff")
# process.load("Configuration.StandardSequences.MagneticField_35T_cff")
process.load("Configuration.StandardSequences.MagneticField_38T_cff")
# process.load("Configuration.StandardSequences.MagneticField_40T_cff") 

# Monte Carlo gun
process.load("IOMC.FlatProtonLogKsiLogTGun.Beta2_cfi")
process.FlatProtonLogKsiLogTGun.Verbosity = 0

# Smearing
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")
process.SmearingGenerator.verbosity = 0

# Oscar - G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")

# redirect the reference:
# BeamProtTransportSetup (@ Configuration.TotemCommon.g4SimHits_cfi)
# --> BeamProtTransportSetup (@ Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_2_cfi)
process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup

process.g4SimHits.Generator.HepMCProductLabel = 'source'    # energy+vertex smearing

process.g4SimHits.Physics.BeamProtTransportSetup.Verbosity = False

process.load("SimGeneral.MixingModule.mixNoPU_cfi")

# RP Strip digitization
process.load("SimTotem.RPDigiProducer.RPSiDetConf_cfi")

process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
process.RPClustProd.Verbosity = 0

process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")
process.RPHecoHitProd.Verbosity = 0

process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")
process.RPSinglTrackCandFind.Verbosity = 0

process.load("RecoTotemRP.RPTrackCandidateCollectionFitter.RPSingleTrackCandCollFitted_cfi")
process.RPSingleTrackCandCollFit.Verbosity = 0

# 2sided reconstruction
process.load("RecoTotemRP.RPInelasticReconstruction.Rec_beta_2_220_2_Arm_cfi")

# redirect the reference:
# BeamProtTransportSetup (@ RecoTotemRP.RPInelasticReconstruction.Rec_beta_2_220_2_Arm_cfi)
# --> BeamProtTransportSetup (@ Configuration.TotemOpticsConfiguration.OpticsConfig_7000GeV_2_cfi)
process.RP220Reconst.BeamProtTransportSetup = process.BeamProtTransportSetup

process.RP220Reconst.Verbosity = 1

process.load("TotemRPValidation.InelasticReconstructionValidation.Inelastic2ArmRecValidation_2_cfi")
process.RPInelProtRecVal.HistogramFileName = 'RPInelasticRecValHistsBeta2DiffractiveMass.root'
process.RPInelProtRecVal.Verbosity = 1

process.HitDists = cms.EDAnalyzer("HitDistributions",
    verbosity = cms.untracked.uint32(0),
    outputFile = cms.string('hitDistributions_beta_2_2_sided_rec.root'),
    rpList = cms.vuint32(100, 101, 102, 103, 104, 
        105, 120, 121, 122, 123, 
        124, 125, 0, 1, 2, 
        3, 4, 5, 20, 21, 
        22, 23, 24, 25)
)

process.load("TotemRPValidation.RPReconstructedTracksValidation.ReconstTracksVal_beta_2_cfi")
process.RPRecTracksVal.Verbosity = 0

process.o1 = cms.OutputModule("PoolOutputModule",
    fileName = cms.untracked.string('inelastic_2_smeared_sim_2Arm.root')
)

process.p1 = cms.Path(process.SmearingGenerator*process.g4SimHits*process.mix*process.RPSiDetDigitizer*process.RPClustProd*process.RPHecoHitProd*process.RPSinglTrackCandFind*process.RPSingleTrackCandCollFit*process.HitDists*process.RP220Reconst*process.RPInelProtRecVal*process.RPRecTracksVal)

process.outpath = cms.EndPath(process.o1)

process.Tracer = cms.Service("Tracer")
