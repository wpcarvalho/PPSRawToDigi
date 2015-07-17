import FWCore.ParameterSet.Config as cms
#RPElasticbeta90energy3p5TeV220October2011a3
process = cms.Process("RPElasticbeta90energy3p5TeV220October2011a3")

# process.load("Configuration.TotemCommon.LoggerMax_cfi")

process.load("Configuration.TotemCommon.RandomNumbers_cfi")

# Load particle table
process.load("SimGeneral.HepPDTESSource.pdt_cfi")

import IOMC.Elegent.ElegentSource_cfi

# Smearing: vertex+energy
process.load("IOMC.SmearingGenerator.SmearingGenerator_cfi")

# process.load("Configuration.StandardSequences.Geometry_cff")  

process.load("Configuration.TotemCommon.geometryRP_real_cfi")




process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")

# declare optics parameters
process.load("Configuration.TotemOpticsConfiguration.OpticsConfig_3500GeV_90_cfi")

# Magnetic Field
# by default we have 3.8T
process.load("Configuration.StandardSequences.MagneticField_cff")


# Oscar - G4 simulation & proton transport
process.load("Configuration.TotemCommon.g4SimHits_cfi")

# No pile up for the mixing module
process.load("Configuration.TotemCommon.mixNoPU_cfi")


process.load("SimTotem.RPDigiProducer.RPSiDetConf_cfi")


process.load("SimTransport.HectorProducer.HectorTransport_cfi")

########################### COMMON PART ##########################################

# Specify the maximum event to simulate
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
)

# minimum of logs
# process.load("Configuration.TotemCommon.LoggerMin_cfi")



# Configure the output module (save the result in a file -- RPT1T2phojetDPEbeta3.5energy3.5TeV.root)
process.o1 = cms.OutputModule("PoolOutputModule",
    outputCommands = cms.untracked.vstring('keep *',
        'drop *_*_TrackerHits*_*',
        'drop *_*_Muon*_*',
        'drop *_*_Ecal*_*',
        'drop *_*_Hcal*_*',
        'drop *_*_Calo*_*',
        'drop *_*_Castor*_*',
        'drop *_*_FP420SI_*',
        'drop *_*_ZDCHITS_*'),
    fileName = cms.untracked.string('RPElastic.root')
)

# Configure if you want to detail or simple log information.
# LoggerMax -- detail log info output including:
#   - errors.log
#   - warnings.log
#   - infos.log
#   - debugs.log
# LoggerMin -- simple log info output to the standard output (e.g. screen)
#process.load("Configuration.TotemCommon.LoggerMax_cfi")
#process.load("Configuration.TotemCommon.LoggerMin_cfi")


########################### SIMULATION ##########################################



process.source = cms.Source("EmptySource")


############## Generator
# Use random number generator service


energy = "3500"

process.generator = IOMC.Elegent.ElegentSource_cfi.generator
process.generator.fileName = IOMC.Elegent.ElegentSource_cfi.ElegentDefaultFileName(energy)
process.generator.t_min = '0.001'  # beta* specific
process.generator.t_max = '10'  # beta* specific


#process.SmearingGenerator.originalLabel = 'original'
#process.SmearingGenerator.modifyLabel = 'source'
# process.SmearingGenerator.verbosity = 1

# Geometry - beta* specific
# this file contains geometry of T1T2+RP



# real geometry
# TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule",
#     verbosity = cms.untracked.uint32(1)
# )



toberemoved = []
for xmlfile in process.XMLIdealGeometryESSource.geomXMLFiles:
    if xmlfile.endswith("RP_Dist_Beam_Cent.xml"):
       toberemoved.append(xmlfile)
for xmlfile in toberemoved:
    process.XMLIdealGeometryESSource.geomXMLFiles.remove(xmlfile)
  
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_10_20_3/RP_Dist_Beam_Cent.xml")


# # process.XMLIdealGeometryESSource.geomXMLFiles.append('Geometry/TotemRPData/data/RP_Param_Beam_Region.xml')
# totemGeomXMLFiles = cms.vstring(
# #         'Geometry/CMSCommonData/data/materials.xml', 
# #         'Geometry/CMSCommonData/data/rotations.xml', 
# #         'Geometry/CMSCommonData/data/normal/cmsextent.xml', 
# #         'Geometry/CMSCommonData/data/cms.xml', 
# #         'Geometry/CMSCommonData/data/beampipe.xml', 
# #         'Geometry/CMSCommonData/data/cmsBeam.xml', 
# #         'Geometry/CMSCommonData/data/cmsMother.xml', 
# #         'Geometry/CMSCommonData/data/mgnt.xml', 
# #         'Geometry/ForwardCommonData/data/forward.xml', 
# #         'Geometry/ForwardCommonData/data/totemRotations.xml', 
# #         'Geometry/ForwardCommonData/data/totemMaterials.xml', 
# #         'Geometry/ForwardCommonData/data/totemt1.xml', 
# #         'Geometry/ForwardCommonData/data/totemt2.xml', 
# #         'Geometry/ForwardCommonData/data/ionpump.xml', 
#         'Geometry/TotemRPData/data/RP_Box.xml', 
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_000.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_001.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_002.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_003.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_004.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_005.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_020.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_021.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_022.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_023.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_024.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_025.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_100.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_101.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_102.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_103.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_104.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_105.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_120.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_121.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_122.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_123.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_124.xml',
#         'Geometry/TotemRPData/data/RP_Box/RP_Box_125.xml',
#         'Geometry/TotemRPData/data/RP_Hybrid.xml', 
#         'Geometry/TotemRPData/data/RP_Materials.xml', 
#         'Geometry/TotemRPData/data/RP_Transformations.xml', 
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly.xml', 
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_000.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_001.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_002.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_003.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_004.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_005.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_020.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_021.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_022.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_023.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_024.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_025.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_100.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_101.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_102.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_103.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_104.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_105.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_120.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_121.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_122.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_123.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_124.xml',
#         'Geometry/TotemRPData/data/RP_Detectors_Assembly/RP_Detectors_Assembly_125.xml',
#         'Geometry/TotemRPData/data/RP_Device.xml', 
#         'Geometry/TotemRPData/data/RP_Vertical_Device.xml', 
#         'Geometry/TotemRPData/data/RP_Horizontal_Device.xml', 
#         'Geometry/TotemRPData/data/RP_220_Right_Station.xml', 
#         'Geometry/TotemRPData/data/RP_220_Left_Station.xml', 
#         'Geometry/TotemRPData/data/RP_147_Right_Station.xml', 
#         'Geometry/TotemRPData/data/RP_147_Left_Station.xml', 
#         'Geometry/TotemRPData/data/RP_Stations_Assembly.xml', 
#         'Geometry/TotemRPData/data/RP_Sensitive_Dets.xml', 
#         'Geometry/TotemRPData/data/RP_Cuts_Per_Region.xml', 
#         'Geometry/TotemRPData/data/TotemRPGlobal.xml', 
#         'Geometry/TotemRPData/data/RP_Param_Beam_Region.xml')           
# process.XMLIdealGeometryESSource.geomXMLFiles += totemGeomXMLFiles
# process.XMLIdealGeometryESSource.rootNodeName = cms.string('TotemRPGlobal:OTOTEM')


# geo_prefere = cms.ESPrefer("XMLIdealGeometryESSource", "XMLIdealGeometryESSource")

#this is configured in geometryRP_real_cfi.
# process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")


#  
# import TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi
#   
# process.TotemRPIncludeAlignments.MisalignedFiles = cms.vstring('TotemAlignment/RPData/LHC/2011_10_20_3/sr+el/45_220.xml',
#         'TotemAlignment/RPData/LHC/2011_10_20_3/sr+el/56_220.xml') 
# process.TotemRPIncludeAlignments.MeasuredFiles = cms.vstring()
#  
#  








process.g4SimHits.Physics.BeamProtTransportSetup = process.BeamProtTransportSetup
#
# NOTE:
#   In CMSSW_3_1_1, the Gun ("Pythia") is now created as a EDProducer and labeled
# as "process.generator". So the data source for Geant 4 simulation module (g4SimHits.Generator.HepMCProductLabel)
# should now be connected to "process.generator".
# -----------------------------------------------------------------------------------
#   In CMSSW_1_7_7. "g4SimHits.Generator.HepMCProductLabel" is conntected to "process.source", because
# "Pythia" is directly used as the source and also labeled "generator.source".
#
#process.g4SimHits.Generator.HepMCProductLabel = 'source'    # energy+vertex smearing
process.g4SimHits.Physics.DefaultCutValue = 100.
# process.g4SimHits.Generator.ApplyPtCuts = False (ApplyPtCuts does not already exist)
process.g4SimHits.Generator.ApplyPCuts = False
process.g4SimHits.Generator.ApplyEtaCuts = False
process.g4SimHits.UseMagneticField = False
process.g4SimHits.Generator.LeaveOnlyScatteredProtons = True
#
# verbosity control
#
# process.g4SimHits.G4EventManagerVerbosity = 1
# process.g4SimHits.G4StackManagerVerbosity = 1
# process.g4SimHits.G4TrackingManagerVerbosity = 1
# process.g4SimHits.MagneticField.Verbosity = True
# process.g4SimHits.Physics.Verbosity = 1
# process.g4SimHits.Physics.BeamProtTransportSetup.Verbosity = True
# process.g4SimHits.Generator.Verbosity = 1
# process.g4SimHits.SteppingAction.Verbosity = 1
# process.g4SimHits.Totem_RP_SD.Verbosity = 1
# process.g4SimHits.TotemSD.Verbosity = 1



########################### DIGI+RECO RP ##########################################


# process.RPSiDetDigitizer.RPVerbosity = 1

process.OptInfo = cms.EDAnalyzer("OpticsInformation")

process.SmearingGenerator.verbosity = 7

process.p1 = cms.Path(
    process.generator
    *process.SmearingGenerator
    *process.OptInfo
    *process.g4SimHits
    *process.mix
    *process.RPSiDetDigitizer
    )

process.outpath = cms.EndPath(process.o1)



