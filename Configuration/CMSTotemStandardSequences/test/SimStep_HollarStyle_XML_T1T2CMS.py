# Auto generated configuration file
# using: 
# Revision: 1.19 
# Source: /local/reps/CMSSW/CMSSW/Configuration/Applications/python/ConfigBuilder.py,v 
# with command line options: --filein file:MinbiasTest_GEN.root --fileout file:MinbiasTest_SIM.root --mc --eventcontent RAWSIM --customise Configuration/StandardSequences/SimWithCastor_cff.customise, SLHCUpgradeSimulations/Configuration/postLS1Customs.customisePostLS1 --datatier SIM --conditions POSTLS170_V5::All --step SIM --magField 38T_PostLS1 --geometry DBExtendedPostLS1 --python_filename SimStep_HollarStyle.py --no_exec -n 1
import FWCore.ParameterSet.Config as cms

process = cms.Process('SIM')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryRecoDB_cff')
#process.load('Configuration.StandardSequences.GeometrySimDB_cff')#merijn changed this to next statement. It seems to be needed to load things locally.
#process.load('Geometry.CMSCommonData.cmsExtendedGeometry2015XML_cfi') Load T1T2CMS  geo
process.load('Configuration.TotemCommon.geometryT1T2CMS_cfi')

process.load('Configuration.StandardSequences.MagneticField_38T_PostLS1_cff')
process.load('Configuration.StandardSequences.SimIdeal_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(20)
)

#merijn add:
process.RandomNumberGeneratorService.generator.initialSeed=456789
process.RandomNumberGeneratorService.g4SimHits.initialSeed=9876
process.RandomNumberGeneratorService.VtxSmeared.initialSeed=123456789

# Input source
process.source = cms.Source("PoolSource",
    secondaryFileNames = cms.untracked.vstring(),
    fileNames = cms.untracked.vstring('file:MinbiasTest_GEN.root')
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.19 $'),
    annotation = cms.untracked.string('--filein nevts:1'),
    name = cms.untracked.string('Applications')
)

# Output definition

process.RAWSIMoutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    eventAutoFlushCompressedSize = cms.untracked.int32(5242880),
    outputCommands = process.RAWSIMEventContent.outputCommands,
    fileName = cms.untracked.string('file:MinbiasTest_SIM.root'),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('SIM')
    )
)

# Additional output definition

# Other statements
#process.XMLFromDBSource.label = cms.string("ExtendedPostLS1") #merijn changed this to next statement. It seems to be needed to load things locally.
from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, 'POSTLS170_V5::All', '')

# Path and EndPath definitions
process.simulation_step = cms.Path(process.psim)
process.endjob_step = cms.EndPath(process.endOfProcess)
process.RAWSIMoutput_step = cms.EndPath(process.RAWSIMoutput)

# Schedule definition
process.schedule = cms.Schedule(process.simulation_step,process.endjob_step,process.RAWSIMoutput_step)

# customisation of the process.

# Automatic addition of the customisation function from Configuration.StandardSequences.SimWithCastor_cff
from Configuration.StandardSequences.SimWithCastor_cff import customise 

#call to customisation function customise imported from Configuration.StandardSequences.SimWithCastor_cff
process = customise(process)

# End of customisation functions
