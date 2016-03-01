# example config, copied from /afs/cern.ch/work/t/totemprd/totem/test_leszek/reco_run9510/workspace/totcsi_rec01/run1/parallel1/submission_files

import FWCore.ParameterSet.Config as cms

process = cms.Process("9998")
process.maxEvents = cms.untracked.PSet(
input = cms.untracked.int32(1000)
)

process.load('TotemRawData.Readers.RawDataSource_cfi')
process.source.fileNames = cms.untracked.vstring()

process.source.fileNames.append('/afs/cern.ch/work/p/polme/public/totemdata/run_9881_EVB13_1.000.srs')


process.load("Configuration.TotemCommon.LoggerMax_cfi")
process.load("Configuration.TotemCommon.geometryGlobal_cfi")


# TODO - check RP distances, check if commited to repository
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2015_08_25_fill4269_2/RP_Dist_Beam_Cent.xml")
# TODO - check RP alignement, check if commited to repository
process.load("TotemAlignment.RPDataFormats.TotemRPIncludeAlignments_cfi")
process.TotemRPIncludeAlignments.RealFiles = cms.vstring(
"TotemAlignment/RPData/alignment_tb_run9509_45_no20.xml",
"TotemAlignment/RPData/alignment_tb_run9509_56.xml"
)

process.load('TotemCondFormats.DAQInformation.DAQMappingSourceXML_cfi')
process.DAQMappingSourceXML.mappingFileNames.append('TotemCondFormats/DAQInformation/data/rp_220_210far.xml')

process.load("Configuration.TotemCommon.RandomNumbers_cfi")

process.load('TotemRawData.RawToDigi.Raw2DigiProducer_cfi')
process.Raw2DigiProducer.rpDataProductLabel = cms.untracked.string("")

process.load("RecoTotemRP.RPClusterizer.RPClusterizationConf_cfi")
process.RPClustProd.DigiLabel = cms.InputTag("Raw2DigiProducer")
process.load("RecoTotemRP.RPRecoHitProducer.RPRecoHitProdConf_cfi")
process.load("RecoTotemRP.RPSingleCandidateTrackFinder.RPSingleTrackCandFindConf_cfi")


process.output = cms.OutputModule(
"PoolOutputModule",
fileName = cms.untracked.string("file:./out.root"),
outputCommands = cms.untracked.vstring('keep *')
)

process.path = cms.Path(process.Raw2DigiProducer
*process.RPClustProd
*process.RPRecoHitProd
*process.RPSinglTrackCandFind
)

process.outpath = cms.EndPath(process.output)
