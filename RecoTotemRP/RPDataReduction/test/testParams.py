import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(100) )


# geometry
process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_10_20_1/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")


process.load("Configuration/TotemOpticsConfiguration/OpticsConfig_3500GeV_90_cfi")

process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use

    fileNames = cms.untracked.vstring(
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6883.000_recoNewT1codeTry2.root'  
            'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6883/run_6883.000-1.root'
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_5657.000.vmea_.root',
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/RPT1T2_CommonReco_run5657/run_run_5657.012.vmea_.root'
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.000-1.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.001-2.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.002-3.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.003-4.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.004-5.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.005-6.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.006-7.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.007-8.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.008-9.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.009-10.root',
#           'rfio:/castor/cern.ch/totem/offline/Analysis/2011/T2RP/inelastic_reco/ver3.2/verguilo/6883/run_6883.020-21.root'

    )
)

process.load("RecoTotemRP.RPDataReduction.RPDataReduction_cfi")

process.demo.verb = cms.uint32(11)

process.demo.readMultiTrk = cms.bool(False)
process.demo.readT2 = cms.bool(False)
process.demo.readT1 = cms.bool(False)
process.demo.readLoNeg = cms.bool(False)
process.demo.readClusters = cms.bool(False)
process.demo.readRecoProt = cms.bool(False)


process.demo.tmcModule=cms.string('Raw2DigiProducer')
process.demo.tmcProd=cms.string('rpCCOutput')
process.demo.whichDiag=cms.uint32(0)

process.demo.fileName = cms.string('newReco-r6883-testAllReadParamsFalse-diag0.root')



process.p = cms.Path(process.demo)
