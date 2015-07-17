import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# geometry
process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_10_20_1/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

process.load("Configuration/TotemOpticsConfiguration/OpticsConfig_3500GeV_90_cfi")


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use

    fileNames = cms.untracked.vstring(
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6895.000_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6895.001_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6895.012_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6895.024_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6895.025_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6895.036_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6895.037_recoNewT1codeTry2.root'

#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_5657.000.vmea_.root',
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/RPT1T2_CommonReco_run5657/run_run_5657.012.vmea_.root'
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.000-1.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.001-2.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.002-3.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.003-4.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.004-5.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.005-6.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.006-7.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.007-8.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.008-9.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.009-10.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.010-11.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.011-12.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.012-13.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.013-14.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.014-15.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.015-16.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.016-17.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.017-18.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.018-19.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.019-20.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.020-21.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.021-22.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.022-23.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.023-24.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.024-25.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.025-26.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.026-27.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.027-28.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.028-29.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.029-30.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.030-31.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.031-32.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.032-33.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.033-34.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.034-35.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.035-36.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6895/run_6895.036-37.root'

    )
)

process.load("RecoTotemRP.RPDataReduction.RPDataReduction_cfi")

process.demo.verb = cms.uint32(4)

process.demo.readMultiTrk = cms.bool(True)
process.demo.readT2 = cms.bool(True)
process.demo.readT1 = cms.bool(True)
process.demo.readLoNeg = cms.bool(True)
process.demo.readClusters = cms.bool(True)
process.demo.readRecoProt = cms.bool(True)


process.demo.tmcModule=cms.string('Raw2DigiProducer')
process.demo.tmcProd=cms.string('rpCCOutput')
process.demo.whichDiag=cms.uint32(1)

process.demo.fileName = cms.string('DataReduction-newReco-r6895-AllSegs-fiver4-diag1.root')

process.demo.bunchText = cms.string('r6895-bunchText-diag1.txt')

process.p = cms.Path(process.demo)
