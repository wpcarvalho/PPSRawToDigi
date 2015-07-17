import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# geometry
process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_10_20_2/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")

process.load("Configuration/TotemOpticsConfiguration/OpticsConfig_3500GeV_90_cfi")


process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use

    fileNames = cms.untracked.vstring(
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_5657.000.vmea_.root',
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/RPT1T2_CommonReco_run5657/run_run_5657.012.vmea_.root'
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6936.000_recoNewT1code.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6936.001_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6936.002_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6936.003_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6936.004_recoNewT1codeTry2.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6936.005_recoNewT1codeTry2.root'
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.000-1.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.001-2.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.002-3.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.003-4.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.004-5.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.005-6.root',
           'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.006-7.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.007-8.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.008-9.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.009-10.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.010-11.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.011-12.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.012-13.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.013-14.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.014-15.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.015-16.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.016-17.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.017-18.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.018-19.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.019-20.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.020-21.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.021-22.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.022-23.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.023-24.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.024-25.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.025-26.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.026-27.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.027-28.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.028-29.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.029-30.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.030-31.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.031-32.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.032-33.root',
          'rfio:/castor/cern.ch/totem/offline/Analysis/2012/T1T2RP/inelastic_reco/ver3.4_test/6936/run_6936.033-34.root'
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

process.demo.fileName = cms.string('DataReduction-r6936-newReco-AllSegs-v2-diag1.root')

process.demo.bunchText = cms.string('r6936-bunchText-diag1-v2.txt')

process.p = cms.Path(process.demo)
