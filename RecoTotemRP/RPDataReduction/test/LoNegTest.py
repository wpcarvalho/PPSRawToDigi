import FWCore.ParameterSet.Config as cms

process = cms.Process("Demo")

process.load("FWCore.MessageService.MessageLogger_cfi")

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )


# geometry
process.load("Configuration.TotemCommon.geometryRP_real_cfi")
process.XMLIdealGeometryESSource.geomXMLFiles.append("Geometry/TotemRPData/data/2011_10_20/RP_Dist_Beam_Cent.xml")
process.TotemRPGeometryESModule = cms.ESProducer("TotemRPGeometryESModule")



process.source = cms.Source("PoolSource",
    # replace 'myfile.root' with the source file you want to use

    fileNames = cms.untracked.vstring(
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_5657.000.vmea_.root',
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/RPT1T2_CommonReco_run5657/run_run_5657.012.vmea_.root'
#          'rfio:/castor/cern.ch/totem/offline/Reco/2011/T2/T1T2_5966_5967_5KeV/run_run_5967.005.vmea_.root'
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6814.000_recoRPT1T2-wrongAlign-June29.root'
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6884.000_recoRPT2-noAlign-correctGeom.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6884.001_recoRPT2-noAlign-correctGeom.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6884.002_recoRPT2-noAlign-correctGeom.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6884.007_recoRPT2-noAlign-correctGeom.root'
           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6883.011_recoRPT2-noAlign-correctGeom.root',
           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6883.021_recoRPT2-noAlign-correctGeom.root',
           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6883.031_recoRPT2-noAlign-correctGeom.root',
           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6883.041_recoRPT2-noAlign-correctGeom.root',
           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6883.051_recoRPT2-noAlign-correctGeom.root',
           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6883.061_recoRPT2-noAlign-correctGeom.root'

#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6917.007_recoRPT2-noAlign-correctGeom.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6917.017_recoRPT2-noAlign-correctGeom.root',
#           'rfio:/castor/cern.ch/totem/offline/Reco/2011/Physics/run_6917.027_recoRPT2-noAlign-correctGeom.root'

    )
)

process.load("RecoTotemRP.RPDataReduction.RPDataReduction_cfi")

process.demo.verb = cms.uint32(4)

process.demo.readMultiTrk = cms.bool(True)
process.demo.readT2 = cms.bool(True)
process.demo.readLoNeg = cms.bool(True)
process.demo.readClusters = cms.bool(True)


process.demo.tmcModule=cms.string('Raw2DigiProducer')
process.demo.tmcProd=cms.string('rpCCOutput')
process.demo.whichDiag=cms.uint32(0)

#process.demo.fileName = cms.string('DataReduction-r5967-005-officialReco-diag0-45tp_56bt-LoNeg-p3.root')
process.demo.fileName = cms.string('DataReduction-r6883-6segs-0x1-RP_OR-diag0-45tp_56bt-LoNegSDCombi-v3-check10of12rpHitNum.root')

#= cms.EDAnalyzer('RPDataReduction',
#    nPl=cms.uint32(3),
#    nMaxPrPl=cms.uint32(10),
#    verb=cms.uint32(4),
#    sigmaCut=cms.double(3.0),andNotOr=cms.uint32(0),
#    elCut56=cms.double(8.4), elCut45t=cms.double(8.5), elCut45b=cms.double(0.0),
#    errX=cms.double(2.0),xLow=cms.untracked.double(-0.6),xHigh=cms.untracked.double(0.0),
#    errY=cms.double(0.3), whichDiag=cms.uint32(0),
   #arm 5-6, far pots are 1mm further from the beam than the near pots
#    refPot=cms.uint32(24),
#    oppositePot=cms.uint32(125),
#    checkPot=cms.uint32(20),
#    fileName = cms.string('r3581-024vs020-UorV-Tracking-ElNoClu-6a-8p5mmCut-newReco-diag4pot-OffIfMultiTS-allEl-0mmCut.root')

#)


process.p = cms.Path(process.demo)
