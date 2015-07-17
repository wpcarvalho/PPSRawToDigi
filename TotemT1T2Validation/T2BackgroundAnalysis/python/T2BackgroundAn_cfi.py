import FWCore.ParameterSet.Config as cms

T2BackgroundAn = cms.EDAnalyzer("T2BackgroundAn",

  verbosity = cms.int32(0),
  outputFileName = cms.string("file:OutputFile.root"),
  VtxClusterDistance = cms.double(5.),
  fastSimulation= cms.bool(True),                              
  MaxPadCluOfFittingCurves= cms.int32(65),
  NameOfGenerator = cms.string("Py8"),#Chose between EPOS and Py8
  UseselectedplanesforAPM= cms.bool(False),
  CluLabel = cms.string("CluLabel"),
  HitLabel = cms.string("HitLabel"),
  RoadLabel = cms.string("RoadLabel"),
  RoadInstanceLabel= cms.string("RoadInstanceLabel"),
  TrackLabel= cms.string("T2TrackColl2"),
  PadRoadFinderAnalysis= cms.bool(True),
  T2_QuarterUsed=  cms.vint32(0),                           
  selected_event= cms.int32(1),
  fastAnalysis= cms.bool(False),
  ZEffiCutImpact= cms.vdouble(5000.0,5000.0,5000.0,5000.0),
  EnergyCutinPrimaryEfficiency= cms.double(0.),
  PtCutinPrimaryEfficiency= cms.double(0.03),
  numhitRequiredFormMatching=cms.int32(1),
  vtxseednumber=cms.int32(1)
)





