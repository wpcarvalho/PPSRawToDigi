import FWCore.ParameterSet.Config as cms

RecoCTPPSFEVT = cms.PSet(
  outputCommands = cms.untracked.vstring(
    'keep TotemFEDInfos_*_*_*',
    'keep TotemTriggerCounters_*_*_*',
    'keep TotemRPDigiedmDetSetVector_*_*_*',
    'keep TotemVFATStatusedmDetSetVector_*_*_*',
    'keep TotemRPClusteredmDetSetVector_*_*_*',
    'keep TotemRPRecHitedmDetSetVector_*_*_*',
    'keep TotemRPUVPatternedmDetSetVector_*_*_*',
    'keep TotemRPLocalTrackedmDetSetVector_*_*_*'
  )
)


RecoCTPPSRECO = cms.PSet(
  outputCommands = cms.untracked.vstring(
    'keep TotemFEDInfos_*_*_*',
    'keep TotemTriggerCounters_*_*_*',
    'keep TotemRPDigiedmDetSetVector_*_*_*',
    'keep TotemVFATStatusedmDetSetVector_*_*_*',
    'keep TotemRPClusteredmDetSetVector_*_*_*',
    'keep TotemRPRecHitedmDetSetVector_*_*_*',
    'keep TotemRPUVPatternedmDetSetVector_*_*_*',
    'keep TotemRPLocalTrackedmDetSetVector_*_*_*'
  )
)


RecoCTPPSAOD = cms.PSet(
  outputCommands = cms.untracked.vstring(
    'keep TotemRPClusteredmDetSetVector_*_*_*',
    'keep TotemRPRecHitedmDetSetVector_*_*_*',
    'keep TotemRPUVPatternedmDetSetVector_*_*_*',
    'keep TotemRPLocalTrackedmDetSetVector_*_*_*'
  )
)
