import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.TrackJetParameters_cfi import *
from RecoJets.JetProducers.AnomalousCellParameters_cfi import *

from RecoJets.JetProducers.ak5GenJets_cfi import ak5GenJets
from RecoJets.JetProducers.ak5TrackJets_cfi import ak5TrackJets

FastjetWithAreaPU = cms.PSet(
    Active_Area_Repeats = cms.int32(5),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    UE_Subtraction = cms.string('no')
)


ueAk5ChgGenJet500 = ak5GenJets.clone(
    src = cms.InputTag("chargeParticles"),
    jetPtMin       = cms.double(1.0),
    inputEtMin     = cms.double(0.5)
)

ueAk5TracksJet500 =  ak5TrackJets.clone(
    src = cms.InputTag("goodTracks"),
    jetPtMin       = cms.double(1.0),
    inputEtMin     = cms.double(0.5)
)

#ueAk5TracksJet.jetType = 'BasicJet'

UEAnalysisJetsAkOnlyMC = cms.Sequence(ueAk5ChgGenJet500)
UEAnalysisJetsAkOnlyReco = cms.Sequence(ueAk5TracksJet500)
UEAnalysisJetsAk = cms.Sequence(UEAnalysisJetsAkOnlyMC*UEAnalysisJetsAkOnlyReco)
