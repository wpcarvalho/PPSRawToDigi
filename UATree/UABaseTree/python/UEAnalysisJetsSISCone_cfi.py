import FWCore.ParameterSet.Config as cms

from RecoJets.JetProducers.FastjetParameters_cfi import *
from RecoJets.JetProducers.sc5GenJets_cfi import sisCone5GenJets
from RecoJets.JetProducers.sc5TrackJets_cfi import sisCone5TrackJets



FastjetWithAreaPU = cms.PSet(
    Active_Area_Repeats = cms.int32(5),
    GhostArea = cms.double(0.01),
    Ghost_EtaMax = cms.double(6.0),
    UE_Subtraction = cms.string('no')
)

#MC jet

ueSisCone5ChgGenJet500 = sisCone5GenJets.clone( 
src = cms.InputTag("chargeParticles"),
jetPtMin       = cms.double(1.0),
inputEtMin     = cms.double(0.5)
 )


#RECO jet Tracks

ueSisCone5TracksJet500 = sisCone5TrackJets.clone(
    src = cms.InputTag("goodTracks"),
jetPtMin       = cms.double(1.0),
inputEtMin     = cms.double(0.5)
)

UEAnalysisNeededJets = cms.Sequence( ueSisCone5ChgGenJet500 * ueSisCone5TracksJet500 )
