#from JetMETCorrections.Configuration.JetCorrectionServicesAllAlgos_cff import *

import FWCore.ParameterSet.Config as cms
from JetMETCorrections.Configuration.DefaultJEC_cff import *

#Declare here your ESSources that are not in JetCorrectionServices_cff ...

#########################
## L1 Offset Correction #
#########################

##-------------------- Disable the CondDB for the L1Offset (until they are included in a new global tag) -------
#ak5CaloL1Offset.useCondDB = False
#ak5PFL1Offset.useCondDB = False

##########################
## L1 FastJet Correction #
##########################

##-------------------- Disable the CondDB for the L1FastJet (until they are included in a new global tag) -------
#ak5CaloL1Fastjet.useCondDB = False
#ak5PFL1Fastjet.useCondDB = False

##-------------------- Import the Jet RECO modules -----------------------
from RecoJets.Configuration.RecoPFJets_cff import *
##-------------------- Turn-on the FastJet density calculation -----------------------
#kt6CaloJets.doRhoFastjet = True
#kt6CaloJets.Rho_EtaMax= cms.double(4.4)
kt6PFJets.doRhoFastjet = True
kt6PFJets.Rho_EtaMax= cms.double(4.4)
##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
#ak5CaloJets.doAreaFastjet = True
#ak5CaloJets.Rho_EtaMax= cms.double(4.5)
ak5PFJets.doAreaFastjet = True
ak5PFJets.Rho_EtaMax= cms.double(4.5)

ak5PFJetsL1Offset   = ak5PFJetsL2L3.clone(correctors = ['ak5PFL1Offset'])
ak5PFJetsL1Area     = ak5PFJetsL2L3.clone(correctors = ['ak5PFL1Fastjet'])


# Residual corrections not in DB ? --------------------------------------------------------------------------------
#ak5CaloResidual.useCondDB = False
#ak5PFResidual.useCondDB = False 

# tests ------------------------------------------------------------------------------------------------------------
#ak5PFJets.Active_Area_Repeats = 5
#ak5PFJets.GhostArea = 0.005 
#kt6PFJets.Active_Area_Repeats = 5
#kt6PFJets.GhostArea = 0.005


L1FastJet = cms.Sequence (kt6PFJets * ak5PFJets) 
#L1FastJet = cms.Sequence (kt6PFJets * ak5PFJets * ak5PFJetsL1Area * ak5PFJetsL1Offset)
#L1FastJet = cms.Sequence (ak5PFJets * kt6PFJets * ak5PFJetsL1Area * ak5PFJetsL1Offset)
