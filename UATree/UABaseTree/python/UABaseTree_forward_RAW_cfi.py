#   For all possibilities that are not here (like some of the input collections),
#   please look at all parameters retrieved in src/GetParameters.cc
#   All the parameters have anyway a default value !

import FWCore.ParameterSet.Config as cms

# General switches --------------------------------------------------------------------
storeZDC         = True
storeFSC         = True

# Standard Parameters For UABaseTree Process   ----------------------------------------
from uabasetree_cfi import *

uabasetree.storeL1Trig    = cms.untracked.bool(False)
uabasetree.storeL1TrigOld = cms.untracked.bool(False)
uabasetree.hlt_paths      = cms.untracked.vstring()
uabasetree.beamspots      = cms.untracked.VInputTag()

#uabasetree.tracks   = cms.untracked.VInputTag()
#uabasetree.vertices = cms.untracked.VInputTag()

#uabasetree.vpfjets   = cms.untracked.VPSet()
#uabasetree.vcalojets = cms.untracked.VPSet()

#uabasetree.basicjets = cms.untracked.VInputTag()
#uabasetree.trackjets = cms.untracked.VInputTag()

#uabasetree.mets = cms.untracked.VInputTag()

#uabasetree.castorrechits = cms.untracked.InputTag()
#uabasetree.castorjetid   = cms.untracked.InputTag()
#uabasetree.castordigis   = cms.untracked.InputTag()

#uabasetree.pfcands = cms.untracked.VInputTag()

#uabasetree.calotowercoll = cms.untracked.InputTag()

#uabasetree.muons = cms.untracked.VInputTag()

#uabasetree.electrons = cms.untracked.VInputTag()

# ZDC
if storeZDC:
    uabasetree.storeZDCInfo  = cms.untracked.bool(False)
    uabasetree.storeZDCHits  = cms.untracked.bool(False)
    uabasetree.storeZDCDigis = cms.untracked.bool(True)
    #uabasetree.zdcrechits    = cms.untracked.InputTag()
    uabasetree.zdcdigis      = cms.untracked.InputTag('hcalDigis')
else:
    uabasetree.storeZDCInfo  = cms.untracked.bool(False)
    uabasetree.storeZDCHits  = cms.untracked.bool(False)
    uabasetree.storeZDCDigis = cms.untracked.bool(False)

# FSC
if storeFSC:
    uabasetree.storeFSCInfo  = cms.untracked.bool(False)
    uabasetree.storeFSCHits  = cms.untracked.bool(False)
    uabasetree.storeFSCDigis = cms.untracked.bool(True)
    #uabasetree.fscrechits    = cms.untracked.InputTag()
    uabasetree.fscdigis      = cms.untracked.InputTag('hcalDigis')
else:
    uabasetree.storeFSCInfo  = cms.untracked.bool(False)
    uabasetree.storeFSCHits  = cms.untracked.bool(False)
    uabasetree.storeFSCDigis = cms.untracked.bool(False)
