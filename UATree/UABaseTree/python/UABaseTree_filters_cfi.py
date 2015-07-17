#HLT Stuff
from HLTrigger.special.hltPhysicsDeclared_cfi import *
hltPhysicsDeclared.L1GtReadoutRecordTag = 'gtDigis'
#mypath=cms.Path(process.hltPhysicsDeclared+<my_other_things_here>)


#noScraping official filter
noscraping = cms.EDFilter("FilterOutScraping",
    applyfilter = cms.untracked.bool(True),
    debugOn = cms.untracked.bool(False), ## Or 'True' to get some per-event info
    numtrack = cms.untracked.uint32(10),
    thresh = cms.untracked.double(0.25)
)

#evtSelData production + MiTFilter 
from UATree.UABaseTree.UABaseTree_tracking_cfi import * 
from UATree.MitEdm.evtSelData_cfi import *
from UATree.MitEdm.FilterEvtSel_cff import *
mitfilter = cms.Sequence(redoSiHits * evtSelData * looseEvtSelFilter)

