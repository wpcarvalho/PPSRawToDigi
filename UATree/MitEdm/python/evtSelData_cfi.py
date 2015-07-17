# $Id: evtSelData_cfi.py,v 1.3 2010/01/05 16:39:44 edwenger Exp $

import FWCore.ParameterSet.Config as cms

evtSelData = cms.EDProducer("ProducerEvtSelData",
    hfRecHits   = cms.untracked.string('hfreco'),
    hbheRecHits = cms.untracked.string('hbhereco'),
    castorRecHits = cms.untracked.string('castorreco'),
    zdcRecHits = cms.untracked.string('zdcreco'),
    pixelRecHits = cms.untracked.string('siPixelRecHits'),
    tracks = cms.untracked.string('generalTracks')
)
