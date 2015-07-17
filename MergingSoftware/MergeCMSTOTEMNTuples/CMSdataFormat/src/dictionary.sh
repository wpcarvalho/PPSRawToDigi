#!/bin/bash

if [ ! -f LinkDef.h ];then mv LinkDefh LinkDef.h;fi

rootcint -f eventdict.cc -c -p \
     MassParticles.h \
     MyBeamSpot.h \
     MyCastorDigi.h \
     MyCastorJet.h \
     MyCastorRecHit.h \
     MyDiJet.h \
     MyEvtId.h \
     MyFwdGap.h \
     MyMet.h \
     MyGenMet.h \
     MyGenKin.h \
     MyGenPart.h \
     MyGenJet.h \
     MyPUSumInfo.h \
     MyHLTrig.h \
     MyL1Trig.h \
     MyPart.h \
     MySimVertex.h \
     MyVertex.h \
     MyTracks.h \
     MyElectron.h \
     MyMuon.h \
     MyBaseJet.h \
     MyJet.h \
     MyCaloJet.h \
     MyTrackJet.h \
     MyPFJet.h \
     MyPFCand.h \
     MyCaloTower.h \
     LinkDef.h \

mv LinkDef.h LinkDefh



