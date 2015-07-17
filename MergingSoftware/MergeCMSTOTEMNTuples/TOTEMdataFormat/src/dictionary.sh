#!/bin/bash
if [ ! -f TotemEvent_LinkDef.h ];then mv TotemEvent_LinkDefh TotemEvent_LinkDef.h;fi
rootcint -f eventdictT.cc -c -p EventMetaData.h TriggerData.h T1Event.h T2Event.h RPRootDumpDigiInfo.h RPEvent.h RPRootDumpPattern.h RPRootDumpPatternInfo.h RPRootDumpReconstructedProton.h RPRootDumpReconstructedProtonPair.h RPRootDumpTrackInfo.h NtupleInfo.h TotemEvent_LinkDef.h 
mv TotemEvent_LinkDef.h TotemEvent_LinkDefh
