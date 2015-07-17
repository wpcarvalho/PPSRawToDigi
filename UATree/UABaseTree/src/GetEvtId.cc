// Description: Function to retrieve Evt Id information

// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool EvtIdDebug = false;

void UABaseTree::GetEvtId(const edm::Event& iEvent)
{
  evtId.Reset();

  evtId.Run       = iEvent.id().run() ;   
  evtId.Evt       = iEvent.id().event() ;
  evtId.LumiSect  = iEvent.luminosityBlock();
  edm::Timestamp Time=iEvent.time();
  evtId.Time      = Time.value(); 
  evtId.IsData    = iEvent.isRealData();
  evtId.ExpType   = iEvent.experimentType();
  evtId.Bunch     = iEvent.bunchCrossing();
  evtId.Orbit     = iEvent.orbitNumber();

  if(EvtIdDebug) evtId.Print();

}

