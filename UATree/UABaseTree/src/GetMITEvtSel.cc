// MIT code
#include "UATree/MitEdm/interface/EvtSelData.h"

// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool MITEvtSelDebug = false;

void UABaseTree::GetMITEvtSel(const edm::Event& iEvent)
{
  using namespace mitedm;

  MITEvtSel.Reset();

  Handle<EvtSelData> EvtSel;
  iEvent.getByLabel("evtSelData",EvtSel);

  MITEvtSel.eHcalNeg = EvtSel->eHcalNeg();
  MITEvtSel.eHcalPos = EvtSel->eHcalPos();
  MITEvtSel.eHfNeg = EvtSel->eHfNeg();
  MITEvtSel.eHfPos = EvtSel->eHfPos();
  MITEvtSel.eHfNegTime = EvtSel->eHfNegTime();
  MITEvtSel.eHfPosTime = EvtSel->eHfPosTime();
  MITEvtSel.eCaNeg = EvtSel->eCastorNeg();
  MITEvtSel.eCaPos = EvtSel->eCastorPos();
  MITEvtSel.eCaNegTime = EvtSel->eCastorNegTime();
  MITEvtSel.eCaPosTime = EvtSel->eCastorPosTime();
  MITEvtSel.eZdcNeg = EvtSel->eZdcNeg();
  MITEvtSel.eZdcPos = EvtSel->eZdcPos();
  MITEvtSel.eZdcNegTime = EvtSel->eZdcNegTime();
  MITEvtSel.eZdcPosTime = EvtSel->eZdcPosTime();
  MITEvtSel.ePxbHits = EvtSel->ePxbHits();
  MITEvtSel.ePxHits = EvtSel->ePxHits();
  MITEvtSel.eClusVtxQual = EvtSel->eClusVtxQual();
  MITEvtSel.eClusVtxDiff = EvtSel->eClusVtxDiff();
  MITEvtSel.nHfNegHits = EvtSel->nHfNegHits();
  MITEvtSel.nHfPosHits = EvtSel->nHfPosHits();
  MITEvtSel.nHfTowersP = EvtSel->nHfTowersP();
  MITEvtSel.nHfTowersN = EvtSel->nHfTowersN();
  MITEvtSel.sumEsubEpPos = EvtSel->sumEsubEpPos();
  MITEvtSel.sumEaddEpPos = EvtSel->sumEaddEpPos();
  MITEvtSel.sumEsubEpNeg = EvtSel->sumEsubEpNeg();
  MITEvtSel.sumEaddEpNeg = EvtSel->sumEaddEpNeg();
  MITEvtSel.sumHfEsubEpPos = EvtSel->sumHfEsubEpPos();
  MITEvtSel.sumHfEaddEpPos = EvtSel->sumHfEaddEpPos();
  MITEvtSel.sumHfEsubEpNeg = EvtSel->sumHfEsubEpNeg();
  MITEvtSel.sumHfEaddEpNeg = EvtSel->sumHfEaddEpNeg();
  MITEvtSel.eHPTrkFrac = EvtSel->eHPTrkFrac();

  if(MITEvtSelDebug) MITEvtSel.Print();

}

