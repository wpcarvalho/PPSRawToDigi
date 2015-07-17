#include "MyMITEvtSel.h"
#include <iostream>
using namespace std;

ClassImp(MyMITEvtSel)

MyMITEvtSel::MyMITEvtSel(){
  this->Reset();
}

MyMITEvtSel::~MyMITEvtSel(){}

void MyMITEvtSel::Print(){
cout << "eHcalNeg	 : " << eHcalNeg       << endl;
cout << "eHcalPos	 : " << eHcalPos       << endl;
cout << "eHfNeg	         : " << eHfNeg         << endl;
cout << "eHfPos	         : " << eHfPos         << endl;
cout << "eHfNegTime	 : " << eHfNegTime     << endl;
cout << "eHfPosTime	 : " << eHfPosTime     << endl;
cout << "eCaNeg	         : " << eCaNeg         << endl;
cout << "eCaPos	         : " << eCaPos         << endl;
cout << "eCaNegTime	 : " << eCaNegTime     << endl;
cout << "eCaPosTime	 : " << eCaPosTime     << endl;
cout << "eZdcNeg	 : " << eZdcNeg        << endl;
cout << "eZdcPos	 : " << eZdcPos        << endl;
cout << "eZdcNegTime	 : " << eZdcNegTime    << endl;
cout << "eZdcPosTime	 : " << eZdcPosTime    << endl;
cout << "ePxbHits	 : " << ePxbHits       << endl;
cout << "ePxHits	 : " << ePxHits        << endl;
cout << "eClusVtxQual    : " << eClusVtxQual   << endl;
cout << "eClusVtxDiff    : " << eClusVtxDiff   << endl;
cout << "nHfNegHits	 : " << nHfNegHits     << endl;
cout << "nHfPosHits	 : " << nHfPosHits     << endl;
cout << "nHfTowersP	 : " << nHfTowersP     << endl;
cout << "nHfTowersN	 : " << nHfTowersN     << endl;
cout << "sumEsubEpPos    : " << sumEsubEpPos   << endl;
cout << "sumEaddEpPos    : " << sumEaddEpPos   << endl;
cout << "sumEsubEpNeg    : " << sumEsubEpNeg   << endl;
cout << "sumEaddEpNeg    : " << sumEaddEpNeg   << endl;
cout << "sumHfEsubEpPos  : " << sumHfEsubEpPos << endl;
cout << "sumHfEaddEpPos  : " << sumHfEaddEpPos << endl;
cout << "sumHfEsubEpNeg  : " << sumHfEsubEpNeg << endl;
cout << "sumHfEaddEpNeg  : " << sumHfEaddEpNeg << endl;
cout << "eHPTrkFrac	 : " << eHPTrkFrac     << endl;

}
void MyMITEvtSel::Reset(){
  eHcalNeg	  = 0;
  eHcalPos	  = 0;
  eHfNeg	  = 0;
  eHfPos	  = 0;
  eHfNegTime	  = 0;
  eHfPosTime	  = 0;
  eCaNeg	  = 0;
  eCaPos	  = 0;
  eCaNegTime	  = 0;
  eCaPosTime	  = 0;
  eZdcNeg	  = 0;
  eZdcPos	  = 0;
  eZdcNegTime	  = 0;
  eZdcPosTime	  = 0;
  ePxbHits	  = 0;
  ePxHits	  = 0;
  eClusVtxQual    = 0;
  eClusVtxDiff    = 0;
  nHfNegHits	  = 0;
  nHfPosHits	  = 0;
  nHfTowersP	  = 0;
  nHfTowersN	  = 0;
  sumEsubEpPos    = 0;
  sumEaddEpPos    = 0;
  sumEsubEpNeg    = 0;
  sumEaddEpNeg    = 0;
  sumHfEsubEpPos  = 0;
  sumHfEaddEpPos  = 0;
  sumHfEsubEpNeg  = 0;
  sumHfEaddEpNeg  = 0;
  eHPTrkFrac	  = 0;
}
