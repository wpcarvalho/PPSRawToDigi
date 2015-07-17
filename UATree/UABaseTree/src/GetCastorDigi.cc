//-- Description: Function to retrieve castor digi information


//-- Castor Digi
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalDigi/interface/CastorDataFrame.h"
#include "DataFormats/HcalDetId/interface/HcalCastorDetId.h"

#include "CondFormats/CastorObjects/interface/CastorQIEShape.h"
#include "CondFormats/CastorObjects/interface/CastorQIECoder.h"
#include "CalibFormats/CastorObjects/interface/CastorCoderDb.h"
#include "CalibFormats/CastorObjects/interface/CastorDbService.h"
#include "CalibFormats/CastorObjects/interface/CastorDbRecord.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool CastorDigiDebug = true;

void UABaseTree::GetCastorDigi(const edm::Event& iEvent ,  const edm::EventSetup& iSetup) {
    
  castorDigis.clear();
  
  Handle<CastorDigiCollection> digicoll;
  iEvent.getByLabel(castordigis_,digicoll);
  
  // QIE coder to convert to fC
  // get conditions
  edm::ESHandle<CastorDbService> conditions;
  iSetup.get<CastorDbRecord>().get(conditions);
  const CastorQIEShape* shape = conditions->getCastorShape ();

  if (CastorDigiDebug) cout<<"number of Castor Digi: "<<digicoll->size()<<endl;
  
  //Initializing the vector
  castorDigis.assign(digicoll->size() , MyCastorDigi() );
  
  //-- loop over the digi collection (224 digis)
  for(size_t i = 0; i < digicoll->size(); i++) { 
  
    CastorDataFrame digi = (*digicoll)[i];
    HcalCastorDetId castorid = digi.id();
    
    const CastorQIECoder* channelCoder = conditions->getCastorCoder (castorid);
    CastorCoderDb coder (*channelCoder, *shape);
    CaloSamples tool;
    coder.adc2fC(digi,tool);

    castorDigis[i].mod = castorid.module();
    castorDigis[i].sec = castorid.sector();
    castorDigis[i].cha = 16*(castorid.module()-1) + castorid.sector();
    
    //-- loop over the 6 or 10 time slices for each digi
    for (int ts = 0; ts < digi.size(); ts++) {   
      castorDigis[i].adc.push_back(digi[ts].adc());
      castorDigis[i].fC.push_back(tool[ts]);
    }
        
    if (CastorDigiDebug) castorDigis[i].Print();   
  }
}
