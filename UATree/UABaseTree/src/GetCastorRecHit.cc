//-- Description: Function to retrieve castor rechit information

//-- ROOT include
#include <TRandom3.h>

//-- Castor RecHit
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/CastorRecHit.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool CastorRecHitDebug = false;

void UABaseTree::GetCastorRecHit(const edm::Event& iEvent) {
  
  Double_t GeVtofC = 62.5;    //-- 1 / 0.016

  Double_t mean_cal = 0;      //-- for the variation of the intercalibration
  Double_t sigma_cal = 0.40;  //-- http://indico.cern.ch/getFile.py/access?contribId=4&resId=1&materialId=slides&confId=76797
  
  castorRecHits.clear();
  MyCastorRecHit mycastorrechit;
  
  Handle<CastorRecHitCollection> rechitcoll;
  iEvent.getByLabel(castorrechits_,rechitcoll);
  
  if (CastorRecHitDebug) cout<<"number of Castor Digi: "<<rechitcoll->size()<<endl;

  //-- random number needs to be initialized for EACH event
  //-- in order to get the SAME random numbers sequence for EACH event
  //-- because we want the SAME calibration for EACH event (do not change calibration event by event)
  
  TRandom3* rnd = new TRandom3();

  //-- loop over the rechit collection (224 rechits)
  if(rechitcoll.isValid()) {
    
    for(size_t i = 0; i < rechitcoll->size(); ++i) {

      const CastorRecHit & rh = (*rechitcoll)[i];
      HcalCastorDetId castorid = rh.id();

      mycastorrechit.mod = castorid.module();
      mycastorrechit.sec = castorid.sector();
      mycastorrechit.cha = 16*(castorid.module()-1) + castorid.sector();
  
      mycastorrechit.energy = rh.energy();
      mycastorrechit.fC = rh.energy()*GeVtofC;
      mycastorrechit.time = rh.time();

      Double_t smearing = rnd->Gaus(mean_cal,sigma_cal);
      mycastorrechit.smearing = smearing;
      mycastorrechit.energy_smeared = rh.energy()*(1+smearing);
      mycastorrechit.fC_smeared = rh.energy()*GeVtofC*(1+smearing);
      
      castorRecHits.push_back(mycastorrechit);
    
      if (CastorRecHitDebug) mycastorrechit.Print();   
    }
  }
}
