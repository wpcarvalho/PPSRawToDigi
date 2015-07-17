// Description: Function to retrieve CaloTowers

#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"

#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

Bool_t CaloTowerDebug = false;

void UABaseTree::GetCaloTower(const edm::Event& iEvent){

   caloTowers.clear();

   Handle<CaloTowerCollection> CaloTowers;
   try {
     iEvent.getByLabel(calotowercoll_,CaloTowers);
   }
   catch ( ... ) {
     cout << "[UABaseTree::GetCaloTower] Can't find the collection " << calotowercoll_ << endl;
   }
   
   caloTowers.assign( CaloTowers->size() , MyCaloTower() );

   Int_t i = 0;
   for (CaloTowerCollection::const_iterator iCT = CaloTowers->begin() ; iCT != CaloTowers->end() ; ++iCT , ++i) {
     
     caloTowers[i].SetPxPyPzE(iCT->px(), iCT->py(), iCT->pz(), iCT->energy());
     caloTowers[i].emEnergy  = iCT->emEnergy();
     caloTowers[i].hadEnergy = iCT->hadEnergy();

     //-- loop over CaloTower constituents
     for(size_t iconst = 0; iconst < iCT->constituentsSize(); iconst++){
     
       DetId detId = iCT->constituent(iconst);

       if(detId.det()==DetId::Ecal){
         EcalSubdetector ecalSubDet = (EcalSubdetector)detId.subdetId();
	 if(ecalSubDet == EcalBarrel) caloTowers[i].hasEB = true;
         else if(ecalSubDet == EcalEndcap) caloTowers[i].hasEE = true;
       }

       else if(detId.det()==DetId::Hcal){
	 HcalDetId hcalDetId(detId);
	 if(hcalDetId.subdet()==HcalBarrel) caloTowers[i].hasHB = true;
	 else if(hcalDetId.subdet()==HcalEndcap) caloTowers[i].hasHE = true;
	 else if(hcalDetId.subdet()==HcalForward) caloTowers[i].hasHF = true;
       } 

     }

     caloTowers[i].zside = iCT->zside();

     if(CaloTowerDebug) caloTowers[i].Print();

   } // end for CaloTowerCollection 

}
