// Description: Function to retrieve FwdGap info (given by Fwd Group)

// user include files
#include "DataFormats/CaloTowers/interface/CaloTower.h" 
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"
#include "DataFormats/CaloTowers/interface/CaloTowerDetId.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalDetId/interface/HcalSubdetector.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"


// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"

bool FwdGapDebug = false;

void UABaseTree::GetFwdGap(const edm::Event& iEvent)
{

  int nTowersHF_plus = 0;
  int nTowersHF_minus = 0;
  int nTowersHE_plus = 0;
  int nTowersHE_minus = 0;
  int nTowersHB_plus = 0;
  int nTowersHB_minus = 0;
  int nTowersEE_plus = 0;
  int nTowersEE_minus = 0;
  int nTowersEB_plus = 0;
  int nTowersEB_minus = 0; 
  //Sum(E)
  double sumEHF_plus = 0.;
  double sumEHF_minus = 0.;
  double sumEHE_plus = 0.;
  double sumEHE_minus = 0.;
  double sumEHB_plus = 0.;
  double sumEHB_minus = 0.;
  double sumEEE_plus = 0.;
  double sumEEE_minus = 0.;
  double sumEEB_plus = 0.;
  double sumEEB_minus = 0.;
  // Sum(ET)
  double sumETHF_plus = 0.;
  double sumETHF_minus = 0.;
  double sumETHE_plus = 0.;
  double sumETHE_minus = 0.;
  double sumETHB_plus = 0.;
  double sumETHB_minus = 0.;
  double sumETEE_plus = 0.;
  double sumETEE_minus = 0.;
  double sumETEB_plus = 0.;
  double sumETEB_minus = 0.;

  fwdGap.Reset();

  // Calo tower collection from event
  edm::Handle<CaloTowerCollection> towerCollectionH;
  iEvent.getByLabel(calotower_,towerCollectionH);
  const CaloTowerCollection& towerCollection = *towerCollectionH;

  // Loop over calo towers
  CaloTowerCollection::const_iterator calotower = towerCollection.begin();
  CaloTowerCollection::const_iterator calotowers_end = towerCollection.end();
  for(; calotower != calotowers_end; ++calotower) {
     //bool hasHCAL = false;
     bool hasHF = false;
     bool hasHE = false;
     bool hasHB = false;
     //bool hasHO = false;
     //bool hasECAL = false;
     bool hasEE = false;
     bool hasEB = false;     
     for(size_t iconst = 0; iconst < calotower->constituentsSize(); iconst++){
        DetId detId = calotower->constituent(iconst);
        if(detId.det()==DetId::Hcal){
           //hasHCAL = true;
           HcalDetId hcalDetId(detId);
           if(hcalDetId.subdet()==HcalForward) hasHF = true;
           else if(hcalDetId.subdet()==HcalEndcap) hasHE = true;
           else if(hcalDetId.subdet()==HcalBarrel) hasHB = true;
           //else if(hcalDetId.subdet()==HcalOuter) hasHO = true;  
        } else if(detId.det()==DetId::Ecal){
           //hasECAL = true;
           EcalSubdetector ecalSubDet = (EcalSubdetector)detId.subdetId();
           if(ecalSubDet == EcalEndcap) hasEE = true;
           else if(ecalSubDet == EcalBarrel) hasEB = true;
        }
     }

     int zside = calotower->zside();
     double caloTowerEnergy = calotower->energy();
     // FIXME
     //double caloTowerET = calotower->et(primVtx.position());
     //double caloTowerET = calotower->et(primVtx.z());
     double caloTowerET = calotower->et();

     // HCAL: Towers made of at least one component from HB,HE,HF
     if( hasHF && !hasHE ){
        if( caloTowerEnergy >= energyThresholdHF_ ){
           if(zside >= 0){
              ++nTowersHF_plus;
              sumEHF_plus += caloTowerEnergy; 
              sumETHF_plus += caloTowerET;
           } else{
              ++nTowersHF_minus;
              sumEHF_minus += caloTowerEnergy;
              sumETHF_minus += caloTowerET;
           } 
        }
     } else if( hasHE && !hasHF && !hasHB ){
        if( caloTowerEnergy >= energyThresholdHE_ ){
           if(zside >= 0){
              ++nTowersHE_plus;
              sumEHE_plus += caloTowerEnergy;
              sumETHE_plus += caloTowerET;
           } else{
              ++nTowersHE_minus;
              sumEHE_minus += caloTowerEnergy;
              sumETHE_minus += caloTowerET;
           }
        }
     } else if( hasHB && !hasHE ){
        if( caloTowerEnergy >= energyThresholdHB_ ){
           if(zside >= 0){
              ++nTowersHB_plus;
              sumEHB_plus += caloTowerEnergy;
              sumETHB_plus += caloTowerET;
           } else{
              ++nTowersHB_minus;
              sumEHB_minus += caloTowerEnergy;
              sumETHB_minus += caloTowerET;
           }
        }
     }

     // ECAL: Towers made of at least one component from EB,EE
     if( hasEE && !hasEB ){
        if( caloTowerEnergy >= energyThresholdEE_ ){
           if(zside >= 0){
              ++nTowersEE_plus;
              sumEEE_plus += caloTowerEnergy;
              sumETEE_plus += caloTowerET;
           } else{
              ++nTowersEE_minus;
              sumEEE_minus += caloTowerEnergy;
              sumETEE_minus += caloTowerET;
           }
        }
     } else if( hasEB && !hasEE ){
        if( caloTowerEnergy >= energyThresholdEB_ ){
           if(zside >= 0){
              ++nTowersEB_plus;
              sumEEB_plus += caloTowerEnergy;
              sumETEB_plus += caloTowerET;
           } else{
              ++nTowersEB_minus;
              sumEEB_minus += caloTowerEnergy;
              sumETEB_minus += caloTowerET;
           }
        }
     }
  }

   fwdGap.nTowersHF_plus = nTowersHF_plus ;
   fwdGap.nTowersHF_minus = nTowersHF_minus ;
   fwdGap.nTowersHE_plus = nTowersHE_plus ;
   fwdGap.nTowersHE_minus = nTowersHE_minus ;
   fwdGap.nTowersHB_plus = nTowersHB_plus ;
   fwdGap.nTowersHB_minus = nTowersHB_minus ;
   fwdGap.nTowersEE_plus = nTowersEE_plus ;
   fwdGap.nTowersEE_minus = nTowersEE_minus ;
   fwdGap.nTowersEB_plus = nTowersEB_plus ;
   fwdGap.nTowersEB_minus = nTowersEB_minus ;
   fwdGap.sumEHF_plus = sumEHF_plus ;
   fwdGap.sumEHF_minus = sumEHF_minus ;
   fwdGap.sumEHE_plus = sumEHE_plus ;
   fwdGap.sumEHE_minus = sumEHE_minus ;
   fwdGap.sumEHB_plus = sumEHB_plus ;
   fwdGap.sumEHB_minus = sumEHB_minus ;
   fwdGap.sumEEE_plus = sumEEE_plus ;
   fwdGap.sumEEE_minus = sumEEE_minus ;
   fwdGap.sumEEB_plus = sumEEB_plus ;
   fwdGap.sumEEB_minus = sumEEB_minus ;
   fwdGap.sumETHF_plus = sumETHF_plus ;
   fwdGap.sumETHF_minus = sumETHF_minus ;
   fwdGap.sumETHE_plus = sumETHE_plus ;
   fwdGap.sumETHE_minus = sumETHE_minus ;
   fwdGap.sumETHB_plus = sumETHB_plus ;
   fwdGap.sumETHB_minus = sumETHB_minus ;
   fwdGap.sumETEE_plus = sumETEE_plus ;
   fwdGap.sumETEE_minus = sumETEE_minus ;
   fwdGap.sumETEB_plus = sumETEB_plus ;
   fwdGap.sumETEB_minus = sumETEB_minus ;

  if(FwdGapDebug) fwdGap.Print();
}

