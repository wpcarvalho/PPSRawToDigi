#include "MyElectron.h"
#include <iostream>

using namespace std;

ClassImp(MyElectron)

MyElectron::MyElectron():MyPart(){
  this->Reset();
}

MyElectron::~MyElectron(){}


void MyElectron::Reset(){


 this->MyPart::Reset();
 GsfTrack.Reset();
 TrackerTrack.Reset();

 eSupClusOverP      = 0 ;
 eSeedClusOverPout  = 0 ;
 PIn		    = 0 ;
 POut		    = 0 ;


 dEtaSupClusTrVtx   = 0 ;
 dEtaSeedClusTrCalo = 0 ;
 dPhiSupClusTrVtx   = 0 ;
 dPhiSeedClusTrCalo = 0 ;
 fbrem  	    = 0 ;

 isBarrel	    = 0 ;
 isEndCap	    = 0 ;

 isEcalDriven	    = 0 ;
 isTrackerDriven    = 0 ;

 isEScaleCorr	    = 0 ;
 isMomentumCorr     = 0 ;
 nClus  	    = 0 ;
 classification     = 0 ;
 expectedInnerHits  = 0 ;

 EcalRecHitIsoDr04	  = 0 ;
 HcalDepth1TowerSumEtDr04 = 0 ;
 HcalDepth2TowerSumEtDr04 = 0 ;
 TrackIsolationDr04	  = 0 ;
 EcalRecHitIsoDr03	  = 0 ;
 HcalTowerSumEtDr03	  = 0 ;
 HcalDepth1TowerSumEtDr03 = 0 ;
 HcalDepth2TowerSumEtDr03 = 0 ;
 TrackIsolationDr03	  = 0 ;
 
 eidRobustLoose       = 0 ; 
 eidRobustTight       = 0 ;
 eidLoose	      = 0 ;
 eidTight	      = 0 ;
 
 E15			= 0 ;  
 E25Max 		= 0 ; 
 E55			= 0 ;
 CovEtaEta		= 0 ;  
 CoviEtaiEta		= 0 ;  
 HadronicOverEm 	= 0 ;  
 HcalDepth1OverEcal	= 0 ;  
 HcalDepth2OverEcal	= 0 ;
 
 dist_conv = 0 ;
 dcot_conv = 0 ;
}

void MyElectron::Print(){

  this->MyPart::Print();
  GsfTrack.Print();
  TrackerTrack.Print();

  cout << "eSupClusOverP      : " << eSupClusOverP      << endl ;
  cout << "eSeedClusOverPout  : " << eSeedClusOverPout  << endl ;
  cout << "PIn		      : " << PIn	        << endl ;
  cout << "POut		      : " << POut	        << endl ;


  cout << "dEtaSupClusTrVtx   : " << dEtaSupClusTrVtx   << endl ;
  cout << "dEtaSeedClusTrCalo : " << dEtaSeedClusTrCalo << endl ;
  cout << "dPhiSupClusTrVtx   : " << dPhiSupClusTrVtx   << endl ;
  cout << "dPhiSeedClusTrCalo : " << dPhiSeedClusTrCalo << endl ;
  cout << "fbrem  	      : " << fbrem	        << endl ;

  cout << "isBarrel	      : " << isBarrel	        << endl ;
  cout << "isEndCap	      : " << isEndCap	        << endl ;

  cout << "isEcalDriven	      : " << isEcalDriven       << endl ;
  cout << "isTrackerDriven    : " << isTrackerDriven    << endl ;

  cout << "isEScaleCorr	      : " << isEScaleCorr       << endl ;
  cout << "isMomentumCorr     : " << isMomentumCorr     << endl ;
  cout << "nClus  	      : " << nClus	        << endl ;
  cout << "classification     : " << classification     << endl ;
  cout << "expectedInnerHits  : " << expectedInnerHits  << endl ;

  cout << "EcalRecHitIsoDr04	    : " << EcalRecHitIsoDr04	    << endl ;
  cout << "HcalDepth1TowerSumEtDr04 : " << HcalDepth1TowerSumEtDr04 << endl ;
  cout << "HcalDepth2TowerSumEtDr04 : " << HcalDepth2TowerSumEtDr04 << endl ;
  cout << "TrackIsolationDr04	    : " << TrackIsolationDr04	    << endl ;
  cout << "EcalRecHitIsoDr03	    : " << EcalRecHitIsoDr03	    << endl ;
  cout << "HcalTowerSumEtDr03	    : " << HcalTowerSumEtDr03	    << endl ;
  cout << "HcalDepth1TowerSumEtDr03 : " << HcalDepth1TowerSumEtDr03 << endl ;
  cout << "HcalDepth2TowerSumEtDr03 : " << HcalDepth2TowerSumEtDr03 << endl ;
  cout << "TrackIsolationDr03	    : " << TrackIsolationDr03	    << endl ;
  
  cout << "eidRobustLoose       : " << eidRobustLoose << endl ; 
  cout << "eidRobustTight       : " << eidRobustTight << endl ;
  cout << "eidLoose	        : " << eidLoose       << endl ;
  cout << "eidTight	        : " << eidTight       << endl ;
  
  cout << "E15			: " << E15	    	  << endl ;  
  cout << "E25Max 		: " << E25Max	    	  << endl ; 
  cout << "E55			: " << E55	    	  << endl ;
  cout << "CovEtaEta		: " << CovEtaEta    	  << endl ;  
  cout << "CoviEtaiEta		: " << CoviEtaiEta  	  << endl ;  
  cout << "HadronicOverEm 	: " << HadronicOverEm 	  << endl ;  
  cout << "HcalDepth1OverEcal	: " << HcalDepth1OverEcal << endl ;  
  cout << "HcalDepth2OverEcal	: " << HcalDepth2OverEcal << endl ;
  
  cout << "dist_conv : " << dist_conv << endl ;
  cout << "dcot_conv : " << dcot_conv << endl ;

}
