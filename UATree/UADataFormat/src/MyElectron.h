#ifndef __MyElectron_H__
#define __MyElectron_H__

#include <vector>
#include "MyPart.h"
#include "MyTracks.h"

class MyElectron : public MyPart {

  public :
    MyElectron();
    ~MyElectron();
    
    void Reset();
    void Print();
    
    MyTracks    GsfTrack;    
    MyTracks    TrackerTrack;

//    MySupClus   SupClus;
    
    // Basic properties            
    Double_t    eSupClusOverP      ;
    Double_t    eSeedClusOverPout  ; 
    Double_t    PIn                ;
    Double_t    POut               ;
    
    
    Double_t    dEtaSupClusTrVtx   ; 
    Double_t    dEtaSeedClusTrCalo ; 
    Double_t    dPhiSupClusTrVtx   ; 
    Double_t    dPhiSeedClusTrCalo ; 
    Double_t    fbrem              ; 
    
    Bool_t      isBarrel           ; 
    Bool_t      isEndCap           ;
    
    Bool_t      isEcalDriven       ;
    Bool_t      isTrackerDriven    ;
   
    Bool_t      isEScaleCorr       ; 
    Bool_t      isMomentumCorr     ; 
    Int_t       nClus              ; 
    Int_t       classification     ; 
    Int_t       expectedInnerHits  ;
    

    // Extra Isolation variables 
    
   
    Double_t    EcalRecHitIsoDr04        ;
    Double_t    HcalDepth1TowerSumEtDr04 ;
    Double_t    HcalDepth2TowerSumEtDr04 ;
    Double_t    TrackIsolationDr04       ;
    Double_t    EcalRecHitIsoDr03        ;
    Double_t    HcalTowerSumEtDr03       ;
    Double_t    HcalDepth1TowerSumEtDr03 ;
    Double_t    HcalDepth2TowerSumEtDr03 ;
    Double_t    TrackIsolationDr03       ;
    
   
    // Id variables
    
    Bool_t    eidRobustLoose	   ; 
    Bool_t    eidRobustTight	   ;
    Bool_t    eidLoose  	   ;
    Bool_t    eidTight  	   ;
    
    // Extra shape variables
    
    Double_t E15                    ;                 
    Double_t E25Max                 ; 
    Double_t E55                    ;
    Double_t CovEtaEta              ;      
    Double_t CoviEtaiEta            ;     
    Double_t HadronicOverEm         ;  
    Double_t HcalDepth1OverEcal     ;    
    Double_t HcalDepth2OverEcal     ;
    
    // conversion variables
    
    Double_t dist_conv ;
    Double_t dcot_conv ;
    
   
  private:



  ClassDef (MyElectron,1)
};

#endif

