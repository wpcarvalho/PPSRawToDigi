// -*- C++ -*-
//
// Package:     Forward
// Class  :     TotemT2OrganizationGem
//
// Implementation:
//     <Notes on implementation>
//
// Original Author: 
//         Created:  Tue May 16 10:14:34 CEST 2006
//

// system include files

// user include files
#include "SimG4CMS/Forward/interface/TotemT2OrganizationGem.h"
#include "SimG4CMS/Forward/interface/TotemNumberMerger.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4VPhysicalVolume.hh"
#include "G4VTouchable.hh" 

//
// constructors and destructor
//

TotemT2OrganizationGem :: TotemT2OrganizationGem() :
  _needUpdateUnitID(false), 
  _needUpdateData(false),
  _currentUnitID(0), 
  _currentDetectorPosition(3), 
  _currentHalfTelescope(3), 
  _currentPlane(7),
  _currentPlaneSide(3) {
  
  edm::LogInfo("ForwardSim") << "Creating TotemT2OrganizationGem";

} // TotemT2OrganizationGem

/**
 *
 */

TotemT2OrganizationGem :: ~TotemT2OrganizationGem() {
} // ~TotemT2OrganizationGem

/**
 *
 */

uint32_t TotemT2OrganizationGem :: GetUnitID(const G4Step* aStep) const {

  return const_cast<TotemT2OrganizationGem *>(this)->GetUnitID(aStep);

} // GetUnitID

/**
 *
 */

uint32_t TotemT2OrganizationGem :: GetUnitID(const G4Step* aStep) {
  
  G4VPhysicalVolume* physVol;
  const G4VTouchable* touch = aStep->GetPreStepPoint()->GetTouchable();
  physVol= touch->GetVolume(0);
  
  
  //  if(physVol->GetName() == "TotemT2gem_driftspace7r") cool = 400 + touch->GetVolume(1)->GetCopyNo()+ (touch->GetVolume(3)->GetCopyNo())*1000;
  if(physVol->GetName() == "TotemT2gem_driftspace7r") {
 
    int sDVolMotherNo = touch->GetVolume(3)->GetCopyNo();
    int sDVolNo = touch->GetVolume(1)->GetCopyNo();
    
    
    if( sDVolMotherNo  == 1 ) {

      _currentDetectorPosition = 0;
      
    } else if( sDVolMotherNo == 2 ) {

      _currentDetectorPosition = 1;
    
    } else std::cout<<"ERROR2: WrongDetectorPosition!"<<std::endl;
    
    
    if((sDVolNo & 1) == 0) {
    
      _currentPlaneSide = 1;
    
    } else _currentPlaneSide = 0;
    
    if(_currentDetectorPosition == 0) {

      if(sDVolNo <= 10) {
	
	_currentHalfTelescope = 0;
	_currentPlane = (sDVolNo - 1) / 2;
	
      }else {
	
	_currentHalfTelescope = 1;
	_currentPlane = (sDVolNo - 11) / 2;
      }
      
    } else {
      
      if(sDVolNo <= 10) {
	
	_currentHalfTelescope = 1;
	_currentPlane = (sDVolNo - 1) / 2;
	
      }else {
	
	_currentHalfTelescope = 0;
	_currentPlane = (sDVolNo - 11) / 2;
      } 
    }
    
    
  } else std::cout<<"ERROR: Physical volume is not T2 driftspace!"<<std::endl;
  
  
  
 
  
  LogDebug("ForwardSim") << _currentDetectorPosition << " " 
			 << _currentHalfTelescope << " " 
			 << _currentPlane << " " 
			 << _currentPlaneSide;
  
  uint32_t rawId = T2DetId::calculateRawId(_currentDetectorPosition,
					   _currentHalfTelescope,
					   _currentPlane,
					   _currentPlaneSide);
  
  LogDebug("ForwardSim") << " ID " << rawId;
  
  return rawId; 
  
} // GetUnitId
