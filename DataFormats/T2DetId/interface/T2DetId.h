#ifndef T2DetId_h
#define T2DetId_h

/**
 * Wild copy of T1DetId.h plus some addons and modifications
 *
 * Author of this piece of plagiarism: Erik Brücken
 * email: brucken@cc.helsinki.fi
 * Updated: 2008-05-27
 */

#include <DataFormats/DetId/interface/DetId.h>
#include <FWCore/Utilities/interface/Exception.h>

#include <iosfwd>
#include <iostream>

class T2DetId :public DetId {
  
 public:
      
  T2DetId();

  // Construct from a packed id. It is required that the Detector part of
  // id is Totem and the SubDet part is T2, otherwise an exception is thrown.
  explicit T2DetId(uint32_t id);


  // Construct from fully qualified identifier.
  T2DetId(unsigned int arm,            // (Left - Right)
	  unsigned int halfTelescope,  // (Inner - Outer)
	  unsigned int plane,          // (1 - 5)
	  unsigned int planeSide);     // (Front - Back)
    
  // Bit 24 = arm: 0==z>0 1==z<0
  // Bit 23 = halfTelescope: 0==Inner 1==Outer
  // Bits [20:22] = plane: [0:4]
  // Bit 19 = planeSide: 0==Front 1==Back

  int arm() const{

    return int((id_>>startArmBit) & 0X1);
  }

  static int arm(uint32_t detId) {

    return int((detId>>startArmBit) & 0X1);
  }

  int halfTelescope() const{
   
    return int((id_>>startHalfTelescopeBit) & 0X1);
  }

  static int halfTelescope(uint32_t detId) {
   
    return int((detId>>startHalfTelescopeBit) & 0X1);
  }

  int plane() const{

    return int((id_>>startPlaneBit) & 0X7);
  }
  
  static int plane(uint32_t detId) {

    return int((detId>>startPlaneBit) & 0X7);
  }
  
  int planeSide() const{
  
    return int((id_>>startPlaneSideBit) & 0X1);
  }

  static int planeSide(uint32_t detId) {
  
    return int((detId>>startPlaneSideBit) & 0X1);
  }
 
  int trIndex() const{
    
    return trind;
  }
  
  void set(unsigned int arm, unsigned int halfTelescope, 
	   unsigned int plane, unsigned int planeSide) {
    
    this->init(arm, halfTelescope, plane, planeSide);
  }
  
  static const int startArmBit = 24;
  static const int startHalfTelescopeBit = 23;  
  static const int startPlaneBit = 20;
  static const int startPlaneSideBit = 19;
  
  static uint32_t calculateRawId(unsigned int arm, unsigned int halfTelescope, 
				 unsigned int plane, unsigned int planeSide);


 private:
  
  void init(unsigned int arm, unsigned int halfTelescope, 
	    unsigned int plane, unsigned planeSide); 
  
  int trind;
  
}; // T2DetId

std::ostream& operator<<( std::ostream& os, const T2DetId& id );

#endif // T2DetId_h
