#ifndef T1DetId_h
#define T1DetId_h


#include <DataFormats/DetId/interface/DetId.h>
#include <FWCore/Utilities/interface/Exception.h>

#include <iosfwd>
#include <iostream>

class T1DetId :public DetId {
  
 public:
      
  T1DetId();

  /// Construct from a packed id. It is required that the Detector part of
  /// id is Totem and the SubDet part is T1, otherwise an exception is thrown.
  explicit T1DetId(uint32_t id);


  /// Construct from fully qualified identifier.
  T1DetId(unsigned int Arm, 
	  unsigned int Plane,
	  unsigned int CSC);
    
  T1DetId(unsigned int Arm, 
	  unsigned int Plane,
	  unsigned int CSC,
	  unsigned int Layer);

  /// Bit 24 = Arm: 0=z>0 1=z<0
  /// Bits [21:23] Plane
  /// Bits [18:20] CSC
  /// Bits [16:17] Layer

  int Arm() const{
    return int((id_>>startArmBit) & 0X1);
  }

  int Plane() const{
    return int((id_>>startPlaneBit) & 0X7);
  }

  int CSC() const{
    return int((id_>>startCSCBit) & 0X7);
  }

  int Layer() const{
    return int((id_>>startLayerBit) & 0X3);
  }

  int trIndex() const{
    return trind;
  }

  void set(unsigned int a, unsigned int b, unsigned int c){
    unsigned int d=0;
    this->init(a,b,c,d);
  }

  void setLayer(unsigned int a, unsigned int b, unsigned int c, unsigned int d){
    this->init(a,b,c,d);
  }

  static const int startArmBit = 24;
  static const int startPlaneBit = 21;
  static const int startCSCBit = 18;
  static const int startLayerBit = 16;
 
 private:
  void init(unsigned int Arm, unsigned int Plane, unsigned int CSC, unsigned intLayer); 
  int trind;

}; // T1DetId

std::ostream& operator<<( std::ostream& os, const T1DetId& id );

#endif
