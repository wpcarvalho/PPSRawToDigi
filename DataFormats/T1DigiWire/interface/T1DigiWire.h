#ifndef T1DigiWire_T1DigiWire_h
#define T1DigiWire_T1DigiWire_h

/*

Digi for T1 CSCs
Author F.Ferro

*/

#include <boost/cstdint.hpp>

class T1DigiWire{

 public:
  explicit T1DigiWire (int wire, int bx, float charge);
  explicit T1DigiWire (int wire, int bx);
  explicit T1DigiWire (int wire, int bx, float charge, int quality);
  explicit T1DigiWire (int wire, int bx, int quality);
  T1DigiWire ();

  bool operator==(const T1DigiWire& digi) const;
  bool operator<(const T1DigiWire& digi) const;

  int wire() const;
  int bx() const;
  void print() const;

  float charge() const;
  int quality() const;

 private:
  uint16_t wire_;
  int32_t  bx_; 

  float charge_;
  int quality_;

};


#include<iostream>
inline std::ostream & operator<<(std::ostream & o, const T1DigiWire& digi) 
{
  return o << " " << digi.wire()
	   << " " << digi.bx()
	   << " " << digi.quality();
  
}
#endif
