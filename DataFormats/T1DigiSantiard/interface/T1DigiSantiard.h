#ifndef T1DigiSantiard_T1DigiSantiard_h
#define T1DigiSantiard_T1DigiSantiard_h

/*

Digi for T1 CSCs
Author F.Ferro

*/

#include <boost/cstdint.hpp>

class T1DigiSantiard{

 public:
  explicit T1DigiSantiard (int strip, int balance, int bx);
  explicit T1DigiSantiard (int strip, int balance, int bx, float charge);
  T1DigiSantiard ();

  bool operator==(const T1DigiSantiard& digi) const;
  bool operator<(const T1DigiSantiard& digi) const;

  int balance() const;
  int strip() const ;
  int bx() const;
  void print() const;
  float charge() const;

 private:
  uint16_t balance_;
  uint16_t strip_;
  int32_t  bx_; 
  float charge_;
};


#include<iostream>
inline std::ostream & operator<<(std::ostream & o, const T1DigiSantiard& digi) 
{
  return o << " " << digi.strip()
	   << " " << digi.balance()
	   << " " << digi.bx();
}
#endif
