#ifndef T1Digi_T1Digi_h
#define T1Digi_T1Digi_h

/*

Digi for T1 CSCs
Author F.Ferro

*/

#include <boost/cstdint.hpp>

class T1Digi{

 public:
  explicit T1Digi (int strip, int balance, int bx);
  T1Digi ();

  bool operator==(const T1Digi& digi) const;
  bool operator<(const T1Digi& digi) const;

  int balance() const;
  int strip() const ;
  int bx() const;
  void print() const;

 private:
  uint16_t balance_;
  uint16_t strip_;
  int32_t  bx_; 
};


#include<iostream>
inline std::ostream & operator<<(std::ostream & o, const T1Digi& digi) 
{
  return o << " " << digi.strip()
	   << " " << digi.balance()
	   << " " << digi.bx();
}
#endif
