#ifndef T1TestElements_T1TestElements_h
#define T1TestElements_T1TestElements_h

/*

Digi for T1 CSCs
Author F.Ferro

*/

#include <boost/cstdint.hpp>

class T1TestElements{

 public:
  explicit T1TestElements (int strip, int type, int bx, float charge);
  T1TestElements ();

  bool operator==(const T1TestElements& digi) const;
  bool operator<(const T1TestElements& digi) const;

  int type() const;
  int element() const ;
  int bx() const;
  void print() const;
  float charge() const;

 private:
  int type_;
  int element_;
  int  bx_; 
  float charge_;
};


#include<iostream>
inline std::ostream & operator<<(std::ostream & o, const T1TestElements& digi) 
{
  return o << " " << digi.element()
	   << " " << digi.type()
	   << " " << digi.bx();
}
#endif
