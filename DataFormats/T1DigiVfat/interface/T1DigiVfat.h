#ifndef T1DigiVfat_T1DigiVfat_h
#define T1DigiVfat_T1DigiVfat_h

/*

Digi for T1 CSCs
Author F.Ferro

*/

#include <boost/cstdint.hpp>

class T1DigiVfat{

 public:
  explicit T1DigiVfat (int strip, int threshold, int bx);
  explicit T1DigiVfat (int strip, int threshold, int bx, float charge);
  explicit T1DigiVfat (int strip, int threshold, int bx, int quality);
  explicit T1DigiVfat (int strip, int threshold, int bx, float charge, int quality);
  T1DigiVfat ();

  bool operator==(const T1DigiVfat& digi) const;
  bool operator<(const T1DigiVfat& digi) const;


  int strip() const ;
  int threshold() const;
  int bx() const;
  int quality() const;
  float charge() const;
  void print() const;

  void setThreshold(int thr){threshold_=thr;}

 private:
  int strip_;
  int threshold_;
  int32_t  bx_; 
  float charge_;
  int quality_;
};


#include<iostream>
inline std::ostream & operator<<(std::ostream & o, const T1DigiVfat& digi) 
{
  return o << " " << digi.strip()
	   << " " << digi.threshold()
	   << " " << digi.bx()
	   << " " << digi.quality();
}
#endif
