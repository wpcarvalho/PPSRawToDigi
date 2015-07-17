/** \file
 * 
 * \author  F.Ferro
 */


#include <DataFormats/T1DigiSantiard/interface/T1DigiSantiard.h>

T1DigiSantiard::T1DigiSantiard (int balance, int strip, int bx) :
  balance_(balance),
  strip_(strip),
  bx_(bx),
  charge_(0)
{}
T1DigiSantiard::T1DigiSantiard (int balance, int strip, int bx,float charge) :
  balance_(balance),
  strip_(strip),
  bx_(bx),
  charge_(charge)
{}

T1DigiSantiard::T1DigiSantiard ():
  balance_(0),
  strip_(0),
  bx_(0),
  charge_(0)
{}


// Comparison
bool
T1DigiSantiard::operator == (const T1DigiSantiard& digi) const {
  if ( strip_ != digi.strip() ||
       bx_    != digi.bx() ) return false;
  return true;
}

///Precedence operator
bool 
T1DigiSantiard::operator<(const T1DigiSantiard& digi) const{

  if(digi.bx() == this->bx())
    return digi.strip()<this->strip();
  else 
    return digi.bx()<this->bx();
}


int T1DigiSantiard::strip() const { return strip_; }

int T1DigiSantiard::bx() const { return bx_; }

int T1DigiSantiard::balance() const {return balance_;}

float T1DigiSantiard::charge() const {return charge_;}

void
T1DigiSantiard::print() const {
  std::cout << "Strip " << strip() 
	    << " bx " << bx() <<std::endl;
}
