/** \file
 * 
 * \author  F.Ferro
 */


#include <DataFormats/T1Digi/interface/T1Digi.h>

T1Digi::T1Digi (int balance, int strip, int bx) :
  balance_(balance),
  strip_(strip),
  bx_(bx)
{}

T1Digi::T1Digi ():
  strip_(0),
  bx_(0) 
{}


// Comparison
bool
T1Digi::operator == (const T1Digi& digi) const {
  if ( strip_ != digi.strip() ||
       bx_    != digi.bx() ) return false;
  return true;
}

///Precedence operator
bool 
T1Digi::operator<(const T1Digi& digi) const{

  if(digi.bx() == this->bx())
    return digi.strip()<this->strip();
  else 
    return digi.bx()<this->bx();
}


int T1Digi::strip() const { return strip_; }

int T1Digi::bx() const { return bx_; }

int T1Digi::balance() const {return balance_;}

void
T1Digi::print() const {
  std::cout << "Strip " << strip() 
	    << " bal " << balance() 
	    << " bx " << bx() <<std::endl;
}
