/** \file
 * 
 * \author  F.Ferro
 */


#include <DataFormats/T1DigiVfat/interface/T1DigiVfat.h>

T1DigiVfat::T1DigiVfat (int strip, int threshold, int bx) :
  strip_(strip),
  threshold_(threshold),
  bx_(bx)
  ,charge_(0), quality_(0)
{}

T1DigiVfat::T1DigiVfat (int strip, int threshold, int bx,float charge) :
  strip_(strip),
  threshold_(threshold),
  bx_(bx)
  ,charge_(charge), quality_(0)
{}

T1DigiVfat::T1DigiVfat (int strip, int threshold, int bx, int quality) :
  strip_(strip),
  threshold_(threshold),
  bx_(bx)
  ,charge_(0), quality_(quality)
{}

T1DigiVfat::T1DigiVfat (int strip, int threshold, int bx,float charge,int quality ) :
  strip_(strip),
  threshold_(threshold),
  bx_(bx)
  ,charge_(charge), quality_(quality)
{}

T1DigiVfat::T1DigiVfat ():
  strip_(0),
  bx_(0) 
  ,charge_(0), quality_(0)
{}


// Comparison
bool
T1DigiVfat::operator == (const T1DigiVfat& digi) const {
  if ( strip_ != digi.strip() ||
       bx_    != digi.bx() ) return false;
  return true;
}

///Precedence operator
bool 
T1DigiVfat::operator<(const T1DigiVfat& digi) const{

  if(digi.bx() == this->bx())
    return digi.strip()<this->strip();
  else 
    return digi.bx()<this->bx();
}


int T1DigiVfat::strip() const { return strip_; }

int T1DigiVfat::bx() const { return bx_; }

int T1DigiVfat::threshold() const {return threshold_;}

int T1DigiVfat::quality() const {return quality_;}

float T1DigiVfat::charge() const {return charge_;}

void
T1DigiVfat::print() const {
  std::cout << "Strip " << strip() 
	    << " bx " << bx() 
	    << " quality "<<quality()<< std::endl;
}
