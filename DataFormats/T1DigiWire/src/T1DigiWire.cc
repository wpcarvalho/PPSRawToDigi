/** \file
 * 
 * \author  F.Ferro
 */


#include <DataFormats/T1DigiWire/interface/T1DigiWire.h>

T1DigiWire::T1DigiWire (int wire, int bx, float charge) :
  wire_(wire),
  bx_(bx)
  ,charge_(charge),
 quality_(0)

{}

T1DigiWire::T1DigiWire (int wire, int bx) :
  wire_(wire),
  bx_(bx)
  ,charge_(0),
  quality_(0)

{}

T1DigiWire::T1DigiWire (int wire, int bx, float charge, int quality) :
  wire_(wire),
  bx_(bx),
  charge_(charge),
  quality_(quality)
{}

T1DigiWire::T1DigiWire (int wire, int bx, int quality) :
  wire_(wire),
  bx_(bx),
  charge_(0),
  quality_(quality)


{}

T1DigiWire::T1DigiWire ():
  wire_(0),
  bx_(0),
  charge_(0),
  quality_(0)
{}


// Comparison
bool
T1DigiWire::operator == (const T1DigiWire& digi) const {
  if ( wire_ != digi.wire() ||
       bx_    != digi.bx() ) return false;
  return true;
}

///Precedence operator
bool 
T1DigiWire::operator<(const T1DigiWire& digi) const{

  if(digi.bx() == this->bx())
    return digi.wire()<this->wire();
  else 
    return digi.bx()<this->bx();
}


int T1DigiWire::wire() const { return wire_; }

int T1DigiWire::bx() const { return bx_; }

float T1DigiWire::charge() const { return charge_; }

int T1DigiWire::quality() const { return quality_;}

void
T1DigiWire::print() const {
  std::cout << "Wire " << wire() 
	    << " bx " << bx() 
	    << " quality " <<quality() << std::endl;
}
