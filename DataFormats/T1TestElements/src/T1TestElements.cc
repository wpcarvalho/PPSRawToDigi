/** \file
 * 
 * \author  F.Ferro
 */


#include <DataFormats/T1TestElements/interface/T1TestElements.h>


T1TestElements::T1TestElements (int type, int element, int bx, float charge) :
  type_(type),
  element_(element),
  bx_(bx)
									    ,charge_(charge)
{}

T1TestElements::T1TestElements ():
  element_(0),
  bx_(0) 
				 ,charge_(0)
{}


// Comparison
bool
T1TestElements::operator == (const T1TestElements& digi) const {
  if ( element_ != digi.element() ||
       bx_    != digi.bx() ) return false;
  return true;
}

///Precedence operator
bool 
T1TestElements::operator<(const T1TestElements& digi) const{

  if(digi.bx() == this->bx())
    return digi.element()<this->element();
  else 
    return digi.bx()<this->bx();
}


int T1TestElements::element() const { return element_; }

int T1TestElements::bx() const { return bx_; }

int T1TestElements::type() const {return type_;}

float T1TestElements::charge() const {return charge_;}

void
T1TestElements::print() const {
  std::cout << "Element " << element() 
	    << " bx " << bx() <<std::endl;
}
