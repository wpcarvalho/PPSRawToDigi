/** \file
 * Impl of RPCDetId
 *
 * \author Ilaria Segoni
 * \version $Id: T1DetId.cc,v 1.2 2008/09/12 16:00:33 lgrzanka Exp $
 * \date 02 Aug 2005
 */

#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <FWCore/Utilities/interface/Exception.h>


// TOTEM =7, T1 = 1

T1DetId::T1DetId():DetId(DetId::Totem,1),trind(0){}


T1DetId::T1DetId(uint32_t id):DetId(id),trind(0) {
  //  std::cout<<" constructor of the T1DetId" <<std::endl;
  if (det()!=DetId::Totem || subdetId()!=1) {
    throw cms::Exception("InvalidDetId") << "T1DetId ctor:"
					 << " det: " << det()
					 << " subdet: " << subdetId()
					 << " is not a valid T1 id";  
  }
}



T1DetId::T1DetId(unsigned int Arm, unsigned int Plane, unsigned int CSC):	      
  DetId(DetId::Totem,1),trind(0)
{
  unsigned int d=0;
  this->init(Arm,Plane,CSC,d);
}


T1DetId::T1DetId(unsigned int Arm, unsigned int Plane, unsigned int CSC, unsigned int Layer):	      
  DetId(DetId::Totem,1),trind(0)
{
  this->init(Arm,Plane,CSC,Layer);
}


void
T1DetId::init(unsigned int Arm, unsigned int Plane, unsigned int CSC,unsigned int Layer)
{
  if ( 
      (Arm != 0 && Arm !=1) ||
      Plane > 4 ||
      CSC > 5  ||
      Layer > 5
      ) {
    throw cms::Exception("InvalidDetId") << "T1DetId ctor:" 
					 << " Invalid parameters: " 
					 << " Arm "<<Arm
					 << " Plane "<<Plane
					 << " CSC "<<CSC
                                         << " Layer " << Layer
					 << std::endl;
  }

  uint32_t ok=0x72000000;
  id_ &= ok;

  id_ |= Layer << 16 |
    CSC   << 18    | 
    Plane    << 21    |
    Arm  << 24 ;
   
}

std::ostream& operator<<( std::ostream& os, const T1DetId& id ){
  os <<  " Arm "<<id.Arm()
     << " Plane "<<id.Plane()
     << " CSC "<<id.CSC()
     << " Layer "<<id.Layer();

  return os;
}
