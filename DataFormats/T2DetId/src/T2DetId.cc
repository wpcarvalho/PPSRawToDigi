#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <FWCore/Utilities/interface/Exception.h>

/**
 *
 */

// TOTEM =7, T2 = 2
T2DetId::T2DetId():DetId(DetId::Totem,2),trind(0){}

/**
 *
 */

T2DetId::T2DetId(uint32_t id):DetId(id),trind(0) {

  //  std::cout<<" constructor of the T2DetId" <<std::endl;
  if (det()!=DetId::Totem || subdetId()!=2) {

    throw cms::Exception("InvalidDetId") << "T2DetId ctor:"
					 << " det: " << det()
					 << " subdet: " << subdetId()
					 << " is not a valid T2 id";  
  }

} // T2DetId

/**
 *
 */

T2DetId::T2DetId(unsigned int arm, unsigned halfTelescope, 
		 unsigned int plane, unsigned int planeSide):	      
  DetId(DetId::Totem,2),trind(0) {
  
  this->init(arm, halfTelescope, plane, planeSide);

} // T2DetId

/**
 *
 */

void T2DetId::init(unsigned int arm, unsigned int halfTelescope, 
		   unsigned int plane, unsigned int planeSide) {
  
  id_ = calculateRawId(arm, halfTelescope, plane, planeSide);
  
} // T2DetId

/**
 *
 */

uint32_t T2DetId::calculateRawId(unsigned int arm, unsigned int halfTelescope, 
				 unsigned int plane, unsigned int planeSide) {
  
  if ( arm>1 || halfTelescope>1 ||  plane>4 || planeSide>1 ) {
    
    throw cms::Exception("InvalidDetId") << "T2DetId ctor:" 
					 << " Invalid parameters: " 
					 << " Arm "<< arm
					 << " HalfTelescope "<< halfTelescope
					 << " Plane "<< plane
                                         << " PlaneSide " << planeSide
					 << std::endl;
    
  }
  
  return 0x74000000 | (planeSide<<19) | (plane<<20) | (halfTelescope<<23) | (arm<<24);
  
} // calculateRawId

/**
 *
 */

std::ostream& operator<<( std::ostream& os, const T2DetId& id ){


  os << " Arm " << id.arm()
     << " HalfTelescope " << id.halfTelescope()
     << " Plane " << id.plane()
     << " PlaneSide " << id.planeSide();
  
  return os;

} // operator<<
