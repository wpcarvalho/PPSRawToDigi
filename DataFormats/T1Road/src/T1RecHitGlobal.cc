#include <DataFormats/T1Road/interface/T1RecHitGlobal.h>
#include <iostream>

T1RecHitGlobal::T1RecHitGlobal() :
  
  theGlobalPosition(0.,0.,0.), 
  theGlobalError(0.,0.,0.,0.,0.,0.)


{}

T1RecHitGlobal::T1RecHitGlobal(
			       const GlobalPoint& pos, 
			       const GlobalError& err
			       ) :  
  theGlobalPosition( pos ), 
  theGlobalError( err )
 
{}

T1RecHitGlobal::T1RecHitGlobal(const T1DetId& id,
			       const GlobalPoint& pos, 
			       const GlobalError& err
			       ) :  
  theGlobalPosition( pos ), 
  theGlobalError( err ),
  theId(id)
 
{}

T1RecHitGlobal::~T1RecHitGlobal() {}

float T1RecHitGlobal::eta() const {
  return Eta(theGlobalPosition.x(),theGlobalPosition.y(),theGlobalPosition.z());
}

float T1RecHitGlobal::phi() const {
  return Phi(theGlobalPosition.x(),theGlobalPosition.y());
}


float T1RecHitGlobal::Eta(float x,float y,float z) const {
  float xyt=0;
  float c=0;
  float eta2=0;
  xyt = sqrt(x*x + y*y);
  //theta
  if(z>0) c = atan(xyt/z);
  if(z<0) c = atan(xyt/z)+3.14159;
  if(z==0) {c = 3.14159;}
  //pseudorapidity
  eta2 = -log(tan(c/2.));
  return eta2;
}
float T1RecHitGlobal::Phi(float x,float y) const {
  float c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}

std::ostream& operator<<(std::ostream& os, const T1RecHitGlobal& rh) {
  os << "T1RecHitGlobal: global x = " << rh.GlobalPosition().x() << " +/- " << sqrt( rh.GlobalPositionError().cxx() ) <<
    " y = " << rh.GlobalPosition().y() << " +/- " << sqrt( rh.GlobalPositionError().cyy() ) 
     <<" z = " <<  rh.GlobalPosition().z() << " +/- " << sqrt( rh.GlobalPositionError().czz() ) << " Eta = " << rh.Eta(rh.GlobalPosition().x(), rh.GlobalPosition().y(),rh.GlobalPosition().z())
     << " Phi = " << rh.Phi(rh.GlobalPosition().x(), rh.GlobalPosition().y() );
   
  return os;
}
