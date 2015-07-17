#include <DataFormats/T1RecHit/interface/T1RecHit2D.h>
#include <iostream>

T1RecHit2D::T1RecHit2D() :
  theDetId(),
  theLocalPosition(0.,0.), 
  theLocalError(0.,0.,0.),

  theChi2( -1. )
  
{}

T1RecHit2D::T1RecHit2D( const T1DetId& id, 
			const LocalPoint& pos, 
			const LocalError& err, 
	
			float chi2
			) :
  theDetId( id ), 
  theLocalPosition( pos ), 
  theLocalError( err ),

  theChi2( chi2 )
 
{}

T1RecHit2D::~T1RecHit2D() {}


std::ostream& operator<<(std::ostream& os, const T1RecHit2D& rh) {
  os << "T1RecHit2D: local x = " << rh.localPosition().x() << " +/- " << sqrt( rh.localPositionError().xx() ) <<
    " y = " << rh.localPosition().y() << " +/- " << sqrt( rh.localPositionError().yy() ) <<
    " chi2 = " << rh.chi2();
  return os;
}

