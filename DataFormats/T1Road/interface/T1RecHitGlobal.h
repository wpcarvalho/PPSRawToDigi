#ifndef DataFormats_T1RecHitGlobal_H
#define DataFormats_T1RecHitGlobal_H

#include "DataFormats/GeometryVector/interface/GlobalPoint.h"
#include "DataFormats/GeometryCommonDetAlgo/interface/GlobalError.h"
#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <vector>
#include <map>
#include <iosfwd>

class T1RecHitGlobal{

 public:



  T1RecHitGlobal();
 
  explicit T1RecHitGlobal( const GlobalPoint& pos, const GlobalError& err );
  explicit T1RecHitGlobal(const T1DetId& id, const GlobalPoint& pos, const GlobalError& err );

  ~T1RecHitGlobal();


  T1RecHitGlobal* clone() const { return new T1RecHitGlobal( *this ); }
  GlobalPoint GlobalPosition() const { return theGlobalPosition; }
  GlobalError GlobalPositionError() const { return theGlobalError; }

  float Eta(float,float,float) const;
  float Phi(float,float) const;
  float eta() const;
  float phi() const;
  T1DetId Id() const {return theId;}

 private:
 
  GlobalPoint theGlobalPosition;
  GlobalError theGlobalError;
  T1DetId theId;

  
};

std::ostream& operator<<(std::ostream& os, const T1RecHitGlobal& rh);

#endif
