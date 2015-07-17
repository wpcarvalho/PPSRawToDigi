#ifndef DataFormats_T1RecHit2D_H
#define DataFormats_T1RecHit2D_H


#include <DataFormats/TrackingRecHit/interface/RecHit2DLocalPos.h>
#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <vector>
#include <map>
#include <iosfwd>

class T1RecHit2D : public RecHit2DLocalPos {

 public:



  T1RecHit2D();
  //  T1RecHit2D( const DetId& id, const GeomDet* det, 
  T1RecHit2D( const T1DetId& id, 
	      const LocalPoint& pos, const LocalError& err,  
	      float chi2 );
  ~T1RecHit2D();

  /// RecHit2DLocalPos base class interface
  T1RecHit2D* clone() const { return new T1RecHit2D( *this ); }
  LocalPoint localPosition() const { return theLocalPosition; }
  LocalError localPositionError() const { return theLocalError; }

  /// TrackingRecHit base class interface
  DetId geographicalId() const { return theDetId; }
  T1DetId t1DetId() const { return theDetId; }

  /// Probability from fit during rechit build
  //  float prob() const { return theProb; }

  /// Chi-squared from fit during rechit build
  float chi2() const { return theChi2; }

  /// Fitted peaking time
  //  float tpeak() const { return theTpeak; }

  /// Container of strip channel numbers comprising the rechit
  //  const ChannelContainer& channels() const {
  //    return theChaCo;
  //  }




  // To handle global values must use DetId to identify Det, hence Surface, which can transform from local
  // GlobalPoint globalPosition() const;

  //  Useful when building segments...
  //  bool nearby(const T1RecHit2D& other, float maxDeltaRPhi);
  //  bool nearby(float otherX, float maxDeltaRPhi);

 private:
  T1DetId theDetId;
  //  const GeomDet* theDet;
  LocalPoint theLocalPosition;
  LocalError theLocalError;
  
  float theChi2;
  
};

/// Output operator for T1RecHit2D
std::ostream& operator<<(std::ostream& os, const T1RecHit2D& rh);

#endif
