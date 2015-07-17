#ifndef DataFormats_T1RecHit2DCollection_H
#define DataFormats_T1RecHit2DCollection_H

#include <vector>
#include <map>
#include <iterator>


#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <DataFormats/T1RecHit/interface/T1RecHit2D.h>

#include <DataFormats/Common/interface/RangeMap.h>
#include <DataFormats/Common/interface/ClonePolicy.h>
#include <DataFormats/Common/interface/OwnVector.h>

typedef edm::RangeMap <T1DetId, edm::OwnVector<T1RecHit2D> > T1RecHit2DCollection;





#endif

