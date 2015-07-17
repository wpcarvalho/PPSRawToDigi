#ifndef DataFormats_T1Road_H
#define DataFormats_T1Road_H

#include <vector>
#include <iterator>

#include <DataFormats/T1Road/interface/T1RecHitGlobal.h>



typedef std::vector<T1RecHitGlobal> T1Road;

typedef std::vector< std::vector<T1RecHitGlobal> > T1RoadCollection;



#endif

