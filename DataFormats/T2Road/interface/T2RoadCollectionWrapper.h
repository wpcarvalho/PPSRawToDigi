#ifndef T2RoadCollectionWrapper_h
#define T2RoadCollectionWrapper_h          //not very clear for me

/**
 * Class T2StripClusterCollectionWrapper
 *
 * author Mirko Berretti
 */

#include <vector>
#include <DataFormats/T2Road/interface/T2Road.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T2Road/interface/T2RoadCollection.h>

typedef edm::Wrapper<std::vector<T2RoadCollection> > T2RoadCollectionWrapper;

#endif

