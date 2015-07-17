#ifndef T2StripClusterCollectionWrapper_h
#define T2StripClusterCollectionWrapper_h          //not very clear for me

/** 
 * Class T2StripClusterCollectionWrapper
 * 
 * author Mirko Berretti
 */
#include <vector>
#include <map>
#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <DataFormats/T2Cluster/interface/T2Cluster.h>
//#include <DataFormats/TotemData/interface/TotemDigiCollection.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T2Cluster/interface/T2StripClusterCollection.h>

typedef edm::Wrapper<std::vector<T2StripClusterCollection> > T2StripClusterCollectionWrapper;


#endif
