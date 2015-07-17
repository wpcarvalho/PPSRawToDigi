#ifndef T2PadClusterCollectionWrapper_h
#define T2PadClusterCollectionWrapper_h          //not very clear for me

/** 
 * Class T2PadClusterCollection
 * 
 * author Mirko Berretti
 */
#include <vector>
#include <map>
#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <DataFormats/T2Cluster/interface/T2Cluster.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T2Cluster/interface/T2PadClusterCollection.h>
//#include <DataFormats/TotemData/interface/TotemDigiCollection.h>


typedef edm::Wrapper<std::vector<T2PadClusterCollection> > T2PadClusterCollectionWrapper;
//typedef TotemDigiCollection<T2DetId,std::vector<T2Cluster> > T2PadClusterCollection
#endif
