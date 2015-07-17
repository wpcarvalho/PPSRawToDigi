#ifndef T2PadClusterCollection_h
#define T2PadClusterCollection_h          //not very clear for me

/** 
 * Class T2PadClusterCollection
 * 
 * author Mirko Berretti
 */
#include <vector>
#include <map>
#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <DataFormats/T2Cluster/interface/T2Cluster.h>
//#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef std::map<T2DetId, std::vector<T2Cluster> > T2PadClusterCollection;
//typedef TotemDigiCollection<T2DetId,std::vector<T2Cluster> > T2PadClusterCollection
#endif
