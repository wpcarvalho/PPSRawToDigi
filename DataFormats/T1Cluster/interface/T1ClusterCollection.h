#ifndef T1Cluster_T1ClusterCollection_h
#define T1Cluster_T1ClusterCollection_h

#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <DataFormats/T1Cluster/interface/T1Cluster.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T1DetId, T1Cluster> T1ClusterCollection;

#endif
