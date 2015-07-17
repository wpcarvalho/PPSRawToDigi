#ifndef DataFormats_T2Cluster_Classes_H
#define DataFormats_T2Cluster_Classes_H

#include <boost/cstdint.hpp>
#include <vector>
#include <string>
#include <map>

#include <DataFormats/Common/interface/Wrapper.h>
#include "DataFormats/T2Cluster/interface/T2StripClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollection.h"
#include "DataFormats/T2Cluster/interface/T2PadClusterCollectionWrapper.h"
#include "DataFormats/T2Cluster/interface/T2StripClusterCollectionWrapper.h"
#include "DataFormats/T2Cluster/interface/cluster_entry.h"
#include "DataFormats/T2Cluster/interface/T2Cluster.h"
#include "DataFormats/T2Cluster/interface/T2ROGeometry.h"




  namespace {
    T2Cluster t2C_;
    T2StripClusterCollection t2SDC_;
    T2PadClusterCollection  t2PDC_;
    T2StripClusterCollectionWrapper t2SCCW_; 
    T2PadClusterCollectionWrapper t2PCCW_;		 
    std::vector<T2Cluster> myvtc;
    std::pair<T2DetId,std::vector<T2Cluster> > t2pc_;
    cluster_entry clentry;
    std::vector<cluster_entry> clentryvt;
    edm::Wrapper<T2StripClusterCollection> t2sccw_;
    edm::Wrapper<T2PadClusterCollection> t2pccw_;
    edm::Wrapper<std::map<T2DetId,std::vector<T2Cluster> > > t2pccwdum_;
    edm::Wrapper<std::pair<T2DetId,std::vector<T2Cluster> > > t2wpc;
  }


#endif // DataFormats_T2Cluster_Classes_H
