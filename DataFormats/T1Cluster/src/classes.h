#include <vector>
#include <map>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T1Cluster/interface/T1Cluster.h>
#include <DataFormats/T1Cluster/interface/T1ClusterCollection.h>

namespace{ 
  namespace {

    T1Cluster d;
    std::vector<T1Cluster>  vv;
    std::vector<std::vector<T1Cluster> >  v1; 
    std::pair<T1DetId, std::vector<T1Cluster> > pv;
    T1ClusterCollection dd;
    
    edm::Wrapper<T1ClusterCollection> dw;
    edm::Wrapper<std::pair<T1DetId, std::vector<T1Cluster> > > pw;

  }
}
