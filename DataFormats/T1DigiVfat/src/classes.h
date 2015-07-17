#include <DataFormats/T1DigiVfat/interface/T1DigiVfat.h>
#include <DataFormats/T1DigiVfat/interface/T1DigiVfatCollection.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>
#include <map>

namespace{ 
  namespace {

    T1DigiVfat d;
    std::vector<T1DigiVfat>  vv;
    std::vector<std::vector<T1DigiVfat> >  v1; 
    std::pair<T1DetId,std::vector<T1DigiVfat> > pv;
    T1DigiVfatCollection dd;
    
    edm::Wrapper<T1DigiVfatCollection> dw;
    edm::Wrapper<std::pair<T1DetId,std::vector<T1DigiVfat> > > pw;
  }
}
