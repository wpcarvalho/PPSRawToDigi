#include <DataFormats/T2DigiVfat/interface/T2DigiVfat.h>
#include <DataFormats/T2DigiVfat/interface/T2DigiVfatCollection.h>
#include <DataFormats/T2DigiVfat/interface/T2VfatInformation.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>
#include <map>

namespace{ 
  namespace {

    T2DigiVfat d;
    std::vector<T2DigiVfat>  vv;
    std::vector<std::vector<T2DigiVfat> >  v1; 
    std::pair<T2DetId,std::vector<T2DigiVfat> > pv;
    T2DigiVfatCollection dd;
    
    edm::Wrapper<T2DigiVfatCollection> dw;
    edm::Wrapper<std::pair<T2DetId,std::vector<T2DigiVfat> > > pw;
    
    T2VfatInformation inf;
    edm::Wrapper<T2VfatInformation> infw;  

  }
}
