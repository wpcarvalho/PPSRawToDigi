#include <DataFormats/T1DigiSantiard/interface/T1DigiSantiard.h>
#include <DataFormats/T1DigiSantiard/interface/T1DigiSantiardCollection.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>
#include <map>

namespace{ 
  namespace {

    T1DigiSantiard d;
    std::vector<T1DigiSantiard>  vv;
    std::vector<std::vector<T1DigiSantiard> >  v1; 
    T1DigiSantiardCollection dd;
    
    edm::Wrapper<T1DigiSantiardCollection> dw;

  }
}
