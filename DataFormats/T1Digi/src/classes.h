#include <DataFormats/T1Digi/interface/T1Digi.h>
#include <DataFormats/T1Digi/interface/T1DigiCollection.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>
#include <map>

namespace{ 
  namespace {

    T1Digi d;
    std::vector<T1Digi>  vv;
    std::vector<std::vector<T1Digi> >  v1; 
    T1DigiCollection dd;
    
    edm::Wrapper<T1DigiCollection> dw;

  }
}
