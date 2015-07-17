#include <DataFormats/T1DigiWire/interface/T1DigiWire.h>
#include <DataFormats/T1DigiWire/interface/T1DigiWireCollection.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>
#include <map>

namespace{ 
  namespace {

    T1DigiWire d;
    std::vector<T1DigiWire>  vv;
    std::vector<std::vector<T1DigiWire> >  v1;
    std::pair<T1DetId,std::vector<T1DigiWire> > pv;
    T1DigiWireCollection dd;
    
    edm::Wrapper<T1DigiWireCollection> dw;
    edm::Wrapper<std::pair<T1DetId,std::vector<T1DigiWire> > > pw;
  }
}
