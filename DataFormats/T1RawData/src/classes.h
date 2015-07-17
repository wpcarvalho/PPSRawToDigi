#include <DataFormats/T1RawData/interface/TotemRawVFATFrame.h>
#include <DataFormats/T1RawData/interface/TotemVFATFrameColl.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>
#include <map>
#include <ext/hash_map>

namespace{ 
  namespace {

    __gnu_cxx::hash_map<int,int> fgh;
    TotemRawVFATFrame d;
    std::vector<TotemRawVFATFrame>  vv;
    std::vector<std::vector<TotemRawVFATFrame> >  v1; 
    TotemVFATFrameColl dd;
    
    edm::Wrapper<TotemVFATFrameColl> dw;

  }
}
