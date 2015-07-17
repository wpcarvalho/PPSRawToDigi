#ifndef DATAFORMATS_T1_DET_ID_HH_
#define DATAFORMATS_T1_DET_ID_HH_

#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include <map>

namespace{ 
  namespace {

    T1DetId d;
    std::pair<T1DetId,std::pair<unsigned int,unsigned int> > dummy1;
    
    edm::Wrapper<std::pair<T1DetId,std::pair<unsigned int,unsigned int> > > dummy2;
  }
}

#endif // DATAFORMATS_T1_DET_ID_HH_