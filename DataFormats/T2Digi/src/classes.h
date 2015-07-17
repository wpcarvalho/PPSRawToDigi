#ifndef DataFormats_T2Digi_Classes_H
#define DataFormats_T2Digi_Classes_H

#include <vector>
#include <string>

#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"

#include <DataFormats/Common/interface/Wrapper.h>
#include <vector>
#include <map>

namespace{ 
  
  namespace {
    
    T2StripDigi tSD_;
    T2PadDigi  tPD_;
    
    std::vector<T2StripDigi>  vtSD_;
    std::vector<T2PadDigi>   vtPD_;

    std::vector<std::vector<T2StripDigi> >  vvtSD_; 
    std::pair<T2DetId,std::vector<T2StripDigi> > pvtSD_;
    std::vector<std::vector<T2PadDigi>  >  vvtPD_;
    std::pair<T2DetId,std::vector<T2PadDigi> > pvtPD_;

    T2StripDigiCollection tlSD_;
    T2PadDigiCollection  tlPD_;

    edm::Wrapper<T2StripDigiCollection> wtSD_;
    edm::Wrapper<std::pair<T2DetId,std::vector<T2StripDigi> > > wpSD_;
    edm::Wrapper<T2PadDigiCollection> wtPD_;
    edm::Wrapper<std::pair<T2DetId,std::vector<T2PadDigi> > > wpPD_;
  }
}

#endif // DataFormats_T2Digi_Classes_H
