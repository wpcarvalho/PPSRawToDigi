#ifndef DataFormats_T2Road_Classes_H
#define DataFormats_T2Road_Classes_H

#include <vector>
#include <string>
#include <DataFormats/Common/interface/Wrapper.h>

#include "DataFormats/T2Road/interface/T2Road.h"
#include "DataFormats/T2Road/interface/T2RoadCollection.h"
#include "DataFormats/T2Road/interface/T2RoadCollectionWrapper.h"



namespace{ 
  namespace {
    T2Road t2road_;
    T2RoadCollection t2roadc_;
    T2RoadCollectionWrapper T2RoadcolW_;
    std::vector<T2Road> dummyt2roadc_;	
    //  std::vector<std::vector<T2Hit> > dummyT2Roadc_;	
    //  std::vector<std::vector<T2Road> > t2roadcoll_;
    // edm::Wrapper<T2Road> t2rW_;
    //edm::Wrapper<std::vector<std::vector<T2Hit> > > t2roadcW_;	
    edm::Wrapper<std::vector<T2Road > > t2roadcW_;
  }
}

#endif // DataFormats_T2Hit_Classes_H
