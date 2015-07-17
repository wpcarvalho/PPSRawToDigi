
#include <DataFormats/T1Road/interface/T1Road.h>



#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/CLHEP/interface/AlgebraicObjects.h>
#include <vector>

namespace{ 
  namespace {

    DeepCopyPointer<AlgebraicSymMatrix> bbbhhjj;
    GlobalError gge;
    T1RecHitGlobal rhG;
    std::vector<T1RecHitGlobal> vrhG;
    std::vector< std::vector<T1RecHitGlobal> > vvrhG;
    T1Road rr;
    T1RoadCollection rrC;

    edm::Wrapper<T1RecHitGlobal> rhgW;
    edm::Wrapper<T1Road> rrW;
    edm::Wrapper<T1RoadCollection> rccW;
  
  }
}

//namespace{ 
//  namespace {
//    T1SegmentCollection seg;    
//    edm::Wrapper<T1SegmentCollection> dwc1;
//  }
//}

