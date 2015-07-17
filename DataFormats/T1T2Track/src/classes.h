#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T1T2Track/interface/T1T2Track.h>
#include <DataFormats/T1T2Track/interface/T1T2TrackCollection.h>
#include <DataFormats/T1Road/interface/T1Road.h>
#include <DataFormats/T2Road/interface/T2Road.h>
#include <DataFormats/T2Road/interface/T2RoadCollection.h>



#include <vector>

namespace{
  namespace {
    T1RecHitGlobal rhG;
    T2Hit t2h;
    T1Road rr;
    T1RoadCollection rrC;
    T2Road rr2;
    T2RoadCollection rrC2;
    T1T2Track tt;
    T1T2TrackCollection ttC;
    edm::Wrapper<T1T2TrackCollection> tccW;
  }
}
