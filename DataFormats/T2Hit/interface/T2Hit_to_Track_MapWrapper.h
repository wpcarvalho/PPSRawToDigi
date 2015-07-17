#ifndef T2Hit_to_Track_MapWrapper_h
#define T2Hit_to_Track_MapWrapper_h          

/** 
 * Class T2StripClusterCollectionWrapper
 * 
 * author Mirko Berretti
 */

#include <vector>
#include <DataFormats/T2Hit/interface/T2Hit.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T2Hit/interface/T2Hit_to_Track_Map.h>

typedef edm::Wrapper<std::vector<T2Hit_to_Track_Map> > T2Hit_to_Track_MapWrapper;


#endif
