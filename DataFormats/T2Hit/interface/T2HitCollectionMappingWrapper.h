#ifndef T2HitCollectionMappingWrapper_h
#define T2HitCollectionMappingWrapper_h          //not very clear for me

/** 
 * Class T2StripClusterCollectionWrapper
 * 
 * author Mirko Berretti
 */

#include <vector>
#include <DataFormats/T2Hit/interface/T2Hit.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T2Hit/interface/T2HitCollectionMapping.h>

typedef edm::Wrapper<std::vector<T2HitCollectionMapping> > T2HitCollectionMappingWrapper;


#endif
