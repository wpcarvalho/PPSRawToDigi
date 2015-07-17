#ifndef T2HitCollectionWrapper_h
#define T2HitCollectionWrapper_h          //not very clear for me

/** 
 * Class T2StripClusterCollectionWrapper
 * 
 * author Mirko Berretti
 */

#include <vector>
#include <DataFormats/T2Hit/interface/T2Hit.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T2Hit/interface/T2HitCollection.h>

typedef edm::Wrapper<std::vector<T2HitCollection> > T2HitCollectionWrapper;


#endif
