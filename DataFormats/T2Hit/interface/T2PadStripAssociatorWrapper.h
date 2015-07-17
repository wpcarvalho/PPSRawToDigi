#ifndef T2PadStripAssociatorWrapper_h
#define T2PadStripAssociatorWrapper_h          

/** 
 * Class T2StripClusterCollectionWrapper
 * 
 * author Mirko Berretti
 */

#include <vector>
#include <DataFormats/T2Hit/interface/T2Hit.h>
#include <DataFormats/Common/interface/Wrapper.h>
#include <DataFormats/T2Hit/interface/T2PadStripAssociator.h>

typedef edm::Wrapper<std::vector<T2PadStripAssociator> > T2PadStripAssociatorWrapper;


#endif
