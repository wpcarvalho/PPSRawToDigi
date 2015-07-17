#ifndef T2HitCollectionMapping_h
#define T2HitCollectionMapping_h          //not very clear for me

/** 
 * Class T2HitCollection
 * 
 * author Mirko Berretti
 */

#include <vector>
#include <map>
#include <DataFormats/T2Hit/interface/T2Hit.h>
//typedef std::map<int,double> T2HitMapIdToZ;
typedef std::map<std::pair<long int,long int>, unsigned long int> T2HitCollectionMapping;
#endif
