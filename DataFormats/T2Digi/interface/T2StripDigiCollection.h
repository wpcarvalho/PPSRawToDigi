#ifndef T2StripDigi_T2StripDigiCollection_h
#define T2StripDigi_T2StripDigiCollection_h

/** 
 * Class T2StripDigiCollection
 * 
 * Author: Erik Brücken / University of Helsinki
 * Email: brucken@cc.helsinki.fi
 */

#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <DataFormats/T2Digi/interface/T2StripDigi.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T2DetId, T2StripDigi> T2StripDigiCollection;

#endif
