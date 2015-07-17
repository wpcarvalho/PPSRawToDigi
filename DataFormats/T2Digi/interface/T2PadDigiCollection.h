#ifndef T2PadDigi_T2PadDigiCollection_h
#define T2PadDigi_T2PadDigiCollection_h

/** 
 * Class T2PadDigiCollection
 * 
 * Author: Erik Brücken / University of Helsinki
 * Email: brucken@cc.helsinki.fi
 */

#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <DataFormats/T2Digi/interface/T2PadDigi.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T2DetId, T2PadDigi> T2PadDigiCollection;

#endif
