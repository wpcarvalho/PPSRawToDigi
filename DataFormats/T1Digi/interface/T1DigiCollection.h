#ifndef T1Digi_T1DigiCollection_h
#define T1Digi_T1DigiCollection_h
/** \class T1DigiCollection
 *  
 *  \author F.Ferro
 */

#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <DataFormats/T1Digi/interface/T1Digi.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T1DetId, T1Digi> T1DigiCollection;

#endif
