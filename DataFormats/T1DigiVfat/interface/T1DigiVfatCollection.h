#ifndef T1DigiVfat_T1DigiVfatCollection_h
#define T1DigiVfat_T1DigiVfatCollection_h
/** \class T1DigiVfatCollection
 *  
 *  \author F.Ferro
 */

#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <DataFormats/T1DigiVfat/interface/T1DigiVfat.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T1DetId, T1DigiVfat> T1DigiVfatCollection;

#endif
