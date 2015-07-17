#ifndef T2DigiVfat_T2DigiVfatCollection_h
#define T2DigiVfat_T2DigiVfatCollection_h


#include <DataFormats/T2DetId/interface/T2DetId.h>
#include <DataFormats/T2DigiVfat/interface/T2DigiVfat.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T2DetId, T2DigiVfat> T2DigiVfatCollection;

#endif
