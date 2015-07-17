#ifndef T1DigiSantiard_T1DigiSantiardCollection_h
#define T1DigiSantiard_T1DigiSantiardCollection_h
/** \class T1DigiSantiardCollection
 *  
 *  \author F.Ferro
 */

#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <DataFormats/T1DigiSantiard/interface/T1DigiSantiard.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T1DetId, T1DigiSantiard> T1DigiSantiardCollection;

#endif
