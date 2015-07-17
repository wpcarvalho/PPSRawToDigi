#ifndef T1DigiWire_T1DigiWireCollection_h
#define T1DigiWire_T1DigiWireCollection_h
/** \class T1DigiWireCollection
 *  
 *  \author F.Ferro
 */

#include <DataFormats/T1DetId/interface/T1DetId.h>
#include <DataFormats/T1DigiWire/interface/T1DigiWire.h>
#include <DataFormats/TotemData/interface/TotemDigiCollection.h>

typedef TotemDigiCollection<T1DetId, T1DigiWire> T1DigiWireCollection;

#endif
