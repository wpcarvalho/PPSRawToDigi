#ifndef PPSDigi_PPSTimingDigiCollection_h
#define PPSDigi_PPSTimingDigiCollection_h

/** \class PPSTimingigiCollection
 *  The collection containingPPSTiming Digis in the event.
 *  Digis are grouped by detId.
 *
 *  \author Seyed Mohsen Etesami
 */ April 2016

#include <DataFormats/PPSDetId/interface/PPSTimingDetId.h>
#include <DataFormats/CTPPSDigi/interface/DiamondDigi.h>
#include <DataFormats/CTPPSDigi/interface/PPSDigiCollection.h>

typedef PPSDigiCollection<detId, DiamondDigi> PPSTimingDigiCollection;

#endif

