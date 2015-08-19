//
// Created by lichon on 18.08.15.
//

#include "SimG4Core/TotemRPProtTransp/interface/LogicalVolumeStoreFix.h"
#include "G4RegionStore.hh"
#include <thread>
#include "FWCore/MessageLogger/interface/MessageLogger.h"


void LogicalVolumeStoreFix::copyToSingleton(G4LogicalVolumeStore* theStore) {
    G4LogicalVolumeStore::const_iterator it;
    for (it = theStore->begin(); it != theStore->end(); it++)
    {
        G4LogicalVolume * v = *it;
        G4LogicalVolumeStore::GetInstance()->Register(v);
    }
}