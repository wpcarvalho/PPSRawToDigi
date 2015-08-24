//
// Created by lichon on 18.08.15.
//

#include "SimG4Core/TotemRPProtTransp/interface/LogicalVolumeStoreFix.h"
#include "G4RegionStore.hh"
#include <thread>
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4RunManagerKernel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

void LogicalVolumeStoreFix::copyToSingleton(G4LogicalVolumeStore* theStore, G4RunManagerKernel* theKernel) {
    G4LogicalVolumeStore::const_iterator it;
    for (it = theStore->begin(); it != theStore->end(); it++)
    {
        G4LogicalVolume * v = *it;
        G4LogicalVolumeStore::GetInstance()->Register(v);
    }
//    edm::LogInfo("PhysicsList") << "q44";
//    G4VPhysicalVolume* world = theKernel->GetCurrentWorld();
//    edm::LogInfo("PhysicsList") << "q55" << world << " " << G4RunManagerKernel::GetRunManagerKernel() ;
//    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetWorldVolume(world);
//    edm::LogInfo("PhysicsList") << "q66";
}