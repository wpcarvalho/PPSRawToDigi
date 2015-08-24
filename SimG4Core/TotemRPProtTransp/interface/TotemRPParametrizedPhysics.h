#ifndef SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__
#define SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__

#include "G4VPhysicsConstructor.hh"
#include "SimG4Core/TotemRPProtTransp/interface/TotemRPParametrizedPhysics.h"

#include "G4LogicalVolumeStore.hh"
#include <thread>
#include "FWCore/ParameterSet/interface/ParameterSet.h"


class TotemRPParametrizedPhysics : public G4VPhysicsConstructor
{
public:
        TotemRPParametrizedPhysics();
        virtual ~TotemRPParametrizedPhysics();
        virtual void ConstructParticle();
        virtual void ConstructProcess();

protected:
        static G4ThreadLocal bool fInitialized;
        std::string processDefFilePath;
};

#endif  //SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__

