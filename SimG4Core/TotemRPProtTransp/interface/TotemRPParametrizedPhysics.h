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
        TotemRPParametrizedPhysics(std::string name, const edm::ParameterSet & p);
        virtual ~TotemRPParametrizedPhysics();
        virtual void ConstructParticle();
        virtual void ConstructProcess();

protected:
        static G4ThreadLocal bool fInitialized;
        std::string processDefFilePath;
        edm::ParameterSet parameters;
};

#endif  //SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__

