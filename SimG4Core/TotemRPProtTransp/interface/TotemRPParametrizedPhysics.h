#ifndef SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__
#define SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__

#include "G4VPhysicsConstructor.hh"

class TotemRPParametrizedPhysics : public G4VPhysicsConstructor
{
public:
        TotemRPParametrizedPhysics(std::string name);
        virtual ~TotemRPParametrizedPhysics();
        virtual void ConstructParticle();
        virtual void ConstructProcess();
        void addParametrisation();

protected:
        static G4ThreadLocal bool fInitialized;
        std::string processDefFilePath;
};

#endif  //SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__

