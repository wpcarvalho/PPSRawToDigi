#ifndef SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__
#define SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__

//#define G4V7

#include "G4VPhysicsConstructor.hh"
 
// Joanna Weng 08.2005
// Physics process for Gflash parameterisation
 
class TotemRPParametrizedPhysics : public G4VPhysicsConstructor
{
        public:
        TotemRPParametrizedPhysics(std::string name);
        virtual ~TotemRPParametrizedPhysics();
 
        protected:
        virtual void ConstructParticle();
        virtual void ConstructProcess();
        void addParametrisation();
};
 

#endif  //SimG4Core_TotemRPProtTransp_TotemRPParametrizedPhysics_h__

