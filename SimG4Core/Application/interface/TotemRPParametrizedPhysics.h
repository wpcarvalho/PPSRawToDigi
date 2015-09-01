#ifndef SimG4Core_TotemRPParametrizedPhysics_H
#define SimG4Core_TotemRPParametrizedPhysics_H

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4VPhysicsConstructor.hh"

class TotemRPParametrizedPhysics : public G4VPhysicsConstructor
{
public:

  TotemRPParametrizedPhysics(std::string name, const edm::ParameterSet & p);
  virtual ~TotemRPParametrizedPhysics();
	
  virtual void ConstructParticle();
  virtual void ConstructProcess();

private:
  edm::ParameterSet theParSet;

};

#endif

