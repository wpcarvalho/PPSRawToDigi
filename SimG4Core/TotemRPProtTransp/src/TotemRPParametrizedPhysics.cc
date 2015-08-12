//#define G4V7

#include "SimG4Core/TotemRPProtTransp/interface/TotemRPParametrizedPhysics.h"
 
#include "G4Electron.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4ProcessManager.hh"
 
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "FWCore/MessageLogger/interface/MessageLogger.h"
 
TotemRPParametrizedPhysics::TotemRPParametrizedPhysics(std::string name) : G4VPhysicsConstructor(name) {}

TotemRPParametrizedPhysics::~TotemRPParametrizedPhysics() {}
 
void TotemRPParametrizedPhysics::ConstructParticle()
{
    G4LeptonConstructor pLeptonConstructor;
    pLeptonConstructor.ConstructParticle();
 
    G4MesonConstructor pMesonConstructor;
    pMesonConstructor.ConstructParticle();
 
    G4BaryonConstructor pBaryonConstructor;
    pBaryonConstructor.ConstructParticle();
 
    G4ShortLivedConstructor pShortLivedConstructor;
    pShortLivedConstructor.ConstructParticle();
     
    G4IonConstructor pConstructor;
    pConstructor.ConstructParticle();
}
 
void TotemRPParametrizedPhysics::ConstructProcess() { addParametrisation(); }
 
void TotemRPParametrizedPhysics::addParametrisation()
{
    edm::LogInfo("TotemRPParametrizedPhysics") << "TotemRPParametrizedPhysics: adding the FastSimulationManagerProcess" << std::endl;

    G4FastSimulationManagerProcess * theFastSimulationManagerProcess =
        new G4FastSimulationManagerProcess("TotemRPParameterisationProcess", fParameterisation);
    theParticleIterator->reset();
    while ((*theParticleIterator)())
    {
        G4ParticleDefinition * particle = theParticleIterator->value();
        G4ProcessManager * pmanager = particle->GetProcessManager();
        pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
    }
}
