#include "SimG4Core/TotemRPProtTransp/interface/TotemRPParametrizedPhysics.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "G4Electron.hh"
#include "G4FastSimulationManagerProcess.hh"
#include "G4ProcessManager.hh"
#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4LogicalVolumeStore.hh"
#include "SimG4Core/TotemRPProtTransp/interface/BeamProtTransportSetup.h"
#include <thread>
#include "FWCore/ParameterSet/interface/ParameterSet.h"


using namespace CLHEP;

//G4ThreadLocal bool TotemRPParametrizedPhysics::fInitialized = false;

TotemRPParametrizedPhysics::TotemRPParametrizedPhysics(std::string name, const edm::ParameterSet & p)
        :  G4VPhysicsConstructor(name),
           parameters(p)
{}

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

void TotemRPParametrizedPhysics::ConstructProcess()
{
//    if(fInitialized) { return; }
//    fInitialized = true;

    BeamProtTransportSetup* _setup = new BeamProtTransportSetup(parameters);

    edm::LogInfo("TotemRPParametrizedPhysics")
    << "TotemRPParametrizedPhysics: adding the FastSimulationManagerProcess"
    << _setup << std::endl;

    edm::LogInfo("TotemRPParametrizedPhysics")<< "j0" << std::endl;
    G4FastSimulationManagerProcess * theFastSimulationManagerProcess =
            new G4FastSimulationManagerProcess("TotemRPParameterisationProcess", fParameterisation);

    edm::LogInfo("TotemRPParametrizedPhysics")<< "j1" << std::endl;

    aParticleIterator->reset();

    while ((*aParticleIterator)())
    {
        G4ParticleDefinition * particle = aParticleIterator->value();
        G4ProcessManager * pmanager = particle->GetProcessManager();
        pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
    }
}
