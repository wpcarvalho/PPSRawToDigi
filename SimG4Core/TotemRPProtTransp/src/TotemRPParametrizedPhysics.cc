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

#include "globals.hh"
#include "G4PhysicsConstructorFactory.hh"
#include "G4TransportationManager.hh"
#include "G4PathFinder.hh"
#include "G4GlobalFastSimulationManager.hh"
#include "G4RunManagerKernel.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Navigator.hh"
#include "G4TransportationManager.hh"

using namespace CLHEP;

G4_DECLARE_PHYSCONSTR_FACTORY(TotemRPParametrizedPhysics);


G4ThreadLocal bool TotemRPParametrizedPhysics::fInitialized = false;

TotemRPParametrizedPhysics::TotemRPParametrizedPhysics()
        :  G4VPhysicsConstructor("totem_parametrised_prot_transp")
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

    edm::LogInfo("TotemRPParametrizedPhysics")
    << "TotemRPParametrizedPhysics: adding the FastSimulationManagerProcess"
    << std::endl;

    edm::LogInfo("TotemRPParametrizedPhysics")<< "j0" << std::endl;

    //todo fix undefined world volume
//    G4VPhysicalVolume* world = QGSPCMS_BERT_EMV::kernel->GetCurrentWorld();
    edm::LogInfo("PhysicsList") << "j0.1" << G4TransportationManager::GetTransportationManager();
    edm::LogInfo("PhysicsList") << "j0.2" << G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking();
//    G4TransportationManager::GetTransportationManager()->GetNavigatorForTracking()->SetWorldVolume(world);

    //todo this crashes
    G4FastSimulationManagerProcess * theFastSimulationManagerProcess =
            new G4FastSimulationManagerProcess();


    edm::LogInfo("TotemRPParametrizedPhysics")<< "j1" << std::endl;

    //todo this is example of PhysicsListHelper usage (a different approach)
// Add Decay Process
//    G4PhysicsListHelper* ph = G4PhysicsListHelper::GetPhysicsListHelper();
//    G4Decay* theDecayProcess = new G4Decay();
//    aParticleIterator->reset();
//    while( (*aParticleIterator)() ){
//        G4ParticleDefinition* particle = aParticleIterator->value();
//        G4ProcessManager* pmanager = particle->GetProcessManager();
//        if (theDecayProcess->IsApplicable(*particle)) {
//            pmanager->AddProcess(theDecayProcess);
//        }
//    }

    aParticleIterator->reset();

    while ((*aParticleIterator)())
    {
        G4ParticleDefinition * particle = aParticleIterator->value();
        G4ProcessManager * pmanager = particle->GetProcessManager();
        pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
    }
}
