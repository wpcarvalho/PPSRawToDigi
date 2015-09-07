//
// Joanna Weng 08.2005
// Physics process for Gflash parameterisation
// modified by Soon Yung Jun, Dongwook Jang
// V.Ivanchenko rename the class, cleanup, and move
//              to SimG4Core/Application - 2012/08/14

#include "SimG4Core/Application/interface/TotemRPParametrizedPhysics.h"
#include "SimG4Core/Application/interface/GFlashEMShowerModel.h"
#include "SimG4Core/Application/interface/GFlashHadronShowerModel.h"
#include "SimG4Core/Application/interface/ElectronLimiter.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4FastSimulationManagerProcess.hh"
#include "G4ProcessManager.hh"

#include "G4LeptonConstructor.hh"
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4ShortLivedConstructor.hh"
#include "G4IonConstructor.hh"
#include "G4RegionStore.hh"
#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4MuonMinus.hh"
#include "G4MuonPlus.hh"
#include "G4PionMinus.hh"
#include "G4PionPlus.hh"
#include "G4KaonMinus.hh"
#include "G4KaonPlus.hh"
#include "G4Proton.hh"
#include "G4AntiProton.hh"

#include "G4MuonNuclearProcess.hh"
#include "G4MuonVDNuclearModel.hh"

#include "G4EmProcessOptions.hh"
#include "G4PhysicsListHelper.hh"
#include "G4SystemOfUnits.hh"
#include "G4UAtomicDeexcitation.hh"
#include "G4LossTableManager.hh"
#include "SimG4Core/Application/interface/BeamProtTransportSetup.h"


G4ThreadLocal G4FastSimulationManagerProcess* theFastSimulationManagerProcess = 0;
G4ThreadLocal BeamProtTransportSetup* beam_prot_transp_setup_ = 0;

TotemRPParametrizedPhysics::TotemRPParametrizedPhysics(std::string name,
					     const edm::ParameterSet & p) 
  : G4VPhysicsConstructor(name), theParSet(p) 
{
}

TotemRPParametrizedPhysics::~TotemRPParametrizedPhysics() {
}

void TotemRPParametrizedPhysics::ConstructParticle()
{
  edm::LogInfo("TotemRPParametrizedPhysics") << "ConstructParticle";
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
  edm::LogInfo("TotemRPParametrizedPhysics") << "ConstructProcess";

  if(beam_prot_transp_setup_ == 0)
    beam_prot_transp_setup_ = new BeamProtTransportSetup(theParSet);

  if(theFastSimulationManagerProcess ==0)
    theFastSimulationManagerProcess =
          new G4FastSimulationManagerProcess("TotemRPParameterisationProcess", fParameterisation);

  aParticleIterator->reset();

  while ((*aParticleIterator)())
  {
    G4ParticleDefinition * particle = aParticleIterator->value();
    G4ProcessManager * pmanager = particle->GetProcessManager();
    pmanager->AddProcess(theFastSimulationManagerProcess, -1, 1, 1);
  }
}
