#include "SimG4Core/PhysicsLists/interface/CMSEmStandardPhysics.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimG4Core/TotemRPProtTransp/interface/TotemRPPhysicsList.h"
#include "SimG4Core/TotemRPProtTransp/interface/BeamProtTransportSetup.h"
#include "SimG4Core/TotemRPProtTransp/interface/TotemRPParametrizedPhysics.h"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4QStoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh" 
#include "G4NeutronTrackingCut.hh"

#include "G4DataQuestionaire.hh"
#include "HadronPhysicsQGSP.hh"

#include "SimG4Core/Physics/interface/PhysicsListFactory.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
TotemRPPhysicsList::TotemRPPhysicsList(G4LogicalVolumeToDDLogicalPartMap & map, const HepPDT::ParticleDataTable * table_,
	      sim::FieldBuilder *fieldBuilder_, const edm::ParameterSet & p) :
	      PhysicsList(map, table_, fieldBuilder_, p), beam_prot_transp_setup_(0) {

  G4DataQuestionaire it(photon);

  edm::LogInfo("PhysicsList") << "You are using the simulation engine: " 
                              << "QGSP_EMV 3.3 "
                              << " + TOTEM Fast Proton Beam Transport "
                              << "\n";

  // EM Physics
    RegisterPhysics( new CMSEmStandardPhysics("standard EM v?",0));

  // Synchroton Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics("extra EM"));

  // Decays
  RegisterPhysics(new G4DecayPhysics("decay"));

  // Hadron Elastic scattering
  RegisterPhysics(new G4HadronElasticPhysics("elastic",0,false)); 

  // Hadron Physics
  G4bool quasiElastic=true;
  RegisterPhysics(new HadronPhysicsQGSP("hadron",quasiElastic));
  //RegisterPhysics(new HadronPhysicsQGSP("hadron"));

  // Stopping Physics
  RegisterPhysics(new G4QStoppingPhysics("stopping"));

  // Ion Physics
  RegisterPhysics(new G4IonPhysics("ion"));

  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut("Neutron tracking cut", 0));

  // Custom Physics
  if (beam_prot_transp_setup_==0) beam_prot_transp_setup_ = new BeamProtTransportSetup(p);
  RegisterPhysics(new TotemRPParametrizedPhysics("totem_parametrised_prot_transp"));
}

TotemRPPhysicsList::~TotemRPPhysicsList()
{ 
  if (beam_prot_transp_setup_!=0)
    delete beam_prot_transp_setup_; 
}

DEFINE_PHYSICSLIST(TotemRPPhysicsList);
