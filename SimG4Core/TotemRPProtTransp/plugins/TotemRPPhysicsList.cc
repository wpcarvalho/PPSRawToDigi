#include "SimG4Core/PhysicsLists/interface/CMSEmStandardPhysics.h"


#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimG4Core/TotemRPProtTransp/interface/TotemRPPhysicsList.h"
#include "SimG4Core/TotemRPProtTransp/interface/BeamProtTransportSetup.h"
#include "SimG4Core/TotemRPProtTransp/interface/TotemRPParametrizedPhysics.h"

#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh" 
#include "G4NeutronTrackingCut.hh"
#include "G4HadronicProcessStore.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_FTFP_BERT.hh"

#include "SimG4Core/Physics/interface/PhysicsListFactory.h"
#include "FWCore/Framework/interface/MakerMacros.h"
 
TotemRPPhysicsList::TotemRPPhysicsList(G4LogicalVolumeToDDLogicalPartMap & map, const HepPDT::ParticleDataTable * table_,
	      sim::ChordFinderSetter *chordFinderSetter_, const edm::ParameterSet & p) :
	      PhysicsList(map, table_, chordFinderSetter_, p), beam_prot_transp_setup_(0) { //todo check sim::FieldBuilder vs sim::ChordFinderSetter

  int  ver     = p.getUntrackedParameter<int>("Verbosity",0);
  G4DataQuestionaire it(photon);

  edm::LogInfo("PhysicsList") << "You are using the simulation engine: " 
                              << "QGSP_EMV 3.3 "
                              << " + TOTEM Fast Proton Beam Transport "
                              << "\n";

  // EM Physics
  RegisterPhysics( new CMSEmStandardPhysics(ver));
  // Synchroton Radiation & GN Physics
  RegisterPhysics(new G4EmExtraPhysics(ver));
  // Decays
  RegisterPhysics(new G4DecayPhysics(ver));
  // Hadron Elastic scattering
  G4HadronicProcessStore::Instance()->SetVerbose(ver);
  RegisterPhysics(new G4HadronElasticPhysics(ver));
  // Hadron Physics
  G4bool quasiElastic=true;
  RegisterPhysics(new G4HadronPhysicsQGSP_FTFP_BERT(ver, quasiElastic));
  // Stopping Physics
  RegisterPhysics(new G4StoppingPhysics(ver));
  // Ion Physics
  RegisterPhysics(new G4IonPhysics(ver));
  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));
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
