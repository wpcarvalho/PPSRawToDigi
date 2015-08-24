#include "QGSPCMS_BERT_EMV.hh"
#include "SimG4Core/PhysicsLists/interface/CMSMonopolePhysics.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "G4EmStandardPhysics_option1.hh"
#include "G4DecayPhysics.hh"
#include "G4EmExtraPhysics.hh"
#include "G4IonPhysics.hh"
#include "G4StoppingPhysics.hh"
#include "G4HadronElasticPhysics.hh"
#include "G4NeutronTrackingCut.hh"
#include "G4HadronicProcessStore.hh"

#include "G4DataQuestionaire.hh"
#include "G4HadronPhysicsQGSP_BERT.hh"

#include "SimG4Core/TotemRPProtTransp/interface/BeamProtTransportSetup.h"
#include "SimG4Core/TotemRPProtTransp/interface/LogicalVolumeStoreFix.h"
#include "SimG4Core/TotemRPProtTransp/interface/TotemRPParametrizedPhysics.h"

#include "G4LogicalVolumeStore.hh"
#include <thread>
#include "G4RunManagerKernel.hh"
#include "G4VPhysicalVolume.hh"

QGSPCMS_BERT_EMV::QGSPCMS_BERT_EMV(G4LogicalVolumeToDDLogicalPartMap& map, 
			   const HepPDT::ParticleDataTable * table_,
			   sim::ChordFinderSetter *chordFinderSetter_, 
			   const edm::ParameterSet & p) : PhysicsList(map, table_, chordFinderSetter_, p){

  G4DataQuestionaire it(photon);

  int  ver     = p.getUntrackedParameter<int>("Verbosity",0);
  bool emPhys  = p.getUntrackedParameter<bool>("EMPhysics",true);
  bool hadPhys = p.getUntrackedParameter<bool>("HadPhysics",true);
  bool tracking= p.getParameter<bool>("TrackingCut");
  edm::LogInfo("PhysicsList") << "You are using the simulation engine: "
			      << "QGSP_BERT_EMV \n Flags for EM Physics "
			      << emPhys << ", for Hadronic Physics "
			      << hadPhys << " and tracking cut " << tracking;
  // EM Physics
  RegisterPhysics( new G4EmStandardPhysics_option1(ver));

  // Synchroton Radiation & GN Physics
  G4EmExtraPhysics* gn = new G4EmExtraPhysics(ver);
  RegisterPhysics(gn);

  // Decays
  this->RegisterPhysics( new G4DecayPhysics(ver) );

  G4HadronicProcessStore::Instance()->SetVerbose(ver);

  // Hadron Elastic scattering
  RegisterPhysics( new G4HadronElasticPhysics(ver));

  // Hadron Physics
  //todo set G4bool quasiElastic=true;
  RegisterPhysics(  new G4HadronPhysicsQGSP_BERT(ver));

  // Stopping Physics
  RegisterPhysics( new G4StoppingPhysics(ver));

  // Ion Physics
  RegisterPhysics( new G4IonPhysics(ver));

  // Neutron tracking cut
  RegisterPhysics( new G4NeutronTrackingCut(ver));

  // Monopoles
//  RegisterPhysics( new CMSMonopolePhysics(table_,chordFinderSetter_,p));

  // Custom Physics
  LogicalVolumeStoreFix::copyToSingleton(G4LogicalVolumeStore::GetInstance(),
                                         G4RunManagerKernel::GetRunManagerKernel());
  if (beam_prot_transp_setup_ == 0) beam_prot_transp_setup_ = new BeamProtTransportSetup(p);
  RegisterPhysics(new TotemRPParametrizedPhysics());
}

QGSPCMS_BERT_EMV::~QGSPCMS_BERT_EMV()
{
  if (beam_prot_transp_setup_!=0)
    delete beam_prot_transp_setup_;
}
