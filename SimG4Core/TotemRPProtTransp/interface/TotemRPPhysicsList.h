#ifndef SimG4Core_TotemRPProtTransp_TotemRPPhysicsList_h__
#define SimG4Core_TotemRPProtTransp_TotemRPPhysicsList_h__

//#define G4V7

//TotemRPPhysicsList


//setup of beam volumes for fast transport by parametrization
 
#include "SimG4Core/Physics/interface/PhysicsList.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "SimG4Core/TotemRPProtTransp/interface/BeamProtTransportSetup.h"
 
class TotemRPPhysicsList : public PhysicsList
{
public:
    TotemRPPhysicsList(G4LogicalVolumeToDDLogicalPartMap & map, const HepPDT::ParticleDataTable * table_,
  	      sim::FieldBuilder *fieldBuilder_, const edm::ParameterSet & p);
    virtual ~TotemRPPhysicsList();
private:
    BeamProtTransportSetup * beam_prot_transp_setup_;
};




#endif  //SimG4Core_TotemRPProtTransp_TotemRPPhysicsList_h__
