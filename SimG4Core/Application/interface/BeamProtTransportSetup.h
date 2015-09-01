/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*    Hubert Niewiadomski
*    Jan Kašpar (jan.kaspar@gmail.com)
*
* $$RCSfile: BeamProtTransportSetup.h,v $: $
* $Revision: 1.1.1.1.6.2 $
* $Date: 2009/11/16 16:55:54 $
*
****************************************************************************/

#ifndef SimG4Core_TotemRPPRotTransp_BeamProtTransportSetup_h_
#define SimG4Core_TotemRPPRotTransp_BeamProtTransportSetup_h_

//#define G4V7

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimG4Core/Application/interface/ProtTranspFastSimModel.h"
#include "FWCore/Framework/interface/ESWatcher.h"

#include "G4String.hh"
#include "G4ThreeVector.hh"

namespace edm {
  class EventSetup;
  class Run;
}

class ProtonTransportRcd;

/**
 * \brief Singleton class to build and maintain proton transport infrastructure.
 * Creates fast simulation model for fast proton transport (ProtTranspFastSimModel)
 * for the following beam pipe segments:
 *  - Beam_IP_150_R
 *  - Beam_IP_150_L
 * defined in detector xml description
**/
class BeamProtTransportSetup
{
public:
    BeamProtTransportSetup(edm::ParameterSet const & p);
    ~BeamProtTransportSetup();

    void UpdateSetup(const edm::EventSetup &);

  private:
    edm::ParameterSet m_pBeamProtTransportSetup;

    // in m
    double model_ip_150_r_zmin;
    double model_ip_150_r_zmax;
    double model_ip_150_l_zmin;
    double model_ip_150_l_zmax;

    ProtTranspFastSimModel *model_ip_150_r;
    ProtTranspFastSimModel *model_ip_150_l;
    
    G4LogicalVolume *Beam_IP_150_R_LV;
    G4LogicalVolume *Beam_IP_150_L_LV;
    
    G4String Beam_IP_150_R_LV_Name;
    G4String Beam_IP_150_L_LV_Name;

    bool verbosity_;

    void BuildTransportModels(const edm::ParameterSet & p);
    void FindLogicalVolumes();
};


#endif  //SimG4Core_TotemRPPRotTransp_BeamProtTransportSetup_h_

