#ifndef TotemRPValidation_RPGeant4Validation_RPSUPPLEMENTARYINFO_H_
#define TotemRPValidation_RPGeant4Validation_RPSUPPLEMENTARYINFO_H_

#include <vector>
#include <map>

#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPPSimHitDebugInfo.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPDetSpaceGeometry.h"

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"

class RPSupplementaryInfo
{
  public:
    RPSupplementaryInfo();
    void SetCurrentEvent(const std::vector<RPPSimHitDebugInfo> &in_sim_vect);
    
//    inline bool OriginalProtonTowardsRightGot() const;
//    inline const RPPSimHitDebugInfo &GetOriginalProtonTowardsRight() const;
//    inline bool OriginalProtonTowardsRight220Got() const;
//    inline const RPPSimHitDebugInfo &GetOriginalProtonTowardsRight220() const;
//    inline bool PrimaryProtonOutgoingFrom220RightStation() const;
//    inline const RPPSimHitDebugInfo &GetOutgoingPrimaryProtonFrom220RightStation() const;
//    inline int NumberofPrimaryDesintegratedProtons() const;
//    inline const std::vector<RPPSimHitDebugInfo> & GetPrimaryDesintegratedProtons();
    
    //station independent interface
    inline bool OriginalProtonTowardsStationGot(RPStationId st_id) {return _primary_protons_at_IP[RPDetSpaceGeometry::IsOnTheRight(st_id)].size()>0;}
    inline const RPPSimHitDebugInfo &GetOriginalProtonTowardsStation(RPStationId st_id) {return _primary_protons_at_IP[RPDetSpaceGeometry::IsOnTheRight(st_id)][0];}
    double Get_Log10t_OfOriginalProtonTowardsStation(RPStationId st_id);
    inline bool OriginalProtonAtStationGot(RPStationId st_id) {return _primary_proton_enters_station[st_id].size()>0;}
    inline const RPPSimHitDebugInfo &GetOriginalProtonAtStation(RPStationId st_id) {return _primary_proton_enters_station[st_id][0];}
    inline bool PrimaryProtonOutgoingFromStation(RPStationId st_id) {return _primary_proton_leaves_station[st_id].size()>0;}
    inline const RPPSimHitDebugInfo &GetOutgoingPrimaryProtonFromStation(RPStationId st_id) {return _primary_proton_leaves_station[st_id][0];}
    inline int NumberofPrimaryDesintegratedProtons(RPStationId st_id) {return _primary_proton_desintegrated_in_station[st_id].size();}
    inline const std::vector<RPPSimHitDebugInfo> & GetPrimaryDesintegratedProtons(RPStationId st_id) {return _primary_proton_desintegrated_in_station[st_id];}
    inline int NumberofRPsSeingPrimaryProton(RPStationId st_id) {return _number_of_RPs_on_the_way_of_prim_proton[st_id];}
        
//    inline int GetOutgoingParticleMultiplicityFrom220RightStation() const;
//    inline int GetOutgoingProtonMultiplicityFrom220RightStation() const;
    
    inline int GetOutgoingParticleMultiplicityFromStation(RPStationId st_id) {return _particle_leaves_station[st_id].size();}
    //inline int GetOutgoingProtonMultiplicityFromStation(RPStationId st_id) {return _primary_proton_leaves_station[st_id].size();}
    
    inline int GetParticleEnteringRPMultiplicity(RPId rp_id) {return _the_particle_enters_pot[rp_id].size();}
    inline int GetParticleLeavingRPMultiplicity(RPId rp_id) {return _the_particle_leaves_pot[rp_id].size();}
    inline int GetPrimProtonEnteringRPMultiplicity(RPId rp_id) {return _the_primary_proton_enters_pot[rp_id].size();}
    inline int GetPrimProtonLeavingRPMultiplicity(RPId rp_id) {return _the_primary_proton_leaves_pot[rp_id].size();}
    inline int GetPrimProtonInelasticRPMultiplicity(RPId rp_id) {return _primary_proton_desintegrated_in_the_pot[rp_id].size();}
    inline int GetParticleLeavingRPFrontWallMult(RPId rp_id) {return _particle_leaves_the_front_wall_of_the_pot[rp_id].size();}
    
    inline const std::vector<RPPSimHitDebugInfo> & GetParticleEnteringRP(RPId rp_id) {return _the_particle_enters_pot[rp_id];}
    inline const std::vector<RPPSimHitDebugInfo> & GetParticleLeavingRP(RPId rp_id) {return _the_particle_leaves_pot[rp_id];}
    inline const std::vector<RPPSimHitDebugInfo> & GetPrimProtonEnteringRP(RPId rp_id) {return _the_primary_proton_enters_pot[rp_id];}
    inline const std::vector<RPPSimHitDebugInfo> & GetPrimProtonLeavingRP(RPId rp_id) {return _the_primary_proton_leaves_pot[rp_id];}
    inline const std::vector<RPPSimHitDebugInfo> & GetPrimProtonInelasticRP(RPId rp_id) {return _primary_proton_desintegrated_in_the_pot[rp_id];}
    inline const std::vector<RPPSimHitDebugInfo> & GetParticlesLeavingRPFrontWall(RPId rp_id) {return _particle_leaves_the_front_wall_of_the_pot[rp_id];}
    
  private:
    typedef std::vector<RPPSimHitDebugInfo> PSimHitVector;
    typedef std::map<RPId, PSimHitVector> RP_Debug_Info_Vector;
    RP_Debug_Info_Vector _the_primary_proton_enters_pot;
    RP_Debug_Info_Vector _the_particle_enters_pot;
    RP_Debug_Info_Vector _the_primary_proton_leaves_pot;
    RP_Debug_Info_Vector _the_particle_leaves_pot;
    RP_Debug_Info_Vector _primary_proton_desintegrated_in_the_pot;
    RP_Debug_Info_Vector _particle_leaves_the_front_wall_of_the_pot;
    
    typedef std::map<RPStationId, PSimHitVector> Station_Debug_Info_Vector;
    Station_Debug_Info_Vector _primary_proton_enters_station;
    Station_Debug_Info_Vector _particle_leaves_station;
    Station_Debug_Info_Vector _primary_proton_leaves_station;
    Station_Debug_Info_Vector _primary_proton_desintegrated_in_station;
    
    typedef std::map<bool, PSimHitVector> Hemisphere_Debug_Info_Vector;
    Hemisphere_Debug_Info_Vector _primary_protons_at_IP;  //true to the right, false to the left
    
    typedef std::map<RPStationId, int> Count_Map;
    Count_Map _number_of_RPs_on_the_way_of_prim_proton;
    
    double mp, p0, E1;
};


#endif /*RPSUPPLEMENTARYINFO_H_*/
