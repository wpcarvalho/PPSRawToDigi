#include "TotemRPValidation/RPGeant4Validation/interface/RPSupplementaryInfo.h"

#include "TMath.h"

RPSupplementaryInfo::RPSupplementaryInfo()
{
  mp = 0.938272029;
  p0 = 7e3;
  E1 = TMath::Sqrt(p0*p0 + mp*mp);
}

void RPSupplementaryInfo::SetCurrentEvent(const std::vector<RPPSimHitDebugInfo> &in_sim_vect)
{
//  std::cout<<"in_sim_vect size="<<in_sim_vect.size()<<std::endl;
  _the_primary_proton_enters_pot.clear();
  _the_particle_enters_pot.clear();
  _the_primary_proton_leaves_pot.clear();
  _the_particle_leaves_pot.clear();
  _primary_proton_desintegrated_in_the_pot.clear();
  
  _primary_proton_enters_station.clear();
  _particle_leaves_station.clear();
  _primary_proton_leaves_station.clear();
  _primary_proton_desintegrated_in_station.clear();
  _particle_leaves_the_front_wall_of_the_pot.clear();
    
  _primary_protons_at_IP.clear();  //true to the right, false to the left
  
  _number_of_RPs_on_the_way_of_prim_proton.clear();
  
  int size = in_sim_vect.size();
  for(int i=0; i<size; ++i)
  {
    //std::cout<<"in_sim_vect[i].detUnitId()="<<(int)in_sim_vect[i].detUnitId()<<std::endl;
    
    if((int)in_sim_vect[i].detUnitId()>=0)
      continue;
      
    if((int)in_sim_vect[i].detUnitId()==-1)
    {
      //get primary protons at diffrent locations
      //std::cout<<"primary proton added="<<std::endl;
      if(in_sim_vect[i].particleType()==2212 && in_sim_vect[i].GetParentId()==0)
      {
        //in IP
        double z_vertex_pos = in_sim_vect[i].GetPrimaryVertex().z();
        if( fabs(z_vertex_pos)<1000 )
        {
          //if to the right; true to the right, false to the left
          _primary_protons_at_IP[ tan(in_sim_vect[i].thetaAtEntry())>0 ].push_back(in_sim_vect[i]);
        }
//        else if( RPDetSpaceGeometry::IsInsideStation(z_vertex_pos) )
//        {
//          RPStationId st_id = RPDetSpaceGeometry::GetStationIdByVertex(z_vertex_pos);
//          _primary_proton_enters_station[st_id].push_back(in_sim_vect[i]);
//        }
      }
    }
    //particle leaves 220 right
    else if((int)in_sim_vect[i].detUnitId()==-2)
    {
      RPStationId st_id = RPDetSpaceGeometry::GetStationIdBySeqNo(3);
      _particle_leaves_station[st_id].push_back(in_sim_vect[i]);
      if(in_sim_vect[i].particleType()==2212 && in_sim_vect[i].GetParentId()==0)
        _primary_proton_leaves_station[st_id].push_back(in_sim_vect[i]);
    }
    //particle enters rp
    else if((int)in_sim_vect[i].detUnitId()==-3)
    {
//      cout<<"particle enters RP"<<endl;
      if(in_sim_vect[i].particleType()==2212 && in_sim_vect[i].GetParentId()==0)
      {
        _the_primary_proton_enters_pot[in_sim_vect[i].GetRPId()].push_back(in_sim_vect[i]);
        RPStationId st_id = RPDetSpaceGeometry::GetStationId(in_sim_vect[i].GetRPId());
        _number_of_RPs_on_the_way_of_prim_proton[st_id]++;
      }
      _the_particle_enters_pot[in_sim_vect[i].GetRPId()].push_back(in_sim_vect[i]);
    }
    //particle leaves rp
    else if((int)in_sim_vect[i].detUnitId()==-4)
    {
//      cout<<"particle leaves RP"<<endl;
      if(in_sim_vect[i].particleType()==2212 && in_sim_vect[i].GetParentId()==0)
      {
        _the_primary_proton_leaves_pot[in_sim_vect[i].GetRPId()].push_back(in_sim_vect[i]);
      }
      _the_particle_leaves_pot[in_sim_vect[i].GetRPId()].push_back(in_sim_vect[i]);
    }
    //particle leaves front_wall of rp
    else if((int)in_sim_vect[i].detUnitId()==-9)
    {
      //std::cout<<"particle leaves front_wall of rp:"<<in_sim_vect[i].GetRPId()<<std::endl;
      _particle_leaves_the_front_wall_of_the_pot[in_sim_vect[i].GetRPId()].push_back(in_sim_vect[i]);
    }
    //particle stopped in rp
    else if((int)in_sim_vect[i].detUnitId()==-5)
    {
//      cout<<"primary proton lost in RP:"<<in_sim_vect[i].GetRPId()
//          <<" part:"<<in_sim_vect[i].GetRPDetPartId()
//          <<" copy no:"<<in_sim_vect[i].GetCopyNo()<<endl;
      _primary_proton_desintegrated_in_the_pot[in_sim_vect[i].GetRPId()].push_back(in_sim_vect[i]);
      _primary_proton_desintegrated_in_station[RPDetSpaceGeometry::GetStationId(in_sim_vect[i].GetRPId())].push_back(in_sim_vect[i]);
    }
    //particle leaves 220 left
    else if((int)in_sim_vect[i].detUnitId()==-6)
    {
      RPStationId st_id = RPDetSpaceGeometry::GetStationIdBySeqNo(-3);
      _particle_leaves_station[st_id].push_back(in_sim_vect[i]);
      if(in_sim_vect[i].particleType()==2212 && in_sim_vect[i].GetParentId()==0)
        _primary_proton_leaves_station[st_id].push_back(in_sim_vect[i]);
    }
    //particle leaves 147 right
    else if((int)in_sim_vect[i].detUnitId()==-7)
    {
      RPStationId st_id = RPDetSpaceGeometry::GetStationIdBySeqNo(1);
      _particle_leaves_station[st_id].push_back(in_sim_vect[i]);
      if(in_sim_vect[i].particleType()==2212 && in_sim_vect[i].GetParentId()==0)
        _primary_proton_leaves_station[st_id].push_back(in_sim_vect[i]);
    }
    //particle leaves 147 left
    else if((int)in_sim_vect[i].detUnitId()==-8)
    {
      RPStationId st_id = RPDetSpaceGeometry::GetStationIdBySeqNo(-1);
      _particle_leaves_station[st_id].push_back(in_sim_vect[i]);
      if(in_sim_vect[i].particleType()==2212 && in_sim_vect[i].GetParentId()==0)
        _primary_proton_leaves_station[st_id].push_back(in_sim_vect[i]);
    }
    //particle enters station
    else if((int)in_sim_vect[i].detUnitId()==-10)
    {
      double z_vertex_pos = in_sim_vect[i].GetPrimaryVertex().z();
      if( RPDetSpaceGeometry::IsInsideStation(z_vertex_pos) )
      {
        RPStationId st_id = RPDetSpaceGeometry::GetStationIdByVertex(z_vertex_pos);
        _primary_proton_enters_station[st_id].push_back(in_sim_vect[i]);
      }
    }
  }
}

double RPSupplementaryInfo::Get_Log10t_OfOriginalProtonTowardsStation(RPStationId st_id)
{
  const RPPSimHitDebugInfo &proton = GetOriginalProtonTowardsStation(st_id);
  LocalVector p_0 = proton.momentumAtEntry();
  double PABS = proton.pabs();
  double E3 = TMath::Sqrt(PABS*PABS+mp*mp);
  double delta_p_z = p0 - p_0.z();
  double t = (E1-E3)*(E1-E3) - p_0.x()*p_0.x() - p_0.y()*p_0.y() - delta_p_z*delta_p_z;
  double log10_t = TMath::Log(-t)/TMath::Log(10.0);
  return log10_t;
}

