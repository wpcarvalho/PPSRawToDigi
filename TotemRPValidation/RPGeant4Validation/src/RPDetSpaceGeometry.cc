#include "TotemRPValidation/RPGeant4Validation/interface/RPDetSpaceGeometry.h"
#include "TMath.h"


RPDetSpaceGeometry::RPDetSpaceGeometry()
{
}


std::vector<RPDetId> RPDetSpaceGeometry::GetRPDetIdList()
{
  /* tempemorarily instantiate det. digitizers explicitly not basing on the geometry
   * information
   * xyzk, x-side, y-station, z-RP, k-det
   * x: 0.1
   * y: 0..2
   * z: 0..6
   * k: 0..9
   */
   
  std::vector<RPDetId> det_ids;
  for(int x=0; x<=1; x++)
    for(int y=0; y<=2; y++)
      for(int z=0; z<=5; z++)
        for(int k=0; k<=9; k++)
        {
          RPDetId det_id = 1000*x+100*y+10*z+k; 
          det_ids.push_back(det_id);
        }
        
  return det_ids;
}


std::vector<RPId> RPDetSpaceGeometry::GetAllRPRomanPotsIdList()
{
  /* tempemorarily instantiate det. digitizers explicitly not basing on the geometry
   * information
   * xyzk, x-side, y-station, z-RP
   * x: 0.1
   * y: 0..2
   * z: 0..6
   */
   
  std::vector<RPId> rp_ids;
  for(int x=0; x<=1; x++)
    for(int y=0; y<=2; y++)
      for(int z=0; z<=5; z++)
      {
        RPId rp_id = 1000*x+100*y+10*z; 
        rp_ids.push_back(rp_id);
      }
        
  return rp_ids;
}


std::vector<RPStationId> RPDetSpaceGeometry::GetRPStationIdList()
{
  /* tempemorarily instantiate det. digitizers explicitly not basing on the geometry
   * information
   * xyzk, x-side, y-station, z-RP, k-det
   * x: 0.1
   * y: 0..2
   */
  
  std::vector<RPStationId> station_ids;
  
  for(int x=0; x<=1; x++)
    for(int y=0; y<=2; y++)
    {
      RPStationId station_id = 1000*x+100*y;
      station_ids.push_back(station_id);
    }
  return station_ids;
}


std::vector<RPId> RPDetSpaceGeometry::GetRPRomanPotsIdList(RPStationId station_id)
{
  /* tempemorarily instantiate det. digitizers explicitly not basing on the geometry
   * information
   * xyzk, x-side, y-station, z-RP, k-det
   * x: 0.1
   * y: 0..2
   * z: 0..6
   * k: 0..9
   */
   
  std::vector<RPId> roman_pot_ids;
  for(int z=0; z<=5; z++)
  {
    RPId rp_id = station_id+10*z;
    roman_pot_ids.push_back(rp_id);
  }
        
  return roman_pot_ids;
}

std::vector<RPDetId> RPDetSpaceGeometry::GetRPDetIdList(RPId rp_id)
{
  /* tempemorarily instantiate det. digitizers explicitly not basing on the geometry
   * information
   * xyzk, x-side, y-station, z-RP, k-det
   * x: 0.1
   * y: 0..2
   * z: 0..6
   * k: 0..9
   */
   
  std::vector<RPDetId> det_ids;
  for(int k=0; k<=9; k++)
  {
    RPDetId det_id = rp_id+k; 
    det_ids.push_back(det_id);
  }
        
  return det_ids;
}

bool RPDetSpaceGeometry::IsInsideStation(double z)
{
  double abs_z = TMath::Abs(z);
  return (130000<abs_z && abs_z<160000) || (210000<abs_z && abs_z<230000);
}

RPStationId RPDetSpaceGeometry::GetStationIdByVertex(double z)
{
  if(210000<z && z<230000)
    return 1200;
  if(130000<z && z<160000)
    return 1000;
  if(-230000<z && z<-210000)
    return 200;
  if(-160000<z && z<-130000)
    return 0;
  return 0;
}


RPStationId RPDetSpaceGeometry::GetStationIdBySeqNo(int seq_no)
{
  switch(seq_no)
  {
    case 1: return 1000;
    case 2: return 1100;
    case 3: return 1200;
    case -1: return 0;
    case -2: return 100;
    case -3: return 200;
  };
  return 0;
}
