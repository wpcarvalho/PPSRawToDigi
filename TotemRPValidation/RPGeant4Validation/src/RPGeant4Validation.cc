#include "TotemRPValidation/RPGeant4Validation/interface/RPGeant4Validation.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPSupplementaryInfo.h"
#include "TotemRPValidation/RPGeant4Validation/interface/StationValidation.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPDetSpaceGeometry.h"
#include "TFile.h"
#include <iostream>


RPGeant4Validation::RPGeant4Validation()
 : rp_coll("/RP_Hists/RP_", edm::ParameterSet()), station_coll("/Station_Hists/Station_", edm::ParameterSet())
{
  rp_ids = det_sp_geom.GetAllRPRomanPotsIdList();
  station_ids = det_sp_geom.GetRPStationIdList();
}


void RPGeant4Validation::FillHistogramms(const std::vector<PSimHit>& hits, const std::vector<RPPSimHitDebugInfo>& debug_hits)
{
  RPSupplementaryInfo debugging_info;
  debugging_info.SetCurrentEvent(debug_hits);
  
//  std::cout<<"Fill rp_coll..."<<std::endl;
  for(unsigned int i=0; i<rp_ids.size(); ++i)
  {
    rp_coll.GetObj(rp_ids[i])->FillHistogramms(debugging_info);
  }
  
//  std::cout<<"Fill station_coll..."<<std::endl;
  for(unsigned int i=0; i<station_ids.size(); ++i)
  {
    station_coll.GetObj(station_ids[i])->FillHistogramms(debugging_info);
  }
  
//  std::cout<<"Fill rp_coll, simhits..."<<std::endl;
  for(unsigned int i=0; i<hits.size(); ++i)
  {
    //TotRPDetId id(hits[i].detUnitId());
    rp_coll.GetObj(det_sp_geom.GetRPId(hits[i].detUnitId()))->FillHistogramms(hits[i], debugging_info);
  }
//  std::cout<<"Finnished."<<std::endl;
}


void RPGeant4Validation::WriteHistograms(const std::string &root_file_name)
{
  TFile *f = TFile::Open(root_file_name.c_str(), "recreate");
  rp_coll.Write(f);
  station_coll.Write(f);
  f->Close();
}
