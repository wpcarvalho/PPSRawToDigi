#include "TotemRPValidation/RPGeant4Validation/interface/StationValidation.h"
#include "TotemRPValidation/RPGeant4Validation/interface/RPSupplementaryInfo.h"
#include "Minuit2/MnMigrad.h"
#include "TMath.h"
#include <iostream>

StationValidation::StationValidation(const std::string &path, RPStationId st_id, const edm::ParameterSet& conf)
 : BaseHistogramManager(path), _station_id(st_id)
{
  Init();
}


void StationValidation::FillHistogramms(RPSupplementaryInfo &sup_info)
{
//  std::cout<<"sup_info.OriginalProtonAtStationGot(_station_id)="<<
//      sup_info.OriginalProtonAtStationGot(_station_id)<<
//      " sup_info.OriginalProtonTowardsStationGot(_station_id)="<<
//      sup_info.OriginalProtonTowardsStationGot(_station_id)<<" _station_id="<<_station_id<<std::endl;
  
  if(sup_info.OriginalProtonAtStationGot(_station_id) && sup_info.OriginalProtonTowardsStationGot(_station_id))
  {
    const RPPSimHitDebugInfo & prot_st = sup_info.GetOriginalProtonAtStation(_station_id);
    const RPPSimHitDebugInfo & prot_ip = sup_info.GetOriginalProtonTowardsStation(_station_id);
    double mom = prot_st.pabs();
    double pos_x_st = prot_st.GetGlobalPosition().x();
    double pos_y_st = prot_st.GetGlobalPosition().y();
    double pos_x_ip = prot_ip.GetGlobalPosition().x();
    double pos_y_ip = prot_ip.GetGlobalPosition().y();
    double ksi = (mom-7e3)/7e3;
    double theta_x_st = TMath::ATan2(prot_st.momentumAtEntry().x(), TMath::Abs(prot_st.momentumAtEntry().z()));
    double theta_y_st = TMath::ATan2(prot_st.momentumAtEntry().y(), TMath::Abs(prot_st.momentumAtEntry().z()));
    double theta_x_ip = TMath::ATan2(prot_ip.momentumAtEntry().x(), TMath::Abs(prot_ip.momentumAtEntry().z()));
    double theta_y_ip = TMath::ATan2(prot_ip.momentumAtEntry().y(), TMath::Abs(prot_ip.momentumAtEntry().z()));
    
//    std::cout<<theta_x_st<<" "<<theta_y_st<<" "<<theta_x_ip<<" "<<theta_y_ip<<std::endl;
//    
    _prim_prot_at_station_mom_dist.Fill(mom);
    _prim_prot_at_station_theta_dist_x.Fill(theta_x_st);
    _prim_prot_at_station_theta_dist_y.Fill(theta_y_st);
    _prim_prot_at_station_pos_dist_x.Fill(pos_x_st);
    _prim_prot_at_station_pos_dist_y.Fill(pos_y_st);
    _prim_prot_at_station_pos_dist_xy.Fill(pos_x_st, pos_y_st);
    _prim_prot_pos_at_station_vs_theta_at_0_x.Fill(theta_x_ip, pos_x_st);
    _prim_prot_pos_at_station_vs_theta_at_0_y.Fill(theta_y_ip, pos_y_st);
    _prim_prot_pos_at_station_vs_pos_at_0_x.Fill(pos_x_ip, pos_x_st);
    _prim_prot_pos_at_station_vs_pos_at_0_y.Fill(pos_y_ip, pos_y_st);
    _prim_prot_pos_at_station_vs_ksi_x.Fill(ksi, pos_x_st);
    _prim_prot_pos_at_station_vs_ksi_y.Fill(ksi, pos_y_st);
    _prim_prot_theta_at_station_vs_theta_at_0_x.Fill(theta_x_ip, theta_x_st);
    _prim_prot_theta_at_station_vs_theta_at_0_y.Fill(theta_y_ip, theta_y_st);
  }
}


void StationValidation::Init()
{
  char name[1024];
  
  sprintf(name, "prim_prot_at_station_mom_dist_%04i", _station_id);
  _prim_prot_at_station_mom_dist = TH1F(name, name, 601, 6.9999e3, 7.00001e3);
  _prim_prot_at_station_mom_dist.SetDirectory(0);
  _prim_prot_at_station_mom_dist.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_at_station_mom_dist);
  
  sprintf(name, "prim_prot_at_station_theta_dist_x_%04i", _station_id);
  _prim_prot_at_station_theta_dist_x = TH1F(name, name, 601, -1e-6, 1e-6);
  _prim_prot_at_station_theta_dist_x.SetDirectory(0);
  _prim_prot_at_station_theta_dist_x.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_at_station_theta_dist_x);
  
  sprintf(name, "prim_prot_at_station_theta_dist_y_%04i", _station_id);
  _prim_prot_at_station_theta_dist_y = TH1F(name, name, 601, -1e-6, 1e-6);
  _prim_prot_at_station_theta_dist_y.SetDirectory(0);
  _prim_prot_at_station_theta_dist_y.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_at_station_theta_dist_y);
  
  sprintf(name, "prim_prot_at_station_pos_dist_x_%04i", _station_id);
  _prim_prot_at_station_pos_dist_x = TH1F(name, name, 601, -1e-6, 1e-6);
  _prim_prot_at_station_pos_dist_x.SetDirectory(0);
  _prim_prot_at_station_pos_dist_x.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_at_station_pos_dist_x);
  
  sprintf(name, "prim_prot_at_station_pos_dist_y_%04i", _station_id);
  _prim_prot_at_station_pos_dist_y = TH1F(name, name, 601, -1e-6, 1e-6);
  _prim_prot_at_station_pos_dist_y.SetDirectory(0);
  _prim_prot_at_station_pos_dist_y.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_at_station_pos_dist_y);
  
  sprintf(name, "prim_prot_at_station_pos_dist_xy_%04i", _station_id);
  _prim_prot_at_station_pos_dist_xy = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_at_station_pos_dist_xy.SetDirectory(0);
  _prim_prot_at_station_pos_dist_xy.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_at_station_pos_dist_xy);
  
  sprintf(name, "prim_prot_pos_at_station_vs_theta_at_0_x_%04i", _station_id);
  _prim_prot_pos_at_station_vs_theta_at_0_x = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_pos_at_station_vs_theta_at_0_x.SetDirectory(0);
  _prim_prot_pos_at_station_vs_theta_at_0_x.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_pos_at_station_vs_theta_at_0_x);
  
  sprintf(name, "prim_prot_pos_at_station_vs_theta_at_0_y_%04i", _station_id);
  _prim_prot_pos_at_station_vs_theta_at_0_y = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_pos_at_station_vs_theta_at_0_y.SetDirectory(0);
  _prim_prot_pos_at_station_vs_theta_at_0_y.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_pos_at_station_vs_theta_at_0_y);
  
  sprintf(name, "prim_prot_pos_at_station_vs_pos_at_0_x_%04i", _station_id);
  _prim_prot_pos_at_station_vs_pos_at_0_x = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_pos_at_station_vs_pos_at_0_x.SetDirectory(0);
  _prim_prot_pos_at_station_vs_pos_at_0_x.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_pos_at_station_vs_pos_at_0_x);
  
  sprintf(name, "prim_prot_pos_at_station_vs_pos_at_0_y_%04i", _station_id);
  _prim_prot_pos_at_station_vs_pos_at_0_y = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_pos_at_station_vs_pos_at_0_y.SetDirectory(0);
  _prim_prot_pos_at_station_vs_pos_at_0_y.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_pos_at_station_vs_pos_at_0_y);
  
  sprintf(name, "prim_prot_pos_at_station_vs_ksi_x_%04i", _station_id);
  _prim_prot_pos_at_station_vs_ksi_x = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_pos_at_station_vs_ksi_x.SetDirectory(0);
  _prim_prot_pos_at_station_vs_ksi_x.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_pos_at_station_vs_ksi_x);
  
  sprintf(name, "prim_prot_pos_at_station_vs_ksi_y_%04i", _station_id);
  _prim_prot_pos_at_station_vs_ksi_y = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_pos_at_station_vs_ksi_y.SetDirectory(0);
  _prim_prot_pos_at_station_vs_ksi_y.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_pos_at_station_vs_ksi_y);
  
  sprintf(name, "prim_prot_theta_at_station_vs_theta_at_0_x_%04i", _station_id);
  _prim_prot_theta_at_station_vs_theta_at_0_x = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_theta_at_station_vs_theta_at_0_x.SetDirectory(0);
  _prim_prot_theta_at_station_vs_theta_at_0_x.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_theta_at_station_vs_theta_at_0_x);
  
  sprintf(name, "prim_prot_theta_at_station_vs_theta_at_0_y_%04i", _station_id);
  _prim_prot_theta_at_station_vs_theta_at_0_y = TH2F(name, name, 601, -1e-6, 1e-6, 601, -1e-6, 1e-6);
  _prim_prot_theta_at_station_vs_theta_at_0_y.SetDirectory(0);
  _prim_prot_theta_at_station_vs_theta_at_0_y.SetBit(TH1::kCanRebin);
  RegisterHistogram(_prim_prot_theta_at_station_vs_theta_at_0_y);
}


void StationValidation::Finalize()
{
}
