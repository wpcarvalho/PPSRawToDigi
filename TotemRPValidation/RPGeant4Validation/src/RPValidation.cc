#include "TotemRPValidation/RPGeant4Validation/interface/RPValidation.h"

RPValidation::RPValidation(const std::string &path, RPId rp_id, const edm::ParameterSet&)
 : BaseHistogramManager(path), _rp_id(rp_id)
{
  InitializeHistograms();
}


void RPValidation::FillHistogramms(RPSupplementaryInfo &sup_info)
{
  if(sup_info.GetPrimProtonInelasticRPMultiplicity(_rp_id)==1)
  {
    //interacting parts of the RP
    //std::cout<<"interacting parts of the RP"<<std::endl;
    _pot_interaction_parts.Fill(sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetRPDetPartId());

//    std::cout<<"sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetRPDetPartId()="<<sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetRPDetPartId()<<
  //    " sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetCopyNo()="<<sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetCopyNo()<<std::endl;            
    //fill histograms for the interaction in the front wall
    if(sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetRPDetPartId()==6 && sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetCopyNo()==0)
    {
      int part_num = sup_info.GetParticleLeavingRPFrontWallMult(_rp_id);
      const std::vector<RPPSimHitDebugInfo> &particles_leaving_fr_wall = sup_info.GetParticlesLeavingRPFrontWall(_rp_id);
    //  std::cout<<"Interaction in the front wall of "<<_rp_id<<" RP"<<std::endl;
      //std::cout<<"all particles: "<<std::endl;
      int track_no_above_1_GeV = 0;
      for(int i=0; i<part_num; ++i)
      {
//        cout<<i<<" "<<RPHepPDTWrapper::GetName(particles_leaving_fr_wall[i].particleType())
////              <<", dir:"<<particles_leaving_fr_wall[i].localDirection()
////              <<", pos:"<<particles_leaving_fr_wall[i].GetGlobalPosition()
//            <<", Theta:"<<particles_leaving_fr_wall[i].thetaAtEntry()*1000.0<<" mrad, Pabs:"
//            <<particles_leaving_fr_wall[i].pabs()*1000.0<<" MeV"<<endl;
        
        _front_wall_inter_PDG_dist.Fill(particles_leaving_fr_wall[i].particleType());
        _front_wall_momentum_theta_dist.Fill(particles_leaving_fr_wall[i].pabs(), particles_leaving_fr_wall[i].thetaAtEntry());
        if(RPHepPDTWrapper::GetCharge(particles_leaving_fr_wall[i].particleType())!=0.0)
        {
          _front_wall_momentum_theta_charged_particle_dist.Fill(particles_leaving_fr_wall[i].pabs(), particles_leaving_fr_wall[i].thetaAtEntry());
        }
        if(particles_leaving_fr_wall[i].pabs()>1.0)
        {
          ++track_no_above_1_GeV;
          _front_wall_inter_PDG_dist_charged_above_1_GeV.Fill(particles_leaving_fr_wall[i].particleType());
          _front_wall_inter_theta_dist_charged_above_1_GeV.Fill(particles_leaving_fr_wall[i].thetaAtEntry());
        }
      }
      _front_wall_inter_track_no_dist.Fill(part_num);
      _front_wall_inter_track_no_dist_charged_above_1_GeV.Fill(track_no_above_1_GeV);
    }
  }
  else if(sup_info.GetPrimProtonEnteringRPMultiplicity(_rp_id)==1 && sup_info.GetPrimProtonLeavingRPMultiplicity(_rp_id)==1)
  {
    _pot_interaction_parts.Fill(-1);
  }
  
  //multiple scattering
  if(sup_info.GetPrimProtonEnteringRPMultiplicity(_rp_id)==1 && sup_info.GetPrimProtonLeavingRPMultiplicity(_rp_id)==1)
  {
    LocalVector scat_vect = sup_info.GetPrimProtonLeavingRP(_rp_id)[0].momentumAtEntry().unit() - 
          sup_info.GetPrimProtonEnteringRP(_rp_id)[0].momentumAtEntry().unit();
    
    _pot_scat_dist_x.Fill(scat_vect.x());
    _pot_scat_dist_y.Fill(scat_vect.y());
  }
  
  //angle and position of entry at RP
  if(sup_info.GetPrimProtonEnteringRPMultiplicity(_rp_id)==1)
  {
    LocalVector prot_enter_direction = sup_info.GetPrimProtonEnteringRP(_rp_id)[0].momentumAtEntry().unit();
    _angle_at_RP_x.Fill(prot_enter_direction.x());
    _angle_at_RP_y.Fill(prot_enter_direction.y());
    
    _position_at_RP_xy.Fill(sup_info.GetPrimProtonEnteringRP(_rp_id)[0].GetGlobalPosition().x(), 
        sup_info.GetPrimProtonEnteringRP(_rp_id)[0].GetGlobalPosition().y());
  }
  
  //fraction of inelastic interactions by RP part 
  if( sup_info.GetPrimProtonInelasticRPMultiplicity(_rp_id)>0 )
  {
    int part_id = sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetRPDetPartId();
    int copy_no = sup_info.GetPrimProtonInelasticRP(_rp_id)[0].GetCopyNo();
    if(part_id==1) //silicon
    {
      _interactions_per_rp_part.Fill(copy_no);
    }
    else if(part_id==6) //wall
    {
      _interactions_per_rp_part.Fill(copy_no+10);
    }
    else if(part_id==13) //bottom foil
    {
      _interactions_per_rp_part.Fill(12);
    }
  }
  
  
  if(sup_info.OriginalProtonTowardsStationGot(_rp_id))
  {
    LocalVector p_0 = sup_info.GetOriginalProtonTowardsStation(_rp_id).momentumAtEntry();
    double PABS = sup_info.GetOriginalProtonTowardsStation(_rp_id).pabs();
    double E3 = TMath::Sqrt(PABS*PABS+mp*mp);
    double delta_p_z = p0 - TMath::Abs(p_0.z());
    double t = (E1-E3)*(E1-E3) - p_0.x()*p_0.x() - p_0.y()*p_0.y() - delta_p_z*delta_p_z;
    double log10_t = TMath::Log(-t)/TMath::Log(10.0);
    _log10_t_dist.Fill(log10_t);
    
    if(sup_info.GetPrimProtonInelasticRPMultiplicity(_rp_id)>0)
      _log10_t_dist_inelastic.Fill(log10_t);
    if(sup_info.GetPrimProtonEnteringRPMultiplicity(_rp_id)>0)
      _log10_t_dist_entered_rp.Fill(log10_t);
  }
}


void RPValidation::FillHistogramms(const PSimHit & sim_hit, RPSupplementaryInfo &supp_info)
{
  int mother_part_id = sim_hit.trackId();
  //double momentum = sim_hit.pabs();
  //double theta = sim_hit.thetaAtEntry();
  double en_dep = sim_hit.energyLoss();
  int pdg_part_type = sim_hit.particleType();
  
  if(pdg_part_type==2212 && mother_part_id==0)
  {
    _en_dep_for_primary_protons.Fill(en_dep);
  }
}


void RPValidation::Finalize()
{
  double sc_factor = 1.0/_pot_interaction_parts.GetEntries();
  _pot_interaction_parts.Scale(sc_factor);
  
  RPHepPDTWrapper::SetBinLabels(_front_wall_inter_PDG_dist);
  RPHepPDTWrapper::SetBinLabels(_front_wall_inter_PDG_dist_charged_above_1_GeV);
  _log10_t_dist_inelastic_ratio.Divide(&_log10_t_dist_inelastic, &_log10_t_dist);
  _log10_t_dist_geometrical_acceptance.Divide(&_log10_t_dist_entered_rp, &_log10_t_dist);
}


void RPValidation::InitializeHistograms()
{
//  mp = 0.938272029;
//  p0 = 7e3;
//  E1 = TMath::Sqrt(p0*p0 + mp*mp);
  
  char name[1024];
  
  sprintf(name, "rp_proton_inelastic_interaction_parts_%04i", _rp_id);
  _pot_interaction_parts = TH1F(name, name, 32, -1.5, 30.5);
  _pot_interaction_parts.SetDirectory(0);
  RegisterHistogram(_pot_interaction_parts);
  
  sprintf(name, "rp_proton_mult_scatter_angle_prim_prot_x_%04i", _rp_id);
  _pot_scat_dist_x = TH1F(name, name, 2000, -1e-6, 1e-6);
  _pot_scat_dist_x.SetBit(TH1::kCanRebin);
  _pot_scat_dist_x.SetDirectory(0);
  RegisterHistogram(_pot_scat_dist_x);
  
  sprintf(name, "rp_proton_mult_scatter_angle_prim_prot_y_%04i", _rp_id);
  _pot_scat_dist_y = TH1F(name, name, 2000, -1e-6, 1e-6);
  _pot_scat_dist_y.SetBit(TH1::kCanRebin);
  _pot_scat_dist_y.SetDirectory(0);
  RegisterHistogram(_pot_scat_dist_y);
  
  sprintf(name, "rp_en_dep_for_primary_protons_%04i", _rp_id);
  _en_dep_for_primary_protons = TH1F(name, name, 1000, 0, 0.0015);
  _en_dep_for_primary_protons.SetDirectory(0);
  RegisterHistogram(_en_dep_for_primary_protons);
  
  sprintf(name, "rp_proton_angle_at_RP_x_%04i", _rp_id);
  _angle_at_RP_x = TH1F(name, name, 600, -1e-6, 1e-6);
  _angle_at_RP_x.SetBit(TH1::kCanRebin);
  _angle_at_RP_x.SetDirectory(0);
  RegisterHistogram(_angle_at_RP_x);
  
  sprintf(name, "rp_proton_angle_at_RP_y_%04i", _rp_id);
  _angle_at_RP_y = TH1F(name, name, 600, -1e-6, 1e-6);
  _angle_at_RP_y.SetBit(TH1::kCanRebin);
  _angle_at_RP_y.SetDirectory(0); 
  RegisterHistogram(_angle_at_RP_y);
  
  sprintf(name, "rp_position_at_RP_xy_%04i", _rp_id);
  _position_at_RP_xy = TH2F(name, name, 600, -1e-6, 1e-6, 600, -1e-6, 1e-6);
  _position_at_RP_xy.SetBit(TH1::kCanRebin);
  _position_at_RP_xy.SetDirectory(0);
  RegisterHistogram(_position_at_RP_xy);
  
  sprintf(name, "rp_interactions_per_rp_part_%04i", _rp_id);
  _interactions_per_rp_part = TH1F(name, name, 13, -0.5, 12.5);
  _interactions_per_rp_part.SetDirectory(0);
  RegisterHistogram(_interactions_per_rp_part);
  
  sprintf(name, "rp_front_wall_inter_track_no_dist_%04i", _rp_id);
  _front_wall_inter_track_no_dist = TH1F(name, name, 301, -0.5, 300.5);
  _front_wall_inter_track_no_dist.SetDirectory(0);
  RegisterHistogram(_front_wall_inter_track_no_dist);
  
  sprintf(name, "rp_front_wall_inter_PDG_dist_%04i", _rp_id);
  _front_wall_inter_PDG_dist = TH1F(name, name, 8401, -4200.5, 4200.5);
  _front_wall_inter_PDG_dist.SetDirectory(0);
  RegisterHistogram(_front_wall_inter_PDG_dist);
  
  sprintf(name, "rp_front_wall_momentum_theta_dist_%04i", _rp_id);
  _front_wall_momentum_theta_dist = TH2F(name, name, 600, -1e-10, 1e-10, 600, -1e-10, 1e-10);
  _front_wall_momentum_theta_dist.SetBit(TH1::kCanRebin);
  _front_wall_momentum_theta_dist.SetDirectory(0);
  RegisterHistogram(_front_wall_momentum_theta_dist);
  
  sprintf(name, "rp_front_wall_momentum_theta_charged_particle_dist_%04i", _rp_id);
  _front_wall_momentum_theta_charged_particle_dist = TH2F(name, name, 600, -1e-10, 1e-10, 600, -1e-10, 1e-10);
  _front_wall_momentum_theta_charged_particle_dist.SetBit(TH1::kCanRebin);
  _front_wall_momentum_theta_charged_particle_dist.SetDirectory(0);
  RegisterHistogram(_front_wall_momentum_theta_charged_particle_dist);
  
  sprintf(name, "rp_front_wall_inter_track_no_dist_dist_charged_above_1_GeV_%04i", _rp_id);
  _front_wall_inter_track_no_dist_charged_above_1_GeV = TH1F(name, name, 301, -0.5, 300.5);
  _front_wall_inter_track_no_dist_charged_above_1_GeV.SetDirectory(0);
  RegisterHistogram(_front_wall_inter_track_no_dist_charged_above_1_GeV);
  
  sprintf(name, "rp_front_wall_inter_PDG_dist_charged_above_1_GeV_%04i", _rp_id);
  _front_wall_inter_PDG_dist_charged_above_1_GeV = TH1F(name, name, 8401, -4200.5, 4200.5);
  _front_wall_inter_PDG_dist_charged_above_1_GeV.SetDirectory(0);
  RegisterHistogram(_front_wall_inter_PDG_dist_charged_above_1_GeV);
  
  sprintf(name, "rp_front_wall_inter_theta_dist_charged_above_1_GeV_%04i", _rp_id);
  _front_wall_inter_theta_dist_charged_above_1_GeV = TH1F(name, name, 600, -1e-10, 1e-10);
  _front_wall_inter_theta_dist_charged_above_1_GeV.SetBit(TH1::kCanRebin);
  _front_wall_inter_theta_dist_charged_above_1_GeV.SetDirectory(0);
  RegisterHistogram(_front_wall_inter_theta_dist_charged_above_1_GeV);
  
  sprintf(name, "rp_log10_t_dist_%04i", _rp_id);
  _log10_t_dist = TH1F(name, name, 100, -5, 0);
  _log10_t_dist.SetDirectory(0);
  RegisterHistogram(_log10_t_dist);
  
  sprintf(name, "rp_log10_t_dist_inelastic_%04i", _rp_id);
  _log10_t_dist_inelastic =  TH1F(name, name, 100, -5, 0);
  _log10_t_dist_inelastic.SetDirectory(0);
  RegisterHistogram(_log10_t_dist_inelastic);
  
  sprintf(name, "rp_log10_t_dist_inelastic_ratio_%04i", _rp_id);
  _log10_t_dist_inelastic_ratio =  TH1F(name, name, 100, -5, 0);
  _log10_t_dist_inelastic_ratio.SetDirectory(0);
  RegisterHistogram(_log10_t_dist_inelastic_ratio);
  
  sprintf(name, "rp_log10_t_dist_entered_rp_%04i", _rp_id);
  _log10_t_dist_entered_rp =  TH1F(name, name, 100, -5, 0);
  _log10_t_dist_entered_rp.SetDirectory(0);
  RegisterHistogram(_log10_t_dist_entered_rp);
  
  sprintf(name, "rp_log10_t_dist_geometrical_acceptance_%04i", _rp_id);
  _log10_t_dist_geometrical_acceptance =  TH1F(name, name, 100, -5, 0);
  _log10_t_dist_geometrical_acceptance.SetDirectory(0);
  RegisterHistogram(_log10_t_dist_geometrical_acceptance);
}

