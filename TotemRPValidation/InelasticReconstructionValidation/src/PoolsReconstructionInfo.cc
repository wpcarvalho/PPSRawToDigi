#include "TotemRPValidation/InelasticReconstructionValidation/interface/PoolsReconstructionInfo.h"
#include "TH1F.h"
#include "TVector2.h"


PoolsReconstructionInfo::PoolsReconstructionInfo(const std::string &path, RPId rp_id, const edm::ParameterSet&)
 : BaseHistogramManager(path), rp_id_(rp_id)
{
  InitializeHistograms();
}


void PoolsReconstructionInfo::InitializeHistograms()
{
  char name[1024];
  sprintf(name, "simulated_vs_theoretical_residuum_x_%04i", rp_id_);
  simulated_vs_theoretical_residuum_x_ = TH1F(name, name, 200, -0.001, 0.001);
  simulated_vs_theoretical_residuum_x_.SetBit(TH1::kCanRebin);
  simulated_vs_theoretical_residuum_x_.SetDirectory(0);
  simulated_vs_theoretical_residuum_x_.SetXTitle("#Delta x [mm]");
  simulated_vs_theoretical_residuum_x_.SetYTitle("Entries");
  RegisterHistogram(simulated_vs_theoretical_residuum_x_);
  
  sprintf(name, "simulated_vs_reconstructed_residuum_x_%04i", rp_id_);
  simulated_vs_reconstructed_residuum_x_ = TH1F(name, name, 200, -0.001, 0.001);
  simulated_vs_reconstructed_residuum_x_.SetBit(TH1::kCanRebin);
  simulated_vs_reconstructed_residuum_x_.SetDirectory(0);
  simulated_vs_reconstructed_residuum_x_.SetXTitle("#Delta x [mm]");
  simulated_vs_reconstructed_residuum_x_.SetYTitle("Entries");
  RegisterHistogram(simulated_vs_reconstructed_residuum_x_);
  
  sprintf(name, "reconstructed_vs_theoretical_residuum_x_%04i", rp_id_);
  reconstructed_vs_theoretical_residuum_x_ = TH1F(name, name, 200, -0.001, 0.001);
  reconstructed_vs_theoretical_residuum_x_.SetBit(TH1::kCanRebin);
  reconstructed_vs_theoretical_residuum_x_.SetDirectory(0);
  reconstructed_vs_theoretical_residuum_x_.SetXTitle("#Delta x [mm]");
  reconstructed_vs_theoretical_residuum_x_.SetYTitle("Entries");
  RegisterHistogram(reconstructed_vs_theoretical_residuum_x_);
  
  sprintf(name, "simulated_vs_theoretical_residuum_y_%04i", rp_id_);
  simulated_vs_theoretical_residuum_y_ = TH1F(name, name, 200, -0.001, 0.001);
  simulated_vs_theoretical_residuum_y_.SetBit(TH1::kCanRebin);
  simulated_vs_theoretical_residuum_y_.SetDirectory(0);
  simulated_vs_theoretical_residuum_y_.SetXTitle("#Delta y [mm]");
  simulated_vs_theoretical_residuum_y_.SetYTitle("Entries");
  RegisterHistogram(simulated_vs_theoretical_residuum_y_);
  
  sprintf(name, "simulated_vs_reconstructed_residuum_y_%04i", rp_id_);
  simulated_vs_reconstructed_residuum_y_ = TH1F(name, name, 200, -0.001, 0.001);
  simulated_vs_reconstructed_residuum_y_.SetBit(TH1::kCanRebin);
  simulated_vs_reconstructed_residuum_y_.SetDirectory(0);
  simulated_vs_reconstructed_residuum_y_.SetXTitle("#Delta y [mm]");
  simulated_vs_reconstructed_residuum_y_.SetYTitle("Entries");
  RegisterHistogram(simulated_vs_reconstructed_residuum_y_);
  
  sprintf(name, "reconstructed_vs_theoretical_residuum_y_%04i", rp_id_);
  reconstructed_vs_theoretical_residuum_y_ = TH1F(name, name, 200, -0.001, 0.001);
  reconstructed_vs_theoretical_residuum_y_.SetBit(TH1::kCanRebin);
  reconstructed_vs_theoretical_residuum_y_.SetDirectory(0);
  reconstructed_vs_theoretical_residuum_y_.SetXTitle("#Delta y [mm]");
  reconstructed_vs_theoretical_residuum_y_.SetYTitle("Entries");
  RegisterHistogram(reconstructed_vs_theoretical_residuum_y_);
  
//  sprintf(name, "reconstructed_vs_theoretical_pool_x_%04i", rp_id_);
//  reconstructed_vs_theoretical_pool_x_ = TH1F(name, name, 2000, -20, 20);
//  //position_residuum_y_.SetBit(TH1::kCanRebin);
//  reconstructed_vs_theoretical_pool_x_.SetDirectory(0);
//  RegisterHistogram(reconstructed_vs_theoretical_pool_x_);
  
//  sprintf(name, "reconstructed_vs_theoretical_pool_y_%04i", rp_id_);
//  reconstructed_vs_theoretical_pool_y_ = TH1F(name, name, 2000, -20, 20);
//  //reconstructed_vs_theoretical_residuum_y_.SetBit(TH1::kCanRebin);
//  reconstructed_vs_theoretical_pool_y_.SetDirectory(0);
//  RegisterHistogram(reconstructed_vs_theoretical_pool_y_);
}

//const TVector2 &theoretical_hit (the generated proton transported without smearing to the RP location) 
//const RP2DHit &simulated_hit (the reconstructed proton transported without smearing to the RP location)
//const TVector2 &reconstructed_hit (the reocnstructed hit in the PRs)

void PoolsReconstructionInfo::Fill(const TVector2 &theoretical_hit, const RP2DHit &simulated_hit, 
    const TVector2 &reconstructed_hit)
{
  double res_sim_theor_x = simulated_hit.X() - theoretical_hit.X();
  double res_sim_theor_y = simulated_hit.Y() - theoretical_hit.Y();
  double res_sim_rec_x = simulated_hit.X() - reconstructed_hit.X();
  double res_sim_rec_y = simulated_hit.Y() - reconstructed_hit.Y();
  double res_rec_theor_x = reconstructed_hit.X() - theoretical_hit.X();
  double res_rec_theor_y = reconstructed_hit.Y() - theoretical_hit.Y();
//  double pool_x = res_rec_theor_x/simulated_hit.Sx();
//  double pool_y = res_rec_theor_y/simulated_hit.Sy();
  
  simulated_vs_theoretical_residuum_x_.Fill(res_sim_theor_x);
  simulated_vs_reconstructed_residuum_x_.Fill(res_sim_rec_x);
  reconstructed_vs_theoretical_residuum_x_.Fill(res_rec_theor_x);
  
  simulated_vs_theoretical_residuum_y_.Fill(res_sim_theor_y);
  simulated_vs_reconstructed_residuum_y_.Fill(res_sim_rec_y);
  reconstructed_vs_theoretical_residuum_y_.Fill(res_rec_theor_y);
  
//  reconstructed_vs_theoretical_pool_x_.Fill(pool_x);
//  reconstructed_vs_theoretical_pool_y_.Fill(pool_y);
}

