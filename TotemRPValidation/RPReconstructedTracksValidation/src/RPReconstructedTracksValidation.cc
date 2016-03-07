#include "TotemRPValidation/RPReconstructedTracksValidation/interface/RPReconstructedTracksValidation.h"
#include "TFile.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
//#include "FWCore/Framework/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProton.h"
#include <boost/shared_ptr.hpp>
#include "TotemRPValidation/RPReconstructedTracksValidation/interface/RPTrackInfo.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDigCluster.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/TotemL1Trigger/interface/RPCCId.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"


RPReconstructedTracksValidation::RPReconstructedTracksValidation(const edm::ParameterSet& conf)
 : rp_histograms_("/RP_Tracks/RP_", edm::ParameterSet()), 
   det_histograms_("/Dets/Det_", edm::ParameterSet()),
   station_histograms_("/Stations/Stat_", edm::ParameterSet()),
   resol_degrad_service_(conf),
   conf_(conf)
{
  gSystem->Load("libHistPainter.so");
  hist_file_name_ = conf.getParameter<std::string>("HistogramFileName");
  OriginalHepMCProductLabel_ = conf.getParameter<std::string>("OriginalHepMCProductLabel");
  OriginalHepMCModuleName_ = conf.getParameter<std::string>("OriginalHepMCModuleName");
  SmearedHepMCProductLabel_ = conf.getParameter<std::string>("SmearedHepMCProductLabel");
  SmearedHepMCModuleName_ = conf.getParameter<std::string>("SmearedHepMCModuleName");
  verbosity_ = conf.getParameter<int>("Verbosity");
  SDValidation_ = conf.getParameter<bool>("SDValidation");
  modulLabelSimu_   = conf.getParameter<std::string>("ModulLabelSimu");
  productLabelSimu_ = conf.getParameter<std::string>("ProductLabelSimu");
  
  rpFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
  rpDigClusterSetLabel = conf.getParameter<edm::InputTag>("RPDigClusterSetLabel");
  rpDetTriggerSetLabel = conf.getParameter<edm::InputTag>("RPDetTriggerSetLabel");


  //init rp ids for the stations
  //order: from ip to outside, then top bottom
  //top first vertical, bottom first vertical, first horizontal, second horizontal, last top vertical, lat bottom vertical 
  for(int i=100; i<=105; ++i)
    right_rp_ids_.push_back(i);
  for(int i=120; i<=125; ++i)
    right_rp_ids_.push_back(i);
  
  for(int i=0; i<=5; ++i)
    left_rp_ids_.push_back(i);
  for(int i=20; i<=25; ++i)
    left_rp_ids_.push_back(i);
     
  for(int i=100; i<=105; ++i)
    station_150_right_ids_.push_back(i);
  for(int i=0; i<=5; ++i)
    station_150_left_ids_.push_back(i);
  for(int i=120; i<=125; ++i)
    station_220_right_ids_.push_back(i);
  for(int i=20; i<=25; ++i)
    station_220_left_ids_.push_back(i);
  
  stations_left_.push_back(0);
  stations_left_.push_back(2);
  stations_right_.push_back(10);
  stations_right_.push_back(12);
}


RPReconstructedTracksValidation::~RPReconstructedTracksValidation()
{
}


void RPReconstructedTracksValidation::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("RPReconstructedTracksValidation") << " edm::ESHandle<BeamOpticsParams> is invalid";
  BOPar_ = *BOParH;
  parameterized_madx_transport_ = std::auto_ptr<ParamMADRefTransport>(new ParamMADRefTransport(conf_, es));
  InitBeamSmearingHistograms();
  InitOverlapsHistograms();
  InitStationFlux();
  
  edm::ESHandle<TotemRPGeometry> Totem_RP_geometry;
  es.get<RealGeometryRecord>().get(Totem_RP_geometry);
  if(!Totem_RP_geometry.isValid()) 
  {
    throw cms::Exception("RPReconstructedTracksValidation") << "edm::ESHandle<TotemRPGeometry> is invalid, RP contours drawing impossible!";
  }
  
  Totem_RP_geometry_ = *Totem_RP_geometry;
}

void RPReconstructedTracksValidation::InitOverlapsHistograms()
{
  char name[1024];
  
  int hist_bins = 160;
  
  sprintf(name, "overlaps_RP_120_121_122_");
  overlaps_RP_120_121_122_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_120_121_122_->SetDirectory(0);
  overlaps_RP_120_121_122_->SetXTitle("x [mm]");
  overlaps_RP_120_121_122_->SetYTitle("y [mm]");

  sprintf(name, "overlaps_RP_100_101_102_");
  overlaps_RP_100_101_102_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_100_101_102_->SetDirectory(0);
  overlaps_RP_100_101_102_->SetXTitle("x [mm]");
  overlaps_RP_100_101_102_->SetYTitle("y [mm]");
  
  sprintf(name, "overlaps_RP_020_021_022_");
  overlaps_RP_020_021_022_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_020_021_022_->SetDirectory(0);
  overlaps_RP_020_021_022_->SetXTitle("x [mm]");
  overlaps_RP_020_021_022_->SetYTitle("y [mm]");
  
  sprintf(name, "overlaps_RP_000_001_002_");
  overlaps_RP_000_001_002_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_000_001_002_->SetDirectory(0);
  overlaps_RP_000_001_002_->SetXTitle("x [mm]");
  overlaps_RP_000_001_002_->SetYTitle("y [mm]");
  
  sprintf(name,
      "overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_");
  overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_->SetDirectory(0);
  overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_->SetXTitle("x [mm]");
  overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_->SetYTitle("y [mm]");
  
  sprintf(name,
      "overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_");
  overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_->SetDirectory(0);
  overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_->SetXTitle("x [mm]");
  overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_->SetYTitle("y [mm]");
  
  sprintf(name,
      "overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_");
  overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_->SetDirectory(0);
  overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_->SetXTitle("x [mm]");
  overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_->SetYTitle("y [mm]");
  
  sprintf(name,
      "overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_");
  overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_->SetDirectory(0);
  overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_->SetXTitle("x [mm]");
  overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_->SetYTitle("y [mm]");
}


void RPReconstructedTracksValidation::InitBeamSmearingHistograms()
{
  vertex_x_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("vertex_x_smearing_dist_", "vertex_x_smearing_dist_", 100,
          BOPar_.GetBeamDisplacementX()*1000, BOPar_.GetBeamDisplacementX()*1000+0.001));
  vertex_x_smearing_dist_->SetDirectory(0);
  vertex_x_smearing_dist_->SetXTitle("#Delta x [mm]");
  vertex_x_smearing_dist_->SetYTitle("Entries");

  vertex_y_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("vertex_y_smearing_dist_", "vertex_y_smearing_dist_", 100,
          BOPar_.GetBeamDisplacementY()*1000, BOPar_.GetBeamDisplacementY()*1000+0.001));
  vertex_y_smearing_dist_->SetDirectory(0);
  vertex_y_smearing_dist_->SetXTitle("#Delta y [mm]");
  vertex_y_smearing_dist_->SetYTitle("Entries");
  
  vertex_z_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("vertex_z_smearing_dist_", "vertex_z_smearing_dist_", 100,
          BOPar_.GetBeamDisplacementZ()*1000, BOPar_.GetBeamDisplacementZ()*1000+0.001));
  vertex_z_smearing_dist_->SetDirectory(0);
  vertex_z_smearing_dist_->SetXTitle("#Delta z [mm]");
  vertex_z_smearing_dist_->SetYTitle("Entries");
  
  px_right_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("px_right_smearing_dist_", "px_right_smearing_dist_", 100,
          BOPar_.GetCrossingAngleX()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetCrossingAngleX()*BOPar_.GetBeamMomentum()+0.001));
  px_right_smearing_dist_->SetDirectory(0);
  px_right_smearing_dist_->SetXTitle("#Delta p_{x} [GeV/c]");
  px_right_smearing_dist_->SetYTitle("Entries");
  
  py_right_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("py_right_smearing_dist_", "py_right_smearing_dist_", 100,
          BOPar_.GetCrossingAngleY()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetCrossingAngleY()*BOPar_.GetBeamMomentum()+0.001));
  py_right_smearing_dist_->SetDirectory(0);
  py_right_smearing_dist_->SetXTitle("#Delta p_{y} [GeV/c]");
  py_right_smearing_dist_->SetYTitle("Entries");
  
  pz_right_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("pz_right_smearing_dist_", "pz_right_smearing_dist_", 100,
          BOPar_.GetMeanXi()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetMeanXi()*BOPar_.GetBeamMomentum()+0.0001));
  pz_right_smearing_dist_->SetDirectory(0);
  pz_right_smearing_dist_->SetXTitle("#Delta p_{z} [GeV/c]");
  pz_right_smearing_dist_->SetYTitle("Entries");
  
  px_left_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("px_left_smearing_dist_", "px_left_smearing_dist_", 100,
          BOPar_.GetCrossingAngleX()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetCrossingAngleX()*BOPar_.GetBeamMomentum()+0.001));
  px_left_smearing_dist_->SetDirectory(0);
  px_left_smearing_dist_->SetXTitle("#Delta p_{x} [GeV/c]");
  px_left_smearing_dist_->SetYTitle("Entries");
  
  py_left_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("py_left_smearing_dist_", "py_left_smearing_dist_", 100,
          BOPar_.GetCrossingAngleY()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetCrossingAngleY()*BOPar_.GetBeamMomentum()+0.001));
  py_left_smearing_dist_->SetDirectory(0);
  py_left_smearing_dist_->SetXTitle("#Delta p_{y} [GeV/c]");
  py_left_smearing_dist_->SetYTitle("Entries");
  
  pz_left_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("pz_left_smearing_dist_", "pz_left_smearing_dist_", 100,
          -BOPar_.GetMeanXi()*BOPar_.GetBeamMomentum(), 
          -BOPar_.GetMeanXi()*BOPar_.GetBeamMomentum()+0.0001));
  pz_left_smearing_dist_->SetDirectory(0);
  pz_left_smearing_dist_->SetXTitle("#Delta p_{z} [GeV/c]");
  pz_left_smearing_dist_->SetYTitle("Entries");
  
  thx_right_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("thx_right_smearing_dist_", "thx_right_smearing_dist_", 100,
          BOPar_.GetCrossingAngleX(), 
          BOPar_.GetCrossingAngleX()+1e-7));
  thx_right_smearing_dist_->SetDirectory(0);
  thx_right_smearing_dist_->SetXTitle("#Delta #Theta_{x} [rad]");
  thx_right_smearing_dist_->SetYTitle("Entries");
  
  thy_right_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("thy_right_smearing_dist_", "thy_right_smearing_dist_", 100,
          BOPar_.GetCrossingAngleY(), 
          BOPar_.GetCrossingAngleY()+1e-7));
  thy_right_smearing_dist_->SetDirectory(0);
  thy_right_smearing_dist_->SetXTitle("#Delta #Theta_{y} [rad]");
  thy_right_smearing_dist_->SetYTitle("Entries");
  
  thx_left_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("thx_left_smearing_dist_", "thx_left_smearing_dist_", 100,
          BOPar_.GetCrossingAngleX(), 
          BOPar_.GetCrossingAngleX()+1e-7));
  thx_left_smearing_dist_->SetDirectory(0);
  thx_left_smearing_dist_->SetXTitle("#Delta #Theta_{x} [rad]");
  thx_left_smearing_dist_->SetYTitle("Entries");
  
  thy_left_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("thy_left_smearing_dist_", "thy_left_smearing_dist_", 100,
          BOPar_.GetCrossingAngleY(), 
          BOPar_.GetCrossingAngleY()+1e-7));
  thy_left_smearing_dist_->SetDirectory(0);
  thy_left_smearing_dist_->SetXTitle("#Delta #Theta_{y} [rad]");
  thy_left_smearing_dist_->SetYTitle("Entries");
  
  xi_right_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("xi_right_smearing_dist_", "xi_right_smearing_dist_", 100,
          BOPar_.GetMeanXi(), 
          BOPar_.GetMeanXi()+1e-5));
  xi_right_smearing_dist_->SetDirectory(0);
  xi_right_smearing_dist_->SetXTitle("#Delta #xi");
  xi_right_smearing_dist_->SetYTitle("Entries");
  
  xi_left_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("xi_left_smearing_dist_", "xi_left_smearing_dist_", 100,
          BOPar_.GetMeanXi(), 
          BOPar_.GetMeanXi()+1e-5));
  xi_left_smearing_dist_->SetDirectory(0);
  xi_left_smearing_dist_->SetXTitle("#Delta #xi");
  xi_left_smearing_dist_->SetYTitle("Entries");
}

void RPReconstructedTracksValidation::InitStationFlux()
{
  char name[1024];
  int hist_bins = 80;
  
  sprintf(name, "station_flux_RP_120_121_122_");
  station_flux_RP_120_121_122_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  station_flux_RP_120_121_122_->SetDirectory(0);
  station_flux_RP_120_121_122_->SetXTitle("x [mm]");
  station_flux_RP_120_121_122_->SetYTitle("y [mm]");
  
  sprintf(name, "station_flux_RP_100_101_102_");
  station_flux_RP_100_101_102_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  station_flux_RP_100_101_102_->SetDirectory(0);
  station_flux_RP_100_101_102_->SetXTitle("x [mm]");
  station_flux_RP_100_101_102_->SetYTitle("y [mm]");
  
  sprintf(name, "station_flux_RP_020_021_022_");
  station_flux_RP_020_021_022_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  station_flux_RP_020_021_022_->SetDirectory(0);
  station_flux_RP_020_021_022_->SetXTitle("x [mm]");
  station_flux_RP_020_021_022_->SetYTitle("y [mm]");
  
  sprintf(name, "station_flux_RP_000_001_002_");
  station_flux_RP_000_001_002_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  station_flux_RP_000_001_002_->SetDirectory(0);
  station_flux_RP_000_001_002_->SetXTitle("x [mm]");
  station_flux_RP_000_001_002_->SetYTitle("y [mm]");
}

void RPReconstructedTracksValidation::WriteBeamSmearingHistograms(TFile *f)
{
  f->cd("/");
  f->mkdir("beam_smearing_validation");
  f->cd("beam_smearing_validation");
  vertex_x_smearing_dist_->Write();
  vertex_y_smearing_dist_->Write();
  vertex_z_smearing_dist_->Write();
  
  px_right_smearing_dist_->Write();
  py_right_smearing_dist_->Write();
  pz_right_smearing_dist_->Write();
  
  px_left_smearing_dist_->Write();
  py_left_smearing_dist_->Write();
  pz_left_smearing_dist_->Write();
  
  thx_right_smearing_dist_->Write();
  thy_right_smearing_dist_->Write();
  
  thx_left_smearing_dist_->Write();
  thy_left_smearing_dist_->Write();
  
  xi_right_smearing_dist_->Write();
  xi_left_smearing_dist_->Write();
  f->cd("/");
}

void RPReconstructedTracksValidation::InitGeneratorDebugHistograms()
{
  char name[1024];
  int hist_bins = 200;
  
  sprintf(name, "fractional_momentum_left_vs_right_");
  fractional_momentum_left_vs_right_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -9, 0, hist_bins, -9, 0));
  fractional_momentum_left_vs_right_->SetDirectory(0);
  fractional_momentum_left_vs_right_->SetYTitle("Log(-#xi_{right})");
  fractional_momentum_left_vs_right_->SetYTitle("Log(-#xi_{left})");
}

void RPReconstructedTracksValidation::WriteGeneratorDebugHistograms(TFile *f)
{
  f->cd("/");
  f->mkdir("generator_debug_histograms");
  f->cd("generator_debug_histograms");
  
  fractional_momentum_left_vs_right_->Write();
  
  f->cd("/");
}

void RPReconstructedTracksValidation::FillGeneratorDebugHistogram_InserProtPair(const PrimaryProton &prim_prot_left_, const PrimaryProton &prim_prot_right_)
{
  if(prim_prot_left_.found && prim_prot_right_.found)
  {
    double pxl = prim_prot_left_.momentum.x();
    double pyl = prim_prot_left_.momentum.y();
    double pzl = prim_prot_left_.momentum.z();
    
    double pxr = prim_prot_right_.momentum.x();
    double pyr = prim_prot_right_.momentum.y();
    double pzr = prim_prot_right_.momentum.z();
    
    double xil = BOPar_.IPNonSmearedProtonMomentumToXi(pxl, pyl, pzl);
    double xir = BOPar_.IPNonSmearedProtonMomentumToXi(pxr, pyr, pzr);
    
    double logxil = TMath::Log10(-xil);
    double logxir = TMath::Log10(-xir);
    
    fractional_momentum_left_vs_right_->Fill(logxil, logxir);
  }
}

void RPReconstructedTracksValidation::WriteOverlapsHistograms(TFile *f)
{
  std::vector<RPId> rpids_st_12; 
  rpids_st_12.push_back(120); rpids_st_12.push_back(121); rpids_st_12.push_back(122);
  
  std::vector<RPId> rpids_st_10;
  rpids_st_10.push_back(100); rpids_st_10.push_back(101); rpids_st_10.push_back(102); 

  std::vector<RPId> rpids_st_02;
  rpids_st_02.push_back(20); rpids_st_02.push_back(21); rpids_st_02.push_back(22);
  
  std::vector<RPId> rpids_st_00;
  rpids_st_00.push_back(0); rpids_st_00.push_back(1); rpids_st_00.push_back(2); 
  
  f->cd("/");
  f->mkdir("overlapping_proton_flux");
  f->cd("overlapping_proton_flux");
  
  overlaps_RP_120_121_122_->Write();
  WriteHistogramAsCanvas((&*overlaps_RP_120_121_122_), rpids_st_12);
  overlaps_RP_100_101_102_->Write();
  WriteHistogramAsCanvas((&*overlaps_RP_100_101_102_), rpids_st_10);
  overlaps_RP_020_021_022_->Write();
  WriteHistogramAsCanvas((&*overlaps_RP_020_021_022_), rpids_st_02);
  overlaps_RP_000_001_002_->Write();
  WriteHistogramAsCanvas((&*overlaps_RP_000_001_002_), rpids_st_00);
  
  overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_->Write();
  WriteHistogramAsCanvas((&*overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_), rpids_st_12);
  overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_->Write();
  WriteHistogramAsCanvas((&*overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_), rpids_st_10);
  overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_->Write();
  WriteHistogramAsCanvas((&*overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_), rpids_st_02);
  overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_->Write();
  WriteHistogramAsCanvas((&*overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_), rpids_st_00);
  
  f->cd("/");
}

void RPReconstructedTracksValidation::WriteStationFlux(TFile *f)
{
  std::vector<RPId> rpids_st_12; 
  rpids_st_12.push_back(120); rpids_st_12.push_back(121); rpids_st_12.push_back(122);
  
  std::vector<RPId> rpids_st_10;
  rpids_st_10.push_back(100); rpids_st_10.push_back(101); rpids_st_10.push_back(102); 

  std::vector<RPId> rpids_st_02;
  rpids_st_02.push_back(20); rpids_st_02.push_back(21); rpids_st_02.push_back(22);
  
  std::vector<RPId> rpids_st_00;
  rpids_st_00.push_back(0); rpids_st_00.push_back(1); rpids_st_00.push_back(2); 
  
  f->cd("/");
  f->mkdir("station_proton_flux");
  f->cd("station_proton_flux");
  
  station_flux_RP_120_121_122_->Write();
  WriteHistogramAsCanvas((&*station_flux_RP_120_121_122_), rpids_st_12);
  
  station_flux_RP_100_101_102_->Write();
  WriteHistogramAsCanvas((&*station_flux_RP_100_101_102_), rpids_st_10);
  
  station_flux_RP_020_021_022_->Write();
  WriteHistogramAsCanvas((&*station_flux_RP_020_021_022_), rpids_st_02);
  
  station_flux_RP_000_001_002_->Write();
  WriteHistogramAsCanvas((&*station_flux_RP_000_001_002_), rpids_st_02);
  
  f->cd("/");
}


void RPReconstructedTracksValidation::WriteHistogramAsCanvas(TH1 *h, std::vector<RPId> rpids)
{
  char buf[1024];
  sprintf(buf, "%s_canv", h->GetName());
  TCanvas * c1 = new TCanvas(buf,buf,0,0,400,400);

  c1->SetTickx(1);
  c1->SetTicky(1);
  c1->cd();
  gStyle->SetPalette(1);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetPadBorderMode(0);
  gStyle->SetPadColor(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetStatColor(0);
  gStyle->SetFillColor(0);
  
  h->SetFillColor(0);
  c1->SetFrameFillColor(0);
  c1->SetFrameFillColor(0);
  
  h->Draw("colz");
  
  for(unsigned int i=0; i<rpids.size(); ++i)
  {
    DrawRPContour(rpids[i]);
  }

  c1->Modified();
  c1->Update();
  c1->Write("");
}

void RPReconstructedTracksValidation::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  //fill track and cluster related histograms 
  edm::Handle< RPFittedTrackCollection > input;
  e.getByLabel(rpFittedTrackCollectionLabel, input);
  if (!input.isValid()) 
    throw cms::Exception("RPReconstructedTracksValidation") << "edm::Handle< RPFittedTrackCollection > is invalid";
  
  edm::Handle< cluster_set > det_clusters;
  e.getByLabel(rpDigClusterSetLabel, det_clusters);
  if (!det_clusters.isValid()) 
    throw cms::Exception("RPReconstructedTracksValidation") << "edm::Handle<edm::det_clusters> det_clusters is invalid";
  
  edm::Handle<edm::DetSetVector<RPDetTrigger> > det_triggers;
  e.getByLabel(rpDetTriggerSetLabel, det_triggers);
//  if (!det_triggers.isValid()) 
//    throw cms::Exception("RPReconstructedTracksValidation") << "edm::Handle<DetSetVector<RPDetTrigger> > det_triggers is invalid";
  
  edm::Handle<std::vector<RPCCBits> > cc_chip_bits;
  e.getByLabel(modulLabelSimu_, productLabelSimu_, cc_chip_bits ); 
    
  //Fill info for all reconstructed tracks (no cuts, all track included)
  FillTrackPositions(input, cc_chip_bits);
  FillDetectorHistogramsForReconstructedTracks(input, det_clusters, det_triggers);
  FillOverlappingAreas(*input);
  FillStationFlux(*input);
  
  //selects hits for stations (protons in 2 units requested per station), reconstruction cut applied
  rec_tracks_collection hits_220_r;
  rec_tracks_collection hits_220_l;
  rec_tracks_collection hits_150_r;
  rec_tracks_collection hits_150_l;   
  int no_of_pots_220_r = Select220RightHits(*input, hits_220_r);
  int no_of_pots_150_r = Select150RightHits(*input, hits_150_r);
  int no_of_pots_220_l = Select220LeftHits(*input, hits_220_l);
  int no_of_pots_150_l = Select150LeftHits(*input, hits_150_l);
  
  //select arm protons, all track included, tracks just grouped per arm
  rec_tracks_collection hits_r;
  rec_tracks_collection hits_l;
  SelectArmHits(*input, hits_l, hits_r);
  
  //fill station protons, at least 2 RP tracks required per station
  FillStationProtonTracks(hits_220_r);
  FillStationProtonTracks(hits_150_r);
  FillStationProtonTracks(hits_220_l);
  FillStationProtonTracks(hits_150_l);
  
  //fill acceptance and smearing information if available
  edm::Handle<edm::HepMCProduct> OriginalHepMCEvt;
  e.getByLabel(OriginalHepMCModuleName_, OriginalHepMCProductLabel_, OriginalHepMCEvt ); 
  if(!OriginalHepMCEvt.isValid())
  {
    if (verbosity_ )
      std::cout<<"Primary protons not found!! Skipping the beam related event validation."<<std::endl;
    return;
  }
  
  edm::Handle<edm::HepMCProduct> SmearedHepMCEvt;
  e.getByLabel(SmearedHepMCModuleName_, SmearedHepMCProductLabel_, SmearedHepMCEvt );
  if (!SmearedHepMCEvt.isValid())
  {
    if (verbosity_ )
      std::cout<<"Smeared protons not found!! Skipping the beam related event validation."<<std::endl;
    return;
  }    
  
  bool prim_found = FindPrimaryProtons(OriginalHepMCEvt);
  if (verbosity_ && !prim_found)
  {
    std::cout<<"Primary protons not found!! Skipping the beam related event validation."<<std::endl;
  }
  if(!prim_found)
  {
    return;
  }
  
  bool smeared_found = FindSmearedProtons(SmearedHepMCEvt);
  if (verbosity_ && !smeared_found)
  {
    std::cout<<"Smeared protons not found!! Skipping the beam related event validation."<<std::endl;
  }
  if(!smeared_found)
  {
    return;
  }

  if(right_prim_prot_.found)
  {
    FillReferenceStationHistograms(stations_right_, right_prim_prot_);
    for(unsigned int i=0; i<right_rp_ids_.size(); ++i)
    {
      FillReferenceHistograms(right_rp_ids_[i], right_prim_prot_);
    }
  }
  if(left_prim_prot_.found)
  {
    FillReferenceStationHistograms(stations_left_, left_prim_prot_);
    for(unsigned int i=0; i<left_rp_ids_.size(); ++i)
    {
      FillReferenceHistograms(left_rp_ids_[i], left_prim_prot_);
    }
  }
  
  //FillGeneratorDebugHistogram_InserProtPair(left_prim_prot_, right_prim_prot_);
  FillEventSmearingInformation();

  //For all reconstructed tracks and all pots, no cuts applied, physics acceptance info included
  FillResidualHistograms(hits_r, right_prim_prot_, right_smeared_prot_, no_of_pots_150_r==2);
  FillResidualHistograms(hits_l, left_prim_prot_, left_smeared_prot_, no_of_pots_150_l==2);

  //All stations included, cuts applied to have at least two tracks per station
  if(no_of_pots_220_r)
  {
    FillStationHistograms(hits_220_r, right_prim_prot_, right_smeared_prot_, no_of_pots_150_r==2);
  }
  if(no_of_pots_220_l)
  {
    FillStationHistograms(hits_220_l, left_prim_prot_, left_smeared_prot_, no_of_pots_150_l==2);
  }
  if(no_of_pots_150_r)
  {
    FillStationHistograms(hits_150_r, right_prim_prot_, right_smeared_prot_, no_of_pots_150_r==2);
  }
  if(no_of_pots_150_l)
  {
    FillStationHistograms(hits_150_l, left_prim_prot_, left_smeared_prot_, no_of_pots_150_l==2);
  }
}


void RPReconstructedTracksValidation::FillOverlappingAreas(const RPFittedTrackCollection &tracks)
{
  RPFittedTrackCollection::const_iterator it1, it2;
  
  //Fill overlaps_RP_120_121_122_
  if((it1=tracks.find(120))!=tracks.end() && (it2=tracks.find(122))!=tracks.end())
    overlaps_RP_120_121_122_->Fill(it1->second.X0(), it1->second.Y0());
  if((it1=tracks.find(121))!=tracks.end() && (it2=tracks.find(122))!=tracks.end())
    overlaps_RP_120_121_122_->Fill(it1->second.X0(), it1->second.Y0());

  //Fill overlaps_RP_100_101_102_
  if((it1=tracks.find(100))!=tracks.end() && (it2=tracks.find(102))!=tracks.end())
    overlaps_RP_100_101_102_->Fill(it1->second.X0(), it1->second.Y0());
  if((it1=tracks.find(101))!=tracks.end() && (it2=tracks.find(102))!=tracks.end())
    overlaps_RP_100_101_102_->Fill(it1->second.X0(), it1->second.Y0());
  
  //Fill overlaps_RP_020_021_022_
  if((it1=tracks.find(20))!=tracks.end() && (it2=tracks.find(22))!=tracks.end())
    overlaps_RP_020_021_022_->Fill(it1->second.X0(), it1->second.Y0());
  if((it1=tracks.find(21))!=tracks.end() && (it2=tracks.find(22))!=tracks.end())
    overlaps_RP_020_021_022_->Fill(it1->second.X0(), it1->second.Y0());
  
  //Fill overlaps_RP_000_001_002_
  if((it1=tracks.find(0))!=tracks.end() && (it2=tracks.find(2))!=tracks.end())
    overlaps_RP_000_001_002_->Fill(it1->second.X0(), it1->second.Y0());
  if((it1=tracks.find(1))!=tracks.end() && (it2=tracks.find(2))!=tracks.end())
    overlaps_RP_000_001_002_->Fill(it1->second.X0(), it1->second.Y0());
  
  //Fill interstation overlaps 
  //overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_
  //overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_
  for(int i = 120; i<=122; ++i)
  {
    it1 = tracks.find(i);
    if(it1!=tracks.end())
      break;
  }
  for(int i = 100; i<=102; ++i)
  {
    it2 = tracks.find(i);
    if(it2!=tracks.end())
      break;
  }
  if(it1!=tracks.end() && it2!=tracks.end())
  {
    overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_->Fill(it1->second.X0(), it1->second.Y0());
    overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_->Fill(it2->second.X0(), it2->second.Y0());
  }
  
  //Fill interstation overlaps 
  //overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_;
  //overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_;
  for(int i = 20; i<=22; ++i)
  {
    it1 = tracks.find(i);
    if(it1!=tracks.end())
      break;
  }
  for(int i = 0; i<=2; ++i)
  {
    it2 = tracks.find(i);
    if(it2!=tracks.end())
      break;
  }
  if(it1!=tracks.end() && it2!=tracks.end())
  {
    overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_->Fill(it1->second.X0(), it1->second.Y0());
    overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_->Fill(it2->second.X0(), it2->second.Y0());
  }
}


void RPReconstructedTracksValidation::FillTrackPositions(edm::Handle< RPFittedTrackCollection > &input, edm::Handle<std::vector<RPCCBits> > &cc_chip_bits)
{
  RPFittedTrackCollection::const_iterator it;
  for(it=(*input).begin(); it!=(*input).end(); ++it)
  {
    rp_histograms_.GetObj(it->first)->FillBasicTrackInfo(it->second);
  }
  
  if(cc_chip_bits.isValid())
  {    //{u,v}
    std::map<RPId, std::pair<RPCCBits,RPCCBits> > cc_bits;
    
    std::vector<RPCCBits>::const_iterator itc;
    for(itc = cc_chip_bits->begin(); itc!=cc_chip_bits->end(); ++itc)
    {
      RPCCId rawid(itc->getId());
      RPId decid = rawid.Arm()*100 + rawid.Station()*10 + rawid.RomanPot(); 
      if(rawid.IsStripsCoordinateUDirection())
      {
        cc_bits[decid].first = *itc;
      }
      else
      {
        cc_bits[decid].second = *itc;
      }
    }
    
    std::map<RPId, std::pair<RPCCBits,RPCCBits> >::iterator u_it;
    for(u_it = cc_bits.begin(); u_it!=cc_bits.end(); ++u_it)
    {
      rp_histograms_.GetObj(u_it->first)->FillRPCoincidenceChipInfo(u_it->second.first, u_it->second.second);      
    }
  }
}


void RPReconstructedTracksValidation::FillEventSmearingInformation()
{
  bool vertex_info_filled = false;
  if(right_prim_prot_.found && right_smeared_prot_.found)
  {
    vertex_x_smearing_dist_->Fill(right_smeared_prot_.vertex.x()-
        right_prim_prot_.vertex.x());
    vertex_y_smearing_dist_->Fill(right_smeared_prot_.vertex.y()-
        right_prim_prot_.vertex.y());
    vertex_z_smearing_dist_->Fill(right_smeared_prot_.vertex.z()-
        right_prim_prot_.vertex.z());
    vertex_info_filled = true;
    
    px_right_smearing_dist_->Fill(right_smeared_prot_.momentum.px()-
        right_prim_prot_.momentum.px());
    py_right_smearing_dist_->Fill(right_smeared_prot_.momentum.py()-
        right_prim_prot_.momentum.py());
    pz_right_smearing_dist_->Fill(right_smeared_prot_.momentum.pz()-
        right_prim_prot_.momentum.pz());
    
    thx_right_smearing_dist_->Fill(
        right_smeared_prot_.momentum.px()/TMath::Abs(right_smeared_prot_.momentum.pz())-
        right_prim_prot_.momentum.px()/TMath::Abs(right_prim_prot_.momentum.pz()));
    
    
    thy_right_smearing_dist_->Fill(
        right_smeared_prot_.momentum.py()/TMath::Abs(right_smeared_prot_.momentum.pz())-
        right_prim_prot_.momentum.py()/TMath::Abs(right_prim_prot_.momentum.pz()));
    xi_right_smearing_dist_->Fill(right_smeared_prot_.momentum.rho()/BOPar_.GetBeamMomentum()-
        right_prim_prot_.momentum.rho()/BOPar_.GetBeamMomentum());
  }

  if(left_prim_prot_.found && left_smeared_prot_.found)
  {
    if(!vertex_info_filled)
    {
      vertex_x_smearing_dist_->Fill(left_smeared_prot_.vertex.x()-
          left_prim_prot_.vertex.x());
      vertex_y_smearing_dist_->Fill(left_smeared_prot_.vertex.y()-
          left_prim_prot_.vertex.y());
      vertex_z_smearing_dist_->Fill(left_smeared_prot_.vertex.z()-
          left_prim_prot_.vertex.z());
    }
    
    px_left_smearing_dist_->Fill(left_smeared_prot_.momentum.px()-
        left_prim_prot_.momentum.px());
    py_left_smearing_dist_->Fill(left_smeared_prot_.momentum.py()-
        left_prim_prot_.momentum.py());
    pz_left_smearing_dist_->Fill(left_smeared_prot_.momentum.pz()-
        left_prim_prot_.momentum.pz());
    
    
    thx_left_smearing_dist_->Fill(
        left_smeared_prot_.momentum.px()/TMath::Abs(left_smeared_prot_.momentum.pz())-
        left_prim_prot_.momentum.px()/TMath::Abs(left_prim_prot_.momentum.pz()));
    thy_left_smearing_dist_->Fill(
        left_smeared_prot_.momentum.py()/TMath::Abs(left_smeared_prot_.momentum.pz())-
        left_prim_prot_.momentum.py()/TMath::Abs(left_prim_prot_.momentum.pz()));
    xi_left_smearing_dist_->Fill(left_smeared_prot_.momentum.rho()/BOPar_.GetBeamMomentum()-
        left_prim_prot_.momentum.rho()/BOPar_.GetBeamMomentum());
  }
}


void RPReconstructedTracksValidation::FillStationFlux(
    const RPFittedTrackCollection &tracks)
{
  RPFittedTrackCollection::const_iterator it1;
  
  //station_flux_RP_120_121_122_
  for (int i = 120; i<=122; ++i)
  {
    it1 = tracks.find(i);
    if (it1!=tracks.end())
    {
      station_flux_RP_120_121_122_->Fill(it1->second.X0(), it1->second.Y0());
      break;
    }
  }
  
  //station_flux_RP_100_101_102_
  for (int i = 100; i<=102; ++i)
  {
    it1 = tracks.find(i);
    if (it1!=tracks.end())
    {
      station_flux_RP_100_101_102_->Fill(it1->second.X0(), it1->second.Y0());
      break;
    }
  }
  
  //station_flux_RP_020_021_022_
  for (int i = 20; i<=22; ++i)
  {
    it1 = tracks.find(i);
    if (it1!=tracks.end())
    {
      station_flux_RP_020_021_022_->Fill(it1->second.X0(), it1->second.Y0());
      break;
    }
  }
  
  //station_flux_RP_000_001_002_
  for (int i = 0; i<=2; ++i)
  {
    it1 = tracks.find(i);
    if (it1!=tracks.end())
    {
      station_flux_RP_000_001_002_->Fill(it1->second.X0(), it1->second.Y0());
      break;
    }
  }
}


void RPReconstructedTracksValidation::FillDetectorHistogramsForReconstructedTracks(
    const edm::Handle<RPFittedTrackCollection> &tracks, const edm::Handle<cluster_set> &clusters, const edm::Handle< edm::DetSetVector<RPDetTrigger> > &det_triggers)
{
  if(tracks.isValid() && clusters.isValid())
  {
    for (RPFittedTrackCollection::const_iterator trit = tracks->begin(); trit!=tracks->end(); ++trit)
    {
      FillDetectorHistograms(trit->second, clusters);
    }
  }
  
  if(det_triggers.isValid())
  {
    int triggeredSectorNo, dec_det_id;
    edm::DetSetVector<RPDetTrigger>::const_iterator inputIterator = det_triggers->begin();
    
    //Fill trigger sectors
    for(; inputIterator != det_triggers->end(); inputIterator++) 
    {
      TotRPDetId simid(inputIterator->id);
      dec_det_id = simid.DetectorDecId();
  
      // inputIterator->data : vector< RPDetTrigger>
      // (inputIterator->data)[i] : RPDetTrigger
      for (unsigned int i=0; i < (inputIterator->data).size(); ++i) 
      {
        triggeredSectorNo = (unsigned int) ((inputIterator->data)[i].GetSector());
        det_histograms_.GetObj(dec_det_id)->FillTriggerSector(triggeredSectorNo);
      }
    }
  }
}


void RPReconstructedTracksValidation::FillReferenceHistograms(RPId id, const PrimaryProton &proton)
{
  boost::shared_ptr<RPTrackInfo> ptr = rp_histograms_.GetObj(id);
  
  double px = proton.momentum.x();
  double py = proton.momentum.y();
  double pz = proton.momentum.z();

  double xi = BOPar_.IPNonSmearedProtonMomentumToXi(px, py, pz);
  double t = BOPar_.IPNonSmearedProtonMomentumTot(px, py, pz);
  
//  double t1=-pz*pz*((px/pz)*(px/pz) + (py/pz)*(py/pz)); 
  ptr->Fill_t_Dist(t);
  ptr->Fill_ksi_t_Dist(t, xi);
  if(verbosity_>0)
  {
    std::cout<<"RPReconstructedTracksValidation::FillReferenceHistograms"<<std::endl;
    std::cout<<"id="<<id<<" px="<<px<<" py="<<py<<" pz="<<pz<<std::endl;
    std::cout<<"xi="<<xi<<" t="<<t<<std::endl;
  }
}

void RPReconstructedTracksValidation::FillReferenceStationHistograms(station_rp_ids_type stations, const PrimaryProton &proton)
{
  double px = proton.momentum.x();
  double py = proton.momentum.y();
  double pz = proton.momentum.z();

  double xi = BOPar_.IPNonSmearedProtonMomentumToXi(px, py, pz);
  double t = BOPar_.IPNonSmearedProtonMomentumTot(px, py, pz);
  
  for(station_rp_ids_type::iterator it = stations.begin(); it!=stations.end(); ++it)
  {
    boost::shared_ptr<RPStationInfo> ptr = station_histograms_.GetObj(*it); 
    ptr->Fill_t_Dist(t);
    ptr->Fill_ksi_t_Dist(t, xi);
    if(verbosity_>0)
    {
      std::cout<<"RPReconstructedTracksValidation::FillReferenceStationHistograms"<<std::endl;
      std::cout<<"st_id="<<*it<<" px="<<px<<" py="<<py<<" pz="<<pz<<std::endl;
      std::cout<<"xi="<<xi<<" t="<<t<<std::endl;
    }    
  }
}



//return number of pots through which the proton goes
int RPReconstructedTracksValidation::SelectHits(const RPFittedTrackCollection &tracks, 
    rec_tracks_collection & coll, const station_rp_ids_type& st_ids)
{
  RPFittedTrackCollection::const_iterator it;
  RPFittedTrackCollection::const_iterator low_it;
  RPFittedTrackCollection::const_iterator high_it;
  
  assert(st_ids.size() == 6);
  
  low_it = tracks.lower_bound(st_ids[0]);
  high_it = tracks.upper_bound(st_ids[5]);
  
  for(it = low_it; it!=high_it; ++it)
  {
    if(!it->second.IsValid())
      continue;
    coll[it->first]=resol_degrad_service_.Create2DHit(it->first, it->second);
  }
  
  rec_tracks_collection::const_iterator end;
  end = coll.end();
  
  bool station_accepted = 
    (coll.find(st_ids[0])!=end || coll.find(st_ids[1])!=end || coll.find(st_ids[2])!=end) && 
    (coll.find(st_ids[3])!=end || coll.find(st_ids[4])!=end || coll.find(st_ids[5])!=end) &&
    !(coll.find(st_ids[0])!=end && coll.find(st_ids[1])!=end) &&
    !(coll.find(st_ids[4])!=end && coll.find(st_ids[5])!=end);
  
  if(!station_accepted)
    coll.clear();
  
  return coll.size();
}

int RPReconstructedTracksValidation::SelectArmHits(const RPFittedTrackCollection &tracks, 
    rec_tracks_collection & coll_l, rec_tracks_collection & coll_r)
{
  RPFittedTrackCollection::const_iterator it;
  RPFittedTrackCollection::const_iterator low_it;
  RPFittedTrackCollection::const_iterator high_it;
  coll_l.clear();
  coll_r.clear();
  
  low_it = tracks.begin();
  high_it = tracks.end();
  for(it = low_it; it!=high_it; ++it)
  {
    if(!it->second.IsValid())
      continue;
    if( TotRPDetId::ArmOfRP(it->first)==0 )
    {
      coll_l[it->first]=resol_degrad_service_.Create2DHit(it->first, it->second);
    }
    else
    {
      coll_r[it->first]=resol_degrad_service_.Create2DHit(it->first, it->second);
    }
  }
  return coll_l.size() + coll_r.size(); 
}


int RPReconstructedTracksValidation::Select220RightHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll)
{
  return SelectHits(tracks, coll, station_220_right_ids_);
}


int RPReconstructedTracksValidation::Select220LeftHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll)
{
  return SelectHits(tracks, coll, station_220_left_ids_);
}


int RPReconstructedTracksValidation::Select150RightHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll)
{
  return SelectHits(tracks, coll, station_150_right_ids_);
}


int RPReconstructedTracksValidation::Select150LeftHits(const RPFittedTrackCollection &tracks, rec_tracks_collection & coll)
{
  return SelectHits(tracks, coll, station_150_left_ids_);
}


void RPReconstructedTracksValidation::FillStationProtonTracks(rec_tracks_collection & coll)
{
  if(coll.size()<2)
    return;
  
  rec_tracks_collection::iterator first_track = coll.begin();
  rec_tracks_collection::reverse_iterator last_track = coll.rbegin();
  
  RPStationId st_id1 = TotRPDetId::StOfRP(first_track->first);
  RPStationId st_id2 = TotRPDetId::StOfRP(last_track->first);
  
  if(st_id1!=st_id2)
    throw cms::Exception("RPReconstructedTracksValidation") << "Tracks wrongly classified to the station, tr1="<<first_track->first
    <<" tr2="<<last_track->first;

  double st_centre = Totem_RP_geometry_.GetStationCentreZPosition(st_id1);
  station_histograms_.GetObj(st_id1)->FillBasicTrackInfo(first_track->second, last_track->second, st_centre);
}


void RPReconstructedTracksValidation::FillStationAcceptance(RPStationId id, double t, double xi)
{
  station_histograms_.GetObj(id)->FillTrackSeen(t, xi);
  station_histograms_.GetObj(id)->FillTrackSeen(t);
}


void RPReconstructedTracksValidation::FillStationHistograms(const rec_tracks_collection &tracks, const PrimaryProton &prim_proton, 
    const PrimaryProton &smeared_proton, bool proton_at_150_station)
{
  if(!prim_proton.found || !smeared_proton.found)
    return;
  
  double opx = prim_proton.momentum.x();
  double opy = prim_proton.momentum.y();
  double opz = prim_proton.momentum.z();
  
  double xi = BOPar_.IPNonSmearedProtonMomentumToXi(opx, opy, opz);
  double t = BOPar_.IPNonSmearedProtonMomentumTot(opx, opy, opz);
  
  //fill station tracks
  rec_tracks_collection::const_iterator first_tr = tracks.begin(); 
  if(first_tr!=tracks.end())
  {
    RPStationId st_id = TotRPDetId::StOfRP(first_tr->first);
    FillStationAcceptance(st_id, t, xi);
  }
}


void RPReconstructedTracksValidation::FillResidualHistograms(
    const rec_tracks_collection &tracks, const PrimaryProton &prim_proton, 
    const PrimaryProton &smeared_proton, bool proton_at_150_station)
{
  if(!prim_proton.found || !smeared_proton.found)
    return;
  
  TVector2 rp_transverse_track_position;
  
  double spx = smeared_proton.momentum.x();
  double spy = smeared_proton.momentum.y();
  double spz = smeared_proton.momentum.z();
  
  double svx = smeared_proton.vertex.x();//+BOPar_.GetBeamDisplacementX()*1000.0;
  double svy = smeared_proton.vertex.y();//+BOPar_.GetBeamDisplacementY()*1000.0;
  double svz = smeared_proton.vertex.z();//+BOPar_.GetBeamDisplacementZ()*1000.0;

  double ovx = prim_proton.vertex.x();//+BOPar_.GetBeamDisplacementX()*1000.0;
  double ovy = prim_proton.vertex.y();//+BOPar_.GetBeamDisplacementY()*1000.0;
  double ovz = prim_proton.vertex.z();//+BOPar_.GetBeamDisplacementZ()*1000.0;
  
  double opx = prim_proton.momentum.x();
  double opy = prim_proton.momentum.y();
  double opz = prim_proton.momentum.z();
  
  double xi = BOPar_.IPNonSmearedProtonMomentumToXi(opx, opy, opz);
  double t = BOPar_.IPNonSmearedProtonMomentumTot(opx, opy, opz);
//    double t1=-opz*opz*((opx/opz)*(opx/opz) + (opy/opz)*(opy/opz));
  
  //fill residuals for each RP
  for(rec_tracks_collection::const_iterator it = tracks.begin(); 
        it!=tracks.end(); ++it)
  { 
    bool res = parameterized_madx_transport_->MADXTeoreticalRPTransversePosition(it->first, 
          svx, svy, svz, spx, spy, spz, rp_transverse_track_position, false);
    
    boost::shared_ptr<RPTrackInfo> ptr = rp_histograms_.GetObj(it->first);
    ptr->FillTrackSeen(t);
    ptr->FillTrackSeen(t, xi);
    
    if(verbosity_>0)
    {
      std::cout<<"RPReconstructedTracksValidation::FillResidualHistograms"<<std::endl;
      std::cout<<"PrimaryProton:"<<std::endl;
      std::cout<<"id="<<it->first<<" px="<<opx<<" py="<<opy<<" pz="<<opz<<std::endl;
      std::cout<<"id="<<it->first<<" vx="<<ovx<<" vy="<<ovy<<" vz="<<ovz<<std::endl;
      std::cout<<"xi="<<xi<<" t="<<t<<std::endl;
      std::cout<<"SmearedProton:"<<std::endl;
      std::cout<<"id="<<it->first<<" px="<<spx<<" py="<<spy<<" pz="<<spz<<std::endl;
      std::cout<<"id="<<it->first<<" vx="<<svx<<" vy="<<svy<<" vz="<<svz<<std::endl;
    }
    
    if(res)
    {
      double x_diff = it->second.X() - rp_transverse_track_position.X();
      double y_diff = it->second.Y() - rp_transverse_track_position.Y();
      ptr->FillRPPositionResiduals(x_diff, y_diff);
      ptr->FillPools(x_diff/it->second.Sx(), y_diff/it->second.Sy());
      
      if(verbosity_>0)
      {
        std::cout<<"x_diff="<<x_diff<<" y_diff="<<y_diff<<std::endl;
      }
            
      if(proton_at_150_station)
      {
        ptr->FillRPPositionResiduals2ProtonsSeenAt150Station(x_diff, y_diff);
      }
    }
  }
}


void RPReconstructedTracksValidation::FillDetectorHistograms(
        const RPFittedTrack &track, const edm::Handle< cluster_set > &clusters)
{
  for(int i = 0; i<track.GetHitEntries(); ++i)
  {
    RPDetId sim_det_id = track.GetHit(i).DetId();
    double residual = track.GetHit(i).Residual();
    TotRPDetId det_raw_id(sim_det_id);
    RPDetId dec_det_id = det_raw_id.DetectorDecId();
    cluster_set::const_iterator it = clusters->find(sim_det_id);
    if(it!=clusters->end())
    {
      det_histograms_.GetObj(dec_det_id)->FillHistograms(*it, residual);
    }
  }
}


bool RPReconstructedTracksValidation::FindPrimaryProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt)
{
  if(!HepMCEvt.isValid())
  {
     throw SimG4Exception("RPReconstructedTracksValidation : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");
  }
  
  const HepMC::GenEvent *evt = HepMCEvt->GetEvent();
  
  int right_count = 0;
  int left_count = 0;
  right_prim_prot_.found = false;
  left_prim_prot_.found = false;
  
  for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); 
      it != evt->particles_end(); ++it )
  {
    HepMC::GenParticle * g = (*it);
    int g_status = g->status();
    int pdg_id = g->pdg_id();
    
    if(verbosity_)
    {
      g->print(std::cout);
      if(g->production_vertex() != NULL){
    	const HepMC::FourVector &v = g->production_vertex()->position();
    	std::cout<< "[" << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.t() << "]"<<std::endl;
      }
    }
    
    // scanning only for particles with status == 1 
    if (g_status == 1 && pdg_id == 2212)
    {
      const HepMC::FourVector &vtx = g->production_vertex()->position();
      const HepMC::FourVector &mom  = g->momentum();
      
      if(verbosity_)
        std::cout<<"Setting the primary protons: mom="<<
        "[" << mom.x() << ", " << mom.y() << ", " << mom.z() << ", " << mom.t() << "]"
        <<" vert="
        <<"[" << vtx.x() << ", " << vtx.y() << ", " << vtx.z() << ", " << vtx.t() << "]"<<std::endl;
      
      if(mom.z()>0)
      {
        ++right_count;
        if(!right_prim_prot_.found || right_prim_prot_.momentum.rho()<mom.rho())
        {
          right_prim_prot_.vertex = vtx; //[mm]
          right_prim_prot_.momentum = mom;  //[GeV]
          right_prim_prot_.found = true;
        }
      }
      if(mom.z()<0)
      {
        ++left_count;
        if(!left_prim_prot_.found || left_prim_prot_.momentum.rho()<mom.rho())
        {
          left_prim_prot_.vertex = vtx; //[mm]
          left_prim_prot_.momentum = mom;  //[GeV]
          left_prim_prot_.found = true; 
        }
      }
    }
  } // end loop on HepMC particles
  
  if(SDValidation_ && right_prim_prot_.found && left_prim_prot_.found)
  {
    if(right_prim_prot_.momentum.rho()>left_prim_prot_.momentum.rho())
    {
      left_prim_prot_.found = false;
      left_count = 0;
    }
    else
    {
      right_prim_prot_.found = false;
      right_count = 0;
    }
  }
  bool result = right_count>0 || left_count>0;
  
  if(verbosity_)
    std::cout<<"right_count="<<right_count<<" left_count="<<left_count
        <<" result="<<result<<std::endl;
  return result;
}


bool RPReconstructedTracksValidation::FindSmearedProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt)
{
  if(!HepMCEvt.isValid())
  {
     throw SimG4Exception("RPReconstructedTracksValidation : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");
  }
  
  const HepMC::GenEvent *evt = HepMCEvt->GetEvent();
  
  int right_count = 0;
  int left_count = 0;
  right_smeared_prot_.found = false;
  left_smeared_prot_.found = false;
  
  for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); 
      it != evt->particles_end(); ++it )
  {
    HepMC::GenParticle * g = (*it);
    int g_status = g->status();
    int pdg_id = g->pdg_id();
    
    if(verbosity_)
    {
      g->print(std::cout);
      if(g->production_vertex() != NULL){
    	const HepMC::FourVector &v = g->production_vertex()->position();
    	std::cout<< "[" << v.x() << ", " << v.y() << ", " << v.z() << ", " << v.t() << "]"<<std::endl;
      }
    }
    
    // scanning only for particles with status == 1 
    if (g_status == 1 && pdg_id == 2212)
    {
      if(verbosity_)
        std::cout<<"Setting the smeared protons."<<std::endl;
      const HepMC::FourVector &vtx = g->production_vertex()->position();
      const HepMC::FourVector &mom  = g->momentum();
      
      if(mom.z()>0)
      {
        ++right_count;
        if(!right_smeared_prot_.found || right_smeared_prot_.momentum.rho()<mom.rho())
        {
          right_smeared_prot_.vertex = vtx; //[mm]
          right_smeared_prot_.momentum = mom;  //[GeV]
          right_smeared_prot_.found = true;
        }
      }
      if(mom.z()<0)
      {
        ++left_count;
        if(!left_smeared_prot_.found || left_smeared_prot_.momentum.rho()<mom.rho())
        {
          left_smeared_prot_.vertex = vtx; //[mm]
          left_smeared_prot_.momentum = mom;  //[GeV]
          left_smeared_prot_.found = true; 
        }
      }
    }
  } // end loop on HepMC particles
  
  if(SDValidation_ && right_smeared_prot_.found && left_smeared_prot_.found)
  {
    if(right_smeared_prot_.momentum.rho()>left_smeared_prot_.momentum.rho())
    {
      left_smeared_prot_.found = false;
      left_count = 0;
    }
    else
    {
      right_smeared_prot_.found = false;
      right_count = 0;
    }
  }
  bool result = right_count>0 || left_count>0;
  
  if(verbosity_)
    std::cout<<"right_count="<<right_count<<" left_count="<<left_count
        <<" result="<<result<<std::endl;
  return result;
}



bool RPReconstructedTracksValidation::IsRightRPId(RPId id)
{
  return TotRPDetId::ArmOfRP(id) == 1;
}


bool RPReconstructedTracksValidation::IsLeftRPId(RPId id)
{
  return TotRPDetId::ArmOfRP(id) == 0;
}


void RPReconstructedTracksValidation::endJob()
{
  WriteHistograms(hist_file_name_);
}

void RPReconstructedTracksValidation::DrawRPContour(RPId dec_rpid)
{
  unsigned int rawid = TotRPDetId::DecToRawId(10*dec_rpid);
  TotRPDetId rpid(rawid);
  
  CLHEP::Hep3Vector trans = Totem_RP_geometry_.GetDetTranslation(rpid);
  
  double xx, xy, xz, yx, yy, yz, zx, zy, zz;
  (Totem_RP_geometry_.GetDetector(rpid)->rotation()).GetComponents(xx, xy, xz, yx, yy, yz, zx, zy, zz);
  CLHEP::HepRep3x3 rot_mat( xx, xy, xz, yx, yy, yz, zx, zy, zz);
  CLHEP::HepRotation rot(rot_mat);
    
  CLHEP::Hep3Vector pA,pB,pC,pD,pE;
  // FIXME: these values need to match the XML geometry description
  // TODO: find a way to do it automatically
  double len = 36.070; // length of detector's edge
  double half_len = len * 0.5;
  double cut = len - 22.276 * sqrt(2.0) * 0.5; // length of the edge adjacent to the cut

  pA.setX(half_len);
  pA.setY(-half_len);

  pB.setX(half_len);
  pB.setY(half_len);

  pC.setX(-half_len);
  pC.setY(half_len);

  pD.setX(-half_len);
  pD.setY(half_len-cut);

  pE.setX(half_len-cut);
  pE.setY(-half_len);

  CLHEP::Hep3Vector pAn = rot * pA + trans;
  CLHEP::Hep3Vector pBn = rot * pB + trans;
  CLHEP::Hep3Vector pCn = rot * pC + trans;
  CLHEP::Hep3Vector pDn = rot * pD + trans;
  CLHEP::Hep3Vector pEn = rot * pE + trans;
  
  TLine *lAB = new TLine(pAn.x(),pAn.y(),pBn.x(),pBn.y());
  lAB->Draw();
  TLine *lBC = new TLine(pBn.x(),pBn.y(),pCn.x(),pCn.y());
  lBC->Draw();
  TLine *lCD = new TLine(pCn.x(),pCn.y(),pDn.x(),pDn.y());
  lCD->Draw();
  TLine *lDE = new TLine(pDn.x(),pDn.y(),pEn.x(),pEn.y());
  lDE->Draw();
  TLine *lEA = new TLine(pEn.x(),pEn.y(),pAn.x(),pAn.y());
  lEA->Draw();
}


void RPReconstructedTracksValidation::WriteHistograms(const std::string &root_file_name)
{
  TFile *f = TFile::Open(root_file_name.c_str(), "recreate");
  if(!f || !f->IsWritable())
  {
    std::cout<<"Output file not opened correctly!!"<<std::endl;
  }
  if(verbosity_)
    std::cout<<"Writting histograms..."<<std::endl;
  WriteBeamSmearingHistograms(f);
  WriteOverlapsHistograms(f);
  WriteStationFlux(f);
  rp_histograms_.Write(f);
  det_histograms_.Write(f);
//  station_histograms_.Write(f);
  if(verbosity_)
    std::cout<<"Writting histograms finnished."<<std::endl;
  f->Close();
}

DEFINE_FWK_MODULE(RPReconstructedTracksValidation);
