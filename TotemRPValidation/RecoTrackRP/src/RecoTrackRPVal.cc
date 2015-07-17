#include <iostream>
#include <string>
#include <boost/shared_ptr.hpp>

#include <TFile.h>
#include <TCanvas.h>
#include <TGraph.h>
#include <TLegend.h>
#include <TLine.h>
#include <TStyle.h>
#include <TH1.h>
#include <TH1D.h>
#include <TH2F.h>

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "HepMC/SimpleVector.h"
#include "HepMC/GenEvent.h"
#include "CLHEP/Vector/ThreeVector.h"
#include "CLHEP/Vector/Rotation.h"

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimG4Core/Notification/interface/SimG4Exception.h"

#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RP2DHit.h"
#include "TotemRPValidation/ParamMADRefTransport/interface/ParamMADRefTransport.h"
#include "TotemRPValidation/RecoTrackRP/interface/RecoTrackRPVal.h"


RecoTrackRPVal::RecoTrackRPVal(const edm::ParameterSet& conf)
 : rp_histograms_("/RP_Tracks/RP_", edm::ParameterSet()), 
   det_histograms_("/Dets/Det_", edm::ParameterSet()),
   resol_degrad_service_(conf),
   conf_(conf)
{
  rPFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
  clusterSetLabel = conf.getParameter<edm::InputTag>("ClusterSetLabel");
  hist_file_name_ = conf.getParameter<std::string>("HistogramFileName");
  OriginalHepMCProductLabel_ = conf.getParameter<std::string>("OriginalHepMCProductLabel");
  OriginalHepMCModuleName_ = conf.getParameter<std::string>("OriginalHepMCModuleName");
  SmearedHepMCProductLabel_ = conf.getParameter<std::string>("SmearedHepMCProductLabel");
  SmearedHepMCModuleName_ = conf.getParameter<std::string>("SmearedHepMCModuleName");
  verbosity_ = conf.getParameter<int>("Verbosity");
  validTracksOnly = conf.getUntrackedParameter<int>("ValidTracksOnly", 1);
  
  //init rp ids for the stations
  //order: from ip to outside, then top bottom
  //top first vertical, bottom first vertical, first horizontal, second horizontal, last top vertical, lat bottom vertical 
  for(int i=100; i<100+NPOTS; ++i)
    station_ids_[rp_150_r].push_back(i);
  for(int i=0; i<NPOTS; ++i)
    station_ids_[rp_150_l].push_back(i);
  for(int i=120; i<120+NPOTS; ++i)
    station_ids_[rp_220_r].push_back(i);
  for(int i=20; i<20+NPOTS; ++i)
    station_ids_[rp_220_l].push_back(i);
}


RecoTrackRPVal::~RecoTrackRPVal()
{
}


void RecoTrackRPVal::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  edm::ESHandle<BeamOpticsParams> BOParH;
  es.get<BeamOpticsParamsRcd>().get(BOParH);
  if(!BOParH.isValid())
    throw cms::Exception("RecoTrackRPVal") << " edm::ESHandle<BeamOpticsParams> is invalid";
  BOPar_ = *BOParH;
  parameterized_madx_transport_ = std::auto_ptr<ParamMADRefTransport>(new ParamMADRefTransport(conf_, es));
  InitBeamSmearingHistograms();
  InitOverlapsHistograms();
  InitHitHistograms();
  InitHitDistHistograms();
  InitAngleHistograms();
  
  edm::ESHandle<TotemRPGeometry> Totem_RP_geometry;
  es.get<RealGeometryRecord>().get(Totem_RP_geometry);
  if(!Totem_RP_geometry.isValid()) 
  {
    throw cms::Exception("RecoTrackRPVal") << "edm::ESHandle<TotemRPGeometry> is invalid, RP contours drawing impossible!";
  }
  
  Totem_RP_geometry_ = *Totem_RP_geometry;
}

void RecoTrackRPVal::InitOverlapsHistograms()
{
  char name[1024];
  
  int hist_bins = 160;
  
  sprintf(name, "overlaps_RP_120_121_122_");
  overlaps_RP_120_121_122_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_120_121_122_->SetDirectory(0);
  
  sprintf(name, "overlaps_RP_100_101_102_");
  overlaps_RP_100_101_102_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_100_101_102_->SetDirectory(0);
  
  sprintf(name, "overlaps_RP_020_021_022_");
  overlaps_RP_020_021_022_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_020_021_022_->SetDirectory(0);
  
  sprintf(name, "overlaps_RP_000_001_002_");
  overlaps_RP_000_001_002_ = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_RP_000_001_002_->SetDirectory(0);
  
  sprintf(name,
      "overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_");
  overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_120_121_122_tracks_when_RP_100_101_102_on_->SetDirectory(0);
  
  sprintf(name,
      "overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_");
  overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_100_101_102_tracks_when_RP_120_121_122_on_->SetDirectory(0);
  
  sprintf(name,
      "overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_");
  overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_020_021_022_tracks_when_RP_000_001_002_on_->SetDirectory(0);
  
  sprintf(name,
      "overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_");
  overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_
      = std::auto_ptr<TH2F>(new TH2F(name, name, hist_bins, -40, 40, hist_bins, -40, 40));
  overlaps_between_stations_RP_000_001_002_tracks_when_RP_020_021_022_on_->SetDirectory(0);
}


void RecoTrackRPVal::InitBeamSmearingHistograms()
{
  vertex_x_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("vertex_x_smearing_dist_", "vertex_x_smearing_dist_", 100,
          BOPar_.GetBeamDisplacementX()*1000, BOPar_.GetBeamDisplacementX()*1000+0.001));
  vertex_x_smearing_dist_->SetBit(TH1::kCanRebin);
  vertex_x_smearing_dist_->SetDirectory(0);

  vertex_y_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("vertex_y_smearing_dist_", "vertex_y_smearing_dist_", 100,
          BOPar_.GetBeamDisplacementY()*1000, BOPar_.GetBeamDisplacementY()*1000+0.001));
  vertex_y_smearing_dist_->SetBit(TH1::kCanRebin);
  vertex_y_smearing_dist_->SetDirectory(0);
  
  vertex_z_smearing_dist_ = std::auto_ptr<TH1D>(
      new TH1D("vertex_z_smearing_dist_", "vertex_z_smearing_dist_", 100,
          BOPar_.GetBeamDisplacementZ()*1000, BOPar_.GetBeamDisplacementZ()*1000+0.001));
  vertex_z_smearing_dist_->SetBit(TH1::kCanRebin);
  vertex_z_smearing_dist_->SetDirectory(0);
  
  for (int i=0; i<2; i++)
  {
	std::string name;
	name="px_";
	name+=(i?"right":"left");
	name+="smearing_dist_";
	px_smearing_dist_[i] = std::auto_ptr<TH1D>(
	  new TH1D(name.c_str(),name.c_str(), 100,
          BOPar_.GetCrossingAngleX()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetCrossingAngleX()*BOPar_.GetBeamMomentum()+0.001));
	px_smearing_dist_[i]->SetBit(TH1::kCanRebin);
	px_smearing_dist_[i]->SetDirectory(0);
  
	name="py_";
	name+=(i?"right":"left");
	name+="smearing_dist_";
	py_smearing_dist_[i] = std::auto_ptr<TH1D>(
	  new TH1D(name.c_str(),name.c_str(), 100,
          BOPar_.GetCrossingAngleY()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetCrossingAngleY()*BOPar_.GetBeamMomentum()+0.001));
	py_smearing_dist_[i]->SetBit(TH1::kCanRebin);
	py_smearing_dist_[i]->SetDirectory(0);
  
	name="pz_";
	name+=(i?"right":"left");
	name+="smearing_dist_";
	pz_smearing_dist_[i] = std::auto_ptr<TH1D>(
	  new TH1D(name.c_str(),name.c_str(), 100,
          BOPar_.GetMeanXi()*BOPar_.GetBeamMomentum(), 
          BOPar_.GetMeanXi()*BOPar_.GetBeamMomentum()+0.0001));
	pz_smearing_dist_[i]->SetBit(TH1::kCanRebin);
	pz_smearing_dist_[i]->SetDirectory(0);

	name="thx_";
	name+=(i?"right":"left");
	name+="smearing_dist_";
	thx_smearing_dist_[i] = std::auto_ptr<TH1D>(
	  new TH1D(name.c_str(),name.c_str(), 100,
          BOPar_.GetCrossingAngleX(), 
          BOPar_.GetCrossingAngleX()+1e-7));
	thx_smearing_dist_[i]->SetBit(TH1::kCanRebin);
	thx_smearing_dist_[i]->SetDirectory(0);
  
	name="thy_";
	name+=(i?"right":"left");
	name+="smearing_dist_";
	thy_smearing_dist_[i] = std::auto_ptr<TH1D>(
	  new TH1D(name.c_str(),name.c_str(), 100,
          BOPar_.GetCrossingAngleY(), 
          BOPar_.GetCrossingAngleY()+1e-7));
	thy_smearing_dist_[i]->SetBit(TH1::kCanRebin);
	thy_smearing_dist_[i]->SetDirectory(0);
  
	name="xi_";
	name+=(i?"right":"left");
	name+="smearing_dist_";
	xi_smearing_dist_[i] = std::auto_ptr<TH1D>(
	  new TH1D(name.c_str(),name.c_str(), 100,
          BOPar_.GetMeanXi(), 
          BOPar_.GetMeanXi()+1e-5));
	xi_smearing_dist_[i]->SetBit(TH1::kCanRebin);
	xi_smearing_dist_[i]->SetDirectory(0);
  }
}

void RecoTrackRPVal::InitHitHistograms()
{
	std::string name;
  int hist_bins = 80;
  
	for (int i=0; i<NSTATIONS; i++) {
		name="hits_"+StIdToName(i);
		hits_rp_[i] = std::auto_ptr<TH2F>(new TH2F(name.c_str(), name.c_str(), hist_bins, -40, 40, hist_bins, -40, 40));
		hits_rp_[i]->SetDirectory(0);
	}
}

void RecoTrackRPVal::InitHitDistHistograms()
{
	for (int s = 0; s < NSTATIONS; s++)
		for (int rp = 0; rp < NPOTS; rp++) {
			int i=station_ids_[s][rp];
			char buf[8];
			sprintf(buf, "det_%i", i);
			hitDists[i] = new TGraph();
			// set name for TGraph
			hitDists[i]->SetName(buf);
			hitDists[i]->SetTitle(buf);
		}
}

std::string RecoTrackRPVal::StIdToName(int st_id)
{
	std::string name;
	switch (st_id) {
		case rp_150_l: name="rp_150_l"; break;
		case rp_150_r: name="rp_150_r"; break;
		case rp_220_l: name="rp_220_l"; break;
		case rp_220_r: name="rp_220_r"; break;
	}
	return name;
}

void RecoTrackRPVal::InitAngleHistograms()
{
	  for (int st_id=0;st_id<NSTATIONS;st_id++) {
		  std::string name="angle_dist_ip_vs_rp_theta_x_"+StIdToName(st_id);
		  angle_dists_[0][st_id]=std::auto_ptr<TH2D>(new TH2D(name.c_str(), name.c_str(), 1000, -1.0e-3, 1.0e-3, 1000, -1.0e-3, 1.0e-3));
		  name="angle_dist_ip_vs_rp_theta_y_"+StIdToName(st_id);
		  angle_dists_[1][st_id]=std::auto_ptr<TH2D>(new TH2D(name.c_str(), name.c_str(), 1000, -1.0e-3, 1.0e-3, 1000, -1.0e-3, 1.0e-3));
		  name="1d_angle_dist_ip_theta_x_"+StIdToName(st_id);
		  angle_dists_1d_[0][0][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
		  name="1d_angle_dist_rp_theta_x_"+StIdToName(st_id);
		  angle_dists_1d_[1][0][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
		  name="1d_angle_dist_ip_theta_y_"+StIdToName(st_id);
		  angle_dists_1d_[0][1][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
		  name="1d_angle_dist_rp_theta_y_"+StIdToName(st_id);
		  angle_dists_1d_[1][1][st_id]=std::auto_ptr<TH1D>(new TH1D(name.c_str(), name.c_str(), 1000, -5.0e-3, 5.0e-3));
	  }
}

void RecoTrackRPVal::WriteBeamSmearingHistograms(TFile *f)
{
  f->cd("/");
  f->mkdir("beam_smearing_validation");
  f->cd("beam_smearing_validation");
  vertex_x_smearing_dist_->Write();
  vertex_y_smearing_dist_->Write();
  vertex_z_smearing_dist_->Write();
  
  for (int i=0; i<2; i++){
	  px_smearing_dist_[i]->Write();
	  py_smearing_dist_[i]->Write();
	  pz_smearing_dist_[i]->Write();

	  thx_smearing_dist_[i]->Write();
	  thy_smearing_dist_[i]->Write();

	  xi_smearing_dist_[i]->Write();
  }
  f->cd("/");
}

void RecoTrackRPVal::WriteOverlapsHistograms(TFile *f)
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
  f->mkdir("overlapping_proton_hits");
  f->cd("overlapping_proton_hits");
  
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

void RecoTrackRPVal::setupCanvas()
{
  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatBorderSize(1);
  gStyle->SetFuncWidth(1);
  gStyle->SetFuncColor(4);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.25);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetOptTitle(0);
  gStyle->ToggleEventStatus();
}

void RecoTrackRPVal::WriteHitDistHistograms(TFile *f)
{
	f->cd("/");
	try {
	f->mkdir("station_hits");
	} catch (...) {};
	f->cd("station_hits");
	// station_ids_[station] = <x00, x01, ..., x05>
	// unitIDs = 000, 003, 020, 023, 100, 103, 120, 123
	for (int st_id=0; st_id<NSTATIONS; st_id++){
	    char buf[15];
	    for(unsigned int unitID=station_ids_[st_id][0]; unitID<station_ids_[st_id][0]+6; unitID+=3) {
	    	sprintf(buf, "hits_unit_%i", unitID);
	        // create one TCanvas per unitID
	        TCanvas * c1 = new TCanvas(buf,buf,200,10,700,500);
	        c1->cd();
	        TH2D *frame = new TH2D(buf, ";x   (mm);y   (mm)", 100, -70., +70., 100, -70., +70.);
	        frame->Draw();

	        setupCanvas();
	        TLegend *l = new TLegend(0.1, 0.1, 0.2, 0.2);

	        int color = 0;
	        int totalNumberOfHits = 0;

	        for (int i=0; i<3; i++) {
				totalNumberOfHits += hitDists.count( unitID+i );

				// draw RP envelopes for all 3 RPs
				DrawRPContour(unitID+i);

				// draw hits in first RP (if present)
				if( hitDists.count( unitID+i ) > 0 ){
					TGraph *g = hitDists[unitID+i];
					if (g->GetN() != 0){
					g->SetMarkerColor(++color);
					g->SetMarkerStyle(20);
					g->SetMarkerSize(0.3);
					g->Draw("P");
					l->AddEntry(g, g->GetName(), "p");
					}
				}
	        }

	        if( totalNumberOfHits > 0 )
	        	l->Draw("same");
			c1->Modified();
			c1->Update();
			c1->Write("");
			// delete l; delete frame; delete c1;
	    }
	}
    f->cd("/");
}

void RecoTrackRPVal::WriteHitHistorgrams(TFile *f)
{
  f->cd("/");
  try {
  f->mkdir("station_hits");
  } catch (...) {};
  f->cd("station_hits");
  
  for (int i=0; i<NSTATIONS; i++) {
		hits_rp_[i]->Write();
		WriteHistogramAsCanvas(&*(hits_rp_[i]), station_ids_[i]);
  }
  
  f->cd("/");
}


void RecoTrackRPVal::WriteHistogramAsCanvas(TH2F *h, const std::vector<RPId> &rpids)
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
  
  // we draw the contours of only half RPs from a station
  for(unsigned int i=0; i<NPOTS/2; ++i)
  {
    DrawRPContour(rpids[i]);
  }

  c1->Modified();
  c1->Update();
  c1->Write("");
}

void RecoTrackRPVal::WriteAngleHistograms(TFile *f){
	f->cd("/");
	f->mkdir("angle");
	f->cd("angle");
	for (int st_id=0;st_id<NSTATIONS;st_id++) {
		for (int arm=0; arm<2; arm++) {
			angle_dists_[arm][st_id]->Write();
			angle_dists_1d_[0][arm][st_id]->Write();
			angle_dists_1d_[1][arm][st_id]->Write();
		}
	}
	f->cd("/");
}

// calculates tangent of the angle of the fitted track in the ZX or ZY plane
int RecoTrackRPVal::calcTrackTheta(std::vector<double> &theta, const rec_tracks_collection &coll, int arm)
{
	rec_tracks_collection::const_iterator it;
	rec_tracks_collection vert, horiz;
	RP2DHit th1, th2;
	double dz;
	if (verbosity_) {
		std::cout<<"RecoTrackRP::calcTrackTheta: arm "<< (arm?"right":"left")
		         << ", got "<<coll.size() << " tracks, RPs ";
		for (rec_tracks_collection::const_iterator i=coll.begin(); i!=coll.end(); i++) {
			std::cout << " " << i->first;
		}
		std::cout << std::endl;
	}
	if (coll.size()<2)
		return 0;

	for (it=coll.begin(); it!=coll.end(); it++) {
		int pot_id = it->first % 20;
		switch (pot_id) {
			case near_top:
			case near_bottom:
			case far_top:
			case far_bottom: vert[it->first]=it->second; break;
			case near_horiz:
			case far_horiz: horiz[it->first]=it->second; break;
		}
	}
	if (verbosity_) std::cout << "RecoTrackRP::calcTrackTheta: "
	                          << vert.size() << " vertical, "
	                          << horiz.size() << " horizontal pots... ";
	if (vert.size()>=2) {
		it=vert.begin();
		if (verbosity_) std::cout << "using vertical." << std::endl;
	} else if (horiz.size()>=2) {
		it=horiz.begin();
		if (verbosity_) std::cout << "using horizontal." << std::endl;
	} else
		return 0;
	th1=it->second;
	it++;
	th2=it->second;
	if (verbosity_) std::cout<<"RecoTrackRP::calcTrackTheta: hit1: "
	                         << th1.X() << " " << th1.Y() << " " << th1.Z() << " hit2: "
	                         << th2.X() << " " << th2.Y() << " " << th2.Z() << std::endl;
	dz=th2.Z()-th1.Z();
	theta.reserve(2);
	theta[0]=(th2.X()-th1.X())/dz;
	theta[1]=(th2.Y()-th1.Y())/dz;
	if (verbosity_) std::cout << "RecoTrackRP::calcTrackTheta: thetaX: " << theta[0]
	                          << " thetaY: " << theta[1] << std::endl;
	return 2;
}

void RecoTrackRPVal::analyze(const edm::Event& e, const edm::EventSetup& es)
{
  //fill track and cluster related histograms 
  edm::Handle< RPFittedTrackCollection > input;
  e.getByLabel(rPFittedTrackCollectionLabel, input);
  if (!input.isValid()) 
    throw cms::Exception("RecoTrackRPVal") << "edm::Handle< RPFittedTrackCollection > is invalid";
  
  edm::Handle< cluster_set > det_clusters;
  e.getByLabel(clusterSetLabel, det_clusters);
  if (!det_clusters.isValid()) 
    throw cms::Exception("RecoTrackRPVal") << "edm::Handle<edm::det_clusters> det_clusters is invalid";
  
  FillTrackPositions(input);
  FillDetectorHistogramsForReconstructedTracks(*input, *det_clusters);
  FillOverlappingAreas(*input);
  FillHits(*input);
  FillHitDists(*input);

  //fill acceptance and smearing information if available
  edm::Handle<edm::HepMCProduct> OriginalHepMCEvt, SmearedHepMCEvt;
  e.getByLabel(OriginalHepMCModuleName_, OriginalHepMCProductLabel_, OriginalHepMCEvt ); 
  e.getByLabel(SmearedHepMCModuleName_, SmearedHepMCProductLabel_, SmearedHepMCEvt );
  if( !(OriginalHepMCEvt.isValid() && FindProtons(OriginalHepMCEvt,0)) )
  {
    if (verbosity_ )
      std::cout<<"Primary protons not found!! Skipping the beam related event validation."<<std::endl;
    return;
  }
  
  if ( !(SmearedHepMCEvt.isValid() && FindProtons(SmearedHepMCEvt,1)) )
  {
    if (verbosity_ )
      std::cout<<"Smeared protons not found!! Skipping the beam related event validation."<<std::endl;
    return;
  }    
  
  FillEventSmearingInformation();

  rec_tracks_collection hits[NSTATIONS];
  int no_of_pots[NSTATIONS];

  for (int i=0; i<NSTATIONS; i++)
  {
	  if (prot_[i/2][0].found)
		  for (size_t j=0; j<station_ids_[i].size(); j++)
			  FillReferenceHistograms(station_ids_[i][j], prot_[i/2][0]);
	  no_of_pots[i] = SelectHits(*input, hits[i], station_ids_[i]);
	  if (no_of_pots[i])
		  FillResidualHistograms(hits[i], prot_[i/2][0], prot_[i/2][1], no_of_pots[i%2]==2);
	  if (prot_[i/2][0].found && prot_[i/2][1].found && no_of_pots[i]>1)
			FillAngleHistograms(hits[i], i);
  }
}

void RecoTrackRPVal::FillAngleHistograms(const rec_tracks_collection &hits, int sid)
{
	std::vector<double> trackTheta;
	if (calcTrackTheta(trackTheta, hits, sid/2)<2) {
		if (verbosity_) std::cout << "RecoTrackRP::analyze: calcTrackTheta returned no value." << std::endl;
		return;
	}
	double xi=BOPar_.IPSmearedProtonMomentumToXi(prot_[sid%2][1].momentum.px(),prot_[sid%2][1].momentum.py(),prot_[sid%2][1].momentum.pz());
	angle_dists_[0][sid]->Fill(prot_[sid/2][1].thetaX,trackTheta[0]);
	angle_dists_[1][sid]->Fill(prot_[sid/2][1].thetaY,trackTheta[1]);
	if (verbosity_) {
		std::cout<<"RecoTrackRP::analyze: thetaX(IP)   thetaX(RP)   thetaY(IP)  thetaY(RP)   xi" << std::endl;
	    std::cout<<"                      " << prot_[sid%2][1].thetaX << " " << trackTheta[0] << " " <<
	                                           prot_[sid%2][1].thetaY << " " << trackTheta[1] << " " <<
	                                           xi << std::endl;
	}
	angle_dists_1d_[0][0][sid]->Fill(prot_[sid/2][1].thetaX);
	angle_dists_1d_[0][1][sid]->Fill(prot_[sid/2][1].thetaY);
	angle_dists_1d_[1][0][sid]->Fill(trackTheta[0]);
	angle_dists_1d_[1][1][sid]->Fill(trackTheta[1]);
}

void RecoTrackRPVal::FillOverlappingAreas(const RPFittedTrackCollection &tracks)
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


void RecoTrackRPVal::FillTrackPositions(edm::Handle< RPFittedTrackCollection > &input)
{
  RPFittedTrackCollection::const_iterator it;
  for(it=(*input).begin(); it!=(*input).end(); ++it)
  {
    rp_histograms_.GetObj(it->first)->FillBasicTrackInfo(it->second);
  }
}


void RecoTrackRPVal::FillEventSmearingInformation()
{
	bool vertex_info_filled = false;
	for (int arm=0; arm<2; arm++)
		if (prot_[arm][0].found && prot_[arm][1].found)
		{
			if(!vertex_info_filled)
			{
				vertex_x_smearing_dist_->Fill(prot_[arm][1].vertex.x()-prot_[arm][0].vertex.x());
				vertex_y_smearing_dist_->Fill(prot_[arm][1].vertex.y()-prot_[arm][0].vertex.y());
				vertex_z_smearing_dist_->Fill(prot_[arm][1].vertex.z()-prot_[arm][0].vertex.z());
				vertex_info_filled = true;
			}

		    px_smearing_dist_[arm]->Fill(prot_[arm][1].momentum.px()-prot_[arm][0].momentum.px());
		    py_smearing_dist_[arm]->Fill(prot_[arm][1].momentum.py()-prot_[arm][0].momentum.py());
		    pz_smearing_dist_[arm]->Fill(prot_[arm][1].momentum.pz()-prot_[arm][0].momentum.pz());

		    thx_smearing_dist_[arm]->Fill(
		        prot_[arm][1].momentum.px()/TMath::Abs(prot_[arm][1].momentum.pz())-
		        prot_[arm][0].momentum.px()/TMath::Abs(prot_[arm][0].momentum.pz()));

		    thy_smearing_dist_[arm]->Fill(
		        prot_[arm][1].momentum.py()/TMath::Abs(prot_[arm][1].momentum.pz())-
		        prot_[arm][0].momentum.py()/TMath::Abs(prot_[arm][0].momentum.pz()));
		    xi_smearing_dist_[arm]->Fill(prot_[arm][1].momentum.rho()/BOPar_.GetBeamMomentum()-
		        prot_[arm][0].momentum.rho()/BOPar_.GetBeamMomentum());
		}
}


void RecoTrackRPVal::FillHits(
    const RPFittedTrackCollection &tracks)
{
  RPFittedTrackCollection::const_iterator it1;
  
  for (int sid=0; sid<NSTATIONS; sid++)
	  for (int i=0; i<NPOTS/2; i++) {
		  it1 = tracks.find(station_ids_[sid][i]);
		  if( it1 != tracks.end() ) {
			  hits_rp_[sid]->Fill(it1->second.X0(), it1->second.Y0());
			  break;
		  }

	  }
}

void RecoTrackRPVal::FillHitDists(const RPFittedTrackCollection &tracks)
{
	for (RPFittedTrackCollection::const_iterator it = tracks.begin(); it != tracks.end(); ++it) {
		// check if track is valid
		if( validTracksOnly )
			if (!it->second.IsValid()) continue;
	    // check if graph corresponding to RPId exists in hitDists map
	    if (hitDists.find(it->first) != hitDists.end()) {
	    	// add point to corresponding TGraph
	    	hitDists[it->first]->SetPoint(hitDists[it->first]->GetN(), it->second.X0(), it->second.Y0());
	    	if(verbosity_)
				std::cout << "Hit in detector " << it->first << "  : ( " << it->second.X0() << " , " << it->second.Y0() << " ) , IsValid()=" << it->second.IsValid() << std::endl;
	    } else {
	    	if(verbosity_)
	    		std::cout << "RPId " << it->first << " not specified in cfg file" << std::endl;
	    }
	}
}

void RecoTrackRPVal::FillDetectorHistogramsForReconstructedTracks(
    const RPFittedTrackCollection &tracks, const cluster_set &clusters)
{
  for (RPFittedTrackCollection::const_iterator trit = tracks.begin(); trit
      !=tracks.end(); ++trit)
  {
    FillDetectorHistograms(trit->second, clusters);
  }
}


void RecoTrackRPVal::FillReferenceHistograms(RPId id, const PrimaryProton &proton)
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
    std::cout<<"RecoTrackRPVal::FillReferenceHistograms"<<std::endl;
    std::cout<<"id="<<id<<" px="<<px<<" py="<<py<<" pz="<<pz<<std::endl;
    std::cout<<"xi="<<xi<<" t="<<t<<std::endl;
  }
}


//return number of pots through which the proton goes
int RecoTrackRPVal::SelectHits(const RPFittedTrackCollection &tracks, 
    rec_tracks_collection & coll, const station_rp_ids_type& st_ids)
{
  RPFittedTrackCollection::const_iterator it;
  RPFittedTrackCollection::const_iterator low_it;
  RPFittedTrackCollection::const_iterator high_it;
  
  assert(st_ids.size() == NPOTS);
  
  low_it = tracks.lower_bound(st_ids[near_top]);
  high_it = tracks.upper_bound(st_ids[far_bottom]);
  
  for(it = low_it; it!=high_it; ++it)
  {
    if(!it->second.IsValid())
      continue;
    coll[it->first]=resol_degrad_service_.Create2DHit(it->first, it->second);
  }
  
  rec_tracks_collection::const_iterator end;
  end = coll.end();
  
  bool station_accepted = 
    (coll.find(st_ids[near_top])!=end || coll.find(st_ids[near_bottom])!=end || coll.find(st_ids[near_horiz])!=end) &&
    (coll.find(st_ids[far_horiz])!=end || coll.find(st_ids[far_top])!=end || coll.find(st_ids[far_bottom])!=end) &&
    !(coll.find(st_ids[near_top])!=end && coll.find(st_ids[near_bottom])!=end) &&
    !(coll.find(st_ids[far_top])!=end && coll.find(st_ids[far_bottom])!=end);
  
  if(!station_accepted)
    coll.clear();
  
  return coll.size();
}



void RecoTrackRPVal::FillResidualHistograms(
    const rec_tracks_collection &tracks, const PrimaryProton &prim_proton, 
    const PrimaryProton &smeared_proton, bool proton_at_150_station)
{
  if(!prim_proton.found || !smeared_proton.found)
    return;
  
  TVector2 rp_transverse_track_position;
  
  for(rec_tracks_collection::const_iterator it = tracks.begin(); 
        it!=tracks.end(); ++it)
  {
    double spx = smeared_proton.momentum.x();
    double spy = smeared_proton.momentum.y();
    double spz = smeared_proton.momentum.z();
    
    double svx = smeared_proton.vertex.x();//+BOPar_.GetBeamDisplacementX()*1000.0;
    double svy = smeared_proton.vertex.y();//+BOPar_.GetBeamDisplacementY()*1000.0;
    double svz = smeared_proton.vertex.z();//+BOPar_.GetBeamDisplacementZ()*1000.0;
    
    bool res = parameterized_madx_transport_->MADXTeoreticalRPTransversePosition(it->first, 
          svx, svy, svz, spx, spy, spz, rp_transverse_track_position, false);
    
    double ovx = prim_proton.vertex.x();//+BOPar_.GetBeamDisplacementX()*1000.0;
    double ovy = prim_proton.vertex.y();//+BOPar_.GetBeamDisplacementY()*1000.0;
    double ovz = prim_proton.vertex.z();//+BOPar_.GetBeamDisplacementZ()*1000.0;
    
    double opx = prim_proton.momentum.x();
    double opy = prim_proton.momentum.y();
    double opz = prim_proton.momentum.z();
    
    double xi = BOPar_.IPNonSmearedProtonMomentumToXi(opx, opy, opz);
    double t = BOPar_.IPNonSmearedProtonMomentumTot(opx, opy, opz);
//    double t1=-opz*opz*((opx/opz)*(opx/opz) + (opy/opz)*(opy/opz));
    
    boost::shared_ptr<RPTrackInfo> ptr = rp_histograms_.GetObj(it->first);
    ptr->FillTrackSeen(t);
    ptr->FillTrackSeen(t, xi);
    
    if(verbosity_>0)
    {
      std::cout<<"RecoTrackRPVal::FillResidualHistograms"<<std::endl;
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


void RecoTrackRPVal::FillDetectorHistograms(
        const RPFittedTrack &track, const cluster_set &clusters)
{
  for(int i = 0; i<track.GetHitEntries(); ++i)
  {
    RPDetId sim_det_id = track.GetHit(i).DetId();
    double residual = track.GetHit(i).Residual();
    TotRPDetId det_raw_id(sim_det_id);
    RPDetId dec_det_id = det_raw_id.DetectorDecId();
    cluster_set::const_iterator it = clusters.find(sim_det_id);
    if(it!=clusters.end())
    {
      det_histograms_.GetObj(dec_det_id)->FillHistograms(*it, residual);
    }
  }
}


bool RecoTrackRPVal::FindProtons(const edm::Handle<edm::HepMCProduct> &HepMCEvt, int smeared)
{
  if(!HepMCEvt.isValid())
  {
     throw SimG4Exception("RecoTrackRPVal : Unable to find HepMCProduct(HepMC::GenEvent) in edm::Event  ");
  }
  
  const HepMC::GenEvent *evt = HepMCEvt->GetEvent();
  
  int count[2] = { 0, 0 };
  prot_[0][smeared].found = prot_[1][smeared].found = false;
  
  for(HepMC::GenEvent::particle_const_iterator it = evt->particles_begin(); 
      it != evt->particles_end(); ++it )
  {
    HepMC::GenParticle * g = (*it);
    int g_status = g->status();
    int pdg_id = g->pdg_id();
    
    if(verbosity_)
    {
      g->print(std::cout);
      //std::cout<<g->production_vertex()->position()<<std::endl;
    }
    
    // scanning only for undecayed particles, status == 1, see http://cepa.fnal.gov/psm/simulation/mcgen/lund/pythia_manual/pythia6.3/pythia6301/node39.html
    // protons have pdg_id == 2212, see http://www-pdg.lbl.gov/2009/reviews/rpp2009-rev-monte-carlo-numbering.pdf
    if (g_status == 1 && pdg_id == 2212)
    {
      const HepMC::FourVector &vtx = g->production_vertex()->position();
      const HepMC::FourVector &mom  = g->momentum();
      if(verbosity_ > 1) {
        std::cout<<"Setting the " << (smeared?"smeared":"primary") << " protons."<<std::endl;
      }
      if(mom.z())
      {
        int arm=(mom.z()>0.0);
        ++count[arm];
        prot_[arm][smeared].vertex = vtx; //[mm]
        prot_[arm][smeared].momentum = mom;  //[GeV]
        prot_[arm][smeared].found = true;
        prot_[arm][smeared].thetaX = BOPar_.ComputeCrossingAngleCorrectedThetaX(mom);
        prot_[arm][smeared].thetaY = BOPar_.ComputeCrossingAngleCorrectedThetaY(mom);
      }
    }
  } // end loop on HepMC particles

  bool result = count[1]<=1 && count[0]<=1 && (count[1]>0 || count[0]>0);

  if(verbosity_ > 1)
    std::cout<<"right_count="<<count[1]<<" left_count="<<count[0]
        <<" result="<<result<<std::endl;
  return result;
}

void RecoTrackRPVal::endJob()
{
	TFile *f = TFile::Open(hist_file_name_.c_str(), "recreate");
	if(!f || !f->IsWritable())
		std::cout<<"Output file not opened correctly!!"<<std::endl;
	if(verbosity_)
		std::cout<<"Writing histograms..."<<std::endl;
	WriteBeamSmearingHistograms(f);
	WriteOverlapsHistograms(f);
	WriteHitHistorgrams(f);
	WriteHitDistHistograms(f);
	WriteAngleHistograms(f);
	rp_histograms_.Write(f);
	det_histograms_.Write(f);
	if(verbosity_)
		std::cout<<"Writing histograms finished."<<std::endl;
	f->Close();
}

void RecoTrackRPVal::DrawRPContour(RPId dec_rpid)
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

DEFINE_FWK_MODULE(RecoTrackRPVal);
