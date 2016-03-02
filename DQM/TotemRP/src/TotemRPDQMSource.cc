#include "DQM/TotemRP/interface/TotemRPDQMSource.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <iostream>
#include <iomanip>
#include <stdio.h>
#include <string>
#include <sstream>
#include <math.h>

using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::DiagonalPlots::DiagonalPlots(DQMStore::IBooker &ibooker, int _id) : id(_id)
{
  bool top45 = id & 2;
  bool top56 = id & 1;
  bool diag = (top45 != top56);

  char name[50];
  sprintf(name, "%s 45%s - 56%s",
    (diag) ? "diagonal" : "antidiagonal",
    (top45) ? "top" : "bot",
    (top56) ? "top" : "bot"
  );

  ibooker.setCurrentFolder(string("TotemRP/") + name);

  h_lrc_x_d = ibooker.book2D("dx left vs right", string(name) + " : dx left vs. right, histogram;#Delta x_{45};#Delta x_{56}", 50, 0., 0., 50, 0., 0.);
  h_lrc_x_n = ibooker.book2D("xn left vs right", string(name) + " : xn left vs. right, histogram;x^{N}_{45};x^{N}_{56}", 50, 0., 0., 50, 0., 0.);
  h_lrc_x_f = ibooker.book2D("xf left vs right", string(name) + " : xf left vs. right, histogram;x^{F}_{45};x^{F}_{56}", 50, 0., 0., 50, 0., 0.);

  h_lrc_y_d = ibooker.book2D("dy left vs right", string(name) + " : dy left vs. right, histogram;#Delta y_{45};#Delta y_{56}", 50, 0., 0., 50, 0., 0.);
  h_lrc_y_n = ibooker.book2D("yn left vs right", string(name) + " : yn left vs. right, histogram;y^{N}_{45};y^{N}_{56}", 50, 0., 0., 50, 0., 0.);
  h_lrc_y_f = ibooker.book2D("yf left vs right", string(name) + " : yf left vs. right, histogram;y^{F}_{45};y^{F}_{56}", 50, 0., 0., 50, 0., 0.);
}

//----------------------------------------------------------------------------------------------------

#if 0

TotemRPDQMSource::ArmPlots::ArmPlots(int _id) : id(_id)
{
  MultiRootPlot *pl_numRPWithTrack = new MultiRootPlot();

  h_numRPWithTrack_top = new TH1D("", "", 5, -0.5, 4.5);
  h_numRPWithTrack_top->SetLineColor(1);
  pl_numRPWithTrack->Add(h_numRPWithTrack_top, "");

  h_numRPWithTrack_hor = new TH1D("", "", 5, -0.5, 4.5);
  h_numRPWithTrack_hor->SetLineColor(2);
  pl_numRPWithTrack->Add(h_numRPWithTrack_hor, "same");

  h_numRPWithTrack_bot = new TH1D("", "", 5, -0.5, 4.5);
  h_numRPWithTrack_bot->SetLineColor(4);
  pl_numRPWithTrack->Add(h_numRPWithTrack_bot, "same");

  PlotManager::RegisterRP(TotRPDetId::lArm, id, "number of RPs with tracks in parallel groups",
    "number of RPs with tracks in parallel groups;number of RPs with tracks: top [black], horiz [red], bottom [blue]", pl_numRPWithTrack);

  h_trackCorr = new TH2D("", "", 13, -0.5, 12.5, 13, -0.5, 12.5);
  TAxis *xa = h_trackCorr->GetXaxis(), *ya = h_trackCorr->GetYaxis();
  xa->SetBinLabel(1, "210, near, top"); ya->SetBinLabel(1, "210, near, top");
  xa->SetBinLabel(2, "bot"); ya->SetBinLabel(2, "bot");
  xa->SetBinLabel(3, "hor"); ya->SetBinLabel(3, "hor");
  xa->SetBinLabel(4, "far, hor"); ya->SetBinLabel(4, "far, hor");
  xa->SetBinLabel(5, "top"); ya->SetBinLabel(5, "top");
  xa->SetBinLabel(6, "bot"); ya->SetBinLabel(6, "bot");
  xa->SetBinLabel(8, "220, near, top"); ya->SetBinLabel(8, "220, near, top");
  xa->SetBinLabel(9, "bot"); ya->SetBinLabel(9, "bot");
  xa->SetBinLabel(10, "hor"); ya->SetBinLabel(10, "hor");
  xa->SetBinLabel(11, "far, hor"); ya->SetBinLabel(11, "far, hor");
  xa->SetBinLabel(12, "top"); ya->SetBinLabel(12, "top");
  xa->SetBinLabel(13, "bot"); ya->SetBinLabel(13, "bot");
  PlotManager::RegisterRP(TotRPDetId::lArm, id, "track RP correlation", "track RP correlation", h_trackCorr, "colz,text");

  h_trackCorr_overlap = new TH2D(*h_trackCorr);
  PlotManager::RegisterRP(TotRPDetId::lArm, id, "track RP correlation, hor-vert overlaps", "track RP correlation, hor-vert overlaps",
    h_trackCorr_overlap, "colz,text");
}

//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::StationPlots::StationPlots(int _id, std::set<unsigned int> planes, 
  bool allocateCorrelationPlots,
  CorrelationPlotsSelector *correlationPlotsSelector, int limit) : 
    id(_id), correlationPlotsSelector(correlationPlotsSelector)
{
  rpHits = new TGraph();
  rpHits->SetMarkerStyle(2);
  rpHits->SetMarkerSize(1.0);

  PlotManager::RegisterRP(TotRPDetId::lStation, id, "RP fits in X,Y", "RP fits in X,Y;x   (mm);y   (mm)", rpHits);

  if (allocateCorrelationPlots)
    Add(planes, limit);
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::StationPlots::Add(std::set<unsigned int> planes, int limit)
{
  int correlationPlotsCounter = 0;
  bool limited;
  if (limit == -1)
    limited = false;
  else
    limited = true;
  for (std::set<unsigned int>::iterator i = planes.begin(); i != planes.end(); i++)
  {
    for (std::set<unsigned int>::iterator j = i; j != planes.end(); j++)
    {
      if (hist[*i][*j] == 0 && hist[*j][*i] == 0 && *i != *j && i != j)
      {
        unsigned int plane1 = std::min(*i, *j);
        unsigned int plane2 = std::max(*i, *j);
        if (correlationPlotsSelector->IfTwoCorrelate(plane1, plane2) && (!limited || correlationPlotsCounter < limit))
        {
          char buf1[200];
          char buf2[200];
          unsigned int RPPlaneId1 = plane1 + 100 * id;
          unsigned int RPPlaneId2 = plane2 + 100 * id;
          std::string RPPlane1 = "";
          std::string RPPlane2 = "";
          RPPlane1 += TotRPDetId::PlaneName(RPPlaneId1, TotRPDetId::nPath);
          RPPlane2 += TotRPDetId::PlaneName(RPPlaneId2, TotRPDetId::nPath);
          size_t pos1 = RPPlane1.rfind('/');
          size_t pos2 = RPPlane2.rfind('/');
          
          if (pos1 == std::string::npos)
          {
            pos1 = 0;
          } else {
            RPPlane1 = RPPlane1.substr(0, pos1);
          }
          
          if (pos2 == std::string::npos) {
            pos2 = 0;
          } else {
            RPPlane2 = RPPlane2.substr(0, pos2);
          }

          pos1 = RPPlane1.rfind('/');
          pos2 = RPPlane2.rfind('/');

          if (pos1 == std::string::npos) {
            pos1 = 0;
          } else {
            RPPlane1 = RPPlane1.substr(pos1 + 1);
          }

          if (pos2 == std::string::npos)
          {
            pos2 = 0;
          } else {
            RPPlane2 = RPPlane2.substr(pos2 + 1);
          }
          RPPlane1 += '_';
          RPPlane2 += '_';
          sprintf(buf1, "%s%u", RPPlane1.c_str(), plane1 % 10);
          sprintf(buf2, "%s%u", RPPlane2.c_str(), plane2 % 10);
          Int_t bins[2] = {512, 512};
          Double_t xmin[2] = {-0.5, -0.5};
          Double_t xmax[2] = {511.5, 511.5};
          hist[plane1][plane2] = hist[plane2][plane1] = new THnSparseD("name", "title", 2, bins, xmin, xmax);
          correlationPlotsCounter++;
          PlotManager::RegisterRP(TotRPDetId::lStation, id, std::string() + "correlation profile  " +  buf1 + " vs " + buf2,
            std::string("correlation profile ") + buf1 + " vs " + buf2 + ";strip in " + buf1 + ";strip in " + buf2, hist[plane1][plane2]);
        }
      }
    }
  }

  if (limited && correlationPlotsCounter >= limit)
    printf("WARNING in TotemRPDQMSource > Number of correlation plots for station %i has been limited to %i.\n", id, limit);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::PotPlots::PotPlots() :
  activity(NULL), activity_u(NULL), activity_v(NULL),
  hit_plane_hist(NULL),
  patterns_u(NULL), patterns_v(NULL),
  event_category(NULL),
  uHitsAll(NULL), vHitsAll(NULL), uHitsSel(NULL), 
  vHitsSel(NULL), uTrack(NULL), vTrack(NULL), uView(NULL), vView(NULL), uvGlobalView(NULL)
{
}

//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::PotPlots::PotPlots(unsigned int id, const TotemRPGeometryLite *geometry)
{
  MultiRootPlot *plane_activity = new MultiRootPlot();
  activity = new TH1D("", "", 11, -0.5, 10.5); activity->SetLineColor(1);
  plane_activity->Add(activity, "");
  activity_u = new TH1D("", "", 11, -0.5, 10.5);activity_u->SetLineColor(2);
  plane_activity->Add(activity_u, "same");
  activity_v = new TH1D("", "", 11, -0.5, 10.5); activity_v->SetLineColor(4);
  plane_activity->Add(activity_v, "same");
  PlotManager::RegisterRP(TotRPDetId::lRP, id,
    "active planes histogram", "active planes;number of active planes [black], U planes only [red], V planes only [blue]", plane_activity);


  hit_plane_hist = new TH2D("", "", 10, -0.5, 9.5, 32, -0.5, 511.5);
  hit_plane_hist->SetStats(false);
  PlotManager::RegisterRP(TotRPDetId::lRP, id,
    "activity in planes (2D)", "activity in planes;plane number;strip number", hit_plane_hist);

  
  MultiRootPlot *patterns_hists = new MultiRootPlot();
  patterns_u = new TH1D("", "", 11, -0.5, 10.5); patterns_u->SetLineColor(2);
  patterns_hists->Add(patterns_u, "");
  patterns_v = new TH1D("", "", 11, -0.5, 10.5); patterns_v->SetLineColor(4);
  patterns_hists->Add(patterns_v, "same");
  PlotManager::RegisterRP(TotRPDetId::lRP, id,
    "numbers of recognized patterns", "recognized patterns;number of recognized patterns in U [red] or V [blue]", patterns_hists);


  MultiRootPlot *pl_planes_in_fit = new MultiRootPlot();
  h_planes_fit_u = new TH1D("", "", 6, -0.5, 5.5); h_planes_fit_u->SetLineColor(2);
  pl_planes_in_fit->Add(h_planes_fit_u, "");
  h_planes_fit_v = new TH1D("", "", 6, -0.5, 5.5); h_planes_fit_v->SetLineColor(4);
  pl_planes_in_fit->Add(h_planes_fit_v, "same");

  PlotManager::RegisterRP(TotRPDetId::lRP, id,
    "planes contributing to fit", "planes contributing to fit;number of planes contributing to fit, U [red], V [blue]", pl_planes_in_fit);


  event_category = new TH1D("", "", 5, -0.5, 4.5);
  event_category->GetXaxis()->SetBinLabel(1, "empty");
  event_category->GetXaxis()->SetBinLabel(2, "insufficient");
  event_category->GetXaxis()->SetBinLabel(3, "single-track");
  event_category->GetXaxis()->SetBinLabel(4, "multi-track");
  event_category->GetXaxis()->SetBinLabel(5, "shower");
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "event category", "event category;", event_category);


  trackHitsCumulative = new TGraph();
  PlotManager::RegisterRP(TotRPDetId::lRP, id,
    "track XY profile (graph)", "track XY profile (graph);x   (mm);y   (mm)", trackHitsCumulative);
  trackHitsCumulativeHist = new TH2D("", "", 100, -18., +18., 100, -18., +18.);
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "track XY profile (histogram)",
    "track XY profile (histogram);x   (mm);y   (mm)", trackHitsCumulativeHist);

  MultiRootPlot *track_profiles = new MultiRootPlot();
  track_u_profile = new TH1D("", "", 512, -256*66E-3, +256*66E-3); track_u_profile->SetLineColor(2);
  track_profiles->Add(track_u_profile, "");
  track_v_profile = new TH1D("", "", 512, -256*66E-3, +256*66E-3); track_v_profile->SetLineColor(4);
  track_profiles->Add(track_v_profile, "same");
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "track u and v profiles",
    "track u and v profiles;u [red] or v [blue]   (mm)", track_profiles);

  PlotManager::RegisterRP(TotRPDetId::lRP, id, "current track 3D", "current track 3D",
    currentTrackInRP = new TrackInRPFigure(id, geometry, "current track 3D"));
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "all tracks in 3D", "all tracks in 3D",
    allTracksInRP = new TrackInRPFigure(id, geometry, "all tracks in 3D"));

  // detector shape in U-V
  double C1 = RPTopology::x_width_ / 2.;
  double C2 = C1 - RPTopology::phys_edge_lenght_ / sqrt(2.);

  double ds_u[6] = {+C1, +C1, -C1, -C1, -C2, +C1};
  double ds_v[6] = {-C1, +C1, +C1, -C2, -C1, -C1};

  // detector shape in X-Y
  CLHEP::Hep3Vector rp_pos = geometry->GetRPPosition(id);
  CLHEP::Hep3Vector rod_U = geometry->GetRPMeanUDirection(id);
  CLHEP::Hep3Vector rod_V = geometry->GetRPMeanVDirection(id);

  double ds_x[6], ds_y[6];
  for (unsigned int i = 0; i < 6; i++)
  {
    ds_x[i] = rod_U.x() * ds_u[i] + rod_V.x() * ds_v[i] + rp_pos.x();
    ds_y[i] = rod_U.y() * ds_u[i] + rod_V.y() * ds_v[i] + rp_pos.y();
  }

  detectorShape = new TGraph(6, ds_x, ds_y);

  currentTrackXY = new TGraph();
  currentTrackXY->SetMarkerStyle(20);
  currentTrackXY->SetMarkerSize(1.);
  
  currentMultiTracksXY = new TGraph();
  currentMultiTracksXY->SetMarkerStyle(24);
  currentMultiTracksXY->SetMarkerSize(1.);

  MultiRootPlot *view2D = new MultiRootPlot();
  view2D->Add(detectorShape, "AL");
  view2D->Add(currentMultiTracksXY, "P");
  view2D->Add(currentTrackXY, "P");
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "current track 2D", "current track 2D;x   (mm);y   (mm)", view2D);
  
  uHitsAll = new TGraphErrors(); uHitsAll->SetMarkerColor(1); uHitsAll->SetMarkerStyle(2); uHitsAll->SetMarkerSize(0.7); uHitsAll->SetName("uHitsAll");
  vHitsAll = new TGraphErrors(); vHitsAll->SetMarkerColor(1); vHitsAll->SetMarkerStyle(2); vHitsAll->SetMarkerSize(0.7); vHitsAll->SetName("vHitsAll");
  uHitsSel = new TGraphErrors(); uHitsSel->SetMarkerColor(8); uHitsSel->SetMarkerStyle(2); uHitsSel->SetMarkerSize(1.2); uHitsSel->SetName("uHitsSel");
  vHitsSel = new TGraphErrors(); vHitsSel->SetMarkerColor(8); vHitsSel->SetMarkerStyle(2); vHitsSel->SetMarkerSize(1.2); vHitsSel->SetName("vHitsSel");
  uTrack = new TF1("uTrack", "[0]*x+[1]", 0, 1); uTrack->SetLineColor(2); uTrack->SetLineWidth(2);
  vTrack = new TF1("vTrack", "[0]*x+[1]", 0, 1); vTrack->SetLineColor(4); vTrack->SetLineWidth(2);

  // typical z
  double z0 = rp_pos.z();
  char z0_str[30];
  sprintf(z0_str, "%.3f", z0 * 1E-3);

  vView = new MultiRootPlot();
  vView->Add(vHitsAll, "AP");
  vView->Add(vHitsSel, "P");
  vView->Add(vTrack, "SAME");
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "actual track - v projection",
    string() + "actual track - v projection;z (mm) - " + z0_str + " m;v   (mm)",  vView);

  uView = new MultiRootPlot();
  uView->Add(uHitsAll, "AP");
  uView->Add(uHitsSel, "P");
  uView->Add(uTrack, "SAME");
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "actual track - u projection",
    string() + "actual track - u projection;z (mm) - " + z0_str + " m;u   (mm)",  uView);
 
  scaleFixer = new TH2F("", "", 100, -25, +25, 100, -25, +25);
  scaleFixer->SetStats(false);
  
  uvGlobalView = new MultiRootPlot("", false);
  uvGlobalView->Add(scaleFixer, "");
  uvGlobalView->Add(uHitsSel, "P");
  uvGlobalView->Add(uTrack, "SAME");
  uvGlobalView->Add(vHitsSel, "P");
  uvGlobalView->Add(vTrack, "SAME");
  PlotManager::RegisterRP(TotRPDetId::lRP, id, "actual track - both projections",
    string() + "actual track - both projections;z (mm) - " + z0_str + " m;u, v   (mm)", uvGlobalView);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::PlanePlots::PlanePlots(unsigned int id)
{
  PlotManager::RegisterRP(TotRPDetId::lPlane, id, "digi profile cumulative", "digi profile cumulative;strip number", digi_profile_cumulative = new TH1D("", "title", 512, -0.5, 511.5));
  PlotManager::RegisterRP(TotRPDetId::lPlane, id, "digi profile one-event", "digi profile one-event;strip number", digi_profile_one_event = new TH1D("", "title", 512, -0.5, 511.5));
  PlotManager::RegisterRP(TotRPDetId::lPlane, id, "cluster profile cumulative", "cluster profile cumulative;cluster center" , cluster_profile_cumulative = new TH1D("", "title", 1024, -0.25, 511.75));
  PlotManager::RegisterRP(TotRPDetId::lPlane, id, "cluster profile one-event", "cluster profile one-event;cluster center" , cluster_profile_one_event = new TH1D("", "title", 1024, -0.25, 511.75));
  PlotManager::RegisterRP(TotRPDetId::lPlane, id, "hit multiplicity", "hit multiplicity;hits/detector/event" , hit_multiplicity = new TH1D("", "title", 6, -0.5, 5.5));
  PlotManager::RegisterRP(TotRPDetId::lPlane, id, "cluster size", "cluster size;hits per cluster" , cluster_size = new TH1D("", "title", 5, 0.5, 5.5));
}

#endif

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::TotemRPDQMSource(const edm::ParameterSet& ps)
  /*
  buildCorrelationPlots(ps.getUntrackedParameter<bool>("buildCorrelationPlots", false)),
  correlationPlotsFilter(ps.getUntrackedParameter<std::string>("correlationPlotsFilter", "")),
  correlationPlotsLimit(ps.getUntrackedParameter<unsigned int>("correlationPlotsLimit", 50)),
  correlationPlotsSelector(correlationPlotsFilter)
  */
{
  edm::LogInfo("TotemRPDQMSource") <<  "Constructor  TotemRPDQMSource::TotemRPDQMSource " << std::endl;
  
  tokenStripDigi = consumes< DetSetVector<RPStripDigi> >(ps.getParameter<edm::InputTag>("tagStripDigi"));
  tokenDigiCluster = consumes< edm::DetSetVector<RPDigCluster> >(ps.getParameter<edm::InputTag>("tagDigiCluster"));
  tokenRecoHit = consumes< edm::DetSetVector<RPRecoHit> >(ps.getParameter<edm::InputTag>("tagRecoHit"));
  tokenPatternColl = consumes< RPRecognizedPatternsCollection >(ps.getParameter<edm::InputTag>("tagPatternColl"));
  tokenTrackCandColl = consumes< RPTrackCandidateCollection >(ps.getParameter<edm::InputTag>("tagTrackCandColl"));
  tokenTrackColl = consumes< RPFittedTrackCollection >(ps.getParameter<edm::InputTag>("tagTrackColl"));
  tokenMultiTrackColl = consumes< RPMulFittedTrackCollection >(ps.getParameter<edm::InputTag>("tagMultiTrackColl"));
}

//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::~TotemRPDQMSource()
{
  edm::LogInfo("TotemRPDQMSource") <<  "Destructor TotemRPDQMSource::~TotemRPDQMSource " << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::dqmBeginRun(edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::beginRun" << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const &, edm::EventSetup const &)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::bookHistograms" << std::endl;
  
  ibooker.cd();
  ibooker.setCurrentFolder("TotemRP");

  //h_test = ibooker_.book1D("test", "some title", 40, -0.5, 39.5);

  // initialize diagonals
  diagonalPlots[1] = DiagonalPlots(ibooker, 1);  // 45 bot - 56 top
  diagonalPlots[2] = DiagonalPlots(ibooker, 2);  // 45 top - 45 bot

  // initialize anti-diagonals
  diagonalPlots[0] = DiagonalPlots(ibooker, 0);  // 45 bot - 56 bot
  diagonalPlots[3] = DiagonalPlots(ibooker, 3);  // 45 top - 56 top

  /*
  // loop over arms
  for (unsigned int arm = 0; arm < 2; arm++)
  {
    armPlots[arm] = ArmPlots(arm);

    // loop over stations
    for (unsigned int st = 0; st < 3; st += 2)
    {
      unsigned int stId = 10*arm + st;

      set<unsigned int> stationPlanes;

      // loop over RPs
      for (unsigned int rp = 0; rp < 6; ++rp)
      {
        unsigned int rpId = 10*stId + rp;

        potPlots[rpId] = PotPlots(rpId, geometry);

        // loop over planes
        for (unsigned int pl = 0; pl < 10; ++pl)
        {
          unsigned int plId = 10*rpId + pl;
          planePlots[plId] = PlanePlots(plId);

		  if (correlationPlotsSelector.IfCorrelate(plId))
		    stationPlanes.insert(plId % 100);
        }
      }

      stationPlots[stId] = StationPlots(stId, stationPlanes,
        buildCorrelationPlots, &correlationPlotsSelector, correlationPlotsLimit);
    }
  }
  */
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, 
                                            edm::EventSetup const& context) 
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::beginLuminosityBlock" << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::analyze(edm::Event const& event, edm::EventSetup const& eSetup)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::analyze" << std::endl;

  Handle< DetSetVector<RPStripDigi> > digi;
  event.getByToken(tokenStripDigi, digi);

  Handle< DetSetVector<RPDigCluster> > digCluster;
  event.getByToken(tokenDigiCluster, digCluster);

  Handle< DetSetVector<RPRecoHit> > hits;
  event.getByToken(tokenRecoHit, hits);

  Handle<RPRecognizedPatternsCollection> patterns;
  event.getByToken(tokenPatternColl, patterns);

  Handle< RPTrackCandidateCollection > trackCanColl;
  event.getByToken(tokenTrackCandColl, trackCanColl);

  Handle< RPFittedTrackCollection > tracks;
  event.getByToken(tokenTrackColl, tracks);

  Handle< RPMulFittedTrackCollection > multiTracks;
  event.getByToken(tokenMultiTrackColl, multiTracks);
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup) 
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::endLuminosityBlock" << std::endl;
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMSource::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
  edm::LogInfo("TotemRPDQMSource") <<  "TotemRPDQMSource::endRun" << std::endl;
}
