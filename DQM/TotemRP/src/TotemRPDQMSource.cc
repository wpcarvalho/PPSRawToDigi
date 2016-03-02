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

  ibooker.setCurrentFolder(string("Totem/RP/") + name);

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

#endif

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDQMSource::PlanePlots::PlanePlots(DQMStore::IBooker &ibooker, unsigned int id)
{
  ibooker.setCurrentFolder(string("Totem/") + TotRPDetId::PlaneName(id, TotRPDetId::nPath));

  digi_profile_cumulative = ibooker.book1D("digi profile", "digi profile;strip number", 512, -0.5, 511.5);
  cluster_profile_cumulative = ibooker.book1D("cluster profile", "cluster profile;cluster center", 1024, -0.25, 511.75);
  hit_multiplicity = ibooker.book1D("hit multiplicity", "hit multiplicity;hits/detector/event", 6, -0.5, 5.5);
  cluster_size = ibooker.book1D("cluster size", "cluster size;hits per cluster", 5, 0.5, 5.5);
}

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

  // loop over arms
  for (unsigned int arm = 0; arm < 2; arm++)
  {
    // TODO
    //armPlots[arm] = ArmPlots(arm);

    // loop over stations
    for (unsigned int st = 0; st < 3; st += 2)
    {
      unsigned int stId = 10*arm + st;

      // TODO
      //set<unsigned int> stationPlanes;

      // loop over RPs
      for (unsigned int rp = 0; rp < 6; ++rp)
      {
        unsigned int rpId = 10*stId + rp;

        // TODO
        //potPlots[rpId] = PotPlots(rpId, geometry);

        // loop over planes
        for (unsigned int pl = 0; pl < 10; ++pl)
        {
          unsigned int plId = 10*rpId + pl;
          planePlots[plId] = PlanePlots(ibooker, plId);

          // TODO
		  //if (correlationPlotsSelector.IfCorrelate(plId))
		    //stationPlanes.insert(plId % 100);
        }
      }

      // TODO
      //stationPlots[stId] = StationPlots(stId, stationPlanes,
        //buildCorrelationPlots, &correlationPlotsSelector, correlationPlotsLimit);
    }
  }
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

  // check validity
  bool valid = true;
  valid &= digi.isValid();
  valid &= digCluster.isValid();
  // TODO
  //valid &= hits.isValid();
  //valid &= trackCanColl.isValid();
  //valid &= tracks.isValid();
  //valid &= patterns.isValid();
  //valid &= multiTracks.isValid();

  if (!valid)
  {
    printf("ERROR in TotemDQMModuleRP::analyze > some of the required inputs are not valid. Skipping this event.\n");
    printf("\tdigi.isValid = %i\n", digi.isValid());
    printf("\tdigCluster.isValid = %i\n", digCluster.isValid());
    // TODO
    //printf("\thits.isValid = %i\n", hits.isValid());
    //printf("\ttrackCanColl.isValid = %i\n", trackCanColl.isValid());
    //printf("\ttracks.isValid = %i\n", tracks.isValid());
    //printf("\tpatterns.isValid = %i\n", patterns.isValid());
    //printf("\tmultiTracks.isValid = %i\n", multiTracks.isValid());

    return;
  }
  
  //-------------------------------------------------------------------------------------------------
  // Plane Plots

  // digi profile cumulative
  for (DetSetVector<RPStripDigi>::const_iterator it = digi->begin(); it != digi->end(); ++it)
  {
    unsigned int DetId = TotRPDetId::RawToDecId(it->detId());
    for (DetSet<RPStripDigi>::const_iterator dit = it->begin(); dit != it->end(); ++dit)
    {
      planePlots[DetId].digi_profile_cumulative->Fill(dit->GetStripNo());
    }
  }

  // cluster profile cumulative
  for (DetSetVector<RPDigCluster>::const_iterator it = digCluster->begin(); it != digCluster->end(); it++)
  {
    unsigned int DetId = TotRPDetId::RawToDecId(it->detId());
    for (DetSet<RPDigCluster>::const_iterator dit = it->begin(); dit != it->end(); ++dit)
    {
      planePlots[DetId].cluster_profile_cumulative->Fill(dit->CentreStripPos());
    }
  }

  // hit multiplicity
  for (DetSetVector<RPDigCluster>::const_iterator it = digCluster->begin(); it != digCluster->end(); it++)
  {
    unsigned int DetId = TotRPDetId::RawToDecId(it->detId());
    planePlots[DetId].hit_multiplicity->Fill(it->size());
  }

  // cluster size
  for (DetSetVector<RPDigCluster>::const_iterator it = digCluster->begin(); it != digCluster->end(); it++)
  {
    unsigned int DetId = TotRPDetId::RawToDecId(it->detId());
    for (DetSet<RPDigCluster>::const_iterator dit = it->begin(); dit != it->end(); ++dit)
      planePlots[DetId].cluster_size->Fill(dit->GetNumberOfStrips());
  }

#if 0
  //---------------------------------------------------------------------------------------------
  // Roman Pots Plots

  // plane activity histogram
  if (CumulativeMode())
  {
    map<unsigned int, set<unsigned int> > planes;
    map<unsigned int, set<unsigned int> > planes_u;
    map<unsigned int, set<unsigned int> > planes_v;
    for (DetSetVector<RPRecoHit>::const_iterator it = hits->begin(); it != hits->end(); ++it)
    {
      unsigned int DetId = TotRPDetId::RawToDecId(it->detId());
      unsigned int RPId = TotRPDetId::RPOfDet(DetId);
      unsigned int planeNum = DetId % 10;
      planes[RPId].insert(planeNum);
      if (TotRPDetId::IsStripsCoordinateUDirection(DetId))
        planes_u[RPId].insert(planeNum);
      else
        planes_v[RPId].insert(planeNum);
    }
    
    for (std::map<unsigned int, PotPlots>::iterator it = potPlots.begin(); it != potPlots.end(); it++)
    {
      it->second.activity->Fill(planes[it->first].size());
      it->second.activity_u->Fill(planes_u[it->first].size());
      it->second.activity_v->Fill(planes_v[it->first].size());
    }
  }
  
  if (CumulativeMode())
  {
    for (DetSetVector<RPDigCluster>::const_iterator it = digCluster->begin(); it != digCluster->end(); it++)
    {
      unsigned int DetId = TotRPDetId::RawToDecId(it->detId());
      unsigned int RPId = TotRPDetId::RPOfDet(DetId);
      unsigned int planeNum = DetId % 10;
      PotPlots &pp = potPlots[RPId];
      for (DetSet<RPDigCluster>::const_iterator dit = it->begin(); dit != it->end(); ++dit)
        pp.hit_plane_hist->Fill(planeNum, dit->CentreStripPos());   
    }
  }

  // recognized pattern histograms and event-category histogram
  if (CumulativeMode())
  {
    for (RPRecognizedPatternsCollection::const_iterator rpit = patterns->begin(); rpit != patterns->end(); ++rpit)
    {
      PotPlots &pp = potPlots[rpit->first];

      unsigned int u = rpit->second.uLines.size(), v = rpit->second.vLines.size();

      pp.patterns_u->Fill(u);
      pp.patterns_v->Fill(v);
    
      // determine category
      unsigned int category = 100;

      if (u == 0 && v == 0) category = 0;                   // empty
      if (u > 0 && v > 0 && u <= 3 && v <= 3) category = 3; // multi-track
      if (u+v == 1) category = 1;                           // insuff
      if (u == 1 && v == 1) category = 2;                   // 1-track
      if (u > 3 || v > 3) category = 4;                     // shower

      pp.event_category->Fill(category);
    }
  }

  // 3D plots cumulative and one-event & 3D graphs
  if (!CumulativeMode())
  {
    // plot all hits
    for (DetSetVector<RPRecoHit>::const_iterator it = hits->begin(); it != hits->end(); ++it)
    {
      unsigned int DetId = TotRPDetId::RawToDecId(it->detId());
      unsigned int RPId = TotRPDetId::RPOfDet(DetId);

      double z0 = geometry->GetRPPosition(RPId).z();

      for (DetSet<RPRecoHit>::const_iterator dit = it->begin(); dit != it->end(); ++dit)
      {
        unsigned int decId = TotRPDetId::RawToDecId(dit->DetId());
        double z = geometry->GetSensorPosition(decId).z();

        if (potPlots[RPId].uHitsAll == 0 && potPlots[RPId].vHitsAll == 0)
          continue;

        if (TotRPDetId::IsStripsCoordinateUDirection(TotRPDetId::RawToDecId(dit->DetId())))
        {
          Int_t pN = potPlots[RPId].uHitsAll->GetN();
          potPlots[RPId].uHitsAll->SetPoint(pN, z - z0, dit->Position());
          potPlots[RPId].uHitsAll->SetPointError(pN, 0., dit->Sigma());
        } else {
          Int_t pN = potPlots[RPId].vHitsAll->GetN();
          potPlots[RPId].vHitsAll->SetPoint(pN, z - z0, dit->Position());
          potPlots[RPId].vHitsAll->SetPointError(pN, 0., dit->Sigma());
        }
      }
    }
  }

  // draw recognized patterns
  if (drawRecognizedPatterns && !CumulativeMode())
  {
    for (RPRecognizedPatternsCollection::const_iterator rit = patterns->begin(); rit != patterns->end(); ++rit)
    {
      unsigned int RPId = rit->first;
      const PotPlots &pp = potPlots[RPId];
      double z0 = geometry->GetRPPosition(RPId).z();

      for (unsigned int i = 0; i < 2; i++)
      {
        MultiRootPlot *view = (i == 0) ? pp.uView : pp.vView;
        const vector<RPRecognizedPatterns::Line> *lines = (i == 0) ? &rit->second.uLines : &rit->second.vLines;
        int color = (i == 0) ? 2 : 4;
        
        for (vector<RPRecognizedPatterns::Line>::const_iterator lit = lines->begin(); lit != lines->end(); ++lit)
        {
          // skip "empty" patterns
          if (lit->hits.size() < 2)
            continue;

          // find z range
          double z_min = +std::numeric_limits<double>::max(), z_max = -std::numeric_limits<double>::max();
          for (RPRecognizedPatterns::Line::HitCollection::const_iterator hit = lit->hits.begin(); hit != lit->hits.end(); ++hit)
          {
            unsigned int decId = TotRPDetId::RawToDecId(hit->DetId());
            double z = geometry->GetSensorPosition(decId).z();
            z_min = min(z_min, z);
            z_max = max(z_max, z);
          }
          
          // create function
          TF1 *f = new TF1("", "[0]*x+[1]", 0, 1);
          f->SetParameters(lit->a, lit->b);
          f->SetRange(z_min - z0, z_max - z0);
          f->SetLineColor(color);
          f->SetLineStyle(2);
          f->SetLineWidth(1);

          view->Add(f, "same");
        }
      }
    }
  }


  // u and v views
  if (!CumulativeMode())
  {
    for (RPFittedTrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
    {
      unsigned int RPId = it->first;
      const RPFittedTrack &ft = it->second;
      double z0 = geometry->GetRPPosition(RPId).z();

      // to set range for u,v Track
      double u_min_z = +std::numeric_limits<double>::max(), u_max_z = -std::numeric_limits<double>::max(), v_min_z = u_min_z, v_max_z = u_max_z;

      // skip invalid tracks
      if (!ft.IsValid())
        continue;

      if (!ft.IsSourceTrackCandidateValid())
        continue;

      // build u and v projection views
      const PotPlots & pp = potPlots[RPId];

#ifdef DEBUG
      printf("------- RPId = %u\n", RPId);
#endif
#ifdef DEBUG
      cout << ">> range0 u: min: " << u_min_z << "\tmax: " << u_max_z << endl;
      cout << ">> range0 v: min: " << v_min_z << "\tmax: " << v_max_z << endl;
#endif

      // plot selected hits and establish track range
      const vector<RPRecoHit> &hits = ft.sourceTrackCandidate.TrackRecoHits();
      for (vector<RPRecoHit>::const_iterator it = hits.begin(); it != hits.end(); ++it)
      {
        unsigned int decId = TotRPDetId::RawToDecId(it->DetId());
        double z = geometry->GetSensorPosition(decId).z();

        if (TotRPDetId::IsStripsCoordinateUDirection(TotRPDetId::RawToDecId(it->DetId())))
        {
          Int_t pN = potPlots[RPId].uHitsSel->GetN();
          pp.uHitsSel->SetPoint(pN, z - z0, it->Position());
          pp.uHitsSel->SetPointError(pN, 0, it->Sigma());
          u_min_z = std::min(u_min_z, z);
          u_max_z = std::max(u_max_z, z);
          //printf("DrawUHit(%u, %.3f, %.3f, %.3f, %.3f)\n", RPId, z0, z - z0, it->Position(), it->Sigma());
        } else {
          Int_t pN = potPlots[RPId].vHitsSel->GetN();
          pp.vHitsSel->SetPoint(pN, z - z0, it->Position());
          pp.vHitsSel->SetPointError(pN, 0, it->Sigma());
          v_min_z = std::min(v_min_z, z);
          v_max_z = std::max(v_max_z, z);
          //printf("DrawVHit(%u, %.3f, %.3f, %.3f, %.3f)\n", RPId, z0, z - z0, it->Position(), it->Sigma());
        }
      }

#ifdef DEBUG
      cout << ">> range u: min: " << u_min_z << "\tmax: " << u_max_z << endl;
      cout << ">> range v: min: " << v_min_z << "\tmax: " << v_max_z << endl;
#endif

      // set the ranges
      pp.uTrack->SetRange(u_min_z - z0, u_max_z - z0);
      pp.vTrack->SetRange(v_min_z - z0, v_max_z - z0);

      // get (x, y) parameters in local system
      double ax = ft.GetTx(), ay = ft.GetTy();
      double bx = ft.X0() + ax * (z0 - ft.Z0());
      double by = ft.Y0() + ay * (z0 - ft.Z0());

      //printf("DrawXYFit(%u, %.3f, %.3f, %.3f, %.3f, %.3f)\n", RPId, z0, ax, ay, bx, by);

#ifdef DEBUG
      printf(">> original fit: Z0 = %f, X0 = %f, Tx = %f, Y0 = %f, Ty = %f\n", ft.Z0(), ft.X0(), ft.GetTx(), ft.Y0(), ft.GetTy());
      printf(">> z0 = %f\n", z0);
      printf(">> x fit: a = %f, b = %f\n", ax, bx);
      printf(">> y fit: a = %f, b = %f\n", ay, by);
      printf(">> track at z0: x = %f, y = %f\n", ft.GetTrackPoint(z0).X(), ft.GetTrackPoint(z0).Y());
#endif

      // mean read-out direction of U and V planes
      CLHEP::Hep3Vector rod_U = geometry->GetRPMeanUDirection(RPId);
      CLHEP::Hep3Vector rod_V = geometry->GetRPMeanVDirection(RPId);

      // mean position of U and V planes
      CLHEP::Hep3Vector mp_U = geometry->GetSensorPosition(RPId*10 + 1);  // one U plane
      CLHEP::Hep3Vector mp_V = geometry->GetSensorPosition(RPId*10 + 0);  // one V plane

      // convert (x, y) to (u, v) parameters
      double au = ax * rod_U.x() + ay * rod_U.y();
      double av = ax * rod_V.x() + ay * rod_V.y();

      double bu = (bx - mp_U.x()) * rod_U.x() + (by - mp_U.y()) * rod_U.y();
      double bv = (bx - mp_V.x()) * rod_V.x() + (by - mp_V.y()) * rod_V.y();
      
      //printf(">> fit in RP %u\n", RPId);
      //printf("\tU: a = %f, b = %f\n", au, bu);
      //printf("\tV: a = %f, b = %f\n", av, bv);

      // set the track parameters
      pp.uTrack->SetParameters(au, bu);
      pp.vTrack->SetParameters(av, bv);
    }
  }

  // 3D track and 2D hit plots
  for (RPFittedTrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
  {
    const RPFittedTrack &ft = it->second;
    
    if (!ft.IsValid())
      continue;
    
    const PotPlots & pp = potPlots[it->first];

    if (!CumulativeMode())
    {
      pp.currentTrackInRP->AddTrack(ft.GetTx(), ft.X0(), ft.GetTy(), ft.Y0());
      pp.currentTrackXY->SetPoint(pp.currentTrackXY->GetN(), ft.X0(), ft.Y0());
    } else {
      pp.allTracksInRP->AddTrack(ft.GetTx(), ft.X0(), ft.GetTy(), ft.Y0());
    }
  }

  if (!CumulativeMode() && multiTracks.isValid())
  {
    for (RPMulFittedTrackCollection::const_iterator rpit = multiTracks->begin(); rpit != multiTracks->end(); ++rpit)
    {
      const PotPlots &pp = potPlots[rpit->first];

      for (vector<RPFittedTrack>::const_iterator trit = rpit->second.begin(); trit != rpit->second.end(); ++trit)
      {
        const RPFittedTrack &ft = *trit;
        if (!ft.IsValid())
          continue;

        pp.currentMultiTracksXY->SetPoint(pp.currentMultiTracksXY->GetN(), ft.X0(), ft.Y0());
      }
    }
  }

  // cumulative RP fit plots
  if (CumulativeMode())
  {
	for (RPFittedTrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
    {
	  const RPFittedTrack &ft = it->second;
	    
	  if (!ft.IsValid())
	    continue;
    
      unsigned int RPId = it->first;
      PotPlots &pp = potPlots[RPId];

      // number of planes contributing to (valid) fits
      unsigned int n_pl_in_fit_u = 0, n_pl_in_fit_v = 0;
      for (int hi = 0; hi < ft.GetHitEntries(); hi++)
      {
        unsigned int rawId = ft.GetHit(hi).DetId();  
        unsigned int decId = TotRPDetId::RawToDecId(rawId);
        if (TotRPDetId::IsStripsCoordinateUDirection(decId))
          n_pl_in_fit_u++;
        else
          n_pl_in_fit_v++;
      }
      pp.h_planes_fit_u->Fill(n_pl_in_fit_u);
      pp.h_planes_fit_v->Fill(n_pl_in_fit_v);

      // mean position of U and V planes
      CLHEP::Hep3Vector mp_U = geometry->GetSensorPosition(RPId*10 + 1);  // one U plane
      CLHEP::Hep3Vector mp_V = geometry->GetSensorPosition(RPId*10 + 0);  // one V plane

      double rp_x = ( mp_U.x() + mp_V.x() ) / 2.;
      double rp_y = ( mp_U.y() + mp_V.y() ) / 2.;

      // mean read-out direction of U and V planes
      CLHEP::Hep3Vector rod_U = geometry->GetRPMeanUDirection(RPId);
      CLHEP::Hep3Vector rod_V = geometry->GetRPMeanVDirection(RPId);

      double x = ft.X0() - rp_x;
      double y = ft.Y0() - rp_y;

      pp.trackHitsCumulative->SetPoint(pp.trackHitsCumulative->GetN(), x, y);
      pp.trackHitsCumulativeHist->Fill(x, y);

      double U = x * rod_U.x() + y * rod_U.y();
      double V = x * rod_V.x() + y * rod_V.y();

      pp.track_u_profile->Fill(U);
      pp.track_v_profile->Fill(V);
    }
  }

  //-----------------------------------------------------------------------------------
  // Station Plots

  // Correlation profile
  if (CumulativeMode() && buildCorrelationPlots)
  {
    for (DetSetVector<RPStripDigi>::const_iterator i = digi->begin(); i != digi->end(); i++)
    {
      for (DetSetVector<RPStripDigi>::const_iterator j = i; j != digi->end(); j++)
      {
        if (i == j)
          continue;

        unsigned int DetId1 = TotRPDetId::RawToDecId(i->detId());
        unsigned int DetId2 = TotRPDetId::RawToDecId(j->detId());
        unsigned int StationId1 = DetId1 / 100;
        unsigned int StationId2 = DetId2 / 100;

        if (StationId1 != StationId2)
          continue;

        unsigned int RPPlaneId1 = DetId1 % 100;
        unsigned int RPPlaneId2 = DetId2 % 100;
        if (stationPlots[StationId1].hist[RPPlaneId1][RPPlaneId2])
        {
          for (DetSet<RPStripDigi>::const_iterator di = i->begin(); di != i->end(); di++)
          {
            for (DetSet<RPStripDigi>::const_iterator dj = j->begin(); dj != j->end(); dj++)
            {
              Double_t temp[2];
              temp[0] = di->GetStripNo();
              temp[1] = dj->GetStripNo();
              stationPlots[StationId1].hist[RPPlaneId1][RPPlaneId2]->Fill(temp);
            }
          }
        }

      }
    }
  }
  
  // plot of hits in a station
  if (!CumulativeMode())
  {
	for (RPFittedTrackCollection::const_iterator it = tracks->begin(); it != tracks->end(); ++it)
    {
	  const RPFittedTrack &ft = it->second;
	  
	  if (!ft.IsValid())
	    continue;
	
	  unsigned int sId = it->first / 10;
	  TGraph *g = stationPlots[sId].rpHits;
	  g->SetPoint(g->GetN(), ft.X0(), ft.Y0());
      //printf("%u: %f, %f\n", it->first, ft.X0(), ft.Y0());
	}
  }
  
  //-----------------------------------------------------------------------------------
  // Arm Plots
  if (CumulativeMode())
  {
    map<unsigned int, unsigned int> mTop, mHor, mBot;

    for (auto p : armPlots)
    {
      mTop[p.first] = 0;
      mHor[p.first] = 0;
      mBot[p.first] = 0;
    }

    for (auto p : *tracks)
    {
      if (!p.second.IsValid())
        continue;
  
      unsigned int armNum = p.first / 100;
      unsigned int rpNum = p.first % 10;

      if (rpNum == 0 || rpNum == 4)
        mTop[armNum]++;
      if (rpNum == 2 || rpNum == 3)
        mHor[armNum]++;
      if (rpNum == 1 || rpNum == 5)
        mBot[armNum]++;
    }

    for (auto &p : armPlots)
    {
      p.second.h_numRPWithTrack_top->Fill(mTop[p.first]);
      p.second.h_numRPWithTrack_hor->Fill(mHor[p.first]);
      p.second.h_numRPWithTrack_bot->Fill(mBot[p.first]);
    }

    // track RP correlation
    for (auto t1 : *tracks)
    {
      if (!t1.second.IsValid())
        continue;

      unsigned int arm1 = t1.first / 100;
      unsigned int stNum1 = (t1.first / 10) % 10;
      unsigned int rpNum1 = t1.first % 10;
      unsigned int idx1 = stNum1/2 * 7 + rpNum1;
      bool hor1 = (rpNum1 == 2 || rpNum1 == 3);

      ArmPlots &ap = armPlots[arm1];

      for (auto t2 : *tracks)
      {
        if (!t2.second.IsValid())
          continue;
      
        unsigned int arm2 = t2.first / 100;
        unsigned int stNum2 = (t2.first / 10) % 10;
        unsigned int rpNum2 = t2.first % 10;
        unsigned int idx2 = stNum2/2 * 7 + rpNum2;
        bool hor2 = (rpNum2 == 2 || rpNum2 == 3);

        if (arm1 != arm2)
          continue;

        ap.h_trackCorr->Fill(idx1, idx2); 
        
        if (hor1 != hor2)
          ap.h_trackCorr_overlap->Fill(idx1, idx2); 
      }
    }
  }
  
  //---------------------------------------------------------------------------------
  // RP-system plots
  // TODO
  if (CumulativeMode())
  {
    for (auto &dp : diagonalPlots)
    {
      unsigned int id = dp.first;
      bool top45 = id & 2;
      bool top56 = id & 1;

      unsigned int id_45_n = (top45) ? 20 : 21;
      unsigned int id_45_f = (top45) ? 24 : 25;
      unsigned int id_56_n = (top56) ? 120 : 121;
      unsigned int id_56_f = (top56) ? 124 : 125;
    
      bool h_45_n = (tracks->find(id_45_n) != tracks->end() && tracks->find(id_45_n)->second.IsValid());
      bool h_45_f = (tracks->find(id_45_f) != tracks->end() && tracks->find(id_45_f)->second.IsValid());
      bool h_56_n = (tracks->find(id_56_n) != tracks->end() && tracks->find(id_56_n)->second.IsValid());
      bool h_56_f = (tracks->find(id_56_f) != tracks->end() && tracks->find(id_56_f)->second.IsValid());
    
      if (! (h_45_n && h_45_f && h_56_n && h_56_f) )
        continue;

      double x_45_n = tracks->find(id_45_n)->second.X0(), y_45_n = tracks->find(id_45_n)->second.Y0();
      double x_45_f = tracks->find(id_45_f)->second.X0(), y_45_f = tracks->find(id_45_f)->second.Y0();
      double x_56_n = tracks->find(id_56_n)->second.X0(), y_56_n = tracks->find(id_56_n)->second.Y0();
      double x_56_f = tracks->find(id_56_f)->second.X0(), y_56_f = tracks->find(id_56_f)->second.Y0();

      double dx_45 = x_45_f - x_45_n;
      double dy_45 = y_45_f - y_45_n;
      double dx_56 = x_56_f - x_56_n;
      double dy_56 = y_56_f - y_56_n;

      DiagonalPlots &pl = dp.second;

      int idx = pl.g_lrc_x_d->GetN();
      pl.g_lrc_x_d->SetPoint(idx, dx_45, dx_56);  
      pl.g_lrc_y_d->SetPoint(idx, dy_45, dy_56);  
      
      pl.g_lrc_x_n->SetPoint(idx, x_45_n, x_56_n);  
      pl.g_lrc_y_n->SetPoint(idx, y_45_n, y_56_n);  
      
      pl.g_lrc_x_f->SetPoint(idx, x_45_f, x_56_f);  
      pl.g_lrc_y_f->SetPoint(idx, y_45_f, y_56_f);  


      pl.h_lrc_x_d->Fill(dx_45, dx_56);  
      pl.h_lrc_y_d->Fill(dy_45, dy_56);  
      
      pl.h_lrc_x_n->Fill(x_45_n, x_56_n);  
      pl.h_lrc_y_n->Fill(y_45_n, y_56_n);  
      
      pl.h_lrc_x_f->Fill(x_45_f, x_56_f);  
      pl.h_lrc_y_f->Fill(y_45_f, y_56_f);  
    }
  }
#endif

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
