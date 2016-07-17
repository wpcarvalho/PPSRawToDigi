/****************************************************************************
*
* This is a part of TotemDQM and TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*   Rafał Leszko (rafal.leszko@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DQMServices/Core/interface/DQMEDAnalyzer.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/TotemRPDetId/interface/TotemRPDetId.h"

//Diamond specific
#include "DataFormats/CTPPSDigi/interface/DiamondDigi.h"
#include "DataFormats/CTPPSDigi/interface/DiamondVFATStatus.h"

// #include "DataFormats/CTPPSReco/interface/TotemRPCluster.h"
// #include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
// #include "DataFormats/CTPPSReco/interface/TotemRPUVPattern.h"
//#include "RecoTotemRP/RPRecoDataFormats/interface/RPMulFittedTrackCollection.h"

//Correlation with tracks in Strips
#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"

#include "Geometry/Records/interface/VeryForwardRealGeometryRecord.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"

#include <string>

//----------------------------------------------------------------------------------------------------

namespace TotemRPDiamondDQM {
  bool TotemRP_isDiamond(const int& detectorSymbolicID) {
    return true;
  }
  
  std::string TotemRP_getDiamondName(const int& detectorSymbolicID) {
    std::string name;
    return name;
  }
}
 
class TotemRPDiamondDQMSource: public DQMEDAnalyzer
{
  public:
    TotemRPDiamondDQMSource(const edm::ParameterSet& ps);
    virtual ~TotemRPDiamondDQMSource();
  
  protected:
    void dqmBeginRun(edm::Run const &, edm::EventSetup const &) override;
    void bookHistograms(DQMStore::IBooker &, edm::Run const &, edm::EventSetup const &) override;
    void analyze(edm::Event const& e, edm::EventSetup const& eSetup);
    void beginLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup);
    void endLuminosityBlock(edm::LuminosityBlock const& lumi, edm::EventSetup const& eSetup);
    void endRun(edm::Run const& run, edm::EventSetup const& eSetup);

  private:
    unsigned int verbosity;

    edm::EDGetTokenT< edm::DetSetVector<DiamondVFATStatus> > tokenStatus;
    edm::EDGetTokenT< edm::DetSetVector<DiamondDigi> > tokenDigi;
    edm::EDGetTokenT< edm::DetSetVector<TotemRPLocalTrack> > tokenLocalTrack;

    /// plots related to the whole system
    struct GlobalPlots
    {
      MonitorElement *h_trackCorr = NULL;

      void Init(DQMStore::IBooker &ibooker);
    };

    GlobalPlots globalPlots;

    /// plots related to one arm
    struct ArmPlots
    {
      int id;

      MonitorElement *h_numRPWithTrack=NULL;

      ArmPlots(){}

      ArmPlots(DQMStore::IBooker &ibooker, int _id);
    };

    std::map<unsigned int, ArmPlots> armPlots;

    /// plots related to one station
    struct StationPlots
    {
      StationPlots() {}
      StationPlots(DQMStore::IBooker &ibooker, int _id);
    };

    std::map<unsigned int, StationPlots> stationPlots;

    /// plots related to one RP
    struct PotPlots
    {
      MonitorElement *frame_problem=NULL, *frame_missing=NULL;

      MonitorElement *activity=NULL;
      MonitorElement *activity_per_bx=NULL, *activity_per_bx_short=NULL;
      MonitorElement *hit_distribution_2d=NULL;

      PotPlots() {}
      PotPlots(DQMStore::IBooker &ibooker, unsigned int id);
    };

    std::map<unsigned int, PotPlots> potPlots;
    
    /// plots related to one Digitizer Board
    struct DigitizerBoardPlots
    {
      MonitorElement *activity=NULL;
      MonitorElement *activity_per_bx=NULL, *activity_per_bx_short=NULL;

      DigitizerBoardPlots() {}
      DigitizerBoardPlots(DQMStore::IBooker &ibooker, unsigned int id);
    };

    std::map<unsigned int, DigitizerBoardPlots> digitizerBoardPlots;

    /// plots related to one RP plane
    struct PlanePlots
    {
      MonitorElement *digi_profile_cumulative = NULL;
      MonitorElement *hit_multiplicity = NULL;
      MonitorElement *threshold_voltage = NULL;
//       MonitorElement *time_over_threshold_cumulative_2d = NULL;
//       MonitorElement *hit_rate_2d = NULL;
      

      PlanePlots() {}
      PlanePlots(DQMStore::IBooker &ibooker, unsigned int id);
    };

    std::map<unsigned int, PlanePlots> planePlots;
    
    /// plots related to one RP plane
    struct ChannelPlots
    {
      MonitorElement *error_flags = NULL;
      MonitorElement *leading_edge_cumulative = NULL;
      MonitorElement *time_over_threshold_cumulative = NULL;
      MonitorElement *leading_trailing_correlation = NULL;
      MonitorElement *hit_rate = NULL;

      ChannelPlots() {}
      ChannelPlots(DQMStore::IBooker &ibooker, unsigned int id);
    };

    std::map<unsigned int, ChannelPlots> channelPlots;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

void TotemRPDiamondDQMSource::GlobalPlots::Init(DQMStore::IBooker &ibooker)
{
  std::map<int,std::string> correlationMap;
  
  h_trackCorr = ibooker.book2D("track correlation RP-allCTPPS", "rp, Track+Timing, hor", 6, -0.5, 5.5, 6, -0.5, 5.5);
  TH2F *hist = h_trackCorr->getTH2F();
  TAxis *xa = hist->GetXaxis(), *ya = hist->GetYaxis();
  xa->SetBinLabel(1, "45, 210, near"); ya->SetBinLabel(1, "45, 210, near");
  xa->SetBinLabel(2, "45, 210, far"); ya->SetBinLabel(2, "45, 210, far");
  xa->SetBinLabel(3, "45, cyl"); ya->SetBinLabel(3, "45, cyl");
  xa->SetBinLabel(4, "56, 210, near"); ya->SetBinLabel(4, "56, 210, near");
  xa->SetBinLabel(5, "56, 210, far"); ya->SetBinLabel(5, "56, 210, far");
  xa->SetBinLabel(6, "56, cyl"); ya->SetBinLabel(6, "56, cyl");
}

//----------------------------------------------------------------------------------------------------

TotemRPDiamondDQMSource::ArmPlots::ArmPlots(DQMStore::IBooker &ibooker, int _id) : id(_id)
{
  string path = TotemRPDetId::armName(id, TotemRPDetId::nPath);
  path.replace(0, 2, "TimingDiamond");
  ibooker.setCurrentFolder(string("CTPPS/") + path);

  string title = TotemRPDetId::armName(id, TotemRPDetId::nFull);

  h_numRPWithTrack = ibooker.book1D("number of hor+cyl RPs with tracks", title+";number of hor+cyl RPs with tracks", 7, -0.5, 6.5);
}

//----------------------------------------------------------------------------------------------------

TotemRPDiamondDQMSource::StationPlots::StationPlots(DQMStore::IBooker &ibooker, int id) 
{
  string path = TotemRPDetId::stationName(id, TotemRPDetId::nPath);
  path.replace(0, 2, "TimingDiamond");
  ibooker.setCurrentFolder(string("CTPPS/") + path);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDiamondDQMSource::PotPlots::PotPlots(DQMStore::IBooker &ibooker, unsigned int id)
{
  string path = TotemRPDetId::rpName(id, TotemRPDetId::nPath);
  path.replace(0, 2, "TimingDiamond");
  ibooker.setCurrentFolder(string("CTPPS/") + path);

  string title = TotemRPDetId::rpName(id, TotemRPDetId::nFull);

  activity = ibooker.book1D("active planes", title+";number of active planes", 5, -0.5, 4.5);


  activity_per_bx = ibooker.book1D("activity per BX", title+";Event.BX", 4002, -1.5, 4000. + 0.5);
  activity_per_bx_short = ibooker.book1D("activity per BX (short)", title+";Event.BX", 102, -1.5, 100. + 0.5);

  hit_distribution_2d = ibooker.book2D("activity in planes (2D)", title+";plane number;channel number", 4, 0.5, 4.5, 12, +0.5, 12.5);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDiamondDQMSource::DigitizerBoardPlots::DigitizerBoardPlots(DQMStore::IBooker &ibooker, unsigned int id)
{
  string path = TotemRPDetId::rpName(id, TotemRPDetId::nPath);
  path.replace(0, 2, "TimingDiamond");
  ibooker.setCurrentFolder(string("CTPPS/") + path);

  string title = TotemRPDetId::rpName(id, TotemRPDetId::nFull);		//TODO: Digitizer Board ID! dbName?

  activity = ibooker.book1D("active planes", title+";number of active planes", 5, -0.5, 4.5);


  activity_per_bx = ibooker.book1D("activity per BX (DB)", title+";Event.BX", 4002, -1.5, 4000. + 0.5);
  activity_per_bx_short = ibooker.book1D("activity per BX (short) (DB)", title+";Event.BX", 102, -1.5, 100. + 0.5);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDiamondDQMSource::PlanePlots::PlanePlots(DQMStore::IBooker &ibooker, unsigned int id)
{
  string path = TotemRPDetId::planeName(id, TotemRPDetId::nPath);
  path.replace(0, 2, "TimingDiamond");
  ibooker.setCurrentFolder(string("CTPPS/") + path);

  string title = TotemRPDetId::planeName(id, TotemRPDetId::nFull);

  digi_profile_cumulative = ibooker.book1D("digi profile", title+";channel number", 12, +0.5, 12.5);
  hit_multiplicity = ibooker.book1D("hit multiplicity", title+";hits/detector/event", 13, -0.5, 12.5);
  
  threshold_voltage = ibooker.book1D("Threshold Voltage", title+";threshold voltage", 12, +0.5, 12.5);
}

TotemRPDiamondDQMSource::ChannelPlots::ChannelPlots(DQMStore::IBooker &ibooker, unsigned int id)
{
  string path = TotemRPDetId::planeName(id, TotemRPDetId::nPath);	//TODO
  path.replace(0, 2, "TimingDiamond");
  ibooker.setCurrentFolder(string("CTPPS/") + path);

  string title = TotemRPDetId::planeName(id, TotemRPDetId::nFull);	//TODO

  error_flags = ibooker.book1D("digi profile", title+";channel number", 16, -0.5, 15.5);
  for (int error_index=1; error_index<16; ++error_index) 
    error_flags->getTH1F()->GetXaxis()->SetBinLabel(error_index, "");
  error_flags->getTH1F()->GetXaxis()->SetBinLabel(16, "MH");
  
  leading_edge_cumulative = ibooker.book1D("Leading Edge", title+";leading edge", 100, -0.5e-9, 19.5e-3);
  time_over_threshold_cumulative = ibooker.book1D("Time over Threshold", title+";time over threshold", 100, -0.5e-9, 19.5e-3);
  leading_trailing_correlation = ibooker.book2D("Leading Trailing Correlation", title+";leading trailing corr", 100, -0.5e-9, 19.5e-3, 100, -0.5e-9, 19.5e-3);
  hit_rate = ibooker.book1D("hit rate", title+";hit rate", 20, -0.5, 1.5);	// Hz?
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

TotemRPDiamondDQMSource::TotemRPDiamondDQMSource(const edm::ParameterSet& ps) :
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0))
{
  tokenStatus = consumes<DetSetVector<DiamondVFATStatus>>(ps.getParameter<edm::InputTag>("tagStatus"));

  tokenDigi = consumes< DetSetVector<DiamondDigi> >(ps.getParameter<edm::InputTag>("tagDigi"));
  tokenLocalTrack = consumes< DetSetVector<TotemRPLocalTrack> >(ps.getParameter<edm::InputTag>("tagLocalTrack"));
  //tokenMultiTrackColl = consumes< RPMulFittedTrackCollection >(ps.getParameter<edm::InputTag>("tagMultiTrackColl"));
}

//----------------------------------------------------------------------------------------------------

TotemRPDiamondDQMSource::~TotemRPDiamondDQMSource()
{
}

//----------------------------------------------------------------------------------------------------

void TotemRPDiamondDQMSource::dqmBeginRun(edm::Run const &, edm::EventSetup const &)
{
}

//----------------------------------------------------------------------------------------------------

void TotemRPDiamondDQMSource::bookHistograms(DQMStore::IBooker &ibooker, edm::Run const &, edm::EventSetup const &)
{
  ibooker.cd();
  ibooker.setCurrentFolder("CTPPS");

  // global plots
  globalPlots.Init(ibooker);

  // loop over arms
  for (unsigned int arm = 0; arm < 2; arm++)
  {
    armPlots[arm] = ArmPlots(ibooker, arm);

    // loop over stations
    for (unsigned int st = 0; st < 3; st += 2)
    {
      unsigned int stId = 10*arm + st;
      
//       stationPlots[stId] = StationPlots(ibooker, stId);	//not used for Diamond

      // loop over RPs with diamonds
      for (unsigned int rp = 0; rp < 1; ++rp)
      {
        unsigned int rpId = 10*stId + rp;

        potPlots[rpId] = PotPlots(ibooker, rpId);

	// loop over Digitizer Boards
	for (unsigned int db = 1; db <= 2; ++db)
	{
	  unsigned int dbId = 10*rpId + db;

	  digitizerBoardPlots[dbId] = DigitizerBoardPlots(ibooker, dbId);
	}

	// loop over planes
	for (unsigned int pl = 1; pl <= 4; ++pl)
	{
	  unsigned int plId = 100*rpId + pl;
	  planePlots[plId] = PlanePlots(ibooker, plId);
	  
	  // loop over channels
	  for (unsigned int ch = 1; ch <= 12; ++ch)
	  {
	    unsigned int chId = 100*plId + ch;
	    cout<<chId<<std::endl;
	    channelPlots[chId] = ChannelPlots(ibooker, chId);
	  }
	  
	}
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------

void TotemRPDiamondDQMSource::beginLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& context) 
{
}

//----------------------------------------------------------------------------------------------------

void TotemRPDiamondDQMSource::analyze(edm::Event const& event, edm::EventSetup const& eventSetup)
{
  
  // get event setup data
  ESHandle<TotemRPGeometry> geometry;
  eventSetup.get<VeryForwardRealGeometryRecord>().get(geometry);

  // get event data
  Handle< DetSetVector<DiamondVFATStatus> > status;
  event.getByToken(tokenStatus, status);

  Handle< DetSetVector<DiamondDigi> > digi;
  event.getByToken(tokenDigi, digi);

  Handle< DetSetVector<TotemRPLocalTrack> > tracks;
  event.getByToken(tokenLocalTrack, tracks);

  //Handle< RPMulFittedTrackCollection > multiTracks;
  //event.getByToken(tokenMultiTrackColl, multiTracks);

  // check validity
  bool valid = true;
  valid &= status.isValid();
  valid &= digi.isValid();
  valid &= tracks.isValid();
  //valid &= multiTracks.isValid();

  if (!valid)
  {
    if (verbosity)
    {
      LogProblem("TotemRPDiamondDQMSource") <<
        "ERROR in TotemDQMModuleRP::analyze > some of the required inputs are not valid. Skipping this event.\n"
        << "    status.isValid = " << status.isValid() << "\n"
        << "    digi.isValid = " << digi.isValid() << "\n"
        << "    tracks.isValid = " << tracks.isValid();
    }

    return;
  }
  
  
  //------------------------------
  // Global Plots

  

  //------------------------------
  // Status Plots
/*
  for (auto &ds : *status)
  {
    unsigned int decId = TotemRPDetId::rawToDecId(ds.detId());
    unsigned int rpId = decId / 10;
    unsigned int plNum = decId % 10;

    auto &plots = potPlots[rpId];

    for (auto &s : ds)
    {
      if (s.isMissing())
      {
        plots.frame_problem->Fill(plNum, s.getChipPosition());
        plots.frame_missing->Fill(plNum, s.getChipPosition());
      }


      if (s.isIDMismatch() || s.isFootprintError() || s.isCRCError())
      {
        plots.frame_problem->Fill(plNum, s.getChipPosition());
      }
    }
  }
  
  //------------------------------
  // Plane Plots

  // digi profile cumulative TODO
  for (DetSetVector<DiamondDigi>::const_iterator it = digi->begin(); it != digi->end(); ++it)
  {
    unsigned int DetId = TotemRPDetId::rawToDecId(it->detId());
    for (DetSet<DiamondDigi>::const_iterator dit = it->begin(); dit != it->end(); ++dit)
    {
      planePlots[dit->getCHID].digi_profile_cumulative->Fill(dit->getCHID);
    }
  }

  // hit multiplicity
  for (DetSetVector<TotemRPCluster>::const_iterator it = digCluster->begin(); it != digCluster->end(); it++)
  {
    unsigned int DetId = TotemRPDetId::rawToDecId(it->detId());
    planePlots[DetId].hit_multiplicity->Fill(it->size());
  }

  //------------------------------
  // Roman Pots Plots

  // determine active planes (from RecHits and VFATStatus)
  map<unsigned int, set<unsigned int> > planes;
  map<unsigned int, set<unsigned int> > planes_u;
  map<unsigned int, set<unsigned int> > planes_v;
  for (const auto &ds : *hits)
  {
    if (ds.empty())
      continue;

    unsigned int DetId = TotemRPDetId::rawToDecId(ds.detId());
    unsigned int RPId = TotemRPDetId::rpOfDet(DetId);
    unsigned int planeNum = DetId % 10;
    planes[RPId].insert(planeNum);
    if (TotemRPDetId::isStripsCoordinateUDirection(DetId))
      planes_u[RPId].insert(planeNum);
    else
      planes_v[RPId].insert(planeNum);
  }

  for (const auto &ds : *status)
  {
    bool activity = false;
    for (const auto &s : ds)
    {
      if (s.isNumberOfClustersSpecified() && s.getNumberOfClusters() > 0)
      {
        activity = true;
        break;
      }
    } 

    if (!activity)
      continue;

    unsigned int DetId = TotemRPDetId::rawToDecId(ds.detId());
    unsigned int RPId = TotemRPDetId::rpOfDet(DetId);
    unsigned int planeNum = DetId % 10;
    planes[RPId].insert(planeNum);
    if (TotemRPDetId::isStripsCoordinateUDirection(DetId))
      planes_u[RPId].insert(planeNum);
    else
      planes_v[RPId].insert(planeNum);
  }

  // plane activity histogram
  for (std::map<unsigned int, PotPlots>::iterator it = potPlots.begin(); it != potPlots.end(); it++)
  {
    it->second.activity->Fill(planes[it->first].size());


    if (planes[it->first].size() >= 6)
    {
      it->second.activity_per_bx->Fill(event.bunchCrossing());
      it->second.activity_per_bx_short->Fill(event.bunchCrossing());
    }
  }
  
  for (DetSetVector<TotemRPCluster>::const_iterator it = digCluster->begin(); it != digCluster->end(); it++)
  {
    unsigned int DetId = TotemRPDetId::rawToDecId(it->detId());
    unsigned int RPId = TotemRPDetId::rpOfDet(DetId);
    unsigned int planeNum = DetId % 10;
    PotPlots &pp = potPlots[RPId];
//     for (DetSet<TotemRPCluster>::const_iterator dit = it->begin(); dit != it->end(); ++dit)
//       pp.hit_plane_hist->Fill(planeNum, dit->getCenterStripPosition());  
  }

  // recognized pattern histograms
  for (auto &ds : *patterns)
  {
    unsigned int rpId = ds.detId();
    PotPlots &pp = potPlots[rpId];

    // count U and V patterns
    unsigned int u = 0, v = 0;
    for (auto &p : ds)
    {
      if (! p.getFittable())
        continue;

      if (p.getProjection() == TotemRPUVPattern::projU)
        u++;

      if (p.getProjection() == TotemRPUVPattern::projV)
        v++;
    }

//     pp.patterns_u->Fill(u);
//     pp.patterns_v->Fill(v);
  }

  // event-category histogram
  for (auto &it : potPlots)
  {
    const unsigned int &id = it.first;
    auto &pp = it.second;

    // process hit data for this plot
    unsigned int pl_u = planes_u[id].size();
    unsigned int pl_v = planes_v[id].size();

    // process pattern data for this pot
    const auto &rp_pat_it = patterns->find(id);

    unsigned int pat_u = 0, pat_v = 0;
    if (rp_pat_it != patterns->end())
    {
      for (auto &p : *rp_pat_it)
      {
        if (! p.getFittable())
          continue;
  
        if (p.getProjection() == TotemRPUVPattern::projU)
          pat_u++;
  
        if (p.getProjection() == TotemRPUVPattern::projV)
          pat_v++;
      }
    }

    // determine category
    signed int category = -1;

    if (pl_u == 0 && pl_v == 0) category = 0;   // empty
    
    if (category == -1 && pat_u + pat_v <= 1)
    {
      if (pl_u + pl_v < 6)
        category = 1;                           // insuff
      else
        category = 4;                           // shower
    }

    if (pat_u == 1 && pat_v == 1) category = 2; // 1-track

    if (category == -1) category = 3;           // multi-track

//     pp.event_category->Fill(category);
  }

  // RP track-fit plots
  for (auto &ds : *tracks)
  {
    unsigned int RPId = ds.detId();
    PotPlots &pp = potPlots[RPId];

    for (auto &ft : ds)
    {
      if (!ft.isValid())
        continue;
     
      // number of planes contributing to (valid) fits
      unsigned int n_pl_in_fit_u = 0, n_pl_in_fit_v = 0;
      for (auto &hds : ft.getHits())
      {
        unsigned int rawId = hds.detId();  
        unsigned int decId = TotemRPDetId::rawToDecId(rawId);
        bool uProj =TotemRPDetId::isStripsCoordinateUDirection(decId);

        for (auto &h : hds)
        {
          h.getPosition();  // just to keep compiler silent
          if (uProj)
            n_pl_in_fit_u++;
          else
            n_pl_in_fit_v++;
        }
      }

//       pp.h_planes_fit_u->Fill(n_pl_in_fit_u);
//       pp.h_planes_fit_v->Fill(n_pl_in_fit_v);
  
      // mean position of U and V planes
      double rp_x = ( geometry->GetDetector(TotemRPDetId::decToRawId(RPId*10 + 0))->translation().x() +
                      geometry->GetDetector(TotemRPDetId::decToRawId(RPId*10 + 1))->translation().x() ) / 2.;
      double rp_y = ( geometry->GetDetector(TotemRPDetId::decToRawId(RPId*10 + 0))->translation().y() +
                      geometry->GetDetector(TotemRPDetId::decToRawId(RPId*10 + 1))->translation().y() ) / 2.;
  
      // mean read-out direction of U and V planes
      CLHEP::Hep3Vector rod_U = geometry->LocalToGlobalDirection(TotemRPDetId::decToRawId(RPId*10 + 1), CLHEP::Hep3Vector(0., 1., 0.));
      CLHEP::Hep3Vector rod_V = geometry->LocalToGlobalDirection(TotemRPDetId::decToRawId(RPId*10 + 0), CLHEP::Hep3Vector(0., 1., 0.));
  
      double x = ft.getX0() - rp_x;
      double y = ft.getY0() - rp_y;
  
//       pp.trackHitsCumulativeHist->Fill(x, y);
  
      double U = x * rod_U.x() + y * rod_U.y();
      double V = x * rod_V.x() + y * rod_V.y();
  
//       pp.track_u_profile->Fill(U);
//       pp.track_v_profile->Fill(V);
    }
  }

  //------------------------------
  // Station Plots

  
  //------------------------------
  // Arm Plots
  {
    map<unsigned int, unsigned int> mTop, mHor, mBot;

    for (auto p : armPlots)
    {
      mTop[p.first] = 0;
      mHor[p.first] = 0;
      mBot[p.first] = 0;
    }

    for (auto &ds : *tracks)
    {
      unsigned int rpId = ds.detId();
      unsigned int armNum = rpId / 100;
      unsigned int rpNum = rpId % 10;

      for (auto &tr : ds)
      {
        if (! tr.isValid())
          continue;
  
        if (rpNum == 0 || rpNum == 4)
          mTop[armNum]++;
        if (rpNum == 2 || rpNum == 3)
          mHor[armNum]++;
        if (rpNum == 1 || rpNum == 5)
          mBot[armNum]++;
      }
    }

    for (auto &p : armPlots)
    {
//       p.second.h_numRPWithTrack_top->Fill(mTop[p.first]);
      p.second.h_numRPWithTrack->Fill(mHor[p.first]);
//       p.second.h_numRPWithTrack_bot->Fill(mBot[p.first]);
    }

    // track RP correlation
    for (auto &ds1 : *tracks)
    {
      for (auto &tr1 : ds1)
      {
        if (! tr1.isValid())
          continue;
  
        unsigned int rpId1 = ds1.detId();
        unsigned int arm1 = rpId1 / 100;
        unsigned int stNum1 = (rpId1 / 10) % 10;
        unsigned int rpNum1 = rpId1 % 10;
        unsigned int idx1 = stNum1/2 * 7 + rpNum1;
        bool hor1 = (rpNum1 == 2 || rpNum1 == 3);
  
        ArmPlots &ap = armPlots[arm1];
  
        for (auto &ds2 : *tracks)
        {
          for (auto &tr2 : ds2)
          {
            if (! tr2.isValid())
              continue;
          
            unsigned int rpId2 = ds2.detId();
            unsigned int arm2 = rpId2 / 100;
            unsigned int stNum2 = (rpId2 / 10) % 10;
            unsigned int rpNum2 = rpId2 % 10;
            unsigned int idx2 = stNum2/2 * 7 + rpNum2;
            bool hor2 = (rpNum2 == 2 || rpNum2 == 3);
    
            if (arm1 != arm2)
              continue;
    
//             ap.h_trackCorr->Fill(idx1, idx2); 
            
//             if (hor1 != hor2)
//               ap.h_trackCorr_overlap->Fill(idx1, idx2); 
          }
        }
      }
    }
  } */
}

//----------------------------------------------------------------------------------------------------

void TotemRPDiamondDQMSource::endLuminosityBlock(edm::LuminosityBlock const& lumiSeg, edm::EventSetup const& eSetup) 
{
}

//----------------------------------------------------------------------------------------------------

void TotemRPDiamondDQMSource::endRun(edm::Run const& run, edm::EventSetup const& eSetup)
{
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRPDiamondDQMSource);
