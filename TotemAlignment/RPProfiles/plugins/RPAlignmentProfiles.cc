/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: RPAlignmentProfiles.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "Geometry/TotemRecords/interface/RealGeometryRecord.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "TotemAlignment/RPTrackBased/interface/LocalTrackFitter.h"
#include "TotemAlignment/RPTrackBased/interface/AlignmentGeometry.h"
#include "TotemAlignment/RPTrackBased/interface/AlignmentTask.h"

#include <map>
#include <fstream>

#include "TGraph.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TF1.h"
#include "TVirtualFitter.h"

//#define DEBUG 1
  

/**
 *\brief Builds alignment-related hit profiles for RPs, units and stations.
 **/
class RPAlignmentProfiles : public edm::EDAnalyzer
{
  public:
    struct LinearFit
    {
      double S1, Sz, Szz, Sx, Sy, Szx, Szy;
      LinearFit() : S1(0.), Sz(0.), Szz(0.), Sx(0.), Sy(0.), Szx(0.), Szy(0.) {}
      void Fill(double x, double y, double z); 
      double Det() const { return Szz*S1 - Sz*Sz; }
      double GetAx() const { return (S1*Szx - Sz*Sx) / Det(); }
      double GetAy() const { return (S1*Szy - Sz*Sy) / Det(); }
    };

    /// set of profiles corresponding to a RP
    struct RPProfileSet {
      TGraph *hits;
      TH1D *hori, *vert;
      TH2D *prof;
      unsigned int x_slices_number;
      double x_slices_start, x_slice_width;
      unsigned int y_slices_number;
      double y_slices_start, y_slice_width;

      std::vector<TH1D*> vertSlices, horSlices;
      RPProfileSet() {}
      RPProfileSet(unsigned int id, const edm::ParameterSet &ps);
      void Fill(double x, double y);
    };

    /// set of profiles corresponding to a RP
    struct UnitProfileSet {
      unsigned int x_slices_number;
      double x_slices_start, x_slice_width;
      TGraph *hits;
      TH1D *vertOfVert;   ///< vertical profile from vertical detectors
      std::vector<TH1D*> vertSlices;
      UnitProfileSet() {}
      UnitProfileSet(unsigned int id, const edm::ParameterSet &ps);
      void Fill(unsigned int rpId, double x, double y);
    };

    /// set of profiles corresponding to a station
    struct StationProfileSet
    {
      TH1D *ax_hist, *ay_hist;
      StationProfileSet() {}
      StationProfileSet(unsigned int id, const edm::ParameterSet &ps);
      void Fill(const LinearFit &);
      void Fill(const LocalTrackFit &);
    };

    RPAlignmentProfiles(const edm::ParameterSet &ps);
    ~RPAlignmentProfiles() {}

  private:
    edm::ParameterSet parameterSet;

    /// verbosity level
    unsigned int verbosity;

    /// choice of pattern-recognition algorithm
    edm::InputTag RPTrackCandidateCollectionLabel;
    
    /// TODO
    edm::InputTag RPFittedTrackCollectionLabel;

    /// the numbers (=last digits) of RPs that will be processed
    std::vector<unsigned int> rpNums;

    /// will accept LocalTrackFitter results with two units only
    bool requireTwoUnits;

    /// output ROOT file anem
    std::string outputFile;

    /// XML file with fit results
    std::string resultFile;
    
    /// threshold-to-maximum ratios for FitPeak method
    double thresholdToMax_main, thresholdToMax_secondary;

    /// station fitter
    LocalTrackFitter fitter;

    /// a flag whether the module has already been initialized (after the first call of beginRun)
    bool initialized;

    /// map: station ID --> alignment geometry
    std::map<unsigned int, AlignmentGeometry> alGeometries;

    /// map: RPId -> set of rp profiles
    std::map<unsigned int, RPProfileSet> rp_1rp_profiles, rp_st_profiles;
    
    /// map: UnitId -> set of unit profiles
    std::map<unsigned int, UnitProfileSet> un_1rp_profiles, un_st_profiles;

    /// map: station ID -> set of profiles
    std::map<unsigned int, StationProfileSet> st_1rp_profiles, st_st_profiles;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event &e, const edm::EventSetup &es);
    virtual void endJob();

    void SaveResultSet(std::map<unsigned int, RPProfileSet>&, 
      std::map<unsigned int, UnitProfileSet> &, std::map<unsigned int, StationProfileSet> &,
      std::fstream &);

    /// analyzes vertical distributions from vertical pots
    void AnalyzeVertOfVert(TH1D *h, std::fstream &rf);

    // static members seem to be the only way to accomplish the fit...
    static TH1D *histToFit;     ///< the histogram to be fitted
    static int gapLeft;         ///< the first (leftmost) bin that is skipped
    static int gapRight;        ///< the last (rightmost) bin that is skipped

    /// user fitting function that skips region between gapLeft and gapRight
    static void MyFitChisquare(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag);

    /// finds and fits the main peak of a distribution
    void FitPeak(TH1D *h, double &pos, double &err, double thresholdToMax);
};

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::LinearFit::Fill(double x, double y, double z)
{
  S1++;
  Sz += z;
  Szz += z*z;
  Sx += x;
  Sy += y;
  Szx += z*x;
  Szy += z*y;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPAlignmentProfiles::RPProfileSet::RPProfileSet(unsigned int id, const edm::ParameterSet &ps)
{
  const ParameterSet &p = ps.getParameterSet("hit_profiles");

  char buf[20];
  hits = new TGraph();  // coordinates in mm
  sprintf(buf, "RP%u: hitmap", id); hits->SetName(buf);
  sprintf(buf, "RP%u: horizontal", id);
  hori = new TH1D(buf, "; x   (mm)", p.getParameter<unsigned int>("bins_1d"),
    p.getParameter<double>("x_min"), p.getParameter<double>("x_max"));

  sprintf(buf, "RP%u: vertical", id);
  vert = new TH1D(buf, "; y   (mm)", p.getParameter<unsigned int>("bins_1d"),
    p.getParameter<double>("y_min"), p.getParameter<double>("y_max"));

  sprintf(buf, "RP%u: 2D profile", id);
  prof = new TH2D(buf, "; x   (mm); y   (mm)",
    p.getParameter<unsigned int>("bins_2d"), p.getParameter<double>("x_min"), p.getParameter<double>("x_max"),
    p.getParameter<unsigned int>("bins_2d"), p.getParameter<double>("y_min"), p.getParameter<double>("y_max")); 

  
  x_slices_start = p.getParameter<double>("x_slices_start");
  x_slice_width = p.getParameter<double>("x_slice_width");
  x_slices_number = p.getParameter<unsigned int>("x_slices_number");
  y_slices_start = p.getParameter<double>("y_slices_start");
  y_slice_width = p.getParameter<double>("y_slice_width");
  y_slices_number = p.getParameter<unsigned int>("y_slices_number");

  for (unsigned int i = 0; i < x_slices_number; i++) {
    sprintf(buf, "RP%u: vertical, %.2f <= x < %.2f", id, x_slices_start+i*x_slice_width,
      x_slices_start+(i+1)*x_slice_width);
    TH1D *vs = new TH1D(buf, "; y   (mm)", p.getParameter<unsigned int>("bins_1d"),
      p.getParameter<double>("y_min"), p.getParameter<double>("y_max"));
    vertSlices.push_back(vs);
  }
   
  for (unsigned int i = 0; i < y_slices_number; i++) {
    sprintf(buf, "RP%u: horizontal, %.2f <= y < %.2f", id, y_slices_start+i*y_slice_width,
      y_slices_start+(i+1)*y_slice_width);
    TH1D *vs = new TH1D(buf, "; y   (mm)", p.getParameter<unsigned int>("bins_1d"),
      p.getParameter<double>("x_min"), p.getParameter<double>("x_max"));
    horSlices.push_back(vs);
  }
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::RPProfileSet::Fill(double x, double y)
{
  hits->SetPoint(hits->GetN(), x, y);
  hori->Fill(x);
  vert->Fill(y);
  prof->Fill(x, y);

  unsigned int x_idx = (unsigned int) ((x - x_slices_start) / x_slice_width);
  if (x_idx < x_slices_number)
    vertSlices[x_idx]->Fill(y);
  
  unsigned int y_idx = (unsigned int) ((y - y_slices_start) / y_slice_width);
  if (y_idx < y_slices_number)
    horSlices[y_idx]->Fill(x);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPAlignmentProfiles::UnitProfileSet::UnitProfileSet(unsigned int id, const edm::ParameterSet &ps)
{
  const ParameterSet &p = ps.getParameterSet("vertOfVert");

  char buf[20];
  hits = new TGraph();  // coordinates in mm
  sprintf(buf, "Unit%u: hitmap", id); hits->SetName(buf);
  sprintf(buf, "Unit%u: vertOfVert", id);
  vertOfVert = new TH1D(buf, "; y   (mm)", p.getParameter<unsigned int>("bins_1d"),
    p.getParameter<double>("y_min"), p.getParameter<double>("y_max"));

  x_slices_start = p.getParameter<double>("x_slices_start");
  x_slice_width = p.getParameter<double>("x_slice_width");
  x_slices_number = p.getParameter<unsigned int>("x_slices_number");

  for (unsigned int i = 0; i < x_slices_number; i++) {
    sprintf(buf, "Unit%u: vertOfVert, %.2f <= x < %.2f", id, x_slices_start+i*x_slice_width,
      x_slices_start+(i+1)*x_slice_width);
    TH1D *vs = new TH1D(buf, "; y   (mm)", p.getParameter<unsigned int>("bins_1d"),
      p.getParameter<double>("y_min"), p.getParameter<double>("y_max"));
    vertSlices.push_back(vs);
  }
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::UnitProfileSet::Fill(unsigned int rpID, double x, double y)
{
  hits->SetPoint(hits->GetN(), x, y);

  if (rpID % 10 != 2 && rpID % 10 != 3)
    vertOfVert->Fill(y);

  unsigned int x_idx = (unsigned int) ((x - x_slices_start) / x_slice_width);
  if (x_idx < x_slices_number)
    vertSlices[x_idx]->Fill(y);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPAlignmentProfiles::StationProfileSet::StationProfileSet(unsigned int id, const edm::ParameterSet &ps)
{
  const ParameterSet &p = ps.getParameterSet("angular_profiles");

  char buf[20];
  sprintf(buf, "station %u: ax_hist", id);
  ax_hist = new TH1D(buf, ";a_{x}   (mrad)", p.getParameter<unsigned int>("bins_1d"),
    p.getParameter<double>("x_min"), p.getParameter<double>("x_max"));
      
  sprintf(buf, "station %u: ay_hist", id);
  ay_hist = new TH1D(buf, ";a_{y}   (mrad)", p.getParameter<unsigned int>("bins_1d"),
    p.getParameter<double>("y_min"), p.getParameter<double>("y_max"));
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::StationProfileSet::Fill(const LinearFit &f)
{
  double det = f.Det();

#ifdef DEBUG
  printf("* RPStationAngularProfiles fit\n\tS1 = %.0f, det = %.1E\n", f.S1, det);
#endif

  if (det == 0.)
    return;

  double ax = f.GetAx();
  double ay = f.GetAy();

  ax_hist->Fill(ax*1E3);
  ay_hist->Fill(ay*1E3);

#ifdef DEBUG
  printf("\tax = %E, ay = %E (rad)\n", ax, ay);
#endif
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::StationProfileSet::Fill(const LocalTrackFit &f)
{
  ax_hist->Fill(f.ax*1E3);
  ay_hist->Fill(f.ay*1E3);
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

RPAlignmentProfiles::RPAlignmentProfiles(const edm::ParameterSet &ps) :
  parameterSet(ps),
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  rpNums(ps.getParameter< vector<unsigned int> >("rpNums")),
  requireTwoUnits(ps.getParameter<bool>("requireTwoUnits")),
  outputFile(ps.getParameter<string>("outputFile")),
  resultFile(ps.getParameter<string>("resultFile")),
  thresholdToMax_main(ps.getParameter<double>("thresholdToMax_main")),
  thresholdToMax_secondary(ps.getParameter<double>("thresholdToMax_secondary")),
  fitter(ps),
  initialized(false)
{
	RPFittedTrackCollectionLabel = ps.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
	RPTrackCandidateCollectionLabel = ps.getParameter<edm::InputTag>("RPTrackCandidateCollectionLabel");
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  if (initialized)
    return;

  // prepare AlignmentGeometry per station
  ESHandle<TotemRPGeometry> geom;
  es.get<RealGeometryRecord>().get(geom);

  for (unsigned int a = 0; a < 2; a++)
    for (unsigned int s = 0; s < 3; s++) {
      if (s == 1)
        continue;
      unsigned int stId = a*10 + s;

      vector<unsigned int> rps;
      for (unsigned int r = 0; r < 6; r++)
        rps.push_back(stId*10 + r);

      AlignmentTask::BuildGeometry(rps, geom.product(), 0., alGeometries[stId]);
    }

  initialized = true;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  ESHandle<TotemRPGeometry> geom;
  eSetup.get<RealGeometryRecord>().get(geom);
  
  Handle< RPTrackCandidateCollection > trackColl;
  event.getByLabel(RPTrackCandidateCollectionLabel, trackColl);

  edm::Handle< RPFittedTrackCollection > fitTrCol;
  event.getByLabel(RPFittedTrackCollectionLabel, fitTrCol);

  //----------------------- 1-RP-fit profiles ---------------------------

  // RP and Unit profiles
  for (std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol->begin(); it != fitTrCol->end(); ++it) {
    if (!it->second.IsValid())
      continue;

    // RP profiles
    map<unsigned int, RPProfileSet>::iterator pit = rp_1rp_profiles.find(it->first);
    if (pit == rp_1rp_profiles.end()) 
      pit = rp_1rp_profiles.insert(pair<unsigned int, RPProfileSet>(it->first, RPProfileSet(it->first, parameterSet))).first;
    
    pit->second.Fill(it->second.X0(), it->second.Y0());

    // unit profiles
    unsigned int stId = it->first / 10;
    unsigned int rpNum = it->first % 10;
    unsigned int unId = stId * 10;
    if (rpNum > 2)
      unId++;

    map<unsigned int, UnitProfileSet>::iterator uit = un_1rp_profiles.find(unId);
    if (uit == un_1rp_profiles.end()) 
      uit = un_1rp_profiles.insert(pair<unsigned int, UnitProfileSet>(unId, UnitProfileSet(unId, parameterSet))).first;
      
      // TODO: without this the code does not work ????
      char buf[20];
      sprintf(buf, "%u, %u", unId, uit->first);

    if (verbosity > 5)
      printf("* event %u, fit in RP%u: X0 = %.1f um, Y0 = %.1f um\n", event.id().event(), it->first, it->second.X0()*1E3, it->second.Y0()*1E3);

    uit->second.Fill(it->first, it->second.X0(), it->second.Y0());
  }

  // fit all stations
  map<unsigned int, LinearFit> fits;  // map: station -> fit
  for (std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol->begin(); it != fitTrCol->end(); ++it) {
    if (!it->second.IsValid())
      continue;

    double z = geom->GetRPDevice(it->first)->translation().z();   // mm
    double x = it->second.X0(), y = it->second.Y0();              // mm

    fits[TotRPDetId::StOfRP(it->first)].Fill(x, y, z);
  }

  // fill fit results in station st_1rp_profiles
  for (map<unsigned int, LinearFit>::iterator it = fits.begin(); it != fits.end(); ++it) {
    map<unsigned int, StationProfileSet>::iterator pit = st_1rp_profiles.find(it->first);
    if (pit == st_1rp_profiles.end())
      pit = st_1rp_profiles.insert(pair<unsigned int, StationProfileSet>(it->first, StationProfileSet(it->first, parameterSet))).first;
    pit->second.Fill(it->second);
  }

  //----------------------- station-fit profiles ---------------------------

  // build hit collections per station
  map<unsigned int, HitCollection> hits;

  for (RPTrackCandidateCollection::const_iterator it = trackColl->begin(); it != trackColl->end(); ++it) {
    // skip non fittable candidates
    if (!it->second.Fittable())
      continue;

    unsigned int stId = it->first/10;
    unsigned int rpNum = it->first % 10;

    // skip if RP not selected by user
    if (find(rpNums.begin(), rpNums.end(), rpNum) == rpNums.end())
      continue;

    const vector<RPRecoHit> &cHits = it->second.TrackRecoHits();
    for (unsigned int i = 0; i < cHits.size(); i++)
      hits[stId].push_back(cHits[i]);
  }

  for (map<unsigned int, HitCollection>::iterator it = hits.begin(); it != hits.end(); ++it) {
    // fit hits per station
    LocalTrackFit ltf;
    if (! fitter.Fit(it->second, alGeometries[it->first], ltf) )
      continue;

    // get list of selected RPs and units
    set<unsigned int> selectedRPs;
    set<unsigned int> selectedUnits;
    for (HitCollection::iterator hit = it->second.begin(); hit != it->second.end(); ++hit) {
      unsigned int rpId = hit->id / 10;
      unsigned int rpNum = rpId % 10;
      unsigned int unId = rpId - rpNum;
      if (rpNum > 2)
        unId++;
      selectedRPs.insert(rpId);
      selectedUnits.insert(unId);
    }

    // has two units?
    if (requireTwoUnits && selectedUnits.size() < 2)
      return;
    
    // fill track angles
    map<unsigned int, StationProfileSet>::iterator sit = st_st_profiles.find(it->first);
    if (sit == st_st_profiles.end())
      sit = st_st_profiles.insert(pair<unsigned int, StationProfileSet>(it->first, StationProfileSet(it->first, parameterSet))).first;
    sit->second.Fill(ltf);

    // fill hit positions in all RPs
    for (set<unsigned int>::iterator rit = selectedRPs.begin(); rit != selectedRPs.end(); ++rit) {
      double z = geom->GetRPDevice(*rit)->translation().z(); 
      double x = 0., y = 0.;
      ltf.Eval(z, x, y);

      // RP plots
      map<unsigned int, RPProfileSet>::iterator pit = rp_st_profiles.find(*rit);
      if (pit == rp_st_profiles.end()) 
        pit = rp_st_profiles.insert(pair<unsigned int, RPProfileSet>(*rit, RPProfileSet(*rit, parameterSet))).first;
      pit->second.Fill(x, y);
      
      // unit plots
      unsigned int rpNum = *rit % 10;
      unsigned int unId = *rit - rpNum;
      if (rpNum > 2)
        unId++;

      map<unsigned int, UnitProfileSet>::iterator uit = un_st_profiles.find(unId);
      if (uit == un_st_profiles.end()) 
        uit = un_st_profiles.insert(pair<unsigned int, UnitProfileSet>(unId, UnitProfileSet(unId, parameterSet))).first;

      // TODO: without this the code does not work ????
      char buf[20];
      sprintf(buf, "%u, %u", unId, uit->first);

      uit->second.Fill(*rit, x, y);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::endJob()
{
  if (outputFile.empty())
    return;

  TFile *f = new TFile(outputFile.c_str(), "recreate");
  if (f->IsZombie())
    return;

  fstream rf;
  if (!resultFile.empty())
    rf.open(resultFile.c_str(), ios_base::out);

  if (!rf.is_open())
    return;

  rf << "<xml DocumentType=\"ProfileAlignmentResult\">" << endl;

  gDirectory = f->mkdir("1-RP fits");
  rf << "  <results fit_type=\"1-RP fits\">" << endl;
  SaveResultSet(rp_1rp_profiles, un_1rp_profiles, st_1rp_profiles, rf);
  rf << "  </results>" << endl;
  
  gDirectory = f->mkdir("station fits");
  rf << "  <results fit_type=\"station fits\">" << endl;
  SaveResultSet(rp_st_profiles, un_st_profiles, st_st_profiles, rf);
  rf << "  </results>" << endl;
 
  rf << "</xml>" << endl;

  delete f;
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::SaveResultSet(std::map<unsigned int, RPProfileSet> &rp_profiles, 
      std::map<unsigned int, UnitProfileSet> &un_profiles, 
      std::map<unsigned int, StationProfileSet> &st_profiles,
      std::fstream &rf)
{
  TDirectory *topDir = gDirectory;

  for (map<unsigned int, RPProfileSet>::iterator it = rp_profiles.begin(); it != rp_profiles.end(); ++it) {
    unsigned int rpIdx = it->first % 10;
    bool horizontal = (rpIdx == 2 || rpIdx == 3);
    double v, ve, h, he;
    FitPeak(it->second.hori, h, he, (horizontal) ? thresholdToMax_main : thresholdToMax_secondary);
    FitPeak(it->second.vert, v, ve, (horizontal) ? thresholdToMax_main : thresholdToMax_secondary);

    char buf[20];
    sprintf(buf, "RP %u", it->first);
    gDirectory = topDir->mkdir(buf);
    it->second.hits->Write("hit map");
    it->second.hori->Write();

    gDirectory = gDirectory->mkdir("horizontal slices");
    vector<TH1D*> &horSlices = it->second.horSlices;
    for (vector<TH1D*>::iterator sit = horSlices.begin(); sit != horSlices.end(); ++sit) {
      double p, pe;
      FitPeak(*sit, p, pe, thresholdToMax_secondary);
      (*sit)->Write();
    }
    gDirectory->cd("..");
    
    it->second.vert->Write();

    gDirectory = gDirectory->mkdir("vertical slices");
    vector<TH1D*> &vertSlices = it->second.vertSlices;
    for (vector<TH1D*>::iterator sit = vertSlices.begin(); sit != vertSlices.end(); ++sit) {
      double p, pe;
      FitPeak(*sit, p, pe, thresholdToMax_secondary);
      (*sit)->Write();
    }
    gDirectory->cd("..");

    it->second.prof->Write();

    rf << "    <rp id=\"" << it->first << 
      "\" peak_x=\"" << h*1E3 << "\" peak_x_error=\"" << he*1E3 << 
      "\" peak_y=\"" << v*1E3 << "\" peak_y_error=\"" << ve*1E3 << 
      "\"/>" << endl;
  }
  
  rf << endl;

  for (map<unsigned int, UnitProfileSet>::iterator it = un_profiles.begin(); it != un_profiles.end(); ++it) {
    printf("** %u\n", it->first);
    char buf[20];
    sprintf(buf, "Unit %u", it->first);
    gDirectory = topDir->mkdir(buf);
    it->second.hits->Write("hit map");
    rf << "    <unit id=\"" << it->first << "\" ";
    AnalyzeVertOfVert(it->second.vertOfVert, rf);
    rf << "/>" << endl;
    it->second.vertOfVert->Write();
    
    gDirectory = gDirectory->mkdir("vertical slices");
    vector<TH1D*> &vertSlices = it->second.vertSlices;
    for (vector<TH1D*>::iterator sit = vertSlices.begin(); sit != vertSlices.end(); ++sit) {
      AnalyzeVertOfVert(*sit, rf);
      (*sit)->Write();
    }
  }

  rf << endl;

  for (map<unsigned int, StationProfileSet>::iterator it = st_profiles.begin(); it != st_profiles.end(); ++it) {
    char buf[20];
    sprintf(buf, "Station %u", it->first);
    gDirectory = topDir->mkdir(buf);
    it->second.ax_hist->Write();
    it->second.ay_hist->Write();

    rf << "    <station id=\"" << it->first << 
      "\" mean_ax=\"" << it->second.ax_hist->GetMean() << "\" mean_ax_error=\"" << it->second.ax_hist->GetMeanError() << 
      "\" mean_ay=\"" << it->second.ay_hist->GetMean() << "\" mean_ay_error=\"" << it->second.ay_hist->GetMeanError() << 
      "\"/>" << endl;
  }
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::AnalyzeVertOfVert(TH1D *h, fstream &rf)
{
  if (verbosity)
    printf(">> RPAlignmentProfiles::AnalyzeVertOfVert > analyzing histogram %s\n", h->GetName());

  // find the largest gap
  pair<int, int> theGap(0, 0);
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    if (h->GetBinContent(i) == 0) {
      int b = i;
      while (i <= h->GetNbinsX() && h->GetBinContent(i) == 0)
        i++;
      int e = i - 1;

      //printf("\t\tgap: %i to %i\n", b, e);

      // skip edge gaps
      if (b == 1 || e == h->GetNbinsX())
        continue;

      int size = theGap.second - theGap.first;
      if (e - b > size)
        theGap = pair<int, int>(b, e);
    }
  }

  // prepare the fitter
  TVirtualFitter::Fitter(h)->SetFCN(MyFitChisquare);
  histToFit = h;
  gapLeft = theGap.first - 2; 
  gapRight = theGap.second + 2;

  if (verbosity)
    printf("\tfit with gap between bins %i to %i (included), i.e. from %E mm to %E mm\n", gapLeft, gapRight, 
        h->GetBinCenter(gapLeft), h->GetBinCenter(gapRight));

  TF1 *fitFun = new TF1("fitFun", "gaus(0)");
  fitFun->SetLineColor(2);
  double a = h->GetBinContent(h->GetMaximumBin())*0.8;
  double mu = h->GetMean();
  double si = h->GetRMS();

  if (verbosity)
    printf("\tpre-fit estimate: a = %.2E, mu = %.3f mm, sigma = %.3f mm\n", a, mu, si);

  fitFun->SetParameters(a, mu, si);
  h->Fit(fitFun, "uq");

  a = fitFun->GetParameter(0);
  mu = fitFun->GetParameter(1);
  si = fitFun->GetParameter(2);

  if (verbosity)
    printf("\tfit results: a = %.2E, mu = %.3f mm, sigma = %.3f mm\n", a, mu, si);

  if (rf.is_open())
    rf << "center_y=\"" << (mu*1E3) << "\" y_dist_sigma=\"" << (si*1E3) << "\"";
}

//----------------------------------------------------------------------------------------------------

TH1D* RPAlignmentProfiles::histToFit = NULL;
int RPAlignmentProfiles::gapLeft = -1, RPAlignmentProfiles::gapRight = -1;

void RPAlignmentProfiles::MyFitChisquare(Int_t &npar, Double_t *gin, Double_t &f, Double_t *u, Int_t flag)
{
  double &a = u[0];
  double &mu = u[1];
  double sisq = u[2]; sisq *= sisq;

  f = 0.;
  for (int i = 1; i <= histToFit->GetNbinsX(); i++) {
    if (i == gapLeft) {
      //printf("* fit of %s: jumping from %i to %i\n", histToFit->GetName(), gapLeft, gapRight);
      i = gapRight;
      continue;
    }
    double x = histToFit->GetBinCenter(i);
    double c = histToFit->GetBinContent(i);
    double fun = x - mu;
    fun = c - a * exp(-fun*fun/2./sisq);
    f += fun*fun;
  }
}

//----------------------------------------------------------------------------------------------------

void RPAlignmentProfiles::FitPeak(TH1D *h, double &pos, double &err, double thresholdToMax)
{
  if (verbosity)
    printf(">> RPAlignmentProfiles::FitPeak > analyzing histogram %s\n", h->GetName());

  // defaults
  pos = err = 0.;

  // find maximum
  double max = 0.;
  int maxBin = 0;
  for (int i = 1; i <= h->GetNbinsX(); i++) {
    if (h->GetBinContent(i) > max) {
      max = h->GetBinContent(i);
      maxBin = i;
    }
  }
  
  if (verbosity)
    printf("\tmaximum: %.3f at bin %u (%.3f)\n", max, maxBin, h->GetBinCenter(maxBin));

  if (max == 0.)
    return;

  double th = thresholdToMax*max;

  if (verbosity)
    printf("\tthresholdToMax: %.3f, threshold: %.3f\n", thresholdToMax, th);

  int leftBin = 0;
  for (int i = maxBin; i > 0; --i)
    if (h->GetBinContent(i) < th) {
      leftBin = i;
      break;
    }

  int rightBin = 0;
  for (int i = maxBin; i <= h->GetNbinsX(); ++i)
    if (h->GetBinContent(i) < th) {
      rightBin = i;
      break;
    }

  double left = h->GetBinCenter(leftBin);
  double right = h->GetBinCenter(rightBin);
  
  if (verbosity) {
    printf("\tleftBin: %u (%.3f)\n", leftBin, left);
    printf("\trightBin: %u (%.3f)\n", rightBin, right);
  }

  TF1 *fitFun = new TF1("fitFun", "gaus(0)");
  fitFun->SetLineColor(2);
  fitFun->SetParameters(max, h->GetBinCenter(maxBin), right-left);
  h->Fit(fitFun, "q", "", left, right);

  double a = fitFun->GetParameter(0);
  double mu = fitFun->GetParameter(1);
  double si = fitFun->GetParameter(2);

  if (verbosity)
    printf("\tfit results: a = %.2E, mu = %.3f mm, sigma = %.3f mm\n", a, mu, si);

  pos = mu;
  err = fitFun->GetParError(1); 
}

DEFINE_FWK_MODULE(RPAlignmentProfiles);

