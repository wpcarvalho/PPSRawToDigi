/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*  Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/ProtonTransportFunctions/interface/ProtonTransportFunctions.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "TotemCondFormats/DataRecord/interface/ProtonTransportRcd.h"
#include "SimG4CMS/TotemRPProtTranspPar/interface/LHCOpticsApproximator.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidate.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPTrackCandidateCollection.h"
#include "Geometry/TotemRecords/interface/MisalignedGeometryRecord.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"

#include "TGraph.h"
#include "TH2D.h"
#include "TFile.h"

//----------------------------------------------------------------------------------------------------

/**
 *\brief Fast (no G4) proton simulation from IP to all RPs at the 200 - 220m regions.
 * This is a temporary module, before the updated geometry will be implemented.
 */
class RPFastFullSimulation200 : public edm::EDProducer
{
  public:
    RPFastFullSimulation200(const edm::ParameterSet &);
    virtual ~RPFastFullSimulation200();

  protected:
    unsigned int verbosity;
    string hepMCLabel;

    /// strip pitch in mm
    double pitch;

    /// widt of double-strip cluster region, in mm
    double dscrWidth;

    /// position (value of v) where double-strip cluster region starts, in mm
    double dscrStart;

    /// position (value of v) where double-strip cluster region eds, in mm
    double dscrEnd;

    /// whether measurement values shall be rounded to the nearest strip
    bool roundToPitch;

    /// variance of box distribution on (-P/2, +P/2)
    double sigma0;

    /// angular limit to distinguish between forward, central and backward particles
    double thetaLim;

    /// v position of strip 0, in mm
    double stripZeroPosition;

    /// minimum number of U (or V) sensors required to consider the RP active
    unsigned int minUVSensorsPerRP;

    /// minimum number of active RPs to consider the track "reconstructable"
    unsigned int minRPsPerStation;

    /// IDs of RPs which will be analysed, if empty no filter is assumed
    vector<unsigned int> RPsAllowed;

    bool savePlots;
    string plotFileName;
    double plot_z;

    TGraph *g_scoringPlane;
    TH2D *h_scoringPlane;
    TH2D *h_ax_ay;

    bool dumpTopRPProtons, dumpBottomRPProtons, dumpHorizontalRPProtons;

    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void produce(edm::Event &, const edm::EventSetup&);
};


using namespace edm;
using namespace std;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

RPFastFullSimulation200::RPFastFullSimulation200(const ParameterSet &ps) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
    hepMCLabel(ps.getParameter<string>("hepMCLabel")),
    pitch(ps.getParameter<double>("pitch")),
    dscrWidth(ps.getParameter<double>("dscrWidth")),
    dscrStart(pitch/2. - dscrWidth/2.),
    dscrEnd(pitch/2. + dscrWidth/2.),
    roundToPitch(ps.getParameter<bool>("roundToPitch")),
    sigma0(pitch/sqrt(12)),
    thetaLim(ps.getParameter<double>("thetaLimit")),
    minUVSensorsPerRP(ps.getParameter<unsigned int>("minUVSensorsPerRP")),
    minRPsPerStation(ps.getParameter<unsigned int>("minRPsPerStation")),
    RPsAllowed(ps.getParameter<vector<unsigned int> >("RPsAllowed")),

    savePlots(ps.getParameter<bool>("savePlots")),
    plotFileName(ps.getParameter<string>("plotFileName")),
    plot_z(ps.getParameter<double>("plot_z")),

    g_scoringPlane(NULL),
    h_scoringPlane(NULL),
    h_ax_ay(NULL),

    dumpTopRPProtons(ps.getParameter<bool>("dumpTopRPProtons")),
    dumpBottomRPProtons(ps.getParameter<bool>("dumpBottomRPProtons")),
    dumpHorizontalRPProtons(ps.getParameter<bool>("dumpHorizontalRPProtons"))
{
  // v position of strip 0
  stripZeroPosition = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_
    - RPTopology::y_width_/2.;

  if (savePlots)
  {
    g_scoringPlane = new TGraph();
    g_scoringPlane->SetName("g_scoringPlane");
    g_scoringPlane->SetTitle(";x   (mm);y   (mm)");

    h_scoringPlane = new TH2D("h_scoringPlane", ";x   (mm);y   (mm)", 100, 0., 0., 100, 0., 0.);

    h_ax_ay = new TH2D("h_ax_ay", ";ax (urad);ay (urad)", 100, 0., 0., 100, 0., 0.);
  }

  produces<RPTrackCandidateCollection>();
  produces< DetSetVector<RPRecoHit> >();
  produces< RPStationTrackFitCollection >();
}

//----------------------------------------------------------------------------------------------------

RPFastFullSimulation200::~RPFastFullSimulation200()
{
  if (savePlots)
  {
    TFile *f_out = new TFile(plotFileName.c_str(), "recreate");

    g_scoringPlane->Write();
    h_scoringPlane->Write();
    h_ax_ay->Write();

    delete f_out;
  }
}

//----------------------------------------------------------------------------------------------------

void RPFastFullSimulation200::beginRun(edm::Run&, edm::EventSetup const& es)
{
  ESHandle<TotemRPGeometry> geom;
  es.get<MisalignedGeometryRecord>().get(geom);
}

//----------------------------------------------------------------------------------------------------

void RPFastFullSimulation200::produce(edm::Event &event, const EventSetup &es)
{
  if (verbosity > 5)
    printf("\n>> RPFastFullSimulation200::produce > event %i\n", event.id().event());

  // get geometry and optics
  ESHandle<TotemRPGeometry> geom;
  es.get<MisalignedGeometryRecord>().get(geom);

  ESHandle<BeamOpticsParams> optPar;
  es.get<BeamOpticsParamsRcd>().get(optPar);

  ESHandle<ProtonTransportFunctions> optFun;
  es.get<ProtonTransportRcd>().get(optFun);

  Handle<HepMCProduct> mcProd;
  event.getByLabel(hepMCLabel, mcProd);
  const GenEvent *mcEvent = mcProd->GetEvent();

  // prepare output
  auto_ptr<RPTrackCandidateCollection> output(new RPTrackCandidateCollection());
  auto_ptr< DetSetVector<RPRecoHit> > recoColl(new DetSetVector<RPRecoHit>());
  auto_ptr< RPStationTrackFitCollection > trackColl(new RPStationTrackFitCollection());

  const FourVector *last_vertex = 0;

  // loop over vertices
  for (GenEvent::vertex_const_iterator vit = mcEvent->vertices_begin(); vit != mcEvent->vertices_end(); ++vit) {
    const FourVector &vertex = (*vit)->position(); // in mm
    // ignore vertex if a duplicate
    if (last_vertex != 0) {
        if (*last_vertex == vertex) {
            continue;
        }
    }
    last_vertex = &vertex;

    // loop over outgoing particles
    for (GenVertex::particles_out_const_iterator pit = (*vit)->particles_out_const_begin(); pit != (*vit)->particles_out_const_end(); ++pit) {
      // keep only protons
      if ((*pit)->pdg_id() != 2212)
        continue;

      const FourVector &p = (*pit)->momentum();

      // determine the arm of action
      double th = p.theta();
      unsigned int arm = 123;
      if (th < thetaLim) arm = 1;
      if (th > (M_PI - thetaLim)) arm = 0;

      // skip central particles
      if (arm > 1)
        continue;

      set<unsigned int> RPs, RPs_p;
      unsigned int ref_RP;

      if (arm == 0)
      {
        ref_RP = 20;
        RPs_p.insert(  0); RPs_p.insert(  1); RPs_p.insert(  2); RPs_p.insert(  3); RPs_p.insert(  4); RPs_p.insert(  5);
        RPs_p.insert( 20); RPs_p.insert( 21); RPs_p.insert( 22); RPs_p.insert( 23); RPs_p.insert( 24); RPs_p.insert( 25);
      } else {
        ref_RP = 120;
        RPs_p.insert(100); RPs_p.insert(101); RPs_p.insert(102); RPs_p.insert(103); RPs_p.insert(104); RPs_p.insert(105);
        RPs_p.insert(120); RPs_p.insert(121); RPs_p.insert(122); RPs_p.insert(123); RPs_p.insert(124); RPs_p.insert(125);
      }

      // leave only the RPs which are in the allowed RPs list
      // if an empty vector is passed assume no filtering
      if (RPsAllowed.size() > 0) {
          for (vector<unsigned int>::const_iterator rp_i = RPsAllowed.begin();
                  rp_i != RPsAllowed.end(); rp_i++) {
              if (RPs_p.find(*rp_i) != RPs_p.end()) {
                RPs.insert(*rp_i);
              }
          }
      } else {
          RPs.insert(RPs_p.begin(), RPs_p.end());
      }

      // center of the reference RP box, in mm
      double z0 = geom->GetRPDevice(ref_RP)->translation().z();

      // optical functions for the RP
      LHCOpticsApproximator *of = optFun->GetFunction(ref_RP);
      if (!of)
        throw cms::Exception("RPFastFullSimulation200::produce") << "No optical function found for ref_RP " << ref_RP;

      // get position and direction of the proton at the center of the RP box
      double xi = p.e() / optPar->GetBeamEnergy() - 1.;
      double parIn[5] = {vertex.x() * 1E-3, p.x()/p.rho(), vertex.y() * 1E-3, p.y()/p.rho(), xi};
      double parOut[5]; // in m, rad and xi dimensionless
      bool transportable = of->Transport(parIn, parOut, true);

      if (verbosity > 5)
        printf("ref_RP%u\tz0 = %E\t p.z = %E\t transportable = %i, beam type = %i\n", ref_RP, z0, p.z(), transportable, of->GetBeamType());

      if (verbosity > 8)
      {
        printf("param input:  %E\t%E\t%E\t%E\t%E\n", parIn[0], parIn[1], parIn[2], parIn[3], parIn[4]);
        printf("param output: %E\t%E\t%E\t%E\t%E\n", parOut[0], parOut[1], parOut[2], parOut[3], parOut[4]);
      }

      // do not continue if the proton is stopped
      if (!transportable)
        continue;

      double ax = parOut[1], bx = parOut[0]*1E3,
             ay = parOut[3], by = parOut[2]*1E3;

      // save hit position at the scoring plane
      if (savePlots && arm == 1)
      {
        double x = bx  + ax * (plot_z - z0); // in mm
        double y = by  + ay * (plot_z - z0); // in mm

        g_scoringPlane->SetPoint(g_scoringPlane->GetN(), x, y);
        h_scoringPlane->Fill(x, y);
        h_ax_ay->Fill(ax*1E6, ay*1E6);
      }

      // counter of RPs with hits
      unsigned int hitTopRPs = 0, hitBottomRPs = 0, hitHorizontalRPs = 0;

      // how many RPs were activated
      unsigned int activeRPs = 0;

      // loop over RPs
      for (set<unsigned int>::const_iterator rp = RPs.begin(); rp != RPs.end(); ++rp)
      {
        // collection of hits
        vector<RPRecoHit> hits;

        // u and v hit counter
        unsigned uHits = 0, vHits = 0;

        // loop over all detectors within the RP
        set<unsigned int> dets = geom->DetsInRP(*rp);
        for (set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit)
        {
          // calculate hit position in global coordinates
          unsigned int rawId = TotRPDetId::DecToRawId(*dit);
          DetSet<RPRecoHit> &hits_plane = recoColl->find_or_insert(rawId);
          double z = geom->GetDetector(rawId)->translation().z(); // in mm
          double x = bx  + ax * (z - z0); // in mm
          double y = by  + ay * (z - z0); // in mm

          if (verbosity > 5)
            printf("\tdet=%u: z=%E,x=%E,y=%E mm", *dit, z, x, y);

          // convert to local coordinates
          CLHEP::Hep3Vector hit(x, y, z);
          hit = geom->GlobalToLocal(rawId, hit);
          double u = hit.x(); // in mm
          double v = hit.y(); // in mm

          // is it within detector?
          if (!RPTopology::IsHit(u, v))
          {
            if (verbosity > 5)
              printf("\n");
            continue;
          }

          if (TotRPDetId::IsStripsCoordinateUDirection(*dit))
            uHits++;
          else
            vHits++;

          // round the measurement
          if (roundToPitch)
          {
            double m = stripZeroPosition - v;
            signed int strip = (int) floor(m / pitch);
            double offset = m - strip*pitch;
            double dstrip = 0.;

            if (offset < dscrStart)
            {
              dstrip = 0.;
            } else {
              if (offset < dscrEnd)
                dstrip = 0.5;
              else
                dstrip = 1.;
            }

            v = stripZeroPosition - pitch * (strip + dstrip);
            if (verbosity > 5)
              printf(" | stripZeroPosition = %+8.4f, strip = %+6.1f", stripZeroPosition, strip+dstrip);
          }

          if (verbosity > 5)
            printf(", v = %E\n", v);

          RPRecoHit rh = RPRecoHit(rawId, v, sigma0);
          hits.push_back(rh);
          hits_plane.push_back(rh);
        }

        if (verbosity > 5)
          printf("\t--> RP %i has %lu hits\n", *rp, hits.size());

        // insert to output collection
        (*output)[*rp] = RPTrackCandidate(hits, (uHits > 1 && vHits > 1));

        // RP hit?
        if (arm == 1 && uHits > 1 && vHits > 1)
        {
          unsigned int rpIdx = *rp % 10;
          if (rpIdx == 0 || rpIdx == 4) hitTopRPs++;
          if (rpIdx == 2 || rpIdx == 3) hitHorizontalRPs++;
          if (rpIdx == 1 || rpIdx == 5) hitBottomRPs++;
        }

        // Reconstructable RP hit?
        if (uHits >= minUVSensorsPerRP && vHits >= minUVSensorsPerRP) {
            activeRPs++;
        }
      }

      RPStationTrackFit trackFit;
      trackFit.ax = ax;
      trackFit.bx = bx;
      trackFit.ay = ay;
      trackFit.by = by;
      trackFit.z0 = z0;
      trackFit.valid = activeRPs >= minRPsPerStation;
      (*trackColl)[arm].push_back(trackFit); // TODO: should be stationID, not arm

      if (verbosity > 5)
        printf("hitTopRPs=%u, hitHorizontalRPs=%u, hitBottomRPs=%u\n", hitTopRPs, hitHorizontalRPs, hitBottomRPs);

      // dump proton
      if ( (dumpTopRPProtons && hitTopRPs > 0) || (dumpBottomRPProtons && hitBottomRPs > 0) || (dumpHorizontalRPProtons && hitHorizontalRPs > 0) )
      {
        printf("%+E\t%+E\t%+E\t%+E\t%+E\t%+E\n", vertex.x(), vertex.y(), vertex.z(), p.x(), p.y(), p.z());
      }
    }
  }

  // save output
  event.put(output);
  event.put(recoColl);
  event.put(trackColl);
}

DEFINE_FWK_MODULE(RPFastFullSimulation200);
