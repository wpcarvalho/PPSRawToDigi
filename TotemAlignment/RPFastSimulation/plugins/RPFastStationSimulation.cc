/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
****************************************************************************/

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGauss.h"
#include "CLHEP/Random/RandExponential.h"

#include "DataFormats/TotemRPDataTypes/interface/RPRecoHit.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/TotemRecords/interface/MisalignedGeometryRecord.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPStationTrackFitCollection.h"

//#define DEBUG 1

/**
 *\brief Fast (no G4) proton simulation in within one station.
 * Uses misaligned geometry.
 */
class RPFastStationSimulation : public edm::EDProducer
{
  public:
	RPFastStationSimulation(const edm::ParameterSet &);
	virtual ~RPFastStationSimulation();

  protected:
    /// whether a HepMC description of the proton shall be saved in the event
    bool makeHepMC;
    
    /// whether the hits of the proton shall be calculated and saved
    bool makeHits;

	unsigned int verbosity;

    /// the list of RPs to simulate
    vector<unsigned int> RPs;

    /// number of particles to generate per event
    unsigned int particlesPerEvent;

    /// the "origin" of tracks, in mm
    double z0;                                  

    /// in mm
    double pitch;
    
    /// widt of double-strip cluster region, in mm
    double dscrWidth;
    
    /// position (value of v) where double-strip cluster region starts, in mm
    double dscrStart;
    
    /// position (value of v) where double-strip cluster region eds, in mm
    double dscrEnd;

    /// if true, double-strip clusters will be assigned the uncertainty of sigma0/2,
    /// if false, sigma0 will be used
    bool dscReduceUncertainty;

    /// size of insensitive margin at sensor's edge facing the beam, in mm
    double insensitiveMargin;

    /// variance of box distribution on (-P/2, +P/2)
    double sigma0;

    /// whether measurement values shall be rounded to the nearest strip
    bool roundToPitch;
  
    /// the index key for output RPStationTrackFitCollection
    unsigned int stationId;

    /// minimum number of U (or V) sensors required to consider the RP active
    unsigned int minUVSensorsPerRP;

    /// minimum number of active RPs to consider the track "reconstructable"
    unsigned int minRPsPerStation;

    struct Distribution
    {
      enum Type { dtBox, dtGauss, dtGaussLimit, dtExpt } type;
      double x_mean, x_width, x_min, x_max;
      double y_mean, y_width, y_min, y_max;
      Distribution(const edm::ParameterSet &);
      void Generate(double &x, double &y);
    };

    /// position parameters in mm
    Distribution position_dist;

    /// angular parameters in rad
    Distribution angular_dist;

    //---------- internal parameters ----------

    /// v position of strip 0, in mm
    double stripZeroPosition;

    /// particle energy and momentum
    double E, p;

    edm::ESHandle<TotemRPGeometry> geom;

    static CLHEP::HepRandomEngine *rndEng;

    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
	virtual void produce(edm::Event &, const edm::EventSetup&);

    void GenerateTrack(unsigned int, HepMC::GenEvent* gEv,
      auto_ptr< edm::DetSetVector<RPRecoHit> > &recoColl,
      auto_ptr< RPStationTrackFitCollection > &trackColl);
};

using namespace edm;
using namespace std;
using namespace CLHEP;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

CLHEP::HepRandomEngine *RPFastStationSimulation::rndEng = NULL;

RPFastStationSimulation::RPFastStationSimulation(const edm::ParameterSet &ps) :
  makeHepMC(ps.getParameter<bool>("makeHepMC")),
  makeHits(ps.getParameter<bool>("makeHits")),
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  RPs(ps.getParameter< vector<unsigned int> >("RPs")),
  particlesPerEvent(ps.getParameter<unsigned int>("particlesPerEvent")),
  z0(ps.getParameter<double>("z0")),
  pitch(ps.getParameter<double>("pitch")),
  dscrWidth(ps.getParameter<double>("dscrWidth")),
  dscrStart(pitch/2. - dscrWidth/2.),
  dscrEnd(pitch/2. + dscrWidth/2.),
  dscReduceUncertainty(ps.getParameter<bool>("dscReduceUncertainty")),
  insensitiveMargin(ps.getParameter<double>("insensitiveMargin")),
  sigma0(pitch/sqrt(12)),
  roundToPitch(ps.getParameter<bool>("roundToPitch")),
  stationId(ps.getParameter<unsigned int>("stationId")),
  minUVSensorsPerRP(ps.getParameter<unsigned int>("minUVSensorsPerRP")),
  minRPsPerStation(ps.getParameter<unsigned int>("minRPsPerStation")),
  position_dist(ps.getParameterSet("position_distribution")),

  angular_dist(ps.getParameterSet("angular_distribution"))

{
  // v position of strip 0
  stripZeroPosition = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_
    - RPTopology::y_width_/2.;

  // initialize random engine
  Service<RandomNumberGenerator> rng;
  rndEng = &(rng->getEngine());
  if (verbosity > 0)
      printf(">> RPFastStationSimulation > seed = %li\n", rndEng->getSeed());

  // register the output
  if (makeHepMC)
    produces<HepMCProduct>();

  if (makeHits)
  {
    produces< DetSetVector<RPRecoHit> >();
    produces< RPStationTrackFitCollection >();
  }
}

//----------------------------------------------------------------------------------------------------

RPFastStationSimulation::Distribution::Distribution(const edm::ParameterSet &ps)
{
  // get type
  string typeName = ps.getParameter<string>("type");
  if (!typeName.compare("box"))
    type = dtBox;
  else if (!typeName.compare("gauss"))
      type = dtGauss;
    else if (!typeName.compare("gauss-limit"))
        type = dtGaussLimit;
      else if (!typeName.compare("expt"))
          type = dtExpt;
        else
          throw cms::Exception("RPFastStationSimulation") << "Unknown distribution type `" << typeName << "'.";

  x_mean = ps.getParameter<double>("x_mean");
  x_width = ps.getParameter<double>("x_width");
  x_min = ps.getParameter<double>("x_min");
  x_max = ps.getParameter<double>("x_max");

  y_mean = ps.getParameter<double>("y_mean");
  y_width = ps.getParameter<double>("y_width");
  y_min = ps.getParameter<double>("y_min");
  y_max = ps.getParameter<double>("y_max");
}

//----------------------------------------------------------------------------------------------------

RPFastStationSimulation::~RPFastStationSimulation()
{
}

//----------------------------------------------------------------------------------------------------

void RPFastStationSimulation::beginRun(edm::Run const&, edm::EventSetup const& es)
{
  // get geometry
  es.get<MisalignedGeometryRecord>().get(geom);

  // get beam energy
  if (makeHepMC)
  {
    ESHandle<BeamOpticsParams> optPar;
    es.get<BeamOpticsParamsRcd>().get(optPar);

    E = optPar->GetBeamEnergy();
    p = optPar->GetBeamMomentum();

    if (verbosity > 5)
      printf(">> RPFastStationSimulation::beginRun > E = %E, p = %E\n", E, p);
  }
}

//----------------------------------------------------------------------------------------------------

void RPFastStationSimulation::Distribution::Generate(double &x, double &y)
{
  switch (type)
  {
    case dtBox:
      x = x_mean + x_width * (rndEng->flat() - 0.5);
      y = y_mean + y_width * (rndEng->flat() - 0.5);
      break;
    case dtGauss:
      x = x_mean + RandGauss::shoot(rndEng) * x_width;
      y = y_mean + RandGauss::shoot(rndEng) * y_width;
      break;
    case dtGaussLimit:
      {
        double u_x = rndEng->flat(), u_y = rndEng->flat();

        double cdf_x_min = (1. + TMath::Erf((x_min - x_mean) / x_width / sqrt(2.))) / 2.;
        double cdf_x_max = (1. + TMath::Erf((x_max - x_mean) / x_width / sqrt(2.))) / 2.;
        double a_x = cdf_x_max - cdf_x_min, b_x = cdf_x_min;
  
        double cdf_y_min = (1. + TMath::Erf((y_min - y_mean) / y_width / sqrt(2.))) / 2.;
        double cdf_y_max = (1. + TMath::Erf((y_max - y_mean) / y_width / sqrt(2.))) / 2.;
        double a_y = cdf_y_max - cdf_y_min, b_y = cdf_y_min;
  
        x = x_mean + x_width * sqrt(2.) * TMath::ErfInverse(2.*(a_x * u_x + b_x) - 1.); 
        y = y_mean + y_width * sqrt(2.) * TMath::ErfInverse(2.*(a_y * u_y + b_y) - 1.); 
      }

      break;
    case dtExpt:
      // TODO
      x = y = 0.;
      break;
    default:
      x = y = 0.;
  }
}

//----------------------------------------------------------------------------------------------------

void RPFastStationSimulation::GenerateTrack(unsigned int idx,
  HepMC::GenEvent* gEv,
  auto_ptr< edm::DetSetVector<RPRecoHit> > &recoColl,
  auto_ptr< RPStationTrackFitCollection > &trackColl)
{
  // generate track
  double bx = 0., by = 0., ax = 0., ay = 0.;
  position_dist.Generate(bx, by);
  angular_dist.Generate(ax, ay);

  if (verbosity > 5)
    printf("\tax = %.3f mrad, bx = %.3f mm, ay = %.3f mrad, by = %.3f mm, z0 = %.3f m\n", ax*1E3, bx, ay*1E3, by, z0*1E-3);

  // add HepMC track description
  if (makeHepMC)
  { 
    GenVertex* gVx = new GenVertex(HepMC::FourVector(bx, by, z0, 0.));
    gEv->add_vertex(gVx);

    GenParticle* gPe; 
    double az = sqrt(1. - ax*ax - ay*ay);
    gPe = new GenParticle(HepMC::FourVector(p*ax, p*ay, p*az, E), 2212, 1); // add a proton in final state
    gPe->suggest_barcode(idx + 1);
    gVx->add_particle_out(gPe);
  }

  if (makeHits)
  {
    unsigned int activeRPs = 0;

    // loop over all RPs
    for (vector<unsigned int>::iterator it = RPs.begin(); it != RPs.end(); ++it)
    {
      // u and v hit counter
      unsigned int uHits = 0, vHits = 0;
  
      // loop over sensors
      const set<unsigned int> &dets = geom->DetsInRP(*it);
      for (set<unsigned int>::iterator dit = dets.begin(); dit != dets.end(); ++dit)
      {
        // collection of hits
        unsigned int rawId = TotRPDetId::DecToRawId(*dit);
        DetSet<RPRecoHit> &hits = recoColl->find_or_insert(rawId);

        // calculate hit in global coordinates
        DetGeomDesc *ginfo = geom->GetDetector(rawId);
        double z = ginfo->translation().z();
        double x = bx + ax * (z - z0);
        double y = by + ay * (z - z0);
  
        // convert to local coordinates
        CLHEP::Hep3Vector hit(x, y, z);
        hit = geom->GlobalToLocal(rawId, hit);
        double u = hit.x();
        double v = hit.y();
  
        if (verbosity > 5)
          printf("\t%4u: u=%+8.4f, v=%+8.4f, x=%+8.4f, y=%+8.4f, z=%+11.3f", *dit, u, v, x, y, z);
  
        // is it within detector?
        if (!RPTopology::IsHit(u, v, insensitiveMargin))
        {
          if (verbosity > 5)
            printf(" | no hit\n");
          continue;
        }
  
        if (TotRPDetId::IsStripsCoordinateUDirection(*dit))
          uHits++;
        else
          vHits++;
  
        // measurement uncertainty
        double sigma = 1E-9;  // value for no rounding, 0 may let some algorithms down
  
        // round the measurement
        if (roundToPitch)
        {
          double m = stripZeroPosition - v;
          signed int strip = (int) floor(m / pitch);
          double offset = m - strip*pitch;
          double dstrip = 0.;
  
          if (offset < dscrStart)
          {
            sigma = sigma0;
            dstrip = 0.;
          } else {
            if (offset < dscrEnd)
            {
              sigma = (dscReduceUncertainty) ? sigma0/2. : sigma0;
              dstrip = 0.5;
            } else {
              sigma = sigma0;
              dstrip = 1.;
            }
          }
  
          v = stripZeroPosition - pitch * (strip + dstrip);
          if (verbosity > 5)
            printf(" | s0=%+8.4f, strip=%+6.1f", stripZeroPosition, strip+dstrip);
        }
  
        if (verbosity > 5)
          printf(" | m=%+8.4f, sigma=%+8.4f\n", v, sigma);
  
        hits.push_back(RPRecoHit(rawId, v, sigma));
      }
      
      if (verbosity > 5)
        printf("\tRP %i has %u U hits and %u V hits\n", *it, uHits, vHits);

      if (uHits >= minUVSensorsPerRP && vHits >= minUVSensorsPerRP)
        activeRPs++;
    }

    RPStationTrackFit tf;
    tf.ax = ax;
    tf.ay = ay;
    tf.bx = bx;
    tf.by = by;
    tf.z0 = z0;
    tf.valid = (activeRPs >= minRPsPerStation);
    (*trackColl)[stationId].push_back(tf);
  } 
}

//----------------------------------------------------------------------------------------------------

void RPFastStationSimulation::produce(edm::Event &event, const edm::EventSetup &es)
{
  if (verbosity > 2)
    printf(">> RPFastStationSimulation::produce > event %i\n", event.id().event());

  if (verbosity > 5)
    printf("\tseed = %li\n", rndEng->getSeed());

  // initialize products
  GenEvent* gEv = new GenEvent();
  gEv->set_event_number(event.id().event());

  auto_ptr< DetSetVector<RPRecoHit> > recoColl(new DetSetVector<RPRecoHit>());
  auto_ptr< RPStationTrackFitCollection > trackColl(new RPStationTrackFitCollection());

  // run particle loop
  for (unsigned int pi = 0; pi < particlesPerEvent; pi++)
  {
    if (verbosity > 5)
      printf("    generating track %u\n", pi);

    GenerateTrack(pi, gEv, recoColl, trackColl);
  }

  // save products
  if (makeHepMC)
  { 
    auto_ptr<HepMCProduct> hepMCoutput(new HepMCProduct());
    hepMCoutput->addHepMCData(gEv);
    event.put(hepMCoutput);
  }

  if (makeHits)
  {
    event.put(recoColl);
    event.put(trackColl);
  }
}

DEFINE_FWK_MODULE(RPFastStationSimulation);
