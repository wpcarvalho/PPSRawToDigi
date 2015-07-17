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

//----------------------------------------------------------------------------------------------------

/**
 *\brief Fast (no G4) proton simulation from IP to all RP stations.
 */
class RPFastFullSimulation : public edm::EDProducer
{
  public:
    RPFastFullSimulation(const edm::ParameterSet &);
    virtual ~RPFastFullSimulation();

    struct BeamMisalignment
    {
      double z0;                ///< in mm
      double shift_x, shift_y;  ///< shift of the beam at z0 from IP, in um
      double tilt_x, tilt_y;    ///< angular shift of the beam, in mrad

      ///< caculates beam shift at position z (in mm), dx and dy in mm
      void CalculateShift(double z, double &dx, double &dy);
    };

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

    /// size of insensitive margin at sensor's edge facing the beam, in mm
    double insensitiveMargin;
    
    /// whether measurement values shall be rounded to the nearest strip
    bool roundToPitch;
    
    /// variance of box distribution on (-P/2, +P/2)
    double sigma0;

    /// angular limit to distinguish between forward, central and backward particles
    double thetaLim;

    /// v position of strip 0, in mm
    double stripZeroPosition;

    /// minimum number of U/V sensors active to tag the RP data as fittable
    unsigned int minUVSensorsPerRP;

    std::map<unsigned int, BeamMisalignment> bms;   ///< map: station id --> beam error

    virtual void beginRun(edm::Run&, edm::EventSetup const&);
    virtual void produce(edm::Event &, const edm::EventSetup&);
};


using namespace edm;
using namespace std;
using namespace HepMC;

//----------------------------------------------------------------------------------------------------

RPFastFullSimulation::RPFastFullSimulation(const ParameterSet &ps) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
    hepMCLabel(ps.getParameter<string>("hepMCLabel")),
    pitch(ps.getParameter<double>("pitch")),
    dscrWidth(ps.getParameter<double>("dscrWidth")),
    dscrStart(pitch/2. - dscrWidth/2.),
    dscrEnd(pitch/2. + dscrWidth/2.),
    insensitiveMargin(ps.getParameter<double>("insensitiveMargin")),
    roundToPitch(ps.getParameter<bool>("roundToPitch")),
    sigma0(pitch/sqrt(12)),
    thetaLim(ps.getParameter<double>("thetaLimit")),
    minUVSensorsPerRP(ps.getParameter<unsigned int>("minUVSensorsPerRP"))
{
  // v position of strip 0
  stripZeroPosition = RPTopology::last_strip_to_border_dist_ + (RPTopology::no_of_strips_-1)*RPTopology::pitch_
    - RPTopology::y_width_/2.;

  const vector<ParameterSet> &be = ps.getParameterSetVector("beam_misalignment");
  for (vector<ParameterSet>::const_iterator it = be.begin(); it != be.end(); ++it)
  {
    BeamMisalignment bm;
    bm.z0 = it->getParameter<double>("z0");
    bm.shift_x = it->getParameter<double>("shift_x");
    bm.shift_y = it->getParameter<double>("shift_y");
    bm.tilt_x = it->getParameter<double>("tilt_x");
    bm.tilt_y = it->getParameter<double>("tilt_y");

    bms[it->getParameter<unsigned int>("station")] = bm;
  }

  if (verbosity)
    printf(">> RPFastFullSimulation > %lu beam misalignments loaded.\n", bms.size());
  
  produces<RPTrackCandidateCollection>();
}

//----------------------------------------------------------------------------------------------------

RPFastFullSimulation::~RPFastFullSimulation()
{
}

//----------------------------------------------------------------------------------------------------

void RPFastFullSimulation::beginRun(edm::Run&, edm::EventSetup const& es)
{
  ESHandle<TotemRPGeometry> geom;
  es.get<MisalignedGeometryRecord>().get(geom);

  if (verbosity > 1)
  {
    printf(">> RPFastFullSimulation::beginJob > beam positions at RPs\n");

    for (TotemRPGeometry::RPDeviceMapType::const_iterator it = geom->beginRP(); it != geom->endRP(); ++it)
    {
      unsigned int rp_id = it->first;
      unsigned int st_id = rp_id / 10;
      double z = it->second->translation().z();

      double x=0, y=0;
      bms[st_id].CalculateShift(z, x, y);

      printf("\tRP=%3u, z=%+.0f mm | beam_x=%+.3f mm, beam_y=%+.3f mm\n", rp_id, z, x, y);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void RPFastFullSimulation::BeamMisalignment::CalculateShift(double z, double &dx, double &dy)
{
  dx = (shift_x + (z - z0)*tilt_x)*1E-3; 
  dy = (shift_y + (z - z0)*tilt_y)*1E-3; 

#ifdef DEBUG
  printf("CalculateShift: z = %E, z0 = %E, dx = %E, dy = %E\n", z, z0, dx, dy);
#endif
}

//----------------------------------------------------------------------------------------------------

void RPFastFullSimulation::produce(edm::Event &event, const EventSetup &es)
{
  if (verbosity > 5)
    printf(">> RPFastFullSimulation::produce > event %i\n", event.id().event());

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

  // loop over vertices
  for (GenEvent::vertex_const_iterator vit = mcEvent->vertices_begin(); vit != mcEvent->vertices_end(); ++vit)
  {
    const FourVector &vertex = (*vit)->position(); // in mm
    
    // loop over outgoing particles
    for (GenVertex::particles_out_const_iterator pit = (*vit)->particles_out_const_begin(); pit != (*vit)->particles_out_const_end(); ++pit)
    {
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
    
      // loop over all stations of the arm
      const set<unsigned int> stations = geom->StationsInArm(arm);
      for (set<unsigned int>::const_iterator st = stations.begin(); st != stations.end(); ++st)
      {
        //printf("station %u\n", *st);
        // beam misalignment for the station
        map<unsigned int, BeamMisalignment>::iterator bm = bms.find(*st);
  
        const set<unsigned int> RPs = geom->RPsInStation(*st);
        for (set<unsigned int>::const_iterator rp = RPs.begin(); rp != RPs.end(); ++rp)
        {
            //printf("RP %u\n", *rp);
            // center of the RP box, in mm
            double z0 = geom->GetRPDevice(*rp)->translation().z();
        
            // optical functions for the RP
            LHCOpticsApproximator *of = optFun->GetFunction(*rp);
            if (!of)
              throw cms::Exception("RPFastFullSimulation::produce") << "No optical function found for RP " << *rp;
        
            // get position and direction of the proton at the center of the RP box
            double xi = p.e() / optPar->GetBeamEnergy() - 1.;
            double parIn[5] = {vertex.x() * 1E-3, p.x()/p.rho(), vertex.y() * 1E-3, p.y()/p.rho(), xi};
            double parOut[5]; // in m, rad and xi dimensionless
            bool transportable = of->Transport(parIn, parOut, true);
        
            if (verbosity > 5)
              printf("RP%u\tz0 = %E\t p.z = %E\t transportable = %i, beam type = %i\n", *rp, z0, p.z(), transportable, of->GetBeamType());
  
            if (verbosity > 8)
            {
              printf("param input:  %E\t%E\t%E\t%E\t%E\n", parIn[0], parIn[1], parIn[2], parIn[3], parIn[4]);
              printf("param output: %E\t%E\t%E\t%E\t%E\n", parOut[0], parOut[1], parOut[2], parOut[3], parOut[4]);
            }
        
            // do not continue if the proton is stopped
            if (!transportable)
              continue;
        
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
              double z = geom->GetDetector(rawId)->translation().z(); // in mm
              double x = parOut[0]*1E3  + parOut[1] * (z - z0); // in mm
              double y = parOut[2]*1E3  + parOut[3] * (z - z0); // in mm
  
              // apply beam misalignment
              if (bm != bms.end())
              {
                double dx = 0., dy = 0.; // in mm
                bm->second.CalculateShift(z, dx, dy);
                x += dx;
                y += dy;
              }
        
              if (verbosity > 5)
                printf("\t%u\t%u\t%E\t%e\t%E\n", *dit, rawId, z, x, y);
        
              // convert to local coordinates
              CLHEP::Hep3Vector hit(x, y, z);
              hit = geom->GlobalToLocal(rawId, hit);
              double u = hit.x(); // in mm
              double v = hit.y(); // in mm
        
              // is it within detector?
              if (!RPTopology::IsHit(u, v, insensitiveMargin))
                continue;
  
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
                printf("\t\tv = %E\n", v);
        
              hits.push_back(RPRecoHit(rawId, v, sigma0));
            }
        
            if (verbosity > 5)
              printf("\tRP %i has %lu hits\n", *rp, hits.size());
        
            // insert to output collection
            bool fittable = (uHits > minUVSensorsPerRP && vHits > minUVSensorsPerRP);
            (*output)[*rp] = RPTrackCandidate(hits, fittable);
        }
      }
    }
  }

  // save output
  event.put(output);
}

DEFINE_FWK_MODULE(RPFastFullSimulation);
