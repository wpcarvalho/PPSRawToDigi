/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*  Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandomEngine.h"

#include <vector>

#include "Math/GSLRndmEngines.h"

namespace RPStationMultiTrackFinderFitter {

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Track source reading events from files
 *
 * It is possible to specify multiple sources of particles. Each source can have
 * a probability of it being used.
 *
 * Number of particles to be returned per event must be specified.
 **/
class CachedMultiTrackSource : public edm::EDProducer
{
public:
    /// caches all the track sources in pd
    CachedMultiTrackSource(const edm::ParameterSet &);
    virtual ~CachedMultiTrackSource();
    virtual void beginJob();

    struct Track {
        /// particle direction of movement
        double vx, vy, vz;
        /// momentum of the particle
        double px, py, pz;
    };

    /// single track source
    struct ProcessData {
        double weight;
        std::vector<Track> tracks;
        unsigned long idx;
    };

protected:
    unsigned int verbosity;
    unsigned int tracksPerEvent;

    /// collection of all sources to be mixed
    std::vector<ProcessData> pd;

    ROOT::Math::GSLRandomEngine gslRndmEng;

    virtual void produce(edm::Event&, const edm::EventSetup&);
};



using namespace HepMC;
using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

CachedMultiTrackSource::CachedMultiTrackSource(const ParameterSet& pSet) :
  verbosity(pSet.getUntrackedParameter<unsigned int>("verbosity", 1)),
  tracksPerEvent(pSet.getParameter<unsigned int>("tracksPerEvent"))
{
  // load input
  if (verbosity > 0)
    printf(">> CachedMultiTrackSource::CachedMultiTrackSource > input:\n");

  const vector<ParameterSet> input = pSet.getParameter< vector<ParameterSet> >("input");
  double S_w = 0.;
  for (unsigned int pi = 0; pi < input.size(); pi++) {
      ProcessData p;
      p.weight = input[pi].getParameter<double>("weight");
      S_w += p.weight;
      p.idx = 0;
      string fileName = input[pi].getParameter<string>("fileName");

      FILE *f_in = fopen(fileName.c_str(), "r");
      if (!f_in)
      {
          throw cms::Exception("CachedMultiTrackSource::CachedMultiTrackSource")
              << "Can't open file `" << fileName << "'.";
      }

      while (!feof(f_in))
      {
          Track t;
          int re = fscanf(f_in, "%lE %lE %lE %lE %lE %lE", &t.vx, &t.vy, &t.vz, &t.px, &t.py, &t.pz);
          if (re == 6)
              p.tracks.push_back(t);
          else
              break;
      }

      fclose(f_in);

      if (verbosity > 0)
          printf("\tprocess %u: weight = %.3f, fileName = %s, tracks loaded = %lu\n", pi, p.weight,
                  fileName.c_str(), p.tracks.size());

      pd.push_back(p);
  }

  // normalise weights
  for (unsigned int pi = 0; pi < pd.size(); pi++)
    pd[pi].weight /= S_w;

  // set random seed
  Service<RandomNumberGenerator> rng;
  CLHEP::HepRandomEngine &rndEng = rng->getEngine();
  unsigned int seed = rndEng.getSeed();

  gslRndmEng.Initialize();
  gslRndmEng.SetSeed(seed);
  if (verbosity > 0)
    printf(">> CachedMultiTrackSource > seed = %u\n", seed);

  produces<HepMCProduct>();
}

//----------------------------------------------------------------------------------------------------

CachedMultiTrackSource::~CachedMultiTrackSource()
{
}

//----------------------------------------------------------------------------------------------------

void CachedMultiTrackSource::beginJob()
{
}

//----------------------------------------------------------------------------------------------------

void CachedMultiTrackSource::produce(edm::Event &e, const edm::EventSetup &es)
{
    if (verbosity > 1)
        printf("\n>> CachedMultiTrackSource::produce(event %i)\n", e.id().event());

    // create event structure
    GenEvent* gEv = new GenEvent();
    gEv->set_event_number(e.id().event());

    // generate number of tracks per process
    vector<double> v_prob;
    for (unsigned int pi = 0; pi < pd.size(); pi++)
        v_prob.push_back(pd[pi].weight);

    vector<unsigned int> v_nTrack = gslRndmEng.Multinomial(tracksPerEvent, v_prob);

    if (verbosity > 1)
        for (unsigned int pi = 0; pi < pd.size(); pi++)
            printf("\tprocess %u: %u tracks\n", pi, v_nTrack[pi]);

    // mix tracks
    for (unsigned int pi = 0; pi < pd.size(); pi++)
    {
        for (unsigned int ti = 0; ti < v_nTrack[pi]; ti++)
        {
            // get track
            if (pd[pi].idx >= pd[pi].tracks.size())
                throw cms::Exception("CachedMultiTrackSource::produce")
                    << "Ran out of tracks in sample of process " << pi << ".";

            const Track t = pd[pi].tracks[pd[pi].idx];
            pd[pi].idx++;

            // add vertex
            GenVertex* gVx = new GenVertex(FourVector(t.vx, t.vy, t.vz, 0.));
            gEv->add_vertex(gVx);

            // add particle
            const int PID = 2212;
            const int FinalState = 1;
            double m = 0.938; // proton mass in GeV
            double E = sqrt(t.px*t.px + t.py*t.py + t.pz*t.pz + m*m);
            GenParticle* gPe;
            gPe = new GenParticle(HepMC::FourVector(t.px, t.py, t.pz, E), PID, FinalState);
            gPe->suggest_barcode(1);
            gVx->add_particle_out(gPe);
        }
    }

    if (verbosity > 1)
        gEv->print();

    // store generator event to the FW event
    auto_ptr<HepMCProduct> output(new HepMCProduct());
    output->addHepMCData(gEv);
    e.put(output);
}

DEFINE_FWK_MODULE(CachedMultiTrackSource);

} // namespace RPStationMultiTrackFinderFitter
