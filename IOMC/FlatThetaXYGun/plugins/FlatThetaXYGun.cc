#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "CLHEP/Random/RandFlat.h"

#include <auto_ptr.h>

//----------------------------------------------------------------------------------------------------

class FlatThetaXYGun : public edm::EDProducer
{
  public:
    FlatThetaXYGun(const edm::ParameterSet &);
    virtual ~FlatThetaXYGun();
  
  private:
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void produce(edm::Event&, const edm::EventSetup&);
  
  protected:
    unsigned int verbosity;

    // range for theta_x and theta_y, in rad
    double th_x_min, th_x_max;
    double th_y_min, th_y_max;

    // energy, mass and momentum of the particles, in GeV
    double energy, mass, momentum;

    // particle code (PDG)
    int particle_id;

    // flags to select arms where the particle should be generated
    bool generate_right_arm, generate_left_arm;

    // random generator
    CLHEP::RandFlat *randomGenerator;
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;
using namespace HepMC;
using namespace CLHEP;
using namespace HepPDT;

//----------------------------------------------------------------------------------------------------

FlatThetaXYGun::FlatThetaXYGun(const edm::ParameterSet &ps) : randomGenerator(NULL)
{
  verbosity = ps.getUntrackedParameter<unsigned int>("verbosity");
  
  th_x_min = ps.getParameter<double>("th_x_min");
  th_x_max = ps.getParameter<double>("th_x_max");
  th_y_min = ps.getParameter<double>("th_y_min");
  th_y_max = ps.getParameter<double>("th_y_max");

  energy = ps.getParameter<double>("energy");

  particle_id = ps.getParameter<int>("particle_id");

  generate_right_arm = ps.getParameter<bool>("generate_right_arm");
  generate_left_arm = ps.getParameter<bool>("generate_left_arm");
  
  // initialise random-number generator
  Service<RandomNumberGenerator> rng;
  if (!rng.isAvailable())
    throw cms::Exception("FlatThetaXYGun::FlatThetaXYGun") << "Random number service not available.";

  HepRandomEngine &rndEng = rng->getEngine();
  randomGenerator = new RandFlat(rndEng);
  unsigned int seed = rndEng.getSeed();

  if (verbosity)
  {
    printf(">> FlatThetaXYGun::FlatThetaXYGun>\n");
    printf("\ttheta_x: min = %.1E, max = %.1E\n", th_x_min, th_x_max);
    printf("\ttheta_y: min = %.1E, max = %.1E\n", th_y_min, th_y_max);
    printf("\tenergy: %.1E GeV\n", energy);
    printf("\tparticle ID: %u\n", particle_id);
    printf("\tgenerate in arm: left=%u, right=%u\n", generate_left_arm, generate_right_arm);
    printf("\tseed = %u\n", seed);
  }

  // declare output
  produces<HepMCProduct>();
}

//----------------------------------------------------------------------------------------------------

FlatThetaXYGun::~FlatThetaXYGun()
{
  if (randomGenerator)
    delete randomGenerator;
}

//----------------------------------------------------------------------------------------------------

void FlatThetaXYGun::beginRun(Run const &, EventSetup const &es)
{
  // calculate mass and momentum
  ESHandle<ParticleDataTable> fPDGTable;
  es.getData(fPDGTable);
  const ParticleData* pData = fPDGTable->particle(ParticleID(particle_id));
  mass = pData->mass().value();
  momentum = sqrt(energy*energy - mass*mass);
}

//----------------------------------------------------------------------------------------------------

void FlatThetaXYGun::produce(edm::Event &ev, const edm::EventSetup &es)
{
  GenEvent *mc_ev = new GenEvent() ;

  GenVertex *vtx = new GenVertex(HepMC::FourVector(0., 0., 0.));
  mc_ev->add_vertex(vtx);

  const int finalState = 1;

  if (generate_right_arm)
  {
    double th_x = randomGenerator->fire(th_x_min, th_x_max);
    double th_y = randomGenerator->fire(th_y_min, th_y_max);
    double th = sqrt(th_x*th_x + th_y*th_y);
    double p_x = momentum * th_x, p_y = momentum * th_y, p_z = momentum * cos(th);
    GenParticle* particle = new GenParticle(FourVector(+p_x, +p_y, +p_z, energy), particle_id, finalState);
    particle->suggest_barcode(1);
    vtx->add_particle_out(particle);

    if (verbosity > 1)
      printf(">> FlatThetaXYGun::produce > right arm, th_x = %.1E, th_y = %.1E\n", th_x, th_y);
  }

  if (generate_left_arm)
  {
    double th_x = randomGenerator->fire(th_x_min, th_x_max);
    double th_y = randomGenerator->fire(th_y_min, th_y_max);
    double th = sqrt(th_x*th_x + th_y*th_y);
    double p_x = momentum * th_x, p_y = momentum * th_y, p_z = momentum * cos(th);
    GenParticle* particle = new GenParticle(FourVector(-p_x, -p_y, -p_z, energy), particle_id, finalState);
    particle->suggest_barcode(2);
    vtx->add_particle_out(particle);

    if (verbosity > 1)
      printf(">> FlatThetaXYGun::produce > left arm, th_x = %.1E, th_y = %.1E\n", th_x, th_y);
  }
  
  if (verbosity > 2)
    mc_ev->print();

  // save generated event
  auto_ptr<HepMCProduct> output(new HepMCProduct());
  output->addHepMCData(mc_ev);
  ev.put(output);
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(FlatThetaXYGun);
