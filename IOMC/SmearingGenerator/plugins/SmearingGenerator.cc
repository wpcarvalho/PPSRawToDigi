/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: SmearingGenerator.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "TRandom2.h"
#include "TLorentzVector.h"
#include "TVector3.h"

#include "IOMC/SmearingGenerator/interface/SmearingGenerator.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemCondFormats/DataRecord/interface/BeamOpticsParamsRcd.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "CLHEP/Random/RandomEngine.h"
#include "CLHEP/Random/RandGauss.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

#include <iostream>

using namespace std;
using namespace HepMC;
using namespace CLHEP;
using namespace edm;

//----------------------------------------------------------------------------------------------------

SmearingGenerator::SmearingGenerator(const edm::ParameterSet& pSet) :
  verbosity(pSet.getUntrackedParameter<unsigned int>("verbosity", 1)),
  modifyLabel(pSet.getParameter<string>("modifyLabel")),
  originalLabel(pSet.getParameter<string>("originalLabel"))
{
  if (verbosity > 0)
    cout << "SmearingGenerator.verbosity="<<verbosity<<endl;

  // register data to consume
  tokenHepMcProduct = consumes<edm::HepMCProduct>(modifyLabel);

  // register fake output
  produces<HepMCProduct>(originalLabel);
}

//----------------------------------------------------------------------------------------------------

SmearingGenerator::~SmearingGenerator()
{
}

//----------------------------------------------------------------------------------------------------

void SmearingGenerator::produce(edm::Event& event, const edm::EventSetup& es)
{
  // initialize random engine
  Service<RandomNumberGenerator> rng;
  HepRandomEngine &rndEng = rng->getEngine(event.streamID());
  rand = new RandGauss(rndEng);

  // retrieve (the only) HepMCEvent from the event
  edm::Handle<edm::HepMCProduct> mcProd;
  event.getByToken(tokenHepMcProduct, mcProd);

  // back up the orignal event
  auto_ptr<HepMCProduct> output(new HepMCProduct(*mcProd));
  event.put(output, originalLabel);

  // apply smearing to the actual event
  // the type cast below is a really aweful practique but it is necessary due to the CMSSW framework
  GenEvent *evt = (GenEvent *) mcProd->GetEvent();
  ApplyBeamSmearing(evt);
  ApplyVertexSmearing(evt);

  return;
}

//----------------------------------------------------------------------------------------------------



//----------------------------------------------------------------------------------------------------

void SmearingGenerator::ApplyVertexSmearing(HepMC::GenEvent *evt)
{
  // generate vertex smearing
  double v_x = rand->fire(MeanX, SiX);
  double v_y = rand->fire(MeanY, SiY);
  double v_z = rand->fire(MeanZ, SiZ);
  double v_t = 0.;

  // shift all vertices
  if (verbosity > 5) cout << ">> SmearingGenerator::ApplyVertexSmearing" << endl;

  for (GenEvent::vertex_const_iterator vt = evt->vertices_begin(); vt != evt->vertices_end(); ++vt) {
	FourVector ov = (*vt)->position();
	if (verbosity > 5) cout << "Vertex\n\told: " << ov.x() << " " << ov.y() << " " << ov.z() << " " << ov.t() << endl;
	(*vt)->set_position(FourVector(ov.x() + v_x, ov.y() + v_y, ov.z() + v_z, ov.t() + v_t));
	ov = (*vt)->position();
	if (verbosity > 5) cout << "\tnew: " << ov.x() << " " << ov.y() << " " << ov.z() << " " << ov.t() << endl;
  }
}

//----------------------------------------------------------------------------------------------------

void SmearingGenerator::ApplyBeamSmearing(HepMC::GenEvent *evt)
{
  // generate energy/angular smearing
  double al1 = Al;
  double al2 = -Al;
  double thx1 = rand->fire(0., SiThX);
  double thy1 = rand->fire(0., SiThY);
  double thx2 = rand->fire(0., SiThX);
  double thy2 = rand->fire(0., SiThY);
  double xi1 = rand->fire(MeanXi, SiXi);
  double xi2 = rand->fire(MeanXi, SiXi);

  // compute transform parameters
  double m = 0.938271;
  double p_nom = sqrt(E_CMS * E_CMS - m*m);

  double thz1 = sqrt(1. - thx1*thx1 - thy1*thy1);
  double thz2 = sqrt(1. - thx2*thx2 - thy2*thy2);

  TVector3 p1(cos(al1)*thx1 + sin(al1)*thz1, thy1, -sin(al1)*thx1 + cos(al1)*thz1);
  p1 *= (1. + xi1) * p_nom;
  double E1 = sqrt(p1.Mag()*p1.Mag() + m*m);
  TLorentzVector P1(p1.x(), p1.y(), p1.z(), E1);

  TVector3 p2(cos(al2)*thx2 + sin(al2)*thz2, thy2, -sin(al2)*thx2 + cos(al2)*thz2);
  p2 *= -(1. + xi2) * p_nom;
  double E2 = sqrt(p2.Mag()*p2.Mag() + m*m);
  TLorentzVector P2(p2.X(), p2.Y(), p2.Z(), E2);

  double factor = (P1 + P2).Mag() / 2. / E_CMS;             // energy correction factor

  TVector3 boost = (p1 + p2) * (1. / (E1 + E2));            // boost vector (direction * beta)
  double beta = (p1 + p2).Mag() / (E1 + E2);                // beta of the boost
  P1.Boost(-boost);
  TVector3 axis(P1.Y(), -P1.X(), 0.);                       // rotation axis
  double angle = -acos( P1.Z() / P1.Vect().Mag() );         // angle of rotation

  if (verbosity > 5) {
    cout << ">> SmearingGenerator::ApplyMomentumSmearing" << endl;
    //cout << "transformation:" << endl;
    //cout << "\tdir = " << boost.Unit() << "\n\tbeta = " << beta << "\n\taxis = " << axis << "\n\tangle " << angle << endl;
  }

  // apply transform to all particles
  for (GenEvent::particle_const_iterator particle = evt->particles_begin(); particle != evt->particles_end(); ++particle) {
    FourVector P = (*particle)->momentum();

    if (verbosity > 5)
      printf("\tbefore: (%+.0f| %+4.2f; %+4.2f; %+.0f)\n", P.t(), P.x(), P.y(), P.z());

    // energy scaling
    TVector3 p(P.x(), P.y(), P.z());
    double E = P.e() * factor, m = P.m();
    p = sqrt(E*E - m*m) / p.Mag() * p;
    TLorentzVector PP(p, E);

    // rotation
    if (fabs(angle) > 1E-8)
      PP.Rotate(angle, axis);

    // boost
    if (fabs(beta) > 1E-8)
      PP.Boost(boost);

    if (verbosity > 5)
      printf("\t after: (%+.0f| %+4.2f; %+4.2f; %+.0f), th_x = %.2E, th_y = %.2E\n\n",
          PP.T(), PP.X(), PP.Y(), PP.Z(), PP.X()/PP.Z(), PP.Y()/PP.Z());

    (*particle)->set_momentum(FourVector(PP.X(), PP.Y(), PP.Z(), PP.T()));
  }
}


DEFINE_FWK_MODULE(SmearingGenerator);

//virtual void beginRun(edm::Run const&iR, edm::EventSetup const&iE)
void SmearingGenerator::beginRun(edm::Run const&iR, edm::EventSetup const&es) {
	if(verbosity > 5){
		cout<<"SmearingGenerator beginRun"<<endl;
	}
  edm::ESHandle<BeamOpticsParams> par;
  es.get<BeamOpticsParamsRcd>().get(par);

  if(!par.isValid()){
      cout<<"<<SmearingGenerator:: Invalid BeamOpticsParams"<<endl;
	  throw cms::Exception("SmearingGenerator") << " edm::ESHandle<BeamOpticsParams> is invalid";
  }
  // angles in rad and positions in m
  Al = par->GetCrossingAngleX();
  SiThX = par->GetBeamDivergenceX();
  SiThY = par->GetBeamDivergenceY();
  MeanXi = par->GetMeanXi();
  SiXi = par->GetSigmaXi();
  E_CMS = par->GetBeamEnergy();
  MeanX = par->GetBeamDisplacementX() * 1000.;
  SiX = par->GetPrimVertSizeX() * 1000.;
  MeanY = par->GetBeamDisplacementY() * 1000.;
  SiY = par->GetPrimVertSizeY() * 1000.;
  MeanZ = par->GetBeamDisplacementZ() * 1000.;
  SiZ = par->GetPrimVertSizeZ() * 1000.;

  if (verbosity) {
    printf(">> SmearingGenerator::beginJob\n");
    printf("Al = %E\n", Al);
    printf("SiThX = %E\n", SiThX);
    printf("SiThY = %E\n", SiThY);
    printf("MeanXi = %E\n", MeanXi);
    printf("SiXi = %E\n", SiXi);
    printf("E_CMS = %E GeV\n", E_CMS);
    printf("MeanX = %E mm\n", MeanX);
    printf("SiX = %E mm\n", SiX);
    printf("MeanY = %E mm\n", MeanY);
    printf("SiY = %E mm\n", SiY);
    printf("MeanZ = %E mm\n", MeanZ);
    printf("SiZ = %E mm\n", SiZ);
  }
}
