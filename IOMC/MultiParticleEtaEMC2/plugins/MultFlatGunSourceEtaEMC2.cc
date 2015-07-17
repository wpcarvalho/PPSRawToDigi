
#include <ostream>

#include "IOMC/MultiParticleEtaEMC2/interface/MultFlatGunSourceEtaEMC2.h"


#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"

// #include "FWCore/Framework/intercface/ESHandle.h"
// #include "FWCore/Framework/interface/EventSetup.h"
#include "SimGeneral/HepPDTRecord/interface/ParticleDataTable.h"

#include "FWCore/Utilities/interface/Exception.h"

#include <iostream>

using namespace edm;
using namespace std;
using namespace CLHEP;

MultFlatGunSourceEtaEMC2::MultFlatGunSourceEtaEMC2( const ParameterSet& pset) :
   fEvt(0)
   // fPDGTable( new DefaultConfig::ParticleDataTable("PDG Table") )
{

   ParameterSet defpset ;
   //ParameterSet pgun_params = pset.getParameter<ParameterSet>("PGunParameters") ;
   ParameterSet pgun_params = 
      pset.getUntrackedParameter<ParameterSet>("PGunParameters", defpset ) ;
  
   // although there's the method ParameterSet::empty(),  
   // it looks like it's NOT even necessary to check if it is,
   // before trying to extract parameters - if it is empty,
   // the default values seem to be taken
   vector<int> defids ;
   defids.push_back(13) ;

   vector<double> defetasmin;
   vector<double> defetasmax;  
   vector<double> defphismin;  
   vector<double> defphismax;
   defetasmin.push_back(-5.5);
   defetasmin.push_back(5.5);
   defphismin.push_back(-3.14);
   defphismin.push_back(3.14);


   fPartIDs    = pgun_params.getUntrackedParameter< vector<int> >("PartID",defids);  
   // fMinEtas     = pgun_params.getUntrackedParameter<vector<double> >("MinEtas",defetasmin);
   // fMaxEtas     = pgun_params.getUntrackedParameter<vector<double> >("MaxEtas",defetasmax);
   fMinPhis     = pgun_params.getUntrackedParameter<vector<double> >("MinPhis",defphismin);
   fMaxPhis     = pgun_params.getUntrackedParameter<vector<double> >("MaxPhis",defphismax);



  fVerbosity = pset.getUntrackedParameter<int>( "Verbosity",0 ) ;

   Service<RandomNumberGenerator> rng;
   long seed = (long)(rng->mySeed()) ;
   fRandomEngine = new CLHEP::HepJamesRandom(seed) ;
   fRandomGenerator = new CLHEP::RandFlat(fRandomEngine) ;
   
   fAddAntiParticle = pset.getUntrackedParameter("AddAntiParticle", false) ;

   InputMCFileName = pgun_params.getUntrackedParameter<std::string>("InputMCFileName");
   InputMCHistoName= pgun_params.getUntrackedParameter<std::string>("InputMCHistoName");

    MCfile = boost::shared_ptr<TFile> (new TFile(InputMCFileName.c_str())); 
   MCHisto=(TH2F*)MCfile->Get(InputMCHistoName.c_str());

   
}

MultFlatGunSourceEtaEMC2::~MultFlatGunSourceEtaEMC2()
{
  
  if ( fRandomGenerator != NULL ) delete fRandomGenerator;
  // do I need to delete the Engine, too ?
  
  // no need to cleanup GenEvent memory - done in HepMCProduct
  // if (fEvt != NULL) delete fEvt ; // double check
  // delete fPDGTable;
  
}


void MultFlatGunSourceEtaEMC2::beginRun(Run& r, EventSetup const& es)
{

   es.getData( fPDGTable ) ;
   //MCfile = boost::shared_ptr<TFile> (new TFile(InputMCFileName.c_str())); 
   //MCfile = TFile::Open(InputMCFileName.c_str()); 
   //MCHisto=(TH2F*)MCfile->Get(InputMCHistoName.c_str());
   return ;

}
