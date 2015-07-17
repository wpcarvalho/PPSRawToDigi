/*
 *  $Date: 2009/11/07 20:17:24 $
 *  $Revision: 1.1.2.1 $ from FlatRandomEGunSource
 *  Modified By Mirko Berretti (TOTEM exp)
 */

#include <ostream>

#include "IOMC/MultiParticleEtaEMC2/interface/MultiParticleGunSourceEtaEMC2.h"
//#include "SimDataFormats/HepMCProduct/interface/HepMCProduct.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

using namespace edm;
using namespace std;

MultiParticleGunSourceEtaEMC2::MultiParticleGunSourceEtaEMC2(const ParameterSet& pset) : MultFlatGunSourceEtaEMC2(pset)
{

   ParameterSet defpset ;
   // ParameterSet pgun_params = pset.getParameter<ParameterSet>("PGunParameters") ;
   ParameterSet pgun_params = 
      pset.getUntrackedParameter<ParameterSet>("PGunParameters",defpset) ;
  
   vector<double> defEsmin;
   defEsmin.push_back(1.0);
   vector<double> defEsmax ;
   defEsmax.push_back(100.0);
 
  // doesn't seem necessary to check if pset is empty - if this
   // is the case, default values will be taken for params
   // fMinEs = pgun_params.getUntrackedParameter<vector<double> >("MinEs",defEsmin);
   
   //fMaxEs = pgun_params.getUntrackedParameter<vector<double> >("MaxEs",defEsmax);

   produces<HepMCProduct>();

   cout << "Internal FlatRandomEGun is initialzed" << endl ;
   
}

MultiParticleGunSourceEtaEMC2::~MultiParticleGunSourceEtaEMC2()
{
   // no need to cleanup fEvt since it's done in HepMCProduct
}

void MultiParticleGunSourceEtaEMC2::produce(edm::Event& e, const edm::EventSetup& es)
{

   if ( fVerbosity > 0 )
   {
      cout << " FlatRandomEGunSource : Begin New Event Generation" << endl ; 
   }
   
   // event loop (well, another step in it...)
          
   // no need to clean up GenEvent memory - done in HepMCProduct

   // here re-create fEvt (memory)
   //
   fEvt = new HepMC::GenEvent() ;
   
   // now actualy, cook up the event from PDGTable and gun parameters
   //

   // 1st, primary vertex
   //
   HepMC::GenVertex* Vtx = new HepMC::GenVertex( HepMC::FourVector(0.,0.,0.));
   
   // loop over particles
   //
   int barcode = 1;
   std::cout<<"Size of particles"<<fPartIDs.size()<<std::endl;
   for (unsigned int ip=0; ip<fPartIDs.size(); ip++)
   {
     //double energy = fRandomGenerator->fire(fMinEs.at(ip), fMaxEs.at(ip)) ;
     // double eta    = fRandomGenerator->fire(fMinEtas.at(ip), fMaxEtas.at(ip)) ;
       //double phi    = fRandomGenerator->fire(fMinPhis.at(ip), fMaxPhis.at(ip)) ;
       double phi    = fRandomGenerator->fire(fMinPhis.at(0), fMaxPhis.at(0)) ;
       Double_t energyr=10.; //default value
       Double_t etar=5.7; //default value
       MCHisto->GetRandom2(etar,energyr);
       double energy= (double) energyr;
       double eta= (double) etar;
       std::cout<<"Eta-Energy:  "<<eta<<"-"<<energy<<std::endl;
       int PartID = fPartIDs[ip] ;
       const HepPDT::ParticleData* 
          PData = fPDGTable->particle(HepPDT::ParticleID(abs(PartID))) ;
       double mass   = PData->mass().value() ;
       double mom2   = energy*energy - mass*mass ;
       double mom    = 0. ;
       if (mom2 > 0.) 
       {
          mom = sqrt(mom2) ;
       }
       else
       {
          mom = 0. ;
       }
       double theta  = 2.*atan(exp(-eta)) ;
       double px     = mom*sin(theta)*cos(phi) ;
       double py     = mom*sin(theta)*sin(phi) ;
       double pz     = mom*cos(theta) ;

       HepMC::FourVector p(px,py,pz,energy) ;
       HepMC::GenParticle* Part = 
           new HepMC::GenParticle(p,PartID,1);
       Part->suggest_barcode( barcode ) ;
       barcode++ ;
       Vtx->add_particle_out(Part);
       
       if ( fAddAntiParticle )
       {
          HepMC::FourVector ap(-px,-py,-pz,energy) ;
	  int APartID = -PartID ;
	  if ( PartID == 22 || PartID == 23 )
	  {
	     APartID = PartID ;
	  }
	  HepMC::GenParticle* APart =
	     new HepMC::GenParticle(ap,APartID,1);
	  APart->suggest_barcode( barcode ) ;
	  barcode++ ;
	  Vtx->add_particle_out(APart) ;
       }
       
   }
   fEvt->add_vertex(Vtx) ;
   fEvt->set_event_number(e.id().event()) ;
   fEvt->set_signal_process_id(20) ;  
   
   
   if ( fVerbosity > 0 )
   {
      fEvt->print() ;  
   }  

   auto_ptr<HepMCProduct> BProduct(new HepMCProduct()) ;
   BProduct->addHepMCData( fEvt );
   e.put(BProduct);
    
   if ( fVerbosity > 0 )
   {
      // for testing purpose only
      //fEvt->print() ;  // for some strange reason, it prints NO event info
      // after it's been put into edm::Event...
      cout << " MultiParticleGunEtaEMC2 : Event Generation Done " << endl;
   }
}

DEFINE_FWK_MODULE(MultiParticleGunSourceEtaEMC2);
