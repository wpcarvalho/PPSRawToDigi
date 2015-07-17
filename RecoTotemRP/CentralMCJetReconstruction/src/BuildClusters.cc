#include <HepMC/GenEvent.h>
#include "RecoTotemRP/CentralMCJetReconstruction/interface/BuildClusters.h"


//#include "RecoTotemRP/CentralMCJetReconstruction/interface/Selector.hh"


void BuildClusters::AddParticles(HepMC::GenEvent *evt)
{
  for (HepMC::GenEvent::particle_const_iterator particle = evt->particles_begin(); particle != evt->particles_end(); ++particle)
  {
    if( (*particle)->is_undecayed())
    {
      HepMC::FourVector P = (*particle)->momentum();
      particles_.push_back( fastjet::PseudoJet(P.px(), P.py(), P.pz(), P.e()) );
      
      int pdg = (*particle)->pdg_id ();
      double pseudo_rap = P.eta();
      
      bool outgoing_proton = (pdg==2212 && fabs(P.pz()/beamEnergy_)>0.5 && fabs(pseudo_rap)> 5);
      if(!outgoing_proton)
      {
        centr_mass_info_.AddParticle(P.px(), P.py(), P.pz(), P.e());
      }
    }
  }
}


vector<fastjet::PseudoJet> BuildClusters::GetInclusiveJets(fastjet::JetAlgorithm jet_alg, 
  double R, double min_pt, double min_eta, double max_eta)
{
  fastjet::JetDefinition jet_def(jet_alg, R);
  fastjet::ClusterSequence cs(particles_, jet_def);
  
  vector<fastjet::PseudoJet> jets = sorted_by_pt( cs.inclusive_jets(min_pt) );
    
  vector<fastjet::PseudoJet> central_jets;
  
  for(unsigned int i=0; i<jets.size(); i++)
  {
    if(jets[i].rap()<=max_eta && jets[i].rap()>=min_eta)
    {
      central_jets.push_back(jets[i]);
    }
  }
  
  
//  return central_jets;
  return central_jets;
}
