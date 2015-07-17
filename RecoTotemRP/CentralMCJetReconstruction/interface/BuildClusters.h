#ifndef RecoTotemRP_CentralMCJetReconstruction_BuildClusters_h_
#define RecoTotemRP_CentralMCJetReconstruction_BuildClusters_h_

#include <vector>
#include "HepMC/GenEvent.h"
#include "fastjet/ClusterSequence.hh"
#include "RecoTotemRP/RPRecoDataFormats/interface/CentralMassInfo.h"

using namespace std;
//using namespace fastjet;

class BuildClusters
{
  public:
    void ClearData()
    {
      particles_.clear();
      centr_mass_info_.clear();
    }
    
    BuildClusters(double beamEnergy) : beamEnergy_(beamEnergy) {}
    void AddParticles(HepMC::GenEvent *evt);
    vector<fastjet::PseudoJet> GetInclusiveJets(fastjet::JetAlgorithm jet_alg, double R, double min_pt, double min_eta, double max_eta);
    CentralMassInfo GetCentralMassInfo() {return centr_mass_info_;}
  
  private:
    vector<fastjet::PseudoJet> particles_;
    CentralMassInfo centr_mass_info_;
    double beamEnergy_;
};


#endif
