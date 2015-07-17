#ifndef TotemRPValidation_InelasticReconstructionValidation_Inelastic2ArmReconstructionValidation_h
#define TotemRPValidation_InelasticReconstructionValidation_Inelastic2ArmReconstructionValidation_h

#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h"
#include "HepMC/GenEvent.h"
#include "TotemRPValidation/BaseValidationClasses/interface/BaseCollectionManager.h"
#include "TotemRPValidation/InelasticReconstructionValidation/interface/RPProtonReconstructionInfo.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TotemRPValidation/InelasticReconstructionValidation/interface/PoolsReconstructionInfo.h"
#include "TotemRPValidation/ParamMADRefTransport/interface/ParamMADRefTransport.h"
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include "TRandom2.h"
#include "TotemRPValidation/ValidationTools/interface/ReconstructionProfile.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPair.h"
#include "TotemRPValidation/InelasticReconstructionValidation/interface/RPProtonPairReconstructionInfo.h"
#include <memory>

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}


class Inelastic2ArmReconstructionValidation : public edm::EDAnalyzer
{
  public:
    explicit Inelastic2ArmReconstructionValidation(const edm::ParameterSet&);
    ~Inelastic2ArmReconstructionValidation();
    
  private:
    typedef std::map<int, const HepMC::GenParticle*> primary_prot_map_type;
    typedef std::map<int, const RPReconstructedProton*> reconstructed_prot_map_type;
    virtual void beginRun(edm::Run const&, edm::EventSetup const&);
    virtual void analyze(const edm::Event&, const edm::EventSetup&);
    virtual void endJob();
    
    bool FindPrimaryProtons(const edm::Event& e);
    bool FindPrimaryVertex(const edm::Event& e);
    bool FindReconstrucedProtons(const edm::Event& e);
    void WriteHistograms(const std::string &root_file_name);
    void FillReferenceHistograms();
    void FillResidualHistograms();
    void FillRPPools();
    
    inline bool RightArm(RPId id) const {return id>=100;}
    inline bool LeftArm(RPId id) const {return id<100;}
    
    std::auto_ptr<RPProtonPairReconstructionInfo> prot_pair_histograms_;
    std::auto_ptr<BaseCollectionManager<PoolsReconstructionInfo, int, edm::ParameterSet> > residual_histograms_;
    std::auto_ptr<ParamMADRefTransport> mad_param_transport_;
    
    primary_prot_map_type primary_protons_;
    primary_prot_map_type primary_vertex_;
    RPReconstructedProtonPair reconstructed_protons_;

    edm::ParameterSet conf_; 
    std::string hist_file_name_;
    std::string HepMCProductLabelVertex_;
    std::string HepMCProductLabelProton_;
    std::string HepMCVertexModuleName_;
    std::string HepMCProtonModuleName_;
    int verbosity_;
    TRandom2 rand_gen_;
    BeamOpticsParams BOPar_;
    bool SDValidation_;
    edm::InputTag rpReconstructedProtonPairCollectionLabel;
};

#endif
