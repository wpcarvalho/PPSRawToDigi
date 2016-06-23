/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Hubert Niewiadomski
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/


#include "TotemAnalysis/TotemNtuplizer/interface/Ntuplizer.h"
#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProton.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPReconstructedProtonPair.h"

#include "DataFormats/TotemDigi/interface/TotemRPDigi.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "DataFormats/Common/interface/Wrapper.h"
#include "fastjet/ClusterSequence.hh"

class PrimaryProton
{
  public: 
    HepMC::FourVector vertex;
    HepMC::FourVector momentum;
    bool found;
};

/**
 *\brief Saves RP reconstruction data into the common ntuple.
 **/
class RPNtuplizer : public Ntuplizer
{
  public:
    RPNtuplizer(const edm::ParameterSet&);
    virtual ~RPNtuplizer() {}
    
    virtual void DeclareConsumes(edm::EDAnalyzer *analyzer) override;
    virtual void CreateBranches(const edm::EventSetup&, TTree *) override;
    virtual void FillEvent(const edm::Event&, const edm::EventSetup&) override;

    void SetOpticsConfig(const BeamOpticsParams & opt_cfg) {BOPar_ = opt_cfg;}
    
  private:
    typedef std::map<int, const RPReconstructedProton*> reconstructed_prot_map_type;
    
    bool FindReconstrucedProtons(const edm::Event& e);

    /// The purpose of this function is to read most forward primary protons for given event
    bool FindSimulatedProtons(const edm::Event& e);

    /// The purpose of this function is to read most forward primary protons' vertices for given event
    bool FindSimulatedProtonsVertex(const edm::Event& e);
    bool FindReconstrucedProtonPair(const edm::Event& e);

    /**
    * This function computes the parameters of simulated proton that will be saved in root file.
    * proton - can be left_prim_prot_ or right_prim_prot_
    * simulatedProtonNumber - left_prim_prot_ - 0, right_prim_prot_ - 1
    */
    void FindParametersOfSimulatedProtons(PrimaryProton proton, int simulatedProtonNumber);
    
    BeamOpticsParams BOPar_;
    bool optics_valid_;
    int Verbosity_;
    std::string RootFileName_;
    
    std::map<unsigned int, RPRootDumpTrackInfo> track_info_;
    std::map<unsigned int, RPRootDumpDigiInfo> digi_info_;
    std::map<unsigned int, RPRootDumpPatternInfo > patterns_info_;

    std::map<unsigned int, std::vector<RPRootDumpTrackInfo> > multi_track_info_;
    std::map<unsigned int, RPRootDumpReconstructedProton> rec_pr_info_;
    std::map<unsigned int, RPRootDumpReconstructedProton> sim_pr_info_;
    RPRootDumpReconstructedProtonPair rec_pr_pair_info_;
    
    RPReconstructedProtonPair reconstructed_proton_pair_;
    reconstructed_prot_map_type reconstructed_protons_;
    int right_rec_prot_fund_;
    int left_rec_prot_fund_;
    int file_number_; 
    int prev_event_;
    
    std::string modulLabelSimu_;
    std::string productLabelSimu_;

    PrimaryProton right_prim_prot_;
    PrimaryProton left_prim_prot_;
    
    /// this flag indicates whether we want to save primary protons information in root file
    bool primaryProtons; 

    /// whether we want to store information about digi or not
    bool includeDigi;
    
    /// whether we want to store the pattern recognition result(s)
    bool includePatterns;
    bool primaryJets_;
    std::string primaryJetsInstance_;
    std::string primaryJetsLabel_;
    RPRootDumpJetInfo MCjets_;
    RPRootDiffMassInfo diff_mass_info_;

    edm::InputTag rpFittedTrackCollectionLabel;
    edm::InputTag rpMulFittedTrackCollectionLabel;
    edm::InputTag rpStripDigiSetLabel;
    edm::InputTag rpDigClusterLabel;
    edm::InputTag rpUVPatternLabel;
    edm::InputTag rpReconstructedProtonCollectionLabel;
    edm::InputTag rpReconstructedProtonPairCollectionLabel;
};
