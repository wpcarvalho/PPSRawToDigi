/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"
#include "TotemAnalysis/TotemNtuplizer/interface/Ntuplizer.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"


/**
 *\brief Saves event metadata (OptoRx BX numbers etc.) into the common ntuple.
 **/
class RawMetaDataNtuplizer : public Ntuplizer
{
  public:
    RawMetaDataNtuplizer(const edm::ParameterSet &ps);

    virtual ~RawMetaDataNtuplizer() {}

    virtual void DeclareConsumes(edm::EDAnalyzer *analyzer) override {}
    virtual void CreateBranches(const edm::EventSetup&, TTree *) override;
    virtual void FillEvent(const edm::Event&, const edm::EventSetup&) override;

  private:
    /// data to populate one `row' of the tree
    EventMetaData data;
    edm::InputTag fedInfosLabel;
};
