/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "FWCore/Utilities/interface/InputTag.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "TotemAnalysis/TotemNtuplizer/interface/Ntuplizer.h"
#include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"

/**
 *\brief Saves the trigger data (from LoneG) to the TOTEM ntuple
 **/
class TriggerDataNtuplizer : public Ntuplizer
{
  public:
    TriggerDataNtuplizer(const edm::ParameterSet &ps);
    virtual ~TriggerDataNtuplizer() {}

    virtual void DeclareConsumes(edm::EDAnalyzer *analyzer) override {}
    virtual void CreateBranches(const edm::EventSetup&, TTree *) override;
    virtual void FillEvent(const edm::Event&, const edm::EventSetup&) override;

  private:
    edm::InputTag triggerCountersLabel;
    TriggerData data;
};
