/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataNtuplizer.h"

#include "DataFormats/TotemDigi/interface/TotemTriggerCounters.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TTree.h"

#include <map>

using namespace std;
using namespace edm;

ClassImp(TriggerData)

//----------------------------------------------------------------------------------------------------

TriggerDataNtuplizer::TriggerDataNtuplizer(const edm::ParameterSet &ps) : Ntuplizer(ps),
    triggerCountersLabel(ps.getParameter<edm::InputTag>("TriggerCountersLabel"))
{
}

//----------------------------------------------------------------------------------------------------

void TriggerDataNtuplizer::CreateBranches(const edm::EventSetup&, TTree *tree)
{
  tree->Branch("trigger_data.", &data);
}

//----------------------------------------------------------------------------------------------------

void TriggerDataNtuplizer::FillEvent(const edm::Event &event, const edm::EventSetup &es)
{
  Handle< TotemTriggerCounters > input;
  event.getByLabel(triggerCountersLabel, input);

  data.type = input->type;
  data.event_num = input->event_num;
  data.bunch_num = input->bunch_num;
  data.src_id = input->src_id;
  data.orbit_num = input->orbit_num;
  data.revision_num = input->revision_num;
  data.run_num = input->run_num;
  data.trigger_num = input->trigger_num;
  data.inhibited_triggers_num = input->inhibited_triggers_num;
  data.input_status_bits = input->input_status_bits;
}
