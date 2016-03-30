/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataNtuplizer.h"

#include "DataFormats/TotemRawData/interface/TotemRawEvent.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "TTree.h"

#include <map>

using namespace std;
using namespace edm;

ClassImp(TriggerData)

//----------------------------------------------------------------------------------------------------

TriggerDataNtuplizer::TriggerDataNtuplizer(const edm::ParameterSet &ps) : Ntuplizer(ps),
    rawEventLabel(ps.getParameter<edm::InputTag>("RawEventLabel"))
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
  Handle< TotemRawEvent > input;
  event.getByLabel(rawEventLabel, input);

  auto &td = input->getTriggerData();
  data.type = td.type;
  data.event_num = td.event_num;
  data.bunch_num = td.bunch_num;
  data.src_id = td.src_id;
  data.orbit_num = td.orbit_num;
  data.revision_num = td.revision_num;
  data.run_num = td.run_num;
  data.trigger_num = td.trigger_num;
  data.inhibited_triggers_num = td.inhibited_triggers_num;
  data.input_status_bits = td.input_status_bits;
}
