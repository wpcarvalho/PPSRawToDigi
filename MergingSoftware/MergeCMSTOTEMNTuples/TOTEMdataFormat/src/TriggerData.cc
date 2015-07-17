#include "TriggerData.h"

ClassImp(TriggerData);

TriggerData::TriggerData()
  : type(0),
    event_num(0),
    bunch_num(0),
    src_id(0),
    orbit_num(0),
    revision_num(0),
    run_num(0),
    trigger_num(0),
    inhibited_triggers_num(0),
    input_status_bits(0) {
}
