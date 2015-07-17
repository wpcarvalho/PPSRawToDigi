#ifndef Totem_TriggerData_include_h_
#define Totem_TriggerData_include_h_

#include <TObject.h>

class TriggerData : public TObject {
 public:
  TriggerData();
  

  unsigned char type;
  unsigned int event_num;
  unsigned int bunch_num;
  unsigned int src_id;
  unsigned int orbit_num;
  unsigned char revision_num;
  unsigned int run_num;
  unsigned int trigger_num;
  unsigned int inhibited_triggers_num;
  unsigned int input_status_bits;

  ClassDef(TriggerData, 10023);
};


#endif
