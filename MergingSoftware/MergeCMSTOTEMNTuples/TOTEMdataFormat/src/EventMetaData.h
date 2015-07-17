#ifndef _include_EventMetaData_h_
#define _include_EventMetaData_h_

#include <TObject.h>
#include <vector>

class EventMetaData : public TObject {

 public:
  EventMetaData();
  
  ULong_t run_no; // ULong_t unsigned long 
  ULong_t event_no;
  ULong_t daq_event_number;
  ULong64_t timestamp; //  ULong64_t, unsigned long long
  std::vector<unsigned int> optoRx_Id;
  std::vector<unsigned int> optoRx_BX;
  std::vector<unsigned int> optoRx_LV1;

  ClassDef(EventMetaData, 10023);
};

#endif
