#ifndef cluster_entry_h
#define cluster_entry_h

#include <vector>

class cluster_entry
{
 public:
  cluster_entry() {rad_coord=0; ang_coord=0;}
  cluster_entry(unsigned short rad, unsigned short ang) {rad_coord=rad; ang_coord=ang;}
    
  unsigned short rad_coord; 
  unsigned short ang_coord;  //{strip/row, sector/col}
  
  //ClassDef(cluster_entry,1)
};  

#endif
