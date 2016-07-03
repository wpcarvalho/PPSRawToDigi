/****************************************************************************
* Seyed Mohsen Etesami    
****************************************************************************/

#include "DataFormats/CTPPSDigi/interface/DiamondVFATStatus.h"

#include <ostream>

std::ostream& operator << (std::ostream& s, const DiamondVFATStatus &st)
{
  return s
      << "miss=" << st.status[0]
      << ",ID=" << st.status[1]
      << ",foot=" << st.status[2]
      << ",CRC=" << st.status[3]
      << ",EC=" << st.status[4]
      << ",BC=" << st.status[5]
      << ",fm=" << st.status[6]
      << ",pm=" << st.status[7];
}
