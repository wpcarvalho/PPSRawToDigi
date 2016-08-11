/****************************************************************************
 *  Syed Mohsen Etesami
 ****************************************************************************/

#include "FWCore/Utilities/interface/typelookup.h"

#include "CondFormats/CTPPSReadoutObjects/interface/TotemDAQMappingDiamond.h"

using namespace std;

//----------------------------------------------------------------------------------------------------

std::ostream& operator << (std::ostream& s, const DiamondVFATInfo &vi)
{
  if (vi.type == DiamondVFATInfo::data)
    s << "type=data, ";
  else
    s << "type=  CC, ";

  s << vi.symbolicID << ", hw id=0x" << hex << vi.hwID << dec;

  return s;
}

//----------------------------------------------------------------------------------------------------

void TotemDAQMappingDiamond::insert(const DiamondFramePosition &fp, const DiamondVFATInfo &vi)
{
  auto it = VFATMapping.find(fp);  
  if (it != VFATMapping.end())
    {
      cerr << "WARNING in DAQMapping::Insert > Overwriting entry at " << fp << ". Previous: " << endl 
	   << "    " << VFATMapping[fp] << "," << endl << "  new: " << endl << "    " << vi << ". " << endl;
    }

  VFATMapping[fp] = vi;
}

//----------------------------------------------------------------------------------------------------

TYPELOOKUP_DATA_REG(TotemDAQMappingDiamond);
