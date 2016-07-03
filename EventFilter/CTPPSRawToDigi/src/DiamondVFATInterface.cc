/****************************************************************************
*   Seyed Mohsen Etesami    
****************************************************************************/


#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATInterface.h"

//----------------------------------------------------------------------------------------------------
    
const DiamondVFATFrame* DiamondVFATInterface::GetFrameByIndexID(DiamondFramePosition index, unsigned int ID)
{
  const DiamondVFATFrame* returnframe = GetFrameByIndex(index);
  if (returnframe == NULL)
    return NULL;
  return (returnframe->getChannelID() == (ID & 0xFFF)) ? returnframe : NULL;
}
