/****************************************************************************
*   Seyed Mohsen Etesami    
****************************************************************************/


#include "EventFilter/CTPPSRawToDigi/interface/DiamondVFATFrameCollection.h"

//----------------------------------------------------------------------------------------------------

using namespace std;

DiamondVFATFrameCollection::DiamondVFATFrameCollection()
{
}

//----------------------------------------------------------------------------------------------------

DiamondVFATFrameCollection::~DiamondVFATFrameCollection()
{
  data.clear();
}

//----------------------------------------------------------------------------------------------------

const DiamondVFATFrame* DiamondVFATFrameCollection::GetFrameByID(unsigned int ID) const
{
  // first convert ID to 12bit form
  ID = ID & 0xFFF;

  for (MapType::const_iterator it = data.begin(); it != data.end(); ++it)
    if (it->second.getChannelID() == ID)
      if (it->second.checkFootprint())
        return &(it->second);

  return NULL;
}

//----------------------------------------------------------------------------------------------------

const DiamondVFATFrame* DiamondVFATFrameCollection::GetFrameByIndex(DiamondFramePosition index) const
{
  MapType::const_iterator it = data.find(index);
  if (it != data.end())
    return &(it->second);
  else
    return NULL;
}

//----------------------------------------------------------------------------------------------------

DiamondVFATInterface::value_type DiamondVFATFrameCollection::BeginIterator() const
{
  MapType::const_iterator it = data.begin();
  return (it == data.end()) ? value_type(DiamondFramePosition(), NULL) : value_type(it->first, &it->second);
}

//----------------------------------------------------------------------------------------------------

DiamondVFATInterface::value_type DiamondVFATFrameCollection::NextIterator(const value_type &value) const
{
  if (!value.second)
    return value;

  MapType::const_iterator it = data.find(value.first);
  it++;

  return (it == data.end()) ? value_type(DiamondFramePosition(), NULL) : value_type(it->first, &it->second);
}

//----------------------------------------------------------------------------------------------------

bool DiamondVFATFrameCollection::IsEndIterator(const value_type &value) const
{
  return (value.second == NULL);
}
