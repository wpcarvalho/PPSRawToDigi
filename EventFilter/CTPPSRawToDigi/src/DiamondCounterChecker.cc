/****************************************************************************
 *  Seyed Mohsen Etesami
 ****************************************************************************/

#include "EventFilter/CTPPSRawToDigi/interface/DiamondCounterChecker.h"

//-------------------------------------------------------------------------------------------------

using namespace std;

//-------------------------------------------------------------------------------------------------

void DiamondCounterChecker::Fill(word counter, DiamondFramePosition fr)
{ 
  pair<CounterMap::iterator, bool> ret;
  
  vector<DiamondFramePosition> list;
  list.push_back(fr);
  ret = relationMap.insert(pair<word, vector<DiamondFramePosition> >(counter, list));
  if (ret.second == false)
    relationMap[counter].push_back(fr);
}
