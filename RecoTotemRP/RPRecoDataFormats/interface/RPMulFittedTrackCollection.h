/****************************************************************************
* This module is directly copied from 
* RPRecoDataFormats/RPFittedTrackCollection and expanded to hold 
* multiple fitted tracks.
*   
* Original Authors:
*   Hubert Niewiadomski (Hubert.Niewiadomski@cern.ch)
* Secondary Authors:
*   Zhang Zhengkui (zhang.zhengkui.fin@gmail.com)
* 
* $Id: RPMulFittedTrackCollection.h 0001 2010-08-20 17:00:00Z zhangz $
* $Revision: 0001 $
* $Date: 2010-08-20 17:00:00 +0100 (Fri, 20 Aug 2010) $
*
****************************************************************************/


#ifndef RecoTotemRP_RPRecoDataFormats_RPMulFittedTrackCollection_h
#define RecoTotemRP_RPRecoDataFormats_RPMulFittedTrackCollection_h


#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrack.h"

#include <vector>
#include <map>


class RPMulFittedTrackCollection : public std::map<RPId, std::vector<RPFittedTrack> >
{
};

class RPMulFittedTrackSetsCollection : public std::map<RPId, std::vector<std::vector<RPFittedTrack> > >
{
};


#endif
