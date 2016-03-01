/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
* $Id: RPNonParallelTrackCandidateFinder.cc,v 1.5.2.3 2009/11/18 17:42:14 jkaspar Exp $
* $Revision: 1.6.2.2 $
* $Date: 2009/12/09 15:34:58 $
*
****************************************************************************/

#ifndef RecoTotemRP_RPRecoDataFormats_RPRecognizedPatternsCollection_h
#define RecoTotemRP_RPRecoDataFormats_RPRecognizedPatternsCollection_h

#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatterns.h"

#include <map>

/**
 *\brief A collection (per RP) of pattern-recognition results.
 **/
class RPRecognizedPatternsCollection : public std::map<RPId, RPRecognizedPatterns>
{
};


#endif

