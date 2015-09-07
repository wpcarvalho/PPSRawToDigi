/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*    
* $$RCSfile: DAQInformationT2.h,v $: $
* $Revision: 4564 $
* $Date: 2011-07-29 13:48:25 +0200 (Pt, 29 lip 2011) $
*
****************************************************************************/
#ifndef DAQInformationT1_h
#define DAQInformationT1_h


#include "TotemRawDataLibrary/DataFormats/interface/FramePosition.h"
#include "TotemCondFormats/DAQInformation/interface/VFATRegisters.h"

#include <map>

using namespace Totem;


/**
 *\brief Container for T1 related DAQ data (mappings etc.).
 */
class DAQInformationT1
{
  public:
    /// 5 digit integer
    typedef unsigned int VFATID;

    /// 3 digit integer
    typedef unsigned int T1ID;

    /// mappings for readout chips
    //readoutPositionToId[position_t1] = symId;   //position_t1=slink
    //readoutIdToRegisters[symId] = VFATRegisters(ID_t1);

    std::map<FramePosition, VFATID> readoutPositionToId;
    std::map<VFATID, VFATRegisters> readoutIdToRegisters;

    /// mapping for coincidence VFATs
//    std::map<FramePosition, T2ID> coincidencePositionToId;
//    std::map<T2ID, VFATRegisters> coincidenceIdToRegisters;

    void Reset();
};

#endif
