/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com)
*    
* $$RCSfile: DAQInformationT2.h,v $: $
* $Revision: 1.1.2.1 $
* $Date: 2009/11/06 16:03:06 $
*
****************************************************************************/
#ifndef DAQInformationT2Braw__
#define DAQInformationT2Braw__

#include "TotemRawDataLibrary/DataFormats/interface/FramePosition.h"
#include "TotemCondFormats/DAQInformation/interface/VFATRegisters.h"

#include <map>

using namespace Totem;


/**
 *\brief Container for T2 related DAQ data (mappings etc.).
 */
class DAQInformationT2_b
{
  public:
    /// 5 digit integer
    typedef unsigned int VFATID;

    /// 3 digit integer
    typedef unsigned int T2ID;

    /// mappings for readout chips
    //readoutPositionToId[position_t2] = symId;   //position_t2=slink
    //readoutIdToRegisters[symId] = VFATRegisters(ID_t2);

    std::map<FramePosition, VFATID> readoutPositionToId;
    std::map<VFATID, VFATRegisters> readoutIdToRegisters;

    /// mapping for coincidence VFATs
    std::map<FramePosition, T2ID> coincidencePositionToId;
    std::map<T2ID, VFATRegisters> coincidenceIdToRegisters;

    //void Reset();
};

#endif
