#ifndef TotemRPValidation_RPGeant4Validation_RP_PSIMHIT_DEBUG_INFO_
#define TotemRPValidation_RPGeant4Validation_RP_PSIMHIT_DEBUG_INFO_

#include "SimDataFormats/TrackingHit/interface/PSimHit.h"
  
class RPPSimHitDebugInfo : public PSimHit
{
  public:
    RPPSimHitDebugInfo(const Local3DPoint& entry, const Local3DPoint& exit, 
        float pabs, float tof, float eloss, int particleType,
        int detId, unsigned int trackId,
        float theta, float phi, unsigned short processType=0)
    : PSimHit(entry, exit, pabs, tof, eloss, particleType, detId, trackId, theta, 
        phi, processType) {}
        
    RPPSimHitDebugInfo(const Local3DPoint& entry, const Local3DPoint& exit, 
        float pabs, float tof, float eloss, int particleType,
        int detId, unsigned int trackId,
        float theta, float phi, unsigned short processType,
        Local3DPoint globalPosition, Local3DPoint primaryVertex, int RPDetPartId, 
        int parentId, int RPId, int copy_no=0)
    : PSimHit(entry, exit, pabs, tof, eloss, particleType, detId, trackId, theta, 
        phi, processType), theGlobalPosition(globalPosition), thePrimaryVertex(primaryVertex), 
        theRPDetPartId(RPDetPartId), theParentId(parentId), theRPId(RPId), 
        theCopyNo(copy_no){}
       
    const Local3DPoint & GetGlobalPosition() const {return theGlobalPosition;}
    const Local3DPoint & GetPrimaryVertex() const {return thePrimaryVertex;}
    int GetRPDetPartId() const {return theRPDetPartId;}
    int GetParentId() const {return theParentId;}
    int GetRPId() const {return theRPId;}
    int GetCopyNo() const {return theCopyNo;}

  private:
    Local3DPoint theGlobalPosition;
    Local3DPoint thePrimaryVertex;
    int theRPDetPartId;
    int theParentId;
    int theRPId;
    int theCopyNo;
};

#endif /*RP_PSIMHIT_DEBUG_INFO_*/
