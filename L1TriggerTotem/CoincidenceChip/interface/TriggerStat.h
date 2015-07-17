#ifndef _L1TriggerTotemTriggerStat_H_
#define _L1TriggerTotemTriggerStat_H_

#include "L1TriggerTotem/CoincidenceChip/interface/PreTriggerStat.h"


class TriggerStat: public PreTriggerStat{
    public:
        TriggerStat(const TString& iname=""): PreTriggerStat(iname){
            fConditionANDTrigger     = 0;
            fConditionANDNoTrigger   = 0;
        }

         ~TriggerStat(){}

        unsigned int GetConditionNEvents() const{ return fConditionANDTrigger+fConditionANDNoTrigger;}

        unsigned int fConditionANDTrigger;
        unsigned int fConditionANDNoTrigger;
    ClassDef(TriggerStat,1)
};
#endif
