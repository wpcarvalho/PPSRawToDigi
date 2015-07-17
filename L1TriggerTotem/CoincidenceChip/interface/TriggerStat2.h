#ifndef _L1TriggerTotemTriggerStat2_H_
#define _L1TriggerTotemTriggerStat2_H_

#include "L1TriggerTotem/CoincidenceChip/interface/PreTriggerStat.h"

#if 1
class TriggerStat2: public PreTriggerStat{
    public:
        TriggerStat2(const TString& iname=""):PreTriggerStat(iname){
            fConditionANDRawSimuCCoutputTheSame            = 0;
            fConditionANDRawSimuCCoutputNotTheSame         = 0;
            fConditionANDRawSimuCCoutputNotTheSameOdd      = 0;
            fConditionANDRawSimuCCoutputNotTheSameEven     = 0;
            fConditionANDRawSimuCCoutputNotTheSameEvenOdd  = 0;
            fConditionANDTriggerRaw0Simu0                  = 0;
            fConditionANDTriggerRaw1Simu0                  = 0;
            fConditionANDTriggerRaw0Simu1                  = 0;
            fConditionANDTriggerRaw1Simu1                  = 0;
        }

         ~TriggerStat2(){}

        // unsigned int GetConditionNEvents() const{ return fConditionANDTrigger+fConditionANDNoTrigger;}
        // unsigned int GetNoConditionNEvents() const{ return fNoConditionANDTrigger+fNoConditionANDNoTrigger;}
        // unsigned int GetNoConditionNEvents() const{ return fNoConditionANDTrigger+fNoConditionANDNoTrigger;}
//        unsigned int GetConditionAndTriggers() const{ return fConditionANDRawSimuCCoutputTheSame + fConditionANDRawSimuCCoutputNotTheSame;}

        unsigned int fConditionANDRawSimuCCoutputTheSame;
        unsigned int fConditionANDRawSimuCCoutputNotTheSame;
        unsigned int fConditionANDRawSimuCCoutputNotTheSameOdd;
        unsigned int fConditionANDRawSimuCCoutputNotTheSameEven;
        unsigned int fConditionANDRawSimuCCoutputNotTheSameEvenOdd;
        unsigned int fConditionANDTriggerRaw0Simu0;
        unsigned int fConditionANDTriggerRaw1Simu0;
        unsigned int fConditionANDTriggerRaw0Simu1;
        unsigned int fConditionANDTriggerRaw1Simu1;
        ClassDef(TriggerStat2,1)
};
#endif
#endif
