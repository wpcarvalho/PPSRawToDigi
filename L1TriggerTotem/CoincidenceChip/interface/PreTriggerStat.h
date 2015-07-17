#ifndef _L1TriggerTotemPreTriggerStat_H_
#define _L1TriggerTotemPreTriggerStat_H_

#include <cassert>

#include "TString.h" 
#include "TNamed.h" 

enum EDetOrientation {kEvenV = 0, kOddU = 1};
enum ERawOrSimu {kRaw, kSimu};

class PreTriggerStat: public TNamed{
//class PreTriggerStat{
    public:
        unsigned int fConditions;

//    private:
//        TString name;
//        bool fPrintConditionFlag;

    public:
       // virtual const char* ClassName(){ return "PreTriggerStat";}
        PreTriggerStat(const TString& iname=""):TNamed(iname,iname){
//            name                = iname;
          //  fPrintConditionFlag = printFlag;
            fConditions         = 0;
        }
         ~PreTriggerStat(){}
//        TString GetName() const{ return name;}
        unsigned int GetConditionNEvents() const{ return fConditions; }
        static TString RawOrSimuString(ERawOrSimu type){
            TString name;
            switch (type) {
                case kRaw:  name = "Raw"; break;
                case kSimu:   name = "Simu"; break;
                default: assert(0);
            }
            return name;
        }
    ClassDef(PreTriggerStat,1)
};

#endif
