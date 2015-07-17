#ifndef _L1TriggerTotemRPTriggerAnalyzerInfoCollector_H_
#define _L1TriggerTotemRPTriggerAnalyzerInfoCollector_H_

#include "TNamed.h"
#include "TString.h"
// #include "L1TriggerTotem/CoincidenceChip/interface/PotCollection.h"

#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStat.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStat2.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStatCollection.h"


class RPTriggerAnalyzerInfoCollector: public TNamed{

   public:
       RPTriggerAnalyzerInfoCollector(const TString& name="_def_name_RPTriggerAnalyzerInfoCollector_"): TNamed(name,name){
           fNEvents                 = 0;
           fNEmptyEvents            = 0;
           runName = "_def_run_name_";
       }
       ~RPTriggerAnalyzerInfoCollector(){}
        unsigned int fNEvents;
        unsigned int fNEmptyEvents;
        TString runName;
      //  void SetTriggerStatCollection(TriggerStatCollection<TriggerStat2>* coll){ assert(coll); fTriggerStatCollectionTriggerStat2 = coll;}
      //   TriggerStatCollection<TriggerStat2>* fTriggerStatCollectionTriggerStat2; //->
//        PotCollection  fPotCollection;
   ClassDef(RPTriggerAnalyzerInfoCollector,1)
};

#endif
