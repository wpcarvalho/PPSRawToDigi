#include "L1TriggerTotem/CoincidenceChip/interface/TwoHalfPot.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStats.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
// #include "L1TriggerTotem/CoincidenceChip/interface/RPTriggerAnalyzer.h"

// ClassImp(TwoHalfPot)

TwoHalfPot::TwoHalfPot(unsigned int irpID, ERawOrSimu iType, const edm::ParameterSet& conf, unsigned int verb): 
           // fName(name),
            fType(iType), 
            fConfig(conf),
            verbosity(verb), 
            fOddHalfPot( *(new HalfPot(irpID,kOddU,iType,conf))),
            fEvenHalfPot( *(new HalfPot(irpID,kEvenV,iType,conf)))
          //  fTrigStat(name,true)
        {
            //   fNEvents               = 0;
            // fUseValidTrackCondition  = useValTrackCondition;
        }


void TwoHalfPot::beginJob(){
}

void TwoHalfPot::analyze(const edm::Event& event, const edm::EventSetup& iSetup){
	// if(1) cout << __func__ << " Sonda 0 type " << GetType() << endl;

    GetHalfPot(kOddU).analyze(event,iSetup);
    GetHalfPot(kEvenV).analyze(event,iSetup);

#if 0    
    fNEvents++;
	// Read real strips: input of CC
	edm::Handle< edm::DetSetVector<RPDetTrigger> > inputTrigger;
	event.getByType( inputTrigger );

    bool emptytest = true;
    for(edm::DetSetVector<RPDetTrigger>::const_iterator it = inputTrigger->begin(); it!= inputTrigger->end(); it++){
      for(edm::DetSet<RPDetTrigger>::const_iterator it2 = it->begin(); it2!= it->end(); it2++){
        if(this->GetDecRPIdFull()==TotRPDetId::RawToDecId(it2->GetDetId())/10) emptytest = false; 
      }
    }
   // cout << "#Event" << fNEvents << endl;
    analyzePotUnderCondition("EmptyPot",emptytest);
    analyzePotUnderCondition("NonEmptyPot",!emptytest);
    analyzePotUnderCondition("True",true);

    if(fConfig.getParameter<bool>("tracks")){
        // Read fitted tracks
        edm::Handle< RPFittedTrackCollection > fitTrCol; 
        event.getByType(fitTrCol);

        bool recoTrackTest = false;
        unsigned int  i=0;
        for(std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol->begin(); it != fitTrCol->end(); ++it) {
            unsigned int rpId=it->first;
            if(rpId!=this->GetDecRPIdFull()) continue; 
            const RPFittedTrack& fitTrack = it->second;
            // analyzePotUnderCondition("RecoTrack", fitTrack.IsValid());
            if(fitTrack.IsValid()) recoTrackTest = true;

            // there should be one valid or non valid track per event but no more 
            i++; assert(i==1);
        }

        analyzePotUnderCondition("RecoTrack", recoTrackTest); 
        analyzePotUnderCondition("NoRecoTrack", !recoTrackTest); 
    }

#endif

}

void TwoHalfPot::endJob(){
    if(verbosity){
        std::cout << "Pot " <<  TotRPDetId::RPName(this->GetDecRPIdFull())  << ", " << this->GetType() << " trigger" << std:: endl;
        for(unsigned int i = 0; i < fConditionalTriggerStat.size(); i++){
            const TriggerStat& trigStat = *fConditionalTriggerStat[i];

            // if(trigStat.fPrintConditionFlag){
            // if(strcmp(trigStat.GetName(),"NoTrue")){
              std::cout << "condition                = " << trigStat.GetName() << std::endl;
              std::cout << "NEvents satisfaing cond. = " << trigStat.GetConditionNEvents()    << std::endl;
              std::cout << "condition and Trigger    = " << trigStat.fConditionANDTrigger     << std::endl;
              std::cout << "condition and NoTrigger  = " << trigStat.fConditionANDNoTrigger   << std::endl;
              std::cout << std::endl;
            // }
        }
    }
}

unsigned int TwoHalfPot::GetDecRPIdFull() const{ assert(fOddHalfPot.GetDecRPIdFull()==fEvenHalfPot.GetDecRPIdFull()); return fOddHalfPot.GetDecRPIdFull();}

HalfPot& TwoHalfPot::GetHalfPot(EDetOrientation orient){
    if(orient == kOddU)
        return fOddHalfPot;
    else
        return fEvenHalfPot;
}

const HalfPot& TwoHalfPot::GetHalfPot(EDetOrientation orient) const{
  if(orient == kOddU)
      return fOddHalfPot;
  else
      return fEvenHalfPot;
}
 //     //  virtual std::string ClassName() const { return "TwoHalfPot"; }
 //     //  virtual const TString& GetName() const { return fName; }

TString TwoHalfPot::GetType() const{
    TString name;
    switch (fType) {
        case kRaw:  name = "Raw"; break;
        case kSimu: name = "Simu"; break;
        default: assert(0);
    }
    return name;
}

bool TwoHalfPot::TriggerSignal_UandV() const{
    // return true if at least one bit is ON in both bitsets
    return (fEvenHalfPot.fCCOutputBits.getBS().count() and fOddHalfPot.fCCOutputBits.getBS().count());
}

bool TwoHalfPot::TriggerSignal_UorV() const{
    // return true if at least one bit is ON in at least one bitset
    return (fEvenHalfPot.fCCOutputBits.getBS().count() or fOddHalfPot.fCCOutputBits.getBS().count());
}

bool TwoHalfPot::TriggerSignal_UorV1sectorOnly() const{
    // return true if and only if 1 bit (sector) is ON in U or V  bitset
    return (fEvenHalfPot.fCCOutputBits.getBS().count()==1 or fOddHalfPot.fCCOutputBits.getBS().count()==1);
} 

bool TwoHalfPot::TriggerSignal_UandV1sectorOnly() const{
    // return true if and only if 1 bit (sector) is ON in U and V  bitset
    return (fEvenHalfPot.fCCOutputBits.getBS().count()==1 and fOddHalfPot.fCCOutputBits.getBS().count()==1);
}

bool TwoHalfPot::TriggerSignal() const{
    TString triggerName = fConfig.getParameter<std::string>("trigger");
    // cout << " trigger = " << triggerName << endl;

    if(!strcmp(triggerName,"UandV")){
        return TriggerSignal_UandV();
    }
    else if(!strcmp(triggerName,"UorV")){
        return TriggerSignal_UorV();
    }
    else if(!strcmp(triggerName,"UandV1sectorOnly")){
        return TriggerSignal_UandV1sectorOnly();
    }
    else if(!strcmp(triggerName,"UorV1sectorOnly")){
        return TriggerSignal_UorV1sectorOnly();
    }
    else if(!strcmp(triggerName,"False")){
        return false;
    }


    cout << "trigger = " << triggerName << " is not defined!!!" << endl;
    assert(0);
    return false;
    // return TriggerSignalOR();
}

void TwoHalfPot::analyzePotUnderCondition(const TString& conditionName, bool condition) {
    // cout << "cond2=" << conditionName << endl;
    fConditionalTriggerStat.AddIfNotExists(conditionName);
    if (condition) {
        if(TriggerSignal())
            fConditionalTriggerStat.Find(conditionName)->fConditionANDTrigger++;
        else
            fConditionalTriggerStat.Find(conditionName)->fConditionANDNoTrigger++;
    }
}

