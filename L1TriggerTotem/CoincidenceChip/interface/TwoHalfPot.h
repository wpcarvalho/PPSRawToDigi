// -*- C++ -*-
//
// Package:    L1TriggerTotem
// Class:      RPTriggerAnalyzer
//
// Original Author:  Jiri Prochazka
//         Created:  Mon Mar  1 15:34:46 CET 2010
// $Id$

#ifndef _L1TriggerTotemTwoHalfPot_H_
#define _L1TriggerTotemTwoHalfPot_H_
// #include "L1TriggerTotem/CoincidenceChip/interface/TriggerStats.h"

// #include "TObject.h"
// TOTEM
// #include "TotemCondFormats/DAQInformation/interface/DAQInformationRP.h"
// #include "L1TriggerTotem/CoincidenceChip/interface/PreTriggerStat.h"

// #include "L1TriggerTotem/CoincidenceChip/interface/TriggerStats.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStat.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStatCollection.h"

namespace edm{
    class EventSetup;
    class Event;
    class ParameterSet;
}

class HalfPot;
// class DAQInformationRP;

class TwoHalfPot{
// class TwoHalfPot: public TObject{
    public:
        // TString fName;
        // TwoHalfPot():fConfig(*(new edm::ParameterSet()) ){} 
        TwoHalfPot(unsigned int irpID, ERawOrSimu iType, const edm::ParameterSet& conf, unsigned int verb = 0);

        ~TwoHalfPot(){};

        unsigned int GetDecRPIdFull() const;
        HalfPot& GetHalfPot(EDetOrientation orient);
        const HalfPot& GetHalfPot(EDetOrientation orient) const;
        TString GetType() const;
        bool TriggerSignal_UandV() const; 
        bool TriggerSignal_UorV() const;
        bool TriggerSignal_UorV1sectorOnly() const;
        bool TriggerSignal_UandV1sectorOnly() const;
        bool TriggerSignal() const;
        void beginJob();
        void analyze(const edm::Event& event, const edm::EventSetup& iSetup);
        void endJob();

        void analyzePotUnderCondition(const TString& conditionName, bool condition);

        ERawOrSimu fType;
       // unsigned int fNEvents;

        TriggerStatCollection<TriggerStat> fConditionalTriggerStat;

        const edm::ParameterSet& fConfig; //!
        unsigned int verbosity;
        // bool         fUseValidTrackCondition;

        HalfPot& fOddHalfPot; //!
        HalfPot& fEvenHalfPot; //!
       // TriggerStat fTrigStat;
     //   ClassDef(TwoHalfPot,1)
};

#endif
