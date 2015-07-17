#ifndef _L1TriggerTotemRawVsSimuPotComparator_H_
#define _L1TriggerTotemRawVsSimuPotComparator_H_

// system

// #include <boost/ptr_container/ptr_vector.hpp>



// ROOT 
//#include "TFile.h"

#include "TH1D.h"
#include "TH2D.h"

//#include "TCanvas.h"
//#include "TSystem.h"
//#include "TTree.h"

// CMSSW 
// #include "DataFormats/Common/interface/DetSetVector.h"
// #include "FWCore/ParameterSet/interface/ParameterSet.h"


// TOTEM
// #include "TotemCondFormats/DAQInformation/interface/DAQInformationRP.h"
// #include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
// #include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"
// #include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"

#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStat2.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStatCollection.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TwoHalfPot.h"
// #include "L1TriggerTotem/CoincidenceChip/interface/TriggerStats.h"

class RawVsSimuPotComparator: public TNamed{
    public:
        // RawVsSimuPotComparator(unsigned int rpID = 0, const TString& iname="_RawVsSimuPotComparator_", int verb=0, const edm::ParameterSet* conf=0);
        RawVsSimuPotComparator(unsigned int rpID = 0, int verb=0, const edm::ParameterSet* conf=0);
        ~RawVsSimuPotComparator(){};
#if 0        
        void AddPotCondition(const TString& condName){
             fPotRaw.AddIfNotExists(new TwoHalfPot(condName,GetDecRPIdFull(),kRaw ,fConfig,verbosity));
            fPotSimu.AddIfNotExists(new TwoHalfPot(condName,GetDecRPIdFull(),kSimu,fConfig,verbosity));
            fConditionalTriggerStat.AddIfNotExists(condName);
        };

        TwoHalfPot& GetPotCondition(ERawOrSimu type, const TString& condName) {
           return *GetPot(type).AssertFind(condName);
        }

        const TwoHalfPot& GetPotCondition(ERawOrSimu type, const TString& condName) const {
           return *GetPot(type).AssertFind(condName);
        }
#endif

#if 1        
        TwoHalfPot& GetPot(ERawOrSimu type) {
            if(type == kRaw)
                return fPotRaw;
            else
                return fPotSimu;
        }

        const TwoHalfPot& GetPot(ERawOrSimu type) const{
            if(type == kRaw)
                return fPotRaw;
            else
                return fPotSimu;
        }
#endif

        // TString GetName() const{ return fName;}

        void beginJob();
        void analyze(const edm::Event& event, const edm::EventSetup& iSetup);
        void endJob();
        void endJob2();

        TString FillHist2D(const TString& conditionName,const TString& name,const edm::Event& event, EDetOrientation orient);
        void PrintEventTwoHalfPot() const;
        void analyzePotUnderCondition(const TString& conditionName, bool condition, const edm::Event& event, EDetOrientation orient);
        void FailingSectors();
        void StripsOnDistrib(const TString& conditionName ,const TString& name,const edm::Event& event, EDetOrientation orient,unsigned int iplane );
        void SectorsOnDistrib(const TString& conditionName ,const TString& name,const edm::Event& event, EDetOrientation orient,unsigned int iplane );
        void AddHist2D(TH2D* hist){
            TString orient=hist->GetName();
            if(orient.Contains("odd"))
                orient = "odd";
            else if(orient.Contains("even"))
                orient = "even";
            else
                orient = "";

            hist->SetTitle((TString)hist->GetName()+";Sector;Detector plane ("+orient+")");
            fHist2D.push_back(hist);
        }

        void AddHist1D(TH1D* hist){
            TString hname=hist->GetName();
            if(hname.Contains("StripsON"))
                hname = "Strip";
            else if(hname.Contains("SectorsON"))
                hname = "Sector";

            hist->SetTitle((TString)hist->GetName()+";"+hname);
            fHist1D.push_back(hist);
        }

        TH2D* FindHist2D(const TString& name){
            for(unsigned int i=0; i < fHist2D.size(); ++i){
                if(!strcmp(fHist2D[i]->GetName(),name)) return fHist2D[i];
            }

            //edm::LogInfo("CClogic") <<  "inputBits" << nevents;
   //          std::cout << "Can not find histogram named \"" << name << "\"" << std::endl;
            //throw cms::Exception("PPBckgAnalyzer2::FindHist2") << "Can not find histogram named " << name << std::endl;

   //         assert(0);
            return 0;
        }

        TH1D* FindHist1D(const TString& name){
            for(unsigned int i=0; i < fHist1D.size(); ++i){
                if(!strcmp(fHist1D[i]->GetName(),name)) return fHist1D[i];
            }

    //      std::cout << "Can not find histogram named \"" << name << "\"" << std::endl;
            //throw cms::Exception("PPBckgAnalyzer2::FindHist2") << "Can not find histogram named " << name << std::endl;

     //       assert(0);
            return 0;
        }

        unsigned int GetDecRPIdFull() const{
           // this implementation whould be better leads to problem with use of this function from python 
           // return value = 0 always because base class TotRPDetId is not saved (no dictionary)
           //  assert(GetPot(kRaw).GetDecRPIdFull() == GetPot(kSimu).GetDecRPIdFull());
           // return GetPot(kRaw).GetDecRPIdFull();

           // use alternative implementation 
          return RPNameToId(this->GetName());
        }
        unsigned int RPNameToId(const TString& name) const{

            if (!strcmp(name,"rp_45_220_nr_tp")) return 20;
            if (!strcmp(name,"rp_45_220_nr_bt")) return 21;
            if (!strcmp(name,"rp_45_220_nr_hr")) return 22;
            if (!strcmp(name,"rp_45_220_fr_hr")) return 23;
            if (!strcmp(name,"rp_45_220_fr_tp")) return 24;
            if (!strcmp(name,"rp_45_220_fr_bt")) return 25;
            if (!strcmp(name,"rp_56_220_nr_tp")) return 120;
            if (!strcmp(name,"rp_56_220_nr_bt")) return 121;
            if (!strcmp(name,"rp_56_220_nr_hr")) return 122;
            if (!strcmp(name,"rp_56_220_fr_hr")) return 123;
            if (!strcmp(name,"rp_56_220_fr_tp")) return 124;
            if (!strcmp(name,"rp_56_220_fr_bt")) return 125;
            std::cout << __func__ << "  bad input name =" << name << std::endl; 
            assert(0);
        }

#if 0
        unsigned int GetNEmptyPots() const{                
            std::cout << "Do we need this function??" << std::endl;
            assert(0);
            TString condName = "EmptyPot";
            assert(GetPot(kRaw).fConditionalTriggerStat.Find(condName));
            assert(GetPot(kRaw).fConditionalTriggerStat.Find(condName)->GetConditionNEvents() == GetPot(kSimu).fConditionalTriggerStat.Find(condName)->GetConditionNEvents() );
            return GetPot(kRaw).fConditionalTriggerStat.Find(condName)->GetConditionNEvents();
        }
#endif
        unsigned int GetRPNPlanesDividedByTwo() const;
        unsigned int GetRPPlanes() const;
        unsigned int GetNSectors() const;
        unsigned int GetNStrips() const;
        bool hasTheSameRawSimuCCoutputPattern() const;
        bool hasTheSameRawSimuCCoutputPattern(EDetOrientation orient) const;

        TriggerStatCollection<TriggerStat2> fConditionalTriggerStat;
    private:
       // TString fName;
        const edm::ParameterSet& fConfig; //!
        TwoHalfPot fPotRaw;
        TwoHalfPot fPotSimu;
        int verbosity;
//        edm::InputTag detTriggerLabel;
//        edm::InputTag fittedTrackCollectionLabel;
//        edm::InputTag stripDigiLabel;
        std::string detTriggerLabel;
        std::string fittedTrackCollectionLabel;
        std::string stripDigiLabel;
    public:
        unsigned int fNEmptyEvents;

        std::vector<TH1D*> fHist1D; 
        std::vector<TH2D*> fHist2D; 


       //  boost::ptr_vector<TH1D> fHist1D; //!
       //  boost::ptr_vector<TH2D> fHist2D; //!
     ClassDef(RawVsSimuPotComparator,1)
};

#endif
