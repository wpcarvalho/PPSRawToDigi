#include "L1TriggerTotem/CoincidenceChip/interface/RawVsSimuPotComparator.h"
#include "L1TriggerTotem/CoincidenceChip/interface/CoincidenceChipConfiguration.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStats.h"

// CMSSW
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/DetSetVector.h"
//#include "FWCore/Framework/interface/EDConsumerBase.h"
#include "FWCore/Utilities/interface/InputTag.h"

// TOTEM
// #include "TotemCondFormats/DAQInformation/interface/DAQInformationRP.h"
#include "RecoTotemRP/RPRecoDataFormats/interface/RPFittedTrackCollection.h"
#include "DataFormats/TotemRPDataTypes/interface/RPStripDigi.h"
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"

ClassImp(RawVsSimuPotComparator)

// RawVsSimuPotComparator::RawVsSimuPotComparator(unsigned int rpID, const TString& iname, int verb, const edm::ParameterSet* conf):
RawVsSimuPotComparator::RawVsSimuPotComparator(unsigned int rpID, int verb, const edm::ParameterSet* conf):
  TNamed(TotRPDetId::RPName(rpID),TotRPDetId::RPName(rpID)),
  fConfig(*conf),
  fPotRaw(rpID, kRaw, *conf, verb),
  fPotSimu(rpID, kSimu, *conf, verb),
  verbosity(verb),
  fNEmptyEvents(0)
{
	 detTriggerLabel =  conf->getParameter<edm::InputTag>("DetTriggerLabel").label();
	 fittedTrackCollectionLabel = conf->getParameter<edm::InputTag>("FittedTrackCollectionLabel").label();
	 stripDigiLabel = conf->getParameter<edm::InputTag>("StripDigiLabel").label();
    // check mapping, consitency test
    for(unsigned int i=20; i <=25; i++)
      assert(RPNameToId(TotRPDetId::RPName(i)) == i);

    for(unsigned int i=120; i <=125; i++)
      assert(RPNameToId(TotRPDetId::RPName(i)) == i);
}

void RawVsSimuPotComparator::beginJob(){

#if 0    
    vector<TString> modif;
    modif.push_back("_2planes");
    modif.push_back("_Not2planes");

    for(unsigned int i = 0; i < modif.size(); i++){
        AddHist2D(new TH2D(name+"RecoTrack_NotTheSame_even"+modif[i] , "" , GetNSectors() , -0.5 , GetNSectors()-0.5 , GetRPNPlanesDividedByTwo() , 0.5 , GetRPNPlanesDividedByTwo()+0.5));
        AddHist2D(new TH2D(name+"RecoTrack_NotTheSame_odd"+modif[i]  , "" , GetNSectors() , -0.5 , GetNSectors()-0.5 , GetRPNPlanesDividedByTwo() , 0.5 , GetRPNPlanesDividedByTwo()+0.5));
        AddHist2D(new TH2D(name+"RecoTrack_TheSame_even"+modif[i]    , "" , GetNSectors() , -0.5 , GetNSectors()-0.5 , GetRPNPlanesDividedByTwo() , 0.5 , GetRPNPlanesDividedByTwo()+0.5));
        AddHist2D(new TH2D(name+"RecoTrack_TheSame_odd"+modif[i]     , "" , GetNSectors() , -0.5 , GetNSectors()-0.5 , GetRPNPlanesDividedByTwo() , 0.5 , GetRPNPlanesDividedByTwo()+0.5));
    }

    AddHist2D(new TH2D(name+"RecoTrack_Raw0Simu1_FailingSectors_even" , "" , GetNSectors() , -0.5 , GetNSectors()-0.5 , GetRPNPlanesDividedByTwo() , 0.5 , GetRPNPlanesDividedByTwo()+0.5));
    AddHist2D(new TH2D(name+"RecoTrack_Raw0Simu1_FailingSectors_odd"  , "" , GetNSectors() , -0.5 , GetNSectors()-0.5 , GetRPNPlanesDividedByTwo() , 0.5 , GetRPNPlanesDividedByTwo()+0.5));

    for(unsigned int plane=1; plane<=GetRPNPlanesDividedByTwo(); plane++){
        TString nplane="";
        nplane += plane;
        fHist1D.push_back(new TH1D(name+"RecoTrack_StripsON_odd_Plane"+nplane   , "" , GetNStrips()  , -0.5 , GetNStrips()-0.5));
        fHist1D.push_back(new TH1D(name+"RecoTrack_StripsON_even_Plane"+nplane  , "" , GetNStrips()  , -0.5 , GetNStrips()-0.5));
        fHist1D.push_back(new TH1D(name+"RecoTrack_SectorsON_odd_Plane"+nplane  , "" , GetNSectors() , -0.5 , GetNSectors()-0.5));
        fHist1D.push_back(new TH1D(name+"RecoTrack_SectorsON_even_Plane"+nplane , "" , GetNSectors() , -0.5 , GetNSectors()-0.5));

        fHist1D.push_back(new TH1D(name+"AllEvents_StripsON_odd_Plane"+nplane   , "" , GetNStrips()  , -0.5 , GetNStrips()-0.5));
        fHist1D.push_back(new TH1D(name+"AllEvents_StripsON_even_Plane"+nplane  , "" , GetNStrips()  , -0.5 , GetNStrips()-0.5));
        fHist1D.push_back(new TH1D(name+"AllEvents_SectorsON_odd_Plane"+nplane  , "" , GetNSectors() , -0.5 , GetNSectors()-0.5));
        fHist1D.push_back(new TH1D(name+"AllEvents_SectorsON_even_Plane"+nplane , "" , GetNSectors() , -0.5 , GetNSectors()-0.5));
    }
#endif    
}

void RawVsSimuPotComparator::analyze(const edm::Event& event, const edm::EventSetup& iSetup){
    //if(verbosity > 2) cout << __func__ << endl;
    GetPot(kRaw).analyze(event, iSetup);
    GetPot(kSimu).analyze(event, iSetup);

    DAQInformationRP::RPID chosenRPID = GetDecRPIdFull();
    if(verbosity > 2) cout << __func__ << "  id=" << chosenRPID << "  label=" << GetName() << endl;

        // Read real strips: input of CC
    edm::Handle< edm::DetSetVector<RPDetTrigger> > inputTrigger;
    event.getByLabel(detTriggerLabel, inputTrigger );
    bool emptytest = true;
    for(edm::DetSetVector<RPDetTrigger>::const_iterator it = inputTrigger->begin(); it!= inputTrigger->end(); it++){
        for(edm::DetSet<RPDetTrigger>::const_iterator it2 = it->begin(); it2!= it->end(); it2++){
            if(this->GetDecRPIdFull()==TotRPDetId::RawToDecId(it2->GetDetId())/10) emptytest = false;
        }
    }
    if(emptytest) fNEmptyEvents++;


    // test whether the activity is less then VoutOfNP (at trigger level, 
	// CC settings e.g. "3 out of 5 planes in coincidence" - VoutOfNP CC block)
    bool lessThanXPlanesTest = false;
    unsigned int lessThanXPlanesTestCounterEvenV = 0;
    unsigned int lessThanXPlanesTestCounterOddU = 0;
    for(edm::DetSetVector<RPDetTrigger>::const_iterator it = inputTrigger->begin(); it!= inputTrigger->end(); it++){
		// cout << "new" << endl;
        if(this->GetDecRPIdFull()==TotRPDetId::RawToDecId((*it)[0].GetDetId())/10) {
		   if(TotRPDetId::RawToDecId((*it)[0].GetDetId()) % 2)
		     lessThanXPlanesTestCounterEvenV++;
		   else
		     lessThanXPlanesTestCounterOddU++;
		}
    }
    CoincidenceChipConfiguration ccConfig;
	ccConfig.configure(fConfig.getParameterSet("params").getParameterSet("coincidenceChipConfig"));
	unsigned int VoutOfNP = ccConfig.getV();
	 cout << "VoutOfNP = " <<  VoutOfNP << endl;
	 cout << "V=" << lessThanXPlanesTestCounterEvenV << "  U=" << lessThanXPlanesTestCounterOddU  << endl;
	if((lessThanXPlanesTestCounterEvenV < VoutOfNP) and (lessThanXPlanesTestCounterOddU < VoutOfNP)) lessThanXPlanesTest = true; 
	// cout << "lessThanThreePlanesTestCounter = " << lessThanXPlanesTestCounter << endl;
    // END OF lessThanXPlanesTest



    bool recoTrackTest = false;
    if(fConfig.getParameter<bool>("tracks")){
        // Read fitted tracks
        edm::Handle< RPFittedTrackCollection > fitTrCol;
        event.getByLabel(fittedTrackCollectionLabel, fitTrCol);

        unsigned int  i=0;
        for(std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol->begin(); it != fitTrCol->end(); ++it) {
            unsigned int rpId=it->first;
            if(rpId!=chosenRPID) continue; //cout << rpId << endl;
            //cout << "  RPID=" << rpId << endl;
            const RPFittedTrack& fitTrack = it->second;
            if(fitTrack.IsValid()) recoTrackTest = true;
            // cout << "RecoTrack" << endl;

            // there should be one and only one valid or non valid track per event but no more
            i++; assert(i==1);
        }

    }

  //  bool trigger_signal_test = GetPot(kRaw).TriggerSignal() and GetPot(kSimu).TriggerSignal();
 
    std::vector<EDetOrientation> orient;
    orient.push_back(kOddU);
    orient.push_back(kEvenV);

#if 1
    for(unsigned int i = 0; i < orient.size();i++){
        analyzePotUnderCondition("EmptyPot",emptytest,event,orient[i]);
        analyzePotUnderCondition("NonEmptyPot",!emptytest,event,orient[i]);
        analyzePotUnderCondition("VoutOfNP",!lessThanXPlanesTest,event,orient[i]);
        analyzePotUnderCondition("!VoutOfNP",lessThanXPlanesTest,event,orient[i]);
       // analyzePotUnderCondition("NonEmptyPot+Triggers",!emptytest and trigger_signal_test,event,orient[i]);
        analyzePotUnderCondition("True",true,event,orient[i]);
//        analyzePotUnderCondition("Triggers",trigger_signal_test, event, orient[i]);

        if(fConfig.getParameter<bool>("tracks")){
            analyzePotUnderCondition("RecoTrack", recoTrackTest,event,orient[i]);
//            analyzePotUnderCondition("NoRecoTrack", !recoTrackTest,event,orient[i]);
//            analyzePotUnderCondition("RecoTrack+Trigger", recoTrackTestandtrigger_signal_test,event,orient[i]);
        }
    }
#endif

    if(verbosity) PrintEventTwoHalfPot();
}

void RawVsSimuPotComparator::endJob(){
    if(verbosity){
        std::cout << "#############################################################" << endl;
        std::cout << "Pot " << TotRPDetId::RPName(GetDecRPIdFull()) << " stats: " << endl;
        std::cout << "NEmptyEvents = " << fNEmptyEvents << endl;
        cout << endl;
        for(unsigned int i = 0; i < fConditionalTriggerStat.size(); i++){
            const TriggerStat2& trigStat = *fConditionalTriggerStat[i];
            // const TString& condition = "condition";
            // std::cout << "Pot " << this->GetDecRPIdFull()  << ", " << this->GetType() << " trigger" << std:: endl;

            // if(strcmp(trigStat.GetName(),"NoTrue")){
            std::cout << "condition                                      = " << trigStat.GetName()                                     << std::endl;
            std::cout << "NEvents satisfaing condition                   = " << trigStat.fConditions                                   << std::endl;
            std::cout << "condition and RawSimuCCoutputTheSame           = " << trigStat.fConditionANDRawSimuCCoutputTheSame           << std::endl;
            std::cout << "condition and RawSimuCCoutputNotTheSame        = " << trigStat.fConditionANDRawSimuCCoutputNotTheSame        << std::endl;
            std::cout << "condition and RawSimuCCoutputNotTheSameEven    = " << trigStat.fConditionANDRawSimuCCoutputNotTheSameEven    << std::endl;
            std::cout << "condition and RawSimuCCoutputNotTheSameOdd     = " << trigStat.fConditionANDRawSimuCCoutputNotTheSameOdd     << std::endl;
            std::cout << "condition and RawSimuCCoutputNotTheSameEvenOdd = " << trigStat.fConditionANDRawSimuCCoutputNotTheSameEvenOdd << std::endl;
            std::cout << "condition and TriggerRaw0Simu0                 = " << trigStat.fConditionANDTriggerRaw0Simu0                 << std::endl;
            std::cout << "condition and TriggerRaw1Simu0                 = " << trigStat.fConditionANDTriggerRaw1Simu0                 << std::endl;
            std::cout << "condition and TriggerRaw0Simu1                 = " << trigStat.fConditionANDTriggerRaw0Simu1                 << std::endl;
            std::cout << "condition and TriggerRaw1Simu1                 = " << trigStat.fConditionANDTriggerRaw1Simu1                 << std::endl;
            std::cout << std::endl;
            //  }
        }

            GetPot(kRaw).endJob();
            std::cout << endl;
            GetPot(kSimu).endJob();
            std::cout << endl;
            std::cout << endl;
    }
}

void RawVsSimuPotComparator::endJob2(){
    assert(0);
#if 0
    if(verbosity){
        const edm::ParameterSet params = fConfig.getParameterSet("params");
        std::string runInfo            = params.getParameter<std::string>("runInfo");
        std::string runName            = params.getParameter<std::string>("runName");

        cout << " Run   POT   cond. N_t      TM_raw/N_t     TM_simu/N_t    TMsame/N_t     TM_notsame/N_t " << endl;
        for(unsigned int i = 0; i < fConditionalTriggerStat.size(); i++){
            const TriggerStat2& trigStat2 = fConditionalTriggerStat[i];
            const TriggerStat& trigStatRaw  = *GetPot( kRaw).fConditionalTriggerStat.AssertFind(trigStat2.GetName());
            const TriggerStat& trigStatSimu = *GetPot(kSimu).fConditionalTriggerStat.AssertFind(trigStat2.GetName());
          //  assert(GetPot(kRaw).fConditionalTriggerStat.Find(trigStat.GetName()));
          //  assert(GetPot(kSimu).fConditionalTriggerStat.Find(trigStat.GetName()));
            assert(trigStat2.GetConditionNEvents() == trigStatRaw.GetConditionNEvents());
            assert(trigStat2.GetConditionNEvents() == trigStatSimu.GetConditionNEvents());

            cout << runName << " &  " << GetDecRPIdFull() << " & " <<   trigStat2.GetName()  << " & " << trigStat2.GetConditionNEvents();
            if(trigStat2.GetConditionNEvents()){
                cout << " & " << trigStatRaw.fConditionANDTrigger/trigStat2.GetConditionNEvents();
                cout << " & " << trigStatSimu.fConditionANDTrigger/trigStat2.GetConditionNEvents();
                cout << " & " << trigStat2.fConditionANDRawSimuCCoutputNotTheSame/trigStat2.GetConditionNEvents() << "  \\\\ "<< endl;
            }else{
                cout << " & -" ;
                cout << " & -" ;
                cout << " & - \\\\" << endl;
            }

        }
    }
#endif
}

void RawVsSimuPotComparator::PrintEventTwoHalfPot() const{
   // if(GetPot(kRaw).TriggerSignal() != GetPot(kSimu).TriggerSignal())
   //     cout << " NOT the same triggers";
   // else
   //     cout << " the same triggers";
    cout << "trigger: Raw = " << GetPot(kRaw).TriggerSignal() << "  Simu = " << GetPot(kSimu).TriggerSignal() << endl;
    cout << "has the same SimuRawCCoutputPattern: " << hasTheSameRawSimuCCoutputPattern() <<  "  (OddU = " << hasTheSameRawSimuCCoutputPattern(kOddU) << "  EvenV = " << hasTheSameRawSimuCCoutputPattern(kEvenV) << ")" << endl;
    cout << endl;
    std::vector<EDetOrientation> orient;
    orient.push_back(kOddU);
    orient.push_back(kEvenV);

    for (unsigned int i = 0; i < orient.size();i++){
        // print input to CC
        std::cout << "CC input " << OrientedRPId::DetOrientation(orient[i]) << endl;
        // std::cout << RPCCInputBits::PrintCCInputBits(fPotCollection.FindHalfPot(fPotCollection[potIndex].id,orient[i])->fCCInputBits.GetBits());
        // std::cout << RPCCInputBits::PrintCCInputBits(GetPot(kRaw).GetHalfPot(orient[i]).fCCInputBits.GetBits());
        std::cout << GetPot(kRaw).GetHalfPot(orient[i]).PrintCCInputBits();

        // print output of CC
        std::cout << "CC output " << OrientedRPId::DetOrientation(orient[i]) << endl;
        std::cout << GetPot(kRaw).GetHalfPot(orient[i]).fCCOutputBits.getBS()  << "  RAW" << endl;

        CoincidenceChip::OutputBits outputBitsSimu  = GetPot(kSimu).GetHalfPot(orient[i]).fCCOutputBits.getBS();
        std::cout << outputBitsSimu  << "  Simu" << endl;
 
#if 0        
        // check if outputBitsSimu is realy equal to outputBitsSimu2 assert(0) otherwise
        // CoincidenceChip::OutputBits outputBitsSimu2 = fCChip.process(fPotCollection.FindHalfPot(fPotCollection[potIndex].id,orient[i])->fCCInputBits.GetBits());
        CoincidenceChip coincidenceChip;
        // configure the coincidence chip
        coincidenceChip.configure(fConfig.getParameterSet("params").getParameterSet("coincidenceChipConfig"));
        CoincidenceChip::OutputBits outputBitsSimu2 = coincidenceChip.process(GetPot(kRaw).GetHalfPot(orient[i]).fCCInputBits.GetBits());
        if(outputBitsSimu!=outputBitsSimu2){
            std::cout << outputBitsSimu2 << "  Simu2" << endl;
            assert(outputBitsSimu==outputBitsSimu2);
        }
#endif        

        std::cout << endl;
    }

}

bool RawVsSimuPotComparator::hasTheSameRawSimuCCoutputPattern(EDetOrientation orient) const{
  return GetPot(kRaw).GetHalfPot(orient).fCCOutputBits.getBS() == GetPot(kSimu).GetHalfPot(orient).fCCOutputBits.getBS();
}

bool RawVsSimuPotComparator::hasTheSameRawSimuCCoutputPattern() const{
  return  hasTheSameRawSimuCCoutputPattern(kEvenV) and hasTheSameRawSimuCCoutputPattern(kOddU);
}

void RawVsSimuPotComparator::analyzePotUnderCondition(const TString& conditionName, bool condition, const edm::Event& event, EDetOrientation orient) {
    // cout << "cond=" << conditionName << endl;

    if(orient==kEvenV){  // dont increase it twice...
      GetPot(kRaw).analyzePotUnderCondition(conditionName,condition);
      GetPot(kSimu).analyzePotUnderCondition(conditionName,condition);
    }

    fConditionalTriggerStat.AddIfNotExists(conditionName);
    if(condition){
        TriggerStat2& trigStat = *fConditionalTriggerStat.AssertFind(conditionName);
        if(orient==kEvenV)  // dont increase it twice...
            trigStat.fConditions++;

        for(unsigned int plane=1; plane <= GetRPNPlanesDividedByTwo(); plane++){
            StripsOnDistrib(conditionName,"",event,orient,plane);
            SectorsOnDistrib(conditionName,"",event,orient,plane);
        }

        // if(GetPot(kRaw).TriggerSignal() and GetPot(kSimu).TriggerSignal()){
        if(hasTheSameRawSimuCCoutputPattern()) {

            if(orient==kEvenV)  // dont increase it twice...
                trigStat.fConditionANDRawSimuCCoutputTheSame++;
            FillHist2D(conditionName,"TheSame_2planes",event,orient);
            FillHist2D(conditionName,"TheSame_Not2planes",event,orient);
        }else{
            if(orient==kEvenV)  // dont increase it twice...
                trigStat.fConditionANDRawSimuCCoutputNotTheSame++;

            // FailingSectors();

            if(verbosity){
                if(orient==kEvenV){  // dont increase it twice...
                    //  std::cout << " ###  " + conditionName+ "ANDTriggerNotTheSame " << trigStat.fConditionANDRawSimuCCoutputNotTheSame << endl;
                    // PrintEventTwoHalfPot(event);
                }
            }

            if(!hasTheSameRawSimuCCoutputPattern(orient) and  hasTheSameRawSimuCCoutputPattern(OrientedRPId::SecondDetOrientation(orient))){
                if(orient==kEvenV)
                    trigStat.fConditionANDRawSimuCCoutputNotTheSameEven++;
                else
                    trigStat.fConditionANDRawSimuCCoutputNotTheSameOdd++;

                FillHist2D(conditionName,"NotTheSame_2planes",event,orient);
                FillHist2D(conditionName,"NotTheSame_Not2Planes",event,orient);
            }

            if(!hasTheSameRawSimuCCoutputPattern(kEvenV) and !hasTheSameRawSimuCCoutputPattern(kOddU)){
                if(orient==kEvenV)  // dont increase it twice...
                    trigStat.fConditionANDRawSimuCCoutputNotTheSameEvenOdd++;
                // if(verbosity) std::cout << "rare event" << endl;
            }

        }

        if(orient==kEvenV){  // dont increase it twice...
          if((GetPot(kRaw).TriggerSignal()==0) and (GetPot(kSimu).TriggerSignal()==0)) trigStat.fConditionANDTriggerRaw0Simu0++;
          if((GetPot(kRaw).TriggerSignal()==1) and (GetPot(kSimu).TriggerSignal()==0)) trigStat.fConditionANDTriggerRaw1Simu0++;
          if((GetPot(kRaw).TriggerSignal()==0) and (GetPot(kSimu).TriggerSignal()==1)) trigStat.fConditionANDTriggerRaw0Simu1++;
          if((GetPot(kRaw).TriggerSignal()==1) and (GetPot(kSimu).TriggerSignal()==1)) trigStat.fConditionANDTriggerRaw1Simu1++;
        }

        // if(!strcmp(conditionName,"NonEmptyEvent") and(GetPot(kRaw).TriggerSignal()==0) and (GetPot(kSimu).TriggerSignal()!=0)){
        if((GetPot(kRaw).TriggerSignal()==0) and (GetPot(kSimu).TriggerSignal()!=0)){
            FillHist2D(conditionName,"Raw0Simu1_FailingSectors",event,orient);
        }

    }

    }

TString RawVsSimuPotComparator::FillHist2D(const TString& conditionName, const TString& name, const edm::Event& event, EDetOrientation orient) {
    const TString& hname = TString(GetName())+"_"+conditionName+"_"+name+"_"+OrientedRPId::DetOrientation(orient);
    TH2D* hist = FindHist2D(hname);
     if(!hist) {
         AddHist2D(new TH2D(hname, "" , GetNSectors() , -0.5 , GetNSectors()-0.5 , GetRPNPlanesDividedByTwo() , 0.5 , GetRPNPlanesDividedByTwo()+0.5));
         hist = fHist2D.back();
     }

    const CoincidenceChip::InputBits bs  = GetPot(kRaw).GetHalfPot(orient).fCCInputBits.GetBits();
    stringstream out;
    // const unsigned int nplanes = GetRPNPlanesDividedByTwo();
    // std::vector<std::bitset<RPCCInputBits::fNPlanes> > coBits(GetNSectors());
    // boost::dynamic_bitset<>  test;
    std::vector<boost::dynamic_bitset<> > coBits(GetNSectors(), boost::dynamic_bitset<>(GetRPNPlanesDividedByTwo()));
    // boost::dynamic_bitset<>
    //std::bitset<GetRPNPlanesDividedByTwo()> > coBits;
    // cout << __func__ << endl;
    coBits.resize(GetNSectors());
    for( unsigned int i = 0 ; i < GetNSectors() ; i++ ){
        for(unsigned int j = 0; j < GetRPNPlanesDividedByTwo(); j++){
            // It is a feature of bitset that:
            //   std::bitset<3> bset;
            //   bset[0]=0;
            //   bset[1]=0;
            //   bset[2]=1;
            //
            //   // first method
            //   cout << bitset << endl;
            //  
            //   // second method
            //   for(unsigned int i = 0; i<3; i++)
            //        cout << bitset[i] << endl;
            //
            // produces:
            //   100  ...  first method
            //   001  ...  second method
            // i.e. REVERSE order!!

            // coBits[(nSectors-1)-i][j] = bs[i+j*nSectors];  // The (nSectors-1) factor is important, otherwise wrong order order...
            coBits[i][j] = bs[i+j*GetNSectors()];
        }
    }

#if 0
    for(unsigned int j=0; j<GetRPNPlanesDividedByTwo(); j++){
        for( unsigned int i = 0; i < GetNSectors() ; i++ ){
            out <<  coBits[i][j];
        }
        out << endl;
    }

    for( unsigned int i = 0 ; i < nSectors ; i++ ){
        out << coBits[i].count();
    }

    out << "  CHECK" << endl;
#endif

    for( unsigned int i = 0 ; i < GetNSectors() ; i++ ){
        for(unsigned int j=0; j < GetRPNPlanesDividedByTwo(); j++){

            if(hname.Contains("Raw0Simu1_FailingSectors")){
                // if(coBits[i][j]  and (GetPot(kRaw).GetHalfPot(orient).fCCOutputBits.getBS()[i] != GetPot(kSimu).GetHalfPot(orient).fCCOutputBits.getBS()[i])){
                if(coBits[i][j]  and (GetPot(kRaw).GetHalfPot(orient).fCCOutputBits.getBS()[i] == 0) and ( GetPot(kSimu).GetHalfPot(orient).fCCOutputBits.getBS()[i] ==1)){
                    hist->Fill(i,j+1);

                    // std::cout << " ########################################################################"  << endl;
                    // std::cout << " ###  " << "  Event " << fNEvents  << endl;

                    // PrintEventTwoHalfPot(event);
                }
            }

            if(hname.Contains("NotTheSame") and !hname.Contains("Not2planes")){
                if (coBits[i].count()==2 and coBits[i][j]  and(GetPot(kRaw).GetHalfPot(orient).fCCOutputBits.getBS()[i] != GetPot(kSimu).GetHalfPot(orient).fCCOutputBits.getBS()[i])){
                    hist->Fill(i,j+1);
                }
            }else if(hname.Contains("NotTheSame") and hname.Contains("Not2planes")){
                if (coBits[i].count()!=2 and coBits[i][j] and(GetPot(kRaw).GetHalfPot(orient).fCCOutputBits.getBS()[i] != GetPot(kSimu).GetHalfPot(orient).fCCOutputBits.getBS()[i])){
                    hist->Fill(i,j+1);
                }
            }else if(!hname.Contains("NotTheSame") and hname.Contains("TheSame") and !hname.Contains("Not2planes")){
                if(coBits[i][j] and coBits[i].count()==2) hist->Fill(i,j+1);

            } else if(!hname.Contains("NotTheSame") and hname.Contains("TheSame") and hname.Contains("Not2planes")){
                if(coBits[i][j] and coBits[i].count()!=2) hist->Fill(i,j+1);
            }
        }
    }

    return out.str();
}


void RawVsSimuPotComparator::FailingSectors() {
    assert(0);
#if 0    
    if(GetPot(kRaw).GetEvenCCOutputBits().getBS().count()==0 or GetPot(kRaw).GetOddCCOutputBits().getBS().count()==0){
        assert(0);
        for(unsigned int i=0; i < GetPot(kRaw).GetEvenCCOutputBits().getBS().size(); i++){

            if((GetPot(kRaw).GetEvenCCOutputBits().getBS()[i] != GetPot(kSimu).GetEvenCCOutputBits().getBS()[i])){
                FindHist1D("NotTheSameRaw0FailingSectorsEven")->Fill(i);
            }

            if(GetPot(kRaw).GetOddCCOutputBits().getBS()[i] != GetPot(kSimu).GetOddCCOutputBits().getBS()[i])
                FindHist1D("NotTheSameRaw0FailingSectorsOdd")->Fill(i);          

        }

    }
#endif    
}

void RawVsSimuPotComparator::StripsOnDistrib(const TString& conditionName, const TString&, const edm::Event& event, EDetOrientation orient,unsigned int iplane ) {
    TString distribName = "_StripsON_";
    TString hname = TString(GetName())+"_"+conditionName+distribName+OrientedRPId::DetOrientation(orient)+"_Plane";
    hname+=iplane;
    TH1D* hist = FindHist1D(hname);
    if(!hist){
        fHist1D.push_back(new TH1D(hname, "" , GetNStrips()  , -0.5 , GetNStrips()-0.5));
        hist = fHist1D.back();
    }

    // cout << hname << endl;

    edm::Handle< edm::DetSetVector<RPStripDigi> >  stripDigiCollection;
    event.getByLabel(stripDigiLabel, stripDigiCollection );

    for(edm::DetSetVector<RPStripDigi>::const_iterator stripDigiIt=stripDigiCollection->begin(); stripDigiIt != stripDigiCollection->end(); ++stripDigiIt) {
        for(edm::DetSet<RPStripDigi>::const_iterator it=stripDigiIt->begin(); it != stripDigiIt->end(); ++it) {
            TotRPDetId rpId(it->GetDetId());
            if(GetDecRPIdFull()!= (unsigned int) rpId.RPCopyNumber()) continue;
            unsigned int plane = rpId.Detector()/2+1;
            if(plane == iplane){// TotRPDetId::RawToDecId(it->GetDetId())
                if((orient == kOddU and rpId.IsStripsCoordinateUDirection()) or (orient == kEvenV and rpId.IsStripsCoordinateVDirection()) ){
                    // std::cout << kOddU << "=orient plane=" << plane << "  detId=" << rpId.Detector() << "  strip=" <<  it->GetStripNo() << "  sector=" << it->GetStripNo()/32 << endl;
                    // assert(0);
                    //if(it->GetStripNo()!= 103 and it->GetStripNo()!= 104)

                   // if(hname.Contains(distribName)){
                       // if(verbosity > 1) cout << "plane=" << plane << "  " << it->GetStripNo() << endl;
                        // hist->Fill(it->GetStripNo(),plane);
                        hist->Fill(it->GetStripNo());
                   //  }
                   // if(name.Contains("SectorsON"))
                    //    hist.Fill(it->GetStripNo()/(GetNStrips()/GetNSectors()),plane);
                }
            }
        }
    }
}

void RawVsSimuPotComparator::SectorsOnDistrib(const TString& conditionName, const TString&, const edm::Event& event, EDetOrientation orient,unsigned int iplane ) {
    TString distribName = "_SectorsON_";
    TString hname = TString(GetName())+"_"+conditionName+distribName+OrientedRPId::DetOrientation(orient)+"_Plane";
    hname+=iplane;
    TH1D* hist = FindHist1D(hname);
    if(!hist){
        fHist1D.push_back(new TH1D(hname, "" , GetNSectors()  , -0.5 , GetNSectors()-0.5));
        hist = fHist1D.back();
    }

    edm::Handle<edm::DetSetVector<RPDetTrigger> > inputTrigger;
    event.getByLabel(detTriggerLabel, inputTrigger );

    for(edm::DetSetVector<RPDetTrigger>::const_iterator rpDetTrigIt=inputTrigger->begin(); rpDetTrigIt != inputTrigger->end(); ++rpDetTrigIt) {
        for(edm::DetSet<RPDetTrigger>::const_iterator it=rpDetTrigIt->begin(); it != rpDetTrigIt->end(); ++it) {
            TotRPDetId rpId(it->GetDetId());
            if(GetDecRPIdFull()!= (unsigned int) rpId.RPCopyNumber()) continue;
            unsigned int plane = rpId.Detector()/2+1;
            if(plane == iplane){
                if((orient == kOddU and rpId.IsStripsCoordinateUDirection()) or (orient == kEvenV and rpId.IsStripsCoordinateVDirection()) ){
                    // std::cout << kOddU << "=orient plane=" << plane << "  detId=" << rpId.Detector() << "  strip=" <<  it->GetStripNo() << "  sector=" << it->GetStripNo()/32 << endl;
                    // assert(0);
                    //if(it->GetStripNo()!= 103 and it->GetStripNo()!= 104)

                    // if(hname.Contains(distribName)){
                     //   assert(0);
                     //   cout << " filling " << endl;
                        hist->Fill(it->GetSector());
                   // }
                }
            }
        }
    }
}

unsigned int RawVsSimuPotComparator::GetRPNPlanesDividedByTwo() const{ return GetPot(kRaw).GetHalfPot(kEvenV).GetNPlanes(); }
unsigned int RawVsSimuPotComparator::GetRPPlanes() const{ return 2*GetPot(kRaw).GetHalfPot(kEvenV).GetNPlanes();}
unsigned int RawVsSimuPotComparator::GetNSectors() const{ return GetPot(kRaw).GetHalfPot(kEvenV).GetNSectors();}
unsigned int RawVsSimuPotComparator::GetNStrips() const{ return GetPot(kRaw).GetHalfPot(kEvenV).GetNStrips();}
