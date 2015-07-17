#include "L1TriggerTotem/CoincidenceChip/interface/PotCollection.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStats.h"

// CMSSW 
#include "DataFormats/Common/interface/DetSetVector.h"

// TOTEM
#include "DataFormats/TotemRPDataTypes/interface/RPDetTrigger.h"

// ClassImp(PotCollection)
 
void PotCollection::CreateNewRPCCInput(const edm::Event& event, const std::string detTriggerLabel) {
    //cout << __func__ << endl;
    // reset all CC input bits first;
    for(unsigned int i=0; i< this->size(); i++){
        (*this)[i]->GetPot(kRaw).GetHalfPot(kOddU).fCCInputBits.Reset();
        (*this)[i]->GetPot(kRaw).GetHalfPot(kEvenV).fCCInputBits.Reset();
        (*this)[i]->GetPot(kSimu).GetHalfPot(kOddU).fCCInputBits.Reset();
        (*this)[i]->GetPot(kSimu).GetHalfPot(kEvenV).fCCInputBits.Reset();
    }

   // cout <<  __func__ << "  sonda1" << endl;
    // Read real strips: input of CC
    edm::Handle<edm::DetSetVector<RPDetTrigger> > detTrig;
    event.getByLabel(detTriggerLabel, detTrig );

    for(edm::DetSetVector<RPDetTrigger>::const_iterator detTrigIt = detTrig->begin(); detTrigIt != detTrig->end(); detTrigIt++) {

      RPCCInputBits::RPDetBits rpDetBits(CoincidenceChip::OutputBits_size);
      rpDetBits.reset();
      for (unsigned int i = 0; i < (detTrigIt->data).size(); ++i) {

    //cout <<  __func__ << "  sonda1a" << endl;
          //cout << "id= " << TotRPDetId::RawToDecId((detTrigIt->data)[i].GetDetId()) << "  sec=" << (detTrigIt->data)[i].GetSector() << endl;
   // cout <<  __func__ << "  sonda1b" << endl;
          unsigned int triggeredSectorNo = (detTrigIt->data)[i].GetSector();
          rpDetBits.set(triggeredSectorNo);
      }

   // cout <<  __func__ << "  sonda2" << endl;
      TotRPDetId detId(detTrigIt->id);
      OrientedRPId orientedRPId(detId);
      //cout << TotRPDetId::RawToDecId(detectorId) << endl;
      //if(verbosity) std::cout << rpDetBits << "  " <<  __func__ << " id=" << detectorId.DetectorDecId() << endl;
      //if(0) cout << "U?=" << detectorId.IsStripsCoordinateUDirection() << endl;
 //     cout << "rpID =" << orientedRPId.RomanPot() << "   ??" << detId.RomanPot() << "  !!" << detId.DetectorDecId()/10 << endl;
 //     cout << "rpID2=" << orientedRPId.RomanPot() << "   ??" << detId.RomanPot() << "  !!" << detId.DetectorDecId()/10 << endl;
      RawVsSimuPotComparator* potInfo =  FindPot(orientedRPId.DetectorDecId()/10);
      if(!potInfo) continue;
     // cout << " asd=" << endl;
      potInfo->GetPot(kRaw).GetHalfPot(orientedRPId.GetDetOrientation()).fCCInputBits.SetBits(detId,rpDetBits);
      potInfo->GetPot(kSimu).GetHalfPot(orientedRPId.GetDetOrientation()).fCCInputBits.SetBits(detId,rpDetBits);
    }
    // cout <<  __func__ << "  sonda3" << endl;

#if 0
    CoincidenceChip::InputBits nullBitSet; nullBitSet.reset();
  for(unsigned int i=0; i< fPotCollection.size(); i++){
      // Check if event is empty (input to CC consists of zeros)
      // if we have no input for given RP then create null bitset

      if(!FindHalfPot(fPotCollection[i].id,kOddU))
          ccInputBitsColl->Add(fPotCollection[i].id, kOddU, nullBitSet);

      if(!ccInputBitsColl->Find(fPotCollection[i].id,kEvenV))
          ccInputBitsColl->Add(fPotCollection[i].id, kEvenV, nullBitSet);

  }
#endif

  // return ccInputBitsColl;
}
