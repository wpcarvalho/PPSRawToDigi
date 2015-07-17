#ifndef _L1TriggerTotemPotCollection_H_
#define _L1TriggerTotemPotCollection_H_

//#include <boost/shared_ptr.hpp>
//#include <boost/ptr_container/ptr_vector.hpp>

// system
#include <vector>
// #include <boost/ptr_container/ptr_vector.hpp>
#include "L1TriggerTotem/CoincidenceChip/interface/RawVsSimuPotComparator.h"


class PotCollection: public TObject, public std::vector<RawVsSimuPotComparator*>{
// class PotCollection: public boost::ptr_vector<RawVsSimuPotComparator>{
  public:
    // PotCollection(): boost::ptr_vector<RawVsSimuPotComparator>(){}
    // PotCollection(): std::vector<RawVsSimuPotComparator*>(){}
    PotCollection(){}
    virtual ~PotCollection(){ /* cout << __func__ << endl;*/ }
    //virtual std::string ClassName() const { return "PotCollection"; }
    const RawVsSimuPotComparator* FindPot(unsigned int rpID) const{
         for(unsigned int i=0; i < this->size(); ++i){
             // std::cout << __func__ << " "  << (*this)[i]->DetectorDecId() << "NULL" << endl;
             if((*this)[i]->GetDecRPIdFull() == rpID){
                 // std::cout << __func__ << "  NOTNULL" << endl;
                 // return &(*this)[i] ;
                 return (*this)[i] ;
             }
         }
       
         return 0;
     }

     RawVsSimuPotComparator* FindPot(unsigned int rpID){
         for(unsigned int i=0; i < this->size(); ++i){
             // std::cout << __func__ << " "  << (*this)[i]->DetectorDecId() << "NULL" << endl;
             if((*this)[i]->GetDecRPIdFull() == rpID){
                 // std::cout << __func__ << "  NOTNULL" << endl;
                 // return &(*this)[i] ;
                 return (*this)[i] ;
             }
         }
       
         return 0;
     }

#if 0    
     HalfPot* FindHalfPot(DAQInformationRP::RPID rpID, EDetOrientation orient, ERawOrSimu type){
         RawVsSimuPotComparator* potInfo = FindPot(rpID);

         if(potInfo)
             return &potInfo->GetHalfPot(orient);
         else
             return 0;

     }
     const HalfPot* FindHalfPot(DAQInformationRP::RPID rpID, EDetOrientation orient, ERawOrSimu type) const{
         const RawVsSimuPotComparator* potInfo = FindPot(rpID);

         if(potInfo)
             return &potInfo->GetHalfPot(orient);
         else
             return 0;
     }
#endif    

#if 0    
     void AddCCInputBits(RPCCInputBits* input){
         HalfPot* halfPot = FindHalfPot(input->RomanPot(),input->GetDetOrientation());
         if(halfPot) halfPot->fCCInputBits.SetBits(input->GetBits());
     }

     void AddCCInputBits(DAQInformationRP::RPID rpID, EDetOrientation orient, CoincidenceChip::InputBits bitSet){
        OrientedRPId totDet(rpID, orient);
        AddCCInputBits(new RPCCInputBits(totDet,bitSet));
     }

#if 0
     void Add(uint32_t id, const CoincidenceChip::InputBits& iCCInputBitsRPCC){
         this->push_back(new RPCCInputBits(id,iCCInputBitsRPCC));
     }
#endif
     
     void AddCCInputBits(const TotRPDetId& detId, const  RPCCInputBits::RPDetBits& rpDetBits){
           // cout << __func__ << " asdf=" << detId.DetectorDecId() << endl;
         OrientedRPId orientedRPId(detId);

         PotInfo* potInfo =  FindPot(orientedRPId.RomanPot());
         if(!potInfo) return;

          // RPCCInputBits* inputBits =
          potInfo->GetHalfPot(orientedRPId.GetDetOrientation()).fCCInputBits.SetBits(detId,rpDetBits);
          // if(inputBits){
          //  inputBits->SetBits(detId,rpDetBits);
          // }else{
          //    AddCCInputBits(new RPCCInputBits(detId,rpDetBits));
          // }
     }

     RPCCInputBits*  FindCCInputBits(const TotRPDetId& id){
         OrientedRPId orientedRPId(id);

         //  std::cout  << __func__ << " asdf=" <<  id.DetectorDecId() << endl;
         // unsigned int orientation = id.IsStripsCoordinateUDirection();
         //  unsigned int rpId = id.RomanPot();
         for(unsigned int i=0; i < this->size(); ++i){
             // std::cout << __func__ << " "  << (*this)[i]->DetectorDecId() << "NULL" << endl;

             PotInfo *potInfo = FindPot(orientedRPId.RomanPot());
             if(potInfo)  return &potInfo->GetHalfPot(orientedRPId.GetDetOrientation()).fCCInputBits;

            // if((id.IsStripsCoordinateUDirection()==(*this)[i].IsStripsCoordinateUDirection()) and(id.RomanPot() == (*this)[i].RomanPot()) ){
                 // std::cout << __func__ << "  NOTNULL" << endl;
           //      return &(*this)[i] ;
            // }
         }
         // std::cout << __func__ << "  NULL" << endl;

         //edm::LogInfo("CClogic") <<  "inputBits" << nevents;
         //cout << "Can not find id \"" <<  id.DetectorDecId() << "\"" << endl;
         //throw cms::Exception("PPBckgAnalyzer2::FindHist2") << "Can not find histogram named " << name << endl;

         //assert(0);
         return 0;
     }

     RPCCInputBits*  FindCCInputBits(unsigned int decRPId, EDetOrientation orientation){
         OrientedRPId id(decRPId,orientation);
         return FindCCInputBits(id);
     }

#endif

     //label as string because InputTag requires c++11.
    void CreateNewRPCCInput(const edm::Event& event, const std::string detTriggerLabel);
    ClassDef(PotCollection,1)
};


#endif
