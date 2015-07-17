// -*- C++ -*-
//
// Package:    L1TriggerTotem
// Class:      RPTriggerAnalyzer
//
// Original Author:  Jiri Prochazka
//         Created:  Mon Mar  1 15:34:46 CET 2010
// $Id$

#ifndef _L1TriggerTotemTriggerStats_H_
#define _L1TriggerTotemTriggerStats_H_

// system 
#include <boost/shared_ptr.hpp>
#include "boost/dynamic_bitset.hpp"
// #include <boost/ptr_container/ptr_vector.hpp>

// ROOT
#include "TString.h"

// CMSSW 
#ifndef __CINT__
#include "FWCore/Framework/interface/Event.h"
#endif

// TOTEM
#include "TotemCondFormats/DAQInformation/interface/DAQInformationRP.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQRecord.h"
#ifndef __CINT__
#include "DataFormats/TotemL1Trigger/interface/RPCCBits.h"
#include "L1TriggerTotem/CoincidenceChip/interface/CoincidenceChip.h"
#endif
#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"

#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStat.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStatCollection.h"


// enum EDetOrientation {kEvenV = 0, kOddU = 1};
// enum ERawOrSimu {kRaw, kSimu};

#if 1

class OrientedRPId: public  TotRPDetId{
    public:
        OrientedRPId(): TotRPDetId(){}
        virtual ~OrientedRPId(){}

        /// Construct from a raw id. It is required that the Detector part of
        /// id is Totem and the SubDet part is RP, otherwise an exception is thrown.
        explicit OrientedRPId(uint32_t id): TotRPDetId(id){};

        /// Construct from fully qualified identifier.
        OrientedRPId(unsigned int Arm, unsigned int Station, unsigned int RomanPot, EDetOrientation iDetOrient):
            TotRPDetId(Arm,Station,RomanPot,iDetOrient == kEvenV?0:1){} ;

        OrientedRPId(unsigned int rpId, EDetOrientation iDetOrient):
            TotRPDetId(rpId/100,(rpId-(rpId/100)*100)/10,rpId - (rpId/10)*10, iDetOrient == kEvenV?0:1){ 
                // cout << rpId/100  << endl;
                // cout << (rpId-(rpId/100)*100)/10 << endl;
                // cout << rpId - (rpId/10)*10 << endl;
                // cout << rpId << "  vs. "  << DetectorDecId() << endl;
                assert(rpId*10+(iDetOrient== kEvenV?0:1) == DetectorDecId());
            } ;

        unsigned int GetDecRPIdFull() const { return  DetectorDecId()/10;}

        OrientedRPId(const OrientedRPId& id): TotRPDetId(id.Arm(),id.Station(),id.RomanPot(),id.IsStripsCoordinateUDirection()?kOddU:kEvenV){}
        OrientedRPId(const TotRPDetId& id): TotRPDetId(id.Arm(),id.Station(),id.RomanPot(),id.Detector()){ }

        virtual std::string ClassName() const { return "OrientedRPId"; }


        EDetOrientation GetDetOrientation() const{
            if(IsStripsCoordinateUDirection())
                return kOddU;
            else
                return kEvenV;       
        }

        static EDetOrientation SecondDetOrientation(EDetOrientation orient){
            if(orient==kEvenV)
                return kOddU;
            else
                return kEvenV;   
        }

        static TString DetOrientation(EDetOrientation orient){
            TString name;
            switch (orient) {
                case kEvenV:  name = "even"; break;
                case kOddU:   name = "odd"; break;
                default: assert(0);
            }
            return name;
        }
};

class RPCCInputBits : public OrientedRPId{
    private:
#if 0        
        static const unsigned int fNPlanes=5;   // 5;  
        static const unsigned int fNSectors=16; // 16;
        static const unsigned int fNPlanes=8;   // 5;  
        static const unsigned int fNSectors=10; // 16;
#endif

    private:
        void Init(){
           // fNPlanes=5;
           // fNSectors=16;
        }

    public:
        // typedef std::bitset<fNSectors> RPDetBits;
        typedef boost::dynamic_bitset<> RPDetBits;
#if 0        
        unsigned int GetNPlanes() const{ return 8; }
        unsigned int GetNSectors() const{ return 10; }
#else
        unsigned int GetNPlanes() const{ return  5; }
        unsigned int GetNSectors() const{ return 16; }
#endif

        RPCCInputBits(): OrientedRPId(){
            fCCInputBits.reset();
            Init();
        }

        RPCCInputBits(DAQInformationRP::RPID rpID, EDetOrientation orient): OrientedRPId(rpID,orient){
            fCCInputBits.reset();
            Init();
        }

        ~RPCCInputBits(){ /*cout << __func__ << endl; */}
        virtual std::string ClassName() const { return "RPCCInputBits"; }

        /// Construct from a raw id. It is required that the Detector part of
        /// id is Totem and the SubDet part is RP, otherwise an exception is thrown.
        explicit RPCCInputBits(uint32_t id, const CoincidenceChip::InputBits& iCCInputBits): OrientedRPId(id), fCCInputBits(iCCInputBits){Init(); /* fCCInputBits.reset();*/};
        explicit RPCCInputBits(uint32_t id): OrientedRPId(id){
           Init(); fCCInputBits.reset();
        };

        /// Construct from fully qualified identifier.
        RPCCInputBits(unsigned int Arm, unsigned int Station, unsigned int RomanPot, EDetOrientation iDetOrient, const CoincidenceChip::InputBits& iCCInputBits):
            OrientedRPId(Arm,Station,RomanPot,iDetOrient), fCCInputBits(iCCInputBits){ Init();} ;
        RPCCInputBits(const TotRPDetId& detId, const CoincidenceChip::InputBits& iCCInputBits): OrientedRPId(detId), fCCInputBits(iCCInputBits){  Init();/* fCCInputBits.reset();*/ }
        RPCCInputBits(const TotRPDetId& detId, const RPCCInputBits::RPDetBits& rpDetBits): OrientedRPId(detId){
            Init();
            fCCInputBits.reset();
            SetBits(detId,rpDetBits);
        }

        void Reset() {
            fCCInputBits.reset();
        }

        CoincidenceChip::InputBits& GetBits() { return fCCInputBits;}
        const CoincidenceChip::InputBits& GetBits() const { return fCCInputBits;}
        void SetBits(const CoincidenceChip::InputBits& iCCInputBits){ fCCInputBits = iCCInputBits;}
        //void SetBits(const RPDetTrigger& detTrig ){ 
        void SetBits(const TotRPDetId&  detId, const  RPCCInputBits::RPDetBits& rpDetBits){ 
            assert(detId.RomanPot() == this->RomanPot());
            assert(rpDetBits.size()==GetNSectors());
            assert(detId.IsStripsCoordinateUDirection() == this->IsStripsCoordinateUDirection());
            unsigned int detPlane = detId.Detector();
            unsigned int index = 0;
            if(detPlane%2==0)
                index = detPlane/2;
            else
                index = (detPlane-1)/2;

            for(unsigned int i=0; i < rpDetBits.size(); i++){
               // std::cout << __func__ << "  index=" << index << " " << fCCInputBits[i+index*rpDetBits.size()] << "  " <<  rpDetBits[i] << endl;
              fCCInputBits[i+index*rpDetBits.size()] = rpDetBits[i]; 
              
            }
        }

    private:
        CoincidenceChip::InputBits fCCInputBits;
};
#endif

class HalfPot: public OrientedRPId{
    public:
        HalfPot(DAQInformationRP::RPID rpId, EDetOrientation orient, ERawOrSimu type, const edm::ParameterSet& conf): 
            OrientedRPId(rpId,orient), 
            fCCInputBits(rpId,orient), 
            fRawOrSimu(type),
            fConfig(conf)
    {
    }
        ~HalfPot(){}

        void FindAndSetCCOutputBits(const edm::Event& event, const edm::EventSetup& iSetup);

        void analyze(const edm::Event& event, const edm::EventSetup& iSetup){
            FindAndSetCCOutputBits(event,iSetup);
        }

        RPCCInputBits fCCInputBits;
        RPCCBits fCCOutputBits;
        ERawOrSimu fRawOrSimu;

        //   boost::ptr_vector<TH1D> fHist1D;
        //   boost::ptr_vector<TH2D> fHist2D;

        // unsigned int GetNPlanes() const{ return 5; }
        // unsigned int GetNSectors() const{ return 16; }
        unsigned int GetNPlanes() const{ return fCCInputBits.GetNPlanes(); }
        unsigned int GetNSectors() const{ return fCCInputBits.GetNSectors(); }
        unsigned int GetNStrips() const{ return 512; }

        TString PrintCCInputBits() const{
            const CoincidenceChip::InputBits& bs = fCCInputBits.GetBits(); 
            // cout << __func__ << " X  sonda 1" << endl;

            // std::vector<std::bitset<fNSectors> > s(GetNPlanes());
            std::vector<boost::dynamic_bitset<> > s(GetNPlanes(),boost::dynamic_bitset<>(GetNSectors()));
             // cout << __func__ << " X  sonda 2" << endl;
             // cout <<"size=" << s.size() << endl;

            for( unsigned int i = 0 ; i < GetNSectors() ; i++ ){
                for(unsigned int j=0; j<GetNPlanes(); j++){
                    s[j][i] = bs[i+j*GetNSectors()];
                }
            }

            stringstream out;
            for(unsigned int j=0; j<GetNPlanes(); j++){
                out << s[j] << endl;
            }

            return out.str();
            //if(0) cout << __func__ << " X  sonda 2" << endl;
        }

    private: 
        const edm::ParameterSet& fConfig;
};

inline void HalfPot::FindAndSetCCOutputBits(const edm::Event& event, const edm::EventSetup& iSetup) {
    // Mask eventualy CC output bits in the case of simulation 
    if(fRawOrSimu==kSimu){
        const edm::ParameterSet& params = fConfig.getParameterSet("params");
        std::vector<std::string> maskCC = params.getParameter<std::vector<std::string> >("maskCC");
        for(unsigned int i=0; i < maskCC.size(); i++){
            if(!strcmp(maskCC[i].c_str(),TotRPDetId::RPName(GetDecRPIdFull()).c_str())){
                // cout << "maskCC[i]=" << maskCC[i] << "  " << fCCOutputBits.getBS() <<  endl; 
                 fCCOutputBits.reset();
                 return;
                // cout << "maskCC[i]=" << maskCC[i] << "  " << fCCOutputBits.getBS() <<  endl; 
            }
        }
    }
    const edm::InputTag tagRaw(fConfig.getParameter<std::string>("modulLabelRaw"),fConfig.getParameter<std::string>("productLabelRaw"));
    const edm::InputTag tagSimu(fConfig.getParameter<std::string>("modulLabelSimu"),fConfig.getParameter<std::string>("productLabelSimu"));
    const edm::InputTag* tag=0;
    if(fRawOrSimu==kRaw)
        tag=&tagRaw;
    else
        tag=&tagSimu;

    // for(unsigned int orient = 0; orient < 2; orient++)
    // GetHalfPot(orient).analyze(event,iSetup, iTag);
  //  cout << __func__ << " 2halfpot   type=" << fRawOrSimu << "  instance=" << tag->instance() << "  label=" << tag->label() <<endl;

	//if(0) cout << __func__ << " X  sonda 1" << endl;
  //   if(verbosity) std::cout << event.id()  << std::endl;

	// Don't forget to reset CC bitset before analyzing event !!!
	fCCOutputBits.reset();

	// edm::InputTag iTag(modulLabel.Data(), productLabel.Data());
	// Read output of CC 
	edm::Handle<std::vector<RPCCBits> > outputCC;
	//event.getByLabel(modulLabel.Data(), productLabel.Data(), outputCC); 
	event.getByLabel(*tag, outputCC); 

	// CoincidenceChip::OutputBits nullBitset;
	// nullBitset.reset();
	unsigned int arm       = this->Arm(); // fRPID/100; // 1;
	unsigned int rpStation = this->Station(); // (fRPID - arm*100)/10; // 2;
	unsigned int rp        = this->RomanPot(); //fRPID - arm*100 - rpStation*10; // 0;
    // cout << "Arm=" << arm << "  rpStation=" << rpStation << "  rp=" << rp << endl;
	RPCCBits ccBits;
	ccBits.setId(TotRPDetId(arm,rpStation,rp,GetDetOrientation()));
	ccBits.reset();
	fCCOutputBits.setId(ccBits.getId());
	
	for(std::vector<RPCCBits>::const_iterator it = outputCC->begin(); it != outputCC->end(); it++) {
		//if(0) cout << __func__ << " X  sonda 2" << endl;
		TotRPDetId detId(it->getId());
        //if(detId.DetectorDecId()/10 != fRPID) continue;
		//cout << " id=" << detId.DetectorDecId()/10 << "  rp2=" << fRPID << "   ccbits = " << it->getBS() << endl;
		//if(0) cout << __func__ << " X  sonda 3" << endl;
		ccBits.setBS(it->getBS());
        // cout << detId.DetectorDecId() << "  ==??  "  << TotRPDetId::RawToDecId(ccBits.getId()) << endl;
		if(detId.DetectorDecId() == TotRPDetId::RawToDecId(ccBits.getId())){
			fCCOutputBits.setBS(ccBits.getBS());
		}
		//else if( detId.DetectorDecId() == TotRPDetId::RawToDecId(ccBits.getId()) and (detId.IsStripsCoordinateVDirection() == kEvenV))
	}


    // cout << "sonda 1" << endl; 

    //cout <<  __func__ << "  " << fCCOutputBits.getBS() << endl; // "   type=" << fRawOrSimu << "  instance=" << iTag.instance() << "  label=" << iTag.label() <<endl;

	//if(0) cout << __func__ << " X  sonda 4" << endl;
	//if(verbosity > 1) PrintCCInputAndOutput(event,orient);
	//if(0) cout << __func__ << " X  sonda 5" << endl;
}

#endif

