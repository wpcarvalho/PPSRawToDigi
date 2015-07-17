/*
 * T1DeadChannelDataAnalyzer.cc
 *
 *  Created on: Sep 2, 2011
 *      Author: Marcin Borratynski (mborratynski@gmail.com)
 */

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"

#include "DataFormats/Common/interface/DetSetVector.h"

#include "TotemRawData/RawToDigi/interface/T1DeadChannelDataAnalyzer.h"

#include <iostream>
#include <fstream>

#define DIS_CONN 56

using namespace std;
using namespace edm;
using namespace Totem;


//----------------------------------------------------------------------------------------------------
T1DeadChannelDataAnalyzer::T1DeadChannelDataAnalyzer(edm::ParameterSet const& conf) :
		  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
		  deadCutValue(conf.getUntrackedParameter<unsigned int>("deadCutValue", 100)),
		  noiseCutValue(conf.getUntrackedParameter<unsigned int>("noiseCutValue", 9)),
		  outputFile(conf.getParameter<string>("outputFile")),
		  OnlyStrips(conf.getParameter<bool>("onlyStrips"))
{
	rawEventLabel = conf.getParameter<InputTag>("RawEventLabel");
	myFile = new ofstream (outputFile.c_str());
	if (myFile->is_open())
		cout << "Opened file" << endl;
	else
		cout << "Unable to open file" << endl;
	if(OnlyStrips)cout << "T1DeadChannelDataAnalyzer: looking ONLY for dead/noisy cathodes (not wires)"<<endl;

}
//----------------------------------------------------------------------------------------------------
T1DeadChannelDataAnalyzer::~T1DeadChannelDataAnalyzer()
{
}
//----------------------------------------------------------------------------------------------------
int numberOfChannels(int plane, int chamber, int channelType){
	// plane:		<0,4>
	// chamber:		<0,5>
	// channeltype:	<1,3>
	const int nWires[5][2] = {{165, 127}, {181, 150}, {196, 158}, {213, 180}, {226, 204}};
	const int nStrips[5][2] = {{118, 97}, {129, 108}, {138, 114}, {151, 127}, {164, 151}};
	if(channelType == 1 || channelType == 2){
		return nStrips[plane][(chamber%3)/2];
	}else if (channelType == 3){
		return nWires[plane][(chamber%3)/2] ;
	}else{
		return 0;
	}
}
//----------------------------------------------------------------------------------------------------

void T1DeadChannelDataAnalyzer::beginRun(edm::Run const& run,edm::EventSetup const& es)
{
  cout << "INSIDE BEGINRUN "<<endl;
t1DataBeginJob(es);
}
//void T1DeadChannelDataAnalyzer::beginJob()
//{
// cout << "INSIDE BEGINJOB 2"<<endl;
//}
void T1DeadChannelDataAnalyzer::t1DataBeginJob(const edm::EventSetup &es) {
  // for t1 - generating connections
  // get mapping
  ESHandle<DAQMapping> mapping;
  es.get<TotemDAQMappingRecord>().get(mapping);

  map<FramePosition, VFATInfo>::const_iterator iter;
  for(iter = mapping->VFATMapping.begin(); iter != mapping->VFATMapping.end(); iter++)
    {     
//      cout << iter->second.symbolicID.subSystem << endl;
      //leaving all non-t1 vfats
      if (iter->second.symbolicID.subSystem!=TotemSymbID::T1 || iter->second.type!=VFATInfo::data) continue;

#ifdef DEBUG
      std::cout << "Position ID:"; 
      std::cout <<iter->first<<" - ";
      std::cout << iter->first.GetRawPosition();
      std::cout << " - symId " << iter->second.symbolicID.symbolicID; 
      std::cout << " - hwId " << hex << iter->second.hwID << dec; 
      std::cout << std::endl;
#endif
    }
  // get mapping
  // Configuration parameters

  const unsigned int nWires[5][2] = {{165, 127}, {181, 150}, {196, 158}, {213, 180}, {226, 204}};
  const unsigned int nStrips[5][2] = {{118, 97}, {129, 108}, {138, 114}, {151, 127}, {164, 151}};
  const unsigned int wireGroupFirst[5][16] = {{0, 6, 13, 22, 30, 39, 50, 60, 72, 84, 96, 108, 119, 131, 143, 155},
                {0, 7, 15, 25, 33, 44, 55, 67, 80, 92, 105, 118, 131, 144, 157, 170},
                {0, 9, 17, 27, 37, 48, 60, 72, 86, 100, 114, 128, 142, 156, 169, 184},
                {0, 9, 19, 29, 40, 51, 64, 78, 92, 107, 122, 137, 152, 166, 181, 197},
                {0, 12, 22, 34, 44, 57, 71, 86, 102, 118, 134, 150, 166, 182, 198, 214}};
  const int cfecSide[5] = {1, 1, -1, -1, -1};  // +1 = CFEC "outside" layer (away from IP); -1 = CFC "inside" layer

  const string connettori[DIS_CONN]={
    "-21C13",
    "-21C22",
    "-23C03",
    "-23C20",
    "-31C10",
    "-31C21",
    "-33C10",
    "-33C21",
    "-41C12",
    "-41C13",
    "-43C12",
    "-43C13",
    "-52C12",
    "-52C13",
    "-20C10",
    "-20C13",
    "-24C03",
    "-24C20",
    "-30C12",
    "-30C13",
    "-34C10",
    "-34C23",
    "-40C10",
    "-40C11",
    "-44C10",
    "-44C11",
    "-55C10",
    "-55C11",
    "+21C10",
    "+21C13",
    "+23C10",
    "+23C13",
    "+31C12",
    "+31C13",
    "+33C12",
    "+33C13",
    "+41C10",
    "+41C11",
    "+43C12",
    "+43C13",
    "+52C12",
    "+52C13",
    "+20C12",
    "+20C13",
    "+24C10",
    "+24C13",
    "+30C12",
    "+30C13",
    "+34C22",
    "+34C23",
    "+40C12",
    "+40C13",
    "+44C10",
    "+44C11",
    "+55C10",
    "+55C11"
  };

  // CFEC connectors configuration
  // Read list of disconnected CFEC connectors
  map<string, int> disconnected;  // empty connectors in CFECs

  for(int iii = 0; iii<DIS_CONN; iii++)
    {
      disconnected[connettori[iii] ] = 1;
    }
  // Construct connectors
  cfecConn* connector[2][5][6][3][4];
  for (int iArm = 0; iArm < 2; iArm++)
    for (int iLayer = 0; iLayer < 5; iLayer++)
      for (int iSext = 0; iSext < 6; iSext++)
  for (int iView = 1; iView <= 2; iView++) 
    {
      unsigned int strip = 0;
      unsigned int lastStrip = nStrips[iLayer][(iSext%3)/2] - 1;
      unsigned int nCfec = lastStrip/64 + 1;
      int cfecIncr = (iView == 1) ? 1 : -1;
      unsigned int cfec = (iView == 1) ? 0 : nCfec - 1;
      for (unsigned int iCfec = 0; iCfec < nCfec; iCfec++) 
        {
    int connIncr = (iView == 1) ? -2*cfecSide[iLayer] : 2*cfecSide[iLayer];
    unsigned int conn = (6 - 3*connIncr)/4;
    for (int iConn = 0; iConn < 2; iConn++) 
      {
        stringstream connIdStr;
        connIdStr << (iArm == 0 ? '+' : '-' ) << iLayer + 1 << iSext << 'C' << cfec << conn;
        if (disconnected.find(connIdStr.str()) == disconnected.end()) 
          {
      unsigned int nextStrip = min(strip + 32, lastStrip + 1);
      connector[iArm][iLayer][iSext][cfec][conn] = new cfecConn(cfecSide[iLayer], iView, strip, nextStrip - 1);
      strip = nextStrip;
          } 
        else 
          {
      connector[iArm][iLayer][iSext][cfec][conn] = new cfecConn(0, 0, 0, 0);
          }
        conn += connIncr;
      }
    cfec += cfecIncr;
        }
    }

  for(iter = mapping->VFATMapping.begin(); iter!=mapping->VFATMapping.end(); iter++)
    {        
      //leaving all non-t1 vfats
      if (iter->second.symbolicID.subSystem!=TotemSymbID::T1 || iter->second.type!=VFATInfo::data) continue;

      unsigned int vfatN = iter->first.GetRawPosition();
      //printf("vfat number %i\n", vfatN);       
      uint32_t id_ = 0x72000000;
      uint32_t symId_ = 0;
      symId_ = iter->second.symbolicID.symbolicID;
      symId_ <<= 15;
      id_ |= symId_;

      T1DetId ID__(id_);
      unsigned int arm = ID__.Arm();
      unsigned int layer = ID__.Plane();
      unsigned int sextant = ID__.CSC();
      unsigned int T1Key ;
      unsigned int vfat = (id_ >> 15) & 0x3;
      unsigned int physCh;

      if(ID__.Layer() == 0){
  // anode chip
  unsigned int lastWire = nWires[layer][(sextant%3)/2] - 1;
  unsigned int wireGroup = 16 + vfat;  // first VFAT channel => last wires
  for (int vfatCh = 0; vfatCh < 128; vfatCh++) {
    if (vfatCh%16 == 0)
      wireGroup -= 2;
    physCh = wireGroupFirst[layer][wireGroup] + (15 - vfatCh%16);
    if (physCh <= lastWire &&
        (wireGroup == 15 || physCh < wireGroupFirst[layer][wireGroup + 1])){
      T1DetId IIIDDD(ID__.Arm(),ID__.Plane(),ID__.CSC(),3);
      T1FisChannel FC(IIIDDD,0,physCh);
      T1Key = vfatN*1000 + vfatCh;

      T1Map[T1Key] = FC;
    }
  }
      }else{
  //cathode chip
  unsigned int cfec = vfat;
  for (int iConn = 0; iConn < 4; iConn++) 
    {
      cfecConn* conn = connector[arm][layer][sextant][cfec][iConn];
      if (conn->isConnected()) 
        {
    int view =  conn->getView();
    for (int iCh = 0; iCh < 32; iCh++) 
      {
        //      unsigned int vfatCh = iCh + 32*iConn;
        unsigned int vfatCh = 127 - (iCh + 32*iConn);  // Try reversing the channel order
        int stripNumber = conn->getStrip(iCh);
        if (stripNumber >= 0)
          {
      T1DetId IIIDDD(ID__.Arm(),ID__.Plane(),ID__.CSC(),view);
      T1FisChannel FC(IIIDDD,0,stripNumber);
      T1Key = vfatN*1000 + vfatCh;

      T1Map[T1Key] = FC;
          }
      }
        }
    }
      } 
    }
      //-----------create T1deadChannelManager
	edm::ESHandle<AnalysisMask> analysisMask;
	es.get<TotemDAQMappingRecord> ().get(analysisMask);
	_T1deadChannelManager = T1DeadChannelManager(analysisMask); //set analysisMask in deadChannelsManager
          
}

//----------------------------------------------------------------------------------------------------
void T1DeadChannelDataAnalyzer::t1DataProduce(const edm::Event& e, VFATFrameCollection::Iterator &fr, const VFATInfo &info) {

  const VFATFrame* frame = fr.Data();

  //    unsigned short ID = frame.getChipID();
  unsigned int i = fr.Position().GetGOHId()*16 + fr.Position().GetIdxInFiber();

  i=fr.Position().GetRawPosition();

  //    if(frame.getChipID()== 0xc26d)
  //  printf(" CIPPA LIPPA %x, %i", frame.getChipID(), i);
  // skip unmapped S-Link positions
    
  // get IDs
  unsigned short symId = info.symbolicID.symbolicID;

  // FILL T1Digi
  uint32_t id_ = 0x72000000;
  uint32_t symId_ = 0;
  
  // TODO: is this correct?
  symId_ = symId;
  symId_ <<= 15;

  id_ |= symId_;

  const vector<unsigned char> activeCh = frame->getActiveChannels();

		/* in this case layer is not defined :
		 * in simulation you have: arm=1bit, plane=3bits, csc=3bits, layer=2bits
		 * here you have:  arm=1bit, plane=3bits, csc=3bits, anode or cathode=1bit, position=2bits
		 */
		int channel;
		int arm;
		int plane;
		int chamber;
		int channelType;
		long long  key;
//		cout << activeCh.size() << endl;
		for (unsigned int j = 0; j < activeCh.size(); j++)
		{
			if(activeCh[j]>=128){
				cout << "Error : activeCh >= 128" << endl;
				throw (const char*)"activeCh >= 128";
			}

			int electronic_index = i*1000 + activeCh[j];

			T1FisChannel FC2 = T1Map.find(electronic_index)->second;
		  T1DetId myDetid = FC2.DetId();
		  if(myDetid.Layer()==3 || myDetid.Layer()==2 || myDetid.Layer()==1)
			if( FC2.Type()==0 || FC2.Type()==-1 ){

				arm = FC2.DetId().Arm();
				plane = FC2.DetId().Plane();
				chamber = FC2.DetId().CSC();
				channelType = FC2.DetId().Layer();
				channel = FC2.Channel()+1;
//				cout << arm << " "<< plane << " "<< chamber << " "<< channelType << " "<< channel << " "<<endl;
				key = 	1000000	* arm +
						100000	* plane +
						10000 	* chamber +
						1000	* channelType +
								  channel;
			       
				channelsHitNumber[key] += 1;
			}
		}

}
void T1DeadChannelDataAnalyzer::endJob()
{
	long long sum = 0;
	int numberOfDead = 0;

	for(int arm = 0; arm < 2; arm ++){
		for(int plane = 0; plane < 5; plane++){

			//-----------avg computation-----------------------------
			long long plane_sum = 0;
			long long plane_ctr = 0;
			long long plane_avg = 0;
			for(int chamber = 0; chamber < 6; chamber++){
				for(int channelType = 1; channelType < 4; channelType++){
					int tmpChannelNumber = numberOfChannels(plane, chamber, channelType);
					for(int channelNumber = 1; channelNumber <= tmpChannelNumber; channelNumber++){
						int key = 	1000000	* arm +
									100000	* plane +
									10000 	* chamber +
									1000	* channelType +
											  channelNumber;
//						if(key<1999)cout << channelsHitNumber[key] << endl;
						sum += channelsHitNumber[key];
						plane_sum += channelsHitNumber[key];
						plane_ctr++;
					}
				}
			}
			plane_avg = plane_sum / plane_ctr;

//			cout << plane_avg <<" "<< plane_sum << " "<< plane_ctr<<endl;;
		

			//-------------------- write out dead/noisy channels --------------------------------
			for(int chamber = 0; chamber < 6; chamber++){
				for(int channelType = 1; channelType < 4; channelType++){
					int tmpChannelNumber = numberOfChannels(plane, chamber, channelType);
					for(int channelNumber = 1; channelNumber <= tmpChannelNumber; channelNumber++){
						int key = 	1000000	* arm +
									100000	* plane +
									10000 	* chamber +
									1000	* channelType +
											  channelNumber;

//						if(arm==0 && plane==2 && chamber==4 && channelType==1)cout <<  channelsHitNumber[key]<< " " << plane_avg << " " << deadCutValue << endl;

						if(channelsHitNumber[key] >= plane_avg * noiseCutValue || channelsHitNumber[key] <= plane_avg / deadCutValue || plane_avg == 0 ){
							numberOfDead++;
							if( (OnlyStrips && channelType != 3) || !OnlyStrips) 
							(*myFile) << arm << ":" << plane << ":" << chamber  << ":" << channelType << ":" << channelNumber << endl;
							if(plane_avg == 0) cout << "NO STATISTICS!!! " <<  arm << ":" << plane << ":" << chamber  << ":" << channelType << ":" << channelNumber << endl;
							else if(channelsHitNumber[key] >= plane_avg * noiseCutValue) cout << "NOISY " <<  arm << ":" << plane << ":" << chamber  << ":" << channelType << ":" << channelNumber << endl;
							else if(channelsHitNumber[key] <= plane_avg / deadCutValue) cout << "DEAD " <<  arm << ":" << plane << ":" << chamber  << ":" << channelType << ":" << channelNumber << endl;
						

						}
					}
				}
			}
		}
	}


	//(*myFile) << "global avg: " << sum / channelsHitNumber.size() << endl;
	//map<long long, long long>::iterator it;
	//for ( it=channelsHitNumber.begin() ; it != channelsHitNumber.end(); it++ ){
	//	(*myFile) << it->first << " -> " << it->second << endl;
	//}
	//(*myFile) << "number of dead: " << numberOfDead  << endl;

	myFile->close();
	delete myFile;
}
//----------------------------------------------------------------------------------------------------
void T1DeadChannelDataAnalyzer::PrintErrorHeader(const edm::EventID &evId, const FramePosition &fp, unsigned short symId, unsigned short vfatId)
{
	if (!eventError) {
		eventError = true;
		if (verbosity > 0)
			cerr << ">> T1DeadChannelDataAnalyzer::produce > event " << evId.run() << ":" << evId.event() << " seems to contain corrupted frames" << endl;
	}

	if (!positionError) {
		positionError = true;
		if (verbosity > 1) {

			cerr << "\tcorrupted frame at position " << fp << " (symId = " << symId;
			if (vfatId != 0)
				cerr << ", vfatId = 0x" << hex << vfatId << dec;
			cerr << ")" << endl;
		}
	}
}
//----------------------------------------------------------------------------------------------------

void T1DeadChannelDataAnalyzer::analyze(const edm::Event& e, const edm::EventSetup& es)
{
// cout << "INSIDE ANALYZE "<<endl;
  stringstream ees;



  // get input
  Handle< RawEvent > input;
  e.getByLabel(rawEventLabel, input);
  if (!input.isValid())
    return;



  // get mapping
  ESHandle<DAQMapping> mapping;
  es.get<TotemDAQMappingRecord>().get(mapping);


  // process all GOOD frames
  for (VFATFrameCollection::Iterator fr(input->frames); !fr.IsEnd(); fr.Next()) {
    map<FramePosition, VFATInfo>::const_iterator mappingIter = mapping->VFATMapping.find(fr.Position());
    if (mappingIter == mapping->VFATMapping.end())
      continue;

 
    if(mappingIter->second.symbolicID.subSystem == TotemSymbID::T1)
      t1DataProduce(e, fr, mappingIter->second);
   
    // decide which method should process that frame

  }


}


DEFINE_FWK_MODULE(T1DeadChannelDataAnalyzer);
