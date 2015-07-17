/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*   Maciej Wróbel (wroblisko@gmail.com)
*   Jan Kašpar (jan.kaspar@gmail.com)
*   Marcin Borratynski (mborratynski@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Utilities/interface/EDMException.h"

#include "DataFormats/TotemRPDetId/interface/TotRPDetId.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "TotemRawData/RawToDigi/interface/CounterChecker.h"
#include "TotemRawData/RawToDigi/interface/Raw2DigiProducer.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawData/RawToDigi/interface/Raw2DigiStatus.h"
#include "TotemCondFormats/DAQInformation/interface/AnalysisMask.h"
#include "TotemCondFormats/DAQInformation/interface/DAQMapping.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQMappingRecord.h"
#include "SimTotem/T1Digitizer/interface/T1DeadChannelManager.h"
 
#include <iostream>
#include <fstream>

using namespace std;
using namespace edm;
using namespace Totem;

//#define DEBUG

//----------------------------------------------------------------------------------------------------

Raw2DigiProducer::Raw2DigiProducer(edm::ParameterSet const& conf) :
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  printErrorSummary(conf.getUntrackedParameter<unsigned int>("printErrorSummary", 1)),
  printUnknownFrameSummary(conf.getUntrackedParameter<unsigned int>("printUnknownFrameSummary", 1)),

  testFootprint(conf.getParameter<unsigned int>("testFootprint")),
  testCRC(conf.getParameter<unsigned int>("testCRC")),
  testID(conf.getParameter<unsigned int>("testID")),
  testECRaw(conf.getParameter<unsigned int>("testECRaw")),
  testECDAQ(conf.getParameter<unsigned int>("testECDAQ")),
  testECMostFrequent(conf.getParameter<unsigned int>("testECMostFrequent")),
  testBCMostFrequent(conf.getParameter<unsigned int>("testBCMostFrequent")),
  
  EC_min(conf.getUntrackedParameter<unsigned int>("EC_min", 10)),
  BC_min(conf.getUntrackedParameter<unsigned int>("BC_min", 10)),
  
  EC_fraction(conf.getUntrackedParameter<double>("EC_fraction", 0.6)),
  BC_fraction(conf.getUntrackedParameter<double>("BC_fraction", 0.6))
{
  rawEventLabel = conf.getParameter<InputTag>("RawEventLabel");
  rpCCProductLabel = conf.getUntrackedParameter<std::string>("rpCCProductLabel", "rpCCOutput");
  rpDataProductLabel = conf.getUntrackedParameter<std::string>("rpDataProductLabel", "rpDataOutput");
  t1DataProductLabel = conf.getUntrackedParameter<std::string>("t1DataProductLabel", "t1DataOutput");
  t2DataProductLabel = conf.getUntrackedParameter<std::string>("t2DataProductLabel", "t2DataOutput");
  conversionStatusLabel = conf.getUntrackedParameter<std::string>("conversionStatusLabel", "conversionStatus");

  // RP data
  produces< edm::DetSetVector<RPStripDigi> > (rpDataProductLabel);

  // RP CC
  produces < std::vector <RPCCBits> > (rpCCProductLabel);

  // T2 data
  produces<T2StripDigiCollection>(t2DataProductLabel);
  produces<T2PadDigiCollection>(t2DataProductLabel);
  produces<T2DigiVfatCollection>(t2DataProductLabel);
  produces<T2VfatInformation>(t2DataProductLabel);
  t2DiscardHighOccupancyVfat = conf.getUntrackedParameter<bool>("discardHighOccupancyVfatverbosity", false); //T2data

  // T1 data
  produces<T1DigiWireCollection>(t1DataProductLabel);
  produces<T1DigiVfatCollection>(t1DataProductLabel);       

  // status
  produces <Raw2DigiStatus>(conversionStatusLabel);
}
  
//----------------------------------------------------------------------------------------------------
 
Raw2DigiProducer::~Raw2DigiProducer()
{
}

//----------------------------------------------------------------------------------------------------
 
void Raw2DigiProducer::beginRun(const edm::Run& run,const edm::EventSetup& es)
{
  // generate connections for t1
  t1DataBeginJob(es);
}

//----------------------------------------------------------------------------------------------------
 
void Raw2DigiProducer::endJob()
{
  if (printErrorSummary) {
    cout << "* Error summary (error signature: number of such events)" << endl;
    for (map<Totem::FramePosition, map<VFATStatus, unsigned int> >::iterator vit = errorSummary.begin();
        vit != errorSummary.end(); ++vit) {
      //map<FramePosition, VFATInfo>::const_iterator mappingIter = mapping->VFATMapping.find(vit->first);
      //cout << "  " << vit->first << " (VFAT " << mappingIter->second.hwId << "):" << endl;
      cout << "  " << vit->first << endl;
      for (map<VFATStatus, unsigned int>::iterator it = vit->second.begin(); it != vit->second.end();
          ++it) {
        cout << "    " << it->first << ": " << it->second << endl;
      }
    }
  }

  if (printUnknownFrameSummary) {
    cout << "* Frames found in data, but not in the mapping (frame position: number of events)" << endl;
    for (map<FramePosition, unsigned int>::iterator it = unknownSummary.begin();
        it != unknownSummary.end(); ++it) {
      cout << "  " << it->first << ":" << it->second << endl;
    }
  }
}

//----------------------------------------------------------------------------------------------------
 
void Raw2DigiProducer::produce(edm::Event& e, const edm::EventSetup& es)
{
  rpDataOutput = auto_ptr<edm::DetSetVector<RPStripDigi> > (new edm::DetSetVector<RPStripDigi>);  
  rpCCOutput = auto_ptr< std::vector<RPCCBits> > (new std::vector<RPCCBits>);
  t2DataTheStripDigis = auto_ptr <T2StripDigiCollection> (new T2StripDigiCollection);
  t2DataThePadDigis = auto_ptr <T2PadDigiCollection> (new T2PadDigiCollection);
  t2DataTheVfats = auto_ptr<T2DigiVfatCollection> (new T2DigiVfatCollection());
  t2DatavfatsStatusMap = auto_ptr <T2VfatInformation> (new T2VfatInformation());  
  t1DigiWireColl = auto_ptr<T1DigiWireCollection> (new T1DigiWireCollection);
  t1DigiVfatColl = std::auto_ptr<T1DigiVfatCollection> (new T1DigiVfatCollection);
  
  stringstream ees;

#ifdef DEBUG
  printf(">> Raw2DigiProducer::produce\n");
  printf("\tevent #%u at %p\n", e.id().event(), (void *) &e);
#endif

  // get input
  Handle< RawEvent >  input;
  e.getByLabel(rawEventLabel, input);
  if (!input.isValid())
    return;

#ifdef DEBUG
  printf("\tevent.frames at %p\n", (void *)input->frames);
  printf("\tevent.dataEventNumber = %lu\n", input->dataEventNumber);
  printf("\tcollection type %s\n", input->frames->GetClassName().c_str());
#endif

  // get mapping
  ESHandle<DAQMapping> mapping;
  es.get<TotemDAQMappingRecord>().get(mapping);

  // get analysis mask to mask channels
  ESHandle<AnalysisMask> analysisMask;
  es.get<TotemDAQMappingRecord>().get(analysisMask);

  // map which will contain FramePositions from mapping, which data is missing in raw event
  map<FramePosition, VFATInfo> missingFrames(mapping->VFATMapping);

  // prepare products
  rpDataBeginJob();  
  t2DataBeginJob();

  // expected EC number for all frames in collection
  // CMSSW events are counted from 1, raw events from 0 => the shift by -1
  unsigned long ECExpected = (e.id().event() - 1) % 0x100; 

  // object containing the information about the conversion
  auto_ptr<Raw2DigiStatus> conversionStatus(new Raw2DigiStatus());
  
  // EC and BC checks (wrt. the most frequent value), BC checks per subsystem
  CounterChecker ECChecker(CounterChecker::ECChecker, "EC", EC_min, EC_fraction, verbosity);
  map<unsigned int, CounterChecker> BCCheckers;
  BCCheckers[TotemSymbID::RP] = CounterChecker(CounterChecker::BCChecker, "BC/RP", BC_min, BC_fraction, verbosity);
  BCCheckers[TotemSymbID::T1] = CounterChecker(CounterChecker::BCChecker, "BC/T1", BC_min, BC_fraction, verbosity);
  BCCheckers[TotemSymbID::T2] = CounterChecker(CounterChecker::BCChecker, "BC/T2", BC_min, BC_fraction, verbosity);

  // aplly tests to all frames
  for (VFATFrameCollection::Iterator fr(input->frames); !fr.IsEnd(); fr.Next()) {
    stringstream fes;
    bool stopProcessing = false;

    // contain information about processed frame
    VFATStatus &actualStatus = (*conversionStatus)[fr.Position()];
    
    // skip unlisted positions (VFATInfo)
    map<FramePosition, VFATInfo>::const_iterator mappingIter = mapping->VFATMapping.find(fr.Position());
    if (mappingIter == mapping->VFATMapping.end()) {
      unknownSummary[fr.Position()]++;
      continue;
    }

    // remove from missingFrames
    map<FramePosition, VFATInfo>::iterator iter = missingFrames.find(fr.Position());
    missingFrames.erase(iter);

    // check footprint
    if (testFootprint != tfNoTest && !fr.Data()->checkFootprint()) {
      fes << "    invalid footprint\n";
      if ((testFootprint == tfErr)) {
        actualStatus.setFootprintError();
        stopProcessing = true;
      }
    }
    
    // check CRC
    if (testCRC != tfNoTest && !fr.Data()->checkCRC()) {
      fes << "    CRC failure\n";
      if (testCRC == tfErr) {
        actualStatus.setCRCError();
        stopProcessing = true;
      }
    }

    // check the id mismatch
    if (testID != tfNoTest && (fr.Data()->getChipID() & 0xFFF) != (mappingIter->second.hwID & 0xFFF)) {
      fes << "    ID mismatch (data: 0x" << hex << fr.Data()->getChipID()
        << ", mapping: 0x" << mappingIter->second.hwID  << dec << ", symbId: " << mappingIter->second.symbolicID.symbolicID << ")\n";
      if (testID == tfErr) {
        actualStatus.setIDMismatch();
        //stopProcessing = true;
      }
    }

    // EC progress check (wrt. raw data event number, i.e. EC expected)
    if (testECRaw != tfNoTest && fr.Data()->getEC() != ECExpected) {
      fes << "    EC (" << fr.Data()->getEC() << ") doesn't match the expectation form event number (" << ECExpected << ")\n";
      if (testECRaw == tfErr) {
        actualStatus.setECProgressError();
        stopProcessing = true;
      }
    }

    // check, if EC number from VFATFrame is the same as obtained from DAQ (only if the there is not ECProgressError)
    if (testECDAQ != tfNoTest && !actualStatus.isECProgressError() && ((fr.Data()->getEC()%0x100) != (input->dataEventNumber%0x100))) {
      fes << "    EC (" << fr.Data()->getEC() << ") doesn't match DAQ event counter (" << input->dataEventNumber << ")\n";
      if (testECDAQ == tfErr) { 
        actualStatus.setECProgressError();
        stopProcessing = true;
      }
    }    

    // if there were errors, put the information to ees buffer
    if (verbosity > 0 && !fes.rdbuf()->str().empty()) {
      string message = (stopProcessing) ? "(and will be dropped)" : "(but will be used though)";
      if (verbosity > 2) {
        ees << "  Frame at " << fr.Position() << " seems corrupted " << message << ":\n";
        ees << fes.rdbuf();
      } else
        ees << "  Frame at " << fr.Position() << " seems corrupted " << message << ".\n";
    }

    // if there were serious errors, do not process this frame
    if (stopProcessing)
      continue;
    
    // fill EC and BC values to the statistics
    ECChecker.Fill(fr.Data()->getEC(), fr.Position());
    BCCheckers[mappingIter->second.symbolicID.subSystem].Fill(fr.Data()->getBC(), fr.Position());
  }

  // analyze EC and BC statistics
  if (testECMostFrequent != tfNoTest)
    ECChecker.Analyze(*conversionStatus, (testECMostFrequent == tfErr), ees);
  if (testBCMostFrequent != tfNoTest) {
    for (map<unsigned int, CounterChecker>::iterator it = BCCheckers.begin(); it != BCCheckers.end(); ++it)
      it->second.Analyze(*conversionStatus, (testBCMostFrequent == tfErr), ees);
  }

  // save the information about missing frames to conversionStatus
  for (map<FramePosition, VFATInfo>::iterator iter = missingFrames.begin(); iter != missingFrames.end(); iter++) {
    VFATStatus &actualStatus = (*conversionStatus)[iter->first];
    actualStatus.setMissing();
    if (verbosity > 1) 
      ees << "Frame for VFAT " << iter->first << " is not present in the data.\n";
  }

  // print error message
  if (verbosity > 0 && !ees.rdbuf()->str().empty()) {
    char buf[20];
    sprintf(buf, "%u.%04u", e.id().run() / 10000, e.id().run() % 10000);
    if (verbosity > 1)
      ERROR("Raw2DigiProducer") << "event " << buf << ":" << e.id().event() <<
        " contains the following problems:\n" << ees.rdbuf() << c_endl;
    else
      ERROR("Raw2DigiProducer") << "event " << buf << ":" << e.id().event() << 
        " contains problems." << c_endl;
  }

  // update error summary
  if (printErrorSummary) {
    for (Raw2DigiStatus::iterator it = conversionStatus->begin(); it != conversionStatus->end(); ++it) {
      if (!it->second.OK()) {
        map<VFATStatus, unsigned int> &m = errorSummary[it->first];
        m[it->second]++;
      }
    }
  }

  // process all GOOD frames
  for (VFATFrameCollection::Iterator fr(input->frames); !fr.IsEnd(); fr.Next()) {
    map<FramePosition, VFATInfo>::const_iterator mappingIter = mapping->VFATMapping.find(fr.Position());
    if (mappingIter == mapping->VFATMapping.end())
      continue;

    VFATStatus &actualStatus = (*conversionStatus)[fr.Position()];

    if (!actualStatus.OK())
      continue;

    // prepare analysis mask class
    map<TotemSymbID, VFATAnalysisMask>::const_iterator analysisIter = analysisMask->analysisMask.find(mappingIter->second.symbolicID);

    // find analysis mask
    VFATAnalysisMask anMa;
    anMa.fullMask = false;
    if (analysisIter != analysisMask->analysisMask.end()) {            
      // if there is some information about masked channels - save it into conversionStatus
      anMa = analysisIter->second;
      if (anMa.fullMask)
        actualStatus.setFullyMaskedOut();
      else
        actualStatus.setPartiallyMaskedOut();
    }
    
    // decide which method should process that frame
    switch (mappingIter->second.symbolicID.subSystem) {
      case TotemSymbID::RP:  
        switch (mappingIter->second.type) {
          case VFATInfo::data:
            rpDataProduce(e, fr, mappingIter->second, anMa);    
            break;
          case VFATInfo::CC:
            rpCCProduce(e, fr, mappingIter->second, anMa);
            break;
        }  
        break;

      case TotemSymbID::T1:
        t1DataProduce(e, fr, mappingIter->second);
        break;

      case TotemSymbID::T2:
        if (mappingIter->second.type==VFATInfo::data)
          t2DataProduce(e, fr, mappingIter->second, anMa);
        break;
    }
  }

  // commit products to event
  e.put(rpDataOutput, rpDataProductLabel);
  e.put(rpCCOutput, rpCCProductLabel);
  e.put(t2DataTheStripDigis, t2DataProductLabel);
  e.put(t2DataThePadDigis, t2DataProductLabel);
  e.put(t2DataTheVfats,t2DataProductLabel);
  e.put(t2DatavfatsStatusMap,t2DataProductLabel);
  e.put(t1DigiWireColl, t1DataProductLabel);
  e.put(t1DigiVfatColl, t1DataProductLabel);
  e.put(conversionStatus, conversionStatusLabel);

  rpDataEndJob(e);   
}

//-------------------------------------------------------------------------------------------------
//                                           RP CC
//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::rpCCBeginJob() 
{
}

//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::rpCCEndJob(edm::Event& e) 
{
}


//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::rpCCProduce(edm::Event& e, VFATFrameCollection::Iterator &fr, const VFATInfo &info, const VFATAnalysisMask &analysisMask)
{
  // cache for error messages
  stringstream fes;

  const VFATFrame *frame = fr.Data();

  // get IDs TODO - construct RPCCId
  unsigned short symId = info.symbolicID.symbolicID;

  const vector<unsigned char> activeCh = frame->getActiveChannels();
  
  std::bitset<16> bs_even;
  std::bitset<16> bs_odd;

  bs_even.reset();
  bs_odd.reset();

  unsigned int stripNo;

  // if all channels are masked out, do not process all frame
  if (!analysisMask.fullMask) 
    for (unsigned int j = 0; j < activeCh.size(); j++)
      {
  //      std::cout << "Active channel " << j << " value " << (int)(activeCh[j]) << std::endl;
  // check, whether j channel is not masked out
  if (analysisMask.maskedChannels.find(j) == analysisMask.maskedChannels.end())
    {
      stripNo = (unsigned int) (activeCh[j]);
      unsigned int ch = stripNo + 2; // TODO check if +2 is necessary
      //  std::cout << "Strip no " << (unsigned int)(activeCh[j]) << std::endl;
      //  std::cout << "Channel no " << ch << std::endl;
      if (ch >= 72 && ch <= 100 && (ch % 4 == 0)) 
        {
    bs_even.set(ch/4-18);
    continue;
        }
      if (ch >= 40 && ch <= 68 && (ch % 4 == 0)) 
        {
    bs_even.set(ch/4 - 2);
    continue;
        }
      if (ch == 38) {
        bs_odd.set(15);
        continue;
      }
      if (ch >= 104 && ch <= 128 && (ch % 4 == 0)) 
        {
    bs_odd.set(ch/4 - 18);
    continue;
        }
      if (ch >= 42 && ch <= 70 && (ch % 4 == 2)) 
        {
    bs_odd.set((ch-2)/4 - 9 - 1);
    continue;
        }
    } // end if
      } // end for
  //     std::cout << "Odd " << bs_odd << std::endl;
  //      std::cout << "Even " << bs_even << std::endl;

  unsigned int evendetId = TotRPDetId::DecToRawId(symId * 10);
  unsigned int odddetId = TotRPDetId::DecToRawId(symId * 10 + 1);
  RPCCBits ccbits_even( evendetId , bs_even );
  RPCCBits ccbits_odd( odddetId, bs_odd);
  
  rpCCOutput->push_back( ccbits_even );
  rpCCOutput->push_back( ccbits_odd );
} 
  


//-------------------------------------------------------------------------------------------------
//                                             RP DATA
//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::rpDataBeginJob() 
{

}

//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::rpDataEndJob(edm::Event& e) 
{
}


//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::rpDataProduce(edm::Event& e, VFATFrameCollection::Iterator &fr,
  const VFATInfo &info, const VFATAnalysisMask &analysisMask)
{
#ifdef DEBUG
  cout <<  "void Raw2DigiProducer::rpDataProduce(VFATFrameCollection::Iterator &fr, hwID=0x" << hex << info.hwID << dec << ")\n";
#endif
  
  const VFATFrame *frame = fr.Data();
    
  // cache for error messages
  stringstream fes;
 
  // get IDs
  unsigned short symId = info.symbolicID.symbolicID;

  // add RPStripDigi for each hit
  unsigned int detId = TotRPDetId::DecToRawId(symId / 10);
  unsigned short offset = (symId % 10) * 128;
  const vector<unsigned char> activeCh = frame->getActiveChannels();
  DetSet<RPStripDigi> &detSet = rpDataOutput->find_or_insert(detId);

#ifdef DEBUG
  cout << "activeCh.size() = " << activeCh.size() << endl;
#endif

  for (unsigned int j = 0; j < activeCh.size(); j++) {
    // skip masked channels
    if (!analysisMask.fullMask && analysisMask.maskedChannels.find(j) == analysisMask.maskedChannels.end()) {
#ifdef DEBUG
      printf("\tdetSet.push_back(RPStripDigi(detId, offset + activeCh[%d]))\n", j);
#endif
      detSet.push_back(RPStripDigi(detId, offset + activeCh[j]));
    }
  }  
}

//-------------------------------------------------------------------------------------------------
//                                               T1
//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::t1DataBeginJob(const edm::EventSetup &es) 
{
  // for t1 - generating connections
  // get mapping
  ESHandle<DAQMapping> mapping;
  es.get<TotemDAQMappingRecord>().get(mapping);

  map<FramePosition, VFATInfo>::const_iterator iter;
  for(iter = mapping->VFATMapping.begin(); iter != mapping->VFATMapping.end(); iter++)
    {     
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

    //    if(ID__.Arm()==0 && ID__.Plane()==4 && ID__.CSC()==5)
//	cout << stripNumber << " " << vfatCh << " "<< view << endl;


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

//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::t1DataEndJob(edm::Event& e)
{
}

//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::t1DataProduce(edm::Event& e, VFATFrameCollection::Iterator &fr, const VFATInfo &info) 
{

  const VFATFrame *frame = fr.Data();

  unsigned short BX = frame->getBC();
  int BX_ = int(BX);
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
  symId_ = symId;
  symId_ <<= 15;
  id_ |= symId_;

  const vector<unsigned char> activeCh = frame->getActiveChannels();
  //  Int_t VFatFlagSP; Int_t convertedIndex;
  //  int iconv;
  std::auto_ptr<T1DigiWire> t1w;
  std::auto_ptr<T1DigiVfat> t1s;
  
  T1DigiWire t1wv;
  T1DigiVfat t1sv;      
  /* in this case layer is not defined :
     in simulation you have: arm=1bit, plane=3bits, csc=3bits, layer=2bits
     here you have:  arm=1bit, plane=3bits, csc=3bits, anode or cathode=1bit, position=2bits
  */
  for (unsigned int j = 0; j < activeCh.size(); j++)
  {
	  if(activeCh[j]>=128)
	  {
		  cout << "activeCh >= 128" << endl;
		  throw (const char*)"activeCh >= 128";
	  }

	  int electronic_index=0;
	  electronic_index = i*1000 + activeCh[j];
	  T1FisChannel FC2 =T1Map.find(electronic_index)->second;
	  if(FC2.Type()==0 || FC2.Type()==-1 )
	  {
		  /* DetId recalled here was created this way
			 Arm 1 bit
			 Plane 3 bits
			 CSC 3 bits
			 A or C 1 bit (0=A, 1=C)
			 Pos 2 bits (0 to 2)
		   */

		  T1DetId myDetid = FC2.DetId();
		  int channelNumber = FC2.Channel()+1;

		  if( ! _T1deadChannelManager.isChannelDead(myDetid, channelNumber) ){
			  if(myDetid.Layer() == 3 && channelNumber > 0 )
			  {
				  T1DigiWire t1digiwire(channelNumber,BX_);
				  t1DigiWireColl->insertDigi(myDetid,t1digiwire);
			  }

			  if( (myDetid.Layer() == 1 || myDetid.Layer() == 2) && channelNumber > 0)
			  {
				  T1DigiVfat t1digivfat(channelNumber ,0,BX_,FC2.Type());
				  t1DigiVfatColl->insertDigi(myDetid,t1digivfat);
			  }
		  }
	  }
  }
}

//-------------------------------------------------------------------------------------------------
//                                              T2
//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::t2DataBeginJob() 
{
  int framestatus=5;//It means the vfat frame has not been  checked , lost in principle.
  for(unsigned int pl=0;pl<40;pl++)
    for(unsigned int iid=0;iid<17;iid++)
      {
  unsigned int key=pl*100+iid;
  t2DatavfatsStatusMap->insert(pair<unsigned int, unsigned int>(key, framestatus));
      }

}

//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::t2DataEndJob(edm::Event& e) 
{
}

//-------------------------------------------------------------------------------------------------

void Raw2DigiProducer::t2DataProduce(edm::Event& e, VFATFrameCollection::Iterator &iter, const VFATInfo &info, const VFATAnalysisMask &analysisMask) 
{
  const VFATFrame *frame = iter.Data();

  //error buffer
  stringstream fes;

  unsigned int framestatus=4;
  bool T2reading=true;
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinformation;
  T2Mapping VFATconvertDigichannell;

  // get IDs
  unsigned short symId = info.symbolicID.symbolicID;
  
  unsigned short dataId = frame->getChipID();

  // FILL T2Digi AND VFATS
  vector<unsigned char> activeCh = frame->getActiveChannels();

  if((t2DiscardHighOccupancyVfat)&&(activeCh.size()>20))
    T2reading=false;


  if(T2reading)
    {
      framestatus=0;
  
      unsigned int SymID= (symId%100); //Vfat Gius. conv Id  
  
      Int_t VFatFlagSP; Int_t convertedIndex;
      int iconv;
      std::auto_ptr<T2PadDigi> t2p;
      std::auto_ptr<T2StripDigi> t2s;
  
      T2PadDigi t2pv;
      T2StripDigi t2sv;
      int col=6;
      int row=6;
      Bool_t flagconv=true;
  
      planeinformation=conv.GetT2Info(symId/100);
      //This should be consistent also in the xml source for non data vfat

      T2DetId currentDet(planeinformation.cmsswid);
  
      std::auto_ptr<T2DigiVfat> actualvfat(new T2DigiVfat(planeinformation.cmsswid, (Int_t)SymID,dataId, activeCh));
      actualvfat->SetFrameStatus(framestatus);
      t2DataTheVfats->insertDigi(currentDet,(*actualvfat)); 
  
      /*
  if(e.id ().event () < 10)
  {
  map<FramePosition, unsigned int>::const_iterator iter;
  for(iter=Vfat_t2_mapping->readoutPositionToId.begin();iter!=Vfat_t2_mapping->readoutPositionToId.end();iter++)
  {
  std::cout<<"Inserted VFAT in plane"<<planeinformation.symb<<"with "<<activeCh.size()<<" channels"<<std::endl;
  }
      
  }
      */

      if (!analysisMask.fullMask)
        for (unsigned int j = 0; j < activeCh.size(); j++)
        {
        // if channel is not masked
          if (analysisMask.maskedChannels.find(activeCh[j]) != analysisMask.maskedChannels.end()) 
          {
            continue;
          }




      //std::cout<<"rawVFATtoCMSplane[rawIds[i]]=  "<<rawVFATtoCMSplane[rawIds[i]]<< std::endl;
      // T2DetId currentDet(cmsIds[i]);
      //std::cout<<"DetRawId= "<<currentDet.calculateRawId(currentDet.arm(), currentDet.halfTelescope(), currentDet.plane(), currentDet.planeSide())<<std::endl;

      flagconv=VFATconvertDigichannell.convertToGeo((Int_t)SymID, (Int_t)activeCh[j], VFatFlagSP, convertedIndex);

      //HistActiveChannel->Fill(activeCh[j]);
      //Convertion from index-Giuseppe_convenction to row;col is necessary in order to have the same MakeCluster 
      //for Digitization-Simulation and the data.
      
      if(flagconv)    
        {
    if(SymID<=16) //Is a data Vfat, Fill DIGI
      {
        if (VFatFlagSP == 1) {                  //enum fStrip,fPad
          if (SymID >= 2 && SymID <= 14 && activeCh[j] >= 4 && activeCh[j] <= 123) 
      {         
        
        col = (convertedIndex-1)/24;
        row = (convertedIndex-1)%24;
          
        /*Only for Monitor check Warning: Monitor use -4!!But Active channel Pad from 3 to 122!!*/
        /*
          int col2= (activeCh[j] - 3) / 24;
          int row2= (activeCh[j] - 3) % 24;
          col2= col2+(SymID-2)*5;

          PadCOLindex->Fill(col2);
          PadROWindex->Fill(row2);

          Padconvertedindex->Fill(convertedIndex);
        */
        /* End on Monitor check */

        iconv=convertedIndex;
        
        t2p=std::auto_ptr<T2PadDigi>(new T2PadDigi ((int)iconv,(int)row,(int)col,0));   //third parameter=adc, set to 0.
        t2pv= *t2p;
      
        t2DataThePadDigis->insertDigi(currentDet,t2pv);//insert(std::make_pair(currentDet,t2p(iconv,(int)row,(int)col,0)));  
      }
          else
      if((SymID >= 2 && SymID <= 14)==false)
        std::cout<<"Warning: VFat-PAD non-data channel on "<<std::endl;
        }   
        // std::cout<<"--"<<std::endl;
        if (VFatFlagSP == 0) 
          {
      if (SymID == 0) 
        {
          //row=activeCh[j];
          row=convertedIndex%256;
          col= 0;
          iconv=convertedIndex;      
          t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row,(int)col,0));
          t2sv=*t2s;
          t2DataTheStripDigis->insertDigi(currentDet,t2sv); 
        }
      if (SymID == 1) 
        {       
          //row=activeCh[j]; Modif on 24 Jan 2010
          //row=activeCh[j]+128;
          row=convertedIndex%256;
          col= 0;
          iconv=convertedIndex;      
          t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row,(int)col,0));
          t2sv=*t2s;
          t2DataTheStripDigis->insertDigi(currentDet,t2sv); 
        } 
      if (SymID == 15)
        {
        
          //row=activeCh[j]; Modif on 24 Jan 2010
          //row=activeCh[j]+128;
          row=convertedIndex%256;
          col= 1;
          iconv=convertedIndex;      
          t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row,(int)col,0));
          t2sv=*t2s;
          t2DataTheStripDigis->insertDigi(currentDet,t2sv); 
        }
      if (SymID == 16) 
        {
          //row=activeCh[j];
          row=convertedIndex%256;
          col= 1;
          iconv=convertedIndex;      
          t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row,(int)col,0));
          t2sv=*t2s;
          t2DataTheStripDigis->insertDigi(currentDet,t2sv); 
        }
      
      if((SymID!=0)&&(SymID!=1)&&(SymID!=15)&&(SymID!=16))
        std::cout<<"Warning: VFat-STRIP non-data channel on "<<std::endl;      
          }
      }
    else //Trigger VFAT
      {
        std::cout<<"Warning: non-data VFAT with ID "<<SymID<<std::endl;
      }
        }
      else
        {
    std::cout<<"Error: T2Mapping::convertToGeo does not work properly (flagconv=false)!!"<<std::endl;    
        }                  
    }

    }//if T2Reading

  map<unsigned int, unsigned int>::iterator itm = t2DatavfatsStatusMap->find(symId);
  if(itm==t2DatavfatsStatusMap->end())
    std::cout<<"symb error: impossible value"<<std::endl;
  else
    itm->second=framestatus;
}

DEFINE_FWK_MODULE(Raw2DigiProducer);
