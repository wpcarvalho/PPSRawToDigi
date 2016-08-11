/****************************************************************************
 * Seyed Mohsen Etesami
 ****************************************************************************/

#include "EventFilter/CTPPSRawToDigi/interface/DiamondRawDataUnpacker.h"

#include "FWCore/MessageLogger/interface/MessageLogger.h"

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

DiamondRawDataUnpacker::DiamondRawDataUnpacker(const edm::ParameterSet &conf)
{
}

//----------------------------------------------------------------------------------------------------

int DiamondRawDataUnpacker::Run(int fedId, const FEDRawData &data, vector<DiamondFEDInfo> &fedInfoColl, DiamondVFATFrameCollection &coll) const
{
  unsigned int size_in_words = data.size() / 8; // bytes -> words
  if (size_in_words < 2)
    {
      LogProblem("Totem") << "Error in RawDataUnpacker::Run > " <<
	"Data in FED " << fedId << " too short (size = " << size_in_words << " words).";
      return 1;
    }

  fedInfoColl.push_back(DiamondFEDInfo(fedId));

  return ProcessOptoRxFrame((const word *) data.data(), size_in_words, fedInfoColl.back(), &coll);
}

//----------------------------------------------------------------------------------------------------

int DiamondRawDataUnpacker::ProcessOptoRxFrame(const word *buf, unsigned int frameSize, DiamondFEDInfo &fedInfo, DiamondVFATFrameCollection *fc) const
{ 
  // get OptoRx metadata
  unsigned long long head = buf[0];
  unsigned long long foot = buf[frameSize-1];

  fedInfo.setHeader(head);
  fedInfo.setFooter(foot);

  unsigned int BOE = (head >> 60) & 0xF;
  unsigned int H0 = (head >> 0) & 0xF;

  //unsigned long LV1 = (head >> 32) & 0xFFFFFF;
  //unsigned long BX = (head >> 20) & 0xFFF;
  unsigned int FEDId = (head >> 8) & 0xFFF;

  unsigned int FOV = (head >> 4) & 0xF;

  unsigned int EOE = (foot >> 60) & 0xF;
  unsigned int F0 = (foot >> 0) & 0xF;
  unsigned int FSize = (foot >> 32) & 0x3FF;

  // check header and footer structure
  if (BOE != 5 || H0 != 0 || EOE != 10 || F0 != 0 || FSize != frameSize)
    {
      LogProblem("Totem") << "Error in DiamondRawDataUnpacker::ProcessOptoRxFrame > " << "Wrong structure of OptoRx header/footer: "
			  << "BOE=" << BOE << ", H0=" << H0 << ", EOE=" << EOE << ", F0=" << F0
			  << ", size (OptoRx)=" << FSize << ", size (DATE)=" << frameSize
			  << ". FEDID=" << FEDId << ". Skipping frame." << endl;
      return 0;
    }

#ifdef DEBUG
  printf(">> DiamondRawDataUnpacker::ProcessOptoRxFrame > OptoRxId = %u, BX = %lu, LV1 = %lu, frameSize = %u, subFrames = %u)\n",
	 OptoRxId, BX, LV1, frameSize, subFrames);
#endif

  // parallel or serial transmission?

  if (FOV == 2)
    return ProcessOptoRxFrameParallel(buf, frameSize, fedInfo, fc);

  LogProblem("Totem") << "Error in DiamondRawDataUnpacker::ProcessOptoRxFrame > " << "Unknown FOV = " << FOV << endl;

  return 0;
}

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

int DiamondRawDataUnpacker::ProcessOptoRxFrameParallel(const word *buf, unsigned int frameSize, DiamondFEDInfo &fedInfo, DiamondVFATFrameCollection *fc) const
{
  // get OptoRx metadata
  unsigned long long head = buf[0];
  //unsigned int OptoRxId = (head >> 8) & 0xFFF;
  unsigned int FEDId = (head >> 8) & 0xFFF;

  // recast data as buffer or 16bit words, skip header
  const uint16_t *payload = (const uint16_t *) (buf + 1);

  // read in OrbitCounter block
  const uint32_t *ocPtr = (const uint32_t *) payload;
  fedInfo.setOrbitCounter(*ocPtr);
  payload += 2;

  // size in 16bit words, without header, footer and orbit counter block
  unsigned int nWords = (frameSize-2) * 4 - 2;

  // process all VFAT data
  for (unsigned int offset = 0; offset < nWords;)
    {
      unsigned int wordsProcessed = ProcessVFATDataFED(payload + offset, FEDId, fc);
   
      offset += wordsProcessed;
    }

  return 0;
}

//----------------------------------------------------------------------------------------------------
int DiamondRawDataUnpacker::ProcessVFATDataFED(const uint16_t *buf, unsigned int FEDId, DiamondVFATFrameCollection *fc) const
{   
  // start counting processed words
  unsigned int wordsProcessed = 1;

  // padding word? skip it
  if (buf[0] == 0xFFFF)
    return wordsProcessed;

  // check header flag
  unsigned int hFlag = (buf[0] >> 8) & 0xFF;
  if (hFlag != vmCluster && hFlag != vmRaw)
    {
      LogProblem("Totem") << "Error in DiamondRawDataUnpacker::ProcessVFATDataParallel > "
			  << "Unknown header flag " << hFlag << ". Skipping this word." << endl; 
      return wordsProcessed;
    }

  // compile frame position
  // NOTE: DAQ group uses terms GOH and fiber in the other way
  unsigned int gohIdx = (buf[0] >> 4) & 0xF;
  unsigned int fiberIdx = (buf[0] >> 0) & 0xF;

  // prepare temporary VFAT frame
  DiamondVFATFrame f;
  DiamondVFATFrame::word *fd = f.getData();

  // copy footprint, BC, EC, Flags, ID, if they exist
  uint8_t presenceFlags = 0;

  if (((buf[wordsProcessed] >> 12) & 0xF) == 0xA)  // BC
    {
      presenceFlags |= 0x1;
      fd[11] = buf[wordsProcessed];
      wordsProcessed++;
    }

  if (((buf[wordsProcessed] >> 12) & 0xF) == 0xC)  // EC, flags
    {
      presenceFlags |= 0x2;
      fd[10] = buf[wordsProcessed];
      wordsProcessed++;
    }

  if (((buf[wordsProcessed] >> 12) & 0xF) == 0xE)  // ID
    {
      presenceFlags |= 0x4;
      fd[9] = buf[wordsProcessed];
      wordsProcessed++;
    }

  // save offset where channel data start
  unsigned int dataOffset = wordsProcessed;

  //Assign diamond frame position 
  DiamondFramePosition fp(FEDId, gohIdx, fiberIdx);


  if (hFlag == vmRaw)
    wordsProcessed += 9;

  // process trailer

  bool skipFrame = false;



  wordsProcessed++;


  if (skipFrame)
    return wordsProcessed;


  // get channel data and CRC - raw mode
  if (hFlag == vmRaw)
    {
      for (unsigned int i = 0; i < 8; i++)
	fd[8 - i] = buf[dataOffset + i];

      // copy CRC
      presenceFlags |= 0x8;
      fd[0] = buf[dataOffset + 8];
    }

  // save frame to output
  f.setPresenceFlags(presenceFlags);
  fc->Insert(fp, f);

  return wordsProcessed;
}

