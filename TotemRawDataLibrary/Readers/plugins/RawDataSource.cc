/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jan Kašpar (jan.kaspar@gmail.com)
 *
 ****************************************************************************/

// this include is a dirty trick
#include "TotemRawDataLibrary/Readers/plugins/Event_hacked.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/InputSourceMacros.h"
#include "DataFormats/Provenance/interface/RunAuxiliary.h"
#include "DataFormats/Provenance/interface/LuminosityBlockAuxiliary.h"
#include "DataFormats/Provenance/interface/EventAuxiliary.h"
#include "DataFormats/Provenance/interface/ProcessHistoryRegistry.h"
#include "FWCore/Framework/interface/EventPrincipal.h"
#include "FWCore/Framework/interface/InputSource.h"

#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/Readers/interface/DataFile.h"

#include <iostream>
#include <iomanip>

#include "TotemRawDataLibrary/Readers/interface/RawDataSource.h"

using namespace std;
using namespace edm;
using namespace Totem;


//#define DEBUG 1

//----------------------------------------------------------------------------------------------------

RawDataSource::RawDataSource(const edm::ParameterSet& pSet, const edm::InputSourceDescription& desc) :
  InputSource(pSet, desc),
  verbosity(pSet.getUntrackedParameter<unsigned int>("verbosity", 0)),
  fileNames(pSet.getUntrackedParameter<vector<string> >("fileNames")),
  setRunNumberFromFileName(pSet.getParameter<bool>("setRunNumberFromFileName")),
  printProgressFrequency(pSet.getUntrackedParameter<unsigned int>("printProgressFrequency", 0)),
  eventID(0, 0, 0),
  previousTimestamp(0)
{
#ifdef DEBUG
  cout << ">> RawDataSource::RawDataSource" << endl;
#endif
  produces<Totem::RawEvent>();
}

//----------------------------------------------------------------------------------------------------

RawDataSource::~RawDataSource()
{
}

//----------------------------------------------------------------------------------------------------
//boost::shared_ptr<RunAuxiliary> RawDataSource::readRunAuxiliary_()
std::shared_ptr<RunAuxiliary> RawDataSource::readRunAuxiliary_()
{
#ifdef DEBUG
  printf(">> RawDataSource::readRunAuxiliary_\n");
#endif

  Timestamp ts_beg(currentRawEvent->timestamp << 32);
  Timestamp ts_end(Timestamp::endOfTime().value() - 0);

  return std::shared_ptr < RunAuxiliary > (new RunAuxiliary(eventID.run(), ts_beg, ts_end));
}

//----------------------------------------------------------------------------------------------------
//boost::shared_ptr<LuminosityBlockAuxiliary> RawDataSource::readLuminosityBlockAuxiliary_()
std::shared_ptr<LuminosityBlockAuxiliary> RawDataSource::readLuminosityBlockAuxiliary_()
{
#ifdef DEBUG
  printf(">> RawDataSource::readLuminosityBlockAuxiliary_\n");
#endif

  // we create luminosity blocks of 1s duration
  Timestamp ts_beg(currentRawEvent->timestamp << 32);
  Timestamp ts_end(((currentRawEvent->timestamp + 1) << 32) - 1);

  return std::shared_ptr <LuminosityBlockAuxiliary> (new LuminosityBlockAuxiliary(eventID.run(),
    eventID.luminosityBlock(), ts_beg, ts_end));
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::readEvent_(EventPrincipal& eventPrincipal)
{
#ifdef DEBUG
  printf(">> RawDataSource::readEvent_ : %lu, %lu | %u\n", currentRawEvent->dataEventNumber, currentRawEvent->timestamp, eventID.event());
#endif
  
  // create Event structure and fill with auxiliary data
  Timestamp ts(currentRawEvent->timestamp << 32);  // conversion from UNIX timestamp: see src/DataFormats/Provenance/interface/Timestamp.h
  bool isRealData = true;
  EventAuxiliary::ExperimentType expType(EventAuxiliary::Undefined);
  EventAuxiliary aux(eventID, processGUID(), ts, isRealData, expType);

  ProcessHistoryRegistry phr;
  eventPrincipal.fillEventPrincipal(aux, phr);
  ModuleCallingContext const* mcc = NULL;
  Event e(eventPrincipal, moduleDescription(), mcc);

  // put raw data into the event
  e.put(currentRawEvent);
  e.commit_();

  // TODO
  //if (printProgressFrequency > 0 && (event() % printProgressFrequency) == 0)
  //  cout << "\033[25Devent " << run() << ":" << event();
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::LoadRawDataEvent()
{
#ifdef DEBUG
  printf(">> RawDataSource::LoadRawDataEvent\n");
#endif

  // prepare structure for the raw event
  currentRawEvent = auto_ptr<RawEvent>(new RawEvent(collections[files[fileIdx].collectionIdx], false));

  // load next raw event
  bool newFile = false;
  while (true)
  {
    unsigned int result = files[fileIdx].file->GetNextEvent(currentRawEvent.get());

    if (result == 0)
      break;
    else {
      // try moving to next file
      fileIdx++;
      newFile = true;

      // stop if there are no more files
      if (fileIdx >= files.size())
      {
        items.push_back(IsStop);
        return;
      }
    }
  }

#ifdef DEBUG
  printf("\t%lu, %lu\n", currentRawEvent->dataEventNumber, currentRawEvent->timestamp);
#endif

  bool beginning = (eventID.run() == 0);

  if (newFile || beginning)
  {
    items.push_back(IsRun);
    items.push_back(IsLumi);
    items.push_back(IsEvent);

    eventID = EventID(files[fileIdx].runNumber, 1, 1);

    previousTimestamp = currentRawEvent->timestamp;

    return;
  } 

  // new luminosity block ??
  if (currentRawEvent->timestamp != previousTimestamp)  // 1s resolution
  {
    items.push_back(IsLumi);
    eventID = EventID(eventID.run(), eventID.luminosityBlock() + 1, eventID.event()); // advance lumi number
  }

  previousTimestamp = currentRawEvent->timestamp;

  eventID = EventID(eventID.run(), eventID.luminosityBlock(), eventID.event() + 1); // advance event number
  items.push_back(IsEvent);
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::skip(int offset)
{
  throw cms::Exception("RawDataSource") << "Method 'skip' not implemented." << std::endl;
  /*
  for (; offset < 0; ++offset)
  {
    retreatToPrevious();
  }

  for (; offset > 0; --offset)
  {
    advanceToNext();
  }
  */
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::setRun(RunNumber_t r)
{
  throw cms::Exception("RawDataSource") << "Method 'setRun' not implemented." << std::endl;

  /*
  // No need to check for invalid (zero) run number,
  // as this is a legitimate way of stopping the job.
  // Do nothing if the run is not changed.
  if (r != eventID_.run())
  {
    eventID_ = EventID(r, origEventID_.luminosityBlock(), 0);
    newRun_ = newLumi_ = true;
  }
  */
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::rewind_()
{
  throw cms::Exception("RawDataSource") << "Method 'rewind_' not implemented." << std::endl;
  /*
  eventID_ = EventID(eventID_.run(), eventID_.luminosityBlock(), origEventID_.event());
  newRun_ = newLumi_ = true;
  */
}

//----------------------------------------------------------------------------------------------------

#ifdef DEBUG
InputSource::ItemType prevItem = InputSource::ItemType::IsInvalid;
#endif

InputSource::ItemType RawDataSource::getNextItemType()
{
  if (items.empty())
    LoadRawDataEvent();

  ItemType item = items.front();
  items.pop_front();

#ifdef DEBUG
  if (prevItem != item)
    printf(">> RawDataSource::getNextItemType > %u --> %u (run=%u, lumi block=%u, event=%u)\n", prevItem, item,
      eventID.run(), eventID.luminosityBlock(), eventID.event());

  prevItem = item;
#endif

  return item;
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::advanceToNext()
{
  throw cms::Exception("RawDataSource") << "Method 'advanceToNext' not implemented." << std::endl;
  //eventID_ = eventID_.next(eventID_.luminosityBlock());
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::retreatToPrevious()
{
  throw cms::Exception("RawDataSource") << "Method 'retreatToPrevious' not implemented." << std::endl;
  //eventID_ = eventID_.previous(eventID_.luminosityBlock());
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::beginJob()
{
#ifdef DEBUG
  printf(">> RawDataSource::beginJob\n");
#endif

  // try to open all files
  for (unsigned int i = 0; i < fileNames.size(); i++)
  {
	if(verbosity)
    {
		cout<< "\n>> RawDataSource::beginJob > Opening file `" << fileNames[i] << "'." << endl;
	}

    // add prefix
    FileInfo fi;
    char *cmsswPath = getenv("CMSSW_BASE");
    fi.fileName = fileNames[i];
    fi.fileName.erase(0, fi.fileName.find_first_not_of(' '));
    if (fi.fileName[0] != '/' && fi.fileName.find("./") != 0 && fi.fileName.find("://") == string::npos)
      if (cmsswPath)
        fi.fileName = string(cmsswPath) + string("/src/") + fi.fileName;

    // run number
    if (setRunNumberFromFileName)
      fi.runNumber = GetRunNumberFromFileName(fi.fileName);
    else
      fi.runNumber = i + 1;

    // open file
    fi.file = Totem::DataFile::OpenStandard(fi.fileName);
    if (!fi.file)
      throw cms::Exception("RawDataSource") << "Cannot open file " << fileNames[i] << std::endl;

    // start indexing
    fi.file->StartIndexing();

    // instantiate appropriete collection object (only once per type)
    bool found = false;
    for (unsigned int j = 0; j < collections.size(); j++)
    {
      if (fi.file->IsCollectionCompatible(collections[j]))
      {
        fi.collectionIdx = j;
        found = true;
        break;
      }
    }

#ifdef DEBUG
    printf("\tcollection found = %u\n", found);
#endif

    if (!found)
    {
      fi.collectionIdx = collections.size();
#ifdef DEBUG
      printf("\tcreating new collection\n");
#endif
      collections.push_back(fi.file->CreateCollection());
    }

    files.push_back(fi);
  }

  if (!files.size())
    throw cms::Exception("RawDataSource") << "No files to read." << std::endl;

  // print file info
  if (verbosity > 0)
  {
    printf(">> RawDataSource::beginJob > files opened:\n");
    printf("idx|    run.file | file name                                                    | file type  | collection type\n");
    for (unsigned int i = 0; i < files.size(); i++)
    {
      unsigned int rn = files[i].runNumber / 10000;
      unsigned int fn = files[i].runNumber % 10000;
      unsigned int len = files[i].fileName.size();
      string shortFileName = (len <= 60) ? files[i].fileName : string("...") + files[i].fileName.substr(len - 57);
      printf("%2u | %6u.%04u | %60s | %10s | %25s", i, rn, fn, shortFileName.c_str(),
          files[i].file->GetClassName().c_str(), collections[files[i].collectionIdx]->GetClassName().c_str());
      cout << endl;
    }
  }

  // initialize the state of reader
  fileIdx = 0;

  items.push_back(IsFile);  // needed for the logic in InputSource::nextItemType 
}

//----------------------------------------------------------------------------------------------------

void RawDataSource::endJob()
{
#ifdef DEBUG
  printf(">> RawDataSource::endJob\n");
#endif

  for (unsigned int j = 0; j < collections.size(); j++)
  {
#ifdef DEBUG
    printf("\t%u\t%p\n", j, collections[j]);
#endif
    delete collections[j];
  }
}

//----------------------------------------------------------------------------------------------------

edm::RunNumber_t RawDataSource::GetRunNumberFromFileName(const string &path)
{
  size_t pos = path.find_last_of('/');
  string s = (pos == string::npos) ? path : path.substr(pos + 1);

  string digits = "0123456789";

  size_t p1 = s.find_first_of(digits);
  if (p1 == string::npos)
    throw cms::Exception("RawDataSource::GetRunNumberFromFileName")
        << "Run number cannot be extracted from file name `" << path << "'." << endl;
  size_t p2 = s.find_first_not_of(digits, p1);
  string runStr = (p2 == string::npos) ? s.substr(p1) : s.substr(p1, p2 - p1);

  unsigned int rn = atoi(runStr.c_str());
  unsigned int fn = 0;

  if (p2 != string::npos)
  {
    p1 = s.find_first_of(digits, p2);
    if (p1 != string::npos)
    {
      p2 = s.find_first_not_of(digits, p1);
      string fileStr = (p2 == string::npos) ? s.substr(p1) : s.substr(p1, p2 - p1);
      fn = atoi(fileStr.c_str());
    }
  }

  if (rn < 1)
    throw cms::Exception("RawDataSource::GetRunNumberFromFileName") << "Run number `" << rn
        << "' must be greater than 1." << endl;

  if (fn > 9999)
    throw cms::Exception("RawDataSource::GetRunNumberFromFileName") << "File number `" << fn
        << "' must be smaller than 10000." << endl;

  RunNumber_t ret = rn * 10000 + fn;
  return ret;
}

DEFINE_FWK_INPUT_SOURCE(RawDataSource);
