/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/Framework/interface/InputSource.h"

#include "TotemRawDataLibrary/Readers/interface/DataFile.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "DataFormats/Provenance/interface/EventID.h"


namespace edm {
    class Event;
    class ParameterSet;
    class RunPrincipal;
    class LuminosityBlockPrincipal;
    class EventPrincipal;
}

namespace Totem {
    class RawEvent;
    class FramePosition;
}


/**
 * \ingroup TotemRawData
 * \brief Reads in TOTEM raw data file(s) and puts Totem::RawEvent into the Event structure.
 *
 * Each input file is treated as one run. If setRunNumberFromFileName is set to false, the runs
 * are enumerated consecutively: 1, 2, 3... If setRunNumberFromFileName is true, run numbers are
 * extracted from file names. The following file name format is expected:
 * \verbatim path/some text<run>some text<file>some text \endverbatim
 * The "/" denotes the last slash in the path. <run> and <file> contain digits only, <run> must
 * be greater than 1, <file> can have 4 digits only (i.e. it must smaller than 10000). <file> is
 * optional, if missing, 0 will be used instead. The (CMSSW) run number is then composed as
 * \verbatim run number = <run> * 1E4 + <file> \endverbatim
 *
 * The module starts a new luminosity block whenever the UNIX timestamp of raw event changes.
 * In CMSSW 7.x it is essential for synchronisation with EventSetup: it can only be updated at
 * luminosity-block boundaries.
**/

using namespace edm;
using namespace std;
using namespace Totem;

class RawDataSource : public InputSource
{
  public:
    struct FileInfo
	{
      std::string fileName;         ///< path to the file
      Totem::DataFile *file;        ///< instance of appropriate reader class
      unsigned int collectionIdx;   ///< index in the 'collections' list
      edm::RunNumber_t runNumber;   ///< associated run number
    };

    RawDataSource(const edm::ParameterSet &, const edm::InputSourceDescription&);
    virtual ~RawDataSource();
    virtual void beginJob();
    virtual void endJob();

  protected:
    unsigned int verbosity;

    std::vector<std::string> fileNames;                     ///< vector of raw data files names
    bool setRunNumberFromFileName;                          ///< whether run number shall be set from file names (see GetRunNumberFromFileName method)
    unsigned int printProgressFrequency;                    ///< frequency with which the progress (i.e. event number) is to be printed

    unsigned int fileIdx;                                   ///< current file index (within files), counted from 0

    friend class ProcessManager;
    std::vector<FileInfo> files;                            ///< to keep information about opened files
    std::vector<Totem::VFATFrameCollection *> collections;  ///< the list of instantiated data collections

    /// \brief extracts run number from a file name
    edm::RunNumber_t GetRunNumberFromFileName(const std::string &str);

  private:
    /// list of the next state items
    deque<ItemType> items;
    
    /// tries to load a next raw event and updates the list of next states 'items'
    void LoadRawDataEvent();

    /// called by the framework to determine the next state (run, lumi, event, stop, ...)
    /// here it simply returns the popped state from the 'items' queue
    virtual ItemType getNextItemType();
    
 /// called by the framework for item type 'IsRun'
    virtual boost::shared_ptr<RunAuxiliary> readRunAuxiliary_();
  
/// called by the framework for item type 'IsLumi'
    virtual boost::shared_ptr<LuminosityBlockAuxiliary> readLuminosityBlockAuxiliary_();
  
  /// called by the framework for item type 'IsEvent'
    /// the cached (pre-loaded) raw data (currentRawEvent) are inserted into the event
    virtual void readEvent_(EventPrincipal& eventPrincipal);

    // TODO: remove unnecessary, describe
    virtual void skip(int offset);
    virtual void setRun(RunNumber_t r);
    virtual void rewind_();
    void advanceToNext() ;
    void retreatToPrevious();

    /// pre-loaded raw event
    auto_ptr<RawEvent> currentRawEvent;

    /// ID of the current (next) event
    EventID eventID;

    /// UNIX timestamp of the previous event
    time_t previousTimestamp;
};
