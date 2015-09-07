/**********************************************************
*
* This is a part of TOTEM testbeam/monitoring software.
* This is a part of TOTEM offline software.
* Authors:
*   Michal Zmuda (m.zmuda@aol.com)
*
**********************************************************/

/**
 * \defgroup TotemRawDataLibrary TOTEM raw data package
 * \brief The classes for reading TOTEM raw data files and streams.
 * 
 * 
 * \section Design
 * The classes are designed for the following concept.
 * \li A \b file is a sequence of events. Space between events may be filled with rubbish.
 * \li An \b event consists of a list of VFAT frames and additional information.
 * 
 * There might be several data formats. For each one there shall be a class reading corresponding data files. Common interface for those
 * classes is defined in Totem::DataFile.
 * 
 * An event is represented by Totem::RawEvent class. As already said it includes collection of VFAT frames (class Totem::OldVFATFrameColl) and
 * additional information.
 * 
 * \section RecommendedUsage Recommended usage
 * The following code opens a file (with one of the standard data formats), reads it event by event and for each event it prints ID of the
 * first VFAT frame in the event.
 * \code
 *   DataFile* input = DataFile::OpenStandard("some file");
 *   if (!input) {
 *     printf("Error in opening file\n");
 *   }
 * 
 *   // This creates the event object with a VFAT frame
 *   // collection that is compatible with the file reader.
 *   RawEvent *event = input->CreateEvent();
 * 
 *   input->StartIndexing();
 *   while (!input->GetNextEvent(event)) {
 *     // event processing here
 *   }
 * 
 *   delete event;
 * \endcode
 * The virtue of using Totem::DataFile::OpenStandard is that it first determines type of the file and then it creates instance of appropriate
 * class (e.g. Totem::SlinkFile, Totem::TTPFile ...). The command <c>input->GetNextEvent()</c> finds the next event in the file and
 * reads data (and performs necessary transformations) into \c event.
 * 
 * Strictly speaking, the <c>input->StartIndexing();</c> is not necessary. However, it is useful because when the loop finishes, the file is
 * indexed. And hence one can use Totem::DataFile::GetEvent function.
 * 
**/

#ifndef _Totem_VFATFrame_h_
#define _Totem_VFATFrame_h_

#include <vector>

/**
 * \ingroup TotemRawDataLibrary
 * Interface to data frames from VFAT chip.
 **/
namespace Totem {

class VFATFrame
{
  public:
    typedef unsigned short word;

    virtual ~VFATFrame() {};

    virtual word getBC() const = 0;                                    ///< Returns Bunch Crossing number (BC<11:0>)
    virtual word getEC() const = 0;                                    ///< Returns Event Counter (EV<7:0>)
    virtual word getFlags() const = 0;                                 ///< Returns flags
    virtual word getChipID() const = 0;                                ///< Returns ChipID (ChipID<11:0>)
    virtual bool checkFootprint() const = 0;                           ///< checks the dummy bits, returns true if they're as they should be
    virtual bool checkCRC() const = 0;                                 ///< checks the if CRC is ok (or if there is none)
    virtual bool channelActive(unsigned char channel) const = 0;       ///< Checks if channel number 'channel' was active, returns positive number if it was active, 0 otherwise
    virtual std::vector<unsigned char> getActiveChannels() const = 0;  ///<\brief Returns array of active channels
                                                                       ///< it's MORE EFFICIENT than the previous method for events with low channel occupancy

    virtual void Print(bool binary = false) const = 0;                 ///< Prints the frame. If binary is true, binary format is used.
};                                                                     

}
#endif
