/****************************************************************************
*
* This is a part of the TOTEM testbeam/monitoring software.
* This is a part of the TOTEM offline software.
* Authors: 
*   Jan Kašpar (jan.kaspar@gmail.com) 
*
****************************************************************************/

#include "TotemRawDataLibrary/DataFormats/interface/CommonDef.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrameCollection.h"
#include "TotemRawDataLibrary/Readers/interface/DataFile.h"
#include "TotemRawDataLibrary/Utilities/interface/CablingAnalyzer.h"

#ifdef _MONITOR_
  #include "VFATtest/interface/VFAT2RegistersXML.h"
#else
  // TODO: solve the problem with missing VFAT2RegistersXML
#endif

#include <iostream>
#include <map>

using namespace Totem;
using namespace std;


#ifdef _MONITOR_
unsigned short CablingAnalyzer::run(Totem::DataFile *data, VFAT2RegistersXML *conf) const
{
  if (!data) return 1;

  Totem::RawEvent *event = data->CreateEvent();
  unsigned short prevEC = 0;
  record result, prevResult;

  FramePosition masterIdx;
  bool masterFound = false;

  unsigned char mode = 0;  // 0... waiting for first configuration block, 1...first conf. block, 2...successive conf. blocks

  unsigned short ID = 0;  // current 16bit ID
  signed int eventCount = -1;

  unsigned int chipEvents = 0;
  unsigned int chipErrors = 0;
  unsigned int chipWarnings = 0;

  map<FramePosition, unsigned short> results;

  while ( !data->GetNextEvent(event) ) {
    eventCount++;

    // get active chips
    vector<record> active;
    for (VFATFrameCollection::Iterator it(event->frames); !it.IsEnd(); it.Next()) {
      if (!it.Data()->checkFootprint()) continue;
      if (!it.Data()->checkCRC()) continue;

      active.push_back(record(it.Data()->getChipID(), it.Position()));
    }

    // set master position (if not known yet)
    if (!masterFound) {
      if (active.size() == 0) {
        WARN("CablingAnalyzer::run") << "There should be exactly one active (master) VFAT frame in data. I read " <<  active.size() << " frames." << c_endl;
        continue;
      }
      if (active.size() > 1) {
        WARN("CablingAnalyzer::run") << "There should be exactly one active (master) VFAT frame in data. I read " <<  active.size() << " frames." << c_endl;
        return 2;
      }
      masterIdx = active[0].index;
      masterFound = true;
    }

    // get master EC
    const VFATFrame *masterFrame = event->frames->GetFrameByIndex(masterIdx);
    if (!masterFrame || !masterFrame->checkFootprint() || !masterFrame->checkCRC()) {
      WARN("CablingAnalyzer::run") << "Master VFAT frame not found (event " << eventCount << ")." << c_endl;
      continue;
    }
    unsigned short EC = masterFrame->getEC();

    // EC processing
    bool confBegin = false;
    if (EC == 0) {
      // print result
      if (mode > 0) {
        printf(">> ID %04x: events analyzed %i, errors %i, warnings %i", ID, chipEvents, chipErrors, chipWarnings);
        if (chipErrors == 0 && (chipEvents > chipErrors + chipWarnings)) {
          results[result.index] = ID;
          cout << ", position " << result.index;
        }
        printf("\n");
      }

      // shift mode
      if (mode == 1) mode = 2;
      if (mode == 0) mode = 1;

      // load new configuration
      if (conf->LoadNextConfiguration()) {
        ERROR("CablingAnalyzer::run") << "There is no more configuration as requested by data file (number "
          << conf->getConfigurationNumber() << ")." << c_endl;
        return 6;
      }

      // get the ID for VFAT being searched
      const std::vector<VFAT2Registers> &reg = conf->GetRegisters();
      if (reg.size() != 1) {
        ERROR("CablingAnalyzer::run") << "Wrong number (" << reg.size() << ") of VFAT records in I2C file." << c_endl;
        return 7;
      }
      ID = reg[0].GetFullChipID();
      printf("> now reading data for VFAT ID %04x\n", ID);

      // reset counters, flag
      confBegin = true;
      chipEvents = 0;
      chipErrors = 0;
      chipWarnings = 0;
    } else {
      // check EC progress
      if (prevEC + 1 != EC) {
        WARN("CablingAnalyzer::run") << "Wrong EC progress. Previous EC = " << prevEC << ", current EC = " << EC << "." << c_endl;
        //return 3;
      }
    }
    prevEC = EC;

    chipEvents++;

    // analyze VFATs data
    unsigned int hits = 0;  
    unsigned int ID12 = ID & 0xFFF;
    for (unsigned int i = 0; i < active.size(); i++) {
      if (mode > 1 && active[i].index == masterIdx) continue;
      if (active[i].ID != ID12) continue;

      result = active[i];
      hits++;
    }

    if (hits == 0) {
      WARN("CablingAnalyzer::run") << "No VFAT frame corresponding to ID " << hex << ID << dec << " (event " << eventCount << ")." << c_endl;
      chipWarnings++;
      continue;
    }

    if (hits > 1) {
      WARN("CablingAnalyzer::run") << hits << " VFAT frame corresponding to ID " << hex << ID << dec << " (event " << eventCount << ")." << c_endl;
      chipWarnings++;
      continue;
    }

    if (hits == 1) {
      if (!confBegin && !(result == prevResult)) {
        WARN("CablingAnalyzer::run") << "Previous and current positions for the same " << hex << ID << dec << 
          " are different (event " << eventCount << ")." << c_endl;
        chipErrors++;
        continue;
      }
      /*
      if (confBegin || !(result == prevResult)) {
        printf("ID = %04x\tposition = ", ID);
        cout << result.index << endl;
      }
      */
      
      prevResult = result;
    }
  }

  // print results in XML format
  printf("\n\n");
  for (map<FramePosition, unsigned short>::iterator it = results.begin(); it != results.end(); ++it) {
    printf("<vfat id=\"0x%04x\" idxinfiber=\"%i\" gohid=\"%i\" optorxid=\"%i\" totfedid=\"%i\" subsystemid=\"%i\" />", it->second,
      it->first.GetIdxInFiber(),it->first.GetGOHId(),it->first.GetOptoRxId(),it->first.GetTOTFEDId(),it->first.GetSubSystemId());
  }
  
  delete event;
  return 0;
}

#endif
