#include "TotemAnalysis/TotemNtuplizer/interface/RawDataFormats.h"
#include "TotemAnalysis/TotemNtuplizer/interface/TriggerDataFormats.h"
#include "TotemAnalysis/TotemNtuplizer/interface/RPRootTrackInfo.h"
//#include "TotemAnalysis/T2MakeNtples/interface/T2Event.h"
#include "TotemAnalysis/TotemNtuplizer/interface/T1Event.h"
#include "TotemAnalysis/TotemNtuplizer/interface/T2Event.h"


#ifdef __CINT__

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ class EventMetaData+;
#pragma link C++ class TriggerData+;

#pragma link C++ class RPRootDumpTrackInfo+;
#pragma link C++ class std::vector<RPRootDumpTrackInfo>+;

#pragma link C++ class RPRootDumpDigiInfo+;
#pragma link C++ class std::vector<RPRootDumpDigiInfo>+;

#pragma link C++ class RPRootDumpReconstructedProton+;
#pragma link C++ class RPRootDumpReconstructedProtonPair+;

#pragma link C++ class RPRootDumpPattern;
#pragma link C++ class std::vector<RPRootDumpPattern>;
#pragma link C++ class RPRootDumpPatternInfo;

#pragma link C++ class RPRootDumpJet;
#pragma link C++ class std::vector<RPRootDumpJet>;
#pragma link C++ class RPRootDumpJetInfo;
#pragma link C++ class RPRootDiffMassInfo;

#pragma link C++ class T2Event+;
#pragma link C++ class T1Event+;

#endif
