#include <vector>

#include "TH1D.h"
#include "TH2D.h"



#include "L1TriggerTotem/CoincidenceChip/interface/PreTriggerStat.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStat.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStat2.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TriggerStatCollection.h"
#include "L1TriggerTotem/CoincidenceChip/interface/TwoHalfPot.h"
#include "L1TriggerTotem/CoincidenceChip/interface/RawVsSimuPotComparator.h"
#include "L1TriggerTotem/CoincidenceChip/interface/PotCollection.h"
#include "L1TriggerTotem/CoincidenceChip/interface/RPTriggerAnalyzerInfoCollector.h"



// #include <boost/ptr_container/ptr_vector.hpp>

#ifdef __CINT__
#ifndef L1TriggerTotemCoincidenceChip_LINKDEF
#define L1TriggerTotemCoincidenceChip_LINKDEF


#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;

#pragma link C++ nestedclasses;
#pragma link C++ nestedtypedefs; 
// #pragma link C++ class std::vector<double>+;
// #pragma link C++ class std::vector< std::vector<double> >+;
// #pragma link C++ class std::vector< std::vector< std::vector<double> > >+;
// #pragma link C++ class std::vector< std::vector< std::vector< std::vector<double> > > >+;
#pragma link C++ enum ERawOrSimu+;
#pragma link C++ enum EDetOrientation+;
#pragma link C++ class PreTriggerStat+;
#pragma link C++ class TriggerStat+;
#pragma link C++ class TriggerStat2+;
#pragma link C++ class std::vector< TriggerStat*>+;
#pragma link C++ class TriggerStatCollection<TriggerStat>+;
#pragma link C++ class std::vector<TriggerStat2*>+;
#pragma link C++ class TriggerStatCollection<TriggerStat2>+;
#pragma link C++ class TwoHalfPot+;
#pragma link C++ class std::vector<TH1D*>+;
#pragma link C++ class std::vector<TH2D*>+;

#pragma link C++ class RawVsSimuPotComparator+;
// #pragma link C++ class boost::ptr_vector<RawVsSimuPotComparator>+;
#pragma link C++ class std::vector<RawVsSimuPotComparator*>+;
#pragma link C++ class PotCollection+;
#pragma link C++ class RPTriggerAnalyzerInfoCollector+;

#endif
#endif
