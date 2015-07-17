/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $Id: TotemRawMetaDataDumper.cc 3290 2010-09-07 11:54:29Z jkaspar $
* $Revision: 2145 $
* $Date: 2010-02-17 15:01:22 +0100 (Wed, 17 Feb 2010) $
*
****************************************************************************/

#ifndef _TotemAnalysis_TotemNtuplizer_Ntuplizer_h_
#define _TotemAnalysis_TotemNtuplizer_Ntuplizer_h_

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

class TTree;

/**
 * A template for ntuple-making classes
 **/
class Ntuplizer
{
 public:
  /// constructor
  Ntuplizer(const edm::ParameterSet&);

  /// book your branches here
  virtual void CreateBranches(const edm::EventSetup&, TTree *) = 0;

  /// fill the event data here
  virtual void FillEvent(const edm::Event&, const edm::EventSetup&) = 0;
};

#endif

