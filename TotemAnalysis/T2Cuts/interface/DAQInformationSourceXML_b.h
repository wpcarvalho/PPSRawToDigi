/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*   Mirko Berretti (mirko.berretti@gmail.com) 
*
* $Id: Element.cc 2133 2010-02-16 13:08:48Z jkaspar $
* $Revision: 1.1.2.7 $
* $Date: 2009/11/28 17:14:04 $
*
* IMPORTANT:
* Please keep the documentation at
* https://twiki.cern.ch/twiki/bin/view/TOTEM/CompTotemRawData
* uptodate.
*
****************************************************************************/

//#include "FWCore/Framework/interface/MakerMacros.h"
//#include "FWCore/Framework/interface/SourceFactory.h"
//#include "FWCore/Framework/interface/ModuleFactory.h"
//#include "FWCore/Framework/interface/ESHandle.h"
//#include "FWCore/ParameterSet/interface/ParameterSet.h"
//#include "FWCore/Framework/interface/ESProducer.h"
//#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
//#include "FWCore/Framework/interface/ESProducts.h"
//#include "FWCore/Framework/interface/SourceFactory.h"
#ifndef DAQInformationSourceB_
#define DAQInformationSourceB_

#include "FWCore/MessageLogger/interface/MessageLogger.h"


#include "TotemRawDataLibrary/DataFormats/interface/FramePosition.h"
#include "TotemAnalysis/T2Cuts/interface/DAQInformationT2_b.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQRecord.h"

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <memory>
#include <boost/shared_ptr.hpp>

//#define DEBUG 1


/**
 * \brief Provides   DAQInformationT2_b 
 *
 * Two input sources are supported: directly from cfg file and from Monitor XML.
 *
 * TODO: support coincidence chips for T2
 *
 * T1 information as well
 **/
class DAQInformationSourceXML_b
{
  public:
    /// Common position tags
    static const std::string tagArm;
    
    /// RP XML tags
    static const std::string tagRPStation;
    static const std::string tagRPPot;
    static const std::string tagRPPlane;
    
    /// T2 XML tags
    static const std::string tagT2;
    static const std::string tagT2Half;
    static const std::string tagT2detector;
    
    /// T1 XML tags
    static const std::string tagT1;
    static const std::string tagT1Arm;
    static const std::string tagT1Plane;
    static const std::string tagT1CSC;

    /// COMMON Chip XML tags 
    static const std::string tagChip1;
    static const std::string tagChip2;
    static const std::string tagTriggerVFAT1;

  DAQInformationSourceXML_b(std::string nameXMLfile);
  ~DAQInformationSourceXML_b();
   boost::shared_ptr<DAQInformationT2_b> produceMap();

  
  private:
    unsigned int verbosity;
    std::string xmlFileName;

    /// node types
    enum nodeType { nUnknown, nTop, nArm, nRPStation, nRPPot, nRPPlane, nChip, nTriggerVFAT,
      nT2, nT2Half, nT2Det, nT1, nT1Arm, nT1Plane, nT1CSC};                    
        
    /// parses XML file
    void ParseXML(std::string xmlfilename, const boost::shared_ptr<DAQInformationT2_b>&);     
  
   
    /// recursive method to extract T2-related information from the DOM tree
    void ParseTreeT2(xercesc::DOMNode *, nodeType type, unsigned int ID, const boost::shared_ptr<DAQInformationT2_b>&); 

   
    
    /// returns true iff the node is of the given name
    bool Test(xercesc::DOMNode *node, const std::string &name);

    /// determines node type
    nodeType GetNodeType(xercesc::DOMNode *);

    /// extracts VFAT's DAQ channel from XML attributes
    FramePosition ChipFramePosition(xercesc::DOMNode *chipnode);

    /// returns true if the node type is accepted by ParseTreeRP
    bool RPNode(nodeType type);        

    bool T2Node(nodeType type);        
  
    bool T1Node(nodeType type);
    
    bool CommonNode(nodeType type);

    // protected:
    /// sets infinite validity of this information
    //virtual void setIntervalFor(const edm::eventsetup::EventSetupRecordKey&, const edm::IOVSyncValue&, edm::ValidityInterval&);
};




#endif

//DEFINE_ANOTHER_FWK_EVENTSETUP_SOURCE(DAQInformationSourceXML_b);

