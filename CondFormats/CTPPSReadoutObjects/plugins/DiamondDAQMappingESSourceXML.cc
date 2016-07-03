/****************************************************************************
* Seyed Mohsen Etesami
****************************************************************************/

#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/Framework/interface/ModuleFactory.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESProducer.h"
#include "FWCore/Framework/interface/EventSetupRecordIntervalFinder.h"
#include "FWCore/Framework/interface/ESProducts.h"
#include "FWCore/Framework/interface/SourceFactory.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "CondFormats/DataRecord/interface/DiamondReadoutRcd.h"

#include "CondFormats/CTPPSReadoutObjects/interface/TotemDAQMappingDiamond.h"
#include "CondFormats/CTPPSReadoutObjects/interface/DiamondAnalysisMask.h"
#include "CondFormats/CTPPSReadoutObjects/interface/DiamondFramePosition.h"

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/dom/DOM.hpp>
#include <xercesc/sax/HandlerBase.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>

#include <memory>
#include <sstream>

//#define DEBUG 1

//----------------------------------------------------------------------------------------------------

using namespace std;

/**
 * \brief Loads TotemDAQMapping and DiamondAnalysisMask from two XML files.
 **/
class DiamondDAQMappingESSourceXML: public edm::ESProducer, public edm::EventSetupRecordIntervalFinder
{
public:
  static const std::string tagVFAT;
  static const std::string tagChannel;
  static const std::string tagAnalysisMask;

  /// Common position tags
  static const std::string tagArm;

  /// RP XML tags
  static const std::string tagRPStation;
  static const std::string tagRPPot;
  static const std::string tagRPPlane;


  /// COMMON Chip XML tags
  static const std::string tagChip1;

  DiamondDAQMappingESSourceXML(const edm::ParameterSet &);
  ~DiamondDAQMappingESSourceXML();

  edm::ESProducts< boost::shared_ptr<TotemDAQMappingDiamond>, boost::shared_ptr<DiamondAnalysisMask> > produce( const DiamondReadoutRcd & );

private:
  unsigned int verbosity;

  /// the mapping files
  std::vector<std::string> mappingFileNames;

  /// the mask files
  std::vector<std::string> maskFileNames;

  /// enumeration of XML node types
  enum NodeType { nUnknown, nTop, nArm, nRPStation, nRPPot, nRPPlane, nChip, nChannel };

  /// whether to parse a mapping of a mask XML
  enum ParseType { pMapping, pMask };

  /// parses XML file
  void ParseXML(ParseType, const string &file, const boost::shared_ptr<TotemDAQMappingDiamond>&, const boost::shared_ptr<DiamondAnalysisMask>&);

  /// recursive method to extract RP-related information from the DOM tree
  void ParseTreeDiamond(ParseType, xercesc::DOMNode *, NodeType, unsigned int parentID,
    const boost::shared_ptr<TotemDAQMappingDiamond>&, const boost::shared_ptr<DiamondAnalysisMask>&);


private:
  /// adds the path prefix, if needed
  string CompleteFileName(const string &fn);

  /// returns the top element from an XML file
  xercesc::DOMDocument* GetDOMDocument(string file);

  /// returns true iff the node is of the given name
  bool Test(xercesc::DOMNode *node, const std::string &name)
  {
    return !(name.compare(xercesc::XMLString::transcode(node->getNodeName())));
  }

  /// determines node type
  NodeType GetNodeType(xercesc::DOMNode *);

  /// returns the content of the node
  string GetNodeContent(xercesc::DOMNode *parent)
  {
    return string(xercesc::XMLString::transcode(parent->getTextContent()));
  }

  /// returns the value of the node
  string GetNodeValue(xercesc::DOMNode *node)
  {
    return string(xercesc::XMLString::transcode(node->getNodeValue()));
  }

  DiamondFramePosition ChipFramePosition(xercesc::DOMNode *chipnode);

/// extracts VFAT's DAQ channel from XML attributes
  DiamondFramePosition DiamondchFramePosition(xercesc::DOMNode *chipnode);

  void GetChannels(xercesc::DOMNode *n, std::set<unsigned char> &channels);

  bool RPNode(NodeType type)
  {
    return ((type == nArm)||(type == nRPStation)||(type == nRPPot)||(type == nRPPlane)||(type == nChip));
  }


  bool CommonNode(NodeType type)
  {
    return ((type==nChip)||(type==nArm));
  }

protected:
  /// sets infinite validity of this data
  virtual void setIntervalFor(const edm::eventsetup::EventSetupRecordKey&, const edm::IOVSyncValue&, edm::ValidityInterval&);
};

//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;
using namespace xercesc;

//const string DiamondDAQMappingESSourceXML::tagVFAT="vfat";
//const string DiamondDAQMappingESSourceXML::tagVFAT="diamond_channel";
 
const string DiamondDAQMappingESSourceXML::tagChannel="channel";
const string DiamondDAQMappingESSourceXML::tagAnalysisMask="analysisMask";

// common XML position tags
const string DiamondDAQMappingESSourceXML::tagArm = "arm";

// common XML Chip tags
const string DiamondDAQMappingESSourceXML::tagChip1 = "diamond_channel";


// specific RP XML tags
const string DiamondDAQMappingESSourceXML::tagRPStation = "station";
const string DiamondDAQMappingESSourceXML::tagRPPot = "rp_detector_set";
const string DiamondDAQMappingESSourceXML::tagRPPlane = "rp_plane_diamond";


//----------------------------------------------------------------------------------------------------

DiamondDAQMappingESSourceXML::DiamondDAQMappingESSourceXML(const edm::ParameterSet& conf) :
  verbosity(conf.getUntrackedParameter<unsigned int>("verbosity", 0)),
  mappingFileNames(conf.getUntrackedParameter< vector<string> >("mappingFileNames")),
  maskFileNames(conf.getUntrackedParameter< vector<string> >("maskFileNames"))
{
  setWhatProduced(this);
  findingRecord<DiamondReadoutRcd>();
}

//----------------------------------------------------------------------------------------------------

string DiamondDAQMappingESSourceXML::CompleteFileName(const string &fn)
{
  FileInPath fip(fn);
  return fip.fullPath();
}

//----------------------------------------------------------------------------------------------------

void DiamondDAQMappingESSourceXML::setIntervalFor(const edm::eventsetup::EventSetupRecordKey&,
  const edm::IOVSyncValue& iosv, edm::ValidityInterval& oValidity)
{
  ValidityInterval infinity(iosv.beginOfTime(), iosv.endOfTime());
  oValidity = infinity;
}

//----------------------------------------------------------------------------------------------------

DOMDocument* DiamondDAQMappingESSourceXML::GetDOMDocument(string file)
{
  XercesDOMParser* parser = new XercesDOMParser();
  parser->parse(file.c_str());

  DOMDocument* xmlDoc = parser->getDocument();

  if (!xmlDoc)
    throw cms::Exception("DiamondDAQMappingESSourceXML::GetDOMDocument") << "Cannot parse file `" << file
      << "' (xmlDoc = NULL)." << endl;

  return xmlDoc;
}

//----------------------------------------------------------------------------------------------------

DiamondDAQMappingESSourceXML::~DiamondDAQMappingESSourceXML()
{ 
}

//----------------------------------------------------------------------------------------------------

edm::ESProducts< boost::shared_ptr<TotemDAQMappingDiamond>, boost::shared_ptr<DiamondAnalysisMask> >
  DiamondDAQMappingESSourceXML::produce( const DiamondReadoutRcd & )
{
  boost::shared_ptr<TotemDAQMappingDiamond> mapping(new TotemDAQMappingDiamond());
  boost::shared_ptr<DiamondAnalysisMask> mask(new DiamondAnalysisMask());
  // initialize Xerces
  try
  {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch)
  {
    char* message = XMLString::transcode(toCatch.getMessage());
    throw cms::Exception("DiamondDAQMappingESSourceXML") << "An XMLException caught with message: " << message << ".\n";
    XMLString::release(&message);
  }

  // load mapping files
  for (unsigned int i = 0; i < mappingFileNames.size(); ++i)
    ParseXML(pMapping, CompleteFileName(mappingFileNames[i]), mapping, mask);

  // load mask files
  for (unsigned int i = 0; i < maskFileNames.size(); ++i)
    ParseXML(pMask, CompleteFileName(maskFileNames[i]), mapping, mask);

  // release Xerces
  XMLPlatformUtils::Terminate();

  // commit the products
  return edm::es::products(mapping, mask);
}

//----------------------------------------------------------------------------------------------------

void DiamondDAQMappingESSourceXML::ParseXML(ParseType pType, const string &file,
  const boost::shared_ptr<TotemDAQMappingDiamond> &mapping, const boost::shared_ptr<DiamondAnalysisMask> &mask)
{

  DOMDocument* domDoc = GetDOMDocument(file);

  DOMElement* elementRoot = domDoc->getDocumentElement();

  if (!elementRoot)
    throw cms::Exception("DiamondDAQMappingESSourceXML::ParseMappingXML") << "File `" <<
      file << "' is empty." << endl;

  ParseTreeDiamond(pType, elementRoot, nTop, 0, mapping, mask);
}

//-----------------------------------------------------------------------------------------------------------

void DiamondDAQMappingESSourceXML::ParseTreeDiamond(ParseType pType, xercesc::DOMNode * parent, NodeType parentType,
  unsigned int parentID, const boost::shared_ptr<TotemDAQMappingDiamond>& mapping,
  const boost::shared_ptr<DiamondAnalysisMask>& mask)
{
#ifdef DEBUG
  printf(">> DiamondDAQMappingESSourceXML::ParseTreeDiamond(%s, %u, %u)\n", XMLString::transcode(parent->getNodeName()),
    parentType, parentID);
#endif

  DOMNodeList *children = parent->getChildNodes();

  for (unsigned int i = 0; i < children->getLength(); i++)
  {
    DOMNode *n = children->item(i);
    if (n->getNodeType() != DOMNode::ELEMENT_NODE)
      continue;
  
    NodeType type = GetNodeType(n);
#ifdef DEBUG
    printf("\tname = %s, type = %u\n", XMLString::transcode(n->getNodeName()), type);
#endif

    // structure control
    if (!RPNode(type))
      continue;

    if ((type != parentType + 1)&&parentType != nRPPot)
    {
      if (parentType == nTop && type == nRPPot)
      {
	    LogPrint("DiamondDAQMappingESSourceXML") << ">> DiamondDAQMappingESSourceXML::ParseTreeDiamond > Warning: tag `" << tagRPPot
					<< "' found in global scope, assuming station ID = 12.";
	    parentID = 12;
      } else {
        throw cms::Exception("DiamondDAQMappingESSourceXML") << "Node " << XMLString::transcode(n->getNodeName())
          << " not allowed within " << XMLString::transcode(parent->getNodeName()) << " block.\n";
      }
    }

    // parse tag attributes
    unsigned int id =0,hw_id = 0; 
    bool id_set = false,hw_id_set = false; 
    bool fullMask = false;
    DOMNamedNodeMap* attr = n->getAttributes();

    for (unsigned int j = 0; j < attr->getLength(); j++)
    {
      DOMNode *a = attr->item(j);



      if (!strcmp(XMLString::transcode(a->getNodeName()), "id"))
      {
        sscanf(XMLString::transcode(a->getNodeValue()), "%u", &id);
        id_set = true;
      }

      if (!strcmp(XMLString::transcode(a->getNodeName()), "hw_id"))
      {
        sscanf(XMLString::transcode(a->getNodeValue()), "%x", &hw_id);
        hw_id_set = true;
      }

      if (!strcmp(XMLString::transcode(a->getNodeName()), "full_mask"))
        fullMask = (strcmp(XMLString::transcode(a->getNodeValue()), "no") != 0);
    }

    // content control
    if (!id_set) 
      throw cms::Exception("DiamondDAQMappingESSourceXML::ParseTreeDiamond") << "id not given for element `"
       << XMLString::transcode(n->getNodeName()) << "'" << endl;


      if (!hw_id_set && type == nChip && pType == pMapping)
      throw cms::Exception("TotemDAQMappingESSourceXML::ParseTreeDiamond") << "hw_id not given for element `"
     << XMLString::transcode(n->getNodeName()) << "'" << endl;
 
    if (type == nRPPlane && id > 9)
      throw cms::Exception("DiamondDAQMappingESSourceXML::ParseTreeDiamond") <<
        "Plane IDs range from 0 to 9. id = " << id << " is invalid." << endl;

#ifdef DEBUG
    printf("\tID found: 0x%x\n", id);
#endif






    // store mapping data
    if (pType == pMapping &&type == nChip)
    {

       const DiamondFramePosition &framepos = ChipFramePosition(n);

      DiamondVFATInfo vfatInfo;
      vfatInfo.hwID = hw_id;     
      vfatInfo.symbolicID.subSystem = CTPPSTimingSymbID::Diamond;

      if (type == nChip)
      {
     vfatInfo.symbolicID.symbolicID = parentID * 100 + id;
     
   vfatInfo.type = DiamondVFATInfo::data;
      }


      mapping->insert(framepos, vfatInfo);

      continue;
    }

    // store mask data
    if (pType == pMask && type == nChip)
    {
      CTPPSTimingSymbID symbId;
      symbId.subSystem = CTPPSTimingSymbID::Diamond;
      symbId.symbolicID = parentID * 10 + id;

      DiamondVFATAnalysisMask am;
      am.fullMask = fullMask;
      GetChannels(n, am.maskedChannels);

      mask->insert(symbId, am);

      continue;
    }
    // recursion (deeper in the tree)
    ParseTreeDiamond(pType, n, type,  parentID * 10 + id, mapping, mask);
  }
}

//----------------------------------------------------------------------------------------------------
DiamondFramePosition DiamondDAQMappingESSourceXML::DiamondchFramePosition(xercesc::DOMNode *chipnode)
{
  DiamondFramePosition fp;
  unsigned char attributeFlag = 0;
  unsigned int Values=0;
  std::string values;

  DOMNamedNodeMap* attr = chipnode->getAttributes();
  for (unsigned int j = 0; j < attr->getLength(); j++)
  {
    DOMNode *a = attr->item(j);

    Values=0;
    sscanf(XMLString::transcode(a->getNodeValue()), "%x", &Values);
    std::ostringstream converts;
    converts << Values;
    values = converts.str();
     if(fp.setXMLAttribute(XMLString::transcode(a->getNodeName()),values , attributeFlag)>1)
   {
      throw cms::Exception("DiamondDAQMappingESSourceXML") <<
        "Unrecognized tag `" << XMLString::transcode(a->getNodeName()) <<
        "' or incompatible value `" << XMLString::transcode(a->getNodeValue()) <<
        "'." << endl;
    }
  }

  if (!fp.checkXMLAttributeFlag(attributeFlag))
  {
    throw cms::Exception("DiamondDAQMappingESSourceXML") <<
      "Wrong/incomplete DAQ channel specification (attributeFlag = " << attributeFlag << ")." << endl;
  }

  return fp;
}


//----------------------------------------------------------------------------------------------------

DiamondFramePosition DiamondDAQMappingESSourceXML::ChipFramePosition(xercesc::DOMNode *chipnode)
{
  DiamondFramePosition fp;
  unsigned char attributeFlag = 0;


  DOMNamedNodeMap* attr = chipnode->getAttributes();
  for (unsigned int j = 0; j < attr->getLength(); j++)
  {
    DOMNode *a = attr->item(j);
   

 if (fp.setXMLAttribute(XMLString::transcode(a->getNodeName()), XMLString::transcode(a->getNodeValue()), attributeFlag) > 1)
   {
      throw cms::Exception("DiamondDAQMappingESSourceXML") <<
        "Unrecognized tag `" << XMLString::transcode(a->getNodeName()) <<
        "' or incompatible value `" << XMLString::transcode(a->getNodeValue()) <<
        "'." << endl;
    }
  }

  if (!fp.checkXMLAttributeFlag(attributeFlag))
  {
    throw cms::Exception("DiamondDAQMappingESSourceXML") <<
      "Wrong/incomplete DAQ channel specification (attributeFlag = " << attributeFlag << ")." << endl;
  }

  return fp;
}

//----------------------------------------------------------------------------------------------------

DiamondDAQMappingESSourceXML::NodeType DiamondDAQMappingESSourceXML::GetNodeType(xercesc::DOMNode *n)
{
  // common node types
  if (Test(n, tagArm)) return nArm;
  if (Test(n, tagChip1)) return nChip;

  // RP node types
  if (Test(n, tagRPStation)) return nRPStation;
  if (Test(n, tagRPPot)) return nRPPot;
  if (Test(n, tagRPPlane)) return nRPPlane;



  throw cms::Exception("DiamondAQMappingESSourceXML::GetNodeType") << "Unknown tag `"
    << XMLString::transcode(n->getNodeName()) << "'.\n";
}

//----------------------------------------------------------------------------------------------------

void DiamondDAQMappingESSourceXML::GetChannels(xercesc::DOMNode *n, set<unsigned char> &channels)
{
  DOMNodeList *children = n->getChildNodes();
  for (unsigned int i = 0; i < children->getLength(); i++)
  {
    DOMNode *n = children->item(i);
    if (n->getNodeType() != DOMNode::ELEMENT_NODE || !Test(n, "channel"))
      continue;

    DOMNamedNodeMap* attr = n->getAttributes();
    bool idSet = false;
    for (unsigned int j = 0; j < attr->getLength(); j++)
    {
      DOMNode *a = attr->item(j);

      if (!strcmp(XMLString::transcode(a->getNodeName()), "id"))
      {
        unsigned int id = 0;
        sscanf(XMLString::transcode(a->getNodeValue()), "%u", &id);
        channels.insert(id);
        idSet = true;
        break;
      }
    }

    if (!idSet)
    {
      throw cms::Exception("DiamondDAQMappingESSourceXML::GetChannels") <<
        "Channel tags must have an `id' attribute.";
    }
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_EVENTSETUP_SOURCE(DiamondDAQMappingESSourceXML);
