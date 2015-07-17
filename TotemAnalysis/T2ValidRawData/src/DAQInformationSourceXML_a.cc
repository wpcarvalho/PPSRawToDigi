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
//#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "TotemAnalysis/T2ValidRawData/interface/DAQInformationSourceXML_a.h"



//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;
using namespace xercesc;

// common XML position tags
const string DAQInformationSourceXML_a::tagArm = "arm";

// common XML Chip tags
const string DAQInformationSourceXML_a::tagChip1 = "vfat";
const string DAQInformationSourceXML_a::tagChip2 = "test_vfat";
const string DAQInformationSourceXML_a::tagTriggerVFAT1 = "trigger_vfat";

// specific RP XML tags
const string DAQInformationSourceXML_a::tagRPStation = "station";
const string DAQInformationSourceXML_a::tagRPPot = "rp_detector_set";
const string DAQInformationSourceXML_a::tagRPPlane = "rp_plane";

// specific T2 XML tags
const string DAQInformationSourceXML_a::tagT2="t2_detector_set";
const string DAQInformationSourceXML_a::tagT2detector="t2_detector";
const string DAQInformationSourceXML_a::tagT2Half="t2_half";

// specific T1 XML tags
const string DAQInformationSourceXML_a::tagT1="t1_detector_set";
const string DAQInformationSourceXML_a::tagT1Arm="t1_arm";
const string DAQInformationSourceXML_a::tagT1Plane="t1_plane";
const string DAQInformationSourceXML_a::tagT1CSC="t1_csc";

//----------------------------------------------------------------------------------------------------

DAQInformationSourceXML_a::DAQInformationSourceXML_a(std::string nameXMLfile) : xmlFileName("")
{
  xmlFileName = nameXMLfile;
}




//----------------------------------------------------------------------------------------------------

DAQInformationSourceXML_a::~DAQInformationSourceXML_a()
{ 
}

//----------------------------------------------------------------------------------------------------

boost::shared_ptr<DAQInformationT2_a> DAQInformationSourceXML_a::produceMap()
{
   
   boost::shared_ptr<DAQInformationT2_a> dataT2 (new DAQInformationT2_a());
  
  if (xmlFileName.empty()) { 
    throw cms::Exception("DAQInformationSourceXML_a::produce") << "Xml file not found" << endl;
  } else {
    ParseXML(xmlFileName,  dataT2);
  }

  return dataT2;
}

//----------------------------------------------------------------------------------------------------


void DAQInformationSourceXML_a::ParseXML(std::string filename, const  boost::shared_ptr<DAQInformationT2_a> &dataT2)
{
#ifdef DEBUG
  printf(">> DAQInformationSourceXML_a::ParseXML(%s)\n", filename.c_str());
#endif

  // load DOM tree from the file
  try {
    XMLPlatformUtils::Initialize();
  }
  catch (const XMLException& toCatch) {
    char* message = XMLString::transcode(toCatch.getMessage());
    throw cms::Exception("DAQInformationSourceXML_a") << "An XMLException caught with message: " << message << ".\n";
    XMLString::release(&message);
  }

  XercesDOMParser* parser = new XercesDOMParser();
  parser->setValidationScheme(XercesDOMParser::Val_Always);
  parser->setDoNamespaces(true);

  try {
    parser->parse(filename.c_str());
  }
  catch (...) {
    throw cms::Exception("DAQInformationSourceXML_a") << "Cannot parse file 1`" << filename << "'." << endl;
  }

  if (!parser)
    throw cms::Exception("DAQInformationSourceXML_a") << "Cannot parse file 2`" << filename << "'." << endl;
  
  DOMDocument* xmlDoc = parser->getDocument();

  if (!xmlDoc)
    throw cms::Exception("DAQInformationSourceXML_a") << "Cannot parse file 3`" << filename << "'." << endl;

 
  DOMElement* elementRootT2 = xmlDoc->getDocumentElement();
  

  //if (!elementRoot)
  //throw cms::Exception("DAQInformationSourceXML_a") << "File `" << filename << "' is empty." << endl;
  
  // extract useful information form the DOM tree
 
  //dataT2->Reset();
 
  ParseTreeT2(elementRootT2, nTop, 0, dataT2);
 
  XMLPlatformUtils::Terminate();
}

//----------------------------------------------------------------------------------------------------

FramePosition DAQInformationSourceXML_a::ChipFramePosition(xercesc::DOMNode *chipnode)
{
  FramePosition fp;
  unsigned char attributeFlag = 0;

  DOMNamedNodeMap* attr = chipnode->getAttributes();
  for (unsigned int j = 0; j < attr->getLength(); j++) {    
    DOMNode *a = attr->item(j);
    if (fp.SetXMLAttribute(XMLString::transcode(a->getNodeName()), XMLString::transcode(a->getNodeValue()), attributeFlag) > 1)
      throw cms::Exception("DAQInformationSourceXML_a") <<
        "Unrecognized tag `" << XMLString::transcode(a->getNodeName()) <<
        "' or incompatible value `" << XMLString::transcode(a->getNodeValue()) <<
        "'." << endl;
  }

  if (!fp.CheckXMLAttributeFlag(attributeFlag))
    throw cms::Exception("DAQInformationSourceXML_a") <<
      "Wrong/incomplete DAQ channel specification (attributeFlag = " << attributeFlag << ")." << endl;

  return fp;
}

//----------------------------------------------------------------------------------------------------

bool DAQInformationSourceXML_a::RPNode(nodeType type)
{
  return ((type == nArm)||(type == nRPStation)||(type == nRPPot)||(type == nRPPlane)||(type == nChip)||(type == nTriggerVFAT));
}

//----------------------------------------------------------------------------------------------------

bool DAQInformationSourceXML_a::T2Node(nodeType type)
{
  return ((type==nT2)||(type==nT2Det)|| (type==nT2Half));
}

//----------------------------------------------------------------------------------------------------

bool DAQInformationSourceXML_a::T1Node(nodeType type)
{
  return ((type==nT1)||(type==nT1Arm)|| (type==nT1Plane) || (type==nT1CSC));
}

//----------------------------------------------------------------------------------------------------

bool DAQInformationSourceXML_a::CommonNode(nodeType type)
{
  return ((type==nChip)||(type==nArm));
}

//----------------------------------------------------------------------------------------------------

bool DAQInformationSourceXML_a::Test(DOMNode *node, const std::string &name)
{
  return !(name.compare(XMLString::transcode(node->getNodeName())));
}

//----------------------------------------------------------------------------------------------------

DAQInformationSourceXML_a::nodeType DAQInformationSourceXML_a::GetNodeType(xercesc::DOMNode *n) {
 
  // common node types
  if (Test(n, tagArm)) return nArm;
  if (Test(n, tagChip1)) return nChip;
  if (Test(n, tagChip2)) return nChip;
  if (Test(n, tagTriggerVFAT1)) return nTriggerVFAT;

  // RP node types
  if (Test(n, tagRPStation)) return nRPStation;
  if (Test(n, tagRPPot)) return nRPPot;
  if (Test(n, tagRPPlane)) return nRPPlane;

  // T2 node types
  if (Test(n, tagT2)) return nT2;
  if (Test(n, tagT2detector)) return nT2Det;  
  if (Test(n, tagT2Half)) return nT2Half;

  // T1 node types
  if (Test(n, tagT1)) return nT1;
  if (Test(n, tagT1Arm)) return nT1Arm;  
  if (Test(n, tagT1Plane)) return nT1Plane;
  if (Test(n, tagT1CSC)) return nT1CSC;
  

  throw cms::Exception("DAQInformationSourceXML_a") << "Unknown tag `" << XMLString::transcode(n->getNodeName()) << "'.\n";
}



//----------------------------------------------------------------------------------------------------

void DAQInformationSourceXML_a::ParseTreeT2(DOMNode *parent, nodeType parentType, unsigned int parentID, const  boost::shared_ptr<DAQInformationT2_a> &dataT2) {

  DOMNodeList *children = parent->getChildNodes();

#ifdef DEBUG
  printf(">> ParseTreeT2(parent,parentType,parentID)=(%p, %i, %u)\n", parent, parentType, parentID);
  printf("\tchildren: Numero children: %li\n", children->getLength());
#endif

  
  for (unsigned int i = 0; i < children->getLength(); i++) {
    DOMNode *n = children->item(i);

    if (n->getNodeType() != DOMNode::ELEMENT_NODE) continue; //bypass iteration 

    // get node type for RP or T2
    nodeType type = GetNodeType(n);
   
#ifdef DEBUG
    printf("\t\tchildren #%i: is a %s, (of type %i) \n", i, XMLString::transcode(n->getNodeName()), type);
#endif

    if ((type == nUnknown)) 
      {
	continue;
#ifdef DEBUG
	printf("Found Unknown tag during T2 reading.. EXIT ");
#endif 	
      }
    if((T2Node(type)==false)&&(CommonNode(type)==false)) 
      {
#ifdef DEBUG
	printf("Found Non-T2 tag during T2 reading.. EXIT ");
	printf("\t The tag is:  %s \n", XMLString::transcode(n->getNodeName()));
#endif
	continue;
      }
     
    // get ID_t2 and position
    unsigned int ID_t2 = 0;
    //id  for T2 plane goes from 0..9; for chip is the 16 bit ID
    unsigned int position_t2 = 0; 
    //position_t2 was the S-link for chip and for the plane should be a number compatible with arm,ht,pl,pls or HS position

    unsigned int arm=0,ht=0,pl=0,pls=0;
    unsigned int VFPOS0_16=0; 
      
    bool idSet_t2 = false;      
    int attribcounter_t2planedescript=0;
    unsigned int toaddForParentID=0;

    DOMNamedNodeMap* attr = n->getAttributes();
    
    //    Begin loop for save T2 element attriute  ------------------------------------------------------------------
    for (unsigned int j = 0; j < attr->getLength(); j++) {
	  
      DOMNode *a = attr->item(j);
      if (!strcmp(XMLString::transcode(a->getNodeName()), "id")) {
	sscanf(XMLString::transcode(a->getNodeValue()), "%u", &ID_t2);
	idSet_t2 = true;
      }
      
      if (!strcmp(XMLString::transcode(a->getNodeName()), "position")) {
	position_t2 = atoi(XMLString::transcode(a->getNodeValue()));
      }
      
      if(type == nArm) {
	if (!strcmp(XMLString::transcode(a->getNodeName()), "id")) {
	  //arm is the top node and should be reset to 0.
	  parentID=0;
	  unsigned int id_arm = atoi(XMLString::transcode(a->getNodeValue()));
	  toaddForParentID=20*id_arm;
	}
      }
      
      if(type == nT2Half) {
	if (!strcmp(XMLString::transcode(a->getNodeName()), "id")) {
	  unsigned int id_half = atoi(XMLString::transcode(a->getNodeValue()));	
	  toaddForParentID=10*id_half;	
	}
      }
      
      //This is needed in principle only for the old formats
      if(type == nT2Det) {
	if (!strcmp(XMLString::transcode(a->getNodeName()), "arm")) {
	  sscanf(XMLString::transcode(a->getNodeValue()), "%u", &arm);
	  attribcounter_t2planedescript++;
	}
	
	if (!strcmp(XMLString::transcode(a->getNodeName()), "ht")) {
	  sscanf(XMLString::transcode(a->getNodeValue()), "%u", &ht);
	  attribcounter_t2planedescript++;
	}
	
	if (!strcmp(XMLString::transcode(a->getNodeName()), "pl")) {
	  sscanf(XMLString::transcode(a->getNodeValue()), "%u", &pl);
	  attribcounter_t2planedescript++;
	}
	
	if (!strcmp(XMLString::transcode(a->getNodeName()), "pls")) {
	  sscanf(XMLString::transcode(a->getNodeValue()), "%u", &pls);
	  attribcounter_t2planedescript++;
	}
	
	//remember id in monitor goes from 0 -- 39
	if (!strcmp(XMLString::transcode(a->getNodeName()), "id")) {
	  //Id saved another time ... just to increment attribcounter
	  sscanf(XMLString::transcode(a->getNodeValue()), "%u", &ID_t2);
	  attribcounter_t2planedescript++;
	}
	
	if (!strcmp(XMLString::transcode(a->getNodeName()), "position")) {		  
	  sscanf(XMLString::transcode(a->getNodeValue()), "%u", &position_t2);
	  attribcounter_t2planedescript++;
	  //Just another indication for further checking. This attribute was not compulsory in monitor.
	  attribcounter_t2planedescript=attribcounter_t2planedescript+20;
	  //20 is just a "big number"
	}
      }
      
      if(type == nChip)
	if (!strcmp(XMLString::transcode(a->getNodeName()), "iid")) 
	  sscanf(XMLString::transcode(a->getNodeValue()), "%u", &VFPOS0_16);
      
    }  
    // End loop for save T2 element attriute  ------------------------------------------------------------------
    
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // When a plane tag is found, calculate Parent-Id and allows different xml formats:
    // Provide compatibility with old monitor formats.
    
    //Note:
    //plane position and id foreseen in  final ip5 mapping.
    //plane position NOT foreseen in Monitor.
    
    if(type == nT2Det) {
      //Calculate the parent-id from attributes or xml STRUCTURE
      
      if(attribcounter_t2planedescript>=(21)) {
	
	//there is already in xml plane the  "position" attribute + other attributes. It is assumed to utilize the parent info
	
#ifdef DEBUG  
	printf("DAQInformationSourceXML_a: attribcounter_t2planedescript: %i \n",attribcounter_t2planedescript);
#endif
	
	if(attribcounter_t2planedescript>=(25)) {
	  
	  edm::LogVerbatim("DAQInformationSourceXML_a")<<"T2-Plane attribute utilezed for parentID: position+info from parent ";
	  //Plane Seems fully specified
	  //all T2 plane attribute read correctly. Check if it is consitent
	  unsigned int test_position_t2=arm*20+ht*10+pl*2+pls;
	  unsigned int testHS=pl*2+pls;
	  if(testHS!=position_t2) {
	    edm::LogPrint("DAQInformationSourceXML_a") <<"T2 Xml inconsistence in pl-pls attributes and position. Only 'position attribute' taken ";
	    testHS=position_t2;		    
	  }
	  
	  // For plane, ID_t2 should go from 0..39 position_t2 from 0..9
	  ID_t2=parentID+position_t2;
	  toaddForParentID=position_t2;
	  if(ID_t2!=test_position_t2)
	    edm::LogPrint("DAQInformationSourceXML_a") <<"T2 Xml inconsistence in plane attributes and xml parents structure. Plane attributes ignored";
	  
	}
	else { 
	  //Case where arm-ht-pl-pls are NOT specified	      
	  edm::LogVerbatim("DAQInformationSourceXML_a")<<"T2 Plane have parentID: "<<parentID<<" for its VFATs. Plane Position read: "<<position_t2;
	  
	  if(attribcounter_t2planedescript==21) {
	    //You have put in XML only position and not Plane id (the last is obligatory)
	    
	    ID_t2=parentID+position_t2;
	    toaddForParentID=position_t2;
	    idSet_t2=true;
	  }
	} 
      }   
      else {
	//Construct plane position from other attributes cause "position" is not inserted;
	//Ex- monitor:    <t2_detector id="0" z="13871.3" arm="0" ht="0" pl="0" pls="0" >
	
	if(attribcounter_t2planedescript>=1) {
	  //Remember, Z attribute is not counted
	    
	  if(attribcounter_t2planedescript>=5) {
	    int test_position_t2=arm*20+ht*10+pl*2+pls;
	    
	    //case for xml from monitor		   
	    ID_t2=test_position_t2;
	    toaddForParentID=test_position_t2;
	    
	    if (parentID!=ID_t2) {
	      edm::LogPrint("DAQInformationSourceXML_a") <<"T2 Inconsistence between plane 'id' and position from attributes. Id ignored";
	      edm::LogPrint("DAQInformationSourceXML_a") <<" T2-Parent = "<<parentID;		 			  
	    }		  				    
	  }
	  else {		   
	    toaddForParentID=ID_t2;		
	    edm::LogVerbatim("DAQInformationSourceXML_a")<<" Number of T2 plane attributes: "<< attribcounter_t2planedescript<<" T2-Plane attribute utilezed for parentID: plane 'id' only";
	  }
	} 
	else {
	  edm::LogProblem ("DAQInformationSourceXML_a") << "T2 plane not enough specified from its attribute!";	      	      
	}
      }	   
    }
    
    // ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    // content control 
    // Note: each element has an id!! However if the element is a plane, it could be enough to use position 0..9
    
    if(idSet_t2==false){
      throw cms::Exception("DAQInformationSourceXML_a::ParseTree") << "ID_t2 not given for element `" << XMLString::transcode(n->getNodeName()) << "'" << endl;
      edm::LogProblem ("DAQInformationSourceXML_a") <<"ID_t2 not given for element `"<<XMLString::transcode(n->getNodeName()) << "'";
    }
    
    if (type == nT2Det && position_t2 > 39) {
      throw cms::Exception("DAQInformationSourceXML_a::ParseTree") << "Plane position_t2 range from 0 to 39. position_t2 = " << position_t2 << " is invalid." << endl;
      edm::LogProblem ("DAQInformationSourceXML_a") <<"Plane position_t2 range from 0 to 39. position_t2 = "<<position_t2<< " is invalid.";
      
    }
    // if (!position_t2Set && type == nChip) {
    // throw cms::Exception("DAQInformationSourceXML_a::ParseTree") << "Position not given for T2 chip in detector " << parentID << endl;
    // edm::LogProblem ("DAQInformationSourceXML_a") <<"Position not given for T2 chip in parent element-id: "<< parentID;
    // }
    
    if (type == nChip) {
      // save information      
#ifdef DEBUG
      printf("T2 Vfat in plane (parentID): %i || GeomPosition %i \n",parentID,VFPOS0_16);
      printf("\t\t\tID_t2 = 0x%x\n", ID_t2);
      printf("\t\t\tpos = %i\n", position_t2);
#endif     
      unsigned int symId=0; 
      //Check if it is a special chip
      if(!tagT2detector.compare(XMLString::transcode((n->getParentNode()->getNodeName()))))
	symId = parentID * 100 + VFPOS0_16; // same conv = symbplaneNumber*100 +iid used in DQM 
      else {
	//It is a special VFAT and the special number is set directly in the XML file
	symId=VFPOS0_16;      //17,18,19,20
#ifdef DEBUG
	printf("DAQInformationSourceXML_a Found T2 special Vfat ChId-SLink-Symb  0x%x - %i - %i \n",ID_t2,position_t2,symId );
#endif
      }
      
      FramePosition framepos = ChipFramePosition(n);
      dataT2->readoutPositionToId[framepos] = symId;      
      //Assign a contanaier for the register of that VFAT with ChipId (Hex)
      dataT2->readoutIdToRegisters[symId] = VFATRegisters(ID_t2);   
    } 
    else {
      //Look for the children of n (recursion)
      //3° argument=parentId  is needed for calculate VFAT-id startintg from the parent plane 
      ParseTreeT2(n, type, parentID+toaddForParentID, dataT2);
    }
            
  }//Go to the next children

}








//DEFINE_ANOTHER_FWK_EVENTSETUP_SOURCE(DAQInformationSourceXML_a);

