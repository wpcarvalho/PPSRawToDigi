/****************************************************************************
*  Seyed Mohsen Etesami
****************************************************************************/

#include "CondFormats/CTPPSReadoutObjects/interface/DiamondFramePosition.h"

#include <iomanip>
#include <cstdlib>

using namespace std;

//----------------------------------------------------------------------------------------------------

void DiamondFramePosition::printXML()
{
   cout << "\" FEDId=\"" << getFEDId()
    << "\" GOHId=\"" << getGOHId()
    << "\" IdxInFiber=\"" << getIdxInFiber()
    << "\"";
}

//----------------------------------------------------------------------------------------------------

unsigned char DiamondFramePosition::setXMLAttribute(const std::string &attribute, const std::string &value,
    unsigned char &flag)
{
  unsigned int v = atoi(value.c_str());

  if (attribute == "FEDId")
  {
    setFEDId(v);
    flag |= 0x1c; 
    return 0;
  }


  if (attribute == "GOHId")
  {
    setGOHId(v);
    flag |= 0x2;
    return 0;
  }

  if (attribute == "IdxInFiber")
  {
    setIdxInFiber(v);
    flag |= 0x1;
    return 0;
  }

  return 1;
}
