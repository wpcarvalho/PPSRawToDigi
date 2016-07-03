/****************************************************************************
*  Seyed Mohsen Etesami
****************************************************************************/

#ifndef CondFormats_DiamondReadoutObjects_DiamondFramePosition
#define CondFormats_DiamondReadoutObjects_DiamondFramePosition

#include <iostream>
#include <string>

/**
 * Uniquely identifies the DAQ channel through which a VFAT frame has been received.
 * /verbatim
 * |          32 bits raw position                    |
 * | 14 bits |  10 bits | 4 bits |       4 bits       |
 * | empty   |  FED ID  | GOH ID | index within fiber |
 * \endverbatim
 */
class DiamondFramePosition
{
  public:
    static const unsigned int offsetIdxInFiber = 0, maskIdxInFiber = 0xF;
    static const unsigned int offsetGOHId = 4, maskGOHId = 0xF;
    static const unsigned int offsetFEDId = 8, maskFEDId = 0x3FF;


    DiamondFramePosition(unsigned short FEDId, unsigned short GOHId, unsigned short IdxInFiber) :
      rawPosition(IdxInFiber<<offsetIdxInFiber | GOHId<<offsetGOHId | FEDId<<offsetFEDId)
    {
    }

    /// don't use this constructor unless you have a good reason
    DiamondFramePosition(unsigned int pos = 0) : rawPosition(pos)
    {
    }

    ~DiamondFramePosition()
    {
    }


    unsigned short getFEDId() const
    {
      return (rawPosition >> offsetFEDId) & maskFEDId;
    }

    void setFEDId(unsigned short v)
    { v &= maskFEDId; rawPosition &= 0xFFFFFFFF - (maskFEDId << offsetFEDId); rawPosition |= (v << offsetFEDId); }

    unsigned short getGOHId() const       { return (rawPosition >> offsetGOHId) & maskGOHId; }

    void setGOHId(unsigned short v)
    { v &= maskGOHId; rawPosition &= 0xFFFFFFFF - (maskGOHId << offsetGOHId); rawPosition |= (v << offsetGOHId); }

    unsigned short getIdxInFiber() const  { return (rawPosition >> offsetIdxInFiber) & maskIdxInFiber; }
    
    void setIdxInFiber(unsigned short v)
    { v &= maskIdxInFiber; rawPosition &= 0xFFFFFFFF - (maskIdxInFiber << offsetIdxInFiber); rawPosition |= (v << offsetIdxInFiber); }


    /// don't use this method unless you have a good reason
    unsigned int getRawPosition() const
    {
      return rawPosition;
    }

    bool operator < (const DiamondFramePosition &pos) const
    {
      return (rawPosition < pos.rawPosition);
    }

    bool operator == (const DiamondFramePosition &pos) const
    {
      return (rawPosition == pos.rawPosition);
    }

   
    /// prints 5-digit hex number, the digits correspond to SubSystem, FED ID,GOH ID, index within fiber respectively
    friend std::ostream& operator << (std::ostream& s, const DiamondFramePosition &fp)
    {

    return s<<fp.getFEDId() << ":"
    << fp.getGOHId() << ":"
    << fp.getIdxInFiber();
    }

 
    /// prints XML formatted DAQ channel to stdout
    void printXML();

    /// Sets attribute with XML name 'attribute' and value 'value'.
    /// Also turns on attribute presents bit in the flag parameter
    /// returns 0 if the attribute is known, non-zero value else
    unsigned char setXMLAttribute(const std::string &attribute, const std::string &value, unsigned char &flag);

    /// returns true if all attributes have been set
    static bool checkXMLAttributeFlag(unsigned char flag)
    {
      return (flag == 0x1f);
    }

  protected:
    unsigned int rawPosition;
};

#endif

