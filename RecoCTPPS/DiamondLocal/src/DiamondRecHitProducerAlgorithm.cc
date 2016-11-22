/****************************************************************************
*
* CTPPS diamond timing detector reconstruction 
*
* Based on RecoCTPPS/TotemRPLocal/src/TotemRPRecHitProducerAlgorithm.cc 
*
* Author:
*   Wagner Carvalho (wcarvalh@cern.ch)
*
****************************************************************************/

#include "RecoCTPPS/DiamondLocal/interface/DiamondRecHitProducerAlgorithm.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

//----------------------------------------------------------------------------------------------------

DiamondRecHitProducerAlgorithm::DiamondRecHitProducerAlgorithm(const edm::ParameterSet &param) 
  :param_(param)
{
}

//----------------------------------------------------------------------------------------------------

DiamondRecHitProducerAlgorithm::~DiamondRecHitProducerAlgorithm() 
{
}

//----------------------------------------------------------------------------------------------------

void DiamondRecHitProducerAlgorithm::buildRecHit(unsigned int detId, 
    const std::vector<DiamondDigi> &input, std::vector<DiamondRecHit> &output)
{
  for (std::vector<DiamondDigi>::const_iterator it = input.begin(); it!=input.end(); ++it)
  {
    edm::LogInfo("CTPPS") << " DiamondDigi::chid = " << it->getChannelId() << " , " 
                          << " DiamondDigi::ledgt = " << it->getLeadingEdge() << std::endl ; 
    DetId id = it->getChannelId() ;
    float time = it->getLeadingEdge() ;
    float time_error = 0 ;   // Still to be implemented
    float x = 0 ;         // Still to be implemented
    float x_error = 0 ;   // Still to be implemented
    float y = 0 ;         // Still to be implemented
    float y_error = 0 ;   // Still to be implemented
    uint32_t flag = 0 ;   // Still to be implemented based on HPTDCErrorFlags
    output.push_back(DiamondRecHit(id, time, time_error, x, x_error, y, y_error, flag));
  }  
}
