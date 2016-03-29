/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "RecoTotemRP/RPClusterizer/interface/RPClusterizerAlgorithm.h"

#include <iostream>

//----------------------------------------------------------------------------------------------------

RPClusterizerAlgorithm::RPClusterizerAlgorithm(const edm::ParameterSet &param)
 :param_(param)
{
  verbosity_ = param_.getParameter<int>("Verbosity");
}

//----------------------------------------------------------------------------------------------------

RPClusterizerAlgorithm::~RPClusterizerAlgorithm()
{
}

//----------------------------------------------------------------------------------------------------

int RPClusterizerAlgorithm::BuildClusters(unsigned int detId, const std::vector<TotemRPDigi> &digi, std::vector<TotemRPCluster> &clusters)
{
  clusters.clear();

  strip_digi_set_.clear();
  strip_digi_set_.insert(digi.begin(), digi.end());

  if (strip_digi_set_.size() == 0)
    return 0;
    
  bool iter_beg=true;
  int cluster_beg=-16;
  int cluster_end;
  int prev_strip=-16;
  int cur_strip;
  
  for (TotemRPDigiSet::const_iterator i=strip_digi_set_.begin(); i!=strip_digi_set_.end(); ++i)
  {
    cur_strip=i->GetStripNo();
    bool non_continuity = (cur_strip!=prev_strip+1);
    
    if (iter_beg)
    {
      cluster_beg=cur_strip;
      iter_beg=false;
    }
    else if (non_continuity)
    {
      cluster_end=prev_strip;
      clusters.push_back(TotemRPCluster(detId, (unsigned short)cluster_beg, (unsigned short) cluster_end));
      
      cluster_beg=cur_strip;
    }
    
    prev_strip=cur_strip;
  }
  
  if (!iter_beg)
  {
    cluster_end=prev_strip;
    clusters.push_back(TotemRPCluster(detId, (unsigned short)cluster_beg, (unsigned short) cluster_end));
  }

  if (verbosity_)
  {
    for(unsigned int i=0; i<clusters.size(); ++i)
    {
      std::cout<<"Cluster, DetId="<<clusters[i].DetId()<<" StrBeg="<<clusters[i].StrBeg()
          <<" StrEnd="<<clusters[i].StrEnd()<<" Pos="<<clusters[i].CentreStripPos()<<std::endl;
    }
  }
    
  return clusters.size();
}
