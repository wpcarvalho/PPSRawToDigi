

#include <DataFormats/T2DigiVfat/interface/T2DigiVfat.h>
using namespace std;

T2DigiVfat::T2DigiVfat (uint32_t cmsswt2detid, unsigned int pos,unsigned short chipID):
  pos_(pos),
  cmsswt2detid_(cmsswt2detid),chipID_(chipID)
{
  Efficiency_=0.8;
}
//  Efficiency_=0.8;

//T2DigiVfat::T2DigiVfat(uint32_t cmsswt2detid,unsigned int pos):pos_(pos),cmsswt2detid_(cmsswt2detid)
//{}



T2DigiVfat::T2DigiVfat(uint32_t cmsswt2detid,unsigned int pos,unsigned short chipID,const vector<unsigned char> activeCh):pos_(pos),cmsswt2detid_(cmsswt2detid),chipID_(chipID)
{
  Efficiency_=0.8;

  for(unsigned int i=0;i<activeCh.size();i++)
    {
      if((int)activeCh.at(i)>128)
	std::cout<<"Warning in T2DigiVfat- Vfat channel>128"<<std::endl;
      else
	{	  
	    SetChannel((int)activeCh.at(i),1);	  
	}
    }

  //Put 0 to the unactive channels

  std::map<unsigned int,unsigned int>::const_iterator it;
  //Look if it is dead
  std::map<unsigned int,unsigned int>::const_iterator iter;
	
  for (unsigned int i=0;i<128;i++)
    { 
      iter=DeadChannel.find(i);  
      it=ChActMap.find(i); 
      if((it==ChActMap.end())||(iter!=DeadChannel.end())) 
	{
	  SetChannel(i,0);
	}
    }
}




T2DigiVfat::T2DigiVfat ()
{}


int T2DigiVfat::GetChannel(unsigned int channel)                       //return 1 or 0 
{
  int out=-1;
  //if(find(ChActMap.begin(), ChActMap.end(), channel) != ChActMap.end())

  std::map<unsigned int,unsigned int>::const_iterator it=ChActMap.find(channel); 
  if(it!=ChActMap.end())
    {
      out=ChActMap[channel];
    }
  return out;
}

int T2DigiVfat::GetThreshold(unsigned int channel)  
{
  int out=-1;
   std::map<unsigned int,int>::const_iterator it=ChThrMap.find(channel); 
   if(it!=ChThrMap.end())
     {
       out=ChThrMap[channel];
     }
   return out;
}

void T2DigiVfat::SetChannel(int channel,unsigned int value)
{
  std::map<unsigned int,unsigned int>::const_iterator test;
  
  if(channel==-1)     //set all channels
    {      
      for(unsigned int k=0;k<128;k++)
	{
	  test=ChActMap.find(k);
	  if(test==ChActMap.end())
	    {
	      ChActMap.insert(std::pair<unsigned int,int>(k,value));
	    }
	  else{
	    ChActMap[k]=value;
	  }
	}
    }
  else
    {
      std::map<unsigned int,unsigned int>::const_iterator itd=DeadChannel.find(channel);
      unsigned int val=value;
      if(itd!=DeadChannel.end())
	val=0;

      std::map<unsigned int,unsigned int>::const_iterator it=ChActMap.find(channel);
      if(it==ChActMap.end()){
	ChActMap.insert(std::pair<unsigned int,unsigned int>(channel,val));    //Channel-Active
	//ChThrMap.insert(std::pair<unsigned int,int>(channel,thr));
      }
      else
	ChActMap[channel]=val;
    }

    
}


bool T2DigiVfat::IsChannelDead(unsigned int channel)
{
  bool retval=false;
  if(DeadChannel.find(channel)!=DeadChannel.end())
    retval=true;

  return retval;

}

void T2DigiVfat::SetAllThresholds(std::vector<double> thr)
{
   for(unsigned int k=0;k<thr.size();k++)
     {
       if(k>128)
	 std::cout<<"T2DigiVfat warning!!! channel number exceed 128!!"<<std::endl;
       else
	 ChThrMap.insert(std::pair<unsigned int,int>(k,(int)thr.at(k)));
     }
}


void T2DigiVfat::SetThreshold(int channel, int thr)
{
  std::map<unsigned int,int>::const_iterator test;
  if(channel==-1)     //set all channels
    {      
      for(unsigned int k=0;k<128;k++)
	{
	  test=ChThrMap.find(k);
	  if(test==ChThrMap.end())
	    {
	      ChThrMap.insert(std::pair<unsigned int,int>(k,thr));
	    }
	  else{
	    ChThrMap[k]=thr;
	  }
	}
    }
  else
    {
      std::map<unsigned int,int>::const_iterator it=ChThrMap.find(channel);
      if(it==ChThrMap.end())
	ChThrMap.insert(std::pair<unsigned int,int>(channel,thr));
      else
	ChThrMap[channel]=thr;
    }
}


void T2DigiVfat::SetDeadChannelsFromEff(std::vector<double> channels)
{
  for(unsigned int k=0;k<channels.size();k++)
    {
      DeadChannel.insert(std::pair<unsigned int,unsigned int>((unsigned int)channels.at(k),1));
    }
}

void T2DigiVfat::SetDeadChannels(std::vector<unsigned int> channels)
{
  for(unsigned int k=0;k<channels.size();k++)
    {
      DeadChannel.insert(std::pair<unsigned int,unsigned int>(channels.at(k),1));
    }
}

void T2DigiVfat::SetNoiseChannel(unsigned int channel, int noiselev)
{
  NoisyChannel.insert(std::pair<unsigned int,int>(channel,noiselev));
}

int T2DigiVfat::GetNoiseChannel(unsigned int channel)
{
  int out=-1;
   std::map<unsigned int,unsigned int>::const_iterator it=NoisyChannel.find(channel); 
   if(it!=NoisyChannel.end())
     {
       out=NoisyChannel[channel];
     }
   return out;
}


void T2DigiVfat::print() const {
  std::cout << "POS " << pos_ 
	    << " CMSSWID " << cmsswt2detid_ <<std::endl;
}


