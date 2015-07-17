#ifndef T2DigiVfat_T2DigiVfat_h
#define T2DigiVfat_T2DigiVfat_h



#include <boost/cstdint.hpp>
#include <map>
#include <vector>

class T2DigiVfat{

 public:
  explicit T2DigiVfat (uint32_t cmsswt2detid, unsigned int pos,unsigned short chipID);
  // T2DigiVfat (uint32_t cmsswt2detid,unsigned int pos);
  T2DigiVfat (uint32_t cmsswt2detid,unsigned int pos, unsigned short chipID,const std::vector<unsigned char> activeCh);
  //explicit T2DigiVfat (int channel, int threshold, int pos, float charge);

  T2DigiVfat ();
  virtual ~T2DigiVfat() {};

  //bool operator==(const T2DigiVfat& digi) const;
  //bool operator<(const T2DigiVfat& digi) const;
  //int threshold() const;
  //int bx() const;
  //float charge() const;

  bool IsPad()
  {
    bool flag=false;
    if(pos_ >= 2 && pos_ <= 14)
      flag=true;
    
    return flag;
  };

  bool IsStrip()
  {
    bool flag=false;
    if(pos_ <= 1)
      flag=true;
    if(pos_ <= 16 && pos_ >= 15)
      flag=true;

    return flag;
  };

  void SetFrameStatus(unsigned int FrameStatus) //0:Frame Ok, 1:CRC fails, 2:FootPrint Fail 
    {
      FrameStatus_=FrameStatus;
    };
  
  void print() const;
  
  int GetChannel(unsigned int channel);                       //return 1 or 0 
  int GetThreshold(unsigned int channel);  
  void SetChannel(int channel,unsigned int value);
  //void SetAllChannel(unsigned int value);
  //void SetThreshold(int thr);
  void SetThreshold(int channel, int thr);
  void SetAllThresholds(std::vector<double> thr);
  int GetNoiseChannel(unsigned int channel); 
  void SetNoiseChannel(unsigned int channel, int noise);
  void SetDeadChannels (std::vector<unsigned int> channels);
  void SetDeadChannelsFromEff(std::vector<double> channels);
  bool IsChannelDead(unsigned int channel);


  unsigned int GetVfatPos() const  {return pos_;}
  uint32_t GetVfatDetId(){return cmsswt2detid_;}
  unsigned short GetChipID(){return chipID_;}
  unsigned int GetFrameStatus(){return FrameStatus_;}
  std::map<unsigned int,unsigned int> ChActMap;    //Channel-Active
  std::map<unsigned int,int> ChThrMap;             //Channel-Threshold  
  std::map<unsigned int,unsigned int> DeadChannel; //Associate 0 for working and 1 for dead channels     
  std::map<unsigned int,unsigned int> NoisyChannel;   //Associate a number 0-5 for the noise level     
  
  double Efficiency_;
  
 private:
  unsigned int  pos_; 
  uint32_t cmsswt2detid_;
  unsigned short chipID_;  
  unsigned int FrameStatus_;   //0:Ok, 1:CRC fails, 2:FootPrint Fail 3:Other problem, not setted.
};


#include<iostream>
/*
inline std::ostream & operator<<(std::ostream & o, const T2DigiVfat& digi) 
{
  return o << " " << digi.strip()
	   << " " << digi.threshold()
	   << " " << digi.bx();
}
*/
#endif
