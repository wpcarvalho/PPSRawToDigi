#define IOMC_DPEProtons_h101_cxx

#include "IOMC/DPEProtons/interface/h101.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <iostream>


bool h101::GetForwardProton(int direction, MADProton &proton)
{
  bool found = false;
  for(int i=0; i<Nhep; i++)
  {
    int ISTHEP = Jsdhep[i]/16000000*100 + Jsmhep[i]/16000000;
    if(ISTHEP==1 && Idhep[i]==2212 && Jsdhep[i]==0 && Phep[i][3]>3500 && Phep[i][2]*direction>0)
    {
      proton.x = 0;
      proton.y = 0;
      proton.px = Phep[i][0];
      proton.py = Phep[i][1];
      proton.pz = Phep[i][2];
      found = true;
      break;
    }
  }
  
  return found;
}


Long64_t h101::GetEntries()
{
  if(fChain)
    return fChain->GetEntries();
  else return 0;
}

MADProtonPair h101::GetProtonPair(Long64_t jentry)
{
  MADProtonPair pair;
  
  Long64_t ientry = LoadTree(jentry);
  if (ientry < 0)
  {
    std::cout<<"DPEProton : fatal input problem, exit(0)"<<std::cout;
  }
  
  fChain->GetEntry(jentry);
  
  GetForwardProton(-1, pair.l);
  GetForwardProton(1, pair.r);
  
  return pair;
}

