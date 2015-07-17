// Description: Function to retrieve Generated Kinematic et al.

// DataFormats
#include "DataFormats/Common/interface/Handle.h"
#include "DataFormats/Candidate/interface/Candidate.h"
#include "DataFormats/Candidate/interface/CandidateFwd.h"
#include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"


// UABaseTree Analysis class decleration
#include "UATree/UABaseTree/interface/UABaseTree.h"


bool PUSumInfoDebug = false ;

void UABaseTree::GetPUSumInfo(const edm::Event& iEvent)
{

   pusuminfo.Reset();
   
   edm::Handle<std::vector<PileupSummaryInfo> > puInfoH;
   std::vector<PileupSummaryInfo>::const_iterator PVI;
   try {
     iEvent.getByLabel(pusuminfo_,puInfoH);
     for(PVI = puInfoH->begin(); PVI != puInfoH->end(); ++PVI) {
	int BX = PVI->getBunchCrossing();
	if(BX == 0){
	   pusuminfo.nPU = PVI->getPU_NumInteractions();
	   for(int i=0; i < PVI->getPU_NumInteractions(); ++i){
	      pusuminfo.zposition.push_back(PVI->getPU_zpositions()[i]);
	      pusuminfo.sumpT_lowpT.push_back(PVI->getPU_sumpT_lowpT()[i]);
	      pusuminfo.sumpT_highpT.push_back(PVI->getPU_sumpT_highpT()[i]);
	      pusuminfo.ntrks_lowpT.push_back(PVI->getPU_ntrks_lowpT()[i]);
	      pusuminfo.ntrks_highpT.push_back(PVI->getPU_ntrks_highpT()[i]);
	   }
	}
        //FIXME: Add out-of-time pile-up 
     }
   }
   catch (...){
     cout << "[UABaseTree::GetPUSumInfo] Was not able to retrieve PU Summary" << endl;
   }
   
   if(PUSumInfoDebug) pusuminfo.Print();
}

