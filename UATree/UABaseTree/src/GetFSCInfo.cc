//-- Description: Function to retrieve FSC rec. hit information

//-- Dataformats
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool FSCInfoDebug = false;

void UABaseTree::GetFSCInfo(const edm::Event& event, const edm::EventSetup& setup) {

   fscInfo.Reset();
   fscHits.clear();
   fscDigis.clear();

   // Get the ZDC rechits collection from the event
   edm::Handle<ZDCRecHitCollection> zdcHitsH;
   event.getByLabel(fscrechits_, zdcHitsH);

   std::map<int,int> nHitsFSCPerChannel;
   std::map<int,double> sumEnergyFSCPerChannel;

   if( zdcHitsH.isValid() ){ 

      const ZDCRecHitCollection* zdchits = zdcHitsH.product();
      ZDCRecHitCollection::const_iterator zdchit;

      MyFSCHit myfschit;
      /*int nFSChitCand = 0;
	double FSCsumHADminus = 0.;
	double FSCsumEMminus = 0.;
	double FSCsumHADplus = 0.;
	double FSCsumEMplus = 0.;*/
      for ( zdchit = zdchits->begin(); zdchit != zdchits->end(); ++zdchit )
      {
	 HcalZDCDetId id(zdchit->id());
	 int side      = id.zside();
	 int section   = id.section();
	 int channel   = id.channel();
	 int channelId = (section-1)*5+(side+1)/2*9+(channel-1);
	 double energy = zdchit->energy();
	 double time = zdchit->time();

	 /*if((section == 1) && (side == 1))  FSCsumEMplus   += energy;
	   if((section == 1) && (side == -1)) FSCsumEMminus  += energy;
	   if((section == 2) && (side == 1))  FSCsumHADplus  += energy;
	   if((section == 2) && (side == -1)) FSCsumHADminus += energy;
	   ++nFSChitCand;*/
         if(nHitsFSCPerChannel.find(channelId) == nHitsFSCPerChannel.end()){
            nHitsFSCPerChannel[channelId] = 0; 
            sumEnergyFSCPerChannel[channelId] = 0.;
         }
         ++nHitsFSCPerChannel[channelId];
         sumEnergyFSCPerChannel[channelId] += energy;

	 if(storeFSCHits_){
	    myfschit.Reset();
	    myfschit.side = side;
	    myfschit.section = section;
	    myfschit.channel = channel;
	    myfschit.channelId = channelId;
	    myfschit.energy = energy;
	    myfschit.time = time;
	    fscHits.push_back(myfschit);
	    if (FSCInfoDebug) fscHits.back().Print();   
	 }
      }
   }

   edm::Handle<ZDCDigiCollection> zdcDigisH;
   event.getByLabel(fscdigis_, zdcDigisH);

   if( zdcDigisH.isValid() ){

      const ZDCDigiCollection* zdcdigis = zdcDigisH.product();
      ZDCDigiCollection::const_iterator zdcdigi;

      edm::ESHandle<HcalDbService> conditions;
      setup.get<HcalDbRecord>().get(conditions);

      MyFSCDigi myfscdigi;
      for(zdcdigi = zdcdigis->begin(); zdcdigi != zdcdigis->end(); ++zdcdigi){
         /*
         const ZDCDataFrame digi = (const ZDCDataFrame)(*j);		
	 int iSide      = 1+2*(12-digi.elecId().spigot());
	 int iSection   = digi.id().section();
	 int iChannel   = digi.id().channel();
	 bool isFsc=(digi[0].fiber()<2||digi[0].fiberChan()<6);

	 int fscCh = isFsc?(9+(iSide+1)/2*17+(digi[0].fiber()-4)*3+digi[0].fiberChan()):0;
	 int chid = iSection?((iSection-1)*5+(iSide+1)/2*17+(iChannel-1)):fscCh;
         */
	 const ZDCDataFrame digi = (const ZDCDataFrame)(*zdcdigi);
	 //int side      = digi.id().zside();
	 int side      = 1+2*(12-digi.elecId().spigot());
	 int section   = digi.id().section();
	 int channel   = digi.id().channel();
	 //int channelId = (section-1)*5+(side+1)/2*9+(channel-1);

	 bool isFsc    = (digi[0].fiber()<2||digi[0].fiberChan()<6);
	 int fscCh     = isFsc?(9+(side+1)/2*17+(digi[0].fiber()-4)*3+digi[0].fiberChan()):0;
         int channelId = fscCh;

         /*
         const HcalQIEShape* qieshape=conditions->getHcalShape();
  	 const HcalQIECoder* qiecoder=iSection?conditions->getHcalCoder(digi.id()):0;
	 CaloSamples caldigi;

	 if(iSection){ 
	 	HcalCoderDb coder(*qiecoder,*qieshape);
		coder.adc2fC(digi,caldigi);
	 }
	
	 int fTS = digi.size();
			
	 if(iSection||isFsc){
	   for (int i = 0; i < fTS; ++i) {
		   DigiDatafC[i+chid*10] = iSection?caldigi[i]:digi[i].nominal_fC();
		   DigiDataADC[i+chid*10] = digi[i].adc();
	   }
	 }
         */

         if( section || !isFsc ) continue;

	 /*const HcalQIEShape* qieshape=conditions->getHcalShape();
	 const HcalQIECoder* qiecoder=conditions->getHcalCoder(digi.id());
	 CaloSamples caldigi;
	 HcalCoderDb coder(*qiecoder,*qieshape);

	 coder.adc2fC(digi,caldigi);*/

	 int fTS = digi.size();
         std::vector<int> digiADC(fTS);
         std::vector<float> digifC(fTS);
	 for(int iTS = 0; iTS < fTS; ++iTS) {
            digiADC[iTS] = digi[iTS].adc();
            //digifC[iTS]  = caldigi[iTS];
            digifC[iTS]  = digi[iTS].nominal_fC();
	    /*DigiDatafC[i+chid*10] = caldigi[i];
	    DigiDataADC[i+chid*10] = digi[i].adc();*/
	 }

	 if(storeFSCDigis_){
	    myfscdigi.Reset();
	    myfscdigi.side = side;
	    myfscdigi.section = section;
	    myfscdigi.channel = channel;
	    myfscdigi.channelId = channelId;
	    myfscdigi.digiADC = digiADC;
	    myfscdigi.digifC = digifC;
	    fscDigis.push_back(myfscdigi);
	    if (FSCInfoDebug) fscDigis.back().Print();   
	 }
      }
   }
 
   if(storeFSCInfo_){
      /*fscInfo.nHits = nFSChitCand;
      fscInfo.sumEnergyEMPlus = FSCsumEMplus;
      fscInfo.sumEnergyEMMinus = FSCsumEMminus;
      fscInfo.sumEnergyHADPlus = FSCsumHADplus;
      fscInfo.sumEnergyHADMinus = FSCsumHADminus;*/
      fscInfo.nHitsPerChannel = nHitsFSCPerChannel;
      fscInfo.sumEnergyPerChannel = sumEnergyFSCPerChannel;
      if (FSCInfoDebug) fscInfo.Print();   
   }
}
