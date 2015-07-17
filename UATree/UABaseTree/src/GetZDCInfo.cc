//-- Description: Function to retrieve ZDC information

//-- Dataformats
#include "DataFormats/HcalDigi/interface/HcalDigiCollections.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/HcalRecHit/interface/ZDCRecHit.h"
#include "CalibFormats/HcalObjects/interface/HcalCoderDb.h"
#include "CalibFormats/HcalObjects/interface/HcalDbService.h"
#include "CalibFormats/HcalObjects/interface/HcalDbRecord.h"

#include "UATree/UABaseTree/interface/UABaseTree.h"

bool ZDCInfoDebug = false;

void UABaseTree::GetZDCInfo(const edm::Event& event, const edm::EventSetup& setup) {

   zdcInfo.Reset();
   zdcHits.clear();
   zdcDigis.clear();

   // Get the ZDC rechits collection from the event
   edm::Handle<ZDCRecHitCollection> zdcHitsH;
   event.getByLabel(zdcrechits_, zdcHitsH);

   std::map<int,int> nHitsZDCPerChannel;
   std::map<int,double> sumEnergyZDCPerChannel;

   if( zdcHitsH.isValid() ){ 

      const ZDCRecHitCollection* zdchits = zdcHitsH.product();
      ZDCRecHitCollection::const_iterator zdchit;

      MyZDCHit myzdchit;
      /*int nZDChitCand = 0;
	double ZDCsumHADminus = 0.;
	double ZDCsumEMminus = 0.;
	double ZDCsumHADplus = 0.;
	double ZDCsumEMplus = 0.;*/
      for ( zdchit = zdchits->begin(); zdchit != zdchits->end(); ++zdchit )
      {
	 HcalZDCDetId id(zdchit->id());
	 int side      = id.zside();
	 int section   = id.section();
	 int channel   = id.channel();
	 //int channelId = (section-1)*5+(side+1)/2*9+(channel-1);
	 int channelId = (section-1)*5+(side+1)/2*17+(channel-1);

	 double energy = zdchit->energy();
	 double time = zdchit->time();

	 /*if((section == 1) && (side == 1))  ZDCsumEMplus   += energy;
	   if((section == 1) && (side == -1)) ZDCsumEMminus  += energy;
	   if((section == 2) && (side == 1))  ZDCsumHADplus  += energy;
	   if((section == 2) && (side == -1)) ZDCsumHADminus += energy;
	   ++nZDChitCand;*/
         if(nHitsZDCPerChannel.find(channelId) == nHitsZDCPerChannel.end()){
            nHitsZDCPerChannel[channelId] = 0; 
            sumEnergyZDCPerChannel[channelId] = 0.;
         }
         ++nHitsZDCPerChannel[channelId];
         sumEnergyZDCPerChannel[channelId] += energy;

	 if(storeZDCHits_){
	    myzdchit.Reset();
	    myzdchit.side = side;
	    myzdchit.section = section;
	    myzdchit.channel = channel;
	    myzdchit.channelId = channelId;
	    myzdchit.energy = energy;
	    myzdchit.time = time;
	    zdcHits.push_back(myzdchit);
	    if (ZDCInfoDebug) zdcHits.back().Print();   
	 }
      }
   }

   edm::Handle<ZDCDigiCollection> zdcDigisH;
   event.getByLabel(zdcdigis_, zdcDigisH);

   if( zdcDigisH.isValid() ){

      const ZDCDigiCollection* zdcdigis = zdcDigisH.product();
      ZDCDigiCollection::const_iterator zdcdigi;

      edm::ESHandle<HcalDbService> conditions;
      setup.get<HcalDbRecord>().get(conditions);

      MyZDCDigi myzdcdigi;
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
	 int channelId = (section-1)*5+(side+1)/2*17+(channel-1);

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

	 if( !section ) continue; 

	 const HcalQIEShape* qieshape=conditions->getHcalShape(digi.id());
	 const HcalQIECoder* qiecoder=conditions->getHcalCoder(digi.id());
	 CaloSamples caldigi;
	 HcalCoderDb coder(*qiecoder,*qieshape);

	 coder.adc2fC(digi,caldigi);

	 int fTS = digi.size();
         std::vector<int> digiADC(fTS);
         std::vector<float> digifC(fTS);
	 for(int iTS = 0; iTS < fTS; ++iTS) {
            digiADC[iTS] = digi[iTS].adc();
            digifC[iTS]  = caldigi[iTS];
	 }

	 if(storeZDCDigis_){
	    myzdcdigi.Reset();
	    myzdcdigi.side = side;
	    myzdcdigi.section = section;
	    myzdcdigi.channel = channel;
	    myzdcdigi.channelId = channelId;
	    myzdcdigi.digiADC = digiADC;
	    myzdcdigi.digifC = digifC;
	    zdcDigis.push_back(myzdcdigi);
	    if (ZDCInfoDebug) zdcDigis.back().Print();   
	 }
      }
   }
 
   if(storeZDCInfo_){
      zdcInfo.nHitsPerChannel = nHitsZDCPerChannel;
      zdcInfo.sumEnergyPerChannel = sumEnergyZDCPerChannel;
      if (ZDCInfoDebug) zdcInfo.Print();   
   }
}
