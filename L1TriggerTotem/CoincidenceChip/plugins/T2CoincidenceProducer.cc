#include "L1TriggerTotem/CoincidenceChip/interface/T2CoincidenceProducer.h"
#define ZMax 0

//
// constructors and destructor
//
T2CoincidenceProducer::T2CoincidenceProducer(const edm::ParameterSet& iConfig) :
    src_(iConfig.getParameter<edm::InputTag> ("src")), verbose_(iConfig.getParameter<unsigned int> ("verbose")), quarter_(iConfig.getParameter<unsigned int> ("quarter")) {
    //src_(iConfig.getParameter<edm::InputTag> ("src")), verbose_(iConfig.getParameter<bool> ("verbose")) {
    //verbose_(iConfig.getParameter<bool> ("verbose")) {

  // configure the coincidence chip
  cchip.configure(iConfig.getParameterSet("coincidenceChipConfig"));

  // now do what ever other initialization is needed
  produces < std::vector < T2TriggerBits > > ("T2TriggerBits");

  T2CoincidenceProducer::createMetaPadMap();
  hh = 1;// if =0, am debugging :: FO14dec2009 no longer in use, use verbose_ > 2 instead


}

T2CoincidenceProducer::~T2CoincidenceProducer() {
  // do anything here that needs to be done at destruction time
  // (e.g. close files, deallocate resources etc.)
}

//
// member functions
//

// ------------ method called to produce the data  ------------
void T2CoincidenceProducer::produce(edm::Event& iEvent, const edm::EventSetup& iSetup) {
  using namespace edm;
  using namespace std;
  edm::Handle<T2PadDigiCollection> input;
  iEvent.getByLabel(src_, input);

  if (verbose_)
    edm::LogPrint("T2CoincidenceProducer")<<"new event FOFOFO"<<std::endl;


  // reset the metaPadMap
  T2CoincidenceProducer::resetMetaPadMap();


  T2PadDigiCollection::DigiRangeIterator detUnitIt;
  std::bitset<16> processedBits;
  std::bitset<8> T2outputBits;
  if( verbose_ > 2) edm::LogPrint("T2CoincidenceProducer") << "CoincidenceChipLoop coming up";

  //  int counter=0;
  std::bitset<80> logicBits[4][numOfMetaPadColumns];
  for (int jj=0;jj<4;jj++) {for (int ii=0;ii<numOfMetaPadColumns;ii++) logicBits[jj][ii].reset();}
 


 const unsigned int zEff=6;
 const int quarterLook=quarter_; // arm=1 = z<0 = minus & halfTele=inner = near=0
//run 2137, minusNear should be=1, but try the other 3 possibilities too


 // access and iterate over digicontainer

  for (detUnitIt=input->begin(); detUnitIt != input->end(); ++detUnitIt) {
    const T2PadDigiCollection::Range& range = (*detUnitIt).second;
 
    for (T2PadDigiCollection::const_iterator digiIt = range.first; digiIt != range.second; ++digiIt) {



      // filling the MetaPadDetector with information from Digitisation
      int mPadRow = digiIt->getRow()/3;
      int mPadCol = digiIt->getCol()/5;

      // for simplification of TriggerEmulation Back-planes will be numbered such that
      // rows and columns considate with Front-planes
      if ( ( ((*detUnitIt).first).planeSide() ) == 0 ) {

      metaPadMap[((*detUnitIt).first).rawId()][mPadRow + mPadCol * numOfMetaPadRows] = 1;
      
      } else metaPadMap[((*detUnitIt).first).rawId()]
	       [mPadRow + (numOfMetaPadColumns-1 - mPadCol) * numOfMetaPadRows] = 1;
      
      //  std::cout<<"Row: "<<digiIt->getRow()<<", Col: "<<digiIt->getCol()<<std::endl;




    }
  }

  //Fredrik
  std::vector<int> ccZ;
  ccZ.clear();

  for(int i=0;i<4;i++){
    
    simplePaths[i]=0;
  }

  simplePath(simplePaths);

  //  //  if (verbose_) std::cout << "FT2-first-arm  "<<std::endl;
 if( verbose_) 
   edm::LogPrint("T2CoincidenceProducer") << "T2-CC-triggerbits  ";

  std::map<int, int*>::iterator FT2Ind;
  T2DetId FT2_x;
  
  for(int k=0;k<numOfMetaPadColumns; k++){
    //Loop over the 13 sectors
    int counterT2[4];
    for (int ii=0;ii<4;ii++) counterT2[ii]=0;
  
    for (FT2Ind=metaPadMap.begin(); FT2Ind != metaPadMap.end(); ++FT2Ind) {
      
      //    std::cout << (*FT2Ind).first << std::endl;
      //    if ((*FT2Ind).first.arm()==0) {
      FT2_x=(T2DetId) (*FT2Ind).first;
      int tempQuadrant=2*FT2_x.arm()+FT2_x.halfTelescope();
      //changed to correct H0...H3 numbering, from arm+2*halftele

      for(int l=0;l<numOfMetaPadRows; l++) {
	
	//	if ((*FT2Ind).second[l+k*numOfMetaPadRows]&&counterT2[FT2_x.arm()+2*FT2_x.halfTelescope()]<numOfMetaPadColumns) {

	//	if (counterT2[tempQuadrant]<numOfMetaPadColumns) {
	  //	  int temp1=counterT2[tempQuadrant];

	if (counterT2[tempQuadrant]<10) { //numOfPlanes, with 8bits/plane
	  if ((*FT2Ind).second[l+k*numOfMetaPadRows]) logicBits[tempQuadrant][k].set(l+counterT2[tempQuadrant]*numOfMetaPadRows);
	} else {
	  //	  std::cout<<"ERROR: MetaPadInfo-More than 10 of them: "<<counterT2[tempQuadrant]<<" "<<std::endl;
	  edm::LogWarning("T2CoincidenceProducer") << "More than 10 MetaPadPlanes in a T2 quadrant";

	}
	//	std::cout << "RawPad: "<< (*FT2Ind).second[l+k*numOfMetaPadRows] <<std::endl;
	// as I assume above, .second is only= 0 or 1
	
	//	} else {
	
	//	}
	//	if ((*FT2Ind).second[l+k*numOfMetaPadRows]) logicBits[k].set(l+k*numOfMetaPadRows);
	//	if ((*FT2Ind).second[l+k*numOfMetaPadRows]) logicBits[FT2_x.arm()+2*FT2_x.halfTelescope()].set();
	//	if ((*FT2Ind).second[l+k*numOfMetaPadRows]) logicBits[(*FT2Ind).first.arm()+2*(*FT2Ind).first.halfTelescope()].set();
	
	//std::cout << "RawId: "<< (*FT2Ind).first << " T2Cast: " <<   ((T2DetId) (*FT2Ind).first)<<std::endl; //.arm();
      }
      counterT2[tempQuadrant]++;
    }
  }
  std::auto_ptr<std::vector<T2TriggerBits> > output(new std::vector<T2TriggerBits>);
    
  const int size=numOfMetaPadColumns;
  output->reserve(size);
  if( verbose_ > 2) edm::LogPrint("T2CoincidenceProducer") << "T2CC debug output follows";

 
  
  for (int jj=0;jj<4;jj++) {
    for (int ii=0;ii<numOfMetaPadColumns;ii++) {
      processedBits=cchip.process(logicBits[jj][ii]);
      if (logicBits[jj][ii].count()+processedBits.count()){
	if( hh&&(verbose_>2))
	  edm::LogPrint("T2CoincidenceProducer") << "NumSetBits L "<<logicBits[jj][ii].count()<<"  P "<<processedBits.count();

	//	if( verbose_&&hh==0 ) std::cout<<"NumSetBits L "<<logicBits[jj][ii].count()<<"  P "<<processedBits.count()<<std::endl;
      }
      //      const uint32_t rawChipId= (uint32_t) ii+jj*numOfMetaPadColumns;
      const T2DetId rawChipId=T2DetId(jj/2, jj%2,0,0); 

      //FO:10dec08 debug
      //     if( verbose_ ) std::cout<<"CC Hit! L:"<<processedBits.to__string()<<std::endl;
      
      if (processedBits.count()||(logicBits[jj][ii].count())) {
	//      if (processedBits.count()||(logicBits[jj][ii].count()>5)) {
	//for-loop necessary since (bitset a).to_long() works, but not a.to_string()
	if( hh&&(verbose_>2))
	  edm::LogPrint("T2CoincidenceProducer") << "CC Hit? L:";

	for (int u=0;u<10;u++) {
	  std::bitset<8> temp1b8;
	  
	  for (int u2=0;u2<8;u2++) {
	    //	    if( verbose_&&hh==0 ) std::cout<<(logicBits[jj][ii].test(u*8+u2)&& (int) 1);
	    temp1b8[u2] = logicBits[jj][ii].test(u*8+u2);
	  }
	  if(( hh)&&(jj==quarterLook)&&(verbose_>2))
	    edm::LogPrint("T2CoincidenceProducer") << temp1b8;
	  
	  //	  if( verbose_&&hh==0 ) std::cout<<std::endl;
	}
	if( (hh)&&(jj==quarterLook)&&(verbose_>2))
	  edm::LogPrint("T2CoincidenceProducer") << "  P:";
	
	//	if( verbose_&&hh==0 ) std::cout<<"  P:";
	/*	for (int u=0;u<16;u++) {
		if( verbose_&&hh==0 ) std::cout<<(processedBits.test(u)&& (int) 1);
		}
	*/
	if( (hh)&&(jj==quarterLook)&&(verbose_>2))
	  edm::LogPrint("T2CoincidenceProducer") << processedBits;

	//	if( verbose_&&hh==0 ) std::cout<<std::endl;
      }
      for (int u=0;u<8;u++) T2outputBits[u]=processedBits[u];
      //save the first 8 bits of the CC's 16bit output; need change only if
      //CC code is materially changed, or parameter OV is set nonzero

      if ((T2outputBits.count()+processedBits.count()>ZMax)&&(jj==quarterLook)&&verbose_)
	edm::LogPrint("T2CoincidenceProducer")<<"FO number of bits on in T2outputBits:"<<T2outputBits.count()<<" and in processedBits(16bits, take half for T2):"<<processedBits.count()<<std::endl;
      if ((T2outputBits.count()+processedBits.count()>ZMax)&&(jj!=quarterLook)&&verbose_)
        edm::LogPrint("T2CoincidenceProducer")<<"FOFO CC triggered in nontriggering quarter: number of bits on in T2outputBits:"<<T2outputBits.count()<<" and in processedBits:"<<processedBits.count()<<std::endl;

      if (jj==quarterLook)
//      if (1)
	{
	  if (T2outputBits.count()>zEff) 
	    ccZ.push_back(2);
	  if ((T2outputBits.count()<zEff+1) &&(T2outputBits.count()>0))
	    ccZ.push_back(1);

	  if (!T2outputBits.count())
	    ccZ.push_back(0);
	      
	    
	}
    
      const T2TriggerBits coincidenceBits(rawChipId, ii, T2outputBits);
      output->push_back(coincidenceBits);
    }
  }



  //if (hh==0)  printoutStuff(); // debugging


  if (ccZ.size()>52) 
    edm::LogPrint("T2CoincidenceProducer")<<"More than 52 CC chips in event!!"<<std::endl;
  if (ccZ.size()<1) 
    edm::LogPrint("T2CoincidenceProducer")<<"Zero CC chips in event!!"<<std::endl;
  bool onlyZeroes=true;
  bool onlyZFailAnd0=false;

  if (verbose_>2) 
    edm::LogPrint("T2CoincidenceProducer")<<"FO CC classify:";
 
  std::vector<int>::const_iterator thisCC;
  for (thisCC=ccZ.begin();thisCC!=ccZ.end();thisCC++)
    {
      if (verbose_>2) 
	edm::LogPrint("T2CoincidenceProducer")<<(*thisCC);
      if ((*thisCC)==1)
	{
	  onlyZeroes=false;
	  onlyZFailAnd0=false;
	}
      if (((*thisCC)==2)&&(onlyZeroes)) 
	{
	  onlyZFailAnd0=true;
	  onlyZeroes=false;
	}
    }
  if (verbose_>2)
    edm::LogPrint("T2CoincidenceProducer")<<std::endl;

  if (onlyZeroes&&verbose_) 
    edm::LogPrint("T2CoincidenceProducer")<<"FOFO Event is all zeroes"<<std::endl;
  
  if (onlyZFailAnd0&&verbose_) 
    edm::LogPrint("T2CoincidenceProducer")<<"FOFO Event has only zeroed CCs and Z-filtered CCs. Z="<<zEff<<std::endl;


  if( hh&&verbose_ ) 
    LogPrint("CoincidenceChip") << "output size : " << output->size();
  iEvent.put(output, "T2TriggerBits");
}

// ------------ method called once each job just before starting event loop  ------------
void T2CoincidenceProducer::beginJob() {
  if( verbose_ ) edm::LogPrint("T2CoincidenceProducer") << "T2CoincidenceProducer: beginJob";

  /*
  cchip.config().setV(minNumOfPlanes);
  cchip.config().setNP(0x0); //8 groups of 10 planes
  cchip.config().setOV(0x0);
  cchip.config().setW(0x0);
  cchip.config().setZ(0xF);
  cchip.config().setO2(0x0);
  cchip.config().setLI(0x0);
  cchip.config().setLO(0x0);
  cchip.config().setAO(0x0);
  */

  if( verbose_ ) edm::LogPrint("CoincidenceChip") << "configuration : " << cchip.getControlReg();
  if( verbose_ ) edm::LogPrint("CoincidenceChip") << "cfg summary   : " << cchip.configSummary();
}

// ------------ method called once each job just after ending the event loop  ------------
void T2CoincidenceProducer::endJob() {
  if( verbose_ ) edm::LogPrint("T2CoincidenceProducer") << "T2CoincidenceProducer: endJob";
}










/**
 *
 */

void T2CoincidenceProducer::createMetaPadMap() {
  
  for(int i=0; i<2; i++) {
    
    for(int j=0; j<2; j++) {
      
      for(int k=0; k<5; k++) {
        
        for(int l=0; l<2; l++) {
          
          int* tmp = new int[numOfMetaPadColumns * numOfMetaPadRows];

          int detId = T2DetId::calculateRawId(i, j, k, l);

          metaPadMap[detId] = tmp;

        }
      }
    }
  }

} // createMetaPadMap

/**
 *
 */

void T2CoincidenceProducer::resetMetaPadMap() {
  
  std::map<int, int*>::iterator i;

  for (i=metaPadMap.begin(); i != metaPadMap.end(); ++i) {

    for (int j=0; j < numOfMetaPadColumns * numOfMetaPadRows; j++) {
      
      (*i).second[j] = 0;
    }  
  }  

} // resetMetaPadMap

/**
 *
 */

 void T2CoincidenceProducer::deleteMetaPadMap() {
  
  std::map<int, int*>::iterator i;
  
  for (i=metaPadMap.begin(); i != metaPadMap.end(); ++i) {
    
    delete[] (*i).second;
  }  

  metaPadMap.clear();  

} // resetMetaPadMap

/**
 * to be deleted !
 */
/*
void T2CoincidenceProducer::printoutStuff() {
  
  std::map<int, int*>::iterator i;
  int counter=0;
  for (i=metaPadMap.begin(); i != metaPadMap.end(); ++i) {
    
    std::cout<<std::endl;
    std::cout << (*i).first << std::endl;
    
    for(int k=0;k<numOfMetaPadColumns; k++){
      for(int l=0;l<numOfMetaPadRows; l++) {
	
 	std::cout << " " << (*i).second[l+k*numOfMetaPadRows];
      }
      
      std::cout<<std::endl;
    }
    
    std::cout << counter << std::endl;
    counter++;

  }
  
} // printoutStuff
*/

/**
 *  Simple algorithm to scan for straight path pattern in T2 MetaPadDetector
 */

void T2CoincidenceProducer::simplePath(int* count) {
  
  std::map<int, int*>::iterator it;
  
  for (int row=0; row<numOfMetaPadRows; ++row) {
    
    for (int col=0; col<numOfMetaPadColumns; ++col) {
      
      int counter = 0;
      int planeNr = 0;
      
      for ( it=metaPadMap.begin(); it != metaPadMap.end(); ++it ) {
	
	if ( (*it).second[row + col * numOfMetaPadRows] == 1) counter++;
	
	if (planeNr == 9) {
	  
	  if (counter >= minNumOfPlanes) count[0]++;
	  counter = 0;
	}
	
	if (planeNr == 19) {
	  
	  if (counter >= minNumOfPlanes) count[1]++;
	  counter = 0;
	}
	
	if (planeNr == 29) {
	  
	  if (counter >= minNumOfPlanes) count[2]++;
	  counter = 0;
	}
	
	if (planeNr == 39) {
	  
	  if (counter >= minNumOfPlanes) count[3]++;
	  counter = 0;
	}
	
	planeNr++;
      }
    }
  }
  
} // simplePath




//define this as a plug-in
DEFINE_FWK_MODULE(T2CoincidenceProducer);
