#include "TotemT1T2Validation/T2GeometryValidation/interface/T2GeometryAnalyzer.h"
//#include "TotemT1T2Validation/T2GeometryValidation/interface/DAQInformationT2_a.h"
#include "TotemT1T2Validation/T2GeometryValidation/interface/DAQInformationSourceXML_a.h"
#include <memory>
#include <string>
#include <iostream>
#include <vector>
#include <map>
#include <cmath>

#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "TotemRawData/RawToDigi/interface/T2Mapping.h"
#include "TotemRawDataLibrary/DataFormats/interface/RawEvent.h"
#include "TotemRawDataLibrary/DataFormats/interface/VFATFrame.h"
#include "TotemRawDataLibrary/Utilities/interface/CommonDef.h"
#include "TotemCondFormats/DAQInformation/interface/DAQInformationT2.h"
#include "TotemCondFormats/DataRecord/interface/TotemDAQRecord.h"

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/T2Digi/interface/T2StripDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2StripDigi.h"
#include "DataFormats/T2Digi/interface/T2PadDigiCollection.h"
#include "DataFormats/T2Digi/interface/T2PadDigi.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "DataFormats/GeometryVector/interface/LocalPoint.h"

#include "TFile.h"
#include "TNtuple.h"
#include "TH1.h"
#include "TH2.h"
#include "TGraph.h"
#include "TTree.h"

#include "TMath.h"
#include "CLHEP/Vector/LorentzVector.h"


T2GeometryAnalyzer::T2GeometryAnalyzer(const edm::ParameterSet& iConfig){

	T2StripDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2StripDigiCollectionLabel");
T2PadDigiCollectionLabel = iConfig.getParameter<edm::InputTag>("T2PadDigiCollectionLabel");

  Chicut = iConfig.getParameter<double>("Chicut");
  PhiChiProbCut = iConfig.getParameter<double>("PhiChiProbCut");
  RChiProbCut = iConfig.getParameter<double>("RChiProbCut");  
  energy=iConfig.getParameter<double>("FitdzEnergy");
  DZScale=iConfig.getParameter<double>("DZScale");
  tracketamin=iConfig.getParameter<double>("TrkEtamin");
  tracketamax=iConfig.getParameter<double>("TrkEtaMAX");
  singleparticle=iConfig.getParameter<bool>("singleparticle");

  outputFileName = iConfig.getUntrackedParameter<std::string>("OutputFile");
  HepMCProductLabel = iConfig.getParameter<std::string>("HepMCProductLabel");
  CluLabel = iConfig.getParameter<std::string>("CluLabel");
  HitLabel = iConfig.getParameter<std::string>("HitLabel");
  RoadLabel = iConfig.getParameter<std::string>("RoadLabel");
  TrackLabel= iConfig.getParameter<std::string>("TrackLabel");
  DoClusterFromVFatChannel=iConfig.getParameter<bool>("DoClusterFromVFatChannel");
  PrintOnlyFinalOutput=iConfig.getParameter<bool>("PrintOnlyFinalOutput");
  PrintChannelNumerationAsInVfatController=iConfig.getParameter<bool>("PrintChannelNumerationAsInVfatController");
  OnlyVfatANDchannel_Finalprinting=iConfig.getParameter<bool>("OnlyVfatANDchannel_Finalprinting");
}


T2GeometryAnalyzer::~T2GeometryAnalyzer()
{
}


void T2GeometryAnalyzer::InfoForVfatController(int* CCuAddrInt, std::string* Fec,  std::string* CCuAddr, std::string* I2CChan, int absvfatiid)
{

  int plane=absvfatiid/100;
  int iid=absvfatiid%100;

  if(plane>=20)
    {
      (*Fec)="19,7";
    }
  else
    {
      (*Fec)="19,8";
    }


  if(plane<10)
    {
      if((plane<2)&&(plane>=0))
	(*CCuAddr)= "0x5b";
      if((plane<4)&&(plane>=2))
	(*CCuAddr)= "0x5c";
      if((plane<6)&&(plane>=4))
	(*CCuAddr)= "0x5d";
      if((plane<8)&&(plane>=6))
	(*CCuAddr)= "0x5e";
      if((plane<10)&&(plane>=8))
	(*CCuAddr)= "0x5f";
    }
    
  if((plane>=10)&&(plane<20))
    {
      if((plane<12)&&(plane>=10))
	(*CCuAddr)= "0x61";
      if((plane<14)&&(plane>=12))
	(*CCuAddr)= "0x75";
      if((plane<16)&&(plane>=14))
	(*CCuAddr)= "0x66";
      if((plane<18)&&(plane>=16))
	(*CCuAddr)= "0x6b";
      if((plane<20)&&(plane>=18))
	(*CCuAddr)= "0x6c";
    }

  
  //0xb 0xc 0xd 0xe 0x1a
  if((plane>=20)&&(plane<30))
    {
      if((plane<22)&&(plane>=20))
	(*CCuAddr)= "0xb";
      if((plane<24)&&(plane>=22))
	(*CCuAddr)= "0xc";
      if((plane<26)&&(plane>=24))
	(*CCuAddr)= "0xd";
      if((plane<28)&&(plane>=26))
	(*CCuAddr)= "0xe";
      if((plane<30)&&(plane>=28))
	(*CCuAddr)= "0x1a";
    }
 
  // 0x70 ..  0x74
  if((plane>=30)&&(plane<40))
    {
      if((plane<32)&&(plane>=30))
	(*CCuAddr)= "0x70";
      if((plane<34)&&(plane>=32))
	(*CCuAddr)= "0x71";
      if((plane<36)&&(plane>=34))
	(*CCuAddr)= "0x72";
      if((plane<38)&&(plane>=36))
	(*CCuAddr)= "0x73";
      if((plane<40)&&(plane>=38))
	(*CCuAddr)= "0x74";
    }


  // std::cout<<"HEre"<<std::endl;
  std::string mah= (*CCuAddr);
  const char *s = mah.c_str();
  //  std::cout<<mah.c_str()<<std::endl;
  //ul = strtoul (szInput,NULL,0x);

  (*CCuAddrInt) = strtol(s,NULL,16); 

  //(*CCuAddrInt)= boost::lexical_cast<unsigned int>(mah);//("0xdeadbeef");
  //std::cout<<"string "<<(*CCuAddr).c_str()<<" converted as"<<(*CCuAddrInt)<<std::endl;

if((plane%2)==0)
{
  if((iid<=7)&&(iid>=0))
    (*I2CChan)="1a";

  if(iid==8)
    (*I2CChan)="11";
  
  if((iid<=16)&&(iid>=9))
    (*I2CChan)="12";
}
else
{
  if((iid<=7)&&(iid>=0))
    (*I2CChan)="14";

  if(iid==8)
    (*I2CChan)="15";
  
  if((iid<=16)&&(iid>=9))
    (*I2CChan)="16";

}

  /*

fec 19,7 -
fec 19,8 +



  Quarto "D3" (Arm plus, Side Near )

ccu 5b  HS 0-1
ccu 5c  HS 2-3
ccu 5d  HS 4-5
ccu 5e  HS 6-7
ccu 5f   HS 8-9


Quarto "D2" (Arm plus, Side Far)
ccu 61  HS 0-1
ccu 75  HS 2-3
ccu 66  HS 4-5
ccu 6b  HS 6-7
ccu 6c   HS 8-9


Nota: qui soto ci sono alcuni puntini per il ccu address. Questi address non li ho trascritti perche'
gli indirizzi sono in ordine nella tendina del vfat: ad esempio l'indirizzo che troverai sotto 0xb del D4 gestisce le hs 2-3

Quarto "D4" (Arm minus, Side Near)
ccu 0xb  HS 0-1
ccu ..  HS 2-3
ccu ..  HS 4-5
ccu ..  HS 6-7
ccu 0x1a   HS 8-9


Quarto "D1" (Arm minus, Side Far)
ccu 0x70  HS 0-1
ccu ..  HS 2-3
ccu ..  HS 4-5
ccu ..  HS 6-7
ccu 0x74  HS 8-9

Fredrik:

 ccu 0xb  HS 0-1 
 ccu 0xc  HS 2-3
 ccu 0xd HS 4-5
 ccu 0xe  HS 6-7
 ccu 0x1a   HS 8-9

 ccu 0x70  HS 0-1
 ccu 0x71  HS 2-3
 ccu 0x72  HS 4-5
 ccu 0x73  HS 6-7
 ccu 0x74  HS 8-9



Le Hs pari hanno I2C:

1a: contiene vfat iid da 0 a 7 (importante:sono posizionati in ordine nella tendina del vfat controller: cioe' il primo vfat che trovi per questo ccu i2c e' lo 0)
11:  contiene vfat iid 8
12:  contiene vfat iid da 9 a16



Le Hs dispari hanno I2C:
 
14 : contiene vfat iid da 0 a 7
15 : contiene vfat iid 8
16 : contiene vfat iid da 9 a16

  */

}



void T2GeometryAnalyzer::MakeEvCluster(std::auto_ptr<T2PadDigiCollection> PadDigiptr,std::auto_ptr<T2StripDigiCollection> StripDigiptr)
{


  std::vector<int> myvp;
  std::vector<int> myvs;


  DigiContainerIterator<T2DetId, T2PadDigi> itp;   
  std::vector<T2Cluster> myclusters; 
  
  for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)

    {
          
      T2DetId mydet=(*itp).first;  
      //mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;

	for(std::vector<T2PadDigi>::const_iterator itpad =(*itp).second.first; itpad !=(*itp).second.second; itpad++)
	{
	  //mypadsdigi=(*itpad).second; 

	  if((*itpad).getPadNr()==0)
	    myvp.push_back(((*itpad).getRow()+1)+24*((*itpad).getCol()) );   //Prima del dataraw
	  else
	    myvp.push_back(((*itpad).getPadNr())); 

	}
     

      theT2Clusterizer.SetDetId(mydet.rawId());

      theT2Clusterizer.SetPadHits(myvp);

      theT2Clusterizer.BuildClusters();

      myvp.clear();
      myclusters=theT2Clusterizer.GetPadClusters();
      
     
      thePadClusters->insert(std::make_pair(mydet,myclusters));
       
       
       
      myclusters.clear();
    }
  
 
  DigiContainerIterator<T2DetId, T2StripDigi> its;

  for(its= StripDigiptr->begin(); its!=StripDigiptr->end(); ++its)

    {
      T2DetId mydet=(*its).first;
      //symb=mydet.arm()*20+mydet.plane()*4+mydet.halfTelescope()*2+mydet.planeSide();
  //    mydet.plane()*2 + mydet.planeSide() + mydet.arm()*20 + mydet.halfTelescope()*10;
      

     
      for(std::vector<T2StripDigi>::const_iterator itstrip =(*its).second.first; itstrip !=(*its).second.second; itstrip++)
	{
	
	  if((*itstrip).getStripNr()==0)
	    myvs.push_back(((*itstrip).getRow()+1)+256*((*itstrip).getCol()) ); //PrimaDataraw
	  else
	    myvs.push_back(((*itstrip).getStripNr())); 
	       
	  //     cout<<" Detector pl: "<< mydet.plane() <<". Strip digitizzati di Erik: Row:  "<< (*itstrip).getRow() <<" , Col: "<< (*itstrip).getCol() <<endl;
    
	}
      theT2Clusterizer.SetDetId(mydet.rawId());
      theT2Clusterizer.SetStripHits(myvs); 
      theT2Clusterizer.BuildClusters();

      //Part added for noise studies      
      //StrCluStrRawID_=theT2Clusterizer.StrCluStrRawID;
      //StrCluStrColID_=theT2Clusterizer.StrCluStrColID;
      //PadCluPadRawID_=theT2Clusterizer.PadCluPadRawID;
      //PadCluPadColID_=theT2Clusterizer.PadCluPadColID;


      myvs.clear();
     
      myclusters=theT2Clusterizer.GetStripClusters();
	
      
      theStripClusters->insert(std::make_pair(mydet,myclusters));

    }
}





//
// member functions
//
// ------------ method called to for each event  ------------
void T2GeometryAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  using namespace std;
  using namespace edm;
  using namespace HepMC;

  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  /* LOADING OF ALL THE RECORDS FROM THE EVENT */
  /* @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@ */
  ESHandle<DAQInformationT2> Vfat_t2_mapping;
  iSetup.get<TotemDAQRecord>().get(Vfat_t2_mapping);

    /* :::::::::::::TakeDigi::::::::::::*/

  edm::Handle<T2PadDigiCollection> t2paddigicoll;
  iEvent.getByLabel(T2PadDigiCollectionLabel, t2paddigicoll);

  t2paddigicoll.product();
  edm::Handle<T2StripDigiCollection> t2stripdigicoll;
  iEvent.getByLabel(T2StripDigiCollectionLabel, t2stripdigicoll);
  t2stripdigicoll.product();
  DigiContainerIterator<T2DetId, T2PadDigi> itp;
  DigiContainerIterator<T2DetId, T2StripDigi> its;

  /* :::::::::::::Take The Clusters::::::::::::*/

  //std::cout << " check --------- A " <<  std::endl;
  Handle<T2StripClusterCollection> t2strclcoll;
  iEvent.getByLabel(CluLabel,"T2StripClusters",t2strclcoll);
  Handle<T2PadClusterCollection> t2padclcoll;
  iEvent.getByLabel(CluLabel,"T2PadClusters",t2padclcoll);

  //::::::Take particle generated with Gun::::::
  
  Handle<HepMCProduct> EvtHandle ;
  iEvent.getByLabel(HepMCProductLabel, EvtHandle ) ;
  //iEvent.getByLabel("source", EvtHandle ) ;
  EvtHandle->GetEvent();

 // ::::::Take Geant local hit on T2::::::
  Handle<CrossingFrame<PSimHit> > cFrame;
  iEvent.getByLabel("mix", "g4SimHitsTotemHitsT2Gem", cFrame);
  // iEvent.getByLabel("mix", "TotemHitsT2Gem", cFrame);
  // get hits from G4Sim
  const string nameOfHits("TotemHitsT2Gem");
  auto_ptr<MixCollection<PSimHit> >  T2SimHits( new MixCollection<PSimHit>( cFrame.product() ) );

  // map hits to T2 plane
  map<int, PSimHitContainer> hitMap;
  for(MixCollection<PSimHit>::MixItr hitItr = T2SimHits->begin(); hitItr != T2SimHits->end(); ++hitItr) {
    hitMap[hitItr->detUnitId()].push_back(*hitItr);
  }
  
  //::::::Take  T2  Hits::::::
  Handle<T2HitCollection> t2hitcoll;

  iEvent.getByLabel(HitLabel,"T2Hits",t2hitcoll);
  //::::::Take  T2  Roads::::::

  Handle<T2RoadCollection> t2roadcoll;
  iEvent.getByLabel(RoadLabel,"T2RoadColl",t2roadcoll);

  //:::::: Take T2 tracks ::::::
  Handle<T1T2TrackCollection> trackCollection;
  iEvent.getByLabel(TrackLabel,"T2TrackColl",trackCollection);
  












  // +++++++++++++++++++++
  // Begin of test code 
  // +++++++++++++++++++++

  T2ROGeometry geo; 
  std::string infostr;
  double padr1=0.;
  double padr2=0.; 
  double strr=0.;
  double Mypadr1=0.;
  double Mypadr2=0.; 
  double Mystrr1=0.;
  double Mystrr2=0.;
  
  for(unsigned int symb=0;symb<10;symb++)
    {
      
      T2GeometryUtil conv;
      T2GeometryUtil::T2DetInfo planeinfo;
      planeinfo=conv.GetT2Info(symb);
      
      geo=T2ROGeometry(planeinfo.cmsswid);
     
      infostr=conv.DetInfoFromRawId(planeinfo.cmsswid);
      int row=3;
      padr1=geo.GetPadRMin(row,0);
      padr2=geo.GetPadRMax(row,0);
     
      if(symb==0)
	{
	  Mypadr1=geo.GetPadRMin(4,0);
	  Mypadr2=geo.GetPadRMax(4,0);
	  Mystrr1=geo.GetStripRMin(18,0);
	  Mystrr2=geo.GetStripRMax(18,0);
	  if(PrintOnlyFinalOutput==false)
	    std::cout<<"Detector 0 has Strip(18,0) Min-Max= "<<Mystrr1<<"-"<<Mystrr2<<" and  Pad(5,0) Min-Max= "<<Mypadr1<<"-"<<Mypadr2<<std::endl;
	}

 
      strr=-1.;
      int introwstrip=0;
      while (((strr<padr2)&&(strr>padr1))==false)
	{
	  strr=geo.GetStripRMin(introwstrip,0);
	  introwstrip++;
	}

       /*
      T2GeometryUtil::vfatid_channel  vfchstr= conv.StripVfatsIdFromRowCol(introwstrip, 0, planeinfo.cmsswid);
      T2GeometryUtil::vfatid_channel  vfchpad= conv.PadVfatsIdFromRowCol(row, 0, planeinfo.cmsswid);
     
      std::cout<<"The overlap for pad row=|"<<row<<"| is given by strip row=|"<<introwstrip<<"|.  Strip R="<<strr<<" Pad R min-max="<<padr1<<" "<<padr2<<std::endl;

      std::cout<<"Strip Vfat Number="<<vfchstr.vfatiid<<".  Channel="<<vfchstr.channel<<std::endl;
      std::cout<<"Pad Vfat Number="<<vfchpad.vfatiid<<".  Channel="<<vfchpad.channel<<std::endl;
      std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
      std::cout<<" "<<std::endl;
      */
      if(PrintOnlyFinalOutput==false)
	{
	  std::cout<<" "<<std::endl;
	  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
	  std::cout<<"Plane Info: "<<infostr.c_str()<<std::endl;
	}

      T1T2TrackCollection::const_iterator TrkCit;
      for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++)
	{
	   for(unsigned int l=0; l<(*TrkCit).GetHitEntries();l++)
	     if((*TrkCit).GetHitT2(l).GetHitDetRawId()==planeinfo.cmsswid)
	       {
		 if(PrintOnlyFinalOutput==false)
		   std::cout<<"Hit R-Phi: "<< ((*TrkCit).GetHitT2(l)).GetHitR()<<"   "<<((*TrkCit).GetHitT2(l)).GetHitPhi()<<std::endl;
		 
		
		 std::vector<T2GeometryUtil::vfatid_channel> vectStripvf=conv.StripVfatsIdsFromStripVect((*TrkCit).GetHitT2(l)); 
		 std::vector<cluster_entry> entriesstripcl= ((*TrkCit).GetHitT2(l)).ClusterStrip_entries;
		 //	 std::cout<<"Strip (row,col)="<<std::endl;
		 for(unsigned int iu=0;iu<entriesstripcl.size();iu++)
		    {
		      entriesstripcl.at(iu);


		      //   std::cout<<"("<<row<<" , "<<col<<")";
		    }

		 if(PrintOnlyFinalOutput==false)
		   {
		     std::cout<<" "<<std::endl;
		 
		     for(unsigned g=0;g<vectStripvf.size();g++)
		       {
			 std::cout<<"Strip Vfat Number="<<vectStripvf.at(g).vfatiid<<".  Channel="<<vectStripvf.at(g).channel<<std::endl;
		       }
		   }
		 
		 std::vector<T2GeometryUtil::vfatid_channel> vectPadvf=conv.PadVfatsIdsFromPadVect((*TrkCit).GetHitT2(l)); 
		 std::vector<cluster_entry> entriespadcl= ((*TrkCit).GetHitT2(l)).ClusterPad_entries;
		 //std::cout<<"Pad (row,col)="<<std::endl;
		 for(unsigned int iu=0;iu<entriespadcl.size();iu++)
		    {
		      entriespadcl.at(iu);


		      //   std::cout<<"("<<row<<" , "<<col<<")";
		    }

		 if(PrintOnlyFinalOutput==false)
		   {
		     std::cout<<" "<<std::endl;

		     for(unsigned g=0;g<vectPadvf.size();g++)
		       {
			 std::cout<<"Pad Vfat Number="<<vectPadvf.at(g).vfatiid<<".  Channel="<<vectPadvf.at(g).channel<<std::endl;
		       }
		   }



	       }
	}
      if(PrintOnlyFinalOutput==false)
	{
	  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
	  std::cout<<" "<<std::endl; std::cout<<" "<<std::endl; std::cout<<" "<<std::endl; std::cout<<" "<<std::endl; std::cout<<" "<<std::endl;
	}

    }

  // +++++++++++++++++++++
  // End of test code 
  // +++++++++++++++++++++












  // ++++++++++++++++++++++++++++++++++
  // Begin of pattern generator code 
  // ++++++++++++++++++++++++++++++++++


  std::vector<std::map<unsigned int, T2GeometryUtil::vfatid_channel> >  Vect_Map_vfsymb_to_vfatidch;
  //24x65 pad
  for(unsigned col=0;col<=64;col+=5)//20 or 5
    {
      for(unsigned row=3;row<=23;row+=15)//25 or 12
	{
	   std::map<unsigned int, T2GeometryUtil::vfatid_channel>   Map_vfsymb_to_vfatidch;
	   if(PrintOnlyFinalOutput==false)
	     std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;

	   //reference (plane 0 variables)
	   double padrmin;double padrmax;double padphimin;
	   double strr=-1.;
	   int introwstrip=-1;
	   int stripcol=0;

	   for(unsigned int plsymb=0;plsymb<=39;plsymb++)
	    {
	      if(PrintOnlyFinalOutput==false){
	      std::cout<<"|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|-|"<<std::endl;
	      std::cout<<"Plane symb: "<<plsymb<<std::endl;
	      }

	      T2GeometryUtil conv;
	      T2GeometryUtil::T2DetInfo planeinfo;
	      padphimin = 0;
	      if(plsymb%10==0) //Get Reference Geometrical coordinate R-Phi on plane 0.
		{		  
		  //unsigned int symb=0;
		  if(PrintOnlyFinalOutput==false)
		    std::cout<<"Reference Pad row-col: "<<row<<"-"<<col<<std::endl;
		  
		  planeinfo=conv.GetT2Info(plsymb);
		  geo=T2ROGeometry(planeinfo.cmsswid);
		  infostr=conv.DetInfoFromRawId(planeinfo.cmsswid);
	   
	   
		   padrmin=geo.GetPadRMin(row,col);
		   padrmax=geo.GetPadRMax(row,col);
		   padphimin=geo.GetPadPhiMin(row,col);
		   geo.GetPadPhiMax(row,col);


		  //Find overlapping strip row;
		  strr=-1.;
		  introwstrip=-1;	   
		  while (((strr<padrmax)&&(strr>padrmin))==false)
		    {
		      introwstrip++;
		      strr=geo.GetStripRMin(introwstrip,0);	       
		    }
	   
		  //Find overlapping strip col;
		 
		  if((padphimin<=geo.GetStripPhiMax(0,0))&&(padphimin>=geo.GetStripPhiMin(0,0)))
		    {
		      stripcol=0;
		    }
		  else
		    stripcol=1;
		}
	   
	     
	     
	      planeinfo=conv.GetT2Info(plsymb);
	      geo=T2ROGeometry(planeinfo.cmsswid);
	      int mycolpad=0;
	      double actpadrmin=-1.;
	      double actpadrmax=-1.;double actpadphimin;double actpadphimax;
	       //Find the collinear pad of the reference pad in symb 0
	       
	      actpadrmin=geo.GetPadRMin(row,mycolpad);
	      actpadrmax=geo.GetPadRMax(row,mycolpad);
	      actpadphimin=geo.GetPadPhiMin(row,mycolpad);
	      actpadphimax=geo.GetPadPhiMax(row,mycolpad);

	      if(PrintOnlyFinalOutput==false){
	      std::cout<<"Before find: actpadphimin:"<<actpadphimin<<" padphimin:"<<padphimin<<" Actual pad r min-max:"<<actpadrmin<<"-"<<actpadrmax<<".  Reference strip r:"<<strr<<std::endl;
	      }
	      double dphi=fabs(actpadphimin-padphimin);
	     

	      bool stophere=( (dphi<0.1)  && (actpadrmin<=strr) && (actpadrmax>=strr) );
	      /*
	      if(stophere)
		{
		   if(dphi<0.1)		
		     std::cout<<"a";		  
	
		   if(actpadrmin<=strr)		    
		     std::cout<<"b";

		   if(actpadrmax>=strr)
		     std::cout<<"c";
		}
	      */
	       while(stophere==false)	      
		{
		  
		  //std::cout<<dphi<<" "<<actpadrmin<<" "<<strr<<" "<<actpadrmax<<std::endl;;
		  //std::cout<<actpadphimin<<padphimin<<std::endl;
		   // do
		   // myrowpad++;
		   //while(); 
		   mycolpad++;
		   actpadrmin=geo.GetPadRMin(row,mycolpad);
		   actpadrmax=geo.GetPadRMax(row,mycolpad);
		   actpadphimin=geo.GetPadPhiMin(row,mycolpad);
		   actpadphimax=geo.GetPadPhiMax(row,mycolpad);
		   dphi=fabs(actpadphimin-padphimin);
		   stophere=( (dphi<0.1)  && (actpadrmin<=strr) && (actpadrmax>=strr) );
		   
		}
		 
	      
	       //Find the corresponding overlappig strip
	       double actstrr=0.;
	       int introwactstrip=0;	   
	       actstrr=geo.GetStripRMin(introwactstrip,0);
	       while (((actstrr<actpadrmax)&&(actstrr>actpadrmin))==false)
		 {
		   introwactstrip++;
		   actstrr=geo.GetStripRMin(introwactstrip,0);		   
		 }

	       int actstripcol=0;
	       
	       //To put a condition compatible with phi conv.
	       double actpadphicentre=(actpadphimin+actpadphimax)/2.;
	       //double dphi1=(fabs(actpadphimin-geo.GetStripPhiMax(0,0)));
	       double dphi1=(fabs(actpadphicentre-geo.GetStripPhiMax(0,0)));
	       if (dphi1>180)
		 dphi1=360.-dphi1;
	       // double dphi2=(fabs(actpadphimin-geo.GetStripPhiMin(0,0)));
	        double dphi2=(fabs(actpadphicentre-geo.GetStripPhiMin(0,0)));
		if(PrintOnlyFinalOutput==false){
		std::cout<<"Plane pad Phi min: "<<actpadphimin<<" Phi centre: "<<actpadphicentre<<std::endl;
		std::cout<<" Strip(0,0) phimax:"<<geo.GetStripPhiMax(0,0)<<" Strip(0,0) phimin:"<<geo.GetStripPhiMin(0,0)<<std::endl;
		std::cout<<" Strip(0,1) phimax:"<<geo.GetStripPhiMax(0,1)<<" Strip(0,1) phimin:"<<geo.GetStripPhiMin(0,1)<<std::endl;
		std::cout<<"Strip distance from col0: "<<dphi1<<" "<<dphi2<<std::endl;
		}

	       if (dphi2>180)
		 dphi2=360.-dphi2;

	       if((dphi1<(192./2.))&&(dphi2<(192./2.)))
		 {		   
		   actstripcol=0;
		 }
	       else
		 actstripcol=1;        

	       if(PrintOnlyFinalOutput==false){
	       std::cout<<"Reference     pad r-c:"<<row<<" "<<col<<" strip r-c:"<<introwstrip<<"-"<<stripcol<<std::endl;
	       std::cout<<"Correspond to pad r-c:"<<row<<" "<<mycolpad<<" strip r-c:"<<introwactstrip<<"-"<<actstripcol<<std::endl;
	       }
	       //
	       double actpadrmin2=geo.GetPadRMin(row,mycolpad);
	       double actpadrmax2=geo.GetPadRMax(row,mycolpad);
	       double actstriprmax2=geo.GetStripRMin(introwactstrip,actstripcol);
	       double actstriprmin2=geo.GetStripRMax(introwactstrip,actstripcol);
	       if(PrintOnlyFinalOutput==false){
	       std::cout<<"After find: Actual pad r min-max:"<<actpadrmin2<<"-"<<actpadrmax2<<".  Strip min-max: "<<actstriprmin2<<"-"<<actstriprmax2<<". Pad r-c="<<row<<"-"<<mycolpad<<"  Strip r-c="<<introwactstrip<<"-"<<actstripcol<<std::endl;
	      //
	       }


	       T2GeometryUtil::vfatid_channel  vfchstr= conv.StripVfatsIdFromRowCol(introwactstrip, actstripcol, planeinfo.cmsswid);
	       T2GeometryUtil::vfatid_channel  vfchpad= conv.PadVfatsIdFromRowCol(row, mycolpad, planeinfo.cmsswid);     
	       if(PrintOnlyFinalOutput==false){
	       std::cout<<"Strip Vfat iid-channel : "<<vfchstr.vfatiid<<" - "<<vfchstr.channel<<std::endl;
	       std::cout<<"Pad Vfat iid-channel   : "<<vfchpad.vfatiid<<" - "<<vfchpad.channel<<std::endl;
	       }

	       unsigned int abssymbpad=plsymb*100+vfchpad.vfatiid;
	       Map_vfsymb_to_vfatidch.insert(pair<unsigned int,T2GeometryUtil::vfatid_channel>(abssymbpad,vfchpad));
	       unsigned int abssymbstrip=plsymb*100+vfchstr.vfatiid;
	       Map_vfsymb_to_vfatidch.insert(pair<unsigned int,T2GeometryUtil::vfatid_channel>(abssymbstrip,vfchstr));
	       if(PrintOnlyFinalOutput==false)
		 std::cout<<plsymb<<" "<<abssymbpad<<" "<<vfchpad.vfatiid<<" "<<abssymbstrip<<" "<<vfchstr.vfatiid<<" inserted in the map"<<std::endl;
	    }
	  Vect_Map_vfsymb_to_vfatidch.push_back(Map_vfsymb_to_vfatidch);
	  if(PrintOnlyFinalOutput==false){
	  std::cout<<"@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@"<<std::endl;
	  std::cout<<std::endl;
	  }
	}            
    }

  
  // ++++++++++++++++++++++++++++++++++
  // End of pattern generator code 
  // ++++++++++++++++++++++++++++++++++








  // +++++++++++++++++++++++++++++++++++++++++++
  // Begin of final output format generator code 
  // +++++++++++++++++++++++++++++++++++++++++++


  std::cout<<std::endl;    std::cout<<std::endl;
  std::cout<<"^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^"<<std::endl;


  std::map<int, std::vector<std::string> >  Mapfinaloutputformat;
  std::map<int, std::vector<std::string> >::iterator Mapfinaloutputformat_iterator;


  for (unsigned int trknumb=0;trknumb<Vect_Map_vfsymb_to_vfatidch.size();trknumb++)
    {

    std::map<unsigned int, T2GeometryUtil::vfatid_channel>   Map_vfsymb_to_vfatidch=Vect_Map_vfsymb_to_vfatidch.at(trknumb);
    /*
      std::map<unsigned int, T2GeometryUtil::vfatid_channel>   Map_vfsymb_to_vfatidch=Vect_Map_vfsymb_to_vfatidch.at(0);
      std::cout<<"Final Printing"<<std::endl;
      map<unsigned int, T2GeometryUtil::vfatid_channel>::const_iterator ditx;
      for(ditx= Map_vfsymb_to_vfatidch.begin();ditx!=Map_vfsymb_to_vfatidch.end();ditx++)
      {
      std::cout<<"plane: "<<(ditx->first)/100<<"      Plane GeoId: "<<ditx->second.vfatiid<<"("<<(ditx->first)%100<<")"<<"   Channel: "<<ditx->second.channel<<std::endl;
      }
      std::cout<<std::endl;    std::cout<<std::endl;
    */

    if(PrintOnlyFinalOutput==false){
      std::cout<<std::endl;    std::cout<<std::endl;  std::cout<<std::endl;
      std::cout<<"@@@@@@@@@@@@@ NEW Track Type:  4 tracks (one for each quarter) will be printed: @@@@@@@@@@@@@ "<<std::endl;
      std::cout<<std::endl;    std::cout<<std::endl;
    }

    for(unsigned int plane=0;plane<=39;plane++)
     { 
      if((plane%10)==0)
	{
	  if(PrintOnlyFinalOutput==false)
	  printf("\n\n\n******* *********** ***********   NEW QUARTER ********* ******** ******** \n\n");
	}

	  T2GeometryUtil converter;
	  T2GeometryUtil::T2DetInfo planeinfo2;
	  planeinfo2=converter.GetT2Info(plane);
	  if(PrintOnlyFinalOutput==false)
	    std::cout<<"(t2_detector_set="<<(plane/10)<<")"<<" plane Info: "<<converter.DetInfoFromRawId(planeinfo2.cmsswid).c_str()<<std::endl;
	  
      for(unsigned int VFPOS0_16=0;VFPOS0_16<17;VFPOS0_16++)
	{
	
	  unsigned int vfatsymb= plane * 100 + VFPOS0_16;
	  map<unsigned int, VFATRegisters>::const_iterator dit = Vfat_t2_mapping->readoutIdToRegisters.find(vfatsymb);
	
	  unsigned short printId = 0;
	  if (dit != Vfat_t2_mapping->readoutIdToRegisters.end()) {
	    
	    printId = dit->second.GetFullChipID();  
	    map<unsigned int, T2GeometryUtil::vfatid_channel>::const_iterator dit2 = Map_vfsymb_to_vfatidch.find(vfatsymb);
	    
	    if (dit2 != Map_vfsymb_to_vfatidch.end()) {

	      T2GeometryUtil::vfatid_channel thevfat=dit2->second;//Map_vfsymb_to_vfatidch(vfatsymb);
	      if(PrintOnlyFinalOutput==false)
		printf("Vfat Position:  ChipID=0x%x  VfatIId=%d  Channel=%d \n",printId,thevfat.vfatiid,thevfat.channel);
	    


	      char buffer[555];char buffer2[555];

	      int CCuAddrInt=0;
	      std::string Fec,CCuAddr,I2CChan; 
	      InfoForVfatController(&CCuAddrInt, &Fec,  &CCuAddr, &I2CChan, vfatsymb);
	      //nunused=sprintf(buffer,"detector_set: %d  HS: %d   Fec:%s   CCuAddr:%s   I2C:%s   VfatIId=%d  ChipID=0x%x  Channel=%d \n",(plane/10),(plane%10),Fec,CCuAddr,I2CChan,thevfat.vfatiid,printId,thevfat.channel);
	      // nunused=sprintf(buffer,"detector_set: %d  HS: %d      VfatIId=%d  ChipID=0x%x  Channel=%d \n",(plane/10),(plane%10),thevfat.vfatiid,printId,thevfat.channel);
	      
	      int thechanneltoPrint=thevfat.channel;
	      if(PrintChannelNumerationAsInVfatController)
		thechanneltoPrint++;

	      if(OnlyVfatANDchannel_Finalprinting==false)
		{
		  sprintf(buffer,"detector_set: %d  HS: %d      VfatIId=%d",(plane/10),(plane%10),thevfat.vfatiid);
		  sprintf(buffer2," ChipID=0x%x  Channel=%d",printId,thechanneltoPrint);
		}
	      else
		{

		  // nunused=sprintf(buffer2," ChipID=0x%x  Channel=%d",printId,thechanneltoPrint);
		  sprintf(buffer2," 0x%x,%d",printId,thechanneltoPrint);
		}

	      std::string Thenewstr1(buffer);
	      std::string Thenewstr2(buffer2);
	      std::string Thenewstr;

	      if(OnlyVfatANDchannel_Finalprinting==false)
		Thenewstr=Thenewstr1+"  Fec: "+Fec+"  CCuAddr: "+CCuAddr+"  I2C: "+I2CChan+"  "+Thenewstr2;
	      else
		Thenewstr=Thenewstr1+Thenewstr2;


	      Mapfinaloutputformat_iterator=Mapfinaloutputformat.find(CCuAddrInt);//vfatsymb
	      if (Mapfinaloutputformat_iterator != Mapfinaloutputformat.end())//Make a  bigger string vector 
		{
		  std::vector<string> initialVect=Mapfinaloutputformat_iterator->second;
		  //std::cout<<"!£"<<(Mapfinaloutputformat_iterator->second).size()<<"  "<<std::endl;
		  initialVect.push_back(Thenewstr);
		  
		  Mapfinaloutputformat_iterator->second=initialVect;
		  //std::cout<<" "<<Mapfinaloutputformat[vfatsymb].size()<<"  "<<std::endl;
		  if(PrintOnlyFinalOutput==false)
		  std::cout<<"VfatAbsID: "<<vfatsymb<<" with String: "<<Thenewstr.c_str()<<" inserted."<<std::endl;
		}
	      else //create a new vector and assign a string
		{
		  std::vector<string> initialVect;
		  initialVect.push_back(Thenewstr);
		  Mapfinaloutputformat.insert(pair<int,std::vector<string> >(CCuAddrInt,initialVect));//vfatsymb
		  if(PrintOnlyFinalOutput==false)
		  std::cout<<"VfatAbsID: "<<vfatsymb<<" with String: "<<Thenewstr.c_str()<<" created."<<std::endl;
		}
	      
	    }
	    else
	      {
		//	std::cout<<"Vfat with symbolic absolute id: "<<vfatsymb<<"not included in the simulated active vfats"<<std::endl;
		//std::cout<<""<<std::endl;
	      }
	  }
	  else
	    {
	      
	      std::cout<<"Vfat with symbolic absolute id: "<<vfatsymb<<"has not been parsed in XML file"<<std::endl;
	    }

	}
      std::cout<<" "<<std::endl;
    }


    } //next track printing
  
  
  std::map<int, std::vector<std::string> >::iterator  pos;
  for (pos = Mapfinaloutputformat.begin(); pos != Mapfinaloutputformat.end(); ++pos) {

    if(OnlyVfatANDchannel_Finalprinting==false)
      std::cout<<"----------------------------------------------------------------------"<<std::endl;
    int CCUint= pos->first;
    char outStr[256];
    
    sprintf(outStr,"%x",CCUint);

    if(OnlyVfatANDchannel_Finalprinting==false)
      printf("CCU 0x%s\n",outStr);

    //cout << "CCU: " << pos->first << std::endl;
    std::vector<std::string> vectstr=pos->second;
    for(unsigned int u=0;u<vectstr.size();u++)
      std::cout<<vectstr.at(u).c_str()<<std::endl;

    if(OnlyVfatANDchannel_Finalprinting==false)
      std::cout<<" "<<std::endl;
  }
  

  // +++++++++++++++++++++++++++++++++++++++++++
  // End of final output format generator code 
  // +++++++++++++++++++++++++++++++++++++++++++


      
  /*
for (unsigned int trknumb=0;trknumb<Vect_Map_vfsymb_to_vfatidch.size();trknumb++)
  {    
    std::map<unsigned int, T2GeometryUtil::vfatid_channel>   Map_vfsymb_to_vfatidch=Vect_Map_vfsymb_to_vfatidch.at(trknumb);
    for(unsigned int plane=0;plane<=39;plane++)
      { 
	for(unsigned int VFPOS0_16=0;VFPOS0_16<17;VFPOS0_16++)
	  {	    
	    unsigned int vfatsymb= plane * 100 + VFPOS0_16;
	    map<unsigned int, T2GeometryUtil::vfatid_channel>::const_iterator dit2 = Map_vfsymb_to_vfatidch.find(vfatsymb);
	    
	    if (dit2 != Map_vfsymb_to_vfatidch.end()) {
	      
	      T2GeometryUtil::vfatid_channel thevfat=dit2->second;//Map_vfsymb_to_vfatidch(vfatsymb);
	      printf("Vfat Position:  Abssymb%d  VfatIId=%d  Channel=%d \n",vfatsymb,thevfat.vfatiid,thevfat.channel);

	      int CCuAddrInt=0;
	      std::string Fec,CCuAddr,I2CChan;
 
	      InfoForVfatController(&CCuAddrInt, &Fec,  &CCuAddr, &I2CChan, vfatsymb);
	      
	    }
	  }
      }
  }
  */


  // ++++++++++++++++++++++++++++++++++++++++++++++
  // Begin of Validation code on the vfats pattern
  // ++++++++++++++++++++++++++++++++++++++++++++++

  if(DoClusterFromVFatChannel)
    {
     

      //Now Start from Vfats channel and Look Reconstruct the PAD-STRIP Digi.
      std::auto_ptr<T2StripDigiCollection> theStripDigis(new T2StripDigiCollection());
      std::auto_ptr<T2PadDigiCollection> thePadDigis(new T2PadDigiCollection());
      Int_t VFatFlagSP; Int_t convertedIndex;
      int iconv;
      std::auto_ptr<T2PadDigi> t2p;
      std::auto_ptr<T2StripDigi> t2s;
	
      T2PadDigi t2pv;
      T2StripDigi t2sv;
      Bool_t flagconv=true;
      T2Mapping VFATconvertDigichannell;
  
  
      for (unsigned int trknumb=0;trknumb<Vect_Map_vfsymb_to_vfatidch.size();trknumb++)
	{

	  std::map<unsigned int, T2GeometryUtil::vfatid_channel> Map_vfsymb_to_vfatidch=Vect_Map_vfsymb_to_vfatidch.at(trknumb);
	  map<unsigned int, T2GeometryUtil::vfatid_channel>::const_iterator ditx;

	  for(ditx= Map_vfsymb_to_vfatidch.begin();ditx!=Map_vfsymb_to_vfatidch.end();ditx++)
	    {    
	      unsigned int SymID= ((ditx->first)%100); //Vfat Gius. conv Id
	      //Nota: convert to Geo vuole indici a partire da 0, come nei MONITOR-daTA 
	      flagconv=VFATconvertDigichannell.convertToGeo((Int_t)SymID, (Int_t)(ditx->second.channel), VFatFlagSP, convertedIndex);
	      int actchj=ditx->second.channel;
	      T2GeometryUtil convs;
	      T2GeometryUtil::T2DetInfo planeinfos;
	      planeinfos=convs.GetT2Info(((ditx->first)/100));
	      T2DetId currentDet(planeinfos.cmsswid);
	      std::cout<<"Digi Process plane: "<<planeinfos.symb<<std::endl;
	      int row_,col_;
	      if(flagconv)    
		{
		  if(SymID<=16) //Is a data Vfat, Fill DIGI
		    {
		      if (VFatFlagSP == 1) {                  //enum fStrip,fPad
			if (SymID >= 2 && SymID <= 14 && actchj >= 4 && actchj <= 123) 
			  {			   
			    col_ = (convertedIndex-1)/24;
			    row_ = (convertedIndex-1)%24;
			    iconv=convertedIndex;		
			    t2p=std::auto_ptr<T2PadDigi>(new T2PadDigi ((int)iconv,(int)row_,(int)col_,0));   //third parameter=adc, set to 0.
			    t2pv= *t2p;		
			    thePadDigis->insertDigi(currentDet,t2pv);//insert(std::make_pair(currentDet,t2p(iconv,(int)row_,(int)col_,0)));  
			  }
			else
			  std::cout<<"Warning: VFat-PAD non-data channel on "<<std::endl;
		      }   
		  
		      if (VFatFlagSP == 0) {
			if (SymID == 0) 
			  {			  
			    row_=convertedIndex%256;
			    col_= 0;
			    iconv=convertedIndex;			
			    t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row_,(int)col_,0));
			    t2sv=*t2s;
			    theStripDigis->insertDigi(currentDet,t2sv); 
			  }
			if (SymID == 1) 
			  {			 
			    row_=convertedIndex%256;
			    col_= 0;
			    iconv=convertedIndex;			
			    t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row_,(int)col_,0));
			    t2sv=*t2s;
			    theStripDigis->insertDigi(currentDet,t2sv); 
			  } 
			if (SymID == 15)
			  {
			    
			    row_=convertedIndex%256;
			    col_= 1;
			    iconv=convertedIndex;			
			    t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row_,(int)col_,0));
			    t2sv=*t2s;
			    theStripDigis->insertDigi(currentDet,t2sv); 
			  }
			if (SymID == 16) 
			  {
			    row_=convertedIndex%256;
			    col_= 1;
			    iconv=convertedIndex;			
			    t2s=std::auto_ptr<T2StripDigi>(new T2StripDigi ((int)iconv,(int)row_,(int)col_,0));
			    t2sv=*t2s;
			    theStripDigis->insertDigi(currentDet,t2sv); 
			  }
	    
			if((SymID!=0)&&(SymID!=1)&&(SymID!=15)&&(SymID!=16))
			  std::cout<<"Warning: VFat-STRIP non-data channel on "<<std::endl;
	    
		      }
		    }
		  else //Trigger VFAT
		    {
		      std::cout<<"Warning: non-data VFAT with ID "<<SymID<<std::endl;
		    }
		  
		}

	    }
	}

 

      theStripClusters= auto_ptr<T2StripClusterCollection>(new T2StripClusterCollection());
      thePadClusters=auto_ptr<T2PadClusterCollection>( new T2PadClusterCollection());
      MakeEvCluster( thePadDigis,theStripDigis);


 
      std::cout<<"Strip Cluster:"<<std::endl;

      for(T2StripClusterCollection::const_iterator itstrip = theStripClusters->begin(); itstrip != theStripClusters->end(); itstrip++){
	vector<T2Cluster> stripClv = itstrip->second;
	for(unsigned int k=0;k<stripClv.size();k++){
	  clusterstripentries->Fill(stripClv[k].GetNoOfEntries());
     
	  T2GeometryUtil convs;
	  T2GeometryUtil::T2DetInfo planeinfos;
	  planeinfos=convs.GetT2Info(stripClv[k].GetDetID());
      
	  std::cout<<convs.DetInfoFromRawId(planeinfos.cmsswid).c_str()<<"  CLU-Phi:"<<(stripClv[k].GetClusterPhi())<<"    CLU-R:"<<(stripClv[k].GetClusterR())<<std::endl;
	}
      }

      std::cout<<"Pad Cluster:"<<std::endl;
      for(T2PadClusterCollection::const_iterator itpad = thePadClusters->begin(); itpad != thePadClusters->end(); itpad++){

	
	std::cout<<"Arm:"<<itpad->first.arm()<<"  Half:"<<itpad->first.halfTelescope()<<"  Plane:"<<itpad->first.plane()<<"  PlaneSide:"<<itpad->first.planeSide()<<std::endl;

	vector<T2Cluster> padClv = itpad->second;
	T2GeometryUtil conver_;
	T2GeometryUtil::T2DetInfo planeinfos_;

	for(unsigned int k=0;k<padClv.size();k++){
	  clusterpadentries->Fill(padClv[k].GetNoOfEntries());
	  
	  planeinfos_=conver_.GetT2Info(itpad->first.rawId());
	  // std::cout<<"$"<<itpad->first.rawId()<<" "<<padClv[k].GetDetID()<<" "<< planeinfos_.cmsswid<<std::endl;
	  std::string descr=conver_.DetInfoFromRawId(planeinfos_.cmsswid);
	 
	  std::cout<<descr.c_str()<<"  CLU-Phi:"<<(padClv[k].GetClusterPhi())<<"    CLU-R:"<<(padClv[k].GetClusterR())<<std::endl;
	}
      }
   



    }

  //std::cout<<"Detector 0 has Strip(18,0) Min-Max= "<<Mystrr1<<"-"<<Mystrr2<<" and  Pad(4,0) Min-Max= "<<Mypadr1<<"-"<<Mypadr2<<std::endl;


   // ++++++++++++++++++++++++++++++++++++++++++++++
   // End of Validation code on the vfats pattern
   // ++++++++++++++++++++++++++++++++++++++++++++++
































































































  //auto_ptr<T2PadClusterCollection>
  //Trketaetacut = std::auto_ptr<TH1F>(new TH1F("Trketaetacut","Reconstructed Track #eta (#eta cut)",1800,-9.0,9.0));




  /*
 
      std::cout<<"   Big test: "<<std::endl;
 int CCuAddrInt=0;
 std::string Fec,CCuAddr,I2CChan;
 int absvfatiid=1310;
  for(unsigned int j=0;j<40;j++) 
    for(unsigned int i=0;i<=16;i++)
      {
       absvfatiid=j*100+i;
       InfoForVfatController(&CCuAddrInt, &Fec,  &CCuAddr, &I2CChan, absvfatiid);
       std::cout<<"Plane "<<j<<"  Vfat "<< i <<" Fec "<<Fec<<" Ccuaddress: "<<CCuAddr.c_str()<<" converted to "<<CCuAddrInt<<"  I2C="<<I2CChan.c_str()<<std::endl;
     }
  */





 


















  //Handle< RawEvent >  input;
  // iEvent.getByType(input);
  //std::cout << " check --------- B " <<  std::endl;

  //--***************************************************         Reco TRACK        ************************************************************--/
  //--***************************************************        Reco TRACK        ************************************************************--/
  //--***************************************************        Reco TRACK        ************************************************************--/
  //Minimum hit distance
  /*
  T1T2TrackCollection::const_iterator TrkCit;
  double trketa=0.;
  double trkphi=0.;
  unsigned int trackcounter=0;
  unsigned int trackcountergood=0;
  unsigned int trackcounterchi=0;
  unsigned int trackcounterDZ=0;
  unsigned int trackcountereta=0;
  double DZgood;
  double numrecotrackcutmag0=0.;
  double chiRProb=0.;
  double chiPhiProb=0.;
  bool chi2condition;
  bool Rmincondition;
  bool Multcondition;
  double phitraccia;

  double polarangle;
  double lastParticleeta=0.;
  double chp12mag0=0.;
  double chp12=0.;
  double nump04575=0.;
  double numchout=0.;
  unsigned int numpmeno5070=0;
  
  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    lastParticleeta=(*p)->momentum().eta();
    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()<6.5)&&((*p)->momentum().eta()>5.3))
      {
	chp12mag0++;
	ChEnergyinT2->Fill(((*p)->momentum().e()));
      }

    if((*p)->momentum().eta()>4.5)
      if((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&((*p)->momentum().eta()>6.5)||((*p)->momentum().eta()<5.3))
	{
	  numchout++;

	  ChOutEtaEnergy->Fill(((*p)->momentum().eta()),((*p)->momentum().e()));;
	}

    if(((*p)->pdg_id()==111)&&((*p)->momentum().eta()<7.5)&&((*p)->momentum().eta()>4.5))
      {
	
	nump04575++;
	energyP04575->Fill((*p)->momentum().e());
	etaP04575->Fill((*p)->momentum().eta());
	if(((*p)->momentum().e())>1.)
	  P0EtaEnergy->Fill(((*p)->momentum().eta()),((*p)->momentum().e()));
	P0EtaEnergycorr->Fill(((*p)->momentum().eta()),((*p)->momentum().e()));
	
	
      }

    if(((*p)->pdg_id()==-211)&&((*p)->momentum().eta()<7.0)&&((*p)->status()==1)&&((*p)->momentum().eta()>5.0))
      {
	numpmeno5070++;
      }


    if((*p)->pdg_id()==-13)
       muEtaEnergy->Fill(((*p)->momentum().eta()),((*p)->momentum().e()));

    if((*p)->status()==1)
      stablepdg->Fill((*p)->pdg_id());
    if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&(fabs((*p)->momentum().eta()<6.5))&&(fabs((*p)->momentum().eta()>5.3)))
      {chp12++;}
  }



  for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
    if(numpmeno5070<15)
      PimenoEtaEnergy[numpmeno5070]->Fill(((*p)->momentum().eta()),((*p)->momentum().e()));
  }
  
  bool trkinH0;
  bool trkinH1;
unsigned int  roadentries;
roadentries=(*TrkCit).GetHitEntries();
 if(trackCollection->size()==0)
 std::cout<<"size: "<<trackCollection->size()<<std::endl;

  for(TrkCit=trackCollection->begin(); TrkCit!=trackCollection->end(); TrkCit++){
    chiRdistro->Fill((*TrkCit).ChiSquaredR());
    chiPhidistro->Fill((*TrkCit).ChiSquaredPhi());
    
    polarangle=(*TrkCit).GetTy();
    
    trketa= (*TrkCit).Eta();
    trkphi=(*TrkCit).Phi()*180/3.14159265;
    Trketa->Fill(trketa);
    Trkphi->Fill(trkphi);
    trkinH0=false;
    trkinH1=false;

    if((trkphi<80)||(trkphi>280))
      trkinH0=true;

    if((trkphi>100)&&(trkphi<260))
      trkinH1=true;

    if(trkinH0)
      {
	TrkEtaH0->Fill(trketa);
      }  
    if(trkinH1)
      {
	TrkEtaH1->Fill(trketa);  
      } 

    
    chiRProb=TMath::Prob((*TrkCit).ChiSquaredR(),((*TrkCit).GetHitEntries()-2));
    chiPhiProb=TMath::Prob((*TrkCit).ChiSquaredPhi(),((*TrkCit).GetHitEntries()-1));  
    Chi2RProb->Fill(chiRProb);
    Chi2PhiProb->Fill(chiPhiProb);
    phitraccia= (*TrkCit).Phi()*180./3.14159265;
    if(phitraccia>90)	
      TrKPhi0degstudyg->Fill(phitraccia-360.);	 
    else
      TrKPhi0degstudyg->Fill(phitraccia);


     MeanEtaResvsPhi->Fill(trkphi,trketa); 
     if(trketa>0)
       MeanEtaResvsPhiHem1->Fill(trkphi,trketa);
     if(trketa<0)
      MeanEtaResvsPhiHem2->Fill(trkphi,trketa);
      
      
    
    if((chiPhiProb<PhiChiProbCut)&&(chiRProb<RChiProbCut))  
      chi2condition=false;
    else
      {
	Trketachicut->Fill(trketa);
	Z0distro->Fill((*TrkCit).Z_at_Rmin());
	R0distro->Fill((*TrkCit).GetTx());
	chi2condition=true;
	trackcounterchi++;
      }	
    
    Rmincondition=false;
    if((*TrkCit).GetTx()<25.0)
      Rmincondition=true;
	

    Multcondition=true;
    NumhitinTrack->Fill((*TrkCit).GetHitEntries());
    
    
    	
    DZgood=1000.;
    if (energy==10.)
      DZgood=sigmaZfunctE10(fabs(trketa));
    
    if (energy==30.)
      DZgood=sigmaZfunctE30(fabs(trketa));

    if (energy==50.)
      DZgood=sigmaZfunctE50(fabs(trketa));
    
    DZgood=DZScale*DZgood;
   
       

    if(fabs((*TrkCit).Z_at_Rmin())<DZgood)
      {
	trackcounterDZ++;
	Trketadzcut->Fill(trketa);
      }
    
    if(IsinT2(trketa,tracketamin,tracketamax)==true)
      {
	trackcountereta++;
	Trketaetacut->Fill(trketa);
	if(chi2condition==true)
	  Trketaetachicut->Fill(trketa);
      }

    
    if((IsinT2(trketa,tracketamin,tracketamax)==true)&&(fabs((*TrkCit).Z_at_Rmin())<DZgood))
      Trketadzetacut->Fill(trketa);
      
 

    if((Multcondition)&&(Rmincondition)&&(chi2condition)&&(IsinT2(trketa,tracketamin,tracketamax))&&(fabs((*TrkCit).Z_at_Rmin())<DZgood)){
      trackcountergood++;
      
      Trketagood->Fill(trketa);
      Trkphigood->Fill(trkphi);

      trkinH0=false;
      trkinH1=false;

      if((trkphi<80)||(trkphi>280))
	trkinH0=true;

      if((trkphi>100)&&(trkphi<260))
	trkinH1=true;

      if(trkinH0)
	{
	  TrkEtaGoodH0->Fill(trketa);
	  
	}  
      if(trkinH1)
	{
	  TrkEtaGoodH1->Fill(trketa);  
	}



      tantheta->Fill((*TrkCit).GetTy());
      if(trackcountergood==1)
	 if((trkphi<80)||(trkphi>280))
	   NumhitinTrackGood->Fill((*TrkCit).GetHitEntries());

      MeanEtaResvsPhiCut->Fill(trkphi,trketa); 
      if(trketa>0)
	MeanEtaResvsPhiHem1Cut->Fill(trkphi,trketa);
      if(trketa<0)
	MeanEtaResvsPhiHem2Cut->Fill(trkphi,trketa);

         

      if((trketa>0)&&(polarangle>2.0*atan(exp(-tracketamax)))) //ulteriore taglio per eta mag 0
	{
	  tanthetam0->Fill((*TrkCit).GetTy());
	  numrecotrackcutmag0= numrecotrackcutmag0+1;
	}
      Trketadzetachicut->Fill(trketa);

      if(((*TrkCit).GetHitEntries()<10)&&((*TrkCit).GetHitEntries()>5))
	{
	  for(unsigned int l=0; l<(*TrkCit).GetHitEntries();l++)
		{
		  RecoZHit->Fill((*TrkCit).GetHitT2(l).GetHitZ());
		}
	}

    }
    trackcounter++;
  }
  
  HNumrecotrackcutmag0->Fill(numrecotrackcutmag0);
  HTrackcounter->Fill(trackcounter);  
  
*/
  
  /*

for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
  if(((*p)->pdg_id()==111)&&((*p)->momentum().eta()<7.5)&&((*p)->momentum().eta()>4.5))
      {
	P0NumEtacorr->Fill(nump04575,(*p)->momentum().eta());
	P0NumEnergycorr->Fill(nump04575,(*p)->momentum().e());
      }
}


  
  if((fabs(lastParticleeta)<6.4)&&(fabs(lastParticleeta)>5.7))
    {
      HTrackcounterNonCrit->Fill(trackcounter);      
    }

  if((fabs(lastParticleeta)<5.6)&&(fabs(lastParticleeta)>5.5))
    {
      HTrackcounterCrit->Fill(trackcounter);
    }
  
  chPartinT2->Fill(chp12);

  

  SingleParticleEfficiency->Fill(lastParticleeta,trackcounter);
  SingleParticleEfficiencyAllCuts->Fill(lastParticleeta,trackcountergood);
  SingleParticleEfficiencyCutChi->Fill(lastParticleeta,trackcounterchi);
  SingleParticleEfficiencyCutEta->Fill(lastParticleeta,trackcountereta);  
  SingleParticleEfficiencyCutDZ->Fill(lastParticleeta,trackcounterDZ);  

  if((trackcounter>0)&&(chp12==0))
    SourceofSecondary->Fill(lastParticleeta);

  if((trackcounter>0))
    SourceofReco->Fill(lastParticleeta);
 

  //std::cout << " check --------- C " <<  std::endl;

 
  double SingleTrackEvent_EtaRec=0;
  double PI = 3.1415927;
  double SingleTrackEvent_PhiRec=0;

  //get G4 tracks
  vector< pair<double,double> > PrimarySimTracks;

  Handle<SimTrackContainer> G4TrkContainer;
  iEvent.getByType(G4TrkContainer);
  if (!G4TrkContainer.isValid()) {
    LogError("TrackerHitAnalyzer::analyze") << "Unable to find SimTrack in event!";
    return;
  }
  

  if(singleparticle)
    {
  for (SimTrackContainer::const_iterator itTrk = G4TrkContainer->begin(); itTrk != G4TrkContainer->end(); ++itTrk) {
    double eta =0, phi =0, p =0;
    const HepLorentzVector G4Trk(itTrk->momentum().x(),itTrk->momentum().y(),itTrk->momentum().z(), itTrk->momentum().e() ) ;
    p = sqrt(G4Trk[0]*G4Trk[0]+G4Trk[1]*G4Trk[1]+G4Trk[2]*G4Trk[2]);
    if ( p == 0){
      LogError("TrackerHitAnalyzer::analyze") << "TrackerTest::INFO: Primary has p = 0 ";
    } else {
      double costheta  = G4Trk[2]/p;
      double theta = acos(TMath::Min(TMath::Max(costheta, -1.),1.));
      eta = -log(tan(theta/2));
      if ( G4Trk[0] != 0 || G4Trk[1] != 0)
        phi = atan2(G4Trk[1],G4Trk[0]);

      if(phi<0)
        phi = 2*PI + phi;
      phi=phi*180/PI;

      pair<double,double> *myprimarytrack = new pair<double,double>(eta,phi);
      PrimarySimTracks.push_back(*myprimarytrack);
    }
  }

  // comparison reconstructed track  - GEANT track

  for(T1T2TrackCollection::const_iterator TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++)
    {
      DZgood=1000.;
      if (energy==10.){
        DZgood=sigmaZfunctE10(fabs((*TC_it).Eta()));
      }
      if (energy==50.) {
        DZgood=sigmaZfunctE50(fabs((*TC_it).Eta()));
      }
      DZgood=DZgood*DZScale;

      SingleTrackEvent_EtaRec = (*TC_it).Eta();
      SingleTrackEvent_PhiRec = (*TC_it).Phi()*180.0/3.14159265;

      chiRProb=TMath::Prob((*TC_it).ChiSquaredR(),((*TC_it).GetHitEntries()-2));
      chiPhiProb=TMath::Prob((*TC_it).ChiSquaredPhi(),((*TC_it).GetHitEntries()-1));

   

      for(vector< pair<double,double> >::iterator pair_it=PrimarySimTracks.begin(); pair_it!=PrimarySimTracks.end(); pair_it++){
        double DE_temp = SingleTrackEvent_EtaRec - (*pair_it).first;
        double DF_temp = SingleTrackEvent_PhiRec - (*pair_it).second;

        //std::cout << "Track DE, DF = " << DE_temp<< " , " << DF_temp << std::endl;
        //std::cout << "Phi reco: "<< SingleTrackEvent_PhiRec << "Phi G4: "  << (*pair_it).second << std::endl;
	 if((chiPhiProb>PhiChiProbCut)&&(chiRProb>RChiProbCut)&&(trackcounter==1))
	   {
	     DPhiChiCutOneTrk->Fill(DF_temp);
	     DEtaChiCutOneTrk->Fill(DE_temp);
	   }   
      
	if((chiPhiProb<PhiChiProbCut)&&(chiRProb<RChiProbCut))  
	  chi2condition=false;
	else
	  {
	    chi2condition=true;	    
	  } 
        if((chi2condition)&&(IsinT2((*TC_it).Eta(),tracketamin,tracketamax))){
          double DZgood=1000.;

          if (energy==10.)
            DZgood=sigmaZfunctE10(fabs((*TC_it).Eta()));
          if (energy==50.)
            DZgood=sigmaZfunctE50(fabs((*TC_it).Eta()));
	  
	  DZgood=DZgood*DZScale;

          if(fabs((*TC_it).Z_at_Rmin())<DZgood){
            DEtaGoodTrk->Fill(DE_temp);
            DPhiGoodTrk->Fill(DF_temp);
          }

        }
      }
    }
  }
  //std::cout << " check --------- D " <<  std::endl;

  //::::::::::Local variable Declaration:::::::::
  double lmyx=0.0;
  double lmyy=0.0;
  double lmyr=0.0;
  double phipart=0.0;

  T2DetId myT2Det;
  int myhalftele;
  int myplaneside;
  int myplane;
  double zdetshift;
  double zglobhit;

  //int zmm[2][10]; //2 row for halftelescope. 10 coloumn for planes
  double z1=13828.3;          //14035.605; //first Gem first drift gas zone (mm)
  double planedist= 86.0;
  double btbdist=24.6;   //25.0;
  double ovdist=43.0;
  double zinsidedet=1.5;           //4.5; From first Drift zone to RO board = 9mm
*/

  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //
  // **************************************************** STRIP STUDIES ************************************************ //

  /*
  for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
    vector<T2Cluster> stripClv = itstrip->second;
    for(unsigned int k=0;k<stripClv.size();k++){
      clusterstripentries->Fill(stripClv[k].GetNoOfEntries());
    }
  }
  */

  /*
if(trackcountergood==1)
  if(singleparticle)  
 for(T2StripClusterCollection::const_iterator itstrip = t2strclcoll->begin(); itstrip != t2strclcoll->end(); itstrip++){
      myT2Det = itstrip->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();
      vector<T2Cluster> stripClv = itstrip->second;
      for(unsigned int k=0;k<stripClv.size();k++){
    

      for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin(); hitMapItr != hitMap.end(); ++hitMapItr){

          const PSimHitContainer & planeSimHits = hitMapItr->second;
          T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

          if((theT2DetId->plane()==myplane)&&(theT2DetId->planeSide()==myplaneside)&&(planeSimHits.size()==1)&&(stripClv.size()==1)){
            for (unsigned int l=0; l<planeSimHits.size(); l++){

              PSimHit myplanehits= planeSimHits[l];
              lmyx= myplanehits.localPosition().x();
              lmyy= myplanehits.localPosition().y();
              lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);
	      diffRCluHit->Fill(stripClv[k].GetClusterR()-lmyr);
	      
	    }
	  }
      }

    }
  }
  */
 

  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
  // **************************************************** PAD STUDIES ************************************************ //
/*
 if(singleparticle)
  if(trackcountergood==1)                       //Comupte the differences only if there is no secondaries
    for(T2PadClusterCollection::const_iterator itpad = t2padclcoll->begin(); itpad != t2padclcoll->end(); itpad++){

      myT2Det = itpad->first;
      myhalftele= myT2Det.halfTelescope();
      myplaneside= myT2Det.planeSide();
      myplane= myT2Det.plane();

      vector<T2Cluster> padClv = itpad->second;

      for(unsigned int k=0;k<padClv.size();k++){

        clusterpadentries->Fill(padClv[k].GetNoOfEntries());

        // now we run for each plane the digi-simulation
        for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin(); hitMapItr != hitMap.end(); ++hitMapItr){

          const PSimHitContainer & planeSimHits = hitMapItr->second;
          T2DetId *theT2DetId =new T2DetId(hitMapItr->first);

          if((theT2DetId->plane()==myplane)&&(theT2DetId->planeSide()==myplaneside)&&(planeSimHits.size()==1)&&(padClv.size()==1)){
            for (unsigned int l=0; l<planeSimHits.size(); l++){

              PSimHit myplanehits= planeSimHits[l];
              lmyx= myplanehits.localPosition().x();
              lmyy= myplanehits.localPosition().y();
              lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);

              //strange GEANT x-y sign convention
              if((myplaneside==1)&&(myhalftele==0)){
                phipart=atan2(lmyy,lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==0)&&(myhalftele==0)){
                phipart=atan2(-lmyy,lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==1)&&(myhalftele==1)){
                phipart=atan2(-lmyy,-lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
              }
              if((myplaneside==0)&&(myhalftele==1)){
                phipart=atan2(lmyy,-lmyx);
                phipart=phipart*180.0/3.14159265;
                if (phipart<0)
                  phipart= 360.0 - fabs(phipart);
                //std::cout<<"PS0 HT1  Reco: "<<padClv[k].GetClusterPhi()<<"; Geant: "<<phipart<<";    X-Y: "<<lmyx<<"-"<<lmyy<<std::endl;
              }
              diffphiCluHit->Fill(padClv[k].GetClusterPhi()-phipart);
            }
          }
        }                 // RPhiHitGeant

        for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
          phipart=(*p)->momentum().phi()*180.0/3.14159265;

          if (phipart<0)
            phipart= 360.0 - fabs(phipart);  // From particleGun Phi to Clusterizator Phi
          diffphiCluGun->Fill(padClv[k].GetClusterPhi()-phipart);
        }
      }
    }

  double expsimuz;
 for(unsigned int pln=0;pln<5;pln++)
   for(unsigned int ht=0;ht<2;ht++)
     for(unsigned int pls=0;pls<2;pls++)
       {
	 expsimuz=13828.3+pln*86.0+pls*24.6+(1-ht)*43.0;
	 SimuZHit->Fill(expsimuz);
	 //std::cout << "Plane: "<<pln << "     Planse Side: "<<pls<<"      Half Telescope: "<<ht<<"   ---      Z= "<<expsimuz<<std::endl;
       }
*/
  //std::cout << " check --------- E " <<  std::endl;
  // **************************************************** HIT STUDIES *********************************************** //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //
  // **************************************************** HIT STUDIES ************************************************ //

  // COMPARISON RECOHIT LOCALHIT
  // now we run for each plane the digi-simulation
/*
if(singleparticle)
  if(trackcountergood==1)
    {
    unsigned int  symb;
    for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++){
      symb=(*ithit).GetHitPlane()*4+(*ithit).GetHitHalftele()*2+(*ithit).GetHitPlaneSide();
	
      //mydet.arm()*20+mydet.plane()*4+mydet.halfTelescope()*2+mydet.planeSide()

      if((*ithit).GetHitHalftele()==0)
	SymbIdHT0->Fill(symb);
    }
 
    for(map<int, PSimHitContainer>::const_iterator hitMapItr = hitMap.begin();  hitMapItr != hitMap.end(); ++hitMapItr){
      const PSimHitContainer & planeSimHits = hitMapItr->second;
      T2DetId *theT2DetId =new T2DetId(hitMapItr->first);
      
      myhalftele= theT2DetId->halfTelescope();
      myplaneside= theT2DetId->planeSide();
      myplane= theT2DetId->plane();
      zdetshift=myplane*planedist+myplaneside*btbdist;
      if (myhalftele==0)
        zdetshift=zdetshift+ovdist;
      zglobhit=zdetshift+z1+zinsidedet;

      for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++){
        if (fabs(zglobhit-(ithit->GetHitZ()))<0.5)
          for (unsigned int l=0; l<planeSimHits.size(); l++){
            PSimHit myplanehits= planeSimHits[l];
            lmyx= myplanehits.localPosition().x();
            lmyy= myplanehits.localPosition().y();  //prima entryPoint()
            lmyr= sqrt(lmyx*lmyx + lmyy*lmyy);

	    

            RLocHRecoH->Fill((ithit->GetHitR())-lmyr);
          }
      }
    }
    }
*/
  //Studies on Reco Hit









  /*
   for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
     
     if(((*p)->status()==1)&&(((*p)->pdg_id())!=2212))
       EtaParticleAll->Fill((*p)->momentum().eta());           

     if ((PartCharge((*p)->pdg_id())!=0)&&((*p)->status()==1)&&(fabs((*p)->momentum().eta())<6.5)&&(fabs((*p)->momentum().eta())>5.3))
       {
	 phipar=(*p)->momentum().phi()*180.0/3.14159265;
	 pareta=(*p)->momentum().eta();
	 
	 if (phipar<0)
          phipar= 360.0 + phipar;
	 
	 EtaParticle->Fill(pareta);     
	 if(pareta>0)
	   MeanEtaResvsPhiHem1Part->Fill(phipar,pareta);

	 
	 if(pareta<0)
	   MeanEtaResvsPhiHem2Part->Fill(phipar,pareta);
	 
	 MeanEtaResvsPhiPart->Fill(phipar,pareta);
       }
   }
  */
  /*
if(singleparticle)
  if(trackcountergood==1)
    for(T2HitCollection::const_iterator ithit = t2hitcoll->begin(); ithit != t2hitcoll->end(); ithit++) {
      hitr=ithit->GetHitR();
      hitphi=ithit->GetHitPhi();
      hitdr=ithit->GetHitDR();
      hitdphi=ithit->GetHitDPhi();
      hnumstrip=ithit->GetHitNumStrip();
      hnumpad=ithit->GetHitNumPad();

      for (GenEvent::particle_const_iterator p = evt->particles_begin(); p != evt->particles_end(); p++ ){
        phipar=(*p)->momentum().phi()*180.0/3.14159265;

        if (phipar<0)
          phipar= 360.0 + phipar;  // From particleGun Phi to Clusterizator Phi

        if ((fabs((*p)->momentum().eta())<6.5)&&(fabs((*p)->momentum().eta())>5.3)){
          if (fabs(phipar-hitphi)<25.0){
            //cout<<"fabs(phipar-hitphi)  ="<<fabs(phipar-hitphi)<<"   particle gen phi: "<<phipar<<"   Hit phi"<< hitphi<<endl;
            diffphiGUNHIT->Fill(phipar-hitphi);
          }
        }
      }
      
      double minZdist=10000.;
      bool oneDZfound=false;
      if((ithit->GetHitPlane()==1)&&(ithit->GetHitPlaneSide()==1)&&(ithit->GetHitHalftele()==0))
	{
	  for(T2HitCollection::const_iterator ithit3 = t2hitcoll->begin(); ithit3 != t2hitcoll->end(); ithit3++) {

	    if((ithit->GetHitZ()<ithit3->GetHitZ())&&(ithit3->GetHitHalftele()==0))
	      {
	      if((fabs(ithit3->GetHitZ()-ithit->GetHitZ())<minZdist)&&((ithit3->GetHitZ()-ithit->GetHitZ())>2.))
		{
		  if(((fabs(ithit3->GetHitPhi()-ithit->GetHitPhi())<10.)||(fabs(ithit3->GetHitPhi()-ithit->GetHitPhi())>350.))&&(fabs(ithit3->GetHitR()-ithit->GetHitR())<3.  ) )
		    {
		      minZdist=fabs(ithit3->GetHitZ()-ithit->GetHitZ());
		      oneDZfound=true;
		    }
		}
	      }
	  }
	}
      if(oneDZfound)
	MinDz2Plane->Fill(minZdist);
      
    }


*/
  /*
  if(trackcountergood==1)   
    if((lastParticleeta<6.5)&&(lastParticleeta>5.3))
      if((trkphi<80)||(trkphi>280))
	{
	  unsigned int  Hitentries=0;
	  for(T2HitCollection::const_iterator itss =t2hitcoll->begin(); itss != t2hitcoll->end(); itss++){
	    if(itss->GetHitClass()==1)
	      Hitentries++;
	  }
	  MeanT2HitCl1vsEta->Fill(lastParticleeta,Hitentries);
	  MeanT2HitvsEta->Fill(lastParticleeta,((*t2hitcoll).size()));
	  MeanT2GeantHitvsEta->Fill(lastParticleeta,(hitMap.size()));
	  MeanT2PadvsEta->Fill(lastParticleeta,((*t2padclcoll).size()));
	  MeanT2StripvsEta->Fill(lastParticleeta,((*t2strclcoll).size()));
	  unsigned int pc=0;
	  unsigned int sc=0;
	  for(itp= PadDigiptr->begin(); itp!=PadDigiptr->end(); ++itp)
	     {
	       pc++;
	     }
	  for(its= StripDigiptr->begin(); its!=StripDigiptr->end(); ++its)
	     {
	       sc++;
	     }
	  MeanT2PadDigivsEta->Fill(lastParticleeta,pc);
	  MeanT2StripDigivsEta->Fill(lastParticleeta,sc);
	  //itp= PadDigiptr->begin()
	  unsigned int  roadentries=0;
	  for(T2RoadCollection::const_iterator itss = t2roadcoll->begin(); itss != t2roadcoll->end(); itss++){
	  
	    roadentries=((*itss).thisRoad).size();
	    if(roadentries>=3)
	    MeanT2RoadEntriesvsEta->Fill(lastParticleeta,roadentries);
	    //roadentries++;
	  }


	 
	}

NumP0VsNumTracks->Fill(nump04575,((double)numrecotrackcutmag0));
NumChPartVsNumP0->Fill(((double)chp12mag0),nump04575);
NumChPartInVsNumChPartOut->Fill(((double)chp12mag0),numchout);



if(chp12mag0>0)
  {
 MultiParticleEfficiencyCutNorm->Fill((double)chp12mag0,(double)(((double)numrecotrackcutmag0)/((double)chp12mag0)));
 //cout<<" eta>0 Ch -- Tr:   "<<chp12mag0 <<" --  "<<numrecotrackcutmag0<<endl;
 
}

 MultiParticleEfficiencyCut->Fill(((double)chp12mag0),((double)numrecotrackcutmag0));

  */




   
}






/*

//----------------------------------------------------------------\\
//----------------- Note about the histograms---------------------\\
//----------------------------------------------------------------\

The meaningfull histograms for one particle events are:

  clusterstripentries
  clusterpadentries
  diffphiGUNHIT
  DPhiGoodTrk
  DEtaGoodTrk
  RLocHRecoH
  diffRCluHit
  diffphiCluGun
  diffphiCluHit
  SingleParticleEfficiencyAllCuts
  SingleParticleEfficiency


The meaningfull histograms for one particle or multiparticle events are:
  Trketagood
  Trkphigood
  Trketa
  Trkphi
  MultiParticleEfficiencyCutNorm
  MultiParticleEfficiencyCut
  Chi2PhiProb
  Chi2RProb
  Chi2PhiProbLogy
  Chi2RProbLogy

*/



// ------------ method called once each job just before starting event loop  ------------
void T2GeometryAnalyzer::beginJob()
{

  std::string xmlfilename="/home/mirko/SL/WorkingArea/CMSSW_311SVN2/CMSSW_3_1_1/src/TotemT1T2Validation/T2GeometryValidation/test/T2GeoMapIP5_4quarter_vmea_cmssw.xml";
  DAQInformationSourceXML_a* onemap=new DAQInformationSourceXML_a(xmlfilename);
   boost::shared_ptr<DAQInformationT2_a> testmap;
  testmap=onemap->produceMap();

  for(unsigned int VFPOS0_16=0;VFPOS0_16<17;VFPOS0_16++)
    {
      
      unsigned int vfatsymb= /*plane * 100 +*/ VFPOS0_16;
      map<unsigned int, VFATRegisters>::const_iterator dit = testmap->readoutIdToRegisters.find(vfatsymb);
	
      unsigned short printId = 0;
      if (dit != testmap->readoutIdToRegisters.end()) {
	
	printId = dit->second.GetFullChipID();  
	//	map<unsigned int, T2GeometryUtil::vfatid_channel>::const_iterator dit2 = Map_vfsymb_to_vfatidch.find(vfatsymb);
	    
	//if (dit2 != Map_vfsymb_to_vfatidch.end()) {
	  
	// T2GeometryUtil::vfatid_channel thevfat=dit2->second;//Map_vfsymb_to_vfatidch(vfatsymb);
	  
	printf("Class test: Vfat Position:  ChipID=0x%x  vfatsymb=%d   \n",printId,vfatsymb);
	  //}
      }
    }



  TH1::AddDirectory(kFALSE);
 
  R0distro=std::auto_ptr<TH1F>(new TH1F("R0distro","Track r0 only chi2 cut",32,-100,100));
  R0distro->SetDirectory(0);			       
  Z0distro=std::auto_ptr<TH1F>(new TH1F("Z0distro","Track Z at r min only chi2 cut",50,-8000,8000));
  Z0distro->SetDirectory(0);
  NumhitinTrack=std::auto_ptr<TH1F>(new TH1F("NumhitinTrack","# hit in Track",18,-0.5,17.5));
  NumhitinTrack->SetDirectory(0);
  NumhitinTrackGood=std::auto_ptr<TH1F>(new TH1F("NumhitinTrackGood","# hit in Track",18,-0.5,17.5));
  NumhitinTrackGood->SetDirectory(0);
  
  chPartinT2=std::auto_ptr<TH1F>(new TH1F("chPartinT2","# Charged Particle in T2",16,-0.5,15.5));
  chPartinT2->SetDirectory(0);
  HNumrecotrackcutmag0=std::auto_ptr<TH1F>(new TH1F("HNumrecotrackcutmag0std","# Track reconstructed, #eta >0, all cut",16,-0.5,15.5));
  HNumrecotrackcutmag0->SetDirectory(0);
  HTrackcounter=std::auto_ptr<TH1F>(new TH1F("HTrackcounter","# Track reconstructed",26,-0.5,25.5));
  HTrackcounter->SetDirectory(0);
  HTrackcounterCrit=std::auto_ptr<TH1F>(new TH1F("HTrackcounterCrit","# Track reconstructed 5.5<#eta Ch< 5.6",10,-0.5,9.5));
  HTrackcounterCrit->SetDirectory(0);
  HTrackcounterNonCrit=std::auto_ptr<TH1F>(new TH1F("HTrackcounterNonCrit","# Track reconstructed 5.7<#eta Ch< 6.4",10,-0.5,9.5));     
  HTrackcounterNonCrit->SetDirectory(0);
   
  tantheta=std::auto_ptr<TH1F>(new TH1F("tantheta","tantheta, arz",100,-0.05,0.1));     
  tantheta->SetDirectory(0);

  tanthetam0=std::auto_ptr<TH1F>(new TH1F("tanthetam0","tantheta, arz",100,-0.05,0.1));     
  tanthetam0->SetDirectory(0);

  RecoZHit = std::auto_ptr<TH1F>(new TH1F("RecoZHit","Reco Hit Z",2000, 13500, 14500));
  RecoZHit->SetDirectory(0);
  
  SimuZHit = std::auto_ptr<TH1F>(new TH1F("SimuZHit","Simulated Hit Z",2000, 13500, 14500));
  SimuZHit->SetDirectory(0);
  // Chi2PhiProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2PhiProblog","Azimuthal #chi^{2} probability",400,400));
  // Chi2PhiProbLogy->SetDirectory(0);
  //Chi2RProbLogy= std::auto_ptr<TCanvas>(new TCanvas("Chi2RProblog","Radial #chi^{2} probability",400,400));
  //Chi2RProbLogy->SetDirectory(0);

  clusterstripentries = std::auto_ptr<TH1F>(new TH1F("clusterstripentries","Number of strips in strip-clusters",10, -0.5, 9.5));
  clusterstripentries->SetDirectory(0); 
  clusterpadentries = std::auto_ptr<TH1F>(new TH1F("clusterpadentries","Number of pads in pad-clusters",10, -0.5, 9.5));
  clusterpadentries->SetDirectory(0); 
  diffphiGUNHIT= std::auto_ptr<TH1F>(new TH1F("diffphiGUNHIT","#Phi Generated particle - #Phi Reconstructed Hit",120, -30, 30));
  diffphiGUNHIT->SetDirectory(0); 
  RLocHRecoH= std::auto_ptr<TH1F>(new  TH1F("RLocHRecoH"," Reconstructed Hit R - GEANT4 Hit R",100,-1.5, 1.5));
  RLocHRecoH->SetDirectory(0); 
  diffRCluHit= std::auto_ptr<TH1F>(new  TH1F("diffRCluHit"," Cluster R - GEANT4 Hit R",100,-1.5, 1.5));
  diffRCluHit->SetDirectory(0); 
  SourceofSecondary= std::auto_ptr<TH1F>(new TH1F("SourceofSecondary","#eta which generate secondaries",1800,-9.0,9.0));
  SourceofSecondary->SetDirectory(0);
  SourceofReco= std::auto_ptr<TH1F>(new TH1F("SourceofReco","#eta which generate reconstruction",1800,-9.0,9.0));
  SourceofReco->SetDirectory(0);
  SymbIdHT0= std::auto_ptr<TH1F>(new TH1F("SymbIdHT0","Symb det id",21,-0.5,19.5));
  SymbIdHT0->SetDirectory(0);
  TrkEtaGoodH0= std::auto_ptr<TH1F>(new TH1F("TrkEtaGoodH0","Reconstructed Track #eta (all cut used)",1800,-9.0,9.0));
  TrkEtaGoodH0->SetDirectory(0); 
  TrkEtaGoodH1= std::auto_ptr<TH1F>(new TH1F("TrkEtaGoodH1","Reconstructed Track #eta (all cut used)",1800,-9.0,9.0));
  TrkEtaGoodH1->SetDirectory(0); 
  TrkEtaH0= std::auto_ptr<TH1F>(new TH1F("TrkEtaH0","Reconstructed Track #eta (no cut used)",1800,-9.0,9.0));
  TrkEtaH0->SetDirectory(0); 
  TrkEtaH1= std::auto_ptr<TH1F>(new TH1F("TrkEtaH1","Reconstructed Track #eta (no cut used)",1800,-9.0,9.0));
  TrkEtaH1->SetDirectory(0); 

  Trketa= std::auto_ptr<TH1F>(new TH1F("Trketa","Reconstructed Track #eta (no cut used)",1800,-9.0,9.0));
  Trketa->SetDirectory(0); 
  Trkphi= std::auto_ptr<TH1F>(new TH1F("Trkphi","Reconstructed Track #phi (no cut used)",361,0,360));
  Trkphi->SetDirectory(0); 
  DPhiGoodTrk= std::auto_ptr<TH1F>(new TH1F("DPhiGoodTrk","#phi Track - #phi generated particle (Only primary tracks are considered)",51,-25.5,25.5));
  DPhiGoodTrk->SetDirectory(0); 
  DEtaGoodTrk= std::auto_ptr<TH1F>(new TH1F("DEtaGoodTrk","#eta Track - #eta generated particle (Only primary tracks are considered)",200,-1.0,1.0));
  DEtaGoodTrk->SetDirectory(0); 
  diffphiCluGun = std::auto_ptr<TH1F>(new TH1F("diffphiCluGun","Cluster #Phi - #Phi generated particle",81, -40.5, 40.5));
  diffphiCluGun->SetDirectory(0); 
  diffphiCluHit= std::auto_ptr<TH1F>(new TH1F("diffphiCluHit","Cluster #Phi - GEANT4 Hit #Phi",50, -5, 5));
  diffphiCluHit->SetDirectory(0); 
  Trketagood= std::auto_ptr<TH1F>(new TH1F("Trketagood","Reconstructed Track #eta (all cut used)",1800,-9.0,9.0));
  Trketagood->SetDirectory(0); 
  Trkphigood= std::auto_ptr<TH1F>(new TH1F("Trkphigood","Reconstructed Track #phi (all cut used)",361,0,360));
  Trkphigood->SetDirectory(0); 

  TrKPhi0degstudyg= std::auto_ptr<TH1F>(new TH1F("TrKPhi0degstudyg","trk phi",720, -270.25, 89.75));
  TrKPhi0degstudyg->SetDirectory(0); 

  DEtaChiCutOneTrk= std::auto_ptr<TH1F>(new TH1F("DEtaChiCutOneTrk","#eta Track - #eta generated particle, chi cut",200,-1.0,1.0));
  DEtaChiCutOneTrk->SetDirectory(0); 
  DPhiChiCutOneTrk= std::auto_ptr<TH1F>(new TH1F("DPhiChiCutOneTrk","#phi Track - #phi generated particle, chi cut ",51,-25.5,25.5));
  DPhiChiCutOneTrk->SetDirectory(0); 

  Chi2RProb = std::auto_ptr<TH1F>(new TH1F("Chi2RProb","Radial #chi^{2} probability",105,0.0,1.05));
  Chi2RProb->SetDirectory(0); 
  Chi2PhiProb = std::auto_ptr<TH1F>(new TH1F("Chi2PhiProb","Azimuthal #chi^{2} probability",105,0.0,1.05));
  Chi2PhiProb->SetDirectory(0);

  chiRdistro = std::auto_ptr<TH1F>(new TH1F("chiRdistro","Radial #chi^{2} ",105,0.0,100));
  chiRdistro->SetDirectory(0);
  chiPhidistro = std::auto_ptr<TH1F>(new TH1F("chiPhidistro","Azimuthal #chi^{2}",105,0.0,100));
  chiPhidistro->SetDirectory(0);


  stablepdg = std::auto_ptr<TH1F>(new TH1F("stablepdg","Pdg Id stable particles",200001,-10000.5,10000.5));
  stablepdg->SetDirectory(0);
  energyP04575= std::auto_ptr<TH1F>(new TH1F("energyP04575","#pi^{0} energy",2000,0.,2000.));
  energyP04575->SetDirectory(0);
  etaP04575= std::auto_ptr<TH1F>(new TH1F("etaP04575","#pi^{0} #eta",60,4.5,7.5));
  etaP04575->SetDirectory(0);
  
  P0EtaEnergy= std::auto_ptr<TH2F>(new TH2F("P0EtaEnergy","#eta vs Energy #pi^{0}",30,4.5,7.5,40,1.,400.));
  P0EtaEnergy->SetDirectory(0);

  ChOutEtaEnergy= std::auto_ptr<TH2F>(new TH2F("ChOutEtaEnergy","#eta vs Energy CH PARTICLE OUTSIDE T2",30,4.5,8.5,40,1.,400.));
  ChOutEtaEnergy->SetDirectory(0);



  ChEnergyinT2= std::auto_ptr<TH1F>(new TH1F("ChEnergyinT2","Charged particle energy in T2",100,0,1000));
  ChEnergyinT2->SetDirectory(0);



  muEtaEnergy= std::auto_ptr<TH2F>(new TH2F("muEtaEnergy","#eta vs Energy muons",30,4.5,8.5,40,1.,400.));
  muEtaEnergy->SetDirectory(0);


  P0EtaEnergycorr= std::auto_ptr<TProfile>(new TProfile("P0EtaEnergycorr","#eta vs Energy #pi^{0}",60,4.5,7.5,""));
  P0EtaEnergycorr->SetDirectory(0);

  P0NumEtacorr= std::auto_ptr<TProfile>(new TProfile("P0NumEtacorr"," # #pi^{0} VS #pi^{0} <#eta>",21,-0.5,20.5,""));
  P0NumEtacorr->SetDirectory(0);

  

  P0NumEnergycorr= std::auto_ptr<TProfile>(new TProfile("P0NumEnergycorr"," # #pi^{0} VS #pi^{0} <energy>",21,-0.5,20.5,""));
  P0NumEnergycorr->SetDirectory(0);

  SingleParticleEfficiencyCutEta= std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCutEta","# Tracks reconstructed vs generated particle #eta (cut #eta)",34,4.35,7.75,""));//Prima l'estr. era a 7.75
  SingleParticleEfficiencyCutEta->SetDirectory(0);
  SingleParticleEfficiencyCutDZ= std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCutDZ","# Tracks reconstructed vs generated particle #eta (cut DZ)",34,4.35,7.75,""));
  SingleParticleEfficiencyCutDZ->SetDirectory(0);
  SingleParticleEfficiencyCutChi= std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyCutChi","# Tracks reconstructed vs generated particle #eta (cut #c)",34,4.35,7.75,""));
  SingleParticleEfficiencyCutChi->SetDirectory(0);

  SingleParticleEfficiencyAllCuts = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyAllCuts","# Tracks reconstructed vs generated particle #eta (all cuts)",34,4.35,7.75,""));
  SingleParticleEfficiencyAllCuts->SetDirectory(0); 
  SingleParticleEfficiency = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiency","# Tracks reconstructed vs generated particle #eta (all track)",34,4.35,8.75,""));
  SingleParticleEfficiency->SetDirectory(0); 
  MultiParticleEfficiencyCut= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCut","<# Tracks> vs # Charged Particles (all track cuts)",12,-0.5,11.5,""));
  MultiParticleEfficiencyCut->SetDirectory(0); 
  MultiParticleEfficiencyCutNorm= std::auto_ptr<TProfile>(new TProfile("MultiParticleEfficiencyCutNorm","<# Tracks / # Charged Particles > vs # Charged Particles (all track cuts)",12,-0.5,11.5,"")); 
  MultiParticleEfficiencyCutNorm->SetDirectory(0); 
  MeanT2HitvsEta = std::auto_ptr<TProfile>(new TProfile("MeanT2HitvsEta","<# Hit> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2HitvsEta->SetDirectory(0);

  MeanT2PadDigivsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2PadDigivsEta","<# detector> with at least 1 Digi-Pad vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2PadDigivsEta->SetDirectory(0);
 MeanT2StripDigivsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2StripDigivsEta","<# detector> with at least 1 Digi-Strip vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2StripDigivsEta->SetDirectory(0);

  MeanT2HitCl1vsEta = std::auto_ptr<TProfile>(new TProfile("MeanT2HitCl1vsEta","<# Hit> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2HitCl1vsEta->SetDirectory(0);

  MeanT2GeantHitvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2GeantHitvsEta","<# detector> with at least 1 Geant Hit vs generated particle #eta ",34,4.35,7.75,""));
  MeanT2GeantHitvsEta->SetDirectory(0);

  // SingleParticleEfficiencyAllCutsH0 = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyAllCutsH0","# Tracks reconstructed vs generated particle #eta (all cuts) HalfTele0",34,4.35,7.75,""));
  //SingleParticleEfficiencyAllCutsH0->SetDirectory(0);

  // SingleParticleEfficiencyAllCutsH1 = std::auto_ptr<TProfile>(new TProfile("SingleParticleEfficiencyAllCutsH1","# Tracks reconstructed vs generated particle #eta (all cuts) HalfTele1",34,4.35,7.75,""));
  //SingleParticleEfficiencyAllCutsH1->SetDirectory(0);

  MeanT2StripvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2StripvsEta","<# Strip> reconstructed vs generated particle #eta ",34,4.35,7.75,"")); 
  MeanT2PadvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2PadvsEta","<#  Pad> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  
  MeanT2RoadEntriesvsEta= std::auto_ptr<TProfile>(new TProfile("MeanT2RoadEntriesvsEta","<# Road Entries> reconstructed vs generated particle #eta ",34,4.35,7.75,""));
  
  MeanEtaResvsPhi= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhi","<Track #eta>  vs #phi (no cuts)",24,0.0,360.0,""));
  MeanEtaResvsPhi->SetDirectory(0); 
  MeanEtaResvsPhiHem1= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem1","<Track #eta>  vs #phi (no cuts) hemis. 1",24,0.0,360.0,""));//15 deg
  MeanEtaResvsPhiHem1->SetDirectory(0); 
  MeanEtaResvsPhiHem2= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem2","<Track #eta>  vs #phi (no cuts) hemis. 2",24,0.0,360.0,""));
  MeanEtaResvsPhiHem2->SetDirectory(0); 
  EtaParticle=std::auto_ptr<TH1F>(new TH1F("EtaParticle","Generated Particle #eta",1800,-9.0,9.0));
  EtaParticle->SetDirectory(0); 
  EtaParticleAll=std::auto_ptr<TH1F>(new TH1F("EtaParticleAll","Generated Particle #eta",60,-15.0,15.0));
  EtaParticleAll->SetDirectory(0);
  MeanEtaResvsPhiPart= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiPart","<Generated #eta>  vs #phi ",24,0.0,360.0,""));
  MeanEtaResvsPhiPart->SetDirectory(0); 
  MeanEtaResvsPhiHem1Part= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem1Part","<Generated #eta>  vs #phi  hemis. 1",24,0.0,360.0,""));//15 deg
  MeanEtaResvsPhiHem1Part->SetDirectory(0); 
  MeanEtaResvsPhiHem2Part= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem2Part","<Generated #eta>  vs #phi  hemis. 2",24,0.0,360.0,""));
  MeanEtaResvsPhiHem2Part->SetDirectory(0); 
  MeanEtaResvsPhiCut= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiCut","<Track #eta>  vs #phi (all cuts)",24,0.0,360.0,""));
  MeanEtaResvsPhiCut->SetDirectory(0);   
  MeanEtaResvsPhiHem1Cut= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem1Cut","<Track #eta>  vs #phi (all cuts) hemis. 1",24,0.0,360.0,""));
  MeanEtaResvsPhiHem1Cut->SetDirectory(0); 
  MeanEtaResvsPhiHem2Cut= std::auto_ptr<TProfile>(new TProfile("MeanEtaResvsPhiHem2Cut","<Track #eta>  vs #phi (all cuts) hemis. 2",24,0.0,360.0,""));
  MeanEtaResvsPhiHem2Cut->SetDirectory(0); 

  MinDz2Plane= std::auto_ptr<TH1F>(new TH1F("MinDz2Plane","Minimum Z Distance resp. to second plane",400,0.0,400.0));
  MinDz2Plane->SetDirectory(0); 

  Trketachicut= std::auto_ptr<TH1F>(new TH1F("Trketachicut","Reconstructed Track #eta (#chi^{2}-p cut)",1800,-9.0,9.0));
 Trketachicut->SetDirectory(0);
  Trketaetacut = std::auto_ptr<TH1F>(new TH1F("Trketaetacut","Reconstructed Track #eta (#eta cut)",1800,-9.0,9.0));
  Trketaetacut->SetDirectory(0);
  
  Trketaetachicut= std::auto_ptr<TH1F>(new TH1F("Trketaetachicut","Reconstructed Track #eta (#eta , #chi^{2}-p cut)",1800,-9.0,9.0));
  Trketaetachicut->SetDirectory(0);
 
  Trketadzcut= std::auto_ptr<TH1F>(new TH1F("Trketadzcut","Reconstructed Track #eta (#Delta Z cut)",1800,-9.0,9.0));
  Trketadzcut->SetDirectory(0);
  Trketadzetacut= std::auto_ptr<TH1F>(new TH1F("Trketadzetacut","Reconstructed Track #eta (#eta, #Delta Z  cut)",1800,-9.0,9.0));
  Trketadzetacut->SetDirectory(0); 
  Trketadzetachicut= std::auto_ptr<TH1F>(new TH1F("Trketadzetachicut","Reconstructed Track #eta (#eta, #Delta Z  cut)",1800,-9.0,9.0));
  Trketadzetachicut->SetDirectory(0);
  

  NumChPartVsNumP0 = std::auto_ptr<TProfile>(new TProfile("NumChPartVsNumP0","# P0 in 4.5-7.5 vs # Ch Part in T2",21,-0.5,20.5,""));
  NumChPartInVsNumChPartOut= std::auto_ptr<TProfile>(new TProfile("NumChPartInVsNumChPartOut","# Ch outside T2  vs # Ch Part in T2",21,-0.5,20.5,""));

  NumP0VsNumTracks= std::auto_ptr<TProfile>(new TProfile("NumP0VsNumTracks","# P0 in 4.5-7.5  vs # Tracks (all cuts)",15,-0.5,14.5,""));

  

  char sZname2[1024];
 char sZnamehist[1024];

  for(unsigned int m=0;m<15; m++)
   {
       sprintf(sZname2, "PimenoEtaEnergy%d", m); 
       sprintf(sZnamehist, "PimenoEtaEnergy%d", m);
       PimenoEtaEnergy[m] = std::auto_ptr<TH2F>(new TH2F(sZname2,sZnamehist,20,5.0,7.0,40,1.,400.));
       PimenoEtaEnergy[m]->SetDirectory(0);
       PimenoEtaEnergy[m]->SetOption("lego");
   }

  chPartinT2->SetXTitle("# ch. Particle in T2");
  ChEnergyinT2->SetXTitle("ch. Particle Energy in T2");
  clusterstripentries->SetXTitle("# strips in one strip-cluster");
  clusterstripentries->SetYTitle("# cluster-strip");
  clusterpadentries->SetXTitle("# pad in one pad-cluster");
  clusterpadentries->SetYTitle("# cluster-pad");
  clusterpadentries->GetYaxis()->SetTitleOffset(1.25);
  clusterstripentries->GetYaxis()->SetTitleOffset(1.25);
  diffphiGUNHIT->SetXTitle("#Delta #Phi (deg)");
  RLocHRecoH->SetXTitle("#Delta R (mm)");
  diffRCluHit->SetXTitle("#Delta R (mm)");
  Trketa->SetXTitle("Track #eta");

  Trketachicut->SetXTitle("Track #eta");
  Trketaetacut->SetXTitle("Track #eta");
  Trketaetachicut->SetXTitle("Track #eta");
  Trketadzcut->SetXTitle("Track #eta");
  Trketadzetacut->SetXTitle("Track #eta");
  Trketadzetachicut->SetXTitle("Track #eta");


  Trkphi->SetXTitle("Track #phi (deg)");
  TrKPhi0degstudyg->SetXTitle("Track #phi (deg)");
  DEtaGoodTrk->SetXTitle("#Delta #eta");
  diffphiCluGun->SetXTitle("#Delta #Phi (deg)");
  diffphiCluHit->SetXTitle("#Delta #Phi (deg)");
  DPhiGoodTrk->SetXTitle("#Delta #Phi (deg)");
  Trkphigood->SetXTitle("Track #phi (deg)");
  Trketagood->SetXTitle("Track #eta");
  SingleParticleEfficiency->SetXTitle("Particle #eta");
  SingleParticleEfficiencyAllCuts->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCutEta->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCutChi->SetXTitle("Particle #eta");
  SingleParticleEfficiencyCutDZ ->SetXTitle("Particle #eta"); 
  SingleParticleEfficiency->SetYTitle("Pseudo-Efficiency");
  SingleParticleEfficiencyAllCuts->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetYTitle("Efficiency");
  MultiParticleEfficiencyCutNorm->SetXTitle("ch. particle in T2");  
  MultiParticleEfficiencyCut->SetYTitle("<Track>");
  MultiParticleEfficiencyCut->SetYTitle("ch. particle in T2");    
  Chi2RProb->SetXTitle("Probability");
  Chi2PhiProb->SetXTitle("Probability");
  HTrackcounter->SetXTitle("# Tracks"); 
  HNumrecotrackcutmag0->SetXTitle("# Tracks Reconstructed");
  HTrackcounterNonCrit->SetXTitle("# Tracks Reconstructed");
  HTrackcounterCrit->SetXTitle("# Tracks Reconstructed");  
  stablepdg->SetXTitle("# Particle ID");  
  NumhitinTrack->SetXTitle("# Hit");  
  NumhitinTrackGood->SetXTitle("# Hit");  

  P0EtaEnergy->SetXTitle("#eta"); 
  P0EtaEnergy->SetYTitle("Energy (GeV)"); 
  // P0EtaEnergy->Draw("");
  P0EtaEnergy->SetOption("lego");
  ChOutEtaEnergy->SetOption("lego");

  muEtaEnergy->SetOption("lego");
  //Chi2RProbLogy->cd();
  // Chi2RProb->Draw("");
  //Chi2RProbLogy->SetLogy();
  //Chi2PhiProbLogy->cd();
  // Chi2PhiProb->Draw("");
  //Chi2PhiProbLogy->SetLogy();


}

// ------------ method called once each job just after ending the event loop  ------------

void T2GeometryAnalyzer::endJob()
{
  



TFile *f = TFile::Open(outputFileName.c_str(), "recreate");
 if( !f || !f->IsWritable() ){
   std::cout << "Output file not opened correctly !!" << std::endl;
 }
  clusterstripentries->Write("");
  clusterpadentries->Write("");
  diffphiGUNHIT->Write("");
  DPhiGoodTrk->Write("");
  MeanT2HitvsEta->Write("");
  DEtaGoodTrk->Write("");
  RLocHRecoH->Write("");
  SymbIdHT0->Write("");


 for(unsigned int m=0;m<15; m++)
   {
       PimenoEtaEnergy[m]->Write(""); 
   } 






  //  SingleParticleEfficiencyAllCutsH1->Write("");
  //SingleParticleEfficiencyAllCutsH0->Write("");
MeanT2StripDigivsEta->Write("");
 MeanT2PadDigivsEta->Write("");
  TrkEtaGoodH0->Write("");
  TrkEtaGoodH1->Write("");
  TrkEtaH1->Write("");
  TrkEtaH0->Write("");
  MeanT2HitCl1vsEta->Write("");
  chPartinT2->Write("");
  ChEnergyinT2->Write("");
  SourceofSecondary->Write("");
  SourceofReco->Write("");
  NumhitinTrack->Write("");
  diffRCluHit->Write("");
  diffphiCluGun->Write("");
  diffphiCluHit->Write("");
MeanT2RoadEntriesvsEta->Write("");
  SingleParticleEfficiencyAllCuts->Write("");
  SingleParticleEfficiencyCutEta->Write("");
  SingleParticleEfficiencyCutChi->Write("");
  SingleParticleEfficiencyCutDZ->Write("");
  NumhitinTrackGood->Write("");
  DPhiChiCutOneTrk->Write("");
  DEtaChiCutOneTrk->Write("");
  Z0distro->Write("");
  R0distro->Write("");

  chiRdistro->Write("");
  chiPhidistro->Write("");

  tantheta->Write("");
  tanthetam0->Write("");
  SingleParticleEfficiency->Write("");
  TrKPhi0degstudyg->Write("");
  HTrackcounter->Write("");
  HNumrecotrackcutmag0->Write("");
  Trketagood->Write("");
  Trkphigood->Write("");
  Trketa->Write("");
  Trkphi->Write("");
  
  Trketachicut->Write("");
  Trketaetacut->Write("");
  Trketaetachicut->Write("");
  Trketadzcut->Write("");
  Trketadzetacut->Write("");
  Trketadzetachicut->Write("");
MeanT2PadvsEta->Write(""); 
MeanT2StripvsEta->Write("");
  HTrackcounterNonCrit->Write("");
  HTrackcounterCrit->Write("");
  SimuZHit->Write("");
  RecoZHit->Write("");
  MinDz2Plane->Write("");
  EtaParticle->Write("");
  MeanEtaResvsPhiHem2->Write(""); 
  MeanEtaResvsPhiHem1->Write(""); 
  MeanEtaResvsPhi->Write("");
  MeanEtaResvsPhiHem1Part->Write("");
  MeanEtaResvsPhiHem2Part->Write("");  
  MeanEtaResvsPhiPart->Write("");
  EtaParticleAll->Write("");
  MeanEtaResvsPhiHem1Cut->Write("");
  MeanEtaResvsPhiHem2Cut->Write("");
  MeanEtaResvsPhiCut->Write("");
  muEtaEnergy->Write("");
  P0EtaEnergy->Write("");
  ChOutEtaEnergy->Write("");
  P0EtaEnergycorr->Write("");
  NumChPartVsNumP0->Write("");
  NumChPartInVsNumChPartOut->Write("");
  NumP0VsNumTracks->Write("");
  P0NumEtacorr->Write("");
  P0NumEnergycorr->Write("");
  MeanT2GeantHitvsEta->Write("");
  energyP04575->Write("");
  etaP04575->Write("");
  MultiParticleEfficiencyCutNorm->Write("");
  MultiParticleEfficiencyCut->Write("");
  Chi2PhiProb->Write("");
  Chi2RProb->Write(""); 
  //Chi2PhiProbLogy->Write("");
  //Chi2RProbLogy->Write("");
  stablepdg->Write("");
  f->Close();


  
  //file for ntuple storing
  
   
  
  



}

//define this as a plug-in
DEFINE_FWK_MODULE(T2GeometryAnalyzer);
