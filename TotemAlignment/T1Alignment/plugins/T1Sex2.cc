/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
*/
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "TotemAlignment/T1Alignment/interface/T1Sex2.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include <vector>
#include "TFile.h"
#include "TH1.h"
#include "TVector3.h"
#include "TVector2.h"
#include <fstream>
#include <sstream>


#ifndef PI
#define PI 3.14159
#endif
#define _30DEG 0.5236
#define _25DEG 0.4363
#define _60DEG 1.0472
#define _90DEG 1.5708
#define _120DEG 2.0944
#define _150DEG 2.618
#define _210DEG 3.665
#define _240DEG 4.1888
#define _270DEG 4.7124
#define _300DEG 5.2360
#define _330DEG 5.7596
#define _3DEG 0.05239878

using namespace edm;
using namespace std;

//#define _DEBUG_

T1Sex2::T1Sex2(const edm::ParameterSet& iConfig):_Verbosity(0)
{
 
  _Verbosity =  iConfig.getParameter<int>("Verbosity");

}


T1Sex2::~T1Sex2()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
T1Sex2::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{

 
/*  float Z_temp[5][6]={{7521.4,7569.6,7526.4,7569.6,7521.4,7564.6},
		      {8191.4,8239.4,8196.4,8239.6,8191.4,8234.6},
		      {8798.4,8846.6,8803.4,8846.6,8798.4,8841.6},
		      {9405.4,9453.6,9410.4,9453.6,9405.4,9448.6},
		      {10168.4,10216.6,10173.4,10216.6,10168.4,10211.6}};
*/





  edm::Handle<T1T2TrackCollection> trackCollection;
  iEvent.getByLabel("t1tracks2","T1TrackColl",trackCollection);


  int TracceRicostruite=0;


#ifdef _DEBUG_
  std::cout << " Taglia della track collection: " << trackCollection->size() <<std::endl;
#endif
  if(_Verbosity>0)
    if (trackCollection->size() == 0){ std::cout << " WARNING: No tracks in T1 in the event !      ]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]]] " <<std::endl;
    }
    else{
      std::cout << "                                                                " <<trackCollection->size() << " tracks in the event " << std::endl;
    }
  TracceRicostruite=trackCollection->size();


  
  T1T2TrackCollection::const_iterator TC_it;




  if(TracceRicostruite<10)
    for(TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++){
      unsigned int entries = (*TC_it).GetHitEntries();
      float PhiDEG = ((*TC_it).Phi()/3.14159*180);
      bool traccia_al_confine = false;

      if(  fabs((*TC_it).Z_at_Rmin()) < 1500 && (*TC_it).Rmin()<200 && (*TC_it).ChiSquared()<60 
	   && (*TC_it).Eta()>0 && entries>=4
	   ){

	if(PhiDEG > 24 && PhiDEG < 36){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();

	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _plus.push_back(V3);
	      _plus_prev.push_back(V2);
	      _plus_error.push_back(eV2);


	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_0_0_50->Fill(DDDx);
		hDDDy_0_0_50->Fill(DDDy);
		break;
	      case 1:
		hDDDx_0_1_50->Fill(DDDx);
		hDDDy_0_1_50->Fill(DDDy);
		break;
	      case 2:
		hDDDx_0_2_50->Fill(DDDx);
		hDDDy_0_2_50->Fill(DDDy);
		break;
	      case 3:
		hDDDx_0_3_50->Fill(DDDx);
		hDDDy_0_3_50->Fill(DDDy);
		break;
	      case 4:
		hDDDx_0_4_50->Fill(DDDx);
		hDDDy_0_4_50->Fill(DDDy);
		break;


	      }

	    }

	  }


	}





      





	if(PhiDEG > 84 && PhiDEG < 96){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();

	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _plus.push_back(V3);
	      _plus_prev.push_back(V2);
	      _plus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_0_0_01->Fill(DDDx);
		hDDDy_0_0_01->Fill(DDDy);
		break;
	      case 1:
		hDDDx_0_1_01->Fill(DDDx);
		hDDDy_0_1_01->Fill(DDDy);
		break;
	      case 2:
		hDDDx_0_2_01->Fill(DDDx);
		hDDDy_0_2_01->Fill(DDDy);
		break;
	      case 3:
		hDDDx_0_3_01->Fill(DDDx);
		hDDDy_0_3_01->Fill(DDDy);
		break;
	      case 4:
		hDDDx_0_4_01->Fill(DDDx);
		hDDDy_0_4_01->Fill(DDDy);
		break;


	      }

	    }

	  }


	}















	if(PhiDEG > 144 && PhiDEG < 156){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();

	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _plus.push_back(V3);
	      _plus_prev.push_back(V2);
	      _plus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_0_0_12->Fill(DDDx);
		hDDDy_0_0_12->Fill(DDDy);
		break;
	      case 1:
		hDDDx_0_1_12->Fill(DDDx);
		hDDDy_0_1_12->Fill(DDDy);
		break;
	      case 2:
		hDDDx_0_2_12->Fill(DDDx);
		hDDDy_0_2_12->Fill(DDDy);
		break;
	      case 3:
		hDDDx_0_3_12->Fill(DDDx);
		hDDDy_0_3_12->Fill(DDDy);
		break;
	      case 4:
		hDDDx_0_4_12->Fill(DDDx);
		hDDDy_0_4_12->Fill(DDDy);
		break;


	      }

	    }

	  }


	}








	if(PhiDEG > 204 && PhiDEG < 216){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();

	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _plus.push_back(V3);
	      _plus_prev.push_back(V2);
	      _plus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_0_0_23->Fill(DDDx);
		hDDDy_0_0_23->Fill(DDDy);
		break;
	      case 1:
		hDDDx_0_1_23->Fill(DDDx);
		hDDDy_0_1_23->Fill(DDDy);
		break;
	      case 2:
		hDDDx_0_2_23->Fill(DDDx);
		hDDDy_0_2_23->Fill(DDDy);
		break;
	      case 3:
		hDDDx_0_3_23->Fill(DDDx);
		hDDDy_0_3_23->Fill(DDDy);
		break;
	      case 4:
		hDDDx_0_4_23->Fill(DDDx);
		hDDDy_0_4_23->Fill(DDDy);
		break;


	      }

	    }

	  }


	}








	if(PhiDEG > 264 && PhiDEG < 276){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();

	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _plus.push_back(V3);
	      _plus_prev.push_back(V2);
	      _plus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_0_0_34->Fill(DDDx);
		hDDDy_0_0_34->Fill(DDDy);
		break;
	      case 1:
		hDDDx_0_1_34->Fill(DDDx);
		hDDDy_0_1_34->Fill(DDDy);
		break;
	      case 2:
		hDDDx_0_2_34->Fill(DDDx);
		hDDDy_0_2_34->Fill(DDDy);
		break;
	      case 3:
		hDDDx_0_3_34->Fill(DDDx);
		hDDDy_0_3_34->Fill(DDDy);
		break;
	      case 4:
		hDDDx_0_4_34->Fill(DDDx);
		hDDDy_0_4_34->Fill(DDDy);
		break;


	      }

	    }

	  }


	}












	if(PhiDEG > 324 && PhiDEG < 336){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();

	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _plus.push_back(V3);
	      _plus_prev.push_back(V2);
	      _plus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_0_0_45->Fill(DDDx);
		hDDDy_0_0_45->Fill(DDDy);
		break;
	      case 1:
		hDDDx_0_1_45->Fill(DDDx);
		hDDDy_0_1_45->Fill(DDDy);
		break;
	      case 2:
		hDDDx_0_2_45->Fill(DDDx);
		hDDDy_0_2_45->Fill(DDDy);
		break;
	      case 3:
		hDDDx_0_3_45->Fill(DDDx);
		hDDDy_0_3_45->Fill(DDDy);
		break;
	      case 4:
		hDDDx_0_4_45->Fill(DDDx);
		hDDDy_0_4_45->Fill(DDDy);
		break;


	      }

	    }

	  }


	}












      }
      
     
      if(  fabs((*TC_it).Z_at_Rmin()) < 1500 && (*TC_it).Rmin()<200 && (*TC_it).ChiSquared()<60 
	   && (*TC_it).Eta()<0 && entries>=4
	   ){

	if(PhiDEG > 24 && PhiDEG < 36){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();

	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _minus.push_back(V3);
	      _minus_prev.push_back(V2);
	      _minus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_1_0_50->Fill(DDDx);
		hDDDy_1_0_50->Fill(DDDy);
		break;
	      case 1:
		hDDDx_1_1_50->Fill(DDDx);
		hDDDy_1_1_50->Fill(DDDy);
		break;
	      case 2:
		hDDDx_1_2_50->Fill(DDDx);
		hDDDy_1_2_50->Fill(DDDy);
		break;
	      case 3:
		hDDDx_1_3_50->Fill(DDDx);
		hDDDy_1_3_50->Fill(DDDy);
		break;
	      case 4:
		hDDDx_1_4_50->Fill(DDDx);
		hDDDy_1_4_50->Fill(DDDy);
		break;


	      }

	    }

	  }


	}





      





	if(PhiDEG > 84 && PhiDEG < 96){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();
	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _minus.push_back(V3);
	      _minus_prev.push_back(V2);
	      _minus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_1_0_01->Fill(DDDx);
		hDDDy_1_0_01->Fill(DDDy);
		break;
	      case 1:
		hDDDx_1_1_01->Fill(DDDx);
		hDDDy_1_1_01->Fill(DDDy);
		break;
	      case 2:
		hDDDx_1_2_01->Fill(DDDx);
		hDDDy_1_2_01->Fill(DDDy);
		break;
	      case 3:
		hDDDx_1_3_01->Fill(DDDx);
		hDDDy_1_3_01->Fill(DDDy);
		break;
	      case 4:
		hDDDx_1_4_01->Fill(DDDx);
		hDDDy_1_4_01->Fill(DDDy);
		break;


	      }

	    }

	  }


	}















	if(PhiDEG > 144 && PhiDEG < 156){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();
	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _minus.push_back(V3);
	      _minus_prev.push_back(V2);
	      _minus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_1_0_12->Fill(DDDx);
		hDDDy_1_0_12->Fill(DDDy);
		break;
	      case 1:
		hDDDx_1_1_12->Fill(DDDx);
		hDDDy_1_1_12->Fill(DDDy);
		break;
	      case 2:
		hDDDx_1_2_12->Fill(DDDx);
		hDDDy_1_2_12->Fill(DDDy);
		break;
	      case 3:
		hDDDx_1_3_12->Fill(DDDx);
		hDDDy_1_3_12->Fill(DDDy);
		break;
	      case 4:
		hDDDx_1_4_12->Fill(DDDx);
		hDDDy_1_4_12->Fill(DDDy);
		break;


	      }

	    }

	  }


	}








	if(PhiDEG > 204 && PhiDEG < 216){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();
	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _minus.push_back(V3);
	      _minus_prev.push_back(V2);
	      _minus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_1_0_23->Fill(DDDx);
		hDDDy_1_0_23->Fill(DDDy);
		break;
	      case 1:
		hDDDx_1_1_23->Fill(DDDx);
		hDDDy_1_1_23->Fill(DDDy);
		break;
	      case 2:
		hDDDx_1_2_23->Fill(DDDx);
		hDDDy_1_2_23->Fill(DDDy);
		break;
	      case 3:
		hDDDx_1_3_23->Fill(DDDx);
		hDDDy_1_3_23->Fill(DDDy);
		break;
	      case 4:
		hDDDx_1_4_23->Fill(DDDx);
		hDDDy_1_4_23->Fill(DDDy);
		break;


	      }

	    }

	  }


	}








	if(PhiDEG > 264 && PhiDEG < 276){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();
	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _minus.push_back(V3);
	      _minus_prev.push_back(V2);
	      _minus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_1_0_34->Fill(DDDx);
		hDDDy_1_0_34->Fill(DDDy);
		break;
	      case 1:
		hDDDx_1_1_34->Fill(DDDx);
		hDDDy_1_1_34->Fill(DDDy);
		break;
	      case 2:
		hDDDx_1_2_34->Fill(DDDx);
		hDDDy_1_2_34->Fill(DDDy);
		break;
	      case 3:
		hDDDx_1_3_34->Fill(DDDx);
		hDDDy_1_3_34->Fill(DDDy);
		break;
	      case 4:
		hDDDx_1_4_34->Fill(DDDx);
		hDDDy_1_4_34->Fill(DDDy);
		break;


	      }

	    }

	  }


	}












	if(PhiDEG > 324 && PhiDEG < 336){


	  float x00 = (*TC_it).GetHitT1(0).GlobalPosition().x();
	  float y00 = (*TC_it).GetHitT1(0).GlobalPosition().y();
	  float z00 = (*TC_it).GetHitT1(0).GlobalPosition().z();
	  int sestante = sextant(x00,y00,z00);
	  for(unsigned int i =1 ; i< entries; i++){
	    float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	    float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	    float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	    int sestante_i = sextant(xi,yi,zi);
	    if(sestante != sestante_i)traccia_al_confine = true;
	  }
	  if(traccia_al_confine){
//	  cout << "traccia al confine"<<endl;
	    for(unsigned int i =0 ; i< entries; i++){
	      float xi = (*TC_it).GetHitT1(i).GlobalPosition().x();
	      float yi = (*TC_it).GetHitT1(i).GlobalPosition().y();
	      float zi = (*TC_it).GetHitT1(i).GlobalPosition().z();
	      float x_prev=(*TC_it).GetTx() * zi + (*TC_it).X0();
	      float y_prev=(*TC_it).GetTy() * zi + (*TC_it).Y0();
	      float eeX=(*TC_it).GetTxSigma()*fabs(zi)+(*TC_it).X0Sigma();
	      float eeY=(*TC_it).GetTySigma()*fabs(zi)+(*TC_it).Y0Sigma();
	      TVector3 V3(xi,yi,zi);
	      TVector2 V2(x_prev,y_prev);
	      TVector2 eV2(eeX,eeY);
	      
	      _minus.push_back(V3);
	      _minus_prev.push_back(V2);
	      _minus_error.push_back(eV2);

	      float DDDx = xi - x_prev;
	      float DDDy = yi - y_prev;
	      int piano = plane(zi);
	  
	      switch(piano){
	      case 0:
		hDDDx_1_0_45->Fill(DDDx);
		hDDDy_1_0_45->Fill(DDDy);
		break;
	      case 1:
		hDDDx_1_1_45->Fill(DDDx);
		hDDDy_1_1_45->Fill(DDDy);
		break;
	      case 2:
		hDDDx_1_2_45->Fill(DDDx);
		hDDDy_1_2_45->Fill(DDDy);
		break;
	      case 3:
		hDDDx_1_3_45->Fill(DDDx);
		hDDDy_1_3_45->Fill(DDDy);
		break;
	      case 4:
		hDDDx_1_4_45->Fill(DDDx);
		hDDDy_1_4_45->Fill(DDDy);
		break;


	      }

	    }

	  }


	}












      }


 
    }
}






// ------------ method called once each job just before starting event loop  ------------
void 
T1Sex2::beginJob()
{
  if(_Verbosity>0)
    cout << "T1Sex2::beginJob" << endl;


  assert(  _plus.empty() && 
	   _plus_prev.empty() && 
	   _plus_error.empty() &&
	   
	   _minus.empty() &&
	   _minus_prev.empty() &&
	   _minus_error.empty() );






  hDDDx_0_0_01 = new TH1D("DeltaX_0_0_01","DeltaX_0_0_01",100,-50.,50.);
  hDDDy_0_0_01 = new TH1D("DeltaY_0_0_01","DeltaY_0_0_01",100,-50.,50.);

  hDDDx_0_0_12 = new TH1D("DeltaX_0_0_12","DeltaX_0_0_12",100,-50.,50.);
  hDDDy_0_0_12 = new TH1D("DeltaY_0_0_12","DeltaY_0_0_12",100,-50.,50.);

  hDDDx_0_0_23 = new TH1D("DeltaX_0_0_23","DeltaX_0_0_23",100,-50.,50.);
  hDDDy_0_0_23 = new TH1D("DeltaY_0_0_23","DeltaY_0_0_23",100,-50.,50.);

  hDDDx_0_0_34 = new TH1D("DeltaX_0_0_34","DeltaX_0_0_34",100,-50.,50.);
  hDDDy_0_0_34 = new TH1D("DeltaY_0_0_34","DeltaY_0_0_34",100,-50.,50.);

  hDDDx_0_0_45 = new TH1D("DeltaX_0_0_45","DeltaX_0_0_45",100,-50.,50.);
  hDDDy_0_0_45 = new TH1D("DeltaY_0_0_45","DeltaY_0_0_45",100,-50.,50.);

  hDDDx_0_0_50 = new TH1D("DeltaX_0_0_50","DeltaX_0_0_50",100,-50.,50.);
  hDDDy_0_0_50 = new TH1D("DeltaY_0_0_50","DeltaY_0_0_50",100,-50.,50.);





  hDDDx_0_1_01 = new TH1D("DeltaX_0_1_01","DeltaX_0_1_01",100,-50.,50.);
  hDDDy_0_1_01 = new TH1D("DeltaY_0_1_01","DeltaY_0_1_01",100,-50.,50.);

  hDDDx_0_1_12 = new TH1D("DeltaX_0_1_12","DeltaX_0_1_12",100,-50.,50.);
  hDDDy_0_1_12 = new TH1D("DeltaY_0_1_12","DeltaY_0_1_12",100,-50.,50.);

  hDDDx_0_1_23 = new TH1D("DeltaX_0_1_23","DeltaX_0_1_23",100,-50.,50.);
  hDDDy_0_1_23 = new TH1D("DeltaY_0_1_23","DeltaY_0_1_23",100,-50.,50.);

  hDDDx_0_1_34 = new TH1D("DeltaX_0_1_34","DeltaX_0_1_34",100,-50.,50.);
  hDDDy_0_1_34 = new TH1D("DeltaY_0_1_34","DeltaY_0_1_34",100,-50.,50.);

  hDDDx_0_1_45 = new TH1D("DeltaX_0_1_45","DeltaX_0_1_45",100,-50.,50.);
  hDDDy_0_1_45 = new TH1D("DeltaY_0_1_45","DeltaY_0_1_45",100,-50.,50.);

  hDDDx_0_1_50 = new TH1D("DeltaX_0_1_50","DeltaX_0_1_50",100,-50.,50.);
  hDDDy_0_1_50 = new TH1D("DeltaY_0_1_50","DeltaY_0_1_50",100,-50.,50.);



  hDDDx_0_2_01 = new TH1D("DeltaX_0_2_01","DeltaX_0_2_01",100,-50.,50.);
  hDDDy_0_2_01 = new TH1D("DeltaY_0_2_01","DeltaY_0_2_01",100,-50.,50.);

  hDDDx_0_2_12 = new TH1D("DeltaX_0_2_12","DeltaX_0_2_12",100,-50.,50.);
  hDDDy_0_2_12 = new TH1D("DeltaY_0_2_12","DeltaY_0_2_12",100,-50.,50.);

  hDDDx_0_2_23 = new TH1D("DeltaX_0_2_23","DeltaX_0_2_23",100,-50.,50.);
  hDDDy_0_2_23 = new TH1D("DeltaY_0_2_23","DeltaY_0_2_23",100,-50.,50.);

  hDDDx_0_2_34 = new TH1D("DeltaX_0_2_34","DeltaX_0_2_34",100,-50.,50.);
  hDDDy_0_2_34 = new TH1D("DeltaY_0_2_34","DeltaY_0_2_34",100,-50.,50.);

  hDDDx_0_2_45 = new TH1D("DeltaX_0_2_45","DeltaX_0_2_45",100,-50.,50.);
  hDDDy_0_2_45 = new TH1D("DeltaY_0_2_45","DeltaY_0_2_45",100,-50.,50.);

  hDDDx_0_2_50 = new TH1D("DeltaX_0_2_50","DeltaX_0_2_50",100,-50.,50.);
  hDDDy_0_2_50 = new TH1D("DeltaY_0_2_50","DeltaY_0_2_50",100,-50.,50.);





  hDDDx_0_3_01 = new TH1D("DeltaX_0_3_01","DeltaX_0_3_01",100,-50.,50.);
  hDDDy_0_3_01 = new TH1D("DeltaY_0_3_01","DeltaY_0_3_01",100,-50.,50.);

  hDDDx_0_3_12 = new TH1D("DeltaX_0_3_12","DeltaX_0_3_12",100,-50.,50.);
  hDDDy_0_3_12 = new TH1D("DeltaY_0_3_12","DeltaY_0_3_12",100,-50.,50.);

  hDDDx_0_3_23 = new TH1D("DeltaX_0_3_23","DeltaX_0_3_23",100,-50.,50.);
  hDDDy_0_3_23 = new TH1D("DeltaY_0_3_23","DeltaY_0_3_23",100,-50.,50.);

  hDDDx_0_3_34 = new TH1D("DeltaX_0_3_34","DeltaX_0_3_34",100,-50.,50.);
  hDDDy_0_3_34 = new TH1D("DeltaY_0_3_34","DeltaY_0_3_34",100,-50.,50.);

  hDDDx_0_3_45 = new TH1D("DeltaX_0_3_45","DeltaX_0_3_45",100,-50.,50.);
  hDDDy_0_3_45 = new TH1D("DeltaY_0_3_45","DeltaY_0_3_45",100,-50.,50.);

  hDDDx_0_3_50 = new TH1D("DeltaX_0_3_50","DeltaX_0_3_50",100,-50.,50.);
  hDDDy_0_3_50 = new TH1D("DeltaY_0_3_50","DeltaY_0_3_50",100,-50.,50.);





  hDDDx_0_4_01 = new TH1D("DeltaX_0_4_01","DeltaX_0_4_01",100,-50.,50.);
  hDDDy_0_4_01 = new TH1D("DeltaY_0_4_01","DeltaY_0_4_01",100,-50.,50.);

  hDDDx_0_4_12 = new TH1D("DeltaX_0_4_12","DeltaX_0_4_12",100,-50.,50.);
  hDDDy_0_4_12 = new TH1D("DeltaY_0_4_12","DeltaY_0_4_12",100,-50.,50.);

  hDDDx_0_4_23 = new TH1D("DeltaX_0_4_23","DeltaX_0_4_23",100,-50.,50.);
  hDDDy_0_4_23 = new TH1D("DeltaY_0_4_23","DeltaY_0_4_23",100,-50.,50.);

  hDDDx_0_4_34 = new TH1D("DeltaX_0_4_34","DeltaX_0_4_34",100,-50.,50.);
  hDDDy_0_4_34 = new TH1D("DeltaY_0_4_34","DeltaY_0_4_34",100,-50.,50.);

  hDDDx_0_4_45 = new TH1D("DeltaX_0_4_45","DeltaX_0_4_45",100,-50.,50.);
  hDDDy_0_4_45 = new TH1D("DeltaY_0_4_45","DeltaY_0_4_45",100,-50.,50.);

  hDDDx_0_4_50 = new TH1D("DeltaX_0_4_50","DeltaX_0_4_50",100,-50.,50.);
  hDDDy_0_4_50 = new TH1D("DeltaY_0_4_50","DeltaY_0_4_50",100,-50.,50.);









  hDDDx_1_0_01 = new TH1D("Deltax_1_0_01","Deltax_1_0_01",100,-50.,50.);
  hDDDy_1_0_01 = new TH1D("Deltay_1_0_01","Deltay_1_0_01",100,-50.,50.);

  hDDDx_1_0_12 = new TH1D("Deltax_1_0_12","Deltax_1_0_12",100,-50.,50.);
  hDDDy_1_0_12 = new TH1D("Deltay_1_0_12","Deltay_1_0_12",100,-50.,50.);

  hDDDx_1_0_23 = new TH1D("Deltax_1_0_23","Deltax_1_0_23",100,-50.,50.);
  hDDDy_1_0_23 = new TH1D("Deltay_1_0_23","Deltay_1_0_23",100,-50.,50.);

  hDDDx_1_0_34 = new TH1D("Deltax_1_0_34","Deltax_1_0_34",100,-50.,50.);
  hDDDy_1_0_34 = new TH1D("Deltay_1_0_34","Deltay_1_0_34",100,-50.,50.);

  hDDDx_1_0_45 = new TH1D("Deltax_1_0_45","Deltax_1_0_45",100,-50.,50.);
  hDDDy_1_0_45 = new TH1D("Deltay_1_0_45","Deltay_1_0_45",100,-50.,50.);

  hDDDx_1_0_50 = new TH1D("Deltax_1_0_50","Deltax_1_0_50",100,-50.,50.);
  hDDDy_1_0_50 = new TH1D("Deltay_1_0_50","Deltay_1_0_50",100,-50.,50.);





  hDDDx_1_1_01 = new TH1D("Deltax_1_1_01","Deltax_1_1_01",100,-50.,50.);
  hDDDy_1_1_01 = new TH1D("Deltay_1_1_01","Deltay_1_1_01",100,-50.,50.);

  hDDDx_1_1_12 = new TH1D("Deltax_1_1_12","Deltax_1_1_12",100,-50.,50.);
  hDDDy_1_1_12 = new TH1D("Deltay_1_1_12","Deltay_1_1_12",100,-50.,50.);

  hDDDx_1_1_23 = new TH1D("Deltax_1_1_23","Deltax_1_1_23",100,-50.,50.);
  hDDDy_1_1_23 = new TH1D("Deltay_1_1_23","Deltay_1_1_23",100,-50.,50.);

  hDDDx_1_1_34 = new TH1D("Deltax_1_1_34","Deltax_1_1_34",100,-50.,50.);
  hDDDy_1_1_34 = new TH1D("Deltay_1_1_34","Deltay_1_1_34",100,-50.,50.);

  hDDDx_1_1_45 = new TH1D("Deltax_1_1_45","Deltax_1_1_45",100,-50.,50.);
  hDDDy_1_1_45 = new TH1D("Deltay_1_1_45","Deltay_1_1_45",100,-50.,50.);

  hDDDx_1_1_50 = new TH1D("Deltax_1_1_50","Deltax_1_1_50",100,-50.,50.);
  hDDDy_1_1_50 = new TH1D("Deltay_1_1_50","Deltay_1_1_50",100,-50.,50.);



  hDDDx_1_2_01 = new TH1D("Deltax_1_2_01","Deltax_1_2_01",100,-50.,50.);
  hDDDy_1_2_01 = new TH1D("Deltay_1_2_01","Deltay_1_2_01",100,-50.,50.);

  hDDDx_1_2_12 = new TH1D("Deltax_1_2_12","Deltax_1_2_12",100,-50.,50.);
  hDDDy_1_2_12 = new TH1D("Deltay_1_2_12","Deltay_1_2_12",100,-50.,50.);

  hDDDx_1_2_23 = new TH1D("Deltax_1_2_23","Deltax_1_2_23",100,-50.,50.);
  hDDDy_1_2_23 = new TH1D("Deltay_1_2_23","Deltay_1_2_23",100,-50.,50.);

  hDDDx_1_2_34 = new TH1D("Deltax_1_2_34","Deltax_1_2_34",100,-50.,50.);
  hDDDy_1_2_34 = new TH1D("Deltay_1_2_34","Deltay_1_2_34",100,-50.,50.);

  hDDDx_1_2_45 = new TH1D("Deltax_1_2_45","Deltax_1_2_45",100,-50.,50.);
  hDDDy_1_2_45 = new TH1D("Deltay_1_2_45","Deltay_1_2_45",100,-50.,50.);

  hDDDx_1_2_50 = new TH1D("Deltax_1_2_50","Deltax_1_2_50",100,-50.,50.);
  hDDDy_1_2_50 = new TH1D("Deltay_1_2_50","Deltay_1_2_50",100,-50.,50.);





  hDDDx_1_3_01 = new TH1D("Deltax_1_3_01","Deltax_1_3_01",100,-50.,50.);
  hDDDy_1_3_01 = new TH1D("Deltay_1_3_01","Deltay_1_3_01",100,-50.,50.);

  hDDDx_1_3_12 = new TH1D("Deltax_1_3_12","Deltax_1_3_12",100,-50.,50.);
  hDDDy_1_3_12 = new TH1D("Deltay_1_3_12","Deltay_1_3_12",100,-50.,50.);

  hDDDx_1_3_23 = new TH1D("Deltax_1_3_23","Deltax_1_3_23",100,-50.,50.);
  hDDDy_1_3_23 = new TH1D("Deltay_1_3_23","Deltay_1_3_23",100,-50.,50.);

  hDDDx_1_3_34 = new TH1D("Deltax_1_3_34","Deltax_1_3_34",100,-50.,50.);
  hDDDy_1_3_34 = new TH1D("Deltay_1_3_34","Deltay_1_3_34",100,-50.,50.);

  hDDDx_1_3_45 = new TH1D("Deltax_1_3_45","Deltax_1_3_45",100,-50.,50.);
  hDDDy_1_3_45 = new TH1D("Deltay_1_3_45","Deltay_1_3_45",100,-50.,50.);

  hDDDx_1_3_50 = new TH1D("Deltax_1_3_50","Deltax_1_3_50",100,-50.,50.);
  hDDDy_1_3_50 = new TH1D("Deltay_1_3_50","Deltay_1_3_50",100,-50.,50.);





  hDDDx_1_4_01 = new TH1D("Deltax_1_4_01","Deltax_1_4_01",100,-50.,50.);
  hDDDy_1_4_01 = new TH1D("Deltay_1_4_01","Deltay_1_4_01",100,-50.,50.);

  hDDDx_1_4_12 = new TH1D("Deltax_1_4_12","Deltax_1_4_12",100,-50.,50.);
  hDDDy_1_4_12 = new TH1D("Deltay_1_4_12","Deltay_1_4_12",100,-50.,50.);

  hDDDx_1_4_23 = new TH1D("Deltax_1_4_23","Deltax_1_4_23",100,-50.,50.);
  hDDDy_1_4_23 = new TH1D("Deltay_1_4_23","Deltay_1_4_23",100,-50.,50.);

  hDDDx_1_4_34 = new TH1D("Deltax_1_4_34","Deltax_1_4_34",100,-50.,50.);
  hDDDy_1_4_34 = new TH1D("Deltay_1_4_34","Deltay_1_4_34",100,-50.,50.);

  hDDDx_1_4_45 = new TH1D("Deltax_1_4_45","Deltax_1_4_45",100,-50.,50.);
  hDDDy_1_4_45 = new TH1D("Deltay_1_4_45","Deltay_1_4_45",100,-50.,50.);

  hDDDx_1_4_50 = new TH1D("Deltax_1_4_50","Deltax_1_4_50",100,-50.,50.);
  hDDDy_1_4_50 = new TH1D("Deltay_1_4_50","Deltay_1_4_50",100,-50.,50.);






}

// ------------ method called once each job just after ending the event loop  ------------
void T1Sex2::endJob() {
  if(_Verbosity>0)
    cout << "T1Sex2::endJob" << endl;


  theFile = new TFile("sexFile2.root","RECREATE");





  hDDDx_0_0_01 ->Write(); // TH1D("DeltaX_0_0_01","DeltaX_0_0_01",100,-50.,50.);
  hDDDy_0_0_01 ->Write(); // TH1D("DeltaY_0_0_01","DeltaY_0_0_01",100,-50.,50.);

  hDDDx_0_0_12 ->Write(); // TH1D("DeltaX_0_0_12","DeltaX_0_0_12",100,-50.,50.);
  hDDDy_0_0_12 ->Write(); // TH1D("DeltaY_0_0_12","DeltaY_0_0_12",100,-50.,50.);

  hDDDx_0_0_23 ->Write(); // TH1D("DeltaX_0_0_23","DeltaX_0_0_23",100,-50.,50.);
  hDDDy_0_0_23 ->Write(); // TH1D("DeltaY_0_0_23","DeltaY_0_0_23",100,-50.,50.);

  hDDDx_0_0_34 ->Write(); // TH1D("DeltaX_0_0_34","DeltaX_0_0_34",100,-50.,50.);
  hDDDy_0_0_34 ->Write(); // TH1D("DeltaY_0_0_34","DeltaY_0_0_34",100,-50.,50.);

  hDDDx_0_0_45 ->Write(); // TH1D("DeltaX_0_0_45","DeltaX_0_0_45",100,-50.,50.);
  hDDDy_0_0_45 ->Write(); // TH1D("DeltaY_0_0_45","DeltaY_0_0_45",100,-50.,50.);

  hDDDx_0_0_50 ->Write(); // TH1D("DeltaX_0_0_50","DeltaX_0_0_50",100,-50.,50.);
  hDDDy_0_0_50 ->Write(); // TH1D("DeltaY_0_0_50","DeltaY_0_0_50",100,-50.,50.);





  hDDDx_0_1_01 ->Write(); // TH1D("DeltaX_0_1_01","DeltaX_0_1_01",100,-50.,50.);
  hDDDy_0_1_01 ->Write(); // TH1D("DeltaY_0_1_01","DeltaY_0_1_01",100,-50.,50.);

  hDDDx_0_1_12 ->Write(); // TH1D("DeltaX_0_1_12","DeltaX_0_1_12",100,-50.,50.);
  hDDDy_0_1_12 ->Write(); // TH1D("DeltaY_0_1_12","DeltaY_0_1_12",100,-50.,50.);

  hDDDx_0_1_23 ->Write(); // TH1D("DeltaX_0_1_23","DeltaX_0_1_23",100,-50.,50.);
  hDDDy_0_1_23 ->Write(); // TH1D("DeltaY_0_1_23","DeltaY_0_1_23",100,-50.,50.);

  hDDDx_0_1_34 ->Write(); // TH1D("DeltaX_0_1_34","DeltaX_0_1_34",100,-50.,50.);
  hDDDy_0_1_34 ->Write(); // TH1D("DeltaY_0_1_34","DeltaY_0_1_34",100,-50.,50.);

  hDDDx_0_1_45 ->Write(); // TH1D("DeltaX_0_1_45","DeltaX_0_1_45",100,-50.,50.);
  hDDDy_0_1_45 ->Write(); // TH1D("DeltaY_0_1_45","DeltaY_0_1_45",100,-50.,50.);

  hDDDx_0_1_50 ->Write(); // TH1D("DeltaX_0_1_50","DeltaX_0_1_50",100,-50.,50.);
  hDDDy_0_1_50 ->Write(); // TH1D("DeltaY_0_1_50","DeltaY_0_1_50",100,-50.,50.);



  hDDDx_0_2_01 ->Write(); // TH1D("DeltaX_0_2_01","DeltaX_0_2_01",100,-50.,50.);
  hDDDy_0_2_01 ->Write(); // TH1D("DeltaY_0_2_01","DeltaY_0_2_01",100,-50.,50.);

  hDDDx_0_2_12 ->Write(); // TH1D("DeltaX_0_2_12","DeltaX_0_2_12",100,-50.,50.);
  hDDDy_0_2_12 ->Write(); // TH1D("DeltaY_0_2_12","DeltaY_0_2_12",100,-50.,50.);

  hDDDx_0_2_23 ->Write(); // TH1D("DeltaX_0_2_23","DeltaX_0_2_23",100,-50.,50.);
  hDDDy_0_2_23 ->Write(); // TH1D("DeltaY_0_2_23","DeltaY_0_2_23",100,-50.,50.);

  hDDDx_0_2_34 ->Write(); // TH1D("DeltaX_0_2_34","DeltaX_0_2_34",100,-50.,50.);
  hDDDy_0_2_34 ->Write(); // TH1D("DeltaY_0_2_34","DeltaY_0_2_34",100,-50.,50.);

  hDDDx_0_2_45 ->Write(); // TH1D("DeltaX_0_2_45","DeltaX_0_2_45",100,-50.,50.);
  hDDDy_0_2_45 ->Write(); // TH1D("DeltaY_0_2_45","DeltaY_0_2_45",100,-50.,50.);

  hDDDx_0_2_50 ->Write(); // TH1D("DeltaX_0_2_50","DeltaX_0_2_50",100,-50.,50.);
  hDDDy_0_2_50 ->Write(); // TH1D("DeltaY_0_2_50","DeltaY_0_2_50",100,-50.,50.);





  hDDDx_0_3_01 ->Write(); // TH1D("DeltaX_0_3_01","DeltaX_0_3_01",100,-50.,50.);
  hDDDy_0_3_01 ->Write(); // TH1D("DeltaY_0_3_01","DeltaY_0_3_01",100,-50.,50.);

  hDDDx_0_3_12 ->Write(); // TH1D("DeltaX_0_3_12","DeltaX_0_3_12",100,-50.,50.);
  hDDDy_0_3_12 ->Write(); // TH1D("DeltaY_0_3_12","DeltaY_0_3_12",100,-50.,50.);

  hDDDx_0_3_23 ->Write(); // TH1D("DeltaX_0_3_23","DeltaX_0_3_23",100,-50.,50.);
  hDDDy_0_3_23 ->Write(); // TH1D("DeltaY_0_3_23","DeltaY_0_3_23",100,-50.,50.);

  hDDDx_0_3_34 ->Write(); // TH1D("DeltaX_0_3_34","DeltaX_0_3_34",100,-50.,50.);
  hDDDy_0_3_34 ->Write(); // TH1D("DeltaY_0_3_34","DeltaY_0_3_34",100,-50.,50.);

  hDDDx_0_3_45 ->Write(); // TH1D("DeltaX_0_3_45","DeltaX_0_3_45",100,-50.,50.);
  hDDDy_0_3_45 ->Write(); // TH1D("DeltaY_0_3_45","DeltaY_0_3_45",100,-50.,50.);

  hDDDx_0_3_50 ->Write(); // TH1D("DeltaX_0_3_50","DeltaX_0_3_50",100,-50.,50.);
  hDDDy_0_3_50 ->Write(); // TH1D("DeltaY_0_3_50","DeltaY_0_3_50",100,-50.,50.);





  hDDDx_0_4_01 ->Write(); // TH1D("DeltaX_0_4_01","DeltaX_0_4_01",100,-50.,50.);
  hDDDy_0_4_01 ->Write(); // TH1D("DeltaY_0_4_01","DeltaY_0_4_01",100,-50.,50.);

  hDDDx_0_4_12 ->Write(); // TH1D("DeltaX_0_4_12","DeltaX_0_4_12",100,-50.,50.);
  hDDDy_0_4_12 ->Write(); // TH1D("DeltaY_0_4_12","DeltaY_0_4_12",100,-50.,50.);

  hDDDx_0_4_23 ->Write(); // TH1D("DeltaX_0_4_23","DeltaX_0_4_23",100,-50.,50.);
  hDDDy_0_4_23 ->Write(); // TH1D("DeltaY_0_4_23","DeltaY_0_4_23",100,-50.,50.);

  hDDDx_0_4_34 ->Write(); // TH1D("DeltaX_0_4_34","DeltaX_0_4_34",100,-50.,50.);
  hDDDy_0_4_34 ->Write(); // TH1D("DeltaY_0_4_34","DeltaY_0_4_34",100,-50.,50.);

  hDDDx_0_4_45 ->Write(); // TH1D("DeltaX_0_4_45","DeltaX_0_4_45",100,-50.,50.);
  hDDDy_0_4_45 ->Write(); // TH1D("DeltaY_0_4_45","DeltaY_0_4_45",100,-50.,50.);

  hDDDx_0_4_50 ->Write(); // TH1D("DeltaX_0_4_50","DeltaX_0_4_50",100,-50.,50.);
  hDDDy_0_4_50 ->Write(); // TH1D("DeltaY_0_4_50","DeltaY_0_4_50",100,-50.,50.);









  hDDDx_1_0_01 ->Write(); // TH1D("Deltax_1_0_01","Deltax_1_0_01",100,-50.,50.);
  hDDDy_1_0_01 ->Write(); // TH1D("Deltay_1_0_01","Deltay_1_0_01",100,-50.,50.);

  hDDDx_1_0_12 ->Write(); // TH1D("Deltax_1_0_12","Deltax_1_0_12",100,-50.,50.);
  hDDDy_1_0_12 ->Write(); // TH1D("Deltay_1_0_12","Deltay_1_0_12",100,-50.,50.);

  hDDDx_1_0_23 ->Write(); // TH1D("Deltax_1_0_23","Deltax_1_0_23",100,-50.,50.);
  hDDDy_1_0_23 ->Write(); // TH1D("Deltay_1_0_23","Deltay_1_0_23",100,-50.,50.);

  hDDDx_1_0_34 ->Write(); // TH1D("Deltax_1_0_34","Deltax_1_0_34",100,-50.,50.);
  hDDDy_1_0_34 ->Write(); // TH1D("Deltay_1_0_34","Deltay_1_0_34",100,-50.,50.);

  hDDDx_1_0_45 ->Write(); // TH1D("Deltax_1_0_45","Deltax_1_0_45",100,-50.,50.);
  hDDDy_1_0_45 ->Write(); // TH1D("Deltay_1_0_45","Deltay_1_0_45",100,-50.,50.);

  hDDDx_1_0_50 ->Write(); // TH1D("Deltax_1_0_50","Deltax_1_0_50",100,-50.,50.);
  hDDDy_1_0_50 ->Write(); // TH1D("Deltay_1_0_50","Deltay_1_0_50",100,-50.,50.);





  hDDDx_1_1_01 ->Write(); // TH1D("Deltax_1_1_01","Deltax_1_1_01",100,-50.,50.);
  hDDDy_1_1_01 ->Write(); // TH1D("Deltay_1_1_01","Deltay_1_1_01",100,-50.,50.);

  hDDDx_1_1_12 ->Write(); // TH1D("Deltax_1_1_12","Deltax_1_1_12",100,-50.,50.);
  hDDDy_1_1_12 ->Write(); // TH1D("Deltay_1_1_12","Deltay_1_1_12",100,-50.,50.);

  hDDDx_1_1_23 ->Write(); // TH1D("Deltax_1_1_23","Deltax_1_1_23",100,-50.,50.);
  hDDDy_1_1_23 ->Write(); // TH1D("Deltay_1_1_23","Deltay_1_1_23",100,-50.,50.);

  hDDDx_1_1_34 ->Write(); // TH1D("Deltax_1_1_34","Deltax_1_1_34",100,-50.,50.);
  hDDDy_1_1_34 ->Write(); // TH1D("Deltay_1_1_34","Deltay_1_1_34",100,-50.,50.);

  hDDDx_1_1_45 ->Write(); // TH1D("Deltax_1_1_45","Deltax_1_1_45",100,-50.,50.);
  hDDDy_1_1_45 ->Write(); // TH1D("Deltay_1_1_45","Deltay_1_1_45",100,-50.,50.);

  hDDDx_1_1_50 ->Write(); // TH1D("Deltax_1_1_50","Deltax_1_1_50",100,-50.,50.);
  hDDDy_1_1_50 ->Write(); // TH1D("Deltay_1_1_50","Deltay_1_1_50",100,-50.,50.);



  hDDDx_1_2_01 ->Write(); // TH1D("Deltax_1_2_01","Deltax_1_2_01",100,-50.,50.);
  hDDDy_1_2_01 ->Write(); // TH1D("Deltay_1_2_01","Deltay_1_2_01",100,-50.,50.);

  hDDDx_1_2_12 ->Write(); // TH1D("Deltax_1_2_12","Deltax_1_2_12",100,-50.,50.);
  hDDDy_1_2_12 ->Write(); // TH1D("Deltay_1_2_12","Deltay_1_2_12",100,-50.,50.);

  hDDDx_1_2_23 ->Write(); // TH1D("Deltax_1_2_23","Deltax_1_2_23",100,-50.,50.);
  hDDDy_1_2_23 ->Write(); // TH1D("Deltay_1_2_23","Deltay_1_2_23",100,-50.,50.);

  hDDDx_1_2_34 ->Write(); // TH1D("Deltax_1_2_34","Deltax_1_2_34",100,-50.,50.);
  hDDDy_1_2_34 ->Write(); // TH1D("Deltay_1_2_34","Deltay_1_2_34",100,-50.,50.);

  hDDDx_1_2_45 ->Write(); // TH1D("Deltax_1_2_45","Deltax_1_2_45",100,-50.,50.);
  hDDDy_1_2_45 ->Write(); // TH1D("Deltay_1_2_45","Deltay_1_2_45",100,-50.,50.);

  hDDDx_1_2_50 ->Write(); // TH1D("Deltax_1_2_50","Deltax_1_2_50",100,-50.,50.);
  hDDDy_1_2_50 ->Write(); // TH1D("Deltay_1_2_50","Deltay_1_2_50",100,-50.,50.);





  hDDDx_1_3_01 ->Write(); // TH1D("Deltax_1_3_01","Deltax_1_3_01",100,-50.,50.);
  hDDDy_1_3_01 ->Write(); // TH1D("Deltay_1_3_01","Deltay_1_3_01",100,-50.,50.);

  hDDDx_1_3_12 ->Write(); // TH1D("Deltax_1_3_12","Deltax_1_3_12",100,-50.,50.);
  hDDDy_1_3_12 ->Write(); // TH1D("Deltay_1_3_12","Deltay_1_3_12",100,-50.,50.);

  hDDDx_1_3_23 ->Write(); // TH1D("Deltax_1_3_23","Deltax_1_3_23",100,-50.,50.);
  hDDDy_1_3_23 ->Write(); // TH1D("Deltay_1_3_23","Deltay_1_3_23",100,-50.,50.);

  hDDDx_1_3_34 ->Write(); // TH1D("Deltax_1_3_34","Deltax_1_3_34",100,-50.,50.);
  hDDDy_1_3_34 ->Write(); // TH1D("Deltay_1_3_34","Deltay_1_3_34",100,-50.,50.);

  hDDDx_1_3_45 ->Write(); // TH1D("Deltax_1_3_45","Deltax_1_3_45",100,-50.,50.);
  hDDDy_1_3_45 ->Write(); // TH1D("Deltay_1_3_45","Deltay_1_3_45",100,-50.,50.);

  hDDDx_1_3_50 ->Write(); // TH1D("Deltax_1_3_50","Deltax_1_3_50",100,-50.,50.);
  hDDDy_1_3_50 ->Write(); // TH1D("Deltay_1_3_50","Deltay_1_3_50",100,-50.,50.);





  hDDDx_1_4_01 ->Write(); // TH1D("Deltax_1_4_01","Deltax_1_4_01",100,-50.,50.);
  hDDDy_1_4_01 ->Write(); // TH1D("Deltay_1_4_01","Deltay_1_4_01",100,-50.,50.);

  hDDDx_1_4_12 ->Write(); // TH1D("Deltax_1_4_12","Deltax_1_4_12",100,-50.,50.);
  hDDDy_1_4_12 ->Write(); // TH1D("Deltay_1_4_12","Deltay_1_4_12",100,-50.,50.);

  hDDDx_1_4_23 ->Write(); // TH1D("Deltax_1_4_23","Deltax_1_4_23",100,-50.,50.);
  hDDDy_1_4_23 ->Write(); // TH1D("Deltay_1_4_23","Deltay_1_4_23",100,-50.,50.);

  hDDDx_1_4_34 ->Write(); // TH1D("Deltax_1_4_34","Deltax_1_4_34",100,-50.,50.);
  hDDDy_1_4_34 ->Write(); // TH1D("Deltay_1_4_34","Deltay_1_4_34",100,-50.,50.);

  hDDDx_1_4_45 ->Write(); // TH1D("Deltax_1_4_45","Deltax_1_4_45",100,-50.,50.);
  hDDDy_1_4_45 ->Write(); // TH1D("Deltay_1_4_45","Deltay_1_4_45",100,-50.,50.);

  hDDDx_1_4_50 ->Write(); // TH1D("Deltax_1_4_50","Deltax_1_4_50",100,-50.,50.);
  hDDDy_1_4_50 ->Write(); //     ("Deltay_1_4_50","Deltay_1_4_50",100,-50.,50.);






 

  







  theFile->Close();




  assert( ! _plus.empty()  );
  assert( 	!   _plus_prev.empty() ); 
  assert( !	   _plus_error.empty() );
 
  assert( !	   _minus.empty() );
  assert( !	   _minus_prev.empty() );
  assert( !	   _minus_error.empty() );
  assert (_plus.size() == _plus_prev.size() && _plus.size() == _plus_error.size() && _plus_prev.size() == _plus_error.size());
  assert (_minus.size() == _minus_prev.size() && _minus.size() == _minus_error.size() && _minus_prev.size() == _minus_error.size());
  

// fitting plus arm

// salvo i risultati dell'allineamento

// forcing not to align the following CSC
/*
    temp_Align_x[0][3][0]=0; temp_Align_y[0][3][0]=0;
    temp_Align_x[0][4][0]=0; temp_Align_y[0][4][0]=0;
    temp_Align_x[0][3][1]=0; temp_Align_y[0][3][1]=0;
    temp_Align_x[0][1][1]=0; temp_Align_y[0][1][1]=0;
    temp_Align_x[0][4][2]=0; temp_Align_y[0][4][2]=0;
    temp_Align_x[0][3][4]=0; temp_Align_y[0][3][4]=0;
*/


    FILE *file = fopen("sex_fit_results.txt","a");
    FILE *fileB = fopen("sex_NewT1Align.txt","w");

// leggo precedente file dell'allineamento
    double Align_x[2][6];
    double Align_y[2][6];
 
    for(int i=0; i<2; i++)
	for(int k=0; k<6; k++){
	  Align_x[i][k]=0;
	  Align_y[i][k]=0;

      }
    unsigned int a,c;
    float x_,y_;
    
    FILE *fileA = fopen("T1SexAlign.txt","r");
    if(fileA!=NULL){
      while(!feof(fileA)){
	
	fscanf(fileA,"%u%u%f%f",&a,&c,&x_,&y_);

	Align_x[a][c]=x_;
	Align_y[a][c]=y_;
      
      }
      fclose(fileA);
    }
 
		  

 if(0==0){
   cout << " START FITTING PLUS SIDE WITH " <<_plus.size() << " POINTS"<<endl;

  TFitterMinuit * minuit = new TFitterMinuit();
  MyFCN3 fcn3(_plus,_plus_prev,_plus_error);
  minuit->SetMinuitFCN(&fcn3);
  double startX = 0; 
  double startY = 0;
// SetParameter: name, param value, param error, min, max
// if not limited (vhigh <= vlow) 
  minuit->SetParameter(0,"x0",startX,0.01,-10,10); 
  minuit->SetParameter(1,"y0",startY,0.01,-10,10);
  minuit->SetParameter(2,"x1",startX,0.01,-10,10); 
  minuit->SetParameter(3,"y1",startY,0.01,-10,10);
  minuit->SetParameter(4,"x2",startX,0.01,-10,10); 
  minuit->SetParameter(5,"y2",startY,0.01,-10,10);
  minuit->SetParameter(6,"x3",startX,0.01,-10,10); 
  minuit->SetParameter(7,"y3",startY,0.01,-10,10);
  minuit->SetParameter(8,"x4",startX,0.01,-10,10); 
  minuit->SetParameter(9,"y4",startY,0.01,-10,10);
  minuit->SetParameter(10,"x5",0,0.01,-10,10);
  minuit->SetParameter(11,"y5",0,0.01,-10,10);

//		  minuit->SetParameter(2,"theta",0,0.1,0,0);
  minuit->SetPrintLevel(5);
// create Minimizer (default is Migrad)
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();

  if(iret == 0 ){
    cout << "PLUS: sextant alignment" << endl;
    for(unsigned int j=0; j<6; j++){
      fprintf(file,"%hd %hd",0,j);
      fprintf(fileB,"%hd %hd",0,j);

      fprintf(file," %g %g\n",minuit->GetParameter(2*j),minuit->GetParameter(2*j+1));
      fprintf(fileB," %g %g\n",minuit->GetParameter(2*j)+Align_x[0][j],minuit->GetParameter(2*j+1)+Align_y[0][j]);
    }
  double chi2, edm, errdef; 
  int nvpar, nparx;
  minuit->GetStats(chi2,edm,errdef,nvpar,nparx);  
  std::cout << "                         Chi2 Fit = " << chi2 <<std::endl;
  }else{
    cout << " PLUS: MINIMIZATION FAILED "<<endl; 
  }


 }
// fitting minus arm
 if(0==0){
   cout << " START FITTING MINUS SIDE WITH " <<_minus.size() << " POINTS"<<endl;

  TFitterMinuit * minuit = new TFitterMinuit();
  MyFCN3 fcn3(_minus,_minus_prev,_minus_error);
  minuit->SetMinuitFCN(&fcn3);
  double startX = 0; 
  double startY = 0;
// SetParameter: name, param value, param error, min, max
// if not limited (vhigh <= vlow) 
  minuit->SetParameter(0,"x0",startX,0.01,-10,10); 
  minuit->SetParameter(1,"y0",startY,0.01,-10,10);
  minuit->SetParameter(2,"x1",startX,0.01,-10,10); 
  minuit->SetParameter(3,"y1",startY,0.01,-10,10);
  minuit->SetParameter(4,"x2",startX,0.01,-10,10); 
  minuit->SetParameter(5,"y2",startY,0.01,-10,10);
  minuit->SetParameter(6,"x3",startX,0.01,-10,10); 
  minuit->SetParameter(7,"y3",startY,0.01,-10,10);
  minuit->SetParameter(8,"x4",startX,0.01,-10,10); 
  minuit->SetParameter(9,"y4",startY,0.01,-10,10);
  minuit->SetParameter(10,"x5",0,0.01,-10,10);
  minuit->SetParameter(11,"y5",0,0.01,-10,10);

//		  minuit->SetParameter(2,"theta",0,0.1,0,0);
  minuit->SetPrintLevel(5);
// create Minimizer (default is Migrad)
  minuit->CreateMinimizer();
  int iret = minuit->Minimize();

  if(iret == 0 ){
    cout << "MINUS: sextant alignment" << endl;
    for(unsigned int j=0; j<6; j++){
      fprintf(file,"%hd %hd",1,j);
      fprintf(fileB,"%hd %hd",1,j);

      fprintf(file," %g %g\n",minuit->GetParameter(2*j),minuit->GetParameter(2*j+1));
      fprintf(fileB," %g %g\n",minuit->GetParameter(2*j)+Align_x[1][j],minuit->GetParameter(2*j+1)+Align_y[1][j]);
    }
 
  double chi2, edm, errdef; 
  int nvpar, nparx;
  minuit->GetStats(chi2,edm,errdef,nvpar,nparx);  
  std::cout << "                         Chi2 Fit = " << chi2 <<std::endl;
  }else{
    cout << " MINUS: MINIMIZATION FAILED "<<endl; 
  }


 }

 fclose(fileB);
 fclose(file);


}




float T1Sex2::Eta(float x,float y,float z){
  float xyt;
  float c=0;
  float eta2;
  xyt = sqrt(x*x + y*y);
//theta
  if(z>0) c = atan(xyt/z);
  if(z<0) c = atan(xyt/z)+3.14159;
  if(z==0) {c = 3.14159;}
//pseudorapidity
  eta2 = -log(tan(c/2.));
  return eta2;
}
float T1Sex2::Phi(float x,float y){
  float c=0;
  if(x>0 && y>0) c = atan(y/x);
  if(x<0) c = atan(y/x)+3.14159;
  if(x>0 && y<0) c = atan(y/x)+6.28318;
  return c;
}

int T1Sex2::sextant(float x, float y, float z){
  int sestante = -1;
  float fi = Phi(x,y);
//  cout << " FI " << fi << endl;
  if(z>0){
    if(plane(z)==4){
      if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 5;
      if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 0;
      if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 1;
      if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 2;
      if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 3;
      if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 4;
    } 
    if(plane(z)==3){
      if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 5;
      if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 0;
      if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 1;
      if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 2;
      if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 3;
      if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 4;
    } 
    if(plane(z)==2){
      if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 5;
      if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
      if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
      if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 2;
      if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
      if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
    } 
    if(plane(z)==1){
      if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 5;
      if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 0;
      if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 1;
      if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 2;
      if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 3;
      if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 4;
    } 
    if(plane(z)==0){
      if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 5;
      if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
      if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
      if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 2;
      if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
      if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
    } 
    
  }


  if(z<0){
    if(plane(z)==4){
      if( fi < _30DEG || fi > (2*PI - _30DEG) ) sestante = 2;
      if(fi > (_60DEG -_30DEG) && fi < (_60DEG +_30DEG) ) sestante = 1;
      if(fi > (_120DEG -_30DEG) && fi < (_120DEG +_30DEG) ) sestante = 0;
      if(fi > (PI -_30DEG) && fi < (PI +_30DEG) ) sestante = 5;
      if(fi > (_240DEG -_30DEG) && fi < (_240DEG +_30DEG) ) sestante = 4;
      if(fi > (_300DEG -_30DEG) && fi < (_300DEG +_30DEG) ) sestante = 3;
    } 
    if(plane(z)==3){
      if( fi < (_30DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG) ) sestante = 2;
      if(fi > (_60DEG -_30DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG) ) sestante = 1;
      if(fi > (_120DEG -_30DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG) ) sestante = 0;
      if(fi > (PI -_30DEG+_3DEG) && fi < (PI +_30DEG+_3DEG) ) sestante = 5;
      if(fi > (_240DEG -_30DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG) ) sestante = 4;
      if(fi > (_300DEG -_30DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG) ) sestante = 3;
    } 
    if(plane(z)==2){
      if( fi < (_30DEG+_3DEG+_3DEG) || fi > (2*PI - _30DEG+_3DEG+_3DEG) ) sestante = 2;
      if(fi > (_60DEG -_30DEG+_3DEG+_3DEG) && fi < (_60DEG +_30DEG+_3DEG+_3DEG) ) sestante = 1;
      if(fi > (_120DEG -_30DEG+_3DEG+_3DEG) && fi < (_120DEG +_30DEG+_3DEG+_3DEG) ) sestante = 0;
      if(fi > (PI -_30DEG+_3DEG+_3DEG) && fi < (PI +_30DEG+_3DEG+_3DEG) ) sestante = 5;
      if(fi > (_240DEG -_30DEG+_3DEG+_3DEG) && fi < (_240DEG +_30DEG+_3DEG+_3DEG) ) sestante = 4;
      if(fi > (_300DEG -_30DEG+_3DEG+_3DEG) && fi < (_300DEG +_30DEG+_3DEG+_3DEG) ) sestante = 3;
    } 
    if(plane(z)==1){
      if( fi < (_30DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG) ) sestante = 2;
      if(fi > (_60DEG -_30DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG) ) sestante = 1;
      if(fi > (_120DEG -_30DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG) ) sestante = 0;
      if(fi > (PI -_30DEG-_3DEG) && fi < (PI +_30DEG-_3DEG) ) sestante = 5;
      if(fi > (_240DEG -_30DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG) ) sestante = 4;
      if(fi > (_300DEG -_30DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG) ) sestante = 3;
    } 
    if(plane(z)==0){
      if( fi < (_30DEG-_3DEG-_3DEG) || fi > (2*PI - _30DEG-_3DEG-_3DEG) ) sestante = 2;
      if(fi > (_60DEG -_30DEG-_3DEG-_3DEG) && fi < (_60DEG +_30DEG-_3DEG-_3DEG) ) sestante = 1;
      if(fi > (_120DEG -_30DEG-_3DEG-_3DEG) && fi < (_120DEG +_30DEG-_3DEG-_3DEG) ) sestante = 0;
      if(fi > (PI -_30DEG-_3DEG-_3DEG) && fi < (PI +_30DEG-_3DEG-_3DEG) ) sestante = 5;
      if(fi > (_240DEG -_30DEG-_3DEG-_3DEG) && fi < (_240DEG +_30DEG-_3DEG-_3DEG) ) sestante = 4;
      if(fi > (_300DEG -_30DEG-_3DEG-_3DEG) && fi < (_300DEG +_30DEG-_3DEG-_3DEG) ) sestante = 3;
    } 
    
  }

  
  return sestante;
}
int T1Sex2::plane(float z) {

  int piano = -1;
  if(fabs(z)<8000){
    piano = 0;
 
  }       
  if(fabs(z)<8400 && fabs(z)>8000){
    piano = 1;
  }
  if(fabs(z)<9000 && fabs(z)>8400){
    piano = 2;
  }
  if(fabs(z)<9500 && fabs(z)>9000){
    piano = 3;
  }
  if(fabs(z)>9500){
    piano = 4;
  }
  return piano;
 	
}
