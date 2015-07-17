/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
  Modified by Marcin Borratynski
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackAnalyzerTB.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"
#include "RecoVertex/VertexTools/interface/FsmwModeFinder3d.h"

#include <vector>
#include "TFile.h"
#include "TH1.h"

#ifndef PI
#define PI 3.141592653589793
#endif

//#define _DEBUG_

T1TrackAnalyzerTB::T1TrackAnalyzerTB(const edm::ParameterSet& iConfig):_SeeTracks(0),_SeeHits(0),_ChiOverNCut(10000),_ZRange(1000),_Verbosity(0),eeevvv(0),_realBeamPosX(0),_realBeamPosY(0)
{
	_SeeTracks = iConfig.getParameter<double>("SeeTrack");
	_SeeHits = iConfig.getParameter<double>("SeeHits");
	_ChiOverNCut = iConfig.getParameter<double>("ChiOverNCut");
	_ZRange =  iConfig.getParameter<double>("ZRange");
	_Zmin =  iConfig.getParameter<double>("Zmin");
	_Zmax =  iConfig.getParameter<double>("Zmax");
	_Verbosity =  iConfig.getParameter<int>("Verbosity");
	_realBeamPosX = iConfig.getParameter<double>("realBeamPosX");
	_realBeamPosY = iConfig.getParameter<double>("realBeamPosY");
	_realBeamAngle = iConfig.getParameter<double>("realBeamAngle");


}


T1TrackAnalyzerTB::~T1TrackAnalyzerTB()
{
}

// ------------ method called to for each event  ------------
void
T1TrackAnalyzerTB::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
	using namespace edm;

	double SingleTrackEvent_EtaRec=0;

	int reconstructedTracks=0;
	int reconstructedTracksZRange=0;
	int reconstructedTracksChiCut=0;

	vector< pair<double,double> >::iterator pair_it;

	edm::Handle<T1T2TrackCollection> trackCollection;
	iEvent.getByLabel("t1tracks","T1TrackColl",trackCollection);

	edm::SimTrackContainer::const_iterator itTrk;

#ifdef _DEBUG_
	std::cout << " Taglia della track collection: " << trackCollection->size() <<std::endl;
#endif
	if(_Verbosity>0)
		if (trackCollection->size() == 0){ std::cout << " WARNING: No tracks in T1 in the event !      " <<std::endl;
		}
		else{
			std::cout << "\t\t\t\t\t\t" <<trackCollection->size() << " tracks in event " << eeevvv<<std::endl;
		}
	reconstructedTracks=trackCollection->size();

	T1T2TrackCollection::const_iterator TC_it;
	if(_SeeTracks > 0){

		for(TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++){
			allTracks->push_back((*TC_it));

			bool flag_ZRange = false;
			bool flag_ChiCut = false;

			double ax = tan( atan((*TC_it).GetTx()) + PI*_realBeamAngle/180.);

			double ay = (*TC_it).GetTy();
			double bx = (*TC_it).X0() - _realBeamPosX;
			double by = (*TC_it).Y0() - _realBeamPosY;

			double __Z_at_Rmin = -(ax*bx + ay*by)/(ax*ax + ay*ay);
			double __Rmin = sqrt( (ax*__Z_at_Rmin + bx)*(ax*__Z_at_Rmin + bx) + (ay*__Z_at_Rmin + by)*(ay*__Z_at_Rmin + by) );


			hRZ->Fill(__Rmin,__Z_at_Rmin);

#ifdef _DEBUG_
		int numero=0;
			std::cout << "Track #"<<++numero;
			std::cout << "  Track Size: " << (*TC_it).GetHitEntries() << std::endl;
#endif

			if(fabs(__Z_at_Rmin ) < _ZRange ){
				reconstructedTracksZRange++;
				flag_ZRange = true;
			}

			if((*TC_it).ChiSquaredOverN()<_ChiOverNCut){
				reconstructedTracksChiCut++;
				flag_ChiCut = true;
			}

			SingleTrackEvent_EtaRec= (*TC_it).Eta();

			if(flag_ZRange)
				hAllEtaRecZRange -> Fill(SingleTrackEvent_EtaRec);

			if(flag_ChiCut)
				hAllEtaRecChiCut -> Fill(SingleTrackEvent_EtaRec);

			hAllEtaRec -> Fill(SingleTrackEvent_EtaRec);


#ifdef _DEBUG_
			std::cout << "         XZ: x = ( "<<ax<<" +/- " << (*TC_it).GetTxSigma()<< " ) * z + ( "<< (*TC_it).X0()<<" +/- " << (*TC_it).X0Sigma()<< " )"<<std::endl;
			std::cout << "         YZ: y = ( "<<(*TC_it).GetTy()<<" +/- " << (*TC_it).GetTySigma()<< " ) * z + ( "<< (*TC_it).Y0()<<" +/- " << (*TC_it).Y0Sigma()<< " )"<<std::endl;
			std::cout << "         Eta = "<< (*TC_it).Eta() << " Phi = " << (*TC_it).Phi() << " Rmin = " << __Rmin << " Z_at_Rmin = " << __Z_at_Rmin <<std::endl;
#endif

			hXYatZ5000->Fill(ax*5000.+(*TC_it).X0()- _realBeamPosX,(*TC_it).GetTy()*5000.+(*TC_it).Y0()- _realBeamPosY );
			hXYatZ0->Fill((*TC_it).X0()- _realBeamPosX,(*TC_it).Y0()- _realBeamPosY );
			if(__Z_at_Rmin >4000 && __Z_at_Rmin<7000){
				hXYatZ5000cut->Fill(ax*5000.+(*TC_it).X0()- _realBeamPosX,(*TC_it).GetTy()*5000.+(*TC_it).Y0() - _realBeamPosY);
				hXYatZ5700cut->Fill(ax*5700.+(*TC_it).X0()- _realBeamPosX,(*TC_it).GetTy()*5700.+(*TC_it).Y0() - _realBeamPosY);
				hXYatZ6750cut->Fill(ax*6750.+(*TC_it).X0()- _realBeamPosX,(*TC_it).GetTy()*6750.+(*TC_it).Y0()- _realBeamPosY );
			}

			if(__Z_at_Rmin >-2000 && __Z_at_Rmin<2000){
				hXYatZm500cut->Fill(ax* -500.+(*TC_it).X0()- _realBeamPosX,(*TC_it).GetTy()* -500.+(*TC_it).Y0()- _realBeamPosY );
				hXYatZ0cut->Fill((*TC_it).X0()- _realBeamPosX,(*TC_it).Y0()- _realBeamPosY );


			}
			hChiSquaredOverN->Fill((*TC_it).ChiSquaredOverN() );

			if(_SeeHits > 0)
				for(unsigned int iy = 0; iy < (*TC_it).GetHitEntries(); iy++){
					std::cout << (*TC_it).GetHitT1(iy) <<std::endl;
				}
		}
	}

	hRecTracks->Fill(reconstructedTracks);
	hRecTracksZRange->Fill(reconstructedTracksZRange);
	hRecTracksChiCut->Fill(reconstructedTracksChiCut);

	eeevvv++;
}


// ------------ method called once each job just before starting event loop  ------------
void 
T1TrackAnalyzerTB::beginJob(const edm::EventSetup& iSetup)
{


	theFile = new TFile("tracksFile.root","RECREATE");

	hRecTracks = new TH1D("recTracks","recTracks",30,-0.5,29.5);
	hRecTracksZRange = new TH1D("recTracksZRange","recTracksZRange",30,-0.5,29.5);
	hRecTracksChiCut = new TH1D("recTracksChiCut","recTracksChiCut",30,-0.5,29.5);

	hGoodTracks = new TH1D("goodTracks","goodTracks",30,-0.5,29.5);
	hGoodTracksOverAllReco = new TH1D("GoodReco/All","GoodReco/All",11,-0.05,1.05);
	hNotRecoTracks =  new TH1D("notRecoTracks","notRecoTracks",30,-0.5,29.5);


	hNumOfTracksInEventsWnoGoodTracks = new TH1D("NumOfSimTracksInEventsWnoGoodTracks","NumOfSimTracksInEventsWnoGoodTracks",30,-0.5,29.5);
	hNumOfTracksInLostEvents = new TH1D("NumOfSimTracksInLostEvents","NumOfSimTracksInLostEvents",30,-0.5,29.5);
	hEtaOfTracksInLostEvents = new TH1D("EtaOfSimTracksInLostEvents","EtaOfSimTracksInLostEvents",50,2,6);

	hAllEtaRec = new TH1D("AllChEtaRec","AllChEtaRec",480,-15,15);
	hAllEtaRecZRange = new TH1D("AllChEtaRecZRange","AllChEtaRecZRange",480,-15,15);
	hAllEtaRecChiCut = new TH1D("AllChEtaRecChiCut","AllChEtaRecChiCut",480,-15,15);

	hDEvsCHIrid = new TH2D("DEvsCHIrid","DEvsCHIrid",50,0,10,100,-4,4);

	hChiSquaredOverN = new TH1D("ChiSquaredOverN","ChiSquaredOverN",100,0,30);

	hXYatZ5000 = new TH2D("XYatZ5000","XYatZ5000",500,-1000,1000,500,-1000,1000);
	hXYatZ0 = new TH2D("XYatZ0","XYatZ0",500,-1000,1000,500,-1000,1000);
	hXYatZ5000cut = new TH2D("XYatZ5000cut","XYatZ5000cut",500,-1000,1000,500,-1000,1000);
	hXYatZ5700cut = new TH2D("XYatZ5700cut","XYatZ5700cut",500,-1000,1000,500,-1000,1000);
	hXYatZ6750cut = new TH2D("XYatZ6750cut","XYatZ6750cut",500,-1000,1000,500,-1000,1000);
	hXYatZ0cut = new TH2D("XYatZ0cut","XYatZ0cut",500,-1000,1000,500,-1000,1000);

	hXYatZm500cut = new TH2D("XYatm500cut","XYatZm500cut",500,-1000,1000,500,-1000,1000);

	hRZ = new TH2D("RZ","RZ",500,-1000,1000,500,-20000,20000);

	allTracks = new T1T2TrackCollection();
}

// ------------ method called once each job just after ending the event loop  ------------
void 
T1TrackAnalyzerTB::endJob() {

	std::cout << allTracks->size() << std::endl;

	vector< pair< GlobalPoint, float > > V_PaD;
	T1T2TrackCollection::const_iterator TC_it1;
	T1T2TrackCollection::const_iterator TC_it2;
	int aa = 0;
	for(TC_it1=allTracks->begin(); TC_it1!=allTracks->end(); TC_it1++){
		std::cout << aa++ << std::endl;
		for(TC_it2=TC_it1+1; TC_it2 != allTracks->end(); TC_it2++){
			if( (*TC_it1).Z_at_Rmin() > _Zmin &&  (*TC_it2).Z_at_Rmin() > _Zmin && (*TC_it1).Z_at_Rmin() < _Zmax && (*TC_it2).Z_at_Rmin() < _Zmax)
			{
				pair< GlobalPoint, float >  coppia ;
				TwoTracksMinD((*TC_it1),(*TC_it2),coppia);
				V_PaD.push_back(coppia);
			}
		}
	}

	const vector< pair< GlobalPoint, float > > V_PaD_c = V_PaD;

	FsmwModeFinder3d Finder;
	GlobalPoint Vert = Finder(V_PaD_c);
	math::XYZPointD Vert_(Vert.x(),Vert.y(),Vert.z() );

	std::cout << "        RECO Vert "<<Vert << "  cm" <<std::endl;

	delete allTracks;

	theFile->Write();
	theFile->Close();
}





void T1TrackAnalyzerTB::TwoTracksMinD(T1T2Track t1, T1T2Track t2, pair<GlobalPoint, float> & coppia_){

	//  std::cout << "Inside TwoTracksMinD " <<std::endl;
	long double a1x ;
	long double b1x ;
	long double a1y;
	long double b1y ;

	long double a2x ;
	long double b2x ;
	long double a2y ;
	long double b2y ;

	if(t1.GetDet()==1 && t2.GetDet()==1){
		a1x = tan( atan(t1.GetTx()) + PI *_realBeamAngle/180.);
		b1x = t1.X0()- _realBeamPosX;
		a1y = t1.GetTy();
		b1y = t1.Y0()- _realBeamPosY;

		a2x = tan( atan(t2.GetTx()) + _realBeamAngle);
		b2x = t2.X0()- _realBeamPosX;
		a2y = t2.GetTy();
		b2y = t2.Y0()- _realBeamPosY;
	}else if(t1.GetDet()==2 && t2.GetDet()==2){
		a1x = t1.GetTy() * cos( t1.Phi() );
		b1x = t1.GetTx() * cos( t1.Phi() );
		a1y = t1.GetTy() * sin( t1.Phi() );
		b1y = t1.GetTx() * sin( t1.Phi() );

		a2x = t2.GetTy() * cos( t2.Phi() );
		b2x = t2.GetTx() * cos( t2.Phi() );
		a2y = t2.GetTy() * sin( t2.Phi() );
		b2y = t2.GetTx() * sin( t2.Phi() );

	}else if(t1.GetDet()==1 && t2.GetDet()==2){
		a1x = t1.GetTx();
		b1x = t1.X0();
		a1y = t1.GetTy();
		b1y = t1.Y0();

		a2x = t2.GetTy() * cos( t2.Phi() );
		b2x = t2.GetTx() * cos( t2.Phi() );
		a2y = t2.GetTy() * sin( t2.Phi() );
		b2y = t2.GetTx() * sin( t2.Phi() );

	}else if(t1.GetDet()==2 && t2.GetDet()==1){
		a2x = t2.GetTx();
		b2x = t2.X0();
		a2y = t2.GetTy();
		b2y = t2.Y0();

		a1x = t1.GetTy() * cos( t1.Phi() );
		b1x = t1.GetTx() * cos( t1.Phi() );
		a1y = t1.GetTy() * sin( t1.Phi() );
		b1y = t1.GetTx() * sin( t1.Phi() );
	}
	else{
		std::cout << "ERROR: wrong detector number in track. " << std::endl;
		assert(0);
	}

	long double A = 1 + a1x * a2x + a1y * a2y;
	long double B1 = a1x*b2x + a1y*b2y - a1x*b1x - a1y*b1y;
	long double B2 =  a2x*b1x + a2y*b1y - a2x*b2x - a2y*b2y;
	long double C1 = 1 + a1x*a1x + a1y*a1y;
	long double C2 = 1 + a2x*a2x + a2y*a2y;
	long double D = -2*b1x*b2x -2*b1y*b2y + b1x*b1x + b2x*b2x + b1y*b1y + b2y*b2y;
	long double S = (B2*C1 + A*B1)/(C1*C2-A*A);
	long double T = A/C1*S + B1/C1;

	long double dsquared = (S*S*(C2-(A*A/C1)))-2.*(S*(A*B1/C1+B2))+(D-(B1*B1/C1));

	//-----------------------------
	// new algorithm
/*
	double dsquared_2 = 0;
	double d_min =0 ;
	double AAAA=0; double AAAAAA=0; double mod_v_ort_sq=0;
	if(a1x != a2x){

		AAAA = (a2y - a1y)/(a2x-a1x);
		AAAAAA = (a1x*a2y - a2x*a1y)/(a2x-a1x);

		mod_v_ort_sq = AAAA*AAAA + AAAAAA*AAAAAA + 1;

		d_min = ( (b1x-b2x)*(a2y-a1y)/(a1x-a2x)+(b1y-b2y))/sqrt(mod_v_ort_sq);

		dsquared_2 = d_min * d_min;

	}

	else{
		d_min = sqrt( (b1x-b2x)*(b1x-b2x) )/ (1 + (a1x+a2x)*(a1x+a2x)/4.) ;
		dsquared_2 = d_min*d_min;
	}


	if(dsquared<0 && dsquared> -0.001) dsquared=0;
	assert(dsquared >= 0);
*/
	//----------------------------

	double d = sqrt(dsquared);

	double x1 = b1x+T*a1x;
	double y1 = b1y+T*a1y;
	double z1 = T;
	double x2 = b2x+S*a2x;
	double y2 = b2y+S*a2y;
	double z2 = S;

	if(_Verbosity>1){
		std::cout << dsquared << std::endl;

		std::cout << A << " " << B1<< " " << C1 << std::endl;

		std::cout << T << " " << S << std::endl;

		std::cout << " Lines: " <<std::endl;
		std::cout << "1: x = "<< a1x << " z + " << b1x
				<< "    y = "<< a1y << " z + " << b1y <<std::endl;
		std::cout << "2: x = "<< a2x << " z + " << b2x
				<< "    y = "<< a2y << " z + " << b2y <<std::endl;
		std::cout << "Minimum distance " << d << std::endl;
		std::cout << "between points:" << std::endl;
		std::cout << "1: "<< x1 << ", "<<  y1 << ", "<<  z1 << " "<< std::endl;
		std::cout << "2: "<< x2 << ", "<<  y2 << ", "<<  z2 << " "<< std::endl;
	}

	GlobalPoint InMezzo(x1/20.+x2/20.,y1/20.+y2/20.,z1/20.+z2/20.); //to cm

	coppia_.first = InMezzo;
	coppia_.second = d/10.;

	return;
}
