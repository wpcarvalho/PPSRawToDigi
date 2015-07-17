/*
  Created by Fabrizio Ferro - INFN Genova for TOTEM
 */
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackAnalyzer.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "SimDataFormats/Track/interface/SimTrack.h"
#include "SimDataFormats/Track/interface/SimTrackContainer.h"
#include "SimDataFormats/Vertex/interface/SimVertex.h"
#include "SimDataFormats/Vertex/interface/SimVertexContainer.h"

#include "CLHEP/Vector/LorentzVector.h"

#include <vector>
#include <utility>
#include "TFile.h"
#include "TH1.h"

//#define _DEBUG_
#define T2_TOO
#ifndef PI
#define PI 3.141592653589793
#endif

T1TrackAnalyzer::T1TrackAnalyzer(const edm::ParameterSet& iConfig):_SeeTracks(0),_SeeHits(0),_ChiOverNCut(10000),_ZRange(1000),eeevvv(0),_CUTS(0)
{
	simVertexContainerLabel = iConfig.getParameter<edm::InputTag>("SimVertexContainerLabel");
	simTrackContainerLabel = iConfig.getParameter<edm::InputTag>("SimTrackContainerLabel");
	_SeeTracks = iConfig.getParameter<double>("SeeTrack");
	_SeeHits = iConfig.getParameter<double>("SeeHits");
	_ChiOverNCut = iConfig.getParameter<double>("ChiOverNCut");
	_ZRange =  iConfig.getParameter<double>("ZRange");
	_CUTS = iConfig.getParameter<int>("Cuts");
}


T1TrackAnalyzer::~T1TrackAnalyzer()
{
}


//
// member functions
//

// ------------ method called to for each event  ------------
void
T1TrackAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup )
{
	using namespace edm;

	double SingleTrackEvent_EtaRec=0;

	double SingleTrackEvent_PhiRec=0;

	int simulatedTraces=0;
	int reconstructedTracks=0;
	int reconstructedTracksMinus=0;
	int reconstructedTracksPlus=0;
	int reconstructedTracksT2=0;
	int reconstructedTracksT2Minus=0;
	int reconstructedTracksT2Plus=0;
	int reconstructedTracksZRange=0;
	int reconstructedTracksChiCut=0;
	int notReconstructedTracks=0;
	int GoodTracksNum=0;

	vector< pair<double,double> >::iterator pair_it;


	// get G4 Vertexes
	std::vector<SimVertex> theSimVertexes;

	Handle<SimVertexContainer> SimVtx;
	iEvent.getByLabel(simVertexContainerLabel, SimVtx);
	theSimVertexes.insert(theSimVertexes.end(),SimVtx->begin(),SimVtx->end());

	edm::Handle<T1T2TrackCollection> trackCollection;
	iEvent.getByLabel("t1tracks","T1TrackColl",trackCollection);

	int hjk=0;
	for(std::vector<SimVertex>::iterator isimvtx = theSimVertexes.begin(); isimvtx != theSimVertexes.end(); ++isimvtx){
#ifdef _DEBUG_
		cout << " Z Vertice["<<hjk<<"] = "<<isimvtx->position().z() << "  Parent indec = "<<isimvtx->parentIndex()<<endl;
#endif
		hjk++;
	}

	Handle<T1T2TrackCollection> trackCollectionT2;
	iEvent.getByLabel("T2TrackColl2","T2TrackColl",trackCollectionT2);


	//get G4 tracks
	vector< pair<double,double> > PrimarySimTracks;

	edm::Handle<edm::SimTrackContainer> G4TrkContainer;
	iEvent.getByLabel(simTrackContainerLabel, G4TrkContainer);
	if (!G4TrkContainer.isValid()) {
		edm::LogError("TrackerHitAnalyzer::analyze") << "Unable to find SimTrack in event!";
		return;
	}

	float Proton_p_max=0;

	edm::SimTrackContainer::const_iterator itTrk;



	nTr=0;
	for (itTrk = G4TrkContainer->begin(); itTrk != G4TrkContainer->end(); ++itTrk){

		double eta =0, phi =0, p =0;
		const CLHEP::HepLorentzVector G4Trk(itTrk->momentum().x(),itTrk->momentum().y(),itTrk->momentum().z(), itTrk->momentum().e() ) ;

		p =sqrt(G4Trk[0]*G4Trk[0]+G4Trk[1]*G4Trk[1]+G4Trk[2]*G4Trk[2]);
		if ( p == 0)
			edm::LogError("TrackerHitAnalyzer::analyze") << "TrackerTest::INFO: Primary has p = 0 ";
		else {
			double costheta  = G4Trk[2]/p;
			double theta = acos(TMath::Min(TMath::Max(costheta, -1.),1.));
			eta = -log(tan(theta/2));
			if ( G4Trk[0] != 0 || G4Trk[1] != 0) phi = atan2(G4Trk[1],G4Trk[0]);

			if(phi<0)phi = 2*PI + phi;

			float charge =  itTrk->charge();

			nTr++;
			eeTr[nTr]=eta;
			ffTr[nTr]=phi;
			chTr[nTr]=charge;
			tyTr[nTr]=itTrk->type();
			pTr[nTr]=p;


			evto[nTr] = eeevvv; //iEvent.id().event();

			if(itTrk->type()==2212){
				if(p>Proton_p_max){
					Proton_p_max = p;
				}
			}

		}
	}

	for (itTrk = G4TrkContainer->begin(); itTrk != G4TrkContainer->end(); ++itTrk) {

		double eta =0, phi =0, p =0;
		const CLHEP::HepLorentzVector G4Trk(itTrk->momentum().x(),itTrk->momentum().y(),itTrk->momentum().z(), itTrk->momentum().e() ) ;
		p =sqrt(G4Trk[0]*G4Trk[0]+G4Trk[1]*G4Trk[1]+G4Trk[2]*G4Trk[2]);
		if ( p == 0)
			edm::LogError("TrackerHitAnalyzer::analyze") << "TrackerTest::INFO: Primary has p = 0 ";
		else {
			double costheta  = G4Trk[2]/p;
			double theta = acos(TMath::Min(TMath::Max(costheta, -1.),1.));
			eta = -log(tan(theta/2));
			if ( G4Trk[0] != 0 || G4Trk[1] != 0) phi = atan2(G4Trk[1],G4Trk[0]);

			if(phi<0)phi = 2*PI + phi;


			float carica =  itTrk->charge();

			if(itTrk->charge() != 0  )
				hAllEtaSim -> Fill(eta);

			if(fabs(eta)>3. && fabs(eta)<5 && carica != 0){
				simulatedTraces++;

				hAllEtaSimT1Range -> Fill(eta);

#ifdef _DEBUG_
			int IndiceVerticeitTrk->vertIndex();
				std::cout << "TRACCIA G4: Eta = "<< eta<< "  Phi = " << phi << " Vert num. = "<<IndiceVertice <<  "   V = (" << theSimVertexes[IndiceVertice].position().x()<<","<<  theSimVertexes[IndiceVertice].position().y()<<","<<  theSimVertexes[IndiceVertice].position().z()<<")  PDG = " <<itTrk->type() << "  Charge = "<< itTrk->charge() << "   P = " << p
						<< std::endl;
#endif

				pair<double,double> *myprimarytrack = new pair<double,double>(eta,phi);
				PrimarySimTracks.push_back(*myprimarytrack);
			}
		}
	}


#ifdef _DEBUG_
	std::cout << " Size of track collection: " << trackCollection->size() <<std::endl;
#endif

	if (trackCollection->size() == 0) std::cout << "             WARNING: No tracks in T1 in the event !     " <<std::endl;

	reconstructedTracks=trackCollection->size();

	double DE = 1000;
	double DF = 1000;
	double DR = 1000;

	T1T2TrackCollection::const_iterator TC_it;
	if(_SeeTracks > 0){

		vector<int> *indiceTracce = new vector<int>;
		vector<int>::iterator pointerTrack_it;
		for(TC_it=trackCollection->begin(); TC_it!=trackCollection->end(); TC_it++){


			bool flag_ZRange = false;
			bool flag_ChiCut = false;
			DE = 1000; DF = 1000; DR = 1000;

#ifdef _DEBUG_
		int numero=0;
std::cout << "Track #"<<++numero;
std::cout << "  Track Size: " << (*TC_it).GetHitEntries() << std::endl;
#endif

			if(fabs((*TC_it).Z_at_Rmin() ) < _ZRange ){
				reconstructedTracksZRange++;
				flag_ZRange = true;
			}

			if((*TC_it).ChiSquaredOverN()<_ChiOverNCut){
				reconstructedTracksChiCut++;
				flag_ChiCut = true;
			}

			SingleTrackEvent_EtaRec= (*TC_it).Eta();
			SingleTrackEvent_PhiRec=(*TC_it).Phi();


			if(_CUTS>=1){
				if(SingleTrackEvent_EtaRec > 3 && SingleTrackEvent_EtaRec<5 && fabs((*TC_it).Z_at_Rmin() ) < _ZRange)reconstructedTracksPlus++;
				if(SingleTrackEvent_EtaRec < -3 && SingleTrackEvent_EtaRec>-5 && fabs((*TC_it).Z_at_Rmin() ) < _ZRange)reconstructedTracksMinus++;
			}else{
				if(SingleTrackEvent_EtaRec > 0)reconstructedTracksPlus++;
				if(SingleTrackEvent_EtaRec < 0)reconstructedTracksMinus++;
			}


			if(flag_ZRange)
				hAllEtaRecZRange -> Fill(SingleTrackEvent_EtaRec);

			if(flag_ChiCut)
				hAllEtaRecChiCut -> Fill(SingleTrackEvent_EtaRec);


			hAllEtaRec -> Fill(SingleTrackEvent_EtaRec);



#ifdef _DEBUG_
std::cout << "         XZ: x = ( "<<(*TC_it).GetTx()<<" +/- " << (*TC_it).GetTxSigma()<< " ) * z + ( "<< (*TC_it).X0()<<" +/- " << (*TC_it).X0Sigma()<< " )"<<std::endl;
std::cout << "         YZ: y = ( "<<(*TC_it).GetTy()<<" +/- " << (*TC_it).GetTySigma()<< " ) * z + ( "<< (*TC_it).Y0()<<" +/- " << (*TC_it).Y0Sigma()<< " )"<<std::endl;
std::cout << "         Eta = "<< (*TC_it).Eta() << " Phi = " << (*TC_it).Phi() << " Rmin = " << (*TC_it).Rmin() << " Z_at_Rmin = " << (*TC_it).Z_at_Rmin() <<std::endl;
#endif

			hChiSquaredOverN->Fill((*TC_it).ChiSquaredOverN() );

			if(_SeeHits > 0)
				for(unsigned int iy = 0; iy < (*TC_it).GetHitEntries(); iy++){
					std::cout << (*TC_it).GetHitT1(iy) <<std::endl;
				}

			// check compatibility with a G4 track
			int tmp_pointer=0;
			int my_pointer = 0;
			for(	pair_it=PrimarySimTracks.begin(); pair_it!=PrimarySimTracks.end(); pair_it++){
				tmp_pointer++;
				double DE_temp = SingleTrackEvent_EtaRec - (*pair_it).first;
				double DF_temp = SingleTrackEvent_PhiRec - (*pair_it).second;

				if(DF_temp>PI)DF_temp=2*PI-DF_temp;
				if(DF_temp<-PI)DF_temp=2*PI+DF_temp;
				if(DF_temp<PI && DF_temp>PI/2) DF_temp=-PI+DF_temp;
				if(DF_temp>-PI && DF_temp<-PI/2.)DF_temp=PI+DF_temp;
				double DR_temp = sqrt(DE_temp*DE_temp+DF_temp*DF_temp);

				if(DR_temp < DR){
					DR = DR_temp;
					DE=DE_temp;
					DF=DF_temp;
					my_pointer = tmp_pointer;

				}
			}

			if(DR<1000){
				hDeltaEta->Fill(DE);
				hDeltaPhi->Fill(DF);
				hDeltaR->Fill(DR);
				if(flag_ZRange){
					hDeltaEtaZRange->Fill(DE);
					hDeltaPhiZRange->Fill(DF);
					hDeltaRZRange->Fill(DR);
				}
				hDEvsCHIrid -> Fill ((*TC_it).ChiSquaredOverN(),DE);

				if((*TC_it).ChiSquaredOverN() < _ChiOverNCut){
					hDeltaEtaChiCut->Fill(DE);
					hDeltaPhiChiCut->Fill(DF);
					hDeltaRChiCut->Fill(DR);
				}
			}

			if(DE<0.5 && DF < 0.5){
				GoodTracksNum++;
				bool ffff = true;

				for(pointerTrack_it=indiceTracce->begin();pointerTrack_it!=indiceTracce->end();pointerTrack_it++){
					if( (*pointerTrack_it) == my_pointer && my_pointer != 0) ffff=false;
				}
				if(ffff) indiceTracce->push_back(my_pointer);
			}
		}

		hGoodTracks->Fill(GoodTracksNum);
		if(reconstructedTracks>0)
			hGoodTracksOverAllReco->Fill((float)GoodTracksNum/(float)reconstructedTracks);

		notReconstructedTracks = PrimarySimTracks.size() - indiceTracce->size();
		hNotRecoTracks->Fill(notReconstructedTracks);

		if(simulatedTraces>0){
			hNotRecoTracksRatio->Fill((float)notReconstructedTracks/(float)simulatedTraces);
			if((float)notReconstructedTracks/(float)simulatedTraces == 1)
				hNumOfTracksInEventsWnoGoodTracks->Fill(simulatedTraces);
		}

		if(reconstructedTracks == 0){
			hNumOfTracksInLostEvents->Fill(simulatedTraces);
			for(	pair_it=PrimarySimTracks.begin(); pair_it!=PrimarySimTracks.end(); pair_it++){

				hEtaOfTracksInLostEvents->Fill( fabs((*pair_it).first) );
			}
		}
		delete indiceTracce;
	}



	//T2
#ifdef T2_TOO

	T1T2TrackCollection::const_iterator TrkCit;
	double trketa=0.;
	for(TrkCit=trackCollectionT2->begin(); TrkCit!=trackCollectionT2->end(); TrkCit++){
		trketa= (*TrkCit).Eta();
		reconstructedTracksT2++;

		if(_CUTS>=1){


			double Z0xy=0;
			double R0xy=0;
			double ProbChi2X_xy=0; double ProbChi2Y_xy=0;
			R0xy=sqrt(((*TrkCit).bx_)*((*TrkCit).bx_)+((*TrkCit).by_)*((*TrkCit).by_));
			Z0xy=(*TrkCit).Z_at_Rmin();
			ProbChi2X_xy=TMath::Prob((*TrkCit).ChiSquaredX(),((*TrkCit).GetHitEntries()-2));
			ProbChi2Y_xy=TMath::Prob((*TrkCit).ChiSquaredY(),((*TrkCit).GetHitEntries()-2));

			if(trketa>0  && (fabs(Z0xy)<5000.)&&(R0xy<60.)&&(ProbChi2X_xy>0.01)&&(ProbChi2Y_xy>0.01)&&(TrkCit->GetHitEntries()>=4)


			) reconstructedTracksT2Plus++;

			if(trketa<0  && (fabs(Z0xy)<5000.)&&(R0xy<60.)&&(ProbChi2X_xy>0.01)&&(ProbChi2Y_xy>0.01)&&(TrkCit->GetHitEntries()>=4)


			) reconstructedTracksT2Minus++;

		}else{

			if(trketa>0) reconstructedTracksT2Plus++;
			if(trketa<0) reconstructedTracksT2Minus++;
		}
		hAllEtaRec->Fill(trketa);
	}

#endif

	// Left means Minus; Right means Plus

	// VETO on Plus Side
	if(reconstructedTracksT2Plus==0){

		hNTracksLeftHalfT1->Fill(reconstructedTracksT2Minus+reconstructedTracksMinus);
		hNTracksLeftNoT1->Fill(reconstructedTracksT2Minus);

	}
	if(reconstructedTracksT2Plus==0 && reconstructedTracksPlus==0){
		hNTracksLeft->Fill(reconstructedTracksT2Minus+reconstructedTracksMinus);
	}

	// VETO on Minus Side
	if(reconstructedTracksT2Minus==0 && reconstructedTracksMinus ==0){

		hNTracksRightHalfT1->Fill(reconstructedTracksT2Plus);
		hNTracksRight->Fill(reconstructedTracksT2Plus+reconstructedTracksPlus);

	}
	if(reconstructedTracksT2Minus==0){
		hNTracksRightNoT1->Fill(reconstructedTracksT2Plus);
	}



	hNTracks->Fill(reconstructedTracksT2Plus+reconstructedTracksPlus+reconstructedTracksT2Minus+reconstructedTracksMinus);
	hNTracksHalfT1->Fill(reconstructedTracksT2Plus+reconstructedTracksT2Minus+reconstructedTracksMinus);
	hNTracksNoT1->Fill(reconstructedTracksT2Plus+reconstructedTracksT2Minus);

	int EventType = 0;

	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus == 0) EventType = 0;
	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus > 0) EventType = 1;
	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus == 0) EventType = 2;
	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus > 0) EventType = 3;
	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus == 0) EventType = 4;
	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus > 0) EventType = 5;
	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus == 0) EventType = 6;
	if(reconstructedTracksT2Minus == 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus > 0) EventType = 7;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus == 0) EventType = 8;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus > 0) EventType = 9;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus == 0) EventType = 10;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus == 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus > 0) EventType = 11;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus == 0) EventType = 12;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus == 0 && reconstructedTracksT2Plus > 0) EventType = 13;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus == 0) EventType = 14;
	if(reconstructedTracksT2Minus > 0 &&  reconstructedTracksMinus > 0 && reconstructedTracksPlus > 0 && reconstructedTracksT2Plus > 0) EventType = 15;


	int Tr12 =  reconstructedTracksT2Plus +  reconstructedTracksT2Minus + reconstructedTracksPlus +  reconstructedTracksMinus;
	int Tr2 = reconstructedTracksT2Plus +  reconstructedTracksT2Minus;
	int TrHalf12 = reconstructedTracksT2Plus +  reconstructedTracksT2Minus +  reconstructedTracksMinus;




	hEventType->Fill(EventType);
	if(EventType==0){hTriggerType->Fill(0);}
	else if(EventType==1 ||EventType==2 ||EventType==3 || EventType==4 || EventType==8 ||EventType==12){
		hTriggerType->Fill(1);
		hTrackNumberSingleArmCondition12->Fill(Tr12);
	}
	else {
		hTriggerType->Fill(2);
		hTrackNumberSingleArmCondition12->Fill(Tr12);
		hTrackNumberDoubleArmCondition12->Fill(Tr12);
	}

	if(EventType == 0 || EventType == 2 )EventType = 0;
	if(EventType == 1 || EventType == 3 )EventType = 1;;
	if(EventType == 5 || EventType == 7 )EventType = 5;
	if(EventType == 6 || EventType == 4 )EventType = 4;
	if(EventType == 8 || EventType == 10 || EventType == 12)EventType = 8;
	if(EventType == 9 || EventType == 11 || EventType == 13)EventType = 9;
	if(EventType == 14 )EventType = 12;
	if(EventType == 15 )EventType = 13;
	hEventTypeHalfT1->Fill(EventType);
	if(EventType==0){hTriggerTypeHalfT1->Fill(0);

	}
	else if(EventType==1 ||EventType==2 ||EventType==3 || EventType==4 || EventType==8 ||EventType==12){
		hTriggerTypeHalfT1->Fill(1);
		hTrackNumberSingleArmConditionHalf12->Fill(TrHalf12);
	}
	else {
		hTriggerTypeHalfT1->Fill(2);
		hTrackNumberSingleArmConditionHalf12->Fill(TrHalf12);
		hTrackNumberDoubleArmConditionHalf12->Fill(TrHalf12);
	}

	if(EventType == 0 || EventType == 2 || EventType == 4 || EventType == 6 )EventType = 0;
	if(EventType == 1 || EventType == 3 || EventType == 5 || EventType == 7 )EventType = 1;
	if(EventType == 8 || EventType == 10 || EventType == 12 || EventType == 14 )EventType = 8;
	if(EventType == 9 || EventType == 11 || EventType == 13 || EventType == 15 )EventType = 9;
	hEventTypeNoT1->Fill(EventType);
	if(EventType==0){hTriggerTypeNoT1->Fill(0);}
	else if(EventType==1 ||EventType==2 ||EventType==3 || EventType==4 || EventType==8 ||EventType==12){
		hTriggerTypeNoT1->Fill(1);
		hTrackNumberSingleArmCondition2->Fill(Tr2);
		if(Tr2==0)cout << " AHIA " << EventType << endl;
	}
	else {
		hTriggerTypeNoT1->Fill(2);
		hTrackNumberSingleArmCondition2->Fill(Tr2);
		hTrackNumberDoubleArmCondition2->Fill(Tr2);
		if(Tr2==0)cout << " AHIA " << EventType << endl;
	}

	hSimTracks->Fill(simulatedTraces);
	hRecTracks->Fill(reconstructedTracks);
	hRecTracksZRange->Fill(reconstructedTracksZRange);
	hRecTracksChiCut->Fill(reconstructedTracksChiCut);
	if(simulatedTraces>0){
		hRecSimTrackRatio->Fill((float)reconstructedTracks/(float)simulatedTraces);
		hRecZRangeSimTrackRatio->Fill((float)reconstructedTracksZRange/(float)simulatedTraces);
	}else{
		hRecSimTrackRatio->Fill(1.);
		hRecZRangeSimTrackRatio->Fill(1.);
	}

#ifdef _DEBUG_
	cout << "Reconstructed tracks: " << reconstructedTracks << "    Simulated tracks " << simulatedTraces << endl;
#endif

	eeevvv++;

}


// ------------ method called once each job just before starting event loop  ------------
void 
T1TrackAnalyzer::beginJob()
{

	hDeltaEta = std::auto_ptr<TH1D>(new TH1D("deltaeta","deltaeta",320,-3.7,3.7));
	hDeltaPhi = std::auto_ptr<TH1D>(new TH1D("deltaphi","deltaphi",320,-3.7,3.7));
	hDeltaR = std::auto_ptr<TH1D>(new TH1D("deltaR","deltaR",320,0.,3.7));
	hDeltaEtaZRange = std::auto_ptr<TH1D>(new TH1D("deltaetaZRange","deltaetaZRange",320,-3.7,3.7));
	hDeltaPhiZRange = std::auto_ptr<TH1D>(new TH1D("deltaphiZRange","deltaphiZRange",320,-3.7,3.7));
	hDeltaRZRange = std::auto_ptr<TH1D>(new TH1D("deltaRZRange","deltaRZRange",320,0.,3.7));

	hDeltaEtaChiCut= std::auto_ptr<TH1D>(new TH1D("deltaetaChiCut","deltaetaChiCut",320,-3.7,3.7));
	hDeltaPhiChiCut = std::auto_ptr<TH1D>(new TH1D("deltaphiChiCut","deltaphiChiCut",320,-3.7,3.7));
	hDeltaRChiCut = std::auto_ptr<TH1D>(new TH1D("deltaRChiCut","deltaRChiCut",320,0.,3.7));

	hEventType = std::auto_ptr<TH1D>(new TH1D("EventType","EventType", 20, -0.5,19.5));
	hEventTypeNoT1 = std::auto_ptr<TH1D>(new TH1D("EventTypeNoT1","EventTypeNoT1", 20, -0.5,19.5));
	hEventTypeHalfT1 = std::auto_ptr<TH1D>(new TH1D("EventTypeHalfT1","EventTypeHalfT1", 20, -0.5,19.5));

	hTriggerType = std::auto_ptr<TH1D>(new TH1D("TriggerType","TriggerType", 5, -0.5,4.5));
	hTriggerTypeNoT1 = std::auto_ptr<TH1D>(new TH1D("TriggerTypeNoT1","TriggerTypeNoT1", 5, -0.5,4.5));
	hTriggerTypeHalfT1 = std::auto_ptr<TH1D>(new TH1D("TriggerTypeHalfT1","TriggerTypeHalfT1", 5, -0.5,4.5));

	hNTracks = std::auto_ptr<TH1D>(new TH1D("NTracks","NTracks", 20, -0.5,19.5));
	hNTracksNoT1 = std::auto_ptr<TH1D>(new TH1D("NTracksNoT1","NTracksNoT1", 20, -0.5,19.5));
	hNTracksHalfT1 = std::auto_ptr<TH1D>(new TH1D("NTracksHalfT1","NTracksHalfT1", 20, -0.5,19.5));
	hNTracksLeft = std::auto_ptr<TH1D>(new TH1D("NTracksLeft","NTracksLeft", 20, -0.5,19.5));
	hNTracksLeftNoT1 = std::auto_ptr<TH1D>(new TH1D("NTracksLeftNoT1","NTracksLeftNoT1", 20, -0.5,19.5));
	hNTracksLeftHalfT1 = std::auto_ptr<TH1D>(new TH1D("NTracksLeftHalfT1","NTracksLeftHalfT1", 20, -0.5,19.5));
	hNTracksRight = std::auto_ptr<TH1D>(new TH1D("NTracksRight","NTracksRight", 20, -0.5,19.5));
	hNTracksRightNoT1 = std::auto_ptr<TH1D>(new TH1D("NTracksRightNoT1","NTracksRightNoT1", 20, -0.5,19.5));
	hNTracksRightHalfT1 = std::auto_ptr<TH1D>(new TH1D("NTracksRightHalfT1","NTracksRightHalfT1", 20, -0.5,19.5));

	hSimTracks = std::auto_ptr<TH1D>(new TH1D("simTracks","simTracks",30,-0.5,29.5));
	hRecTracks = std::auto_ptr<TH1D>(new TH1D("recTracks","recTracks",30,-0.5,29.5));
	hRecTracksZRange = std::auto_ptr<TH1D>(new TH1D("recTracksZRange","recTracksZRange",30,-0.5,29.5));
	hRecTracksChiCut = std::auto_ptr<TH1D>(new TH1D("recTracksChiCut","recTracksChiCut",30,-0.5,29.5));
	hRecSimTrackRatio = std::auto_ptr<TH1D>(new TH1D("RecSimTrackRatio","RecSimTrackRatio",101,-0.05,10.05));
	hRecZRangeSimTrackRatio = std::auto_ptr<TH1D>(new TH1D("RecZRangeSimTrackRatio","RecZRangeSimTrackRatio",101,-0.05,10.05));
	hGoodTracks = std::auto_ptr<TH1D>(new TH1D("goodTracks","goodTracks",30,-0.5,29.5));
	hGoodTracksOverAllReco = std::auto_ptr<TH1D>(new TH1D("GoodReco/All","GoodReco/All",11,-0.05,1.05));
	hNotRecoTracks =  std::auto_ptr<TH1D>(new TH1D("notRecoTracks","notRecoTracks",30,-0.5,29.5));
	hNotRecoTracksRatio =  std::auto_ptr<TH1D>(new TH1D("notRecoTracksRatio","notRecoTracksRatio",11,-0.05,1.05));

	hNumOfTracksInEventsWnoGoodTracks = std::auto_ptr<TH1D>(new TH1D("NumOfSimTracksInEventsWnoGoodTracks","NumOfSimTracksInEventsWnoGoodTracks",30,-0.5,29.5));
	hNumOfTracksInLostEvents = std::auto_ptr<TH1D>(new TH1D("NumOfSimTracksInLostEvents","NumOfSimTracksInLostEvents",30,-0.5,29.5));
	hEtaOfTracksInLostEvents = std::auto_ptr<TH1D>(new TH1D("EtaOfSimTracksInLostEvents","EtaOfSimTracksInLostEvents",50,2,6));
	hAllEtaSim = std::auto_ptr<TH1D>(new TH1D("AllChEtaSim","AllChEtaSim",100,-15,15));
	hAllEtaSimT1Range = std::auto_ptr<TH1D>(new TH1D("AllChEtaSimT1Range","AllChEtaSimT1Range",100,-15,15));
	hAllEtaRec = std::auto_ptr<TH1D>(new TH1D("AllChEtaRec","AllChEtaRec",100,-15,15));
	hAllEtaRecZRange = std::auto_ptr<TH1D>(new TH1D("AllChEtaRecZRange","AllChEtaRecZRange",100,-15,15));
	hAllEtaRecChiCut = std::auto_ptr<TH1D>(new TH1D("AllChEtaRecChiCut","AllChEtaRecChiCut",100,-15,15));

	hDEvsCHIrid =  std::auto_ptr<TH2D>(new TH2D("DEvsCHIrid","DEvsCHIrid",50,0,10,100,-4,4));

	hChiSquaredOverN = std::auto_ptr<TH1D>(new TH1D("ChiSquaredOverN","ChiSquaredOverN",100,0,30));


	hTrackNumberSingleArmCondition2= std::auto_ptr<TH1D> (new TH1D("hTrackNumberSingleArmCondition2","hTrackNumberSingleArmCondition2",50,-0.5,49.5));
	hTrackNumberDoubleArmCondition2=  std::auto_ptr<TH1D> (new TH1D("hTrackNumberDoubleArmCondition2","hTrackNumberDoubleArmCondition2",50,-0.5,49.5));
	hTrackNumberSingleArmCondition12= std::auto_ptr<TH1D> (new TH1D("hTrackNumberSingleArmCondition12","hTrackNumberSingleArmCondition12",50,-0.5,49.5));
	hTrackNumberDoubleArmCondition12= std::auto_ptr<TH1D> (new TH1D("hTrackNumberDoubleArmCondition12","hTrackNumberDoubleArmCondition12",50,-0.5,49.5));
	hTrackNumberSingleArmConditionHalf12= std::auto_ptr<TH1D> (new TH1D("hTrackNumberSingleArmConditionHalf12","hTrackNumberSingleArmConditionHalf12",50,-0.5,49.5));
	hTrackNumberDoubleArmConditionHalf12= std::auto_ptr<TH1D> (new TH1D("hTrackNumberDoubleArmConditionHalf12","hTrackNumberDoubleArmConditionHalf12",50,-0.5,49.5));

}

// ------------ method called once each job just after ending the event loop  ------------
void 
T1TrackAnalyzer::endJob() {

	theFile = TFile::Open("tracksFile.root","RECREATE");
	if(!theFile || !theFile->IsWritable())
	{
		std::cout<<"Output file not opened correctly!!"<<std::endl;
	}

	hDeltaEta -> Write();
	hDeltaPhi -> Write();
	hDeltaR -> Write();
	hDeltaEtaZRange -> Write();
	hDeltaPhiZRange -> Write();
	hDeltaRZRange -> Write();

	hDeltaEtaChiCut-> Write();
	hDeltaPhiChiCut -> Write();
	hDeltaRChiCut -> Write();

	hEventType -> Write();
	hEventTypeNoT1 -> Write();
	hEventTypeHalfT1 -> Write();

	hTriggerType -> Write();
	hTriggerTypeNoT1 -> Write();
	hTriggerTypeHalfT1 -> Write();

	hNTracks -> Write();
	hNTracksNoT1 -> Write();
	hNTracksHalfT1 -> Write();
	hNTracksLeft -> Write();
	hNTracksLeftNoT1 -> Write();
	hNTracksLeftHalfT1 -> Write();
	hNTracksRight -> Write();
	hNTracksRightNoT1 -> Write();
	hNTracksRightHalfT1 -> Write();

	hSimTracks -> Write();
	hRecTracks -> Write();
	hRecTracksZRange -> Write();
	hRecTracksChiCut -> Write();
	hRecSimTrackRatio -> Write();
	hRecZRangeSimTrackRatio -> Write();
	hGoodTracks -> Write();
	hGoodTracksOverAllReco -> Write();
	hNotRecoTracks -> Write();
	hNotRecoTracksRatio -> Write();

	hNumOfTracksInEventsWnoGoodTracks -> Write();
	hNumOfTracksInLostEvents -> Write();
	hEtaOfTracksInLostEvents -> Write();
	hAllEtaSim -> Write();
	hAllEtaSimT1Range -> Write();
	hAllEtaRec -> Write();
	hAllEtaRecZRange -> Write();
	hAllEtaRecChiCut -> Write();

	hDEvsCHIrid -> Write();

	hChiSquaredOverN -> Write();

	hTrackNumberSingleArmCondition2->Write();
	hTrackNumberDoubleArmCondition2->Write();
	hTrackNumberSingleArmCondition12->Write();
	hTrackNumberDoubleArmCondition12->Write();
	hTrackNumberSingleArmConditionHalf12->Write();
	hTrackNumberDoubleArmConditionHalf12->Write();

	theFile->Close();
}
