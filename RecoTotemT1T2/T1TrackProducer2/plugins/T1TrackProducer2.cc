/*
	Created by Fabrizio Ferro - INFN Genova for TOTEM
	Optimized by Marcin Borratynski - AGH University of Science and Technology in Cracow
*/

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/ESHandle.h"

#include "DataFormats/T1DigiWire/interface/T1DigiWireCollection.h"
#include "DataFormats/T1Cluster/interface/T1Cluster.h"
#include "DataFormats/T1Cluster/interface/T1ClusterCollection.h"

#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "RecoTotemT1T2/T1TrackProducer2/interface/T1TrackProducer2.h"

#include <string>
#include <sys/time.h>
#include <unistd.h>

#define NUMBER_OF_PLANES 5

using namespace edm;
using namespace std;


T1TrackProducer2::T1TrackProducer2(const ParameterSet& config) :
	theChiRidThr(config.getParameter<double> ("ChiRidThr")),
	theVerbosity(config.getParameter<int> ("Verbosity"))
{
	t1RoadCollectionLabel = config.getParameter<edm::InputTag>("T1RoadCollectionLabel");
	produces<T1T2TrackCollection> ("T1TrackColl");

	if (theVerbosity >= 2) {
		std::cout << "Inside T1TrackProducer2" << std::endl;
	}
}

T1TrackProducer2::~T1TrackProducer2() {
}

void T1TrackProducer2::produce(Event& event, const EventSetup& setup) {

	int track_number = 0;
	auto_ptr <T1T2TrackCollection> trackCollection(new T1T2TrackCollection());
	edm::Handle <T1RoadCollection> roadCollection;

	if (theVerbosity >= 2) {
		std::cout << "*************************************************************" << std::endl;
		std::cout << "\t\tT1TrackProducer2::produce" << std::endl;
		std::cout << "*************************************************************" << std::endl;
	}
	// Get hits from the event
	event.getByLabel(t1RoadCollectionLabel, roadCollection);

	if (theVerbosity >= 2) {
		std::cout << "Road collection size: " << roadCollection->size() << std::endl;
	}


	for (T1RoadCollection::const_iterator RC_it = roadCollection->begin(); RC_it != roadCollection->end(); RC_it++) {
		if (theVerbosity >= 2) {
			std::cout << "\tRoad size: " << (*RC_it).size() << std::endl;
		}

		// For each road should build a selection of 3, 4 or 5 points, all on different levels(planes), fit to reconstruct the track.
		track_number += LookForTracks(*RC_it, *trackCollection);
	}


	event.put(trackCollection, "T1TrackColl");

	if (theVerbosity >= 1) {
		std::cout << "Final number of fitted tracks in event: " << track_number << std::endl;
	}
}

int T1TrackProducer2::LookForTracks(T1Road road,	T1T2TrackCollection &track_collection) {

	int number_of_tracks = 0;

	// one list for arm of the detector.
	list<int> index_per_plane_positive[5];
	list<int> index_per_plane_negative[5];

	for (int plane_index = 0; plane_index < NUMBER_OF_PLANES; plane_index++){
		index_per_plane_positive[plane_index] = list<int>();
		index_per_plane_negative[plane_index] = list<int>();
	}

	// save the indexes of hits in road corresponding to
	// plane 1 in index_per_plane[0],
	// plane 2 in index_per_plane[1],
	// ...
	// plane 5 in index_per_plane[4]

	for (unsigned int index_in_road = 0; index_in_road < road.size(); index_in_road++) {

		float tmp_z = road[index_in_road].GlobalPosition().z();


		// Z values - from T1geometry :
		//	{7521.4, 7569.6, 7526.4, 7569.6, 7521.4, 7564.6},	7400 8000
		//	{8191.4, 8239.4, 8196.4, 8239.6, 8191.4, 8234.6},	8000 8600
		//	{8798.4, 8846.6, 8803.4, 8846.6, 8798.4, 8841.6},	8600 9200
		//	{9405.4, 9453.6, 9410.4, 9453.6, 9405.4, 9448.6},	9200 9800
		//	{10168.4,10216.6,10173.4,10216.6,10168.4,10211.6}};	9800 10400

		for(int plane = 0; plane < NUMBER_OF_PLANES; plane ++){
			if ((tmp_z > 7400 + (plane * 600)) && (tmp_z < 8000 + (plane * 600))  ){

				index_per_plane_positive[plane].push_back(index_in_road);

			}else if ((tmp_z < -7400 - (plane * 600)) && (tmp_z > -8000 - (plane * 600))){

				index_per_plane_negative[plane].push_back(index_in_road);
			}
		}
	}

	//find fits for 5,4 and 3 hits. Order of the execution is important. Do not change it.
	int number_5_plane_tracks_positive = Fit5PlaneTrack(index_per_plane_positive, track_collection, road);
	int number_5_plane_tracks_negative = Fit5PlaneTrack(index_per_plane_negative, track_collection, road);

	int number_4_plane_tracks_positive = Fit4PlaneTrack(index_per_plane_positive, track_collection, road);
	int number_4_plane_tracks_negative = Fit4PlaneTrack(index_per_plane_negative, track_collection, road);

	int number_3_plane_tracks_positive = Fit3PlaneTrack(index_per_plane_positive, track_collection, road);
	int number_3_plane_tracks_negative = Fit3PlaneTrack(index_per_plane_negative, track_collection, road);

	if (theVerbosity >= 2) {
		std::cout << "\t\tNumber of 5 plane tracks - positive z value : " << number_5_plane_tracks_positive << std::endl;
		std::cout << "\t\tNumber of 5 plane tracks - negative z value : " << number_5_plane_tracks_negative << std::endl;
		std::cout << "\t\tNumber of 4 plane tracks - positive z value : " << number_4_plane_tracks_positive << std::endl;
		std::cout << "\t\tNumber of 4 plane tracks - negative z value : " << number_4_plane_tracks_negative << std::endl;
		std::cout << "\t\tNumber of 3 plane tracks - positive z value : " << number_3_plane_tracks_positive << std::endl;
		std::cout << "\t\tNumber of 3 plane tracks - negative z value : " << number_3_plane_tracks_negative << std::endl;
	}

	number_of_tracks =
			number_5_plane_tracks_positive + number_5_plane_tracks_negative +
			number_4_plane_tracks_positive + number_4_plane_tracks_negative +
			number_3_plane_tracks_positive + number_3_plane_tracks_negative ;

	return number_of_tracks;

}


int T1TrackProducer2::Fit5PlaneTrack(list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road){

	double chi_square_over_n;
	list <Hits_pointer> track_candidates;
	T1RecHitGlobal hit_vector[5];
	int hits_pointer_values[5];
	double values[10];

	// check all possible combinations for 5 planes and one hit on every plane.
	for (list<int>:: iterator iter_0 = index_per_plane[0].begin(); iter_0 != index_per_plane[0].end(); iter_0++)
		for (list<int>:: iterator iter_1 = index_per_plane[1].begin(); iter_1 != index_per_plane[1].end(); iter_1++)
			for (list<int>:: iterator iter_2 = index_per_plane[2].begin(); iter_2 != index_per_plane[2].end(); iter_2++)
				for (list<int>:: iterator iter_3 = index_per_plane[3].begin(); iter_3 != index_per_plane[3].end(); iter_3++)
					for (list<int>:: iterator iter_4 = index_per_plane[4].begin(); iter_4 != index_per_plane[4].end(); iter_4++){

						hits_pointer_values[0] = *iter_0;
						hits_pointer_values[1] = *iter_1;
						hits_pointer_values[2] = *iter_2;
						hits_pointer_values[3] = *iter_3;
						hits_pointer_values[4] = *iter_4;

						for(int i = 0; i < 5; ++i){
							hit_vector[i] = road[hits_pointer_values[i]];
						}

						//compute parameters and store result in values array.
						ComputeParametersOfTrack(hit_vector, 5, values);

						//chi_square_over_n is equal (chi_square_x + chi_square_y) / (float)(2*(track_hits_vector.size()-2)
						chi_square_over_n = (values[0] + values[1]) / 6.0;

						if (theVerbosity >= 3) cout << "\t\tchi2 value of 5 hits track = " <<  chi_square_over_n << std::endl;

						//store the candidate
						if (chi_square_over_n < theChiRidThr){
							track_candidates.push_back(	Hits_pointer(hits_pointer_values, 5, chi_square_over_n) );
						}
					}

	return ObtainTracksFromCandidates(track_candidates, index_per_plane, track_collection, road, 5);
}

int T1TrackProducer2::Fit4PlaneTrack(list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road){

	double chi_square_over_n;
	int number_of_combinations = 5;
	int hits_pointer_values[4];
	int combinations[5][4] = {
			{0,1,2,3},
			{0,1,2,4},
			{0,1,3,4},
			{0,2,3,4},
			{1,2,3,4}
	};
	double values[10];

	list <Hits_pointer> track_candidates;
	T1RecHitGlobal hit_vector[4];

	// check all possible combinations for 4 planes and one hit on every one.
	for(int index = 0; index < number_of_combinations; index++)
		for (list<int>:: iterator iter_0 = index_per_plane[combinations[index][0]].begin(); iter_0 != index_per_plane[combinations[index][0]].end(); iter_0++)
			for (list<int>:: iterator iter_1 = index_per_plane[combinations[index][1]].begin(); iter_1 != index_per_plane[combinations[index][1]].end(); iter_1++)
				for (list<int>:: iterator iter_2 = index_per_plane[combinations[index][2]].begin(); iter_2 != index_per_plane[combinations[index][2]].end(); iter_2++)
					for (list<int>:: iterator iter_3 = index_per_plane[combinations[index][3]].begin(); iter_3 != index_per_plane[combinations[index][3]].end(); iter_3++){

						hits_pointer_values[0] = *iter_0;
						hits_pointer_values[1] = *iter_1;
						hits_pointer_values[2] = *iter_2;
						hits_pointer_values[3] = *iter_3;

						for(int i = 0; i < 4; ++i){
							hit_vector[i] = road[hits_pointer_values[i]];
						}

						//compute parameters and store result in values array.
						ComputeParametersOfTrack(hit_vector, 4, values);

						//chi_square_over_n is equal (chi_square_x + chi_square_y) / (float)(2*(track_hits_vector.size()-2)
						chi_square_over_n = (values[0] + values[1]) / 4.0 ;

						if (theVerbosity >= 3) cout << "\t\tchi2 value of 4 hits track = " <<  chi_square_over_n << std::endl;

						//store the candidate
						if (chi_square_over_n < theChiRidThr){
							track_candidates.push_back(	Hits_pointer(hits_pointer_values, 4, chi_square_over_n) );
						}
					}
	return ObtainTracksFromCandidates(track_candidates, index_per_plane, track_collection, road, 4);
}

int T1TrackProducer2::Fit3PlaneTrack(list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road){

	double chi_square_over_n;
	int number_of_combinations = 10;
	int hits_pointer_values[3];
	int combinations[10][3] = {
			{0,1,2},
			{0,1,3},
			{0,1,4},
			{0,2,3},
			{0,2,4},
			{0,3,4},
			{1,2,3},
			{1,2,4},
			{1,3,4},
			{2,3,4},
	};
	double values[10];

	list <Hits_pointer> track_candidates;
	T1RecHitGlobal hit_vector[3];

	// check all possible combinations for 3 planes and one hit on every one.
	for(int index = 0; index < number_of_combinations; index++)
		for (list<int>:: iterator iter_0 = index_per_plane[combinations[index][0]].begin(); iter_0 != index_per_plane[combinations[index][0]].end(); iter_0++)
			for (list<int>:: iterator iter_1 = index_per_plane[combinations[index][1]].begin(); iter_1 != index_per_plane[combinations[index][1]].end(); iter_1++)
				for (list<int>:: iterator iter_2 = index_per_plane[combinations[index][2]].begin(); iter_2 != index_per_plane[combinations[index][2]].end(); iter_2++){

					hits_pointer_values[0] = *iter_0;
					hits_pointer_values[1] = *iter_1;
					hits_pointer_values[2] = *iter_2;

					for(int i = 0; i < 3; ++i){
						hit_vector[i] = road[hits_pointer_values[i]];
					}

					//compute parameters and store result in values array.
					ComputeParametersOfTrack(hit_vector,3, values);

					//chi_square_over_n is equal (chi_square_x + chi_square_y) / (float)(2*(track_hits_vector.size()-2)
					chi_square_over_n = (values[0] + values[1]) / 2.0;

					if (theVerbosity >= 3) std::cout << "\t\tchi2 value of 3 hits track = " <<  chi_square_over_n << std::endl;

					//store the candidate
					if ( chi_square_over_n  < theChiRidThr){
						track_candidates.push_back(	Hits_pointer(hits_pointer_values, 3, chi_square_over_n) );
					}
				}
	return ObtainTracksFromCandidates(track_candidates, index_per_plane, track_collection, road, 3);
}

int T1TrackProducer2::ObtainTracksFromCandidates(list<Hits_pointer> track_candidates, list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road, int size){

	int number_of_tracks = 0;
	set<int> used_values;
	T1RecHitGlobal hits[5];

	// fits with the lowest chi value first.
	track_candidates.sort();

	int ctr = 1;

	if (theVerbosity >= 2 && (!track_candidates.empty())) {
		std::cout <<  "\t\t----------- Best candidates for track size : " << size << " -----------"  << std::endl;
		for(list<Hits_pointer>:: iterator iter=track_candidates.begin(); iter!=track_candidates.end(); iter++){
			std::cout <<  "\t\ttrack candidate: " << ctr << ", chi2:  " <<  (*iter).chi() << std::endl;
			ctr++;
		}
		std::cout <<  "\t\t----------- Fitted candidates -----------" << std::endl;
	}

	for(list<Hits_pointer>:: iterator iter=track_candidates.begin(); iter!=track_candidates.end(); iter++){
		Hits_pointer candidate = *iter;
		bool accept_track = true;

		//check if candidate does not contain already used indexes.
		for(int iter = 0; iter < candidate.size(); iter++ )
			//check if set used_values contains index.
			if (used_values.find(candidate.index(iter)) != used_values.end())
				accept_track = false;

		//if all indexes are unused
		if(accept_track){

			if (theVerbosity >= 2 ){
				std::cout <<  "\t\t\t track chi: " <<  candidate.chi() << std::endl;
			}


			for(int hit_number = 0; hit_number < candidate.size(); hit_number++ ){
				//store result
				hits[hit_number] = road[candidate.index(hit_number)];
				//mark index as used
				used_values.insert(candidate.index(hit_number));

				if (theVerbosity >= 2) std::cout <<  "\t\t\t\t index of fitted hit in road: " << candidate.index(hit_number) << std::endl;

				//remove value from  index_per_plane list
				for(int i = 0; i < 5; ++i){
					index_per_plane[i].remove(candidate.index(hit_number));
				}

			}
			//create new track and add it to track collection
			track_collection.push_back(fitTrack(hits, candidate.size()));
			number_of_tracks++;
		}
	}
	return number_of_tracks;
}

void T1TrackProducer2::ComputeParametersOfTrack(T1RecHitGlobal* hit_vector, int size, double* output_values){

	float x[5], y[5], z[5], ex[5], ey[5], ez[5];
	float Sx = 0;
	float Sy = 0;
	float Sxz = 0;
	float Szz_x = 0;
	float Szz_y = 0;
	float Syz = 0;
	float Sz_x = 0;
	float Sz_y = 0;
	float S0_x = 0;
	float S0_y = 0;

	/**
	 * Chi2 is a variable that gives you an quantitative idea of how well a
	 * certain function matches a set of measured points.
	 * Function is a straight line in the space and the set of
	 * measured points are the subset of hits in the detector.
	 * A straight line in the space can be seen as two straight lines:
	 * the two projections in the xz and yz planes.
	 * The Chi2 of the fit in the yz plane is defined as the
	 * sum for i from 1 to N of (f(z(i)) - y(i))^2 / sigma(i)^2
	 * where N is the number of points to fit (3, 4 or 5),
	 * f(z) is the resulting straight line,
	 * (z(i),y(i)) are the hit coordinates and sigma is the uncertainty on the hit position.
	 * The final Chi2 is the sum of the Chi2's of the two projections.
	 */
	double chi2X = 0;
	double chi2Y = 0;

	for (int index = 0; (index < size) && (index < 5); index++) {
		x[index] = hit_vector[index].GlobalPosition().x();
		y[index] = hit_vector[index].GlobalPosition().y();
		z[index] = hit_vector[index].GlobalPosition().z();
		ex[index] = hit_vector[index].GlobalPositionError().cxx();
		ey[index] = hit_vector[index].GlobalPositionError().cyy();
		ez[index] = hit_vector[index].GlobalPositionError().czz();

		if (theVerbosity >= 3) {
			std::cout << "\t\t\t x = " << x[index] << " +/- " << sqrt(ex[index]) << std::endl;
			std::cout << "\t\t\t y = " << y[index] << " +/- " << sqrt(ey[index]) << std::endl;
			std::cout << "\t\t\t z = " << z[index] << " +/- " << sqrt(ez[index]) << std::endl;
		}

		Sxz += x[index] * z[index] / ex[index];
		Szz_x += z[index] * z[index] / ex[index];
		Sz_x += z[index] / ex[index];
		Szz_y += z[index] * z[index] / ey[index];
		Sz_y += z[index] / ey[index];
		Sx += x[index] / ex[index];
		S0_x += 1.0 / ex[index];
		S0_y += 1.0 / ey[index];
		Syz += y[index] * z[index] / ey[index];
		Sy += y[index] / ey[index];
	}

	// fit in XZ plane
	double a_xz = (Sxz * S0_x - Sz_x * Sx) / (Szz_x * S0_x - Sz_x * Sz_x);
	double b_xz = (Sx * Szz_x - Sz_x * Sxz) / (Szz_x * S0_x - Sz_x * Sz_x);
	double e_a_xz = sqrt(S0_x / (S0_x * Szz_x - Sz_x * Sz_x));
	double e_b_xz = sqrt(Szz_x / (S0_x * Szz_x - Sz_x * Sz_x));

	// fit in YZ plane
	double a_yz = (Syz * S0_y - Sz_y * Sy) / (Szz_y * S0_y - Sz_y * Sz_y);
	double b_yz = (Sy * Szz_y - Sz_y * Syz) / (Szz_y * S0_y - Sz_y * Sz_y);
	double e_a_yz = sqrt(S0_y / (S0_y * Szz_y - Sz_y * Sz_y));
	double e_b_yz = sqrt(Szz_y / (S0_y * Szz_y - Sz_y * Sz_y));

	for (int index = 0; index < size; index++) {
		if (theVerbosity >= 3) {
			std::cout << " Point " << x[index] << " " << y[index] << " " << z[index] << std::endl;
			std::cout << "\t\t\t aX=" << a_xz << " bX=" << b_xz << " eX=" << sqrt(ex[index]) << std::endl;
			std::cout << "\t\t\t aY=" << a_yz << " bY=" << b_yz << " eY=" << sqrt(ey[index]) << std::endl;
		}
		chi2X += (a_xz * z[index] + b_xz - x[index]) * (a_xz * z[index] + b_xz - x[index]) / ex[index];
		chi2Y += (a_yz * z[index] + b_yz - y[index]) * (a_yz * z[index] + b_yz - y[index]) / ey[index];
	}

	if (theVerbosity >= 3) {
		std::cout << " Track fit results: " << "a_x = " << a_xz << " +/- " << e_a_xz << std::endl;
		std::cout << "                    " << "b_x = " << b_xz << " +/- " << e_b_xz << std::endl;
		std::cout << "                    " << "a_y = " << a_yz << " +/- " << e_a_yz << std::endl;
		std::cout << "                    " << "b_y = " << b_yz << " +/- " << e_b_yz << std::endl;
		std::cout << "                    " << "Chi2= " << chi2X + chi2Y   << std::endl;
	}

	output_values[0] = chi2X;
	output_values[1] = chi2Y;

	output_values[2] = b_xz;
	output_values[3] = b_yz;
	output_values[4] = a_xz;
	output_values[5] = a_yz;

	output_values[6] = e_b_xz;
	output_values[7] = e_b_yz;
	output_values[8] = e_a_xz;
	output_values[9] = e_a_yz;
}

T1T2Track T1TrackProducer2::fitTrack(T1RecHitGlobal* hitvec, int size) {

	if (size != 3 && size != 4 && size != 5 ){
		cout <<  "ERROR T1TrackProducer2::fitTrack : Number of hits: " << size << " should be 3,4 or 5 " <<endl  ;
		return T1T2Track();
	}

	double values[10];
	int vector_size = 4;

	//compute parameters and store result in values array.
	ComputeParametersOfTrack(hitvec, size, values);

	// chi = chi_X + chi_Y
	double chi2 = values[0] + values[1];

	TVectorD vect(vector_size);

	vect[0] = values[2];
	vect[1] = values[3];
	vect[2] = values[4];
	vect[3] = values[5];

	TMatrixD mat(vector_size, vector_size);
	for (int i = 0; i < vector_size; i++)
		for (int j = 0; j < vector_size; j++)
			mat[i][j] = 0;

	//values only on diagonal
	mat[0][0] = values[6] * values[6];
	mat[1][1] = values[7] * values[7];
	mat[2][2] = values[8] * values[8];
	mat[3][3] = values[9] * values[9];

	// 1 or -1
	int hemisphere = (int)(hitvec[0].GlobalPosition().z()/fabs(hitvec[0].GlobalPosition().z()));

	T1T2Track track(vect, mat, chi2, values[0], values[1], hemisphere, 1);
	for (int i = 0; i < size; i++)
		track.AddHit(hitvec[i]);

	return track;
}

float T1TrackProducer2::Eta(float x, float y, float z) {
	float xyt;
	float c = 0.0;

	xyt = sqrt((x * x) + (y * y));
	//theta
	if (z > 0)
		c = atan(xyt / z);
	else if (z < 0)
		c = atan(xyt / z) + PI;
	else
		c = PI;
	//pseudorapidity
	return -log(tan(c / 2.));
}

float T1TrackProducer2::Phi(float x, float y) {
	float c = 0.0;

	if (x > 0 && y > 0)
		c = atan(y / x);
	else if (x < 0)
		c = atan(y / x) + PI;
	else if (x > 0 && y < 0)
		c = atan(y / x) + (2 * PI);

	return c;
}
