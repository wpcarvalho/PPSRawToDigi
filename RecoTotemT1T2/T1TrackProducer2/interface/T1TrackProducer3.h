//	Created by Fabrizio Ferro - INFN Genova for TOTEM
//	Modified and documented by Marcin Borratynski - AGH University of Science and Technology in Cracow

#ifndef _T1TrackProducer3_h
#define _T1TrackProducer3_h

#include "FWCore/Framework/interface/EDProducer.h"
#include "Geometry/TotemGeometry/interface/T1Geometry.h"

#include "DataFormats/T1DetId/interface/T1DetId.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2D.h"
#include "DataFormats/T1RecHit/interface/T1RecHit2DCollection.h"
#include "DataFormats/T1Road/interface/T1Road.h"
#include "DataFormats/T1T2Track/interface/T1T2TrackCollection.h"

namespace edm {
	class ParameterSet;
	class Event;
	class EventSetup;
}

/**
 *  Auxiliary class used to store the hit vector indexes in road and chi square value of the fit.
 *  Used by T1TrackProducer3 functions
 */
class Hits_pointer3 {
	public:

		/**
		 * A constructor.
		 * @param indexes is an integer pointer argument to the table of indexes in road of one hit vector.
		 * @param size is a size of the indexes array
		 * @param chi is the chi squared value of the fit
		 */
		Hits_pointer3(int* indexes, int size, double chi) {
			_chi = chi;
			_size = size;
			for(int i = 0; (i < size) && (i < 5); i ++){
				_indexes[i] = indexes[i];
			}
		};

		~Hits_pointer3() {};

		/**
		 * Is required by stl::list.stort() function.
		 */
		bool operator <(const Hits_pointer3& indici) {
			return _chi < indici.chi();
		}

		/**
		 * Function return the value of .
		 * @param i is an index in indexes table.
		 * @return Index of i-th element of the hit vector in road list. Road is a T1TrackProducer3::LookForTracks() member.
		 */
		int index(int i) const { return _indexes[i]; };

		/**
		 * Getter type function.
		 * @return Size of indexes array
		 */
		int size() const { return _size; };
		/**
		 * Getter type function.
		 * @return Chi squared value of the fit
		 */
		double chi() const { return _chi; };

	private:
		int _indexes[5];
		int _size;
		double _chi;
};

/**
 *  Obtains tracks form the roads in road collection.
 */
class T1TrackProducer3: public edm::EDProducer {
	public:
		/**
		 *  Constructor.
		 */
		T1TrackProducer3(const edm::ParameterSet&);

		/**
		 *  Function that compute Eta value form given (x,y,z)
		 *  @return -log(tan(theta / 2.0))
		 */
		float Eta(float x, float y, float z);

		/**
		 *  Function that compute Phi value form given (x, y)
		 */
		float Phi(float x, float y);

		virtual ~T1TrackProducer3();

		/**
		 *  The method which produces the reconstructed hits
		 *  @param event contains road collection. T1TrackProducer3::LookForTracks() is executed for every road in road collection.
		 *  @param setup in unused
		 */
		virtual void produce(edm::Event& event, const edm::EventSetup& setup);

	private:
		edm::InputTag t1RoadCollectionLabel;
		/**
		 *  Executed for every road in road collection. Use Fit5PlaneTrack, Fit4PlaneTrack and Fit3PlaneTrack method.
		 *  @param road in the bunch of hits.
		 *  @param track_collection is an output parameter. Best fitted tracks are stored in it.
		 *  @return number of fitted tracks.
		 */
		virtual int LookForTracks(T1Road road, T1T2TrackCollection& track_collection);
		/**
		 *  Function that looks for tracks on 5 planes of the detector. Function checks every possible fitting of a track and choose best solutions.
		 *  @param index_per_plane is a table of lists. index_per_plane[i] where (0 <= i <= 4) contains indexes of hits in road that are associated with the (i+1)-th plane of the detector.
		 *  @param track_collection is an output parameter. Best fitted tracks are stored in it.
		 *  @param road is a vector of T1RecHitGlobal. Many roads may be connected to one event.
		 *  @return number of fitted tracks. This is a sum of results returned by Fit5PlaneTrack(), Fit4PlaneTrack() and Fit3PlaneTrack().
		 */
		int Fit5PlaneTrack(list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road);
		/**
		 *  Function that looks for tracks on 4 planes of the detector. Function checks every possible fitting of a track and choose best solutions.
		 *  @param index_per_plane is a table of lists. index_per_plane[i] where (0 <= i <= 4) contains indexes of hits in road that are associated with the (i+1)-th plane of the detector.
		 *  @param track_collection is an output parameter. Best fitted tracks are stored in it.
		 *  @param road is a vector of T1RecHitGlobal. Many roads may be connected to one event.
		 *  @return number of fitted tracks.
		 */
		int Fit4PlaneTrack(list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road);
		/**
		 *  Function that looks for tracks on 3 planes of the detector. Function checks every possible fitting of a track and choose best solutions.
		 *  @param index_per_plane is a table of lists. index_per_plane[i] where (0 <= i <= 4) contains indexes of hits in road that are associated with the (i+1)-th plane of the detector.
		 *  @param track_collection is an output parameter. Best fitted tracks are stored in it.
		 *  @param road is a vector of T1RecHitGlobal. Many roads may be connected to one event.
		 *  @return number of fitted tracks.
		 */
		int Fit3PlaneTrack(list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road);
		/**
		 *  Method that picks up the best candidates of fitting.
		 *  @param track_candidates is a list of Hits_pointer3. It contain all tracks with chiSquared / 2(sizeOfVector-2) lower than theChiRidThr set up in config file.
		 *  @param index_per_plane is a table of lists. index_per_plane[i] where (0 <= i <= 4) contains indexes of hits in road that are associated with the (i+1)-th plane of the detector.
		 *  @param track_collection is an output parameter. Best fitted tracks are stored in it.
		 *  @param road is a vector of T1RecHitGlobal. Many roads may be connected to one event.
		 *  @param size is a size of hit vector stored in track_candidates. Used when verbocity > 0.
		 *  @return number of picked up candidates. That is equal to number of fitted tracks.
		 */
		int ObtainTracksFromCandidates(list<Hits_pointer3> track_candidates, list<int> index_per_plane[], T1T2TrackCollection& track_collection, T1Road road, int size);
		/**
		 *  Method that computes chi square for x plane fit, chi squared for y plane fit and other parameters associated with fitting a hit vector.
		 *  @param hit_vector contains hits - one per plane.
		 *  @param size is a size of hit_vector. May be only  3, 4 or 5.
		 *  @param output_values is an output parameter. It is a 10 elements array. output_values[0] contains chi square for x plane, output_values[1] contains chi square for y plane. Other values are used only in fitTrack().
		 */
		void ComputeParametersOfTrack(T1RecHitGlobal* hit_vector, int size, double* output_values);
		/**
		 *  Method that creates T1T2Track form given hit vector. This method is called only for selected tracks.
		 *  @param hit_vector contains hits - one per plane. T
		 *  @param size is a size of hit_vector. May be only  3, 4 or 5.
		 *  @return track with computed chi squared value, hemisphere, track parameters vector, covariance matrix,  chiSquaredX, chiSquaredY and detector.
		 */
		virtual T1T2Track fitTrack(T1RecHitGlobal* hit_vector, int size);

		/**
		 *  value set in python/RecoTotemT1T2/T1TrackProducer/T1TrackProducer_cfi.py
		 */
		double theChiProbThr;
		/**
		 *  value set in python/RecoTotemT1T2/T1TrackProducer/T1TrackProducer_cfi.py
		 */
		int theVerbosity;
};
#endif
