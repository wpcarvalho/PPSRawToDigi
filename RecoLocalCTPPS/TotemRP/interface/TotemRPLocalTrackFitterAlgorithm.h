/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
* 	Hubert Niewiadomski
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef RecoTotemRP_TotemRPLocalTrackFitterAlgorithm_TotemRPLocalTrackFitterAlgorithm_h
#define RecoTotemRP_TotemRPLocalTrackFitterAlgorithm_TotemRPLocalTrackFitterAlgorithm_h

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/CTPPSReco/interface/TotemRPLocalTrack.h"
#include "Geometry/VeryForwardGeometryBuilder/interface/TotemRPGeometry.h"
#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"
#include "Geometry/VeryForwardRPTopology/interface/RPTopology.h"

#include "TVector3.h"
#include "TVector2.h"

// TODO
#include <ext/hash_map>

//----------------------------------------------------------------------------------------------------

/**
 *\brief TODO
 **/
struct RPDetCoordinateAlgebraObjs
{
  //u*nom_pitch = (u0*pitch+u_cor)+(pitch/eff_pitch*rot_cor*v)*(x-x0)
  //             eff_u_middle_edge +    eff_v*(x-middle_of_edge_pos)
  //       (eff_u_middle_edge-eff_v*middle_of_edge_pos) + eff_v*x
  TVector3 centre_of_det_global_position_;
  double rec_u_0_;              ///< in mm, position of det. centre projected on readout direction
  TVector2 readout_direction_;  ///< non paralell projection and rot_cor included
  bool available_;              ///< if det should be included in the reconstruction
};

//----------------------------------------------------------------------------------------------------

/**
 *\brief TODO
 **/
class TotemRPLocalTrackFitterAlgorithm
{
  public:
    TotemRPLocalTrackFitterAlgorithm(const edm::ParameterSet &conf);

    /// performs the track fit, returns true if successful
    bool FitTrack(const vector<const TotemRPRecHit *> &hits, double z_0, const TotemRPGeometry &tot_geom, TotemRPLocalTrack &fitted_track);

    /// Resets the reconstruction-data cache.
    void Reset();

    // TODO: needed?
#if 0
    TVector2 ComputeXYPointInZDir(const TotemRPRecHit& hit_0, const TotemRPRecHit& hit_1,
      const TotemRPGeometry &tot_geom);

    TVector2 ComputeXYPointOfTheGivenLine(const TotemRPRecHit& hit_0, const TotemRPRecHit& hit_1,
      double tx, double ty, double z0, const TotemRPGeometry &tot_geom);
#endif

  private:
    /// A cache of reconstruction data. Must be reset every time the geometry chagnges.
    /// TODO: use unordered_map
    typedef __gnu_cxx::hash_map<unsigned int, RPDetCoordinateAlgebraObjs> DetReconstructionDataMap;
    DetReconstructionDataMap det_data_map_;

    RPTopology rp_topology_;

    /// Returns the reconstruction data for the chosen detector from the cache DetReconstructionDataMap.
    /// If it is not yet in the cache, calls PrepareReconstAlgebraData to make it.
    RPDetCoordinateAlgebraObjs *GetDetAlgebraData(unsigned int det_id, const TotemRPGeometry &tot_rp_geom);

    /// Build the reconstruction data.
    RPDetCoordinateAlgebraObjs PrepareReconstAlgebraData(unsigned int det_id, const TotemRPGeometry &tot_rp_geom);

    /// A matrix multiplication shorthand.
    void MultiplyByDiagonalInPlace(TMatrixD &mt, const TVectorD &diag);
    
    static TVector3 convert3vector(const CLHEP::Hep3Vector & v)
    {
      return TVector3(v.x(),v.y(),v.z()) ;
    }
};


#endif

