/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *    Jan Kaspar <jan.kaspar@gmail.com>
 *    Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RPStationMultiTrackFinderFitter_Cluster
#define RPStationMultiTrackFinderFitter_Cluster

#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h" // RPId
#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecognizedPatterns.h"

#include "TMatrixD.h"
#include "TVectorD.h"

#include <map>
#include <set>
#include <string>
#include <cstdio>

namespace RPStationMultiTrackFinderFitter {

using namespace std;

/**
 * \brief Keeps the information about a track
 *
 * Each track is described in terms of slope and intercept at z = z0.
 * Slope is composed of ax and ay, intercept of bx and by.
 *
 * Track's equation becomes:
 *
 * \f{eqnarray*}{
 * x(z) &=& ax * (z - z0) + bx \\
 * y(z) &=& ay * (z - z0) + by
 * \f}
 */
struct Track
{
    double ax, ///< in μrad
           ay; ///< in μrad
    double bx, ///< in mm
           by; ///< in mm
    double z0; ///< in mm

    Track(double _ax=0, double _ay=0, double _bx=0, double _by=0, double _z0=0) :
        ax(_ax), ay(_ay), bx(_bx), by(_by), z0(_z0)
    {
    }

    void Print(const string &prefix = "") const
    {
        printf("%sax = %+.1f urad, ay = %+.1f urad, bx = %+.3f mm, by = %+.3f mm, z0 = %.3f m\n",
                prefix.c_str(), ax*1E6, ay*1E6, bx, by, z0/1E3);
    }
};


/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Groups track candidates together forming a cluster with uncertainties
 *
 * \author Jan Kaspar <jan.kaspar@gmail.com>
 * \author Jakub Sawicki <jakub.kuba.sawicki@gmail.com>
 */
struct Cluster
{
    /// map: RP ID --> count
    map<RPId, unsigned int> rpIdCount;

    /// hits associated with the Cluster
    set<RPRecoHit> hits;

    /// row and column order: ax, bx, ay, by
    TMatrixD sum_W;   ///< inverse of uncertainty matrix for center (W = V^-1)
    TMatrixD sum_W_i; ///< uncertainty matrix for center
    TVectorD sum_Wv;  ///< sum of W * v = V^-1 * v for each contained track candidate
    TVectorD center;  ///< the track being the result of sum_W_i * sum_Wv

    Cluster() :
        sum_W(4, 4), sum_W_i(4, 4), sum_Wv(4), center(4) {}

    /// Returns a \ref Track corresponding to the \ref Cluster.
    /// No z0 is set though, if needed must be filled later.
    inline Track GetCenterTrack() const
    {
        Track t_center;
        t_center.ax = center(0);
        t_center.bx = center(1);
        t_center.ay = center(2);
        t_center.by = center(3);
        return t_center;
    }

    /**
     * \brief Calculates the distance between the \ref Cluster and a track
     *
     * It calculates the uncertainty matrix for `diff = v - center`:
     * `V_diff_i = (V_v + sum_W_)^-1`.
     * Based on it the distance is calculated by taking a sum over
     * `diff(i) * V_diff_i(i, j) * diff(j)`.
     *
     * \returns distance between the \ref Cluster and the track
     */
    double Distance(
            const TVectorD &v,            ///< [in] vector in format [ax bx ay by]
            const TMatrixD &V_v,          ///< [in] uncertainty matrix for v
            const TMatrixD &V_v_i) const; ///< [in] inverted uncertainty matrix for v

    /**
     * \brief Inserts a track candidate into \ref Cluster
     *
     * Increments the counters for rpId1 and rpId2.
     * Calculates `Wv = W * v` adding it to sum_Wv.
     * Adds W to sum_W.
     */
    void Add(RPId rpId1 /**< [in] */, RPId rpId2 /**< [in] */,
            const TVectorD &v,  ///< [in] vector in format [ax bx ay by]
            const TMatrixD &W); ///< [in] uncertainty matrix for v

    /// Copies the \ref Cluster into copy
    Cluster& Copy(Cluster &copy);
};

} // namespace RPStationMultiTrackFinderFitter

#endif // RPStationMultiTrackFinderFitter_Cluster
