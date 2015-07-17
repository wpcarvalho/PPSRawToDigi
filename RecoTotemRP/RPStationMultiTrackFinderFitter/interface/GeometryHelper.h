/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#ifndef RPStationMultiTrackFinderFitter_GeometryHelper
#define RPStationMultiTrackFinderFitter_GeometryHelper

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/TotemRPDataTypes/interface/RPTypes.h"
#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"

#include <map>

using namespace std;

namespace RPStationMultiTrackFinderFitter {

struct RPInfo;
struct DetInfo;

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Helper for the geometry
 *
 * Provides convenient view of the geometry from the perspective of track finding
 * and fitting.
 *
 * For each RP, two directions are available, U and V. The cut edge's position is
 * represented as a position vector and normal vector facing outwards.
 *
 * The returned RP information is cached for multiple use.
 */
class GeometryHelper
{
  public:
    explicit GeometryHelper(const edm::ParameterSet &ps);
    virtual ~GeometryHelper();

    /// Updates the geometry information, resets the cache of RPInfo objects
    inline void setGeometry(const TotemRPGeometry *geometry)
    {
        this->geometry = geometry;
        rpInfos.clear();
        detInfos.clear();
    }

    /// \returns if the point (in global coordinates) lies within the sensor
    virtual bool isPointOutside(RPId rpID, CLHEP::Hep3Vector point);

    /// \returns if the point (in global coordinates) lies outside the sensor's edge
    virtual bool isPointOutsideEdge(RPId rpID, CLHEP::Hep3Vector point);

    /// RPInfo is created if not have existed or found in rpInfos otherwise
    /// \returns the RPInfo for an RP with ID = RPId
    RPInfo getRPInfo(RPId RPId);

    /// DetInfo is created if not have existed or found in detInfos otherwise
    /// \returns the DetInfo for a sensor with ID = detId
    DetInfo getDetInfo(RPDetId detId);

  protected:
    /// Creates RPInfo object.
    /// The values are averaged information about the angles and position
    /// of the center of the RP. Edge position is also determined.
    RPInfo buildRPInfo(RPId rpId);

    /// Created DetInfo object.
    /// Copies information about the RP from geometry to the object.
    DetInfo buildDetInfo(RPDetId detId);

  private:
    const unsigned int verbosity;
    const TotemRPGeometry *geometry;

    /// error margin for the edge check
    const double outside_edge_th;

    /// map: RPId -> RPInfo
    map<RPId, RPInfo> rpInfos;

    /// map: detId -> DetInfo
    map<RPDetId, DetInfo> detInfos;
};

//----------------------------------------------------------------------------------------------------

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Contains RP information
 */
struct RPInfo
{
    RPId RPid;                        ///< ID of the RP described
    double cx, cy, cz;                ///< position of the RP's center, TODO: units
    double u_dx, u_dy, v_dx, v_dy;    ///< projection of direction vectors of U
                                      ///< and V directions onto X and Y
    CLHEP::Hep3Vector edge_p, edge_n; ///< global position and normal (pointing
                                      ///< outwards the sensor) of the edge
};

//----------------------------------------------------------------------------------------------------

/**
 * \ingroup RPStationMultiTrackFinderFitter
 * \brief Contains single sensor information
 */
struct DetInfo
{
    RPDetId detId;      ///< ID of the sensor described
    double cx, cy, cz;  ///< position of the RP's center, TODO: units
    double dx, dy;      ///< projection of the readout direction of the sensor onto X and Y
};

} // namespace RPStationMultiTrackFinderFitter

#endif // #ifndef RPStationMultiTrackFinderFitter_GeometryHelper
