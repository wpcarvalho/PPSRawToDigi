/****************************************************************************
 *
 * This is a part of TOTEM offline software.
 * Authors:
 *   Jakub Sawicki (jakub.kuba.sawicki@gmail.com)
 *
 ****************************************************************************/

#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/GeometryHelper.h"
#include "Geometry/TotemRPDetTopology/interface/RPTopology.h"

#include <set>

namespace RPStationMultiTrackFinderFitter {

GeometryHelper::GeometryHelper(const edm::ParameterSet &ps) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity")),
    outside_edge_th(ps.getParameter<double>("outside_edge_th"))
{
}

//----------------------------------------------------------------------------------------------------

GeometryHelper::~GeometryHelper()
{
}

//----------------------------------------------------------------------------------------------------

bool GeometryHelper::isPointOutside(RPId rpID, CLHEP::Hep3Vector point)
{
    RPInfo rpinfo = getRPInfo(rpID);

    CLHEP::Hep3Vector point_rp = point - CLHEP::Hep3Vector(rpinfo.cx, rpinfo.cy, rpinfo.cz);

    double u = point_rp.x() * rpinfo.u_dx + point_rp.y() * rpinfo.u_dy,
           v = point_rp.x() * rpinfo.v_dx + point_rp.y() * rpinfo.v_dy;

    return !RPTopology::IsHit(u, v);
}

//----------------------------------------------------------------------------------------------------

bool GeometryHelper::isPointOutsideEdge(RPId rpID, CLHEP::Hep3Vector point)
{
    RPInfo rpinfo = getRPInfo(rpID);
    CLHEP::Hep3Vector d = point - rpinfo.edge_p;

    // a dot product of a normal of the edge (facing outwards) and the position relative
    // to the detector's edge
    // if positive then point lies outside the edge
    return (rpinfo.edge_n.dot(d) > outside_edge_th);
}

//----------------------------------------------------------------------------------------------------

RPInfo GeometryHelper::getRPInfo(RPId rpId)
{
    if (rpInfos.find(rpId) == rpInfos.end())
    {
        rpInfos.insert(pair<RPId, RPInfo>(rpId, buildRPInfo(rpId)));

        // TODO: return result immediately
    }

    return rpInfos.find(rpId)->second;
}

//----------------------------------------------------------------------------------------------------

RPInfo GeometryHelper::buildRPInfo(RPId rpId)
{
    unsigned int Nu = 0, Nv = 0;
    double Su = 0., Sv = 0.;
    CLHEP::Hep3Vector edge_p, edge_n, trans; // defaults to (0,0,0)

    set<unsigned int> DetIds = geometry->DetsInRP(rpId);

    for (set<unsigned int>::iterator det = DetIds.begin(); det != DetIds.end(); ++det)
    {
        unsigned int rawId = TotRPDetId::DecToRawId(*det);

        CLHEP::Hep3Vector d = geometry->LocalToGlobalDirection(rawId, CLHEP::Hep3Vector(0., 1., 0.));
        double phi = atan2(d.y(), d.x());
        if (TotRPDetId::IsStripsCoordinateUDirection(*det))
        {
            Nu++;
            Su += phi;
        } else {
            Nv++;
            Sv += phi;
        }

        edge_p += geometry->GetDetEdgePosition(rawId);
        edge_n += geometry->GetDetEdgeNormalVector(rawId);

        trans += geometry->GetDetTranslation(rawId);
    }

    // the "average" values of dx and dy
    double md_u_x = cos(Su / Nu), md_u_y = sin(Su / Nu);
    double md_v_x = cos(Sv / Nv), md_v_y = sin(Sv / Nv);

    // and the edge vectors
    edge_p /= Nu + Nv;
    edge_n /= Nu + Nv;
    trans /= Nu + Nv;

    RPInfo rpinfo;
    rpinfo.RPid   = rpId;
    rpinfo.cx     = trans.x();
    rpinfo.cy     = trans.y();
    rpinfo.cz     = trans.z();
    rpinfo.u_dx   = md_u_x;
    rpinfo.u_dy   = md_u_y;
    rpinfo.v_dx   = md_v_x;
    rpinfo.v_dy   = md_v_y;
    rpinfo.edge_p = edge_p;
    rpinfo.edge_n = edge_n;

    return rpinfo;
}

//----------------------------------------------------------------------------------------------------

DetInfo GeometryHelper::getDetInfo(RPDetId detID)
{
    if (detInfos.find(detID) == detInfos.end())
    {
        detInfos.insert(pair<RPDetId, DetInfo>(detID, buildDetInfo(detID)));
        // TODO: return immediately
    }

    return detInfos.find(detID)->second;
}

//----------------------------------------------------------------------------------------------------

DetInfo GeometryHelper::buildDetInfo(RPDetId detID)
{
    CLHEP::Hep3Vector d = geometry->LocalToGlobalDirection(detID, CLHEP::Hep3Vector(0., 1., 0.));
    DDTranslation c = geometry->GetDetector(detID)->translation();

    DetInfo newInfo;
    newInfo.detId = detID;
    newInfo.dx    = d.x();
    newInfo.dy    = d.y();
    newInfo.cx    = c.x();
    newInfo.cy    = c.y();
    newInfo.cz    = c.z();

    return newInfo;
}


} // namespace RPStationMultiTrackFinderFitter
