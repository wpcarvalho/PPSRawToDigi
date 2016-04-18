/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#ifndef DataFormats_CTPPSReco_TotemRPUVPattern
#define DataFormats_CTPPSReco_TotemRPUVPattern

#include "DataFormats/CTPPSReco/interface/TotemRPRecHit.h"

#include <vector>

/**
 *\brief A linear pattern in U or V projection.
 * The intercept b is taken at the middle of a RP:
 *     (geometry->GetRPDevice(RPId)->translation().z())
 * The global coordinate system is used (wrt. the beam). This is the same convention
 * as for the 1-RP track fits.
 **/
class TotemRPUVPattern
{
  public:
    enum ProjectionType { projInvalid, projU, projV };

    TotemRPUVPattern() : projection(projInvalid), a(0.), b(0.), w(0.), fittable(false)
    {
    }

    ProjectionType getProjection() const { return projection; }
    void setProjection(ProjectionType p_) { projection = p_; }

    double getA() const { return a; }
    void setA(double a_) { a = a_; }

    double getB() const { return b; }
    void setB(double b_) { b = b_; }

    double getW() const { return w; }
    void setW(double w_) { w = w_; }

    bool getFittable() const { return fittable; }
    void setFittable(bool f_) { fittable = f_; }

    void addHit(const TotemRPRecHit &hit) { hits.push_back(hit); }

    const std::vector<TotemRPRecHit>& getHits() const { return hits; }

    friend bool operator< (const TotemRPUVPattern &l, const TotemRPUVPattern &r);

  private:
    ProjectionType projection;        ///< projection
    double a;                         ///< slope in rad
    double b;                         ///< intercept in mm
    double w;                         ///< weight
    bool fittable;                    ///< whether this pattern is worth including in track fits
    std::vector<TotemRPRecHit> hits;  ///< hits associated with the pattern
};

//----------------------------------------------------------------------------------------------------

extern bool operator< (const TotemRPUVPattern &l, const TotemRPUVPattern &r);

#endif
