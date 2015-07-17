/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
* $Id: TrackData.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#include "TotemCondFormats/BeamOpticsParamsObjects/interface/BeamOpticsParams.h"
#include "TotemRPValidation/ElasticReconstruction/interface/TrackData.h"

#include "TMath.h"

using namespace HepMC;


//----------------------------------------------------------------------------------------------------

TrackData::TrackData(const HepMC::FourVector &p, const BeamOpticsParams &par, bool correctForCrossingAngle)
{
  th_x = p.x() / p.z();
  if (correctForCrossingAngle)
    if (p.z() > 0)
      th_x -= par.GetCrossingAngleX();
    else
      th_x += par.GetCrossingAngleX();
  th_y = p.y() / p.z();
  th = sqrt(th_x*th_x + th_y*th_y);
  double P = p.rho();
  t_x = P*P * th_x*th_x;
  t_y = P*P * th_y*th_y;
  t = t_x + t_y;
  phi = TMath::ATan2(th_y, th_x);
}

//----------------------------------------------------------------------------------------------------

TrackData::TrackData(const RPRecoElasticEvent::fit_type &f, const BeamOpticsParams &par)
{
  th_x = f.th_x;
  th_y = f.th_y;
  th = sqrt(th_x*th_x + th_y*th_y);
  double E = par.GetBeamEnergy();
  t_x = E*E * th_x*th_x;
  t_y = E*E * th_y*th_y;
  t = t_x + t_y;
  phi = TMath::ATan2(th_y, th_x);
}


