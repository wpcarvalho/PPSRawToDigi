/****************************************************************************
*
* This is a part of the TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: TrackData.h 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/


#ifndef _TrackData_h_
#define _TrackData_h_

class BeamOpticsParams;


#include "RecoTotemRP/RPRecoDataFormats/interface/RPRecoElasticEvent.h"

/**
 *\brief precalculated track data: theta, t, phi, etc.
 **/
struct TrackData
{
  double t, t_x, t_y;  ///< absolute value of t
  double th, th_x, th_y;
  double phi;
  TrackData() : t(0.), t_x(0.), t_y(0.), th(0.), th_x(0.), th_y(0.), phi(0.) {}
  TrackData(const HepMC::FourVector &, const BeamOpticsParams &par, bool correctForCrossingAngle = true);
  TrackData(const RPRecoElasticEvent::fit_type &, const BeamOpticsParams &par);
};

#endif 

