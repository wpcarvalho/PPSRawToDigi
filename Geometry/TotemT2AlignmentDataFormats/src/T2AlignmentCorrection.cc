/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	 
*    
* $$RCSfile: T2AlignmentCorrection.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include "Geometry/TotemT2AlignmentDataFormats/interface/T2AlignmentCorrection.h"

#include "Math/GenVector/RotationZ.h"

T2AlignmentCorrection::T2AlignmentCorrection(double sh_x, double sh_y, double rot_z) : 
  translation(sh_x, sh_y, 0.), rotation(ROOT::Math::RotationZ(rot_z))
{
}

//----------------------------------------------------------------------------------------------------

void T2AlignmentCorrection::Add(const T2AlignmentCorrection &a)
{
  rotation = a.rotation * rotation;
  translation = a.translation + translation;
}

//----------------------------------------------------------------------------------------------------

int T2AlignmentCorrection::ExtractStandardParameters(double &sh_x, double &sh_y, double &rot_z) const
{
  // extract shifts
  double limit = 1E-3;  // 1 um
  if (fabs(translation.z()) > limit)
    return 1;

  sh_x = translation.x();
  sh_y = translation.y();

  // extract z rotation
  double xx, xy, xz, yx, yy, yz, zx, zy, zz;
  rotation.GetComponents(xx, xy, xz, yx, yy, yz, zx, zy, zz);//Get all the 3x3 matrix components
  limit = 1E-6;
  if (fabs(zz - 1.) > limit || fabs(zx) > limit || fabs(zy) > limit || fabs(xz) > limit || fabs(yz) > limit)
    return 2;
  rot_z = atan2(xy, xx); // in -pi to +pi range

  return 0;
}

std::vector<double> T2AlignmentCorrection::getValues () const
{
  std::vector<double> result;
  double val1, val2, val3;
  this->ExtractStandardParameters (val1, val2, val3);
  result.push_back (val1);
  result.push_back (val2);
  result.push_back (val3);
  return result;
}



