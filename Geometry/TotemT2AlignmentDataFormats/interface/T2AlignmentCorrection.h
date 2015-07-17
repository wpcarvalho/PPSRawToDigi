/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: T2AlignmentCorrection.h,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#ifndef _T2AlignmentCorrection_h_
#define _T2AlignmentCorrection_h_

#include "DetectorDescription/Base/interface/DDRotationMatrix.h"
#include "DetectorDescription/Base/interface/DDTranslation.h"
#include <vector>

/**
 *\brief Alignment correction of a T2 silicon detector
 * If a detector's original rotation is R and shift is S, then, when a correction is applied, 
 *    the full rotation is: rotation * R and
 *    the full shift: is shift + S
 * see DetGeomDesc::ApplyAlignment .
 *
 * In other words, the transformation is in the global coordinate system x, y, z (or s). The
 * transformation applies TO DETECTOR POSITION and NOT TO HIT POSITION.
 **/
class T2AlignmentCorrection {
  protected:
    DDTranslation translation;        ///< shift in mm; currently implemented as ROOT::Math::DisplacementVector3D
    DDRotationMatrix rotation;        ///< rotation "in rad"; currently implemented as ROOT::Math::Rotation3D

  public:
    T2AlignmentCorrection() {}
    //T2AlignmentCorrection(const DDTranslation &_t, const DDRotationMatrix &_r) : translation(_t), rotation(_r) {}

    /// constructor for simplified case with shifts in x and y and rotation around z only
    /// shifts in mm, rotation in rad
    T2AlignmentCorrection(double sh_x, double sh_y, double rot_z);

    /// merges (cumulates) alignements
    void Add(const T2AlignmentCorrection&);

    const DDTranslation& Translation() const
      { return translation; }

    const DDRotationMatrix& Rotation() const
      { return rotation; }

    /// extracts the standard parameters from general translation and rotation
    /// returns non-zero value if the alignment correction is more general than the standard
    int ExtractStandardParameters(double &sh_x, double &sh_y, double &rot_z) const;

    /**
   * Return a vector<double> describing T2 Alignment Correction.
   * Should return the same values as ExtractStandardParameters
   */
  std::vector<double> getValues () const;


};

#endif

