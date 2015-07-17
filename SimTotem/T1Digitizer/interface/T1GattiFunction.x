#ifndef T1_GATTI_FUNCTION_H
#define T1_GATTI_FUNCTION_H


#include "Geometry/TotemGeometry/interface/T1ChamberSpecs.h"

/** \class T1GattiFunction
 *
 * Represent functional form of charge distribution over strips
 * in Endcap Muon CSC's.
 *
 *
 *
 * This is required in building RecHits from strips in MEClusterizer.
 * It was ported from FORTRAN. <BR>
 *
 *  Function: describes the cathode signal using                       <BR>
 *                the single-parameter Gatti formula:                  <BR>
 *                              1 - tanh(K_2 * lambda)**2              <BR>
 *     Gamma(lambda) = K_1 * -------------------------------           <BR>
 *                           1 + K_3 * tanh (K_2 *lambda)**2           <BR>
 *     lambda = x/h, h is anode cathode spacing                        <BR>
 *                                                                     <BR>
 *     K_2 = pi/2*(1 - 0.5*sqrt(K_3))                                  <BR>
 *                                                                     <BR>
 *              K_2*sqrt(K_3)                                          <BR>
 *      K_1 = -------------------                                      <BR>
 *            4 * atan(sqrt(K_3))                                      <BR>
 *                                                                     <BR>
 *  References  : E.Gatti, A.Longoni, NIM 163 (1979) 82-93.            <BR>
 *                                                                     <BR>
 *  For K_3, "It is used parametrization from Fig.2 from E.Mathieson   <BR>
 *            J.S.Gordon, "Cathode charge distributions in multi-      <BR>
 *            wire chambers", NIM 227 (1984) 277-282"                  <BR>
 *  (comment from GATTI3 in cmsim/src/mc_uty/.)                        <BR>
 */

//class T1ChamberSpecs;

class T1GattiFunction
{
public:
  T1GattiFunction() {};
  /// this routine calculates k1, k2, k3, and h, if needed
  void initChamberSpecs();

  ///  returns the fraction of charge on a strip centered
  ///  a distance of x away from the center of the shower,
  ///  at zero.  Note that the user is responsible for making
  ///  sure the constants have been initialized using the chamber specs.
  double binValue( double x, double stripWidth) const;

private:
  // geometry constants for the detector
  double k1, k2, k3, h;

  double norm, sqrtk3;

};

#endif
