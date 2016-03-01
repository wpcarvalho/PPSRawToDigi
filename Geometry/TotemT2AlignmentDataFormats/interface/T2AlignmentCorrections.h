/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: T2AlignmentCorrections.h,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#ifndef _T2AlignmentCorrections_h_
#define _T2AlignmentCorrections_h_

#include "Geometry/TotemT2AlignmentDataFormats/interface/T2AlignmentCorrection.h"

#include <map>

/**
 *\brief Container for T2 alignment corrections
 * map: raw detector ID --> alignment correction
 **/
class T2AlignmentCorrections : public std::map<unsigned int, T2AlignmentCorrection>
{
  public:
    void Add(unsigned int, const T2AlignmentCorrection&);

  /**
   * Inserts into T2AlignmentCorrections a T2AlignmentCorrection object with detector identifier.
   * \param identifier - detector id
   * \param values vector of arguments for T2AlignmentCorrection constructor
   * \throw cms::Exception if values vector size is not correct (!=3)
   */
  void insertValues (const std::string& identifier, const std::vector<double>& values);
};



#endif

