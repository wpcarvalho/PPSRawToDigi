/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	
*    
* $$RCSfile: T2AlignmentCorrections.cc,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include "Geometry/TotemT2AlignmentDataFormats/interface/T2AlignmentCorrections.h"

#include "FWCore/Utilities/interface/Exception.h"
#include "DataFormats/T2DetId/interface/T2DetId.h"
#include "Geometry/TotemGeometry/interface/T2GeometryUtil.h"

#include "FWCore/Utilities/interface/typelookup.h"


void T2AlignmentCorrections::Add(unsigned int id, const T2AlignmentCorrection &a)
{
  iterator it = find(id);
  if (it == end())
    insert(value_type(id, a));
  else it->second.Add(a);
}

void T2AlignmentCorrections::insertValues (const std::string& identifier, const std::vector<double>& values)
{
  unsigned int key = 0;
  sscanf (identifier.c_str (), "%u", &key);
  T2GeometryUtil conv;
  T2GeometryUtil::T2DetInfo planeinformation;
  uint32_t k2=uint32_t(key);
  planeinformation=conv.GetT2Info(k2);
  unsigned int t2rawid=(unsigned int) planeinformation.cmsswid;

  if (values.size () != 3) throw cms::Exception ("T2AlignmentCorrections::insertValues") << "Input vector size wrong" << std::endl;
  (*this)[ t2rawid /*TotemRPId::DecToRawId (key)*/ ] = T2AlignmentCorrection (values[0], values[1], values[2]);
}

TYPELOOKUP_DATA_REG(T2AlignmentCorrections);

