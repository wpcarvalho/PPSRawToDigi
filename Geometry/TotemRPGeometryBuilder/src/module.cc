/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: module.cc,v $: $
* $Revision: 10226 $
* $Date: 2015-03-18 18:44:24 +0200 (śro, 18 mar 2015) $
*
****************************************************************************/

#include "FWCore/Utilities/interface/typelookup.h"

#include "Geometry/TotemRPGeometryBuilder/interface/TotemRPGeometry.h"
#include "Geometry/TotemRPGeometryBuilder/interface/DetGeomDesc.h"

#include "DataFormats/DetId/interface/DetId.h"

TYPELOOKUP_DATA_REG(TotemRPGeometry);
TYPELOOKUP_DATA_REG(DetGeomDesc);
