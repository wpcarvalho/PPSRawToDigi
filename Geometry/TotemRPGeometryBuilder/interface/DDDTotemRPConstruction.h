/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: DDDTotemRPConstruction.h,v $: $
* $Revision: 10226 $
* $Date: 2015-03-18 17:44:24 +0100 (Wed, 18 Mar 2015) $
*
****************************************************************************/

#ifndef Geometry_DDDTotemRPContruction_H
#define Geometry_DDDTotemRPContruction_H

#include "Geometry/TotemRPGeometryBuilder/interface/DetGeomDesc.h"
#include "Geometry/TotemRPGeometryBuilder/interface/DDDTotemRPCommon.h"

class DDCompactView;
class DDFilteredView;


/**
 * \ingroup TotemRPGeometry
 * \brief Builds scructure of DetGeomDesc objects out of DDCompactView (resp. DDFilteredView).
 *
 * See schema of \ref TotemRPGeometry "TOTEM RP geometry classes"
 *
 * It adds detector IDs (via class TotemRPDetId).
 * intended to be called from: modul TotemRPDetGeomDescESModule.
 **/

class DDDTotemRPContruction {
	public:
		DDDTotemRPContruction();
		const DetGeomDesc* construct(const DDCompactView* cpv);

	protected:
		void buildDetGeomDesc(DDFilteredView *fv, DetGeomDesc *gd);
};

#endif
