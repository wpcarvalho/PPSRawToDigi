/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: RealGeometryRecord.h,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#ifndef RECORDS_REALGEOMETRYRECORD_H
#define RECORDS_REALGEOMETRYRECORD_H

#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "Geometry/TotemRecords/interface/MeasuredGeometryRecord.h"

#include "boost/mpl/vector.hpp"

#include "TotemAlignment/RPRecords/interface/RPRealAlignmentRecord.h"

/**
 * \ingroup TotemRPGeometry
 * \brief Event setup record containing the real (actual) geometry information.
 **/
class RealGeometryRecord : public edm::eventsetup::DependentRecordImplementation
						   <RealGeometryRecord, boost::mpl::vector<MeasuredGeometryRecord, RPRealAlignmentRecord /*, ... */> >
{
};

#endif

