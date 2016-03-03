/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Dominik Mierzejewski (dmierzej@cern.ch)
*    
* $$RCSfile: MeasuredGeometryRecord.h,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#ifndef RECORDS_MEASUREDGEOMETRYRECORD_H
#define RECORDS_MEASUREDGEOMETRYRECORD_H

#include "FWCore/Framework/interface/DependentRecordImplementation.h"
#include "Geometry/Records/interface/IdealGeometryRecord.h"

#include "boost/mpl/vector.hpp"

#include "TotemAlignment/RPRecords/interface/RPMeasuredAlignmentRecord.h"

/**
 * \ingroup TotemRPGeometry
 * \brief Event setup record containing the Measured (measured) geometry information.
 **/
class MeasuredGeometryRecord : public edm::eventsetup::DependentRecordImplementation
						   <MeasuredGeometryRecord, boost::mpl::vector<IdealGeometryRecord, RPMeasuredAlignmentRecord /*, ... */> >
{
};

#endif

