/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: AlignmentAlgorithm.cc 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#include "TotemAlignment/RPTrackBased/interface/AlignmentAlgorithm.h"
#include "TotemAlignment/RPTrackBased/interface/AlignmentTask.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

//----------------------------------------------------------------------------------------------------

AlignmentAlgorithm::AlignmentAlgorithm(const edm::ParameterSet& ps, AlignmentTask *_t) :
  verbosity(ps.getUntrackedParameter<unsigned int>("verbosity", 0)),
  task(_t),
  singularLimit(ps.getParameter<double>("singularLimit")),
  useExternalFitter(ps.getParameter<bool>("useExternalFitter"))
{
}

