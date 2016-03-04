/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*	Jan Kaspar (jan.kaspar@gmail.com) 
*    
* $$RCSfile: GeometryTestModule.h,v $: $
* $Revision: 10226 $
* $Date: 2015-03-18 17:44:24 +0100 (Wed, 18 Mar 2015) $
*
****************************************************************************/

#ifndef Geometry_TestModule_H
#define Geometry_TestModule_H

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"

/**
 * \ingroup TotemRPGeometry
 * \brief Testing module.
 *
 * See schema of \ref TotemRPGeometry "TOTEM RP geometry classes"
 **/
class GeometryTestModule : public edm::EDAnalyzer {
   public:
      explicit GeometryTestModule(const edm::ParameterSet&);
      ~GeometryTestModule();

   private:
      virtual void beginJob();
      virtual void analyze(const edm::Event&, const edm::EventSetup&);
      virtual void endJob() ;
};

#endif
