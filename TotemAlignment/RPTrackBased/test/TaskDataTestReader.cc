/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: analyzeStatisticalResults.cc 3294 2010-09-07 15:15:32Z jkaspar $
* $Revision: 3294 $
* $Date: 2010-09-07 17:15:32 +0200 (Tue, 07 Sep 2010) $
*
****************************************************************************/

#include "TotemAlignment/RPTrackBased/interface/AlignmentGeometry.h"

#include "TFile.h"

using namespace std;


//----------------------------------------------------------------------------------------------------

int main()
{
  printf("bla\n");

  TFile *f = new TFile("task_data.root");

  printf("f = %p\n", (void*) f);

  AlignmentGeometry *g = (AlignmentGeometry *) f->Get("geometry");

  printf("g = %p\n", (void*) g);

  g->Print();

  return 0;
}

