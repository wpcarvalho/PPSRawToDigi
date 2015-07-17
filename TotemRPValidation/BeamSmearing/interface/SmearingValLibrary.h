/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors: 
*  Jan Kašpar (jan.kaspar@gmail.com) 
*    
* $Id: SmearingValLibrary.h 9977 2015-01-12 14:00:26Z tsodzawi $
* $Revision: 9977 $
* $Date: 2015-01-12 16:00:26 +0200 (pon, 12 sty 2015) $
*
****************************************************************************/

#ifndef _TotemRPValidationElasticReconstructionSmearingValLibrary_H_
#define _TotemRPValidationElasticReconstructionSmearingValLibrary_H_

#include <string>
#include <map>

namespace edm {
  class ParameterSet;
  class EventSetup;
  class Event;
}

class TH1D;

/**
 *\brief Code library for validation of SmearingGenerator package.
 * Actually called from SmearingValidation plugin.
**/
class SmearingValLibrary
{
 public:
  struct SmearInfo
  {
    TH1D *xi, *th_x, *th_y, *t;
    SmearInfo();
  };

  enum region {rForw, rCent, rBack};
  
  SmearingValLibrary(const edm::ParameterSet&);
  ~SmearingValLibrary();

  void initialize();
  void analyze(const edm::Event&, const edm::EventSetup&);
  void finalize();
  void writeHistogramsToFile();

 private:
  static const char* regionNames[];

  unsigned char verbosity;
  std::string generatorLabel;
  std::string originalLabel;
  std::string outputFile;
  double thetaLim;                    ///< angular limit to distinguish between forward, central and backward particles

  TH1D *vtx_x, *vtx_y, *vtx_z;
  std::map<region, SmearInfo> hists;
};

#endif 

