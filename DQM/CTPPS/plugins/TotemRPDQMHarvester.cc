/****************************************************************************
*
* This is a part of TotemDQM and TOTEM offline software.
* Authors:
*   Jan Kašpar (jan.kaspar@gmail.com)
*
****************************************************************************/

#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DQMServices/Core/interface/DQMEDHarvester.h"
#include "DQMServices/Core/interface/DQMStore.h"
#include "DQMServices/Core/interface/MonitorElement.h"

#include "DataFormats/TotemRPDetId/interface/TotemRPDetId.h"

//----------------------------------------------------------------------------------------------------
 
class TotemRPDQMHarvester: public DQMEDHarvester
{
  public:
    TotemRPDQMHarvester(const edm::ParameterSet& ps);
    virtual ~TotemRPDQMHarvester();
  
  protected:
    void dqmEndJob(DQMStore::IBooker &, DQMStore::IGetter &) override;

  private:
    void MakePlaneEfficiencyHistograms(unsigned int plId, DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter);
};

//----------------------------------------------------------------------------------------------------
//----------------------------------------------------------------------------------------------------

using namespace std;
using namespace edm;

//----------------------------------------------------------------------------------------------------

TotemRPDQMHarvester::TotemRPDQMHarvester(const edm::ParameterSet& ps)
{
}

//----------------------------------------------------------------------------------------------------

TotemRPDQMHarvester::~TotemRPDQMHarvester()
{
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMHarvester::MakePlaneEfficiencyHistograms(unsigned int id, DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter)
{
  printf(">> TotemRPDQMHarvester::MakePlaneEfficiencyHistograms(%u)\n", id);

  string path = TotemRPDetId::planeName(id, TotemRPDetId::nPath);
  path.replace(0, 2, "TrackingStrip");
  ibooker.setCurrentFolder(string("CTPPS/") + path);


  printf("    num: %s\n", ("CTPPS/" + path + "/efficiency_num").c_str());

  MonitorElement *efficiency_num = igetter.get("CTPPS/" + path + "/efficiency_num");
  MonitorElement *efficiency_den = igetter.get("CTPPS/" + path + "/efficiency_den");

  printf("    num = %p, den = %p\n", efficiency_num, efficiency_den);

  if (!efficiency_num || !efficiency_den)
    return;
  
  string title = TotemRPDetId::planeName(id, TotemRPDetId::nFull);
  
  MonitorElement *efficiency = ibooker.book1D("efficiency", title+";track position   (mm)", 30, -15., 0.);

  for (signed int bi = 1; bi <= efficiency->getNbinsX(); bi++)
  {
    double num = efficiency_num->getBinContent(bi);
    double den = efficiency_den->getBinContent(bi);

    if (den > 0)
    {
      double p = num / den;
      double p_unc = sqrt(den * p * (1. - p));
      efficiency->setBinContent(bi, p);
      efficiency->setBinError(bi, p);
    } else {
      efficiency->setBinContent(bi, 0.);
    }
  }
}

//----------------------------------------------------------------------------------------------------

void TotemRPDQMHarvester::dqmEndJob(DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter)
{
  // loop over arms
  for (unsigned int arm = 0; arm < 2; arm++)
  {
    // loop over stations
    for (unsigned int st = 0; st < 3; st += 2)
    {
      unsigned int stId = 10*arm + st;
      
      // loop over RPs
      for (unsigned int rp = 0; rp < 6; ++rp)
      {
        unsigned int rpId = 10*stId + rp;

        // loop over planes
        for (unsigned int pl = 0; pl < 10; ++pl)
        {
          unsigned int plId = 10*rpId + pl;

          MakePlaneEfficiencyHistograms(plId, ibooker, igetter);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRPDQMHarvester);
