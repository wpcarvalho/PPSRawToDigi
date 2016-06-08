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
    void MakePlaneEfficiencyHistograms(unsigned int plId, DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter,
      MonitorElement* &rp_efficiency);
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

void TotemRPDQMHarvester::MakePlaneEfficiencyHistograms(unsigned int id, DQMStore::IBooker& ibooker, DQMStore::IGetter& igetter,
  MonitorElement* &rp_efficiency)
{
  printf(">> TotemRPDQMHarvester::MakePlaneEfficiencyHistograms(%u)\n", id);

  // get source histograms
  string path = TotemRPDetId::planeName(id, TotemRPDetId::nPath);
  path.replace(0, 2, "TrackingStrip");

  MonitorElement *efficiency_num = igetter.get("CTPPS/" + path + "/efficiency num");
  MonitorElement *efficiency_den = igetter.get("CTPPS/" + path + "/efficiency den");

  if (!efficiency_num || !efficiency_den)
    return;

  // book new histogram
  ibooker.setCurrentFolder(string("CTPPS/") + path);

  string title = TotemRPDetId::planeName(id, TotemRPDetId::nFull);
  
  TAxis *axis = efficiency_den->getTH1()->GetXaxis();
  MonitorElement *efficiency = ibooker.book1D("efficiency", title+";track position   (mm)", axis->GetNbins(), axis->GetXmin(), axis->GetXmax());

  // book new RP histogram (if not yet done)
  if (rp_efficiency == NULL)
  {
    printf("    new RP hist, %i\n", id);

    path = TotemRPDetId::rpName(id/10, TotemRPDetId::nPath);
    path.replace(0, 2, "TrackingStrip");
    title = TotemRPDetId::rpName(id/10, TotemRPDetId::nFull);
    ibooker.setCurrentFolder(string("CTPPS/") + path);
    rp_efficiency = ibooker.book2D("plane efficiency", title+";plane;track position   (mm)",
      10, -0.5, 9.5, axis->GetNbins(), axis->GetXmin(), axis->GetXmax());
  }

  // calculate and fill efficiencies
  for (signed int bi = 1; bi <= efficiency->getNbinsX(); bi++)
  {
    double num = efficiency_num->getBinContent(bi);
    double den = efficiency_den->getBinContent(bi);

    if (den > 0)
    {
      double p = num / den;
      double p_unc = sqrt(p * (1. - p) / den);
      efficiency->setBinContent(bi, p);
      efficiency->setBinError(bi, p_unc);

      int pl_bi = (id%10) + 1;
      rp_efficiency->setBinContent(pl_bi, bi, p);

      printf("        %i, %i, %.2f\n", pl_bi, bi, p);
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

        MonitorElement *rp_efficiency = NULL;

        // loop over planes
        for (unsigned int pl = 0; pl < 10; ++pl)
        {
          unsigned int plId = 10*rpId + pl;

          MakePlaneEfficiencyHistograms(plId, ibooker, igetter, rp_efficiency);
        }
      }
    }
  }
}

//----------------------------------------------------------------------------------------------------

DEFINE_FWK_MODULE(TotemRPDQMHarvester);
