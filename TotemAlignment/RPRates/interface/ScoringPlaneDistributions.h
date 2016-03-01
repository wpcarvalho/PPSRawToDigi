/****************************************************************************
*
* This is a part of TOTEM offline software.
* Authors:
*  Jan Kašpar (jan.kaspar@gmail.com)
*
* $$RCSfile: ScoringPlaneDistributions.h,v $: $
* $Revision: 9977 $
* $Date: 2015-01-12 15:00:26 +0100 (Mon, 12 Jan 2015) $
*
****************************************************************************/

#include <vector>

class TH1D;
class TGraph;
class TotemRPGeometry;

/**
 *\brief 
 **/
class ScoringPlaneDistributions
{
  public:
    struct RP_profiles {
      unsigned long count;
      TH1D *xProfile, *yProfile;
      RP_profiles(signed int Id, const char *label);
      void Fill(double x, double y);                  ///< x and y in mm, in internal detector frame
    };

    typedef std::vector< std::vector< unsigned long > > overlap_statistics;

    ScoringPlaneDistributions() : geometry(NULL) {}
    void Init(const TotemRPGeometry *g, unsigned int N, double from, double to, signed int Id = 0);

    unsigned long eventCounter;

    void Fill(double x, double y);                    ///< x and y in mm
    void Write();                                     ///< writes data in gDirectory

  private:
    const TotemRPGeometry *geometry;
    std::vector<double> offsets;                      ///< offsets (edge to beam distances) in mm

    TGraph *hitDistribution;

    std::vector<RP_profiles> rpTop, rpHor, rpBot;
    overlap_statistics overlapTopHor, overlapBotHor;  ///< first index is for vertical RP while second for horizontal RP offset

    void WriteRPProfiles(const std::vector<RP_profiles> &);
    void WriteOverlapProfiles(const overlap_statistics &);
};

