#include "FWCore/Framework/interface/EDAnalyzer.h"

#include <TGraph.h>

#include <map>
#include <vector>

class InelasticPlots
{
	public:
		explicit InelasticPlots();
		~InelasticPlots();

    		void prepareResPlots(std::string);
		std::vector<TCanvas*> getResolutionPlots();
	private:
		unsigned char verbosity;
		std::vector<TCanvas*> resolutionPlots;

		std::string histogramFile;
};
