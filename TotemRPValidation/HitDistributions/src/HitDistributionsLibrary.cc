#include "TotemRPValidation/HitDistributions/interface/HitDistributionsLibrary.h"

//----------------------------------------------------------------------------------------------------

HitDistributionsLibrary::HitDistributionsLibrary(const edm::ParameterSet& conf)
{
  gSystem->Load("libHistPainter.so");
  using namespace std;
  verbosity = conf.getUntrackedParameter<unsigned int>("verbosity", 0);
  validTracksOnly = conf.getUntrackedParameter<unsigned int>("validTracksOnly", 1);
  outputFile = conf.getParameter<std::string>("outputFile");
  vector< unsigned int > rpList = conf.getParameter< vector<unsigned int> >("rpList");
  rpFittedTrackCollectionLabel = conf.getParameter<edm::InputTag>("RPFittedTrackCollectionLabel");
  // initialize vector of TGraph's and vector of unitIDs
  for (unsigned int i = 0; i < rpList.size(); i++) {
    hitDists[rpList[i]] = new TGraph();
    // set name for TGraph
    char buf[20];
    sprintf(buf, "det_%i", rpList[i]);
    hitDists[rpList[i]]->SetName(buf);
    hitDists[rpList[i]]->SetTitle(buf);

    // calculate unitID based on detector ID
    // if there is at least on detector per unit
    //   present in rpList, then unitID will be added
    //   to unitIDs vector
    unsigned int unitID;
    if( ((rpList[i]%10) == 0) || ((rpList[i]%10) == 3))
      unitID = rpList[i];
    if( ((rpList[i]%10) == 1) || ((rpList[i]%10) == 4))
      unitID = rpList[i]-1;
    if( ((rpList[i]%10) == 2) || ((rpList[i]%10) == 5))
      unitID = rpList[i]-2;
    if( unitIDs.count(unitID) == 0){
      unitIDs.insert(unitID);
    }

    // correspondence between unitIDs and detectorIDs :
    // unit 000 : detectors 000,001,002
    // unit 003 : detectors 003,004,005
    // unit 020 : detectors 020,021,022
    // unit 023 : detectors 023,024,025
    // unit 100 : detectors 100,101,102
    // unit 103 : detectors 103,104,105
    // unit 120 : detectors 120,121,122
    // unit 123 : detectors 123,124,125

  }

}

//----------------------------------------------------------------------------------------------------

HitDistributionsLibrary::~HitDistributionsLibrary()
{
}

//----------------------------------------------------------------------------------------------------

void HitDistributionsLibrary::analyze(const edm::Event& event, const edm::EventSetup& eSetup)
{
  using namespace std;
  using namespace edm;

  // ------------------------------ GET EVENTS
  //for the moment - ideal geometry
  edm::ESHandle<TotemRPGeometry> Totem_RP_geometry;
  eSetup.get<RealGeometryRecord>().get(Totem_RP_geometry);

  edm::Handle< RPFittedTrackCollection > fitTrCol;
  event.getByLabel(rpFittedTrackCollectionLabel, fitTrCol);

  if(verbosity)
    std::cout << " Total fitted tracks: " << (fitTrCol)->size() << std::endl;

  if(verbosity)
    std::cout << " Number of detectors: " << Totem_RP_geometry->NumberOfDetsIncluded() << std::endl;

  for (std::map< unsigned int, TGraph *>::const_iterator it = hitDists.begin(); it != hitDists.end(); ++it) {
    // create id of first detector in RP station
    unsigned int rawid = TotRPDetId::DecToRawId(10*(it->first));
    TotRPDetId rpid(rawid);

    // add translation of detector to list
    if (rpTranslations.find(it->first) == rpTranslations.end()){
      CLHEP::Hep3Vector trans = Totem_RP_geometry->GetDetTranslation(rpid);
      rpTranslations[it->first] = trans;
    }

    // add rotation of detector to list
    if (rpRotations.find(it->first) == rpRotations.end()){
      double xx, xy, xz, yx, yy, yz, zx, zy, zz;
      (Totem_RP_geometry->GetDetector(rpid)->rotation()).GetComponents(xx, xy, xz, yx, yy, yz, zx, zy, zz);
      CLHEP::HepRep3x3 rot_mat( xx, xy, xz, yx, yy, yz, zx, zy, zz);
      CLHEP::HepRotation rot(rot_mat);
      rpRotations[it->first] = rot;
    }
  }


  // ------------------------------ ONE RP FIT
  for (std::map<RPId, RPFittedTrack>::const_iterator it = fitTrCol->begin(); it != fitTrCol->end(); ++it) {

    // check if track is valid
    if( validTracksOnly )
      if (!it->second.IsValid()) continue;

    // check if graph corresponding to RPId exists in hitDists map
    if (hitDists.find(it->first) != hitDists.end()){

      // add point to corresponding TGraph
      hitDists[it->first]->SetPoint(hitDists[it->first]->GetN(), it->second.X0(), it->second.Y0());

      if(verbosity)
        std::cout << "Hit in detector " << it->first << "  : ( " << it->second.X0() << " , " << it->second.Y0() << " ) , IsValid()=" << it->second.IsValid() << std::endl;

    } else {
      if(verbosity)
        std::cout << "RPId " << it->first << " not specified in cfg file" << std::endl;
    }

  }

}

//----------------------------------------------------------------------------------------------------

void HitDistributionsLibrary::writeToFile()
{
  using namespace std;
  using namespace edm;

  // =================================== save
  TFile *of = TFile::Open(outputFile.c_str(), "recreate");
  if(!of || !of->IsWritable())
  {
    std::cout<<"Output file not opened correctly!!"<<std::endl;
  }

  std::vector<TCanvas*>  histGraphs = getHistograms();
  for(unsigned int i=0; i<histGraphs.size(); ++i){
    histGraphs[i]->Write("");
  }

  of->Close();
}

//----------------------------------------------------------------------------------------------------

std::vector<TCanvas*>  HitDistributionsLibrary::getHistograms()
{
  return hitsGraphs;
}

//----------------------------------------------------------------------------------------------------

void HitDistributionsLibrary::ExportAllHistograms()
{

  // loop over units
  for (std::set< unsigned int>::iterator it = unitIDs.begin(); it != unitIDs.end(); ++it) {
    unsigned int unitID = *it;
    char buf[20];
    sprintf(buf, "unit_%i", unitID);

    // create one TCanvas per unitID
    TCanvas * c1 = new TCanvas(buf,buf,200,10,700,500);
    c1->cd();
    TH2D *frame = new TH2D(buf, ";x   (mm);y   (mm)", 100, -70., +70., 100, -70., +70.);
    frame->Draw();

    setupCanvas();
    TLegend *l = new TLegend(0.1, 0.1, 0.2, 0.2);

    int color = 0;
    int totalNumberOfHits = hitDists.count( unitID );
    totalNumberOfHits += hitDists.count( unitID+1 );
    totalNumberOfHits += hitDists.count( unitID+2 );

    // draw RP envelopes for all 3 RPs
    drawRPPlaneEnvelopes(unitID);
    drawRPPlaneEnvelopes(unitID+1);
    drawRPPlaneEnvelopes(unitID+2);


    // draw hits in first RP (if present)
    if( hitDists.count( unitID ) > 0 ){
      TGraph *g = hitDists[unitID];
      if (g->GetN() != 0){
        g->SetMarkerColor(++color);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.3);
        g->Draw("P");
        l->AddEntry(g, g->GetName(), "p");
      }
    }

    // draw hits in second RP (if present)
    if( hitDists.count( unitID + 1 ) > 0 ){
      TGraph *g = hitDists[unitID+1];
      if (g->GetN() != 0){
        g->SetMarkerColor(++color);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.3);
        g->Draw("P");
        l->AddEntry(g, g->GetName(), "p");
      }
    }

    // draw hits in third RP (if present)
    if( hitDists.count( unitID + 2) > 0 ){
      TGraph *g = hitDists[unitID+2];
      if (g->GetN() != 0){
        g->SetMarkerColor(++color);
        g->SetMarkerStyle(20);
        g->SetMarkerSize(0.3);
        g->Draw("P");
        l->AddEntry(g, g->GetName(), "p");
      }
    }

    if( totalNumberOfHits > 0 ){
      l->Draw("same");
    }
    c1->Modified();
    c1->Update();
    hitsGraphs.push_back(c1);
  }
}

void HitDistributionsLibrary::drawRPPlaneEnvelopes(unsigned int rpId)
{

  if ((rpRotations.find(rpId) != rpRotations.end()) &&  (rpTranslations.find(rpId) != rpTranslations.end())){

    CLHEP::Hep3Vector pA,pB,pC,pD,pE;
    // FIXME: these values need to match the XML geometry description
    // TODO: find a way to do it automatically
    double len = 36.070; // length of detector's edge
    double half_len = len * 0.5;
    double cut = len - 22.276 * sqrt(2.0) * 0.5; // length of the edge adjacent to the cut

    pA.setX(half_len);
    pA.setY(-half_len);

    pB.setX(half_len);
    pB.setY(half_len);

    pC.setX(-half_len);
    pC.setY(half_len);

    pD.setX(-half_len);
    pD.setY(half_len-cut);

    pE.setX(half_len-cut);
    pE.setY(-half_len);

    CLHEP::Hep3Vector pAn = rpRotations[rpId] * pA + rpTranslations[rpId];
    CLHEP::Hep3Vector pBn = rpRotations[rpId] * pB + rpTranslations[rpId];
    CLHEP::Hep3Vector pCn = rpRotations[rpId] * pC + rpTranslations[rpId];
    CLHEP::Hep3Vector pDn = rpRotations[rpId] * pD + rpTranslations[rpId];
    CLHEP::Hep3Vector pEn = rpRotations[rpId] * pE + rpTranslations[rpId];

    TLine *lAB = new TLine(pAn.x(),pAn.y(),pBn.x(),pBn.y());
    lAB->Draw();
    TLine *lBC = new TLine(pBn.x(),pBn.y(),pCn.x(),pCn.y());
    lBC->Draw();
    TLine *lCD = new TLine(pCn.x(),pCn.y(),pDn.x(),pDn.y());
    lCD->Draw();
    TLine *lDE = new TLine(pDn.x(),pDn.y(),pEn.x(),pEn.y());
    lDE->Draw();
    TLine *lEA = new TLine(pEn.x(),pEn.y(),pAn.x(),pAn.y());
    lEA->Draw();

    if(verbosity){
      double dx = rpTranslations[rpId].x();
      double dy = rpTranslations[rpId].y();
      std::cout << "rpId = " << rpId << " dx = "  << dx << " , dy = " << dy << std::endl;
    }

  } else {
    if(verbosity)
      std::cout << "Rotation or translation for RP " << rpId << " not found !" << std::endl;
  }

}


void HitDistributionsLibrary::setupCanvas()
{
//  gStyle->SetOptStat(0);
  gStyle->SetCanvasColor(0);
  gStyle->SetCanvasBorderMode(0);
  gStyle->SetFrameBorderMode(0);
  gStyle->SetOptTitle(0);
  gStyle->SetTitleBorderSize(1);
  gStyle->SetStatBorderSize(1);
  gStyle->SetFuncWidth(1);
  gStyle->SetFuncColor(4);
  gStyle->SetStatW(0.25);
  gStyle->SetStatH(0.25);
  gStyle->SetLegendBorderSize(1);
  gStyle->SetTitleYOffset(1.2);
  gStyle->SetTitleXOffset(1.2);
  gStyle->SetPadLeftMargin(0.10);
  gStyle->SetPadRightMargin(0.10);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadTopMargin(0.08);
  gStyle->SetOptTitle(0);
  gStyle->ToggleEventStatus();
}
