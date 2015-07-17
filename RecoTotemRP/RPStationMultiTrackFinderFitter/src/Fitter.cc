#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Fitter.h"

namespace RPStationMultiTrackFinderFitter {

Fitter::Fitter(
        const edm::ParameterSet& ps,
        GeometryHelper &geometry) :
    verbosity(ps.getUntrackedParameter<unsigned int>("verbosity")),
    geometry(geometry),
    enableFitting(ps.getParameter<bool>("enableFitting")),
    z_h(ps.getParameter<double>("z_h")) {}

Fitter::~Fitter() {}

void Fitter::fitCollection( RPStationTrackFitCollection *rpTrColl, const vector<Cluster> *input) {
    for (vector<Cluster>::const_iterator cluIt = input->begin();
            cluIt != input->end(); cluIt++) {
        RPStationTrackFit track;
        if (enableFitting) {
            const vector<RPRecoHit> hits(cluIt->hits.begin(), cluIt->hits.end());
            fit(track, hits);
        } else {
            clusterToTrack(track, *cluIt);
        }
        (*rpTrColl)[1].push_back(track);
    }
}

void Fitter::fitCollection(
        RPStationTrackFitCollection *rpTrColl,
        const RPTrackCandidateCollection *input) {
    if (input->size() < 2) return;

    vector<RPRecoHit> hits;
    for (RPTrackCandidateCollection::const_iterator rpIt = input->begin();
            rpIt != input->end(); rpIt++) {
        if (rpIt->second.Fittable() != true) return;
        const vector<RPRecoHit> &trackHits = rpIt->second.TrackRecoHits();
        hits.insert(hits.end(), trackHits.begin(), trackHits.end());
    }

    RPStationTrackFit track;

    fit(track, hits);
    (*rpTrColl)[1].push_back(track);
}

void Fitter::clusterToTrack(RPStationTrackFit &out, const Cluster &cluster) {
    out.z0 = z_h;
    out.ax = cluster.center(0);
    out.bx = cluster.center(1);
    out.ay = cluster.center(2);
    out.by = cluster.center(3);

    out.covarianceMatrix = cluster.sum_W_i;
    out.valid = true;
}

void Fitter::fit( RPStationTrackFit &out, const vector<RPRecoHit> &hits) {
    int matrixSize = hits.size();
    TVectorD M(matrixSize);
    TMatrixD G(matrixSize, 4);
    TMatrixD Vi(matrixSize, matrixSize);

    for (int i=0; i < matrixSize; i++)
    {
        RPRecoHit hit = hits[i];
        DetInfo detInfo = geometry.getDetInfo(hit.DetId());

        //////////////////////////////////////////////////////////////////////
        // 							matrix M								//
        //////////////////////////////////////////////////////////////////////

        M(i) = hit.Position() + detInfo.cx*detInfo.dx + detInfo.cy*detInfo.dy;

        //////////////////////////////////////////////////////////////////////
        // 							matrix G								//
        //////////////////////////////////////////////////////////////////////

        G(i, 0) = detInfo.dx * (detInfo.cz - z_h);
        G(i, 1) = detInfo.dy * (detInfo.cz - z_h);
        G(i, 2) = detInfo.dx;
        G(i, 3) = detInfo.dy;

        //////////////////////////////////////////////////////////////////////
        // 							matrix Vi								//
        //////////////////////////////////////////////////////////////////////

        Vi(i, i) = 1.0/(hit.Sigma()*hit.Sigma());
    }

    //matrix calculation
    TMatrixD GTVi(G, TMatrixD::kTransposeMult, Vi);
    TMatrixD GTViG(GTVi, TMatrixD::kMult, G);
    TMatrixD GTViG_i(TMatrixD::kInverted, GTViG);
    TVectorD P = GTViG_i * GTVi * M;

    //chi^2
    TVectorD D = M - G * P;
    double chiSq = 0.;
    for (int i = 0; i < D.GetNrows(); i++)
        for (int j = 0; j < D.GetNrows(); j++)
            chiSq += D(i) * Vi(i, j) * D(j);

    //filling output
    out.covarianceMatrix = GTViG_i;
    out.chiSq = chiSq;

    out.ndf= matrixSize - 4;
    out.ax = P(0);
    out.ay = P(1);
    out.bx = P(2);
    out.by = P(3);
    out.z0 = z_h;

}

} // namespace RPStationMultiTrackFinderFitter
