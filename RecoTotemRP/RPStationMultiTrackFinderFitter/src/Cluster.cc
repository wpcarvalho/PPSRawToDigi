#include "RecoTotemRP/RPStationMultiTrackFinderFitter/interface/Cluster.h"

namespace RPStationMultiTrackFinderFitter {

double Cluster::Distance(
        const TVectorD &v, const TMatrixD &V_v, const TMatrixD &V_v_i) const {
    // (inverted) uncertainty matrix for difference v - center
    TMatrixD V_diff_i = V_v + sum_W_i;
    V_diff_i.SetTol(1.e-22);
    V_diff_i.Invert();

    double d2 = 0.;
    unsigned int dim = sum_Wv.GetNrows();
    for (unsigned int i = 0; i < dim; i++)
        for (unsigned int j = 0; j < dim; j++)
            d2 += (v(i) - center(i)) * V_diff_i(i, j) * (v(j) - center(j));

    return d2;
}

void Cluster::Add(RPId rpId1, RPId rpId2, const TVectorD &v, const TMatrixD &W) {
    rpIdCount[rpId1]++;
    rpIdCount[rpId2]++;

    // W = V^-1, V = covariance matrix for v
    TVectorD Wv = W * v;

    unsigned int dim = sum_Wv.GetNrows();
    for (unsigned int i = 0; i < dim; i++)
    {
        sum_Wv(i) += Wv(i);
        for (unsigned int j = 0; j < dim; j++)
            sum_W(i, j) += W(i, j);
    }

    sum_W_i = sum_W;
    sum_W_i.Invert();

    center = sum_W_i * sum_Wv;
}

Cluster& Cluster::Copy(Cluster &copy) {
    copy.rpIdCount.insert(rpIdCount.begin(), rpIdCount.end());

    copy.sum_W = sum_W;
    copy.sum_W_i = sum_W_i;
    copy.sum_Wv = sum_Wv;
    copy.center = center;

    return copy;
}


} // namespace RPStationMultiTrackFinderFitter
