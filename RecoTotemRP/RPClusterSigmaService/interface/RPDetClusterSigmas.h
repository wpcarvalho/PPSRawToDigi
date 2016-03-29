#ifndef RecoTotemRP_RPClusterSigmaService_RPDetClusterSigmas_h
#define RecoTotemRP_RPClusterSigmaService_RPDetClusterSigmas_h

class RPDetClusterSigmas
{
  public:
    RPDetClusterSigmas();

    inline double GetClusterSigma(unsigned int rp_id, int clu_size, double avg_strip_no) const
    {
      return 0.0191;
    }
};

#endif
