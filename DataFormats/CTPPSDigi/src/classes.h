/****************************************************************************
* Seyed Mohsen Etesami
****************************************************************************/

#include "DataFormats/Common/interface/Wrapper.h"
#include "DataFormats/Common/interface/DetSet.h"
#include "DataFormats/Common/interface/DetSetVector.h"

#include "DataFormats/CTPPSDigi/interface/DiamondDigi.h"
#include "DataFormats/CTPPSDigi/interface/DiamondVFATStatus.h"
#include "DataFormats/CTPPSDigi/interface/DiamondFEDInfo.h"

#include <vector>

namespace {
  namespace {
    DiamondDigi rm_diamo_dig;
    edm::DetSet<DiamondDigi> ds_rp_diamo_dig;
    std::vector<DiamondDigi> vec_rp_diamo_dig;
    edm::DetSetVector<DiamondDigi> dsv_rp_diamo_dig;
    std::vector<edm::DetSet<DiamondDigi> > vec_ds_rp_diamo_dig;
    edm::Wrapper<edm::DetSet<DiamondDigi> > wds_rp_diamo_dig;
    edm::Wrapper<edm::DetSetVector<DiamondDigi> > wdsv_rp_diamo_dig;


    std::map<unsigned int, uint64_t> dummy37;

    DiamondVFATStatus dummy40;
    edm::Wrapper< DiamondVFATStatus > dummy41;
    edm::DetSetVector<DiamondVFATStatus> dummy42;
    edm::Wrapper< edm::DetSetVector<DiamondVFATStatus> > dummy43;

    std::bitset<8> dummy60;
    edm::Wrapper< std::bitset<8> > dummy61;

    std::bitset<16> dummy70;
    edm::Wrapper< std::bitset<16> > dummy71;

    DiamondFEDInfo di;
    std::vector<DiamondFEDInfo> v_di;
    edm::Wrapper<std::vector<DiamondFEDInfo>> w_v_di;
  }
}
