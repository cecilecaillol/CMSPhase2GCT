#ifndef __STITCHTOWERS_H__
#define __STITCHTOWERS_H__ 

#include "algo_top.h"
#include <algorithm>
#include <utility>

#include "algo_top_parameters.h"
#include "../../../../include/objects.h"

using namespace std;

ap_uint<1> is_gt_3x5neighbors(TowersInEta stitchedInPhi[TOWERS_IN_PHI*2], ap_uint<6> iphi, ap_uint<6> ieta){
#pragma HLS PIPELINE II=6
  if (ieta<0 or ieta>=TOWERS_IN_ETA) return 0;
  ap_uint<13> central_tower_et = stitchedInPhi[iphi].towers[ieta].tower_et();
  if (iphi-2>=0){ 
     if (stitchedInPhi[iphi-2].towers[ieta].tower_et() > central_tower_et) return 0;
     if (ieta-1>=0 && stitchedInPhi[iphi-2].towers[ieta-1].tower_et() > central_tower_et) return 0;
     if (ieta+1<=TOWERS_IN_ETA<2 && stitchedInPhi[iphi-2].towers[ieta+1].tower_et() > central_tower_et) return 0;
  }
  if (iphi-1>=0){
     if (stitchedInPhi[iphi-1].towers[ieta].tower_et() > central_tower_et) return 0;
     if (ieta-1>=0 && stitchedInPhi[iphi-1].towers[ieta-1].tower_et() > central_tower_et) return 0;
     if (ieta+1<=TOWERS_IN_ETA<2 && stitchedInPhi[iphi-1].towers[ieta+1].tower_et() > central_tower_et) return 0;
  }

  if (ieta-1>=0 && stitchedInPhi[iphi].towers[ieta-1].tower_et() > central_tower_et) return 0;
  if (ieta+1<=TOWERS_IN_ETA<2 && stitchedInPhi[iphi].towers[ieta+1].tower_et() > central_tower_et) return 0;

  if (iphi+1<=TOWERS_IN_PHI*2){
     if (stitchedInPhi[iphi+1].towers[ieta].tower_et() > central_tower_et) return 0;
     if (ieta-1>=0 && stitchedInPhi[iphi+1].towers[ieta-1].tower_et() > central_tower_et) return 0;
     if (ieta+1<=TOWERS_IN_ETA<2 && stitchedInPhi[iphi+1].towers[ieta+1].tower_et() > central_tower_et) return 0;
  }
  if (iphi+2<=TOWERS_IN_PHI*2){
     if (stitchedInPhi[iphi+2].towers[ieta].tower_et() > central_tower_et) return 0;
     if (ieta-1>=0 && stitchedInPhi[iphi+2].towers[ieta-1].tower_et() > central_tower_et) return 0;
     if (ieta+1<=TOWERS_IN_ETA<2 && stitchedInPhi[iphi+2].towers[ieta+1].tower_et() > central_tower_et) return 0;
  }
  return 1;
}

ap_uint<13> sum_five_phi(TowersInEta stitchedInPhi[TOWERS_IN_PHI*2], ap_uint<6> iphi, ap_uint<6> ieta) {
#pragma HLS PIPELINE II=6
  if (ieta<0 or ieta>=TOWERS_IN_ETA) return 0;

  ap_uint<13> e0=0;
  if (iphi-2>=0) e0=stitchedInPhi[iphi-2].towers[ieta].tower_et();
  ap_uint<13> e1=0;
  if (iphi-1>=0) e1=stitchedInPhi[iphi-1].towers[ieta].tower_et();
  ap_uint<13> e2=0;
  if (iphi>=0) e2=stitchedInPhi[iphi].towers[ieta].tower_et();
  ap_uint<13> e3=0;
  if (iphi+1<TOWERS_IN_PHI*2) e3=stitchedInPhi[iphi+1].towers[ieta].tower_et();
  ap_uint<13> e4=0;
  if (iphi+2<TOWERS_IN_PHI*2) e4=stitchedInPhi[iphi+2].towers[ieta].tower_et();

  ap_uint<16> e01 = e0 + e1;
  ap_uint<16> e23 = e2 + e3;
  ap_uint<16> e0123 = e01 + e23;
  ap_uint<16> e = e0123 + e4;
  return e;
}

ap_uint<13> get_3x5_et(TowersInEta stitchedInPhi[TOWERS_IN_PHI*2], CaloTau tau){
    ap_uint<13> left_band_e=0;
    left_band_e=sum_five_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()-1);

    ap_uint<13> central_band_e=0;
    central_band_e=sum_five_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta());

    ap_uint<13> right_band_e=0;
    right_band_e=sum_five_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()+1);

    ap_uint<13> total_e=0;

    total_e=left_band_e+central_band_e+right_band_e;

    return total_e;
}

ap_uint<13> sum_seven_phi(TowersInEta stitchedInPhi[TOWERS_IN_PHI*2], ap_uint<6> iphi, ap_uint<6> ieta) {
#pragma HLS PIPELINE II=6
  if (ieta<0 or ieta>=TOWERS_IN_ETA) return 0;

  ap_uint<13> e0=0;
  if (iphi-3>=0) e0=stitchedInPhi[iphi-3].towers[ieta].tower_et();
  ap_uint<13> e1=0;
  if (iphi-2>=0) e1=stitchedInPhi[iphi-2].towers[ieta].tower_et();
  ap_uint<13> e2=0;
  if (iphi-1>=0) e2=stitchedInPhi[iphi-1].towers[ieta].tower_et();
  ap_uint<13> e3=0;
  if (iphi>=0) e3=stitchedInPhi[iphi].towers[ieta].tower_et();
  ap_uint<13> e4=0;
  if (iphi+1<TOWERS_IN_PHI*2) e4=stitchedInPhi[iphi+1].towers[ieta].tower_et();
  ap_uint<13> e5=0;
  if (iphi+2<TOWERS_IN_PHI*2) e5=stitchedInPhi[iphi+2].towers[ieta].tower_et();
  ap_uint<13> e6=0;
  if (iphi+3<TOWERS_IN_PHI*2) e6=stitchedInPhi[iphi+3].towers[ieta].tower_et();

  ap_uint<13> e01 = e0 + e1;
  ap_uint<13> e23 = e2 + e3;
  ap_uint<13> e0123 = e01 + e23;
  ap_uint<13> e45 = e4 + e5;
  ap_uint<13> e456 = e45 + e6;
  ap_uint<13> e = e0123 + e456;
  return e;
}

ap_uint<13> get_7x7_et(TowersInEta stitchedInPhi[TOWERS_IN_PHI*2], CaloTau tau){
  ap_uint<13> e0=0;
  e0=sum_seven_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()-3);

  ap_uint<13> e1=0;
  e1=sum_seven_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()-2);

  ap_uint<13> e2=0;
  e2=sum_seven_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()-1);

  ap_uint<13> e3=0;
  e3=sum_seven_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta());

  ap_uint<13> e4=0;
  e4=sum_seven_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()+1);

  ap_uint<13> e5=0;
  e5=sum_seven_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()+2);

  ap_uint<13> e6=0;
  e6=sum_seven_phi(stitchedInPhi, tau.peak_phi(), tau.peak_eta()+3);

  ap_uint<13> e01 = e0 + e1;
  ap_uint<13> e23 = e2 + e3;
  ap_uint<13> e0123 = e01 + e23;
  ap_uint<13> e45 = e4 + e5;
  ap_uint<13> e456 = e45 + e6;
  ap_uint<13> e = e0123 + e456;
  return e;

}



#endif /*!__STITCHTOWERS_H__ */
