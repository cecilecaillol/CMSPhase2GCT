#include "algo_top_parameters.h"
#include "algo_top.h"
#include <algorithm>
#include <utility>

#include "../../../../include/objects.h"
#include "tauFunctions.h"

using namespace std;
using namespace algo;


TowersInEta unpackInputLink(hls::stream<algo::axiword576> &link) {
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME
#pragma HLS INTERFACE axis port=link
#pragma HLS INLINE

  TowersInEta tEta_;
  ap_uint<576> word_576b_;

#ifndef __SYNTHESIS__
  // Avoid simulation warnings
  if (link.empty()) return tEta_ ;
#endif

  word_576b_ = link.read().data;

  tEta_.towers[0]  = Tower(word_576b_( 31,   0));
  tEta_.towers[1]  = Tower(word_576b_( 63,  32));
  tEta_.towers[2]  = Tower(word_576b_( 95,  64));
  tEta_.towers[3]  = Tower(word_576b_(127,  96));
  tEta_.towers[4]  = Tower(word_576b_(159, 128));
  tEta_.towers[5]  = Tower(word_576b_(191, 160));
  tEta_.towers[6]  = Tower(word_576b_(223, 192));
  tEta_.towers[7]  = Tower(word_576b_(255, 224));
  tEta_.towers[8]  = Tower(word_576b_(287, 256));
  tEta_.towers[9]  = Tower(word_576b_(319, 288));
  tEta_.towers[10] = Tower(word_576b_(351, 320));
  tEta_.towers[11] = Tower(word_576b_(383, 352));
  tEta_.towers[12] = Tower(word_576b_(415, 384));
  tEta_.towers[13] = Tower(word_576b_(447, 416));
  tEta_.towers[14] = Tower(word_576b_(479, 448));
  tEta_.towers[15] = Tower(word_576b_(511, 480));
  tEta_.towers[16] = Tower(word_576b_(543, 512));

  return tEta_;
}

bool packOutput(TowersInEta tEta_, hls::stream<algo::axiword576> &olink){
#pragma HLS PIPELINE II=N_OUTPUT_WORDS_PER_FRAME
#pragma HLS INTERFACE axis port=link
#pragma HLS INLINE

  ap_uint<576> word_576b_;

  word_576b_( 31,   0) = (ap_uint<32>) tEta_.towers[0].data;
  word_576b_( 63,  32) = (ap_uint<32>) tEta_.towers[1].data;
  word_576b_( 95,  64) = (ap_uint<32>) tEta_.towers[2].data;
  word_576b_(127,  96) = (ap_uint<32>) tEta_.towers[3].data;
  word_576b_(159, 128) = (ap_uint<32>) tEta_.towers[4].data;
  word_576b_(191, 160) = (ap_uint<32>) tEta_.towers[5].data;
  word_576b_(223, 192) = (ap_uint<32>) tEta_.towers[6].data;
  word_576b_(255, 224) = (ap_uint<32>) tEta_.towers[7].data;
  word_576b_(287, 256) = (ap_uint<32>) tEta_.towers[8].data;
  word_576b_(319, 288) = (ap_uint<32>) tEta_.towers[9].data;
  word_576b_(351, 320) = (ap_uint<32>) tEta_.towers[10].data;
  word_576b_(383, 352) = (ap_uint<32>) tEta_.towers[11].data;
  word_576b_(415, 384) = (ap_uint<32>) tEta_.towers[12].data;
  word_576b_(447, 416) = (ap_uint<32>) tEta_.towers[13].data;
  word_576b_(479, 448) = (ap_uint<32>) tEta_.towers[14].data;
  word_576b_(511, 480) = (ap_uint<32>) tEta_.towers[15].data;
  word_576b_(543, 512) = (ap_uint<32>) tEta_.towers[16].data;
  word_576b_(575, 544) = (ap_uint<32>) 0;

  axiword576 r; r.last = 0; r.user = 0;
  r.data = word_576b_;
  olink.write(r);

  return true;
}


void algo_top(hls::stream<axiword576> link_in[N_INPUT_LINKS], hls::stream<axiword576> link_out[N_OUTPUT_LINKS]) {
#pragma HLS INTERFACE axis port=link_in
#pragma HLS INTERFACE axis port=link_out
#pragma HLS PIPELINE II=N_WORDS_PER_FRAME

#pragma HLS ARRAY_PARTITION variable=link_in complete dim=0
#pragma HLS ARRAY_PARTITION variable=link_out complete dim=0


  // Step 1: Unpack links
  // Input is 64 links carrying 32phix34eta towers
  TowersInEta towersInPosEta[TOWERS_IN_PHI];
  TowersInEta towersInNegEta[TOWERS_IN_PHI];
#pragma HLS ARRAY_PARTITION variable=towersInPosEta complete dim=0
#pragma HLS ARRAY_PARTITION variable=towersInNegEta complete dim=0
     
  for (size_t ilink = 0; ilink < N_INPUT_LINKS/2; ilink++) {
#pragma LOOP UNROLL
#pragma HLS latency min=1
    size_t iPosEta = ilink;
    size_t iNegEta = ilink + (N_INPUT_LINKS/2);
    towersInPosEta[ilink] = unpackInputLink(link_in[iPosEta]);
    towersInNegEta[ilink] = unpackInputLink(link_in[iNegEta]);
  }

  // Step 2: Tau Algo goes here

  CaloTau taus[MAX_TAUS];
#pragma HLS ARRAY_PARTITION variable=taus  complete dim=0
  for(int n=0; n<MAX_TAUS; n++){
#pragma LOOP UNROLL
#pragma HLS latency min=1
    taus[n]=CaloTau(0,0,0,0,0);
  }

  // First save highest seed towers with Et and position
  for(int tphi=0; tphi<TOWERS_IN_PHI; tphi++){
#pragma LOOP UNROLL
#pragma HLS latency min=1
    for(int teta=0; teta<TOWERS_IN_ETA/2; teta++){
#pragma LOOP UNROLL
#pragma HLS latency min=1
        if(towersInPosEta[tphi].towers[teta].tower_et() > 0 && is_gt_3x5neighbors(towersInPosEta,tphi,teta)){ // only take the tower that has the largest energy of the 3x5 cluster
	   if (towersInPosEta[tphi].towers[teta].tower_et()>taus[0].tau_et() ){
	      taus[2]=taus[1];
              taus[1]=taus[0];
	      taus[0]=CaloTau(towersInPosEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
           else if(towersInPosEta[tphi].towers[teta].tower_et()>taus[1].tau_et() ){
              taus[2]=taus[1];
              taus[1]=CaloTau(towersInPosEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
           else if(towersInPosEta[tphi].towers[teta].tower_et()>taus[2].tau_et() ){
              taus[2]=CaloTau(towersInPosEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
	   teta++; // No need to check the neighbor tower because it has less energy and will be merged in the same cluster
	}
     }
  }

  // For each seed compute the cluster energy (3x5) and the isolation (7x7)
  if (taus[0].peak_eta()>=0){
    taus[0]=CaloTau(get_3x5_et(towersInPosEta, taus[0]), taus[0].peak_phi(), taus[0].peak_eta(), get_7x7_et(towersInPosEta, taus[0]), 0);
  }

  if (taus[1].peak_eta()>=0){
    taus[1]=CaloTau(get_3x5_et(towersInPosEta, taus[1]), taus[1].peak_phi(), taus[1].peak_eta(), get_7x7_et(towersInPosEta, taus[1]), 0);
  }

  if (taus[2].peak_eta()>=0){
    taus[2]=CaloTau(get_3x5_et(towersInPosEta, taus[2]), taus[2].peak_phi(), taus[2].peak_eta(), get_7x7_et(towersInPosEta, taus[2]), 0);
  }


  // Step 3: Pack the outputs
  for (size_t olink = 0; olink < N_OUTPUT_LINKS/2; olink++) {
#pragma LOOP UNROLL
#pragma HLS latency min=1
    size_t iPosEta = olink;              
    size_t iNegEta = olink + (N_OUTPUT_LINKS/2);
    packOutput(towersInPosEta[olink], link_out[iPosEta]);
    packOutput(towersInNegEta[olink], link_out[iNegEta]);
  }
}
