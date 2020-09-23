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

  return tEta_;
}

bool packOutput(CaloTau taus[5], hls::stream<algo::axiword576> &olink){
#pragma HLS PIPELINE II=N_OUTPUT_WORDS_PER_FRAME
#pragma HLS INTERFACE axis port=link
#pragma HLS INLINE

  ap_uint<576> word_576b_;

  word_576b_( 33,   0) = (ap_uint<34>) taus[0].data;
  word_576b_( 67,  34) = (ap_uint<34>) taus[1].data;
  word_576b_( 101,  68) = (ap_uint<34>) taus[2].data;
  word_576b_(135,  102) = (ap_uint<34>) taus[3].data;
  word_576b_(169, 136) = (ap_uint<34>) taus[4].data;
  word_576b_(575, 170) = 0;

  axiword576 r; r.last = 0; r.user = 0;
  r.data = word_576b_;
  olink.write(r);

  return true;
}

bool packOutputZero(hls::stream<algo::axiword576> &olink){
#pragma HLS PIPELINE II=N_OUTPUT_WORDS_PER_FRAME
#pragma HLS INTERFACE axis port=link
#pragma HLS INLINE

  ap_uint<576> word_576b_;

  word_576b_(575, 0) = 0;

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
  // Input is 60 links carrying 30phix30eta towers x 2 pos/neg
  TowersInEta towersEta[TOWERS_IN_PHI*2];
#pragma HLS ARRAY_PARTITION variable=towersEta complete dim=0
     
  for (size_t ilink = 0; ilink < N_INPUT_LINKS/2; ilink++) {
#pragma LOOP UNROLL
#pragma HLS latency min=1
    towersEta[ilink] = unpackInputLink(link_in[ilink]);
    towersEta[ilink+TOWERS_IN_PHI] = unpackInputLink(link_in[ilink+TOWERS_IN_PHI]);
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
  for(int tphi=0; tphi<TOWERS_IN_PHI*2; tphi++){
#pragma LOOP UNROLL
#pragma HLS latency min=1
    for(int teta=0; teta<TOWERS_IN_ETA; teta++){
#pragma LOOP UNROLL
#pragma HLS latency min=1
        if(towersEta[tphi].towers[teta].tower_et() > 0 && is_gt_3x5neighbors(towersEta,tphi,teta)){ // only take the tower that has the largest energy of the 3x5 cluster
	   if (towersEta[tphi].towers[teta].tower_et()>taus[0].tau_et() ){
              taus[4]=taus[3];
              taus[3]=taus[2];
	      taus[2]=taus[1];
              taus[1]=taus[0];
	      taus[0]=CaloTau(towersEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
           else if(towersEta[tphi].towers[teta].tower_et()>taus[1].tau_et() ){
              taus[4]=taus[3];
              taus[3]=taus[2];
              taus[2]=taus[1];
              taus[1]=CaloTau(towersEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
           else if(towersEta[tphi].towers[teta].tower_et()>taus[2].tau_et() ){
              taus[4]=taus[3];
              taus[3]=taus[2];
              taus[2]=CaloTau(towersEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
           else if(towersEta[tphi].towers[teta].tower_et()>taus[3].tau_et() ){
              taus[4]=taus[3];
              taus[3]=CaloTau(towersEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
           else if(towersEta[tphi].towers[teta].tower_et()>taus[4].tau_et() ){
              taus[4]=CaloTau(towersEta[tphi].towers[teta].tower_et(),tphi,teta,0,0);
           }
	   teta++; // No need to check the neighbor tower because it has less energy and will be merged in the same cluster
	}
     }
  }

 #ifndef __SYNTHESIS__  
   for(int tphi=0; tphi<TOWERS_IN_PHI*2; tphi++){
     for(int teta=0; teta<TOWERS_IN_ETA; teta++){
       if(towersEta[tphi].towers[teta].cluster_et() != 0 )
        cout<<std::dec<<"[tphi, teta] = ["<<tphi<<", "<<teta<<"]: "<<towersEta[tphi].towers[teta].toString()<<endl;
     }
   }
 #endif

 #ifndef __SYNTHESIS__
        cout<<std::dec<<taus[0].tau_et()<<" "<<taus[1].tau_et()<<" "<<taus[2].tau_et()<<" "<<taus[3].tau_et()<<" "<<taus[4].tau_et()<<endl;
        cout<<std::dec<<taus[0].peak_eta()<<" "<<taus[1].peak_eta()<<" "<<taus[2].peak_eta()<<" "<<taus[3].peak_eta()<<" "<<taus[4].peak_eta()<<endl;
        cout<<std::dec<<taus[0].peak_phi()<<" "<<taus[1].peak_phi()<<" "<<taus[2].peak_phi()<<" "<<taus[3].peak_phi()<<" "<<taus[4].peak_phi()<<endl;
 #endif

  // For each seed compute the cluster energy (3x5) and the isolation (7x7)
  if (taus[0].peak_eta()>=0){
    taus[0]=CaloTau(get_3x5_et(towersEta, taus[0]), taus[0].peak_phi(), taus[0].peak_eta(), get_7x7_et(towersEta, taus[0]), 0);
  }

  if (taus[1].peak_eta()>=0){
    taus[1]=CaloTau(get_3x5_et(towersEta, taus[1]), taus[1].peak_phi(), taus[1].peak_eta(), get_7x7_et(towersEta, taus[1]), 0);
  }

  if (taus[2].peak_eta()>=0){
    taus[2]=CaloTau(get_3x5_et(towersEta, taus[2]), taus[2].peak_phi(), taus[2].peak_eta(), get_7x7_et(towersEta, taus[2]), 0);
  }

  if (taus[3].peak_eta()>=0){
    taus[3]=CaloTau(get_3x5_et(towersEta, taus[3]), taus[3].peak_phi(), taus[3].peak_eta(), get_7x7_et(towersEta, taus[3]), 0);
  }

  if (taus[4].peak_eta()>=0){
    taus[4]=CaloTau(get_3x5_et(towersEta, taus[4]), taus[4].peak_phi(), taus[4].peak_eta(), get_7x7_et(towersEta, taus[4]), 0);
  }

 #ifndef __SYNTHESIS__
        cout<<std::dec<<taus[0].tau_et()<<" "<<taus[1].tau_et()<<" "<<taus[2].tau_et()<<" "<<taus[3].tau_et()<<" "<<taus[4].tau_et()<<endl;
        cout<<std::dec<<taus[0].peak_eta()<<" "<<taus[1].peak_eta()<<" "<<taus[2].peak_eta()<<" "<<taus[3].peak_eta()<<" "<<taus[4].peak_eta()<<endl;
        cout<<std::dec<<taus[0].peak_phi()<<" "<<taus[1].peak_phi()<<" "<<taus[2].peak_phi()<<" "<<taus[3].peak_phi()<<" "<<taus[4].peak_phi()<<endl;
 #endif


  // Step 3: Pack the outputs
  /*for (size_t olink = 0; olink < N_OUTPUT_LINKS/2; olink++) {
#pragma LOOP UNROLL
#pragma HLS latency min=1
    size_t iPosEta = olink;              
    size_t iNegEta = olink + (N_OUTPUT_LINKS/2);
    packOutput(towersInPosEta[olink], link_out[iPosEta]);
    packOutput(towersInNegEta[olink], link_out[iNegEta]);
  }*/
  packOutput(taus, link_out[0]); // pack the 5 taus in a single link
  for (size_t olink = 1; olink < N_OUTPUT_LINKS; olink++) {
#pragma LOOP UNROLL
#pragma HLS latency min=1
    packOutputZero(link_out[olink]);
  }
}
