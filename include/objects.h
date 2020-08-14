#ifndef __OBJECTS_H__
#define __OBJECTS_H__

#include <stdint.h>
#include <ap_int.h>
#include <algorithm>
#include <utility>

using namespace std;

/*Calo-tau object definition*/
class CaloTau{

  public:
    CaloTau() : data(0) {;}

    CaloTau(ap_uint<10> tau_et, ap_uint<5> peak_phi, ap_uint<5> peak_eta, ap_uint<10> isolation, ap_uint<2> number_eg) {

      data = (tau_et) |
        (((ap_uint<32>) peak_phi)  << 10) |
        (((ap_uint<32>) peak_eta)  << 15) |
        (((ap_uint<32>) isolation) << 20) |
        (((ap_uint<32>) number_eg) << 30);
    }

    CaloTau(uint32_t i) {data = i;}

    CaloTau& operator=(const CaloTau& rhs) {
      data = rhs.data;
      return *this;
    }

    operator uint32_t() {return (uint32_t) data;}

    ap_uint<10> tau_et(){return (data & 0x3FF);}
    ap_uint<5>  peak_phi()  {return ((data >> 10) & 0x1F);}
    ap_uint<5>  peak_eta()  {return ((data >> 15) & 0x1F);}
    ap_uint<10> isolation() {return ((data >> 20) & 0x3FF);}
    ap_uint<2>  number_eg() {return ((data >> 30) & 0x3);}

#ifndef __SYNTHESIS__
    string toString() {
      return "CaloTau [" +
        to_string(this->tau_et()) + ", (" +
        to_string(this->peak_phi()) + ", " +
        to_string(this->peak_eta()) + "), " +
        to_string(this->isolation()) + ", " +
        to_string(this->number_eg()) + "] ";
    }
#endif


    ap_uint<32> data;

};


/*ECAL tower object definition*/
class Tower{

  public:
    Tower() : data(0) {;}

    Tower(ap_uint<10> cluster_et, ap_uint<10> tower_et, ap_uint<3> peak_phi, ap_uint<3> peak_eta, ap_uint<3> peak_time, ap_uint<3> hOe) {

      data = (cluster_et) | 
	(((ap_uint<32>) tower_et)  << 10) | 
	(((ap_uint<32>) peak_phi)  << 20) | 
	(((ap_uint<32>) peak_eta)  << 23) | 
	(((ap_uint<32>) peak_time) << 26) | 
	(((ap_uint<32>) hOe)       << 29);
    }

    Tower(uint32_t i) {data = i;}

    Tower& operator=(const Tower& rhs) {
      data = rhs.data;
      return *this;
    }

    operator uint32_t() {return (uint32_t) data;}

    ap_uint<10> cluster_et(){return (data & 0x3FF);}
    ap_uint<10> tower_et()  {return ((data >> 10) & 0x3FF);}
    ap_uint<3>  peak_phi()  {return ((data >> 20) & 0x7);}
    ap_uint<3>  peak_eta()  {return ((data >> 23) & 0x7);}
    ap_uint<3>  peak_time() {return ((data >> 26) & 0x7);}
    ap_uint<3>  hOe()       {return ((data >> 29) & 0x7);}

#ifndef __SYNTHESIS__
    string toString() {
      return "Tower [" + 
	to_string(this->cluster_et()) + ", (" + 
	to_string(this->peak_phi()) + ", " + 
	to_string(this->peak_eta()) + "), " + 
	to_string(this->peak_time()) + ", " + 
	to_string(this->hOe()) + ", " + 
	to_string(this->tower_et()) + "] ";
    }
#endif


    ap_uint<32> data;

};

class TowersInEta{
  public:
    TowersInEta(){
#pragma HLS ARRAY_PARTITION variable=towers complete dim=0
      for (size_t eta = 0; eta < 17; eta++) {
#pragma LOOP UNROLL
	this->towers[eta] = Tower();
      }
    }

    Tower towers[17];
};

template<typename T, int N>
void buffering(T in_[N], T out_[N]) {

  T tmp[N];
#pragma HLS ARRAY_PARTITION variable=tmp complete

  for (int i = 0; i < N; ++i) {
#pragma HLS latency min=1
#pragma HLS LOOP UNROLL
    tmp[i] = in_[i];
  }
  for (int i = 0; i < N; ++i) {
#pragma HLS latency min=1
#pragma HLS LOOP UNROLL
    out_[i] = tmp[i];
  }
}

#endif /*!__OBJECTS_H__*/
