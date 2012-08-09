#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

/**
 * standard constructor
 */
dDPM::dDPM() : rxTPM() { }

/**
 * copy constructor 
 * @param ddpm object that will be copied into this.
 */
dDPM::dDPM(const dDPM &ddpm) : rxTPM(ddpm) { }

/**
 * destructor
 */
dDPM::~dDPM(){ }

/**
 * @return the trace of the dDPM object
 */
double dDPM::trace() const {

   return Tools::gL() * rxTPM::trace();

}


/**
 * @return inproduct of (*this) ddpm with ddpm_i, defined as Tr (A B)
 * @param ddpm_i input ddpm
 */
double dDPM::ddot(const dDPM &ddpm_i) const{

   return Tools::gL() * rxTPM::ddot(ddpm_i);

}

/**
 * access the numbers in sp mode
 * @param S dp spin
 * @param S_ab intermediate spin for a and b
 * @param a first spatial index for the row
 * @param b second spatial index for the row
 * @param S_cd intermediate spin for c and d
 * @param c first spatial index for the column
 * @param d second spatial index for the column
 */
double dDPM::operator()(int S,int S_ab,int a,int b,int S_cd,int c,int d) const {

   int phase_i = get_inco(S,S_ab,a,b);

   if(phase_i == 0)
      return 0.0;

   int phase_j = get_inco(S,S_cd,c,d);

   if(phase_j == 0)
      return 0.0;

   return phase_i * phase_j * (*this)(S,rxTPM::gs2t(S,S_ab,a,b),rxTPM::gs2t(S,S_cd,c,d));

}

/**
 * access the numbers in sp mode: using translational invariance for the block index l!
 * @param l the block index
 * @param S dp spin
 * @param S_ab intermediate spin for a and b
 * @param a first spatial index for the row
 * @param b second spatial index for the row
 * @param S_cd intermediate spin for c and d
 * @param c first spatial index for the column
 * @param d second spatial index for the column
 */
double dDPM::operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const {

   return (*this)(S,S_ab, (a - l + Tools::gL())%Tools::gL() , (b - l + Tools::gL())%Tools::gL() ,S_cd, (c - l + Tools::gL())%Tools::gL() , (d - l + Tools::gL())%Tools::gL() );

}

/**
 * @param S dp spin
 * @param S_ab intermediate spin
 * @param a first sp index
 * @param b second sp index
 * @return the right phase for this order of sp indices as a function of my basis.
 */
int dDPM::get_inco(int S,int S_ab,int a,int b){

   //3 indices can never be equal at the same time!
   if(a == 0 && b == 0)
      return 0;

   if(S == 0){//S == 1/2

      if(a == b && S_ab == 1)
         return 0;

      if(a > b)
         return 1 - 2*S_ab;
      else
         return 1;

   }
   else{//S == 3/2

      //here totally antisymmetrical
      if(a == b || a == 0 || b == 0)
         return 0;

      //intermediate spin can never be 0
      if(S_ab == 0)
         return 0;

      if(a > b)
         return -1;
      else
         return 1;

   }

}
