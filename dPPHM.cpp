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
dPPHM::dPPHM() : xTPM() { }

/**
 * copy constructor 
 * @param dpphm object that will be copied into this.
 */
dPPHM::dPPHM(const dPPHM &dpphm) : xTPM(dpphm) { }

/**
 * destructor
 */
dPPHM::~dPPHM(){ }

/**
 * @return the trace of the dPPHM object
 */
double dPPHM::trace() const {

   return Tools::gL() * xTPM::trace();

}

/**
 * @return inproduct of (*this) dpphm with dpphm_i, defined as Tr (A B)
 * @param dpphm_i input dpphm
 */
double dPPHM::ddot(const dPPHM &dpphm_i) const{

   return Tools::gL() * xTPM::ddot(dpphm_i);

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
double dPPHM::operator()(int S,int S_ab,int a,int b,int S_cd,int c,int d) const {

   int phase_i = get_inco(S,S_ab,a,b);

   if(phase_i == 0)
      return 0.0;

   int phase_j = get_inco(S,S_cd,c,d);

   if(phase_j == 0)
      return 0.0;

   return phase_i * phase_j * (*this)(S,gs2t(S,S_ab,a,b),gs2t(S,S_cd,c,d));

}

/**
 * access the numbers in sp mode: with shift for the block index!
 * @param l block-index
 * @param S dp spin
 * @param S_ab intermediate spin for a and b
 * @param a first spatial index for the row
 * @param b second spatial index for the row
 * @param S_cd intermediate spin for c and d
 * @param c first spatial index for the column
 * @param d second spatial index for the column
 */
double dPPHM::operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const {

   return (*this)(S,S_ab,Tools::shift(a,l),Tools::shift(b,l),S_cd,Tools::shift(c,l),Tools::shift(d,l));

}

/**
 * @param S dp spin
 * @param S_ab intermediate spin
 * @param a first sp index
 * @param b second sp index
 * @return the right phase for this order of sp indices as a function of my basis.
 */
int dPPHM::get_inco(int S,int S_ab,int a,int b){

   if(S == 0){//S == 1/2

      if(a == b && S_ab == 1)
         return 0;

      if(a > b)
         return 1 - 2*S_ab;
      else
         return 1;

   }
   else{//S == 3/2

      //intermediate spin can never be 0
      if(S_ab == 0)
         return 0;

      //totally antisymmetrical
      if(a == b)
         return 0;

      if(a > b)
         return -1;
      else
         return 1;

   }

}

/**
 * map a dDPM on a dPPHM using the I2 map
 * @param ddpm input dDPM
 */
void dPPHM::I(const dDPM &ddpm){

   int a,b,c,d;

   int S_ab,S_cd;

   TPM tpm;
   tpm.bar(1.0/(Tools::gN() - 2.0),ddpm);

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         S_ab = gt2s(S,i,0);

         a = gt2s(S,i,1);
         b = gt2s(S,i,2);

         for(int j = i;j < gdim(S);++j){

            S_cd = gt2s(S,j,0);

            c = gt2s(S,j,1);
            d = gt2s(S,j,2);

            (*this)(S,i,j) = 0.0;

            for(int S_ = 0;S_ < 2;++S_)
               (*this)(S,i,j) += (2* (S_ + 0.5) + 1.0) * Tools::g6j(S,S_,S_ab,S_cd) * ddpm(S_,S_ab,a,b,S_cd,c,d);

            if(S_ab == S_cd)
               (*this)(S,i,j) += tpm(S_ab,a,b,c,d);

         }
      }
   }

   this->symmetrize();

}
