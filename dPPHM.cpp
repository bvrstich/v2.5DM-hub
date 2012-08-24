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

/**
 * map a dDPM on a dPPHM using the Q1 map
 * @param ddpm input dDPM
 */
void dPPHM::Q(const dDPM &ddpm){

   int a,b,c,d;

   int S_ab,S_cd;

   int sign_ab,sign_cd;

   double norm_ab,norm_cd;

   TPM tpm;
   tpm.bar(1.0/(Tools::gN() - 2.0),ddpm);

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   double ward = 2.0 * ddpm.trace() / (Tools::gN()*(Tools::gN() - 1.0)*(Tools::gN() - 2.0));

   //first S = 1/2
   for(int i = 0;i < gdim(0);++i){

      S_ab = gt2s(0,i,0);

      a = gt2s(0,i,1);
      b = gt2s(0,i,2);

      norm_ab = 1.0;

      if(a == b)
         norm_ab /= std::sqrt(2.0);

      sign_ab = 1 - 2*S_ab;

      for(int j = i;j < gdim(0);++j){

         S_cd = gt2s(0,j,0);

         c = gt2s(0,j,1);
         d = gt2s(0,j,2);

         sign_cd = 1 - 2*S_cd;

         norm_cd = 1.0;

         if(c == d)
            norm_cd /= std::sqrt(2.0);

         //dp term
         (*this)(0,i,j) = 0.0;

         for(int S_ = 0;S_ < 2;++S_)
            (*this)(0,i,j) -= (2* (S_ + 0.5) + 1.0) * Tools::g6j(0,S_,S_ab,S_cd) * ddpm(S_,S_ab,a,b,S_cd,c,d);

         if(b == d){

            if(a == c){

               //np_a
               if(a == 0)
                  (*this)(0,i,j) += 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

               //np_d
               if(b == 0)
                  (*this)(0,i,j) += 0.5 * sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

               //sp(2)_a
               if(S_ab == S_cd)
                  (*this)(0,i,j) += norm_ab * norm_cd * spm(0,0);

            }

            //sp(1)_d
            if(b == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,c);

            //sp(3)_a first part
            if(c == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,0);

            //sp(3)_c first part
            if(a == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(c,0);

            //tp(2)_d
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_cd) * tpm(Z,a,0,c,0);

            if(a == 0)
               hard *= std::sqrt(2.0);

            if(c == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;

         }

         if(a == d){

            if(b == c){

               //np_b
               if(a == 0)
                  (*this)(0,i,j) += sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

               //np_c
               if(b == 0)
                  (*this)(0,i,j) += sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * ward;

               //sp(2)_b
               if(S_ab == S_cd)
                  (*this)(0,i,j) += sign_ab * norm_ab * norm_cd * spm(0,0);

            }

            //sp(1)_b
            if(a == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(b,c);

            //sp(3)_b second part
            if(c == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(0,b);

            //sp(3)_c second part
            if(b == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(0,c);

            //tp(2)_b
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_cd) * tpm(Z,b,0,c,0);

            if(b == 0)
               hard *= std::sqrt(2.0);

            if(c == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;

         }

         if(b == c){

            //sp(1)_c
            if(b == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_ab * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,d);

            //sp(3)_a second part
            if(d == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(a,0);

            //sp(3)_d second part
            if(a == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(d,0);

            //tp(2)_c
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_cd) * tpm(Z,a,0,d,0);

            if(a == 0)
               hard *= std::sqrt(2.0);

            if(d == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= sign_ab * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;

         }

         if(a == c){

            //sp(1)_a
            if(a == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(b,d);

            //sp(3)_b first part
            if(d == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(b,0);

            //sp(3)_d first part
            if(b == 0)
               (*this)(0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * 0.5 * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * spm(d,0);

            //tp(2)_a
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(0,Z,S_ab,S_cd) * tpm(Z,b,0,d,0);

            if(b == 0)
               hard *= std::sqrt(2.0);

            if(d == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * norm_ab * norm_cd * hard;

         }

         //tp(1)_a
         if(a == 0){

            if(0 == b)
               (*this)(0,i,j) += 0.5 * norm_ab * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,0,b,c,d);
            else
               (*this)(0,i,j) += 0.5 * norm_ab * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,0,b,c,d);

         }

         //tp(1)_b
         if(b == 0){

            if(a == 0)
               (*this)(0,i,j) += 0.5 * sign_ab * sign_cd * norm_ab * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,a,0,c,d);
            else
               (*this)(0,i,j) += 0.5 * sign_ab * sign_cd * norm_ab * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_cd,a,0,c,d);

         }

         //tp(1)_c
         if(c == 0){

            if(d == 0)
               (*this)(0,i,j) += 0.5 * norm_cd * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,0,d);
            else
               (*this)(0,i,j) += 0.5 * norm_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,0,d);

         }

         //tp(1)_d
         if(d == 0){

            if(c == 0)
               (*this)(0,i,j) += 0.5 * norm_cd * sign_ab * sign_cd * std::sqrt( 2.0 * (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,c,0);
            else
               (*this)(0,i,j) += 0.5 * norm_cd * sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * tpm(S_ab,a,b,c,0);

         }

      }
   }

   //S = 3/2 is a lot easier
   for(int i = 0;i < gdim(1);++i){

      S_ab = gt2s(1,i,0);

      a = gt2s(1,i,1);
      b = gt2s(1,i,2);

      for(int j = i;j < gdim(1);++j){

         S_cd = gt2s(1,j,0);

         c = gt2s(1,j,1);
         d = gt2s(1,j,2);

         //dp term
         (*this)(1,i,j) = 0.0;

         for(int S_ = 0;S_ < 2;++S_)
            (*this)(1,i,j) -= (2* (S_ + 0.5) + 1.0) * Tools::g6j(1,S_,S_ab,S_cd) * ddpm(S_,S_ab,a,b,S_cd,c,d);

         //sp(2) full
         if(i == j)
            (*this)(1,i,j) += spm(0,0);

         if(b == d){

            //tp(2)_d
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(1,Z,S_ab,S_cd) * tpm(Z,a,0,c,0);

            if(a == 0)
               hard *= std::sqrt(2.0);

            if(c == 0)
               hard *= std::sqrt(2.0);

            (*this)(1,i,j) -= 3.0 * hard;

         }

         if(a == d){

            //tp(2)_b
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(1,Z,S_ab,S_cd) * tpm(Z,b,0,c,0);

            if(b == 0)
               hard *= std::sqrt(2.0);

            if(c == 0)
               hard *= std::sqrt(2.0);

            (*this)(1,i,j) += 3.0 * hard;

         }

         if(b == c){

            //tp(2)_c
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(1,Z,S_ab,S_cd) * tpm(Z,a,0,d,0);

            if(a == 0)
               hard *= std::sqrt(2.0);

            if(d == 0)
               hard *= std::sqrt(2.0);

            (*this)(1,i,j) += 3.0 * hard;


         }

         if(a == c){

            //tp(2)_a
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(1,Z,S_ab,S_cd) * tpm(Z,b,0,d,0);

            if(b == 0)
               hard *= std::sqrt(2.0);

            if(d == 0)
               hard *= std::sqrt(2.0);

            (*this)(1,i,j) -= 3.0 * hard;

         }

      }
   }

   this->symmetrize();

}

/**
 * a sort of skew trace, see notes for info
 */
double dPPHM::barbreve() const {

   double ward = 0.0;

   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd){

         double hard = 0.0;

         for(int l = 0;l < Tools::gL();++l)
            for(int b = 0;b < Tools::gL();++b){

               if(l == b)
                  hard += 2.0 * (*this)(l,0,S_ab,l,b,S_cd,l,b);
               else
                  hard += (*this)(l,0,S_ab,l,b,S_cd,l,b);

            }
         
         ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * hard;

      }

   return ward;

}
