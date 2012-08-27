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
dPHHM::dPHHM() : rxPHM() { }

/**
 * copy constructor 
 * @param dphhm object that will be copied into this.
 */
dPHHM::dPHHM(const dPHHM &dphhm) : rxPHM(dphhm) { }

/**
 * destructor
 */
dPHHM::~dPHHM(){ }

/**
 * @return the trace of the dPHHM object
 */
double dPHHM::trace() const {

   return Tools::gL() * rxPHM::trace();

}

/**
 * @return inproduct of (*this) dphhm with dphhm_i, defined as Tr (A B)
 * @param dphhm_i input dphhm
 */
double dPHHM::ddot(const dPHHM &dphhm_i) const{

   return Tools::gL() * rxPHM::ddot(dphhm_i);

}

/**
 * access the numbers in sp mode
 * @param S pph spin
 * @param S_bl intermediate spin for b and l
 * @param a first spatial index for the row
 * @param b second spatial index for the row
 * @param S_dl intermediate spin for d and l
 * @param c first spatial index for the column
 * @param d second spatial index for the column
 */
double dPHHM::operator()(int S,int S_bl,int a,int b,int S_dl,int c,int d) const {

   if( S == 1 && ( S_bl == 0 || S_dl == 0) )
      return 0.0;

   if(S_bl == 1 && b == 0)
      return 0.0;

   if(S_dl == 1 && d == 0)
      return 0.0;

   return (*this)(S,gs2ph(S,S_bl,a,b),gs2ph(S,S_dl,c,d));

}

/**
 * access the numbers in sp mode
 * @param l the diagonal third index
 * @param S pph spin
 * @param S_bl intermediate spin for b and l
 * @param a first spatial index for the row
 * @param b second spatial index for the row
 * @param S_dl intermediate spin for d and l
 * @param c first spatial index for the column
 * @param d second spatial index for the column
 */
double dPHHM::operator()(int l,int S,int S_bl,int a,int b,int S_dl,int c,int d) const {

   return (*this)(S,S_bl,Tools::shift(a,l),Tools::shift(b,l),S_dl,Tools::shift(c,l),Tools::shift(d,l));

}

/**
 * Map a dDPM on a dPHHM using the G1 map
 * @param ddpm input dDPM
 */
void dPHHM::G1(const dDPM &ddpm){

   int a,b,c,d;
   int S_bl,S_dl;

   int sign_bl,sign_dl;

   TPM tpm;
   tpm.bar(1.0/(Tools::gN() - 2.0),ddpm);

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   //first S = 1/2
   for(int i = 0;i < gdim(0);++i){

      S_bl = gph2s(0,i,0);

      a = gph2s(0,i,1);
      b = gph2s(0,i,2);

      sign_bl = 1 - 2*S_bl;

      for(int j = i;j < gdim(0);++j){

         S_dl = gph2s(0,j,0);

         c = gph2s(0,j,1);
         d = gph2s(0,j,2);

         sign_dl = 1 - 2*S_dl;

         //let the games begin: first the dp
         (*this)(0,i,j) = 0.0;

         for(int S_ = 0;S_ < 2;S_++)
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  (*this)(0,i,j) += (2*(S_ + 0.5) + 1.0) * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) )

                     * Tools::g6j(S_,0,S_dl,S_ab) * Tools::g6j(S_,0,S_bl,S_cd) * Tools::g6j(S_,0,S_bl,S_dl)

                     * ddpm(S_,S_ab,a,d,S_cd,c,b);

               }

         if(a == d)
            (*this)(0,i,j) *= std::sqrt(2.0);

         if(c == b)
            (*this)(0,i,j) *= std::sqrt(2.0);

         //tp(1)
         if(S_bl == S_dl)
            if(a == c){

               double hard = tpm(S_bl,b,0,d,0);

               if(b == 0)
                  hard *= std::sqrt(2.0);

               if(d == 0)
                  hard *= std::sqrt(2.0);

               (*this)(0,i,j) += hard;

            }

         //part only in S = 1/2 blocks
         if(a == b){

            //sp_d
            if(c == d)
               (*this)(0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * sign_dl * spm(0,0);

            //sp_b
            if(c == 0)
               (*this)(0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * spm(d,0);

            //tp(2)_a
            double hard = tpm(S_dl,c,0,d,0);

            if(c == 0)
               hard *= std::sqrt(2.0);

            if(d == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * sign_dl * hard;

         }

         if(c == d){

            //sp_c
            if(a == 0)
               (*this)(0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_dl * spm(b,0);

            //tp(2)_d
            double hard = tpm(S_bl,a,0,b,0);

            if(a == 0)
               hard *= std::sqrt(2.0);

            if(b == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * sign_dl * hard;

         }

         if(a == 0){

            //sp_a
            if(c == 0)
               (*this)(0,i,j) += 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * spm(b,d);

            //tp(2)_b
            double hard = tpm(S_dl,c,b,d,0);

            if(c == b)
               hard *= std::sqrt(2.0);

            if(d == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_dl * hard;

         }

         if(c == 0){

            //tp(2)_c
            double hard = tpm(S_bl,a,d,b,0);

            if(a == d)
               hard *= std::sqrt(2.0);

            if(b == 0)
               hard *= std::sqrt(2.0);

            (*this)(0,i,j) -= 0.5 * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * sign_bl * hard;

         }

      }
   }

   //then S = 3/2
   for(int i = 0;i < gdim(1);++i){

      S_bl = gph2s(1,i,0);

      a = gph2s(1,i,1);
      b = gph2s(1,i,2);

      for(int j = i;j < gdim(1);++j){

         S_dl = gph2s(1,j,0);

         c = gph2s(1,j,1);
         d = gph2s(1,j,2);

         //dp
         (*this)(1,i,j) = 0.0;

         for(int S_ = 0;S_ < 2;S_++)
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  (*this)(1,i,j) += (2*(S_ + 0.5) + 1.0) * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) )

                     * Tools::g6j(S_,0,S_dl,S_ab) * Tools::g6j(S_,0,S_bl,S_cd) * Tools::g6j(S_,1,S_bl,S_dl)

                     * ddpm(S_,S_ab,a,d,S_cd,c,b);

               }

         if(a == d)
            (*this)(1,i,j) *= std::sqrt(2.0);

         if(c == b)
            (*this)(1,i,j) *= std::sqrt(2.0);

         //tp(1)
         if(a == c)
            (*this)(1,i,j) += tpm(1,b,0,d,0);

      }
   }

   this->symmetrize();

}

/**
 * @return the skew-trace of the this dPHHM object, see notes for info (phd)
 */
double dPHHM::skew_trace() const{

   double ward;

   double dspm = 0.0;

   for(int S_bl = 0;S_bl < 2;++S_bl)
      for(int S_dl = 0;S_dl < 2;++S_dl){

         ward = 0.0;

         for(int a = 0;a < Tools::gL();++a)
            for(int c = 0;c < Tools::gL();++c)
               ward += (*this)(0,S_bl,a,a,S_dl,c,c);

         dspm += std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * (1 - 2*S_bl) * (1 - 2*S_dl) * ward;

      }

   return dspm;

}

/**
 * Map a dDPM on a dPHHM using the G2 map
 * @param ddpm input dDPM
 */
void dPHHM::G2(const dDPM &ddpm){
   
   int a,b,c,d;
   int S_bl,S_dl;

   int sign_bl,sign_dl;

   TPM tpm;
   tpm.bar(1.0/(Tools::gN() - 2.0),ddpm);

   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         S_bl = gph2s(S,i,0);

         a = gph2s(S,i,1);
         b = gph2s(S,i,2);

         sign_bl = 1 - 2*S_bl;

         for(int j = i;j < gdim(S);++j){

            S_dl = gph2s(S,j,0);

            c = gph2s(S,j,1);
            d = gph2s(S,j,2);

            sign_dl = 1 - 2*S_dl;

            //first the dp
            (*this)(S,i,j) = 0.0;

            for(int S_ = 0;S_ < 2;S_++)
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)(S,i,j) -= (2*(S_ + 0.5) + 1.0) * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) )

                        * Tools::g6j(S_,0,S_dl,S_ab) * Tools::g6j(S_,0,S_bl,S_cd) * Tools::g6j(S_,S,S_bl,S_dl)

                        * ddpm(S_,S_ab,a,d,S_cd,c,b);

                  }

            if(a == d)
               (*this)(S,i,j) *= std::sqrt(2.0);

            if(c == b)
               (*this)(S,i,j) *= std::sqrt(2.0);

            //tp_a
            double hard = 0.0;

            for(int Z = 0;Z < 2;++Z)
               hard += (2*Z + 1.0) * Tools::g9j(S,Z,S_dl,S_bl) * tpm(Z,a,d,c,b);

            if(a == d)
               hard *= std::sqrt(2.0);

            if(c == b)
               hard *= std::sqrt(2.0);

            (*this)(S,i,j) -= std::sqrt( (2*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * hard;

            if(b == d){

               //sp
               if(S_bl == S_dl){

                  //sp_a
                  (*this)(S,i,j) += spm(a,c);

                  //sp_b
                  if(b == 0)
                     (*this)(S,i,j) += sign_bl * spm(a,c);

               }

               //tp_d
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::g9j(S,Z,S_dl,S_bl) * tpm(Z,a,0,c,0);

               if(a == 0)
                  hard *= std::sqrt(2.0);

               if(c == 0)
                  hard *= std::sqrt(2.0);

               (*this)(S,i,j) -= sign_bl * sign_dl * std::sqrt( (2*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * hard;

            }

            if(b == 0){

               //tp_b
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::g9j(S,Z,S_dl,S_bl) * tpm(Z,a,d,c,0);

               if(a == d)
                  hard *= std::sqrt(2.0);

               if(c == 0)
                  hard *= std::sqrt(2.0);

               (*this)(S,i,j) -= sign_bl * std::sqrt( (2*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * hard;

            }

            if(d == 0){

               //tp_c
               double hard = 0.0;

               for(int Z = 0;Z < 2;++Z)
                  hard += (2*Z + 1.0) * Tools::g9j(S,Z,S_dl,S_bl) * tpm(Z,a,0,c,b);

               if(a == 0)
                  hard *= std::sqrt(2.0);

               if(c == b)
                  hard *= std::sqrt(2.0);

               (*this)(S,i,j) -= sign_dl * std::sqrt( (2*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * hard;

            }

         }
      }
   }

   this->symmetrize();

}
