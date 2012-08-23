#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

/**
 * constructor, makes matrix of dimension L, there will be two degenerate blocks +1/2,-1/2. So
 * only one matrix is needed of dimension L is needed to represent the SPM.
 */
SPM::SPM() : Matrix(Tools::gL()) { }

/**
 * copy constructor
 * @param spm_copy content of this matrix will be copied into the constructed matrix
 */
SPM::SPM(const SPM &spm_copy) : Matrix(spm_copy) { }

/**
 * destructor
 */
SPM::~SPM(){

}

ostream &operator<<(ostream &output,const SPM &spm_p){

   for(int i = 0;i < spm_p.gn();++i)
      for(int j = 0;j < spm_p.gn();++j)
         output << i << "\t" << j << "\t" << spm_p(i,j) << endl;

   return output;

}

/**
 * Trace out a set of indices to create the "bar" matrix of a TPM
 * @param scale the factor u want the SPM to be scaled with (1/N-1 for normal sp density matrix)
 * @param tpm the TPM out of which the SPM will be filled
 */
void SPM::bar(double scale,const TPM &tpm){

   //hulpvariabele
   double ward;

   for(int a = 0;a < Tools::gL();++a)
      for(int c = a;c < Tools::gL();++c){

         (*this)(a,c) = 0.0;

         for(int b = 0;b < Tools::gL();++b){

            //S = 0 stuk
            ward = tpm(0,a,b,c,b);

            if(a == b)
               ward *= std::sqrt(2.0);

            if(c == b)
               ward *= std::sqrt(2.0);

            (*this)(a,c) += ward;

            //S = 1 stuk: hier kan nooit a = b en c = d wegens antisymmetrie
            (*this)(a,c) += 3.0*tpm(1,a,b,c,b);

         }

         //nog schalen
         (*this)(a,c) *= 0.5*scale;

      }

   this->symmetrize();

}

/**
 * map a dDPM matrix on a SPM by tracing the first indices together with the diagonal third (see symmetry.pdf)
 * @param scale the SPM with this number
 * @param ddpm input dDPM object
 */
void SPM::breve(double scale,const dDPM &ddpm){

   double ward,hard;

   for(int b = 0;b < Tools::gL();++b)
      for(int d = b;d < Tools::gL();++d){

         (*this)(b,d) = 0.0;

         //only S = 1/2 contribution
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int l = 0;l < Tools::gL();++l){

                  hard = ddpm(l,0,S_ab,l,b,S_cd,l,d);

                  if(b == l)
                     hard *= std::sqrt(2.0);

                  if(d == l)
                     hard *= std::sqrt(2.0);

                  ward += hard;

               }

               (*this)(b,d) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd) * ward;

            }

         (*this)(b,d) *= scale;

      }

   this->symmetrize();

}
