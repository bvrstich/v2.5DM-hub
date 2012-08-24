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

vector< vector<int> > *TPM::t2s;
int ***TPM::s2t;

/**
 * initialize the static variables and allocate and initialize the static lists
 */
void TPM::init(){

   //allocate
   t2s = new vector< vector<int> > [2];

   s2t = new int ** [2];

   for(int S = 0;S < 2;++S){

      s2t[S] = new int * [Tools::gL()];

      for(int a = 0;a < Tools::gL();++a)
         s2t[S][a] = new int [Tools::gL()];

   }

   vector<int> v(2);

   int i;

   for(int S = 0;S < 2;++S){

      i = 0;

      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + S;b < Tools::gL();++b){

            v[0] = a;
            v[1] = b;

            t2s[S].push_back(v);

            s2t[S][a][b] = i;
            s2t[S][b][a] = i;

            ++i;

         }

   }

}

/**
 * deallocate the static lists
 */
void TPM::clear(){

   delete [] t2s;

   for(int S = 0;S < 2;++S){

      for(int a = 0;a < Tools::gL();++a)
         delete [] s2t[S][a];

      delete [] s2t[S];

   }

   delete [] s2t;

}

/**
 * standard constructor for a spinsymmetrical tp matrix: constructs BlockMatrix object with 2 blocks, for S = 0 or 1,
 */
TPM::TPM() : BlockMatrix(2) {

   //set the dimension and degeneracy of the two blocks:
   this->setMatrixDim(0,t2s[0].size(),1);
   this->setMatrixDim(1,t2s[1].size(),3);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of matrix tpm_c
 * @param tpm_c object that will be copied into this.
 */
TPM::TPM(const TPM &tpm_c) : BlockMatrix(tpm_c){ }

/**
 * destructor
 */
TPM::~TPM(){ }

/**
 * access the elements of the the blocks in sp mode, the symmetry or antisymmetry of the blocks is automatically accounted for:\n\n
 * Antisymmetrical for S = 1, symmetrical in the sp orbitals for S = 0\n\n
 * @param S The spinquantumnumber that identifies the block
 * @param a first sp index that forms the tp row index i of spin S, together with b
 * @param b second sp index that forms the tp row index i of spin S, together with a
 * @param c first sp index that forms the tp column index j of spin S, together with d
 * @param d second sp index that forms the tp column index j of spin S, together with c
 * @return the number on place TPM(B,i,j) with the right phase.
 */
double TPM::operator()(int S,int a,int b,int c,int d) const{

   if(S == 0){

      int i = s2t[0][a][b];
      int j = s2t[0][c][d];

      return (*this)(S,i,j);

   }
   else{

      if( (a == b) || (c == d) )
         return 0;
      else{

         int i = s2t[1][a][b];
         int j = s2t[1][c][d];

         int phase = 1;

         if(a > b)
            phase *= -1;
         if(c > d)
            phase *= -1;

         return phase*(*this)(S,i,j);

      }

   }

}

ostream &operator<<(ostream &output,const TPM &tpm_p){

   for(int S = 0;S < 2;++S){

      output << S << "\t" << tpm_p.gdim(S) << "\t" << tpm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < tpm_p.gdim(S);++i)
         for(int j = 0;j < tpm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t" << tpm_p.t2s[S][i][0] << "\t" << tpm_p.t2s[S][i][1]

               << "\t" << tpm_p.t2s[S][j][0] << "\t" << tpm_p.t2s[S][j][1] << "\t" << tpm_p(S,i,j) << endl;

         }

      std::cout << std::endl;

   }

   return output;

}

/**
 * @return The expectation value of the total spin for the TPM.
 */
double TPM::spin() const{

   double ward = 0.0;

   for(int i = 0;i < this->gdim(0);++i)
      ward += -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) * (*this)(0,i,i);

   for(int i = 0;i < this->gdim(1);++i)
      ward += 3.0 * ( -1.5 * (Tools::gN() - 2.0)/(Tools::gN() - 1.0) + 2.0 ) * (*this)(1,i,i);

   return ward;

}

/**
 * make a TPM object from a dDPM object by tracing the third index.
 * @param scale factor to scale the TPM with
 * @param ddpm input dDPM
 */
void TPM::bar(double scale,const dDPM &ddpm){

   int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         a = t2s[S][i][0];
         b = t2s[S][i][1];

         for(int j = i;j < gdim(S);++j){

            c = t2s[S][j][0];
            d = t2s[S][j][1];

            (*this)(S,i,j) = 0;

            for(int l = 0;l < Tools::gL();++l)
               for(int Z = 0;Z < 2;++Z)
                  (*this)(S,i,j) += (2.0*(Z + 0.5)+ 1.0) * ddpm(l,Z,S,a,b,S,c,d);

            (*this)(S,i,j) *= scale/(2.0*S + 1.0);

         }
      }

   }

   this->symmetrize();

}

/**
 * make a TPM object from a dPPHM object by tracing the third index.
 * @param scale factor to scale the TPM with
 * @param dpphm input dPPHM
 */
void TPM::bar(double scale,const dPPHM &dpphm){

   int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         a = t2s[S][i][0];
         b = t2s[S][i][1];

         for(int j = i;j < gdim(S);++j){

            c = t2s[S][j][0];
            d = t2s[S][j][1];

            (*this)(S,i,j) = 0;

            for(int l = 0;l < Tools::gL();++l)
               for(int Z = 0;Z < 2;++Z)
                  (*this)(S,i,j) += (2.0*(Z + 0.5)+ 1.0) * dpphm(l,Z,S,a,b,S,c,d);

            (*this)(S,i,j) *= scale/(2.0*S + 1.0);

         }
      }

   }

   this->symmetrize();

}
