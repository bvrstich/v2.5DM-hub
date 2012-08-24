#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cmath>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;

#include "include.h"

vector< vector<int> > PHM::ph2s;
int **PHM::s2ph;

/**
 * initialize the static variables and lists
 */
void PHM::init(){

   //allocate
   s2ph = new int * [Tools::gL()];

   for(int a = 0;a < Tools::gL();++a)
      s2ph[a] = new int [Tools::gL()];

   //initialize
   int i = 0;

   vector<int> v(2);

   for(int a = 0;a < Tools::gL();++a)
      for(int b = 0;b < Tools::gL();++b){

         v[0] = a;
         v[1] = b;

         ph2s.push_back(v);
         s2ph[a][b] = i;

         ++i;

      }

}

/**
 * deallocate the statics
 */
void PHM::clear(){

   for(int a = 0;a < Tools::gL();++a)
      delete [] s2ph[a];

   delete [] s2ph;

}

/**
 * standard constructor: constructs BlockMatrix object with 2 blocks, for S = 0 and 1 of dimension M*M/4.
 */
PHM::PHM() : BlockMatrix(2) {

   //set the dimension of the blocks
   this->setMatrixDim(0,ph2s.size(),1);
   this->setMatrixDim(1,ph2s.size(),3);

}

/**
 * copy constructor: constructs BlockMatrix object with two blocks of dimension M*M/4 and copies the content of phm_c into it,
 * if counter == 0, allocates and constructs the lists containing the relationship between sp and ph basis.
 * @param phm_c PHM to be copied into (*this)
 */
PHM::PHM(const PHM &phm_c) : BlockMatrix(phm_c){ }

/**
 * destructor: if counter == 1 the memory for the static lists ph2s en s2ph will be deleted.
 */
PHM::~PHM(){ }

/**
 * access the elements of the matrix in sp mode, 
 * @param S The spin of the block you want to access
 * @param a first sp index that forms the ph row index i in block S together with b
 * @param b second sp index that forms the ph row index i in block S together with a
 * @param c first sp index that forms the ph column index j in block S together with d
 * @param d second sp index that forms the ph column index j in block S together with c
 * @return the number on place PHM(S,i,j)
 */
double &PHM::operator()(int S,int a,int b,int c,int d){

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(S,i,j);

}

/**
 * access the elements of the matrix in sp mode: read only mode
 * @param S The spin of the block you want to access
 * @param a first sp index that forms the ph row index i in block S together with b
 * @param b second sp index that forms the ph row index i in block S together with a
 * @param c first sp index that forms the ph column index j in block S together with d
 * @param d second sp index that forms the ph column index j in block S together with c
 * @return the number on place PHM(S,i,j)
 */
double PHM::operator()(int S,int a,int b,int c,int d) const{

   int i = s2ph[a][b];
   int j = s2ph[c][d];

   return (*this)(S,i,j);

}

ostream &operator<<(ostream &output,const PHM &phm_p){

   for(int S = 0;S < phm_p.gnr();++S){

      output << S << "\t" << phm_p.gdim(S) << "\t" << phm_p.gdeg(S) << std::endl;
      output << std::endl;

      for(int i = 0;i < phm_p.gdim(S);++i)
         for(int j = 0;j < phm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t" << phm_p.ph2s[i][0] << "\t" << phm_p.ph2s[i][1]

               << "\t" << phm_p.ph2s[j][0] << "\t" << phm_p.ph2s[j][1] << "\t" << phm_p(S,i,j) << endl;

         }

      output << endl;

   }

   return output;

}

/**
 * The G map, maps a TPM object on a PHM object.
 * @param tpm input TPM
 */
void PHM::G(const TPM &tpm){

   //construct the SPM corresponding to the TPM
   SPM spm;
   spm.bar(1.0/(Tools::gN() - 1.0),tpm);

   int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         a = ph2s[i][0];
         b = ph2s[i][1];

         for(int j = i;j < gdim(S);++j){

            c = ph2s[j][0];
            d = ph2s[j][1];

            //tp part
            (*this)(S,i,j) = -Tools::g6j(0,0,S,0)*tpm(0,a,d,c,b) - 3.0*Tools::g6j(0,0,S,1)*tpm(1,a,d,c,b);

            //norm
            if(a == d)
               (*this)(S,i,j) *= std::sqrt(2.0);

            if(c == b)
               (*this)(S,i,j) *= std::sqrt(2.0);

            //sp part
            if(b == d)
               (*this)(S,i,j) += spm(a,c);

         }
      }

   }

   this->symmetrize();

}

/** 
 * map a dDPM object on a PHM object by summing over the spin of elements with one index diagonal.
 * see symmetry.pdf for more info. Watch out, not symmetrical in general!
 * @param scale the number you scale the PHM with
 * @param ddpm input dDPM
 */
void PHM::spinsum(double scale,const dDPM &ddpm){

  int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         a = ph2s[i][0];
         b = ph2s[i][1];

         for(int j = 0;j < gdim(S);++j){

            c = ph2s[j][0];
            d = ph2s[j][1];

            (*this)(S,i,j) = 0.0;

            //only S = 1/2 contribution
            for(int S_ab = 0;S_ab < 2;++S_ab)
               (*this)(S,i,j) += 2.0 * std::sqrt( (2.0*S_ab + 1.0) / (2.0*S + 1.0) ) * Tools::g6j(0,0,S,S_ab) * ddpm(a,0,S_ab,a,b,S,c,d);

            (*this)(S,i,j) *= scale;

         }
      }
   }

}

/** 
 * map a dPPHM object on a PHM object by summing over the spin of elements with one index diagonal.
 * see symmetry.pdf for more info. Watch out, not symmetrical in general!
 * @param scale the number you scale the PHM with
 * @param dpphm input dPPHM
 */
void PHM::spinsum(double scale,const dPPHM &dpphm){

  int a,b,c,d;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < this->gdim(S);++i){

         a = ph2s[i][0];
         b = ph2s[i][1];

         for(int j = 0;j < this->gdim(S);++j){

            c = ph2s[j][0];
            d = ph2s[j][1];

            (*this)(S,i,j) = 0.0;

            //only S = 1/2 contribution
            for(int S_ab = 0;S_ab < 2;++S_ab)
               (*this)(S,i,j) += std::sqrt( (2.0*S_ab + 1.0) / (2.0*S + 1.0) ) * dpphm(a,0,S_ab,a,b,S,c,d);

            (*this)(S,i,j) *= scale;

         }
      }
   }

}
