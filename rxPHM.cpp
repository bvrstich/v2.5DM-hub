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

vector< vector<int> > *rxPHM::ph2s;
int ****rxPHM::s2ph;

/**
 * initialize the static variables and allocate and initialize the static lists
 */
void rxPHM::init(){

   ph2s = new vector< vector<int> > [2];

   s2ph = new int *** [2];

   for(int S = 0;S < 2;++S){

      s2ph[S] = new int ** [2];

      for(int S_bl = 0;S_bl < 2;++S_bl){

         s2ph[S][S_bl] = new int * [Tools::gL()];

         for(int a = 0;a < Tools::gL();++a)
            s2ph[S][S_bl][a] = new int [Tools::gL()];

      }

   }

   vector<int> v(3);

   int ph = 0;

   //S == 1/2

   //S_bl = 0;
   for(int a = 0;a < Tools::gL();++a)
      for(int b = 0;b < Tools::gL();++b){

         v[0] = 0;//S_bl
         v[1] = a;
         v[2] = b;

         ph2s[0].push_back(v);

         s2ph[0][0][a][b] = ph;

         ++ph;

      }

   //S_bl = 1
   for(int a = 0;a < Tools::gL();++a)
      for(int b = 1;b < Tools::gL();++b){

         v[0] = 1;//S_bl
         v[1] = a;
         v[2] = b;

         ph2s[0].push_back(v);

         s2ph[0][1][a][b] = ph;

         ++ph;

      }

   //S == 3/2: only S_bl = 1
   ph = 0;

   for(int a = 0;a < Tools::gL();++a)
      for(int b = 1;b < Tools::gL();++b){

         v[0] = 1;//S_bl
         v[1] = a;
         v[2] = b;

         ph2s[1].push_back(v);

         s2ph[1][1][a][b] = ph;

         ++ph;

      }

}

/**
 * deallocate the static lists
 */
void rxPHM::clear(){

   delete [] ph2s;

   for(int S = 0;S < 2;++S){

      for(int S_bl = 0;S_bl < 2;++S_bl){

         for(int a = 0;a < Tools::gL();++a)
            delete [] s2ph[S][S_bl][a];

         delete [] s2ph[S][S_bl];

      }

      delete [] s2ph[S];

   }

   delete [] s2ph;

}

/**
 * standard constructor, constructs two blocks, the S = 1/2 and S = 3/2
 * @param l the parameter l that determines which spatial sp indices are blocked out
 */
rxPHM::rxPHM() : BlockMatrix(2) {

   this->setMatrixDim(0,ph2s[0].size(),2);
   this->setMatrixDim(1,ph2s[1].size(),4);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of blockmatrix rxphm_c
 * @param rxphm_c object that will be copied into this.
 */
rxPHM::rxPHM(const rxPHM &rxphm_c) : BlockMatrix(rxphm_c){ }

/**
 * destructor
 */
rxPHM::~rxPHM(){ }

ostream &operator<<(ostream &output,const rxPHM &rxphm_p){

   for(int S = 0;S < 2;++S){

      output << endl;
      output << "S = " << 2*S + 1 << "/" << 2 << endl;
      output << endl;

      for(int i = 0;i < rxphm_p.gdim(S);++i)
         for(int j = 0;j < rxphm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t(" << rxphm_p.ph2s[S][i][0] << ")\t" << 

               rxphm_p.ph2s[S][i][1] << "\t" << rxphm_p.ph2s[S][i][2]

               << "\t;\t(" << rxphm_p.ph2s[S][j][0] << ")\t" << rxphm_p.ph2s[S][j][1]

               << "\t" << rxphm_p.ph2s[S][j][2] << "\t" 

               << rxphm_p[S](i,j) << endl;

         }

   }

   return output;

}

/**
 * @param S the dp spin
 * @param i the tp index
 * @param option == 0 return a, == 1 return b
 * @return the sp indices corresponding to the tp index i with parameter l
 */
int rxPHM::gph2s(int S,int i,int option){

   return ph2s[S][i][option];

}

/**
 * @param S the dp spin
 * @param S_ab intermediate spin
 * @param a the first sp index
 * @param b the second sp index
 * @return the tp index corresponding to the sp indices a and b with parameter l
 */
int rxPHM::gs2ph(int S,int S_ab,int a,int b){

   return s2ph[S][S_ab][a][b];

}
