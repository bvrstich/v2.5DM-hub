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

vector< vector<int> > *xTPM::t2s;
int ****xTPM::s2t;

/**
 * initialize the static variables and allocate and initialize the static lists
 */
void xTPM::init(){

   //allocate
   t2s = new vector< vector<int> > [2];

   s2t = new int *** [2];

   for(int S = 0;S < 2;++S){

      s2t[S] = new int ** [2];

      for(int S_ab = 0;S_ab < 2;++S_ab){

         s2t[S][S_ab] = new int * [Tools::gL()];

         for(int a = 0;a < Tools::gL();++a)
            s2t[S][S_ab][a] = new int [Tools::gL()];

      }

   }

   vector<int> v(3);

   int tp;

   tp = 0;

   //S == 1/2
   for(int S_ab = 0;S_ab < 2;++S_ab){

      for(int a = 0;a < Tools::gL();++a)
         for(int b = a + S_ab;b < Tools::gL();++b){

            v[0] = S_ab;//S_ab
            v[1] = a;
            v[2] = b;

            t2s[0].push_back(v);

            s2t[0][S_ab][a][b] = tp;
            s2t[0][S_ab][b][a] = tp;

            ++tp;

         }

   }

   tp = 0;

   //S == 3/2
   for(int a = 0;a < Tools::gL();++a)
      for(int b = a + 1;b < Tools::gL();++b){

         v[0] = 1;//S_ab
         v[1] = a;
         v[2] = b;

         t2s[1].push_back(v);

         s2t[1][1][a][b] = tp;
         s2t[1][1][b][a] = tp;

         ++tp;

      }

}

/**
 * deallocate the static lists
 */
void xTPM::clear(){

   delete [] t2s;

   for(int S = 0;S < 2;++S){

      for(int S_ab = 0;S_ab < 2;++S_ab){

         for(int a = 0;a < Tools::gL();++a)
            delete [] s2t[S][S_ab][a];

         delete [] s2t[S][S_ab];

      }

      delete [] s2t[S];

   }

   delete [] s2t;

}

/**
 * standard constructor, constructs two blocks, the S = 1/2 and S = 3/2
 */
xTPM::xTPM() : BlockMatrix(2) {

   this->setMatrixDim(0,t2s[0].size(),2);
   this->setMatrixDim(1,t2s[1].size(),4);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of blockmatrix xtpm_c
 * @param xtpm_c object that will be copied into this.
 */
xTPM::xTPM(const xTPM &xtpm_c) : BlockMatrix(xtpm_c){ }

/**
 * destructor
 */
xTPM::~xTPM(){ }

ostream &operator<<(ostream &output,const xTPM &xtpm_p){

   for(int S = 0;S < 2;++S){

      output << endl;
      output << "S = " << 2*S + 1 << "/" << 2 << endl;
      output << endl;

      for(int i = 0;i < xtpm_p.gdim(S);++i)
         for(int j = 0;j < xtpm_p.gdim(S);++j){

            output << i << "\t" << j << "\t|\t(" << xtpm_p.t2s[S][i][0] << ")\t" << xtpm_p.t2s[S][i][1] << "\t" << xtpm_p.t2s[S][i][2]

               << "\t;\t(" << xtpm_p.t2s[S][j][0] << ")\t" << xtpm_p.t2s[S][j][1] << "\t" << xtpm_p.t2s[S][j][2] << "\t" 

               << xtpm_p[S](i,j) << endl;

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
int xTPM::gt2s(int S,int i,int option){

   return t2s[S][i][option];

}

/**
 * @param S the dp spin
 * @param S_ab intermediate spin
 * @param a the first sp index
 * @param b the second sp index
 * @return the tp index corresponding to the sp indices a and b with parameter l
 */
int xTPM::gs2t(int S,int S_ab,int a,int b){

   return s2t[S][S_ab][a][b];

}

/**
 * test if basis is correctly constructed
 */
void xTPM::print_basis(){

   for(int S = 0;S < 2;++S){

      cout << endl;
      cout << "S =\t" << 2*S + 1 << "/" << 2 << endl;
      cout << endl;

      for(unsigned int i = 0;i < t2s[S].size();++i)
         cout << i << "\t|\t(" << t2s[S][i][0] << ")\t" << t2s[S][i][1] << "\t" << t2s[S][i][2] << endl;

   }

}
