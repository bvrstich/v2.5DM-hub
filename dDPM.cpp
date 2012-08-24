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

   return (*this)(S,S_ab, Tools::shift(a,l) , Tools::shift(b,l) ,S_cd, Tools::shift(c,l) , Tools::shift(d,l) );

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

/**
 * project (*this) on a dDPM object with correct symmetry with the third index
 */
void dDPM::proj_W(){

   int i,j;

   int phase_i,phase_j;

   //storage for symmetries with itsself
   double mat[2][2];
   double vec[2];

   //first the S = 1/2 , the difficult part

   //all equal: W^l|aa;aa: three equalities!
   for(int a = 1;a < Tools::gL();++a){

      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int S_cd = 0;S_cd < 2;++S_cd)
            mat[S_ab][S_cd] = (*this)(0,0,S_ab,0,a,S_cd,0,a);

      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int S_cd = S_ab;S_cd < 2;++S_cd){

            //1
            double ward = mat[S_ab][S_cd];

            //2
            for(int S_lb = 0;S_lb < 2;++S_lb)
               ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * Tools::g6j(0,0,S_ab,S_lb) * mat[S_lb][S_cd];

            //3
            for(int S_ld = 0;S_ld < 2;++S_ld)
               ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * Tools::g6j(0,0,S_cd,S_ld) * mat[S_ab][S_ld];

            //4
            for(int S_lb = 0;S_lb < 2;++S_lb)
               for(int S_ld = 0;S_ld < 2;++S_ld){

                  ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) 

                     * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * mat[S_lb][S_ld];

               }

            //5
            ward += 2.0*std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 
            
               * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,0) * (1 - 2*S_ab) * (1 - 2*S_cd) * (*this)(a,0,0,0,0,0,0,0);

            i = rxTPM::gs2t(0,S_ab,0,a);
            j = rxTPM::gs2t(0,S_cd,0,a);

            (*this)(0,i,j) = 0.2 * ward;
            (*this)(0,j,i) = (*this)(0,i,j);

         }

      i = rxTPM::gs2t(0,0,Tools::par(a),Tools::par(a));

      (*this)(0,i,i) = 0.0;

      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int S_cd = 0;S_cd < 2;++S_cd){

            (*this)(0,i,i) += 0.5*std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,0) 

               * (1 - 2*S_ab) * (1 - 2*S_cd) * (*this)(0,S_ab,0,a,S_cd,0,a);

         }

   }

   //other with 3 equalities: a=b, b=c and l = d
   for(int a = 1;a < Tools::gL();++a){

      if(a == Tools::gL()/2){//if a == L/2, par(a) = a and all the symmetries involve only the same indices

         for(int S_cd = 0;S_cd < 2;++S_cd)
            vec[S_cd] = (*this)(0,0,a,a,S_cd,a,0);

         //first the "in block symmetries"
         for(int S_cd = 0;S_cd < 2;++S_cd){

            //1)
            double ward = vec[S_cd];

            //2)
            for(int S_cl = 0;S_cl < 2;++S_cl)
               ward += std::sqrt( (2.0*S_cl + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl) * (1 - 2*S_cd) * Tools::g6j(0,0,S_cd,S_cl) * vec[S_cl];

             //3)
            for(int S_lb = 0;S_lb < 2;++S_lb)
               ward += 2.0*std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) ) * Tools::g6j(0,0,0,S_lb) * Tools::g6j(0,0,S_cd,0) * vec[S_lb];

            i = rxTPM::gs2t(0,0,a,a);
            j = rxTPM::gs2t(0,S_cd,a,0);

            phase_j = 1 - 2*S_cd;

            (*this)(0,i,j) = 0.25 * phase_j * ward;
            (*this)(0,j,i) = (*this)(0,i,j);

         }

      }
      else{

         for(int S_cd = 0;S_cd < 2;++S_cd)
            vec[S_cd] = (*this)(0,0,a,a,S_cd,a,0);

         //first the "in block symmetries"
         for(int S_cd = 0;S_cd < 2;++S_cd){

            //1)
            double ward = vec[S_cd];

            //2)
            for(int S_cl = 0;S_cl < 2;++S_cl)
               ward += std::sqrt( (2.0*S_cl + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl) * (1 - 2*S_cd) * Tools::g6j(0,0,S_cd,S_cl) * vec[S_cl];

            //3)
            for(int S_lb = 0;S_lb < 2;++S_lb)
               ward += 2.0*std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) ) * Tools::g6j(0,0,0,S_lb) * Tools::g6j(0,0,S_cd,0) * (*this)(a,0,S_lb,0,a,0,0,0);

            i = rxTPM::gs2t(0,0,a,a);
            j = rxTPM::gs2t(0,S_cd,0,a);

            (*this)(0,i,j) = 0.25 * (1 - 2*S_cd) * ward;
            (*this)(0,j,i) = (*this)(0,i,j);

         }

         //then the "a" block: which is of course projected back onto in the '0' block now
         for(int S_lb = 0;S_lb < 2;++S_lb){

            i = rxTPM::gs2t(0,S_lb,0,Tools::par(a));
            j = rxTPM::gs2t(0,0,Tools::par(a),Tools::par(a));

            phase_i = 1 - 2*S_lb;

            (*this)(0,i,j) = 0.0;

            for(int S_cd = 0;S_cd < 2;++S_cd){

               (*this)(0,i,j) += phase_i * std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) ) 

                  * Tools::g6j(0,0,0,S_lb) * Tools::g6j(0,0,S_cd,0) * (*this)(0,0,a,a,S_cd,a,0);

            }

            (*this)(0,j,i) = (*this)(0,i,j);

         }

      }

   }

   //then 2 equalities, first ( a = c ; b = d ) W^l_{ab;ab} <--> W^a_{lb;lb} <--> W^b_{al;al}
   for(int a = 1;a < Tools::gL();++a){

      for(int b = a + 1;b < Tools::gL();++b){

         if(Tools::par(a) == b && Tools::shift(b,a) == a){//then the elements are mappen onto themselves

            //first store the elements
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd)
                  mat[S_ab][S_cd] = (*this)(0,S_ab,a,b,S_cd,a,b);

            //first calculate the avarage:
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = S_ab;S_cd < 2;++S_cd){

                  double ward = mat[S_ab][S_cd];

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) )

                           * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (1 - 2*S_lb) * (1 - 2*S_ld) * mat[S_lb][S_ld];

                     }

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_al)

                           * Tools::g6j(0,0,S_cd,S_cl)  * (1 - 2*S_ab) * (1 - 2*S_cd) * mat[S_al][S_cl];

                     }

                  i = rxTPM::gs2t(0,S_ab,a,b);
                  j = rxTPM::gs2t(0,S_cd,a,b);

                  (*this)(0,i,j) = ward/3.0;
                  (*this)(0,j,i) = (*this)(0,i,j);

               }

         }
         else{

            //first calculate the avarage:
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = S_ab;S_cd < 2;++S_cd){

                  double ward = (*this)(0,S_ab,a,b,S_cd,a,b);

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) )

                           * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(a,0,S_lb,0,b,S_ld,0,b);

                     }

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_al)

                           * Tools::g6j(0,0,S_cd,S_cl) * (1 - 2*S_al) * (1 - 2*S_cl) * (1 - 2*S_ab) * (1 - 2*S_cd) * (*this)(b,0,S_al,a,0,S_cl,a,0);

                     }

                  i = rxTPM::gs2t(0,S_ab,a,b);
                  j = rxTPM::gs2t(0,S_cd,a,b);

                  (*this)(0,i,j) = ward/3.0;
                  (*this)(0,j,i) = (*this)(0,i,j);

               }

            //then make the rest symmetric
            for(int S_lb = 0;S_lb < 2;++S_lb)
               for(int S_ld = S_lb;S_ld < 2;++S_ld){

                  i = rxTPM::gs2t(0,S_lb,Tools::par(a),Tools::shift(b,a));
                  j = rxTPM::gs2t(0,S_ld,Tools::par(a),Tools::shift(b,a));

                  if(Tools::par(a) > Tools::shift(b,a)){

                     phase_i = 1 - 2*S_lb;
                     phase_j = 1 - 2*S_ld;

                  }
                  else{

                     phase_i = 1;
                     phase_j = 1;

                  }

                  (*this)(0,i,j) = 0.0;

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     for(int S_cd = 0;S_cd < 2;++S_cd){

                        (*this)(0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) ) 

                           * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(0,S_ab,a,b,S_cd,a,b);

                     }

                  (*this)(0,j,i) = (*this)(0,i,j);

               }

            for(int S_al = 0;S_al < 2;++S_al)
               for(int S_cl = S_al;S_cl < 2;++S_cl){

                  i = rxTPM::gs2t(0,S_al,Tools::shift(a,b),Tools::par(b));
                  j = rxTPM::gs2t(0,S_cl,Tools::shift(a,b),Tools::par(b));

                  if(Tools::shift(a,b) > Tools::par(b)){

                     phase_i = 1 - 2*S_al;
                     phase_j = 1 - 2*S_cl;

                  }
                  else{

                     phase_i = 1;
                     phase_j = 1;

                  }

                  (*this)(0,i,j) = 0.0;

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     for(int S_cd = 0;S_cd < 2;++S_cd){

                        (*this)(0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) * 

                           (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al) * (1 - 2*S_cl) 

                           * Tools::g6j(0,0,S_ab,S_al) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(0,S_ab,a,b,S_cd,a,b);

                     }

                  (*this)(0,j,i) = (*this)(0,i,j);

               }

         }

      }
   }

   //next with two equalities:
   for(int a = 1;a < Tools::gL();++a){

      //a = d and  b = l: W^l_{al;ca}
      for(int c = 1;c < a;++c){

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               mat[S_ab][S_cd] = (*this)(0,S_ab,a,0,S_cd,c,a);

         //first set the "0" block right
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               //1)
               double ward = mat[S_ab][S_cd];

               //2)
               for(int S_al = 0;S_al < 2;++S_al)
                  ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * Tools::g6j(0,0,S_al,S_ab) * mat[S_al][S_cd];

               //3)
               for(int S_cl = 0;S_cl < 2;++S_cl){

                  ward += std::sqrt( 2.0 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) 

                     * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(a,0,0,0,0,S_cl,c,0);

               }

               i = rxTPM::gs2t(0,S_ab,a,0);
               j = rxTPM::gs2t(0,S_cd,c,a);

               phase_i = 1 - 2*S_ab;

               (*this)(0,i,j) = phase_i * ward / 3.0;
               (*this)(0,j,i) = (*this)(0,i,j);

            }

         //then the "a" block
         for(int S_cl = 0;S_cl < 2;++S_cl){

            i = rxTPM::gs2t(0,0,Tools::par(a),Tools::par(a));
            j = rxTPM::gs2t(0,S_cl,Tools::shift(c,a),Tools::par(a));

            if(Tools::shift(c,a) > Tools::par(a))
               phase_j = 1 - 2*S_cl;
            else
               phase_j = 1;

            (*this)(0,i,j) = 0.0;

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  (*this)(0,i,j) += phase_j * std::sqrt( 0.5 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cl) * (1 - 2*S_cd)

                     * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cl,S_cd) * (*this)(0,S_ab,a,0,S_cd,c,a);

               }

            (*this)(0,j,i) = (*this)(0,i,j);

         }

      }

      //a = c and  b = l: W^l_{al;ad}
      for(int d = a + 1;d < Tools::gL();++d){

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               mat[S_ab][S_cd] = (*this)(0,S_ab,a,0,S_cd,a,d);

         //first set the "l" block right
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               //1)
               double ward = mat[S_ab][S_cd];

               for(int S_al = 0;S_al < 2;++S_al)
                  ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * Tools::g6j(0,0,S_al,S_ab) * mat[S_al][S_cd];

               for(int S_ld = 0;S_ld < 2;++S_ld){

                  ward += std::sqrt( 2.0 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) )

                     * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(a,0,0,0,0,S_ld,0,d);

               }

               i = rxTPM::gs2t(0,S_ab,a,0);
               j = rxTPM::gs2t(0,S_cd,a,d);

               phase_i = 1 - 2*S_ab;

               (*this)(0,i,j) = phase_i * ward / 3.0;
               (*this)(0,j,i) = (*this)(0,i,j);

            }

         //then the "a" block
         for(int S_ld = 0;S_ld < 2;++S_ld){

            i = rxTPM::gs2t(0,0,Tools::par(a),Tools::par(a));
            j = rxTPM::gs2t(0,S_ld,Tools::par(a),Tools::shift(d,a));

            if(Tools::par(a) > Tools::shift(d,a))
               phase_j = 1 - 2*S_ld;
            else
               phase_j = 1;

            (*this)(0,i,j) = 0.0;

            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  (*this)(0,i,j) += phase_j * std::sqrt( 0.5 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) 

                     * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_ld,S_cd) * (*this)(0,S_ab,a,0,S_cd,a,d);

               }

            (*this)(0,j,i) = (*this)(0,i,j);

         }

      }

   }
   
   //another one with 2 equalities: a = c = l
   for(int b = 1;b < Tools::gL();++b)
      for(int d = b + 1;d < Tools::gL();++d){

         //save the numbers first in the matrix mat
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               mat[S_ab][S_cd] = (*this)(0,S_ab,0,b,S_cd,0,d);

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               i = rxTPM::gs2t(0,S_ab,0,b);
               j = rxTPM::gs2t(0,S_cd,0,d);

               double ward = mat[S_ab][S_cd];

               for(int S_lb = 0;S_lb < 2;++S_lb)
                  ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * Tools::g6j(0,0,S_ab,S_lb) * mat[S_lb][S_cd];

               for(int S_ld = 0;S_ld < 2;++S_ld)
                  ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * Tools::g6j(0,0,S_cd,S_ld) * mat[S_ab][S_ld];

               for(int S_lb = 0;S_lb < 2;++S_lb)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) )

                        * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * mat[S_lb][S_ld];

                  }

               (*this)(0,i,j) = 0.25 * ward;
               (*this)(0,j,i) = (*this)(0,i,j);

            }

      }
   
   //one equality: start with a == c
   for(int a = 1;a < Tools::gL();++a){

      for(int b = a + 1;b < Tools::gL();++b){

         for(int d = b + 1;d < Tools::gL();++d){

            //first take the average
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  //1)
                  double ward = (*this)(0,S_ab,a,b,S_cd,a,d);

                  //2)
                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) ) 

                           * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(a,0,S_lb,0,b,S_ld,0,d);

                     }

                  i = rxTPM::gs2t(0,S_ab,a,b);
                  j = rxTPM::gs2t(0,S_cd,a,d);

                  (*this)(0,i,j) = 0.5*ward;
                  (*this)(0,j,i) = (*this)(0,i,j);

               }

            //then symmetrize the other term
            for(int S_lb = 0;S_lb < 2;++S_lb)
               for(int S_ld = 0;S_ld < 2;++S_ld){

                  i = rxTPM::gs2t(0,S_lb,Tools::par(a),Tools::shift(b,a));
                  j = rxTPM::gs2t(0,S_ld,Tools::par(a),Tools::shift(d,a));

                  if(Tools::par(a) > Tools::shift(b,a))
                     phase_i = 1 - 2*S_lb;
                  else
                     phase_i = 1;

                  if(Tools::par(a) > Tools::shift(d,a))
                     phase_j = 1 - 2*S_ld;
                  else
                     phase_j = 1;

                  (*this)(0,i,j) = 0.0;

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     for(int S_cd = 0;S_cd < 2;++S_cd){

                        (*this)(0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) )

                           * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(0,S_ab,a,b,S_cd,a,d);

                     }

                  (*this)(0,j,i) = (*this)(0,i,j);

               }

         }
      }
   }

   //then b == c: ab;bd
   for(int a = 1;a < Tools::gL();++a){

      for(int b = a + 1;b < Tools::gL();++b){

         for(int d = b + 1;d < Tools::gL();++d){

            if(b == Tools::par(b) && Tools::shift(a,b) == d && Tools::shift(d,b) == a){//then it is projected onto itsself!

               //store original elements
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd)
                     mat[S_ab][S_cd] = (*this)(0,S_ab,a,b,S_cd,b,d);

               //first take average
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     //1)
                     double ward = mat[S_ab][S_cd];

                     //2)
                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_ld)

                              * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cd,S_ld) * mat[S_ld][S_al];

                        }

                     i = rxTPM::gs2t(0,S_ab,a,b);
                     j = rxTPM::gs2t(0,S_cd,b,d);

                     (*this)(0,i,j) = 0.5 * ward;
                     (*this)(0,j,i) = (*this)(0,i,j);

                  }

            }
            else{

               //first take average
               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     //1)
                     double ward = (*this)(0,S_ab,a,b,S_cd,b,d);

                     //2)
                     for(int S_al = 0;S_al < 2;++S_al)
                        for(int S_ld = 0;S_ld < 2;++S_ld){

                           ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al)

                              * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(b,0,S_al,a,0,S_ld,0,d);

                        }

                     i = rxTPM::gs2t(0,S_ab,a,b);
                     j = rxTPM::gs2t(0,S_cd,b,d);

                     (*this)(0,i,j) = 0.5 * ward;
                     (*this)(0,j,i) = (*this)(0,i,j);

                  }

               //then symmetrize rest
               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     i = rxTPM::gs2t(0,S_al,Tools::shift(a,b),Tools::par(b));
                     j = rxTPM::gs2t(0,S_ld,Tools::par(b),Tools::shift(d,b));

                     if(Tools::shift(a,b) > Tools::par(b))
                        phase_i = 1 - 2*S_al;
                     else
                        phase_i = 1;

                     if(Tools::par(b) > Tools::shift(d,b))
                        phase_j = 1 - 2*S_ld;
                     else
                        phase_j = 1;

                     (*this)(0,i,j) = 0.0;

                     for(int S_ab = 0;S_ab < 2;++S_ab)
                        for(int S_cd = 0;S_cd < 2;++S_cd){

                           (*this)(0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) 

                              * (1 - 2*S_ab) * (1 - 2*S_al) * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(0,S_ab,a,b,S_cd,b,d);

                        }

                     (*this)(0,j,i) = (*this)(0,i,j);

                  }

            }

         }
      }
   }

   //last regular: b == d --> W^l_{ab;cb}
   for(int a = 1;a < Tools::gL();++a){

      for(int c = a + 1;c < Tools::gL();++c){

         for(int b = c + 1;b < Tools::gL();++b){

            //first average out
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = (*this)(0,S_ab,a,b,S_cd,c,b);

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_cd) 

                           * (1 - 2*S_al) * (1 - 2*S_cl) * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cl,S_cd) * (*this)(b,0,S_al,a,0,S_cl,c,0);

                     }

                  i = rxTPM::gs2t(0,S_ab,a,b);
                  j = rxTPM::gs2t(0,S_cd,c,b);

                  (*this)(0,i,j) = 0.5*ward;
                  (*this)(0,j,i) = (*this)(0,i,j);

               }

            //then make the rest symmetric
            for(int S_al = 0;S_al < 2;++S_al)
               for(int S_cl = 0;S_cl < 2;++S_cl){

                  i = rxTPM::gs2t(0,S_al,Tools::shift(a,b),Tools::par(b));
                  j = rxTPM::gs2t(0,S_cl,Tools::shift(c,b),Tools::par(b));

                  if(Tools::shift(a,b) > Tools::par(b))
                     phase_i = 1 - 2*S_al;
                  else
                     phase_i = 1;

                  if(Tools::shift(c,b) > Tools::par(b))
                     phase_j = 1 - 2*S_cl;
                  else
                     phase_j = 1;

                  (*this)(0,i,j) = 0.0;

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     for(int S_cd = 0;S_cd < 2;++S_cd){

                        (*this)(0,i,j) += phase_i*phase_j * std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                           * (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al ) * (1 - 2*S_cl) 

                           * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cl,S_cd) * (*this)(0,S_ab,a,b,S_cd,c,b);

                     }

                  (*this)(0,j,i) = (*this)(0,i,j);

               }

         }
      }
   }

   //and at last: a == l
   for(int b = 1;b < Tools::gL();++b){

      for(int S_cd = 0;S_cd < 2;++S_cd){

         for(int c = 1;c < Tools::gL();++c){

            if(c != b){

               for(int d = c + S_cd;d < Tools::gL();++d){

                  if(d != b){

                     for(int S_ab = 0;S_ab < 2;++S_ab)
                        vec[S_ab] = (*this)(0,S_ab,0,b,S_cd,c,d);

                     for(int S_ab = 0;S_ab < 2;++S_ab){

                        double ward = vec[S_ab];

                        for(int S_lb = 0;S_lb < 2;++S_lb)
                           ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ab + 1.0) ) * Tools::g6j(0,0,S_lb,S_ab) * vec[S_lb];

                        i = rxTPM::gs2t(0,S_ab,0,b);
                        j = rxTPM::gs2t(0,S_cd,c,d);

                        (*this)(0,i,j) = 0.5 * ward;
                        (*this)(0,j,i) = (*this)(0,i,j);

                     }

                  }

               }

            }

         }

      }

   }

   //easy part: S = 3/2

   //first diagonal part: ab;ab
   for(int a = 1;a < Tools::gL();++a){

      for(int b = a + 1;b < Tools::gL();++b){

         i = rxTPM::gs2t(1,1,a,b);

         (*this)(1,i,i) = 1.0/3.0 * ( (*this)(1,i,i) + (*this)(a,1,1,0,b,1,0,b) + (*this)(b,1,1,a,0,1,a,0) );

         //rest is symmetric
         i = rxTPM::gs2t(1,1,Tools::par(a),Tools::shift(b,a));

         (*this)(1,i,i) = (*this)(0,1,1,a,b,1,a,b);

         i = rxTPM::gs2t(1,1,Tools::shift(a,b),Tools::par(b));

         (*this)(1,i,i) = (*this)(0,1,1,a,b,1,a,b);

      }
   }

   //three terms with one diagonal index:

   //1) b == d: ab;cb
   for(int a = 1;a < Tools::gL();++a){

      for(int c = a + 1;c < Tools::gL();++c){

         for(int b = c + 1;b < Tools::gL();++b){

            i = rxTPM::gs2t(1,1,a,b);
            j = rxTPM::gs2t(1,1,c,b);

            (*this)(1,i,j) = 0.5 * ( (*this)(1,i,j) + (*this)(b,1,1,a,0,1,c,0) );
            (*this)(1,j,i) = (*this)(1,i,j);

            i = rxTPM::gs2t(1,1,Tools::shift(a,b),Tools::par(b));
            j = rxTPM::gs2t(1,1,Tools::shift(c,b),Tools::par(b));

            if(Tools::shift(a,b) > Tools::par(b))
               phase_i = -1;
            else
               phase_i = 1;

            if(Tools::shift(c,b) > Tools::par(b))
               phase_j = -1;
            else
               phase_j = 1;

            (*this)(1,i,j) = phase_i*phase_j*(*this)(0,1,1,a,b,1,c,b);
            (*this)(1,j,i) = (*this)(1,i,j);

         }
      }
   }

   //2) b == c: ab;bc
   for(int a = 1;a < Tools::gL();++a){

      for(int b = a + 1;b < Tools::gL();++b){

         for(int c = b + 1;c < Tools::gL();++c){

            i = rxTPM::gs2t(1,1,a,b);
            j = rxTPM::gs2t(1,1,b,c);

            (*this)(1,i,j) = 0.5 * ( (*this)(1,i,j) + (*this)(b,1,1,a,0,1,0,c) );
            (*this)(1,j,i) = (*this)(1,i,j);

            i = rxTPM::gs2t(1,1,Tools::shift(a,b),Tools::par(b));
            j = rxTPM::gs2t(1,1,Tools::par(b),Tools::shift(c,b));

            if(Tools::shift(a,b) > Tools::par(b))
               phase_i = -1;
            else
               phase_i = 1;

            if(Tools::par(b) > Tools::shift(c,b))
               phase_j = -1;
            else
               phase_j = 1;

            (*this)(1,i,j) = phase_i*phase_j*(*this)(0,1,1,a,b,1,b,c);
            (*this)(1,j,i) = (*this)(1,i,j);

         }
      }
   }

   //3) a == c: ab;ac
   for(int a = 1;a < Tools::gL();++a){

      for(int b = a + 1;b < Tools::gL();++b){

         for(int c = b + 1;c < Tools::gL();++c){

            i = rxTPM::gs2t(1,1,a,b);
            j = rxTPM::gs2t(1,1,a,c);

            (*this)(1,i,j) = 0.5 * ( (*this)(1,i,j) + (*this)(a,1,1,0,b,1,0,c) );
            (*this)(1,j,i) = (*this)(1,i,j);

            i = rxTPM::gs2t(1,1,Tools::par(a),Tools::shift(b,a));
            j = rxTPM::gs2t(1,1,Tools::par(a),Tools::shift(c,a));

            if(Tools::par(a) > Tools::shift(b,a))
               phase_i = -1;
            else
               phase_i = 1;

            if(Tools::par(a) > Tools::shift(c,a))
               phase_j = -1;
            else
               phase_j = 1;

            (*this)(1,i,j) = phase_i * phase_j * (*this)(0,1,1,a,b,1,a,c);
            (*this)(1,j,i) = (*this)(1,i,j);

         }
      }
   }

}

/**
 * test if the projection is correct, number 1 -> print all
 */
void dDPM::test_proj_1() const {

   cout << "S = 1/2" << endl;
   cout << endl;

   for(int a = 0;a < Tools::gL();++a)
      for(int b = 0;b < Tools::gL();++b)
         for(int c = 0;c < Tools::gL();++c){

            cout << endl;

            //1) b = d
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "ab;cb\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,c,b) << "\t";

                  double ward = 0.0;

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_cl)

                           * (1 - 2*S_ab) * (1 - 2*S_cd) * Tools::g6j(0,0,S_ab,S_al) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(b,0,S_al,a,0,S_cl,c,0);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == b)
                           hard /= std::sqrt(2.0);

                        if(a == 0)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  cout << ward << endl;

               }

            //2) b = c
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "ab;bc\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,b,c) << "\t";

                  double ward = 0.0;

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * 

                           (1 - 2*S_ab) * Tools::g6j(0,0,S_ab,S_al) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(b,0,S_al,a,0,S_ld,0,c);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == b)
                           hard /= std::sqrt(2.0);

                        if(a == 0)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  cout << ward << endl;

               }

            //3) a = d
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "ab;ca\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,c,a) << "\t";

                  double ward = 0.0;

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl)

                           * (1 - 2*S_cd) * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(a,0,S_lb,0,b,S_cl,c,0);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == a)
                           hard /= std::sqrt(2.0);

                        if(0 == b)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  cout << ward << endl;

               }

            //4) a = c
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "ab;ac\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,a,c) << "\t";

                  double ward = 0.0;

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                           * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(a,0,S_lb,0,b,S_ld,0,c);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == a)
                           hard /= std::sqrt(2.0);

                        if(0 == b)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  cout << ward << endl;

               }

            cout << endl;

            //)5 a = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "lb;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,0,a,S_cd,b,c) << "\t";

                  double ward = 0.0;

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * Tools::g6j(0,0,S_lb,S_ab) * (*this)(0,0,S_lb,0,a,S_cd,b,c);

                  cout << ward << endl;

               }

            //6) b = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "al;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,0,S_cd,b,c) << "\t";

                  double ward = 0.0;

                  for(int S_al = 0;S_al < 2;++S_al){

                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al) 

                        * Tools::g6j(0,0,S_al,S_ab) * (*this)(0,0,S_al,a,0,S_cd,b,c);

                  }

                  cout << ward << endl;

               }

            //7) c = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "ab;ld\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,0,c) << "\t";

                  double ward = 0.0;

                  for(int S_ld = 0;S_ld < 2;++S_ld)
                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(0,0,S_ab,a,b,S_ld,0,c);

                  cout << ward << endl;

               }

            //8) d = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  cout << "ab;cl\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,c,0) << "\t";

                  double ward = 0.0;

                  for(int S_cl = 0;S_cl < 2;++S_cl)
                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(0,0,S_ab,a,b,S_cl,c,0);

                  cout << ward << endl;

               }

         }

   cout << endl;


   cout << endl;
   cout << "S = 3/2" << endl;
   cout << endl;

   for(int a = 0;a < Tools::gL();++a)
      for(int b = 0;b < Tools::gL();++b)
         for(int c = 0;c < Tools::gL();++c){

            cout << endl;
            cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,c,b) << "\t" << (*this)(b,1,1,a,0,1,c,0) << endl;
            cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,b,c) << "\t" << (*this)(b,1,1,a,0,1,0,c) << endl;
            cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,c,a) << "\t" << (*this)(a,1,1,0,b,1,c,0) << endl;
            cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,a,c) << "\t" << (*this)(a,1,1,0,b,1,0,c) << endl;
            cout << endl;

         }

}

/**
 * test if the projection is correct, number 2 -> print only the errors
 */
void dDPM::test_proj_2() const {

   cout << "S = 1/2" << endl;
   cout << endl;

   for(int a = 0;a < Tools::gL();++a)
      for(int b = 0;b < Tools::gL();++b)
         for(int c = 0;c < Tools::gL();++c){

            //1) b = d
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_cl)

                           * (1 - 2*S_ab) * (1 - 2*S_cd) * Tools::g6j(0,0,S_ab,S_al) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(b,0,S_al,a,0,S_cl,c,0);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == b)
                           hard /= std::sqrt(2.0);

                        if(a == 0)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  if(fabs(ward - (*this)(0,0,S_ab,a,b,S_cd,c,b)) > 1.0e-12)
                     cout << "ab;cb\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,c,b) << "\t" << ward << endl;

               }

            //2) b = c
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_al = 0;S_al < 2;++S_al)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        double hard = std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_al) * 

                           (1 - 2*S_ab) * Tools::g6j(0,0,S_ab,S_al) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(b,0,S_al,a,0,S_ld,0,c);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == b)
                           hard /= std::sqrt(2.0);

                        if(a == 0)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  if(fabs(ward - (*this)(0,0,S_ab,a,b,S_cd,b,c)) > 1.0e-12)
                     cout << "ab;bc\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,b,c) << "\t" << ward << endl;

               }

            //3) a = d
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_cl = 0;S_cl < 2;++S_cl){

                        double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_cl)

                           * (1 - 2*S_cd) * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(a,0,S_lb,0,b,S_cl,c,0);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == a)
                           hard /= std::sqrt(2.0);

                        if(0 == b)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  if(fabs(ward - (*this)(0,0,S_ab,a,b,S_cd,c,a)) > 1.0e-12)
                     cout << "ab;ca\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,c,a) << "\t" << ward << endl;

               }

            //4) a = c
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     for(int S_ld = 0;S_ld < 2;++S_ld){

                        double hard = std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                           * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(a,0,S_lb,0,b,S_ld,0,c);

                        if(a == b)
                           hard /= std::sqrt(2.0);

                        if(c == a)
                           hard /= std::sqrt(2.0);

                        if(0 == b)
                           hard *= std::sqrt(2.0);

                        if(c == 0)
                           hard *= std::sqrt(2.0);

                        ward += hard;

                     }

                  if(fabs(ward - (*this)(0,0,S_ab,a,b,S_cd,a,c))>1.0e-12)
                     cout << "ab;ac\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,a,c) << "\t" << ward << endl;

               }

            //)5 a = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_lb = 0;S_lb < 2;++S_lb)
                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_lb + 1.0) ) * Tools::g6j(0,0,S_lb,S_ab) * (*this)(0,0,S_lb,0,a,S_cd,b,c);

                  if(fabs(ward - (*this)(0,0,S_ab,0,a,S_cd,b,c)) > 1.0e-12)
                     cout << "lb;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,0,a,S_cd,b,c) << "\t" << ward << endl;


               }

            //6) b = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_al = 0;S_al < 2;++S_al){

                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al)

                        * Tools::g6j(0,0,S_al,S_ab) * (*this)(0,0,S_al,a,0,S_cd,b,c);

                  }

                  if(fabs(ward - (*this)(0,0,S_ab,a,0,S_cd,b,c)) > 1.0e-12)
                     cout << "al;cd\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,0,S_cd,b,c) << "\t" << ward << endl;

               }

            //7) c = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_ld = 0;S_ld < 2;++S_ld)
                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(0,0,S_ab,a,b,S_ld,0,c);

                  if(fabs(ward - (*this)(0,0,S_ab,a,b,S_cd,0,c)) > 1.0e-12)
                     cout << "ab;ld\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,0,c) << "\t" << ward << endl;

               }

            //8) d = l
            for(int S_ab = 0;S_ab < 2;++S_ab)
               for(int S_cd = 0;S_cd < 2;++S_cd){

                  double ward = 0.0;

                  for(int S_cl = 0;S_cl < 2;++S_cl){

                     ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) 

                        * Tools::g6j(0,0,S_cd,S_cl) * (*this)(0,0,S_ab,a,b,S_cl,c,0);

                  }


                  if(fabs(ward - (*this)(0,0,S_ab,a,b,S_cd,c,0)) > 1.0e-12)
                     cout << "ab;cl\t" << a << "\t" << b << "\t" << c << "\t(" << S_ab << ")\t(" << S_cd << ")\t" << (*this)(0,0,S_ab,a,b,S_cd,c,0) << "\t" << ward << endl;


               }

         }

   cout << endl;


   cout << endl;
   cout << "S = 3/2" << endl;
   cout << endl;

   for(int a = 0;a < Tools::gL();++a)
      for(int b = 0;b < Tools::gL();++b)
         for(int c = 0;c < Tools::gL();++c){

            if(fabs((*this)(0,1,1,a,b,1,c,b) - (*this)(b,1,1,a,0,1,c,0)) > 1.0e-12)
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,c,b) << "\t" << (*this)(b,1,1,a,0,1,c,0) << endl;

            if(fabs((*this)(0,1,1,a,b,1,b,c) - (*this)(b,1,1,a,0,1,0,c)) > 1.0e-12)
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,b,c) << "\t" << (*this)(b,1,1,a,0,1,0,c) << endl;

            if(fabs((*this)(0,1,1,a,b,1,c,a) - (*this)(a,1,1,0,b,1,c,0)) > 1.0e-12)
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,c,a) << "\t" << (*this)(a,1,1,0,b,1,c,0) << endl;

            if(fabs((*this)(0,1,1,a,b,1,a,c) - (*this)(a,1,1,0,b,1,0,c)) > 1.0e-12)
               cout << a << "\t" << b << "\t" << c << "\t|\t" << (*this)(0,1,1,a,b,1,a,c) << "\t" << (*this)(a,1,1,0,b,1,0,c) << endl;

         }

}

/**
 * Fill the dDPM object with the Hubbard Hamiltonian with on-site repulsion U
 * @param U the on-site repulsion
 */
void dDPM::hubbard(double U) {

   int a,b,c,d;
   int S_ab,S_cd;

   int sign;

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

            if(S_ab == S_cd){

               sign = 1 - 2*S_ab;

               //first hopping
               if( (a == c) && ( ( (b + 1)%Tools::gL() == d ) || ( b == (d + 1)%Tools::gL() ) ) )
                  (*this)(S,i,j) -= 1.0/( (Tools::gN() - 1.0) * (Tools::gN() - 2.0) );

               if( (b == c) && ( ( (a + 1)%Tools::gL() == d ) || ( a == (d + 1)%Tools::gL() ) ) )
                  (*this)(S,i,j) -= sign/( (Tools::gN() - 1.0) * (Tools::gN() - 2.0) );

               if( (a == d) && ( ( (b + 1)%Tools::gL() == c ) || ( b == (c + 1)%Tools::gL() ) ) )
                  (*this)(S,i,j) -= sign/( (Tools::gN() - 1.0) * (Tools::gN() - 2.0) );

               if( (b == d) && ( ( (a + 1)%Tools::gL() == c ) || ( a == (c + 1)%Tools::gL() ) ) )
                  (*this)(S,i,j) -= 1.0/( (Tools::gN() - 1.0) * (Tools::gN() - 2.0) );

               //only on-site interaction for singlet tp states:
               if(S_ab == 0)
                  if(i == j && a == b)
                     (*this)(S,i,j) += 2.0*U/(Tools::gN() - 2.0);

               if(a == b)
                  (*this)(S,i,j) /= std::sqrt(2.0);

               if(c == d)
                  (*this)(S,i,j) /= std::sqrt(2.0);

            }

         }

      }

   }

   this->symmetrize();

}

/**
 * initialize (*this) on the correctly normalized unitmatrix (so that the trace is N(N-1)(N-2)/2)
 */
void dDPM::unit(){

   int S_ab,S_cd;
   int a,b,c,d;

   double norm;

   for(int S = 0;S < 2;++S){

      for(int i = 0;i < gdim(S);++i){

         S_ab = gt2s(S,i,0);

         a = gt2s(S,i,1);
         b = gt2s(S,i,2);

         for(int j = i;j < gdim(S);++j){

            S_cd = gt2s(S,j,0);

            c = gt2s(S,j,1);
            d = gt2s(S,j,2);

            //set the norm
            norm = 1.0;

            if(a == b)
               norm /= std::sqrt(2.0);

            if(c == d)
               norm /= std::sqrt(2.0);

            (*this)(S,i,j) = 0.0;

            //set the unitmatrix
            if(a == c && b == d){

               if(S_ab == S_cd)
                  (*this)(S,i,j) += 1.0;

               if(a == 0)
                  (*this)(S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd);

               if(b == 0)
                  (*this)(S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd) * (1 - 2*S_ab) * (1 - 2*S_cd);

            }

            if(a == d && b == c){

               if(S_ab == S_cd)
                  (*this)(S,i,j) += (1 - 2*S_ab);

               if(a == 0)
                  (*this)(S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd) * (1 - 2*S_cd);

               if(b == 0)
                  (*this)(S,i,j) += std::sqrt( (2*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd) * (1 - 2*S_ab);

            }

            (*this)(S,i,j) *= norm * (Tools::gN()*(Tools::gN() - 1.0)*(Tools::gN() - 2.0)/(2*Tools::gL()*(2*Tools::gL() - 1.0)*(2*Tools::gL() - 2.0)));

            (*this)(S,j,i) = (*this)(S,i,j);

         }

      }

   }

}


/**
 * Construct the right hand side of the Newton equation for the determination of the search direction, 
 * the negative gradient of the potential:
 * @param t scaling factor of the potential
 * @param ham Hamiltonian of the current problem
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 */
void dDPM::constr_grad(double t,const dDPM &ham,const SUP &P){

   *this = P.gI1();

   dDPM hulp;

#ifdef __Q2_CON

   hulp.Q('D',P.gQ2());

   *this += hulp;

#endif

#ifdef __I2_CON

   hulp.I(P.gI2());

   *this += hulp;

#endif

#ifdef __Q1_CON

   hulp.Q(P.gQ1());

   *this += hulp;

#endif

#ifdef __G1_CON

   hulp.G1(P.gG1());

   *this += hulp;

#endif

#ifdef __G2_CON

   hulp.G2(P.gG2());

   *this += hulp;

#endif

   this->dscal(t);

   *this -= ham;

   this->proj();

}

/**
 * project *this on the traceless space but staying in the right symmetry area.
 */
void dDPM::proj_Tr(){

   double ward = this->ddot(Tools::gunit())/Tools::gunit().ddot(Tools::gunit());

   this->daxpy(-ward,Tools::gunit());

}


/**
 * total projection, both on traceless space and good third index symmetry
 */
void dDPM::proj(){

   this->proj_W();
   this->proj_Tr();

}

/**
 * solve the Newton equations for the determination of the search direction,
 * @param t scaling factor of the potential
 * @param P SUP matrix containing the inverse of the constraint matrices (carrier space matrices).
 * @param b right hand side (the gradient constructed int dDPM::constr_grad)
 * @return nr of iterations needed to converge to the desired accuracy
 */
int dDPM::solve(double t,const SUP &P,dDPM &b){

   int iter = 0;

   //delta = 0
   *this = 0;

   //residu:
   dDPM r(b);

   //norm van het residu
   double rr = r.ddot(r);

   //enkele variabelen
   double rr_old,ward;

   dDPM Hb;

   while(rr > 1.0e-7){ 

      Hb.H(t,b,P);

      ward = rr/b.ddot(Hb);

      //delta += ward*b
      this->daxpy(ward,b);

      //r -= ward*Hb
      r.daxpy(-ward,Hb);

      //nieuwe variabelen berekenen en oude overdragen
      rr_old = rr;
      rr = r.ddot(r);

      //nieuwe b nog:
      b.dscal(rr/rr_old);

      b += r;

      ++iter;

   }

   return iter;

}

/**
 * The hessian-map of the Newton system:
 * @param t potential scaling factor
 * @param b the dDPM on which the hamiltonian will work, the image will be put in (*this)
 * @param P the SUP matrix containing the constraints, (can be seen as the metric).
 */
void dDPM::H(double t,const dDPM &b,const SUP &P){

   this->L_map(P.gI1(),b);

   dDPM ddpm;

#ifdef __Q2_CON

   dDPM Q2b;
   Q2b.Q('U',b);

   dDPM hulp_dp;
   hulp_dp.L_map(P.gQ2(),Q2b);

   ddpm.Q('D',hulp_dp);

   *this += ddpm;

#endif

#ifdef __I2_CON
   
   dPPHM Ib;
   Ib.I(b);

   dPPHM hulp_pph;
   hulp_pph.L_map(P.gI2(),Ib);

   ddpm.I(hulp_pph);

   *this += ddpm;

#endif

#ifdef __Q1_CON
   
   dPPHM Q1b;
   Q1b.Q(b);

   hulp_pph.L_map(P.gQ1(),Q1b);

   ddpm.Q(hulp_pph);

   *this += ddpm;

#endif

#ifdef __G1_CON
   
   dPHHM G1b;
   G1b.G1(b);

   dPHHM hulp_phh;
   hulp_phh.L_map(P.gG1(),G1b);

   ddpm.G1(hulp_phh);

   *this += ddpm;

#endif

#ifdef __G2_CON
   
   dPHHM G2b;
   G2b.G2(b);

   hulp_phh.L_map(P.gG2(),G2b);

   ddpm.G2(hulp_phh);

   *this += ddpm;

#endif

   this->dscal(t);

   this->proj();

}

/**
 * perform a line search to determine what step size in along the Newton direction is ideal.
 * @param t potential scaling factor
 * @param P SUP matrix containing the inverse of the constraints (carrier space matrices)
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double dDPM::line_search(double t,SUP &P,const dDPM &ham){

   double tolerance = 1.0e-5*t;

   if(tolerance < 1.0e-12)
      tolerance = 1.0e-12;

   //neem de wortel uit P
   P.sqrt(1);

   //maak eerst een SUP van delta
   SUP S_delta;

   S_delta.fill(*this);

   //hulpje om dingskes in te steken:
   SUP hulp;
   hulp.L_map(P,S_delta);

   EIG eigen(hulp);

   double a = 0;

   double b = -1.0/eigen.min();

   double c(0);

   double ham_delta = ham.ddot(*this);

   while(b - a > tolerance){

      c = (b + a)/2.0;

      if( (ham_delta - t*eigen.lsfunc(c)) < 0.0)
         a = c;
      else
         b = c;

   }

   return c;

}

/**
 * perform a line search what step size in along a certain direction (*this) minimizes the potential, this one is used for extrapolation.
 * @param t potential scaling factor
 * @param W dDPM containing the current approximation of the W object
 * @param ham Hamiltonian of the problem
 * @return the steplength
 */
double dDPM::line_search(double t,const dDPM &W,const dDPM &ham){

   SUP P;

   P.fill(W);

   P.invert();

   return this->line_search(t,P,ham);

}

/**
 * The spincoupled Q2 map: maps a dDPM object onto itself
 * @param option if == 'U' up, if == 'D' down
 * @param ddpm_i input TPM
 */
void dDPM::Q(char option,const dDPM &ddpm_i){

   if(option == 'U'){

      TPM tpm;
      tpm.bar(1.0/(Tools::gN() - 2.0),ddpm_i);

      SPM spm;
      spm.bar(1.0/(Tools::gN() - 1.0),tpm);

      double ward = 2.0 * ddpm_i.trace()/( Tools::gN()*(Tools::gN() - 1.0)*(Tools::gN() - 2.0) );

      int a,b,c,d;
      int S_ab,S_cd;

      int sign_ab,sign_cd;

      double norm_ab,norm_cd;

      double hard;

      //start with the S = 1/2 block, this is the most difficult one:
      for(int i = 0;i < gdim(0);++i){

         S_ab = gt2s(0,i,0);

         a = gt2s(0,i,1);
         b = gt2s(0,i,2);

         sign_ab = 1 - 2*S_ab;

         norm_ab = 1.0;

         if(a == b)
            norm_ab /= std::sqrt(2.0);

         for(int j = i;j < gdim(0);++j){

            S_cd = gt2s(0,j,0);

            c = gt2s(0,j,1);
            d = gt2s(0,j,2);

            sign_cd = 1 - 2*S_cd;

            norm_cd = 1.0;

            if(c == d)
               norm_cd /= std::sqrt(2.0);

            hard = std::sqrt( (2*S_ab + 1.0) * (2*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd);

            //dp part
            (*this)(0,i,j) = -ddpm_i(0,i,j);

            //np(1)
            if(i == j)
               (*this)(0,i,j) += ward - spm(0,0);

            //terms that contribute when the spin is diagonal:
            if(S_ab == S_cd){

               //tp(1)
               (*this)(0,i,j) += tpm(S_ab,a,b,c,d);

               //sp(1) first term
               if(b == d)
                  (*this)(0,i,j) -= norm_ab * norm_cd * spm(a,c);

               //sp(2) first term
               if(a == d)
                  (*this)(0,i,j) -= sign_ab * norm_ab * norm_cd * spm(b,c);

               //sp(4) first term
               if(b == c)
                  (*this)(0,i,j) -= sign_cd * norm_ab * norm_cd * spm(a,d);

               //sp(5) first term
               if(a == c)
                  (*this)(0,i,j) -= norm_ab * norm_cd * spm(b,d);

            }

            if(b == 0){

               //tp(2)
               if(a == 0)
                  (*this)(0,i,j) += std::sqrt(2.0) * norm_ab * sign_ab * sign_cd * hard * tpm(S_cd,a,0,c,d);
               else
                  (*this)(0,i,j) += norm_ab * sign_ab * sign_cd * hard * tpm(S_cd,a,0,c,d);

               //sp(1) second term
               if(d == 0)
                  (*this)(0,i,j) -= sign_ab * sign_cd * norm_ab * norm_cd * hard * spm(a,c);

               //sp(3)
               if(a == d)
                  (*this)(0,i,j) -= sign_ab * norm_ab * norm_cd * hard * spm(c,0);

               //sp(4) second term
               if(0 == c)
                  (*this)(0,i,j) -= sign_ab * norm_ab * norm_cd * hard * spm(a,d);

               //sp(6)
               if(a == c)
                  (*this)(0,i,j) -= sign_ab * sign_cd * norm_ab * norm_cd * hard * spm(0,d);

               //np(4)
               if(c == 0 && a == d)
                  (*this)(0,i,j) += sign_ab * norm_ab * norm_cd * hard * ward;

               //np(6)
               if(d == 0 && a == c)
                  (*this)(0,i,j) += sign_ab * sign_cd * norm_ab * norm_cd * hard * ward;

            }

            if(a == 0){

               //tp(3)
               if(b == 0)
                  (*this)(0,i,j) += std::sqrt(2.0) * norm_ab * hard * tpm(S_cd,0,b,c,d);
               else
                  (*this)(0,i,j) += norm_ab * hard * tpm(S_cd,0,b,c,d);

               //sp(2) second term
               if(d == 0)
                  (*this)(0,i,j) -= sign_cd * norm_ab * norm_cd * hard * spm(b,c);

               //sp(5) second term
               if(c == 0)
                  (*this)(0,i,j) -= norm_ab * norm_cd * hard * spm(b,d);

               //sp(3) second part
               if(b == d)
                  (*this)(0,i,j) -= norm_ab * norm_cd * hard * spm(c,0);

               //sp(6) second part
               if(b == c)
                  (*this)(0,i,j) -= sign_cd * norm_ab * norm_cd * hard * spm(d,0);

               //np(3)
               if(c == 0 && b == d)
                  (*this)(0,i,j) += norm_ab * norm_cd * hard * ward;

               //np(5)
               if(d == 0 && b == c)
                  (*this)(0,i,j) += sign_cd * norm_ab * norm_cd * hard * ward;

            }

            if(0 == d){

               //tp(4)
               if(c == 0)
                  (*this)(0,i,j) += std::sqrt(2.0) * norm_cd * sign_ab * sign_cd * hard * tpm(S_ab,a,b,c,0);
               else
                  (*this)(0,i,j) += norm_cd * sign_ab * sign_cd * hard * tpm(S_ab,a,b,c,0);

               //sp(7) first term
               if(b == c)
                  (*this)(0,i,j) -= norm_ab * norm_cd * sign_cd * hard * spm(a,0);

               //sp(8) first term
               if(a == c)
                  (*this)(0,i,j) -= norm_ab * norm_cd * sign_ab * sign_cd * hard * spm(b,0);

            }

            if(b == d){

               //tp(5)
               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_cd) * tpm(Z,a,0,c,0);

               //correct for norms of the tpm
               if(a == 0)
                  hulp *= std::sqrt(2.0);

               if(c == 0)
                  hulp *= std::sqrt(2.0);

               (*this)(0,i,j) += norm_ab * norm_cd * sign_ab * sign_cd * std::sqrt( (2*S_ab + 1.0) * (2*S_cd + 1.0) ) * hulp;

               //sp(7) second term
               if(c == 0)
                  (*this)(0,i,j) -= norm_ab * norm_cd * hard * spm(a,0);

            }

            if(a == d){

               //tp(6)
               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_cd) * tpm(Z,b,0,c,0);

               if(b == 0)
                  hulp *= std::sqrt(2.0);

               if(c == 0)
                  hulp *= std::sqrt(2.0);

               (*this)(0,i,j) += sign_cd * std::sqrt( (2*S_ab + 1) * (2*S_cd + 1.0) ) * norm_ab * norm_cd * hulp;

               //sp(8) second term
               if(c == 0)
                  (*this)(0,i,j) -= sign_ab * norm_ab * norm_cd * hard * spm(b,0);

            }

            if(c == 0){

               //tp(7)
               if(d == 0)
                  (*this)(0,i,j) += std::sqrt(2.0) * norm_cd * hard * tpm(S_ab,a,b,0,d);
               else
                  (*this)(0,i,j) += norm_cd * hard * tpm(S_ab,a,b,0,d);

            }

            if(b == c){

               //tp(8)
               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_cd) * tpm(Z,a,0,d,0);

               if(a == 0)
                  hulp *= std::sqrt(2.0);

               if(d == 0)
                  hulp *= std::sqrt(2.0);

               (*this)(0,i,j) += sign_ab * std::sqrt( (2*S_ab + 1) * (2*S_cd + 1.0) ) * norm_ab * norm_cd * hulp;

            }

            if(a == c){

               //tp(8)
               double hulp = 0.0;

               //sum over intermediate spin
               for(int Z = 0;Z < 2;++Z)
                  hulp += (2*Z + 1.0) * Tools::g6j(0,0,Z,S_ab) * Tools::g6j(0,0,Z,S_cd) * tpm(Z,b,0,d,0);

               if(b == 0)
                  hulp *= std::sqrt(2.0);

               if(d == 0)
                  hulp *= std::sqrt(2.0);

               (*this)(0,i,j) += std::sqrt( (2*S_ab + 1) * (2*S_cd + 1.0) ) * norm_ab * norm_cd * hulp;

            }

         }
      }

      //then the S = 3/2 block, this should be easy, totally antisymmetrical 
      for(int i = 0;i < gdim(1);++i){

         a = gt2s(1,i,1);
         b = gt2s(1,i,2);

         for(int j = i;j < gdim(1);++j){

            c = gt2s(1,j,1);
            d = gt2s(1,j,2);

            (*this)(1,i,j) = tpm(1,a,b,c,d) - ddpm_i(1,i,j);

            if(i == j)
               (*this)(1,i,j) += ward - spm(0,0);

            if(b == d)
               (*this)(1,i,j) += tpm(1,a,0,c,0) - spm(a,c);

            if(b == c)
               (*this)(1,i,j) -= tpm(1,a,0,d,0) - spm(a,d);

            if(a == c)
               (*this)(1,i,j) += tpm(1,b,0,d,0) - spm(b,d);

         }
      }

   }
   else{

      TPM tpm;
      tpm.bar(1.0/(Tools::gN() - 2.0),ddpm_i);

      SPM spm;
      spm.bar(1.0/(Tools::gN() - 1.0),tpm);

      double dspm = ddpm_i.trace()/((Tools::gN() - 1.0)*(Tools::gN() - 2.0) * Tools::gL() );

      SPM breve;
      breve.breve(1.0/((Tools::gN() - 1.0)*(Tools::gN() - 2.0)),ddpm_i);

      ssdTPM ssdtpm;
      ssdtpm.bar(1.0/((Tools::gN() - 1.0)*(Tools::gN() - 2.0)),ddpm_i);

      PHM phm;
      phm.spinsum(1.0/(Tools::gN() - 2.0),ddpm_i);

      dTPM dtpm;
      dtpm.bar(1.0/(Tools::gN() - 2.0),ddpm_i);

      double ward = 2.0 * ddpm_i.dotunit()/( Tools::gN()*(Tools::gN() - 1.0)*(Tools::gN() - 2.0) );

      int a,b,c,d;
      int S_ab,S_cd;

      int sign_ab,sign_cd;

      double norm_ab,norm_cd;

      //start with the S = 1/2 block, this is the most difficult one:
      for(int i = 0;i < gdim(0);++i){

         S_ab = gt2s(0,i,0);

         a = gt2s(0,i,1);
         b = gt2s(0,i,2);

         sign_ab = 1 - 2*S_ab;

         norm_ab = 1.0;

         if(a == b)
            norm_ab /= std::sqrt(2.0);

         for(int j = i;j < gdim(0);++j){

            S_cd = gt2s(0,j,0);

            c = gt2s(0,j,1);
            d = gt2s(0,j,2);

            sign_cd = 1 - 2*S_cd;

            norm_cd = 1.0;

            if(c == d)
               norm_cd /= std::sqrt(2.0);

            //dp part
            (*this)(0,i,j) = -ddpm_i(0,i,j);

            if(i == j)
               (*this)(0,i,j) += ward - dspm;

            if(S_ab == S_cd){

               (*this)(0,i,j) += tpm(S_ab,a,b,c,d) + phm(S_ab,a,b,c,d) + (1 - 2*S_ab)*phm(S_ab,b,a,c,d) 
               
                  + (1 - 2*S_cd) * phm(S_cd,d,c,a,b) + phm(S_cd,c,d,a,b);

               if(b == d)
                  (*this)(0,i,j) -= norm_ab * norm_cd * ( spm(a,c) + breve(a,c) + ssdtpm(a,a,c) + ssdtpm(c,a,c) - dtpm(b,S_ab,a,c) );

               if(a == d)
                  (*this)(0,i,j) -= sign_ab * norm_ab * norm_cd * ( spm(b,c) + breve(b,c) + ssdtpm(b,b,c) + ssdtpm(c,b,c) - dtpm(a,S_ab,b,c) );

               if(b == c)
                  (*this)(0,i,j) -= sign_cd * norm_ab * norm_cd * ( spm(a,d) + breve(a,d) + ssdtpm(a,a,d) + ssdtpm(d,a,d) - dtpm(b,S_ab,a,d));

               if(a == c)
                  (*this)(0,i,j) -= norm_ab * norm_cd * ( spm(b,d) + breve(b,d) + ssdtpm(b,b,d) + ssdtpm(d,b,d) - dtpm(a,S_ab,b,d) );

            }

         }
      }

      //then the S = 3/2 block, this should be easy, totally antisymmetrical 
      for(int i = 0;i < gdim(1);++i){

         a = gt2s(1,i,1);
         b = gt2s(1,i,2);

         for(int j = i;j < gdim(1);++j){

            c = gt2s(1,j,1);
            d = gt2s(1,j,2);

            (*this)(1,i,j) = tpm(1,a,b,c,d) - ddpm_i(1,i,j) + phm(1,a,b,c,d) - phm(1,b,a,c,d) - phm(1,d,c,a,b) + phm(1,c,d,a,b);

            if(i == j)
               (*this)(1,i,j) += ward - dspm;

            if(b == d)
               (*this)(1,i,j) -= spm(a,c) + breve(a,c) + ssdtpm(a,a,c) + ssdtpm(c,a,c) - dtpm(b,1,a,c);

            if(b == c)
               (*this)(1,i,j) += spm(a,d) + breve(a,d) + ssdtpm(a,a,d) + ssdtpm(d,a,d) - dtpm(b,1,a,d);

            if(a == c)
               (*this)(1,i,j) -= spm(b,d) + breve(b,d) + ssdtpm(b,b,d) + ssdtpm(d,b,d) - dtpm(a,1,b,d);

         }
      }

   }

   this->symmetrize();

}

/**
 * @return the trace of *this with the unitmatrix in dDP space.
 */
double dDPM::dotunit() const {

   double ward = this->trace();

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

         ward += 2.0 * hard * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd);

      }

   return ward;

}

/**
 * map a dPPHM on a dDPM object with the I2 down map.
 * @param dpphm input dPPHM matrix
 */
void dDPM::I(const dPPHM &dpphm){

   int a,b,c,d;

   int S_ab,S_cd;

   TPM tpm;
   tpm.bar(1.0/(Tools::gN() - 2.0),dpphm);

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
               (*this)(S,i,j) += (2* (S_ + 0.5) + 1.0) * Tools::g6j(S,S_,S_ab,S_cd) * dpphm(S_,S_ab,a,b,S_cd,c,d);

            if(S_ab == S_cd)
               (*this)(S,i,j) += tpm(S_ab,a,b,c,d);

         }
      }
   }

   this->symmetrize();

}
