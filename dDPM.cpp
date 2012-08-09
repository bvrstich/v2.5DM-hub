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

   /*
   //next with two equalities:
   for(int a = 0;a < M;++a){

   if(a == l)
   ++a;

   if(a == M)
   break;

//a = d and  b = l: W^l_{al;ca}
for(int c = 0;c < a;++c){

if(c == l)
c++;

if(c == a)
break;

for(int S_ab = 0;S_ab < 2;++S_ab)
for(int S_cd = 0;S_cd < 2;++S_cd)
mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,c,a);

//first set the "l" block right
for(int S_ab = 0;S_ab < 2;++S_ab)
for(int S_cd = 0;S_cd < 2;++S_cd){

double ward = mat[S_ab][S_cd];

for(int S_al = 0;S_al < 2;++S_al)
ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * Tools::g6j(0,0,S_al,S_ab) * mat[S_al][S_cd];

for(int S_cl = 0;S_cl < 2;++S_cl){

ward += std::sqrt( 2.0 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) 

    * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,S_cl) * (*this)(a,0,0,l,l,S_cl,c,l);

    }

    i = rxTPM::gs2t(l,0,S_ab,a,l);
    j = rxTPM::gs2t(l,0,S_cd,c,a);

    if(a > l)
    phase_i = 1 - 2*S_ab;
    else
    phase_i = 1;

    (*this)[l](0,i,j) = phase_i * ward / 3.0;
    (*this)[l](0,j,i) = (*this)[l](0,i,j);

    }

//then the "a" block
for(int S_cl = 0;S_cl < 2;++S_cl){

i = rxTPM::gs2t(a,0,0,l,l);
j = rxTPM::gs2t(a,0,S_cl,c,l);

if(c > l)
phase_j = 1 - 2*S_cl;
else
phase_j = 1;

(*this)[a](0,i,j) = 0.0;

for(int S_ab = 0;S_ab < 2;++S_ab)
for(int S_cd = 0;S_cd < 2;++S_cd){

(*this)[a](0,i,j) += phase_j * std::sqrt( 0.5 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cl) * (1 - 2*S_cd)

* Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cl,S_cd) * (*this)(l,0,S_ab,a,l,S_cd,c,a);

}

(*this)[a](0,j,i) = (*this)[a](0,i,j);

}

}

//a = c and  b = l: W^l_{al;ad}
for(int d = a + 1;d < M;++d){

   if(d == l)
      d++;

   if(d == M)
      break;

   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd)
         mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,a,d);

   //first set the "l" block right
   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd){

         double ward = mat[S_ab][S_cd];

         for(int S_al = 0;S_al < 2;++S_al)
            ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * Tools::g6j(0,0,S_al,S_ab) * mat[S_al][S_cd];

         for(int S_ld = 0;S_ld < 2;++S_ld){

            ward += std::sqrt( 2.0 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) )

               * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(a,0,0,l,l,S_ld,l,d);

         }

         i = rxTPM::gs2t(l,0,S_ab,a,l);
         j = rxTPM::gs2t(l,0,S_cd,a,d);

         if(a > l)
            phase_i = 1 - 2*S_ab;
         else
            phase_i = 1;

         (*this)[l](0,i,j) = phase_i * ward / 3.0;
         (*this)[l](0,j,i) = (*this)[l](0,i,j);

      }

   //then the "a" block
   for(int S_ld = 0;S_ld < 2;++S_ld){

      i = rxTPM::gs2t(a,0,0,l,l);
      j = rxTPM::gs2t(a,0,S_ld,l,d);

      if(l > d)
         phase_j = 1 - 2*S_ld;
      else
         phase_j = 1;

      (*this)[a](0,i,j) = 0.0;

      for(int S_ab = 0;S_ab < 2;++S_ab)
         for(int S_cd = 0;S_cd < 2;++S_cd){

            (*this)[a](0,i,j) += phase_j * std::sqrt( 0.5 * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) 

               * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_ld,S_cd) * (*this)(l,0,S_ab,a,l,S_cd,a,d);

         }

      (*this)[a](0,j,i) = (*this)[a](0,i,j);

   }

}

}

//another one with 2 equalities: a = c = l
for(int b = l + 1;b < M;++b)
for(int d = b + 1;d < M;++d){

   //save the numbers first in the matrix mat
   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd)
         mat[S_ab][S_cd] = (*this)(l,0,S_ab,l,b,S_cd,l,d);

   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd){

         i = rxTPM::gs2t(l,0,S_ab,l,b);
         j = rxTPM::gs2t(l,0,S_cd,l,d);

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

         (*this)[l](0,i,j) = 0.25 * ward;
         (*this)[l](0,j,i) = (*this)[l](0,i,j);

      }

}

//next one with 2 equalities: b = c = l
for(int a = 0;a < l;++a)
for(int d = l + 1;d < M;++d){

   //save the numbers first in the matrix mat
   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd)
         mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,l,d);

   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd){

         i = rxTPM::gs2t(l,0,S_ab,a,l);
         j = rxTPM::gs2t(l,0,S_cd,l,d);

         double ward = mat[S_ab][S_cd];

         for(int S_al = 0;S_al < 2;++S_al)
            ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * Tools::g6j(0,0,S_ab,S_al) * mat[S_al][S_cd];

         for(int S_ld = 0;S_ld < 2;++S_ld)
            ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_ld + 1.0) ) * Tools::g6j(0,0,S_cd,S_ld) * mat[S_ab][S_ld];

         for(int S_al = 0;S_al < 2;++S_al)
            for(int S_ld = 0;S_ld < 2;++S_ld){

               ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab)

                  * Tools::g6j(0,0,S_ab,S_al) * Tools::g6j(0,0,S_cd,S_ld) * mat[S_al][S_ld];

            }

         (*this)[l](0,i,j) = 0.25*ward;
         (*this)[l](0,j,i) = (*this)[l](0,i,j);

      }

}

//next one with 2 equalities: b = d = l
for(int a = 0;a < l;++a)
for(int c = a + 1;c < l;++c){

   //save the numbers first in the matrix mat
   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd)
         mat[S_ab][S_cd] = (*this)(l,0,S_ab,a,l,S_cd,c,l);

   for(int S_ab = 0;S_ab < 2;++S_ab)
      for(int S_cd = 0;S_cd < 2;++S_cd){

         i = rxTPM::gs2t(l,0,S_ab,a,l);
         j = rxTPM::gs2t(l,0,S_cd,c,l);

         double ward = mat[S_ab][S_cd];

         for(int S_al = 0;S_al < 2;++S_al)
            ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_al + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * Tools::g6j(0,0,S_ab,S_al) * mat[S_al][S_cd];

         for(int S_cl = 0;S_cl < 2;++S_cl)
            ward += std::sqrt( (2.0*S_cd + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_cd) * (1 - 2*S_cl) * Tools::g6j(0,0,S_cd,S_cl) * mat[S_ab][S_cl];

         for(int S_al = 0;S_al < 2;++S_al)
            for(int S_cl = 0;S_cl < 2;++S_cl){

               ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab)

                  * (1 - 2*S_cd) * (1 - 2*S_cl) * Tools::g6j(0,0,S_ab,S_al) * Tools::g6j(0,0,S_cd,S_cl) * mat[S_al][S_cl];

            }

         (*this)[l](0,i,j) = 0.25*ward;
         (*this)[l](0,j,i) = (*this)[l](0,i,j);

      }

}

//one equality: start with a == c
for(int a = 0;a < M;++a){

   if(a == l)
      a++;

   if(a == M)
      break;

   for(int b = a + 1;b < M;++b){

      if(b == l)
         b++;

      if(b == M)
         break;

      for(int d = b + 1;d < M;++d){

         if(d == l)
            d++;

         if(d == M)
            break;

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               //first take the average
               double ward = (*this)(l,0,S_ab,a,b,S_cd,a,d);

               for(int S_lb = 0;S_lb < 2;++S_lb)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) ) 

                        * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(a,0,S_lb,l,b,S_ld,l,d);

                  }

               i = rxTPM::gs2t(l,0,S_ab,a,b);
               j = rxTPM::gs2t(l,0,S_cd,a,d);

               (*this)[l](0,i,j) = 0.5*ward;
               (*this)[l](0,j,i) = (*this)[l](0,i,j);

            }

         //then symmetrize the other term
         for(int S_lb = 0;S_lb < 2;++S_lb)
            for(int S_ld = 0;S_ld < 2;++S_ld){

               i = rxTPM::gs2t(a,0,S_lb,l,b);
               j = rxTPM::gs2t(a,0,S_ld,l,d);

               if(l > b)
                  phase_i = 1 - 2*S_lb;
               else
                  phase_i = 1;

               if(l > d)
                  phase_j = 1 - 2*S_ld;
               else
                  phase_j = 1;

               (*this)[a](0,i,j) = 0.0;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[a](0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_lb + 1.0) * (2.0*S_ld + 1.0) )

                        * Tools::g6j(0,0,S_ab,S_lb) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(l,0,S_ab,a,b,S_cd,a,d);

                  }

               (*this)[a](0,j,i) = (*this)[a](0,i,j);

            }

      }
   }
}


//then b == c: ab;bd
for(int a = 0;a < M;++a){

   if(a == l)
      a++;

   if(a == M)
      break;

   for(int b = a + 1;b < M;++b){

      if(b == l)
         b++;

      if(b == M)
         break;

      for(int d = b + 1;d < M;++d){

         if(d == l)
            d++;

         if(d == M)
            break;

         //first take average
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               double ward = (*this)(l,0,S_ab,a,b,S_cd,b,d);

               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_ld = 0;S_ld < 2;++S_ld){

                     ward += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_al)

                        * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(b,0,S_al,a,l,S_ld,l,d);

                  }

               i = rxTPM::gs2t(l,0,S_ab,a,b);
               j = rxTPM::gs2t(l,0,S_cd,b,d);

               (*this)[l](0,i,j) = 0.5 * ward;
               (*this)[l](0,j,i) = (*this)[l](0,i,j);

            }

         //then symmetrize rest
         for(int S_al = 0;S_al < 2;++S_al)
            for(int S_ld = 0;S_ld < 2;++S_ld){

               i = rxTPM::gs2t(b,0,S_al,a,l);
               j = rxTPM::gs2t(b,0,S_ld,l,d);

               if(a > l)
                  phase_i = 1 - 2*S_al;
               else
                  phase_i = 1;

               if(l > d)
                  phase_j = 1 - 2*S_ld;
               else
                  phase_j = 1;

               (*this)[b](0,i,j) = 0.0;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[b](0,i,j) += phase_i * phase_j * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) * (2.0*S_al + 1.0) * (2.0*S_ld + 1.0) ) 

                        * (1 - 2*S_ab) * (1 - 2*S_al) * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cd,S_ld) * (*this)(l,0,S_ab,a,b,S_cd,b,d);

                  }

               (*this)[b](0,j,i) = (*this)[b](0,i,j);

            }

      }
   }
}

//last regular: b == d --> W^l_{ab;cb}
for(int a = 0;a < M;++a){

   if(a == l)
      a++;

   if(a == M)
      break;

   for(int c = a + 1;c < M;++c){

      if(c == l)
         c++;

      if(c == M)
         break;

      for(int b = c + 1;b < M;++b){

         if(b == l)
            b++;

         if(b == M)
            break;

         //first average out
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               double ward = (*this)(l,0,S_ab,a,b,S_cd,c,b);

               for(int S_al = 0;S_al < 2;++S_al)
                  for(int S_cl = 0;S_cl < 2;++S_cl){

                     ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * (1 - 2*S_ab) * (1 - 2*S_cd) 

                        * (1 - 2*S_al) * (1 - 2*S_cl) * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cl,S_cd) * (*this)(b,0,S_al,a,l,S_cl,c,l);

                  }

               i = rxTPM::gs2t(l,0,S_ab,a,b);
               j = rxTPM::gs2t(l,0,S_cd,c,b);

               (*this)[l](0,i,j) = 0.5*ward;
               (*this)[l](0,j,i) = (*this)[l](0,i,j);

            }

         //then make the rest symmetric
         for(int S_al = 0;S_al < 2;++S_al)
            for(int S_cl = 0;S_cl < 2;++S_cl){

               i = rxTPM::gs2t(b,0,S_al,a,l);
               j = rxTPM::gs2t(b,0,S_cl,c,l);

               if(a > l)
                  phase_i = 1 - 2*S_al;
               else
                  phase_i = 1;

               if(c > l)
                  phase_j = 1 - 2*S_cl;
               else
                  phase_j = 1;

               (*this)[b](0,i,j) = 0.0;

               for(int S_ab = 0;S_ab < 2;++S_ab)
                  for(int S_cd = 0;S_cd < 2;++S_cd){

                     (*this)[b](0,i,j) += phase_i*phase_j * std::sqrt( (2.0*S_al + 1.0) * (2.0*S_cl + 1.0) * (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) 

                        * (1 - 2*S_ab) * (1 - 2*S_cd) * (1 - 2*S_al ) * (1 - 2*S_cl) 

                        * Tools::g6j(0,0,S_al,S_ab) * Tools::g6j(0,0,S_cl,S_cd) * (*this)(l,0,S_ab,a,b,S_cd,c,b);

                  }

               (*this)[b](0,j,i) = (*this)[b](0,i,j);

            }

      }
   }
}

//two more terms: b == l
for(int a = 0;a < l;++a){

   for(int S_cd = 0;S_cd < 2;++S_cd){

      for(int c = 0;c < M;++c){

         if(c != l && c != a){

            for(int d = c + S_cd;d < M;++d){

               if(d != l && d != a){

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     vec[S_ab] = (*this)(l,0,S_ab,a,l,S_cd,c,d);

                  for(int S_ab = 0;S_ab < 2;++S_ab){

                     double ward = vec[S_ab];

                     for(int S_al = 0;S_al < 2;++S_al)
                        ward += std::sqrt( (2.0*S_al + 1.0) * (2.0*S_ab + 1.0) ) * (1 - 2*S_al) * (1 - 2*S_ab) * Tools::g6j(0,0,S_al,S_ab) * vec[S_al];

                     i = rxTPM::gs2t(l,0,S_ab,a,l);
                     j = rxTPM::gs2t(l,0,S_cd,c,d);

                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

               }

            }

         }

      }

   }

}

//and at last: a == l
for(int b = l + 1;b < M;++b){

   for(int S_cd = 0;S_cd < 2;++S_cd){

      for(int c = 0;c < M;++c){

         if(c != l && c != b){

            for(int d = c + S_cd;d < M;++d){

               if(d != l && d != b){

                  for(int S_ab = 0;S_ab < 2;++S_ab)
                     vec[S_ab] = (*this)(l,0,S_ab,l,b,S_cd,c,d);

                  for(int S_ab = 0;S_ab < 2;++S_ab){

                     double ward = vec[S_ab];

                     for(int S_lb = 0;S_lb < 2;++S_lb)
                        ward += std::sqrt( (2.0*S_lb + 1.0) * (2.0*S_ab + 1.0) ) * Tools::g6j(0,0,S_lb,S_ab) * vec[S_lb];

                     i = rxTPM::gs2t(l,0,S_ab,l,b);
                     j = rxTPM::gs2t(l,0,S_cd,c,d);

                     (*this)[l](0,i,j) = 0.5 * ward;
                     (*this)[l](0,j,i) = (*this)[l](0,i,j);

                  }

               }

            }

         }

      }

   }

}

//easy part: S = 3/2

//first diagonal part: ab;ab
for(int a = 0;a < M;++a){

   if(a == l)
      ++a;

   if(a == M)
      break;

   for(int b = a + 1;b < M;++b){

      if(b == l)
         ++b;

      if(b == M)
         break;

      i = rxTPM::gs2t(l,1,1,a,b);

      (*this)[l](1,i,i) = 1.0/3.0 * ( (*this)[l](1,i,i) + (*this)(a,1,1,l,b,1,l,b) + (*this)(b,1,1,a,l,1,a,l) );

      //rest is symmetric
      i = rxTPM::gs2t(a,1,1,l,b);

      (*this)[a](1,i,i) = (*this)(l,1,1,a,b,1,a,b);

      i = rxTPM::gs2t(b,1,1,a,l);

      (*this)[b](1,i,i) = (*this)(l,1,1,a,b,1,a,b);

   }
}

//three terms with one diagonal index:

//1) b == d: ab;cb
for(int a = 0;a < M;++a){

   if(a == l)
      a++;

   if(a == M)
      break;

   for(int c = a + 1;c < M;++c){

      if(c == l)
         c++;

      if(c == M)
         break;

      for(int b = c + 1;b < M;++b){

         if(b == l)
            b++;

         if(b == M)
            break;

         i = rxTPM::gs2t(l,1,1,a,b);
         j = rxTPM::gs2t(l,1,1,c,b);

         (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(b,1,1,a,l,1,c,l) );
         (*this)[l](1,j,i) = (*this)[l](1,i,j);

         i = rxTPM::gs2t(b,1,1,a,l);
         j = rxTPM::gs2t(b,1,1,c,l);

         if(a > l)
            phase_i = -1;
         else
            phase_i = 1;

         if(c > l)
            phase_j = -1;
         else
            phase_j = 1;

         (*this)[b](1,i,j) = phase_i*phase_j*(*this)(l,1,1,a,b,1,c,b);
         (*this)[b](1,j,i) = (*this)[b](1,i,j);

      }
   }
}

//2) b == c: ab;bc
for(int a = 0;a < M;++a){

   if(a == l)
      a++;

   if(a == M)
      break;

   for(int b = a + 1;b < M;++b){

      if(b == l)
         b++;

      if(b == M)
         break;

      for(int c = b + 1;c < M;++c){

         if(c == l)
            c++;

         if(c == M)
            break;

         i = rxTPM::gs2t(l,1,1,a,b);
         j = rxTPM::gs2t(l,1,1,b,c);

         (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(b,1,1,a,l,1,l,c) );
         (*this)[l](1,j,i) = (*this)[l](1,i,j);

         i = rxTPM::gs2t(b,1,1,a,l);
         j = rxTPM::gs2t(b,1,1,l,c);

         if(a > l)
            phase_i = -1;
         else
            phase_i = 1;

         if(l > c)
            phase_j = -1;
         else
            phase_j = 1;

         (*this)[b](1,i,j) = phase_i*phase_j*(*this)(l,1,1,a,b,1,b,c);
         (*this)[b](1,j,i) = (*this)[b](1,i,j);

      }
   }
}

//3) a == c: ab;ac
for(int a = 0;a < M;++a){

   if(a == l)
      a++;

   if(a == M)
      break;

   for(int b = a + 1;b < M;++b){

      if(b == l)
         b++;

      if(b == M)
         break;

      for(int c = b + 1;c < M;++c){

         if(c == l)
            c++;

         if(c == M)
            break;

         i = rxTPM::gs2t(l,1,1,a,b);
         j = rxTPM::gs2t(l,1,1,a,c);

         (*this)[l](1,i,j) = 0.5 * ( (*this)[l](1,i,j) + (*this)(a,1,1,l,b,1,l,c) );
         (*this)[l](1,j,i) = (*this)[l](1,i,j);

         i = rxTPM::gs2t(a,1,1,l,b);
         j = rxTPM::gs2t(a,1,1,l,c);

         if(l > b)
            phase_i = -1;
         else
            phase_i = 1;

         if(l > c)
            phase_j = -1;
         else
            phase_j = 1;

         (*this)[a](1,i,j) = phase_i * phase_j * (*this)(l,1,1,a,b,1,a,c);
         (*this)[a](1,j,i) = (*this)[a](1,i,j);

      }
   }
}
*/
}
