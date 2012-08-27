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
ssdTPM::ssdTPM() : SPM() { }

/**
 * copy constructor 
 * @param ssdtpm object that will be copied into this.
 */
ssdTPM::ssdTPM(const ssdTPM &ssdtpm) : SPM(ssdtpm) { }

/**
 * destructor
 */
ssdTPM::~ssdTPM(){ }

/**
 * @return the trace of the ssdTPM object
 */
double ssdTPM::trace() const {

   return Tools::gL() * SPM::trace();

}


/**
 * @return inproduct of (*this) ssdtpm with ssdtpm_i, defined as Tr (A B)
 * @param ssdtpm_i input ssdtpm
 */
double ssdTPM::ddot(const ssdTPM &ssdtpm_i) const{

   return Tools::gL() * SPM::ddot(ssdtpm_i);

}

/**
 * access the numbers in sp mode
 * @param l block index
 * @param a row spatial index 
 * @param b column spatial index 
 */
double ssdTPM::operator()(int l,int a,int b) const {

   return (*this)(Tools::shift(a,l),Tools::shift(b,l));

}

/**
 * map a dDPM matrix on a ssdTPM by using a special "bar" function, the different blocks of a dDPM matrix.
 * special because the intermediate spin is not required to be diagonal, instead the usual [][]{...} is used.
 * @param scale the ssdTPM with this number
 * @param ddpm input dDPM object
 */
void ssdTPM::bar(double scale,const dDPM &ddpm){

   double ward,hard;

   for(int a = 0;a < Tools::gL();++a)
      for(int c = a;c < Tools::gL();++c){

         (*this)(a,c) = 0.0;

         //first S = 1/2
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int b = 0;b < Tools::gL();++b){

                  hard = ddpm(0,0,S_ab,a,b,S_cd,c,b);

                  if(a == b)
                     hard *= std::sqrt(2.0);

                  if(c == b)
                     hard *= std::sqrt(2.0);

                  ward += hard;

               }

               (*this)(a,c) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,S_cd) * ward;

            }

         //then S = 3/2
         for(int b = 0;b < Tools::gL();++b)
            (*this)(a,c) -= 2.0 * ddpm(0,1,1,a,b,1,c,b);


         //finally scale
         (*this)(a,c) *= scale;

      }

   this->symmetrize();

}


/**
 * map a dPPHM matrix on a ssdTPM by using a special "bar" function, the different blocks of a dPPHM matrix.
 * special because the intermediate spin is not required to be diagonal, instead the usual [][] is used.
 * @param scale the ssdTPM with this number
 * @param dpphm input dPPHM object
 */
void ssdTPM::bar(double scale,const dPPHM &dpphm){

   double ward,hard;

   for(int a = 0;a < Tools::gL();++a)
      for(int c = a;c < Tools::gL();++c){

         (*this)(a,c) = 0.0;

         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int b = 0;b < Tools::gL();++b){

                  hard = dpphm(0,S_ab,a,b,S_cd,c,b);

                  if(a == b)
                     hard *= std::sqrt(2.0);

                  if(c == b)
                     hard *= std::sqrt(2.0);

                  ward += hard;

               }

               (*this)(a,c) += std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * ward;

            }

         //finally scale
         (*this)(a,c) *= scale;

      }

   this->symmetrize();

}

/**
 * map a dPHHM matrix on a ssdTPM by using a spinsummed "skew-bar" function on the different blocks of a dPHHM matrix.
 * spinsummed because the intermediate spin is not required to be diagonal, instead the usual [][] is used.
 * BEWARE: matrix not symmetric in general!
 * @param scale the ssdTPM with this number
 * @param dphhm input dPHHM object
 */
void ssdTPM::skew_bar(double scale,const dPHHM &dphhm){

   double ward;

   for(int a = 0;a < Tools::gL();++a)
      for(int c = 0;c < Tools::gL();++c){

         (*this)(a,c) = 0.0;

         for(int S_bl = 0;S_bl < 2;++S_bl)
            for(int S_dl = 0;S_dl < 2;++S_dl){

               ward = 0.0;

               for(int b = 0;b < Tools::gL();++b)
                  ward += dphhm(0,S_bl,b,b,S_dl,a,c);

               (*this)(a,c) += (1 - 2*S_bl) * std::sqrt( (2.0*S_bl + 1.0) * (2.0*S_dl + 1.0) ) * ward;

            }

         //finally scale
         (*this)(a,c) *= scale;

      }

}
