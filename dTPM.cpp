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
dTPM::dTPM() : xSPM() { }

/**
 * copy constructor 
 * @param dtpm object that will be copied into this.
 */
dTPM::dTPM(const dTPM &dtpm) : xSPM(dtpm) { }

/**
 * destructor
 */
dTPM::~dTPM(){ }

/**
 * @return the trace of the dTPM object
 */
double dTPM::trace() const {

   return Tools::gL() * xSPM::trace();

}

/**
 * @return inproduct of (*this) dtpm with dtpm_i, defined as Tr (A B)
 * @param dtpm_i input dtpm
 */
double dTPM::ddot(const dTPM &dtpm_i) const{

   return Tools::gL() * xSPM::ddot(dtpm_i);

}

/**
 * access the numbers in sp mode
 * @param l block index
 * @param S spin index
 * @param a row spatial index 
 * @param b column spatial index 
 */
double dTPM::operator()(int l,int S,int a,int b) const {

   return (*this)(S,Tools::shift(a,l),Tools::shift(b,l));

}

/**
 * special "bar" function that maps a dDPM on a dTPM object, see notes for info.
 * @param scale the factor you scale the dTPM with
 * @param ddpm input dDPM object
 */
void dTPM::bar(double scale,const dDPM &ddpm){

   double ward,hard;

   for(int b = 0;b < Tools::gL();++b)
      for(int d = b;d < Tools::gL();++d){

         //first Z = 0
         (*this)(0,b,d) = 0.0;

         //only S = 1/2 contribution
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int a = 0;a < Tools::gL();++a){

                  hard = ddpm(0,0,S_ab,a,b,S_cd,a,d);

                  if(a == b)
                     hard *= std::sqrt(2.0);

                  if(a == d)
                     hard *= std::sqrt(2.0);

                  ward += hard;

               }

               (*this)(0,b,d) += 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,0) * Tools::g6j(0,0,S_cd,0) * ward;

            }

         (*this)(0,b,d) *= scale;

         //then Z = 1
         (*this)(1,b,d) = 0.0;

         //first S = 1/2 contribution
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd){

               ward = 0.0;

               for(int a = 0;a < Tools::gL();++a){

                  hard = ddpm(0,0,S_ab,a,b,S_cd,a,d);

                  if(a == b)
                     hard *= std::sqrt(2.0);

                  if(a == d)
                     hard *= std::sqrt(2.0);

                  ward += hard;

               }

               (*this)(1,b,d) += 2.0 * std::sqrt( (2.0*S_ab + 1.0) * (2.0*S_cd + 1.0) ) * Tools::g6j(0,0,S_ab,1) * Tools::g6j(0,0,S_cd,1) * ward;

            }

         //finally S = 3/2 contribution
         for(int a = 0;a < Tools::gL();++a)
            (*this)(1,b,d) += 4.0/3.0 * ddpm(0,1,1,a,b,1,a,d);

         (*this)(1,b,d) *= scale;

      }

   this->symmetrize();

}
