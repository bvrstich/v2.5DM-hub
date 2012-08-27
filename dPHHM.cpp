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
