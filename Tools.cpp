#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_sf_coupling.h>

#include "include.h"

using std::cout;
using std::endl;

dDPM *Tools::unit;

double *Tools::x6j;
double *Tools::x9j;

int Tools::L;
int Tools::N;


/**
 * initialize and allocate the static variables
 * @param L_in nr of sites
 * @param N_in nr of particles
 */
void Tools::init(int L_in,int N_in){

   L = L_in;
   N = N_in;

   x6j = new double [16];

   x9j = new double [16];

   //fill the x6j list
   for(int S_ = 0;S_ < 2;++S_)
      for(int S = 0;S < 2;++S)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x6j[8*S_ + 4*S + 2*S_ab + S_cd] = gsl_sf_coupling_6j(2*S_ + 1,1,2*S_ab,2*S + 1,1,2*S_cd);

   //fill the x9j list
   for(int S = 0;S < 2;++S)
      for(int Z = 0;Z < 2;++Z)
         for(int S_ab = 0;S_ab < 2;++S_ab)
            for(int S_cd = 0;S_cd < 2;++S_cd)
               x9j[8*S + 4*Z + 2*S_ab + S_cd] = gsl_sf_coupling_9j(1,1,2*Z,2*S + 1,2*S_ab,1,2*S_cd,1,1);

}

/**
 * This needs a separate function because I'm stupid
 */
void Tools::sunit(){

   unit = new dDPM();

   *unit = 0.0;

   for(int S = 0;S < 2;++S)
      for(int i = 0;i < unit->gdim(S);++i)
         (*unit)(S,i,i) = 1.0;

   unit->proj_W();

}

/**
 * deallocate the static lists and objects
 */
void Tools::clear(){

   delete [] x6j;

   delete [] x9j;

   delete unit;

}

/**
 * store only the six j symbols needed
 * @param S_ 0 or 1 means 1/2 or 3/2
 * @param S 0 or 1 means 1/2 or 3/2
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g6j(int S_,int S,int S_ab,int S_cd){

   return x6j[8*S_ + 4*S + 2*S_ab + S_cd];

}

/**
 * store only the nine j symbols needed
 * @param S 0 or 1 means 1/2 or 3/2
 * @param Z 0 or 1 means 0 or 1
 * @param S_ab 0 or 1 means 0 or 1
 * @param S_cd 0 or 1 means 0 or 1
 * @return the 6j symbol
 */
double Tools::g9j(int S,int Z,int S_ab,int S_cd){

   return x9j[8*S + 4*Z + 2*S_ab + S_cd];

}

/**
 * @return the nr of sites
 */
int Tools::gL(){

   return L;

}

/**
 * @return the nr of particles
 */
int Tools::gN(){

   return N;

}

/**
 * @param index a
 * @param index b
 * @return the index a shifted by b, used for translational invariance
 */
int Tools::shift(int a,int b){

   return (a - b + L)%L;

}

/**
 * @param index a
 * @return the parity inverse of the index a
 */
int Tools::par(int a){

   return (L - a)%L;

}

/**
 * @return the projected unit matrix
 */
dDPM &Tools::gunit(){

   return *unit;

}
