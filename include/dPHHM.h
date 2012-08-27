#ifndef dPHHM_H
#define dPHHM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "rxPHM.h"

class SUP;

/**
 * @author Brecht Verstichel
 * @date 27-08-2012\n\n
 * This class dPHHM is a class written for the G1 and G2 conditions, which is a PHHM matrix with the spatial third index diagonal.
 * The v2.5DM program will use this as the central variable. Translational invariance is exploited, which means only one block is
 * to be stored. It enherits from the rxPHM class and overloads the functions ddot and trace.
 */
class dPHHM : public rxPHM {

   public:

      //constructor
      dPHHM();

      //copy constructor
      dPHHM(const dPHHM &);

      //destructor
      virtual ~dPHHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      double operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      double trace() const;

      double ddot(const dPHHM &) const;

      void G1(const dDPM &);

      double skew_trace() const;

      void G2(const dDPM &);

   private:

};

#endif
