#ifndef dPPHM_H
#define dPPHM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "xTPM.h"

class SUP;

/**
 * @author Brecht Verstichel
 * @date 08-08-2012\n\n
 * This class dPPHM is a class written for the I2 and Q1 conditions, which is a PPHM matrix with the spatial third index diagonal.
 * The v2.5DM program will use this as the central variable. Translational invariance is exploited, which means only one block is
 * to be stored. It enherits from the xTPM class and overloads the functions ddot and trace.
 */
class dPPHM : public xTPM {


   public:

      //constructor
      dPPHM();

      //copy constructor
      dPPHM(const dPPHM &);

      //destructor
      virtual ~dPPHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      double operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      static int get_inco(int S,int S_ab,int a,int b);

      double trace() const;

      double ddot(const dPPHM &) const;

      void I(const dDPM &);

      void Q(const dDPM &);

      double barbreve() const;

   private:

};

#endif
