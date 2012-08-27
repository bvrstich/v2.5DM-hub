#ifndef dDPM_H
#define dDPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "rxTPM.h"

class SUP;
class dPPHM;

/**
 * @author Brecht Verstichel
 * @date 08-08-2012\n\n
 * This class dDPM is a class written for the I1 and Q2 conditions, which is a DPM matrix with the spacial third index diagonal.
 * The v2.5DM program will use this as the central variable. Translational invariance is exploited, which means only one block is
 * to be stored. It enherits from the rxTPM class and overloads the functions ddot and trace.
 */
class dDPM : public rxTPM {

   public:

      //constructor
      dDPM();

      //copy constructor
      dDPM(const dDPM &);

      //destructor
      virtual ~dDPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      double operator()(int l,int S,int S_ab,int a,int b,int S_cd,int c,int d) const;

      static int get_inco(int S,int S_ab,int a,int b);

      double trace() const;

      double ddot(const dDPM &) const;

      void proj_W();

      void test_proj_1() const;

      void test_proj_2() const;

      void hubbard(double);

      void unit();

      void constr_grad(double,const dDPM &,const SUP &);

      void proj_Tr();

      void proj();

      int solve(double,const SUP &,dDPM &);

      void H(double,const dDPM &,const SUP &);

      double line_search(double,SUP &,const dDPM &);

      double line_search(double,const dDPM &,const dDPM &);

      void Q(char option,const dDPM &);

      double dotunit() const;

      void I(const dPPHM &);

      void Q(const dPPHM &);

   private:

};

#endif
