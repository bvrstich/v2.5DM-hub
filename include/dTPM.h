#ifndef dTPM_H
#define dTPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "xSPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-08-2012\n\n
 * This class dTPM is an array of xSPM matrices, it will contain a special "bar of the different blocks of the dDPM and consorts.
 * Because translational invariance is exploited, only one xSPM block is actually stored.
 */
class dTPM : public xSPM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param dtpm_p the dTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const dTPM &dtpm_p);

   public:

      //constructor
      dTPM();

      //copy constructor
      dTPM(const dTPM &);

      //destructor
      virtual ~dTPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      double operator()(int,int,int,int) const;

      double trace() const;

      double ddot(const dTPM &) const;

      void bar(double,const dDPM &);

      void bar(double,const dPPHM &);

   private:

};

#endif
