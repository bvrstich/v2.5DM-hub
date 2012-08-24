#ifndef ssdTPM_H
#define ssdTPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "SPM.h"

/**
 * @author Brecht Verstichel
 * @date 23-08-2012\n\n
 * This class ssdTPM is an array of SPM matrices, it will contain a special "bar of the different blocks of the dDPM and consorts.
 * It is a called ssdTPM because it has no spin dependence anymore, unlike the true dTPM, so spinsummed (ss) dTPM.
 * Because translational invariance is exploited, only one SPM block is actually stored.
 */
class ssdTPM : public SPM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ssdtpm_p the ssdTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const ssdTPM &ssdtpm_p);

   public:

      //constructor
      ssdTPM();

      //copy constructor
      ssdTPM(const ssdTPM &);

      //destructor
      virtual ~ssdTPM();

      using Matrix::operator=;

      using Matrix::operator();

      double operator()(int,int,int) const;

      double trace() const;

      double ddot(const ssdTPM &) const;

      void bar(double,const dDPM &);

      void bar(double,const dPPHM &);

   private:

};

#endif
