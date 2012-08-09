#ifndef dDPM_H
#define dDPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "rxTPM.h"

/**
 * @author Brecht Verstichel
 * @date 08-08-2012\n\n
 * This class dDPM is a class written for the I1 and Q2 conditions, which is a DPM matrix with the spacial third index diagonal.
 * The v2.5DM program will use this as the central variable. Translational invariance is exploited, which means only one block is
 * to be stored. It enherits from the rxTPM class and overloads the functions ddot and trace.
 */
class dDPM : public rxTPM {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param ddpm_p the dDPM you want to print
    */
   friend ostream &operator<<(ostream &output,const dDPM &ddpm_p);

   public:

      //constructor
      dDPM();

      //copy constructor
      dDPM(const dDPM &);

      //destructor
      virtual ~dDPM();

      double trace() const;

      double ddot(const dDPM &) const;

   private:

};

#endif
