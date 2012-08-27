#ifndef SPM_H
#define SPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "Matrix.h"

class TPM;
class dDPM;
class dPPHM;
class dPHHM;

/**
 * @author Brecht Verstichel
 * @date 20-04-2010\n\n
 * This class SPM was written for single particle matrices in a spinsymmetrical system. It inherits from the class Matrix and expands it with
 * specific memberfunction and a knowledge of the nr of sp orbitals and particles.
 */

class SPM : public Matrix {

   /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << spm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << spm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param spm_p de SPM you want to print
    */
   friend ostream &operator<<(ostream &output,const SPM &spm_p);

   public:
      
      //constructor
      SPM();

      //copy constructor
      SPM(const SPM &);

      //destructor
      virtual ~SPM();

      using Matrix::operator=;

      void bar(double,const TPM &);

      void breve(double,const dDPM &);

      void breve(double,const dPPHM &);

      void breve(double,const dPHHM &);

      void barbar(double,const dPHHM &);

      void breve_si(double,const dPHHM &);

   private:

};

#endif
