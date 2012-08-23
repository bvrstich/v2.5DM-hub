#ifndef xSPM_H
#define xSPM_H

#include <iostream>
#include <fstream>

using std::ostream;

#include "BlockMatrix.h"

/**
 * @author Brecht Verstichel
 * @date 19-04-2010\n\n
 * This class xSPM is a class written for the blocks of the dTPM object, which is a TPM object with the second index diagonal
 * It consists of two blocks with equal dimensions and inherits from BlockMatrix.
 */
class xSPM : public BlockMatrix {

   public:
      
      xSPM();

      //copy constructor
      xSPM(const xSPM &);

      //destructor
      virtual ~xSPM();

      void constr_lists();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

   private:

};

#endif
