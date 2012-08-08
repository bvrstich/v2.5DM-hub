#ifndef rxTPM_H
#define rxTPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

/**
 * @author Brecht Verstichel
 * @date 08-08-2012\n\n
 * This class rxTPM is a class written for the blocks of the dDPM matrices, it is called rxTPM because it is a reduced expansion of a TPM matrix.
 * expanded for the same reason as xTPM and reduced because there is an extra parameter l which the sp indices cannot be equal to in the S = 3/2 block.
 * This is written for translationally invariant dDPM matrices, so the only the l = 0 part is stored.
 */
class rxTPM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param rxtpm_p the rxTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const rxTPM &rxtpm_p);

   public:
      
      //constructor
      rxTPM();

      //copy constructor
      rxTPM(const rxTPM &);

      //destructor
      virtual ~rxTPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      void pseudo_invert(int,int);

      void pseudo_sqrt(int,int,int);

      static void init();

      static void clear();

      static int gt2s(int,int,int);

      static int gs2t(int,int,int,int);

   private:

      //!static list that takes a dp-spinindex S and a tp index i and returns two sp indices a and b and intermediate spin S_ab
      static vector< vector<int> > *t2s;

      //!static list that takes in a dp-spinindex S, an intermediate spin S_ab and two sp indices a,b and returns a tp index i
      static int ****s2t;

};

#endif
