#ifndef rxTPM_H
#define rxTPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class dDPM;

/**
 * @author Brecht Verstichel
 * @date 29-09-2011\n\n
 * This class rxTPM is a class written for the blocks of the dDPM matrices, it is called rxTPM because it is a reduced expansion of a TPM matrix.
 * expanded for the same reason as xTPM and reduced because there is an extra parameter l which the sp indices cannot be equal to in the S = 3/2 block.
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
      rxTPM(int);

      //copy constructor
      rxTPM(const rxTPM &);

      //destructor
      virtual ~rxTPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //geef N terug
      int gN() const;

      //geef M terug
      int gM() const;

      //geef l terug
      int gl() const;

      void pseudo_invert();

      void pseudo_sqrt(int);

      static void init(int,int);

      static void clear();

      static int gt2s(int,int,int,int);

      static int gs2t(int,int,int,int,int);

      static void print_basis();

   private:

      //!static list that takes in a paramater l, a dp-spinindex S and a tp index i and returns two sp indices a and b and intermediate spin S_ab
      static vector< vector<int> > **t2s;

      //!static list that takes in a parameter l, a dp-spinindex S, intermediate spin S_ab and two sp indices a,b and returns a tp index i
      static int *****s2t;

      //!parameter l, sp index that is blocked out from the a = b part.
      int l;

      //!nr of particles
      static int N;

      //!dimension of sp hilbert space
      static int M;

};

#endif
