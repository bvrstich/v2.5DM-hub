#ifndef xTPM_H
#define xTPM_H

#include <iostream>
#include <fstream>
#include <vector>

using std::ostream;
using std::vector;

#include "Matrix.h"

class dDPM;

/**
 * @author Brecht Verstichel
 * @date 04-08-2011\n\n
 * This class xTPM is a class written for spinsymmetrical the blocks of the dPPHM matrices. It is called xTPM for expanded TPM.
 * because it is a two particle matrix, but it is expanded because it is assumed that it can couple with a third index, so it has
 * more spinfreedom then a regular TPM.
 */
class xTPM : public BlockMatrix {

   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param xtpm_p the xTPM you want to print
    */
   friend ostream &operator<<(ostream &output,const xTPM &xtpm_p);

   public:
      
      //constructor
      xTPM();

      //copy constructor
      xTPM(const xTPM &);

      //destructor
      virtual ~xTPM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      static void init();

      static void clear();

      static int gt2s(int,int,int);

      static int gs2t(int,int,int,int);

      static void print_basis();

   private:

      //!static list that takes in a dp-spinindex S and a tp index i and returns two sp indices a and b and intermediate spin S_ab
      static vector< vector<int> > *t2s;

      //!static list that takes in a dp-spinindex S, intermediate spin S_ab and two sp indices a,b and returns a tp index i
      static int ****s2t;

};

#endif
