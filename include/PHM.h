#ifndef PHM_H
#define PHM_H

#include <iostream>

using std::ostream;

#include "BlockMatrix.h"
#include "TPM.h"

class dPHHM;

/**
 * @author Brecht Verstichel
 * @date 23-04-2010\n\n
 * This class, PHM, is a class written for spinsymmetrical particle-hole matrices, it inherits all the functions from its mother class
 * BlockMatrix, some special member functions and two lists that give the relationship between the sp and the ph basis.
 */
class PHM : public BlockMatrix {

    /**
    * Output stream operator overloaded, the usage is simple, if you want to print to a file, make an
    * ifstream object and type:\n\n
    * object << phm_p << endl;\n\n
    * For output onto the screen type: \n\n
    * cout << phm_p << endl;\n\n
    * @param output The stream to which you are writing (e.g. cout)
    * @param phm_p the PHM you want to print
    */
   friend ostream &operator<<(ostream &output,const PHM &phm_p);

   public:
      
      //constructor
      PHM();

      //copy constructor
      PHM(const PHM &);

      //destructor
      virtual ~PHM();

      using BlockMatrix::operator=;

      using BlockMatrix::operator();

      //change the numbers in sp mode
      double &operator()(int S,int a,int b,int c,int d);

      //change the numbers in sp mode: read mode
      double operator()(int S,int a,int b,int c,int d) const;

      void G(const TPM &);

      void spinsum(double,const dDPM &);

      void spinsum(double,const dPPHM &);

      void spinsum(double,const dPHHM &);

      void bar(double,const dPHHM &);

      void spinsum_si(double,const dPHHM &);

      static void init();

      static void clear();

   private:

      //!static list of dimension [n_ph][2] that takes in a ph index i and returns two sp indices: a = ph2s[i][0] and b = ph2s[i][1]
      static vector< vector<int> > ph2s;

      //!static list of dimension [M/2][M/2] that takes two sp indices a,b and returns a ph index i: i = s2ph[a][b]
      static int **s2ph;

};

#endif
