#ifndef SUP_H
#define SUP_H

#include <iostream> 
#include <fstream> 

using std::ostream;

class dDPM;
class dPPHM;
class dPHHM;
class EIG;

/**
 * @author Brecht Verstichel
 * @date 22-08-2012\n\n
 * This class, SUP is a blockmatrix over the carrierspace's of active N-representability conditions. 
 * This class contains two dDPM objects containing the I1 and Q2 conditions, two dPPHM objects containing I2 and Q1
 * and two dPHHM objects containing G1 and G2.
 */
class SUP{
  
   /**
    * Output stream operator overloaded
    * @param output The stream to which you are writing (e.g. cout)
    * @param X_p the SUP you want to print
    */
   friend ostream &operator<<(ostream &output,const SUP &X_p);

   public:

      //constructor
      SUP();

      //copy constructor
      SUP(const SUP &);

      //destructor
      ~SUP();

      //overload += operator
      SUP &operator+=(const SUP &);

      //overload -= operator
      SUP &operator-=(const SUP &);

      //overload equality operator
      SUP &operator=(const SUP &);

      //overload equality operator
      SUP &operator=(double);

      dDPM &gI1();

      const dDPM &gI1() const;

#ifdef __Q2_CON
      dDPM &gQ2();

      const dDPM &gQ2() const;
#endif

#ifdef __I2_CON
      dPPHM &gI2();

      const dPPHM &gI2() const;
#endif

#ifdef __Q1_CON
      dPPHM &gQ1();

      const dPPHM &gQ1() const;
#endif

#ifdef __G1_CON
      dPHHM &gG1();

      const dPHHM &gG1() const;
#endif

#ifdef __G2_CON
      dPHHM &gG2();

      const dPHHM &gG2() const;
#endif

      int gdim() const;

      double ddot(const SUP &) const;

      void invert();

      void dscal(double alpha);

      //positieve of negatieve vierkantswortel uit een supermatrix
      void sqrt(int option);

      void L_map(const SUP &,const SUP &);

      void daxpy(double alpha,const SUP &);

      double trace() const;

      SUP &mprod(const SUP &,const SUP &);

      void fill(const dDPM &);

      void fill();

      void fill_Random();

      void sep_pm(SUP &,SUP &);

      static void init();

   private:

      //!pointer to a dDPM object, containts the W matrix that has to be positive semidefinite: I1 condition
      dDPM *I1;

#ifdef __Q2_CON
      //!pointer to a dDPM object, containts the Q2(W) matrix that has to be positive semidefinite: Q2 condition
      dDPM *Q2;
#endif

#ifdef __I2_CON
      //!pointer to a dPPHM object, containts the I2(W) matrix that has to be positive semidefinite: I2 condition
      dPPHM *I2;
#endif

#ifdef __Q1_CON
      //!pointer to a dPPHM object, containts the Q1(W) matrix that has to be positive semidefinite: Q1 condition
      dPPHM *Q1;
#endif

#ifdef __G1_CON
      //!pointer to a dPHHM object, containts the G1(W) matrix that has to be positive semidefinite: G1 condition
      dPHHM *G1;
#endif

#ifdef __G2_CON
      //!pointer to a dPHHM object, containts the G2(W) matrix that has to be positive semidefinite: G2 condition
      dPHHM *G2;
#endif

      //!total dimension of the SUP matrix
      static int dim;

};

#endif
