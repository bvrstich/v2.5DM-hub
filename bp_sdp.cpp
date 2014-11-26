/**
 * @mainpage 
 * This is an implementation of the dual only, potential reduction interior point method
 * for optimizing the 2.5DM using the I1, I2, Q1, Q2, G1 and G2 N-representability conditions.
 * Compiling can be done with the options I1 , I, I1Q2, IQ, IQG1, IQG2 and IQG with logical consequences for the program.
 * @author Brecht Verstichel, Ward Poelmans
 * @date 10-06-2011
 */

#include <iostream>
#include <fstream>
#include <cmath>
#include <gsl/gsl_sf_coupling.h>

using std::cout;
using std::endl;
using std::ifstream;
using std::ofstream;

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

/**
 * In the main the actual program is run.\n 
 * We start from the unity 2.5DM normed on the particle number and minimize the ojective function:\n\n
 * Tr (Gamma H) - t * ln(det P(Gamma)) \n\n
 * Once the minimum is found the parameter t is reduced and a new search is initiated,
 * this goes on until convergence is reached.\n
 * The potential is minimized using the Newton-Raphson method and the resulting linear system
 * is solved via the linear conjugate gradient method.
 */
int main(void){

   srand(time(NULL));

   cout.precision(10);

   const int L = 4;//nr of sites;
   const int N = 4;//nr of particles

   Tools::init(L,N);

   rxTPM::init();
   TPM::init();
   PHM::init();
   xTPM::init();
   rxPHM::init();

   SUP::init();
   EIG::init();

   //make the unit matrix for the proj_Tr
   Tools::sunit();

   //hamiltoniaan
   dDPM ham;
   ham.hubbard(1.0);

   dDPM ham_copy(ham);

   //only symmetrical, traceless hamiltonian needed in program.
   ham.proj();

   //primal
   SUP X;

   //dual
   SUP Z;

   //Lagrange multiplier
   SUP V;

   //just dubya
   SUP W;

   SUP u_0;

   //little help
   dDPM hulp;

   u_0.gI1().set_u_0();

   u_0.fill();

   X = 0.0;
   Z = 0.0;

   //what does this do?
   double sigma = 1.0;

   double tolerance = 1.0e-5;

   double D_conv(1.0),P_conv(1.0),convergence(1.0);

   // mazziotti uses 1.6 for this
   double mazzy = 2.0;

   int iter;
   int max_iter = 1;

   while(D_conv > tolerance || P_conv > tolerance || fabs(convergence) > tolerance){

      D_conv = 1.0;

      iter = 0;

      while(D_conv > tolerance && iter <= max_iter){

         ++iter;

         //solve system
         SUP B(Z);

         B -= u_0;

         B.daxpy(mazzy/sigma,X);

         dDPM b;

         b.collaps(B);

         b.daxpy(-mazzy/sigma,ham);

         hulp.solve(b);

         //construct W
         W.fill(hulp);

         W += u_0;

         W.daxpy(-1.0/sigma,X);

         //update Z and V with eigenvalue decomposition:
         W.sep_pm(Z,V);

         V.dscal(-sigma);

         //check infeasibility of the dual problem:
         dDPM v;

         v.collaps(V);

         v -= ham;

         D_conv = sqrt(v.ddot(v));

     }

      //update primal:
      X = V;

      //check primal feasibility (W is a helping variable now)
      W.fill(hulp);

      W += u_0;

      W -= Z;

      P_conv = sqrt(W.ddot(W));

      convergence = Z.gI1().ddot(ham) + X.ddot(u_0);

      cout << Z.gI1().trace() << "\t" << P_conv << "\t" << D_conv << "\t" << sigma << "\t" << convergence << "\t" << Z.gI1().ddot(ham_copy) << endl;

      if(D_conv < P_conv)
         sigma *= 1.01;
      else
         sigma /= 1.01;

   }

   cout << endl;
   cout << "Trace" << "\t" << Z.gI1().trace() << endl;
   cout << "Energy: " << ham_copy.ddot(Z.gI1()) << endl;
   cout << "pd gap: " << Z.ddot(X) << endl;
   cout << "dual conv: " << D_conv << endl;
   cout << "primal conv: " << P_conv << endl;

   cout << Z.gI1() << endl;

   Tools::clear();

   rxPHM::clear();
   rxTPM::clear();
   TPM::clear();
   PHM::clear();
   xTPM::clear();

   return 0;

}
