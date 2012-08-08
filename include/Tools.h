#ifndef TOOLS_H
#define TOOLS_H

#include <iostream>
#include <cstdlib>

/**
 * @author Brecht Verstichel
 * @date 26-09-2011\n\n
 * This is a class that contains some static objects (lists, matrices) needed by all the classes,
 * and it would be stupid to store it in every class.
 */

class Tools{

   public:

      static void init(int,int);

      static void clear();

      static double g6j(int,int,int,int);

      static double g9j(int,int,int,int);

      static int gL();

      static int gN();

   private:

      //!array in which I store 6j symbols needed
      static double *x6j;

      //!array in which I store 9j symbols needed
      static double *x9j;

      //!nr of sites
      static int L;

      //!nr of particles
      static int N;

};

#endif
