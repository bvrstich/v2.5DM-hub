#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;
using std::ofstream;
using std::ifstream;
using std::cout;
using std::endl;
using std::ios;

#include "include.h"

/**
 * standard constructor
 */
xSPM::xSPM() : BlockMatrix(2) {

   //set the dimension and degeneracy of the two blocks:
   this->setMatrixDim(0,Tools::gL(),1);
   this->setMatrixDim(1,Tools::gL(),3);

}

/**
 * copy constructor: constructs Matrix object and fills it with the content of matrix xspm_c
 * @param xspm_c object that will be copied into this.
 */
xSPM::xSPM(const xSPM &xspm_c) : BlockMatrix(xspm_c){ }

/**
 * destructor
 */
xSPM::~xSPM(){ }
