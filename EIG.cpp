#include <iostream>
#include <fstream>
#include <cmath>

using std::ostream;
using std::endl;

#include "include.h"

int EIG::dim;

/**
 * initialize the statics
 */
void EIG::init(){

   dim = Tools::gL()*(2*Tools::gL() - 2)*(2*Tools::gL() - 1);

#ifdef __Q2_CON
   dim += Tools::gL()*(2*Tools::gL() - 2)*(2*Tools::gL() - 1);
#endif

#ifdef __I2_CON
   dim += 2*Tools::gL()*Tools::gL()*(2*Tools::gL() - 1);
#endif

#ifdef __Q1_CON
   dim += 2*Tools::gL()*Tools::gL()*(2*Tools::gL() - 1);
#endif

#ifdef __G1_CON
   dim += 4*Tools::gL()*Tools::gL()*(2*Tools::gL() - 1);
#endif

#ifdef __G2_CON
   dim += 4*Tools::gL()*Tools::gL()*(2*Tools::gL() - 1);
#endif

}

/**
 * standard constructor with initialization on the eigenvalues of a SUP object.
 * @param X input SUP object that will be destroyed after this function is called. The eigenvectors
 * of the matrix will be stored in the columns of the original SUP matrix.
 */
EIG::EIG(SUP &X){

   v_I1 = new BlockVector<dDPM>(X.gI1());

#ifdef __Q2_CON
   v_Q2 = new BlockVector<dDPM>(X.gQ2());
#endif 

#ifdef __I2_CON
   v_I2 = new BlockVector<dPPHM>(X.gI2());
#endif 

#ifdef __Q1_CON
   v_Q1 = new BlockVector<dPPHM>(X.gQ1());
#endif 

#ifdef __G1_CON
   v_G1 = new BlockVector<dPHHM>(X.gG1());
#endif 

#ifdef __G2_CON
   v_G2 = new BlockVector<dPHHM>(X.gG2());
#endif 

}

/**
 * Copy constructor\n
 * allocates the memory for the eigenvalues of a SUP object and copies the content of eig_c into it.
 * @param eig_c The input EIG that will be copied into this.
 */
EIG::EIG(const EIG &eig_c){

   v_I1 = new BlockVector<dDPM>(eig_c.gv_I1());

#ifdef __Q2_CON
   v_Q2 = new BlockVector<dDPM>(eig_c.gv_Q2());
#endif 

#ifdef __I2_CON
   v_I2 = new BlockVector<dPPHM>(eig_c.gv_I2());
#endif 

#ifdef __Q1_CON
   v_Q1 = new BlockVector<dPPHM>(eig_c.gv_Q1());
#endif 

#ifdef __G1_CON
   v_G1 = new BlockVector<dPHHM>(eig_c.gv_G1());
#endif 

#ifdef __G2_CON
   v_G2 = new BlockVector<dPHHM>(eig_c.gv_G2());
#endif 

}

/**
 * overload equality operator
 * @param eig_c object that will be copied into this.
 */
EIG &EIG::operator=(const EIG &eig_c){

   *v_I1 = eig_c.gv_I1();

#ifdef __Q2_CON
   *v_Q2 = eig_c.gv_Q2();
#endif 

#ifdef __I2_CON
   *v_I2 = eig_c.gv_I2();
#endif 

#ifdef __Q1_CON
   *v_Q1 = eig_c.gv_Q1();
#endif 

#ifdef __G1_CON
   *v_G1 = eig_c.gv_G1();
#endif 

#ifdef __G2_CON
   *v_G2 = eig_c.gv_G2();
#endif 

   return *this;

}

/**
 * Destructor, deallocation of the memory
 */
EIG::~EIG(){

   delete v_I1;

#ifdef __Q2_CON
   delete v_Q2;
#endif 

#ifdef __I2_CON
   delete v_I2;
#endif 

#ifdef __Q1_CON
   delete v_Q1;
#endif 

#ifdef __G1_CON
   delete v_G1;
#endif 

#ifdef __G2_CON
   delete v_G2;
#endif 

}

ostream &operator<<(ostream &output,const EIG &eig_p){

   std::cout << eig_p.gv_I1() << std::endl;

#ifdef __Q2_CON
   std::cout << eig_p.gv_Q2() << std::endl;
#endif 

#ifdef __I2_CON
   std::cout << eig_p.gv_I2() << std::endl;
#endif 

#ifdef __Q1_CON
   std::cout << eig_p.gv_Q1() << std::endl;
#endif 

#ifdef __G1_CON
   std::cout << eig_p.gv_G1() << std::endl;
#endif 

#ifdef __G2_CON
   std::cout << eig_p.gv_G2() << std::endl;
#endif 

   return output;

}

/** 
 * get the BlockVector<dDPM> object containing the eigenvalues of the dDPM block I1
 * @return a BlockVector<dDPM> object containing the desired eigenvalues
 */
BlockVector<dDPM> &EIG::gv_I1(){

   return *v_I1;

}

/** 
 * the const version
 * get the BlockVector<dDPM> object containing the eigenvalues of the dDPM block I1
 * @return a BlockVector<dDPM> object containing the desired eigenvalues
 */
const BlockVector<dDPM> &EIG::gv_I1() const{

   return *v_I1;

}

#ifdef __Q2_CON

/** 
 * get the BlockVector<dDPM> object containing the eigenvalues of the dDPM block Q2
 * @return a BlockVector<dDPM> object containing the desired eigenvalues
 */
BlockVector<dDPM> &EIG::gv_Q2(){

   return *v_Q2;

}

/** 
 * the const version
 * get the BlockVector<dDPM> object containing the eigenvalues of the dDPM block Q2
 * @return a BlockVector<dDPM> object containing the desired eigenvalues
 */
const BlockVector<dDPM> &EIG::gv_Q2() const{

   return *v_Q2;

}

#endif

#ifdef __I2_CON

/** 
 * get the BlockVector<dPPHM> object containing the eigenvalues of the dDPM block I2
 * @return a BlockVector<dPPHM> object containing the desired eigenvalues
 */
BlockVector<dPPHM> &EIG::gv_I2(){

   return *v_I2;

}

/** 
 * the const version
 * get the BlockVector<dPPHM> object containing the eigenvalues of the dDPM block I2
 * @return a BlockVector<dPPHM> object containing the desired eigenvalues
 */
const BlockVector<dPPHM> &EIG::gv_I2() const{

   return *v_I2;

}

#endif

#ifdef __Q1_CON

/** 
 * get the BlockVector<dPPHM> object containing the eigenvalues of the dDPM block Q1
 * @return a BlockVector<dPPHM> object containing the desired eigenvalues
 */
BlockVector<dPPHM> &EIG::gv_Q1(){

   return *v_Q1;

}

/** 
 * the const version
 * get the BlockVector<dPPHM> object containing the eigenvalues of the dDPM block Q1
 * @return a BlockVector<dPPHM> object containing the desired eigenvalues
 */
const BlockVector<dPPHM> &EIG::gv_Q1() const{

   return *v_Q1;

}

#endif

#ifdef __G1_CON

/** 
 * get the BlockVector<dPHHM> object containing the eigenvalues of the dDPM block G1
 * @return a BlockVector<dPHHM> object containing the desired eigenvalues
 */
BlockVector<dPHHM> &EIG::gv_G1(){

   return *v_G1;

}

/** 
 * the const version
 * get the BlockVector<dPHHM> object containing the eigenvalues of the dDPM block G1
 * @return a BlockVector<dPHHM> object containing the desired eigenvalues
 */
const BlockVector<dPHHM> &EIG::gv_G1() const{

   return *v_G1;

}

#endif

#ifdef __G2_CON

/** 
 * get the BlockVector<dPHHM> object containing the eigenvalues of the dDPM block G2
 * @return a BlockVector<dPHHM> object containing the desired eigenvalues
 */
BlockVector<dPHHM> &EIG::gv_G2(){

   return *v_G2;

}

/** 
 * the const version
 * get the BlockVector<dPHHM> object containing the eigenvalues of the dDPM block G2
 * @return a BlockVector<dPHHM> object containing the desired eigenvalues
 */
const BlockVector<dPHHM> &EIG::gv_G2() const{

   return *v_G2;

}

#endif

/**
 * @return total dimension of the EIG object
 */
int EIG::gdim() const{

   return dim;

}

/**
 * @return the minimal element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::min() const{

   //lowest eigenvalue of P block
   double ward = v_I1->min();

#ifdef __Q2_CON
   if(v_Q2->min() < ward)
      ward = v_Q2->min();
#endif

#ifdef __I2_CON
   if(v_I2->min() < ward)
      ward = v_I2->min();
#endif

#ifdef __Q1_CON
   if(v_Q1->min() < ward)
      ward = v_Q1->min();
#endif

#ifdef __G1_CON
   if(v_G1->min() < ward)
      ward = v_G1->min();
#endif

#ifdef __G2_CON
   if(v_G2->min() < ward)
      ward = v_G2->min();
#endif

   return ward;

}

/**
 * @return the maximum element present in this EIG object.
 * watch out, only works when EIG is filled with the eigenvalues of a diagonalized SUP matrix
 */
double EIG::max() const{

   //highest eigenvalue of P block
   double ward = v_I1->max();

#ifdef __Q2_CON
   if(v_Q2->max() > ward)
      ward = v_Q2->max();
#endif

#ifdef __I2_CON
   if(v_I2->max() > ward)
      ward = v_I2->max();
#endif

#ifdef __Q1_CON
   if(v_Q1->max() > ward)
      ward = v_Q1->max();
#endif

#ifdef __G1_CON
   if(v_G1->max() > ward)
      ward = v_G1->max();
#endif

#ifdef __G2_CON
   if(v_G2->max() > ward)
      ward = v_G2->max();
#endif

   return ward;

}

/**
 * @param alpha step length along the Newton direction
 * @return The line search function, gradient of the potential in the Newton direction as a function of the step length alpha
 */
double EIG::lsfunc(double alpha) const{

   double ward = v_I1->lsfunc(alpha);

#ifdef __Q2_CON
   ward += v_Q2->lsfunc(alpha);
#endif

#ifdef __I2_CON
   ward += v_I2->lsfunc(alpha);
#endif

#ifdef __Q1_CON
   ward += v_Q1->lsfunc(alpha);
#endif

#ifdef __G1_CON
   ward += v_G1->lsfunc(alpha);
#endif

#ifdef __G2_CON
   ward += v_G2->lsfunc(alpha);
#endif

   return Tools::gL() * ward;

}
