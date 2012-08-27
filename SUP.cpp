#include <iostream>
#include <cmath>
#include <fstream>

using std::ostream;

#include "include.h"

int SUP::dim;

/**
 * initialize the statics
 */
void SUP::init(){

   dim = Tools::gL()*(2*Tools::gL() - 2)*(2*Tools::gL() - 1);

#ifdef __Q2_CON
   dim += Tools::gL()*(2*Tools::gL() - 2)*(2*Tools::gL() - 1);
#endif

#ifdef __I2_CON
   dim += 2 * Tools::gL() * Tools::gL() * (2*Tools::gL() - 1);
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
 * standard constructor\n
 * Allocates two dDPM matrices, I1 and Q2 and one dPPHM matrix I2
 */
SUP::SUP(){

   I1 = new dDPM();

#ifdef __Q2_CON
   Q2 = new dDPM();
#endif

#ifdef __I2_CON
   I2 = new dPPHM();
#endif

#ifdef __Q1_CON
   Q1 = new dPPHM();
#endif

#ifdef __G1_CON
   G1 = new dPHHM();
#endif

#ifdef __G2_CON
   G2 = new dPHHM();
#endif

}

/**
 * standard constructor\n
 * Allocates two dDPM's and a dPPHM, then copies the content of input SUP X into it.
 * @param X_c input SUP
 */
SUP::SUP(const SUP &X_c){

   I1 = new dDPM(X_c.gI1());

#ifdef __Q2_CON
   Q2 = new dDPM(X_c.gQ2());
#endif

#ifdef __I2_CON
   I2 = new dPPHM(X_c.gI2());
#endif

#ifdef __Q1_CON
   Q1 = new dPPHM(X_c.gQ1());
#endif

#ifdef __G1_CON
   G1 = new dPHHM(X_c.gG1());
#endif

#ifdef __G2_CON
   G2 = new dPHHM(X_c.gG2());
#endif

}

/**
 * Destructor
 */
SUP::~SUP(){

   delete I1;

#ifdef __Q2_CON
   delete Q2;
#endif

#ifdef __I2_CON
   delete I2;
#endif

#ifdef __Q1_CON
   delete Q1;
#endif

#ifdef __G1_CON
   delete G1;
#endif

#ifdef __G2_CON
   delete G2;
#endif

}

/**
 * Overload += operator
 * @param X_pl The SUP matrix that has to be added to this
 */
SUP &SUP::operator+=(const SUP &X_pl){

   *I1 += X_pl.gI1();

#ifdef __Q2_CON
   *Q2 += X_pl.gQ2();
#endif

#ifdef __I2_CON
   *I2 += X_pl.gI2();
#endif

#ifdef __Q1_CON
   *Q1 += X_pl.gQ1();
#endif

#ifdef __G1_CON
   *G1 += X_pl.gG1();
#endif

#ifdef __G2_CON
   *G2 += X_pl.gG2();
#endif

   return *this;

}

/**
 * Overload -= operator
 * @param X_pl The SUP that will be deducted from this
 */
SUP &SUP::operator-=(const SUP &X_pl){

   *I1 -= X_pl.gI1();

#ifdef __Q2_CON
   *Q2 -= X_pl.gQ2();
#endif

#ifdef __I2_CON
   *I2 -= X_pl.gI2();
#endif

#ifdef __Q1_CON
   *Q1 -= X_pl.gQ1();
#endif

#ifdef __G1_CON
   *G1 -= X_pl.gG1();
#endif

#ifdef __G2_CON
   *G2 -= X_pl.gG2();
#endif

   return *this;

}

/**
 * Overload equality operator, copy X_c into this
 * @param X_c SUP to be copied into this
 */
SUP &SUP::operator=(const SUP &X_c){

   *I1 = X_c.gI1();

#ifdef __Q2_CON
   *Q2 = X_c.gQ2();
#endif

#ifdef __I2_CON
   *I2 = X_c.gI2();
#endif

#ifdef __Q1_CON
   *Q1 = X_c.gQ1();
#endif

#ifdef __G1_CON
   *G1 = X_c.gG1();
#endif

#ifdef __G2_CON
   *G2 = X_c.gG2();
#endif

   return *this;

}

/**
 * overload operator = number, all the blockmatrices in SUP are put equal to the number a.
 * e.g. X = 0 makes all the Matrix elements zero.
 * @param a the number
 */
SUP &SUP::operator=(double a){

   *I1 = a;

#ifdef __Q2_CON
   *Q2 = a;
#endif

#ifdef __I2_CON
   *I2 = a;
#endif

#ifdef __Q1_CON
   *Q1 = a;
#endif

#ifdef __G1_CON
   *G1 = a;
#endif

#ifdef __G2_CON
   *G2 = a;
#endif

   return *this;

}

/**
 * @return pointer to the dDPM object containing the I1 block
 */
dDPM &SUP::gI1(){

   return *I1;

}

/**
 * the const version
 * @return pointer to the dDPM object containing the I1 block
 */
const dDPM &SUP::gI1() const{

   return *I1;

}

#ifdef __Q2_CON

/**
 * @return pointer to the dDPM object containing the Q2 block
 */
dDPM &SUP::gQ2(){

   return *Q2;

}

/**
 * the const version
 * @return pointer to the dDPM object containing the Q2 block
 */
const dDPM &SUP::gQ2() const{

   return *Q2;

}

#endif

#ifdef __I2_CON

/**
 * @return pointer to the dPPHM object containing the I2 block
 */
dPPHM &SUP::gI2(){

   return *I2;

}

/**
 * the const version
 * @return pointer to the dPPHM object containing the I2 block
 */
const dPPHM &SUP::gI2() const{

   return *I2;

}

#endif

#ifdef __Q1_CON

/**
 * @return pointer to the dPPHM object containing the Q1 block
 */
dPPHM &SUP::gQ1(){

   return *Q1;

}

/**
 * the const version
 * @return pointer to the dPPHM object containing the Q1 block
 */
const dPPHM &SUP::gQ1() const{

   return *Q1;

}

#endif

#ifdef __G1_CON

/**
 * @return pointer to the dPHHM object containing the G1 block
 */
dPHHM &SUP::gG1(){

   return *G1;

}

/**
 * the const version
 * @return pointer to the dPHHM object containing the G1 block
 */
const dPHHM &SUP::gG1() const{

   return *G1;

}

#endif

#ifdef __G2_CON

/**
 * @return pointer to the dPHHM object containing the G2 block
 */
dPHHM &SUP::gG2(){

   return *G2;

}

/**
 * the const version
 * @return pointer to the dPHHM object containing the G2 block
 */
const dPHHM &SUP::gG2() const{

   return *G2;

}

#endif

ostream &operator<<(ostream &output,const SUP &X_p){

   output << (rxTPM)(X_p.gI1()) << std::endl;

#ifdef __Q2_CON
   output << (rxTPM)(X_p.gQ2()) << std::endl;
#endif

#ifdef __I2_CON
   output << (xTPM)(X_p.gI2()) << std::endl;
#endif

#ifdef __Q1_CON
   output << (xTPM)(X_p.gQ1()) << std::endl;
#endif

#ifdef __G1_CON
   output << (X_p.gG1()) << std::endl;
#endif

#ifdef __G2_CON
   output << (X_p.gG2()) << std::endl;
#endif

   return output;

}

/**
 * Fill the SUP matrix with random elements, watch out, not necessarily positive definite
 */
void SUP::fill_Random(){

   I1->fill_Random();

#ifdef __Q2_CON
   Q2->fill_Random();
#endif

#ifdef __I2_CON
   I2->fill_Random();
#endif

#ifdef __Q1_CON
   Q1->fill_Random();
#endif

#ifdef __G1_CON
   G1->fill_Random();
#endif

#ifdef __G2_CON
   G2->fill_Random();
#endif

}

/**
 * @return total dimension of SUP (carrier) space
 */
int SUP::gdim() const{

   return dim;

}

/**
 * @param X_i input SUP_PQ X_i
 * @return inproduct between this and input matrix X_i, defined as Tr(this X_i)
 */
double SUP::ddot(const SUP &X_i) const{

   double ward = 0.0;

   ward += I1->ddot(X_i.gI1());

#ifdef __Q2_CON
   ward += Q2->ddot(X_i.gQ2());
#endif

#ifdef __I2_CON
   ward += I2->ddot(X_i.gI2());
#endif

#ifdef __Q1_CON
   ward += Q1->ddot(X_i.gQ1());
#endif

#ifdef __G1_CON
   ward += G1->ddot(X_i.gG1());
#endif

#ifdef __G2_CON
   ward += G2->ddot(X_i.gG2());
#endif

   return ward;

}

/**
 * Invert all the Matrices in SUP and put it in this, watch out, destroys original matrices.
 * Makes use of cholesky decomposition, so only positive definite matrices can be used as input!
 */
void SUP::invert(){

   I1->pseudo_invert(Tools::gL() - 1,0);

#ifdef __Q2_CON
   Q2->pseudo_invert(Tools::gL() - 1,0);
#endif

#ifdef __I2_CON
   I2->invert();
#endif

#ifdef __Q1_CON
   Q1->invert();
#endif

#ifdef __G1_CON
   G1->invert();
#endif

#ifdef __G2_CON
   G2->invert();
#endif

}

/**
 * Scale all the elements in the blockmatrix with parameter alpha
 * @param alpha the scalefactor
 */
void SUP::dscal(double alpha){

   I1->dscal(alpha);

#ifdef __Q2_CON
   Q2->dscal(alpha);
#endif

#ifdef __I2_CON
   I2->dscal(alpha);
#endif

#ifdef __Q1_CON
   Q1->dscal(alpha);
#endif

#ifdef __G1_CON
   G1->dscal(alpha);
#endif

#ifdef __G2_CON
   G2->dscal(alpha);
#endif

}

/**
 * Take the square root out of the positive semidefinite SUP matrix this and put it in this, watch out, original matrix is destroyed.
 * @param option = 1 positive square root, = -1 negative square root.
 */
void SUP::sqrt(int option){

   I1->pseudo_sqrt(option,Tools::gL() - 1,0);

#ifdef __Q2_CON
   Q2->pseudo_sqrt(option,Tools::gL() - 1,0);
#endif

#ifdef __I2_CON
   I2->sqrt(option);
#endif

#ifdef __Q1_CON
   Q1->sqrt(option);
#endif

#ifdef __G1_CON
   G1->sqrt(option);
#endif

#ifdef __G2_CON
   G2->sqrt(option);
#endif

}

/**
 * Multiply symmetric SUP blockmatrix object left en right with symmetric SUP blockmatrix map to 
 * form another symmetric matrix and put it in (*this): this = map*object*map
 * @param map SUP that will be multiplied to the left en to the right of matrix object
 * @param object central SUP
 */
void SUP::L_map(const SUP &map,const SUP &object){

   I1->L_map(map.gI1(),object.gI1());

#ifdef __Q2_CON
   Q2->L_map(map.gQ2(),object.gQ2());
#endif

#ifdef __I2_CON
   I2->L_map(map.gI2(),object.gI2());
#endif

#ifdef __Q1_CON
   Q1->L_map(map.gQ1(),object.gQ1());
#endif

#ifdef __G1_CON
   G1->L_map(map.gG1(),object.gG1());
#endif

#ifdef __G2_CON
   G2->L_map(map.gG2(),object.gG2());
#endif

}

/**
 * add the SUP X_p times the constant alpha to this
 * @param alpha the constant to multiply the X_p with
 * @param X_p the SUP to be multiplied by alpha and added to (*this)
 */
void SUP::daxpy(double alpha,const SUP &X_p){

   I1->daxpy(alpha,X_p.gI1());

#ifdef __Q2_CON
   Q2->daxpy(alpha,X_p.gQ2());
#endif

#ifdef __I2_CON
   I2->daxpy(alpha,X_p.gI2());
#endif

#ifdef __Q1_CON
   Q1->daxpy(alpha,X_p.gQ1());
#endif

#ifdef __G1_CON
   G1->daxpy(alpha,X_p.gG1());
#endif

#ifdef __G2_CON
   G2->daxpy(alpha,X_p.gG2());
#endif

}

/**
 * @return trace of the SUP matrix, defined as sum of the traces of the separate carrierspace matrices
 */
double SUP::trace() const{

   double ward = I1->trace();

#ifdef __Q2_CON
   ward += Q2->trace();
#endif

#ifdef __I2_CON
   ward += I2->trace();
#endif

#ifdef __Q1_CON
   ward += Q1->trace();
#endif

#ifdef __G1_CON
   ward += G1->trace();
#endif

#ifdef __G2_CON
   ward += G2->trace();
#endif

   return ward;

}

/**
 * General matrixproduct between two SUP matrices, act with Matrix::mprod on every block
 * @param A left hand matrix
 * @param B right hand matrix
 * @return The product AB
 */
SUP &SUP::mprod(const SUP &A,const SUP &B){

   I1->mprod(A.gI1(),B.gI1());

#ifdef __Q2_CON
   Q2->mprod(A.gQ2(),B.gQ2());
#endif

#ifdef __I2_CON
   I2->mprod(A.gI2(),B.gI2());
#endif

#ifdef __Q1_CON
   Q1->mprod(A.gQ1(),B.gQ1());
#endif

#ifdef __G1_CON
   G1->mprod(A.gG1(),B.gG1());
#endif

#ifdef __G2_CON
   G2->mprod(A.gG2(),B.gG2());
#endif

   return *this;

}

/**
 * Fill the SUP matrix (*this) with two dDPM matrices like: this = diag[I^l(W)  Q^l(W)]
 * @param W input dDPM
 */
void SUP::fill(const dDPM &W){

   *I1 = W;

#ifdef __Q2_CON
   Q2->Q('U',W);
#endif

#ifdef __I2_CON
   I2->I(W);
#endif

#ifdef __Q1_CON
   Q1->Q(W);
#endif

#ifdef __G1_CON
   G1->G1(W);
#endif

#ifdef __G2_CON
   G2->G2(W);
#endif

}

/**
 * Fill the SUP matrix (*this) with the matrix in I1 as input
 */
void SUP::fill(){

#ifdef __Q2_CON
   Q2->Q('U',*I1);
#endif

#ifdef __I2_CON
   I2->I(*I1);
#endif

#ifdef __Q1_CON
   Q1->Q(*I1);
#endif

#ifdef __G1_CON
   G1->G1(*I1);
#endif

#ifdef __G2_CON
   G2->G2(*I1);
#endif

}

/**
 * Seperate SUP into two SUP's, a positive and negative semidefinite part.
 * @param p positive (plus) output part
 * @param m negative (minus) output part
 */
void SUP::sep_pm(SUP &p,SUP &m){

   I1->sep_pm(p.gI1(),m.gI1());

#ifdef __Q2_CON
   Q2->sep_pm(p.gQ2(),m.gQ2());
#endif

#ifdef __I2_CON
   I2->sep_pm(p.gI2(),m.gI2());
#endif

#ifdef __Q1_CON
   Q1->sep_pm(p.gQ1(),m.gQ1());
#endif

#ifdef __G1_CON
   G1->sep_pm(p.gG1(),m.gG1());
#endif

#ifdef __G2_CON
   G2->sep_pm(p.gG2(),m.gG2());
#endif

}
