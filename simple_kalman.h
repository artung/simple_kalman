#ifndef SIMPLE_KALMAN_H
#define SIMPLE_KALMAN_H

#include <assert.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class SimpleKalman {
public:
  SimpleKalman();
  ~SimpleKalman();
  void InitSystem (const mat& A, const mat& B, const mat& H, const mat& Q, const mat& R);
  void Kalmanf(mat& x, mat& x_m, const mat& u);
private:
  
  // System
  mat A_;
  mat B_;
  mat H_;
  mat Q_;
  mat R_;
  mat x_;
  mat z_;
  mat v_;
  mat w_;
 
  // Variables for the Kalman
  mat x_m_;
  mat x_p_;
  
  mat P_p_;
  mat P_m_;
  
  
  int nIter_;
};

#endif