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
  
  /*
   * Define a system: x = Ax + Bu + v
   *                  z = Hx + w
   * where: process noise v ~ N(0,Q)
   *        measurement noise w ~ N(0,R)
   */
  void InitSystem (const mat& A, const mat& B, const mat& H, const mat& Q, const mat& R);
  
  /*
   * Do Kalman iteration step-by-step
   * x is the true states
   * x_m is the estimates states
   * u is the applied input to the system
   */
  void Kalmanf(colvec& x, colvec& x_m, const colvec& u);
private:
  
  /// Variables to hold the system parameters
  mat A_;
  mat B_;
  mat H_;
  mat Q_;
  mat R_;
  colvec v_;
  colvec w_;
  
  mat sqrt_Q_;
  mat sqrt_R_;
  
  /// System states and outputs
  colvec x_;
  colvec z_;
  
 
  /*
   * Variables for the Kalman
   * subscript m means measurement
   * subscript p means priori
   */
  colvec x_m_;
  colvec x_p_;
  
  mat P_p_;
  mat P_m_;
  
};

#endif