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
   * @param A System matrix
   * @param B Input matrix
   * @param H Output matrix
   * @param Q Process noise covariance
   * @param R Measurement noise covariance
   */
  void InitSystem (const mat& A, const mat& B, const mat& H, const mat& Q, const mat& R);
  
  /*
   * Do Kalman iteration step-by-step
   * @return x the true states, this is a returned value
   * @return x_m the estimated states, this is a returned value
   * @param u the applied input to the system
   */
  void Kalmanf(colvec& x, colvec& x_m, const colvec& u);
  
private:
  
  mat A_;      ///< System matrix
  mat B_;      ///< Input matrix
  mat H_;      ///< Output matrix
  mat Q_;      ///< Process noise covariance
  mat R_;      ///< Measurement noise covariance
  colvec v_;   ///< Gaussian process noise
  colvec w_;   ///< Gaussian measurement noise
  
  mat sqrt_Q_; ///< Process noise stdev
  mat sqrt_R_; ///< Measurement noise stdev
  
  colvec x_;   ///< State vector
  colvec z_;   ///< Output matrix
 
  colvec x_m_; ///< State vector after measurement update
  colvec x_p_; ///< State vector after a priori update
  
  mat P_p_;    ///< State covariance after a priori update
  mat P_m_;    ///< State covariance after measurement update
};

#endif