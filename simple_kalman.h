/**
 * @file simple_kalman.h
 * @author Auralius Manurung
 * @date 18 Apr 2015
 * @brief Header file for the Kalman filter implementation.
 * 
 * @section DESCRIPTION
 * Define a system:
 * \f[x_k = Ax_{k-1} + Bu_{k-1} + v_{k-1}\f]
 * \f[z_k = Hx_k + w_k\f]
 * where:\n
 * \f$v\f$ is the process noise (Gaussian with covariance Q)\n  
 * \f$w\f$ is the measurement noise (Gaussian with covariance R)\n 
 * \f$A\f$ is the system matrix\n
 * \f$B\f$ is the input matrix\n
 * \f$H\f$ is the output matrix\n
 * \f$x\f$ is the state vector\n
 * \f$z\f$ is the output vector\n
 * \f$u\f$ is the input vector\n
 */

#ifndef SIMPLE_KALMAN_H
#define SIMPLE_KALMAN_H

#include <assert.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

/*!
 * @brief Implemetation of the Kalman filter. 
 */
class SimpleKalman {
public:
  SimpleKalman();
  ~SimpleKalman();
  
  /*!
   * @brief Define the system.
   * @param A System matrix
   * @param B Input matrix
   * @param H Output matrix
   * @param Q Process noise covariance
   * @param R Measurement noise covariance
   */
  void InitSystem (const mat& A, const mat& B, const mat& H, const mat& Q, const mat& R);
  
  /*!
   * @brief Initialize the system states.
   * Must be called after InitSystem.
   * If not, called, system states are initialized to zero.
   * @param x0 Inital value for the system state
   */
  void InitSystemState(const colvec& x0);
  
  /*!
   * @brief Initialize the state covariance.
   * Must be called after InitSystem.
   * If not called, covariance state is Initialized to an identity matrix.
   * @param P0 Inital value for the stat covariance
   */
  void InitStateCovariance(const mat& P0);
  
  /*!
   * @brief Do Kalman iteration step-by-step
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