/**
 * @file ekf.h
 * @author Auralius Manurung
 * @date 18 Apr 2015
 * @brief Header file for the Extended Kalman filter implementation.
 * 
 * @section DESCRIPTION
 * Define a non-linear discrete-time process: 
 * \f[x_k = f(x_{k-1}, u_{k-1}, k) + v_{k-1}\f]
 * \f[z_k = h(x_k, u_k, k) + w_k\f]
 * where:\n 
 * \f$f\f$ is the dynamic model of the system\n
 * \f$h\f$ is the measurement model of the system\n
 * \f$v\f$ is the process noise (Gaussian with covariance Q)\n  
 * \f$w\f$ is the measurement noise (Gaussian with covariance R)\n 
 * \f$x\f$ is the state vector\n
 * \f$z\f$ is the output vector\n
 * \f$u\f$ is the input vector\n
 */

#ifndef EKF_H
#define EKF_H

#include <assert.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

/*!
 * @brief Implemetation of the extended Kalman filter. 
 * This class needs to be derived.
 */
class EKF {
public:  
  /*!
   * \brief Constructor, nothing happens here.
   */
  EKF();
  
  /*!
   * \brief Destructur, nothing happens here.
   */
  ~EKF();
  
  /*!
   * \brief Tell me how many states and outputs you have!
   * @param n_states Number of the elements on the input vector x
   * @param n_output Number of the elements on the otput vector z
   * @param Q Process noise covariance
   * @param R Measurement noise covariance
   */
  void InitSystem(int n_states, int n_outputs, const mat& Q, const mat& R);
  /*!
   * \brief Define model of your system.
   * @param x System states
   * @param u System inputs
   * @param k k-th iteration
   */
  virtual colvec f(const colvec &x, const colvec &u, const int k);
  
  /*!
   * \brief Define the output model of your system.
   * @param x System states
   * @param u Input vector
   * @param k k-th iteration
   */
  virtual colvec h(const colvec &x, const colvec &u, const int k);
  
   
  /*!
   * \brief Initialize the system states.
   * Must be called after InitSystem.
   * If not called, covariance state is Initialized to an identity matrix.
   * @param x0 Inital value for the system state
   * @param k k-th iteration
   */
  void InitSystemState(const colvec& x0);
 
  /*!
   * \brief Do the extended Kalman iteration step-by-step.
   * @return x the true states, this is a returned value
   * @return x_m the estimated states, this is a returned value
   * @param u the applied input to the system
   * @param k k-th iteration
   */
  void EKalmanf(colvec& x, colvec& x_m, const colvec& u, const int k);
  
private:
  /*!
   * \brief Compute the Jacobian of f numerically using  a  small  
   * finite-difference perturbation magnitude. 
   * @param x System states
   * @param u Input vector
   * @param k k-th iteration
   */
  void CalcF(const colvec &x, const colvec &u, const int k);
  
  /*!
   * \brief Compute the Jacobian of h numerically using  a  small  
   * finite-difference perturbation magnitude.
   * @param x System states
   * @param u Input vector
   * @param k k-th iteration 
   */
  void CalcH(const colvec &x, const colvec &u, const int k);
  
 ;

  mat F_;          ///< Jacobian of F	
  mat H_;          ///< Jacobian of H
  mat Q_;          ///< Process noise covariance
  mat R_;          ///< Measurement noise covariance
  colvec v_;       ///< Gaussian process noise
  colvec w_;       ///< Gaussian measurement noise
  
  mat sqrt_Q_;     ///< Process noise stdev
  mat sqrt_R_;     ///< Measurement noise stdev
  

 
  colvec x_m_;     ///< State vector after measurement update
  colvec x_p_;     ///< State vector after a priori update
  
  mat P_p_;        ///< State covariance after a priori update
  mat P_m_;        ///< State covariance after measurement update
  
  double epsilon_; ///< Very small number
  
protected:
  
  int nStates_;   ///< Number of the states
  int nOutputs_;  ///< Number of outputs
  
  colvec x_;      ///< State vector
  colvec z_;      ///< Output matrix
};


#endif
