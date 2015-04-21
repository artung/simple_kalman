/**
 * @file ekf.h
 * @author Auralius Manurung
 * @date 18 Apr 2015
 * @brief Header file for the Extended Kalman filter implementation.
 * 
 * @section DESCRIPTION
 * Define a non-linear discrete-time process: 
 * \f[x_k = f(x_(k-1), u_k, k) + v_(k-1)\f]
 * \f[z_k = h(x_k, u_k, k) + w_k\f]
 * where:\n 
 * \f$f\f$ is the dynamic model of the system\n
 * \f$h\f$ is the measurement model of the system\n
 * \f$v ~ N(0,Q)\f$ is the process noise\n  
 * \f$w ~ N(0,R)\f$ is the measurement noise\n
 */

#ifndef EKF_H
#define EKF_H

#include <assert.h>
#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

class EKF {
public:  
  /*!
   * \brief Tell me how many states and outputs you have!
   * @param n_states Number of the elements on the input vector x
   * @param n_output Number of the elements on the otput vector z
   */
  EKF(int n_states, int n_outputs, const mat& Q, const mat& R);
  
  /*!
   * \brief Destructur, nothing happens here.
   */
  ~EKF();
  
  /*!
   * \brief Define model of your system.
   */
  virtual colvec f(const colvec &x, const colvec &u, const int k);
  
  /*!
   * \brief Define the output model of your system.
   */
  virtual colvec h(const colvec &x, const colvec &u, const int k);
  
   
  /*!
   * \brief Initialize the system states.
   */
  void InitSystemState(const colvec& x0);
 
  /*!
   * \brief Do the extended Kalman iteration step-by-step.
   * @return x the true states, this is a returned value
   * @return x_m the estimated states, this is a returned value
   * @param u the applied input to the system
   */
  void EKalmanf(colvec& x, colvec& x_m, const colvec& u, const int k);
  
private:
  /*!
   * \brief Compute the Jacobian of f numerically using  a  small  
   * finite-difference perturbation magnitude. 
   */
  void CalcF(const colvec &x, const colvec &u, const int k);
  
  /*!
   * \brief Compute the Jacobian of h numerically using  a  small  
   * finite-difference perturbation magnitude. 
   */
  void CalcH(const colvec &x, const colvec &u, const int k);
  
  int nStates_;
  int nOutputs_;

  mat F_;
  mat H_;
  mat Q_;      ///< Process noise covariance
  mat R_;      ///< Measurement noise covariance
  colvec v_;   ///< Gaussian process noise
  colvec w_;   ///< Gaussian measurement noise
  
  mat sqrt_Q_; ///< Process noise stdev
  mat sqrt_R_; ///< Measurement noise stdev
  

 
  colvec x_m_; ///< State vector after measurement update
  colvec x_p_; ///< State vector after a priori update
  
  mat P_p_;    ///< State covariance after a priori update
  mat P_m_;    ///< State covariance after measurement update
  
  double epsilon_;
  
protected:
    
  colvec f_;
  colvec h_;
  
  colvec x_;   ///< State vector
  colvec z_;   ///< Output matrix
};


#endif
