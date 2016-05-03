/**
 * @file main4.cpp
 * @author Auralius Manurung
 * @date 18 Apr 2015
 * @brief Example for the Extended Kalman FIlter Implementation.
 * 
 * The following example is taken from:\n
 * http://ch.mathworks.com/matlabcentral/fileexchange/38302-kalman-filter-package/content//Kalman%20Filter%20Package/Examples/ExtendedKalmanFilterDemo.m
 */

#include <iostream>
#include <fstream>

#include "ekf.h"

/*!
 * @brief Class EKF needs to be derived, two virtual functions are provided in 
 * which system model and output model are described.
 */
class MyEKF: public EKF
{
public:  
  virtual colvec f(const colvec& x, const colvec& u, const int k) {
    colvec xk(nOutputs_);
    xk(0) = sin(x(1) * k);
    xk(1) = x(1);
    return xk;
  }
  
  virtual colvec h(const colvec& x, const colvec& u, const int k) {
    colvec zk(nOutputs_);
    zk(0) = x(0);
    zk(1) = x(1);
    return zk;
  }
};

/////////////////////////////////////////////////////////////////////////////

int main(int argc, char** argv)
{
  /* 
   * Log the result into a tab delimitted file, later we can open 
   * it with Matlab. Use: plot_data4.m to plot the results.
   */
  ofstream log_file;
  log_file.open("log_file.txt");
  
  int n_states = 2;
  int n_outputs = 2;
  mat Q(2, 2);
  mat R(2, 2);
  
  Q << 0.001 << 0    << endr
    << 0     << 0    << endr;
    
  R << 0.1   << 0    << endr
    <<   0   << 0.01 << endr;
  
  colvec x0(3);
  x0 << 0 << 3 * M_PI / 500;
  
  colvec x, x_m;
  colvec u;
  
  MyEKF myekf;
  myekf.InitSystem(n_states, n_outputs, Q, R);
  myekf.InitSystemState(x0);
  
  for (int k = 0; k < 1000; k ++) {
    myekf.EKalmanf(x, x_m, u, k);
    log_file << k << '\t' << x(0) << '\t' << x_m(0) << '\n'; 
  }
  
  log_file.close();
  
  return 0;
}