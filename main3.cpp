#include <iostream>
#include <fstream>

#include "simple_kalman.h"


int main(int argc, char** argv)
  {
    /* 
     * Log the result into a tab delimitted file, later we can open 
     * it with Matlab. Use: plot_data2.m to plot the results.
     */
    ofstream log_file;
    log_file.open("log_file.txt");
   
    /*
     * Define kinematic system, with position and
     * velocity as the states
     * 4 input states
     * 2 output state
     */
    mat A(4,4), B(4,1), H(4,2), Q(4,4), R(2,2);
    
    A << 1 << 0 << 1 << 0 << endr
      << 0 << 1 << 0 << 1 << endr
      << 0 << 0 << 1 << 0 << endr
      << 0 << 0 << 0 << 1 << endr;
    
    B = B.zeros();
    
    H << 1 << 0 << 0 << 0 << endr
      << 0 << 1 << 0 << 0 << endr;
    
    Q = Q.eye();
      
    R = 10 * R.eye();
    
    colvec x0(4);
    x0 << 10 << 10 << 1 << 0;
    
    mat P0(4,4);
    P0 = 10 * P0.eye();
      
    SimpleKalman kalman;
    kalman.InitSystem(A, B, H, Q, R);
    kalman.InitSystemState(x0);
    kalman.InitStateCovariance(P0);
    
    colvec z(4);
    colvec x(4), x_m(4);
    colvec u(1);
    
    // No inputs, system subjects only to random perturbation
    u = u.zeros(); 
  
    for (int i = 0; i < 100 ; i ++) {  
      kalman.Kalmanf(x, x_m, u);
      log_file << x(0) << '\t' << x(1) << '\t' << x(2) << '\t' << x(3) << '\t' 
               << x_m(0) << '\t' << x_m(1) << '\t' << x_m(2) << '\t' << x_m(3) << '\t' << '\n';
    }
    
    log_file.close();
    
    return 0;
  }