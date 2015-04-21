/**
 * @file main2.cpp
 * @author Auralius Manurung
 * @date 18 Apr 2015
 * @brief Example for the Kalman Filter Implementation.
 */

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
     */
    mat A(2,2), B(2,1), H(1,1), Q(2,2), R(1,1);
    
    A << 1 << 1 << endr
      << 0 << 1 << endr;
    
    B << 0 << endr
      << 1 << endr;
    
    H << 1 << 0;
    
    Q << 0 << 0 << endr
      << 0 << 0.1 << endr;
      
    R << 0.5;
      
    SimpleKalman kalman;
    kalman.InitSystem(A, B, H, Q, R);
    
    colvec z(2);
    colvec x(2), x_m(2);
    colvec u(1);
  
    for (int i = 0; i < 30 ; i ++) {
      if (i < 10)
	u << 1;
      else if (i >= 10 && i < 20)
	u << -1;
      else 
	u << 0;
      
      kalman.Kalmanf(x, x_m, u);
      log_file << i << '\t' << x(0) << '\t' << x(1) << '\t' << x_m(0) << '\t' << x_m(1) << '\t' << u(0) << '\n';
      
    }
    
    log_file.close();
    
    return 0;
  }