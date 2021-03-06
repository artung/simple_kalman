/**
 * @file main1.cpp
 * @author Auralius Manurung
 * @date 18 Apr 2015
 * @brief Example for the Kalman Filter Implementation.
 */

#include <iostream>
#include <fstream>

#include "simple_kalman.h"


/*!
 * Example case is taken from:\n
 * http://www.mathworks.com/matlabcentral/fileexchange/18465-learning-the-kalman-filter-in-simulink-v2-1/content/html/runkalmanfilter.html\n
 * Constant voltage of 12 V is measured as the output of a system.
 */
int main(int argc, char** argv)
  {
    /* 
     * Log the result into a tab delimitted file, later we can open 
     * it with Matlab. Use: plot_data1.m to plot the results.
     */
    ofstream log_file;
    log_file.open("log_file.txt");
   
    /*
     * Define the system
     */
    mat A(1,1), B(1,1), H(1,1), Q(1,1), R(1,1);
    
    A(0, 0) = 0;
    B(0, 0) = 1;
    Q(0, 0) = 2;
    H(0, 0) = 1;
    R(0, 0) = 2;
    
    SimpleKalman kalman;
    kalman.InitSystem(A, B, H, Q, R);
    
    colvec x(1), x_m(1);
    colvec u(1);
  
    for (int k = 0; k < 100; k++) {
      u(0,0) = 12;
      kalman.Kalmanf(x, x_m, u);
      log_file << k << '\t' << x.at(0,0) << '\t' << x_m.at(0,0) << '\t' << u.at(0,0) << '\n';      
    }
    
    log_file.close();
    
    return 0;
  }