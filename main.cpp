#include <iostream>
#include <fstream>

#include "simple_kalman.h"


int main(int argc, char** argv)
  {
    /* 
     * Log the result into a tab delimitted file, later we can open 
     * it with Matlab
     */
    ofstream log_file;
    log_file.open("log_file.txt");
   
    /*
     * Define the system 
     * Sinusoid voltage is measured as an output of a system
     */
    int nStates = 1;
    mat A(1,1), B(1,1), H(1,1), Q(1,1), R(1,1);
    
    A(0, 0) = 0;
    B(0, 0) = 1;
    Q(0, 0) = 2;
    H(0, 0) = 1;
    R(0, 0) = 2;
    
    SimpleKalman kalman;
    kalman.InitSystem(A, B, H, Q, R);
    
    rowvec z(1);
    rowvec x(1), x_m(1);
    rowvec u(1);
  
    double t = 0;

    do {
      u(0,0) = 12 * sin(2 * M_PI * 0.1 * t);
      kalman.Kalmanf(x, x_m, u);
      log_file << t << '\t' << x.at(0,0) << '\t' << x_m.at(0,0) << '\t' << u.at(0,0) << '\n';
      t = t + 0.1;
    } while (t < 20);
    
    log_file.close();
    
    return 0;
  }