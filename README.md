# simple_kalman

## Kalman filter implementation in C++ 

This work is motivated by Lecture #7 of the [this]: http://www.idsc.ethz.ch/education/lectures/recursive-estimation.html course.

The Kalman filter implementation is located in: simple_kalman.h and simple_kalman.cpp while the extended Kalman filter implementation is located in ekf.h and ekf.cpp.

Please use cmake to build the codes by first modifying the main file name in CMakeList.txt.

The steps to compile are:

```
cd build
cmake ..
make
```

## Dependencies

Please install Armadillo (http://arma.sourceforge.net), take it from your distribution repository. Armadillo is very easy to use. Have a look on this [link]: http://arma.sourceforge.net/docs.html for quick tutorial on using Armadillo.

