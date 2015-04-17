#include "simple_kalman.h"

SimpleKalman::SimpleKalman()
{
  nIter_ = 0;
}

SimpleKalman::~SimpleKalman()
{
  
}

void SimpleKalman::InitSystem(const mat& A, const mat& B, const mat& H, const mat& Q, const mat& R)
{
  A_ = A;
  B_ = B;
  H_ = H;
  Q_ = Q;
  R_ = R;

  x_.set_size(A_.n_rows, 1);
  x_.zeros();
  
  x_p_.set_size(A_.n_rows, 1);
  x_m_.set_size(A_.n_rows, 1);
  
  P_p_.set_size(A_.n_rows, A_.n_cols);
  P_m_.set_size(A_.n_rows, A_.n_cols);
  
  P_m_.eye();
  x_m_.zeros();
  
}

void SimpleKalman::Kalmanf(mat &x, mat& x_m, const mat& u)
{
  // True system:
  v_.randn(A_.n_rows, 1);
  w_.randn(H_.n_rows, 1);
  v_ = Q_ * v_;
  w_ = R_ * w_;
  x_ = A_ * x_ + B_ * u + v_;
  z_ = H_ * x_ + w_;
  
  // Prior update:
  x_p_ = A_ * x_m_ + B_ * u;
  P_p_ = A_ * P_m_ * trans(A_) + Q_;
  
  // Measurement update
  mat K = P_p_ * trans(H_) * inv(H_ * P_p_ * trans(H_) + R_);
  x_m = x_p_ + K * (z_ - H_ * x_p_);
  P_m_ = P_p_ - K * H_ * P_p_;

  x = x_;
}

