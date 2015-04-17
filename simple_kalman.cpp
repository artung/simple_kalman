#include "simple_kalman.h"

SimpleKalman::SimpleKalman()
{

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
    
  sqrt_Q_ = sqrt(Q_);
  sqrt_R_ = chol(R_);
        
  int n_states = A.n_cols;
  int n_outputs = H.n_rows;
  
  assert(A.is_square() && "Whoops, A must be a square matrix");
  assert(B.n_rows == A.n_rows && "Whoops, B has wrong dimension");
  assert(Q.is_square() && "Whoops, Q must be a square matrix");
  assert(R.n_rows == n_outputs && "Whoops, R must be a row vector with element numbers equal to number of the outputs");
  
  x_.resize(n_states);
  x_.zeros();
  
  x_p_.resize(n_states);
  x_m_.resize(n_states);
  
  P_p_.resize(n_states, n_states);
  P_m_.resize(n_states, n_states);
  
  P_m_.eye();
  x_m_.zeros();
  
  v_.resize(n_states);
  w_.resize(n_outputs);
}

void SimpleKalman::Kalmanf(colvec &x, colvec& x_m, const colvec& u)
{
  // True system:
  v_.randn(A_.n_rows);
  w_.randn(H_.n_rows);
  v_ = sqrt_Q_ * v_;
  w_ = sqrt_R_ * w_;
  x_ = A_ * x_ + B_ * u + v_;
  z_ = H_ * x_ + w_;
  
  // Prior update:
  x_p_ = A_ * x_m_ + B_ * u;
  P_p_ = A_ * P_m_ * trans(A_) + Q_;
  
  // Measurement update
  mat K = P_p_ * trans(H_) * inv(H_ * P_p_ * trans(H_) + R_);
  x_m_ = x_p_ + K * (z_ - H_ * x_p_);
  P_m_ = P_p_ - K * H_ * P_p_;

  // Return the results
  x = x_;
  x_m = x_m_;
}

