#include "ekf.h"

EKF::EKF(int n_states, int n_outputs)
{
  nStates_ = n_states;
  nOutputs_ = n_outputs;

  x_.resize(n_states);
  F_.resize(n_states, n_states);
  H_.resize(n_outputs, n_states);
  
  Q_.resize(n_outputs, n_outputs);
  R_.resize(n_outputs, n_outputs);
}


EKF::~EKF()
{

}

void EKF::ModelF()
{

}

void EKF::ModelH()
{

}

void EKF::CalcJaccF()
{

}

void EKF::CalcJaccH()
{

}

void EKF::Kalmanf(colvec& x, colvec& x_m, const colvec& u)
{

}
