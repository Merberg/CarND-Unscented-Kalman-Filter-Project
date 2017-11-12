#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  VectorXd rmse(4);
  rmse << 0, 0, 0, 0;
  unsigned int e_size = estimations.size();

  //Parameter checking
  if (e_size == 0 || e_size != ground_truth.size()) {
    cout << "Tools::CalculateRMSE - parameter error" << endl;
    return rmse;
  }

  //accumulate squared residuals
  for (unsigned int i = 0; i < e_size; ++i) {

    VectorXd residual = estimations[i] - ground_truth[i];

    //coefficient-wise multiplication
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  //calculate the mean
  rmse = rmse / e_size;

  //calculate the squared root
  rmse = rmse.array().sqrt();

  cout << "rmse: " << rmse << endl;
  cout << endl;

  return rmse;
}
