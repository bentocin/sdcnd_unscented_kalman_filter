#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

Tools::Tools() {}

Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Calculate the RMSE here.
  VectorXd rmse(4);
  rmse.fill(0.0);

  // Sanity checks before calculation
  if (estimations.size() != ground_truth.size())
  {
    std::cout << "ERROR: Estimations and ground truth have different length" << endl;
    return rmse;
  }

  if (estimations.size() == 0)
  {
    std::cout << "ERROR: No values in estimations" << endl;
    return rmse;
  }

  // Calculate sum of squared residuals
  for (int i = 0; i < estimations.size(); i++)
  {
    // Calculate error and add square to the sum
    VectorXd residual = estimations[i] - ground_truth[i];
    residual = residual.array() * residual.array();
    rmse += residual;
  }

  // Average and root of the squared errors
  rmse /= estimations.size();
  rmse = rmse.array().sqrt();

  return rmse;

}