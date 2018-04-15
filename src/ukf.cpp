#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.001

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;
using namespace std;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // If this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // If this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // Initial state vector
  x_ = VectorXd(5);

  // Initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  // Can be tuned to improve NIS
  std_a_ = 1.0;

  // Process noise standard deviation yaw acceleration in rad/s^2
  // Can be tuned to improve NIS
  std_yawdd_ = 0.5;
  
  //DO NOT MODIFY measurement noise values below these are provided by the sensor manufacturer.
  // Laser measurement noise standard deviation position1 in m
  std_laspx_ = 0.15;

  // Laser measurement noise standard deviation position2 in m
  std_laspy_ = 0.15;

  // Radar measurement noise standard deviation radius in m
  std_radr_ = 0.3;

  // Radar measurement noise standard deviation angle in rad
  std_radphi_ = 0.03;

  // Radar measurement noise standard deviation radius change in m/s
  std_radrd_ = 0.3;
  //DO NOT MODIFY measurement noise values above these are provided by the sensor manufacturer.

  // Not initialized until first measurement
  is_initialized_ = false;

  // Set dimension of state vector
  n_x_ = 5;

  // Set dimension of augmented vector
  n_aug_ = 7;

  // Set spread parameter
  lambda_ = 3 - n_x_;

  // Initial matrix for sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);

  // Initial time
  time_us_ = 0;

  // Initialize vector for weights
  weights_ = VectorXd(2 * n_aug_ + 1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < weights_.size(); i++)
  {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  // Initial matrix for radar noise 
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_, 0, 0,
              0, std_radphi_*std_radphi_, 0,
              0, 0, std_radrd_*std_radrd_;

  // Initial matrix for laser noise
  R_laser_ = MatrixXd(2,2);
  R_laser_ << std_laspx_*std_laspx_, 0,
              0, std_laspy_*std_laspy_;

  // Set NIS for laser
  nis_laser_ = 0;

  // Set NIS for radar
  nis_radar_ = 0;

}

UKF::~UKF() {}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {

  /********************************************************************
   * Initialization
  ********************************************************************/
  if(!is_initialized_)
  {
    if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      // Initialize state
      x_ << meas_package.raw_measurements_(0), meas_package.raw_measurements_(1), 0.0, 0.0, 0.0;
    } 
    else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      double rho = meas_package.raw_measurements_(0);
      double phi = meas_package.raw_measurements_(1);
      double rhod = meas_package.raw_measurements_(2);

      // Initialize state
      x_ << rho * cos(phi), rho * sin(phi), 0.0, 0.0, 0.0;
    }

    // Check whether values for p_x and p_y are too small
    if (x_(0) < EPS)
    {
      x_(0) = EPS;
    } 

    if (x_(1) < EPS)
    {
      x_(1) = EPS;
    }

    // Initialize covariance matrix
    P_ << 1, 0, 0, 0, 0,
          0, 1, 0, 0, 0,
          0, 0, 1, 0, 0,
          0, 0, 0, 1, 0,
          0, 0, 0, 0, 1;

    // Set time to timespamp of the current measurement
    time_us_ = meas_package.timestamp_;

    // Done initializing
    is_initialized_ = true;

    return;
  }

  /********************************************************************
   * Prediction
  ********************************************************************/

  // Calculate time delta between previous and current measurement
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;
cout << "x (prä) = " << x_ << endl;
  Prediction(dt);

  /********************************************************************
   * Update
  ********************************************************************/

  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    UpdateLidar(meas_package);
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    UpdateRadar(meas_package);
  }
  

  cout << "x = " << x_ << endl;
  //cout << "P = " << P_ << endl;

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {

  /********************************************************************
   * Generate Augmented Sigma Points
  ********************************************************************/
  // Augmented state vector
  VectorXd x_aug = VectorXd(n_aug_);
  // nu_a and nu_phi are zero as the mean of the noise is zero
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  // Augmented covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;
  // Square root of augmented state covariance
  MatrixXd A = P_aug.llt().matrixL();

  // Augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);

  // Calculate equal part of sigma point equation
  MatrixXd A_equal = sqrt(lambda_ + n_aug_) * A;
  // Calculate augmented sigma points
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; i++)
  {
    Xsig_aug.col(i+1) = x_aug + A_equal.col(i);
    Xsig_aug.col(i+1+n_aug_) = x_aug - A_equal.col(i);
  }

  /********************************************************************
   * Predict Sigma Points
  ********************************************************************/
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    // Extract values for readibility
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // Precalculate values for optimization
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);
    double dt_square = delta_t * delta_t;

    // Predict state values
    double px_p, py_p;
    if (fabs(yawd) > EPS)
    {
      px_p = p_x + v / yawd * (sin(yaw + yawd * delta_t) - sin_yaw);
      py_p = p_y + v / yawd * (cos_yaw - cos(yaw + yawd * delta_t));
    } else {
      px_p = p_x + v * cos_yaw * delta_t;
      py_p = p_y + v * sin_yaw * delta_t;
    }

    double v_p = v;
    double yaw_p = yaw + yawd * delta_t;
    double yawd_p = yawd;

    // Add noise
    px_p += 0.5 * nu_a * dt_square * cos_yaw;
    py_p += 0.5 * nu_a * dt_square * sin_yaw;
    v_p += delta_t * nu_a;
    yaw_p += 0.5 * nu_yawdd * dt_square;
    yawd_p += delta_t * nu_yawdd;

    // Write predicted sigma point into matrix
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  /********************************************************************
   * Predicted State & Covariance
  ********************************************************************/

  // Mean of the predicted sigma points
  x_ = Xsig_pred_ * weights_;

  // Predicted state covariance
  P_.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;

    // Normalize the angle
    NormalizeAngle(&x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {

  // Measurement space dimension
  int n_z = 2;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = Xsig_pred_.block(0, 0, n_z, 2 * n_aug_ + 1);

  // Update the state
  UpdateCommon(meas_package, Zsig, n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {

  // Measurement space dimension
  int n_z = 3;

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_z, 2 * n_aug_ + 1);

  // Transform sigma points into measurement space
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    const double px = max(Xsig_pred_(0,i), EPS);
    const double py = max(Xsig_pred_(1,i), EPS);
    const double v = Xsig_pred_(2,i);
    const double yaw = Xsig_pred_(3,i);
    const double yawd = Xsig_pred_(4,i);

    const double r = max(sqrt(px*px + py*py), EPS);

    Zsig(0,i) = r;
    Zsig(1,i) = atan2(py, px);
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / r;
  }

  // Update the state
  UpdateCommon(meas_package, Zsig, n_z);
}

/**
 * Normalizes the angle to the range of -pi to pi
 * @param {double} *angle pointer to the angle that should be normalized
 */ 
void UKF::NormalizeAngle(double *angle) {
  
  while (*angle > M_PI) *angle -= 2.0 * M_PI;
  while (*angle < -M_PI) *angle += 2.0 * M_PI;
}

/**
 * Common part of the lidar and radar update function
 * @param {MeasurementPackage} meas_package The measurement at k+1
 * @param {MatrixXd} Zsig Sigma points in measurement space
 * @param {int} n_z Dimension of the measurement space
 */
void UKF::UpdateCommon(MeasurementPackage meas_package, MatrixXd Zsig, int n_z) {

  // Extract the current measurement
  VectorXd z = meas_package.raw_measurements_;

  // Calculate the mean of the predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred = Zsig * weights_;

  // Calculate innovation covariance matrix
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    VectorXd z_diff = Zsig.col(i) - z_pred;
    S += weights_(i) * z_diff * z_diff.transpose();
  }

  // Add the measurement noise uncertainty
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    S += R_laser_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    S += R_radar_;
  }

  // Calculate cross correlation function
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);

  for (int i = 0; i < 2 * n_aug_ + 1; i++)
  {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    // Normalize the angle
    NormalizeAngle(&z_diff(1));

    // State difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    
    // Normalize the angle
    NormalizeAngle(&x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  // Kalman gain
  MatrixXd K = Tc * S.inverse();

  // Residual
  VectorXd z_diff = z - z_pred;

  // Normalize the angle
  NormalizeAngle(&z_diff(1));

  // Update state mean and covariance
  x_ += K * z_diff;
  P_ -= K * S * K.transpose();

  // Calculate Normalized Innovation Squared (NIS)
  if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  {
    nis_laser_ = z_diff.transpose() * S.inverse() * z_diff;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    nis_radar_ = z_diff.transpose() * S.inverse() * z_diff;
  }
}
