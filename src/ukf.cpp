#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

//******************************************************************************
UKF::UKF()
    : is_initialized_(false),
      use_laser_(true),
      use_radar_(true),
      time_us_(0),
      std_a_(8),  //TODO adjust
      std_yawdd_(8),  //TODO adjust
      std_laspx_(0.15),
      std_laspy_(0.15),
      std_radr_(0.3),
      std_radphi_(0.03),
      std_radrd_(0.3),
      n_x_(5),
      n_aug_(7),
      n_sig_(2 * n_aug_ + 1),
      lambda_(3 - n_aug_) {
  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  weights_ = VectorXd(n_sig_);
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);
}

UKF::~UKF() {
}

//******************************************************************************
VectorXd UKF::ProcessMeasurement(const MeasurementPackage &measurement) {

  //Split depending on the measurement type
  if (measurement.sensor_type_ == MeasurementPackage::LASER) {
    ProcessLidarMeasurement(measurement);
  } else {
    ProcessRadarMeasurement(measurement);
  }

  //Return the updated state
  return x_;
}

//******************************************************************************
void UKF::ProcessLidarMeasurement(const MeasurementPackage &measurement) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    float px = measurement.raw_measurements_(0);
    float py = measurement.raw_measurements_(1);
    Init(px, py, measurement.timestamp_);

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  Prediction(measurement.timestamp_);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  UpdateLidar(measurement);

}

//******************************************************************************
void UKF::ProcessRadarMeasurement(const MeasurementPackage &measurement) {

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    float ro = measurement.raw_measurements_(0);
    float theta = measurement.raw_measurements_(1);
    Init(ro * cos(theta), ro * sin(theta), measurement.timestamp_);

    is_initialized_ = true;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/
  Prediction(measurement.timestamp_);

  /*****************************************************************************
   *  Update
   ****************************************************************************/
  UpdateRadar(measurement);
}

//******************************************************************************
void UKF::Init(float px, float py, long long ts) {
  /**
   Initialize state.
   */
  x_ = VectorXd(n_x_);
  x_.fill(0.0);
  x_(0) = px;
  x_(1) = py;

  // record the time
  time_us_ = ts;
}

//******************************************************************************
void UKF::Prediction(double timestamp) {

  double delta_t = (timestamp - time_us_);
  time_us_ = timestamp;

  /*
   * Generate Sigma Points
   */
  //create augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug.head(5) = x_;
  x_aug(5) = 0;
  x_aug(6) = 0;

  //create augmented state covariance
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug.fill(0.0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(n_x_, n_x_) = std_a_ * std_a_;
  P_aug(n_x_ + 1, n_x_ + 1) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);
  Xsig_aug.col(0) = x_aug;

  for (int i = 0; i < n_aug_; i++) {
    Xsig_aug.col(i + 1) = x_aug + sqrt(lambda_ + n_aug_) * L.col(i);
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt(lambda_ + n_aug_) * L.col(i);
  }

  /*
   * Predict Sigma Points
   */
  double half_delta_t_sqr = 0.5 * delta_t * delta_t;

  for (int i = 0; i < n_sig_; i++) {
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    if (yawd != 0) {
      Xsig_pred_(0, i) = p_x
          + v / yawd * (sin(yaw + yawd * delta_t) - sin(yaw));
      Xsig_pred_(1, i) = p_y
          + v / yawd * (-cos(yaw + yawd * delta_t) + cos(yaw));
    } else {
      Xsig_pred_(0, i) = p_x + v * cos(yaw) * delta_t;
      Xsig_pred_(1, i) = p_y + v * sin(yaw) * delta_t;
    }
    Xsig_pred_(0, i) += half_delta_t_sqr * cos(yaw) * nu_a;
    Xsig_pred_(1, i) += half_delta_t_sqr * sin(yaw) * nu_a;
    Xsig_pred_(2, i) = v + delta_t * nu_a;
    Xsig_pred_(3, i) = yaw + yawd * delta_t + half_delta_t_sqr * nu_yawdd;
    Xsig_pred_(4, i) = yawd + delta_t * nu_yawdd;
  }

  /*
   * Predict Mean and Covariance
   */
  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Normalize(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }

}

//******************************************************************************
void UKF::UpdateLidar(const MeasurementPackage &measurement) {
  MatrixXd H = MatrixXd(2, n_x_);
  H << 1, 0, 0, 0, 0, 0, 1, 0, 0, 0;

  MatrixXd R = MatrixXd(2, 2);
  R << 0.0225, 0, 0, 0.0225;

  VectorXd y = measurement.raw_measurements_ - H * x_;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R;
  MatrixXd Si = S.inverse();
  MatrixXd K = P_ * Ht * Si;

  //new state
  x_ = x_ + (K * y);
  long x_size = x_.size();
  MatrixXd I = MatrixXd::Identity(x_size, x_size);
  P_ = (I - K * H) * P_;
}

//******************************************************************************
void UKF::UpdateRadar(const MeasurementPackage &measurement) {
  /**
   TODO: calculate the radar NIS.
   */
  //transform sigma points into measurement space
  static const int n_z = 3;  //r, phi, and r_dot
  MatrixXd Zsig = MatrixXd(n_z, n_sig_);
  for (int i = 0; i < n_sig_; i++) {
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    Zsig(0, i) = hypot(p_x, p_y);
    Zsig(1, i) = atan2(p_y, p_x);
    Zsig(2, i) = (p_x * cos(yaw) * v + p_y * sin(yaw) * v) / Zsig(0, i);
  }

  /*
   * Predict Measurement
   */
  //calculate mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = Normalize(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  MatrixXd R = MatrixXd(n_z, n_z);
  R << std_radr_ * std_radr_, 0, 0, 0, std_radphi_ * std_radphi_, 0, 0, 0, std_radrd_
      * std_radrd_;
  S = S + R;

  /*
   * Update State
   */
  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = Normalize(z_diff(1));

    // state difference
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Normalize(x_diff(3));

    Tc += weights_(i) * x_diff * z_diff.transpose();
  }

  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();

  //update state mean and covariance matrix
  //residual
  VectorXd z_diff = measurement.raw_measurements_ - z_pred;
  z_diff(1) = Normalize(z_diff(1));
  x_ += K * z_diff;
  P_ += -K * S * K.transpose();

}
