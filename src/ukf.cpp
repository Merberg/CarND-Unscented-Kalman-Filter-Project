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
      std_a_(2.5),  //TODO adjust
      std_yawdd_(1),  //TODO adjust
      n_x_(5),
      n_z_(3),
      n_aug_(7),
      n_sig_(2 * n_aug_ + 1),
      lambda_(3 - n_aug_) {
  static const int NOISE_PX = 0.15;
  static const int NOISE_PY = 0.15;
  static const int NOISE_RHO = 0.3;
  static const int NOISE_PHI = 0.03;
  static const int NOISE_RD = 0.3;

  // initial state vector
  x_ = VectorXd(n_x_);

  // initial covariance matrix
  P_ = MatrixXd::Identity(n_x_, n_x_);
  P_ *= 0.15;

  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  weights_ = VectorXd(n_sig_);
  weights_.fill(1 / (2 * (lambda_ + n_aug_)));
  weights_(0) = lambda_ / (lambda_ + n_aug_);

  R_lidar_ = MatrixXd(2, 2);
  R_lidar_ << NOISE_PX, 0,  //
  0, NOISE_PY;

  R_radar_ = MatrixXd(n_z_, n_z_);
  R_radar_ << NOISE_RHO * NOISE_RHO, 0, 0,  //
  0, NOISE_PHI * NOISE_PHI, 0,  //
  0, 0, NOISE_RD * NOISE_RD;

  test_ProcessRadarMeasurement();
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

  return x_;
}

//******************************************************************************
void UKF::ProcessLidarMeasurement(const MeasurementPackage &measurement) {

  if (!use_laser_) {
    cout << "Lidar NOT IN USE" << endl;
    return;
  }
  cout << "Process Lidar" << endl;

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

  //Return the updated state
  cout << "x: " << x_ << endl;
  cout << "P:" << P_ << endl;
  cout << endl;

}

//******************************************************************************
void UKF::ProcessRadarMeasurement(const MeasurementPackage &measurement) {

  if (!use_radar_) {
    cout << "Radar NOT IN USE" << endl;
    return;
  }
  cout << "Process Radar" << endl;

  /*****************************************************************************
   *  Initialization
   ****************************************************************************/
  if (!is_initialized_) {
    float rho = measurement.raw_measurements_(0);
    float phi = measurement.raw_measurements_(1);
    Init(rho * cos(phi), rho * sin(phi), measurement.timestamp_);

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

  //Return the updated state
  cout << "x: " << x_ << endl;
  cout << "P:" << P_ << endl;
  cout << endl;
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

  double dt = (timestamp - time_us_) / 1000000.0;  //dt in seconds
  time_us_ = timestamp;

  cout << "dt:" << dt << endl;

  MatrixXd Xsig_aug = Prediction_GenerateSigmaPoints();
  Xsig_pred_ = Prediction_PredictSigmaPoints(dt, Xsig_aug);
  P_ = Prediction_PredictMeanAndCovariance(weights_, Xsig_pred_);

}

//******************************************************************************
MatrixXd UKF::Prediction_GenerateSigmaPoints() {
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
  return Xsig_aug;
}

//******************************************************************************
MatrixXd UKF::Prediction_PredictSigmaPoints(double dt, MatrixXd Xsig_aug) {
  double half_dt_sqr = 0.5 * dt * dt;
  MatrixXd Xsig_pred = MatrixXd(n_x_, n_sig_);

  for (int i = 0; i < n_sig_; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    if (fabs(yawd) > 0.001) {
      Xsig_pred(0, i) = px + v / yawd * (sin(yaw + yawd * dt) - sin(yaw));
      Xsig_pred(1, i) = py + v / yawd * (-cos(yaw + yawd * dt) + cos(yaw));
    } else {
      Xsig_pred(0, i) = px + v * cos(yaw) * dt;
      Xsig_pred(1, i) = py + v * sin(yaw) * dt;
    }
    Xsig_pred(0, i) += half_dt_sqr * cos(yaw) * nu_a;
    Xsig_pred(1, i) += half_dt_sqr * sin(yaw) * nu_a;
    Xsig_pred(2, i) = v + dt * nu_a;
    Xsig_pred(3, i) = yaw + yawd * dt + half_dt_sqr * nu_yawdd;
    Xsig_pred(4, i) = yawd + dt * nu_yawdd;
  }
  return Xsig_pred;
}

//******************************************************************************
MatrixXd UKF::Prediction_PredictMeanAndCovariance(VectorXd weights,
                                                  MatrixXd Xsig_pred) {
  //predict state mean
  VectorXd x_pred = VectorXd(n_x_);
  x_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_pred += weights(i) * Xsig_pred.col(i);
  }

  //predict state covariance matrix
  MatrixXd P_pred = MatrixXd(n_x_, n_x_);
  P_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd x_diff = Xsig_pred.col(i) - x_pred;
    x_diff(3) = Normalize(x_diff(3));
    P_pred += weights(i) * x_diff * x_diff.transpose();
  }
  return P_pred;
}

//******************************************************************************
void UKF::UpdateLidar(const MeasurementPackage &measurement) {
  MatrixXd H = MatrixXd(2, n_x_);
  H << 1, 0, 0, 0, 0,  //
  0, 1, 0, 0, 0;

  VectorXd y = measurement.raw_measurements_ - H * x_;
  MatrixXd Ht = H.transpose();
  MatrixXd S = H * P_ * Ht + R_lidar_;
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
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    Zsig(0, i) = hypot(px, py);
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / Zsig(0, i);
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
  S += R_radar_;

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

//******************************************************************************
float UKF::Normalize(float angle_radians) {
  //float int_part = 0.0;
  //return modff((angle_radians + M_PI)/2.*M_PI, &int_part) - M_PI;
  while (angle_radians > M_PI)
    angle_radians -= 2. * M_PI;
  while (angle_radians < -M_PI)
    angle_radians += 2. * M_PI;
  return angle_radians;
}


//******************************************************************************
void UKF::test_ProcessRadarMeasurement()
{

}
