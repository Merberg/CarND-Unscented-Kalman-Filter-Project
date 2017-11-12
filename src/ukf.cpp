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
      std_a_(1.8),
      std_yawdd_(0.7),
      n_x_(5),
      n_z_(3),
      n_aug_(7),
      n_sig_(2 * n_aug_ + 1),
      lambda_(3 - n_aug_),
      NIS_radar_(7.815)
{
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

  //Uncomment for unit tests
  //test_Prediction(false);
  //test_UpdateRadar(true);
}

UKF::~UKF()
{
}

//******************************************************************************
VectorXd UKF::ProcessMeasurement(const MeasurementPackage &measurement)
{

  //Split depending on the measurement type
  if (measurement.sensor_type_ == MeasurementPackage::LASER) {
    ProcessLidarMeasurement(measurement);
  } else {
    ProcessRadarMeasurement(measurement);
  }

  return x_;
}

//******************************************************************************
void UKF::ProcessLidarMeasurement(const MeasurementPackage &measurement)
{

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
void UKF::ProcessRadarMeasurement(const MeasurementPackage &measurement)
{

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
void UKF::Init(float px, float py, long long ts)
{
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
void UKF::Prediction(double timestamp)
{

  double dt = (timestamp - time_us_) / 1000000.0;  //dt in seconds
  time_us_ = timestamp;

  cout << "dt:" << dt << endl;

  MatrixXd Xsig_aug = Prediction_GenerateSigmaPoints();
  Prediction_PredictSigmaPoints(dt, Xsig_aug);
  Prediction_PredictMeanAndCovariance();

}

//******************************************************************************
MatrixXd UKF::Prediction_GenerateSigmaPoints()
{
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
void UKF::Prediction_PredictSigmaPoints(const double dt,
                                        const MatrixXd& Xsig_aug)
{
  double half_dt_sqr = 0.5 * dt * dt;

  for (int i = 0; i < n_sig_; i++) {
    double px = Xsig_aug(0, i);
    double py = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    if (fabs(yawd) > 0.001) {
      Xsig_pred_(0, i) = px + v / yawd * (sin(yaw + yawd * dt) - sin(yaw));
      Xsig_pred_(1, i) = py + v / yawd * (-cos(yaw + yawd * dt) + cos(yaw));
    } else {
      Xsig_pred_(0, i) = px + v * cos(yaw) * dt;
      Xsig_pred_(1, i) = py + v * sin(yaw) * dt;
    }
    Xsig_pred_(0, i) += half_dt_sqr * cos(yaw) * nu_a;
    Xsig_pred_(1, i) += half_dt_sqr * sin(yaw) * nu_a;
    Xsig_pred_(2, i) = v + dt * nu_a;
    Xsig_pred_(3, i) = yaw + yawd * dt + half_dt_sqr * nu_yawdd;
    Xsig_pred_(4, i) = yawd + dt * nu_yawdd;
  }
}

//******************************************************************************
void UKF::Prediction_PredictMeanAndCovariance()
{
  //predict state mean
  x_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    x_ += weights_(i) * Xsig_pred_.col(i);
  }

  //predict state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    x_diff(3) = Normalize(x_diff(3));
    P_ += weights_(i) * x_diff * x_diff.transpose();
  }
}

//******************************************************************************
void UKF::UpdateLidar(const MeasurementPackage &measurement)
{
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
void UKF::UpdateRadar(const MeasurementPackage &measurement)
{
  MatrixXd Zsig = MatrixXd(n_z_, n_sig_);
  VectorXd z = measurement.raw_measurements_;
  VectorXd z_pred = VectorXd(n_z_);
  MatrixXd S = MatrixXd(n_z_, n_z_);

  UpdateRadar_PredictMeanAndCovariance(Zsig, z_pred, S);
  UpdateRadar_UpdateState(z, Zsig, z_pred, S);
  UpdateRadar_CalculateNIS(z, z_pred, S);
}

//******************************************************************************
void UKF::UpdateRadar_PredictMeanAndCovariance(MatrixXd& Zsig, VectorXd& z_pred,
                                               MatrixXd& S)
{

  //transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++) {
    double px = Xsig_pred_(0, i);
    double py = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);

    Zsig(0, i) = hypot(px, py);
    Zsig(1, i) = atan2(py, px);
    Zsig(2, i) = (px * cos(yaw) * v + py * sin(yaw) * v) / Zsig(0, i);
  }

  //calculate mean predicted measurement
  z_pred.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    z_pred += weights_(i) * Zsig.col(i);
  }

  //calculate measurement covariance matrix S
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++) {
    // state difference
    VectorXd z_diff = Zsig.col(i) - z_pred;
    z_diff(1) = Normalize(z_diff(1));

    S += weights_(i) * z_diff * z_diff.transpose();
  }

  //add measurement noise covariance matrix
  S += R_radar_;
}

//******************************************************************************
void UKF::UpdateRadar_UpdateState(const VectorXd& z, const MatrixXd& Zsig,
                                  const VectorXd& z_pred, const MatrixXd& S)
{
  //calculate cross correlation matrix
  MatrixXd Tc = MatrixXd(n_x_, n_z_);
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
  VectorXd z_diff = z - z_pred;
  z_diff(1) = Normalize(z_diff(1));
  x_ += K * z_diff;
  P_ += -K * S * K.transpose();

}

//******************************************************************************
void UKF::UpdateRadar_CalculateNIS(const VectorXd& z, const VectorXd& z_pred,
                                   const MatrixXd& S)
{
  VectorXd z_diff = z_pred - z;
  z_diff(1) = Normalize(z_diff(1));
  float nis = z_diff.transpose() * S.inverse() * z_diff;
  NIS_radar_ = (NIS_radar_ + nis) / 2.0;
  cout << "Radar NIS:" << nis << " (" << NIS_radar_ << ")" << endl;
}

//******************************************************************************
float UKF::Normalize(float angle_radians)
{
  //float int_part = 0.0;
  //return modff((angle_radians + M_PI)/2.*M_PI, &int_part) - M_PI;
  while (angle_radians > M_PI)
    angle_radians -= 2. * M_PI;
  while (angle_radians < -M_PI)
    angle_radians += 2. * M_PI;
  return angle_radians;
}

//******************************************************************************
void UKF::test_Prediction(bool exitAtEnd)
{
  //***Test 1 Prediction_GenerateSigmaPoints************************************
  x_ << 5.7441, 1.3800, 2.2049, 0.5015, 0.3528;
  std_a_ = 0.2;
  std_yawdd_ = 0.2;
  P_ <<  //
  0.0043, -0.0013, 0.0030, -0.0022, -0.0020,  //
  -0.0013, 0.0077, 0.0011, 0.0071, 0.0060,  //
  0.0030, 0.0011, 0.0054, 0.0007, 0.0008,  //
  -0.0022, 0.0071, 0.0007, 0.0098, 0.0100,  //
  -0.0020, 0.0060, 0.0008, 0.0100, 0.0123;

  MatrixXd actual_Xsig_aug = Prediction_GenerateSigmaPoints();

  MatrixXd expect_Xsig_aug = MatrixXd(n_aug_, n_sig_);
  expect_Xsig_aug <<  //
      5.7441, 5.85768, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.63052, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441, 5.7441,  //
  1.38, 1.34566, 1.52806, 1.38, 1.38, 1.38, 1.38, 1.38, 1.41434, 1.23194, 1.38, 1.38, 1.38, 1.38, 1.38,  //
  2.2049, 2.28414, 2.24557, 2.29582, 2.2049, 2.2049, 2.2049, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.2049, 2.2049,  //
  0.5015, 0.44339, 0.631886, 0.516923, 0.595227, 0.5015, 0.5015, 0.5015, 0.55961, 0.371114, 0.486077, 0.407773, 0.5015, 0.5015, 0.5015,  //
  0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.3528, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.3528,  //
  0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641, 0,  //
  0, 0, 0, 0, 0, 0, 0, 0.34641, 0, 0, 0, 0, 0, 0, -0.34641;

  test_CompareResults(expect_Xsig_aug, actual_Xsig_aug,
                      "Prediction_GenerateSigmaPoints");

  //***Test 2 Prediction_PredictSigmaPoints*************************************
  double test_dt = 0.1;
  MatrixXd expect_Xsig_pred = MatrixXd(n_x_, n_sig_);
  expect_Xsig_pred <<  //
      5.93553, 6.06251, 5.92217, 5.9415, 5.92361, 5.93516, 5.93705, 5.93553, 5.80832, 5.94481, 5.92935, 5.94553, 5.93589, 5.93401, 5.93553,  //
  1.48939, 1.44673, 1.66484, 1.49719, 1.508, 1.49001, 1.49022, 1.48939, 1.5308, 1.31287, 1.48182, 1.46967, 1.48876, 1.48855, 1.48939,  //
  2.2049, 2.28414, 2.24557, 2.29582, 2.2049, 2.2049, 2.23954, 2.2049, 2.12566, 2.16423, 2.11398, 2.2049, 2.2049, 2.17026, 2.2049,  //
  0.53678, 0.473387, 0.678098, 0.554557, 0.643644, 0.543372, 0.53678, 0.538512, 0.600173, 0.395462, 0.519003, 0.429916, 0.530188, 0.53678, 0.535048,  //
  0.3528, 0.299973, 0.462123, 0.376339, 0.48417, 0.418721, 0.3528, 0.387441, 0.405627, 0.243477, 0.329261, 0.22143, 0.286879, 0.3528, 0.318159;

  Prediction_PredictSigmaPoints(test_dt, expect_Xsig_aug);

  test_CompareResults(expect_Xsig_pred, Xsig_pred_,
                      "Prediction_PredictSigmaPoints");

  //***Test 3 Prediction_PredictMeanAndCovariance*******************************
  Xsig_pred_ <<  //
      5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,  //
  1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,  //
  2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,  //
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,  //
  0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

  MatrixXd expect_P = MatrixXd(n_x_, n_x_);
  expect_P <<  //
      0.00543425, -0.0024053, 0.00341576, -0.00348196, -0.00299378,  //
  -0.0024053, 0.010845, 0.0014923, 0.00980182, 0.00791091,  //
  0.00341576, 0.0014923, 0.00580129, 0.000778632, 0.000792973,  //
  -0.00348196, 0.00980182, 0.000778632, 0.0119238, 0.0112491,  //
  -0.00299378, 0.00791091, 0.000792973, 0.0112491, 0.0126972;

  Prediction_PredictMeanAndCovariance();

  test_CompareResults(expect_P, P_, "Prediction_PredictMeanAndCovariance");

  if (exitAtEnd)
    exit(0);
}

//******************************************************************************
void UKF::test_UpdateRadar(bool exitAtEnd)
{
//***Test 1 UpdateRadar_PredictMeanAndCovariance******************************
  R_radar_ <<  //
      0.3 * 0.3, 0, 0,  //
  0, 0.0175 * 0.0175, 0,  //
  0, 0, 0.1 * 0.1;
  Xsig_pred_ <<  //
      5.9374, 6.0640, 5.925, 5.9436, 5.9266, 5.9374, 5.9389, 5.9374, 5.8106, 5.9457, 5.9310, 5.9465, 5.9374, 5.9359, 5.93744,  //
  1.48, 1.4436, 1.660, 1.4934, 1.5036, 1.48, 1.4868, 1.48, 1.5271, 1.3104, 1.4787, 1.4674, 1.48, 1.4851, 1.486,  //
  2.204, 2.2841, 2.2455, 2.2958, 2.204, 2.204, 2.2395, 2.204, 2.1256, 2.1642, 2.1139, 2.204, 2.204, 2.1702, 2.2049,  //
  0.5367, 0.47338, 0.67809, 0.55455, 0.64364, 0.54337, 0.5367, 0.53851, 0.60017, 0.39546, 0.51900, 0.42991, 0.530188, 0.5367, 0.535048,  //
  0.352, 0.29997, 0.46212, 0.37633, 0.4841, 0.41872, 0.352, 0.38744, 0.40562, 0.24347, 0.32926, 0.2214, 0.28687, 0.352, 0.318159;

  MatrixXd Zsig = MatrixXd(n_z_, n_sig_);
  VectorXd actual_z_pred = VectorXd(n_z_);
  VectorXd expect_z_pred = VectorXd(n_z_);
  expect_z_pred << 6.12155, 0.245993, 2.10313;

  MatrixXd actual_S = MatrixXd(n_z_, n_z_);
  MatrixXd expect_S = MatrixXd(n_z_, n_z_);
  expect_S <<  //
      0.0946171, -0.000139448, 0.00407016,  //
  -0.000139448, 0.000617548, -0.000770652,  //
  0.00407016, -0.000770652, 0.0180917,

  UpdateRadar_PredictMeanAndCovariance(Zsig, actual_z_pred, actual_S);

  test_CompareResults(expect_z_pred, actual_z_pred,
                      "UpdateRadar_PredictMeanAndCovariance z_pred");
  test_CompareResults(expect_S, actual_S,
                      "UpdateRadar_PredictMeanAndCovariance S");

//***Test 2 UpdateRadar_UpdateState*******************************************
  x_ << 5.93637, 1.49035, 2.20528, 0.536853, 0.353577;
  VectorXd expect_x = VectorXd(n_x_);
  expect_x << 5.92276, 1.41823, 2.15593, 0.489274, 0.321338;
  P_ <<  //
  0.0054342, -0.002405, 0.0034157, -0.0034819, -0.00299378,  //
  -0.002405, 0.01084, 0.001492, 0.0098018, 0.00791091,  //
  0.0034157, 0.001492, 0.0058012, 0.00077863, 0.000792973,  //
  -0.0034819, 0.0098018, 0.00077863, 0.011923, 0.0112491,  //
  -0.0029937, 0.0079109, 0.00079297, 0.011249, 0.0126972;
  MatrixXd expect_P = MatrixXd(n_x_, n_x_);
  expect_P <<  //
      0.00361579, -0.000357881, 0.00208316, -0.000937196, -0.00071727,  //
  -0.000357881, 0.00539867, 0.00156846, 0.00455342, 0.00358885,  //
  0.00208316, 0.00156846, 0.00410651, 0.00160333, 0.00171811,  //
  -0.000937196, 0.00455342, 0.00160333, 0.00652634, 0.00669436,  //
  -0.00071719, 0.00358884, 0.00171811, 0.00669426, 0.00881797;

  Zsig <<  //
      6.1190, 6.2334, 6.1531, 6.1283, 6.1143, 6.1190, 6.1221, 6.1190, 6.0079, 6.0883, 6.1125, 6.1248, 6.1190, 6.1188, 6.12057,  //
  0.24428, 0.2337, 0.27316, 0.24616, 0.24846, 0.24428, 0.24530, 0.24428, 0.25700, 0.21692, 0.24433, 0.24193, 0.24428, 0.24515, 0.245239,  //
  2.1104, 2.2188, 2.0639, 2.187, 2.0341, 2.1061, 2.1450, 2.1092, 2.0016, 2.129, 2.0346, 2.1651, 2.1145, 2.0786, 2.11295;
  VectorXd z = VectorXd(n_z_);
  z << 5.9214, 0.2187, 2.0062;

  UpdateRadar_UpdateState(z, Zsig, expect_z_pred, expect_S);

  test_CompareResults(expect_x, x_, "UpdateRadar_UpdateState x");
  test_CompareResults(expect_P, P_, "UpdateRadar_UpdateState P");

  if (exitAtEnd)
    exit(0);
}

//******************************************************************************
void UKF::test_CompareResults(MatrixXd expected, MatrixXd actual,
                              std::string testString)
{
  const float EPSILON = 1e-5;
  MatrixXd diff = (expected - actual).cwiseAbs();

  if ((expected - actual).norm() < EPSILON || diff.maxCoeff() < EPSILON) {
    cout << "Test " << testString << " PASSED" << endl;
  } else {
    cout << "Test " << testString << " FAILED" << endl;
    cout << "Actual:" << actual << endl;
    cout << endl;
  }
}

