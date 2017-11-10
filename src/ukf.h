#ifndef UKF_H
#define UKF_H

#include "measurement_package.h"
#include "Eigen/Dense"
#include <vector>
#include <string>
#include <fstream>

using Eigen::MatrixXd;
using Eigen::VectorXd;

class UKF {
 public:

  /**
   * Constructor
   */
  UKF();

  /**
   * Destructor
   */
  virtual ~UKF();

  /**
   * Run the Unscented Kalman Filter process and return the state
   */
  VectorXd ProcessMeasurement(const MeasurementPackage &);

 private:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  const double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  const double std_yawdd_;

  ///* Laser measurement noise covariance matrix
  MatrixXd R_lidar_;

  ///* Radar measurement noise covariance matrix
  MatrixXd R_radar_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  const int n_x_;

  ///* Radar dimension (rho, phi, and rho_dot)
  const int n_z_;

  ///* Augmented state dimension
  const int n_aug_;

  ///* Number of sigma points
  const int n_sig_;

  ///* Sigma point spreading parameter
  const double lambda_;

  /**
   * Process the Measurements for Lidar and Radar (splitting into methods to
   * enable change to Template Method in the event another sensor type is added
   * later.
   */
  void ProcessLidarMeasurement(const MeasurementPackage &);
  void ProcessRadarMeasurement(const MeasurementPackage &);

  /**
   * Initialize the state and time
   */
  void Init(float, float, long long);

  /**
   * Prediction methods predict sigma points, the state, and the state covariance
   * matrix.  The functionality is split into methods to make it easier for
   * testing with known data points (i.e. from the lectures).
   */
  void Prediction(double);
  MatrixXd Prediction_GenerateSigmaPoints();
  MatrixXd Prediction_PredictSigmaPoints(double dt, MatrixXd Xsig_aug);
  MatrixXd Prediction_PredictMeanAndCovariance(VectorXd weights,
                                               MatrixXd Xsig_pred);

  /**
   * Updates the state and the state covariance matrix for laser and radar
   * measurements
   */
  void UpdateLidar(const MeasurementPackage &);
  void UpdateRadar(const MeasurementPackage &);

  /**
   * Normalizes the angles if needed
   */
  float Normalize(float angle_radians);

  /**
   * Method to test implementation against the lectures
   */
  void test_ProcessRadarMeasurement();

};

#endif /* UKF_H */
