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

  ///* Laser measurement noise standard deviation position1 in m
  const double std_laspx_;

  ///* Laser measurement noise standard deviation position2 in m
  const double std_laspy_;

  ///* Radar measurement noise standard deviation radius in m
  const double std_radr_;

  ///* Radar measurement noise standard deviation angle in rad
  const double std_radphi_;

  ///* Radar measurement noise standard deviation radius change in m/s
  const double std_radrd_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  const int n_x_;

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
   * Prediction Predicts sigma points, the state, and the state covariance
   * matrix
   */
  void Prediction(double);

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


};

#endif /* UKF_H */
