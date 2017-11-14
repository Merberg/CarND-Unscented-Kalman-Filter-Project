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
  void ProcessMeasurement(const MeasurementPackage);

  ///* state vector: [pos1 pos2 vel_abs yaw_angle yaw_rate] in SI units and rad
  VectorXd x_;

 private:
  ///* initially set to false, set to true in first call of ProcessMeasurement
  bool is_initialized_;

  ///* if this is false, laser measurements will be ignored (except for init)
  bool use_laser_;

  ///* if this is false, radar measurements will be ignored (except for init)
  bool use_radar_;

  ///* state covariance matrix
  MatrixXd P_;

  ///* predicted sigma points matrix
  MatrixXd Xsig_pred_;

  ///* time when the state is true, in us
  long long time_us_;

  ///* Process noise standard deviation longitudinal acceleration in m/s^2
  double std_a_;

  ///* Process noise standard deviation yaw acceleration in rad/s^2
  double std_yawdd_;

  ///* Laser measurement noise covariance matrix
  MatrixXd R_lidar_;

  ///* Radar measurement noise covariance matrix
  MatrixXd R_radar_;

  ///* Weights of sigma points
  VectorXd weights_;

  ///* State dimension
  const int N_x_;

  ///* Radar dimension (rho, phi, and rho_dot)
  const int N_z_;

  ///* Augmented state dimension
  const int N_AUG_;

  ///* Number of sigma points
  const int N_SIG_;

  ///* Sigma point spreading parameter
  const double LAMBDA_;

  ///* Lidar NIS running average
  float NIS_lidar_;

  ///* Radar NIS running average
  float NIS_radar_;

  /**
   * Process the Measurements for Lidar and Radar (splitting into methods to
   * enable change to Template Method in the event another sensor type is added
   * later.
   */
  void ProcessLidarMeasurement(const MeasurementPackage);
  void ProcessRadarMeasurement(const MeasurementPackage);

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
  void Prediction_PredictSigmaPoints(const double, const MatrixXd&);
  void Prediction_PredictMeanAndCovariance();

  /**
   * Updates the state and the state covariance matrix for laser and radar
   * measurements.  The functionality is split into methods to make it easier for
   * testing with known data points.
   */
  void UpdateLidar(const MeasurementPackage);
  void UpdateLidar_CalculateNIS(const VectorXd&, const VectorXd&,
                                const MatrixXd&);

  void UpdateRadar(const MeasurementPackage);
  void UpdateRadar_PredictMeanAndCovariance(MatrixXd&, VectorXd&, MatrixXd&);
  void UpdateRadar_UpdateState(const VectorXd&, const MatrixXd&,
                               const VectorXd&, const MatrixXd&);
  void UpdateRadar_CalculateNIS(const VectorXd&, const VectorXd&,
                                const MatrixXd&);

  /**
   * Normalizes the angles if needed
   */
  float Normalize(float);

  /**
   * Method to test implementation against the lectures
   */
  void test_Prediction(bool);
  void test_UpdateRadar(bool);
  void test_CompareResults(MatrixXd, MatrixXd, std::string);

};

#endif /* UKF_H */
