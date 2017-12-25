#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

static const string NIS_dir = "../NIS/";
static const string radar_NIS_path = NIS_dir + "radar_NIS.txt";
static const string laser_NIS_path = NIS_dir + "laser_NIS.txt";

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF()
{
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 2;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.3;

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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */

  is_initialized_ = false;

  // Open files to store the NIS
  radar_NIS_file_.open(radar_NIS_path, ios::out);
  if (radar_NIS_file_.fail())
  {
    cerr << "Error opening" << radar_NIS_path << endl;
    exit(1);
  }

  laser_NIS_file_.open(laser_NIS_path, ios::out);
  if (laser_NIS_file_.fail())
  {
    cerr << "Error opening" << laser_NIS_path << endl;
    exit(1);
  }

  // State dimension
  n_x_ = x_.size();

  // Augmented state dimension
  n_aug_ = n_x_ + 2; // We will create 2 * n_aug_ + 1 sigma points

  // Number of sigma points
  n_sig_ = 2 * n_aug_ + 1;

  // Sigma point spreading parameter
  lambda_ = 3 - n_x_;

  // predicted sigma points matrix
  Xsig_pred_ = MatrixXd(n_x_, n_sig_);

  //create vector for weights
  weights_ = VectorXd(n_sig_);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i = 1; i < 2 * n_aug_ + 1; i++)
  {
    weights_(i) = 0.5 / (n_aug_ + lambda_);
  }

  // init covariance matrix
  P_ << 1.0, 0, 0, 0, 0,
      0, 1.0, 0, 0, 0,
      0, 0, 1.0, 0, 0,
      0, 0, 0, 1.0, 0,
      0, 0, 0, 0, 1.0;

  // the current NIS for radar
  NIS_radar_ = 0.0;

  // the current NIS for laser
  NIS_laser_ = 0.0;

  n_radar_ = 3;

  n_laser_ = 2;

  // Measurement noise covariance matrices initialization
  R_radar_ = MatrixXd(n_radar_, n_radar_);
  R_radar_ << std_radr_ * std_radr_, 0, 0,
      0, std_radphi_ * std_radphi_, 0,
      0, 0, std_radrd_ * std_radrd_;
  R_laser_ = MatrixXd(n_laser_, n_laser_);
  R_laser_ << std_laspx_ * std_laspx_, 0,
      0, std_laspy_ * std_laspy_;
}

UKF::~UKF()
{
  radar_NIS_file_.close();
  laser_NIS_file_.close();
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Make sure you switch between laser and radar
  measurements.
  */

  if (!is_initialized_)
  {
    // Initialize the state ekf_.x_ with the first measurement.
    // Create the covariance matrix.

    cout << "Initializing unscented Kalman filter" << endl;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      // Convert radar from polar to cartesian coordinates and initialize state.
      // Convert radar from polar to cartesian coordinates and initialize state.
      float rho = meas_package.raw_measurements_[0];     // range
      float phi = meas_package.raw_measurements_[1];     // bearing
      float rho_dot = meas_package.raw_measurements_[2]; // velocity of rho

      // Coordinates convertion from polar to cartesian
      float px = rho * cos(phi);
      float py = rho * sin(phi);
      float vx = rho_dot * cos(phi);
      float vy = rho_dot * sin(phi);
      float v = sqrt(vx * vx + vy * vy);

      x_ << px, py, v, 0, 0;
    }
    else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
    {
      // Initialize state.
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0., 0., 0.;
    }

    time_us_ = meas_package.timestamp_;

    // done initializing, no need to predict or update
    is_initialized_ = true;
    return;
  }

  // Compute the time from the previous measurement in seconds.
  double dt = (meas_package.timestamp_ - time_us_) / 1000000.0;
  time_us_ = meas_package.timestamp_;

  // First predict then update
  Prediction(dt);

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    // Radar update
    UpdateRadar(meas_package);
  }
  else
  {
    // Laser update
    UpdateLidar(meas_package);
  }
}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t)
{
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */

  double delta_t2 = delta_t * delta_t;
  // Augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);
  // Augmented state covarience matrix
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  // Sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_, n_sig_);

  // Fill the matrices
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_, n_x_) = P_;
  P_aug(5, 5) = std_a_ * std_a_;
  P_aug(6, 6) = std_yawdd_ * std_yawdd_;

  //create square root matrix
  MatrixXd L = P_aug.llt().matrixL();

  //create augmented sigma points
  Xsig_aug.col(0) = x_aug;
  double sqrt_lambda_n_aug = sqrt(lambda_ + n_aug_); // Save some computations

  VectorXd sqrt_lambda_n_aug_L;
  for (int i = 0; i < n_aug_; i++)
  {
    sqrt_lambda_n_aug_L = sqrt_lambda_n_aug * L.col(i);
    Xsig_aug.col(i + 1) = x_aug + sqrt_lambda_n_aug_L;
    Xsig_aug.col(i + 1 + n_aug_) = x_aug - sqrt_lambda_n_aug_L;
  }

  // Run each augmented sigma points through the process model with noise
  // (prediction step)

  for (int i = 0; i < n_sig_; i++)
  {
    // Extract values for better readability
    double p_x = Xsig_aug(0, i);
    double p_y = Xsig_aug(1, i);
    double v = Xsig_aug(2, i);
    double yaw = Xsig_aug(3, i);
    double yawd = Xsig_aug(4, i);
    double nu_a = Xsig_aug(5, i);
    double nu_yawdd = Xsig_aug(6, i);

    // Precalculate sin and cos for optimization
    double sin_yaw = sin(yaw);
    double cos_yaw = cos(yaw);
    double arg = yaw + yawd * delta_t;

    // Predicted state values
    double px_p, py_p;
    // Avoid division by zero
    if (yawd > 0.0)
    {
      double v_yawd = v / yawd;
      px_p = p_x + v_yawd * (sin(arg) - sin_yaw);
      py_p = p_y + v_yawd * (cos_yaw - cos(arg));
    }
    else
    {
      double v_delta_t = v * delta_t;
      px_p = p_x + v_delta_t * cos_yaw;
      py_p = p_y + v_delta_t * sin_yaw;
    }
    double v_p = v;
    double yaw_p = arg;
    double yawd_p = yawd;

    // Add noise
    px_p += 0.5 * nu_a * delta_t2 * cos_yaw;
    py_p += 0.5 * nu_a * delta_t2 * sin_yaw;
    v_p += nu_a * delta_t;
    yaw_p += 0.5 * nu_yawdd * delta_t2;
    yawd_p += nu_yawdd * delta_t;

    // Write predicted sigma point into right column
    Xsig_pred_(0, i) = px_p;
    Xsig_pred_(1, i) = py_p;
    Xsig_pred_(2, i) = v_p;
    Xsig_pred_(3, i) = yaw_p;
    Xsig_pred_(4, i) = yawd_p;
  }

  // Predicted state mean
  x_ = Xsig_pred_ * weights_; // vectorised sum
  // Predicted state covariance matrix
  P_.fill(0.0);
  for (int i = 0; i < n_sig_; i++)
  { //iterate over sigma points
    // State difference
    VectorXd deltax = Xsig_pred_.col(i) - x_;

    //Normalize the angle
    while (deltax(3) > M_PI)
    {
      deltax(3) -= 2. * M_PI;
    }
    while (deltax(3) < -M_PI)
    {
      deltax(3) += 2. * M_PI;
    }

    // outer product
    // https://eigen.tuxfamily.org/dox/group__QuickRefPage.html#matrixonly
    P_ = P_ + weights_(i) * deltax * deltax.transpose();
  }
}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */

  // Create matrix for sigma points in measurement space
  // Transform sigma points into measurement space
  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_laser_, n_sig_);

  // Transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++)
  {
    Zsig(0, i) = Xsig_pred_(0, i);
    Zsig(1, i) = Xsig_pred_(1, i);
  }
  NIS_laser_ = UpdateCommon(meas_package, Zsig, n_laser_);

  laser_NIS_file_ << NIS_laser_ << endl;
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package)
{
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */

  // Create matrix for sigma points in measurement space
  MatrixXd Zsig = MatrixXd(n_radar_, n_sig_);
  // Transform sigma points into measurement space
  for (int i = 0; i < n_sig_; i++)
  {
    // extract values for better readibility
    double p_x = Xsig_pred_(0, i);
    double p_y = Xsig_pred_(1, i);
    double v = Xsig_pred_(2, i);
    double yaw = Xsig_pred_(3, i);
    double v1 = cos(yaw) * v;
    double v2 = sin(yaw) * v;
    // Measurement model
    Zsig(0, i) = sqrt(p_x * p_x + p_y * p_y);        //r
    Zsig(1, i) = atan2(p_y, p_x);                    //phi
    Zsig(2, i) = (p_x * v1 + p_y * v2) / Zsig(0, i); //r_dot
  }

  NIS_radar_ = UpdateCommon(meas_package, Zsig, n_radar_);
  radar_NIS_file_ << NIS_radar_ << endl;
}

// Common update function
double UKF::UpdateCommon(MeasurementPackage meas_package, MatrixXd Zsig, int n_z)
{
  // Mean predicted measurement
  VectorXd z_pred = VectorXd(n_z);

  //calculate mean predicted measurement
  z_pred.fill(0.);
  for (int i = 0; i < n_sig_; i++)
  {
    z_pred = z_pred + weights_(i) * Zsig.col(i);
  }

  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  {
    while (z_pred(1) > M_PI)
    {
      z_pred(1) -= 2. * M_PI;
    }
    while (z_pred(1) < -M_PI)
    {
      z_pred(1) += 2. * M_PI;
    }
  }

  //measurement covariance matrix S
  MatrixXd S = MatrixXd(n_z, n_z);
  S.fill(0.0);
  for (int i = 0; i < n_sig_; i++)
  {
    // Residual
    VectorXd z_diff = Zsig.col(i) - z_pred;

    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    {
      //Normalize the angle
      while (z_diff(1) > M_PI)
      {
        z_diff(1) -= 2. * M_PI;
      }
      while (z_diff(1) < -M_PI)
      {
        z_diff(1) += 2. * M_PI;
      }
    }
    S = S + weights_(i) * z_diff * z_diff.transpose();
  }
  // Add measurement noise covariance matrix
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  { // Radar
    S = S + R_radar_;
  }
  else if (meas_package.sensor_type_ == MeasurementPackage::LASER)
  { // Lidar
    S = S + R_laser_;
  }

  // Create matrix for cross correlation Tc
  MatrixXd Tc = MatrixXd(n_x_, n_z);
  // Calculate cross correlation matrix
  Tc.fill(0.0);
  for (int i = 0; i < n_sig_; i++)
  {
    //residual
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
    { // Radar
      //Normalize the angle
      while (z_diff(1) > M_PI)
      {
        z_diff(1) -= 2. * M_PI;
      }
      while (z_diff(1) < -M_PI)
      {
        z_diff(1) += 2. * M_PI;
      }
    }
    // State difference
    VectorXd deltax = Xsig_pred_.col(i) - x_;

    //Normalize the angle
    while (deltax(3) > M_PI)
    {
      deltax(3) -= 2. * M_PI;
    }
    while (deltax(3) < -M_PI)
    {
      deltax(3) += 2. * M_PI;
    }

    Tc = Tc + weights_(i) * deltax * z_diff.transpose();
  }
  // Measurements
  VectorXd z = meas_package.raw_measurements_;
  //Kalman gain K;
  MatrixXd K = Tc * S.inverse();
  // Residual
  VectorXd z_diff = z - z_pred;
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR)
  { // Radar
    //Normalize the angle
    while (z_diff(1) > M_PI)
    {
      z_diff(1) -= 2. * M_PI;
    }
    while (z_diff(1) < -M_PI)
    {
      z_diff(1) += 2. * M_PI;
    }
  }
  // Update state mean and covariance matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();
  // Calculate NIS

  double NIS = z_diff.transpose() * S.inverse() * z_diff;

  return NIS;
}