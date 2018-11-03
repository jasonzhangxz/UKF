#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>

#define EPS 0.00001 //very small positive number

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  is_initialized_ = false;

  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);
  P_ << 1, 0, 0, 0, 0,
        0, 1, 0, 0, 0,
        0, 0, 1, 0, 0,
        0, 0, 0, 1, 0,
        0, 0, 0, 0, 1;

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1.5;

  // Process noise standard deviation yaw acceleration in rad/s^2
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

  /**
  TODO:

  Complete the initialization. See ukf.h for other member properties.

  Hint: one or more values initialized above might be wildly off...
  */
  //state dimension
  n_x_ = x_.size();

  //Augmented state dimension
  n_aug_ = n_x_ + 2;

  //number of sigma points
  n_sig_ = 2*n_aug_+1;

  //set predicted sigma points matrix dimension
  Xsig_pred_ = MatrixXd(n_x_,n_sig_);

  //sigma points spreading parameter
  lambda_ = 3-n_x_;

  //weights of sigma points
  weights_ = VectorXd(n_sig_);

  //measurement noise covariance matrix initialization
  R_radar_ = MatrixXd(3,3);
  R_radar_ << std_radr_*std_radr_,0,0,
              0, std_radphi_*std_radphi_, 0,
              0,0,std_radrd_*std_radrd_;

  R_lidar_ = MatrixXd(2,2);
  R_lidar_ << std_laspx_*std_laspx_,0,
              0,std_laspy_*std_laspy_;

}

UKF::~UKF() {}

/**
 *  Angle normalization to [-Pi, Pi]
 */
void UKF::AngNorm(double *angle){
  while (*angle>M_PI) {*angle -= 2.0*M_PI;}
  while (*angle<-M_PI){*angle += 2.0*M_PI;}
}

/**
 * @param {MeasurementPackage} meas_package The latest measurement data of
 * either radar or laser.
 */
void UKF::ProcessMeasurement(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Make sure you switch between lidar and radar
  measurements.
  */
  //CTRV model, x_ si [px,py,vel,angle,angle_rate]
  if(!is_initialized_){
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      float rho=meas_package.raw_measurements_[0];
      float phi=meas_package.raw_measurements_[1];
      float rho_dot=meas_package.raw_measurements_[2];
      //convert from polar to cartesian
      float px=rho*cos(phi);
      float py=rho*sin(phi);
      float vx=rho_dot*cos(phi);
      float vy=rho_dot*sin(phi);
      float v =sqrt(vx*vx + vy*vy);
      x_ << px,py,v,0,0;
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER){
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1],0,0,0;
      if (fabs(x_(0)) < EPS and fabs(x_(1)) < EPS){
        x_(0) = EPS;
        x_(1) = EPS;
      }
    }

    //initialize weights
    weights_(0) = lambda_/(lambda_+n_aug_);
    for(int i=1;i<weights_.size();i++){
      weights_(i) = 0.5/(lambda_+n_aug_);
    }

    //save initial time stamp
    time_us_ = meas_package.timestamp_;

    is_initialized_ = true;

    return;
  }

  //calculate time difference btw measurements
  double dt = (meas_package.timestamp_ - time_us_);
  dt /= 1000000.0;//from microsseconds to sencods
  time_us_ = meas_package.timestamp_;

  Prediction(dt) ;

  //Update
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR && use_radar_){
    UpdateRadar(meas_package);
  }

  if(meas_package.sensor_type_ == MeasurementPackage::LASER and use_laser_){
    UpdateLidar(meas_package);
  }

}

/**
 * Predicts sigma points, the state, and the state covariance matrix.
 * @param {double} delta_t the change in time (in seconds) between the last
 * measurement and this one.
 */
void UKF::Prediction(double delta_t) {
  /**
  TODO:

  Complete this function! Estimate the object's location. Modify the state
  vector, x_. Predict sigma points, the state, and the state covariance matrix.
  */
  double dt_2 = delta_t * delta_t;

  //Augmented mean vector
  VectorXd x_aug = VectorXd(n_aug_);

  //Augmented state covariance Matrix
  MatrixXd P_aug = MatrixXd(n_aug_,n_aug_);

  //sigma point matrix
  MatrixXd Xsig_aug = MatrixXd(n_aug_,n_sig_);

  //Fill the matrices
  x_aug.fill(0.0);
  x_aug.head(n_x_) = x_;
  P_aug.fill(0);
  P_aug.topLeftCorner(n_x_,n_x_) = P_;
  P_aug(5,5) = std_a_*std_a_;
  P_aug(6,6) = std_yawdd_*std_yawdd_;

  MatrixXd L = P_aug.llt().matrixL();

  //create sigma points
  Xsig_aug.col(0) = x_aug;
  double sqrt_lambda_n_aug = sqrt(lambda_+n_aug_);
  VectorXd sqrt_lambda_n_aug_L;
  for(int i=0;i<n_aug_;i++){
    sqrt_lambda_n_aug_L = sqrt_lambda_n_aug * L.col(i);
    Xsig_aug.col(i+1) = x_aug + sqrt_lambda_n_aug_L;
    Xsig_aug.col(i+1+n_aug_) = x_aug - sqrt_lambda_n_aug_L;
  }

  //Predict sigma points
  double px,py,v,yaw,yawd,nu_a,nu_yawdd,sin_yaw,cos_yaw,arg;
  double px_p,py_p,v_p,yaw_p,yawd_p;
  for(int i=0;i<n_sig_;i++){
    //extract values for better writting
    px = Xsig_aug(0,i);
    py = Xsig_aug(1,i);
    v  = Xsig_aug(2,i);
    yaw= Xsig_aug(3,i);
    yawd=Xsig_aug(4,i);
    nu_a=Xsig_aug(5,i);
    nu_yawdd=Xsig_aug(6,i);

    //precalculate some values
    sin_yaw = sin(yaw);
    cos_yaw = cos(yaw);
    arg = yaw + yawd*delta_t;

    //predicted state values
    if(fabs(yawd)>EPS){
      double v_yawd = v/yawd;
      px_p = px + v_yawd *(sin(arg) - sin_yaw);
      py_p = py + v_yawd *(cos_yaw - cos(arg));
    } else {
      px_p = px + v*delta_t*cos_yaw;
      py_p = py + v*delta_t*sin_yaw;
    }

    v_p = v;
    yaw_p = arg;
    yawd_p= yawd;

    //Add noise
    px_p += 0.5*nu_a*dt_2*cos_yaw;
    py_p += 0.5*nu_a*dt_2*sin_yaw;
    v_p  += nu_a*delta_t;
    yaw_p+= 0.5*nu_yawdd*dt_2;
    yawd_p+=nu_yawdd*delta_t;

    //write predicted sigma point into the right column
    Xsig_pred_(0,i) = px_p;
    Xsig_pred_(1,i) = py_p;
    Xsig_pred_(2,i) = v_p;
    Xsig_pred_(3,i) = yaw_p;
    Xsig_pred_(4,i) = yawd_p;
  }

  //Predicted state mean
  x_ = Xsig_pred_ * weights_;

  //predicted state covariance matrix
  P_.fill(0.0);
  for(int i=0;i<n_sig_;i++){
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    AngNorm(&(x_diff(3)));
    P_ = P_ + weights_(i) * x_diff * x_diff.transpose();
  }

}

/**
 * Updates the state and the state covariance matrix using a laser measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateLidar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use lidar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the lidar NIS.
  */
  int n_z = 2; //Lidar measurement dimension - x,y
  MatrixXd Zsig = Xsig_pred_.block(0,0,n_z,n_sig_); //sigma points in Lidar measurement space

  UpdateUKF(meas_package,Zsig,n_z);
}

/**
 * Updates the state and the state covariance matrix using a radar measurement.
 * @param {MeasurementPackage} meas_package
 */
void UKF::UpdateRadar(MeasurementPackage meas_package) {
  /**
  TODO:

  Complete this function! Use radar data to update the belief about the object's
  position. Modify the state vector, x_, and covariance, P_.

  You'll also need to calculate the radar NIS.
  */
  int n_z = 3; //radar measurement dimension - rho, phi,rho_dot
  MatrixXd Zsig = MatrixXd(n_z,n_sig_);//sigma points in measurement

  double px,py,v,yaw,v1,v2;
  for(int i=0;i<n_sig_;i++){
    px = Xsig_pred_(0,i);
    py = Xsig_pred_(1,i);
    v  = Xsig_pred_(2,i);
    yaw= Xsig_pred_(3,i);
    v1 = cos(yaw)*v;
    v2 = sin(yaw)*v;

    Zsig(0,i) = sqrt(px*px + py*py); //rho
    Zsig(1,i) = atan2(py,px);  //phi
    Zsig(2,i) = (px*v1 + py*v2)/Zsig(0,i); //rho_dot

  }

  UpdateUKF(meas_package,Zsig,n_z);

}

void UKF::UpdateUKF(MeasurementPackage meas_package, MatrixXd Zsig, int n_z){
  //predicted measurement mean
  VectorXd z_pred = VectorXd(n_z);
  z_pred = Zsig * weights_;

  //measurement covariance matrix
  MatrixXd S = MatrixXd(n_z,n_z);
  S.fill(0.0);
  for(int i=0;i<n_sig_;i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    AngNorm(&(z_diff(1)));
    S = S + weights_(i) * z_diff *z_diff.transpose();
  }

  //add noise
  MatrixXd R = MatrixXd(n_z,n_z);
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    R = R_radar_;
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    R = R_lidar_;
  }
  S = S + R;

  //correlation matrix
  MatrixXd Tc = MatrixXd(n_x_,n_z);
  Tc.fill(0.0);
  for(int i=0;i<n_sig_;i++){
    VectorXd z_diff = Zsig.col(i) - z_pred;
    if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
      AngNorm(&(z_diff(1)));
    }
    VectorXd x_diff = Xsig_pred_.col(i) - x_;
    AngNorm(&(x_diff(3)));
    Tc = Tc + weights_(i) * x_diff * z_diff.transpose();
  }

  VectorXd z = meas_package.raw_measurements_;
  MatrixXd K = Tc * S.inverse();
  VectorXd z_diff = z - z_pred;
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    AngNorm(&(z_diff(1)));
  }

  //Update state mean and covariance Matrix
  x_ = x_ + K * z_diff;
  P_ = P_ - K * S * K.transpose();

  //calculate NIS
  if(meas_package.sensor_type_ == MeasurementPackage::RADAR){
    NIS_radar_ = z.transpose() * S.inverse() * z;
  } else if(meas_package.sensor_type_ == MeasurementPackage::LASER){
    NIS_laser_ = z.transpose() * S.inverse() * z;
  }

}
