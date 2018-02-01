#include "ukf.h"
#include "Eigen/Dense"
#include <iostream>
#include "tools.h"
using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/**
 * Initializes Unscented Kalman filter
 * This is scaffolding, do not modify
 */
UKF::UKF() {
  // if this is false, laser measurements will be ignored (except during init)
  use_laser_ = true;

  // if this is false, radar measurements will be ignored (except during init)
  use_radar_ = true;

  // state dim
  n_x_= 5;
  // initial state vector
  x_ = VectorXd(5);

  // initial covariance matrix
  P_ = MatrixXd(5, 5);

  // Process noise standard deviation longitudinal acceleration in m/s^2
  std_a_ = 1;

  // Process noise standard deviation yaw acceleration in rad/s^2
  std_yawdd_ = 0.07;
  
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
  

  // Set augmentated state dimension
  n_aug_ = 7;
  
  // Set sigma point spraeding parameter
  lambda_ = 3 - n_aug_;
  cout << "T:contruction" << endl;
  weights_ = VectorXd(2 * n_aug_ +1);
  weights_(0) = lambda_ / (lambda_ + n_aug_);
  for (int i=1; i< 2 * n_aug_ + 1; i++){
      weights_(i) = 1/(2 * (lambda_ + n_aug_));
  }
    
  
  is_initialized_ = false;

}

UKF::~UKF() {}

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
  /*********Initialization**********/    
  if (!is_initialized_) {   
    counter_debug_=0;
    cout << "Initialize" << endl;  
	  x_ = VectorXd(n_x_);
    P_ = MatrixXd(n_x_, n_x_);
    P_ = MatrixXd::Identity(n_x_,n_x_)*1000;  	    
    if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {
      // Convert RADA measurement to x_
      // code borrowed from EKF project      
      float ro =  meas_package.raw_measurements_[0];
      float theta =  meas_package.raw_measurements_[1];
      float ro_dot =  meas_package.raw_measurements_[2];		
      x_ << ro * cos(theta), ro * sin(theta), 0, 0, 0;
      MatrixXd Hj_pxpy_inv(2,3);
      Hj_pxpy_inv << cos(theta), -ro * sin(theta), 0,
        sin(theta), ro * cos(theta),0;      
      MatrixXd R_radar(3, 3);
      R_radar << std_radr_ * std_radr_, 0, 0,
        0, std_radphi_ * std_radphi_, 0,
        0, 0, std_radrd_ * std_radrd_;	  
      MatrixXd R_pxy_radar = Hj_pxpy_inv * R_radar * Hj_pxpy_inv.transpose();      
      P_(0, 0) = R_pxy_radar(0, 0);
      P_(0, 1) = R_pxy_radar(0, 1);
      P_(1, 0) = R_pxy_radar(1, 0);
      P_(1, 1) = R_pxy_radar(1, 1);	
      time_us_ = meas_package.timestamp_;
      is_initialized_ = true;           
    } else if (meas_package.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */	      
      x_ << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], 0, 0, 0;	  	                   
      P_(0, 0) = std_laspx_ * std_laspx_;	  	 
      P_(1, 1) = std_laspy_ * std_laspy_;
      time_us_ = meas_package.timestamp_;	       
      is_initialized_ = true;      
    }

      // done initializing, no need to predict or update    
    return;
  }
  counter_debug_++;
  cout << "counter " << counter_debug_ <<endl;
  /*********Prediction**********/      
  Prediction((meas_package.timestamp_ - time_us_) / 1000000.0);  
  //cout << "Predict : x_" << x_ << endl;
  //cout << "Predict: P_" << P_ << endl;
  /*********Update**********/
  
  if (meas_package.sensor_type_ == MeasurementPackage::RADAR) {     
    cout << "RADAR" << endl;
    UpdateRadar(meas_package);
    time_us_ = meas_package.timestamp_;    
  } else {    
    cout << "LIDAR" << endl;
    UpdateLidar(meas_package);
    time_us_ = meas_package.timestamp_;
  }  
  
  cout << "Update: x_" << x_ << endl;
  cout << "Update: P_" << P_ << endl;
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
  VectorXd x_aug = VectorXd(n_aug_);
  x_aug << x_, 0, 0;
  //create augmented covariance matrix
  MatrixXd Q = MatrixXd(2, 2);
  Q << std_a_ * std_a_, 0,
    0, std_yawdd_ * std_yawdd_;
  MatrixXd P_aug = MatrixXd(n_aug_, n_aug_);
  P_aug << P_, MatrixXd::Zero(n_x_, 2),
        MatrixXd::Zero(2, n_x_), Q;
  //create square root matrix
  MatrixXd A = P_aug.llt().matrixL();
  //create augmented sigma points
  MatrixXd Xsig_aug = MatrixXd(n_aug_, 2 * n_aug_ + 1);
  Xsig_aug << x_aug,
    x_aug.replicate(1, n_aug_) + sqrt(lambda_ + n_aug_) * A,
    x_aug.replicate(1, n_aug_) - sqrt(lambda_ + n_aug_) *  A;  	 
	
  //predict sigma points
  Xsig_pred_ = MatrixXd(n_x_, 2 * n_aug_ + 1);  
  Xsig_pred_.fill(0);
  
  for (int i = 0; i < 2 * n_aug_ +1; i++)
  {
    float v = Xsig_aug(2, i);
    float psi = Xsig_aug(3, i);
    float psi_dot = Xsig_aug(4, i);
    float va = Xsig_aug(5, i);
    float vpsi = Xsig_aug(6, i);
    if (fabs(psi_dot)<=0.001) {
      Xsig_pred_(0, i) = Xsig_aug(0, i) + v * cos(psi) * delta_t +
        0.5*delta_t*delta_t * cos(psi) * va;
      Xsig_pred_(1, i) = Xsig_aug(1, i) + v * sin(psi) * delta_t +
        0.5*delta_t*delta_t * sin(psi) * va;
    } else {
      Xsig_pred_(0, i) = Xsig_aug(0, i) + 
        v / psi_dot * (sin(psi + psi_dot *delta_t) - sin(psi)) + 
        0.5*delta_t*delta_t * cos(psi) * va;
      Xsig_pred_(1, i) = Xsig_aug(1, i) + 
        v / psi_dot * (-cos(psi + psi_dot *delta_t) + cos(psi)) +
        0.5*delta_t*delta_t * sin(psi) * va;
    }
    Xsig_pred_(2, i) = Xsig_aug(2, i) + delta_t * va;            
    Xsig_pred_(3, i) = Xsig_aug(3, i) + psi_dot * delta_t + 0.5 * delta_t * delta_t * vpsi; 
    Xsig_pred_(4, i) = Xsig_aug(4, i) + delta_t * vpsi;	
  }
  
  //Xsig_pred_(3,0) = Tools::CalculateAngle(Xsig_pred_(3,0), 0);
  //for (int i = 1; i < 2 * n_aug_ +1; i++) {
  // Xsig_pred_(3, i) = Tools::CalculateAngle(Xsig_pred_(3, i), Xsig_pred_(3,0));    
  //}     
    
  //predict state mean
  x_ = Xsig_pred_ * weights_;  
  //predict state covariance matrix
  MatrixXd Xsig_pred_residue = Xsig_pred_ - x_.replicate(1, 2*n_aug_+1);  
  for (int i = 0; i < 2 * n_aug_ +1; i++) {
    Xsig_pred_residue(3, i) = Tools::FoldAngle(Xsig_pred_residue(3, i), 0);    
  }    
  x_(3) = Tools::FoldAngle(x_(3), 0);  
  P_ = Xsig_pred_residue * weights_.asDiagonal() * Xsig_pred_residue.transpose();          
  return;
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
	MatrixXd H = MatrixXd(2,5);
  H << 1, 0, 0, 0, 0,
	  0, 1, 0, 0, 0;  
	VectorXd z = VectorXd(2);
	z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1];		    
  MatrixXd R(2,2);
  R << std_laspx_ * std_laspx_, 0,
    0, std_laspy_ * std_laspy_;
  VectorXd y = z - H * x_;
  MatrixXd S = H * P_ * H.transpose() + R;  
  MatrixXd K = P_ * H.transpose() * S.inverse();
  x_ = x_ + K * y;    
  x_(3) = Tools::FoldAngle(x_(3), 0);    
  
  P_ = P_ - K * H * P_;
  float NIS = y.transpose() * S.inverse() * y;
  
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
  // Get measurement  
  VectorXd z = VectorXd(3);
	z << meas_package.raw_measurements_[0], meas_package.raw_measurements_[1], meas_package.raw_measurements_[2];      
  //transform sigma points into measurement space    
  MatrixXd Zsig = MatrixXd(3, 2 * n_aug_ + 1);
  Zsig.fill(0);  
  for (int i=0; i< 2 * n_aug_ + 1; i++){
    float px = Xsig_pred_(0, i);
    float py = Xsig_pred_(1, i);
    float v = Xsig_pred_(2, i);
    float psi = Xsig_pred_(3, i);
    float psi_dot = Xsig_pred_(4, i);
    float pxy = sqrt(px*px + py*py);
    Zsig(1, i) = atan2(py, px);
    Zsig(0, i) = pxy;
    if (fabs(pxy)<=0.001){
      cout << "*****************" <<endl;
      Zsig(2, i) = v;//0;    
    } else {      
      Zsig(2, i) = (px * cos(psi) * v + py * sin(psi) *v) / pxy;
    }
  }    
  //calculate mean predicted measurement
  //cout << "z " << z << endl;
  //cout << "Zsig " << Zsig << endl;
  VectorXd z_pred = Zsig * weights_;  
  //calculate innovation covariance matrix S
  //cout << "z_pred " << z_pred << endl;
  MatrixXd Z_residue = Zsig - z_pred.replicate(1, 2 * n_aug_ + 1);

  for (int i = 0; i < 2 * n_aug_ +1; i++) {  
    Z_residue(1, i) = Tools::FoldAngle(Z_residue(1, i), 0);     
  }   
  
  MatrixXd R(3,3);
  R << std_radr_*std_radr_, 0, 0,
    0, std_radphi_*std_radphi_, 0,
    0, 0, std_radrd_*std_radrd_;  
  //cout << "Z_residue " << Z_residue << endl;
  MatrixXd S = Z_residue * weights_.asDiagonal() * Z_residue.transpose() + R;    
  
  //calculate cross correlation matrix
  MatrixXd X_residue = Xsig_pred_ - x_.replicate(1, 2 * n_aug_ + 1);  
  for (int i = 0; i < 2 * n_aug_ +1; i++) {    
    X_residue(3, i) = Tools::FoldAngle(X_residue(3, i), 0);  
  }        
  //cout << "Xsig_pred_ " << Xsig_pred_ << endl;
  //cout << "X_residue " << X_residue << endl;
  MatrixXd Tc = X_residue  *  weights_.asDiagonal() * Z_residue.transpose();
  //cout << "Tc " << Tc << endl;
  //calculate Kalman gain K;
  MatrixXd K = Tc * S.inverse();  
  cout << "K " << K << endl;
  //update state mean and covariance matrix
  VectorXd y = z - z_pred;
  y(1) = Tools::FoldAngle(y(1), 0);    
  x_ = x_ + K * y;
  x_(3) = Tools::FoldAngle(x_(3), 0);  
  P_ = P_ - K * S * K.transpose();    
  float NIS = y.transpose() * S.inverse() * y;  
  cout << "NIS: " << NIS << endl;
}
