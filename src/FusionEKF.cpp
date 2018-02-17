#include "FusionEKF.h"
#include "tools.h"
#include "Eigen/Dense"
#include <iostream>

using namespace std;
using Eigen::MatrixXd;
using Eigen::VectorXd;
using std::vector;

/*
 * Constructor.
 */
FusionEKF::FusionEKF() {
  is_initialized_ = false;

  previous_timestamp_ = 0;

  // initializing matrices those who will be different for radar and laser
  R_laser_ = MatrixXd(2, 2);
  R_radar_ = MatrixXd(3, 3);
  H_laser_ = MatrixXd(2, 4);
  Hj_ = MatrixXd(3, 4);

  //measurement covariance matrix - laser
  R_laser_ << 0.0225, 0,
        0, 0.0225;

  //measurement covariance matrix - radar
  R_radar_ << 0.09, 0, 0,
        0, 0.0009, 0,
        0, 0, 0.09;


  // Measurement matrix - Laser
  H_laser_ << 1,0,0,0,
              0,1,0,0;

  // Hj cannot be initialized as it depends upon px, py and needs to be updated every steps

  // Initialize ekf_ parameters which are not going to change
  //create a 4D state vector, we don't know yet the values of the x state
	ekf_.x_ = VectorXd(4);

	//state covariance matrix P
	ekf_.P_ = MatrixXd(4, 4);

	//the initial transition matrix F_
	ekf_.F_ = MatrixXd(4, 4);

}

/**
* Destructor.
*/
FusionEKF::~FusionEKF() {}

void FusionEKF::ProcessMeasurement(const MeasurementPackage &measurement_pack) {


  /*****************************************************************************
   *  Initialization
   ****************************************************************************/

  //acceleration noise components
  float noise_ax;
  float noise_ay;

  noise_ax = 9;
  noise_ay = 9;

  if (!is_initialized_) {

    // first measurement
    if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
      /**
      Convert radar from polar to cartesian coordinates and initialize state.
      */
      float ro;
      float theta;
      float ro_dot;
      ro = measurement_pack.raw_measurements_[0];
      theta = measurement_pack.raw_measurements_[1];
      ro_dot = measurement_pack.raw_measurements_[2];

      //ro_dot =

      // why not velocity
      float px;
      float py;

      px = ro * cos(theta);
      py = ro * sin(theta);

      // updating the x
      ekf_.x_ << px, py, 5, 0;
      previous_timestamp_ = measurement_pack.timestamp_;

    }
    else if (measurement_pack.sensor_type_ == MeasurementPackage::LASER) {
      /**
      Initialize state.
      */
      //set the state with the initial location and zero velocity
      float px;
      float py;
      px = measurement_pack.raw_measurements_[0];
      py = measurement_pack.raw_measurements_[1];
		  ekf_.x_ << px, py, 5, 0;
      previous_timestamp_ = measurement_pack.timestamp_;

    }

    // done initializing, no need to predict or update
    is_initialized_ = true;
    // covariance matrix remain the same intitlize with high variance value in velocity
    //ekf_.P_.setRandom();
    ekf_.P_ << 1, 0, 0, 0,
  			  0, 1, 0, 0,
  			  0, 0, 10, 0,
  			  0, 0, 0, 10;
   // State Transition function initialization, only change where dt is required
    ekf_.F_ << 1, 0, 1, 0,
  			  0, 1, 0, 1,
  			  0, 0, 1, 0,
  			  0, 0, 0, 1;
    //cout << "x_ = " << ekf_.x_ << endl;
    //cout << "P_ = " << ekf_.P_ << endl;
    return;
  }

  /*****************************************************************************
   *  Prediction
   ****************************************************************************/

  /**
     * Update the state transition matrix F according to the new elapsed time.
      - Time is measured in seconds.
     * Update the process noise covariance matrix.
     * Use noise_ax = 9 and noise_ay = 9 for your Q matrix.
   */

   //compute the time elapsed between the current and previous measurements
 	float dt = (measurement_pack.timestamp_ - previous_timestamp_) / 1000000.0;	//dt - expressed in seconds
 	previous_timestamp_ = measurement_pack.timestamp_;

 	float dt_2 = dt * dt;
 	float dt_3 = dt_2 * dt;
 	float dt_4 = dt_3 * dt;

 	//Modify the F matrix so that the time is integrated
 	ekf_.F_(0, 2) = dt;
 	ekf_.F_(1, 3) = dt;

 	//set the process covariance matrix Q
 	ekf_.Q_ = MatrixXd(4, 4);
 	ekf_.Q_ <<  dt_4/4*noise_ax, 0, dt_3/2*noise_ax, 0,
 			   0, dt_4/4*noise_ay, 0, dt_3/2*noise_ay,
 			   dt_3/2*noise_ax, 0, dt_2*noise_ax, 0,
 			   0, dt_3/2*noise_ay, 0, dt_2*noise_ay;

  ekf_.Predict();

  /*****************************************************************************
   *  Update
   ****************************************************************************/

  /**
     * Use the sensor type to perform the update step.
     * Update the state and covariance matrices.
   */

  if (measurement_pack.sensor_type_ == MeasurementPackage::RADAR) {
    // Radar updates

    ekf_.R_ = R_radar_;
    ekf_.H_ = tools.CalculateJacobian(ekf_.x_);

    //measurement update
	  ekf_.UpdateEKF(measurement_pack.raw_measurements_);

  } else {
    // Laser updates
    ekf_.R_ = R_laser_;
    ekf_.H_ = H_laser_;

    //measurement update
	  ekf_.Update(measurement_pack.raw_measurements_);
  }

  // print the output
  //cout << "x_ = " << ekf_.x_ << endl;
  //cout << "P_ = " << ekf_.P_ << endl;
}
