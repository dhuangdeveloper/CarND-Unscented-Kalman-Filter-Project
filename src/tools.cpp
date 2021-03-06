#include <iostream>
#include "tools.h"

using Eigen::VectorXd;
using Eigen::MatrixXd;
using std::vector;

//Tools::Tools() {}

//Tools::~Tools() {}

VectorXd Tools::CalculateRMSE(const vector<VectorXd> &estimations,
                              const vector<VectorXd> &ground_truth) {
  // Calculate the RMSE here
  // Code is from the one I submitted in quiz, with some portion coming from prefilled-code in the quiz
  VectorXd rmse(4);
  rmse << 0,0,0,0;
  // check the validity of the following inputs:
  // * the estimation vector size should not be zero
  //  * the estimation vector size should equal ground truth vector size
  if ( (estimations.size()==0) || (ground_truth.size()==0) || (ground_truth.size()!= estimations.size()))
  {
    return rmse;
  }
  //accumulate squared residuals
  
  for(int i=0; i < estimations.size(); ++i){
	  // ... your code here
	  VectorXd residual(4);
	  residual = estimations[i] - ground_truth[i];
	  residual = residual.array()*residual.array(); 
	  rmse += residual;
  }
  rmse = rmse / estimations.size();
  rmse = rmse.array().sqrt();
  return rmse;
}


float Tools::FoldAngle(float angle, float center_angle){
  // Fold angle into center_angle-pi, center_angle+pi range
  const float pi = 3.14159265;
  if (angle > center_angle+pi) {    
    //cout << angle;
    angle = angle - ceil((angle-pi-center_angle) / pi/2) * 2 * pi;    
    //cout << "->" << angle << endl;
  }
  
  if (angle < center_angle-pi) {        
    //cout << angle;
    angle = angle + ceil((center_angle-pi-angle) / pi/2) * 2 * pi;    
    //cout << "->" << angle << endl;
  }      
  return angle;
}