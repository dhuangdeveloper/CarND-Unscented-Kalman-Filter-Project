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