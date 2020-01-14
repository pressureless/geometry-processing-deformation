#include "biharmonic_solve.h"
#include <igl/min_quad_with_fixed.h>
#include "biharmonic_precompute.h"

void biharmonic_solve(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & D)
{
    // REPLACE WITH YOUR CODE
    D = Eigen::MatrixXd::Zero(data.n,3);
    Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0,0); 
    Eigen::MatrixXd B = Eigen::MatrixXd::Zero(data.n, 1); 
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, D); 
}

