#include "arap_single_iteration.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/polar_svd3x3.h>

void arap_single_iteration(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::SparseMatrix<double> & K,
  const Eigen::MatrixXd & bc,
  Eigen::MatrixXd & U)
{
    // REPLACE WITH YOUR CODE 
    Eigen::MatrixXd Beq = Eigen::MatrixXd::Zero(0,0); 
    Eigen::MatrixXd B; 

    Eigen::VectorXi known = data.known;
    Eigen::MatrixXd UU = U;
    for (int i = 0; i < known.size(); ++i)
    {
        // std::cout<<"previous i:"<<i<<", value:"<<UU.row(known(i))<<", after:"<<bc.row(i)<<std::endl;
        UU.row(known(i)) = bc.row(i);
    }
    // std::cout<<"UU:"<<UU<<std::endl;

    Eigen::MatrixXd Ct = U.transpose()*K;
    Eigen::MatrixXd C = Ct.transpose();
    Eigen::Matrix3d Ri(3,3);
    Eigen::MatrixXd R(3*data.n,3); 
    for (int i = 0; i < data.n; ++i)
    { 
    	Eigen::Matrix3d Ci = C.block(3*i, 0, 3, 3); 
    	igl::polar_svd3x3<Eigen::Matrix3d>(Ci, Ri); 
    	R.block(3*i, 0, 3, 3) = Ri;
    }
    B = K * R;
    // 
    Eigen::MatrixXd BB = Eigen::MatrixXd::Zero(data.n, 1); 
    igl::min_quad_with_fixed_solve(data, B, bc, Beq, U); 
}
