#include "arap_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h> 
#include <igl/cotmatrix_entries.h> 
#include <igl/edge_lengths.h>

void arap_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data,
  Eigen::SparseMatrix<double> & K)
{
    // REPLACE WITH YOUR CODE
    Eigen::SparseMatrix<double> L;
	igl::cotmatrix(V, F, L);
    
    Eigen::MatrixXd angle;
    igl::cotmatrix_entries(V, F, angle);

	Eigen::MatrixXd l;
	igl::edge_lengths(V,F,l); 

    Eigen::MatrixXd KK(V.rows(), 3*V.rows());

    for (int i = 0; i < F.rows(); ++i)
    { 
        for (int j = 0; j < 3; ++j)
        {
            int s = (j+1) % 3;
            int t = (j+2) % 3;

            double cos = (l(i,j)*l(i,j) + l(i,s)*l(i,s) - l(i,t)*l(i,t)) / (2 * l(i,j) * l(i,s)); //cos
            double cot = cos/std::sqrt(1 - cos*cos);  //cot
            Eigen::RowVectorXd length = cot * (V.row(F(i,j)) - V.row(F(i,s)));

            for (int k = 0; k < 3; ++k)
            {
                KK(F(i,j), 3*F(i,j) + k) += length(k);
                KK(F(i,s), 3*F(i,j) + k) -= length(k);

                KK(F(i,j), 3*F(i,s) + k) += length(k);
                KK(F(i,s), 3*F(i,s) + k) -= length(k);

                KK(F(i,j), 3*F(i,t) + k) += length(k);
                KK(F(i,s), 3*F(i,t) + k) -= length(k);
            }
        } 
    }
    K = KK.sparseView();

	Eigen::SparseMatrix<double> LL = L * 2;

	Eigen::SparseMatrix<double> AeqSparse(0, 0); 
	bool result = igl::min_quad_with_fixed_precompute(LL, b, AeqSparse, false, data);
	std::cout<<"arap_precompute result:"<<result<<std::endl;
	// std::cout<<"arap_precompute K:"<<K<<std::endl;
}
