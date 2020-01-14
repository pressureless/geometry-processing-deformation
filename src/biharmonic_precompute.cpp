#include "biharmonic_precompute.h"
#include <igl/min_quad_with_fixed.h>
#include <igl/cotmatrix.h> 
#include <igl/massmatrix.h>  

void biharmonic_precompute(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const Eigen::VectorXi & b,
  igl::min_quad_with_fixed_data<double> & data)
{
    // REPLACE WITH YOUR CODE
    Eigen::SparseMatrix<double> L;
	  igl::cotmatrix(V, F, L);

	  Eigen::SparseMatrix<double> M;
	  igl::massmatrix(V, F, igl::MASSMATRIX_TYPE_DEFAULT, M);

    Eigen::SparseMatrix<double> MM(M.rows(), M.cols());
    MM.setZero();

	  std::vector<Eigen::Triplet<double>> tripletList;
    for (int i = 0; i < M.rows(); ++i)
    {
    	tripletList.push_back(Eigen::Triplet<double>(i, i, 1/(M.diagonal()(i)) ));
    } 
	  MM.setFromTriplets(tripletList.begin(), tripletList.end());
    

	  Eigen::SparseMatrix<double> Q = L.transpose()*MM*L;
  
	  Eigen::SparseMatrix<double> AeqSparse(0, 0); 

	  bool result = igl::min_quad_with_fixed_precompute(Q, b, AeqSparse, false, data);  
}

