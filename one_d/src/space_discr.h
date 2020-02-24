#ifndef MATRICES_ASSEMBLY_HPP
#define MATRICES_ASSEMBLY_HPP

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <Eigen/IterativeLinearSolvers>
#include <numeric>
#include "params.h"
#include <iostream>

namespace Mesh{ /* namespace Mesh */
  /**
   * @class meshUniformRefinement
   *
   * @brief This class holds a sequence of refining 1D meshes according
   *	    to the following structure: N0 is fixed. L is the levels that
   *	    can be the set {0,1,...,L-1}. The mesh sizes are then given by
   *	    1/N0*2^{-L} and the total number of nodes are give by 1+N0*2^{L}. 
   *	    Note that N0+1 is the total number of nodes in the coarsest mesh.
   *        The constructed meshes cover the interval [a,b].
   *
   */
  class meshUniformRefinement{ 
    public: 
      meshUniformRefinement(size_t N_0, size_t numOfLevels, double a, double b){ 
	h_.resize(numOfLevels); 
	N_dofs_.resize(numOfLevels); 
	meshes_.resize(numOfLevels);
	std::vector<double> L(numOfLevels); 
	std::iota(L.begin(), L.end(), 0); /// Refinement levels
	/// mesh sizes: 1/N0 * 2^{-L}
	std::transform(L.begin(), L.end(), h_.begin(), [&](double level){
	    return 1.0/N_0 * std::pow(2,-level);} );	
        /// Total number of nodes: 1+N0*2^{L}
	std::transform(L.begin(), L.end(), N_dofs_.begin(), [&](double level){
	    return 1 + N_0 * std::pow(2,level);} );
	std::transform(N_dofs_.begin(), N_dofs_.end(), meshes_.begin(), [&](size_t Ndofs){
    	    return Eigen::VectorXd::LinSpaced(Ndofs, a, b); } ); 
      }
      
      /// Returns a pointer to the mesh at a given level 
      std::shared_ptr<Eigen::VectorXd> getMesh(size_t level){
	  return std::make_shared<Eigen::VectorXd>(meshes_.at(level));
	}
      
      std::vector<double> h_; 
      std::vector<size_t> N_dofs_; 
   
    private: 
      std::vector<Eigen::VectorXd> meshes_; 
  };

} /* namespace Mesh */ 

/**
 * @brief Imposing 0 Dirichlet boundary conditions on the galerkin matrices 
 *	  by setting the corresponding rows/columns to 0, and the diagonal entries to 1
 *
 * @tparam SparseMatrix_t the type of the sparse matrix to be modified
 * @param A the matrix to be modified
 * @param N the number of rows/columns of the matrix A
 *
 */
template<typename SparseMatrix_t> 
void imposeZeroDirichletBoundaryConditions(SparseMatrix_t& A, size_t N){
  /* Using coeffRef is not expensive since entries already exist. */ 
  /* Exploit the tridiagonal structure of the galerkin matrices.  */
    A.coeffRef(N-2,N-1) = 0; A.coeffRef(N-1,N-2) = 0; A.coeffRef(N-1,N-1) = 1; 
    A.coeffRef(0,0) = 1; A.coeffRef(0,1) = 0; A.coeffRef(1,0) = 0; 
 
}
/* */

/** 
 * @function ConjGradDescentSolver
 *
 * @brief Solves a linear system of equations
 *        Ax=b using the iterative CGD
 *
 * @tparam A the lhs matrix A
 * @tparam b the rhs vector b
 * @tparam x0 the initial guess for x
 * @return the solution x
 *
 */ 
inline Eigen::VectorXd ConjGradDescentSolver(Eigen::SparseMatrix<double>& A, 
    Eigen::VectorXd& b, const Eigen::VectorXd& x0){

  const double tol = 1e-10; 
  const double maxIterations = 1000; 
  size_t n=0; 
  
  Eigen::VectorXd x = Eigen::VectorXd::Zero(x0.size()); 
  Eigen::VectorXd old_resid = b - A*x; 

  if (old_resid.norm() < tol)
    return x;
  else{
    Eigen::VectorXd p = old_resid;
    Eigen::VectorXd new_resid = old_resid; 
    while (new_resid.norm() > tol){
      if (new_resid.norm() > 1e3)
	throw std::overflow_error("The residual is diverging in CGD!");
      if (n > maxIterations)
	throw std::overflow_error("Maximum CGD iterations reached!");
      old_resid.swap(new_resid);
      double alpha = old_resid.dot(old_resid) / 
			 old_resid.dot(A*p);
      x += alpha * p; 
      new_resid = old_resid - alpha * A*p; 
      double beta_tmp = new_resid.dot(new_resid) / old_resid.dot(old_resid); 
      p = new_resid + beta_tmp * p; 
      n++;
    }
    return x;
  }
}

/** 
 * @function SparseDiagonal
 *
 * @brief Constructs a diagonal matrix in sparse format
 *	  of type Eigen::SparseMatrix from a diagonal 
 *	  vector
 *
 * @tparam M the sparse matrix to be constructed
 * @tparam diag the diagonal vector
 *
 */ 
inline void SparseDiagonal(Eigen::SparseMatrix<double>& M, 
    				Eigen::VectorXd& diag){
  //std::vector<Eigen::Triplet<double>> triplets(diag.size()); 
  for (int i=0; i<diag.size(); i++)
    M.coeffRef(i,i) = diag(i);
}

/** 
 * @function PrecondConjGradDescentSolver
 *
 * @brief Solves a linear system of equations
 *        Ax=b using the iterative PrecondCGD
 *	  where the precondioner is D^{-1} and
 *	  D is the diagonal of A. 
 *
 * @tparam A the lhs matrix A
 * @tparam b the rhs vector b
 * @tparam x0 the initial guess for x
 * @return the solution x
 *
 */ 
inline Eigen::VectorXd PrecondConjGradDescentSolver(Eigen::SparseMatrix<double>& A, 
    Eigen::VectorXd& b, const Eigen::VectorXd& x0){

  const double tol = 1e-10; 
  const double maxIterations = 1000; 
  size_t n=0; 

  Eigen::SparseMatrix<double> P(b.size(), b.size()); // The preconditionner 
  Eigen::VectorXd Pdiag = 1 / A.diagonal().array();
  SparseDiagonal(P, Pdiag); 

  Eigen::SparseMatrix<double> Amod = P*A; 
  Eigen::VectorXd bmod = P*b; 
  
  Eigen::VectorXd x = Eigen::VectorXd::Zero(x0.size()); 
  Eigen::VectorXd old_resid = bmod - Amod*x; 

  if (old_resid.norm() < tol)
    return x;
  else{
    Eigen::VectorXd p = old_resid;
    Eigen::VectorXd new_resid = old_resid; 
    while (new_resid.norm() > tol){
      if (new_resid.norm() > 1e3)
	throw std::overflow_error("The residual is diverging in CGD!");
      if (n > maxIterations)
	throw std::overflow_error("Maximum CGD iterations reached!");
      old_resid.swap(new_resid);
      double alpha = old_resid.dot(old_resid) / 
			 old_resid.dot(Amod*p);
      x += alpha * p; 
      new_resid = old_resid - alpha * Amod*p; 
      double beta_tmp = new_resid.dot(new_resid) / old_resid.dot(old_resid); 
      p = new_resid + beta_tmp * p; 
      n++;
    }
    return x;
  }
}


/**
 * @function getStiffnessMatrixWeighted_TripletForm
 * @function getCrossMatrixWeighted_TripletForm
 * @function getMassMatrixWeighted_TripletForm
 * @function getLoadVector
 *
 * @brief Assembly of the Galerkin matrices: stiffness, cross and mass matrices
 * 	  for general weights alpha(x), beta(x), gamma(x), and the load vector
 * 	  for a general function a load function f(x)
 *
 * Local trapezoidal rule is used for the cross and mass matrices. 
 * Local midpoint rule is used for the stiffness matrix.
 * Local trapezoidal rule is used for the load vector. 
 *
 * @tparam FUNCTOR the type of the coefficient objects
 * @param mesh the 1d mesh
 * @param alpha,beta,gamma,f the coefficients
 * @return the galerkin matrices in triplet formet of type TripletForm
 *
 */
using SparseMatrix =  Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;
using TripletForm = std::vector<Triplet>; 

namespace GalerkinMatricesAssembly{ /* namespace GalerkinMatricesAssembly */

template <typename FUNCTOR>
TripletForm getStiffnessMatrixWeighted_TripletForm(const Eigen::VectorXd &mesh, FUNCTOR alpha) {
  
  TripletForm triplets;
  unsigned M = mesh.size() - 1;
  triplets.reserve(3 * M + 1);
  
  /* Diagonal Entries */
  double diagEntry = 1.0/(mesh(1)-mesh(0)) * alpha(0.5*(mesh(0)+mesh(1))); 
  triplets.push_back(Triplet(0,0,diagEntry));
  diagEntry = 1.0/(mesh(M)-mesh(M-1)) * alpha(0.5*(mesh(M)+mesh(M-1))); 
  triplets.push_back(Triplet(M,M,diagEntry));  
  for (int i=1; i<M; i++){
	  double dx_left = mesh(i) - mesh(i-1);
	  double dx_right = mesh(i+1) - mesh(i); 
	  diagEntry = 1.0/dx_left * alpha(0.5*(mesh(i)+mesh(i-1)))
		  	+ 1.0/dx_right * alpha(0.5*(mesh(i+1)+mesh(i)));
	  triplets.push_back(Triplet(i,i,diagEntry)); 
  }
  
  /* Off-Diagonal Entries */
  for (int i=0; i<M; i++){ 
	  double dx_right = mesh(i+1) - mesh(i); 
	  double offDiag = -1/dx_right 
		  		* alpha(0.5*(mesh(i)+mesh(i+1))); 
	  triplets.push_back(Triplet(i,i+1,offDiag)); 
	  triplets.push_back(Triplet(i+1,i,offDiag));
  }
  return triplets;
}

template<typename FUNCTOR> 
TripletForm getCrossMatrixWeighted_TripletForm(const Eigen::VectorXd& mesh, FUNCTOR beta){
  TripletForm triplets; 
  unsigned M = mesh.size() - 1; 
  triplets.reserve(2 * M + 1); 

  /* Diagonal Entries */
  double diagEntry = -0.5 * beta(mesh(0)); 
  triplets.push_back(Triplet(0, 0, diagEntry));
  diagEntry = 0.5 * beta(mesh(M)); 
  triplets.push_back(Triplet(M, M, diagEntry)); 

  /* Off-Diagonal Entries */
  for (int i=0; i<M; i++){ 
    triplets.push_back(Triplet(i, i+1, 0.5*beta(mesh(i))));
    triplets.push_back(Triplet(i+1, i, -0.5*beta(mesh(i+1))));
  }
  return triplets;
}

template <typename FUNCTOR>
TripletForm getMassMatrixWeighted_TripletForm(const Eigen::VectorXd &mesh, FUNCTOR gamma) {
  TripletForm triplets;
  unsigned M = mesh.size() - 1;
  triplets.reserve(M + 1);

  /* Diagonal Entries*/
  double diagEntry = 0.5 * gamma(mesh(0)) * (mesh(1)-mesh(0)); 
  triplets.push_back(Triplet(0,0,diagEntry));
  
  diagEntry = 0.5 * gamma(mesh(M)) * (mesh(M) - mesh(M-1)); 
  triplets.push_back(Triplet(M,M,diagEntry));
  
  for (int i=1; i<M; i++){
	  double dx_right = mesh(i+1) - mesh(i);
	  double dx_left = mesh(i) - mesh(i-1); 
	  diagEntry = 0.5 * gamma(mesh(i)) * (dx_right+dx_left);
	  triplets.push_back(Triplet(i,i,diagEntry));
  }
  return triplets;
}

template<typename FUNCTOR>
Eigen::VectorXd getLoadVector(const Eigen::VectorXd& mesh, FUNCTOR f){
  size_t N = mesh.size();
  Eigen::VectorXd res(N);
 
  res[0] = 0.5 * (mesh(1) - mesh(0)) * f(mesh(0)); 
  for (int i=1; i<N-1; i++)
    res[i] = 0.5 * (mesh(i+1) - mesh(i-1)) * f(mesh(i)); 
  res[N-1] = 0.5 * (mesh(N-1) - mesh(N-2)) * f(mesh(N-1)); 

  return res;
}

} /* namespace GalerkinMatricesAssembly */


#endif
