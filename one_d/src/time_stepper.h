#ifndef TIME_STEPPER_HPP
#define TIME_STEPPER_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <cassert>
#include <iostream>
#include "params.h" 

using SparseMatrix = Eigen::SparseMatrix<double>; 

/**
 * @class thetaSchemeTimeStepper
 *
 * @brief Performing the \theta scheme time stepping.
 *
 * @tparam SincMatType the type of sinc matrix
 * @tparam is_beta_zero a boolean to distinguish between 
 *         fractional and non-fractional cases
 *
 */
template<typename SincMatType, bool is_beta_zero> 
class thetaSchemeTimeStepper{
  public: 
    thetaSchemeTimeStepper() = delete;

    /// Move and copy constructors 
    thetaSchemeTimeStepper(SincMatType&& Q, SparseMatrix&& A, SparseMatrix&& Ma, 
	double dt, double h, double N_dofs);
    thetaSchemeTimeStepper(const SincMatType& Q, const SparseMatrix& A, const SparseMatrix& Ma,
	double dt, double h, double N_dofs);

    /// \Psi^{\tau} (u(t)) = u(t+\tau)
    template<bool it_is = is_beta_zero>
      typename std::enable_if<it_is,Eigen::VectorXd>::type
      discreteEvolutionOperator(Eigen::VectorXd&) const; 
    template<bool it_is = is_beta_zero>
      typename std::enable_if<!it_is,Eigen::VectorXd>::type 
      discreteEvolutionOperator(Eigen::VectorXd&) const;

    virtual ~thetaSchemeTimeStepper() = default;
  private: 
    double dt_, h_; int N_dofs_;  
    SparseMatrix A_; /// The stiffness matrix w/o MQL
    SparseMatrix L_; /// FEM approx. of \mathcal{L}
    SparseMatrix M_; /// The mass matrix
    SparseMatrix Ma_;
    SincMatType Q_; /// FEM approx. of \mathcal{L}^{-beta = -(1-\beta)}
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_; /// for non-fractional case
};

template<typename SincMatType,bool is_beta_zero>
thetaSchemeTimeStepper<SincMatType, is_beta_zero>::thetaSchemeTimeStepper(
    SincMatType&& Q, SparseMatrix&& A, SparseMatrix&& Ma,
    double dt, double h, double N_dofs):
  dt_(dt), h_(h), N_dofs_(N_dofs), Q_(std::move(Q)), A_(std::move(A)), 
  L_(std::move(Q.L_)), M_(std::move(Q.M_)), Ma_(std::move(Ma)){   
    solver_.compute(M_+dt_*theta*(A_+L_)); 
  }

template<typename SincMatType, bool is_beta_zero>
thetaSchemeTimeStepper<SincMatType, is_beta_zero>::thetaSchemeTimeStepper(
    const SincMatType& Q, const SparseMatrix& A, const SparseMatrix& Ma, 
    double dt, double h, double N_dofs):
   dt_(dt), h_(h), N_dofs_(N_dofs), Q_(Q), A_(A), 
   L_(Q.L_), M_(Q.M_), Ma_(Ma){
     solver_.compute(M_+dt_*theta*(A_+L_)); 
  }

/** 
 * @member discreteEvolutionOperator
 *
 * @brief Performing the time stepping of the solution
 *
 * @tparam SincMatType the type of the sinc matrix
 * @tparam is_beta_zero indicating fractional/non-fractional cases
 * @param u_prev the solution at time t
 * @return the solution at time t+\tau
 *
 */
template<typename SincMatType, bool is_beta_zero> 
template<bool it_is> typename std::enable_if<it_is,Eigen::VectorXd>::type
thetaSchemeTimeStepper<SincMatType, 
  			is_beta_zero>::discreteEvolutionOperator(
     Eigen::VectorXd& u_prev) const{ 
  Eigen::VectorXd u_next(u_prev.size());
  u_next = solver_.solve( (M_-dt_*(1-theta)*(A_+L_))*u_prev );
  return u_next; 
}

/** 
 * @member discreteEvolutionOperator
 *
 * @overload
 *
 * Conjugate gradient descent is implemented for fractional case
 *
 * @throws std::overflow_error if residuals in the CGD diverge
 *
 */
template<typename SincMatType, bool is_beta_zero> 
template<bool it_is> typename std::enable_if<!it_is,Eigen::VectorXd>::type
thetaSchemeTimeStepper<SincMatType,
  			is_beta_zero>::discreteEvolutionOperator(
     Eigen::VectorXd& u_prev) const{
  /// Conjugate gradient descent algorithm for the non-fractional case
  size_t n = 0; 
  Eigen::VectorXd u_next = Eigen::VectorXd::Zero(u_prev.size()); 
  Eigen::VectorXd tmp = L_ * u_prev;
  Eigen::VectorXd rhs_vec = M_*u_prev - 
    		dt_*(1-theta) * (A_*u_prev + Ma_*(Q_*tmp) ); 
  Eigen::VectorXd tmpp = L_ * u_next;
  Eigen::VectorXd old_resid = rhs_vec - (M_*u_next + theta*dt_*
      			(A_*u_next+Ma_*(Q_*tmpp) ) );
  if (old_resid.norm() < tol)
    return u_next; 
  else{
    Eigen::VectorXd p = old_resid; 
    Eigen::VectorXd new_resid = old_resid;
    while (new_resid.norm() > tol){
      if (new_resid.norm() > 1e3)
	throw std::overflow_error("The residual is diverging in CGD!");
      if (n > maxIterations)
	throw std::overflow_error("Maximum CGD iterations reached!");
      old_resid.swap(new_resid);  
      tmp = L_ * p; 
      Eigen::VectorXd Stiff_times_p = M_*p + 
		theta * dt_ * (A_*p + Ma_*(Q_*tmp)); 
      double alpha = old_resid.dot(old_resid) / 
			old_resid.dot(Stiff_times_p);
      u_next += alpha * p; 
      new_resid = old_resid - alpha * Stiff_times_p;
      double beta_tmp = new_resid.dot(new_resid) / old_resid.dot(old_resid); 
      p = new_resid + beta_tmp * p;
      n++; 
    }
    return u_next; 
  }
}
#endif

