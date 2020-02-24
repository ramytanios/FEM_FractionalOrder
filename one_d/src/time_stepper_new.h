#ifndef TIME_STEPPER_HPP
#define TIME_STEPPER_HPP

#include <Eigen/Sparse>
#include <Eigen/Dense>
#include <Eigen/SparseCholesky>
#include <cassert>
#include <iostream>
#include "params.h"
#include "wrapper.h"

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
    const SparseMatrix A_; /// The stiffness matrix w/o MQL
    const SparseMatrix L_; /// FEM approx. of \mathcal{L}
    const SparseMatrix M_; /// The mass matrix
    const SparseMatrix Ma_; /// Mass matrix weighted by a(x): See report
    const SincMatType Q_; /// FEM approx. of \mathcal{L}^{-beta = -(1-\beta)}
    const MatrixReplacement Lhs_Wrapped_; /// Lhs matrix of time stepping scheme, wrapped. 
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver_; /// for non-fractional case
    Eigen::ConjugateGradient<MatrixReplacement, Eigen::Lower|Eigen::Upper,
        Eigen::IdentityPreconditioner> IterativeSolver_;  /// for fractional case
};

template<typename SincMatType,bool is_beta_zero>
thetaSchemeTimeStepper<SincMatType, is_beta_zero>::thetaSchemeTimeStepper(
    SincMatType&& Q, SparseMatrix&& A, SparseMatrix&& Ma,
                    double dt, double h, double N_dofs):
  dt_(dt), h_(h), N_dofs_(N_dofs),
  Q_(std::move(Q)),
  A_(std::move(A)),
  L_(std::move(Q.L_)),
  M_(std::move(Q.M_)),
  Ma_(std::move(Ma)),
  Lhs_Wrapped_(M_, A_, Ma_, L_, Q_, dt_), IterativeSolver_(Lhs_Wrapped_)
{
    IterativeSolver_.setTolerance(1e-5);
    IterativeSolver_.setMaxIterations(1000);
    solver_.compute(M_+dt_*theta*(A_+L_));
}

template<typename SincMatType, bool is_beta_zero>
thetaSchemeTimeStepper<SincMatType, is_beta_zero>::thetaSchemeTimeStepper(                                             
    const SincMatType& Q, const SparseMatrix& A, const SparseMatrix& Ma,
        		double dt, double h, double N_dofs):
  dt_(dt), h_(h), N_dofs_(N_dofs),
  Q_(Q), 
  A_(A),
  L_(Q.L_),
  M_(Q.M_), 
  Ma_(Ma), 
  Lhs_Wrapped_(M_, A_, Ma_, L_, Q_, dt_), IterativeSolver_(Lhs_Wrapped_)
{
    IterativeSolver_.setTolerance(1e-5);
    IterativeSolver_.setMaxIterations(1000);
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
    Eigen::VectorXd u_next(u_prev.size());
    Eigen::VectorXd tmpp = L_*u_prev;
    u_next = IterativeSolver_.solve(M_*u_prev - (dt_*(1-theta))*
                                    (A_*u_prev+M_*(Q_*tmpp)));
    return u_next;
}
#endif

