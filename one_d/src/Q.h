#ifndef SINC_MATRIX_HPP
#define SINC_MATRIX_HPP

#include <Eigen/Sparse>
#include <Eigen/SparseCholesky>
#include <Eigen/SVD>
#include <Eigen/IterativeLinearSolvers>
#include <numeric>
#include "params.h"
#include <iostream>

/**
 * @class SincMatrix
 *
 * @brief Representing the sinc matrix Q, without storing its elements explicitly.
 *        Accessing the matrix Q through matrix-vector product by overloading the * operator.
 *        Given f \in \mathbb{R}^n, the * method returns Q*f.
 *
 */
class SincMatrix{
  public:
    SincMatrix() = delete; 
    SincMatrix(const SincMatrix&);
    SincMatrix(double, const Eigen::SparseMatrix<double>& L, const Eigen::SparseMatrix<double>& M);
    SincMatrix(double, Eigen::SparseMatrix<double>&& L, Eigen::SparseMatrix<double>&& M);
    Eigen::VectorXd operator* (Eigen::VectorXd& f) const; 
  private:
    double k_; 
    Eigen::SparseMatrix<double> L_, M_;
    Eigen::VectorXd K_; 

    template<typename SincMatType, bool is_theta_zero>
    friend class thetaSchemeTimeStepper;
};

SincMatrix::SincMatrix(const SincMatrix& other):
k_(other.k_), L_(other.L_),
M_(other.M_),
K_(other.K_){}

SincMatrix::SincMatrix(double h, const Eigen::SparseMatrix<double>& L,
                       const Eigen::SparseMatrix<double>& M):
M_(M),
L_(L) {
    k_ = -1/(std::log(h));
    K_ = Eigen::VectorXd::LinSpaced(
                                    ceil(M_PI*M_PI / (4.*beta*k_*k_)) + ceil(M_PI*M_PI / (4.*(1-beta)*k_*k_)) + 1,
                                    -ceil(M_PI*M_PI / (4.*beta*k_*k_)),  ceil(M_PI*M_PI / (4.*(1-beta)*k_*k_)) );
}

SincMatrix::SincMatrix(double h, Eigen::SparseMatrix<double>&& L,
                       Eigen::SparseMatrix<double>&& M):
M_(std::move(M)),
L_(std::move(L)) {
    k_ = -1/(std::log(h));
    K_ = Eigen::VectorXd::LinSpaced(
                                    ceil(M_PI*M_PI / (4.*beta*k_*k_)) + ceil(M_PI*M_PI / (4.*(1-beta)*k_*k_)) + 1,
                                    -ceil(M_PI*M_PI / (4.*beta*k_*k_)),  ceil(M_PI*M_PI / (4.*(1-beta)*k_*k_)) );
}
/**
 * @member operator*
 *
 * @brief Performing the matrix-vector product Q*f
 *
 * @param f the vector f
 * @return the product Q*f
 *
 */
Eigen::VectorXd SincMatrix::operator* (Eigen::VectorXd& f) const{
    if (M_.cols() != f.size())
        throw std::length_error("Matrix-vector sizes do not match!");
    
    Eigen::VectorXd res = Eigen::VectorXd::Zero(f.size());
    
    double Cb = M_PI/(2*std::sin(M_PI*beta));
    
    Eigen::ConjugateGradient<Eigen::SparseMatrix<double>, Eigen::Lower|Eigen::Upper,
    Eigen::IdentityPreconditioner> PCGD_Solver;
    PCGD_Solver.setMaxIterations(1000);
    PCGD_Solver.setTolerance(1e-10);
    
    Eigen::SparseLU<Eigen::SparseMatrix<double>> SparseLU_Solver;
    
    std::for_each(K_.data(), K_.data()+K_.size(), [&](double arg){
        SparseLU_Solver.compute(M_+std::exp(2*arg*k_)*L_);
        res = res + (k_/Cb)*std::exp(2*beta*arg*k_)*SparseLU_Solver.solve(f); });
    
    //  std::for_each(K_.data(), K_.data()+K_.size(), [&](double arg){
    //       SparseLU_Solver.compute(M_+std::exp(2*arg*k_)*L_);
    //       res = res + (k_/Cb)*std::exp(2*beta*arg*k_)*SparseLU_Solver.solve(f); });
    
    //  std::for_each(K_.data(), K_.data()+K_.size(), [&](double arg){
    //      Eigen::SparseMatrix<double> lhs = M_+std::exp(2*arg*k_)*L_;
    //      res = res + (k_/Cb)*std::exp(2*beta*arg*k_)
    //       *ConjGradDescentSolver(lhs,f,Eigen::VectorXd::Zero(f.size())); });
    
    //  std::for_each(K_.data(), K_.data()+K_.size(), [&](double arg){
    //      Eigen::SparseMatrix<double> lhs = M_+std::exp(2*arg*k_)*L_;
    //      res = res + (k_/Cb)*std::exp(2*beta*arg*k_)
    //       *PrecondConjGradDescentSolver(lhs,f,Eigen::VectorXd::Zero(f.size())); });
    
    return res;
}

#endif 
