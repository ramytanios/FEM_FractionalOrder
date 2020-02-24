#include "galerkin_mat.h"
#include "pde_weights.h"
#include "space_discr.h"
#include <iostream>

Galerkin_Mat::Galerkin_Mat(size_t s, std::shared_ptr<Eigen::VectorXd> mesh_p): s_(s){ 
      M_.resize(s,s); L_.resize(s,s); 
      A_.resize(s,s); Ma_.resize(s,s);
      mesh_p_ = mesh_p;
    }

void Galerkin_Mat::setup_galerkin_mat(){ 

  TripletForm L_Triplet =
    GalerkinMatricesAssembly::getStiffnessMatrixWeighted_TripletForm(*mesh_p_, PdeWeights::alpha_x); 
  TripletForm Mr_Triplet = 
  	GalerkinMatricesAssembly::getMassMatrixWeighted_TripletForm(*mesh_p_, PdeWeights::gamma_x); 
  TripletForm M_Triplet = 
  	GalerkinMatricesAssembly::getMassMatrixWeighted_TripletForm(*mesh_p_, [](double x){return 1;}); 
  TripletForm B_Triplet =
    GalerkinMatricesAssembly::getCrossMatrixWeighted_TripletForm(*mesh_p_, PdeWeights::beta_x); 
  TripletForm Ma_Triplet = 
    GalerkinMatricesAssembly::getMassMatrixWeighted_TripletForm(*mesh_p_, [&](double x){return 1;});
 
  L_.setFromTriplets(L_Triplet.begin(), L_Triplet.end());
  A_.setFromTriplets(B_Triplet.begin(), B_Triplet.end());
  A_.setFromTriplets(Mr_Triplet.begin(), Mr_Triplet.end());
  M_.setFromTriplets(M_Triplet.begin(), M_Triplet.end());
  Ma_.setFromTriplets(Ma_Triplet.begin(), Ma_Triplet.end());
}
 
void Galerkin_Mat::impose_dirichlet(){
   imposeZeroDirichletBoundaryConditions(M_, s_);
   imposeZeroDirichletBoundaryConditions(L_, s_);
   imposeZeroDirichletBoundaryConditions(A_, s_);
   imposeZeroDirichletBoundaryConditions(Ma_, s_);
}

std::shared_ptr<SP_MAT> Galerkin_Mat::M(){ return std::make_shared<SP_MAT>(M_);}
std::shared_ptr<SP_MAT> Galerkin_Mat::L(){ return std::make_shared<SP_MAT>(L_);}
std::shared_ptr<SP_MAT> Galerkin_Mat::A(){ return std::make_shared<SP_MAT>(A_);}
std::shared_ptr<SP_MAT> Galerkin_Mat::Ma(){ return std::make_shared<SP_MAT>(Ma_);}
