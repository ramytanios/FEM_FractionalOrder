#include "space_discr.h"
#include "time_stepper.h"
#include "pde_weights.h"
#include "routines.h" 
#include <algorithm> 
#include <iostream>

int main(int argc, char** argv) {
  
  double T = 1; 
  size_t N_0 = 4; /// (N_0+1) nodes on the coarsest mesh
  size_t Levels = 6; /// Number of levels in mesh hierarchy 
  double a=0, b=1; /// Spatial domain [a,b]; 

  Mesh::meshUniformRefinement mesh_structure(N_0, Levels, a, b);

  std::vector<double> l2_error(Levels); 
  std::vector<double> linf_error(Levels); 

  std::string file_name = "error_"+std::to_string(theta)+".txt"; 
  std::ofstream file(file_name); 
  
  /// Looping over mesh hierarchy levels
  for (int curr_level=0; curr_level<Levels; curr_level++){

    /// Get mesh information: mesh, mesh_size, number of nodes. 
    std::shared_ptr<Eigen::VectorXd> mesh_p(mesh_structure.getMesh(curr_level));
    double h = mesh_structure.h_.at(curr_level); 
    size_t N_dofs = mesh_structure.N_dofs_.at(curr_level); 
   
    std::cout << "Solving with N_dofs = " << N_dofs <<
      ", h = " << h <<
      ", theta = " << theta <<
      " and beta = " << beta << std::endl;
    
    const double dt = std::pow(h,1); 
    const int num_steps = T/dt; 

    /* Galerkin matrices assembly and transformation to Sparse format */ 
    std::cout << "Assembling Galerkin matrices. " << std::endl;

    TripletForm L_Triplet =
      GalerkinMatricesAssembly::getStiffnessMatrixWeighted_TripletForm(*mesh_p, PdeWeights::alpha_x); 
    TripletForm M_Triplet = 
      GalerkinMatricesAssembly::getMassMatrixWeighted_TripletForm(*mesh_p, PdeWeights::gamma_x); 
    TripletForm B_Triplet =
      GalerkinMatricesAssembly::getCrossMatrixWeighted_TripletForm(*mesh_p, PdeWeights::beta_x); 
     TripletForm Ma_Triplet =                                                                                              
       GalerkinMatricesAssembly::getMassMatrixWeighted_TripletForm(*mesh_p, [&](double x){return sigma*std::sqrt(x);});  
    SparseMatrix A(N_dofs, N_dofs);
    SparseMatrix L(N_dofs, N_dofs);
    SparseMatrix M(N_dofs, N_dofs);
    SparseMatrix Ma(N_dofs, N_dofs); 
    L.setFromTriplets(L_Triplet.begin(), L_Triplet.end());
    A.setFromTriplets(B_Triplet.begin(), B_Triplet.end());
    M.setFromTriplets(M_Triplet.begin(), M_Triplet.end());
    Ma.setFromTriplets(Ma_Triplet.begin(), Ma_Triplet.end());
    /* ----END---- */ 

    /* Projecting initial condition on FEM space */ 
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M); 
    Eigen::VectorXd g = GalerkinMatricesAssembly::getLoadVector(*mesh_p, PdeWeights::initial_cond); 
    Eigen::VectorXd u0 = solver.solve(g); 
    /* ----END----*/ 

    Eigen::VectorXd u_h(u0.size()); 
    u_h = u0; 

    /* Dirichlet boundary condtions */
    imposeZeroDirichletBoundaryConditions(M, N_dofs); 
    imposeZeroDirichletBoundaryConditions(A, N_dofs); 
    imposeZeroDirichletBoundaryConditions(L, N_dofs);
    imposeZeroDirichletBoundaryConditions(Ma, N_dofs);
    /* ----END---- */

    SincMatrix Q(h, std::move(L), std::move(M)); 
   
    /* Time stepping and solving  */ 
    thetaSchemeTimeStepper<SincMatrix,testbool> 
      time_stepper(Q, std::move(A), std::move(Ma), dt, h, N_dofs); 
    std::cout << "Stepping in time." << std::endl;
    for (int m=0;m<num_steps; m++){
      u_h = time_stepper.discreteEvolutionOperator(u_h);
    }
    /* ----END---- */

    Eigen::VectorXd u_ref = ReadNthLine("uref_new_beta1o2.txt",curr_level, N_dofs);
    Eigen::VectorXd e_h = u_ref - u_h; 
    l2_error[curr_level] = l2_norm(e_h, h); 
    linf_error[curr_level] = e_h.lpNorm<Eigen::Infinity>(); 

    /* Writing results to text files */
    std::cout << "Writing data to text file." << std::endl;
    file << h << "," << l2_norm(e_h, h) 
      	      << "," << e_h.lpNorm<Eigen::Infinity>() << std::endl; 
    /* ----END---- */
    
    std::cout << "L2-norm of error " << l2_norm(e_h, h) << " and " 
              << "Loo-norm of error " << e_h.lpNorm<Eigen::Infinity>() << std::endl;
    std::cout << "Solving ended!" << std::endl;
    std::cout << "--------------" << std::endl;
}

  std::vector<double> l2_error_log(Levels);
  std::transform(l2_error.begin(), l2_error.end(), l2_error_log.begin(), log_op()); 
  
  std::vector<double> linf_error_log(Levels); 
  std::transform(linf_error.begin(), linf_error.end(), linf_error_log.begin(), log_op()); 
  
  std::vector<double> h_log(Levels); 
  std::transform(mesh_structure.h_.begin(), mesh_structure.h_.end(), h_log.begin(), log_op());

  std::pair<double, double> coeff_l2 = linear_fit(l2_error_log, h_log); 
  std::pair<double, double> coeff_linf = linear_fit(linf_error_log, h_log); 

  std::cout << "The algebraic rate of convergence in the L2 norm is " << coeff_l2.second << std::endl; 
  std::cout << "The algebraic rate of convergence in the Linf norm  is " << coeff_linf.second << std::endl; 

  file << coeff_l2.second << "," << coeff_linf.second << std::endl;
  return 0;
}
