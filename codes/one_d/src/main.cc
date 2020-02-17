#include "space_discr.h"
#include "time_stepper.h"
#include "pde_weights.h"
#include "routines.h" 
#include <algorithm> 
#include <iostream>
#include <chrono>

int main(int argc, char** argv) {
  
  double T = 1; 
  size_t N_0 = 50; /// (N_0+1) nodes on the coarsest mesh
  size_t Levels = 1; /// Number of levels in mesh hierarchy 
  double a=0, b=70; /// Spatial domain [a,b]; 

  Mesh::meshUniformRefinement mesh_structure(N_0, Levels, a, b);

  std::string file_name = "solution.txt"; 
  std::ofstream file(file_name); 
  std::ofstream file_tmp("50pmodmodmod.txt"); 
  
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
    
    const double dt = std::pow(h, 2); 
    const int num_steps = T/dt; 

    /* Galerkin matrices assembly and transformation to Sparse format */ 
    std::cout << "Assembling Galerkin matrices. " << std::endl;

    TripletForm L_Triplet =
      GalerkinMatricesAssembly::getStiffnessMatrixWeighted_TripletForm(*mesh_p, PdeWeights::alpha_x); 
    TripletForm Mr_Triplet = 
      GalerkinMatricesAssembly::getMassMatrixWeighted_TripletForm(*mesh_p, PdeWeights::gamma_x); 
    TripletForm M_Triplet = 
      GalerkinMatricesAssembly::getMassMatrixWeighted_TripletForm(*mesh_p, [](double x){return 1;}); 
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
    A.setFromTriplets(Mr_Triplet.begin(), Mr_Triplet.end());
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
    
    /* Boundary condtions */
    imposeZeroDirichletBoundaryConditions(M, N_dofs); 
    imposeZeroDirichletBoundaryConditions(A, N_dofs); 
    imposeZeroDirichletBoundaryConditions(Ma, N_dofs); 
    imposeZeroDirichletBoundaryConditions(L, N_dofs);
    /* ----END---- */

    SincMatrix Q(h, std::move(L), std::move(M)); 	
   
    /* Time stepping and solving  */ 
    thetaSchemeTimeStepper<SincMatrix,testbool>
      time_stepper(Q, std::move(A), std::move(Ma), dt, h, N_dofs); 
    std::cout << "Stepping in time." << std::endl;
    auto start_time = std::chrono::high_resolution_clock::now(); 
    for (int m=0;m<num_steps; m++){
      u_h = time_stepper.discreteEvolutionOperator(u_h);
    }
    auto end_time = std::chrono::high_resolution_clock::now(); 
    std::cout << std::chrono::duration<double>(end_time-start_time).count()
      << std::endl; 
    /* ----END---- */

    Eigen::VectorXd BS_price = DigitalBS(*mesh_p, K, r, sigma, 0, T);
    /* Writing results to text files */
    std::cout << "Writing data to text file." << std::endl;
    for (int i=0; i<u_h.size() ; i++)
    	file << (*mesh_p)[i] << "," << u0[i] << "," << u_h[i] <<
	  "," << BS_price[i] << std::endl; 
    /* ----END---- */
    
    std::cout << "Solving ended!" << std::endl;
    std::cout << "--------------" << std::endl;
//    std::for_each(Q.kappa_.begin(), Q.kappa_.end(), [&](double arg){
//	file_tmp << arg << std::endl; }); 
}

  return 0;
}
