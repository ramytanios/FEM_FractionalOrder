#include "space_discr.h"
#include "pde_weights.h"
#include "time_stepper_new.h"
#include "Q.h"
#include "routines.h" 
#include "galerkin_mat.h" 
#include <algorithm> 
#include <iostream>
#include <chrono>

int main(int argc, char** argv) {
  
  double T = 1; 
  size_t N_0 = std::stoi(argv[1]); /// (N_0+1) nodes on the coarsest mesh
  size_t Levels = 1; /// Number of levels in mesh hierarchy 
  double a=0, b=1; /// Spatial domain [a,b]; 

  Mesh::meshUniformRefinement mesh_structure(N_0, Levels, a, b);

  std::string file_name = "solution.txt"; 
  std::ofstream file(file_name); 
  
  /// Looping over mesh hierarchy levels
  for (int curr_level=0; curr_level<Levels; curr_level++){

    /// Get mesh information: mesh, mesh_size, number of nodes. 
    std::shared_ptr<Eigen::VectorXd> mesh_p(mesh_structure.getMesh(curr_level));
    double h = mesh_structure.h_.at(curr_level); 
    size_t N_dofs = mesh_structure.N_dofs_.at(curr_level); 
   
    std::cout << "Solving with N_dofs = " << N_dofs <<
      ", h = " << h << ", theta = " << theta << " and beta = " << beta << std::endl;
    
    const double dt = std::pow(h, 2); 
    const int num_steps = T/dt; 

    /* Galerkin matrices assembly, transformation to Sparse format
     and imposing Dirichlet boundary conditions. */ 
    std::cout << "Assembling Galerkin matrices. " << std::endl;
    Galerkin_Mat GalerkinMatSetup(N_dofs, mesh_p);
    GalerkinMatSetup.setup_galerkin_mat(); 
    GalerkinMatSetup.impose_dirichlet();
    Eigen::SparseMatrix<double> M(*GalerkinMatSetup.M()); 
    Eigen::SparseMatrix<double> L(*GalerkinMatSetup.L()); 
    Eigen::SparseMatrix<double> Ma(*GalerkinMatSetup.Ma()); 
    Eigen::SparseMatrix<double> A(*GalerkinMatSetup.A()); 
    /* ----END---- */ 

    /* Projecting initial condition on FEM space */ 
    Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
    solver.compute(M); 
    Eigen::VectorXd g = GalerkinMatricesAssembly::getLoadVector(*mesh_p, PdeWeights::initial_cond); 
    Eigen::VectorXd u0 = solver.solve(g); 
    /* ----END----*/

    Eigen::VectorXd u_h(u0.size()); 
    u_h = u0;
    
    SincMatrix Q(h, L, M); 	
   
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

    /* Writing results to text files */
    std::cout << "Writing data to text file." << std::endl;
    for (int i=0; i<u_h.size() ; i++)
    	file << (*mesh_p)[i] << "," << u0[i] << "," << u_h[i] << std::endl; 
    /* ----END---- */
    
    std::cout << "Solving ended!" << std::endl;
    std::cout << "--------------" << std::endl;
//    std::for_each(Q.kappa_.begin(), Q.kappa_.end(), [&](double arg){
//	file_tmp << arg << std::endl; });
}
  return 0;
}
