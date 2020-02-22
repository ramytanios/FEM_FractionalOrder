#include "space_discr.h"
#include "time_stepper.h"
#include "my_io.h"
#include "analytical_nonfrac.h"
//#include "routines.h"
#include "pde_weights.h"
#include "params.h"
#include <lf/base/base.h>
#include <lf/io/io.h>
#include <lf/mesh/hybrid2d/hybrid2d.h>
#include <boost/filesystem.hpp>

int main(){
    int N0 = 4; 
//  for (int l=0; l<5; l++){
  /* Discrete finite element space */
  std::cout << "Setting up mesh and FEM space. " << std::endl;
 // boost::filesystem::path here = __FILE__;
 // auto mesh_path = here.parent_path().parent_path() / "mesh/square.msh";
 // auto mesh_factory = std::make_unique<lf::mesh::hybrid2d::MeshFactory>(2);
 // const lf::io::GmshReader reader(std::move(mesh_factory), mesh_path.string());
 // auto mesh_p = reader.mesh();
  std::shared_ptr<lf::mesh::hybrid2d::MeshFactory> mesh_factory_ptr =
  std::make_shared<lf::mesh::hybrid2d::MeshFactory>(2);
// Triangular tensor product mesh
  lf::mesh::hybrid2d::TPQuadMeshBuilder builder(mesh_factory_ptr);
// Set mesh parameters following the Builder pattern
// Domain is the unit square
  int NoXCells = N0 * std::pow(2,3); 
  builder.setBottomLeftCorner(Eigen::Vector2d{0, 0})
    .setTopRightCorner(Eigen::Vector2d{1, 1})
    .setNoXCells(NoXCells)
    .setNoYCells(NoXCells);
  std::shared_ptr<lf::mesh::Mesh> mesh_p{builder.Build()};

  std::shared_ptr<lf::uscalfe::FeSpaceLagrangeO1<double>> fe_space =
    	std::make_shared<lf::uscalfe::FeSpaceLagrangeO1<double>>(mesh_p); 
  const lf::assemble::DofHandler &dofh{fe_space->LocGlobMap()}; 
  const lf::base::size_type N_dofs(dofh.NoDofs()); 
  /* ----END---- */ 

  double T = 0.1; 
  double mesh_size = 1.0/NoXCells * std::sqrt(2);
  double dt = std::pow(mesh_size,2); 
  const int num_steps = T/dt;
  std::cout << "The number of time steps: " << num_steps << std::endl;
  
  /* Galerkin matrices and vectors assembly */
  std::cout << "Assembly of galerkin matrices. " << std::endl;
  
  ParametricElementMatrices::alphaElementMatrixProvider<decltype(PdeCoefficients::alpha_x)> 
    	stiffnessElemMatrixProvider(fe_space, PdeCoefficients::alpha_x);
  ParametricElementMatrices::gammaElementMatrixProvider<decltype(PdeCoefficients::gamma_x)>
    	massElemMatrixProvider(fe_space, [](Eigen::Vector2d x)->double{(void)x; return 1;});
  ParametricElementMatrices::gammaElementMatrixProvider<decltype(PdeCoefficients::gamma_x)>
    	weightedMassElemMatrixProvider(fe_space, PdeCoefficients::gamma_x);
  ParametricElementMatrices::betaElementMatrixProvider<decltype(PdeCoefficients::beta_x)>
    	crossElemMatrixProvider(fe_space, PdeCoefficients::beta_x);
  ParametricElementMatrices::elementVectorProvider
  	loadElemVectorProvider(fe_space, PdeCoefficients::f_x); 
  
  Eigen::VectorXd phi(N_dofs);
  phi.setZero();   
  lf::assemble::COOMatrix<double> L_COO(N_dofs, N_dofs); 
  lf::assemble::COOMatrix<double> A_COO(N_dofs, N_dofs);
  lf::assemble::COOMatrix<double> M_COO(N_dofs, N_dofs);

  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, stiffnessElemMatrixProvider, L_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, weightedMassElemMatrixProvider, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, crossElemMatrixProvider, A_COO);
  lf::assemble::AssembleMatrixLocally(0, dofh, dofh, massElemMatrixProvider, M_COO);
  lf::assemble::AssembleVectorLocally(0, dofh, loadElemVectorProvider, phi); 
  /* ----END---- */ 

  
  /* Initial condition - L2 projected on the discrete FEM space  */
  Eigen::VectorXd u0 = Eigen::VectorXd::Zero(N_dofs);
  lf::uscalfe::MeshFunctionGlobal mf_initial_cond(PdeCoefficients::initial_cond); 
  u0 = lf::uscalfe::NodalProjection(*fe_space, mf_initial_cond); 

  lf::io::VtkWriter vtk_writer0(mesh_p, "payoff.vtk");
  auto nodal_data0 = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2); 
  for (lf::base::size_type i=0; i<N_dofs; i++){ 
    nodal_data0->operator()(dofh.Entity(i)) = u0[i];
  };
  vtk_writer0.WritePointData("payoff.vtk", *nodal_data0);
  std::cout << "Initial condition  written to payoff.vtk ." << std::endl;

  my_io::pyWriter py_writer(mesh_p); 
  py_writer.writePointData("payoff.txt", u0);
  std::cout << "Initial condition  written to payoff.txt ." << std::endl;
  /* ----END---- */ 
 

  /* Enforce dirichlet BC, 0 on the matrices and vectors */
  lf::mesh::utils::CodimMeshDataSet<bool>
    		bd_flags = lf::mesh::utils::flagEntitiesOnBoundary(mesh_p,2); 
  std::function<bool(unsigned int dof_idx)> bd_selector_d=
      [&] (unsigned int dof_idx){
      return bd_flags(dofh.Entity(dof_idx));
  };
  dropMatrixRowsColumns(bd_selector_d ,L_COO);
  dropMatrixRowsColumns(bd_selector_d ,M_COO);
  dropMatrixRowsColumns(bd_selector_d ,A_COO);
  for (lf::base::size_type i=0; i<N_dofs; i++){ // for Dirichlet BC. 
    if (bd_selector_d(i) == 1)
      u0[i] = 0; 
  }
  /* ----END---- */ 

  Eigen::SparseMatrix<double> L = L_COO.makeSparse();
  Eigen::SparseMatrix<double> A = A_COO.makeSparse(); 
  Eigen::SparseMatrix<double> M = M_COO.makeSparse();
  
  SincMatrix Q(mesh_size, std::move(L), std::move(M));
  /* Time stepper and solver */
  thetaSchemeTimeStepper<SincMatrix,testbool> 
        time_stepper(Q, std::move(A), dt, mesh_size, N_dofs);
  Eigen::VectorXd u_h = u0;
  std::cout << "Stepping in time " << std::endl;
  for (int i=0; i<num_steps; i++){
    u_h = time_stepper.discreteEvolutionOperator(u_h); 
  }
  /* ----END---- */ 


  /* Writing results to data files */
  lf::io::VtkWriter vtk_writer(mesh_p, "discrete_solution.vtk");
  auto nodal_data = lf::mesh::utils::make_CodimMeshDataSet<double>(mesh_p, 2); 
  for (lf::base::size_type i=0; i<N_dofs; i++){ 
    nodal_data->operator()(dofh.Entity(i)) = u_h[i];
  }
  vtk_writer.WritePointData("discrete_solution.vtk", *nodal_data);
  std::cout << "Solution written to vtk discrete_solution.vtk ." << std::endl;

  py_writer.writePointData("discrete_solution.txt", u_h);
  std::cout << "Solution written to discrete_solution.txt ." << std::endl;

  return 0;

}
