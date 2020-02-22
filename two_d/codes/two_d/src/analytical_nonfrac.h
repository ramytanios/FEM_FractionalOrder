#ifndef ANALYTICAL_NON_FRACTIONAL
#define ANALYTICAL_NON_FRACTIONAL

#include <Eigen/Eigen> 

template<typename mesh_ptr_t> 
Eigen::VectorXd exact(double T, double initial_cond_const, mesh_ptr_t mesh_p){
  int K = 1000; 
  Eigen::VectorXd sol = Eigen::VectorXd::Zero(mesh_p->NumEntities(2)); 
  for (const lf::mesh::Entity &ent: mesh_p->Entities(2)){
    const lf::geometry::Geometry* geo_p = ent.Geometry(); 
    Eigen::MatrixXd node = lf::geometry::Corners(*geo_p);
    double x = node(0,0); double y = node(1,0); 
    double eval = 0; 
    for (int m=1; m<K; m++){
      for (int n=1; n<K; n++){
	double lambda = std::pow(M_PI,2)*(std::pow(m,2) + std::pow(n,2));
	eval += 4 * initial_cond_const / (M_PI*M_PI*m*n) * 
	  (std::cos(m*M_PI)-1) * (std::cos(n*M_PI)-1) * 
	  std::sin(m*M_PI*x) * std::sin(n*M_PI*y) * 
	  std::exp(-lambda*T); 
      }
    }
    sol(mesh_p->Index(ent)) = eval; 
  }
  return sol; 
}

#endif 
