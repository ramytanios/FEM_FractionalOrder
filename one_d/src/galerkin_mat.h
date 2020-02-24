#ifndef GALERKIN_MAT_HPP
#define GALERKIN_MAT_HPP

#include <memory>
#include <Eigen/Sparse>

using SP_MAT = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;
using TripletForm = std::vector<Triplet>; 

class Galerkin_Mat{
  public: 
    Galerkin_Mat(size_t s, std::shared_ptr<Eigen::VectorXd> mesh_p); 
    
    void setup_galerkin_mat(); 
   
    void impose_dirichlet();

    std::shared_ptr<SP_MAT> M();
   
    std::shared_ptr<SP_MAT> L();
    
    std::shared_ptr<SP_MAT> A();
    
    std::shared_ptr<SP_MAT> Ma();
    
    ~Galerkin_Mat() = default; 
  private:
    SP_MAT M_; 
    SP_MAT L_;
    SP_MAT A_; 
    SP_MAT Ma_; 
    size_t s_;
    std::shared_ptr<Eigen::VectorXd> mesh_p_; 
};

#endif 
