#ifndef GALERKIN_MATRICES_HPP
#define GALERKIN_MATRICES_HPP

#include <lf/assemble/assemble.h>
#include <lf/geometry/geometry.h>
#include <lf/mesh/mesh.h>
#include <lf/mesh/test_utils/test_meshes.h>
#include <lf/quad/quad.h> 
#include <lf/uscalfe/uscalfe.h> 

#include <iostream> 
#include <stdexcept> 

#include "params.h"
/**
 * @function dropMatrixRowColumns
 *
 * @brief Enforcing 0 Dirichlet boundary conditions on the galerkin matrices
          by setting the rows/columns corresponding to a vertex on the boundary
          to 0 and the diagonal entry to 1.
 *
 * @tparam SCALAR the type of entries of the matrix
 * @tparam SELECTOR the type of selector identifying the boundary indices
 * @param selectvals the selector that identifies the boundary vertices
 * @param A the galerkin matrix to be modified in triplet format
 *
 */
template <typename SCALAR, typename SELECTOR>
void dropMatrixRowsColumns(SELECTOR &&selectvals,
                           lf::assemble::COOMatrix<SCALAR> &A) {
  const lf::assemble::size_type N(A.cols());
  LF_ASSERT_MSG(A.rows() == N, "Matrix must be square!");
  A.setZero(
      [&selectvals](lf::assemble::gdof_idx_t i, lf::assemble::gdof_idx_t j) {
        return (selectvals(i) || selectvals(j));
      });
  for (lf::assemble::gdof_idx_t dofnum = 0; dofnum < N; ++dofnum) {
    const auto selval{selectvals(dofnum)};
    if (selval) {
      A.AddToEntry(dofnum, dofnum, 1.0);
    }
  }
} /* Copyright ETH Zurich */


namespace ParametricElementMatrices{ 
  
  using elemMatrix = Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>; 
  
  /** 
   * @base_class ElementMatrixProvider
   *
   * @brief Computing the element matrix corresponding to a cell 
   *
   * @tparam MAP the type of the corresponding coefficient functor
   *
   */
  template<class MAP> 
  class ElementMatrixProvider{ 
    public: 
      ElementMatrixProvider() {}
      bool isActive(const lf::mesh::Entity& ){ return true; }
      virtual const elemMatrix Eval(const lf::mesh::Entity&) = 0; 
     virtual ~ElementMatrixProvider(){}
    protected:
	lf::quad::QuadRule myQuadRule_Tria_;
      	lf::quad::QuadRule myQuadRule_Quad_;
	MAP coeff_mapping_; 
        std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>, 5> fe_precomp_; 
  };

  /** 
   * @derived_class alphaElementMatrixProvider
   *
   * @brief Computing the stiffness element matrix corresponding
            to the laplace operator, for a given mesh cell
   *
   * @tparam MAP the type of the corresponding coefficient functor alpha(x)
   *
   */
  template<class MAP> 
  class alphaElementMatrixProvider: public ElementMatrixProvider<MAP>{
    public:
      alphaElementMatrixProvider(
	  std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>,
	  		MAP);
      const elemMatrix Eval(const lf::mesh::Entity&) override; 
      ~alphaElementMatrixProvider() override {}
 };

  /** 
   * @derived_class betaElementMatrixProvider
   *
   * @brief Computing the cross element matrix
            for a given mesh cell
   *
   * @tparam MAP the type of the corresponding coefficient functor beta(x)
   *
   */
  template<class MAP> 
  class betaElementMatrixProvider: public ElementMatrixProvider<MAP>{
    public: 
      betaElementMatrixProvider(
	  std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>,
	  		MAP);
      const elemMatrix Eval(const lf::mesh::Entity&) override; 
      ~betaElementMatrixProvider() override {}
 };

  /** 
   * @derived_class gammaElementMatrixProvider
   *
   * @brief Computing the mass element matrix
            for a given mesh cell
   *
   * @tparam MAP the type of the corresponding coefficient functor gamma(x)
   *
   */
  template<class MAP> 
  class gammaElementMatrixProvider: public ElementMatrixProvider<MAP>{
    public: 
      gammaElementMatrixProvider(
	  std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>,
	  		MAP);
      const elemMatrix Eval(const lf::mesh::Entity&) override; 
      ~gammaElementMatrixProvider() override {}
 };

  /** 
   * @class elementVectorProvider
   *
   * @brief Computing the element vector
            for a given mesh cell
   *
   */
  class elementVectorProvider{
    public:
      elementVectorProvider(
	  std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>>,
	  		std::function<double(Eigen::Vector2d)>);
      const Eigen::VectorXd Eval(const lf::mesh::Entity&);
      bool isActive(const lf::mesh::Entity& ){ return true; }
      ~elementVectorProvider() {}
    private:
	lf::quad::QuadRule myQuadRule_Tria_;
      	lf::quad::QuadRule myQuadRule_Quad_;
	std::function<double(Eigen::Vector2d)> coeff_mapping_; 
        std::array<lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>, 5> fe_precomp_; 
 
  };
      
  /* --------------------Implementation of classes----------------------*/


  /* ----------------------Stiffness matrix---------------------------- */
 
  /**
   * @member alphaElementMatrixProvider
   *
   * @brief Constructing the quadrature rules for the local integral
            computations. Local edge midpoint rule for triangles and 
            quadrilaterals 
            Precomputing reference shape functions and their gradients
            at the quadrature points
   *
   * @tparam MAP the type of the functor coefficient
   * @param fe_space a pointer to the finite element space
   * @param alpha the corresponding coefficient
   *
   */
  template<class MAP> 
  alphaElementMatrixProvider<MAP>::alphaElementMatrixProvider(
      std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space
      ,MAP alpha){
    
    this->coeff_mapping_ = alpha;
    this->myQuadRule_Tria_ = lf::quad::make_TriaQR_EdgeMidpointRule(); 
    this->myQuadRule_Quad_ = lf::quad::make_QuadQR_EdgeMidpointRule(); 
    /* Construction of quadrature rule for local precomputations */
    /* local composite midpoint rule */
    lf::base::RefEl ref_el = lf::base::RefEl::kQuad(); 
    auto fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, lf::quad::make_QuadQR_EdgeMidpointRule()); 
    ref_el = lf::base::RefEl::kTria();
    fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, lf::quad::make_TriaQR_EdgeMidpointRule());
  }

    /**
     * @member Eval 
     *
     * @brief Calculating the element matrix for a given cell 
     *        using the predefined local quadrature rules.
     *
     * @tparam MAP the type of the functor coefficient
     * @param cell the cell to calculate the element matrix
     * @return the element matrix corresponding to 'cell'
     *
     */ 
    template<class MAP>
    const elemMatrix alphaElementMatrixProvider<MAP>::Eval(const lf::mesh::Entity& cell){
   
    lf::geometry::Geometry* geom = cell.Geometry(); 
    lf::base::RefEl ref_elem = cell.RefEl();
    const lf::base::size_type num_nodes = ref_elem.NumNodes(); 
    lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>& pfe=
	this->fe_precomp_[ref_elem.Id()];      
    Eigen::MatrixXd result; 
    result = Eigen::MatrixXd::Zero(num_nodes, num_nodes); 
    Eigen::MatrixXd midPointsRef(2, num_nodes); 
    Eigen::MatrixXd gradRef = 
			pfe.PrecompGradientsReferenceShapeFunctions(); 
 
    if (num_nodes==3){
      midPointsRef = lf::quad::make_TriaQR_EdgeMidpointRule().Points();
      Eigen::VectorXd weights = lf::quad::make_TriaQR_EdgeMidpointRule().Weights(); 
      auto globalMidPoints = geom->Global(midPointsRef);
      const Eigen::MatrixXd DphiInvT = 
		geom->JacobianInverseGramian(midPointsRef);
      const Eigen::VectorXd detPhi = geom->IntegrationElement(midPointsRef);

      /* loop over quadrature nodes */
      for (int n=0; n<num_nodes; n++){ 
	Eigen::MatrixXd glob_grad(2,num_nodes);
	glob_grad = DphiInvT.block<2,2>(0,2*n) *
	  		gradRef.block(0,2*n,num_nodes,2).transpose();
	Eigen::Matrix2d alphaAtx = this->coeff_mapping_(globalMidPoints.col(n)); 
	result += weights(n) * detPhi(n) * glob_grad.transpose() * 
	      alphaAtx * glob_grad; 
      }
    } else if (num_nodes==4){
      midPointsRef = lf::quad::make_QuadQR_EdgeMidpointRule().Points();  
      Eigen::VectorXd weights = lf::quad::make_QuadQR_EdgeMidpointRule().Weights(); 
      auto globalMidPoints = geom->Global(midPointsRef); 
      const Eigen::MatrixXd DphiInvT = 
	geom->JacobianInverseGramian(midPointsRef); 
      const Eigen::VectorXd detPhi = geom->IntegrationElement(midPointsRef); 
      
      /* loop over quadrature nodes */
      for (int n=0; n<num_nodes; n++){ 
	Eigen::MatrixXd glob_grad(2,num_nodes);
	glob_grad = DphiInvT.block<2,2>(0,2*n) *
	  	gradRef.block(0,2*n,num_nodes,2).transpose(); 
	Eigen::Matrix2d alphaAtx = this->coeff_mapping_(globalMidPoints.col(n)); 
	result += weights(n) * detPhi(n) * glob_grad.transpose() * 
	      alphaAtx * glob_grad;
      }
    }
    return result; 
  }

  /* ----------------------Cross matrix---------------------------- */
 
  /**
   * @member betaElementMatrixProvider
   *
   * @brief Constructing the quadrature rules for the local integral
            computations. Composite local trapezoidal rule for triangles and 
            quadrilaterals 
            Precomputing reference shape functions and their gradients
            at the quadrature points
   *
   * @tparam MAP the type of the functor coefficient
   * @param fe_space a pointer to the finite element space
   * @param beta the corresponding coefficient
   *
   */
  template<class MAP> 
  betaElementMatrixProvider<MAP>::betaElementMatrixProvider(
      std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space
      ,MAP beta){
 
    this->coeff_mapping_ = beta;
    lf::base::RefEl ref_el = lf::base::RefEl::kQuad(); 
    Eigen::MatrixXd QP_quad(2,4);
    QP_quad << 0.0, 1.0, 1.0, 0.0,
      	       0.0, 0.0, 1.0, 1.0; 
    Eigen::VectorXd QW_quad(4);
    QW_quad << 0.25, 0.25, 0.25, 0.25; 
    lf::quad::QuadRule myQuadRule_kQuad(ref_el, QP_quad, QW_quad, 3);
    this->myQuadRule_Quad_ = myQuadRule_kQuad; 
    auto fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, myQuadRule_kQuad); 
    
    ref_el = lf::base::RefEl::kTria();
    Eigen::MatrixXd QP_tria(2,3);
    QP_tria << 0.0, 1.0, 0.0, 
      	       0.0, 0.0, 1.0; 
    Eigen::VectorXd QW_tria(3);
    QW_tria << 1./6., 1./6., 1./6.; 
    lf::quad::QuadRule myQuadRule_kTria(ref_el, QP_tria, QW_tria, 3);
    this->myQuadRule_Tria_ = myQuadRule_kTria; 
    fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, myQuadRule_kTria);
  }

    /**
     * @member Eval 
     *
     * @brief Calculating the element matrix for a given cell 
     *        using the predefined local quadrature rules.
     *
     * @tparam MAP the type of the functor coefficient
     * @param cell the cell to calculate the element matrix
     * @return the element matrix corresponding to 'cell'
     *
     */ 
    template<class MAP>
    const elemMatrix betaElementMatrixProvider<MAP>::Eval(const lf::mesh::Entity& cell){
   
    lf::geometry::Geometry* geom = cell.Geometry(); 
    lf::base::RefEl ref_elem = cell.RefEl();
    const lf::base::size_type num_nodes = ref_elem.NumNodes(); 
    lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>& pfe=
	this->fe_precomp_[ref_elem.Id()];      
    Eigen::MatrixXd result; 
    result = Eigen::MatrixXd::Zero(num_nodes, num_nodes); 
    Eigen::MatrixXd PointsRef(2, num_nodes); 
    Eigen::MatrixXd shapeRef = 
		pfe.PrecompReferenceShapeFunctions();
    Eigen::MatrixXd gradRef = 
		pfe.PrecompGradientsReferenceShapeFunctions(); 
 
    if (num_nodes==3){
      PointsRef = this->myQuadRule_Tria_.Points();
      Eigen::VectorXd weights = this->myQuadRule_Tria_.Weights(); 
      auto globalPoints = geom->Global(PointsRef);
      const Eigen::VectorXd detPhi = geom->IntegrationElement(PointsRef);
      const Eigen::MatrixXd DphiInvT = 
		geom->JacobianInverseGramian(PointsRef);
      /* loop over quadrature nodes */
      for (int n=0; n<num_nodes; n++){ 	
	Eigen::MatrixXd glob_grad(2, num_nodes);
	glob_grad = DphiInvT.block<2,2>(0,2*n) *
	  		gradRef.block(0,2*n,num_nodes,2).transpose();
	Eigen::Vector2d betaAtx = this->coeff_mapping_(globalPoints.col(n));
	result += weights(n) * detPhi(n) * 
	       shapeRef.col(n) * betaAtx.transpose() * glob_grad; 
      }

    } else if (num_nodes==4){
      PointsRef = this->myQuadRule_Quad_.Points();  
      Eigen::VectorXd weights = this->myQuadRule_Quad_.Weights(); 
      auto globalPoints = geom->Global(PointsRef); 
      const Eigen::VectorXd detPhi = geom->IntegrationElement(PointsRef); 
      /* loop over quadrature nodes */
      const Eigen::MatrixXd DphiInvT = 
		geom->JacobianInverseGramian(PointsRef);
 	for (int n=0; n<num_nodes; n++){
	Eigen::MatrixXd glob_grad(2,num_nodes);
	glob_grad = DphiInvT.block<2,2>(0,2*n) *
	  		gradRef.block(0,2*n,num_nodes,2).transpose();
	Eigen::Vector2d betaAtx = this->coeff_mapping_(globalPoints.col(n)); 	
	result += weights(n) * detPhi(n) *
	  	shapeRef.col(n) * betaAtx.transpose() * glob_grad; 
      }
    }
    return result; 
  }


  /* ----------------------Mass matrix---------------------------- */

 /**
   * @member gammaElementMatrixProvider
   *
   * @brief Constructing the quadrature rules for the local integral
            computations. Composite local trapezoidal rule for triangles and 
            quadrilaterals 
            Precomputing reference shape functions and their gradients
            at the quadrature points
   *
   * @tparam MAP the type of the functor coefficient
   * @param fe_space a pointer to the finite element space
   * @param gamma the corresponding coefficient
   *
   */
  template<class MAP> 
  gammaElementMatrixProvider<MAP>::gammaElementMatrixProvider(
      std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space
      ,MAP gamma){
    
    this->coeff_mapping_ = gamma;
    lf::base::RefEl ref_el = lf::base::RefEl::kQuad(); 
    Eigen::MatrixXd QP_quad(2,4);
    QP_quad << 0.0, 1.0, 1.0, 0.0,
      	       0.0, 0.0, 1.0, 1.0; 
    Eigen::VectorXd QW_quad(4);
    QW_quad << 0.25, 0.25, 0.25, 0.25; 
    lf::quad::QuadRule myQuadRule_kQuad(ref_el, QP_quad, QW_quad, 3);
    this->myQuadRule_Quad_ = myQuadRule_kQuad; 
    auto fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, myQuadRule_kQuad); 
    
    ref_el = lf::base::RefEl::kTria();
    Eigen::MatrixXd QP_tria(2,3);
    QP_tria << 0.0, 1.0, 0.0, 
      	       0.0, 0.0, 1.0; 
    Eigen::VectorXd QW_tria(3);
    QW_tria << 1./6., 1./6., 1./6.; 
    lf::quad::QuadRule myQuadRule_kTria(ref_el, QP_tria, QW_tria, 3);
    this->myQuadRule_Tria_ = myQuadRule_kTria; 
    fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, myQuadRule_kTria);
  }

    /**
     * @member Eval 
     *
     * @brief Calculating the element matrix for a given cell 
     *        using the predefined local quadrature rules.
     *
     * @tparam MAP the type of the functor coefficient
     * @param cell the cell to calculate the element matrix
     * @return the element matrix corresponding to 'cell'
     *
     */ 
    template<class MAP>
    const elemMatrix gammaElementMatrixProvider<MAP>::Eval(const lf::mesh::Entity& cell){
   
    lf::geometry::Geometry* geom = cell.Geometry(); 
    lf::base::RefEl ref_elem = cell.RefEl();
    const lf::base::size_type num_nodes = ref_elem.NumNodes(); 
    lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>& pfe=
	this->fe_precomp_[ref_elem.Id()];      
    Eigen::MatrixXd result; 
    result = Eigen::MatrixXd::Zero(num_nodes, num_nodes); 
    Eigen::MatrixXd PointsRef(2, num_nodes); 
    Eigen::MatrixXd shapeRef = 
			pfe.PrecompReferenceShapeFunctions();
      
    if (num_nodes==3){
      PointsRef = this->myQuadRule_Tria_.Points();
      Eigen::VectorXd weights = this->myQuadRule_Tria_.Weights(); 
      auto globalPoints = geom->Global(PointsRef);
      const Eigen::VectorXd detPhi = geom->IntegrationElement(PointsRef);
      /* loop over quadrature nodes */
      for (int n=0; n<num_nodes; n++){ 	
	double gammaAtx = this->coeff_mapping_(globalPoints.col(n));
	result += weights(n) * detPhi(n) * gammaAtx * 
	  shapeRef.col(n) * shapeRef.col(n).transpose(); 
      }

    } else if (num_nodes==4){
      PointsRef = this->myQuadRule_Quad_.Points();  
      Eigen::VectorXd weights = this->myQuadRule_Quad_.Weights(); 
      auto globalPoints = geom->Global(PointsRef); 
      const Eigen::VectorXd detPhi = geom->IntegrationElement(PointsRef); 
      /* loop over quadrature nodes */
      for (int n=0; n<num_nodes; n++){
	double gammaAtx = this->coeff_mapping_(globalPoints.col(n)); 	
	result += weights(n) * detPhi(n) * gammaAtx * 
	  shapeRef.col(n) * shapeRef.col(n).transpose(); 
      }
    }
    return result; 
  }
    
    
  /* ----------------------Load Vector---------------------------- */

 /**
   * @member elementVectorProvider
   *
   * @brief Constructing the quadrature rules for the local integral
            computations. Composite local trapezoidal rule for triangles and 
            quadrilaterals 
            Precomputing reference shape functions and their gradients
            at the quadrature points
   *
   * @param fe_space a pointer to the finite element space
   * @param f the corresponding load functor
   *
   */
  elementVectorProvider::elementVectorProvider(
      std::shared_ptr<const lf::uscalfe::UniformScalarFESpace<double>> fe_space
      ,std::function<double(Eigen::Vector2d)> f){
    
    this->coeff_mapping_ = f;
    lf::base::RefEl ref_el = lf::base::RefEl::kQuad(); 
    Eigen::MatrixXd QP_quad(2,4);
    QP_quad << 0.0, 1.0, 1.0, 0.0,
      	       0.0, 0.0, 1.0, 1.0; 
    Eigen::VectorXd QW_quad(4);
    QW_quad << 0.25, 0.25, 0.25, 0.25; 
    lf::quad::QuadRule myQuadRule_kQuad(ref_el, QP_quad, QW_quad, 3);
    this->myQuadRule_Quad_ = myQuadRule_kQuad; 
    auto fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, myQuadRule_kQuad); 
    
    ref_el = lf::base::RefEl::kTria();
    Eigen::MatrixXd QP_tria(2,3);
    QP_tria << 0.0, 1.0, 0.0, 
      	       0.0, 0.0, 1.0; 
    Eigen::VectorXd QW_tria(3);
    QW_tria << 1./6., 1./6., 1./6.; 
    lf::quad::QuadRule myQuadRule_kTria(ref_el, QP_tria, QW_tria, 3);
    this->myQuadRule_Tria_ = myQuadRule_kTria; 
    fe = fe_space->ShapeFunctionLayout(ref_el); 
    if (fe!=nullptr)
      this->fe_precomp_[ref_el.Id()] = 
	  lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>(
		  fe, myQuadRule_kTria);
  }


    /**
     * @member Eval 
     *
     * @brief Calculating the element vector for a given cell 
     *        using the predefined local quadrature rules.
     *
     * @param cell the cell to calculate the element matrix
     * @return the element vector corresponding to 'cell'
     *
     */ 
    const Eigen::VectorXd elementVectorProvider::Eval(const lf::mesh::Entity& cell){
   
    lf::geometry::Geometry* geom = cell.Geometry(); 
    lf::base::RefEl ref_elem = cell.RefEl();
    const lf::base::size_type num_nodes = ref_elem.NumNodes(); 
    lf::uscalfe::PrecomputedScalarReferenceFiniteElement<double>& pfe=
	this->fe_precomp_[ref_elem.Id()];      
    Eigen::VectorXd result; 
    result = Eigen::VectorXd::Zero(num_nodes);
    Eigen::MatrixXd PointsRef(2, num_nodes); 
    Eigen::MatrixXd shapeRef = 
			pfe.PrecompReferenceShapeFunctions();
    if (num_nodes==3){
      PointsRef = this->myQuadRule_Tria_.Points();
      Eigen::VectorXd weights = this->myQuadRule_Tria_.Weights(); 
      auto globalPoints = geom->Global(PointsRef);
      const Eigen::VectorXd detPhi = geom->IntegrationElement(PointsRef);
      /* loop over quadrature nodes */
      for (int n=0; n<num_nodes; n++){ 	
	double fAtx = this->coeff_mapping_(globalPoints.col(n));
	for (int i=0; i<num_nodes; i++){ 
	    result(i) += weights(n) * detPhi(n) * shapeRef(i,n) * 
	      fAtx;	
	} 
      }

    } else if (num_nodes==4){
      PointsRef = this->myQuadRule_Quad_.Points();  
      Eigen::VectorXd weights = this->myQuadRule_Quad_.Weights(); 
      auto globalPoints = geom->Global(PointsRef); 
      const Eigen::VectorXd detPhi = geom->IntegrationElement(PointsRef); 
      /* loop over quadrature nodes */
      for (int n=0; n<num_nodes; n++){
	double fAtx = this->coeff_mapping_(globalPoints.col(n)); 	
	for (int i=0; i<num_nodes; i++){ 
	    result(i) += weights(n) * detPhi(n) * shapeRef(i,n) * 
	      fAtx;
	}
      }
    }
    return result; 
  }


} // namespace ParametricElementMatrices

using SparseMatrix = Eigen::SparseMatrix<double>; 
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
    template<typename T> 
    SincMatrix(T&&); /// Move and copy constructors 
    template<typename T> 
    SincMatrix(double, T&& L, T&& M); 
    Eigen::VectorXd operator* (Eigen::VectorXd& f) const; 
  private:
    double k_; 
    SparseMatrix L_, M_;
    Eigen::VectorXd K_; 
    std::vector<double> K_STL_;  
    template<typename SincMatType, bool is_theta_zero>
    friend class thetaSchemeTimeStepper; 
};

template<typename T> 
SincMatrix::SincMatrix(T&& other): 
	k_(other.k_), L_(std::forward<T>(other.L_)),
	    M_(std::forward<T>(other.M_)),
	  K_(std::forward<T>(other.K_)),
	  K_STL_(std::forward<T>(other.K_STL_)) {}


template<typename T> 
SincMatrix::SincMatrix(double h, T&& L, 
    T&& M):
      M_(std::forward<T>(M)), 
	L_(std::forward<T>(L)) {
  	k_ = -1/(std::log(h));
	K_ = Eigen::VectorXd::LinSpaced(
	    ceil(M_PI*M_PI / (4.*beta*k_*k_)) + ceil(M_PI*M_PI / (4.*(1-beta)*k_*k_)) + 1, 
	    -ceil(M_PI*M_PI / (4.*beta*k_*k_)),  ceil(M_PI*M_PI / (4.*(1-beta)*k_*k_)) );
	K_STL_.resize(K_.size());
	std::cout << K_.size() << std::endl;
	for (int i=0; i<K_.size(); i++)
	  K_STL_.at(i) = K_(i); 
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
  Eigen::SparseLU<Eigen::SparseMatrix<double>> solver;
  
  double Cb = M_PI/(3*std::sin(M_PI*beta));

  std::for_each(K_STL_.begin(), K_STL_.end(), [&](double arg){
       solver.compute(M_+std::exp(2*arg*k_)*L_);
       res = res + (k_/Cb)*std::exp(2*beta*arg*k_)*solver.solve(f);});
 
  return res; 
}

#endif 
