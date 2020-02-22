#include <Eigen/Core>
#include <Eigen/Dense>
#include "Q.h"
#include <memory>
#include <type_traits>
class MatrixReplacement;
namespace Eigen {
    namespace internal {
        // MatrixReplacement looks-like a SparseMatrix, so let's inherits its traits:
        template<>
        struct traits<MatrixReplacement> :  public Eigen::internal::traits<Eigen::SparseMatrix<double> >
        {};
    }
}

class MatrixReplacement : public Eigen::EigenBase<MatrixReplacement> {
public:
    // Required typedefs, constants, and method:
    typedef double Scalar;
    typedef double RealScalar;
    typedef int StorageIndex;
    enum {
        ColsAtCompileTime = Eigen::Dynamic,
        MaxColsAtCompileTime = Eigen::Dynamic,
        IsRowMajor = false
    };
    Index rows() const { return M_->rows(); }
    Index cols() const { return M_->cols(); }
    template<typename Rhs>
    Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct> operator*(const Eigen::MatrixBase<Rhs>& x) const {
        return Eigen::Product<MatrixReplacement,Rhs,Eigen::AliasFreeProduct>(*this, x.derived());
    }
    
    MatrixReplacement(const Eigen::SparseMatrix<double>& M, const Eigen::SparseMatrix<double>& A,
                        const Eigen::SparseMatrix<double>& Ma, const Eigen::SparseMatrix<double>& L,
                        const SincMatrix& Q, const double& dt) {
      	
        M_ = std::make_shared<Eigen::SparseMatrix<double>>(M);
	A_ = std::make_shared<Eigen::SparseMatrix<double>>(A);
	L_ = std::make_shared<Eigen::SparseMatrix<double>>(L);
	Ma_ = std::make_shared<Eigen::SparseMatrix<double>>(Ma);
	Q_ = std::make_shared<SincMatrix>(Q);
	dt_= std::make_shared<double>(dt);
    }
    
     std::shared_ptr< Eigen::SparseMatrix<double> > M_;
     std::shared_ptr< Eigen::SparseMatrix<double> > Ma_;
     std::shared_ptr< Eigen::SparseMatrix<double> > L_;
     std::shared_ptr< Eigen::SparseMatrix<double> > A_;
     std::shared_ptr<SincMatrix> Q_;
     std::shared_ptr<double> dt_;
};
// Implementation of MatrixReplacement * Eigen::DenseVector though a specialization of internal::generic_product_impl:
namespace Eigen {
    namespace internal {
        template<typename Rhs>
        struct generic_product_impl<MatrixReplacement, Rhs, SparseShape, DenseShape, GemvProduct> // GEMV stands for matrix-vector
        : generic_product_impl_base<MatrixReplacement,Rhs,generic_product_impl<MatrixReplacement,Rhs> >
        {
            typedef typename Product<MatrixReplacement,Rhs>::Scalar Scalar;
            template<typename Dest>
            static void scaleAndAddTo(Dest& dst, const MatrixReplacement& lhs, const Rhs& rhs, const Scalar& alpha)
            {
                // This method should implement "dst += alpha * lhs * rhs" inplace,
                // however, for iterative solvers, alpha is always equal to 1, so let's not bother about it.
                assert(alpha==Scalar(1) && "scaling is not implemented");
                EIGEN_ONLY_USED_FOR_DEBUG(alpha);
                // Here we could simply call dst.noalias() += lhs.my_matrix() * rhs,
                // but let's do something fancier (and less efficient):
                Eigen::VectorXd tmpp = *lhs.L_ * rhs;
                dst = *lhs.M_*rhs + (theta**lhs.dt_) * (*lhs.A_*rhs +
                                                              *lhs.M_*(*lhs.Q_*tmpp));
            }
        };
    }
}
