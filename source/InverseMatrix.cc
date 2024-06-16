#include "InverseMatrix.h"

namespace project{
    using namespace dealii;

    template <class MatrixType, class PreconditionerType>
    InverseMatrix<MatrixType, PreconditionerType>::InverseMatrix(
            const MatrixType &        m,
            const PreconditionerType &preconditioner)
            : matrix(&m)
            , preconditioner(&preconditioner)
    {}


    template<class MatrixType, class PreconditionerType>
    void InverseMatrix<MatrixType, PreconditionerType>::vmult(
            Vector<double> &dst,
            const Vector<double> &src) const {
        SolverControl solver_control(src.size(), 1e-6 * src.l2_norm());
        SolverCG<Vector<double>> cg(solver_control);

        dst = 0;

        cg.solve(*matrix, dst, src, *preconditioner);
    }


    template class InverseMatrix<SparseMatrix<double>, SparseDirectUMFPACK>;
    template class InverseMatrix<SparseMatrix<double>, SparseILU<double>>;



}