//
// Created by Florian on 16.06.2024.
//

#ifndef dealii_SCHURCOMPLEMENT_H
#define dealii_SCHURCOMPLEMENT_H
#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/logstream.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor_function.h>

#include <deal.II/lac/block_vector.h>
#include <deal.II/lac/full_matrix.h>
#include <deal.II/lac/block_sparse_matrix.h>
#include <deal.II/lac/solver_cg.h>
#include <deal.II/lac/precondition.h>
#include <deal.II/lac/affine_constraints.h>

#include <deal.II/grid/tria.h>
#include <deal.II/grid/grid_generator.h>
#include <deal.II/grid/grid_tools.h>
#include <deal.II/grid/grid_refinement.h>

#include <deal.II/dofs/dof_handler.h>
#include <deal.II/dofs/dof_renumbering.h>
#include <deal.II/dofs/dof_tools.h>


#include <deal.II/fe/fe_q.h>
#include <deal.II/fe/fe_system.h>
#include <deal.II/fe/fe_values.h>

#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/matrix_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>

#include <deal.II/lac/sparse_direct.h>

#include <deal.II/lac/sparse_ilu.h>

#include <iostream>
#include <fstream>
#include <memory>
#include "InverseMatrix.h"



namespace project{
    using namespace dealii;

    template <class PreconditionerType>
    class SchurComplement : public Subscriptor
    {
    public:
        SchurComplement(
                const BlockSparseMatrix<double> &system_matrix,
                const InverseMatrix<SparseMatrix<double>, PreconditionerType> &A_inverse);

        void vmult(Vector<double> &dst, const Vector<double> &src) const;

    private:
        const SmartPointer<const BlockSparseMatrix<double>> system_matrix;
        const SmartPointer<
                const InverseMatrix<SparseMatrix<double>, PreconditionerType>>
                A_inverse;

        mutable Vector<double> tmp1, tmp2;
    };

}



#endif //LAPLACIAN_SCHURCOMPLEMENT_H
