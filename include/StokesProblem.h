//
// Created by Florian on 16.06.2024.
//

#ifndef dealii_STOKES_H
#define dealii_STOKES_H

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

namespace project {
    using namespace dealii;


    template<int dim>
    struct InnerPreconditioner;

    template<>
    struct InnerPreconditioner<2> {
        using type = SparseDirectUMFPACK;
    };

    template<>
    struct InnerPreconditioner<3> {
        using type = SparseILU<double>;
    };


    template<int dim>
    class StokesProblem {
    public:
        StokesProblem(const unsigned int degree);

        void run();

    private:
        void setup_dofs();

        void assemble_system();

        void solve();

        void output_results(const unsigned int refinement_cycle) const;

        void refine_mesh();

        const unsigned int degree;

        Triangulation<dim> triangulation;
        FESystem<dim> fe;
        DoFHandler<dim> dof_handler;

        AffineConstraints<double> constraints;

        BlockSparsityPattern sparsity_pattern;
        BlockSparseMatrix<double> system_matrix;

        BlockSparsityPattern preconditioner_sparsity_pattern;
        BlockSparseMatrix<double> preconditioner_matrix;

        BlockVector<double> solution;
        BlockVector<double> system_rhs;

        std::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
    };
}

#endif //dealii_STOKES_H
