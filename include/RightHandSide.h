//
// Created by Florian on 16.06.2024.
//

#ifndef dealii_RIGHTHANDSIDE_H
#define dealii_RIGHTHANDSIDE_H
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


namespace project{
    using namespace dealii;
    template <int dim>
    class RightHandSide : public TensorFunction<1, dim>
    {
    public:
        RightHandSide()
                : TensorFunction<1, dim>()
        {}

        virtual Tensor<1, dim> value(const Point<dim> &p) const override;

        virtual void value_list(const std::vector<Point<dim>> &p,
                                std::vector<Tensor<1, dim>> &value) const override;
    };

}


#endif //LAPLACIAN_RIGHTHANDSIDE_H
