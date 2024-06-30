//
// Created by Florian on 16.06.2024.
//

#ifndef dealii_STOKES_H
#define dealii_STOKES_H

#include <deal.II/base/quadrature_lib.h>
#include <deal.II/base/function.h>
#include <deal.II/base/utilities.h>
#include <deal.II/base/tensor_function.h>
#include <deal.II/base/tensor.h>
 
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
 #include <deal.II/fe/fe_nothing.h>
 #include <deal.II/fe/fe_interface_values.h>
 #include <deal.II/fe/fe_update_flags.h>
 
#include <deal.II/numerics/vector_tools.h>
#include <deal.II/numerics/data_out.h>
#include <deal.II/numerics/error_estimator.h>
 
#include <deal.II/lac/sparse_direct.h>
 
#include <deal.II/lac/sparse_ilu.h>

#include <deal.II/non_matching/mesh_classifier.h>
#include <deal.II/non_matching/fe_immersed_values.h>
#include <deal.II/non_matching/fe_values.h>

#include <deal.II/base/function_signed_distance.h>

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

   template <int dim>
  class StokesProblem
  {
  public:
    StokesProblem(const unsigned int degree);
    void run();
 
  private:
    void setup_dofs();
    void assemble_system();
    void solve();
    void output_results(const unsigned int refinement_cycle) const;
    void refine_mesh();
    void setup_discrete_level_set();
    

    const unsigned int degree;
 
    Triangulation<dim>  triangulation;
    
    
    //2
    bool face_has_ghost_penalty(
      const typename Triangulation<dim>::active_cell_iterator &cell,
      const unsigned int face_index) const;
       
    hp::FECollection<dim> fe_collection;

    const FE_Q<dim> fe_level_set;
    DoFHandler<dim> level_set_dof_handler;
    Vector<double>  level_set;
    
    NonMatching::MeshClassifier<dim> mesh_classifier;

    
    //2end

    //removed const FESystem<dim> fe; 
    DoFHandler<dim>     dof_handler;
 
    AffineConstraints<double> constraints;
 
    BlockSparsityPattern      sparsity_pattern;
    BlockSparseMatrix<double> system_matrix;
 
    BlockSparsityPattern      preconditioner_sparsity_pattern;
    BlockSparseMatrix<double> preconditioner_matrix;
 
    BlockVector<double> solution;
    BlockVector<double> system_rhs;
 
    std::shared_ptr<typename InnerPreconditioner<dim>::type> A_preconditioner;
  };
}

#endif //dealii_STOKES_H
