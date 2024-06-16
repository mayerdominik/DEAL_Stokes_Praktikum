//
// Created by Florian on 16.06.2024.
//
#include "StokesProblem.h"
#include "SchurComplement.h"
#include "RightHandSide.h"
#include "BoundaryValues.h"



namespace project {
    using namespace dealii;


    template<int dim>
    StokesProblem<dim>::StokesProblem(const unsigned int degree)
            : degree(degree),
            triangulation(Triangulation<dim>::maximum_smoothing),
            fe(FE_Q<dim>(degree + 1), dim, FE_Q<dim>(degree), 1),
            dof_handler(triangulation) {}

    template<int dim>
    void StokesProblem<dim>::setup_dofs() {
        A_preconditioner.reset();
        system_matrix.clear();
        preconditioner_matrix.clear();

        dof_handler.distribute_dofs(fe);
        DoFRenumbering::Cuthill_McKee(dof_handler);

        std::vector<unsigned int> block_component(dim + 1, 0);
        block_component[dim] = 1;
        DoFRenumbering::component_wise(dof_handler, block_component);

        {
            constraints.clear();

            const FEValuesExtractors::Vector velocities(0);
            DoFTools::make_hanging_node_constraints(dof_handler, constraints);
            VectorTools::interpolate_boundary_values(dof_handler,
                                                     1,
                                                     BoundaryValues<dim>(),
                                                     constraints,
                                                     fe.component_mask(velocities));
        }

        constraints.close();
        const std::vector<types::global_dof_index> dofs_per_block =
                DoFTools::count_dofs_per_fe_block(dof_handler, block_component);
        const types::global_dof_index n_u = dofs_per_block[0];
        const types::global_dof_index n_p = dofs_per_block[1];

        std::cout << "   Number of active cells: " << triangulation.n_active_cells()
                  << std::endl
                  << "   Number of degrees of freedom: " << dof_handler.n_dofs()
                  << " (" << n_u << '+' << n_p << ')' << std::endl;

        {
            BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);

            Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1);
            for (unsigned int c = 0; c < dim + 1; ++c)
                for (unsigned int d = 0; d < dim + 1; ++d)
                    if (!((c == dim) && (d == dim)))
                        coupling[c][d] = DoFTools::always;
                    else
                        coupling[c][d] = DoFTools::none;

            DoFTools::make_sparsity_pattern(
                    dof_handler, coupling, dsp, constraints, false);

            sparsity_pattern.copy_from(dsp);
        }

        {
            BlockDynamicSparsityPattern preconditioner_dsp(dofs_per_block,
                                                           dofs_per_block);

            Table<2, DoFTools::Coupling> preconditioner_coupling(dim + 1, dim + 1);
            for (unsigned int c = 0; c < dim + 1; ++c)
                for (unsigned int d = 0; d < dim + 1; ++d)
                    if (((c == dim) && (d == dim)))
                        preconditioner_coupling[c][d] = DoFTools::always;
                    else
                        preconditioner_coupling[c][d] = DoFTools::none;

            DoFTools::make_sparsity_pattern(dof_handler,
                                            preconditioner_coupling,
                                            preconditioner_dsp,
                                            constraints,
                                            false);

            preconditioner_sparsity_pattern.copy_from(preconditioner_dsp);
        }
        system_matrix.reinit(sparsity_pattern);
        preconditioner_matrix.reinit(preconditioner_sparsity_pattern);

        solution.reinit(dofs_per_block);
        system_rhs.reinit(dofs_per_block);
    }


    template <int dim>
    void StokesProblem<dim>::assemble_system()
    {
        system_matrix         = 0;
        system_rhs            = 0;
        preconditioner_matrix = 0;

        QGauss<dim> quadrature_formula(degree + 2);

        FEValues<dim> fe_values(fe,
                                quadrature_formula,
                                update_values | update_quadrature_points |
                                update_JxW_values | update_gradients);

        const unsigned int dofs_per_cell = fe.n_dofs_per_cell();

        const unsigned int n_q_points = quadrature_formula.size();

        FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
        FullMatrix<double> local_preconditioner_matrix(dofs_per_cell,
                                                       dofs_per_cell);
        Vector<double>     local_rhs(dofs_per_cell);

        std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);

        const RightHandSide<dim>    right_hand_side;
        std::vector<Tensor<1, dim>> rhs_values(n_q_points, Tensor<1, dim>());

        const FEValuesExtractors::Vector velocities(0);
        const FEValuesExtractors::Scalar pressure(dim);

        std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
        std::vector<double>                  div_phi_u(dofs_per_cell);
        std::vector<Tensor<1, dim>>          phi_u(dofs_per_cell);
        std::vector<double>                  phi_p(dofs_per_cell);

        for (const auto &cell : dof_handler.active_cell_iterators())
        {
            fe_values.reinit(cell);
            local_matrix                = 0;
            local_preconditioner_matrix = 0;
            local_rhs                   = 0;

            right_hand_side.value_list(fe_values.get_quadrature_points(),
                                       rhs_values);

            for (unsigned int q = 0; q < n_q_points; ++q)
            {
                for (unsigned int k = 0; k < dofs_per_cell; ++k)
                {
                    symgrad_phi_u[k] =
                            fe_values[velocities].symmetric_gradient(k, q);
                    div_phi_u[k] = fe_values[velocities].divergence(k, q);
                    phi_u[k]     = fe_values[velocities].value(k, q);
                    phi_p[k]     = fe_values[pressure].value(k, q);
                }

                for (unsigned int i = 0; i < dofs_per_cell; ++i)
                {
                    for (unsigned int j = 0; j <= i; ++j)
                    {
                        local_matrix(i, j) +=
                                (2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) // (1)
                                 - div_phi_u[i] * phi_p[j]                 // (2)
                                 - phi_p[i] * div_phi_u[j])                // (3)
                                * fe_values.JxW(q);                        // * dx

                        local_preconditioner_matrix(i, j) +=
                                (phi_p[i] * phi_p[j]) // (4)
                                * fe_values.JxW(q);   // * dx
                    }
                    local_rhs(i) += phi_u[i]            // phi_u_i(x_q)
                                    * rhs_values[q]     // * f(x_q)
                                    * fe_values.JxW(q); // * dx
                }
            }

            for (unsigned int i = 0; i < dofs_per_cell; ++i)
                for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
                {
                    local_matrix(i, j) = local_matrix(j, i);
                    local_preconditioner_matrix(i, j) =
                            local_preconditioner_matrix(j, i);
                }

            cell->get_dof_indices(local_dof_indices);
            constraints.distribute_local_to_global(local_matrix,
                                                   local_rhs,
                                                   local_dof_indices,
                                                   system_matrix,
                                                   system_rhs);
            constraints.distribute_local_to_global(local_preconditioner_matrix,
                                                   local_dof_indices,
                                                   preconditioner_matrix);
        }

        std::cout << "   Computing preconditioner..." << std::endl << std::flush;

        A_preconditioner =
                std::make_shared<typename InnerPreconditioner<dim>::type>();
        A_preconditioner->initialize(
                system_matrix.block(0, 0),
                typename InnerPreconditioner<dim>::type::AdditionalData());
    }
    template <int dim>
    void StokesProblem<dim>::solve()
    {
        const InverseMatrix<SparseMatrix<double>,
                typename InnerPreconditioner<dim>::type>
                A_inverse(system_matrix.block(0, 0), *A_preconditioner);
        Vector<double> tmp(solution.block(0).size());

        {
            Vector<double> schur_rhs(solution.block(1).size());
            A_inverse.vmult(tmp, system_rhs.block(0));
            system_matrix.block(1, 0).vmult(schur_rhs, tmp);
            schur_rhs -= system_rhs.block(1);

            SchurComplement<typename InnerPreconditioner<dim>::type> schur_complement(
                    system_matrix, A_inverse);

            SolverControl            solver_control(solution.block(1).size(),
                                                    1e-6 * schur_rhs.l2_norm());
            SolverCG<Vector<double>> cg(solver_control);

            SparseILU<double> preconditioner;
            preconditioner.initialize(preconditioner_matrix.block(1, 1),
                                      SparseILU<double>::AdditionalData());

            InverseMatrix<SparseMatrix<double>, SparseILU<double>> m_inverse(
                    preconditioner_matrix.block(1, 1), preconditioner);

            cg.solve(schur_complement, solution.block(1), schur_rhs, m_inverse);

            constraints.distribute(solution);

            std::cout << "  " << solver_control.last_step()
                      << " outer CG Schur complement iterations for pressure"
                      << std::endl;
        }

        {
            system_matrix.block(0, 1).vmult(tmp, solution.block(1));
            tmp *= -1;
            tmp += system_rhs.block(0);

            A_inverse.vmult(solution.block(0), tmp);

            constraints.distribute(solution);
        }
    }




    template <int dim>
    void
    StokesProblem<dim>::output_results(const unsigned int refinement_cycle) const
    {
        std::vector<std::string> solution_names(dim, "velocity");
        solution_names.emplace_back("pressure");

        std::vector<DataComponentInterpretation::DataComponentInterpretation>
                data_component_interpretation(
                dim, DataComponentInterpretation::component_is_part_of_vector);
        data_component_interpretation.push_back(
                DataComponentInterpretation::component_is_scalar);

        DataOut<dim> data_out;
        data_out.attach_dof_handler(dof_handler);
        data_out.add_data_vector(solution,
                                 solution_names,
                                 DataOut<dim>::type_dof_data,
                                 data_component_interpretation);
        data_out.build_patches();

        std::ofstream output(
                "solution-" + Utilities::int_to_string(refinement_cycle, 2) + ".vtk");
        data_out.write_vtk(output);
    }



    template <int dim>
    void StokesProblem<dim>::refine_mesh()
    {
        Vector<float> estimated_error_per_cell(triangulation.n_active_cells());

        const FEValuesExtractors::Scalar pressure(dim);
        KellyErrorEstimator<dim>::estimate(
                dof_handler,
                QGauss<dim - 1>(degree + 1),
                std::map<types::boundary_id, const Function<dim> *>(),
                solution,
                estimated_error_per_cell,
                fe.component_mask(pressure));

        GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                        estimated_error_per_cell,
                                                        0.3,
                                                        0.0);
        triangulation.execute_coarsening_and_refinement();
    }



    template <int dim>
    void StokesProblem<dim>::run()
    {
        {
            std::vector<unsigned int> subdivisions(dim, 1);
            subdivisions[0] = 4;

            const Point<dim> bottom_left = (dim == 2 ?
                                            Point<dim>(-2, -1) :    // 2d case
                                            Point<dim>(-2, 0, -1)); // 3d case

            const Point<dim> top_right = (dim == 2 ?
                                          Point<dim>(2, 0) :    // 2d case
                                          Point<dim>(2, 1, 0)); // 3d case

            GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                      subdivisions,
                                                      bottom_left,
                                                      top_right);
        }

        for (const auto &cell : triangulation.active_cell_iterators())
            for (const auto &face : cell->face_iterators())
                if (face->center()[dim - 1] == 0)
                    face->set_all_boundary_ids(1);


        triangulation.refine_global(4 - dim);

        for (unsigned int refinement_cycle = 0; refinement_cycle < 6;
             ++refinement_cycle)
        {
            std::cout << "Refinement cycle " << refinement_cycle << std::endl;

            if (refinement_cycle > 0)
                refine_mesh();

            setup_dofs();

            std::cout << "   Assembling..." << std::endl << std::flush;
            assemble_system();

            std::cout << "   Solving..." << std::flush;
            solve();

            output_results(refinement_cycle);

            std::cout << std::endl;
        }
    }



    template class StokesProblem<2>;
    template class StokesProblem<3>;

}