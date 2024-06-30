//
// Created by Florian on 16.06.2024.
//
#include "StokesProblem.h"
#include "SchurComplement.h"
#include "RightHandSide.h"
#include "BoundaryValues.h"



namespace project {
    using namespace dealii;


  template <int dim>
  StokesProblem<dim>::StokesProblem(const unsigned int degree)
    : degree(degree)
    , triangulation(Triangulation<dim>::maximum_smoothing)
    //removed , fe(FE_Q<dim>(degree + 1) ^ dim, FE_Q<dim>(degree))
    , dof_handler(triangulation)
    , level_set_dof_handler(triangulation)
    , fe_level_set(degree)
    , mesh_classifier(level_set_dof_handler, level_set)
  {}
  
   template<int dim>
   double double_contraction(dealii::Tensor<2,dim>& t1, dealii::Tensor<2,dim>& t2) {
       double sum = 0;
       //for(int j = 0; j < 2; j++){
       for(int i = 0; i < 2; i++){
           sum += t1[i] * t2[i];
      // }
       }
       return sum;
   }
 
  
   enum ActiveFEIndex
  {
    lagrange = 0,
    nothing  = 1
  };

    template <int dim>
    void StokesProblem<dim>::setup_discrete_level_set()
    {
        std::cout << "Setting up discrete level set function" << std::endl;

        level_set_dof_handler.distribute_dofs(fe_level_set);
        level_set.reinit(level_set_dof_handler.n_dofs());

        const Functions::SignedDistance::Sphere<dim> signed_distance_sphere;
        VectorTools::interpolate(level_set_dof_handler,
                                 signed_distance_sphere,
                                 level_set);
    }
 

  template <int dim>
  void StokesProblem<dim>::setup_dofs()
  {
    A_preconditioner.reset();
    system_matrix.clear();
    preconditioner_matrix.clear();
 
    //removed dof_handler.distribute_dofs(fe);

    
    const auto face_has_flux_coupling = [&](const auto        &cell,
                                            const unsigned int face_index) {
      return this->face_has_ghost_penalty(cell, face_index);
    };
 
    
    //2
    level_set_dof_handler.distribute_dofs(fe_level_set);
    level_set.reinit(level_set_dof_handler.n_dofs());
 
    Point<dim> center(0,0);
    const Functions::SignedDistance::Sphere<dim> signed_distance_sphere(center,1);
    VectorTools::interpolate(level_set_dof_handler,
                             signed_distance_sphere,
                             level_set);
    
    fe_collection.push_back(FESystem<dim>(FE_Q<dim>(degree + 1) ^ dim, FE_Q<dim>(degree)));
    fe_collection.push_back(FESystem<dim>(FE_Nothing<dim>() ^ dim, FE_Nothing<dim>()));
    
    mesh_classifier.reclassify();

    for (const auto &cell : dof_handler.active_cell_iterators())
      {
        const NonMatching::LocationToLevelSet cell_location = mesh_classifier.location_to_level_set(cell);
 
        if (cell_location == NonMatching::LocationToLevelSet::outside)
          cell->set_active_fe_index(ActiveFEIndex::nothing);
        else
          cell->set_active_fe_index(ActiveFEIndex::lagrange);
      }
 
    dof_handler.distribute_dofs(fe_collection);   
        DoFRenumbering::Cuthill_McKee(dof_handler);                       
    //2end*/
 
    std::vector<unsigned int> block_component(dim + 1, 0);
    block_component[dim] = 1;
    DoFRenumbering::component_wise(dof_handler, block_component);
 
    {
      constraints.clear();
 
      const FEValuesExtractors::Vector velocities(0);
      /*DoFTools::make_hanging_node_constraints(dof_handler, constraints);*/
      /*VectorTools::interpolate_boundary_values(dof_handler,
                                               1,
                                               BoundaryValues<dim>(),
                                               constraints,
                                               fe_collection[0].component_mask(velocities));*/
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
      /*const auto face_has_flux_coupling = [&](const auto        &cell,
                                              const unsigned int face_index) {
        return this->face_has_ghost_penalty(cell, face_index);
      };*/
  
      BlockDynamicSparsityPattern dsp(dofs_per_block, dofs_per_block);
 
      Table<2, DoFTools::Coupling> coupling(dim + 1, dim + 1); //solution = |(u_x, ..., v)| = dim + 1
      Table<2, DoFTools::Coupling> face_coupling(dim + 1, dim + 1);
      for (unsigned int c = 0; c < dim + 1; ++c)
        for (unsigned int d = 0; d < dim + 1; ++d)
          if (!((c == dim) && (d == dim))) {
            coupling[c][d] = DoFTools::always;
            face_coupling[c][d] = DoFTools::always;
            }
          else {
            coupling[c][d] = DoFTools::none;
            face_coupling[c][d] = DoFTools::always;
         }
					
     const bool                      keep_constrained_dofs = true;
 
    DoFTools::make_flux_sparsity_pattern(dof_handler,
                                         dsp,
                                         constraints,
                                         keep_constrained_dofs,
                                         coupling,
                                         face_coupling,
                                         numbers::invalid_subdomain_id,
                                         face_has_flux_coupling);
                                         
	  
        /*DoFTools::make_sparsity_pattern(
          dof_handler, coupling, dsp, constraints, false);*/
 
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
  bool StokesProblem<dim>::face_has_ghost_penalty(
    const typename Triangulation<dim>::active_cell_iterator &cell,
    const unsigned int                                       face_index) const
  {
    if (cell->at_boundary(face_index))
      return false;
 
    const NonMatching::LocationToLevelSet cell_location =
      mesh_classifier.location_to_level_set(cell);
 
    const NonMatching::LocationToLevelSet neighbor_location =
      mesh_classifier.location_to_level_set(cell->neighbor(face_index));
 
    if (cell_location == NonMatching::LocationToLevelSet::intersected &&
        neighbor_location != NonMatching::LocationToLevelSet::outside)
      return true;
 
    if (neighbor_location == NonMatching::LocationToLevelSet::intersected &&
        cell_location != NonMatching::LocationToLevelSet::outside)
      return true;
 
    return false;
  }

  template <int dim>
  void StokesProblem<dim>::assemble_system()
  {
    system_matrix         = 0;
    system_rhs            = 0;
    preconditioner_matrix = 0;
 
    const QGauss<dim> quadrature_formula(degree + 2);
 
    /*removedFEValues<dim> fe_values(fe,
                            quadrature_formula,
                            update_values | update_quadrature_points |
                              update_JxW_values | update_gradients);*/
 
    const unsigned int dofs_per_cell = fe_collection[0].n_dofs_per_cell();
 
    const unsigned int n_q_points = quadrature_formula.size();
 
    FullMatrix<double> local_matrix(dofs_per_cell, dofs_per_cell);
    FullMatrix<double> local_preconditioner_matrix(dofs_per_cell,
                                                   dofs_per_cell);
    Vector<double>     local_rhs(dofs_per_cell);
 
    std::vector<types::global_dof_index> local_dof_indices(dofs_per_cell);
 
    const RightHandSide<dim>    right_hand_side;
    std::vector<Tensor<1, dim>> rhs_values(n_q_points, Tensor<1, dim>());
    
    //Added{
    const double ghost_parameter   = 0.5;
    const double nitsche_parameter = 5 * (degree + 1) * degree;
    
    const QGauss<dim - 1>  face_quadrature(degree + 1);
    FEInterfaceValues<dim> fe_interface_values(fe_collection[0],
                                               face_quadrature,
                                               update_gradients |
                                                 update_JxW_values | update_hessians |
                                                 update_normal_vectors);
 
 
    const QGauss<1> quadrature_1D(degree + 1);
 
    NonMatching::RegionUpdateFlags region_update_flags;
    region_update_flags.inside = update_values | update_gradients |
                                 update_JxW_values | update_quadrature_points;
    region_update_flags.surface = update_values | update_gradients |
                                  update_JxW_values | update_quadrature_points |
                                  update_normal_vectors;
 
    NonMatching::FEValues<dim> non_matching_fe_values(fe_collection,
                                                      quadrature_1D,
                                                      region_update_flags,
                                                      mesh_classifier,
                                                      level_set_dof_handler,
                                                      level_set);
     //Added}
 
    const FEValuesExtractors::Vector velocities(0);
    const FEValuesExtractors::Scalar pressure(dim);
 
     std::vector<SymmetricTensor<2, dim>> symgrad_phi_u(dofs_per_cell);
    std::vector<Tensor<2, dim, double>> grad_phi_u(dofs_per_cell);
    std::vector<double>                  div_phi_u(dofs_per_cell);
    std::vector<Tensor<1, dim, double>>          phi_u(dofs_per_cell);
    std::vector<double>                  phi_p(dofs_per_cell);
    
    std::vector<Tensor<1, dim, double>>  grad_phi_p(dofs_per_cell);

 	int count = 0;
    for (const auto &cell : dof_handler.active_cell_iterators() | IteratorFilters::ActiveFEIndexEqualTo(ActiveFEIndex::lagrange))
      {
 	//std::cout << "Cell number " << count++ << std::endl;
        local_matrix                = 0;
        local_preconditioner_matrix = 0;
        local_rhs                   = 0;
                                   
        //added 
        const double cell_side_length = cell->minimum_vertex_distance();
        
        non_matching_fe_values.reinit(cell);
        
          const std::optional<FEValues<dim>> &inside_fe_values = non_matching_fe_values.get_inside_fe_values();
          

 	if(inside_fe_values) {
 	   right_hand_side.value_list(inside_fe_values->get_quadrature_points(), rhs_values);

          for (const unsigned int q :inside_fe_values->quadrature_point_indices())
          {
            for (unsigned int k = 0; k < dofs_per_cell; ++k)
              {
                         symgrad_phi_u[k] =
                  (*inside_fe_values)[velocities].symmetric_gradient(k, q);
                grad_phi_u[k] = (*inside_fe_values)[velocities].gradient(k, q);
                div_phi_u[k] = (*inside_fe_values)[velocities].divergence(k, q);
                phi_u[k]     = (*inside_fe_values)[velocities].value(k, q);
                phi_p[k]     = (*inside_fe_values)[pressure].value(k, q);
                grad_phi_p[k] = (*inside_fe_values)[pressure].gradient(k, q);
              }
            for (unsigned int i = 0; i < dofs_per_cell; ++i)
              {
                for (unsigned int j = 0; j <= i; ++j)
                  {
                 // std::cout << "innerloop " << i << " " << j << " " <<dofs_per_cell  << std::endl;
                    local_matrix(i, j) +=
                      ((double_contraction(grad_phi_u[i], grad_phi_u[j])) // (1)
                       //(2 * (symgrad_phi_u[i] * symgrad_phi_u[j]) // (1)
                       - (div_phi_u[i] * phi_p[j])                 // (2)
                       - phi_p[i] * div_phi_u[j] )                // (3)
                       //- phi_u[j] * grad_phi_p[i]
                       //+ phi_u[i] * grad_phi_p[j])
                      * inside_fe_values->JxW(q);                        // * dx
                      
                      //if((symgrad_phi_u[i] * symgrad_phi_u[j]) == 0 && phi_p[i] * div_phi_u[j] == 0 && div_phi_u[i] * phi_p[j] == 0) std::cout << "Both 0 " << std::endl;
                      //if((symgrad_phi_u[i] * symgrad_phi_u[j]) == 0 && phi_p[i] * div_phi_u[j] == 0 && div_phi_u[i] * phi_p[j] != 0) std::cout << "B " << std::endl;
                      //if(((symgrad_phi_u[i] * symgrad_phi_u[j]) != 0 || phi_p[i] * div_phi_u[j] != 0) && div_phi_u[i] * phi_p[j] == 0) std::cout << "A " << std::endl;
                      //if(double_contraction(grad_phi_u[i] , grad_phi_u[j]) != 0 && (phi_p[i] * div_phi_u[j] != 0 || div_phi_u[i] * phi_p[j] != 0)) std::cout << "Both non Zero " << std::endl;
                    
                      
 
                    local_preconditioner_matrix(i, j) +=
                      (phi_p[i] * phi_p[j]) // (4)
                      * inside_fe_values->JxW(q);   // * dx
                  }

                local_rhs(i) += phi_u[i]            // phi_u_i(x_q)
                                * rhs_values[q]     // * f(x_q)
                                * inside_fe_values->JxW(q); // * dx
              }
          }
          }
          
          for (unsigned int i = 0; i < dofs_per_cell; ++i)  
            for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
            {
              local_matrix(i, j) = local_matrix(j, i);
              local_preconditioner_matrix(i, j) =
                local_preconditioner_matrix(j, i);
                
                
            }
	

         // std::cout << "t" << std::endl;
          BoundaryValues<dim> boundary_condition;
          
          const std::optional<NonMatching::FEImmersedSurfaceValues<dim>> &surface_fe_values = non_matching_fe_values.get_surface_fe_values();
          
          if (surface_fe_values)
          {
            for (const unsigned int q :
                 surface_fe_values->quadrature_point_indices())
              {
                const Point<dim> &point = surface_fe_values->quadrature_point(q);
                const Tensor<1, dim> &normal = surface_fe_values->normal_vector(q);
                
                //std::cout << point << std::endl;
                
                //std::vector<double> b();
                Vector<double> b(dim + 1);
                std::vector<double> v(dim);

                boundary_condition.vector_value(point, b);
                for(int i = 0; i < dim; i++){
                	v[i] = b[i];
                	//std::cout << v[i] << std::endl;
                }
                
                //std::cout << "surface_fe_values "<< count << std::endl;
                
		dealii::ArrayView<double> b1(v);
		dealii::Tensor<1, dim> g(b1);
                
                for (const unsigned int i : surface_fe_values->dof_indices())
                  {
                    for (const unsigned int j :
                         surface_fe_values->dof_indices())
                      {
                        local_matrix(i, j) +=
                          (-normal * (*surface_fe_values)[velocities].gradient(i, q) *
                             (*surface_fe_values)[velocities].value(j, q) +
                           -normal * (*surface_fe_values)[velocities].gradient(j, q) *
                             (*surface_fe_values)[velocities].value(i, q) +
                           nitsche_parameter / cell_side_length *
                             (*surface_fe_values)[velocities].value(i, q) *
                             (*surface_fe_values)[velocities].value(j, q)) *
                          surface_fe_values->JxW(q);
                      }
                    local_rhs(i) +=
                      g *
                      (nitsche_parameter / cell_side_length *
                         (*surface_fe_values)[velocities].value(i, q) -
                       normal * (*surface_fe_values)[velocities].gradient(i, q)
                      // + normal * (*surface_fe_values)[pressure].value(i, q)
                       ) *
                      surface_fe_values->JxW(q);
                  }
              }
          }
        cell->get_dof_indices(local_dof_indices);
        constraints.distribute_local_to_global(local_matrix,
                                               local_rhs,
                                               local_dof_indices,
                                               system_matrix,
                                               system_rhs);
                                               //std::cout<<"const" <<std::endl;
        constraints.distribute_local_to_global(local_preconditioner_matrix,
                                               local_dof_indices,
                                               preconditioner_matrix);

/*

        for (unsigned int i = 0; i < dofs_per_cell; ++i)
          for (unsigned int j = i + 1; j < dofs_per_cell; ++j)
            {
              local_matrix(i, j) = local_matrix(j, i);
              local_preconditioner_matrix(i, j) =
                local_preconditioner_matrix(j, i);
            }*/
            
          /*for(int i = 0; i < dofs_per_cell; ++i){
          	for(int j = 0; j < dofs_per_cell; ++j){
          		std::cout << std::setprecision(3) << local_matrix(i,j) << " ";
          	}
          	std::cout << std::endl;
          }*/
                               
        for (const unsigned int f : cell->face_indices())
          if (face_has_ghost_penalty(cell, f))
            {
              const unsigned int invalid_subface =
                numbers::invalid_unsigned_int;
 
              fe_interface_values.reinit(cell,
                                         f,
                                         invalid_subface,
                                         cell->neighbor(f),
                                         cell->neighbor_of_neighbor(f),
                                         invalid_subface);
 
              const unsigned int n_interface_dofs =
                fe_interface_values.n_current_interface_dofs();
              FullMatrix<double> local_stabilization(n_interface_dofs,
                                                     n_interface_dofs);
              for (unsigned int q = 0;
                   q < fe_interface_values.n_quadrature_points;
                   ++q)
                {
                  const Tensor<1, dim> normal = fe_interface_values.normal(q);
                  for (unsigned int i = 0; i < n_interface_dofs; ++i)
                    for (unsigned int j = 0; j < n_interface_dofs; ++j)
                      {
                        local_stabilization(i, j) +=
                          1 * pow(cell_side_length, 3) * normal *
                          fe_interface_values[pressure].jump_in_gradients(i, q) *
                          normal *
                          fe_interface_values[pressure].jump_in_gradients(j, q) *
                          fe_interface_values.JxW(q);
                       
                       local_stabilization(i, j) +=
                          0.5 * ghost_parameter * cell_side_length * (normal *
                          fe_interface_values[velocities].jump_in_gradients(i, q)) * 
                          (normal *
                          fe_interface_values[velocities].jump_in_gradients(j, q)) *
                          fe_interface_values.JxW(q) 
                          +   
                           0.5 * ghost_parameter * cell_side_length * cell_side_length * cell_side_length * 0.25 * (
                          normal * (normal * fe_interface_values[velocities].jump_in_hessians(i, q))) * 
                          (
                          normal * (normal * fe_interface_values[velocities].jump_in_hessians(j, q))) *
                          fe_interface_values.JxW(q);
                      }
                }
 
              const std::vector<types::global_dof_index>
                local_interface_dof_indices =
                  fe_interface_values.get_interface_dof_indices();
 
 		//constraints.distribute_local_to_global(local_stabilization, local_interface_dof_indices, system_matrix);
              system_matrix.add(local_interface_dof_indices,
                                   local_stabilization);
            }
      }
 //std::cout << system_matrix.block(0,0).n() << std::endl;
 // std::cout << system_matrix.block(0,0).n_actually_nonzero_elements() << std::endl;
    std::cout << "Computing preconditioner..." << std::endl << std::flush;
 
    A_preconditioner =
      std::make_shared<typename InnerPreconditioner<dim>::type>();
            //std::cout << "done prec1" << std::endl;
    A_preconditioner->initialize(
      system_matrix.block(0, 0),
      typename InnerPreconditioner<dim>::type::AdditionalData());
      
      //std::cout << "done prec" << std::endl;
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

            SolverControl            solver_control(solution.block(1).size() + 150,
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

	data_out.set_cell_selection(
	      [this](const typename Triangulation<dim>::cell_iterator &cell) {
		return cell->is_active() &&
		       mesh_classifier.location_to_level_set(cell) !=
		         NonMatching::LocationToLevelSet::outside;
	      });

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
                 fe_collection[0].component_mask(pressure));

        GridRefinement::refine_and_coarsen_fixed_number(triangulation,
                                                        estimated_error_per_cell,
                                                        0.3,
                                                        0.0);
        triangulation.execute_coarsening_and_refinement();
    }



    template <int dim>
    void StokesProblem<dim>::run() {
        {
            std::vector<unsigned int> subdivisions(dim, 1);
            subdivisions[0] = 4;

            const Point<dim> bottom_left = (dim == 2 ?
                                            Point<dim>(-2, -1) :    // 2d case
                                            Point<dim>(-2, 0, -1)); // 3d case

            const Point<dim> top_right = (dim == 2 ?
                                          Point<dim>(2, 0) :    // 2d case
                                          Point<dim>(2, 1, 0)); // 3d case

            GridGenerator::hyper_cube(triangulation, -1.21, 1.21);
            /* GridGenerator::subdivided_hyper_rectangle(triangulation,
                                                       subdivisions,
                                                       bottom_left,
                                                       top_right);*/
        }

        for (const auto &cell: triangulation.active_cell_iterators())
            for (const auto &face: cell->face_iterators())
                if (face->center()[dim - 1] >= 1.20)
                    face->set_all_boundary_ids(1);


        triangulation.refine_global(1);

        double n_refinements = 5;

        for (unsigned int refinement_cycle = 0; refinement_cycle < n_refinements;
             ++refinement_cycle) {
            std::cout << "Refinement cycle " << refinement_cycle << std::endl;
            triangulation.refine_global(1);
            setup_discrete_level_set();
            std::cout << "Classifying cells" << std::endl;
            mesh_classifier.reclassify();
            setup_dofs();
            assemble_system();
            solve();
            //refine_mesh();
            output_results(refinement_cycle);

        }

    }




    template class StokesProblem<2>;
    template class StokesProblem<3>;

}
