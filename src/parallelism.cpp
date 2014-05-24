#include "parallelism.hpp"















//
//int process::main_loop(){
//	std::cout << "the basic main loop for a generic process does nothing.  create a derived class." << std::endl;
//	
//	return SUCCESSFUL;
//}









int ubermaster_process::main_loop()
{
	
	boost::timer::auto_cpu_timer t;
	
	
	
	
	
	program_options.splash_screen();
	
	
	//parse the options
	program_options.startup(); // tests for existence of necessary files, etc.
	
	//if desired, display the options
	if (program_options.verbose_level>=3)
		program_options.display_current_options();
	
	
	
	
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	int num_vars = get_num_vars_PPD(solve_options.PPD);
	
	
	witness_data data_mc_data;
	data_mc_data.populate();
	
	data_mc_data.print();

	
	
	witness_set W = data_mc_data.choose(program_options);
	
	if (W.num_points==0) {
		return 1;
	}
	
	program_options.MPType = MPType;
	
	
	W.get_variable_names(num_vars);
	W.input_filename = program_options.input_filename;
	
	W.print_to_screen();
	
	
	
	int incidence_number = get_incidence_number(W.pts_mp[0], program_options, program_options.input_filename);
	
	W.incidence_number = incidence_number;
	
	vertex_set V(num_vars);
	
	
	
	
	
	vec_mp *pi = (vec_mp *) br_malloc(W.dim*sizeof(vec_mp ));
	for (int ii=0; ii<W.dim; ii++) {
		init_vec_mp2(pi[ii],W.num_variables, solve_options.T.AMP_max_prec);
		pi[ii]->size = W.num_variables;
	}
	get_projection(pi, program_options, solve_options, W.num_variables, W.dim);
	
    for (int ii=0; ii<W.dim; ii++) {
        V.add_projection(pi[ii]);
    }
    
	
	
	switch (W.dim) {
		case 1:
			{
				curve_decomposition C;
				
				C.component_num = W.comp_num;
				
				
				std::stringstream converter;
				converter << "_dim_" << C.dimension << "_comp_" << C.component_num;
				program_options.output_dir += converter.str();
				
				// curve
				C.main(V, W, pi, program_options, solve_options);
				
				if (program_options.verbose_level>=2)
					printf("outputting data\n");
				
				
				
				
				C.output_main(program_options.output_dir);
				
				V.print(program_options.output_dir/ "V.vertex");
				
			}
			break;
			
			
		case 2:
		{
			
			surface_decomposition S;
			S.component_num = W.comp_num;
			
			std::stringstream converter;
			converter << "_dim_" << S.dimension << "_comp_" << S.component_num;
			program_options.output_dir += converter.str();
			
			
			// surface
			S.main(V, W, pi, program_options, solve_options);
			
			
			
			
			S.output_main(program_options.output_dir);
			
			V.print(program_options.output_dir/ "V.vertex");
		}
			break;
			
		default:
		{
			std::cout << "bertini_real not programmed for components of dimension " << W.dim << std::endl;
		}
			break;
	}
	
	
	

	
	
	for (int ii=0; ii<W.dim; ii++)
		clear_vec_mp(pi[ii]);
	
	// dismiss the workers
	int sendme = TERMINATE;
	MPI_Bcast(&sendme, 1, MPI_INT, 0, MPI_COMM_WORLD);
	return SUCCESSFUL;
}





























int worker_process::main_loop()
{
	int single_int_buffer = 0;
	
	int solver_choice = INITIAL_STATE;
	
	
	
	while (solver_choice != TERMINATE) {
		
		MPI_Bcast(&solver_choice, 1, MPI_INT, solve_options.head(), MPI_COMM_WORLD);
		
		if ( (solver_choice!=0) && (solve_options.id()==1) && (solve_options.verbose_level>=1)) {
			std::cout << "received call for help for solver " << enum_lookup(solver_choice) << std::endl;
		}
		
		switch (solver_choice) {
			case NULLSPACE:
				nullspace_slave_entry_point(this->solve_options);
				break;
				
			case MIDPOINT_SOLVER:
			{
				surface_decomposition S;
				S.worker_connect(solve_options, program_options);
			}
				break;
				
			case PARSING:
				MPI_Bcast(&single_int_buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
				break;
				
			case MULTILIN:
				multilin_slave_entry_point(this->solve_options);
				break;
				
			case SPHERE_SOLVER:
				sphere_slave_entry_point(this->solve_options);
				break;
				
			case TERMINATE:
				break;
				
			default:
				break;
		}
	}
	
	
	
	
	
	return SUCCESSFUL;
}




