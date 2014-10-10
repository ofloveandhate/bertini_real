#include "parallelism.hpp"
























int ubermaster_process::main_loop()
{
	
	boost::timer::auto_cpu_timer t;
	
	
	program_options.MPType = MPType;
	
	
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
	

	witness_set W = data_mc_data.choose(program_options);
	
	if (W.num_points()==0) {
		return 1;
	}
	
	
	
	
	W.get_variable_names(num_vars);
	W.set_input_filename(program_options.input_filename);
	
	
		
	
	if (program_options.verbose_level>=1) {
		W.print_to_screen();
	}
	
	
	
	
	
	
	vertex_set V(num_vars);
	
	
	
	
	
	vec_mp *pi = (vec_mp *) br_malloc(W.dimension()*sizeof(vec_mp ));
	for (int ii=0; ii<W.dimension(); ii++) {
		init_vec_mp2(pi[ii],W.num_variables(), solve_options.T.AMP_max_prec);
		pi[ii]->size = W.num_variables();
	}
	get_projection(pi, program_options, solve_options, W.num_variables(), W.dimension());
	
    for (int ii=0; ii<W.dimension(); ii++) {
        V.add_projection(pi[ii]);
    }
    
	
	if (program_options.primary_mode==BERTINIREAL) {
		
		
		
		W.set_incidence_number(get_incidence_number( *(W.point(0)), program_options, program_options.input_filename));
		
		
		
		switch (W.dimension()) {
			case 1:
			{
				curve_decomposition C;
				
				C.component_num = W.component_number();
				
				
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
				S.component_num = W.component_number();
				
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
				std::cout << "bertini_real not programmed for components of dimension " << W.dimension() << std::endl;
			}
				break;
		}
		
	}
	else if(program_options.primary_mode==CRIT)
	{
		critreal(W,pi,V);
	}
	
	
	

	
	
	for (int ii=0; ii<W.dimension(); ii++)
		clear_vec_mp(pi[ii]);
	
	// dismiss the workers
	int sendme = TERMINATE;
	MPI_Bcast(&sendme, 1, MPI_INT, 0, MPI_COMM_WORLD);
	return SUCCESSFUL;
}



void ubermaster_process::critreal(witness_set & W, vec_mp *pi, vertex_set & V)
{
	
	
	
	system_randomizer randomizer;
	
	randomizer.setup(W.num_variables()-1-W.dimension(),solve_options.PPD.num_funcs);
	
	
	
	print_point_to_screen_matlab(pi[0],"pi");

	FILE *OUT;
	OUT = safe_fopen_write("most_recent_pi");
	
	print_vec_out_mp(OUT, pi[0]); /* Print vector to file */
	fclose(OUT);
	
	nullspace_config ns_config; // this is set up in the nullspace call.
	
	solver_output solve_out;
	
	
	
	
	compute_crit_nullspace(solve_out, // the returned value
						   W,            // input the original witness set
						   &randomizer,
						   pi,
						   W.dimension(),  // dimension of ambient complex object
						   W.dimension(),   //  target dimension to find
						   W.dimension(),   // COdimension of the critical set to find.
						   program_options,
						   solve_options,
						   &ns_config);
	
	
	witness_set W_crit, W_nonsing, W_sing;
	
	solve_out.get_noninfinite_w_mult_full(W_crit);
	solve_out.get_nonsing_finite_multone(W_nonsing);
	solve_out.get_sing_finite(W_sing);
	
	W_crit.write_dehomogenized_coordinates("crit_data");
	W_nonsing.write_dehomogenized_coordinates("nonsingular");
	W_sing.write_dehomogenized_coordinates("singular");
	
	
	
	return;
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




