#include "parallelism.hpp"
















int process::main_loop(){
	std::cout << "the basic main loop for a generic process does nothing.  create a derived class." << std::endl;
	
	return SUCCESSFUL;
}









int ubermaster_process::main_loop()
{
	witness_set W;

	
	program_options.splash_screen();
	
	
	//parse the options
	program_options.startup(); // tests for existence of necessary files, etc.
	
	//if desired, display the options
	if (program_options.verbose_level>=3)
		program_options.display_current_options();
	
	
	
	
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	int num_vars = get_num_vars_PPD(solve_options.PPD);
	
	
	W.witnessSetParse(program_options.witness_set_filename,num_vars);
	
	program_options.MPType = MPType;
	
	
	W.get_variable_names();
	W.input_filename = program_options.input_filename;
	
	int incidence_number = get_incidence_number(W,num_vars,
																							program_options,
																							program_options.input_filename);
	
	W.incidence_number = incidence_number;
	
	vertex_set V(num_vars);
	curve_decomposition C;
	surface_decomposition S;
	
	
	
	vec_mp *pi = (vec_mp *) br_malloc(W.dim*sizeof(vec_mp ));
	for (int ii=0; ii<W.dim; ii++) {
		init_vec_mp2(pi[ii],W.num_variables, solve_options.T.AMP_max_prec);
		pi[ii]->size = W.num_variables;
	}
	get_projection(pi, program_options, solve_options, W.num_variables, W.dim);
	

	
	switch (W.dim) {
		case 1:
			// curve
			C.main(V, W, pi, program_options, solve_options);
			
			if (program_options.verbose_level>=2)
				printf("outputting data\n");
			
			
			C.output_main(program_options, V);
			
			break;
			
			
		case 2:
			// surface
			S.main(V, W, pi, program_options, solve_options);

			S.output_main(program_options, V);
			
			break;
			
		default:
			std::cout << "bertini_real not programmed for components of dimension " << W.dim << std::endl;
			break;
	}
	
	
	for (int ii=0; ii<W.dim; ii++)
		clear_vec_mp(pi[ii]);
	
	int sendme = TERMINATE;
	MPI_Bcast(&sendme, 1, MPI_INT, 0, MPI_COMM_WORLD);
	return SUCCESSFUL;
}




























int worker_process::main_loop()
{
	
	int single_int_buffer = 0;
	
	int solver_choice = INITIAL_STATE;
	
	while (solver_choice != TERMINATE) {
		
		MPI_Bcast(&solver_choice, 1, MPI_INT, program_options.head(), MPI_COMM_WORLD);
		std::cout << "worker" << program_options.id() << " received call for help for solver " << enum_lookup(solver_choice) << std::endl;
		
		switch (solver_choice) {
			case NULLSPACE:
				nullspace_slave_entry_point(this->solve_options);
				// call the blabla here
				break;
				
			case MIDPOINT_SOLVER:
//				midpoint_slave_entry_point(this->solve_options);
				// call the blabla here
				break;
				
			case TERMINATE:
				break;
				
			case PARSING:
				MPI_Bcast(&single_int_buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
			default:
				break;
		}
	}
	
	
	
	
	
	return SUCCESSFUL;
}







