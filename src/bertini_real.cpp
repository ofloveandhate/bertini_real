#include "bertini_real.hpp"

int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{

	
	witness_set W;
	
	
	int MPType;
	
	
	////
	//  INITIALIZATION
	////
	MPI_Init(&argC,&args);
	
	
	//instantiate options
	BR_configuration program_options;
	program_options.splash_screen();
	
	program_options.parse_commandline(argC, args);
	
	//parse the options
	program_options.startup(); // tests for existence of necessary files, etc.
	
	//if desired, display the options
	if (program_options.verbose_level>=3)
		program_options.display_current_options();
	
	
	srand(time(NULL));
	
	
	// split the input_file.  this must be called before setting up the solver config.
	parse_input_file(program_options.input_filename, &MPType);
	
	
	// set up the solver configuration
	solver_configuration solve_options;  solver_init_config(&solve_options);
	get_tracker_config(&solve_options,MPType);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	solve_options.verbose_level = program_options.verbose_level;
	solve_options.T.ratioTol = 0.999; // manually assert to be more permissive.
	solve_options.use_gamma_trick = program_options.use_gamma_trick;
	
	if (program_options.verbose_level>=2)
		solve_options.show_status_summary = 1;
	else
		solve_options.show_status_summary = 0;
	
	
	initMP(solve_options.T.Precision); // set up some globals.
	
	
	int num_vars = get_num_vars_PPD(solve_options.PPD);
	
	
	W.witnessSetParse(program_options.witness_set_filename,num_vars);

	W.MPType = program_options.MPType = MPType;
	

	W.get_variable_names();
	
	
	
	int incidence_number = get_incidence_number(W,num_vars,
																							program_options.input_filename,
																							program_options.stifle_text);
	
	W.incidence_number = incidence_number;
	
	vertex_set V;
	curve_decomposition C;
	surface_decomposition S;
	

		
	vec_mp *pi = (vec_mp *) br_malloc(W.dim*sizeof(vec_mp ));
	for (int ii=0; ii<W.dim; ii++) {
		init_vec_mp2(pi[ii],W.num_variables, solve_options.T.AMP_max_prec);
		pi[ii]->size = W.num_variables;
	}
	get_projection(pi, program_options, solve_options, W.num_variables, W.dim);
	
//	for (int ii=0; ii<W.dim; ii++) {
//		print_point_to_screen_matlab(pi[ii],"pi");
//	}
	
	switch (W.dim) {
		case 1:
			// curve
			curve_main(V, C, W, pi, &program_options, &solve_options);
			
			if (program_options.verbose_level>=2)
				printf("outputting data\n");
			
			Output_Main(program_options, W, C, V);
			
			break;
			

		case 2:
			// surface
			surface_main(V, S, W, pi, &program_options, &solve_options);

			Output_Main(program_options, W, S.crit_curve, V);// temporarily using the curve output.
			
			break;
			
		default:
			std::cout << "bertini_real not programmed for components of dimension " << W.dim << std::endl;
			break;
	}
	
	
	for (int ii=0; ii<W.dim; ii++)
		clear_vec_mp(pi[ii]);
	
	
	
	MPI_Finalize();
	
	
	return 0;
}



