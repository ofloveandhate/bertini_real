#include "surface.hpp"



void surface_main(vertex_set & V,
									surface_decomposition & surf,
									witness_set & W,
									vec_mp *pi,
									BR_configuration *program_options,
									solver_configuration *solve_options)
{

	int ambient_dim = 2;
	
	int num_vars = W.num_variables; //HERE CHANGE THIS
	
	
	
	// perform an isosingular deflation
	// note: you do not need witness_data to perform isosingular deflation
	if (program_options->verbose_level>=2)
		printf("performing isosingular deflation\n");
	
	W.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
	int num_deflations, *deflation_sequence = NULL;
	int rV = isosingular_deflation(&num_deflations, &deflation_sequence, program_options->current_working_filename,
																 "witness_points_dehomogenized", "bertini", "matlab -nosplash",
																 program_options->max_deflations, W.dim, W.comp_num);
	free(deflation_sequence);
	
	
	program_options->input_deflated_filename = program_options->current_working_filename;
	
	std::stringstream converter;
	converter << "_dim_" << W.dim << "_comp_" << W.comp_num << "_deflated";
	program_options->input_deflated_filename += converter.str();
	converter.clear(); converter.str("");
	
	// this wraps around a bertini routine
	parse_input_file(program_options->input_deflated_filename);
	
	preproc_data_clear(&solve_options->PPD);
	parse_preproc_data("preproc_data", &solve_options->PPD);
	
	
	
	if (program_options->verbose_level>=2) {
		printf("checking if component is self-conjugate\n");
	}
	int self_conjugate = checkSelfConjugate(W,num_vars,program_options->current_working_filename, program_options->stifle_text);  //later:  could be passed in from user, if we want
	
	//regenerate the various files, since we ran bertini since then.
	parse_input_file(program_options->input_deflated_filename);
	
	
	
	


	//create the matrix
	mat_mp randomizer_matrix;
	init_mat_mp2(randomizer_matrix,W.num_variables-1-ambient_dim,solve_options->PPD.num_funcs,solve_options->T.AMP_max_prec);
	
	//create the array of integers
	int *randomized_degrees = (int *)bmalloc((W.num_variables-1-ambient_dim)*sizeof(int));
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(randomizer_matrix, &randomized_degrees, W.num_variables-1-ambient_dim, solve_options->PPD.num_funcs);
	
	
	

	
	
	// 2(b) compute NID of critical curve, deflating if necessary.
	witness_set W_crit;
	
	nullspace_config ns_config;
	compute_crit_nullspace(&W_crit, // the returned value
												 W,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 1,   // COdimension of the critical set to find.
												 program_options,
												 solve_options,
												 &ns_config);
	
	W_crit.write_dehomogenized_coordinates("W_crit"); // write the points to file
	

	program_options->move_to_temp();
	
	create_nullspace_system("nullspace_jac_system",
													boost::filesystem::path(program_options->called_dir) /= program_options->input_filename,
													program_options, &ns_config);
	ns_config.clear();
	
	
	// W_crit contains a witness SUPERset for the critical curve.
	// somehow supposed to perform NID here on the crit curve.  there may be more than one separate component, some of which need to be deflated, and some which do not.
	
	
//	W_crit.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
//	rV = isosingular_deflation(&num_deflations, &deflation_sequence,
//														 "nullspace_jac_system",
//														 "witness_points_dehomogenized", "bertini", "matlab -nosplash", program_options->max_deflations, W.dim, W.comp_num);
//	free(deflation_sequence);
	
	program_options->current_working_filename = "nullspace_jac_system";
	
	vec_mp temp_proj; init_vec_mp(temp_proj,0);
	vec_cp_mp(temp_proj, pi[1]);
	std::cout << temp_proj->size << std::endl;
	change_size_vec_mp(temp_proj, W_crit.num_variables); temp_proj->size = W_crit.num_variables;
	for (int ii=W.num_variables; ii<W_crit.num_variables; ii++) {
		set_zero_mp(&temp_proj->coord[ii]);
	}
	
	for (int ii=0; ii<W_crit.num_linears; ii++) {
		increase_size_vec_mp(W_crit.L_mp[ii], W_crit.num_variables); W_crit.L_mp[ii]->size = W_crit.num_variables;
		
		for (int jj=W.num_variables; jj<W_crit.num_variables; jj++) {
			set_zero_mp(&W_crit.L_mp[ii]->coord[jj]);
		}
	}
	

	W_crit.print_to_screen();
	std::cout << "about to enter curve decomp for critical curve\n";
//	mypause();
	
	curve_main(V,
						 surf.crit_curve,
						 W_crit,
						 &temp_proj,
						 program_options,
						 solve_options);
	
	program_options->move_to_called();
	
	// now we have a file named 'nullspace_jac_system' which can be used to generate the curve decomposition
	

	std::cout << "done with critical curve" << std::endl;
	exit(0);
	
		
	// 2(c) Find singular points \emph{not} on the critical curve, and keep them in the set of vertices.
	
	
	witness_set W_crit_real2;
	
	compute_crit_nullspace(&W_crit_real2, // the returned value
												 W,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 2,   // COdimension of the critical set to find.
												 program_options,
												 solve_options,
												 &ns_config);
	
	W_crit_real2.write_dehomogenized_coordinates("W_crit_real2"); // write the points to file
	
	program_options->move_to_temp();
	
	create_nullspace_system("nullspace_jac_system2", boost::filesystem::path(program_options->called_dir) /= program_options->input_filename, program_options, &ns_config);
	ns_config.clear();
	
	
	program_options->move_to_called();
	
	
	
	
}


