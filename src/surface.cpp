#include "surface.hpp"



void surface_main(vertex_set & V,
									surface_decomposition & surf,
									witness_set & W_surf,
									vec_mp *pi,
									BR_configuration & program_options,
									solver_configuration solve_options)
{

#ifdef usetempfolders
	program_options.move_to_temp();
#endif
	
	
	int ambient_dim = 2;
	
	int num_vars = W_surf.num_variables; //HERE CHANGE THIS
	
	
	
	// perform an isosingular deflation
	// note: you do not need witness_data to perform isosingular deflation
	if (program_options.verbose_level>=2)
		printf("performing isosingular deflation\n");
	
	W_surf.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
	int num_deflations, *deflation_sequence = NULL;
	
	isosingular_deflation(&num_deflations, &deflation_sequence, program_options,
																 program_options.current_working_filename,
																 "witness_points_dehomogenized",
																 program_options.max_deflations, W_surf.dim, W_surf.comp_num);
	free(deflation_sequence);
	
	
	program_options.input_deflated_filename = program_options.current_working_filename;
	
	std::stringstream converter;
	converter << "_dim_" << W_surf.dim << "_comp_" << W_surf.comp_num << "_deflated";
	program_options.input_deflated_filename += converter.str();
	converter.clear(); converter.str("");
	
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename);
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	if (program_options.verbose_level>=2) {
		printf("checking if component is self-conjugate\n");
	}
	
	
	checkSelfConjugate(W_surf,num_vars,program_options, program_options.current_working_filename);  //later:  could be passed in from user, if we want
	
	//regenerate the various files, since we ran bertini since then.
	parse_input_file(program_options.input_deflated_filename);
	
	
	
	


	//create the matrix
	mat_mp randomizer_matrix;
	init_mat_mp2(randomizer_matrix,W_surf.num_variables-1-ambient_dim,solve_options.PPD.num_funcs,solve_options.T.AMP_max_prec);
	
	//create the array of integers
	int *randomized_degrees = (int *)bmalloc((W_surf.num_variables-1-ambient_dim)*sizeof(int));
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(randomizer_matrix, &randomized_degrees, W_surf.num_variables-1-ambient_dim, solve_options.PPD.num_funcs);
	
	
	

	
	// find critical points on the critical curve.
	
	
	witness_set W_curve_crit;
	nullspace_config ns_config;
	compute_crit_nullspace(&W_curve_crit, // the returned value
												 W_surf,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 2,   // COdimension of the critical set to find.
												 program_options,
												 solve_options,
												 &ns_config);
	
	W_curve_crit.only_first_vars(num_vars);
	W_curve_crit.sort_for_real(solve_options.T);
	W_curve_crit.write_dehomogenized_coordinates("W_curve_crit"); // write the points to file
	
	ns_config.clear();
	
	
	
	
	
	// compute witness points for the critical curve.  
	witness_set W_curve;
	
	
	compute_crit_nullspace(&W_curve, // the returned value
												 W_surf,            // input the original witness set
												 randomizer_matrix,
												 pi,
												 randomized_degrees,
												 ambient_dim,  // dimension of ambient complex object
												 2,   //  target dimension to find
												 1,   // COdimension of the critical set to find.
												 program_options,
												 solve_options,
												 &ns_config);
	
//	W_curve.only_first_vars(num_vars);
	W_curve.write_dehomogenized_coordinates("W_crit"); // write the points to file

	
	
	// this system describes the system for the critical curve
	create_nullspace_system("input_critical_curve",
													boost::filesystem::path(program_options.called_dir) / program_options.input_filename,
													program_options, &ns_config);
	ns_config.clear();
	
	program_options.current_working_filename = boost::filesystem::absolute("input_critical_curve");
	
	

	
	vec_mp temp_proj; init_vec_mp(temp_proj,0);
	vec_cp_mp(temp_proj, pi[0]);
	std::cout << temp_proj->size << std::endl;
	increase_size_vec_mp(temp_proj, W_curve.num_variables); temp_proj->size = W_curve.num_variables;
	for (int ii=W_surf.num_variables; ii<W_curve.num_variables; ii++) {
		set_zero_mp(&temp_proj->coord[ii]);
	}
	
	for (int ii=0; ii<W_curve.num_linears; ii++) {
		increase_size_vec_mp(W_curve.L_mp[ii], W_curve.num_variables); W_curve.L_mp[ii]->size = W_curve.num_variables;
		
		for (int jj=W_surf.num_variables; jj<W_curve.num_variables; jj++) {
			set_zero_mp(&W_curve.L_mp[ii]->coord[jj]);
		}
	}
	
		

	

	
	slice_and_dice(W_curve,
								 W_curve_crit,
								 randomizer_matrix,
								 &temp_proj,
								 program_options,
								 solve_options,
								 surf.crit_curve,
								 V);
	
	surf.crit_curve.add_projection(pi[0]);
	
	surf.crit_curve.num_variables = num_vars;
	
	// DONE COMPUTING THE CRITICAL CURVE NOW.
	
#ifdef usetempfolders
	program_options.move_to_called();
#endif
	
	
	
}


