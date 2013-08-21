#include "surface.hpp"



void surface_main(vertex_set & V,
									surface_decomposition & surf,
									witness_set & W_surf,
									vec_mp *pi,
									BR_configuration & program_options,
									solver_configuration & solve_options)
{

#ifdef usetempfolders
	program_options.move_to_temp();
#endif
	
	surf.input_filename = W_surf.input_filename;
	int ambient_dim = 2;
	
	surf.num_variables = W_surf.num_variables;
	surf.component_num = W_surf.comp_num;
	surf.add_projection(pi[0]);
	surf.add_projection(pi[1]);
	
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
	int *randomized_degrees = (int *)br_malloc((W_surf.num_variables-1-ambient_dim)*sizeof(int));
	
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

	W_curve.input_filename = "input_critical_curve";
	
	// this system describes the system for the critical curve
	create_nullspace_system("input_critical_curve",
													boost::filesystem::path(program_options.called_dir) / program_options.input_filename,
													program_options, &ns_config);
	ns_config.clear();
	
	// this should be deprecated
	program_options.current_working_filename = boost::filesystem::absolute("input_critical_curve");
	
	

	
	vec_mp temp_proj; init_vec_mp2(temp_proj,0,1024);
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
	
	
	interslice(W_curve,
								 W_curve_crit,
								 randomizer_matrix,
								 &temp_proj,
								 program_options,
								 solve_options,
								 surf.crit_curve,
								 V);
	
	clear_vec_mp(temp_proj);
	surf.crit_curve.add_projection(pi[0]);
	
	surf.crit_curve.num_variables = num_vars;
	
	// DONE COMPUTING THE CRITICAL CURVE NOW.
	std::cout << "done decomposing critical curve" << std::endl;
	
	
	
	vec_mp *multilin_linears = (vec_mp *) br_malloc(2*sizeof(vec_mp));
	
	
	init_vec_mp2(multilin_linears[0], W_surf.num_variables,1024);
	init_vec_mp2(multilin_linears[1], W_surf.num_variables,1024);
	vec_cp_mp(multilin_linears[0],pi[0]);
	vec_cp_mp(multilin_linears[1],W_surf.L_mp[1]);
	
	vec_mp midpoints_downstairs, crit_downstairs;
	init_vec_mp(midpoints_downstairs,0); init_vec_mp(crit_downstairs,0);
	std::vector< int > index_tracker;
	
	W_curve_crit.compute_downstairs_crit_midpts(crit_downstairs, midpoints_downstairs, index_tracker, pi[0]);
	
	print_point_to_screen_matlab(crit_downstairs,"crit_down");
	print_point_to_screen_matlab(midpoints_downstairs,"middown");
	std::vector< witness_set > midpoint_witness_sets;
	midpoint_witness_sets.resize(midpoints_downstairs->size);
	
	surf.midpoint_slices.resize(midpoints_downstairs->size);
	
	

	
	
	for (int ii=0; ii<midpoints_downstairs->size; ii++){
		std::stringstream converter;
		converter << ii;
		
		int blabla, *declarations;
		parse_input_file(W_surf.input_filename, &blabla);
		partition_parse(&declarations, W_surf.input_filename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		set_mp(&multilin_linears[0]->coord[0], &midpoints_downstairs->coord[ii]);
		
		multilintolin_solver_main(solve_options.T.MPType,
															W_surf,         // witness_set
															randomizer_matrix,
															multilin_linears, //  the set of linears we will solve at.
															&midpoint_witness_sets[ii], // the new data is put here!
															solve_options); // already a pointer
		
		boost::filesystem::path slicename = W_surf.input_filename;
		slicename += "_slice_";
		slicename += converter.str();
		create_sliced_system(W_surf.input_filename, slicename, &multilin_linears[0], 1, W_surf);
		
		
		parse_input_file(slicename, &blabla);
		partition_parse(&declarations, slicename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		
		midpoint_witness_sets[ii].input_filename = slicename;
		midpoint_witness_sets[ii].add_linear(multilin_linears[1]);
		midpoint_witness_sets[ii].add_patch(W_surf.patch_mp[0]);
		midpoint_witness_sets[ii].variable_names = W_surf.variable_names;
		midpoint_witness_sets[ii].dim = 1;
		
		// we already know the component is self-conjugate (by entry condition), so we are free to call this function
		computeCurveSelfConj(midpoint_witness_sets[ii],
												 &pi[1],
												 surf.midpoint_slices[ii],V,
												 num_vars,
												 program_options, solve_options);
		
		surf.midpoint_slices[ii].add_projection(pi[1]);
		surf.midpoint_slices[ii].num_variables = num_vars;
		
		std::cout << "done decomposing the " << ii << "th midpoint slice" << std::endl;

	}

	
	
	
	std::vector< witness_set > critpoint_witness_sets;
	critpoint_witness_sets.resize(crit_downstairs->size);
	
	surf.critpoint_slices.resize(crit_downstairs->size);
	
	for (int ii=0; ii<crit_downstairs->size; ii++){
		std::stringstream converter;
		converter << ii;
		
		int blabla, *declarations;
		parse_input_file(W_surf.input_filename, &blabla);
		partition_parse(&declarations, W_surf.input_filename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		set_mp(&multilin_linears[0]->coord[0], &crit_downstairs->coord[ii]);
		
		multilintolin_solver_main(solve_options.T.MPType,
															W_surf,         // witness_set
															randomizer_matrix,
															multilin_linears, //  the set of linears we will solve at.
															&critpoint_witness_sets[ii], // the new data is put here!
															solve_options); // already a pointer
		
		boost::filesystem::path slicename = W_surf.input_filename;
		slicename += "_critslice_";
		slicename += converter.str();
		create_sliced_system(W_surf.input_filename, slicename, &multilin_linears[0], 1, W_surf);
		
		
		parse_input_file(slicename, &blabla);
		partition_parse(&declarations, slicename, "func_input", "config", 0); // the 0 means not self conjugate.
		free(declarations);																													 // i would like to move this.
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		
		critpoint_witness_sets[ii].input_filename = slicename;
		critpoint_witness_sets[ii].add_linear(multilin_linears[1]);
		critpoint_witness_sets[ii].add_patch(W_surf.patch_mp[0]);
		critpoint_witness_sets[ii].variable_names = W_surf.variable_names;
		critpoint_witness_sets[ii].dim = 1;
		
		// we already know the component is self-conjugate (by entry condition), so we are free to call this function
		computeCurveSelfConj(critpoint_witness_sets[ii],
												 &pi[1],
												 surf.critpoint_slices[ii],V,
												 num_vars,
												 program_options, solve_options);
		
		surf.critpoint_slices[ii].add_projection(pi[1]);
		surf.critpoint_slices[ii].num_variables = num_vars;
		
		std::cout << "done decomposing the " << ii << "th critpoint slice" << std::endl;
		
	}
	
	
	//connect the dots
	connect_the_dots(V,
									 surf,
									 pi,
									 program_options,
									 solve_options);
	
	

#ifdef usetempfolders
	program_options.move_to_called();
#endif
	
	
	
}







void connect_the_dots(vertex_set & V,
											surface_decomposition & surf,
											vec_mp *pi,
											BR_configuration & program_options,
											solver_configuration & solve_options)
{
	
	midpoint_config mid_config;
	
	for (int ii=0; ii<surf.midpoint_slices.size(); ii++) {
		for (int jj=0; jj<surf.midpoint_slices[ii].num_edges; jj++) {
			face F;
			
			//create the face here
			
			
			surf.add_face(F);
		}
	}
	
	return;
}


















void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
													vec_mp * linears, int num_to_add,
													const witness_set & W)
{
	
	if (W.variable_names.size()==0) {
		std::cout << "trying to create a sliced system, but witness set does not have the variable names." << std::endl;
		deliberate_segfault();
	}
	int *declarations = NULL;
	
	partition_parse(&declarations, input_file, "func_input", "config", 0); // the 0 means not self conjugate.
	free(declarations);
	
	FILE *OUT = safe_fopen_write(output_file);
	FILE *IN = safe_fopen_read("func_input");
	
	
	
	
	fprintf(OUT,"INPUT\n\n");
	copyfile(IN,OUT);
	fclose(IN);
	
	std::vector< int > indicators;
	for (int ii=0; ii<num_to_add; ii++) {
		indicators.push_back(rand());
	}
	
	for (int ii=0; ii<num_to_add; ii++) {
		std::stringstream linname;
		linname << "supp_lin_" << indicators[ii];
		write_vector_as_constants(linears[ii], linname.str(), OUT);
		
		linname.clear();  linname.str("");
	}
	
	for (int ii=0; ii<num_to_add; ii++) {
		fprintf(OUT,"function supp_lin_%d;\n",indicators[ii]);
	}
	for (int ii=0; ii<num_to_add; ii++) {
		fprintf(OUT,"supp_lin_%d = supp_lin_%d_1",indicators[ii],indicators[ii]);
		for (int jj=1; jj<W.num_variables; jj++) {
			fprintf(OUT," + %s*supp_lin_%d_%d",W.variable_names[jj].c_str(), indicators[ii],jj+1);
		}
		fprintf(OUT, ";\n\n");
	}
	fprintf(OUT,"END\n");
	
	fclose(OUT);
	
	
	
	
}

