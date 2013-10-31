#include "surface.hpp"












void surface_decomposition::main(vertex_set & V,
                                 const witness_set & W_surf,
                                 vec_mp *pi,
                                 BR_configuration & program_options,
                                 solver_configuration & solve_options)
{
    
	
    
	
#ifdef usetempfolders
	program_options.move_to_temp();
#endif
	
	// set some member information.
	this->input_filename = W_surf.input_filename;
	this->num_variables = W_surf.num_variables;
	this->component_num = W_surf.comp_num;
    
	add_projection(pi[0]); // add to this
	add_projection(pi[1]); // add to this
	
	
    
	//copy the patch over into this object
	for (int ii=0; ii<W_surf.num_patches; ii++)
		this->add_patch(W_surf.patch_mp[ii]);
	
	
	this->beginning_stuff( W_surf, program_options, solve_options);
	
	
    
    // get the witness points for the critical curve.
    witness_set W_critcurve;
    compute_critcurve_witness_set(W_critcurve,
                                  W_surf,
                                  program_options,
                                  solve_options);
    
    
    // get the critical points and the sphere intersection points for the critical curve
    witness_set W_critcurve_crit;
    compute_critcurve_critpts(W_critcurve_crit, // the computed value
                              W_surf, // input witness set
                              W_critcurve,
                              program_options,
                              solve_options);
    
    
    this->crit_curve.add_witness_set(W_critcurve_crit,CRITICAL,V);
    
    
    
    
    
    
    
    
    
    
    
    // now we get the critical points for the sphere intersection curve.

    // make the input file
	create_sphere_system(W_surf.input_filename,
                         "input_surf_sphere",
                         sphere_radius,
                         sphere_center,
                         W_surf);
    
    // get witness points
	witness_set W_intersection_sphere;
	compute_sphere_witness_set(W_surf,
                               W_intersection_sphere, // output
                               program_options,
                               solve_options);
	
    
    // compute critical points
    witness_set W_sphere_crit;
	compute_sphere_crit(W_intersection_sphere,
                        W_sphere_crit, // output
                        program_options,
                        solve_options);
    
    this->sphere_curve.add_witness_set(W_sphere_crit,CRITICAL,V);
    
    
    // merge together the critical points from both the critical curve and the sphere intersection curve.
    witness_set W_total_crit;
    
	W_total_crit.merge(W_critcurve_crit);
	W_total_crit.only_first_vars( this->num_variables); // chop down the number of variables
    
    
	W_sphere_crit.only_first_vars(this->num_variables);
    W_total_crit.merge(W_sphere_crit);
    
    
    
    W_total_crit.input_filename = "total_crit___there-is-a-problem";
	W_total_crit.sort_for_unique(solve_options.T);
    
    
    
    
    
    
    
    
    
    
    
    

	compute_critical_curve(W_critcurve, // all input.
                           W_total_crit,
                           V,
                           program_options,
                           solve_options);
	
	// DONE COMPUTING THE CRITICAL CURVE NOW.
	this->output_main(program_options, V);
	
	
	
    
	
	
	
	

	
	
	
    // actually perform the interslice on the bounding sphere curve.
	compute_bounding_sphere(W_intersection_sphere, // the witness points we will track from
                            W_total_crit,          // the critical points for the both sphere and critical curve.
                            V,                     // vertex set.  it goes almost everywhere.
                            program_options, solve_options); // configuration
    
	
    this->output_main(program_options, V);
    

    
    
    
    
    
    
	
	
	
    
	
	
	
	
	// compute the downstairs crit and midpoints for slicing
	vec_mp midpoints_downstairs, crit_downstairs;  std::vector< int > index_tracker;
	init_vec_mp(midpoints_downstairs,0); init_vec_mp(crit_downstairs,0);
    
    V.compute_downstairs_crit_midpts(W_total_crit, crit_downstairs, midpoints_downstairs, index_tracker, pi[0]);
	
	
	
	
	
	
	if (program_options.verbose_level>=0) {
        std::cout << color::green() << "the pi[0] projection values at which we will be slicing:\n\n" << color::console_default();
		print_point_to_screen_matlab(crit_downstairs,"crit_down");
		print_point_to_screen_matlab(midpoints_downstairs,"middown");
        W_total_crit.print_to_screen();
	}
	
    
	
	
	
	
	
	
	// get the midpoint slices
	program_options.merge_edges=true;
	bool rerun_empty = false;
	compute_slices(W_surf, V,
                   midpoints_downstairs, this->mid_slices,
                   program_options, solve_options, rerun_empty,"mid");
	
	
	
	//incremental output
	this->output_main(program_options, V);
	
	// get the critical slices
	
	program_options.merge_edges=true;
	rerun_empty = true;
	compute_slices(W_surf, V,
                   crit_downstairs, this->crit_slices,
                   program_options, solve_options, rerun_empty,"crit");
	
    
	//incremental output
	this->output_main(program_options, V);
	
	
	//connect the dots - the final routine
	connect_the_dots(V, pi, program_options, solve_options);
	
	
	
    
#ifdef usetempfolders
	// this is likely extremely broken.  i'd like to use temp folders, but it's not there yet.
	program_options.move_to_called();
#endif
	
	
	
}


void surface_decomposition::beginning_stuff(const witness_set & W_surf,
                                            BR_configuration & program_options,
                                            solver_configuration & solve_options)
{
	// perform an isosingular deflation
	// note: you do not need witness_data to perform isosingular deflation
	if (program_options.verbose_level>=2)
		printf("performing isosingular deflation\n");
	
	W_surf.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
	int num_deflations, *deflation_sequence = NULL;
	
	isosingular_deflation(&num_deflations, &deflation_sequence, program_options,
                          program_options.input_filename,
                          "witness_points_dehomogenized",
                          program_options.max_deflations, W_surf.dim, W_surf.comp_num);
	free(deflation_sequence);
	
	
	program_options.input_deflated_filename = program_options.input_filename;
	
	std::stringstream converter;
	converter << "_dim_" << W_surf.dim << "_comp_" << W_surf.comp_num << "_deflated";
	program_options.input_deflated_filename += converter.str();
	converter.clear(); converter.str("");
	
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename);
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	if (program_options.verbose_level>=2)
		printf("checking if component is self-conjugate\n");
	
	
	
	checkSelfConjugate(W_surf,num_variables,program_options, W_surf.input_filename);
	
	//regenerate the various files, since we ran bertini since then and many files were deleted.
	parse_input_file(program_options.input_deflated_filename);
	
	
	
	
	
	
	//create the matrix
	init_mat_mp2(this->randomizer_matrix,W_surf.num_variables-1-this->dimension,solve_options.PPD.num_funcs,solve_options.T.AMP_max_prec);
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(this->randomizer_matrix, this->randomized_degrees, W_surf.num_variables-1-dimension, solve_options.PPD.num_funcs);
	
	
	
	
	if (!verify_projection_ok(W_surf,
                              this->randomizer_matrix,
                              this->pi,
                              solve_options))
	{
		std::cout << "the projections being used appear to suffer rank deficiency with Jacobian matrix..." << std::endl;
		mypause();
	}
    
	
}














void surface_decomposition::compute_critcurve_witness_set(witness_set & W_critcurve,
                                                          const witness_set & W_surf,
                                                          BR_configuration & program_options,
                                                          solver_configuration & solve_options)
{
  	
	// find witness points on the critical curve.
	
	
	
	
	nullspace_config ns_config;
    
	
	std::cout << color::bold("m") << "computing witness points for the critical curve" << color::console_default() << std::endl;
    
	compute_crit_nullspace(&W_critcurve, // the returned value
                           W_surf,            // input the original witness set
                           this->randomizer_matrix,
                           this->pi,
                           this->randomized_degrees,
                           2,  // dimension of ambient complex object
                           2,   //  target dimension to find
                           1,   // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	//this uses both pi[0] and pi[1] to compute witness points

	W_critcurve.write_dehomogenized_coordinates("W_curve"); // write the points to file
	
	W_critcurve.input_filename = "input_critical_curve";
	
	// this system describes the system for the critical curve
	create_nullspace_system("input_critical_curve",
                            boost::filesystem::path(program_options.called_dir) / program_options.input_deflated_filename,
                            program_options, &ns_config);
    
    
    
    
    // adjust the size of the linears to match the number of variables.  this should be a method in the witness set data type.
    
    for (int ii=0; ii<W_critcurve.num_linears; ii++) {
		increase_size_vec_mp(W_critcurve.L_mp[ii], W_critcurve.num_variables); W_critcurve.L_mp[ii]->size = W_critcurve.num_variables;
		
		for (int jj=W_surf.num_variables; jj<W_critcurve.num_variables; jj++) {
			set_zero_mp(&W_critcurve.L_mp[ii]->coord[jj]);
		}
	}
    
    
    
    return;
    
}


void surface_decomposition::compute_critcurve_critpts(witness_set & W_critcurve_crit, // the computed value
                                                      const witness_set & W_surf, // input witness set
                                                      const witness_set & W_critcurve,
                                                      BR_configuration & program_options,
                                                      solver_configuration & solve_options)
{
  
    
	nullspace_config ns_config;
	
	std::cout << color::bold("m") << "\ncomputing critical points of the critical curve" <<  color::console_default() << std::endl;
	
	
	compute_crit_nullspace(&W_critcurve_crit, // the returned value
                           W_surf,            // input the original witness set
                           this->randomizer_matrix,
                           this->pi,
                           this->randomized_degrees,
                           2,  // dimension of ambient complex object
                           2,   //  target dimension to find
                           2,   // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	// this will use pi[0] to compute critical points
	
	W_critcurve_crit.only_first_vars(num_variables);
//	W_critcurve_crit.write_dehomogenized_coordinates("W_critcurve_crit_all"); // write the points to file
	W_critcurve_crit.sort_for_real(solve_options.T);
//	W_critcurve_crit.write_dehomogenized_coordinates("W_critcurve_crit"); // write the points to file
	
	ns_config.clear();
	
	
	this->compute_sphere_bounds(W_critcurve_crit); // sets the radius and center in this decomposition.  Must propagate to the constituent decompositions as well.   fortunately, i have a method for that!!!
	crit_curve.copy_sphere_bounds(*this); // copy the bounds into the critcurve.
	
    
    std::cout << color::bold("m") << "intersecting critical curve with sphere" << color::console_default() << std::endl;
    
    W_critcurve_crit.input_filename = "input_critical_curve";
    
    
	witness_set W_sphere_intersection;
    // now get the sphere intersection critical points and ends of the interval
	crit_curve.get_additional_critpts(&W_sphere_intersection,  // the returned value
                                      W_critcurve,       // all else here is input
                                      this->randomizer_matrix,  //
                                      program_options,
                                      solve_options);
    
    
    W_sphere_intersection.only_first_vars(W_surf.num_variables);
    
    W_sphere_intersection.sort_for_real(solve_options.T);
	W_sphere_intersection.sort_for_unique(solve_options.T);
    
    W_critcurve_crit.merge(W_sphere_intersection);
    
    
    
    
    return;
}

void surface_decomposition::compute_critical_curve(const witness_set & W_critcurve,
                                                   const witness_set & W_critpts,
                                                   vertex_set & V,
                                                   BR_configuration & program_options,
                                                   solver_configuration & solve_options)
{
	
    std::cout << color::bold("m") << "interslicing critical curve" << color::console_default() << std::endl;
    
	
    
    
	program_options.merge_edges=false;
	
    
	//fluff up the projection to have 0 entries for all the synthetic variables.
	vec_mp temp_proj; init_vec_mp2(temp_proj,0,1024);
	vec_cp_mp(temp_proj, pi[0]);
	std::cout << temp_proj->size << std::endl;
	increase_size_vec_mp(temp_proj, W_critcurve.num_variables); // nondestructive increase in size
    temp_proj->size = W_critcurve.num_variables;
    
	for ( int ii=this->num_variables; ii<W_critcurve.num_variables; ii++) {
		set_zero_mp(&temp_proj->coord[ii]);
	}
	

    

	
	crit_curve.interslice(W_critcurve,
                          W_critpts,
                          this->randomizer_matrix,
                          &temp_proj,
                          program_options,
                          solve_options,
                          V);
	
	clear_vec_mp(temp_proj);
	
	crit_curve.add_projection(pi[0]);
	
	crit_curve.num_variables = num_variables + num_variables-1;
	
    std::cout << color::magenta() << "done decomposing critical curve" << color::console_default() << std::endl;
	return;
}


















void surface_decomposition::compute_sphere_witness_set(const witness_set & W_surf,
                                                       witness_set & W_intersection_sphere,
                                                       BR_configuration & program_options,
                                                       solver_configuration & solve_options)
{
	//build up the start system
	solve_options.robust = true;
	
	
	int blabla;
	
	
	
	parse_input_file(W_surf.input_filename, &blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
    
    
	
	
	make_randomization_matrix_based_on_degrees(randomizer_matrix, sphere_curve.randomized_degrees,
                                               W_surf.num_variables-W_surf.num_patches-W_surf.num_linears,
                                               solve_options.PPD.num_funcs);
	sphere_curve.randomized_degrees.push_back(2);
	
	
	multilin_config ml_config(solve_options); // copies in the randomizer matrix and sets up the SLP & globals.
	
	ml_config.set_randomizer(randomizer_matrix);
	
	vec_mp *multilin_linears = (vec_mp *) br_malloc(2*sizeof(vec_mp));
	init_vec_mp2(multilin_linears[0],W_surf.num_variables,solve_options.T.AMP_max_prec);
	multilin_linears[0]->size = W_surf.num_variables;
	init_vec_mp2(multilin_linears[1],0,solve_options.T.AMP_max_prec);
	vec_cp_mp(multilin_linears[1],W_surf.L_mp[1]);
	
	witness_set W_sphere;
	sphere_config sp_config(randomizer_matrix);
	
	
	for (int ii=0; ii<2; ii++) {
		
		for (int jj=0; jj<W_surf.num_variables; jj++) {
			get_comp_rand_mp(&multilin_linears[0]->coord[jj]);
		}
		
		vec_cp_mp(sp_config.starting_linear[ii], multilin_linears[0]);
		
		witness_set W_temp;
		
		
		
		solve_options.allow_singular = 0;
		solve_options.complete_witness_set = 0;
		solve_options.allow_multiplicity = 0;
		solve_options.allow_unsuccess = 0;
		
		
		
		multilin_solver_master_entry_point(W_surf,         // witness_set
                                           &W_temp, // the new data is put here!
                                           multilin_linears,
                                           ml_config,
                                           solve_options);
		
		W_sphere.merge(W_temp);
		
	}
	
	clear_vec_mp(multilin_linears[0]);
	clear_vec_mp(multilin_linears[1]);
	free(multilin_linears);
	
	
	W_sphere.add_patch(W_surf.patch_mp[0]);
	W_sphere.add_linear(W_surf.L_mp[1]);
	
	
	
	
//	W_sphere.print_to_screen();
	
	
	
	// need to actually move to the sphere system now.
	std::cout << "getting witness points for sphere intersection curve" << std::endl;
	
	
	sp_config.set_memory(solve_options);
	sp_config.set_center(this->sphere_center);
	sp_config.set_radius(this->sphere_radius);
	
	
	
	solve_options.allow_singular = 0;
	solve_options.complete_witness_set = 1;
	solve_options.allow_multiplicity = 0;
	solve_options.allow_unsuccess = 0;
	
	
	
	
	
	sphere_solver_master_entry_point(W_sphere, &W_intersection_sphere,
                                     sp_config,
                                     solve_options);
	
	
	W_intersection_sphere.input_filename = "input_surf_sphere";
	
}




void surface_decomposition::compute_sphere_crit(const witness_set & W_intersection_sphere,
                                                witness_set & W_sphere_crit,
                                                BR_configuration & program_options,
                                                solver_configuration & solve_options)
{
	int blabla;
	parse_input_file("input_surf_sphere", &blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
//	prog_t SLP;
//	setupProg(&SLP, solve_options.T.Precision, solve_options.T.MPType);
//	
//	comp_mp temp;  init_mp(temp);
//	for (int ii=0; ii<W_intersection_sphere.num_pts; ii++) {
//		dot_product_mp(temp, W_intersection_sphere.pts_mp[ii], W_intersection_sphere.L_mp[0]);
//		print_comp_matlab(temp, "result");
//		
//		
//		set_one_mp(temp);
//		eval_struct_mp e; init_eval_struct_mp(e, 0, 0, 0);
//		evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, W_intersection_sphere.pts_mp[ii], temp, &SLP);
//		
//		print_point_to_screen_matlab(e.funcVals,"func_vals");
//		clear_eval_struct_mp(e);
//	}
//	clear_mp(temp);
//	
//	clearProg(&SLP, solve_options.T.MPType, 1); // 1 means call freeprogeval()
//	mypause();
    
    
    
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(sphere_curve.randomizer_matrix, sphere_curve.randomized_degrees, W_intersection_sphere.num_variables-W_intersection_sphere.num_patches-1, solve_options.PPD.num_funcs);
    
    
	
	
	
	nullspace_config ns_config;
	
	std::cout << "computing critical points of sphere curve" << std::endl;
	
	compute_crit_nullspace(&W_sphere_crit,                   // the returned value
                           W_intersection_sphere,            // input the witness set with linears
                           sphere_curve.randomizer_matrix,
                           &(this->pi[0]),
                           sphere_curve.randomized_degrees,
                           1,                                // dimension of ambient complex object
                           1,                                //  target dimension to find
                           1,                                // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	
	ns_config.clear();
	
	W_sphere_crit.sort_for_unique(solve_options.T);
	W_sphere_crit.sort_for_real(solve_options.T);
//	W_sphere_crit.sort_for_inside_sphere()
    
	W_sphere_crit.input_filename = "input_surf_sphere";
	
	return;
}



void surface_decomposition::compute_bounding_sphere(const witness_set & W_intersection_sphere,
                                                    const witness_set & W_crit,
                                                    vertex_set & V,
                                                    BR_configuration & program_options,
                                                    solver_configuration & solve_options)
{
	
    
	program_options.merge_edges=false;
	
	this->sphere_curve.compute_sphere_bounds(W_crit); // inflate around this thing because the interslice method requires having bounds in place.
	
	std::cout << "entering interslice for sphere" << std::endl;
	// then feed it to the interslice algorithm
	this->sphere_curve.interslice(W_intersection_sphere, // the witness set with a single linear and some patches.
                                  W_crit, // the critical points
                                  sphere_curve.randomizer_matrix,
                                  &pi[0],
                                  program_options,
                                  solve_options,
                                  V);
	
	
	this->sphere_curve.add_projection(pi[0]);
	this->sphere_curve.num_variables = num_variables;
	
	
	
	return;
}















void surface_decomposition::compute_slices(const witness_set W_surf,
                                           vertex_set & V,
                                           vec_mp projection_values_downstairs,
                                           std::vector< curve_decomposition > & slices,
                                           BR_configuration & program_options,
                                           solver_configuration & solve_options,
                                           bool rerun_empty,
                                           std::string kindofslice)
{
	
	
	
	
	vec_mp *multilin_linears = (vec_mp *) br_malloc(2*sizeof(vec_mp));
	
	
	init_vec_mp2(multilin_linears[0], W_surf.num_variables,1024);
	init_vec_mp2(multilin_linears[1], W_surf.num_variables,1024);
	vec_cp_mp(multilin_linears[0],pi[0]);
	vec_cp_mp(multilin_linears[1],W_surf.L_mp[1]);
	
	
	comp_mp rand_perterb;  init_mp(rand_perterb); comp_mp(h);  init_mp(h); comp_d jalk;
	jalk->r = 10*solve_options.T.final_tolerance; jalk->i = 0.0;
	d_to_mp(h, jalk);                 // h = 1e-10
	
	
	
	slices.resize(projection_values_downstairs->size);
	
	int blabla;
	parse_input_file(W_surf.input_filename, & blabla);
	
	solve_options.robust = true;
	
	multilin_config ml_config(solve_options); // copies in the randomizer matrix and sets up the SLP & globals.
	
	for (int ii=0; ii<projection_values_downstairs->size; ii++){
		
		
		std::cout << color::magenta() << "decomposing the " << ii << "th " << kindofslice << " slice" << color::console_default() << std::endl;
		print_comp_matlab(&projection_values_downstairs->coord[ii], "target_proj");
		
		solve_options.backup_tracker_config();
		
		int iterations=0;
		
		witness_set slice_witness_set; // deliberately scoped variable
		
		curve_decomposition new_slice;
		
		while ((new_slice.num_edges==0) && (iterations<10)) {// try not to allow an empty edge, if rerun_empty is true.  break below.
			
			std::stringstream converter;
			converter << ii;
			
			int blabla;
			parse_input_file(W_surf.input_filename, &blabla);
			//			partition_parse(&declarations, W_surf.input_filename, "func_input", "config", 0); // the 0 means not self conjugate.
			//			free(declarations);																													 // i would like to move this->
			preproc_data_clear(&solve_options.PPD); // ugh this sucks
			parse_preproc_data("preproc_data", &solve_options.PPD);
			
			make_randomization_matrix_based_on_degrees(randomizer_matrix, randomized_degrees, W_surf.num_variables-W_surf.num_patches-W_surf.num_linears, solve_options.PPD.num_funcs);
			
			ml_config.set_randomizer(randomizer_matrix);
			neg_mp(&multilin_linears[0]->coord[0], &projection_values_downstairs->coord[ii]);
			
			
			
			solve_options.allow_singular = 0;
			solve_options.complete_witness_set = 0;
			solve_options.allow_multiplicity = 0;
			solve_options.allow_unsuccess = 0;
			
			
			//			print_matrix_to_screen_matlab(ml_config.randomizer_matrix,"rand_surf_slice_314");
			multilin_solver_master_entry_point(W_surf,         // witness_set
                                               &slice_witness_set, // the new data is put here!
                                               multilin_linears,
                                               ml_config,
                                               solve_options);
			
			boost::filesystem::path slicename = W_surf.input_filename;
			slicename += "_"; slicename += kindofslice; slicename += "slice_"; slicename += converter.str();
			create_sliced_system(W_surf.input_filename, slicename, &multilin_linears[0], 1, W_surf);
			
			
			parse_input_file(slicename, &blabla);
			//			partition_parse(&declarations, slicename, "func_input", "config", 0); // the 0 means not self conjugate.
			//			free(declarations);																													 // i would like to move this->
			preproc_data_clear(&solve_options.PPD);
			parse_preproc_data("preproc_data", &solve_options.PPD);
			
			
			slice_witness_set.num_variables = W_surf.num_variables;
			slice_witness_set.input_filename = slicename;
			slice_witness_set.add_linear(multilin_linears[1]);
			slice_witness_set.add_patch(W_surf.patch_mp[0]);
			slice_witness_set.variable_names = W_surf.variable_names;
			slice_witness_set.dim = 1;
			
			
			
			
			new_slice.reset();
			
			solve_options.complete_witness_set = 1;
			
			new_slice.copy_sphere_bounds(*this);
			
			// we already know the component is self-conjugate (by entry condition), so we are free to call this function
			// the memory for the multilin system will get erased in this call...
			new_slice.computeCurveSelfConj(slice_witness_set,
                                           &pi[1],
                                           V,
                                           num_variables,
                                           program_options, solve_options);
			
			
			if (iterations<3) {
				solve_options.T.securityMaxNorm = 2*solve_options.T.securityMaxNorm;
			}
			else if (iterations< 7){
				//				solve_options.T.final_tolerance = 0.5*solve_options.T.final_tolerance;
				solve_options.T.securityMaxNorm = 10*solve_options.T.securityMaxNorm;
				
				jalk->r = 1e-17;
				d_to_mp(h, jalk);                 // h = 1e-10
				
				get_comp_rand_real_mp(rand_perterb); //  rand_perterb = real rand
				
				
				
				mul_mp(rand_perterb,rand_perterb,h); // rand_perterb = small real rand
				add_mp(&projection_values_downstairs->coord[ii],&projection_values_downstairs->coord[ii],rand_perterb); // perturb the projection value
				
				
			}
			else
			{
				jalk->r = 1e-16;
				d_to_mp(h, jalk);                 // h = 1e-10
				
				get_comp_rand_real_mp(rand_perterb); //  rand_perterb = real rand
				
				mul_mp(rand_perterb,rand_perterb,h); // rand_perterb = small real rand
				add_mp(&projection_values_downstairs->coord[ii],&projection_values_downstairs->coord[ii],rand_perterb); // perturb the projection value
				
				//				solve_options.T.final_tolerance = 0.1*solve_options.T.final_tolerance;
				solve_options.T.securityMaxNorm = 100*solve_options.T.securityMaxNorm;
			}
			
			
			
			
			
			if (rerun_empty == false)
				break;
			
			iterations++;
			
			slice_witness_set.reset();
			
			
		} // re: while loop
		solve_options.reset_tracker_config();
		
		new_slice.add_projection(pi[1]);
		new_slice.num_variables = num_variables;
		
		slices[ii] = new_slice;
		
		
		this->output_main(program_options, V);
		std::cout << color::magenta() << "DONE decomposing the " << ii << "th " << kindofslice << " slice" << color::console_default() << std::endl;
		
	} // re: for loop
	
	clear_mp(rand_perterb);  clear_mp(h);
	
	return;
}




















void surface_decomposition::connect_the_dots(vertex_set & V,
                                             vec_mp *pi,
                                             BR_configuration & program_options,
                                             solver_configuration & solve_options)
{
	
	
	std::cout << color::bold("m") << "***\n\nCONNECT THE DOTS\n\n***" << color::console_default() << std::endl;
    
	midpoint_config md_config;
	
	
	md_config.setup(*this, solve_options);
	
	
	
	comp_mp temp, temp2, temp3;
	init_mp(temp);
	init_mp(temp2);
	init_mp(temp3);
	
	comp_mp numer, denom;
	init_mp(numer);
	init_mp(denom);
	
	comp_mp proj_top, proj_bottom, proj_mid;
	init_mp(proj_mid);
	init_mp(proj_bottom);
	init_mp(proj_top);
	
	comp_mp found_proj_0, found_proj_1;  init_mp(found_proj_0); init_mp(found_proj_1);
	
	
	vec_mp found_point;
    init_vec_mp(found_point, this->num_variables);
    found_point->size = this->num_variables;
	
	
	witness_set W_midtrack;
	
    
	vec_mp blank_point;  init_vec_mp2(blank_point, 0,1024);
	W_midtrack.add_point(blank_point);
	
	
    
	
	
	for (int ii=0; ii<mid_slices.size(); ii++) {
		
		
		for (int jj=0; jj<mid_slices[ii].num_edges; jj++) {
			this->output_main(program_options, V);
			std::cout << color::magenta() << "\n\n\n*****************************\nface " << this->num_faces << ", midslice " << ii << " edge " << jj << color::console_default() << std::endl;
			
			if (mid_slices[ii].edges[jj].is_degenerate()) {
				std::cout << "no computation necessary -- midslice edge is degenerate" << std::endl;
				continue;
			}
			
			
			
			
            
			//create the face
			face F;
			
			F.index = ii; // the index of which midslice this face came from.
			F.midpt = mid_slices[ii].edges[jj].midpt; // index the point
			
			
			// perform a search to find the top and bottom edges in the appropriate curve.
			
			std::cout << color::brown() << "current midpoint: " <<  mid_slices[ii].edges[jj].midpt  << " " << color::console_default() << std::endl;
            
			
			
			// get the type of system for the top and bottom edges.  this is determined by reading the system name for the midpoints.
			// info on the files from which the points came from.
			int bottom_index = V.vertices[mid_slices[ii].edges[jj].left].input_filename_index;
			int top_index = V.vertices[mid_slices[ii].edges[jj].right].input_filename_index;
			
			if (V.filenames[bottom_index].string().compare(sphere_curve.input_filename.string())==0) {// note that this comparison will fail if using temp directories.
				std::cout << "bottom point is a sphere point" << std::endl;
				md_config.system_type_bottom = SYSTEM_SPHERE;
			}
			else if (V.filenames[bottom_index].string().compare(crit_curve.input_filename.string())==0) {
				std::cout << "bottom point is a crit point" << std::endl;
				md_config.system_type_bottom = SYSTEM_CRIT;
			}
			else{
                md_config.system_type_bottom = UNSET;
				std::cout << color::red() << "bottom point does not appear to be a point of either sphere or crit type." << std::endl <<
				"instead, it appears to have came from " << V.filenames[bottom_index].string() << color::console_default() << std::endl;
			}
			
			
			if (V.filenames[top_index].string().compare(sphere_curve.input_filename.string())==0) {
				std::cout << "top point is a sphere point" << std::endl;
				md_config.system_type_top = SYSTEM_SPHERE;
			}
			else if (V.filenames[top_index].string().compare(crit_curve.input_filename.string())==0) {
				std::cout << "top point is a crit point" << std::endl;
				md_config.system_type_top = SYSTEM_CRIT;
			}
			else{
                md_config.system_type_top = UNSET;
				std::cout << color::red() << "top point does not appear to be a point of either sphere or crit type." << std::endl <<
				"instead, it appears to have came from " << V.filenames[top_index].string() << color::console_default() << std::endl;
			}
			
			
			
			
			
			
			
			// get the bottom and top edges for this face.
            
			
			int num_bottom_vars;
			if (md_config.system_type_bottom == SYSTEM_CRIT) {
				num_bottom_vars = md_config.num_crit_vars;
				F.bottom = crit_curve.edge_w_midpt(mid_slices[ii].edges[jj].left); // index the *edge*
				F.system_type_bottom = SYSTEM_CRIT;
			}
			else if (md_config.system_type_bottom == SYSTEM_SPHERE){
				num_bottom_vars = md_config.num_sphere_vars;
				F.bottom = sphere_curve.edge_w_midpt(mid_slices[ii].edges[jj].left); // index the *edge*
				F.system_type_bottom = SYSTEM_SPHERE;
			}
			else {
				F.bottom = -10;
				F.system_type_bottom = UNSET;
			}
			
			
			int num_top_vars;
			if (md_config.system_type_top == SYSTEM_CRIT) {
				num_top_vars = md_config.num_crit_vars;
				F.top = crit_curve.edge_w_midpt(mid_slices[ii].edges[jj].right); // index the *edge*
				F.system_type_top = SYSTEM_CRIT;
			}
			else if (md_config.system_type_top == SYSTEM_SPHERE){
				num_top_vars = md_config.num_sphere_vars;
				F.top = sphere_curve.edge_w_midpt(mid_slices[ii].edges[jj].right); // index the *edge*
				F.system_type_top = SYSTEM_SPHERE;
			}
			else {
				F.top = -10;
				F.system_type_top = UNSET;
			}
			
			
			bool bottom_ok = ((md_config.system_type_bottom == SYSTEM_CRIT) || (md_config.system_type_bottom == SYSTEM_SPHERE));
			bool top_ok = ((md_config.system_type_top == SYSTEM_CRIT) || (md_config.system_type_top == SYSTEM_SPHERE));
			
			if (! (bottom_ok && top_ok) ) {
				std::cout << color::red() << "one of the points was neither sphere nor crit in origin." << color::console_default() << std::endl;
				
                std::cout << "tracking from these point indices:" << std::endl;
                std::cout <<  mid_slices[ii].edges[jj].left  << " " << mid_slices[ii].edges[jj].midpt  << " "  << mid_slices[ii].edges[jj].right << std::endl;
                
				add_face(F);
//				mypause();
				continue;
			}
            
			
			
			//copy in the start point as three points concatenated.
			
			W_midtrack.num_variables = this->num_variables + num_bottom_vars + num_top_vars;
			
			change_size_vec_mp(W_midtrack.pts_mp[0], W_midtrack.num_variables); // destructive resize
			W_midtrack.pts_mp[0]->size = W_midtrack.num_variables;
			
			// mid
			int var_counter = 0;
			for (int kk=0; kk<this->num_variables; kk++) {
				set_mp(&W_midtrack.pts_mp[0]->coord[kk], &V.vertices[mid_slices[ii].edges[jj].midpt].pt_mp->coord[kk]);
				var_counter++;
			}
			
			// bottom
			int offset = var_counter;
			for (int kk=0; kk<num_bottom_vars; kk++) {
				set_mp(&W_midtrack.pts_mp[0]->coord[kk+offset], &V.vertices[mid_slices[ii].edges[jj].left].pt_mp->coord[kk]); // y0
				var_counter++;
			}
			
			// top
			offset = var_counter;
			for (int kk=0; kk<num_top_vars; kk++) {
				set_mp(&W_midtrack.pts_mp[0]->coord[kk+offset], &V.vertices[mid_slices[ii].edges[jj].right].pt_mp->coord[kk]); // y2
				var_counter++;
			}
			
			
			
			
			
			
			
			
			
			//copy in the patches appropriate for the systems we will be tracking on.  this could be improved.
			W_midtrack.reset_patches();
			
			for (int qq = 0; qq<this->num_patches; qq++) {
				W_midtrack.add_patch(this->patch[qq]);
			}
			
			if (md_config.system_type_bottom == SYSTEM_CRIT) {
				for (int qq = 0; qq<crit_curve.num_patches; qq++) {
					W_midtrack.add_patch(crit_curve.patch[qq]);
				}
			}
			else{
				for (int qq = 0; qq<sphere_curve.num_patches; qq++) {
					W_midtrack.add_patch(sphere_curve.patch[qq]);
				}
			}
			
			if (md_config.system_type_top == SYSTEM_CRIT) {
				for (int qq = 0; qq<crit_curve.num_patches; qq++) {
					W_midtrack.add_patch(crit_curve.patch[qq]);
				}
			}
			else{
				for (int qq = 0; qq<sphere_curve.num_patches; qq++) {
					W_midtrack.add_patch(sphere_curve.patch[qq]);
				}
			}
			
			
			
			
			
			
			
			// make u, v target values.
            
			projection_value_homogeneous_input(md_config.crit_val_left,  V.vertices[ crit_slices[ii].edges[0].midpt ].pt_mp,pi[0]);
			projection_value_homogeneous_input(md_config.crit_val_right, V.vertices[crit_slices[ii+1].edges[0].midpt].pt_mp,pi[0]);
			
			
			// the u direction corresponds to pi[0].
			for (int zz=0; zz<2; zz++) { // go left (zz=0) and right (zz=1)
				if (zz==0) {
					std::cout << "\n      <<=========   going left" << std::endl;
                    std::cout << "tracking from these point indices:" << std::endl;
                    std::cout <<  mid_slices[ii].edges[jj].left  << " " << mid_slices[ii].edges[jj].midpt  << " "  << mid_slices[ii].edges[jj].right << std::endl;
				}
				else{
					std::cout << "\n\n           going right   =======>> " << std::endl;
                    std::cout << "tracking from these point indices:" << std::endl;
                    std::cout <<  mid_slices[ii].edges[jj].left  << " " << mid_slices[ii].edges[jj].midpt  << " "  << mid_slices[ii].edges[jj].right << std::endl;
				}
				
				
                
				//track
				int final_top_ind, final_bottom_ind; // indexes in V of the bottom and top points of the left or right edge.
				
				
				if (zz==0) {
					if (md_config.system_type_bottom == SYSTEM_CRIT) {
						final_bottom_ind = crit_curve.edges[F.bottom].left;
					}
					else{
						final_bottom_ind = sphere_curve.edges[F.bottom].left;
					}
					
					if (md_config.system_type_top == SYSTEM_CRIT) {
						final_top_ind = crit_curve.edges[F.top].left;
					}
					else{
						final_top_ind = sphere_curve.edges[F.top].left;
					}
					
					set_zero_mp(md_config.u_target);
				}
				else{ // zz==1, and going right
					if (md_config.system_type_bottom == SYSTEM_CRIT) {
						final_bottom_ind = crit_curve.edges[F.bottom].right;
					}
					else{
						final_bottom_ind = sphere_curve.edges[F.bottom].right;
					}
					
					if (md_config.system_type_top == SYSTEM_CRIT) {
						final_top_ind = crit_curve.edges[F.top].right;
					}
					else{
						final_top_ind = sphere_curve.edges[F.top].right;
					}
					
					set_one_mp(md_config.u_target);
				}
				
                
				
				// check for degeneracy
				if (final_bottom_ind==final_top_ind) {
					// can simply set the top or bottom edge to be this one.  know it goes there.
					std::cout << "current edge is degenerate" << std::endl;
					
					
					int current_edge = crit_slices[ii+zz].edge_w_midpt(final_bottom_ind);
					
					if (current_edge<0) {
						std::cout << "unable to find an edge in crit_slices[" << ii+zz << "] with midpoint " << final_bottom_ind << std::endl;
						std::cout << "checking for edges with point as left or right point" << std::endl;
						
						current_edge = crit_slices[ii+zz].edge_w_left(final_bottom_ind);
						
						if (current_edge < 0) {
							current_edge = crit_slices[ii+zz].edge_w_right(final_bottom_ind);
						}
						
						if (current_edge < 0) {
							std::cout << "unable to find any edges in crit_slice with point" << std::endl;
							continue;
						}
					}
					
					
					if (zz==0){
						F.left.push_back(current_edge);
						F.num_left++;
					}
					else {
						F.right.push_back(current_edge);
						F.num_right++;
					}
					continue; // go to next zz value, or next midslice edge, or whatever
				}
				
				std::cout << "final top: " << final_top_ind << ", final bottom:	" << final_bottom_ind << std::endl;
				
				// get the projection values of top and bottom final points.
				projection_value_homogeneous_input(proj_top,    V.vertices[ final_top_ind ].pt_mp,   pi[1]); //w2
				projection_value_homogeneous_input(proj_bottom, V.vertices[ final_bottom_ind ].pt_mp,pi[1]); //w0
				
				//initialize current index trackers.
				int current_bottom_ind = final_bottom_ind;
				int current_top_ind = -12131; // initialize to impossible value;
				
				
				
				// candidates
				std::set< int > found_edges;
				std::set< int > possible_edges;
				
                
				for (int rr = 0; rr< crit_slices[ii+zz].num_edges; rr++){
					possible_edges.insert(rr);}
				
				
				while ((current_top_ind != final_top_ind) && (possible_edges.size()>0)) // int yy=0; yy<crit_slices[ii+zz].num_edges; yy++
				{
					
					std::cout << "target bottom: " << final_bottom_ind << " current bottom: " << current_bottom_ind << " current top: " << current_top_ind << " final top: " << final_top_ind << std::endl;
					
					
					
					std::vector< int > candidates; // indices of candidates for next one.
					
					
					
					int candidate_counter = 0;
					std::cout << "\nfinding candidates for bottom index " << current_bottom_ind << std::endl;
					for (int qq=0; qq< crit_slices[ii+zz].num_edges; qq++) {
						
						
						bool correct_interval = false;
						bool matches_end = ((crit_slices[ii+zz].edges[qq].left == current_bottom_ind) || (crit_slices[ii+zz].edges[qq].right == current_bottom_ind));
						bool havent_found_yet = (possible_edges.find(qq)!=possible_edges.end());
                        
						
						// we gotta be moving from lower to higher...  so temp > temp2 is required
						if (matches_end) {
							projection_value_homogeneous_input(temp, V.vertices[ crit_slices[ii+zz].edges[qq].midpt].pt_mp,pi[1]);
							projection_value_homogeneous_input(temp2, V.vertices[ final_bottom_ind].pt_mp,pi[1]);
							projection_value_homogeneous_input(temp3, V.vertices[ final_top_ind].pt_mp,pi[1]);
							
							correct_interval =  ( mpf_get_d(temp3->r) > mpf_get_d(temp->r)) && (mpf_get_d(temp->r) > mpf_get_d(temp2->r)) ;
						}
                        
						
						if (havent_found_yet && matches_end && correct_interval) {
							candidates.push_back(qq);
							
							std::cout << "candidate [" << candidate_counter << "] = " <<
                            crit_slices[ii+zz].edges[qq].left << " " << crit_slices[ii+zz].edges[qq].midpt << " " << crit_slices[ii+zz].edges[qq].right <<  std::endl;
							
							candidate_counter++;
						}
						else
						{
                            //							std::cout << "edge " << qq << " excluded: " << correct_interval << " dir, (fabs=" << fabs( mpf_get_d(temp->r) - mpf_get_d(temp2->r)) << ") " << matches_end << " matches, " << havent_found_yet << "  ~found yet" << std::endl;
							
						}
                        
					}
					
					std::cout << std::endl;
					
					
					if (candidate_counter==0) {
						std::cout << "found 0 candidates for left endpoint, bottom index " << current_bottom_ind << std::endl;
						break; // out of the while loop
					}
					
					for (int qq=0; qq<candidate_counter; qq++)
					{
						int current_edge = candidates[qq];
						
						std::cout << "face #: " << this->num_faces << ", zz: " << zz << ", current_edge: " << current_edge << std::endl;
						std::cout << "tracking to these indices: " << final_bottom_ind << " " << crit_slices[ii+zz].edges[current_edge].midpt << " " << final_top_ind << std::endl;
                        
						
						
						//target midpoint e.w from paper.
						projection_value_homogeneous_input(proj_mid, V.vertices[ crit_slices[ii+zz].edges[current_edge].midpt ].pt_mp,pi[1]);
                        
						
						if (solve_options.use_real_thresholding) {
							if (fabs(mpf_get_d(proj_top->i))<1e-10)
								mpf_set_str(proj_top->i,"0.0",10);
							
							if (fabs(mpf_get_d(proj_bottom->i))<1e-10)
								mpf_set_str(proj_bottom->i,"0.0",10);
							
							if (fabs(mpf_get_d(proj_mid->i))<1e-10)
								mpf_set_str(proj_mid->i,"0.0",10);
						}
						
						
						
						sub_mp(denom, proj_top, proj_bottom); // p2(w2) - p2(w0);
						sub_mp(numer, proj_mid, proj_bottom); // p2(e.w) - p2(w0);
						div_mp(md_config.v_target, numer, denom); // [p2(e.w) - p2(w0)] / [p2(w2) - p2(w0)]
						
						//print some display to screen
						if (solve_options.verbose_level >= 3)
						{
							print_point_to_screen_matlab(V.vertices[mid_slices[ii].edges[jj].right].pt_mp,"top_start");
							print_point_to_screen_matlab(V.vertices[mid_slices[ii].edges[jj].left].pt_mp,"bottom_start");
							print_point_to_screen_matlab(V.vertices[ crit_slices[ii+zz].edges[current_edge].midpt ].pt_mp,"midpoint_target");
						}
						
						if (solve_options.verbose_level >= 0)
						{
							print_comp_matlab(proj_top,"upper");
							print_comp_matlab(proj_bottom,"lower");
							print_comp_matlab(proj_mid,"mid");
							
							print_comp_matlab(numer,"numer");
							print_comp_matlab(denom,"denom");
							
							print_comp_matlab(md_config.u_target,"u_target");
							print_comp_matlab(md_config.v_target,"v_target");
							
							print_comp_matlab(md_config.crit_val_left,"crit_val_left");
							print_comp_matlab(md_config.crit_val_right,"crit_val_right");
							
							
							set_one_mp(temp);
							sub_mp(temp2, temp, md_config.v_target);
							mul_mp(temp, temp2, proj_bottom);
							
							mul_mp(temp2, md_config.v_target, proj_top);
							
							add_mp(temp3, temp, temp2);
							
							print_comp_matlab(temp3,"proj_1_target_mid");
						}
						
						
						
						solve_options.allow_multiplicity = 1;
						solve_options.allow_singular = 1;
						solve_options.robust = true;
						
						
						witness_set W_new;
						midpoint_solver_master_entry_point(W_midtrack, // carries with it the start points, and the linears.
                                                           &W_new, // new data goes in here
                                                           *this,
                                                           md_config,
                                                           solve_options);
						
						// should get a single point back from this solver.
						
						if (W_new.num_pts==0) {
							std::cout << color::red() << "midpoint tracker did not return any points :(" << color::console_default() << std::endl;
                            possible_edges.erase(current_edge);
							continue;
						}
						
						
						// get only the midpoint coordinates out of the returned point
						for (int tt = 0; tt<this->num_variables; tt++) {
							set_mp(&found_point->coord[tt], &W_new.pts_mp[0]->coord[tt]);
						}
						
						
						//need to look the found point up in vertex set V
						int found_index = index_in_vertices(V,
                                                            found_point);
						
						
						
						
						// find projection value?  remove this?
						projection_value_homogeneous_input(found_proj_1, found_point, pi[1]);
						
						
						if (solve_options.verbose_level>=0) {
							vec_mp tempvec; init_vec_mp(tempvec,0);
							dehomogenize(&tempvec, found_point);
							print_point_to_screen_matlab(tempvec, "found_point");
							clear_vec_mp(tempvec);
							
							std::cout << "found_index of point: " << found_index << std::endl;
							
						}
						
						
						
						
						
						//search among the current edge possibilities for the one containing the found (mid) point
						int index_in_set = -1;
						for (std::set<int>::iterator possibility_iter=possible_edges.begin(); possibility_iter!=possible_edges.end(); possibility_iter++) {
							if (found_index == crit_slices[ii+zz].edges[*possibility_iter].midpt) {
								index_in_set = *possibility_iter;
								break;
							}
						}
						
						
						if (index_in_set < 0) {
							index_in_set = crit_slices[ii+zz].edge_w_removed(found_index);
							if (index_in_set>=0) {
								std::cout << color::green() << "found point as removed point from edge " << index_in_set << color::console_default() << std::endl;
							}
						}
						
						
						if (index_in_set < 0) {
							std::cout << color::red() << "did not find the indexed point as the midpoint of any current possibilitiy." << color::console_default() << std::endl;
                            
							// search edges for the found point as a removed point.
						}
						
						
						
						if (index_in_set>=0)
						{
							
							int next_edge = index_in_set; // index the *edge*
							
							std::cout << "next_edge " << next_edge << ", l m r: " << crit_slices[ii+zz].edges[next_edge].left << " " << crit_slices[ii+zz].edges[next_edge].midpt << " " << crit_slices[ii+zz].edges[next_edge].right << std::endl;
							
							
							
							if ( (next_edge<0) || !(  (crit_slices[ii+zz].edges[next_edge].left!=current_bottom_ind) ||  crit_slices[ii+zz].edges[next_edge].right!=current_bottom_ind))  {
								continue;
							}
							
							if (crit_slices[ii+zz].edges[next_edge].left==current_bottom_ind) {
								current_bottom_ind = current_top_ind = crit_slices[ii+zz].edges[next_edge].right; // the upper value
							}
							else {
								current_bottom_ind = current_top_ind = crit_slices[ii+zz].edges[next_edge].left; // the lower value
							}
							
							// keep track of those edges we found.
							found_edges.insert(next_edge);
							
							possible_edges.erase(current_edge);
							possible_edges.erase(next_edge);
                            ////TODO:  check this for possibility of removing too many edges.  namely those we haven't actually found or whatever.
                            //							for (int ww=0; ww<candidate_counter; ww++) {
                            //								possible_edges.erase(candidates[ww]);
                            //							}
							
							// add the next edge to the set we can connect together.
							if (zz==0) {
								F.left.push_back(next_edge);
								F.num_left++;
							}
							else
							{
								F.right.push_back(next_edge);
								F.num_right++;
							}
							
							break;
						}
						else
						{
							possible_edges.erase(current_edge);
						}
						
					}//re: for qq
					
				} // re: while...
				
                
                
                
			}
			
			
			std::cout << "F.top " << F.top << std::endl;
			std::cout << "F.bottom " << F.bottom << std::endl;
			std::cout << "F.num_left " << F.num_left << std::endl;
			std::cout << "F.num_right " << F.num_right << std::endl;
			std::cout << "%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%\n\n";
			add_face(F);
		}
	}
	
	clear_vec_mp(blank_point);
	clear_mp(temp);
	clear_mp(temp2);
	clear_mp(numer);
	clear_mp(denom);
	
	
	clear_mp(found_proj_0);
	clear_mp(found_proj_1);
	return;
}









void surface_decomposition::print(boost::filesystem::path base)
{
	
	//	std::cout << "printing surface decomposition to folder " << base << std::endl;
	decomposition::print(base);
	
	
	
	boost::filesystem::path summaryname = base;
	summaryname /= "S.surf";
	FILE *OUT = safe_fopen_write(summaryname);
	fprintf(OUT,"%d %d %ld %ld\n\n", num_faces, num_edges, mid_slices.size(), crit_slices.size());
	// what more to print here?
	fclose(OUT);
	
	summaryname = base;
	summaryname /= "F.faces";
	print_faces(summaryname);
	
	
	boost::filesystem::path curve_location = base;
	curve_location /= "curve";
	
	for (int ii=0; ii<mid_slices.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_midslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		mid_slices[ii].print(specific_loc);
	}
	
	for (int ii=0; ii<crit_slices.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_critslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		crit_slices[ii].print(specific_loc);
	}
	
	boost::filesystem::path specific_loc = curve_location;
	specific_loc += "_crit";
	crit_curve.print(specific_loc);
	
	specific_loc = curve_location;
	specific_loc += "_sphere";
	sphere_curve.print(specific_loc);
}







void surface_decomposition::print_faces(boost::filesystem::path outputfile)
{
    //	std::cout << "printing faces to file " << outputfile << std::endl;
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d\n\n",num_faces);
	
	for(int ii=0;ii<num_faces;ii++){
		fprintf(OUT,"%d %d\n%d %d\n", faces[ii].midpt, faces[ii].index, faces[ii].top, faces[ii].bottom);
		fprintf(OUT,"%d %d\n",faces[ii].system_type_top,faces[ii].system_type_bottom);
		fprintf(OUT,"%ld\n",faces[ii].left.size());
		for (int jj=0; jj<faces[ii].left.size(); jj++) {
			fprintf(OUT,"%d ",faces[ii].left[jj]);
		}
		fprintf(OUT,"\n");
		fprintf(OUT,"%ld\n",faces[ii].right.size());
		for (int jj=0; jj<faces[ii].right.size(); jj++) {
			fprintf(OUT,"%d ",faces[ii].right[jj]);
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	fclose(OUT);
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
	fprintf(OUT,"END\n\n\n\n");
	
    
    for (int ii=0; ii<W.num_patches; ii++) {
		std::stringstream linname;
		linname << "patch_" << ii;
		write_vector_as_constants(W.patch_mp[ii], linname.str(), OUT);
		linname.clear();  linname.str("");
	}
    
	fclose(OUT);
	
	
	
	
}









void create_sphere_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
                          comp_mp sphere_radius,
                          vec_mp sphere_center,
                          const witness_set & W)
{
	// a bit of error checking
    
	
	if (W.variable_names.size()==0) {
		std::cout << color::red() << "trying to create a sphere system, but witness set does not have the variable names." << color::console_default() << std::endl;
		br_exit(-9846);
	}
	
	
	// got here, so ok to continue.
	
	
	int *declarations = NULL;
	
	partition_parse(&declarations, input_file, "func_input", "config", 0); // the 0 means not self conjugate.
	free(declarations);
	
	FILE *OUT = safe_fopen_write(output_file);
	
	
	
	// copy the original input file as parsed above.
	FILE *IN = safe_fopen_read("func_input");
	fprintf(OUT,"INPUT\n\n");
	copyfile(IN,OUT);
	fclose(IN);
	
	
	
	// now put in the sphere equations
	
	int rand_index = rand(); // what if there are multiple sphere equations???  put in this random number to be safe.
	
	fprintf(OUT,"function sphere_%d;\n",rand_index);
    
	fprintf(OUT,"sphere_%d = ",rand_index);
	for (int jj=1; jj<W.num_variables; jj++) { // start at 1 to omit the homogenizing variable
		fprintf(OUT," (%s-(", W.variable_names[jj].c_str());
		mpf_out_str(OUT,10,0,sphere_center->coord[jj-1].r);
		fprintf(OUT,"+I*");
		mpf_out_str(OUT,10,0,sphere_center->coord[jj-1].i);
		fprintf(OUT,"))^2");
		
		if (jj!=W.num_variables-1) {
			fprintf(OUT," + ");
		}
		
	}
    
	fprintf(OUT," - (");
	mpf_out_str(OUT,10,0,sphere_radius->r);
	fprintf(OUT,")^2 ");
    
	fprintf(OUT, ";\n\nEND;\n\n\n\n\n");
	
    for (int ii=0; ii<W.num_patches; ii++) {
		std::stringstream linname;
		linname << "patch_" << ii;
		write_vector_as_constants(W.patch_mp[ii], linname.str(), OUT);
		linname.clear();  linname.str("");
	}
    
    
	fclose(OUT);
	
	
	
}





