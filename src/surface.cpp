#include "surface.hpp"












void surface_decomposition::main(vertex_set & V,
                                 const witness_set & W_surf,
                                 vec_mp *pi,
                                 BR_configuration & program_options,
                                 solver_configuration & solve_options)
{
    
	
#ifdef functionentry_output
	std::cout << "surface_main" << std::endl;
#endif
	
#ifdef usetempfolders
	program_options.move_to_temp();
#endif
	
	// set some member information.
	this->input_filename = W_surf.input_filename;
	this->num_variables = W_surf.num_variables;
	this->component_num = W_surf.comp_num;
    
	add_projection(pi[0]); // add to *this
	add_projection(pi[1]); // add to *this
	
	
	
	
	
	
    
	//copy the patch over into this object
	for (int ii=0; ii<W_surf.num_patches; ii++)
		this->add_patch(W_surf.patch_mp[ii]);
	
	
	this->beginning_stuff( W_surf, program_options, solve_options);
	
	
    
    // get the witness points for the critical curve.
    witness_set W_critcurve;
	std::map< int, witness_set> higher_multiplicity_witness_sets;
	
	
    compute_critcurve_witness_set(W_critcurve,
								  higher_multiplicity_witness_sets,
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
    
	
	
	
	
	
	///////////////////////////////
	
	// the bounding sphere must be set before here, or this will eat it.
	
	
	std::map< std::pair<int,int>, witness_set > split_sets;
	deflate_and_split(split_sets,
					  higher_multiplicity_witness_sets,
					  program_options,
					  solve_options);
	
	
	
	witness_set W_singular_crit;
	compute_singular_crit(W_singular_crit,
						  split_sets,
						  V,
						  program_options,
						  solve_options);
	///////////////////////////////
	
	
	
	

	
    
    // merge together the critical points from both the critical curve and the sphere intersection curve.
    witness_set W_total_crit;
    
	W_total_crit.merge(W_critcurve_crit);

    W_total_crit.merge(W_sphere_crit);
    
	

    W_total_crit.merge(W_singular_crit);
    
	
    W_total_crit.input_filename = "total_crit___there-is-a-problem";
	W_total_crit.sort_for_unique(&solve_options.T);
    

    
    
    
	
	compute_critical_curve(W_critcurve, // all input.
                           W_total_crit,
                           V,
                           program_options,
                           solve_options);
	
		
	
	this->output_main(program_options.output_dir);
	V.print(program_options.output_dir/ "V.vertex");
	


	
    // actually perform the interslice on the bounding sphere curve.
	compute_bounding_sphere(W_intersection_sphere, // the witness points we will track from
                            W_total_crit,          // the critical points for the both sphere and critical curve.
                            V,                     // vertex set.  it goes almost everywhere.
                            program_options, solve_options); // configuration
    
	
	this->output_main(program_options.output_dir);
	V.print(program_options.output_dir/ "V.vertex");
	
	
	
	compute_singular_curves(W_total_crit,
							split_sets,
							V,
							program_options, solve_options);
	
	
	
    
    
    
    
    
    
	
	
	
	
	
	
	
	// compute the downstairs crit and midpoints for slicing
	vec_mp midpoints_downstairs, crit_downstairs;  std::vector< int > index_tracker;
	init_vec_mp(midpoints_downstairs,0); init_vec_mp(crit_downstairs,0);
    
    V.compute_downstairs_crit_midpts(W_total_crit, crit_downstairs, midpoints_downstairs, index_tracker, pi[0]);
	
	
	
	
	
	
	if (program_options.verbose_level>=0) {
        std::cout << color::green() << "the pi[0] projection values at which we will be slicing:\n\n" << color::console_default();
		print_point_to_screen_matlab(crit_downstairs,"crit_down");
		print_point_to_screen_matlab(midpoints_downstairs,"middown");
        
	}
	
    if (program_options.verbose_level>=2) {
		W_total_crit.print_to_screen();
	}
	
	
	
	
	
	
	// get the midpoint slices
	program_options.merge_edges=true;
	bool rerun_empty = false;
	compute_slices(W_surf, V,
                   midpoints_downstairs, this->mid_slices,
                   program_options, solve_options, rerun_empty,"mid");
	
	
	
	//incremental output
	this->output_main(program_options.output_dir);
	V.print(program_options.output_dir/ "V.vertex");
	
	// get the critical slices
	
	program_options.merge_edges=true;
	rerun_empty = true;
	compute_slices(W_surf, V,
                   crit_downstairs, this->crit_slices,
                   program_options, solve_options, rerun_empty,"crit");
	
    
	//incremental output
	this->output_main(program_options.output_dir);
	V.print(program_options.output_dir/ "V.vertex");
	
	//connect the dots - the final routine
	connect_the_dots(V, pi, program_options, solve_options);
	
	
	this->output_main(program_options.output_dir);
	V.print(program_options.output_dir/ "V.vertex");
	
	clear_vec_mp(crit_downstairs); clear_vec_mp(midpoints_downstairs);
    
#ifdef usetempfolders
	// this is likely extremely broken.  i'd like to use temp folders, but it's not there yet.
	program_options.move_to_called();
#endif
	
	
	
}








void surface_decomposition::beginning_stuff(const witness_set & W_surf,
                                            BR_configuration & program_options,
                                            solver_configuration & solve_options)
{
	
#ifdef functionentry_output
	std::cout << "surface::beginning_stuff" << std::endl;
#endif
	
	if (1) {
		// perform an isosingular deflation
		// note: you do not need witness_data to perform isosingular deflation
		if (program_options.verbose_level>=2)
			printf("performing isosingular deflation\n");
		
		
		program_options.input_deflated_filename = program_options.input_filename;
		
		std::stringstream converter;
		converter << "_dim_" << W_surf.dim << "_comp_" << W_surf.comp_num << "_deflated";
		program_options.input_deflated_filename += converter.str();
		
		
		W_surf.write_dehomogenized_coordinates("witness_points_dehomogenized"); // write the points to file
		int num_deflations, *deflation_sequence = NULL;
		
		isosingular_deflation(&num_deflations, &deflation_sequence, program_options,
							  program_options.input_filename,
							  "witness_points_dehomogenized",
							  program_options.input_deflated_filename, // the desired output name
							  program_options.max_deflations);
		free(deflation_sequence);
		
		
		
		converter.clear(); converter.str("");
	}
	else {
		program_options.input_deflated_filename = program_options.input_filename;
	}
	
	
	
	
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename);
	
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	if (0) {
		if (program_options.verbose_level>=2)
			printf("checking if component is self-conjugate\n");
		checkSelfConjugate(W_surf.pts_mp[0],program_options, W_surf.input_filename);
		
		//regenerate the various files, since we ran bertini since then and many files were deleted.
		parse_input_file(program_options.input_deflated_filename);
	}
	
	
	
	if (program_options.user_sphere) {
		read_sphere(program_options.bounding_sphere_filename);
	}
	
	
	randomizer->setup(W_surf.num_variables-1-this->dimension,solve_options.PPD.num_funcs);
	
	

	
	
	
	
	if (!verify_projection_ok(W_surf,
                              this->randomizer,
                              this->pi,
                              solve_options))
	{
		std::cout << "the projections being used appear to suffer rank deficiency with Jacobian matrix..." << std::endl;
		mypause();
	}
    
	
}














void surface_decomposition::compute_critcurve_witness_set(witness_set & W_critcurve,
														  std::map<int, witness_set> & higher_multiplicity_witness_sets,
                                                          const witness_set & W_surf,
                                                          BR_configuration & program_options,
                                                          solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_critcurve_witness_set" << std::endl;
#endif
	
	
	// find witness points on the critical curve.
	
	
	bool prev_quick_state = program_options.quick_run;
	program_options.quick_run = false;
	
	
	nullspace_config ns_config;
    
	solve_options.use_gamma_trick = 0;
	std::cout << color::bold("m") << "computing witness points for the critical curve" << color::console_default() << std::endl;
    
	solver_output solve_out;
	
	compute_crit_nullspace(solve_out,	// the returned value
                           W_surf,            // input the original witness set
                           this->randomizer,
                           this->pi,
                           2,  // dimension of ambient complex object
                           2,   //  target dimension to find
                           1,   // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	
	
	solve_out.get_nonsing_finite_multone(W_critcurve);
	solve_out.get_patches_linears(W_critcurve);
	solve_out.cp_names(W_critcurve);
	
	
	solve_out.get_multpos_full(higher_multiplicity_witness_sets);
	//now we have a map of multiplicities and witness sets.  for each point of each multiplicity, we need to perform iso. defl.  this will enable us to get ahold of the singular curves.
	
//	W_critcurve.print_to_screen();
//	//
//	for (auto iter = higher_multiplicity_witness_sets.begin(); iter!=higher_multiplicity_witness_sets.end(); iter++) {
//		std::cout << "found " << iter->second.num_points << " points of multiplicity " << iter->first << std::endl;
//
//		iter->second.print_to_screen();
//	}

	
	
	
	
	
	
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
    
	
    program_options.quick_run = prev_quick_state;
	
    return;
    
}


void surface_decomposition::compute_critcurve_critpts(witness_set & W_critcurve_crit,  // the computed value
                                                      const witness_set & W_surf,      // input witness set
                                                      witness_set & W_critcurve, // input witness set
                                                      BR_configuration & program_options, //as usual, goes everywhere.
                                                      solver_configuration & solve_options) // wtb: a way to make this global
{
#ifdef functionentry_output
	std::cout << "surface::compute_critcurve_critpts" << std::endl;
#endif
	
	std::cout << color::bold("m") << "\ncomputing critical points of the critical curve" <<  color::console_default() << std::endl;
	
	
	
	
	
	
	
	solve_options.use_gamma_trick = 0;
	
	nullspace_config ns_config; // this is set up in the nullspace call.
	
	solver_output solve_out;
	

	if (0) {
		
		int blabla;
		parse_input_file(W_critcurve.input_filename,&blabla);
		preproc_data_clear(&solve_options.PPD); // ugh this sucks
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		
		crit_curve.randomizer->setup(W_critcurve.num_variables - W_critcurve.num_patches - 1,solve_options.PPD.num_funcs);
		
		
		
		
		compute_crit_nullspace(solve_out, // the returned value
							   W_critcurve,            // input the original witness set
							   crit_curve.randomizer,
							   &this->pi[0],
							   1,  // dimension of ambient complex object
							   1,   //  target dimension to find
							   1,   // COdimension of the critical set to find.
							   program_options,
							   solve_options,
							   &ns_config);
	}


	//get crit points of the surface.
	int blabla;
	parse_input_file(W_surf.input_filename,&blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
	compute_crit_nullspace(solve_out, // the returned value
						   W_surf,            // input the original witness set
						   this->randomizer,
						   &this->pi[0],
						   2,  // dimension of ambient complex object
						   2,   //  target dimension to find
						   2,   // COdimension of the critical set to find.
						   program_options,
						   solve_options,
						   &ns_config);
	
	solve_out.get_noninfinite_w_mult_full(W_critcurve_crit);
	ns_config.clear();
	solve_out.reset();
	
	
	
	//now get those critpoints which lie on the crit curve (there may be some intersection with the above)
	
	W_critcurve.only_natural_vars();// trim off the synthetic variables for the input witness set for computing critical points
	parse_input_file(W_critcurve.input_filename,&blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	std::ifstream fin("deg.out");
	int max_degree = 0;
	int temp_degree;
	for (int ii=0; ii<solve_options.PPD.num_funcs; ii++) {
		fin >> temp_degree;
		if (temp_degree>max_degree) {
			max_degree = temp_degree;
		}
	}
	fin.close();
	
	temp_degree = solve_options.T.AMP_bound_on_degree;
	solve_options.T.AMP_bound_on_degree = max_degree;
	
	crit_curve.randomizer->setup(W_critcurve.num_variables - W_critcurve.num_patches - 1,solve_options.PPD.num_funcs);
	
	
	compute_crit_nullspace(solve_out, // the returned value
						   W_critcurve,            // input the original witness set
						   crit_curve.randomizer,
						   &this->pi[0],
						   1,  // dimension of ambient complex object
						   1,   //  target dimension to find
						   1,   // COdimension of the critical set to find.
						   program_options,
						   solve_options,
						   &ns_config);
	
	
	
	
	witness_set W_temp;
	solve_out.get_noninfinite_w_mult_full(W_temp);
	ns_config.clear();
	solve_out.reset();
	
	
	W_critcurve_crit.merge(W_temp);
	
	W_critcurve_crit.sort_for_real(&solve_options.T);
	W_critcurve_crit.sort_for_unique(&solve_options.T);
	
//	W_critcurve_crit.print_to_screen();
//	
//	
//	
//	
	
	
	
	if (have_sphere_radius) {
		W_critcurve_crit.sort_for_inside_sphere(sphere_radius, sphere_center);
	}
	else
	{
		this->compute_sphere_bounds(W_critcurve_crit); // sets the radius and center in this decomposition.  Must propagate to the constituent decompositions as well.   fortunately, i have a method for that!!!
	}
	
	
	crit_curve.copy_sphere_bounds(*this); // copy the bounds into the critcurve.
	
	
	
	std::cout << color::bold("m") << "intersecting critical curve with sphere" << color::console_default() << std::endl;
	
	W_critcurve_crit.input_filename = "input_critical_curve";
	
	
	witness_set W_sphere_intersection;
	// now get the sphere intersection critical points and ends of the interval
	crit_curve.get_additional_critpts(&W_sphere_intersection,  // the returned value
									  W_critcurve,       // all else here is input
									  program_options,
									  solve_options);
	

	W_sphere_intersection.sort_for_real(&solve_options.T);
	W_sphere_intersection.sort_for_unique(&solve_options.T);
	
	W_critcurve_crit.merge(W_sphere_intersection);
	
	solve_options.T.AMP_bound_on_degree = temp_degree;
    
    return;
}

void surface_decomposition::compute_critical_curve(const witness_set & W_critcurve,
                                                   const witness_set & W_critpts,
                                                   vertex_set & V,
                                                   BR_configuration & program_options,
                                                   solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_critical_curve" << std::endl;
#endif
	
	
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
                          &temp_proj,
                          program_options,
                          solve_options,
                          V);
	
	clear_vec_mp(temp_proj);
	
	crit_curve.add_projection(pi[0]);
	
	crit_curve.num_variables = W_critcurve.num_variables;
	
    std::cout << color::magenta() << "done decomposing critical curve" << color::console_default() << std::endl;
	return;
}













void surface_decomposition::deflate_and_split(std::map< std::pair<int,int>, witness_set > & split_sets,
											  std::map<int, witness_set > & higher_multiplicity_witness_sets,
											  BR_configuration & program_options,
											  solver_configuration & solve_options)
{
	
	for (auto iter = higher_multiplicity_witness_sets.begin(); iter!=higher_multiplicity_witness_sets.end(); ++iter) {
		std::cout << "multiplicity " << iter->first << std::endl;
	}
	
	for (auto iter = higher_multiplicity_witness_sets.begin(); iter!=higher_multiplicity_witness_sets.end(); ++iter) {
		std::cout << std::endl << color::magenta() << "splitting points for multiplicity " << iter->first << " singular curve" << color::console_default() << std::endl;
		
		int num_this_multiplicity = 0;
		
		
		
		
		
		
		
		witness_set active_set = iter->second; // seed the loop/  this sucks because it duplicates data
		active_set.only_first_vars(this->num_variables);
		
		witness_set W_only_one_witness_point;
		W_only_one_witness_point.copy_skeleton(active_set);
		
		
		while (active_set.has_points()) {
			
			witness_set W_reject; // these will be populated in the find matching call.
			std::pair<int,int> curr_index(iter->first,num_this_multiplicity);
			std::stringstream converter;
			
			converter << "_singcurve_mult_" << iter->first << "_" << num_this_multiplicity;
			
			boost::filesystem::path singcurve_filename = program_options.input_filename;
			singcurve_filename += converter.str(); converter.clear(); converter.str("");
			
			
			W_only_one_witness_point.reset_points();
			W_only_one_witness_point.add_point(active_set.pts_mp[0]); // exists by entrance condition
			W_only_one_witness_point.write_dehomogenized_coordinates("singular_witness_points_dehomogenized"); // write the points to file
			
			
			
			int num_deflations, *deflation_sequence = NULL;
			isosingular_deflation(&num_deflations, &deflation_sequence, program_options,
								  program_options.input_filename, // start from the beginning.
								  "singular_witness_points_dehomogenized",
								  singcurve_filename,
								  program_options.max_deflations);
			free(deflation_sequence);
			
			active_set.input_filename = singcurve_filename;
			active_set.dim = 1; //why again does a witness set need a dimension?
			
			int blabla;
			parse_input_file(singcurve_filename,&blabla);
			preproc_data_clear(&solve_options.PPD); // ugh this sucks
			parse_preproc_data("preproc_data", &solve_options.PPD);
			
			std::cout << "testing points for deflation validity" << std::endl;
			
			
			find_matching_singular_witness_points(split_sets[curr_index],
												  W_reject, //W_reject contains the points of this multiplicity which DO NOT satisfy this deflation.  they must be deflated again.
												  active_set,//input witness set
												  solve_options);
			
			split_sets[curr_index].input_filename = singcurve_filename;
			
			
			if (W_reject.has_points()) {
				std::cout << "found that current singular witness set had " << W_reject.num_points << " non-deflated points" << std::endl;
			}
			
			//TODO: this is an ideal place for a swap operator.
			active_set = W_reject;
			num_this_multiplicity++;
			
		}
		
		
		
		
		
		
	}
	
	for (auto iter = split_sets.begin(); iter!=split_sets.end(); ++iter) {
		std::cout << "multiplicity " << iter->first.first << " " << iter->first.second << std::endl;
	}
	
}



void surface_decomposition::compute_singular_crit(witness_set & W_singular_crit,
												  const std::map<std::pair<int,int>, witness_set> & split_sets,
												  vertex_set & V,
												  BR_configuration & program_options,
												  solver_configuration & solve_options)
{
	
	
	
	W_singular_crit.num_variables = this->num_variables;
	W_singular_crit.num_natural_vars = this->num_variables;
	W_singular_crit.copy_patches(*this);
	for (auto iter = split_sets.begin(); iter!=split_sets.end(); ++iter) {
		
		std::cout << std::endl << color::magenta() << "getting critical points for multiplicity " << iter->first.first << " " << iter->first.second << " singular curve" << color::console_default() << std::endl;
		
		
		
		
		int blabla;
		parse_input_file(iter->second.input_filename,&blabla);
		preproc_data_clear(&solve_options.PPD); // ugh this sucks
		parse_preproc_data("preproc_data", &solve_options.PPD);
		
		std::ifstream fin("deg.out");
		int max_degree = 0;
		int temp_degree;
		for (int ii=0; ii<solve_options.PPD.num_funcs; ii++) {
			fin >> temp_degree;
			if (temp_degree>max_degree) {
				max_degree = temp_degree;
			}
		}
		fin.close();
		
		temp_degree = solve_options.T.AMP_bound_on_degree;
		solve_options.T.AMP_bound_on_degree = max_degree;
		
		
		
		
		singular_curves[iter->first].randomizer->setup(iter->second.num_variables-iter->second.num_patches-1, solve_options.PPD.num_funcs);

		
		
		nullspace_config ns_config;
		solver_output solve_out;
		
		
		compute_crit_nullspace(solve_out,                   // the returned value
							   iter->second,            // input the witness set with linears
							   singular_curves[iter->first].randomizer,
							   &(this->pi[0]),
							   1,                                // dimension of ambient complex object
							   1,                                //  target dimension to find
							   1,                                // COdimension of the critical set to find.
							   program_options,
							   solve_options,
							   &ns_config);
		
		
		solve_options.T.AMP_bound_on_degree = temp_degree;
		
		witness_set W_this_round;
		solve_out.get_noninfinite_w_mult_full(W_this_round);
		
		ns_config.clear();
		
//		W_this_round.only_first_vars(this->num_variables);
		W_this_round.sort_for_unique(&solve_options.T); // this could be made to be unnecessary, after rewriting a bit of solverout
		W_this_round.sort_for_real(&solve_options.T);
		W_this_round.sort_for_inside_sphere(sphere_radius,sphere_center);
		W_this_round.input_filename = iter->second.input_filename;
		
		singular_curves[iter->first].add_witness_set(W_this_round,CRITICAL,V); // creates the curve decomposition for this multiplicity
		
		
		
		W_singular_crit.merge(W_this_round);
		
	}
	
}





void surface_decomposition::compute_singular_curves(const witness_set & W_total_crit,
													const std::map< std::pair<int,int>, witness_set> & split_sets,
													vertex_set & V,
													BR_configuration & program_options,
													solver_configuration & solve_options)
{
	program_options.merge_edges=false;
	
	for (auto iter = split_sets.begin(); iter!=split_sets.end(); ++iter) {
		
		singular_curves[iter->first].input_filename = iter->second.input_filename;
		singular_curves[iter->first].copy_sphere_bounds(*this);
		
		
		std::cout << "getting sphere intersection with singular curve " << iter->first.first << " " << iter->first.second << std::endl;
		witness_set W_sphere_intersection;
		// now get the sphere intersection critical points and ends of the interval
		singular_curves[iter->first].get_additional_critpts(&W_sphere_intersection,  // the returned value
															iter->second,       // all else here is input
															program_options,
															solve_options);
		
		W_sphere_intersection.input_filename = iter->second.input_filename;
		
		
//		W_sphere_intersection.only_first_vars(this->num_variables); // throw out the extra variables.
		W_sphere_intersection.sort_for_real(&solve_options.T);
		W_sphere_intersection.sort_for_unique(&solve_options.T);
		
		singular_curves[iter->first].add_witness_set(W_sphere_intersection,CRITICAL,V); // creates the curve decomposition for this multiplicity
		
		W_sphere_intersection.merge(W_total_crit);
		W_sphere_intersection.input_filename = "should_have_already_been_added_elsewhere";
		
		
		
		
		
		
		std::cout << "interslicing singular curve " << std::endl;
		
		// we already know the component is self-conjugate (by entry condition), so we are free to call this function
		singular_curves[iter->first].interslice(iter->second,
												W_sphere_intersection,
												&pi[0],
												program_options,
												solve_options,
												V);
		singular_curves[iter->first].add_projection(this->pi[0]);
		
		num_singular_curves++;
		
		
		
	}
	
}



// will compute a randomizer matrix since you don't provide one. must have current PPD in solve_options for this to work correctly
int find_matching_singular_witness_points(witness_set & W_match,
										  witness_set & W_reject,
										  const witness_set & W,
										  solver_configuration & solve_options)
{
	
	if (W.has_no_points()) {
		std::cout << color::red() << "input witness set for find_matching_  has NO points, but hypothetically it does..." << color::console_default() << std::endl;
		return TOLERABLE_FAILURE;
	}
	
	
	// assumes input file for W is already parsed.
	
	
	
	prog_t SLP;
	setupProg(&SLP, solve_options.T.Precision, 2);
	
	
	comp_mp zerotime; init_mp(zerotime);
	set_zero_mp(zerotime);
	
	
	
	eval_struct_mp ED; init_eval_struct_mp(ED, 0, 0, 0);
	
	tracker_config_t * T = &solve_options.T;
	double tol = MAX(T->final_tol_times_mult, T->sing_val_zero_tol);
	
	mat_mp U, E, V; init_mat_mp(U, 0, 0); init_mat_mp(E, 0, 0); init_mat_mp(V, 0, 0);
	
	
	
	
	evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, W.pts_mp[0], zerotime, &SLP);
	int hypothesis_corank = svd_jacobi_mp_prec(U, E, V, ED.Jv, tol, T->Precision); // this wraps around svd_jacobi_mp.
	
	
	
	std::vector< bool > validity_flag;
	validity_flag.push_back(true);
	for (int zz = 1;zz<W.num_points;++zz) // by hypothesis, the first (0th) point satisfies the deflation
	{
		
		evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, W.pts_mp[zz], zerotime, &SLP);
		
		// first, check that the point satifies the system.
		if (d_vec_abs_mp(ED.funcVals)>T->final_tol_times_mult) { // is this the correct measure of the vector to compare?
			validity_flag.push_back(false);
			continue;
		}
		
		
		// now, check the rank and make sure is same as first (0th) point.
		int corank = svd_jacobi_mp_prec(U, E, V, ED.Jv, tol, T->Precision); // this wraps around svd_jacobi_mp.
		
		if (corank != hypothesis_corank) {
			
			validity_flag.push_back(false);
			continue;
		}
		else{
			validity_flag.push_back(true);
		}
	}
	
	
	for (int zz=0; zz<W.num_points; zz++) {
		if (validity_flag[zz]==true) { // trivially true for first point -- it generated the deflation!
			W_match.add_point(W.pts_mp[zz]);
		}
		else
		{
			std::cout << "adding reject point" << std::endl;
			W_reject.add_point(W.pts_mp[zz]);
		}
	}
	
	
	
	W_match.copy_skeleton(W);
	W_reject.copy_skeleton(W);
	
	W_match.copy_linears(W);
	W_reject.copy_linears(W);
	
	W_match.copy_patches(W);
	W_reject.copy_patches(W);
	
	
	
	clear_mat_mp(U); clear_mat_mp(E); clear_mat_mp(V);
	clear_mp(zerotime);
	clear_eval_struct_mp(ED); clearProg(&SLP, solve_options.T.MPType, 1);
	
	
	return SUCCESSFUL;
}







void surface_decomposition::compute_sphere_witness_set(const witness_set & W_surf,
                                                       witness_set & W_intersection_sphere,
                                                       BR_configuration & program_options,
                                                       solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_sphere_witness_set" << std::endl;
#endif
    
	
    
    
	//build up the start system
	solve_options.robust = true;
	solve_options.use_gamma_trick = 0;
	
	int blabla;
	
	
	
	parse_input_file(W_surf.input_filename, &blabla);
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
    
    
	if (!this->randomizer->is_ready()) {
		std::cout << "randomizer not ready" << std::endl;
		br_exit(19355);
	}
	
	
	sphere_curve.randomizer->setup(W_surf.num_variables-W_surf.num_patches-W_surf.num_linears, solve_options.PPD.num_funcs);
	
//	make_randomization_matrix_based_on_degrees(randomizer_matrix, sphere_curve.randomized_degrees,
//                                               W_surf.num_variables-W_surf.num_patches-W_surf.num_linears,
//                                               solve_options.PPD.num_funcs);
	// what do do about this??????????????????????????????????????
//	sphere_curve.randomized_degrees.push_back(2);
	
	
	multilin_config ml_config(solve_options,this->randomizer); // copies in the randomizer matrix and sets up the SLP & globals.
	
    
	vec_mp *multilin_linears = (vec_mp *) br_malloc(2*sizeof(vec_mp));
    
	init_vec_mp2(multilin_linears[0],W_surf.num_variables,solve_options.T.AMP_max_prec);
	multilin_linears[0]->size = W_surf.num_variables;
    
	init_vec_mp2(multilin_linears[1],0,solve_options.T.AMP_max_prec);
	vec_cp_mp(multilin_linears[1],W_surf.L_mp[1]);
	
    
    
	witness_set W_sphere;
	sphere_config sp_config(sphere_curve.randomizer);
	
	
	for (int ii=0; ii<2; ii++) {
		
		for (int jj=0; jj<W_surf.num_variables; jj++) {
			get_comp_rand_mp(&multilin_linears[0]->coord[jj]);
		}
		
		vec_cp_mp(sp_config.starting_linear[ii], multilin_linears[0]); // copy in the first multilinear to the new witness set we are computing.
		
		witness_set W_temp;
		
		
		
		solve_options.allow_singular = 0;
		solve_options.complete_witness_set = 0;
		solve_options.allow_multiplicity = 0;
		solve_options.allow_unsuccess = 0;
		
		
		solver_output fillme;
		multilin_solver_master_entry_point(W_surf,         // witness_set
                                           fillme, // the new data is put here!
                                           multilin_linears,
                                           ml_config,
                                           solve_options);
		
		//get stuff from fillme into W_temp, or whatever
		fillme.get_noninfinite_w_mult(W_temp);
		
		W_sphere.merge(W_temp);
		
	}
	
	clear_vec_mp(multilin_linears[0]);
	clear_vec_mp(multilin_linears[1]);
	free(multilin_linears);
	
	
	
	W_sphere.add_linear(W_surf.L_mp[1]);
	W_sphere.copy_patches(W_surf); //.patch_mp[0]
	
	
	
	// need to actually move to the sphere system now.
	
	
	sp_config.set_memory(solve_options);
	sp_config.set_center(this->sphere_center);
	sp_config.set_radius(this->sphere_radius);
	
	
	
	solve_options.allow_singular = 0;
	solve_options.complete_witness_set = 1;
	solve_options.allow_multiplicity = 0;
	solve_options.allow_unsuccess = 0;
	
	
	
	
	
	solver_output fillme;
	sphere_solver_master_entry_point(W_sphere,
									 fillme,
                                     sp_config,
                                     solve_options);
	
	fillme.get_noninfinite_w_mult_full(W_intersection_sphere);
	
	W_intersection_sphere.input_filename = "input_surf_sphere";
	
}




void surface_decomposition::compute_sphere_crit(const witness_set & W_intersection_sphere,
                                                witness_set & W_sphere_crit,
                                                BR_configuration & program_options,
                                                solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_sphere_crit" << std::endl;
#endif
	
	
	int blabla;
	parse_input_file("input_surf_sphere", &blabla); // having already been written to disk
	preproc_data_clear(&solve_options.PPD); // ugh this sucks
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	
    
    
    
	sphere_curve.randomizer->setup(W_intersection_sphere.num_variables-W_intersection_sphere.num_patches-1, solve_options.PPD.num_funcs);
    
	
	
	solve_options.use_gamma_trick = 0;
	
	nullspace_config ns_config;
	
	std::cout << "computing critical points of sphere curve" << std::endl;
	
	solver_output solve_out;
	
	compute_crit_nullspace(solve_out,                   // the returned value
                           W_intersection_sphere,            // input the witness set with linears
                           sphere_curve.randomizer,
                           &(this->pi[0]),
                           1,                                // dimension of ambient complex object
                           1,                                //  target dimension to find
                           1,                                // COdimension of the critical set to find.
                           program_options,
                           solve_options,
                           &ns_config);
	
	solve_out.get_noninfinite_w_mult_full(W_sphere_crit);
	
	ns_config.clear();
	
	W_sphere_crit.sort_for_unique(&solve_options.T);
	W_sphere_crit.sort_for_real(&solve_options.T);
    
	W_sphere_crit.input_filename = "input_surf_sphere";
	
	return;
}



void surface_decomposition::compute_bounding_sphere(const witness_set & W_intersection_sphere,
                                                    const witness_set & W_crit,
                                                    vertex_set & V,
                                                    BR_configuration & program_options,
                                                    solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::compute_bounding_sphere" << std::endl;
#endif
    
	program_options.merge_edges=false;
	
	this->sphere_curve.compute_sphere_bounds(W_crit); // inflate around this thing because the interslice method requires having bounds in place.
	
	std::cout << color::magenta() << "entering interslice for sphere" << color::console_default() << std::endl;
	// then feed it to the interslice algorithm
	this->sphere_curve.interslice(W_intersection_sphere, // the witness set with a single linear and some patches.
                                  W_crit, // the critical points
                                  &pi[0],
                                  program_options,
                                  solve_options,
                                  V);
	
	sphere_curve.input_filename = W_intersection_sphere.input_filename;
	this->sphere_curve.add_projection(pi[0]);
	this->sphere_curve.num_variables = num_variables;
	
	
	std::cout << color::magenta() << "done decomposing sphere curve" << color::console_default() << std::endl;
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
#ifdef functionentry_output
	std::cout << "surface::compute_slices" << std::endl;
#endif
	
	
	
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
	
	
	
	multilin_config ml_config(solve_options); // copies in the randomizer matrix and sets up the SLP & globals.
	
	for (int ii=0; ii<projection_values_downstairs->size; ii++){
		
		
		std::cout << color::magenta() << "decomposing the " << ii << "th " << kindofslice << " slice" << color::console_default() << std::endl;
		print_comp_matlab(&projection_values_downstairs->coord[ii], "target_proj");
		
		solve_options.backup_tracker_config();
		
		int iterations=0;
		
		witness_set slice_witness_set; // deliberately scoped variable
		
		
		
		std::stringstream converter;
		converter << ii;
		
		int blabla;
		parse_input_file(W_surf.input_filename, &blabla);
		
		preproc_data_clear(&solve_options.PPD); // ugh this sucks
		parse_preproc_data("preproc_data", &solve_options.PPD);
		

		ml_config.set_randomizer(this->randomizer);
		neg_mp(&multilin_linears[0]->coord[0], &projection_values_downstairs->coord[ii]);
		
		
		
		solve_options.robust = true;
		
		
		solver_output fillme;
		multilin_solver_master_entry_point(W_surf,         // witness_set
										   fillme, // the new data is put here!
										   multilin_linears,
										   ml_config,
										   solve_options);
		
		fillme.get_noninfinite_w_mult(slice_witness_set);
		
		fillme.reset();
		
		boost::filesystem::path slicename = W_surf.input_filename;
		slicename += "_"; slicename += kindofslice; slicename += "slice_"; slicename += converter.str();
		create_sliced_system(W_surf.input_filename, slicename, &multilin_linears[0], 1, W_surf);
		
		
		parse_input_file(slicename, &blabla);
		preproc_data_clear(&solve_options.PPD);
		parse_preproc_data("preproc_data", &solve_options.PPD);
		

		
		
		slice_witness_set.num_variables = W_surf.num_variables;
		slice_witness_set.input_filename = slicename;
		slice_witness_set.add_linear(multilin_linears[1]);
		slice_witness_set.add_patch(W_surf.patch_mp[0]);
		slice_witness_set.variable_names = W_surf.variable_names;
		slice_witness_set.dim = 1;
		
		
		solve_options.complete_witness_set = 1;
		
		
		slices[ii].input_filename = slicename;
		
		slices[ii].copy_sphere_bounds(*this);
		
		// we already know the component is self-conjugate (by entry condition), so we are free to call this function
		// the memory for the multilin system will get erased in this call...
		bool prev_quick_state = program_options.quick_run;
//		program_options.quick_run = false;
		
		slices[ii].computeCurveSelfConj(slice_witness_set,
									   &pi[1],
									   V,
									   program_options, solve_options);
		
		program_options.quick_run = prev_quick_state;
		
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
		
		
		
		slice_witness_set.reset();
			

		solve_options.reset_tracker_config();
		
		slices[ii].add_projection(pi[1]);
		slices[ii].num_variables = num_variables;
		
		
		
        // does it matter speedwise whether i do this before or after the copy immediately above?  i think the answer is no.
        V.assert_projection_value(slices[ii].all_edge_indices(), &projection_values_downstairs->coord[ii], 0); // the 0 is an index into the number of projections.
		
		this->output_main(program_options.output_dir);
		V.print(program_options.output_dir/ "V.vertex");
		
		std::cout << color::magenta() << "DONE decomposing the " << ii << "th " << kindofslice << " slice" << color::console_default() << std::endl;
		
	} // re: for loop
	
	clear_mp(rand_perterb);  clear_mp(h);
	clear_vec_mp(multilin_linears[0]);
	clear_vec_mp(multilin_linears[1]);
	free(multilin_linears);
	return;
}




















void surface_decomposition::connect_the_dots(vertex_set & V,
                                             vec_mp *pi,
                                             BR_configuration & program_options,
                                             solver_configuration & solve_options)
{
#ifdef functionentry_output
	std::cout << "surface::connect_the_dots" << std::endl;
#endif
	
	std::cout << color::bold("m") << "***\n\nCONNECT THE DOTS\n\n***" << color::console_default() << std::endl;
    
	
	
	
	
	
	
	
	
	midpoint_config md_config;
	md_config.setup(*this, solve_options); // yep, pass 'this' object into another call. brilliant.
	
	
	
	
	if (solve_options.use_parallel()) {
		master_connect(V, md_config, solve_options, program_options);
	}
	else
	{
		serial_connect(V, md_config, solve_options, program_options);
	}
	
    
    
	
	
	
	
	
	
	return;
}


void surface_decomposition::serial_connect(vertex_set & V, midpoint_config & md_config, solver_configuration & solve_options, BR_configuration & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::serial_connect" << std::endl;
#endif
	
    
	this->output_main(program_options.output_dir);
	V.print(program_options.output_dir/ "V.vertex");
	
	
	for (unsigned int ii=0; ii!=mid_slices.size(); ii++) { // each edge of each midslice will become a face.  degenerate edge => degenerate face.
		
		
		for (int jj=0; jj<mid_slices[ii].num_edges; jj++) {
			
			
			
			//make face
			
			face F = make_face(ii,jj, V, md_config, solve_options, program_options);
			
			//			std::cout << "F.top " << F.top << std::endl;
			//			std::cout << "F.bottom " << F.bottom << std::endl;
			//			std::cout << "F.num_left " << F.num_left << std::endl;
			//			std::cout << "F.num_right " << F.num_right << std::endl;
			
			
			if (!F.is_degenerate())
			{
				add_face(F);
				
//					boost::filesystem::path::rename(program_options.output_dir / "F.faces");
//					this->print_faces(program_options.output_dir + "_bak" / "F.faces");
				this->print_faces(program_options.output_dir / "F.faces");
//					this->output_main(program_options.output_dir);
//					V.print(program_options.output_dir/ "V.vertex");
			}
			

		}
	}
	
	return;
}



void surface_decomposition::master_connect(vertex_set & V, midpoint_config & md_config, solver_configuration & solve_options, BR_configuration & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::master_connect" << std::endl;
#endif
	
	
	boost::timer::auto_cpu_timer t;
    
	MPI_Status statty_mc_gatty;
	
	solve_options.call_for_help(MIDPOINT_SOLVER); // sets available workers, too
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	md_config.bcast_send(solve_options);
	
	
	
	if (V.num_vertices > 1e5) {
		std::cout << color::red() << "attempting to transmit over 1e5 points to all workers..." << color::console_default() << std::endl;
	}
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	//seed the workers
	for (int ii=1; ii<solve_options.numprocs; ii++) {
		this->send(ii, solve_options);
		V.send(ii, solve_options);
	}
	
	
	this->output_main(program_options.output_dir);
	V.print(program_options.output_dir/ "V.vertex");
	
	// this loop is semi-self-seeding
	for (unsigned int ii=0; ii!=mid_slices.size(); ii++) { // each edge of each midslice will become a face.  degenerate edge => degenerate face.
		
		
		for (int jj=0; jj<mid_slices[ii].num_edges; jj++) {
			
			if (mid_slices[ii].edges[jj].is_degenerate()) {
				continue;
			}
			
			int next_worker = solve_options.activate_next_worker();
			
			int send_num_faces = 1;// num_faces doubles as the keep_going signal.  if 0, the worker halts.
			MPI_Send(&send_num_faces, 1, MPI_INT, next_worker, NUMPACKETS, MPI_COMM_WORLD);
			
			
			master_face_requester(ii,jj, next_worker, solve_options);
			
			
			
			if (solve_options.have_available()) {
				
				continue;
			}
			else
			{
				//perform blocking receive of the face.
				int recv_num_faces;
				MPI_Recv(&recv_num_faces, 1, MPI_INT, MPI_ANY_SOURCE, NUMPACKETS, MPI_COMM_WORLD, &statty_mc_gatty);
				bool added_face = false;
				for (int qq = 0; qq<recv_num_faces; qq++) {
					face F;
					F.receive(statty_mc_gatty.MPI_SOURCE, solve_options);
					if (!F.is_degenerate()) {
						add_face(F);
						added_face = true;
					}
					
				}
				
				solve_options.deactivate(statty_mc_gatty.MPI_SOURCE);
				
				
				
				//TODO: this needs to be inside of a better protector, in the sense that we shouldn't do it after EVERY face, but rather at thoughtful times.
				if (added_face) { // only the faces file really needs to be updated.  the rest stay the same.  this is very wasteful.
//					boost::filesystem::path::rename(program_options.output_dir / "F.faces");
//					this->print_faces(program_options.output_dir + "_bak" / "F.faces");
					this->print_faces(program_options.output_dir / "F.faces");
//					this->output_main(program_options.output_dir);
//					V.print(program_options.output_dir/ "V.vertex");
				}
				
				
			}
		}
	}
	
	
	//cleanup
	
	
	while (solve_options.have_active()) {
		int recv_num_faces;
		MPI_Recv(&recv_num_faces, 1, MPI_INT, MPI_ANY_SOURCE, NUMPACKETS, MPI_COMM_WORLD, &statty_mc_gatty);
		for (int ii=0; ii<recv_num_faces; ii++) {
			face F;
			F.receive(statty_mc_gatty.MPI_SOURCE, solve_options);
			add_face(F);
			solve_options.deactivate(statty_mc_gatty.MPI_SOURCE);
		}
	}
	
	solve_options.send_all_available(0);
	
	
	return;
}





void surface_decomposition::worker_connect(solver_configuration & solve_options, BR_configuration & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::worker_connect" << std::endl;
#endif
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	MPI_Barrier(MPI_COMM_WORLD);
	solve_options.robust = true;
	
//	std::cout << "worker getting md_congid" << std::endl;
	midpoint_config md_config;
	md_config.bcast_receive(solve_options);
	//receive the md_config from the master.  it holds the three SLP's, as well as everything needed to run the system except:
	//	• system types
	//	• patches
	//	• starting point
	//which are all updated from another call later.
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
//	std::cout << "worker getting surface" << std::endl;
	this->receive(solve_options.head(), solve_options);
	
	
	
	vertex_set V;
	V.receive(solve_options.head(), solve_options);
	
	
	
	
	
	while (1) {
		
		// get the continue or discontinue signal.
		int keep_going;
		MPI_Status statty_mc_gatty;
		MPI_Recv(&keep_going, 1, MPI_INT, solve_options.head(), NUMPACKETS, MPI_COMM_WORLD, &statty_mc_gatty);
		if (keep_going==0) {
			break;
		}
		
		
		
		int ii, jj;
		// get the indices of the face to make.
		worker_face_requester(ii,jj,solve_options);
		
		
		//make the face
		face F = make_face(ii, jj, V, md_config, solve_options, program_options);
		
		
		//send the face back to master.
		int send_num_faces = 1;
		MPI_Send(&send_num_faces, 1, MPI_INT, solve_options.head(), NUMPACKETS, MPI_COMM_WORLD);
		
		F.send(solve_options.head(), solve_options);
	}
	
	
	
	return;
}




void surface_decomposition::master_face_requester(int ii, int jj, int next_worker, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "surface::master_face_requester" << std::endl;
#endif
	
	
	
	
	int * buffer = new int[2];
	buffer[0] = ii;
	buffer[1] = jj;
	MPI_Send(&buffer[0], 2, MPI_INT, next_worker, DATA_TRANSMISSION, MPI_COMM_WORLD);
	delete [] buffer;
}


void surface_decomposition::worker_face_requester(int & ii, int & jj, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "surface::worker_face_requester" << std::endl;
#endif
	
	
	
	
	int * buffer = new int[2];
	MPI_Status statty_mc_gatty;
	
	MPI_Recv(&buffer[0], 2, MPI_INT, mpi_config.head(), DATA_TRANSMISSION, MPI_COMM_WORLD, &statty_mc_gatty);
	ii = buffer[0];
	jj = buffer[1];
	
	delete [] buffer;
	return;
}

face surface_decomposition::make_face(int ii, int jj, vertex_set & V,
									  midpoint_config & md_config,
									  solver_configuration & solve_options, BR_configuration & program_options)
{
#ifdef functionentry_output
	std::cout << "surface::make_face" << std::endl;
#endif
	
	
	
	
	// assert some solver options
	solve_options.use_gamma_trick = 0;
    solve_options.robust = true;
	
	//create the face
	face F;
	
	std::cout << color::magenta() << "\n\n\n*****************************\nmidslice " << ii << " / " << this->mid_slices.size() <<  ", edge " << jj << " / " << mid_slices[ii].edges.size() << color::console_default() << "\n";
	
	if (mid_slices[ii].edges[jj].is_degenerate()) {
//		std::cout << "no computation necessary -- midslice edge is degenerate" << std::endl;
		return F;
	}
	
	
	
	comp_mp temp, temp2, temp3; init_mp2(temp,1024); init_mp2(temp2,1024); init_mp2(temp3,1024);
	comp_mp numer, denom; init_mp2(numer,1024); init_mp2(denom,1024);
	comp_mp proj_top, proj_bottom, proj_mid; init_mp2(proj_mid,1024); init_mp2(proj_bottom,1024); init_mp2(proj_top,1024);
	vec_mp found_point; init_vec_mp2(found_point, this->num_variables,1024); found_point->size = this->num_variables;
	
	
	
	
	
	
	
	F.crit_slice_index = ii; // the index of which midslice this face came from.
	F.midpt = mid_slices[ii].edges[jj].midpt; // index the point
	
	std::cout << color::brown() << "current midpoint: " <<  mid_slices[ii].edges[jj].midpt  << " " << color::console_default() << "\n";
	
	
	
	// perform a search to find the top and bottom edges in the appropriate curve.
	
	
	// get the type of system for the top and bottom edges.  this is determined by reading the system name for the midpoints.
	// info on the files from which the points came from.
	int bottom_input_index = V.vertices[mid_slices[ii].edges[jj].left].input_filename_index;
	int top_input_index = V.vertices[mid_slices[ii].edges[jj].right].input_filename_index;
	
	bool bail_out = false;
	if (md_config.systems.find(V.filenames[bottom_input_index].string())==md_config.systems.end()) {
		std::cout << "bottom system is " << V.filenames[bottom_input_index] << ", which is not in md_config" << std::endl;
		bail_out = true;
	}
	if (md_config.systems.find(V.filenames[top_input_index].string())==md_config.systems.end()) {
		std::cout << "top system is " << V.filenames[top_input_index] << ", which is not in md_config" << std::endl;
		bail_out = true;
	}
	
	
	md_config.system_name_bottom = V.filenames[bottom_input_index].filename().string();
	md_config.system_name_top = V.filenames[top_input_index].filename().string();
	md_config.system_name_mid = this->input_filename.filename().string();
	
	

	std::cout << md_config.system_name_top << " " << md_config.system_name_bottom << std::endl;
	
	
	curve_decomposition * top_curve = curve_with_name(md_config.system_name_top);
	curve_decomposition * bottom_curve = curve_with_name(md_config.system_name_bottom);
	
	//man, i hate checking for null...
	
	if (top_curve==NULL) {
		std::cout << "did not find matching top curve: " << md_config.system_name_top << std::endl;
		bail_out = true;
	}
	
	if (bottom_curve==NULL) {
		std::cout << "did not find matching bottom curve: " << md_config.system_name_bottom << std::endl;
		bail_out = true;
	}
	
	
	
	if (crit_slices[ii].num_edges == 0) {
		std::cout << "critslice " << ii << " has no edges!" << std::endl;
		bail_out = true;
	}
	if (crit_slices[ii+1].num_edges == 0) {
		std::cout << "critslice " << ii+1 << " has no edges!" << std::endl;
		bail_out = true;
	}
	
	
	// get the bottom and top edges for this face.
	
	int num_bottom_vars = md_config.num_bottom_vars();
	int num_top_vars = md_config.num_top_vars();
	
	F.bottom = bottom_curve->edge_w_midpt(mid_slices[ii].edges[jj].left); // index the *edge*
	F.system_name_bottom = md_config.system_name_bottom;
	
	F.top= top_curve->edge_w_midpt(mid_slices[ii].edges[jj].right); // index the *edge*
	F.system_name_top = md_config.system_name_top;
	
	if (F.bottom < 0) { // this would happen because the bottom point was of type PROBLEMATIC
		std::cout << "F.bottom is set to negative" << std::endl;
		bail_out = true;
	}
	
	if (F.top < 0) { // this would happen because the bottom point was of type PROBLEMATIC
		std::cout << "F.top is set to negative" << std::endl;
		bail_out = true;
	}
		
	
	if (num_bottom_vars==0) {
		std::cout << "0 bottom variables" << std::endl;
		bail_out = true;
	}
	
	if (num_top_vars==0) {
		std::cout << "0 top variables" << std::endl;
		bail_out = true;
	}
	
	
	
	if (bail_out) {
		std::cout << color::red() << "bailing out " << ii << " " << jj << "." << std::endl;
		
		std::cout << "tracking from these point indices:" << std::endl;
		std::cout <<  mid_slices[ii].edges[jj].left  << " " << mid_slices[ii].edges[jj].midpt  << " "  << mid_slices[ii].edges[jj].right << color::console_default() << std::endl;
		
		return F;
	}
	
	
	//copy in the start point as three points concatenated.
	witness_set W_midtrack;
	vec_mp blank_point;  init_vec_mp2(blank_point, 0,1024);
	W_midtrack.add_point(blank_point);
	clear_vec_mp(blank_point);
	
	W_midtrack.num_variables = this->num_variables + num_bottom_vars + num_top_vars;
	W_midtrack.num_natural_vars = this->num_variables;
	change_size_vec_mp(W_midtrack.pts_mp[0], W_midtrack.num_variables); W_midtrack.pts_mp[0]->size = W_midtrack.num_variables; // destructive resize
	
	
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
		W_midtrack.add_patch(this->patch_mp[qq]);
	}
	
	for (int qq = 0; qq<bottom_curve->num_patches; qq++) {
		W_midtrack.add_patch(bottom_curve->patch_mp[qq]);
	}
	
	for (int qq = 0; qq<top_curve->num_patches; qq++) {
		W_midtrack.add_patch(top_curve->patch_mp[qq]);
	}
	
	

	
	
	
	
	
	
	// make u, v target values.
	
	
//	print_comp_matlab(&V.vertices[ crit_slices[ii].edges[0].midpt ].projection_values->coord[0],"a");
//	print_comp_matlab(&V.vertices[ crit_slices[ii+1].edges[0].midpt ].projection_values->coord[0],"b");
	
	set_mp(md_config.crit_val_left,   &V.vertices[ crit_slices[ii].edges[0].midpt ].projection_values->coord[0]);
	set_mp(md_config.crit_val_right,  &V.vertices[ crit_slices[ii+1].edges[0].midpt ].projection_values->coord[0]);
	
	
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
		int final_top_ind = -1, final_bottom_ind = -2; // indexes in V of the bottom and top points of the left or right edge.
		
		
		if (zz==0) {
			
			set_zero_mp(md_config.u_target);
			final_bottom_ind = bottom_curve->edges[F.bottom].left;
			final_top_ind = top_curve->edges[F.top].left;
			
		}
		else{ // zz==1, and going right
			
			set_one_mp(md_config.u_target);
			final_bottom_ind = bottom_curve->edges[F.bottom].right;
			final_top_ind = top_curve->edges[F.top].right;
			
		}
		
		
		
		
		// check for degeneracy
		if (final_bottom_ind==final_top_ind) {
			// can simply set the top or bottom edge to be this one.  know it goes there.
			std::cout << "current edge is degenerate, " << final_bottom_ind <<"=="<<final_top_ind << std::endl;
			
			
			int current_edge = crit_slices[ii+zz].edge_w_midpt(final_bottom_ind);
			
			if (current_edge<0) {
				std::cout << "unable to find a degenerate edge in crit_slices[" << ii+zz << "] with midpoint " << final_bottom_ind << std::endl;
//				std::cout << "making new degenerate edge" << std::endl;
//				edge E(final_bottom_ind,final_bottom_ind,final_bottom_ind);
//				current_edge = crit_slices[ii+zz].add_edge(E);
				continue;
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
		set_mp(proj_top, &V.vertices[ final_top_ind ].projection_values->coord[1]);
		set_mp(proj_bottom, &V.vertices[ final_bottom_ind ].projection_values->coord[1]);
		
		// i think the projection values have already been thresholded
		if (solve_options.use_real_thresholding) {
			if (fabs(mpf_get_d(proj_top->i))<1e-10)
				mpf_set_str(proj_top->i,"0.0",10);
			
			if (fabs(mpf_get_d(proj_bottom->i))<1e-10)
				mpf_set_str(proj_bottom->i,"0.0",10);
			
		}
		
		
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
			for (std::set<int>::iterator poss_iter=possible_edges.begin(); poss_iter != possible_edges.end(); poss_iter++) {
				
				int qq = *poss_iter;
				
				bool matches_end = ((crit_slices[ii+zz].edges[qq].left == current_bottom_ind) || (crit_slices[ii+zz].edges[qq].right == current_bottom_ind));
				bool already_found = (found_edges.find(qq)!=found_edges.end());
				
				
				// we gotta be moving from lower to higher...  so temp > temp2 is required
				bool correct_interval = false;
				if (matches_end) {
					set_mp(temp , &V.vertices[ crit_slices[ii+zz].edges[qq].midpt].projection_values->coord[1]);
					set_mp(temp2, &V.vertices[ final_bottom_ind].projection_values->coord[1]);
					set_mp(temp3, &V.vertices[ final_top_ind].projection_values->coord[1]);
					correct_interval =  ( mpf_get_d(temp3->r) > mpf_get_d(temp->r)) && (mpf_get_d(temp->r) > mpf_get_d(temp2->r)) ;
				}
				
				bool degenerate = crit_slices[ii+zz].edges[qq].is_degenerate();
				
				if ( (!already_found) && matches_end && correct_interval && (!degenerate) ) { //
					candidates.push_back(qq);
					
					if (program_options.verbose_level>=1) {
						std::cout << "candidate [" << candidate_counter << "] = " <<
						crit_slices[ii+zz].edges[qq].left << " " << crit_slices[ii+zz].edges[qq].midpt << " " << crit_slices[ii+zz].edges[qq].right <<  std::endl;
					}
					
					
					candidate_counter++;
				}
				else
				{
					//                            if (!correct_interval) {
					//                                print_comp_matlab(temp3,"final_top_proj_1");
					//                                print_comp_matlab(temp ,"critslices[].proj_val_1");
					//                                print_comp_matlab(temp2,"final_bottom_proj_1");
					//                            }
					//							std::cout << "edge " << qq << " excluded: " << correct_interval << " direction, " << matches_end << " matches, " << already_found << "  already found" << std::endl;
					//
				}
				
			}
			
			std::cout << std::endl;
			
			
			if (candidate_counter==0) {
				std::cout << "found 0 candidates for bottom index " << current_bottom_ind << std::endl;
				break; // out of the while loop
			}
			
			for (int qq=0; qq<candidate_counter; qq++)
			{
				int current_edge = candidates[qq];
				
				if (program_options.verbose_level>=-1) {
					//					std::cout << "face #: " << this->num_faces << ", zz: " << zz << ", current_edge: " << current_edge << std::endl;
					std::cout << "tracking to these indices: " << final_bottom_ind << " " << crit_slices[ii+zz].edges[current_edge].midpt << " " << final_top_ind << std::endl;
				}
				
				
				
				
				//target midpoint e.w from paper.
				
				// this line has three index references in it.
				set_mp(proj_mid, &V.vertices[ crit_slices[ii+zz].edges[current_edge].midpt ].projection_values->coord[1]);
				
				
				
				if (solve_options.use_real_thresholding) {
					real_threshold(proj_mid,1e-10);
				}
				
				
				
				sub_mp(denom, proj_top, proj_bottom); // p2(w2) - p2(w0);
				sub_mp(numer, proj_mid, proj_bottom); // p2(e.w) - p2(w0);
				div_mp(md_config.v_target, numer, denom); // [p2(e.w) - p2(w0)] / [p2(w2) - p2(w0)]
				
				
				
				
				
				
				
				if (solve_options.verbose_level >= 0)
				{

					
					print_comp_matlab(md_config.u_target,"u_target");
					print_comp_matlab(md_config.v_target,"v_target");
					

				}
				
				
				if (solve_options.verbose_level >= 1){
					print_comp_matlab(proj_top,"upper");
					print_comp_matlab(proj_bottom,"lower");
					print_comp_matlab(proj_mid,"mid");
					
					print_comp_matlab(numer,"numer");
					print_comp_matlab(denom,"denom");
					
					print_comp_matlab(md_config.crit_val_left,"proj_val_left");
					print_comp_matlab(md_config.crit_val_right,"proj_val_right");
					
					
					set_one_mp(temp);
					sub_mp(temp2, temp, md_config.v_target);
					mul_mp(temp, temp2, proj_bottom);
					
					mul_mp(temp2, md_config.v_target, proj_top);
					
					add_mp(temp3, temp, temp2);
					
					print_comp_matlab(temp3,"proj_1_target_mid");
					
				}
				
				
				solver_output fillme;
				witness_set W_new;
				midpoint_solver_master_entry_point(W_midtrack, // carries with it the start points, and the linears.
												   fillme, // new data goes in here
												   md_config,
												   solve_options);
				
				fillme.get_noninfinite_w_mult_full(W_new);
				
				
				
				//print some display to screen
				if (solve_options.verbose_level >= 3)
				{
					print_point_to_screen_matlab(V.vertices[mid_slices[ii].edges[jj].right].pt_mp,"top_start");
					print_point_to_screen_matlab(V.vertices[mid_slices[ii].edges[jj].left].pt_mp,"bottom_start");
					print_point_to_screen_matlab(V.vertices[ crit_slices[ii+zz].edges[current_edge].midpt ].pt_mp,"midpoint_target");
				}
				
				
				
				
				
				// should get a single point back from this solver.
				
				if (W_new.num_points==0) {
					std::cout << color::red() << "midpoint tracker did not return any points :(" << color::console_default() << std::endl;
					possible_edges.erase(current_edge);
					continue;
				}
				
				vec_mp top_found, bottom_found;  init_vec_mp(top_found,md_config.num_top_vars());
				init_vec_mp(bottom_found,md_config.num_bottom_vars());
				
				// get only the midpoint coordinates out of the returned point
				for (int tt = 0; tt<this->num_variables; tt++) {
					set_mp(&found_point->coord[tt], &W_new.pts_mp[0]->coord[tt]);
				}
				
				int offset = md_config.num_mid_vars();
				// get only the bottom coordinates out of the returned point
				for (int tt = 0; tt<md_config.num_bottom_vars(); tt++) {
					set_mp(&bottom_found->coord[tt], &W_new.pts_mp[0]->coord[offset+tt]);
				}
				
				offset += md_config.num_bottom_vars();
				// get only the top coordinates out of the returned point
				for (int tt = 0; tt<md_config.num_top_vars(); tt++) {
					set_mp(&top_found->coord[tt], &W_new.pts_mp[0]->coord[offset+tt]);
				}
				
				//need to look the found point up in vertex set V
				int found_index = index_in_vertices(V, found_point);
				std::cout << index_in_vertices(V, bottom_found) << " bottom_found" << std::endl;
				std::cout << index_in_vertices(V, top_found) << " top_found" << std::endl;
				
				clear_vec_mp(bottom_found);
				clear_vec_mp(top_found);
				
				
				if (solve_options.verbose_level>=0) {
					std::cout << "found_index of point: " << found_index << std::endl;
				}
				
				if (solve_options.verbose_level>=1) {
					
					projection_value_homogeneous_input(temp, found_point, pi[1]);
					print_comp_matlab(temp, "found_point_proj_val");
					vec_mp tempvec; init_vec_mp(tempvec,0);
					dehomogenize(&tempvec, found_point);
					print_point_to_screen_matlab(tempvec, "found_point");
					clear_vec_mp(tempvec);
				}
				
				
				
				//search among the current edge possibilities for the one containing the found (mid) point
				int index_in_set = -1;
				for (std::set<int>::iterator possibility_iter=possible_edges.begin(); possibility_iter!=possible_edges.end(); possibility_iter++) {
					if (found_index == crit_slices[ii+zz].edges[*possibility_iter].midpt) {
						index_in_set = *possibility_iter;
						break;
					}
				}
				
				
				// search edges for the found point as a removed point.
				if (index_in_set < 0) {
					index_in_set = crit_slices[ii+zz].edge_w_removed(found_index);
					
					if (index_in_set>=0 && solve_options.verbose_level>=1) {
						std::cout << color::green() << "found point as removed point from edge " << index_in_set << color::console_default() << std::endl;
					}
				}
				
				
				if (index_in_set < 0 && solve_options.verbose_level>=0) {
					std::cout << color::red() << "did not find the indexed point as the midpoint of any current possibilitiy." << color::console_default() << std::endl;
				}
				
				//perhaps check for the point as left or right point of an edge?
				
				if (index_in_set>=0 && found_edges.find(index_in_set)==found_edges.end())
				{
					int next_edge = index_in_set; // index the *edge*
					
					if (program_options.verbose_level>=0) {
						std::cout << "added_edge " << next_edge << ", l m r: " << crit_slices[ii+zz].edges[next_edge].left << " " << crit_slices[ii+zz].edges[next_edge].midpt << " " << crit_slices[ii+zz].edges[next_edge].right << std::endl;
					}
					
					
					
					
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
					
					// erase both the currently testing edge, and the found one, from the list of possibilities.
					possible_edges.erase(current_edge);
					possible_edges.erase(next_edge);
					
					
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
					//didn't find, so simply remove from the list of possibilities.
					possible_edges.erase(current_edge);
				}
				
			}//re: for qq
			
		} // re: while...
		
		
		
		
	}
	
	clear_vec_mp(found_point);
	
	clear_mp(temp); clear_mp(temp2); clear_mp(temp3);
	clear_mp(numer); clear_mp(denom);
	
	clear_mp(proj_top); clear_mp(proj_bottom); clear_mp(proj_mid);
	
	return F;
}






















































void surface_decomposition::send(int target, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "surface::send" << std::endl;
#endif
	
	
	
	
	decomposition::send(target, mpi_config);
	int * buffer = new int[5];
	
	buffer[0] = num_edges;
	buffer[1] = num_faces;
	
	buffer[2] = mid_slices.size();
	buffer[3] = crit_slices.size();
	buffer[4] = num_singular_curves;
	
	MPI_Send(buffer, 5, MPI_INT, target, SURFACE, MPI_COMM_WORLD);
	delete [] buffer;
	
	
	
	
	
	for (unsigned int ii=0; ii<num_edges; ii++) {
		edges[ii].send(target, mpi_config);
	}
	
	for (unsigned int ii=0; ii<num_faces; ii++) {
		faces[ii].send(target, mpi_config);
	}
	
	
	for (unsigned int ii=0; ii!=mid_slices.size(); ii++) {
		mid_slices[ii].send(target, mpi_config);
	}
	
	for (unsigned int ii=0; ii!=crit_slices.size(); ii++) {
		crit_slices[ii].send(target, mpi_config);
	}
	
	crit_curve.send(target, mpi_config);
	sphere_curve.send(target, mpi_config);
	
	
	buffer = new int[2*num_singular_curves];
	int counter = 0;
	for (auto iter = singular_curves.begin(); iter!= singular_curves.end(); ++iter) {
		buffer[counter] = iter->first.first;
		buffer[counter+1] = iter->first.second;
		counter+=2;
	}
	MPI_Send(buffer, 2*num_singular_curves, MPI_INT, target, SURFACE, MPI_COMM_WORLD);
	delete [] buffer;
	
	for (auto iter = singular_curves.begin(); iter!= singular_curves.end(); ++iter) {
		iter->second.send(target, mpi_config);
	}
	
	
	return;
}


void surface_decomposition::receive(int source, parallelism_config & mpi_config)
{
	
#ifdef functionentry_output
	std::cout << "surface::receive" << std::endl;
#endif
	
	
	
	
	MPI_Status statty_mc_gatty;
	
	decomposition::receive(source, mpi_config);
	
	
	int * buffer = new int[5];
	MPI_Recv(buffer, 5, MPI_INT, source, SURFACE, MPI_COMM_WORLD, &statty_mc_gatty);
	int a, b, c, d;
	a = buffer[0];
	b = buffer[1];
	c = buffer[2];
	d = buffer[3];
	num_singular_curves = buffer[4];
	delete [] buffer;
	
	
	
	for (int ii=0; ii<a; ii++) {
		edge E;
		E.receive(source, mpi_config);
		add_edge(E);
	}
	
	for (int ii=0; ii<b; ii++) {
		face F;
		F.receive(source, mpi_config);
		add_face(F);
	}
	
	//TODO: rewrite this to prevent unnecessary copy operation.
	mid_slices.resize(c);
	for (int ii=0; ii<c; ii++) {
		mid_slices[ii].receive(source,mpi_config);
	}
	
	
	crit_slices.resize(d);
	for (int ii=0; ii<d; ii++) {
		crit_slices[ii].receive(source, mpi_config);
	}
	
	crit_curve.receive(source, mpi_config);
	sphere_curve.receive(source, mpi_config);
	
	
	buffer = new int[2*num_singular_curves];//
	MPI_Recv(buffer, 2*num_singular_curves, MPI_INT, source, SURFACE, MPI_COMM_WORLD, &statty_mc_gatty);
	
	// receive and unpack the buffer at same time
	int counter = 0;
	for (int ii=0; ii<num_singular_curves; ii++) {
		singular_curves[std::pair<int,int>(buffer[counter],buffer[counter+1])].receive(source, mpi_config);
		counter+=2;
	}
	
	delete [] buffer;
	return;
}




















void surface_decomposition::print(boost::filesystem::path base)
{
	
#ifdef functionentry_output
	std::cout << "surface::print" << std::endl;
#endif
	
	
	
	
	//	std::cout << "printing surface decomposition to folder " << base << std::endl;
	decomposition::print(base);
	
	
	
	boost::filesystem::path summaryname = base;
	summaryname /= "S.surf";
	FILE *OUT = safe_fopen_write(summaryname);
	fprintf(OUT,"%d %d %ld %ld\n\n", num_faces, num_edges, mid_slices.size(), crit_slices.size());
	fprintf(OUT,"%ld\n",singular_curves.size());
	for (auto iter = singular_curves.begin(); iter!=singular_curves.end(); ++iter) {
		fprintf(OUT,"%d %d ",iter->first.first,iter->first.second); // TODO:  clean up this first.first nonsense.  it is terrible.
	}
	fprintf(OUT,"\n");
	// what more to print here?
	fclose(OUT);
	
	summaryname = base;
	summaryname /= "F.faces";
	print_faces(summaryname);
	
	
	boost::filesystem::path curve_location = base;
	curve_location /= "curve";
	
	for (unsigned int ii=0; ii!=mid_slices.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_midslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		mid_slices[ii].output_main(specific_loc);
	}
	
	for (unsigned int ii=0; ii!=crit_slices.size(); ii++) {
		
		std::stringstream converter;
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_critslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		crit_slices[ii].output_main(specific_loc);
	}
	
	boost::filesystem::path specific_loc = curve_location;
	specific_loc += "_crit";
	crit_curve.output_main(specific_loc);
	
	specific_loc = curve_location;
	specific_loc += "_sphere";
	sphere_curve.output_main(specific_loc);
	
	for (auto iter = singular_curves.begin(); iter!=singular_curves.end(); ++iter) {
		
		std::stringstream converter;
		converter << curve_location.string() << "_singular_mult_" << iter->first.first << "_" << iter->first.second;
		
		specific_loc = converter.str();
		converter.clear(); converter.str("");
		
		iter->second.output_main(specific_loc);
	}
	
	
}







void surface_decomposition::print_faces(boost::filesystem::path outputfile)
{
#ifdef functionentry_output
	std::cout << "surface::print_faces" << std::endl;
#endif
	
	
	
	
    //	std::cout << "printing faces to file " << outputfile << std::endl;
	
	
	std::ofstream fout(outputfile.c_str());
	fout << num_faces << std::endl << std::endl;
	for (unsigned int ii=0; ii<num_faces; ++ii) {
		fout << faces[ii] << std::endl;
	}
	fout.close();
	
}






void surface_decomposition::read_faces(boost::filesystem::path load_from_me)
{
#ifdef functionentry_output
	std::cout << "surface::read_faces" << std::endl;
#endif
	
	
	
	std::ifstream fin(load_from_me.c_str());
	int temp_num_faces;
	fin >> temp_num_faces;
	for (int ii=0; ii<temp_num_faces; ii++) {
		face F;
		
		fin >> F;
		add_face(F);
	}
	
	fin.close();

	
	
	return;
}






void surface_decomposition::setup(boost::filesystem::path base)
{
	decomposition::setup(base / "decomp");
	
	std::vector<std::pair<int,int>> singular_multiplicities;
	int temp_num_crit, temp_num_mid;
	
	read_summary(singular_multiplicities,temp_num_mid, temp_num_crit, base / "S.surf");
	
	read_faces(base / "F.faces");
	
	mid_slices.resize(temp_num_mid);
	crit_slices.resize(temp_num_crit);
	
	
	
	boost::filesystem::path curve_location = base;
	curve_location /= "curve";
	
	std::stringstream converter;
	
	for (int ii=0; ii<temp_num_mid; ii++) {
		
		
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_midslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		mid_slices[ii].setup(specific_loc);
	}
	
	for (int ii=0; ii<temp_num_crit; ii++) {
		
		converter << ii;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_critslice_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		crit_slices[ii].setup(specific_loc);
	}
	
	
	for (auto iter = singular_multiplicities.begin(); iter!=singular_multiplicities.end(); ++iter) {
		
		converter << iter->first << "_" << iter->second;
		
		boost::filesystem::path specific_loc = curve_location;
		specific_loc += "_singular_mult_";
		specific_loc += converter.str();
		converter.clear(); converter.str("");
		
		singular_curves[*iter].setup(specific_loc);
	}
	
	boost::filesystem::path specific_loc = curve_location;
	specific_loc += "_crit";
	crit_curve.setup(specific_loc);
	
	specific_loc = curve_location;
	specific_loc += "_sphere";
	sphere_curve.setup(specific_loc);
	
	//	singular_curve.setup(specific_loc);
	
	
	
	
	return;
	
}





void read_summary(std::vector<std::pair<int,int>> & singular_multiplicities, int & temp_num_mid, int & temp_num_crit, boost::filesystem::path INfile)
{
	FILE *IN = safe_fopen_read(INfile);
	int temp_num_faces, temp_num_edges;
	
	fscanf(IN,"%d %d %d %d", &temp_num_faces, &temp_num_edges, &temp_num_mid, &temp_num_crit);
	
	int temp_num_sing;
	fscanf(IN,"%d",&temp_num_sing);
	int temp1,temp2;
	for (int ii=0; ii<temp_num_sing; ii++) {
		fscanf(IN,"%d %d",&temp1, &temp2);
		singular_multiplicities.push_back(std::pair<int,int>(temp1,temp2));
	}
	fclose(IN);
	
	
	return;
}


void create_sliced_system(boost::filesystem::path input_file, boost::filesystem::path output_file,
                          vec_mp * linears, int num_to_add,
                          const witness_set & W)
{
#ifdef functionentry_output
	std::cout << "create_sliced_system" << std::endl;
#endif
	
	
	
	
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
	fprintf(OUT,"END;\n\n\n\n");
	
    
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
	
#ifdef functionentry_output
	std::cout << "create_sphere_system" << std::endl;
#endif
	
	
	
	
	
	
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




















void face::send(int target, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "face::send" << std::endl;
#endif
	
	
	
	
	cell::send(target,mpi_config);
	
	int * buffer = (int *) br_malloc(10*sizeof(int));
	
	buffer[0] = left.size();
	buffer[1] = right.size();
	buffer[3] = top;
	buffer[4] = bottom;
	buffer[5] = num_left;
	buffer[6] = num_right;
	buffer[7] = system_name_bottom.size();
	buffer[8] = system_name_top.size();
	buffer[9] = crit_slice_index;
	
	MPI_Send(buffer, 10, MPI_INT, target, FACE, MPI_COMM_WORLD);
	free(buffer);
	
	if (num_left != left.size()) {
		std::cout << "left sizes for face DO NOT match" << std::endl;
	}
	if (num_right != right.size()) {
		std::cout << "left sizes for face DO NOT match" << std::endl;
	}
	
	if (num_left>0) {
		buffer = (int *) br_malloc(num_left*sizeof(int));
		for (unsigned int ii=0; ii<num_left; ii++) {
			buffer[ii] = left[ii];
		}
		MPI_Send(buffer, num_left, MPI_INT, target, FACE, MPI_COMM_WORLD);
		free(buffer);
	}
	
	
	
	if (num_right>0) {
		buffer = (int *) br_malloc(num_right*sizeof(int));
		for (unsigned int ii=0; ii<num_right; ii++) {
			buffer[ii] = right[ii];
		}
		MPI_Send(buffer, num_right, MPI_INT, target, FACE, MPI_COMM_WORLD);
		free(buffer);
	}
	
	
	
	std::string sendme = system_name_bottom;
	sendme.append(system_name_top);
	
	int num_to_send = sendme.size()+1;
	char * charbuff = new char[num_to_send];
	strcpy(charbuff, sendme.c_str());
	charbuff[num_to_send-1] = '\0';
	MPI_Send(&charbuff[0], num_to_send, MPI_CHAR, target, FACE, MPI_COMM_WORLD);
	delete [] charbuff;
	
	
	send_comp_mp(left_crit_val, target);
	send_comp_mp(right_crit_val, target);
	
	
	return;
	
}

void face::receive(int source, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "face::receive" << std::endl;
#endif
	
	
	int nchars_name_bottom, nchars_name_top;
	
	cell::receive(source,mpi_config);
	
	MPI_Status statty_mc_gatty;
	int * buffer= (int *) br_malloc(10*sizeof(int));
	
	MPI_Recv(buffer, 10, MPI_INT, source, FACE, MPI_COMM_WORLD, &statty_mc_gatty);
	
	int tmp_size_left = buffer[0];
	int tmp_size_right = buffer[1];
	top = buffer[3];
	bottom = buffer[4];
	num_left = buffer[5];
	num_right = buffer[6];
	nchars_name_bottom = buffer[7];
	nchars_name_top = buffer[8];
	crit_slice_index = buffer[9];
	
	free(buffer);
	
	
	if (tmp_size_left>0) {
		int * buffer2 = (int *) br_malloc(tmp_size_left*sizeof(int));
		MPI_Recv(buffer2, tmp_size_left, MPI_INT, source, FACE, MPI_COMM_WORLD, &statty_mc_gatty);
		for (int ii=0; ii<tmp_size_left; ii++) {
			left.push_back(buffer2[ii]);
		}
		free(buffer2);
	}
	
	
	
	if (tmp_size_right>0) {
		int * buffer3 = (int *) br_malloc(tmp_size_right*sizeof(int));
		MPI_Recv(buffer3, tmp_size_right, MPI_INT, source, FACE, MPI_COMM_WORLD, &statty_mc_gatty);
		for (int ii=0; ii<tmp_size_right; ii++) {
			right.push_back(buffer3[ii]);
		}
		free(buffer3);
	}
	
	
	
	
	
	receive_comp_mp(left_crit_val,source);
	receive_comp_mp(right_crit_val,source);
	
	
	
	char * charbuff = new char[nchars_name_bottom+nchars_name_top+1];
	
	MPI_Recv(charbuff, nchars_name_bottom+nchars_name_top+1, MPI_CHAR, source, FACE, MPI_COMM_WORLD, &statty_mc_gatty);
	
	std::stringstream converter;
	for (int jj=0; jj<nchars_name_bottom; ++jj) {
		converter << charbuff[jj];
	}
	system_name_bottom = converter.str();
	converter.clear();
	converter.str("");
	
	int offset = nchars_name_bottom;
	for (int jj=0; jj<nchars_name_top; ++jj) {
		converter << charbuff[jj+offset];
	}
	system_name_top = converter.str();
	
	delete [] charbuff;
	
	
	
	return;
	
}










