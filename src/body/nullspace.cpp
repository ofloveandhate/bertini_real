#include "nullspace.hpp"




int compute_crit_nullspace(SolverOutput & solve_out, // the returned value
						   const WitnessSet & W,
						   std::shared_ptr<SystemRandomizer> randomizer,
						   vec_mp *pi,
						   int ambient_dim,
						   int target_dim, // this should also be the number of vectors in the *pi entry
						   int target_crit_dim,
						   BertiniRealConfig & program_options,
						   SolverConfiguration & solve_options,
						   NullspaceConfiguration *ns_config)
{
	if (1) { // target_dim==target_crit_codim
		return compute_crit_nullspace_right(solve_out, // the returned value
											W,            // input the original witness set
											randomizer,
											pi,
											ambient_dim,  // dimension of ambient complex object
											program_options,
											solve_options,
											ns_config);
	}
	else{
		return compute_crit_nullspace_left(solve_out, // the returned value
										   W,            // input the original witness set
										   randomizer,
										   pi,
										   ambient_dim,  // dimension of ambient complex object
										   target_dim,   //  target dimension to find
										   target_crit_dim,   // COdimension of the critical set to find.
										   program_options,
										   solve_options,
										   ns_config);
	}
	
}




int compute_crit_nullspace_left(SolverOutput & solve_out, // the returned value
								const WitnessSet & W,
								std::shared_ptr<SystemRandomizer> randomizer,
								vec_mp *pi, // an array of projections, the number of which is the target dimensions
								int ambient_dim,
								int target_dim,
								int target_crit_codim,
								BertiniRealConfig & program_options,
								SolverConfiguration & solve_options,
								NullspaceConfiguration *ns_config)
{
	//many of the 1's here should be replaced by the number of patch equations, or the number of variable_groups
	
	int offset;
	WitnessSet Wtemp, Wtemp2;
	
	
	
	// get the max degree of the derivative functions.  this is unique to the left nullspace formulation, as the functions become mixed together
	int max_degree;
	
	
	nullspace_config_setup_left(ns_config,
								pi,
								ambient_dim,
								target_dim,
								target_crit_codim,
								&max_degree,
								randomizer,
								W,
								solve_options);
	
	
	solve_options.T.AMP_bound_on_degree = (double) max_degree+1;
	
	
	//
	///
	/////        end setup
	////////
	//////////////
	///////////////////////////
	
	
	
	
	if (max_degree==0) {
		// this will probably need tweaking when the dimension is higher than 1.  but who really is going to decompose a linear surface?
		solve_out.copy_patches(W);
		ns_concluding_modifications(solve_out, W, ns_config);
		
		
		std::cout << "the highest degree of any derivative equation is 0.  Returning empty SolverOutput." << std::endl;
		//then there cannot possibly be any critical points, with respect to ANY projection.  simply return an empty but complete set.
		
		return 0;
	}
	
	
	
	
	//  2.  Do a bunch of homotopies in $x$, each set of which will be followed by a single linear solve in $v$.
	if (program_options.verbose_level()>=3) {
		std::cout << "building up linprod start system for left nullspace" << std::endl;
	}
	
	// setup for the multilin moves
	//  these are for feeding into the multilin solver -- and that's it.  the majority will be overridden in the while loop as the start x linears
	vec_mp *multilin_linears = (vec_mp *) br_malloc(ambient_dim*sizeof(vec_mp)); // target dim is the number of linears in the input witness set
	for (int ii=0; ii<ambient_dim; ii++) {
		init_vec_mp2(multilin_linears[ii],W.num_variables(), solve_options.T.AMP_max_prec);
		multilin_linears[ii]->size = W.num_variables();
		vec_cp_mp(multilin_linears[ii], W.linear(ii));
	}
	
	MultilinConfiguration ml_config(solve_options,randomizer);
	
	
	
	// this is for performing the matrix inversion to get ahold of the $v$ values corresponding to $x$
	mat_mp tempmat;  init_mat_mp2(tempmat, ns_config->num_v_vars, ns_config->num_v_vars,solve_options.T.AMP_max_prec);
	tempmat->rows = tempmat->cols = ns_config->num_v_vars;
	
	offset = ns_config->num_v_vars-1;
	for (int jj=0; jj<ns_config->num_v_vars; jj++)
		set_mp(&tempmat->entry[offset][jj], &ns_config->v_patch->coord[jj]);
	
	// for holding the result of the matrix inversion
	vec_mp result; init_vec_mp2(result,ns_config->num_v_vars,solve_options.T.AMP_max_prec);  result->size = ns_config->num_v_vars;
	
	// use this for the matrix inversion
	vec_mp invert_wrt_me;
	init_vec_mp2(invert_wrt_me,ns_config->num_v_vars,solve_options.T.AMP_max_prec); invert_wrt_me->size = ns_config->num_v_vars;
	for (int ii=0; ii<ns_config->num_v_vars-1; ii++)
		set_zero_mp(&invert_wrt_me->coord[ii]); // set zero
	
	set_one_mp(&invert_wrt_me->coord[ns_config->num_v_vars-1]); // the last entry is set to 1 for the patch equation
	
	
	
	
	
	
	vec_mp temppoint;  init_vec_mp2(temppoint, ns_config->num_natural_vars +ns_config->num_synth_vars + ns_config->num_v_vars,solve_options.T.AMP_max_prec);
	temppoint->size = ns_config->num_natural_vars + ns_config->num_synth_vars + ns_config->num_v_vars;
	
	
	WitnessSet W_step_one;
	W_step_one.set_num_variables(W.num_variables());
	W_step_one.set_num_natural_variables(W.num_natural_variables());
	W_step_one.copy_patches(W);
	W_step_one.copy_names(W);
	
	
	WitnessSet W_linprod;
	W_linprod.set_num_variables(ns_config->num_natural_vars + ns_config->num_v_vars + ns_config->num_synth_vars);
	W_linprod.set_num_natural_variables(W.num_natural_variables());
	
	
	if (program_options.quick_run()<=1)
		solve_options.robust = true;
	else
		solve_options.robust = false;
	
	
	
	DoubleOdometer odo(ns_config->num_jac_equations, target_crit_codim, max_degree);
	
	int increment_status = 0;
	while (increment_status!=-1) { // current_absolute_index incremented at the bottom of loop
		
		if (program_options.verbose_level()>=5) {
			odo.print();
		}
		
		
		//copy in the linears for the solve
		for (int ii=0; ii<target_crit_codim; ii++) {
			vec_cp_mp(multilin_linears[ii], ns_config->starting_linears[odo.act_reg(ii)][odo.reg_val(ii)]);
		}
		// the remainder of the linears are left alone (stay stationary).
		
		if (program_options.verbose_level()>=6) {
			std::cout << "moving FROM this set:\n";
			for (unsigned int ii=0; ii<W.num_linears(); ii++) {
				print_point_to_screen_matlab(W.linear(ii),"L");
			}
			std::cout << "\nTO this set:\n";
			for (unsigned int ii=0; ii<W.num_linears(); ii++) {
				print_point_to_screen_matlab(multilin_linears[ii],"ELL");
			}
		}
		
		
		
		// actually solve WRT the linears
		
		
		SolverOutput fillme;
		multilin_solver_master_entry_point(W,         // WitnessSet
										   fillme, // the new data is put here!
										   multilin_linears,
										   ml_config,
										   solve_options);
		
		WitnessSet Wtemp;
		fillme.get_noninfinite_w_mult_full(Wtemp); // should be ordered
		
		W_step_one.merge(Wtemp, &solve_options.T);
		
		Wtemp.reset();
		
		
		
		
		//set the v_linears (M_i).
		
		for (int ii=0; ii<ns_config->num_v_vars-1; ii++) { // subtract one from upper limit because of the patch equation
			
			if (program_options.verbose_level()>=7)
			{
				std::cout << "copy into tempmat v_linears[" << odo.inact_reg(ii) << "]\n";
				print_point_to_screen_matlab(ns_config->v_linears[odo.inact_reg(ii)], "v_linears");
			}
			for (int jj=0; jj<ns_config->num_v_vars; jj++)
				set_mp(&tempmat->entry[ii][jj], &ns_config->v_linears[odo.inact_reg(ii)]->coord[jj]);
		}
		
		
		
		increment_status = odo.increment();  // increment the tracking indices
		
		
		// if it's time to move on to next function combo, must commit the points
		if (increment_status!=0)
		{
			
			// invert the matrix for the v variables.
			matrixSolve_mp(result, tempmat,  invert_wrt_me);
			
			
			offset = ns_config->num_natural_vars+ns_config->num_synth_vars;
			for (int jj=0; jj<ns_config->num_v_vars; jj++)
				set_mp(&temppoint->coord[jj+offset], &result->coord[jj]);
			
			
			for (unsigned int ii=0; ii<W_step_one.num_points(); ii++) {
				for (int jj=0; jj<ns_config->num_natural_vars+ns_config->num_synth_vars; jj++) {
					set_mp(&temppoint->coord[jj], & W_step_one.point(ii)->coord[jj]);
				}
				W_linprod.add_point(temppoint);
			}
			
			
			
			W_step_one.reset();
			W_step_one.set_num_variables(W.num_variables());
			W_step_one.set_num_natural_variables(W.num_natural_variables());
			W_step_one.copy_patches(W);  // necessary?
			W_step_one.copy_names(W); // necessary?
			
		}
		
	}
	
	
	
	for (int ii=0; ii<ambient_dim; ii++) {
		clear_vec_mp(multilin_linears[ii]);
	}
	free(multilin_linears);
	
	clear_vec_mp(temppoint);
	
	int num_before = W_linprod.num_points();
	W_linprod.sort_for_unique(&solve_options.T);
	if (num_before - W_linprod.num_points()>0) {
		std::cout << "there were non-unique start points" << std::endl;
		mypause();
	}
	
	W_linprod.set_num_natural_variables(ns_config->num_v_vars+ns_config->num_synth_vars);
	W_linprod.copy_patches(W);
	W_linprod.copy_names(W);
	
	//set some solver options
	
	if (program_options.quick_run()<=0)
		solve_options.robust = true;
	else
		solve_options.robust = false;
	
	
	
	solve_options.use_midpoint_checker = 0;
	
	
	
	if (program_options.verbose_level()>=6)
		ns_config->print();
	
	
	
	if (program_options.verbose_level()>=3) {
		std::cout << "running nullspace method" << std::endl;
	}
	
	
	nullspacejac_solver_master_entry_point(solve_options.T.MPType,
										   W_linprod, // carries with it the start points, but not the linears.
										   solve_out,   // the created data goes in here.
										   ns_config,
										   solve_options);
	
	
	ns_concluding_modifications(solve_out, W, ns_config);
	
	
	clear_mat_mp(tempmat);
	clear_vec_mp(invert_wrt_me);
	clear_vec_mp(result);
	
	
	return SUCCESSFUL;
}
















void nullspace_config_setup_left(NullspaceConfiguration *ns_config,
								 vec_mp *pi, // an array of projections, the number of which is the target dimensions
								 int ambient_dim,
								 int target_dim,
								 int target_crit_codim,
								 int *max_degree, // a pointer to the value
								 std::shared_ptr<SystemRandomizer> randomizer,
								 const WitnessSet & W,
								 SolverConfiguration & solve_options)
{
	
	ns_config->set_side(nullspace_handedness::LEFT);
	
	int toss;
	parse_input_file(W.input_filename(), &toss); // re-create the parsed files for the stuffs (namely the SLP).
	
	ns_config->set_randomizer(randomizer); // set the pointer.  this randomizer is for the underlying system.
	
	
	
	*max_degree = randomizer->max_degree()-1; // minus one for differentiated degree
	ns_config->max_degree = *max_degree;
	
	
	
	// set some integers
	
	ns_config->num_projections = ambient_dim - target_crit_codim + 1;
	
	
	ns_config->num_v_vars = (W.num_natural_variables()-1) - ambient_dim + ns_config->num_projections;
	
	
	ns_config->num_synth_vars = W.num_synth_vars(); // this may get a little crazy if we chain into this more than once.  this code is written to be called into only one time beyond the first.
	ns_config->num_natural_vars = W.num_natural_variables();
	
	ns_config->ambient_dim = ambient_dim;
	ns_config->target_dim = target_dim;
	ns_config->target_crit_codim = target_crit_codim;
	
	
	ns_config->target_projection = (vec_mp *) br_malloc(ns_config->num_projections * sizeof(vec_mp));
	for (int ii=0; ii<ns_config->num_projections; ii++) {
		init_vec_mp2(ns_config->target_projection[ii], W.num_variables(),solve_options.T.AMP_max_prec);
		ns_config->target_projection[ii]->size = W.num_variables();
		vec_cp_mp(ns_config->target_projection[ii], pi[ii]);
	}
	
	
	
	ns_config->num_jac_equations = (ns_config->num_natural_vars - 1);// N-1;  the subtraction of 1 is for the 1 hom-var.
																	 // me must omit any previously added synthetic vars.
	
	ns_config->num_additional_linears = target_dim - target_crit_codim;
	
	ns_config->num_v_linears = ns_config->num_jac_equations;
	
	
	
	// this check is correct.
	int check_num_f = randomizer->num_rand_funcs() + ns_config->num_jac_equations + ns_config->num_additional_linears + W.num_patches() + 1; // +1 for v patch from this incoming computation
	int check_num_v = ns_config->num_natural_vars + ns_config->num_synth_vars + ns_config->num_v_vars;
	if (check_num_f != check_num_v) {
		std::cout << color::red();
		std::cout << "mismatch in number of equations...\n" << std::endl;
		std::cout << "left: " << check_num_f << " right " << check_num_v << std::endl;
		std::cout << color::console_default();
		throw std::logic_error("logic error in nullspace_left");
	}
	else{
		//		std::cout << "passed the square check in nullspace_left setup" << std::endl;
	}
	
	
	
	
	
	
	
	// set up the linears in $v$  ( the M_i linears)
	ns_config->v_linears = (vec_mp *)br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
	for (int ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_mp2(ns_config->v_linears[ii],ns_config->num_v_vars,solve_options.T.AMP_max_prec);
		ns_config->v_linears[ii]->size = ns_config->num_v_vars;
		for (int jj=0; jj<ns_config->num_v_vars; jj++){
			get_comp_rand_mp(&ns_config->v_linears[ii]->coord[jj]); // should this be real? no.
		}
	}
	
	
	// the last of the linears will be used for the slicing, and passed on to later routines
	int offset = ambient_dim - target_dim + target_crit_codim;
	
	
	ns_config->additional_linears_terminal = (vec_mp *)br_malloc((ns_config->num_additional_linears)*sizeof(vec_mp));
	ns_config->additional_linears_starting = (vec_mp *)br_malloc((ns_config->num_additional_linears)*sizeof(vec_mp));
	
	for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
		init_vec_mp2(ns_config->additional_linears_terminal[ii],W.num_variables(),solve_options.T.AMP_max_prec);
		ns_config->additional_linears_terminal[ii]->size = W.num_variables();
		if (1) {
			for (int jj=0; jj<W.num_natural_variables(); jj++){
				get_comp_rand_mp(&ns_config->additional_linears_terminal[ii]->coord[jj]); // should this be real?  no.
			}
			for (int jj=W.num_natural_variables(); jj<W.num_variables(); jj++) {
				set_zero_mp(&ns_config->additional_linears_terminal[ii]->coord[jj]);
			}
		}
		else{
			for (int jj=0; jj<W.num_variables(); jj++){
				get_comp_rand_mp(&ns_config->additional_linears_terminal[ii]->coord[jj]); // should this be real?  no.
			}
		}
		
		init_vec_mp2(ns_config->additional_linears_starting[ii],W.num_variables(),solve_options.T.AMP_max_prec);
		ns_config->additional_linears_starting[ii]->size = W.num_variables();
		vec_cp_mp(ns_config->additional_linears_starting[ii], W.linear(ii+offset));
	}
	
	
	// set up the patch in $v$.  we will include this in an inversion matrix to get the starting $v$ values.
	init_vec_mp2(ns_config->v_patch,ns_config->num_v_vars,solve_options.T.AMP_max_prec);
	ns_config->v_patch->size = ns_config->num_v_vars;
	for (int ii=0; ii<ns_config->num_v_vars; ii++) {
		get_comp_rand_mp(&ns_config->v_patch->coord[ii]);
	}
	
	
	
	mat_mp temp_getter;  init_mat_mp2(temp_getter,*max_degree, W.num_variables(),solve_options.T.AMP_max_prec);
	temp_getter->rows = *max_degree;
	temp_getter->cols = W.num_variables();
	
	//the 'ns_config->starting_linears' will be used for the x variables.  we will homotope to these $k-\ell$ at a time
	ns_config->starting_linears = (vec_mp **)br_malloc( ns_config->num_jac_equations*sizeof(vec_mp *));
	for (int ii=0; ii<ns_config->num_jac_equations; ii++) {
		ns_config->starting_linears[ii] = (vec_mp *) br_malloc((*max_degree)*sizeof(vec_mp));
		
		if (1) {
			make_matrix_random_mp(temp_getter,*max_degree, W.num_natural_variables(), solve_options.T.AMP_max_prec); // this matrix is nearly orthogonal
			
			for (int jj=0; jj<(*max_degree); jj++) {
				init_vec_mp2(ns_config->starting_linears[ii][jj],W.num_variables(),solve_options.T.AMP_max_prec);
				ns_config->starting_linears[ii][jj]->size = W.num_variables();
				
				for (int kk=0; kk<W.num_natural_variables(); kk++) {
					set_mp(&ns_config->starting_linears[ii][jj]->coord[kk], &temp_getter->entry[jj][kk]);
				}
				
				for (int kk=W.num_natural_variables(); kk<W.num_variables(); kk++) {
					set_zero_mp(&ns_config->starting_linears[ii][jj]->coord[kk]);
				}
			}
		}
		else{
			make_matrix_random_mp(temp_getter,*max_degree, W.num_variables(), solve_options.T.AMP_max_prec); // this matrix is nearly orthogonal
			
			for (int jj=0; jj<(*max_degree); jj++) {
				init_vec_mp2(ns_config->starting_linears[ii][jj],W.num_variables(),solve_options.T.AMP_max_prec);
				ns_config->starting_linears[ii][jj]->size = W.num_variables();
				
				for (int kk=0; kk<W.num_variables(); kk++) {
					set_mp(&ns_config->starting_linears[ii][jj]->coord[kk], &temp_getter->entry[jj][kk]);
				}
			}
		}
	}
	
	
	clear_mat_mp(temp_getter);
	
	
	
	return;
}
















int compute_crit_nullspace_right(SolverOutput & solve_out, // the returned value
						   const WitnessSet & W,
						   std::shared_ptr<SystemRandomizer> randomizer,
						   vec_mp *pi, // an array of projections, the number of which is the target dimensions
						   int ambient_dim,
						   BertiniRealConfig & program_options,
						   SolverConfiguration & solve_options,
						   NullspaceConfiguration *ns_config)
{
	//many of the 1's here should be replaced by the number of patch equations, or the number of variable_groups
	
	
	
	
	
	// get the max degree of the derivative functions.  this is unique to the left nullspace formulation, as the functions become mixed together
	int max_degree;
	
	
	nullspace_config_setup_right(ns_config,
						   pi,
						   ambient_dim,
						   &max_degree,
						   randomizer,
						   W,
						   solve_options);
	
	
	
	
	
	if (max_degree==0) {
		// this will probably need tweaking when the dimension is higher than 1.  but who really is going to decompose a linear surface?
		solve_out.copy_patches(W);
		ns_concluding_modifications(solve_out, W, ns_config);
		
		
		std::cout << "the highest degree of any derivative equation is 0.  Returning empty SolverOutput." << std::endl;
		//then there cannot possibly be any critical points, with respect to ANY projection.  simply return an empty but complete set.
		
		return 0;
	}
	
	
	
	
	
	
	//
	///
	/////        end setup
	////////
	//////////////
	///////////////////////////
	
	
	
	
	
	
	int offset;
	WitnessSet Wtemp, Wtemp2;
	
	
	//  2.  Do a bunch of homotopies in $x$, each set of which will be followed by a single linear solve in $v$.
	if (program_options.verbose_level()>=3) {
		std::cout << "building up linprod start system for left nullspace" << std::endl;
	}
	
	// setup for the multilin moves
	//  these are for feeding into the multilin solver -- and that's it.  the majority will be overridden in the while loop as the start x linears
	vec_mp *multilin_linears = (vec_mp *) br_malloc(W.num_linears()*sizeof(vec_mp)); // target dim is the number of linears in the input witness set
	for (unsigned int ii=0; ii<W.num_linears(); ii++) {
		init_vec_mp2(multilin_linears[ii],W.num_variables(), solve_options.T.AMP_max_prec);
		multilin_linears[ii]->size = W.num_variables();
		vec_cp_mp(multilin_linears[ii], W.linear(ii));
	}
	
	MultilinConfiguration ml_config(solve_options,randomizer);
	
	
	
	// this is for performing the matrix inversion to get ahold of the $v$ values corresponding to $x$
	mat_mp tempmat;  init_mat_mp2(tempmat, ns_config->num_v_vars, ns_config->num_v_vars,solve_options.T.AMP_max_prec);
	tempmat->rows = tempmat->cols = ns_config->num_v_vars;
	
	offset = ns_config->num_v_vars-1;
	for (int jj=0; jj<ns_config->num_v_vars; jj++)
		set_mp(&tempmat->entry[offset][jj], &ns_config->v_patch->coord[jj]);
	
	
	// for holding the result of the matrix inversion
	vec_mp result; init_vec_mp2(result,ns_config->num_v_vars,solve_options.T.AMP_max_prec);  result->size = ns_config->num_v_vars;
	
	// use this for the matrix inversion
	vec_mp invert_wrt_me;
	init_vec_mp2(invert_wrt_me,ns_config->num_v_vars,solve_options.T.AMP_max_prec); invert_wrt_me->size = ns_config->num_v_vars;
	for (int ii=0; ii<ns_config->num_v_vars-1; ii++)
		set_zero_mp(&invert_wrt_me->coord[ii]); // set zero
	
	set_one_mp(&invert_wrt_me->coord[ns_config->num_v_vars-1]); // the last entry is set to 1 for the patch equation
	
	
	
	
	
	
	vec_mp temppoint;  init_vec_mp2(temppoint, ns_config->num_natural_vars +ns_config->num_synth_vars + ns_config->num_v_vars,solve_options.T.AMP_max_prec);
	temppoint->size = ns_config->num_natural_vars + ns_config->num_synth_vars + ns_config->num_v_vars;
	
	
	WitnessSet W_step_one;
	W_step_one.set_num_variables(W.num_variables());
	W_step_one.set_num_natural_variables(W.num_natural_variables());
	W_step_one.copy_patches(W);
	W_step_one.copy_names(W);
	
	
	WitnessSet W_linprod;
	W_linprod.set_num_variables(ns_config->num_natural_vars + ns_config->num_v_vars + ns_config->num_synth_vars);
	W_linprod.set_num_natural_variables(W.num_natural_variables());
	
	
	if (program_options.quick_run()<=1)
		solve_options.robust = true;
	else
		solve_options.robust = false;
	
	
	
	
	
	for (int ii=0; ii<randomizer->num_rand_funcs(); ii++) {
		
		int differentiated_degree = randomizer->randomized_degree(ii)-1; // the -1 is for differentiating.  this could be 0.
		
		if (differentiated_degree==0) {
			continue;
		}
		else{
			for (int jj=0; jj<differentiated_degree; jj++) {
				
				//copy in the linear for the solve
				vec_cp_mp(multilin_linears[0], ns_config->starting_linears[ii][jj]);
				
				// the remainder of the linears are left alone (stay stationary).
				
				if (program_options.verbose_level()>=6) {
					std::cout << "moving FROM this set:\n";
					for (unsigned int ii=0; ii<W.num_linears(); ii++) {
						print_point_to_screen_matlab(W.linear(ii),"L");
					}
					std::cout << "\nTO this set:\n";
					for (unsigned int ii=0; ii<W.num_linears(); ii++) {
						print_point_to_screen_matlab(multilin_linears[ii],"ELL");
					}
				}
				
				
				
				// actually solve WRT the linears
				
				
				SolverOutput fillme;
				multilin_solver_master_entry_point(W,         // WitnessSet
												   fillme, // the new data is put here!
												   multilin_linears,
												   ml_config,
												   solve_options);
				
				WitnessSet Wtemp;
				fillme.get_noninfinite_w_mult_full(Wtemp); // should be ordered
				
				W_step_one.merge(Wtemp, &solve_options.T);
				
				Wtemp.reset();
				
				
				
				
				//set the v_linears (M_i).
				
				
				
				
			}
			
			
			
			int curr_index = 0;
			for (int kk=0; kk<ns_config->num_v_vars; kk++) { // subtract one from upper limit because of the patch equation
				if (kk!=ii) {
					
					for (int mm=0; mm<ns_config->num_v_vars; mm++){
						set_mp(&tempmat->entry[curr_index][mm], &ns_config->v_linears[kk]->coord[mm]);}
					
					curr_index++;
				}
				
			}
			
			
			// invert the matrix for the v variables.
			matrixSolve_mp(result, tempmat,  invert_wrt_me);
			
			//set the bottom part of the temppoint, which will be a startpoint for the nullspace call later.
			offset = ns_config->num_natural_vars+ns_config->num_synth_vars;
			for (int mm=0; mm<ns_config->num_v_vars; mm++)
				set_mp(&temppoint->coord[mm+offset], &result->coord[mm]);
			
			
			//set the top part, x, of the start point, and copy it in.
			for (unsigned int kk=0; kk<W_step_one.num_points(); kk++) {
				for (int mm=0; mm<ns_config->num_natural_vars+ns_config->num_synth_vars; mm++) {
					set_mp(&temppoint->coord[mm], & W_step_one.point(kk)->coord[mm]);
				}
				W_linprod.add_point(temppoint);
			}
			
			
			
			W_step_one.reset();
			W_step_one.set_num_variables(W.num_variables());
			W_step_one.set_num_natural_variables(W.num_natural_variables());
			W_step_one.copy_patches(W);  // necessary?
			W_step_one.copy_names(W); // necessary?
		}
		
		
		
		
		
		
		
	}
	
	
	for (int ii=0; ii<1; ii++) {
		clear_vec_mp(multilin_linears[ii]);
	}
	free(multilin_linears);
	
	clear_vec_mp(temppoint);
	
	
	
	int num_before = W_linprod.num_points();
	W_linprod.sort_for_unique(&solve_options.T);
	if (num_before - W_linprod.num_points()>0) {
		std::cout << "there were non-unique start points" << std::endl;
		mypause();
	}
	
	
	
	W_linprod.set_num_natural_variables(ns_config->num_v_vars+ns_config->num_synth_vars);
	W_linprod.copy_patches(W);
	W_linprod.copy_names(W);
	
	//set some solver options
	
	if (program_options.quick_run()<=0)
		solve_options.robust = true;
	else
		solve_options.robust = false;
	
	
	
	solve_options.use_midpoint_checker = 0;
	
	
	
	if (program_options.verbose_level()>=6)
		ns_config->print();
	
	
	
	if (program_options.verbose_level()>=3) {
		std::cout << "running nullspace right method" << std::endl;
	}
	
	
	
	
	nullspacejac_solver_master_entry_point(solve_options.T.MPType,
										   W_linprod, // carries with it the start points, but not the linears.
										   solve_out,   // the created data goes in here.
										   ns_config,
										   solve_options);
	
	
	ns_concluding_modifications(solve_out, W, ns_config);
	
	
	clear_mat_mp(tempmat);
	clear_vec_mp(invert_wrt_me);
	clear_vec_mp(result);
	
	
	return SUCCESSFUL;
}









void nullspace_config_setup_right(NullspaceConfiguration *ns_config,
								  vec_mp *pi, // an array of projections, the number of which is the target dimensions
								  int ambient_dim,
								  int *max_degree, // a pointer to the value
								  std::shared_ptr<SystemRandomizer> randomizer,
								  const WitnessSet & W,
								  SolverConfiguration & solve_options)
{
	
	ns_config->set_side(nullspace_handedness::RIGHT);
	
	
	
	int toss;
	parse_input_file(W.input_filename(), &toss); // re-create the parsed files for the stuffs (namely the SLP).
	
	ns_config->set_randomizer(randomizer); // set the pointer.  this randomizer is for the underlying system.
	
	*max_degree = randomizer->max_degree()-1; // minus one for differentiated degree
	ns_config->max_degree = *max_degree;
	
	
	
	// set some integers
	
	ns_config->num_projections = ambient_dim;
	
	ns_config->num_v_vars = W.num_natural_variables()-1;
	
	
	ns_config->num_synth_vars = W.num_synth_vars(); // this may get a little crazy if we chain into this more than once.  this code is written to be called into only one time beyond the first.
	ns_config->num_natural_vars = W.num_natural_variables();
	
	ns_config->ambient_dim = ambient_dim;
	
	
	ns_config->target_projection = (vec_mp *) br_malloc(ns_config->num_projections * sizeof(vec_mp));
	for (int ii=0; ii<ns_config->num_projections; ii++) {
		init_vec_mp2(ns_config->target_projection[ii], W.num_variables(),solve_options.T.AMP_max_prec);
		ns_config->target_projection[ii]->size = W.num_variables();
		vec_cp_mp(ns_config->target_projection[ii], pi[ii]);
	}
	
	
	
	ns_config->num_jac_equations = (ns_config->num_natural_vars - 1);// N-1;  the subtraction of 1 is for the 1 hom-var.
																	 // me must omit any previously added synthetic vars.
	
	ns_config->num_additional_linears = ambient_dim-1;
	
	ns_config->num_v_linears = ns_config->num_jac_equations;
	
	
	
	// this check is correct.
	int check_num_func = randomizer->num_rand_funcs() + ns_config->num_jac_equations + ns_config->num_additional_linears + W.num_patches() + 1; // +1 for v patch from this incoming computation
	int check_num_vars = ns_config->num_natural_vars + ns_config->num_synth_vars + ns_config->num_v_vars;
	if (check_num_func != check_num_vars) {
		std::cout << color::red();
		std::cout << "mismatch in number of equations...\n" << std::endl;
		std::cout << "left: " << check_num_func << " right " << check_num_vars << std::endl;
		std::cout << color::console_default();
		throw std::logic_error("logic error in nullspace_left");
	}
//	else{
//		std::cout << "passed the square check in nullspace_right setup" << std::endl;
//	}
	
	
	
	
	
	
	
	// set up the linears in $v$  ( the M_i linears)
	ns_config->v_linears = (vec_mp *)br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
	for (int ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_mp2(ns_config->v_linears[ii],ns_config->num_v_vars,solve_options.T.AMP_max_prec);
		ns_config->v_linears[ii]->size = ns_config->num_v_vars;
		for (int jj=0; jj<ns_config->num_v_vars; jj++){
			get_comp_rand_mp(&ns_config->v_linears[ii]->coord[jj]); // should this be real? no.
		}
	}
	
	
	// the last of the linears will be used for the slicing, and passed on to later routines
	int offset = 1;
	
	
	ns_config->additional_linears_terminal = (vec_mp *)br_malloc((ns_config->num_additional_linears)*sizeof(vec_mp));
	ns_config->additional_linears_starting = (vec_mp *)br_malloc((ns_config->num_additional_linears)*sizeof(vec_mp));
	
	for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
		init_vec_mp2(ns_config->additional_linears_terminal[ii],W.num_variables(),solve_options.T.AMP_max_prec);
		ns_config->additional_linears_terminal[ii]->size = W.num_variables();
		
		for (int jj=0; jj<W.num_natural_variables(); jj++){
			get_comp_rand_mp(&ns_config->additional_linears_terminal[ii]->coord[jj]); // should this be real?  no.
		}
		for (int jj=W.num_natural_variables(); jj<W.num_variables(); jj++) {
			set_zero_mp(&ns_config->additional_linears_terminal[ii]->coord[jj]);
		}
		
		
		init_vec_mp2(ns_config->additional_linears_starting[ii],W.num_variables(),solve_options.T.AMP_max_prec);
		ns_config->additional_linears_starting[ii]->size = W.num_variables();
		vec_cp_mp(ns_config->additional_linears_starting[ii], W.linear(ii+offset));
	}
	
	
	// set up the patch in $v$.  we will include this in an inversion matrix to get the starting $v$ values.
	init_vec_mp2(ns_config->v_patch,ns_config->num_v_vars,solve_options.T.AMP_max_prec);
	ns_config->v_patch->size = ns_config->num_v_vars;
	for (int ii=0; ii<ns_config->num_v_vars; ii++) {
		get_comp_rand_mp(&ns_config->v_patch->coord[ii]);
	}
	
	
	
	mat_mp temp_getter;  init_mat_mp2(temp_getter,0, 0,solve_options.T.AMP_max_prec);
	temp_getter->rows = 0; temp_getter->cols = 0;
	
	
	
	//the 'ns_config->starting_linears' will be used for the x variables.  we will homotope to these 1 at a time
	ns_config->starting_linears = (vec_mp **)br_malloc( randomizer->num_rand_funcs()*sizeof(vec_mp *));
	for (int ii=0; ii<randomizer->num_rand_funcs(); ii++) {
		
		unsigned int curr_degree = randomizer->randomized_degree(ii)-1;
		
		
		ns_config->starting_linears[ii] = (vec_mp *) br_malloc(curr_degree*sizeof(vec_mp));
		
		make_matrix_random_mp(temp_getter,curr_degree, W.num_natural_variables(), solve_options.T.AMP_max_prec); // this matrix is nearly orthogonal
		
		for (unsigned int jj=0; jj<curr_degree; jj++) {
			init_vec_mp2(ns_config->starting_linears[ii][jj],W.num_variables(),solve_options.T.AMP_max_prec);
			ns_config->starting_linears[ii][jj]->size = W.num_variables();
			
			for (int kk=0; kk<W.num_natural_variables(); kk++) {
				set_mp(&ns_config->starting_linears[ii][jj]->coord[kk], &temp_getter->entry[jj][kk]);
			}
			
			for (int kk=W.num_natural_variables(); kk<W.num_variables(); kk++) {
				set_zero_mp(&ns_config->starting_linears[ii][jj]->coord[kk]);
			}
		}
		
	}
	
	
	clear_mat_mp(temp_getter);
	
	
	
	return;
}





void ns_concluding_modifications(SolverOutput & solve_out,
								 const WitnessSet & W,
								 NullspaceConfiguration * ns_config)
{
	solve_out.num_variables  = ns_config->num_natural_vars + ns_config->num_v_vars;
	solve_out.num_natural_vars = W.num_natural_variables();
	
	solve_out.add_patch(ns_config->v_patch);
	
	for (int ii=0; ii< ns_config->num_additional_linears; ii++) {
		solve_out.add_linear(ns_config->additional_linears_terminal[ii]);
	}
}










void create_nullspace_system(boost::filesystem::path output_name,
							 boost::filesystem::path input_name,
							 BertiniRealConfig & program_options,
							 NullspaceConfiguration *ns_config)
/***************************************************************\
 * USAGE: setup input file for one deflation iteration           *
 * ARGUMENTS: number of declaration statments, name of file,     *
 *  command to run Matlab, and dimension of nullspace            *
 * RETURN VALUES: creates new file                               *
 * NOTES:                                                        *
 \***************************************************************/
{
	
	int *declarations = NULL;
	partition_parse(&declarations, input_name, "func_input_real" , "config_real" ,0); // the 0 means not self conjugate mode
	
	
	
	int ii, numVars = 0, numFuncs = 0, numConstants = 0;
	int *lineVars = NULL, *lineFuncs = NULL, *lineConstants = NULL;
	char ch, *str = NULL, **vars = NULL, **funcs = NULL, **consts = NULL;
	FILE *IN = NULL, *OUT = NULL;
	
	
	// move the file & open it
	IN = safe_fopen_read("func_input_real");
	// setup variables
	if (declarations[0] > 0)
	{ // using variable_group
		parse_names(&numVars, &vars, &lineVars, IN,const_cast< char *>("variable_group"), declarations[0]);
	}
	else
	{ // using hom_variable_group
		parse_names(&numVars, &vars, &lineVars, IN,const_cast< char *>("hom_variable_group"), declarations[1]);
	}
	
	// setup constants
	rewind(IN);
	parse_names(&numConstants, &consts, &lineConstants, IN,const_cast< char *>("constant"), declarations[8]);
	
	// setup functions
	rewind(IN);
	parse_names(&numFuncs, &funcs, &lineFuncs, IN, const_cast< char *>("function"), declarations[9]);
	
	fclose(IN);
	
	// setup Matlab script
	
	if (0) {
		createMatlabDerivative("matlab_nullspace_system.m", "func_input_real",
							   ns_config,
							   numVars, vars, lineVars, numConstants, consts, lineConstants, numFuncs, funcs, lineFuncs);
	}
	else{
		create_matlab_determinantal_system("matlab_nullspace_system.m", "func_input_real",
										   ns_config,
										   numVars, vars, lineVars, numConstants, consts, lineConstants, numFuncs, funcs, lineFuncs);
	}
	
	
	
	
	
	// run Matlab script
	std::stringstream converter;
	converter << program_options.matlab_command() << "matlab_nullspace_system";
	system(converter.str().c_str());
	converter.clear(); converter.str("");
	
	
	
	// setup new file
	
	OUT = safe_fopen_write(output_name.c_str());
	
	
	//	fprintf(OUT, "CONFIG\n");
	//
	//	IN = safe_fopen_read("config_real");
	//	copyfile(IN,OUT);
	//	fclose(IN);
	//	fprintf(OUT, "\nEND;");
	fprintf(OUT, "\n\nINPUT\n");
	
	
	std::string constants = just_constants("func_input_real",numConstants,consts,lineConstants);
	
	fprintf(OUT, "%s\n\n",constants.c_str());
	//	IN = safe_fopen_read("func_input_real");
	//	copyfile(IN,OUT);
	//	fclose(IN);
	//	fprintf(OUT,"\n");
	
	//	if (numConstants>0) {
	//
	//		fprintf(OUT,"constant ");
	//		for (ii=0;ii<numConstants; ii++){
	//			fprintf(OUT,"%s",consts[ii]);
	//			if (ii!=numConstants-1)
	//				fprintf(OUT,", ");
	//			else
	//				fprintf(OUT,";\n");
	//		}
	//	}
	
	
	// this commented block produces TWO variables groups.
	//	fprintf(OUT,"variable_group ");
	//	for (ii=0;ii<numVars; ii++){
	//		fprintf(OUT,"%s",vars[ii]);
	//		(ii==numVars-1) ? fprintf(OUT,";\n") : fprintf(OUT,", ");
	//	}
	//
	//	fprintf(OUT,"hom_variable_group ");
	//	for (ii=0; ii<(ns_config->num_v_vars); ii++) {
	//		fprintf(OUT,"synth%d",ii+1);
	//		(ii==ns_config->num_v_vars-1) ? fprintf(OUT,";\n") : fprintf(OUT,", ");
	//	}
	
	if (0) {
		fprintf(OUT,"variable_group ");
		for (ii=0;ii<numVars; ii++){
			fprintf(OUT,"%s",vars[ii]);
			(ii==numVars-1) ? fprintf(OUT,";\n") : fprintf(OUT,", ");
		}
		
		fprintf(OUT,"hom_variable_group ");
		for (ii=0; ii<(ns_config->num_v_vars); ii++) {
			fprintf(OUT,"synth%d",ii+1);
			(ii==ns_config->num_v_vars-1) ? fprintf(OUT,";\n") : fprintf(OUT,", ");
		}
	}
	else if (0)
	{
		fprintf(OUT,"variable_group ");
		for (ii=0;ii<numVars; ii++){
			fprintf(OUT,"%s, ",vars[ii]);
			//			(ii==numVars-1) ? fprintf(OUT,";\n") : fprintf(OUT,", ");
		}
		
		for (ii=0; ii<(ns_config->num_v_vars); ii++) {
			fprintf(OUT,"synth%d",ii+1);
			(ii==ns_config->num_v_vars-1) ? fprintf(OUT,";\n") : fprintf(OUT,", ");
		}
	}
	else
	{
		fprintf(OUT,"variable_group ");
		for (ii=0;ii<numVars; ii++){
			fprintf(OUT,"%s",vars[ii]);
			(ii==numVars-1) ? fprintf(OUT,";\n") : fprintf(OUT,", ");
		}
	}
	
	
	for (ii=0; ii<ns_config->num_projections; ii++) {
		std::stringstream projname;
		projname << "pi_" << ii+1;
		write_vector_as_constants(ns_config->target_projection[ii], projname.str(), OUT);
	}
	
	if (!ns_config->randomizer()->is_square()) {
		write_matrix_as_constants( *(ns_config->randomizer()->get_mat_full_prec()), "r", OUT);
	}
	
	WaitOnGeneratedFile("derivative_polynomials_declaration");
	
	IN = safe_fopen_read("derivative_polynomials_declaration");
	while ((ch = fgetc(IN)) != EOF)
		fprintf(OUT, "%c", ch);
	fclose(IN);
	
	// END; written in the above transciption
	fprintf(OUT,"END;\n\n\n");
	
	write_vector_as_constants(ns_config->v_patch, "synthetic_var_patch", OUT);
	fprintf(OUT,"function synth_var_patch;\n");
	fprintf(OUT,"synth_var_patch = ");
	for (ii=0; ii<(ns_config->num_v_vars); ii++) {
		fprintf(OUT,"synthetic_var_patch_%d * synth%d",ii+1,ii+1);
		(ii==ns_config->num_v_vars-1) ? fprintf(OUT,";\n") : fprintf(OUT," + ");
	}
	
	
	fprintf(OUT,"\n\nTODO: add natural variable patch");
	
	fclose(OUT);
	// clear memory
	free(str);
	
	for (ii = 0; ii < numVars; ii++)
		free(vars[ii]);
	free(vars);
	free(lineVars);
	for (ii = 0; ii < numConstants; ii++)
		free(consts[ii]);
	free(consts);
	free(lineConstants);
	for (ii = 0; ii < numFuncs; ii++)
		free(funcs[ii]);
	free(funcs);
	free(lineFuncs);
	
	return;
}

void createMatlabDerivative(boost::filesystem::path output_name,
							boost::filesystem::path input_name,
							NullspaceConfiguration *ns_config,
							int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs)
/***************************************************************\
 * USAGE: setup a Matlab script to perform the deflation         *
 * ARGUMENTS: input file, declaration name, and number of lines  *
 * RETURN VALUES: Matlab script                                  *
 * NOTES:                                                        *
 \***************************************************************/
{
	int ii, lineNumber = 1, cont = 1, declares = 0, strSize = 1;
	char *str = (char *)br_malloc(strSize * sizeof(char));
	
	std::ifstream IN(input_name.c_str());
	std::ofstream OUT(output_name.c_str());
	// setup Bertini constants in Matlab
	OUT << "syms I Pi;\nI = 1i;\n";
	
	// setup variables
	OUT << "syms";
	for (ii = 0; ii < numVars; ii++)
		OUT << " " << vars[ii];
	OUT << ";\nmy_var_names = [";
	for (ii = 0; ii < numVars; ii++)
		OUT << " " << vars[ii];
	OUT << "];\n\n";
	
	OUT << "syms";
	for (ii = 0; ii < ns_config->num_v_vars; ii++)
		OUT << " synth" << ii+1;
	
	OUT << ";\nsynth_vars = [";
	for (ii = 0; ii < ns_config->num_v_vars; ii++)
		OUT << " " << "synth" << ii+1 << ";";
	OUT << "];\n\n";
	
	
	
	OUT << "syms";
	for (ii = 1; ii <= ns_config->num_projections; ii++)
		for (int jj = 1; jj< ns_config->num_natural_vars; jj++)
			OUT << " pi_" << ii << "_" << jj+1;
	OUT << "\n";
	
	OUT << "proj = [";
	for (ii = 1; ii <= ns_config->num_projections; ii++){
		OUT << "[";
		for (int jj = 1; jj< ns_config->num_natural_vars; jj++){
			OUT << " pi_" << ii << "_" << jj+1 << ";";
		}
		OUT << " ] ";
	}
	OUT << "];\n\n";
	
	// setup constants
	if (numConstants > 0)
	{
		OUT << "syms";
		for (ii = 0; ii < numConstants; ii++)
			OUT << " " << consts[ii];
		OUT << ";\n";
	}
	
	
	OUT << "num_jac_equations = " << ns_config->num_jac_equations << ";\n";
	OUT << "num_randomized_eqns = " << ns_config->randomizer()->num_rand_funcs() << ";\n";
	OUT << "target_crit_codim = " << ns_config->target_crit_codim << ";\n";
	
	// copy lines which do not declare items or define constants (keep these as symbolic objects)
	while (cont)
	{ // see if this line number declares items
		declares = 0;
		for (ii = 0; ii < numVars; ii++)
			if (lineNumber == lineVars[ii])
				declares = 1;
		for (ii = 0; ii < numConstants; ii++)
			if (lineNumber == lineConstants[ii])
				declares = 1;
		for (ii = 0; ii < numFuncs; ii++)
			if (lineNumber == lineFuncs[ii])
				declares = 1;
		
		std::string current_line;
		getline(IN,current_line);
		//		std::cout << "curr line: " << current_line << std::endl;
		
		if (declares)
		{ // move past this line
			
		}
		else
		{ // check to see if this defines a constant - line must be of the form '[NAME]=[EXPRESSION];' OR EOF
			
			std::string name;
			
			size_t found=current_line.find('=');
			if (found!=std::string::npos) {
				name = current_line.substr(0,found);
			}
			
			declares = 0;
			// compare against constants
			for (ii = 0; ii < numConstants; ii++)
				if (strcmp(name.c_str(), consts[ii]) == 0)
					declares = 1;
			
			if (!declares)// print line
				OUT << current_line << std::endl;
		}
		
		// increment lineNumber
		lineNumber++;
		
		// test for EOF
		if (IN.eof())
			cont = 0;
	}
	
	OUT << "syms ";
	for (ii=0; ii<ns_config->randomizer()->num_rand_funcs(); ii++) {
		for (int jj=0; jj<ns_config->randomizer()->num_base_funcs(); jj++) {
			OUT << "r_" << ii+1 << "_" << jj+1 << " "	;
		}
	}
	OUT << ";\n";
	
	
	
	
	
	
	
	// setup functions
	OUT << "\nF_orig = [";
	for (ii = 0; ii < numFuncs; ii++)
		OUT << " " << funcs[ii] << ";";
	OUT << "]; %collect the functions into a single matrix\n";
	
	// put in the randomization matrices.
	if (!ns_config->randomizer()->is_square()) {
		OUT << "R = [";
		for (ii=0; ii<ns_config->randomizer()->num_rand_funcs(); ii++) {
			for (int jj=0; jj<ns_config->randomizer()->num_base_funcs(); jj++) {
				OUT << "r_" << ii+1 << "_" << jj+1 << " "	;
			}
			OUT << ";\n";
		}
		OUT << "];\n\n";
		
		OUT << "F_rand = R*F_orig; % randomize\n";
	}
	else{
		
		OUT << "F_rand = F_orig; % no need to randomize\n";
	}
	
	
	// compute the jacobian
	OUT << "J = [transpose(jacobian(F_rand,my_var_names)) proj];  %compute the transpose of the jacobian\n";
	OUT << "                                           %concatenate the projections\n\n";
	
	OUT << "new_eqns = J*synth_vars; % multiply the synth vars\n\n";
	
	OUT << "OUT = fopen('derivative_polynomials_declaration','w'); %open the file to write\n";
	
	OUT << "fprintf(OUT, 'function ');\n";
	OUT << "for ii=1:num_randomized_eqns\n";
	OUT << "  fprintf(OUT, 'f%i',ii);\n";
	OUT << "  if ii~=num_randomized_eqns\n";
	OUT << "    fprintf(OUT,', ');\n";
	OUT << "  else\n";
	OUT << "    fprintf(OUT,';\\n');\n";
	OUT << "  end %re: if\n";
	OUT << "end\n\n";
	
	OUT << "for ii=1:num_randomized_eqns\n";
	OUT << "  fprintf(OUT,'f%i = %s;\\n',ii,char(F_rand(ii)));\n";
	OUT << "end\n\n";
	
	OUT << "fprintf(OUT, 'function ');\n";
	OUT << "for jj = 1:num_jac_equations\n";
	OUT << "  fprintf(OUT, 'der_func_%i', jj);\n";
	OUT << "  if jj~=num_jac_equations\n    fprintf(OUT,', ');\n  end %re: if\n";
	OUT << "end %re: jj\n";
	OUT << "fprintf(OUT,';\\n');\n\n";
	
	OUT << "for jj = 1:num_jac_equations\n";
	
	OUT << "	curr_eqn = char(new_eqns(jj));\n";
	OUT << "	i_locations = regexp(curr_eqn,'[\\W\\s]i[\\W\\s]');\n";
	OUT << "	curr_eqn(i_locations+1) = 'I';\n";
	OUT << "  fprintf(OUT, 'der_func_%i = %s;  \\n',jj,curr_eqn);\n";
	//	OUT << "  for kk = 1:num_jac_equations\n";
	//	OUT << "    fprintf(OUT,'%s\n',);\n";
	//	OUT << "  end %re: kk\n";
	//	OUT << "  for kk = 1:target_crit_codim\n";
	//	OUT << "    fprintf(OUT,'pi_%i_%i*synth%i',kk,jj+1,kk+num_randomized_eqns);\n";
	//	OUT << "    if kk~=target_crit_codim\n";
	//	OUT << "      fprintf(OUT,' + ');\n";
	//	OUT << "    else\n";
	//	OUT << "      fprintf(OUT,';\\n');\n";
	//	OUT << "    end %re: if\n";
	//	OUT << "  end %re: kk\n";
	OUT << "end %re: jj\n\n";
	OUT << "\nexit %exit the script\n";
	
	
	
	// clear memory
	free(str);
	
	return;
}














void create_matlab_determinantal_system(boost::filesystem::path output_name,
										boost::filesystem::path input_name,
										NullspaceConfiguration *ns_config,
										int numVars, char **vars, int *lineVars, int numConstants, char **consts, int *lineConstants, int numFuncs, char **funcs, int *lineFuncs)
/***************************************************************\
 * USAGE: setup a Matlab script to perform the deflation         *
 * ARGUMENTS: input file, declaration name, and number of lines  *
 * RETURN VALUES: Matlab script                                  *
 * NOTES:                                                        *
 \***************************************************************/
{
	int ii, lineNumber = 1, declares = 0, strSize = 1;
	char *str = (char *)br_malloc(strSize * sizeof(char));
	
	std::ifstream IN(input_name.c_str());
	std::ofstream OUT(output_name.c_str());
	// setup Bertini constants in Matlab
	OUT << "syms I Pi;\n";
	
	// setup variables
	OUT << "syms";
	for (ii = 0; ii < numVars; ii++)
		OUT << " " << vars[ii];
	OUT << ";\nmy_var_names = [";
	for (ii = 0; ii < numVars; ii++)
		OUT << " " << vars[ii];
	OUT << "];\n\n";
	
	
	
	
	
	OUT << "syms";
	for (ii = 1; ii <= ns_config->num_projections; ii++)
		for (int jj = 1; jj< ns_config->num_natural_vars; jj++)
			OUT << " pi_" << ii << "_" << jj+1;
	OUT << "\n";
	
	OUT << "proj = [";
	for (ii = 1; ii <= ns_config->num_projections; ii++){
		OUT << "[";
		for (int jj = 1; jj< ns_config->num_natural_vars; jj++){
			OUT << " pi_" << ii << "_" << jj+1 << "";
		}
		OUT << " ]; ";
	}
	OUT << "];\n\n";
	
	// setup constants
	if (numConstants > 0)
	{
		OUT << "syms";
		for (ii = 0; ii < numConstants; ii++)
			OUT << " " << consts[ii];
		OUT << ";\n";
	}
	
	OUT << "num_projections = " << ns_config->num_projections << ";\n";
	OUT << "num_jac_equations = " << ns_config->num_jac_equations << ";\n";
	OUT << "num_randomized_eqns = " << ns_config->randomizer()->num_rand_funcs() << ";\n";
	OUT << "target_crit_codim = " << ns_config->target_crit_codim << ";\n";
	
	int cont = 1;
	// copy lines which do not declare items or define constants (keep these as symbolic objects)
	while (cont)
	{ // see if this line number declares items
		declares = 0;
		for (ii = 0; ii < numVars; ii++)
			if (lineNumber == lineVars[ii])
				declares = 1;
		for (ii = 0; ii < numConstants; ii++)
			if (lineNumber == lineConstants[ii])
				declares = 1;
		for (ii = 0; ii < numFuncs; ii++)
			if (lineNumber == lineFuncs[ii])
				declares = 1;
		
		std::string current_line;
		getline(IN,current_line);
		//		std::cout << "curr line: " << current_line << std::endl;
		
		if (declares)
		{ // move past this line
			
		}
		else
		{ // check to see if this defines a constant - line must be of the form '[NAME]=[EXPRESSION];' OR EOF
			
			std::string name;
			
			size_t found=current_line.find('=');
			if (found!=std::string::npos) {
				name = current_line.substr(0,found);
			}
			
			declares = 0;
			// compare against constants
			for (ii = 0; ii < numConstants; ii++)
				if (strcmp(name.c_str(), consts[ii]) == 0)
					declares = 1;
			
			if (!declares)// print line
				OUT << current_line << std::endl;
		}
		
		// increment lineNumber
		lineNumber++;
		
		// test for EOF
		if (IN.eof())
			cont = 0;
	}
	
	OUT << "syms ";
	for (ii=0; ii<ns_config->randomizer()->num_rand_funcs(); ii++) {
		for (int jj=0; jj<ns_config->randomizer()->num_base_funcs(); jj++) {
			OUT << "r_" << ii+1 << "_" << jj+1 << " "	;
		}
	}
	OUT << ";\n";
	
	
	
	
	
	
	
	// setup functions
	OUT << "\nF_orig = [";
	for (ii = 0; ii < numFuncs; ii++)
		OUT << " " << funcs[ii] << ";";
	OUT << "]; %collect the functions into a single matrix\n";
	
	// put in the randomization matrices.
	if (!ns_config->randomizer()->is_square()) {
		OUT << "R = [";
		for (ii=0; ii<ns_config->randomizer()->num_rand_funcs(); ii++) {
			for (int jj=0; jj<ns_config->randomizer()->num_base_funcs(); jj++) {
				OUT << "r_" << ii+1 << "_" << jj+1 << " "	;
			}
			OUT << ";\n";
		}
		OUT << "];\n\n";
		
		OUT << "F_rand = R*F_orig; % randomize\n";
	}
	else{
		
		OUT << "F_rand = F_orig; % no need to randomize\n";
	}
	
	
	// compute the jacobian
	OUT << "J = [jacobian(F_rand,my_var_names); proj];  %compute the transpose of the jacobian\n";
	OUT << "                                           %concatenate the projections\n\n";
	
	OUT << "new_eqn = det(J); \n\n";
	
	OUT << "OUT = fopen('derivative_polynomials_declaration','w'); %open the file to write\n";
	
	OUT << "fprintf(OUT, 'function ');\n";
	OUT << "for ii=1:num_randomized_eqns\n";
	OUT << "  fprintf(OUT, 'f%i',ii);\n";
	OUT << "  if ii~=num_randomized_eqns\n";
	OUT << "    fprintf(OUT,', ');\n";
	OUT << "  else\n";
	OUT << "    fprintf(OUT,';\\n');\n";
	OUT << "  end %re: if\n";
	OUT << "end\n\n";
	
	OUT << "for ii=1:num_randomized_eqns\n";
	OUT << "  fprintf(OUT,'f%i = %s;\\n',ii,char(F_rand(ii)));\n";
	OUT << "end\n\n";
	
	OUT << "fprintf(OUT, 'function der_func;\\n');\n";
	//	OUT << "for jj = 1:num_jac_equations\n";
	//	OUT << "  fprintf(OUT, 'der_func_%i', jj);\n";
	//	OUT << "  if jj~=num_jac_equations\n    fprintf(OUT,', ');\n  end %re: if\n";
	//	OUT << "end %re: jj\n";
	//	OUT << "fprintf(OUT,';\\n');\n\n";
	
	
	OUT << "fprintf(OUT,'der_func = %s;\\n',char(new_eqn));\n\n";
	//	OUT << "for jj = 1:num_jac_equations\n";
	//	OUT << "  fprintf(OUT, 'der_func_%i = %s;  \\n',jj,char(new_eqns(jj)));\n";
	//	OUT << "  for kk = 1:num_jac_equations\n";
	//	OUT << "    fprintf(OUT,'%s\n',);\n";
	//	OUT << "  end %re: kk\n";
	//	OUT << "  for kk = 1:target_crit_codim\n";
	//	OUT << "    fprintf(OUT,'pi_%i_%i*synth%i',kk,jj+1,kk+num_randomized_eqns);\n";
	//	OUT << "    if kk~=target_crit_codim\n";
	//	OUT << "      fprintf(OUT,' + ');\n";
	//	OUT << "    else\n";
	//	OUT << "      fprintf(OUT,';\\n');\n";
	//	OUT << "    end %re: if\n";
	//	OUT << "  end %re: kk\n";
	//	OUT << "end %re: jj\n\n";
	
	OUT << "exit %quit the script\n";
	
	
	// clear memory
	free(str);
	
	return;
}
















