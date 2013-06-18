#include "nullspace_left.hpp"


int compute_crit_nullspace(witness_set *W_crit_real, // the returned value
													 witness_set W,
													 mat_mp randomizer_matrix,
													 vec_mp *pi, // an array of projections, the number of which is the target dimensions
													 int *randomized_degrees, // an array of integers holding the degrees of the randomized equations.
													 int ambient_dim,
													 int target_dim,
													 int target_crit_codim,
													 program_configuration *program_options,
													 solver_configuration *solve_options)
{
	//many of the 1's here should be replaced by the number of patch equations, or the number of variable_groups
	
	int ii, jj;
	int offset;
	witness_set Wtemp, Wtemp2;
	
	nullspace_config ns_config;
	
	// get the max degree of the derivative functions.  this is unique to the left formulation
	int max_degree;
	
	
	nullspace_config_setup(&ns_config,
												 pi,
												 ambient_dim,
												 target_dim,
												 target_crit_codim,
												 &max_degree,
												 randomized_degrees,
												 randomizer_matrix,
												 W,
												 solve_options);
	
	
	
	printf("max_degree: %d\nnum_jac_equations: %d\n",max_degree,ns_config.num_jac_equations);
	

	//
	///
	/////        end setup
	////////
	//////////////
	///////////////////////////
	
	
	

	
	
	
	
	
	//  2.  Do a bunch of homotopies in $x$, each set of which will be followed by a single linear solve in $v$.
	
	
	// setup for the multilin moves
	//  these are for feeding into the multilin solver -- and that's it.  they'll be set in the while loop
	vec_mp *multilin_linears = (vec_mp *) br_malloc(ambient_dim*sizeof(vec_mp)); // target dim is the number of linears in the input witness set
	for (ii=0; ii<ambient_dim; ii++) {
		init_vec_mp(multilin_linears[ii],W.num_variables); multilin_linears[ii]->size = W.num_variables;
		vec_cp_mp(multilin_linears[ii], W.L_mp[ii]);
	}

	
	int num_funcs_atatime = target_dim; // subtract 1 for the patch equation
	
	// create and seed the function indices -- keep track of which functions we are working on
	int *function_indices = NULL;
	function_indices = (int *)br_malloc( num_funcs_atatime*sizeof(int) );
	for (ii=0; ii<num_funcs_atatime; ii++) {
		function_indices[ii] = ii;
	}
	
	
	//create and seed the unused_function_indices -- keep track of the linears we are NOT using
	int *unused_function_indices = NULL;
	unused_function_indices = (int *)br_malloc( (ns_config.num_jac_equations - num_funcs_atatime) *sizeof(int));
	for (ii=num_funcs_atatime; ii<ns_config.num_jac_equations; ii++) {
		unused_function_indices[ii-num_funcs_atatime] = ii; // shift index for assignment to start at 0
	}
	
	int *degrees_for_counter = (int *) br_malloc(ns_config.num_jac_equations * sizeof(int));
	for (ii=0; ii<ns_config.num_jac_equations; ii++)
		degrees_for_counter[ii] = max_degree;
	
	
	//create and seed the subindices -- keep track of the specific linears we are working on.
	int *subindices = (int *)br_malloc(num_funcs_atatime*sizeof(int));
	memset(subindices, 0, (num_funcs_atatime)*sizeof(int)); // initialize to 0
	
	
	
	
	
	
	// this is for performing the matrix inversion to get ahold of the $v$ values corresponding to $x$
	mat_mp tempmat;  init_mat_mp(tempmat, ns_config.num_v_vars, ns_config.num_v_vars);
	tempmat->rows = tempmat->cols = ns_config.num_v_vars;
	
	offset = ns_config.num_v_vars-1;
	for (jj=0; jj<ns_config.num_v_vars; jj++)
		set_mp(&tempmat->entry[offset][jj], &ns_config.v_patch->coord[jj]);
	
	// for holding the result of the matrix inversion
	vec_mp result; init_vec_mp(result,ns_config.num_v_vars);  result->size = ns_config.num_v_vars;
	
	// use this for the matrix inversion
	vec_mp invert_wrt_me;
	init_vec_mp(invert_wrt_me,ns_config.num_v_vars); invert_wrt_me->size = ns_config.num_v_vars;
	for (ii=0; ii<ns_config.num_v_vars-1; ii++)
		set_zero_mp(&invert_wrt_me->coord[ii]); // set zero
	
	set_one_mp(&invert_wrt_me->coord[ns_config.num_v_vars-1]); // the last entry is set to 1 for the patch equation
	
	
	
	
	
	
	int ticked_new_functions; // for telling whether we are done with a set of functions.  should be a bool.
	
	
	
	vec_mp temppoint;  init_vec_mp(temppoint, ns_config.num_x_vars + ns_config.num_v_vars);
	temppoint->size = ns_config.num_x_vars + ns_config.num_v_vars;
	
	
	witness_set W_step_one;  init_witness_set(&W_step_one);
	W_step_one.num_variables = W.num_variables;
	cp_patches(&W_step_one, W);
	cp_names(&W_step_one, W);
	
	
	witness_set W_linprod;  init_witness_set(&W_linprod);
	W_linprod.num_variables = ns_config.num_x_vars + ns_config.num_v_vars;
	
//	int current_absolute_index = 0; // the while loop termination condition\//current_absolute_index<terminationint
	int reached_end = 0;
	while (!reached_end) { // current_absolute_index incremented at the bottom of loop
		
		if (program_options->verbose_level>=2)
		{  //some display
			printf("functions: ");
			for (ii=0; ii<num_funcs_atatime; ii++)
				printf("%d ", function_indices[ii]);
			printf("\t\tsubindices: ");
			for (ii=0; ii<num_funcs_atatime; ii++)
				printf("%d ", subindices[ii]);
			printf("\t\tunused_functions: ");
			for (ii=0; ii<ns_config.num_jac_equations - num_funcs_atatime; ii++)
				printf("%d ", unused_function_indices[ii]);
			
			printf("\n");
			
			for (ii=0; ii<ambient_dim; ii++) 
				print_point_to_screen_matlab_mp(multilin_linears[ii],"newlin");
		}
		
		
		//copy in the linears for the solve
		for (ii=0; ii<target_dim; ii++) {
			vec_cp_mp(multilin_linears[ii], ns_config.starting_linears[function_indices[ii]][subindices[ii]]);
		}
		// the remainder of the linears are left alone (stay stationary).
		
		
		{
			solve_options->allow_singular = 1;
			solve_options->complete_witness_set = 1;
			init_witness_set(&Wtemp); // intialize for holding the data
			// actually solve WRT the linears
			multilintolin_solver_main(W.MPType,
																W,         // witness_set
																randomizer_matrix,
																multilin_linears, //  the set of linears we will solve at.
																&Wtemp, // the new data is put here!
																solve_options); // already a pointer
			

			// merge into previously obtained data
			init_witness_set(&Wtemp2);
			cp_witness_set(&Wtemp2, W_step_one);   // i dislike this pattern of copy and merge.  it seems wasteful
			clear_witness_set(W_step_one); init_witness_set(&W_step_one);
			
			
			merge_witness_sets(&W_step_one, Wtemp2, Wtemp);
			clear_witness_set(Wtemp);  clear_witness_set(Wtemp2);
		}
		
		
		// increment the tracking indices
		ticked_new_functions = increment_subindices(&function_indices, &subindices,
																								degrees_for_counter, num_funcs_atatime,
																								ns_config.num_jac_equations);
		
		// if it's time to move on to next function combo, must commit the points
		if (ticked_new_functions==1)
		{
			//set the v_linears (M_i).
			for (ii=0; ii<ns_config.num_v_vars-1; ii++) // deficient one because of the patch equation
				for (jj=0; jj<ns_config.num_v_vars; jj++)
					set_mp(&tempmat->entry[ii][jj], &ns_config.v_linears[unused_function_indices[ii]]->coord[jj]);
			
			// invert the matrix for the v variables.
			matrixSolve_LU_mp(result, tempmat,  invert_wrt_me, 1e-14, 1e9);
					
			
			
			
			
			offset = ns_config.num_x_vars;
			for (jj=0; jj<ns_config.num_v_vars; jj++) 
				set_mp(&temppoint->coord[jj+offset], &result->coord[jj]);
			
			
			for (ii=0; ii<W_step_one.num_pts; ii++) {
				for (jj=0; jj<ns_config.num_x_vars; jj++) {
					set_mp(&temppoint->coord[jj], &W_step_one.pts_mp[ii]->coord[jj]);
				}
				add_point_to_witness_set(&W_linprod, temppoint);
			}
			
			
			
			
			clear_witness_set(W_step_one);
			init_witness_set(&W_step_one);
			W_step_one.num_variables = W.num_variables;
			cp_patches(&W_step_one, W);  // necessary?
			cp_names(&W_step_one, W); // necessary?
			
			
			
			reached_end = increment_function_indices(&function_indices, &unused_function_indices,
																							 num_funcs_atatime,
																							 ns_config.num_jac_equations);
		}
	}
	
	
	clear_vec_mp(temppoint);
	
	
	cp_patches(&W_linprod, W);
	cp_names(&W_linprod, W);
	
	//set some solver options
	
	solve_options->complete_witness_set = 0;
	solve_options->allow_multiplicity = 1;
	solve_options->allow_singular = 1;
	solve_options->use_midpoint_checker = 0;
	
	
	printf("entering nullspace solver\n");
	init_witness_set(&Wtemp);
	nullspacejac_solver_main(W.MPType,
													 W_linprod, // carries with it the start points, but not the linears.
													 &Wtemp,   // the created data goes in here.
													 &ns_config,
													 solve_options);
	

	
	
	
	init_witness_set(&Wtemp2);
	Wtemp2.num_variables = W.num_variables;
	cp_patches(&Wtemp2, W);
	cp_names(&Wtemp2, W);
	vec_mp tempvec;  init_vec_mp(tempvec, W.num_variables);  tempvec->size = W.num_variables;
	for (ii=0; ii<Wtemp.num_pts; ii++) {
		
		for (jj=0; jj<W.num_variables; jj++) {
			set_mp(&tempvec->coord[jj], &Wtemp.pts_mp[ii]->coord[jj]);
		}
		add_point_to_witness_set(&Wtemp2,tempvec);
	}

	clear_witness_set(Wtemp);
	
	
	init_witness_set(&Wtemp);
	
	
	sort_for_real(&Wtemp, Wtemp2,solve_options->T); // get only the real solutions.
	clear_witness_set(Wtemp2);

	W_crit_real->num_variables = W.num_variables;
	cp_patches(W_crit_real, W);
	cp_names(W_crit_real, W);
	sort_for_unique(W_crit_real, Wtemp,solve_options->T); // get only the real solutions.
	
	clear_witness_set(Wtemp);
	

	
	
	
	
	clear_mat_mp(tempmat);
	clear_vec_mp(invert_wrt_me);
	clear_vec_mp(result);
	
	
	
	
	
	
	
	
	return SUCCESSFUL;
}










int increment_subindices(int **function_indices,
												 int **subindices,
												 int *degrees,
												 int num_inner_indices, //these should be implicitly available.  switch to c++ plx.
												 int num_functions)
{
	
	
	int ii;
	int carry = 1; // seed carry so that it causes addition of at least the last entry of the odometer
	for (ii=num_inner_indices-1; ii>=0; ii--) { // count down from the end of the indexes
		
		if (carry==1)
			(*subindices)[ii]++;
		
		if ( (*subindices)[ii]>=(degrees[ (*function_indices)[ii]]) ) {
			(*subindices)[ii] = 0;
			carry = 1;
		}
		else{
			carry = 0;
			break;
		}
		
		
	}
	return  carry;  // if return 1, then need to increment the functions.
}








int increment_function_indices(int **function_indices,
															 int **unused_function_indices,
															 int num_inner_indices,
															 int num_functions)
{
	int ii, jj;
	int carry = 1; // seed 1
	// increment the function indices
	
	int local_counter = 1;
	for (ii=num_inner_indices-1; ii>=0; ii--) {
		
		if (carry==1){
			(*function_indices)[ii]++;
			for (jj=ii+1; jj<num_inner_indices; jj++) {
				(*function_indices)[jj] = (*function_indices)[jj]+1;
			}
		}
		
		if ( (*function_indices)[ii]>num_functions-local_counter)
		{
			carry = 1;
		}
		else{
			break; // all done!
		}
		
		local_counter++;
	}
	
	local_counter = 0;
	for (ii=0; ii<num_functions; ii++) {
		
		// a crappy find function
		int ok = 1;
		for (jj=0; jj<num_inner_indices; jj++) {
			if ((*function_indices)[jj]==ii) {
				ok=0;
			}
		}
		
		if ( (ok==1) && (num_functions != num_inner_indices)){ // if didn't find current index in the list, put it in the unused list.
			(*unused_function_indices)[local_counter] = ii;
			local_counter++;
		}
	}
	

	if ( (*function_indices)[num_inner_indices-1]>=num_functions) { // if the last index is greater than allowed
		return 1;
	}
	else{
		return 0;
	}
	
}



void nullspace_config_setup(nullspace_config *ns_config,
														vec_mp *pi, // an array of projections, the number of which is the target dimensions
														int ambient_dim,
														int target_dim,
														int target_crit_codim,
														int *max_degree, // a pointer to the value
														int *randomized_degrees, // an array of randomized degrees
														mat_mp randomizer_matrix,
														witness_set W,
														solver_configuration *solve_options)
{

	int ii, jj, kk;
	
	
	int maxiii = 0;
	ns_config->randomized_degrees = (int *)br_malloc(randomizer_matrix->rows * sizeof(int));
	for (ii=0; ii<randomizer_matrix->rows; ii++	) {
		ns_config->randomized_degrees[ii] = randomized_degrees[ii]; // store the full degree (not derivative).
		if ( (randomized_degrees[ii]-1) > maxiii)
			maxiii = randomized_degrees[ii]-1; // minus one for the derivative
	}
	*max_degree = maxiii;
	
	
	
	// set some integers
	ns_config->num_v_vars = W.num_variables-1 - ambient_dim + target_crit_codim; //  N-k+l
	ns_config->num_x_vars = W.num_variables;
	
	ns_config->ambient_dim = ambient_dim;
	ns_config->target_dim = target_dim;
	ns_config->target_crit_codim = target_crit_codim;
	ns_config->max_degree = (*max_degree);

	ns_config->target_projection = (vec_mp *) br_malloc(target_crit_codim * sizeof(vec_mp));
	for (ii=0; ii<target_crit_codim; ii++) {
		init_vec_mp2(ns_config->target_projection[ii], W.num_variables,solve_options->T.AMP_max_prec);
		ns_config->target_projection[ii]->size = W.num_variables;
		vec_cp_mp(ns_config->target_projection[ii], pi[ii]);
	}
	
	
	// make the post-randomizer matrix $S$
	int num_jac_equations = ns_config->num_v_vars + target_dim - 1;//N-k+l+r-1;  the subtraction of 1 is for the 1 hom-var
	
	init_mat_mp2(ns_config->post_randomizer_matrix,
							 num_jac_equations, W.num_variables-1,
							 solve_options->T.AMP_max_prec);
	
	ns_config->post_randomizer_matrix->rows = num_jac_equations;
	ns_config->post_randomizer_matrix->cols = W.num_variables-1;
	if (ns_config->post_randomizer_matrix->rows == ns_config->post_randomizer_matrix->cols) {
		make_matrix_ID_mp(ns_config->post_randomizer_matrix,num_jac_equations,num_jac_equations);
	}
	else{
	make_matrix_random_real_mp(ns_config->post_randomizer_matrix,
														 num_jac_equations,W.num_variables-1,
														 solve_options->T.AMP_max_prec);
	}
	ns_config->num_jac_equations = num_jac_equations;
	
	ns_config->num_randomized_eqns = W.num_variables-1-ambient_dim;
	
	ns_config->num_v_linears = ns_config->num_jac_equations;
	
	
	// set up the linears in $v$  ( the M_i linears)
	ns_config->v_linears = (vec_mp *)br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
	for (ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_mp2(ns_config->v_linears[ii],ns_config->num_v_vars,solve_options->T.AMP_max_prec);
		ns_config->v_linears[ii]->size = ns_config->num_v_vars;
		for (jj=0; jj<ns_config->num_v_vars; jj++){
			get_comp_rand_mp(&ns_config->v_linears[ii]->coord[jj]); // should this be real?
		}
	}
	
	
	
	int offset = target_dim;  //  CHECK ME!!!
	ns_config->num_additional_linears = ambient_dim - target_dim;
	ns_config->additional_linears_terminal = (vec_mp *)br_malloc((ns_config->num_additional_linears)*sizeof(vec_mp));
	ns_config->additional_linears_starting = (vec_mp *)br_malloc((ns_config->num_additional_linears)*sizeof(vec_mp));
	for (ii=0; ii<ns_config->num_additional_linears; ii++) {
		init_vec_mp2(ns_config->additional_linears_terminal[ii],W.num_variables,solve_options->T.AMP_max_prec);
		ns_config->additional_linears_terminal[ii]->size = W.num_variables;
		for (jj=0; jj<W.num_variables; jj++){
			get_comp_rand_real_mp(&ns_config->additional_linears_terminal[ii]->coord[jj]); // should this be real?
		}
		
		
		init_vec_mp2(ns_config->additional_linears_starting[ii],W.num_variables,solve_options->T.AMP_max_prec);
		ns_config->additional_linears_starting[ii]->size = W.num_variables;
		vec_cp_mp(ns_config->additional_linears_starting[ii], W.L_mp[ii+offset]);
	}
	
	
	// set up the patch in $v$.  we will include this in an inversion matrix to get the starting $v$ values.
	init_vec_mp2(ns_config->v_patch,ns_config->num_v_vars,solve_options->T.AMP_max_prec);
	ns_config->v_patch->size = ns_config->num_v_vars;
	for (ii=0; ii<ns_config->num_v_vars; ii++) {
		get_comp_rand_mp(&ns_config->v_patch->coord[ii]);
	}
	
	
	
	
	//the 'ns_config->starting_linears' will be used for the x variables.  we will homotope to these $r$ at a time
	ns_config->starting_linears = (vec_mp **)br_malloc( num_jac_equations*sizeof(vec_mp *));
	for (ii=0; ii<num_jac_equations; ii++) {
		ns_config->starting_linears[ii] = (vec_mp *) br_malloc((*max_degree)*sizeof(vec_mp)); //subtract 1 for differentiation
		for (jj=0; jj<(*max_degree); jj++) {
			init_vec_mp2(ns_config->starting_linears[ii][jj],W.num_variables,solve_options->T.AMP_max_prec);
			ns_config->starting_linears[ii][jj]->size = W.num_variables;
			//			set_zero_mp(&ns_config->starting_linears[ii][jj]->coord[0]);  // maybe? but prolly not
			for (kk=0; kk<W.num_variables; kk++) {
				get_comp_rand_mp(&ns_config->starting_linears[ii][jj]->coord[kk]);
			}
		}
	}
	
	// copy the main randomizer matrix
	init_mat_mp2(ns_config->randomizer_matrix,randomizer_matrix->rows, randomizer_matrix->cols,solve_options->T.AMP_max_prec);
	ns_config->randomizer_matrix->rows = randomizer_matrix->rows;
	ns_config->randomizer_matrix->cols = randomizer_matrix->cols;
	mat_cp_mp(ns_config->randomizer_matrix, randomizer_matrix);

	

	
	return;
}












