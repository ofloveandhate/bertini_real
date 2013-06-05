#include "nullspace.h"


int compute_crit_nullspace(witness_set *W_crit_real, // the returned value
													 witness_set W,
													 mat_mp randomizer_matrix,
													 vec_mp *random_complex_projections,
													 vec_mp *pi, // an array of projections, the number of which is the target dimensions
													 int *randomized_degrees, // an array of integers holding the degrees of the randomized equations.
													 int target_dim,
													 program_configuration *program_options,
													 solver_configuration *solve_options)
{
	//many of the 1's here should be replaced by the number of patch equations, or the number of variable_groups
	
	
	int ii, jj, kk;
	witness_set Wtemp, Wtemp2;
	
	nullspace_config ns_config;
	
	nullspace_config_setup(&ns_config,
												 target_dim,
												 randomized_degrees,
												 randomizer_matrix,
												 W);
		//
	///
	/////        end setup
	////////
	//////////////
	///////////////////////////
	
	
	
	
	
	//  1.  find witness set for
	//      \[  [ f(x) \\ L_1(x) \\ \vdots \\ L_r(x) ] \]
	//
	//      and extend all witness points to the new variables $v$ via a single linear solve.
	
	
	// A:  is this not given to us with the problem statement?
	
	
	
	
	
	
	
	
	
	
	//  2.  Use (1) above to do a bunch of homotopies in $x$ each set of which will be followed by a single linear solve in $v$.
	
	//  these are for feeding into the multilin solver -- and that's it.  they'll be set in the while loop
	vec_mp *multilin_linears = br_malloc(target_dim*sizeof(vec_mp));
	for (ii=0; ii<target_dim; ii++) {
		init_vec_mp(multilin_linears[ii],W.num_variables);
		multilin_linears[ii]->size = W.num_variables;
	}

	// this is for terminating the while loop used for building up the start system.
	int terminationint = 1;
	for (ii=0; ii<randomizer_matrix->rows; ii++) 
		terminationint *= MAX((randomized_degrees[ii])-1,1);
	
	printf("%d terminationint\n",terminationint);
	
	// create and seed the function indices -- keep track of which functions we are working on
	int *function_indices = NULL;
	function_indices = (int *)br_malloc( target_dim*sizeof(int));
	
	
	
	//create and seed the unused_function_indices -- keep track of the linears we are NOT using
	int *unused_function_indices = NULL;
	unused_function_indices = (int *)br_malloc( (ns_config.num_v_linears - target_dim) *sizeof(int));
	for (ii=target_dim; ii<randomizer_matrix->rows; ii++) {
		unused_function_indices[ii-target_dim] = ii;
	}
	
	
	//create and seed the subindices -- keep track of the specific linears we are working on.
	int *subindices = (int *)br_malloc(target_dim*sizeof(int));
	memset(subindices, 0, (target_dim)*sizeof(int)); // initialize to 0
	
	
	
	
	
	
	// this is for performing the matrix inversion to get ahold of the $v$ values corresponding to $x$
	mat_mp tempmat;  init_mat_mp(tempmat, W.num_variables-1, W.num_variables-1);
	tempmat->rows = tempmat->cols = W.num_variables-1;
	
	// for holding the result of the matrix inversion
	vec_mp result; init_vec_mp(result,W.num_variables-1);  result->size = W.num_variables-1;
	
	// use this for the matrix inversion
	vec_mp invert_wrt_me;
	init_vec_mp(invert_wrt_me,W.num_variables-1); invert_wrt_me->size = W.num_variables-1;
	for (ii=0; ii<W.num_variables-2; ii++) 
		set_zero_mp(&invert_wrt_me->coord[ii]); // set zero
	
	set_one_mp(&invert_wrt_me->coord[W.num_variables-2]); // the last entry is 1 for the patch equation
	
	
	
	
	int offset;
	int current_absolute_index = 0; // the while loop termination condition
	int ticked_new_functions; // for telling whether we are done with a set of functions.  should be a bool.
	
	
	
	witness_set W_step_one;  init_witness_set(&W_step_one);
	W_step_one.num_variables = W.num_variables;
	cp_patches(&W_step_one, W);
	cp_names(&W_step_one, W);
	
	
	witness_set W_linprod;  init_witness_set(&W_linprod);
	W_linprod.num_variables = W.num_variables + W.num_variables-1;
	
	
	while (current_absolute_index<terminationint) { // current_absolute_index incremented at the bottom of loop
		
		{
			//some display
			for (ii=0; ii<target_dim; ii++) {
				printf("functions: %d ", function_indices[ii]);
			}
			printf("\t\t");
			for (ii=0; ii<randomizer_matrix->rows - target_dim; ii++) {
				printf("unused_functions: %d ", unused_function_indices[ii]);
			}
			printf("\t\t");
			for (ii=0; ii<target_dim; ii++) {
				printf("subindices: %d ", subindices[ii]);
			}
			printf("\n");
			
			//copy in the linears for the solve
			for (ii=0; ii<target_dim; ii++) {
				vec_cp_mp(multilin_linears[ii], ns_config.starting_linears[function_indices[ii]][subindices[ii]]);
			}
		}
		
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
		cp_witness_set(&Wtemp2, W_step_one);
		merge_witness_sets(&W_step_one, Wtemp2, Wtemp);
		clear_witness_set(Wtemp);  clear_witness_set(Wtemp2);
		
		// increment the tracking indices
		ticked_new_functions = increment_subindices(&function_indices, &subindices, randomized_degrees, target_dim, randomizer_matrix->rows);
		
		// if it's time to move on to next function combo, must commit the points
		if (ticked_new_functions==1)
		{
			//perform the matrix inversion to obtain the v values
			
			
			// invert the matrix for the v variables.
			
			
			offset = 0;
			for (ii=0; ii<ns_config.num_v_linears - target_dim; ii++) 
				for (jj=0; jj<W.num_variables-1; jj++) 
					set_mp(&tempmat->entry[ii+offset][jj], &ns_config.v_linears[function_indices[ii]]->coord[jj]);
			
			
			offset = W.num_variables-target_dim-2;
			for (ii=0; ii<target_dim; ii++) 
				for (jj=0; jj<W.num_variables-1; jj++) 
					set_mp(&tempmat->entry[ii+offset][jj], &random_complex_projections[ii]->coord[jj+1]);
			
			
			offset = W.num_variables-target_dim-1;
			for (jj=0; jj<W.num_variables-1; jj++) 
				set_mp(&tempmat->entry[offset][jj], &ns_config.v_patch->coord[jj]);
			
			
			matrixSolve_LU_mp(result, tempmat,  invert_wrt_me, 1e-14, 1e9);
					
			
			vec_mp temppoint;  init_vec_mp(temppoint, W.num_variables+result->size);
			temppoint->size = W.num_variables + result->size;
			
			
			
			for (ii=0; ii<W_step_one.num_pts; ii++) {
				for (jj=0; jj<W.num_variables; jj++) {
					set_mp(&temppoint->coord[jj], &W_step_one.pts_mp[ii]->coord[jj]);
				}
				offset = W.num_variables;
				for (jj=0; jj<result->size; jj++) {
					set_mp(&temppoint->coord[jj+offset], &result->coord[jj]);
				}
				
				add_point_to_witness_set(&W_linprod, temppoint);
			}
			
			clear_vec_mp(temppoint);
			
			
			clear_witness_set(W_step_one);
			init_witness_set(&W_step_one);
			W_step_one.num_variables = W.num_variables;
			cp_patches(&W_step_one, W);
			cp_names(&W_step_one, W);
			
			
			
			increment_function_indices(&function_indices, &unused_function_indices, target_dim, randomizer_matrix->rows);
		}

		
		current_absolute_index++;
	}
	
//	print_witness_set_to_screen(W_linprod);
//	
//	mypause();
	
	
//TODO: set the config data
	
	printf("entering nullspace solver\n");
	init_witness_set(&Wtemp);
	nullspacejac_solver_main(W.MPType,
													 W_linprod, // carries with it the start points, but not the linears.
													 &Wtemp,   // the created data goes in here.
													 &ns_config,
													 solve_options);
	
	printf("done with nullspace solver\n");
//TODO: stuff with Wtemp here
	
	clear_witness_set(Wtemp);
	

	
	
	
	
	//  3.  (2) Gives the set of solutions to a polynomial system with products of linears of the appropriate degrees in the appropriate spots.  So, you just follow all of those points in a straight-line homotopy to this fairly complicated system.
	
	
	
	
	
	
	
	clear_mat_mp(tempmat);
	clear_vec_mp(invert_wrt_me);
	clear_vec_mp(result);
	
	
	
	
	
	
	return SUCCESSFUL;
}










int increment_subindices(int **function_indices,
											 int **subindices,
											 int * randomized_degrees,
											 int num_inner_indices, //these should be implicitly available.  switch to c++ plx.
											 int num_functions)
{
	
	
	int ii;
	int carry = 1; // seed carry so that it causes addition of at least the last entry of the odometer
	for (ii=num_inner_indices-1; ii>=0; ii--) { // count down from the end of the indexes
		
		if (carry==1)
			(*subindices)[ii]++;
		
		if ( (*subindices)[ii]>=(randomized_degrees[ (*function_indices)[ii]]-1) ) {
			(*subindices)[ii] = 0;
			carry = 1;
		}
		else{
			carry = 0;
			break;
		}
		
		
	}
	return  carry;
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
	

	
	
	return 0;
}



void nullspace_config_setup(nullspace_config *ns_config,
														int target_dim,
														int *randomized_degrees,
														mat_mp randomizer_matrix,
														witness_set W)
{

	int ii, jj, kk;
	
	ns_config->target_dim = target_dim;
	
	ns_config->jac_degrees = (int *)br_malloc(randomizer_matrix->rows*sizeof(int));
	for (ii=0; ii<randomizer_matrix->rows; ii++)
		ns_config->jac_degrees[ii] = randomized_degrees[ii]-1;
	
	
	printf("randomized degrees of jacobian functions:\n");
	for (ii=0; ii<randomizer_matrix->rows; ii++) {
		printf("%d ",ns_config->jac_degrees[ii]);
	}
	printf("\n");
	
	
	
	
	
	ns_config->num_v_linears = W.num_variables-1 - target_dim;
	
	
	// set up the linears in $v$
	ns_config->v_linears = br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
	for (ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_mp(ns_config->v_linears[ii],W.num_variables-1); ns_config->v_linears[ii]->size = W.num_variables-1;
		for (jj=0; jj<W.num_variables-1; jj++){
			get_comp_rand_mp(&ns_config->v_linears[ii]->coord[jj]); // should this be real?
		}
	}
	
	
	
	// set up the patch in $v$
	init_vec_mp(ns_config->v_patch,W.num_variables-1);  ns_config->v_patch->size = W.num_variables-1;
	for (ii=0; ii<W.num_variables-1; ii++) {
		get_comp_rand_mp(&ns_config->v_patch->coord[ii]);
	}
	
	
	
	
	//the 'ns_config->starting_linears' will be used for the x variables
	ns_config->starting_linears = NULL;
	ns_config->starting_linears = (vec_mp **)br_malloc( (W.num_variables - target_dim)*sizeof(vec_mp *));
	for (ii=0; ii<target_dim; ii++) {
		ns_config->starting_linears[ii] = (vec_mp *) br_malloc((randomized_degrees[ii]-1)*sizeof(vec_mp)); //subtract 1 for differentiation
		for (jj=0; jj<randomized_degrees[ii]-1; jj++) {
			init_vec_mp(ns_config->starting_linears[ii][jj],W.num_variables);
			ns_config->starting_linears[ii][jj]->size = W.num_variables;
			//			set_zero_mp(&ns_config->starting_linears[ii][jj]->coord[0]);
			for (kk=0; kk<W.num_variables; kk++) {
				get_comp_rand_mp(&ns_config->starting_linears[ii][jj]->coord[kk]);
			}
		}
	}
	
	
	init_mat_mp(ns_config->randomizer_matrix,randomizer_matrix->rows, randomizer_matrix->cols);
	ns_config->randomizer_matrix->rows = randomizer_matrix->rows;
	ns_config->randomizer_matrix->cols = randomizer_matrix->cols;
	mat_cp_mp(ns_config->randomizer_matrix, randomizer_matrix);

	
	return;
}














