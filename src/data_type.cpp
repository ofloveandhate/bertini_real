#include "data_type.hpp"
#include "programConfiguration.hpp"



std::string enum_lookup(int flag)
{
	switch (flag) {
		case SUCCESSFUL:
			return "SUCCESSFUL";
			break;
			
		case CRITICAL_FAILURE:
			return "CRITICAL_FAILURE";
			break;
			
		case TOLERABLE_FAILURE:
			return "TOLERABLE_FAILURE";
			break;
			
		case UNSET:
			return "UNSET";
			break;
			
		case CRITICAL:
			return "CRITICAL";
			break;
			
		case NEW:
			return "NEW";
			break;
			
		case MIDPOINT:
			return "MIDPOINT";
			break;
			
		case ISOLATED:
			return "ISOLATED";
			break;			
			
		case NULLSPACE:
			return "NULLSPACE";
			break;
			
		case LINPRODTODETJAC:
			return "LINPRODTODETJAC";
			break;
			
		case DETJACTODETJAC:
			return "DETJACTODETJAC";
			break;
			
		case LINTOLIN:
			return "LINTOLIN";
			break;
			
		case MULTILIN:
			return "MULTILIN";
			break;
		
		case MIDPOINT_SOLVER:
			return "MIDPOINT_SOLVER";
			break;
			
		case SPHERE_SOLVER:
			return "SPHERE_SOLVER";
			break;
						
			
		case TERMINATE:
			return "TERMINATE";
			break;
			
		case INITIAL_STATE:
			return "INITIAL_STATE";
			break;
			
			
		case PARSING:
			return "PARSING";
			break;
			
		case TYPE_CONFIRMATION:
			return "TYPE_CONFIRMATION";
			break;
			
		case DATA_TRANSMISSION:
			return "DATA_TRANSMISSION";
			break;
			
		case NUMPACKETS:
			return "NUMPACKETS";
			break;
			
		case INACTIVE:
			return "INACTIVE";
			break;
			
		case VEC_MP:
			return "VEC_MP";
			break;
			
		case VEC_D:
			return "VEC_D";
			break;
			
		case MAT_MP:
			return "MAT_MP";
			break;
			
		case MAT_D:
			return "MAT_D";
			break;
			
		case COMP_MP:
			return "COMP_MP";
			break;
			
		case COMP_D:
			return "COMP_D";
			break;
			
		case INDICES:
			return "INDICES";
			break;
			
			
		default:
			break;
	}
	
	return "unknown...  check out data_type.cpp";
}




















void system_randomizer::randomize(vec_d randomized_func_vals, mat_d randomized_jacobian,
								  vec_d func_vals, mat_d jacobian_vals,
								  comp_d hom_var)
{
	bool bail_out = false;
	
	if (!is_ready()) {
		std::cout << "trying to randomize_d, but is not set up!" << std::endl;
		bail_out = true;
	}
	
	if (func_vals->size != this->num_original_funcs) {
		std::cout << "mismatch in number of expected input functions (" << num_original_funcs << ") and actual inputted-number of functions (" << func_vals->size << ")." << std::endl;
		bail_out = true;
	}
	
	if (jacobian_vals->rows != this->num_original_funcs) {
		std::cout << "mismatch in size of expected input jacobian (" << num_original_funcs << ") and actual number of rows in jacobian (" << jacobian_vals->rows << ")." << std::endl;
		bail_out = true;
	}
	
	
	if (bail_out)
		br_exit(-9509);
	
	
	if (is_square()) {
		mat_cp_d(randomized_jacobian,jacobian_vals);
		vec_cp_d(randomized_func_vals,func_vals);
		return;
	}
	
	
	
	//ensure outputs are of correct size.
	increase_size_mat_d(randomized_jacobian, num_randomized_funcs, jacobian_vals->cols);
	randomized_jacobian->rows = num_randomized_funcs; randomized_jacobian->cols = jacobian_vals->cols;
	
	increase_size_vec_d(randomized_func_vals, num_randomized_funcs);
	randomized_func_vals->size = num_randomized_funcs;
	
	// do a little precomputation
	//0th entry is 1, having been set previously.
	for (int ii=1; ii<=max_degree_deficiency; ii++) {
		mul_d(&temp_homogenizer_d->coord[ii],&temp_homogenizer_d->coord[ii-1],hom_var);
	}
	
	increase_size_mat_d(temp_jac_d,num_original_funcs,jacobian_vals->cols);
	temp_jac_d->rows = num_original_funcs; temp_jac_d->cols = jacobian_vals->cols;
	
	// we do both function and jacobian randomization in the same loop for optimization.
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		// this must be done in a loop because each output (randomized function) has a different homogeneous structure.
		//
		//TODO: optimization could be done by looking at the previous function's structure and omitting previously done calculations (for the second and subsequent functions)
		
		
		// copy the current randomizer coefficients into a single row matrix
		for (int jj=0; jj<num_original_funcs; jj++) {
			set_d(&single_row_input_d->entry[0][jj], &randomizer_matrix_d->entry[ii][jj]); //TODO this is repetitive, and wasteful.  optimize away.
		} // at least is is just setting, not multiplying
		
		
		
		
		/////////////////
		//
		//  functions
		//
		//////////////
		increase_size_vec_d(temp_funcs_d, num_original_funcs);  temp_funcs_d->size = num_original_funcs;
		
		for (int jj=0; jj<num_original_funcs; jj++) {
			// structure_matrix[ii][jj] gives the degree deficiency
			if (structure_matrix[ii][jj]>0) { // if must homogenize at least one degree.
				mul_d(&temp_funcs_d->coord[jj],&func_vals->coord[jj], &temp_homogenizer_d->coord[ structure_matrix[ii][jj] ]);
			}
			else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
			{
				set_d(&temp_funcs_d->coord[jj], &func_vals->coord[jj]);
			}
			else // yeah...  maybe optimize away?
			{
				set_d(&temp_funcs_d->coord[jj], &func_vals->coord[jj]);
			}
			
		}
		
		mul_mat_vec_d(temp_vec_d, single_row_input_d, temp_funcs_d);
		set_d(&randomized_func_vals->coord[ii], &temp_vec_d->coord[0]);
		
		
		
		
		/////////////////
		//
		//  jacobian
		//
		//////////////
		
		
		for (int kk=0; kk<jacobian_vals->cols; kk++) { // kk indexes the variables in jacobian_vals (columns in jacobian)
			for (int jj=0; jj<num_original_funcs; jj++) { // jj indexes the original functions (rows in jacobian)
													  // structure_matrix[ii][jj] gives the degree deficiency
				if (structure_matrix[ii][jj]>0) { // must homogenize at least one degree.
					mul_d(&temp_jac_d->entry[jj][kk],&jacobian_vals->entry[jj][kk], &temp_homogenizer_d->coord[ structure_matrix[ii][jj] ]);
				}
				else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
				{
					set_d(&temp_jac_d->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				else
				{
					set_d(&temp_jac_d->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				
				
			}
		}
		
		
		// actually randomize here
		mat_mul_d(temp_mat_d, single_row_input_d, temp_jac_d);
		
		
		//copy the output into the returned value.
		for (int kk=0; kk<jacobian_vals->cols; kk++) {
			set_d(&randomized_jacobian->entry[ii][kk], &temp_mat_d->entry[0][kk]);
		}
		
		
		// for the last step, the first variable is the hom_var, and it must use the product rule.
		for (int jj=0; jj<num_original_funcs; jj++) {
			if (structure_matrix[ii][jj]>0) { // only if we need to actually homogenize this function
											  // we abuse the temp_funcs here and use as a temp storage area.
				mul_d(&temp_funcs_d->coord[jj], //d•h^(d-1)
					  &integer_coeffs_d->coord[ structure_matrix[ii][jj] ], &temp_homogenizer_d->coord[ structure_matrix[ii][jj]-1 ]); // this could be optimized if degree deficiency == 1.
				mul_d(&temp_funcs_d->coord[jj], &temp_funcs_d->coord[jj], &func_vals->coord[jj]);
				//temp_funcs_ = d•h^(d-1)•f_{jj}
			}
			else
			{
				set_zero_d(&temp_funcs_d->coord[jj]);
			}
		}
		
		//randomize
		mul_mat_vec_d(temp_vec_d, single_row_input_d, temp_funcs_d); // this produces a single number as output.
																	 // M = R•(hommed_f)
																	 // combine for power rule
		add_d(&randomized_jacobian->entry[ii][0], &randomized_jacobian->entry[ii][0], &temp_vec_d->coord[0]);
		
	}
	
	
	
}


void system_randomizer::randomize(vec_mp randomized_func_vals, mat_mp randomized_jacobian,
								  vec_mp func_vals, mat_mp jacobian_vals,
								  comp_mp hom_var)
{
	bool bail_out = false;
	
	if (!is_ready()) {
		std::cout << "trying to randomize_d, but is not set up!" << std::endl;
		bail_out = true;
	}
	
	if (func_vals->size != this->num_original_funcs) {
		std::cout << "mismatch in number of expected input functions (" << num_original_funcs << ") and actual inputted-number of functions (" << func_vals->size << ")." << std::endl;
		bail_out = true;
	}
	
	if (jacobian_vals->rows != this->num_original_funcs) {
		std::cout << "mismatch in size of expected input jacobian (" << num_original_funcs << ") and actual number of rows in jacobian (" << jacobian_vals->rows << ")." << std::endl;
		bail_out = true;
	}
	
	
	if (bail_out)
		br_exit(-9510);
	
	
	
	if (is_square()) {
		mat_cp_mp(randomized_jacobian,jacobian_vals);
		vec_cp_mp(randomized_func_vals,func_vals);
		return;
	}
	
	//ensure outputs are of correct size.
	increase_size_mat_mp(randomized_jacobian, num_randomized_funcs, jacobian_vals->cols);
	randomized_jacobian->rows = num_randomized_funcs; randomized_jacobian->cols = jacobian_vals->cols;
	
	increase_size_vec_mp(randomized_func_vals, num_randomized_funcs);
	randomized_func_vals->size = num_randomized_funcs;
	
	// do a little precomputation
	//0th entry is 1, having been set previously.
	for (int ii=1; ii<=max_degree_deficiency; ii++) {
		mul_mp(&temp_homogenizer_mp->coord[ii],&temp_homogenizer_mp->coord[ii-1],hom_var);
	}
	
	increase_size_mat_mp(temp_jac_mp,num_original_funcs,jacobian_vals->cols);
	temp_jac_mp->rows = num_original_funcs;
	temp_jac_mp->cols = jacobian_vals->cols;
	
	// we do both function and jacobian randomization in the same loop for optimization.
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		// this must be done in a loop because each output (randomized function) has a different homogeneous structure.
		//
		// optimization could be done by looking at the previous function's structure and omitting previously done calculations (for the second and subsequent functions)
		
		
		// copy the current randomizer coefficients into a single row matrix
		for (int jj=0; jj<num_original_funcs; jj++) {
			set_mp(&single_row_input_mp->entry[0][jj], &randomizer_matrix_mp->entry[ii][jj]); // this is repetitive, and wasteful.  optimize away.
		} // at least is is just setting, not multiplying
		
		
		
		
		/////////////////
		//
		//  functions
		//
		//////////////
		
		increase_size_vec_mp(temp_funcs_mp, num_original_funcs);  temp_funcs_mp->size = num_original_funcs;
		
		for (int jj=0; jj<num_original_funcs; jj++) {
			// structure_matrix[ii][jj] gives the degree deficiency
			if (structure_matrix[ii][jj]>0) { // if must homogenize at least one degree.
				mul_mp(&temp_funcs_mp->coord[jj],&func_vals->coord[jj], &temp_homogenizer_mp->coord[ structure_matrix[ii][jj] ]);
			}
			else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
			{
				set_mp(&temp_funcs_mp->coord[jj], &func_vals->coord[jj]);
			}
			else // yeah...  maybe optimize away?
			{
				set_mp(&temp_funcs_mp->coord[jj], &func_vals->coord[jj]);
			}
			
		}
		
		mul_mat_vec_mp(temp_vec_mp, single_row_input_mp, temp_funcs_mp);
		set_mp(&randomized_func_vals->coord[ii], &temp_vec_mp->coord[0]);
		
		
		
		
		/////////////////
		//
		//  jacobian
		//
		//////////////
		
		
		for (int kk=0; kk<jacobian_vals->cols; kk++) { // kk indexes the variables in jacobian_vals (columns in jacobian)
			for (int jj=0; jj<num_original_funcs; jj++) { // jj indexes the original functions (rows in jacobian)
													  // structure_matrix[ii][jj] gives the degree deficiency
				if (structure_matrix[ii][jj]>0) { // must homogenize at least one degree.
					mul_mp(&temp_jac_mp->entry[jj][kk],&jacobian_vals->entry[jj][kk], &temp_homogenizer_mp->coord[ structure_matrix[ii][jj] ]);
				}
				else if (structure_matrix[ii][jj]==0) // no need to additionally homogenize -- already at the correct degree
				{
					set_mp(&temp_jac_mp->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				else
				{
					set_mp(&temp_jac_mp->entry[jj][kk],&jacobian_vals->entry[jj][kk]);
				}
				
				
			}
		}
		
		
		// actually randomize here
		mat_mul_mp(temp_mat_mp, single_row_input_mp, temp_jac_mp);
		
		
		//copy the output into the returned value.
		for (int kk=0; kk<jacobian_vals->cols; kk++) {
			set_mp(&randomized_jacobian->entry[ii][kk], &temp_mat_mp->entry[0][kk]);
		}
		
		
		// for the last step, the first variable is the hom_var, and it must use the product rule.
		for (int jj=0; jj<num_original_funcs; jj++) {
			if (structure_matrix[ii][jj]>0) { // only if we need to actually homogenize this function
											  // we abuse the temp_funcs here and use as a temp storage area.
				mul_mp(&temp_funcs_mp->coord[jj], //d•h^(d-1)
					   &integer_coeffs_mp->coord[ structure_matrix[ii][jj] ], &temp_homogenizer_mp->coord[ structure_matrix[ii][jj]-1 ]); // this could be optimized if degree deficiency == 1.
				mul_mp(&temp_funcs_mp->coord[jj], &temp_funcs_mp->coord[jj], &func_vals->coord[jj]);
				//temp_funcs_ = d•h^(d-1)•f_{jj}
			}
			else
			{
				set_zero_mp(&temp_funcs_mp->coord[jj]);
			}
		}
		
		//randomize
		mul_mat_vec_mp(temp_vec_mp, single_row_input_mp, temp_funcs_mp);
		// M = R•(hommed_f)
		// combine for power rule
		add_mp(&randomized_jacobian->entry[ii][0], &randomized_jacobian->entry[ii][0], &temp_vec_mp->coord[0]);
		
	}
}


void system_randomizer::change_prec(int new_prec)
{
	if (!is_ready()) {
		std::cout << "trying to change precision when not set up!" << std::endl;
		br_exit(65621);
	}
	change_prec_mat_mp(randomizer_matrix_mp, new_prec);
	mat_cp_mp(randomizer_matrix_mp,randomizer_matrix_full_prec);
	
	change_prec_vec_mp(temp_homogenizer_mp,new_prec);
	change_prec_vec_mp(temp_funcs_mp,new_prec);
	change_prec_mat_mp(temp_jac_mp,new_prec);
	
	
	change_prec_mat_mp(single_row_input_mp,new_prec);
	
	change_prec_vec_mp(integer_coeffs_mp,new_prec);
	change_prec_mat_mp(temp_mat_mp,new_prec);
	change_prec_vec_mp(temp_vec_mp,new_prec);
	
	
}

void system_randomizer::setup(int num_desired_rows, int num_funcs)
{
	if (num_desired_rows<=0) {
		std::cout << "requested a randomizer for <= 0 output functions..." << std::endl;
		br_exit(80146);
	}
	if (num_funcs<=0) {
		std::cout << "requested a randomizer for <= 0 input functions..." << std::endl;
		br_exit(80147);
	}
	
	if (num_desired_rows>num_funcs) {
		std::cout << "requested randomizer which would UP-randomize." << std::endl;
		std::cout << num_desired_rows << " > " << num_funcs << ", which is not allowed" << std::endl;
		br_exit(80148);
	}
	setup_indicator = false; // reset;
	randomized_degrees.resize(0);
	original_degrees.resize(0);
	max_base_degree = 0;
	
	//get unique degrees
	int *unique_degrees = (int *) br_malloc(num_funcs*sizeof(int));
	
	
	FILE *IN = safe_fopen_read("deg.out"); //open the deg.out file for reading.
	int num_unique_degrees = 0;
	int occurrence_counter;
	int tempdegree;
	//TODO: this doesn't read correctly if there is more than one variable group.
	for (int ii=0; ii<num_funcs; ++ii) {
		fscanf(IN,"%d\n",&tempdegree); // read data
		original_degrees.push_back(tempdegree);
		
		occurrence_counter = 0; // set the counter for how many times the current degree has already been found.
		for (int jj=0; jj<ii; jj++) {
			if (original_degrees[jj]==tempdegree) { // if previously stored degree is same as current one
				occurrence_counter++; // increment counter
			}
		}
		
		if (occurrence_counter==0) { // if did not find already in list
			unique_degrees[num_unique_degrees] = tempdegree; // add to list of unique degrees.
			num_unique_degrees++; // have one more unique degree
		} // re: jj
		
		if (tempdegree>max_base_degree) {
			max_base_degree = tempdegree;
		}
		
	}// re: ii
	fclose(IN);
	
	
	
	num_randomized_funcs = num_desired_rows;
	num_original_funcs = num_funcs;
	
	

	square_indicator = false;
	
	//sort the unique degrees into decreasing order
	qsort(unique_degrees, num_unique_degrees, sizeof(int), compare_integers_decreasing);
	
	//count how many of each unique degree there are.
	int *num_of_each_degree = (int *) br_malloc(num_unique_degrees*sizeof(int));
	for (int ii=0; ii<num_unique_degrees; ii++) {
		num_of_each_degree[ii] = 0;
		for (int jj=0; jj<num_funcs; ++jj) {
			if (unique_degrees[ii]==original_degrees[jj]) {
				num_of_each_degree[ii]++;
			}
		}
	}
	
	
	
	
	//resize the matrix
	change_size_mat_mp(randomizer_matrix_full_prec,num_desired_rows,num_funcs);
	randomizer_matrix_full_prec->rows = num_desired_rows;
	randomizer_matrix_full_prec->cols = num_funcs;
	
	structure_matrix.resize(num_desired_rows);
	for (auto iter = structure_matrix.begin(); iter!=structure_matrix.end(); ++iter) {
		iter->resize(num_funcs);
	}
	
	
	
	
	if (num_desired_rows==num_funcs) {
		for (int ii=0; ii<num_desired_rows; ii++) {
			randomized_degrees.push_back(original_degrees[ii]);
			
			for (int jj=0; jj<num_funcs; jj++) {
				if (ii==jj) {
					structure_matrix[ii][jj] = max_base_degree - original_degrees[ii];
				}
				else{
					structure_matrix[ii][jj] = 0;
				}
			}
		}
		make_matrix_ID_mp(randomizer_matrix_full_prec,num_funcs,num_funcs);
		square_indicator = true;
	}
	else{
		
		int counter = 0;
		int current_degree_index = 0; // start at the end
		
		for (int ii=0; ii<num_desired_rows; ii++) {
		
			counter++;
			if (counter>num_of_each_degree[current_degree_index]) {
				current_degree_index++;
				counter = 1;
			}
			
			int current_degree = unique_degrees[current_degree_index];
			randomized_degrees.push_back(current_degree);
			
			int encountered_current_degree = 0;
			for (int jj=0; jj<num_funcs; jj++) {
				if ( original_degrees[jj]<= current_degree ) {
					encountered_current_degree++;
					if (encountered_current_degree >= counter){
						get_comp_rand_real_mp(&randomizer_matrix_full_prec->entry[ii][jj]);
						structure_matrix[ii][jj] = current_degree - original_degrees[jj];  // deficiency level
					}
					else{
						set_zero_mp(&randomizer_matrix_full_prec->entry[ii][jj]);
						structure_matrix[ii][jj] = 0;
					}
				}
				else
				{
					set_zero_mp(&randomizer_matrix_full_prec->entry[ii][jj]);
					structure_matrix[ii][jj] = 0;
				}
			}
			

		}
		
		square_indicator = false;
	}
	
	
	max_degree_deficiency = unique_degrees[0] - unique_degrees[num_unique_degrees-1];
	
	
	
	mat_cp_mp(randomizer_matrix_mp,randomizer_matrix_full_prec);
	mat_mp_to_d(randomizer_matrix_d,randomizer_matrix_full_prec);
	
	
	free(num_of_each_degree);
	free(unique_degrees);
	
	
	
	
	
	
	setup_temps();
	
	
}


void system_randomizer::setup_temps()
{
	change_size_vec_d(integer_coeffs_d,max_degree_deficiency+1);  integer_coeffs_d->size = max_degree_deficiency+1;
	change_size_vec_mp(integer_coeffs_mp,max_degree_deficiency+1); integer_coeffs_mp->size = max_degree_deficiency+1;
	
	for (int ii=0; ii<=max_degree_deficiency; ii++) {
		integer_coeffs_d->coord[ii].r = ii; integer_coeffs_d->coord[ii].i = 0;
		set_zero_mp(&integer_coeffs_mp->coord[ii]); // initialize
		mpf_set_d(integer_coeffs_mp->coord[ii].r,ii); // set to integer value
	}
	
	
	change_size_mat_d(single_row_input_d,1,num_original_funcs);  single_row_input_d->rows = 1; single_row_input_d->cols = num_original_funcs;
	change_size_mat_mp(single_row_input_mp,1,num_original_funcs);  single_row_input_mp->rows = 1; single_row_input_mp->cols = num_original_funcs;
	
	change_size_vec_d(temp_homogenizer_d,max_degree_deficiency+1);  temp_homogenizer_d->size = max_degree_deficiency+1;
	change_size_vec_mp(temp_homogenizer_mp,max_degree_deficiency+1); temp_homogenizer_mp->size = max_degree_deficiency+1;
	set_one_d(&temp_homogenizer_d->coord[0]);
	set_one_mp(&temp_homogenizer_mp->coord[0]);
	
	
	
	
	setup_indicator = true;
}




void system_randomizer::send(int target, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "system_randomizer::send" << std::endl;
#endif
	int sendme = setup_indicator;
	MPI_Send(&sendme,1,MPI_INT,target,UNUSED,MPI_COMM_WORLD);

	if (!setup_indicator) {
		std::cout << "bailing on sending upsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	buffer[0] = square_indicator;
	buffer[1] = num_randomized_funcs;
	buffer[2] = num_original_funcs;
	buffer[3] = max_base_degree;
	buffer[4] = max_degree_deficiency;
	MPI_Send(buffer,5,MPI_INT,target, SYSTEM_RANDOMIZER ,mpi_config.comm());
	
	delete[] buffer;
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_send = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_send];
	
	
	
	if (randomized_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = randomized_degrees.begin(); iter!=randomized_degrees.end(); iter++) {
			buffer2[cnt] = *iter;
			cnt++;
		}
		
	}
	
	int offset = num_randomized_funcs;
	if (original_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = original_degrees.begin(); iter!=original_degrees.end(); iter++) {
			buffer2[cnt+offset] = *iter;
			cnt++;
		}
	}
	
	
	offset += num_original_funcs;
	
	int cnt = 0;
	for (auto iter = structure_matrix.begin(); iter!=structure_matrix.end(); iter++) {
		for (auto jter = iter->begin(); jter != iter->end(); jter++) {
			buffer2[cnt+offset] = *jter;
			cnt++;
		}
		
	}

	
	
	
	MPI_Send(buffer2, size_to_send, MPI_INT, target, SYSTEM_RANDOMIZER, mpi_config.comm());
	delete [] buffer2;
	
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		send_mat_mp(randomizer_matrix_full_prec, target);
	}
	
	
	
}


void system_randomizer::receive(int source, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "system_randomizer::receive" << std::endl;
#endif
	
	MPI_Status statty_mc_gatty;
	
	int recvme;
	MPI_Recv(&recvme,1,MPI_INT,source,UNUSED,MPI_COMM_WORLD,&statty_mc_gatty);
	setup_indicator = recvme;
	if (!setup_indicator) {
		std::cout << "bailing from getting unsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	MPI_Recv(buffer,5,MPI_INT,source,SYSTEM_RANDOMIZER,mpi_config.comm(),&statty_mc_gatty);
	
	square_indicator = buffer[0];
	num_randomized_funcs = buffer[1];
	num_original_funcs = buffer[2];
	max_base_degree = buffer[3];
	max_degree_deficiency = buffer[4];
	
	
	delete[] buffer;
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_receive = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_receive];
	
	
	
	
	MPI_Recv(buffer2, size_to_receive, MPI_INT, source, SYSTEM_RANDOMIZER, mpi_config.comm(),&statty_mc_gatty);
	
	int cnt = 0;
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		randomized_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	for (int ii=0; ii<num_original_funcs; ii++) {
		original_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	structure_matrix.resize(num_randomized_funcs);
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		for (int jj=0; jj<num_original_funcs; jj++) {
			structure_matrix[ii].push_back(buffer2[cnt]);
			cnt++;
		}
	}
	
	
	delete [] buffer2;
	
	change_size_mat_mp(randomizer_matrix_full_prec,num_randomized_funcs,num_original_funcs);
	randomizer_matrix_full_prec->rows = num_randomized_funcs;
	randomizer_matrix_full_prec->cols = num_original_funcs;
	
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		receive_mat_mp(randomizer_matrix_full_prec, source);
	}
	
	
	setup_temps();
	// temps and things set up after you get the matrix and a few parameters.

	
}



void system_randomizer::bcast_send(parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "system_randomizer::bcast_send" << std::endl;
#endif
	
	
	int sendme = setup_indicator;
	MPI_Bcast(&sendme,1,MPI_INT,mpi_config.head(),MPI_COMM_WORLD);
	
	if (!setup_indicator) {
		std::cout << "bailing on sending unsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	buffer[0] = square_indicator;
	buffer[1] = num_randomized_funcs;
	buffer[2] = num_original_funcs;
	buffer[3] = max_base_degree;
	buffer[4] = max_degree_deficiency;
	MPI_Bcast(buffer,5,MPI_INT,mpi_config.head(),MPI_COMM_WORLD);
	
	delete[] buffer;
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_send = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_send];
	
	
	
	if (randomized_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = randomized_degrees.begin(); iter!=randomized_degrees.end(); iter++) {
			buffer2[cnt] = *iter;
			cnt++;
		}
		
	}
	
	int offset = num_randomized_funcs;
	if (original_degrees.size()>0) {
		int cnt = 0;
		for (auto iter = original_degrees.begin(); iter!=original_degrees.end(); iter++) {
			buffer2[cnt+offset] = *iter;
			cnt++;
		}
	}
	
	
	offset += num_original_funcs;
	
	int cnt = 0;
	for (auto iter = structure_matrix.begin(); iter!=structure_matrix.end(); iter++) {
		for (auto jter = iter->begin(); jter != iter->end(); jter++) {
			buffer2[cnt+offset] = *jter;
			cnt++;
		}
		
	}
	
	

	MPI_Bcast(buffer2, size_to_send, MPI_INT, mpi_config.head(), MPI_COMM_WORLD);
	delete [] buffer2;
	
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		bcast_mat_mp(randomizer_matrix_full_prec, mpi_config.id(), mpi_config.head());
	}
	
	

	
	
	
	
}


void system_randomizer::bcast_receive(parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "system_randomizer::bcast_receive" << std::endl;
#endif
	
	
	
	
	
	int recvme;
	MPI_Bcast(&recvme,1,MPI_INT,mpi_config.head(),MPI_COMM_WORLD);
	setup_indicator = recvme;
	if (!setup_indicator) {
		std::cout << "bailing from getting unsetup randomizer" << std::endl;
		return;
	}
	
	int *buffer = new int[1+2+2]; // square + size_of_randomizer + max_base + max_degree
	MPI_Bcast(buffer,5,MPI_INT,mpi_config.head(),MPI_COMM_WORLD);
	
	square_indicator = buffer[0];
	num_randomized_funcs = buffer[1];
	num_original_funcs = buffer[2];
	max_base_degree = buffer[3];
	max_degree_deficiency = buffer[4];
	
	
	delete[] buffer;

	
	
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_receive = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_receive];
	MPI_Bcast(buffer2, size_to_receive, MPI_INT, mpi_config.head(), MPI_COMM_WORLD);
	
	int cnt = 0;
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		randomized_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	for (int ii=0; ii<num_original_funcs; ii++) {
		original_degrees.push_back(buffer2[cnt]);
		cnt++;
	}
	
	structure_matrix.resize(num_randomized_funcs);
	for (int ii=0; ii<num_randomized_funcs; ii++) {
		for (int jj=0; jj<num_original_funcs; jj++) {
			structure_matrix[ii].push_back(buffer2[cnt]);
			cnt++;
		}
	}
	
	
	delete [] buffer2;
	
	change_size_mat_mp(randomizer_matrix_full_prec,num_randomized_funcs,num_original_funcs);
	randomizer_matrix_full_prec->rows = num_randomized_funcs;
	randomizer_matrix_full_prec->cols = num_original_funcs;
	
	if ( (randomizer_matrix_full_prec->rows != 0) || (randomizer_matrix_full_prec->cols != 0)) {
		bcast_mat_mp(randomizer_matrix_full_prec, mpi_config.id(), mpi_config.head());
		mat_cp_mp(randomizer_matrix_mp,randomizer_matrix_full_prec);
		mat_mp_to_d(randomizer_matrix_d,randomizer_matrix_full_prec);
	}
	
	
	setup_temps();
	// temps and things set up after you get the matrix and a few parameters.

}





int point_holder::add_point(vec_mp new_point)
{
	
	if (num_pts_!=0 && this->pts_mp_==NULL) {
		printf("trying to add point to point_holder with non-zero num_points and NULL container!\n");
		br_exit(9713);
	}
	
	if (num_pts_==0 && this->pts_mp_!=NULL) {
		printf("trying to add point to point_holder with num_points==0 and non-NULL container!\n");
		br_exit(9713);
	}
	
	
	if (num_pts_==0) {
		pts_mp_ = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		pts_mp_ = (vec_mp *)br_realloc(pts_mp_, (num_pts_+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(pts_mp_[num_pts_], new_point->size, new_point->curr_prec);
	pts_mp_[num_pts_]->size = new_point->size;
	vec_cp_mp(pts_mp_[num_pts_], new_point);
	
	num_pts_++;
	
	return num_pts_-1;
}

int patch_holder::add_patch(vec_mp new_patch)
{
	
	if (num_patches_!=0 && patch_mp_==NULL) {
		printf("trying to add patch to witness set with non-zero num_patches and NULL container!\n");
		deliberate_segfault();
	}
	
	if (num_patches_==0 && patch_mp_!=NULL) {
		printf("trying to add point to witness set with num_points==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (num_patches_==0) {
		patch_mp_ = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		patch_mp_ = (vec_mp *)br_realloc(patch_mp_, (num_patches_+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(patch_mp_[num_patches_], new_patch->size,new_patch->curr_prec);
	patch_mp_[num_patches_]->size = new_patch->size;
	vec_cp_mp(patch_mp_[num_patches_], new_patch);
	
	this->num_patches_++;
	
	
	int burnsome = 0;
	for (unsigned int ii=0; ii<this->num_patches_; ii++) {
		burnsome += this->patch_mp_[ii]->size;
	}
	

	return num_patches_-1;
}


int linear_holder::add_linear(vec_mp new_linear)
{
	
	if (num_linears_!=0 && L_mp_==NULL) {
		printf("trying to add linear to linear holder with non-zero num_linears and NULL container!\n");
		br_exit(9711);
	}
	
	if (num_linears_==0 && L_mp_!=NULL) {
		printf("trying to add linear to linear holder with num_linears==0 and non-NULL container!\n");
		br_exit(9711);
	}
	
	
	if (num_linears_==0) {
		L_mp_ = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		L_mp_ = (vec_mp *)br_realloc(L_mp_, (num_linears_+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(L_mp_[num_linears_], new_linear->size,new_linear->curr_prec);
	L_mp_[num_linears_]->size = new_linear->size;
	vec_cp_mp(L_mp_[num_linears_], new_linear);
	
	this->num_linears_++;
	
	
	return num_linears_-1;
}







int witness_set::witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars)
{
	
	
	
	
	FILE *IN = safe_fopen_read(witness_set_file);
	
	
	int temp_num_patches, patch_size, temp_num_linears, temp_num_points, num_vars_in_linears;
	
	
	fscanf(IN, "%d %d %d", &temp_num_points, &dim_, &comp_num_); scanRestOfLine(IN);
	
	
	
	this->num_vars_ = num_vars;
	this->num_natty_vars_ = num_vars;
	
	
	vec_mp temp_vec;  init_vec_mp2(temp_vec, num_vars,1024); temp_vec->size = num_vars;
	
	for (int ii=0; ii < temp_num_points; ii++) {
		
		//read the witness points into memory
		for (int jj=0; jj < num_vars; ++jj) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10); // 10 is the base
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
			
			scanRestOfLine(IN);
		}
		
		add_point(temp_vec);
	}
	
	
	
	
	
	fscanf(IN, "%d %d", &temp_num_linears, &num_vars_in_linears);  scanRestOfLine(IN);
	
	
	for (int ii=0; ii < temp_num_linears; ii++) {
		change_size_vec_mp(temp_vec,num_vars_in_linears);
		temp_vec->size = num_vars_in_linears;
		
		//read the witness linears into memory
		for (int jj=0; jj < num_vars_in_linears; jj++) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
		}
		
		add_linear(temp_vec);
	}
	
	
	
	
	fscanf(IN, "%d %d", &temp_num_patches, &patch_size); scanRestOfLine(IN);
	
	if (temp_num_patches>1) {
		std::cerr << temp_num_patches << " patches detected.  this probably indicates a problem." << std::endl;
		std::cerr << "the file being read: " << witness_set_file << std::endl;
		std::cerr << "trying to read " << num_vars << " variables." << std::endl;
		mypause();
	}
	
	
	for (int ii=0; ii < temp_num_patches; ii++) {
		
		change_size_vec_mp(temp_vec,patch_size);
		temp_vec->size = patch_size;
		
		//read the patch into memory
		for (int jj=0; jj < patch_size; jj++) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
		}
		
		add_patch(temp_vec);
	}
	
	fclose(IN);
	
	
	clear_vec_mp(temp_vec);
	
	return 0;
}










void witness_set::write_homogeneous_coordinates(boost::filesystem::path filename) const
{

	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%zu\n\n",num_points()); // print the header line
	
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		vec_mp &curr_point = point(ii);
		for ( int jj=0; jj< curr_point->size; jj++) {
			print_mp(OUT,0,&( curr_point->coord[jj]));
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void witness_set::write_dehomogenized_coordinates(boost::filesystem::path filename) const
{

	
	vec_mp result; init_vec_mp(result,1);
	
	
	FILE *OUT = safe_fopen_write(filename.c_str()); // open the output file.
	
	fprintf(OUT,"%zu\n\n",num_points()); // print the header line
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		if (this->num_synth_vars()>0) {
			dehomogenize(&result,point(ii), num_natty_vars_);
		}
		else{
			dehomogenize(&result,point(ii));
		}
		
		for (int jj=0; jj<num_natty_vars_-1; jj++) {
			print_mp(OUT, 0, &result->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
	clear_vec_mp(result);
	return;
}



void witness_set::write_dehomogenized_coordinates(boost::filesystem::path filename,std::set<unsigned int> indices) const
{
	
	for (auto ii=indices.begin(); ii!=indices.end(); ++ii) {
		if (*ii >= this->num_points()) {
			std::cout << "requested to print out-of-range point index " << *ii << " to a dehomogenized file." << std::endl;
			std::cout << "[this witness_set contains " << num_points() << " points.]" << std::endl;
			br_exit(66190);
		}
	}
	
	
	
	
	vec_mp result; init_vec_mp(result,1);
	
	
	FILE *OUT = safe_fopen_write(filename.c_str()); // open the output file.
	
	fprintf(OUT,"%lu\n\n",indices.size()); // print the header line
	for (auto ii=indices.begin(); ii!=indices.end(); ++ii) {
		if (this->num_synth_vars()>0) {
			dehomogenize(&result,this->point(*ii), num_natty_vars_);
		}
		else{
			dehomogenize(&result,this->point(*ii));
		}
		
		for (int jj=0; jj<num_natty_vars_-1; jj++) {
			print_mp(OUT, 0, &result->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
	clear_vec_mp(result);
	return;
}



void witness_set::write_linears(boost::filesystem::path filename) const
{

	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%zu\n\n",num_linears()); // print the header line
	
	for (unsigned int ii=0; ii<num_linears(); ++ii) {
		for (int jj=0; jj<linear(ii)->size; jj++) {
			print_mp(OUT, 0, &linear(ii)->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}



void witness_set::print_patches(boost::filesystem::path filename) const
{

	
	FILE *OUT  = safe_fopen_write(filename); // open the output file
	
	
	fprintf(OUT,"%zu\n\n",num_patches()); // print the header line
	
	for (unsigned int ii=0; ii<num_patches(); ++ii) {
		fprintf(OUT,"%d\n",patch(ii)->size);
		for (int jj=0; jj<patch(ii)->size; jj++) {
			print_mp(OUT, 0, &patch(ii)->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void witness_set::read_patches_from_file(boost::filesystem::path filename)
{
	
	FILE *IN = safe_fopen_read(filename);
	
	witness_set::reset_patches();
	int curr_num_patches;
	fscanf(IN,"%d\n",&curr_num_patches);
	
	vec_mp temp_patch; init_vec_mp2(temp_patch,1,1024); temp_patch->size = 1;
	for (int ii=0; ii<curr_num_patches; ii++) {
		int curr_size;
		fscanf(IN,"%d\n",&curr_size);
		change_size_vec_mp(temp_patch,curr_size); temp_patch->size = curr_size;
		
		for (int jj=0; jj<curr_size; jj++) {
			mpf_inp_str(temp_patch->coord[jj].r, IN, 10);
			mpf_inp_str(temp_patch->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
		}
		
		witness_set::add_patch(temp_patch);
	}
	
	clear_vec_mp(temp_patch);
	fclose(IN);
	
	return;
}



void witness_set::print_to_screen() const
{
	
	vec_mp dehom;  init_vec_mp(dehom,1); dehom->size = 1;
	
	std::stringstream varname;
	
	std::cout << "witness set has " << num_vars_ << " total variables, " << num_natty_vars_ << " natural variables." << std::endl;
	
	
	std::cout << "dim " << dim_ << ", comp " << comp_num_ << std::endl;
	std::cout << "input file name " << input_filename_ << std::endl;

	
	printf("******\n%zu points\n******\n",num_points());
	std::cout << color::green();
	for (unsigned ii=0; ii<num_points(); ii++) {
		
		dehomogenize(&dehom, point(ii), num_natty_vars_);
		
		varname << "point_" << ii;
		
		print_point_to_screen_matlab(dehom,varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::blue();
	printf("******\n%zu linears\n******\n",num_linears());
	
	for (unsigned ii=0; ii<num_linears(); ii++) {
		varname << "linear_" << ii;
		print_point_to_screen_matlab(linear(ii),varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::cyan();
	printf("******\n%zu patches\n******\n",num_patches());
	
	for (unsigned ii=0; ii<num_patches(); ii++) {
		varname << "patch_" << ii;
		print_point_to_screen_matlab(patch(ii),varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << "variable names:\n";
	for (unsigned ii=0; ii< num_var_names(); ii++) {
		std::cout << name(ii) << "\n";
	}
	printf("\n\n");
	
	
	clear_vec_mp(dehom);
}


void witness_set::print_to_file(boost::filesystem::path filename) const
{
	// print back into the same format we parse from.
	
	
	FILE *OUT = safe_fopen_write(filename);
	
	
	fprintf(OUT, "%zu %d %d\n\n", num_points(), dim_, comp_num_);
	


	for (unsigned int ii=0; ii < num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		//print the witness points into file
		for (int jj=0; jj < curr_point->size; ++jj) {
			mpf_out_str(OUT,10,0,curr_point->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,curr_point->coord[jj].i);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	
	
	
	fprintf(OUT, "%zu %d\n", num_linears(), num_vars_);
	
	
	for (unsigned int ii=0; ii < num_linears(); ii++) {

		for (int jj=0; jj < linear(ii)->size; jj++) {
			mpf_out_str(OUT,10,0,linear(ii)->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,linear(ii)->coord[jj].i);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	
	
	fprintf(OUT, "%zu %d\n", num_patches(), num_vars_);// TODO: this is incorrect
	
	for (unsigned int ii=0; ii < num_patches(); ii++) {
		
		vec_mp & curr_patch = patch(ii);
		
		//read the patch into memory
		for (int jj=0; jj < curr_patch->size; jj++) {
			mpf_out_str(OUT,10,0,curr_patch->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,curr_patch->coord[jj].i);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	fclose(OUT);
	
	
	
	
	
	return;
}


void witness_set::only_natural_vars()
{
	witness_set::only_first_vars(num_natty_vars_);
}


void witness_set::only_first_vars(int num_vars)
{
	
	vec_mp tempvec;  init_vec_mp2(tempvec, num_vars, 1024);
	tempvec->size = num_vars;
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		
		for (int jj=0; jj<num_vars; jj++) {
			set_mp(&tempvec->coord[jj], &curr_point->coord[jj]);
		}
		
		change_size_vec_mp(curr_point, num_vars);  curr_point->size = num_vars;
		vec_cp_mp(curr_point, tempvec);
	}
	
	
	this->num_vars_ = num_vars;
	
	int patch_size_counter = 0, trim_from_here =  0;
	for (unsigned int ii=0; ii<num_patches(); ii++) {
		patch_size_counter += patch(ii)->size;
		if (patch_size_counter == num_vars)
		{
			trim_from_here = ii+1;
		}
	}
	
	if (trim_from_here==0) {
		std::cerr << "problem: the sum of the patch sizes never equalled the number of variables to trim to...\nhence, the trimming operation could not complete." << std::endl;
		this->print_to_screen();
		deliberate_segfault();
	}
	
	for (unsigned int ii=trim_from_here; ii<num_patches(); ii++) {
		clear_vec_mp(patch(ii));
	}
	
	patch_mp_ = (vec_mp *) br_realloc(patch_mp_, trim_from_here* sizeof(vec_mp));
	num_patches_ = trim_from_here;
	
	for (unsigned int ii=0; ii<num_linears(); ii++) {
		linear(ii)->size = num_vars;
	}
	
	clear_vec_mp(tempvec);
	return;
}



void witness_set::sort_for_real(tracker_config_t * T)
{
	

	
	
	int *real_indicator = new int[num_points()];
	int counter = 0;
	
	vec_mp result; init_vec_mp(result,num_natty_vars_-1);
	result->size = num_natty_vars_-1;
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		for (int jj=1; jj<num_natty_vars_; jj++) {
			div_mp(&result->coord[jj-1], &curr_point->coord[jj], &curr_point->coord[0]);
		}
		real_indicator[ii] = checkForReal_mp(result, T->real_threshold);
		
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	vec_mp *tempvec = (vec_mp *)br_malloc(counter * sizeof(vec_mp));
	
	counter = 0;  // reset
	for (unsigned int ii=0; ii<num_points(); ii++) {
		vec_mp & curr_point = point(ii);
		if (real_indicator[ii]==1) {
			
			init_vec_mp2(tempvec[counter],this->num_vars_,1024); tempvec[counter]->size = this->num_vars_;
			vec_cp_mp(tempvec[counter],curr_point);
			counter++;
		}
		else{
			
		}
	}
	
	clear_vec_mp(result);
	
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		clear_vec_mp(point(ii));
	}
	free(pts_mp_);
	
	pts_mp_ = tempvec;
	num_pts_ = counter;
	
	delete[] real_indicator;
	return;
}






// T is necessary for the tolerances.
void witness_set::sort_for_unique(tracker_config_t * T)
{
	
	if (num_vars_==0) {
		throw std::logic_error("sorting witness set with 0 variables for uniqueness");
	}
	
	int curr_uniqueness;
	int num_good_pts = 0;
	std::vector<int> is_unique;
	
	for (unsigned int ii = 0; ii<num_points(); ++ii) {
		vec_mp &curr_point = point(ii);
		
		curr_uniqueness = 1;
		
		int prev_size_1 = curr_point->size;  curr_point->size = num_natty_vars_; // cache and change to natural number
		
		for (unsigned int jj=ii+1; jj<num_points(); ++jj) {
			vec_mp & inner_point = point(jj);
			int prev_size_2 = inner_point->size; inner_point->size = num_natty_vars_; // cache and change to natural number
			if ( isSamePoint_homogeneous_input(curr_point,inner_point) ){
				curr_uniqueness = 0;
			}
			inner_point->size = prev_size_2; // restore
		}
		curr_point->size = prev_size_1; // restore
		
		if (curr_uniqueness==1) {
			is_unique.push_back(1);
			num_good_pts++;
		}
		else {
			is_unique.push_back(0);
		}
		
	}
	
	
	
	vec_mp *transferme = (vec_mp *)br_malloc(num_good_pts*sizeof(vec_mp));
	int counter = 0;
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		if (is_unique[ii]==1) {
			
			init_vec_mp2(transferme[counter],num_vars_,1024);  transferme[counter]->size = num_vars_;
			vec_cp_mp(transferme[counter], point(ii));
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		std::logic_error("counter mismatch");
	}
	
	if (num_points()>0) {
		for (unsigned int ii=0; ii<num_points(); ii++) {
			clear_vec_mp(point(ii));
		}
		free(pts_mp_);
	}
	
	
	num_pts_ = num_good_pts;
	pts_mp_ = transferme;
	
	return;
}





void witness_set::sort_for_inside_sphere(comp_mp radius, vec_mp center)
{
	
	
	
	
	int num_good_pts = 0;
	
	
	std::vector<int> is_ok;
	
	vec_mp temp_vec; init_vec_mp(temp_vec,0);
	comp_mp temp; init_mp(temp);
	
	for (unsigned int ii = 0; ii<num_points(); ++ii) {
		
		
		
		dehomogenize(&temp_vec, point(ii),num_natty_vars_);
		temp_vec->size = center->size;
		
		norm_of_difference(temp->r, temp_vec, center);
		if ( mpf_cmp(temp->r, radius->r) < 0   ){
			is_ok.push_back(1);
			num_good_pts++;
		}
		else
		{
			is_ok.push_back(0);
		}
		
		
		
	}
	
	
	
	vec_mp *transferme = (vec_mp *)br_malloc(num_good_pts*sizeof(vec_mp));
	int counter = 0;
	for (unsigned int ii=0; ii<num_points(); ++ii) {
		if (is_ok[ii]==1) {
			init_vec_mp2(transferme[counter],0,1024);  transferme[counter]->size = 0;
			vec_cp_mp(transferme[counter], point(ii));
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		printf("counter mismatch\n");
		br_exit(271);
	}
	
	for (unsigned int ii=0; ii<num_points(); ii++) {
		clear_vec_mp(point(ii));
	}
	free(pts_mp_);
	
	
	num_pts_ = num_good_pts;
	pts_mp_ = transferme;
	
	clear_vec_mp(temp_vec);
	clear_mp(temp);
	return;
}




void witness_set::merge(const witness_set & W_in)
{
	
	//error checking first
	if ( (num_vars_==0) && (W_in.num_vars_!=0) && (num_points()==0)) {
		num_vars_ = W_in.num_vars_;
		num_natty_vars_ = W_in.num_natty_vars_;
	}

	
	
	
	if (W_in.num_natty_vars_ != this->num_natty_vars_) {
		printf("merging two witness sets with differing numbers of natural variables. %d merging set, %d existing\n",
			   W_in.num_natural_variables(), this->num_natural_variables());
		br_exit(95);
	}
	
	//just mindlessly add the linears.  up to user to ensure linears get merged correctly.  no way to know what they want...
	for (unsigned int ii = 0; ii<W_in.num_linears(); ii++) {
		witness_set::add_linear(W_in.linear(ii));
	}
	
	
	for (unsigned int ii = 0; ii<W_in.num_points(); ii++) {
		int is_new = 1;
		vec_mp & in_point = W_in.point(ii);
		for (unsigned int jj = 0; jj<num_points(); jj++){
			vec_mp & curr_point = this->point(jj);
			if ( curr_point->size == in_point->size) {
				if (isSamePoint_inhomogeneous_input(curr_point, in_point)) {
					is_new = 0;
					break;
				}
			}
			
		}
		
		if (is_new==1)
			witness_set::add_point( (in_point) );
	}
	
	
	for (unsigned int ii = 0; ii<W_in.num_patches(); ii++) {
		int is_new = 1;
		
		for (unsigned int jj = 0; jj<this->num_patches(); jj++){
			if ( this->patch(jj)->size ==  W_in.patch(ii)->size) {
				if (isSamePoint_inhomogeneous_input(patch(jj), W_in.patch(ii))) {
					is_new = 0;
					break;
				}
			}
		}
		
		if (is_new==1)
			witness_set::add_patch(W_in.patch(ii));
	}
	
	
	
	return;
}//re: merge_witness_sets
















void witness_set::send(parallelism_config & mpi_config, int target) const
{
    //send all the data to the other party
    

//
//	
//	std::vector< std::string > variable_names;
//	
//	boost::filesystem::path input_filename;
//	function input_file;
//	// end data members
    
    int *buffer = (int *) br_malloc(8*sizeof(int));
    buffer[0] = dimension();
    buffer[1] = component_number();
    buffer[2] = incidence_number();
    buffer[3] = num_variables();
    buffer[4] = num_natural_variables();
    buffer[5] = num_points();
    buffer[6] = num_linears();
    buffer[7] = num_patches();
    
    MPI_Send(buffer, 8, MPI_INT, WITNESS_SET, target, mpi_config.comm());
    
    free(buffer);
    
    for (unsigned int ii=0; ii<num_linears(); ii++) {
        send_vec_mp( linear(ii),target);
    }
    for (unsigned int ii=0; ii<num_patches(); ii++) {
        send_vec_mp( patch(ii),target);
    }
    for (unsigned int ii=0; ii<num_points(); ii++) {
        send_vec_mp( point(ii) ,target);
    }
    
    char * namebuffer = (char *) br_malloc(1024*sizeof(char));
    
	
    free(namebuffer);
    return;
}

void witness_set::receive(int source, parallelism_config & mpi_config)
{
    MPI_Status statty_mc_gatty;
    
	int *buffer = new int[8];
    
    
    MPI_Recv(buffer, 8, MPI_INT, WITNESS_SET, source, mpi_config.comm(), &statty_mc_gatty);
    
    set_dimension(buffer[0]);
    comp_num_ = buffer[1];
    incid_num_ = buffer[2];
    num_vars_ = buffer[3];
    num_natty_vars_ = buffer[4];
    unsigned int temp_num_pts = buffer[5];
    unsigned int temp_num_linears = buffer[6];
    unsigned int temp_num_patches = buffer[7];
    
	delete [] buffer;
	
    vec_mp tempvec; init_vec_mp2(tempvec,0,1024);
    
    for (unsigned int ii=0; ii<temp_num_linears; ii++) {
        receive_vec_mp(tempvec,source);
        add_linear(tempvec);
    }
    
    for (unsigned int ii=0; ii<temp_num_patches; ii++) {
        receive_vec_mp(tempvec,source);
        add_patch(tempvec);
    }
    
	for (unsigned int ii=0; ii<temp_num_pts; ii++) {
        receive_vec_mp(tempvec,source);
        add_point(tempvec);
    }
    
    clear_vec_mp(tempvec);
    return;
}






void witness_data::populate()
{
	FILE *IN;
	

	
	IN = safe_fopen_read("witness_data");
	
	
	int num_nonempty_codims;

	
	fscanf(IN,"%d %d", &num_variables_, &num_nonempty_codims);
	
	
	
	vec_mp prev_approx;
	init_vec_mp(prev_approx,num_variables_); prev_approx->size = num_variables_;
	
	vec_mp temp_vec;
	init_vec_mp(temp_vec,num_variables_); temp_vec->size = num_variables_;
	
	
	//create counters and initialize to zeros
	std::vector<int> component_counter(num_nonempty_codims,0),codim_indicator(num_nonempty_codims,0);
	
	
	for (int ii=0; ii<num_nonempty_codims; ii++) {
		int current_codimension,num_points_this_dim;
		fscanf(IN,"%d %d",&current_codimension,&num_points_this_dim);
		
		
		//std::cout << "getting codimension " << current_codimension << " with " << num_points_this_dim << " points " << std::endl;
		
		codim_indicator[ii] = current_codimension;
		
		int current_dimension = num_variables_-1-current_codimension;
		
		nonempty_dimensions.push_back(current_dimension);
		
		for (int jj=0; jj<num_points_this_dim; jj++) {
			
			
			int precision;
			// last approximation
			fscanf(IN,"%d",&precision);
			change_prec_vec_mp(temp_vec,precision);
			for (int kk=0; kk<num_variables_; kk++) {
				mpf_inp_str(temp_vec->coord[kk].r, IN, 10); // 10 is the base
				mpf_inp_str(temp_vec->coord[kk].i, IN, 10);
			}
			
			
			 
			// previous approximation
			fscanf(IN,"%d",&precision);
			change_prec_vec_mp(temp_vec,precision);
			for (int kk=0; kk<num_variables_; kk++) {
				mpf_inp_str(prev_approx->coord[kk].r, IN, 10); // 10 is the base
				mpf_inp_str(prev_approx->coord[kk].i, IN, 10);
			}
			
			
			witness_point_metadata meta(current_dimension);
			
			meta.set_from_file(IN);
			

			int index = add_solution(temp_vec, meta);
			
			int num_existing_pts_this_comp = map_lookup_with_default( dimension_component_counter[current_dimension], meta.component_number(), 0 ); // the right hand 0 is the default if not found
			
			dimension_component_counter[current_dimension][meta.component_number()] = num_existing_pts_this_comp+1;
			index_tracker[current_dimension][meta.component_number()].push_back(index);
			
			
			if (meta.component_number()+1 > component_counter[ii]) {
				component_counter[ii]  = meta.component_number()+1;  // keeps track of the number of components per dimension
			}
			
			
			
		}
	}
	
	
	
	int should_be_minus_one;
	fscanf(IN,"%d",&should_be_minus_one);
	
	if (should_be_minus_one!=(-1)) {
		std::cerr << "did not parse top of file correctly.  got " << should_be_minus_one << " instead of -1\n";
		br_exit(-1);
	}
	
	int numbertype;
	fscanf(IN,"%d",&numbertype);
	
	if (numbertype==2 || numbertype==1) {
		change_prec_vec_mp(temp_vec,1024); // TODO: make this adapt to settings.
	}
	else {
		change_prec_vec_mp(temp_vec,64);
	}
	
	
	mpq_t *temp_rat = (mpq_t *) br_malloc(2*sizeof(mpq_t)); // create and allocate two rationals.
	
	
	for (int mm=0; mm<num_nonempty_codims; mm++) {
		//BEGIN REPEATING BLOCK, one per nonempty codimension
		int current_codimension;
		current_codimension = codim_indicator[mm];
		int current_dimension = num_variables_-1-current_codimension;
		
		int num_rows_randomization, num_cols_randomization;
		fscanf(IN,"%d %d",&num_rows_randomization,&num_cols_randomization);
		
		mat_mp randomization_matrix;  init_mat_mp(randomization_matrix,num_rows_randomization,num_cols_randomization);
		for (int ii = 0; ii < num_rows_randomization; ii++) {
			for (int jj=0; jj<num_cols_randomization; jj++) {
				if (numbertype==2) {
					setup_comp_in_rat(temp_rat,IN);
					rat_to_mp(&randomization_matrix->entry[ii][jj],temp_rat); clear_rat(temp_rat);
				}
				else{
					mpf_inp_str(randomization_matrix->entry[ii][jj].r, IN, 10);
					mpf_inp_str(randomization_matrix->entry[ii][jj].i, IN, 10);
				}
			}
		}
		clear_mat_mp(randomization_matrix);
		
		
		//MATRIX W FOR HOMOGENIZATION
		//  same length as the randomization matrix
		homogenization_matrix_.resize(num_rows_randomization);
		for (int ii=0; ii<num_rows_randomization; ii++){
			homogenization_matrix_[ii].resize(num_cols_randomization);
			for (int jj=0; jj<num_cols_randomization; jj++){
				fscanf(IN,"%d",&homogenization_matrix_[ii][jj]);}
		
		}
		
		
		int num_entries_hom_patch_eqn;
		fscanf(IN,"%d",&num_entries_hom_patch_eqn);
		vec_mp homogenization_patch_eqn;  init_vec_mp(homogenization_patch_eqn,num_entries_hom_patch_eqn); homogenization_patch_eqn->size = num_entries_hom_patch_eqn;
		//		//VECTOR H FOR HOMOGENIZATION
		for (int ii = 0; ii<num_entries_hom_patch_eqn; ii++) {
			if (numbertype==2) {
				setup_comp_in_rat(temp_rat,IN);
				rat_to_mp(&temp_vec->coord[ii],temp_rat);  clear_rat(temp_rat);
			}
			else{
				mpf_inp_str(temp_vec->coord[ii].r, IN, 10);
				mpf_inp_str(temp_vec->coord[ii].i, IN, 10);
			}
		}
		clear_vec_mp(homogenization_patch_eqn);
		
		
		
		
		//   HOMVARCONST
		comp_mp hom_variable_constant;  init_mp(hom_variable_constant);
		setup_comp_in_rat(temp_rat,IN);
		rat_to_mp(hom_variable_constant,temp_rat); clear_rat(temp_rat);
		clear_mp(hom_variable_constant);
		
		
		
		//MATRIX B FOR LINEAR SLICE COEFFICIENTS
		int num_linears, num_lin_entries;
		fscanf(IN,"%d %d",&num_linears,&num_lin_entries);
		
		change_size_vec_mp(temp_vec,num_lin_entries);
		temp_vec->size = num_lin_entries;
		
		for (int ii = 0; ii < num_linears; ii++) {
			for (int jj=0; jj< num_lin_entries; jj++) {
				if (numbertype==2) {
					setup_comp_in_rat(temp_rat,IN);
					rat_to_mp(&temp_vec->coord[jj],temp_rat); clear_rat(temp_rat);
				}
				else{
					mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
					mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
				}
			}
			
			add_linear_w_meta(temp_vec, witness_linear_metadata(current_dimension));
			
		}
		
		
		
		
		int num_patch_coefficients;
		fscanf(IN,"%d",&num_patch_coefficients);
		
		
		
		change_size_vec_mp(temp_vec,num_patch_coefficients);
		temp_vec->size = num_patch_coefficients;
		
		//		// PATCH COEFFICIENTS
		for (int ii = 0; ii<num_patch_coefficients; ii++) {
			if (numbertype==2) {
				setup_comp_in_rat(temp_rat,IN);
				rat_to_mp(&temp_vec->coord[ii],temp_rat);  clear_rat(temp_rat);
			}
			else{
				mpf_inp_str(temp_vec->coord[ii].r, IN, 10);
				mpf_inp_str(temp_vec->coord[ii].i, IN, 10);
			}
			
		}
		
		add_patch_w_meta(temp_vec, witness_patch_metadata(current_dimension));
		
		//END REPEATED BLOCK (ONE FOR EACH NONEMPTY CODIM).
	}
	
	fclose(IN);
	
	free(temp_rat);
	
	
	
	
	clear_vec_mp(temp_vec);
	clear_vec_mp(prev_approx);
	
	
	
	
	
	
	
	
	
	
	
	
	return;
}






witness_set witness_data::choose(BR_configuration & options)
{
#ifdef functionentry_output
	std::cout << "witness_data::choose" << std::endl;
#endif

	int target_dimension = options.target_dimension;
	int target_component = options.target_component;
	
	
	
	if (target_dimension == -1) {
		//try to get it by inference.  this is default behaviour.
		
		if (nonempty_dimensions.size()==1) {
			
			
			target_dimension = nonempty_dimensions[0];
			if (target_component==-1) { // want all components
				return best_possible_automatic_set(options);
			}
			else if(target_component==-2) // this is default
			{
				if (dimension_component_counter[target_dimension].size()==1) {
					return form_specific_witness_set(target_dimension,0);
				}
				else
				{
					return choose_set_interactive(options);
				}
			}
			else{ // want a specific witness set
				
//				std::cout << map_lookup_with_default(dimension_component_counter[target_dimension],target_component,0) << " " << dimension_component_counter[target_dimension][target_component] << std::endl;
				
				
				if (map_lookup_with_default(dimension_component_counter[target_dimension],target_component,0)==0) {
					std::cout << "you asked for component " << target_component << " (of dimension " << nonempty_dimensions[0] << ") which does not exist" << std::endl;
					witness_set W(num_variables());
					W.set_dimension(nonempty_dimensions[0]);
					W.set_component_number(-1);
					return W;
				}
				else{
					return form_specific_witness_set(nonempty_dimensions[0],target_component);
				}
			}
		}
		else{
			return choose_set_interactive(options);
		}
	}
	else{
		//want a specific dimension
		
		
		if (dimension_component_counter.find(target_dimension)==dimension_component_counter.end()) {
			std::cout << "there are no components of dimension " << options.target_dimension << std::endl;
			
			witness_set W(num_variables());
			W.set_dimension(options.target_dimension);
			W.set_component_number(-1);
			return W;
		}
		else{
			if (target_component==-2) {
				if (dimension_component_counter[target_dimension].size()==1) { // this line is fragile because if something is looked up and doesn't exist, and entry is added.  this is why there is the function map_lookup_with_default.
					return form_specific_witness_set(target_dimension,0);
				}
				else{
					return choose_set_interactive(options);
				}
			}
			else if (target_component==-1) {
				return best_possible_automatic_set(options); // this may eventually call the interactive chooser
			}
			else{
				// want both a specific dimension and component number
				if (map_lookup_with_default(dimension_component_counter[target_dimension],target_component,0)==0) {
					std::cout << "you asked for a component (" << target_component << ") which does not exist" << std::endl;
					witness_set W(num_variables());
					W.set_dimension(target_dimension);
					W.set_component_number(-1);
					return W;
				}
				else{
					return form_specific_witness_set(target_dimension,target_component);
				}
			}
		}
		
		
		
	}
	
	
}













witness_set witness_data::best_possible_automatic_set(BR_configuration & options)
{
#ifdef functionentry_output
	std::cout << "witness_data::best_possible_automatic_set" << std::endl;
#endif
	int target_dimension = options.target_dimension;
	
	witness_set W(num_variables()); // create blank witness set
	W.set_dimension(target_dimension);
	W.set_component_number(-1);
	
	
	if (index_tracker.find(target_dimension)==index_tracker.end()) {
		std::cout << "must return empty set from auto constructor.  have no components of dimension " << target_dimension << std::endl;
		return W;
	}
	
	
	
	std::vector<int> components_with_no_deflations_needed;
	for (auto iter=index_tracker[target_dimension].begin(); iter!=index_tracker[target_dimension].end(); ++iter) {
		// iterate over components for target dimension
		if (iter->second.size()==0) {
			std::cout << "detected a witness set with no points.  this should be impossible...  by definition a witness set has at least one point in it." << std::endl;
			mypause();
		}
		else{
			if (point_metadata[iter->second[0]].num_deflations_needed()==0) {
				components_with_no_deflations_needed.push_back(iter->first);
			}
		}
		
	}
	
	if (components_with_no_deflations_needed.size()==0) { // everything needs deflation
		return choose_set_interactive(options);
	}
	
	std::cout << "checking for self-conjugate components" << std::endl;
	
	int sc_counter = 0;
	// need only copy the points, as all witness sets of a dimension have the same linears and patches.
	for (auto iter=components_with_no_deflations_needed.begin(); iter!=components_with_no_deflations_needed.end(); ++iter) {
		// iterate over deflation-free components for target dimension
		
		int current_index = index_tracker[target_dimension][*iter][0]; // guaranteed to exist, b/c nonempty.  already checked.
		
		if (checkSelfConjugate(point(current_index), options, options.input_filename)==true) {
			std::cout << "dim " << target_dimension << ", comp " << *iter << " is self-conjugate" << std::endl;
			for (int ii=0; ii<dimension_component_counter[target_dimension][*iter]; ++ii) {
				W.add_point( point(index_tracker[target_dimension][*iter][ii]) );
			}
			sc_counter++;
		}
	}
	
	
	
	
	
	if (sc_counter==0) { // found 0 self-conjugate components.
		std::cout << color::green() << "found 0 self-conjugate components, which do not need deflation." << color::console_default() << std::endl;
		return choose_set_interactive(options);
	}
	else{
		
		for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
			if (linear_metadata[ii].dimension() == target_dimension) {
				W.add_linear(linear(ii));
			}
		}
		
		for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
			if (patch_metadata[ii].dimension() == target_dimension) {
				W.add_patch(patch(ii));
			}
		}
		
		
		return W;
	}
		

	
}







// lets the user choose a set, and returns a copy of it.
witness_set witness_data::choose_set_interactive(BR_configuration & options)
{
#ifdef functionentry_output
	std::cout << "witness_data::choose_set_interactive" << std::endl;
#endif
	int target_dimension = options.target_dimension;
	
	std::cout << "the nonempty dimensions:" << std::endl;
	for (auto iter=nonempty_dimensions.begin(); iter!=nonempty_dimensions.end(); ++iter) {
		std::cout << *iter << " ";
	}
	std::cout << std::endl;
	
	
		
	if (target_dimension==-1) { // display all dimensions, which is default
		
		for (auto iter=index_tracker.begin(); iter!=index_tracker.end(); ++iter) {
			std::cout << "dimension " << iter->first << std::endl;
			for (auto jter = iter->second.begin(); jter!=iter->second.end(); ++jter) {
				
				
				int first_index = *(jter->second.begin());
				std::cout << "\tcomponent " << jter->first << ": multiplicity " << point_metadata[first_index].multiplicity() << ", deflations needed: " << point_metadata[first_index].num_deflations_needed() << std::endl;
			}
			std::cout << std::endl;
		}
		
		std::set<int> valid_choices;
		for (auto iter=nonempty_dimensions.begin(); iter!=nonempty_dimensions.end(); ++iter) {
			valid_choices.insert(*iter);
		}
		
		target_dimension = get_int_choice("\n\nchoose a single dimension:\n",valid_choices);
		options.target_dimension = target_dimension;
	}
		
	
	std::cout << "dimension " << target_dimension << std::endl;
	for (auto jter = index_tracker[target_dimension].begin(); jter!= index_tracker[target_dimension].end(); ++jter) {
		int first_index = *(jter->second.begin());
		std::cout << "\tcomponent " << jter->first << ": multiplicity " << point_metadata[first_index].multiplicity() << ", deflations needed: " << point_metadata[first_index].num_deflations_needed() << ", degree " << dimension_component_counter[target_dimension][jter->first] << std::endl;
		
	}
	std::cout << std::endl;
	
	
	
	int num_components = get_int_choice("how many components would you like to decompose? (0 for all mult-one components)\n",0,dimension_component_counter[target_dimension].size());
	
	
	std::set<int> chosen_few;
	
	if (num_components==0) {
		for (unsigned int ii=0; ii<dimension_component_counter[target_dimension].size(); ++ii) {
			if (point_metadata[index_tracker[target_dimension][ii][0]].multiplicity()==1) {
				chosen_few.insert(ii);
			}
			
		}
	}
	else{
		std::set<int> valid_choices;
		for (unsigned int ii=0; ii<dimension_component_counter[target_dimension].size(); ++ii) {
			valid_choices.insert(ii);
		}
		
		for (int ii=0; ii<num_components; ii++) {
			int new_choice = get_int_choice("choose a component:\n",valid_choices);
			chosen_few.insert(new_choice);
			valid_choices.erase(new_choice);
		}
		
	}
	
	
	
	
	//3. return the set via a copy.
	witness_set W(num_variables());
	
	W.set_dimension(target_dimension);
	if (chosen_few.size()==1) {
		W.set_component_number(*chosen_few.begin());
	}
	else{
		W.set_component_number(-1);
	}
	
	for (auto iter=chosen_few.begin(); iter!=chosen_few.end(); ++iter) {
		
		//iterate over each point in the component
		for (auto jter=index_tracker[target_dimension][*iter].begin(); jter!=index_tracker[target_dimension][*iter].end(); ++jter) {
			W.add_point( point(*jter) );
		}
		
	}
	
	for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
		if (linear_metadata[ii].dimension() == target_dimension) {
			W.add_linear(linear(ii));
		}
	}
	
	for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
		if (patch_metadata[ii].dimension() == target_dimension) {
			W.add_patch(patch(ii));
		}
	}
	
	
	return W;
}





witness_set witness_data::form_specific_witness_set(int dim, int comp)
{
	witness_set W(num_variables());
	
	W.set_dimension(dim);
	W.set_component_number(comp);
	
	
	// check to make sure actually have the set we need.
	if (dimension_component_counter.find(dim)==dimension_component_counter.end()) {
		std::cout << "trying to construct witness set of dimension " << dim << " but there are no components of that dimension." << std::endl;
		return W;
	}
	
	
	if (dimension_component_counter[dim].find(comp) == dimension_component_counter[dim].end()) {
		std::cout << "trying to constuct witness set for dimension " << dim << ", component " << comp << ", but that component does not exist" << std::endl;
		return W;
	}
	
	// ok, it exists.  lets form it.
	
	for (auto iter = index_tracker[dim][comp].begin(); iter!= index_tracker[dim][comp].end(); ++iter) {
		//iter points to an index into the vertices stored in the vertex set.
		W.add_point( point(*iter) );
	}
	
	
	for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
		if (linear_metadata[ii].dimension() == dim) {
			W.add_linear(linear(ii));
		}
	}
	
	for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
		if (patch_metadata[ii].dimension() == dim) {
			W.add_patch(patch(ii));
		}
	}
	
	
	return W;
}


void vertex::set_point(const vec_mp new_point)
{
	change_prec_vec_mp(this->pt_mp_, new_point->curr_prec);
	vec_cp_mp(this->pt_mp_, new_point);
}




void vertex::send(int target, parallelism_config & mpi_config)
{
	
	send_vec_mp(pt_mp_, target);
	
	send_vec_mp(projection_values_, target);
	
	int * buffer = (int *) br_malloc(3*sizeof(int));
	buffer[0] = type_;
	buffer[1] = removed_;
	buffer[2] = input_filename_index_;
	
	MPI_Send(buffer, 3, MPI_INT, target, VERTEX, mpi_config.comm());
	free(buffer);
	
}


void vertex::receive(int source, parallelism_config & mpi_config)
{
	MPI_Status statty_mc_gatty;
	int * buffer = (int *) br_malloc(3*sizeof(int));
	
	
	receive_vec_mp(pt_mp_, source);
	receive_vec_mp(projection_values_, source);
	
	MPI_Recv(buffer, 3, MPI_INT, source, VERTEX, mpi_config.comm(), &statty_mc_gatty);
	
	type_ = buffer[0];
	removed_ = buffer[1];
	input_filename_index_ = buffer[2];
//	print_point_to_screen_matlab(pt_mp,"recvpt");
	free(buffer);
}




int vertex_set::search_for_point(vec_mp testpoint)
{
    int index = -1;
	
    index = search_for_active_point(testpoint);
    
    if (index==-1) {
        index = search_for_removed_point(testpoint);
    }
    
    return index;
}


int vertex_set::search_for_active_point(vec_mp testpoint)
{
    int index = -1; // initialize the index we will return

    // dehomogenize the testpoint into the internal temp container.
    for (int jj=1; jj<num_natural_variables_; jj++) {
		div_mp(&checker_1_->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
	}
    
    
	//	WTB: a faster comparison search.
	for (unsigned int ii=0; ii<num_vertices_; ii++) {
		
		int current_index = ii;
		
		if (vertices_[current_index].is_removed()==0) {
			vec_mp& current_point = vertices_[current_index].point();
			
			// dehomogenize the current point under investigation
			for (int jj=1; jj<num_natural_variables_; jj++) {
				div_mp(&checker_2_->coord[jj-1], &(current_point)->coord[jj], &(current_point)->coord[0]);
			}
			
			if (isSamePoint_inhomogeneous_input(checker_1_, checker_2_)){
				index = current_index;
				break;
			}
			
			if (index!=-1)
				break;
			
		}
	}
    
    return index;
}



int vertex_set::search_for_removed_point(vec_mp testpoint)
{
    
    int index = -1; // initialize the index we will return
    
    
    // dehomogenize the testpoint into the internal temp container.
    
    for (int jj=1; jj<num_natural_variables_; jj++) {
		div_mp(&checker_1_->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
	}
    
    
    //	WTB: a faster comparison search.
    for (unsigned int ii=0; ii<num_vertices_; ii++) {
		int current_index = ii;
		
		if (vertices_[current_index].is_removed()) {
			
			vec_mp & current_point = vertices_[current_index].point();
			
			for (int jj=1; jj<num_natural_variables_; jj++) {
				div_mp(&checker_2_->coord[jj-1], &(current_point)->coord[jj], &(current_point)->coord[0]);
			}
			
			if (isSamePoint_inhomogeneous_input(checker_1_, checker_2_)){
				index = current_index;
				break;
			}
			
			if (index!=-1)
				break;
			
		}
	}
    
    return index;
}





int vertex_set::compute_downstairs_crit_midpts(const witness_set & W,
                                               vec_mp crit_downstairs,
                                               vec_mp midpoints_downstairs,
                                               std::vector< int > & index_tracker,
                                               vec_mp pi)
{
    
    
	
	int retVal = SUCCESSFUL;
	
    int proj_index = this->get_proj_index(pi);
    
	vec_mp projection_values; init_vec_mp2(projection_values, int(W.num_points()) ,1024);
	projection_values->size = W.num_points();
	
	for (unsigned int ii=0; ii<W.num_points(); ii++){
		vec_mp & curr_point = W.point(ii);
		
		for (int jj=0; jj<curr_point->size; jj++) {
			
			if (!(mpfr_number_p(curr_point->coord[jj].r) && mpfr_number_p(curr_point->coord[jj].i))) {
				std::cout << color::red();
				std::cout << "there was NAN in a coordinate for a point in the projections to sort :(" << std::endl;
				print_point_to_screen_matlab(curr_point, "bad_point");
				
				print_point_to_screen_matlab(pi,"pi");
				std::cout << color::console_default();
				return CRITICAL_FAILURE;
			}
			
		}
		
        
        int curr_index = search_for_point(curr_point);
        
        if (curr_index < 0) {
            std::cout << color::red() << "trying to retrieve projection value from a non-stored point" << color::console_default() << std::endl;
            mypause();
        }
        

        
        set_mp(&projection_values->coord[ii], & (vertices_[curr_index].projection_values())->coord[proj_index]);
        
		
		
		if (!(mpfr_number_p(projection_values->coord[ii].r) && mpfr_number_p(projection_values->coord[ii].i)))
		{
			print_comp_matlab(&projection_values->coord[ii],"not_a_number");
			print_point_to_screen_matlab(pi,"pi");
			
		}
		
	}
	
	
#ifdef thresholding
	real_threshold(projection_values, 1e-10);
#endif
	
	
	change_size_vec_mp(crit_downstairs,1); // destructive resize
	crit_downstairs->size = 1;
	
	retVal = sort_increasing_by_real(crit_downstairs, index_tracker, projection_values);
	
	clear_vec_mp(projection_values); // done with this data.  clear it.
	
	int num_midpoints = crit_downstairs->size - 1;
	
	if (num_midpoints<1) {
        // work is done, simply return
		return retVal;
	}
	
	
	comp_d h;  h->r = 0.5; h->i = 0.0;
	comp_mp half; init_mp2(half,1024); d_to_mp(half, h);
	
	change_size_vec_mp(midpoints_downstairs, num_midpoints);
	midpoints_downstairs->size = num_midpoints;
	
	comp_mp temp; init_mp2(temp,1024);
	for (int ii=0; ii<num_midpoints; ii++){
		add_mp(temp, &crit_downstairs->coord[ii], &crit_downstairs->coord[ii+1]);
		mul_mp(&midpoints_downstairs->coord[ii], temp, half);
	}
	
	clear_mp(temp);
	clear_mp(half);
	
	
	return retVal;
}




std::vector<int> vertex_set::assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value)
{
    if (this->curr_projection_<0) {
        std::cout << color::red() << "trying to assert projection value (current index) without having index set" << color::console_default() << std::endl;
        br_exit(-91621);
    }
    return assert_projection_value(relevant_indices,new_value,this->curr_projection_);
}

std::vector<int> vertex_set::assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value, int proj_index)
{
	std::vector<int> bad_indices;
    
    comp_mp temp; init_mp(temp);
    
    if ( proj_index > num_projections_ ) {
        std::cout << color::red() << "trying to assert projection value, but index of projection is larger than possible" << color::console_default() << std::endl;
		br_exit(3091); // throw?
    }
    
    for (std::set<int>::iterator ii=relevant_indices.begin(); ii!=relevant_indices.end(); ii++) {
        //*ii
        
        sub_mp(temp, &(vertices_[*ii].projection_values())->coord[proj_index], new_value);
        if (fabs(mpf_get_d(temp->r))>0.0001) {
            std::cout << "trying to assert projection value of " << mpf_get_d(new_value->r)
				      << " but original value is " << mpf_get_d((vertices_[*ii].projection_values())->coord[proj_index].r) << std::endl;
            std::cout << "point index is " << *ii << std::endl;
			bad_indices.push_back(*ii);
			continue;
        }
        
        set_mp(&(vertices_[*ii].projection_values())->coord[proj_index], new_value);
    }
    
    
    clear_mp(temp);
    return bad_indices;
}





int vertex_set::add_vertex(const vertex source_vertex)
{
	
	
//	if (this->num_projections <= 0) {
//		std::cout << color::red() << "vertex set has no projections!!!" << color::console_default() << std::endl;
//		br_exit(-9644);
//	}
//	
//	if (curr_input_index<0 && source_vertex.input_filename_index == -1) {
//		std::cout << color::red() << "adding points to vertex set, but input file index unset." << color::console_default() << std::endl;
//		br_exit(6711);
//	}
	
	
	vertices_.push_back(source_vertex);
	
	
	if ((vertices_[num_vertices_].projection_values())->size < num_projections_) {
		increase_size_vec_mp((vertices_[num_vertices_].projection_values()), num_projections_ );
		(vertices_[num_vertices_].projection_values())->size = num_projections_;
	}
	
	for (int ii=0; ii<this->num_projections_; ii++){
		bool compute_proj_val = false;
		
		if (! (mpfr_number_p( (vertices_[num_vertices_].projection_values())->coord[ii].r) && mpfr_number_p( (vertices_[num_vertices_].projection_values())->coord[ii].i)  ) )
		{
			compute_proj_val = true;
		}
		//		else if ( false )//yeah, i dunno what else right yet.
		//		{
		//			compute_proj_val = true;
		//		}
		
		
		if (compute_proj_val==true) {
			projection_value_homogeneous_input(&(vertices_[num_vertices_].projection_values())->coord[ii],
											   vertices_[num_vertices_].point(),
											   projections_[ii]);
#ifdef thresholding
            real_threshold(&(vertices_[num_vertices_].projection_values())->coord[ii],1e-13);
#endif
		}
		
	}
	
    
	
	if (vertices_[num_vertices_].input_filename_index() == -1)
	{
		vertices_[num_vertices_].set_input_filename_index(curr_input_index_);
	}
		
	
	this->num_vertices_++;
	return this->num_vertices_-1;
}


void vertex_set::print_to_screen()
{
	printf("vertex set has %zu vertices:\n\n",num_vertices_);
	for (unsigned int ii=0; ii<this->num_vertices_; ++ii) {
		print_point_to_screen_matlab(vertices_[ii].point(),"vert");
		print_point_to_screen_matlab(vertices_[ii].projection_values(),"projection_values");
		printf("type: %d\n", vertices_[ii].type());
	}
}


int vertex_set::setup_vertices(boost::filesystem::path INfile)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	unsigned int temp_num_vertices;
	int num_vars;
	int tmp_num_projections;
	int tmp_num_filenames;
	fscanf(IN, "%u %d %d %d\n\n", &temp_num_vertices, &tmp_num_projections, &num_natural_variables_, &tmp_num_filenames);
	
	
	vec_mp temp_vec; init_vec_mp2(temp_vec,num_natural_variables_,1024);
	for (int ii=0; ii<tmp_num_projections; ii++) {
		for (int jj=0; jj<num_natural_variables_; jj++) {
			mpf_inp_str(temp_vec->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vec->coord[jj].i, IN, 10);
		}
		add_projection(temp_vec);
	}
	clear_vec_mp(temp_vec);
	
	scanRestOfLine(IN);
	scanRestOfLine(IN);
	
	
	for (int ii=0; ii<tmp_num_filenames; ii++) {
		int tmp_size;
		fscanf(IN,"%d\n",&tmp_size);
		
		char * buffer = new char[tmp_size];
		fgets(buffer, tmp_size, IN);
		boost::filesystem::path temppath = buffer;
		this->filenames_.push_back(temppath);
		
		delete [] buffer;
	}
	
	
	
	vertex temp_vertex;
	
	for (unsigned int ii=0; ii<temp_num_vertices; ii++)
	{
		fscanf(IN, "%d\n", &num_vars);
		if ((temp_vertex.point())->size != num_vars) {
			change_size_vec_mp(temp_vertex.point(),num_vars); (temp_vertex.point())->size = num_vars;
		}
		
		for (int jj=0; jj<num_vars; jj++)
		{
			mpf_inp_str((temp_vertex.point())->coord[jj].r, IN, 10);
			mpf_inp_str((temp_vertex.point())->coord[jj].i, IN, 10);
		}
		
		int temp_num;
		fscanf(IN,"%d\n",&temp_num);
		increase_size_vec_mp(temp_vertex.projection_values(),temp_num);
		(temp_vertex.projection_values())->size = temp_num;
		for (int jj=0; jj<temp_num; jj++) {
			mpf_inp_str((temp_vertex.projection_values())->coord[jj].r, IN, 10);
			mpf_inp_str((temp_vertex.projection_values())->coord[jj].i, IN, 10);
		}
		
		int temp_int;
		fscanf(IN,"%d\n",&temp_int);
		temp_vertex.set_input_filename_index(temp_int);
		
		fscanf(IN,"%d\n",&temp_int);
	   temp_vertex.set_type(temp_int);
		
		vertex_set::add_vertex(temp_vertex);
	}
	
	
	
	fclose(IN);
	
	if (this->num_vertices_!=temp_num_vertices) {
		printf("parity error in num_vertices.\n\texpected: %zu\tactual: %u\n",num_vertices_,temp_num_vertices); // this is totally impossible.
		br_exit(25943);
	}
	
	return num_vertices_;
}







void vertex_set::print(boost::filesystem::path outputfile) const
{
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%zu %d %d %lu\n\n",num_vertices_,num_projections_, num_natural_variables_, filenames_.size());
	
	
	
	for (int ii=0; ii<num_projections_; ii++) {
		for (int jj=0; jj<num_natural_variables_; jj++) {
			print_mp(OUT, 0, &projections_[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	
	
	for (unsigned int ii=0; ii!=filenames_.size(); ii++) {
		int strleng = filenames_[ii].string().size() + 1; // +1 for the null character
		char * buffer = new char[strleng];
		memcpy(buffer, filenames_[ii].c_str(), strleng);
		fprintf(OUT,"%d\n",strleng);
		fprintf(OUT,"%s\n",buffer);
		delete [] buffer;
	}
	
	for (unsigned int ii = 0; ii < num_vertices_; ii++)
	{ // output points
		fprintf(OUT,"%d\n", (vertices_[ii].get_point())->size);
		for(int jj=0;jj<(vertices_[ii].get_point())->size;jj++) {
			print_mp(OUT, 0, &(vertices_[ii].get_point())->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"%d\n",(vertices_[ii].get_projection_values())->size);
		for(int jj=0;jj<(vertices_[ii].get_projection_values())->size;jj++) {
			print_mp(OUT, 0, &(vertices_[ii].get_projection_values())->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"%d\n",vertices_[ii].input_filename_index());
		
		fprintf(OUT,"\n");
		if (vertices_[ii].is_removed()) {
			fprintf(OUT,"%d\n\n",REMOVED);
		}
		else{
			fprintf(OUT,"%d\n\n",vertices_[ii].type());
		}
	}
	
	
	
	fclose(OUT);
	
}





void vertex_set::send(int target, parallelism_config & mpi_config)
{
	

	int num_filenames = filenames_.size();
	
	
	int * buffer2 = new int[6];
	buffer2[0] = num_natural_variables_;
	buffer2[1] = num_projections_;
	buffer2[2] = curr_projection_;
	buffer2[3] = num_filenames;
	buffer2[4] = curr_input_index_;
	buffer2[5] = num_vertices_;
	
	MPI_Send(buffer2, 6, MPI_INT, target, VERTEX_SET, MPI_COMM_WORLD);
	
	delete [] buffer2;
//	std::cout << "sending " << num_projections << " projections" << std::endl;
	
	buffer2 = new int[num_projections_];
	for (int ii=0; ii<num_projections_; ii++) {
		buffer2[ii] = projections_[ii]->size;
	}
	MPI_Send(buffer2, num_projections_, MPI_INT, target, VERTEX_SET, MPI_COMM_WORLD);
	delete [] buffer2;
	
	for (int ii=0; ii<num_projections_; ii++) {
		send_vec_mp(projections_[ii],target);
	}

//	std::cout << "sending " << num_filenames << " filenames" << std::endl;
	for (int ii=0; ii<num_filenames; ii++) {
		char * buffer;
		
		
		int strleng = filenames_[ii].string().size()+1;
		buffer = new char[strleng];
		memcpy(buffer, filenames_[ii].string().c_str(), strleng-1); // this sucks
		buffer[strleng-1] = '\0';
		
		MPI_Send(&strleng, 1, MPI_INT, target, VERTEX_SET, MPI_COMM_WORLD);
//		std::cout << "sending filename length " << strleng << " " << filenames[ii].string() << std::endl;
		MPI_Send(&buffer[0], strleng, MPI_CHAR, target, VERTEX_SET, MPI_COMM_WORLD);
		
		delete [] buffer;
		
	}
	
	
	for (unsigned int ii=0; ii<num_vertices_; ii++) {
		vertices_[ii].send(target, mpi_config);
	}
	
	
	

	
	
	
	
	return;
}


void vertex_set::receive(int source, parallelism_config & mpi_config)
{
	MPI_Status statty_mc_gatty;
	
	int * buffer2 = new int[6];
	MPI_Recv(buffer2, 6, MPI_INT, source, VERTEX_SET, MPI_COMM_WORLD, &statty_mc_gatty);
	
	
	
	int temp_num_natural_variables = buffer2[0];
	int temp_num_projections = buffer2[1];
	curr_projection_ = buffer2[2];
	int temp_num_filenames = buffer2[3];
	curr_input_index_ = buffer2[4];
	unsigned int temp_num_vertices = buffer2[5];
	
	delete [] buffer2;
	
	set_num_vars(temp_num_natural_variables);
	
	
//	std::cout << "receiving " << temp_num_projections << " projections" << std::endl;
	
	buffer2 = new int[temp_num_projections];
	
	MPI_Recv(buffer2, temp_num_projections, MPI_INT, source, VERTEX_SET, MPI_COMM_WORLD, &statty_mc_gatty);

	
	
	vec_mp tempvec; init_vec_mp2(tempvec, 0, 1024);
	for (int ii=0; ii<temp_num_projections; ii++) {
//		std::cout << "recving " << ii << "th proj" << std::endl;
		change_size_vec_mp(tempvec,buffer2[ii]); tempvec->size = buffer2[ii];
		receive_vec_mp(tempvec,source);
		add_projection(tempvec);
//		print_point_to_screen_matlab(tempvec,"tempvec_recvd_proj");
		
	}
	clear_vec_mp(tempvec);
	delete [] buffer2;
	
	if (num_projections_!=temp_num_projections) {
		std::cout << "num_projections doesn't match!" << std::endl;
	}
	
//	std::cout << "receiving " << temp_num_filenames << " filenames" << std::endl;
	for (int ii=0; ii<temp_num_filenames; ii++) {
		char * buffer; int strleng;
		
		MPI_Recv(&strleng, 1, MPI_INT, source, VERTEX_SET, MPI_COMM_WORLD, &statty_mc_gatty);
		
		buffer = new char[strleng];
//		std::cout << "recving filename length " << strleng << std::endl;
		MPI_Recv(&buffer[0], strleng, MPI_CHAR, source, VERTEX_SET, MPI_COMM_WORLD, &statty_mc_gatty);
		filenames_.push_back(boost::filesystem::path(std::string(buffer)));
		
		delete [] buffer;
		
	}
	
	
	
	

	
	for (unsigned int ii=0; ii<temp_num_vertices; ii++) {
		vertex tempvert;
		tempvert.receive(source, mpi_config);
		add_vertex(tempvert);
	}
	
	if (num_vertices_ != temp_num_vertices) {
		std::cout << "logical inconsistency.  do not have correct num vertices." << std::endl;
	}
	
	
	

	
	return;
}










































int decomposition::add_witness_set(const witness_set & W, int add_type, vertex_set & V)
{
#ifdef functionentry_output
	std::cout << "decomposition::add_witness_set" << std::endl;
#endif
	
	
    V.set_curr_input(W.input_filename());
    
    vertex temp_vertex;
    temp_vertex.set_type(add_type);
    
    for (unsigned int ii=0; ii<W.num_points(); ii++) {
        vec_cp_mp(temp_vertex.point(), W.point(ii));
        this->index_in_vertices_with_add(V, temp_vertex);
    }
    
    return 0;
}


int decomposition::add_vertex(vertex_set & V, vertex source_vertex)
{
#ifdef functionentry_output
	std::cout << "decomposition::add_vertex" << std::endl;
#endif
	
	int current_index = V.add_vertex(source_vertex);
	
	
	return current_index;
}






int decomposition::index_in_vertices(vertex_set & V,
									 vec_mp testpoint)
{
#ifdef functionentry_output
	std::cout << "decomposition::index_in_vertices" << std::endl;
#endif
	int index = -1;
	
	
	// first we search the non-removed points.
	index = V.search_for_point(testpoint);
    
    
	return index;
}





//TODO: make T here a pointer

int decomposition::index_in_vertices_with_add(vertex_set &V,
											  vertex vert)
{
	int index = decomposition::index_in_vertices(V, vert.point());
	
	if (index==-1) {
		index = decomposition::add_vertex(V, vert);
	}
	
	return index;
	
}











int decomposition::setup(boost::filesystem::path INfile)
{
	boost::filesystem::path directoryName = INfile.parent_path();

	std::stringstream converter;
	std::string tempstr;
	std::ifstream fin(INfile.c_str());
	
	
	getline(fin, tempstr);
	set_input_filename(directoryName / tempstr);
	
	getline(fin, tempstr);
	converter << tempstr;
	converter >> this->num_variables_ >> this->dim_;
	converter.clear(); converter.str("");
	


	vec_mp tempvec; init_vec_mp2(tempvec, this->num_variables_,1024);
	tempvec->size = this->num_variables_;
	

	for (int ii=0; ii<dimension(); ii++) {
		int temp_size;
		fin >> temp_size;
		
		change_size_vec_mp(tempvec,temp_size);  tempvec->size = temp_size;
		for (int jj=0;jj<temp_size;jj++)
		{
			std::string re, im;
			fin >> re >> im ;
			mpf_set_str(tempvec->coord[jj].r, const_cast<char *>(re.c_str()), 10);
			mpf_set_str(tempvec->coord[jj].i, const_cast<char *>(im.c_str()), 10);
		}
		
		decomposition::add_projection(tempvec);
	}
	
	
	
	
	int curr_num_patches = -1;
	fin >> curr_num_patches;
	
	
	vec_mp temp_patch; init_vec_mp2(temp_patch,1,1024); temp_patch->size = 1;
	for (int ii=0; ii<curr_num_patches; ii++) {
		

		int curr_size = -1;
		fin >> curr_size;

		change_size_vec_mp(temp_patch,curr_size); temp_patch->size = curr_size;
		
		for (int jj=0; jj<curr_size; jj++) {
			std::string re, im;
			fin >> re >> im; // this line is correct
		
			mpf_set_str(temp_patch->coord[jj].r, const_cast<char *>(re.c_str()), 10);
			mpf_set_str(temp_patch->coord[jj].i, const_cast<char *>(im.c_str()), 10);
		}
		
		decomposition::add_patch(temp_patch);
	}
	
	
	clear_vec_mp(tempvec);
	clear_vec_mp(temp_patch);
	
	fin.close();
	
	
	
	return 0;
}






void decomposition::print(boost::filesystem::path base)
{
	
#ifdef functionentry_output
	std::cout << "decomposition::print" << std::endl;
#endif
	
	if (dimension() != num_curr_projections()) {
//		std::cout << "decomposition was short projections\nneeded	" << this->dimension << " but had " << num_curr_projections << std::endl;;
	}
	
	
	int ii;
	
	boost::filesystem::create_directory(base);
	
	FILE *OUT = safe_fopen_write(base / "decomp");
	
	fprintf(OUT,"%s\n",input_filename().filename().c_str());
	
	fprintf(OUT,"%d %d\n\n",num_variables(), dimension());
	
	
	for (ii=0; ii<num_curr_projections(); ii++) {
		fprintf(OUT,"%d\n",pi_[ii]->size);
		for(int jj=0;jj<pi_[ii]->size;jj++)
		{
			print_mp(OUT, 0, &pi_[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	fprintf(OUT,"%zu\n\n",this->num_patches()); // print the header line
	
	for (unsigned ii=0; ii<this->num_patches(); ++ii) {
		vec_mp & curr_patch = patch(ii);
		fprintf(OUT,"%d\n",curr_patch->size);
		for (int jj=0; jj<curr_patch->size; jj++) {
			print_mp(OUT, 0, &curr_patch->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
    
    fprintf(OUT,"\n\n");
    
    print_mp(OUT, 0, this->sphere_radius_);
    fprintf(OUT, "\n%d\n",this->sphere_center_->size);
    
    for (int jj=0; jj<this->sphere_center_->size; jj++) {
        print_mp(OUT, 0, &this->sphere_center_->coord[jj]);
        fprintf(OUT, "\n");
    }
    fprintf(OUT,"\n");
    
	fclose(OUT);
	
	
}




int decomposition::read_sphere(const boost::filesystem::path & bounding_sphere_filename)
{
#ifdef functionentry_output
	std::cout << "decomposition::read_sphere" << std::endl;
#endif
	
	if (num_variables_<2) {
		std::cout << "during read of sphere, decomposition of dimension	" << dimension() << " has " << num_variables() << " variables!" << std::endl;
		mypause();
	}

	change_size_vec_mp(this->sphere_center_, num_variables()-1); //destructive resize
	sphere_center_->size = num_variables()-1;
	
	
	FILE *IN = safe_fopen_read(bounding_sphere_filename);
	
	mpf_inp_str(sphere_radius_->r, IN, 10);
	mpf_set_str(sphere_radius_->i,"0",10);
	
	
	for (int jj=1; jj<num_variables(); jj++) {
		mpf_inp_str(sphere_center_->coord[jj-1].r, IN, 10);
		mpf_set_str(sphere_center_->coord[jj-1].i,"0",10);
	}
	
	//TODO: check this line for validity
	int tmp_num_vars;
	fscanf(IN,"%d",&tmp_num_vars); scanRestOfLine(IN);
	
	
	fclose(IN);
	
	
	have_sphere_ = true;
	
	
	return SUCCESSFUL;
}



void decomposition::compute_sphere_bounds(const witness_set & W_crit)
{
	
#ifdef functionentry_output
	std::cout << "decomposition::compute_sphere_bounds" << std::endl;
#endif

	
	int num_vars = W_crit.num_natural_variables()-1;
	
	change_size_vec_mp(this->sphere_center_, num_vars); //destructive resize
	sphere_center_->size = num_vars;
	
	if (W_crit.num_points() == 0) {
		set_zero_mp(sphere_radius_);
		mpf_set_str(sphere_radius_->r, "3.0",10);
		for (int ii=0; ii<num_vars; ii++) {
			set_zero_mp(&sphere_center_->coord[ii]);
		}
		this->have_sphere_ = true;
		return;
	}
	
	vec_mp(temp_vec); init_vec_mp2(temp_vec,0,1024);
	if (W_crit.num_points()==1)
	{
		set_zero_mp(sphere_radius_);
		mpf_set_str(sphere_radius_->r, "3.0",10);
		dehomogenize(&temp_vec, W_crit.point(0));
		
		for (int ii=0; ii<num_vars; ii++) {
			set_mp(&sphere_center_->coord[ii], &temp_vec->coord[ii]);
		}
#ifdef thresholding
		real_threshold(sphere_center_, 1e-13);
#endif
		clear_vec_mp(temp_vec);
		this->have_sphere_ = true;
		return;
	}
	
	
	//	W_crit.print_to_screen();
	
	
	
	
	comp_mp temp_rad; init_mp2(temp_rad,1024);
	set_zero_mp(temp_rad);
	
	set_one_mp(sphere_radius_);
	neg_mp(sphere_radius_,sphere_radius_); // set to impossibly low value.
	
	
	comp_mp temp_mp;  init_mp2(temp_mp,1024);
	
	
	vec_mp(cumulative_sum); init_vec_mp2(cumulative_sum,num_vars,1024);
	cumulative_sum->size = num_vars;
	
	for (int ii=0; ii<num_vars; ii++) {
		set_zero_mp(&cumulative_sum->coord[ii]);
	}
	
	
	for (unsigned int ii=0; ii<W_crit.num_points(); ii++) {
		dehomogenize(&temp_vec, W_crit.point(ii), num_vars+1);
		temp_vec->size = num_vars;
		add_vec_mp(cumulative_sum, cumulative_sum, temp_vec);
	}
	
	set_zero_mp(temp_mp);
	mpf_set_d(temp_mp->r, double(W_crit.num_points()));
	
	
	
	for (int ii=0; ii<num_vars; ii++) {
		div_mp(&sphere_center_->coord[ii], &cumulative_sum->coord[ii], temp_mp);
	}
	
	
	
	for (unsigned int ii=0; ii<W_crit.num_points(); ii++) {
		dehomogenize(&temp_vec, W_crit.point(ii), num_vars+1);
		temp_vec->size = num_vars;
		vec_sub_mp(temp_vec, temp_vec, sphere_center_);
		
		
		twoNormVec_mp(temp_vec, temp_mp);
		mpf_abs_mp(temp_rad->r, temp_mp);
		
		
		if (mpf_cmp(sphere_radius_->r, temp_rad->r) < 0){
			set_mp(sphere_radius_, temp_rad);
		}
	}
	
	
	
	mpf_set_str(temp_rad->r,"2.0",10);
	mpf_set_str(temp_rad->i,"0.0",10);
	mul_mp(sphere_radius_,temp_rad,sphere_radius_);  // double the radius to be safe.
	
	
	clear_mp(temp_mp); clear_mp(temp_rad);
	clear_vec_mp(temp_vec); clear_vec_mp(cumulative_sum);
	
	
	
	this->have_sphere_ = true;
#ifdef thresholding
	real_threshold(sphere_center_,1e-13);
#endif

	return;
}



void decomposition::copy_data_from_witness_set(const witness_set & W)
{
	// set some member information.
	set_input_filename(W.input_filename());
	set_num_variables(W.num_variables());
	set_component_number(W.component_number());
	
	this->W_ = W;
}

void decomposition::output_main(const boost::filesystem::path base)
{
#ifdef functionentry_output
	std::cout << "decomposition::output_main" << std::endl;
#endif

	FILE *OUT;
	boost::filesystem::path backupdir = base;
	backupdir += "_bak";
	if (boost::filesystem::exists(base)) {
		
		if (boost::filesystem::exists(backupdir)) {
			boost::filesystem::remove_all(backupdir);
		}
		boost::filesystem::rename(base, backupdir);
	}
	boost::filesystem::create_directory(base);
	
	
	copyfile("witness_data",base / "witness_data"); // this is wrong for nested decompositions.
	
	W_.print_to_file(base / "witness_set");
	
	
// TODO:  this should be a write call, not a copy!
	if (input_filename().filename().string().compare("unset")) {
		copyfile(input_filename(), base / input_filename().filename());
	}
	else{
//		std::cout << "not copying inputfile because name was unset -- " << input_filename << std::endl;
	}
	
	
	this->print(base); // using polymorphism and virtualism here!
	
	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%s\n",base.c_str());
	fprintf(OUT,"%d\n",2);//remove this
	fprintf(OUT,"%d\n",dimension());
	fclose(OUT);
	
	
//	if (boost::filesystem::exists(backupdir)) {
//		boost::filesystem::remove_all(backupdir);
//	}
}



void decomposition::send(int target, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "decomposition::send" << std::endl;
#endif

	
	int * buffer2;
	
	
	//pack and send numbers of things.
	buffer2 = new int[7];
	buffer2[0] = num_variables();
	buffer2[1] = dimension();
	buffer2[2] = component_number();
	buffer2[3] = num_curr_projections();

	buffer2[4] = num_patches();
	buffer2[5] = int(have_sphere_);
	int strleng = input_filename().string().size() + 1;
	buffer2[6] = strleng;

	MPI_Send(buffer2, 7, MPI_INT, target, 6, MPI_COMM_WORLD);
	delete [] buffer2;
	
	
	
	
	
	if (num_curr_projections()>0) {
		for (int ii=0; ii<num_curr_projections(); ii++) {
			send_vec_mp(pi_[ii],target);
		}
	}
	
	
	
	
	if ( num_patches()>0) {

		for (unsigned int ii=0; ii<num_patches(); ii++) {
			send_vec_mp(patch(ii),target);
		}
	}
	
	
	
	
	if (have_sphere_) {

		send_vec_mp(sphere_center_,target);
		send_comp_mp(sphere_radius_,target);
	}
	
	
	
	if (strleng>1) {
		char * buffer = new char[strleng];
		memcpy(buffer, input_filename().c_str(), strleng);
		MPI_Send(buffer, strleng, MPI_CHAR, target, 7, MPI_COMM_WORLD);
		delete [] buffer;
	}
	
	
	
	
	randomizer_->send(target,mpi_config);



	
	
	
	return;
}



void decomposition::receive(int source, parallelism_config & mpi_config)
{
#ifdef functionentry_output
	std::cout << "decomposition::receive" << std::endl;
#endif

	
	
	MPI_Status statty_mc_gatty;
	
	
	
	vec_mp tempvec;  init_vec_mp2(tempvec, 0, 1024);
	
	int * buffer2;
	
	
	
	
	
	
	buffer2 = new int[7];
	
	MPI_Recv(buffer2, 7, MPI_INT, source, 6, MPI_COMM_WORLD, &statty_mc_gatty);
	set_num_variables(buffer2[0]);
	dim_ = buffer2[1];
	comp_num_ = buffer2[2];
	int temp_num_projections = buffer2[3];
	int temp_num_patches = buffer2[4];
	have_sphere_ = bool(buffer2[5]);
	int strleng = buffer2[6];

	delete [] buffer2;
	
	
//	std::cout << "receieved:" << std::endl;
//	std::cout << num_variables << std::endl;
//	std::cout << dimension << std::endl;
//	std::cout << component_num << std::endl;
//	std::cout << temp_num_projections << std::endl;
//	std::cout << num_rand_degrees << std::endl;
//	std::cout << rand_rows << std::endl;
//	std::cout << rand_cols << std::endl;
//	std::cout << temp_num_patches << std::endl;
//	std::cout << have_sphere_radius << std::endl;
//	std::cout << strleng << std::endl;
//	std::cout << num_counters << std::endl;
//	std::cout << num_indices << std::endl;
	
	

	
	
	
	if (temp_num_projections>0) {

		for (int ii=0; ii<temp_num_projections; ii++) {
			receive_vec_mp(tempvec,source);
			add_projection(tempvec);
		}
	}
	
	
	


	
	
	

	
	
	if (temp_num_patches>0) {
		for (int ii=0; ii<temp_num_patches; ii++) {
			receive_vec_mp(tempvec,source);
			add_patch(tempvec);
		}
	}
	
	
	
	if (have_sphere_) {
		receive_vec_mp(sphere_center_,source);
		receive_comp_mp(sphere_radius_,source);
	}
	
	if (strleng>1) {
		char * buffer = new char[strleng];
		MPI_Recv(buffer, strleng, MPI_CHAR, source, 7, MPI_COMM_WORLD, &statty_mc_gatty);
		
		set_input_filename(buffer);
		delete [] buffer;
	}
	
	

	

	clear_vec_mp(tempvec);
	
	randomizer_->receive(source,mpi_config);
	
	return;
}



















int compare_integers_decreasing(const void * left_in, const void * right_in)
{
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left<right) {
		return 1;
	}
	else if (left>right){
		return -1;
	}
	else{
		return 0;
	}
	
}

int compare_integers_increasing(const void * left_in, const void * right_in)
{
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left>right) {
		return 1;
	}
	else if(left < right){
		return -1;
	}
	else{
		return 0;
	}
	
}








namespace color {
	int color_to_int(const std::string c)
	{
		
		if(c.compare("k") == 0){
			return 30;
		}
		else if(c.compare("r") == 0){
			return 31;
		}
		else if(c.compare("g") == 0){
			return 32;
		}
		else if(c.compare("y") == 0){
			return 33;
		}
		else if (c.compare("b") == 0) {
			return 34;
		}
		else if(c.compare("m") == 0){
			return 35;
		}
		else if(c.compare("c") == 0){
			return 36;
		}
		else if(c.compare("l") == 0){
			return 37;
		}
		
		return 30;
		
	}
	
	
	std::string bold(std::string new_color)
	{
		std::stringstream ss;
		ss << "\033[1;" << color_to_int(new_color) << "m" ;     //;
		return ss.str();
	}
	
	std::string dark(std::string new_color)
	{
		std::stringstream ss;
		ss << "\033[2;" << color_to_int(new_color) << "m" ;     //;
		return ss.str();
	}
	
	
	std::string underline(std::string new_color)
	{
		std::stringstream ss;
		ss << "\033[4;" << color_to_int(new_color) << "m" ;     //;
		return ss.str();
	}
	
	
	std::string background(std::string new_color)
	{
		std::stringstream ss;
		ss << "\033[7;" << color_to_int(new_color) << "m" ;     //;
		return ss.str();
	}
	
	
	std::string strike(std::string new_color)
	{
		std::stringstream ss;
		ss << "\033[9;" << color_to_int(new_color) << "m" ;     //;
		return ss.str();
	}
	
	std::string console_default(){
		return "\033[0m";
	}
	
	std::string black(){
		return "\033[0;30m";
	}
	
	std::string red(){
		return "\033[0;31m";
	}
	
	std::string green(){
		return "\033[0;32m";
	}
	
	
	std::string brown(){
		return "\033[0;33m";
	}
	
	std::string blue(){
		return "\033[0;34m";
	}
	
	std::string magenta(){
		return "\033[0;35m";
	}
	
	std::string cyan(){
		return "\033[0;36m";
	}
	
	std::string gray(){
		return "\033[0;37m";
	}
	
	
	//black - 30
	//red - 31
	//green - 32
	//brown - 33
	//blue - 34
	//magenta - 35
	//cyan - 36
	//lightgray - 37
}








