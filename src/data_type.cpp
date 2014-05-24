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
	if (!is_ready()) {
		std::cout << "trying to randomize_d, but is not set up!" << std::endl;
		br_exit(-9509);
	}
	
	
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
		// optimization could be done by looking at the previous function's structure and omitting previously done calculations (for the second and subsequent functions)
		
		
		// copy the current randomizer coefficients into a single row matrix
		for (int jj=0; jj<num_original_funcs; jj++) {
			set_d(&single_row_input_d->entry[0][jj], &randomizer_matrix_d->entry[ii][jj]); // this is repetitive, and wasteful.  optimize away.
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
	if (!is_ready()) {
		std::cout << "trying to randomize_mp, but is not set up!" << std::endl;
		br_exit(-9510);
	}
	
	
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
		std::cout << "requested randomizer which would UP-randomize" << std::endl;
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
	
	int counter = 0;
	int current_degree_index = 0; // start at the end
	int current_degree;
	for (int ii=0; ii<num_desired_rows; ii++) {
		
		counter++;
		if (counter>num_of_each_degree[current_degree_index]) {
			current_degree_index++;
			counter = 1;
		}
		
		current_degree = unique_degrees[current_degree_index];
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
	
	
	max_degree_deficiency = unique_degrees[0] - unique_degrees[num_unique_degrees-1];
	
	
	
	if (num_desired_rows==num_funcs) {
		make_matrix_ID_mp(randomizer_matrix_full_prec,num_funcs,num_funcs);
		square_indicator = true;
	}
	else{
		square_indicator = false;
	}
	
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
	MPI_Send(buffer,5,MPI_INT,target,UNUSED,MPI_COMM_WORLD);
	
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

	
	
	
	MPI_Send(buffer2, size_to_send, MPI_INT, target, UNUSED, MPI_COMM_WORLD);
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
	MPI_Recv(buffer,5,MPI_INT,source,UNUSED,MPI_COMM_WORLD,&statty_mc_gatty);
	
	square_indicator = buffer[0];
	num_randomized_funcs = buffer[1];
	num_original_funcs = buffer[2];
	max_base_degree = buffer[3];
	max_degree_deficiency = buffer[4];
	
	
	delete[] buffer;
	
	// need to send original_degrees, randomized_degrees, randomizer_matrix_full_prec, max_degree_deficiency, and what else?
	int size_to_receive = num_randomized_funcs+num_original_funcs + num_randomized_funcs*num_original_funcs;
	int *buffer2 = new int[size_to_receive];
	
	
	
	
	MPI_Recv(buffer2, size_to_receive, MPI_INT, source, UNUSED, MPI_COMM_WORLD,&statty_mc_gatty);
	
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
	
	if (this->num_points!=0 && this->pts_mp==NULL) {
		printf("trying to add point to point_holder with non-zero num_points and NULL container!\n");
		br_exit(9713);
	}
	
	if (this->num_points==0 && this->pts_mp!=NULL) {
		printf("trying to add point to point_holder with num_points==0 and non-NULL container!\n");
		br_exit(9713);
	}
	
	
	if (this->num_points==0) {
		this->pts_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		this->pts_mp = (vec_mp *)br_realloc(pts_mp, (num_points+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(pts_mp[num_points], new_point->size, new_point->curr_prec);
	this->pts_mp[num_points]->size = new_point->size;
	vec_cp_mp(pts_mp[num_points], new_point);
	
	num_points++;
	
	return num_points-1;
}

int patch_holder::add_patch(vec_mp new_patch)
{
	
	if (this->num_patches!=0 && this->patch_mp==NULL) {
		printf("trying to add patch to witness set with non-zero num_patches and NULL container!\n");
		deliberate_segfault();
	}
	
	if (this->num_patches==0 && this->patch_mp!=NULL) {
		printf("trying to add point to witness set with num_points==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (this->num_patches==0) {
		this->patch_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		this->patch_mp = (vec_mp *)br_realloc(this->patch_mp, (this->num_patches+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(this->patch_mp[this->num_patches], new_patch->size,new_patch->curr_prec);
	this->patch_mp[this->num_patches]->size = new_patch->size;
	vec_cp_mp(this->patch_mp[this->num_patches], new_patch);
	
	this->num_patches++;
	
	
	int burnsome = 0;
	for (int ii=0; ii<this->num_patches; ii++) {
		burnsome += this->patch_mp[ii]->size;
	}
	

	return num_patches-1;
}


int linear_holder::add_linear(vec_mp new_linear)
{
	
	if (this->num_linears!=0 && this->L_mp==NULL) {
		printf("trying to add linear to linear holder with non-zero num_linears and NULL container!\n");
		br_exit(9711);
	}
	
	if (this->num_linears==0 && this->L_mp!=NULL) {
		printf("trying to add linear to linear holder with num_linears==0 and non-NULL container!\n");
		br_exit(9711);
	}
	
	
	if (this->num_linears==0) {
		this->L_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		this->L_mp = (vec_mp *)br_realloc(L_mp, (num_linears+1) * sizeof(vec_mp));
	}
	
	init_vec_mp2(L_mp[num_linears], new_linear->size,new_linear->curr_prec);
	L_mp[this->num_linears]->size = new_linear->size;
	vec_cp_mp(L_mp[num_linears], new_linear);
	
	this->num_linears++;
	
	
	return num_linears-1;
}





//use:  call to parse the file witness_set_file, into the struct W.


int witness_set::witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars)
{
	
	
	
	
	FILE *IN = safe_fopen_read(witness_set_file);
	
	
	int temp_num_patches, patch_size, temp_num_linears, temp_num_points, num_vars_in_linears;
	
	
	fscanf(IN, "%d %d %d", &temp_num_points, &this->dim, &this->comp_num); scanRestOfLine(IN);
	
	
	
	this->num_variables = num_vars;
	this->num_natural_vars = num_vars;
	
	
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
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%d\n\n",this->num_points); // print the header line
	
	for (ii=0; ii<this->num_points; ++ii) {
		for (jj=0; jj<this->pts_mp[ii]->size; jj++) {
			print_mp(OUT,0,&this->pts_mp[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void witness_set::write_dehomogenized_coordinates(boost::filesystem::path filename) const
{
	int ii,jj;
	
	vec_mp result; init_vec_mp(result,1);
	
	
	FILE *OUT = safe_fopen_write(filename.c_str()); // open the output file.
	
	fprintf(OUT,"%d\n\n",this->num_points); // print the header line
	for (ii=0; ii<this->num_points; ++ii) {
		//		print_point_to_screen_matlab(pts_mp[ii],"hompt");
		if (this->num_synth_vars()>0) {
			dehomogenize(&result,this->pts_mp[ii], num_natural_vars);
		}
		else{
			dehomogenize(&result,this->pts_mp[ii]);
		}
		
		for (jj=0; jj<num_natural_vars-1; jj++) {
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
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%d\n\n",this->num_linears); // print the header line
	
	for (ii=0; ii<this->num_linears; ++ii) {
		for (jj=0; jj<this->L_mp[ii]->size; jj++) {
			print_mp(OUT, 0, &this->L_mp[ii]->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}



void witness_set::print_patches(boost::filesystem::path filename) const
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename); // open the output file
	
	
	fprintf(OUT,"%d\n\n",this->num_patches); // print the header line
	
	for (ii=0; ii<this->num_patches; ++ii) {
		fprintf(OUT,"%d\n",this->patch_mp[ii]->size);
		for (jj=0; jj<this->patch_mp[ii]->size; jj++) {
			print_mp(OUT, 0, &this->patch_mp[ii]->coord[jj]);
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
	
	std::cout << "witness set has " << num_variables << " total variables, " << num_natural_vars << " natural variables." << std::endl;
	
	
	std::cout << "dim " << dim << ", comp " << comp_num << std::endl;
	std::cout << "input file name " << input_filename << std::endl;

	int ii;
	printf("******\n%d points\n******\n",this->num_points);
	std::cout << color::green();
	for (ii=0; ii<this->num_points; ii++) {
		
		dehomogenize(&dehom, this->pts_mp[ii], num_natural_vars);
		
		varname << "point_" << ii;
		
		print_point_to_screen_matlab(dehom,varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::blue();
	printf("******\n%d linears\n******\n",this->num_linears);
	
	for (ii=0; ii<this->num_linears; ii++) {
		varname << "linear_" << ii;
		print_point_to_screen_matlab(this->L_mp[ii],varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::cyan();
	printf("******\n%d patches\n******\n",this->num_patches);
	
	for (ii=0; ii<this->num_patches; ii++) {
		varname << "patch_" << ii;
		print_point_to_screen_matlab(this->patch_mp[ii],varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << "variable names:\n";
	for (ii=0; ii< int(this->variable_names.size()); ii++) {
		std::cout << this->variable_names[ii] << "\n";
	}
	printf("\n\n");
	
	
	clear_vec_mp(dehom);
}


void witness_set::print_to_file(boost::filesystem::path filename) const
{
	// print back into the same format we parse from.
	
	
	FILE *OUT = safe_fopen_write(filename);
	
	
	fprintf(OUT, "%d %d %d\n\n", this->num_points, this->dim, this->comp_num);
	


	for (int ii=0; ii < num_points; ii++) {
		
		//read the witness points into memory
		for (int jj=0; jj < pts_mp[ii]->size; ++jj) {
			mpf_out_str(OUT,10,0,pts_mp[ii]->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,pts_mp[ii]->coord[jj].i);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	
	
	
	fprintf(OUT, "%d %d\n", num_linears, num_variables);
	
	
	for (int ii=0; ii < num_linears; ii++) {

		for (int jj=0; jj < L_mp[ii]->size; jj++) {
			mpf_out_str(OUT,10,0,L_mp[ii]->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,L_mp[ii]->coord[jj].i);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	
	
	fprintf(OUT, "%d %d\n", num_patches, num_variables);// TODO: this is incorrect
	
	for (int ii=0; ii < num_patches; ii++) {
		
		
		//read the patch into memory
		for (int jj=0; jj < patch_mp[ii]->size; jj++) {
			mpf_out_str(OUT,10,0,patch_mp[ii]->coord[jj].r); // 10 is the base
			fprintf(OUT," ");
			mpf_out_str(OUT,10,0,patch_mp[ii]->coord[jj].i);
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
	witness_set::only_first_vars(num_natural_vars);
}


void witness_set::only_first_vars(int num_vars)
{
	
	vec_mp tempvec;  init_vec_mp2(tempvec, num_vars, 1024);
	tempvec->size = num_vars;
	
	for (int ii=0; ii<this->num_points; ii++) {
		
		for (int jj=0; jj<num_vars; jj++) {
			set_mp(&tempvec->coord[jj], &this->pts_mp[ii]->coord[jj]);
		}
		
		change_size_vec_mp(this->pts_mp[ii], num_vars);  this->pts_mp[ii]->size = num_vars;
		vec_cp_mp(this->pts_mp[ii], tempvec);
	}
	
	
	this->num_variables = num_vars;
	
	int patch_size_counter = 0, trim_from_here =  0;
	for (int ii=0; ii<this->num_patches; ii++) {
		patch_size_counter += this->patch_mp[ii]->size;
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
	
	for (int ii=trim_from_here; ii<this->num_patches; ii++) {
		clear_vec_mp(this->patch_mp[ii]);
	}
	
	this->patch_mp = (vec_mp *) br_realloc(this->patch_mp, trim_from_here* sizeof(vec_mp));
	this->num_patches = trim_from_here;
	
	for (int ii=0; ii<num_linears; ii++) {
		L_mp[ii]->size = num_vars;
	}
	
	clear_vec_mp(tempvec);
	return;
}



void witness_set::sort_for_real(tracker_config_t * T)
{
	
	
	int ii;
	
	
	int *real_indicator = new int[this->num_points];
	int counter = 0;
	
	vec_mp result; init_vec_mp(result,num_natural_vars-1);
	result->size = num_natural_vars-1;
	
	for (ii=0; ii<this->num_points; ii++) {
		for (int jj=1; jj<num_natural_vars; jj++) {
			div_mp(&result->coord[jj-1], &this->pts_mp[ii]->coord[jj], &this->pts_mp[ii]->coord[0]);
		}
		real_indicator[ii] = checkForReal_mp(result, T->real_threshold);
		
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	vec_mp *tempvec = (vec_mp *)br_malloc(counter * sizeof(vec_mp));
	
	counter = 0;  // reset
	for (ii=0; ii<this->num_points; ii++) {
		if (real_indicator[ii]==1) {
			
			init_vec_mp2(tempvec[counter],this->num_variables,1024); tempvec[counter]->size = this->num_variables;
			vec_cp_mp(tempvec[counter],this->pts_mp[ii]);
			counter++;
		}
		else{
			
		}
	}
	
	clear_vec_mp(result);
	
	
	for (ii=0; ii<this->num_points; ii++) {
		clear_vec_mp(this->pts_mp[ii]);
	}
	free(this->pts_mp);
	
	this->pts_mp = tempvec;
	this->num_points = counter;
	
	delete[] real_indicator;
	return;
}






// T is necessary for the tolerances.
void witness_set::sort_for_unique(tracker_config_t * T)
{
	
	if (num_natural_vars==0) {
		std::cout << "blabla" << std::endl;
	}
	
	int curr_uniqueness;
	int num_good_pts = 0;
	std::vector<int> is_unique;
	
	for (int ii = 0; ii<this->num_points; ++ii) {
		curr_uniqueness = 1;
		
		int prev_size_1 = pts_mp[ii]->size;  pts_mp[ii]->size = num_natural_vars; // cache and change to natural number
		
		for (int jj=ii+1; jj<this->num_points; ++jj) {
			int prev_size_2 = pts_mp[jj]->size; pts_mp[jj]->size = num_natural_vars; // cache and change to natural number
			if ( isSamePoint_homogeneous_input(this->pts_mp[ii],this->pts_mp[jj]) ){
				curr_uniqueness = 0;
			}
			pts_mp[jj]->size = prev_size_2; // restore
		}
		pts_mp[ii]->size = prev_size_1; // restore
		
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
	for (int ii=0; ii<this->num_points; ++ii) {
		if (is_unique[ii]==1) {
			init_vec_mp2(transferme[counter],this->num_variables,1024);  transferme[counter]->size = this->num_variables;
			vec_cp_mp(transferme[counter], this->pts_mp[ii]);
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		printf("counter mismatch\n");
		exit(270);
	}
	
	for (int ii=0; ii<this->num_points; ii++) {
		clear_vec_mp(this->pts_mp[ii]);
	}
	free(this->pts_mp);
	
	
	this->num_points = num_good_pts;
	this->pts_mp = transferme;
	
	return;
}





void witness_set::sort_for_inside_sphere(comp_mp radius, vec_mp center)
{
	
	
	
	
	int num_good_pts = 0;
	
	
	std::vector<int> is_ok;
	
	vec_mp temp_vec; init_vec_mp(temp_vec,0);
	comp_mp temp; init_mp(temp);
	
	for (int ii = 0; ii<this->num_points; ++ii) {
		
		
		
		dehomogenize(&temp_vec, this->pts_mp[ii],num_natural_vars);
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
	for (int ii=0; ii<this->num_points; ++ii) {
		if (is_ok[ii]==1) {
			init_vec_mp2(transferme[counter],0,1024);  transferme[counter]->size = 0;
			vec_cp_mp(transferme[counter], this->pts_mp[ii]);
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		printf("counter mismatch\n");
		br_exit(271);
	}
	
	for (int ii=0; ii<this->num_points; ii++) {
		clear_vec_mp(this->pts_mp[ii]);
	}
	free(this->pts_mp);
	
	
	this->num_points = num_good_pts;
	this->pts_mp = transferme;
	
	clear_vec_mp(temp_vec);
	clear_mp(temp);
	return;
}




void witness_set::merge(const witness_set & W_in)
{
	
	//error checking first
	if ( (this->num_variables==0) && (W_in.num_variables!=0) && (this->num_points==0)) {
		this->num_variables = W_in.num_variables;
		this->num_natural_vars = W_in.num_natural_vars;
	}
//	else if ( (this->num_variables!=0) && (W_in.num_variables==0) ){
//		
//	}
	else if ( (W_in.num_variables!=this->num_variables) ) {//&& (W_in.num_points!=0) && (this->num_points!=0)
		printf("merging two witness sets with differing numbers of variables.\n");
		std::cout << "existing: " << this->num_variables << ", input: " << W_in.num_variables << std::endl;
//		deliberate_segfault();
	}
	
	
	if (W_in.num_natural_vars != this->num_natural_vars) {
		printf("merging two witness sets with differing numbers of natural variables. %d merging set, %d existing\n",
			   W_in.num_natural_vars, this->num_natural_vars);
		deliberate_segfault();
//		br_exit(95);
	}
	
	//just mindlessly add the linears.  up to user to ensure linears get merged correctly.  no way to know what they want...
	for (int ii = 0; ii<W_in.num_linears; ii++) {
		witness_set::add_linear(W_in.L_mp[ii]);
	}
	
	
	for (int ii = 0; ii<W_in.num_points; ii++) {
		int is_new = 1;
		for (int jj = 0; jj<this->num_points; jj++){
			if (this->pts_mp[jj]->size == W_in.pts_mp[ii]->size) {
				if (isSamePoint_inhomogeneous_input(this->pts_mp[jj], W_in.pts_mp[ii])) {
					is_new = 0;
					break;
				}
			}
			
		}
		
		if (is_new==1)
			witness_set::add_point(W_in.pts_mp[ii]);
	}
	
	
	for (int ii = 0; ii<W_in.num_patches; ii++) {
		int is_new = 1;
		
		for (int jj = 0; jj<this->num_patches; jj++){
			if (this->patch_mp[jj]->size == W_in.patch_mp[ii]->size) {
				if (isSamePoint_inhomogeneous_input(this->patch_mp[jj], W_in.patch_mp[ii])) {
					is_new = 0;
					break;
				}
			}
		}
		
		if (is_new==1)
			witness_set::add_patch(W_in.patch_mp[ii]);
	}
	
	
	
	return;
}//re: merge_witness_sets
















void witness_set::send(parallelism_config & mpi_config, int target)
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
    buffer[0] = dim;
    buffer[1] = comp_num;
    buffer[2] = incidence_number;
    buffer[3] = num_variables;
    buffer[4] = num_natural_vars;
    buffer[5] = num_points;
    buffer[6] = num_linears;
    buffer[7] = num_patches;
    
    MPI_Send(buffer, 8, MPI_INT, WITNESS_SET, target, mpi_config.my_communicator);
    
    free(buffer);
    
    for (int ii=0; ii<num_linears; ii++) {
        send_vec_mp(L_mp[ii],target);
    }
    for (int ii=0; ii<num_patches; ii++) {
        send_vec_mp(patch_mp[ii],target);
    }
    for (int ii=0; ii<num_points; ii++) {
        send_vec_mp(pts_mp[ii],target);
    }
    
    char * namebuffer = (char *) br_malloc(1024*sizeof(char));
    
	
    free(namebuffer);
    return;
}

void witness_set::receive(int source, parallelism_config & mpi_config)
{
    MPI_Status statty_mc_gatty;
    
    int *buffer = (int *) br_malloc(8*sizeof(int));
    
    
    MPI_Recv(buffer, 8, MPI_INT, WITNESS_SET, source, mpi_config.my_communicator, &statty_mc_gatty);
    
    dim = buffer[0];
    comp_num = buffer[1];
    incidence_number = buffer[2];
    num_variables = buffer[3];
    num_natural_vars = buffer[4];
    num_points = buffer[5];
    num_linears = buffer[6];
    num_patches = buffer[7];
    
    free(buffer);
    
    vec_mp tempvec; init_vec_mp2(tempvec,0,1024);
    
    for (int ii=0; ii<num_linears; ii++) {
        receive_vec_mp(tempvec,source);
        add_linear(tempvec);
    }
    
    for (int ii=0; ii<num_patches; ii++) {
        receive_vec_mp(tempvec,source);
        add_patch(tempvec);
    }
    
    for (int ii=0; ii<num_points; ii++) {
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

	
	fscanf(IN,"%d %d", &num_variables, &num_nonempty_codims);
	
	
	
	vec_mp prev_approx;
	init_vec_mp(prev_approx,num_variables); prev_approx->size = num_variables;
	
	vec_mp temp_vec;
	init_vec_mp(temp_vec,num_variables); temp_vec->size = num_variables;
	
	
	//create counters and initialize to zeros
	std::vector<int> component_counter(num_nonempty_codims,0),codim_indicator(num_nonempty_codims,0);
	
	
	for (int ii=0; ii<num_nonempty_codims; ii++) {
		int current_codimension,num_points_this_dim;
		fscanf(IN,"%d %d",&current_codimension,&num_points_this_dim);
		
		codim_indicator[ii] = current_codimension;
		
		int current_dimension = num_variables-1-current_codimension;
		
		nonempty_dimensions.push_back(current_dimension);
		
		for (int jj=0; jj<num_points_this_dim; jj++) {
			
			
			int precision;
			// last approximation
			fscanf(IN,"%d",&precision);
			change_prec_vec_mp(temp_vec,precision);
			for (int kk=0; kk<num_variables; kk++) {
				mpf_inp_str(temp_vec->coord[kk].r, IN, 10); // 10 is the base
				mpf_inp_str(temp_vec->coord[kk].i, IN, 10);
			}
			
			
			
			// previous approximation
			fscanf(IN,"%d",&precision);
			change_prec_vec_mp(temp_vec,precision);
			for (int kk=0; kk<num_variables; kk++) {
				mpf_inp_str(prev_approx->coord[kk].r, IN, 10); // 10 is the base
				mpf_inp_str(prev_approx->coord[kk].i, IN, 10);
			}
			
			
			witness_point_metadata meta(current_dimension);
			
			meta.set_from_file(IN);

			int index = add_solution(temp_vec, meta);
			
			int num_existing_pts_this_comp = map_lookup_with_default( dimension_component_counter[current_dimension], meta.component_number, 0 );
			
			dimension_component_counter[current_dimension][meta.component_number] = num_existing_pts_this_comp+1;
			index_tracker[current_dimension][meta.component_number].push_back(index);
			
			
			if (meta.component_number+1 > component_counter[ii]) {
				component_counter[ii]  = meta.component_number+1;  // keeps track of the number of components per dimension
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
	
	
	mpq_t *temp_rat = new mpq_t[2];  init_rat(temp_rat);
	
	
	for (int mm=0; mm<num_nonempty_codims; mm++) {
		//BEGIN REPEATING BLOCK, one per nonempty codimension
		int current_codimension;
		current_codimension = codim_indicator[mm];
		int current_dimension = num_variables-1-current_codimension;
		
		int num_rows_randomization, num_cols_randomization;
		fscanf(IN,"%d %d",&num_rows_randomization,&num_cols_randomization);
		
		mat_mp randomization_matrix;  init_mat_mp(randomization_matrix,num_rows_randomization,num_cols_randomization);
		for (int ii = 0; ii < num_rows_randomization; ii++) {
			for (int jj=0; jj<num_cols_randomization; jj++) {
				if (numbertype==2) {
					setup_comp_in_rat(temp_rat,IN);
					rat_to_mp(&randomization_matrix->entry[ii][jj],temp_rat);
				}
				else{
					mpf_inp_str(randomization_matrix->entry[ii][jj].r, IN, 10);
					mpf_inp_str(randomization_matrix->entry[ii][jj].i, IN, 10);
				}
			}
		}
		
		
		
		//MATRIX W FOR HOMOGENIZATION
		//  same length as the randomization matrix
		int hom_mat[num_rows_randomization][num_cols_randomization];
		for (int ii=0; ii<num_rows_randomization; ii++)
			for (int jj=0; jj<num_cols_randomization; jj++)
				fscanf(IN,"%d",&hom_mat[ii][jj]);
		
		
		
		int num_entries_hom_patch_eqn;
		fscanf(IN,"%d",&num_entries_hom_patch_eqn);
		vec_mp homogenization_patch_eqn;  init_vec_mp(homogenization_patch_eqn,num_entries_hom_patch_eqn); homogenization_patch_eqn->size = num_entries_hom_patch_eqn;
		//		//VECTOR H FOR HOMOGENIZATION
		for (int ii = 0; ii<num_entries_hom_patch_eqn; ii++) {
			if (numbertype==2) {
				setup_comp_in_rat(temp_rat,IN);
				rat_to_mp(&temp_vec->coord[ii],temp_rat);
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
		rat_to_mp(hom_variable_constant,temp_rat);
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
					rat_to_mp(&temp_vec->coord[jj],temp_rat);
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
				rat_to_mp(&temp_vec->coord[ii],temp_rat);
			}
			else{
				mpf_inp_str(temp_vec->coord[ii].r, IN, 10);
				mpf_inp_str(temp_vec->coord[ii].i, IN, 10);
			}
			
		}
		
		add_patch_w_meta(temp_vec, witness_patch_metadata(current_dimension));
		
		//END REPEATED BLOCK (ONE FOR EACH NONEMPTY CODIM).
	}
	clear_rat(temp_rat);
	delete [] temp_rat;
	
	
	fclose(IN);
	
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
					witness_set W(num_variables);
					W.dim = nonempty_dimensions[0];
					W.comp_num = -1;
					deliberate_segfault();
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
			
			witness_set W(num_variables);
			W.dim = options.target_dimension;
			W.comp_num = -1;
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
					witness_set W(num_variables);
					W.dim = target_dimension;
					W.comp_num = -1;
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
	
	witness_set W(num_variables); // create blank witness set
	W.dim = target_dimension;
	W.comp_num = -1;
	
	
	if (index_tracker.find(target_dimension)==index_tracker.end()) {
		std::cout << "must return empty set from auto constructor.  have no components of dimension " << target_dimension << std::endl;
		return W;
	}
	
	
	
	std::vector<int> components_with_no_deflations_needed;
	for (auto iter=index_tracker[target_dimension].begin(); iter!=index_tracker[target_dimension].end(); ++iter) {
		// iterate over components for target dimension
		if (iter->second.size()==0) {
			std::cout << "what the hell?" << std::endl;
			mypause();
		}
		else{
			if (point_metadata[iter->second[0]].deflations_needed==0) {
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
		
		if (checkSelfConjugate(pts_mp[current_index], options, options.input_filename)==true) {
			std::cout << "dim " << target_dimension << ", comp " << *iter << " is self-conjugate" << std::endl;
			for (int ii=0; ii<dimension_component_counter[target_dimension][*iter]; ++ii) {
				W.add_point(pts_mp[index_tracker[target_dimension][*iter][ii]]);
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
			if (linear_metadata[ii].dimension == target_dimension) {
				W.add_linear(L_mp[ii]);
			}
		}
		
		for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
			if (patch_metadata[ii].dimension == target_dimension) {
				W.add_patch(patch_mp[ii]);
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
	
	
		
	if (target_dimension==-1) { // all dimensions
		
		for (auto iter=index_tracker.begin(); iter!=index_tracker.end(); ++iter) {
			std::cout << "dimension " << iter->first << std::endl;
			for (auto jter = iter->second.begin(); jter!=iter->second.end(); ++jter) {
				
				
				int first_index = *(jter->second.begin());
				std::cout << "\tcomponent " << jter->first << ": multiplicity " << point_metadata[first_index].multiplicity << ", deflations needed: " << point_metadata[first_index].deflations_needed << std::endl;
			}
			std::cout << std::endl;
		}
		
		std::set<int> valid_choices;
		for (auto iter=nonempty_dimensions.begin(); iter!=nonempty_dimensions.end(); ++iter) {
			valid_choices.insert(*iter);
		}
		
		target_dimension = get_int_choice("\n\nchoose a dimension:\n",valid_choices);
		
	}
		
	
	std::cout << "dimension " << options.target_dimension << std::endl;
	for (auto jter = index_tracker[target_dimension].begin(); jter!= index_tracker[target_dimension].end(); ++jter) {
		int first_index = *(jter->second.begin());
		std::cout << "\tcomponent " << jter->first << ": multiplicity " << point_metadata[first_index].multiplicity << ", deflations needed: " << point_metadata[first_index].deflations_needed << ", degree " << dimension_component_counter[target_dimension][jter->first] << std::endl;
		
	}
	std::cout << std::endl;
	
	
	
	int num_components = get_int_choice("how many components would you like to decompose? (0 for all mult-one components)\n",0,dimension_component_counter[target_dimension].size());
	
	
	std::set<int> chosen_few;
	
	if (num_components==0) {
		for (unsigned int ii=0; ii<dimension_component_counter[target_dimension].size(); ++ii) {
			if (point_metadata[index_tracker[target_dimension][ii][0]].multiplicity==1) {
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
	witness_set W(num_variables);
	
	W.dim = target_dimension;
	if (chosen_few.size()==1) {
		W.comp_num = *chosen_few.begin();
	}
	else{
		W.comp_num = -1;
	}
	
	for (auto iter=chosen_few.begin(); iter!=chosen_few.end(); ++iter) {
		
		//iterate over each point in the component
		for (auto jter=index_tracker[target_dimension][*iter].begin(); jter!=index_tracker[target_dimension][*iter].end(); ++jter) {
			std::cout << "adding point " << *jter << std::endl;
			print_point_to_screen_matlab(pts_mp[*jter],"p");
			W.add_point(pts_mp[*jter]);
		}
		
	}
	
	for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
		if (linear_metadata[ii].dimension == target_dimension) {
			W.add_linear(L_mp[ii]);
		}
	}
	
	for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
		if (patch_metadata[ii].dimension == target_dimension) {
			W.add_patch(patch_mp[ii]);
		}
	}
	
	return W;
}





witness_set witness_data::form_specific_witness_set(int dim, int comp)
{
	witness_set W(num_variables);
	
	W.dim = dim;
	W.comp_num = comp;
	
	
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
		W.add_point(pts_mp[*iter]);
	}
	
	
	for (unsigned int ii=0; ii<linear_metadata.size(); ii++) {
		if (linear_metadata[ii].dimension == dim) {
			W.add_linear(L_mp[ii]);
		}
	}
	
	for (unsigned int ii=0; ii<patch_metadata.size(); ii++) {
		if (patch_metadata[ii].dimension == dim) {
			W.add_patch(patch_mp[ii]);
		}
	}
	
	
	return W;
}


void vertex::set_point(const vec_mp new_point)
{
	change_prec_vec_mp(this->pt_mp, new_point->curr_prec);
	vec_cp_mp(this->pt_mp, new_point);
}




void vertex::send(int target, parallelism_config & mpi_config)
{
	
	send_vec_mp(pt_mp, target);
	
	send_vec_mp(projection_values, target);
	
	int * buffer = (int *) br_malloc(3*sizeof(int));
	buffer[0] = type;
	buffer[1] = removed;
	buffer[2] = input_filename_index;
	
	MPI_Send(buffer, 3, MPI_INT, target, VERTEX, MPI_COMM_WORLD);
	free(buffer);
	
}


void vertex::receive(int source, parallelism_config & mpi_config)
{
	MPI_Status statty_mc_gatty;
	int * buffer = (int *) br_malloc(3*sizeof(int));
	
	
	receive_vec_mp(pt_mp, source);
	receive_vec_mp(projection_values, source);
	
	MPI_Recv(buffer, 3, MPI_INT, source, VERTEX, MPI_COMM_WORLD, &statty_mc_gatty);
	
	type = buffer[0];
	removed = buffer[1];
	input_filename_index = buffer[2];
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
    for (int jj=1; jj<num_natural_variables; jj++) {
		div_mp(&checker_1->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
	}
    
    
	//	WTB: a faster comparison search.
	for (int ii=0; ii<num_vertices; ii++) {
		
		int current_index = ii;
		
		if (vertices[current_index].removed==0) {
			
			// dehomogenize the current point under investigation
			for (int jj=1; jj<num_natural_variables; jj++) {
				div_mp(&checker_2->coord[jj-1], &vertices[current_index].pt_mp->coord[jj], &vertices[current_index].pt_mp->coord[0]);
			}
			
			if (isSamePoint_inhomogeneous_input(checker_1, checker_2)){
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
    
    for (int jj=1; jj<num_natural_variables; jj++) {
		div_mp(&checker_1->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
	}
    
    
    //	WTB: a faster comparison search.
    for (int ii=0; ii<num_vertices; ii++) {
		int current_index = ii;
		
		if (vertices[current_index].removed==1) {
			
			for (int jj=1; jj<num_natural_variables; jj++) {
				div_mp(&checker_2->coord[jj-1], &vertices[current_index].pt_mp->coord[jj], &vertices[current_index].pt_mp->coord[0]);
			}
			
			if (isSamePoint_inhomogeneous_input(checker_1, checker_2)){
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
    
	vec_mp projection_values; init_vec_mp2(projection_values,W.num_points,1024);
	projection_values->size = W.num_points;
	
	for (int ii=0; ii<W.num_points; ii++){
		
		for (int jj=0; jj<W.pts_mp[ii]->size; jj++) {
			
			if (!(mpfr_number_p(W.pts_mp[ii]->coord[jj].r) && mpfr_number_p(W.pts_mp[ii]->coord[jj].i))) {
				std::cout << color::red();
				std::cout << "there was NAN in a coordinate for a point in the projections to sort :(" << std::endl;
				print_point_to_screen_matlab(W.pts_mp[ii], "bad_point");
				
				print_point_to_screen_matlab(pi,"pi");
				std::cout << color::console_default();
				return CRITICAL_FAILURE;
			}
			
		}
		
        
        int curr_index = search_for_point(W.pts_mp[ii]);
        
        if (curr_index < 0) {
            std::cout << color::red() << "trying to retrieve projection value from a non-stored point" << color::console_default() << std::endl;
            mypause();
        }
        

        
        set_mp(&projection_values->coord[ii], &vertices[curr_index].projection_values->coord[proj_index]);
        
		
		
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
    if (this->curr_projection<0) {
        std::cout << color::red() << "trying to assert projection value (current index) without having index set" << color::console_default() << std::endl;
        br_exit(-91621);
    }
    return assert_projection_value(relevant_indices,new_value,this->curr_projection);
}

std::vector<int> vertex_set::assert_projection_value(const std::set< int > & relevant_indices, comp_mp new_value, int proj_index)
{
	std::vector<int> bad_indices;
    
    comp_mp temp; init_mp(temp);
    
    if ( proj_index > num_projections ) {
        std::cout << color::red() << "trying to assert projection value, but index of projection is larger than possible" << color::console_default() << std::endl;
		br_exit(3091); // throw?
    }
    
    for (std::set<int>::iterator ii=relevant_indices.begin(); ii!=relevant_indices.end(); ii++) {
        //*ii
        
        sub_mp(temp, &vertices[*ii].projection_values->coord[proj_index], new_value);
        if (fabs(mpf_get_d(temp->r))>0.0001) {
            std::cout << "trying to assert projection value of " << mpf_get_d(new_value->r)
				      << " but original value is " << mpf_get_d(vertices[*ii].projection_values->coord[proj_index].r) << std::endl;
            std::cout << "point index is " << *ii << std::endl;
			bad_indices.push_back(*ii);
			continue;
        }
        
        set_mp(&vertices[*ii].projection_values->coord[proj_index], new_value);
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
	
	
	this->vertices.push_back(source_vertex);
	
	
	if (this->vertices[num_vertices].projection_values->size < this->num_projections) {
		increase_size_vec_mp(this->vertices[num_vertices].projection_values, this->num_projections );
		this->vertices[num_vertices].projection_values->size = this->num_projections;
	}
	
	for (int ii=0; ii<this->num_projections; ii++){
		bool compute_proj_val = false;
		
		if (! (mpfr_number_p( vertices[num_vertices].projection_values->coord[ii].r) && mpfr_number_p( vertices[num_vertices].projection_values->coord[ii].i)  ) )
		{
			compute_proj_val = true;
		}
		//		else if ( false )//yeah, i dunno what else right yet.
		//		{
		//			compute_proj_val = true;
		//		}
		
		
		if (compute_proj_val==true) {
			projection_value_homogeneous_input(&vertices[num_vertices].projection_values->coord[ii],
											   vertices[num_vertices].pt_mp,
											   projections[ii]);
#ifdef thresholding
            real_threshold(&vertices[num_vertices].projection_values->coord[ii],1e-13);
#endif
		}
		
	}
	
    
	
	if (vertices[num_vertices].input_filename_index == -1)
	{
		vertices[num_vertices].input_filename_index = curr_input_index;
	}
		
	
	this->num_vertices++;
	return this->num_vertices-1;
}


void vertex_set::print_to_screen()
{
	printf("vertex set has %d vertices:\n\n",this->num_vertices);
	for (int ii=0; ii<this->num_vertices; ++ii) {
		print_point_to_screen_matlab(this->vertices[ii].pt_mp,"vert");
		print_point_to_screen_matlab(this->vertices[ii].projection_values,"projection_values");
		printf("type: %d\n", this->vertices[ii].type);
	}
}


int vertex_set::setup_vertices(boost::filesystem::path INfile)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices;
	int num_vars;
	int tmp_num_projections;
	int tmp_num_filenames;
	fscanf(IN, "%d %d %d %d\n\n", &num_vertices, &tmp_num_projections, &this->num_natural_variables, &tmp_num_filenames);
	
	
	vec_mp temp_vec; init_vec_mp2(temp_vec,num_natural_variables,1024);
	for (int ii=0; ii<tmp_num_projections; ii++) {
		for (int jj=0; jj<num_natural_variables; jj++) {
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
		this->filenames.push_back(temppath);
		
		delete [] buffer;
	}
	
	
	
	vertex temp_vertex;
	
	for (int ii=0; ii<num_vertices; ii++)
	{
		fscanf(IN, "%d\n", &num_vars);
		if (temp_vertex.pt_mp->size != num_vars) {
			change_size_vec_mp(temp_vertex.pt_mp,num_vars); temp_vertex.pt_mp->size = num_vars;
		}
		
		for (int jj=0; jj<num_vars; jj++)
		{
			mpf_inp_str(temp_vertex.pt_mp->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vertex.pt_mp->coord[jj].i, IN, 10);
		}
		
		int temp_num;
		fscanf(IN,"%d\n",&temp_num);
		increase_size_vec_mp(temp_vertex.projection_values,temp_num);
		temp_vertex.projection_values->size = temp_num;
		for (int jj=0; jj<temp_num; jj++) {
			mpf_inp_str(temp_vertex.projection_values->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vertex.projection_values->coord[jj].i, IN, 10);
		}
		
//		print_point_to_screen_matlab(temp_vertex.projection_values,"p");
		fscanf(IN,"%d\n",&temp_vertex.input_filename_index);
		fscanf(IN,"%d\n",&temp_vertex.type);
		
		vertex_set::add_vertex(temp_vertex);
	}
	
	
	
	fclose(IN);
	
	if (this->num_vertices!=num_vertices) {
		printf("parity error in num_vertices.\n\texpected: %d\tactual: %d\n",this->num_vertices,num_vertices); // this is totally impossible.
		br_exit(25943);
	}
	
	return num_vertices;
}







void vertex_set::print(boost::filesystem::path outputfile) const
{
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d %d %d %lu\n\n",num_vertices,num_projections, num_natural_variables, filenames.size());
	
	
	
	for (int ii=0; ii<num_projections; ii++) {
		for (int jj=0; jj<num_natural_variables; jj++) {
			print_mp(OUT, 0, &projections[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	
	
	for (unsigned int ii=0; ii!=filenames.size(); ii++) {
		int strleng = filenames[ii].string().size() + 1; // +1 for the null character
		char * buffer = new char[strleng];
		memcpy(buffer, filenames[ii].c_str(), strleng);
		fprintf(OUT,"%d\n",strleng);
		fprintf(OUT,"%s\n",buffer);
		delete [] buffer;
	}
	
	for (int ii = 0; ii < num_vertices; ii++)
	{ // output points
		fprintf(OUT,"%d\n", vertices[ii].pt_mp->size);
		for(int jj=0;jj<vertices[ii].pt_mp->size;jj++) {
			print_mp(OUT, 0, &vertices[ii].pt_mp->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"%d\n",vertices[ii].projection_values->size);
		for(int jj=0;jj<vertices[ii].projection_values->size;jj++) {
			print_mp(OUT, 0, &vertices[ii].projection_values->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"%d\n",vertices[ii].input_filename_index);
		
		fprintf(OUT,"\n");
		if (vertices[ii].removed==1) {
			fprintf(OUT,"%d\n\n",REMOVED);
		}
		else{
			fprintf(OUT,"%d\n\n",vertices[ii].type);
		}
	}
	
	
	
	fclose(OUT);
	
}





void vertex_set::send(int target, parallelism_config & mpi_config)
{
	

	int num_filenames = filenames.size();
	
	
	int * buffer2 = new int[6];
	buffer2[0] = num_natural_variables;
	buffer2[1] = num_projections;
	buffer2[2] = curr_projection;
	buffer2[3] = num_filenames;
	buffer2[4] = curr_input_index;
	buffer2[5] = num_vertices;
	
	MPI_Send(buffer2, 6, MPI_INT, target, VERTEX_SET, MPI_COMM_WORLD);
	
	delete [] buffer2;
//	std::cout << "sending " << num_projections << " projections" << std::endl;
	
	buffer2 = new int[num_projections];
	for (int ii=0; ii<num_projections; ii++) {
		buffer2[ii] = projections[ii]->size;
	}
	MPI_Send(buffer2, num_projections, MPI_INT, target, VERTEX_SET, MPI_COMM_WORLD);
	delete [] buffer2;
	
	for (int ii=0; ii<num_projections; ii++) {
		send_vec_mp(projections[ii],target);
	}

//	std::cout << "sending " << num_filenames << " filenames" << std::endl;
	for (int ii=0; ii<num_filenames; ii++) {
		char * buffer;
		
		
		int strleng = filenames[ii].string().size()+1;
		buffer = new char[strleng];
		memcpy(buffer, filenames[ii].string().c_str(), strleng-1); // this sucks
		buffer[strleng-1] = '\0';
		
		MPI_Send(&strleng, 1, MPI_INT, target, VERTEX_SET, MPI_COMM_WORLD);
//		std::cout << "sending filename length " << strleng << " " << filenames[ii].string() << std::endl;
		MPI_Send(&buffer[0], strleng, MPI_CHAR, target, VERTEX_SET, MPI_COMM_WORLD);
		
		delete [] buffer;
		
	}
	
	
	for (int ii=0; ii<num_vertices; ii++) {
		vertices[ii].send(target, mpi_config);
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
	curr_projection = buffer2[2];
	int temp_num_filenames = buffer2[3];
	curr_input_index = buffer2[4];
	int temp_num_vertices = buffer2[5];
	
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
	
	if (num_projections!=temp_num_projections) {
		std::cout << "num_projections doesn't match!" << std::endl;
	}
	
//	std::cout << "receiving " << temp_num_filenames << " filenames" << std::endl;
	for (int ii=0; ii<temp_num_filenames; ii++) {
		char * buffer; int strleng;
		
		MPI_Recv(&strleng, 1, MPI_INT, source, VERTEX_SET, MPI_COMM_WORLD, &statty_mc_gatty);
		
		buffer = new char[strleng];
//		std::cout << "recving filename length " << strleng << std::endl;
		MPI_Recv(&buffer[0], strleng, MPI_CHAR, source, VERTEX_SET, MPI_COMM_WORLD, &statty_mc_gatty);
		filenames.push_back(boost::filesystem::path(std::string(buffer)));
		
		delete [] buffer;
		
	}
	
	
	
	

	
	for (int ii=0; ii<temp_num_vertices; ii++) {
		vertex tempvert;
		tempvert.receive(source, mpi_config);
		add_vertex(tempvert);
	}
	
	if (num_vertices != temp_num_vertices) {
		std::cout << "logical inconsistency.  do not have correct num vertices." << std::endl;
	}
	
	
	

	
	return;
}










































int decomposition::add_witness_set(const witness_set & W, int add_type, vertex_set & V)
{
#ifdef functionentry_output
	std::cout << "decomposition::add_witness_set" << std::endl;
#endif
	
	
    V.set_curr_input(W.input_filename);
    
    vertex temp_vertex;
    temp_vertex.type = add_type;
    
    for (int ii=0; ii<W.num_points; ii++) {
        vec_cp_mp(temp_vertex.pt_mp, W.pts_mp[ii]);
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
	
	if (this->counters.find(source_vertex.type) == this->counters.end()) {
		this->counters[source_vertex.type] = 0;
	}
	
	this->counters[source_vertex.type]++;
	this->indices[source_vertex.type].push_back(current_index);
	
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
	int index = decomposition::index_in_vertices(V, vert.pt_mp);
	
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
	input_filename = directoryName / tempstr;
	
	getline(fin, tempstr);
	converter << tempstr;
	converter >> this->num_variables >> this->dimension;
	converter.clear(); converter.str("");
	

	int num_types;
	fin >> num_types;
	
	
	for (int ii =0; ii<num_types; ii++) {
		int current_type, num_this_type, current_index;
		fin >> current_type >> num_this_type;

		this->counters[current_type] = num_this_type;
		
		for (int jj=0; jj<num_this_type; jj++) {
			fin >> current_index;
			this->indices[current_type].push_back(current_index);
		}
	}
	
	vec_mp tempvec; init_vec_mp2(tempvec, this->num_variables,1024);
	tempvec->size = this->num_variables;
	

	for (int ii=0; ii<dimension; ii++) {
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
	
	
	
	clear_vec_mp(temp_patch);
	
	fin.close();
	
	
	
	return 0;
}






void decomposition::print(boost::filesystem::path base)
{
	
#ifdef functionentry_output
	std::cout << "decomposition::print" << std::endl;
#endif
	
	if (dimension != num_curr_projections) {
//		std::cout << "decomposition was short projections\nneeded	" << this->dimension << " but had " << num_curr_projections << std::endl;;
	}
	
	
	int ii;
	
	boost::filesystem::create_directory(base);
	
	FILE *OUT = safe_fopen_write(base / "decomp");
	
	fprintf(OUT,"%s\n",input_filename.filename().c_str());
	
	fprintf(OUT,"%d %d\n\n",num_variables, dimension);
	
	fprintf(OUT, "%d\n", int(this->counters.size()));
	
	std::map< int, int>::iterator type_iter;
	for (type_iter = this->counters.begin(); type_iter!= this->counters.end(); type_iter++) {
		fprintf(OUT, "%d %d\n",type_iter->first, type_iter->second);  // print the number corresponding to the type, and the number of that type.
		
		for (int jj = 0; jj<type_iter->second; jj++) {
			fprintf(OUT, "%d\n", indices[type_iter->first][jj]);
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	
	for (ii=0; ii<num_curr_projections; ii++) {
		fprintf(OUT,"%d\n",pi[ii]->size);
		for(int jj=0;jj<pi[ii]->size;jj++)
		{
			print_mp(OUT, 0, &pi[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	fprintf(OUT,"%d\n\n",this->num_patches); // print the header line
	
	for (ii=0; ii<this->num_patches; ++ii) {
		fprintf(OUT,"%d\n",this->patch_mp[ii]->size);
		for (int jj=0; jj<this->patch_mp[ii]->size; jj++) {
			print_mp(OUT, 0, &this->patch_mp[ii]->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
    
    fprintf(OUT,"\n\n");
    
    print_mp(OUT, 0, this->sphere_radius);
    fprintf(OUT, "\n%d\n",this->sphere_center->size);
    
    for (int jj=0; jj<this->sphere_center->size; jj++) {
        print_mp(OUT, 0, &this->sphere_center->coord[jj]);
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
	
	if (num_variables<2) {
		std::cout << "during read of sphere, decomposition of dimension	" << dimension << " has " << num_variables << " variables!" << std::endl;
		mypause();
	}

	change_size_vec_mp(this->sphere_center, num_variables-1); //destructive resize
	sphere_center->size = num_variables-1;
	
	
	FILE *IN = safe_fopen_read(bounding_sphere_filename);
	
	mpf_inp_str(sphere_radius->r, IN, 10);
	mpf_set_str(sphere_radius->i,"0",10);
	
	
	for (int jj=1; jj<num_variables; jj++) {
		mpf_inp_str(sphere_center->coord[jj-1].r, IN, 10);
		mpf_set_str(sphere_center->coord[jj-1].i,"0",10);
	}
	
	//TODO: check this line for validity
	int tmp_num_vars;
	fscanf(IN,"%d",&tmp_num_vars); scanRestOfLine(IN);
	
	
	fclose(IN);
	
	
	have_sphere_radius = true;
	
	
	return SUCCESSFUL;
}



void decomposition::compute_sphere_bounds(const witness_set & W_crit)
{
	
#ifdef functionentry_output
	std::cout << "decomposition::compute_sphere_bounds" << std::endl;
#endif

	
	int num_vars = W_crit.num_natural_vars-1;
	
	change_size_vec_mp(this->sphere_center, num_vars); //destructive resize
	sphere_center->size = num_vars;
	
	if (W_crit.num_points == 0) {
		set_zero_mp(sphere_radius);
		mpf_set_str(sphere_radius->r, "3.0",10);
		for (int ii=0; ii<num_vars; ii++) {
			set_zero_mp(&sphere_center->coord[ii]);
		}
		return;
	}
	
	vec_mp(temp_vec); init_vec_mp2(temp_vec,0,1024);
	if (W_crit.num_points==1)
	{
		set_zero_mp(sphere_radius);
		mpf_set_str(sphere_radius->r, "3.0",10);
		dehomogenize(&temp_vec, W_crit.pts_mp[0]);
		
		for (int ii=0; ii<num_vars; ii++) {
			set_mp(&sphere_center->coord[ii], &temp_vec->coord[ii]);
		}
#ifdef thresholding
		real_threshold(sphere_center, 1e-13);
#endif
		clear_vec_mp(temp_vec);
		return;
	}
	
	
	//	W_crit.print_to_screen();
	
	
	
	
	comp_mp temp_rad; init_mp2(temp_rad,1024);
	set_zero_mp(temp_rad);
	
	set_one_mp(sphere_radius);
	neg_mp(sphere_radius,sphere_radius); // set to impossibly low value.
	
	
	comp_mp temp_mp;  init_mp2(temp_mp,1024);
	
	
	vec_mp(cumulative_sum); init_vec_mp2(cumulative_sum,num_vars,1024);
	cumulative_sum->size = num_vars;
	
	for (int ii=0; ii<num_vars; ii++) {
		set_zero_mp(&cumulative_sum->coord[ii]);
	}
	
	
	for (int ii=0; ii<W_crit.num_points; ii++) {
		dehomogenize(&temp_vec, W_crit.pts_mp[ii], num_vars+1);
		temp_vec->size = num_vars;
		add_vec_mp(cumulative_sum, cumulative_sum, temp_vec);
	}
	
	set_zero_mp(temp_mp);
	mpf_set_d(temp_mp->r, double(W_crit.num_points));
	
	
	
	for (int ii=0; ii<num_vars; ii++) {
		div_mp(&sphere_center->coord[ii], &cumulative_sum->coord[ii], temp_mp);
	}
	
	
	
	for (int ii=0; ii<W_crit.num_points; ii++) {
		dehomogenize(&temp_vec, W_crit.pts_mp[ii], num_vars+1);
		temp_vec->size = num_vars;
		vec_sub_mp(temp_vec, temp_vec, sphere_center);
		
		
		twoNormVec_mp(temp_vec, temp_mp);
		mpf_abs_mp(temp_rad->r, temp_mp);
		
		
		if (mpf_cmp(sphere_radius->r, temp_rad->r) < 0){
			set_mp(sphere_radius, temp_rad);
		}
	}
	
	
	
	mpf_set_str(temp_rad->r,"2.0",10);
	mpf_set_str(temp_rad->i,"0.0",10);
	mul_mp(sphere_radius,temp_rad,sphere_radius);  // double the radius to be safe.
	
	
	clear_mp(temp_mp); clear_mp(temp_rad);
	clear_vec_mp(temp_vec); clear_vec_mp(cumulative_sum);
	
	
	
	this->have_sphere_radius = true;
#ifdef thresholding
	real_threshold(sphere_center,1e-13);
#endif

	return;
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
	
	W.print_to_file(base / "witness_set");
	
	
// TODO:  this should be a write call, not a copy!
	if (input_filename.filename().string().compare("unset")) {
		copyfile(input_filename, base / input_filename.filename());
	}
	else{
//		std::cout << "not copying inputfile because name was unset -- " << input_filename << std::endl;
	}
	
	
	this->print(base); // using polymorphism and virtualism here!
	
	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%s\n",base.c_str());
	fprintf(OUT,"%d\n",2);//remove this
	fprintf(OUT,"%d\n",dimension);
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
	buffer2 = new int[9];
	buffer2[0] = num_variables;
	buffer2[1] = dimension;
	buffer2[2] = component_num;
	buffer2[3] = num_curr_projections;

	buffer2[4] = num_patches;
	buffer2[5] = have_sphere_radius;
	int strleng = input_filename.string().size() + 1;
	buffer2[6] = strleng;
	buffer2[7] = counters.size();
	buffer2[8] = indices.size();
	MPI_Send(buffer2, 9, MPI_INT, target, 6, MPI_COMM_WORLD);
	delete [] buffer2;
	
	
	
	if (counters.size()>0) {

		//pack and send the counters.
		int * buffer3 = new int[2*counters.size()];
		int cnt = 0;
		for (auto iter = counters.begin(); iter!=counters.begin(); iter++) {
			buffer3[2*cnt] = iter->first;
			buffer3[2*cnt+1] = iter->second;
		}
		MPI_Send(buffer3, 2*counters.size(), MPI_INT, target, 2, MPI_COMM_WORLD);
		delete [] buffer3;
	}
	
	
	

	int * intbuff = new int[2];

	//pack and send the indices
	
	if (indices.size()>0) {
		for (auto iter = indices.begin(); iter!=indices.end(); iter++) { // a std::map<int, std::vectors<ints>>
			
			intbuff[0] = iter->first;
			int num_these_indices = iter->second.size();
			
			intbuff[1] = num_these_indices;
			MPI_Send(intbuff, 2, MPI_INT, target, 4, MPI_COMM_WORLD);
			
			

			if (num_these_indices>0) {
				int * buffer4 = new int[num_these_indices];
				int cnt = 0;
				for (auto jter = iter->second.begin(); jter != iter->second.end(); jter++) {
					buffer4[cnt] = *jter;
					cnt++;
				}
				MPI_Send(buffer4, num_these_indices, MPI_INT, target, 5, MPI_COMM_WORLD);
				delete [] buffer4;
			}
			
		}
	}
	delete [] intbuff;
	
	
	if (num_curr_projections>0) {
		for (int ii=0; ii<num_curr_projections; ii++) {
			send_vec_mp(pi[ii],target);
		}
	}
	
	
	
	
	if ( num_patches>0) {

		for (int ii=0; ii<num_patches; ii++) {
			send_vec_mp(patch_mp[ii],target);
		}
	}
	
	
	
	
	if (have_sphere_radius) {

		send_vec_mp(sphere_center,target);
		send_comp_mp(sphere_radius,target);
	}
	
	
	
	if (strleng>1) {
		char * buffer = new char[strleng];
		memcpy(buffer, input_filename.c_str(), strleng);
		MPI_Send(buffer, strleng, MPI_CHAR, target, 7, MPI_COMM_WORLD);
		delete [] buffer;
	}
	
	
	
	
	randomizer->send(target,mpi_config);



	
	
	
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
	
	
	
	
	
	
	buffer2 = new int[9];
	
	MPI_Recv(buffer2, 9, MPI_INT, source, 6, MPI_COMM_WORLD, &statty_mc_gatty);
	num_variables = buffer2[0];
	dimension = buffer2[1];
	component_num = buffer2[2];
	int temp_num_projections = buffer2[3];
	int temp_num_patches = buffer2[4];
	have_sphere_radius = buffer2[5];
	int strleng = buffer2[6];
	
	int num_counters = buffer2[7];
	int num_indices = buffer2[8];
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
	
	
	if (num_counters>0) {
		int * buffer3 = new int[2*num_counters];
		MPI_Recv(buffer3, 2*num_counters, MPI_INT, source, 2, MPI_COMM_WORLD, &statty_mc_gatty);
		for (int ii=0; ii<num_counters; ii++) {
			counters[buffer3[2*ii]] = buffer3[2*ii+1]; // counters is a map, this is ok.
		}
		delete [] buffer3;
	}
	
	
	int * intbuff = new int[2];
	

	if (num_indices>0) {
		for (int ii=0; ii<num_indices; ii++) {
			
			MPI_Recv(intbuff, 2, MPI_INT, source, 4, MPI_COMM_WORLD, &statty_mc_gatty);
			int indices_index_lol = intbuff[0];
			int num_these_indices = intbuff[1];
			
			
			
			std::vector<int> tempind;
			
			if (num_these_indices>0) {

				int * buffer3 = new int[num_these_indices];
				MPI_Recv(buffer3, num_these_indices, MPI_INT, source, 5, MPI_COMM_WORLD, &statty_mc_gatty);
				for (int jj=0; jj<num_these_indices; jj++) {
					tempind.push_back(buffer3[jj]);
				}
				delete [] buffer3;
			}
			
			indices[indices_index_lol] = tempind;
		}
	}
	
	
	
	
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
	
	
	
	if (have_sphere_radius) {
		receive_vec_mp(sphere_center,source);
		receive_comp_mp(sphere_radius,source);
	}
	
	if (strleng>1) {
		char * buffer = new char[strleng];
		MPI_Recv(buffer, strleng, MPI_CHAR, source, 7, MPI_COMM_WORLD, &statty_mc_gatty);
		
		input_filename = buffer;
		delete [] buffer;
	}
	
	
	
	delete [] intbuff;
	

	clear_vec_mp(tempvec);
	
	randomizer->receive(source,mpi_config);
	
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








