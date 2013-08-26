#include "solver_nullspace_left.hpp"







void nullspace_config::clear()
{
	
	free(randomized_degrees); // the degrees of the randomized functions (not derivatives)
	
	for (int ii=0; ii<ambient_dim-target_dim; ii++) {
		for (int jj=0; jj<max_degree; jj++) {
			clear_vec_mp(starting_linears[ii][jj]);
		}
		free(starting_linears[ii]);
	}
	free(starting_linears);
	
	for (int ii=0; ii<num_additional_linears; ii++){
		clear_vec_mp(additional_linears_starting[ii]);
		clear_vec_mp(additional_linears_terminal[ii]);
	}
	free(additional_linears_terminal);
	free(additional_linears_starting);
	

	for (int ii=0; ii<num_v_linears; ii++) 
		clear_vec_mp(v_linears[ii]);
	free(v_linears);
	
	clear_vec_mp(v_patch);
	
	for (int ii=0; ii<num_projections; ii++)
		clear_vec_mp(target_projection[ii]);
	free(target_projection);
	
	clear_mat_mp(randomizer_matrix);
	clear_mat_mp(post_randomizer_matrix);
}




void nullspace_config::print()
{
	std::stringstream converter;
	
	std::cout << "*******************\n\tns_config:\n\n";

	std::cout << "#x: " << num_x_vars << "\n#v: " << num_v_vars << std::endl;
	
	std::cout << "# additional linears: " << this->num_additional_linears << std::endl;
	std::cout << "# jacobian equations: " << this->num_jac_equations << std::endl;
	std::cout << "# randomized equations: " << this->num_randomized_eqns << std::endl;
	std::cout << "max_degree:\t" << max_degree << std::endl;
	
	print_matrix_to_screen_matlab(this->randomizer_matrix,"R");
	print_matrix_to_screen_matlab(this->post_randomizer_matrix,"S");
	
	
	for (int ii=0; ii<num_additional_linears; ii++) {
		converter << "additional_linear_starting_" << ii+1;
		print_point_to_screen_matlab(additional_linears_starting[ii],converter.str());
		converter.str("");
	}
	
	
	for (int ii=0; ii<num_additional_linears; ii++) {
		converter << "additional_linear_terminal_" << ii+1;
		print_point_to_screen_matlab(additional_linears_terminal[ii],converter.str());
		converter.str("");
	}
	
	
	
	for (int ii=0; ii<num_jac_equations; ii++) {
		for (int jj=0; jj<this->max_degree; jj++) {
			converter << "starting_linear_" << ii+1 << "_" << jj+1;
			print_point_to_screen_matlab(starting_linears[ii][jj],converter.str());
			converter.str("");
		}
	}
	
	for (int ii=0; ii<num_jac_equations; ii++) {
		converter << "v_linear_" << ii+1;
			print_point_to_screen_matlab(v_linears[ii],converter.str());
		converter.str("");
	}
	
	std::cout << "\n********************\n" << std::endl;
	
}


///////////////
//
//   end the nullspace_config
//
/////////////






///////////////
//
//   begin nullspace_eval_data_mp
//
/////////////
void nullspacejac_eval_data_mp::init()
{
	this->is_solution_checker_d = &check_issoln_nullspacejac_d;
	this->is_solution_checker_mp = &check_issoln_nullspacejac_mp;
	this->evaluator_function_d = &nullspacejac_eval_d;
	this->evaluator_function_mp = &nullspacejac_eval_mp;
	this->precision_changer = &change_nullspacejac_eval_prec;
	this->dehomogenizer = &nullspacejac_dehom;
	
	additional_linears_terminal = NULL;
	additional_linears_terminal_full_prec = NULL;
	
	additional_linears_starting = NULL;
	additional_linears_starting_full_prec = NULL;
	
	
	init_mat_mp(post_randomizer_matrix,0,0);  // S, for randomizing the jacobian subsystem down to N-k+\ell-1 equations
	
	starting_linears = NULL; // outer layer should have as many as there are randomized equations
													 // inside layer has number corresponding to randomized_degrees
	starting_linears_full_prec = NULL; // outer layer should have as many as there are randomized equations
																		 // inside layer has number corresponding to randomized_degrees
	
	v_linears = NULL;         // should be as many in here as there are randomized equations
	v_linears_full_prec = NULL;         // should be as many in here as there are randomized equations
	
	init_vec_mp(v_patch,0);
	
	
	
	init_mat_mp(jac_with_proj,0,0);
	
	
	init_mp(perturbation);
	
	init_mp(half);
	
	
	comp_d h;
	
	h->r = 0.5;
	h->i = 0.0;
	
	d_to_mp(half, h);
	
	
	comp_d p; p->r = PERTURBATION_VALUE_mp; p->i = PERTURBATION_VALUE_mp;
	d_to_mp(perturbation,p);
	
	
	target_projection = NULL; //
	target_projection_full_prec = NULL; //
	
	this->randomized_degrees.clear();
	
	if (this->MPType==2) {
		init_mat_mp2(post_randomizer_matrix_full_prec,0,0,1024);  // S, for randomizing the jacobian subsystem down to N-k+\ell-1 equations
		init_vec_mp2(v_patch_full_prec,0,1024);
		init_mat_mp2(jac_with_proj_full_prec,0,0,1024);
		init_mp2(perturbation_full_prec,1024);
		init_mp2(half_full_prec,1024);
		d_to_mp(perturbation_full_prec,p);
		d_to_mp(half_full_prec, h);
	}
	
	
}


int nullspacejac_eval_data_mp::send(parallelism_config & mpi_config)
{

	int solver_choice = NULLSPACE;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
	// send the confirmation integer, to ensure that we are sending the correct type.
	
	//send the base class stuff.
	solver_mp::send(mpi_config);
	
	std::cout << "master done sending solver_d" << std::endl;
	
	int *buffer = new int[12];
	
	buffer[0] = num_additional_linears;
	buffer[1] = num_jac_equations;
	buffer[2] = max_degree;
	buffer[3] = num_v_linears;
	buffer[4] = num_projections;
	
	buffer[5] = num_x_vars;
	buffer[6] = num_v_vars;
	
	buffer[7] = target_dim;
	buffer[8] = ambient_dim;
	buffer[9] = target_crit_codim;
	
	buffer[10] = num_randomized_eqns;
	buffer[11] = num_natural_vars;
	// now can actually send the data.
	
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	delete[] buffer;
	buffer = new int[num_randomized_eqns];
	for (int ii=0; ii<num_randomized_eqns; ii++) {
		buffer[ii] = randomized_degrees[ii];
	}
	MPI_Bcast(buffer, num_randomized_eqns, MPI_INT, 0, mpi_config.my_communicator);
	delete[] buffer;
	
	
	if (this->MPType==2){
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				// receive the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal_full_prec[ii], mpi_config.id(), mpi_config.head());
				
				// receive the full precision starting additional linear
				bcast_vec_mp(this->additional_linears_starting_full_prec[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else {} // num_additional_linears == 0
		
		
		// recieve the post-randomizer-matrix
		bcast_mat_mp(post_randomizer_matrix_full_prec, mpi_config.id(), mpi_config.head());
		
		
		if (num_jac_equations>0) {
			for (int ii=0; ii<num_jac_equations; ii++) {
				for (int jj=0; jj<max_degree; jj++) {
					bcast_vec_mp(starting_linears_full_prec[ii][jj], mpi_config.id(), mpi_config.head());
				}
			}
		}
		else{}
		
		
		
		
		if (num_v_linears>0) {
			for (int ii=0; ii<num_v_linears; ii++) {
				bcast_vec_mp(v_linears_full_prec[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else{}
		
		bcast_vec_mp(v_patch_full_prec, mpi_config.id(), mpi_config.head());
		bcast_mat_mp(jac_with_proj_full_prec, mpi_config.id(), mpi_config.head());
		bcast_comp_mp(perturbation_full_prec, mpi_config.id(), mpi_config.head());
		bcast_comp_mp(half_full_prec, mpi_config.id(), mpi_config.head());

		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				bcast_vec_mp(target_projection_full_prec[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else{}
	}
	else {
		
		if (num_additional_linears>0) {
			for (int ii=0; ii<num_additional_linears; ii++) {
				// receive the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
				
				// receive the full precision starting additional linear
				bcast_vec_mp(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else {} // num_additional_linears == 0
		
		
		// recieve the post-randomizer-matrix
		bcast_mat_mp(post_randomizer_matrix, mpi_config.id(), mpi_config.head());
		
	
		for (int ii=0; ii<num_jac_equations; ii++) {
			for (int jj=0; jj<max_degree; jj++) {
				bcast_vec_mp(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
			}
		}
	
		for (int ii=0; ii<num_v_linears; ii++) {
			bcast_vec_mp(v_linears[ii], mpi_config.id(), mpi_config.head());
		}

		
		bcast_vec_mp(v_patch, mpi_config.id(), mpi_config.head());
		bcast_mat_mp(jac_with_proj, mpi_config.id(), mpi_config.head());
		bcast_comp_mp(perturbation, mpi_config.id(), mpi_config.head());
		bcast_comp_mp(half, mpi_config.id(), mpi_config.head());
		
		if (num_projections>0) {
			for (int ii=0; ii<num_projections; ii++) {
				bcast_vec_mp(target_projection[ii], mpi_config.id(), mpi_config.head());
			}
		}
		else{}
	}
	
	
	return SUCCESSFUL;
}




int nullspacejac_eval_data_mp::receive(parallelism_config & mpi_config)
{
	int *buffer = new int[12];
	MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (buffer[0] != NULLSPACE) {
		std::cout << "worker failed to confirm it is receiving the NULLSPACE type eval data" << std::endl;
		mpi_config.abort(777);
	}
	
	solver_mp::receive(mpi_config);

	// now can actually receive the data from whoever.
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	num_additional_linears = buffer[0];
	num_jac_equations = buffer[1];
	max_degree = buffer[2];
	num_v_linears = buffer[3];
	num_projections = buffer[4];
	
	num_x_vars = buffer[5];
	num_v_vars = buffer[6];
	
	target_dim = buffer[7];
	ambient_dim = buffer[8];
	target_crit_codim = buffer[9];
	
	num_randomized_eqns = buffer[10];
	num_natural_vars = buffer[11];
	
	delete[] buffer;
	buffer = new int[num_randomized_eqns];
	
	MPI_Bcast(buffer, num_randomized_eqns, MPI_INT, 0, mpi_config.my_communicator);
	for (int ii=0; ii<num_randomized_eqns; ii++) {
		randomized_degrees.push_back(buffer[ii]);
	}
	delete[] buffer;
	
	
	
	if (this->MPType==2) {
		if (num_additional_linears>0) {
			
			this->additional_linears_terminal						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_terminal_full_prec = (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting_full_prec = (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				init_vec_mp2(this->additional_linears_terminal_full_prec[ii],0,1024);
				init_vec_mp(this->additional_linears_terminal[ii],0);
				// receive the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->additional_linears_terminal[ii],this->additional_linears_terminal_full_prec[ii]);

				// receive the full precision starting additional linear
				init_vec_mp2(this->additional_linears_starting_full_prec[ii],0,1024);
				init_vec_mp(this->additional_linears_starting[ii],0);
				bcast_vec_mp(this->additional_linears_starting_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->additional_linears_starting[ii],this->additional_linears_starting_full_prec[ii]);
				
			}
		}
		else {} // num_additional_linears == 0
		
		
		// recieve the post-randomizer-matrix
		// recieve the post-randomizer-matrix
		init_mat_mp(post_randomizer_matrix,num_jac_equations,num_x_vars-1);
		post_randomizer_matrix->rows = num_jac_equations;
		post_randomizer_matrix->cols = num_x_vars-1;
		
		init_mat_mp(post_randomizer_matrix_full_prec,num_jac_equations,num_x_vars-1);
		post_randomizer_matrix_full_prec->rows = num_jac_equations;
		post_randomizer_matrix_full_prec->cols = num_x_vars-1;
		
		bcast_mat_mp(post_randomizer_matrix_full_prec, mpi_config.id(), mpi_config.head());
		mat_cp_mp(this->post_randomizer_matrix, post_randomizer_matrix_full_prec);
		
		
		if (num_jac_equations>0) {
			this->starting_linears_full_prec = (vec_mp **) br_malloc(num_jac_equations*sizeof(vec_mp *));
			this->starting_linears = (vec_mp **) br_malloc(num_jac_equations*sizeof(vec_mp *));
			
			for (int ii=0; ii<num_jac_equations; ii++) {
				
				this->starting_linears_full_prec[ii] = (vec_mp *) br_malloc(max_degree*sizeof(vec_mp ));
				this->starting_linears[ii] = (vec_mp *) br_malloc(max_degree*sizeof(vec_mp ));
				for (int jj=0; jj<max_degree; jj++) {
					init_vec_mp(this->starting_linears[ii][jj],0);
					init_vec_mp2(this->starting_linears_full_prec[ii][jj],0,1024);
					
					// recieve the starting linearsin full prec and convert
					bcast_vec_mp(starting_linears_full_prec[ii][jj], mpi_config.id(), mpi_config.head());
					vec_cp_mp(this->starting_linears[ii][jj],starting_linears_full_prec[ii][jj]);
				}
			}
		}
		else{}
		
		if (num_v_linears>0) {
			this->v_linears = (vec_mp *) br_malloc(num_v_linears*sizeof(vec_mp));
			this->v_linears_full_prec = (vec_mp *) br_malloc(num_v_linears*sizeof(vec_mp));
			for (int ii=0; ii<num_v_linears; ii++) {
				init_vec_mp(this->v_linears[ii],0);
				init_vec_mp2(this->v_linears_full_prec[ii],0,1024);
				// receive the full precision v_linears
				bcast_vec_mp(v_linears_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->v_linears[ii],v_linears_full_prec[ii]);
			}
		}
		else{}
		
		init_vec_mp2(v_patch_full_prec,0,1024);
		init_vec_mp(v_patch,0);
		bcast_vec_mp(v_patch_full_prec, mpi_config.id(), mpi_config.head());
		vec_cp_mp(this->v_patch, v_patch_full_prec);
		
		init_mat_mp2(jac_with_proj_full_prec, num_x_vars-1, num_randomized_eqns+num_projections,1024);
		
		
		init_mat_mp(jac_with_proj, num_x_vars-1, num_randomized_eqns+num_projections);
		jac_with_proj->rows = jac_with_proj_full_prec->rows = num_x_vars-1;
		jac_with_proj->cols = jac_with_proj_full_prec->cols = num_randomized_eqns+num_projections;
		
		bcast_mat_mp(jac_with_proj_full_prec, mpi_config.id(), mpi_config.head());
		mat_cp_mp(this->jac_with_proj, jac_with_proj_full_prec);
		
		init_mp2(perturbation_full_prec,1024);
		init_mp(perturbation);
		bcast_comp_mp(perturbation_full_prec, mpi_config.id(), mpi_config.head());
		set_mp(this->perturbation,perturbation_full_prec);
		
		init_mp2(half_full_prec,1024);
		init_mp(half);
		bcast_comp_mp(half_full_prec, mpi_config.id(), mpi_config.head());
		set_mp(this->half, half_full_prec);
		
		
		
		if (num_projections>0) {
			this->target_projection = (vec_mp *) br_malloc(num_projections*sizeof(vec_mp));
			this->target_projection_full_prec = (vec_mp *) br_malloc(num_projections*sizeof(vec_mp));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_mp(this->target_projection[ii],0);
				init_vec_mp2(this->target_projection[ii],0,1024);
				
				bcast_vec_mp(target_projection_full_prec[ii], mpi_config.id(), mpi_config.head());
				vec_cp_mp(this->target_projection[ii],target_projection_full_prec[ii]);
			}
		}
		else{}
	
	
	}
	else{ // MPType == 1
		if (num_additional_linears>0) {
			
			this->additional_linears_terminal						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			this->additional_linears_starting						= (vec_mp *) br_malloc(num_additional_linears*sizeof(vec_mp));
			
			for (int ii=0; ii<num_additional_linears; ii++) {
				init_vec_mp2(this->additional_linears_terminal[ii],0,1024);
				// receive the full precision terminal additional linear
				bcast_vec_mp(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
				
				// receive the full precision starting additional linear
				init_vec_mp2(this->additional_linears_starting[ii],0,1024);
				bcast_vec_mp(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
				
			}
		}
		else {} // num_additional_linears == 0
		
		std::cout << "worker got additional linears" << std::endl;
		
		// recieve the post-randomizer-matrix
		init_mat_mp(post_randomizer_matrix,num_jac_equations,num_x_vars-1);
		post_randomizer_matrix->rows = num_jac_equations;
		post_randomizer_matrix->cols = num_x_vars-1;
		
		bcast_mat_mp(post_randomizer_matrix, mpi_config.id(), mpi_config.head());
		
		std::cout << "worker got post-randomizer matrix" << std::endl;
		
		std::cout << "num_jac " << num_jac_equations << std::endl;
		this->starting_linears = (vec_mp **) br_malloc(num_jac_equations*sizeof(vec_mp *));
		
		for (int ii=0; ii<num_jac_equations; ii++) {
			this->starting_linears[ii] = (vec_mp *) br_malloc(max_degree*sizeof(vec_mp ));
			for (int jj=0; jj<max_degree; jj++) {
				init_vec_mp(this->starting_linears[ii][jj],num_x_vars);  starting_linears[ii][jj]->size = num_x_vars;
				std::cout << "worker getting starting_linear[" << ii << "][" << jj << "]" << std::endl;
				// recieve the starting linearsin full prec and convert
				bcast_vec_mp(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
			}
		}

		this->v_linears = (vec_mp *) br_malloc(num_v_linears*sizeof(vec_mp));
		for (int ii=0; ii<num_v_linears; ii++) {
			init_vec_mp(this->v_linears[ii],0);
			// receive the full precision v_linears
			bcast_vec_mp(v_linears[ii], mpi_config.id(), mpi_config.head());
		}
		
		init_vec_mp(v_patch,0);
		bcast_vec_mp(v_patch, mpi_config.id(), mpi_config.head());
		std::cout << "a";
		
		init_mat_mp(jac_with_proj,num_x_vars-1,num_randomized_eqns+num_projections);
		jac_with_proj->rows = num_x_vars-1;
		jac_with_proj->cols = num_randomized_eqns+num_projections;
		bcast_mat_mp(jac_with_proj, mpi_config.id(), mpi_config.head());
		std::cout << "b";
		bcast_comp_mp(perturbation, mpi_config.id(), mpi_config.head());
		std::cout << "c";
		bcast_comp_mp(half, mpi_config.id(), mpi_config.head());
		std::cout << "d";
		
		if (num_projections>0) {
			this->target_projection = (vec_mp *) br_malloc(num_projections*sizeof(vec_mp));
			
			for (int ii=0; ii<num_projections; ii++) {
				init_vec_mp(this->target_projection[ii],0);
				bcast_vec_mp(target_projection[ii], mpi_config.id(), mpi_config.head());
			}
		}

	}

	
	return SUCCESSFUL;
}

int nullspacejac_eval_data_mp::setup(prog_t * _SLP,
																		 nullspace_config *ns_config,
																		 witness_set & W,
																		 solver_configuration & solve_options)
{

	verbose_level = solve_options.verbose_level;
	
	solver_mp::setup(_SLP);
	
	generic_setup_patch(&patch,W);

	
	
	
	if (solve_options.use_gamma_trick==1)
		get_comp_rand_mp(this->gamma); // set gamma to be random complex value
	else{
		set_one_mp(this->gamma);
	}
	
	comp_d temp;
	if (this->MPType==2) {
		if (solve_options.use_gamma_trick==1){
			get_comp_rand_rat(temp, this->gamma, this->gamma_rat, 64, solve_options.T.AMP_max_prec, 0, 0);
		}
		else{
			set_one_mp(this->gamma);
			set_one_rat(this->gamma_rat);
		}
	}
	
	num_jac_equations = ns_config->num_jac_equations;
	target_dim = ns_config->target_dim;
	ambient_dim = ns_config->ambient_dim;
	target_crit_codim = ns_config->target_crit_codim;
	
	num_natural_vars = W.num_variables - W.num_synth_vars - ns_config->num_v_vars;
	num_v_vars = ns_config->num_v_vars;
	num_x_vars = ns_config->num_x_vars;
	
	num_v_linears = ns_config->num_v_linears;   //
	num_variables = ns_config->num_x_vars + ns_config->num_v_vars;
	
	num_projections = ns_config->num_projections;
	
	num_additional_linears = ns_config->num_additional_linears;
	
	num_randomized_eqns = ns_config->num_randomized_eqns;
	max_degree = ns_config->max_degree;
	randomized_degrees.resize(ns_config->num_randomized_eqns);
	for (int ii=0; ii<ns_config->randomizer_matrix->rows; ii++)
		randomized_degrees[ii] = ns_config->randomized_degrees[ii]; // store the full degree (not derivative).
	
	
	
	
	
	mat_cp_mp(randomizer_matrix,
						ns_config->randomizer_matrix);
	
	
	mat_cp_mp(post_randomizer_matrix,
						ns_config->post_randomizer_matrix);
	
	
	
	// set up the vectors to hold the linears.

	target_projection = (vec_mp *)br_malloc(ns_config->num_projections * sizeof(vec_mp));

	init_mat_mp(jac_with_proj, ns_config->num_x_vars-1,ns_config->num_v_vars);
	jac_with_proj->rows = ns_config->num_x_vars-1;
	jac_with_proj->cols = ns_config->num_v_vars;
	
	int offset = ns_config->num_randomized_eqns;
	for (int ii=0; ii<ns_config->num_projections; ii++) {
		init_vec_mp(target_projection[ii],ns_config->num_x_vars);
		target_projection[ii]->size =  ns_config->num_x_vars;
		vec_cp_mp(target_projection[ii], ns_config->target_projection[ii]);
		
		for (int jj=1; jj<ns_config->num_x_vars; jj++) {
			set_mp(&jac_with_proj->entry[jj-1][ii+offset], &ns_config->target_projection[ii]->coord[jj]);
		}
	}
	

	vec_cp_mp(v_patch, ns_config->v_patch);
	

	
	
	this->v_linears = (vec_mp *) br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
	for (int ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_mp(v_linears[ii],ns_config->num_v_vars);
		v_linears[ii]->size = ns_config->num_v_vars;
		vec_cp_mp(v_linears[ii], ns_config->v_linears[ii]);
	}
	
	
	additional_linears_terminal = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
	additional_linears_starting = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
	
	for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
		init_vec_mp(additional_linears_terminal[ii], ns_config->num_x_vars);
		additional_linears_terminal[ii]->size = ns_config->num_x_vars;
		vec_cp_mp(additional_linears_terminal[ii],ns_config->additional_linears_terminal[ii]);
		
		init_vec_mp(additional_linears_starting[ii], ns_config->num_x_vars);
		additional_linears_starting[ii]->size = ns_config->num_x_vars;
		vec_cp_mp(additional_linears_starting[ii],ns_config->additional_linears_starting[ii]);
	}
	
	starting_linears = (vec_mp **)br_malloc(ns_config->num_jac_equations*sizeof(vec_mp *));
	
	for (int ii=0; ii<ns_config->num_jac_equations; ++ii) {
		starting_linears[ii] = (vec_mp *)br_malloc(ns_config->max_degree*sizeof(vec_mp));
		for (int jj=0; jj<ns_config->max_degree; jj++) {
			init_vec_mp(starting_linears[ii][jj],W.num_variables);
			starting_linears[ii][jj]->size = W.num_variables;
			
			vec_cp_mp(starting_linears[ii][jj], ns_config->starting_linears[ii][jj]);
		}
	}
	
	
	
	if (this->MPType==2) {
		mat_cp_mp(randomizer_matrix_full_prec, ns_config->randomizer_matrix);
		
		
		mat_cp_mp(post_randomizer_matrix_full_prec, ns_config->post_randomizer_matrix);
		
		
		
		// set up the vectors to hold the linears.
		
		target_projection_full_prec = (vec_mp *)br_malloc(ns_config->num_projections * sizeof(vec_mp));
		
		init_mat_mp2(jac_with_proj_full_prec, ns_config->num_x_vars-1,ns_config->num_v_vars,solve_options.T.AMP_max_prec);
		jac_with_proj_full_prec->rows = ns_config->num_x_vars-1;
		jac_with_proj_full_prec->cols = ns_config->num_v_vars;
		
		int offset = ns_config->num_randomized_eqns;
		for (int ii=0; ii<ns_config->num_projections; ii++) {
			init_vec_mp2(target_projection_full_prec[ii],ns_config->num_x_vars,solve_options.T.AMP_max_prec);
			target_projection_full_prec[ii]->size =  ns_config->num_x_vars;
			vec_cp_mp(target_projection_full_prec[ii], ns_config->target_projection[ii]);
			
			for (int jj=1; jj<ns_config->num_x_vars; jj++) {
				set_mp(&jac_with_proj_full_prec->entry[jj-1][ii+offset], &ns_config->target_projection[ii]->coord[jj]);
			}
		}
		
		
		vec_cp_mp(v_patch_full_prec, ns_config->v_patch);
		
		
		
		
		this->v_linears_full_prec = (vec_mp *) br_malloc(ns_config->num_v_linears*sizeof(vec_mp));
		for (int ii=0; ii<ns_config->num_v_linears; ii++) {
			init_vec_mp2(v_linears_full_prec[ii],ns_config->num_v_vars,solve_options.T.AMP_max_prec);
			v_linears_full_prec[ii]->size = ns_config->num_v_vars;
			vec_cp_mp(v_linears_full_prec[ii], ns_config->v_linears[ii]);
		}
		
		
		additional_linears_terminal_full_prec = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
		additional_linears_starting_full_prec = (vec_mp *) br_malloc(ns_config->num_additional_linears*sizeof(vec_mp));
		
		for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
			init_vec_mp2(additional_linears_terminal_full_prec[ii], ns_config->num_x_vars,solve_options.T.AMP_max_prec);
			additional_linears_terminal_full_prec[ii]->size = ns_config->num_x_vars;
			vec_cp_mp(additional_linears_terminal_full_prec[ii],ns_config->additional_linears_terminal[ii]);
			
			init_vec_mp2(additional_linears_starting_full_prec[ii], ns_config->num_x_vars,solve_options.T.AMP_max_prec);
			additional_linears_starting_full_prec[ii]->size = ns_config->num_x_vars;
			vec_cp_mp(additional_linears_starting_full_prec[ii],ns_config->additional_linears_starting[ii]);
		}
		
		starting_linears_full_prec = (vec_mp **)br_malloc(ns_config->num_jac_equations*sizeof(vec_mp *));
		
		for (int ii=0; ii<ns_config->num_jac_equations; ++ii) {
			starting_linears_full_prec[ii] = (vec_mp *)br_malloc(ns_config->max_degree*sizeof(vec_mp));
			for (int jj=0; jj<ns_config->max_degree; jj++) {
				init_vec_mp2(starting_linears_full_prec[ii][jj],W.num_variables,solve_options.T.AMP_max_prec);
				starting_linears_full_prec[ii][jj]->size = W.num_variables;
				
				vec_cp_mp(starting_linears_full_prec[ii][jj], ns_config->starting_linears[ii][jj]);
			}
		}
	}
	
							
	return SUCCESSFUL;
}

///////////////
//
//   end nullspace_eval_data_mp
//
/////////////



















///////////////
//
//   begin nullspace_eval_data_d
//
/////////////

void nullspacejac_eval_data_d::init()
{

	if (this->MPType==2)
		this->BED_mp = new nullspacejac_eval_data_mp(2);
	else
		this->BED_mp = NULL;
	
	this->is_solution_checker_d = &check_issoln_nullspacejac_d;
	this->is_solution_checker_mp = &check_issoln_nullspacejac_mp;
	this->evaluator_function_d = &nullspacejac_eval_d;
	this->evaluator_function_mp = &nullspacejac_eval_mp;
	this->precision_changer = &change_nullspacejac_eval_prec;
	this->dehomogenizer = &nullspacejac_dehom;
	
	additional_linears_terminal = NULL;
	
	additional_linears_starting = NULL;
	
	
	init_mat_d(post_randomizer_matrix,0,0);  // S, for randomizing the jacobian subsystem down to N-k+\ell-1 equations

	starting_linears = NULL; // outer layer should have as many as there are randomized equations
													 // inside layer has number corresponding to randomized_degrees
	
	
	
	v_linears = NULL;         // should be as many in here as there are randomized equations
	
	init_vec_d(v_patch,0);

	init_mat_d(jac_with_proj,0,0);

	int id;
	MPI_Comm_rank(MPI_COMM_WORLD, &id);

	half->r = 0.5; half->i = 0.0;
	
	perturbation->r = PERTURBATION_VALUE; perturbation->i = PERTURBATION_VALUE;

	
	target_projection = NULL; //

	
	this->randomized_degrees.clear();
}


int nullspacejac_eval_data_d::send(parallelism_config & mpi_config)
{
	
	int solver_choice = NULLSPACE;
	std::cout << "master bcasting solver choice " << solver_choice << std::endl;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
	// send the confirmation integer, to ensure that we are sending the correct type.
	
	//send the base class stuff.
	solver_d::send(mpi_config);

	std::cout << "master done sending solver_d" << std::endl;
	
	int *buffer = new int[12];
	
	buffer[0] = num_additional_linears;
	buffer[1] = num_jac_equations;
	buffer[2] = max_degree;
	buffer[3] = num_v_linears;
	buffer[4] = num_projections;
	
	buffer[5] = num_x_vars;
	buffer[6] = num_v_vars;
	
	buffer[7] = target_dim;
	buffer[8] = ambient_dim;
	buffer[9] = target_crit_codim;
	
	buffer[10] = num_randomized_eqns;
	buffer[11] = num_natural_vars;
	// now can actually send the data.
	
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	delete[] buffer;
	buffer = new int[num_randomized_eqns];
	for (int ii=0; ii<num_randomized_eqns; ii++) {
		buffer[ii] = randomized_degrees[ii];
	}
	MPI_Bcast(buffer, num_randomized_eqns, MPI_INT, 0, mpi_config.my_communicator);
	delete[] buffer;
	
	if (num_additional_linears>0) {
		for (int ii=0; ii<num_additional_linears; ii++) {
			// receive the full precision terminal additional linear
			bcast_vec_d(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
			
			// receive the full precision starting additional linear
			bcast_vec_d(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else {} // num_additional_linears == 0
	
	
	// recieve the post-randomizer-matrix
	bcast_mat_d(post_randomizer_matrix, mpi_config.id(), mpi_config.head());
	
	
	if (num_jac_equations>0) {
		for (int ii=0; ii<num_jac_equations; ii++) {
			for (int jj=0; jj<max_degree; jj++) {
				bcast_vec_d(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
			}
		}
	}
	else{}
	
	
	
	
	if (num_v_linears>0) {
		for (int ii=0; ii<num_v_linears; ii++) {
			bcast_vec_d(v_linears[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	bcast_vec_d(v_patch, mpi_config.id(), mpi_config.head());
	bcast_mat_d(jac_with_proj, mpi_config.id(), mpi_config.head());

	
	if (num_projections>0) {
		for (int ii=0; ii<num_projections; ii++) {
			bcast_vec_d(target_projection[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	if (this->MPType==2) {
		this->BED_mp->send(mpi_config);
	}
	
	return SUCCESSFUL;
}

int nullspacejac_eval_data_d::receive(parallelism_config & mpi_config)
{
	int *buffer = new int[12];
	MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (buffer[0] != NULLSPACE){
		std::cout << "worker failed to confirm it is receiving the nullspace type eval data" << std::endl;
		mpi_config.abort(777);
	}

	solver_d::receive(mpi_config);
	
	// now can actually receive the data from whoever.
	
	
	
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	
	num_additional_linears = buffer[0];
	num_jac_equations = buffer[1];
	max_degree = buffer[2];
	num_v_linears = buffer[3];
	num_projections = buffer[4];
	
	num_x_vars = buffer[5];
	num_v_vars = buffer[6];
	
	target_dim = buffer[7];
	ambient_dim = buffer[8];
	target_crit_codim = buffer[9];
	
	num_randomized_eqns = buffer[10];
	num_natural_vars = buffer[11];
	
	
	
	delete[] buffer;
	buffer = new int[num_randomized_eqns];
	
	MPI_Bcast(buffer, num_randomized_eqns, MPI_INT, 0, mpi_config.my_communicator);
	for (int ii=0; ii<num_randomized_eqns; ii++) {
		randomized_degrees.push_back(buffer[ii]);
	}
	delete[] buffer;
	
	
	if (num_additional_linears>0) {
		
		this->additional_linears_terminal						= (vec_d *) br_malloc(num_additional_linears*sizeof(vec_d));
		this->additional_linears_starting						= (vec_d *) br_malloc(num_additional_linears*sizeof(vec_d));
		
		for (int ii=0; ii<num_additional_linears; ii++) {
			
			// receive the full precision terminal additional linear
			bcast_vec_d(this->additional_linears_terminal[ii], mpi_config.id(), mpi_config.head());
			
			// receive the full precision starting additional linear
			bcast_vec_d(this->additional_linears_starting[ii], mpi_config.id(), mpi_config.head());
			
		}
	}
	else {} // num_additional_linears == 0
	
	
	// recieve the post-randomizer-matrix
	bcast_mat_d(post_randomizer_matrix, mpi_config.id(), mpi_config.head());
	
	
	if (num_jac_equations>0) {
		this->starting_linears = (vec_d **) br_malloc(num_jac_equations*sizeof(vec_d *));
		
		for (int ii=0; ii<num_jac_equations; ii++) {
			
			this->starting_linears[ii] = (vec_d *) br_malloc(max_degree*sizeof(vec_d ));
			for (int jj=0; jj<max_degree; jj++) {
				init_vec_d(this->starting_linears[ii][jj],0);
				
				// recieve the starting linearsin full prec and convert
				bcast_vec_d(starting_linears[ii][jj], mpi_config.id(), mpi_config.head());
			}
		}
	}
	else{}
	
	if (num_v_linears>0) {
		this->v_linears = (vec_d *) br_malloc(num_v_linears*sizeof(vec_d));
		for (int ii=0; ii<num_v_linears; ii++) {
			init_vec_d(this->v_linears[ii],0);
			// receive the full precision v_linears
			bcast_vec_d(v_linears[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	bcast_vec_d(v_patch, mpi_config.id(), mpi_config.head());
	
	bcast_mat_d(jac_with_proj, mpi_config.id(), mpi_config.head());
	

	
	
	if (num_projections>0) {
		this->target_projection = (vec_d *) br_malloc(num_projections*sizeof(vec_d));
		
		for (int ii=0; ii<num_projections; ii++) {
			init_vec_d(this->target_projection[ii],0);
			
			bcast_vec_d(target_projection[ii], mpi_config.id(), mpi_config.head());
		}
	}
	else{}
	
	std::cout << "worker entering mp receive" << std::endl;
	if (this->MPType==2) {
		this->BED_mp->receive(mpi_config);
	}
	
	return SUCCESSFUL;
}




int nullspacejac_eval_data_d::setup(prog_t * _SLP,
																		nullspace_config *ns_config,
																		witness_set & W,
																		solver_configuration & solve_options)
{
	
	solver_d::setup(_SLP);
	
	verbose_level = solve_options.verbose_level;
	
	generic_setup_patch(&patch,W);
	

	
	
	
	if (solve_options.use_gamma_trick==1)
		get_comp_rand_d(this->gamma); // set gamma to be random complex value
	else
		set_one_d(this->gamma);
	
	
	num_jac_equations = ns_config->num_jac_equations;
	target_dim = ns_config->target_dim;
	ambient_dim = ns_config->ambient_dim;
	target_crit_codim = ns_config->target_crit_codim;
	
	num_natural_vars = W.num_variables - W.num_synth_vars - ns_config->num_v_vars;
	num_v_vars = ns_config->num_v_vars;
	num_x_vars = ns_config->num_x_vars;
	
	num_randomized_eqns = ns_config->num_randomized_eqns;
	max_degree = ns_config->max_degree;
	randomized_degrees.resize(ns_config->num_randomized_eqns);
	for (int ii=0; ii<ns_config->randomizer_matrix->rows; ii++)
		randomized_degrees[ii] = ns_config->randomized_degrees[ii]; // store the full degree (not derivative).
	
	
	num_v_linears = ns_config->num_v_linears;   //
	
	
	num_variables = ns_config->num_x_vars + ns_config->num_v_vars;
	
  
	mat_mp_to_d(randomizer_matrix,
						ns_config->randomizer_matrix);
	
	
	//  THE STUFF PROPRIETARY TO THIS METHOD
	mat_mp_to_d(post_randomizer_matrix,
						ns_config->post_randomizer_matrix);
	
	
	
	// set up the vectors to hold the linears.
	num_projections = ns_config->num_projections;
	target_projection = (vec_d *)br_malloc(ns_config->num_projections * sizeof(vec_d));
	
	init_mat_d(jac_with_proj,ns_config->num_x_vars-1,ns_config->num_v_vars);
	jac_with_proj->rows = ns_config->num_x_vars-1;
	jac_with_proj->cols = ns_config->num_v_vars;
	
	int offset = ns_config->num_randomized_eqns;
	for (int ii=0; ii<ns_config->num_projections; ii++) {
		init_vec_d(target_projection[ii],ns_config->num_x_vars);
		target_projection[ii]->size =  ns_config->num_x_vars;
		vec_mp_to_d(target_projection[ii], ns_config->target_projection[ii]);
		
		for (int jj=1; jj<ns_config->num_x_vars; jj++) {
			mp_to_d(&jac_with_proj->entry[jj-1][ii+offset], &ns_config->target_projection[ii]->coord[jj]);
		}
	}
	
	
	vec_mp_to_d(v_patch, ns_config->v_patch);
	
	
	
	
	this->v_linears = (vec_d *) br_malloc(ns_config->num_v_linears*sizeof(vec_d));
	for (int ii=0; ii<ns_config->num_v_linears; ii++) {
		init_vec_d(v_linears[ii],ns_config->num_v_vars);
		v_linears[ii]->size = ns_config->num_v_vars;
		vec_mp_to_d(v_linears[ii], ns_config->v_linears[ii]);
	}
	
	num_additional_linears = ns_config->num_additional_linears;
	additional_linears_terminal = (vec_d *) br_malloc(ns_config->num_additional_linears*sizeof(vec_d));
	additional_linears_starting = (vec_d *) br_malloc(ns_config->num_additional_linears*sizeof(vec_d));
	
	for (int ii=0; ii<ns_config->num_additional_linears; ii++) {
		init_vec_d(additional_linears_terminal[ii], ns_config->num_x_vars);
		additional_linears_terminal[ii]->size = ns_config->num_x_vars;
		vec_mp_to_d(additional_linears_terminal[ii],ns_config->additional_linears_terminal[ii]);
		
		init_vec_d(additional_linears_starting[ii], ns_config->num_x_vars);
		additional_linears_starting[ii]->size = ns_config->num_x_vars;
		vec_mp_to_d(additional_linears_starting[ii],ns_config->additional_linears_starting[ii]);
	}
	
	
	
	starting_linears = (vec_d **)br_malloc(ns_config->num_jac_equations*sizeof(vec_d *));
	
	for (int ii=0; ii<ns_config->num_jac_equations; ++ii) {
		starting_linears[ii] = (vec_d *)br_malloc(ns_config->max_degree*sizeof(vec_d));
		for (int jj=0; jj<ns_config->max_degree; jj++) {
			init_vec_d(starting_linears[ii][jj],W.num_variables);
			starting_linears[ii][jj]->size = W.num_variables;
			
			vec_mp_to_d(starting_linears[ii][jj], ns_config->starting_linears[ii][jj]);
		}
	}
	
	if (this->MPType==2)
	{
		this->BED_mp->setup(_SLP,ns_config, W, solve_options);
		rat_to_d(this->gamma, this->BED_mp->gamma_rat);
	}
	
	return SUCCESSFUL;
}

///////////////
//
//   end nullspace_eval_data_d
//
/////////////






















int nullspacejac_solver_master_entry_point(int										MPType,
																					 witness_set						&W, // carries with it the start points, and the linears.
																					 witness_set						*W_new, // new data goes in here
																					 nullspace_config				*ns_config,
																					 solver_configuration		& solve_options)
{

	
	if (solve_options.use_parallel()) {
		solve_options.call_for_help(NULLSPACE);
	}
	
	
	
	
	
	//	solve_options.T.numVars = setupProg(this->SLP, solve_options.T.Precision, solve_options.T.MPType);
	
	int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
	
	prog_t SLP;
	//	// setup a straight-line program, using the file(s) created by the parser
  solve_options.T.numVars = setupProg_count(&SLP, solve_options.T.Precision, solve_options.T.MPType,
																						&startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv,
																						&subFuncsBelow);
	
	nullspacejac_eval_data_d *ED_d = NULL;
	nullspacejac_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new nullspacejac_eval_data_d(0);
			
			ED_d->setup(&SLP,
									ns_config,
									W,
									solve_options);
			break;
			
		case 1:
			ED_mp = new nullspacejac_eval_data_mp(1);
			
			ED_mp->setup(&SLP,
									 ns_config,
									 W,
									 solve_options);
			// initialize latest_newton_residual_mp
			mpf_init(solve_options.T.latest_newton_residual_mp);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new nullspacejac_eval_data_d(2);
			
			ED_mp = ED_d->BED_mp;
			
			
			ED_d->setup(&SLP,
									ns_config,
									W,
									solve_options);
			
			
			
			adjust_tracker_AMP(& (solve_options.T), W.num_variables); // & initialize latest_newton_residual_mp
			
			break;
		default:
			break;
	}

	
	
	generic_solver_master(W_new, W,
												ED_d, ED_mp,
												solve_options);
	
		
  return SUCCESSFUL;

}






void nullspace_slave_entry_point(solver_configuration & solve_options)
{
	
	
	// already received the flag which indicated that this worker is going to be performing the nullspace calculation.
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	nullspacejac_eval_data_d *ED_d = NULL;
	nullspacejac_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new nullspacejac_eval_data_d(0);
			ED_d->receive(solve_options);
			break;
			
		case 1:
			ED_mp = new nullspacejac_eval_data_mp(1);
			ED_mp->receive(solve_options);
			
			std::cout << "worker done receiving mp type" << std::endl;
			// initialize latest_newton_residual_mp
			mpf_init(solve_options.T.latest_newton_residual_mp);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new nullspacejac_eval_data_d(2);
			
			ED_d->receive(solve_options);
			
			ED_mp = ED_d->BED_mp;
			std::cout << "worker done receiving double_mp type" << std::endl;
			
			
			
			// initialize latest_newton_residual_mp
			mpf_init(solve_options.T.latest_newton_residual_mp);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		default:
			break;
	}
	
	

	// call the file setup function
	FILE *OUT = NULL, *midOUT = NULL;
	
	generic_setup_files(&OUT, "output",
											&midOUT, "midpath_data");
	
	trackingStats trackCount; init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	std::cout << "slave entering tracker loop" << std::endl;
	generic_tracker_loop_worker(&trackCount, OUT, midOUT,
															ED_d, ED_mp,
															solve_options);
	
	
	// close the files
	fclose(midOUT);   fclose(OUT);
	
	
	//clear data
}






int nullspacejac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real

	
//	std::cout << "++++++++++++++++++++++++++\n";
	
	
  nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED; // to avoid having to cast every time
	
	
	
  int ii, jj, kk, mm;
	int offset;
  comp_d one_minus_s, gamma_s;
	
	set_one_d(one_minus_s);
  sub_d(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
  mul_d(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	

	vec_d curr_x_vars; init_vec_d(curr_x_vars, BED->num_x_vars);
	curr_x_vars->size = BED->num_x_vars;
	for (ii=0; ii<BED->num_x_vars; ii++) 
		set_d(&curr_x_vars->coord[ii], &current_variable_values->coord[ii]);
	
	vec_d curr_v_vars; init_vec_d(curr_v_vars, BED->num_v_vars);
	curr_v_vars->size = BED->num_v_vars;
	for (ii=0; ii<BED->num_v_vars; ii++)
		set_d(&curr_v_vars->coord[ii], &current_variable_values->coord[ii+BED->num_x_vars]);
	
	vec_d patchValues; init_vec_d(patchValues, 0);
	vec_d temp_function_values; init_vec_d(temp_function_values,0);
	

	vec_d AtimesF;  init_vec_d(AtimesF,0);
	vec_d linprod_x;  init_vec_d(linprod_x, BED->num_jac_equations);
		linprod_x->size = BED->num_jac_equations;
	vec_d linprod_times_gamma_s; init_vec_d(linprod_times_gamma_s,BED->num_jac_equations);
		linprod_times_gamma_s->size = BED->num_jac_equations;
	vec_d tempvec; init_vec_d(tempvec,0);
	vec_d tempvec2; init_vec_d(tempvec2,0);
	
	
	
	mat_d Jv_Patch; init_mat_d(Jv_Patch, 0, 0);
	mat_d tempmat; init_mat_d(tempmat,BED->num_variables-1,BED->num_variables-1);
		tempmat->rows = tempmat->cols = BED->num_variables-1; // change the size indicators

	mat_d lin_func_vals; init_mat_d(lin_func_vals,BED->num_jac_equations, BED->max_degree);
		lin_func_vals->rows = BED->num_jac_equations; lin_func_vals->cols = BED->max_degree;
	
	
	
	
	
	mat_d AtimesJ; init_mat_d(AtimesJ,1,1); AtimesJ->rows = AtimesJ->cols = 1;
	
	mat_d Jv_jac; init_mat_d(Jv_jac,0,0);
	mat_d temp_jacobian_functions, temp_jacobian_parameters;
		init_mat_d(temp_jacobian_functions,0,0); init_mat_d(temp_jacobian_parameters,0,0);
	
	mat_d linprod_derivative_wrt_x;
		init_mat_d(linprod_derivative_wrt_x, BED->num_jac_equations, BED->num_x_vars);
		linprod_derivative_wrt_x->rows = BED->num_jac_equations; linprod_derivative_wrt_x->cols = BED->num_x_vars;
	
	
	comp_d running_prod;
	comp_d temp, temp2, temp3;
	
	
	
	mat_d S_times_Jf_pi;  init_mat_d(S_times_Jf_pi, BED->post_randomizer_matrix->rows, BED->num_v_vars); // set up temp matrix
		S_times_Jf_pi->rows = BED->post_randomizer_matrix->rows; S_times_Jf_pi->cols = BED->num_v_vars;
	vec_d target_function_values;  init_vec_d(target_function_values,0);
	vec_d target_function_values_times_oneminus_s;  init_vec_d(target_function_values_times_oneminus_s,0);
	
	vec_d start_function_values;
	init_vec_d(start_function_values,BED->num_jac_equations); start_function_values->size = BED->num_jac_equations;
	
	//initialize some containers, for the unused stuff from the called evaluators.
	point_d unused_function_values, unused_parVals;
		init_vec_d(unused_function_values,0);init_vec_d(unused_parVals,0);
	vec_d unused_parDer; init_vec_d(unused_parDer,0);
	mat_d unused_Jp; init_mat_d(unused_Jp,0,0);
	mat_d perturbed_Jv; init_mat_d(perturbed_Jv,0,0);
	
	
	mat_d perturbed_AtimesJ, tempmat3; // create matrices
		init_mat_d(perturbed_AtimesJ,1,1); init_mat_d(tempmat3,0,0);
	perturbed_AtimesJ->rows = perturbed_AtimesJ->cols = 1;
	
	mat_d tempmat1,tempmat2; // create temp matrices
		init_mat_d(tempmat1,BED->num_jac_equations,BED->num_v_vars);
		tempmat1->rows = BED->num_jac_equations; tempmat1->cols = BED->num_v_vars;
	
		init_mat_d(tempmat2,BED->num_jac_equations,BED->num_v_vars);
		tempmat2->rows = BED->num_jac_equations; tempmat2->cols = BED->num_v_vars;
	
	mat_d jac_homogenizing_matrix; init_mat_d(jac_homogenizing_matrix,0,0);
	
	
	//initialize the jacobians we will work with.
	point_d perturbed_forward_variables, perturbed_backward_variables;
	init_vec_d(perturbed_forward_variables,0); init_vec_d(perturbed_backward_variables,0);
	change_size_vec_d(perturbed_forward_variables, BED->num_x_vars);
	change_size_vec_d(perturbed_backward_variables,BED->num_x_vars);
	perturbed_backward_variables->size = perturbed_forward_variables->size = BED->num_x_vars;
	
	
	// the main evaluations for $x$
	evalProg_d(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, curr_x_vars, pathVars, BED->SLP);
	
	
  // evaluate the patch
  patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, curr_x_vars, pathVars, &BED->patch);  // Jp is ignored
	

	
	//resize output variables to correct size
	change_size_vec_d(funcVals,BED->num_variables);
  change_size_mat_d(Jv, BED->num_variables, BED->num_variables);
  change_size_mat_d(Jp, BED->num_variables, 1);
	
	
	//////
	// initialize stuff to all 0's
	///////
	
  funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
  Jv->cols = BED->num_variables;  //  <-- this must be square
  Jp->cols = 1;
	
	for (ii=0; ii<Jv->rows; ii++) 
		for (jj=0; jj<Jv->cols; jj++) 
			set_zero_d(&Jv->entry[ii][jj]);
	
	for (ii = 0; ii<BED->num_variables; ii++) 
		set_zero_d(&Jp->entry[ii][0]);  // initialize entire matrix to 0
	
	
	// orig eqns
	
	
	// randomize
	mul_mat_vec_d(AtimesF,BED->randomizer_matrix, temp_function_values); // set values of AtimesF (A is randomization matrix)

	// set func vals
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real) randomization
		set_d(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	
	

	// the jacobian equations for orig
	
	//  randomize the original functions and jacobian

	mat_mul_d(AtimesJ,BED->randomizer_matrix,temp_jacobian_functions);
	
	// copy the jacobian into the return value for the evaluator
	for (ii=0; ii< AtimesJ->rows; ii++) 
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_d(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]); 

	// copy in the transpose of the (randomized) jacobian, omitting the homogenizing variable
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=1; jj<BED->num_x_vars; jj++)
			set_d(&BED->jac_with_proj->entry[jj-1][ii], &AtimesJ->entry[ii][jj]);

	
	
	// the additional linears.  there are $r-\ell$ of them.
	offset = BED->num_randomized_eqns;
	for (ii=0; ii< BED->num_additional_linears; ii++) {
		dot_product_d(temp, BED->additional_linears_terminal[ii], curr_x_vars);
		mul_d(temp3, temp, one_minus_s);
		neg_d(&Jp->entry[offset+ii][0], temp);  // Jp = -terminal
		
		dot_product_d(temp,  BED->additional_linears_starting[ii], curr_x_vars);
		
		mul_d(temp2, temp, BED->gamma);
		add_d(&Jp->entry[offset+ii][0], &Jp->entry[offset+ii][0], temp2);   // Jp = -terminal + gamma*start
		
		mul_d(temp2, temp, gamma_s);
		
		add_d(&funcVals->coord[offset+ii],temp2, temp3); // (gamma*s)*start(x) + (1-s)*terminal(x)
		
		for (jj=0; jj<BED->num_x_vars; jj++) {
			mul_d(temp, gamma_s,    &BED->additional_linears_starting[ii]->coord[jj]); 
			mul_d(temp2,one_minus_s,&BED->additional_linears_terminal[ii]->coord[jj]);
			
			add_d(&Jv->entry[ii+offset][jj], temp, temp2);
		}
	} // √
	
	
	
	// NOW WE WILL WORK ON THE TARGET SYSTEM'S FUNCTION VALUES
	
	
	
	// make the homogenizing matrix for the $x$ variables
	make_matrix_ID_d(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
	
	for (ii=0; ii<BED->num_randomized_eqns; ii++)
		for (jj=0; jj<(BED->max_degree - (BED->randomized_degrees[ii]-1)); jj++)
			mul_d(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
	
	for (ii=BED->num_randomized_eqns; ii<BED->num_v_vars; ii++)
		for (jj=0; jj<(BED->max_degree); jj++) // these are all degree 1
			mul_d(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
	

	mat_mul_d(tempmat, BED->post_randomizer_matrix, BED->jac_with_proj); // jac with proj having been set above with unperturbed values
	mat_mul_d(S_times_Jf_pi, tempmat, jac_homogenizing_matrix); // jac with proj having been set above with unperturbed values


	mul_mat_vec_d(target_function_values, S_times_Jf_pi, curr_v_vars);
	vec_mulcomp_d(target_function_values_times_oneminus_s, target_function_values, one_minus_s);
	
	
	
	//  THE LINPROD START SYSTEM FUNCTION VALUES
	offset = BED->num_randomized_eqns + BED->num_additional_linears;  
	// the product of the linears
	for (jj=0; jj<BED->num_jac_equations; jj++) {
		
		//perform the $x$ evaluation
		set_one_d(&linprod_x->coord[jj]); // initialize to 1 for multiplication
		for (ii=0; ii<BED->max_degree; ++ii) {
			dot_product_d(&lin_func_vals->entry[jj][ii],BED->starting_linears[jj][ii],curr_x_vars); // save into a buffer for calculating the derivative later.
			mul_d(&linprod_x->coord[jj], &linprod_x->coord[jj], &lin_func_vals->entry[jj][ii]);// multiply linprod_x times the value we just created
		}
		
		//perform the $v$ evaluation
		dot_product_d(temp, BED->v_linears[jj], curr_v_vars);
		
		//now set the combined $x,v$ value
		mul_d(&start_function_values->coord[jj], &linprod_x->coord[jj], temp); // start_function_values = linprod_x(x)*linear(v)
		mul_d(&linprod_times_gamma_s->coord[jj], gamma_s, &start_function_values->coord[jj]); // sets the value linprod_x*gamma*s
		
		
		neg_d( &Jp->entry[jj + offset][0], &target_function_values->coord[jj]);  // temp = -target
		mul_d( temp, BED->gamma, &start_function_values->coord[jj]);  // temp = gamma*linprod_x(x)*linear(v)
		add_d( &Jp->entry[jj + offset][0], &Jp->entry[jj + offset][0], temp); // Jp = -target + gamma*start
	}
	
//	print_point_to_screen_matlab(target_function_values,"target_f");
//	print_point_to_screen_matlab(target_function_values_times_oneminus_s,"target_f_one_minus_s");
//	print_point_to_screen_matlab(start_function_values,"start_f");
//	print_point_to_screen_matlab(linprod_times_gamma_s,"start_f_gamma_s");
//	print_comp_matlab(gamma_s,"gamma_s");
	
	
	offset = BED->num_randomized_eqns + BED->num_additional_linears; // N-r+k
	for (ii=0; ii<BED->num_jac_equations; ii++) {
		add_d(&funcVals->coord[ii+offset], &target_function_values_times_oneminus_s->coord[ii], &linprod_times_gamma_s->coord[ii]);
	}
	
	
	
	
	
	
	
	
	// DERIVATIVES OF THE HOMOTOPY WRT V
	offset = BED->num_randomized_eqns + BED->num_additional_linears; // N-k+l
	for (ii=0; ii<BED->num_jac_equations; ii++) {
		for (jj=0; jj<BED->num_v_vars; jj++) {
			mul_d(temp, &BED->v_linears[ii]->coord[jj], &linprod_x->coord[ii]);  // temp = M_ij * linprod_x(x)
			mul_d(temp2, temp, gamma_s);                                      // temp2 = gamma*s*M_ij * linprod_x(x)
			
//			std::stringstream converter;
//			converter << "linprod_v_der_" << ii << "_" << jj;
//			print_comp_matlab(temp, converter.str());
//			converter.str("");
			
			mul_d(temp, one_minus_s, &S_times_Jf_pi->entry[ii][jj]);           // temp = (1-s)*(S*[Jf^T pi^T])_ij
			
			add_d(&Jv->entry[ii+offset][jj+BED->num_x_vars], temp, temp2);    // Jv = temp + temp2
		}
	}
	// √ for sphere

//	print_matrix_to_screen_matlab(S_times_Jf_pi,"S_times_Jf_pi");
	
	// now the x derivatives corresponding to the linprod start system
	
	// an implementation of the product rule
	for (mm=0; mm< BED->num_jac_equations; mm++) {
		for (kk=0; kk<BED->num_x_vars; kk++) { // for each variable
			set_zero_d(&linprod_derivative_wrt_x->entry[mm][kk]); // initialize to 0 for the sum
			
			for (ii=0; ii<BED->max_degree; ++ii) { //  for each linear
				
				set_d(running_prod, &BED->starting_linears[mm][ii]->coord[kk]);// initialize the product
				for (jj=0; jj<BED->max_degree; jj++) {
					if (jj!=ii) {
						mul_d(running_prod,running_prod,&lin_func_vals->entry[mm][jj]); // the linear evaluated at curr_var_vals
					}
				}//re: jj
				add_d(&linprod_derivative_wrt_x->entry[mm][kk],&linprod_derivative_wrt_x->entry[mm][kk],running_prod);
				
			}// re:ii
			
			dot_product_d(temp, BED->v_linears[mm], curr_v_vars);  // these two lines multiply by  (v_linear •	v)
			mul_d(&linprod_derivative_wrt_x->entry[mm][kk], &linprod_derivative_wrt_x->entry[mm][kk], temp);
		} // re: kk
	} // re: mm
	
	
	
	
	
	// NUMERICALLY DIFFERENTIATE THE derivative of the target jacobian system wrt $x$.
	
	offset = BED->num_randomized_eqns + BED->num_additional_linears;
	for (ii=0; ii<BED->num_x_vars; ii++) {
		
		//go forward
		vec_cp_d(perturbed_forward_variables, curr_x_vars);
		add_d( &perturbed_forward_variables->coord[ii], &perturbed_forward_variables->coord[ii],BED->perturbation);
		
		
		make_matrix_ID_d(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
		for (kk=0; kk<BED->num_randomized_eqns; kk++)
			for (jj=0; jj<(BED->max_degree - (BED->randomized_degrees[kk]-1)); jj++) // -1 because differentiation
				mul_d(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_forward_variables->coord[0]);
		
		for (kk=BED->num_randomized_eqns; kk<BED->num_v_vars; kk++)
			for (jj=0; jj<(BED->max_degree); jj++)
				mul_d(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_forward_variables->coord[0]);
		
		// evaluate perturbed forwards
		evalProg_d(unused_function_values, unused_parVals, unused_parDer,  //  unused output
							 perturbed_Jv,  // <---- the output we need
							 unused_Jp, //unused output
							 perturbed_forward_variables, pathVars, BED->SLP); // input
		
		
		mat_mul_d(perturbed_AtimesJ,BED->randomizer_matrix,perturbed_Jv);
		
		mat_cp_d(tempmat1, BED->jac_with_proj);
		// copy in the transpose of the (randomized) jacobian
		for (mm=0; mm< (BED->num_randomized_eqns); mm++) {
			for (jj=1; jj<BED->num_x_vars; jj++) {
				set_d(&tempmat1->entry[jj - 1][mm],&perturbed_AtimesJ->entry[mm][jj]); 
			}
		}

		mat_mul_d(tempmat2, BED->post_randomizer_matrix, tempmat1);
		mat_mul_d(tempmat3, tempmat2, jac_homogenizing_matrix);
		mul_mat_vec_d(tempvec, tempmat3, curr_v_vars);
		
		
		
		//go backward
		
		vec_cp_d(perturbed_backward_variables, curr_x_vars);
		sub_d(&perturbed_backward_variables->coord[ii], &perturbed_backward_variables->coord[ii],BED->perturbation);
		
		make_matrix_ID_d(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
		for (kk=0; kk<BED->num_randomized_eqns; kk++)
			for (jj=0; jj<(BED->max_degree - (BED->randomized_degrees[kk]-1)); jj++)
				mul_d(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_backward_variables->coord[0]);
		for (kk=BED->num_randomized_eqns; kk<BED->num_v_vars; kk++)
			for (jj=0; jj<(BED->max_degree); jj++)
				mul_d(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_backward_variables->coord[0]);
		

		// evaluate perturbed back
		evalProg_d(unused_function_values, unused_parVals, unused_parDer,  //  unused output
							 perturbed_Jv,  // <---- the output we need
							 unused_Jp, //unused output
							 perturbed_backward_variables, pathVars, BED->SLP); // input


		mat_cp_d(tempmat1, BED->jac_with_proj); // probably unnecessary
		
		// copy in the transpose of the (randomized) jacobian
		mat_mul_d(perturbed_AtimesJ,BED->randomizer_matrix,perturbed_Jv);
		for (mm=0; mm< (BED->num_randomized_eqns); mm++) {
			for (jj=1; jj<BED->num_x_vars; jj++) {
				set_d(&tempmat1->entry[jj - 1][mm], &perturbed_AtimesJ->entry[mm][jj]);
			}
		}

		mat_mul_d(tempmat2, BED->post_randomizer_matrix, tempmat1);
		mat_mul_d(tempmat3, tempmat2, jac_homogenizing_matrix);
		mul_mat_vec_d(tempvec2, tempmat3, curr_v_vars);
		
		
		vec_sub_d(tempvec, tempvec, tempvec2); // tempvec = forward - backward

		div_d(temp, BED->half, BED->perturbation);   //this is repetitively wasteful
		
		vec_mulcomp_d(tempvec, tempvec, temp); //tempvec = (forward-backward)/(2h)
		// √ this is verified correct for sphere
		
//		print_point_to_screen_matlab(tempvec,"Jv_target");
		
		vec_mulcomp_d(tempvec, tempvec, one_minus_s);// now tempvec has the numerical derivatives wrt $x$ variable ii. (for all jac eqns)
		
		// now, combine this numerical derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
		
		for (mm=0; mm<BED->num_jac_equations; mm++) {
			mul_d(temp, gamma_s, &linprod_derivative_wrt_x->entry[mm][ii]);
			add_d(&Jv->entry[mm+offset][ii], &tempvec->coord[mm], temp);
//			std::cout << "setting Jv[" << mm+offset << "][" << ii << "]\n";
		}
		
	}//re: ii for numerical diff

	// END NUMERICAL DIFF wrt x
	
	
	
	//set the X PATCH values
	offset = BED->num_randomized_eqns + BED->num_additional_linears + BED->num_jac_equations;
	
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_d(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
		
		// Jv = Jv_Patch
		for (jj = 0; jj<BED->num_x_vars; jj++) // for each variable
			set_d(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	

	offset = BED->num_randomized_eqns + BED->num_additional_linears + BED->num_jac_equations + BED->patch.num_patches;
	if (offset != BED->num_variables-1) {
		std::cout << "mismatch in number of blabla, line 1006;\n" << offset << " " << BED->num_variables-1 << std::endl;
		print_matrix_to_screen_matlab(Jv,"Jv");
		deliberate_segfault();
	}
	
	// V patch
	set_one_d(temp2);
	dot_product_d(temp, BED->v_patch, curr_v_vars);
	sub_d(&funcVals->coord[BED->num_variables-1], temp, temp2);  // f = patch*v-1

	for (ii=0; ii<BED->num_v_vars; ii++) 
		set_d(&Jv->entry[BED->num_variables-1][BED->num_x_vars+ii], &BED->v_patch->coord[ii]);
	
	
	
	// finally, set parVals & parDer correctly
	
  change_size_point_d(parVals, 1);  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
	
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	if (BED->verbose_level==10) {
	printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
//	print_matrix_to_screen_matlab(jac_homogenizing_matrix,"jac_hom_1044");
//	print_matrix_to_screen_matlab(BED->post_randomizer_matrix,"S");
//	print_matrix_to_screen_matlab(BED->randomizer_matrix,"R");
//
//	
//	print_matrix_to_screen_matlab( AtimesJ,"jac");
//	print_point_to_screen_matlab(curr_x_vars,"currxvars");
//	print_point_to_screen_matlab(current_variable_values,"curr_vars");
	print_point_to_screen_matlab(funcVals,"F");
	print_matrix_to_screen_matlab(Jv,"Jv");
//	print_matrix_to_screen_matlab(Jp,"Jp");
//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jacwithproj");
//			//these values are set in this function:  point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp
//	print_matrix_to_screen_matlab(BED->randomizer_matrix,"randomizer_matrix");

//	std::cout << "\n\n**************\n\n";
//	mypause();
		if (BED->verbose_level==10)
			mypause();
	}
	
	
	
	clear_vec_d(curr_x_vars);
	clear_vec_d(curr_v_vars);
	clear_vec_d(patchValues);
	clear_vec_d(temp_function_values);

	
	clear_vec_d(AtimesF);
	clear_vec_d(linprod_x);
	clear_vec_d(linprod_times_gamma_s);
	clear_vec_d(tempvec);
	clear_vec_d(tempvec2);
	

	clear_mat_d(Jv_Patch);
	clear_mat_d(tempmat);
	clear_mat_d(lin_func_vals);
	

	
	clear_mat_d(AtimesJ);
	clear_mat_d(Jv_jac);
	clear_mat_d(temp_jacobian_functions);
	clear_mat_d(temp_jacobian_parameters);
	clear_mat_d(linprod_derivative_wrt_x);
	
	clear_vec_d(start_function_values);
	clear_vec_d(target_function_values);
	clear_vec_d(target_function_values_times_oneminus_s);
	clear_mat_d(S_times_Jf_pi);
	
	
	clear_vec_d(unused_function_values);
	clear_vec_d(unused_parVals);
	clear_vec_d(unused_parDer);

	
	clear_mat_d(unused_Jp);
	clear_mat_d(perturbed_Jv);
	clear_mat_d(perturbed_AtimesJ);
	clear_mat_d(tempmat3);
	clear_mat_d(tempmat1);
	clear_mat_d(tempmat2);

	clear_vec_d(perturbed_forward_variables);
	clear_vec_d(perturbed_backward_variables);

	
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_d dehommed; init_vec_d(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,curr_x_vars);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_d(dehommed);
#endif
	
	//	printf("exiting eval\n");
  return 0;
}




int nullspacejac_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	//	printf("entering eval_mp\n");
  nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED; // to avoid having to cast every time
	
  int ii, jj, kk, mm;
	int offset;
  comp_mp one_minus_s, gamma_s;  init_mp(one_minus_s); init_mp(gamma_s);
	
	set_one_mp(one_minus_s);
  sub_mp(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
  mul_mp(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
  // we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	vec_mp curr_x_vars; init_vec_mp(curr_x_vars, BED->num_x_vars);
	curr_x_vars->size = BED->num_x_vars;
	for (ii=0; ii<BED->num_x_vars; ii++)
		set_mp(&curr_x_vars->coord[ii], &current_variable_values->coord[ii]);
	
	vec_mp curr_v_vars; init_vec_mp(curr_v_vars, BED->num_v_vars);
	curr_v_vars->size = BED->num_v_vars;
	for (ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&curr_v_vars->coord[ii], &current_variable_values->coord[ii+BED->num_x_vars]);
	
	
	
	
	
	
	
	
	vec_mp patchValues; init_vec_mp(patchValues, 0);
	vec_mp temp_function_values; init_vec_mp(temp_function_values,0);
	
	
	vec_mp AtimesF;  init_vec_mp(AtimesF,0);
	vec_mp linprod_x;  init_vec_mp(linprod_x, BED->num_jac_equations);
	linprod_x->size = BED->num_jac_equations;
	vec_mp linprod_times_gamma_s; init_vec_mp(linprod_times_gamma_s,BED->num_jac_equations);
	linprod_times_gamma_s->size = BED->num_jac_equations;
	vec_mp tempvec; init_vec_mp(tempvec,0);
	vec_mp tempvec2; init_vec_mp(tempvec2,0);
	
	
	
	mat_mp Jv_Patch; init_mat_mp(Jv_Patch, 0, 0);
	mat_mp tempmat; init_mat_mp(tempmat,BED->num_variables-1,BED->num_variables-1);
	tempmat->rows = tempmat->cols = BED->num_variables-1; // change the size indicators
	
	mat_mp lin_func_vals; init_mat_mp(lin_func_vals,BED->num_jac_equations, BED->max_degree);
	lin_func_vals->rows = BED->num_jac_equations; lin_func_vals->cols = BED->max_degree;
	
	
	
	
	
	mat_mp AtimesJ; init_mat_mp(AtimesJ,1,1);
	mat_mp Jv_jac; init_mat_mp(Jv_jac,0,0);
	mat_mp temp_jacobian_functions, temp_jacobian_parameters;
	init_mat_mp(temp_jacobian_functions,0,0); init_mat_mp(temp_jacobian_parameters,0,0);
	
	mat_mp linprod_derivative_wrt_x;
	init_mat_mp(linprod_derivative_wrt_x, BED->num_jac_equations, BED->num_x_vars);
	linprod_derivative_wrt_x->rows = BED->num_jac_equations; linprod_derivative_wrt_x->cols = BED->num_x_vars;
	
	
	comp_mp running_prod;  init_mp(running_prod);
	comp_mp temp, temp2, temp3; init_mp(temp); init_mp(temp2); init_mp(temp3);
	
	
	
	mat_mp S_times_Jf_pi;  init_mat_mp(S_times_Jf_pi, BED->post_randomizer_matrix->rows, BED->num_v_vars); // set up temp matrix
	S_times_Jf_pi->rows = BED->post_randomizer_matrix->rows; S_times_Jf_pi->cols = BED->num_v_vars;
	vec_mp target_function_values;  init_vec_mp(target_function_values,0);
	vec_mp target_function_values_times_oneminus_s;  init_vec_mp(target_function_values_times_oneminus_s,0);
	
	vec_mp start_function_values;
	init_vec_mp(start_function_values,BED->num_jac_equations); start_function_values->size = BED->num_jac_equations;
	
	//initialize some containers, for the unused stuff from the called evaluators.
	point_mp unused_function_values, unused_parVals;
	init_vec_mp(unused_function_values,0);init_vec_mp(unused_parVals,0);
	vec_mp unused_parDer; init_vec_mp(unused_parDer,0);
	mat_mp unused_Jp; init_mat_mp(unused_Jp,0,0);
	mat_mp perturbed_Jv; init_mat_mp(perturbed_Jv,0,0);
	
	
	mat_mp perturbed_AtimesJ, tempmat3; // create matrices
	init_mat_mp(perturbed_AtimesJ,0,0); init_mat_mp(tempmat3,0,0);
	
	
	mat_mp tempmat1,tempmat2; // create temp matrices
	init_mat_mp(tempmat1,BED->num_jac_equations,BED->num_v_vars);
	tempmat1->rows = BED->num_jac_equations; tempmat1->cols = BED->num_v_vars;
	
	init_mat_mp(tempmat2,BED->num_jac_equations,BED->num_v_vars);
	tempmat2->rows = BED->num_jac_equations; tempmat2->cols = BED->num_v_vars;
	
	mat_mp jac_homogenizing_matrix; init_mat_mp(jac_homogenizing_matrix,0,0);
	
	//initialize the jacobians we will work with.
	point_mp perturbed_forward_variables, perturbed_backward_variables;
	init_vec_mp(perturbed_forward_variables,0); init_vec_mp(perturbed_backward_variables,0);
	change_size_vec_mp(perturbed_forward_variables,BED->num_x_vars);
	change_size_vec_mp(perturbed_backward_variables,BED->num_x_vars);
	perturbed_backward_variables->size = perturbed_forward_variables->size = BED->num_x_vars;
	
	
	// the main evaluations for $x$
	evalProg_mp(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, curr_x_vars, pathVars, BED->SLP);
	
	
  patch_eval_mp(    patchValues, parVals, parDer, Jv_Patch, Jp, curr_x_vars, pathVars, &BED->patch);  // Jp is ignored
	
	
	
	
	//resize output variables to correct size
	change_size_vec_mp(funcVals,BED->num_variables);
  change_size_mat_mp(Jv, BED->num_variables, BED->num_variables);
  change_size_mat_mp(Jp, BED->num_variables, 1);
	
	// initialize stuff to all 0's
  funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
  Jv->cols = BED->num_variables;  //  <-- this must be square
  Jp->cols = 1;
	
	for (ii=0; ii<Jv->rows; ii++)
		for (jj=0; jj<Jv->cols; jj++)
			set_zero_mp(&Jv->entry[ii][jj]);
	
	for (ii = 0; ii<BED->num_variables; ii++)
		set_zero_mp(&Jp->entry[ii][0]);  // initialize entire matrix to 0
	
	// orig eqns
	
	
	// randomize
	mul_mat_vec_mp(AtimesF, BED->randomizer_matrix, temp_function_values); // set values of AtimesF (A is randomization matrix)
	
	// set func vals
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real) randomization
		set_mp(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	
	
	
	// the jacobian equations for orig
	
	//  randomize the original functions and jacobian
	mat_mul_mp(AtimesJ,BED->randomizer_matrix,temp_jacobian_functions);
	
	// copy the jacobian into the return value for the evaluator
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_mp(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
	
	// copy in the transpose of the (randomized) jacobian, omitting the homogenizing variable
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=1; jj<BED->num_x_vars; jj++)
			set_mp(&BED->jac_with_proj->entry[jj-1][ii], &AtimesJ->entry[ii][jj]);
	//√
	
	
	// the additional linears.  there are $r-\ell$ of them.
	offset = BED->num_randomized_eqns;
	for (ii=0; ii< BED->num_additional_linears; ii++) {
		dot_product_mp(temp, BED->additional_linears_terminal[ii], curr_x_vars);
		mul_mp(temp3, temp, one_minus_s);
		neg_mp(&Jp->entry[offset+ii][0], temp);  // Jp = -terminal
		
		dot_product_mp(temp,  BED->additional_linears_starting[ii], curr_x_vars);
		
		mul_mp(temp2, temp, BED->gamma);
		add_mp(&Jp->entry[offset+ii][0], &Jp->entry[offset+ii][0], temp2);   // Jp = -terminal + gamma*start
		
		mul_mp(temp2, temp, gamma_s);
		
		add_mp(&funcVals->coord[offset+ii],temp2, temp3); // (gamma*s)*start(x) + (1-s)*terminal(x)
		
		for (jj=0; jj<BED->num_x_vars; jj++) {
			mul_mp(temp, gamma_s,    &BED->additional_linears_starting[ii]->coord[jj]);
			mul_mp(temp2,one_minus_s,&BED->additional_linears_terminal[ii]->coord[jj]);
			
			add_mp(&Jv->entry[ii+offset][jj], temp, temp2);
		}
	} // √
	
	
	
	// NOW WE WILL WORK ON THE TARGET SYSTEM'S FUNCTION VALUES
	
	
	
	// make the homogenizing matrix for the $x$ variables
	make_matrix_ID_mp(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
	
	for (ii=0; ii<BED->num_randomized_eqns; ii++)
		for (jj=0; jj<(BED->max_degree - (BED->randomized_degrees[ii]-1)); jj++)
			mul_mp(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
	
	for (ii=BED->num_randomized_eqns; ii<BED->num_v_vars; ii++)
		for (jj=0; jj<(BED->max_degree); jj++) // these are all degree 1
			mul_mp(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
	
	
	mat_mul_mp(tempmat, BED->post_randomizer_matrix, BED->jac_with_proj); // jac with proj having been set above with unperturbed values
	mat_mul_mp(S_times_Jf_pi, tempmat, jac_homogenizing_matrix); // jac with proj having been set above with unperturbed values
	
	
	mul_mat_vec_mp(target_function_values, S_times_Jf_pi, curr_v_vars);
	vec_mulcomp_mp(target_function_values_times_oneminus_s, target_function_values, one_minus_s);
	
	
	
	//  THE LINPROD START SYSTEM FUNCTION VALUES
	offset = BED->num_randomized_eqns + BED->num_additional_linears;
	// the product of the linears
	for (jj=0; jj<BED->num_jac_equations; jj++) {
		
		//perform the $x$ evaluation
		set_one_mp(&linprod_x->coord[jj]); // initialize to 1 for multiplication
		for (ii=0; ii<BED->max_degree; ++ii) {
			dot_product_mp(&lin_func_vals->entry[jj][ii],BED->starting_linears[jj][ii],curr_x_vars); // save into a buffer for calculating the derivative later.
			mul_mp(&linprod_x->coord[jj], &linprod_x->coord[jj], &lin_func_vals->entry[jj][ii]);// multiply linprod_x times the value we just created
		}
		
		//perform the $v$ evaluation
		dot_product_mp(temp, BED->v_linears[jj], curr_v_vars);
		
		//now set the combined $x,v$ value
		mul_mp(&start_function_values->coord[jj], &linprod_x->coord[jj], temp); // start_function_values = linprod_x(x)*linear(v)
		mul_mp(&linprod_times_gamma_s->coord[jj], gamma_s, &start_function_values->coord[jj]); // sets the value linprod_x*gamma*s
		
		
		neg_mp( &Jp->entry[jj + offset][0], &target_function_values->coord[jj]);  // temp = -target
		mul_mp( temp, BED->gamma, &start_function_values->coord[jj]);  // temp = gamma*linprod_x(x)*linear(v)
		add_mp( &Jp->entry[jj + offset][0], &Jp->entry[jj + offset][0], temp); // Jp = -target + gamma*start
	}
	
	//	print_point_to_screen_matlab(target_function_values,"target_f");
	//	print_point_to_screen_matlab(target_function_values_times_oneminus_s,"target_f_one_minus_s");
	//	print_point_to_screen_matlab(start_function_values,"start_f");
	//	print_point_to_screen_matlab(linprod_times_gamma_s,"start_f_gamma_s");
	//	print_comp_matlab(gamma_s,"gamma_s");
	
	
	offset = BED->num_randomized_eqns + BED->num_additional_linears; // N-r+k
	for (ii=0; ii<BED->num_jac_equations; ii++) {
		add_mp(&funcVals->coord[ii+offset], &target_function_values_times_oneminus_s->coord[ii], &linprod_times_gamma_s->coord[ii]);
	}
	
	
	
	
	
	
	
	
	// DERIVATIVES OF THE HOMOTOPY WRT V
	offset = BED->num_randomized_eqns + BED->num_additional_linears; // N-k+l
	for (ii=0; ii<BED->num_jac_equations; ii++) {
		for (jj=0; jj<BED->num_v_vars; jj++) {
			mul_mp(temp, &BED->v_linears[ii]->coord[jj], &linprod_x->coord[ii]);  // temp = M_ij * linprod_x(x)
			mul_mp(temp2, temp, gamma_s);                                      // temp2 = gamma*s*M_ij * linprod_x(x)
			
			//			std::stringstream converter;
			//			converter << "linprod_v_der_" << ii << "_" << jj;
			//			print_comp_matlab(temp, converter.str());
			//			converter.str("");
			
			mul_mp(temp, one_minus_s, &S_times_Jf_pi->entry[ii][jj]);           // temp = (1-s)*(S*[Jf^T pi^T])_ij
			
			add_mp(&Jv->entry[ii+offset][jj+BED->num_x_vars], temp, temp2);    // Jv = temp + temp2
		}
	}
	// √ for sphere
	
	//	print_matrix_to_screen_matlab(S_times_Jf_pi,"S_times_Jf_pi");
	
	// now the x derivatives corresponding to the linprod start system
	
	// an implementation of the product rule
	for (mm=0; mm< BED->num_jac_equations; mm++) {
		for (kk=0; kk<BED->num_x_vars; kk++) { // for each variable
			set_zero_mp(&linprod_derivative_wrt_x->entry[mm][kk]); // initialize to 0 for the sum
			
			for (ii=0; ii<BED->max_degree; ++ii) { //  for each linear
				
				set_mp(running_prod, &BED->starting_linears[mm][ii]->coord[kk]);// initialize the product
				for (jj=0; jj<BED->max_degree; jj++) {
					if (jj!=ii) {
						mul_mp(running_prod,running_prod,&lin_func_vals->entry[mm][jj]); // the linear evaluated at curr_var_vals
					}
				}//re: jj
				add_mp(&linprod_derivative_wrt_x->entry[mm][kk],&linprod_derivative_wrt_x->entry[mm][kk],running_prod);
				
			}// re:ii
			
			dot_product_mp(temp, BED->v_linears[mm], curr_v_vars);  // these two lines multiply by  (v_linear •	v)
			mul_mp(&linprod_derivative_wrt_x->entry[mm][kk], &linprod_derivative_wrt_x->entry[mm][kk], temp);
		} // re: kk
	} // re: mm
	
	
	
	
	
	// NUMERICALLY DIFFERENTIATE THE derivative of the target jacobian system wrt $x$.
	
	offset = BED->num_randomized_eqns + BED->num_additional_linears;
	for (ii=0; ii<BED->num_x_vars; ii++) {
		
		//go forward
		vec_cp_mp(perturbed_forward_variables, curr_x_vars);
		add_mp( &perturbed_forward_variables->coord[ii], &perturbed_forward_variables->coord[ii],BED->perturbation);
		
		
		make_matrix_ID_mp(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
		for (kk=0; kk<BED->num_randomized_eqns; kk++)
			for (jj=0; jj<(BED->max_degree - (BED->randomized_degrees[kk]-1)); jj++) // -1 because differentiation
				mul_mp(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_forward_variables->coord[0]);
		
		for (kk=BED->num_randomized_eqns; kk<BED->num_v_vars; kk++)
			for (jj=0; jj<(BED->max_degree); jj++)
				mul_mp(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_forward_variables->coord[0]);
		
		// evaluate perturbed forwards
		evalProg_mp(unused_function_values, unused_parVals, unused_parDer,  //  unused output
								perturbed_Jv,  // <---- the output we need
								unused_Jp, //unused output
								perturbed_forward_variables, pathVars, BED->SLP); // input
		
		
		mat_mul_mp(perturbed_AtimesJ,BED->randomizer_matrix,perturbed_Jv);
		
		mat_cp_mp(tempmat1, BED->jac_with_proj);
		// copy in the transpose of the (randomized) jacobian
		for (mm=0; mm< (BED->num_randomized_eqns); mm++) {
			for (jj=1; jj<BED->num_x_vars; jj++) {
				set_mp(&tempmat1->entry[jj - 1][mm],&perturbed_AtimesJ->entry[mm][jj]);
			}
		}
		
		mat_mul_mp(tempmat2, BED->post_randomizer_matrix, tempmat1);
		mat_mul_mp(tempmat3, tempmat2, jac_homogenizing_matrix);
		mul_mat_vec_mp(tempvec, tempmat3, curr_v_vars);
		
		
		
		//go backward
		
		vec_cp_mp(perturbed_backward_variables, curr_x_vars);
		sub_mp(&perturbed_backward_variables->coord[ii], &perturbed_backward_variables->coord[ii],BED->perturbation);
		
		make_matrix_ID_mp(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
		for (kk=0; kk<BED->num_randomized_eqns; kk++)
			for (jj=0; jj<(BED->max_degree - (BED->randomized_degrees[kk]-1)); jj++)
				mul_mp(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_backward_variables->coord[0]);
		for (kk=BED->num_randomized_eqns; kk<BED->num_v_vars; kk++)
			for (jj=0; jj<(BED->max_degree); jj++)
				mul_mp(&jac_homogenizing_matrix->entry[kk][kk], &jac_homogenizing_matrix->entry[kk][kk], &perturbed_backward_variables->coord[0]);
		
		
		// evaluate perturbed back
		evalProg_mp(unused_function_values, unused_parVals, unused_parDer,  //  unused output
								perturbed_Jv,  // <---- the output we need
								unused_Jp, //unused output
								perturbed_backward_variables, pathVars, BED->SLP); // input
		
		
		mat_cp_mp(tempmat1, BED->jac_with_proj); // probably unnecessary
		
		// copy in the transpose of the (randomized) jacobian
		mat_mul_mp(perturbed_AtimesJ,BED->randomizer_matrix,perturbed_Jv);
		for (mm=0; mm< (BED->num_randomized_eqns); mm++) {
			for (jj=1; jj<BED->num_x_vars; jj++) {
				set_mp(&tempmat1->entry[jj - 1][mm], &perturbed_AtimesJ->entry[mm][jj]);
			}
		}
		
		mat_mul_mp(tempmat2, BED->post_randomizer_matrix, tempmat1);
		mat_mul_mp(tempmat3, tempmat2, jac_homogenizing_matrix);
		mul_mat_vec_mp(tempvec2, tempmat3, curr_v_vars);
		
		
		vec_sub_mp(tempvec, tempvec, tempvec2); // tempvec = forward - backward
		
		div_mp(temp, BED->half, BED->perturbation);   //this is repetitively wasteful
		
		
		
		vec_mulcomp_mp(tempvec, tempvec, temp); //tempvec = (forward-backward)/(2h)
																						// √ this is verified correct for sphere
		
		//		print_point_to_screen_matlab(tempvec,"Jv_target");
		
		vec_mulcomp_mp(tempvec, tempvec, one_minus_s);// now tempvec has the numerical derivatives wrt $x$ variable ii. (for all jac eqns)
		
		// now, combine this numerical derivative and the linprod derivative to get the exact values for the $x$ portion of the jacobian matrix to return
		
		for (mm=0; mm<BED->num_jac_equations; mm++) {
			mul_mp(temp, gamma_s, &linprod_derivative_wrt_x->entry[mm][ii]);
			add_mp(&Jv->entry[mm+offset][ii], &tempvec->coord[mm], temp);
			//			std::cout << "setting Jv[" << mm+offset << "][" << ii << "]\n";
		}
		
	}//re: ii for numerical diff
	
	// END NUMERICAL DIFF wrt x
	
	
	
	//set the X PATCH values
	offset = BED->num_randomized_eqns + BED->num_additional_linears + BED->num_jac_equations;
	
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_mp(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
		
		// Jv = Jv_Patch
		for (jj = 0; jj<BED->num_x_vars; jj++) // for each variable
			set_mp(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	
	offset = BED->num_randomized_eqns + BED->num_additional_linears + BED->num_jac_equations + BED->patch.num_patches;
	if (offset != BED->num_variables-1) {
		std::cout << "mismatch in number of blabla, line 2701;\n" << offset << " " << BED->num_variables-1 << std::endl;
		print_matrix_to_screen_matlab(Jv,"Jv");
		deliberate_segfault();
	}
	
	// V patch
	set_one_mp(temp2);
	dot_product_mp(temp, BED->v_patch, curr_v_vars);
	sub_mp(&funcVals->coord[BED->num_variables-1], temp, temp2);  // f = patch*v-1
	
	for (ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&Jv->entry[BED->num_variables-1][BED->num_x_vars+ii], &BED->v_patch->coord[ii]);
	
	
	// finally, set parVals & parDer correctly
	
  change_size_point_mp(parVals, 1);  change_size_vec_mp(parDer, 1);
  parVals->size = parDer->size = 1;
	
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	
	if (BED->verbose_level==10) {
		//	print_matrix_to_screen_matlab( AtimesJ,"jac");
		//	print_point_to_screen_matlab(curr_x_vars,"currxvars");
		print_point_to_screen_matlab(funcVals,"F_mp");
		//	print_point_to_screen_matlab(parVals,"parVals");
		//	print_point_to_screen_matlab(parDer,"parDer");
		print_matrix_to_screen_matlab(Jv,"Jv_mp");
		//	print_matrix_to_screen_matlab(Jp,"Jp");
		
		//	print_matrix_to_screen_matlab(BED->jac_with_proj,"jacwithproj");
		//these values are set in this function:  point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp
		//	print_matrix_to_screen_matlab(BED->randomizer_matrix,"randomizer_matrix");
		
		if (BED->verbose_level==10)
			mypause();
	}
	
	
	
	
	clear_vec_mp(curr_x_vars);
	clear_vec_mp(curr_v_vars);
	clear_vec_mp(patchValues);
	clear_vec_mp(temp_function_values);
	
	
	clear_vec_mp(AtimesF);
	clear_vec_mp(linprod_x);
	clear_vec_mp(linprod_times_gamma_s);
	clear_vec_mp(tempvec);
	clear_vec_mp(tempvec2);
	
	
	clear_mat_mp(Jv_Patch);
	clear_mat_mp(tempmat);
	clear_mat_mp(lin_func_vals);
	
	
	
	clear_mat_mp(AtimesJ);
	clear_mat_mp(Jv_jac);
	clear_mat_mp(temp_jacobian_functions);
	clear_mat_mp(temp_jacobian_parameters);
	clear_mat_mp(linprod_derivative_wrt_x);
	
	clear_mp(running_prod);
	clear_mp(temp);
	clear_mp(temp2);
	clear_mp(temp3);
	
	
	clear_vec_mp(target_function_values);
	clear_vec_mp(target_function_values_times_oneminus_s);
	clear_mat_mp(S_times_Jf_pi);
	
	
	clear_vec_mp(unused_function_values);
	clear_vec_mp(unused_parVals);
	clear_vec_mp(unused_parDer);
	
	
	clear_mat_mp(unused_Jp);
	clear_mat_mp(perturbed_Jv);
	clear_mat_mp(perturbed_AtimesJ);
	clear_mat_mp(tempmat3);
	clear_mat_mp(tempmat1);
	clear_mat_mp(tempmat2);
	
	clear_vec_mp(perturbed_forward_variables);
	clear_vec_mp(perturbed_backward_variables);
	
	
	
	
	
	
	
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_mp dehommed; init_vec_mp(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize_mp(&dehommed,curr_x_vars);
	mpf_out_str (BED->FOUT, 10, 15, pathVars->r);
	fprintf(BED->FOUT," ");
	mpf_out_str (BED->FOUT, 10, 15, pathVars->i);
	fprintf(BED->FOUT," ");
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		mpf_out_str (BED->FOUT, 10, 15, dehommed->coord[ii].r);
		fprintf(BED->FOUT," ");
		mpf_out_str (BED->FOUT, 10, 15, dehommed->coord[ii].i);
		fprintf(BED->FOUT," ");
		
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_mp(dehommed);
#endif
	
	
  return 0;
} 





int nullspacejac_dehom(point_d out_d, point_mp out_mp,
											 int *out_prec,
											 point_d in_d, point_mp in_mp,
											 int in_prec,
											 void const *ED_d, void const *ED_mp)
{
  
  
	
  *out_prec = in_prec;
	
	
	
  if (in_prec < 64)
  { // compute out_d
		nullspacejac_eval_data_d *BED_d = (nullspacejac_eval_data_d *)ED_d;
		
		change_size_vec_d(out_d,BED_d->num_natural_vars-1);
		out_d->size = BED_d->num_natural_vars-1;
		
		
		for (int ii=0; ii<BED_d->num_natural_vars-1; ++ii) {
			div_d(&out_d->coord[ii],&in_d->coord[ii+1],&in_d->coord[0]); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		
		BED_d = NULL;
		
  }
  else
  { // compute out_mp
		
		nullspacejac_eval_data_mp *BED_mp = (nullspacejac_eval_data_mp *)ED_mp;
		// set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);
		
		change_size_vec_mp(out_mp,BED_mp->num_natural_vars-1);
		out_mp->size = BED_mp->num_natural_vars-1;
		
		for (int ii=0; ii<BED_mp->num_natural_vars-1; ++ii) {
			div_mp(&out_mp->coord[ii],&in_mp->coord[ii+1],&in_mp->coord[0]); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		
		BED_mp = NULL;
		
	}
	

	
  
	
  return 0;
}









int change_nullspacejac_eval_prec(void const *ED, int new_prec)
{
	nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED; // to avoid having to cast every time
	
	int ii, jj;
	
	if (new_prec != BED->curr_prec){
		// change the precision for the patch
		changePatchPrec_mp(new_prec, &BED->patch);
		
		if (BED->verbose_level >=4)
			printf("prec  %d\t-->\t%d\n",BED->curr_prec, new_prec);

		BED->SLP->precision = new_prec;
		
		BED->curr_prec = new_prec;
		
		setprec_mp(BED->gamma, new_prec);
		mpf_set_q(BED->gamma->r, BED->gamma_rat[0]);
		mpf_set_q(BED->gamma->i, BED->gamma_rat[1]);
		
		change_prec_mat_mp(BED->randomizer_matrix,new_prec);
		mat_cp_mp(BED->randomizer_matrix,BED->randomizer_matrix_full_prec);
		
		change_prec_mat_mp(BED->post_randomizer_matrix,new_prec);
		mat_cp_mp(BED->post_randomizer_matrix,BED->post_randomizer_matrix_full_prec);
		
		
		for (ii=0; ii<BED->num_additional_linears; ii++) {
			change_prec_point_mp(BED->additional_linears_terminal[ii],new_prec);
			vec_cp_mp(BED->additional_linears_terminal[ii],BED->additional_linears_terminal_full_prec[ii]);
			
			
			change_prec_point_mp(BED->additional_linears_starting[ii],new_prec);
			vec_cp_mp(BED->additional_linears_starting[ii],BED->additional_linears_starting_full_prec[ii]);
		}
		
		
		
		for (ii=0; ii<BED->num_jac_equations; ++ii) {
			for (jj=0; jj<BED->max_degree; jj++) {
				change_prec_point_mp(BED->starting_linears[ii][jj],new_prec);
				vec_cp_mp(BED->starting_linears[ii][jj], BED->starting_linears_full_prec[ii][jj]);
			}
		}
		
		
		
		for (ii=0; ii<BED->num_v_linears; ii++) {
			change_prec_point_mp(BED->v_linears[ii],new_prec);
			vec_cp_mp(BED->v_linears[ii],BED->v_linears_full_prec[ii]);
		}

		
		change_prec_point_mp(BED->v_patch,new_prec);
		vec_cp_mp(BED->v_patch,BED->v_patch_full_prec);
		
		change_prec_mat_mp(BED->jac_with_proj,new_prec);
		mat_cp_mp(BED->jac_with_proj,BED->jac_with_proj_full_prec);
		
		change_prec_mp(BED->perturbation,new_prec);
		change_prec_mp(BED->half,new_prec);
		
		
	}
	
	
  return 0;
}






int check_issoln_nullspacejac_d(endgame_data_t *EG,
																tracker_config_t *T,
																void const *ED)
{
  nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED; // to avoid having to cast every time
	
	
	int ii;
	
	vec_d curr_v_vars; init_vec_d(curr_v_vars, BED->num_v_vars);
	curr_v_vars->size = BED->num_v_vars;
	
	vec_d curr_x_vars; init_vec_d(curr_x_vars, BED->num_x_vars);
	curr_x_vars->size = BED->num_x_vars;
	
	double n1, n2, max_rat;
	point_d f;
	eval_struct_d e;
	//
	//	mpf_init(n1); mpf_init(n2); mpf_init(zero_thresh); mpf_init(max_rat);
	init_point_d(f, 1);
	init_eval_struct_d(e,0, 0, 0);
	
	max_rat = T->ratioTol;
	
	// setup threshold based on given threshold and precision
	//	if (num_digits > 300)
	//		num_digits = 300;
	//	num_digits -= 2;
	double tol = MAX(T->funcResTol, 1e-10);
	
	
	if (EG->prec>=64){
		vec_d terminal_pt;  init_vec_d(terminal_pt,1);
		vec_mp_to_d(terminal_pt,EG->PD_mp.point);
		evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, terminal_pt, EG->PD_d.time, BED->SLP);
		//		lin_to_lin_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, terminal_pt, EG->PD_d.time, ED);
		

		for (ii=0; ii<BED->num_x_vars; ii++)
			set_d(&curr_x_vars->coord[ii], &terminal_pt->coord[ii]);
		
		for (ii=0; ii<BED->num_v_vars; ii++)
			set_d(&curr_v_vars->coord[ii], &terminal_pt->coord[ii+BED->num_x_vars]);
		
		clear_vec_d(terminal_pt);
	}
	else{
		evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_d.point, EG->PD_d.time, BED->SLP);
		//		lin_to_lin_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_d.point, EG->PD_d.time, ED);
		
		for (ii=0; ii<BED->num_x_vars; ii++)
			set_d(&curr_x_vars->coord[ii], &EG->PD_d.point->coord[ii]);
		
		for (ii=0; ii<BED->num_v_vars; ii++)
			set_d(&curr_v_vars->coord[ii], &EG->PD_d.point->coord[ii+BED->num_x_vars]);
		
	}
	
	
	if (EG->last_approx_prec>=64) {
		vec_d prev_pt;  init_vec_d(prev_pt,1);
		vec_mp_to_d(prev_pt,EG->PD_mp.point);
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, prev_pt, EG->PD_d.time, BED->SLP);
		clear_vec_d(prev_pt);}
	else{
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_d, EG->PD_d.time, BED->SLP);
	}
	
	
	
	
	

	
	
//	print_point_to_screen_matlab(EG->PD_d.point,"soln");
//	print_point_to_screen_matlab(e.funcVals,"howfaroff");	// compare the function values
	int isSoln = 1;
	for (ii = 0; (ii < BED->SLP->numFuncs) && isSoln; ii++)
	{
		n1 = d_abs_d( &e.funcVals->coord[ii]); // corresponds to final point
		n2 = d_abs_d( &f->coord[ii]); // corresponds to the previous point
		
		
		if (tol <= n1 && n1 <= n2)
		{ // compare ratio
			if (n1 > max_rat * n2){ // seriously what is the point of this
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (d) 1 coord %d\n",ii);
			}
		}
		else if (tol <= n2 && n2 <= n1)
		{ // compare ratio
			if (n2 > max_rat * n1){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (d) 2 coord %d\n",ii);
			}
		}
	}
	

	
	
	if (!isSoln) {
		
		print_point_to_screen_matlab(e.funcVals,"terminal");
		print_point_to_screen_matlab(f,"prev");

		printf("tol was %le\nmax_rat was %le\n",tol,max_rat);
	}
	
	
	
	
//	mat_d AtimesJ; init_mat_d(AtimesJ,0,0);  AtimesJ->rows = AtimesJ->cols = 0;
//	mat_d jac_homogenizing_matrix;  init_mat_d(jac_homogenizing_matrix,0,0);
//	jac_homogenizing_matrix->rows = jac_homogenizing_matrix->cols = 0;
//	mat_d tempmat; init_mat_d(tempmat,0,0); tempmat->rows = tempmat->cols = 0;
//	
//	mat_d Jf_pi; init_mat_d(Jf_pi,0,0); Jf_pi->rows = Jf_pi->cols = 0;
//	vec_d target_function_values; init_vec_d(target_function_values,0); target_function_values->size = 0;
//	
//	mat_mul_d(AtimesJ,BED->randomizer_matrix,e.Jv);
//	for (ii=0; ii< AtimesJ->rows; ii++)
//		for (int jj=1; jj<BED->num_x_vars; jj++)
//			set_d(&BED->jac_with_proj->entry[jj - 1][ii],&AtimesJ->entry[ii][jj]); // copy in the transpose of the (randomized) jacobian, omitting the homogenizing variables
//	
//	make_matrix_ID_d(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
//	
//	for (ii=0; ii<BED->num_randomized_eqns; ii++)
//		for (int jj=0; jj<(BED->max_degree - (BED->randomized_degrees[ii]-1)); jj++)
//			mul_d(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
//	
//	for (ii=BED->num_randomized_eqns; ii<BED->num_v_vars; ii++)
//		for (int jj=0; jj<(BED->max_degree); jj++) // these are all degree 1
//			mul_d(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
//	
//	mat_mul_d(Jf_pi, BED->jac_with_proj, jac_homogenizing_matrix);
//	
//	mul_mat_vec_d(target_function_values, Jf_pi, curr_v_vars);
//	
//	mul_mat_vec_d(target_function_values, BED->post_randomizer_matrix, target_function_values);
//
//	
//	for (ii=0; ii<target_function_values->size; ii++) {
//		if (d_abs_d(&target_function_values->coord[ii]) > tol) {
//			isSoln = 0;
//			print_point_to_screen_matlab(target_function_values,"target_func_vals");
//			break;
//		}
//	}
	
	clear_eval_struct_d(e);
	clear_vec_d(f);
	
//	std::cout << isSoln << std::endl;
	
	
	return isSoln;
	
}


int check_issoln_nullspacejac_mp(endgame_data_t *EG,
																 tracker_config_t *T,
																 void const *ED)
{
  nullspacejac_eval_data_mp *BED = (nullspacejac_eval_data_mp *)ED; // to avoid having to cast every time
	
	int ii;
	
	vec_mp curr_v_vars; init_vec_mp(curr_v_vars, BED->num_v_vars);
	curr_v_vars->size = BED->num_v_vars;
	
	vec_mp curr_x_vars; init_vec_mp(curr_x_vars, BED->num_x_vars);
	curr_x_vars->size = BED->num_x_vars;
	
	
	
	for (ii = 0; ii < T->numVars; ii++)
	{
    if (!(mpfr_number_p(EG->PD_mp.point->coord[ii].r) && mpfr_number_p(EG->PD_mp.point->coord[ii].i)))
		{
			printf("got not a number\n");
			print_point_to_screen_matlab(EG->PD_mp.point,"bad solution");
      return 0;
		}
	}
	
	
	
	mpf_t n1, n2, zero_thresh, max_rat;
	mpf_init(n1); mpf_init(n2); mpf_init(zero_thresh); mpf_init(max_rat);
	
	point_mp f; init_point_mp(f, 1);f->size = 1;
	eval_struct_mp e; init_eval_struct_mp(e, 0, 0, 0);
	
	mpf_set_d(max_rat, T->ratioTol);
	
	
	int num_digits = prec_to_digits((int) mpf_get_default_prec());
	// setup threshold based on given threshold and precision
	if (num_digits > 300)
		num_digits = 300;
	num_digits -= 4;
	double tol = MAX(T->funcResTol, pow(10,-num_digits));
	mpf_set_d(zero_thresh, tol);
	
	
	for (ii=0; ii<BED->num_x_vars; ii++)
		set_mp(&curr_x_vars->coord[ii], &EG->PD_mp.point->coord[ii]);
	
	for (ii=0; ii<BED->num_v_vars; ii++)
		set_mp(&curr_v_vars->coord[ii], &EG->PD_mp.point->coord[ii+BED->num_x_vars]);
	
	
	
	//this one guaranteed by entry condition
	//	lin_to_lin_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, ED);
	evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, BED->SLP);
	
//	print_point_to_screen_matlab(EG->PD_mp.point,"soln");
//	print_point_to_screen_matlab(e.funcVals,"howfaroff");
	
	if (EG->last_approx_prec < 64) { // copy to _mp
		point_d_to_mp(EG->last_approx_mp, EG->last_approx_d);
	}
	
	evalProg_mp(f, e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_mp, EG->PD_mp.time, BED->SLP);
	//	lin_to_lin_eval_mp(f,          e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_mp, EG->PD_mp.time, ED);
	// compare the function values
	int isSoln = 1;
	for (ii = 0; ii < BED->SLP->numFuncs && isSoln; ii++)
	{
		mpf_abs_mp(n1, &e.funcVals->coord[ii]);
		mpf_abs_mp(n2, &f->coord[ii]);
		
		//		mpf_out_str(NULL,10,9,n1);
		
		if ( (mpf_cmp(zero_thresh, n1) <= 0) &&  (mpf_cmp(n1, n2) <= 0) )
		{ // compare ratio
			mpf_mul(n2, max_rat, n2);
			if (mpf_cmp(n1, n2) > 0){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (mp) 1\n");
			}
		}
		else if ( (mpf_cmp(zero_thresh, n2) <= 0) &&  (mpf_cmp(n2, n1) <= 0) )
		{ // compare ratio
			mpf_mul(n1, max_rat, n1);
			if (mpf_cmp(n2, n1) > 0){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (mp) 2\n");
			}
		}
	}
	
	
//	mat_mp AtimesJ; init_mat_mp(AtimesJ,0,0);  AtimesJ->rows = AtimesJ->cols = 0;
//	mat_mp jac_homogenizing_matrix;  init_mat_mp(jac_homogenizing_matrix,0,0);
//	jac_homogenizing_matrix->rows = jac_homogenizing_matrix->cols = 0;
//	mat_mp tempmat; init_mat_mp(tempmat,0,0); tempmat->rows = tempmat->cols = 0;
//	mat_mp S_times_Jf_pi; init_mat_mp(S_times_Jf_pi,0,0); S_times_Jf_pi->rows = S_times_Jf_pi->cols = 0;
//	vec_mp target_function_values; init_vec_mp(target_function_values,0); target_function_values->size = 0;
//	
//	mat_mul_mp(AtimesJ,BED->randomizer_matrix,e.Jv);
//	for (ii=0; ii< AtimesJ->rows; ii++)
//		for (int jj=1; jj<BED->num_x_vars; jj++)
//			set_mp(&BED->jac_with_proj->entry[jj - 1][ii],&AtimesJ->entry[ii][jj]); // copy in the transpose of the (randomized) jacobian, omitting the homogenizing variables
//	make_matrix_ID_mp(jac_homogenizing_matrix,BED->num_v_vars,BED->num_v_vars);
//	
//	for (ii=0; ii<BED->num_randomized_eqns; ii++)
//		for (int jj=0; jj<(BED->max_degree - (BED->randomized_degrees[ii]-1)); jj++)
//			mul_mp(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
//	
//	for (ii=BED->num_randomized_eqns; ii<BED->num_v_vars; ii++)
//		for (int jj=0; jj<(BED->max_degree); jj++) // these are all degree 1
//			mul_mp(&jac_homogenizing_matrix->entry[ii][ii], &jac_homogenizing_matrix->entry[ii][ii], &curr_x_vars->coord[0]);
//	
//	mat_mul_mp(tempmat, BED->post_randomizer_matrix, BED->jac_with_proj); // jac with proj having been set above with unperturbed values
//	mat_mul_mp(S_times_Jf_pi, tempmat, jac_homogenizing_matrix); // jac with proj having been set above with unperturbed values
//	mul_mat_vec_mp(target_function_values, S_times_Jf_pi, curr_v_vars);
//	
//	
//	for (ii=0; ii<target_function_values->size; ii++) {
//		mpf_abs_mp(n1, &target_function_values->coord[ii]);
//		if (mpf_cmp(n1,zero_thresh)>0) {
//			isSoln = 0;
//			std::cout << T->funcResTol << " " <<  pow(10,-num_digits) << std::endl;
//			mpf_out_str (NULL, 10, 6, zero_thresh); // base 10, 6 digits
//			std::cout << std::endl;
//			print_point_to_screen_matlab(target_function_values,"target_func_vals");
//			break;
//		}
//	}
	
	
	
	mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);
	
	
	clear_eval_struct_mp(e);
	clear_vec_mp(f);
	
//	std::cout << isSoln << std::endl;
	
	return isSoln;
	
}







int check_isstart_nullspacejac_d(point_d testpoint,
																tracker_config_t *T,
																void const *ED)
{

	eval_struct_d e;
	init_eval_struct_d(e,0, 0, 0);
		
	comp_d time;
	set_one_d(time);
	
	
	double tol = (1e-9);
	
	nullspacejac_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, testpoint, time, ED);
	
	int isSoln = 1;
	
	for (int ii = 0; (ii < e.funcVals->size) && isSoln; ii++) // function by function
	{		
		if (tol <= d_abs_d( &e.funcVals->coord[ii])){ // compare
			isSoln = 0;
			print_point_to_screen_matlab(testpoint,"invalid_startpoint");
			print_point_to_screen_matlab(e.funcVals,"start_residual");
		}
		
	}
	
	
	clear_eval_struct_d(e);
	
	return isSoln;
	
}



void check_nullspace_evaluator(point_mp current_values,
															 void const *ED)
{
	int ii;
	printf("checking homogeneousness of double evaluator\n");
  nullspacejac_eval_data_d *BED = (nullspacejac_eval_data_d *)ED; // to avoid having to cast every time
	//initialize
	eval_struct_d e_d; init_eval_struct_d(e_d, 0, 0, 0);
	eval_struct_d e_d2; init_eval_struct_d(e_d2, 0, 0, 0);
	
	

	
	comp_d zerotime; set_zero_d(zerotime);
	

	
	
	point_d tempvec;  init_point_d(tempvec,0);
	vec_mp_to_d(tempvec, current_values);
	
	nullspacejac_eval_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, zerotime, ED);

	
	comp_d lambda; get_comp_rand_d(lambda);
	
	for (ii=0; ii<BED->num_x_vars; ii++) {
		mul_d(&tempvec->coord[ii],&tempvec->coord[ii],lambda);
	}
	
	nullspacejac_eval_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec, zerotime, ED);
	

	printf("lambda = %lf+1i*%lf\n",lambda->r, lambda->i);
	print_point_to_screen_matlab(e_d.funcVals,"f");
	print_point_to_screen_matlab(e_d2.funcVals,"f2");
	
	
	mypause();

	return;
	
}















