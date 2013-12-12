#include "solver_midpoint_tracker.hpp"
#include "surface.hpp"


void midpoint_config::setup(const surface_decomposition & surf,
                            solver_configuration & solve_options)
{
	
	this->MPType = solve_options.T.MPType;
    
    add_projection(surf.pi[0]);
	add_projection(surf.pi[1]);
    
	num_mid_vars = surf.num_variables;
	num_crit_vars = surf.crit_curve.num_variables;
	num_sphere_vars = surf.sphere_curve.num_variables;
	
	
	int blabla;  // i would like to move this.
	parse_input_file(surf.input_filename, &blabla);
    
	
	//	// setup a straight-line program, using the file(s) created by the parser
	solve_options.T.numVars = setupProg(this->SLP_mid, solve_options.T.Precision, solve_options.T.MPType);
    
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	//create the array of integers
	std::vector<int> randomized_degrees;
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(this->randomizer_matrix, randomized_degrees, surf.num_variables-surf.num_patches-2, solve_options.PPD.num_funcs);
	//W.num_variables-W.num_patches-W.num_linears
    
	
	this->mid_memory.capture_globals();
	this->mid_memory.set_globals_null();
	
    
    
	//do the necessary parsing
	
    
	parse_input_file(surf.crit_curve.input_filename, &blabla);
    
	
	
	//	// setup a straight-line program, using the file(s) created by the parser
	solve_options.T.numVars = setupProg(this->SLP_crit, solve_options.T.Precision, solve_options.T.MPType);
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	//get the matrix and the degrees of the resulting randomized functions.
	randomized_degrees.clear();
	make_randomization_matrix_based_on_degrees(this->randomizer_matrix_crit, randomized_degrees, surf.crit_curve.num_variables-surf.crit_curve.num_patches-1, solve_options.PPD.num_funcs);
	
	this->crit_memory.capture_globals();
	this->crit_memory.set_globals_null();
	
	
	
	//do the necessary parsing
	
	
	parse_input_file(surf.sphere_curve.input_filename, &blabla);
    
	
	//	// setup a straight-line program, using the file(s) created by the parser
	setupProg(this->SLP_sphere, solve_options.T.Precision, solve_options.T.MPType);
	preproc_data_clear(&solve_options.PPD);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	
	//get the matrix and the degrees of the resulting randomized functions.
	randomized_degrees.clear();
	make_randomization_matrix_based_on_degrees(this->randomizer_matrix_sph, randomized_degrees, surf.sphere_curve.num_variables-surf.sphere_curve.num_patches-1, solve_options.PPD.num_funcs);
	
	this->sphere_memory.capture_globals();
	this->sphere_memory.set_globals_null();
	
	
    //	// are there other things in T which need to be set?
}


void midpoint_config::init()
{
	num_projections = 0;
    
	this->MPType = -1;
	
	init_mat_mp2(randomizer_matrix_sph,1,1,1024);randomizer_matrix_sph->rows = randomizer_matrix_sph->cols = 1;
	init_mat_mp2(randomizer_matrix_crit,1,1,1024);randomizer_matrix_crit->rows = randomizer_matrix_crit->cols = 1;
	init_mat_mp2(randomizer_matrix,1,1,1024);randomizer_matrix->rows = randomizer_matrix->cols = 1;
	
	
	init_mp2(v_target,1024);
	init_mp2(u_target,1024);
	
	init_mp2(crit_val_left,1024);
	init_mp2(crit_val_right,1024);
	
	this->SLP_crit =  new prog_t;
	this->SLP_mid = new prog_t;//(prog_t *) br_malloc(sizeof(prog_t));
	this->SLP_sphere = new prog_t;//(prog_t *) br_malloc(sizeof(prog_t));
	
	system_type_top = UNSET;
	system_type_bottom = UNSET;
	
    
}



void midpoint_config::initial_send(parallelism_config & mpi_config)
{
	
	int * buffer = new int[7];
	
	buffer[0] = randomizer_matrix_crit->rows;
	buffer[1] = randomizer_matrix_crit->cols;
	buffer[2] = randomizer_matrix_sph->rows;
	buffer[3] = randomizer_matrix_sph->cols;
	buffer[4] = randomizer_matrix->rows;
	buffer[5] = randomizer_matrix->cols;
	buffer[6] = MPType;
	
	MPI_Bcast(buffer, 7, MPI_INT, mpi_config.id(), mpi_config.my_communicator);
	
	delete [] buffer;
	
	
	
	mid_memory.set_globals_to_this();
    bcast_prog_t(SLP_mid, this->MPType, 0, 0); // last two arguments are: myid, headnode
	
	sphere_memory.set_globals_to_this();
    bcast_prog_t(SLP_sphere, this->MPType, 0, 0); // last two arguments are: myid, headnode
	
    crit_memory.set_globals_to_this();
    bcast_prog_t(SLP_crit, this->MPType, 0, 0); // last two arguments are: myid, headnode

    

    
	
	
	
	
	
    bcast_mat_mp(randomizer_matrix_crit,0,0);
    bcast_mat_mp(randomizer_matrix_sph,0,0);
    bcast_mat_mp(randomizer_matrix,0,0);
    
	
	
    
    MPI_Bcast(&num_mid_vars, 1, MPI_INT, mpi_config.id(), mpi_config.my_communicator);
    MPI_Bcast(&num_crit_vars, 1, MPI_INT, mpi_config.id(), mpi_config.my_communicator);
    MPI_Bcast(&num_sphere_vars, 1, MPI_INT, mpi_config.id(), mpi_config.my_communicator);
	
    MPI_Bcast(&num_projections, 1, MPI_INT, mpi_config.id(), mpi_config.my_communicator);
    for (int ii=0; ii<num_projections; ii++) {
        bcast_vec_mp(pi[ii],0,0);
    }
	
	
}

void midpoint_config::initial_receive(parallelism_config & mpi_config)
{
	
	
	int * buffer = new int[7];
	
	MPI_Bcast(buffer, 7, MPI_INT, mpi_config.id(), mpi_config.my_communicator);
	
	change_size_mat_mp(randomizer_matrix_crit,buffer[0],buffer[1]);
	change_size_mat_mp(randomizer_matrix_sph,buffer[2],buffer[3]);
	change_size_mat_mp(randomizer_matrix,buffer[4],buffer[5]);
	
	
	randomizer_matrix_crit->rows = buffer[0];
	randomizer_matrix_crit->cols = buffer[1];
	randomizer_matrix_sph->rows = buffer[2];
	randomizer_matrix_sph->cols = buffer[3];
	randomizer_matrix->rows = buffer[4];
	randomizer_matrix->cols = buffer[5];
	
	MPType = buffer[6];
	
	delete [] buffer;
	
	
	
	
	
    bcast_prog_t(SLP_mid, this->MPType, 1, 0); // last two arguments are: myid, headnode
	initEvalProg(this->MPType);
    mid_memory.capture_globals();
    mid_memory.set_globals_null();
	
    sphere_memory.set_globals_null();
	bcast_prog_t(SLP_sphere, this->MPType, 1, 0); // last two arguments are: myid, headnode
	initEvalProg(this->MPType);
    sphere_memory.capture_globals();
    sphere_memory.set_globals_null();
	
	
    bcast_prog_t(SLP_crit, this->MPType, 1, 0); // last two arguments are: myid, headnode
	initEvalProg(this->MPType);
    crit_memory.capture_globals();
	crit_memory.set_globals_null();
	
	
    
	
	
	
	
	
    bcast_mat_mp(randomizer_matrix_crit,1,0);
    bcast_mat_mp(randomizer_matrix_sph,1,0);
    bcast_mat_mp(randomizer_matrix,1,0);
    
    
    MPI_Bcast(&num_mid_vars, 1, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
    MPI_Bcast(&num_crit_vars, 1, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
    MPI_Bcast(&num_sphere_vars, 1, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
	int temp_num_projections;
	MPI_Bcast(&temp_num_projections, 1, MPI_INT, mpi_config.id(), mpi_config.my_communicator);
    vec_mp tempvec;  init_vec_mp2(tempvec, 1, 1024);
    for (int ii=0; ii<temp_num_projections; ii++) {
        bcast_vec_mp(tempvec,1,0);
        add_projection(tempvec);
    }
    clear_vec_mp(tempvec);
	
	
}

//void midpoint_config::update_send(parallelism_config & mpi_config,int target)
//{
//    send_comp_mp(v_target,target); // set during the loop in connect the dots
//    send_comp_mp(u_target,target); // set during the loop in connect the dots
//    
//    send_comp_mp(crit_val_left,target); // set during the loop in connect the dots
//    send_comp_mp(crit_val_right,target); // set during the loop in connect the dots
//    
//
//	MPI_Send(&system_type_bottom, 1, MPI_INT, target, UNUSED, MPI_COMM_WORLD);
//    MPI_Send(&system_type_top, 1, MPI_INT, target, UNUSED, MPI_COMM_WORLD);
//}

//void midpoint_config::update_receive(parallelism_config & mpi_config)
//{
//    MPI_Status statty_mc_gatty;
//    
//    receive_comp_mp(v_target,MPI_ANY_SOURCE); // set during the loop in connect the dots
//    receive_comp_mp(u_target,MPI_ANY_SOURCE); // set during the loop in connect the dots
//    
//    receive_comp_mp(crit_val_left,MPI_ANY_SOURCE); // set during the loop in connect the dots
//    receive_comp_mp(crit_val_right,MPI_ANY_SOURCE); // set during the loop in connect the dots
//    
//    
//    MPI_Recv(&system_type_bottom, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &statty_mc_gatty);
//    MPI_Recv(&system_type_top, 1, MPI_INT, MPI_ANY_SOURCE, MPI_ANY_TAG, MPI_COMM_WORLD, &statty_mc_gatty);
//}

///////////////
//
//   begin midpoint_eval_data_mp
//
/////////////


void midpoint_eval_data_mp::init()
{
	this->is_solution_checker_d = &check_issoln_midpoint_d;
	this->is_solution_checker_mp = &check_issoln_midpoint_mp;
	this->evaluator_function_d = &midpoint_eval_d;
	this->evaluator_function_mp = &midpoint_eval_mp;
	this->precision_changer = &change_midpoint_eval_prec;
	this->dehomogenizer = &midpoint_dehom;
	
	num_projections = 0;
	pi = NULL;
	SLP_bottom = NULL;
	SLP_mid = NULL;
	SLP_top = NULL;
	
	init_mp(half);
	init_mp(one);
	init_mp(zero);
	init_mp(u_start);
	init_mp(v_start);
	init_mp(u_target);
	init_mp(v_target);
	init_mp(crit_val_right);
	init_mp(crit_val_left);
	
	
	comp_d h;
	
	h->r = 0.5;
	h->i = 0.0;
	
	d_to_mp(half, h);
    
	set_one_mp(this->one);
	set_zero_mp(this->zero);
	
	
	set_mp(this->u_start,half);
	set_mp(this->v_start,half);
	
	
	
	this->num_projections = 0;
	this->num_variables = -1;
	this->num_mid_vars = -1;
	this->num_bottom_vars = -1;
	this->num_top_vars = -1;
	
    
    
	
	
	
	init_mat_mp(randomizer_matrix_bottom,1,1);randomizer_matrix_bottom->rows = randomizer_matrix_bottom->cols = 1;
	init_mat_mp(randomizer_matrix_top,1,1);randomizer_matrix_top->rows = randomizer_matrix_top->cols = 1;
	
	
	
	if (this->MPType==2) {
		init_mat_mp2(randomizer_matrix_bottom_full_prec, 0, 0,1024);
		init_mat_mp2(randomizer_matrix_top_full_prec, 0, 0,1024);
		
		init_mp2(half_full_prec,1024);
		set_zero_mp(half_full_prec);
		mpf_set_str(half_full_prec->r,"0.5",10);
		
		
		init_mp2(crit_val_right_full_prec,1024);
		init_mp2(crit_val_left_full_prec,1024);
		
		
		init_mp2(one_full_prec,1024);
		init_mp2(zero_full_prec,1024);
		init_mp2(u_start_full_prec,1024);
		init_mp2(v_start_full_prec,1024);
		init_mp2(u_target_full_prec,1024);
		init_mp2(v_target_full_prec,1024);
		
		set_one_mp(this->one_full_prec);
		set_zero_mp(this->zero_full_prec);
		
		
		set_mp(this->u_start_full_prec,half_full_prec);
		set_mp(this->v_start_full_prec,half_full_prec);
		
	}
	
	
	
    
	
	
	
}


int midpoint_eval_data_mp::send(parallelism_config & mpi_config)
{
	
	int solver_choice = MIDPOINT_SOLVER;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
	// send the confirmation integer, to ensure that we are sending the correct type.
	
	//send the base class stuff.
	solver_mp::send(mpi_config);
	
	int *buffer = new int[12];
	
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	delete[] buffer;
    
	return SUCCESSFUL;
}




int midpoint_eval_data_mp::receive(parallelism_config & mpi_config)
{
    
	int *buffer = new int[12]; // allocate 12 here because we will be sending 12 integers
	MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (buffer[0] != MIDPOINT_SOLVER) {
		std::cout << "worker failed to confirm it is receiving the midpoint solver type eval data" << std::endl;
		mpi_config.abort(778);
	}
	// now can actually receive the data from whomever.
    
    //the base class receive
	solver_mp::receive(mpi_config);
	
	
    
    
    
    
    
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	
	
	delete[] buffer;
    
	
	return SUCCESSFUL;
}



int midpoint_eval_data_mp::setup(const midpoint_config & md_config,
                                 const witness_set & W,
                                 solver_configuration & solve_options)
{
	
	verbose_level = solve_options.verbose_level;
	
	solver_mp::setup();
	
	generic_setup_patch(&patch,W);
	
	
	
	this->mid_memory = md_config.mid_memory;
	this->SLP_mid = md_config.SLP_mid;
	this->num_mid_vars = md_config.num_mid_vars;
	mat_cp_mp(randomizer_matrix, md_config.randomizer_matrix);
	
	
	if (md_config.system_type_top == SYSTEM_CRIT) {
		this->top_memory = md_config.crit_memory; // copy the pointers to the memory
		this->SLP_top = md_config.SLP_crit;
		this->num_top_vars = md_config.num_crit_vars;
		mat_cp_mp(randomizer_matrix_top, md_config.randomizer_matrix_crit);
		if (MPType == 2) {
			mat_cp_mp(randomizer_matrix_top_full_prec, md_config.randomizer_matrix_crit);
		}
	}
	else if (md_config.system_type_top == SYSTEM_SPHERE){
		this->top_memory = md_config.sphere_memory; // copy the pointers to the memory
		this->SLP_top = md_config.SLP_sphere;
		this->num_top_vars = md_config.num_sphere_vars;
		mat_cp_mp(randomizer_matrix_top, md_config.randomizer_matrix_sph);
		if (MPType == 2) {
			mat_cp_mp(randomizer_matrix_top_full_prec, md_config.randomizer_matrix_sph);
		}
	}
	else{
		br_exit(98711);
	}
	
	
	
	
	if (md_config.system_type_bottom == SYSTEM_CRIT) {
		this->bottom_memory = md_config.crit_memory; // copy the pointers to the memory
		this->SLP_bottom = md_config.SLP_crit;
		this->num_bottom_vars = md_config.num_crit_vars;
		mat_cp_mp(randomizer_matrix_bottom, md_config.randomizer_matrix_crit);
		if (MPType == 2) {
			mat_cp_mp(randomizer_matrix_bottom_full_prec, md_config.randomizer_matrix_crit);
		}
	}
	else if (md_config.system_type_bottom == SYSTEM_SPHERE){
		this->bottom_memory = md_config.sphere_memory; // copy the pointers to the memory
		this->SLP_bottom = md_config.SLP_sphere;
		this->num_bottom_vars = md_config.num_sphere_vars;
		mat_cp_mp(randomizer_matrix_bottom, md_config.randomizer_matrix_sph);
		if (MPType == 2) {
			mat_cp_mp(randomizer_matrix_bottom_full_prec, md_config.randomizer_matrix_sph);
		}
	}
	else{
		br_exit(98711);
	}
	
	
    
	this->num_variables = num_mid_vars + num_bottom_vars + num_top_vars;
	
    
	
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
	
	
    
	set_mp(this->crit_val_left, md_config.crit_val_left);
	set_mp(this->crit_val_right, md_config.crit_val_right);
	
	set_mp(this->u_target, md_config.u_target);
	set_mp(this->v_target, md_config.v_target);
	
    
	
	add_projection(md_config.pi[0]);
	add_projection(md_config.pi[1]);
    
	
	
	
	
	
	if (this->MPType==2) {
		mat_cp_mp(randomizer_matrix_full_prec, md_config.randomizer_matrix);
        
		set_mp(this->crit_val_left_full_prec, md_config.crit_val_left);
		set_mp(this->crit_val_right_full_prec, md_config.crit_val_right);
		
		set_mp(this->u_target_full_prec, md_config.u_target);
		set_mp(this->v_target_full_prec, md_config.v_target);
	}
	
	
	return SUCCESSFUL;
}

///////////////
//
//   end midpoint_eval_data_mp
//
/////////////



















///////////////
//
//   begin midpoint_eval_data_d
//
/////////////

void midpoint_eval_data_d::init()
{
	
	if (this->MPType==2)
    {
		this->BED_mp = new midpoint_eval_data_mp(2);
        solver_d::BED_mp = this->BED_mp;
    }
	else{
		this->BED_mp = NULL;
    }
	
	this->is_solution_checker_d = &check_issoln_midpoint_d;
	this->is_solution_checker_mp = &check_issoln_midpoint_mp;
	this->evaluator_function_d = &midpoint_eval_d;
	this->evaluator_function_mp = &midpoint_eval_mp;
	this->precision_changer = &change_midpoint_eval_prec;
	this->dehomogenizer = &midpoint_dehom;
	
	half->r = 0.5;
	half->i = 0.0;
    
	set_one_d(this->one);
	set_zero_d(this->zero);
	
	set_d(this->u_start,half);
	set_d(this->v_start,half);
	
	this->pi = NULL;
	this->num_projections = 0;
	this->num_variables = -1;
	this->num_mid_vars = -1;
	this->num_bottom_vars = -1;
	this->num_top_vars = -1;
	
    
	
	init_mat_d(randomizer_matrix_bottom,1,1);randomizer_matrix_bottom->rows = randomizer_matrix_bottom->cols = 1;
	init_mat_d(randomizer_matrix_top,1,1); randomizer_matrix_top->rows = randomizer_matrix_top->cols = 1;
}


int midpoint_eval_data_d::send(parallelism_config & mpi_config)
{
	
	int solver_choice = MIDPOINT_SOLVER;
	MPI_Bcast(&solver_choice, 1, MPI_INT, mpi_config.head(), mpi_config.my_communicator);
	// send the confirmation integer, to ensure that we are sending the correct type.
	
	//send the base class stuff.
	solver_d::send(mpi_config);
	
	
	
	int *buffer = new int[12];
	
	// now can actually send the data.
	
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	delete[] buffer;
    return SUCCESSFUL;
}

int midpoint_eval_data_d::receive(parallelism_config & mpi_config)
{
	int *buffer = new int[12];
	MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (buffer[0] != MIDPOINT_SOLVER){
		std::cout << "worker failed to confirm it is receiving the midpoint solver type eval data" << std::endl;
		mpi_config.abort(777);
	}
	
	solver_d::receive(mpi_config);
	
	// now can actually receive the data from whoever.
	
	
	
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	
    
	
	
	delete[] buffer;
    
	
	if (this->MPType==2) {
		this->BED_mp->receive(mpi_config);
	}
	
	return SUCCESSFUL;
}




int midpoint_eval_data_d::setup(const midpoint_config & md_config,
                                const witness_set & W,
                                solver_configuration & solve_options)
{
	
	solver_d::setup();
	
	
	
	verbose_level = solve_options.verbose_level;
	
	generic_setup_patch(&patch,W);
	
	
	
	this->mid_memory = md_config.mid_memory;
	this->SLP_mid = md_config.SLP_mid;
	this->num_mid_vars = md_config.num_mid_vars;
	mat_mp_to_d(randomizer_matrix, md_config.randomizer_matrix);
	
	
	if (md_config.system_type_top == SYSTEM_CRIT) {
		this->top_memory = md_config.crit_memory; // copy the pointers to the memory
		this->SLP_top = md_config.SLP_crit;
		this->num_top_vars = md_config.num_crit_vars;
		mat_mp_to_d(randomizer_matrix_top, md_config.randomizer_matrix_crit);
	}
	else if (md_config.system_type_top == SYSTEM_SPHERE){
		this->top_memory = md_config.sphere_memory; // copy the pointers to the memory
		this->SLP_top = md_config.SLP_sphere;
		this->num_top_vars = md_config.num_sphere_vars;
		mat_mp_to_d(randomizer_matrix_top, md_config.randomizer_matrix_sph);
	}
	else{
		std::cout << "in double setup for midpoint eval, top system was neither sphere nor crit." << std::endl;
		br_exit(98711);
	}
	
	
	
	
	if (md_config.system_type_bottom == SYSTEM_CRIT) {
		this->bottom_memory = md_config.crit_memory; // copy the pointers to the memory
		this->SLP_bottom = md_config.SLP_crit;
		this->num_bottom_vars = md_config.num_crit_vars;
		mat_mp_to_d(randomizer_matrix_bottom, md_config.randomizer_matrix_crit);
	}
	else if (md_config.system_type_bottom == SYSTEM_SPHERE){
		this->bottom_memory = md_config.sphere_memory; // copy the pointers to the memory
		this->SLP_bottom = md_config.SLP_sphere;
		this->num_bottom_vars = md_config.num_sphere_vars;
		mat_mp_to_d(randomizer_matrix_bottom, md_config.randomizer_matrix_sph);
	}
	else{
		std::cout << "in double setup for midpoint eval, bottom system was neither sphere nor crit." << std::endl;
		br_exit(98712);
	}
	
	
	
	this->num_variables = num_mid_vars + num_bottom_vars + num_top_vars;
	
    
	
	
	
	
	add_projection(md_config.pi[0]);
	add_projection(md_config.pi[1]);
	
	mp_to_d(this->crit_val_left, md_config.crit_val_left);
	mp_to_d(this->crit_val_right, md_config.crit_val_right);
	
	mp_to_d(this->u_target, md_config.u_target);
	mp_to_d(this->v_target, md_config.v_target);
	
	if (solve_options.use_gamma_trick==1)
		get_comp_rand_d(this->gamma); // set gamma to be random complex value
	else
		set_one_d(this->gamma);
	
	
	
	mat_mp_to_d(randomizer_matrix,
                md_config.randomizer_matrix);
	
	
	
	if (this->MPType==2)
	{
		this->BED_mp->setup(md_config, W, solve_options); // must be called before the gamma line below.
		rat_to_d(this->gamma, this->BED_mp->gamma_rat);
	}
	
	return SUCCESSFUL;
}

///////////////
//
//   end midpoint_eval_data_d
//
/////////////






















int midpoint_solver_master_entry_point(const witness_set						&W, // carries with it the start points, and the linears.
                                       witness_set						*W_new, // new data goes in here
                                       midpoint_config & md_config,
                                       solver_configuration		& solve_options)
{
    
	bool prev_state = solve_options.force_no_parallel;
	solve_options.force_no_parallel = true;
	
	
	if (solve_options.use_parallel()) {
		solve_options.call_for_help(MIDPOINT);
	}
	
	midpoint_eval_data_d *ED_d = NULL;
	midpoint_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new midpoint_eval_data_d(0);
			
			ED_d->setup(md_config,
                        W,
                        solve_options);
			break;
			
		case 1:
			ED_mp = new midpoint_eval_data_mp(1);
			
			ED_mp->setup(md_config,
                         W,
                         solve_options);
			// initialize latest_newton_residual_mp
			mpf_init(solve_options.T.latest_newton_residual_mp);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new midpoint_eval_data_d(2);
			
			ED_mp = ED_d->BED_mp;
			
			
			ED_d->setup(md_config,
                        W,
                        solve_options);
			
			
			
			adjust_tracker_AMP( &(solve_options.T), W.num_variables);
			// initialize latest_newton_residual_mp
			break;
		default:
			break;
	}
	
	print_point_to_screen_matlab(W.pts_mp[0],"start_pt");
	ED_d->print();
	
	master_solver(W_new, W,
                  ED_d, ED_mp,
                  solve_options);
    
	solve_options.force_no_parallel = prev_state;
	
	switch (solve_options.T.MPType) {
		case 0:
			delete ED_d;
			break;
			
		case 1:
			delete ED_mp;
			break;
		case 2:
			delete ED_d;
			
		default:
			break;
	}
	
    return SUCCESSFUL;
	
}






void midpoint_slave_entry_point(solver_configuration & solve_options)
{
	
	
	// already received the flag which indicated that this worker is going to be performing the nullspace calculation.
	bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
	
	midpoint_eval_data_d *ED_d = NULL;
	midpoint_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new midpoint_eval_data_d(0);
			ED_d->receive(solve_options);
			break;
			
		case 1:
			ED_mp = new midpoint_eval_data_mp(1);
			
			
			ED_mp->receive(solve_options);
			// initialize latest_newton_residual_mp
			mpf_init(solve_options.T.latest_newton_residual_mp);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new midpoint_eval_data_d(2);
			
			ED_d->receive(solve_options);
			
			ED_mp = ED_d->BED_mp;
            
			
			
			
			// initialize latest_newton_residual_mp
			mpf_init2(solve_options.T.latest_newton_residual_mp,solve_options.T.AMP_max_prec);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		default:
			break;
	}
	
	
    
	// call the file setup function
	FILE *OUT = NULL, *midOUT = NULL;
	
	generic_setup_files(&OUT, "output",
						&midOUT, "midpath_data");
	
	trackingStats trackCount; init_trackingStats(&trackCount); // initialize trackCount to all 0
	
    
	worker_tracker_loop(&trackCount, OUT, midOUT,
                        ED_d, ED_mp,
                        solve_options);
	
	
	// close the files
	fclose(midOUT);   fclose(OUT);
	
	
	//clear data
}






int midpoint_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	
	
	
    midpoint_eval_data_d *BED = (midpoint_eval_data_d *)ED; // to avoid having to cast every time
	
	
	
    int ii, jj;
	int offset;
    comp_d one_minus_s, gamma_s;
	
	set_one_d(one_minus_s);
    sub_d(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
    mul_d(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	//parse out the variables into proper segments.
	vec_d curr_mid_vars; init_vec_d(curr_mid_vars, BED->num_mid_vars);
	curr_mid_vars->size = BED->num_mid_vars;
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_d(&curr_mid_vars->coord[ii], &current_variable_values->coord[ii]);
	
	offset = BED->num_mid_vars; // y0
	vec_d curr_bottom_vars; init_vec_d(curr_bottom_vars, BED->num_bottom_vars);
	curr_bottom_vars->size = BED->num_bottom_vars;
	for (ii=0; ii<BED->num_bottom_vars; ii++)
		set_d(&curr_bottom_vars->coord[ii], &current_variable_values->coord[ii+offset]);
	
	offset = BED->num_mid_vars + BED->num_bottom_vars; // y2
	vec_d curr_top_vars; init_vec_d(curr_top_vars, BED->num_top_vars);
	curr_top_vars->size = BED->num_top_vars;
	for (ii=0; ii<BED->num_top_vars; ii++)
		set_d(&curr_top_vars->coord[ii], &current_variable_values->coord[ii+offset]);
	
    
	
	
	//create the variables to hold temp output
	vec_d patchValues; init_vec_d(patchValues, 0);
	vec_d temp_function_values; init_vec_d(temp_function_values,0);
	
	vec_d AtimesF;  init_vec_d(AtimesF,0);
    
	vec_d tempvec; init_vec_d(tempvec,0);
	vec_d tempvec2; init_vec_d(tempvec2,0);
	
	
	mat_d temp_jacobian_functions, temp_jacobian_parameters;
	init_mat_d(temp_jacobian_functions,0,0); init_mat_d(temp_jacobian_parameters,0,0);
	
	mat_d AtimesJ; init_mat_d(AtimesJ,1,1); AtimesJ->rows = AtimesJ->cols = 1;
	
	mat_d Jv_jac; init_mat_d(Jv_jac,0,0);
	mat_d Jv_Patch; init_mat_d(Jv_Patch,0,0);
    
	comp_d temp, temp2, temp3;
	comp_d proj_bottom, proj_top, proj_mid;
	
	comp_d u, v;
	mul_d(temp, one_minus_s, BED->u_target);
	mul_d(temp2, pathVars, BED->u_start);
	add_d(u,temp, temp2);
	
	mul_d(temp, one_minus_s, BED->v_target);
	mul_d(temp2, pathVars, BED->v_start);
	add_d(v,temp, temp2);
	
	
	comp_d one_minus_u, one_minus_v;
	set_one_d(one_minus_u);
	sub_d(one_minus_u, one_minus_u, u);
	
	set_one_d(one_minus_v);
	sub_d(one_minus_v, one_minus_v, v);
	
	//initialize some more containers, for the unused stuff from the called evaluators.
	point_d unused_function_values, unused_parVals;
	init_vec_d(unused_function_values,0);init_vec_d(unused_parVals,0);
	vec_d unused_parDer; init_vec_d(unused_parDer,0);
	mat_d unused_Jp; init_mat_d(unused_Jp,0,0);
    
	
	
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
	
	for (ii = 0; ii<BED->num_variables; ii++)
		set_zero_d(&funcVals->coord[ii]);  // initialize entire matrix to 0
	
	// the main evaluations for $x$
	
	BED->mid_memory.set_globals_to_this();
	
    
//	print_point_to_screen_matlab(curr_mid_vars,"curr_mid_vars");
	offset = 0;
	evalProg_d(temp_function_values, parVals, parDer,
               temp_jacobian_functions, unused_Jp, curr_mid_vars, pathVars, BED->SLP_mid);
	
    
	
	// randomize
	mul_mat_vec_d(AtimesF,BED->randomizer_matrix, temp_function_values);
	mat_mul_d(AtimesJ,BED->randomizer_matrix,temp_jacobian_functions);
	
	// for midpoint functions
	for (ii=0; ii<AtimesF->size; ii++)
		set_d(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	// the jacobian equations for midpoint
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_d(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
	
	
	
	BED->bottom_memory.set_globals_to_this();
	
	
	offset = BED->randomizer_matrix->rows; //y0
	int offset_horizontal = BED->num_mid_vars;
	evalProg_d(temp_function_values, parVals, parDer,
               temp_jacobian_functions, unused_Jp, curr_bottom_vars, pathVars, BED->SLP_bottom);
    
    //	std::cout << offset << " " << offset_horizontal << std::endl;
	// randomize
	mul_mat_vec_d(AtimesF,BED->randomizer_matrix_bottom, temp_function_values); // set values of AtimesF (A is randomization matrix)
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real) randomization
		set_d(&funcVals->coord[ii+offset], &AtimesF->coord[ii]);
    
	
	mat_mul_d(AtimesJ,BED->randomizer_matrix_bottom,temp_jacobian_functions);
    //	print_matrix_to_screen_matlab(AtimesJ,"AtimesJy2");
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_d(&Jv->entry[ii+offset][jj+offset_horizontal],&AtimesJ->entry[ii][jj]);
	
	
	
	
	
	offset = BED->randomizer_matrix->rows + BED->randomizer_matrix_bottom->rows; // y2
	offset_horizontal = BED->num_mid_vars + BED->num_bottom_vars;
	
	BED->top_memory.set_globals_to_this();
	
	evalProg_d(temp_function_values, parVals, parDer,
               temp_jacobian_functions, unused_Jp, curr_top_vars, pathVars, BED->SLP_top);
    
	
	// randomize
	mul_mat_vec_d(AtimesF,BED->randomizer_matrix_top, temp_function_values); // set values of AtimesF (A is randomization matrix)
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real) randomization
		set_d(&funcVals->coord[ii+offset], &AtimesF->coord[ii]);
	// the jacobian equations for orig
	//  randomize the original functions and jacobian
	mat_mul_d(AtimesJ,BED->randomizer_matrix_top,temp_jacobian_functions);
	
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_d(&Jv->entry[ii+offset][jj+offset_horizontal],&AtimesJ->entry[ii][jj]);
	
	
	// done with evaluation, so set these to NULL to prevent accidental deletion of the data.
	BED->top_memory.set_globals_null();
	
	
	
	
	
	
	
	
	//
	// now three equations involving the projections, causing movement in the u direction.
	//
	
	
	
	offset = BED->randomizer_matrix->rows + BED->randomizer_matrix_bottom->rows + BED->randomizer_matrix_top->rows;
	
	projection_value_homogeneous_input(proj_mid, curr_mid_vars,BED->pi[0]);
	projection_value_homogeneous_input(proj_top, curr_top_vars,BED->pi[0]);
	projection_value_homogeneous_input(proj_bottom, curr_bottom_vars,BED->pi[0]);
	
    
	
	
	mul_d(temp, one_minus_u, BED->crit_val_left);
	mul_d(temp2, u, BED->crit_val_right);
    
	
	add_d(temp3, temp, temp2); // temp3 = (1-u)*p2(w0) + u*p2(w2)
    
	sub_d(&funcVals->coord[offset+0], proj_mid, temp3);
	sub_d(&funcVals->coord[offset+1], proj_bottom, temp3);
	sub_d(&funcVals->coord[offset+2], proj_top, temp3);
	
	
	// now the derivatives of the supplemental \pi[0] equations
	// d/dt
	sub_d(temp, BED->u_start, BED->u_target);
	mul_d(&Jp->entry[offset][0],BED->crit_val_left, temp);
	
	neg_d(temp, temp);
	mul_d(temp2, BED->crit_val_right, temp);
	add_d(&Jp->entry[offset+0][0],&Jp->entry[offset][0],temp2); // c_l(u_0 - u_t) + c_r(u_t - u_0)
	
	set_d(&Jp->entry[offset+1][0],&Jp->entry[offset][0]);
	set_d(&Jp->entry[offset+2][0],&Jp->entry[offset][0]);
	
	
	
	// d/dx
	offset_horizontal = 0;
	div_d(&Jv->entry[offset][offset_horizontal], proj_mid, &curr_mid_vars->coord[0]);
	neg_d(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -*proj_mid / x[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++)  // mid
		set_d(&Jv->entry[offset][ii], &BED->pi[0]->coord[ii]);
	
	offset_horizontal = BED->num_mid_vars; // y2
	div_d(&Jv->entry[offset+1][offset_horizontal], proj_bottom, &curr_bottom_vars->coord[0]);
	neg_d(&Jv->entry[offset+1][offset_horizontal], &Jv->entry[offset+1][offset_horizontal]); //  -*proj_bottom / y0[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++)
		set_d(&Jv->entry[offset+1][offset_horizontal+ii], &BED->pi[0]->coord[ii]);
	
	offset_horizontal = BED->num_mid_vars + BED->num_bottom_vars; // y2
	div_d(&Jv->entry[offset+2][offset_horizontal], proj_top, &curr_top_vars->coord[0]);
	neg_d(&Jv->entry[offset+2][offset_horizontal], &Jv->entry[offset+2][offset_horizontal]); //  -*proj_mid / y2[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++)
		set_d(&Jv->entry[offset+2][offset_horizontal+ii], &BED->pi[0]->coord[ii]);
	
    
	
	
	
	
	
	
	// finally, the fourth equation, which forces the midpoint to be a midpoint
	
	offset+=3;
	
	projection_value_homogeneous_input(proj_mid, curr_mid_vars,BED->pi[1]); // mid
	projection_value_homogeneous_input(proj_bottom, curr_bottom_vars,BED->pi[1]); // y0
	projection_value_homogeneous_input(proj_top, curr_top_vars,BED->pi[1]); // y2
	
	
	
	//function value
	mul_d(temp, one_minus_v, proj_bottom);
	mul_d(temp2, v, proj_top);
	add_d(temp3, temp, temp2);
	sub_d(&funcVals->coord[offset], proj_mid, temp3);
	
	
	//Jv for the last equation
	// midpoint
	
	offset_horizontal = 0;
	div_d(&Jv->entry[offset][offset_horizontal], proj_mid, &curr_mid_vars->coord[0]);
	neg_d(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -proj_mid / x0[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++)
		div_d(&Jv->entry[offset][ii], &BED->pi[1]->coord[ii], &curr_mid_vars->coord[0]);
	
	
	// y0
	offset_horizontal = BED->num_mid_vars;
	div_d(&Jv->entry[offset][offset_horizontal], proj_bottom, &curr_bottom_vars->coord[0]);
	mul_d(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal], one_minus_v);
    //	neg_d(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -(1-v)*proj_bottom / y0[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++) // only go this far, because all remaining entries are 0
	{
		mul_d(&Jv->entry[offset][ii+offset_horizontal], one_minus_v, &BED->pi[1]->coord[ii]);
		neg_d(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal]);
		div_d(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal], &curr_bottom_vars->coord[0]);
	}
	
	
	// y2
	offset_horizontal = BED->num_mid_vars+BED->num_bottom_vars;
	div_d(&Jv->entry[offset][offset_horizontal], proj_top, &curr_top_vars->coord[0]);
	mul_d(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal], v);
    //	neg_d(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -(v)*proj_top / y2[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++)
	{
		mul_d(&Jv->entry[offset][ii+offset_horizontal], v, &BED->pi[1]->coord[ii]);
		neg_d(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal]);
		div_d(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal], &curr_top_vars->coord[0]);
	}
	
	
    //	the Jp entry for this last equation.
    //should be correct up to perhaps switching left and right proj val
	
	sub_d(temp, BED->v_start, BED->v_target);
	mul_d(&Jp->entry[offset][0], proj_bottom, temp);
	
	neg_d(temp, temp);
	mul_d(temp2, proj_top, temp);
	
	add_d(&Jp->entry[offset][0], &Jp->entry[offset][0], temp2);
	
	offset++;
	
	
	
	
	if (offset != BED->num_variables - BED->patch.num_patches) {
		std::cout << "appear to have offset " << offset << " but should be " << BED->num_variables - BED->patch.num_patches << std::endl;
        //		print_matrix_to_screen_matlab(Jv,"Jv");
		mypause();
	}
    // evaluate the patch
    patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, unused_Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	offset = BED->num_variables - BED->patch.num_patches;
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_d(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
		
		// Jv = Jv_Patch
		for (jj = 0; jj<BED->num_variables; jj++) // for each variable
			set_d(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	
    
	
    
	
	// finally, set parVals & parDer correctly
	
    change_size_point_d(parVals, 1);  change_size_vec_d(parDer, 1);
    parVals->size = parDer->size = 1;
	
    set_d(&parVals->coord[0], pathVars); // s = t
    set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
    
	if ( (BED->verbose_level==14) || (BED->verbose_level == -14)) {
		std::cout << color::blue();
		printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
		std::cout << color::console_default();
		
		print_point_to_screen_matlab(current_variable_values,"curr_vars");
        //		print_point_to_screen_matlab(funcVals,"F_d");
        //		print_matrix_to_screen_matlab(Jv,"Jv");
        //		print_matrix_to_screen_matlab(Jp,"Jp");
		
	}
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_d dehommed; init_vec_d(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,curr_mid_vars);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_d(dehommed);
#endif
	
	
	
	clear_vec_d(curr_mid_vars);
	clear_vec_d(curr_top_vars);
	clear_vec_d(curr_bottom_vars);
	clear_vec_d(patchValues);
	clear_vec_d(temp_function_values);
	
	
	clear_vec_d(AtimesF);
	
	clear_vec_d(tempvec);
	clear_vec_d(tempvec2);
	
	
	clear_mat_d(Jv_Patch);
	
	
	
	clear_mat_d(AtimesJ);
	clear_mat_d(Jv_jac);
	clear_mat_d(temp_jacobian_functions);
	clear_mat_d(temp_jacobian_parameters);
	
	
	
	
	clear_vec_d(unused_function_values);
	clear_vec_d(unused_parVals);
	clear_vec_d(unused_parDer);
	
	
	clear_mat_d(unused_Jp);
	
    
    return 0;
}




int midpoint_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
    
	
    midpoint_eval_data_mp *BED = (midpoint_eval_data_mp *)ED; // to avoid having to cast every time
	
	
	
    int ii, jj;
	int offset;
    
	
	
	
	comp_mp one_minus_s, gamma_s;  init_mp(one_minus_s); init_mp(gamma_s);
	vec_mp curr_mid_vars; init_vec_mp(curr_mid_vars, BED->num_mid_vars);
	vec_mp curr_bottom_vars; init_vec_mp(curr_bottom_vars, BED->num_bottom_vars);
	curr_bottom_vars->size = BED->num_bottom_vars;
	vec_mp curr_top_vars; init_vec_mp(curr_top_vars, BED->num_top_vars);
	curr_top_vars->size = BED->num_top_vars;
	
	//create the variables to hold temp output
	vec_mp patchValues; init_vec_mp(patchValues, 0);
	vec_mp temp_function_values; init_vec_mp(temp_function_values,0);
	
	vec_mp AtimesF;  init_vec_mp(AtimesF,0);
	
	vec_mp tempvec; init_vec_mp(tempvec,0);
	vec_mp tempvec2; init_vec_mp(tempvec2,0);
	
	
	mat_mp temp_jacobian_functions, temp_jacobian_parameters;
	init_mat_mp(temp_jacobian_functions,0,0); init_mat_mp(temp_jacobian_parameters,0,0);
	
	mat_mp AtimesJ; init_mat_mp(AtimesJ,1,1); AtimesJ->rows = AtimesJ->cols = 1;
	
	mat_mp Jv_jac; init_mat_mp(Jv_jac,0,0);
	mat_mp Jv_Patch; init_mat_mp(Jv_Patch,0,0);
	
	
	comp_mp temp, temp2, temp3;  init_mp(temp);
	init_mp(temp2); init_mp(temp3);
	comp_mp proj_bottom, proj_top, proj_mid;
	init_mp(proj_bottom); init_mp(proj_top); init_mp(proj_mid);
	point_mp unused_function_values, unused_parVals;
	comp_mp one_minus_u, one_minus_v;
	comp_mp u, v;  init_mp(u); init_mp(v);
	
	
	
	
	
	
	
	set_one_mp(one_minus_s);
    sub_mp(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
    mul_mp(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	//parse out the variables into proper segments.
	
	curr_mid_vars->size = BED->num_mid_vars;
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_mp(&curr_mid_vars->coord[ii], &current_variable_values->coord[ii]);
	
    
	offset = BED->num_mid_vars;// y0
	
	curr_bottom_vars->size = BED->num_bottom_vars;
	for (ii=0; ii<BED->num_bottom_vars; ii++)
		set_mp(&curr_bottom_vars->coord[ii], &current_variable_values->coord[ii+offset]);
	
	
	offset = BED->num_mid_vars + BED->num_bottom_vars; // y2
    
	curr_top_vars->size = BED->num_top_vars;
	for (ii=0; ii<BED->num_top_vars; ii++)
		set_mp(&curr_top_vars->coord[ii], &current_variable_values->coord[ii+offset]);
	
	
	
	
	
	mul_mp(temp, one_minus_s, BED->u_target);
	mul_mp(temp2, pathVars, BED->u_start);
	add_mp(u,temp, temp2);
	
	mul_mp(temp, one_minus_s, BED->v_target);
	mul_mp(temp2, pathVars, BED->v_start);
	add_mp(v,temp, temp2);
	
	
	
	init_mp(one_minus_u); init_mp(one_minus_v);
	set_one_mp(one_minus_u);
	sub_mp(one_minus_u, one_minus_u, u);
	
	set_one_mp(one_minus_v);
	sub_mp(one_minus_v, one_minus_v, v);
	
	//initialize some more containers, for the unused stuff from the called evaluators.
	
	
	
	init_vec_mp(unused_function_values,0);init_vec_mp(unused_parVals,0);
	vec_mp unused_parDer; init_vec_mp(unused_parDer,0);
	mat_mp unused_Jp; init_mat_mp(unused_Jp,0,0);
	
	
	
	//resize output variables to correct size
	change_size_vec_mp(funcVals,BED->num_variables);
    change_size_mat_mp(Jv, BED->num_variables, BED->num_variables);
    change_size_mat_mp(Jp, BED->num_variables, 1);
	
	
	//////
	// initialize stuff to all 0's
	///////
	
    funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
    Jv->cols = BED->num_variables;  //  <-- this must be square
    Jp->cols = 1;
	
	for (ii=0; ii<Jv->rows; ii++)
		for (jj=0; jj<Jv->cols; jj++)
			set_zero_mp(&Jv->entry[ii][jj]);
	
	for (ii = 0; ii<BED->num_variables; ii++)
		set_zero_mp(&Jp->entry[ii][0]);  // initialize entire matrix to 0
	
	
	// the main evaluations for $x$
	
	BED->mid_memory.set_globals_to_this();
	
	evalProg_mp(temp_function_values, parVals, parDer,
                temp_jacobian_functions, unused_Jp, curr_mid_vars, pathVars, BED->SLP_mid);
	
	
	
	// randomize
	mul_mat_vec_mp(AtimesF,BED->randomizer_matrix, temp_function_values);
	mat_mul_mp(AtimesJ,BED->randomizer_matrix,temp_jacobian_functions);
	
	// for midpoint functions
	for (ii=0; ii<AtimesF->size; ii++)
		set_mp(&funcVals->coord[ii], &AtimesF->coord[ii]);
	
	// the jacobian equations for midpoint
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_mp(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
	
	
	
	BED->bottom_memory.set_globals_to_this();
	
	offset = BED->randomizer_matrix->rows; //y0
	int offset_horizontal = BED->num_mid_vars;
	evalProg_mp(temp_function_values, parVals, parDer,
                temp_jacobian_functions, unused_Jp, curr_bottom_vars, pathVars, BED->SLP_bottom);
	
	//	std::cout << offset << " " << offset_horizontal << std::endl;
	// randomize
	mul_mat_vec_mp(AtimesF,BED->randomizer_matrix_bottom, temp_function_values); // set values of AtimesF (A is randomization matrix)
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real) randomization
		set_mp(&funcVals->coord[ii+offset], &AtimesF->coord[ii]);
	
	
	mat_mul_mp(AtimesJ,BED->randomizer_matrix_bottom,temp_jacobian_functions);
	//	print_matrix_to_screen_matlab(AtimesJ,"AtimesJy2");
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_mp(&Jv->entry[ii+offset][jj+offset_horizontal],&AtimesJ->entry[ii][jj]);
	
	
	
	
	
	offset = BED->randomizer_matrix->rows + BED->randomizer_matrix_bottom->rows; // y2
	offset_horizontal = BED->num_mid_vars + BED->num_bottom_vars;
	
	BED->top_memory.set_globals_to_this();
    
	
	evalProg_mp(temp_function_values, parVals, parDer,
                temp_jacobian_functions, unused_Jp, curr_top_vars, pathVars, BED->SLP_top);
    
	
	// randomize
	mul_mat_vec_mp(AtimesF,BED->randomizer_matrix_top, temp_function_values); // set values of AtimesF (A is randomization matrix)
	for (ii=0; ii<AtimesF->size; ii++)  // for each function, after (real) randomization
		set_mp(&funcVals->coord[ii+offset], &AtimesF->coord[ii]);
	
	// the jacobian equations for orig
	//  randomize the original functions and jacobian
	mat_mul_mp(AtimesJ,BED->randomizer_matrix_top,temp_jacobian_functions);
	for (ii=0; ii< AtimesJ->rows; ii++)
		for (jj=0; jj< AtimesJ->cols; jj++)
			set_mp(&Jv->entry[ii+offset][jj+offset_horizontal],&AtimesJ->entry[ii][jj]);
	
	
	// done with evaluation, so set these to NULL to prevent accidental deletion of the data.
	
	BED->top_memory.set_globals_null();
	
	
	
	
	
	
	
	//
	// now three equations involving the projections, causing movement in the u direction.
	//
	
	
	
	offset = BED->randomizer_matrix->rows + BED->randomizer_matrix_bottom->rows + BED->randomizer_matrix_top->rows;
	
	projection_value_homogeneous_input(proj_mid, curr_mid_vars,BED->pi[0]);
	projection_value_homogeneous_input(proj_top, curr_top_vars,BED->pi[0]);
	projection_value_homogeneous_input(proj_bottom, curr_bottom_vars,BED->pi[0]);
	
    
	
	
	mul_mp(temp, one_minus_u, BED->crit_val_left);
	mul_mp(temp2, u, BED->crit_val_right);
	
	
	add_mp(temp3, temp, temp2); // temp3 = (1-u)*p2(w0) + u*p2(w2)
	
	sub_mp(&funcVals->coord[offset+0], proj_mid, temp3);
	sub_mp(&funcVals->coord[offset+1], proj_bottom, temp3);
	sub_mp(&funcVals->coord[offset+2], proj_top, temp3);
	
	
	// now the derivatives of the supplemental \pi[0] equations
	// d/dt
	sub_mp(temp, BED->u_start, BED->u_target);
	mul_mp(&Jp->entry[offset][0],BED->crit_val_left, temp);
	
	neg_mp(temp, temp);
	mul_mp(temp2, BED->crit_val_right, temp);
	add_mp(&Jp->entry[offset+0][0],&Jp->entry[offset][0],temp2); // c_l(u_0 - u_t) + c_r(u_t - u_0)
	
	set_mp(&Jp->entry[offset+1][0],&Jp->entry[offset][0]);
	set_mp(&Jp->entry[offset+2][0],&Jp->entry[offset][0]);
	
	
	
	// d/dx
	offset_horizontal = 0;
	div_mp(&Jv->entry[offset][offset_horizontal], proj_mid, &curr_mid_vars->coord[0]);
	neg_mp(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -*proj_mid / x[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++){  // mid
		div_mp(&Jv->entry[offset][ii], &BED->pi[0]->coord[ii], &curr_mid_vars->coord[0]);
	}
	
	
	offset_horizontal = BED->num_mid_vars; // y2
	div_mp(&Jv->entry[offset+1][offset_horizontal], proj_bottom, &curr_bottom_vars->coord[0]);
	neg_mp(&Jv->entry[offset+1][offset_horizontal], &Jv->entry[offset+1][offset_horizontal]); //  -*proj_bottom / y0[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++){
		div_mp(&Jv->entry[offset+1][offset_horizontal+ii], &BED->pi[0]->coord[ii], &curr_bottom_vars->coord[0]);
	}
	
	offset_horizontal = BED->num_mid_vars + BED->num_bottom_vars; // y2
	div_mp(&Jv->entry[offset+2][offset_horizontal], proj_top, &curr_top_vars->coord[0]);
	neg_mp(&Jv->entry[offset+2][offset_horizontal], &Jv->entry[offset+2][offset_horizontal]); //  -*proj_mid / y2[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++){
		div_mp(&Jv->entry[offset+2][offset_horizontal+ii], &BED->pi[0]->coord[ii], &curr_top_vars->coord[0]);
	}
	
	
	
	
	
	
	
	// finally, the fourth equation, which forces the midpoint to be a midpoint
	
	offset+=3;
	
	projection_value_homogeneous_input(proj_mid, curr_mid_vars,BED->pi[1]); // mid
	projection_value_homogeneous_input(proj_bottom, curr_bottom_vars,BED->pi[1]); // y0
	projection_value_homogeneous_input(proj_top, curr_top_vars,BED->pi[1]); // y2
	
	
	
	//function value
	mul_mp(temp, one_minus_v, proj_bottom);
	mul_mp(temp2, v, proj_top);
	add_mp(temp3, temp, temp2);
	sub_mp(&funcVals->coord[offset], proj_mid, temp3);
	
	
	//Jv for the last equation
	// midpoint
	
	offset_horizontal = 0;
	div_mp(&Jv->entry[offset][offset_horizontal], proj_mid, &curr_mid_vars->coord[0]);
	neg_mp(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -*proj_mid / x0[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++){
		div_mp(&Jv->entry[offset][ii], &BED->pi[1]->coord[ii], &curr_mid_vars->coord[0]);
	}
	
	
	// y0
	offset_horizontal = BED->num_mid_vars;
	div_mp(&Jv->entry[offset][offset_horizontal], proj_bottom, &curr_bottom_vars->coord[0]);
	mul_mp(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal], one_minus_v);
    //	neg_mp(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -(1-v)*proj_bottom / y0[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++) // only go this far, because all remaining entries are 0
	{
		mul_mp(&Jv->entry[offset][ii+offset_horizontal], one_minus_v, &BED->pi[1]->coord[ii]);
		neg_mp(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal]);
		div_mp(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal], &curr_bottom_vars->coord[0]);
	}
	
	
	// y2
	offset_horizontal = BED->num_mid_vars+BED->num_bottom_vars;
	div_mp(&Jv->entry[offset][offset_horizontal], proj_top, &curr_top_vars->coord[0]);
	mul_mp(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal], v);
    //	neg_mp(&Jv->entry[offset][offset_horizontal], &Jv->entry[offset][offset_horizontal]); //  -(v)*proj_top / y2[0]
	for (int ii=1; ii<BED->num_mid_vars; ii++)
	{
		mul_mp(&Jv->entry[offset][ii+offset_horizontal], v, &BED->pi[1]->coord[ii]);
		neg_mp(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal]);
		div_mp(&Jv->entry[offset][ii+offset_horizontal], &Jv->entry[offset][ii+offset_horizontal], &curr_top_vars->coord[0]);
	}
	
	
	//	the Jp entry for this last equation.
	//should be correct up to perhaps switching left and right proj val
	
	sub_mp(temp, BED->v_start, BED->v_target);
	mul_mp(&Jp->entry[offset][0], proj_bottom, temp);
	
	neg_mp(temp, temp);
	mul_mp(temp2, proj_top, temp);
	
	add_mp(&Jp->entry[offset][0], &Jp->entry[offset][0], temp2);
	
	offset++;
	
	
	
	
	
	if (offset != BED->num_variables - BED->patch.num_patches) {
		std::cout << "appear to have offset " << offset << " but should be " << BED->num_variables - BED->patch.num_patches << std::endl;
		print_matrix_to_screen_matlab(Jv,"Jv");
		mypause();
	}
    // evaluate the patch
    patch_eval_mp(    patchValues, parVals, parDer, Jv_Patch, unused_Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	offset = BED->num_variables - BED->patch.num_patches;
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_mp(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
		
		// Jv = Jv_Patch
		for (jj = 0; jj<BED->num_variables; jj++) // for each variable
			set_mp(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	
	
	
	
	
	// finally, set parVals & parDer correctly
	
    change_size_point_mp(parVals, 1);  change_size_vec_mp(parDer, 1);
    parVals->size = parDer->size = 1;
	
    set_mp(&parVals->coord[0], pathVars); // s = t
    set_one_mp(&parDer->coord[0]);       // ds/dt = 1
	
    
	
	if ( (BED->verbose_level==14) || (BED->verbose_level == -14)) {
		std::cout << color::blue();
		print_comp_matlab(pathVars,"t");
		std::cout << color::console_default();
		
		
		print_point_to_screen_matlab(current_variable_values,"curr_vars");
        //		print_point_to_screen_matlab(funcVals,"F_mp");
        //		print_matrix_to_screen_matlab(Jv,"Jv");
        //		print_matrix_to_screen_matlab(Jp,"Jp");
		
        //		vec_mp result; init_vec_mp(result,0);
        //		dehomogenize(&result, curr_mid_vars);
        //		print_point_to_screen_matlab(result,"mid");
        //
        //		//	std::cout << "\n\n**************\n\n";
        //		clear_vec_mp(result);
        //		mypause();
	}
	
	
	
	
	
	
	
#ifdef printpathnullspace_left
	BED->num_steps++;
	vec_mp dehommed; init_vec_mp(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,curr_mid_vars);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_mp(dehommed);
#endif
    
    
	clear_mp(temp);
	clear_mp(temp2);
	clear_mp(temp3);
	
	clear_mp(proj_bottom);
	clear_mp(proj_top);
	clear_mp(proj_mid);
	
	clear_mp(u);
	clear_mp(v);
	clear_mp(one_minus_u);
	clear_mp(one_minus_v);
	clear_mp(one_minus_s);
	clear_mp(gamma_s);
	
	clear_vec_mp(curr_mid_vars);
	clear_vec_mp(curr_top_vars);
	clear_vec_mp(curr_bottom_vars);
	clear_vec_mp(patchValues);
	clear_vec_mp(temp_function_values);
	
	
	clear_vec_mp(AtimesF);
	
	clear_vec_mp(tempvec);
	clear_vec_mp(tempvec2);
	
	
	clear_mat_mp(Jv_Patch);
	
	
	
	clear_mat_mp(AtimesJ);
	clear_mat_mp(Jv_jac);
	clear_mat_mp(temp_jacobian_functions);
	clear_mat_mp(temp_jacobian_parameters);
	
	
	
	
	clear_vec_mp(unused_function_values);
	clear_vec_mp(unused_parVals);
	clear_vec_mp(unused_parDer);
	
	
	clear_mat_mp(unused_Jp);
	
	
	
	
    return 0;
}






int midpoint_dehom(point_d out_d, point_mp out_mp,
                   int *out_prec,
                   point_d in_d, point_mp in_mp,
                   int in_prec,
                   void const *ED_d, void const *ED_mp)
{
    
    
	
    *out_prec = in_prec;
	
	
	
    if (in_prec < 64)
    { // compute out_d
		midpoint_eval_data_d *BED_d = (midpoint_eval_data_d *)ED_d;
		
		change_size_vec_d(out_d,BED_d->num_mid_vars-1);
		out_d->size = BED_d->num_mid_vars-1;
		
		
		for (int ii=0; ii<BED_d->num_mid_vars-1; ++ii) {
			div_d(&out_d->coord[ii],&in_d->coord[ii+1],&in_d->coord[0]); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
        //
        //		for (int ii=BED_d->num_natural_vars-1; ii<in_d->size-1; ++ii) {
        //			set_d( &out_d->coord[ii],&in_d->coord[ii+1]);
        //		}
        //
		
		
		BED_d = NULL;
		
    }
    else
    { // compute out_mp
		midpoint_eval_data_mp *BED_mp = (midpoint_eval_data_mp *)ED_mp;
		
        
		change_size_vec_mp(out_mp,BED_mp->num_mid_vars-1);
		out_mp->size = BED_mp->num_mid_vars-1;
		
		for (int ii=0; ii<BED_mp->num_mid_vars-1; ++ii) {
			div_mp(&out_mp->coord[ii],&in_mp->coord[ii+1],&in_mp->coord[0]); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		

		
        // set prec on out_mp
        setprec_point_mp(out_mp, *out_prec);
		
		BED_mp = NULL;
		
	}

	
    return 0;
}









int change_midpoint_eval_prec(void const *ED, int new_prec)
{
	midpoint_eval_data_mp *BED = (midpoint_eval_data_mp *)ED; // to avoid having to cast every time
	
    
	BED->SLP_mid->precision = new_prec;
	BED->SLP_bottom->precision = new_prec;
	BED->SLP_top->precision = new_prec;
	
	changePatchPrec_mp(new_prec, &BED->patch);
	
	if (new_prec != BED->curr_prec){
		// change the precision for the patch
		if (BED->verbose_level >=5)
		{
			std::cout << color::brown();
			printf("prec  %d\t-->\t%d\n",BED->curr_prec, new_prec);
			std::cout << color::console_default();
		}
		
		BED->curr_prec = new_prec;
		
		
		
		
		
		
		
		setprec_mp(BED->gamma, new_prec);
		mpf_set_q(BED->gamma->r, BED->gamma_rat[0]);
		mpf_set_q(BED->gamma->i, BED->gamma_rat[1]);
		
		change_prec_mat_mp(BED->randomizer_matrix,new_prec);
		mat_cp_mp(BED->randomizer_matrix,BED->randomizer_matrix_full_prec);
		
		change_prec_mat_mp(BED->randomizer_matrix_bottom,new_prec);
		mat_cp_mp(BED->randomizer_matrix_bottom,BED->randomizer_matrix_bottom_full_prec);
        
		change_prec_mat_mp(BED->randomizer_matrix_top,new_prec);
		mat_cp_mp(BED->randomizer_matrix_top,BED->randomizer_matrix_top_full_prec);
		
		
		for (int ii=0; ii<BED->num_projections; ii++) {
			change_prec_vec_mp(BED->pi[ii],new_prec);
			vec_cp_mp(BED->pi[ii],BED->pi_full_prec[ii]);
		}
		
		setprec_mp(BED->v_target,new_prec);  set_mp(BED->v_target, BED->v_target_full_prec);
		setprec_mp(BED->u_target,new_prec);	set_mp(BED->u_target, BED->u_target_full_prec);
		setprec_mp(BED->v_start,new_prec); set_mp(BED->v_start, BED->v_start_full_prec);
		setprec_mp(BED->u_start,new_prec); set_mp(BED->u_start, BED->u_start_full_prec);
		
		setprec_mp(BED->half, new_prec); set_mp(BED->half, BED->half_full_prec);
		setprec_mp(BED->one, new_prec);  set_mp(BED->one, BED->one_full_prec);
		setprec_mp(BED->zero, new_prec); set_mp(BED->zero, BED->zero_full_prec);
		
		setprec_mp(BED->crit_val_left, new_prec); set_mp(BED->crit_val_left, BED->crit_val_left_full_prec);
		setprec_mp(BED->crit_val_right, new_prec); set_mp(BED->crit_val_right, BED->crit_val_right_full_prec);
        
		
	}
	
	
    return 0;
}






int check_issoln_midpoint_d(endgame_data_t *EG,
                            tracker_config_t *T,
                            void const *ED)
{
    midpoint_eval_data_d *BED = (midpoint_eval_data_d *)ED; // to avoid having to cast every time
	
	
	int ii;
	int offset;
	
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
	
	
	
	
	vec_d f_terminal; init_vec_d(f_terminal, 1); f_terminal->size =1;
	vec_d f_prev; init_vec_d(f_prev, 1); f_prev->size =1;
	
	vec_d curr_mid_vars; init_vec_d(curr_mid_vars, BED->num_mid_vars);
	curr_mid_vars->size = BED->num_mid_vars;
	
	vec_d curr_bottom_vars; init_vec_d(curr_bottom_vars, BED->num_bottom_vars);
	curr_bottom_vars->size = BED->num_bottom_vars;
	
	vec_d curr_top_vars; init_vec_d(curr_top_vars, BED->num_top_vars);
	curr_top_vars->size = BED->num_top_vars;
	
	vec_d temp_function_values; init_vec_d(temp_function_values,1);
	temp_function_values->size = 1;
	
	
	
	vec_d terminal_pt;  init_vec_d(terminal_pt,1); terminal_pt->size = 1;
	vec_d prev_pt;  init_vec_d(prev_pt,1); prev_pt->size = 1;
	
	if (EG->prec>=64){
		vec_mp_to_d(terminal_pt,EG->PD_mp.point);
	}
	else{
		vec_cp_d(terminal_pt,EG->PD_d.point);
	}
	
	
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_d(&curr_mid_vars->coord[ii], &terminal_pt->coord[ii]);
	
	for (ii=0; ii<BED->num_bottom_vars; ii++)
		set_d(&curr_bottom_vars->coord[ii], &terminal_pt->coord[ii+BED->num_mid_vars]);
	
	for (ii=0; ii<BED->num_top_vars; ii++)
		set_d(&curr_top_vars->coord[ii], &terminal_pt->coord[ii+BED->num_mid_vars+BED->num_bottom_vars]);
	
	
	// the main evaluations for $x$
	
	BED->mid_memory.set_globals_to_this();
	
	
	// for midpoint functions
	offset = 0;
	evalProg_d(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_mid_vars, EG->PD_d.time, BED->SLP_mid);
	
	//resize output variables to correct size
	increase_size_vec_d(f_terminal,temp_function_values->size);
	f_terminal->size = temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_d(&f_terminal->coord[ii], &temp_function_values->coord[ii]);
	
	
	BED->bottom_memory.set_globals_to_this();
	
	
	offset += temp_function_values->size; //y0
	evalProg_d(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_bottom_vars, EG->PD_d.time, BED->SLP_bottom);
	
	//resize output variables to correct size
	increase_size_vec_d(f_terminal,f_terminal->size + temp_function_values->size);
	f_terminal->size = f_terminal->size + temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_d(&f_terminal->coord[ii+offset], &temp_function_values->coord[ii]);
	
	
	BED->top_memory.set_globals_to_this();
	
	offset += temp_function_values->size; //y2
	evalProg_d(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_top_vars, EG->PD_d.time, BED->SLP_top);
	
	//resize output variables to correct size
	increase_size_vec_d(f_terminal,f_terminal->size + temp_function_values->size);
	f_terminal->size = f_terminal->size + temp_function_values->size;
	for (ii=0; ii<temp_function_values->size; ii++)
		set_d(&f_terminal->coord[ii+offset], &temp_function_values->coord[ii]);
	
	
	
	
	
	if (EG->last_approx_prec>=64) {
		vec_mp_to_d(prev_pt,EG->last_approx_mp);
	}
	else{
		vec_cp_d(prev_pt, EG->last_approx_d);
	}
	
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_d(&curr_mid_vars->coord[ii], &prev_pt->coord[ii]);
	
	for (ii=0; ii<BED->num_bottom_vars; ii++)
		set_d(&curr_bottom_vars->coord[ii], &prev_pt->coord[ii+BED->num_mid_vars]);
	
	for (ii=0; ii<BED->num_top_vars; ii++)
		set_d(&curr_top_vars->coord[ii], &prev_pt->coord[ii+BED->num_mid_vars+BED->num_bottom_vars]);
	
	
	// the main evaluations for $x$
	
	BED->mid_memory.set_globals_to_this();
	
	
	// for midpoint functions
	offset = 0;
	evalProg_d(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_mid_vars, EG->PD_d.time, BED->SLP_mid);
	
	//resize output variables to correct size
	increase_size_vec_d(f_prev,temp_function_values->size);
	f_prev->size = temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_d(&f_prev->coord[ii], &temp_function_values->coord[ii]);
	
	
	BED->bottom_memory.set_globals_to_this();
	
	
	offset += temp_function_values->size; //y0
	evalProg_d(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_bottom_vars, EG->PD_d.time, BED->SLP_bottom);
	increase_size_vec_d(f_prev, f_prev->size + temp_function_values->size);
	f_prev->size = f_prev->size + temp_function_values->size;
    
	for (ii=0; ii<temp_function_values->size; ii++)
		set_d(&f_prev->coord[ii+offset], &temp_function_values->coord[ii]);
	
	
	BED->top_memory.set_globals_to_this();
	
	offset += temp_function_values->size; //y2
	evalProg_d(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_top_vars, EG->PD_d.time, BED->SLP_top);
	
	increase_size_vec_d(f_prev, f_prev->size + temp_function_values->size);
	f_prev->size = f_prev->size + temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_d(&f_prev->coord[ii+offset], &temp_function_values->coord[ii]);
    
    
	BED->top_memory.set_globals_null();
	
	
	
	
	//	print_point_to_screen_matlab(EG->PD_d.point,"soln");
	//	print_point_to_screen_matlab(e.funcVals,"howfaroff");	// compare the function values
	int isSoln = 1;
	for (ii = 0; (ii < f_terminal->size) && isSoln; ii++)
	{
		n1 = d_abs_d( &f_terminal->coord[ii]); // corresponds to final point
		n2 = d_abs_d( &f_prev->coord[ii]); // corresponds to the previous point
		
		
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
		
		print_point_to_screen_matlab(f_terminal,"terminal");
		print_point_to_screen_matlab(f_prev,"prev");
		
		printf("tol was %le\nmax_rat was %le\n",tol,max_rat);
	}
	
	
	
	
	
	clear_eval_struct_d(e);
	clear_vec_d(f_prev);
	clear_vec_d(f_terminal);
    
	clear_vec_d(prev_pt);
	clear_vec_d(terminal_pt);
	
	
	clear_vec_d(curr_mid_vars);
	clear_vec_d(temp_function_values);
	clear_vec_d(curr_bottom_vars);
	clear_vec_d(curr_top_vars);
	
	
	return isSoln;
	
}


int check_issoln_midpoint_mp(endgame_data_t *EG,
                             tracker_config_t *T,
                             void const *ED)
{
    midpoint_eval_data_mp *BED = (midpoint_eval_data_mp *)ED; // to avoid having to cast every time
	
	int ii;
	int offset;
	
	mpf_t n1, n2, zero_thresh, max_rat;
	mpf_init(n1); mpf_init(n2); mpf_init(zero_thresh); mpf_init(max_rat);
	
	eval_struct_mp e; init_eval_struct_mp(e, 0, 0, 0);
	
	mpf_set_d(max_rat, T->ratioTol);
	
	
	for (ii = 0; ii < T->numVars; ii++)
	{
        if (!(mpfr_number_p(EG->PD_mp.point->coord[ii].r) && mpfr_number_p(EG->PD_mp.point->coord[ii].i)))
		{
			printf("got not a number\n");
			print_point_to_screen_matlab(EG->PD_mp.point,"bad solution");
            return 0;
		}
	}
	
	int num_digits = prec_to_digits((int) mpf_get_default_prec());
	// setup threshold based on given threshold and precision
	if (num_digits > 300)
		num_digits = 300;
	num_digits -= 4;
	double tol = MAX(T->funcResTol, pow(10,-num_digits));
	mpf_set_d(zero_thresh, tol);
	
	
	
	
	
	
	
	vec_mp f_terminal; init_vec_mp(f_terminal, 1); f_terminal->size =1;
	vec_mp f_prev; init_vec_mp(f_prev, 1); f_prev->size =1;
	
	vec_mp curr_mid_vars; init_vec_mp(curr_mid_vars, BED->num_mid_vars);
	curr_mid_vars->size = BED->num_mid_vars;
	
	vec_mp curr_bottom_vars; init_vec_mp(curr_bottom_vars, BED->num_bottom_vars);
	curr_bottom_vars->size = BED->num_bottom_vars;
	
	vec_mp curr_top_vars; init_vec_mp(curr_top_vars, BED->num_top_vars);
	curr_top_vars->size = BED->num_top_vars;
	
	vec_mp temp_function_values; init_vec_mp(temp_function_values,1);
	temp_function_values->size = 1;
	
	
	
	vec_mp terminal_pt;  init_vec_mp(terminal_pt,1); terminal_pt->size = 1;
	vec_mp prev_pt;  init_vec_mp(prev_pt,1); prev_pt->size = 1;
	
	vec_cp_mp(terminal_pt,EG->PD_mp.point);
	
	
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_mp(&curr_mid_vars->coord[ii], &terminal_pt->coord[ii]);
	
	for (ii=0; ii<BED->num_bottom_vars; ii++)
		set_mp(&curr_bottom_vars->coord[ii], &terminal_pt->coord[ii+BED->num_mid_vars]);
	
	for (ii=0; ii<BED->num_top_vars; ii++)
		set_mp(&curr_top_vars->coord[ii], &terminal_pt->coord[ii+BED->num_mid_vars+BED->num_bottom_vars]);
	
	
	// the main evaluations for $x$
	
	BED->mid_memory.set_globals_to_this();
	
	
	
    
	// for midpoint functions
	offset = 0;
	evalProg_mp(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_mid_vars, EG->PD_mp.time, BED->SLP_mid);
	
	//resize output variables to correct size
	increase_size_vec_mp(f_terminal,temp_function_values->size);
	f_terminal->size = temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_mp(&f_terminal->coord[ii], &temp_function_values->coord[ii]);
	
	
	BED->bottom_memory.set_globals_to_this();
	
	offset += temp_function_values->size; //y0
	evalProg_mp(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_bottom_vars, EG->PD_mp.time, BED->SLP_bottom);
	
	//resize output variables to correct size
	increase_size_vec_mp(f_terminal,f_terminal->size + temp_function_values->size);
	f_terminal->size = f_terminal->size + temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_mp(&f_terminal->coord[ii+offset], &temp_function_values->coord[ii]);
	
	
	
	BED->top_memory.set_globals_to_this();
    
	offset += temp_function_values->size; //y2
	evalProg_mp(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_top_vars, EG->PD_mp.time, BED->SLP_top);
	
	//resize output variables to correct size
	increase_size_vec_mp(f_terminal,f_terminal->size + temp_function_values->size);
	f_terminal->size = f_terminal->size + temp_function_values->size;
	for (ii=0; ii<temp_function_values->size; ii++)
		set_mp(&f_terminal->coord[ii+offset], &temp_function_values->coord[ii]);
	
	
	
	
	
	if (EG->last_approx_prec < 64) { // copy to _mp
		vec_d_to_mp(prev_pt,EG->last_approx_d);
	}
	else{
		vec_cp_mp(prev_pt,EG->last_approx_mp);
	}
    
    
    
    
    
    
    
	
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_mp(&curr_mid_vars->coord[ii], &prev_pt->coord[ii]);
	
	for (ii=0; ii<BED->num_bottom_vars; ii++)
		set_mp(&curr_bottom_vars->coord[ii], &prev_pt->coord[ii+BED->num_mid_vars]);
	
	for (ii=0; ii<BED->num_top_vars; ii++)
		set_mp(&curr_top_vars->coord[ii], &prev_pt->coord[ii+BED->num_mid_vars+BED->num_bottom_vars]);
	
	
	// the main evaluations for $x$
	
	BED->mid_memory.set_globals_to_this();
	
	
	// for midpoint functions
	offset = 0;
	evalProg_mp(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_mid_vars, EG->PD_mp.time, BED->SLP_mid);
	
	//resize output variables to correct size
	increase_size_vec_mp(f_prev,temp_function_values->size);
	f_prev->size = temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_mp(&f_prev->coord[ii], &temp_function_values->coord[ii]);
	
	
	BED->bottom_memory.set_globals_to_this();
	
	
	offset += temp_function_values->size; //y0
	evalProg_mp(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_bottom_vars, EG->PD_mp.time, BED->SLP_bottom);
	increase_size_vec_mp(f_prev, f_prev->size + temp_function_values->size);
	f_prev->size = f_prev->size + temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_mp(&f_prev->coord[ii+offset], &temp_function_values->coord[ii]);
	
    
    
	BED->top_memory.set_globals_to_this();
	
	
	offset += temp_function_values->size; //y2
	evalProg_mp(temp_function_values, e.parVals, e.parDer, e.Jv, e.Jp, curr_top_vars, EG->PD_mp.time, BED->SLP_top);
	
	increase_size_vec_mp(f_prev, f_prev->size + temp_function_values->size);
	f_prev->size = f_prev->size + temp_function_values->size;
	
	for (ii=0; ii<temp_function_values->size; ii++)
		set_mp(&f_prev->coord[ii+offset], &temp_function_values->coord[ii]);
	
	
	BED->top_memory.set_globals_null();
    
	
	
	
	
	
	
	
	
	
	// compare the function values
	int isSoln = 1;
	for (ii = 0; ii < f_terminal->size && isSoln; ii++)
	{
		mpf_abs_mp(n1, &f_terminal->coord[ii]);
		mpf_abs_mp(n2, &f_prev->coord[ii]);
		
		//		mpf_out_str(NULL,10,9,n1);
		
		if ( (mpf_cmp(zero_thresh, n1) <= 0) &&  (mpf_cmp(n1, n2) <= 0) )
		{ // compare ratio
			mpf_mul(n2, max_rat, n2);
			if (mpf_cmp(n1, n2) > 0){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (mp) 1\nmax_rat: ");
                mpf_out_str(NULL,10,10,max_rat);
                printf("\n");
			}
		}
		else if ( (mpf_cmp(zero_thresh, n2) <= 0) &&  (mpf_cmp(n2, n1) <= 0) )
		{ // compare ratio
			mpf_mul(n1, max_rat, n1);
			if (mpf_cmp(n2, n1) > 0){
				isSoln = 0;
				printf("labeled as non_soln due to max_rat (mp) 2\nmax_rat: ");
                mpf_out_str(NULL,10,10,max_rat);
                printf("\n");
			}
		}
	}
	
	
    
	
	
	mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);
	
	
	clear_eval_struct_mp(e);
	
	clear_vec_mp(f_terminal);
	clear_vec_mp(f_prev);
	
	clear_vec_mp(prev_pt);
	clear_vec_mp(terminal_pt);
    
	clear_vec_mp(curr_mid_vars);
	clear_vec_mp(temp_function_values);
	clear_vec_mp(curr_bottom_vars);
	clear_vec_mp(curr_top_vars);
	
	

	

	
	
	
	
	return isSoln;
	
}







int check_isstart_midpoint_d(point_d testpoint,
                             tracker_config_t *T,
                             void const *ED)
{
	
	eval_struct_d e;
	init_eval_struct_d(e,0, 0, 0);
	
	comp_d time;
	set_one_d(time);
	
	
	double tol = (1e-9);
	
	midpoint_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, testpoint, time, ED);
	
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



void check_midpoint_evaluator(point_mp current_values,
                              void const *ED)
{
	int ii;
	printf("checking homogeneousness of double evaluator\n");
    midpoint_eval_data_d *BED = (midpoint_eval_data_d *)ED; // to avoid having to cast every time
    //initialize
	eval_struct_d e_d; init_eval_struct_d(e_d, 0, 0, 0);
	eval_struct_d e_d2; init_eval_struct_d(e_d2, 0, 0, 0);
	
	
	
	
	comp_d zerotime; set_zero_d(zerotime);
	
	
	
	
	point_d tempvec;  init_point_d(tempvec,0);
	vec_mp_to_d(tempvec, current_values);
	
	midpoint_eval_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, tempvec, zerotime, ED);
	
	
	comp_d lambda; get_comp_rand_d(lambda);
	
	for (ii=0; ii<BED->num_mid_vars; ii++) {
		mul_d(&tempvec->coord[ii],&tempvec->coord[ii],lambda);
	}
	
	midpoint_eval_d(e_d2.funcVals, e_d2.parVals, e_d2.parDer, e_d2.Jv, e_d2.Jp, tempvec, zerotime, ED);
	
	
	printf("lambda = %lf+1i*%lf\n",lambda->r, lambda->i);
	print_point_to_screen_matlab(e_d.funcVals,"f");
	print_point_to_screen_matlab(e_d2.funcVals,"f2");
	
	
	mypause();
	
	return;
	
}















