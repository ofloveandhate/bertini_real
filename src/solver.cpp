#include "solver.hpp"






void adjust_tracker_AMP(tracker_config_t * T, int num_variables)
{
		T->AMP_eps = (double) num_variables * num_variables;  //According to Demmel (as in the AMP paper), n^2 is a very reasonable bound for \epsilon.
		T->AMP_Phi = T->AMP_bound_on_degree*(T->AMP_bound_on_degree-1.0)*T->AMP_bound_on_abs_vals_of_coeffs;  //Phi from the AMP paper.
		T->AMP_Psi = T->AMP_bound_on_degree*T->AMP_bound_on_abs_vals_of_coeffs;  //Psi from the AMP paper.
																																					// initialize latest_newton_residual_mp to the maximum precision
		mpf_init2(T->latest_newton_residual_mp, T->AMP_max_prec);
}






int solver::send(parallelism_config & mpi_config)
{
	
	
	int *buffer = new int[4];
	
	buffer[0] = this->num_variables;
	buffer[1] = this->num_steps;
	buffer[2] = this->verbose_level;
	buffer[3] = this->MPType;
	
	MPI_Bcast(buffer, 4, MPI_INT, 0, MPI_COMM_WORLD);
	
	send_preproc_data(&this->preProcData);
	
	int num_SLP;
	if (this->have_SLP) 
		num_SLP = 1;
	else
		num_SLP = 0;
	
	MPI_Bcast(&num_SLP, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	for (int ii=0; ii<num_SLP; ii++) {
		std::cout << "master bcasting the SLP, MPType" << this->MPType << std::endl;
		bcast_prog_t(this->SLP, MPType, 0, 0);
	}
	
	delete(buffer);
	return SUCCESSFUL;
}

int solver::receive(parallelism_config & mpi_config)
{
	
	int *buffer = new int[4];
	

	MPI_Bcast(buffer, 4, MPI_INT, 0, MPI_COMM_WORLD);
	
	this->num_variables = buffer[0];
	this->num_steps = buffer[1];
	this->verbose_level = buffer[2];
	this->MPType = buffer[3];
	
	
	receive_preproc_data(&this->preProcData);
	
	if (this->have_SLP) {
		clearProg(this->SLP, this->MPType, 1);
		this->have_SLP = false;
	}
//	this->SLP = (prog_t *) br_malloc(sizeof(prog_t));
	
	
	int num_SLP = 1;
	MPI_Bcast(&num_SLP, 1, MPI_INT, 0, MPI_COMM_WORLD); // get the number of SLP's to receieve
	
	if (num_SLP>0) {
		prog_t * _SLP = (prog_t *) br_malloc(num_SLP*sizeof(prog_t));
		for (int ii=0; ii<num_SLP; ii++) {
			std::cout << "worker bcasting the SLP, MPType" << this->MPType << std::endl;
			bcast_prog_t(&_SLP[ii], this->MPType, 1, 0); // last two arguments are: myid, headnode
			std::cout << "worker copying the SLP" << std::endl;
			this->SLP = &_SLP[ii];
//			cp_prog_t(this->SLP, &_SLP[ii]);
			std::cout << "worker copied the SLP" << std::endl;
		}
		
		this->have_SLP = true;
		initEvalProg(this->MPType);
	}
	
	
	
	
	delete(buffer);
	return SUCCESSFUL;
}




//////////////
//
//	SOLVER MP
//
////////////

int solver_mp::send(parallelism_config & mpi_config)
{
	solver::send(mpi_config);
	
	std::cout << "master sending patch " << std::endl;
	print_matrix_to_screen_matlab(this->patch.patchCoeff,"patchCoeff");
	std::cout << patch.num_patches << " " << patch.curr_prec << std::endl;
	send_patch_mp(&this->patch);
	
	std::cout << "master sent patch_mp" << std::endl;
	
	bcast_comp_mp(this->gamma, 0,0);
	
	int *buffer = new int[2];
	buffer[0] = randomizer_matrix->rows;
	buffer[1] = randomizer_matrix->cols;
	
	MPI_Bcast(buffer, 2, MPI_INT, 0, MPI_COMM_WORLD);
	if (this->MPType==2) {
		bcast_mat_mp(randomizer_matrix_full_prec, 0, 0);
	}
	else{
		print_matrix_to_screen_matlab(randomizer_matrix,"R");
		bcast_mat_mp(randomizer_matrix, 0, 0);
	}
	
	delete[] buffer;
	std::cout << "master sent randomizer matrix" << std::endl;
	buffer = new int[1];
	
	
	buffer[0] = this->curr_prec;
	
	
	MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	delete[] buffer;
	
	return SUCCESSFUL;
}


int solver_mp::receive(parallelism_config & mpi_config)
{
	
	solver::receive(mpi_config);
	
	
	receive_patch_mp(&this->patch);
	
	bcast_comp_mp(this->gamma, 1,0);
	
	int *buffer = new int[2];
	MPI_Bcast(buffer, 2, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (this->MPType==2) {
		init_mat_mp2(randomizer_matrix_full_prec,buffer[0],buffer[1],1024);
		randomizer_matrix_full_prec->rows = buffer[0];
		randomizer_matrix_full_prec->cols = buffer[1];
		
		bcast_mat_mp(randomizer_matrix_full_prec, 1, 0);
		
		init_mat_mp(randomizer_matrix,0,0);
		mat_cp_mp(randomizer_matrix, randomizer_matrix_full_prec);

	}
	else{
		init_mat_mp(randomizer_matrix,buffer[0],buffer[1]);
		randomizer_matrix->rows = buffer[0];
		randomizer_matrix->cols = buffer[1];
		bcast_mat_mp(randomizer_matrix, 1, 0);
	}
	
	delete[] buffer;
	
	
	buffer = new int[1];
	
	MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	this->curr_prec = buffer[0];
	
	delete[] buffer;
	return SUCCESSFUL;
}







//////////////
//
//	SOLVER D
//
////////////

int solver_d::send(parallelism_config & mpi_config)
{
	solver::send(mpi_config);
	
	send_patch_d(&(this->patch));

	bcast_comp_d(this->gamma, 0,0);
	
	bcast_mat_d(randomizer_matrix, 0, 0);
	
	return SUCCESSFUL;
}


int solver_d::receive(parallelism_config & mpi_config)
{
	
	solver::receive(mpi_config);

	receive_patch_d(&this->patch); // the receiving part of the broadcast
	
	bcast_comp_d(this->gamma, 1,0);
	
	bcast_mat_d(randomizer_matrix, 1, 0);
	
	return SUCCESSFUL;
}








//////





void get_tracker_config(solver_configuration & solve_options,int MPType)
{
	
	//necessary for the setupConfig call
	double intrinsicCutoffMultiplier;
	int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, supersetOnly = 0, paramHom = 0;
	//end necessaries for the setupConfig call.
	int constructWitnessSet = 0;
	
  setupConfig(&solve_options.T,
							&solve_options.midpoint_tol,
							&userHom,
							&useRegen,
							&regenStartLevel,
							&maxCodim,
							&specificCodim,
							&solve_options.path_number_modulus,
							&intrinsicCutoffMultiplier,
							&reducedOnly,
							&constructWitnessSet,
							&supersetOnly,
							&paramHom,
							MPType);
	
//	if (solve_options.T.final_tolerance > solve_options.T.real_threshold)
		solve_options.T.real_threshold = 1e-6;
	
	
	cp_tracker_config_t(&solve_options.T_orig,&solve_options.T);
	
	
	return;
}






void solver_clear_config(solver_configuration & options){
	//has no fields which require clearing.
	return;
}






void generic_solver_master(witness_set * W_new, const witness_set & W,
													 solver_d * ED_d, solver_mp * ED_mp,
													 solver_configuration & solve_options)
{
	
	
	
	
  int num_crossings = 0;
	
	W_new->num_variables = W.num_variables;
	W_new->num_synth_vars = W.num_synth_vars;
	
	if (solve_options.complete_witness_set==1){
		
		W_new->cp_patches(W); // copy the patches over from the original witness set
		W_new->cp_names(W);
	}
	
	
	if (solve_options.use_parallel()) {
		
		bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
		
		switch (solve_options.T.MPType) {
			case 1:
				
				ED_mp->send(solve_options);
//				std::cout << "master done sending mp type" << std::endl;
				break;
				
			default:
				
				ED_d->send(solve_options);
//				std::cout << "master done sending double type " << solve_options.T.MPType <<  std::endl;
				break;
		}
	}
	
	
	trackingStats trackCount; init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	post_process_t *endPoints = (post_process_t *)br_malloc(W.num_pts * sizeof(post_process_t)); //overallocate, expecting full
	
	
	// call the file setup function
	FILE *OUT = NULL, *midOUT = NULL;
	
	generic_setup_files(&OUT, "output",
											&midOUT, "midpath_data");
	
	if (solve_options.use_parallel()) {
		
		
		generic_tracker_loop_master(&trackCount, OUT, midOUT,
																W,
																endPoints,
																ED_d, ED_mp,
																solve_options);
	}
	else{
		generic_tracker_loop(&trackCount, OUT, midOUT,
												 W,
												 endPoints,
												 ED_d, ED_mp,
												 solve_options);
	}
	
	
	
	// close the files
	fclose(midOUT);   fclose(OUT);
	
	
	
	
	// check for path crossings
	if (solve_options.use_midpoint_checker==1) {
		midpoint_checker(trackCount.numPoints, solve_options.T.numVars,solve_options.midpoint_tol, &num_crossings);
	}
	
	// post process
	switch (solve_options.T.MPType) {
		case 0:
			BRpostProcessing(endPoints, W_new, trackCount.successes, &ED_d->preProcData, &solve_options.T, solve_options);
			break;
			
		default:
			BRpostProcessing(endPoints, W_new, trackCount.successes, &ED_mp->preProcData, &solve_options.T, solve_options);
			break;
	}
	
  //clear the endpoints here
	
}






/**
 sets the start_pts structure to hold all points in W
 
 \param startPts the value being set.  should be NULL input.
 \param W the witness_set input
 
 */
void generic_set_start_pts(point_data_d ** startPts,
													 const witness_set & W)
{
	int ii; // counters
	
	*startPts = (point_data_d *)br_malloc(W.num_pts * sizeof(point_data_d));
	
	for (ii = 0; ii < W.num_pts; ii++)
	{ // setup startPts[ii]
		init_point_data_d(&(*startPts)[ii], W.num_variables); // also performs initialization on the point inside startPts
		change_size_vec_d((*startPts)[ii].point,W.num_variables);
		(*startPts)[ii].point->size = W.num_variables;
		
		//1 set the coordinates
		vec_mp_to_d((*startPts)[ii].point, W.pts_mp[ii]);
		
		//2 set the start time to 1.
		set_one_d((*startPts)[ii].time);
	}
}



void generic_set_start_pts(point_data_mp ** startPts,
													 const witness_set & W)
{
	int ii; // counters
	
	(*startPts) = (point_data_mp *)br_malloc(W.num_pts * sizeof(point_data_mp));
	
	for (ii = 0; ii < W.num_pts; ii++)
	{ // setup startPts[ii]
		init_point_data_mp(&(*startPts)[ii], W.num_variables); // also performs initialization on the point inside startPts
		change_size_vec_mp((*startPts)[ii].point,W.num_variables);
		(*startPts)[ii].point->size = W.num_variables;
		
		//1 set the coordinates
		vec_cp_mp((*startPts)[ii].point, W.pts_mp[ii]);
		
		//2 set the start time to 1.
		set_one_mp((*startPts)[ii].time);
	}
}






void generic_tracker_loop(trackingStats *trackCount,
										 FILE * OUT, FILE * MIDOUT,
										 const witness_set & W,  // was the startpts file pointer.
										 post_process_t *endPoints,
										 solver_d * ED_d, solver_mp * ED_mp,
										 solver_configuration & solve_options)
{
	
	
	
	

	
	int (*curr_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*curr_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	int (*change_prec)(void const *, int) = NULL;
	int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;
	
	switch (solve_options.T.MPType) {
		case 0:
			curr_eval_d = ED_d->evaluator_function_d;
			curr_eval_mp = ED_d->evaluator_function_mp;
			change_prec = ED_d->precision_changer;
			find_dehom = ED_d->dehomogenizer;
			break;

		default:
			curr_eval_d = ED_mp->evaluator_function_d;
			curr_eval_mp = ED_mp->evaluator_function_mp;
			change_prec = ED_mp->precision_changer;
			find_dehom = ED_mp->dehomogenizer;
			break;
			
	}

	
	
	
	point_data_d *startPts_d = NULL;
	generic_set_start_pts(&startPts_d, W);
	
	point_data_mp *startPts_mp = NULL;
	generic_set_start_pts(&startPts_mp, W);
	
	solve_options.T.endgameOnly = 0;
	
	
  // setup the rest of the structures
	endgame_data_t EG; //this will hold the temp solution data produced for each individual track
	init_endgame_data(&EG, solve_options.T.Precision);
	
	
	
	
	trackCount->numPoints = W.num_pts;
	int solution_counter = 0;
	
	
	
	// track each of the start points
	
	for (int ii = 0; ii < W.num_pts; ii++)
	{
		if ((solve_options.verbose_level>=0) && ((ii%(solve_options.path_number_modulus))==0))
			printf("tracking path %d of %d\n",ii,W.num_pts);
		
		if (solve_options.T.MPType==2) {
			ED_mp->curr_prec = 64;
		}
		
		if (solve_options.robust==true) {
			robust_track_path(solution_counter, &EG,
												&startPts_d[ii], &startPts_mp[ii],
												OUT, MIDOUT,
												solve_options, ED_d, ED_mp,
												curr_eval_d, curr_eval_mp, change_prec, find_dehom);
		}
		else{
		// track the path
			generic_track_path(solution_counter, &EG,
											 &startPts_d[ii], &startPts_mp[ii],
											 OUT, MIDOUT,
											 &solve_options.T, ED_d, ED_mp,
											 curr_eval_d, curr_eval_mp, change_prec, find_dehom);
		
		}
		
		// check to see if it should be sharpened
		if (EG.retVal == 0 && solve_options.T.sharpenDigits > 0)
		{ // use the sharpener for after an endgame
			sharpen_endpoint_endgame(&EG, &solve_options.T, OUT, ED_d, ED_mp, curr_eval_d, curr_eval_mp, change_prec);
		}
		
		
		
		int issoln;
		
		switch (solve_options.T.MPType) {
			case 0:
					issoln = ED_d->is_solution_checker_d(&EG,  &solve_options.T, ED_d);

				break;
				
			default:
				
				if (EG.prec<64){
					issoln = ED_mp->is_solution_checker_d(&EG,  &solve_options.T, ED_d); }
				else {
					issoln = ED_mp->is_solution_checker_mp(&EG, &solve_options.T, ED_mp); }
				break;
		}
		
		
		
		//get the terminal time in double form
		comp_d time_to_compare;
		if (EG.prec < 64) {
			set_d(time_to_compare,EG.PD_d.time);}
		else {
			mp_to_d(time_to_compare, EG.PD_mp.time); }
		
		
		if ((EG.retVal != 0 && time_to_compare->r > solve_options.T.minTrackT) || !issoln) {  // <-- this is the real indicator of failure...
			
			trackCount->failures++;
			
			if (solve_options.verbose_level>=3) {

				printf("\nthere was a path failure nullspace_left tracking witness point %d\nretVal = %d; issoln = %d\n",ii,EG.retVal, issoln);
				
				print_path_retVal_message(EG.retVal);
				
				if (solve_options.verbose_level >= 5) {
					if (EG.prec < 64)
						print_point_to_screen_matlab(EG.PD_d.point,"bad_terminal_point");
					else
						print_point_to_screen_matlab(EG.PD_mp.point,"bad_terminal_point");
				}
			}
		}
		else
		{
			//otherwise converged, but may have still had non-zero retval due to other reasons.
			endgamedata_to_endpoint(&endPoints[solution_counter], &EG);
			trackCount->successes++;
			solution_counter++; // probably this could be eliminated
		}
		
	}// re: for (ii=0; ii<W.num_pts ;ii++)
	
	
	
	//clear the data structures.
  for (int ii = 0; ii >W.num_pts; ii++)
  { // clear startPts[ii]
    clear_point_data_d(&startPts_d[ii]);
		clear_point_data_mp(&startPts_mp[ii]);
  }
  free(startPts_d);
	free(startPts_mp);

	
	
}



void generic_tracker_loop_master(trackingStats *trackCount,
																 FILE * OUT, FILE * MIDOUT,
																 const witness_set & W,  // was the startpts file pointer.
																 post_process_t *endPoints,
																 solver_d * ED_d, solver_mp * ED_mp,
																 solver_configuration & solve_options)
{
	
	
	point_data_d *startPts_d = NULL;
	point_data_mp *startPts_mp = NULL;
	
	switch (solve_options.T.MPType) {
		case 1:
			generic_set_start_pts(&startPts_mp, W); 
			break;
			
		default:
			generic_set_start_pts(&startPts_d, W);
			break;
	}
	
	
	
	
	
	solve_options.T.endgameOnly = 0;
	
	
	// setup the rest of the structures
	endgame_data_t *EG_receives = (endgame_data_t *) br_malloc(1*sizeof(endgame_data_t)); //this will hold the temp solution data produced for each individual track
	init_endgame_data(&EG_receives[0], solve_options.T.Precision);
	
	
	int *indices_outgoing= (int *) br_malloc(sizeof(int));
	int max_outgoing = 1;
	int max_incoming = 1;
	
	
	
	trackCount->numPoints = W.num_pts;
	int solution_counter = 0;
	
	
	
	// track each of the start points
	
	
	int next_index = 0;
	
	int num_packets = 1;
	
	// seed the workers
	for (int ii=1; ii<solve_options.numprocs; ii++) {
		int next_worker = solve_options.activate_next_worker();
//		std::cout << "master sending first packet to worker" << ii << std::endl;
		send_start_points(next_worker, num_packets,
													startPts_d,
													startPts_mp,
													next_index,
													solve_options);
		
	}
	
	

	while (1)
	{
//		std::cout << "master waiting to receive" << std::endl;
		int source = receive_endpoints(trackCount,
											EG_receives, max_incoming,
											solution_counter,
											endPoints,
											ED_d, ED_mp,
											solve_options);
		
		
		
			int next_worker = solve_options.activate_next_worker();
			
//		std::cout << "master sending" << std::endl;
		// this loop currently assumes the numpackets==1
			send_start_points(next_worker, num_packets,
											 startPts_d,
											 startPts_mp,
											 next_index,
											 solve_options);
		
		if (next_index==W.num_pts)
			break;
		
	}// re: for (ii=0; ii<W.num_pts ;ii++)
	
	while (solve_options.have_active()) {
		
		int source = receive_endpoints(trackCount,
																	 EG_receives, max_incoming,
																	 solution_counter,
																	 endPoints,
																	 ED_d, ED_mp,
																	 solve_options);
	}
	
	
	solve_options.send_all_available(0);
	
	//clear the data structures.
	switch (solve_options.T.MPType) {
		case 1:
			for (int ii = 0; ii >W.num_pts; ii++)
				clear_point_data_mp(&startPts_mp[ii]);
			free(startPts_mp);
			break;
			
		default:
			for (int ii = 0; ii >W.num_pts; ii++)
				clear_point_data_d(&startPts_d[ii]);
			free(startPts_d);
			break;
	}
  
}


void send_start_points(int next_worker, int num_packets,
											point_data_d *startPts_d,
											point_data_mp *startPts_mp,
											int & next_index,
											solver_configuration & solve_options)
{
	MPI_Send(&num_packets, 1, MPI_INT, next_worker, NUMPACKETS, MPI_COMM_WORLD);
	
	int *indices_outgoing = new int[num_packets];
	
	for (int ii=0; ii<num_packets; ii++) {
		indices_outgoing[ii] = next_index;
		next_index++;
	}
	

	MPI_Send(indices_outgoing, num_packets, MPI_INT, next_worker, INDICES, MPI_COMM_WORLD);
	
	
	
	for (int ii=indices_outgoing[0]; ii<=indices_outgoing[num_packets-1]; ii++) {
		if (solve_options.T.MPType==1) {
			send_vec_mp( startPts_mp[ii].point, next_worker);
			
		}
		else
		{
			send_vec_d( startPts_d[ii].point, next_worker);
			
		}
	}
	
	delete[] indices_outgoing;
	return;
}


int receive_endpoints(trackingStats *trackCount,
											endgame_data_t *EG_receives, int & max_incoming,
											int & solution_counter,
											post_process_t *endPoints,
											solver_d * ED_d, solver_mp * ED_mp,
											solver_configuration & solve_options)
{

	//now to receive data
	int num_incoming;
	MPI_Status statty_mc_gatty;
	MPI_Recv(&num_incoming, 1, MPI_INT, MPI_ANY_SOURCE, NUMPACKETS, MPI_COMM_WORLD, &statty_mc_gatty);
	
	if (num_incoming > max_incoming) {
		EG_receives = (endgame_data_t *) br_realloc(EG_receives, num_incoming * sizeof(endgame_data_t));
		for (int ii=max_incoming; ii<num_incoming; ii++) {
			init_endgame_data(&EG_receives[ii], solve_options.T.Precision);
		}
		max_incoming = num_incoming;
	}
	
	int incoming_id = send_recv_endgame_data_t(&EG_receives, &num_incoming, solve_options.T.MPType, statty_mc_gatty.MPI_SOURCE, 0); // the trailing 0 indicates receiving
	

	solve_options.deactivate(statty_mc_gatty.MPI_SOURCE);
	
	for (int ii=0; ii<num_incoming; ii++) {
		int issoln;
		
		switch (solve_options.T.MPType) {
			case 0:
				issoln = ED_d->is_solution_checker_d(&EG_receives[ii],  &solve_options.T, ED_d);
				
				break;
				
			default:
				
				if (EG_receives[ii].prec<64){
					issoln = ED_mp->is_solution_checker_d(&EG_receives[ii],  &solve_options.T, ED_d); } // this function call is a reference!
				else {
					issoln = ED_mp->is_solution_checker_mp(&EG_receives[ii], &solve_options.T, ED_mp); } // this function call is a reference!
				break;
		}
		
		
		
		//get the terminal time in double form
		comp_d time_to_compare;
		if (EG_receives[ii].prec < 64) {
			set_d(time_to_compare,EG_receives[ii].PD_d.time);}
		else {
			mp_to_d(time_to_compare, EG_receives[ii].PD_mp.time); }
		
		
		if ((EG_receives[ii].retVal != 0 && time_to_compare->r > solve_options.T.minTrackT) || !issoln) {  // <-- this is the real indicator of failure...
			
			trackCount->failures++;
			
			if (solve_options.verbose_level>=1) {
				printf("\nthere was a path failure nullspace_left tracking witness point %d\nretVal = %d; issoln = %d\n",EG_receives[ii].pathNum, EG_receives[ii].retVal, issoln);
				
				print_path_retVal_message(EG_receives[ii].retVal);
				
				if (solve_options.verbose_level >= 5) {
					if (EG_receives[ii].prec < 64)
						print_point_to_screen_matlab(EG_receives[ii].PD_d.point,"bad_terminal_point");
					else
						print_point_to_screen_matlab(EG_receives[ii].PD_mp.point,"bad_terminal_point");
				}
			}
	
		}
		else
		{
			//otherwise converged, but may have still had non-zero retval due to other reasons.
			endgamedata_to_endpoint(&endPoints[solution_counter], &EG_receives[ii]);
			trackCount->successes++;
			solution_counter++; // probably this could be eliminated
		}
	}
	
	
	
	
	
	return incoming_id;
}

void generic_tracker_loop_worker(trackingStats *trackCount,
																 FILE * OUT, FILE * MIDOUT,
																 solver_d * ED_d, solver_mp * ED_mp,
																 solver_configuration & solve_options)
{
	
	
//	
//	
//	std::cout << "worker" << solve_options.id() << " entering loop" << std::endl;
//	
	int (*curr_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*curr_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	int (*change_prec)(void const *, int) = NULL;
	int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;
	
	switch (solve_options.T.MPType) {
		case 0:
			curr_eval_d = ED_d->evaluator_function_d;
			curr_eval_mp = ED_d->evaluator_function_mp;
			change_prec = ED_d->precision_changer;
			find_dehom = ED_d->dehomogenizer;
			break;
			
		default:
			curr_eval_d = ED_mp->evaluator_function_d;
			curr_eval_mp = ED_mp->evaluator_function_mp;
			change_prec = ED_mp->precision_changer;
			find_dehom = ED_mp->dehomogenizer;
			break;
			
	}
	
	
	
	
	point_data_d *startPts_d;
	point_data_mp *startPts_mp;
	
	switch (solve_options.T.MPType) {
		case 1:
			startPts_mp = (point_data_mp *) br_malloc(1*sizeof(point_data_mp));
			init_point_data_mp(&startPts_mp[0], ED_mp->num_variables);
			break;
			
		default:
			startPts_d = (point_data_d *) br_malloc(1*sizeof(point_data_d));
			init_point_data_d(&startPts_d[0], ED_d->num_variables);
			break;
	}
	
	
  // setup the rest of the structures
	endgame_data_t * EG = (endgame_data_t *) br_malloc(1*sizeof(endgame_data_t)); //this will hold the temp solution data produced for each individual track
	init_endgame_data(&EG[0], solve_options.T.Precision);
	
	int *indices_incoming = (int *) br_malloc(1*sizeof(int));
	
	MPI_Status statty_mc_gatty;
	int max_num_allocated = 1;
	
	int numStartPts = 1;
	
	
	while (1)
	{
//		std::cout << "worker" << solve_options.id() << " receiving work" << std::endl;
		MPI_Recv(&numStartPts, 1, MPI_INT, solve_options.head(), NUMPACKETS, MPI_COMM_WORLD, &statty_mc_gatty);
		// recv next set of start points
		
		if (numStartPts==0) {
			break;
		}
		
		if (numStartPts<max_num_allocated) {
			switch (solve_options.T.MPType) {
				case 1:
					startPts_mp = (point_data_mp *) br_realloc(startPts_mp, numStartPts*sizeof(point_data_mp));
					break;
					
				default:
					startPts_d = (point_data_d *) br_realloc(startPts_d, numStartPts*sizeof(point_data_d));
					break;
			}
			
			indices_incoming = (int *) br_realloc(indices_incoming, numStartPts*sizeof(int));
			
			max_num_allocated = numStartPts;
		}
		
		MPI_Recv(indices_incoming, numStartPts, MPI_INT, solve_options.head(), INDICES, MPI_COMM_WORLD, &statty_mc_gatty);
		
		for (int ii=0; ii<numStartPts; ii++) {
			switch (solve_options.T.MPType) {
				case 1:
					receive_vec_mp(startPts_mp[ii].point, solve_options.head());
					set_one_mp(startPts_mp[ii].time);
					break;
					
				default:
					receive_vec_d(startPts_d[ii].point, solve_options.head());
					set_one_d(startPts_d[ii].time);
					break;
			}
		}

		
		// track each of the start points
		for (int ii = 0; ii < numStartPts; ii++)
		{
			int current_index = indices_incoming[ii];;
			
			if (solve_options.verbose_level>=0 && (current_index%solve_options.path_number_modulus==0) )
				printf("tracking path %d, worker %d\n", current_index, solve_options.id());
			
			
			if (solve_options.T.MPType==2) {
				ED_mp->curr_prec = 64;
			}
			
			
			// track the path
			generic_track_path(indices_incoming[ii], &EG[ii],
												 &startPts_d[ii], &startPts_mp[ii],
												 OUT, MIDOUT,
												 &solve_options.T, ED_d, ED_mp,
												 curr_eval_d, curr_eval_mp, change_prec, find_dehom);
			
			
			
			// check to see if it should be sharpened
			if (EG[ii].retVal == 0 && solve_options.T.sharpenDigits > 0)
			{ // use the sharpener for after an endgame
				sharpen_endpoint_endgame(&EG[ii], &solve_options.T, OUT, ED_d, ED_mp, curr_eval_d, curr_eval_mp, change_prec);
			}
			
		}// re: for (ii=0; ii<W.num_pts ;ii++)

		std::cout << "worker" << solve_options.id() << " waiting to send point " << indices_incoming[0] << std::endl;
		MPI_Send(&numStartPts, 1, MPI_INT, solve_options.head(), NUMPACKETS, MPI_COMM_WORLD);
		send_recv_endgame_data_t(&EG, &numStartPts, solve_options.T.MPType, solve_options.head(), 1);
		
		
	}

	
	switch (solve_options.T.MPType) {
		case 1:
			for (int ii=0; ii<max_num_allocated; ii++) {
				clear_point_data_mp(&startPts_mp[ii]);
				clear_endgame_data(&EG[ii]);
			}
			break;
			
		default:
			for (int ii=0; ii<max_num_allocated; ii++) {
				clear_point_data_d(&startPts_d[ii]);
				clear_endgame_data(&EG[ii]);
			}
			break;
	}
	
	free(indices_incoming);
}





void generic_track_path(int pathNum, endgame_data_t *EG_out,
												point_data_d *Pin, point_data_mp *Pin_mp,
												FILE *OUT, FILE *MIDOUT,
												tracker_config_t *T,
												void const *ED_d, void const *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
{
	
	EG_out->pathNum = pathNum;
	EG_out->codim = 0; // this is ignored
	
		T->first_step_of_path = 1;
		
		if (T->MPType == 2)
		{ // track using AMP
			EG_out->prec = EG_out->last_approx_prec = 52;
			
			EG_out->retVal = endgame_amp(T->endgameNumber, EG_out->pathNum, &EG_out->prec, &EG_out->first_increase, &EG_out->PD_d, &EG_out->PD_mp, &EG_out->last_approx_prec, EG_out->last_approx_d, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
			
			if (EG_out->prec == 52)
			{ // copy over values in double precision
				EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
				EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
				EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
				findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
			}
			else
			{ // make sure that the other MP things are set to the correct precision
				mpf_clear(EG_out->function_residual_mp);
				mpf_init2(EG_out->function_residual_mp, EG_out->prec);
				
				mpf_clear(EG_out->latest_newton_residual_mp);
				mpf_init2(EG_out->latest_newton_residual_mp, EG_out->prec);
				
				mpf_clear(EG_out->t_val_at_latest_sample_point_mp);
				mpf_init2(EG_out->t_val_at_latest_sample_point_mp, EG_out->prec);
				
				mpf_clear(EG_out->error_at_latest_sample_point_mp);
				mpf_init2(EG_out->error_at_latest_sample_point_mp, EG_out->prec);
				
				// copy over the values
				mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
				mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
				mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
				findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
			}
		}
		else if (T->MPType == 0)
		{ // track using double precision
			EG_out->prec = EG_out->last_approx_prec = 52;
			
			EG_out->retVal = endgame_d(T->endgameNumber, EG_out->pathNum, &EG_out->PD_d, EG_out->last_approx_d, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);  // WHERE THE ACTUAL TRACKING HAPPENS
			EG_out->first_increase = 0;
			// copy over values in double precision
			EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
			EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
			EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
			findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
		}
		else if (T->MPType == 1)
		{
			EG_out->pathNum = pathNum;
			EG_out->codim = 0; // zero dimensional - this is ignored
			
			T->first_step_of_path = 1;
			
			// track using MP
			EG_out->retVal = endgame_mp(T->endgameNumber, EG_out->pathNum, &EG_out->PD_mp, EG_out->last_approx_mp, Pin_mp, T, OUT, MIDOUT, ED_mp, eval_func_mp, find_dehom);
			
			
			EG_out->prec = EG_out->last_approx_prec = T->Precision;
			EG_out->first_increase = 0;
			
			// copy over the values
			mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
			mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
			mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
			findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
			
		}
		
	
	return;
}



void robust_track_path(int pathNum, endgame_data_t *EG_out,
												point_data_d *Pin, point_data_mp *Pin_mp,
												FILE *OUT, FILE *MIDOUT,
												solver_configuration & solve_options,
												solver_d *ED_d, solver_mp *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
{
	
//	std::cout << "using robust tracker" << std::endl;
	EG_out->pathNum = pathNum;
	EG_out->codim = 0; // this is ignored
	
	
	tracker_config_t * T = &solve_options.T;
	
	int iterations=0, max_iterations = 10;
	
	solve_options.backup_tracker_config();
	
	
	
	std::map <int,int> setting_increments;

	
	
	EG_out->retVal = -876; // set to bad return value
	while ((iterations<max_iterations) && (EG_out->retVal!=0)) {
		
		// reset a few things here
		
		EG_out->retVal = 0;
		T->first_step_of_path = 1;
		
		if (T->MPType == 2)
		{ // track using AMP
			
			if (solve_options.T.MPType==2) {
				ED_mp->curr_prec = 64; // reset!
			}
			
			
			EG_out->prec = EG_out->last_approx_prec = 52;
			
			EG_out->retVal = endgame_amp(T->endgameNumber, EG_out->pathNum, &EG_out->prec, &EG_out->first_increase, &EG_out->PD_d, &EG_out->PD_mp, &EG_out->last_approx_prec, EG_out->last_approx_d, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
			
			if (EG_out->prec == 52)
			{ // copy over values in double precision
				EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
				EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
				EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
				findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
			}
			else
			{ // make sure that the other MP things are set to the correct precision
				mpf_clear(EG_out->function_residual_mp);
				mpf_init2(EG_out->function_residual_mp, EG_out->prec);
				
				mpf_clear(EG_out->latest_newton_residual_mp);
				mpf_init2(EG_out->latest_newton_residual_mp, EG_out->prec);
				
				mpf_clear(EG_out->t_val_at_latest_sample_point_mp);
				mpf_init2(EG_out->t_val_at_latest_sample_point_mp, EG_out->prec);
				
				mpf_clear(EG_out->error_at_latest_sample_point_mp);
				mpf_init2(EG_out->error_at_latest_sample_point_mp, EG_out->prec);
				
				// copy over the values
				mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
				mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
				mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
				findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
			}
		}
		else if (T->MPType == 0)
		{ // track using double precision
			EG_out->prec = EG_out->last_approx_prec = 52;
			
			EG_out->retVal = endgame_d(T->endgameNumber, EG_out->pathNum, &EG_out->PD_d, EG_out->last_approx_d, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);  // WHERE THE ACTUAL TRACKING HAPPENS
			EG_out->first_increase = 0;
			// copy over values in double precision
			EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
			EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
			EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
			findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
		}
		else if (T->MPType == 1)
		{
			EG_out->pathNum = pathNum;
			EG_out->codim = 0; // zero dimensional - this is ignored
			
			T->first_step_of_path = 1;
			
			// track using MP
			EG_out->retVal = endgame_mp(T->endgameNumber, EG_out->pathNum, &EG_out->PD_mp, EG_out->last_approx_mp, Pin_mp, T, OUT, MIDOUT, ED_mp, eval_func_mp, find_dehom);
			
			
			EG_out->prec = EG_out->last_approx_prec = T->Precision;
			EG_out->first_increase = 0;
			
			// copy over the values
			mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
			mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
			mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
			findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
			
		}
		
		
		
//TODO: PUT the is_solution check right here.
		
		// get how many times we have changed settings due to this type of failure.
		int current_retval_counter = map_lookup_with_default( setting_increments, EG_out->retVal, 0 );
		
		if (!(EG_out->retVal==0 ||   EG_out->retVal==-50 )) {
			std::cout << "solution had non-zero retVal " << EG_out->retVal << " (" << current_retval_counter << ")th occurrence on iteration " << iterations << "." << std::endl;
			switch (EG_out->retVal) {
					
//				case -50:
////					retVal_reached_minTrackT
////					NBHDRADIUS
//					solve_options.T.minTrackT *= 1e-50;
//					break;
				case -100:
					solve_options.T.AMP_max_prec +=256;
					
					
					//this is higher precision needed.
					break;
					
				case 100:
					solve_options.T.AMP_max_prec +=256; // linear increase by 256.
					break;
					
				case -10:
					solve_options.T.maxNumSteps *=2; // double each time
					
					break;
					
				
					
				case -1:
					
					break;
					
				case -3:
					solve_options.T.minStepSizeBeforeEndGame *= 1e-5;
					solve_options.T.minStepSizeDuringEndGame *= 1e-5;
					solve_options.T.minStepSize *=  1e-5;
					break;
					
				case -4:
					if (current_retval_counter<3) {
						solve_options.T.securityMaxNorm *= 10;  // exponential increase by 10's
						std::cout << "increasing securityMaxNorm to " << solve_options.T.securityMaxNorm << std::endl;
					}
					else
					{
						// on the third try, go to security level 1.
						solve_options.T.securityLevel = 1; // just turn on security level 1 
						std::cout << "setting securityLevel to 1" << std::endl;
					}
						
					//	solve_options.T.final_tolerance = 0.1*solve_options.T.final_tolerance;
					
					break;
					
				default:
					solve_options.T.odePredictor ++;
					break;
			}
		}
		
		// increment the counter for how many times we have changed this type of setting.
		setting_increments[EG_out->retVal] = current_retval_counter + 1;
		
		
		iterations++;
	} // re: while
	
	
	solve_options.reset_tracker_config();
	
	return;
} // re: robust_track_path





// this is to be deprecated and removed shortly.  next time you see this, do it.
void generic_track_path_d(int pathNum, endgame_data_t *EG_out,
													point_data_d *Pin,
													FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
													void const *ED_d, void const *ED_mp,
													int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
													int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
													int (*change_prec)(void const *, int),
													int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
{
	
	
	
	EG_out->pathNum = pathNum;
	EG_out->codim = 0; // this is ignored
	
	T->first_step_of_path = 1;
	
	if (T->MPType == 2)
	{ // track using AMP
		EG_out->prec = EG_out->last_approx_prec = 52;
		
		EG_out->retVal = endgame_amp(T->endgameNumber, EG_out->pathNum, &EG_out->prec, &EG_out->first_increase, &EG_out->PD_d, &EG_out->PD_mp, &EG_out->last_approx_prec, EG_out->last_approx_d, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED_d, ED_mp, eval_func_d, eval_func_mp, change_prec, find_dehom);
		
		if (EG_out->prec == 52)
		{ // copy over values in double precision
			EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
			EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
			EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
			findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
		}
		else
		{ // make sure that the other MP things are set to the correct precision
			mpf_clear(EG_out->function_residual_mp);
			mpf_init2(EG_out->function_residual_mp, EG_out->prec);
			
			mpf_clear(EG_out->latest_newton_residual_mp);
			mpf_init2(EG_out->latest_newton_residual_mp, EG_out->prec);
			
			mpf_clear(EG_out->t_val_at_latest_sample_point_mp);
			mpf_init2(EG_out->t_val_at_latest_sample_point_mp, EG_out->prec);
			
			mpf_clear(EG_out->error_at_latest_sample_point_mp);
			mpf_init2(EG_out->error_at_latest_sample_point_mp, EG_out->prec);
			
			// copy over the values
			mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
			mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
			mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
			findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED_mp, eval_func_mp);
		}
	}
	else if (T->MPType == 0)
	{ // track using double precision
		EG_out->prec = EG_out->last_approx_prec = 52;
		
		EG_out->retVal = endgame_d(T->endgameNumber, EG_out->pathNum, &EG_out->PD_d, EG_out->last_approx_d, Pin, T, OUT, MIDOUT, ED_d, eval_func_d, find_dehom);  // WHERE THE ACTUAL TRACKING HAPPENS
		EG_out->first_increase = 0;
		// copy over values in double precision
		EG_out->latest_newton_residual_d = T->latest_newton_residual_d;
		EG_out->t_val_at_latest_sample_point_d = T->t_val_at_latest_sample_point;
		EG_out->error_at_latest_sample_point_d = T->error_at_latest_sample_point;
		findFunctionResidual_conditionNumber_d(&EG_out->function_residual_d, &EG_out->condition_number, &EG_out->PD_d, ED_d, eval_func_d);
	}
	
	
	
	return;
}



// this is to be deprecated and removed shortly.  next time you see this, do it.
void generic_track_path_mp(int pathNum, endgame_data_t *EG_out,
													 point_data_mp *Pin,
													 FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
													 void const *ED,
													 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
													 int (*change_prec)(void const *, int),
													 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
{
	
	EG_out->pathNum = pathNum;
  EG_out->codim = 0; // zero dimensional - this is ignored
	
  T->first_step_of_path = 1;
	
  // track using MP
  EG_out->retVal = endgame_mp(T->endgameNumber, EG_out->pathNum, &EG_out->PD_mp, EG_out->last_approx_mp, Pin, T, OUT, MIDOUT, ED, eval_func_mp, find_dehom);
	
	
  EG_out->prec = EG_out->last_approx_prec = T->Precision;
  EG_out->first_increase = 0;
	
  // copy over the values
  mpf_set(EG_out->latest_newton_residual_mp, T->latest_newton_residual_mp);
  mpf_set_d(EG_out->t_val_at_latest_sample_point_mp, T->t_val_at_latest_sample_point);
  mpf_set_d(EG_out->error_at_latest_sample_point_mp, T->error_at_latest_sample_point);
  findFunctionResidual_conditionNumber_mp(EG_out->function_residual_mp, &EG_out->condition_number, &EG_out->PD_mp, ED, eval_func_mp);
	
  return;
}





void generic_setup_patch(patch_eval_data_d *P, const witness_set & W)
{
	int ii;
	P->num_patches = W.num_patches;
	init_mat_d(P->patchCoeff, W.num_patches, W.num_variables);
	P->patchCoeff->rows = W.num_patches; P->patchCoeff->cols = W.num_variables;
	
	int varcounter = 0;
	for (int jj=0; jj<W.num_patches; jj++) {
		for (ii=0; ii<varcounter; ii++) {
			set_zero_d(&P->patchCoeff->entry[jj][ii]);
		}
		
		int offset = varcounter;
		for (ii = 0; ii < W.patch_mp[jj]->size ; ii++){
			mp_to_d(&P->patchCoeff->entry[jj][ii+offset],&W.patch_mp[jj]->coord[ii]);
			varcounter++;
		}
		
		for (ii=varcounter; ii<W.num_variables; ii++) {
			set_zero_d(&P->patchCoeff->entry[jj][ii]);
		}
	}
}


void generic_setup_patch(patch_eval_data_mp *P, const witness_set & W)
{
	int ii;
	init_mat_rat(P->patchCoeff_rat, W.num_patches, W.num_variables);
	
	init_mat_mp2(P->patchCoeff, W.num_patches, W.num_variables, mpf_get_default_prec());
	
	P->curr_prec = mpf_get_default_prec();
	P->num_patches = W.num_patches;
	P->patchCoeff->rows = W.num_patches;
	P->patchCoeff->cols = W.num_variables;
	

	
	
	if (W.num_patches==0) {
		std::cerr << "the number of patches in input W is 0.  this is not allowed, the number must be positive.\n" << std::endl;
		deliberate_segfault();
	}
	
//	for (int jj=0; jj<W.num_patches; jj++) {
//		print_point_to_screen_matlab(W.patch_mp[jj],"patch_jj");
//	}
//	std::cout << "W.num_variables " << W.num_variables << std::endl;
//	print_matrix_to_screen_matlab(P->patchCoeff, "patch_mat");
	
	
	int varcounter = 0;
	for (int jj=0; jj<W.num_patches; jj++) {
		
		for (ii=0; ii<varcounter; ii++) {
			set_zero_mp(&P->patchCoeff->entry[jj][ii]);
			mp_to_rat(P->patchCoeff_rat[jj][ii], &P->patchCoeff->entry[jj][ii]);
		}
		
		int offset = varcounter;
		for (ii = 0; ii < W.patch_mp[jj]->size; ii++){
			set_mp(&P->patchCoeff->entry[jj][ii+offset],&W.patch_mp[jj]->coord[ii]);
			mp_to_rat(P->patchCoeff_rat[jj][ii+offset],
								&P->patchCoeff->entry[jj][ii+offset]);
			varcounter++;
		}
		
		for (ii=varcounter; ii<W.num_variables; ii++) {
			set_zero_mp(&P->patchCoeff->entry[jj][ii]);
			mp_to_rat(P->patchCoeff_rat[jj][ii], &P->patchCoeff->entry[jj][ii]);
		}
	}
}



int generic_setup_files(FILE ** OUT, boost::filesystem::path outname,
												FILE ** MIDOUT, boost::filesystem::path midname)
{
	
	*OUT = safe_fopen_write(outname);  // open the main output files.
  *MIDOUT = safe_fopen_write(midname);
	
	return SUCCESSFUL;
}

void get_projection(vec_mp *pi,
										BR_configuration program_options,
										const solver_configuration & solve_options,
										int num_vars,
										int num_projections)
{
	
	int ii,jj;
	for (ii=0; ii<num_projections; ii++) {
		change_size_vec_mp(pi[ii], num_vars);  pi[ii]->size = num_vars;
	}
	
	
	
	//assumes the vector pi is already initialized
	if (program_options.user_projection==1) {
		FILE *IN = safe_fopen_read(program_options.projection_filename.c_str()); // we are already assured this file exists, but safe fopen anyway.
		int tmp_num_vars;
		fscanf(IN,"%d",&tmp_num_vars); scanRestOfLine(IN);
		if (tmp_num_vars!=num_vars-1) {
			printf("the number of variables appearing in the projection\nis not equal to the number of non-homogeneous variables in the problem\n");
			printf("please modify file to have %d coordinate pairs.\n",num_vars-1);
			abort();
		}
		
		for (ii=0; ii<num_projections; ii++) {
			set_zero_mp(&pi[ii]->coord[0]);
			for (jj=1; jj<num_vars; jj++) {
				mpf_inp_str(pi[ii]->coord[jj].r, IN, 10);
				mpf_inp_str(pi[ii]->coord[jj].i, IN, 10);
				scanRestOfLine(IN);
			}
		}
		fclose(IN);
	}
	else{
//		for (ii=0; ii<num_projections; ii++) {
//			set_zero_mp(&pi[ii]->coord[0]);
//			for (jj=1; jj<num_vars; jj++)
//				get_comp_rand_real_mp(&pi[ii]->coord[jj]);//, &temp_getter->entry[ii][jj-1]);
//
//		}
		mat_mp temp_getter;
		init_mat_mp2(temp_getter,0,0,1024);
		make_matrix_random_real_mp(temp_getter,num_projections, num_vars-1, 1024); // this matrix is ~orthogonal
		
		for (ii=0; ii<num_projections; ii++) {
			set_zero_mp(&pi[ii]->coord[0]);
			for (jj=1; jj<num_vars; jj++)
				set_mp(&pi[ii]->coord[jj], &temp_getter->entry[ii][jj-1]);
			
		}

		clear_mat_mp(temp_getter);
	}
	
	return;
}






