#include "solver_midpoint_tracker.hpp"
#include "surface.hpp"



void midpoint_config::setup(const surface_decomposition & surf)
{
	//copy the projections
	for (int ii=0; ii< surf.dimension; ii++)
		add_projection(surf.pi[ii]);
	
	
}






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
	
	
	if (this->MPType==2) {

	}
	
	
}


int midpoint_eval_data_mp::send(parallelism_config & mpi_config)
{
	
	int solver_choice = MIDPOINT;
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
	int *buffer = new int[12];
	MPI_Bcast(buffer, 1, MPI_INT, 0, MPI_COMM_WORLD);
	
	if (buffer[0] != MIDPOINT) {
		mpi_config.abort(777);
	}
	
	solver_mp::receive(mpi_config);
	
	// now can actually receive the data from whoever.
	MPI_Bcast(buffer,12,MPI_INT, 0, mpi_config.my_communicator);
	
	
	
	delete[] buffer;
		
	
	return SUCCESSFUL;
}

int midpoint_eval_data_mp::setup(prog_t * _SLP,
																 midpoint_config *ns_config,
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
	
	mat_cp_mp(randomizer_matrix,
						ns_config->randomizer_matrix);
	
	if (this->MPType==2) {
		mat_cp_mp(randomizer_matrix_full_prec, ns_config->randomizer_matrix);
		
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
		this->BED_mp = new midpoint_eval_data_mp(2);
	else
		this->BED_mp = NULL;
	
	this->is_solution_checker_d = &check_issoln_midpoint_d;
	this->is_solution_checker_mp = &check_issoln_midpoint_mp;
	this->evaluator_function_d = &midpoint_eval_d;
	this->evaluator_function_mp = &midpoint_eval_mp;
	this->precision_changer = &change_midpoint_eval_prec;
	this->dehomogenizer = &midpoint_dehom;
	

	
	
}


int midpoint_eval_data_d::send(parallelism_config & mpi_config)
{
	
	int solver_choice = MIDPOINT;
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
	
	if (buffer[0] != MIDPOINT){
		std::cout << "worker failed to confirm it is receiving the nullspace type eval data" << std::endl;
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




int midpoint_eval_data_d::setup(prog_t * _SLP,
																		midpoint_config *ns_config,
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
	
	
	
	mat_mp_to_d(randomizer_matrix,
							ns_config->randomizer_matrix);
	
	
	if (this->MPType==2)
	{
		this->BED_mp->setup(_SLP,ns_config, W, solve_options);
		rat_to_d(this->gamma, this->BED_mp->gamma_rat);
	}
	
	return SUCCESSFUL;
}

///////////////
//
//   end midpoint_eval_data_d
//
/////////////






















int midpoint_solver_master_entry_point(int										MPType,
																					 witness_set						&W, // carries with it the start points, and the linears.
																					 witness_set						*W_new, // new data goes in here
																					 midpoint_config				*ns_config,
																					 solver_configuration		& solve_options)
{
	
	
	if (solve_options.use_parallel()) {
		solve_options.call_for_help(MIDPOINT);
	}
	
	W_new->num_variables = W.num_variables;
	W_new->num_synth_vars = W.num_synth_vars;
	
	if (solve_options.complete_witness_set==1){
		
		W_new->cp_patches(W); // copy the patches over from the original witness set
		W_new->cp_names(W);
	}
	
	
	
  int num_crossings = 0;
	
  trackingStats trackCount; init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	
	//	solve_options.T.numVars = setupProg(this->SLP, solve_options.T.Precision, solve_options.T.MPType);
	
	int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
	
	prog_t SLP;
	//	// setup a straight-line program, using the file(s) created by the parser
  solve_options.T.numVars = setupProg_count(&SLP, solve_options.T.Precision, solve_options.T.MPType,
																						&startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv,
																						&subFuncsBelow);
	
	midpoint_eval_data_d *ED_d = NULL;
	midpoint_eval_data_mp *ED_mp = NULL;
	
	
	switch (solve_options.T.MPType) {
		case 0:
			ED_d = new midpoint_eval_data_d(0);
			
			ED_d->setup(&SLP,
									ns_config,
									W,
									solve_options);
			break;
			
		case 1:
			ED_mp = new midpoint_eval_data_mp(1);
			
			ED_mp->setup(&SLP,
									 ns_config,
									 W,
									 solve_options);
			// initialize latest_newton_residual_mp
			mpf_init(solve_options.T.latest_newton_residual_mp);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
			break;
		case 2:
			ED_d = new midpoint_eval_data_d(2);
			
			ED_mp = ED_d->BED_mp;
			
			
			ED_d->setup(&SLP,
									ns_config,
									W,
									solve_options);
			
			
			
			adjust_tracker_AMP(& (solve_options.T), W.num_variables);
			// initialize latest_newton_residual_mp
			break;
		default:
			break;
	}
	
	
	
	
	if (solve_options.use_parallel()) {
		
		bcast_tracker_config_t(&solve_options.T, solve_options.id(), solve_options.head() );
		
		switch (solve_options.T.MPType) {
			case 1:
				std::cout << "master sending mp type" << std::endl;
				ED_mp->send(solve_options);
				break;
				
			default:
				std::cout << "master sending double type " << solve_options.T.MPType <<  std::endl;
				ED_d->send(solve_options);
				break;
		}
	}
	
	
	
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
	
	
	
  //clear the endopints here
	
	
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
	

	generic_tracker_loop_worker(&trackCount, OUT, midOUT,
															ED_d, ED_mp,
															solve_options);
	
	
	// close the files
	fclose(midOUT);   fclose(OUT);
	
	
	//clear data
}






int midpoint_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	
	
	
  midpoint_eval_data_d *BED = (midpoint_eval_data_d *)ED; // to avoid having to cast every time
	
	
	
  int ii, jj, kk, mm;
	int offset;
  comp_d one_minus_s, gamma_s;
	
	set_one_d(one_minus_s);
  sub_d(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
  mul_d(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	vec_d curr_x_vars; init_vec_d(curr_x_vars, BED->num_mid_vars);
	curr_x_vars->size = BED->num_mid_vars;
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_d(&curr_x_vars->coord[ii], &current_variable_values->coord[ii]);
	
	offset = BED->num_mid_vars;
	vec_d curr_y_vars; init_vec_d(curr_y_vars, BED->num_y_vars);
	curr_y_vars->size = BED->num_y_vars;
	for (ii=0; ii<BED->num_y_vars; ii++)
		set_d(&curr_y_vars->coord[ii], &current_variable_values->coord[ii+offset]);
	
	offset = BED->num_mid_vars + BED->num_y_vars;
	vec_d curr_z_vars; init_vec_d(curr_z_vars, BED->num_z_vars);
	curr_z_vars->size = BED->num_z_vars;
	for (ii=0; ii<BED->num_z_vars; ii++)
		set_d(&curr_z_vars->coord[ii], &current_variable_values->coord[ii+offset]);
	
	
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
	


	//initialize some containers, for the unused stuff from the called evaluators.
	point_d unused_function_values, unused_parVals;
	init_vec_d(unused_function_values,0);init_vec_d(unused_parVals,0);
	vec_d unused_parDer; init_vec_d(unused_parDer,0);
	mat_d unused_Jp; init_mat_d(unused_Jp,0,0);

	
	// the main evaluations for $x$
	evalProg_d(temp_function_values, parVals, parDer, temp_jacobian_functions, unused_Jp, curr_x_vars, pathVars, BED->SLP);
	
	
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
	
	
	
	
	
	
	
	

	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		set_d(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
		
		// Jv = Jv_Patch
		for (jj = 0; jj<BED->num_mid_vars; jj++) // for each variable
			set_d(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
	}
	
	

	

	
	// finally, set parVals & parDer correctly
	
  change_size_point_d(parVals, 1);  change_size_vec_d(parDer, 1);
  parVals->size = parDer->size = 1;
	
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	if (BED->verbose_level>=5) {
		printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
		//
		//
		//	print_matrix_to_screen_matlab( AtimesJ,"jac");
		//	print_point_to_screen_matlab(curr_x_vars,"currxvars");
		//	print_point_to_screen_matlab(current_variable_values,"curr_vars");
		//	print_point_to_screen_matlab(funcVals,"F");
		//	print_matrix_to_screen_matlab(Jv,"Jv");
		//	print_matrix_to_screen_matlab(Jp,"Jp");

		
		//	std::cout << "\n\n**************\n\n";
		if (BED->verbose_level==10)
			mypause();
	}
	
	
	
	clear_vec_d(curr_x_vars);
	clear_vec_d(curr_y_vars);
	clear_vec_d(curr_z_vars);
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




int midpoint_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
	printf("entering eval_mp\n");
  midpoint_eval_data_mp *BED = (midpoint_eval_data_mp *)ED; // to avoid having to cast every time
	
  	
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
		
		comp_d denom;
		change_size_vec_d(out_d,in_d->size-1);
		out_d->size = in_d->size-1;
		
		set_d(denom, &in_d->coord[0]);
		
		for (int ii=0; ii<BED_d->num_mid_vars-1; ++ii) {
			set_d(&out_d->coord[ii],&in_d->coord[ii+1]);
			div_d(&out_d->coord[ii],&out_d->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
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
		
		comp_mp denom; init_mp(denom);
		change_size_vec_mp(out_mp,in_mp->size-1);
		out_mp->size = in_mp->size-1;
		
		set_mp(denom, &in_mp->coord[0]);
		
		for (int ii=0; ii<BED_mp->num_mid_vars-1; ++ii) {
			set_mp(&out_mp->coord[ii],&in_mp->coord[ii+1]);
			div_mp(&out_mp->coord[ii],&out_mp->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
		}
		
//		for (int ii=BED_mp->num_natural_vars-1; ii<in_mp->size-1; ++ii) {
//			set_mp( &out_mp->coord[ii],&in_mp->coord[ii+1]);
//		}
		
		clear_mp(denom);
		
		
    // set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);
		
		BED_mp = NULL;
		
	}
	
	
//	if (in_prec < 64) {
//		print_point_to_screen_matlab(in_d,"in");
//		print_point_to_screen_matlab(out_d,"out");
//	}
//	else{
//		print_point_to_screen_matlab(in_mp,"in");
//		print_point_to_screen_matlab(out_mp,"out");
//	}
//	mypause();

	
  
	
  return 0;
}









int change_midpoint_eval_prec(void const *ED, int new_prec)
{
	midpoint_eval_data_mp *BED = (midpoint_eval_data_mp *)ED; // to avoid having to cast every time
	
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
		
				
	}
	
	
  return 0;
}






int check_issoln_midpoint_d(endgame_data_t *EG,
																tracker_config_t *T,
																void const *ED)
{
  midpoint_eval_data_d *BED = (midpoint_eval_data_d *)ED; // to avoid having to cast every time
	
	
	int ii;
	

	
	vec_d curr_x_vars; init_vec_d(curr_x_vars, BED->num_mid_vars);
	curr_x_vars->size = BED->num_mid_vars;
	
	vec_d curr_y_vars; init_vec_d(curr_y_vars, BED->num_y_vars);
	curr_y_vars->size = BED->num_y_vars;
	
	
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
		
		
		for (ii=0; ii<BED->num_mid_vars; ii++)
			set_d(&curr_x_vars->coord[ii], &terminal_pt->coord[ii]);
		
		for (ii=0; ii<BED->num_y_vars; ii++)
			set_d(&curr_y_vars->coord[ii], &terminal_pt->coord[ii+BED->num_mid_vars]);
		
		clear_vec_d(terminal_pt);
	}
	else{
		evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_d.point, EG->PD_d.time, BED->SLP);
		//		lin_to_lin_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_d.point, EG->PD_d.time, ED);
		
		for (ii=0; ii<BED->num_mid_vars; ii++)
			set_d(&curr_x_vars->coord[ii], &EG->PD_d.point->coord[ii]);
		
		for (ii=0; ii<BED->num_y_vars; ii++)
			set_d(&curr_y_vars->coord[ii], &EG->PD_d.point->coord[ii+BED->num_mid_vars]);
		
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
	
	
	
	
	
	clear_eval_struct_d(e);
	clear_vec_d(f);
	

	return isSoln;
	
}


int check_issoln_midpoint_mp(endgame_data_t *EG,
																 tracker_config_t *T,
																 void const *ED)
{
  midpoint_eval_data_mp *BED = (midpoint_eval_data_mp *)ED; // to avoid having to cast every time
	
	int ii;
	
	
	
	vec_mp curr_x_vars; init_vec_mp(curr_x_vars, BED->num_mid_vars);
	curr_x_vars->size = BED->num_mid_vars;
	
	vec_mp curr_y_vars; init_vec_mp(curr_y_vars, BED->num_y_vars);
	curr_y_vars->size = BED->num_y_vars;
	
	vec_mp curr_z_vars; init_vec_mp(curr_z_vars, BED->num_z_vars);
	curr_z_vars->size = BED->num_z_vars;
	
	
	
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
	
	
	for (ii=0; ii<BED->num_mid_vars; ii++)
		set_mp(&curr_x_vars->coord[ii], &EG->PD_mp.point->coord[ii]);
	
	for (ii=0; ii<BED->num_y_vars; ii++)
		set_mp(&curr_y_vars->coord[ii], &EG->PD_mp.point->coord[ii+BED->num_mid_vars]);
	
	
	
	//this one guaranteed by entry condition
	evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, BED->SLP);
	

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
	
	
		
	
	
	mpf_clear(n1); mpf_clear(n2); mpf_clear(zero_thresh); mpf_clear(max_rat);
	
	
	clear_eval_struct_mp(e);
	clear_vec_mp(f);
	

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















