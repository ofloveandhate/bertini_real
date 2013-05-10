#include "lintolinSolver.h"

//the main wrapper function for chaining into the lin_to_lin solver.
// pass in:
//     • the MPType, witness_set with the linear and points to move FROM (as well as the patch equation),
//     • max_precision randomizer matrix and linears,
//	   • initialized W_new
//
// you will get out: W_new populated in both mp and double types, with the linears we moved TO, and the new points,
//  in corresponding order.

//solve_options has both a tracker_config_t and a preproc_data.
int lin_to_lin_solver_main(int MPType,
													 witness_set W,
													 mat_mp n_minusone_randomizer_matrix_full_prec,
													 vec_mp *new_linears_full_prec, int num_new_linears,
													 witness_set *W_new,
													 solver_configuration *solve_options){
	
	cp_patches(W_new,W); // copy the patches over from the original witness set
	W_new->num_variables = W.num_variables;
	W_new->MPType = MPType;
	cp_names(W_new,W);
	
	if (num_new_linears==0) {
		W_new->num_linears = 0;
		printf("\nno new linears at which to solve.  returning out of lin_to_lin\n");
		return 0;
	}
	
	if (W.num_linears==0) {
		printf("input witness set had 0 linears!\n");
		exit(01);
	}
	
	W_new->num_linears = (num_new_linears);
	W_new->L = (vec_d *)bmalloc((num_new_linears)*sizeof(vec_d));
	W_new->L_mp = (vec_mp *)bmalloc((num_new_linears)*sizeof(vec_mp));
	int ii;
	for (ii=0; ii<num_new_linears; ii++) {
		//copy the linears into the new witness_set
		init_vec_d(W_new->L[ii],0); init_vec_mp2(W_new->L_mp[ii],0,1024); //the 1024 here is incorrect
		vec_mp_to_d(   W_new->L[ii],new_linears_full_prec[ii]);
		vec_cp_mp(W_new->L_mp[ii],new_linears_full_prec[ii]);
	}
	
	

	if (MPType==1){
		lin_to_lin_solver_mp(MPType,W,n_minusone_randomizer_matrix_full_prec,new_linears_full_prec,num_new_linears,W_new,solve_options);
	}
	else{
		lin_to_lin_solver_d( MPType,W,n_minusone_randomizer_matrix_full_prec,new_linears_full_prec,num_new_linears,W_new,solve_options);
	}
		
	
	return 0;
}


int lin_to_lin_solver_d(int MPType, //, double parse_time, unsigned int currentSeed
												witness_set W,  // includes the initial linear.
												mat_mp n_minusone_randomizer_matrix_full_prec,  // for randomizing down to N-1 equations.
												vec_mp *new_linears_full_prec,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												int num_new_linears,
												witness_set *W_new,
												solver_configuration *solve_options)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	double parse_time = 0;
  FILE *OUT = NULL, *FAIL = safe_fopen_write("failed_paths"), *midOUT = NULL, *rawOUT = safe_fopen_write("raw_data");
  tracker_config_t T;
  prog_t dummyProg;
  bclock_t time1, time2;
  int num_variables = 0, num_sols = 0;
	
	
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	
  lintolin_eval_data_d ED;
  trackingStats trackCount;
  double track_time;
	
  bclock(&time1); // initialize the clock.
  init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	//necessary for later whatnot
	int userHom = 0, useRegen = 0, pathMod = 0, paramHom = 0;
	
	cp_tracker_config_t(&T, &solve_options->T);
	
	
	
	
	//  // call the setup function
	num_variables = lin_to_lin_setup_d(&OUT, "output",
																		 &midOUT, "midpath_data",
																		 &T, &ED,
																		 &dummyProg,  //arg 7
																		 &startSub, &endSub, &startFunc, &endFunc,
																		 &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow,
																		 &ptr_to_eval_d, &ptr_to_eval_mp,
																		 "preproc_data", "deg.out",
																		 !useRegen, "nonhom_start", "start",
																		 n_minusone_randomizer_matrix_full_prec,W,
																		 solve_options);
  
	int (*change_prec)(void const *, int) = NULL;
	change_prec = &change_lintolin_eval_prec;
	
	int (*dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;
	dehom = &lintolin_dehom;
	
	
	
  // error checking
  if (userHom <= 0 && paramHom != 2)
  { // no pathvariables or parameters allowed!
    if (dummyProg.numPathVars > 0)
    { // path variable present
      printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
    if (dummyProg.numPars > 0)
    { // parameter present
      printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
  }
	
  if (T.MPType == 2)  //If we are doing adaptive precision path-tracking, we must set up AMP_eps, AMP_Phi, AMP_Psi based on config settings.
  {
    T.AMP_eps = (double) num_variables * num_variables;  //According to Demmel (as in the AMP paper), n^2 is a very reasonable bound for \epsilon.
    T.AMP_Phi = T.AMP_bound_on_degree*(T.AMP_bound_on_degree-1.0)*T.AMP_bound_on_abs_vals_of_coeffs;  //Phi from the AMP paper.
    T.AMP_Psi = T.AMP_bound_on_degree*T.AMP_bound_on_abs_vals_of_coeffs;  //Psi from the AMP paper.
    // initialize latest_newton_residual_mp to the maximum precision
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);
  }
	
	
	post_process_t *endPoints = (post_process_t *)bmalloc(W.num_pts*(num_new_linears) * sizeof(post_process_t)); //overallocate, expecting full number of solutions.
	
	
	if (T.endgameNumber == 3)
	{ // use the track-back endgame
		//        zero_dim_trackBack_d(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ED.BED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_basic_eval_prec, zero_dim_dehom);
		printf("bertini_real not equipped to deal with endgameNumber 3\nexiting\n");
		exit(-99);
	}
	else
	{ // use regular endgame
		lin_to_lin_track_d(&trackCount, OUT, rawOUT, midOUT,
											 W,  // was the startpts file pointer.
											 new_linears_full_prec,
											 num_new_linears,
											 endPoints,
											 FAIL, pathMod,
											 &T, &ED, ED.BED_mp,
											 ptr_to_eval_d, ptr_to_eval_mp,
											 change_prec, dehom,
											 solve_options);
	}
	
	fclose(midOUT);
	

	// finish the output to rawOUT
	fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
	

	// check for path crossings
	int num_crossings = 0;
	if (solve_options->use_midpoint_checker==1) {
		midpoint_checker(trackCount.numPoints, num_variables,solve_options->midpoint_tol, &num_crossings);
	}
	
	// setup num_sols
	num_sols = trackCount.successes;
	
  // we report how we did with all paths:
  bclock(&time2);
  totalTime(&track_time, time1, time2);
  if (T.screenOut)
  {
    printf("Number of failures:  %d\n", trackCount.failures);
    printf("Number of successes:  %d\n", trackCount.successes);
    printf("Number of paths:  %d\n", trackCount.numPoints);
    printf("Parse Time = %fs\n", parse_time);
    printf("Track Time = %fs\n", track_time);
  }
  fprintf(OUT, "Number of failures:  %d\n", trackCount.failures);
  fprintf(OUT, "Number of successes:  %d\n", trackCount.successes);
  fprintf(OUT, "Number of paths:  %d\n", trackCount.numPoints);
  fprintf(OUT, "Parse Time = %fs\n", parse_time);
  fprintf(OUT, "Track Time = %fs\n", track_time);
	
  // print the system to rawOUT
	//  printlintolinRelevantData(&ED, ED.BED_mp, T.MPType, usedEq, rawOUT);
	//	printZeroDimRelevantData(&ED, ED.BED_mp, T.MPType, usedEq, rawOUT);// legacy
	
  // close all of the files
  fclose(OUT);
  fclose(rawOUT);
  fprintf(FAIL, "\n");
  fclose(FAIL);
	

	
	BRpostProcessing(endPoints, W_new, trackCount.successes, &ED.preProcData, &T,solve_options);
	
	

	
//TODO: DAB is there other stuff which should be cleared here?
	
    free(startSub);
    free(endSub);
    free(startFunc);
    free(endFunc);
    free(startJvsub);
    free(endJvsub);
    free(startJv);
    free(endJv);

	
	
  lintolin_eval_clear_d(&ED, userHom, T.MPType);
  tracker_config_clear(&T);

  return 0;
}




void lin_to_lin_track_d(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set W,
												vec_mp *new_linears_full_prec,
												int num_new_linears,
												post_process_t *endPoints,
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												lintolin_eval_data_d *ED_d,
												lintolin_eval_data_mp *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
												solver_configuration *solve_options)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does standard zero dimensional tracking                *
 *  in either double precision or adaptive precision             *
 \***************************************************************/
{
	
  int ii,jj,kk, oid, startPointIndex, max = max_threads();
  tracker_config_t *T_copy = NULL;
  lintolin_eval_data_d *BED_copy = NULL;
  trackingStats *trackCount_copy = NULL;
  FILE **OUT_copy = NULL, **MIDOUT_copy = NULL, **RAWOUT_copy = NULL, **FAIL_copy = NULL, *NONSOLN = NULL, **NONSOLN_copy = NULL;
	
  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);
	
	
	int (*curr_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*curr_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	
	curr_eval_d = &lin_to_lin_eval_d;   //DAB
  curr_eval_mp = &lin_to_lin_eval_mp; // DAB 
	

	
	point_data_d *startPts = NULL;
	startPts = (point_data_d *)bmalloc(W.num_pts * sizeof(point_data_d));
	
	for (ii = 0; ii < W.num_pts; ii++)
	{ // setup startPts[ii]
		init_point_data_d(&startPts[ii], W.num_variables); // also performs initialization on the point inside startPts
		change_size_vec_d(startPts[ii].point,W.num_variables);
		startPts[ii].point->size = W.num_variables;
		
		//NEED TO COPY IN THE WITNESS POINT
	
		//1 set the coordinates
		for (jj = 0; jj<W.num_variables; jj++) {
			startPts[ii].point->coord[jj].r = W.pts_d[ii]->coord[jj].r;
			startPts[ii].point->coord[jj].i = W.pts_d[ii]->coord[jj].i;
		}
		//2 set the start time to 1.
		set_one_d(startPts[ii].time);
	}
	T->endgameOnly = 0;
	
	
  // setup the rest of the structures
	endgame_data_t *EG = NULL; //this will hold the temp solution data produced for each individual track
  setup_lin_to_lin_omp_d(max,
												 &EG, &trackCount_copy, trackCount,
												 &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN,
												 &T_copy, T,
												 &BED_copy, ED_d, ED_mp);

	
	trackCount->numPoints = W.num_pts*(num_new_linears);
	int solution_counter = 0;
	
	
	for (kk = 0; kk< num_new_linears; kk++)
	{
		if (solve_options->verbose_level>=1)
			printf("solving for linear %d\n",kk);
		
		// we pass the particulars of the information for this solve mode via the ED.
		
		//set current linear in the evaluator data's
		vec_mp_to_d(     ED_d->current_linear,new_linears_full_prec[kk]);
		
		
		if (T->MPType==2) {
			//q: should i reset the precision here?
			vec_cp_mp(ED_d->BED_mp->current_linear,new_linears_full_prec[kk]);
			vec_cp_mp(ED_d->BED_mp->current_linear_full_prec,new_linears_full_prec[kk]);
		}
		


		
		// track each of the start points
#ifdef _OPENMP
#pragma omp parallel for private(ii, oid, startPointIndex) schedule(runtime)
#endif
		for (ii = 0; ii < W.num_pts; ii++)
		{ // get current thread number
			oid = thread_num();
			
			if (solve_options->verbose_level>=1)
				printf("\t\tpoint %d\n",ii);
			
			startPointIndex = ii;

//TODO: eliminate this header stuff
			// print the header of the path to OUT
			printPathHeader_d(OUT_copy[oid], &startPts[startPointIndex], &T_copy[oid], ii, &BED_copy[oid], eval_func_d);
			
			
#ifdef printpathlintolin
			BED_copy[oid].num_steps = 0;
#endif
			
			// track the path
			lin_to_lin_track_path_d(solution_counter, &EG[oid], &startPts[startPointIndex], OUT_copy[oid], MIDOUT_copy[oid], &T_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, curr_eval_d, curr_eval_mp, change_prec, find_dehom);
			
#ifdef printpathlintolin
			fprintf(BED_copy[oid].FOUT,"-100 %d ",BED_copy[oid].num_steps);
			for (mm=0; mm<BED_copy[oid].num_variables-1; ++mm) {
				fprintf(BED_copy[oid].FOUT,"0 0 ");
			}
			fprintf(BED_copy[oid].FOUT,"\n%d\n\n",EG->retVal);
#endif
			
			
			// check to see if it should be sharpened
			if (EG[oid].retVal == 0 && T_copy[oid].sharpenDigits > 0)
			{ // use the sharpener for after an endgame
				sharpen_endpoint_endgame(&EG[oid], &T_copy[oid], OUT_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, curr_eval_d, curr_eval_mp, change_prec);
			}
			

			//get the terminal time in double form
			comp_d time_to_compare;
			if (EG->prec < 64) {
				set_d(time_to_compare,EG->PD_d.time);}
			else {
				mp_to_d(time_to_compare, EG->PD_mp.time); }
			
			
			int issoln = 0;
			if (EG->last_approx_prec!=1) {
				if (EG->prec<64){
					issoln = check_issoln_lintolin_d(&EG[oid],  &T_copy[oid], &BED_copy[oid]); }
				else {
					issoln = check_issoln_lintolin_mp(&EG[oid], &T_copy[oid], BED_copy[oid].BED_mp); }
			}
			else{
				
				if (EG->prec<64){
					print_point_to_screen_matlab(EG->PD_d.point,"solution");}
				else {
					print_point_to_screen_matlab_mp(EG->PD_mp.point,"solution");}
				
				printf("the last approximation was of precision %d\n",EG->last_approx_prec);
				printf("this is probably a problem\n");
				mypause();
			}




			
			
			if ( (EG->retVal != 0 && time_to_compare->r > T->minTrackT) || !issoln) {  // <-- this is the real indicator of failure...
				trackCount->failures++;
				printf("\nretVal = %d\nthere was a fatal path failure tracking linear %d, witness point %d\n\n",EG->retVal,kk,ii);
				print_path_retVal_message(EG->retVal);
				if (!issoln) {
					printf("the following was labeled as not a solution\n");
					if (EG->prec < 64) {
						print_point_to_screen_matlab(EG->PD_d.point,"variablevalues");
					}
					else{
						print_point_to_screen_matlab_mp(EG->PD_mp.point,"variablevalues");
					}
				}
				print_point_to_screen_matlab(BED_copy->old_linear,"old");
				print_point_to_screen_matlab(BED_copy->current_linear,"new");
				exit(EG->retVal); //failure intolerable in this solver.
			}
			else
			{
				//otherwise converged, but may have still had non-zero retval due to other reasons.
				endgamedata_to_endpoint(&endPoints[solution_counter], EG);
				trackCount->successes++;
				solution_counter++; // probably this could be eliminated
			}
		}// re: for (ii=0; ii<W.num_pts ;ii++)
	} // for each new linear
	
	
	//clear the data structures.
	
  for (ii = 0; ii >W.num_pts; ii++)
  { // clear startPts[ii]
    clear_point_data_d(&startPts[ii]);
  }
  free(startPts);
	
  // clear the structures
  clear_lintolin_omp_d(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, &BED_copy);
	
	
//TODO: CLEAR MORE MEMORY
	
	
	
  return;
}




// derived from zero_dim_track_path_d
void lin_to_lin_track_path_d(int pathNum, endgame_data_t *EG_out,
														 point_data_d *Pin,
														 FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
														 void const *ED_d, void const *ED_mp,
														 int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
														 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 int (*change_prec)(void const *, int),
														 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: actually does the zero-dimensional tracking and sets   *
 *  up EG_out                                                    *
 \***************************************************************/
{

	
	
	EG_out->pathNum = pathNum;
  EG_out->codim = 0; // zero dimensional - this is ignored
	
  T->first_step_of_path = 1;
	
//	printf("mytype %d\n",T->MPType);
  if (T->MPType == 2)
  { // track using AMP
		
		//verification code:  
//		lintolin_eval_data_mp *LtLED = (lintolin_eval_data_mp *)ED_mp; // to avoid having to cast every time
//		print_matrix_to_screen_matlab_mp(LtLED->n_minusone_randomizer_matrix,"nminus1");
//		print_point_to_screen_matlab_mp(LtLED->old_linear,"old");
//		print_point_to_screen_matlab_mp(LtLED->current_linear,"curr");
//		print_matrix_to_screen_matlab_mp(LtLED->patch.patchCoeff,"patch");
//		printf("patch precision %d\n",LtLED->patch.curr_prec);
		
    EG_out->prec = EG_out->last_approx_prec = 1;
		
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
//			printf("%d prec in track_path\n",EG_out->prec);
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
  else
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




// derived from zero_dim_basic_setup_d
int lin_to_lin_setup_d(FILE **OUT, char *outName,
											 FILE **midOUT, char *midName,
											 tracker_config_t *T,
											 lintolin_eval_data_d *ED,
											 prog_t *dummyProg,
											 int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
											 int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
											 int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
											 char *preprocFile, char *degreeFile,
											 int findStartPts, char *pointsIN, char *pointsOUT,
											 mat_mp n_minusone_randomizer_matrix_full_prec,
											 witness_set W,
											 solver_configuration *solve_options)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES: number of original variables                   *
 * NOTES: setup for zero dimensional tracking                    *
 \***************************************************************/
{ // need to create the homotopy
//	printf("entering lin_to_lin_setup_d 606\n");
  int rank, patchType, ssType, numOrigVars, adjustDegrees, numGps;
	
  *eval_d = &lin_to_lin_eval_d;   //DAB
  *eval_mp = &lin_to_lin_eval_mp; // DAB  // lol send them to the same place for now.
	
	//  *eval_mp = &lintolin_eval_mp; // DAB
	
  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");
	
  if (T->MPType == 2) // using AMP - need to allocate space to store BED_mp
    ED->BED_mp = (lintolin_eval_data_mp *)bmalloc(1 * sizeof(lintolin_eval_data_mp));
  else
    ED->BED_mp = NULL;
	
	
  // setup a straight-line program, using the file(s) created by the parser
  T->numVars = numOrigVars = setupProg_count(dummyProg, T->Precision, T->MPType, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);
	
  // setup preProcData
  setupPreProcData(preprocFile, &ED->preProcData);
	
	
	
  numGps = ED->preProcData.num_var_gp + ED->preProcData.num_hom_var_gp;//this should probably be removed
  // find the rank
  rank = rank_finder_d(&ED->preProcData, dummyProg, T, T->numVars); //this should probably be removed
	

    patchType = 2; // 1-hom patch
    ssType = 0;    // with 1-hom, we use total degree start system  // in BR, this is irrelevant
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay  // irrelevant in BR
    setuplintolinEval_d(T,preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->MPType, &T->numVars, NULL, NULL, NULL, ED, adjustDegrees, n_minusone_randomizer_matrix_full_prec, W,solve_options);

	
	
#ifdef printpathlintolin
	int ii;
	ED->FOUT = safe_fopen_append("pathtrack_lintolin");
	fprintf(ED->FOUT,"%d ",W.num_variables);
	for (ii=0; ii<W.num_variables; ii++) {
		fprintf(ED->FOUT,"%s ",W.variable_names[ii]);
	}
	fprintf(ED->FOUT,"\n%d ",W.num_pts);
	fprintf(ED->FOUT,"%d %d %d ",T->MPType, T->odePredictor, T->endgameNumber);
	fprintf(ED->FOUT,"\n");
#endif
	
	
  return numOrigVars;
}





//this derived from basic_eval_d
int lin_to_lin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, build for bertini_real
	
	// uncomment to see the time at each step.
//	printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);

  lintolin_eval_data_d *BED = (lintolin_eval_data_d *)ED; // to avoid having to cast every time
	
  int ii, jj; // counters
  comp_d one_minus_s, gamma_s;

	
  set_one_d(one_minus_s);
  sub_d(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
  mul_d(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	comp_d temp, temp2;
	
	vec_d patchValues; init_vec_d(patchValues, 0);
	vec_d temp_function_values; init_vec_d(temp_function_values,0);
	vec_d AtimesF; init_vec_d(AtimesF,BED->n_minusone_randomizer_matrix->rows); AtimesF->size = BED->n_minusone_randomizer_matrix->rows;// declare  // initialize
	vec_d vars_times_curr_linear; init_vec_d(vars_times_curr_linear,BED->num_variables);  vars_times_curr_linear->size = BED->num_variables;
	vec_d vars_times_old_linear; init_vec_d(vars_times_old_linear,BED->num_variables);  vars_times_old_linear->size =  BED->num_variables;
	vec_d one_minus_s_times_current_linear; init_vec_d(one_minus_s_times_current_linear,BED->num_variables); one_minus_s_times_current_linear->size = BED->num_variables;
	vec_d gamma_s_times_old_linear; init_vec_d(gamma_s_times_old_linear,BED->num_variables); gamma_s_times_old_linear->size = BED->num_variables;
	
	mat_d temp_jacobian_functions; init_mat_d(temp_jacobian_functions,BED->n_minusone_randomizer_matrix->cols,BED->num_variables);
		temp_jacobian_functions->rows = BED->n_minusone_randomizer_matrix->cols; temp_jacobian_functions->cols = BED->num_variables;
	mat_d temp_jacobian_parameters; init_mat_d(temp_jacobian_parameters,0,0);
	mat_d Jv_Patch; init_mat_d(Jv_Patch, 0, 0);
	mat_d AtimesJ; init_mat_d(AtimesJ,BED->n_minusone_randomizer_matrix->rows,BED->num_variables);
		AtimesJ->rows = BED->n_minusone_randomizer_matrix->rows; AtimesJ->cols = BED->num_variables;
	
	
	//set the sizes
	change_size_vec_d(funcVals,BED->num_variables); funcVals->size = BED->num_variables;
  change_size_mat_d(Jv, BED->num_variables, BED->num_variables); Jv->rows = Jv->cols = BED->num_variables; //  -> this should be square!!!

	
	

	
	// evaluate the SLP to get the system's whatnot.
	evalProg_d(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	
	
  // evaluate the patch
  patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	
	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	// note that you can only really do this after you are done calling other evaluators.
  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
	change_size_mat_d(Jp, BED->num_variables, 1); Jp->rows = BED->num_variables; Jp->cols = 1;

  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	///////// / / / /  /   /
  // combine everything
	///////// / / / /  /   /
	
	//perform the randomization multiplications
	mat_mul_d(AtimesJ,BED->n_minusone_randomizer_matrix,temp_jacobian_functions);
	mul_mat_vec_d(AtimesF,BED->n_minusone_randomizer_matrix, temp_function_values ); // set values of AtimesF (A is randomization matrix)
	
	for (ii=0; ii<AtimesF->size; ii++) { // for each function, after (real orthogonal) randomization
		set_d(&funcVals->coord[ii], &AtimesF->coord[ii]);
	}
	


	//
	// function values for the linears
	//
	
	// multiply vars times the new linear, with (1-s)
	for (ii=0; ii<BED->num_variables; ii++) { // for each variable, including the homogeneous ones.
		mul_d(&one_minus_s_times_current_linear->coord[ii], &BED->current_linear->coord[ii], one_minus_s);
		mul_d(&vars_times_curr_linear->coord[ii], &one_minus_s_times_current_linear->coord[ii], &current_variable_values->coord[ii]);
	}
	// multiply vars times the old linear, with gamma*s
	for (ii=0; ii<BED->num_variables; ii++) { // for each variable, including the homogeneous ones.
		mul_d(&gamma_s_times_old_linear->coord[ii], &BED->old_linear->coord[ii], gamma_s);
		mul_d(&vars_times_old_linear->coord[ii],&gamma_s_times_old_linear->coord[ii],&current_variable_values->coord[ii]);
	}
	// add the old and the new
	set_zero_d(temp);  //  initialize to 0
	for (ii=0; ii<BED->num_variables; ii++) {  // for each variable, including the homogenizing ones.
		add_d(temp2, &vars_times_old_linear->coord[ii], &vars_times_curr_linear->coord[ii]);
		add_d(temp,temp,temp2); //  tested correct, does not overwrite values in temp before accessed for read.
	}
	//finally set the value of the linears we are homotoping.
	set_d(&funcVals->coord[BED->num_variables-2],temp); // this is the entry for the linear's homotopy.
	
	
	
	//set the function PATCH values
	int offset = BED->num_variables-1;
	for (ii=0; ii<BED->patch.num_patches; ii++) {
		set_d(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
	}
	
	
	

	//////////////
	//
	// set the JACOBIAN values.
	//
	////////////////////
	
	//first, the entries related to the functions
	
  for (ii = 0; ii < BED->n_minusone_randomizer_matrix->rows; ii++)
  {
		for (jj = 0; jj < BED->num_variables; jj++)
		{
			set_d(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
		}
  }
	

	

	
	// second, the entries for the linears we are homotoping
	
	
	offset = BED->num_variables-2;
	for (<#initialization#>; <#condition#>; <#increment#>) {
		for (ii=0; ii<BED->num_variables; ii++) {
			add_d(temp,&gamma_s_times_old_linear[]->coord[ii], &one_minus_s_times_current_linear[]->coord[ii]);
			set_d(&Jv->entry[offset][ii], temp); // HERE DAB!!!
		}
	}

	
	offset = BED->num_variables-1;
	for (ii=0; ii<BED->num_variables; ii++) {
		set_d(&Jv->entry[offset][ii],&Jv_Patch->entry[0][ii]);
	}
	

	
	
	for (ii = 0; ii<BED->num_variables-2; ii++) {
		set_zero_d(&Jp->entry[ii][0]);  // no parameter dependence means zero derivative for these functions
	}
	

	// Jp = -current_linear_times_vars + gamma*old_linear_times_vars
	set_zero_d(&Jp->entry[BED->num_variables-2][0]);
//TODO:  this goes in to a loop over the number of linears
	dot_product_d(temp,BED->current_linear,current_variable_values);
	neg_d(temp,temp);
	dot_product_d(temp2,BED->old_linear,current_variable_values);
	mul_d(temp2,temp2,BED->gamma);
	add_d(&Jp->entry[BED->num_variables-2][0],temp,temp2);
	
	
	
	
	// the entries in the jacobian for the patch equations.
	offset = BED->num_variables-1;
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		// Jp = 0
		set_zero_d(&Jp->entry[ii+offset][0]);
		// Jv = Jv_Patch
		for (jj = 0; jj<BED->num_variables; jj++) // for each variable
		{
			set_d(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
		}
	}
	
	// done!  yay!
	
	
	//uncomment to see screen output of important variables at each solve step.
//	printf("gamma = %lf+1i*%lf;\n", BED->gamma->r, BED->gamma->i);
//	printf("time = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
//	print_matrix_to_screen_matlab( temp_jacobian_functions,"jac");
//	print_point_to_screen_matlab(current_variable_values,"currvars");
//	print_point_to_screen_matlab(BED->current_linear,"new");
//	print_point_to_screen_matlab(BED->old_linear,"old");
//	print_point_to_screen_matlab(funcVals,"F");
//	print_matrix_to_screen_matlab(Jv,"Jv");
//	print_matrix_to_screen_matlab(Jp,"Jp");
//	print_matrix_to_screen_matlab(BED->n_minusone_randomizer_matrix,"n_minusone_randomizer_matrix");
//
//	mypause();
//

//	fprintf(BED->FOUT,"jamesbrown\n");
	
	clear_vec_d(patchValues);
	clear_vec_d(temp_function_values);
	clear_vec_d(AtimesF); 
	clear_vec_d(vars_times_curr_linear);
	clear_vec_d(vars_times_old_linear);
	clear_vec_d(one_minus_s_times_current_linear);
	clear_vec_d(gamma_s_times_old_linear);
	
	clear_mat_d(temp_jacobian_functions);
	clear_mat_d(temp_jacobian_parameters);
	clear_mat_d(Jv_Patch);
	clear_mat_d(AtimesJ);
	
#ifdef printpathlintolin
	BED->num_steps++;
	vec_d dehommed; init_vec_d(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize(&dehommed,current_variable_values);
	fprintf(BED->FOUT,"%.15lf %.15lf ", pathVars->r, pathVars->i);
	for (ii=0; ii<BED->num_variables-1; ++ii) {
		fprintf(BED->FOUT,"%.15lf %.15lf ",dehommed->coord[ii].r,dehommed->coord[ii].i);
	}
	fprintf(BED->FOUT,"\n");
	clear_vec_d(dehommed);
#endif
	
	
  return 0;
}




void printlintolinRelevantData(lintolin_eval_data_d *ED_d, lintolin_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: prints the relevant data to FP so that we can recover  *
 * the system if needed for sharpening                           *
 \***************************************************************/
{
  // print the MPType and if an eq-by-eq method (diagona/regen) was used
//  fprintf(FP, "%d %d\n", MPType, eqbyeqMethod);
	//
	//  // print the patch
//  printPatchCoeff(FP, MPType, ED_d, ED_mp);
	//
	//  // print the start system
	//  printStartSystem(FP, MPType, ED_d, ED_mp);
	//
	//  // print the square system
	//  printSquareSystem(FP, MPType, ED_d, ED_mp);
	
  return;
}



void lintolin_eval_clear_d(lintolin_eval_data_d *ED, int clearRegen, int MPType)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: clear ED                                               *
 \***************************************************************/
{
	//clear the patch
  patch_eval_data_clear_d(&ED->patch);
  preproc_data_clear(&ED->preProcData);
	
	
	
  if (MPType == 2)
  { // clear the MP stuff
    lintolin_eval_clear_mp(ED->BED_mp, 0, 0);
  }
	
	//specifics for the lintolin method.
	clear_vec_d(ED->current_linear);
	clear_vec_d(ED->old_linear);
	clear_mat_d(ED->n_minusone_randomizer_matrix);
	
	
#ifdef printpathlintolin
	fclose(ED->FOUT);
#endif
  return;
}






void setup_lin_to_lin_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														FILE ***NONSOLN_copy, FILE *NONSOLN,
														tracker_config_t **T_copy, tracker_config_t *T,
														lintolin_eval_data_d **BED_copy, lintolin_eval_data_d *ED_d, lintolin_eval_data_mp *ED_mp)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: setup everything needed to do zero dimensional tracking*
 *  using OpenMP                                                 *
 \***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int ii;
  // error checking
  if (max_threads <= 0)
  {
    printf("\n\nERROR: The number of threads (%d) needs to be positive when setting up for tracking!\n", max_threads);
    bexit(ERROR_CONFIGURATION);
  }
	
	
	
	
  // allocate space for EG
  *EG = (endgame_data_t *)bmalloc(max_threads * sizeof(endgame_data_t));
	
	// initialize
  for (ii = 0; ii < max_threads; ii++)
	{
    if (T->MPType == 2)
    { // initialize for AMP tracking
      init_endgame_data(&(*EG)[ii], 64);
    }
    else
    { // initialize for double precision tracking
      init_endgame_data(&(*EG)[ii], 52);
    }
	}
	
  // allocate space to hold pointers to the files
  *OUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *MIDOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *RAWOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *FAIL_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *NONSOLN_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
	
  if (max_threads == 1)
  { // setup the pointers
    *trackCount_copy = trackCount;
    *T_copy = T;
    *BED_copy = ED_d;
    (*BED_copy)->BED_mp = ED_mp; // make sure that this is pointed to inside of ED_d
		
    (*OUT_copy)[0] = OUT;
    (*RAWOUT_copy)[0] = RAWOUT;
    (*MIDOUT_copy)[0] = MIDOUT;
    (*FAIL_copy)[0] = FAIL;
    (*NONSOLN_copy)[0] = NONSOLN;
  }
  else // max_threads > 1
  { // allocate memory
    *trackCount_copy = (trackingStats *)bmalloc(max_threads * sizeof(trackingStats));
    *T_copy = (tracker_config_t *)bmalloc(max_threads * sizeof(tracker_config_t));
    *BED_copy = (lintolin_eval_data_d *)bmalloc(max_threads * sizeof(lintolin_eval_data_d));
		
    // copy T, ED_d, ED_mp, & trackCount
    for (ii = 0; ii < max_threads; ii++)
    { // copy T
      cp_tracker_config_t(&(*T_copy)[ii], T);
      // copy ED_d & ED_mp
      cp_lintolin_eval_data_d(&(*BED_copy)[ii], ED_d, ED_mp, T->MPType);
      // initialize trackCount_copy
      init_trackingStats(&(*trackCount_copy)[ii]);
      (*trackCount_copy)[ii].numPoints = trackCount->numPoints;
    }
		
    // setup the files
    char *str = NULL;
    int size;
    for (ii = 0; ii < max_threads; ii++)
    {
      size = 1 + snprintf(NULL, 0, "output_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", ii);
      (*OUT_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "midout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", ii);
      (*MIDOUT_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "rawout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", ii);
      (*RAWOUT_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "fail_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", ii);
      (*FAIL_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", ii);
      (*NONSOLN_copy)[ii] = fopen(str, "w+");
    }
    free(str);
  }
	
	
  return;
}

void clear_lintolin_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, lintolin_eval_data_d **BED_copy)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: copy the relevant data back to the standard spot and   *
 *  clear the allocated data that was used by OpenMP             *
 \***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int ii;
	
  // clear EG
  for (ii = max_threads - 1; ii >= 0; ii--)
  {
    clear_endgame_data(&(*EG)[ii]);
  }
  free(*EG);
	
  if (max_threads == 1)
  { // set the pointers to NULL since they just pointed to the actual values
    *trackCount_copy = NULL;
    *T_copy = NULL;
    *BED_copy = NULL;
		
    *OUT_copy[0] = NULL;
    *RAWOUT_copy[0] = NULL;
    *MIDOUT_copy[0] = NULL;
    *FAIL_copy[0] = NULL;
		
    // free the memory of the file pointers
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }
  else if (max_threads > 1)
  {
    // combine trackCount_copy
    add_trackingStats(trackCount, *trackCount_copy, max_threads);
		
    // clear the copies T, ED_d & ED_mp
    for (ii = max_threads - 1; ii >= 0; ii--)
    { // clear BED_copy - 0 since not using regeneration
      lintolin_eval_clear_d(&(*BED_copy)[ii], 0, (*T_copy)[ii].MPType);
			// clear T_copy
      tracker_config_clear(&(*T_copy)[ii]);
    }
		
    // free the memory
    free(*trackCount_copy);
    free(*T_copy);
    free(*BED_copy);
		
    // copy all of the files to the appropriate place
    char ch, *str = NULL;
    int size;
    for (ii = 0; ii < max_threads; ii++)
    {
      // rewind to beginning for reading
      rewind((*OUT_copy)[ii]);
      // copy over to OUT
      ch = fgetc((*OUT_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(OUT, "%c", ch);
        ch = fgetc((*OUT_copy)[ii]);
      }
      // close file & delete
      fclose((*OUT_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "output_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*MIDOUT_copy)[ii]);
      // copy over to MIDOUT
      ch = fgetc((*MIDOUT_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(MIDOUT, "%c", ch);
        ch = fgetc((*MIDOUT_copy)[ii]);
      }
      // close file & delete
      fclose((*MIDOUT_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "midout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*RAWOUT_copy)[ii]);
      // copy over to RAWOUT
      ch = fgetc((*RAWOUT_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(RAWOUT, "%c", ch);
        ch = fgetc((*RAWOUT_copy)[ii]);
      }
      // close file & delete
      fclose((*RAWOUT_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "rawout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*FAIL_copy)[ii]);
      // copy over to FAIL
      ch = fgetc((*FAIL_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(FAIL, "%c", ch);
        ch = fgetc((*FAIL_copy)[ii]);
      }
      // close file & delete
      fclose((*FAIL_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "fail_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*NONSOLN_copy)[ii]);
      // copy over to NONSOLN
      ch = fgetc((*NONSOLN_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(NONSOLN, "%c", ch);
        ch = fgetc((*NONSOLN_copy)[ii]);
      }
      // close file & delete
      fclose((*NONSOLN_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", ii);
      remove(str);
    }
    free(str);
    // free file memory
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }
	
  return;
}




void setuplintolinEval_d(tracker_config_t *T,char preprocFile[], char degreeFile[], prog_t *dummyProg,
												 int squareSize, int patchType, int ssType, int MPType,
												 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,// what are these supposed to point to?
												 lintolin_eval_data_d *BED, int adjustDegrees,
												 mat_mp n_minusone_randomizer_matrix_full_prec,
												 witness_set W,
												 solver_configuration *solve_options)
{
  int ii;
	BED->num_variables = W.num_variables;
	
	
  setupPreProcData(preprocFile, &BED->preProcData);

  setupPatch_d(patchType, &BED->patch, ptr1, ptr2);
	for (ii = 0; ii < BED->num_variables ; ii++)
		set_d(&BED->patch.patchCoeff->entry[0][ii],&W.patch[0]->coord[ii]);
	


	BED->SLP = dummyProg;

	
	init_mat_d(BED->n_minusone_randomizer_matrix,0,0);
	mat_mp_to_d(BED->n_minusone_randomizer_matrix,
					 n_minusone_randomizer_matrix_full_prec);
	
	
	// set up the vectors to hold the two linears.
	
	
	init_vec_d(BED->current_linear,0);
	change_size_vec_d(BED->current_linear, W.num_variables);
	BED->current_linear->size =  W.num_variables;
	
	
	
	
	
	init_vec_d(BED->old_linear,0);
	change_size_vec_d(BED->old_linear, dummyProg->numVars);
	BED->old_linear->size =  dummyProg->numVars; // seriously, why is this not in the bertini method?

	vec_cp_d(BED->old_linear,W.L[0]);

	
	if (solve_options->use_gamma_trick==1)
		get_comp_rand_d(BED->gamma); // set gamma to be random complex value
	else
		set_one_d(BED->gamma);
	
	
	if (MPType == 2)
  { // using AMP - initialize using 16 digits & 64-bit precison
    int digits = 16, prec = 64;
		initMP(prec);
		BED->BED_mp->curr_prec = prec;
		
		
#ifdef printpathlintolin
		BED->BED_mp->FOUT = BED->FOUT;
#endif
		
		
		BED->BED_mp->gamma_rat = (mpq_t *)bmalloc(2 * sizeof(mpq_t));
		if (solve_options->use_gamma_trick==1){
			get_comp_rand_rat(BED->gamma, BED->BED_mp->gamma, BED->BED_mp->gamma_rat, prec, T->AMP_max_prec, 1, 1);
		}
		else{
			init_rat(BED->BED_mp->gamma_rat);
			init_mp(BED->BED_mp->gamma);
			set_one_d(BED->gamma);
			set_one_mp(BED->BED_mp->gamma);
			set_one_rat(BED->BED_mp->gamma_rat);
		}
		
		
		BED->BED_mp->num_variables = W.num_variables;

		
		BED->BED_mp->SLP = BED->SLP; // assign the SLP pointer

		init_vec_mp2(BED->BED_mp->old_linear,W.num_variables,prec);
		init_vec_mp2(BED->BED_mp->old_linear_full_prec,W.num_variables,T->AMP_max_prec);
		vec_cp_mp(BED->BED_mp->old_linear,W.L_mp[0]); //copy in the old linear
		vec_cp_mp(BED->BED_mp->old_linear_full_prec,W.L_mp[0]); //copy in the old linear

		init_vec_mp2(BED->BED_mp->current_linear,W.num_variables,prec); //initialize the new linear container
		init_vec_mp2(BED->BED_mp->current_linear_full_prec,W.num_variables,T->AMP_max_prec); //initialize the new linear container
		
		init_mat_mp2(BED->BED_mp->n_minusone_randomizer_matrix,0,0,prec); // initialize the randomizer matrix
		init_mat_mp2(BED->BED_mp->n_minusone_randomizer_matrix_full_prec,0,0,T->AMP_max_prec); // initialize the randomizer matrix
		
		mat_cp_mp(BED->BED_mp->n_minusone_randomizer_matrix_full_prec,n_minusone_randomizer_matrix_full_prec);
		mat_cp_mp(BED->BED_mp->n_minusone_randomizer_matrix,n_minusone_randomizer_matrix_full_prec);
		
		
    // setup preProcData
    setupPreProcData(preprocFile, &BED->BED_mp->preProcData);
		
		
		
		// setup the patch
    setupPatch_d_to_mp(&BED->patch, &BED->BED_mp->patch, digits, prec, patchType, dummyProg->numVars, ptr3, ptr4);
		for (ii = 0; ii < BED->num_variables ; ii++)
		{
			set_mp(&BED->BED_mp->patch.patchCoeff->entry[0][ii],&W.patch_mp[0]->coord[ii]);
			mp_to_rat(BED->BED_mp->patch.patchCoeff_rat[0][ii],&W.patch_mp[0]->coord[ii]);
		}
		
		for (ii = 0; ii < BED->num_variables ; ii++)
			mp_to_d(&BED->patch.patchCoeff->entry[0][ii],&W.patch_mp[0]->coord[ii]);
		
		
  }//re: if mptype==2

	
	
	//	//     this commented code is for checking that the first point solves the patch equation.   it damn well should!
//	print_matrix_to_screen_matlab(BED->patch.patchCoeff,"userpatch");
//	
//	vec_d patchValues;
//	init_vec_d(patchValues,0);
//	point_d parVals;
//	init_vec_d(parVals,0);
//	vec_d parDer;
//	init_vec_d(parDer,0);
//	mat_d Jv_Patch;
//	init_mat_d(Jv_Patch,0,0);
//	mat_d Jp;
//	init_mat_d(Jp,0,0);
//	comp_d pathVars;
//	set_one_d(pathVars);
//	patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W.pts_d[0], pathVars, &BED->patch);  // Jp is ignored
//
//	print_point_to_screen_matlab(patchValues,"patchvalues");
//	print_point_to_screen_matlab(W.pts_d[0],"initialpoint");
//	mypause();
	
	
	
//	BED->FOUT = safe_fopen_write("pathhistory_lintolin");
//	printf("opened pathhistory\n");
//	mypause();
	
	
  return;
}



void cp_lintolin_eval_data_d(lintolin_eval_data_d *BED, lintolin_eval_data_d *BED_d_input, lintolin_eval_data_mp *BED_mp_input, int MPType)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: stores a copy of BED_(t)_input to BED                  *
 \***************************************************************/
{
	printf("entering cp lintolin_eval_data_d\nthis function needs much attention, as things which should be copied are not!\n");
	exit(-1);
  cp_preproc_data(&BED->preProcData, &BED_d_input->preProcData);
	//  cp_square_system_d(&BED->squareSystem, &BED_d_input->squareSystem);
  cp_patch_d(&BED->patch, &BED_d_input->patch);
	//  cp_start_system_d(&BED->startSystem, &BED_d_input->startSystem);
	BED->SLP = (prog_t *)bmalloc(1 * sizeof(prog_t));
  cp_prog_t(BED->SLP, BED_d_input->SLP);
	
	
	//HERE COPY THE MATRICES  DAB !!!
	
	set_d(BED->gamma, BED_d_input->gamma);
  if (MPType == 2)
  { // need to also setup MP versions since using AMP
    BED->BED_mp = (lintolin_eval_data_mp *)bmalloc(1 * sizeof(lintolin_eval_data_mp));
		
    cp_preproc_data(&BED->BED_mp->preProcData, &BED_mp_input->preProcData);
    // simply point to the SLP that was setup in BED
		//    cp_square_system_mp(&BED->BED_mp->squareSystem, &BED_mp_input->squareSystem, 0, BED->squareSystem.Prog);
    cp_patch_mp(&BED->BED_mp->patch, &BED_mp_input->patch);
		//    cp_start_system_mp(&BED->BED_mp->startSystem, &BED_mp_input->startSystem);
  }
  else
    BED->BED_mp = NULL;
	
  return;
}




int lintolin_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: compute the dehom point                                *
 \***************************************************************/
{
//  basic_eval_data_d *BED_d = NULL;
//  basic_eval_data_mp *BED_mp = NULL;
	
  *out_prec = in_prec;
	
	
	
  if (in_prec < 64)
  { // compute out_d
//    BED_d = (basic_eval_data_d *)ED_d;
		point_cp_d(out_d,in_d);
//    getDehomPoint_d(out_d, in_d, in_d->size, &BED_d->preProcData);
  }
  else
  { // compute out_mp
//    BED_mp = (basic_eval_data_mp *)ED_mp;
    // set prec on out_mp
    setprec_point_mp(out_mp, *out_prec);
		point_cp_mp(out_mp,in_mp);
//    getDehomPoint_mp(out_mp, in_mp, in_mp->size, &BED_mp->preProcData);
  }
	
//  BED_d = NULL;
//  BED_mp = NULL;
	
  return 0;
}









int change_lintolin_eval_prec(void const *ED, int prec)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: change precision for standard zero dimensional solving *
 \***************************************************************/
{
  change_lintolin_eval_prec_mp(prec, (lintolin_eval_data_mp *)ED);
  return 0;
}

void change_lintolin_eval_prec_mp(int new_prec, lintolin_eval_data_mp *BED)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: change precision for the main part of BED              *
 \***************************************************************/
{
	

	
  // change the precision for the patch
  changePatchPrec_mp(new_prec, &BED->patch);
	

	if (new_prec != BED->curr_prec){
		BED->curr_prec = new_prec;
		
		setprec_mp(BED->gamma, new_prec);
		mpf_set_q(BED->gamma->r, BED->gamma_rat[0]);
		mpf_set_q(BED->gamma->i, BED->gamma_rat[1]);
		
		change_prec_point_mp(BED->current_linear,new_prec);
		vec_cp_mp(BED->current_linear, BED->current_linear_full_prec);
		
		change_prec_point_mp(BED->old_linear,new_prec);
		vec_cp_mp(BED->old_linear, BED->old_linear_full_prec);
		
		change_prec_mat_mp(BED->n_minusone_randomizer_matrix,new_prec);
		mat_cp_mp(BED->n_minusone_randomizer_matrix,BED->n_minusone_randomizer_matrix_full_prec);
		
	}
	

	


	
  return;
}














int lin_to_lin_solver_mp(int MPType,
												 witness_set W,  // includes the initial linear.
												 mat_mp n_minusone_randomizer_matrix,  // for randomizing down to N-1 equations.
												 vec_mp *new_linears,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												 int num_new_linears,
												 witness_set *W_new,
												 solver_configuration *solve_options)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	double parse_time = 0;
  FILE *OUT = NULL, *FAIL = fopen("failed_paths", "w"), *midOUT = NULL, *rawOUT = fopen("raw_data", "w");
  tracker_config_t T;
  prog_t dummyProg;
  bclock_t time1, time2;
  int num_variables = 0, num_sols = 0;
	
	
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	
  lintolin_eval_data_mp ED;  // was basic_eval_data_d  DAB
  trackingStats trackCount;
  double track_time;
	
  bclock(&time1); // initialize the clock.
  init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	//necessary for later whatnot
	int userHom = 0, useRegen = 0, pathMod = 0, paramHom = 0;
	
	cp_tracker_config_t(&T, &solve_options->T);
	

//	
	// initialize latest_newton_residual_mp
  mpf_init(T.latest_newton_residual_mp);   //<------ THIS LINE IS ABSOLUTELY CRITICAL TO CALL
	
	
	//  // call the setup function
	// setup for standard tracking - 'useRegen' is used to determine whether or not to setup 'start'
	num_variables = lin_to_lin_setup_mp(&OUT, "output",
																		 &midOUT, "midpath_data",
																		 &T, &ED,
																		 &dummyProg,  //arg 7
																		 &startSub, &endSub, &startFunc, &endFunc,
																		 &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow,
																		 &ptr_to_eval_d, &ptr_to_eval_mp,  //args 17,18
																		 "preproc_data", "deg.out",
																		 !useRegen, "nonhom_start", "start",
																		 n_minusone_randomizer_matrix,W,
																			solve_options);
  
	int (*change_prec)(void const *, int) = NULL;
	change_prec = &change_basic_eval_prec;
	
	int (*dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;
	dehom = &lintolin_dehom;
	
	
	
  // error checking
  if (userHom <= 0 && paramHom != 2)
  { // no pathvariables or parameters allowed!
    if (dummyProg.numPathVars > 0)
    { // path variable present
      printf("ERROR: Bertini does not expect path variables when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
    if (dummyProg.numPars > 0)
    { // parameter present
      printf("ERROR: Bertini does not expect parameters when user-defined homotopies are not being used!\n");
      bexit(ERROR_INPUT_SYSTEM);
    }
  }
	
	
	post_process_t *endPoints = (post_process_t *)bmalloc(W.num_pts*(num_new_linears) * sizeof(post_process_t)); //overallocate, expecting full number of solutions.
	
	
	
	if (T.endgameNumber == 3)
	{ // use the track-back endgame
		//        zero_dim_trackBack_d(&trackCount, OUT, rawOUT, midOUT, StartPts, FAIL, pathMod, &T, &ED, ED.BED_mp, ptr_to_eval_d, ptr_to_eval_mp, change_basic_eval_prec, zero_dim_dehom);
		printf("bertini_real not equipped to deal with endgameNumber 3\nexiting\n");
		exit(-99);
	}
	else
	{ // use regular endgame
		lin_to_lin_track_mp(&trackCount, OUT, rawOUT, midOUT,
											 W,  // was the startpts file pointer.
											 new_linears,
											 num_new_linears,
											 endPoints,
											 FAIL, pathMod,
											 &T, &ED,
											 ptr_to_eval_mp, //ptr_to_eval_d,
												change_prec, dehom,
												solve_options);
	}
	
	
	
	
	fclose(midOUT);
	
	
	
	
	// finish the output to rawOUT
	fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
	
	// check for path crossings
	int num_crossings = 0;
	if (solve_options->use_midpoint_checker==1) {
		midpoint_checker(trackCount.numPoints, num_variables,solve_options->midpoint_tol, &num_crossings);
	}
	
	// setup num_sols
	num_sols = trackCount.successes;
	
  // we report how we did with all paths:
  bclock(&time2);
  totalTime(&track_time, time1, time2);
  if (T.screenOut)
  {
    printf("Number of failures:  %d\n", trackCount.failures);
    printf("Number of successes:  %d\n", trackCount.successes);
    printf("Number of paths:  %d\n", trackCount.numPoints);
    printf("Parse Time = %fs\n", parse_time);
    printf("Track Time = %fs\n", track_time);
  }
  fprintf(OUT, "Number of failures:  %d\n", trackCount.failures);
  fprintf(OUT, "Number of successes:  %d\n", trackCount.successes);
  fprintf(OUT, "Number of paths:  %d\n", trackCount.numPoints);
  fprintf(OUT, "Parse Time = %fs\n", parse_time);
  fprintf(OUT, "Track Time = %fs\n", track_time);
	

  // close all of the files
  fclose(OUT);
  fclose(rawOUT);
  fprintf(FAIL, "\n");
  fclose(FAIL);
	

	BRpostProcessing(endPoints, W_new, trackCount.successes, &ED.preProcData, &T,solve_options);
	
	


	free(startSub);
	free(endSub);
	free(startFunc);
	free(endFunc);
	free(startJvsub);
	free(endJvsub);
	free(startJv);
	free(endJv);

	
	
  lintolin_eval_clear_mp(&ED, userHom, T.MPType);
  tracker_config_clear(&T);
	
  return 0;
}




void lin_to_lin_track_mp(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set W,
												vec_mp *new_linears,
												int num_new_linears,
												post_process_t *endPoints,
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												lintolin_eval_data_mp *ED,
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *),
												 solver_configuration *solve_options)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does standard zero dimensional tracking                *
 *  in either double precision or adaptive precision             *
 \***************************************************************/
{
	
  int ii,kk, oid, max = max_threads();
  tracker_config_t *T_copy = NULL;
  lintolin_eval_data_mp *BED_copy = NULL;
  trackingStats *trackCount_copy = NULL;
  FILE **OUT_copy = NULL, **MIDOUT_copy = NULL, **RAWOUT_copy = NULL, **FAIL_copy = NULL, *NONSOLN = NULL, **NONSOLN_copy = NULL;
	
  // top of RAWOUT - number of variables and that we are doing zero dimensional
  fprintf(RAWOUT, "%d\n%d\n", T->numVars, 0);
	
	
	
	// setup NONSOLN
	//  if (!ED_d->squareSystem.noChanges)
	//  { // setup NONSOLN
	//    NONSOLN = fopen("nonsolutions", "w");
	//    fprintf(NONSOLN, "                                    \n\n");
	//  }
	

  int (*curr_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
  curr_eval_mp = &lin_to_lin_eval_mp; // custom evaluator for this method
	
	
	point_data_mp *startPts = NULL;
	startPts = (point_data_mp *)bmalloc(W.num_pts * sizeof(point_data_mp));
	
	
	for (ii = 0; ii < W.num_pts; ii++)
	{ // setup startPts[ii]
		init_point_data_mp2(&startPts[ii], W.num_variables, T->Precision);
		startPts[ii].point->size = W.num_variables;
		
		//NEED TO COPY IN THE WITNESS POINT
		
		//1 set the coordinates
		vec_cp_mp(startPts[ii].point, W.pts_mp[ii] );

		//2 set the start time to 1.
		set_one_mp(startPts[ii].time);
	}
	
	
	T->endgameOnly = 0;
	
	
  // setup the rest of the structures
	endgame_data_t *EG = NULL;
  setup_lin_to_lin_omp_mp(max,
												 &EG, &trackCount_copy, trackCount,
												 &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN,
												 &T_copy, T,
												 &BED_copy, ED);
	
	
	
	trackCount->numPoints = W.num_pts*(num_new_linears);
	int solution_counter = 0;
	
	
	for (kk = 0; kk< num_new_linears; kk++)
	{
		if (solve_options->verbose_level>=1)
			printf("solving for linear %d\n",kk);
		
//		//set current linear in the evaluator data
		vec_cp_mp(ED->current_linear,new_linears[kk]);


		
		// track each of the start points
#ifdef _OPENMP
#pragma omp parallel for private(ii, oid, startPointIndex) schedule(runtime)
#endif
		for (ii = 0; ii < W.num_pts; ii++)
		{ // get current thread number
			oid = thread_num();
			
			if (solve_options->verbose_level>=1)
				printf("path %d\n",ii);
			
#ifdef printpathlintolin
			BED_copy[oid].num_steps = 0;
#endif
			
			
			// print the header of the path to OUT
//			printPathHeader_mp(OUT_copy[oid], &startPts[startPointIndex], &T_copy[oid], solution_counter, &BED_copy[oid], eval_func_mp);
			lin_to_lin_track_path_mp(solution_counter, &EG[oid], &startPts[ii], OUT_copy[oid], MIDOUT_copy[oid], &T_copy[oid], &BED_copy[oid], curr_eval_mp, change_prec, find_dehom); //curr_eval_d,
			

			
#ifdef printpathlintolin
			int mm;
			fprintf(BED_copy[oid].FOUT,"-100 %d ",BED_copy[oid].num_steps);
			for (mm=0; mm<BED_copy[oid].num_variables-1; ++mm) {
				fprintf(BED_copy[oid].FOUT,"0 0 ");
			}
			fprintf(BED_copy[oid].FOUT,"\n%d\n\n",EG->retVal);
#endif
			
			
			// check to see if it should be sharpened
			if (EG[oid].retVal == 0 && T_copy[oid].sharpenDigits > 0)
			{ // use the sharpener for after an endgame
				sharpen_endpoint_endgame(&EG[oid], &T_copy[oid], OUT_copy[oid], NULL, &BED_copy[oid], NULL, curr_eval_mp, NULL);
			}
			

			
			
			//get the terminal time in double form
			comp_d time_to_compare;
			if (EG->prec < 64) {
				set_d(time_to_compare,EG->PD_d.time);}
			else {
				mp_to_d(time_to_compare, EG->PD_mp.time); }
			
			
			int issoln = 0;
			if (EG->last_approx_prec!=1) {
//				if (EG->prec<64){
//					issoln = check_issoln_lintolin_d(&EG[oid],  &T_copy[oid], &BED_copy[oid]); }
//				else {
					issoln = check_issoln_lintolin_mp(&EG[oid], &T_copy[oid], &BED_copy[oid]);
//			}
			}
			else{
				
//				if (EG->prec<64){
//					print_point_to_screen_matlab(EG->PD_d.point,"solution");}
//				else {
					print_point_to_screen_matlab_mp(EG->PD_mp.point,"solution");
//				}
				
				printf("the last approximation was of precision %d\n",EG->last_approx_prec);
				printf("this is probably a problem\n");
				mypause();
			}
			
			
			
			if ((EG->retVal != 0 && time_to_compare->r > T->minTrackT) || !issoln) {  // <-- this is the real indicator of failure...
				trackCount->failures++;
				printf("\nretVal = %d\nthere was a fatal path failure tracking linear %d, witness point %d\n\n",EG->retVal,kk,ii);
				print_path_retVal_message(EG->retVal);
				print_point_to_screen_matlab_mp(BED_copy->old_linear,"old");
				print_point_to_screen_matlab_mp(BED_copy->current_linear,"new");
				exit(EG->retVal); //failure intolerable in this solver.
			}
			else
			{
				//otherwise converged, but may have still had non-zero retval due to other reasons.
				endgamedata_to_endpoint(&endPoints[solution_counter], EG);
				trackCount->successes++;
				solution_counter++; // probably this could be eliminated
			}
			
			
			

			
			
			
		}// re: for (ii=0; ii<W.num_pts ;ii++)
	} // for each new linear
	
	

	
	
	//clear the data structures.
	
  for (ii = 0; ii >W.num_pts; ii++)
  { // clear startPts[ii]
    clear_point_data_mp(&startPts[ii]);
  }
  free(startPts);
	
  // clear the structures
  clear_lintolin_omp_mp(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, &BED_copy);
	
  return;
}



// derived from zero_dim_track_path_d
void lin_to_lin_track_path_mp(int pathNum, endgame_data_t *EG_out,
														 point_data_mp *Pin,
														 FILE *OUT, FILE *MIDOUT, tracker_config_t *T,
														 void const *ED,
														 int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
														 int (*change_prec)(void const *, int),
														 int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: actually does the zero-dimensional tracking and sets   *
 *  up EG_out                                                    *
 \***************************************************************/
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




// derived from zero_dim_basic_setup_d
int lin_to_lin_setup_mp(FILE **OUT, char *outName,
											 FILE **midOUT, char *midName,
											 tracker_config_t *T,
											 lintolin_eval_data_mp *ED,
											 prog_t *dummyProg,
											 int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
											 int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
											 int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
											 char *preprocFile, char *degreeFile,
											 int findStartPts, char *pointsIN, char *pointsOUT,
											 mat_mp n_minusone_randomizer_matrix,
												witness_set W,
												solver_configuration *solve_options)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES: number of original variables                   *
 * NOTES: setup for zero dimensional tracking                    *
 \***************************************************************/
{ // need to create the homotopy
	
  int rank, patchType, ssType, numOrigVars, adjustDegrees, numGps;
	
  *eval_d = &lin_to_lin_eval_d;   
  *eval_mp = &lin_to_lin_eval_mp;
	
  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");
	

	
	
  // setup a straight-line program, using the file(s) created by the parser
  T->numVars = numOrigVars = setupProg_count(dummyProg, T->Precision, T->MPType, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);
	
	
  // setup preProcData
  setupPreProcData(preprocFile, &ED->preProcData);
	
	
	
  numGps = ED->preProcData.num_var_gp + ED->preProcData.num_hom_var_gp;
	
	
  // now that we know the rank, we can setup the rest of ED
  if (numGps == 1)// this should ALWAYS be the case in this solver.
  { // 1-hom
    patchType = 2; // 1-hom patch
    ssType = 0;    // with 1-hom, we use total degree start system
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
    setuplintolinEval_mp(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->Precision, &T->numVars, NULL, NULL, NULL, ED, adjustDegrees, n_minusone_randomizer_matrix, W,solve_options);
  }
  else
  { // m-hom, m > 1
		printf("there cannot be more than one variable group in bertini_real\n");
		exit(-11011);
  }
	
#ifdef printpathlintolin
	int ii;
	ED->FOUT = safe_fopen_append("pathtrack_lintolin");
	fprintf(ED->FOUT,"%d ",W.num_variables);
	for (ii=0; ii<W.num_variables; ii++) {
		fprintf(ED->FOUT,"%s ",W.variable_names[ii]);
	}
	fprintf(ED->FOUT,"\n%d ",W.num_pts);
	fprintf(ED->FOUT,"%d %d %d ",T->MPType, T->odePredictor, T->endgameNumber);
	fprintf(ED->FOUT,"\n");
#endif
	
	
  return numOrigVars;
}





//this derived from basic_eval_d
int lin_to_lin_eval_mp(point_mp funcVals, point_mp parVals, vec_mp parDer, mat_mp Jv, mat_mp Jp, point_mp current_variable_values, comp_mp pathVars, void const *ED)
{ // evaluates a special homotopy type, build for bertini_real

	//	print_comp_mp_matlab(pathVars,"pathvars");

	
  lintolin_eval_data_mp *BED = (lintolin_eval_data_mp *)ED; // to avoid having to cast every time
	


  int ii, jj;
  comp_mp one_minus_s, gamma_s; init_mp(one_minus_s); init_mp(gamma_s);
	
	comp_mp temp, temp2; init_mp(temp); init_mp(temp2);
	
	vec_mp one_minus_s_times_current_linear;
	init_vec_mp(one_minus_s_times_current_linear,BED->num_variables);
	one_minus_s_times_current_linear->size = BED->num_variables;
	
	vec_mp patchValues; init_vec_mp(patchValues, 0);
	
	vec_mp vars_times_curr_linear; init_vec_mp(vars_times_curr_linear,BED->num_variables);
	vec_mp vars_times_old_linear; init_vec_mp(vars_times_old_linear,BED->num_variables);
	vars_times_old_linear->size = vars_times_curr_linear->size = BED->num_variables;
	
	
	vec_mp temp_function_values; init_vec_mp(temp_function_values,0);
	vec_mp AtimesF; init_vec_mp(AtimesF,0);  // declare  // initialize
	vec_mp gamma_s_times_old_linear; init_vec_mp(gamma_s_times_old_linear,BED->num_variables);
	
	
	mat_mp temp_jacobian_functions; init_mat_mp(temp_jacobian_functions,0,0);
	mat_mp temp_jacobian_parameters; init_mat_mp(temp_jacobian_parameters,0,0);
	mat_mp Jv_Patch; init_mat_mp(Jv_Patch, 0, 0);
	mat_mp AtimesJ; init_mat_mp(AtimesJ,1,1);
	
	
	init_vec_mp(funcVals,BED->num_variables); funcVals->size = BED->num_variables;

	
	//set the sizes
	gamma_s_times_old_linear->size = BED->num_variables;
	
	
	
	

	
	change_size_mat_mp(Jv, BED->num_variables, BED->num_variables); // this is an input
	Jv->rows = Jv->cols = BED->num_variables;  //  -> this should be square!!!
	
	
  set_one_mp(one_minus_s);
  sub_mp(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
  mul_mp(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	

	//the SLP evaluation
	evalProg_mp(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	
	
  // evaluate the patch
  patch_eval_mp(patchValues, parVals, parDer, Jv_Patch, Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	

	// we assume that the only parameter is s = t and setup parVals & parDer accordingly.
	
  // set parVals & parDer correctly
	change_size_mat_mp(Jp, BED->num_variables, 1);
	Jp->rows = BED->num_variables; Jp->cols = 1;
	
  change_size_point_mp(parVals, 1);
  change_size_vec_mp(parDer, 1);
	
  parVals->size = parDer->size = 1;
  set_mp(&parVals->coord[0], pathVars); // s = t
  set_one_mp(&parDer->coord[0]);       // ds/dt = 1
	
  // combine everything
	
	
	
  
  
	
	
  
	
	
	
	
	mat_mul_mp(AtimesJ,BED->n_minusone_randomizer_matrix,temp_jacobian_functions);
	
	
	mul_mat_vec_mp(AtimesF,BED->n_minusone_randomizer_matrix, temp_function_values ); // set values of AtimesF (A is randomization matrix)
	
	


	for (ii=0; ii<AtimesF->size; ii++) { // for each function, after (real orthogonal) randomization
		set_mp(&funcVals->coord[ii], &AtimesF->coord[ii]);
	}
	
	
	
	

	
	// multiply vars times the new linear, with (1-s)
	for (ii=0; ii<BED->num_variables; ii++) { // for each variable, including the homogeneous ones.
		mul_mp(&one_minus_s_times_current_linear->coord[ii], &BED->current_linear->coord[ii], one_minus_s);
		mul_mp(&vars_times_curr_linear->coord[ii], &one_minus_s_times_current_linear->coord[ii], &current_variable_values->coord[ii]);
	}
	


	
	
	// multiply vars times the old linear, with gamma*s
	for (ii=0; ii<BED->num_variables; ii++) { // for each variable, including the homogeneous ones.
		mul_mp(&gamma_s_times_old_linear->coord[ii], &BED->old_linear->coord[ii], gamma_s);
		mul_mp(&vars_times_old_linear->coord[ii],&gamma_s_times_old_linear->coord[ii],&current_variable_values->coord[ii]);
	}
	
 
	// add the old and the new
	set_zero_mp(temp);  //  initialize to 0
	for (ii=0; ii<BED->num_variables; ii++) {  // for each variable, including the homogenizing ones.
		add_mp(temp2, &vars_times_old_linear->coord[ii], &vars_times_curr_linear->coord[ii]);
		add_mp(temp,temp,temp2); //  tested correct, does not overwrite values in temp before accessed for read.
	}
	
	
	
	//set the value of the linears we are homotoping.
	set_mp(&funcVals->coord[BED->num_variables-2],temp); // this is the entry for the linear's homotopy.

	
	
	//set the PATCH values
	int offset = BED->num_variables-1;
	for (ii=0; ii<BED->patch.num_patches; ii++) {
		set_mp(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
	}
	
	
	
	
	//////////////
	//
	// set the JACOBIAN values.
	//
	////////////////////
	
	//first, the entries related to the functions
	
  for (ii = 0; ii < BED->n_minusone_randomizer_matrix->rows; ii++)
  {
		for (jj = 0; jj < BED->num_variables; jj++)
		{
			set_mp(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
		}
  }
	
	
	
	
	
	// the entries for the linears we are homotoping
	
	
	offset = BED->num_variables-2;
	for (ii=0; ii<BED->num_variables; ii++) {
		add_mp(temp,&gamma_s_times_old_linear->coord[ii], &one_minus_s_times_current_linear->coord[ii]);
		set_mp(&Jv->entry[offset][ii], temp); // HERE DAB!!!
	}
	
	
	
	offset = BED->num_variables-1;
	for (ii=0; ii<BED->num_variables; ii++) {
		set_mp(&Jv->entry[offset][ii],&Jv_Patch->entry[0][ii]);
	}
	
	
	
	
	for (ii = 0; ii<BED->num_variables-2; ii++) {
		set_zero_mp(&Jp->entry[ii][0]);  // no parameter dependence means zero derivative for these functions
	}
	
	
	// Jp = -current_linear_times_vars + gamma*old_linear_times_vars
	set_zero_mp(&Jp->entry[BED->num_variables-2][0]);
	dot_product_mp(temp,BED->current_linear,current_variable_values);
	neg_mp(temp,temp);
	dot_product_mp(temp2,BED->old_linear,current_variable_values);
	mul_mp(temp2,temp2,BED->gamma);
	add_mp(&Jp->entry[BED->num_variables-2][0],temp,temp2);

	
	
	
	
	// the entries in the jacobian for the patch equations.
	offset = BED->num_variables-1;
	for (ii = 0; ii<BED->patch.num_patches; ii++)  // for each patch equation
	{ // funcVals = patchValues
		// Jp = 0
		set_zero_mp(&Jp->entry[ii+offset][0]);
		// Jv = Jv_Patch
		for (jj = 0; jj<BED->num_variables; jj++) // for each variable
		{
			set_mp(&Jv->entry[ii+offset][jj], &Jv_Patch->entry[ii][jj]);
		}
		
	}
	
	
	//done!  yay!
	
	//	printf("gamma = %lf+1i*%lf;\n", BED->gamma->r, BED->gamma->i);
	//

//	print_point_to_screen_matlab_mp(parVals,"parVals");
//	print_matrix_to_screen_matlab_mp( temp_jacobian_functions,"jac");
//	print_point_to_screen_matlab_mp(current_variable_values,"currvars");
//	print_point_to_screen_matlab_mp(BED->current_linear,"new");
//	print_point_to_screen_matlab_mp(BED->old_linear,"old");
//	print_point_to_screen_matlab_mp(funcVals,"F");
//	print_matrix_to_screen_matlab_mp(Jv,"Jv");
//	print_matrix_to_screen_matlab_mp(Jp,"Jp");
//	print_matrix_to_screen_matlab_mp(BED->n_minusone_randomizer_matrix,"n_minusone_randomizer_matrix");
//
	//
	
	
	clear_mp(one_minus_s);
	clear_mp(gamma_s);
	clear_mp(temp);
	clear_mp(temp2);
	clear_vec_mp(one_minus_s_times_current_linear);
	
	
	clear_vec_mp(patchValues);
	clear_vec_mp(vars_times_curr_linear);
	clear_vec_mp(vars_times_old_linear);
	clear_vec_mp(temp_function_values);
	clear_vec_mp(AtimesF);  // declare  // initialize
	clear_vec_mp(gamma_s_times_old_linear);
	
	
	clear_mat_mp(temp_jacobian_functions);
	clear_mat_mp(temp_jacobian_parameters);
	clear_mat_mp(Jv_Patch);
	clear_mat_mp(AtimesJ);
	
	
#ifdef printpathlintolin
	BED->num_steps++;
	vec_mp dehommed; init_vec_mp(dehommed,BED->num_variables-1); dehommed->size = BED->num_variables-1;
	dehomogenize_mp(&dehommed,current_variable_values);
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
	clear_vec_d(dehommed);
#endif

  return 0;
}







void lintolin_eval_clear_mp(lintolin_eval_data_mp *ED, int clearRegen, int MPType)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: clear ED                                               *
 \***************************************************************/
{

  patch_eval_data_clear_mp(&ED->patch);
  preproc_data_clear(&ED->preProcData);
	
	
	//SOME OTHER THINGS OUGHT TO BE CLEARED HERE  DAB
	
	
	//specifics for the lintolin method.
	clear_mp(ED->gamma);
	clear_vec_mp(ED->current_linear);
	clear_vec_mp(ED->old_linear);
	clear_mat_mp(ED->n_minusone_randomizer_matrix);
	
#ifdef printpathlintolin
	fclose(ED->FOUT);
#endif
	
  return;
}






void setup_lin_to_lin_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														FILE ***NONSOLN_copy, FILE *NONSOLN,
														tracker_config_t **T_copy, tracker_config_t *T,
														lintolin_eval_data_mp **BED_copy, lintolin_eval_data_mp *ED_mp)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: setup everything needed to do zero dimensional tracking*
 *  using OpenMP                                                 *
 \***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int ii;
  // error checking
  if (max_threads <= 0)
  {
    printf("\n\nERROR: The number of threads (%d) needs to be positive when setting up for tracking!\n", max_threads);
    bexit(ERROR_CONFIGURATION);
  }
	
	
// allocate space for EG
	*EG = (endgame_data_t *)bmalloc(max_threads * sizeof(endgame_data_t));
  for (ii = 0; ii < max_threads; ii++)
  {
    init_endgame_data(&(*EG)[ii], T->Precision);
  }
	
	
//  *EG = (endgame_data_t *)bmalloc(max_threads * sizeof(endgame_data_t));
//	
//	// initialize
//  for (ii = 0; ii < max_threads; ii++)
//	{
//    if (T->MPType == 2)
//    { // initialize for AMP tracking
//      init_endgame_data(&(*EG)[ii], 64);
//    }
//    else
//    { // initialize for double precision tracking
//      init_endgame_data(&(*EG)[ii], 52);
//    }
//	}
	
  // allocate space to hold pointers to the files
  *OUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *MIDOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *RAWOUT_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *FAIL_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
  *NONSOLN_copy = (FILE **)bmalloc(max_threads * sizeof(FILE *));
	
  if (max_threads == 1)
  { // setup the pointers
    *trackCount_copy = trackCount;
    *T_copy = T;
    *BED_copy = ED_mp;
//    (*BED_copy)->BED_mp = ED_mp; // make sure that this is pointed to inside of ED_d
		
    (*OUT_copy)[0] = OUT;
    (*RAWOUT_copy)[0] = RAWOUT;
    (*MIDOUT_copy)[0] = MIDOUT;
    (*FAIL_copy)[0] = FAIL;
    (*NONSOLN_copy)[0] = NONSOLN;
  }
  else // max_threads > 1
  { // allocate memory
		printf("more than one thread? prolly gonna crap out here in a sec\n");
		
    *trackCount_copy = (trackingStats *)bmalloc(max_threads * sizeof(trackingStats));
    *T_copy = (tracker_config_t *)bmalloc(max_threads * sizeof(tracker_config_t));
    *BED_copy = (lintolin_eval_data_mp *)bmalloc(max_threads * sizeof(lintolin_eval_data_mp));
		
    // copy T, ED_d, ED_mp, & trackCount
    for (ii = 0; ii < max_threads; ii++)
    { // copy T
      cp_tracker_config_t(&(*T_copy)[ii], T);
      // copy ED_d & ED_mp
      cp_lintolin_eval_data_mp(&(*BED_copy)[ii], ED_mp, T->MPType);
      // initialize trackCount_copy
      init_trackingStats(&(*trackCount_copy)[ii]);
      (*trackCount_copy)[ii].numPoints = trackCount->numPoints;
    }
		
    // setup the files
    char *str = NULL;
    int size;
    for (ii = 0; ii < max_threads; ii++)
    {
      size = 1 + snprintf(NULL, 0, "output_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", ii);
      (*OUT_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "midout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", ii);
      (*MIDOUT_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "rawout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", ii);
      (*RAWOUT_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "fail_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", ii);
      (*FAIL_copy)[ii] = fopen(str, "w+");
			
      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", ii);
      (*NONSOLN_copy)[ii] = fopen(str, "w+");
    }
    free(str);
  }
	
	
  return;
}

void clear_lintolin_omp_mp(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, lintolin_eval_data_mp **BED_copy)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: copy the relevant data back to the standard spot and   *
 *  clear the allocated data that was used by OpenMP             *
 \***************************************************************/
// if max_threads == 1, things are only pointers to the actual values,
// otherwise, they are copies
{
  int ii;
	
  // clear EG
  for (ii = max_threads - 1; ii >= 0; ii--)
  {
    clear_endgame_data(&(*EG)[ii]);
  }
  free(*EG);
	
  if (max_threads == 1)
  { // set the pointers to NULL since they just pointed to the actual values
    *trackCount_copy = NULL;
    *T_copy = NULL;
    *BED_copy = NULL;
		
    *OUT_copy[0] = NULL;
    *RAWOUT_copy[0] = NULL;
    *MIDOUT_copy[0] = NULL;
    *FAIL_copy[0] = NULL;
		
    // free the memory of the file pointers
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }
  else if (max_threads > 1)
  {
    // combine trackCount_copy
    add_trackingStats(trackCount, *trackCount_copy, max_threads);
		
    // clear the copies T, ED_d & ED_mp
    for (ii = max_threads - 1; ii >= 0; ii--)
    { // clear BED_copy - 0 since not using regeneration
      lintolin_eval_clear_mp(&(*BED_copy)[ii], 0, (*T_copy)[ii].MPType);
			// clear T_copy
      tracker_config_clear(&(*T_copy)[ii]);
    }
		
    // free the memory
    free(*trackCount_copy);
    free(*T_copy);
    free(*BED_copy);
		
    // copy all of the files to the appropriate place
    char ch, *str = NULL;
    int size;
    for (ii = 0; ii < max_threads; ii++)
    {
      // rewind to beginning for reading
      rewind((*OUT_copy)[ii]);
      // copy over to OUT
      ch = fgetc((*OUT_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(OUT, "%c", ch);
        ch = fgetc((*OUT_copy)[ii]);
      }
      // close file & delete
      fclose((*OUT_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "output_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "output_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*MIDOUT_copy)[ii]);
      // copy over to MIDOUT
      ch = fgetc((*MIDOUT_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(MIDOUT, "%c", ch);
        ch = fgetc((*MIDOUT_copy)[ii]);
      }
      // close file & delete
      fclose((*MIDOUT_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "midout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "midout_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*RAWOUT_copy)[ii]);
      // copy over to RAWOUT
      ch = fgetc((*RAWOUT_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(RAWOUT, "%c", ch);
        ch = fgetc((*RAWOUT_copy)[ii]);
      }
      // close file & delete
      fclose((*RAWOUT_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "rawout_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "rawout_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*FAIL_copy)[ii]);
      // copy over to FAIL
      ch = fgetc((*FAIL_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(FAIL, "%c", ch);
        ch = fgetc((*FAIL_copy)[ii]);
      }
      // close file & delete
      fclose((*FAIL_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "fail_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "fail_%d", ii);
      remove(str);
			
      // rewind to beginning for reading
      rewind((*NONSOLN_copy)[ii]);
      // copy over to NONSOLN
      ch = fgetc((*NONSOLN_copy)[ii]);
      while (ch != EOF)
      {
        fprintf(NONSOLN, "%c", ch);
        ch = fgetc((*NONSOLN_copy)[ii]);
      }
      // close file & delete
      fclose((*NONSOLN_copy)[ii]);
      size = 1 + snprintf(NULL, 0, "nonsolutions_%d", ii);
      str = (char *)brealloc(str, size * sizeof(char));
      sprintf(str, "nonsolutions_%d", ii);
      remove(str);
    }
    free(str);
    // free file memory
    free(*OUT_copy);
    free(*MIDOUT_copy);
    free(*RAWOUT_copy);
    free(*FAIL_copy);
    free(*NONSOLN_copy);
  }
	
  return;
}




void setuplintolinEval_mp(char preprocFile[], char degreeFile[], prog_t *dummyProg,
												 int squareSize, int patchType, int ssType, int prec,
												 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
												 lintolin_eval_data_mp *BED, int adjustDegrees,
												 mat_mp n_minusone_randomizer_matrix,
													witness_set W,
													solver_configuration *solve_options)
{
  int ii;
	int digits = prec_to_digits(mpf_get_default_prec());
  setupPreProcData(preprocFile, &BED->preProcData);

	setupPatch_mp(patchType, &BED->patch, digits, prec, ptr1, ptr2);
//

	
	BED->num_variables = W.num_variables;
	

	
	for (ii = 0; ii < BED->num_variables ; ii++)
		set_mp(&BED->patch.patchCoeff->entry[0][ii],&W.patch_mp[0]->coord[ii]);
	
	
	
	BED->SLP = dummyProg;

	
	init_mat_mp(BED->n_minusone_randomizer_matrix,0,0);
	mat_cp_mp(BED->n_minusone_randomizer_matrix,
					 n_minusone_randomizer_matrix);
	
	
	// set up the vectors to hold the two linears.
	
	
	init_vec_mp(BED->current_linear,0);
	change_size_vec_mp(BED->current_linear, n_minusone_randomizer_matrix->cols);
	BED->current_linear->size =  n_minusone_randomizer_matrix->cols;
	
	
	
	
	
	init_vec_mp(BED->old_linear,0);
	change_size_vec_mp(BED->old_linear, dummyProg->numVars);
	BED->old_linear->size =  dummyProg->numVars; // seriously, why is this not in the bertini method?
	

	vec_cp_mp(BED->old_linear,W.L_mp[0]);
	

	
	init_mp2(BED->gamma,prec);
	if (solve_options->use_gamma_trick==1)
		get_comp_rand_mp(BED->gamma); // set gamma to be random complex value
	else
		set_one_mp(BED->gamma);
	
	
  return;
}



void cp_lintolin_eval_data_mp(lintolin_eval_data_mp *BED, lintolin_eval_data_mp *BED_mp_input, int MPType)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: stores a copy of BED_(t)_input to BED                  *
 \***************************************************************/
{
	printf("entering cp_lintolin_eval_data_mp\nthis function is likely broken\n");
	exit(-1);
  cp_preproc_data(&BED->preProcData, &BED_mp_input->preProcData);
	//  cp_square_system_d(&BED->squareSystem, &BED_d_input->squareSystem);
  cp_patch_mp(&BED->patch, &BED_mp_input->patch);
	//  cp_start_system_d(&BED->startSystem, &BED_d_input->startSystem);
	BED->SLP = (prog_t *)bmalloc(1 * sizeof(prog_t));
  cp_prog_t(BED->SLP, BED_mp_input->SLP);
	
	
	//HERE COPY THE MATRICES  DAB !!!
	init_mat_mp(BED->n_minusone_randomizer_matrix,0,0);
	mat_cp_mp(BED->n_minusone_randomizer_matrix,BED_mp_input->n_minusone_randomizer_matrix);
	
	set_mp(BED->gamma, BED_mp_input->gamma);

	
  return;
}







int check_issoln_lintolin_d(endgame_data_t *EG,
													 tracker_config_t *T,
													 void const *ED)
{
  lintolin_eval_data_d *BED = (lintolin_eval_data_d *)ED; // to avoid having to cast every time
	
	
	int ii;
	
	
	
	
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
		clear_vec_d(terminal_pt);
	}
	else{
		evalProg_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_d.point, EG->PD_d.time, BED->SLP);
	}
	
	
	if (EG->last_approx_prec>=64) {
		vec_d prev_pt;  init_vec_d(prev_pt,1);
		vec_mp_to_d(prev_pt,EG->PD_mp.point);
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, prev_pt, EG->PD_d.time, BED->SLP);
		clear_vec_d(prev_pt);}
	else{
		evalProg_d(f, e.parVals, e.parDer, e.Jv, e.Jp, EG->last_approx_d, EG->PD_d.time, BED->SLP);
	}
	
	// compare the function values
	int isSoln = 1;
	for (ii = 0; (ii < BED->SLP->numFuncs) && isSoln; ii++)
	{
		n1 = d_abs_d( &e.funcVals->coord[ii]); // corresponds to final point
		n2 = d_abs_d( &f->coord[ii]); // corresponds to the previous point
		
		if (tol <= n1 && n1 <= n2)
		{ // compare ratio
			if (n1 > max_rat * n2){ // seriously what is the point of this?
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


int check_issoln_lintolin_mp(endgame_data_t *EG,
														tracker_config_t *T,
														void const *ED)
{
  lintolin_eval_data_mp *BED = (lintolin_eval_data_mp *)ED; // to avoid having to cast every time
	
	int ii;
	
	for (ii = 0; ii < T->numVars; ii++)
	{
    if (!(mpfr_number_p(EG->PD_mp.point->coord[ii].r) && mpfr_number_p(EG->PD_mp.point->coord[ii].i)))
		{
			printf("got not a number\n");
			print_point_to_screen_matlab_mp(EG->PD_mp.point,"bad solution");
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
	
	//this one guaranteed by entry condition
//	lin_to_lin_eval_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, ED);
	evalProg_mp(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, EG->PD_mp.point, EG->PD_mp.time, BED->SLP);
//	print_point_to_screen_matlab_mp(e.funcVals,"howfaroff");
	
	if (EG->last_approx_prec < 64)
	{ // copy to _mp
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





