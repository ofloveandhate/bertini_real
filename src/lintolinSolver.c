#include "lintolinSolver.h"




int lin_to_lin_solver_d(int MPType, double parse_time, unsigned int currentSeed,
												witness_set_d W,  // includes the initial linear.
												mat_d n_minusone_randomizer_matrix,  // for randomizing down to N-1 equations.
												vec_d *new_linears,   // collection of random complex linears.  for setting up the regeneration for V(f\\g)
												int num_new_linears,
												witness_set_d *W_new,
												int my_id, int num_processes, int headnode)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	
  FILE *OUT = NULL, *FAIL = fopen("failed_paths", "w"), *midOUT = NULL, *rawOUT = fopen("raw_data", "w");
  tracker_config_t T;
  prog_t dummyProg;
  bclock_t time1, time2;
  int num_variables = 0, userHom = 0, paramHom = 0, pathMod = 0, convergence_failures = 0, sharpening_failures = 0, sharpening_singular = 0, num_crossings = 0, num_sols = 0;
  int useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, supersetOnly = 0, usedEq = 0;
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  int (*ptr_to_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*ptr_to_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	
  lintolin_eval_data_d ED;  // was basic_eval_data_d  DAB
  trackingStats trackCount;
  char inputName[] = "func_input";
  double midpoint_tol, track_time, intrinsicCutoffMultiplier;
	
  bclock(&time1); // initialize the clock.
  init_trackingStats(&trackCount); // initialize trackCount to all 0
	
  // setup T
  setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &supersetOnly, &paramHom, MPType);
	
	T.MPType = 0;
	T.Precision = 53;
	
  // setup the precision structures
  initMP(T.Precision); // initialize MP based on T.Precision
	
#ifdef _OPENMP
#pragma omp parallel
#endif
  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
    initMP(T.Precision);
  }
	
	//  // call the setup function
	// setup for standard tracking - 'useRegen' is used to determine whether or not to setup 'start'
	num_variables = lin_to_lin_setup_d(&OUT, "output",
																		 &midOUT, "midpath_data",
																		 &T, &ED,
																		 &dummyProg,  //arg 7
																		 &startSub, &endSub, &startFunc, &endFunc,
																		 &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow,
																		 &ptr_to_eval_d, &ptr_to_eval_mp,  //args 17,18
																		 "preproc_data", "deg.out",
																		 !useRegen, "nonhom_start", "start",
																		 n_minusone_randomizer_matrix,W);
  
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
	
  if (T.MPType == 2)  //If we are doing adaptive precision path-tracking, we must set up AMP_eps, AMP_Phi, AMP_Psi based on config settings.
  {
    T.AMP_eps = (double) num_variables * num_variables;  //According to Demmel (as in the AMP paper), n^2 is a very reasonable bound for \epsilon.
    T.AMP_Phi = T.AMP_bound_on_degree*(T.AMP_bound_on_degree-1.0)*T.AMP_bound_on_abs_vals_of_coeffs;  //Phi from the AMP paper.
    T.AMP_Psi = T.AMP_bound_on_degree*T.AMP_bound_on_abs_vals_of_coeffs;  //Psi from the AMP paper.
    // initialize latest_newton_residual_mp to the maximum precision
    mpf_init2(T.latest_newton_residual_mp, T.AMP_max_prec);
  }
	
	
	
	
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
											 new_linears,
											 num_new_linears,
											 W_new,
											 FAIL, pathMod,
											 &T, &ED, ED.BED_mp,
											 ptr_to_eval_d, ptr_to_eval_mp,
											 change_prec, dehom);
	}
	
	

	
	fclose(midOUT);
	
	
	

	// finish the output to rawOUT
	fprintf(rawOUT, "%d\n\n", -1);  // bottom of rawOUT
	
	// check for path crossings
	midpoint_checker(trackCount.numPoints, num_variables, midpoint_tol, &num_crossings);
	
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
  printlintolinRelevantData(&ED, ED.BED_mp, T.MPType, usedEq, rawOUT);
	//	printZeroDimRelevantData(&ED, ED.BED_mp, T.MPType, usedEq, rawOUT);// legacy
	
  // close all of the files
  fclose(OUT);
  fclose(rawOUT);
  fprintf(FAIL, "\n");
  fclose(FAIL);
	
  // reproduce the input file needed for this run
  reproduceInputFile(inputName, "func_input", &T, 0, 0, currentSeed, pathMod, userHom, useRegen, regenStartLevel, maxCodim, specificCodim, intrinsicCutoffMultiplier, reducedOnly, supersetOnly, paramHom);
	
	//  // print the output
	//  if (userHom == -59)
	//  { // print the eq-by-eq output chart to the screen
	//    eqbyeqOutputChart_d(ED.EqD, stdout, T.regen_remove_inf);
	//  }
	
  // do the standard post-processing
  sort_points(num_crossings, &convergence_failures, &sharpening_failures, &sharpening_singular, inputName, num_sols, num_variables, midpoint_tol, T.final_tol_times_mult, &T, &ED.preProcData, useRegen == 1 && userHom == 0, userHom == -59);
	

	
  // print the failure summary
  printFailureSummary(&trackCount, convergence_failures, sharpening_failures, sharpening_singular);
	
//  // clear memory
//  if (userHom == 0 && paramHom == 0)
//  {
//		//    if (ED.squareSystem.Prog->numSubfuncs > 0)
//		//    { // clear subFuncsBelow
//		//      for (ii = ED.squareSystem.Prog->numFuncs - 1; ii >= 0; ii--)
//		//        free(subFuncsBelow[ii]);
//		//      free(subFuncsBelow);
//		//    }
    free(startSub);
    free(endSub);
    free(startFunc);
    free(endFunc);
    free(startJvsub);
    free(endJvsub);
    free(startJv);
    free(endJv);
//  }
	
	
  lintolin_eval_clear_d(&ED, userHom, T.MPType);
  tracker_config_clear(&T);
  clearMP();
	
  return 0;
}




//												vec_d *new_witness_points,
//												int *num_total_witness_points, // this will be set in this function

void lin_to_lin_track_d(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set_d W,
												vec_d *new_linears,
												int num_new_linears,
												witness_set_d *W_new,
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												lintolin_eval_data_d *ED_d,
												basic_eval_data_mp *ED_mp,
												int (*eval_func_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
												int (*eval_func_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
												int (*change_prec)(void const *, int),
												int (*find_dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *))
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
	
	
	
	// setup NONSOLN
	//  if (!ED_d->squareSystem.noChanges)
	//  { // setup NONSOLN
	//    NONSOLN = fopen("nonsolutions", "w");
	//    fprintf(NONSOLN, "                                    \n\n");
	//  }
	
	int (*curr_eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *) = NULL;
  int (*curr_eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *) = NULL;
	
	curr_eval_d = &lin_to_lin_eval_d;   //DAB
  curr_eval_mp = &basic_eval_mp; // DAB  // lol send them to the same place for now.
	

	
	point_data_d *startPts = NULL;
	startPts = (point_data_d *)bmalloc(W.W.num_pts * sizeof(point_data_d));
	
	
	
	for (ii = 0; ii < W.W.num_pts; ii++)
	{ // setup startPts[ii]
		init_point_data_d(&startPts[ii], W.num_variables); // also performs initialization on the point inside startPts
		change_size_vec_d(startPts[ii].point,W.num_variables);
		startPts[ii].point->size = W.num_variables;
		
		//NEED TO COPY IN THE WITNESS POINT
		
		//1 homogenize
//		set_one_d(&startPts[ii].point->coord[0]);
		
		//2 set the coordinates
		for (jj = 0; jj<W.num_variables; jj++) {
			startPts[ii].point->coord[jj].r = W.W.pts[ii]->coord[jj].r;
			startPts[ii].point->coord[jj].i = W.W.pts[ii]->coord[jj].i;
		}
		//3 set the start time to 1.
		set_one_d(startPts[ii].time);
	}
	
	
	T->endgameOnly = 0;
	
	
  // setup the rest of the structures
	endgame_data_t *EG = NULL;
  setup_lin_to_lin_omp_d(max,
												 &EG, &trackCount_copy, trackCount,
												 &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN,
												 &T_copy, T,
												 &BED_copy, ED_d, ED_mp);
	
	
	//initialize the structure for holding the produced data
	W_new->num_linears = (num_new_linears+1);
	W_new->L = (vec_d *)bmalloc((num_new_linears+1)*sizeof(vec_d));
	
	W_new->W.num_pts=W.W.num_pts*(num_new_linears+1);
  W_new->W.pts=(point_d *)bmalloc(W_new->W.num_pts*sizeof(point_d));
  W_new->num_variables = W.num_variables;
	
	
	
	int witness_point_counter = 0;
	
	
	for (ii=0; ii<W.W.num_pts; ii++) {
		init_vec_d(W_new->W.pts[ii],0);
		vec_cp_d(W_new->W.pts[ii],
						 startPts[ii].point)
		witness_point_counter++;
	}
	
	//	point_d* solutions_to_lintolin;
	//	solutions_to_lintolin = (point_d *)bmalloc( W.W.num_pts*num_new_linears * sizeof(point_d));
	
	
	init_vec_d(W_new->L[0],0);
	change_size_vec_d(W_new->L[0],W_new->num_variables);
	vec_cp_d(W_new->L[0],W.L[0]);
	
	
	for (kk = 0; kk< num_new_linears; kk++)
	{
#ifdef verbose
		printf("solving for linear %d\n",kk);
#endif
		
		//set current linear
		vec_cp_d(ED_d->current_linear,new_linears[kk]);
		
		init_vec_d(W_new->L[kk+1],0);
		change_size_vec_d(W_new->L[kk+1],W_new->num_variables);
		W_new->L[kk+1]->size = W_new->num_variables;
		vec_cp_d(W_new->L[kk+1],new_linears[kk]);
		
		
		// we pass the particulars of the information for this solve mode via the ED.
		
		
		// track each of the start points
#ifdef _OPENMP
#pragma omp parallel for private(ii, oid, startPointIndex) schedule(runtime)
#endif
		for (ii = 0; ii < W.W.num_pts; ii++)
		{ // get current thread number
			oid = thread_num();
			
			// print the path number if needed
//			if (pathMod > 0 && !(ii % pathMod))
//				printf("Tracking path %d of %d\n", ii, W.num_variables);
			
			//    if (eachStartPt)
			//    { // setup the next start point
			//      startPointIndex = 0;
			//      setupStart_d(&T_copy[oid], &startPts[startPointIndex], START);
			//    }
			//    else
			//    { // next start point is setup at index ii
			startPointIndex = ii;
			//    }
			
			// print the header of the path to OUT
			printPathHeader_d(OUT_copy[oid], &startPts[startPointIndex], &T_copy[oid], ii, &BED_copy[oid], eval_func_d);
			
			
			// track the path
			lin_to_lin_track_path_d(ii, &EG[oid], &startPts[startPointIndex], OUT_copy[oid], MIDOUT_copy[oid], &T_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, curr_eval_d, curr_eval_mp, change_prec, find_dehom);
			
			
			
			// check to see if it should be sharpened
			if (EG[oid].retVal == 0 && T_copy[oid].sharpenDigits > 0)
			{ // use the sharpener for after an endgame
				sharpen_endpoint_endgame(&EG[oid], &T_copy[oid], OUT_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, curr_eval_d, curr_eval_mp, change_prec);
			}
			
			
			init_vec_d(W_new->W.pts[witness_point_counter],0);
			change_size_vec_d(W_new->W.pts[witness_point_counter],W.num_variables);
			W_new->W.pts[witness_point_counter]->size = W.num_variables;

			
			vec_cp_d( W_new->W.pts[witness_point_counter],
							 EG->PD_d.point);
			
			comp_d result;
			dot_product_d(result, W_new->W.pts[witness_point_counter], ED_d->current_linear);
			
			
//			print_point_to_screen_matlab( W_new->W.pts[witness_point_counter],"newwitnesspoint");
//			print_point_to_screen_matlab( ED_d->current_linear,"newlinear");
//			
//			printf("dot_prod=%le+1i*%le;\n",result->r,result->i);
			
			
			if (EG->retVal!=0) {
				printf("retVal = %d\nthere was a path failure tracking linear %d, witness point %d\n",EG->retVal,kk,ii);
				print_path_retVal_message(EG->retVal);
				exit(EG->retVal);
			}

			
			witness_point_counter++;

			
			//				print_point_to_screen_matlab(startPts[startPointIndex].point,"start");
			//				print_point_to_screen_matlab(ED_d->old_linear,"old");
			//
			//				print_point_to_screen_matlab(EG->PD_d.point,"vend"); // the actual solution!!!
			//				print_point_to_screen_matlab(ED_d->current_linear,"curr");
			
			
			//    if (EG[oid].prec < 64)
			//    { // print footer in double precision
			//      printPathFooter_d(&trackCount_copy[oid], &EG[oid], &T_copy[oid], OUT_copy[oid], RAWOUT_copy[oid], FAIL_copy[oid], NONSOLN_copy[oid], &BED_copy[oid]);
			//    }
			//    else
			//    { // print footer in multi precision
			//      printPathFooter_mp(&trackCount_copy[oid], &EG[oid], &T_copy[oid], OUT_copy[oid], RAWOUT_copy[oid], FAIL_copy[oid], NONSOLN_copy[oid], BED_copy[oid].BED_mp);
			//    }
			
			
			//copy the solutions out of EG.
			
			
			
			
		}// re: for (ii=0; ii<W.W.num_pts ;ii++)
	} // for each new linear
	
	
	
	
	
	
	
	//clear the data structures.
	
  for (ii = 0; ii >W.W.num_pts; ii++)
  { // clear startPts[ii]
    clear_point_data_d(&startPts[ii]);
  }
  free(startPts);
	
  // clear the structures
  clear_lintolin_omp_d(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, &BED_copy);
	
	//  if (!ED_d->squareSystem.noChanges)
	//  { // complete NONSOLN
	//    rewind(NONSOLN);
	//    fprintf(NONSOLN, "%d", trackCount->junkCount);
	//    fclose(NONSOLN);
	//  }
	
	
	
	
	
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
											 mat_d n_minusone_randomizer_matrix,
											 witness_set_d W)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES: number of original variables                   *
 * NOTES: setup for zero dimensional tracking                    *
 \***************************************************************/
{ // need to create the homotopy
	printf("setup 569\n");
	
	
	
  int rank, patchType, ssType, numOrigVars, adjustDegrees, numGps;
	
  *eval_d = &lin_to_lin_eval_d;   //DAB
  *eval_mp = &basic_eval_mp; // DAB  // lol send them to the same place for now.
	
	//  *eval_mp = &lintolin_eval_mp; // DAB
	
  *OUT = fopen(outName, "w");  // open the main output files.
  *midOUT = fopen(midName, "w");
	
  if (T->MPType == 2) // using AMP - need to allocate space to store BED_mp
    ED->BED_mp = (basic_eval_data_mp *)bmalloc(1 * sizeof(basic_eval_data_mp));
  else
    ED->BED_mp = NULL;
	
	
  // setup a straight-line program, using the file(s) created by the parser
  T->numVars = numOrigVars = setupProg_count(dummyProg, T->Precision, T->MPType, startSub, endSub, startFunc, endFunc, startJvsub, endJvsub, startJv, endJv, subFuncsBelow);
	
	//	T.numVars = setupProg_count(&SLP, T.Precision, T.MPType, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow);
	
  // setup preProcData
  setupPreProcData(preprocFile, &ED->preProcData);
	
	
	
  numGps = ED->preProcData.num_var_gp + ED->preProcData.num_hom_var_gp;
  // find the rank
  rank = rank_finder_d(&ED->preProcData, dummyProg, T, T->numVars);
  // check to make sure that it is possible to have a zero dimensional component
  if (T->numVars > rank + numGps)
		//  {
		//    printf("The system has no zero dimensional solutions based on its rank!\n");
		//    printf("The rank of the system including the patches is %d while the total number of variables is %d.\n\n", rank + numGps, T->numVars);
		//    bexit(ERROR_INPUT_SYSTEM);
		//  }
		
		//AM I ACUTALLY SUPPOSED TO DO THIS? !!! HERE
		// adjust the number of variables based on the rank
		T->numVars = rank + numGps;
	
	
	
	
  // now that we know the rank, we can setup the rest of ED
  if (numGps == 1)
  { // 1-hom
    patchType = 2; // 1-hom patch
    ssType = 0;    // with 1-hom, we use total degree start system
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
    setuplintolinEval_d(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->MPType, &T->numVars, NULL, NULL, NULL, ED, adjustDegrees, n_minusone_randomizer_matrix, W);
  }
  else
  { // m-hom, m > 1
    patchType = 0; // random patch based on m-hom variable structure
    ssType = 1;    // with m-hom, we use the mhom structure for start system
    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
    setuplintolinEval_d(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->MPType, &ED->preProcData, NULL, NULL, NULL, ED, adjustDegrees, n_minusone_randomizer_matrix, W);
  }
	
	
	
	//  // create the start points, if needed
	//  if (findStartPts)
	//  {
	//    if (ssType == 0)
	//    { // setup total degree start points
	//      setupTD_startPoints_d(pointsIN, pointsOUT, ED->startSystem.size_r, ED->startSystem.degrees, &ED->patch);
	//    }
	//    else
	//    { // setup m-hom start points
	//      MHstartMaker_d(&ED->preProcData, ED->squareSystem.P, ED->startSystem.coeff, ED->patch.patchCoeff);
	//    }
	//  }
	
	
  return numOrigVars;
}





//this derived from basic_eval_d
int lin_to_lin_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, build for bertini_real
//	printf("entering eval\n");
	
	
//	printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);

  lintolin_eval_data_d *BED = (lintolin_eval_data_d *)ED; // to avoid having to cast every time
	
  int ii, jj;
  comp_d one_minus_s, gamma_s;
	vec_d patchValues;
	mat_d Jv_Patch;
	
  // we assume that the only parameter is s = t and setup parVals & parDer accordingly.


	init_vec_d(patchValues, 0);
	init_mat_d(Jv_Patch, 0, 0);
	
  set_one_d(one_minus_s);
  sub_d(one_minus_s, one_minus_s, pathVars);  // one_minus_s = (1 - s)
  mul_d(gamma_s, BED->gamma, pathVars);       // gamma_s = gamma * s
	
	
	vec_d temp_function_values;
	init_vec_d(temp_function_values,0);
	mat_d temp_jacobian_functions, temp_jacobian_parameters;
	init_mat_d(temp_jacobian_functions,0,0);
	init_mat_d(temp_jacobian_parameters,0,0);
	
	evalProg_d(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	
	
  // evaluate the patch
  patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, current_variable_values, pathVars, &BED->patch);  // Jp is ignored
	
	
	
  // set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
	
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
  // combine everything
	

	change_size_vec_d(funcVals,BED->num_variables);
  change_size_mat_d(Jv, BED->num_variables, BED->num_variables);
  change_size_mat_d(Jp, BED->num_variables, 1);
	
	
  funcVals->size = Jv->rows = Jp->rows = BED->num_variables;
  Jv->cols = BED->num_variables;  //  -> this should be square!!!
  Jp->cols = 1;
	
	mat_d AtimesJ;
	init_mat_d(AtimesJ,1,1);
	mat_mul_d(AtimesJ,BED->n_minusone_randomizer_matrix,temp_jacobian_functions);
	
	vec_d AtimesF;  // declare
	init_vec_d(AtimesF,0); // initialize
	mul_mat_vec_d(AtimesF,BED->n_minusone_randomizer_matrix, temp_function_values ); // set values of AtimesF (A is randomization matrix)
	
	
	vec_d vars_times_curr_linear, vars_times_old_linear;
	init_vec_d(vars_times_curr_linear,0);
	init_vec_d(vars_times_old_linear,0);
	change_size_vec_d(vars_times_curr_linear,BED->num_variables);
	change_size_vec_d(vars_times_old_linear,BED->num_variables);
	//set the sizes
	vars_times_old_linear->size = vars_times_curr_linear->size = BED->num_variables;
	
	
	for (ii=0; ii<AtimesF->size; ii++) { // for each function, after (real orthogonal) randomization
		set_d(&funcVals->coord[ii], &AtimesF->coord[ii]);
	}
	
	

	
	comp_d temp, temp2;
	vec_d one_minus_s_times_current_linear;
	init_vec_d(one_minus_s_times_current_linear,0);
	change_size_vec_d(one_minus_s_times_current_linear,BED->num_variables);
	one_minus_s_times_current_linear->size = BED->num_variables;
	
	// multiply vars times the new linear, with (1-s)
	for (ii=0; ii<BED->num_variables; ii++) { // for each variable, including the homogeneous ones.
		mul_d(&one_minus_s_times_current_linear->coord[ii], &BED->current_linear->coord[ii], one_minus_s);
		//		printf("%le+1i*%le = %le+1i*%le * %le+1i*%le\n",
		//					 temp->r, temp->i,
		//					 BED->current_linear->coord[ii].r,BED->current_linear->coord[ii].i,
		//					 one_minus_s->r, one_minus_s->i);
		mul_d(&vars_times_curr_linear->coord[ii], &one_minus_s_times_current_linear->coord[ii], &current_variable_values->coord[ii]);
	}
	
	//	print_point_to_screen_matlab(vars_times_curr_linear,"vars_times_curr_linear");
	
	vec_d gamma_s_times_old_linear;
	init_vec_d(gamma_s_times_old_linear,0);
	change_size_vec_d(gamma_s_times_old_linear,BED->num_variables);
	gamma_s_times_old_linear->size = BED->num_variables;
	
	
	// multiply vars times the old linear, with gamma*s
	for (ii=0; ii<BED->num_variables; ii++) { // for each variable, including the homogeneous ones.
		mul_d(&gamma_s_times_old_linear->coord[ii], &BED->old_linear->coord[ii], gamma_s);
		//		printf("%le+1i*%le = %le+1i*%le * %le+1i*%le\n",
		//					 temp->r, temp->i,
		//					 BED->old_linear->coord[ii].r,BED->old_linear->coord[ii].i,
		//					 gamma_s->r, gamma_s->i);
		mul_d(&vars_times_old_linear->coord[ii],&gamma_s_times_old_linear->coord[ii],&current_variable_values->coord[ii]);
	}
	
	// add the old and the new
	set_zero_d(temp);  //  initialize to 0
	for (ii=0; ii<BED->num_variables; ii++) {  // for each variable, including the homogenizing ones.
		add_d(temp2, &vars_times_old_linear->coord[ii], &vars_times_curr_linear->coord[ii]);
		add_d(temp,temp,temp2); //  tested correct, does not overwrite values in temp before accessed for read.
	}
	
	
	
	//set the value of the linears we are homotoping.
	set_d(&funcVals->coord[BED->num_variables-2],temp); // this is the entry for the linear's homotopy.
	
	
	
	
	
	//set the PATCH values
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
	

	

	
	// the entries for the linears we are homotoping
	
	
	offset = BED->num_variables-2;
	for (ii=0; ii<BED->num_variables; ii++) {
		add_d(temp,&gamma_s_times_old_linear->coord[ii], &one_minus_s_times_current_linear->coord[ii]);
		set_d(&Jv->entry[offset][ii], temp); // HERE DAB!!!
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
	
	for (ii=0; ii<BED->num_variables; ii++) {
		
		
		mul_d(temp,&BED->current_linear->coord[ii],&current_variable_values->coord[ii]);
		neg_d(temp,temp);
		
		
		mul_d(temp2,&BED->old_linear->coord[ii],&current_variable_values->coord[ii]);
		mul_d(temp2,BED->gamma,temp2);
		
		add_d(temp,temp,temp2); // add the two temps together.
		
		//add to the Jp entry
		add_d(&Jp->entry[BED->num_variables-2][0],&Jp->entry[BED->num_variables-2][0],temp);
	}
	
	
	
	
	
	
	
	
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
	

	
	
//	printf("gamma = %lf+1i*%lf;\n", BED->gamma->r, BED->gamma->i);
//
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
	//	print_matrix_to_screen_matlab(Jp,"Jp");
	
	
	//  clear_mat_d(Jv_F); clear_mat_d(Jv_start);
	clear_mat_d(Jv_Patch);
	//  clear_vec_d(F); clear_vec_d(startSys);
	clear_vec_d(patchValues);
	
//	printf("exiting eval\n");
  return 0;
}




void printlintolinRelevantData(lintolin_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP)
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
	//  square_system_eval_data_clear_d(&ED->squareSystem, MPType);
	//  start_system_eval_data_clear_d(&ED->startSystem);
  patch_eval_data_clear_d(&ED->patch);
  preproc_data_clear(&ED->preProcData);
	
	
	
  if (MPType == 2)
  { // ED->BED_mp->RD is just a pointer to ED->RD, so it will be cleared below - and thus 0 for not regen and 0 for not clear Prog since this was done above
    basic_eval_clear_mp(ED->BED_mp, 0, 0);
  }
	
	//  if (clearRegen == -59)
	//  {
	//    eqbyeq_clear_d(ED->EqD, MPType);
	//  }
	
	
	//specifics for the lintolin method.
	clear_vec_d(ED->current_linear);
	clear_vec_d(ED->old_linear);
	clear_mat_d(ED->n_minusone_randomizer_matrix);
	
	
	
  return;
}






void setup_lin_to_lin_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														FILE ***NONSOLN_copy, FILE *NONSOLN,
														tracker_config_t **T_copy, tracker_config_t *T,
														lintolin_eval_data_d **BED_copy, lintolin_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
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




void setuplintolinEval_d(char preprocFile[], char degreeFile[], prog_t *dummyProg,
												 int squareSize, int patchType, int ssType, int MPType,
												 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
												 lintolin_eval_data_d *BED, int adjustDegrees,
												 mat_d n_minusone_randomizer_matrix,
												 witness_set_d W)
{
  int ii;
	
  setupPreProcData(preprocFile, &BED->preProcData);
	//  setupSquareSystem_d(dummyProg, squareSize, &BED->preProcData, degreeFile, &BED->squareSystem, adjustDegrees); // NOTE: squareSystem must be setup before the patch!!!
  setupPatch_d(patchType, &BED->patch, ptr1, ptr2);
	
	
	BED->num_variables = W.num_variables;
	
	//  HERE!!!!!!!!!
//	set_zero_d(&BED->patch.patchCoeff->entry[0][0]);
	
	print_matrix_to_screen_matlab(BED->patch.patchCoeff,"stockpatch");
//	print_point_to_screen_matlab(
//	printf("%d\n",W.patch_size);
	
	for (ii = 0; ii < BED->num_variables ; ii++)
	{
		set_d(&BED->patch.patchCoeff->entry[0][ii],&W.patch[0]->coord[ii]);
	}

	print_matrix_to_screen_matlab(BED->patch.patchCoeff,"userpatch");


	vec_d patchValues;
	init_vec_d(patchValues,0);
	point_d parVals;
	init_vec_d(parVals,0);
	vec_d parDer;
	init_vec_d(parDer,0);
	mat_d Jv_Patch;
	init_mat_d(Jv_Patch,0,0);
	mat_d Jp;
	init_mat_d(Jp,0,0);
	comp_d pathVars;
	set_one_d(pathVars);
	patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W.W.pts[0], pathVars, &BED->patch);  // Jp is ignored
	
	print_point_to_screen_matlab(patchValues,"patchvalues");
	print_point_to_screen_matlab(W.W.pts[0],"initialpoint");
//	mypause();
	
	
	BED->SLP = dummyProg;
  // find the degrees of the 'square system'
	//  int *deg = (int *)bmalloc(BED->squareSystem.size_r * sizeof(int));
	//  for (ii = 0; ii < BED->squareSystem.size_r; ii++)
	//    deg[ii] = BED->squareSystem.orig_degrees[BED->squareSystem.P[ii]];
	
	//  setupStartSystem_d(ssType, BED->squareSystem.size_r, deg, BED->squareSystem.P, &BED->preProcData, degreeFile, &BED->startSystem);
	
  if (MPType == 2)
  { // using AMP - initialize using 16 digits & 64-bit precison
    int digits = 16, prec = 64;
		
    // setup preProcData
    setupPreProcData(preprocFile, &BED->BED_mp->preProcData);
    // setup the square system
		//    setupSquareSystem_d_to_mp(&BED->squareSystem, &BED->BED_mp->squareSystem, digits, prec);
    // setup the patch
    setupPatch_d_to_mp(&BED->patch, &BED->BED_mp->patch, digits, prec, patchType, dummyProg->numVars, ptr3, ptr4);
    // setup the start system
		//    setupStartSystem_d_to_mp(&BED->startSystem, &BED->BED_mp->startSystem, digits, prec);
  }
	
	init_mat_d(BED->n_minusone_randomizer_matrix,0,0);
	mat_cp_d(BED->n_minusone_randomizer_matrix,
					 n_minusone_randomizer_matrix);
	
	
	// set up the vectors to hold the two linears.
	
	
	init_vec_d(BED->current_linear,0);
	change_size_vec_d(BED->current_linear, n_minusone_randomizer_matrix->cols);
	BED->current_linear->size =  n_minusone_randomizer_matrix->cols;
	
	
	
	
	
	init_vec_d(BED->old_linear,0);
	change_size_vec_d(BED->old_linear, dummyProg->numVars);
	BED->old_linear->size =  dummyProg->numVars; // seriously, why is this not in the bertini method?
	
//	int ii;
//	for (ii = 0; ii<W.L->size; ii++) {
//		set_d(&BED->old_linear->coord[ii],&W.L->coord[ii]);
//	}
	vec_cp_d(BED->old_linear,W.L[0]);
//	printf("%d\n",BED->old_linear->size);
	
	
	
	get_comp_rand_d(BED->gamma); // set gamma to be random complex value
	
	
	
	
  return;
}



void cp_lintolin_eval_data_d(lintolin_eval_data_d *BED, lintolin_eval_data_d *BED_d_input, basic_eval_data_mp *BED_mp_input, int MPType)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: stores a copy of BED_(t)_input to BED                  *
 \***************************************************************/
{
	printf("entering cp lintolin\n");
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
    BED->BED_mp = (basic_eval_data_mp *)bmalloc(1 * sizeof(basic_eval_data_mp));
		
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



