

#include "linprodtodetjacSolver.h"




int linprod_to_detjac_solver_d(int MPType, double parse_time, unsigned int currentSeed,
												witness_set_d W,  // includes the initial linear.
												mat_d n_minusone_randomizer_matrix,  // for randomizing down to N-1 equations.
												vec_d projection,
												witness_set_d *W_new, // for passing the data back out of this function tree
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
	
  linprodtodetjac_eval_data_d ED;  // was basic_eval_data_d  DAB
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
	num_variables = linprod_to_detjac_setup_d(&OUT, "output",
																		 &midOUT, "midpath_data",
																		 &T, &ED,
																		 &dummyProg,  //arg 7
																		 &startSub, &endSub, &startFunc, &endFunc,
																		 &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow,
																		 &ptr_to_eval_d, &ptr_to_eval_mp,  //args 17,18
																		 "preproc_data", "deg.out",
																		 !useRegen, "nonhom_start", "start",
																		 n_minusone_randomizer_matrix,W,
																						projection);
  
	int (*change_prec)(void const *, int) = NULL;
	change_prec = &change_basic_eval_prec;
	
	int (*dehom)(point_d, point_mp, int *, point_d, point_mp, int, void const *, void const *) = NULL;
	dehom = &linprodtodetjac_dehom;
	
	
	
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
		printf("bertini_real not equipped to deal with endgameNumber 3\n");
		exit(-99);
	}
	else
	{ // use regular endgame
		linprod_to_detjac_track_d(&trackCount, OUT, rawOUT, midOUT,
											 W,  // was the startpts file pointer.
											 W_new,
											 FAIL, pathMod,
											 &T, &ED, ED.BED_mp,
											 ptr_to_eval_d, ptr_to_eval_mp,
											 change_prec, dehom);
	}
	
	
	
	
	fclose(midOUT);
	
	
	
	//    // delete the file containing the start points
	//    if (paramHom == 2 || userHom == 2)
	//      remove("start_param_hom");
	
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
  printlinprodtodetjacRelevantData(&ED, ED.BED_mp, T.MPType, usedEq, rawOUT);
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
	
  // clear memory
  if (userHom == 0 && paramHom == 0)
  {

    free(startSub);
    free(endSub);
    free(startFunc);
    free(endFunc);
    free(startJvsub);
    free(endJvsub);
    free(startJv);
    free(endJv);
  }
	
	
  linprodtodetjac_eval_clear_d(&ED, userHom, T.MPType);
  tracker_config_clear(&T);
  clearMP();
	
  return 0;
}




//												vec_d *new_witness_points,
//												int *num_total_witness_points, // this will be set in this function

void linprod_to_detjac_track_d(trackingStats *trackCount,
												FILE *OUT, FILE *RAWOUT, FILE *MIDOUT,
												witness_set_d W,
												witness_set_d *W_new,  // for holding the produced data.
												FILE *FAIL,
												int pathMod, tracker_config_t *T,
												linprodtodetjac_eval_data_d *ED_d,
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
	
  int ii,jj, oid, startPointIndex, max = max_threads();
  tracker_config_t *T_copy = NULL;
  linprodtodetjac_eval_data_d *BED_copy = NULL;
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
	
	curr_eval_d = &linprod_to_detjac_eval_d;   //DAB
  curr_eval_mp = &basic_eval_mp; // DAB  // lol send them to the same place for now.
	
	
	
	point_data_d *startPts = NULL;
	startPts = (point_data_d *)bmalloc(W.W.num_pts * sizeof(point_data_d));
	
	
	
	for (ii = 0; ii < W.W.num_pts; ii++)
	{ // setup startPts[ii]
		init_point_data_d(&startPts[ii], W.num_variables); // also performs initialization on the point inside startPts
		change_size_vec_d(startPts[ii].point,W.num_variables);
		startPts[ii].point->size = W.num_variables;
		
		//NEED TO COPY IN THE WITNESS POINT
		
		//1 do not homogenize
		
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
  setup_linprod_to_detjac_omp_d(max,
												 &EG, &trackCount_copy, trackCount,
												 &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN,
												 &T_copy, T,
												 &BED_copy, ED_d, ED_mp);
	
	
	//initialize the structure for holding the produced data
	W_new->num_linears = 1;
	W_new->L = (vec_d *)bmalloc(1*sizeof(vec_d));
	
	init_vec_d(W_new->L[0],0);
	change_size_vec_d(W_new->L[0],W_new->num_variables);
	W_new->L[0]->size = W_new->num_variables;
	vec_cp_d(W_new->L[0],W.L[0]);
	
	W_new->W.num_pts=0;  // initialize to 0.  will increment down below
  W_new->W.pts=(point_d *)bmalloc(W.W.num_pts*sizeof(point_d));
  W_new->num_variables = W.num_variables;
	
	
	int solution_counter = 0;

		
		// we pass the particulars of the information for this solve mode via the ED.
		
		
		// track each of the start points
#ifdef _OPENMP
#pragma omp parallel for private(ii, oid, startPointIndex) schedule(runtime)
#endif
//		for (ii = 0; ii < 1; ii++)
			for (ii = 0; ii < W.W.num_pts; ii++)
		{ // get current thread number
			oid = thread_num();
			

			startPointIndex = ii;
			//    }
			
			// print the header of the path to OUT
			printPathHeader_d(OUT_copy[oid], &startPts[startPointIndex], &T_copy[oid], ii, &BED_copy[oid], eval_func_d);
			
			
			// track the path
			linprod_to_detjac_track_path_d(ii, &EG[oid], &startPts[startPointIndex], OUT_copy[oid], MIDOUT_copy[oid], &T_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, curr_eval_d, curr_eval_mp, change_prec, find_dehom);
			
			
			
			// check to see if it should be sharpened
			if (EG[oid].retVal == 0 && T_copy[oid].sharpenDigits > 0)
			{ // use the sharpener for after an endgame
				sharpen_endpoint_endgame(&EG[oid], &T_copy[oid], OUT_copy[oid], &BED_copy[oid], BED_copy[oid].BED_mp, curr_eval_d, curr_eval_mp, change_prec);
			}
			
			
			

			
			if (EG->retVal!=0) {
				printf("retVal = %d\nthere was a path failure in linprod_to_detjac, tracking witness point %d\n",EG->retVal,ii);
				print_path_retVal_message(EG->retVal);
//				print_point_to_screen_matlab(EG->PD_d.point,"the_solution");
//				printf("%le\n",EG->t_val_at_latest_sample_point_d);
//				printf("%le\n",EG->condition_number);
//				mypause();
//				exit(EG->retVal);
			}
			else{
				//copy the solutions out of EG.
				
				init_vec_d(W_new->W.pts[solution_counter],0);
				change_size_vec_d(W_new->W.pts[solution_counter],W.num_variables);
				W_new->W.pts[solution_counter]->size = W.num_variables;
				vec_cp_d( W_new->W.pts[solution_counter],
								 EG->PD_d.point);
				solution_counter++;
			}
			

			
			//    if (EG[oid].prec < 64)
			//    { // print footer in double precision
			//      printPathFooter_d(&trackCount_copy[oid], &EG[oid], &T_copy[oid], OUT_copy[oid], RAWOUT_copy[oid], FAIL_copy[oid], NONSOLN_copy[oid], &BED_copy[oid]);
			//    }
			//    else
			//    { // print footer in multi precision
			//      printPathFooter_mp(&trackCount_copy[oid], &EG[oid], &T_copy[oid], OUT_copy[oid], RAWOUT_copy[oid], FAIL_copy[oid], NONSOLN_copy[oid], BED_copy[oid].BED_mp);
			//    }
			

			
			
		}// re: for (ii=0; ii<W.W.num_pts ;ii++)
	
	
	W_new->W.num_pts=solution_counter;
	printf("found %d solutions after the regeneration step\n",solution_counter);

	
	
	//clear the data structures.
	
  for (ii = 0; ii >W.W.num_pts; ii++)
  { // clear startPts[ii]
    clear_point_data_d(&startPts[ii]);
  }
  free(startPts);
	
  // clear the structures
  clear_linprodtodetjac_omp_d(max, &EG, &trackCount_copy, trackCount, &OUT_copy, OUT, &RAWOUT_copy, RAWOUT, &MIDOUT_copy, MIDOUT, &FAIL_copy, FAIL, &NONSOLN_copy, NONSOLN, &T_copy, &BED_copy);
	
	//  if (!ED_d->squareSystem.noChanges)
	//  { // complete NONSOLN
	//    rewind(NONSOLN);
	//    fprintf(NONSOLN, "%d", trackCount->junkCount);
	//    fclose(NONSOLN);
	//  }
	
	
	
	
	
  return;
}




// derived from zero_dim_track_path_d
void linprod_to_detjac_track_path_d(int pathNum, endgame_data_t *EG_out,
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
int linprod_to_detjac_setup_d(FILE **OUT, char *outName,
											 FILE **midOUT, char *midName,
											 tracker_config_t *T,
											 linprodtodetjac_eval_data_d *ED,
											 prog_t *dummyProg,
											 int **startSub, int **endSub, int **startFunc, int **endFunc, int **startJvsub, int **endJvsub, int **startJv, int **endJv, int ***subFuncsBelow,
											 int (**eval_d)(point_d, point_d, vec_d, mat_d, mat_d, point_d, comp_d, void const *),
											 int (**eval_mp)(point_mp, point_mp, vec_mp, mat_mp, mat_mp, point_mp, comp_mp, void const *),
											 char *preprocFile, char *degreeFile,
											 int findStartPts, char *pointsIN, char *pointsOUT,
											 mat_d n_minusone_randomizer_matrix,
											 witness_set_d W,
															vec_d projection)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES: number of original variables                   *
 * NOTES: setup for zero dimensional tracking                    *
 \***************************************************************/
{ // need to create the homotopy
	
	
	
  int rank, patchType, ssType, numOrigVars, adjustDegrees, numGps;
	
  *eval_d = &linprod_to_detjac_eval_d;   //DAB
  *eval_mp = &basic_eval_mp; // DAB  // lol send them to the same place for now.
	
	//  *eval_mp = &linprodtodetjac_eval_mp; // DAB
	
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
    setuplinprodtodetjacEval_d(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->MPType, &T->numVars, NULL, NULL, NULL, ED, adjustDegrees, n_minusone_randomizer_matrix, W,projection);
  }
  else
  { // m-hom, m > 1
		printf("more than one variable group. this is impossible in BertiniReal.\n");
		exit(-12321);
//    patchType = 0; // random patch based on m-hom variable structure
//    ssType = 1;    // with m-hom, we use the mhom structure for start system
//    adjustDegrees = 0; // if the system does not need its degrees adjusted, then that is okay
//    setuplinprodtodetjacEval_d(preprocFile, degreeFile, dummyProg, rank, patchType, ssType, T->MPType, &ED->preProcData, NULL, NULL, NULL, ED, adjustDegrees, n_minusone_randomizer_matrix, W,projection);
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
int linprod_to_detjac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d current_variable_values, comp_d pathVars, void const *ED)
{ // evaluates a special homotopy type, built for bertini_real
//	printf("t = %lf+1i*%lf;\n", pathVars->r, pathVars->i);
	
  linprodtodetjac_eval_data_d *BED = (linprodtodetjac_eval_data_d *)ED; // to avoid having to cast every time
	
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
	
	for (ii=0; ii<AtimesF->size; ii++) { // for each function, after (real orthogonal) randomization
		set_d(&funcVals->coord[ii], &AtimesF->coord[ii]);
	}
	
	
	
	
	comp_d temp, temp2;

	
	
	
	vec_d lin_func_vals;
	init_vec_d(lin_func_vals,BED->num_linears);
	lin_func_vals->size = BED->num_linears;
	
	// the product of the linears
	comp_d linprod, linprod_times_gamma_s;
	set_one_d(linprod); // initialize to 1 for multiplication
	for (ii=0; ii<BED->num_linears; ++ii) {
		
		dot_product_d(&lin_func_vals->coord[ii],BED->linears[ii],current_variable_values); // save into a buffer for calculating the derivative later.
		mul_d(linprod,linprod,&lin_func_vals->coord[ii]);// multiply.
	}
	//done with the loop over the linears.  now set the value
	
	mul_d( linprod_times_gamma_s,gamma_s,linprod); // sets the value linprod*gamma*s


	
	
	//the determinant of the jacobian, with the projection.
	
	mat_d tempmat;
	init_mat_d(tempmat,0,0);
	mat_cp_d(tempmat,AtimesJ);// copy into the matrix of which we will take the determinant
	
	increase_size_mat_d(tempmat,BED->num_variables,BED->num_variables); // make it bigger to accomodate more entries.  make square
	tempmat->rows = tempmat->cols = BED->num_variables; // change the size indicators
	
	//copy in the projection from BED
	for (ii=0; ii<BED->num_variables; ++ii) {
		set_d(&tempmat->entry[BED->num_variables-2][ii],&BED->projection->coord[ii]); // copy in the projection
	}
	
	//copy in the jacocian of the patch equation
	for (ii=0; ii<BED->num_variables; ++ii) {
		set_d(&tempmat->entry[BED->num_variables-1][ii],&Jv_Patch->entry[0][ii]); //  copy in the patch jacobian
	}
	
	//now TAKE THE DETERMINANT of tempmat.
	comp_d detjac;
	take_determinant(detjac,tempmat); // the determinant goes into detjac
	
	set_d(&funcVals->coord[BED->num_variables-2],detjac);
	mul_d(&funcVals->coord[BED->num_variables-2],&funcVals->coord[BED->num_variables-2],one_minus_s);//  (1-s)*detjac
	add_d(&funcVals->coord[BED->num_variables-2],&funcVals->coord[BED->num_variables-2],linprod_times_gamma_s); 
	

	

	
	//set the PATCH values
	int offset = BED->num_variables-1;
	for (ii=0; ii<BED->patch.num_patches; ii++) {
		set_d(&funcVals->coord[ii+offset], &patchValues->coord[ii]);
	}
	
	
	
	
	
	//now for the numerical derivative.
	
	mat_d Jv_detjac;
	init_mat_d(Jv_detjac,0,0);
	detjac_numerical_derivative_d(Jv_detjac, // return value
																current_variable_values, pathVars, BED->projection, BED); // input parameters for the method
	
	
	
	//////////////
	//
	// SET THE FINAL RETURNED JACOBIAN ENTRIES.
	//
	////////////////////
	// note: matrix was sized appropriately above.
	
	
	//first, the entries related to the functions
  for (ii = 0; ii < BED->num_variables-2; ii++)
  {
		for (jj = 0; jj < BED->num_variables; jj++)
		{
			set_d(&Jv->entry[ii][jj],&AtimesJ->entry[ii][jj]);
		}
  }
	
	
	
	
	
// the entries in the jacobian for the linprod and detjac terms.
	
	
	//&&&&&&&&&
	
	vec_d linprod_derivative;
	init_vec_d(linprod_derivative,BED->num_variables);
	linprod_derivative->size = BED->num_variables;
	
	
	// an implementation of the product rule
	comp_d running_prod;
	int kk;
	for (kk=0; kk<BED->num_variables; kk++) {
		set_zero_d(&linprod_derivative->coord[kk]); // initialize to 0 for the sum
		for (ii=0; ii<BED->num_linears; ++ii) { //  for each linear
			
			set_one_d(running_prod);// initialize to 1 for the product
			for (jj=0; jj<BED->num_linears; jj++) {
				if (jj==ii) {
					mul_d(running_prod,running_prod,&BED->linears[ii]->coord[kk]);// the coordinate IS the derivative of the linear function
				}
				else{
					mul_d(running_prod,running_prod,&lin_func_vals->coord[jj]); // the linear evaluated at curr_var_vals
				}
			}//re: jj
			
			add_d(&linprod_derivative->coord[kk],&linprod_derivative->coord[kk],running_prod);
		}// re:ii
		
	}//re: kk
	
	
	offset = BED->num_variables-2;
	for (ii=0; ii<BED->num_variables; ii++) {
		mul_d(temp,&linprod_derivative->coord[ii],gamma_s); // temp = d/dx linprod * gamma * s
		mul_d(temp2,&Jv_detjac->entry[0][ii],one_minus_s);
		add_d(&Jv->entry[offset][ii],temp, temp2);  // Jv = (1-s)d/dx(detjac) + s*gamma*d/dx linprod
	}
	//&&&&&&&&&&&
	
	
	// the last entries come from the patch equation(s)
	offset = BED->num_variables-1;
	for (ii=0; ii<BED->num_variables; ii++) {
		set_d(&Jv->entry[offset][ii],&Jv_Patch->entry[0][ii]);
	}
	
	
	
	
	for (ii = 0; ii<BED->num_variables-2; ii++) {
		set_zero_d(&Jp->entry[ii][0]);  // no parameter dependence means zero derivative for these functions
	}
	

	//initialize
		//Jp[N-2][0] = -detjac + gamma*linprod
	
	set_d(&Jp->entry[BED->num_variables-2][0],detjac);
	neg_d(&Jp->entry[BED->num_variables-2][0],&Jp->entry[BED->num_variables-2][0]);
	
	mul_d(temp,BED->gamma,linprod);

	add_d(&Jp->entry[BED->num_variables-2][0],&Jp->entry[BED->num_variables-2][0],temp);
	
	

	
	
	
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
	
	
	
	
	// set parVals & parDer correctly
  change_size_point_d(parVals, 1);
  change_size_vec_d(parDer, 1);
	
  parVals->size = parDer->size = 1;
  set_d(&parVals->coord[0], pathVars); // s = t
  set_one_d(&parDer->coord[0]);       // ds/dt = 1
	
	
	
	

//	print_matrix_to_screen_matlab( AtimesJ,"jac");
	//	print_point_to_screen_matlab(current_variable_values,"currvars");
	
	
//	print_matrix_to_screen_matlab(Jv_detjac,"Jv_detjac");
//	print_point_to_screen_matlab(linprod_derivative,"linprod_derivative");
//	print_point_to_screen_matlab(funcVals,"F");
//	print_point_to_screen_matlab(parVals,"parVals");
//	print_point_to_screen_matlab(parDer,"parDer");
//	print_matrix_to_screen_matlab(Jv,"Jv");
//	print_matrix_to_screen_matlab(Jp,"Jp");
	
	
	//these values are set in this function:  point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp
	//	print_matrix_to_screen_matlab(BED->n_minusone_randomizer_matrix,"n_minusone_randomizer_matrix");
	//
//	mypause();
	//

	

	clear_mat_d(temp_jacobian_functions);
	clear_mat_d(temp_jacobian_parameters);
	clear_mat_d(AtimesJ);clear_mat_d(tempmat);
	clear_mat_d(Jv_detjac);
	clear_mat_d(Jv_Patch);
	
	clear_vec_d(patchValues);
	clear_vec_d(temp_function_values);
	clear_vec_d(AtimesF);
	clear_vec_d(lin_func_vals);
	
	
	
//	printf("exiting eval\n");
  return 0;
}


int take_determinant(comp_d determinant, mat_d source_matrix){
	
	
	
	if (source_matrix->cols!=source_matrix->rows) {
		printf("source matrix is not square! (%d rows, %d columns)\n",source_matrix->rows,source_matrix->cols);
		exit(-108);
	}
	if (source_matrix->cols==0) {
		printf("source matrix has 0 entries!");
		exit(-109);
	}

	int num_variables = source_matrix->cols;
	int ii;
	
	mat_d intermediate;
	init_mat_d(intermediate,0,0);
	
	vec_d zerovec;
	init_vec_d(zerovec,0);
	change_size_vec_d(zerovec,num_variables);
	zerovec->size = num_variables;
	for (ii=0; ii<num_variables; ii++) {
		set_zero_d(&zerovec->coord[ii]);
	}
	
	int *rwnm = NULL;
	rwnm = NULL;
	vec_d garbage;
	init_vec_d(garbage,0);
	
	double tol = 1e-14; //  these should be for realsies
	double largeChange = 1e11;
	
	// returns x, intermediate.
	
	//	print_matrix_to_screen_matlab(tempmat,"tempmat");
	int sign;
	
	int retval = LU_matrixSolve_d(garbage, intermediate, &rwnm, &sign, source_matrix, zerovec,tol,largeChange);
	//the solution is in intermediate.
	//error check.  solution failed if retval!=0
	if (retval!=0) {
		printf("LU decomposition failed\n");
		print_matrix_to_screen_matlab(source_matrix,"source_matrix");
		exit(retval);
	}
	
	
	
	
	//compute the determinant
	set_one_d(determinant); // initialize
	for (ii=0; ii<num_variables; ii++) {
		mul_d(determinant,determinant,&intermediate->entry[rwnm[ii]][ii]);
	}
	determinant->r = determinant->r*sign;
	determinant->i = determinant->i*sign;
	
//	print_matrix_to_screen_matlab(source_matrix,"detme");
//	printf("candidate=%lf+1i*%lf;det(detme)\n",determinant->r,determinant->i);
//	mypause();
	// this verified correct via 20 samples in matlab.  dab.
	return 0;
}


// output: Jv.
// input: current_variable_values, pathvars, ED
int detjac_numerical_derivative_d(mat_d Jv, //  the returned value
																	point_d current_variable_values, comp_d pathVars, vec_d projection, void const *ED) // inputs
{
	int ii,jj,kk,mm; // counters
	linprodtodetjac_eval_data_d *BED = (linprodtodetjac_eval_data_d *)ED; // cast the ED
	change_size_mat_d(Jv, 1, BED->num_variables); // one row, num_variables columns
	Jv->rows = 1;
	Jv->cols = BED->num_variables;
	
	//initialize some containers, for the unused stuff from the called evaluators.
	point_d unused_function_values, unused_parVals; init_vec_d(unused_function_values,0);init_vec_d(unused_parVals,0);
	vec_d unused_parDer; init_vec_d(unused_parDer,0);
	mat_d unused_Jp; init_mat_d(unused_Jp,0,0);
	
	
	//set up the perturbed jacobians
	mat_d perturbed_forward_Jv, perturbed_forward_Jv_Patch, perturbed_backward_Jv, perturbed_backward_Jv_Patch; //, base_Jv, base_Jv_Patch
	init_mat_d(perturbed_backward_Jv,0,0);init_mat_d(perturbed_forward_Jv,0,0);
	init_mat_d(perturbed_backward_Jv_Patch,0,0);init_mat_d(perturbed_forward_Jv_Patch,0,0);
//	init_mat_d(base_Jv_Patch,0,0);init_mat_d(base_Jv,0,0);
	//the sizes for these variables will be set in the evalmethod calls below.  only need to init them here.
	
	//initialize the jacobians we will work with.
	point_d perturbed_forward_variables, perturbed_backward_variables;
	init_vec_d(perturbed_forward_variables,0); init_vec_d(perturbed_backward_variables,0);
	change_size_vec_d(perturbed_forward_variables,BED->num_variables);
	change_size_vec_d(perturbed_backward_variables,BED->num_variables);
	perturbed_backward_variables->size = perturbed_forward_variables->size = BED->num_variables;
	
	comp_d perturbation;
	
	perturbation->r = PERTURBATION_VALUE;
	perturbation->i = PERTURBATION_VALUE;
	
	mat_d AtimesJ;
	init_mat_d(AtimesJ,1,1); // size will be properly set later.
	
	comp_d det_backward, det_forward;
//	//get the baseline values at the current variable values.
//	evalProg_d(  unused_function_values, unused_parVals, unused_parDer,  //  unused output
//						 base_Jv,  // <---- the output we need
//						 unused_Jp, //unused output
//						 current_variable_values, pathVars, BED->SLP); // input
//	
//	patch_eval_d(unused_function_values, unused_parVals, unused_parDer, //  unused output
//							 base_Jv_Patch, // <---- the output we need
//							 unused_Jp, // unused output
//							 current_variable_values, pathVars, &BED->patch);  // input
	
	
	mat_d perturbed_forward_detme,perturbed_backward_detme; // create matrices
	init_mat_d(perturbed_forward_detme,BED->num_variables,BED->num_variables);
	init_mat_d(perturbed_backward_detme,BED->num_variables,BED->num_variables);
	
	perturbed_forward_detme->rows = perturbed_forward_detme->cols = perturbed_backward_detme->rows = perturbed_backward_detme->cols = BED->num_variables; // set the size indicators
	
	//for each variable, we need the derivative.  we will put them in the Jv matrix.
	for (ii=0; ii<BED->num_variables; ++ii) {
		
		//first, we assign the perturbed variables.  some of these calls could be eliminated
		for (jj=0; jj<BED->num_variables; ++jj) {
			set_d(&perturbed_forward_variables->coord[jj],&current_variable_values->coord[jj]);
		}
		for (jj=0; jj<BED->num_variables; ++jj) {
			set_d(&perturbed_backward_variables->coord[jj],&current_variable_values->coord[jj]);
		}
		add_d( &perturbed_forward_variables->coord[ii], &perturbed_forward_variables->coord[ii],perturbation);
		sub_d(&perturbed_backward_variables->coord[ii] ,&perturbed_backward_variables->coord[ii],perturbation);
		
		
		
		evalProg_d(unused_function_values, unused_parVals, unused_parDer,  //  unused output
							 perturbed_forward_Jv,  // <---- the output we need
							 unused_Jp, //unused output
							 perturbed_forward_variables, pathVars, BED->SLP); // input
		
		patch_eval_d(unused_function_values, unused_parVals, unused_parDer, //  unused output
								 perturbed_forward_Jv_Patch, // <---- the output we need
								 unused_Jp, // unused output
								 perturbed_forward_variables, pathVars, &BED->patch);  // input

		mat_mul_d(AtimesJ,BED->n_minusone_randomizer_matrix,perturbed_forward_Jv);
	
		for (kk=0; kk<BED->num_variables-2; ++kk) { // for each (randomized) equation
			for (mm=0; mm<BED->num_variables; ++mm) { //  for each variable
				set_d(&perturbed_forward_detme->entry[kk][mm],&AtimesJ->entry[kk][mm]);
			}
		}
		for (kk=0; kk<BED->num_variables;++kk) { // second to last row corresponds to the projection
			set_d(&perturbed_forward_detme->entry[BED->num_variables-2][kk],&projection->coord[kk]);
		}
		for (kk=0; kk<BED->num_variables;++kk) {// the last row corresponds to the patch
			set_d(&perturbed_forward_detme->entry[BED->num_variables-1][kk],&perturbed_forward_Jv_Patch->entry[0][kk]);
		}
		
		//now take the determinant;
		
		take_determinant(det_forward, perturbed_forward_detme);
		
		
		evalProg_d(  unused_function_values, unused_parVals, unused_parDer,  //  unused output
							 perturbed_backward_Jv,  // <---- the output we need
							 unused_Jp, //unused output
							 perturbed_backward_variables, pathVars, BED->SLP); // input
		
		patch_eval_d(unused_function_values, unused_parVals, unused_parDer, //  unused output
								 perturbed_backward_Jv_Patch, // <---- the output we need
								 unused_Jp, // unused output
								 perturbed_backward_variables, pathVars, &BED->patch);  // input
		//now have the pieces of the puzzle.
		
		mat_mul_d(AtimesJ,BED->n_minusone_randomizer_matrix,perturbed_backward_Jv);

		//copy in the jacobian WRT variables
		for (kk=0; kk<BED->num_variables-2; ++kk) { // for each (randomized) equation
			for (mm=0; mm<BED->num_variables; ++mm) { // for each variable
				set_d(&perturbed_backward_detme->entry[kk][mm],&AtimesJ->entry[kk][mm]);
			}
		}
		//copy in the projection
		for (kk=0; kk<BED->num_variables;++kk) {
			set_d(&perturbed_backward_detme->entry[BED->num_variables-2][kk],&projection->coord[kk]);
		}
		//copy in the patch
		for (kk=0; kk<BED->num_variables;++kk) {
			set_d(&perturbed_backward_detme->entry[BED->num_variables-1][kk],&perturbed_backward_Jv_Patch->entry[0][kk]);
		}
		
		//now take the determinant;
		
		take_determinant(det_backward, perturbed_backward_detme);
		
		sub_d(&Jv->entry[0][ii],det_forward,det_backward);  // Jv = (forward-backward)/(2h)
		div_d(&Jv->entry[0][ii],&Jv->entry[0][ii],perturbation);
		Jv->entry[0][ii].r /= 2.0;
		Jv->entry[0][ii].i /= 2.0;
		
//		print_matrix_to_screen_matlab(perturbed_forward_detme,"forwarddetme");
//		print_matrix_to_screen_matlab(perturbed_backward_detme,"backwarddetme");
//		printf("h=%lf+1i*%lf;\n",perturbation->r,perturbation->i);
//		printf("%lf+1i*%lf\n",Jv->entry[0][ii].r,Jv->entry[0][ii].i);
//		mypause();
	}

	clear_mat_d(AtimesJ);
	clear_mat_d(perturbed_backward_detme);clear_mat_d(perturbed_forward_detme);
	//HERE forgetting to clear some vectors.  losing memory.
	return 0;
}

void printlinprodtodetjacRelevantData(linprodtodetjac_eval_data_d *ED_d, basic_eval_data_mp *ED_mp, int MPType, int eqbyeqMethod, FILE *FP)
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



void linprodtodetjac_eval_clear_d(linprodtodetjac_eval_data_d *ED, int clearRegen, int MPType)
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
	
	
	//specifics for the linprodtodetjac method.
	clear_vec_d(ED->projection);
//	clear_vec_d(ED->old_linear);
	clear_mat_d(ED->n_minusone_randomizer_matrix);
	
	
	
  return;
}






void setup_linprod_to_detjac_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount,
														FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT,
														FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL,
														FILE ***NONSOLN_copy, FILE *NONSOLN,
														tracker_config_t **T_copy, tracker_config_t *T,
														linprodtodetjac_eval_data_d **BED_copy, linprodtodetjac_eval_data_d *ED_d, basic_eval_data_mp *ED_mp)
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
    *BED_copy = (linprodtodetjac_eval_data_d *)bmalloc(max_threads * sizeof(linprodtodetjac_eval_data_d));
		
    // copy T, ED_d, ED_mp, & trackCount
    for (ii = 0; ii < max_threads; ii++)
    { // copy T
      cp_tracker_config_t(&(*T_copy)[ii], T);
      // copy ED_d & ED_mp
      cp_linprodtodetjac_eval_data_d(&(*BED_copy)[ii], ED_d, ED_mp, T->MPType);
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

void clear_linprodtodetjac_omp_d(int max_threads, endgame_data_t **EG, trackingStats **trackCount_copy, trackingStats *trackCount, FILE ***OUT_copy, FILE *OUT, FILE ***RAWOUT_copy, FILE *RAWOUT, FILE ***MIDOUT_copy, FILE *MIDOUT, FILE ***FAIL_copy, FILE *FAIL, FILE ***NONSOLN_copy, FILE *NONSOLN, tracker_config_t **T_copy, linprodtodetjac_eval_data_d **BED_copy)
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
      linprodtodetjac_eval_clear_d(&(*BED_copy)[ii], 0, (*T_copy)[ii].MPType);
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




void setuplinprodtodetjacEval_d(char preprocFile[], char degreeFile[], prog_t *dummyProg,
												 int squareSize, int patchType, int ssType, int MPType,
												 void const *ptr1, void const *ptr2, void const *ptr3, void const *ptr4,
												 linprodtodetjac_eval_data_d *BED, int adjustDegrees,
												 mat_d n_minusone_randomizer_matrix,
												 witness_set_d W,
																vec_d projection)
{
  int ii;
	
  setupPreProcData(preprocFile, &BED->preProcData);
	//  setupSquareSystem_d(dummyProg, squareSize, &BED->preProcData, degreeFile, &BED->squareSystem, adjustDegrees); // NOTE: squareSystem must be setup before the patch!!!
  setupPatch_d(patchType, &BED->patch, ptr1, ptr2);
	
	
	BED->num_variables = W.num_variables;
	
	//  HERE!!!!!!!!!
	//	set_zero_d(&BED->patch.patchCoeff->entry[0][0]);
	
//	print_matrix_to_screen_matlab(BED->patch.patchCoeff,"stockpatch");
	//	print_point_to_screen_matlab(
	//	printf("%d\n",W.patch_size);
	
	for (ii = 0; ii < BED->num_variables ; ii++)
	{
		set_d(&BED->patch.patchCoeff->entry[0][ii],&W.patch[0]->coord[ii]);
	}
	
//	print_matrix_to_screen_matlab(BED->patch.patchCoeff,"userpatch");
	
	
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
//	
//	patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W.W.pts[0], pathVars, &BED->patch);  // Jp is ignored
	
//	print_point_to_screen_matlab(patchValues,"patchvalues");
//	print_point_to_screen_matlab(W.W.pts[0],"initialpoint");
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
	
	
	init_vec_d(BED->projection,0);
	change_size_vec_d(BED->projection, BED->num_variables);
	BED->projection->size = BED->num_variables;
	vec_cp_d(BED->projection,projection);
	
//	BED->current_linear->size =  n_minusone_randomizer_matrix->cols;
	
	
	BED->num_linears = W.num_linears;

	BED->linears = (vec_d *)malloc(W.num_linears*sizeof(vec_d));
	for (ii=0; ii<W.num_linears; ++ii) {
		init_vec_d(BED->linears[ii],0);
		change_size_vec_d(BED->linears[ii],W.L[ii]->size);
		BED->linears[ii]->size = W.L[ii]->size;
		vec_cp_d(BED->linears[ii], W.L[ii]);
	}
//	init_vec_d(BED->old_linear,0);
//	change_size_vec_d(BED->old_linear, dummyProg->numVars);
//	BED->old_linear->size =  dummyProg->numVars; // seriously, why is this not in the bertini method?
	
	//	int ii;
	//	for (ii = 0; ii<W.L->size; ii++) {
	//		set_d(&BED->old_linear->coord[ii],&W.L->coord[ii]);
	//	}
//	vec_cp_d(BED->old_linear,W.L[0]);
	//	printf("%d\n",BED->old_linear->size);
	
	
	
	get_comp_rand_d(BED->gamma); // set gamma to be random complex value
	
	
	
	
  return;
}



void cp_linprodtodetjac_eval_data_d(linprodtodetjac_eval_data_d *BED, linprodtodetjac_eval_data_d *BED_d_input, basic_eval_data_mp *BED_mp_input, int MPType)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: stores a copy of BED_(t)_input to BED                  *
 \***************************************************************/
{
	printf("entering cp_linprodtodetjac\n");
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




int linprodtodetjac_dehom(point_d out_d, point_mp out_mp, int *out_prec, point_d in_d, point_mp in_mp, int in_prec, void const *ED_d, void const *ED_mp)
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





//input : square matrix m
void determinant_jacobian_d(mat_d m, int * rank, double * det ){
	det= 0;  // initialize
	rank = 0;
	
	//will perform this via gaussian elimination and product of diagonal
	
	
	
	return;
}





//
//int detjac_to_detjac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
//
//
//
//int linprod_to_detjac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED);
//
//
//
////this derived from basic_eval_d
//int detjac_to_detjac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
//{ // evaluates the function as described in basic_eval_data_d: 'top' squareSystem*s + gamma*(1-s)*startSystem and 'bottom' patch
//	
//	printf("entering detjac_to_detjac_eval_d 568\n");
//	
//	
//  detjac_eval_data_d *BED = (detjac_eval_data_d *)ED; // to avoid having to cast every time
//  int i, j, size, top;
//  comp_d one_minus_s, gamma_s;
//  vec_d F, startSys, patchValues;
//  mat_d Jv_F, Jv_start, Jv_Patch;
//  // we assume that the only parameter is s = t and setup parVals & parDer accordingly. Also, F, startSys & patchValues do not depend on a parameter
//	
//	//this line verifies that we do indeed have the randomizer matrix
//	//	print_matrix_to_screen_matlab(BED->n_minusone_randomizer_matrix,"n_minusone_randomizer_matrix");
//	//	exit(-1);
//	
//  // initialize
//  init_vec_d(F, 0); init_vec_d(startSys, 0); init_vec_d(patchValues, 0);
//  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_start, 0, 0); init_mat_d(Jv_Patch, 0, 0);
//	
//  set_one_d(one_minus_s);
//  sub_d(one_minus_s, one_minus_s, pathVars);
//  mul_d(gamma_s, BED->startSystem.gamma, pathVars);
//	//	printf("the three evals 708\n");
//  // evalute the 'square' system F
//  square_system_eval_d(      F, parVals, parDer, Jv_F, Jp, vars, pathVars, &BED->squareSystem); // Jp is ignored
//	//	printf("F->size %d 711\n",F->size);
//  // evaluate the patch
//  patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);  // Jp is ignored
//	//	print_point_to_screen_matlab(patchValues, "patchValues");
//  // evaluate the start system
//	//  start_system_eval_d(startSys, parVals, parDer, Jv_start, Jp, vars, pathVars, &BED->startSystem); // Jp is ignored
//	//	printf("717\n");
//  // set parVals & parDer correctly
//  change_size_point_d(parVals, 1);
//  change_size_vec_d(parDer, 1);
//  parVals->size = parDer->size = 1;
//  set_d(&parVals->coord[0], pathVars); // s = t
//  set_one_d(&parDer->coord[0]);       // ds/dt = 1
//	
//  // combine everything
//	printf("gamma is %le+1i*%le\ncombining 726\n",BED->startSystem.gamma->r,BED->startSystem.gamma->i);
//  // find funcVals = 'top' F*(1-s) + gamma*s*startSys and 'bottom' patchValues
//  // find Jv = 'top' Jv_F*(1-s) + gamma*s*Jv_start and 'bottom' Jv_Patch
//  // find Jp = 'top' -F + gamma*startSys and 'bottom' all zeros
//  size = F->size + BED->patch.num_patches;
//  change_size_mat_d(Jv, size, vars->size);
//  change_size_mat_d(Jp, size, 1);
//  change_size_vec_d(funcVals, size);
//  funcVals->size = Jv->rows = Jp->rows = size;
//  top = F->size;
//  Jv->cols = vars->size;  // this should be square!!!
//  Jp->cols = 1;
//  for (i = 0; i < size; i++)
//  {
//    if (i < top)
//    { // funcVals = F*(1-s) + gamma*s*startSys
//      mul_d(&funcVals->coord[i], &F->coord[i], one_minus_s);
//      sum_mul_d(&funcVals->coord[i], &startSys->coord[i], gamma_s);
//      // Jp = -F + gamma*startSys
//      mul_d(&Jp->entry[i][0], BED->startSystem.gamma, &startSys->coord[i]);
//      sub_d(&Jp->entry[i][0], &Jp->entry[i][0], &F->coord[i]);
//      // Jv = Jv_F*(1-s) + gamma*s*Jv_start
//      for (j = 0; j < vars->size; j++)
//      {
//        mul_d(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_s);
//        sum_mul_d(&Jv->entry[i][j], &Jv_start->entry[i][j], gamma_s);
//      }
//    }
//    else
//    { // funcVals = patchValues
//      set_d(&funcVals->coord[i], &patchValues->coord[i - top]);
//      // Jp = 0
//      set_zero_d(&Jp->entry[i][0]);
//      // Jv = Jv_Patch
//      for (j = 0; j < vars->size; j++)
//      {
//        set_d(&Jv->entry[i][j], &Jv_Patch->entry[i - top][j]);
//      }
//    }
//  }
//	
//  clear_mat_d(Jv_F); clear_mat_d(Jv_start); clear_mat_d(Jv_Patch);
//  clear_vec_d(F); clear_vec_d(startSys); clear_vec_d(patchValues);
//	
//  return 0;
//}
//
//
//
//
//
//
//
////this derived from basic_eval_d
//int linprod_to_detjac_eval_d(point_d funcVals, point_d parVals, vec_d parDer, mat_d Jv, mat_d Jp, point_d vars, comp_d pathVars, void const *ED)
//{ // evaluates the function as described in basic_eval_data_d: 'top' squareSystem*s + gamma*(1-s)*startSystem and 'bottom' patch
//	
//	printf("entering linprod_to_detjac_eval_d 755\n");
//	
//	
//  detjac_eval_data_d *BED = (detjac_eval_data_d *)ED; // to avoid having to cast every time
//  int i, j, size, top;
//  comp_d one_minus_s, gamma_s;
//  vec_d F, startSys, patchValues;
//  mat_d Jv_F, Jv_start, Jv_Patch;
//  // we assume that the only parameter is s = t and setup parVals & parDer accordingly. Also, F, startSys & patchValues do not depend on a parameter
//	
//	//this line verifies that we do indeed have the randomizer matrix
//	//	print_matrix_to_screen_matlab(BED->n_minusone_randomizer_matrix,"n_minusone_randomizer_matrix");
//	//	exit(-1);
//	
//  // initialize
//  init_vec_d(F, 0); init_vec_d(startSys, 0); init_vec_d(patchValues, 0);
//  init_mat_d(Jv_F, 0, 0); init_mat_d(Jv_start, 0, 0); init_mat_d(Jv_Patch, 0, 0);
//	
//  set_one_d(one_minus_s);
//  sub_d(one_minus_s, one_minus_s, pathVars);
//  mul_d(gamma_s, BED->startSystem.gamma, pathVars);
//	//	printf("the three evals 708\n");
//  // evalute the 'square' system F
//  square_system_eval_d(      F, parVals, parDer, Jv_F, Jp, vars, pathVars, &BED->squareSystem); // Jp is ignored
//	//	printf("F->size %d 711\n",F->size);
//  // evaluate the patch
//  patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, vars, pathVars, &BED->patch);  // Jp is ignored
//	//	print_point_to_screen_matlab(patchValues, "patchValues");
//  // evaluate the start system
//	//  start_system_eval_d(startSys, parVals, parDer, Jv_start, Jp, vars, pathVars, &BED->startSystem); // Jp is ignored
//	//	printf("717\n");
//  // set parVals & parDer correctly
//  change_size_point_d(parVals, 1);
//  change_size_vec_d(parDer, 1);
//  parVals->size = parDer->size = 1;
//  set_d(&parVals->coord[0], pathVars); // s = t
//  set_one_d(&parDer->coord[0]);       // ds/dt = 1
//	
//  // combine everything
//	printf("gamma is %le+1i*%le\ncombining 726\n",BED->startSystem.gamma->r,BED->startSystem.gamma->i);
//  // find funcVals = 'top' F*(1-s) + gamma*s*startSys and 'bottom' patchValues
//  // find Jv = 'top' Jv_F*(1-s) + gamma*s*Jv_start and 'bottom' Jv_Patch
//  // find Jp = 'top' -F + gamma*startSys and 'bottom' all zeros
//  size = F->size + BED->patch.num_patches;
//  change_size_mat_d(Jv, size, vars->size);
//  change_size_mat_d(Jp, size, 1);
//  change_size_vec_d(funcVals, size);
//  funcVals->size = Jv->rows = Jp->rows = size;
//  top = F->size;
//  Jv->cols = vars->size;  // this should be square!!!
//  Jp->cols = 1;
//  for (i = 0; i < size; i++)
//  {
//    if (i < top)
//    { // funcVals = F*(1-s) + gamma*s*startSys
//      mul_d(&funcVals->coord[i], &F->coord[i], one_minus_s);
//      sum_mul_d(&funcVals->coord[i], &startSys->coord[i], gamma_s);
//      // Jp = -F + gamma*startSys
//      mul_d(&Jp->entry[i][0], BED->startSystem.gamma, &startSys->coord[i]);
//      sub_d(&Jp->entry[i][0], &Jp->entry[i][0], &F->coord[i]);
//      // Jv = Jv_F*(1-s) + gamma*s*Jv_start
//      for (j = 0; j < vars->size; j++)
//      {
//        mul_d(&Jv->entry[i][j], &Jv_F->entry[i][j], one_minus_s);
//        sum_mul_d(&Jv->entry[i][j], &Jv_start->entry[i][j], gamma_s);
//      }
//    }
//    else
//    { // funcVals = patchValues
//      set_d(&funcVals->coord[i], &patchValues->coord[i - top]);
//      // Jp = 0
//      set_zero_d(&Jp->entry[i][0]);
//      // Jv = Jv_Patch
//      for (j = 0; j < vars->size; j++)
//      {
//        set_d(&Jv->entry[i][j], &Jv_Patch->entry[i - top][j]);
//      }
//    }
//  }
//	
//  clear_mat_d(Jv_F); clear_mat_d(Jv_start); clear_mat_d(Jv_Patch);
//  clear_vec_d(F); clear_vec_d(startSys); clear_vec_d(patchValues);
//	
//  return 0;
//}
//
//
//
//
//
//
