#include "solver.hpp"




void get_tracker_config(solver_configuration *solve_options,int MPType)
{
	
	//necessary for the setupConfig call
	double intrinsicCutoffMultiplier;
	int userHom = 0, useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, pathMod = 0, reducedOnly = 0, supersetOnly = 0, paramHom = 0;
	//end necessaries for the setupConfig call.
	
	
  setupConfig(&solve_options->T, &solve_options->midpoint_tol, &userHom, &useRegen, &regenStartLevel, &maxCodim, &specificCodim, &pathMod, &intrinsicCutoffMultiplier, &reducedOnly, &supersetOnly, &paramHom, MPType);
	
	return;
}




void solver_init_config(solver_configuration *options){
	options->allow_multiplicity = 0;
	options->allow_singular = 0;
	options->allow_infinite = 0;
	options->allow_unsuccess = 0;
	options->use_midpoint_checker = 1;
	options->show_status_summary = 0;
	options->verbose_level = 0; // default to 0.  higher is more verbose
	options->use_gamma_trick = 0;
	
	options->complete_witness_set = 1;
}


void solver_clear_config(solver_configuration *options){
	//has no fields which require clearing.
	return;
}




/**
 sets the start_pts structure to hold all points in W
 
 \param startPts the value being set.  should be NULL input.
 \param W the witness_set input
 
 */
void generic_set_start_pts(point_data_d ** startPts,
							 witness_set & W)
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
													 witness_set & W)
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


void generic_track_path_d(int pathNum, endgame_data_t *EG_out,
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




void generic_track_path_mp(int pathNum, endgame_data_t *EG_out,
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

