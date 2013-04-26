#include "postProcessing.h"



void endgamedata_to_endpoint(post_process_t *endPoint, endgame_data_t *EG){
	
	int num_vars, ii;
	endPoint->path_num = EG->pathNum;
	
	endPoint->sol_prec = EG->prec;
	endPoint->cond_est = EG->condition_number;
	endPoint->final_t = EG->t_val_at_latest_sample_point_d;//???
	endPoint->first_increase = EG->first_increase;
	
	
	if (EG->prec==52) {
		num_vars = EG->PD_d.point->size;
		endPoint->sol_d  = (comp_d *)bmalloc(num_vars * sizeof(comp_d));
		endPoint->sol_mp = NULL;
		
		for (ii=0; ii<num_vars; ii++) {
			endPoint->sol_d[ii]->r = EG->PD_d.point->coord[ii].r;
			endPoint->sol_d[ii]->i = EG->PD_d.point->coord[ii].i;
		}
		
		endPoint->size_sol = num_vars;
		endPoint->function_resid_d = EG->function_residual_d;  // the function residual
		endPoint->newton_resid_d = EG->latest_newton_residual_d;
		endPoint->cycle_num = EG->PD_d.cycle_num;
		endPoint->accuracy_estimate = EG->error_at_latest_sample_point_d;//
		
	}
	else
	{
		num_vars = EG->PD_mp.point->size;
		
		endPoint->sol_d  = NULL;
		endPoint->sol_mp = (comp_mp *)bmalloc(num_vars * sizeof(comp_mp));
		for (ii=0; ii<num_vars; ii++) {
			init_mp2(endPoint->sol_mp[ii],EG->prec);
			mpf_set(endPoint->sol_mp[ii]->r,EG->PD_mp.point->coord[ii].r);
			mpf_set(endPoint->sol_mp[ii]->i,EG->PD_mp.point->coord[ii].i);
		}
		endPoint->size_sol = num_vars;
		mpf_init2(endPoint->function_resid_mp, EG->prec); mpf_init2(endPoint->newton_resid_mp, EG->prec);
		
		mpf_set(endPoint->function_resid_mp,EG->function_residual_mp); //this is undoubtedly incorrect
		mpf_set(endPoint->newton_resid_mp,EG->latest_newton_residual_mp);
		endPoint->cycle_num = EG->PD_mp.cycle_num;
		endPoint->accuracy_estimate = mpf_get_d(EG->error_at_latest_sample_point_mp);
	}
	
	
	if (EG->retVal==0) {
		endPoint->success = 1;
	}
	else if (EG->retVal == retVal_sharpening_failed){
		endPoint->success = retVal_sharpening_failed;
	}
	else if (EG->retVal == retVal_sharpening_singular_endpoint){
		endPoint->success = retVal_sharpening_singular_endpoint;
	}
	else{
		endPoint->success = -1;
	}
	
	
	
	endPoint->sol_num = 0; // set up for post-processing
	endPoint->multiplicity = 1;
	endPoint->isFinite = 0;
	//	// the post_process_t structure is used in post-processing //
	//	typedef struct
	//	{
	//		int path_num;     // path number of the solution
	//		int sol_num;      // solution number
	//		comp_d  *sol_d;   // solution
	//		comp_mp *sol_mp;
	//		int sol_prec;     // precision of the solution
	//		int size_sol;     // the number of entries in sol
	//		double function_resid_d;  // the function residual
	//		mpf_t  function_resid_mp;
	//		double cond_est;  // the estimate of the condition number
	//		double newton_resid_d;    // the newton residual
	//		mpf_t  newton_resid_mp;
	//		double final_t;   // the final value of time
	//		double accuracy_estimate; // accuracy estimate between extrapolations
	//		double first_increase;    // time value of the first increase in precision
	//		int cycle_num;    // cycle number used in extrapolations
	//		int success;      // success flag
	//		int multiplicity; // multiplicity
	//		int isReal;       // real flag:  0 - not real, 1 - real
	//		int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
	//		int isSing;       // singular flag: 0 - non-sigular, 1 - singular
	//	} post_process_t;
	
	
	
	//	typedef struct
	//	{
	//		int prec;
	//		point_data_d PD_d;
	//		point_data_mp PD_mp;
	//
	//		int last_approx_prec;       // precision of the last approximation
	//		point_d last_approx_d;      // last approximation to the end point
	//		point_mp last_approx_mp;    // last approximation to the end point
	//
	//		int retVal;
	//		int pathNum;
	//		int codim;
	//		double first_increase;
	//		double condition_number;
	//		double function_residual_d;
	//		mpf_t  function_residual_mp;
	//		double latest_newton_residual_d;
	//		mpf_t  latest_newton_residual_mp;
	//		double t_val_at_latest_sample_point_d;
	//		mpf_t  t_val_at_latest_sample_point_mp;
	//		double error_at_latest_sample_point_d;
	//		mpf_t  error_at_latest_sample_point_mp;
	//	} endgame_data_t;
	
	
	
	
}

//void findSingSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxCondNum, double finalTol, int regenToggle)

int BRfindSingularSolns(post_process_t *endPoints, int num_sols, int num_vars,
												tracker_config_t *T ){
	int ii, sing_count=0;
	
	for (ii = 0; ii < num_sols; ii++){
		if ( (endPoints[ii].cond_est >  T->cond_num_threshold) || (endPoints[ii].cond_est < 0.0) )
			endPoints[ii].isSing = 1;
		else
			endPoints[ii].isSing = 0;
		
		if (endPoints[ii].isSing)
		{
			sing_count++;
		}
	}
	
	return sing_count;
}


int BRfindFiniteSolns(post_process_t *endPoints, int num_sols, int num_vars,
											tracker_config_t *T ){
	int ii, jj, finite_count=0;
	
	//initialize temp stuffs
	comp_d dehom_coord_recip_d;
	comp_mp dehom_coord_recip_mp; init_mp(dehom_coord_recip_mp);
	vec_d dehom_d;   init_vec_d(dehom_d,num_vars-1);   dehom_d->size = num_vars-1;
	vec_mp dehom_mp; init_vec_mp(dehom_mp,num_vars-1); dehom_mp->size = num_vars-1;
	
	
	
	for (ii = 0; ii < num_sols; ii++){
		if (endPoints[ii].sol_prec<64) {
			set_d(dehom_coord_recip_d,endPoints[ii].sol_d[0]);
			recip_d(dehom_coord_recip_d,dehom_coord_recip_d);
			for (jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_d(&dehom_d->coord[jj],dehom_coord_recip_d,endPoints[ii].sol_d[jj])
			}
			
			if (infNormVec_d(dehom_d) < T->finiteThreshold){
				endPoints[ii].isFinite = 1;
				finite_count++;
			}
			else{
				endPoints[ii].isFinite = 0;
			}
			
		}
		else // high precision, do mp
		{
			change_prec_point_mp(dehom_mp,endPoints[ii].sol_prec);
			setprec_mp(dehom_coord_recip_mp,endPoints[ii].sol_prec);
			set_mp(dehom_coord_recip_mp,endPoints[ii].sol_mp[0]);
			for (jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_mp(&dehom_mp->coord[jj],dehom_coord_recip_mp,endPoints[ii].sol_mp[jj])
			}
			
			if (infNormVec_mp(dehom_mp) < T->finiteThreshold){
				endPoints[ii].isFinite = 1;
				finite_count++;
			}
			else{
				endPoints[ii].isFinite = 0;
			}
			
		}
	}
	
	clear_vec_d(dehom_d);
	clear_vec_mp(dehom_mp);
	clear_mp(dehom_coord_recip_mp);
	
	return finite_count;
}




int is_acceptable_solution(post_process_t endPoint, solver_configuration *solve_options){
	int indicator = 1;
	
	if ( (endPoint.multiplicity!=1) && (solve_options->allow_multiplicity==0) ) {
		indicator = 0;
	}
	
	if ( (endPoint.isSing==1) && (solve_options->allow_singular==0) ) {
		indicator = 0;
	}
	
	if ( (endPoint.isFinite!=1) && (solve_options->allow_infinite==0) ) {
		indicator = 0;
	}
	
	if ( (endPoint.success!=1) && (solve_options->allow_unsuccess==0) ) {
		indicator = 0;
	}
	
	
	return indicator;
}


//assumes that W has the number of variables already set, and the pts NOT allocated yet.  should be NULL
void BRpostProcessing(post_process_t *endPoints, witness_set *W_new, int num_pts,
											preproc_data *preProcData, tracker_config_t *T,
											solver_configuration *solve_options)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does the actual post processing for a zero dim run     *
 \***************************************************************/
{
	
	//	options->allow_multiplicity = 0;
	//	options->allow_singular = 0;
	//	options->allow_infinite = 0;
	//	options->allow_unsuccess = 0;
	
	int ii, jj;
	
	
	// sets the multiplicity and solution number in the endPoints data
	//direct from the bertini library:
	findMultSol(endPoints, num_pts, W_new->num_variables, preProcData, T->final_tol_times_mult);
	
	//sets the singularity flag in endPoints.
	//custom, derived from bertini's analagous call.
	int num_singular_solns = BRfindSingularSolns(endPoints, num_pts, W_new->num_variables, T);
	
	//sets the finite flag in endPoints.
	//custom, derived from bertini's analagous call.
	int num_finite_solns = BRfindFiniteSolns(endPoints, num_pts, W_new->num_variables, T);
	
	
	
	if (solve_options->show_status_summary==1) {
		printf("%d singular solutions\n",num_singular_solns);
		printf("%d finite solutions\n",num_finite_solns);
		
		for (ii=0; ii<num_pts; ++ii) {
			//		int success;      // success flag
			//		int multiplicity; // multiplicity
			//		int isReal;       // real flag:  0 - not real, 1 - real
			//		int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
			//		int isSing;       // singular flag: 0 - non-sigular, 1 - singular
			printf("solution %d, success %d, multi %d, isFinite %d, isSing %d\n",ii,endPoints[ii].success,endPoints[ii].multiplicity,endPoints[ii].isFinite,endPoints[ii].isSing);
		}
	}
	int num_actual_solns = 0;
	int *actual_solns_indices;
	actual_solns_indices = (int *)bmalloc(num_pts*sizeof(int));
	
	
	
	for (ii=0; ii<num_pts; ii++) {
		if ( is_acceptable_solution(endPoints[ii],solve_options) )//determine if acceptable based on current configuration
		{
			actual_solns_indices[num_actual_solns] = ii;
			num_actual_solns++;
		}
	}
	
	//initialize the structures for holding the produced data
	W_new->num_pts=num_actual_solns; W_new->num_pts=num_actual_solns;
  W_new->pts_d=(point_d *)bmalloc(num_actual_solns*sizeof(point_d));
  W_new->pts_mp=(point_mp *)bmalloc(num_actual_solns*sizeof(point_mp));
	
	
	
	for (ii=0; ii<num_actual_solns; ++ii) {
		
		init_vec_d(W_new->pts_d[ii],W_new->num_variables); init_vec_mp(W_new->pts_mp[ii],W_new->num_variables);
		W_new->pts_d[ii]->size = W_new->pts_mp[ii]->size = W_new->num_variables;
		
		if (endPoints[actual_solns_indices[ii]].sol_prec<64) {
			//copy out of the double structure.
			for (jj=0; jj<W_new->num_variables; jj++) {
				set_d(&W_new->pts_d[ii]->coord[jj],endPoints[actual_solns_indices[ii]].sol_d[jj]);
			}
			vec_d_to_mp(W_new->pts_mp[ii],W_new->pts_d[ii]);
		}
		else{
			//copy out of the mp structure.
			for (jj=0; jj<W_new->num_variables; jj++) {
				set_mp(&W_new->pts_mp[ii]->coord[jj],endPoints[actual_solns_indices[ii]].sol_mp[jj]);
			}
			vec_mp_to_d(W_new->pts_d[ii],W_new->pts_mp[ii]);
		}
	}
	
	free(actual_solns_indices);
	
	
	
	
  return;
}




