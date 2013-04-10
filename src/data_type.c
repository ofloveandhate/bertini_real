#include "data_type.h"



void init_curveDecomp_d(curveDecomp_d *C){
	
	C->num_V0=C->num_V1=C->num_E=0;
	C->V0=C->V1=C->E=NULL;
	return;
}
void clear_curveDecomp_d(curveDecomp_d *C){
	
	int ii;
	for(ii=0;ii<C->num_V0;ii++)
		clear_vec_d(C->V0[ii].pt);
	free(C->V0);
	for(ii=0;ii<C->num_V1;ii++)
		clear_vec_d(C->V1[ii].pt);
	free(C->V1);
	for(ii=0;ii<C->num_E;ii++)
	{
		clear_vec_d(C->E[ii].midpt);
		clear_vec_d(C->E[ii].pi);
		clear_witness_set(C->E[ii].W);
	}
	free(C->E);
	
}
void clear_sample_d(sample_d *S){
	
	int ii;
	for(ii=0;ii<S->num_E;ii++)
	{
		clear_vec_d(S->pV[ii]);
		clear_mat_d(S->V[ii]);
	}
	free(S->pV);
	free(S->V);
	
}





void merge_witness_sets(witness_set_d *W_out,witness_set_d W_left,witness_set_d W_right){
	
	//error checking first
	if (W_left.num_variables != W_right.num_variables) {
		printf("merging two witness sets with differing numbers of variables.\n");
		exit(-341);
	}
	
	if (W_left.num_patches != W_right.num_patches) {
		printf("merging two witness sets with differing numbers of patch equations.\n");
		exit(-342);
	}
	
	if (W_left.patch_size != W_right.patch_size) {
		printf("merging two witness sets with differing sizes of patch(es).\n");
		exit(-343);
	}
	
	
//	if (W_left.num_linears == 0) { // should this be allowable?  i think yes...
//		printf("left witness_set has no linears.\n");
//		exit(-344);
//	}
//	
//	if (W_right.num_linears == 0) { // should this be allowable?  i think yes...
//		printf("right witness_set has no linears.\n");
//		exit(-345);
//	}
	
	int ii;
	//initialize the structure for holding the produced data
	W_out->num_variables = W_left.num_variables;
	
	
	W_out->num_linears = W_left.num_linears+W_right.num_linears;
	W_out->L = (vec_d *)bmalloc((W_left.num_linears+W_right.num_linears)*sizeof(vec_d));
	W_out->L_mp = (vec_mp *)bmalloc((W_left.num_linears + W_right.num_linears)*sizeof(vec_mp));
	// merge the left and right linears into the output.
	int counter = 0;
	for (ii=0; ii<W_left.num_linears; ii++) {
		init_vec_d(W_out->L[counter],W_out->num_variables); W_out->L[counter]->size = W_out->num_variables;
		vec_cp_d(W_out->L[counter],W_left.L[ii]);
		
		init_vec_mp(W_out->L_mp[counter],W_out->num_variables); W_out->L_mp[counter]->size = W_out->num_variables;
		vec_cp_mp(W_out->L_mp[counter],W_left.L_mp[ii]);
		counter++;
	}
	for (ii=0; ii<W_right.num_linears; ii++) {
		init_vec_d(W_out->L[counter],W_out->num_variables);
		vec_cp_d(W_out->L[counter],W_right.L[ii]);
		
		init_vec_mp(W_out->L_mp[counter],W_out->num_variables);
		vec_cp_mp(W_out->L_mp[counter],W_right.L_mp[ii]);
		counter++;
	}
	
	
	
	//set the number of points
	W_out->W.num_pts  = W_left.W.num_pts+W_right.W.num_pts;
	W_out->W_mp.num_pts = W_left.W.num_pts+W_right.W.num_pts;
	

  W_out->W.pts=(point_d *)bmalloc(W_out->W.num_pts*sizeof(point_d));
	W_out->W_mp.pts=(point_mp *)bmalloc(W_out->W.num_pts*sizeof(point_mp));
	
	counter = 0;
	for (ii=0; ii<W_left.W.num_pts; ii++) {
		init_vec_d(W_out->W.pts[counter],W_out->num_variables); W_out->W.pts[counter]->size = W_out->num_variables;
		vec_cp_d(W_out->W.pts[counter],W_left.W.pts[ii]);
		
		init_vec_mp(W_out->W_mp.pts[counter],W_out->num_variables); W_out->W_mp.pts[counter]->size = W_out->num_variables;
		vec_cp_mp(  W_out->W_mp.pts[counter],W_left.W_mp.pts[ii]);
		counter++;
	}
	
	for (ii=0; ii<W_right.W.num_pts; ii++) {
		init_vec_d(W_out->W.pts[counter],W_out->num_variables);
		vec_cp_d(W_out->W.pts[counter],W_right.W.pts[ii]);
		
		init_vec_mp(W_out->W_mp.pts[counter],W_out->num_variables);
		vec_cp_mp(  W_out->W_mp.pts[counter],W_right.W_mp.pts[ii]);
		counter++;
	}

	cp_patches(W_out,W_left); // copy the patches over from the original witness set
	
	return;
}//re: merge_witness_sets


void init_variable_names(witness_set_d *W, int num_vars){
	int ii;
	
	if (W->variable_names!=NULL) {
		printf("attempting to initialize non-null variable names\n");
		exit(-1);
	}
	
	W->variable_names = (char **)bmalloc(num_vars*sizeof(char*));
	for (ii=0; ii<num_vars; ++ii) {
		W->variable_names[ii] = (char*) bmalloc(64*sizeof(char));
	}
	
	return;
}

//copies names from old to new
void cp_names(witness_set_d *W_out, witness_set_d W_in){
	int ii;
	
	if (W_in.num_variables==0 || W_in.variable_names==NULL) {
		printf("attempting to copy variable names from witness_set with no variables");
		exit(1333);
	}
	
	if (W_out->variable_names==NULL) {
		init_variable_names(W_out,W_in.num_variables);
	}
	
	for (ii=0; ii<W_in.num_variables; ++ii) {
		strcpy(W_out->variable_names[ii],  W_in.variable_names[ii]);
	}
	
//	printf("copied these variable names:\n");
//	for (ii=0; ii<W_in.num_variables; ++ii) {
//		printf("%s\n",W_out->variable_names[ii]);
//	}
//	mypause();
}




//copies the mp and d linears from in to out.
void cp_linears(witness_set_d *W_out, witness_set_d W_in){
	int ii;
	
	
	W_out->num_linears = W_in.num_linears;
	printf("%d linears to copy.\n",W_in.num_linears);


	//
	// DOUBLE COPIES
	//
	if (W_out->L==NULL) {
		W_out->L = (vec_d *)bmalloc(W_in.num_linears * sizeof(vec_d));
	}
	else
	{
		W_out->L = (vec_d *)brealloc(W_out->L, W_in.num_linears * sizeof(vec_d));
	}
	
	for (ii=0; ii<W_in.num_linears; ++ii) {
		init_vec_d(W_out->L[ii],W_in.L[ii]->size);
		vec_cp_d(W_out->L[ii],W_in.L[ii]);
		W_out->L[ii]->size = W_in.L[ii]->size;
	}
	
	//
	// MP COPIES
	//
	
	if (W_out->L_mp==NULL) {
		W_out->L_mp = (vec_mp *)bmalloc(W_in.num_linears * sizeof(vec_mp));
	}
	else
	{
		W_out->L_mp = (vec_mp *)brealloc(W_out->L_mp, W_in.num_linears * sizeof(vec_mp));
	}
	
	for (ii=0; ii<W_in.num_linears; ++ii) {
		init_vec_mp(W_out->L_mp[ii],W_in.L_mp[ii]->size);
		vec_cp_mp(W_out->L_mp[ii],W_in.L_mp[ii]);
		W_out->L_mp[ii]->size = W_in.L_mp[ii]->size;
	}
	
	
	return;
}


void cp_patches(witness_set_d *W_out, witness_set_d W_in){
	int ii;
	
	
	W_out->patch_size = W_in.patch_size;
	W_out->num_patches = W_in.num_patches;
	
	
	if (W_out->patch==NULL) {
		W_out->patch = (vec_d *)bmalloc(W_in.num_patches * sizeof(vec_d));
	}
	else
	{
		W_out->patch = (vec_d *)brealloc(W_out->patch, W_in.num_patches * sizeof(vec_d));
	}
	
	for (ii=0; ii<W_in.num_patches; ++ii) {
		init_vec_d(W_out->patch[ii],0);
		change_size_vec_d(W_out->patch[ii],W_in.patch[ii]->size);
		vec_cp_d(W_out->patch[ii],W_in.patch[ii]);
		W_out->patch[ii]->size = W_in.patch[ii]->size;
	}
	
	
	if (W_out->patch_mp==NULL) {
		W_out->patch_mp = (vec_mp *)bmalloc(W_in.num_patches * sizeof(vec_mp));
	}
	else
	{
		W_out->patch_mp = (vec_mp *)brealloc(W_out->patch_mp, W_in.num_patches * sizeof(vec_mp));
	}
	
	for (ii=0; ii<W_in.num_patches; ++ii) {
		init_vec_mp(W_out->patch_mp[ii],0);
		change_size_vec_mp(W_out->patch_mp[ii],W_in.patch_mp[ii]->size);
		vec_cp_mp(W_out->patch_mp[ii],W_in.patch_mp[ii]);
		W_out->patch_mp[ii]->size = W_in.patch_mp[ii]->size;
	}
	
	
	return;
}


// initializes witness set, both the mp data and double.
//only call this on new witness sets, otherwise you will leak memory when you set the pointers to NULL.
void init_witness_set_d(witness_set_d *W){
	W->num_variables = W->num_patches = W->num_linears = 0;
	W->patch = W->L = W->W.pts = NULL;
	W->patch_mp = W->L_mp = W->W_mp.pts = NULL;
	W->variable_names = NULL;
	W->W.num_pts = W->W_mp.num_pts = 0;
	W->incidence_number = -1;
	return;
}


void dot_product_d(comp_d result, vec_d one, vec_d two){
	if (one->size!=two->size) {
		printf("attempting to dot two vectors not of the same size! (%d!=%d)\n",one->size,two->size);
		exit(-78);
	}
	set_zero_d(result);
	
	comp_d temp;
	int ii;
	for (ii=0; ii<one->size; ++ii) {
		mul_d(temp,&one->coord[ii],&two->coord[ii]);
		add_d(result,result,temp);
	}
}
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two){
	if (one->size!=two->size) {
		printf("attempting to dot two vectors not of the same size! (%d!=%d)\n",one->size,two->size);
		exit(-78);
	}
	set_zero_mp(result);
	
	comp_mp temp; init_mp(temp);
	comp_mp temp2; init_mp(temp2);
	int ii;
	for (ii=0; ii<one->size; ++ii) {
		set_mp(temp2,result)
		mul_mp(temp,&one->coord[ii],&two->coord[ii]);
		add_mp(result,temp2,temp);
	}
	clear_mp(temp);clear_mp(temp2);
}


void write_dehomogenized_coordinates(witness_set_d W, char filename[]){
	
	FILE *OUT;
	OUT = safe_fopen_write(filename);

	int ii,jj;
	
	
	fprintf(OUT,"%d\n\n",W.W.num_pts);
	
	for (ii=0; ii<W.W.num_pts; ++ii) {
		if (W.MPType==1){
			
			vec_mp result;
			init_vec_mp(result,1);
			dehomogenize_mp(&result,W.W_mp.pts[ii]);
			for (jj=0; jj<W.num_variables-1; jj++) {
				print_mp(OUT, 0, &result->coord[jj]);
				fprintf(OUT, "\n");
			}
		}
		else{
			vec_d result;
			init_vec_d(result,0);
			dehomogenize(&result,W.W.pts[ii]);
			for (jj=0; jj<W.num_variables-1; jj++) {
				fprintf(OUT,"%.15le %.15le\n",result->coord[jj].r,result->coord[jj].i);
			}

		}
		
		
		fprintf(OUT,"\n");

	}
	
	fclose(OUT);
}


void dehomogenize(vec_d *result, vec_d dehom_me){
	comp_d denom;
	change_size_vec_d(*result,dehom_me->size-1);
	(*result)->size = dehom_me->size-1;
	set_d(denom, &dehom_me->coord[0]);
	
	int ii;
	for (ii=0; ii<dehom_me->size-1; ++ii) {
		set_d( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
		div_d(&(*result)->coord[ii],&(*result)->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	return;
}
void dehomogenize_mp(vec_mp *result, vec_mp dehom_me){
	comp_mp denom;
	init_mp(denom);
	change_size_vec_mp((*result),dehom_me->size-1);

	(*result)->size = dehom_me->size-1;
	
//	mpf_out_str (NULL, 10, 10, dehom_me->coord[0].r);
	
//	gmp
//	gmp_printf ("dehom_me[0].r  %Zd\n", dehom_me->coord[0].r);
//	
	
	set_mp(denom, &dehom_me->coord[0]);
	int ii;
	for (ii=0; ii<dehom_me->size-1; ++ii) {
		set_mp( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
		div_mp(&(*result)->coord[ii],&(*result)->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	return;
}


void clear_witness_set(witness_set_d W){
	int ii;
	
	
//	printf("freeing %d patches\n",W.num_patches);
	for (ii=0; ii<W.num_patches; ii++) {
		clear_vec_d(W.patch[ii]);
	}
	if (W.num_patches>0) {
		free(W.patch);
	}
	
	
//	printf("freeing %d linears\n",W.num_linears);
	for (ii=0; ii<W.num_linears; ii++) {
		clear_vec_d(W.L[ii]);
	}
	if (W.num_linears>0) {
		free(W.L);
	}
	
//	printf("freeing %d points\n",W.W.num_pts);
	for (ii=0; ii<W.W.num_pts; ii++) {
		clear_point_d(W.W.pts[ii]);
	}
	free(W.W.pts);
	
	
	//now clear the mp structures.
	
//	printf("freeing %d patches\n",W.num_patches);
	for (ii=0; ii<W.num_patches; ii++) {
		clear_vec_mp(W.patch_mp[ii]);
	}
	if (W.num_patches>0) {
		free(W.patch_mp);
	}
	
	
//	printf("freeing %d linears\n",W.num_linears);
	for (ii=0; ii<W.num_linears; ii++) {
		clear_vec_d(W.L_mp[ii]);
	}
	if (W.num_linears>0) {
		free(W.L_mp);
	}
	
//	printf("freeing %d points\n",W.W_mp.num_pts);
	for (ii=0; ii<W.W_mp.num_pts; ii++) {
		clear_point_mp(W.W_mp.pts[ii]);
	}
	free(W.W_mp.pts);
	
	
	
	return;
}

void print_witness_set_to_screen(witness_set_d W){
	int ii;
	printf("******\n%d points in double, %d points in mp\n******\n",W.W.num_pts,W.W_mp.num_pts);
	for (ii=0; ii<W.W.num_pts; ii++) {
		printf("the%dth",ii);
		print_point_to_screen_matlab(W.W.pts[ii],"point");
	}
		
	printf("******\n%d linears\n******\n",W.num_linears);
	
	for (ii=0; ii<W.num_linears; ii++) {
		printf("the%dth",ii);
		print_point_to_screen_matlab(W.L[ii],"linear");
	}
	
	printf("\n\n");
}


void print_point_to_screen_matlab(vec_d M, char name[])
{
	int kk;
	
	if (M->size==0) {
		printf("requested to print a vector '%s' which had size==0.  exiting\n",name);
		exit(-1);
	}
	
	printf("%s = [...\n",name);
	for (kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		printf(" %.10le+1i*%.10le;\n",M->coord[kk].r,M->coord[kk].i);
	}
	printf("];\n\n");
}
void print_point_to_screen_matlab_mp(vec_mp M, char name[])
{
	int kk;
	
	if (M->size==0) {
		printf("requested to print a vector '%s' which had size==0.  exiting\n",name);
		exit(-1);
	}
	
	printf("%s = [...\n",name);
	for (kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		mpf_out_str (NULL, 10, 16, M->coord[kk].r);
		printf("+1i*");
		mpf_out_str (NULL, 10, 16, M->coord[kk].i);
		printf(";\n");
//		printf(" %le+1i*%le;\n",M->coord[kk].r,M->coord[kk].i);
	}
	printf("];\n\n");
}



void print_matrix_to_screen_matlab(mat_d M, char name[])
{
	int jj,kk;
	
	printf("%%matrix '%s' has dimensions %dx%d\n",name, M->rows,M->cols);
	printf("%s = [...\n",name);
	for (kk = 0; kk < M->rows; kk++)
	{ // print kth row
		for (jj = 0; jj < M->cols; jj++)
		{
			printf(" %.15le+1i*%.15le",M->entry[kk][jj].r,M->entry[kk][jj].i);
		}
		printf(";\n");
	}
	printf("];\n\n");
}
void print_matrix_to_screen_matlab_mp(mat_mp M, char name[])
{
	int jj,kk;
	
	printf("%%matrix '%s' has dimensions %dx%d\n",name, M->rows,M->cols);
	printf("%s = [...\n",name);
	for (kk = 0; kk < M->rows; kk++)
	{ // print kth row
		for (jj = 0; jj < M->cols; jj++)
		{
			
			mpf_out_str (NULL, 10, 7, M->entry[kk][jj].r);
			printf("+1i*");
			mpf_out_str (NULL, 10, 7, M->entry[kk][jj].i); // base 10 , 7 digits
			printf(" ");
			
//			printf(" %.15le+1i*%.15le",M->entry[kk][jj].r,M->entry[kk][jj].i);
		}
		printf(";\n");
	}
	printf("];\n\n");
}


void print_comp_mp_matlab(comp_mp M,char name[]){
	printf("%s=",name);
	mpf_out_str (NULL, 10, 6, M->r);
	printf("+1i*");
	mpf_out_str (NULL, 10, 6, M->i); // base 10, 6 digits
	printf("\n");
	return;
}


void print_path_retVal_message(int retVal){
	
	
	if (retVal==100) {
		printf("max_prec_reached\n");
	}
	if (retVal==-50) {
		printf("reached_minTrackT\n");
	}
	else if (retVal==-200){
		printf("cycle_num_too_high\n");
	}
	else if (retVal==-15){
		printf("PSEG_failed\n");
	}
	else if (retVal==-2){
		printf("going_to_infinity\n");
	}
	else if (retVal==-4){
		printf("security_max\n");
	}
	else if (retVal==-3){
		printf("step_size_too_small\n");
	}
	else if (retVal==-10){
		printf("too_many_steps\n");
	}
	else if (retVal==-20){
		printf("refining_failed\n");
	}
	else if (retVal==-100){
		printf("higher_prec_needed\n");
	}
	else if (retVal==-99){
		printf("retVal_NAN\n");
	}
	else if (retVal==-98){
		printf("retVal_Bertini_Junk\n");
	}
	else if (retVal==-97){
		printf("Failed_to_converge (used in newton iterations)\n");
	}
	else if (retVal==-22){
		printf("sharpening_singular_endpoint\nthis is used when the sharpening sees that the endpoint is singular and cannot sharpen it\n");
	}
	else if (retVal==-21){
		printf("sharpening_failed\nthis is used when the sharpening of an endpoint does not reach the desired tolerance\n");
	}
	else if (retVal==-22){
		printf("higher_dim\nthis is used in regeneration when an endpoint lies on a higher dimensional component\n");
	}

	
	return;
}



void endgamedate_to_endpoint(post_process_t *endPoint, endgame_data_t *EG){
	
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
		printf("endpoint %d has multiplicity %d\n",ii,endPoints[ii].multiplicity);
    if (endPoints[ii].multiplicity > 0)
			
    { // regeneration endpoints are always non-singular
      if ( (endPoints[ii].cond_est >  T->cond_num_threshold) || (endPoints[ii].cond_est < 0.0) || (endPoints[ii].multiplicity > 1))
        endPoints[ii].isSing = 1;
      else
        endPoints[ii].isSing = 0;
			
      if (endPoints[ii].isSing)
      { // print info to SINGU
				
        // increment the number of singular solutions printed
        sing_count++;
			}
		}
	}
	
	return sing_count;
}

//void BRPostProcess(FILE *OUT, post_process_t *endPoints, int num_sols, int num_vars, double final_tol, tracker_config_t *T, preproc_data *PPD, int num_crossings, int convergence_failures, char *inputName, int regenToggle, int eqbyeqToggle)
///***************************************************************\
// * USAGE:                                                        *
// * ARGUMENTS:                                                    *
// * RETURN VALUES:                                                *
// * NOTES: does the actual post processing for a zero dim run     *
// \***************************************************************/
//{
//  int i, num_good_sols = 0, size = 0;
//  FILE *NAMES = NULL;
//  char ch, **name_table = (char **)bmalloc(num_vars * sizeof(char *));
//  int *origErrorIsInf = (int *)bmalloc(num_sols * sizeof(int));
//  double *origErrorEst = (double *)bmalloc(num_sols * sizeof(double));
//  point_d *dehomPoints_d = T->MPType == 1 ? NULL : (point_d *)bmalloc(num_sols * sizeof(point_d));
//  point_mp *dehomPoints_mp = T->MPType == 0 ? NULL : (point_mp *)bmalloc(num_sols * sizeof(point_mp));
//	
//  // Reading in the names of the variables.
//  NAMES = fopen("names.out", "r");
//  if (NAMES == NULL)
//  {
//    printf("ERROR: 'names.out' does not exist!\n");
//    bexit(ERROR_FILE_NOT_EXIST);
//  }
//  for (i = 0; i < num_vars; i++)
//  { // initial allocation
//    size = 1;
//    name_table[i] = (char *)bmalloc(size * sizeof(char));
//    // read in name
//    while ((name_table[i][size - 1] = fgetc(NAMES)) != '\n')
//    {
//      size++;
//      name_table[i] = (char *)brealloc(name_table[i], size * sizeof(char));
//    }
//    name_table[i][size - 1] = '\0';
//  }
//  fclose(NAMES);
//	
//  // setup the dehomogenized points and error estimates
//  for (i = 0; i < num_sols; i++)
//  { // setup dehom
//    if (endPoints[i].sol_prec < 64)
//    { // setup _d
//      init_point_d(dehomPoints_d[i], 0);
//      getDehomPoint_comp_d(dehomPoints_d[i], &origErrorIsInf[i], &origErrorEst[i], endPoints[i].sol_d, num_vars, PPD, endPoints[i].accuracy_estimate);
//    }
//    else
//    { // setup _mp
//      initMP(endPoints[i].sol_prec);
//      init_point_mp(dehomPoints_mp[i], 0);
//      getDehomPoint_comp_mp(dehomPoints_mp[i], &origErrorIsInf[i], &origErrorEst[i], endPoints[i].sol_mp, num_vars, PPD, endPoints[i].accuracy_estimate);
//    }
//  }
//	
//  // print top of main_data
//  printFileHeader(OUT, num_vars, PPD, name_table);
//	
//  // create various output from this data
//  zeroDimOutputChart(endPoints, dehomPoints_d, dehomPoints_mp, name_table, num_sols, num_vars, PPD, T->finiteThreshold, T->real_threshold, final_tol, T->cond_num_threshold, regenToggle, eqbyeqToggle);
//  createRawSoln(endPoints, dehomPoints_d, dehomPoints_mp, num_sols, num_vars, PPD);
//	
//  // Now we sort the points - and write the body of main_data
//  if (regenToggle)
//    fprintf(OUT, "\nNOTE: Since regeneration is being used, only non-singular solutions are printed.\n\n");
//  else if (eqbyeqToggle)
//    fprintf(OUT, "\nNOTE: Since equation-by-equation is being used, only non-singular solutions are printed.\n\n");
//	
//  for (i = 0; i < num_sols; i++)
//    if (endPoints[i].multiplicity > 0 && (!(regenToggle || eqbyeqToggle) || ((regenToggle || eqbyeqToggle) && endPoints[i].isSing == 0)))  // print only for regen/eqbyeq if non-singular
//    {
//      if (endPoints[i].sol_prec < 64)
//      { // print header for the solution
//        printMainDataPointHeader_d(OUT, endPoints[i].sol_num, endPoints[i].path_num, endPoints[i].cond_est, endPoints[i].function_resid_d, endPoints[i].newton_resid_d, endPoints[i].final_t, endPoints[i].accuracy_estimate, endPoints[i].sol_prec, endPoints[i].first_increase, endPoints[i].cycle_num, endPoints[i].success, origErrorIsInf[i], origErrorEst[i]);
//        // print the point
//        printRefinedPoint_d(OUT, dehomPoints_d[i]);
//      }
//      else
//      { // print header for the solution
//        printMainDataPointHeader_mp(OUT, endPoints[i].sol_num, endPoints[i].path_num, endPoints[i].cond_est, endPoints[i].function_resid_mp, endPoints[i].newton_resid_mp, endPoints[i].final_t, endPoints[i].accuracy_estimate, endPoints[i].sol_prec, endPoints[i].first_increase, endPoints[i].cycle_num, endPoints[i].success, origErrorIsInf[i], origErrorEst[i]);
//        // print the point
//        printRefinedPoint_mp(OUT, dehomPoints_mp[i]);
//      }
//			
//      // print footer for the solution
//      printMainDataPointFooter(OUT, i, num_sols, endPoints);
//			
//      num_good_sols++;
//    }
//	
//  fprintf(OUT, "-------------------------\n");
//	
//  if (num_good_sols == 1)
//    fprintf(OUT, "At tol=%.12e, there appears to be a unique solution.\n", final_tol);
//  else
//    fprintf(OUT, "At tol=%.12e, there appear to be %d solutions.\n", final_tol, num_good_sols);
//	
//  if (num_crossings > 0)
//  {
//    printf("\nIt appears that %d path crossing(s) occurred prior to t=tEndgame\n\n", num_crossings);
//    fprintf(OUT, "\nIt appears that %d path crossing(s) occurred prior to t=tEndgame\n\n", num_crossings);
//    printf("   Try adjusting TRACKTOLBEFOREEG in the input file\n\n");
//  }
//	
//  if (convergence_failures == 1)
//  {
//    fprintf(OUT, "\nThere is 1 path that was tracked near the target time which did not converge. Try adjusting FinalTol or NbhdRadius in the input file.\n");
//    fprintf(OUT, "This path is marked with '#' after its path number.\n");
//  }
//  else if (convergence_failures > 1)
//  {
//    fprintf(OUT, "\nThere are %d paths that were tracked near the target time which did not converge. Try adjusting FinalTol or NbhdRadius in the input file.\n", convergence_failures);
//    fprintf(OUT, "These paths are marked with '#' after the path numbers.\n");
//  }
//	
//  // print the input onto the bottom of refined solutions
//  NAMES = fopen(inputName, "r");
//  if (NAMES == NULL)
//  {
//    printf("ERROR: '%s' does not exist!\n", inputName);
//    bexit(ERROR_FILE_NOT_EXIST);
//  }
//	
//  while((ch = fgetc(NAMES)) != EOF)
//    fputc(ch, OUT);
//  fclose(NAMES);
//  remove(inputName);
//	
//  // clear memory
//  for (i = num_vars - 1; i >= 0; i--)
//    free(name_table[i]);
//  free(name_table);
//  for (i = 0; i < num_sols; i++)
//  {
//    if (endPoints[i].sol_prec < 64)
//    {
//      clear_point_d(dehomPoints_d[i]);
//    }
//    else
//    {
//      clear_point_mp(dehomPoints_mp[i]);
//    }
//  }
//  if (T->MPType != 1)
//    free(dehomPoints_d);
//  if (T->MPType != 0)
//    free(dehomPoints_mp);
//  free(origErrorEst);
//  free(origErrorIsInf);
//	
//  return;
//}














