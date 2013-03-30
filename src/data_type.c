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
		init_vec_d(W_out->L[counter],W_out->num_variables);
		vec_cp_d(W_out->L[counter],W_left.L[ii]);
		
		init_vec_mp(W_out->L_mp[counter],W_out->num_variables);
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
		init_vec_d(W_out->W.pts[counter],W_out->num_variables);
		vec_cp_d(W_out->W.pts[counter],W_left.W.pts[ii]);
		
		init_vec_mp(W_out->W_mp.pts[counter],W_out->num_variables);
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
	W->W.num_pts = W->W_mp.num_pts = 0;
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
	int ii;
	for (ii=0; ii<one->size; ++ii) {
		mul_mp(temp,&one->coord[ii],&two->coord[ii]);
		add_mp(result,result,temp);
	}
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
	printf("******\n%d points\n******\n",W.W.num_pts);
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
		printf(" %le+1i*%le;\n",M->coord[kk].r,M->coord[kk].i);
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
		mpf_out_str (NULL, 10, 7, M->coord[kk].r);
		printf("+1i*");
		mpf_out_str (NULL, 10, 7, M->coord[kk].i);
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
	mpf_out_str (NULL, 10, 6, M->r);
	printf(" ");
	mpf_out_str (NULL, 10, 6, M->i); // base 10, 6 digits
	printf("\n");
	return;
}


void print_path_retVal_message(int retVal){
	
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
	else if (retVal==100){
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
