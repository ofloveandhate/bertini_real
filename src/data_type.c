#include "data_type.h"

void cp_patches_d(witness_set_d *W_out, witness_set_d W_in){
	int ii;
	
	
	
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
	
	return;
}

void init_witness_set_d(witness_set_d *W){
	W->num_variables = W->num_patches = W->num_linears = 0;
	W->patch = W->L = W->W.pts = NULL;
	
	return;
}


void dot_product_d(comp_d result, vec_d one, vec_d two){

	if (one->size!=two->size) {
		printf("attempting to dot two vectors not of the same size! (%d!=%d)\n",one->size,two->size);
		exit(-78);
	}
//	comp_d result;
	set_zero_d(result);
	
	
	comp_d temp;
	int ii;
	for (ii=0; ii<one->size; ++ii) {
		
		mul_d(temp,&one->coord[ii],&two->coord[ii]);
		add_d(result,result,temp);
	}

	
}

void write_dehomogenized_coordinates(witness_set_d W, char filename[]){
	
	FILE *OUT;
	OUT = safe_fopen_write(filename);

	int ii,jj;
	
//	printf("%%dehomogenizing %d points, writing to %s",W.W.num_pts,filename);
	
	fprintf(OUT,"%d\n\n",W.W.num_pts);
	
	for (ii=0; ii<W.W.num_pts; ++ii) {
		
		vec_d result;
		init_vec_d(result,0);
		dehomogenize(&result,W.W.pts[ii]);
		for (jj=0; jj<W.num_variables-1; jj++) {
			fprintf(OUT,"%.15le %.15le\n",result->coord[jj].r,result->coord[jj].i);
		}
		fprintf(OUT,"\n");
		
//		print_point_to_screen_matlab(W.W.pts[ii],"witnesspoint");
//		print_point_to_screen_matlab(result,"result");
		
		
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

void clear_witness_set(witness_set_d W){
	int ii;
	
	
	printf("freeing %d patches\n",W.num_patches);
	for (ii=0; ii<W.num_patches; ii++) {
		clear_vec_d(W.patch[ii]);
	}
	if (W.num_patches>0) {
		free(W.patch);
	}
	
	
	printf("freeing %d linears\n",W.num_linears);
	for (ii=0; ii<W.num_linears; ii++) {
		clear_vec_d(W.L[ii]);
	}
	if (W.num_linears>0) {
		free(W.L);
	}
	
	printf("freeing %d points\n",W.W.num_pts);
	for (ii=0; ii<W.W.num_pts; ii++) {
		clear_point_d(W.W.pts[ii]);
	}
	free(W.W.pts);
	
	
	
	
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
