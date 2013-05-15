#include "witnessSet.h"




//use:  call to parse the file witness_set_file, into the struct W.  fills the values W, and L, patch, as well as their _mp counterparts
	


int witnessSetParse(witness_set *W, char *witness_set_file, const int num_vars){
  

  int ii, jj;
	
  FILE *IN = safe_fopen_read(witness_set_file);
  
  
	int codim, comp_num;
	int num_patches, patch_size, num_linears, num_pts, num_vars_in_linears;

	
	
  fscanf(IN, "%d %d %d", &num_pts, &codim, &comp_num); scanRestOfLine(IN);
	W->codim = codim;
	W->comp_num = comp_num;
  W->num_pts = W->num_pts = num_pts;
  W->pts_d=(point_d *)bmalloc(W->num_pts*sizeof(point_d));
	W->pts_mp=(point_mp *)bmalloc(W->num_pts*sizeof(point_mp));

	
  W->num_variables = num_vars;

  for (ii=0; ii < num_pts; ii++) {
    init_point_mp2(W->pts_mp[ii],num_vars,1024); init_point_d(W->pts_d[ii],num_vars);
		W->pts_mp[ii]->size = W->pts_d[ii]->size = num_vars;
    //read the witness points into memory
    for (jj=0; jj < num_vars; jj=jj+1) {
			mpf_inp_str(W->pts_mp[ii]->coord[jj].r, IN, 10); // 10 is the base
			mpf_inp_str(W->pts_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
		vec_mp_to_d(W->pts_d[ii],W->pts_mp[ii]);
  }
  
	

	fscanf(IN, "%d %d", &num_linears, &num_vars_in_linears);  scanRestOfLine(IN);
//	printf("%d %d\n",num_linears, num_vars_in_linears);
	W->num_linears = num_linears;
	W->L = (vec_d *)bmalloc(num_linears*sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(num_linears*sizeof(vec_mp));
	
  for (ii=0; ii < num_linears; ii++) {
		init_vec_mp2(W->L_mp[ii],num_vars_in_linears,1024); init_vec_d(W->L[ii],num_vars_in_linears);
		
		W->L[ii]->size =  W->L_mp[ii]->size = num_vars_in_linears;
		
    //read the witness linears into memory
    for (jj=0; jj < num_vars_in_linears; jj++) {
			mpf_inp_str(W->L_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->L_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
		vec_mp_to_d(W->L[ii],W->L_mp[ii]);
  }
  
	
	fscanf(IN, "%d %d", &num_patches, &patch_size); scanRestOfLine(IN);
	W->patch_size = patch_size;
	W->num_patches = num_patches;
	
	W->patch = (vec_d *)bmalloc(num_patches*sizeof(vec_d));
	W->patch_mp = (vec_mp *)bmalloc(num_patches*sizeof(vec_mp));
	
  for (ii=0; ii < num_patches; ii++) {
		init_vec_mp2(W->patch_mp[ii],patch_size,1024);//default max_prec is 1024
		init_vec_d(W->patch[ii],patch_size);
		
		W->patch_mp[ii]->size = W->patch[ii]->size = patch_size;
		
    //read the patch into memory
    for (jj=0; jj < patch_size; jj++) {
			mpf_inp_str(W->patch_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->patch_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
		vec_mp_to_d(W->patch[ii],W->patch_mp[ii]);
  }
	
  fclose(IN);
	
	


  return 0;
}





void merge_witness_sets(witness_set *W_out,witness_set W_left,witness_set W_right){
	
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
	
	cp_names(W_out, W_right);
	
	
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
	W_out->num_pts  = W_left.num_pts+W_right.num_pts;
	
	
  W_out->pts_d =  (point_d *)bmalloc(W_out->num_pts*sizeof(point_d));
	W_out->pts_mp = (point_mp *)bmalloc(W_out->num_pts*sizeof(point_mp));
	
	counter = 0;
	for (ii=0; ii<W_left.num_pts; ii++) {
		init_vec_d(W_out->pts_d[counter],W_out->num_variables); W_out->pts_d[counter]->size = W_out->num_variables;
		vec_cp_d(W_out->pts_d[counter],W_left.pts_d[ii]);
		
		init_vec_mp(W_out->pts_mp[counter],W_out->num_variables); W_out->pts_mp[counter]->size = W_out->num_variables;
		vec_cp_mp(  W_out->pts_mp[counter],W_left.pts_mp[ii]);
		counter++;
	}
	
	for (ii=0; ii<W_right.num_pts; ii++) {
		init_vec_d(W_out->pts_d[counter],W_out->num_variables);
		vec_cp_d(W_out->pts_d[counter],W_right.pts_d[ii]);
		
		init_vec_mp(W_out->pts_mp[counter],W_out->num_variables);
		vec_cp_mp(  W_out->pts_mp[counter],W_right.pts_mp[ii]);
		counter++;
	}
	
	cp_patches(W_out,W_left); // copy the patches over from the original witness set
	
	return;
}//re: merge_witness_sets


void init_variable_names(witness_set *W, int num_vars){
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
void cp_names(witness_set *W_out, witness_set W_in){
	int ii;
	
	if (W_in.num_variables==0) {
		printf("\nattempting to copy variable names from witness_set with no variables\n");
		exit(1333);
	}
	
	if (W_in.variable_names==NULL) {
		printf("\nattempting to copy variable names from witness_set unset variable names\n");
		exit(1334);
	}
	
	if (W_out->variable_names==NULL) {
		init_variable_names(W_out,W_in.num_variables);
	}
	
	for (ii=0; ii<W_in.num_variables; ++ii) {
		strcpy(W_out->variable_names[ii],  W_in.variable_names[ii]);
	}
}




//copies the mp and d linears from in to out.
void cp_linears(witness_set *W_out, witness_set W_in){
	int ii;
	
	
	W_out->num_linears = W_in.num_linears;
	
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


void cp_patches(witness_set *W_out, witness_set W_in){
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

void cp_witness_set(witness_set *W_out, witness_set W_in){
	int ii;
	
	W_out->codim = W_in.codim;
	W_out->comp_num = W_in.comp_num;
	W_out->incidence_number = W_in.incidence_number;
	W_out->num_variables = W_in.num_variables;
	W_out->num_var_gps = W_in.num_var_gps;
	W_out->num_linears = W_in.num_linears;
	W_out->num_patches = W_in.num_patches;
	W_out->patch_size = W_in.patch_size;
	W_out->MPType = W_in.MPType;
	
	cp_patches(W_out,W_in);
	cp_linears(W_out,W_in);
	cp_names(W_out,W_in);
	
	//set the number of points
	W_out->num_pts  = W_in.num_pts;
	
	
  W_out->pts_d=(point_d *)bmalloc(W_out->num_pts*sizeof(point_d));
	W_out->pts_mp=(point_mp *)bmalloc(W_out->num_pts*sizeof(point_mp));
	
	for (ii=0; ii<W_in.num_pts; ii++) {
		init_vec_d(W_out->pts_d[ii],W_out->num_variables); W_out->pts_d[ii]->size = W_out->num_variables;
		vec_cp_d(W_out->pts_d[ii],W_in.pts_d[ii]);
		
		init_vec_mp(W_out->pts_mp[ii],W_out->num_variables); W_out->pts_mp[ii]->size = W_out->num_variables;
		vec_cp_mp(  W_out->pts_mp[ii],W_in.pts_mp[ii]);
	}
	
	return;
	
}

// initializes witness set, both the mp data and double.
//only call this on new witness sets, otherwise you will leak memory when you set the pointers to NULL.
void init_witness_set_d(witness_set *W){
	W->num_variables = W->num_patches = W->num_linears = 0;
	W->patch = W->L = W->pts_d = NULL;
	W->patch_mp = W->L_mp = W->pts_mp = NULL;
	W->variable_names = NULL;
	W->num_pts = W->num_pts = 0;
	W->incidence_number = -1;
	return;
}



//requires the witness set to have set the number of variables.
void get_variable_names(witness_set *W){
	int ii;
	
	
	FILE *IN = safe_fopen_read("names.out");
	
	init_variable_names(W,W->num_variables);
	for (ii=0; ii<W->num_variables; ++ii){
		fscanf(IN,"%s\n",W->variable_names[ii]);
	}
	
	fclose(IN);
	
	
}





void write_homogeneous_coordinates(witness_set W, char filename[])
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename); // open the output file
	
	
	fprintf(OUT,"%d\n\n",W.num_pts); // print the header line
	
	for (ii=0; ii<W.num_pts; ++ii) {
		for (jj=0; jj<W.num_variables; jj++) {
			fprintf(OUT,"%.15le %.15le\n",W.pts_d[ii]->coord[jj].r,W.pts_d[ii]->coord[jj].i);
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void write_dehomogenized_coordinates(witness_set W, char filename[]){
	int ii,jj;
	
	FILE *OUT = safe_fopen_write(filename); // open the output file.
	
	fprintf(OUT,"%d\n\n",W.num_pts); // print the header line
	for (ii=0; ii<W.num_pts; ++ii) {
		if (W.MPType==1){ // both fields should be populated anyway?
			vec_mp result;
			init_vec_mp(result,1);
			dehomogenize_mp(&result,W.pts_mp[ii]);
			for (jj=0; jj<W.num_variables-1; jj++) {
				print_mp(OUT, 0, &result->coord[jj]);
				fprintf(OUT, "\n");
			}
		}
		else{
			vec_d result;
			init_vec_d(result,0);
			dehomogenize_d(&result,W.pts_d[ii]);
			for (jj=0; jj<W.num_variables-1; jj++) {
				fprintf(OUT,"%.15le %.15le\n",result->coord[jj].r,result->coord[jj].i);
			}
			
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
	return;
}




void write_linears(witness_set W, char filename[])
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename); // open the output file
	
	
	fprintf(OUT,"%d\n\n",W.num_linears); // print the header line
	
	for (ii=0; ii<W.num_linears; ++ii) {
		for (jj=0; jj<W.num_variables; jj++) {
			fprintf(OUT,"%.15le %.15le\n",W.L[ii]->coord[jj].r,W.L[ii]->coord[jj].i);
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}



void clear_witness_set(witness_set W){
	int ii;
	
	
	for (ii=0; ii<W.num_patches; ii++) {
		clear_vec_d(W.patch[ii]);
		clear_vec_mp(W.patch_mp[ii]);
	}
	if (W.num_patches>0) {
		free(W.patch);
		free(W.patch_mp);
	}
	
	
	for (ii=0; ii<W.num_linears; ii++) {
		clear_vec_d(W.L[ii]);
		clear_vec_mp(W.L_mp[ii]);
	}
	if (W.num_linears>0) {
		free(W.L);
		free(W.L_mp);
	}
	
	for (ii=0; ii<W.num_pts; ii++) {
		clear_point_d(W.pts_d[ii]);
		clear_point_mp(W.pts_mp[ii]);
	}
	if (W.num_pts>0) {
		free(W.pts_d);
		free(W.pts_mp);
	}
	
	
	if (W.variable_names!=NULL) {
		for (ii=0; ii<W.num_variables; ++ii) {
			free(W.variable_names[ii]);
		}
		free(W.variable_names);
		W.variable_names=NULL;
	}
	
	return;
}

void print_witness_set_to_screen(witness_set W){
	int ii;
	printf("******\n%d points in double, %d points in mp\n******\n",W.num_pts,W.num_pts);
	for (ii=0; ii<W.num_pts; ii++) {
		printf("the%dth",ii);
		print_point_to_screen_matlab(W.pts_d[ii],"point");
	}
	
	printf("******\n%d linears\n******\n",W.num_linears);
	
	for (ii=0; ii<W.num_linears; ii++) {
		printf("the%dth",ii);
		print_point_to_screen_matlab(W.L[ii],"linear");
	}
	
	printf("\n\n");
}






