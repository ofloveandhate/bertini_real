#include "witnessSet.hpp"




//use:  call to parse the file witness_set_file, into the struct W.	


int witnessSetParse(witness_set *W, boost::filesystem::path witness_set_file, const int num_vars){
  

  int ii, jj;
	
  FILE *IN = safe_fopen_read(witness_set_file.c_str());
  
  
	int dim, comp_num;
	int num_patches, patch_size, num_linears, num_pts, num_vars_in_linears;

	
	
  fscanf(IN, "%d %d %d", &num_pts, &dim, &comp_num); scanRestOfLine(IN);
	W->dim = dim;
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

	W->num_linears = num_linears;
	W->L_d = (vec_d *)bmalloc(num_linears*sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(num_linears*sizeof(vec_mp));
	
  for (ii=0; ii < num_linears; ii++) {
		init_vec_mp2(W->L_mp[ii],num_vars_in_linears,1024); init_vec_d(W->L_d[ii],num_vars_in_linears);
		
		W->L_d[ii]->size =  W->L_mp[ii]->size = num_vars_in_linears;
		
    //read the witness linears into memory
    for (jj=0; jj < num_vars_in_linears; jj++) {
			mpf_inp_str(W->L_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->L_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
		vec_mp_to_d(W->L_d[ii],W->L_mp[ii]);
  }
  
	
	fscanf(IN, "%d %d", &num_patches, &patch_size); scanRestOfLine(IN);
//	W->patch_size = patch_size;
	W->num_patches = num_patches;
	
	W->patch_d = (vec_d *)bmalloc(num_patches*sizeof(vec_d));
	W->patch_mp = (vec_mp *)bmalloc(num_patches*sizeof(vec_mp));
	
  for (ii=0; ii < num_patches; ii++) {
		init_vec_mp2(W->patch_mp[ii],patch_size,1024);//default max_prec is 1024
		init_vec_d(W->patch_d[ii],patch_size);
		
		W->patch_mp[ii]->size = W->patch_d[ii]->size = patch_size;
		
    //read the patch into memory
    for (jj=0; jj < patch_size; jj++) {
			mpf_inp_str(W->patch_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->patch_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
		vec_mp_to_d(W->patch_d[ii],W->patch_mp[ii]);
  }
	
  fclose(IN);
	
	


  return 0;
}



void add_patch_to_witness_set(witness_set *W, vec_mp new_patch){
	
	if (W->num_patches!=0 && W->patch_mp==NULL) {
		printf("trying to add patch to witness set with non-zero num_patches and NULL container!\n");
		deliberate_segfault();
	}
	
	if (W->num_patches==0 && W->patch_mp!=NULL) {
		printf("trying to add point to witness set with num_pts==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (W->num_patches==0) {
		W->patch_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		W->patch_mp = (vec_mp *)brealloc(W->patch_mp, (W->num_patches+1) * sizeof(vec_mp));
	}
	
	init_vec_mp(W->patch_mp[W->num_patches], new_patch->size);
	W->patch_mp[W->num_patches]->size = new_patch->size;
	vec_cp_mp(W->patch_mp[W->num_patches], new_patch);
	
	W->num_patches++;
	
	return;
}



void add_point_to_witness_set(witness_set *W, vec_mp new_point){
	
	if (W->num_pts!=0 && W->pts_mp==NULL) {
		printf("trying to add point to witness set with non-zero num_pts and NULL container!\n");
		deliberate_segfault();
	}
	
	if (W->num_pts==0 && W->pts_mp!=NULL) {
		printf("trying to add point to witness set with num_pts==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (W->num_pts==0) {
		W->pts_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		W->pts_mp = (vec_mp *)brealloc(W->pts_mp, (W->num_pts+1) * sizeof(vec_mp));
	}
	
	init_vec_mp(W->pts_mp[W->num_pts], new_point->size);
	W->pts_mp[W->num_pts]->size = new_point->size;
	vec_cp_mp(W->pts_mp[W->num_pts], new_point);
	
	W->num_pts++;
	
	return;
}


void add_linear_to_witness_set(witness_set *W, vec_mp new_linear){
	
	if (W->num_linears!=0 && W->L_mp==NULL) {
		printf("trying to add linear to witness set with non-zero num_linears and NULL container!\n");
		deliberate_segfault();
	}
	
	if (W->num_linears==0 && W->L_mp!=NULL) {
		printf("trying to add linear to witness set with num_linears==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (W->num_linears==0) {
		W->L_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		W->L_mp = (vec_mp *)brealloc(W->L_mp, (W->num_linears+1) * sizeof(vec_mp));
	}
	
	init_vec_mp(W->L_mp[W->num_linears], new_linear->size);
	W->L_mp[W->num_linears]->size = new_linear->size;
	vec_cp_mp(W->L_mp[W->num_linears], new_linear);
	
	W->num_linears++;
	
	return;
}


void merge_witness_sets(witness_set *W_out,witness_set W_left,witness_set W_right){
	
	//error checking first
	if (W_left.num_variables != W_right.num_variables) {
		printf("merging two witness sets with differing numbers of variables.\n");
		deliberate_segfault();
	}
	
	if (W_left.num_patches != W_right.num_patches) {
		printf("merging two witness sets with differing numbers of patch equations. %d left, %d right\n",
					 W_left.num_patches, W_right.num_patches);
		deliberate_segfault();
	}
	

	
	cp_names(W_out, W_right);
	
	
	int ii;
	//initialize the structure for holding the produced data
	W_out->num_variables = W_left.num_variables;
	
	
	W_out->num_linears = W_left.num_linears+W_right.num_linears;
	W_out->L_d = (vec_d *)bmalloc((W_left.num_linears+W_right.num_linears)*sizeof(vec_d));
	W_out->L_mp = (vec_mp *)bmalloc((W_left.num_linears + W_right.num_linears)*sizeof(vec_mp));
	// merge the left and right linears into the output.
	int counter = 0;
	for (ii=0; ii<W_left.num_linears; ii++) {
		init_vec_d(W_out->L_d[counter],W_out->num_variables); W_out->L_d[counter]->size = W_out->num_variables;
		vec_cp_d(W_out->L_d[counter],W_left.L_d[ii]);
		
		init_vec_mp(W_out->L_mp[counter],W_out->num_variables); W_out->L_mp[counter]->size = W_out->num_variables;
		vec_cp_mp(W_out->L_mp[counter],W_left.L_mp[ii]);
		counter++;
	}
	for (ii=0; ii<W_right.num_linears; ii++) {
		init_vec_d(W_out->L_d[counter],W_out->num_variables);
		vec_cp_d(W_out->L_d[counter],W_right.L_d[ii]);
		
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
	if (W_out->L_d==NULL) {
		W_out->L_d = (vec_d *)bmalloc(W_in.num_linears * sizeof(vec_d));
	}
	else
	{
		W_out->L_d = (vec_d *)brealloc(W_out->L_d, W_in.num_linears * sizeof(vec_d));
	}
	
	for (ii=0; ii<W_in.num_linears; ++ii) {
		init_vec_d(W_out->L_d[ii],W_in.L_d[ii]->size);
		vec_cp_d(W_out->L_d[ii],W_in.L_d[ii]);
		W_out->L_d[ii]->size = W_in.L_d[ii]->size;
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
	
	
	W_out->num_patches = W_in.num_patches;
	
	
	if (W_out->patch_d==NULL) {
		W_out->patch_d = (vec_d *)bmalloc(W_in.num_patches * sizeof(vec_d));
	}
	else
	{
		W_out->patch_d = (vec_d *)brealloc(W_out->patch_d, W_in.num_patches * sizeof(vec_d));
	}
	
	for (ii=0; ii<W_in.num_patches; ++ii) {
		init_vec_d(W_out->patch_d[ii],0);
		change_size_vec_d(W_out->patch_d[ii],W_in.patch_d[ii]->size);
		vec_cp_d(W_out->patch_d[ii],W_in.patch_d[ii]);
		W_out->patch_d[ii]->size = W_in.patch_d[ii]->size;
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
	
	W_out->dim = W_in.dim;
	W_out->comp_num = W_in.comp_num;
	W_out->incidence_number = W_in.incidence_number;
	W_out->num_variables = W_in.num_variables;
	W_out->num_var_gps = W_in.num_var_gps;
	W_out->num_linears = W_in.num_linears;
	W_out->num_patches = W_in.num_patches;
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
void init_witness_set(witness_set *W){
	W->num_variables = W->num_patches = W->num_linears = 0;
	W->patch_d = W->L_d = W->pts_d = NULL;
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





void write_homogeneous_coordinates(witness_set W, boost::filesystem::path filename)
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
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

void write_dehomogenized_coordinates(witness_set W, boost::filesystem::path filename){
	int ii,jj;
	
	FILE *OUT = safe_fopen_write(filename.c_str()); // open the output file.
	
	fprintf(OUT,"%d\n\n",W.num_pts); // print the header line
	for (ii=0; ii<W.num_pts; ++ii) {
			vec_mp result;
			init_vec_mp(result,1);
			dehomogenize_mp(&result,W.pts_mp[ii]);
			for (jj=0; jj<W.num_variables-1; jj++) {
				print_mp(OUT, 0, &result->coord[jj]);
				fprintf(OUT, "\n");
			}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
	return;
}




void write_linears(witness_set W, boost::filesystem::path filename)
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%d\n\n",W.num_linears); // print the header line
	
	for (ii=0; ii<W.num_linears; ++ii) {
		for (jj=0; jj<W.num_variables; jj++) {
			print_mp(OUT, 0, &W.L_mp[ii]->coord[jj]);
			fprintf(OUT, "\n");
//			fprintf(OUT,"%.15le %.15le\n",W.L_mp[ii]->coord[jj].r,W.L_mp[ii]->coord[jj].i);
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}



void clear_witness_set(witness_set W){
	int ii;
	
	
	for (ii=0; ii<W.num_patches; ii++) {
//		clear_vec_d(W.patch_d[ii]);
		clear_vec_mp(W.patch_mp[ii]);
	}
	if (W.num_patches>0) {
//		free(W.patch_d);
		free(W.patch_mp);
	}
	
	
	for (ii=0; ii<W.num_linears; ii++) {
//		clear_vec_d(W.L_d[ii]);
		clear_vec_mp(W.L_mp[ii]);
	}
	if (W.num_linears>0) {
//		free(W.L_d);
		free(W.L_mp);
	}
	
	for (ii=0; ii<W.num_pts; ii++) {
//		clear_point_d(W.pts_d[ii]);
		clear_point_mp(W.pts_mp[ii]);
	}
	if (W.num_pts>0) {
//		free(W.pts_d);
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
	printf("******\n%d points\n******\n",W.num_pts);
	for (ii=0; ii<W.num_pts; ii++) {
		printf("the%dth",ii);
		print_point_to_screen_matlab_mp(W.pts_mp[ii],"point");
	}
	
	printf("******\n%d linears\n******\n",W.num_linears);
	
	for (ii=0; ii<W.num_linears; ii++) {
		printf("the%dth",ii);
		print_point_to_screen_matlab(W.L_d[ii],"linear");
	}
	
	printf("\n\n");
}






//send in an initialized but empty witness_set.
// T is necessary for the tolerances.
void sort_for_real(witness_set *W_out,
									 witness_set W_in,
									 tracker_config_t T)
{
	//	printf("sorting points for real-ness; %d points\n",W_in.num_pts);
	int ii;
	
	//	if (W_in.incidence_number==-1) {
	//		printf("input witness_set has unset incidence_number for comparison.\n");
	//		exit(-1112);
	//	}
	
	W_out->MPType = W_in.MPType;
	W_out->incidence_number = W_in.incidence_number;
	W_out->num_variables = W_in.num_variables;
	W_out->num_var_gps = W_in.num_var_gps;
	
	//copy the linears, patches from old to new
	cp_patches(W_out, W_in);
	cp_linears(W_out, W_in);
	cp_names(W_out, W_in);
	//
	
	int real_indicator[W_in.num_pts];
	int counter = 0;
	vec_mp result; init_vec_mp(result,1);
	for (ii=0; ii<W_in.num_pts; ii++) {
		dehomogenize_mp(&result, W_in.pts_mp[ii]);
		real_indicator[ii] = checkForReal_mp(result, T.real_threshold);
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	W_out->num_pts = counter;
	
	W_out->pts_mp = (point_mp *)bmalloc(counter*sizeof(point_mp));
	W_out->pts_d  = (point_d  *)bmalloc(counter*sizeof(point_d));
	
	if ( (W_out->pts_d == NULL) && (counter!=0) )
		printf("WTF\n");
	
	//	printf("here, W_in.numpts %d, counter:%d\n",W_in.num_pts,counter);
	
	counter = 0;  // reset
	for (ii=0; ii<W_in.num_pts; ii++) {
		//		printf("%d",ii);
		if (real_indicator[ii]==1) {
			//			printf("\treal\n");
			
			
			init_vec_d( W_out->pts_d[counter],   W_in.num_variables);
			//			printf("b!\n");
			W_out->pts_d[counter]->size = W_in.num_variables;
			//			printf("b!\n");
			
			
			init_vec_mp(W_out->pts_mp[counter],W_in.num_variables);
			//			printf("a!\n");
			W_out->pts_mp[counter]->size = W_in.num_variables;
			//			printf("a!\n");
			
			
			
			vec_cp_mp(W_out->pts_mp[counter],W_in.pts_mp[ii]);
			//			printf("c!\n");
			vec_mp_to_d(W_out->pts_d[counter],W_in.pts_mp[ii]);
			//			printf("c!\n");
			counter++;
			
		}
		else{
			//			printf("\tnot real \n");
		}
		
		//		printf("counter %d\n",counter);
		
	}
	
	//	printf("here2\n");
	clear_vec_mp(result);
	
	
	return;
}






//send in an initialized but empty witness_set.
// T is necessary for the tolerances.
void sort_for_unique(witness_set *W_out,
										 witness_set W_in,
										 tracker_config_t T)
{
	//	printf("sorting points for unique-ness\n%d points in \n",W_in.num_pts);
	int ii, jj;
	
	//	if (W_in.incidence_number==-1) {
	//		printf("input witness_set has unset incidence_number for comparison.\n");
	//		exit(-1112);
	//	}
	
	W_out->MPType = W_in.MPType;
	W_out->incidence_number = W_in.incidence_number;
	W_out->num_variables = W_in.num_variables;
	W_out->num_var_gps = W_in.num_var_gps;
	
	//copy the linears, patches from old to new
	cp_patches(W_out, W_in);
	cp_linears(W_out, W_in);
	cp_names(W_out, W_in);
	//
	int curr_uniqueness;
	int num_good_pts = 0;
	int *is_unique  = (int *)bmalloc(W_in.num_pts*sizeof(int));
	
	for (ii = 0; ii<W_in.num_pts; ++ii) {
		curr_uniqueness = 1;
		if (ii!= (W_in.num_pts-1) ) { // the last point is unique...  always
																	//			printf("a");
			for (jj=ii+1; jj<W_in.num_pts; ++jj) {
				if ( isSamePoint_homogeneous_input_d(W_in.pts_d[ii],W_in.pts_d[jj]) ){
					curr_uniqueness = 0;
					
					//					printf("the following two points are not distinct:\n");
					//					print_point_to_screen_matlab(W_in.pts_d[ii],"left");
					//					print_point_to_screen_matlab(W_in.pts_d[jj],"right");
				}
			}
		}
		
		//		printf("%d curr_uniqueness\n",curr_uniqueness);
		
		if (curr_uniqueness==1) {
			is_unique[ii] = 1;
			num_good_pts++;
		}
		else
		{
			is_unique[ii] = 0;
		}
		
	}
	
	
	W_out->num_pts = num_good_pts;
	
	W_out->pts_mp = (vec_mp *)bmalloc(num_good_pts*sizeof(vec_mp));
	W_out->pts_d = (vec_d *)bmalloc(num_good_pts*sizeof(vec_d));
	int counter = 0;
	for (ii=0; ii<W_in.num_pts; ++ii) {
		if (is_unique[ii]==1) {
			init_vec_d(W_out->pts_d[counter],W_in.num_variables); W_out->pts_d[counter]->size = W_in.num_variables;
			init_vec_mp2(W_out->pts_mp[counter],W_in.num_variables,1024);  W_out->pts_mp[counter]->size = W_in.num_variables;
			
			vec_cp_d(W_out->pts_d[counter], W_in.pts_d[ii]);
			vec_cp_mp(W_out->pts_mp[counter], W_in.pts_mp[ii]);
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		printf("counter mismatch\n");
		exit(270);
	}
	
	free(is_unique);
	
	
	//	printf("%d points out\n",W_out->num_pts);
	return;
}







