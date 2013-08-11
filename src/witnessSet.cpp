#include "witnessSet.hpp"




//use:  call to parse the file witness_set_file, into the struct W.	


int witness_set::witnessSetParse(const boost::filesystem::path witness_set_file, const int num_vars)
{
  

  int ii, jj;
	
  FILE *IN = safe_fopen_read(witness_set_file.c_str());
  
  
	int dim, comp_num;
	int num_patches, patch_size, num_linears, num_pts, num_vars_in_linears;

	
	
  fscanf(IN, "%d %d %d", &num_pts, &dim, &comp_num); scanRestOfLine(IN);
	this->dim = dim;
	this->comp_num = comp_num;
  this->num_pts = this->num_pts = num_pts;
	this->pts_mp=(point_mp *)bmalloc(this->num_pts*sizeof(point_mp));

	
  this->num_variables = num_vars;
	this->num_synth_vars = 0;
	
  for (ii=0; ii < num_pts; ii++) {
    init_point_mp2(this->pts_mp[ii],num_vars,1024);
		this->pts_mp[ii]->size = num_vars;
    //read the witness points into memory
    for (jj=0; jj < num_vars; jj=jj+1) {
			mpf_inp_str(this->pts_mp[ii]->coord[jj].r, IN, 10); // 10 is the base
			mpf_inp_str(this->pts_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
  }
  
	fscanf(IN, "%d %d", &num_linears, &num_vars_in_linears);  scanRestOfLine(IN);

	this->num_linears = num_linears;
	this->L_mp = (vec_mp *)bmalloc(num_linears*sizeof(vec_mp));
	
  for (ii=0; ii < num_linears; ii++) {
		init_vec_mp2(this->L_mp[ii],num_vars_in_linears,1024);
		this->L_mp[ii]->size = num_vars_in_linears;
		
    //read the witness linears into memory
    for (jj=0; jj < num_vars_in_linears; jj++) {
			mpf_inp_str(this->L_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(this->L_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
  }
  
	
	fscanf(IN, "%d %d", &num_patches, &patch_size); scanRestOfLine(IN);
	this->num_patches = num_patches;
	
	if (this->num_patches>3) {
		std::cerr << num_patches << " patches detected.  this probably indicates a problem." << std::endl;
		std::cerr << "the file being read: " << witness_set_file << std::endl;
		std::cerr << "trying to read " << num_vars << " variables." << std::endl;
		mypause();
	}
	
	this->patch_mp = (vec_mp *)bmalloc(num_patches*sizeof(vec_mp));
	
  for (ii=0; ii < num_patches; ii++) {
		init_vec_mp2(this->patch_mp[ii],patch_size,1024);//default max_prec is 1024
		
		this->patch_mp[ii]->size = patch_size;
    //read the patch into memory
    for (jj=0; jj < patch_size; jj++) {
			mpf_inp_str(this->patch_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(this->patch_mp[ii]->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
    }
  }
	
  fclose(IN);
	
	


  return 0;
}



void witness_set::add_patch(vec_mp new_patch)
{
	
	if (this->num_patches!=0 && this->patch_mp==NULL) {
		printf("trying to add patch to witness set with non-zero num_patches and NULL container!\n");
		deliberate_segfault();
	}
	
	if (this->num_patches==0 && this->patch_mp!=NULL) {
		printf("trying to add point to witness set with num_pts==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (this->num_patches==0) {
		this->patch_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		this->patch_mp = (vec_mp *)br_realloc(this->patch_mp, (this->num_patches+1) * sizeof(vec_mp));
	}
	
	init_vec_mp(this->patch_mp[this->num_patches], new_patch->size);
	this->patch_mp[this->num_patches]->size = new_patch->size;
	vec_cp_mp(this->patch_mp[this->num_patches], new_patch);
	
	this->num_patches++;
	
	return;
}



void witness_set::add_point(vec_mp new_point)
{
	
	if (this->num_pts!=0 && this->pts_mp==NULL) {
		printf("trying to add point to witness set with non-zero num_pts and NULL container!\n");
		deliberate_segfault();
	}
	
	if (this->num_pts==0 && this->pts_mp!=NULL) {
		printf("trying to add point to witness set with num_pts==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (this->num_pts==0) {
		this->pts_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		this->pts_mp = (vec_mp *)br_realloc(this->pts_mp, (this->num_pts+1) * sizeof(vec_mp));
	}
	
	init_vec_mp(this->pts_mp[this->num_pts], new_point->size);
	this->pts_mp[this->num_pts]->size = new_point->size;
	vec_cp_mp(this->pts_mp[this->num_pts], new_point);
	
	this->num_pts++;
	
	return;
}


void witness_set::add_linear(vec_mp new_linear)
{
	
	if (this->num_linears!=0 && this->L_mp==NULL) {
		printf("trying to add linear to witness set with non-zero num_linears and NULL container!\n");
		deliberate_segfault();
	}
	
	if (this->num_linears==0 && this->L_mp!=NULL) {
		printf("trying to add linear to witness set with num_linears==0 and non-NULL container!\n");
		deliberate_segfault();
	}
	
	
	if (this->num_linears==0) {
		this->L_mp = (vec_mp *)br_malloc(sizeof(vec_mp));
	}
	else{
		this->L_mp = (vec_mp *)br_realloc(this->L_mp, (this->num_linears+1) * sizeof(vec_mp));
	}
	
	init_vec_mp(this->L_mp[this->num_linears], new_linear->size);
	this->L_mp[this->num_linears]->size = new_linear->size;
	vec_cp_mp(this->L_mp[this->num_linears], new_linear);
	
	this->num_linears++;
	
	return;
}






//requires the witness set to have set the number of variables.
void witness_set::get_variable_names()
{
	
	this->variable_names.resize(this->num_variables);
	
	std::ifstream fin("names.out");
	
	for (int ii=0; ii<this->num_variables; ++ii){
		fin >> this->variable_names[ii];
	}
	fin.close();

}





void witness_set::write_homogeneous_coordinates(boost::filesystem::path filename)
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%d\n\n",this->num_pts); // print the header line
	
	for (ii=0; ii<this->num_pts; ++ii) {
		for (jj=0; jj<this->pts_mp[ii]->size; jj++) {
			print_mp(OUT,0,&this->pts_mp[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void witness_set::write_dehomogenized_coordinates(boost::filesystem::path filename)
{
	int ii,jj;
	
//	std::cout << "dehomogenizing " << num_variables << " vars, " << num_synth_vars << " synth vars" << std::endl;
	vec_mp result; init_vec_mp(result,1);
	
	
	FILE *OUT = safe_fopen_write(filename.c_str()); // open the output file.
	
	fprintf(OUT,"%d\n\n",this->num_pts); // print the header line
	for (ii=0; ii<this->num_pts; ++ii) {
//		print_point_to_screen_matlab(pts_mp[ii],"hompt");
		if (this->num_synth_vars>0) {
			dehomogenize(&result,this->pts_mp[ii], num_variables-num_synth_vars);
		}
		else{
			dehomogenize(&result,this->pts_mp[ii]);
		}
		
		
//		print_point_to_screen_matlab(result,"soln");
		
		for (jj=0; jj<this->num_variables-1; jj++) {
			print_mp(OUT, 0, &result->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
	
	clear_vec_mp(result);
	return;
}




void witness_set::write_linears(boost::filesystem::path filename)
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename.c_str()); // open the output file
	
	
	fprintf(OUT,"%d\n\n",this->num_linears); // print the header line
	
	for (ii=0; ii<this->num_linears; ++ii) {
		for (jj=0; jj<this->L_mp[ii]->size; jj++) {
			print_mp(OUT, 0, &this->L_mp[ii]->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}



void witness_set::print_patches(boost::filesystem::path filename)
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename); // open the output file
	
	
	fprintf(OUT,"%d\n\n",this->num_patches); // print the header line
	
	for (ii=0; ii<this->num_patches; ++ii) {
		fprintf(OUT,"%d\n",this->patch_mp[ii]->size);
		for (jj=0; jj<this->patch_mp[ii]->size; jj++) {
			print_mp(OUT, 0, &this->patch_mp[ii]->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void witness_set::read_patches_from_file(boost::filesystem::path filename)
{

	FILE *IN = safe_fopen_read(filename);
	
	witness_set::reset_patches();
	int curr_num_patches;
	fscanf(IN,"%d\n",&curr_num_patches);
	
	vec_mp temp_patch; init_vec_mp2(temp_patch,1,1024); temp_patch->size = 1;
	for (int ii=0; ii<curr_num_patches; ii++) {
		int curr_size;
		fscanf(IN,"%d\n",&curr_size);
		change_size_vec_mp(temp_patch,curr_size); temp_patch->size = curr_size;
		
		for (int jj=0; jj<curr_size; jj++) {
			mpf_inp_str(temp_patch->coord[jj].r, IN, 10);
			mpf_inp_str(temp_patch->coord[jj].i, IN, 10);
			scanRestOfLine(IN);
		}
		
		witness_set::add_patch(temp_patch);
	}
	
	clear_vec_mp(temp_patch);
	fclose(IN);
	
	return;
}



void witness_set::print_to_screen()
{
	
	std::stringstream varname;
	
	std::cout << "witness set has " << num_variables << " total variables, " << num_synth_vars << " synthetic vars." << std::endl;
	int ii;
	printf("******\n%d points\n******\n",this->num_pts);
	for (ii=0; ii<this->num_pts; ii++) {
		varname << "point_" << ii;
		print_point_to_screen_matlab(this->pts_mp[ii],varname.str());
		varname.str("");
	}
	
	printf("******\n%d linears\n******\n",this->num_linears);
	
	for (ii=0; ii<this->num_linears; ii++) {
		varname << "linear_" << ii;
		print_point_to_screen_matlab(this->L_mp[ii],varname.str());
		varname.str("");
	}
	
	printf("******\n%d patches\n******\n",this->num_patches);
	
	for (ii=0; ii<this->num_patches; ii++) {
		varname << "patch_" << ii;
		print_point_to_screen_matlab(this->patch_mp[ii],varname.str());
		varname.str("");
	}
	
	std::cout << "variable names:\n";
	for (ii=0; ii<this->variable_names.size(); ii++) {
		std::cout << this->variable_names[ii] << "\n";
	}
	printf("\n\n");
}




void witness_set::only_natural_vars()
{
	witness_set::only_first_vars(this->num_variables - this->num_synth_vars);
}


void witness_set::only_first_vars(int num_vars)
{
	
	vec_mp tempvec;  init_vec_mp2(tempvec, num_vars, 1024);
	tempvec->size = num_vars;
	
	for (int ii=0; ii<this->num_pts; ii++) {
		
		for (int jj=0; jj<num_vars; jj++) {
			set_mp(&tempvec->coord[jj], &this->pts_mp[ii]->coord[jj]);
		}
		
		change_size_vec_mp(this->pts_mp[ii], num_vars);  this->pts_mp[ii]->size = num_vars;
		vec_cp_mp(this->pts_mp[ii], tempvec);
	}
	

	this->num_synth_vars = this->num_synth_vars - (this->num_variables - num_vars); // is this line correct?
	this->num_variables = num_vars;
	
	int patch_size_counter = 0, trim_from_here =  0;
	for (int ii=0; ii<this->num_patches; ii++) {
		patch_size_counter += this->patch_mp[ii]->size;
		if (patch_size_counter == num_vars)
		{
			trim_from_here = ii+1;
		}
	}
	
	if (trim_from_here==0) {
		std::cerr << "problem: the sum of the patche sizes never equalled the number of variables to trim to...\nhence, the trimming operation could not complete." << std::endl;
		deliberate_segfault();
	}
	
	for (int ii=0; ii<this->num_patches; ii++) {
		clear_vec_mp(this->patch_mp[ii]);
	}
	
	this->patch_mp = (vec_mp *) br_realloc(this->patch_mp, trim_from_here* sizeof(vec_mp));
	this->num_patches = trim_from_here;
	
	clear_vec_mp(tempvec);
	return;
}



void witness_set::sort_for_real(tracker_config_t T)
{
	
	
	int ii;
	
	
	int real_indicator[this->num_pts];
	int counter = 0;
	
//	std::cout << "sorting for real\n\n" << num_synth_vars << " synth vars " << num_variables-num_synth_vars << " is diff" << std::endl;
	
	
	vec_mp result; init_vec_mp(result,num_variables-num_synth_vars-1);
	result->size = num_variables-num_synth_vars-1;
	
	for (ii=0; ii<this->num_pts; ii++) {
		for (int jj=1; jj<this->num_variables-this->num_synth_vars; jj++) {
			div_mp(&result->coord[jj-1], &this->pts_mp[ii]->coord[jj], &this->pts_mp[ii]->coord[0]);
		}
//		dehomogenize(&result,this->pts_mp[ii], num_variables-num_synth_vars);
		
//		print_point_to_screen_matlab(result,"isthisreal???");
//		dehomogenize(&result, this->pts_mp[ii]);
		real_indicator[ii] = checkForReal_mp(result, T.real_threshold);
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	vec_mp *tempvec = (vec_mp *)br_malloc(counter * sizeof(vec_mp));
	
	counter = 0;  // reset
	for (ii=0; ii<this->num_pts; ii++) {
		if (real_indicator[ii]==1) {
			
			init_vec_mp(tempvec[counter],this->num_variables); tempvec[counter]->size = this->num_variables;
			vec_cp_mp(tempvec[counter],this->pts_mp[ii]);
			counter++;
		}
		else{
			
		}
	}
	
	clear_vec_mp(result);
	
	
	for (ii=0; ii<this->num_pts; ii++) {
		clear_vec_mp(this->pts_mp[ii]);
	}
	free(this->pts_mp);
	
	this->pts_mp = tempvec;
	this->num_pts = counter;
	
	return;
}






// T is necessary for the tolerances.
void witness_set::sort_for_unique(tracker_config_t T)
{


	
	int curr_uniqueness;
	int num_good_pts = 0;
	std::vector<int> is_unique;
	
	for (int ii = 0; ii<this->num_pts; ++ii) {
		curr_uniqueness = 1;

		for (int jj=ii+1; jj<this->num_pts; ++jj) {
			if ( isSamePoint_homogeneous_input(this->pts_mp[ii],this->pts_mp[jj]) ){
				curr_uniqueness = 0;
			}
		}
		
		if (curr_uniqueness==1) {
			is_unique.push_back(1);
			num_good_pts++;
		}
		else {
			is_unique.push_back(0);
		}
		
	}
	
	
	
	vec_mp *transferme = (vec_mp *)bmalloc(num_good_pts*sizeof(vec_mp));
	int counter = 0;
	for (int ii=0; ii<this->num_pts; ++ii) {
		if (is_unique[ii]==1) {
			init_vec_mp2(transferme[counter],this->num_variables,1024);  transferme[counter]->size = this->num_variables;
			vec_cp_mp(transferme[counter], this->pts_mp[ii]);
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		printf("counter mismatch\n");
		exit(270);
	}
	
	for (int ii=0; ii<this->num_pts; ii++) {
		clear_vec_mp(this->pts_mp[ii]);
	}
	free(this->pts_mp);
	
	
	this->num_pts = num_good_pts;
	this->pts_mp = transferme;
	
	return;
}










void witness_set::merge(const witness_set & W_in){
	
	//error checking first
	if ( (this->num_variables==0) && (W_in.num_variables!=0) && (this->num_pts==0)) {
		this->num_variables = W_in.num_variables;
	}
	else if ( (this->num_variables!=0) && (W_in.num_variables==0) ){
		
	}
	else if ( (W_in.num_variables!=this->num_variables) && (W_in.num_pts!=0) && (this->num_pts!=0) ) {
		printf("merging two witness sets with differing numbers of variables.\n");
		deliberate_segfault();
	}
	
//	if (W_in.num_patches != this->num_patches) {
//		printf("merging two witness sets with differing numbers of patch equations. %d left, %d right\n",
//					 W_in.num_patches, this->num_patches);
//		deliberate_segfault();
//	}
	
	if (W_in.num_synth_vars != this->num_synth_vars) {
		printf("merging two witness sets with differing numbers of synthetic variables. %d merging set, %d existing\n",
					 W_in.num_synth_vars, this->num_synth_vars);
		deliberate_segfault();
	}
	

	for (int ii = 0; ii<W_in.num_linears; ii++) {
		int is_new = 1;
		for (int jj = 0; jj<this->num_linears; jj++){
			if (isSamePoint_inhomogeneous_input(this->L_mp[jj], W_in.L_mp[ii])) {
				is_new = 0;
				break;
			}
			else{
				
			}
		}
		
		if (is_new==1)
			witness_set::add_linear(W_in.L_mp[ii]);
	}
	
	for (int ii = 0; ii<W_in.num_pts; ii++) {
		int is_new = 1;
		for (int jj = 0; jj<this->num_pts; jj++){
			if (isSamePoint_inhomogeneous_input(this->pts_mp[jj], W_in.pts_mp[ii])) {
				is_new = 0;
				break;
			}
		}
		
		if (is_new==1)
			witness_set::add_point(W_in.pts_mp[ii]);
	}
	
	
	for (int ii = 0; ii<W_in.num_patches; ii++) {
		int is_new = 1;
		
		for (int jj = 0; jj<this->num_patches; jj++){
			if (this->patch_mp[jj]->size == W_in.patch_mp[ii]->size) {
				if (isSamePoint_inhomogeneous_input(this->patch_mp[jj], W_in.patch_mp[ii])) {
					is_new = 0;
					break;
				}
			}
		}
		
		if (is_new==1)
			witness_set::add_patch(W_in.patch_mp[ii]);
	}
	
	

	return;
}//re: merge_witness_sets


//copies names from old to new
void cp_names(witness_set *W_out, witness_set & W_in)
{
	
	if (W_in.num_variables==0) {
		printf("\nattempting to copy variable names from witness_set with no variables\n");
		exit(1333);
	}
	
	
	W_out->variable_names.clear();
	
	for (int ii=0; ii<W_in.variable_names.size(); ii++)
		W_out->variable_names.push_back(W_in.variable_names[ii]);
	
}




//copies the mp and d linears from in to out.
void cp_linears(witness_set *W_out, witness_set & W_in)
{
	int ii;
	
	
	W_out->num_linears = W_in.num_linears;
	
	
	if (W_out->L_mp==NULL) {
		W_out->L_mp = (vec_mp *)bmalloc(W_in.num_linears * sizeof(vec_mp));
	}
	else
	{
		W_out->L_mp = (vec_mp *)br_realloc(W_out->L_mp, W_in.num_linears * sizeof(vec_mp));
	}
	
	for (ii=0; ii<W_in.num_linears; ++ii) {
		init_vec_mp(W_out->L_mp[ii],W_in.L_mp[ii]->size);
		vec_cp_mp(W_out->L_mp[ii],W_in.L_mp[ii]);
		W_out->L_mp[ii]->size = W_in.L_mp[ii]->size;
	}
	
	
	return;
}


void cp_patches(witness_set *W_out, witness_set & W_in)
{
	
	W_out->num_patches = W_in.num_patches;
	
	
	if (W_out->patch_mp==NULL)
		W_out->patch_mp = (vec_mp *)bmalloc(W_in.num_patches * sizeof(vec_mp));
	else{
		W_out->patch_mp = (vec_mp *)br_realloc(W_out->patch_mp, W_in.num_patches * sizeof(vec_mp));
	}
	
	for (int ii=0; ii<W_in.num_patches; ++ii) {
		init_vec_mp(W_out->patch_mp[ii],W_in.patch_mp[ii]->size);
		vec_cp_mp(W_out->patch_mp[ii],W_in.patch_mp[ii]);
		W_out->patch_mp[ii]->size = W_in.patch_mp[ii]->size;
	}
	
	
	return;
}

//
////send in an initialized but empty witness_set.
//// T is necessary for the tolerances.
//void sort_for_real(witness_set *W_out,
//									 witness_set & W_in,
//									 tracker_config_t T)
//{
//	int ii;
//	
//	
//	W_out->MPType = W_in.MPType;
//	W_out->incidence_number = W_in.incidence_number;
//	W_out->num_variables = W_in.num_variables;
//
//	
//	//copy the linears, patches from old to new
//	cp_patches(W_out, W_in);
//	cp_linears(W_out, W_in);
//	cp_names(W_out, W_in);
//	
//	int real_indicator[W_in.num_pts];
//	int counter = 0;
//	vec_mp result; init_vec_mp(result,1);
//	for (ii=0; ii<W_in.num_pts; ii++) {
//		dehomogenize(&result, W_in.pts_mp[ii]);
//		real_indicator[ii] = checkForReal_mp(result, T.real_threshold);
//		if (real_indicator[ii]==1) {
//			counter++;
//		}
//	}
//	
//	
//	W_out->num_pts = counter;
//	
//	W_out->pts_mp = (point_mp *)bmalloc(counter*sizeof(point_mp));
//
//	
//	
//
//	counter = 0;  // reset
//	for (ii=0; ii<W_in.num_pts; ii++) {
//		if (real_indicator[ii]==1) {
//
//			init_vec_mp(W_out->pts_mp[counter],W_in.num_variables);
//			W_out->pts_mp[counter]->size = W_in.num_variables;
//
//			vec_cp_mp(W_out->pts_mp[counter],W_in.pts_mp[ii]);
//
//			counter++;
//			
//		}
//		else{
//
//		}
//	}
//	
//	clear_vec_mp(result);
//	
//	
//	return;
//}
//
//
//
//
//
//
////send in an initialized but empty witness_set.
//// T is necessary for the tolerances.
//void sort_for_unique(witness_set *W_out,
//										 witness_set & W_in,
//										 tracker_config_t T)
//{
//
//	int ii, jj;
//	
//
//	
//	W_out->MPType = W_in.MPType;
//	W_out->incidence_number = W_in.incidence_number;
//	W_out->num_variables = W_in.num_variables;
//	
//	//copy the linears, patches from old to new
//	cp_patches(W_out, W_in);
//	cp_linears(W_out, W_in);
//	cp_names(W_out, W_in);
//
//	int curr_uniqueness;
//	int num_good_pts = 0;
//	int *is_unique  = (int *)bmalloc(W_in.num_pts*sizeof(int));
//	
//	for (ii = 0; ii<W_in.num_pts; ++ii) {
//		curr_uniqueness = 1;
//		if (ii!= (W_in.num_pts-1) ) { // the last point is unique...  always
//																	//			printf("a");
//			for (jj=ii+1; jj<W_in.num_pts; ++jj) {
//				if ( isSamePoint_homogeneous_input(W_in.pts_mp[ii],W_in.pts_mp[jj]) ){
//					curr_uniqueness = 0;
//				}
//			}
//		}
//		
//		if (curr_uniqueness==1) {
//			is_unique[ii] = 1;
//			num_good_pts++;
//		}
//		else
//		{
//			is_unique[ii] = 0;
//		}
//		
//	}
//	
//	
//	W_out->num_pts = num_good_pts;
//	
//	W_out->pts_mp = (vec_mp *)bmalloc(num_good_pts*sizeof(vec_mp));
//	int counter = 0;
//	for (ii=0; ii<W_in.num_pts; ++ii) {
//		if (is_unique[ii]==1) {
//			init_vec_mp2(W_out->pts_mp[counter],W_in.num_variables,1024);  W_out->pts_mp[counter]->size = W_in.num_variables;
//			vec_cp_mp(W_out->pts_mp[counter], W_in.pts_mp[ii]);
//			counter++;
//		}
//	}
//	
//	if (counter!= num_good_pts) {
//		printf("counter mismatch\n");
//		exit(270);
//	}
//	
//	free(is_unique);
//	return;
//}
//





