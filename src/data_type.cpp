#include "data_type.hpp"
#include "programConfiguration.hpp"



std::string enum_lookup(int flag)
{
	switch (flag) {
		case SUCCESSFUL:
			return "SUCCESSFUL";
			break;
		
		case CRITICAL_FAILURE:
			return "CRITICAL_FAILURE";
			break;
			
		case TOLERABLE_FAILURE:
			return "TOLERABLE_FAILURE";
			break;
			
		case UNSET:
			return "UNSET";
			break;
			
		case CRITICAL:
			return "CRITICAL";
			break;
			
		case NEW:
			return "NEW";
			break;
			
		case MIDPOINT:
			return "MIDPOINT";
			break;
			
		case ISOLATED:
			return "ISOLATED";
			break;
			
			
		
			
		case NULLSPACE:
			return "NULLSPACE";
			break;
			
		case LINPRODTODETJAC:
			return "LINPRODTODETJAC";
			break;
			
		case DETJACTODETJAC:
			return "DETJACTODETJAC";
			break;
			
		case LINTOLIN:
			return "LINTOLIN";
			break;
			
		case MULTILIN:
			return "MULTILIN";
			break;
			
		case TERMINATE:
			return "TERMINATE";
			break;
			
		case INITIAL_STATE:
			return "INITIAL_STATE";
			break;
			
			
		case PARSING:
			return "PARSING";
			break;
			
		case TYPE_CONFIRMATION:
			return "TYPE_CONFIRMATION";
			break;
			
		case DATA_TRANSMISSION:
			return "DATA_TRANSMISSION";
			break;
			
		case NUMPACKETS:
			return "NUMPACKETS";
			break;
			
		case INACTIVE:
			return "INACTIVE";
			break;
			
		case VEC_MP:
			return "VEC_MP";
			break;
			
		case VEC_D:
			return "VEC_D";
			break;
			
		case MAT_MP:
			return "MAT_MP";
			break;
			
		case MAT_D:
			return "MAT_D";
			break;
			
		case COMP_MP:
			return "COMP_MP";
			break;
		
		case COMP_D:
			return "COMP_D";
			break;
			
		case INDICES:
			return "INDICES";
			break;
			
			
		default:
			break;
	}
	
	return "unknown...  check out data_type.cpp";
}

void * br_malloc(size_t size)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does malloc with error checking                        *
 \***************************************************************/
{
  if (size <= 0)
  { // nothing to allocate
    return NULL;
  }
  else
  { // try to allocate memory
    void *x = malloc(size);
    if (x == NULL)
    {
//			raise(SIGINT);
      printf("ERROR: bertini_real's malloc was unable to allocate memory (%d)!\n", (int) size);
			deliberate_segfault();
//      br_exit(ERROR_MEMORY_ALLOCATION);
    }
    return x;
  }
}

void *br_realloc(void *ptr, size_t size)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does realloc with error checking                       *
 \***************************************************************/
{
  if (size <= 0)
  { // nothing to allocate - free memory and return NULL
    free(ptr);
    ptr = NULL;
  }
  else
  { // try to reallocate memory
    ptr = realloc(ptr, size);
    if (ptr == NULL)
    {
      printf("ERROR: bertini_real's realloc was unable to re-allocate memory!\n");
			deliberate_segfault();
//      br_exit(ERROR_MEMORY_ALLOCATION);
    }
  }
  return ptr;
}


void deliberate_segfault(){
	printf("the following segfault is deliberate\n");
	int *faulty = NULL;
	faulty[-10] = faulty[10]+faulty[0];
	
}


void br_exit(int errorCode)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: exits Bertini - either standard or using MPI           *
 \***************************************************************/
{
  if (errorCode == 0)
    errorCode = ERROR_OTHER;
	
	
  printf("%s\n", "bertini_real quitting\n");
	
	deliberate_segfault();
	
#ifdef _HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, errorCode);
#else
  exit(errorCode);
#endif
}
























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
	this->pts_mp=(point_mp *)br_malloc(this->num_pts*sizeof(point_mp));
	
	
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
	this->L_mp = (vec_mp *)br_malloc(num_linears*sizeof(vec_mp));
	
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
	
	this->patch_mp = (vec_mp *)br_malloc(num_patches*sizeof(vec_mp));
	
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
	
	init_vec_mp2(this->patch_mp[this->num_patches], new_patch->size,1024);
	this->patch_mp[this->num_patches]->size = new_patch->size;
	vec_cp_mp(this->patch_mp[this->num_patches], new_patch);
	
	this->num_patches++;
	
	
	int numasdf = 0;
	for (int ii=0; ii<this->num_patches; ii++) {
		numasdf += this->patch_mp[ii]->size;
	}
	
	if (numasdf>this->num_variables) {
		std::cout << "added patch " << this->num_patches << ", but put the number of patch vars (" << numasdf << ") higher than set number (" << this->num_variables << ")" << std::endl;
		deliberate_segfault();
	}
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
	
	init_vec_mp2(this->pts_mp[this->num_pts], new_point->size,1024);
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
	
	init_vec_mp2(this->L_mp[this->num_linears], new_linear->size,1024);
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





void witness_set::write_homogeneous_coordinates(boost::filesystem::path filename) const
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

void witness_set::write_dehomogenized_coordinates(boost::filesystem::path filename) const
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




void witness_set::write_linears(boost::filesystem::path filename) const
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



void witness_set::print_patches(boost::filesystem::path filename) const
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



void witness_set::print_to_screen() const
{
	
	vec_mp dehom;  init_vec_mp(dehom,1); dehom->size = 1;
	
	std::stringstream varname;
	
	std::cout << "witness set has " << num_variables << " total variables, " << num_synth_vars << " synthetic vars." << std::endl;
	int ii;
	printf("******\n%d points\n******\n",this->num_pts);
	std::cout << color::green();
	for (ii=0; ii<this->num_pts; ii++) {
		
		dehomogenize(&dehom, this->pts_mp[ii], num_variables - num_synth_vars);
		
		varname << "point_" << ii;
		
		print_point_to_screen_matlab(dehom,varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::blue();
	printf("******\n%d linears\n******\n",this->num_linears);
	
	for (ii=0; ii<this->num_linears; ii++) {
		varname << "linear_" << ii;
		print_point_to_screen_matlab(this->L_mp[ii],varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << color::cyan();
	printf("******\n%d patches\n******\n",this->num_patches);
	
	for (ii=0; ii<this->num_patches; ii++) {
		varname << "patch_" << ii;
		print_point_to_screen_matlab(this->patch_mp[ii],varname.str());
		varname.str("");
	}
	std::cout << color::console_default();
	
	std::cout << "variable names:\n";
	for (ii=0; ii< int(this->variable_names.size()); ii++) {
		std::cout << this->variable_names[ii] << "\n";
	}
	printf("\n\n");
	
	
	clear_vec_mp(dehom);
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
		std::cerr << "problem: the sum of the patch sizes never equalled the number of variables to trim to...\nhence, the trimming operation could not complete." << std::endl;
		deliberate_segfault();
	}
	
	for (int ii=trim_from_here; ii<this->num_patches; ii++) {
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
	
	
	int *real_indicator = new int[this->num_pts];
	int counter = 0;
	
	vec_mp result; init_vec_mp(result,num_variables-num_synth_vars-1);
	result->size = num_variables-num_synth_vars-1;
	
	for (ii=0; ii<this->num_pts; ii++) {
		for (int jj=1; jj<this->num_variables-this->num_synth_vars; jj++) {
			div_mp(&result->coord[jj-1], &this->pts_mp[ii]->coord[jj], &this->pts_mp[ii]->coord[0]);
		}
		real_indicator[ii] = checkForReal_mp(result, T.real_threshold);
		
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	vec_mp *tempvec = (vec_mp *)br_malloc(counter * sizeof(vec_mp));
	
	counter = 0;  // reset
	for (ii=0; ii<this->num_pts; ii++) {
		if (real_indicator[ii]==1) {
			
			init_vec_mp2(tempvec[counter],this->num_variables,1024); tempvec[counter]->size = this->num_variables;
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
	
	delete[] real_indicator;
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
	
	
	
	vec_mp *transferme = (vec_mp *)br_malloc(num_good_pts*sizeof(vec_mp));
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










void witness_set::merge(const witness_set & W_in)
{
	
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



int witness_set::compute_downstairs_crit_midpts(vec_mp crit_downstairs,
																								 vec_mp midpoints_downstairs,
																								 std::vector< int > & index_tracker,
																								 vec_mp pi) const
{
	
	int retVal = 0;
	
	
	vec_mp projection_values; init_vec_mp2(projection_values,this->num_pts,1024);
	projection_values->size = this->num_pts;
	
	for (int ii=0; ii<this->num_pts; ii++){
		
		for (int jj=0; jj<this->pts_mp[ii]->size; jj++) {
			
			if (!(mpfr_number_p(this->pts_mp[ii]->coord[jj].r) && mpfr_number_p(this->pts_mp[ii]->coord[jj].i))) {
				std::cout << color::red();
				std::cout << "there was NAN in a point in the projections to sort :(" << std::endl;
				print_point_to_screen_matlab(this->pts_mp[ii], "bad_point");
				
				print_point_to_screen_matlab(pi,"pi");
				std::cout << color::console_default();
				return -14;
			}
			
		}
	
		
		
		projection_value_homogeneous_input(&projection_values->coord[ii], this->pts_mp[ii], pi); // set projection value
		
		print_comp_matlab(&projection_values->coord[ii],"proj_val");
		
		
		if (!(mpfr_number_p(projection_values->coord[ii].r) && mpfr_number_p(projection_values->coord[ii].i)))
		{
			print_comp_matlab(&projection_values->coord[ii],"not_a_number");
			print_point_to_screen_matlab(pi,"pi");
			
		}
			
//		if (fabs(mpf_get_d(projection_values->coord[ii].i))<1e-13) {
//			mpf_set_str(projection_values->coord[ii].i,"0",10);
//		}
	}
	
	
	change_size_vec_mp(crit_downstairs,1); // destructive resize
	crit_downstairs->size = 1;
	
	retVal = sort_increasing_by_real(crit_downstairs, index_tracker, projection_values);
	
	clear_vec_mp(projection_values);
	
	int num_midpoints = crit_downstairs->size - 1;
	
	if (num_midpoints<1) {
		return retVal;
	}
	
	
	comp_d h;  h->r = 0.5; h->i = 0.0;
	comp_mp half; init_mp2(half,1024); d_to_mp(half, h);
	
	change_size_vec_mp(midpoints_downstairs, num_midpoints);
	midpoints_downstairs->size = num_midpoints;
	
	comp_mp temp; init_mp2(temp,1024);
	for (int ii=0; ii<num_midpoints; ii++){
		add_mp(temp, &crit_downstairs->coord[ii], &crit_downstairs->coord[ii+1]);
		mul_mp(&midpoints_downstairs->coord[ii], temp, half);
	}
	
	clear_mp(temp);
	clear_mp(half);
	
	
	return retVal;
}




//copies names from old to new
void witness_set::cp_names(const witness_set & W_in)
{
	
	if (W_in.num_variables==0) {
		printf("\nattempting to copy variable names from witness_set with no variables\n");
		deliberate_segfault();
		exit(1333);
	}
	
	this->variable_names = W_in.variable_names;

}




//copies the mp and d linears from in to out.
void witness_set::cp_linears(const witness_set & W_in)
{
	int ii;
	
	

	
	if (this->num_linears==0)
		this->L_mp = (vec_mp *)br_malloc(W_in.num_linears * sizeof(vec_mp));
	else
		this->L_mp = (vec_mp *)br_realloc(this->L_mp, W_in.num_linears * sizeof(vec_mp));
	
	for (ii=0; ii<W_in.num_linears; ++ii) {
		init_vec_mp2(this->L_mp[ii],W_in.L_mp[ii]->size,1024);
		vec_cp_mp(this->L_mp[ii],W_in.L_mp[ii]);
		this->L_mp[ii]->size = W_in.L_mp[ii]->size;
	}
	
	this->num_linears = W_in.num_linears;
	
	return;
}


void witness_set::cp_patches(const witness_set & W_in)
{
	if (this->num_patches==0)
		this->patch_mp = (vec_mp *)br_malloc(W_in.num_patches * sizeof(vec_mp));
	else
		this->patch_mp = (vec_mp *)br_realloc(this->patch_mp, W_in.num_patches * sizeof(vec_mp));
	
	
	for (int ii=0; ii<W_in.num_patches; ++ii) {
		init_vec_mp2(this->patch_mp[ii],W_in.patch_mp[ii]->size,1024);
		vec_cp_mp(this->patch_mp[ii],W_in.patch_mp[ii]);
		this->patch_mp[ii]->size = W_in.patch_mp[ii]->size;
	}
	
	this->num_patches = W_in.num_patches;
	return;
}














int vertex_set::add_vertex(const vertex source_vertex){
	
	this->num_vertices++;
	this->vertices.push_back(source_vertex);
	
	return this->num_vertices-1;
}


void vertex_set::print_to_screen()
{
	printf("vertex set has %d vertices:\n\n",this->num_vertices);
	for (int ii=0; ii<this->num_vertices; ++ii) {
		print_point_to_screen_matlab(this->vertices[ii].pt_mp,"vert");
		print_point_to_screen_matlab(this->vertices[ii].projection_values,"projection_values");
		printf("type: %d\n", this->vertices[ii].type);
	}
}


int vertex_set::setup_vertices(boost::filesystem::path INfile)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices;
	int num_vars;
	fscanf(IN, "%d %d\n\n", &num_vertices, &this->num_projections);
	
	vertex temp_vertex;
	
	for(int ii=0;ii<num_vertices;ii++)
	{
		fscanf(IN, "%d\n", &num_vars);
		if (temp_vertex.pt_mp->size != num_vars) {
			change_size_vec_mp(temp_vertex.pt_mp,num_vars); temp_vertex.pt_mp->size = num_vars;
		}
		
		for (int jj=0;jj<num_vars;jj++)
		{
			mpf_inp_str(temp_vertex.pt_mp->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vertex.pt_mp->coord[jj].i, IN, 10);
		}
		
		increase_size_vec_mp(temp_vertex.projection_values,num_projections);
		for (int jj=0; jj<num_projections; jj++) {
			mpf_inp_str(temp_vertex.projection_values->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vertex.projection_values->coord[jj].i, IN, 10);
		}
		
		
		fscanf(IN,"%d\n",&temp_vertex.type);
		
		
		vertex_set::add_vertex(temp_vertex);
	}
	
	
	fclose(IN);
	
	if (this->num_vertices!=num_vertices) {
		printf("parity error in num_vertices.\n\texpected: %d\tactual: %d\n",this->num_vertices,num_vertices);
		br_exit(25943);
	}
	
	return num_vertices;
}






/**
 
 //assumes all vertices have the same number of variables in them.
 
 Output vertex structure as follows:
 # pts
 pt.1
 
 pt.2
 
 **/
void vertex_set::print(boost::filesystem::path outputfile)
{
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d %d\n\n",num_vertices,num_projections);
	for (int ii = 0; ii < num_vertices; ii++)
	{ // output points
		fprintf(OUT,"%d\n", vertices[ii].pt_mp->size);
		for(int jj=0;jj<vertices[ii].pt_mp->size;jj++) {
			print_mp(OUT, 0, &vertices[ii].pt_mp->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		for(int jj=0;jj<vertices[ii].projection_values->size;jj++) {
			print_mp(OUT, 0, &vertices[ii].projection_values->coord[jj]);
			fprintf(OUT,"\n");
		}
		
		fprintf(OUT,"\n");
		fprintf(OUT,"%d\n\n",vertices[ii].type);
	}
	
	fclose(OUT);
	
}





int decomposition::add_vertex(vertex_set & V, vertex source_vertex)
{
	
	
	int current_index = V.add_vertex(source_vertex);
	
	if (this->counters.find(source_vertex.type) == this->counters.end()) {
		this->counters[source_vertex.type] = 0;
	}

	this->counters[source_vertex.type]++;
	this->indices[source_vertex.type].push_back(current_index);

	return current_index;
}






int decomposition::index_in_vertices(vertex_set &V,
																		 vec_mp testpoint,
																		 comp_mp proj_val,
																		 tracker_config_t T)
{
	int ii;
	int index = -1;
		
	
	
	for (int jj=1; jj<V.num_natural_variables; jj++) {
		div_mp(&V.checker_1->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
	}
	
//	WTB: a faster comparison search.
	for (ii=0; ii<V.num_vertices; ii++) {
		
		int current_index = ii;
		
		if (V.vertices[current_index].removed!=1) {
		
		//it would be desirable to perform the below comparison, but the projection value depends on the projection used to compute it!
		// hence, it has been commented out.  
//			sub_mp(V.diff, projection_value, V.vertices[current_index].projVal_mp);
//			mpf_abs_mp(V.abs, V.diff);
//			
//			if (mpf_cmp(V.abs, V.zerothresh) > 0){ // i think this is opposite
//				continue;
//			}
		
		for (int jj=1; jj<V.num_natural_variables; jj++) {
			div_mp(&V.checker_2->coord[jj-1],&V.vertices[current_index].pt_mp->coord[jj], &V.vertices[current_index].pt_mp->coord[0]);
		}
		
		if (isSamePoint_inhomogeneous_input(V.checker_1, V.checker_2)){				
			index = current_index;
			break;
		}

		if (index!=-1)
			break;
			
		}
	}
		
	
	return index;	
}





//TODO: make T here a pointer

int decomposition::index_in_vertices_with_add(vertex_set &V,
																							vertex vert,
																							comp_mp proj_val,
																							tracker_config_t T)
{
	int index = decomposition::index_in_vertices(V, vert.pt_mp, proj_val, T);
	
	if (index==-1) {
		index = decomposition::add_vertex(V, vert);
	}
	
	return index;
	
}











int decomposition::setup(boost::filesystem::path INfile,
												 boost::filesystem::path & inputName,
												 boost::filesystem::path directoryName)
{
	

	boost::filesystem::path input_deflated_Name;
	
	std::stringstream converter;
	std::string tempstr;
	std::ifstream fin(INfile.c_str());
	
	int num_vertices = 0;
	
	getline(fin, tempstr);
	input_deflated_Name = tempstr;
	inputName = directoryName / input_deflated_Name;

	getline(fin, tempstr);
	converter << tempstr;
	converter >> this->num_variables >> this->dimension;
	converter.clear(); converter.str("");
	
	int num_types;
	fin >> num_types;
	
	for (int ii =0; ii<num_types; ii++) {
		int current_type, num_this_type, current_index;
		fin >> current_type >> num_this_type;
		std::cout << "vertex type: " << current_type << ", # vertices: " << num_this_type << std::endl;
		this->counters[current_type] = num_this_type;
		
		for (int jj=0; jj<num_this_type; jj++) {
			fin >> current_index;
			this->indices[current_type].push_back(current_index);
		}
	}
	
	vec_mp tempvec; init_vec_mp2(tempvec, this->num_variables,1024);
	tempvec->size = this->num_variables;
	
	getline(fin,tempstr); // burn two lines as a consequence of using getline.
	getline(fin,tempstr);
	for (int ii=0; ii<dimension; ii++) {
		for (int jj=0;jj<this->num_variables;jj++)
		{
			getline(fin,tempstr);
			converter << tempstr;
			std::string re, im;
			converter >> re >> im; // this line is correct
			converter.clear(); converter.str("");
			
			mpf_set_str(tempvec->coord[jj].r, const_cast<char *>(re.c_str()), 10);
			mpf_set_str(tempvec->coord[jj].i, const_cast<char *>(im.c_str()), 10);
		}
		
		decomposition::add_projection(tempvec);
		
	}

	
	
	
	int curr_num_patches;
	
	getline(fin,tempstr); std::cout << tempstr << std::endl;
	getline(fin,tempstr); std::cout << tempstr << std::endl;
	getline(fin,tempstr); std::cout << tempstr << std::endl; // i hate this
	converter << tempstr;
	converter >> curr_num_patches;
	converter.clear(); converter.str("");
	
	std::cout << "should be loading " << curr_num_patches << " patches" << std::endl;
	
	
	vec_mp temp_patch; init_vec_mp2(temp_patch,1,1024); temp_patch->size = 1;
	for (int ii=0; ii<curr_num_patches; ii++) {
		
		getline(fin,tempstr); 
		getline(fin,tempstr); // <--- supposed to have the size in it
		
		int curr_size;
		
		converter << tempstr;
		converter >> curr_size;
		converter.clear(); converter.str("");
		
		change_size_vec_mp(temp_patch,curr_size); temp_patch->size = curr_size;
		
		for (int jj=0; jj<curr_size; jj++) {
			getline(fin,tempstr);
			converter << tempstr;
			std::string re, im;
			converter >> re >> im; // this line is correct
			converter.clear(); converter.str("");
			
			mpf_set_str(temp_patch->coord[jj].r, const_cast<char *>(re.c_str()), 10);
			mpf_set_str(temp_patch->coord[jj].i, const_cast<char *>(im.c_str()), 10);
		}
		
		decomposition::add_patch(temp_patch);
		print_point_to_screen_matlab(temp_patch,"temp_patch");
	}
	

	
	clear_vec_mp(temp_patch);
	
	fin.close();
	
	
	
	return num_vertices;
}





/**
 Output curve overall info as follows:
 
 **/
void decomposition::print(boost::filesystem::path base)
{
	int ii;
	
	boost::filesystem::create_directory(base);
	
	FILE *OUT = safe_fopen_write(base / "decomp");
	
	fprintf(OUT,"%s\n",input_filename.filename().c_str());
	
	fprintf(OUT,"%d %d\n",num_variables, dimension);
	
	fprintf(OUT, "%d\n", int(this->counters.size()));
	
	std::map< int, int>::iterator type_iter;
	for (type_iter = this->counters.begin(); type_iter!= this->counters.end(); type_iter++) {
		fprintf(OUT, "%d %d\n",type_iter->first, type_iter->second);  // print the number corresponding to the type, and the number of that type.
		
		for (int jj = 0; jj<type_iter->second; jj++) {
			fprintf(OUT, "%d\n", indices[type_iter->first][jj]);
		}
		fprintf(OUT,"\n");
	}
	
	if (dimension != num_curr_projections) {
//		std::cout << "decomposition was short projections\nneeded	" << this->dimension << " but had " << num_curr_projections << std::endl;;
	}
	
	for (ii=0; ii<num_curr_projections; ii++) {
		for(int jj=0;jj<pi[ii]->size;jj++)
		{
			print_mp(OUT, 0, &pi[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
		fprintf(OUT,"\n");
	}
	fprintf(OUT,"\n");
	
	fprintf(OUT,"%d\n\n",this->num_patches); // print the header line
	
	for (ii=0; ii<this->num_patches; ++ii) {
		fprintf(OUT,"%d\n",this->patch[ii]->size);
		for (int jj=0; jj<this->patch[ii]->size; jj++) {
			print_mp(OUT, 0, &this->patch[ii]->coord[jj]);
			fprintf(OUT, "\n");
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT);
}






void decomposition::output_main(const BR_configuration & program_options, vertex_set & V)
{
	
	FILE *OUT;
	std::stringstream converter;
	converter << "_dim_" << this->dimension << "_comp_" << this->component_num;
	boost::filesystem::path base = program_options.output_dir;
	base += converter.str();
	converter.clear(); converter.str("");
	
	boost::filesystem::remove_all(base);
	boost::filesystem::create_directory(base);
	
	
	copyfile("witness_data",base / "witness_data");
	
	copyfile(program_options.witness_set_filename, base / "witness_set");
	
	
	// TODO:  this should be a write call, not a copy ?
	copyfile(program_options.input_deflated_filename, base / program_options.input_deflated_filename.filename());
	
	
	V.print(base / "V.vertex");
	
	this->print(base);
	
	OUT = safe_fopen_write("Dir_Name");
	fprintf(OUT,"%s\n",base.c_str());
	fprintf(OUT,"%d\n",program_options.MPType);
	fclose(OUT);
	
}





void clear_sample(sample_data *S, int MPType)
{
	
	int ii;
	for(ii=0;ii<S->num_edges;ii++)
	{
		free(S->sample_indices[ii]);
	}
	free(S->sample_indices);
}



void norm_of_difference(mpf_t result, vec_mp left, vec_mp right)
{
	if (left->size!=right->size || left->size == 0) {
		printf("attempting to take difference of two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		
		print_point_to_screen_matlab(left,"left");
		print_point_to_screen_matlab(right,"right");
		deliberate_segfault();
	}
	
	int ii;
	
	vec_mp difference;  init_vec_mp2(difference, left->size,1024);difference->size = left->size;
	comp_mp temp; init_mp2(temp,1024);
	
	for (ii = 0;  ii< left->size; ++ii) {
		sub_mp(&difference->coord[ii], &left->coord[ii], &right->coord[ii]);
	}
	
	twoNormVec_mp(difference, temp);
	
	mpf_abs_mp(result, temp);
	clear_vec_mp(difference);
	clear_mp(temp);
	return;
}


void dehomogenize(point_d *result, point_d dehom_me)
{
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

void dehomogenize(point_mp *result, point_mp dehom_me)
{
	if (dehom_me->size==0 || dehom_me->size==1) {
		printf("attempting to dehomogenize a vector of length 0 or 1\n");
		exit(977);
	}
	
	comp_mp denom; init_mp2(denom,1024);
	change_size_vec_mp((*result),dehom_me->size-1);
	
	(*result)->size = dehom_me->size-1;
	
	set_mp(denom, &dehom_me->coord[0]);
	int ii;
	for (ii=0; ii<dehom_me->size-1; ++ii) {
		set_mp( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
		div_mp(&(*result)->coord[ii],&(*result)->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	clear_mp(denom);
	return;
}

void dehomogenize(point_mp *result, point_mp dehom_me, int num_variables)
{
	if (dehom_me->size==0 || dehom_me->size==1) {
		printf("attempting to dehomogenize a vector of length 0 or 1\n");
		exit(977);
	}
	
	comp_mp denom; init_mp2(denom,1024);
	change_size_vec_mp((*result),dehom_me->size-1);
	
	(*result)->size = dehom_me->size-1;
	
	set_mp(denom, &dehom_me->coord[0]);
	int ii;
	for (ii=0; ii<num_variables-1; ++ii) {
		set_mp( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
		div_mp(&(*result)->coord[ii],&(*result)->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	for (ii=num_variables-1; ii<dehom_me->size-1; ++ii) {
		set_mp( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
	}
	
	clear_mp(denom);
	return;
}




void dehomogenize(point_d *result, point_d dehom_me, int num_variables)
{
	if (dehom_me->size==0 || dehom_me->size==1) {
		printf("attempting to dehomogenize a vector of length 0 or 1\n");
		exit(977);
	}
	
	comp_d denom;
	change_size_point_d((*result),dehom_me->size-1);
	
	(*result)->size = dehom_me->size-1;
	
	set_d(denom, &dehom_me->coord[0]);
	int ii;
	for (ii=0; ii<num_variables-1; ++ii) {
		set_d( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
		div_d(&(*result)->coord[ii],&(*result)->coord[ii],denom); //  result[ii] = dehom_me[ii+1]/dehom_me[0].
	}
	
	for (ii=num_variables-1; ii<dehom_me->size-1; ++ii) {
		set_d( &(*result)->coord[ii],&(dehom_me)->coord[ii+1]);
	}
	
	return;
}



void dot_product_d(comp_d result, vec_d left, vec_d right)
{
	if (left->size!=right->size) {
		printf("attempting to dot_d two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		deliberate_segfault();
	}
	
	set_zero_d(result);
	
	comp_d temp;
	int ii;
	for (ii=0; ii<left->size; ++ii) {
		mul_d(temp,&left->coord[ii],&right->coord[ii]);
		add_d(result,result,temp);
	}
}

void dot_product_mp(comp_mp result, vec_mp left, vec_mp right)
{
	if (left->size!=right->size) {
		printf("attempting to dot_mp two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		deliberate_segfault();
	}
	
	set_zero_mp(result);
	comp_mp temp; init_mp2(temp,1024);
	int ii;
	for (ii=0; ii<left->size; ++ii) {
		mul_mp(temp,&left->coord[ii],&right->coord[ii]);
		add_mp(result,result,temp);
	}
	clear_mp(temp);
}


/**
 computes the projection value given a homogeneous input.
 
 double type
 */
void projection_value_homogeneous_input(comp_d result, vec_d input, vec_d projection)
{
	set_zero_d(result);
	comp_d temp;
	for (int ii=0; ii<projection->size; ii++) {
		mul_d(temp, &input->coord[ii], &projection->coord[ii]);
		add_d(result, result, temp);
	}
	set_d(temp, result);
	div_d(result, temp, &input->coord[0]);
//
//	if (result->i < 1e-14) {
//		result->i = 0.0;
//	}
	
//	if (result->r < 1e-13) {
//		result->r = 0.0;
//	}
	
	return;
}

/**
 computes the projection value given a homogeneous input.
 
 mp type
*/
void projection_value_homogeneous_input(comp_mp result, vec_mp input, vec_mp projection)
{
	set_zero_mp(result);
	comp_mp temp; init_mp2(temp,1024);
	for (int ii=0; ii<projection->size; ii++) {
		mul_mp(temp, &input->coord[ii], &projection->coord[ii]);
		add_mp(result, result, temp);
	}
	set_mp(temp, result);
	div_mp(result, temp, &input->coord[0]);
	
	clear_mp(temp);
	
//	if (mpf_get_d(result->i) < 1e-14) {
//		mpf_set_d(result->i, 0.0);
//	}
	
//	mp_to_d(temp2, result);
//	if (temp2->r < 1e-13) {
//		mpf_set_d(result->r, 0.0);
//	}
	
}



int isSamePoint_inhomogeneous_input(point_d left, point_d right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_inhom_d with disparate sized points.\n");
		std::cout << "left: " << left->size << "\t right: " << right->size << std::endl;
		deliberate_segfault();
//		exit(-287);
	}
	
	
	int indicator = isSamePoint(left,NULL,52,right,NULL,52,SAMEPOINTTOL);
	
	
	return indicator;
}


int isSamePoint_inhomogeneous_input(point_mp left, point_mp right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_inhom_mp with disparate sized points.\n");
		std::cout << "left: " << left->size << "\t right: " << right->size << std::endl;
		deliberate_segfault();
//		exit(-287);
	}
	
	int indicator = isSamePoint(NULL,left,65,NULL,right,65,SAMEPOINTTOL); // make the bertini library call
	

	return indicator;
}



int isSamePoint_homogeneous_input(point_d left, point_d right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_hom_d with disparate sized points.\n");
		std::cout << "left: " << left->size << "\t right: " << right->size << std::endl;
		deliberate_segfault();
//		exit(-287);
	}
	
	vec_d dehom_left;  init_vec_d(dehom_left,left->size-1);  dehom_left->size = left->size-1;
	vec_d dehom_right; init_vec_d(dehom_right,right->size-1); dehom_right->size = right->size-1;
	
	dehomogenize(&dehom_left,left);
	dehomogenize(&dehom_right,right);
	
	int indicator = isSamePoint(dehom_left,NULL,52,dehom_right,NULL,52,SAMEPOINTTOL);
	
	clear_vec_d(dehom_left); clear_vec_d(dehom_right);
	
	return indicator;
}


int isSamePoint_homogeneous_input(point_mp left, point_mp right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_hom_mp with disparate sized points.\n");
		std::cout << "left: " << left->size << "\t right: " << right->size << std::endl;
		deliberate_segfault();
//		exit(-287);
	}
	
	vec_mp dehom_left;  init_vec_mp(dehom_left,left->size-1);  dehom_left->size = left->size-1;
	vec_mp dehom_right; init_vec_mp(dehom_right,right->size-1); dehom_right->size = right->size-1;
	
	dehomogenize(&dehom_left,left);
	dehomogenize(&dehom_right,right);
	
	int indicator = isSamePoint(NULL,dehom_left,65,NULL,dehom_right,65,SAMEPOINTTOL); // make the bertini library call
	
	clear_vec_mp(dehom_left); clear_vec_mp(dehom_right);
	
	return indicator;
}




void print_point_to_screen_matlab(vec_d M, std::string name)
{
	
	printf("%s = [...\n",name.c_str());
	for (int kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		printf(" %.15le+1i*%.15le;\n",M->coord[kk].r,M->coord[kk].i);
	}
	printf("];\n\n");
}

void print_point_to_screen_matlab(vec_mp M, std::string name)
{
	
	printf("%s = [...\n",name.c_str());
	for (int kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		mpf_out_str (NULL, 10, 16, M->coord[kk].r);
		printf("+1i*");
		mpf_out_str (NULL, 10, 16, M->coord[kk].i);
		printf(";\n");
	}
	printf("];\n\n");
}



void print_matrix_to_screen_matlab(mat_d M, std::string name)
{
	int jj,kk;
	
	printf("%%matrix '%s' has dimensions %dx%d\n", name.c_str(), M->rows,M->cols);
	printf("%s = [...\n",name.c_str());
	for (kk = 0; kk < M->rows; kk++)
	{ // print kth row
		for (jj = 0; jj < M->cols; jj++)
		{
//			if (abs(M->entry[kk][jj].r)<1e-12) {
//				printf("0");
//			}
//			else{
//			printf("%.8le",M->entry[kk][jj].r);
//			}
//			
//			if (abs(M->entry[kk][jj].r)>=1e-12 && abs(M->entry[kk][jj].i)>=1e-12) {
//				printf("+");
//			}
//			else{
//				printf(" ");
//			}
//			
//			if (abs(M->entry[kk][jj].i)<1e-12) {
//				
//			}
//			else{
//				printf("1i*%.8le ",M->entry[kk][jj].i);
//			}
			printf("%.4le+1i*%.4le ",M->entry[kk][jj].r,M->entry[kk][jj].i );
		}
		if (kk!= M->rows-1) {
			printf(";...\n");
		}
		
	}
	printf("];\n\n");
}
void print_matrix_to_screen_matlab(mat_mp M, std::string name)
{
	int jj,kk;
	
	printf("%%matrix '%s' has dimensions %dx%d\n",name.c_str(), M->rows,M->cols);
	printf("%s = [...\n",name.c_str());
	for (kk = 0; kk < M->rows; kk++)
	{ // print kth row
		for (jj = 0; jj < M->cols; jj++)
		{
			
			mpf_out_str (NULL, 10, 15, M->entry[kk][jj].r);
			printf("+1i*");
			mpf_out_str (NULL, 10, 15, M->entry[kk][jj].i); // base 10 , 7 digits
			printf("\t");
		}
		printf(";\n");
	}
	printf("];\n\n");
}


void print_comp_matlab(comp_mp M, std::string name){
	printf("%s=",name.c_str());
	mpf_out_str (NULL, 10, 4, M->r);
	printf("+1i*");
	mpf_out_str (NULL, 10, 4, M->i); // base 10, 6 digits
	printf("\n");
	return;
}

void print_comp_matlab(comp_d M, std::string name){
	printf("%s=%.5le+1i*%.5le\n",name.c_str(),M->r,M->i);
	return;
}

void print_path_retVal_message(int retVal){
	
	
	if (retVal==100) {
		printf("max_prec_reached\n");
	}
	if (retVal==-50) {
		printf("reached_minTrackT\nrelevant setting name is 'NbhdRadius'\n");
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



int get_num_vars_PPD(preproc_data PPD){
	int num_vars = 0; // initialize
	
	int ii;
	
	//run through each variable group
	for (ii=0; ii<(PPD.num_var_gp+PPD.num_hom_var_gp); ii++) {
		num_vars+= PPD.size[ii];
		num_vars+= PPD.type[ii];
	}
	
	return num_vars;
	
}



void cp_patch_mp(patch_eval_data_mp *PED, patch_eval_data_mp PED_input)
{
  PED->num_patches = PED_input.num_patches;
	
  // set the current precision
  PED->curr_prec = PED_input.curr_prec;
	
  // initialize patchCoeff to this preicision
  init_mat_mp2(PED->patchCoeff, PED_input.patchCoeff->rows, PED_input.patchCoeff->cols, PED->curr_prec);
  init_mat_rat(PED->patchCoeff_rat, PED_input.patchCoeff->rows, PED_input.patchCoeff->cols);
	
  // setup patchCoeff
  mat_cp_mp_rat(PED->patchCoeff, PED->patchCoeff_rat, PED_input.patchCoeff, PED_input.patchCoeff_rat);
	
  return;
}


void cp_patch_d(patch_eval_data_d *PED, patch_eval_data_d PED_input)
{
  PED->num_patches = PED_input.num_patches;
	

  // initialize patchCoeff to this preicision
  init_mat_d(PED->patchCoeff, PED_input.patchCoeff->rows, PED_input.patchCoeff->cols);

  // setup patchCoeff
  mat_cp_d(PED->patchCoeff, PED_input.patchCoeff);
	
  return;
}

void cp_preproc_data(preproc_data *PPD, const preproc_data & PPD_input)
{
	
  PPD->num_funcs = PPD_input.num_funcs;
  PPD->num_hom_var_gp = PPD_input.num_hom_var_gp;
  PPD->num_var_gp = PPD_input.num_var_gp;
	
  int total_gp = PPD->num_hom_var_gp + PPD->num_var_gp;
	
	std::cout << PPD_input.num_hom_var_gp << " " << PPD_input.num_var_gp  << " "  << total_gp << std::endl;

	
	if (total_gp==0) {
		return;
	}
	else{
		PPD->size = (int *)br_malloc(total_gp * sizeof(int));
		PPD->type = (int *)br_malloc(total_gp * sizeof(int));
		
		for (int i = 0; i < total_gp; i++)
		{
			PPD->size[i] = PPD_input.size[i];
			PPD->type[i] = PPD_input.type[i];
		}
	}
		
  
	
  return;
}


void clear_post_process_t(post_process_t * endPoint, int num_vars)
{

	if (endPoint->sol_prec >1)
	{
		if (endPoint->sol_prec >= 64)
		{ // clear _mp
			mpf_clear(endPoint->function_resid_mp);
			mpf_clear(endPoint->newton_resid_mp);
			for (int j = 0; j < num_vars; j++)
			{
				clear_mp(endPoint->sol_mp[j]);
			}
		}
		free(endPoint->sol_mp);
		free(endPoint->sol_d);
	}
}

// this sort should be optimized.  it is sloppy and wasteful right now.
int sort_increasing_by_real(vec_mp projections_sorted, std::vector< int > & index_tracker, vec_mp projections_input){
	

	
	
	for (int ii=0; ii<projections_input->size; ii++) {
		if (!(mpfr_number_p(projections_input->coord[ii].r) && mpfr_number_p(projections_input->coord[ii].i))) {
			std::cout << "there was NAN in the projections to sort :(" << std::endl;
			print_point_to_screen_matlab(projections_input, "projections_input");
			
			return -51;
		}
	}
	
	
	
	
	
	std::vector< int > index_tracker_non_unique;
	std::vector< double > projvals_as_doubles;
	
	
	vec_mp projections_sorted_non_unique;
	init_vec_mp2(projections_sorted_non_unique,projections_input->size,1024);
	projections_sorted_non_unique->size = projections_input->size;


	
	std::set<int> unsorted_indices;
	for (int ii=0; ii<projections_input->size; ii++) {
		unsorted_indices.insert(ii);
	}
	
	
	
	
	//sort by size
	for (int ii=0; ii<projections_input->size; ii++) { // for each of the projection values input
		double min = 1e20; // reset this bogus value
		
		int indicator = -1;
		
		// this loop finds the minimum projection value
		for (std::set<int>::iterator set_iter = unsorted_indices.begin(); set_iter!=unsorted_indices.end(); set_iter++) {
			
			double curr = mpf_get_d(projections_input->coord[*set_iter].r); // convert projection value to a double for comparison
			if ( curr < min) { // compare
				indicator = *set_iter;
				min = curr;
			}
		}
		
		if (indicator==-1) { // if min value was larger than a huge number
			printf("min projection value was *insanely* large\n");
			exit(1111);
		}
		
		unsorted_indices.erase(indicator);
		
		projvals_as_doubles.push_back(min);
		index_tracker_non_unique.push_back(indicator);
		set_mp( &projections_sorted_non_unique->coord[ii],&projections_input->coord[indicator]);
	}
	
	
	
	
	// filter for uniqueness
	
	
	double distinct_thresh = 1e-5;  // reasonable?
	
	change_size_vec_mp(projections_sorted,1); projections_sorted->size = 1;
	
	index_tracker.push_back(index_tracker_non_unique[0]);
	set_mp(&projections_sorted->coord[0],&projections_sorted_non_unique->coord[0])
	int unique_counter = 1;
	for (int ii=1; ii<projections_input->size; ii++) {
		if ( fabs( projvals_as_doubles[ii-1]-projvals_as_doubles[ii]) < distinct_thresh) {
			continue;
		}
		else
		{
			increase_size_vec_mp(projections_sorted,unique_counter+1); projections_sorted->size = unique_counter+1;
			set_mp(&projections_sorted->coord[unique_counter],&projections_sorted_non_unique->coord[ii]);
			unique_counter++;
			
			index_tracker.push_back(index_tracker_non_unique[ii]);
		}
	}
	

	return 0;
}



//input the raw number of variables including the homogeneous variable (of which there must be one)
// assume the array of integers 'randomized_degrees' is already initialized to the correct size.
void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, std::vector< int > & randomized_degrees, int num_desired_rows, int num_funcs)
{
	int ii,jj;
	
	
	
	//get unique degrees
	int *degrees = (int *) br_malloc(num_funcs*sizeof(int));
	int *unique_degrees = (int *) br_malloc(num_funcs*sizeof(int));
	
	
	FILE *IN = safe_fopen_read("deg.out"); //open the deg.out file for reading.
	int num_unique_degrees = 0;
	int occurrence_counter;
	for (ii=0; ii<num_funcs; ++ii) {
		fscanf(IN,"%d\n",&degrees[ii]); // read data
		occurrence_counter = 0; // set the counter for how many timmes the current degree has already been found.
		for (jj=0; jj<ii; jj++) {
			if (degrees[jj]==degrees[ii]) { // if previously stored degree is same as current one
				occurrence_counter++; // increment counter
			}
		}
		
		if (occurrence_counter==0) { // if did not find already in list
			unique_degrees[num_unique_degrees] = degrees[ii]; // add to list of unique degrees.
			num_unique_degrees++; // have one more unique degree
		} // re: jj
	}// re: ii
	fclose(IN);
	
	
	if (num_desired_rows==num_funcs) {
		make_matrix_ID_mp(randomization_matrix,num_funcs,num_funcs);
		for (ii=0; ii<num_desired_rows; ++ii) {
			randomized_degrees.push_back(degrees[ii]);
		}
		free(degrees);
		free(unique_degrees);
		return;
	}
	
	//sort the unique degrees into decreasing order
	qsort(unique_degrees, num_unique_degrees, sizeof(int), compare_integers_decreasing);
	
	//count how many of each unique degree there are.
	int *num_of_each_degree = (int *) br_malloc(num_unique_degrees*sizeof(int));
	for (ii=0; ii<num_unique_degrees; ii++) {
		num_of_each_degree[ii] = 0;
		for (jj=0; jj<num_funcs; ++jj) {
			if (unique_degrees[ii]==degrees[jj]) {
				num_of_each_degree[ii]++;
			}
		}
	}
	
	
	
	//	for (ii=0; ii<num_unique_degrees; ii++) {
	//		printf("unique_degrees[%d]=%d; num_of_each_degree=%d\n",ii,unique_degrees[ii],num_of_each_degree[ii]);
	//	}
	
	//resize the matrix
	change_size_mat_mp(randomization_matrix,num_desired_rows,num_funcs);
	randomization_matrix->rows = num_desired_rows; randomization_matrix->cols = num_funcs;
	
	
	
	int counter = 0;
	int current_degree_index = 0; // start at the end
	int current_degree;
	for (ii=0; ii<num_desired_rows; ii++) {
		
		counter++;
		if (counter>num_of_each_degree[current_degree_index]) {
			current_degree_index++;
			counter = 1;
		}
		
		current_degree = unique_degrees[current_degree_index];
		randomized_degrees.push_back(current_degree);
		
		int encountered_current_degree = 0;
		for (jj=0; jj<num_funcs; jj++) {
			if ( (degrees[jj]<= current_degree)  ) {
				encountered_current_degree++;
				if (encountered_current_degree >= counter){
					get_comp_rand_real_mp(&randomization_matrix->entry[ii][jj]);
				}
				else{
					set_zero_mp(&randomization_matrix->entry[ii][jj]);
				}
			}
			else
			{
				set_zero_mp(&randomization_matrix->entry[ii][jj]);
			}
		}
		
		
	}
	
	free(num_of_each_degree);
	free(degrees);
	free(unique_degrees);
	
	return;
}


int compare_integers_decreasing(const void * left_in, const void * right_in){
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left<right) {
		return 1;
	}
	else if(right > left){
		return -1;
	}
	else{
		return 0;
	}
	
}


void send_patch_d(patch_eval_data_d * patch)
{
	comp_d *patch_coeff = NULL;
	patch_eval_data_d_int PED_int;
	MPI_Datatype mpi_comp_d, mpi_patch_d_int;
	
	// setup mpi_comp_d & mpi_patch_d_int
	create_comp_d(&mpi_comp_d);
	create_patch_eval_data_d_int(&mpi_patch_d_int);
	// setup PED_int
	cp_patch_d_int(&PED_int, patch, &patch_coeff, 0);
	
	// broadcast patch structures
	MPI_Bcast(&PED_int, 1, mpi_patch_d_int, 0, MPI_COMM_WORLD);
	MPI_Bcast(patch_coeff, PED_int.patchCoeff_rows * PED_int.patchCoeff_cols, mpi_comp_d, 0, MPI_COMM_WORLD);
	
	// free memory
	MPI_Type_free(&mpi_comp_d);
	MPI_Type_free(&mpi_patch_d_int);
	free(patch_coeff);
}


void receive_patch_d(patch_eval_data_d * patch)
{
	comp_d *patch_coeff = NULL;
	patch_eval_data_d_int PED_int;
	MPI_Datatype mpi_comp_d, mpi_patch_d_int;
	
	// setup mpi_comp_d & mpi_patch_d_int
	create_comp_d(&mpi_comp_d);
	create_patch_eval_data_d_int(&mpi_patch_d_int);
	
	// recv patch structures
	MPI_Bcast(&PED_int, 1, mpi_patch_d_int, 0, MPI_COMM_WORLD);
	// setup patch_coeff
	patch_coeff = (comp_d *)br_malloc(PED_int.patchCoeff_rows * PED_int.patchCoeff_cols * sizeof(comp_d));
	MPI_Bcast(patch_coeff, PED_int.patchCoeff_rows * PED_int.patchCoeff_cols, mpi_comp_d, 0, MPI_COMM_WORLD);
	
	
	// setup patch
	cp_patch_d_int(patch, &PED_int, &patch_coeff, 1);  // patch_coeff is freed in here
	
	// free mpi_comp_d & mpi_patch_d_int
	MPI_Type_free(&mpi_comp_d);
	MPI_Type_free(&mpi_patch_d_int);
}


void send_patch_mp(patch_eval_data_mp * patch)
{
	char *patchStr = NULL;
	patch_eval_data_mp_int PED_int;
	MPI_Datatype mpi_patch_int;
	
	// setup mpi_patch_int
	create_patch_eval_data_mp_int(&mpi_patch_int);
	// setup PED_int
	cp_patch_mp_int(&PED_int, patch, &patchStr, 0, 0);
	
	// send PED_int
	MPI_Bcast(&PED_int, 1, mpi_patch_int, 0, MPI_COMM_WORLD);
	// send patchStr
	MPI_Bcast(patchStr, PED_int.totalLength, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// clear memory
	free(patchStr);
	MPI_Type_free(&mpi_patch_int);
}

void receive_patch_mp(patch_eval_data_mp * patch)
{
	char *patchStr = NULL;
	patch_eval_data_mp_int PED_int;
	MPI_Datatype mpi_patch_int;
	
	// setup mpi_patch_int
	create_patch_eval_data_mp_int(&mpi_patch_int);
	// recv PED_int
	MPI_Bcast(&PED_int, 1, mpi_patch_int, 0, MPI_COMM_WORLD);
	
	// setup patchStr
	patchStr = (char *)br_malloc(PED_int.totalLength * sizeof(char));
	// recv patchStr
	MPI_Bcast(patchStr, PED_int.totalLength, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// setup _mp patch
	cp_patch_mp_int(patch, &PED_int, &patchStr, 1, 1);
	// last 3 arguments:   ,    ,   intype -- 0 if sender, !0, else
	
	
	// free mpi_patch_int (patchStr is freed in cp_patch_mp_int)
	MPI_Type_free(&mpi_patch_int);
}



void send_preproc_data(preproc_data *PPD){
	

	int *buffer = new int[3];
	
	buffer[0] = PPD->num_funcs;
	buffer[1] = PPD->num_hom_var_gp;
	buffer[2] = PPD->num_var_gp;
	
	MPI_Bcast(buffer, 3, MPI_INT, 0, MPI_COMM_WORLD);
	
	int size = PPD->num_hom_var_gp + PPD->num_var_gp;
	MPI_Bcast(PPD->type, size, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(PPD->size, size, MPI_INT, 0, MPI_COMM_WORLD);
	
	delete(buffer);
}

void receive_preproc_data(preproc_data *PPD){
	
	
	int *buffer = new int[3];
	
	MPI_Bcast(buffer, 3, MPI_INT, 0, MPI_COMM_WORLD);
	
	PPD->num_funcs = buffer[0];
  PPD->num_hom_var_gp = buffer[1];
	PPD->num_var_gp = buffer[2];
	
	
	
	int num_groups = PPD->num_hom_var_gp + PPD->num_var_gp;
	
	PPD->type = (int *) br_malloc(num_groups * sizeof(int));
	PPD->size = (int *) br_malloc(num_groups * sizeof(int));
	MPI_Bcast(PPD->type, num_groups, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(PPD->size, num_groups, MPI_INT, 0, MPI_COMM_WORLD);
	
	delete(buffer);
}










void send_vec_mp(vec_mp b, int target)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: broadcasts b                                           *
 \***************************************************************/
{
  MPI_Datatype mpi_vec_mp_int;
  point_mp_int b_int;
  char *bstr = NULL;
	
  // create the datatypes mpi_vec_mp_int
  create_point_mp_int(&mpi_vec_mp_int);
	
	cp_point_mp_int(&b_int, b, &bstr, 0, 0, 0);
	
	// send b_int and bstr
	MPI_Send(&b_int, 1, mpi_vec_mp_int, target, VEC_MP, MPI_COMM_WORLD);
	MPI_Send(bstr, b_int.totalLength, MPI_CHAR, target,  VEC_MP, MPI_COMM_WORLD);
	
	// clear bstr
	free(bstr);


  // clear mpi_vec_mp_int
  MPI_Type_free(&mpi_vec_mp_int);
	
  return;
}



void receive_vec_mp(vec_mp b, int source)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: broadcasts b                                           *
 \***************************************************************/
{
  MPI_Datatype mpi_vec_mp_int;
  point_mp_int b_int;
  char *bstr = NULL;
	
  // create the datatypes mpi_vec_mp_int
  create_point_mp_int(&mpi_vec_mp_int);
	
	MPI_Status statty_mc_gatty;

	MPI_Recv(&b_int, 1, mpi_vec_mp_int, source,  VEC_MP, MPI_COMM_WORLD, &statty_mc_gatty);
	bstr = (char *)br_malloc(b_int.totalLength * sizeof(char));
	MPI_Recv(bstr, b_int.totalLength, MPI_CHAR, source,  VEC_MP, MPI_COMM_WORLD, &statty_mc_gatty);
	
	// setup b and clear bstr
	cp_point_mp_int(b, &b_int, &bstr, 1, 1, 1); 
	
  // clear mpi_vec_mp_int
  MPI_Type_free(&mpi_vec_mp_int);
	
  return;
}




void send_vec_d(vec_d b, int target)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: broadcasts b                                           *
 \***************************************************************/
{
	MPI_Datatype mpi_point_d_int, mpi_comp_d;
  point_d_int b_int;
  comp_d *entries = NULL;
	
  // create the datatype mpi_point_d_int & mpi_comp_d
  create_point_d_int(&mpi_point_d_int);
  create_comp_d(&mpi_comp_d);
	
	cp_point_d_int(&b_int, b, &entries, 0, 0, 0);
	
	// send b_int
	MPI_Send(&b_int, 1, mpi_point_d_int, target, VEC_D, MPI_COMM_WORLD);
	// send entries
	MPI_Send(entries, b_int.size, mpi_comp_d, target, VEC_D, MPI_COMM_WORLD);
	
	// clear entries
	free(entries);

	
  // clear mpi_point_d_int & mpi_comp_d
  MPI_Type_free(&mpi_point_d_int);
  MPI_Type_free(&mpi_comp_d);
	

	
//  MPI_Datatype mpi_vec_d_int;
//  point_d_int b_int;
//  char *bstr = NULL;
//	
//  // create the datatypes mpi_vec_d_int
//  create_point_d_int(&mpi_vec_d_int);
//	
//	cp_point_d_int(&b_int, b, &bstr, 0, 0, 0);
//	
//	// send b_int and bstr
//	MPI_Send(&b_int, 1, mpi_vec_d_int, target, VEC_MP, MPI_COMM_WORLD);
//	MPI_Send(bstr, b_int.totalLength, MPI_CHAR, target,  VEC_MP, MPI_COMM_WORLD);
//	
//	// clear bstr
//	free(bstr);
//	
//	
//  // clear mpi_vec_d_int
//  MPI_Type_free(&mpi_vec_d_int);
	
  return;
}



void receive_vec_d(vec_d b, int source)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: broadcasts b                                           *
 \***************************************************************/
{
	MPI_Datatype mpi_point_d_int, mpi_comp_d;
  point_d_int b_int;
  comp_d *entries = NULL;
	
  // create the datatype mpi_point_d_int & mpi_comp_d
  create_point_d_int(&mpi_point_d_int);
  create_comp_d(&mpi_comp_d);

	MPI_Status statty_mc_gatty;
	
	MPI_Recv(&b_int, 1, mpi_point_d_int, source, VEC_D, MPI_COMM_WORLD, &statty_mc_gatty);
	
	entries = (comp_d *)br_malloc(b_int.size * sizeof(comp_d));
	// recv entries
	MPI_Recv(entries, b_int.size, mpi_comp_d, source, VEC_D, MPI_COMM_WORLD, &statty_mc_gatty);
	
	// setup b
	cp_point_d_int(b, &b_int, &entries, 1, 1, 1);
	
  // clear mpi_point_d_int & mpi_comp_d
  MPI_Type_free(&mpi_point_d_int);
  MPI_Type_free(&mpi_comp_d);
	

	
	
//  MPI_Datatype mpi_vec_d_int;
//  point_d_int b_int;
//  comp_d *bstr = NULL;
//	
//  // create the datatypes mpi_vec_d_int
//  create_point_d_int(&mpi_vec_d_int);
//	
//	MPI_Status statty_mc_gatty;
//	
//	MPI_Recv(&b_int, 1, mpi_vec_d_int, source,  VEC_MP, MPI_COMM_WORLD, &statty_mc_gatty);
//	bstr = (char *)br_malloc(b_int.totalLength * sizeof(char));
//	MPI_Recv(bstr, b_int.totalLength, MPI_CHAR, source,  VEC_MP, MPI_COMM_WORLD, &statty_mc_gatty);
//	
//	// setup b and clear bstr
//	cp_point_d_int(b, &b_int, &bstr, 1, 1, 1);
//	
//  // clear mpi_vec_d_int
//  MPI_Type_free(&mpi_vec_d_int);
	
  return;
}



namespace color {
	
	
	
  std::string console_default(){
		return "\033[0m";
	}
	
	std::string black(){
		return "\033[0;30m";
	}
	
	std::string red(){
		return "\033[0;31m";
	}
	
	std::string green(){
		return "\033[0;32m";
	}
	
	std::string brown(){
		return "\033[0;33m";
	}
	
	std::string blue(){
		return "\033[0;34m";
	}
	
	std::string magenta(){
		return "\033[0;35m";
	}
	std::string cyan(){
		return "\033[0;36m";
	}
	
	std::string gray(){
		return "\033[0;37m";
	}
	
	
}


//black - 30
//red - 31
//green - 32
//brown - 33
//blue - 34
//magenta - 35
//cyan - 36
//lightgray - 37
