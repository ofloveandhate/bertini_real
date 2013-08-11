#include "data_type.hpp"


//enum {SUCCESSFUL=0, CRITICAL_FAILURE=-10, TOLERABLE_FAILURE=-1};
////The following lets us use words instead of numbers to indicate vertex type.
//enum {UNSET=-10, CRITICAL=0, NEW=1, MIDPOINT=2, ISOLATED=-1, SAMPLE_POINT=3};
//// enum for worker mode choice



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
			
			
			
			//enum {NULLSPACE = 3000, LINPRODTODETJAC, DETJACTODETJAC, LINTOLIN, MULTILIN};
			//enum {TERMINATE = 2000, INITIAL_STATE};

			
			
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
			
			//enum {PARSING = 1000, TYPE_CONFIRMATION, DATA_TRANSMISSION, NUMPACKETS};
			//enum {INACTIVE = 0, ACTIVE};
			//enum {VEC_MP = 4000, VEC_D, MAT_MP, MAT_D, COMP_MP, COMP_D, INDICES};
			
			
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
      br_exit(ERROR_MEMORY_ALLOCATION);
    }
    return x;
  }
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
#ifdef _HAVE_MPI
  MPI_Abort(MPI_COMM_WORLD, errorCode);
#else
  exit(errorCode);
#endif
}



void deliberate_segfault(){
	printf("the following segfault is deliberate\n");
	int *faulty = NULL;
	faulty[-10] = faulty[10]+faulty[0];
	
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
		print_comp_matlab(this->vertices[ii].projVal_mp,"proj");
		printf("type: %d\n", this->vertices[ii].type);
	}
}


int vertex_set::setup_vertices(boost::filesystem::path INfile)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices;
	int num_vars;
	fscanf(IN, "%d\n\n", &num_vertices);
	
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
		
		mpf_inp_str(temp_vertex.projVal_mp->r, IN, 10);
		mpf_inp_str(temp_vertex.projVal_mp->i, IN, 10);
		
		fscanf(IN,"%d\n",&temp_vertex.type);
		
		
		vertex_set::add_vertex(temp_vertex);
	}
	
	
	fclose(IN);
	
	if (this->num_vertices!=num_vertices) {
		printf("parity error in num_vertices.\n\texpected: %d\tactual: %d\n",this->num_vertices,num_vertices);
		exit(25943);
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
	fprintf(OUT,"%d\n\n",num_vertices);
	for (int ii = 0; ii < num_vertices; ii++)
	{ // output points
		fprintf(OUT,"%d\n", vertices[ii].pt_mp->size);
		for(int jj=0;jj<vertices[ii].pt_mp->size;jj++) {
			print_mp(OUT, 0, &vertices[ii].pt_mp->coord[jj]);
			fprintf(OUT,"\n");
		}
		print_mp(OUT, 0, vertices[ii].projVal_mp);
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

	print_point_to_screen_matlab(source_vertex.pt_mp,"new_added_point");
	
	this->counters[source_vertex.type]++;
	this->indices[source_vertex.type].push_back(current_index);

	return current_index;
}






int decomposition::index_in_vertices(vertex_set &V,
																		 vec_mp testpoint, comp_mp projection_value,
																		 tracker_config_t T)
{
	int ii;
	int index = -1;
		
//	V.print_to_screen();
	
	std::map< int , int >::iterator type_iter;
	
	for (type_iter = this->counters.begin(); type_iter!= this->counters.end(); type_iter++) {
		for (ii=0; ii<type_iter->second; ii++) {
			int current_index = this->indices[type_iter->first][ii];
			
			
			sub_mp(V.diff, projection_value, V.vertices[current_index].projVal_mp);
			mpf_abs_mp(V.abs, V.diff);
			
			if (mpf_cmp(V.abs, V.zerothresh) > 0){
				continue;
			}
			
			for (int jj=1; jj<V.num_natural_variables; jj++) {
				div_mp(&V.checker_1->coord[jj-1], &testpoint->coord[jj],  &testpoint->coord[0]);
				div_mp(&V.checker_2->coord[jj-1],&V.vertices[current_index].pt_mp->coord[jj], &V.vertices[current_index].pt_mp->coord[0]);
			}
			
			if (isSamePoint_inhomogeneous_input(V.checker_1, V.checker_2)){
//				std::cout << "these two points are the SAME:\n";
//				print_point_to_screen_matlab(V.vertices[current_index].pt_mp,"candidate");
//				print_point_to_screen_matlab(testpoint,"testpoint");
//				mypause();
				
				index = current_index;
				break;
			}
//			else{
//				std::cout << "these two points are DIFFERENT:\n";
//				print_point_to_screen_matlab(V.vertices[current_index].pt_mp,"candidate");
//				print_point_to_screen_matlab(testpoint,"testpoint");
//				mypause();
//			}
			
		}
		if (index!=-1) {
			break;
		}
	}

		
	
	return index;	
}



int decomposition::index_in_vertices_with_add(vertex_set &V,
																							vec_mp testpoint, comp_mp projection_value,
																							tracker_config_t T)
{
	int index = decomposition::index_in_vertices(V, testpoint, projection_value, T);
	
	if (index==-1) {
		std::cout << "adding vertex to set.  type: " << NEW << "\n\n";
		
		vertex temp_vertex;
		
		set_mp(temp_vertex.projVal_mp,  projection_value);
		vec_cp_mp(temp_vertex.pt_mp, testpoint);
		
		
		temp_vertex.type = NEW;
		
		index = decomposition::add_vertex(V, temp_vertex);
		
	}
	
	return index;

}












int decomposition::setup(boost::filesystem::path INfile,
												 boost::filesystem::path & inputName,
												 boost::filesystem::path directoryName)
//setup the vertex structure
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
	
	vec_mp tempvec; init_vec_mp(tempvec, this->num_variables);
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
	converter << tempstr;
	converter >> curr_num_patches;
	converter.clear(); converter.str("");
	
//	std::cout << "during setup, num_patches is to be " << curr_num_patches << std::endl;
//	mypause();
	
	vec_mp temp_patch; init_vec_mp2(temp_patch,1,1024); temp_patch->size = 1;
	for (int ii=0; ii<curr_num_patches; ii++) {
		
		getline(fin,tempstr); //std::cout << tempstr << std::endl;
		getline(fin,tempstr); //std::cout << tempstr << " <--- supposed to have the size in it" << std::endl;
		
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
	}
	
//	std::cout << "correct number patches: " << curr_num_patches << " vs memory-num_patches " << this->num_patches << std::endl;
//	mypause();
	clear_vec_mp(temp_patch);
	
	fin.close();
	
	
	
	return num_vertices;
}





/**
 Output curve overall info as follows:
 
 **/
void decomposition::print(boost::filesystem::path input_deflated_Name, boost::filesystem::path outputfile)
{
	int ii;
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	fprintf(OUT,"%s\n",input_deflated_Name.filename().c_str());
	
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
		std::cerr << "decomposition was short projections\n";
		deliberate_segfault();
	}
	for (ii=0; ii<num_curr_projections; ii++) {
		for(int jj=0;jj<pi_mp[ii]->size;jj++)
		{
			print_mp(OUT, 0, &pi_mp[ii]->coord[jj]);
			fprintf(OUT,"\n");
		}
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




int curve_decomposition::setup_edges(boost::filesystem::path INfile)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	
	fscanf(IN, "%d\n", &this->num_edges);
	int left, midpt, right;
	for(int ii=0;ii<this->num_edges;ii++) {
		fscanf(IN,"%d %d %d",&left, &midpt, &right); scanRestOfLine(IN);
		this->edges.push_back(edge(left, midpt, right));
	}
	
	fclose(IN);
	return this->num_edges;
}





void curve_decomposition::add_edge(edge new_edge)
{
	this->num_edges++;
	this->edges.push_back(new_edge);
	return;
}

void curve_decomposition::print_edges(boost::filesystem::path outputfile)
/**Output edge structure as follows:
 # variables
 # edges
 name of input file
 edge 1
 
 edge 2

 
 for each edge, output the following information:
 index to left vertex in vertices
 index to right vertex in vertices
 index to midpoint vertex in vertices
 
 **/
{
	int ii;
	FILE *OUT = safe_fopen_write(outputfile.c_str());
	
	// output the number of vertices
	fprintf(OUT,"%d\n\n",num_edges);
	
	for(ii=0;ii<num_edges;ii++)
		fprintf(OUT,"%d %d %d \n",edges[ii].left,
						edges[ii].midpt,
						edges[ii].right);
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
	
	vec_mp difference;  init_vec_mp(difference, left->size);difference->size = left->size;
	comp_mp temp; init_mp(temp);
	
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
	
	comp_mp denom; init_mp(denom);
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
	
	comp_mp denom; init_mp(denom);
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
	comp_mp temp; init_mp(temp);
	comp_mp temp2; init_mp(temp2);
	int ii;
	for (ii=0; ii<left->size; ++ii) {
		set_mp(temp2,result)
		mul_mp(temp,&left->coord[ii],&right->coord[ii]);
		add_mp(result,temp2,temp);
	}
	clear_mp(temp);clear_mp(temp2);
}


/**
 computes the projection value given a homogeneous input.
 
 double type
 */
void projection_value_homogeneous_input(comp_d result, vec_d input, vec_d projection){
	set_zero_d(result);
	comp_d temp;
	for (int ii=0; ii<projection->size; ii++) {
		mul_d(temp, &input->coord[ii], &projection->coord[ii]);
		add_d(result, result, temp);
	}
	set_d(temp, result);
	div_d(result, temp, &input->coord[0])
	return;
}

/**
 computes the projection value given a homogeneous input.
 
 mp type
*/
void projection_value_homogeneous_input(comp_mp result, vec_mp input, vec_mp projection){
	set_zero_mp(result);
	comp_mp temp; init_mp(temp);
	for (int ii=0; ii<projection->size; ii++) {
		mul_mp(temp, &input->coord[ii], &projection->coord[ii]);
		add_mp(result, result, temp);
	}
	set_mp(temp, result);
	div_mp(result, temp, &input->coord[0])
	return;
}



int isSamePoint_inhomogeneous_input(point_d left, point_d right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_inhom_d with disparate sized points.\n");
		std::cout << "left: " << left->size << "\t right: " << right->size << std::endl;
		deliberate_segfault();
//		exit(-287);
	}
	
	
	int indicator = isSamePoint(left,NULL,52,right,NULL,52,1e-6);
	
	
	return indicator;
}


int isSamePoint_inhomogeneous_input(point_mp left, point_mp right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_inhom_mp with disparate sized points.\n");
		std::cout << "left: " << left->size << "\t right: " << right->size << std::endl;
		deliberate_segfault();
//		exit(-287);
	}
	
	int indicator = isSamePoint(NULL,left,65,NULL,right,65,1e-6); // make the bertini library call
	

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
	
	int indicator = isSamePoint(dehom_left,NULL,52,dehom_right,NULL,52,1e-6);
	
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
	
	int indicator = isSamePoint(NULL,dehom_left,65,NULL,dehom_right,65,1e-6); // make the bertini library call
	
	clear_vec_mp(dehom_left); clear_vec_mp(dehom_right);
	
	return indicator;
}




void print_point_to_screen_matlab(vec_d M, std::string name)
{
	int kk;
	
	if (M->size==0) {
		printf("requested to print a vector '%s' which had size==0.  exiting\n",name.c_str());
	}
	
	printf("%s = [...\n",name.c_str());
	for (kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		printf(" %.15le+1i*%.15le;\n",M->coord[kk].r,M->coord[kk].i);
	}
	printf("];\n\n");
}

void print_point_to_screen_matlab(vec_mp M, std::string name)
{
	int kk;
	
	if (M->size==0) {
		printf("requested to print a vector '%s' which had size==0.  exiting\n",name.c_str());
	}
	
	printf("%s = [...\n",name.c_str());
	for (kk = 0; kk < M->size; kk++)
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
			printf("%.8le+1i*%.8le\t",M->entry[kk][jj].r,M->entry[kk][jj].i);
		}
		printf(";\n");
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
	mpf_out_str (NULL, 10, 6, M->r);
	printf("+1i*");
	mpf_out_str (NULL, 10, 6, M->i); // base 10, 6 digits
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

void cp_preproc_data(preproc_data *PPD, preproc_data PPD_input)
{
  int i, total_gp;
	
  PPD->num_funcs = PPD_input.num_funcs;
  PPD->num_hom_var_gp = PPD_input.num_hom_var_gp;
  PPD->num_var_gp = PPD_input.num_var_gp;
	
  total_gp = PPD->num_hom_var_gp + PPD->num_var_gp;
	
  PPD->size = (int *)bmalloc(total_gp * sizeof(int));
  PPD->type = (int *)bmalloc(total_gp * sizeof(int));
	
  for (i = 0; i < total_gp; i++)
  {
    PPD->size[i] = PPD_input.size[i];
    PPD->type[i] = PPD_input.type[i];
  }
	
  return;
}





void sort_increasing_by_real(vec_mp *projections_sorted, int **index_tracker, vec_mp projections_input){
	
	comp_mp large; init_mp(large);
	comp_d l; l->r = 1e9; l->i = 0;
	d_to_mp(large,l);
	
	
	vec_mp raw; init_vec_mp(raw,1);
	vec_cp_mp(raw,projections_input);
	
	change_size_vec_mp( (*projections_sorted), raw->size);
	(*projections_sorted)->size = raw->size;
	int ii,jj;
	//	comp_mp temp1, temp2; init_mp(temp1); init_mp(temp2);
	
	double min;
	double curr;
	int indicator = -1;
	for (ii=0; ii<raw->size; ii++) {
		min = 1e10;
		
		for (jj=0; jj<raw->size; jj++) {
			curr = mpf_get_d(raw->coord[jj].r);
			if ( curr < min) {
				indicator = jj;
				min = curr;
			}
		}
		if (indicator==-1) {
			printf("min projection value was *insanely* large\n");
			exit(1111);
		}
		
		
		
		(*index_tracker)[ii] = indicator;
		set_mp( &(*projections_sorted)->coord[ii],&raw->coord[indicator]);
		set_mp( &raw->coord[indicator],large);
	}
	return;
}



//input the raw number of variables including the homogeneous variable (of which there must be one)
// assume the array of integers 'randomized_degrees' is already initialized to the correct size.
void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, int ** randomized_degrees, int num_desired_rows, int num_funcs)
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
			(*randomized_degrees)[ii] = degrees[ii];
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
		(*randomized_degrees)[ii] = current_degree;
		
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



void send_patch_mp(patch_eval_data_mp * patch)
{
	char *patchStr = NULL;
	patch_eval_data_mp_int PED_int;
	MPI_Datatype mpi_patch_int;
	
	// setup mpi_patch_int
	create_patch_eval_data_mp_int(&mpi_patch_int);
	// setup PED_int
	cp_patch_mp_int(&PED_int, &patch, &patchStr, 0, 0);
	
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
	patchStr = (char *)bmalloc(PED_int.totalLength * sizeof(char));
	// recv patchStr
	MPI_Bcast(patchStr, PED_int.totalLength, MPI_CHAR, 0, MPI_COMM_WORLD);
	
	// setup _mp patch
	cp_patch_mp_int(patch, &PED_int, &patchStr, 1, 1);
	
	//	// setup _d patch
//	for (i = 0; i < BED->patch.patchCoeff->rows; i++)
//		for (j = 0; j < BED->patch.patchCoeff->cols; j++)
//		{
//			mp_to_d(&patch.patchCoeff->entry[i][j], &BED->BED_mp->patch.patchCoeff->entry[i][j]);
//		}
	
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
	patch_coeff = (comp_d *)bmalloc(PED_int.patchCoeff_rows * PED_int.patchCoeff_cols * sizeof(comp_d));
	MPI_Bcast(patch_coeff, PED_int.patchCoeff_rows * PED_int.patchCoeff_cols, mpi_comp_d, 0, MPI_COMM_WORLD);
	
	
// setup patch
	cp_patch_d_int(patch, &PED_int, &patch_coeff, 1);  // patch_coeff is freed in here
	
	// free mpi_comp_d & mpi_patch_d_int
	MPI_Type_free(&mpi_comp_d);
	MPI_Type_free(&mpi_patch_d_int);
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
	bstr = (char *)bmalloc(b_int.totalLength * sizeof(char));
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
	
	entries = (comp_d *)bmalloc(b_int.size * sizeof(comp_d));
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
//	bstr = (char *)bmalloc(b_int.totalLength * sizeof(char));
//	MPI_Recv(bstr, b_int.totalLength, MPI_CHAR, source,  VEC_MP, MPI_COMM_WORLD, &statty_mc_gatty);
//	
//	// setup b and clear bstr
//	cp_point_d_int(b, &b_int, &bstr, 1, 1, 1);
//	
//  // clear mpi_vec_d_int
//  MPI_Type_free(&mpi_vec_d_int);
	
  return;
}




