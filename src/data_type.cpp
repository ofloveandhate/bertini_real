#include "data_type.hpp"

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
	fscanf(IN, "%d %d\n\n", &num_vertices, &num_vars);
	
	vertex temp_vertex;
	change_size_vec_mp(temp_vertex.pt_mp,num_vars); temp_vertex.pt_mp->size = num_vars;
	
	for(int ii=0;ii<num_vertices;ii++)
	{
		for(int jj=0;jj<num_vars;jj++)
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
 
 .
 .
 .
 **/
void vertex_set::print(boost::filesystem::path outputfile)
{
	
	FILE *OUT = safe_fopen_write(outputfile);
	
	// output the number of vertices
	fprintf(OUT,"%d %d\n\n",num_vertices, vertices[0].pt_mp->size);
	for (int ii = 0; ii < num_vertices; ii++)
	{ // output points
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





int decomposition::add_vertex(vertex_set & V, vertex source_vertex){
	
	
	int current_index = V.add_vertex(source_vertex);
	
	std::cout << "vertex type to add: " << source_vertex.type << std::endl;
	
	
	if (this->counters.find(source_vertex.type) == this->counters.end()) {
		this->counters[source_vertex.type] = 0;
	}
	
		
	std::cout << "counters[vertex_type]: " << this->counters[source_vertex.type] << std::endl;
	
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
//			std::cout << "current_index: " << current_index << " type: " << type_iter->first << " " << ii << std::endl;
			if (isSamePoint_homogeneous_input(V.vertices[current_index].pt_mp, testpoint)){
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
	
	for (ii=0; ii<dimension; ii++) {
		for(int jj=0;jj<num_variables;jj++)
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
 .
 .
 .
 
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
		printf("attempting to dot two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		exit(-78);
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



void dot_product_d(comp_d result, vec_d left, vec_d right){
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

void dot_product_mp(comp_mp result, vec_mp left, vec_mp right){
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


