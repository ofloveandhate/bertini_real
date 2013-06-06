#include "data_type.h"

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
      printf("ERROR: malloc was unable to allocate memory (%d)!\n", (int) size);
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

	int faulty[2];
	printf("%d\n",faulty[-1]);
	
}


void init_vertex_set(vertex_set *V){
	V->vertices = NULL;  //Isolated real points.
	V->num_vertices = 0;
}

int add_vertex(vertex_set *V, vertex source_vertex){
	
//	printf("adding vertex\n");
	
	if ( (V->num_vertices==0) && (V->vertices!=NULL) ) {
		printf("initalization error in add_vertex\n");
		printf("V->num_vertices==%d && V->vertices!=NULL\n",V->num_vertices);
		exit(-1);
	}
	
	
	if
		( (V->num_vertices!=0) && (V->vertices==NULL) ) {
		printf("initalization error in add_vertex\n");
		printf("NULL vertices with positive num_vertices\n");
		exit(-1);
	}
	
	
	
	
	V->num_vertices++;
	
	//	printf("\n**************\nadding vertex:\n");
	//	print_point_to_screen_matlab_mp(source_vertex.pt_mp,"adding");
	//	printf("\n---------------------\n");
	
	
	if (V->num_vertices==1) {
		V->vertices = (vertex *)bmalloc( sizeof(vertex));
	}
	else{
		V->vertices = (vertex *)brealloc(V->vertices, V->num_vertices*sizeof(vertex));
	}
	
	int current_index = V->num_vertices-1;
	
	init_vertex(&V->vertices[current_index]);
	cp_vertex(&V->vertices[current_index],source_vertex);
	
		
	return current_index;
}




int curve_add_vertex(curveDecomp_d *C, vertex_set *V, vertex source_vertex){
	
	
	int current_index = add_vertex(V, source_vertex);
	
	
	switch (source_vertex.type) { // if you add a type to the enum, ensure it gets a case here!
		case CRITICAL:
			C->num_V1++; // increment number of points of this type
			
			if (C->num_V1==1)
				C->V1_indices = (int *)bmalloc(C->num_V1*sizeof(int));
			else
				C->V1_indices = (int *)brealloc(C->V1_indices, (C->num_V1)*sizeof(int));
			
			C->V1_indices[C->num_V1-1] = current_index; // the number in C having already been incremented, we -1;
			
			break;
			
		case NEW:
			C->num_new++; // increment number of points of this type
			C->num_V1++;  // increment number of points of this type
			
			
			if (C->num_V1==1)
				C->V1_indices = (int *)bmalloc(C->num_V1*sizeof(int));
			else
				C->V1_indices = (int *)brealloc(C->V1_indices, (C->num_V1)*sizeof(int));
			
			C->V1_indices[C->num_V1-1] = current_index; // the number in C having already been incremented, we -1;
			
			
			if (C->num_new==1)
				C->new_indices = (int *)bmalloc(C->num_new*sizeof(int));
			else
				C->new_indices = (int *)brealloc(C->new_indices, (C->num_new)*sizeof(int));
			
			C->new_indices[C->num_new-1] = current_index; // the number in C having already been incremented, we -1;
			
			break;
			
		case ISOLATED:
			C->num_V0++; // increment number of points of this type
			
			if (C->num_V0==1)
				C->V0_indices = (int *)bmalloc(C->num_V0*sizeof(int)); // if the first, allocate
			else
				C->V0_indices = (int *)brealloc(C->V0_indices, (C->num_V0)*sizeof(int));  // reallocate
			
			C->V0_indices[C->num_V0-1] = current_index; // the number in C having already been incremented, we -1;
			
			
			break;
			
		case MIDPOINT:
			C->num_midpts++; // increment number of points of this type
			
			
			if (C->num_midpts==1)
				C->midpt_indices = (int *)br_malloc(C->num_midpts*sizeof(int)); // if the first, allocate
			else
				C->midpt_indices = (int *)brealloc(C->midpt_indices, (C->num_midpts)*sizeof(int));  // reallocate
			
			C->midpt_indices[C->num_midpts-1] = current_index; // the number in C having already been incremented, we -1;
			
			break;
			
		case SAMPLE_POINT:
			C->num_samples++; // increment number of points of this type
			
			
			if (C->num_samples==1)
				C->sample_indices = (int *)br_malloc(sizeof(int)); // if the first, allocate
			else
				C->sample_indices = (int *)brealloc(C->sample_indices, (C->num_samples)*sizeof(int));  // reallocate
			
			C->sample_indices[C->num_samples-1] = current_index; // the number in C having already been incremented, we -1;
			
			break;
			
		default:
			
			printf("attempting to add a vertex of unset or invalid type (%d)\n",source_vertex.type);
			exit(-220);
			break;
	}

	return current_index;
}






int curve_index_in_vertices(curveDecomp_d *C, vertex_set *V,
														vec_mp testpoint, comp_mp projection_value,
														tracker_config_t T)
{
	int ii;
	
	
	int index = -1;
	
	
	for (ii=0; ii<C->num_V1; ii++) {
		
		int current_index = C->V1_indices[ii];

		if (isSamePoint_homogeneous_input_mp(V->vertices[current_index].pt_mp, testpoint)){
			index = ii;
			break;
		}
		
	}
	
	for (ii=0; ii<C->num_V0; ii++) {
		
		int current_index = C->V0_indices[ii];

		if (isSamePoint_homogeneous_input_mp(V->vertices[current_index].pt_mp, testpoint)){
			index = ii;
			break;
		}
		
	}
	
	
	for (ii=0; ii<C->num_isolated; ii++) {
		
		int current_index = C->isolated_indices[ii];
		if (isSamePoint_homogeneous_input_mp(V->vertices[current_index].pt_mp, testpoint)){
			index = ii;
			break;
		}
		
	}
	
	
	for (ii=0; ii<C->num_new; ii++) {
		
		int current_index = C->new_indices[ii];
		if (isSamePoint_homogeneous_input_mp(V->vertices[current_index].pt_mp, testpoint)){
			index = ii;
			break;
		}
	}
	

	return index;	
}



int curve_index_in_vertices_with_add(curveDecomp_d *C, vertex_set *V,
																		 vec_mp testpoint, comp_mp projection_value,
																		 tracker_config_t T)
{
	int index = curve_index_in_vertices(C, V, testpoint, projection_value, T);
	
	if (index==-1) {
		vertex temp_vertex; init_vertex(&temp_vertex);
		
		set_mp(temp_vertex.projVal_mp,  projection_value);
		vec_cp_mp(temp_vertex.pt_mp, testpoint);
		
		
		temp_vertex.type = NEW;
		index = add_vertex(V,temp_vertex);
	}
	
	return index;

}

void vertex_set_print_to_screen(vertex_set *V){
	int ii;
	printf("vertex set has %d vertices:\n\n",V->num_vertices);
	for (ii=0; ii<V->num_vertices; ++ii) {
		print_point_to_screen_matlab_mp(V->vertices[ii].pt_mp,"vert");
		print_comp_mp_matlab(V->vertices[ii].projVal_mp,"proj");
		printf("type: %d\n", V->vertices[ii].type);
	}
	
}




/**
 copy in a vertex
 */
void cp_vertex(vertex *target_vertex, vertex source_vertex){

	vec_cp_mp(target_vertex->pt_mp,source_vertex.pt_mp);
	set_mp(target_vertex->projVal_mp, source_vertex.projVal_mp);

	target_vertex->type = source_vertex.type;
	
	return;
}



void init_vertex(vertex *curr_vertex){
	
	init_vec_mp(curr_vertex->pt_mp,1); curr_vertex->pt_mp->size = 1; // this is the line giving us trouble.
	init_mp(curr_vertex->projVal_mp);
	
	curr_vertex->type = -100;
	return;
}


void init_edge(edge *curr_edge){
	
	curr_edge->left = curr_edge->right = curr_edge->midpt = -1; // initialize to impossible value.
	return;
}



void cp_edge(edge *new_edge, edge source_edge){
	
	new_edge->left = source_edge.left;
	new_edge->right= source_edge.right;
	new_edge->midpt = source_edge.midpt;
	return;
	
}

void add_edge(curveDecomp_d *C, edge new_edge){
	
	if ( ( C->num_edges==0 && C->edges!=NULL ) || ( C->num_edges!=0 && C->edges==NULL ) ) {
		printf("intialization error in curve decomposition\n");
		exit(-1);
	}
	
	C->num_edges++;
	
	if (C->num_edges==1) {
		C->edges = (edge *)bmalloc(C->num_edges*sizeof(edge));
	}
	else{
		C->edges = (edge *)brealloc(C->edges, C->num_edges*sizeof(edge));
	}
	
	
	init_edge(&C->edges[C->num_edges-1]); // intialize the new edge

	cp_edge(&C->edges[C->num_edges-1], new_edge);
	
	return;
}




void init_curveDecomp_d(curveDecomp_d *C){
	
//	C->MPType = MPType;
	C->num_V0 = C->num_V1 = C->num_new = C->num_midpts = C->num_isolated = 0;
	
	C->num_edges = 0;
	
	C->edges							= NULL;
	C->V0_indices					= NULL;
	C->V1_indices					= NULL;
	C->midpt_indices			= NULL;
	C->new_indices				= NULL;
	C->isolated_indices		= NULL;
	
	
	return;
}



void clear_vertex_set_d(vertex_set *V){
	int ii;
	for(ii=0;ii<V->num_vertices;ii++)
	{
//		if(C->MPType == 0)
//		{
//			clear_vec_d(V->vertices[ii].pt_d);
//		}
//		else
//		{
			clear_vec_mp(V->vertices[ii].pt_mp);
//		}
	}
}
void clear_curveDecomp_d(curveDecomp_d *C){
	
	free(C->edges);
	
	
	free(C->isolated_indices);
	free(C->V0_indices);
	free(C->V1_indices);
	free(C->new_indices);
}



void clear_sample_d(sample_data *S, int MPType){
	
	int ii;
	for(ii=0;ii<S->num_edges;ii++)
	{
		free(S->sample_indices[ii]);
	}
	free(S->sample_indices);
}



void norm_of_difference(mpf_t result, vec_mp left, vec_mp right){
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


void dehomogenize_d(vec_d *result, vec_d dehom_me){
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


void dot_product_d(comp_d result, vec_d left, vec_d right){
	if (left->size!=right->size) {
		printf("attempting to dot_d two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		raise(-78);
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
		raise(-79);
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


void projection_value_homogeneous_input_d(comp_d result, vec_d input, vec_d projection){
	
	if (projection->size != input->size) {
		printf("computing projection values of incompatibly sized vectors\n");
		exit(1024);
	}
	
	vec_d dehom;
	init_vec_d(dehom,input->size-1); dehom->size = input->size-1;
	
	vec_d pi_no_hom;
	init_vec_d(pi_no_hom,input->size-1); pi_no_hom->size = input->size-1;
	int ii;
	for (ii=0; ii<input->size-1; ii++) {
		set_d(&pi_no_hom->coord[ii],&projection->coord[ii+1]);
	}
	
	dehomogenize_d(&dehom,input);
	dot_product_d(result,dehom,pi_no_hom); // set projection value
	
	clear_vec_d(dehom); clear_vec_d(pi_no_hom);
	
	return;
}

// MP version
void projection_value_homogeneous_input(comp_mp result, vec_mp input, vec_mp projection){
	
	if (projection->size != input->size) {
		printf("computing projection values of incompatibly sized vectors\n");
		exit(1024);
	}
	
	vec_mp dehom;
	init_vec_mp(dehom,input->size-1); dehom->size = input->size-1;
	
	vec_mp pi_no_hom;
	init_vec_mp(pi_no_hom,input->size-1); pi_no_hom->size = input->size-1;
	int ii;
	for (ii=0; ii<input->size-1; ii++) {
		set_mp(&pi_no_hom->coord[ii],&projection->coord[ii+1]);
	}
	
	dehomogenize_mp(&dehom,input);
	dot_product_mp(result,dehom,pi_no_hom); // set projection value
	
	clear_vec_mp(dehom); clear_vec_mp(pi_no_hom);
	return;
}



int isSamePoint_inhomogeneous_input_d(point_d left, point_d right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_hom_d with disparate sized points.\n");
		exit(-287);
	}
	
	
	int indicator = isSamePoint(left,NULL,52,right,NULL,52,1e-6);
	
	
	return indicator;
}


int isSamePoint_inhomogeneous_input_mp(point_mp left, point_mp right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_hom_mp with disparate sized points.\n");
		exit(-287);
	}
	
	int indicator = isSamePoint(NULL,left,65,NULL,right,65,1e-6); // make the bertini library call
	

	return indicator;
}



int isSamePoint_homogeneous_input_d(point_d left, point_d right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_hom_d with disparate sized points.\n");
		exit(-287);
	}
	
	vec_d dehom_left;  init_vec_d(dehom_left,left->size-1);  dehom_left->size = left->size-1;
	vec_d dehom_right; init_vec_d(dehom_right,right->size-1); dehom_right->size = right->size-1;
	
	dehomogenize_d(&dehom_left,left);
	dehomogenize_d(&dehom_right,right);
	
	int indicator = isSamePoint(dehom_left,NULL,52,dehom_right,NULL,52,1e-6);
	
	clear_vec_d(dehom_left); clear_vec_d(dehom_right);
	
	return indicator;
}
int isSamePoint_homogeneous_input_mp(point_mp left, point_mp right){
	
	if (left->size!=right->size) {
		printf("attempting to isSamePoint_hom_mp with disparate sized points.\n");
		exit(-287);
	}
	
	vec_mp dehom_left;  init_vec_mp(dehom_left,left->size-1);  dehom_left->size = left->size-1;
	vec_mp dehom_right; init_vec_mp(dehom_right,right->size-1); dehom_right->size = right->size-1;
	
	dehomogenize_mp(&dehom_left,left);
	dehomogenize_mp(&dehom_right,right);
	
	int indicator = isSamePoint(NULL,dehom_left,65,NULL,dehom_right,65,1e-6); // make the bertini library call
	
	clear_vec_mp(dehom_left); clear_vec_mp(dehom_right);
	
	return indicator;
}




void print_point_to_screen_matlab(vec_d M, char name[])
{
	int kk;
	
	if (M->size==0) {
		printf("requested to print a vector '%s' which had size==0.  exiting\n",name);
//		exit(-1);
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
//		exit(-1);
	}
	
	printf("%s = [...\n",name);
	for (kk = 0; kk < M->size; kk++)
	{ // print kth coordinate
		mpf_out_str (NULL, 10, 16, M->coord[kk].r);
		printf("+1i*");
		mpf_out_str (NULL, 10, 16, M->coord[kk].i);
		printf(";\n");
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
			printf(" %.7le+1i*%.7le",M->entry[kk][jj].r,M->entry[kk][jj].i);
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


