#include "data_type.h"


int add_vertex(curveDecomp_d *C, vertex source_vertex){
	
	printf("adding vertex\n");
	
	if ( ( (C->num_vertices==0) && (C->vertices!=NULL) ) || ( (C->num_vertices!=0) && (C->vertices==NULL) ) ) {
		printf("intialization error in add_vertex\n");
		printf("C->num_vertices = %d\n",C->num_vertices);
		exit(-1);
	}
	
	C->num_vertices++;
	
	//	printf("\n**************\nadding vertex:\n");
	//	print_point_to_screen_matlab_mp(source_vertex.pt_mp,"adding");
	//	printf("\n---------------------\n");
	
	
	if (C->num_vertices==1) {
		C->vertices = (vertex *)bmalloc( C->num_vertices*sizeof(vertex));
	}
	else{
		C->vertices = (vertex *)brealloc(C->vertices, C->num_vertices*sizeof(vertex));
	}
	
	int current_index = C->num_vertices-1;
	
	init_vertex_mp(&C->vertices[current_index]);
	
	cp_vertex_mp(&C->vertices[current_index],source_vertex);
	
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
				C->midpt_indices = (int *)bmalloc(C->num_midpts*sizeof(int)); // if the first, allocate
			else
				C->midpt_indices = (int *)brealloc(C->midpt_indices, (C->num_midpts)*sizeof(int));  // reallocate
			
			C->midpt_indices[C->num_midpts-1] = current_index; // the number in C having already been incremented, we -1;
			
			break;
			
		default:
			
			printf("attempting to add a vertex of unset or invalid type (%d)\n",source_vertex.type);
			exit(-220);
			break;
	}
	
	return current_index;
}






int index_in_vertices(curveDecomp_d *C, vec_mp testpoint, comp_mp projection_value, tracker_config_t T, int sidedness){
	int ii;
	
	
	//	sidedness = -1 if from left, 1 if from right.
	int index = -1;
	
	
	
	
	for (ii=0; ii<C->num_vertices; ii++) {
		if (isSamePoint_homogeneous_input_mp(C->vertices[ii].pt_mp, testpoint)){
//			printf("these two points are the same when dehomogenized\n");
//			print_point_to_screen_matlab_mp(testpoint,"testpoint");
//			print_point_to_screen_matlab_mp(C->vertices[ii].pt_mp,"candidate");
			index = ii;
			break;
		}
	}
	
	
	
	if (index==-1) {
//		printf("vertex not found; adding to vertex set\n");
		vertex temp_vertex;  init_vertex_mp(&temp_vertex);
		set_mp(temp_vertex.projVal_mp,  projection_value);
		vec_cp_mp(temp_vertex.pt_mp, testpoint);
		temp_vertex.type = NEW;
		add_vertex(C,temp_vertex);
		index = C->num_vertices-1;
	}
	
	switch (sidedness) {
		case -1:
			C->vertices[index].num_left++; // increase the number of times this vertex has been a left vertex
			break;
		
		case 1:
			C->vertices[index].num_right++;// increase the number of times this vertex has been a right vertex
			break;
		
		case 0: // midpoint

			break;
			
		default:
			break;
	}
	
	
	return index;	
}



/**
 copy in a vertex
 */
void cp_vertex_mp(vertex *target_vertex, vertex source_vertex){
//assume the vertex is initialized
//	printf("copying vertex\n");
	vec_cp_mp(target_vertex->pt_mp,source_vertex.pt_mp);
	set_mp(target_vertex->projVal_mp, source_vertex.projVal_mp);
	target_vertex->num_left = source_vertex.num_left; target_vertex->num_right = source_vertex.num_right;
	target_vertex->type = source_vertex.type;
	
	return;
}



void init_vertex_d(vertex *curr_vertex){
	
	init_vec_d(curr_vertex->pt,1); curr_vertex->pt->size = 1; // this is the line giving us trouble.
	
	curr_vertex->type = -100;
	curr_vertex->num_left = 0;
	curr_vertex->num_right = 0;
	return;
}


void init_vertex_mp(vertex *curr_vertex){
	
	init_vec_mp(curr_vertex->pt_mp,1); curr_vertex->pt_mp->size = 1; // this is the line giving us trouble.
	init_mp(curr_vertex->projVal_mp);
	
	curr_vertex->type = -100;
	curr_vertex->num_left = 0;
	curr_vertex->num_right = 0;
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
	
	C->num_V0 = C->num_V1 = C->num_new = C->num_midpts = 0;
	
	C->num_edges = C->num_vertices = 0;
	
	C->vertices=NULL;
	C->edges=NULL;
	
	C->V0_indices = NULL;
	C->V1_indices = NULL;
	C->midpt_indices = NULL;
	C->new_indices = NULL;
	
	
	return;
}


void clear_curveDecomp_d(curveDecomp_d *C, int MPType){
	
	int ii;
	for(ii=0;ii<C->num_vertices;ii++)
	{
		if(MPType == 0)
		{
			clear_vec_d(C->vertices[ii].pt);
		}
		else
		{
			clear_vec_mp(C->vertices[ii].pt_mp);
		}
	}
	free(C->vertices);

	free(C->edges);
	
}



void clear_sample_d(sample_d *S, int MPType){
	
	int ii, jj;
	for(ii=0;ii<S->num_edges;ii++)
	{
		if(MPType == 0)
		{
			clear_vec_d(S->projection_values[ii]);
		}
		else
		{
			clear_vec_mp(S->projection_values_mp[ii]);
		}
		for(jj=0; jj<S->num_pts[ii];jj++)
		{
			if(MPType == 0)
			{
		    clear_vec_d(S->vertices[ii][jj]);
			}
			else
			{
				clear_vec_mp(S->vertices_mp[ii][jj]);
			}
		}
		free(S->refine[ii]);
	}
	if(MPType == 0)
	{
		free(S->projection_values);
		free(S->vertices);
	}
	else
	{
		free(S->projection_values_mp);
		free(S->vertices_mp);
	}
	free(S->refine);
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
		printf("attempting to dot two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		exit(-78);
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
		printf("attempting to dot two vectors not of the same size! (%d!=%d)\n",left->size,right->size);
		exit(-78);
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


