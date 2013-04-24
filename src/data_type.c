#include "data_type.h"



void init_configuration(program_configuration *options) {
	
	
	options->user_projection = 0;
	options->projection_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->projection_filename = "";
	
	options->user_randomization = 0;
	options->randomization_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->randomization_filename = "";
	
	options->input_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->input_filename = "input\0";
	
	options->witness_set_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	options->witness_set_filename = "witness_set\0";
	
	options->input_deflated_filename = (char *)bmalloc(MAX_STRLEN*sizeof(char));
	//this must be set after deflation is run.
	return;
}

void clear_configuration(program_configuration *options) {
	
	free(options->projection_filename);
	free(options->randomization_filename);
	free(options->input_filename);
	free(options->input_deflated_filename);
	
	return;
}




void add_vertex_to_V0(curveDecomp_d *C, vertex_d new_vertex){
	
	if ( ( C->num_V0==0 && C->V0!=NULL ) || ( C->num_V0!=0 && C->V0==NULL ) ) {
		printf("intialization error in add_vertex_to_V0\n");
		exit(-1);
	}
	
	C->num_V0++;
	
	if (C->num_V0==1) {
		C->V0 = (vertex_d *)bmalloc( C->num_V0*sizeof(vertex_d));
	}
	else{
		C->V0 = (vertex_d *)brealloc(C->V0, C->num_V0*sizeof(vertex_d));
	}
	
	init_vertex_mp(&C->V0[C->num_V0-1]);
	cp_vertex_mp(&C->V0[C->num_V0-1],new_vertex);
	
	return;
}



void add_vertex_to_V1(curveDecomp_d *C, vertex_d new_vertex){
	
	if ( ( (C->num_V1==0) && (C->V1!=NULL) ) || ( (C->num_V1!=0) && (C->V1==NULL) ) ) {
		printf("intialization error in add_vertex_to_V1\n");
		printf("C->num_V1 = %d\n",C->num_V1);
		exit(-1);
	}
	
	C->num_V1++;
	
	if (C->num_V1==1) {
		printf("initializing V1\n");
		C->V1 = (vertex_d *)bmalloc( C->num_V1*sizeof(vertex_d));
	}
	else{
		printf("adding point to already initialized V1\n");
		printf("%d\n",C->num_V1);
		C->V1 = (vertex_d *)brealloc(C->V1, C->num_V1*sizeof(vertex_d));
		printf("after realloc\n");
	}
	
	init_vertex_mp(&C->V1[C->num_V1-1]);
	printf("after init\n");
	cp_vertex_mp(&C->V1[C->num_V1-1],new_vertex);
	
	return;
}



int index_in_V1(curveDecomp_d *C, vec_mp testpoint, comp_mp projection_value, tracker_config_t T, int sidedness){
	int ii;
	
	
	//	sidedness = -1 if from left, 1 if from right.
	int index = -1;
	
	
	
	
	for (ii=0; ii<C->num_V1; ii++) {
		if (isSamePoint(NULL,C->V1[ii].pt_mp,65,NULL,testpoint,65,1e-6)){
			print_point_to_screen_matlab_mp(testpoint,"testpoint");
			print_point_to_screen_matlab_mp(C->V1[ii].pt_mp,"candidate");
			index = ii;
			break;
		}
	}
	
	
	
	
	
	if (index==-1) {
		printf("vertex not found; adding to V1\n");
		vertex_d temp_vertex;  init_vertex_mp(&temp_vertex);
		set_mp(temp_vertex.projVal_mp,  projection_value);
		vec_cp_mp(temp_vertex.pt_mp, testpoint);
		
		add_vertex_to_V1(C,temp_vertex);
		index = C->num_V1-1;
	}
	
	if (sidedness==-1) {
		C->V1[index].num_left++;
	}
	else{
		C->V1[index].num_right++;
	}
	
	//	printf("index = %d\n",index);
	//	mypause();
	
	return index;
	//	typedef struct
	//	{
	//		vertex_d *V0;  //Isolated real points.
	//		vertex_d *V1;  //Critical points AND new non-critical endpoints of edges.
	//		//  vertex_d *midPts;  //Midpoints of edges.
	//		edge_d *edges;
	//		int      num_V0;
	//		int      num_V1;
	//		//  int      num_midPts;
	//		int      num_edges;
	//	}curveDecomp_d;
	
}



/**
 copy in a vertex
 */
void cp_vertex_mp(vertex_d *target_vertex, vertex_d new_vertex){
//assume the vertex is initialized
//	printf("copying vertex\n");
	vec_cp_mp(target_vertex->pt_mp,new_vertex.pt_mp);
	set_mp(target_vertex->projVal_mp, new_vertex.projVal_mp);
	target_vertex->num_left = new_vertex.num_left; target_vertex->num_right = new_vertex.num_right;
	target_vertex->type = new_vertex.type;
	
	return;
}




void init_vertex_mp(vertex_d *curr_vertex){
	
	printf("inside init_vertex_mp\n");
	init_vec_mp(curr_vertex->pt_mp,1); curr_vertex->pt_mp->size = 1; // this is the line giving us trouble.
	printf("here\n");
	init_mp(curr_vertex->projVal_mp);
	
	curr_vertex->type = -1;
	curr_vertex->num_left = 0;
	curr_vertex->num_right = 0;
	printf("exiting init_vertex_mp\n");
	return;
}


void init_edge_mp(edge_d *curr_edge, int num_variables){
	curr_edge->left = curr_edge->right = -1; // initialize to impossible value.
	init_vec_mp(curr_edge->midpt_mp,num_variables); curr_edge->midpt_mp->size = num_variables;
	init_vec_mp(curr_edge->pi_mp,num_variables); curr_edge->pi_mp->size = num_variables;
	return;
}

void init_edge_d(edge_d *curr_edge){
	

	curr_edge->left = curr_edge->right = -1; // initialize to impossible value.
	
	init_vec_d(curr_edge->midpt,1); curr_edge->midpt->size = 1;
	init_vec_d(curr_edge->pi,1); curr_edge->pi->size = 1;
	
	return;
}



void add_edge_mp(curveDecomp_d *C, edge_d new_edge){
	
	if ( ( C->num_edges==0 && C->edges!=NULL ) || ( C->num_edges!=0 && C->edges==NULL ) ) {
		printf("intialization error in curve decomposition\n");
		exit(-1);
	}
	
	if (C->num_edges==0) {
		C->edges = (edge_d *)bmalloc(1*sizeof(edge_d));
	}
	else{
		C->edges = (edge_d *)brealloc(C->edges, C->num_edges+1*sizeof(edge_d));
	}
	
	C->num_edges++;
	init_edge_mp(&C->edges[C->num_edges-1],new_edge.midpt_mp->size);
//	C->edges[C->num_edges-1].size = new_edge.midpt_mp->size;
	vec_cp_mp(C->edges[C->num_edges-1].midpt_mp,new_edge.midpt_mp);
	vec_cp_mp(C->edges[C->num_edges-1].pi_mp,new_edge.pi_mp);
	C->edges[C->num_edges-1].left = new_edge.left;
	C->edges[C->num_edges-1].right = new_edge.right;
	return;
}




void init_curveDecomp_d(curveDecomp_d *C){
	
	C->num_V0=0;
	C->num_V1=0;
	C->num_edges=0;
	
	C->V0=NULL;
	C->V1=NULL;
	
	C->edges=NULL;
	return;
}


void clear_curveDecomp_d(curveDecomp_d *C, int MPType){
	
	int ii;
	for(ii=0;ii<C->num_V0;ii++)
	{
		if(MPType == 0)
		{
			clear_vec_d(C->V0[ii].pt);
		}
		else
		{
			clear_vec_mp(C->V0[ii].pt_mp);
		}
	}
	free(C->V0);
	for(ii=0;ii<C->num_V1;ii++)
	{
		if(MPType == 0)
		{
			clear_vec_d(C->V1[ii].pt);
		}
		else
		{
			clear_vec_mp(C->V1[ii].pt_mp);
		}
	}
	free(C->V1);
	for(ii=0;ii<C->num_edges;ii++)
	{
		if(MPType == 0)
		{
			clear_vec_d(C->edges[ii].midpt);
			clear_vec_d(C->edges[ii].pi);
		}
		else
		{
			clear_vec_mp(C->edges[ii].midpt_mp);
			clear_vec_mp(C->edges[ii].pi_mp);
		}
//		clear_witness_set(C->edges[ii].W);
	}
	free(C->edges);
	
}



void clear_sample_d(sample_d *S, int MPType){
	
	int ii, jj;
	for(ii=0;ii<S->num_edges;ii++)
	{
		if(MPType == 0)
		{
			clear_vec_d(S->proj_vertices[ii]);
		}
		else
		{
			clear_vec_mp(S->proj_vertices_mp[ii]);
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
		free(S->proj_vertices);
		free(S->vertices);
	}
	else
	{
		free(S->proj_vertices_mp);
		free(S->vertices_mp);
	}
	free(S->refine);
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
	W_out->W.num_pts  = W_left.W.num_pts+W_right.W.num_pts;
	W_out->W_mp.num_pts = W_left.W.num_pts+W_right.W.num_pts;
	

  W_out->W.pts=(point_d *)bmalloc(W_out->W.num_pts*sizeof(point_d));
	W_out->W_mp.pts=(point_mp *)bmalloc(W_out->W.num_pts*sizeof(point_mp));
	
	counter = 0;
	for (ii=0; ii<W_left.W.num_pts; ii++) {
		init_vec_d(W_out->W.pts[counter],W_out->num_variables); W_out->W.pts[counter]->size = W_out->num_variables;
		vec_cp_d(W_out->W.pts[counter],W_left.W.pts[ii]);
		
		init_vec_mp(W_out->W_mp.pts[counter],W_out->num_variables); W_out->W_mp.pts[counter]->size = W_out->num_variables;
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


void init_variable_names(witness_set_d *W, int num_vars){
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
void cp_names(witness_set_d *W_out, witness_set_d W_in){
	int ii;
	
	if (W_in.num_variables==0 || W_in.variable_names==NULL) {
		printf("attempting to copy variable names from witness_set with no variables");
		exit(1333);
	}
	
	if (W_out->variable_names==NULL) {
		init_variable_names(W_out,W_in.num_variables);
	}
	
	for (ii=0; ii<W_in.num_variables; ++ii) {
		strcpy(W_out->variable_names[ii],  W_in.variable_names[ii]);
	}
	
//	printf("copied these variable names:\n");
//	for (ii=0; ii<W_in.num_variables; ++ii) {
//		printf("%s\n",W_out->variable_names[ii]);
//	}
//	mypause();
}




//copies the mp and d linears from in to out.
void cp_linears(witness_set_d *W_out, witness_set_d W_in){
	int ii;
	
	
	W_out->num_linears = W_in.num_linears;
	printf("%d linears to copy.\n",W_in.num_linears);


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

void cp_witness_set(witness_set_d *W_out, witness_set_d W_in){
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
	
	
	
	
	W_out->L = (vec_d *)bmalloc(W_in.num_linears*sizeof(vec_d));
	W_out->L_mp = (vec_mp *)bmalloc(W_in.num_linears*sizeof(vec_mp));
	// merge the left and right linears into the output.
	for (ii=0; ii<W_in.num_linears; ii++) {
		init_vec_d(W_out->L[ii],W_out->num_variables); W_out->L[ii]->size = W_out->num_variables;
		vec_cp_d(W_out->L[ii],W_in.L[ii]);
		
		init_vec_mp(W_out->L_mp[ii],W_out->num_variables); W_out->L_mp[ii]->size = W_out->num_variables;
		vec_cp_mp(W_out->L_mp[ii],W_in.L_mp[ii]);
	}

	
	
	
	//set the number of points
	W_out->W.num_pts  = W_in.W.num_pts;
	W_out->W_mp.num_pts = W_in.W.num_pts;
	
	
  W_out->W.pts=(point_d *)bmalloc(W_out->W.num_pts*sizeof(point_d));
	W_out->W_mp.pts=(point_mp *)bmalloc(W_out->W.num_pts*sizeof(point_mp));
	
	for (ii=0; ii<W_in.W.num_pts; ii++) {
		init_vec_d(W_out->W.pts[ii],W_out->num_variables); W_out->W.pts[ii]->size = W_out->num_variables;
		vec_cp_d(W_out->W.pts[ii],W_in.W.pts[ii]);
		
		init_vec_mp(W_out->W_mp.pts[ii],W_out->num_variables); W_out->W_mp.pts[ii]->size = W_out->num_variables;
		vec_cp_mp(  W_out->W_mp.pts[ii],W_in.W_mp.pts[ii]);
	}
	
	return;
	
}

// initializes witness set, both the mp data and double.
//only call this on new witness sets, otherwise you will leak memory when you set the pointers to NULL.
void init_witness_set_d(witness_set_d *W){
	W->num_variables = W->num_patches = W->num_linears = 0;
	W->patch = W->L = W->W.pts = NULL;
	W->patch_mp = W->L_mp = W->W_mp.pts = NULL;
	W->variable_names = NULL;
	W->W.num_pts = W->W_mp.num_pts = 0;
	W->incidence_number = -1;
	return;
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
	printf("dehomogenizing\n");
	print_point_to_screen_matlab_mp(dehom_me,"dehom_me");
	
	if (dehom_me->size==0 || dehom_me->size==1) {
		printf("attempting to dehomogenize a vector of length 0 or 1\n");
		exit(977);
	}
	
	comp_mp denom; init_mp(denom);
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
	
	clear_mp(denom);
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
	comp_mp temp2; init_mp(temp2);
	int ii;
	for (ii=0; ii<one->size; ++ii) {
		set_mp(temp2,result)
		mul_mp(temp,&one->coord[ii],&two->coord[ii]);
		add_mp(result,temp2,temp);
	}
	clear_mp(temp);clear_mp(temp2);
}

void write_homogeneous_coordinates(witness_set_d W, char filename[])
{
	int ii,jj;
	
	FILE *OUT  = safe_fopen_write(filename); // open the output file
	
	
	fprintf(OUT,"%d\n\n",W.W_mp.num_pts); // print the header line
	
	for (ii=0; ii<W.W_mp.num_pts; ++ii) {
		for (jj=0; jj<W.num_variables; jj++) {
			fprintf(OUT,"%.15le %.15le\n",W.W.pts[ii]->coord[jj].r,W.W.pts[ii]->coord[jj].i);
		}
		fprintf(OUT,"\n");
	}
	
	fclose(OUT); // close the output file
	
	return;
}

void write_dehomogenized_coordinates(witness_set_d W, char filename[]){
	int ii,jj;
	
	FILE *OUT = safe_fopen_write(filename); // open the output file.
	
	fprintf(OUT,"%d\n\n",W.W.num_pts); // print the header line
	for (ii=0; ii<W.W.num_pts; ++ii) {
		if (W.MPType==1){ // both fields should be populated anyway?
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
	
	return;
}



void write_linears(witness_set_d W, char filename[])
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



void clear_witness_set(witness_set_d W){
	int ii;
	
	
//	printf("freeing %d patches\n",W.num_patches);
	for (ii=0; ii<W.num_patches; ii++) {
		clear_vec_d(W.patch[ii]);
		clear_vec_mp(W.patch_mp[ii]);
	}
	if (W.num_patches>0) {
		free(W.patch);
		free(W.patch_mp);
	}
	
	
//	printf("freeing %d linears\n",W.num_linears);
	for (ii=0; ii<W.num_linears; ii++) {
		clear_vec_d(W.L[ii]);
		clear_vec_mp(W.L_mp[ii]);
	}
	if (W.num_linears>0) {
		free(W.L);
		free(W.L_mp);
	}
	
	if (W.W.num_pts!=W.W_mp.num_pts) {
		printf("there was a mismatch in the number of points in the witness set being cleared\n");
	}
	
//	printf("freeing %d points\n",W.W.num_pts);
	for (ii=0; ii<W.W.num_pts; ii++) {
		clear_point_d(W.W.pts[ii]);
		clear_point_mp(W.W_mp.pts[ii]);
	}
	free(W.W.pts);
	free(W.W_mp.pts);
	

//TODO: clear variable names
	return;
}

void print_witness_set_to_screen(witness_set_d W){
	int ii;
	printf("******\n%d points in double, %d points in mp\n******\n",W.W.num_pts,W.W_mp.num_pts);
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



void endgamedata_to_endpoint(post_process_t *endPoint, endgame_data_t *EG){
	
	int num_vars, ii;
	endPoint->path_num = EG->pathNum;
	
	endPoint->sol_prec = EG->prec;
	endPoint->cond_est = EG->condition_number;
	endPoint->final_t = EG->t_val_at_latest_sample_point_d;//???
	endPoint->first_increase = EG->first_increase;
	
	
	if (EG->prec==52) {
		num_vars = EG->PD_d.point->size;
		endPoint->sol_d  = (comp_d *)bmalloc(num_vars * sizeof(comp_d));
		endPoint->sol_mp = NULL;
		
		for (ii=0; ii<num_vars; ii++) {
			endPoint->sol_d[ii]->r = EG->PD_d.point->coord[ii].r;
			endPoint->sol_d[ii]->i = EG->PD_d.point->coord[ii].i;
		}

		endPoint->size_sol = num_vars;
		endPoint->function_resid_d = EG->function_residual_d;  // the function residual
		endPoint->newton_resid_d = EG->latest_newton_residual_d;
		endPoint->cycle_num = EG->PD_d.cycle_num;
		endPoint->accuracy_estimate = EG->error_at_latest_sample_point_d;//

	}
	else
	{
		num_vars = EG->PD_mp.point->size;
		
		endPoint->sol_d  = NULL;
		endPoint->sol_mp = (comp_mp *)bmalloc(num_vars * sizeof(comp_mp));
		for (ii=0; ii<num_vars; ii++) {
			init_mp2(endPoint->sol_mp[ii],EG->prec);
			mpf_set(endPoint->sol_mp[ii]->r,EG->PD_mp.point->coord[ii].r);
			mpf_set(endPoint->sol_mp[ii]->i,EG->PD_mp.point->coord[ii].i);
		}
		endPoint->size_sol = num_vars;
		mpf_init2(endPoint->function_resid_mp, EG->prec); mpf_init2(endPoint->newton_resid_mp, EG->prec);
		
		mpf_set(endPoint->function_resid_mp,EG->function_residual_mp); //this is undoubtedly incorrect
		mpf_set(endPoint->newton_resid_mp,EG->latest_newton_residual_mp);
		endPoint->cycle_num = EG->PD_mp.cycle_num;
		endPoint->accuracy_estimate = mpf_get_d(EG->error_at_latest_sample_point_mp);
	}
	
	
	if (EG->retVal==0) {
		endPoint->success = 1;
	}
	else if (EG->retVal == retVal_sharpening_failed){
		endPoint->success = retVal_sharpening_failed;
	}
	else if (EG->retVal == retVal_sharpening_singular_endpoint){
		endPoint->success = retVal_sharpening_singular_endpoint;
	}
	else{
		endPoint->success = -1;
	}
		

	
	endPoint->sol_num = 0; // set up for post-processing
	endPoint->multiplicity = 1;
	endPoint->isFinite = 0;
//	// the post_process_t structure is used in post-processing //
//	typedef struct
//	{
//		int path_num;     // path number of the solution
//		int sol_num;      // solution number
//		comp_d  *sol_d;   // solution
//		comp_mp *sol_mp;
//		int sol_prec;     // precision of the solution
//		int size_sol;     // the number of entries in sol
//		double function_resid_d;  // the function residual
//		mpf_t  function_resid_mp;
//		double cond_est;  // the estimate of the condition number
//		double newton_resid_d;    // the newton residual
//		mpf_t  newton_resid_mp;
//		double final_t;   // the final value of time
//		double accuracy_estimate; // accuracy estimate between extrapolations
//		double first_increase;    // time value of the first increase in precision
//		int cycle_num;    // cycle number used in extrapolations
//		int success;      // success flag
//		int multiplicity; // multiplicity
//		int isReal;       // real flag:  0 - not real, 1 - real
//		int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
//		int isSing;       // singular flag: 0 - non-sigular, 1 - singular
//	} post_process_t;
	
	
	
//	typedef struct
//	{
//		int prec;
//		point_data_d PD_d;
//		point_data_mp PD_mp;
//		
//		int last_approx_prec;       // precision of the last approximation
//		point_d last_approx_d;      // last approximation to the end point
//		point_mp last_approx_mp;    // last approximation to the end point
//		
//		int retVal;
//		int pathNum;
//		int codim;
//		double first_increase;
//		double condition_number;
//		double function_residual_d;
//		mpf_t  function_residual_mp;
//		double latest_newton_residual_d;
//		mpf_t  latest_newton_residual_mp;
//		double t_val_at_latest_sample_point_d;
//		mpf_t  t_val_at_latest_sample_point_mp;
//		double error_at_latest_sample_point_d;
//		mpf_t  error_at_latest_sample_point_mp;
//	} endgame_data_t;
	
	
	
	
}

//void findSingSol(post_process_t *endPoints, point_d *dehomPoints_d, point_mp *dehomPoints_mp, int num_sols, int num_vars, preproc_data *PPD, double maxCondNum, double finalTol, int regenToggle)

int BRfindSingularSolns(post_process_t *endPoints, int num_sols, int num_vars,
												tracker_config_t *T ){
	int ii, sing_count=0;
	
	for (ii = 0; ii < num_sols; ii++){
//    if (endPoints[ii].multiplicity > 0) { 
      if ( (endPoints[ii].cond_est >  T->cond_num_threshold) || (endPoints[ii].cond_est < 0.0) )//|| (endPoints[ii].multiplicity > 1)
        endPoints[ii].isSing = 1;
      else
        endPoints[ii].isSing = 0;
			
      if (endPoints[ii].isSing)
      {
        sing_count++;
			}
//		}
	}
	
	return sing_count;
}


int BRfindFiniteSolns(post_process_t *endPoints, int num_sols, int num_vars,
												tracker_config_t *T ){
	int ii, jj, finite_count=0;
	
	//initialize temp stuffs
	comp_d dehom_coord_recip_d;
	comp_mp dehom_coord_recip_mp; init_mp(dehom_coord_recip_mp);	
	vec_d dehom_d;   init_vec_d(dehom_d,num_vars-1);   dehom_d->size = num_vars-1;
	vec_mp dehom_mp; init_vec_mp(dehom_mp,num_vars-1); dehom_mp->size = num_vars-1;
	
	
	
	for (ii = 0; ii < num_sols; ii++){
		if (endPoints[ii].sol_prec<64) {
			set_d(dehom_coord_recip_d,endPoints[ii].sol_d[0]);
			recip_d(dehom_coord_recip_d,dehom_coord_recip_d);
			for (jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_d(&dehom_d->coord[jj],dehom_coord_recip_d,endPoints[ii].sol_d[jj])
			}
			
			if (infNormVec_d(dehom_d) < T->finiteThreshold){
				endPoints[ii].isFinite = 1;
				finite_count++;
			}
			else{
				endPoints[ii].isFinite = 0;
			}
			
		}
		else // high precision, do mp
		{
			change_prec_point_mp(dehom_mp,endPoints[ii].sol_prec);
			setprec_mp(dehom_coord_recip_mp,endPoints[ii].sol_prec);
			set_mp(dehom_coord_recip_mp,endPoints[ii].sol_mp[0]);
			for (jj=0; jj<num_vars-1; ++jj) {
				//do the division.
				mul_mp(&dehom_mp->coord[jj],dehom_coord_recip_mp,endPoints[ii].sol_mp[jj])
			}
			
			if (infNormVec_mp(dehom_mp) < T->finiteThreshold){
				endPoints[ii].isFinite = 1;
				finite_count++;
			}
			else{
				endPoints[ii].isFinite = 0;
			}
			
		}
	}
	
	clear_vec_d(dehom_d);
	clear_vec_mp(dehom_mp);
	clear_mp(dehom_coord_recip_mp);
	
	return finite_count;
}






//assumes that W has the number of variables already set, and the pts NOT allocated yet.  should be NULL
void BRpostProcessing_AllowDuplicates(post_process_t *endPoints, witness_set_d *W_new, int num_pts, preproc_data preProcData, tracker_config_t *T)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does the actual post processing for a zero dim run     *
 \***************************************************************/
{
	
	int ii, jj;
	//	printf("%lf\n",T_copy[oid].finiteThreshold);
	//direct from the bertini library:
	findMultSol(endPoints, num_pts, W_new->num_variables, preProcData, T->final_tol_times_mult);
	// sets the multiplicity and solution number in the endPoints data
	
	//custom, derived from bertini's analagous call.
	int num_singular_solns = BRfindSingularSolns(endPoints, num_pts, W_new->num_variables, T);
	//sets the singularity flag in endPoints.
	
	int num_finite_solns = BRfindFiniteSolns(endPoints, num_pts, W_new->num_variables, T);
	
	printf("%d singular solutions\n",num_singular_solns);
	printf("%d finite solutions\n",num_finite_solns);
	
	for (ii=0; ii<num_pts; ++ii) {
		//		int success;      // success flag
		//		int multiplicity; // multiplicity
		//		int isReal;       // real flag:  0 - not real, 1 - real
		//		int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
		//		int isSing;       // singular flag: 0 - non-sigular, 1 - singular
		
		printf("solution %d, success %d, multi %d, isFinite %d, isSing %d\n",ii,endPoints[ii].success,endPoints[ii].multiplicity,endPoints[ii].isFinite,endPoints[ii].isSing);
	}
	
	int *actual_solns_indices;
	actual_solns_indices = (int *)bmalloc(num_pts*sizeof(int));
	
	
	int num_actual_solns = 0;
	for (ii=0; ii<num_pts; ii++) {
		if (endPoints[ii].isFinite &&  (endPoints[ii].success==1) ) { //(!endPoints[ii].isSing) && (endPoints[ii].multiplicity==1) &&
			actual_solns_indices[num_actual_solns] = ii;
			num_actual_solns++;
		}
	}
	

	//initialize the structures for holding the produced data
	
	W_new->W.num_pts=num_actual_solns; W_new->W_mp.num_pts=num_actual_solns;
  W_new->W.pts=(point_d *)bmalloc(num_actual_solns*sizeof(point_d));
  W_new->W_mp.pts=(point_mp *)bmalloc(num_actual_solns*sizeof(point_mp));
	
	
	
	for (ii=0; ii<num_actual_solns; ++ii) {
		printf("setting the actual solution %d\n",ii);
		init_vec_d(W_new->W.pts[ii],W_new->num_variables); init_vec_mp(W_new->W_mp.pts[ii],W_new->num_variables);
		W_new->W.pts[ii]->size = W_new->W_mp.pts[ii]->size = W_new->num_variables;
		
		if (endPoints[actual_solns_indices[ii]].sol_prec<64) {
			//copy out of the double structure.
			for (jj=0; jj<W_new->num_variables; jj++) {
				set_d(&W_new->W.pts[ii]->coord[jj],endPoints[actual_solns_indices[ii]].sol_d[jj]);
			}
			vec_d_to_mp(W_new->W_mp.pts[ii],W_new->W.pts[ii]);
		}
		else{
			//copy out of the mp structure.
			for (jj=0; jj<W_new->num_variables; jj++) {
				set_mp(&W_new->W_mp.pts[ii]->coord[jj],endPoints[actual_solns_indices[ii]].sol_mp[jj]);
			}
			vec_mp_to_d(W_new->W.pts[ii],W_new->W_mp.pts[ii]);
		}
	}
	
	free(actual_solns_indices);
	
  return;
}





//assumes that W has the number of variables already set, and the pts NOT allocated yet.  should be NULL
void BRpostProcessing(post_process_t *endPoints, witness_set_d *W_new, int num_pts, preproc_data preProcData, tracker_config_t *T)
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES: does the actual post processing for a zero dim run     *
 \***************************************************************/
{
	
	int ii, jj;
	//	printf("%lf\n",T_copy[oid].finiteThreshold);
	//direct from the bertini library:
	findMultSol(endPoints, num_pts, W_new->num_variables, preProcData, T->final_tol_times_mult);
	// sets the multiplicity and solution number in the endPoints data
	
	//custom, derived from bertini's analagous call.
	int num_singular_solns = BRfindSingularSolns(endPoints, num_pts, W_new->num_variables, T);
	//sets the singularity flag in endPoints.
	
	int num_finite_solns = BRfindFiniteSolns(endPoints, num_pts, W_new->num_variables, T);
	
	printf("%d singular solutions\n",num_singular_solns);
	printf("%d finite solutions\n",num_finite_solns);
	
	for (ii=0; ii<num_pts; ++ii) {
		//		int success;      // success flag
		//		int multiplicity; // multiplicity
		//		int isReal;       // real flag:  0 - not real, 1 - real
		//		int isFinite;     // finite flag: -1 - no finite/infinite distinction, 0 - infinite, 1 - finite
		//		int isSing;       // singular flag: 0 - non-sigular, 1 - singular
		
		printf("solution %d, success %d, multi %d, isFinite %d, isSing %d\n",ii,endPoints[ii].success,endPoints[ii].multiplicity,endPoints[ii].isFinite,endPoints[ii].isSing);
	}
	
	int num_actual_solns = 0;
	int *actual_solns_indices;
	actual_solns_indices = (int *)bmalloc(num_pts*sizeof(int));
	
	
	for (ii=0; ii<num_pts; ii++) {
		if (endPoints[ii].isFinite && (!endPoints[ii].isSing) && (endPoints[ii].multiplicity==1) && (endPoints[ii].success==1) ) {
			actual_solns_indices[num_actual_solns] = ii;
			num_actual_solns++;
		}
	}
	
//TODO: here, get the good points, so can copy them
	
	//initialize the structures for holding the produced data

	W_new->W.num_pts=num_actual_solns; W_new->W_mp.num_pts=num_actual_solns;
  W_new->W.pts=(point_d *)bmalloc(num_actual_solns*sizeof(point_d));
  W_new->W_mp.pts=(point_mp *)bmalloc(num_actual_solns*sizeof(point_mp));
	
	
	
	for (ii=0; ii<num_actual_solns; ++ii) {
		printf("setting the actual solution %d\n",ii);
		init_vec_d(W_new->W.pts[ii],W_new->num_variables); init_vec_mp(W_new->W_mp.pts[ii],W_new->num_variables);
		W_new->W.pts[ii]->size = W_new->W_mp.pts[ii]->size = W_new->num_variables;
		
		if (endPoints[actual_solns_indices[ii]].sol_prec<64) {
			//copy out of the double structure.
			for (jj=0; jj<W_new->num_variables; jj++) {
				set_d(&W_new->W.pts[ii]->coord[jj],endPoints[actual_solns_indices[ii]].sol_d[jj]);
			}
			vec_d_to_mp(W_new->W_mp.pts[ii],W_new->W.pts[ii]);
		}
		else{
			//copy out of the mp structure.
			for (jj=0; jj<W_new->num_variables; jj++) {
				set_mp(&W_new->W_mp.pts[ii]->coord[jj],endPoints[actual_solns_indices[ii]].sol_mp[jj]);
			}
			vec_mp_to_d(W_new->W.pts[ii],W_new->W_mp.pts[ii]);
		}
	}
	
	free(actual_solns_indices);

	
	
	
  return;
}







void insert_randomization_matrix_witness_data(int rows, int cols, int codim_index){
	
	
	size_t len = 0;
	int ii,jj,kk;
	
	
	//make the matrix
	mat_mp randomizer_matrix; init_mat_mp2(randomizer_matrix,rows,cols,1024);
	randomizer_matrix->rows = rows; randomizer_matrix->cols = cols;
	make_matrix_random_mp(randomizer_matrix,rows,cols,1024);
	
	
	
	//rename the witness_data file for re-use later.
	rename("witness_data","witness_data_predeflation");
	
	FILE *IN = NULL, *OUT = NULL;
	
	IN = safe_fopen_read("witness_data_predeflation");
	
	OUT = safe_fopen_write("witness_data_newrand");
	
	
	int *nonempty_codimensions;
	
	char ch;
	int num_variables, num_nonempty_codims, num_pts_this_codim, current_codimension;
	int bytes_read;
	fscanf(IN,"%d\n",&num_variables); fprintf(OUT,"%d\n",num_variables);
	fscanf(IN,"%d\n",&num_nonempty_codims);fprintf(OUT,"%d\n",num_nonempty_codims);
	char *line = NULL;  // getline will alloc
	
	
	for (ii=0; ii<num_nonempty_codims; ii++) {
		fscanf(IN,"%d\n",&current_codimension); fprintf(OUT,"%d\n",current_codimension);
		fscanf(IN,"%d\n",&num_pts_this_codim); fprintf(OUT,"%d\n",num_pts_this_codim);
		
		
		for (jj=0; jj<num_pts_this_codim; jj++) {
			bytes_read = getline(&line, &len, IN); fprintf(OUT,"%s",line);  free(line); line=NULL; // precision of soln
			
			//copy the solution itself.
			for (kk=0; kk<num_variables; kk++) {
				bytes_read = getline(&line, &len, IN); fprintf(OUT,"%s",line);  free(line); line=NULL; // precision of soln
			}
			
			bytes_read = getline(&line, &len, IN); fprintf(OUT,"%s",line);  free(line); line=NULL; // precision of soln
			
			//copy the previous approximation
			for (kk=0; kk<num_variables; kk++) {
				bytes_read = getline(&line, &len, IN); fprintf(OUT,"%s",line);  free(line); line=NULL; // precision of soln
			}
			
			for (kk=0; kk<8; kk++) { //  copy the next 8 lines.
				bytes_read = getline(&line, &len, IN); fprintf(OUT,"%s",line);  free(line); line=NULL; // precision of soln
			}
			
			// then do it again!
		}
	}
	
	int checker;
	fscanf(IN,"%d\n",&checker); fprintf(OUT,"-1\n\n");
	if (checker!=-1) {
		printf("the check integer was not -1\n");
		exit(-1);
	}
	
	int num_type;
	fscanf(IN,"%d\n",&num_type);
	
//	for (ii=0; ii<num_nonempty_codims; ++ii) {
//		//for each nonempty codimension, copy the rest of the info, unless codimen matches our target.
//		if (codim_index == ) {
//			
//		}
//	}

	
	fclose(IN); fclose(OUT);
	clear_mat_mp(randomizer_matrix);
	
	//rename the witness_data file for re-use later.
	rename("witness_data_predeflation","witness_data");
	
	
	return;
}





