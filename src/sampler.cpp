#include "sampler.hpp"


int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	
	MPI_Init(&argC,&args);
	
	
	int num_vars=0;
	boost::filesystem::path inputName, RandMatName, witnessSetName, samplingNamenew;
	
	vertex_set V;
	curve_decomposition C;   //new data type; stores vertices, edges, etc.
	
	sample_data   S_old,S_new;
	mat_mp randomizer_matrix;
	
	
	srand(time(NULL));
	
	
	sampler_configuration sampler_options;
	
	sampler_options.splash_screen();
	sampler_options.parse_commandline(argC, args);

	int MPType;
	
	boost::filesystem::path Dir_Name;
	get_dir_mptype( Dir_Name, &MPType);
	
	
	
	
	int successful_startup = curve_sampler_startup(Dir_Name,
																								 inputName, witnessSetName,RandMatName,samplingNamenew,
																								 C,V);
	if (successful_startup!=1)
		return -445;
	
	
	
	parse_input_file(inputName);
	

	
	// set up the solver configuration
	solver_configuration solve_options;  solver_init_config(&solve_options);
	get_tracker_config(&solve_options,MPType);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	initMP(solve_options.T.Precision);
	
	num_vars = get_num_vars_PPD(solve_options.PPD);		
	
	


	
	if (solve_options.verbose_level>=1)
		printf("loading curve decomposition\n");
	
	witness_set W;
	W.num_variables = num_vars;
	W.MPType = MPType;
	W.get_variable_names();
	
//	W.witnessSetParse(witnessSetName,num_vars);
	
	
	W.reset_patches();
	for (int ii=0; ii<C.num_patches; ii++) {
		W.add_patch(C.patch[ii]);
	}
//	W.read_patches_from_file(Dir_Name / "patches");
	read_matrix(RandMatName, randomizer_matrix);
	
	set_initial_sample_data(&S_old,C,V,
											 num_vars);

	solve_options.verbose_level = sampler_options.verbose_level;
	if (solve_options.verbose_level>=2) {
		solve_options.show_status_summary=1;
	}
	solve_options.T.ratioTol = 1; // manually assert to be more permissive.  i don't really like this.
	solve_options.use_midpoint_checker = 0;
	solve_options.use_gamma_trick = sampler_options.use_gamma_trick;
	solve_options.allow_singular = 1;
	
	/////////
	////////
	//////
	////
	//
	//  Generate new sampling data
	//

	
	generate_new_sampling_pts(&S_new,
														randomizer_matrix,
														S_old,
														C,
														V,
														W,
														MPType,
														&sampler_options,
														&solve_options);
	
	

	//
	//   Done with the main call
	////
	/////
	///////
	////////
	
	
	
	//output
	output_sampling_data(S_new,V,
											 samplingNamenew,num_vars,MPType);
	


	
	clear_mat_mp(randomizer_matrix);
		
	clear_sample(&S_new, MPType);
	
	
	clearMP();
	
	MPI_Finalize();
	return 0;
}





//
//
//void generate_new_sampling_pts(sample_data *S_new,
//															 mat_mp randomizer_matrix,
//															 sample_data S_old,
//															 curve_decomposition C,
//															 vertex_set *V,
//															 witness_set W,
//															 int MPType ,
//															 sampler_configuration *sampler_options,
//															 solver_configuration *solve_options)
//{
//	
//	if(MPType==0)
//		generate_new_sampling_pts_d(S_new, randomizer_matrix, S_old, C, V, W, MPType, sampler_options,solve_options);
//	else
//		generate_new_sampling_pts_d(S_new, randomizer_matrix, S_old, C, V, W, MPType, sampler_options, solve_options);
//	
//	
//	return;
//}


void generate_new_sampling_pts(sample_data *S_new,
																 mat_mp randomizer_matrix,
																 sample_data S_old,
																 curve_decomposition C,
																 vertex_set &V,
																 witness_set & W,
																 int MPType,
																 sampler_configuration *sampler_options,
																 solver_configuration *solve_options)
{

	
	int				num_vars = W.num_variables;
	
	
	witness_set  Wnew;
	
	
	vec_mp				target_projection;
	init_vec_mp(target_projection,num_vars); target_projection->size = num_vars;
	vec_cp_mp(target_projection,C.pi_mp[0]); // copy the projection into target_projection
	
	
	vec_mp startpt;
	init_vec_mp(startpt,num_vars); startpt->size = num_vars;
	
	
	vec_mp				start_projection;
	init_vec_mp(start_projection,num_vars);  start_projection->size = num_vars;
	vec_cp_mp(start_projection,C.pi_mp[0]); // grab the projection, copy it into start_projection
	
	
	int ii,jj;
	
	
	comp_mp				temp, temp1, target_projection_value;
	init_mp(temp);  init_mp(temp1); init_mp(target_projection_value);

	int				prev_num_samp, sample_counter;
	int				*refine_current = NULL, *refine_next = NULL;
	
	
	
	
	S_new->num_edges = S_old.num_edges;
	S_new->num_samples_each_edge = (int *) bmalloc(S_new->num_edges * sizeof(int));
	S_new->sample_indices = (int **) br_malloc(S_new->num_edges * sizeof(int *));
	
	
	vertex temp_vertex;
//	init_vertex(&temp_vertex);
	mpf_t dist_away; mpf_init(dist_away);
	int interval_counter;
	int num_refinements;
	int *current_indices;
	
	
	for(ii=0;ii<S_old.num_edges;ii++) // for each of the edges
	{

		set_initial_refinement_flags(&num_refinements,
																 &refine_current,
																 &current_indices,
																 S_old, V,
																 ii, sampler_options);
		
		prev_num_samp = S_old.num_samples_each_edge[ii]; // grab the number of points from the array of integers
		
		
		
		int pass_number  = 0;//this should be the only place this is reset.
		while(1) // breaking condition is all samples being less than TOL away from each other (in the infty norm sense).
		{

			refine_next = (int * )br_malloc((prev_num_samp+num_refinements-1) * sizeof(int)); // refinement flag	

			int *new_indices = (int *) br_malloc((prev_num_samp+num_refinements)*sizeof(int));
			
			new_indices[0] = current_indices[0];
			sample_counter = 1; // this will be incremented every time we put a point into new_points
			// starts at 1 because already committed one.
			// this should be the only place this counter is reset.
			
			
			if (sampler_options->verbose_level>=0)
				printf("edge %d, pass %d, %d refinements\n",ii,pass_number,num_refinements);
			
			if (sampler_options->verbose_level>=1) {
				printf("the current indices:\n");
				for (jj=0; jj<prev_num_samp; jj++) 
					printf("%d ",current_indices[jj]);
				printf("\n\n");
			}
			
			if (sampler_options->verbose_level>=1) {
				printf("refine_flag:\n");
				for (jj=0; jj<prev_num_samp-1; jj++) {
						printf("%d ",refine_current[jj]);
				}
				printf("\n\n");
			}
			
			if (sampler_options->verbose_level>=1) {
				printf("the current projection values are:\n");
				for (jj=0; jj<prev_num_samp; jj++) {
					print_comp_matlab(V.vertices[current_indices[jj]].projVal_mp,"proj");
				}
				printf("\n\n");
			}
			
			num_refinements = 0; // reset this counter.  this should be the only place this is reset
			interval_counter = 0;
			for(jj=0;jj<prev_num_samp-1;jj++) // for each interval in the previous set
			{
				
				
				
				if (sampler_options->verbose_level>=2) 
					printf("interval %d of %d\n",jj,prev_num_samp-1);
				
				
				
				int startpt_index; int left_index; int right_index;
				
				// set the starting projection and point.
				if(jj==0){// if on the first sample, go the right
					startpt_index = current_indices[1]; //right!
				}
				else{ //go to the left
					startpt_index = current_indices[jj]; // left!
				}
				
				left_index = current_indices[jj];
				right_index = current_indices[jj+1];
				
				
				if (new_indices[sample_counter-1]!=left_index) {
					printf("index mismatch\n");
					exit(1);
				}
				

				
				if(refine_current[jj]==1) //
				{
			
					vec_cp_mp(startpt,V.vertices[startpt_index].pt_mp);
					set_mp(&(start_projection->coord[0]),V.vertices[startpt_index].projVal_mp);
					neg_mp(&(start_projection->coord[0]),&(start_projection->coord[0]));
					
					
//					print_point_to_screen_matlab(C.pi_mp[0],"C.pi_mp[0]");
					
					
					estimate_new_projection_value(target_projection_value,				// the new value
																				V.vertices[left_index].pt_mp,	//
																				V.vertices[right_index].pt_mp, // two points input
																				C.pi_mp[0]);												// projection (in homogeneous coordinates)
					
					
					
//					print_comp_matlab(V.vertices[left_index].projVal_mp,"left");
//					print_comp_matlab(V.vertices[right_index].projVal_mp,"right");
//					print_comp_matlab(target_projection_value,"target");
//					
//					printf("\n\n");
					
					
					neg_mp(&target_projection->coord[0],target_projection_value); // take the opposite :)
					
					
					set_witness_set_mp(&W, start_projection,startpt,num_vars); // set the witness point and linear in the input for the lintolin solver.
					
					if (sampler_options->verbose_level>=3) {
						print_point_to_screen_matlab(W.pts_mp[0],"startpt");
						print_comp_matlab(&W.L_mp[0]->coord[0],"initial_projection_value");
						print_comp_matlab(target_projection_value,"target_projection_value");
						
//						W.print_to_screen();
					}
					
					
					
					multilintolin_solver_main(MPType,
																		W,         // witness_set
																		randomizer_matrix,
																		&target_projection, //  the set of linears we will solve at.
																		&Wnew, // the new data is put here!
																		solve_options); // already a pointer
					
					
					
					
					
					if (sampler_options->verbose_level>=3) 
						print_point_to_screen_matlab(Wnew.pts_mp[0], "new_solution");
					
					
					
					
					
					// check how far away we were from the LEFT interval point
					norm_of_difference(dist_away,
														 Wnew.pts_mp[0], // the current new point
														 V.vertices[left_index].pt_mp);// jj is left, jj+1 is right
					
					if ( mpf_cmp(dist_away, sampler_options->TOL )>0 ){
						refine_next[interval_counter] = 1; 
						num_refinements++;
					}
					else{
					 refine_next[interval_counter] = 0;
					}
					interval_counter++;
					
					
					
					
					// check how far away we were from the RIGHT interval point
					norm_of_difference(dist_away,
														 Wnew.pts_mp[0], // the current new point
														 V.vertices[right_index].pt_mp);
					
					if (mpf_cmp(dist_away, sampler_options->TOL ) > 0){
						refine_next[interval_counter] = 1; 
						num_refinements++;
					}
					else{
						refine_next[interval_counter] = 0;
					}
					interval_counter++;
					
					
					vec_cp_mp(temp_vertex.pt_mp,Wnew.pts_mp[0]);
					set_mp(temp_vertex.projVal_mp,target_projection_value);
					temp_vertex.type = SAMPLE_POINT;
					new_indices[sample_counter] = V.add_vertex(temp_vertex);
					sample_counter++;
					
					new_indices[sample_counter] = right_index;
					sample_counter++;

					Wnew.reset();
					
				}
				else {
					if (sampler_options->verbose_level>=2) 
						printf("adding sample %d\n",sample_counter);
					
					refine_next[interval_counter] = 0;
					new_indices[sample_counter] = right_index;
					interval_counter++;
					sample_counter++;
				}
			}

			
			if (sampler_options->verbose_level>=1) // print by default
				printf("\n\n");
			
			if( (num_refinements == 0) || (pass_number >= sampler_options->maximum_num_iterations) ) // if have no need for new samples
			{
				
				if (sampler_options->verbose_level>=1) // print by default
					printf("breaking\nsample_counter = %d\n",sample_counter);
				
				
				S_new->sample_indices[ii] = new_indices;
				S_new->num_samples_each_edge[ii] = sample_counter;
		
				free(refine_current);
				refine_current = NULL;
				break; // BREAKS THE WHILE LOOP
			}
			else{
				free(refine_current);
				free(current_indices);
				
				refine_current = refine_next; // reassign this pointer
				current_indices = new_indices;

				prev_num_samp=sample_counter; // update the number of samples
				pass_number++;
			}

		}//while loop
		if (sampler_options->verbose_level>=2) {
			printf("exiting while loop\n");
		}
		printf("\n");
	}  // re: ii (for each edge)
	
	
	
	
	clear_mp(temp); clear_mp(temp1);
	clear_vec_mp(start_projection);
	clear_vec_mp(target_projection);
}








void set_witness_set_mp(witness_set *W, vec_mp new_linear,vec_mp pts,int num_vars)
{
	
	W->num_variables = num_vars;
	
	W->num_pts=1;
	W->pts_mp=(point_mp *)bmalloc(sizeof(point_mp)); // apparently you can only pass in a single point to copy in.
	
	//initialize the memory
	init_point_mp(W->pts_mp[0],num_vars); W->pts_mp[0]->size = num_vars;
	point_cp_mp(W->pts_mp[0],pts);
	
	
	W->num_linears = 1;
	W->L_mp = (vec_mp *)bmalloc(sizeof(vec_mp));
	init_vec_mp(W->L_mp[0],num_vars); W->L_mp[0]->size = num_vars;
	vec_cp_mp(W->L_mp[0],new_linear);
	
}


void  output_sampling_data(sample_data S, vertex_set V,
													 boost::filesystem::path samplingName,int num_vars,int MPType)
{
	FILE *OUT =  safe_fopen_write(samplingName);
	int ii,jj,kk;
	// output the number of vertices
	fprintf(OUT,"%d\n\n",S.num_edges);
	for (ii=0; ii<S.num_edges; ii++) {
		
		
		fprintf(OUT,"%d\n\n",S.num_samples_each_edge[ii]);
		for (jj=0; jj<S.num_samples_each_edge[ii]; jj++) {
			print_mp(OUT,0,V.vertices[S.sample_indices[ii][jj]].projVal_mp);
			fprintf(OUT,"\n");
			for(kk=0;kk<num_vars;kk++) {
				print_mp(OUT, 0, &V.vertices[S.sample_indices[ii][jj]].pt_mp->coord[kk]);
				fprintf(OUT,"\n");
			}
			fprintf(OUT,"\n");
		}
	}
}


int  curve_sampler_startup(boost::filesystem::path directoryName,
													 boost::filesystem::path &inputName,
													 boost::filesystem::path &witnessSetName,
													 boost::filesystem::path &RandMatName,
													 boost::filesystem::path &samplingNamenew,
													 curve_decomposition &C,
													 vertex_set &V)
/***************************************************************\
 * USAGE:    curve_sampler_startup structure and inputname
 * ARGUMENTS:                                                    *
 * RETURN VALUES:           an integer                           *
 * NOTES:                                                        *
 \***************************************************************/
{


	int successful_startup = 1;
	//read in the gross information from the summary file.
	C.setup(directoryName / "C.curve", inputName, directoryName);
	
	//setup E structure from E.edge
	C.setup_edges(directoryName / "E.edge");
	
	//setup V structure from V.vertex
	
	V.setup_vertices(directoryName / "V.vertex");
	

	witnessSetName = directoryName / "witness_set";
	

	RandMatName = directoryName / "Rand_Matrix";
	
	samplingNamenew = directoryName / "samp.dat";

	
	if(C.num_edges==0)
	{
		printf("sampler will not generate sampling data since there are no edges\n");
		successful_startup = 0;
	}
	
	return successful_startup;
}


int get_dir_mptype(boost::filesystem::path & Dir_Name, int * MPType){

	std::string tempstr;
	std::ifstream fin("Dir_Name");
	fin >> tempstr;
	fin >> *MPType;
	
	Dir_Name = tempstr;
	return *MPType;
}












//dehomogenizes, takes the average, computes the projection.
//takes in the full projection \pi, including the homogenizing coordinate.
void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi){
	int ii;
	
	if (left->size != right->size) {
		printf("attempting to estimate_new_projection_value on vectors of different size\n");
		exit(8776);
	}
	
	
	
	
	vec_mp dehom_left, dehom_right;
	init_vec_mp(dehom_left,left->size-1);   dehom_left->size  = left->size-1;
	init_vec_mp(dehom_right,right->size-1); dehom_right->size = right->size-1;
	
	dehomogenize(&dehom_left,left);
	dehomogenize(&dehom_right,right);
	
	
	
//	print_point_to_screen_matlab(dehom_left,"left_pt");
//	print_point_to_screen_matlab(dehom_right,"right_pt");
//	print_point_to_screen_matlab(pi,"pi");
	
	comp_mp temp, temp2, half; init_mp(temp); init_mp(temp2);  init_mp(half);
	
	mpf_set_d(half->r, 0.5); mpf_set_d(half->i, 0.0);
	
	 
	set_zero_mp(result);                                           // result = 0;  initialize
	
	for (ii = 0; ii<dehom_left->size; ii++) {
		add_mp(temp,&dehom_left->coord[ii],&dehom_right->coord[ii]); //  a = (x+y)
		mul_mp(temp2, temp, half);                                   //  b = a/2
		mul_mp(temp,&pi->coord[ii+1],temp2);                          //  a = b.pi
		set_mp(temp2,result);                                        //  b = result
		add_mp(result, temp, temp2);                                  //  result = a+b
	}

	
	// in other words, result += (x+y)/2 \cdot pi

	mpf_t zerothresh; mpf_init(zerothresh);
	mpf_set_d(zerothresh, 1e-9);
	if (mpf_cmp(result->i, zerothresh) < 0){
		mpf_set_d(result->i, 0.0);
	}
	
	clear_mp(temp); clear_mp(temp2); clear_mp(half);
	mpf_clear(zerothresh);
	clear_vec_mp(dehom_right);clear_vec_mp(dehom_left);
	
	return;
}




int  set_initial_sample_data(sample_data *S, curve_decomposition C, vertex_set V,
													int num_vars)
{
	int ii;
	
	
	comp_mp temp_mp;  init_mp(temp_mp);
	
	
	S->num_variables = num_vars;
	S->num_edges = C.num_edges;
	S->num_samples_each_edge = (int *)bmalloc(C.num_edges * sizeof(int));
	for(ii=0;ii<C.num_edges;ii++){
		S->num_samples_each_edge[ii]=3; // left, right, mid
	}
	
	
	S->sample_indices = (int **)br_malloc(C.num_edges*sizeof(int *));
	
	
	for(ii=0;ii<S->num_edges;ii++)
	{
		
		S->sample_indices[ii] = (int *) br_malloc(3*sizeof(int));
		S->sample_indices[ii][0] = C.edges[ii].left;
		S->sample_indices[ii][1] = C.edges[ii].midpt;
		S->sample_indices[ii][2] = C.edges[ii].right;
		
	}
	return 0;
}




void set_initial_refinement_flags(int *num_refinements, int **refine_flags, int ** current_indices,
																 sample_data S, vertex_set &V,
																 int current_edge, sampler_configuration *sampler_options)
{
	
	
	
//	printf("setting refinement flags\n");
	*num_refinements = 0;
	
	if (*refine_flags==NULL)
		(* refine_flags) = (int *)br_malloc((S.num_samples_each_edge[current_edge]-1)*sizeof(int));
	else
		(* refine_flags) = (int *)brealloc((refine_flags), (S.num_samples_each_edge[current_edge]-1)*sizeof(int));
	
	(* current_indices) = (int *)br_malloc(S.num_samples_each_edge[current_edge]*sizeof(int));
	
	
	memset (* refine_flags, 0, (S.num_samples_each_edge[current_edge]-1) *sizeof(int));
	
	int ii;
	mpf_t dist_away;  mpf_init(dist_away);
	
	(* current_indices)[0] = S.sample_indices[current_edge][0];
	for (ii=0; ii<(S.num_samples_each_edge[current_edge]-1); ii++) {
		(* current_indices)[ii+1] = S.sample_indices[current_edge][ii+1];
		
		norm_of_difference(dist_away, V.vertices[S.sample_indices[current_edge][ii]].pt_mp,
																	V.vertices[S.sample_indices[current_edge][ii+1]].pt_mp); // get the distance between the two adjacent points.
		if ( mpf_cmp(dist_away, sampler_options->TOL)>0 ){
			(*refine_flags)[ii] = 1;
			(*num_refinements)++;
		}
	}
//	printf("done setting refinement flags\n");
	mpf_clear(dist_away);
}




