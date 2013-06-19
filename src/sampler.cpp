#include "sampler.hpp"


int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	int num_vars=0;
	char *inputName = NULL, *RandMatName = NULL, *witnessSetName = NULL, *samplingNamenew = NULL;
	
	vertex_set V;
	curveDecomp_d C;   //new data type; stores vertices, edges, etc.
	witness_set W; init_witness_set(&W);
	sample_data   S_old,S_new;
	mat_mp n_minusone_randomizer_matrix;
	
	
	srand(time(NULL));
	
	
	
	sampler_splash_screen();
	sampler_configuration sampler_options;  sampler_init_config(&sampler_options);
	sampler_parse_commandline(argC, args, &sampler_options);

	int MPType;
	char *Dir_Name = NULL;
	get_dir_mptype( &Dir_Name, &MPType);
	
	
	init_curveDecomp_d(&C);
	init_vertex_set(&V);
	
	
	int successful_startup = curve_sampler_startup(argC, args,
																								 Dir_Name,
																								 &inputName, &witnessSetName,&RandMatName,&samplingNamenew,
																								 &C,&V,
																								 &num_vars, MPType);
	if (successful_startup!=1)
		return -445;
	
	
	printf("%s\n",inputName);
	parse_input_file(inputName, &MPType);
	

	
	// set up the solver configuration
	solver_configuration solve_options;  solver_init_config(&solve_options);
	get_tracker_config(&solve_options,MPType);
	parse_preproc_data("preproc_data", &solve_options.PPD);
	
	initMP(solve_options.T.Precision);
	
	num_vars = get_num_vars_PPD(solve_options.PPD);		
	
	


	
	if (solve_options.verbose_level>=1)
		printf("loading curve decomposition\n");
	
	witnessSetParse(&W, witnessSetName,num_vars);
	W.num_var_gps = solve_options.PPD.num_var_gp;
	W.MPType = MPType;
	get_variable_names(&W);

	read_rand_matrix(RandMatName, n_minusone_randomizer_matrix);
	
	set_initial_sample_data(&S_old,C,V,
											 num_vars);
	

	solve_options.verbose_level = sampler_options.verbose_level-1;
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

	if (solve_options.verbose_level>=1)
		printf("generating new_sample points\n");

	
	generate_new_sampling_pts(&S_new,
														n_minusone_randomizer_matrix,
														S_old,
														C,
														&V,
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
	
	
	
	
	
	// clear memory
	free(inputName);
	free(witnessSetName);
	
	clear_mat_mp(n_minusone_randomizer_matrix);
	clear_witness_set(W);
	
	
	clear_sample_d(&S_new, MPType);
	
//	clear_sample_d(&S_old, MPType);
	
	clear_curveDecomp_d(&C);
	clearMP();
	return 0;
}





//
//
//void generate_new_sampling_pts(sample_data *S_new,
//															 mat_mp n_minusone_randomizer_matrix,
//															 sample_data S_old,
//															 curveDecomp_d C,
//															 vertex_set *V,
//															 witness_set W,
//															 int MPType ,
//															 sampler_configuration *sampler_options,
//															 solver_configuration *solve_options)
//{
//	
//	if(MPType==0)
//		generate_new_sampling_pts_d(S_new, n_minusone_randomizer_matrix, S_old, C, V, W, MPType, sampler_options,solve_options);
//	else
//		generate_new_sampling_pts_d(S_new, n_minusone_randomizer_matrix, S_old, C, V, W, MPType, sampler_options, solve_options);
//	
//	
//	return;
//}


void generate_new_sampling_pts(sample_data *S_new,
																 mat_mp n_minusone_randomizer_matrix,
																 sample_data S_old,
																 curveDecomp_d C,
																 vertex_set *V,
																 witness_set W,
																 int MPType,
																 sampler_configuration *sampler_options,
																 solver_configuration *solve_options)
{

	
	int				num_vars = W.num_variables;
	
	
	witness_set  Wnew;
	
	
	vec_mp				target_projection;
	init_vec_mp(target_projection,num_vars); target_projection->size = num_vars;
	vec_cp_mp(target_projection,C.pi_mp); // copy the projection into target_projection
	
	
	vec_mp startpt;
	init_vec_mp(startpt,num_vars); startpt->size = num_vars;
	
	
	vec_mp				start_projection;
	init_vec_mp(start_projection,num_vars);  start_projection->size = num_vars;
	vec_cp_mp(start_projection,C.pi_mp); // grab the projection, copy it into start_projection
	
	
	int ii,jj;
	
	
	comp_mp				temp, temp1, target_projection_value;
	init_mp(temp);  init_mp(temp1); init_mp(target_projection_value);

	int				prev_num_samp, sample_counter;
	int				*refine_current = NULL, *refine_next = NULL;
	
	
	
	
	S_new->num_edges = S_old.num_edges;
	S_new->num_samples_each_edge = (int *) bmalloc(S_new->num_edges * sizeof(int));
	S_new->sample_indices = (int **) br_malloc(S_new->num_edges * sizeof(int *));
	
	
	vertex temp_vertex;  init_vertex(&temp_vertex);
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
//			memset (refine_next, 0, (prev_num_samp+num_refinements-1) *sizeof(int)); // uncomment me and remove set=0 statements later
			

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
					print_comp_mp_matlab(V->vertices[current_indices[jj]].projVal_mp,"proj");
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
					
					

			
					vec_cp_mp(startpt,V->vertices[startpt_index].pt_mp);
					set_mp(&(start_projection->coord[0]),V->vertices[startpt_index].projVal_mp);
					neg_mp(&(start_projection->coord[0]),&(start_projection->coord[0]));
					
					
					
					estimate_new_projection_value(target_projection_value,				// the new value
																				V->vertices[left_index].pt_mp,	//
																				V->vertices[right_index].pt_mp, // two points input
																				C.pi_mp);												// projection (in homogeneous coordinates)
					
					
//					print_comp_mp_matlab(V->vertices[left_index].projVal_mp,"left");
//					print_comp_mp_matlab(V->vertices[right_index].projVal_mp,"right");
//					print_comp_mp_matlab(target_projection_value,"target");
//					
//					printf("\n\n");
					
					
					neg_mp(&target_projection->coord[0],target_projection_value); // take the opposite :)
					
					
					set_witness_set_mp(&W, start_projection,startpt,num_vars); // set the witness point and linear in the input for the lintolin solver.
					
					if (sampler_options->verbose_level>=3) {
						print_point_to_screen_matlab_mp(W.pts_mp[0],"startpt");
						print_comp_mp_matlab(&W.L_mp[0]->coord[0],"initial_projection_value");
						print_comp_mp_matlab(target_projection_value,"target_projection_value");
					}
					
					
					init_witness_set(&Wnew);
					
					
//					lintolin_solver_main(MPType,
//															 W,         // witness_set
//															 n_minusone_randomizer_matrix,
//															 &target_projection, 1,//  the set of linears we will solve at.
//															 &Wnew, // the new data is put here!
//															 solve_options); // already a pointer
					
					

					
					
					multilintolin_solver_main(MPType,
																		W,         // witness_set
																		n_minusone_randomizer_matrix,
																		&target_projection, //  the set of linears we will solve at.
																		&Wnew, // the new data is put here!
																		solve_options); // already a pointer
					
					
					
					
					
					if (sampler_options->verbose_level>=3) 
						print_point_to_screen_matlab_mp(Wnew.pts_mp[0], "new_solution");
					
					
					
					
					
					// check how far away we were from the LEFT interval point
					norm_of_difference(dist_away,
														 Wnew.pts_mp[0], // the current new point
														 V->vertices[left_index].pt_mp);// jj is left, jj+1 is right
					//					mpf_out_str (NULL, 10, 9, dist_away); printf("\n");
					
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
														 V->vertices[right_index].pt_mp);
					
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
					new_indices[sample_counter] = add_vertex( V, temp_vertex);
					sample_counter++;
					
					new_indices[sample_counter] = right_index;
					sample_counter++;

					clear_witness_set(Wnew); // clear the temporary witness set
					
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

//
//
//void generate_new_sampling_pts_mp(sample_data *S_new,
//																	mat_mp n_minusone_randomizer_matrix,
//																	sample_data S_old,
//																	curveDecomp_d C,
//																	vertex_set *V,
//																	witness_set W,
//																	int MPType,
//																	sampler_configuration *sampler_options,
//																	solver_configuration *solve_options)
//{
//	witness_set  Wnew;
//	vec_mp         target_projection;
//	vec_mp          start_projection,startpt;
//	int            ii,jj;
//	comp_mp         temp, temp1, target_projection_value;
//	vec_mp          *previous_points = NULL, *new_points = NULL, previous_projection_values, new_projection_values;
//	int            prev_num_samp, sample_counter;
//	int            *refine_current = NULL, *refine_next = NULL;
//	
//	int num_vars = W.num_variables;
//	
//	init_vec_mp(target_projection,num_vars); target_projection->size = num_vars;
//	vec_cp_mp(target_projection,C.pi_mp); // copy the projection into target_projection
//	
//	init_vec_mp(start_projection,num_vars);  start_projection->size = num_vars; 	
//	vec_cp_mp(start_projection,C.pi_mp); // grab the projection, copy it into start_projection
//	
//	init_vec_mp(new_projection_values,1); new_projection_values->size = 1;
//	init_vec_mp(startpt,num_vars); startpt->size = num_vars;
//	
//	init_mp(temp);  init_mp(temp1); init_mp(target_projection_value);
//
//	S_new->num_edges = S_old.num_edges;
//	
//	
////	S_new->vertices_mp = (vec_mp **)bmalloc(S_new->num_edges * sizeof(vec_mp*));
////	S_new->projection_values_mp = (vec_mp *)bmalloc(S_new->num_edges * sizeof(vec_mp));
////	S_new->refine = (int **)bmalloc(S_new->num_edges * sizeof(int*));
////	S_new->num_pts = (int *) bmalloc(S_new->num_edges * sizeof(int));
//	
//	
//	int num_refinements;
//	for(ii=0;ii<S_old.num_edges;ii++) // for each of the edges
//	{
//		//copy the initial projection values
//		init_vec_mp(previous_projection_values,S_old.num_pts[ii]);
//		previous_projection_values->size = S_old.num_pts[ii]; // inititalize
//		vec_cp_mp(previous_projection_values, S_old.projection_values_mp[ii]); // set projection values
//		
//		
//
//		
//		previous_points = (vec_mp *)bmalloc(S_old.num_pts[ii]*sizeof(vec_mp));
//		for (jj=0; jj<S_old.num_pts[ii]; ++jj) {
//			init_vec_mp(previous_points[jj],1);  previous_points[jj]->size = 1;
//			vec_cp_mp(previous_points[jj], S_old.vertices_mp[ii][jj]);
//		}
//		
//		
//		set_initial_refinement_flags(&num_refinements,&refine_current, S_old.vertices_mp[ii], S_old.num_pts[ii], sampler_options);
//
//		prev_num_samp = S_old.num_pts[ii]; // grab the number of points from the array of integers
//		
//		
//		int pass_number  = 0;//this should be the only place this is reset.
//		while(1) // breaking condition is all samples being less than TOL away from each other (in the infty norm sense).
//		{
//			
//			if (sampler_options->verbose_level>=0){ // print by default
//				printf("edge %d, pass %d, %d refinements\n",ii,pass_number,num_refinements);
//			}
//			
//			new_points = (vec_mp *)bmalloc((prev_num_samp+num_refinements)* sizeof(vec_mp));
//			// will hold the collected data
//			
//			refine_next = (int * )bmalloc((prev_num_samp+num_refinements-1) * sizeof(int)); // refinement flag 
//			memset (refine_next, 0, (prev_num_samp+num_refinements-1) *sizeof(int));
//			
//			change_size_vec_mp(new_projection_values,prev_num_samp+num_refinements);
//			new_projection_values->size=(prev_num_samp+num_refinements); // hold the values of the projections in new_projection_values
//			
//			
//			
//			//copy left point and projection value
//			set_mp(&(new_projection_values->coord[0]),&(previous_projection_values->coord[0]));//leftmost point always the same?
//			init_vec_mp(new_points[0],num_vars); new_points[0]->size = num_vars; // initialize the first of the samples???
//			vec_cp_mp(new_points[0],previous_points[0]);
//			
//			//now we can always set the point to the right as well as the new point
//			
//			
//			sample_counter = 1; // this will be incremented every time we put a point into new_points 
//			// starts at 1 because already committed one.
//			// this should be the only place this counter is reset.
//	
//			
//			if (sampler_options->verbose_level>=2) {
//				printf("the current projection values are:\n");
//				for (jj=0; jj<prev_num_samp; jj++) {
//					print_comp_mp_matlab(&previous_projection_values->coord[jj],"proj");
//				}
//				printf("\n\n");
//			}
//			
//			if (sampler_options->verbose_level>=1) {
//				printf("will refine these at these interval indices:\n");
//				for (jj=0; jj<prev_num_samp-1; jj++) {
//					if (refine_current[jj]) {
//						printf("%d ",jj);
//					}
//				}
//				printf("\n\n");
//			}
//			
//			
//			num_refinements = 0; // reset this counter.  this should be the only place this is reset
//			for(jj=0;jj<prev_num_samp-1;jj++) // for each sample in the previous set
//			{
//				
//				if (sampler_options->verbose_level>=2) {
//					printf("interval %d of %d\n",jj,prev_num_samp-1);
//				}
//				
//				if(refine_current[jj]) //
//				{
//
//					int startpt_index;
//					// set the starting projection and point.
//					if(jj==0)// if on the first sample, go the right
//					{
//						startpt_index = 1; //right!
//					}
//					else //go to the left
//					{
//						startpt_index = jj; // left!
//					}
//
//					
//					
////					vec_cp_mp(startpt,S_old.vertices_mp[ii][1]);
////					set_mp(&(start_projection->coord[0]),&S_old.projection_values_mp[ii]->coord[1]);
////					neg_mp(&(start_projection->coord[0]),&(start_projection->coord[0]));
//					
//					vec_cp_mp(startpt,previous_points[startpt_index]);
//					set_mp(&(start_projection->coord[0]),&(previous_projection_values->coord[startpt_index]));
//					neg_mp(&(start_projection->coord[0]),&(start_projection->coord[0]));
//					
//					
//					estimate_new_projection_value(target_projection_value, // the new value
//																				previous_points[jj], previous_points[jj+1], // two points input
//																				C.pi_mp); // projection (in homogeneous coordinates)
//					
//					neg_mp(&target_projection->coord[0],target_projection_value); // take the opposite :)
//			
//					
//					set_witness_set_mp(&W, start_projection,startpt,num_vars); // set the witness point and linear in the input for the lintolin solver.
//					
//					if (sampler_options->verbose_level>=3) {
//						print_point_to_screen_matlab_mp(W.pts_mp[0],"startpt");
//						print_comp_mp_matlab(&W.L_mp[0]->coord[0],"initial_projection_value");
//						print_comp_mp_matlab(target_projection_value,"target_projection_value");
//						
//					}
//				
//					printf("W.num_linears = %d\n", W.num_linears);
//					init_witness_set(&Wnew);
//					multilintolin_solver_main(MPType,
//																 W,         // witness_set
//																 n_minusone_randomizer_matrix,
//																 &target_projection, //  the set of linears we will solve at.
//																 &Wnew, // the new data is put here!
//																 solve_options);
//					
//					
//
//					
//					
//					if (sampler_options->verbose_level>=3) {
//						print_point_to_screen_matlab_mp(Wnew.pts_mp[0], "new_solution");
//					}
//					mpf_t dist_away; mpf_init(dist_away);
//					
//					
//					// check how far away we were from the left interval point
//					norm_of_difference(dist_away,
//														 Wnew.pts_mp[0], // the current new point
//														 previous_points[jj]);// jj is left, jj+1 is right
////					mpf_out_str (NULL, 10, 9, dist_away); printf("\n");
//					
//					if ( mpf_cmp(dist_away, sampler_options->TOL )>0 ){
//						refine_next[sample_counter-1] = 1; // we started at 1, so to get 0 must -1
//						num_refinements++;
//					}
//					
//					
//					
//					// check how far away we were from the right interval point
//					norm_of_difference(dist_away,
//														 Wnew.pts_mp[0], // the current new point
//														 previous_points[jj+1]);// jj is left, jj+1 is right
//					
//					if (mpf_cmp(dist_away, sampler_options->TOL ) > 0){
//						refine_next[sample_counter] = 1; // we started at 1, so to get 0 must -1
//						num_refinements++;
//					}
//					
//
//					if (sampler_options->verbose_level>=2) {
//						printf("adding sample %d\n",sample_counter);
//					}
//					
//					set_mp(&(new_projection_values->coord[sample_counter]),target_projection_value);
//					init_vec_mp(new_points[sample_counter],num_vars); // add the midpoint we just solved at
//					vec_cp_mp(new_points[sample_counter],Wnew.pts_mp[0]);
//					
//					if (sampler_options->verbose_level>=2) {
//						printf("adding sample %d\n",sample_counter+1);
//					}
//					
//					set_mp(&(new_projection_values->coord[sample_counter+1]),&(previous_projection_values->coord[jj+1])); //copy the right point
//					init_vec_mp(new_points[sample_counter+1],num_vars);
//					vec_cp_mp(new_points[sample_counter+1],previous_points[jj+1]); // copy in the point
//					clear_witness_set(Wnew); // clear the temporary witness set
//					
//					
//					sample_counter++;
//					sample_counter++;
//					
//				}
//				else
//				{
//					if (sampler_options->verbose_level>=2) {
//						printf("adding sample %d\n",sample_counter);
//					}
//					//simply copy in the right point
//					set_mp(&(new_projection_values->coord[sample_counter]),&(previous_projection_values->coord[jj+1]));
//					init_vec_mp(new_points[sample_counter],num_vars);
//					vec_cp_mp(new_points[sample_counter],previous_points[jj+1]);
//					sample_counter++;
//				}
//				
//				
//			}
//			
//			
//			
////			printf("prev_num_samp %d\n",prev_num_samp);
//			//free some memory
//			for (jj=0; jj<prev_num_samp; jj++) {
////				printf("clearing previous_points[%d]\n",jj);
//				clear_vec_mp(previous_points[jj]);
//			}
//			free(previous_points);//frees the malloc above
//			
//			// copy previous_points as a pointer
//			
//			
//
//			
//			if (sampler_options->verbose_level>=1) // print by default
//				printf("\n\n");
//			
//			if( (num_refinements == 0) || (pass_number >= sampler_options->maximum_num_iterations) ) // if have no need for new samples
//			{
//				
//				if (sampler_options->verbose_level>=1) // print by default
//					printf("breaking\nsample_counter = %d\n",sample_counter);
//				
//				S_new->vertices_mp[ii] = (vec_mp *)bmalloc(sample_counter*sizeof(vec_mp));
//				S_new->vertices_mp[ii] = new_points;  // this should be a pointer reassignment
//				S_new->num_pts[ii] = sample_counter;
//				
//				init_vec_mp(S_new->projection_values_mp[ii],new_projection_values->size);
//				vec_cp_mp(S_new->projection_values_mp[ii],new_projection_values);
//				
//				S_new->refine[ii] = refine_next; // will get freed later
//				free(refine_current);
//				refine_current = NULL;
//				break; // BREAKS THE WHILE LOOP
//			}
//			else{
//				refine_current = refine_next; // reassign this pointer
//				vec_cp_mp(previous_projection_values,new_projection_values);
//				previous_points = new_points; // this has gotta fail???  leaking memory?
//				prev_num_samp=sample_counter; // update the number of samples
//				pass_number++;
//			}
//			
//		}//while loop
//		if (sampler_options->verbose_level>=2) {
//			printf("exiting while loop\n");
//		}
//		printf("\n");
//	}  // re: ii (for each edge)
//	
//	
//
//	
//	
//	
//	clear_mp(temp); clear_mp(temp1);
//	clear_vec_mp(start_projection);
//	clear_vec_mp(target_projection);
//}


void read_rand_matrix(char *INfile, mat_mp n_minusone_randomizer_matrix)
{
	int ii,jj,rows,cols,cur_precision;
	FILE *IN= safe_fopen_read(INfile);
	fscanf(IN,"%d %d",&rows,&cols); scanRestOfLine(IN);
	fscanf(IN,"%d",&cur_precision); scanRestOfLine(IN);
	init_mat_mp2(n_minusone_randomizer_matrix,rows,cols,cur_precision);
	n_minusone_randomizer_matrix->rows=rows; n_minusone_randomizer_matrix->cols=cols;
	
	for(ii=0;ii<n_minusone_randomizer_matrix->rows;ii++)
		for(jj=0;jj<n_minusone_randomizer_matrix->cols;jj++)
		{
			mpf_inp_str(n_minusone_randomizer_matrix->entry[ii][jj].r, IN, 10);
			mpf_inp_str(n_minusone_randomizer_matrix->entry[ii][jj].i, IN, 10);
			scanRestOfLine(IN);
		}
	
	fclose(IN);
}

void set_witness_set_d(witness_set *W, vec_d new_linear,vec_d pts,int num_vars)
{
	
	W->num_variables = num_vars;
	
	W->num_pts=1;
	W->pts_mp=(point_mp *)bmalloc(sizeof(point_mp)); // apparently you can only pass in a single point to copy in.
	
	
	W->num_pts=1;
	W->pts_d=(point_d *)bmalloc(sizeof(point_d)); // apparently you can only pass in a single point to copy in.
	
	
	//initialize the memory
	init_point_d(W->pts_d[0],num_vars); W->pts_d[0]->size = num_vars;
	init_point_mp(W->pts_mp[0],num_vars); W->pts_mp[0]->size = num_vars;
	
	
	point_cp_d(W->pts_d[0],pts);
	point_d_to_mp(W->pts_mp[0],pts);
	
	W->num_linears = 1;
	W->L_d = (vec_d *)bmalloc(sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(sizeof(vec_mp));
	
	init_vec_d( W->L_d[0],   num_vars); W->L_d[0]->size = num_vars;
	init_vec_mp(W->L_mp[0],num_vars); W->L_mp[0]->size = num_vars;
	
	vec_cp_d(W->L_d[0],new_linear);
	vec_d_to_mp(W->L_mp[0],new_linear);
	
}

void set_witness_set_mp(witness_set *W, vec_mp new_linear,vec_mp pts,int num_vars)
{
	
	
	W->num_variables = num_vars;
	
	
	W->num_pts=1;
	W->pts_mp=(point_mp *)bmalloc(sizeof(point_mp)); // apparently you can only pass in a single point to copy in.
	
	
	W->num_pts=1;
	W->pts_d=(point_d *)bmalloc(sizeof(point_d)); // apparently you can only pass in a single point to copy in.
	
	
	//initialize the memory
	init_point_d(W->pts_d[0],num_vars); W->pts_d[0]->size = num_vars;
	init_point_mp(W->pts_mp[0],num_vars); W->pts_mp[0]->size = num_vars;
	
	
	point_cp_mp(W->pts_mp[0],pts);
	point_mp_to_d(W->pts_d[0],pts);
	
	W->num_linears = 1;
	W->L_d = (vec_d *)bmalloc(sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(sizeof(vec_mp));
	
	init_vec_d( W->L_d[0],   num_vars); W->L_d[0]->size = num_vars;
	init_vec_mp(W->L_mp[0],num_vars); W->L_mp[0]->size = num_vars;
	
	vec_cp_mp(W->L_mp[0],new_linear);
	vec_mp_to_d(W->L_d[0],new_linear);
	
}


void  output_sampling_data(sample_data S, vertex_set V,
													 char *samplingName,int num_vars,int MPType)
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


int  curve_sampler_startup(int argC, char *args[], char directoryName[],
													 char **inputName,
													 char **witnessSetName,
													 char **RandMatName,
													 char **samplingNamenew,
													 curveDecomp_d *C,
													 vertex_set *V,
													 int *num_vars,
													 int MPType)
/***************************************************************\
 * USAGE:    curve_sampler_startup structure and inputname
 * ARGUMENTS:                                                    *
 * RETURN VALUES:           an integer                           *
 * NOTES:                                                        *
 \***************************************************************/
{
	int strLength;
	char tmp_file[1024];

	int successful_startup = 1;
	//read in the gross information from the summary file.
	sprintf(tmp_file,  "%s/C.curve", directoryName);
	setup_curve(C,tmp_file,MPType,inputName,directoryName);
	
	//setup E structure from E.edge
	sprintf(tmp_file,  "%s/E.edge", directoryName);
	curve_setup_edges(C,tmp_file);
	
	//setup V structure from V.vertex
	sprintf(tmp_file,  "%s/V.vertex", directoryName);
	setup_vertices(V,tmp_file);
	
	strLength = 1 + snprintf(NULL, 0, "%s/witness_set", directoryName);
	*witnessSetName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*witnessSetName, "%s/witness_set", directoryName);
	
	
	strLength = 1 + snprintf(NULL, 0, "%s/Rand_Matrix", directoryName);
	*RandMatName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*RandMatName, "%s/Rand_Matrix", directoryName);
	
	strLength = 1 + snprintf(NULL, 0, "%s/samp.dat", directoryName);
	*samplingNamenew = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*samplingNamenew, "%s/samp.dat", directoryName);
	if(C->num_edges==0)
	{
		printf("sampler will not generate sampling data since there are no edges\n");
		successful_startup = 0;
	}
	
	return successful_startup;
}


int get_dir_mptype(char ** Dir_Name, int * MPType){

	FILE *IN;
	int strLength;
	//setup the name of directory
	IN = safe_fopen_read("Dir_Name");
	fscanf(IN, "%d\n", &strLength);
	*Dir_Name = (char *)br_malloc(strLength * sizeof(char));
	fgets(*Dir_Name, strLength, IN);
	fscanf(IN, "%d\n", MPType);
	fclose(IN);
	
	return *MPType;
}


int curve_setup_edges(curveDecomp_d *C,char *INfile)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	
	
	fscanf(IN, "%d\n", &C->num_edges);
	
	
	C->edges=(edge*) bmalloc(C->num_edges*sizeof(edge));
	int ii;
	for(ii=0;ii<C->num_edges;ii++) {
		fscanf(IN,"%d %d %d",&C->edges[ii].left, &C->edges[ii].midpt, &C->edges[ii].right); scanRestOfLine(IN);
	}
	
	
	printf("%d\n", C->num_edges);
	for(ii=0;ii<C->num_edges;ii++) {
		printf("%d %d %d\n",C->edges[ii].left, C->edges[ii].midpt, C->edges[ii].right);
	}
	
	fclose(IN);
	return C->num_edges;
}


int setup_vertices(vertex_set *V,
									 char *INfile)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices,ii,jj;
	int num_vars;
	fscanf(IN, "%d %d\n\n", &num_vertices, &num_vars);
	
	printf("%d variables\n",num_vars);
	vertex temp_vertex;
	init_vertex(&temp_vertex);
	change_size_vec_mp(temp_vertex.pt_mp,num_vars);temp_vertex.pt_mp->size = num_vars;
	
	for(ii=0;ii<num_vertices;ii++)
	{
		for(jj=0;jj<num_vars;jj++)
		{
			mpf_inp_str(temp_vertex.pt_mp->coord[jj].r, IN, 10);
			mpf_inp_str(temp_vertex.pt_mp->coord[jj].i, IN, 10);
		}
		
		mpf_inp_str(temp_vertex.projVal_mp->r, IN, 10);
		mpf_inp_str(temp_vertex.projVal_mp->i, IN, 10);
		
		fscanf(IN,"%d\n",&temp_vertex.type);

		
		add_vertex(V, temp_vertex);		
	}
	
	
	fclose(IN);
	
	if (V->num_vertices!=num_vertices) {
		printf("parity error in num_vertices.\n\texpected: %d\tactual: %d\n",V->num_vertices,num_vertices);
		exit(25943);
	}
	
	return num_vertices;
}






int setup_curve(curveDecomp_d *C,char *INfile, int MPType, char **inputName, char *directoryName)
//setup the vertex structure
{
	
	int ii,strLength;
	char *input_deflated_Name=NULL;
	
	
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices = 0;
	fscanf(IN, "%d\n", &strLength);
	
	input_deflated_Name = (char *)bmalloc((strLength+1) * sizeof(char));
	fgets(input_deflated_Name, (strLength+1), IN);
	
	strLength = 1 + snprintf(NULL, 0, "%s/%s", directoryName,input_deflated_Name);
	*inputName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*inputName,  "%s/%s", directoryName,input_deflated_Name);
	
	
	
	
	fscanf(IN, "%d %d", &C->num_variables, &C->num_edges);  scanRestOfLine(IN);
	
	fscanf(IN, "%d %d %d %d %d", &C->num_V0, &C->num_V1, &C->num_midpts, &C->num_new, &C->num_isolated);  scanRestOfLine(IN);
	
	C->V0_indices = (int *)bmalloc(C->num_V0*sizeof(int));
	C->V1_indices = (int *)bmalloc(C->num_V1*sizeof(int));
	C->midpt_indices = (int *)bmalloc(C->num_midpts*sizeof(int));
	C->new_indices = (int *)bmalloc(C->num_new*sizeof(int));

	for (ii=0; ii<C->num_V0; ++ii) {
		fscanf(IN,"%d\n",&C->V0_indices[ii]);
	}
	
	for (ii=0; ii<C->num_V1; ++ii) {
		fscanf(IN,"%d\n",&C->V1_indices[ii]);
	}
	
	for (ii=0; ii<C->num_midpts; ++ii) {
		fscanf(IN,"%d\n",&C->midpt_indices[ii]);
	}
	
	for (ii=0; ii<C->num_new; ++ii) {
		fscanf(IN,"%d\n",&C->new_indices[ii]);
	}
	
	
	init_point_mp(C->pi_mp,C->num_variables); C->pi_mp->size=C->num_variables;

	for(ii=0;ii<C->num_variables;ii++)
	{
		mpf_inp_str(C->pi_mp->coord[ii].r, IN, 10);
		mpf_inp_str(C->pi_mp->coord[ii].i, IN, 10);
	}
	


	
	fclose(IN);
	return num_vertices;
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
	
	dehomogenize_mp(&dehom_left,left);
	dehomogenize_mp(&dehom_right,right);
	
	
	
//	print_point_to_screen_matlab_mp(dehom_left,"left_pt");
//	print_point_to_screen_matlab_mp(dehom_right,"right_pt");
//	print_point_to_screen_matlab_mp(pi,"pi");
	
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




int  set_initial_sample_data(sample_data *S, curveDecomp_d C, vertex_set V,
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
																 sample_data S, vertex_set *V,
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
		
		norm_of_difference(dist_away, V->vertices[S.sample_indices[current_edge][ii]].pt_mp,
																	V->vertices[S.sample_indices[current_edge][ii+1]].pt_mp); // get the distance between the two adjacent points.
		if ( mpf_cmp(dist_away, sampler_options->TOL)>0 ){
			(*refine_flags)[ii] = 1;
			(*num_refinements)++;
		}
	}
//	printf("done setting refinement flags\n");
	mpf_clear(dist_away);
}




