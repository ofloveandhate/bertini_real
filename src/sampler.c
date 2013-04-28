#include "sampler.h"


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
	curveDecomp_d C;  //new data type; stores vertices, edges, etc.
	witness_set W;
	sample_d   S_old,S_new;

	sampler_splash_screen();
	
	sampler_configuration sampler_options;  init_sampler_config(&sampler_options);
	
	sampler_parse_commandline(argC, args, &sampler_options);
	
	////
	//  begin the actual program
	////
	
	if (setup_curveDecomp(argC, args, &inputName, &witnessSetName,&RandMatName,&samplingNamenew,&C,&num_vars))
		return 1;
	
	

	srand(time(NULL));
	
	int MPType;


	
	mat_mp n_minusone_randomizer_matrix;
	
	//end parser-bertini essentials
	
	parse_input_file(inputName, &MPType);

	
	// set up the solver configuration
	solver_configuration solve_options;  init_solver_config(&solve_options);
	get_tracker_config(&solve_options,MPType);
	setupPreProcData("preproc_data", &solve_options.PPD);
	
	initMP(solve_options.T.Precision);
	
	num_vars = get_num_vars_PPD(solve_options.PPD);		
	
	
//	printf("parsing witness set\n");
	witnessSetParse(&W, witnessSetName,num_vars);
	W.num_var_gps = solve_options.PPD.num_var_gp;
	W.MPType = MPType;
	get_variable_names(&W);

	
//	printf("getting randomizer matrix\n");
	read_rand_matrix(RandMatName, n_minusone_randomizer_matrix);
	
//	printf("loading sampling data\n");
	Load_sampling_data(&S_old,C,num_vars,MPType);
	
	solve_options.T.ratioTol = 1;
	
	/////////
	////////
	//////
	////
	//
	//  Generate new sampling data
	//

	
//	printf("generate_new_sampling_pts\n");
	generate_new_sampling_pts(&S_new,
														n_minusone_randomizer_matrix,
														S_old,
														C,
														W,
														MPType,
														&sampler_options,
														&solve_options);
	
	

	//
	//   done with the main call
	////
	/////
	///////
	////////
	
	
	
	//output
	output_sampling_data(S_new,samplingNamenew,num_vars,MPType);
	
	
	
	
	
	// clear memory
	free(inputName);
	free(witnessSetName);
	
	clear_mat_mp(n_minusone_randomizer_matrix);
	clear_witness_set(W);
	
	
	clear_sample_d(&S_new, MPType);
	
//	clear_sample_d(&S_old, MPType);
	
	clear_curveDecomp_d(&C, MPType);
	clearMP();
	return 0;
}







void generate_new_sampling_pts(sample_d *S_new,
															 mat_mp n_minusone_randomizer_matrix,
															 sample_d S_old,
															 curveDecomp_d C,
															 witness_set W,
															 int MPType ,
															 sampler_configuration *sampler_options,
															 solver_configuration *solve_options)
{
	
	
	solve_options->verbose_level = sampler_options->verbose_level;
	
	if(MPType==0)
		generate_new_sampling_pts_d(S_new, n_minusone_randomizer_matrix, S_old, C, W, MPType, sampler_options,solve_options);
	else
		generate_new_sampling_pts_mp(S_new, n_minusone_randomizer_matrix, S_old, C, W, MPType, sampler_options, solve_options);
	return;
}


void generate_new_sampling_pts_d(sample_d *S_new,
																 mat_mp n_minusone_randomizer_matrix,
																 sample_d S_old,
																 curveDecomp_d C
																 ,witness_set W,
																 int MPType,
																 sampler_configuration *sampler_options,
																 solver_configuration *solve_options)
{
	witness_set		Wnew;
	vec_mp         target_projection;
	vec_d          start_projection,startpt,mid_pt;
	int            ii,jj,k;
	comp_d         pi_end,temp,temp1;
	vec_d          *previous_points, *new_points, previous_projection_values, new_projection_values;
	int            prev_num_samp, sample_counter;
	int            *refine_current, *refine_next;
	double         max_norm, TOL=1e-1;
	
	int num_vars =  W.num_variables;

	
	init_vec_mp(target_projection,num_vars); target_projection->size = num_vars;
	
	init_vec_d(start_projection,num_vars);  start_projection->size = num_vars;
	
	init_vec_d(startpt,num_vars);
	startpt->size = num_vars;
	
	init_vec_d(mid_pt,num_vars);
	mid_pt->size = num_vars;
	
	
	S_new->num_edges = S_old.num_edges;
	S_new->vertices = (vec_d **)bmalloc(S_new->num_edges * sizeof(vec_d*));
	S_new->projection_values = (vec_d *)bmalloc(S_new->num_edges * sizeof(vec_d));
	S_new->refine = (int **)bmalloc(S_new->num_edges * sizeof(int*));
	S_new->num_pts = (int *) bmalloc(S_new->num_edges * sizeof(int));
	
	
	
	
	for(ii=0;ii<S_old.num_edges;ii++) // for each edge, refine it.
	{
		previous_points = S_old.vertices[ii]; // copy pointer
		prev_num_samp = S_old.num_pts[ii]; // copy number of points
		
		
		init_vec_d(previous_projection_values,S_old.num_pts[ii]);
		previous_projection_values->size = S_old.num_pts[ii];
		
		vec_cp_d(previous_projection_values, S_old.projection_values[ii]);
		refine_current = S_old.refine[ii];
		
		
		while(1) {
			new_points = (vec_d *)bmalloc(2*prev_num_samp * sizeof(vec_d)); // why the 2?
			init_vec_d(new_projection_values,2*prev_num_samp); // what is projection_values?
			//new_projection_values->size=2*prev_num_samp;
			refine_next = (int * )bmalloc(2*prev_num_samp * sizeof(int));
			//copy left point
			set_d(&(new_projection_values->coord[0]),&(previous_projection_values->coord[0]));
			sample_counter = 1;
			init_vec_d(new_points[0],num_vars);
			vec_cp_d(new_points[0],previous_points[0]);
			for(jj=0;jj<prev_num_samp-1;jj++)
			{
				if(refine_current[jj])
				{
//TODO: resolve this discrepancy
					vec_cp_d(start_projection,C.pi_d);
					point_d_to_mp(target_projection,start_projection);
					if(jj==0)
					{
						vec_cp_d(startpt,previous_points[jj+1]);
						sub_d(&(start_projection->coord[0]),&(start_projection->coord[0]),&(previous_projection_values->coord[jj+1]));
					}
					else
					{
						vec_cp_d(startpt,previous_points[jj]);
						sub_d(&(start_projection->coord[0]),&(start_projection->coord[0]),&(previous_projection_values->coord[jj]));
					}
					//tracking end point pi_end = \pi (si+s_{ii+1})/2
					//w=2w1w2/(w1+w2)
					add_d(&(mid_pt->coord[0]),&(previous_points[jj]->coord[0]),
								&(previous_points[jj+1]->coord[0]));
					mul_d(temp,&(previous_points[jj]->coord[0]),
								&(previous_points[jj+1]->coord[0]));
					div_d(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					temp->r=2.0; temp->i=0.0;
					mul_d(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mul_d(pi_end,&(mid_pt->coord[0]),&(C.pi_d->coord[0]));
					
					for(k=1;k<num_vars;k++)
					{
						div_d(&(mid_pt->coord[k]),&(previous_points[jj]->coord[k]),
									&(previous_points[jj]->coord[0]));
						div_d(temp,&(previous_points[jj+1]->coord[k]),
									&(previous_points[jj+1]->coord[0]));
						add_d(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						temp->r=0.5; temp->i=0.0;
						mul_d(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						mul_d(&(mid_pt->coord[k]),&(mid_pt->coord[0]),&(mid_pt->coord[k]));
						mul_d(temp,&(mid_pt->coord[k]),&(C.pi_d->coord[k]));
						add_d(pi_end,pi_end,temp);
					}
					div_d(pi_end,pi_end,&(mid_pt->coord[0]));
					sub_d(temp,&C.pi_d->coord[0],pi_end);
					d_to_mp(&(target_projection->coord[0]),temp);
					
					// copy in the data to the source witness set?
					set_witness_set_d(&W, start_projection,startpt,num_vars);
					//printVec_d(stdout,0,mid_pt);
					//print_d(stdout,0,pi_end);
					//printf("start_projection\n");
					//printVec_d(stdout,0,start_projection);
					//printf("newlinears\n");
					//printVec_mp(stdout,0,target_projection);
					
					lin_to_lin_solver_main(MPType,
																 W,         // witness_set
																 n_minusone_randomizer_matrix,
																 &target_projection, //  the set of linears we will solve at.
																 1, // the number of new linears.
																 &Wnew,// the new data is put here!
																 solve_options);
					//printVec_d(stdout,0,mid_pt);
					
					set_zero_d(&(new_projection_values->coord[sample_counter]));
					init_vec_d(new_points[sample_counter],num_vars);
					vec_cp_d(new_points[sample_counter],Wnew.pts_d[0]);
					
					//printVec_d(stdout,0,Wnew.pts_d[0]);
					//printVec_d(stdout,0,mid_pt);
					max_norm=0.0;
					for(k=1;k<num_vars;k++)
					{
						div_d(temp1,&(new_points[sample_counter]->coord[k]),
									&(new_points[sample_counter]->coord[0]));
						div_d(temp,&(mid_pt->coord[k]),&(mid_pt->coord[0]));
						sub_d(temp,temp1,temp);
						if(max_norm<d_abs_d(temp))
							max_norm = d_abs_d(temp);
					}
					//printf("pi_end");
					//print_d(stdout,0,pi_end);
					set_d(&(new_projection_values->coord[sample_counter]),pi_end);
					
					if(max_norm<TOL)
					{
						refine_next[sample_counter-1] = 0;
						refine_next[sample_counter++] = 0;
					}
					else
					{
						refine_next[sample_counter-1] = 1;
						refine_next[sample_counter++] = 1;
					}
					set_d(&(new_projection_values->coord[sample_counter]),&(previous_projection_values->coord[jj+1]));
					init_vec_d(new_points[sample_counter],num_vars);
					vec_cp_d(new_points[sample_counter],previous_points[jj+1]);
					sample_counter++;
					
				}
				else
				{
					refine_next[sample_counter-1] = 0;
					set_d(&(new_projection_values->coord[sample_counter]),&(previous_projection_values->coord[jj+1]));
					init_vec_d(new_points[sample_counter],num_vars);
					vec_cp_d(new_points[sample_counter],previous_points[jj+1]);
					sample_counter++;
				}
			}
			new_points = (vec_d *)brealloc(new_points,sample_counter*sizeof(vec_d));
			new_projection_values->size=sample_counter;
			refine_next = (int * )brealloc(refine_next,sample_counter * sizeof(int));
			
			previous_points = new_points;
			vec_cp_d(previous_projection_values,new_projection_values);
			clear_vec_d(new_projection_values);
			refine_current = refine_next;
			if(sample_counter == prev_num_samp)
			{
				S_new->vertices[ii] = previous_points;
				S_new->num_pts[ii] = prev_num_samp;
				init_vec_d(S_new->projection_values[ii],prev_num_samp);
				vec_cp_d(S_new->projection_values[ii],previous_projection_values);
				clear_vec_d(previous_projection_values);
				S_new->refine[ii] = refine_current;
				break;
			}
			else
				prev_num_samp=sample_counter;
		}
	}
	
	clear_witness_set(Wnew);
	clear_vec_d(start_projection);
	clear_vec_d(target_projection);
}



void generate_new_sampling_pts_mp(sample_d *S_new,
																	mat_mp n_minusone_randomizer_matrix,
																	sample_d S_old,
																	curveDecomp_d C,
																	witness_set W,
																	int MPType,
																	sampler_configuration *sampler_options,
																	solver_configuration *solve_options)
{
	witness_set  Wnew;
	vec_mp         target_projection;
	vec_mp          start_projection,startpt;
	int            ii,jj;
	comp_mp         temp, temp1, target_projection_value;
	vec_mp          *previous_points, *new_points, previous_projection_values, new_projection_values;
	int            prev_num_samp, sample_counter;
	int            *refine_current, *refine_next;
	
	int num_vars = W.num_variables;
	
	init_vec_mp(target_projection,num_vars); target_projection->size = num_vars;
	vec_cp_mp(target_projection,C.pi); // copy the projection into target_projection
	
	init_vec_mp(start_projection,num_vars);  start_projection->size = num_vars; 	
	vec_cp_mp(start_projection,C.pi); // grab the projection, copy it into start_projection
	
	init_vec_mp(new_projection_values,1); new_projection_values->size = 1;
	init_vec_mp(startpt,num_vars); startpt->size = num_vars;
	
	init_mp(temp);  init_mp(temp1); init_mp(target_projection_value);

	S_new->num_edges = S_old.num_edges;
	S_new->vertices_mp = (vec_mp **)bmalloc(S_new->num_edges * sizeof(vec_mp*));
	S_new->projection_values_mp = (vec_mp *)bmalloc(S_new->num_edges * sizeof(vec_mp));
	S_new->refine = (int **)bmalloc(S_new->num_edges * sizeof(int*));
	S_new->num_pts = (int *) bmalloc(S_new->num_edges * sizeof(int));
	
	
	int num_refinements;
	for(ii=0;ii<S_old.num_edges;ii++) // for each of the edges
	{
		//copy the initial projection values
		init_vec_mp(previous_projection_values,S_old.num_pts[ii]); previous_projection_values->size = S_old.num_pts[ii]; // inititalize
		
		num_refinements = 2; // one for each, right and left
		
		vec_cp_mp(previous_projection_values, S_old.projection_values_mp[ii]); // set projection values
		refine_current = S_old.refine[ii]; // change the pointer of refine_current to be S_old.refine[ii]
		previous_points = S_old.vertices_mp[ii]; // does this copy a pointer?
		prev_num_samp = S_old.num_pts[ii]; // grab the number of points from the array of integers
		
		int pass_number  = 0;
		while(1) // breaking condition is all samples being less than TOL away from each other (in the infty norm sense).
		{
			
			if (sampler_options->verbose_level>=0){ // print by default
				printf("edge %d,\tpass %d,\t%d refinements\n",ii,pass_number,num_refinements);
			}
			
			new_points = (vec_mp *)bmalloc((prev_num_samp+num_refinements)* sizeof(vec_mp));
			// will hold the collected data
			
			refine_next = (int * )bmalloc((prev_num_samp+num_refinements-1) * sizeof(int)); // refinement flag 
			memset (refine_next, 0, (prev_num_samp+num_refinements-1) *sizeof(int));
			
			change_size_vec_mp(new_projection_values,prev_num_samp+num_refinements);
			new_projection_values->size=(prev_num_samp+num_refinements); // hold the values of the projections in new_projection_values
			
			
			
			//copy left point and projection value
			set_mp(&(new_projection_values->coord[0]),&(previous_projection_values->coord[0]));//leftmost point always the same?
			init_vec_mp(new_points[0],num_vars); new_points[0]->size = num_vars; // initialize the first of the samples???
			vec_cp_mp(new_points[0],previous_points[0]);
			
			//now we can always set the point to the right as well as the new point
			
			
			sample_counter = 1; // this will be incremented every time we put a point into new_points 
			// starts at 1 because already committed one.
			// this should be the only place this counter is reset.
	
			
			if (sampler_options->verbose_level>=2) {
				printf("the current projection values are:\n");
				for (jj=0; jj<prev_num_samp; jj++) {
					print_comp_mp_matlab(&previous_projection_values->coord[jj],"proj");
				}
				printf("\n\n");
			}
			
			if (sampler_options->verbose_level>=1) {
				printf("will refine these at these interval indices:\n");
				for (jj=0; jj<prev_num_samp-1; jj++) {
					if (refine_current[jj]) {
						printf("%d ",jj);
					}
				}
				printf("\n\n");
			}
			
			
			num_refinements = 0; // reset this counter.  this should be the only place this is reset
			for(jj=0;jj<prev_num_samp-1;jj++) // for each sample in the previous set
			{
				
				if (sampler_options->verbose_level>=2) {
					printf("interval %d of %d\n",jj,prev_num_samp-1);
				}
				
				if(refine_current[jj]) //
				{

					int index_of_current_startpoint;
					// set the starting projection and point.
					if(jj==0)// if on the first sample, go the right
					{
						index_of_current_startpoint = 1; //right!
					}
					else //go to the left
					{
						index_of_current_startpoint = jj; // left!
					}

					vec_cp_mp(startpt,previous_points[index_of_current_startpoint]);
					set_mp(&(start_projection->coord[0]),&(previous_projection_values->coord[index_of_current_startpoint]));
					neg_mp(&(start_projection->coord[0]),&(start_projection->coord[0]));
					
					
					estimate_new_projection_value(target_projection_value, // the new value
																				previous_points[jj], previous_points[jj+1], // two points input
																				C.pi); // projection (in homogeneous coordinates)
					
					neg_mp(&target_projection->coord[0],target_projection_value); // take the opposite :)
			
					
					set_witness_set_mp(&W, start_projection,startpt,num_vars); // set the witness point and linear in the input for the lintolin solver.
					
					if (sampler_options->verbose_level>=3) {
						print_comp_mp_matlab(&W.L_mp[0]->coord[0],"initial_projection_value");
						print_comp_mp_matlab(target_projection_value,"target_projection_value");
						print_point_to_screen_matlab_mp(W.pts_mp[0],"startpt");
					}
				
					
					init_witness_set_d(&Wnew);
					lin_to_lin_solver_main(MPType,
																 W,         // witness_set
																 n_minusone_randomizer_matrix,
																 &target_projection, //  the set of linears we will solve at.
																 1, // the number of new linears.
																 &Wnew, // the new data is put here!
																 solve_options);
					
					

					
					
					if (sampler_options->verbose_level>=3) {
						print_point_to_screen_matlab_mp(Wnew.pts_mp[0], "new_solution");
					}
					mpf_t dist_away; mpf_init(dist_away);
					
					
					// check how far away we were from the left interval point
					norm_of_difference(dist_away,
														 Wnew.pts_mp[0], // the current new point
														 previous_points[jj]);// jj is left, jj+1 is right
//					mpf_out_str (NULL, 10, 9, dist_away); printf("\n");
					
					if ( mpf_cmp(dist_away, sampler_options->TOL )>0 ){
						refine_next[sample_counter-1] = 1; // we started at 1, so to get 0 must -1
						num_refinements++;
					}
					
					
					
					// check how far away we were from the right interval point
					norm_of_difference(dist_away,
														 Wnew.pts_mp[0], // the current new point
														 previous_points[jj+1]);// jj is left, jj+1 is right
					
					if (mpf_cmp(dist_away, sampler_options->TOL ) > 0){
						refine_next[sample_counter] = 1; // we started at 1, so to get 0 must -1
						num_refinements++;
					}
					

					if (sampler_options->verbose_level>=2) {
						printf("adding sample %d\n",sample_counter);
					}
					
					set_mp(&(new_projection_values->coord[sample_counter]),target_projection_value);
					init_vec_mp(new_points[sample_counter],num_vars); // add the midpoint we just solved at
					vec_cp_mp(new_points[sample_counter],Wnew.pts_mp[0]);
					
					if (sampler_options->verbose_level>=2) {
						printf("adding sample %d\n",sample_counter+1);
					}
					
					set_mp(&(new_projection_values->coord[sample_counter+1]),&(previous_projection_values->coord[jj+1])); //copy the right point
					init_vec_mp(new_points[sample_counter+1],num_vars);
					vec_cp_mp(new_points[sample_counter+1],previous_points[jj+1]); // copy in the point
					clear_witness_set(Wnew); // clear the temporary witness set
					
					
					sample_counter++;
					sample_counter++;
					
				}
				else
				{
					if (sampler_options->verbose_level>=2) {
						printf("adding sample %d\n",sample_counter);
					}
					//simply copy in the right point
					set_mp(&(new_projection_values->coord[sample_counter]),&(previous_projection_values->coord[jj+1]));
					init_vec_mp(new_points[sample_counter],num_vars);
					vec_cp_mp(new_points[sample_counter],previous_points[jj+1]);
					sample_counter++;
				}
				
				
			}
			
			
			
//			printf("prev_num_samp %d\n",prev_num_samp);
			//free some memory
			for (jj=0; jj<prev_num_samp; jj++) {
//				printf("clearing previous_points[%d]\n",jj);
				clear_vec_mp(previous_points[jj]);
			}
			free(previous_points);//frees the malloc above
			
			// copy previous_points as a pointer
			
			

			
			if (sampler_options->verbose_level>=0) // print by default
				printf("\n\n");
			
			if( (num_refinements == 0) || (pass_number >= sampler_options->maximum_num_iterations) ) // if have no need for new samples
			{
				if (sampler_options->verbose_level>=1) // print by default
					printf("breaking\nsample_counter = %d\n",sample_counter);
				
				S_new->vertices_mp[ii] = (vec_mp *)bmalloc(sample_counter*sizeof(vec_mp));
				S_new->vertices_mp[ii] = new_points;  // this should be a pointer reassignment
				S_new->num_pts[ii] = sample_counter;
				
				init_vec_mp(S_new->projection_values_mp[ii],new_projection_values->size);
				vec_cp_mp(S_new->projection_values_mp[ii],new_projection_values);
				
				S_new->refine[ii] = refine_next; // will get freed later
				break; // BREAKS THE WHILE LOOP
			}
			else{
				refine_current = refine_next; // reassign this pointer
				vec_cp_mp(previous_projection_values,new_projection_values);
				previous_points = new_points; // this has gotta fail???  leaking memory?
				prev_num_samp=sample_counter; // update the number of samples
				pass_number++;
			}
			
		}//while loop
		if (sampler_options->verbose_level>=2) {
			printf("exiting while loop\n");
		}
		
	}
	
	

	
	
	
	clear_mp(temp); clear_mp(temp1);
	clear_vec_mp(start_projection);
	clear_vec_mp(target_projection);
}


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
	W->L = (vec_d *)bmalloc(sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(sizeof(vec_mp));
	
	init_vec_d( W->L[0],   num_vars); W->L[0]->size = num_vars;
	init_vec_mp(W->L_mp[0],num_vars); W->L_mp[0]->size = num_vars;
	
	vec_cp_d(W->L[0],new_linear);
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
	W->L = (vec_d *)bmalloc(sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(sizeof(vec_mp));
	
	init_vec_d( W->L[0],   num_vars); W->L[0]->size = num_vars;
	init_vec_mp(W->L_mp[0],num_vars); W->L_mp[0]->size = num_vars;
	
	vec_cp_mp(W->L_mp[0],new_linear);
	vec_mp_to_d(W->L[0],new_linear);
	
}


void  output_sampling_data(sample_d S,char *samplingName,int num_vars,int MPType)
{
	FILE *OUT =  fopen(samplingName, "w");
	int ii,jj,k;
	// output the number of vertices
	fprintf(OUT,"%d\n\n",S.num_edges);
	for(ii=0;ii<S.num_edges;ii++)
	{
		fprintf(OUT,"%d\n\n",S.num_pts[ii]);
		for(jj=0;jj<S.num_pts[ii];jj++)
		{
			if(MPType==0)
				print_d(OUT,0,&(S.projection_values[ii]->coord[jj]));
			else
				print_mp(OUT,0,&(S.projection_values_mp[ii]->coord[jj]));
			fprintf(OUT,"\n");
			for(k=0;k<num_vars;k++)
			{
				if(MPType==0)
					print_d(OUT, 0, &(S.vertices[ii][jj]->coord[k]));
				else
					print_mp(OUT, 0, &(S.vertices_mp[ii][jj]->coord[k]));
				
				fprintf(OUT,"\n");
			}
			fprintf(OUT,"\n");
		}
	}
}

int  Load_sampling_data(sample_d *S, curveDecomp_d C,int num_vars,int MPType)
{
	int ii,jj,index;
	

	comp_mp temp_mp;
	S->num_edges = C.num_edges;
	S->num_pts = (int *)bmalloc(C.num_edges * sizeof(int));
	
	for(ii=0;ii<C.num_edges;ii++)
		S->num_pts[ii]=3; // left, right, mid
	
	if(MPType==0) {
		S->vertices = (vec_d **)bmalloc(C.num_edges * sizeof(vec_d*));
		S->refine = (int **)bmalloc(C.num_edges * sizeof(int*));
		S->projection_values = (vec_d *)bmalloc(C.num_edges * sizeof(vec_d));
	}
	else {
		init_mp(temp_mp);
		S->vertices_mp = (vec_mp **)bmalloc(C.num_edges * sizeof(vec_mp*));
		S->refine = (int **)bmalloc(C.num_edges * sizeof(int*));
		S->projection_values_mp = (vec_mp *)bmalloc(C.num_edges * sizeof(vec_mp));
	}
	
	for(ii=0;ii<S->num_edges;ii++)
	{
		if(MPType==0)
		{
			S->vertices[ii] = (vec_d *)bmalloc(S->num_pts[ii] * sizeof(vec_d));
			S->refine[ii] = (int *)bmalloc(S->num_pts[ii] * sizeof(int));
			
			//three points: left mid right
			for(jj=0;jj<S->num_pts[ii];jj++) {
				init_vec_d(S->vertices[ii][jj],num_vars);
				S->refine[ii][jj]=1;  //initially set for subsampling refinement
			}
			
			//initialize the thing which holds the projection values
			init_vec_d(S->projection_values[ii],S->num_pts[ii]);
		
			
			//left point
			index = C.edges[ii].left;
			vec_cp_d(S->vertices[ii][0],C.vertices[index].pt);
			projection_value_homogeneous_input_d(&(S->projection_values[ii]->coord[0]), // destination
																					 C.vertices[index].pt, C.pi_d); //source
			
			
			//mid point
			index = C.edges[ii].midpt;
			vec_cp_d(S->vertices[ii][1],C.vertices[index].pt);
			projection_value_homogeneous_input_d(&(S->projection_values[ii]->coord[1]), // destination
																				 C.vertices[index].pt, C.pi_d); //source
			
			
			
			//right point
			index = C.edges[ii].right;
			vec_cp_d(S->vertices[ii][2],C.vertices[index].pt);
			projection_value_homogeneous_input_d(&(S->projection_values[ii]->coord[2]), // destination
																					 C.vertices[index].pt, C.pi_d); //source
			
		}
		else
		{

			
			S->vertices_mp[ii] = (vec_mp *)bmalloc(S->num_pts[ii] * sizeof(vec_mp));
			S->refine[ii] = (int *)bmalloc(S->num_pts[ii] * sizeof(int));
			
			for(jj=0;jj<S->num_pts[ii];jj++)
			{
				init_vec_mp(S->vertices_mp[ii][jj],num_vars);//three points: left mid right
				S->refine[ii][jj]=1; // initially set the refine flag
			}
			
			init_vec_mp(S->projection_values_mp[ii],S->num_pts[ii]); S->projection_values_mp[ii]->size=S->num_pts[ii];
			
			
			//left point
			index = C.edges[ii].left;
			vec_cp_mp(S->vertices_mp[ii][0],C.vertices[index].pt_mp);
			set_mp(&(S->projection_values_mp[ii]->coord[0]),C.vertices[index].projVal_mp);
//			projection_value_homogeneous_input(&(S->projection_values_mp[ii]->coord[0]), // destination
//																					 C.vertices[index].pt_mp, C.pi); //source
//			
//			
//			print_comp_mp_matlab(&S->projection_values_mp[ii]->coord[0],"computed");
//			print_comp_mp_matlab(C.vertices[index].projVal_mp,"stored");									 
			
			//mid point
			index = C.edges[ii].midpt;
			vec_cp_mp(S->vertices_mp[ii][1],C.vertices[index].pt_mp);
			set_mp(&(S->projection_values_mp[ii]->coord[1]),C.vertices[index].projVal_mp);
//			projection_value_homogeneous_input(&(S->projection_values_mp[ii]->coord[1]), // destination
//																					 C.vertices[index].pt_mp, C.pi); //source
			
			
//			print_comp_mp_matlab(&S->projection_values_mp[ii]->coord[1],"computed");
//			print_comp_mp_matlab(C.vertices[index].projVal_mp,"stored");
			
			//right point
			index = C.edges[ii].right;
			vec_cp_mp(S->vertices_mp[ii][2],C.vertices[index].pt_mp);
			set_mp(&(S->projection_values_mp[ii]->coord[2]),C.vertices[index].projVal_mp);
//			projection_value_homogeneous_input(&(S->projection_values_mp[ii]->coord[2]), // destination
//																					 C.vertices[index].pt_mp, C.pi); //source
//			
//			print_comp_mp_matlab(&S->projection_values_mp[ii]->coord[2],"computed");
//			print_comp_mp_matlab(C.vertices[index].projVal_mp,"stored");
//			mypause();
		}
		
	}
	return 0;
}


int  setup_curveDecomp(int argC, char *args[], char **inputName, char **witnessSetName, char **RandMatName,  char **samplingNamenew, curveDecomp_d *C,int *num_vars)
/***************************************************************\
 * USAGE:    setup curveDecomp structure and inputname
 * ARGUMENTS:                                                    *
 * RETURN VALUES:           an integer                           *
 * NOTES:                                                        *
 \***************************************************************/
{
	FILE *IN;
	int strLength, MPType;
	char *directoryName=NULL,tmp_file[1000];
	

	
	//setup the name of directory
	IN = safe_fopen_read("Dir_Name");
	fscanf(IN, "%d\n", &strLength);
	directoryName = (char *)bmalloc(strLength * sizeof(char));
	fgets(directoryName, strLength, IN);
	fscanf(IN, "%d\n", &MPType);
	fclose(IN);
	
	//read in the gross information from the summary file.
	sprintf(tmp_file,  "%s/C.curve", directoryName);
	setup_curve(C,tmp_file,MPType);
	
	//setup E structure from E.edge
	sprintf(tmp_file,  "%s/E.edge", directoryName);
	C->num_edges = setup_edges(&(C->edges),tmp_file,num_vars,inputName,directoryName,MPType);
	
//	//setup V0 structure from V0.vert
//	sprintf(tmp_file,  "%s/V0.vert", directoryName);
//	C->num_V0 = setup_vertices(&(C->V0),tmp_file,*num_vars,MPType);
	
	//setup V1 structure from V1.vert
	sprintf(tmp_file,  "%s/V.vert", directoryName);
	C->num_vertices = setup_vertices(&(C->vertices),tmp_file,*num_vars,MPType);
	
	
	
	strLength = 1 + snprintf(NULL, 0, "%s/witness_set", directoryName);
	*witnessSetName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*witnessSetName, "%s/witness_set", directoryName);
	
	strLength = 1 + snprintf(NULL, 0, "%s/Rand_Matrix", directoryName);
	*RandMatName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*RandMatName, "%s/Rand_Matrix", directoryName);
	
	strLength = 1 + snprintf(NULL, 0, "%s/samp.dat", directoryName);
	*samplingNamenew = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*samplingNamenew, "%s/samp.dat", directoryName);
	if(! C->num_edges)
	{
		printf("sampler will not generate sampling data since there are no edges\n");
		return 1;
	}
	
	return 0;
}

int setup_edges(edge **edges,char *INfile,int *num_vars, char **inputName, char *directoryName, int MPType)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_edges,ii,strLength;
	char *input_deflated_Name=NULL;
	
	fscanf(IN, "%d\n", num_vars);
	fscanf(IN, "%d\n", &num_edges);
	fscanf(IN, "%d\n", &strLength);
	
	input_deflated_Name = (char *)bmalloc((strLength+1) * sizeof(char));
	fgets(input_deflated_Name, (strLength+1), IN);
	
	strLength = 1 + snprintf(NULL, 0, "%s/%s", directoryName,input_deflated_Name);
	*inputName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*inputName,  "%s/%s", directoryName,input_deflated_Name);
	
	*edges=(edge*) bmalloc(num_edges*sizeof(edge));
	for(ii=0;ii<num_edges;ii++)
	{
		fscanf(IN,"%d %d %d",&((*edges)[ii].left),&((*edges)[ii].midpt),&((*edges)[ii].right));
		scanRestOfLine(IN);
			
//		printf("edge %d of %d: = ",ii, num_edges);
//		printf("(%d, %d, %d) \n", ((*edges)[ii].left), ((*edges)[ii].midpt), ((*edges)[ii].right));
		
	}
	
	fclose(IN);
	return num_edges;
}


int setup_vertices(vertex **vertices,char *INfile,int num_vars, int MPType)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices,ii,jj;
	fscanf(IN, "%d\n\n", &num_vertices);
	*vertices=(vertex* )bmalloc(num_vertices*sizeof(vertex));
	for(ii=0;ii<num_vertices;ii++)
	{
		if(MPType==0)
		{
			init_point_d((*vertices)[ii].pt,num_vars);
			(*vertices)[ii].pt->size=num_vars;
			for(jj=0;jj<num_vars;jj++)
			{
				fscanf(IN, "%lf %lf", &((*vertices)[ii].pt->coord[jj].r), &((*vertices)[ii].pt->coord[jj].i));
			}
			fscanf(IN, "%lf %lf", &((*vertices)[ii].projVal->r), &((*vertices)[ii].projVal->i));
		}
		else
		{
			init_point_mp((*vertices)[ii].pt_mp,num_vars); (*vertices)[ii].pt_mp->size=num_vars;
			
			for(jj=0;jj<num_vars;jj++)
			{
				mpf_inp_str((*vertices)[ii].pt_mp->coord[jj].r, IN, 10);
				mpf_inp_str((*vertices)[ii].pt_mp->coord[jj].i, IN, 10);
			}
			init_mp((*vertices)[ii].projVal_mp);
			mpf_inp_str((*vertices)[ii].projVal_mp->r, IN, 10);
			mpf_inp_str((*vertices)[ii].projVal_mp->i, IN, 10);
		}
		
		//regardless of mptype, read the type of the vertex
		fscanf(IN,"%d\n",&(*vertices)[ii].type);
	}
	fclose(IN);
	return num_vertices;
}






int setup_curve(curveDecomp_d *C,char *INfile, int MPType)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices,ii;
	
	fscanf(IN, "%d %d %d", &C->num_variables, &C->num_vertices, &C->num_edges);  scanRestOfLine(IN);
	
	fscanf(IN, "%d %d %d %d", &C->num_V0, &C->num_V1, &C->num_midpts, &C->num_new);  scanRestOfLine(IN);
	
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
	
	
	if(MPType==0)
	{
		init_point_d(C->pi_d,C->num_variables); C->pi_d->size=C->num_variables;
		for(ii=0;ii<C->num_variables;ii++)
		{
			fscanf(IN, "%lf %lf", &(C->pi_d->coord[ii].r), &(C->pi_d->coord[ii].i));
		}
	}
	else
	{
		init_point_mp(C->pi,C->num_variables); C->pi->size=C->num_variables;
		for(ii=0;ii<C->num_variables;ii++)
		{
			mpf_inp_str(C->pi->coord[ii].r, IN, 10);
			mpf_inp_str(C->pi->coord[ii].i, IN, 10);
		}
		
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

	
	// in other words, result += (x+y)/2 . pi

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








