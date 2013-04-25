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

	
	sampler_configuration sampler_options;  init_sampler_config(&sampler_options);
	sampler_parse_options(argC, args, &sampler_options);
	
	////
	//  begin the actual program
	////
	
	if (setup_curveDecomp(argC, args, &inputName, &witnessSetName,&RandMatName,&samplingNamenew,&C,&num_vars))
		return 1;
	
	
	
	//< prints the welcome message, also gets the inputName, witnessSetName and C
	srand(time(NULL));
	// essentials for using the bertini parser

	
	unsigned int currentSeed;
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
	
	
	printf("parsing witness set\n");
	witnessSetParse(&W, witnessSetName,num_vars);
	W.num_var_gps = solve_options.PPD.num_var_gp;
	W.MPType = MPType;
	get_variable_names(&W);
//	W.incidence_number = get_component_number(W,num_vars,inputName, sampler_options.stifle_text);
	
	
	printf("getting randomizer matrix\n");
	read_rand_matrix(RandMatName, n_minusone_randomizer_matrix);
	
	printf("loading sampling data\n");
	Load_sampling_data(&S_old,C,num_vars,MPType);
	
	
	
	/////////////
	//
	//  Generate new sampling data
	//
	//////////////////
	
	printf("generate_new_sampling_pts\n");
	generate_new_sampling_pts(&S_new,n_minusone_randomizer_matrix,S_old, C,W, num_vars,currentSeed,MPType, &solve_options);
	
	//output
	output_sampling_data(S_new,samplingNamenew,num_vars,MPType);
	
	
	
	
	
	// clear memory
	free(inputName);
	free(witnessSetName);
	
	clear_mat_mp(n_minusone_randomizer_matrix);
	clear_witness_set(W);
	
	clear_sample_d(&S_old, MPType);
	clear_sample_d(&S_new, MPType);
	clear_curveDecomp_d(&C, MPType);
	clearMP();
	return 0;
}







void generate_new_sampling_pts(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,sample_d S_old, curveDecomp_d C,witness_set W, int num_vars,unsigned int currentSeed, int MPType , solver_configuration *solve_options)
{
	if(MPType==0)
		generate_new_sampling_pts_d(S_new, n_minusone_randomizer_matrix, S_old, C, W, num_vars, currentSeed, MPType, solve_options);
	else
		generate_new_sampling_pts_mp(S_new, n_minusone_randomizer_matrix, S_old, C, W, num_vars, currentSeed, MPType, solve_options);
	return;
}


void generate_new_sampling_pts_d(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,
																 sample_d S_old, curveDecomp_d C,witness_set W, int num_vars,
																 unsigned int currentSeed, int MPType, solver_configuration *solve_options)
{
	witness_set		Wnew;
	vec_mp         new_linears;
	vec_d          L,startpt,mid_pt;
	int            ii,jj,k;
	comp_d         pi_end,temp,temp1;
	vec_d          *previous_points, *new_points, previous_projection_values, new_projection_values;
	int            num_samp_old, num_samp_new;
	int            *refine_old, *refine_new;
	double         max_norm, TOL=1e-1;
	

	
	init_vec_mp(new_linears,num_vars); new_linears->size = num_vars;
	
	init_vec_d(L,num_vars);  L->size = num_vars; // what is L?
	
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
		num_samp_old = S_old.num_pts[ii]; // copy number of points
		
		
		init_vec_d(previous_projection_values,S_old.num_pts[ii]);
		previous_projection_values->size = S_old.num_pts[ii];
		
		vec_cp_d(previous_projection_values, S_old.projection_values[ii]);
		refine_old = S_old.refine[ii];
		
		
		while(1) {
			new_points = (vec_d *)bmalloc(2*num_samp_old * sizeof(vec_d)); // why the 2?
			init_vec_d(new_projection_values,2*num_samp_old); // what is projection_values?
			//new_projection_values->size=2*num_samp_old;
			refine_new = (int * )bmalloc(2*num_samp_old * sizeof(int));
			//copy left point
			set_d(&(new_projection_values->coord[0]),&(previous_projection_values->coord[0]));
			num_samp_new = 1;
			init_vec_d(new_points[0],num_vars);
			vec_cp_d(new_points[0],previous_points[0]);
			for(jj=0;jj<num_samp_old-1;jj++)
			{
				if(refine_old[jj])
				{
//TODO: resolve this discrepancy
					vec_cp_d(L,C.pi_d);
					point_d_to_mp(new_linears,L);
					if(jj==0)
					{
						vec_cp_d(startpt,previous_points[jj+1]);
						sub_d(&(L->coord[0]),&(L->coord[0]),&(previous_projection_values->coord[jj+1]));
					}
					else
					{
						vec_cp_d(startpt,previous_points[jj]);
						sub_d(&(L->coord[0]),&(L->coord[0]),&(previous_projection_values->coord[jj]));
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
					d_to_mp(&(new_linears->coord[0]),temp);
					
					// copy in the data to the source witness set?
					set_witness_set_d(&W, L,startpt,num_vars);
					//printVec_d(stdout,0,mid_pt);
					//print_d(stdout,0,pi_end);
					//printf("L\n");
					//printVec_d(stdout,0,L);
					//printf("newlinears\n");
					//printVec_mp(stdout,0,new_linears);
					
					lin_to_lin_solver_main(MPType,
																 W,         // witness_set
																 n_minusone_randomizer_matrix,
																 &new_linears, //  the set of linears we will solve at.
																 1, // the number of new linears.
																 &Wnew,// the new data is put here!
																 solve_options);
					//printVec_d(stdout,0,mid_pt);
					
					set_zero_d(&(new_projection_values->coord[num_samp_new]));
					init_vec_d(new_points[num_samp_new],num_vars);
					vec_cp_d(new_points[num_samp_new],Wnew.pts_d[0]);
					
					//printVec_d(stdout,0,Wnew.pts_d[0]);
					//printVec_d(stdout,0,mid_pt);
					max_norm=0.0;
					for(k=1;k<num_vars;k++)
					{
						div_d(temp1,&(new_points[num_samp_new]->coord[k]),
									&(new_points[num_samp_new]->coord[0]));
						div_d(temp,&(mid_pt->coord[k]),&(mid_pt->coord[0]));
						sub_d(temp,temp1,temp);
						if(max_norm<d_abs_d(temp))
							max_norm = d_abs_d(temp);
					}
					//printf("pi_end");
					//print_d(stdout,0,pi_end);
					set_d(&(new_projection_values->coord[num_samp_new]),pi_end);
					
					if(max_norm<TOL)
					{
						refine_new[num_samp_new-1] = 0;
						refine_new[num_samp_new++] = 0;
					}
					else
					{
						refine_new[num_samp_new-1] = 1;
						refine_new[num_samp_new++] = 1;
					}
					set_d(&(new_projection_values->coord[num_samp_new]),&(previous_projection_values->coord[jj+1]));
					init_vec_d(new_points[num_samp_new],num_vars);
					vec_cp_d(new_points[num_samp_new],previous_points[jj+1]);
					num_samp_new++;
					
				}
				else
				{
					refine_new[num_samp_new-1] = 0;
					set_d(&(new_projection_values->coord[num_samp_new]),&(previous_projection_values->coord[jj+1]));
					init_vec_d(new_points[num_samp_new],num_vars);
					vec_cp_d(new_points[num_samp_new],previous_points[jj+1]);
					num_samp_new++;
				}
			}
			new_points = (vec_d *)brealloc(new_points,num_samp_new*sizeof(vec_d));
			new_projection_values->size=num_samp_new;
			refine_new = (int * )brealloc(refine_new,num_samp_new * sizeof(int));
			
			previous_points = new_points;
			vec_cp_d(previous_projection_values,new_projection_values);
			clear_vec_d(new_projection_values);
			refine_old = refine_new;
			if(num_samp_new == num_samp_old)
			{
				S_new->vertices[ii] = previous_points;
				S_new->num_pts[ii] = num_samp_old;
				init_vec_d(S_new->projection_values[ii],num_samp_old);
				vec_cp_d(S_new->projection_values[ii],previous_projection_values);
				clear_vec_d(previous_projection_values);
				S_new->refine[ii] = refine_old;
				break;
			}
			else
				num_samp_old=num_samp_new;
		}
	}
	
	clear_witness_set(Wnew);
	clear_vec_d(L);
	clear_vec_d(new_linears);
}



void generate_new_sampling_pts_mp(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,
																	sample_d S_old, curveDecomp_d C,witness_set W,
																	int num_vars,unsigned int currentSeed, int MPType, solver_configuration *solve_options)
{
	witness_set  Wnew;
	vec_mp         new_linears;
	vec_mp          L,startpt,mid_pt;
	int            ii,jj,k;
	comp_mp         pi_end,temp,temp1;
	vec_mp          *previous_points, *new_points, previous_projection_values, new_projection_values;
	int            num_samp_old, num_samp_new;
	int            *refine_old, *refine_new;
	double         max_norm, TOL=1e-1;
	
	
	
	init_vec_mp(new_linears,num_vars); new_linears->size = num_vars;
	
	init_vec_mp(L,num_vars); // what is the purpose of this?
	L->size = num_vars; // what is L?
	
	init_vec_mp(startpt,num_vars); startpt->size = num_vars;
	
	init_vec_mp(mid_pt,num_vars); mid_pt->size = num_vars;
	
	init_mp(temp); init_mp(pi_end); init_mp(temp1);
	S_new->num_edges = S_old.num_edges;
	S_new->vertices_mp = (vec_mp **)bmalloc(S_new->num_edges * sizeof(vec_mp*));
	S_new->projection_values_mp = (vec_mp *)bmalloc(S_new->num_edges * sizeof(vec_mp));
	S_new->refine = (int **)bmalloc(S_new->num_edges * sizeof(int*));
	S_new->num_pts = (int *) bmalloc(S_new->num_edges * sizeof(int));
	
	
	
	for(ii=0;ii<S_old.num_edges;ii++) // for each of the edges
	{
		previous_points = S_old.vertices_mp[ii]; // does this copy a pointer?
		num_samp_old = S_old.num_pts[ii]; // grab the number of points from the array of integers
		
		init_vec_mp(previous_projection_values,S_old.num_pts[ii]);
		previous_projection_values->size = S_old.num_pts[ii]; //
		
		
		
		vec_cp_mp(previous_projection_values, S_old.projection_values_mp[ii]); // set the edge_proj_old to be the old set of projections???
		
		refine_old = S_old.refine[ii]; // change the pointer of refine_old to be S_old.refine[ii]
		
		
		while(1) // breaking condition is all samples being less than TOL away from each other (in the infty norm sense).
		{
			new_points = (vec_mp *)bmalloc(2*num_samp_old * sizeof(vec_mp));
			
			init_vec_mp(new_projection_values,2*num_samp_old); new_projection_values->size=2*num_samp_old; // hold the values of the projections in new_projection_values???
			
			refine_new = (int * )bmalloc(2*num_samp_old * sizeof(int)); // what is this?
			
			//copy left point and projection value
			set_mp(&(new_projection_values->coord[0]),&(previous_projection_values->coord[0]));//leftmost point always the same?
			init_vec_mp(new_points[0],num_vars); new_points[0]->size = num_vars; // initialize the first of the samples???
			vec_cp_mp(new_points[0],previous_points[0]);
			
			
			num_samp_new = 1; // initialize a counter???
			print_point_to_screen_matlab_mp(previous_projection_values,"old_projections");
			
	
			for(jj=0;jj<num_samp_old-1;jj++) // for each sample in the previous set
			{
				if(refine_old[jj]) // 
				{
					vec_cp_mp(L,C.pi); // grab the projection, copy it into L
					vec_cp_mp(new_linears,C.pi); // copy the projection into new_linears
					
					if(jj==0)// if on the first sample???
					{
						vec_cp_mp(startpt,previous_points[jj+1]);
						set_mp(&(L->coord[0]),&(previous_projection_values->coord[jj+1]));
						neg_mp(&(L->coord[0]),&(L->coord[0]));
//						sub_mp(&(L->coord[0]),&(L->coord[0]),&(previous_projection_values->coord[jj+1]));// L[0]= L[0]-projection???
					}
					else
					{
						vec_cp_mp(startpt,previous_points[jj]);
						set_mp(&(L->coord[0]),&(previous_projection_values->coord[jj+1]));
						neg_mp(&(L->coord[0]),&(L->coord[0]));
					}
					
					 //HERE!!!
					//DAB -- I have no idea what is going on here...???
					add_mp(&(mid_pt->coord[0]),&(previous_points[jj]->coord[0]),
								 &(previous_points[jj+1]->coord[0]));
					mul_mp(temp,&(previous_points[jj]->coord[0]),
								 &(previous_points[jj+1]->coord[0]));
					div_mp(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mpf_set_d(temp->r, 2.0); mpf_set_d(temp->i, 0.0);
					mul_mp(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mul_mp(pi_end,&(mid_pt->coord[0]),&(C.pi->coord[0]));
					for(k=1;k<num_vars;k++)
					{
						div_mp(&(mid_pt->coord[k]),&(previous_points[jj]->coord[k]),
									 &(previous_points[jj]->coord[0]));
						div_mp(temp,&(previous_points[jj+1]->coord[k]),
									 &(previous_points[jj+1]->coord[0]));
						add_mp(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						
						mpf_set_d(temp->r, 0.5); mpf_set_d(temp->i, 0.0);
						mul_mp(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						mul_mp(&(mid_pt->coord[k]),&(mid_pt->coord[0]),&(mid_pt->coord[k]));
						mul_mp(temp,&(mid_pt->coord[k]),&(C.pi->coord[k]));
						add_mp(pi_end,pi_end,temp);
					}
					div_mp(pi_end,pi_end,&(mid_pt->coord[0]));
					sub_mp(temp,&(C.pi->coord[0]),pi_end);
					
					set_mp(&(new_linears->coord[0]),temp); // set the value of the projection we want to move to.
					set_witness_set_mp(&W, L,startpt,num_vars); // set the witness point and linear in the input for the lintolin solver.
					
					
					init_witness_set_d(&Wnew);
					
//					print_point_to_screen_matlab_mp(W.L_mp[0],"initial_projection");
//					print_point_to_screen_matlab_mp(new_linears,"destination");
//					mypause(); 
					
					
					lin_to_lin_solver_main(MPType,
																 W,         // witness_set
																 n_minusone_randomizer_matrix,
																 &new_linears, //  the set of linears we will solve at.
																 1, // the number of new linears.
																 &Wnew, // the new data is put here!
																 solve_options);
					set_zero_mp(&(new_projection_values->coord[num_samp_new]));
					init_vec_mp(new_points[num_samp_new],num_vars);
					vec_cp_mp(new_points[num_samp_new],Wnew.pts_mp[0]);
					clear_witness_set(Wnew);
					
					
					max_norm=0.0; // initialize
					for(k=1;k<num_vars;k++)
					{
						div_mp(temp1,&(new_points[num_samp_new]->coord[k]),
									 &(new_points[num_samp_new]->coord[0]));
						div_mp(temp,&(mid_pt->coord[k]),&(mid_pt->coord[0]));
						sub_mp(temp,temp1,temp);
						if(max_norm<d_abs_mp(temp))
							max_norm = d_abs_mp(temp);
					}
					set_mp(&(new_projection_values->coord[num_samp_new]),pi_end);
					//printf("pi_end=");print_mp(stdout,0,pi_end);
					if(max_norm<TOL)
					{
						refine_new[num_samp_new-1] = 0;
						refine_new[num_samp_new++] = 0;
					}
					else
					{
						refine_new[num_samp_new-1] = 1;
						refine_new[num_samp_new++] = 1;
					}
					set_mp(&(new_projection_values->coord[num_samp_new]),&(previous_projection_values->coord[jj+1]));
					//printf("new_projection_values=");print_mp(stdout,0,&(new_projection_values->coord[num_samp_new]));
					init_vec_mp(new_points[num_samp_new],num_vars);
					vec_cp_mp(new_points[num_samp_new],previous_points[jj+1]);
					num_samp_new++;
					
				}
				else
				{
					refine_new[num_samp_new-1] = 0;
					set_mp(&(new_projection_values->coord[num_samp_new]),&(previous_projection_values->coord[jj+1]));
					init_vec_mp(new_points[num_samp_new],num_vars);
					vec_cp_mp(new_points[num_samp_new],previous_points[jj+1]);
					num_samp_new++;
				}
			}
			new_points = (vec_mp *)brealloc(new_points,num_samp_new*sizeof(vec_mp));
			//printVec_mp(stdout,0,new_projection_values);printf("==%d==%d\n",num_samp_new,num_samp_old);
			//change_size_vec_mp(new_projection_values,num_samp_new); // what is pV?
			new_projection_values->size=num_samp_new;
			refine_new = (int * )brealloc(refine_new,num_samp_new * sizeof(int));
			
			//printVec_mp(stdout,0,new_projection_values);
			previous_points = new_points;
			vec_cp_mp(previous_projection_values,new_projection_values);
			//printVec_mp(stdout,0,previous_projection_values);
			clear_vec_mp(new_projection_values);
			refine_old = refine_new;
			if(num_samp_new == num_samp_old) // if had no new samples???
			{
				S_new->vertices_mp[ii] = previous_points;
				S_new->num_pts[ii] = num_samp_old;
				init_vec_mp(S_new->projection_values_mp[ii],num_samp_old);
				vec_cp_mp(S_new->projection_values_mp[ii],previous_projection_values);
				clear_vec_mp(previous_projection_values);
				S_new->refine[ii] = refine_old;
				break; // BREAKS THE WHILE LOOP
			}
			else
				num_samp_old=num_samp_new; // update the number of samples
		}
	}
	clear_mp(temp); clear_mp(pi_end); clear_mp(temp1);
	
	
	clear_vec_mp(L);
	clear_vec_mp(new_linears);
}


void read_rand_matrix(char *INfile, mat_mp n_minusone_randomizer_matrix)
{
	int ii,jj,rows,cols,cur_precision;
	FILE *IN= safe_fopen_read(INfile);
	fscanf(IN,"%d\n",&rows);
	fscanf(IN,"%d\n",&cols);
	fscanf(IN,"%d\n",&cur_precision);
	init_mat_mp2(n_minusone_randomizer_matrix,0,0,cur_precision);
	increase_size_mat_mp(n_minusone_randomizer_matrix, rows, cols);
	n_minusone_randomizer_matrix->rows=rows;
	n_minusone_randomizer_matrix->cols=cols;
	for(ii=0;ii<n_minusone_randomizer_matrix->rows;ii++)
		for(jj=0;jj<n_minusone_randomizer_matrix->cols;jj++)
		{
			fscanf(IN,"<");
			mpf_inp_str(n_minusone_randomizer_matrix->entry[ii][jj].r, IN, 10);
			mpf_inp_str(n_minusone_randomizer_matrix->entry[ii][jj].i, IN, 10);
			fscanf(IN,">");
		}
	
	fclose(IN);
}

void set_witness_set_d(witness_set *W, vec_d L,vec_d pts,int num_vars)
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
	
	vec_cp_d(W->L[0],L);
	vec_d_to_mp(W->L_mp[0],L);
	
}

void set_witness_set_mp(witness_set *W, vec_mp L,vec_mp pts,int num_vars)
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
	
	vec_cp_mp(W->L_mp[0],L);
	vec_mp_to_d(W->L[0],L);
	
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
				S->refine[ii][jj]=1;
			}
			
			init_vec_mp(S->projection_values_mp[ii],S->num_pts[ii]); S->projection_values_mp[ii]->size=S->num_pts[ii];
			
			
			vec_mp pi_nohom;  init_vec_mp(pi_nohom,num_vars-1); pi_nohom->size = num_vars-1;
			int mm;
			for (mm=0; mm<num_vars-1; mm++) {
				set_mp(&pi_nohom->coord[mm],& C.pi->coord[mm]);
			}
			
			vec_mp dehom; init_vec_mp(dehom,num_vars-1); dehom->size = num_vars-1;
			
			
			//left point
			index = C.edges[ii].left;
			vec_cp_mp(S->vertices_mp[ii][0],C.vertices[index].pt_mp);
			dehomogenize_mp(&dehom,C.vertices[index].pt_mp);
			dot_product_mp(&(S->projection_values_mp[ii]->coord[0]), dehom,pi_nohom);

			
			//mid point
			vec_cp_mp(S->vertices_mp[ii][1],C.vertices[C.edges[ii].midpt].pt_mp);
			dehomogenize_mp(&dehom,S->vertices_mp[ii][1]);
			dot_product_mp(&(S->projection_values_mp[ii]->coord[1]), dehom,pi_nohom);
			
			
			//right point
			index = C.edges[ii].right;
			vec_cp_mp(S->vertices_mp[ii][2],C.vertices[index].pt_mp);
			dehomogenize_mp(&dehom,C.vertices[index].pt_mp);
			dot_product_mp(&(S->projection_values_mp[ii]->coord[1]), dehom,pi_nohom);

			clear_vec_mp(pi_nohom);
			clear_vec_mp(dehom);
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



