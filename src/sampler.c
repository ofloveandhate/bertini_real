#include "sampler.h"


int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
	int num_vars=0;  //1=self-conjugate; 0=not
	char *inputName = NULL, *RandMatName = NULL, *witnessSetName = NULL, *samplingNamenew = NULL;
	curveDecomp_d C;  //new data type; stores vertices, edges, etc.
	witness_set_d W;
	sample_d   S_old,S_new;
	//mat_d n_minusone_randomizer_matrix;
	
	////
	//  begin the actual program
	////
	
	if (setup_curveDecomp(argC, args, &inputName, &witnessSetName,&RandMatName,&samplingNamenew,&C,&num_vars))
		return 1;
	
	
	
	//< prints the welcome message, also gets the inputName, witnessSetName and C
	srand(time(NULL));
	// essentials for using the bertini parser
	prog_t SLP;
	
	unsigned int currentSeed;
	int trackType, genType = 0, MPType,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0;
	int my_id=0, num_processes=1, headnode = 0; // headnode is always 0

	int num_var_gps = 0, userHom = 0;
	mat_mp n_minusone_randomizer_matrix;
	
	//end parser-bertini essentials
	
	parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
	preproc_data PPD;
	setupPreProcData("preproc_data", &PPD);
	num_var_gps = PPD.num_var_gp;

	tracker_config_t T;
	get_tracker_config(&T,MPType);
	
	initMP(T.Precision);
	
	
	num_vars = setupProg(&SLP, T.Precision, MPType); // num_vars includes the number of homogeneous coordinates.
	// the number of homogeneous coordinates is the num_var_gps.
	
	
	
	printf("parsing witness set\n");
	witnessSetParse(&W, witnessSetName,num_vars);
	read_rand_matrix(RandMatName, n_minusone_randomizer_matrix);
	printf("Loading sampling data\n");
	Load_sampling_data(&S_old,C,num_vars,MPType);
	
	//Generate new sampling data
	printf("generate_new_sampling_pts\n");
	generate_new_sampling_pts(&S_new,n_minusone_randomizer_matrix,S_old, C,W, num_vars,currentSeed,MPType);
	
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







void generate_new_sampling_pts(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,sample_d S_old, curveDecomp_d C,witness_set_d W, int num_vars,unsigned int currentSeed, int MPType )
{
	if(MPType==0)
		generate_new_sampling_pts_d(S_new, n_minusone_randomizer_matrix, S_old, C, W, num_vars, currentSeed, MPType);
	else
		generate_new_sampling_pts_mp(S_new, n_minusone_randomizer_matrix, S_old, C, W, num_vars, currentSeed, MPType);
	return;
}


void generate_new_sampling_pts_d(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,sample_d S_old, curveDecomp_d C,witness_set_d W, int num_vars,unsigned int currentSeed, int MPType)
{
	witness_set_d Wnew;
	vec_mp         new_linears;
	vec_d          L,startpt,mid_pt;
	int            ii,jj,k;
	comp_d         pi_end,temp,temp1;
	vec_d          *Edge_samp_old, *Edge_samp_new, Edge_proj_old, Edge_proj_new;
	int            num_samp_old, num_samp_new;
	int            *refine_old, *refine_new;
	double         max_norm, TOL=1e-1;
	
	init_witness_set_d(&Wnew);
	cp_patches(&Wnew,W);
	
	init_vec_mp(new_linears,num_vars);
	new_linears->size = num_vars;
	
	init_vec_d(L,num_vars); // what is the purpose of this?
	L->size = num_vars; // what is L?
	
	init_vec_d(startpt,num_vars);
	startpt->size = num_vars;
	
	init_vec_d(mid_pt,num_vars);
	mid_pt->size = num_vars;
	
	
	S_new->num_edges = S_old.num_edges;
	S_new->vertices = (vec_d **)bmalloc(S_new->num_edges * sizeof(vec_d*));
	S_new->proj_vertices = (vec_d *)bmalloc(S_new->num_edges * sizeof(vec_d));
	S_new->refine = (int **)bmalloc(S_new->num_edges * sizeof(int*));
	S_new->num_pts = (int *) bmalloc(S_new->num_edges * sizeof(int));
	
	for(ii=0;ii<S_old.num_edges;ii++)
	{
		Edge_samp_old = S_old.vertices[ii];
		num_samp_old = S_old.num_pts[ii];
		init_vec_d(Edge_proj_old,S_old.num_pts[ii]);
		vec_cp_d(Edge_proj_old, S_old.proj_vertices[ii]);
		refine_old = S_old.refine[ii];
		while(1)
		{
			Edge_samp_new = (vec_d *)bmalloc(2*num_samp_old * sizeof(vec_d));
			init_vec_d(Edge_proj_new,2*num_samp_old); // what is proj_vertices?
			//Edge_proj_new->size=2*num_samp_old;
			refine_new = (int * )bmalloc(2*num_samp_old * sizeof(int));
			//copy left point
			set_d(&(Edge_proj_new->coord[0]),&(Edge_proj_old->coord[0]));
			num_samp_new = 1;
			init_vec_d(Edge_samp_new[0],num_vars);
			vec_cp_d(Edge_samp_new[0],Edge_samp_old[0]);
			for(jj=0;jj<num_samp_old-1;jj++)
			{
				if(refine_old[jj])
				{
					vec_cp_d(L,C.edges[ii].pi);
					point_d_to_mp(new_linears,L);
					if(jj==0)
					{
						vec_cp_d(startpt,Edge_samp_old[jj+1]);
						sub_d(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[jj+1]));
					}
					else
					{
						vec_cp_d(startpt,Edge_samp_old[jj]);
						sub_d(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[jj]));
					}
					//tracking end point pi_end = \pi (si+s_{ii+1})/2
					//w=2w1w2/(w1+w2)
					add_d(&(mid_pt->coord[0]),&(Edge_samp_old[jj]->coord[0]),
								&(Edge_samp_old[jj+1]->coord[0]));
					mul_d(temp,&(Edge_samp_old[jj]->coord[0]),
								&(Edge_samp_old[jj+1]->coord[0]));
					div_d(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					temp->r=2.0; temp->i=0.0;
					mul_d(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mul_d(pi_end,&(mid_pt->coord[0]),&(C.edges[ii].pi->coord[0]));
					
					for(k=1;k<num_vars;k++)
					{
						div_d(&(mid_pt->coord[k]),&(Edge_samp_old[jj]->coord[k]),
									&(Edge_samp_old[jj]->coord[0]));
						div_d(temp,&(Edge_samp_old[jj+1]->coord[k]),
									&(Edge_samp_old[jj+1]->coord[0]));
						add_d(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						temp->r=0.5; temp->i=0.0;
						mul_d(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						mul_d(&(mid_pt->coord[k]),&(mid_pt->coord[0]),&(mid_pt->coord[k]));
						mul_d(temp,&(mid_pt->coord[k]),&(C.edges[ii].pi->coord[k]));
						add_d(pi_end,pi_end,temp);
					}
					div_d(pi_end,pi_end,&(mid_pt->coord[0]));
					sub_d(temp,&(C.edges[ii].pi->coord[0]),pi_end);
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
																 &Wnew); // the new data is put here!
					//printVec_d(stdout,0,mid_pt);
					
					set_zero_d(&(Edge_proj_new->coord[num_samp_new]));
					init_vec_d(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_d(Edge_samp_new[num_samp_new],Wnew.W.pts[0]);
					
					//printVec_d(stdout,0,Wnew.W.pts[0]);
					//printVec_d(stdout,0,mid_pt);
					max_norm=0.0;
					for(k=1;k<num_vars;k++)
					{
						div_d(temp1,&(Edge_samp_new[num_samp_new]->coord[k]),
									&(Edge_samp_new[num_samp_new]->coord[0]));
						div_d(temp,&(mid_pt->coord[k]),&(mid_pt->coord[0]));
						sub_d(temp,temp1,temp);
						if(max_norm<d_abs_d(temp))
							max_norm = d_abs_d(temp);
					}
					//printf("pi_end");
					//print_d(stdout,0,pi_end);
					set_d(&(Edge_proj_new->coord[num_samp_new]),pi_end);
					
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
					set_d(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[jj+1]));
					init_vec_d(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_d(Edge_samp_new[num_samp_new],Edge_samp_old[jj+1]);
					num_samp_new++;
					
				}
				else
				{
					refine_new[num_samp_new-1] = 0;
					set_d(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[jj+1]));
					init_vec_d(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_d(Edge_samp_new[num_samp_new],Edge_samp_old[jj+1]);
					num_samp_new++;
				}
			}
			Edge_samp_new = (vec_d *)brealloc(Edge_samp_new,num_samp_new*sizeof(vec_d));
			Edge_proj_new->size=num_samp_new;
			refine_new = (int * )brealloc(refine_new,num_samp_new * sizeof(int));
			
			Edge_samp_old = Edge_samp_new;
			vec_cp_d(Edge_proj_old,Edge_proj_new);
			clear_vec_d(Edge_proj_new);
			refine_old = refine_new;
			if(num_samp_new == num_samp_old)
			{
				S_new->vertices[ii] = Edge_samp_old;
				S_new->num_pts[ii] = num_samp_old;
				init_vec_d(S_new->proj_vertices[ii],num_samp_old);
				vec_cp_d(S_new->proj_vertices[ii],Edge_proj_old);
				clear_vec_d(Edge_proj_old);
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
																	sample_d S_old, curveDecomp_d C,witness_set_d W,
																	int num_vars,unsigned int currentSeed, int MPType)
{
	witness_set_d  Wnew;
	vec_mp         new_linears;
	vec_mp          L,startpt,mid_pt;
	int            ii,jj,k;
	comp_mp         pi_end,temp,temp1;
	vec_mp          *Edge_samp_old, *Edge_samp_new, Edge_proj_old, Edge_proj_new;
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
	S_new->proj_vertices_mp = (vec_mp *)bmalloc(S_new->num_edges * sizeof(vec_mp));
	S_new->refine = (int **)bmalloc(S_new->num_edges * sizeof(int*));
	S_new->num_pts = (int *) bmalloc(S_new->num_edges * sizeof(int));
	for(ii=0;ii<S_old.num_edges;ii++) // for each of the edges
	{
		Edge_samp_old = S_old.vertices_mp[ii];
		num_samp_old = S_old.num_pts[ii];
		init_vec_mp(Edge_proj_old,S_old.num_pts[ii]); Edge_proj_old->size = S_old.num_pts[ii];
		vec_cp_mp(Edge_proj_old, S_old.proj_vertices_mp[ii]); // set the edge_proj_old to be the old set of projections???
		refine_old = S_old.refine[ii]; // change the pointer of refine_old to be S_old.refine[ii]
		
		
		while(1) // breaking condition is all samples being less than TOL away from each other (in the infty norm sense).
		{
			Edge_samp_new = (vec_mp *)bmalloc(2*num_samp_old * sizeof(vec_mp));
			
			init_vec_mp(Edge_proj_new,2*num_samp_old); Edge_proj_new->size=2*num_samp_old; // hold the values of the projections in Edge_proj_new???
			
			refine_new = (int * )bmalloc(2*num_samp_old * sizeof(int)); // what is this?
			
			//copy left point
			set_mp(&(Edge_proj_new->coord[0]),&(Edge_proj_old->coord[0]));
			num_samp_new = 1; // initialize a counter???
			init_vec_mp(Edge_samp_new[0],num_vars); Edge_samp_new[0]->size = num_vars; // initialize the first of the samples???
			vec_cp_mp(Edge_samp_new[0],Edge_samp_old[0]);
			
			print_point_to_screen_matlab_mp(Edge_proj_old,"old_projections");
			
	
			for(jj=0;jj<num_samp_old-1;jj++) // for each sample in the previous set
			{
				if(refine_old[jj]) // 
				{
					vec_cp_mp(L,C.edges[ii].pi_mp); // grab the projection for the iith edge, copy it into L
					vec_cp_mp(new_linears,L); // copy the projection into new_linears
					
					if(jj==0)// if on the first sample???
					{
						vec_cp_mp(startpt,Edge_samp_old[jj+1]);
						set_mp(&(L->coord[0]),&(Edge_proj_old->coord[jj+1]));
						neg_mp(&(L->coord[0]),&(L->coord[0]));
//						sub_mp(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[jj+1]));// L[0]= L[0]-projection???
					}
					else
					{
						vec_cp_mp(startpt,Edge_samp_old[jj]);
						sub_mp(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[jj]));
					}
					
					
					//DAB -- I have no idea what is going on here...???
					add_mp(&(mid_pt->coord[0]),&(Edge_samp_old[jj]->coord[0]),
								 &(Edge_samp_old[jj+1]->coord[0]));
					mul_mp(temp,&(Edge_samp_old[jj]->coord[0]),
								 &(Edge_samp_old[jj+1]->coord[0]));
					div_mp(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mpf_set_d(temp->r, 2.0); mpf_set_d(temp->i, 0.0);
					mul_mp(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mul_mp(pi_end,&(mid_pt->coord[0]),&(C.edges[ii].pi_mp->coord[0]));
					for(k=1;k<num_vars;k++)
					{
						div_mp(&(mid_pt->coord[k]),&(Edge_samp_old[jj]->coord[k]),
									 &(Edge_samp_old[jj]->coord[0]));
						div_mp(temp,&(Edge_samp_old[jj+1]->coord[k]),
									 &(Edge_samp_old[jj+1]->coord[0]));
						add_mp(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						
						mpf_set_d(temp->r, 0.5); mpf_set_d(temp->i, 0.0);
						mul_mp(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						mul_mp(&(mid_pt->coord[k]),&(mid_pt->coord[0]),&(mid_pt->coord[k]));
						mul_mp(temp,&(mid_pt->coord[k]),&(C.edges[ii].pi_mp->coord[k]));
						add_mp(pi_end,pi_end,temp);
					}
					div_mp(pi_end,pi_end,&(mid_pt->coord[0]));
					sub_mp(temp,&(C.edges[ii].pi_mp->coord[0]),pi_end);
					
					set_mp(&(new_linears->coord[0]),temp); // set the value of the projection we want to move to.
					set_witness_set_mp(&W, L,startpt,num_vars); // set the witness point and linear in the input for the lintolin solver.
					
					
					init_witness_set_d(&Wnew);
					
					print_point_to_screen_matlab_mp(W.L_mp[0],"initial_projection");
					print_point_to_screen_matlab_mp(new_linears,"destination");
					
					mypause();  //HERE
					lin_to_lin_solver_main(MPType,
																 W,         // witness_set
																 n_minusone_randomizer_matrix,
																 &new_linears, //  the set of linears we will solve at.
																 1, // the number of new linears.
																 &Wnew); // the new data is put here!
					set_zero_mp(&(Edge_proj_new->coord[num_samp_new]));
					init_vec_mp(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_mp(Edge_samp_new[num_samp_new],Wnew.W_mp.pts[0]);
					clear_witness_set(Wnew);
					
					
					max_norm=0.0; // initialize
					for(k=1;k<num_vars;k++)
					{
						div_mp(temp1,&(Edge_samp_new[num_samp_new]->coord[k]),
									 &(Edge_samp_new[num_samp_new]->coord[0]));
						div_mp(temp,&(mid_pt->coord[k]),&(mid_pt->coord[0]));
						sub_mp(temp,temp1,temp);
						if(max_norm<d_abs_mp(temp))
							max_norm = d_abs_mp(temp);
					}
					set_mp(&(Edge_proj_new->coord[num_samp_new]),pi_end);
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
					set_mp(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[jj+1]));
					//printf("Edge_proj_new=");print_mp(stdout,0,&(Edge_proj_new->coord[num_samp_new]));
					init_vec_mp(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_mp(Edge_samp_new[num_samp_new],Edge_samp_old[jj+1]);
					num_samp_new++;
					
				}
				else
				{
					refine_new[num_samp_new-1] = 0;
					set_mp(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[jj+1]));
					init_vec_mp(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_mp(Edge_samp_new[num_samp_new],Edge_samp_old[jj+1]);
					num_samp_new++;
				}
			}
			Edge_samp_new = (vec_mp *)brealloc(Edge_samp_new,num_samp_new*sizeof(vec_mp));
			//printVec_mp(stdout,0,Edge_proj_new);printf("==%d==%d\n",num_samp_new,num_samp_old);
			//change_size_vec_mp(Edge_proj_new,num_samp_new); // what is pV?
			Edge_proj_new->size=num_samp_new;
			refine_new = (int * )brealloc(refine_new,num_samp_new * sizeof(int));
			
			//printVec_mp(stdout,0,Edge_proj_new);
			Edge_samp_old = Edge_samp_new;
			vec_cp_mp(Edge_proj_old,Edge_proj_new);
			//printVec_mp(stdout,0,Edge_proj_old);
			clear_vec_mp(Edge_proj_new);
			refine_old = refine_new;
			if(num_samp_new == num_samp_old) // if had no new samples???
			{
				S_new->vertices_mp[ii] = Edge_samp_old;
				S_new->num_pts[ii] = num_samp_old;
				init_vec_mp(S_new->proj_vertices_mp[ii],num_samp_old);
				vec_cp_mp(S_new->proj_vertices_mp[ii],Edge_proj_old);
				clear_vec_mp(Edge_proj_old);
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

void set_witness_set_d(witness_set_d *W, vec_d L,vec_d pts,int num_vars)
{
	
	
	W->num_variables = num_vars;
	
	
	W->W_mp.num_pts=1;
	W->W_mp.pts=(point_mp *)bmalloc(sizeof(point_mp)); // apparently you can only pass in a single point to copy in.
	
	
	W->W.num_pts=1;
	W->W.pts=(point_d *)bmalloc(sizeof(point_d)); // apparently you can only pass in a single point to copy in.
	
	
	//initialize the memory
	init_point_d(W->W.pts[0],num_vars); W->W.pts[0]->size = num_vars;
	init_point_mp(W->W_mp.pts[0],num_vars); W->W_mp.pts[0]->size = num_vars;
	
	
	point_cp_d(W->W.pts[0],pts);
	point_d_to_mp(W->W_mp.pts[0],pts);
	
	W->num_linears = 1;
	W->L = (vec_d *)bmalloc(sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(sizeof(vec_mp));
	
	init_vec_d( W->L[0],   num_vars); W->L[0]->size = num_vars;
	init_vec_mp(W->L_mp[0],num_vars); W->L_mp[0]->size = num_vars;
	
	vec_cp_d(W->L[0],L);
	vec_d_to_mp(W->L_mp[0],L);
	
}

void set_witness_set_mp(witness_set_d *W, vec_mp L,vec_mp pts,int num_vars)
{
	
	
	W->num_variables = num_vars;
	
	
	W->W_mp.num_pts=1;
	W->W_mp.pts=(point_mp *)bmalloc(sizeof(point_mp)); // apparently you can only pass in a single point to copy in.
	
	
	W->W.num_pts=1;
	W->W.pts=(point_d *)bmalloc(sizeof(point_d)); // apparently you can only pass in a single point to copy in.
	
	
	//initialize the memory
	init_point_d(W->W.pts[0],num_vars); W->W.pts[0]->size = num_vars;
	init_point_mp(W->W_mp.pts[0],num_vars); W->W_mp.pts[0]->size = num_vars;
	
	
	point_cp_mp(W->W_mp.pts[0],pts);
	point_mp_to_d(W->W.pts[0],pts);
	
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
				print_d(OUT,0,&(S.proj_vertices[ii]->coord[jj]));
			else
				print_mp(OUT,0,&(S.proj_vertices_mp[ii]->coord[jj]));
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
	
	comp_d temp;
	comp_mp temp_mp;
	S->num_edges = C.num_edges;
	S->num_pts = (int *)bmalloc(C.num_edges * sizeof(int));
	for(ii=0;ii<C.num_edges;ii++)
		S->num_pts[ii]=3;
	if(MPType==0)
	{
		S->vertices = (vec_d **)bmalloc(C.num_edges * sizeof(vec_d*));
		S->refine = (int **)bmalloc(C.num_edges * sizeof(int*));
		S->proj_vertices = (vec_d *)bmalloc(C.num_edges * sizeof(vec_d));
	}
	else
	{
		init_mp(temp_mp);
		S->vertices_mp = (vec_mp **)bmalloc(C.num_edges * sizeof(vec_mp*));
		S->refine = (int **)bmalloc(C.num_edges * sizeof(int*));
		S->proj_vertices_mp = (vec_mp *)bmalloc(C.num_edges * sizeof(vec_mp));
	}
	for(ii=0;ii<S->num_edges;ii++)
	{
		if(MPType==0)
		{
			S->vertices[ii] = (vec_d *)bmalloc(S->num_pts[ii] * sizeof(vec_d));
			S->refine[ii] = (int *)bmalloc(S->num_pts[ii] * sizeof(int));
			
			for(jj=0;jj<S->num_pts[ii];jj++)
			{
				init_vec_d(S->vertices[ii][jj],num_vars);//three points: left mid right
				S->refine[ii][jj]=1;
			}
			
			//left point
			init_vec_d(S->proj_vertices[ii],S->num_pts[ii]);
			
			
			vec_d pi_nohom;  init_vec_d(pi_nohom,num_vars-1); pi_nohom->size = num_vars-1;
			int mm;
			for (mm=0; mm<num_vars-1; mm++) {
				set_d(&pi_nohom->coord[mm],& C.edges[ii].pi->coord[mm]);
			}
			
			vec_d dehom; init_vec_d(dehom,num_vars-1); dehom->size = num_vars-1;
			
			
			//left point
			index = C.edges[ii].left;
			vec_cp_d(S->vertices[ii][0],C.V1[index].pt);
			dehomogenize(&dehom,C.V1[index].pt);
			dot_product_d(&(S->proj_vertices[ii]->coord[0]), dehom,pi_nohom);
			
			
			//mid point
			vec_cp_d(S->vertices[ii][1],C.edges[ii].midpt);
			dehomogenize(&dehom,C.edges[ii].midpt);
			dot_product_d(&(S->proj_vertices[ii]->coord[1]), dehom,pi_nohom);
			
			
			//right point
			index = C.edges[ii].right;
			vec_cp_d(S->vertices[ii][2],C.V1[index].pt);
			dehomogenize(&dehom,C.V1[index].pt);
			dot_product_d(&(S->proj_vertices[ii]->coord[1]), dehom,pi_nohom);
			
			clear_vec_d(pi_nohom);
			clear_vec_d(dehom);
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
			
			init_vec_mp(S->proj_vertices_mp[ii],S->num_pts[ii]); S->proj_vertices_mp[ii]->size=S->num_pts[ii];
			
			
			vec_mp pi_nohom;  init_vec_mp(pi_nohom,num_vars-1); pi_nohom->size = num_vars-1;
			int mm;
			for (mm=0; mm<num_vars-1; mm++) {
				set_mp(&pi_nohom->coord[mm],& C.edges[ii].pi_mp->coord[mm]);
			}
			
			vec_mp dehom; init_vec_mp(dehom,num_vars-1); dehom->size = num_vars-1;
			
			
			//left point
			index = C.edges[ii].left;
			vec_cp_mp(S->vertices_mp[ii][0],C.V1[index].pt_mp);
			dehomogenize_mp(&dehom,C.V1[index].pt_mp);
			dot_product_mp(&(S->proj_vertices_mp[ii]->coord[0]), dehom,pi_nohom);

			
			//mid point
			vec_cp_mp(S->vertices_mp[ii][1],C.edges[ii].midpt_mp);
			dehomogenize_mp(&dehom,C.edges[ii].midpt_mp);
			dot_product_mp(&(S->proj_vertices_mp[ii]->coord[1]), dehom,pi_nohom);
			
			
			//right point
			index = C.edges[ii].right;
			vec_cp_mp(S->vertices_mp[ii][2],C.V1[index].pt_mp);
			dehomogenize_mp(&dehom,C.V1[index].pt_mp);
			dot_product_mp(&(S->proj_vertices_mp[ii]->coord[1]), dehom,pi_nohom);

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
	
	printf("\n Sampler module for BertiniReal(TM) v%s\n\n", BERTINI_REAL_VERSION_STRING);
	printf(" D.J. Bates, D. Brake,\n W. Hao, J.D. Hauenstein,\n A.J. Sommese, C.W. Wampler\n\n");
	printf("(using GMP v%d.%d.%d, MPFR v%s)\n\n",
				 __GNU_MP_VERSION, __GNU_MP_VERSION_MINOR, __GNU_MP_VERSION_PATCHLEVEL, mpfr_get_version());
	
	// check for write privilege
	if (checkWritePrivilege())
	{
		printf("ERROR: BertiniReal does not have write privileges!\n");
		bexit(ERROR_WRITE_PRIVILEGE);
	}
	
	if (argC > 1 && args[1] != NULL && (!strcmp(args[1], "--help") || !strcmp(args[1], "-help"))) // help
	{ // print information about Bertini
		printf("\nThis is a sampler module for BertiniReal v%s, developed by \n Dan J. Bates, Daniel Brake,\n Wenrui Hao, Jonathan D. Hauenstein, \n Andrew J. Sommmese, and Charles W. Wampler.\n\n", BERTINI_REAL_VERSION_STRING);
		printf("See ??? for details about BertiniReal.\n\n");
	}
	
	
	//setup the name of directory
	IN = safe_fopen_read("Dir_Name");
	fscanf(IN, "%d\n", &strLength);
	directoryName = (char *)bmalloc(strLength * sizeof(char));
	fgets(directoryName, strLength, IN);
	fscanf(IN, "%d\n", &MPType);
	fclose(IN);
	
	//setup E structure from E.edge
	sprintf(tmp_file,  "%s/E.edge", directoryName);
	C->num_edges = setup_edges(&(C->edges),tmp_file,num_vars,inputName,directoryName,MPType);
	
	//setup V0 structure from V0.vert
	sprintf(tmp_file,  "%s/V0.vert", directoryName);
	C->num_V0 = setup_vertices(&(C->V0),tmp_file,*num_vars,MPType);
	
	//setup V1 structure from V1.vert
	sprintf(tmp_file,  "%s/V1.vert", directoryName);
	C->num_V1 = setup_vertices(&(C->V1),tmp_file,*num_vars,MPType);
	
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

int setup_edges(edge_d **edges,char *INfile,int *num_vars, char **inputName, char *directoryName, int MPType)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_edges,ii,jj,strLength;
	char *input_deflated_Name=NULL;
	
	fscanf(IN, "%d\n", num_vars);
	fscanf(IN, "%d\n", &num_edges);
	fscanf(IN, "%d\n", &strLength);
	
	input_deflated_Name = (char *)bmalloc((strLength+1) * sizeof(char));
	fgets(input_deflated_Name, (strLength+1), IN);
	
	strLength = 1 + snprintf(NULL, 0, "%s/%s", directoryName,input_deflated_Name);
	*inputName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*inputName,  "%s/%s", directoryName,input_deflated_Name);
	
	*edges=(edge_d*) bmalloc(num_edges*sizeof(edge_d));
	for(ii=0;ii<num_edges;ii++)
	{
		fscanf(IN,"%d\n",&((*edges)[ii].left));
		fscanf(IN,"%d\n",&((*edges)[ii].right));
		if(MPType==0)
		{
			init_point_d((*edges)[ii].midpt,*num_vars);
			(*edges)[ii].midpt->size=*num_vars;
			for(jj=0;jj<*num_vars;jj++)
				fscanf(IN, "%lf %lf", &((*edges)[ii].midpt->coord[jj].r), &((*edges)[ii].midpt->coord[jj].i));
			
			init_point_d((*edges)[ii].pi,*num_vars);
			(*edges)[ii].pi->size=*num_vars;
			for(jj=0;jj<*num_vars;jj++)
				fscanf(IN, "%lf %lf", &((*edges)[ii].pi->coord[jj].r), &((*edges)[ii].pi->coord[jj].i));
		}
		else
		{
			init_point_mp((*edges)[ii].midpt_mp,*num_vars);
			(*edges)[ii].midpt_mp->size=*num_vars;
			for(jj=0;jj<*num_vars;jj++)
			{
				mpf_inp_str((*edges)[ii].midpt_mp->coord[jj].r, IN, 10);
				mpf_inp_str((*edges)[ii].midpt_mp->coord[jj].i, IN, 10);
			}
			init_point_mp((*edges)[ii].pi_mp,*num_vars);
			(*edges)[ii].pi_mp->size=*num_vars;
			for(jj=0;jj<*num_vars;jj++)
			{
				mpf_inp_str((*edges)[ii].pi_mp->coord[jj].r, IN, 10);
				mpf_inp_str((*edges)[ii].pi_mp->coord[jj].i, IN, 10);
			}
			
			
			printf("ii = %d\n",ii);
			printVec_mp(stdout,0,(*edges)[ii].pi_mp);
			printVec_mp(stdout,0,(*edges)[ii].midpt_mp);
		}
	}
	fclose(IN);
	return num_edges;
}


int setup_vertices(vertex_d **vertices,char *INfile,int num_vars, int MPType)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num_vertices,ii,jj;
	fscanf(IN, "%d\n\n", &num_vertices);
	*vertices=(vertex_d* )bmalloc(num_vertices*sizeof(vertex_d));
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
		}
		else
		{
			init_point_mp((*vertices)[ii].pt_mp,num_vars); (*vertices)[ii].pt_mp->size=num_vars;
			
			for(jj=0;jj<num_vars;jj++)
			{
				mpf_inp_str((*vertices)[ii].pt_mp->coord[jj].r, IN, 10);
				mpf_inp_str((*vertices)[ii].pt_mp->coord[jj].i, IN, 10);
			}
			
		}
	}
	fclose(IN);
	return num_vertices;
}
