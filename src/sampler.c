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
	
	if(setup_curveDecomp(argC, args, &inputName, &witnessSetName,&RandMatName,&samplingNamenew,&C,&num_vars))
		return 1;
	//< prints the welcome message, also gets the inputName, witnessSetName and C
	srand(time(NULL));
	// essentials for using the bertini parser
	prog_t SLP;
	
	unsigned int currentSeed;
	int trackType, genType = 0, MPType,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0;
	int my_id, num_processes, headnode = 0; // headnode is always 0
	int precision = 53;
	num_processes = 1;
	int num_var_gps = 0, userHom = 0;
	mat_mp n_minusone_randomizer_matrix;
	
	//end parser-bertini essentials
	
	parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
	preproc_data PPD;
	setupPreProcData("preproc_data", &PPD);
	num_var_gps = PPD.num_var_gp;
	num_vars = setupProg(&SLP, precision, MPType); // num_vars includes the number of homogeneous coordinates.
	// the number of homogeneous coordinates is the num_var_gps.
	
	tracker_config_t T;
	get_tracker_config(&T,MPType);
	
	initMP(T.Precision);
	
	printf("parsing witness set\n");
	witnessSetParse(&W, witnessSetName,num_vars);
	read_rand_matrix(RandMatName, n_minusone_randomizer_matrix);
	printf("Loading sampling data\n");
	Load_sampling_data(&S_old,C,num_vars,MPType);
	
	//Generate new sampling data
	generate_new_sampling_pts(&S_new,n_minusone_randomizer_matrix,S_old, C,W, num_vars,currentSeed,MPType);
	printf("generate_new_sampling_pts\n");
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
	//TMP END
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
	int            i,j,k;
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
	
	for(i=0;i<S_old.num_edges;i++)
	{
		Edge_samp_old = S_old.vertices[i];
		num_samp_old = S_old.num_pts[i];
		init_vec_d(Edge_proj_old,S_old.num_pts[i]);
		vec_cp_d(Edge_proj_old, S_old.proj_vertices[i]);
		refine_old = S_old.refine[i];
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
			for(j=0;j<num_samp_old-1;j++)
			{
				if(refine_old[j])
				{
					vec_cp_d(L,C.edges[i].pi);
					point_d_to_mp(new_linears,L);
					if(j==0)
					{
						vec_cp_d(startpt,Edge_samp_old[j+1]);
						sub_d(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[j+1]));
					}
					else
					{
						vec_cp_d(startpt,Edge_samp_old[j]);
						sub_d(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[j]));
					}
					//tracking end point pi_end = \pi (si+s_{i+1})/2
					//w=2w1w2/(w1+w2)
					add_d(&(mid_pt->coord[0]),&(Edge_samp_old[j]->coord[0]),
								&(Edge_samp_old[j+1]->coord[0]));
					mul_d(temp,&(Edge_samp_old[j]->coord[0]),
								&(Edge_samp_old[j+1]->coord[0]));
					div_d(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					temp->r=2.0; temp->i=0.0;
					mul_d(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mul_d(pi_end,&(mid_pt->coord[0]),&(C.edges[i].pi->coord[0]));
					
					for(k=1;k<num_vars;k++)
					{
						div_d(&(mid_pt->coord[k]),&(Edge_samp_old[j]->coord[k]),
									&(Edge_samp_old[j]->coord[0]));
						div_d(temp,&(Edge_samp_old[j+1]->coord[k]),
									&(Edge_samp_old[j+1]->coord[0]));
						add_d(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						temp->r=0.5; temp->i=0.0;
						mul_d(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						mul_d(&(mid_pt->coord[k]),&(mid_pt->coord[0]),&(mid_pt->coord[k]));
						mul_d(temp,&(mid_pt->coord[k]),&(C.edges[i].pi->coord[k]));
						add_d(pi_end,pi_end,temp);
					}
					div_d(pi_end,pi_end,&(mid_pt->coord[0]));
					sub_d(temp,&(C.edges[i].pi->coord[0]),pi_end);
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
					set_d(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[j+1]));
					init_vec_d(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_d(Edge_samp_new[num_samp_new],Edge_samp_old[j+1]);
					num_samp_new++;
					
				}
				else
				{
					refine_new[num_samp_new-1] = 0;
					set_d(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[j+1]));
					init_vec_d(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_d(Edge_samp_new[num_samp_new],Edge_samp_old[j+1]);
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
				S_new->vertices[i] = Edge_samp_old;
				S_new->num_pts[i] = num_samp_old;
				init_vec_d(S_new->proj_vertices[i],num_samp_old);
				vec_cp_d(S_new->proj_vertices[i],Edge_proj_old);
				clear_vec_d(Edge_proj_old);
				S_new->refine[i] = refine_old;
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



void generate_new_sampling_pts_mp(sample_d *S_new,mat_mp n_minusone_randomizer_matrix,sample_d S_old, curveDecomp_d C,witness_set_d W, int num_vars,unsigned int currentSeed, int MPType)
{
	witness_set_d  Wnew;
	vec_mp         new_linears;
	vec_mp          L,startpt,mid_pt;
	int            i,j,k;
	comp_mp         pi_end,temp,temp1;
	vec_mp          *Edge_samp_old, *Edge_samp_new, Edge_proj_old, Edge_proj_new;
	int            num_samp_old, num_samp_new;
	int            *refine_old, *refine_new;
	double         max_norm, TOL=1e-1;
	
	init_witness_set_d(&Wnew);
	cp_patches(&Wnew,W);
	
	init_vec_mp(new_linears,num_vars);
	new_linears->size = num_vars;
	
	init_vec_mp(L,num_vars); // what is the purpose of this?
	L->size = num_vars; // what is L?
	
	init_vec_mp(startpt,num_vars);
	startpt->size = num_vars;
	
	init_vec_mp(mid_pt,num_vars);
	mid_pt->size = num_vars;
	
	init_mp(temp); init_mp(pi_end); init_mp(temp1);
	S_new->num_edges = S_old.num_edges;
	S_new->vertices_mp = (vec_mp **)bmalloc(S_new->num_edges * sizeof(vec_mp*));
	S_new->proj_vertices_mp = (vec_mp *)bmalloc(S_new->num_edges * sizeof(vec_mp));
	S_new->refine = (int **)bmalloc(S_new->num_edges * sizeof(int*));
	S_new->num_pts = (int *) bmalloc(S_new->num_edges * sizeof(int));
	for(i=0;i<S_old.num_edges;i++)
	{
		Edge_samp_old = S_old.vertices_mp[i];
		num_samp_old = S_old.num_pts[i];
		init_vec_mp(Edge_proj_old,S_old.num_pts[i]);
		vec_cp_mp(Edge_proj_old, S_old.proj_vertices_mp[i]);
		refine_old = S_old.refine[i];
		while(1)
		{
			Edge_samp_new = (vec_mp *)bmalloc(2*num_samp_old * sizeof(vec_mp));
			init_vec_mp(Edge_proj_new,2*num_samp_old); // what is pV?
			Edge_proj_new->size=2*num_samp_old;
			refine_new = (int * )bmalloc(2*num_samp_old * sizeof(int));
			//copy left point
			set_mp(&(Edge_proj_new->coord[0]),&(Edge_proj_old->coord[0]));
			num_samp_new = 1;
			init_vec_mp(Edge_samp_new[0],num_vars);
			vec_cp_mp(Edge_samp_new[0],Edge_samp_old[0]);
			for(j=0;j<num_samp_old-1;j++)
			{
				if(refine_old[j])
				{
					vec_cp_mp(L,C.edges[i].pi_mp);
					vec_cp_mp(new_linears,L);
					if(j==0)
					{
						vec_cp_mp(startpt,Edge_samp_old[j+1]);
						sub_mp(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[j+1]));
					}
					else
					{
						vec_cp_mp(startpt,Edge_samp_old[j]);
						sub_mp(&(L->coord[0]),&(L->coord[0]),&(Edge_proj_old->coord[j]));
					}
					
					add_mp(&(mid_pt->coord[0]),&(Edge_samp_old[j]->coord[0]),
								 &(Edge_samp_old[j+1]->coord[0]));
					mul_mp(temp,&(Edge_samp_old[j]->coord[0]),
								 &(Edge_samp_old[j+1]->coord[0]));
					div_mp(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mpf_set_d(temp->r, 2.0); mpf_set_d(temp->i, 0.0);
					mul_mp(&(mid_pt->coord[0]),temp,&(mid_pt->coord[0]));
					mul_mp(pi_end,&(mid_pt->coord[0]),&(C.edges[i].pi_mp->coord[0]));
					for(k=1;k<num_vars;k++)
					{
						div_mp(&(mid_pt->coord[k]),&(Edge_samp_old[j]->coord[k]),
									 &(Edge_samp_old[j]->coord[0]));
						div_mp(temp,&(Edge_samp_old[j+1]->coord[k]),
									 &(Edge_samp_old[j+1]->coord[0]));
						add_mp(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						
						mpf_set_d(temp->r, 0.5); mpf_set_d(temp->i, 0.0);
						mul_mp(&(mid_pt->coord[k]),temp,&(mid_pt->coord[k]));
						mul_mp(&(mid_pt->coord[k]),&(mid_pt->coord[0]),&(mid_pt->coord[k]));
						mul_mp(temp,&(mid_pt->coord[k]),&(C.edges[i].pi_mp->coord[k]));
						add_mp(pi_end,pi_end,temp);
					}
					div_mp(pi_end,pi_end,&(mid_pt->coord[0]));
					sub_mp(temp,&(C.edges[i].pi_mp->coord[0]),pi_end);
					
					set_mp(&(new_linears->coord[0]),temp);
					set_witness_set_mp(&W, L,startpt,num_vars);
					lin_to_lin_solver_main(MPType,
																 W,         // witness_set
																 n_minusone_randomizer_matrix,
																 &new_linears, //  the set of linears we will solve at.
																 1, // the number of new linears.
																 &Wnew); // the new data is put here!
					set_zero_mp(&(Edge_proj_new->coord[num_samp_new]));
					init_vec_mp(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_mp(Edge_samp_new[num_samp_new],Wnew.W_mp.pts[0]);
					max_norm=0.0;
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
					set_mp(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[j+1]));
					//printf("Edge_proj_new=");print_mp(stdout,0,&(Edge_proj_new->coord[num_samp_new]));
					init_vec_mp(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_mp(Edge_samp_new[num_samp_new],Edge_samp_old[j+1]);
					num_samp_new++;
					
				}
				else
				{
					refine_new[num_samp_new-1] = 0;
					set_mp(&(Edge_proj_new->coord[num_samp_new]),&(Edge_proj_old->coord[j+1]));
					init_vec_mp(Edge_samp_new[num_samp_new],num_vars);
					vec_cp_mp(Edge_samp_new[num_samp_new],Edge_samp_old[j+1]);
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
			if(num_samp_new == num_samp_old)
			{
				S_new->vertices_mp[i] = Edge_samp_old;
				S_new->num_pts[i] = num_samp_old;
				init_vec_mp(S_new->proj_vertices_mp[i],num_samp_old);
				vec_cp_mp(S_new->proj_vertices_mp[i],Edge_proj_old);
				clear_vec_mp(Edge_proj_old);
				S_new->refine[i] = refine_old;
				break;
			}
			else
				num_samp_old=num_samp_new;
		}
	}
	clear_mp(temp); clear_mp(pi_end); clear_mp(temp1);
	
	clear_witness_set(Wnew);
	clear_vec_mp(L);
	clear_vec_mp(new_linears);
}


void read_rand_matrix(char *INfile, mat_mp n_minusone_randomizer_matrix)
{
	int i,j,rows,cols,cur_precision;
	FILE *IN= safe_fopen_read(INfile);
	fscanf(IN,"%d\n",&rows);
	fscanf(IN,"%d\n",&cols);
	fscanf(IN,"%d\n",&cur_precision);
	init_mat_mp2(n_minusone_randomizer_matrix,0,0,cur_precision);
	increase_size_mat_mp(n_minusone_randomizer_matrix, rows, cols);
	n_minusone_randomizer_matrix->rows=rows;
	n_minusone_randomizer_matrix->cols=cols;
	for(i=0;i<n_minusone_randomizer_matrix->rows;i++)
		for(j=0;j<n_minusone_randomizer_matrix->cols;j++)
		{
			fscanf(IN,"<");
			mpf_inp_str(n_minusone_randomizer_matrix->entry[i][j].r, IN, 10);
			mpf_inp_str(n_minusone_randomizer_matrix->entry[i][j].i, IN, 10);
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
	int i,j,k;
	// output the number of vertices
	fprintf(OUT,"%d\n\n",S.num_edges);
	for(i=0;i<S.num_edges;i++)
	{
		fprintf(OUT,"%d\n\n",S.num_pts[i]);
		for(j=0;j<S.num_pts[i];j++)
		{
			if(MPType==0)
				print_d(OUT,0,&(S.proj_vertices[i]->coord[j]));
			else
				print_mp(OUT,0,&(S.proj_vertices_mp[i]->coord[j]));
			fprintf(OUT,"\n");
			for(k=0;k<num_vars;k++)
			{
				if(MPType==0)
					print_d(OUT, 0, &(S.vertices[i][j]->coord[k]));
				else
					print_mp(OUT, 0, &(S.vertices_mp[i][j]->coord[k]));
				
				fprintf(OUT,"\n");
			}
			fprintf(OUT,"\n");
		}
	}
}

int  Load_sampling_data(sample_d *S, curveDecomp_d C,int num_vars,int MPType)
{
	int i,j,index;
	
	comp_d temp;
	comp_mp temp_mp;
	S->num_edges = C.num_edges;
	S->num_pts = (int *)bmalloc(C.num_edges * sizeof(int));
	for(i=0;i<C.num_edges;i++)
		S->num_pts[i]=3;
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
	for(i=0;i<S->num_edges;i++)
	{
		if(MPType==0)
		{
			S->vertices[i] = (vec_d *)bmalloc(S->num_pts[i] * sizeof(vec_d));
			S->refine[i] = (int *)bmalloc(S->num_pts[i] * sizeof(int));
			for(j=0;j<S->num_pts[i];j++)
			{
				init_vec_d(S->vertices[i][j],num_vars);//three points: left mid right
				S->refine[i][j]=1;
			}
			index = C.edges[i].left;
			//left point
			init_vec_d(S->proj_vertices[i],S->num_pts[i]);
			set_zero_d(&(S->proj_vertices[i]->coord[0]));
			vec_cp_d(S->vertices[i][0],C.V1[index].pt);
			for(j=0;j<num_vars;j++)
			{
				mul_d(temp,&(S->vertices[i][0]->coord[j]),&(C.edges[i].pi->coord[j]));
				add_d(&(S->proj_vertices[i]->coord[0]),&(S->proj_vertices[i]->coord[0]),temp);
			}
			//mid point
			set_zero_d(&(S->proj_vertices[i]->coord[1]));
			vec_cp_d(S->vertices[i][1],C.edges[i].midpt);
			
			for(j=0;j<num_vars;j++)
			{
				mul_d(temp,&(S->vertices[i][1]->coord[j]),&(C.edges[i].pi->coord[j]));
				add_d(&(S->proj_vertices[i]->coord[1]),&(S->proj_vertices[i]->coord[1]),temp);
			}
			//right point
			set_zero_d(&(S->proj_vertices[i]->coord[2]));
			index = C.edges[i].right;
			vec_cp_d(S->vertices[i][2],C.V1[index].pt);
			
			for(j=0;j<num_vars;j++)
			{
				mul_d(temp,&(S->vertices[i][2]->coord[j]),&(C.edges[i].pi->coord[j]));
				add_d(&(S->proj_vertices[i]->coord[2]),&(S->proj_vertices[i]->coord[2]),temp);
			}
		}
		else
		{
			S->vertices_mp[i] = (vec_mp *)bmalloc(S->num_pts[i] * sizeof(vec_mp));
			S->refine[i] = (int *)bmalloc(S->num_pts[i] * sizeof(int));
			for(j=0;j<S->num_pts[i];j++)
			{
				init_vec_mp(S->vertices_mp[i][j],num_vars);//three points: left mid right
				S->refine[i][j]=1;
			}
			index = C.edges[i].left;
			//left point
			init_vec_mp(S->proj_vertices_mp[i],S->num_pts[i]);
			S->proj_vertices_mp[i]->size=S->num_pts[i];
			set_zero_mp(&(S->proj_vertices_mp[i]->coord[0]));
			vec_cp_mp(S->vertices_mp[i][0],C.V1[index].pt_mp);
			for(j=0;j<num_vars;j++)
			{
				mul_mp(temp_mp,&(S->vertices_mp[i][0]->coord[j]),&(C.edges[i].pi_mp->coord[j]));
				add_mp(&(S->proj_vertices_mp[i]->coord[0]),&(S->proj_vertices_mp[i]->coord[0]),temp_mp);
			}
			//mid point
			set_zero_mp(&(S->proj_vertices_mp[i]->coord[1]));
			vec_cp_mp(S->vertices_mp[i][1],C.edges[i].midpt_mp);
			
			for(j=0;j<num_vars;j++)
			{
				mul_mp(temp_mp,&(S->vertices_mp[i][1]->coord[j]),&(C.edges[i].pi_mp->coord[j]));
				add_mp(&(S->proj_vertices_mp[i]->coord[1]),&(S->proj_vertices_mp[i]->coord[1]),temp_mp);
			}
			//right point
			set_zero_mp(&(S->proj_vertices_mp[i]->coord[2]));
			index = C.edges[i].right;
			vec_cp_mp(S->vertices_mp[i][2],C.V1[index].pt_mp);
			
			for(j=0;j<num_vars;j++)
			{
				mul_mp(temp_mp,&(S->vertices_mp[i][2]->coord[j]),&(C.edges[i].pi_mp->coord[j]));
				add_mp(&(S->proj_vertices_mp[i]->coord[2]),&(S->proj_vertices_mp[i]->coord[2]),temp_mp);
			}
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
	int num,i,j,strLength;
	char *input_deflated_Name=NULL;
	
	fscanf(IN, "%d\n", num_vars);
	fscanf(IN, "%d\n", &num);
	fscanf(IN, "%d\n", &strLength);
	
	input_deflated_Name = (char *)bmalloc((strLength+1) * sizeof(char));
	fgets(input_deflated_Name, (strLength+1), IN);
	
	strLength = 1 + snprintf(NULL, 0, "%s/%s", directoryName,input_deflated_Name);
	*inputName = (char *)bmalloc(strLength * sizeof(char));
	sprintf(*inputName,  "%s/%s", directoryName,input_deflated_Name);
	
	*edges=(edge_d*) bmalloc(num*sizeof(edge_d));
	for(i=0;i<num;i++)
	{
		fscanf(IN,"%d\n",&((*edges)[i].left));
		fscanf(IN,"%d\n",&((*edges)[i].right));
		if(MPType==0)
		{
			init_point_d((*edges)[i].midpt,*num_vars);
			(*edges)[i].midpt->size=*num_vars;
			for(j=0;j<*num_vars;j++)
				fscanf(IN, "%lf %lf", &((*edges)[i].midpt->coord[j].r), &((*edges)[i].midpt->coord[j].i));
			
			init_point_d((*edges)[i].pi,*num_vars);
			(*edges)[i].pi->size=*num_vars;
			for(j=0;j<*num_vars;j++)
				fscanf(IN, "%lf %lf", &((*edges)[i].pi->coord[j].r), &((*edges)[i].pi->coord[j].i));
		}
		else
		{
			init_point_mp((*edges)[i].midpt_mp,*num_vars);
			(*edges)[i].midpt_mp->size=*num_vars;
			for(j=0;j<*num_vars;j++)
			{
				mpf_inp_str((*edges)[i].midpt_mp->coord[j].r, IN, 10);
				mpf_inp_str((*edges)[i].midpt_mp->coord[j].i, IN, 10);
			}
			init_point_mp((*edges)[i].pi_mp,*num_vars);
			(*edges)[i].pi_mp->size=*num_vars;
			for(j=0;j<*num_vars;j++)
			{
				mpf_inp_str((*edges)[i].pi_mp->coord[j].r, IN, 10);
				mpf_inp_str((*edges)[i].pi_mp->coord[j].i, IN, 10);
			}
			printf("i = %d",i);
			printVec_mp(stdout,0,(*edges)[i].pi_mp);
			printVec_mp(stdout,0,(*edges)[i].midpt_mp);
		}
	}
	fclose(IN);
	return num;
}


int setup_vertices(vertex_d **vertices,char *INfile,int num_vars, int MPType)
//setup the vertex structure
{
	FILE *IN = safe_fopen_read(INfile);
	int num,i,j;
	fscanf(IN, "%d\n\n", &num);
	*vertices=(vertex_d* )bmalloc(num*sizeof(vertex_d));
	for(i=0;i<num;i++)
	{
		if(MPType==0)
		{
			init_point_d((*vertices)[i].pt,num_vars);
			(*vertices)[i].pt->size=num_vars;
			for(j=0;j<num_vars;j++)
			{
				fscanf(IN, "%lf %lf", &((*vertices)[i].pt->coord[j].r), &((*vertices)[i].pt->coord[j].i));
			}
		}
		else
		{
			init_point_mp((*vertices)[i].pt_mp,num_vars);
			(*vertices)[i].pt_mp->size=num_vars;
			for(j=0;j<num_vars;j++)
			{
				mpf_inp_str((*vertices)[i].pt_mp->coord[j].r, IN, 10);
				mpf_inp_str((*vertices)[i].pt_mp->coord[j].i, IN, 10);
			}
			
		}
	}
	fclose(IN);
	return num;
}
