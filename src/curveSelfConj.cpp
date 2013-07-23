#include "curveSelfConj.hpp"










void computeCurveSelfConj(boost::filesystem::path inputFile,
													witness_set & W,
													vec_mp *pi,
													curve_decomposition &C,
													vertex_set &V,
													int num_vars,
													BR_configuration *program_options,
													solver_configuration *solve_options)
{
	//IN DEVELOPMENT
  
	
	
	int ambient_dim = 1;
	
	int ii,kk;
	
	witness_set Wtemp,Wtemp2;
	
	
	int *declarations = NULL;
	partition_parse(&declarations, inputFile, "func_input", "config", 0); // the 0 means not self conjugate.
	
	
	

	
	
	
	
	// 2) randomize down to N-1 equations
	// to get a square system for the homotopies in the following steps.
	
	
	//create the matrix
	mat_mp randomizer_matrix;
	init_mat_mp2(randomizer_matrix,W.num_variables-W.num_patches-ambient_dim,solve_options->PPD.num_funcs,solve_options->T.AMP_max_prec);
	
	//create the array of integers
	int *randomized_degrees = (int *)br_malloc((W.num_variables-W.num_patches-ambient_dim)*sizeof(int));
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(randomizer_matrix, &randomized_degrees, W.num_variables-W.num_patches-ambient_dim, solve_options->PPD.num_funcs);
	
	if (program_options->verbose_level>=2)
		print_matrix_to_screen_matlab(randomizer_matrix,"randomization");
	
	
	
	
	
	
	
	FILE *OUT = safe_fopen_write("Rand_matrix");
	fprintf(OUT,"%d %d\n",randomizer_matrix->rows,randomizer_matrix->cols);
	fprintf(OUT,"%d\n",solve_options->T.AMP_max_prec);
	fprintf(OUT,"\n");
	print_matrix_to_file_mp(OUT, 0, randomizer_matrix);
	fclose(OUT);
	
	
		
  // 4) solve for critical conditions for random complex projection
	witness_set W_crit_real;
	cp_names(&W_crit_real, W);
	cp_patches(&W_crit_real, W);
	
	if (program_options->crit_solver == LINPRODTODETJAC)
	{
#ifdef printpathlintolin
		remove("pathtrack_lintolin");
#endif
		
		
		//compute the number of linears we will need to move to, to perform the regeneration for linprodtodetjac method.
		int num_new_linears = 0;
		for (ii=0; ii<W.num_variables-W.num_patches-ambient_dim; ii++) {
			num_new_linears+= randomized_degrees[ii]-1; // subtract 1 for differentiation
		}
		num_new_linears -= 1; // subtract 1 because started with a linear (in the witness_set)
		
		if (program_options->verbose_level>=1) {
			printf("the number of *new* linears will be %d\n",num_new_linears);
		}
		
		
		compute_crit_linprodtodetjac(&W_crit_real, // the returned value
																 W,            // input the original witness set
																 randomizer_matrix,
																 *pi,
																 num_new_linears,
																 program_options,
																 solve_options);
	}
	else	{
		
		nullspace_config ns_config;
		compute_crit_nullspace(&W_crit_real, // the returned value
													 W,            // input the original witness set
													 randomizer_matrix, 
													 pi,
													 randomized_degrees,
													 1,  // dimension of ambient complex object
													 1,   //  target dimension to find
													 1,   // COdimension of the critical set to find.
													 program_options,
													 solve_options,
													 &ns_config);
		ns_config.clear();
		
		cp_patches(&W_crit_real, W);
		cp_linears(&W_crit_real, W);
		W_crit_real.only_first_vars(num_vars);
		W_crit_real.sort_for_real(solve_options->T);
		
		
		W_crit_real.write_dehomogenized_coordinates("W_before_additional");
		
		// now get the bounding box critical points and ends of the interval
		curve_get_additional_critpts(&W_crit_real,
																 W,
																 randomizer_matrix,
																 pi[0],
																 randomized_degrees,
																 program_options,
																 solve_options);
		
	}
	
	
	
	
	W_crit_real.write_dehomogenized_coordinates("W_after_additional");
	
	
	
  // 5) compute h_o, get critical points
	

	
	
	
	vec_mp projection_values; init_vec_mp(projection_values,W_crit_real.num_pts); projection_values->size = W_crit_real.num_pts;
	
	
	for (ii=0; ii<W_crit_real.num_pts; ii++) 
		projection_value_homogeneous_input(&projection_values->coord[ii],W_crit_real.pts_mp[ii], *pi); // set projection value
	
	
	
	if (program_options->verbose_level>=1)
		printf("sorting projection values\n");
	
	
	vec_mp projections_sorted;
	init_vec_mp(projections_sorted,W_crit_real.num_pts); projections_sorted->size = W_crit_real.num_pts;
	
	int *index_tracker = (int *)br_malloc(W_crit_real.num_pts*sizeof(vec_mp)*sizeof(int));
	sort_increasing_by_real(&projections_sorted, &index_tracker, projection_values);
	
	
	
	
	
	///////
	//
	//   actually form crit.
	//
	/////////
	
	
	
	double *crit_downstairs; int num_crit;
	
	num_crit = W_crit_real.num_pts;
	crit_downstairs = (double *)br_malloc(num_crit*sizeof(double));
	
	int num_midpoints = W_crit_real.num_pts-1;
	
	if (num_midpoints<1) {
		printf("no midpoints to work with :(\n");
		printf("please program exiting setting C (start at line 449 in curveSelfConj() )\n");
		exit(-1);
	}
	else{
		printf("%d midpoints downstairs\n",num_midpoints);
	}
	
	
	double *midpoints_downstairs = (double *) br_malloc(num_midpoints*sizeof(double));
	
	
	for (ii=0; ii<projections_sorted->size-1; ii++) 
		midpoints_downstairs[ii] = (mpf_get_d(projections_sorted->coord[ii].r)+mpf_get_d(projections_sorted->coord[ii+1].r))/2;
	
	for (ii=0; ii<W_crit_real.num_pts; ii++) 
		crit_downstairs[ii] = mpf_get_d(projections_sorted->coord[ii].r);
	
	
	if (program_options->verbose_level>=1) {
		printf("midpoints:\n");
		for (ii=0; ii<num_midpoints; ii++) {
			printf("%lf\n",midpoints_downstairs[ii]);
		}
		
		printf("crit_downstairs:\n");
		//print out the crit set of projection values.
		for (ii=0; ii<num_crit; ii++) {
			printf("%lf ",crit_downstairs[ii]);
		}
		printf("\n");
	}
	
	
	
	
	
	//copy the points from witness sets into curve decomposition V0;
	
	
	vertex temp_vertex;
	
	for (ii=0; ii< projections_sorted->size; ii++){
		if (program_options->verbose_level>=2)
			printf("adding point %d of %d from crit_real to vertices\n",ii,projections_sorted->size);
		
		set_mp(temp_vertex.projVal_mp,  &projections_sorted->coord[ii]); // set projection value
		vec_cp_mp(temp_vertex.pt_mp,W_crit_real.pts_mp[index_tracker[ii]]);// set point
		temp_vertex.type = CRITICAL; // set type
		C.add_vertex(V,temp_vertex);
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
  // 6) find points P
	// T = \pi(crit)
	// S = midpoints of T
  
	// for ( s \in S)
	//find P
	
	solve_options->allow_multiplicity=1;
	solve_options->allow_singular=1;
	solve_options->use_midpoint_checker = 0;
	
	
	
	int edge_counter = 0; // set the counter
	
	
	std::vector< witness_set> midpoint_witness_sets;
	midpoint_witness_sets.resize(num_midpoints);
	
	
//	witness_set *midpoint_witness_sets;
//	midpoint_witness_sets = (witness_set *)br_malloc(num_midpoints*sizeof(witness_set));
	
	comp_d temp;
  vec_mp particular_projection;  init_vec_mp(particular_projection,W.num_variables); particular_projection->size = W.num_variables;
	vec_cp_mp(particular_projection,*pi);
	
	for (ii=0; ii<num_midpoints; ii++) {
		if (program_options->verbose_level>=2) {
			printf("solving midpoints upstairs %d, projection value %lf\n",ii,midpoints_downstairs[ii]);
		}
		
		

		temp->i = 0;
		temp->r = -midpoints_downstairs[ii];
		d_to_mp(&particular_projection->coord[0],temp);
		
		//make the projection we will solve at.
		
		solve_options->allow_multiplicity = 1;
		solve_options->allow_singular = 1;
		
		lintolin_solver_main(solve_options->T.MPType,
												 W,
												 randomizer_matrix,
												 &particular_projection, 1,
												 &Wtemp,
												 solve_options); //sets this value
		
		
		if (program_options->verbose_level>=2) {
			printf("sorting midpoint witness set %d for realness\n",ii);
		}
		
		sort_for_real(&Wtemp2, Wtemp, solve_options->T);
		sort_for_unique(&midpoint_witness_sets[ii], Wtemp2, solve_options->T);
		edge_counter += midpoint_witness_sets[ii].num_pts;
		
		Wtemp.reset(); Wtemp2.reset();
	}
	
	if (program_options->verbose_level>=1) {
		printf("done finding midpoints upstairs\n");
	}
	
	

	
  // 7) find edge endpoints
  
	solve_options->allow_multiplicity = 1;
	solve_options->allow_singular = 1;
	
	
	witness_set Wleft, Wright;
	
	edge temp_edge;
	comp_mp left_check; init_mp(left_check);
	comp_mp right_check; init_mp(right_check);
	
	for (ii=0; ii<num_midpoints; ++ii) {
		if (program_options->verbose_level>=1) {
			printf("moving from downstairs midpoint %d, to left and right\n",ii);
		}
		
		
		temp->r = crit_downstairs[ii]; temp->i = 0;
		d_to_mp(left_check,temp);
		d_to_mp(&particular_projection->coord[0],temp);
		neg_mp(&particular_projection->coord[0],&particular_projection->coord[0]);
		
		lintolin_solver_main(solve_options->T.MPType,
												 midpoint_witness_sets[ii], //the input
												 randomizer_matrix,
												 &particular_projection, 1,
												 &Wleft,
												 solve_options); // the output
		
		
		temp->r = crit_downstairs[ii+1]; temp->i = 0;
		d_to_mp(right_check,temp);
		d_to_mp(&particular_projection->coord[0],temp);
		neg_mp(&particular_projection->coord[0],&particular_projection->coord[0]);
		
		lintolin_solver_main(solve_options->T.MPType,
												 midpoint_witness_sets[ii], //the input
												 randomizer_matrix,
												 &particular_projection, 1,
												 &Wright,
												 solve_options); // the output
		//each member of Wtemp should real.  if a member of V1 already, mark index.  else, add to V1, and mark.
		
		
		int keep_going = 1;
		if (Wright.num_pts!=midpoint_witness_sets[ii].num_pts) {
			printf("had a critical failure\n moving right was deficient a point\n");
			keep_going = 0;
		}
		if (Wleft.num_pts!=midpoint_witness_sets[ii].num_pts) {
			printf("had a critical failure\n moving left was deficient a point\n");
			keep_going = 0;
		}
		if (!keep_going) {
			exit(3999);
		}
		
		
		for (kk=0; kk<midpoint_witness_sets[ii].num_pts; kk++) {
			neg_mp(temp_vertex.projVal_mp, &midpoint_witness_sets[ii].L_mp[0]->coord[0]  ); // set projection value
			vec_cp_mp(temp_vertex.pt_mp,midpoint_witness_sets[ii].pts_mp[kk]);// set point
			
			
			
			temp_vertex.type = MIDPOINT; // set type
			
			temp_edge.midpt = C.add_vertex(V,temp_vertex); // gets the index of the new midpoint as it is added
			temp_edge.left  = C.index_in_vertices_with_add(V,Wleft.pts_mp[kk],  left_check,  solve_options->T); 
			temp_edge.right = C.index_in_vertices_with_add(V,Wright.pts_mp[kk], right_check, solve_options->T);
			
			C.add_edge(temp_edge);
			
			if (program_options->verbose_level>=3) {
				printf("upstairs midpoint %d\n",kk);
				print_comp_mp_matlab(V.vertices[temp_edge.left].projVal_mp,"left");
				print_comp_mp_matlab(V.vertices[temp_edge.midpt].projVal_mp,"mid");
				print_comp_mp_matlab(V.vertices[temp_edge.right].projVal_mp,"right");
				printf("indices of left, mid, right: %d %d %d\n",temp_edge.left,temp_edge.midpt,temp_edge.right);
				printf("\n\n");
			}
		}
		Wleft.reset();
		Wright.reset();

	}//re: for ii
	clear_mp(left_check); clear_mp(right_check);
	
	if (program_options->verbose_level>=0) {
		printf("C.num_edges = %d\n",C.num_edges);
	}
	
	
	
	
	return;
} // re: computeCurveSelfConj


//subfunctions





int curve_get_additional_critpts(witness_set *W_crit_real,
																 witness_set & W,
																 mat_mp randomizer_matrix,
																 vec_mp pi,
																 int *randomized_degrees,
																 BR_configuration *program_options,
																 solver_configuration *solve_options)
{
	
	
	solve_options->complete_witness_set=1;
	
	int ii, jj;
	
	witness_set Wtemp, Wtemp2;
	
	/////now compute the bounds, whether bounding box or bounding the interval
	
	comp_mp one; init_mp(one);
	set_one_mp(one);
	
	vec_mp variable_projection; init_vec_mp(variable_projection, W.num_variables); variable_projection->size = W.num_variables;
	
	if (!program_options->use_bounding_box) {// do not want bounding box
		
		if (program_options->verbose_level>=0)
			printf("in no-box mode.  computing bounds of interval\n");
		
		
		vec_cp_mp(variable_projection,pi);
		
		// compute projections of W_crit_real
		vec_mp projection_values; init_vec_mp(projection_values,W_crit_real->num_pts); projection_values->size=W_crit_real->num_pts;
		for (jj=0; jj<W_crit_real->num_pts; jj++) {
			projection_value_homogeneous_input(&projection_values->coord[jj], W_crit_real->pts_mp[jj],pi);
		}
		
		vec_mp projections_sorted;
		init_vec_mp(projections_sorted,W_crit_real->num_pts); projections_sorted->size = W_crit_real->num_pts;
		int *index_tracker = (int *)br_malloc(W_crit_real->num_pts*sizeof(int));
		sort_increasing_by_real(&projections_sorted, &index_tracker, projection_values);
		
		printf("%d points in W_crit_real\n", W_crit_real->num_pts);
		print_point_to_screen_matlab(projections_sorted,"proj_sorted");
		
		
		comp_mp left_proj_val, right_proj_val;
		init_mp(left_proj_val); init_mp(right_proj_val);
		
		if (W_crit_real->num_pts==0) {
			set_one_mp(right_proj_val);
			set_one_mp(left_proj_val);
			neg_mp(left_proj_val,left_proj_val);
		}
		else{
			set_mp(left_proj_val, &projections_sorted->coord[0]);
			set_mp(right_proj_val, &projections_sorted->coord[W_crit_real->num_pts-1]);
		}
		
		
		comp_mp relevant_proj_val; init_mp(relevant_proj_val);
		//this loop presupposes that w_crit_real has at least one member
		for (jj=0; jj<2; jj++) {
			
			
			//set the value for the projection to move to
			if (jj==0){ // if on the first of the two iterations
				set_mp(relevant_proj_val,left_proj_val);}
			else{
				set_mp(relevant_proj_val,right_proj_val);}
			
			set_mp(&variable_projection->coord[0],relevant_proj_val);
			
			if (jj==0){
				sub_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);}
			else{
				add_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);}
			
			neg_mp(&variable_projection->coord[0],&variable_projection->coord[0]);
			//done setting the value of the projection.
			
			//run the solver
			lintolin_solver_main(solve_options->T.MPType,
													 W,
													 randomizer_matrix,
													 &variable_projection, 1,
													 &Wtemp,
													 solve_options);
			
			Wtemp.sort_for_real(solve_options->T);
			Wtemp.sort_for_unique(solve_options->T);
			W_crit_real->merge(Wtemp);
			
			Wtemp.reset();
		}
		
		clear_mp(relevant_proj_val);
	}
	else //want a bounding box
	{
		
		if (program_options->verbose_level>=0) {
			printf("computing bounding box\n");
		}
		// for each variable direction, compute projections of crit points onto that axis.
		
		
		for (ii=0; ii<W.num_variables-1; ++ii)
		{
			if (program_options->verbose_level>=1) {
				printf("projecting onto variable %d\n",ii);
			}
			//make the projection
			for (jj=0; jj<W.num_variables; jj++) {
				set_zero_mp(&variable_projection->coord[jj]);
			}
			set_one_mp(&variable_projection->coord[ii+1]);
			
			
			while (verify_projection_ok(W,
																	randomizer_matrix,
																	variable_projection,
																	solve_options)!=1){
				// the projection onto variable ii is bad.
				if (ii==0) {
					get_comp_rand_real_mp(&variable_projection->coord[2]);
				}
				else{
					get_comp_rand_real_mp(&variable_projection->coord[ii]);
				}
			}
			

			nullspace_config ns_config;
			compute_crit_nullspace(&Wtemp, // the returned value
														 W,            // input the original witness set
														 randomizer_matrix,
														 &variable_projection,
														 randomized_degrees,
														 1,  // dimension of ambient complex object
														 1,   //  target dimension to find
														 1,   // COdimension of the critical set to find.
														 program_options,
														 solve_options,
														 &ns_config);
			ns_config.clear();
			
			sort_for_real(&Wtemp2, Wtemp, solve_options->T); // have the real points now.
			
			
			if (Wtemp2.num_pts>0)
			{
				//compute the projections, go outside a bit, find points.
				vec_mp projection_values; init_vec_mp(projection_values,Wtemp2.num_pts); projection_values->size=Wtemp2.num_pts;
				for (jj=0; jj<Wtemp2.num_pts; jj++) {
					projection_value_homogeneous_input(&projection_values->coord[jj], Wtemp2.pts_mp[jj],variable_projection);
				}
				
				
				vec_mp projections_sorted;
				init_vec_mp(projections_sorted,Wtemp2.num_pts); projections_sorted->size = Wtemp2.num_pts;
				int *index_tracker = (int *)br_malloc(Wtemp2.num_pts*sizeof(int));
				sort_increasing_by_real(&projections_sorted, &index_tracker, projection_values);
				
				Wtemp.reset(); Wtemp2.reset();

				
				set_mp(&variable_projection->coord[0],&projections_sorted->coord[0]); // the first one
				sub_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);
				neg_mp(&variable_projection->coord[0],&variable_projection->coord[0]);
				
				lintolin_solver_main(solve_options->T.MPType,
														 W,
														 randomizer_matrix,
														 &variable_projection, 1,
														 &Wtemp,
														 solve_options);
				
				Wtemp.sort_for_real(solve_options->T);
				Wtemp.sort_for_unique(solve_options->T);
				W_crit_real->merge(Wtemp);
				
				Wtemp.reset();

				
				set_mp(&variable_projection->coord[0],&projections_sorted->coord[projections_sorted->size-1]); // the last one
				add_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);
				neg_mp(&variable_projection->coord[0],&variable_projection->coord[0]);
				
				lintolin_solver_main(solve_options->T.MPType,
														 W,
														 randomizer_matrix,
														 &variable_projection, 1,
														 &Wtemp,
														 solve_options);
				
				Wtemp.sort_for_real(solve_options->T);
				Wtemp.sort_for_unique(solve_options->T);
				W_crit_real->merge(Wtemp);
				
				Wtemp.reset();
				
				clear_vec_mp(projection_values);
				clear_vec_mp(projections_sorted);
				free(index_tracker);
			}
			
			Wtemp.reset();
			Wtemp2.reset();
		}
		
		
		if (program_options->verbose_level>=0)
			printf("completed detjac2detjac solve, have %d real solutions to work with\n\n", W_crit_real->num_pts);
		
		Wtemp2 = *W_crit_real;
//		cp_witness_set(&Wtemp2, *W_crit_real);
		W_crit_real->reset();
//		clear_witness_set(*W_crit_real); init_witness_set(W_crit_real); // prep for the merge
		sort_for_unique(W_crit_real, Wtemp2, solve_options->T);
		Wtemp2.reset();
		// if have box, intersect box with C
	} // re: bounding box calc
	
	clear_vec_mp(variable_projection);
	clear_mp(one);
	
	
	
	
	
	return 0;
}















//the linprodtodetjac method for getting the critical points
int compute_crit_linprodtodetjac(witness_set *W_crit_real, // the returned value
																 witness_set & W,
																 mat_mp randomizer_matrix,
																 vec_mp pi,
																 int num_new_linears,
																 BR_configuration *program_options,
																 solver_configuration *solve_options)
{
	
	int ii,jj,kk;
	
	int success = 1;
	
	
	//get the random complex projection, which we will use to find the crit pts
	vec_mp random_complex_projection; init_vec_mp2(random_complex_projection,W.num_variables,solve_options->T.AMP_max_prec); random_complex_projection->size = W.num_variables;
	
	set_zero_mp(&random_complex_projection->coord[0]); // first coordinate is 0
	for (ii=1; ii<W.num_variables; ii++) {
		get_comp_rand_mp(&random_complex_projection->coord[ii]); // all other coordinates random complex
	}
	
	
	//////////////
	//
	//  first of three steps for obtaining the critical points.
	//
	/////////////
	
	
	vec_mp *new_linears_full_prec = (vec_mp *) malloc(num_new_linears*sizeof(vec_mp));
	
	for (ii=0; ii<num_new_linears; ii++) {
		init_vec_mp2(new_linears_full_prec[ii],W.num_variables,solve_options->T.AMP_max_prec); new_linears_full_prec[ii]->size = W.num_variables;
		for (kk=0; kk<W.num_variables; kk++) {
			get_comp_rand_mp(&new_linears_full_prec[ii]->coord[kk]);
		}
	}
	
	
	witness_set Wtemp2;
	witness_set Wtemp;
	witness_set W_lintolin = W;
	
	//set some solver options.  these will have to be reset later;
	solve_options->allow_multiplicity = 1;
	solve_options->allow_singular = 1;
	solve_options->use_midpoint_checker = 0;
	
	
	
	lintolin_solver_main(solve_options->T.MPType,
											 W,
											 randomizer_matrix,
											 new_linears_full_prec, num_new_linears,
											 &Wtemp,
											 solve_options);
	
	
	
	
	//add the W to Wnew.
	W_lintolin.merge(Wtemp);
	Wtemp.reset();
	
	// these three write calls optional.
	W_lintolin.write_homogeneous_coordinates("lintolin_points_homogeneous");
	W_lintolin.write_dehomogenized_coordinates("lintolin_points");
	W_lintolin.write_linears("lintolin_linears");
	
	if (program_options->verbose_level>=1)
		printf("done with initial use of lintolin solver\n");
	
	
	W_lintolin.print_to_screen();
	
	
	//////////////
	//
	//  second of three steps for obtaining the critical points.
	//
	/////////////
	
	
	
	solve_options->allow_multiplicity = 1;
	solve_options->allow_singular = 1;
	
	double initial_bound_on_degree = solve_options->T.AMP_bound_on_degree;
	solve_options->T.AMP_bound_on_degree = (double) num_new_linears+1;
	linprod_to_detjac_solver_main(solve_options->T.MPType,
																W_lintolin,
																randomizer_matrix,
																random_complex_projection,
																&Wtemp,
																solve_options);
	solve_options->T.AMP_bound_on_degree = initial_bound_on_degree;
	
	
	//
	//	print_witness_set_to_screen(W_linprod);
	Wtemp.write_dehomogenized_coordinates("member_points");
	Wtemp.write_dehomogenized_coordinates("linprod_solns");
	
	
	
	//	W_linprod_good should have only the points from W_in which lie on the correct component.
	witness_set W_linprod; 
	Wtemp.incidence_number = W.incidence_number; // copy over from the input.  note that the incidence number may be different from the component number
	
	sort_for_unique(&W_linprod, Wtemp, solve_options->T);
	W_linprod.write_dehomogenized_coordinates("linprod_solns_postmembership");
	
	Wtemp.reset();
	
	
	if (program_options->verbose_level>=1)
		printf("done with linprod_to_detjac\n");
	
	
	
	//////////////
	//
	//  third of three steps for obtaining the critical points.
	//
	/////////////
	witness_set W_detjacdetjac;	
	cp_names(&W_detjacdetjac,W);
	
	solve_options->allow_multiplicity = 1;
	solve_options->allow_singular = 1;
	
	initial_bound_on_degree = solve_options->T.AMP_bound_on_degree;
	solve_options->T.AMP_bound_on_degree = (double) num_new_linears+1;
	
	detjac_to_detjac_solver_main(solve_options->T.MPType,
															 W_linprod,
															 randomizer_matrix,
															 random_complex_projection, // move from
															 pi,  // move to
															 &W_detjacdetjac,
															 solve_options);
	
	solve_options->T.AMP_bound_on_degree = initial_bound_on_degree;
	W_detjacdetjac.write_dehomogenized_coordinates("detjac_solns");
	
	
	sort_for_real(W_crit_real, W_detjacdetjac,solve_options->T); // get only the real solutions.
	
	
	
	
	
	// NOW WE HAVE THE CRITICAL POINTS
	
	
	
	
	
	
	
	/////now compute the bounds, whether bounding box or bounding the interval
	
	comp_mp one; init_mp(one);
	set_one_mp(one);
	
	vec_mp variable_projection; init_vec_mp(variable_projection, W.num_variables); variable_projection->size = W.num_variables;
	
	if (!program_options->use_bounding_box) {// do not want bounding box
		
		if (program_options->verbose_level>=0)
			printf("in no-box mode.  computing bounds of interval\n");
		
		
		vec_cp_mp(variable_projection,pi);
		
		// compute projections of W_crit_real
		vec_mp projection_values; init_vec_mp(projection_values,W_crit_real->num_pts); projection_values->size=W_crit_real->num_pts;
		for (jj=0; jj<W_crit_real->num_pts; jj++) {
			projection_value_homogeneous_input(&projection_values->coord[jj], W_crit_real->pts_mp[jj],pi);
		}
		
		vec_mp projections_sorted;
		init_vec_mp(projections_sorted,W_crit_real->num_pts); projections_sorted->size = W_crit_real->num_pts;
		int *index_tracker = (int *)br_malloc(W_crit_real->num_pts*sizeof(int));
		sort_increasing_by_real(&projections_sorted, &index_tracker, projection_values);
		
		printf("%d points in W_crit_real\n", W_crit_real->num_pts);
		print_point_to_screen_matlab(projections_sorted,"proj_sorted");
		
		
		comp_mp left_proj_val, right_proj_val;
		init_mp(left_proj_val); init_mp(right_proj_val);
		
		if (W_crit_real->num_pts==0) {
			set_one_mp(right_proj_val);
			set_one_mp(left_proj_val);
			neg_mp(left_proj_val,left_proj_val);
		}
		else{
			set_mp(left_proj_val, &projections_sorted->coord[0]);
			set_mp(right_proj_val, &projections_sorted->coord[W_crit_real->num_pts-1]);
		}
		
		
		comp_mp relevant_proj_val; init_mp(relevant_proj_val);
		//this loop presupposes that w_crit_real has at least one member
		for (jj=0; jj<2; jj++) {
			
			
			//set the value for the projection to move to
			if (jj==0){ // if on the first of the two iterations
				set_mp(relevant_proj_val,left_proj_val);}
			else{
				set_mp(relevant_proj_val,right_proj_val);}
			
			set_mp(&variable_projection->coord[0],relevant_proj_val);
			
			if (jj==0){
				sub_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);}
			else{
				add_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);}
			
			neg_mp(&variable_projection->coord[0],&variable_projection->coord[0]);
			//done setting the value of the projection.
			
			//run the solver
			lintolin_solver_main(solve_options->T.MPType,
													 W,
													 randomizer_matrix,
													 &variable_projection, 1,
													 &Wtemp,
													 solve_options);
			
			Wtemp.sort_for_real(solve_options->T);
			Wtemp.sort_for_unique(solve_options->T);
			W_crit_real->merge(Wtemp);
			
			Wtemp.reset();
		}
		
		clear_mp(relevant_proj_val);
	}
	else //want a bounding box
	{
		
		if (program_options->verbose_level>=0) {
			printf("computing bounding box\n");
		}
		// for each variable direction, compute projections of crit points onto that axis.
		
		
		for (ii=0; ii<W.num_variables-1; ++ii)
		{
			if (program_options->verbose_level>=1) {
				printf("projecting onto variable %d\n",ii);
			}
			//make the projection
			for (jj=0; jj<W.num_variables; jj++) {
				set_zero_mp(&variable_projection->coord[jj]);
			}
			set_one_mp(&variable_projection->coord[ii+1]);
			
			
			while (verify_projection_ok(W,
																	randomizer_matrix,
																	variable_projection,
																	solve_options)!=1){
				// the projection onto variable ii is bad.
				if (ii==0) {
					get_comp_rand_real_mp(&variable_projection->coord[2]);
				}
				else{
					get_comp_rand_real_mp(&variable_projection->coord[ii]);
				}
			}
			
			
			initial_bound_on_degree = solve_options->T.AMP_bound_on_degree;
			solve_options->T.AMP_bound_on_degree = (double) num_new_linears+1;
			detjac_to_detjac_solver_main(solve_options->T.MPType,
																	 W_linprod,
																	 randomizer_matrix,
																	 random_complex_projection, // move from
																	 variable_projection,  // move to
																	 &Wtemp,
																	 solve_options);
			solve_options->T.AMP_bound_on_degree = initial_bound_on_degree;
			
			sort_for_real(&Wtemp2, Wtemp, solve_options->T); // have the real points now.
			
			
			if (Wtemp2.num_pts>0)
			{
				//compute the projections, go outside a bit, find points.
				vec_mp projection_values; init_vec_mp(projection_values,Wtemp2.num_pts); projection_values->size=Wtemp2.num_pts;
				for (jj=0; jj<Wtemp2.num_pts; jj++) {
					projection_value_homogeneous_input(&projection_values->coord[jj], Wtemp2.pts_mp[jj],variable_projection);
				}
				
				
				vec_mp projections_sorted;
				init_vec_mp(projections_sorted,Wtemp2.num_pts); projections_sorted->size = Wtemp2.num_pts;
				int *index_tracker = (int *)br_malloc(Wtemp2.num_pts*sizeof(int));
				sort_increasing_by_real(&projections_sorted, &index_tracker, projection_values);
				
				Wtemp.reset(); Wtemp2.reset();

				
				set_mp(&variable_projection->coord[0],&projections_sorted->coord[0]); // the first one
				sub_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);
				neg_mp(&variable_projection->coord[0],&variable_projection->coord[0]);
				
				lintolin_solver_main(solve_options->T.MPType,
														 W,
														 randomizer_matrix,
														 &variable_projection, 1,
														 &Wtemp,
														 solve_options);
				
				Wtemp.sort_for_real(solve_options->T);
				Wtemp.sort_for_unique(solve_options->T);
				W_crit_real->merge(Wtemp);
				
				Wtemp.reset();

				
				set_mp(&variable_projection->coord[0],&projections_sorted->coord[projections_sorted->size-1]); // the last one
				add_mp(&variable_projection->coord[0],&variable_projection->coord[0],one);
				neg_mp(&variable_projection->coord[0],&variable_projection->coord[0]);
				
				lintolin_solver_main(solve_options->T.MPType,
														 W,
														 randomizer_matrix,
														 &variable_projection, 1,
														 &Wtemp,
														 solve_options);
				
				Wtemp.sort_for_real(solve_options->T);
				Wtemp.sort_for_unique(solve_options->T);
				W_crit_real->merge(Wtemp);
				
				Wtemp.reset();
				
				clear_vec_mp(projection_values);
				clear_vec_mp(projections_sorted);
				free(index_tracker);
			}
			
			Wtemp.reset();
			Wtemp2.reset();
		}
		
		
		if (program_options->verbose_level>=0)
			printf("completed detjac2detjac solve, have %d real solutions to work with\n\n", W_crit_real->num_pts);
		
		Wtemp2 = *W_crit_real;
//		cp_witness_set(&Wtemp2, *W_crit_real);
		W_crit_real->reset();
		sort_for_unique(W_crit_real, Wtemp2, solve_options->T);
		Wtemp2.reset();
		// if have box, intersect box with C
	} // re: bounding box calc
	
	clear_vec_mp(variable_projection);
	clear_mp(one);
	
	W_linprod.reset();
	W_lintolin.reset();
	W_detjacdetjac.reset();

	return success;
}




int get_sum_degrees(char filename[], int num_funcs){
	int degsum = 0, tmpdeg, ii;
	
	FILE *IN;
	
	IN =  safe_fopen_read(filename);
	
	for (ii = 0; ii<num_funcs; ii++) {
		fscanf(IN,"%d",&tmpdeg);
		degsum += tmpdeg;
	}
	
	fclose(IN);
	
	
	return degsum;
}


//
//
//void sort_for_membership(char * input_file,
//												 witness_set *W_out,
//												 witness_set W_in,
//												 char *stifle_text){
////	printf("sorting points for membership\n");
//
//
//	if (W_in.incidence_number==-1) {
//		printf("input witness_set has unset incidence_number for comparison.\n");
//		exit(-1112);
//	}
//
//	W_out->MPType = W_in.MPType;
//	W_out->incidence_number = W_in.incidence_number;
//	W_out->num_variables = W_in.num_variables;
//	W_out->num_var_gps = W_in.num_var_gps;
//
//	//copy the linears, patches from old to new
//	cp_patches(W_out, W_in);
//	cp_linears(W_out, W_in);
//	cp_names(W_out, W_in);
//
//	//more important files out of the way.  this will probably crash the program if it is called without these files being present.
//	rename_bertini_files_dotbak();
//
//
//	int strLength = 0, digits = 15, *declarations = NULL;
//  char *SysStr = NULL,*fmt = NULL, *bertini_command="bertini";
//  FILE *IN = NULL;
//
//
//
//	//make the command string to run
//	strLength = 1 + snprintf(NULL, 0, "%s input_membership_test %s ", bertini_command, stifle_text);
//  SysStr = (char *)br_malloc(strLength * sizeof(char));
//  sprintf(SysStr, "%s input_membership_test %s ", bertini_command, stifle_text);
//
//
//
//	// make the format to write the member_points file
//  strLength = 1 + snprintf(NULL, 0, "%%.%dle %%.%dle\n", digits, digits);
//  // allocate size
//  fmt = (char *)br_malloc(strLength * sizeof(char));
//  // setup fmt
//  sprintf(fmt, "%%.%dle %%.%dle\n", digits, digits);
//
//
//  // setup input file
//	IN = safe_fopen_read(input_file);
//  partitionParse(&declarations, IN, "func_input_real", "config_real",0); // the 0 means not self conjugate
//	fclose(IN);
//
//
//	//check existence of the required witness_data file.
//	IN = safe_fopen_read("witness_data");
//	fclose(IN);
//
//
//
//
//
//
//	//only need to do this once.  we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
//  membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
//
//
//
//
//	// Do membership test
//	printf("*\n%s\n*\n",SysStr);
//  system(SysStr);
//
//
//	int on_component_indicator[W_in.num_pts];
//	read_incidence_matrix_wrt_number(on_component_indicator,W_in.incidence_number);
//
////	printf("done reading incidence_matrix\n");
//
//	int *is_on_component;
//	is_on_component = (int *)br_malloc(W_in.num_pts*sizeof(int));
//
//
//
//	int num_pts_on_component =0;
//	int ii,jj;
//	for (ii=0; ii<W_in.num_pts; ii++) {
////		printf("%d\n",on_component_indicator[ii]);
//
//		if (on_component_indicator[ii]==1) {
//			num_pts_on_component++;
//			is_on_component[ii]=1;
//		}
//		else{
//			is_on_component[ii]=0;
//		}
//
//	}
//
//	if (num_pts_on_component==0) {
//		printf("found 0 points out of %d candidates, on component.\n",W_in.num_pts);
//		restore_bertini_files_dotbak();
//		return;
//	}
//
//
//	vec_d *points_on_component;
//	points_on_component =(point_d *)br_malloc(num_pts_on_component*sizeof(point_d));
//
//	int *indices_on_component;
//	indices_on_component = (int *)br_malloc(num_pts_on_component*sizeof(int));
//
//	int counter =0;
//	for (ii=0; ii<W_in.num_pts; ii++) {
//		if (on_component_indicator[ii]==1) {
//			init_vec_d(points_on_component[counter],W_in.pts_d[ii]->size);
//			points_on_component[counter]->size = W_in.pts_d[ii]->size;
//			vec_cp_d(points_on_component[counter],W_in.pts_d[ii]);
//
//			indices_on_component[counter] = ii;
//
//			counter++;
//		}
//	}
//
//
//
//	//now remove duplicates
//
//
//
//	int *is_unique; int curr_uniqueness; int num_good_pts = 0;
//	is_unique  = (int *)br_malloc(num_pts_on_component*sizeof(int));
//	for (ii = 0; ii<num_pts_on_component; ++ii) {
//		curr_uniqueness = 1;
//		if (ii!= (num_pts_on_component-1) ) { // the last point is unique...  always
//			for (jj=ii+1; jj<num_pts_on_component; ++jj) {
//				if (isSamePoint_homogeneous_input_d(points_on_component[ii],points_on_component[jj])){
//					//formerly (isSamePoint(points_on_component[ii],NULL,52,points_on_component[jj],NULL,52,1e-8)){
//					curr_uniqueness = 0;
////					printf("the following two points are not distinct:\n");
////					print_point_to_screen_matlab(points_on_component[ii],"left");
////					print_point_to_screen_matlab(points_on_component[jj],"right");
//				}
//			}
//		}
//
//
//		if (curr_uniqueness==1) {
//			is_unique[ii] = 1;
//			num_good_pts++;
//		}
//		else
//		{
//			is_unique[ii] = 0;
//		}
//
//	}
//
//
//
//
//
//
//
//
//	//allocate the memory
//	W_out->pts_d=(point_d *)br_malloc(num_good_pts*sizeof(point_d));
//	W_out->pts_mp=(point_mp *)br_malloc(num_good_pts*sizeof(point_mp));
//	W_out->num_pts = num_good_pts;
//	W_out->num_pts = num_good_pts;
//
//
//
//	//copy in the distinct points lying on the correct component
//	counter = 0;
//	for (ii=0; ii<num_pts_on_component; ++ii) {
//		if (is_unique[ii]==1) {
//			init_vec_d(W_out->pts_d[counter],W_in.num_variables); W_out->pts_d[counter]->size = W_in.num_variables;
//			init_vec_mp2(W_out->pts_mp[counter],W_in.num_variables,1024);  W_out->pts_mp[counter]->size = W_in.num_variables;
//
//			vec_cp_d(W_out->pts_d[counter], W_in.pts_d[indices_on_component[ii]]);
//			vec_cp_mp(W_out->pts_mp[counter], W_in.pts_mp[indices_on_component[ii]]);
//			counter++;
//		}
//	}
//
//	if (!counter==num_good_pts){
//		printf("counter mismatch; counter!=num_good_pts in sort_for_membership\n");
//		exit(169);
//	}
//
//		//clear the memory
//
//	for (ii=0; ii<num_pts_on_component; ++ii) {
//		clear_vec_d(points_on_component[ii]);
//	}
//	free(points_on_component);
//
//
//	free(is_unique);
//  free(declarations);
//
//  // delete temporary files
//  remove("func_input_real");
//  remove("config_real");
//	remove("incidence_matrix");
////	remove("member_points");
//
//	//move files back into place
//	restore_bertini_files_dotbak();
//	return;
//}
//
//




void sort_increasing_by_real(vec_mp *projections_sorted, int **index_tracker, vec_mp projections_input){
	
	comp_mp large; init_mp(large);
	comp_d l; l->r = 1e9; l->i = 0;
	d_to_mp(large,l);
	
	
	vec_mp raw; init_vec_mp(raw,1);
	vec_cp_mp(raw,projections_input);
	
	change_size_vec_mp( (*projections_sorted), raw->size);
	(*projections_sorted)->size = raw->size;
	int ii,jj;
	//	comp_mp temp1, temp2; init_mp(temp1); init_mp(temp2);
	
	double min;
	double curr;
	int indicator = -1;
	for (ii=0; ii<raw->size; ii++) {
		min = 1e10;
		
		for (jj=0; jj<raw->size; jj++) {
			curr = mpf_get_d(raw->coord[jj].r);
			if ( curr < min) {
				indicator = jj;
				min = curr;
			}
		}
		if (indicator==-1) {
			printf("min projection value was *insanely* large\n");
			exit(1111);
		}
		
		
		
		(*index_tracker)[ii] = indicator;
		set_mp( &(*projections_sorted)->coord[ii],&raw->coord[indicator]);
		set_mp( &raw->coord[indicator],large);
	}
	return;
}



//input the raw number of variables including the homogeneous variable (of which there must be one)
// assume the array of integers 'randomized_degrees' is already initialized to the correct size.
void make_randomization_matrix_based_on_degrees(mat_mp randomization_matrix, int ** randomized_degrees, int num_desired_rows, int num_funcs)
{
	int ii,jj;
	
	
	
	//get unique degrees
	int *degrees = (int *) br_malloc(num_funcs*sizeof(int));
	int *unique_degrees = (int *) br_malloc(num_funcs*sizeof(int));
	
	
	FILE *IN = safe_fopen_read("deg.out"); //open the deg.out file for reading.
	int num_unique_degrees = 0;
	int occurrence_counter;
	for (ii=0; ii<num_funcs; ++ii) {
		fscanf(IN,"%d\n",&degrees[ii]); // read data
		occurrence_counter = 0; // set the counter for how many timmes the current degree has already been found.
		for (jj=0; jj<ii; jj++) {
			if (degrees[jj]==degrees[ii]) { // if previously stored degree is same as current one
				occurrence_counter++; // increment counter
			}
		}
		
		if (occurrence_counter==0) { // if did not find already in list
			unique_degrees[num_unique_degrees] = degrees[ii]; // add to list of unique degrees.
			num_unique_degrees++; // have one more unique degree
		} // re: jj
	}// re: ii
	fclose(IN);
	
	
	if (num_desired_rows==num_funcs) {
		make_matrix_ID_mp(randomization_matrix,num_funcs,num_funcs);
		for (ii=0; ii<num_desired_rows; ++ii) {
			(*randomized_degrees)[ii] = degrees[ii];
		}
		free(degrees);
		free(unique_degrees);
		return;
	}
	
	//sort the unique degrees into decreasing order
	qsort(unique_degrees, num_unique_degrees, sizeof(int), compare_integers_decreasing);
	
	//count how many of each unique degree there are.
	int *num_of_each_degree = (int *) br_malloc(num_unique_degrees*sizeof(int));
	for (ii=0; ii<num_unique_degrees; ii++) {
		num_of_each_degree[ii] = 0;
		for (jj=0; jj<num_funcs; ++jj) {
			if (unique_degrees[ii]==degrees[jj]) {
				num_of_each_degree[ii]++;
			}
		}
	}
	
	
	
	//	for (ii=0; ii<num_unique_degrees; ii++) {
	//		printf("unique_degrees[%d]=%d; num_of_each_degree=%d\n",ii,unique_degrees[ii],num_of_each_degree[ii]);
	//	}
	
	//resize the matrix
	change_size_mat_mp(randomization_matrix,num_desired_rows,num_funcs);
	randomization_matrix->rows = num_desired_rows; randomization_matrix->cols = num_funcs;
	
	
	
	int counter = 0;
	int current_degree_index = 0; // start at the end
	int current_degree;
	for (ii=0; ii<num_desired_rows; ii++) {
		
		counter++;
		if (counter>num_of_each_degree[current_degree_index]) {
			current_degree_index++;
			counter = 1;
		}
		
		current_degree = unique_degrees[current_degree_index];
		(*randomized_degrees)[ii] = current_degree;
		
		int encountered_current_degree = 0;
		for (jj=0; jj<num_funcs; jj++) {
			if ( (degrees[jj]<= current_degree)  ) {
				encountered_current_degree++;
				if (encountered_current_degree >= counter){
					get_comp_rand_real_mp(&randomization_matrix->entry[ii][jj]);
				}
				else{
					set_zero_mp(&randomization_matrix->entry[ii][jj]);
				}
			}
			else
			{
				set_zero_mp(&randomization_matrix->entry[ii][jj]);
			}
		}
		
		
	}
	
	free(num_of_each_degree);
	free(degrees);
	free(unique_degrees);
	
	return;
}


int compare_integers_decreasing(const void * left_in, const void * right_in){
	
	int left = *(const int*)left_in;
	int right = *(const int*)right_in;
	
	
	if (left<right) {
		return 1;
	}
	else if(right > left){
		return -1;
	}
	else{
		return 0;
	}
	
}






int verify_projection_ok(witness_set & W,
												 vec_mp projection,
												 solver_configuration *solve_options)
{
	

	//create a matrix
	mat_mp randomizer_matrix;
	init_mat_mp2(randomizer_matrix,W.num_variables-W.num_patches-W.dim,solve_options->PPD.num_funcs,solve_options->T.AMP_max_prec);
	
	//create the array of integers
	int *randomized_degrees = (int *)br_malloc((W.num_variables-W.num_patches-W.dim)*sizeof(int));
	
	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(randomizer_matrix, &randomized_degrees, W.num_variables-W.num_patches-W.dim, solve_options->PPD.num_funcs);
	
	int invalid_flag = verify_projection_ok(W, randomizer_matrix, projection, solve_options);
	
	clear_mat_mp(randomizer_matrix);
	free(randomized_degrees);
	
	return invalid_flag;
}



int verify_projection_ok(witness_set & W,
												 mat_mp randomizer_matrix,
												 vec_mp projection,
												 solver_configuration *solve_options)
{
	int ii,jj;
	
	
	int invalid_flag;
	
	
	vec_mp temp_rand_point;  init_vec_mp(temp_rand_point,W.num_variables); temp_rand_point->size = W.num_variables;
	set_one_mp(&temp_rand_point->coord[0]); // first coordinate must be 1
	for (ii=1; ii<W.num_variables; ++ii) {
		get_comp_rand_mp(&temp_rand_point->coord[ii]);
	}
	
	prog_t SLP;
	setupProg(&SLP, solve_options->T.Precision, 2);
	
	
	comp_mp zerotime; init_mp(zerotime);
	set_zero_mp(zerotime);
	
	
	std::cout << "in verify_projection, SLP has " << SLP.numVars << " variables\n";
	std::cout << solve_options->T.MPType << " is MPType\n";
	
	print_point_to_screen_matlab(temp_rand_point,"temp_rand_point");
	print_comp_mp_matlab(zerotime,"zerotime");
	
	
	
	eval_struct_mp ED; init_eval_struct_mp(ED, 0, 0, 0);
	evalProg_mp(ED.funcVals, ED.parVals, ED.parDer, ED.Jv, ED.Jp, temp_rand_point, zerotime, &SLP);
	
	mat_mp AtimesJ; init_mat_mp(AtimesJ, 1, 1); AtimesJ->rows = AtimesJ->cols = 1;
	
	mat_mul_mp(AtimesJ, randomizer_matrix, ED.Jv);
	
	//inflate the matrix
	
	mat_mp detme;  init_mat_mp(detme, W.num_variables-1, W.num_variables-1);
	detme->cols = detme->rows= W.num_variables-1;
	for (ii=0; ii<W.num_variables-2; ++ii) {
		for (jj=0; jj<W.num_variables-1; ++jj) {
			set_mp(&detme->entry[ii][jj],&AtimesJ->entry[ii][jj+1]); // omit the homogeneous coordinate's columns
		}
	}
	
	for (ii=0; ii<W.num_variables-1; ii++) {
		set_mp(&detme->entry[W.num_variables-2][ii],&projection->coord[ii+1]);
	}
	
	comp_mp determinant; init_mp(determinant);
	take_determinant_mp(determinant,detme);
	
	if ( d_abs_mp(determinant) < 1e-12){
		invalid_flag = 0;
		std::cout << d_abs_mp(determinant) << "\n";
		print_matrix_to_screen_matlab(ED.Jv,"Jv");
		print_matrix_to_screen_matlab(detme,"detme");
	}
	else
		invalid_flag = 1;
	
	return invalid_flag;
	
}


