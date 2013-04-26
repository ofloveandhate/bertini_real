#include "curveSelfConj.h"










void computeCurveSelfConj(char * inputFile,
													witness_set W,
													vec_mp pi,
													curveDecomp_d *C,
													int num_vars,
													int num_var_gps,
													program_configuration *program_options,
													solver_configuration *solve_options)
{
	//IN DEVELOPMENT
  

	int ii,kk,jj;
	
	witness_set Wtemp,Wtemp2;

	FILE *IN = safe_fopen_read(inputFile);
	int *declarations = NULL;
	partitionParse(&declarations, IN, "func_input", "config",0); // the 0 means not self conjugate.
	
	preproc_data PPD;
	setupPreProcData("preproc_data", &PPD);
	
	

		
	
	//get the random complex projection, which we will use to find the crit pts
	vec_mp random_complex_projection; init_vec_mp2(random_complex_projection,W.num_variables,solve_options->T.AMP_max_prec); random_complex_projection->size = W.num_variables;
	
	set_zero_mp(&random_complex_projection->coord[0]); // first coordinate is 0
	for (ii=1; ii<W.num_variables; ii++) {
		get_comp_rand_mp(&random_complex_projection->coord[ii]); // all other coordinates random complex
	}
	
	
	
	
	
	
	
	
	
	
	
	
	// 2) randomize down to N-1 equations
	// to get a square system for the homotopies in the following steps.
	
	//actually, just assign a random real nearly orthogonal matrix.
	
	//create the matrix
	mat_mp n_minusone_randomizer_matrix_full_prec;
	init_mat_mp2(n_minusone_randomizer_matrix_full_prec,W.num_variables-2,PPD.num_funcs,solve_options->T.AMP_max_prec);
	
	//create the array of integers
	int *randomized_degrees = (int *)bmalloc((W.num_variables-2)*sizeof(int));

	//get the matrix and the degrees of the resulting randomized functions.
	make_randomization_matrix_based_on_degrees(n_minusone_randomizer_matrix_full_prec, &randomized_degrees, W.num_variables-2, PPD.num_funcs);
	


	
	//compute the number of linears we will need to move to, to perform the regeneration for linprodtodetjac method.
	
	int deg_g = 0;
	for (ii=0; ii<W.num_variables-2; ii++) {
		deg_g+= randomized_degrees[ii]-1;
	}
	deg_g -= 1;
#ifdef verbose
	printf("the number of new linears will be %d\n",deg_g);
#endif
	
	
	
	
	vec_mp *new_linears_full_prec; new_linears_full_prec = (vec_mp *) malloc(deg_g*sizeof(vec_mp));
	
	for (ii=0; ii<deg_g; ii++) {
		init_vec_mp2(new_linears_full_prec[ii],W.num_variables,solve_options->T.AMP_max_prec); new_linears_full_prec[ii]->size = W.num_variables;
		for (kk=0; kk<W.num_variables; kk++) {
			get_comp_rand_mp(&new_linears_full_prec[ii]->coord[kk]);
		}
	}
	
	
	
	FILE *OUT=NULL;
	OUT = safe_fopen_write("Rand_matrix");
	fprintf(OUT,"%d %d\n",n_minusone_randomizer_matrix_full_prec->rows,n_minusone_randomizer_matrix_full_prec->cols);
	fprintf(OUT,"%d\n",solve_options->T.AMP_max_prec);
	fprintf(OUT,"\n");
	print_matrix_to_file_mp(OUT, 0, n_minusone_randomizer_matrix_full_prec);
	fclose(OUT);
	
	
	
	// 1) check rank
	// a) use jacobian
	//TODO: RANK CHECK WITH PROJECTION
	
	
	
	// b) use INCREASE_MATRIX from bertini to make more stuff on the bottom
	// c) put 'pi' underneath
	// d) check rank
	
	// a little test code:
	//	increase_size_mat_d(Jv, Jv->rows+1, Jv->cols);
	//	Jv->rows = Jv->rows+1;
	//	for (ii = 0; ii<Jv->cols; ii++) {
	//		Jv->entry[Jv->rows-1][ii].r = random_real_linear->coord[ii].r;
	//		Jv->entry[Jv->rows-1][ii].i = random_real_linear->coord[ii].i;
	//	}
	//	print_matrix_to_screen_matlab(Jv, "jacobian_with_random_real_linear");
	
	
	
	//	// now need to get determinant of jacobian.
	
	
  

	
	
	
	
  // 3) initalize data structures
	//  this initialization is performed before -this- function?
  
	
	
  // 4) solve for critical conditions for random complex projection
	// if no box, determine box
	//  int whatever = compute_box();
	
	// for each variable direction, compute projections of crit points onto that axis.
	
	
	//////////////
	//
	//  first of three steps for obtaining the critical points.
	//
	/////////////
	
	

	init_witness_set_d(&Wtemp);
	cp_patches(&Wtemp,W); // copy the patches over from the original witness set
	
	witness_set W_lintolin;
	init_witness_set_d(&W_lintolin);
	cp_patches(&W_lintolin,W); // copy the patches over from the original witness set
	
	// if have box, intersect box with C
	
	lin_to_lin_solver_main(solve_options->T.MPType,
												 W,
												 n_minusone_randomizer_matrix_full_prec,
												 new_linears_full_prec, deg_g,
												 &Wtemp,
												 solve_options);
	
	printf("made it out of the lin_to_lin solver\n");
	

	//add the W to Wnew.
	merge_witness_sets(&W_lintolin,Wtemp,W);

	clear_witness_set(Wtemp);
	
//	check_patch_values(W_lintolin);
	write_homogeneous_coordinates(W_lintolin, "lintolin_points_homogeneous");
	write_dehomogenized_coordinates(W_lintolin,"lintolin_points");
	write_linears(W_lintolin,"lintolin_linears");
	
	
	
	
	
	
	
	
	//////////////
	//
	//  second of three steps for obtaining the critical points.
	//
	/////////////
	
	witness_set W_linprod; init_witness_set_d(&W_linprod);
	cp_names(&W_lintolin,W);
	
	
	linprod_to_detjac_solver_main(solve_options->T.MPType,
																W_lintolin,
																n_minusone_randomizer_matrix_full_prec,
																random_complex_projection,
																&W_linprod,
																solve_options);
	
	
	
	//
	//	print_witness_set_to_screen(W_linprod);
	write_dehomogenized_coordinates(W_linprod,"member_points");
	write_dehomogenized_coordinates(W_linprod,"linprod_solns");
	witness_set W_linprod_good;
	init_witness_set_d(&W_linprod_good);
	W_linprod.incidence_number = W.incidence_number; // copy over from the input.  note that the incidence number may be different from the component number
	
	
	//	W_linprod_good should have only the points from W_in which lie on the correct component.
	sort_for_membership(inputFile, &W_linprod_good, W_linprod, program_options->stifle_text);
	write_dehomogenized_coordinates(W_linprod_good,"linprod_solns_postmembership");
	
	printf("done sorting membership\n");
	
//	print_witness_set_to_screen(W_linprod_good);
	
	
	

	printf("done with linprod_to_detjac\n");	
	
	
	if (W_linprod_good.num_pts==0) {
		printf("found 0 solutions to feed into the detjactodetjac solver.  exiting.\n");
		exit(-11);
	}

	
	
	//////////////
	//
	//  third of three steps for obtaining the critical points.
	//
	/////////////
	witness_set W_detjacdetjac; init_witness_set_d(&W_detjacdetjac);
	cp_names(&W_detjacdetjac,W);
	
	
	detjac_to_detjac_solver_main(solve_options->T.MPType,
															 W_linprod_good,
															 n_minusone_randomizer_matrix_full_prec,
															 random_complex_projection, // move from
															 pi,  // move to
															 &W_detjacdetjac,
															 solve_options);
	
//	print_witness_set_to_screen(W_detjacdetjac);
	write_dehomogenized_coordinates(W_detjacdetjac,"detjac_solns");

	

	witness_set W_crit_real; init_witness_set_d(&W_crit_real);
	sort_for_real(&W_crit_real, W_detjacdetjac,solve_options->T);
	
	


	printf("completed detjac2detjac solve, have %d real solutions to work with\n\n", W_crit_real.num_pts);

	

	
	
  // 5) compute h_o, get critical points
	
	
	printf("***\n computing crit\n***\n");
	
	solve_options->allow_multiplicity=1;
	solve_options->allow_singular=1;
	
	vec_mp projection_values; init_vec_mp(projection_values,W_crit_real.num_pts); projection_values->size = W_crit_real.num_pts;
	vec_mp dehom; 
	
	vec_mp pi_no_hom; init_vec_mp(pi_no_hom,W_crit_real.num_variables-1); pi_no_hom->size = W_crit_real.num_variables-1;
	for (ii=0; ii<W_crit_real.num_variables-1; ii++) {
		set_mp(&pi_no_hom->coord[ii],&pi->coord[ii+1]);  // omit the first coordinate
	}
	
	
	for (ii=0; ii<W_crit_real.num_pts; ii++) {
		projection_value_homogeneous_input(&projection_values->coord[ii],W_crit_real.pts_mp[ii], pi); // set projection value
	}

	
	int *index_tracker = (int *)bmalloc(W_crit_real.num_pts*sizeof(vec_mp)*sizeof(int));
	
	vec_mp projections_sorted;
	init_vec_mp(projections_sorted,W_crit_real.num_pts); projections_sorted->size = W_crit_real.num_pts;
	
	
	printf("sorting projection values\n");
	sort_increasing_by_real(&projections_sorted, &index_tracker, projection_values);
	
	
//	vec_mp *projections_sorted, vec_mp **points_sorted, int **index_tracker, vec_mp projections_input
//	if (W_crit_real.num_pts>0) {
//		print_point_to_screen_matlab_mp(projection_values,"projection_values");
//		print_point_to_screen_matlab_mp(projections_sorted,"sorted_projection_values");
//	}
	// projections_sorted
	
	
	
	
	
	////////
	//
	//  go right and left to get endpoints for interval of interest
	//
	///////////

	
	double left;
	double right;
	
	if (W_crit_real.num_pts>0) {
		left = mpf_get_d(projections_sorted->coord[0].r)-1;
		right = mpf_get_d(projections_sorted->coord[W_crit_real.num_pts-1].r)+1;
	}
	else{
		left = -1;
		right = 1;
	}
//TODO: make these '1' better;

	
	
	witness_set left_endpoint; init_witness_set_d(&left_endpoint);
	witness_set right_endpoint; init_witness_set_d(&right_endpoint);
	
	comp_d temp;
  vec_mp particular_projection;  init_vec_mp(particular_projection,W.num_variables); particular_projection->size = W.num_variables;
	vec_cp_mp(particular_projection,pi);

	
	/////
	//
	// go left
	//
	////////
	init_witness_set_d(&Wtemp);init_witness_set_d(&Wtemp2);
	temp->r = -left; temp->i = 0;
	d_to_mp(&particular_projection->coord[0],temp);
	printf("going left\n");
	lin_to_lin_solver_main(solve_options->T.MPType,
												 W,
												 n_minusone_randomizer_matrix_full_prec,
												 &particular_projection, 1,
												 &Wtemp,
												 solve_options);
	
	sort_for_real(&Wtemp2,Wtemp,solve_options->T);
	sort_for_unique(&left_endpoint, Wtemp2, solve_options->T);
	int left_has_real_soln;
	if (left_endpoint.num_pts>0){
		left_has_real_soln = 1;
//		merge_witness_sets(mcreallerson);
	}
	else{
		left_has_real_soln = 0;
	}
	clear_witness_set(Wtemp);clear_witness_set(Wtemp2);
	
	// end going left
	
	
	/////
	//
	// go right
	//
	////////
	
	//set the value of the projection to which we wish to homotope.
	temp->r = -right; temp->i = 0;
	d_to_mp(&particular_projection->coord[0],temp);
	init_witness_set_d(&Wtemp);init_witness_set_d(&Wtemp2);
	printf("going right\n");
	lin_to_lin_solver_main(solve_options->T.MPType,
												 W,
												 n_minusone_randomizer_matrix_full_prec,
												 &particular_projection, 1,
												 &Wtemp,
												 solve_options);
	
	sort_for_real(&Wtemp2, Wtemp, solve_options->T);
	sort_for_unique(&right_endpoint, Wtemp2, solve_options->T);
	int right_has_real_soln;
	if (right_endpoint.num_pts>0){
		right_has_real_soln = 1;
	}
	else{
		right_has_real_soln = 0;
	}
	
	clear_witness_set(Wtemp);clear_witness_set(Wtemp2);
	
	//end going right

	
	
	
	
	///////
	//
	//   actually form crit.
	//
	/////////
	

	/// dehomogenize the points to perform the projections using \pi.  the points we
	vec_d mcdehom;  init_vec_d(mcdehom,W.num_variables-1); mcdehom->size = W.num_variables-1;

	
	double *crit_downstairs; int num_crit;
	
	num_crit = W_crit_real.num_pts + left_has_real_soln + right_has_real_soln;
	crit_downstairs = (double *)bmalloc(num_crit*sizeof(double));
	
	int num_midpoints = W_crit_real.num_pts-1 + left_has_real_soln + right_has_real_soln;
	double *midpoints_downstairs = (double *) bmalloc(num_midpoints*sizeof(double));
	
	if (num_midpoints<1) {
		printf("no midpoints to work with :(\n");
		printf("please program exiting setting C\n");
		exit(-1);
	}
	else{
		printf("%d midpoints\n",num_midpoints);
	}


	if (left_has_real_soln==1) {
		if (W_crit_real.num_pts==0) {
			midpoints_downstairs[0] = (left+right)/2;
		}
		else{
			midpoints_downstairs[0] = (left+mpf_get_d(projections_sorted->coord[0].r))/2; //halfway between the endpoint and the leftmost critical point.
		}
		crit_downstairs[0] = left;
	}

	
	int offset = left_has_real_soln;
	for (ii=0; ii<projections_sorted->size-1; ii++) {
		midpoints_downstairs[ii+offset] = (mpf_get_d(projections_sorted->coord[ii].r)+mpf_get_d(projections_sorted->coord[ii+1].r))/2;
	}
	for (ii=0; ii<W_crit_real.num_pts; ii++) {
		crit_downstairs[ii+offset] = mpf_get_d(projections_sorted->coord[ii].r);
	}

	if (right_has_real_soln==1) {
		
		if (W_crit_real.num_pts>0) {
			midpoints_downstairs[num_midpoints-1] = (right+ mpf_get_d(projections_sorted->coord[W_crit_real.num_pts-1].r) )/2; // halfway between the rightmost crit point found earlier, and the right enpoint (which has at least one solution).
		}
		crit_downstairs[num_crit-1] = right;
	}
	
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

	
	
	
	
	
	
	
	//copy the points from witness sets into curve decomposition V0;
	
	
	vertex temp_vertex;  init_vertex_mp(&temp_vertex);
	
	for (ii=0; ii<left_endpoint.num_pts; ii++) {
		printf("adding left endpoint %d of %d\n",ii, right_endpoint.num_pts);
		
		projection_value_homogeneous_input(temp_vertex.projVal_mp, left_endpoint.pts_mp[ii], pi); // set projection value
		vec_cp_mp(temp_vertex.pt_mp,left_endpoint.pts_mp[ii]);// set point
		temp_vertex.type = CRITICAL;// set type
		add_vertex(C,temp_vertex);
	}
	
	
	for (ii=0; ii<right_endpoint.num_pts; ii++) {
		printf("adding right endpoint %d of %d\n",ii, right_endpoint.num_pts);
		
		projection_value_homogeneous_input(temp_vertex.projVal_mp, right_endpoint.pts_mp[ii], pi); // set projection value

		vec_cp_mp(temp_vertex.pt_mp,right_endpoint.pts_mp[ii]);// set point
		temp_vertex.type = CRITICAL;// set type
		add_vertex(C,temp_vertex);
		printf("added point\n");
	}
	
	
	for (ii=0; ii< projections_sorted->size; ii++){
		printf("adding point %d of %d from crit_real to vertices\n",ii,projections_sorted->size);

		set_mp(temp_vertex.projVal_mp,  &projections_sorted->coord[ii]); // set projection value
		vec_cp_mp(temp_vertex.pt_mp,W_crit_real.pts_mp[index_tracker[ii]]);// set point
		temp_vertex.type = CRITICAL; // set type

		add_vertex(C,temp_vertex);
	}
	
	
	printf("there are a total of %d vertices\n",C->num_vertices);

	
	
	
	
	
	
	
	
	
	
	
	
	
  // 6) find points P
	// T = \pi(crit)
	// S = midpoints of T
  
	// for ( s \in S)
	//find P
	
	int edge_counter = 0; // set the counter
	
	OUT = safe_fopen_write("midpoints_upstairs");
	fprintf(OUT,"                                            \n\n%lf %d %d\n",midpoints_downstairs[ii],ii+1,num_midpoints);
	witness_set *midpoint_witness_sets;
	midpoint_witness_sets = (witness_set *)bmalloc(num_midpoints*sizeof(witness_set));
	
	for (ii=0; ii<num_midpoints; ii++) {
		printf("solving midpoints upstairs %d\n",ii);
		
		init_witness_set_d(&Wtemp);
		
		temp->i = 0;
		temp->r = -midpoints_downstairs[ii];
		d_to_mp(&particular_projection->coord[0],temp);
		
		//make the projection we will solve at.
		
		
		lin_to_lin_solver_main(solve_options->T.MPType,
													 W,
													 n_minusone_randomizer_matrix_full_prec,
													 &particular_projection, 1,
													 &Wtemp,
													 solve_options); //sets this value
		
		
		printf("sorting for realness\n");
		init_witness_set_d(&midpoint_witness_sets[ii]);
		sort_for_real(&midpoint_witness_sets[ii], Wtemp, solve_options->T);
		edge_counter += midpoint_witness_sets[ii].num_pts;
		
		//print the midpoints upstairs to the file.
		for (jj=0; jj<midpoint_witness_sets[ii].num_pts; jj++) { // iterate over all current new midpoint witness points
			
			
			
//			
//			add_vertex(C,temp_vertex);
			
			//print them to the screen
			init_vec_mp(dehom,W.num_variables);  dehom->size = W.num_variables;
			dehomogenize_mp(&dehom,midpoint_witness_sets[ii].pts_mp[jj]);
			for (kk=0; kk<W.num_variables-1; kk++) {
				mpf_out_str (OUT, 10, 16, dehom->coord[kk].r); fprintf(OUT," ");
				mpf_out_str (OUT, 10, 16, dehom->coord[kk].i); fprintf(OUT,"\n");
			}
			fprintf(OUT,"\n"); // print the line separating values
			clear_vec_mp(dehom);
		}
		fprintf(OUT,"\n");// print the line separating solutions
		clear_witness_set(Wtemp);
	}
	rewind(OUT);
	fprintf(OUT,"%d",edge_counter);
	fclose(OUT);
													 
	
	
	printf("done finding midpoints upstairs\n");

	
	
	//set C->V_0, C->V_1, C->E
	
	
	
  // 7) find edge endpoints
  
	
	//we assume the particular_projection is set above.
	
	witness_set Wleft, Wright;
	
	edge temp_edge;  init_edge(&temp_edge);

	for (ii=0; ii<num_midpoints; ++ii) {
		

		
		init_witness_set_d(&Wleft); init_witness_set_d(&Wright);
		
		
		temp->r = -crit_downstairs[ii]; temp->i = 0;
		d_to_mp(&particular_projection->coord[0],temp);
		
		comp_mp left_check; init_mp(left_check);
		comp_mp right_check; init_mp(right_check);
		
		d_to_mp(left_check,temp);
		
		lin_to_lin_solver_main(solve_options->T.MPType,
													 midpoint_witness_sets[ii], //the input
													 n_minusone_randomizer_matrix_full_prec,
													 &particular_projection, 1,
													 &Wleft,
													 solve_options); // the output
		
		
		temp->r = -crit_downstairs[ii+1]; temp->i = 0;
		d_to_mp(right_check,temp);
		
		d_to_mp(&particular_projection->coord[0],temp);
		
		
		lin_to_lin_solver_main(solve_options->T.MPType,
													 midpoint_witness_sets[ii], //the input
													 n_minusone_randomizer_matrix_full_prec,
													 &particular_projection, 1,
													 &Wright,
													 solve_options); // the output
		//each member of Wtemp should real.  if a member of V0 already, mark index.  else, add to V0, and mark.
	
//		printf("found %d points left\n",Wleft.num_pts);
//		printf("found %d points right\n",Wright.num_pts);
		
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
//		print_witness_set_to_screen(Wleft);
//		print_witness_set_to_screen(Wright);
		
		for (kk=0; kk<midpoint_witness_sets[ii].num_pts; kk++) {
			
			
			set_mp(temp_vertex.projVal_mp, &midpoint_witness_sets[ii].L_mp[0]->coord[0]  ); // set projection value
			vec_cp_mp(temp_vertex.pt_mp,midpoint_witness_sets[ii].pts_mp[kk]);// set point
			temp_vertex.type = MIDPOINT; // set type
			
			temp_edge.midpt = add_vertex(C,temp_vertex);
			
			printf("midpoint is vertex %d\n",temp_edge.midpt);
			
			
//			print_point_to_screen_matlab_mp(Wleft.pts_mp[kk],"left");
//			vec_cp_mp(temp_edge.midpt_mp,midpoint_witness_sets[ii].pts_mp[kk]);
			temp_edge.left  = index_in_vertices(C,Wleft.pts_mp[kk], left_check, solve_options->T,-1);
			temp_edge.right = index_in_vertices(C,Wright.pts_mp[kk],right_check,solve_options->T, 1);
			
			add_edge(C, temp_edge);
		}
		
		clear_mp(left_check); clear_mp(right_check);
		clear_witness_set(Wleft); clear_witness_set(Wright);
	}//for ii
	
	
	printf("C->num_edges = %d\n",C->num_edges);
	// track left/right
	// start at (x,t) = (p,1)
	// call endpoint q.
	//perform check on q \in crit,
  

	
	
	
	

	
	
	

	clear_witness_set(W_linprod);
	clear_witness_set(W_lintolin);
	clear_witness_set(W_detjacdetjac);
	
//TODO: clear memory
	return;
}



//code for checking that points satisfy the patch equation.
void check_patch_values(witness_set W){//, prog_t SLP, tracker_config_t T, mat_d n_minusone_randomizer_matrix, vec_d projection
	int ii;
	
	mat_d Jv, AtimesJ, tempmat;
	init_mat_d(Jv,1,1);init_mat_d(AtimesJ,1,1);init_mat_d(tempmat,W.num_variables,W.num_variables);
	
	
	tempmat->rows = tempmat->cols = W.num_variables; // change the size indicators
	
	mat_d Jv_Patch; init_mat_d(Jv_Patch,1,1);
	vec_d patchValues, parVals,parDer;  init_vec_d(patchValues,0);init_vec_d(parVals,0);init_vec_d(parDer,0);
	mat_d Jp; init_mat_d(Jp,1,1);
	
	comp_d pathVars;
	
	
	int patchType = 2;
	patch_eval_data_d patch;
	setupPatch_d(patchType, &patch, &W.num_variables, NULL);
	for (ii = 0; ii < W.num_variables ; ii++)
	{
		set_d(&patch.patchCoeff->entry[0][ii],&W.patch[0]->coord[ii]);
	}
	
	
	set_one_d(pathVars);
	
	printf("found %d points\n",W.num_pts);
	for (ii=0; ii<W.num_pts; ii++) {
		
		
		patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W.pts_d[ii], pathVars, &patch);  // Jp is ignored
		
		
		printf("point %d\n",ii);
		print_point_to_screen_matlab(patchValues,"patchvals[ii]");
		mypause();
		
		
	}
	
	return;
}


void check_detjac(witness_set W, prog_t SLP, tracker_config_t T, mat_d n_minusone_randomizer_matrix, vec_d projection)
{
	int ii, jj, kk;
	
	mat_d Jv, AtimesJ, tempmat;
	init_mat_d(Jv,1,1);init_mat_d(AtimesJ,1,1);init_mat_d(tempmat,W.num_variables,W.num_variables);
	
	
	tempmat->rows = tempmat->cols = W.num_variables; // change the size indicators
	
	comp_d detjac;
	mat_d Jv_Patch; init_mat_d(Jv_Patch,1,1);
	vec_d patchValues, parVals,parDer;  init_vec_d(patchValues,0);init_vec_d(parVals,0);init_vec_d(parDer,0);
	mat_d Jp; init_mat_d(Jp,1,1);
	
	comp_d pathVars;
	
	
	int patchType = 2;
	patch_eval_data_d patch;
	setupPatch_d(patchType, &patch, &W.num_variables, NULL);
	for (ii = 0; ii < W.num_variables ; ii++)
	{
		set_d(&patch.patchCoeff->entry[0][ii],&W.patch[0]->coord[ii]);
	}
	
	
	set_zero_d(pathVars);
	
	for (ii=0; ii<W.num_pts; ii++) {
		
		
		get_jacobian(W.pts_d[ii],T.MPType,W.num_var_gps,SLP,T,Jv);
		mat_mul_d(AtimesJ,n_minusone_randomizer_matrix,Jv); // randomize down.
		
		for (jj=0; jj<W.num_variables-2; jj++) {
			for (kk=0; kk<W.num_variables; kk++) {
				set_d(&tempmat->entry[jj][kk],&AtimesJ->entry[jj][kk]);
			}
		}
		//		increase_size_mat_d(tempmat,W.num_variables,W.num_variables); // make it bigger to accomodate more entries.  make square
		
		patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W.pts_d[ii], pathVars, &patch);  // Jp is ignored
		
		
		//copy in the projection from BED
		for (jj=0; jj<W.num_variables; ++jj) {
			set_d(&tempmat->entry[W.num_variables-2][jj],&projection->coord[jj]); // copy in the projection
		}
		
		//copy in the jacobian of the patch equation
		for (jj=0; jj<W.num_variables; ++jj) {
			set_d(&tempmat->entry[W.num_variables-1][jj],&Jv_Patch->entry[0][jj]); //  copy in the patch jacobian
		}
		
		//now TAKE THE DETERMINANT of tempmat.
		printf("point %d\n",ii);
		print_matrix_to_screen_matlab( tempmat,"detjac is det of");
		take_determinant_d(detjac,tempmat); // the determinant goes into detjac
		printf("detjac = %le+1i*%le\n",detjac->r,detjac->i);
		
	}
	return;
	
}








void get_jacobian(point_d current_values,
									int MPType,
									int num_var_gps,
									prog_t SLP,
									tracker_config_t T,
									mat_d jacobian)
{
	
	//initialize
	eval_struct_d e_d;
	init_eval_struct_d(e_d, 0, 0, 0);
	
	comp_d zerotime;
	set_zero_d(zerotime);
	
	
	evalProg_d(e_d.funcVals, e_d.parVals, e_d.parDer, jacobian, e_d.Jp, current_values, zerotime, &SLP);
	// the result of this, most importantly, is e_d.Jv, which contains the (complex) jacobian Jv for the system.
	// this jacobian Jv = \frac{\partial f_i}{\partial x_j} ( current_witness_point, zerotime)
	
	//TODO: NEED TO CLEAR THE EVALDATA, or leak memory.
	
	return;
	
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




void sort_for_membership(char * input_file,
												 witness_set *W_out,
												 witness_set W_in,
												 char *stifle_text){
//	printf("sorting points for membership\n");
	
	
	if (W_in.incidence_number==-1) {
		printf("input witness_set has unset incidence_number for comparison.\n");
		exit(-1112);
	}
	
	W_out->MPType = W_in.MPType;
	W_out->incidence_number = W_in.incidence_number;
	W_out->num_variables = W_in.num_variables;
	W_out->num_var_gps = W_in.num_var_gps;
	
	//copy the linears, patches from old to new
	cp_patches(W_out, W_in);
	cp_linears(W_out, W_in);
	
	//more important files out of the way.  this will probably crash the program if it is called without these files being present.
	rename_bertini_files_dotbak();
	
	
	int strLength = 0, digits = 15, *declarations = NULL;
  char *SysStr = NULL,*fmt = NULL, *bertini_command="bertini";
  FILE *IN = NULL;
	

	
	//make the command string to run
	strLength = 1 + snprintf(NULL, 0, "%s input_membership_test %s ", bertini_command, stifle_text);
  SysStr = (char *)bmalloc(strLength * sizeof(char));
  sprintf(SysStr, "%s input_membership_test %s ", bertini_command, stifle_text);
	
	
	
	// make the format to write the member_points file
  strLength = 1 + snprintf(NULL, 0, "%%.%dle %%.%dle\n", digits, digits);
  // allocate size
  fmt = (char *)bmalloc(strLength * sizeof(char));
  // setup fmt
  sprintf(fmt, "%%.%dle %%.%dle\n", digits, digits);
	
	
  // setup input file
	IN = safe_fopen_read(input_file);
  partitionParse(&declarations, IN, "func_input_real", "config_real",0); // the 0 means not self conjugate
	fclose(IN);
	
	
	//check existence of the required witness_data file.
	IN = safe_fopen_read("witness_data");
	fclose(IN);
	
	
	
	
	
	
	//only need to do this once.  we put the point and its conjugate into the same member points file and run the membership test simultaneously with one bertini call.
  membership_test_input_file("input_membership_test", "func_input_real", "config_real",3);
	

	
	
	// Do membership test
	printf("*\n%s\n*\n",SysStr);
  system(SysStr);
	
	
	int on_component_indicator[W_in.num_pts];
	read_incidence_matrix_wrt_number(on_component_indicator,W_in.incidence_number);
	
	printf("done reading incidence_matrix\n");
	
	int *is_on_component;
	is_on_component = (int *)bmalloc(W_in.num_pts*sizeof(int));
	

	
	int num_pts_on_component =0;
	int ii,jj;
	for (ii=0; ii<W_in.num_pts; ii++) {
//		printf("%d\n",on_component_indicator[ii]);
		
		if (on_component_indicator[ii]==1) {
			num_pts_on_component++;
			is_on_component[ii]=1;
		}
		else{
			is_on_component[ii]=0;
		}
		
	}
	
	if (num_pts_on_component==0) {
		printf("found 0 points out of %d candidates, on component.\n",W_in.num_pts);
		restore_bertini_files_dotbak();
		return;
	}
	
	
	vec_d *points_on_component;
	points_on_component =(point_d *)bmalloc(num_pts_on_component*sizeof(point_d));

	int *indices_on_component;
	indices_on_component = (int *)bmalloc(num_pts_on_component*sizeof(int));
	
	int counter =0;
	for (ii=0; ii<W_in.num_pts; ii++) {
		if (on_component_indicator[ii]==1) {
			init_vec_d(points_on_component[counter],W_in.pts_d[ii]->size);
			points_on_component[counter]->size = W_in.pts_d[ii]->size;
			vec_cp_d(points_on_component[counter],W_in.pts_d[ii]);
			
			indices_on_component[counter] = ii;

			counter++;
		}
	}
	
	
	
	//now remove duplicates

	
	
	int *is_unique; int curr_uniqueness; int num_good_pts = 0;
	is_unique  = (int *)bmalloc(num_pts_on_component*sizeof(int));
	for (ii = 0; ii<num_pts_on_component; ++ii) {
		curr_uniqueness = 1;
		if (ii!= (num_pts_on_component-1) ) { // the last point is unique...  always
			for (jj=ii+1; jj<num_pts_on_component; ++jj) {
				if (isSamePoint_homogeneous_input_d(points_on_component[ii],points_on_component[jj])){
					//formerly (isSamePoint(points_on_component[ii],NULL,52,points_on_component[jj],NULL,52,1e-8)){
					curr_uniqueness = 0;
//					printf("the following two points are not distinct:\n");
					print_point_to_screen_matlab(points_on_component[ii],"left");
					print_point_to_screen_matlab(points_on_component[jj],"right");
				}
			}
		}

				
		if (curr_uniqueness==1) {
			is_unique[ii] = 1;
			num_good_pts++;
		}
		else
		{
			is_unique[ii] = 0;
		}
		
	}

	
	

	
	
	
	
	//allocate the memory
	W_out->pts_d=(point_d *)bmalloc(num_good_pts*sizeof(point_d));
	W_out->pts_mp=(point_mp *)bmalloc(num_good_pts*sizeof(point_mp));
	W_out->num_pts = num_good_pts;
	W_out->num_pts = num_good_pts;
	
	
	
	//copy in the distinct points lying on the correct component
	counter = 0;
	for (ii=0; ii<num_pts_on_component; ++ii) {
		if (is_unique[ii]==1) {
			init_vec_d(W_out->pts_d[counter],W_in.num_variables); W_out->pts_d[counter]->size = W_in.num_variables;
			init_vec_mp2(W_out->pts_mp[counter],W_in.num_variables,1024);  W_out->pts_mp[counter]->size = W_in.num_variables;
			
			vec_cp_d(W_out->pts_d[counter], W_in.pts_d[indices_on_component[ii]]);
			vec_cp_mp(W_out->pts_mp[counter], W_in.pts_mp[indices_on_component[ii]]);
			counter++;
		}
	}
	
	if (!counter==num_good_pts){
		printf("counter mismatch; counter!=num_good_pts in sort_for_membership\n");
		exit(169);
	}

		//clear the memory
	
	for (ii=0; ii<num_pts_on_component; ++ii) {
		clear_vec_d(points_on_component[ii]);
	}
	free(points_on_component);
	

	free(is_unique);
  free(declarations);
	
  // delete temporary files
  remove("func_input_real");
  remove("config_real");
	remove("incidence_matrix");
//	remove("member_points");

	//move files back into place
	restore_bertini_files_dotbak();
	return;
}




//send in an initialized but empty witness_set.
// T is necessary for the tolerances.
void sort_for_unique(witness_set *W_out,
									 witness_set W_in,
									 tracker_config_t T)
{
	printf("sorting points for unique-ness\n%d points in \n",W_in.num_pts);
	int ii, jj;
	
	//	if (W_in.incidence_number==-1) {
	//		printf("input witness_set has unset incidence_number for comparison.\n");
	//		exit(-1112);
	//	}
	
	W_out->MPType = W_in.MPType;
	W_out->incidence_number = W_in.incidence_number;
	W_out->num_variables = W_in.num_variables;
	W_out->num_var_gps = W_in.num_var_gps;
	
	//copy the linears, patches from old to new
	cp_patches(W_out, W_in);
	cp_linears(W_out, W_in);
	
	//
	int curr_uniqueness;
	int num_good_pts = 0;
	int *is_unique  = (int *)bmalloc(W_in.num_pts*sizeof(int));
	
	for (ii = 0; ii<W_in.num_pts; ++ii) {
		curr_uniqueness = 1;
		if (ii!= (W_in.num_pts-1) ) { // the last point is unique...  always
			printf("a");
			for (jj=ii+1; jj<W_in.num_pts; ++jj) {
				if ( isSamePoint_homogeneous_input_d(W_in.pts_d[ii],W_in.pts_d[jj]) ){
					curr_uniqueness = 0;
					
					printf("the following two points are not distinct:\n");
					print_point_to_screen_matlab(W_in.pts_d[ii],"left");
					print_point_to_screen_matlab(W_in.pts_d[jj],"right");
				}
			}
		}
		
		printf("%d curr_uniqueness\n",curr_uniqueness);
		
		if (curr_uniqueness==1) {
			is_unique[ii] = 1;
			num_good_pts++;
		}
		else
		{
			is_unique[ii] = 0;
		}
		
	}

	
	W_out->num_pts = num_good_pts;
	
	W_out->pts_mp = (vec_mp *)bmalloc(num_good_pts*sizeof(vec_mp));
	W_out->pts_d = (vec_d *)bmalloc(num_good_pts*sizeof(vec_d));
	int counter = 0;
	for (ii=0; ii<W_in.num_pts; ++ii) {
		if (is_unique[ii]==1) {
			init_vec_d(W_out->pts_d[counter],W_in.num_variables); W_out->pts_d[counter]->size = W_in.num_variables;
			init_vec_mp2(W_out->pts_mp[counter],W_in.num_variables,1024);  W_out->pts_mp[counter]->size = W_in.num_variables;
			
			vec_cp_d(W_out->pts_d[counter], W_in.pts_d[ii]);
			vec_cp_mp(W_out->pts_mp[counter], W_in.pts_mp[ii]);
			counter++;
		}
	}
	
	if (counter!= num_good_pts) {
		printf("counter mismatch\n");
		exit(270);
	}
	
	free(is_unique);
	
	
	printf("%d points out\n",W_out->num_pts);
	return;
}




//send in an initialized but empty witness_set.
// T is necessary for the tolerances.
void sort_for_real(witness_set *W_out,
									 witness_set W_in,
									 tracker_config_t T)
{
//	printf("sorting points for real-ness\n%d points\n",W_in.num_pts);
	int ii;
	
//	if (W_in.incidence_number==-1) {
//		printf("input witness_set has unset incidence_number for comparison.\n");
//		exit(-1112);
//	}
	
	W_out->MPType = W_in.MPType;
	W_out->incidence_number = W_in.incidence_number;
	W_out->num_variables = W_in.num_variables;
	W_out->num_var_gps = W_in.num_var_gps;
	
	//copy the linears, patches from old to new
	cp_patches(W_out, W_in);
	cp_linears(W_out, W_in);
	
//	
	
	int real_indicator[W_in.num_pts];	
	int counter = 0;
	vec_mp result; init_vec_mp(result,1);
	for (ii=0; ii<W_in.num_pts; ii++) {
		dehomogenize_mp(&result, W_in.pts_mp[ii]);
		real_indicator[ii] = checkForReal_mp(result, T.real_threshold);
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	W_out->num_pts = counter;
	
	W_out->pts_mp = (point_mp *)bmalloc(counter*sizeof(point_mp));
	W_out->pts_d  = (point_d  *)bmalloc(counter*sizeof(point_d));
	
	
	
	counter = 0;
	for (ii=0; ii<W_in.num_pts; ii++) {
		if (real_indicator[ii]==1) {
			init_vec_mp(W_out->pts_mp[counter],W_in.num_variables); W_out->pts_mp[ii]->size = W_in.num_variables;
			init_vec_d( W_out->pts_d[counter],   W_in.num_variables); W_out->pts_d[ii]->size = W_in.num_variables;
			
			vec_cp_mp(W_out->pts_mp[counter],W_in.pts_mp[ii]);
			vec_mp_to_d(W_out->pts_d[counter],W_in.pts_mp[ii]);
			counter++;
		}
	}
	

	clear_vec_mp(result);
	
	
	return;
}





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
	int *degrees = (int *) bmalloc(num_funcs*sizeof(int));
	int *unique_degrees = (int *) bmalloc(num_funcs*sizeof(int));
	
	
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
	int *num_of_each_degree = (int *) bmalloc(num_unique_degrees*sizeof(int));
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

	print_matrix_to_screen_matlab_mp(randomization_matrix,"randomization");
	
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

