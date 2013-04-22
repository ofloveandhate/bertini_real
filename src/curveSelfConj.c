#include "curveSelfConj.h"










void computeCurveSelfConj(char * inputFile,
													witness_set_d W,
													vec_mp pi,
													curveDecomp_d *C,
													int num_vars,
													int num_var_gps,
													unsigned int currentSeed)
{
	//IN DEVELOPMENT
  
	prog_t SLP;
	
	int userHom = 0, pathMod = 0, paramHom = 0;
	double midpoint_tol = 0, intrinsicCutoffMultiplier = 0;
	
	int useRegen = 0, regenStartLevel = 0, maxCodim = 0, specificCodim = 0, reducedOnly = 0, supersetOnly = 0;
	//unused:  int useParameters = 0, timeIn = 0, num_variables = 0,
  int *startSub = NULL, *endSub = NULL, *startFunc = NULL, *endFunc = NULL, *startJvsub = NULL, *endJvsub = NULL, *startJv = NULL, *endJv = NULL, **subFuncsBelow = NULL;
  trackingStats trackCount;
	
	int ii,kk,jj;
	tracker_config_t T;
	
	
	init_trackingStats(&trackCount); // initialize trackCount to all 0
	
	FILE *IN = safe_fopen_read(inputFile);
	int *declarations = NULL;
	partitionParse(&declarations, IN, "func_input", "config",0); // the 0 means not self conjugate.
	
	// setup T
	setupConfig(&T, &midpoint_tol, &userHom, &useRegen, &regenStartLevel,
							&maxCodim, &specificCodim, &pathMod,
							&intrinsicCutoffMultiplier, &reducedOnly, &supersetOnly,
							&paramHom, T.MPType); // Opens the file "config" and stores data in tracker_config_t.
	
	preproc_data PPD;
	setupPreProcData("preproc_data", &PPD);
	
	
	// setup a SLP
	T.numVars = setupProg_count(&SLP, T.Precision, T.MPType, &startSub, &endSub, &startFunc, &endFunc, &startJvsub, &endJvsub, &startJv, &endJv, &subFuncsBelow);
	int SLP_num_vars = T.numVars;
#ifdef verbose
	printf("%d is precision\n%d is MPTYPE\n",T.Precision,T.MPType);
	printf("there are %d variables in the problem\nthere are %d vargrps\n", T.numVars,num_var_gps);
#endif
	
	
	
	
	//get the random complex projection, which we will use to find the crit pts
	vec_mp b_mp; init_vec_mp2(b_mp,SLP_num_vars,T.AMP_max_prec); b_mp->size = SLP_num_vars;
	vec_d b; init_vec_d(b,SLP_num_vars); b->size = SLP_num_vars;
	for (ii=0; ii<SLP_num_vars; ii++) {
		get_comp_rand_mp(&b_mp->coord[ii]);
	}
	vec_mp_to_d(b,b_mp);
	
	
	
	
	
	//need degree of g...  get from deg.out file
	int degsum = get_sum_degrees("deg.out",PPD.num_funcs);
	int deg_g = degsum - PPD.num_funcs;
	deg_g = deg_g - 1; // subtract one because already have one linear.
#ifdef verbose
	printf("the number of new linears will be %d\n",deg_g);
#endif
	
	
	
	
	
	vec_mp *new_linears_full_prec; new_linears_full_prec = (vec_mp *) malloc(deg_g*sizeof(vec_mp));
	
	for (ii=0; ii<deg_g; ii++) {
		init_vec_mp2(new_linears_full_prec[ii],SLP_num_vars,T.AMP_max_prec); new_linears_full_prec[ii]->size = SLP_num_vars;
		for (kk=0; kk<SLP_num_vars; kk++) {
			get_comp_rand_mp(&new_linears_full_prec[ii]->coord[kk]);
		}
	}
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	// 2) randomize down to N-1 equations
	// to get a square system for the homotopies in the following steps.
	
	//actually, just assign a random real nearly orthogonal matrix.
	
	
	mat_mp n_minusone_randomizer_matrix_full_prec;
	init_mat_mp2(n_minusone_randomizer_matrix_full_prec,W.num_variables-1-num_var_gps,PPD.num_funcs,T.AMP_max_prec);
	// note: AMP_max_prec = 1024;
	if (PPD.num_funcs>num_vars-1-num_var_gps){
		//  prepare a matrix, which we will use for left multiplication to get the correct dimensions
		make_matrix_random_real_mp(n_minusone_randomizer_matrix_full_prec,num_vars-1-num_var_gps,PPD.num_funcs,T.AMP_max_prec);
	}
	else if(PPD.num_funcs<num_vars-1-num_var_gps){
		printf("system does not have any one-dimensional components based on dimensions!\n");
		bexit(ERROR_CONFIGURATION);
	}
	else
	{ //no need, already have correct dimensions
		make_matrix_ID_mp(n_minusone_randomizer_matrix_full_prec, W.num_variables-1-num_var_gps, PPD.num_funcs);
	}
	
	// copy the mp into the d
	mat_d n_minusone_randomizer_matrix;
	init_mat_d(n_minusone_randomizer_matrix,n_minusone_randomizer_matrix_full_prec->rows,n_minusone_randomizer_matrix_full_prec->cols);
	n_minusone_randomizer_matrix->rows = n_minusone_randomizer_matrix_full_prec->rows;
	n_minusone_randomizer_matrix->cols = n_minusone_randomizer_matrix_full_prec->cols;
	mat_mp_to_d(n_minusone_randomizer_matrix,n_minusone_randomizer_matrix_full_prec);
	
	FILE *OUT=NULL;
	OUT = safe_fopen_write("Rand_matrix");
	fprintf(OUT,"%d\n",n_minusone_randomizer_matrix_full_prec->rows);
	fprintf(OUT,"%d\n",n_minusone_randomizer_matrix_full_prec->cols);
	fprintf(OUT,"%d\n",T.AMP_max_prec);
	printMat_mp(OUT, 0, n_minusone_randomizer_matrix_full_prec);
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
	
	witness_set_d Wtemp;
	init_witness_set_d(&Wtemp);
	cp_patches(&Wtemp,W); // copy the patches over from the original witness set
	
	witness_set_d W_lintolin;
	init_witness_set_d(&W_lintolin);
	cp_patches(&W_lintolin,W); // copy the patches over from the original witness set
	
	// if have box, intersect box with C
	
	lin_to_lin_solver_main(T.MPType,
												 W,
												 n_minusone_randomizer_matrix_full_prec,
												 new_linears_full_prec, deg_g,
												 &Wtemp);
	
	printf("made it out of the lin_to_lin solver\n");
	

	//add the W to Wnew.
	merge_witness_sets(&W_lintolin,Wtemp,W);

	clear_witness_set(Wtemp);
	
//	check_patch_values(W_lintolin);
	write_dehomogenized_coordinates(W_lintolin,"lintolin_points");
	
	
	witness_set_d W_linprod; init_witness_set_d(&W_linprod);

	witness_set_d W_detjacdetjac; init_witness_set_d(&W_detjacdetjac);
	

	cp_names(&W_lintolin,W);
	
	linprod_to_detjac_solver_main(T.MPType,
																W_lintolin,
																n_minusone_randomizer_matrix_full_prec,
																b_mp,
																&W_linprod);
	
	
	printf("made it out of linprod_to_detjac\n");
	//
	//	print_witness_set_to_screen(W_linprod);
	write_dehomogenized_coordinates(W_linprod,"member_points");
	write_dehomogenized_coordinates(W_linprod,"linprod_solns");
	witness_set_d W_linprod_good;
	init_witness_set_d(&W_linprod_good);
	W_linprod.incidence_number = W.incidence_number; // copy over from the input.  note that the incidence number may be different from the component number
	
	
	//	W_linprod_good should have only the points from W_in which lie on the correct component.
	sort_for_membership(inputFile, &W_linprod_good, W_linprod);
	write_dehomogenized_coordinates(W_linprod_good,"linprod_solns_postmembership");
	printf("done sorting membership\n");
	
//	print_witness_set_to_screen(W_linprod_good);
	
	
	


	
	
	
	detjac_to_detjac_solver_main(T.MPType,
															 W_linprod_good,
															 n_minusone_randomizer_matrix_full_prec,
															 b_mp,
															 pi,
															 &W_detjacdetjac);
	
//	print_witness_set_to_screen(W_detjacdetjac);
	write_dehomogenized_coordinates(W_detjacdetjac,"detjac_solns");

	
//	check_detjac( W_detjacdetjac, SLP, T, n_minusone_randomizer_matrix, pi_d);

	witness_set_d W_crit_real; init_witness_set_d(&W_crit_real);
	sort_for_real(&W_crit_real, W_detjacdetjac,T);
	
	

	
	printf("completed detjac2detjac solve\n\n");
	
	

	
  // 5) compute h_o, get critical points
	
	
	vec_mp projection_values; init_vec_mp(projection_values,W_crit_real.W_mp.num_pts); projection_values->size = W_crit_real.W_mp.num_pts;
	vec_mp dehom; init_vec_mp(dehom,W_crit_real.num_variables-1); dehom->size = W_crit_real.num_variables-1;
	
	vec_mp pi_no_hom; init_vec_mp(pi_no_hom,W_crit_real.num_variables-1); pi_no_hom->size = W_crit_real.num_variables-1;
	for (ii=0; ii<W_crit_real.num_variables-1; ii++) {
		set_mp(&pi_no_hom->coord[ii],&pi->coord[ii+1]);
	}
	
	
	for (ii=0; ii<W_crit_real.W_mp.num_pts; ii++) {
		dehomogenize_mp(&dehom,W_crit_real.W_mp.pts[ii]);
		dot_product_mp(&projection_values->coord[ii],dehom,pi_no_hom);
	}
	
	
		
	int *index_tracker = (int *)bmalloc(W_crit_real.W_mp.num_pts*sizeof(vec_mp)*sizeof(int));
	
	vec_mp projections_sorted;
	init_vec_mp(projections_sorted,W_crit_real.W_mp.num_pts); projections_sorted->size = W_crit_real.W_mp.num_pts;
	
	sort_increasing_by_real(&projections_sorted, &index_tracker, projection_values);
	
	
//	vec_mp *projections_sorted, vec_mp **points_sorted, int **index_tracker, vec_mp projections_input
	print_point_to_screen_matlab_mp(projection_values,"projection_values");
	print_point_to_screen_matlab_mp(projections_sorted,"sorted_projection_values");
	
	// projections_sorted
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	//5 find crit (the projection values)

	printf("***\n computing crit\n***");
	
	double left = mpf_get_d(projections_sorted->coord[0].r)-1;
	double right = mpf_get_d(projections_sorted->coord[W_crit_real.W_mp.num_pts-1].r)+1;
//TODO: make these '1' better;

	
	
	
	witness_set_d left_endpoint; init_witness_set_d(&left_endpoint);
	witness_set_d right_endpoint; init_witness_set_d(&right_endpoint);
	
	comp_d temp;
  vec_mp particular_projection;  init_vec_mp(particular_projection,W.num_variables); particular_projection->size = W.num_variables;
	vec_cp_mp(particular_projection,pi);

	
	/////
	//
	// go left
	//
	////////
	init_witness_set_d(&Wtemp);
	temp->r = -left; temp->i = 0;
	d_to_mp(&particular_projection->coord[0],temp);
	
	lin_to_lin_solver_main(T.MPType,
												 W,
												 n_minusone_randomizer_matrix_full_prec,
												 &particular_projection, 1,
												 &Wtemp);
	
	sort_for_real(&left_endpoint,Wtemp,T);
	int left_has_real_soln;
	if (left_endpoint.W.num_pts>0){
		left_has_real_soln = 1;
//		merge_witness_sets(mcreallerson);
	}
	else{
		left_has_real_soln = 0;
	}
	clear_witness_set(Wtemp); 
	
	/////
	//
	// go right
	//
	////////
	
	//set the value of the projection to which we wish to homotope.
	temp->r = -right; temp->i = 0;
	d_to_mp(&particular_projection->coord[0],temp);
	init_witness_set_d(&Wtemp);
	lin_to_lin_solver_main(T.MPType,
												 W,
												 n_minusone_randomizer_matrix_full_prec,
												 &particular_projection, 1,
												 &Wtemp);
	
	sort_for_real(&right_endpoint, Wtemp, T);
	int right_has_real_soln;
	if (right_endpoint.W.num_pts>0){
		right_has_real_soln = 1;
	}
	else{
		right_has_real_soln = 0;
	}
	
	clear_witness_set(Wtemp);
	
	

	
	
	

	/// dehomogenize the points to perform the projections using \pi.  the points we
	vec_d mcdehom;  init_vec_d(mcdehom,W.num_variables-1); mcdehom->size = W.num_variables-1;

	
	double *crit_downstairs; int num_crit;
	
	num_crit = W_crit_real.W_mp.num_pts + left_has_real_soln + right_has_real_soln;
	crit_downstairs = (double *)bmalloc(num_crit*sizeof(double));
	

	int num_midpoints = W_crit_real.W_mp.num_pts-1 + left_has_real_soln + right_has_real_soln;
	double *midpoints_downstairs = (double *) bmalloc(num_midpoints*sizeof(double));
	


	if (left_has_real_soln==1) {
		midpoints_downstairs[0] = (left+mpf_get_d(projections_sorted->coord[0].r))/2; //halfway between the endpoint and the leftmost critical point.
		crit_downstairs[0] = left;
	}

	
	int offset = left_has_real_soln;
	for (ii=0; ii<projections_sorted->size-1; ii++) {
		midpoints_downstairs[ii+offset] = (mpf_get_d(projections_sorted->coord[ii].r)+mpf_get_d(projections_sorted->coord[ii+1].r))/2;
	}
	for (ii=0; ii<W_crit_real.W_mp.num_pts; ii++) {
		crit_downstairs[ii+offset] = mpf_get_d(projections_sorted->coord[ii].r);
	}

	if (right_has_real_soln==1) {
		midpoints_downstairs[num_midpoints-1] = (right+ mpf_get_d(projections_sorted->coord[W_crit_real.W_mp.num_pts-1].r) )/2;
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
	

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
  // 6) find points P
	// T = \pi(crit)
	// S = midpoints of T
  
	// for ( s \in S)
	//find P
	
	int edge_counter = 0; // set the counter
	
	OUT = safe_fopen_write("midpoints_upstairs");
	fprintf(OUT,"                                            \n\n",midpoints_downstairs[ii],ii+1,num_midpoints);
	witness_set_d *midpoint_witness_sets;
	midpoint_witness_sets = (witness_set_d *)bmalloc(num_midpoints*sizeof(witness_set_d));
	
	for (ii=0; ii<num_midpoints; ii++) {
		init_witness_set_d(&midpoint_witness_sets[ii]);
		init_witness_set_d(&Wtemp);
		
		temp->i = 0;
		temp->r = -midpoints_downstairs[ii];
		d_to_mp(&particular_projection->coord[0],temp);
		
		//make the projection we will solve at.
		
		
		lin_to_lin_solver_main(T.MPType,
													 W,
													 n_minusone_randomizer_matrix_full_prec,
													 &particular_projection, 1,
													 &Wtemp); //sets this value
		
		
		
		sort_for_real(&midpoint_witness_sets[ii], Wtemp, T);
		edge_counter += midpoint_witness_sets[ii].W.num_pts;
		
		for (jj=0; jj<midpoint_witness_sets[ii].W.num_pts; jj++) {
			dehomogenize_mp(&dehom,midpoint_witness_sets[ii].W_mp.pts[jj]);
			for (kk=0; kk<W.num_variables-1; kk++) {
				mpf_out_str (OUT, 10, 16, dehom->coord[kk].r); fprintf(OUT," ");
				mpf_out_str (OUT, 10, 16, dehom->coord[kk].i); fprintf(OUT,"\n");
			}
			fprintf(OUT,"\n"); // print the line separating values
		}
		fprintf(OUT,"\n");// print the line separating solutions
		clear_witness_set(Wtemp);
	}
	rewind(OUT);
	fprintf(OUT,"%d",edge_counter);
	fclose(OUT);
													 
	
	
	
	
	
	
	
	//copy the points from witness sets into curve decomposition V0;
	
	
	vertex_d temp_vertex;  init_vertex_mp(&temp_vertex);
	
	for (ii=0; ii<left_endpoint.W.num_pts; ii++) {
		printf("adding left endpoints\n");
		dehomogenize_mp(&dehom,left_endpoint.W_mp.pts[ii]);
		dot_product_mp(temp_vertex.projVal_mp,dehom,pi_no_hom);
		
		vec_cp_mp(temp_vertex.pt_mp,left_endpoint.W_mp.pts[ii]);
		
		add_point_to_V1(C,temp_vertex);
	}
	
	
	for (ii=0; ii<right_endpoint.W.num_pts; ii++) {
		printf("adding right endpoints\n");
		dehomogenize_mp(&dehom,right_endpoint.W_mp.pts[ii]);
		dot_product_mp(temp_vertex.projVal_mp,dehom,pi_no_hom);
		
		vec_cp_mp(temp_vertex.pt_mp,right_endpoint.W_mp.pts[ii]);
		
		add_point_to_V1(C,temp_vertex);
	}

	
	for (ii=0; ii< projections_sorted->size; ii++){
		printf("adding point %d from crit_real to V1\n",ii);
		set_mp(temp_vertex.projVal_mp,  &projections_sorted->coord[ii]);
		
		
		vec_cp_mp(temp_vertex.pt_mp,W_crit_real.W_mp.pts[index_tracker[ii]]);
		add_point_to_V1(C,temp_vertex);
	}
	
	
	printf("there are a total of %d points to vertex set\n",C->num_V1);

	//set C->V_0, C->V_1, C->E
	
	
	
  // 7) find edge endpoints
  
	
	//we assume the particular_projection is set above.
	
	witness_set_d Wleft, Wright;
	
	edge_d temp_edge;  init_edge_mp(&temp_edge, W.num_variables);
	for (ii=0; ii<num_midpoints; ++ii) {
		

		
		init_witness_set_d(&Wleft); init_witness_set_d(&Wright);
		
		
		temp->r = -crit_downstairs[ii]; temp->i = 0;
		d_to_mp(&particular_projection->coord[0],temp);
		
		comp_mp left_check; init_mp(left_check);
		comp_mp right_check; init_mp(right_check);
		
		d_to_mp(left_check,temp);
		
		lin_to_lin_solver_main(T.MPType,
													 midpoint_witness_sets[ii], //the input
													 n_minusone_randomizer_matrix_full_prec,
													 &particular_projection, 1,
													 &Wleft); // the output
		
		temp->r = -crit_downstairs[ii+1]; temp->i = 0;
		d_to_mp(right_check,temp);
		
		d_to_mp(&particular_projection->coord[0],temp);
		
		
		lin_to_lin_solver_main(T.MPType,
													 midpoint_witness_sets[ii], //the input
													 n_minusone_randomizer_matrix_full_prec,
													 &particular_projection, 1,
													 &Wright); // the output
		//each member of Wtemp should real.  if a member of V0 already, mark index.  else, add to V0, and mark.
	
//		printf("found %d points left\n",Wleft.W_mp.num_pts);
//		printf("found %d points right\n",Wright.W_mp.num_pts);
		
		print_witness_set_to_screen(Wleft);
		print_witness_set_to_screen(Wright);
		
		vec_cp_mp(temp_edge.pi_mp, pi);
		for (kk=0; kk<midpoint_witness_sets[ii].W_mp.num_pts; kk++) {
//			print_point_to_screen_matlab_mp(Wleft.W_mp.pts[kk],"left");
			vec_cp_mp(temp_edge.midpt_mp,midpoint_witness_sets[ii].W_mp.pts[kk]);
			temp_edge.left  = index_in_V1(C,Wleft.W_mp.pts[kk], left_check, T,-1);
			temp_edge.right = index_in_V1(C,Wright.W_mp.pts[kk],right_check,T, 1);
			
			add_edge_mp(C, temp_edge);
		}
		
		clear_mp(left_check); clear_mp(right_check);
		clear_witness_set(Wleft); clear_witness_set(Wright);
	}//for ii
	
	
	printf("C->num_edges = %d\n",C->num_edges);
	// track left/right
	// start at (x,t) = (p,1)
	// call endpoint q.
	//perform check on q \in crit,
  

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
//	
//	
//	
//	vec_d pi; init_vec_d(pi,num_vars); pi->size = num_vars;
//	vec_mp_to_d(pi,pi_mp);
//	
//	C.num_V1=2;
//	C.V1=(vertex_d *)bmalloc(C.num_V1*sizeof(vertex_d));
//	
//	comp_d temp;
//	
//	//point = (1 0) for x^2+y^2-a^=0 x=a,y=0;
//	init_point_d(C.V1[0].pt,num_vars);
//	set_zero_d(&(C.V1[0].pt->coord[2]));
//	set_one_d(&(C.V1[0].pt->coord[0]));
//	add_d(temp,&(Wuser.patch[0]->coord[0]),&(Wuser.patch[0]->coord[1]));
//	div_d(&(C.V1[0].pt->coord[0]),&(C.V1[0].pt->coord[0]),temp);
//	set_d(&(C.V1[0].pt->coord[1]),&(C.V1[0].pt->coord[0]));
//	
//	//point=(-1 0) for x=-a,y=0
//	init_point_d(C.V1[1].pt,num_vars);
//	set_zero_d(&(C.V1[1].pt->coord[2]));
//	set_one_d(&(C.V1[1].pt->coord[0]));
//	sub_d(temp,&(Wuser.patch[0]->coord[0]),&(Wuser.patch[0]->coord[1]));
//	div_d(&(C.V1[1].pt->coord[0]),&(C.V1[1].pt->coord[0]),temp);
//	neg_d(&(C.V1[1].pt->coord[1]),&(C.V1[1].pt->coord[0]));
//	C.V1[1].pt->size=num_vars;
//	C.V1[0].pt->size=num_vars;
//	C.num_edges=2;
//	
//	C.E=(edge_d *)bmalloc(C.num_edges*sizeof(edge_d));
//	C.E[0].left=0;        C.E[0].right=1;
//	
//	//mid point (0 1)
//	init_point_d(C.E[0].midpt,num_vars);
//	C.E[0].midpt->size=num_vars;
//	set_zero_d(&(C.E[0].midpt->coord[1]));
//	set_one_d(&(C.E[0].midpt->coord[0]));
//	add_d(temp,&(Wuser.patch[0]->coord[0]),&(Wuser.patch[0]->coord[2]));
//	div_d(&(C.E[0].midpt->coord[0]),&(C.E[0].midpt->coord[0]),temp);
//	set_d(&(C.E[0].midpt->coord[2]),&(C.E[0].midpt->coord[0]));
//	init_point_d(C.E[0].pi,num_vars);
//	point_cp_d(C.E[0].pi,pi);
//	
//	C.E[1].left=0;        C.E[1].right=1;
//	//mid point (0 -1)
//	init_point_d(C.E[1].midpt,num_vars);
//	C.E[1].midpt->size=num_vars;
//	set_zero_d(&(C.E[1].midpt->coord[1]));
//	set_one_d(&(C.E[1].midpt->coord[0]));
//	sub_d(temp,&(Wuser.patch[0]->coord[0]),&(Wuser.patch[0]->coord[2]));
//	div_d(&(C.E[1].midpt->coord[0]),&(C.E[1].midpt->coord[0]),temp);
//	neg_d(&(C.E[1].midpt->coord[2]),&(C.E[1].midpt->coord[0]));
//	
//	init_point_d(C.E[1].pi,num_vars);
//	point_cp_d(C.E[1].pi,pi);
//	
//	
//	init_point_mp(C.V1[0].pt_mp,num_vars);
//	C.V1[0].pt_mp->size=num_vars;
//	point_d_to_mp(C.V1[0].pt_mp,C.V1[0].pt);
//	init_point_mp(C.V1[1].pt_mp,num_vars);
//	C.V1[1].pt_mp->size=num_vars;
//	point_d_to_mp(C.V1[1].pt_mp,C.V1[1].pt);
//	
//	
//	init_point_mp(C.E[0].midpt_mp,num_vars);
//	C.E[0].midpt_mp->size=num_vars;
//	point_d_to_mp(C.E[0].midpt_mp,C.E[0].midpt);
//	init_point_mp(C.E[0].pi_mp,num_vars);
//	C.E[0].pi_mp->size=num_vars;
//	point_d_to_mp(C.E[0].pi_mp,C.E[0].pi);
//	
//	init_point_mp(C.E[1].midpt_mp,num_vars);
//	C.E[1].midpt_mp->size=num_vars;
//	point_d_to_mp(C.E[1].midpt_mp,C.E[1].midpt);
//	C.E[1].midpt_mp->size=num_vars;
//	init_point_mp(C.E[1].pi_mp,num_vars);
//	C.E[1].pi_mp->size=num_vars;
//	point_d_to_mp(C.E[1].pi_mp,C.E[1].pi);
//	
//	
	
	
	

	clear_witness_set(W_linprod);
	clear_witness_set(W_lintolin);
	clear_witness_set(W_detjacdetjac);
	
//TODO: clear memory
	return;
}



//code for checking that points satisfy the patch equation.
void check_patch_values(witness_set_d W){//, prog_t SLP, tracker_config_t T, mat_d n_minusone_randomizer_matrix, vec_d projection
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
	
	
	set_one_d(pathVars);
	
	printf("found %d points\n",W.W.num_pts);
	for (ii=0; ii<W.W.num_pts; ii++) {
		
		
		patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W.W.pts[ii], pathVars, &patch);  // Jp is ignored
		
		
		printf("point %d\n",ii);
		print_point_to_screen_matlab(patchValues,"patchvals[ii]");
		mypause();
		
		
	}
	
	return;
}


void check_detjac(witness_set_d W, prog_t SLP, tracker_config_t T, mat_d n_minusone_randomizer_matrix, vec_d projection)
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
	
	for (ii=0; ii<W.W.num_pts; ii++) {
		
		
		get_jacobian(W.W.pts[ii],T.MPType,W.num_var_gps,SLP,T,Jv);
		mat_mul_d(AtimesJ,n_minusone_randomizer_matrix,Jv); // randomize down.
		
		for (jj=0; jj<W.num_variables-2; jj++) {
			for (kk=0; kk<W.num_variables; kk++) {
				set_d(&tempmat->entry[jj][kk],&AtimesJ->entry[jj][kk]);
			}
		}
		//		increase_size_mat_d(tempmat,W.num_variables,W.num_variables); // make it bigger to accomodate more entries.  make square
		
		patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W.W.pts[ii], pathVars, &patch);  // Jp is ignored
		
		
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
//		mypause();
		
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
												 witness_set_d *W_out,
												 witness_set_d W_in){
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
	rename("arr.out","arr.out.bak");
	rename("num.out","num.out.bak");
	rename("deg.out","deg.out.bak");
	rename("config","config.bak");
	rename("preproc_data","preproc_data.bak");
	
	
	int strLength = 0, digits = 15, *declarations = NULL;
  char *SysStr = NULL,*fmt = NULL, *bertini_command="bertini";
  FILE *IN = NULL, *OUT=NULL;
	

	
	//make the command string to run
	strLength = 1 + snprintf(NULL, 0, "%s input_membership_test", bertini_command);
  SysStr = (char *)bmalloc(strLength * sizeof(char));
  sprintf(SysStr, "%s input_membership_test", bertini_command);
	
	
	
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
	
	
	int on_component_indicator[W_in.W.num_pts];
	read_incidence_matrix_wrt_number(on_component_indicator,W_in.incidence_number);
	
	printf("done reading incidence_matrix\n");
	
	int *is_on_component;
	is_on_component = (int *)bmalloc(W_in.W.num_pts*sizeof(int));
	

	
	int num_pts_on_component =0;
	int ii,jj;
	for (ii=0; ii<W_in.W.num_pts; ii++) {
		printf("%d\n",on_component_indicator[ii]);
		
		if (on_component_indicator[ii]==1) {
			num_pts_on_component++;
			is_on_component[ii]=1;
		}
		else{
			is_on_component[ii]=0;
		}
		
	}
	
	if (num_pts_on_component==0) {
		printf("found 0 points out of %d candidates, on component.\n",W_in.W.num_pts);
		mypause();
	}
	
	
	vec_d *points_on_component;
	points_on_component =(point_d *)bmalloc(num_pts_on_component*sizeof(point_d));

	int *indices_on_component;
	indices_on_component = (int *)bmalloc(num_pts_on_component*sizeof(int));
	
	int counter =0;
	for (ii=0; ii<W_in.W.num_pts; ii++) {
		if (on_component_indicator[ii]==1) {
			init_vec_d(points_on_component[counter],W_in.W.pts[ii]->size);
			points_on_component[counter]->size = W_in.W.pts[ii]->size;
			vec_cp_d(points_on_component[counter],W_in.W.pts[ii]);
			
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
				if (isSamePoint(points_on_component[ii],NULL,52,points_on_component[jj],NULL,52,1e-8)){
					curr_uniqueness = 0;
					printf("the following two points are not unique:\n");
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
	W_out->W.pts=(point_d *)bmalloc(num_good_pts*sizeof(point_d));
	W_out->W_mp.pts=(point_mp *)bmalloc(num_good_pts*sizeof(point_mp));
	W_out->W.num_pts = num_good_pts;
	W_out->W_mp.num_pts = num_good_pts;
	
	
	
	//copy in the distinct points lying on the correct component
	counter = 0;
	for (ii=0; ii<num_pts_on_component; ++ii) {
		if (is_unique[ii]==1) {
			init_vec_d(W_out->W.pts[counter],W_in.num_variables); W_out->W.pts[counter]->size = W_in.num_variables;
			init_vec_mp2(W_out->W_mp.pts[counter],W_in.num_variables,1024);  W_out->W_mp.pts[counter]->size = W_in.num_variables;
			
			vec_cp_d(W_out->W.pts[counter], W_in.W.pts[indices_on_component[ii]]);
			vec_cp_mp(W_out->W_mp.pts[counter], W_in.W_mp.pts[indices_on_component[ii]]);
			counter++;
		}
	}
	
	printf("check: %d\n",counter==num_good_pts);

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
	rename("config.bak","config");
	rename("arr.out.bak","arr.out");
	rename("num.out.bak","num.out");
	rename("deg.out.bak","deg.out");
	rename("preproc_data.bak","preproc_data");
	return;
}








void sort_for_real(witness_set_d *W_out,
									 witness_set_d W_in,
									 tracker_config_t T)
{
	printf("sorting points for real-ness\n%d points mp %d points d\n",W_in.W_mp.num_pts,W_in.W.num_pts);
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
	
	int real_indicator[W_in.W.num_pts];	
	int counter = 0;
	vec_mp result; init_vec_mp(result,1);
	for (ii=0; ii<W_in.W_mp.num_pts; ii++) {
		dehomogenize_mp(&result, W_in.W_mp.pts[ii]);
//		print_point_to_screen_matlab_mp(result,"isreal?");
		real_indicator[ii] = checkForReal_mp(result, T.real_threshold);
//		printf("isreal is %d\n",real_indicator[ii]);
		if (real_indicator[ii]==1) {
			counter++;
		}
	}
	
	
	W_out->W_mp.pts = (point_mp *)bmalloc(counter*sizeof(point_mp));
	W_out->W.pts    = (point_d  *)bmalloc(counter*sizeof(point_d));
	W_out->W_mp.num_pts = W_out->W.num_pts = counter;
	for (ii=0; ii<W_in.W_mp.num_pts; ii++) {
		if (real_indicator[ii]==1) {
			init_vec_mp(W_out->W_mp.pts[ii],1); init_vec_d(W_out->W.pts[ii],1);
			vec_cp_mp(W_out->W_mp.pts[ii],W_in.W_mp.pts[ii]);
			vec_mp_to_d(W_out->W.pts[ii],W_in.W_mp.pts[ii]);
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
		min = 1e7;
		
		for (jj=0; jj<raw->size; jj++) {
			curr = mpf_get_d(raw->coord[jj].r);
			if ( curr < min) {
				printf("old min %lf, new min %lf\n",min,curr);
				indicator = jj;
				min = curr;
			}
		}
		if (indicator==-1) {
			printf("min projection value was insanely large\n");
			exit(1111);
		}
		
		
		
		(*index_tracker)[ii] = indicator;
		set_mp( &(*projections_sorted)->coord[ii],&raw->coord[indicator]);
		set_mp( &raw->coord[indicator],large);
	}
	return;
}





