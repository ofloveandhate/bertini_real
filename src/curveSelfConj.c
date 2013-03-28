#include "curveSelfConj.h"










void computeCurveSelfConj(char * inputFile,
													witness_set_d W,
													vec_d pi,
													curveDecomp_d *C,
													int num_vars,
													int num_var_gps,
													unsigned int currentSeed){
	//IN DEVELOPMENT
  
	prog_t SLP;
	
	int MPType;
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
							&paramHom, MPType); // Opens the file "config" and stores data in tracker_config_t.
	
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
	vec_mp b_mp; init_vec_mp(b_mp,SLP_num_vars); b_mp->size = SLP_num_vars;
	vec_d b; init_vec_d(b,SLP_num_vars); b->size = SLP_num_vars;
	for (ii=0; ii<SLP_num_vars; ii++) {
		get_comp_rand_mp(&b_mp->coord[ii]);
	}
	vec_mp_to_d(b,b_mp);


	
	
	
	//need degree of g...  get from deg.out file
	int degprod = get_prod_degrees("deg.out",PPD.num_funcs);
	int deg_g = degprod - PPD.num_funcs;
	deg_g = deg_g - 1;
#ifdef verbose
	printf("the number of new linears will be %d\n",deg_g);
#endif
	
	


	
	vec_mp *new_linears_mp; new_linears_mp = (vec_mp *) malloc(deg_g*sizeof(vec_mp));
	vec_d *new_linears; new_linears = (vec_d *) malloc(deg_g*sizeof(vec_d));

	for (ii=0; ii<deg_g; ii++) {
		init_vec_mp(new_linears_mp[ii],SLP_num_vars);
		change_size_vec_mp(new_linears_mp[ii],SLP_num_vars); new_linears_mp[ii]->size = SLP_num_vars;
		for (kk=0; kk<SLP_num_vars; kk++) {
			get_comp_rand_mp(&new_linears_mp[ii]->coord[kk]);
		}
	}
	// copy the mp into the d
	for (ii=0; ii<deg_g; ii++) {
		init_vec_d(new_linears[ii],SLP_num_vars);
		change_size_vec_d(new_linears[ii],SLP_num_vars); new_linears[ii]->size = SLP_num_vars;
		vec_mp_to_d(new_linears[ii],new_linears_mp[ii]);
	}
	
	
//	print_point_to_screen_matlab(new_linears[0],"new_linears[0]");
	
	
	
	
	
	//end initialize
	
	
	
	
	
	
	//
	//
	//    actually start the computations
	//
	//
	
	
	
	

	
// 2) randomize down to N-1 equations
// to get a square system for the homotopies in the following steps.

//actually, just assign a random real nearly orthogonal matrix.

	
	mat_mp n_minusone_randomizer_matrix_mp;
	init_mat_mp(n_minusone_randomizer_matrix_mp,1,1);
	if (PPD.num_funcs>num_vars-1-num_var_gps){
		//  prepare a matrix, which we will use for left multiplication to get the correct dimensions
		make_matrix_random_real_mp(n_minusone_randomizer_matrix_mp,num_vars-1-num_var_gps,PPD.num_funcs,mpf_get_default_prec());
	}
	else if(PPD.num_funcs<num_vars-1-num_var_gps){
		printf("system has no one-dimensional components based on dimensions!\n");
		bexit(ERROR_CONFIGURATION);
	}
	else
	{ //no need, already have correct dimensions
		make_matrix_ID_mp(n_minusone_randomizer_matrix_mp, num_vars-1-num_var_gps, num_vars-1);
	}
	// copy the mp into the d
	mat_d n_minusone_randomizer_matrix;
	init_mat_d(n_minusone_randomizer_matrix,1,1);
	mat_mp_to_d(n_minusone_randomizer_matrix,n_minusone_randomizer_matrix_mp);
	
	
	
	
// 1) check rank
// a) use jacobian
	
	
	
	
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
    // is this initialization performed before -this- function?
  
	
	
  // 4) solve for critical conditions for random complex projection
    // if no box, determine box
//  int whatever = compute_box();
	
			// for each variable direction, compute projections of crit points onto that axis.
	
	witness_set_d Wtemp;
	init_witness_set_d(&Wtemp);
	cp_patches_d(&Wtemp,W); // copy the patches over from the original witness set
	
	witness_set_d W_lintolin;
	init_witness_set_d(&W_lintolin);
	cp_patches_d(&W_lintolin,W); // copy the patches over from the original witness set

    // if have box, intersect box with C
	if (T.MPType==1) {
		printf("running mptype1 solver");
		lin_to_lin_solver_mp(MPType, 0, currentSeed,  //mptype parse_time, current_seed
												 W,         // witness_set
												 n_minusone_randomizer_matrix_mp,
												 new_linears_mp, //  the set of linears we will solve at.
												 deg_g, // the number of new linears.
												 &Wtemp,
												 0, 1, 0);  //these numbers represent: myid, num_processes, headnode.

		
	}
	else {
		printf("running mptype0 solver");

		lin_to_lin_solver_d(MPType, 0, currentSeed,  //mptype parse_time, current_seed
												W,         // witness_set
												n_minusone_randomizer_matrix,
												new_linears, //  the set of linears we will solve at.
												deg_g, // the number of new linears.
												&Wtemp,
												0, 1, 0);  //these numbers represent: myid, num_processes, headnode.
	}
	printf("made it out of the lin_to_lin solver\n");
	
	//add the W to Wnew.
	merge_witness_sets(&W_lintolin,Wtemp,W);
	clear_witness_set(Wtemp);
	
	
	witness_set_d W_linprod;
	init_witness_set_d(&W_linprod);

	print_witness_set_to_screen(W_lintolin);

	if (T.MPType==2) {
		printf("bertini_real not fully equipped for mptype 2 yet.  please change to  0 or 1.\nyou may continue, but the program will probably crash.\n\n");
		mypause();
	}
	
	
	linprod_to_detjac_solver_d(MPType, 0, currentSeed,  //mptype parse_time, current_seed
														 W_lintolin,         // witness_set
														 n_minusone_randomizer_matrix,
														 b,  // the random complex vector we are using as a projection.
														 &W_linprod,
														 0, 1, 0);  //these numbers represent: myid, num_processes, headnode.
	
	
//	printf("made it out of linprod_to_detjac\n");
//	
//	print_witness_set_to_screen(W_linprod);
	write_dehomogenized_coordinates(W_linprod,"postdetjacsolns");
	mypause();
	
	
	
	mat_d Jv, AtimesJ, tempmat;
	init_mat_d(Jv,1,1);init_mat_d(AtimesJ,1,1);
	init_mat_d(tempmat,W.num_variables,W.num_variables);
	tempmat->rows = tempmat->cols = W.num_variables; // change the size indicators

	int patchType = 2;
	patch_eval_data_d patch;
	setupPatch_d(patchType, &patch, &T.numVars, NULL);
	for (ii = 0; ii < W.num_variables ; ii++)
	{
		set_d(&patch.patchCoeff->entry[0][ii],&W.patch[0]->coord[ii]);
	}
	
	comp_d detjac;
	
	mat_d Jv_Patch;
	init_mat_d(Jv_Patch,1,1);
	
	vec_d patchValues, parVals,parDer;
	init_vec_d(patchValues,0);init_vec_d(parVals,0);init_vec_d(parDer,0);
	mat_d Jp;
	init_mat_d(Jp,1,1);
	
	comp_d pathVars;
	set_zero_d(pathVars);
	
	printf("found %d points\n",W_linprod.W.num_pts);
	for (ii=0; ii<W_linprod.W.num_pts; ii++) {

		
		get_jacobian(W_linprod.W.pts[ii],MPType,num_var_gps,SLP,T,Jv);
		mat_mul_d(AtimesJ,n_minusone_randomizer_matrix,Jv); // randomize down.
		
		for (jj=0; jj<W.num_variables-2; jj++) {
			for (kk=0; kk<W.num_variables; kk++) {
				set_d(&tempmat->entry[jj][kk],&AtimesJ->entry[jj][kk]);
			}
		}
//		increase_size_mat_d(tempmat,W.num_variables,W.num_variables); // make it bigger to accomodate more entries.  make square
		
		patch_eval_d(    patchValues, parVals, parDer, Jv_Patch, Jp, W_linprod.W.pts[ii], pathVars, &patch);  // Jp is ignored
	
		
		//copy in the projection from BED
		for (jj=0; jj<W.num_variables; ++jj) {
			set_d(&tempmat->entry[W.num_variables-2][jj],&b->coord[jj]); // copy in the projection
		}
		
		//copy in the jacocian of the patch equation
		for (jj=0; jj<W.num_variables; ++jj) {
			set_d(&tempmat->entry[W.num_variables-1][jj],&Jv_Patch->entry[0][jj]); //  copy in the patch jacobian
		}
		
		//now TAKE THE DETERMINANT of tempmat.
		printf("point %d\n",ii);
		print_matrix_to_screen_matlab( tempmat,"detjac");
		take_determinant(detjac,tempmat); // the determinant goes into detjac
		printf("%le+1i*%le\n",detjac->r,detjac->i);
		mypause();
//#ifdef verbose
//		print_point_to_screen_matlab(W.W.pts[0], "W0");
//		
//		printf("\n\njacobian values at witness point 0:\n");
//		print_matrix_to_screen_matlab(Jv, "jacobian");
//		//working out of function_main.c inside function_eval_main
//		
//		print_point_to_screen_matlab(pi, "pi");
//		
//#endif

	}
	
//	linprod_to_detjac_solver_d();
//	printf("\nthere are now %d total witness points\n\nhere they are:\n",num_total_witness_points);
//	for (ii=0; ii<num_total_witness_points; ii++) {
//		printf("point %d has %d coordinates\n",ii,new_witness_points[ii]->size);
////		print_point_to_screen_matlab(new_witness_points[ii],"witnesspoint");
//	}
	
	//legacy notes:
//  real_curve_crit_pt_main_d(MPType, parse_time, currentSeed, startName, my_id, num_processes, headnode);
//	zero_dim_main(MPType, parse_time, currentSeed, startName, my_id, num_processes, headnode);
	
	clear_witness_set(W_linprod);
	clear_witness_set(W_lintolin);
	
	
  // 5) compute h_o, get critical points
  // isn't this what we are doing in 4)?

	
	

  // 6) find points P
    // T = \pi(crit)
    // S = midpoints of T
  
    // for ( s \in S)
      //find P
  
  
  // 7) find edge endpoints
  
    // track left/right
    // start at (x,t) = (p,1)
    // call endpoint q.
      //perform check on q \in crit, 
  
  //set C->V_0, C->V_1, C->E
	
	return;
}



//	point_d funcVals, point_d parVals;
//	vec_d parDer, mat_d Jv, mat_d Jp, point_d vars;
//	comp_d pathVars;
//	void const *ED;

//	*eval_func_d = &basic_eval_d;
//

//	eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, ED, eval_func_d);


//	eval_struct_d e;
//	init_eval_struct_d(e, 0, 0, 0);
////
//	basic_eval_data_d ED;
////
//	comp_d time;
//	time->r = 0;
//	time->i = 0;
//////	basic_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, PD->point, PD->time, void const *ED);
////
//	basic_eval_d(e.funcVals, e.parVals, e.parDer, e.Jv, e.Jp, Wnew.W.pts[0], time , ED);
////
//	clear_eval_struct_d(e);
//
//
//





void get_jacobian(point_d current_values,
									int MPType,
									int num_var_gps,
									prog_t SLP,
									tracker_config_t T,
									mat_d jacobian)
{
	
	int jj;  // counters.  i always use double letter for counters.
	
	
	
	// the prototype for eval_struct_d
//	// these contain all of the data structures that are returned when doing function evaluation
//	typedef struct
//	{
//		point_d funcVals;
//		point_d parVals;
//		vec_d parDer;
//		mat_d Jv;
//		mat_d Jp;
//	} eval_struct_d;
	
	
	//initialize
	eval_struct_d e_d;
	init_eval_struct_d(e_d, 0, 0, 0);
	
	comp_d zerotime;
	set_zero_d(zerotime);
	
	
	//initialize temp data
	point_d current_witness_point;
	init_point_d(current_witness_point,T.numVars);
	
//	//homogenize
//	for (jj=0; jj<num_var_gps; jj++) {
//		set_one_d(&current_witness_point->coord[jj]);
//	}
	
	for (jj=0; jj<T.numVars; jj++) {
		current_witness_point->coord[jj].r = current_values->coord[jj].r;
		current_witness_point->coord[jj].i = current_values->coord[jj].i;
	}
	
//	evalProg_d(temp_function_values, parVals, parDer, temp_jacobian_functions, temp_jacobian_parameters, current_variable_values, pathVars, BED->SLP);
	
	evalProg_d(e_d.funcVals, e_d.parVals, e_d.parDer, e_d.Jv, e_d.Jp, current_witness_point, zerotime, &SLP);
	// the result of this, most importantly, is e_d.Jv, which contains the (complex) jacobian Jv for the system.
	// this jacobian Jv = \frac{\partial f_i}{\partial x_j} ( current_witness_point, zerotime)
	
//	printf("returned from evalProg_d() Jv has %d rows, %d cols\n",e_d.Jv->rows,e_d.Jv->cols);
	
//	change_size_mat_d(jacobian,e_d.funcVals->size,T.numVars-num_var_gps);//initialize the jacobian to be the correct size
//	jacobian->cols = T.numVars-num_var_gps;
//	jacobian->rows = e_d.funcVals->size;     // why aren't these set in the change_size_mat_d???
//	
	
//	printf("%d funcvals\n",e_d.funcVals->size);
	
	
	//put into place from temp structures.
	mat_cp_d(jacobian,e_d.Jv);
//	for (kk = 0; kk < e_d.Jv->rows; kk++)
//	{ // print kth row
//		for (jj = num_var_gps; jj < e_d.Jv->cols; jj++)
//		{
//			{ // assign jth column
//				jacobian->entry[kk][jj-num_var_gps].r = e_d.Jv->entry[kk][jj].r;
//				jacobian->entry[kk][jj-num_var_gps].i = e_d.Jv->entry[kk][jj].i;
//			}
//		}
//	}
	
	return;
	
}






int get_prod_degrees(char filename[], int num_funcs){
	int degprod = 0, tmpdeg, ii;
	
	FILE *IN;
	
	IN =  safe_fopen_read(filename);
	
	for (ii = 0; ii<num_funcs; ii++) {
		fscanf(IN,"%d",&tmpdeg);
		degprod += tmpdeg;
	}

	fclose(IN);
	
	
	return degprod;
}





