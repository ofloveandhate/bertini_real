#include "determinant_derivative.h"








int take_determinant_d(comp_d determinant, mat_d source_matrix)
{
	
	
	
	if (source_matrix->cols!=source_matrix->rows) {
		printf("source matrix is not square! (%d rows, %d columns)\n",source_matrix->rows,source_matrix->cols);
		exit(-108);
	}
	if (source_matrix->cols==0) {
		printf("source matrix has 0 entries!");
		exit(-109);
	}
	
	int num_variables = source_matrix->cols;
	int ii;
	
	
	mat_d intermediate; init_mat_d(intermediate,0,0);
	
	int *rwnm = NULL; rwnm = NULL;
	vec_d garbage; init_vec_d(garbage,0);
	vec_d zerovec; init_vec_d(zerovec,0); change_size_vec_d(zerovec,num_variables); zerovec->size = num_variables;
	
	
	for (ii=0; ii<num_variables; ii++) {
		set_zero_d(&zerovec->coord[ii]);
	}
	
	
	
	double tol = 1e-20; //  these should be for realsies
	double largeChange = 1e20;
	
	// returns x, intermediate.
	
	//	print_matrix_to_screen_matlab(tempmat,"tempmat");
	int sign;
	
	int retval = LU_matrixSolve_d(garbage, intermediate, &rwnm, &sign, source_matrix, zerovec,tol,largeChange);
	//the solution is in intermediate.
	//error check.  solution failed if retval!=0
	if (retval!=0) {
		printf("LU decomposition failed (d)\n");
		print_matrix_to_screen_matlab(source_matrix,"source_matrix");
		exit(retval);
	}
	
	
	
	
	//compute the determinant
	set_one_d(determinant); // initialize
	for (ii=0; ii<num_variables; ii++) {
		mul_d(determinant,determinant,&intermediate->entry[rwnm[ii]][ii]);
	}
	determinant->r = determinant->r*sign;
	determinant->i = determinant->i*sign;
	
	//	print_matrix_to_screen_matlab(source_matrix,"detme");
	//	printf("candidate=%lf+1i*%lf;det(detme)\n",determinant->r,determinant->i);
	//	mypause();
	// this verified correct via 20 samples in matlab.  dab.
	
	
	clear_vec_d(garbage);
	clear_vec_d(zerovec);
	clear_mat_d(intermediate);
	
	return 0;
}





int take_determinant_mp(comp_mp determinant, mat_mp source_matrix){
	
	
	
	if (source_matrix->cols!=source_matrix->rows) {
		printf("source matrix is not square! (%d rows, %d columns)\n",source_matrix->rows,source_matrix->cols);
		exit(-108);
	}
	if (source_matrix->cols==0) { // no need to check rows -- know they are equal already
		printf("source matrix has 0 entries!");
		exit(-109);
	}
	
	int num_variables = source_matrix->cols;
	int ii;
	
	mat_mp intermediate; init_mat_mp(intermediate,num_variables,num_variables);
	intermediate->rows = intermediate->cols = num_variables;
	vec_mp zerovec; init_vec_mp(zerovec,source_matrix->cols); change_size_vec_mp(zerovec,num_variables); zerovec->size = num_variables;
	vec_mp garbage; init_vec_mp(garbage,source_matrix->cols); garbage->size = source_matrix->cols;
	
	int sign;
	int *rwnm = NULL; rwnm = NULL;
	
	mpf_t tol;  mpfr_init(tol); // = 1e-14
	mpf_t largeChange; mpfr_init(largeChange); //  = 1e11
	
	mpf_set_d(tol, 1e-20); // tol is the minimum acceptable 1-norm for each row during decomposition
	mpf_set_d(largeChange, 1e20);
	
	
	
	
	for (ii=0; ii<num_variables; ii++) {
		set_zero_mp(&zerovec->coord[ii]);
	}
	
	//  these should be for realsies
	
	// returns x, intermediate.
	
	//	print_matrix_to_screen_matlab(tempmat,"tempmat");
	
	
	int retval = LU_matrixSolve_mp(garbage, intermediate, &rwnm, &sign, source_matrix, zerovec,tol,largeChange);
	//the solution is in intermediate.
	//error check.  solution failed if retval!=0
	if (retval!=0) {
		printf("LU decomposition failed (mp)\n");
		print_matrix_to_screen_matlab_mp(source_matrix,"source_matrix");
		exit(retval);
	}
	
	
	
	
	//compute the determinant
	set_one_mp(determinant); // initialize
	for (ii=0; ii<num_variables; ii++) {
		mul_mp(determinant,determinant,&intermediate->entry[rwnm[ii]][ii]);
	}
	if (sign==-1)
	{
		neg_mp(determinant,determinant);
	}
	//	print_matrix_to_screen_matlab(source_matrix,"detme");
	//	printf("candidate=%lf+1i*%lf;det(detme)\n",determinant->r,determinant->i);
	//	mypause();
	// this verified correct via 20 samples in matlab.  dab.
	
	mpf_clear(tol);
	mpf_clear(largeChange);
	clear_mat_mp(intermediate);
	clear_vec_mp(zerovec);
	clear_vec_mp(garbage);
	return 0;
}




//																	patch_eval_data_d patch,


// output: Jv.
// input: current_variable_values, pathvars, ED
int detjac_numerical_derivative_d(mat_d Jv, //  the returned value
																	point_d current_variable_values, comp_d pathVars, vec_d projection,
																	int num_variables,
																	prog_t *SLP,
																	mat_d n_minusone_randomizer_matrix) // inputs
{
	int ii,jj,kk,mm; // counters
//	linprodtodetjac_eval_data_d *BED = (linprodtodetjac_eval_data_d *)ED; // cast the ED
	change_size_mat_d(Jv, 1, num_variables); // one row, num_variables columns
	Jv->rows = 1;
	Jv->cols = num_variables;
	
	//initialize some containers, for the unused stuff from the called evaluators.
	point_d unused_function_values, unused_parVals; init_vec_d(unused_function_values,0);init_vec_d(unused_parVals,0);
	vec_d unused_parDer; init_vec_d(unused_parDer,0);
	mat_d unused_Jp; init_mat_d(unused_Jp,0,0);
	
	
	//set up the perturbed jacobians
	mat_d perturbed_forward_Jv, perturbed_forward_Jv_Patch, perturbed_backward_Jv, perturbed_backward_Jv_Patch; //, base_Jv, base_Jv_Patch
	init_mat_d(perturbed_backward_Jv,0,0);init_mat_d(perturbed_forward_Jv,0,0);
	init_mat_d(perturbed_backward_Jv_Patch,0,0);init_mat_d(perturbed_forward_Jv_Patch,0,0);
	//	init_mat_d(base_Jv_Patch,0,0);init_mat_d(base_Jv,0,0);
	//the sizes for these variables will be set in the evalmethod calls below.  only need to init them here.
	
	//initialize the jacobians we will work with.
	point_d perturbed_forward_variables, perturbed_backward_variables;
	init_vec_d(perturbed_forward_variables,0); init_vec_d(perturbed_backward_variables,0);
	change_size_vec_d(perturbed_forward_variables,num_variables);
	change_size_vec_d(perturbed_backward_variables,num_variables);
	perturbed_backward_variables->size = perturbed_forward_variables->size = num_variables;
	
	comp_d perturbation;
	
	perturbation->r = PERTURBATION_VALUE;
	perturbation->i = PERTURBATION_VALUE;
	
	mat_d AtimesJ; init_mat_d(AtimesJ,1,1); // size will be properly set later.
	
	comp_d det_backward, det_forward;
	//	//get the baseline values at the current variable values.
	//	evalProg_d(  unused_function_values, unused_parVals, unused_parDer,  //  unused output
	//						 base_Jv,  // <---- the output we need
	//						 unused_Jp, //unused output
	//						 current_variable_values, pathVars, SLP); // input
	//
	//	patch_eval_d(unused_function_values, unused_parVals, unused_parDer, //  unused output
	//							 base_Jv_Patch, // <---- the output we need
	//							 unused_Jp, // unused output
	//							 current_variable_values, pathVars, &patch);  // input
	
	
	mat_d perturbed_forward_detme,perturbed_backward_detme; // create matrices
	init_mat_d(perturbed_forward_detme,num_variables-1,num_variables-1);
	init_mat_d(perturbed_backward_detme,num_variables-1,num_variables-1);
	
	perturbed_forward_detme->rows = perturbed_forward_detme->cols = perturbed_backward_detme->rows = perturbed_backward_detme->cols = num_variables-1; // set the size indicators
	set_zero_d(&Jv->entry[0][0]);
	//for each variable, we need the derivative.  we will put them in the Jv matrix.
	for (ii=0; ii<num_variables; ++ii) {
		
		//first, we assign the perturbed variables.  some of these calls could be eliminated
		for (jj=0; jj<num_variables; ++jj) {
			set_d(&perturbed_forward_variables->coord[jj],&current_variable_values->coord[jj]);
		}
		for (jj=0; jj<num_variables; ++jj) {
			set_d(&perturbed_backward_variables->coord[jj],&current_variable_values->coord[jj]);
		}
		add_d( &perturbed_forward_variables->coord[ii], &perturbed_forward_variables->coord[ii],perturbation);
		sub_d(&perturbed_backward_variables->coord[ii] ,&perturbed_backward_variables->coord[ii],perturbation);
		
		
		
		evalProg_d(unused_function_values, unused_parVals, unused_parDer,  //  unused output
							 perturbed_forward_Jv,  // <---- the output we need
							 unused_Jp, //unused output
							 perturbed_forward_variables, pathVars, SLP); // input
		
//		patch_eval_d(unused_function_values, unused_parVals, unused_parDer, //  unused output
//								 perturbed_forward_Jv_Patch, // <---- the output we need
//								 unused_Jp, // unused output
//								 perturbed_forward_variables, pathVars, &patch);  // input
		
		mat_mul_d(AtimesJ,n_minusone_randomizer_matrix,perturbed_forward_Jv);
		
		for (kk=0; kk<num_variables-2; ++kk) { // for each (randomized) equation
			for (mm=0; mm<num_variables-1; ++mm) { //  for each variable
				set_d(&perturbed_forward_detme->entry[kk][mm],&AtimesJ->entry[kk][mm+1]);
			}
		}
		for (kk=0; kk<num_variables-1;++kk) { // second to last row corresponds to the projection
			set_d(&perturbed_forward_detme->entry[num_variables-2][kk],&projection->coord[kk+1]);
		}
//		for (kk=0; kk<num_variables-1;++kk) {// the last row corresponds to the patch
//			set_d(&perturbed_forward_detme->entry[num_variables-1][kk],&perturbed_forward_Jv_Patch->entry[0][kk]);
//		}
		
		//now take the determinant;
		
		take_determinant_d(det_forward, perturbed_forward_detme);
		
		
		evalProg_d(  unused_function_values, unused_parVals, unused_parDer,  //  unused output
							 perturbed_backward_Jv,  // <---- the output we need
							 unused_Jp, //unused output
							 perturbed_backward_variables, pathVars, SLP); // input
		
//		patch_eval_d(unused_function_values, unused_parVals, unused_parDer, //  unused output
//								 perturbed_backward_Jv_Patch, // <---- the output we need
//								 unused_Jp, // unused output
//								 perturbed_backward_variables, pathVars, &patch);  // input
		//now have the pieces of the puzzle.
		
		mat_mul_d(AtimesJ,n_minusone_randomizer_matrix,perturbed_backward_Jv);
		
		//copy in the jacobian WRT variables
		for (kk=0; kk<num_variables-2; ++kk) { // for each (randomized) equation
			for (mm=0; mm<num_variables-1; ++mm) { // for each variable
				set_d(&perturbed_backward_detme->entry[kk][mm],&AtimesJ->entry[kk][mm+1]);
			}
		}
		//copy in the projection
		for (kk=0; kk<num_variables-1;++kk) {
			set_d(&perturbed_backward_detme->entry[num_variables-2][kk],&projection->coord[kk+1]);
		}
//		//copy in the patch
//		for (kk=0; kk<num_variables;++kk) {
//			set_d(&perturbed_backward_detme->entry[num_variables-1][kk],&perturbed_backward_Jv_Patch->entry[0][kk]);
//		}
		
		//now take the determinant;
		
		take_determinant_d(det_backward, perturbed_backward_detme);
		
		sub_d(&Jv->entry[0][ii],det_forward,det_backward);  // Jv = (forward-backward)/(2h)
		div_d(&Jv->entry[0][ii],&Jv->entry[0][ii],perturbation);
		Jv->entry[0][ii].r /= 2.0;
		Jv->entry[0][ii].i /= 2.0;
		
		//		print_matrix_to_screen_matlab(perturbed_forward_detme,"forwarddetme");
		//		print_matrix_to_screen_matlab(perturbed_backward_detme,"backwarddetme");
		//		printf("h=%lf+1i*%lf;\n",perturbation->r,perturbation->i);
		//		printf("%lf+1i*%lf\n",Jv->entry[0][ii].r,Jv->entry[0][ii].i);
		//		mypause();
	}
	
	clear_mat_d(AtimesJ);
	clear_mat_d(perturbed_backward_detme);clear_mat_d(perturbed_forward_detme);
//TODO: HERE forgetting to clear some vectors.  losing memory.
	return 0;
}


//																	 patch_eval_data_mp patch,

// output: Jv.
// input: current_variable_values, pathvars, ED
int detjac_numerical_derivative_mp(mat_mp Jv, //  the returned value
																	 point_mp current_variable_values, comp_mp pathVars, vec_mp projection,
																	 int num_variables,
																	 prog_t *SLP,
																	 mat_mp n_minusone_randomizer_matrix) // inputs
{
	
	int ii,jj,kk,mm; // counters
//	linprodtodetjac_eval_data_mp *BED = (linprodtodetjac_eval_data_mp *)ED; // cast the ED
	
	change_size_mat_mp(Jv, 1, num_variables); // one row, num_variables columns
	Jv->rows = 1;
	Jv->cols = num_variables;
	
	//initialize some containers, for the unused stuff from the called evaluators.
	point_mp unused_function_values, unused_parVals; init_vec_mp(unused_function_values,0);init_vec_mp(unused_parVals,0);
	vec_mp unused_parDer; init_vec_mp(unused_parDer,1);
	mat_mp unused_Jp; init_mat_mp(unused_Jp,0,0); // these 0's should be made what they will be
	
	
	//set up the perturbed jacobians
	mat_mp perturbed_forward_Jv, perturbed_forward_Jv_Patch, perturbed_backward_Jv, perturbed_backward_Jv_Patch; //, base_Jv, base_Jv_Patch
	init_mat_mp(perturbed_backward_Jv,0,0);init_mat_mp(perturbed_forward_Jv,0,0); // change these sizes
	init_mat_mp(perturbed_backward_Jv_Patch,0,0);init_mat_mp(perturbed_forward_Jv_Patch,0,0);
	//	init_mat_mp(base_Jv_Patch,0,0);init_mat_mp(base_Jv,0,0);
	//the sizes for these variables will be set in the evalmethod calls below.  only need to init them here.
	
	//initialize the jacobians we will work with.
	point_mp perturbed_forward_variables, perturbed_backward_variables;
	init_vec_mp(perturbed_forward_variables,num_variables); init_vec_mp(perturbed_backward_variables,num_variables);
	perturbed_backward_variables->size = perturbed_forward_variables->size = num_variables;
	
	comp_mp perturbation; init_mp(perturbation);
	comp_d p; p->r = PERTURBATION_VALUE_mp; p->i = PERTURBATION_VALUE_mp;
	d_to_mp(perturbation,p);
	
	comp_mp temp, temp2; init_mp(temp); init_mp(temp2);
	comp_mp two; init_mp(two);
	comp_d t; t->r = 2.0; t->i = 0;
	d_to_mp(two,t);
	
	mat_mp AtimesJ; init_mat_mp(AtimesJ,1,1); // size will be properly set later.
	
	comp_mp det_backward, det_forward; init_mp(det_backward); init_mp(det_forward);
	//	//get the baseline values at the current variable values.
	//	evalProg_mp(  unused_function_values, unused_parVals, unused_parDer,  //  unused output
	//						 base_Jv,  // <---- the output we need
	//						 unused_Jp, //unused output
	//						 current_variable_values, pathVars, SLP); // input
	//
	//	patch_eval_mp(unused_function_values, unused_parVals, unused_parDer, //  unused output
	//							 base_Jv_Patch, // <---- the output we need
	//							 unused_Jp, // unused output
	//							 current_variable_values, pathVars, &patch);  // input
	
	mat_mp perturbed_forward_detme,perturbed_backward_detme; // create matrices
	init_mat_mp(perturbed_forward_detme,num_variables-1,num_variables-1);
	init_mat_mp(perturbed_backward_detme,num_variables-1,num_variables-1);
	
	perturbed_forward_detme->rows = perturbed_forward_detme->cols = perturbed_backward_detme->rows = perturbed_backward_detme->cols = num_variables-1; // set the size indicators
	
//	set_zero_mp(&Jv->entry[0][0]);
	//for each variable, we need the derivative.  we will put them in the Jv matrix.
	for (ii=0; ii<num_variables; ++ii) {
		
		//first, we assign the perturbed variables.  some of these calls could be eliminated
		for (jj=0; jj<num_variables; ++jj) {
			set_mp(&perturbed_forward_variables->coord[jj],&current_variable_values->coord[jj]);
		}
		for (jj=0; jj<num_variables; ++jj) {
			set_mp(&perturbed_backward_variables->coord[jj],&current_variable_values->coord[jj]);
		}
		add_mp( &perturbed_forward_variables->coord[ii], &perturbed_forward_variables->coord[ii],perturbation);
		sub_mp(&perturbed_backward_variables->coord[ii] ,&perturbed_backward_variables->coord[ii],perturbation);
		
		
		
		evalProg_mp(unused_function_values, unused_parVals, unused_parDer,  //  unused output
								perturbed_forward_Jv,  // <---- the output we need
								unused_Jp, //unused output
								perturbed_forward_variables, pathVars, SLP); // input
		
//		patch_eval_mp(unused_function_values, unused_parVals, unused_parDer, //  unused output
//									perturbed_forward_Jv_Patch, // <---- the output we need
//									unused_Jp, // unused output
//									perturbed_forward_variables, pathVars, &patch);  // input
		
		mat_mul_mp(AtimesJ,n_minusone_randomizer_matrix,perturbed_forward_Jv);
		
		for (kk=0; kk<num_variables-2; ++kk) { // for each (randomized) equation
			for (mm=0; mm<num_variables-1; ++mm) { //  for each variable
				set_mp(&perturbed_forward_detme->entry[kk][mm],&AtimesJ->entry[kk][mm+1]);
			}
		}
		for (kk=0; kk<num_variables-1;++kk) { // second to last row corresponds to the projection
			set_mp(&perturbed_forward_detme->entry[num_variables-2][kk],&projection->coord[kk+1]);
		}
//		for (kk=0; kk<num_variables;++kk) {// the last row corresponds to the patch
//			set_mp(&perturbed_forward_detme->entry[num_variables-1][kk],&perturbed_forward_Jv_Patch->entry[0][kk]);
//		}
		
		//now take the determinant;
		
		take_determinant_mp(det_forward, perturbed_forward_detme);
		
		
		evalProg_mp(  unused_function_values, unused_parVals, unused_parDer,  //  unused output
								perturbed_backward_Jv,  // <---- the output we need
								unused_Jp, //unused output
								perturbed_backward_variables, pathVars, SLP); // input
		
//		patch_eval_mp(unused_function_values, unused_parVals, unused_parDer, //  unused output
//									perturbed_backward_Jv_Patch, // <---- the output we need
//									unused_Jp, // unused output
//									perturbed_backward_variables, pathVars, &patch);  // input
		//now have the pieces of the puzzle.
		
		mat_mul_mp(AtimesJ,n_minusone_randomizer_matrix,perturbed_backward_Jv);
		
		//copy in the jacobian WRT variables
		for (kk=0; kk<num_variables-2; ++kk) { // for each (randomized) equation
			for (mm=0; mm<num_variables-1; ++mm) { // for each variable
				set_mp(&perturbed_backward_detme->entry[kk][mm],&AtimesJ->entry[kk][mm+1]);
			}
		}
		//copy in the projection
		for (kk=0; kk<num_variables-1;++kk) {
			set_mp(&perturbed_backward_detme->entry[num_variables-2][kk],&projection->coord[kk+1]);
		}
		//copy in the patch
//		for (kk=0; kk<num_variables;++kk) {
//			set_mp(&perturbed_backward_detme->entry[num_variables-1][kk],&perturbed_backward_Jv_Patch->entry[0][kk]);
//		}
		
		//now take the determinant;
		
		take_determinant_mp(det_backward, perturbed_backward_detme);
		sub_mp(temp,det_forward,det_backward);  // Jv = (forward-backward)/(2h)
		div_mp(temp2,temp,perturbation);
		div_mp(&Jv->entry[0][ii],temp2,two);
		
		//		print_matrix_to_screen_matlab(perturbed_forward_detme,"forwarddetme");
		//		print_matrix_to_screen_matlab(perturbed_backward_detme,"backwarddetme");
		//		printf("h=%lf+1i*%lf;\n",perturbation->r,perturbation->i);
		//		printf("%lf+1i*%lf\n",Jv->entry[0][ii].r,Jv->entry[0][ii].i);
		//		mypause();
	}
	
	clear_point_mp(unused_function_values);
	clear_point_mp(unused_parVals);
	clear_vec_mp(unused_parDer);
	clear_mat_mp(unused_Jp);
	
	clear_mat_mp(perturbed_forward_Jv);
	clear_mat_mp(perturbed_forward_Jv_Patch);
	clear_mat_mp(perturbed_backward_Jv);
	clear_mat_mp(perturbed_backward_Jv_Patch);
	
	clear_point_mp(perturbed_forward_variables);
	clear_point_mp(perturbed_backward_variables);
	
	clear_mp(perturbation);
	clear_mp(two);
	
	clear_mat_mp(AtimesJ);
	clear_mp(det_backward);
	clear_mp(det_forward);
	clear_mp(temp); clear_mp(temp2);
	
	clear_mat_mp(perturbed_backward_detme);clear_mat_mp(perturbed_forward_detme);
	
	
	return 0;
}







