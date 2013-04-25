#include "bertini_real.h"

int main(int argC, char *args[])
/***************************************************************\
 * USAGE:                                                        *
 * ARGUMENTS:                                                    *
 * RETURN VALUES:                                                *
 * NOTES:                                                        *
 \***************************************************************/
{
  int rV,num_vars=0,sc;  //1=self-conjugate; 0=not

  curveDecomp_d C;  //new data type; stores vertices, edges, etc.
  vec_mp pi_mp;  //random projection
  witness_set Wuser;
	int max_deflations = 10,strLength;
	int num_deflations, *deflation_sequence = NULL;
	int MPType;
	

	////
	//  INITIALIZATION
	////
	
	splash_screen();
	
	//instantiate options
	program_configuration program_options;  init_program_config(&program_options);
	parse_options(argC, args, &program_options);
	
	//parse the options
	startup(program_options); // tests for existence of necessary files, etc.
	
	//if desired, display the options
//	display_current_options(program_options);
	
  
	srand(0);
//	srand(time(NULL));
	// essentials for using the bertini parser

	
	
	//split the input_file
	parse_input_file(program_options.input_filename, &MPType);
	
	
	// set up the solver configuration
	solver_configuration solve_options;  init_solver_config(&solve_options);
	get_tracker_config(&solve_options,MPType);
	setupPreProcData("preproc_data", &solve_options.PPD);
	
	initMP(solve_options.T.Precision);
	

	
	
	num_vars = get_num_vars_PPD(solve_options.PPD);
	

	init_witness_set_d(&Wuser);
	witnessSetParse(&Wuser,program_options.witness_set_filename,num_vars);
	Wuser.num_var_gps = solve_options.PPD.num_var_gp;
	Wuser.MPType = MPType;
	
	get_variable_names(&Wuser);


	Wuser.incidence_number = get_component_number(Wuser,num_vars,
																								program_options.input_filename,
																								program_options.stifle_text);
	
	

	printf("performing isosingular deflation\n");
	// perform an isosingular deflation
	write_dehomogenized_coordinates(Wuser, "witness_points_dehomogenized"); // write the points to file
	
	rV = isosingular_deflation(&num_deflations, &deflation_sequence, program_options.input_filename, "witness_points_dehomogenized", "bertini", "matlab -nosplash", max_deflations);
  
	//this should be made a function
	strLength = 1 + snprintf(NULL, 0, "%s_comp_%d_deflated", program_options.input_filename,deflation_sequence[num_deflations]);
	sprintf(program_options.input_deflated_filename, "%s_comp_%d_deflated", program_options.input_filename,deflation_sequence[num_deflations]);
	
	
	
	
	// this wraps around a bertini routine
	parse_input_file(program_options.input_deflated_filename, &MPType);
	
	preproc_data_clear(&solve_options.PPD);
	setupPreProcData("preproc_data", &solve_options.PPD);

	
	

	
	
	printf("checking if component is self-conjugate\n");
	sc = checkSelfConjugate(Wuser,num_vars,program_options.input_filename, program_options.stifle_text);  //later:  could be passed in from user, if we want

	
	
	
	//regenerate the various files, since we ran bertini since then.  
	parse_input_file(program_options.input_deflated_filename, &MPType);
	
	
	//get the random projection \pi
	init_vec_mp(pi_mp,num_vars); pi_mp->size = num_vars; // should include the homogeneous variable
	get_projection(pi_mp, program_options, solve_options, num_vars);
	

	
	

	

	//initialize the data structure which collets the output
	init_curveDecomp_d(&C);
	
	

	
	if (sc==0)  //C is not self-conjugate
	{
		//Call non-self-conjugate case code
		printf("\n\nentering not-self-conjugate case\n\n");
		computeCurveNotSelfConj(Wuser, pi_mp, &C,num_vars-1,program_options.input_deflated_filename);//This is Wenrui's !!!
		printf("Bertini_real found %d vertices (vertex)\n",C.num_V0);
	}
	else
	{
		//Call self-conjugate case code
		printf("\n\nentering self-conjugate case\n\n");
		computeCurveSelfConj(program_options.input_deflated_filename,Wuser,pi_mp,&C,
												 num_vars,Wuser.num_var_gps,
												 &program_options, &solve_options);  //This is Dans', at least at first !!!
	}
	
	
	
	printf("\n*\ndone with case\n*\n");

		
	
	C.num_variables = num_vars;
	if (MPType==0) {
		init_vec_d(C.pi_d, Wuser.num_variables); C.pi_d->size = Wuser.num_variables;
		vec_mp_to_d(C.pi_d, pi_mp);
	}
	else
	{
		init_vec_mp(C.pi, Wuser.num_variables); C.pi->size = Wuser.num_variables;
		vec_cp_mp(C.pi, pi_mp);
	}

	printf("outputting data\n");
	Output_Main("output", program_options.input_deflated_filename,Wuser.comp_num, num_vars, C, 2);
	
	

	printf("clearing witness_set\n");
	clear_witness_set(Wuser);

	printf("clearing C\n");
	clear_curveDecomp_d(&C,2);
	
//	printf("clearing program_options\n");
//	clear_program_config(&program_options);
	
  //TMP END
  return 0;
}



