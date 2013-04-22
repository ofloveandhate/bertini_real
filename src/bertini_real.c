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
  char *inputName = NULL, *witnessSetName = NULL,*input_deflated_Name=NULL;
  curveDecomp_d C;  //new data type; stores vertices, edges, etc.
  vec_mp pi_mp;  //random projection
  witness_set_d Wuser;
	int max_deflations = 10,strLength;
	int num_deflations, *deflation_sequence = NULL;
	int ii;  // counters
	
	

	////
	//  begin the actual program
	////
	
	
  startup(argC, args, &inputName, &witnessSetName);  //< prints the welcome message,
	//also gets the inputName, witnessSetName
	//default inputName = "input"
	//default witnessSetName = "witness_set"
  
	srand(0);
//	srand(time(NULL));
	// essentials for using the bertini parser
	prog_t SLP;
	
	
	unsigned int currentSeed;
	int trackType, genType = 0, MPType,  sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0;
  int my_id, num_processes, headnode = 0; // headnode is always 0
//	int precision = 53;
	num_processes = 1;
	int num_var_gps = 0, userHom = 0;
	
	//end parser-bertini essentials
	
	
	
	parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
	
	preproc_data PPD;
	setupPreProcData("preproc_data", &PPD);
	
	tracker_config_t T;
	get_tracker_config(&T,MPType);
	
	initMP(T.Precision);
	
//#ifdef _OPENMP
//#pragma omp parallel
//#endif
//  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
//    initMP(T.Precision);
//  }
	
	
	num_var_gps = PPD.num_var_gp;
	int orig_num_func = PPD.num_funcs;
	num_vars = setupProg(&SLP, T.Precision, MPType); // num_vars includes the number of homogeneous coordinates.
	// the number of homogeneous coordinates is the num_var_gps.
	
	
	printf("parsing witness set\n");
	init_witness_set_d(&Wuser);
	witnessSetParse(&Wuser,witnessSetName,num_vars);
	Wuser.num_var_gps = num_var_gps;
	Wuser.MPType = MPType;
	
	get_variable_names(&Wuser);

	
//	witnessSetParse(&Wnew, witnessSetName,num_vars);  // Wnew same as Wuser, except for functions. (needs to be updated)
//	Wnew.MPType = MPType;
//	Wnew.num_var_gps = num_var_gps;
	
	printf("getting component number\n");
	Wuser.incidence_number = get_component_number(Wuser,num_vars,inputName);
//	Wnew.incidence_number = Wuser.incidence_number;
	
	
//	printf("inserting new randomization matrix\n");
//	insert_randomization_matrix_witness_data(1,1,1);// move me to after isosingular deflation.
//	exit(0);
	
	
	
	write_dehomogenized_coordinates(Wuser, "witness_points_dehomogenized");
	

	printf("performing isosingular deflation\n");
	// perform an isosingular deflation
	rV = isosingular_deflation(&num_deflations, &deflation_sequence, inputName, "witness_points_dehomogenized", "bertini", "matlab -nosplash", max_deflations);
  
	strLength = 1 + snprintf(NULL, 0, "%s_comp_%d_deflated", inputName,deflation_sequence[num_deflations]);
	input_deflated_Name = (char *)bmalloc(strLength * sizeof(char));
	sprintf(input_deflated_Name, "%s_comp_%d_deflated", inputName,deflation_sequence[num_deflations]);
	
	
	
	
	parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	preproc_data_clear(&PPD);
	setupPreProcData("preproc_data", &PPD);
	int defl_num_func = PPD.num_funcs;
	
	



	
	
	
	printf("checking if component is self-conjugate\n");
	sc = checkSelfConjugate(Wuser,num_vars,inputName);  //later:  could be passed in from user, if we want
	if (sc==0) {
		printf("component is NOT self conjugate\n");
	}
	else{
		printf("component IS self conjugate\n");
	}
	
	
	
	
	
	
	
	
	
	//regenerate the various files, since we ran bertini since then
	parse_input(input_deflated_Name, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
	
	
	//now need to get new system produced by isosingular_deflation into bertini_real's memory.
	
	
	
	//need to convert Wuser to Wnew here
	
	//q: what could change?
	
	//temp answer:  the functions, but not the points, or the slices.
	
	
	// initialize the projection pi.  for now, get random projection.  would prefer to get it from a file.
	init_vec_mp(pi_mp,num_vars);
	pi_mp->size = num_vars; // should include the homogeneous variable
	
	set_zero_mp(&pi_mp->coord[0]);
	for (ii=1; ii<num_vars; ii++) {
		set_one_mp(&pi_mp->coord[ii]);
	}

	


	
	
	printf("init C\n");
	init_curveDecomp_d(&C);
	if (sc==0)  //C is not self-conjugate
	{
		//Call non-self-conjugate case code
		printf("\n\nentering not-self-conjugate case\n\n");
		computeCurveNotSelfConj(Wuser, pi_mp, &C,num_vars-1,"input_deflated");//This is Wenrui's !!!
		printf("Bertini_real found %d vertices (vertex)\n",C.num_V0);
	}
	else
	{
		//Call self-conjugate case code
		printf("\n\nentering self-conjugate case\n\n");
		computeCurveSelfConj(input_deflated_Name,Wuser,pi_mp,&C,num_vars,num_var_gps,currentSeed);  //This is Dans', at least at first !!!
	}
	
	
	
	printf("\n*\ndone with case\n*\n");

	printf("your projection was:\n");
	print_point_to_screen_matlab_mp(pi_mp,"pi");
//        set_zero_d(&(pi->coord[0]));
//        set_one_d( &(pi->coord[1]));
//        set_zero_d(&(pi->coord[2]));
		
	

	
	Output_Main("output", input_deflated_Name,Wuser.comp_num, num_vars, C, 2);
	
	
	printf("clearing data\n");
	// clear memory
	free(inputName);
	free(witnessSetName);
	
	clear_witness_set(Wuser);
//	clear_witness_set(Wnew);
	clear_curveDecomp_d(&C,2);
	
	
  //TMP END
  return 0;
}



