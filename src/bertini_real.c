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
  vec_d pi;  //random projection
  witness_set_d Wuser, Wnew;
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
  
	
	srand(time(NULL));
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
#ifdef _OPENMP
#pragma omp parallel
#endif
  { // set precision for each thread - all threads will execute this and set the precision correctly on each thread
    initMP(T.Precision);
  }
	
	
	num_var_gps = PPD.num_var_gp;
	num_vars = setupProg(&SLP, T.Precision, MPType); // num_vars includes the number of homogeneous coordinates.
	// the number of homogeneous coordinates is the num_var_gps.
	
	
	printf("parsing witness set\n");
	witnessSetParse(&Wuser,witnessSetName,num_vars);
	Wuser.MPType =MPType;
	witnessSetParse(&Wnew, witnessSetName,num_vars);  // Wnew same as Wuser, except for functions. (needs to be updated)
	Wnew.MPType =MPType;
	
	write_dehomogenized_coordinates(Wuser, "witness_points_dehomogenized");
	

	printf("performing isosingular deflation\n");
	// perform an isosingular deflation
	rV = isosingular_deflation(&num_deflations, &deflation_sequence, inputName, "witness_points_dehomogenized", "bertini", "matlab", max_deflations);
  
	 strLength = 1 + snprintf(NULL, 0, "%s_comp_%d_deflated", inputName,deflation_sequence[0]);
         input_deflated_Name = (char *)bmalloc(strLength * sizeof(char));
         sprintf(input_deflated_Name, "%s_comp_%d_deflated", inputName,deflation_sequence[0]);
	
	
	
	
	// initialize the projection pi.  for now, get random projection.  would prefer to get it from a file.
	init_vec_d(pi,num_vars);
	pi->size = num_vars;
	for (ii=0; ii<num_vars; ii++) {
		get_comp_rand_real_d(&pi->coord[ii]);
	}
	
	
	

	
	
	printf("checking if component is self-conjugate\n");
	sc = checkSelfConjugate(Wuser,num_vars,input_deflated_Name);  //later:  could be passed in from user, if we want
	if (sc==0) {
		printf("component is NOT self conjugate\n");
	}
	else{
		printf("component IS self conjugate\n");
	}
	
	
	
	
	
	
	
	
	
	
	parse_input(input_deflated_Name, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
	
	
	//now need to get new system produced by isosingular_deflation into bertini_real's memory.
	
	
	
	//need to convert Wuser to Wnew here
	
	//q: what could change?
	
	//temp answer:  the functions, but not the points, or the slices.
	
	
	printf("init C\n");
	init_curveDecomp_d(&C);
	if (sc==0)  //C is not self-conjugate
	{
		//Call non-self-conjugate case code
		printf("\n\nentering not-self-conjugate case\n\n");
		computeCurveNotSelfConj(Wnew, pi, &C,num_vars-1,"input_deflated");//This is Wenrui's !!!
		printf("Bertini_real found %d vertices (vertex)\n",C.num_V0);
	}
	else
	{
		//Call self-conjugate case code
		printf("\n\nentering self-conjugate case\n\n");
		computeCurveSelfConj(input_deflated_Name,Wnew,pi,&C,num_vars,num_var_gps,currentSeed);  //This is Dans', at least at first !!!
	}
	printf("\n*\ndone with case\n*\n");

	
        set_zero_d(&(pi->coord[0]));
        set_one_d(&(pi->coord[1]));
        set_zero_d(&(pi->coord[2]));
  
	
        printf("before V1"); fflush(stdout);
        C.num_V1=2;
        C.V1=(vertex_d *)bmalloc(C.num_V1*sizeof(vertex_d));

        init_point_d(C.V1[0].pt,num_vars);
        set_one_d(&(C.V1[0].pt->coord[0]));
        set_one_d(&(C.V1[0].pt->coord[1]));
        set_zero_d(&(C.V1[0].pt->coord[2]));

        init_point_d(C.V1[1].pt,num_vars);
        set_one_d(&(C.V1[1].pt->coord[0]));
        set_neg_one_d(&(C.V1[1].pt->coord[1]));
        set_zero_d(&(C.V1[1].pt->coord[2]));

        printf("after V1");fflush(stdout);
        C.num_E=2;

        C.E=(edge_d *)bmalloc(C.num_E*sizeof(edge_d));
        C.E[0].left=0;        C.E[0].right=1;
        init_point_d(C.E[0].midpt,num_vars);
        set_one_d(&(C.E[0].midpt->coord[0]));
        set_zero_d(&(C.E[0].midpt->coord[1]));
        set_one_d(&(C.E[0].midpt->coord[2]));
        init_point_d(C.E[1].pi,num_vars);
        point_cp_d(C.E[0].pi,pi);


        C.E[1].left=0;        C.E[1].right=1;
        init_point_d(C.E[1].midpt,num_vars);
        set_one_d(&(C.E[1].midpt->coord[0]));
        set_zero_d(&(C.E[1].midpt->coord[1]));
        set_neg_one_d(&(C.E[1].midpt->coord[2]));
        init_point_d(C.E[1].pi,num_vars);
        point_cp_d(C.E[1].pi,pi);

	
	Output_Main(inputName, input_deflated_Name,deflation_sequence, num_vars, C);
	
	
	// clear memory
	free(inputName);
	free(witnessSetName);
	
	clear_witness_set(Wuser);
	clear_witness_set(Wnew);
	
	
	
  //TMP END
  return 0;
}



