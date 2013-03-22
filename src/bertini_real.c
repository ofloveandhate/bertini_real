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
  char *inputName = NULL, *witnessSetName = NULL;
  curveDecomp_d C;  //new data type; stores vertices, edges, etc.
  vec_d pi;  //random projection
  witness_set_d Wuser, Wnew;
	int max_deflations = 10;
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
	int precision = 53;
	num_processes = 1;
	int num_var_gps = 0, userHom = 0;
	
	//end parser-bertini essentials
	
	
	
	parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
	if (MPType==2) {
		printf("bertini_real is not equipped for adaptive multiple precision (AMP).\ndefaulting to fixed non-double precision\n");
		mypause();
		MPType=1;
	}
	
	preproc_data PPD;
	setupPreProcData("preproc_data", &PPD);
	
	
	num_var_gps = PPD.num_var_gp;
	num_vars = setupProg(&SLP, precision, MPType); // num_vars includes the number of homogeneous coordinates.
	// the number of homogeneous coordinates is the num_var_gps.
	
	
	printf("parsing witness set\n");
	witnessSetParse(&Wuser,witnessSetName,num_vars);
	witnessSetParse(&Wnew, witnessSetName,num_vars);  // Wnew same as Wuser, except for functions.
	
	write_dehomogenized_coordinates(Wuser, "witness_points_dehomogenized");
	
	
	printf("performing isosingular deflation\n");
	// perform an isosingular deflation
	rV = isosingular_deflation(&num_deflations, &deflation_sequence, inputName, "witness_points_dehomogenized", "bertini", "matlab", max_deflations);
  
	
	
	
	
	
	
	
	
	
	
	// initialize the projection pi.  for now, get random projection.
	init_vec_d(pi,num_vars);
	pi->size = num_vars;
	for (ii=0; ii<num_vars; ii++) {
		get_comp_rand_real_d(&pi->coord[ii]);
	}
	
	
	
	
  
	
	
	
	
	
	printf("checking if component is self-conjugate\n");
	sc = checkSelfConjugate(Wuser,num_vars,"input_deflated");  //later:  could be passed in from user, if we want
	if (sc==0) {
		printf("component is NOT self conjugate\n");
	}
	else{
		printf("component IS self conjugate\n");
	}
	
	
	
	
	
	
	
	
	
	
	parse_input("input_deflated", &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
	
	
	
	//now need to get new system produced by isosingular_deflation into bertini_real's memory.
	
	
	
	//need to convert Wuser to Wnew here
	
	//q: what could change?
	
	//temp answer:  the functions, but not the points, or the slices.
	
	
	
	
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
		computeCurveSelfConj("input_deflated",
												 Wnew,pi,&C,num_vars,num_var_gps,currentSeed);  //This is Dans', at least at first !!!
	}
  
	
	printf("\n*\ndone with case\n*\n");
	
	
	
	// clear memory
	free(inputName);
	free(witnessSetName);
	
	clear_witness_set(Wuser);
	clear_witness_set(Wnew);
	
	
	
  //TMP END
  return 0;
}



