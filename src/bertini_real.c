#include "bertini_real.h"

int main(int argC, char *args[])
/***************************************************************\
* USAGE:                                                        *
* ARGUMENTS:                                                    *
* RETURN VALUES:                                                *
* NOTES:                                                        *
\***************************************************************/
{
  int rV,num_vars,sc;  //1=self-conjugate; 0=not
  char *inputName = NULL, *witnessSetName = NULL;
  curveDecomp_d C;  //new data type; stores vertices, edges, etc.
  vec_d pi;  //random projection
  witness_set_d Wuser, Wnew;
	int max_deflations = 10;
	int num_deflations, *deflation_sequence = NULL;
	int ii;  // counters
	
	
	// essentials for using the bertini parser
	prog_t SLP;
  unsigned int currentSeed;
	int trackType, genType = 0, MPType, userHom, sharpenOnly, needToDiff, remove_temp, useParallelDiff = 0;
  int my_id, num_processes, headnode = 0; // headnode is always 0
	int precision = 53;
	num_processes = 1;
	//end parser-bertini essentials
	
	
	
	////
	//  begin the actual program
	////
	
	
  startup(argC, args, &inputName, &witnessSetName);  //prints the welcome message,
	   //also gets the inputName, witnessSetName
		//default inputName = "input"
		//default witnessSetName = "witness_set"
  



	

	
		parse_input(inputName, &trackType, &MPType, &genType, &userHom, &currentSeed, &sharpenOnly, &needToDiff, &remove_temp, useParallelDiff, my_id, num_processes, headnode);
		MPType = 0;
	
	num_vars = setupProg(&SLP, precision, MPType); // numbers change for MP.
	printf("there are %d variables in the problem\n", num_vars);
	
  
    printf("parsing witness set\n");
    witnessSetParse(&Wuser,witnessSetName,num_vars);
    witnessSetParse(&Wnew,witnessSetName,num_vars);  // Wnew same as Wuser, except for functions.
	
//		getFunctions(&Wuser,inputName);
    printf("done reading witness set\n");
    
    

	printf("performing isosingular deflation\n");
	// perform an isosingular deflation
	rV = isosingular_deflation(&num_deflations, &deflation_sequence, inputName, witnessSetName, "bertini", "matlab", max_deflations);
  
	
	
 printf("checking if component is self-conjugate\n");
 sc = checkSelfConjugate(Wuser,num_vars,inputName);  //later:  could be passed in from user, if we want
	 
	
	if (sc==0) {
		printf("component is NOT self conjugate\n");
	}
	else{
		printf("component IS self conjugate\n");

	}

	

	//need to convert Wuser to Wnew here
	
	//q: what could change?
	

	//temp answer:  the functions, but not the points, or the slices.
//		getFunctions(&Wnew,"input_deflated");
	
	
	//now need to get new system produced by isosingular_deflation into bertini_real's memory.
	
	
	 
     if (sc==0)  //C is not self-conjugate
     {
       //Call non-self-conjugate case code
       printf("entering not-self-conjugate case\n");
       computeCurveNotSelfConj(Wnew, pi, &C,num_vars,"input");//This is Wenrui's !!! 
       printf("Bertini_real found %d vertices (vertex)\n",C.num_V0);
     }
     else
     {
       //Call self-conjugate case code
       printf("entering self-conjugate case\n");
       computeCurveSelfConj(Wnew,pi,&C);  //This is Dans', at least at first !!!
     }    
  
	
	
	
    // clear memory
	free(inputName);
	free(witnessSetName);
	
	

	clear_vec_d(Wuser.L);
	for (ii=0; ii<Wuser.W.num_pts; ii++) {
		clear_point_d(Wuser.W.pts[ii]);
	}
  free(Wuser.W.pts);
	
	
	
  clear_vec_d(Wnew.L);
	for (ii=0; ii<Wnew.W.num_pts; ii++) {
		clear_point_d(Wnew.W.pts[ii]);
	}
  free(Wnew.W.pts);

  //TMP END
  return 0;
}


