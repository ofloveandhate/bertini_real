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
	
	
  startup(argC, args, &inputName, &witnessSetName);  //prints the welcome message,
	   //also gets the inputName, witnessSetName
		//default inputName = "input"
		//default witnessSetName = "witness_set"
  

    num_vars=2; //hardcoded here until we program to parse the function file.
    
    printf("parsing witness set\n");
    witnessSetParse(&Wuser,witnessSetName,num_vars);
    witnessSetParse(&Wnew,witnessSetName,num_vars);  // Wnew same as Wuser, except for functions.
	
//		getFunctions(&Wuser);
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
//		getFunctions(&Wnew);
	
	
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


