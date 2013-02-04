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
  char *inputName = NULL, *startName = NULL;
  curveDecomp_d C;  //new data type; stores vertices, edges, etc.
  vec_d pi;  //random projection
  witness_set_d Wuser, Wnew;
	int max_deflations = 10;
	int num_deflations, *deflation_sequence = NULL;
	
  startup(argC, args, &inputName, &startName);  //prints the welcome message,
	   //also gets the inputName, startName
  

    num_vars=2; //hardcoded here until we program to parse the function file.
    
    printf("parsing witness set\n");
    witnessSetParse(&Wuser,startName,num_vars);
    
    printf("done reading witness set\n");
    
    


    
     printf("checking if component is self-conjugate\n");
     sc = checkSelfConjugate(Wuser,num_vars,"input");  //later:  could be passed in from user, if we want
	 
	
	if (sc==0) {
		printf("component is NOT self conjugate\n");
	}
	else{
		printf("component IS self conjugate\n");

	}
	return 0;
	
	printf("performing isosingular deflation\n");
	// perform an isosingular deflation
	rV = isosingular_deflation(&num_deflations, &deflation_sequence, inputName, startName, "bertini", "matlab", max_deflations);
	//need to convert Wuser to Wnew here
	
	 
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
	free(startName);
  
  //TMP: clear L W
	
	clear_vec_d(Wuser.L);
  clear_point_d(Wuser.W.pts[0]);
  clear_point_d(Wuser.W.pts[1]);
  free(Wuser.W.pts);
	
	
  clear_vec_d(Wnew.L);
  clear_point_d(Wnew.W.pts[0]);
  clear_point_d(Wnew.W.pts[1]);
  free(Wnew.W.pts);

  //TMP END
  return 0;
}


