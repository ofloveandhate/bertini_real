#include "witnessSetParse.h"




//use:  call to parse the file witness_set_file, into the struc W.  fills the values W and L
// the functions will be parsed elsewhere.


int witnessSetParse(witness_set_d *W, char *witness_set_file, int num_vars){
  
  double i, r;

  int ii=0, jj=0, num_linears, num_pts;
  FILE *IN = NULL;
  
  IN = fopen(witness_set_file, "r");
  if (IN == NULL)
  {
    printf("\n\nERROR: BertiniReal was unable to open witness set file named %s.\n",witness_set_file);
    bexit(ERROR_CONFIGURATION);
  }
  
  
  
  
  fscanf(IN, "%d", &num_pts);
  W->W.num_pts=num_pts;
  W->W.pts=(point_d *)bmalloc(W->W.num_pts*sizeof(point_d));\
  
  


  for (ii=0; ii < num_pts; ii=ii+1) {
   
    //initialize the memory
    init_point_d(W->W.pts[ii],num_vars);
    //read the witness points into memory
    for (jj=0; jj < num_vars; jj=jj+1) {
      fscanf(IN,"%lf %lf",&r,&i); //  assign the real and imaginary part
      W->W.pts[ii]->coord[jj].r = r;
      W->W.pts[ii]->coord[jj].i = i;
      
    }
  }
  
  
  
  init_vec_d(W->L,num_vars);
  W->L->size=num_vars;

  fscanf(IN, "%d", &num_linears);
  for (ii=0; ii < num_linears; ii=ii+1) {
    //read the witness linears into memory
    for (jj=0; jj < num_vars; jj=jj+1) {
      fscanf(IN,"%lf %lf",&r,&i);
      W->L->coord[jj].r = r;
      W->L->coord[jj].i = i;
    }
  }
  
  fclose(IN);
  

  return 0;
}




