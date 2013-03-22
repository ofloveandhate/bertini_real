#include "witnessSetParse.h"




//use:  call to parse the file witness_set_file, into the struc W.  fills the values W and L
// the functions will be parsed elsewhere.

// this function should be rewritten to use MP
int witnessSetParse(witness_set_d *W, char *witness_set_file, const int num_vars){
  
  double i, r;

  int ii=0, jj=0, num_linears, num_pts;
	
  FILE *IN = NULL;
  
  IN = safe_fopen_read(witness_set_file);
  
  
  int num_vars_in_linears;
	
	
  fscanf(IN, "%d", &num_pts);
  W->W.num_pts=num_pts;
  W->W.pts=(point_d *)bmalloc(W->W.num_pts*sizeof(point_d));
  W->num_variables = num_vars;
  


  for (ii=0; ii < num_pts; ii++) {
   
    //initialize the memory
    init_point_d(W->W.pts[ii],num_vars);
		change_size_vec_d(W->W.pts[ii],num_vars);
		W->W.pts[ii]->size = num_vars;
    //read the witness points into memory
    for (jj=0; jj < num_vars; jj=jj+1) {
      fscanf(IN,"%le %le",&r,&i); //  assign the real and imaginary part
      W->W.pts[ii]->coord[jj].r = r;
      W->W.pts[ii]->coord[jj].i = i;
      
    }
  }
  
//	printf("done reading in the %d points\n", W->W.num_pts);
	
	fscanf(IN, "%d %d", &num_linears, &num_vars_in_linears);
	W->num_linears = num_linears;
	W->L = (vec_d *)bmalloc(num_linears*sizeof(vec_d));


  for (ii=0; ii < num_linears; ii++) {
		init_vec_d(W->L[ii],num_vars_in_linears);

		change_size_vec_d(W->L[ii],num_vars_in_linears);
		W->L[ii]->size = num_vars_in_linears;
    //read the witness linears into memory
    for (jj=0; jj < num_vars_in_linears; jj++) {
      fscanf(IN,"%lf %lf",&r,&i);
      W->L[ii]->coord[jj].r = r;
      W->L[ii]->coord[jj].i = i;
    }
  }
  
	
	int num_patches, patch_size;
	fscanf(IN, "%d %d", &num_patches, &patch_size);

	W->patch_size = patch_size;
	W->patch = (vec_d *)bmalloc(num_patches*sizeof(vec_d));
	W->num_patches = num_patches;
	
  for (ii=0; ii < num_patches; ii++) {
		init_vec_d(W->patch[ii],patch_size);
		
		change_size_vec_d(W->patch[ii],patch_size);
		W->patch[ii]->size = patch_size;
    //read the patch into memory
    for (jj=0; jj < patch_size; jj++) {
      fscanf(IN,"%lf %lf",&r,&i);
      W->patch[ii]->coord[jj].r = r;
      W->patch[ii]->coord[jj].i = i;
    }
  }
	
	
  fclose(IN);
  
	
//	print_point_to_screen_matlab(W->W.pts[0],"firstpoint");
//	print_point_to_screen_matlab(W->L[0],"L");
//	print_point_to_screen_matlab(W->patch[0],"patch");
//
//	mypause();
	
  return 0;
}




