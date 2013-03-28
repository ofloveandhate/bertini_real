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
	
	rewind(IN);
	
	
	
	///////////////
	//
	//
	//   READ IN SAME FILE AS THE MP DATA
	//
	//
	/////////
	
	

	
  fscanf(IN, "%d", &num_pts);
  W->W_mp.num_pts=num_pts;
  W->W_mp.pts=(point_mp *)bmalloc(W->W_mp.num_pts*sizeof(point_mp));
//  W->num_variables = num_vars;
  
//	size_t retval;
//	printf("%d pts\n",num_pts);
  for (ii=0; ii < num_pts; ii++) {
		
    //initialize the memory
    init_point_mp(W->W_mp.pts[ii],num_vars);
		change_size_vec_mp(W->W_mp.pts[ii],num_vars);
		W->W_mp.pts[ii]->size = num_vars;
    //read the witness points into memory
    for (jj=0; jj < num_vars; jj=jj+1) {
			mpf_inp_str(W->W_mp.pts[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->W_mp.pts[ii]->coord[jj].i, IN, 10);
			
//			mpf_out_str (NULL, 10, retval, W->W_mp.pts[ii]->coord[jj].r);
//			printf(" ");
//			
//			mpf_out_str (NULL, 10, retval, W->W_mp.pts[ii]->coord[jj].i);
//			printf("\n");
//			retval = gmp_printf ("W->W_mp.pts[%d]->coord[%d].r = %Zd\n",ii,jj, W->W_mp.pts[ii]->coord[jj].r);
//			printf("%d is retval\n",retval);
//      fscanf(IN,"%le %le",&r,&i); //  assign the real and imaginary part
//      W->W_mp.pts[ii]->coord[jj].r = r;
//      W->W_mp.pts[ii]->coord[jj].i = i;
      
    }
//		printf("\n");
  }
  
	//	printf("done reading in the %d points\n", W->W_mp.num_pts);
	
	fscanf(IN, "%d %d", &num_linears, &num_vars_in_linears);
//	W->num_linears = num_linears;
	W->L_mp = (vec_mp *)bmalloc(num_linears*sizeof(vec_mp));
	
//	printf("%d %d\n",num_linears,num_vars_in_linears);
  for (ii=0; ii < num_linears; ii++) {
		init_vec_mp(W->L_mp[ii],num_vars_in_linears);
		
		change_size_vec_mp(W->L_mp[ii],num_vars_in_linears);
		W->L_mp[ii]->size = num_vars_in_linears;
    //read the witness linears into memory
    for (jj=0; jj < num_vars_in_linears; jj++) {
			mpf_inp_str(W->L_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->L_mp[ii]->coord[jj].i, IN, 10);
			//      fscanf(IN,"%lf %lf",&r,&i);
//      W->L_mp[ii]->coord[jj].r = r;
//      W->L_mp[ii]->coord[jj].i = i;
    }
  }
  
	
	fscanf(IN, "%d %d", &num_patches, &patch_size);
	
//	W->patch_size = patch_size;
	W->patch_mp = (vec_mp *)bmalloc(num_patches*sizeof(vec_mp));
//	W->num_patches = num_patches;
	
  for (ii=0; ii < num_patches; ii++) {
		init_vec_mp(W->patch_mp[ii],patch_size);
		
		change_size_vec_mp(W->patch_mp[ii],patch_size);
		W->patch_mp[ii]->size = patch_size;
    //read the patch into memory
    for (jj=0; jj < patch_size; jj++) {
			mpf_inp_str(W->patch_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->patch_mp[ii]->coord[jj].i, IN, 10);
//      fscanf(IN,"%lf %lf",&r,&i);
//      W->patch[ii]->coord[jj].r = r;
//      W->patch[ii]->coord[jj].i = i;
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


//
//
//int witnessSetParse_mp(witness_set_mp *W, char *witness_set_file, int num_vars){
//  
////  double i, r;
//	
//  int ii=0, jj=0, num_linears, num_pts;
//	
//  FILE *IN = NULL;
//  
//  IN = safe_fopen_read(witness_set_file);
//  
//  
//  int num_vars_in_linears;
//	
//	
//  fscanf(IN, "%d", &num_pts);
//  W->W.num_pts=num_pts;
//  W->W.pts=(point_mp *)bmalloc(W->W.num_pts*sizeof(point_mp));
//  W->num_variables = num_vars;
//  
//
//	
//  for (ii=0; ii < num_pts; ii++) {
//		
//    //initialize the memory
//    init_point_mp(W->W.pts[ii],num_vars);
//		change_size_vec_mp(W->W.pts[ii],num_vars);
//		W->W.pts[ii]->size = num_vars;
//    //read the witness points into memory
//    for (jj=0; jj < num_vars; jj=jj+1) {
//			mpf_inp_str(W->W.pts[ii]->coord[jj].r, IN, 10);
//			mpf_inp_str(W->W.pts[ii]->coord[jj].i, IN, 10); // reads the input buffer for a base 10 number
//			
////      fscanf(IN,"%le %le",&r,&i); //  assign the real and imaginary part
////      W->W.pts[ii]->coord[jj].r = r;
////      W->W.pts[ii]->coord[jj].i = i;
//      
//    }
//  }
//  
//	//	printf("done reading in the %d points\n", W->W.num_pts);
//	
//	fscanf(IN, "%d %d", &num_linears, &num_vars_in_linears);
//	W->num_linears = num_linears;
//	W->L = (vec_mp *)bmalloc(num_linears*sizeof(vec_mp));
//	
//	
//  for (ii=0; ii < num_linears; ii++) {
//		init_vec_mp(W->L[ii],num_vars_in_linears);
//		
//		change_size_vec_mp(W->L[ii],num_vars_in_linears);
//		W->L[ii]->size = num_vars_in_linears;
//    //read the witness linears into memory
//    for (jj=0; jj < num_vars_in_linears; jj++) {
//			
//			
//			mpf_inp_str(W->L[ii]->coord[jj].r, IN, 10);
//			mpf_inp_str(W->L[ii]->coord[jj].i, IN, 10);
//			
////      fscanf(IN,"%lf %lf",&r,&i);
////      W->L[ii]->coord[jj].r = r;
////      W->L[ii]->coord[jj].i = i;
//    }
//  }
//  
//	
//	int num_patches, patch_size;
//	fscanf(IN, "%d %d", &num_patches, &patch_size);
//	
//	W->patch_size = patch_size;
//	W->patch = (vec_mp *)bmalloc(num_patches*sizeof(vec_mp));
//	W->num_patches = num_patches;
//	
//  for (ii=0; ii < num_patches; ii++) {
//		init_vec_mp(W->patch[ii],patch_size);
//		
//		change_size_vec_mp(W->patch[ii],patch_size);
//		W->patch[ii]->size = patch_size;
//    //read the patch into memory
//    for (jj=0; jj < patch_size; jj++) {
//			mpf_inp_str(W->patch[ii]->coord[jj].r, IN, 10);
//			mpf_inp_str(W->patch[ii]->coord[jj].i, IN, 10);
//
//			
////      fscanf(IN,"%lf %lf",&r,&i);
////      W->patch[ii]->coord[jj].r = r;
////      W->patch[ii]->coord[jj].i = i;
//    }
//  }
//	
//	
//  fclose(IN);
//  
//	
//	//	print_point_to_screen_matlab(W->W.pts[0],"firstpoint");
//	//	print_point_to_screen_matlab(W->L[0],"L");
//	//	print_point_to_screen_matlab(W->patch[0],"patch");
//	//
//	//	mypause();
//	
//  return 0;
//}



