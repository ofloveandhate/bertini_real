#include "witnessSetParse.h"




//use:  call to parse the file witness_set_file, into the struct W.  fills the values W, and L, patch, as well as their _mp counterparts

// the functions will be parsed elsewhere.
// this function has been be rewritten to use MP


int witnessSetParse(witness_set_d *W, char *witness_set_file, const int num_vars){
  
  double i, r;

  int ii=0, jj=0, num_linears, num_pts;
	
  FILE *IN = NULL;
  
  IN = safe_fopen_read(witness_set_file);
  
  
  int num_vars_in_linears;
	
	int codim, comp_num;
	int num_patches, patch_size;

	
	
		

	
  fscanf(IN, "%d %d %d", &num_pts, &codim, &comp_num);
	W->codim = codim;
	W->comp_num = comp_num;
  W->W.num_pts = W->W_mp.num_pts = num_pts;
  W->W.pts=(point_d *)bmalloc(W->W.num_pts*sizeof(point_d));
	W->W_mp.pts=(point_mp *)bmalloc(W->W_mp.num_pts*sizeof(point_mp));

	
  W->num_variables = num_vars;

  for (ii=0; ii < num_pts; ii++) {
    init_point_mp2(W->W_mp.pts[ii],num_vars,1024); init_point_d(W->W.pts[ii],num_vars);
		W->W_mp.pts[ii]->size = W->W.pts[ii]->size = num_vars;
    //read the witness points into memory
    for (jj=0; jj < num_vars; jj=jj+1) {
			mpf_inp_str(W->W_mp.pts[ii]->coord[jj].r, IN, 10); // 10 is the base
			mpf_inp_str(W->W_mp.pts[ii]->coord[jj].i, IN, 10);      
    }
		vec_mp_to_d(W->W.pts[ii],W->W_mp.pts[ii]);
  }
  
	

	fscanf(IN, "%d %d", &num_linears, &num_vars_in_linears);  scanRestOfLine(IN);
//	printf("%d %d\n",num_linears, num_vars_in_linears);
	W->num_linears = num_linears;
	W->L = (vec_d *)bmalloc(num_linears*sizeof(vec_d));
	W->L_mp = (vec_mp *)bmalloc(num_linears*sizeof(vec_mp));
	
  for (ii=0; ii < num_linears; ii++) {
		init_vec_mp2(W->L_mp[ii],num_vars_in_linears,1024); init_vec_d(W->L[ii],num_vars_in_linears);
		
		W->L[ii]->size =  W->L_mp[ii]->size = num_vars_in_linears;
		
    //read the witness linears into memory
    for (jj=0; jj < num_vars_in_linears; jj++) {
			mpf_inp_str(W->L_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->L_mp[ii]->coord[jj].i, IN, 10);
    }
		vec_mp_to_d(W->L[ii],W->L_mp[ii]);
  }
  
	
	fscanf(IN, "%d %d", &num_patches, &patch_size); scanRestOfLine(IN);
	W->patch_size = patch_size;
	W->num_patches = num_patches;
	
	W->patch = (vec_d *)bmalloc(num_patches*sizeof(vec_d));
	W->patch_mp = (vec_mp *)bmalloc(num_patches*sizeof(vec_mp));
	
  for (ii=0; ii < num_patches; ii++) {
		init_vec_mp2(W->patch_mp[ii],patch_size,1024);//default max_prec is 1024
		init_vec_d(W->patch[ii],patch_size);
		
		W->patch_mp[ii]->size = W->patch[ii]->size = patch_size;
		
    //read the patch into memory
    for (jj=0; jj < patch_size; jj++) {
			mpf_inp_str(W->patch_mp[ii]->coord[jj].r, IN, 10);
			mpf_inp_str(W->patch_mp[ii]->coord[jj].i, IN, 10);
    }
		vec_mp_to_d(W->patch[ii],W->patch_mp[ii]);
  }
	
  fclose(IN);
	
	
//TODO: check for mem leaks here.


  return 0;
}


//requires the witness set to have set the number of variables.
void get_variable_names(witness_set_d *W){
	int ii;
	
	
	FILE *IN = safe_fopen_read("names.out");
	
	init_variable_names(W,W->num_variables);
	for (ii=0; ii<W->num_variables; ++ii){
	fscanf(IN,"%s\n",W->variable_names[ii]);
	}
	
	fclose(IN);
	
	
}

