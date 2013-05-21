
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#ifndef WITNESS_SET_H_
#define WITNESS_SET_H_

#include "data_type.h"
#include "polysolve.h"
#include "fileops.h"





typedef struct
{
	//  system sys;
  vec_d *L_d;
	vec_d *patch_d;
  point_d *pts_d;
	
	vec_mp *L_mp;
	vec_mp *patch_mp;
  point_mp *pts_mp;
	
	
	
	int codim;
  int comp_num;
	
	
	int incidence_number;
	int num_pts;
	int num_var_gps;
	int num_variables;
	int num_linears;
	int num_patches;
//	int *patch_sizes;
	int MPType; //  will indicate the type of the solve.  both fields will have data, but one will be used.
	
	char ** variable_names;
} witness_set;  //For a single irreducible component!  Need num. irred. decomp. type later.
// end the double types





int witnessSetParse(witness_set *W, char *witness_set_file, int num_vars);


void add_patch_to_witness_set(witness_set *W, vec_mp new_patch);


void add_point_to_witness_set(witness_set *W, vec_mp new_point);

void add_linear_to_witness_set(witness_set *W, vec_mp new_linear);


void merge_witness_sets(witness_set *W_out,witness_set W_left,witness_set W_right);

void init_variable_names(witness_set *W, int num_vars);

void cp_names(witness_set *W_out, witness_set W_in);
void cp_linears(witness_set *W_out, witness_set W_in);
void cp_patches(witness_set *W_out, witness_set W_in);
void cp_witness_set(witness_set *W_out, witness_set W_in);


void init_witness_set(witness_set *W);


void get_variable_names(witness_set *W);

void clear_witness_set(witness_set W);
void print_witness_set_to_screen(witness_set W);




void write_homogeneous_coordinates(witness_set W, char filename[]);
void write_dehomogenized_coordinates(witness_set W, char filename[]);
void write_linears(witness_set W, char filename[]);



#endif
