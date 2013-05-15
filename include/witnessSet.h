
#include <stdio.h>
#include <gmp.h>
#include <mpfr.h>

#ifndef WITNESS_SET_H_
#define WITNESS_SET_H_

#include "data_type.h"
#include "polysolve.h"
#include "fileops.h"


int witnessSetParse(witness_set *W, char *witness_set_file, int num_vars);



void merge_witness_sets(witness_set *W_out,witness_set W_left,witness_set W_right);

void init_variable_names(witness_set *W, int num_vars);

void cp_names(witness_set *W_out, witness_set W_in);
void cp_linears(witness_set *W_out, witness_set W_in);
void cp_patches(witness_set *W_out, witness_set W_in);
void cp_witness_set(witness_set *W_out, witness_set W_in);
void init_witness_set_d(witness_set *W);
void get_variable_names(witness_set *W);

void clear_witness_set(witness_set W);
void print_witness_set_to_screen(witness_set W);




void write_homogeneous_coordinates(witness_set W, char filename[]);
void write_dehomogenized_coordinates(witness_set W, char filename[]);
void write_linears(witness_set W, char filename[]);



#endif
