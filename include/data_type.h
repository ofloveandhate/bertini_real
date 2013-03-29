#include <stdio.h>
#include <stdlib.h>
#include <stdarg.h>
#include <string.h>
#include <math.h>
#include <gmp.h>
#include <time.h>
#include <float.h>
#include <limits.h>
#include <mpfr.h>
#include <mpf2mpfr.h>




#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

#include "fileops.h"
#include "cascade.h"
#include "polysolve.h"


/*** low-level data types. ***/




typedef struct
{
  int num_pts;
  point_d *pts;
  int dim;
  int deg;
} witness_point_set_d;


typedef struct
{
  char *func;  //symbolic representation of function (straight from input file).
} function;


//typedef struct
//{
//  function *funcs; //probably not used for now.
//  //prog_t slp; //SLP -- we'll use this primarily, at least at first.
//  char *fName; //system can be given in a file -- this is its name.
//}system;

////begin mp types
typedef struct
{
  int num_pts;
  point_mp *pts;
  int dim;
  int deg;
} witness_point_set_mp;
//end mp types

typedef struct
{
//  system sys;
  vec_d *L;
	vec_d *patch;
  witness_point_set_d W;
	
	vec_mp *L_mp;
	vec_mp *patch_mp;
  witness_point_set_mp W_mp;
	
	int num_variables;
	int num_linears;
	int num_patches;
	int patch_size;
	int MPType; //  will indicate the type of the solve.  both fields will have data, but one will be used.
} witness_set_d;  //For a single irreducible component!  Need num. irred. decomp. type later.
// end the double types



//
//
//typedef struct
//{
//  char *func;  //symbolic representation of function (straight from input file).
//}function_mp;
//
//
//typedef struct
//{
//  function_mp *funcs; //probably not used for now.
//  //prog_t slp; //SLP -- we'll use this primarily, at least at first.
//  char *fName; //system can be given in a file -- this is its name.
//}system_mp;
//
//typedef struct
//{
//  system_mp sys;
//  vec_mp *L;
//	vec_mp *patch;
//  witness_point_set_mp W;
//	int num_variables;
//	int num_linears;
//	int num_patches;
//	int patch_size;
//}witness_set_mp;  //For a single irreducible component!  Need num. irred. decomp. type later.
////end mp types



// CURVE CELL DECOMP DATA TYPES
typedef struct
{
  point_d pt;
  comp_d projVal; //Value of projection pi applied to pt; used for easy comparison of projected values later....
  int type;  //See enum below.
  int num_left;  //this and next line track how often this vertex appears in
  int num_right; //edges;  good for sanity check AND determining V0.
}vertex_d;

typedef struct
{
  int left;  //index from V1
  int right; //index from V1
  point_d midpt; //the midPts
  witness_set_d W; //contains functions; IS THIS OVERKILL????  Could just be a system_d....
  vec_d pi;  //projection
}edge_d;

typedef struct
{
  vertex_d *V0;  //Isolated real points.
  vertex_d *V1;  //Critical points AND new non-critical endpoints of edges.
//  vertex_d *midPts;  //Midpoints of edges.
  edge_d *E;
  int      num_V0;
  int      num_V1;
//  int      num_midPts;
  int      num_E;
}curveDecomp_d;

//The following lets us use words instead of numbers to indicate vertex types.
enum {CRITICAL=0, NEW=1, MIDPOINT=2};









//function prototypes for bertini_real data clearing etc.
void merge_witness_sets(witness_set_d *W_out,witness_set_d W_left,witness_set_d W_right);

void cp_patches(witness_set_d *W_out, witness_set_d W_in);
void init_witness_set_d(witness_set_d *W);

void dot_product_d(comp_d result, vec_d one, vec_d two);
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two);

void write_dehomogenized_coordinates(witness_set_d W, char filename[]);
void dehomogenize(vec_d *result, vec_d dehom_me);
void dehomogenize_mp(vec_mp *result, vec_mp dehom_me);

void clear_witness_set(witness_set_d W);
void print_witness_set_to_screen(witness_set_d W);
void print_point_to_screen_matlab(vec_d M, char name[]);
void print_point_to_screen_matlab_mp(vec_mp M, char name[]);
void print_matrix_to_screen_matlab(mat_d M, char name[]);
void print_matrix_to_screen_matlab_mp(mat_mp M, char name[]);

void print_comp_mp_matlab(comp_mp M,char name[]);


void print_path_retVal_message(int retVal);

#endif

