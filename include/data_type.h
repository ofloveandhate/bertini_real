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
	
	int codim;
  int comp_num;
	
	char ** variable_names;
	int incidence_number;
	int num_var_gps;
	int num_variables;
	int num_linears;
	int num_patches;
	int patch_size;
	int MPType; //  will indicate the type of the solve.  both fields will have data, but one will be used.
} witness_set_d;  //For a single irreducible component!  Need num. irred. decomp. type later.
// end the double types




// CURVE CELL DECOMP DATA TYPES
typedef struct
{
  point_d pt;
  point_mp pt_mp;
  comp_d projVal; //Value of projection pi applied to pt; used for easy comparison of projected values later....
  comp_mp projVal_mp;
  int type;  //See enum below.
  int num_left;  //this and next line track how often this vertex appears in
  int num_right; //edges;  good for sanity check AND determining V0.
}vertex_d;

typedef struct
{
  int left;  //index from V1
  int right; //index from V1
  point_d midpt; //the midpoint
  point_mp midpt_mp;//the midpoint in multiple precision
//  witness_set_d W; //contains functions; IS THIS OVERKILL????  Could just be a system_d....
  vec_d pi;  //projection
  vec_mp pi_mp;
}edge_d;

typedef struct
{
  vertex_d *V0;  //Isolated real points.
  vertex_d *V1;  //Critical points AND new non-critical endpoints of edges.
	//  vertex_d *midPts;  //Midpoints of edges.
  edge_d *edges;
	
  int      num_V0;
  int      num_V1;
  int      num_midpts;
	
	int *V0_indices;
	int *V1_indices;
	int *midpt_indices;
	
  int      num_edges;
}curveDecomp_d;


typedef struct
{
  vec_d    **vertices;  //points
  vec_d    *proj_vertices; //projection of points
  vec_mp   **vertices_mp;
  vec_mp   *proj_vertices_mp;
  int      **refine;
  int      num_edges;// number of edges
  int      *num_pts;// number of points
}sample_d;






#define MAX_STRLEN 200

typedef struct
{
	int user_projection;
	char *projection_filename;
	
	int user_randomization;
	char *randomization_filename;
	
	char *input_filename;
	
	char *witness_set_filename;
	
	char *input_deflated_filename;
} program_configuration;


void init_configuration(program_configuration *options);
void clear_configuration(program_configuration *options);


//The following lets us use words instead of numbers to indicate vertex types.
enum {CRITICAL=0, NEW=1, MIDPOINT=2};





int index_in_V1(curveDecomp_d *C, vec_mp testpoint, comp_mp projection_value, tracker_config_t T, int sidedness);
void add_point_to_V0(curveDecomp_d *C, vertex_d new_vertex);
void add_point_to_V1(curveDecomp_d *C, vertex_d new_vertex);

void init_vertex_mp(vertex_d *curr_vertex);
void cp_vertex_mp(vertex_d *target_vertex, vertex_d new_vertex);
void add_edge_mp(curveDecomp_d *C, edge_d new_edge);


//function prototypes for bertini_real data clearing etc.
void init_curveDecomp_d(curveDecomp_d *C);
void merge_witness_sets(witness_set_d *W_out,witness_set_d W_left,witness_set_d W_right);

void init_variable_names(witness_set_d *W, int num_vars);

void cp_names(witness_set_d *W_out, witness_set_d W_in);
void cp_linears(witness_set_d *W_out, witness_set_d W_in);
void cp_patches(witness_set_d *W_out, witness_set_d W_in);
void cp_witness_set(witness_set_d *W_out, witness_set_d W_in);
void init_witness_set_d(witness_set_d *W);


void dehomogenize(vec_d *result, vec_d dehom_me);
void dehomogenize_mp(vec_mp *result, vec_mp dehom_me);

void dot_product_d(comp_d result, vec_d one, vec_d two);
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two);


void write_dehomogenized_coordinates(witness_set_d W, char filename[]);
void write_linears(witness_set_d W, char filename[]);

void clear_curveDecomp_d(curveDecomp_d *C, int MPType);
void clear_sample_d(sample_d *S, int MPType);
void clear_witness_set(witness_set_d W);
void print_witness_set_to_screen(witness_set_d W);
void print_point_to_screen_matlab(vec_d M, char name[]);
void print_point_to_screen_matlab_mp(vec_mp M, char name[]);
void print_matrix_to_screen_matlab(mat_d M, char name[]);
void print_matrix_to_screen_matlab_mp(mat_mp M, char name[]);

void print_comp_mp_matlab(comp_mp M,char name[]);


void print_path_retVal_message(int retVal);



void endgamedata_to_endpoint(post_process_t *endPoint, endgame_data_t *EG);

int BRfindSingularSolns(post_process_t *endPoints,
												int num_sols, int num_vars,
												tracker_config_t *T );
int BRfindFiniteSolns(post_process_t *endPoints, int num_sols, int num_vars,
											tracker_config_t *T );

void BRpostProcessing(post_process_t *endPoints, witness_set_d *W_new, int num_pts, preproc_data preProcData, tracker_config_t *T);


void insert_randomization_matrix_witness_data(int rows, int cols, int codim_index);



#endif

