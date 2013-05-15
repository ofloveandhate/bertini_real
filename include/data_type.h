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
#include <signal.h>



#ifndef _DATA_TYPE_H
#define _DATA_TYPE_H

#include "fileops.h"
#include "cascade.h"
#include "polysolve.h"


/*** low-level data types. ***/


//The following lets us use words instead of numbers to indicate vertex type.
enum {CRITICAL=0, NEW=1, MIDPOINT=2, ISOLATED=-1, SAMPLE_POINT=3};



void * br_malloc(size_t size);

void br_exit(int errorCode);

void deliberate_segfault();

typedef struct
{
  char *func;  //symbolic representation of function (straight from input file).
} function;






typedef struct
{
//  system sys;
  vec_d *L;
	vec_d *patch;
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
	int patch_size;
	int MPType; //  will indicate the type of the solve.  both fields will have data, but one will be used.
	
	char ** variable_names;
} witness_set;  //For a single irreducible component!  Need num. irred. decomp. type later.
// end the double types




// CURVE CELL DECOMP DATA TYPES
typedef struct
{
  point_d  pt_d;
  point_mp pt_mp;
	
  comp_d   projVal_d; //Value of projection pi applied to pt;
  comp_mp  projVal_mp;
	
  int      type;  //See enum above.
} vertex;


typedef struct
{
	vertex *vertices;  //Isolated real points.
	
	int num_vertices;
	
	
} vertex_set;

typedef struct
{
  int left;  //index into vertices
  int right; //index into vertices
	int midpt; // index into vertices
} edge;

typedef struct
{
	int num_variables;
	
  edge *edges;
	
	//these counters keep track of the number of things
	
	
	int      num_edges;
	
	
	int			 num_isolated;
  int      num_V0;
  int      num_V1;
  int      num_midpts;
	int			 num_new;
	int			 num_samples;
	

	// these arrays of integers index the vertices.
	int			*isolated_indices;
	int			*V0_indices;
	int			*V1_indices;
	int			*midpt_indices;
	int			*new_indices;
	int			*sample_indices;
	
	vec_mp	pi_mp; // the projection
	vec_d		pi_d; // the projection
	mat_mp n_minus_one_randomizer_matrix;
	
	int MPType;
	
} curveDecomp_d;


typedef struct
{
	int num_variables;
	
	int num_edges;
	
	int *num_samples_each_edge;
	
	int **sample_indices;

}
sample_data;




void init_vertex_set(vertex_set *V);
int add_vertex(vertex_set *V, vertex new_vertex);

int curve_add_vertex(curveDecomp_d *C, vertex_set *V, vertex source_vertex);

int curve_index_in_vertices(curveDecomp_d *C,  vertex_set *V,
														vec_mp testpoint, comp_mp projection_value,
														tracker_config_t T);

int curve_index_in_vertices_with_add(curveDecomp_d *C, vertex_set *V,
																		 vec_mp testpoint, comp_mp projection_value,
																		 tracker_config_t T);

void vertex_set_print_to_screen(vertex_set *V);


void init_vertex(vertex *curr_vertex);

/**
 copy in a vertex
 */
void cp_vertex(vertex *target_vertex, vertex source_vertex);


void init_edge(edge *curr_edge);
//void init_edge_mp(edge *curr_edge, int num_variables);
//void init_edge_d(edge *curr_edge, int num_variables);

void add_edge(curveDecomp_d *C, edge new_edge);


//function prototypes for bertini_real data clearing etc.
void init_curveDecomp_d(curveDecomp_d *C, int MPType);
void norm_of_difference(mpf_t result, vec_mp left, vec_mp right);

void dehomogenize_d(vec_d *result, vec_d dehom_me);
void dehomogenize_mp(vec_mp *result, vec_mp dehom_me);

void dot_product_d(comp_d result, vec_d one, vec_d two);
void dot_product_mp(comp_mp result, vec_mp one, vec_mp two);

void projection_value_homogeneous_input_d(comp_d result, vec_d input, vec_d projection);
void projection_value_homogeneous_input(comp_mp result, vec_mp input, vec_mp projection);


int isSamePoint_inhomogeneous_input_d(point_d left, point_d right);
int isSamePoint_inhomogeneous_input_mp(point_mp left, point_mp right);



int isSamePoint_homogeneous_input_d(point_d left, point_d right);

int isSamePoint_homogeneous_input_mp(point_mp left, point_mp right);

void clear_curveDecomp_d(curveDecomp_d *C);
void clear_sample_d(sample_data *S, int MPType);

void print_point_to_screen_matlab(vec_d M, char name[]);
void print_point_to_screen_matlab_mp(vec_mp M, char name[]);
void print_matrix_to_screen_matlab(mat_d M, char name[]);
void print_matrix_to_screen_matlab_mp(mat_mp M, char name[]);

void print_comp_mp_matlab(comp_mp M,char name[]);


void print_path_retVal_message(int retVal);

/**
retrieves the number of variables from the PPD by taking the sum of the sizes, plus the sum of the types.
*/
int get_num_vars_PPD(preproc_data PPD);




#endif

