

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



#include "fileops.h"
#include "polysolve.h"
#include "data_type.h"
#include "programStartup.h"
#include "witnessSet.h"
#include "checkSelfConjugate.h"

#include "lintolinSolver.h"
#include "multilintolinSolver.h"

#ifndef _SAMPLER_H
#define _SAMPLER_H






int  curve_sampler_startup(int argC, char *args[],
													 char directoryName[], char **inputName,
													 char **witnessSetName, char **RandMatName,
													 char **samplingNamenew,
													 curveDecomp_d *C,
													 vertex_set *V,
													 int *num_vars ,
													 int MPType);





int get_dir_mptype(char ** Dir_Name, int * MPType);

int curve_setup_edges(curveDecomp_d *C,
											char *INfile);

int setup_vertices(vertex_set *V,char *INfile);

int  set_initial_sample_data(sample_data *S, curveDecomp_d C, vertex_set V,
													int num_vars);

void  output_sampling_data(sample_data S,  vertex_set V,
													 char *samplingName,int num_vars, int MPType);
void  set_witness_set_d(witness_set *W, vec_d L,vec_d pts,int num_vars);
void set_witness_set_mp(witness_set *W, vec_mp L,vec_mp pts,int num_vars);

void generate_new_sampling_pts(sample_data *S_new,
															 mat_mp n_minusone_randomizer_matrix,
															 sample_data S_old,
															 curveDecomp_d C,
															 vertex_set *V,
															 witness_set W,
															 int  MPType,
															 sampler_configuration *sampler_options,
															 solver_configuration *solve_options);

//void generate_new_sampling_pts_d(sample_data *S_new,
//																 mat_mp n_minusone_randomizer_matrix,
//																 sample_data S_old,
//																 curveDecomp_d C,
//																 vertex_set *V,
//																 witness_set W,
//																 int  MPType,
//																 sampler_configuration *sampler_options,
//																 solver_configuration *solve_options);

//void generate_new_sampling_pts_mp(sample_data *S_new,
//																	mat_mp n_minusone_randomizer_matrix,
//																	sample_data S_old,
//																	curveDecomp_d C,
//																	vertex_set *V,
//																	witness_set W,
//																	int  MPType,
//																	sampler_configuration *sampler_options,
//																	solver_configuration *solve_options);

void read_rand_matrix(char *INfile, mat_mp n_minusone_randomizer_matrix);


int setup_curve(curveDecomp_d *C,char *INfile, int MPType, char **inputName, char *directoryName);


void estimate_new_projection_value(comp_mp result, vec_mp left, vec_mp right, vec_mp pi);
#endif


void set_initial_refinement_flags(int *num_refinements, int **refine_flags, int **current_indices,
																 sample_data S, vertex_set *V,
																 int current_edge, sampler_configuration *sampler_options);


